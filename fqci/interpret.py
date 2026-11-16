"""Pair reports, run profile rules, and derive fastp parameters."""
from __future__ import annotations

import math
import re
import shlex
from dataclasses import asdict, dataclass, field

from . import metrics as M
from .parser import FastQCReport, read_label_is_explicit, sample_and_read
from .rules import SEVERITY_ORDER, Finding, Profile, mark_secondary, render, run_rules

ADAPTERS = {
    "truseq": ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
    "nextera": ("CTGTCTCTTATACACATCT", "CTGTCTCTTATACACATCT"),
}
TAKE_MAX = {"length_required", "cut_mean_quality", "cut_window_size"}
TAKE_MIN = {"n_base_limit", "unqualified_percent_limit"}


@dataclass
class SampleResult:
    sample: str
    files: dict
    metrics: dict
    findings: list[Finding]
    fastp_params: dict
    fastp_command: str
    param_reasons: dict = field(default_factory=dict)
    undetermined: bool = False

    def primary(self) -> list[Finding]:
        return [f for f in self.findings
                if not f.expected and not f.secondary_to and f.severity != "info"]

    @property
    def verdict(self) -> str:
        if self.undetermined:
            return "n/a (unassigned reads)"
        prim = self.primary()
        if not prim:
            return "PASS"
        w = self.worst.upper()
        stuck = sorted({f"{f.rule_id} ({f.severity})" for f in prim if not f.fixable})
        if not stuck:
            return f"{w}: fixable by trimming"
        return f"{w}; not fixable by trimming: {', '.join(stuck)}"

    @property
    def worst(self) -> str:
        sev = [f.severity for f in self.primary()]
        return min(sev, key=lambda s: SEVERITY_ORDER[s]) if sev else "pass"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["worst"] = self.worst
        d["verdict"] = self.verdict
        return d


def natural_key(s: str):
    """Vc2 < Vc10 < Vc_1786: numbers compare as numbers."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def depth_dup_expectation(n_reads: int, genome_size: float) -> tuple[float, float]:
    """Expected FastQC deduplicated % from random start-position collisions alone.

    Reads land on 2*G possible (start, strand) positions; with Poisson rate
    lam = n / 2G, the fraction of reads that are distinct is (1 - e^-lam) / lam.
    """
    lam = n_reads / (2 * genome_size)
    return lam, (100.0 if lam == 0 else 100 * (1 - math.exp(-lam)) / lam)


def group(reports: list[FastQCReport]) -> dict[str, dict[int, FastQCReport]]:
    """Pair R1/R2 reports by sample name.

    Bare _1/_2 names (no R) are only treated as mates when the read counts agree;
    otherwise they are two single-end samples (e.g. isolates named iso_1 and iso_2)."""
    out: dict[str, dict[int, FastQCReport]] = {}
    for r in reports:
        s, rd = sample_and_read(r.filename)
        slot = out.setdefault(s, {})
        key = rd or 1
        if key in slot:
            raise ValueError(f"duplicate read {key} for sample {s}: {slot[key].source} and {r.source}")
        slot[key] = r
    for s in list(out):
        slot = out[s]
        if 1 in slot and 2 in slot and not (read_label_is_explicit(slot[1].filename)
                                            and read_label_is_explicit(slot[2].filename)):
            n1, n2 = (M.compute(slot[k])["n_reads"] for k in (1, 2))
            if n1 != n2:
                del out[s]
                for k in (1, 2):
                    out[_stem(slot[k].filename)] = {1: slot[k]}
    return out


def _stem(filename: str) -> str:
    return re.sub(r"(\.(fastq|fq))?(\.gz)?$", "", filename, flags=re.I)


def adapter_family(*ms: dict) -> str:
    """Dimer sequences are direct evidence and win; Adapter Content columns decide otherwise."""
    seqs = [s["seq"] for m in ms for s in m["overrep"]]
    if any(q.startswith(M.DIMER_PREFIXES[:2]) for q in seqs):
        return "truseq"
    if any(q.startswith(M.DIMER_PREFIXES[2:]) for q in seqs):
        return "nextera"
    nextera = sum((m.get("adapter_last") or {}).get("Nextera Transposase Sequence", 0) for m in ms)
    truseq = sum((m.get("adapter_last") or {}).get("Illumina Universal Adapter", 0) for m in ms)
    return "nextera" if nextera > truseq else "truseq"


def _num(v):
    if isinstance(v, str):
        try:
            f = float(v)
            return int(f) if f.is_integer() else f
        except ValueError:
            return v
    return v


def interpret_sample(sample: str, reads: dict[int, FastQCReport], profile: Profile,
                     genome_size: float | None = None, expected_gc: float | None = None,
                     threads: int = 4) -> SampleResult:
    ms = {k: M.compute(r) for k, r in sorted(reads.items())}
    paired = 1 in ms and 2 in ms
    p = dict(profile.params)
    p["expected_gc"] = expected_gc
    first = ms[min(ms)]
    min_len = max(p["min_len_floor"],
                  min(p.get("min_len_cap", 10**9), round(p["min_len_frac"] * first["read_len"])))
    ceiling = max(1, int(0.8 * first["read_len"]))  # never discard (nearly) every read
    p["min_len_capped"] = min_len > ceiling
    p["min_len"] = min(min_len, ceiling)
    fam = adapter_family(*ms.values())
    p["adapter_r1"], p["adapter_r2"] = ADAPTERS[fam]
    p["adapter_family"] = fam

    p["no_genome_size"] = not genome_size
    p["purpose"] = profile.purpose or "n/a"
    p["unit"] = "read pairs" if paired else "reads"
    findings: list[Finding] = []
    for k, m in ms.items():
        if genome_size and m.get("dedup_pct") is not None:
            m["reads_per_start"], m["dedup_expected"] = depth_dup_expectation(m["n_reads"], genome_size)
            m["dedup_ratio"] = m["dedup_pct"] / m["dedup_expected"]
        findings += run_rules(profile, {**p, **m}, "read", read=f"R{k}")

    pair_ns = dict(p)
    for k, m in ms.items():
        pair_ns.update({f"r{k}_{kk}": v for kk, v in m.items()})
    dimer = max(m["dimer_pct"] for m in ms.values())
    pair_ns["pair_dimer_pct"] = dimer
    pair_ns["adapter_max_any"] = max(m.get("adapter_max_pct") or 0 for m in ms.values())
    if not paired:
        pair_ns["r1_dimer_pct"] = first["dimer_pct"]
        pair_ns["r2_dimer_pct"] = 0.0
        pair_ns["r2_n_reads"] = first["n_reads"]
    pair_ns.setdefault("r1_observed_i7", None)
    pair_ns.setdefault("r1_observed_index_name", None)
    pairs = first["n_reads"]
    pair_ns["usable_pairs"] = pairs * (1 - dimer / 100)
    if genome_size:
        lens = [m.get("mean_len") or m["read_len"] for m in ms.values()]
        pair_ns["mean_len_used"] = sum(lens) / len(lens)
        bases = pair_ns["usable_pairs"] * sum(lens)
        raw = bases / genome_size
        f = p["depth_factor_low"] if raw < p["depth_low_cutoff"] else p["depth_factor"]
        pair_ns["est_coverage_raw"] = raw
        pair_ns["depth_factor_used"] = f
        pair_ns["est_coverage"] = raw * f
        pair_ns["genome_size_mb"] = genome_size / 1e6
    else:
        pair_ns["est_coverage"] = None
    findings += run_rules(profile, pair_ns, "pair")
    mark_secondary(profile, findings)
    findings.sort(key=lambda f: (f.expected, bool(f.secondary_to), SEVERITY_ORDER[f.severity], f.rule_id))

    # --- merge actions
    params: dict = {}
    reasons: dict[str, list[str]] = {}

    def add(k, v, why):
        v = _num(v)
        if not paired and k in ("adapter_sequence_r2", "detect_adapter_for_pe", "trim_tail2"):
            return
        if k in params and k in TAKE_MAX:
            params[k] = max(params[k], v)
        elif k in params and k in TAKE_MIN:
            params[k] = min(params[k], v)
        else:
            params[k] = v
        reasons.setdefault(k, [])
        if why not in reasons[k]:
            reasons[k].append(why)

    for k, v in profile.base_actions.items():
        add(k, render(v, pair_ns, strict=True), f"profile default ({profile.name})")
    for f in findings:
        if f.expected or not f.fixable:
            continue
        for k, v in f.actions.items():
            add(k, v, f.rule_id + (" (secondary)" if f.secondary_to else ""))
    if "length_required" not in params:
        add("length_required", p["min_len"], "minimum length")
    if p["min_len_capped"]:
        reasons["length_required"].append("capped at 80% of the read length")

    cmd = build_fastp(sample, reads, params, paired, threads)
    return SampleResult(sample=sample, files={f"R{k}": r.filename for k, r in sorted(reads.items())},
                        metrics={f"R{k}": {kk: v for kk, v in m.items()} for k, m in ms.items()} |
                        {"pair": {k: v for k, v in pair_ns.items() if not k.startswith(("r1_", "r2_"))}},
                        findings=findings, fastp_params=params, fastp_command=cmd, param_reasons=reasons,
                        undetermined=sample.lower().startswith("undetermined"))


ORDER = ["phred64", "adapter_sequence", "adapter_sequence_r2", "detect_adapter_for_pe",
         "trim_poly_g", "trim_poly_x", "trim_tail1", "trim_tail2", "cut_tail", "cut_right", "cut_window_size", "cut_mean_quality",
         "qualified_quality_phred", "unqualified_percent_limit", "n_base_limit",
         "low_complexity_filter", "length_required"]


FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)


def build_fastp(sample, reads, params, paired, threads: int = 4) -> str:
    """fastp command for one sample. Returns a '# ...' comment instead when fastp cannot
    read the input (e.g. BAM). Output names never equal an input name."""
    keys = sorted(reads)
    inputs = [reads[k].filename for k in keys]
    bad = [f for f in inputs if not FASTQ_RE.search(f)]
    if bad:
        return f"# {sample}: skipped - fastp reads FASTQ (optionally .gz); input is {', '.join(bad)}"
    outs = []
    for k in keys:
        o = f"{sample}_R{k}.trim.fastq.gz"
        if o in inputs:  # re-trimming fqci's own output: never overwrite the input
            o = f"{sample}_R{k}.retrim.fastq.gz"
        outs.append(o)
    parts = ["fastp", "-i", inputs[0], "-o", outs[0]]
    if paired:
        parts += ["-I", inputs[1], "-O", outs[1]]
    for k in sorted(params, key=lambda k: ORDER.index(k) if k in ORDER else 99):
        v = params[k]
        if v is True:
            parts.append(f"--{k}")
        elif v not in (False, None):
            parts += [f"--{k}", str(v)]
    parts += ["--json", f"{sample}.fastp.json", "--html", f"{sample}.fastp.html", "--thread", str(threads)]
    return " ".join(shlex.quote(x) for x in parts)

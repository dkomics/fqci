"""Before/after comparison of two FastQC report sets.

Besides per-sample tables, it predicts how many reads trimming should have removed
(mostly adapter dimers) and flags samples that lost more - the signature of an
over-aggressive setting such as --cut_right or a high --length_required.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from . import metrics as M
from .interpret import group, interpret_sample, natural_key
from .parser import FastQCReport
from .rules import Profile

KEYS = [("n_reads", "Reads", "{:,}"), ("dimer_pct", "Dimers %", "{:.1f}"),
        ("q_drop_lq20", "LQ<20 from", "{}"), ("pct_reads_meanq_lt20", "meanQ<20 %", "{:.1f}"),
        ("dedup_pct", "Dedup %", "{:.1f}"), ("adapter_max_pct", "Adapter %", "{:.2f}"),
        ("polya_pct", "PolyA %", "{:.1f}"), ("comp_bias_body", "Comp. bias", "{:.0f}")]
MODULES = ["Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores",
           "Per base sequence content", "Per sequence GC content", "Sequence Duplication Levels",
           "Overrepresented sequences", "Adapter Content"]


def _f(fmt, v):
    return "n/a" if v is None else fmt.format(v)


@dataclass
class LossCheck:
    sample: str
    kept_pct: float
    expected_loss: float
    actual_loss: float
    excess: float
    flagged: bool
    cause: str


def fastp_cause(js: dict, diagnose: bool = True) -> str:
    """Explain where fastp's removed reads went, from its JSON report.

    With diagnose=False only the breakdown is given: for a sample whose loss matches its
    dimers, 'too short' is mostly dimers trimmed to nothing, and advice would mislead."""
    fr = js.get("filtering_result", {})
    removed = {k: v for k, v in fr.items() if k != "passed_filter_reads" and v}
    top = max(removed, key=lambda k: removed[k]) if removed else None
    total = sum(fr.values()) or 1
    if not removed:
        return "fastp removed no reads."
    parts = [f"{k.replace('_reads', '').replace('_', ' ')} {100 * v / total:.1f}%"
             for k, v in sorted(removed.items(), key=lambda kv: -kv[1])]
    cmd = js.get("command", "")
    msg = "fastp removed: " + ", ".join(parts) + "."
    if not diagnose:
        return msg
    if top == "too_short_reads":
        msg += " Most were too short after trimming."
        if "--cut_right" in cmd:
            msg += (" The command used --cut_right, which cuts each read at its first low-quality"
                    " window; --cut_tail trims only the poor 3' end (on a 2x151 MiSeq run: 59% of"
                    " pairs kept with --cut_right and a 76 bp minimum, 95% with --cut_tail and 50 bp).")
        else:
            msg += " Check --length_required against the read-length distribution."
    elif top == "low_quality_reads":
        msg += " Most failed the per-read quality filter (--qualified_quality_phred / --unqualified_percent_limit)."
    elif top == "adapter_dimer_reads":
        msg += " Most were adapter dimers."
    return msg


def loss_check(sample, before: dict, after: dict, profile: Profile, fastp_dir: Path | None) -> LossCheck:
    mb = {k: M.compute(r) for k, r in before.items()}
    ma = {k: M.compute(r) for k, r in after.items()}
    r = min(set(mb) & set(ma))  # caller guarantees at least one read in both sets
    kept = 100 * ma[r]["n_reads"] / mb[r]["n_reads"] if mb[r]["n_reads"] else 0.0
    expected = max(m["dimer_pct"] for m in mb.values())
    actual = 100 - kept
    excess = actual - expected
    flagged = excess > profile.params.get("trim_loss_allowance", 8)
    cause = ""
    if fastp_dir is not None:
        jf = fastp_dir / f"{sample}.fastp.json"
        if jf.exists():
            try:
                cause = fastp_cause(json.loads(jf.read_text()), diagnose=flagged)
            except (ValueError, OSError) as e:
                cause = f"could not read {jf.name}: {e}"
        else:
            cause = f"no {jf.name} in {fastp_dir}"
    if flagged and not cause:
        cause = ("Likely causes: --length_required too high for the read lengths, or --cut_right "
                 "quality trimming. Pass --fastp-dir to see fastp's breakdown.")
    return LossCheck(sample, kept, expected, actual, excess, flagged, cause)


def compare(before: list[FastQCReport], after: list[FastQCReport], profile: Profile,
            fastp_dir: str | Path | None = None, **kw) -> str:
    gb, ga = group(before), group(after)
    fd = Path(fastp_dir) if fastp_dir else None
    common = sorted(set(gb) & set(ga), key=natural_key)
    shared = [s for s in common if set(gb[s]) & set(ga[s])]
    no_reads = [s for s in common if s not in shared]
    L = ["# Before/after trimming", ""]
    missing = sorted(set(gb) ^ set(ga), key=natural_key)
    if missing:
        L += [f"Samples present in only one set (skipped): {', '.join(missing)}", ""]
    if no_reads:
        L += [f"Samples with no read (R1/R2) in both sets (skipped): {', '.join(no_reads)}", ""]

    checks = {s: loss_check(s, gb[s], ga[s], profile, fd) for s in shared}
    results = {s: interpret_sample(s, ga[s], profile, **kw) for s in shared}
    remaining = {s: [f for f in results[s].findings if not f.expected and f.severity in ("fail", "warn")
                     and f.rule_id not in ("COVERAGE", "INDEX_SEEN", "LOW_INPUT_GC_BIAS")]
                 for s in shared}

    allow = profile.params.get("trim_loss_allowance", 8)
    flagged = [c for c in checks.values() if c.flagged]
    L += ["## Run summary", ""]
    if shared:
        kept_all = [c.kept_pct for c in checks.values()]
        L.append(f"{len(shared)} samples compared. Reads kept: {min(kept_all):.1f}-{max(kept_all):.1f}%.")
        done = sum(1 for s in shared if not remaining[s] and not checks[s].flagged)
        L.append(f"{done} of {len(shared)} have nothing left to fix"
                 + (" (over-trimmed samples are not counted as fixed)." if flagged else "."))
    if flagged:
        L += ["", f"**Trimming removed more than expected in {len(flagged)} of {len(shared)} samples** "
              f"(more than {allow} points beyond the predicted dimer loss). Review the trimming "
              "settings before using these reads."]
    L += ["", "| Sample | Reads kept % | Expected loss % | Actual loss % | Excess | Check |",
          "|---|---:|---:|---:|---:|---|"]
    for s in shared:
        c = checks[s]
        L.append(f"| {s} | {c.kept_pct:.1f} | {c.expected_loss:.1f} | {c.actual_loss:.1f} | "
                 f"{c.excess:+.1f} | {'TOO MUCH REMOVED' if c.flagged else 'ok'} |")
    L += ["", "Expected loss is the adapter-dimer share of the raw reads. Excess is actual minus expected.", ""]

    for s in shared:
        L += [f"## {s}", ""]
        c = checks[s]
        if c.flagged:
            L += [f"**[WARN] Trimming removed {c.actual_loss:.1f}% of reads; about {c.expected_loss:.1f}% "
                  f"was expected.** {c.cause}", ""]
        elif c.cause:
            L += [c.cause, ""]
        L += ["| Read | Metric | Before | After |", "|---|---|---:|---:|"]
        for rd in sorted(set(gb[s]) & set(ga[s])):
            mb, ma = M.compute(gb[s][rd]), M.compute(ga[s][rd])
            for k, lab, fmt in KEYS:
                L.append(f"| R{rd} | {lab} | {_f(fmt, mb.get(k))} | {_f(fmt, ma.get(k))} |")
            kept = 100 * ma["n_reads"] / mb["n_reads"] if mb["n_reads"] else 0
            L.append(f"| R{rd} | Reads kept % | | {kept:.1f} |")
            for mod in MODULES:
                L.append(f"| R{rd} | {mod} | {mb['status'].get(mod, 'n/a')} | {ma['status'].get(mod, 'n/a')} |")
        rem = remaining[s]
        if rem:
            L += ["", "Still flagged after trimming:", ""]
        elif c.flagged:
            L += ["", "No FastQC findings remain, but too many reads were removed (see above): "
                  "re-trim before using these reads.", ""]
        else:
            L += ["", "Nothing left to fix after trimming.", ""]
        L += [f"- [{f.severity.upper()}] {f.title} ({f.read or 'pair'}): {f.explain}" for f in rem]
        L.append("")
    return "\n".join(L)

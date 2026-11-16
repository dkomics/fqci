"""Turn a parsed FastQC report into numbers that rules can reason about.

Every metric is a plain float/int/str/bool or None (None = not computable,
and any rule referring to it is skipped rather than guessed).
"""
from __future__ import annotations

import re
from typing import Any
from collections import Counter

from .parser import FastQCReport, pos_range

# Sequences that begin at the adapter junction = zero-length insert (adapter dimer).
# FastQC reports them without the leading 'A' of AGATCGGAAGAGC (the A-tail base),
# so its Adapter Content module, which searches for AGATCGGAAGAG, largely misses them.
DIMER_PREFIXES = ("AGATCGGAAGAGC", "GATCGGAAGAGC", "CTGTCTCTTATACACATCT",
                  "TGTCTCTTATACACATCT")
TRUSEQ_I7_ANCHOR = "CCAGTCAC"  # index read follows this in the R1 TruSeq adapter
ADAPTER_SOURCE_RE = re.compile(r"adapter|primer|nextera|transposase", re.I)


def _f(x: str) -> float:
    return float(x)


def _sustained_drop(rows, col: int, thr: float) -> int | None:
    """First position from which the value stays below thr in >=80% of remaining bins."""
    vals = [(pos_range(r[0])[0], _f(r[col])) for r in rows]
    for i, (pos, v) in enumerate(vals):
        if v < thr:
            rest = [w for _, w in vals[i:]]
            if sum(w < thr for w in rest) >= 0.8 * len(rest):
                return pos
    return None


def compute(rep: FastQCReport) -> dict:
    m: dict = {}
    b = rep.basic
    m["fastqc_version"] = rep.version
    m["filename"] = rep.filename
    m["encoding"] = b.get("Encoding", "")
    m["phred64"] = "1.5" in m["encoding"] or "1.3" in m["encoding"]
    m["n_reads"] = int(float(b.get("Total Sequences", 0)))
    lens = b.get("Sequence length", "0")
    m["read_len"] = pos_range(lens)[1]
    m["read_len_min"] = pos_range(lens)[0]
    try:
        m["gc"] = float(b["%GC"])
    except (KeyError, ValueError):
        m["gc"] = None  # NaN would break JSON output and every GC rule
    m["status"] = {k: v.status for k, v in rep.modules.items()}

    # --- per-base quality
    pbq = rep.modules.get("Per base sequence quality")
    if pbq and pbq.rows:
        h = pbq.header
        lq, p10, med = h.index("Lower Quartile"), h.index("10th Percentile"), h.index("Median")
        m["q_drop_lq20"] = _sustained_drop(pbq.rows, lq, 20)
        m["q_drop_p10_20"] = _sustained_drop(pbq.rows, p10, 20)
        m["q_drop_median20"] = _sustained_drop(pbq.rows, med, 20)
        m["q_median_last"] = _f(pbq.rows[-1][med])
        m["q_mean_min"] = min(_f(r[1]) for r in pbq.rows)
        m["q_lq_first10_min"] = min(_f(r[lq]) for r in pbq.rows if pos_range(r[0])[0] <= 10)
    # --- per-sequence quality
    psq = rep.modules.get("Per sequence quality scores")
    if psq and psq.rows:
        tot = sum(_f(r[1]) for r in psq.rows) or 1
        m["pct_reads_meanq_lt20"] = 100 * sum(_f(r[1]) for r in psq.rows if _f(r[0]) < 20) / tot
        m["meanq_mode"] = int(max(psq.rows, key=lambda r: _f(r[1]))[0])
    # --- per-tile
    pt = rep.modules.get("Per tile sequence quality")
    if pt and pt.rows:
        worst: dict[str, float] = {}
        for t, _, dev in pt.rows:
            worst[t] = min(worst.get(t, 0.0), _f(dev))
        m["n_tiles"] = len(worst)
        m["pct_tiles_dev5"] = 100 * sum(v <= -5 for v in worst.values()) / len(worst)
        m["pct_tiles_dev10"] = 100 * sum(v <= -10 for v in worst.values()) / len(worst)
        m["tile_worst_dev"] = min(worst.values())
    # --- per-base content
    pbc = rep.modules.get("Per base sequence content")
    if pbc and pbc.rows:
        h = pbc.header
        ig, ia, it, ic = (h.index(x) for x in ("G", "A", "T", "C"))

        def bias(r):
            return max(abs(_f(r[ia]) - _f(r[it])), abs(_f(r[ig]) - _f(r[ic])))
        # Classify FastQC bins by their first position. "End" is only a bin that starts in the
        # last two cycles (e.g. "150-151" of a 2x151 run, "150" of 2x150); a 5-wide bin such as
        # "145-149" belongs to the read body, so real bias there is not mistaken for a
        # last-cycle artefact. "Start" includes a bin that begins within the first 12 bases.
        end_from = m["read_len"] - 1
        start: list[float] = []
        body: list[float] = []
        end: list[float] = []
        for r in pbc.rows:
            first = pos_range(r[0])[0]
            (end if first >= end_from and first > 12 else start if first <= 12 else body).append(bias(r))
        m["comp_bias_start"] = max(start) if start else None
        m["comp_bias_body"] = max(body) if body else None
        m["comp_bias_end"] = max(end) if end else None
    # --- per-sequence GC
    gcm = rep.modules.get("Per sequence GC content")
    if gcm and gcm.rows:
        ys = [_f(r[1]) for r in gcm.rows]
        sm = [sum(ys[max(0, i - 2):i + 3]) / len(ys[max(0, i - 2):i + 3]) for i in range(len(ys))]
        peak = max(sm) or 1
        modes = [i for i in range(1, len(sm) - 1)
                 if sm[i] >= sm[i - 1] and sm[i] > sm[i + 1] and sm[i] >= 0.10 * peak]
        # merge modes closer than 5 %GC
        merged: list[int] = []
        for i in modes:
            if merged and i - merged[-1] < 5:
                if sm[i] > sm[merged[-1]]:
                    merged[-1] = i
            else:
                merged.append(i)
        m["gc_modes"] = len(merged)
        m["gc_mode_positions"] = [int(gcm.rows[i][0]) for i in merged]
        m["gc_peak"] = int(gcm.rows[max(range(len(ys)), key=ys.__getitem__)][0])
    # --- N content
    nc = rep.modules.get("Per base N content")
    if nc and nc.rows:
        m["n_max_pct"] = max(_f(r[1]) for r in nc.rows)
    # --- length distribution
    sld = rep.modules.get("Sequence Length Distribution")
    if sld and sld.rows:
        m["length_uniform"] = len([r for r in sld.rows if _f(r[1]) > 0]) == 1
        tot = sum(_f(r[1]) for r in sld.rows)
        if tot:
            m["mean_len"] = sum(sum(pos_range(r[0])) / 2 * _f(r[1]) for r in sld.rows) / tot
    # --- duplication
    dup = rep.modules.get("Sequence Duplication Levels")
    if dup:
        v = dup.extra.get("Total Deduplicated Percentage")
        m["dedup_pct"] = float(v) if v else None
        # FastQC 0.11/0.12 write "Duplication Level, Percentage of total"; some layouts add
        # "Percentage of deduplicated" first. Find the column by name, not position.
        col = dup.header.index("Percentage of total") if "Percentage of total" in dup.header else 1
        hi = [r for r in dup.rows if r[0].startswith(">") and r[0] not in (">10", ">50")]
        m["pct_seq_dup_gt100"] = sum(_f(r[col]) for r in hi) if dup.rows else None
    # --- overrepresented
    ov = rep.modules.get("Overrepresented sequences")
    seqs: list[dict[str, Any]] = []
    if ov:
        for r in ov.rows:
            seqs.append({"seq": r[0], "count": int(float(r[1])), "pct": _f(r[2]),
                         "source": r[3] if len(r) > 3 else ""})
    m["overrep"] = seqs
    m["dimer_pct"] = sum(s["pct"] for s in seqs if s["seq"].startswith(DIMER_PREFIXES))
    # each overrepresented sequence is counted in exactly one class: poly-G, adapter, other
    def _polyg(s):
        return bool(s["seq"]) and s["seq"].count("G") >= 0.9 * len(s["seq"])

    def _adapter(s):
        return bool(ADAPTER_SOURCE_RE.search(s["source"])) or s["seq"].startswith(DIMER_PREFIXES)
    m["overrep_polyg_pct"] = sum(s["pct"] for s in seqs if _polyg(s))
    m["overrep_adapter_pct"] = sum(s["pct"] for s in seqs if not _polyg(s) and _adapter(s))
    m["overrep_other_pct"] = sum(s["pct"] for s in seqs if not _polyg(s) and not _adapter(s))
    idx: Counter[str] = Counter()
    for s in seqs:
        k = s["seq"].find(TRUSEQ_I7_ANCHOR)
        if s["seq"].startswith(DIMER_PREFIXES[:2]) and k >= 0 and len(s["seq"]) >= k + 14:
            idx[s["seq"][k + 8:k + 14]] += s["count"]
    m["observed_i7"] = idx.most_common(1)[0][0] if idx else None
    srcs: Counter[str] = Counter()
    for s in seqs:
        mm = re.search(r"(TruSeq Adapter, Index \d+)", s["source"])
        if mm:
            srcs[mm.group(1)] += s["count"]
    m["observed_index_name"] = srcs.most_common(1)[0][0] if srcs else None
    # --- adapter content
    ac = rep.modules.get("Adapter Content")
    if ac and ac.rows:
        cols = ac.header[1:]
        last = ac.rows[-1]
        m["adapter_last"] = {c: _f(last[i + 1]) for i, c in enumerate(cols)}
        seq_cols = [c for c in cols if c not in ("PolyA", "PolyG")]
        m["adapter_max_pct"] = max((m["adapter_last"][c] for c in seq_cols), default=0.0)
        m["adapter_top"] = max(seq_cols, key=lambda c: m["adapter_last"][c]) if seq_cols else None

        def onset(thr):
            for r in ac.rows:
                if max(_f(r[i + 1]) for i, c in enumerate(cols) if c in seq_cols) >= thr:
                    return pos_range(r[0])[0]
            return None
        m["adapter_onset_1pct"] = onset(1.0)
        m["adapter_onset_5pct"] = onset(5.0)
        m["polya_pct"] = m["adapter_last"].get("PolyA")
        m["polyg_pct"] = m["adapter_last"].get("PolyG")
    return m

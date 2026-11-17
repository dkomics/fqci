"""Markdown and JSON reporting."""
from __future__ import annotations

import json

from .interpret import SampleResult, natural_key
from .rules import Profile

SEV_LABEL = {"fail": "FAIL", "warn": "WARN", "info": "info"}


def _md(text: str) -> str:
    """Escape characters that break Markdown tables."""
    return str(text).replace("|", "\\|")


def _fmt(v, nd=1):
    if v is None:
        return "n/a"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def _n(k: int, word: str = "sample") -> str:
    return f"{k} {word}" + ("" if k == 1 else "s")


def run_summary(results: list[SampleResult]) -> list[str]:
    """Run-level view: verdict counts, findings shared by most samples, and subsets."""
    real = [r for r in results if not r.undetermined]
    n = len(real)
    if not n:
        return []
    L = ["## Run summary", ""]
    counts: dict[str, int] = {}
    for r in real:
        k = r.verdict.split(";")[0].split(":")[0]
        counts[k] = counts.get(k, 0) + 1
    unassigned = len(results) - n
    L.append(f"{n} samples" + (f" (plus {unassigned} unassigned-read set)" if unassigned else "") + ": "
             + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + ".")
    by_rule: dict[str, tuple[str, set]] = {}
    for r in real:
        for f in r.findings:
            if f.expected or f.severity == "info" or f.secondary_to:
                continue
            by_rule.setdefault(f.rule_id, (f.title, set()))[1].add(r.sample)
    wide = {k: v for k, v in by_rule.items() if len(v[1]) >= max(3, 0.5 * n)}
    some = {k: v for k, v in by_rule.items() if k not in wide}
    if wide:
        L += ["", "Found in most samples, so probably a property of the run or library prep rather than of "
              "individual samples:", ""]
        L += [f"- {t}: {len(s)} of {n}" for k, (t, s) in sorted(wide.items(), key=lambda kv: -len(kv[1][1]))]
    if some:
        L += ["", "Found in a subset of samples:", ""]
        for _k, (t, s) in sorted(some.items(), key=lambda kv: -len(kv[1][1])):
            names = sorted(s, key=natural_key)
            shown = ", ".join(names) if len(names) <= 12 else ", ".join(names[:12]) + f" and {len(names) - 12} more"
            L.append(f"- {t}: {len(s)} of {n} ({shown})")
    dimers = [(r.metrics["pair"].get("pair_dimer_pct") or 0, r.sample) for r in real]
    hi = [d for d in dimers if d[0] >= 10]
    lo = [d for d in dimers if d[0] < 1]
    mid = sorted((d for d in dimers if 1 <= d[0] < 10), key=lambda d: natural_key(d[1]))
    if hi and lo:
        a, b = min(d for d, _ in hi), max(d for d, _ in hi)
        rng = f"{a:.0f}%" if f"{a:.0f}" == f"{b:.0f}" else f"{a:.0f}-{b:.0f}%"
        msg = f"Adapter dimers split the run in two: {_n(len(hi))} at {rng} and {len(lo)} below 1%"
        if mid:
            shown = ", ".join(f"{name} {d:.1f}%" for d, name in mid[:5])
            if len(mid) > 5:
                shown += f" and {len(mid) - 5} more"
            msg += f"; {len(mid)} in between ({shown})"
        L += ["", msg + ". That pattern usually means a library-prep batch, not individual samples."]
    cov = sorted((r.metrics["pair"]["est_coverage"], r.sample) for r in real
                 if r.metrics["pair"].get("est_coverage") is not None)
    if cov:
        L += ["", "Lowest expected depth: " + ", ".join(f"{s} {c:.0f}x" for c, s in cov[:5]) + "."]
    return L + [""]


def to_markdown(results: list[SampleResult], profile: Profile, title: str = "FastQC interpretation") -> str:
    L = [f"# {title}", "",
         f"Profile: **{profile.name}**"
         + (f", purpose: **{profile.purpose}**" if profile.purpose else "")
         + f". {profile.description}", ""]
    L += run_summary(results)
    L += ["## Summary", "",
          "| Sample | Reads/pairs | Dimers % | R1 Q<20 from | R2 Q<20 from | R2 reads meanQ<20 % | Dedup % (R1) | Est. depth | Verdict |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in results:
        reads = [k for k in r.metrics if k != "pair"]
        m1 = r.metrics.get("R1") or (r.metrics[reads[0]] if reads else {})
        m2, pm = r.metrics.get("R2", {}), r.metrics["pair"]
        cov = pm.get("est_coverage")
        L.append("| {} | {:,} | {} | {} | {} | {} | {} | {} | {} |".format(
            _md(r.sample), m1.get("n_reads", 0), _fmt(pm.get("pair_dimer_pct")),
            _fmt(m1.get("q_drop_lq20")), _fmt(m2.get("q_drop_lq20")),
            _fmt(m2.get("pct_reads_meanq_lt20")), _fmt(m1.get("dedup_pct")),
            f"{cov:.0f}x" if cov is not None else "n/a", r.verdict))
    L += ["", "Q<20 from: the position from which the lower quartile stays below Q20. "
          "Est. depth: expected median depth after mapping (read-count estimate x depth factor). "
          "Verdict: the worst primary finding, and whether trimming resolves it.", ""]

    for r in results:
        L += [f"## {r.sample}", "", "Files: " + ", ".join(f"{k} `{v}`" for k, v in r.files.items()), ""]
        if r.undetermined:
            m1 = r.metrics.get("R1", {})
            L += ["These are reads the demultiplexer could not assign to a sample (index errors, "
                  "PhiX spike-in, index hopping). They are not a sample, so no verdict or trimming "
                  "is given. The findings are listed for information only.", "",
                  f"- {m1.get('n_reads', 0):,} reads; GC peaks at {m1.get('gc_mode_positions', 'n/a')}% "
                  "(a PhiX spike-in peaks at about 44%)."]
            L += [f"- {f.title} ({f.read or 'pair'}): {f.explain}" for f in r.findings
                  if f.severity in ("fail", "warn")] + [""]
            continue
        primary = r.primary()
        notes = [f for f in r.findings if not f.expected and not f.secondary_to and f.severity == "info"]
        secondary = [f for f in r.findings if not f.expected and f.secondary_to]
        expected = [f for f in r.findings if f.expected]
        if primary:
            L += ["### Primary findings", ""]
            for f in primary:
                where = f" ({f.read})" if f.read else ""
                L.append(f"**[{SEV_LABEL[f.severity]}] {f.title}{where}.** {f.explain}")
                L.append("")
                if f.cause:
                    L.append(f"Likely cause: {f.cause}")
                    L.append("")
                if f.advice:
                    L.append(f"Advice: {f.advice}")
                    L.append("")
        if notes:
            L += ["### Notes", ""] + [f"- {f.title}: {f.explain}" for f in notes] + [""]
        if secondary:
            L += ["### Probably secondary (re-check after trimming)", "",
                  "These findings fit the primary problem above and are likely to resolve once it is "
                  "trimmed, so they do not count toward the verdict.", ""]
            for f in secondary:
                L.append(f"- [{SEV_LABEL[f.severity]}] {f.title} ({f.read or 'pair'}): {f.explain} "
                         f"*Explained by {', '.join(f.secondary_to)}.*")
            L.append("")
        if expected:
            L += ["### Expected for this library type", ""]
            for f in expected:
                L.append(f"- {f.title} ({f.read or 'pair'}): {f.expected_note}")
            L.append("")
        L += ["### Recommended trimming", "", "```bash", r.fastp_command.replace(" --", " \\\n  --"), "```", "",
              "| Parameter | Value | Why |", "|---|---|---|"]
        for k, v in r.fastp_params.items():
            L.append(f"| `--{k}` | {'on' if v is True else v} | {', '.join(r.param_reasons.get(k, []))} |")
        L.append("")
    L += ["## Next step", "",
          "Run FastQC on the trimmed reads, then `fqci compare --before <raw reports> --after "
          "<trimmed reports>`. Secondary findings that persist after trimming are real "
          "problems in their own right.", ""]
    if profile.notes:
        L += ["## Profile notes", ""] + [f"- {n}" for n in profile.notes] + [""]
    return "\n".join(L)


def to_json(results: list[SampleResult], profile: Profile) -> str:
    return json.dumps({"profile": profile.name, "samples": [r.to_dict() for r in results]},
                      indent=2, default=str)

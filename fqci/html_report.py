"""Self-contained HTML report: one file, no external resources, works offline and by email.

Layout: header, run overview (verdict tiles, run summary, depth and dimer charts), a sortable
sample table, and one collapsible section per sample with findings and the fastp command.
Status colours always come with an icon and a word, never colour alone.
"""
from __future__ import annotations

import datetime as _dt
import html
import json

from . import __version__
from typing import Any

from .interpret import SampleResult
from .report import run_summary
from .rules import Profile

DASH = "\u2013"
STATUS = {  # verdict level -> (css class, icon, label)
    "FAIL": ("critical", "✕", "FAIL"),
    "WARN": ("warning", "!", "WARN"),
    "PASS": ("good", "✓", "PASS"),
    "INFO": ("neutral", "i", "info"),
    "N/A": ("neutral", "–", "n/a"),
}

CSS = """
.fqci{color-scheme:light;--surface-0:#f7f7f5;--surface-1:#fcfcfb;--border:#e2e1dc;--grid:#ecebe7;
--text-primary:#0b0b0b;--text-secondary:#52514e;--text-muted:#77766f;--series-1:#2a78d6;
--good:#0ca30c;--warning:#fab219;--critical:#d03b3b;--neutral:#8a8983;--code-bg:#f0efec;
font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;background:var(--surface-0);
color:var(--text-primary);margin:0}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])) .fqci{color-scheme:dark;
--surface-0:#121211;--surface-1:#1a1a19;--border:#34332f;--grid:#2a2a27;--text-primary:#f0efec;
--text-secondary:#c3c2b7;--text-muted:#9a998f;--series-1:#3987e5;--code-bg:#24241f}}
.fqci *{box-sizing:border-box}
.fqci main{max-width:1180px;margin:0 auto;padding:28px 24px 64px}
.fqci h1{font-size:24px;margin:0 0 4px}
.fqci h2{font-size:18px;margin:32px 0 12px}
.fqci h3{font-size:15px;margin:18px 0 8px}
.fqci .meta{color:var(--text-secondary);font-size:13px}
.fqci .card{background:var(--surface-1);border:1px solid var(--border);border-radius:10px;padding:16px 18px;margin:12px 0}
.fqci .tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:16px 0}
.fqci .tile{background:var(--surface-1);border:1px solid var(--border);border-radius:10px;padding:12px 14px}
.fqci .tile .n{font-size:28px;font-weight:650;line-height:1.1}
.fqci .tile .l{color:var(--text-secondary);font-size:13px}
.fqci .pill{display:inline-flex;align-items:center;gap:6px;font-size:12px;font-weight:600;
white-space:nowrap;padding:1px 9px 1px 3px;border-radius:999px;border:1px solid var(--border);background:var(--surface-1)}
.fqci .pill i{font-style:normal;display:inline-grid;place-items:center;width:17px;height:17px;border-radius:50%;
color:#fff;font-size:11px;font-weight:700}
.fqci .critical i{background:var(--critical)} .fqci .warning i{background:var(--warning);color:#0b0b0b}
.fqci .good i{background:var(--good)} .fqci .neutral i{background:var(--neutral)}
.fqci table{border-collapse:collapse;width:100%;font-size:13.5px}
.fqci th,.fqci td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--grid);vertical-align:top}
.fqci th{color:var(--text-secondary);font-weight:600;cursor:pointer;user-select:none;white-space:nowrap}
.fqci th[aria-sort]::after{content:" \\2195";color:var(--text-muted)}
.fqci td.num,.fqci th.num{text-align:right;font-variant-numeric:tabular-nums}
.fqci .scroll{overflow-x:auto}
.fqci details.sample{background:var(--surface-1);border:1px solid var(--border);border-radius:10px;margin:10px 0}
.fqci details.sample>summary{list-style:none;cursor:pointer;padding:12px 16px;display:flex;gap:12px;align-items:center;flex-wrap:wrap}
.fqci details.sample>summary::-webkit-details-marker{display:none}
.fqci details.sample>summary::before{content:"\\25B8";color:var(--text-muted)}
.fqci details.sample[open]>summary::before{content:"\\25BE"}
.fqci details.sample>summary .name{font-weight:650}
.fqci details.sample>summary .why{color:var(--text-secondary);font-size:13px}
.fqci .body{padding:0 18px 16px}
.fqci .finding{border-left:3px solid var(--border);padding:6px 0 6px 12px;margin:10px 0}
.fqci .finding.critical{border-color:var(--critical)} .fqci .finding.warning{border-color:var(--warning)}
.fqci .finding .t{font-weight:600}
.fqci .finding p{margin:4px 0;color:var(--text-secondary)}
.fqci .finding p b{color:var(--text-primary);font-weight:600}
.fqci ul.plain{margin:6px 0;padding-left:20px;color:var(--text-secondary)}
.fqci pre{background:var(--code-bg);border-radius:8px;padding:12px 14px;overflow-x:auto;font-size:12.5px;margin:8px 0}
.fqci code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.fqci button.copy{font:inherit;font-size:12px;padding:3px 10px;border-radius:6px;border:1px solid var(--border);
background:var(--surface-1);color:var(--text-primary);cursor:pointer}
.fqci .charts{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:12px}
.fqci .chart h3{margin:0 0 2px} .fqci .chart .sub{color:var(--text-secondary);font-size:12.5px;margin-bottom:8px}
.fqci .chart svg{display:block;max-width:620px}
.fqci svg text{fill:var(--text-secondary);font-size:11px}
.fqci svg text.val{paint-order:stroke;stroke:var(--surface-1);stroke-width:3px;stroke-linejoin:round}
.fqci a{color:var(--series-1)}
.fqci svg .bar{fill:var(--series-1)} .fqci svg .bar:hover{opacity:.8}
.fqci svg .ref{stroke:var(--text-muted);stroke-dasharray:3 3;stroke-width:1}
.fqci svg .axis{stroke:var(--grid)}
.fqci #tip{position:fixed;pointer-events:none;background:var(--surface-1);color:var(--text-primary);
border:1px solid var(--border);border-radius:6px;padding:4px 8px;font-size:12px;display:none;box-shadow:0 2px 8px rgba(0,0,0,.12)}
.fqci .muted{color:var(--text-muted)}
.fqci .ai{border:1px dashed var(--series-1);border-radius:8px;padding:10px 14px;margin:10px 0;background:var(--surface-0)}
.fqci .ai .lab{display:inline-block;font-size:11px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;
color:var(--series-1);margin-bottom:4px}
.fqci .ai p{margin:4px 0}
.fqci .ai .chk{color:var(--text-secondary);font-size:12.5px}
.fqci .ai .rec{color:var(--text-muted);font-size:12px}
.fqci .toolbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:4px 0 10px}
.fqci .toolbar button{font:inherit;font-size:12.5px;cursor:pointer;border:1px solid var(--border);
background:var(--surface-1);color:var(--text-primary);border-radius:999px;padding:3px 10px}
.fqci .toolbar button.chip{padding:2px 6px}
.fqci .toolbar button.chip .pill{border:0;padding-right:4px;background:transparent}
.fqci .toolbar button[aria-pressed="true"]{border-color:var(--series-1);box-shadow:0 0 0 1px var(--series-1) inset}
.fqci .toolbar .sep{width:1px;height:20px;background:var(--border);margin:0 4px}
.fqci details.sample[hidden]{display:none}
@media print{.fqci{background:#fff}.fqci button.copy,.fqci #tip{display:none}.fqci details.sample{break-inside:avoid}}
"""

JS = """
(function(){
  var tip=document.getElementById('tip');
  document.querySelectorAll('[data-tip]').forEach(function(el){
    el.addEventListener('mousemove',function(e){tip.textContent=el.getAttribute('data-tip');
      tip.style.display='block';tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY+12)+'px';});
    el.addEventListener('mouseleave',function(){tip.style.display='none';});
  });
  document.querySelectorAll('button.copy').forEach(function(b){
    b.addEventListener('click',function(){var t=document.getElementById(b.getAttribute('data-for')).innerText;
      if(navigator.clipboard){navigator.clipboard.writeText(t).then(function(){b.textContent='Copied';
        setTimeout(function(){b.textContent='Copy';},1500);});}});
  });
  document.querySelectorAll('table.sortable th').forEach(function(th,i){
    th.addEventListener('click',function(){var tb=th.closest('table').tBodies[0];
      var rows=Array.prototype.slice.call(tb.rows);var asc=th.getAttribute('data-dir')!=='asc';
      th.parentNode.querySelectorAll('th').forEach(function(x){x.removeAttribute('data-dir');});
      th.setAttribute('data-dir',asc?'asc':'desc');
      rows.sort(function(a,b){var x=a.cells[i].getAttribute('data-v')||a.cells[i].innerText,
        y=b.cells[i].getAttribute('data-v')||b.cells[i].innerText;var nx=parseFloat(x),ny=parseFloat(y);
        var c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y,undefined,{numeric:true});return asc?c:-c;});
      rows.forEach(function(r){tb.appendChild(r);});});
  });
  var secs=Array.prototype.slice.call(document.querySelectorAll('details.sample'));
  var shown=document.getElementById('shown');
  function count(){var n=secs.filter(function(d){return !d.hidden;}).length;
    if(shown){shown.textContent=n+' of '+secs.length+' shown';}}
  document.querySelectorAll('.toolbar [data-act]').forEach(function(b){
    b.addEventListener('click',function(){var open=b.getAttribute('data-act')==='expand';
      secs.forEach(function(d){if(!d.hidden){d.open=open;}});});
  });
  document.querySelectorAll('.toolbar .chip').forEach(function(c){
    c.addEventListener('click',function(){var f=c.getAttribute('data-filter');
      document.querySelectorAll('.toolbar .chip').forEach(function(x){x.setAttribute('aria-pressed',x===c?'true':'false');});
      secs.forEach(function(d){d.hidden=!(f==='ALL'||d.getAttribute('data-level')===f);});count();});
  });
  document.querySelectorAll('a[href^="#s"]').forEach(function(a){
    a.addEventListener('click',function(){var d=document.querySelector(a.getAttribute('href'));
      if(d){d.hidden=false;d.open=true;count();}});
  });
  count();
  window.addEventListener('beforeprint',function(){document.querySelectorAll('details.sample').forEach(function(d){d.open=true;});});
})();
"""


def _e(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def _level(verdict: str) -> str:
    v = verdict.upper()
    if v.startswith("N/A"):
        return "N/A"
    for k in ("FAIL", "WARN", "PASS"):
        if v.startswith(k):
            return k
    return "INFO"


def _pill(level: str, text: str | None = None) -> str:
    cls, icon, label = STATUS[level]
    return f'<span class="pill {cls}"><i aria-hidden="true">{icon}</i>{_e(text or label)}</span>'


def _num(v, nd=1, suffix="") -> str:
    if v is None:
        return '<span class="muted">n/a</span>'
    if isinstance(v, float):
        return f"{v:,.{nd}f}{suffix}"
    if isinstance(v, int):
        return f"{v:,}{suffix}"
    return _e(v)


def _bar_chart(rows: list[tuple[str, float]], unit: str, refs: list[tuple[float, str]],
               label_fmt: str = "{:.0f}") -> str:
    """Horizontal single-series bar chart with optional dashed reference lines.
    One hue, no legend (the title names the series), labels on bars, hover tooltips."""
    if not rows:
        return '<p class="muted">No data.</p>'
    names_w = min(170, 7 * max(len(n) for n, _ in rows) + 12)
    w, bar_h, gap, top = 560, 16, 6, 18
    plot_w = w - names_w - 54
    vmax = max([v for _, v in rows] + [r for r, _ in refs] + [1e-9]) * 1.08
    h = top + len(rows) * (bar_h + gap) + 4
    x = lambda v: names_w + plot_w * v / vmax  # noqa: E731 - local scale
    out = [f'<svg viewBox="0 0 {w} {h}" width="100%" role="img" aria-label="{_e(unit)} per sample">']
    for r, lab in refs:
        out.append(f'<line class="ref" x1="{x(r):.1f}" x2="{x(r):.1f}" y1="{top - 6}" y2="{h - 2}"/>'
                   f'<text x="{x(r) + 3:.1f}" y="{top - 8}" font-size="10">{_e(lab)}</text>')
    out.append(f'<line class="axis" x1="{names_w}" x2="{names_w}" y1="{top - 4}" y2="{h - 2}"/>')
    for i, (name, v) in enumerate(rows):
        y = top + i * (bar_h + gap)
        bw = max(0.0, x(v) - names_w)
        short = name if len(name) <= 24 else name[:23] + "…"
        out.append(f'<text x="{names_w - 6}" y="{y + bar_h - 4}" text-anchor="end">{_e(short)}</text>')
        out.append(f'<path class="bar" data-tip="{_e(name)}: {_e(label_fmt.format(v))} {_e(unit)}" '
                   f'd="M{names_w},{y} h{max(bw - 4, 0):.1f} q4,0 4,4 v{bar_h - 8} q0,4 -4,4 h{-max(bw - 4, 0):.1f} z"/>')
        out.append(f'<text class="val" x="{names_w + bw + 4:.1f}" y="{y + bar_h - 4}">{_e(label_fmt.format(v))}</text>')
    out.append("</svg>")
    return "".join(out)


def _finding_html(f, cls: str) -> str:
    where = f" ({_e(f.read)})" if f.read else ""
    parts = [f'<div class="finding {cls}"><div class="t">{_pill(f.severity.upper() if f.severity != "info" else "INFO")} '
             f'{_e(f.title)}{where}</div><p>{_e(f.explain)}</p>']
    if f.cause:
        parts.append(f"<p><b>Likely cause:</b> {_e(f.cause)}</p>")
    if f.advice:
        parts.append(f"<p><b>Advice:</b> {_e(f.advice)}</p>")
    parts.append("</div>")
    return "".join(parts)


def _ai_box(s: Any) -> str:
    """Per-sample AI box; s is an ai_summary.SampleSummary."""
    parts = ['<div class="ai"><span class="lab">AI summary \u00b7 generated, advisory</span>',
             f"<p>{_e(s.summary)}</p>"]
    if s.next_steps:
        parts.append('<ul class="plain">' + "".join(f"<li>{_e(x)}</li>" for x in s.next_steps) + "</ul>")
    if s.rule_ids:
        parts.append(f'<p class="chk">Based on: {_e(", ".join(s.rule_ids))}</p>')
    parts += [f'<p class="chk">Check: {_e(w)}</p>' for w in s.warnings]
    parts.append("</div>")
    return "".join(parts)


def _sample_section(r: SampleResult, idx: int, ai_sample: Any = None) -> str:
    level = _level(r.verdict)
    prim = r.primary()
    first = prim[0].title if prim else ("unassigned reads" if r.undetermined else "no problems found")
    head = (f'<summary>{_pill(level)}<span class="name">{_e(r.sample)}</span>'
            f'<span class="why">{_e(r.verdict)} · {_e(first)}</span></summary>')
    body = ['<div class="body">', '<p class="meta">Files: '
            + ", ".join(f"{_e(k)} <code>{_e(v)}</code>" for k, v in r.files.items()) + "</p>"]
    if ai_sample is not None:
        body.append(_ai_box(ai_sample))
    if r.undetermined:
        body.append("<p>Reads the demultiplexer could not assign to a sample. Not a sample: no verdict "
                    "or trimming is given; findings are listed for information.</p>")
    if prim:
        body.append("<h3>Primary findings</h3>")
        body += [_finding_html(f, "critical" if f.severity == "fail" else "warning") for f in prim]
    notes = [f for f in r.findings if not f.expected and not f.secondary_to and f.severity == "info"]
    if notes:
        body.append("<h3>Notes</h3><ul class=\"plain\">"
                    + "".join(f"<li><b>{_e(f.title)}:</b> {_e(f.explain)}</li>" for f in notes) + "</ul>")
    sec = [f for f in r.findings if not f.expected and f.secondary_to]
    if sec:
        body.append("<h3>Probably secondary (re-check after trimming)</h3><ul class=\"plain\">"
                    + "".join(f"<li>{_e(f.title)} ({_e(f.read or 'pair')}): {_e(f.explain)} "
                              f"<i>Explained by {_e(', '.join(f.secondary_to))}.</i></li>" for f in sec) + "</ul>")
    exp = [f for f in r.findings if f.expected]
    if exp:
        body.append("<h3>Expected for this library type</h3><ul class=\"plain\">"
                    + "".join(f"<li>{_e(f.title)} ({_e(f.read or 'pair')}): {_e(f.expected_note)}</li>"
                              for f in exp) + "</ul>")
    if not r.undetermined:
        cid = f"cmd{idx}"
        cmd = r.fastp_command
        shown = cmd if cmd.startswith("#") else cmd.replace(" --", " \\\n  --")
        body.append(f'<h3>Recommended trimming <button class="copy" data-for="{cid}">Copy</button></h3>'
                    f'<pre><code id="{cid}">{_e(shown)}</code></pre>')
        if r.fastp_params:
            body.append('<table><thead><tr><th>Parameter</th><th>Value</th><th>Why</th></tr></thead><tbody>'
                        + "".join(f"<tr><td><code>--{_e(k)}</code></td><td>{'on' if v is True else _e(v)}</td>"
                                  f"<td>{_e(', '.join(r.param_reasons.get(k, [])))}</td></tr>"
                                  for k, v in r.fastp_params.items()) + "</tbody></table>")
    body.append("</div>")
    return f'<details class="sample" id="s{idx}" data-level="{level}">{head}{"".join(body)}</details>'


def to_html(results: list[SampleResult], profile: Profile, title: str = "FastQC interpretation",
            ai: Any = None, ai_error: str | None = None) -> str:
    """ai: an ai_summary.AISummary (optional); ai_error: why a requested AI summary is missing."""
    real = [r for r in results if not r.undetermined]
    counts = {"FAIL": 0, "WARN": 0, "PASS": 0}
    for r in real:
        lv = _level(r.verdict)
        if lv in counts:
            counts[lv] += 1
    fixable = sum(1 for r in real if r.verdict.endswith(": fixable by trimming"))
    gen = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    purpose = f" · purpose <b>{_e(profile.purpose)}</b>" if profile.purpose else ""

    parts = [f'<header><h1>{_e(title)}</h1><div class="meta">Profile <b>{_e(profile.name)}</b>{purpose} '
             f'· {len(real)} samples · generated {gen} by fqci {__version__}</div>'
             f'<p class="meta">{_e(profile.description)}</p></header>']

    # tiles
    parts.append('<div class="tiles">'
                 f'<div class="tile"><div class="n">{len(real)}</div><div class="l">samples</div></div>'
                 + "".join(f'<div class="tile"><div class="n">{n}</div><div class="l">{_pill(k)}</div></div>'
                           for k, n in counts.items())
                 + f'<div class="tile"><div class="n">{fixable}</div><div class="l">fixable by trimming</div></div>'
                 + "</div>")

    # run summary (reuse the Markdown logic; drop its heading)
    rs = [ln for ln in run_summary(results) if ln and not ln.startswith("## ")]
    if rs:
        items, para = [], []
        for ln in rs:
            if ln.startswith("- "):
                items.append(f"<li>{_e(ln[2:])}</li>")
            else:
                if items:
                    para.append('<ul class="plain">' + "".join(items) + "</ul>")
                    items = []
                para.append(f"<p>{_e(ln)}</p>")
        if items:
            para.append('<ul class="plain">' + "".join(items) + "</ul>")
        parts.append('<section class="card"><h2 style="margin-top:0">Run summary</h2>' + "".join(para) + "</section>")
    if ai is not None:
        rec = ai.record
        box = ['<section class="ai"><span class="lab">AI summary \u00b7 generated, advisory</span>',
               f"<p>{_e(ai.run_summary)}</p>"]
        box += [f'<p class="chk">Check: {_e(w)}</p>' for w in ai.warnings]
        box.append(f'<p class="rec">{_e(rec["provider"])} / {_e(rec["model"])} \u00b7 prompt {_e(rec["prompt"])} '
                   f'({_e(rec["prompt_sha256"])}) \u00b7 {_e(rec["generated"])}'
                   f'{" · sample names anonymised before sending" if rec.get("anonymised") else ""}. '
                   "The rule-based findings and verdicts are authoritative; each sample\u2019s AI summary "
                   "is under Findings by sample.</p></section>")
        parts.append("".join(box))
    elif ai_error:
        parts.append(f'<section class="ai"><span class="lab">AI summary</span><p class="chk">Requested but not '
                     f"generated: {_e(ai_error)}. The rule-based report below is complete.</p></section>")

    # charts
    depth = [(r.sample, r.metrics["pair"]["est_coverage"]) for r in real
             if r.metrics["pair"].get("est_coverage") is not None]
    dimer = [(r.sample, r.metrics["pair"].get("pair_dimer_pct") or 0.0) for r in real]
    charts = []
    if depth:
        p = profile.params
        refs = [(p["min_cov_warn"], f"warn {p['min_cov_warn']}x"), (p["min_cov_fail"], f"fail {p['min_cov_fail']}x")]
        charts.append('<div class="card chart"><h3>Expected median depth</h3><div class="sub">'
                      'Read-count estimate x depth factor; dashed lines are the warn and fail thresholds.</div>'
                      + _bar_chart(depth, "x", refs) + "</div>")
    if any(v > 0 for _, v in dimer):
        charts.append('<div class="card chart"><h3>Adapter dimers (% of read pairs)</h3><div class="sub">'
                      f'Dashed line: {profile.params.get("dimer_fail", 10)}% fail threshold.</div>'
                      + _bar_chart(dimer, "% dimers", [(profile.params.get("dimer_fail", 10), "fail")], "{:.1f}")
                      + "</div>")
    if charts:
        parts.append('<h2>Overview</h2><div class="charts">' + "".join(charts) + "</div>")

    # sample table
    rows = []
    for i, r in enumerate(results):
        reads = [k for k in r.metrics if k != "pair"]
        m1 = r.metrics.get("R1") or (r.metrics[reads[0]] if reads else {})
        m2, pm = r.metrics.get("R2", {}), r.metrics["pair"]
        lv = _level(r.verdict)
        order = {"FAIL": 0, "WARN": 1, "PASS": 2}.get(lv, 3)
        cov = pm.get("est_coverage")
        r2q = _num(m2.get("pct_reads_meanq_lt20")) if m2 else DASH
        rows.append(
            f'<tr><td data-v="{_e(r.sample)}"><a href="#s{i}">{_e(r.sample)}</a></td>'
            f'<td data-v="{order}">{_pill(lv)}</td>'
            f'<td class="num" data-v="{m1.get("n_reads", 0)}">{_num(m1.get("n_reads"))}</td>'
            f'<td class="num" data-v="{pm.get("pair_dimer_pct") or 0}">{_num(pm.get("pair_dimer_pct"))}</td>'
            f'<td class="num" data-v="{m2.get("pct_reads_meanq_lt20", -1) if m2 else -1}">{r2q}</td>'
            f'<td class="num" data-v="{m1.get("dedup_pct") or 0}">{_num(m1.get("dedup_pct"))}</td>'
            f'<td class="num" data-v="{cov if cov is not None else -1}">{_num(cov, 0, "x") if cov is not None else _num(None)}</td>'
            f'<td>{_e(r.verdict)}</td></tr>')
    parts.append('<h2>Samples</h2><div class="card scroll"><table class="sortable"><thead><tr>'
                 '<th aria-sort="none">Sample</th><th aria-sort="none">Verdict</th>'
                 '<th class="num" aria-sort="none">Reads/pairs</th><th class="num" aria-sort="none">Dimers %</th>'
                 '<th class="num" aria-sort="none">R2 meanQ&lt;20 %</th><th class="num" aria-sort="none">Dedup % (R1)</th>'
                 '<th class="num" aria-sort="none">Exp. depth</th><th aria-sort="none">Detail</th></tr></thead><tbody>'
                 + "".join(rows) + "</tbody></table></div>"
                 '<p class="meta">Click a column heading to sort; click a sample to jump to its findings.</p>')

    # per-sample sections
    present = [lv for lv in ("FAIL", "WARN", "PASS", "N/A") if any(_level(r.verdict) == lv for r in results)]
    chips = "".join(f'<button type="button" class="chip" data-filter="{lv}" aria-pressed="false">'
                    f'{_pill(lv)}</button>' for lv in present)
    parts.append('<h2>Findings by sample</h2><div class="toolbar" role="toolbar" aria-label="Findings controls">'
                 '<button type="button" class="tb" data-act="expand">Expand all</button>'
                 '<button type="button" class="tb" data-act="collapse">Collapse all</button>'
                 '<span class="sep" aria-hidden="true"></span><span class="meta">Show:</span>'
                 '<button type="button" class="chip" data-filter="ALL" aria-pressed="true">All</button>'
                 f'{chips}<span class="meta" id="shown"></span></div>'
                 + "".join(_sample_section(r, i, ai.samples.get(r.sample) if ai is not None else None)
                           for i, r in enumerate(results)))
    if profile.notes:
        parts.append('<h2>Profile notes</h2><ul class="plain">'
                     + "".join(f"<li>{_e(n)}</li>" for n in profile.notes) + "</ul>")
    parts.append('<p class="meta">Next step: run FastQC on the trimmed reads, then '
                 '<code>fqci compare --before &lt;raw&gt; --after &lt;trimmed&gt; --fastp-dir .</code></p>')

    data = json.dumps({"fqci_version": __version__, "profile": profile.name,
                       "samples": [{"sample": r.sample, "verdict": r.verdict} for r in results]})
    # encode markup characters so embedded JSON can never close or open an element
    safe_data = data.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    return ("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<title>{_e(title)}</title><style>{CSS}</style></head>"
            f'<body class="fqci"><main>{"".join(parts)}</main><div id="tip" role="tooltip"></div>'
            f'<script type="application/json" id="fqci-data">{safe_data}</script>'
            f"<script>{JS}</script></body></html>\n")

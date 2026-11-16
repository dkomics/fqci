"""Optional AI summary of fqci findings (off by default).

The rule-based findings stay authoritative. This module sends a compact JSON description of
them (metrics, findings, verdicts, fastp options; never reads or file contents) together with a
fixed, versioned prompt to a language model, and checks the reply before it is shown:

- the reply must match a JSON schema (one entry per sample id);
- rule IDs it cites must exist in that sample's findings (others are removed);
- numbers it quotes must appear in the input (others are listed as unverified).

Sample names can be replaced by codes (S01, S02, ...) before anything leaves the machine.
Providers are pluggable; Anthropic is implemented first.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
from dataclasses import asdict, dataclass, field
from importlib import resources
from typing import Any, Protocol

from . import __version__
from .interpret import SampleResult
from .report import run_summary
from .rules import Profile

DEFAULT_PROMPT = "findings_v1"
DEFAULT_MODEL = "claude-sonnet-5-5"

OUTPUT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "run_summary": {"type": "string"},
        "samples": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "summary": {"type": "string"},
                    "rule_ids": {"type": "array", "items": {"type": "string"}},
                    "next_steps": {"type": "array", "items": {"type": "string"}},
                },
                "required": ["id", "summary", "rule_ids", "next_steps"],
                "additionalProperties": False,
            },
        },
    },
    "required": ["run_summary", "samples"],
    "additionalProperties": False,
}


class AIError(RuntimeError):
    """The AI summary could not be produced (missing SDK or key, API error, bad reply)."""


class Provider(Protocol):
    name: str
    model: str

    def generate(self, system: str, user: str, schema: dict, max_tokens: int) -> tuple[str, dict]:
        """Return (JSON text, metadata such as the exact model and token usage)."""
        ...


class AnthropicProvider:
    """Claude via the Anthropic API. Reads ANTHROPIC_API_KEY (and ANTHROPIC_BASE_URL if set)."""

    name = "anthropic"

    def __init__(self, model: str = DEFAULT_MODEL, timeout: float = 120.0):
        self.model = model
        try:
            import anthropic  # optional dependency: pip install "fqci[ai]"
        except ImportError:
            raise AIError("the Anthropic SDK is not installed; install it with: pip install anthropic") from None
        if not os.environ.get("ANTHROPIC_API_KEY"):
            raise AIError("ANTHROPIC_API_KEY is not set")
        self._client = anthropic.Anthropic(timeout=timeout, max_retries=2)
        self._anthropic = anthropic

    def generate(self, system: str, user: str, schema: dict, max_tokens: int) -> tuple[str, dict]:
        try:
            resp = self._client.messages.create(
                model=self.model,
                max_tokens=max_tokens,
                system=system,
                messages=[{"role": "user", "content": user}],
                # structured output; sent via extra_body so older SDK versions also work
                extra_body={"output_config": {"format": {"type": "json_schema", "schema": schema}}},
            )
        except self._anthropic.APIError as e:
            raise AIError(f"Anthropic API error: {e}") from None
        if resp.stop_reason != "end_turn":
            raise AIError(f"the model stopped early (stop_reason={resp.stop_reason}); no summary produced")
        # replies can also hold thinking blocks; take the first text block only
        text = next((getattr(b, "text", None) for b in resp.content if getattr(b, "type", "") == "text"), None)
        if not text:
            raise AIError("the reply contained no text")
        usage = getattr(resp, "usage", None)
        meta = {"model": getattr(resp, "model", self.model),
                "input_tokens": getattr(usage, "input_tokens", None),
                "output_tokens": getattr(usage, "output_tokens", None)}
        return text, meta


PROVIDERS = {"anthropic": AnthropicProvider}


def get_provider(name: str, model: str | None = None) -> Provider:
    if name not in PROVIDERS:
        raise AIError(f"unknown AI provider {name!r} (available: {', '.join(PROVIDERS)})")
    return PROVIDERS[name](model=model or DEFAULT_MODEL)  # type: ignore[return-value]


def load_prompt(path: str | None = None) -> tuple[str, str, str]:
    """Return (prompt text, prompt name, short sha256). A custom file overrides the built-in."""
    if path:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        name = os.path.basename(path)
    else:
        text = resources.files("fqci").joinpath(f"prompts/{DEFAULT_PROMPT}.md").read_text(encoding="utf-8")
        name = DEFAULT_PROMPT
    return text, name, hashlib.sha256(text.encode()).hexdigest()[:12]


def _r(v):
    return round(v, 1) if isinstance(v, float) else v


def _category(f) -> str:
    if f.expected:
        return "expected"
    if f.secondary_to:
        return "secondary"
    return "note" if f.severity == "info" else "primary"


def build_payload(results: list[SampleResult], profile: Profile, genome_size: float | None,
                  expected_gc: float | None, anonymise: bool = False) -> tuple[dict, dict[str, str]]:
    """Compact, model-facing description of the run. Returns (payload, id -> real sample name)."""
    ids: dict[str, str] = {}
    samples = []
    width = max(2, len(str(len(results))))
    for i, r in enumerate(results, 1):
        sid = f"S{i:0{width}d}" if anonymise else r.sample
        ids[sid] = r.sample
        reads = [k for k in r.metrics if k != "pair"]
        m1 = r.metrics.get("R1") or (r.metrics[reads[0]] if reads else {})
        m2, pm = r.metrics.get("R2", {}), r.metrics["pair"]
        metrics = {
            "reads_or_pairs": m1.get("n_reads"), "read_length": m1.get("read_len"),
            "paired": "R2" in r.metrics, "adapter_dimer_pct": _r(pm.get("pair_dimer_pct")),
            "expected_depth_x": _r(pm.get("est_coverage")), "sequence_depth_x": _r(pm.get("est_coverage_raw")),
            "dedup_pct_R1": _r(m1.get("dedup_pct")), "gc_pct_R1": _r(m1.get("gc")),
            "pct_reads_meanQ_below_20_R2": _r(m2.get("pct_reads_meanq_lt20")) if m2 else None,
        }
        findings = [{"rule_id": f.rule_id, "severity": f.severity, "category": _category(f),
                     "read": f.read or "pair", "title": f.title,
                     "explanation": f.expected_note if f.expected else f.explain,
                     "explained_by": f.secondary_to}
                    for f in r.findings]
        samples.append({"id": sid, "verdict": r.verdict, "unassigned_reads": r.undetermined,
                        "metrics": {k: v for k, v in metrics.items() if v is not None},
                        "findings": findings,
                        "fastp_options": {} if r.undetermined else r.fastp_params})
    p = profile.params
    payload = {
        "tool": f"fqci {__version__}", "profile": profile.name, "purpose": profile.purpose,
        "genome_size_bp": genome_size, "expected_gc_pct": expected_gc,
        "thresholds": {"depth_warn_x": p.get("min_cov_warn"), "depth_fail_x": p.get("min_cov_fail"),
                       "adapter_dimer_fail_pct": p.get("dimer_fail")},
        "run_summary_from_rules": [ln for ln in run_summary(results) if ln and not ln.startswith("## ")],
        "samples": samples,
    }
    if anonymise:  # the rule-based run summary names samples; replace them with codes
        text = json.dumps(payload["run_summary_from_rules"])
        for sid, name in sorted(ids.items(), key=lambda kv: -len(kv[1])):
            text = text.replace(name, sid)
        payload["run_summary_from_rules"] = json.loads(text)
    return payload, ids


# a number: comma-grouped thousands (4,052,587) or a plain decimal (42.5); not part of a word (R1, S01, Q20)
_NUM_RE = re.compile(r"(?<![A-Za-z_\d.])-?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?![\d])")


def _num(s: str) -> float | None:
    try:
        return float(s.replace(",", ""))
    except ValueError:
        return None


def _numbers_in(obj) -> set[str]:
    out: set[str] = set()
    if isinstance(obj, dict):
        for v in obj.values():
            out |= _numbers_in(v)
    elif isinstance(obj, list):
        for v in obj:
            out |= _numbers_in(v)
    elif isinstance(obj, bool) or obj is None:
        pass
    elif isinstance(obj, (int, float)):
        for f in (obj, round(obj), round(obj, 1), round(obj, 2)):
            out.add(f"{f:g}")
            out.add(f"{float(f):.1f}")
    elif isinstance(obj, str):
        for m in _NUM_RE.findall(obj):
            v = _num(m)
            if v is not None:
                out |= _numbers_in(v)
    return out


def _quoted_numbers(text: str) -> list[str]:
    return [m for m in _NUM_RE.findall(text)]


def _known(n: str, allowed: set[str]) -> bool:
    v = _num(n)
    return v is None or bool(_numbers_in(v) & allowed)


@dataclass
class SampleSummary:
    id: str
    sample: str
    summary: str
    rule_ids: list[str]
    next_steps: list[str]
    warnings: list[str] = field(default_factory=list)


@dataclass
class AISummary:
    run_summary: str
    samples: dict[str, SampleSummary]
    record: dict
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"record": self.record, "warnings": self.warnings, "run_summary": self.run_summary,
                "samples": {k: asdict(v) for k, v in self.samples.items()}}


def validate(reply_text: str, payload: dict, ids: dict[str, str]) -> tuple[str, dict[str, SampleSummary], list[str]]:
    """Check the model's reply against the input. Raises AIError if it is unusable."""
    try:
        data = json.loads(reply_text)
    except json.JSONDecodeError as e:
        raise AIError(f"the reply is not valid JSON: {e}") from None
    if not isinstance(data, dict) or not isinstance(data.get("samples"), list) \
            or not isinstance(data.get("run_summary"), str):
        raise AIError("the reply does not match the expected structure")
    allowed = _numbers_in(payload) | {str(i) for i in range(0, 11)}
    run_warn: list[str] = []
    rules_by_id = {s["id"]: {f["rule_id"] for f in s["findings"]} for s in payload["samples"]}
    out: dict[str, SampleSummary] = {}
    for item in data["samples"]:
        sid = str(item.get("id", ""))
        if sid not in ids:
            run_warn.append(f"reply mentioned an unknown sample id {sid!r}; ignored")
            continue
        if sid in out:
            run_warn.append(f"reply had two entries for {sid}; kept the first")
            continue
        warns = []
        cited = [str(x) for x in item.get("rule_ids", [])]
        bad = [x for x in cited if x not in rules_by_id[sid]]
        if bad:
            warns.append(f"removed rule IDs not among this sample's findings: {', '.join(bad)}")
        steps = [str(s) for s in item.get("next_steps", [])][:3]
        text = str(item.get("summary", ""))
        unverified = sorted({n for n in _quoted_numbers(" ".join([text, *steps])) if not _known(n, allowed)})
        if unverified:
            warns.append("numbers not found in the input (verify before use): " + ", ".join(unverified))
        out[sid] = SampleSummary(id=sid, sample=ids[sid], summary=text,
                                 rule_ids=[x for x in cited if x not in bad], next_steps=steps, warnings=warns)
    missing = [s for s in ids if s not in out]
    if missing:
        run_warn.append(f"no AI summary returned for: {', '.join(ids[s] for s in missing)}")
    rs = data["run_summary"]
    unverified = sorted({n for n in _quoted_numbers(rs) if not _known(n, allowed)})
    if unverified:
        run_warn.append("run summary quotes numbers not found in the input (verify before use): "
                        + ", ".join(unverified))
    return rs, out, run_warn


def summarise(results: list[SampleResult], profile: Profile, provider: Provider, *,
              genome_size: float | None = None, expected_gc: float | None = None,
              anonymise: bool = False, prompt_path: str | None = None) -> tuple[AISummary, dict]:
    """Build the payload, call the provider once, validate the reply. Returns (summary, payload sent)."""
    prompt, prompt_name, prompt_hash = load_prompt(prompt_path)
    payload, ids = build_payload(results, profile, genome_size, expected_gc, anonymise)
    user = ("Here is the fqci output for one sequencing run. Follow the instructions and return "
            "JSON only.\n\n" + json.dumps(payload, indent=1))
    max_tokens = min(16000, 1500 + 400 * len(payload["samples"]))
    text, meta = provider.generate(prompt, user, OUTPUT_SCHEMA, max_tokens)
    rs, samples, warns = validate(text, payload, ids)
    by_name = {s.sample: s for s in samples.values()}
    if anonymise:  # put real names back into the text shown locally
        for s in by_name.values():
            s.summary = _deanon(s.summary, ids)
            s.next_steps = [_deanon(x, ids) for x in s.next_steps]
        rs = _deanon(rs, ids)
    record = {"provider": provider.name, "model": meta.get("model", provider.model),
              "prompt": prompt_name, "prompt_sha256": prompt_hash, "fqci_version": __version__,
              "anonymised": anonymise, "generated": _dt.datetime.now().isoformat(timespec="seconds"),
              "input_tokens": meta.get("input_tokens"), "output_tokens": meta.get("output_tokens")}
    return AISummary(run_summary=rs, samples=by_name, record=record, warnings=warns), payload


def _deanon(text: str, ids: dict[str, str]) -> str:
    for sid in sorted(ids, key=len, reverse=True):
        text = re.sub(rf"\b{re.escape(sid)}\b", ids[sid], text)
    return text


def to_markdown(ai: AISummary) -> str:
    L = ["## AI summary (generated; verify against the findings above)", "",
         f"Generated by {ai.record['provider']} / {ai.record['model']} with prompt "
         f"{ai.record['prompt']} ({ai.record['prompt_sha256']}) on {ai.record['generated']}. "
         "Advisory only: the rule-based findings and verdicts are authoritative.", ""]
    if ai.warnings:
        L += ["Checks on the reply:", ""] + [f"- {w}" for w in ai.warnings] + [""]
    L += ["### Run", "", ai.run_summary, ""]
    for s in ai.samples.values():
        L += [f"### {s.sample}", "", s.summary, ""]
        if s.next_steps:
            L += ["Next steps:", ""] + [f"- {x}" for x in s.next_steps] + [""]
        if s.rule_ids:
            L += [f"Based on: {', '.join(s.rule_ids)}", ""]
        L += [f"*Check: {w}*" for w in s.warnings] + ([""] if s.warnings else [])
    return "\n".join(L)

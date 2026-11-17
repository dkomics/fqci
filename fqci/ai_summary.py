"""Optional AI summary of fqci findings (off by default).

The rule-based findings stay authoritative. This module sends a compact JSON description of
them (metrics, findings, verdicts, fastp options; never reads or file contents) together with a
fixed, versioned prompt to a language model, and checks the reply before it is shown:

- the reply must match a JSON schema (one entry per sample id);
- rule IDs it cites must exist in that sample's findings (others are removed);
- numbers it quotes must appear in the input (others are listed as unverified).

Sample names can be replaced by codes (S01, S02, ...) before anything leaves the machine.
Two routes: an API provider (Anthropic first; pluggable), or a manual route with no API key:
--ai-export writes a prompt to paste into a chat, --ai-import reads the saved reply back.
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

# List prices, USD per million tokens (input, output), from
# https://platform.claude.com/docs/en/about-claude/pricing. Used only to estimate the cost of a run;
# the Claude Console usage page is authoritative. Haiku 5.5 costs more above 100k prompt tokens.
PRICES_AS_OF = "2026-11-17"
PRICES_USD_PER_MTOK: dict[str, tuple[float, float]] = {
    "claude-sonnet-5-5": (2.0, 10.0),
    "claude-opus-5-5": (4.0, 20.0),
    "claude-fable-5-2": (8.0, 40.0),
    "claude-haiku-5-5": (0.10, 0.50),
}
_HAIKU_LONG = (0.50, 2.50)  # Haiku 5.5, prompts over 100k tokens

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
                "output_tokens": getattr(usage, "output_tokens", None),
                "cache_read_input_tokens": getattr(usage, "cache_read_input_tokens", None),
                "cache_creation_input_tokens": getattr(usage, "cache_creation_input_tokens", None)}
        return text, meta


PROVIDERS = {"anthropic": AnthropicProvider}


def get_provider(name: str, model: str | None = None) -> Provider:
    if name not in PROVIDERS:
        raise AIError(f"unknown AI provider {name!r} (available: {', '.join(PROVIDERS)})")
    return PROVIDERS[name](model=model or DEFAULT_MODEL)  # type: ignore[return-value]


def _price(model: str, input_tokens: int) -> tuple[float, float] | None:
    """List price for a model ID; dated IDs (claude-sonnet-5-5-20261001) match their base name."""
    for base in sorted(PRICES_USD_PER_MTOK, key=len, reverse=True):
        if model == base or model.startswith(base + "-"):
            if base == "claude-haiku-5-5" and input_tokens > 100_000:
                return _HAIKU_LONG
            return PRICES_USD_PER_MTOK[base]
    return None


def estimate_cost(model: str, input_tokens: int | None, output_tokens: int | None) -> float | None:
    """Estimated USD cost of one call at list price, or None if the model or token counts are unknown."""
    if input_tokens is None or output_tokens is None:
        return None
    price = _price(model, input_tokens)
    if price is None:
        return None
    return (input_tokens * price[0] + output_tokens * price[1]) / 1e6


def _usd(v: float) -> str:
    return f"${v:.2f}" if v >= 1 else f"${v:.3f}"


def cost_line(record: dict) -> str | None:
    """One line on tokens and estimated cost for API runs; None for manual (chat) replies."""
    tin, tout = record.get("input_tokens"), record.get("output_tokens")
    if tin is None or tout is None:
        return None
    line = f"Tokens: {tin:,} input + {tout:,} output"
    cache = (record.get("cache_read_input_tokens") or 0) + (record.get("cache_creation_input_tokens") or 0)
    if cache:
        line += f" (+ {cache:,} cache tokens, not priced)"
    cost = record.get("cost_usd_estimate")
    if cost is None:
        return line + f". No price on file for {record.get('model')}; see the Claude Console for the cost."
    pin, pout = _price(record["model"], tin) or (0.0, 0.0)
    return (line + f"; estimated cost {_usd(cost)} at list price (${pin:g} / ${pout:g} per million "
            f"tokens, prices as of {PRICES_AS_OF}). The Claude Console usage page is authoritative.")


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


_SCALED_RE = re.compile(r"(?<![A-Za-z_\d.])(-?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?)\s*(million|thousand|billion)\b",
                        re.I)
_SCALE = {"thousand": 1e3, "million": 1e6, "billion": 1e9}


def _quoted_numbers(text: str) -> list[str]:
    """Numbers as written; '1.9 million' becomes '1,900,000' so it is checked at full value."""
    scaled = []
    for m in _SCALED_RE.finditer(text):
        v = _num(m.group(1))
        if v is not None:
            scaled.append(f"{v * _SCALE[m.group(2).lower()]:,.0f}")
    rest = _SCALED_RE.sub(" ", text)
    return scaled + _NUM_RE.findall(rest)


def _known(n: str, allowed: set[str]) -> bool:
    v = _num(n)
    if v is None:
        return True
    if v >= 1e5:  # large counts are quoted rounded ("about 2.3 million"): accept within 2%
        return any(abs(v - a) <= 0.02 * a for a in _floats(allowed) if a >= 1e5)
    # the quoted value must equal an input figure or one of its rounded forms ("93" for 93.4);
    # rounding the quoted value as well would let 78.8 pass because some other figure is 79
    return f"{v:g}" in allowed or f"{v:.1f}" in allowed


def _floats(allowed: set[str]) -> list[float]:
    return [v for v in (_num(a) for a in allowed) if v is not None]


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
    # numbers a sample's text may quote: that sample's own data plus run-level figures
    run_level = {k: v for k, v in payload.items() if k != "samples"}
    small = {str(i) for i in range(0, 11)} | {f"{float(i):.1f}" for i in range(0, 11)}
    base = _numbers_in(run_level) | small
    per_sample = {s["id"]: base | _numbers_in(s) for s in payload["samples"]}
    allowed = base | _numbers_in(payload)  # the run summary may quote any sample's figures
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
        unverified = sorted({n for n in _quoted_numbers(" ".join([text, *steps])) if not _known(n, per_sample[sid])})
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


@dataclass
class Prepared:
    """Everything needed to ask a model, and to check and restore its answer later."""
    system: str
    user: str
    payload: dict
    ids: dict[str, str]
    prompt_name: str
    prompt_sha256: str
    anonymise: bool


def prepare(results: list[SampleResult], profile: Profile, *, genome_size: float | None = None,
            expected_gc: float | None = None, anonymise: bool = False,
            prompt_path: str | None = None) -> Prepared:
    prompt, prompt_name, prompt_hash = load_prompt(prompt_path)
    payload, ids = build_payload(results, profile, genome_size, expected_gc, anonymise)
    user = ("Here is the fqci output for one sequencing run. Follow the instructions and return "
            "JSON only.\n\n" + json.dumps(payload, indent=1))
    return Prepared(prompt, user, payload, ids, prompt_name, prompt_hash, anonymise)


def finish(prep: Prepared, text: str, provider: str, model: str, meta: dict | None = None) -> AISummary:
    """Validate a reply against the prepared request and restore sample names."""
    meta = meta or {}
    rs, samples, warns = validate(text, prep.payload, prep.ids)
    by_name = {s.sample: s for s in samples.values()}
    if prep.anonymise:  # put real names back into the text shown locally
        for s in by_name.values():
            s.summary = _deanon(s.summary, prep.ids)
            s.next_steps = [_deanon(x, prep.ids) for x in s.next_steps]
        rs = _deanon(rs, prep.ids)
    record = {"provider": provider, "model": meta.get("model", model),
              "prompt": prep.prompt_name, "prompt_sha256": prep.prompt_sha256, "fqci_version": __version__,
              "anonymised": prep.anonymise, "generated": _dt.datetime.now().isoformat(timespec="seconds"),
              "input_tokens": meta.get("input_tokens"), "output_tokens": meta.get("output_tokens")}
    for k in ("cache_read_input_tokens", "cache_creation_input_tokens"):
        if meta.get(k):
            record[k] = meta[k]
    if record["input_tokens"] is not None:
        record["cost_usd_estimate"] = estimate_cost(record["model"], record["input_tokens"],
                                                    record["output_tokens"])
        record["prices_as_of"] = PRICES_AS_OF
    return AISummary(run_summary=rs, samples=by_name, record=record, warnings=warns)


def summarise(results: list[SampleResult], profile: Profile, provider: Provider, *,
              genome_size: float | None = None, expected_gc: float | None = None,
              anonymise: bool = False, prompt_path: str | None = None) -> tuple[AISummary, dict]:
    """Build the request, call the provider once, validate the reply. Returns (summary, payload sent)."""
    prep = prepare(results, profile, genome_size=genome_size, expected_gc=expected_gc,
                   anonymise=anonymise, prompt_path=prompt_path)
    max_tokens = min(16000, 1500 + 400 * len(prep.payload["samples"]))
    text, meta = provider.generate(prep.system, prep.user, OUTPUT_SCHEMA, max_tokens)
    return finish(prep, text, provider.name, provider.model, meta), prep.payload


# ---------------------------------------------------------------- manual route (no API key)
EXPORT_FILE = "fqci_ai_export.json"     # kept locally: payload, id mapping, prompt identity
PROMPT_FILE = "fqci_ai_prompt.txt"      # pasted into a chat


def chat_prompt(prep: Prepared) -> str:
    """One self-contained text for a chat window: instructions, reply format, data."""
    return (prep.system.rstrip() + "\n\n"
            "Reply format: return ONLY a JSON object (no code fences, no text before or after) that "
            "matches this JSON Schema:\n" + json.dumps(OUTPUT_SCHEMA, indent=1) + "\n\n" + prep.user + "\n")


def export(prep: Prepared, outdir) -> tuple[str, str]:
    """Write the chat prompt and the local record needed to import the reply. Returns their paths."""
    from pathlib import Path
    out = Path(outdir)
    (out / PROMPT_FILE).write_text(chat_prompt(prep), encoding="utf-8")
    (out / EXPORT_FILE).write_text(json.dumps({
        "fqci_version": __version__, "prompt": prep.prompt_name, "prompt_sha256": prep.prompt_sha256,
        "anonymise": prep.anonymise, "ids": prep.ids, "payload": prep.payload}, indent=1), encoding="utf-8")
    return str(out / PROMPT_FILE), str(out / EXPORT_FILE)


def extract_json(text: str) -> str:
    """Accept a reply pasted from a chat: strip code fences and any text around the JSON object."""
    t = text.strip()
    m = re.search(r"```(?:json)?\s*(.*?)```", t, re.S)
    if m:
        t = m.group(1).strip()
    if not t.startswith("{"):
        i, j = t.find("{"), t.rfind("}")
        if i == -1 or j <= i:
            raise AIError("no JSON object found in the reply file")
        t = t[i:j + 1]
    return t


def import_reply(reply_path: str, export_path: str, current: Prepared, model: str | None = None) -> AISummary:
    """Check a chat reply against the export it answers and the current run, then build the summary."""
    try:
        with open(export_path, encoding="utf-8") as fh:
            exp = json.load(fh)
    except FileNotFoundError:
        raise AIError(f"{export_path} not found: run with --ai-export first, using the same -o folder") from None
    if exp.get("payload") != current.payload or exp.get("prompt_sha256") != current.prompt_sha256:
        raise AIError("the reply was exported from different data, options or prompt than this run; "
                      "re-run --ai-export with the same inputs and options, then ask again")
    with open(reply_path, encoding="utf-8") as fh:
        text = extract_json(fh.read())
    prep = Prepared(current.system, current.user, exp["payload"], exp["ids"], exp["prompt"],
                    exp["prompt_sha256"], exp["anonymise"])
    return finish(prep, text, "manual (chat)", model or "not recorded")


def _deanon(text: str, ids: dict[str, str]) -> str:
    for sid in sorted(ids, key=len, reverse=True):
        text = re.sub(rf"\b{re.escape(sid)}\b", ids[sid], text)
    return text


def to_markdown(ai: AISummary) -> str:
    L = ["## AI summary (generated; verify against the findings above)", "",
         f"Generated by {ai.record['provider']} / {ai.record['model']} with prompt "
         f"{ai.record['prompt']} ({ai.record['prompt_sha256']}) on {ai.record['generated']}. "
         "Advisory only: the rule-based findings and verdicts are authoritative.", ""]
    cl = cost_line(ai.record)
    if cl:
        L += [cl, ""]
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

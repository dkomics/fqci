"""YAML rule engine.

A rule:
  id: ADAPTER_DIMER
  scope: pair | read           # pair rules see r1_*, r2_*, and pair-level metrics
  title: Adapter dimers
  levels:                      # evaluated in order; first true wins
    fail: "dimer_pct >= 10"
    warn: "dimer_pct >= 1"
  explain: "text with {metric:.1f} placeholders"
  cause: "..."
  fixable: true|false
  explains: [OTHER_RULE_IDS]   # when this rule fires at fail, those become 'secondary'
  actions: {fastp_param: value, ...}   # values may be "{metric}" templates
  advice: "non-trimming follow-up"

Profiles may extend another profile, override params, and patch rules
(e.g. mark a finding as expected for that library type).
Expressions are restricted to names, numbers, comparisons, + - * /, and/or/not,
and min/max/abs. They are interpreted by walking the AST (no eval). A metric that
could not be computed is None ('unknown'): comparisons with it are unknown, and
'unknown or True' is True, so one missing metric does not disable a whole rule.
"""
from __future__ import annotations

import ast
import copy
import re
from collections.abc import Callable
from typing import Any
from dataclasses import dataclass, field
from pathlib import Path

import yaml

PROFILE_DIR = Path(__file__).parent / "profiles"
SEVERITY_ORDER = {"fail": 0, "warn": 1, "info": 2}
FASTP_OPTIONS = {
    "phred64", "adapter_sequence", "adapter_sequence_r2", "adapter_fasta", "detect_adapter_for_pe",
    "disable_adapter_trimming", "trim_poly_g", "poly_g_min_len", "trim_poly_x", "poly_x_min_len",
    "cut_front", "cut_tail", "cut_right", "cut_window_size", "cut_mean_quality",
    "cut_front_mean_quality", "cut_tail_mean_quality", "cut_right_mean_quality",
    "qualified_quality_phred", "unqualified_percent_limit", "average_qual", "n_base_limit",
    "low_complexity_filter", "complexity_threshold", "length_required", "length_limit",
    "trim_front1", "trim_front2", "trim_tail1", "trim_tail2", "max_len1", "max_len2",
    "correction", "overlap_len_require", "overlap_diff_limit", "dedup", "dup_calc_accuracy",
}

_FUNCS: dict[str, Callable[..., Any]] = {"min": min, "max": max, "abs": abs}
_CMP = {ast.Gt: lambda a, b: a > b, ast.GtE: lambda a, b: a >= b, ast.Lt: lambda a, b: a < b,
        ast.LtE: lambda a, b: a <= b, ast.Eq: lambda a, b: a == b, ast.NotEq: lambda a, b: a != b}
_BIN = {ast.Add: lambda a, b: a + b, ast.Sub: lambda a, b: a - b,
        ast.Mult: lambda a, b: a * b, ast.Div: lambda a, b: a / b}


class RuleError(ValueError):
    pass


def _check(node, expr: str, names: set) -> None:
    """Reject anything outside the small expression language; collect referenced names."""
    if isinstance(node, ast.Expression):
        return _check(node.body, expr, names)
    if isinstance(node, ast.Constant):
        if not (node.value is None or isinstance(node.value, (bool, int, float))):
            raise RuleError(f"only numbers, True/False/None allowed as constants in: {expr}")
        return None
    if isinstance(node, ast.Name):
        if node.id not in _FUNCS:
            names.add(node.id)
        return None
    if isinstance(node, ast.BoolOp):
        for v in node.values:
            _check(v, expr, names)
        return None
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.Not, ast.USub)):
        return _check(node.operand, expr, names)
    if isinstance(node, ast.BinOp) and type(node.op) in _BIN:
        _check(node.left, expr, names)
        return _check(node.right, expr, names)
    if isinstance(node, ast.Compare):
        if not all(type(o) in _CMP or isinstance(o, (ast.Is, ast.IsNot)) for o in node.ops):
            raise RuleError(f"unsupported comparison in: {expr}")
        _check(node.left, expr, names)
        for c in node.comparators:
            _check(c, expr, names)
        return None
    if isinstance(node, ast.Call):
        if not (isinstance(node.func, ast.Name) and node.func.id in _FUNCS) or node.keywords \
                or not node.args:
            raise RuleError(f"only min/max/abs calls with positional arguments allowed in: {expr}")
        for a in node.args:
            _check(a, expr, names)
        return None
    raise RuleError(f"disallowed syntax {type(node).__name__} in: {expr}")


def compile_expr(expr: str):
    """Parse and validate an expression; returns (tree, referenced names)."""
    try:
        tree = ast.parse(str(expr), mode="eval")
    except SyntaxError as e:
        raise RuleError(f"invalid expression {expr!r}: {e.msg}") from None
    names: set = set()
    _check(tree, expr, names)
    return tree, names


def _ev(node, ns: dict):
    """Three-valued evaluation: None means 'unknown' (a metric could not be computed)."""
    if isinstance(node, ast.Expression):
        return _ev(node.body, ns)
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return ns.get(node.id)
    if isinstance(node, ast.BoolOp):
        vals = [_ev(v, ns) for v in node.values]
        if isinstance(node.op, ast.And):
            if any(v is not None and not v for v in vals):
                return False
            return None if any(v is None for v in vals) else True
        if any(v is not None and v for v in vals):
            return True
        return None if any(v is None for v in vals) else False
    if isinstance(node, ast.UnaryOp):
        v = _ev(node.operand, ns)
        if v is None:
            return None
        return (not v) if isinstance(node.op, ast.Not) else -v
    if isinstance(node, ast.BinOp):
        a, b = _ev(node.left, ns), _ev(node.right, ns)
        if a is None or b is None or isinstance(a, str) or isinstance(b, str):
            return None
        try:
            return _BIN[type(node.op)](a, b)
        except (ZeroDivisionError, TypeError, OverflowError):
            return None
    if isinstance(node, ast.Compare):
        left = _ev(node.left, ns)
        result = True
        for op, comp in zip(node.ops, node.comparators, strict=True):
            right = _ev(comp, ns)
            if isinstance(op, ast.Is):
                ok = left is right
            elif isinstance(op, ast.IsNot):
                ok = left is not right
            elif left is None or right is None:
                return None
            else:
                try:
                    ok = _CMP[type(op)](left, right)
                except TypeError:
                    return None
            if not ok:
                result = False
            left = right
        return result
    if isinstance(node, ast.Call):
        args = [_ev(a, ns) for a in node.args]
        if any(a is None for a in args):
            return None
        fn = _FUNCS[node.func.id] if isinstance(node.func, ast.Name) else None
        if fn is None:
            return None
        try:
            return fn(*args)
        except (TypeError, ValueError):
            return None
    return None


def evaluate(expr: str, ns: dict):
    """Evaluate a rule expression. Returns None when the answer is unknown."""
    tree, _ = compile_expr(expr)
    return _ev(tree, ns)


@dataclass
class Profile:
    name: str
    description: str
    params: dict
    rules: list[dict]
    base_actions: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    purposes: dict = field(default_factory=dict)
    purpose: str | None = None
    own_params: set = field(default_factory=set)  # set by the profile itself, not inherited from base

    def set_purpose(self, name: str | None) -> None:
        """Apply purpose thresholds. An explicit name (from --purpose) overrides the profile;
        the default purpose only fills values the profile did not set itself."""
        explicit = name is not None
        name = name or self.purpose
        if not self.purposes:
            if explicit:
                raise RuleError(f"profile {self.name} has no purposes")
            return
        if name not in self.purposes:
            raise RuleError(f"unknown purpose {name!r} (available: {', '.join(sorted(self.purposes))})")
        for k, v in self.purposes[name].items():
            if explicit or k not in self.own_params:
                self.params[k] = v
        self.purpose = name


def _resolve(name: str, relative_to: Path | None = None) -> Path:
    """Built-in profile names win; otherwise a path (relative to the extending file, if any)."""
    if "/" not in name and not name.endswith((".yaml", ".yml")):
        p = PROFILE_DIR / f"{name}.yaml"
        if p.is_file():
            return p
        raise RuleError(f"profile not found: {name} (available: {', '.join(available_profiles())}, "
                        "or a path to a .yaml file)")
    p = Path(name)
    if not p.is_absolute() and relative_to is not None:
        p = relative_to.parent / p
    if not p.is_file():
        raise RuleError(f"profile file not found: {p}")
    return p


def _read(path: Path) -> dict:
    try:
        d = yaml.safe_load(path.read_text())
    except yaml.YAMLError as e:
        raise RuleError(f"{path}: invalid YAML: {e}") from None
    if not isinstance(d, dict) or "name" not in d:
        raise RuleError(f"{path}: not a profile (needs at least a 'name:' key)")
    return d


def available_profiles() -> list[str]:
    return sorted(p.stem for p in PROFILE_DIR.glob("*.yaml") if p.stem != "base")


def load_profile(name: str) -> Profile:
    path = _resolve(name)
    d = _read(path)
    chain, seen = [d], {path.resolve()}
    while d.get("extends"):
        path = _resolve(str(d["extends"]), relative_to=path)
        if path.resolve() in seen:
            raise RuleError(f"profile inheritance loop at {path}")
        seen.add(path.resolve())
        d = _read(path)
        chain.append(d)
    params: dict[str, Any] = {}
    rules: dict[str, dict] = {}
    actions: dict[str, Any] = {}
    notes: list[str] = []
    purposes: dict[str, dict] = {}
    default_purpose: str | None = None
    own = {k for layer in chain if layer.get("name") != "base" for k in (layer.get("params") or {})}
    order: list[str] = []
    for layer in reversed(chain):
        params.update(layer.get("params", {}))
        actions.update(layer.get("base_actions", {}))
        notes += layer.get("notes", [])
        for k, v in (layer.get("purposes") or {}).items():
            purposes.setdefault(k, {}).update(v)
        default_purpose = layer.get("default_purpose", default_purpose)
        for r in layer.get("rules", []):
            if r["id"] not in rules:
                order.append(r["id"])
            rules[r["id"]] = copy.deepcopy(r)
        for rid, patch in (layer.get("patch") or {}).items():
            if rid not in rules:
                raise RuleError(f"{layer.get('name')}: patch for unknown rule {rid}")
            rules[rid].update(copy.deepcopy(patch))  # keys are replaced wholesale
        for rid in layer.get("disable", []):
            rules.pop(rid, None)
    final = [rules[i] for i in order if i in rules]
    bad = [k for k in actions if k not in FASTP_OPTIONS]
    bad += [f"{r['id']}.{k}" for r in final for k in (r.get("actions") or {}) if k not in FASTP_OPTIONS]
    if bad:
        raise RuleError(f"{chain[0]['name']}: not fastp options: {', '.join(bad)} "
                        "(check indentation - thresholds belong under params:)")
    for r in final:  # validate early
        for e in list(r.get("levels", {}).values()):
            compile_expr(e)
    prof = Profile(name=chain[0]["name"], description=chain[0].get("description", ""),
                   params=params, rules=final, base_actions=actions, notes=notes,
                   purposes=purposes, purpose=default_purpose, own_params=own)
    prof.set_purpose(None)
    return prof


_FIELD_RE = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)(?::([^{}]*))?\}")


def render(template, ns: dict, strict: bool = False):
    """Fill {name} / {name:fmt} from ns.

    For prose, missing values render as 'n/a'. With strict=True (used for fastp options)
    a missing value raises RuleError, so a typo can never reach a generated command."""
    if not isinstance(template, str):
        return template

    def sub(m):
        v = ns.get(m.group(1))
        if v is None:
            if strict:
                raise RuleError(f"no value for {{{m.group(1)}}} in action {template!r}")
            return "n/a"
        try:
            return format(v, m.group(2) or "")
        except (ValueError, TypeError):
            return str(v)
    return _FIELD_RE.sub(sub, template)


@dataclass
class Finding:
    rule_id: str
    title: str
    severity: str
    scope: str
    read: str | None
    explain: str
    cause: str
    fixable: bool
    expected: bool
    expected_note: str
    advice: str
    actions: dict
    secondary_to: list[str] = field(default_factory=list)


def run_rules(profile: Profile, ns: dict, scope: str, read: str | None = None) -> list[Finding]:
    out = []
    for r in profile.rules:
        if r.get("scope", "read") != scope:
            continue
        sev = None
        for level, expr in r.get("levels", {}).items():
            res = evaluate(expr, ns)
            if res:
                sev = level
                break
        if sev is None:
            continue
        acts = {k: render(v, ns, strict=True) for k, v in (r.get("actions") or {}).items()}
        out.append(Finding(
            rule_id=r["id"], title=r.get("title", r["id"]), severity=sev, scope=scope, read=read,
            explain=render(r.get("explain", ""), ns), cause=render(r.get("cause", ""), ns),
            fixable=bool(r.get("fixable", False)), expected=bool(r.get("expected", False)),
            expected_note=render(r.get("expected_note", ""), ns),
            advice=render(r.get("advice", ""), ns), actions=acts))
    return out


def mark_secondary(profile: Profile, findings: list[Finding]) -> None:
    explains = {r["id"]: r.get("explains", []) for r in profile.rules}
    primaries = [f for f in findings if f.severity == "fail" and explains.get(f.rule_id)]
    for f in findings:
        for p in primaries:
            if f.rule_id in explains[p.rule_id] and p is not f and p.rule_id not in f.secondary_to:
                f.secondary_to.append(p.rule_id)

"""fqci - profile-aware interpretation of FastQC reports."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .compare import compare
from .interpret import group, interpret_sample, natural_key
import zipfile

from .parser import discover
from .html_report import to_html
from .report import to_json, to_markdown
from .rules import RuleError, available_profiles, load_profile


def _common(p):
    p.add_argument("-p", "--profile", default="wgs_isolate",
                   help=f"library profile name or YAML path ({', '.join(available_profiles())})")
    p.add_argument("--genome-size", type=float, help="genome size in bp, e.g. 4.03e6 (not Mb)")
    p.add_argument("--expected-gc", type=float, help="expected mean %%GC, e.g. 47.5")
    p.add_argument("--purpose", help="what the reads are for; sets depth thresholds "
                   "(assembly: 50x warn/20x fail; mapping: 30x/15x)")
    p.add_argument("--set", action="append", default=[], metavar="PARAM=VALUE",
                   help="override a profile parameter, e.g. --set min_cov_warn=30 (repeatable)")


def _apply_sets(prof, sets):
    """--set PARAM=VALUE, converted to the type of the profile's existing value."""
    for item in sets:
        k, sep, v = item.partition("=")
        if not sep or k not in prof.params:
            raise SystemExit(f"--set {item}: unknown parameter. Known: {', '.join(sorted(prof.params))}")
        old = prof.params[k]
        try:
            if isinstance(old, bool):
                if v.lower() not in ("true", "false"):
                    raise ValueError
                val = v.lower() == "true"
            elif isinstance(old, (int, float)):
                f = float(v)
                val = int(f) if isinstance(old, int) and f.is_integer() and "." not in v else f
            else:
                val = v
        except ValueError:
            kind = {bool: "true/false", int: "a number", float: "a number"}.get(type(old), "a value")
            raise SystemExit(f"--set {item}: {k} needs {kind} (currently {old!r})") from None
        prof.params[k] = val


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="fqci", description=__doc__)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("interpret", help="interpret FastQC reports and recommend fastp parameters")
    a.add_argument("inputs", nargs="+", help="*_fastqc.zip, *_fastqc/ folders, or directories containing them")
    _common(a)
    a.add_argument("-o", "--outdir", default="fqci_out")
    a.add_argument("--title", default="FastQC interpretation")

    c = sub.add_parser("compare", help="compare raw vs trimmed FastQC reports")
    c.add_argument("--before", nargs="+", required=True)
    c.add_argument("--after", nargs="+", required=True)
    _common(c)
    c.add_argument("-o", "--out", default="fqci_compare.md")
    c.add_argument("--fastp-dir", help="folder with <sample>.fastp.json files, to explain read losses")

    sub.add_parser("profiles", help="list available profiles")
    a.add_argument("--threads", type=int, default=4, help="fastp --thread value in the generated commands")

    args = ap.parse_args(argv)
    try:
        return _run(args)
    except (RuleError, ValueError, OSError, zipfile.BadZipFile) as e:
        print(f"fqci: error: {e}", file=sys.stderr)
        return 2


def _run(args) -> int:
    if args.cmd == "profiles":
        for n in available_profiles():
            print(f"{n}\t{load_profile(n).description.strip()}")
        return 0
    if args.genome_size is not None and args.genome_size < 1e4:
        raise ValueError(f"--genome-size {args.genome_size:g} looks like Mb or kb; give base pairs "
                         f"(e.g. {args.genome_size:g}e6 for {args.genome_size:g} Mb)")
    prof = load_profile(args.profile)
    prof.set_purpose(args.purpose)
    _apply_sets(prof, args.set)
    kw = {"genome_size": args.genome_size, "expected_gc": args.expected_gc}
    if args.cmd == "interpret":
        reps = discover(args.inputs)
        if not reps:
            print("no FastQC reports found", file=sys.stderr)
            return 2
        res = [interpret_sample(s, r, prof, threads=args.threads, **kw)
               for s, r in sorted(group(reps).items(), key=lambda kv: natural_key(kv[0]))]
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "fqci_report.md").write_text(to_markdown(res, prof, args.title))
        (out / "fqci_report.json").write_text(to_json(res, prof))
        (out / "fqci_report.html").write_text(to_html(res, prof, args.title), encoding="utf-8")
        real = [r for r in res if not r.undetermined]
        (out / "fastp_commands.sh").write_text(
            "#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(r.fastp_command for r in real) + "\n")
        (out / "fastp_params.json").write_text(json.dumps({r.sample: r.fastp_params for r in real}, indent=2))
        for r in res:
            print(f"{r.sample}\t{'UNASSIGNED' if r.undetermined else r.worst.upper()}")
        print(f"report: {out / 'fqci_report.md'} and {out / 'fqci_report.html'}", file=sys.stderr)
        return 0
    if args.cmd == "compare":
        Path(args.out).write_text(compare(discover(args.before), discover(args.after), prof,
                                          fastp_dir=args.fastp_dir, **kw))
        print(f"written: {args.out}", file=sys.stderr)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())

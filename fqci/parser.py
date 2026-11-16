"""Parse FastQC `fastqc_data.txt` from a .zip, an extracted folder, or the file itself.

The zip is read in memory; nothing is extracted to disk.
"""
from __future__ import annotations

import io
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Module:
    name: str
    status: str
    header: list[str] = field(default_factory=list)
    rows: list[list[str]] = field(default_factory=list)
    extra: dict[str, str] = field(default_factory=dict)  # '#Key\tvalue' lines


@dataclass
class FastQCReport:
    source: str
    version: str
    modules: dict[str, Module]

    @property
    def basic(self) -> dict[str, str]:
        m = self.modules.get("Basic Statistics")
        return {r[0]: r[1] for r in m.rows} if m else {}

    @property
    def filename(self) -> str:
        return self.basic["Filename"]

    def status(self, module: str) -> str | None:
        m = self.modules.get(module)
        return m.status if m else None


def parse_text(text: str, source: str = "<text>") -> FastQCReport:
    version = ""
    modules: dict[str, Module] = {}
    cur: Module | None = None
    pending_hash: list[str] = []
    for raw in text.splitlines():
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("##FastQC"):
            version = line.split("\t")[1] if "\t" in line else ""
        elif line.startswith(">>END_MODULE"):
            if cur is not None:
                if pending_hash and not cur.header:
                    cur.header = pending_hash[-1].lstrip("#").split("\t")
                modules[cur.name] = cur
            cur, pending_hash = None, []
        elif line.startswith(">>"):
            name, _, status = line[2:].partition("\t")
            cur, pending_hash = Module(name=name, status=status.strip()), []
        elif cur is not None and line.startswith("#"):
            pending_hash.append(line)
        elif cur is not None:
            if pending_hash:
                # last '#' line is the column header; earlier ones are key/value
                cur.header = pending_hash[-1].lstrip("#").split("\t")
                for kv in pending_hash[:-1]:
                    k, _, v = kv.lstrip("#").partition("\t")
                    cur.extra[k] = v
                pending_hash = []
            cur.rows.append(line.split("\t"))
    if not modules:
        raise ValueError(f"{source}: no FastQC modules found - is this fastqc_data.txt?")
    rep = FastQCReport(source=source, version=version, modules=modules)
    if not rep.basic.get("Filename"):
        # without it we cannot know the FASTQ name and would write a wrong fastp command
        raise ValueError(f"{source}: FastQC report has no 'Filename' in Basic Statistics")
    return rep


def load(path: str | Path) -> FastQCReport:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{p}: no such file or directory")
    if p.is_dir():
        f = p / "fastqc_data.txt"
        if not f.exists():
            raise FileNotFoundError(f"{p}: no fastqc_data.txt inside")
        return parse_text(f.read_text(), str(p))
    if p.suffix == ".zip":
        with zipfile.ZipFile(p) as z:
            names = [n for n in z.namelist() if n.endswith("fastqc_data.txt")]
            if not names:
                raise FileNotFoundError(f"{p}: no fastqc_data.txt in zip")
            with z.open(sorted(names, key=len)[0]) as fh:
                return parse_text(io.TextIOWrapper(fh, encoding="utf-8").read(), str(p))
    return parse_text(p.read_text(), str(p))


def discover(paths: list[str | Path]) -> list[FastQCReport]:
    """Accept files and directories; directories are searched for *_fastqc.zip
    and *_fastqc/ folders (a zip and its extracted folder are not double-counted)."""
    found: dict[str, FastQCReport] = {}
    for path in map(Path, paths):
        cands: list[Path]
        if path.is_dir() and not (path / "fastqc_data.txt").exists():
            cands = sorted(path.glob("*_fastqc.zip")) + sorted(
                d for d in path.glob("*_fastqc") if d.is_dir())
        else:
            cands = [path]
        for c in cands:
            rep = load(c)
            prev = found.get(rep.filename)
            if prev is None:
                found[rep.filename] = rep
            elif _report_stem(prev.source) != _report_stem(rep.source):
                # same FASTQ name from two different reports (e.g. a top-up run): never drop silently
                raise ValueError(f"two FastQC reports for the same file name {rep.filename!r}: "
                                 f"{prev.source} and {rep.source}. Analyse them separately or rename.")
    return list(found.values())


def _report_stem(source: str) -> str:
    """X_fastqc.zip and its extracted folder X_fastqc/ count as the same report."""
    p = Path(source)
    return str(p.with_suffix("")) if p.suffix == ".zip" else str(p)


_READ_RE = re.compile(r"^(?P<sample>.+?)[._-](?P<r>[Rr])?(?P<read>[12])(?:_001)?$")
_SUFFIX_RE = re.compile(r"[._-](trim(med)?|clean(ed)?|filt(ered)?|qc|val_[12]|paired)$", re.I)
_EXT_RE = re.compile(r"(\.(fastq|fq|txt|sam|bam))?(\.(gz|bz2|zst))?$", re.I)


def sample_and_read(filename: str) -> tuple[str, int | None]:
    stem = _EXT_RE.sub("", Path(filename).name)
    prev = None
    while prev != stem:  # Vc1_R1.trim / Vc1_R1_trimmed / Vc1_R1.clean -> Vc1_R1
        prev, stem = stem, _SUFFIX_RE.sub("", stem)
    m = _READ_RE.match(stem)
    if m:
        return m.group("sample"), int(m.group("read"))
    return stem, None


def read_label_is_explicit(filename: str) -> bool:
    """True for _R1/_r2 style names; False for bare _1/_2 (which may be two unrelated files)."""
    stem = _EXT_RE.sub("", Path(filename).name)
    prev = None
    while prev != stem:
        prev, stem = stem, _SUFFIX_RE.sub("", stem)
    m = _READ_RE.match(stem)
    return bool(m and m.group("r"))


def pos_range(label: str) -> tuple[int, int]:
    """'10-11' -> (10, 11); '7' -> (7, 7)."""
    a, _, b = label.partition("-")
    return int(a), int(b or a)

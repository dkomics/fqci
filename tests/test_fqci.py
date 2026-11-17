import json as _json
import os
import shutil
import subprocess
import subprocess as _sp
import sys as _sys
from pathlib import Path

import pytest

from fqci import ai_summary as A
from fqci import metrics, parser
from fqci.cli import main
from fqci.compare import compare
from fqci.interpret import depth_dup_expectation, group, interpret_sample, natural_key
from fqci.rules import RuleError, available_profiles, evaluate, load_profile

D = Path(__file__).parent / "data"
REAL = os.environ.get("FQCI_REAL_DATA")  # optional: directory of real FastQC zips


@pytest.fixture(scope="module")
def raw():
    return parser.discover([D / "raw"])


@pytest.fixture(scope="module")
def trimmed():
    return parser.discover([D / "trimmed"])


# ---------------------------------------------------------------- parser
@pytest.mark.parametrize("fn,exp", [
    ("Vc1_Nepal_R1.fastq.gz", ("Vc1_Nepal", 1)),
    ("S1_S1_L001_R2_001.fastq.gz", ("S1_S1_L001", 2)),
    ("ERR123_1.fq.gz", ("ERR123", 1)),
    ("sample.R2.fastq", ("sample", 2)),
    ("Vc1_Nepal_R1.trim.fastq.gz", ("Vc1_Nepal", 1)),
    ("x_R1_val_1.fq.gz", ("x", 1)),
    ("single.fastq.gz", ("single", None)),
])
def test_sample_and_read(fn, exp):
    assert parser.sample_and_read(fn) == exp


def test_parse_zip_and_modules(raw):
    assert len(raw) == 2
    r = raw[0]
    assert r.version.startswith("0.")
    for mod in ("Basic Statistics", "Per base sequence quality", "Sequence Duplication Levels",
                "Overrepresented sequences", "Adapter Content"):
        assert mod in r.modules
    dup = r.modules["Sequence Duplication Levels"]
    assert "Total Deduplicated Percentage" in dup.extra
    assert dup.header[0] == "Duplication Level"


def test_parse_extracted_folder(tmp_path, raw):
    import zipfile
    with zipfile.ZipFile(D / "raw" / "sim_R1_fastqc.zip") as z:
        z.extractall(tmp_path)
    folder = next(tmp_path.glob("*_fastqc"))
    rep = parser.load(folder)
    assert rep.filename == "sim_R1.fastq.gz"
    # zip + extracted folder in the same directory are not double counted
    shutil.copy(D / "raw" / "sim_R1_fastqc.zip", tmp_path)
    assert len(parser.discover([tmp_path])) == 1


def test_parse_rejects_non_fastqc():
    with pytest.raises(ValueError):
        parser.parse_text("hello\nworld\n")


# ---------------------------------------------------------------- metrics
def test_metrics_detect_dimers(raw):
    ms = {parser.sample_and_read(r.filename)[1]: metrics.compute(r) for r in raw}
    # simulation truth: 60,000 / 180,000 pairs are dimers (33.3%)
    assert 30 <= ms[1]["dimer_pct"] <= 34
    assert ms[1]["observed_i7"] == "ATCACG"
    assert 30 <= ms[2]["pct_reads_meanq_lt20"] <= 36  # dimer mates in R2
    assert ms[1]["n_reads"] == 180000 and ms[1]["read_len"] == 101


def test_metrics_clean_after_trim(trimmed):
    for r in trimmed:
        m = metrics.compute(r)
        assert m["dimer_pct"] == 0
        assert m["q_drop_lq20"] is None


# ---------------------------------------------------------------- rules
def test_safe_eval():
    assert evaluate("a > 1 and b < 2", {"a": 2, "b": 1}) is True
    assert evaluate("a > 1", {"a": None}) is None          # missing metric -> not applicable
    assert evaluate("abs(a - b) > 3", {"a": 1, "b": 5}) is True
    for bad in ("__import__('os').system('true')", "a.__class__", "open('x')", "[x for x in a]"):
        with pytest.raises(RuleError):
            evaluate(bad, {"a": 1})


@pytest.mark.parametrize("name", available_profiles())
def test_profiles_load_and_run(name, raw):
    prof = load_profile(name)
    for s, reads in group(raw).items():
        res = interpret_sample(s, reads, prof, genome_size=4e5)
        assert res.fastp_command.startswith("fastp ")


def test_unknown_profile():
    with pytest.raises(RuleError):
        load_profile("nonexistent_profile")


# ---------------------------------------------------------------- interpretation
def test_wgs_interpretation(raw):
    res = interpret_sample("sim", group(raw)["sim"], load_profile("wgs_isolate"), genome_size=4e5)
    ids = {f.rule_id: f for f in res.findings}
    assert ids["ADAPTER_DIMER"].severity == "fail" and not ids["ADAPTER_DIMER"].secondary_to
    for sec in ("QUAL_LOW_READS", "COMP_BODY", "DUPLICATION_DEPTH"):
        assert any(f.rule_id == sec and f.secondary_to == ["ADAPTER_DIMER"] for f in res.findings)
    assert not any(f.rule_id == "DUPLICATION" for f in res.findings)  # fixed rule off with genome size
    p = res.fastp_params
    assert p["adapter_sequence"].startswith("AGATCGGAAGAGCACACG")
    assert p["adapter_sequence_r2"].startswith("AGATCGGAAGAGCGTCG")
    assert p["detect_adapter_for_pe"] is True and p["length_required"] == 50
    assert "dedup" not in p
    # 60x of sequence -> 48x expected depth after calibration: below the 50x assembly threshold
    assert res.verdict == "FAIL; not fixable by trimming: COVERAGE (warn)"
    pm = res.metrics["pair"]
    assert 55 <= pm["est_coverage_raw"] <= 65
    assert abs(pm["est_coverage"] - pm["est_coverage_raw"] * 0.8) < 1e-9  # >=30x uses depth_factor


def test_amplicon_marks_expected(raw):
    res = interpret_sample("sim", group(raw)["sim"], load_profile("amplicon"))
    exp = {f.rule_id for f in res.findings if f.expected}
    assert "DUPLICATION" in exp
    assert "cut_tail" not in res.fastp_params and "cut_right" not in res.fastp_params
    assert not any(f.rule_id == "COVERAGE" for f in res.findings)


def test_single_end(raw):
    r1 = [r for r in raw if "_R1" in r.filename]
    res = interpret_sample("sim", {1: r1[0]}, load_profile("wgs_isolate"))
    assert " -I " not in res.fastp_command
    assert "adapter_sequence_r2" not in res.fastp_params
    assert "detect_adapter_for_pe" not in res.fastp_params


def test_compare_resolves(raw, trimmed):
    md = compare(raw, trimmed, load_profile("wgs_isolate"), genome_size=4e5)
    assert "Nothing left to fix after trimming." in md
    assert "| R1 | Reads kept % | | 66.7 |" in md


def test_cli(tmp_path):
    out = tmp_path / "o"
    assert main(["interpret", str(D / "raw"), "-o", str(out), "--genome-size", "4e5"]) == 0
    for f in ("fqci_report.md", "fqci_report.json", "fastp_commands.sh", "fastp_params.json"):
        assert (out / f).stat().st_size > 0
    assert main(["compare", "--before", str(D / "raw"), "--after", str(D / "trimmed"),
                 "-o", str(tmp_path / "c.md")]) == 0


@pytest.mark.skipif(shutil.which("fastp") is None, reason="fastp not on PATH")
def test_generated_fastp_command_runs(tmp_path, raw):
    """Regenerate the simulated FASTQ, run the recommended command verbatim, check dimers are gone."""
    subprocess.run([_sys.executable, str(D / "simulate.py"), "3000", "1500"], cwd=tmp_path, check=True)
    res = interpret_sample("sim", group(raw)["sim"], load_profile("wgs_isolate"))
    subprocess.run(res.fastp_command, shell=True, cwd=tmp_path, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    import gzip
    with gzip.open(tmp_path / "sim_R1.trim.fastq.gz", "rt") as fh:
        n = sum(1 for i, _ in enumerate(fh) if i % 4 == 0)
    assert n == 3000  # every genomic pair kept, every dimer pair removed


@pytest.mark.skipif(not REAL, reason="set FQCI_REAL_DATA to a directory of real FastQC zips")
def test_real_data_smoke():
    reps = parser.discover([REAL])
    for s, reads in group(reps).items():
        interpret_sample(s, reads, load_profile("wgs_isolate"), genome_size=4.03e6)


# ---------------------------------------------------------------- duplication vs depth, sorting


def _haiti_like(n_reads=26839007, dedup=32.4):
    """Clean single-end report (trimmed fixture) relabelled as a deep SE run."""
    import re as _re
    import zipfile
    with zipfile.ZipFile(D / "trimmed" / "sim_R1.trim_fastqc.zip") as z:
        name = next(n for n in z.namelist() if n.endswith("fastqc_data.txt"))
        txt = z.read(name).decode()
    txt = _re.sub(r"Total Sequences\t\d+", f"Total Sequences\t{n_reads}", txt)
    txt = _re.sub(r"#Total Deduplicated Percentage\t[\d.]+", f"#Total Deduplicated Percentage\t{dedup}", txt)
    txt = txt.replace("sim_R1.trim.fastq.gz", "Vc_1786_Haiti_SE.fastq.gz")
    return parser.parse_text(txt, "haiti_like")


def test_depth_expectation_formula():
    lam, exp = depth_dup_expectation(26839007, 4.03e6)
    assert round(lam, 2) == 3.33 and round(exp, 1) == 29.0
    assert depth_dup_expectation(0, 4e6)[1] == 100.0


def test_deep_se_duplication_is_explained_by_depth():
    rep = _haiti_like()
    res = interpret_sample("Vc_1786_Haiti_SE", {1: rep}, load_profile("wgs_isolate"), genome_size=4.03e6)
    dup = [f for f in res.findings if f.rule_id.startswith("DUPLICATION")]
    assert [(f.rule_id, f.severity) for f in dup] == [("DUPLICATION_DEPTH", "info")]
    cov = next(f for f in res.findings if f.rule_id == "COVERAGE")
    assert " reads remain" in cov.explain and "pairs" not in cov.explain


def test_fixed_duplication_used_without_genome_size():
    res = interpret_sample("x", {1: _haiti_like()}, load_profile("wgs_isolate"))
    assert any(f.rule_id == "DUPLICATION" and f.severity == "fail" for f in res.findings)
    assert not any(f.rule_id == "DUPLICATION_DEPTH" for f in res.findings)


def test_excess_duplication_still_fails():
    res = interpret_sample("x", {1: _haiti_like(dedup=10.0)}, load_profile("wgs_isolate"), genome_size=4.03e6)
    assert any(f.rule_id == "DUPLICATION_DEPTH" and f.severity == "fail" for f in res.findings)


def test_natural_sort():
    names = ["Vc10_Nepal", "Vc1_Nepal", "Vc_1786_Haiti_SE", "Vc2_Nepal", "Vc24_Nepal"]
    assert sorted(names, key=natural_key) == ["Vc1_Nepal", "Vc2_Nepal", "Vc10_Nepal", "Vc24_Nepal",
                                             "Vc_1786_Haiti_SE"]


def test_profile_rejects_non_fastp_action(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: bad\nextends: wgs_isolate\nbase_actions:\n  dup_excess_fail: 30\n")
    with pytest.raises(RuleError, match="not fastp options"):
        load_profile(str(bad))


# ---------------------------------------------------------------- MiSeq 2x151: end bias, undetermined, --set, tiles
def _bins(read_len):
    """FastQC 0.11/0.12 default grouping: 1-9 single, then 5-wide (10-14 ...), last bin shorter."""
    b = [(i, i) for i in range(1, 10)]
    i = 10
    while i <= read_len:
        b.append((i, min(i + 4, read_len)))
        i += 5
    return b


def _miseq151(fname, n=1007380, end_bias=True, start_bias=True, lens=((151, 0.8), (120, 0.2)),
              tile_dev=0.0, gc_modes=(65,), read_len=151, adapter_cols=True):
    """FastQC 0.12-style report for a 2x151 MiSeq library (H37Rv-like): 2-column duplication
    module, 1-9 then 5-wide position bins ending '150-151' for 151 bp, as in real output."""
    L = ["##FastQC\t0.12.1", ">>Basic Statistics\tpass", "#Measure\tValue", f"Filename\t{fname}",
         "Encoding\tSanger / Illumina 1.9", f"Total Sequences\t{n}", f"Sequence length\t35-{read_len}",
         "%GC\t65", ">>END_MODULE", ">>Per base sequence quality\tpass",
         "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile"]
    bins = _bins(read_len)
    def lab(a, b):
        return str(a) if a == b else f"{a}-{b}"
    L += [f"{lab(a, b)}\t36\t37\t34\t38\t30\t39" for a, b in bins] + [">>END_MODULE"]
    L += [">>Per tile sequence quality\tpass", "#Tile\tBase\tMean"]
    L += [f"{t}\t1\t{tile_dev if t < 1105 else 0.0}" for t in range(1101, 1120)] + [">>END_MODULE"]
    L += [">>Per sequence quality scores\tpass", "#Quality\tCount", "15\t10", "36\t1000", ">>END_MODULE"]
    L += [">>Per base sequence content\tfail", "#Base\tG\tA\tT\tC"]
    for a, b in bins:
        if start_bias and b <= 9:
            row = (49.8, 15.2, 7.9, 27.1)
        elif end_bias and a >= read_len - 1:
            row = (41.5, 22.2, 17.8, 18.5)
        else:
            row = (32.6, 17.4, 17.4, 32.6)
        L.append(f"{lab(a, b)}\t" + "\t".join(map(str, row)))
    L += [">>END_MODULE", ">>Per sequence GC content\tpass", "#GC Content\tCount"]
    import math as _m
    L += [f"{g}\t{sum(1000 * _m.exp(-((g - mu) / 4) ** 2) for mu in gc_modes):.1f}" for g in range(0, 101)]
    L += [">>END_MODULE", ">>Per base N content\tpass", "#Base\tN-Count", "1\t0.0", ">>END_MODULE",
          ">>Sequence Length Distribution\twarn", "#Length\tCount"]
    L += [f"{ln}\t{n * fr}" for ln, fr in lens] + [">>END_MODULE"]
    L += [">>Sequence Duplication Levels\tpass", "#Total Deduplicated Percentage\t81.3",
          "#Duplication Level\tPercentage of total", "1\t81.0", ">100\t0.5", ">1k\t0.2", ">>END_MODULE",
          ">>Overrepresented sequences\tpass", ">>END_MODULE"]
    if adapter_cols:
        L += [">>Adapter Content\tpass",
              "#Position\tIllumina Universal Adapter\tNextera Transposase Sequence\tPolyA\tPolyG",
              "1\t0.0\t0.0\t0.0\t0.0", f"{read_len}\t0.0\t0.1\t0.0\t0.0", ">>END_MODULE"]
    else:  # FastQC 0.11.x: no PolyA/PolyG columns
        L += [">>Adapter Content\tpass", "#Position\tIllumina Universal Adapter\tNextera Transposase Sequence",
              "1\t0.0\t0.0", f"{read_len}\t0.0\t0.1", ">>END_MODULE"]
    return parser.parse_text("\n".join(L) + "\n", fname)


def _pair(sample, **kw):
    return {1: _miseq151(f"{sample}_R1_001.fastq.gz", **kw), 2: _miseq151(f"{sample}_R2_001.fastq.gz", **kw)}


def test_last_cycle_bias_is_not_body_bias():
    res = interpret_sample("H37Rv_S18_L001", _pair("H37Rv_S18_L001"), load_profile("wgs_isolate"),
                           genome_size=4.41e6, expected_gc=65.6)
    ids = {(f.rule_id, f.read): f.severity for f in res.findings}
    assert not any(r == "COMP_BODY" for r, _ in ids)
    assert ids[("COMP_END", "R1")] == "warn" and ids[("COMP_START", "R1")] == "info"
    assert res.fastp_params["trim_tail1"] == 1 and res.fastp_params["trim_tail2"] == 1
    assert not any(f.rule_id.startswith("GC_") for f in res.findings)
    assert res.verdict == "WARN: fixable by trimming"


def test_real_body_bias_still_fails():
    rep = _miseq151("x_R1_001.fastq.gz")
    rows = rep.modules["Per base sequence content"].rows
    rows[20][1:] = ["50", "10", "10", "30"]  # a mid-read bin
    m = metrics.compute(rep)
    assert m["comp_bias_body"] >= 20


def test_coverage_uses_mean_length():
    res = interpret_sample("s", _pair("s"), load_profile("wgs_isolate"), genome_size=4.41e6)
    pm = res.metrics["pair"]
    assert abs(pm["mean_len_used"] - (0.8 * 151 + 0.2 * 120)) < 0.01
    assert abs(pm["est_coverage_raw"] - 1007380 * 2 * pm["mean_len_used"] / 4.41e6) < 0.01
    assert pm["depth_factor_used"] == 0.8 and abs(pm["est_coverage"] - 0.8 * pm["est_coverage_raw"]) < 1e-9


def test_single_end_drops_trim_tail2():
    res = interpret_sample("s", {1: _miseq151("s_R1_001.fastq.gz")}, load_profile("wgs_isolate"))
    assert "trim_tail1" in res.fastp_params and "trim_tail2" not in res.fastp_params


def test_undetermined_excluded(tmp_path):
    import zipfile
    d = tmp_path / "qc"
    d.mkdir()
    for s, kw in (("H37Rv_S18_L001", {}), ("Undetermined_S0_L001", {"gc_modes": (44, 65)})):
        for r in (1, 2):
            txt = _to_text(_miseq151(f"{s}_R{r}_001.fastq.gz", **kw))
            with zipfile.ZipFile(d / f"{s}_R{r}_001_fastqc.zip", "w") as z:
                z.writestr(f"{s}_R{r}_001_fastqc/fastqc_data.txt", txt)
    out = tmp_path / "o"
    assert main(["interpret", str(d), "-o", str(out), "--genome-size", "4.41e6"]) == 0
    sh = (out / "fastp_commands.sh").read_text()
    assert "H37Rv_S18_L001" in sh and "Undetermined" not in sh
    md = (out / "fqci_report.md").read_text()
    assert "n/a (unassigned reads)" in md and "could not assign" in md


def _to_text(rep):
    L = [f"##FastQC\t{rep.version}"]
    for m in rep.modules.values():
        L.append(f">>{m.name}\t{m.status}")
        L += [f"#{k}\t{v}" for k, v in m.extra.items()]
        if m.header:
            L.append("#" + "\t".join(m.header))
        L += ["\t".join(r) for r in m.rows]
        L.append(">>END_MODULE")
    return "\n".join(L) + "\n"


def test_set_overrides_threshold(tmp_path, capsys):
    out = tmp_path / "o"
    assert main(["interpret", str(D / "raw"), "-o", str(out), "--genome-size", "4e5",
                 "--set", "min_cov_warn=100", "--set", "min_cov_fail=80"]) == 0
    import json
    d = json.loads((out / "fqci_report.json").read_text())
    cov = next(f for f in d["samples"][0]["findings"] if f["rule_id"] == "COVERAGE")
    assert cov["severity"] == "fail"  # ~61x < 80
    with pytest.raises(SystemExit):
        main(["interpret", str(D / "raw"), "-o", str(out), "--set", "no_such_param=1"])


def test_tile_explained_by_dimers_in_profile():
    dimer = next(r for r in load_profile("wgs_isolate").rules if r["id"] == "ADAPTER_DIMER")
    assert "TILE" in dimer["explains"]


def test_tile_primary_without_dimers():
    res = interpret_sample("s", _pair("s", tile_dev=-12.0), load_profile("wgs_isolate"))
    tile = [f for f in res.findings if f.rule_id == "TILE"]
    assert tile and tile[0].severity == "fail" and not tile[0].secondary_to


def test_compare_shows_tile(raw, trimmed):
    md = compare(raw, trimmed, load_profile("wgs_isolate"))
    assert "Per tile sequence quality" in md


# ---------------------------------------------------------------- minimum length and cut_tail
@pytest.mark.parametrize("read_len,expected", [(151, 50), (101, 50), (76, 38), (50, 36)])
def test_min_length_capped(read_len, expected):
    rep = _miseq151("s_R1_001.fastq.gz")
    rep.modules["Basic Statistics"].rows = [r if r[0] != "Sequence length" else ["Sequence length", f"35-{read_len}"]
                                            for r in rep.modules["Basic Statistics"].rows]
    res = interpret_sample("s", {1: rep}, load_profile("wgs_isolate"))
    assert res.fastp_params["length_required"] == expected


def test_quality_trimming_uses_cut_tail():
    for prof in ("wgs_isolate", "shotgun_metagenome"):
        res = interpret_sample("s", _pair("s"), load_profile(prof))
        assert res.fastp_params.get("cut_tail") is True
        assert "cut_right" not in res.fastp_params
        assert "--cut_tail" in res.fastp_command


# ---------------------------------------------------------------- depth calibration, purpose, GC bias, compare loss check, run summary


def test_depth_factor_low_below_cutoff():
    res = interpret_sample("s", _pair("s", n=150000), load_profile("wgs_isolate"), genome_size=4.41e6)
    pm = res.metrics["pair"]
    assert pm["est_coverage_raw"] < 30 and pm["depth_factor_used"] == 0.7
    assert abs(pm["est_coverage"] - 0.7 * pm["est_coverage_raw"]) < 1e-9


def test_purpose_mapping_relaxes_thresholds():
    prof_a, prof_m = load_profile("wgs_isolate"), load_profile("wgs_isolate")
    prof_m.set_purpose("mapping")
    assert prof_a.purpose == "assembly" and prof_m.params["min_cov_warn"] == 30
    pair = _pair("s", n=600000)  # ~39x raw -> ~31x expected
    def cov(prof):
        return next(f for f in interpret_sample("s", pair, prof, genome_size=4.41e6).findings
                    if f.rule_id == "COVERAGE")
    assert cov(prof_a).severity == "warn" and cov(prof_m).severity == "info"
    assert "for mapping (30x warn" in cov(prof_m).explain
    with pytest.raises(RuleError):
        prof_a.set_purpose("nonsense")


def test_cli_purpose(tmp_path):
    out = tmp_path / "o"
    assert main(["interpret", str(D / "raw"), "-o", str(out), "--genome-size", "4e5", "--purpose", "mapping"]) == 0
    d = _json.loads((out / "fqci_report.json").read_text())
    cov = next(f for f in d["samples"][0]["findings"] if f["rule_id"] == "COVERAGE")
    assert cov["severity"] == "info"  # 48x expected >= 30x mapping threshold
    assert "purpose: **mapping**" in (out / "fqci_report.md").read_text()


def _gc_pair(sample, gc, n):
    pair = _pair(sample, n=n)
    for rep in pair.values():
        rep.modules["Basic Statistics"].rows = [r if r[0] != "%GC" else ["%GC", str(gc)]
                                                for r in rep.modules["Basic Statistics"].rows]
    return pair


def test_low_input_gc_bias():
    prof = load_profile("wgs_isolate")
    hit = interpret_sample("s", _gc_pair("s", 62, 156701), prof, genome_size=4.41e6, expected_gc=65.6)
    assert any(f.rule_id == "LOW_INPUT_GC_BIAS" and f.severity == "warn" for f in hit.findings)
    deep = interpret_sample("s", _gc_pair("s", 62, 1500000), prof, genome_size=4.41e6, expected_gc=65.6)
    assert not any(f.rule_id == "LOW_INPUT_GC_BIAS" for f in deep.findings)
    ongc = interpret_sample("s", _gc_pair("s", 65, 156701), prof, genome_size=4.41e6, expected_gc=65.6)
    assert not any(f.rule_id == "LOW_INPUT_GC_BIAS" for f in ongc.findings)
    nogc = interpret_sample("s", _gc_pair("s", 62, 156701), prof, genome_size=4.41e6)
    assert not any(f.rule_id == "LOW_INPUT_GC_BIAS" for f in nogc.findings)


def _set_reads(pair, n):
    for rep in pair.values():
        rep.modules["Basic Statistics"].rows = [r if r[0] != "Total Sequences" else ["Total Sequences", str(n)]
                                                for r in rep.modules["Basic Statistics"].rows]
    return pair


def test_compare_flags_excess_loss(tmp_path):
    before = list(_pair("s1", n=1000000).values()) + list(_pair("s2", n=1000000).values())
    after = (list(_set_reads(_pair("s1", end_bias=False), 590000).values())
             + list(_set_reads(_pair("s2", end_bias=False), 970000).values()))
    fd = tmp_path / "fp"
    fd.mkdir()
    (fd / "s1.fastp.json").write_text(_json.dumps({
        "command": "fastp -i a -I b --cut_right --length_required 76",
        "filtering_result": {"passed_filter_reads": 1180000, "low_quality_reads": 0,
                             "too_many_N_reads": 0, "too_short_reads": 820000, "too_long_reads": 0}}))
    md = compare(before, after, load_profile("wgs_isolate"), fastp_dir=fd)
    assert "| s1 | 59.0 | 0.0 | 41.0 | +41.0 | TOO MUCH REMOVED |" in md
    assert "| s2 | 97.0 | 0.0 | 3.0 | +3.0 | ok |" in md
    assert "Trimming removed more than expected in 1 of 2 samples" in md
    assert "--cut_right" in md and "too short 41.0%" in md
    assert "no s2.fastp.json" in md
    # over-trimmed s1 must not be counted as fixed
    assert "1 of 2 have nothing left to fix (over-trimmed samples are not counted as fixed)." in md
    s1 = md.split("## s1")[1].split("## s2")[0]
    assert "Nothing left to fix" not in s1 and "re-trim before using these reads" in s1


def test_no_advice_when_loss_is_expected(raw, trimmed):
    md = compare(raw, trimmed, load_profile("wgs_isolate"), fastp_dir=D)
    assert "fastp removed:" in md and "--cut_right" not in md


def test_compare_dimers_are_expected_loss(raw, trimmed):
    md = compare(raw, trimmed, load_profile("wgs_isolate"))
    # sim: 33% dimer pairs; 66.7% kept -> excess within allowance
    row = next(line for line in md.splitlines() if line.startswith("| sim |"))
    assert row.endswith("| ok |")


def test_fastp_cause_real_json():
    from fqci.compare import fastp_cause
    js = _json.loads((D / "sim.fastp.json").read_text())
    msg = fastp_cause(js)
    assert msg.startswith("fastp removed:") and "%" in msg


def test_run_summary_batch_and_runwide():
    from fqci.report import run_summary
    prof = load_profile("wgs_isolate")
    res = [interpret_sample(f"m{i}", _pair(f"m{i}"), prof, genome_size=4.41e6) for i in range(4)]
    L = "\n".join(run_summary(res))
    assert "4 samples" in L and "Composition bias in the last cycles only: 4 of 4" in L
    assert "Lowest expected depth" in L


def test_run_summary_dimer_batch(raw):
    from fqci.report import run_summary
    prof = load_profile("wgs_isolate")
    res = [interpret_sample("sim", group(raw)["sim"], prof)]
    res += [interpret_sample(f"m{i}", _pair(f"m{i}"), prof) for i in range(3)]
    L = "\n".join(run_summary(res))
    assert "Adapter dimers split the run in two: 1 sample at 33% and 3 below 1%." in L
    assert "in between" not in L


def test_run_summary_dimer_split_names_samples_in_between():
    """Real Nepal/Haiti cholera run: 12 samples at 18-43%, Vc13 at 1.3%, 14 below 1%
    (42.53: the real Vc1 value, shown as 42.5 elsewhere).
    The sentence used to drop Vc13 (12 + 14 = 26 of 27)."""
    from fqci.report import run_summary
    prof = load_profile("wgs_isolate")
    hi = [42.53, 31.2, 17.8, 22.5, 29.1, 23.5, 30.3, 27.3, 41.3, 20.6, 35.7, 25.3]
    lo = [0.5, 0.2, 0.4, 0.6, 0.2, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    names = [f"Vc{i}_Nepal" for i in range(1, 25)] + ["Vc_1786_Haiti_SE", "Vc_1792_Haiti_SE", "Vc_1798_Haiti_SE"]
    res = []
    for name, d in zip(names, hi + [1.3] + lo):
        r = interpret_sample(name, _pair(name), prof)
        r.metrics["pair"]["pair_dimer_pct"] = d
        res.append(r)
    L = "\n".join(run_summary(res))
    assert ("Adapter dimers split the run in two: 12 samples at 18-43% and 14 below 1%; "
            "1 in between (Vc13_Nepal 1.3%). That pattern usually means a library-prep batch") in L
    from fqci.ai_summary import build_payload
    payload, _ = build_payload(res, prof, 4.03e6, 47.5, anonymise=True)
    line = next(x for x in payload["run_summary_from_rules"] if "split the run" in x)
    assert "1 in between (S13 1.3%)" in line and "Vc" not in line
    # many in between: list five, count the rest
    for r, d in zip(res[12:20], [1.5, 2, 3, 4, 5, 6, 7, 8]):
        r.metrics["pair"]["pair_dimer_pct"] = d
    L = "\n".join(run_summary(res))
    assert "7 below 1%; 8 in between (Vc13_Nepal 1.5%, Vc14_Nepal 2.0%, Vc15_Nepal 3.0%, " in L
    assert "Vc17_Nepal 5.0% and 3 more)." in L


def test_prompt_asks_for_secondary_rule_ids():
    from fqci.ai_summary import load_prompt
    text, _, _ = load_prompt()
    assert "including secondary findings" in text


# ================================================================ code-review fixes


def _write_zip(d, fname, text):
    import zipfile
    stem = fname.replace(".fastq.gz", "").replace(".fq.gz", "")
    with zipfile.ZipFile(d / f"{stem}_fastqc.zip", "w") as z:
        z.writestr(f"{stem}_fastqc/fastqc_data.txt", text)


def _txt(rep):
    return _to_text(rep)


# --- generated commands
def test_output_never_overwrites_input(trimmed):
    res = interpret_sample("sim", group(trimmed)["sim"], load_profile("wgs_isolate"))
    cmd = res.fastp_command
    assert "-i sim_R1.trim.fastq.gz -o sim_R1.retrim.fastq.gz" in cmd
    assert "-I sim_R2.trim.fastq.gz -O sim_R2.retrim.fastq.gz" in cmd


def test_min_length_never_near_read_length():
    for prof, rl in (("amplicon", 45), ("wgs_isolate", 36)):
        rep = _miseq151("s_R1_001.fastq.gz", read_len=rl)
        res = interpret_sample("s", {1: rep}, load_profile(prof))
        assert res.fastp_params["length_required"] <= int(0.8 * rl)
        assert "capped at 80% of the read length" in res.param_reasons["length_required"]


@pytest.mark.parametrize("names,expect", [
    (("s1_r1.fastq.gz", "s1_r2.fastq.gz"), {"s1": [1, 2]}),      # lowercase r
    (("S1_R1_001.fastq.gz", "S1_R2_001.fastq.gz"), {"S1": [1, 2]}),
])
def test_pairing_names(names, expect):
    reps = [_miseq151(n) for n in names]
    assert {k: sorted(v) for k, v in group(reps).items()} == expect


def test_bare_digits_with_different_counts_are_not_mates():
    reps = [_miseq151("iso_1.fastq.gz", n=500000), _miseq151("iso_2.fastq.gz", n=730000)]
    g = group(reps)
    assert sorted(g) == ["iso_1", "iso_2"] and all(list(v) == [1] for v in g.values())
    same = group([_miseq151("ERR1_1.fastq.gz", n=5), _miseq151("ERR1_2.fastq.gz", n=5)])
    assert list(same) == ["ERR1"] and sorted(same["ERR1"]) == [1, 2]


def test_truseq_dimers_win_over_adapter_columns(raw):
    reps = group(raw)["sim"]
    for r in reps.values():  # Nextera column slightly above Illumina, as can happen
        ac = r.modules["Adapter Content"]
        i_ill, i_nx = ac.header.index("Illumina Universal Adapter"), ac.header.index("Nextera Transposase Sequence")
        for row in ac.rows:
            row[i_ill], row[i_nx] = "0.0", "0.05"
    res = interpret_sample("sim", reps, load_profile("wgs_isolate"))
    assert res.fastp_params["adapter_sequence"].startswith("AGATCGGAAGAGC")


def test_nextera_family_selected():
    rep = _miseq151("n_R1_001.fastq.gz")
    rep.modules["Overrepresented sequences"].header = ["Sequence", "Count", "Percentage", "Possible Source"]
    rep.modules["Overrepresented sequences"].rows = [
        ["CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", "5000", "2.5", "No Hit"]]
    res = interpret_sample("n", {1: rep}, load_profile("wgs_isolate"))
    assert res.fastp_params["adapter_sequence"] == "CTGTCTCTTATACACATCT"


def test_odd_filenames_are_quoted():
    import shlex
    pair = {1: _miseq151("my sample;rm -rf $(x)_R1_001.fastq.gz"),
            2: _miseq151("my sample;rm -rf $(x)_R2_001.fastq.gz")}
    cmd = interpret_sample("my sample;rm -rf $(x)", pair, load_profile("wgs_isolate")).fastp_command
    args = shlex.split(cmd)
    assert args[args.index("-i") + 1] == "my sample;rm -rf $(x)_R1_001.fastq.gz"


def test_non_fastq_input_is_skipped():
    res = interpret_sample("b", {1: _miseq151("b.bam")}, load_profile("wgs_isolate"))
    assert res.fastp_command.startswith("# b: skipped")


def test_r2_only_keeps_its_read_number():
    res = interpret_sample("p", {2: _miseq151("p_R2_001.fastq.gz")}, load_profile("wgs_isolate"))
    assert "-i p_R2_001.fastq.gz -o p_R2.trim.fastq.gz" in res.fastp_command


def test_threads_option():
    res = interpret_sample("s", _pair("s"), load_profile("wgs_isolate"), threads=12)
    assert res.fastp_command.endswith("--thread 12")


def test_cut_q_override_reaches_command():
    prof = load_profile("wgs_isolate")
    prof.params["cut_q"] = 15
    res = interpret_sample("s", _pair("s"), prof)
    assert res.fastp_params["cut_mean_quality"] == 15
    assert load_profile("rnaseq").params["cut_q"] == 15


def test_default_only_options_not_emitted():
    for k in ("qualified_quality_phred", "unqualified_percent_limit", "n_base_limit"):
        for prof in ("wgs_isolate", "amplicon", "rnaseq", "shotgun_metagenome"):
            for r in load_profile(prof).rules:
                assert k not in (r.get("actions") or {}), (prof, r["id"], k)


# --- metrics and rules
def test_body_bias_at_145_149_is_not_last_cycle():
    rep = _miseq151("x_R1_001.fastq.gz", read_len=150, end_bias=False)
    rows = rep.modules["Per base sequence content"].rows
    i = [r[0] for r in rows].index("145-149")
    rows[i][1:] = ["50", "10", "10", "30"]
    m = metrics.compute(rep)
    assert m["comp_bias_body"] >= 20 and m["comp_bias_end"] < 10


def test_polyg_fires_without_polyg_column():
    rep = _miseq151("g_R1_001.fastq.gz", adapter_cols=False)
    rep.modules["Overrepresented sequences"].header = ["Sequence", "Count", "Percentage", "Possible Source"]
    rep.modules["Overrepresented sequences"].rows = [["G" * 50, "50000", "5.0", "No Hit"]]
    res = interpret_sample("g", {1: rep}, load_profile("wgs_isolate"))
    assert any(f.rule_id == "POLYG" for f in res.findings) and res.fastp_params.get("trim_poly_g") is True


def test_dup_gt100_uses_percentage_of_total():
    m = metrics.compute(_miseq151("d_R1_001.fastq.gz"))
    assert abs(m["pct_seq_dup_gt100"] - 0.7) < 1e-9  # '>100' 0.5 + '>1k' 0.2 of total


def test_real_report_dup_column():
    if not REAL:
        pytest.skip("set FQCI_REAL_DATA")
    reps = parser.discover([REAL])
    for r in reps:
        m = metrics.compute(r)
        assert 0 <= m["pct_seq_dup_gt100"] <= 100


def test_overrep_classes_never_negative():
    rep = _miseq151("o_R1_001.fastq.gz")
    rep.modules["Overrepresented sequences"].header = ["Sequence", "Count", "Percentage", "Possible Source"]
    rep.modules["Overrepresented sequences"].rows = [["G" * 50, "4000", "2.0", "PCR Primer Index 1"]]
    m = metrics.compute(rep)
    assert m["overrep_polyg_pct"] == 2.0 and m["overrep_adapter_pct"] == 0 and m["overrep_other_pct"] == 0


def test_strict_action_templates():
    from fqci.rules import Profile, run_rules
    prof = Profile("t", "", {}, [{"id": "X", "scope": "read", "levels": {"warn": "a > 0"},
                                  "actions": {"adapter_sequence": "{adapter_seq}"}}])
    with pytest.raises(RuleError, match="no value for"):
        run_rules(prof, {"a": 1}, "read")


def test_evaluator_three_valued_and_safe():
    assert evaluate("a > 1 or b > 1", {"a": None, "b": 2}) is True
    assert evaluate("a > 1 and b > 1", {"a": None, "b": 0}) is False
    assert evaluate("a / 0 > 1", {"a": 1}) is None
    for bad in ("'x' * 10", "a ** 2", "min(a, key=abs)", "lambda: 1", "a[0]", "(a := 1)"):
        with pytest.raises(RuleError):
            evaluate(bad, {"a": 1})


def test_all_rule_names_are_known_metrics():
    """Guard against a misspelt metric silently disabling a rule."""
    from fqci.rules import compile_expr
    pair = _pair("s")
    for prof_name in available_profiles():
        prof = load_profile(prof_name)
        res = interpret_sample("s", pair, prof, genome_size=4.41e6, expected_gc=65.6)
        read_ns = set(prof.params) | set(res.metrics["R1"]) | {
            "min_len", "adapter_r1", "adapter_r2", "adapter_family", "expected_gc", "no_genome_size",
            "unit", "purpose", "min_len_capped", "reads_per_start", "dedup_expected", "dedup_ratio"}
        pair_ns = read_ns | set(res.metrics["pair"]) | {f"r{k}_{x}" for k in (1, 2) for x in res.metrics["R1"]}
        for r in prof.rules:
            known = pair_ns if r.get("scope") == "pair" else read_ns
            for e in r.get("levels", {}).values():
                missing = compile_expr(e)[1] - known
                assert not missing, (prof_name, r["id"], missing)


@pytest.mark.parametrize("rule,mutate", [
    ("ENCODING", lambda r: _set_basic(r, "Encoding", "Illumina 1.5")),
    ("N_CONTENT", lambda r: _set_rows(r, "Per base N content", [["1", "0.0"], ["50", "8.0"]])),
    ("POLYA", lambda r: _set_adapter_last(r, "PolyA", "7.0")),
    ("POLYG", lambda r: _set_adapter_last(r, "PolyG", "3.0")),
    ("ADAPTER_READTHROUGH", lambda r: _set_adapter_last(r, "Illumina Universal Adapter", "12.0")),
])
def test_rule_fires(rule, mutate):
    rep = _miseq151("f_R1_001.fastq.gz")
    mutate(rep)
    res = interpret_sample("f", {1: rep}, load_profile("wgs_isolate"))
    assert any(f.rule_id == rule for f in res.findings), [f.rule_id for f in res.findings]


def test_pair_count_mismatch_fires():
    pair = _pair("s")
    _set_basic(pair[2], "Total Sequences", "1000")
    res = interpret_sample("s", pair, load_profile("wgs_isolate"))
    assert any(f.rule_id == "PAIR_COUNT_MISMATCH" and f.severity == "fail" for f in res.findings)


def _set_basic(rep, key, val):
    rep.modules["Basic Statistics"].rows = [r if r[0] != key else [key, val] for r in rep.modules["Basic Statistics"].rows]


def _set_rows(rep, module, rows):
    rep.modules[module].rows = rows


def _set_adapter_last(rep, col, val):
    ac = rep.modules["Adapter Content"]
    ac.rows[-1][ac.header.index(col)] = val


# --- files, profiles, CLI
def test_missing_filename_is_an_error():
    txt = _txt(_miseq151("x_R1_001.fastq.gz")).replace("Filename\tx_R1_001.fastq.gz\n", "")
    with pytest.raises(ValueError, match="no 'Filename'"):
        parser.parse_text(txt, "x")


def test_duplicate_filenames_across_runs_raise(tmp_path):
    for run in ("run1", "run2"):
        (tmp_path / run).mkdir()
        _write_zip(tmp_path / run, "s_R1_001.fastq.gz", _txt(_miseq151("s_R1_001.fastq.gz")))
    with pytest.raises(ValueError, match="two FastQC reports for the same file name"):
        parser.discover([tmp_path / "run1", tmp_path / "run2"])


def test_profile_resolution_ignores_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "base").write_text("garbage: [")
    (tmp_path / "amplicon").mkdir()
    assert load_profile("rnaseq").name == "rnaseq" and load_profile("amplicon").name == "amplicon"


def test_profile_extends_relative_and_loops(tmp_path, monkeypatch):
    d = tmp_path / "profs"
    d.mkdir()
    (d / "parent.yaml").write_text("name: parent\nextends: wgs_isolate\nparams: {min_cov_warn: 77}\n")
    (d / "child.yaml").write_text("name: child\nextends: parent.yaml\n")
    monkeypatch.chdir(tmp_path)  # extends must resolve next to child.yaml, not in cwd
    child = load_profile(str(d / "child.yaml"))
    assert child.params["min_cov_warn"] == 77          # default purpose must not override it
    child.set_purpose("mapping")
    assert child.params["min_cov_warn"] == 30          # an explicit --purpose does
    (d / "loop.yaml").write_text("name: loop\nextends: loop.yaml\n")
    with pytest.raises(RuleError, match="loop"):
        load_profile(str(d / "loop.yaml"))


@pytest.mark.parametrize("args,msg", [
    (["interpret", "/nonexistent"], "no such file"),
    (["interpret", str(D / "raw"), "-p", "nosuch"], "profile not found"),
    (["interpret", str(D / "raw"), "--genome-size", "4.03"], "looks like Mb"),
    (["interpret", str(D / "raw"), "--purpose", "x"], "unknown purpose"),
])
def test_cli_errors_exit_2(args, msg, capsys, tmp_path):
    assert main(args + ["-o", str(tmp_path / "o")]) == 2
    assert msg in capsys.readouterr().err


def test_cli_set_type_error(tmp_path):
    with pytest.raises(SystemExit, match="needs a number"):
        main(["interpret", str(D / "raw"), "-o", str(tmp_path / "o"), "--set", "cut_q=abc"])


def test_compare_partial_pair(tmp_path):
    before = list(_pair("s").values())
    after = [_set_reads(_pair("s", end_bias=False), 980000)[2]]  # only R2 trimmed
    md = compare(before, after, load_profile("wgs_isolate"))
    assert "| s | 97.3 |" in md  # 980,000 of 1,007,380, from the R2 present in both sets


def test_python_m_fqci():
    out = _sp.run([_sys.executable, "-m", "fqci", "--version"], capture_output=True, text=True)
    assert out.returncode == 0 and out.stdout.strip() == __import__("fqci").__version__


# ================================================================ HTML report
def _html_of(results, prof=None):
    from fqci.html_report import to_html
    return to_html(results, prof or load_profile("wgs_isolate"), "t")


def test_html_written_by_cli(tmp_path):
    out = tmp_path / "o"
    assert main(["interpret", str(D / "raw"), "-o", str(out), "--genome-size", "4e5"]) == 0
    page = (out / "fqci_report.html").read_text(encoding="utf-8")
    assert page.startswith("<!DOCTYPE html>") and "sim" in page and "Adapter dimers" in page


def test_html_is_well_formed_and_offline(raw):
    from html.parser import HTMLParser
    res = [interpret_sample("sim", group(raw)["sim"], load_profile("wgs_isolate"), genome_size=4e5)]
    page = _html_of(res)

    class P(HTMLParser):
        def __init__(self):
            super().__init__()
            self.stack, self.errors, self.ext = [], [], []

        def handle_starttag(self, tag, attrs):
            if tag not in ("meta", "br", "hr", "img", "input", "link", "path", "line"):
                self.stack.append(tag)
            for k, v in attrs:
                if k in ("src", "href") and v and v.startswith(("http:", "https:", "//")):
                    self.ext.append(v)

        def handle_endtag(self, tag):
            if tag in ("path", "line"):
                return
            if not self.stack or self.stack.pop() != tag:
                self.errors.append(tag)
    p = P()
    p.feed(page)
    assert not p.errors and not p.stack, (p.errors, p.stack[-5:])
    assert not p.ext  # no external resources: works offline and by email


def test_html_escapes_sample_names():
    evil = '<script>alert(1)</script>"x'
    pair = {1: _miseq151(f"{evil}_R1_001.fastq.gz"), 2: _miseq151(f"{evil}_R2_001.fastq.gz")}
    page = _html_of([interpret_sample(evil, pair, load_profile("wgs_isolate"))])
    assert "<script>alert(1)" not in page and "&lt;script&gt;alert(1)" in page


def test_html_fixable_tile_counts_only_fixable():
    import re as _re
    prof = load_profile("wgs_isolate")
    res = [interpret_sample("lo", _pair("lo", n=100000), prof, genome_size=4.41e6),      # depth: not fixable
           interpret_sample("ok", _pair("ok", n=1500000), prof, genome_size=4.41e6)]     # end bias only: fixable
    assert res[0].verdict.startswith("FAIL; not fixable") and res[1].verdict.endswith(": fixable by trimming")
    page = _html_of(res, prof)
    n = _re.search(r'<div class="n">(\d+)</div><div class="l">fixable by trimming', page).group(1)
    assert n == "1"


def test_html_single_end_and_undetermined():
    prof = load_profile("wgs_isolate")
    res = [interpret_sample("se", {1: _miseq151("se.fastq.gz")}, prof, genome_size=4.41e6),
           interpret_sample("Undetermined_S0_L001", _pair("Undetermined_S0_L001"), prof)]
    page = _html_of(res, prof)
    assert "Reads the demultiplexer could not assign" in page and "se" in page


def test_html_findings_controls():
    prof = load_profile("wgs_isolate")
    res = [interpret_sample("lo", _pair("lo", n=100000), prof, genome_size=4.41e6),
           interpret_sample("ok", _pair("ok", n=1500000, end_bias=False), prof, genome_size=4.41e6)]
    page = _html_of(res, prof)
    assert 'data-act="expand"' in page and 'data-act="collapse"' in page
    assert 'data-filter="ALL"' in page and 'data-filter="FAIL"' in page and 'data-filter="PASS"' in page
    assert 'data-filter="WARN"' not in page          # only verdicts present in the run get a chip
    assert 'data-level="FAIL"' in page and 'data-level="PASS"' in page


# ================================================================ AI summary (no network: fake provider)
class FakeProvider:
    name, model = "fake", "fake-model-1"

    def __init__(self, reply=None):
        self.reply, self.calls = reply, []

    def generate(self, system, user, schema, max_tokens):
        self.calls.append({"system": system, "user": user, "schema": schema, "max_tokens": max_tokens})
        payload = _json.loads(user.split("\n\n", 1)[1])
        if self.reply is not None:
            return (self.reply(payload) if callable(self.reply) else self.reply), {"model": self.model}
        samples = [{"id": s["id"], "summary": f"{s['id']} verdict is {s['verdict']}.",
                    "rule_ids": [f["rule_id"] for f in s["findings"] if f["category"] == "primary"][:2],
                    "next_steps": ["Run the recommended fastp command."]} for s in payload["samples"]]
        return _json.dumps({"run_summary": "One run.", "samples": samples}), {"model": self.model,
                                                                               "input_tokens": 1, "output_tokens": 1}


def _two_samples():
    prof = load_profile("wgs_isolate")
    return [interpret_sample("lowdepth_S1", _pair("lowdepth_S1", n=100000), prof, genome_size=4.41e6),
            interpret_sample("clean_S2", _pair("clean_S2", n=1500000, end_bias=False), prof, genome_size=4.41e6)], prof


def test_ai_prompt_ships_and_is_versioned():
    text, name, h = A.load_prompt()
    assert name == "findings_v1" and len(h) == 12 and "authoritative" in text


def test_ai_payload_contains_findings_not_files():
    res, prof = _two_samples()
    payload, ids = A.build_payload(res, prof, 4.41e6, None)
    s = payload["samples"][0]
    assert s["id"] == "lowdepth_S1" and s["findings"] and "fastp_options" in s
    assert "_R1_001.fastq.gz" not in _json.dumps(payload)  # file names are never sent


def test_ai_anonymise_hides_names_everywhere():
    res, prof = _two_samples()
    payload, ids = A.build_payload(res, prof, 4.41e6, None, anonymise=True)
    blob = _json.dumps(payload)
    assert "lowdepth_S1" not in blob and "clean_S2" not in blob
    assert ids == {"S01": "lowdepth_S1", "S02": "clean_S2"}


def test_ai_summarise_with_fake_provider_and_deanonymise():
    res, prof = _two_samples()
    fp = FakeProvider()
    ai, payload = A.summarise(res, prof, fp, genome_size=4.41e6, anonymise=True)
    assert set(ai.samples) == {"lowdepth_S1", "clean_S2"}
    assert ai.samples["lowdepth_S1"].summary.startswith("lowdepth_S1 verdict")  # names restored locally
    assert "lowdepth_S1" not in fp.calls[0]["user"]                             # but never sent
    assert ai.record["prompt"] == "findings_v1" and ai.record["anonymised"] is True
    assert fp.calls[0]["schema"] == A.OUTPUT_SCHEMA


def test_ai_validation_flags_and_filters():
    res, prof = _two_samples()
    payload, ids = A.build_payload(res, prof, 4.41e6, None)
    reply = _json.dumps({"run_summary": "Depth was 999x overall.", "samples": [
        {"id": "lowdepth_S1", "summary": "Only 4 x coverage, about 12345 reads.", "rule_ids": ["COVERAGE", "MADE_UP"],
         "next_steps": ["a", "b", "c", "d"]},
        {"id": "ghost", "summary": "x", "rule_ids": [], "next_steps": []}]})
    rs, samples, warns = A.validate(reply, payload, ids)
    s = samples["lowdepth_S1"]
    assert s.rule_ids == ["COVERAGE"] and len(s.next_steps) == 3
    assert any("MADE_UP" in w for w in s.warnings) and any("12345" in w for w in s.warnings)
    assert any("unknown sample id 'ghost'" in w for w in warns)
    assert any("no AI summary returned for: clean_S2" in w for w in warns)
    assert any("999" in w for w in warns)


def test_ai_numbers_from_input_are_accepted():
    res, prof = _two_samples()
    payload, ids = A.build_payload(res, prof, 4.41e6, None)
    depth = payload["samples"][1]["metrics"]["expected_depth_x"]
    reads = payload["samples"][1]["metrics"]["reads_or_pairs"]
    reply = _json.dumps({"run_summary": "ok", "samples": [
        {"id": "clean_S2", "summary": f"About {depth:.0f}x from {reads:,} pairs; 2 of 2 checks.", "rule_ids": [],
         "next_steps": []}, {"id": "lowdepth_S1", "summary": "Low.", "rule_ids": [], "next_steps": []}]})
    _, samples, warns = A.validate(reply, payload, ids)
    assert samples["clean_S2"].warnings == [] and not any("numbers" in w for w in warns)


def test_ai_bad_reply_raises():
    res, prof = _two_samples()
    for bad in ("not json", _json.dumps({"samples": "x"})):
        with pytest.raises(A.AIError):
            A.summarise(res, prof, FakeProvider(reply=bad))


def test_ai_html_and_markdown_sections():
    from fqci.html_report import to_html
    res, prof = _two_samples()
    ai, _ = A.summarise(res, prof, FakeProvider(reply=lambda p: _json.dumps({
        "run_summary": "<b>run</b>", "samples": [{"id": s["id"], "summary": "<script>x</script>",
                                                  "rule_ids": [], "next_steps": []} for s in p["samples"]]})))
    page = to_html(res, prof, "t", ai=ai)
    assert page.count("AI summary · generated, advisory") == 3   # run box + one per sample
    assert "<script>x</script>" not in page and "&lt;script&gt;x" in page
    assert "fake / fake-model-1" in page
    md = A.to_markdown(ai)
    assert md.startswith("## AI summary (generated; verify") and "### lowdepth_S1" in md
    page2 = to_html(res, prof, "t", ai_error="ANTHROPIC_API_KEY is not set")
    assert "Requested but not generated: ANTHROPIC_API_KEY is not set" in page2


def test_ai_cli_without_key_still_writes_reports(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    out = tmp_path / "o"
    assert main(["interpret", str(D / "raw"), "-o", str(out), "--ai-summary"]) == 0
    assert "AI summary not generated" in capsys.readouterr().err
    assert (out / "fqci_report.html").exists() and not (out / "fqci_ai_summary.json").exists()
    assert "Requested but not generated" in (out / "fqci_report.md").read_text()


def test_ai_cli_with_fake_provider(tmp_path, monkeypatch):
    monkeypatch.setitem(A.PROVIDERS, "fake", lambda model=None: FakeProvider())
    out = tmp_path / "o"
    assert main(["interpret", str(D / "raw"), "-o", str(out), "--ai-summary", "--ai-provider", "fake",
                 "--anonymise"]) == 0
    sent = _json.loads((out / "fqci_ai_request.json").read_text())
    assert sent["samples"][0]["id"] == "S01"
    summary = _json.loads((out / "fqci_ai_summary.json").read_text())
    assert summary["record"]["provider"] == "fake" and "sim" in summary["samples"]
    assert "AI summary · generated, advisory" in (out / "fqci_report.html").read_text()


def test_ai_unknown_provider():
    with pytest.raises(A.AIError, match="unknown AI provider"):
        A.get_provider("nope")


# ---- real Anthropic SDK against a local stand-in server (checks the request we send)
def test_anthropic_provider_request_shape(monkeypatch):
    pytest.importorskip("anthropic")
    import http.server
    import threading
    seen = {}

    class H(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_POST(self):
            n = int(self.headers["content-length"])
            auth = self.headers.get("authorization", "")  # newer SDKs: Bearer; older: x-api-key
            seen["path"] = self.path
            seen["key"] = self.headers.get("x-api-key") or (auth[7:] if auth.startswith("Bearer ") else None)
            seen["body"] = _json.loads(self.rfile.read(n))
            p = _json.loads(seen["body"]["messages"][0]["content"].split("\n\n", 1)[1])
            reply = {"run_summary": "ok", "samples": [{"id": s["id"], "summary": "fine", "rule_ids": [],
                                                       "next_steps": []} for s in p["samples"]]}
            body = _json.dumps({"id": "msg_1", "type": "message", "role": "assistant", "model": "claude-sonnet-5-5",
                                "content": [{"type": "text", "text": _json.dumps(reply)}],
                                "stop_reason": "end_turn", "stop_sequence": None,
                                "usage": {"input_tokens": 123, "output_tokens": 45}}).encode()
            self.send_response(200)
            self.send_header("content-type", "application/json")
            self.send_header("content-length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    srv = http.server.HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")
    monkeypatch.setenv("ANTHROPIC_BASE_URL", f"http://127.0.0.1:{srv.server_port}")
    for k in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy", "ALL_PROXY", "all_proxy"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("NO_PROXY", "127.0.0.1")
    try:
        res, prof = _two_samples()
        ai, _ = A.summarise(res, prof, A.get_provider("anthropic"))
    finally:
        srv.shutdown()
    b = seen["body"]
    assert seen["path"] == "/v1/messages" and seen["key"] == "test-key"
    assert b["model"] == "claude-sonnet-5-5" and "authoritative" in b["system"]
    assert b["output_config"]["format"]["type"] == "json_schema"
    assert b["output_config"]["format"]["schema"] == A.OUTPUT_SCHEMA
    assert "temperature" not in b  # Sonnet 5 models reject non-default sampling settings
    assert ai.record["input_tokens"] == 123 and set(ai.samples) == {"lowdepth_S1", "clean_S2"}


def test_anthropic_provider_stop_reason(monkeypatch):
    pytest.importorskip("anthropic")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "k")
    prov = A.AnthropicProvider()

    class R:
        stop_reason, content, model, usage = "max_tokens", [], "m", None
    monkeypatch.setattr(prov._client.messages, "create", lambda **kw: R())
    with pytest.raises(A.AIError, match="stopped early"):
        prov.generate("s", "u", {}, 10)


# ---- manual route: --ai-export / --ai-import (no API key)
def _fake_reply_for(payload, fence=True):
    body = _json.dumps({"run_summary": "Run ok.", "samples": [
        {"id": s["id"], "summary": f"{s['id']}: {s['verdict']}.", "rule_ids": [], "next_steps": []}
        for s in payload["samples"]]})
    return f"Here is the JSON:\n```json\n{body}\n```\nHope this helps." if fence else body


def test_ai_export_writes_paste_ready_prompt(tmp_path):
    res, prof = _two_samples()
    prep = A.prepare(res, prof, genome_size=4.41e6, anonymise=True)
    prompt_file, export_file = A.export(prep, tmp_path)
    text = Path(prompt_file).read_text(encoding="utf-8")
    assert "authoritative" in text and '"run_summary"' in text and "return ONLY a JSON object" in text
    assert "lowdepth_S1" not in text and '"id": "S01"' in text          # anonymised in what is pasted
    exp = _json.loads(Path(export_file).read_text(encoding="utf-8"))
    assert exp["ids"] == {"S01": "lowdepth_S1", "S02": "clean_S2"}      # mapping stays local


def test_extract_json_from_chat_reply():
    assert A.extract_json('```json\n{"a": 1}\n```') == '{"a": 1}'
    assert A.extract_json('Sure! {"a": {"b": 2}} Done.') == '{"a": {"b": 2}}'
    with pytest.raises(A.AIError):
        A.extract_json("no json here")


def test_ai_import_roundtrip_and_mismatch(tmp_path):
    res, prof = _two_samples()
    prep = A.prepare(res, prof, genome_size=4.41e6, anonymise=True)
    A.export(prep, tmp_path)
    reply = tmp_path / "reply.txt"
    reply.write_text(_fake_reply_for(prep.payload))
    ai = A.import_reply(str(reply), str(tmp_path / A.EXPORT_FILE), prep, model="Claude (chat)")
    assert ai.record["provider"] == "manual (chat)" and ai.record["model"] == "Claude (chat)"
    assert ai.samples["lowdepth_S1"].summary.startswith("lowdepth_S1:")  # names restored
    other = A.prepare(res[:1], prof, genome_size=4.41e6, anonymise=True)  # different run
    with pytest.raises(A.AIError, match="different data"):
        A.import_reply(str(reply), str(tmp_path / A.EXPORT_FILE), other)
    with pytest.raises(A.AIError, match="run with --ai-export first"):
        A.import_reply(str(reply), str(tmp_path / "missing.json"), prep)


def test_ai_export_import_cli(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    out = tmp_path / "o"
    base = ["interpret", str(D / "raw"), "-o", str(out), "--genome-size", "4e5"]
    assert main(base + ["--ai-export", "--anonymise"]) == 0
    assert "Paste it into a Claude chat" in capsys.readouterr().err
    exp = _json.loads((out / A.EXPORT_FILE).read_text())
    reply = tmp_path / "reply.json"
    reply.write_text(_fake_reply_for(exp["payload"]))
    assert main(base + ["--ai-import", str(reply)]) == 0                 # --anonymise read from the export
    page = (out / "fqci_report.html").read_text()
    assert "manual (chat)" in page and "AI summary · generated, advisory" in page
    bad = tmp_path / "bad.json"
    bad.write_text("I could not do this.")
    assert main(base + ["--ai-import", str(bad)]) == 2
    assert "no JSON object found" in capsys.readouterr().err


def test_ai_routes_are_exclusive(tmp_path):
    with pytest.raises(SystemExit):
        main(["interpret", str(D / "raw"), "-o", str(tmp_path), "--ai-summary", "--ai-export"])


def test_ai_numbers_checked_per_sample_and_scaled():
    res, prof = _two_samples()
    payload, ids = A.build_payload(res, prof, 4.41e6, None)
    own = A._numbers_in(payload["samples"][0]) | A._numbers_in({k: v for k, v in payload.items() if k != "samples"})
    other_depth = next(v for v in payload["samples"][1]["metrics"].values()   # a figure only clean_S2 has
                       if isinstance(v, float) and f"{v:g}" not in own and v < 1e5)
    own_reads = payload["samples"][0]["metrics"]["reads_or_pairs"]           # belongs to lowdepth_S1
    reply = _json.dumps({"run_summary": f"Depths up to {other_depth}x.", "samples": [
        {"id": "lowdepth_S1", "summary": f"Depth {other_depth}x; about {own_reads / 1e6:.1f} million reads; "
                                          "roughly 1.9 million pairs lost.", "rule_ids": [], "next_steps": []},
        {"id": "clean_S2", "summary": "Fine.", "rule_ids": [], "next_steps": []}]})
    _, samples, warns = A.validate(reply, payload, ids)
    w = " ".join(samples["lowdepth_S1"].warnings)
    assert str(other_depth) in w             # another sample's number is not accepted for this one
    assert "1,900,000" in w                  # '1.9 million' checked at full value
    assert f"{own_reads:,}" not in w and "100,000" not in w   # own count, rounded, is accepted
    assert not any("run summary" in x for x in warns)          # run summary may quote any sample

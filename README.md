# fqci: profile-aware FastQC interpretation

`fqci` reads FastQC output (`*_fastqc.zip` or extracted folders), computes quantitative metrics
from each module, applies library-type-specific rules written in YAML, and produces:

- a Markdown report that separates primary findings, findings that are probably secondary to a
  primary problem, and findings that are expected for the library type;
- a fastp command for each sample, with the reason for every parameter;
- JSON (full metrics and findings) plus `fastp_params.json`, for use in pipelines;
- a before/after comparison once you have re-run FastQC on the trimmed reads, which also checks
  that trimming did not remove more reads than expected;
- a run-level summary separating run-wide or batch problems from sample-specific ones.

It interprets the numbers behind FastQC's PASS/WARN/FAIL flags instead of the flags themselves.
Those flags are calibrated for random-shotgun genomic libraries.

## Install

```bash
git clone https://github.com/dkomics/fqci.git && cd fqci
pip install .            # needs Python >= 3.9 and PyYAML
pip install '.[test]' && pytest -q
```

fqci reads FastQC reports; FastQC and fastp themselves are installed separately
(e.g. `conda install -c bioconda fastqc fastp`).

## Use

```bash
# 1. interpret raw-read QC
fqci interpret qc_raw/ -p wgs_isolate --purpose assembly --genome-size 4.03e6 --expected-gc 47.5 -o fqci_out

# 2. trim (run from the directory with the raw FASTQ; the file names come from the FastQC reports)
bash fqci_out/fastp_commands.sh

# 3. re-QC and close the loop
fastqc -o qc_trim *.trim.fastq.gz
fqci compare --before qc_raw/ --after qc_trim/ -p wgs_isolate --fastp-dir . -o fqci_compare.md

fqci profiles            # list available library profiles
```

Profiles: `wgs_isolate`, `shotgun_metagenome`, `amplicon`, `rnaseq`. Pass a path to your own
YAML file to `-p` to use a custom profile.

Useful options: `--purpose assembly|mapping` (depth thresholds: 50x/20x or 30x/15x warn/fail),
`--genome-size` in bp and `--expected-gc` (enable depth, duplication-vs-depth and GC checks),
`--set PARAM=VALUE` (override any profile parameter), `--threads` (fastp `--thread`), and for
`compare`, `--fastp-dir` (read `<sample>.fastp.json` to explain read losses).

## How interpretation works

1. **Metrics** (`metrics.py`): for example, the position from which the per-base lower quartile
   *stays* below Q20, the % of reads with mean quality below Q20, the % of tiles deviating by 5
   or more, composition bias at the read start, body and last cycle, the number of GC modes,
   the % of adapter-dimer reads (overrepresented sequences that *start* at the adapter junction,
   which FastQC's Adapter Content module largely misses), adapter onset, and the i7 index read
   out of dimers.
2. **Rules** (`profiles/base.yaml`): each rule has `levels` (fail/warn/info expressions over the
   metrics and profile parameters), a templated explanation, a likely cause, whether trimming can
   fix it, `actions` (fastp options), and `explains:` (rules it accounts for when it fails).
3. **Profiles** extend `base`, override thresholds in `params`, and `patch` rules. For example,
   `amplicon` marks duplication as expected and drops sliding-window trimming.
4. **Secondary marking**: when a primary rule fails (adapter dimers, for instance), the findings
   it explains are listed as provisional, and the verdict counts only primary findings. The
   before/after comparison then shows whether they resolved.
5. **Depth**: with `--genome-size`, expected median depth is the read-count estimate x 0.8
   (x 0.7 below 30x), and duplication is judged against what depth alone would produce, so deep
   runs of small genomes are not flagged for unavoidable duplication.

Rule expressions use a small language (names, numbers, comparisons, `+ - * /`, `and/or/not`,
`min/max/abs`), interpreted from the syntax tree without `eval`. A metric that could not be
computed (e.g. a column older FastQC versions do not write) is "unknown": comparisons with it
are unknown, and `unknown or true` is true, so one missing metric does not disable a rule.
Templated fastp options fail with an error if a value is missing.

## Profile keys

`name`, `extends` (a built-in profile name or a path relative to this file), `description`,
`params`, `base_actions` (fastp options always applied; values may be `{param}` templates),
`purposes` / `default_purpose`, `rules` (new rules), `patch` (replace keys of an inherited rule),
`disable` (list of inherited rule ids to drop), `notes`.

## Adding or changing a rule

```yaml
# my_profile.yaml
name: my_profile
extends: wgs_isolate
params:
  min_cov_warn: 80
patch:
  POLYG:
    levels: {fail: "polyg_pct >= 0.5"}
rules:
  - id: SHORT_READS
    scope: read
    title: Reads shorter than expected
    levels: {warn: "read_len < 100"}
    explain: "Read length {read_len} bp."
    fixable: false
```

## Validation

**Simulated data.** `tests/data/simulate.py` generates 101 bp TruSeq paired-end reads with 33%
adapter-dimer pairs (R1 dimers running off into Q2 poly-A, R2 dimer mates at low quality). The
tool estimates 32.9% dimers. The generated fastp command, run unchanged with fastp 1.4.0, keeps
every genomic pair and removes every dimer pair, and every secondary finding clears afterwards.

**Real data** (FastQC 0.12.1; fastp 0.22.0 and 1.4.0):

- *Vibrio cholerae*, 27 samples: 24 GAIIx 101 bp TruSeq paired-end and 3 Illumina single-end.
  Dimers of 18-42% were found in one library batch (12 samples) and 0-1.3% in the other, while
  FastQC's Adapter Content reported about 0.5% in the dimer samples checked. After trimming one
  sample with the generated command, every FastQC
  module passed and the R2 tile warning, which the dimers had caused, disappeared. Deep
  single-end samples (215-257x) were correctly not flagged for duplication explained by depth.
- *Mycobacterium tuberculosis*, 24 samples, MiSeq 2x151 Nextera-type. Last-cycle composition
  bias in every sample was identified as a run artefact and fixed by trimming one base.
  Comparing `--cut_right` with a 76 bp minimum against `--cut_tail` with 50 bp kept 59.0% and
  95.3% of pairs; `fqci compare` flags the first setting in all 24 samples and none with the
  second. The depth factor was calibrated against TB-Profiler's genome median depth (median
  ratio 0.78).

**Tests.** 93 tests cover the parser, metrics, the rule evaluator, every profile and rule, command
generation (quoting, paired/single-end, input/output names), the CLI and its error handling, and
`compare`. With fastp on `PATH` a generated command is run end to end. Set
`FQCI_REAL_DATA=/path/to/zips` to add tests on your own reports. An independent code review was
done before this release; its findings are fixed and tested.

## Limits

- Illumina short reads only. Nanopore (NanoPlot/pycoQC) input is not supported yet.
- FastQC's duplication and overrepresented-sequence modules use a subset of reads, so the dimer
  fraction is an estimate (within about 1% in simulation).
- Depth is estimated from read counts. Real depth over GC-rich or otherwise hard regions can be
  much lower than the median, especially in low-yield libraries; check depth over the regions
  your analysis depends on after mapping.
- Thresholds were tuned on the datasets above; check them on yours with `--set` and `compare`.

## Release history

See [CHANGELOG.md](CHANGELOG.md).

## Author and licence

Davis Kuchaka (ORCID [0000-0002-0991-703X](https://orcid.org/0000-0002-0991-703X)),
Kilimanjaro Clinical Research Institute, Moshi, Tanzania. Released under the MIT licence; see
[LICENSE](LICENSE). Bug reports and suggestions: https://github.com/dkomics/fqci/issues

# Changelog

This project follows [Semantic Versioning](https://semver.org/). Versions below 1.0 may change
the command-line interface, profile format or outputs between minor releases.

## Unreleased (0.3.0.dev0)

- Optional AI summary (`--ai-summary`): plain-language summaries per sample and for the run, from
  a fixed, versioned prompt (`findings_v1`), using the Anthropic API (`pip install "fqci[ai]"`).
  Advisory only: rule-based findings stay authoritative. Replies are schema-checked, cited rule
  IDs are verified, unverifiable numbers are flagged, `--anonymise` keeps sample names local, and
  the exact request is saved. Without the flag fqci makes no network calls.
- No API key needed: `--ai-export` writes a prompt to paste into a Claude chat; `--ai-import
  REPLY` checks the saved reply against the export and the current run and adds it to the
  reports. Number checks are per sample and read "million"/"thousand" at full value.

## 0.2.0 (2026-11-16)

- HTML report (`fqci_report.html`), written by `fqci interpret` next to the Markdown report:
  verdict tiles, run summary, depth and adapter-dimer charts with thresholds, a sortable sample
  table, and collapsible findings per sample with Expand all / Collapse all and a FAIL / WARN /
  PASS filter. Self-contained (no external resources), dark-mode aware, print-friendly; sample
  names and embedded data are escaped.
- Requires Python 3.10 or newer (3.9 reached end of life in October 2025). Tests run on 3.10-3.13.

## 0.1.0 (2026-11-16)

Initial public release.

- Interprets FastQC reports (`*_fastqc.zip` or extracted folders) for Illumina short reads using
  YAML rule profiles: `wgs_isolate`, `shotgun_metagenome`, `amplicon`, `rnaseq`.
- Separates primary, secondary (explained by another finding) and library-type-expected
  findings, and gives a verdict saying whether trimming resolves the primary problem.
- Detects adapter dimers from overrepresented sequences, which FastQC's Adapter Content module
  under-reports, judged per read pair; reads the sample index out of dimers.
- Generates fastp commands (`cut_tail` 3' trimming, explicit adapters, poly-G/poly-X, last-cycle
  trimming, capped minimum length) with the reason for each option; never overwrites inputs.
- Depth estimates calibrated against measured depth; duplication judged against depth;
  `--purpose assembly|mapping` depth thresholds; low-input GC-bias warning.
- `fqci compare`: before/after tables, a check that trimming removed no more than expected, and
  optional fastp JSON breakdown of where reads went.
- Run-level summaries separating run-wide, batch and sample-specific findings.
- Validated on simulated data and on real *V. cholerae* (27 samples) and *M. tuberculosis*
  (24 samples) runs; independent code review before release. See README, "Validation".

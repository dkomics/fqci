# Changelog

This project follows [Semantic Versioning](https://semver.org/). Versions below 1.0 may change
the command-line interface, profile format or outputs between minor releases.

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

# Changelog

This project follows [Semantic Versioning](https://semver.org/). Versions below 1.0 may change
the command-line interface, profile format or outputs between minor releases.

## 0.3.0 (2026-11-17)

- Optional AI summary (`--ai-summary`): plain-language summaries per sample and for the run, from
  a fixed, versioned prompt (`findings_v1`), using the Anthropic API (`pip install "fqci[ai]"`,
  `ANTHROPIC_API_KEY`; default model `claude-sonnet-5-5`). Advisory only: the rule-based findings,
  verdicts and fastp commands stay authoritative. Replies are schema-checked; cited rule IDs must
  exist in that sample's findings; numbers in a sample's summary must appear in that sample's data
  or the run-level figures (others are flagged); `--anonymise` replaces sample names with codes
  before anything is sent; the exact request is saved. Without an AI option fqci makes no network
  calls, and an API failure never stops the rule-based reports.
- No API key needed: `--ai-export` writes a prompt to paste into a Claude chat; `--ai-import
  REPLY` checks the saved reply against the export and the current run (same data, options and
  prompt) and adds it to the reports.
- Cost: after an API run fqci prints the tokens used and an estimated cost at list price, and adds
  the same line to the reports and `fqci_ai_summary.json`. Validated against the Claude Console:
  four runs estimated at about $0.55-0.60, billed $0.56.
- Run summary: the adapter-dimer split sentence now names samples between 1% and 10% (it used to
  leave them out, so the counts did not add up), and says "1 sample at 33%" instead of "1 samples
  at 33-33%".
- Validated on real data with both routes: a 27-sample *V. cholerae* run (Nepal paired-end and
  Haiti single-end) and a 24-sample *M. tuberculosis* MiSeq run. All automatic checks passed, and
  the summaries were checked by hand against the data sent. Known limit: rankings across samples
  written in prose ("the most pronounced in the run") cannot be checked automatically.

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

# core-unit-bioinformatics/exomeqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Added`

- Synthetic DRAGEN CNV test data (`tests/fixtures/cnvqc/`, with its generator script) and nf-tests for `DRAGEN_CNV_QC` (thresholds passed from the params, no thresholds, purity given in percent → `ERROR`, stub); more self-tests in `analyze_cnv_metrics.R --run-self-tests` (purity range and boundary, unknown options).
- Param `--total_excluded_bp_max` for the threshold on the total excluded bp, which `analyze_cnv_metrics.R` already supported but the pipeline could not set. Not checked if unset (default), so results are unchanged.
- Documentation of the CNV QC thresholds (`docs/usage.md`), of the CNV QC output (`docs/output.md`) and descriptions of the CNV QC params in `nextflow_schema.json`.
- New parameter `--multiqc_report_name` _(Implemented with support from Claude (Anthropic).)_: the file name of the MultiQC report without `.html`, e.g. the run ID. The report is written to `multiqc/<name>.html` (data in `<name>_data/`, plots in `<name>_plots/`). Validated against `^(?!.*\.html$)[A-Za-z0-9][A-Za-z0-9._-]*$`. Not set by default: the report keeps its current name (`multiqc_report.html`, or `<title>_multiqc_report.html` with `--multiqc_title`), so existing runs are unchanged. Report content is not affected.
- Pipeline nf-test `tests/multiqc_report_name.nf.test` (named report, default name, invalid name) with minimal fixtures in `tests/fixtures/multiqc_report_name/`.

### `Fixed`

- DRAGEN CNV QC ignored all thresholds: `modules/local/dragen_cnv_qc/main.nf` built `task.ext.args` from the threshold params in `conf/modules.config` but never passed it to `bin/analyze_cnv_metrics.R`. Only the hard-coded check (`EstimatedTumorPurity=NA` → `FAIL`) was active, so e.g. a sample with a purity below `--purity_min` still got `PASS`. The thresholds are now passed to the script.
  - **Results change:** with thresholds set, samples that violate them now get `QC_STATUS=FAIL` with the reason in `QC_FAIL_REASONS` (e.g. `EstimatedTumorPurity:0.25<0.3`) in `cnv-metrics/cnv-metrics.csv` and the MultiQC section "CNV QC Metrics". Samples within the thresholds, and runs without thresholds, are unchanged.
  - `conf/modules.config`: a threshold is passed if it is set (`!= null`). Before, a threshold of `0` was dropped silently.
  - `bin/analyze_cnv_metrics.R` stops with an `ERROR` on unknown options (e.g. `--purity-min` instead of `--purity_min`), which would otherwise leave that check silently disabled, and prints all errors as `ERROR: <message>` with exit status 1.
  - The stub of `DRAGEN_CNV_QC` now creates the same files as a real run with the default `ext.args` (`cnv-metrics.csv`, `purcov_maxmodel_lineplot.json`, `karyoplot-*_mqc.png`, with `ext.prefix`).
- `DRAGEN_CNV_QC` failed with a missing-output error if the DRAGEN output had no `*.cnv.purity.coverage.models.tsv` (no `*.json`) or no `*.tn.tsv.gz`/`*.tumor.baf.bedgraph.gz` (no `*.png`). These outputs are now optional, as documented in `docs/usage.md`. No change in results.
- `bin/analyze_cnv_metrics.R` was committed without the executable bit (mode 100644), so `DRAGEN_CNV_QC` failed with `Permission denied` (exit 126) when run from a fresh clone. It is now executable (100755). No change in results.

## v0.1dev - [date]

Initial release of core-unit-bioinformatics/exomeqc, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

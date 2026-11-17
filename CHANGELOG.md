# core-unit-bioinformatics/exomeqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Fixed`

- DRAGEN CNV QC ignored all thresholds _(Implemented with support from Claude (Anthropic).)_: `modules/local/dragen_cnv_qc/main.nf` built `task.ext.args` from the threshold params in `conf/modules.config` but never passed it to `bin/analyze_cnv_metrics.R`. Only the hard-coded check (`EstimatedTumorPurity=NA` → `FAIL`) was active, so e.g. a sample with a purity below `--purity_min` still got `PASS`. The thresholds are now passed to the script.
  - **Results change:** with thresholds set, samples that violate them now get `QC_STATUS=FAIL` with the reason in `QC_FAIL_REASONS` (e.g. `EstimatedTumorPurity:0.25<0.3`) in `cnv-metrics/cnv-metrics.csv` and the MultiQC section "CNV QC Metrics". Samples within the thresholds, and runs without thresholds, are unchanged. The production parameters (`pmad_max=0.2`, `outlier_baf_fraction_max=0.1`, `purity_min=0.3`) take effect for the first time.
  - `conf/modules.config`: a threshold is passed if it is set (`!= null`). Before, a threshold of `0` was dropped silently.
  - `bin/analyze_cnv_metrics.R` stops with an `ERROR` on unknown options (e.g. `--purity-min` instead of `--purity_min`), which would otherwise leave that check silently disabled, and prints all errors as `ERROR: <message>` with exit status 1.
  - The stub of `DRAGEN_CNV_QC` now creates the same files as a real run with the default `ext.args` (`cnv-metrics.csv`, `purcov_maxmodel_lineplot.json`, `karyoplot-*_mqc.png`, with `ext.prefix`).
- `DRAGEN_CNV_QC` failed with a missing-output error if the DRAGEN output had no `*.cnv.purity.coverage.models.tsv` (no `*.json`) or no `*.tn.tsv.gz`/`*.tumor.baf.bedgraph.gz` (no `*.png`). These outputs are now optional, as documented in `docs/usage.md`. No change in results.
- `bin/analyze_cnv_metrics.R` was committed without the executable bit (mode 100644), so `DRAGEN_CNV_QC` failed with `Permission denied` (exit 126) when run from a fresh clone. It is now executable (100755). No change in results.

### `Added`

- Synthetic DRAGEN CNV test data (`tests/fixtures/cnvqc/`, with its generator script) and nf-tests for `DRAGEN_CNV_QC` (thresholds passed from the params, no thresholds, purity given in percent → `ERROR`, stub); more self-tests in `analyze_cnv_metrics.R --run-self-tests` (purity range and boundary, unknown options).
- Documentation of the CNV QC thresholds (`docs/usage.md`), of the CNV QC output (`docs/output.md`) and descriptions of the CNV QC params in `nextflow_schema.json`.

## v0.1dev - [date]

Initial release of core-unit-bioinformatics/exomeqc, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

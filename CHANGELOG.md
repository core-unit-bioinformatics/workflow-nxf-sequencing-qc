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
- VAF histogram in the MultiQC report _(Implemented with support from Claude (Anthropic).)_: new local module `VAF_HISTOGRAM` running `bin/vaf_histogram.py` (Python 3.12, standard library) on the VCFs given in the new samplesheet column `vcf` (any variant caller). New MultiQC section "Variant allele fraction" (order 860, after the CNV sections): one histogram per sample, x-axis VAF, y-axis number of variants, drawn as a step outline in a MultiQC line graph (MultiQC bar graphs are horizontal only). Outputs in `vaf-histogram/`: `vaf_histogram.tsv`, `vaf_histogram_summary.tsv`, `vaf_histogram_mqc.json`.
  - Samplesheet: new optional column `vcf` (`.vcf`/`.vcf.gz`). All rows of a sample must give the same `vcf` or none, otherwise the pipeline stops with an `ERROR`.
  - New parameters: `--enable_vaf_histogram` (default `true`; stops with an `ERROR` if no sample has a `vcf`), `--vaf_histogram_allele_fraction` (`auto`, `FORMAT_AF` or `FORMAT_AD`, default `auto`), `--vaf_histogram_bin_width` (default `0.05`), `--vaf_histogram_sample_regex` (default: all VCF columns), `--vaf_histogram_include_filtered` (default `false`, only PASS variants are counted).
  - VAF as in cio-abcd/variantinterpretation: FORMAT/AF, or FORMAT/AD[alt] / FORMAT/DP; `auto` uses FORMAT/AF if the VCF header declares it, else FORMAT/AD + DP. The source per VCF is logged and written to the summary. A VCF without the needed FORMAT fields, or a counted variant without a value, stops the pipeline with an `ERROR` (no fallback). Each ALT allele counts as one variant; with FORMAT/GT only the ALT alleles in the sample's genotype count (hom-ref/no-call skipped and counted); records with FILTER `.` (no filters applied) count like PASS and are logged; bin edges use decimal arithmetic.
  - Report change: the MultiQC report gains this section; no existing section or output changes.
  - Synthetic test VCFs in DRAGEN, Mutect2 and FreeBayes format (`tests/fixtures/vaf_histogram/`, with their generator script), pytest unit tests (`tests/python/test_vaf_histogram.py`), nf-tests for the module and a pipeline test (`tests/vaf_histogram.nf.test`).

### `Fixed`

- DRAGEN CNV QC ignored all thresholds _(Implemented with support from Claude (Anthropic).)_: `modules/local/dragen_cnv_qc/main.nf` built `task.ext.args` from the threshold params in `conf/modules.config` but never passed it to `bin/analyze_cnv_metrics.R`. Only the hard-coded check (`EstimatedTumorPurity=NA` → `FAIL`) was active, so e.g. a sample with a purity below `--purity_min` still got `PASS`. The thresholds are now passed to the script.
  - **Results change:** with thresholds set, samples that violate them now get `QC_STATUS=FAIL` with the reason in `QC_FAIL_REASONS` (e.g. `EstimatedTumorPurity:0.25<0.3`) in `cnv-metrics/cnv-metrics.csv` and the MultiQC section "CNV QC Metrics". Samples within the thresholds, and runs without thresholds, are unchanged. The production parameters (`pmad_max=0.2`, `outlier_baf_fraction_max=0.1`, `purity_min=0.3`) take effect for the first time.
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

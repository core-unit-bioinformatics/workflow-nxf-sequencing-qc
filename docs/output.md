# core-unit-bioinformatics/exomeqc: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Fingerprint check](#fingerprint-check) - Sample identity checks within the run and against the fingerprint database of previous runs
- [FastQ Screen](#fastq-screen) - Read-level contamination screening against multiple reference genomes
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Fingerprint check

<details markdown="1">
<summary>Output files</summary>

- `fingerprintcheck/`
  - `crosscheck.crosscheck_metrics.txt`: Picard CrosscheckFingerprints of all BAMs of this run against each other; shown as table and heatmap in MultiQC.
  - `fingerprints/<sample>_<runid>.fingerprint.vcf.gz` (+ `.tbi`): fingerprint of each sample (Picard ExtractFingerprint), VCF sample name `<sample>@<runid>`, FORMAT fields `AD`, `DP`, `PL` at the haplotype map SNPs.
  - `database/`
    - `exomeqc_fingerprint_database.vcf.gz` (+ `.tbi`): updated fingerprint database, i.e. the given `--fingerprint_database` without the fingerprints replaced by this run, plus all fingerprints of this run. Use it as `--fingerprint_database` of the next run.
    - `fingerprint_database.crosscheck_metrics.txt`: Picard CrosscheckFingerprints of this run's fingerprints against the updated database (all pairs). Not read by MultiQC.
    - `fingerprint_database_matches.tsv`: reported comparisons with fingerprints of other runs (see below).
    - `fingerprint_database_summary.tsv`: counts (`metric`, `value`): `lod_threshold`, `run_samples`, `previous_database_samples`, `replaced_fingerprints`, `database_samples`, `previous_runs_database_samples`, `comparisons_read`, `comparisons_same_run_dropped`, `comparisons_other_runs`, `unexpected_match`, `unexpected_mismatch`, `expected_match`, `inconclusive_same_sample`, `inconclusive_other_sample`, `expected_mismatch`, `reported_comparisons`.
    - `fingerprint_database_matches_mqc.tsv` or, when nothing is reported, `fingerprint_database_matches_mqc.html`: MultiQC section "Crosscheck Fingerprints: Database matches".

</details>

Columns of `fingerprint_database_matches.tsv`:

| Column                   | Description                                                                                         |
| ------------------------ | --------------------------------------------------------------------------------------------------- |
| `sample`, `runid`        | Sample of this run                                                                                  |
| `database_sample`        | Sample ID of the database fingerprint                                                               |
| `database_runid`         | Run of the database fingerprint                                                                     |
| `status`                 | `UNEXPECTED_MATCH`, `UNEXPECTED_MISMATCH`, `INCONCLUSIVE` (same sample ID only) or `EXPECTED_MATCH` |
| `lod_score`              | LOD score: log10 likelihood ratio of same vs different individual                                   |
| `lod_score_tumor_normal` | Tumor-aware LOD score with the sample of this run as tumor                                          |
| `lod_score_normal_tumor` | Tumor-aware LOD score with the sample of this run as normal                                         |

Rows are sorted by status (unexpected results first), sample and descending LOD. See [usage](usage.md#fingerprint-check-and-fingerprint-database) for the classification.

### FastQ Screen

<details markdown="1">
<summary>Output files</summary>

- `fastqscreen/`
  - `*_screen.txt`: tabular FastQ Screen result for each input FASTQ
  - `*_screen.html`: per-sample FastQ Screen report
  - `*_screen.png`: per-sample FastQ Screen summary plot

</details>

FastQ Screen is executed when `--enable_fastqscreen` is set. The `*_screen.txt` files are automatically collected by MultiQC.

### DRAGEN CNV QC

<details markdown="1">
<summary>Output files</summary>

- `cnv-metrics/`
  - `cnv-metrics.csv`: one row per sample with `QC_STATUS` (`PASS`/`FAIL`), `QC_FAIL_REASONS` (`none`, or the failed checks separated by `;`, e.g. `EstimatedTumorPurity:0.25<0.3`) and the DRAGEN CNV metrics the checks use. Shown in the MultiQC section "CNV QC Metrics".
  - `purcov_maxmodel_lineplot.json`: maximum log-likelihood per purity of the DRAGEN purity/coverage models, for the MultiQC section "CNV: Purity/Coverage Model".
  - `karyoplot-<sample>-tn_mqc.png`, `karyoplot-<sample>-tumor.baf_mqc.png`: genome-wide tumor/normal ratio and tumor B-allele frequency.

</details>

The thresholds are described in [usage: DRAGEN CNV QC](usage.md#dragen-cnv-qc).

### VAF histogram

<details markdown="1">
<summary>Output files</summary>

- `vaf-histogram/`
  - `vaf_histogram.tsv`: variants per histogram and VAF bin. Columns: `sample` (histogram name: the samplesheet sample, or `<sample>_<column>` if several VCF columns are counted), `vcf_sample` (VCF column), `vcf` (file name), `bin_start`, `bin_end` (bin [start, end); the last bin includes 1.0), `variants`.
  - `vaf_histogram_summary.tsv`: per histogram, the VAF source and what was read and counted. Columns: `sample`, `vcf_sample`, `vcf`, `vaf_source` (`FORMAT_AF` or `FORMAT_AD`), `records_read`, `records_filtered` (skipped, FILTER not PASS or `.`), `alleles_counted`, `alleles_non_variant` (skipped, ALT `.`, `*`, `<NON_REF>` or `<*>`), `alleles_not_in_genotype` (skipped, ALT allele not in the sample's GT, e.g. hom-ref or no-call), `records_filter_missing` (records with FILTER `.`, no filters applied, counted like `PASS`).
  - `vaf_histogram_mqc.json`: the MultiQC section "Variant allele fraction".

</details>

The parameters are described in [usage: VAF histogram](usage.md#vaf-histogram).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

With `--multiqc_report_name <name>` (e.g. the run ID), these are named `<name>.html`, `<name>_data/` and `<name>_plots/`. Without it but with `--multiqc_title <title>`, MultiQC names them `<title>_multiqc_report.html` etc., with spaces in the title replaced by `-`.

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools (including FastQ Screen when enabled). The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

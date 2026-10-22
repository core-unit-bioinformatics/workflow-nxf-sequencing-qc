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
    - `exomeqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz` (+ `.tbi`): updated fingerprint database, named with the start time of the pipeline run: the kept fingerprints of the given `--fingerprint_database` (all except those replaced by this run, see the warning in [usage](usage.md#updated-fingerprint-database)) and all fingerprints of this run; the sample order follows the input file names of bcftools merge. Use it as `--fingerprint_database` of the next run.
    - `previous_database_check.tsv` (only with `--fingerprint_database`): result of the check of the given database against the haplotype map (`metric`, `value`): `database_file`, `database_compression` (`plain`, `gzip` or `bgzip`), `database_samples`, `database_runs`, `database_sites`, `database_missing_genotypes`, `haplotype_map_snps`, `haplotype_map_blocks`, `replaced_fingerprints`, `kept_fingerprints`.
    - `previous_database_samples.tsv` (only with `--fingerprint_database`): every fingerprint of the given database: `fingerprint_id`, `sample`, `runid`, `action` (`keep`, or `replace` when this run has a fingerprint with the same `<sample>@<runid>`).
    - `fingerprint_database.crosscheck_metrics.txt` (only when the given database has kept fingerprints): Picard CrosscheckFingerprints of this run's fingerprints against the kept database fingerprints, before this run was added (all pairs). Not read by MultiQC's Picard module.
    - `fingerprint_database_fingerprints.tsv`: fingerprint of each sample of this run (see below).
    - `fingerprint_database_matches.tsv`: reported comparisons with fingerprints of other runs (see below).
    - `fingerprint_database_summary.tsv`: counts (`metric`, `value`): `lod_threshold`, `haplotype_map_file`, `haplotype_map_snps`, `haplotype_map_blocks`, `run_samples`, `run_samples_suitable`, `run_samples_insufficient`, `run_samples_no_coverage`, `database_given` (`yes`/`no`), `database_file`, `database_compression`, `previous_database_samples`, `replaced_fingerprints`, `checked_database_samples`, `checked_database_runs`, `updated_database_file`, `database_samples` (in the updated database), `comparisons`, `unexpected_match`, `unexpected_mismatch`, `expected_match`, `inconclusive_same_sample`, `inconclusive_other_sample`, `expected_mismatch`, `reported_comparisons`.
    - `fingerprint_database_overview_mqc.html`: MultiQC section "Fingerprint check: overview".
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

Columns of `fingerprint_database_fingerprints.tsv` (one row per sample of this run):

| Column                       | Description                                                                                       |
| ---------------------------- | ------------------------------------------------------------------------------------------------- |
| `fingerprint_id`             | `<sample>@<runid>`                                                                                |
| `sample`, `runid`            | Sample and run                                                                                    |
| `haplotype_blocks`           | Haplotype blocks in the fingerprint VCF (= blocks of the haplotype map)                           |
| `blocks_with_reads`          | Blocks with at least one read (`DP` > 0)                                                          |
| `fraction_blocks_with_reads` | `blocks_with_reads` / `haplotype_blocks`                                                          |
| `mean_depth`                 | Mean `DP` over all blocks                                                                         |
| `self_lod`                   | LOD of the sample's BAM compared with itself in `crosscheck.crosscheck_metrics.txt`               |
| `status`                     | `SUITABLE` (self LOD ≥ T), `INSUFFICIENT` (self LOD < T) or `NO_COVERAGE` (no block with reads)   |

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

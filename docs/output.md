# core-unit-bioinformatics/workflow-nxf-sequencing-qc: Output

## Introduction

This page describes the output of the pipeline. All paths are relative to `--outdir`.
The output layout is the same for gene panels, exomes and genomes. Most of the plots are taken
from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Fingerprint check](#fingerprint-check) - Sample identity checks within the run and against the fingerprint database of previous runs, with fingerprint QC (any assay, BAM input)
- [FastQ Screen](#fastq-screen) - Read-level contamination / species screening against multiple reference genomes (any assay, single- or paired-end)
- [DRAGEN CNV QC](#dragen-cnv-qc) - QC of **DRAGEN** CNV caller output (DRAGEN only)
- [VAF histogram](#vaf-histogram) - Variant allele fraction histogram per sample (any caller)
- [Summary](#summary) - Compact per-case summary (CNV QC, purity, ploidy, TMB, MSI, HRD) on top of the report (DRAGEN output)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Fingerprint check

<details markdown="1">
<summary>Output files</summary>

- `samtools/`
  - `*.bai`: BAM indexes created for the fingerprint check.
- `fingerprintcheck/`
  - `crosscheck.crosscheck_metrics.txt`: Picard CrosscheckFingerprints of all BAMs of this run against each other (`CROSSCHECK_BY FILE`, with tumor-aware LOD scores); shown as table and heatmap in MultiQC. Main columns see below.
  - `fingerprints/<sample>_<runid>.fingerprint.vcf.gz` (+ `.tbi`): fingerprint of each sample after the [fingerprint QC](usage.md#fingerprint-qc) (Picard ExtractFingerprint, then bcftools annotate), VCF sample name `<sample>@<runid>`, one SNP per haplotype block. FORMAT fields `AD`, `DP` (Picard), `PL` (Picard; `0,0,0` at failed sites), `FT` (`PASS` or the failed filters `FpLowDepth`, `FpLowBaseQuality`, `FpNoLikelihoods`), `FPDP`, `FPAD`, `FPBQ` (depth, REF/ALT depth and mean base quality of the QC). The `##SAMPLE` header line holds the QC status, the values and the thresholds.
  - `qc/<sample>_<runid>.fingerprint_qc_summary.tsv`: fingerprint QC of the sample, one row with the columns of `fingerprint_database_qc.tsv` (see below).
  - `qc/<sample>_<runid>.fingerprint_qc_sites.tsv`: QC of each haplotype site: `#CHROM`, `POS`, `REF`, `ALT`, `FT`, `FPDP`, `FPAD`, `FPBQ`, `PL` (as written into the fingerprint).
  - `database/`
    - `exomeqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz` (+ `.tbi`): updated fingerprint database, named with the start time of the pipeline run: the kept fingerprints of the given `--fingerprint_database` (all except those replaced by this run, see the warning in [usage](usage.md#updated-fingerprint-database)) and all fingerprints of this run; the sample order follows the input file names of bcftools merge. Each fingerprint keeps its site QC fields and has exactly one `##SAMPLE` header line with its QC status, values and thresholds (`ExomeqcFingerprintQC=NOT_RECORDED` for fingerprints added before the fingerprint QC existed). Use it as `--fingerprint_database` of the next run.
    - `previous_database_check.tsv` (only with `--fingerprint_database`): result of the check of the given database against the haplotype map (`metric`, `value`): `database_file`, `database_compression` (`plain`, `gzip` or `bgzip`), `database_samples`, `database_runs`, `database_sites`, `database_missing_genotypes`, `haplotype_map_snps`, `haplotype_map_blocks`, `replaced_fingerprints`, `kept_fingerprints`.
    - `previous_database_samples.tsv` (only with `--fingerprint_database`): every fingerprint of the given database: `fingerprint_id`, `sample`, `runid`, `action` (`keep`, or `replace` when this run has a fingerprint with the same `<sample>@<runid>`).
    - `fingerprint_database.crosscheck_metrics.txt` (only when the given database has kept fingerprints): Picard CrosscheckFingerprints of this run's fingerprints against the kept database fingerprints, before this run was added (all pairs). Not read by MultiQC's Picard module.
    - `fingerprint_database_fingerprints.tsv`: fingerprint of each sample of this run (see below).
    - `fingerprint_database_qc.tsv`: fingerprint QC of each sample of this run (see below).
    - `fingerprint_database_qc_mqc.tsv`: MultiQC section "Fingerprint QC".
    - `fingerprint_database_matches.tsv`: reported comparisons with fingerprints of other runs (see below).
    - `fingerprint_database_summary.tsv`: counts (`metric`, `value`): `lod_threshold`, `haplotype_map_file`, `haplotype_map_snps`, `haplotype_map_blocks`, `run_samples`, `run_samples_suitable`, `run_samples_insufficient`, `run_samples_no_coverage`, `run_samples_qc_pass`, `run_samples_qc_fail`, `database_qc_pass`, `database_qc_fail`, `database_qc_not_recorded` (QC status of the fingerprints in the updated database), `database_given` (`yes`/`no`), `database_file`, `database_compression`, `previous_database_samples`, `replaced_fingerprints`, `checked_database_samples`, `checked_database_runs`, `updated_database_file`, `database_samples` (in the updated database), `comparisons`, `unexpected_match`, `unexpected_mismatch`, `expected_match`, `inconclusive_same_sample`, `inconclusive_other_sample`, `expected_mismatch`, `reported_comparisons`.
    - `fingerprint_database_overview_mqc.html`: MultiQC section "Fingerprint check: overview".
    - `fingerprint_database_matches_mqc.tsv` or, when nothing is reported, `fingerprint_database_matches_mqc.html`: MultiQC section "Crosscheck Fingerprints: Database matches".

</details>

Main columns of `crosscheck.crosscheck_metrics.txt`:

| Column                                              | Meaning                                                                                             |
| --------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| `LEFT_FILE` / `RIGHT_FILE`                          | The two BAM files compared.                                                                         |
| `RESULT`                                            | `EXPECTED_MATCH`, `EXPECTED_MISMATCH`, `UNEXPECTED_MATCH`, `UNEXPECTED_MISMATCH` or `INCONCLUSIVE`. |
| `LOD_SCORE`                                         | Log-odds that both files come from the same individual (positive = match).                          |
| `LOD_SCORE_TUMOR_NORMAL` / `LOD_SCORE_NORMAL_TUMOR` | Tumor-aware LOD scores (loss of heterozygosity taken into account).                                 |

`INCONCLUSIVE` usually means too few haplotype-map SNPs were covered. This is expected for small
gene panels with a genome-wide haplotype map; see [usage: assay-type notes](usage.md#assay-type-notes).

Columns of `fingerprint_database_matches.tsv`:

| Column                   | Description                                                                                         |
| ------------------------ | --------------------------------------------------------------------------------------------------- |
| `sample`, `runid`        | Sample of this run                                                                                  |
| `database_sample`        | Sample ID of the database fingerprint                                                               |
| `database_runid`         | Run of the database fingerprint                                                                     |
| `fingerprint_qc`         | Fingerprint QC of the sample of this run (`PASS`, `FAIL`)                                           |
| `database_fingerprint_qc` | Fingerprint QC of the database fingerprint (`PASS`, `FAIL`, `NOT_RECORDED`)                        |
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
| `fingerprint_qc`             | Fingerprint QC status (`PASS`, `FAIL`)                                                            |

Columns of `fingerprint_database_qc.tsv` (one row per sample of this run; the same values are in the `##SAMPLE` header line of the fingerprint):

| Column                                                                                  | Description                                                                                                                  |
| --------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `fingerprint_id`, `sample`, `runid`                                                     | `<sample>@<runid>`, sample and run                                                                                           |
| `status`                                                                                | `PASS` or `FAIL`                                                                                                             |
| `failed_criteria`                                                                       | `FractionPassingSites` and/or `Contamination` separated by `\|`, `.` if none                                                 |
| `qc_date`, `exomeqc_version`                                                            | Start date of the run (`YYYY-MM-DD`) and pipeline version                                                                    |
| `haplotype_sites`, `passing_sites`, `fraction_passing_sites`                            | Haplotype sites of the fingerprint, sites that pass, and their fraction                                                      |
| `sites_low_depth`, `sites_low_base_quality`, `sites_no_likelihoods`                     | Sites with the filters `FpLowDepth`, `FpLowBaseQuality`, `FpNoLikelihoods` (a site can have several)                         |
| `mean_site_depth`, `mean_site_base_quality`                                             | Mean `FPDP` over all sites; mean `FPBQ` over the sites with reads                                                            |
| `contamination`, `contamination_sites`                                                  | Estimated contamination and the number of passing sites it is based on; `.` and `0` without passing sites or without `--fingerprint_qc_contamination_check` |
| `min_base_quality`, `min_mapping_quality`, `min_site_depth`, `min_site_mean_base_quality` | Site thresholds used (`--fingerprint_qc_*`)                                                                                  |
| `min_fraction_passing_sites`, `max_contamination`                                       | Sample thresholds used (`--fingerprint_qc_*`); `max_contamination` is `.` when the contamination was not checked |

### FastQ Screen

<details markdown="1">
<summary>Output files</summary>

- `fastqscreen/`
  - `*_screen.txt`: tabular FastQ Screen result for each input FASTQ
  - `*_screen.html`: per-sample FastQ Screen report
  - `*_screen.png`: per-sample FastQ Screen summary plot

</details>

Written when `--enable_fastqscreen` is true (default). Paired-end samples are screened with `--paired`,
single-end samples read by read. The `*_screen.txt` files are automatically collected by MultiQC.

### DRAGEN CNV QC

<details markdown="1">
<summary>Output files</summary>

- `cnv-metrics/`
  - `cnv-metrics.csv`: one row per sample with `QC_STATUS` (`PASS`/`FAIL`), `QC_FAIL_REASONS` (`none`, or the failed checks separated by `;`, e.g. `EstimatedTumorPurity:0.25<0.3`) and the DRAGEN CNV metrics the checks use. Shown in the MultiQC section "DRAGEN CNV QC Metrics".
  - `purcov_maxmodel_lineplot.json`: maximum log-likelihood per purity of the DRAGEN purity/coverage models, for the MultiQC section "DRAGEN CNV: Purity/Coverage Model".
  - `karyoplot-<sample>-tn_mqc.png`, `karyoplot-<sample>-tumor.baf_mqc.png`: genome-wide tumor/normal ratio and tumor B-allele frequency.

</details>

Written only when `--enable_cnvqc` is true (default). This step is **specific to DRAGEN CNV caller
output** and is not run for other callers. The thresholds are described in
[usage: DRAGEN CNV QC](usage.md#dragen-cnv-qc).

### VAF histogram

<details markdown="1">
<summary>Output files</summary>

- `vaf-histogram/`
  - `vaf_histogram.tsv`: variants per histogram and VAF bin. Columns: `sample` (histogram name: the samplesheet sample, or `<sample>_<column>` if several VCF columns are counted), `vcf_sample` (VCF column), `vcf` (file name), `bin_start`, `bin_end` (bin [start, end); the last bin includes 1.0), `variants`.
  - `vaf_histogram_summary.tsv`: per histogram, the VAF source and what was read and counted. Columns: `sample`, `vcf_sample`, `vcf`, `vaf_source` (`FORMAT_AF` or `FORMAT_AD`), `records_read`, `records_filtered` (skipped, FILTER not PASS or `.`), `alleles_counted`, `alleles_non_variant` (skipped, ALT `.`, `*`, `<NON_REF>` or `<*>`), `alleles_not_in_genotype` (skipped, ALT allele not in the sample's GT, e.g. hom-ref or no-call), `records_filter_missing` (records with FILTER `.`, no filters applied, counted like `PASS`).
  - `vaf_histogram_mqc.json`: the MultiQC section "Variant allele fraction".

</details>

The parameters are described in [usage: VAF histogram](usage.md#vaf-histogram).

### Summary

<details markdown="1">
<summary>Output files</summary>

- `summary/`
  - `exomeqc_summary.html`: the **Summary** shown on top of the MultiQC report (standalone fragment).
  - `exomeqc_summary_multiqc_config.yml`: MultiQC config that puts the summary below the report title (`intro_text`), above General Statistics.
  - `exomeqc_summary.tsv`: compact table, one row per case.
  - `biomarker.tsv`: all DRAGEN TMB, MSI and HRD values per case.

</details>

`exomeqc_summary.tsv` columns:

| Column   | Description                                                                                    |
| -------- | ---------------------------------------------------------------------------------------------- |
| `sample` | DRAGEN case (output file prefix)                                                               |
| `cnv_qc` | `PASS`/`FAIL` from CNV QC (`n/a` if CNV QC is disabled or the case is missing)                 |
| `purity` | DRAGEN `EstimatedTumorPurity`, in %                                                                 |
| `ploidy` | DRAGEN `OverallPloidy`                                                                         |
| `tmb`    | DRAGEN `TMB` (mut/Mb)                                                                          |
| `msi`    | DRAGEN `PercentageUnstableSites`; suffix `(invalid)` if `ResultIsValid` is false               |
| `hrd`    | DRAGEN `HRD_sum` (LOH + TAI + LST)                                                             |

`biomarker.tsv` has the same columns as `biomarker.csv` from the diagnostic somatic-exome workflow (`05-process-tsv-biomarker.py`): all TMB metrics, all MSI fields except `Settings`, then `LOH_score`, `TAI_score`, `LST_score`, `HRD_sum`. Unlike that script, HRD rows are keyed by file prefix, so all three biomarkers of one case share a row. Missing values are `n/a`.

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

Results generated by MultiQC collate pipeline QC from supported tools (including FastQ Screen when enabled) and
everything MultiQC recognises in `--multiqc_inputs` and `--dragen_dirs`. For DRAGEN this includes mapping metrics,
WGS coverage (genomes) or target-region coverage (panels/exomes), and variant-calling metrics. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - `sequencingqc_software_mqc_versions.yml`: versions of all tools used.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

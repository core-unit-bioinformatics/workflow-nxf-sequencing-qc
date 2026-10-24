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
  - New parameters: `--enable_vaf_histogram` (default `false`, so existing samplesheets without `vcf` keep working; if set, stops with an `ERROR` if no sample has a `vcf`), `--vaf_histogram_allele_fraction` (`auto`, `FORMAT_AF` or `FORMAT_AD`, default `auto`), `--vaf_histogram_bin_width` (default `0.05`), `--vaf_histogram_sample_regex` (default: all VCF columns), `--vaf_histogram_include_filtered` (default `false`, only variants with FILTER PASS or `.` are counted).
  - VAF as in cio-abcd/variantinterpretation: FORMAT/AF, or FORMAT/AD[alt] / FORMAT/DP; `auto` uses FORMAT/AF if the VCF header declares it, else FORMAT/AD + DP. The source per VCF is logged and written to the summary. A VCF without the needed FORMAT fields, or a counted variant without a value, stops the pipeline with an `ERROR` (no fallback). Each ALT allele counts as one variant; with FORMAT/GT only the ALT alleles in the sample's genotype count (hom-ref/no-call skipped and counted); records with FILTER `.` (no filters applied) count like PASS and are logged; bin edges use decimal arithmetic.
  - Report change: the MultiQC report gains this section; no existing section or output changes.
  - Synthetic test VCFs in DRAGEN, Mutect2 and FreeBayes format (`tests/fixtures/vaf_histogram/`, with their generator script), pytest unit tests (`tests/python/test_vaf_histogram.py`), nf-tests for the module and a pipeline test (`tests/vaf_histogram.nf.test`).
- Fingerprint database across sequencing runs _(Implemented with support from Claude (Anthropic).)_:
  - Picard ExtractFingerprint writes one fingerprint VCF per sample, named `<sample>@<runid>` (`fingerprintcheck/fingerprints/`).
  - Fingerprint database of previous runs (new parameter `--fingerprint_database`, optional; plain VCF, gzip or bgzip, with or without index). It is first checked against `--fingerprint_haplotype_map` by the new local module `CHECK_FINGERPRINT_DATABASE` (`bin/check_fingerprint_database.py`, PyVCF): fingerprint IDs, `##contig` names and lengths against the map's sequence dictionary, `FORMAT/PL`, every record an SNP of the map with the map's alleles, every haplotype block present exactly once. Any mismatch stops the pipeline with an `ERROR` (before, a database of another reference failed late with an unclear Picard error, and one of another haplotype map on the same reference was not detected).
  - `BCFTOOLS_VIEW_FINGERPRINT_DATABASE` keeps the database fingerprints that are not replaced by this run and writes them as bgzip with a new tabix index; the index of the given file is not needed.
  - Fingerprints of this run already in the given database (same sample **and** run ID, `<sample>@<runid>`) are dropped and replaced by the new ones; this is logged as `WARNING`, shown in the MultiQC overview, listed in `previous_database_samples.tsv` and counted. The docs state this as a warning.
  - Picard CrosscheckFingerprints (`PICARD_CROSSCHECKFINGERPRINTS_DATABASE`) compares every sample of this run with every kept database fingerprint, before this run is added to the database, so no sample is compared with its own new fingerprint. Then `BCFTOOLS_MERGE_FINGERPRINT_DATABASE` appends this run's fingerprints: updated database `fingerprintcheck/database/exomeqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz` (start time of the pipeline run), to be passed to the next run.
  - `bin/fingerprint_database_report.py` (PyVCF) classifies the comparisons with other runs by LOD score (new parameter `--fingerprint_database_lod_threshold`, default 5) and writes `fingerprint_database_matches.tsv`, `fingerprint_database_fingerprints.tsv`, `fingerprint_database_summary.tsv` and two MultiQC sections: "Fingerprint check: overview" (haplotype map SNPs and blocks, database fingerprints checked, replaced fingerprints, comparison counts, updated database, and per sample of this run the haplotype blocks with reads, mean depth, self LOD and whether the fingerprint is `SUITABLE`, `INSUFFICIENT` (self LOD < T) or has `NO_COVERAGE`), placed before the Picard heatmap, and "Crosscheck Fingerprints: Database matches" (unexpected matches, mismatches and inconclusive comparisons of the same sample ID, and expected matches, with LOD scores), placed right after the Picard sections.
  - Optional alignment of samples without BAM with bwa-mem2 (new parameters `--align_reads_for_fingerprint`, default `false`, and `--bwamem2_index`).
  - New parameters `--fasta_fai` and `--fasta_dict`; both are created from `--fasta` if not given.
  - Synthetic test data for the fingerprint check (`tests/fixtures/fingerprint/`): made-up sample IDs (`SYN001_3-N` ...) and run IDs, a random reference and simulated reads; the generator script `make_fingerprint_fixtures.py` also makes the module test inputs with Picard and bcftools and removes paths and dates from their headers. pytest unit tests for both scripts and nf-tests for the new modules and the subworkflow.

### `Fixed`

<<<<<<< HEAD
- DRAGEN CNV QC ignored all thresholds: `modules/local/dragen_cnv_qc/main.nf` built `task.ext.args` from the threshold params in `conf/modules.config` but never passed it to `bin/analyze_cnv_metrics.R`. Only the hard-coded check (`EstimatedTumorPurity=NA` → `FAIL`) was active, so e.g. a sample with a purity below `--purity_min` still got `PASS`. The thresholds are now passed to the script.
  - **Results change:** with thresholds set, samples that violate them now get `QC_STATUS=FAIL` with the reason in `QC_FAIL_REASONS` (e.g. `EstimatedTumorPurity:0.25<0.3`) in `cnv-metrics/cnv-metrics.csv` and the MultiQC section "CNV QC Metrics". Samples within the thresholds, and runs without thresholds, are unchanged.
=======
- DRAGEN CNV QC ignored all thresholds _(Implemented with support from Claude (Anthropic).)_: `modules/local/dragen_cnv_qc/main.nf` built `task.ext.args` from the threshold params in `conf/modules.config` but never passed it to `bin/analyze_cnv_metrics.R`. Only the hard-coded check (`EstimatedTumorPurity=NA` → `FAIL`) was active, so e.g. a sample with a purity below `--purity_min` still got `PASS`. The thresholds are now passed to the script.
  - **Results change:** with thresholds set, samples that violate them now get `QC_STATUS=FAIL` with the reason in `QC_FAIL_REASONS` (e.g. `EstimatedTumorPurity:0.25<0.3`) in `cnv-metrics/cnv-metrics.csv` and the MultiQC section "CNV QC Metrics". Samples within the thresholds, and runs without thresholds, are unchanged. The production parameters (`pmad_max=0.2`, `outlier_baf_fraction_max=0.1`, `purity_min=0.3`) take effect for the first time.
>>>>>>> a13667f (Add parameter multiqc_report_name for the MultiQC report file name)
  - `conf/modules.config`: a threshold is passed if it is set (`!= null`). Before, a threshold of `0` was dropped silently.
  - `bin/analyze_cnv_metrics.R` stops with an `ERROR` on unknown options (e.g. `--purity-min` instead of `--purity_min`), which would otherwise leave that check silently disabled, and prints all errors as `ERROR: <message>` with exit status 1.
  - The stub of `DRAGEN_CNV_QC` now creates the same files as a real run with the default `ext.args` (`cnv-metrics.csv`, `purcov_maxmodel_lineplot.json`, `karyoplot-*_mqc.png`, with `ext.prefix`).
- `DRAGEN_CNV_QC` failed with a missing-output error if the DRAGEN output had no `*.cnv.purity.coverage.models.tsv` (no `*.json`) or no `*.tn.tsv.gz`/`*.tumor.baf.bedgraph.gz` (no `*.png`). These outputs are now optional, as documented in `docs/usage.md`. No change in results.
- `bin/analyze_cnv_metrics.R` was committed without the executable bit (mode 100644), so `DRAGEN_CNV_QC` failed with `Permission denied` (exit 126) when run from a fresh clone. It is now executable (100755). No change in results.

### `Changed`

- Samplesheet: new column `runid`, required with `--enable_fingerprintcheck`. With the fingerprint check, `sample` and `runid` may only contain `A-Z`, `a-z`, `0-9`, `.`, `_`, `-`.
- `--fasta` (the reference the BAMs are aligned to) is now required with `--enable_fingerprintcheck`, because Picard ExtractFingerprint needs it. It must be uncompressed.
- Rows without `bam` now stop the pipeline with an `ERROR` when the fingerprint check is on and `--align_reads_for_fingerprint` is not set. Before, such rows were silently left out of the fingerprint check (the check for the `bam` column never failed).
- Rows of the same `sample` that differ in `bam` or `runid` now stop the pipeline with an `ERROR`. Before, the values of the first row were used silently.
- Results: the crosscheck of this run's BAMs (`fingerprintcheck/crosscheck.crosscheck_metrics.txt`, Picard 3.4.0, same arguments) has the same LOD scores and results as v1.1.1 for samples with a BAM (compared on the synthetic test data). As before, the order of the samples in the heatmap is set by Picard and can differ between runs. Samples aligned with `--align_reads_for_fingerprint` are added to it. The MultiQC report gains the new database section.
- MultiQC: `fingerprint_database.crosscheck_metrics.txt` is excluded from the Picard module (`sp.picard/crosscheckfingerprints.exclude_fn`), so database samples never enter the heatmap.

### `Dependencies`

- New nf-core modules: `picard/extractfingerprint` (Picard 3.5.0), `bcftools/view` and `bcftools/merge` (bcftools 1.23.1), `bwamem2/index` and `bwamem2/mem` (bwa-mem2 2.3, whose binary reports version 2.2.1; samtools 1.22.1), `samtools/faidx` (samtools 1.24), `picard/createsequencedictionary` (Picard 3.5.0).
- New local modules `check_fingerprint_database` and `fingerprint_database_report` (Python 3.12, PyVCF3 1.0.4, the maintained fork of PyVCF; biocontainer `pyvcf3:1.0.4--py312h0fa9677_0`).

## v0.1dev - [date]

Initial release of core-unit-bioinformatics/exomeqc, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

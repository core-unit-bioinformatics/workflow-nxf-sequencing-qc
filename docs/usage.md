# core-unit-bioinformatics/exomeqc: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

Columns for the fingerprint check (`--enable_fingerprintcheck`, on by default):

| Column  | Description                                                                                                                                                                                                                                                                                      |
| ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `bam`   | Full path to the BAM/CRAM file of the sample, aligned to `--fasta`. Required for every row unless `--align_reads_for_fingerprint` is set; then rows without `bam` are aligned from `fastq_1`/`fastq_2` with bwa-mem2.                                                                            |
| `runid` | RUNID of the sequencing run the sample comes from, e.g. `20260101_LH00001_0001_AFIXTURE01`. Required for every row. Together with `sample` it identifies a fingerprint as `<sample>@<runid>`, because the same sample can be sequenced in several runs. Only `A-Z`, `a-z`, `0-9`, `.`, `_`, `-`. |

With the fingerprint check, `sample` may only contain `A-Z`, `a-z`, `0-9`, `.`, `_` and `-` (it becomes part of VCF sample and file names). Each sample may appear in only one `bam`/`runid` combination per samplesheet: rows of the same `sample` are merged and must agree on `bam` and `runid`.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,bam,runid
SYN001_1-T,SYN001_1-T_R1.fq.gz,SYN001_1-T_R2.fq.gz,SYN001_1-T_tumor.bam,20260101_LH00001_0001_AFIXTURE01
SYN001_3-N,SYN001_3-N_R1.fq.gz,SYN001_3-N_R2.fq.gz,SYN001_3-N.bam,20260101_LH00001_0001_AFIXTURE01
```

Optional column for the VAF histogram:

- `vcf`: Full path to a small-variant VCF (`.vcf` or `.vcf.gz`, any caller) of the sample, used for the [VAF histogram](#vaf-histogram). All rows of a sample must give the same `vcf` or none.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run core-unit-bioinformatics/exomeqc --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

To enable FastQ Screen and include its metrics in MultiQC, add:

```bash
--enable_fastqscreen --fastqscreen_database /path/to/FastQ_Screen_Genomes
```

The `--fastqscreen_database` directory must contain `fastq_screen.conf` and the referenced genome index folders. It can also be given as a `.tar.gz` of that directory (local path or URL); the pipeline unpacks it into a directory named like the archive without `.tar.gz`. `DATABASE` paths in `fastq_screen.conf` must be absolute or relative to the parent of the database directory (`<database_dir_name>/<genome>/<index_prefix>`).

### Fingerprint check and fingerprint database

The fingerprint check (`--enable_fingerprintcheck`, on by default) needs the Picard haplotype map and the reference FASTA the BAMs are aligned to:

```bash
--enable_fingerprintcheck \
--fingerprint_haplotype_map /path/to/haplotype_map.txt \
--fasta /path/to/genome.fasta \
--fingerprint_database /path/to/previous/exomeqc_fingerprint_database.vcf.gz
```

It runs two checks:

1. **This run.** Picard CrosscheckFingerprints compares the BAMs of all samples of this run with each other (heatmap in the MultiQC report, unchanged from earlier versions). `--fingerprint_allow_duplicate_reads` (default `true`) controls whether duplicate reads are used here.
2. **Fingerprint database.** Picard ExtractFingerprint writes one fingerprint VCF per sample with the sample name `<sample>@<runid>`. Each fingerprint passes the [fingerprint QC](#fingerprint-qc) before it is used; failed haplotype sites are masked and the QC result is stored in the fingerprint. If a database of previous runs is given (`--fingerprint_database`, optional), the steps are:
   1. **Check.** The database is checked against `--fingerprint_haplotype_map` (`CHECK_FINGERPRINT_DATABASE`, see below). Fingerprints with the same `<sample>@<runid>` as a sample of this run are marked for replacement.
   2. **Prepare.** `BCFTOOLS_VIEW_FINGERPRINT_DATABASE` keeps all other fingerprints and writes them as bgzipped VCF with a new tabix index. The database may be given as plain VCF, gzip or bgzip, with or without index; the index of the given file is not used. The given file itself is never modified.
   3. **Crosscheck.** Every sample of this run is compared with every kept database fingerprint (`PICARD_CROSSCHECKFINGERPRINTS_DATABASE`: CrosscheckFingerprints, `CROSSCHECK_BY SAMPLE`, `CROSSCHECK_MODE CHECK_ALL_OTHERS`, this run as `INPUT`, the database as `SECOND_INPUT`). This happens _before_ this run is added, so no sample is compared with its own new fingerprint; comparisons within this run are covered by check 1.
   4. **Update.** This run's fingerprints are added to the kept fingerprints (`BCFTOOLS_MERGE_FINGERPRINT_DATABASE`), see [Updated database](#updated-fingerprint-database).

   Without `--fingerprint_database`, steps 1–3 are skipped and the updated database holds only this run's fingerprints. Comparisons with other runs are classified with the threshold T = `--fingerprint_database_lod_threshold` (default `5`):

   | LOD          | Same sample ID (expected match) | Different sample ID |
   | ------------ | ------------------------------- | ------------------- |
   | LOD ≥ T      | `EXPECTED_MATCH`                | `UNEXPECTED_MATCH`  |
   | −T < LOD < T | `INCONCLUSIVE`                  | `INCONCLUSIVE`      |
   | LOD ≤ −T     | `UNEXPECTED_MISMATCH`           | `EXPECTED_MISMATCH` |

   The MultiQC table "Crosscheck Fingerprints: Database matches" shows `UNEXPECTED_MATCH`, `UNEXPECTED_MISMATCH`, `EXPECTED_MATCH` and the `INCONCLUSIVE` comparisons of the same sample ID, each with its LOD score. `EXPECTED_MISMATCH` and `INCONCLUSIVE` comparisons between different sample IDs are only counted (`fingerprint_database_summary.tsv`). "Expected" refers only to identical sample IDs; whether two different sample IDs belong to the same patient is not judged by the pipeline.

   Picard's own `RESULT` column is not used for this table: it treats every database sample as a different sample (the `@<runid>` part always differs) and calls any LOD above −T between different samples a match.

The MultiQC section "Fingerprint check: overview" (`bin/fingerprint_database_report.py`) comes first and summarises the check: the haplotype map (number of SNPs and haplotype blocks; one SNP per block is compared), the database (file, number of fingerprints checked, number of runs), the replaced fingerprints, the counts of the comparison classes and the updated database. It also shows for every sample of this run whether its fingerprint is usable:

| Fingerprint    | Meaning                                                                                                                                                                                                       |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SUITABLE`     | Self LOD ≥ T: the LOD of the sample's BAM compared with itself (diagonal of the crosscheck heatmap) reaches the threshold, so a match with another fingerprint of the same person can be called.             |
| `INSUFFICIENT` | Self LOD < T: too few informative reads at the haplotype map SNPs. Comparisons of this sample stay `INCONCLUSIVE`; a missing match is no evidence of a sample swap.                                           |
| `NO_COVERAGE`  | No read at any haplotype block of the map, e.g. when the target regions do not cover the map's SNPs or the BAM is empty or of another reference.                                                             |

The number of haplotype blocks with reads and the mean depth per sample are listed as well. Fingerprints that are not `SUITABLE` are logged as `WARNING`; they are still added to the database.

#### Fingerprint QC

Each fingerprint is checked before it is compared with the database and added to it (`FINGERPRINT_QC`, `bin/fingerprint_qc.py`). samtools mpileup lists every read at the SNPs of the haplotype map with its base quality, mapping quality and read name (no BAQ, no quality filter, duplicates as set by `--fingerprint_allow_duplicate_reads`); the script applies all thresholds, so the database records exactly what was used. The fingerprint has one SNP per haplotype block, its _haplotype sites_.

| Parameter                                     | Default | Applies to                                                                                                                                                                                                            |
| --------------------------------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--fingerprint_qc_min_base_quality`           | `20`    | Each base: bases below are not counted in the site depth, the allele depths and the contamination estimate.                                                                                                          |
| `--fingerprint_qc_min_mapping_quality`        | `20`    | Each read: reads below are not counted at all.                                                                                                                                                                        |
| `--fingerprint_qc_min_site_depth`             | `10`    | Each site: `FPDP`, the number of bases passing both quality thresholds, must reach it. Each read name counts once per site (both reads of a pair are one DNA fragment), the read with the higher base quality is used. |
| `--fingerprint_qc_min_site_mean_base_quality` | `25`    | Each site: `FPBQ`, the mean base quality of all its reads with sufficient mapping quality, must reach it.                                                                                                            |
| `--fingerprint_qc_min_fraction_passing_sites` | `0.5`   | Each sample: fraction of haplotype sites that pass.                                                                                                                                                                   |
| `--fingerprint_qc_contamination_check`       | `false` | Each sample: estimate the contamination and use it as criterion. Off by default: not suited for tumour samples (see below).                                                                                        |
| `--fingerprint_qc_max_contamination`          | `0.03`  | Each sample, only with `--fingerprint_qc_contamination_check`: estimated contamination (fraction of reads from another person) must not exceed it.                                                                |

A site that fails gets the filters `FpLowDepth`, `FpLowBaseQuality` or `FpNoLikelihoods` in `FORMAT/FT`, and its genotype likelihoods are set to `PL=0,0,0`. Flat likelihoods carry no information, so Picard CrosscheckFingerprints gives the same LOD scores as if the site were missing (a missing `PL` makes Picard fail). Failed sites therefore never contribute to the comparisons with the database or to later runs.

**Contamination check (optional, off by default).** With `--fingerprint_qc_contamination_check`, the contamination is estimated by maximum likelihood from the high-quality bases at the passing sites, with the model of VerifyBamID's `FREEMIX` (Jun et al. 2012): a fraction of the reads comes from an unrelated person of the same population, the genotypes of the sample and of the contaminant follow Hardy-Weinberg proportions of the allele frequency (`MAF`) of the haplotype map, and base errors follow the base qualities. Without any passing site there is no estimate and the sample fails.

> [!WARNING]
> **Not suited for tumour samples.** The model assumes a diploid, copy-neutral sample. In tumours, copy-number-driven allelic imbalance looks like contamination: on TCGA glioblastoma data, VerifyBamID highly overestimated the contamination of the majority of the tumour samples (Bergmann et al. 2016, Conpair). Methods that work for tumours need a matched normal (Conpair, ContEst) or a copy number segmentation (GATK CalculateContamination). The check is therefore off by default; turn it on only for germline (normal) samples. Without it, contamination is no criterion and is recorded as not checked (`Contamination=.`, `MaxContamination=.`).

A sample is `PASS` when the fraction of passing sites (and, with the contamination check, the contamination) meet the thresholds, otherwise `FAIL` with the failed criteria (`FractionPassingSites`, `Contamination`). Failed samples are logged as `WARNING` and **still added to the database**, marked `FAIL`, so they remain visible when the database is reviewed. Set `--fingerprint_qc_min_fraction_passing_sites` for your haplotype map: with a map whose SNPs lie outside the targets of the exome kit, few sites reach the depth.

The QC result is stored with the fingerprint and travels with it into every later database:

- per site in the fingerprint VCF: `FT`, `FPDP`, `FPAD` (REF and ALT bases counted in `FPDP`) and `FPBQ`;
- per fingerprint as VCF header line `##SAMPLE=<ID=<sample>@<runid>,ExomeqcFingerprintQC=PASS|FAIL,FailedCriteria=...,QCDate=YYYY-MM-DD,ExomeqcVersion=...,HaplotypeSites=...,PassingSites=...,FractionPassingSites=...,SitesLowDepth=...,SitesLowBaseQuality=...,SitesNoLikelihoods=...,MeanSiteDepth=...,MeanSiteBaseQuality=...,Contamination=...,ContaminationSites=...,MinBaseQuality=...,MinMappingQuality=...,MinSiteDepth=...,MinSiteMeanBaseQuality=...,MinFractionPassingSites=...,MaxContamination=...>`, i.e. the status, the observed values and the thresholds used, so fingerprints of runs with different settings can be compared.

The updated database has exactly one such line per fingerprint (`FINGERPRINT_DATABASE_HEADER` and `BCFTOOLS_REHEADER_FINGERPRINT_DATABASE`; bcftools merge alone would keep outdated lines of replaced fingerprints). Fingerprints added before the fingerprint QC existed get `ExomeqcFingerprintQC=NOT_RECORDED`, logged as `WARNING`. The MultiQC section "Fingerprint QC" shows status, failed criteria, values and thresholds of each sample of this run; the database matches show the QC status of both fingerprints of a comparison (`QC`, `Database QC`).

#### Updated fingerprint database

The updated database is written to `fingerprintcheck/database/exomeqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz` (with `.tbi`). The timestamp is the start time of the pipeline run (system time), so the databases of successive runs do not overwrite each other. Pass it as `--fingerprint_database` to the next run so the database grows with every run. Keep the haplotype map and reference the same for all runs of one database.

> [!WARNING]
> **Fingerprints with the same sample AND run ID are dropped and replaced.** If the given database already contains a fingerprint with the same `<sample>@<runid>` as a sample of this run (e.g. when a run is analysed again), that old fingerprint is **dropped** from the updated database and **replaced** by the new fingerprint of this run. The old fingerprint is not kept anywhere in the updated database and is not compared with this run. This is logged as a `WARNING`, listed in the overview section of the MultiQC report, marked `replace` in `fingerprintcheck/database/previous_database_samples.tsv` and counted as `replaced_fingerprints`. Fingerprints of the same sample from _other_ runs are kept. Keep the given database file if the old fingerprints may be needed again.

#### Check of the given database

`CHECK_FINGERPRINT_DATABASE` (`bin/check_fingerprint_database.py`, PyVCF) reads the whole database and stops the pipeline with an `ERROR` that names the file, the record and the expected and observed values unless:

- all sample names are unique and of the form `<sample>@<runid>` (i.e. the database was created by exomeqc);
- every `##contig` of the VCF header is in the sequence dictionary (`@SQ` lines) of the haplotype map, with the same length;
- `FORMAT/PL` (genotype likelihoods, read by Picard) is defined and every genotype has 3 PL values or is missing;
- every record is an SNP of the haplotype map, at the same position and with the same two alleles (`REF`/`ALT` = `MAJOR_ALLELE`/`MINOR_ALLELE` in either order);
- every haplotype block of the map is present exactly once (Picard ExtractFingerprint writes one representative SNP per block, also when it has no reads).

A database built with another haplotype map, genome build or contig naming therefore stops the run before any comparison is made. Start a new database (omit `--fingerprint_database`) when the haplotype map or the reference changes.

Samples without `bam` can be aligned with bwa-mem2 for the fingerprint check (off by default):

```bash
--align_reads_for_fingerprint --bwamem2_index /path/to/bwamem2/index/dir
```

Without `--bwamem2_index`, the index is built from `--fasta`; for a human genome this needs about 64 GB of memory and several hours. The aligned BAMs are only used for the fingerprints and are not published. `--fasta_fai` and `--fasta_dict` are created from `--fasta` if not given. The FASTA must be uncompressed.

Tested inputs: Picard 3.4.0 (CrosscheckFingerprints) and 3.5.0 (ExtractFingerprint), bcftools 1.23.1 (view, merge) and 1.25 (annotate, reheader), samtools 1.25 (mpileup), HTSlib 1.25, GNU Awk 5.3.1, bwa-mem2 2.3 (conda package; the binary reports 2.2.1), PyVCF3 1.0.4 on the synthetic test data in `tests/fixtures/fingerprint/` (see `make_fingerprint_fixtures.py`), with the database given as bgzip with index and as plain VCF without index (nf-test), and as plain VCF, gzip and bgzip (pytest). The real haplotype map (e.g. Broad's `Homo_sapiens_assembly38.haplotype_database.txt`) was not tested here. The contamination estimate of the fingerprint QC was tested on simulated pileups (pytest: 0 to 20 % contamination recovered within 0.01 at 400 sites) and on the fixture BAM of `SYN001_3-N` mixed with 5 % and 11 % reads of another synthetic person (true fractions 0.048 and 0.099; estimates 0.047 and 0.096 at 100 sites); it was not validated on real exomes.

To name the MultiQC report after the run, add e.g.:

```bash
--multiqc_report_name 20261023_LH00001_0001_AHXXXXXXX
```

The report is then written to `<OUTDIR>/multiqc/20261023_LH00001_0001_AHXXXXXXX.html`. The name may contain `A-Z`, `a-z`, `0-9`, `.`, `_` and `-`, must start with a letter or digit and must not end in `.html`; other names stop the pipeline at parameter validation. `--multiqc_title` still sets the title printed in the report.

### Report summary

With `--dragen_dirs` set, the MultiQC report starts with a compact **Summary** section (`--enable_summary`, default `true`): one row per DRAGEN case with CNV QC status, purity, ploidy, TMB, MSI and HRD, plus run-level chips for CNV QC and fingerprint results. All further values are in collapsed `<details>` blocks below the table.

| Parameter          | Default | Description                                                                    |
| ------------------ | ------- | ------------------------------------------------------------------------------ |
| `enable_summary`   | `true`  | Add the summary section. Skipped with a warning if `--dragen_dirs` is not set. |
| `summary_tmb_high` | unset   | Highlight TMB (mut/Mb) at or above this value.                                 |
| `summary_msi_high` | unset   | Highlight MSI (% unstable sites) at or above this value.                       |
| `summary_hrd_high` | unset   | Highlight HRD sum (LOH + TAI + LST) at or above this value.                    |

The pipeline ships no clinical cutoffs: without a `summary_*_high` value the numbers are shown without highlighting. Set your lab's validated cutoffs in a params file.

Input requirements (searched recursively in all `--dragen_dirs`, case = file prefix):

- `<case>.tmb.metrics.csv` (DRAGEN TMB, must contain a `TMB` metric)
- `<case>.microsat_output.json` (DRAGEN MSI, must contain `PercentageUnstableSites`; `ResultIsValid: false` is shown as `(invalid)` and never highlighted)
- `<case>.hrdscore.csv` (DRAGEN HRD, exactly one data row with `HRD_sum`)

A case found twice across `--dragen_dirs`, or a malformed file, stops the pipeline with an error naming the file. A case with only some of the files is shown with `n/a` and listed under "Incomplete cases". Tested with synthetic files modelled on DRAGEN output (`tests/fixtures/summary/`); not yet with a real DRAGEN run.

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run core-unit-bioinformatics/exomeqc -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### DRAGEN CNV QC

With `--enable_cnvqc` (default `true`), `bin/analyze_cnv_metrics.R` reads the DRAGEN CNV output of all samples under `--dragen_dirs` (`*.cnv_metrics.csv`, `*.cnv.vcf.gz`, and if present `*.cnv.purity.coverage.models.tsv`, `*.cnv.excluded_intervals.bed.gz`, `*.tn.tsv.gz`, `*.tumor.baf.bedgraph.gz`) and sets a `QC_STATUS` of `PASS` or `FAIL` per sample.

A sample always fails if DRAGEN found no purity model (`EstimatedTumorPurity=NA`) or `ModelSource` is missing. All other checks are optional: a check runs only if its threshold is given. A metric that is missing while its threshold is set fails the sample (`<metric>:MISSING`).

| Parameter                       | Metric (source)                                              | Sample fails if  |
| ------------------------------- | ------------------------------------------------------------ | ---------------- |
| `pmad_max`                      | `PMAD` (cnv_metrics.csv)                                     | value > max      |
| `coverage_mad_max`              | `Coverage MAD` (cnv_metrics.csv)                             | value > max      |
| `outlier_baf_fraction_max`      | `OutlierBafFraction` (cnv_metrics.csv)                       | value > max      |
| `diploid_coverage_min`          | `DiploidCoverage` (CNV VCF header)                           | value < min      |
| `purity_min`, `purity_max`      | `EstimatedTumorPurity` (CNV VCF header), a fraction in [0, 1] | outside min..max |
| `ploidy_min`, `ploidy_max`      | `OverallPloidy` (CNV VCF header)                             | outside min..max |
| `post_norm_sigma_max`           | `Post-Normalization Bin Count Sigma` (cnv_metrics.csv)       | value > max      |
| `num_segments_max`              | `Number of segments` (cnv_metrics.csv)                       | value > max      |
| `filtered_records_relative_max` | `Number of filtered records (total)` percentage (cnv_metrics.csv) | value > max |
| `homozygosity_index_max`        | `HomozygosityIndex` (CNV VCF header)                         | value > max      |
| `total_excluded_bp_max`         | total bp in `*.cnv.excluded_intervals.bed.gz` (`total_excluded_bp`) | value > max |

The limits themselves pass (e.g. a purity of exactly `purity_min` passes). Purity is a fraction: `--purity_min 0.4` means 40 %. A value outside [0, 1] (e.g. `40`) stops the pipeline with an `ERROR`, as does any non-numeric threshold or an unknown option passed to the script.

### VAF histogram

With `--enable_vaf_histogram` (default `false`), `bin/vaf_histogram.py` counts the variants per variant allele fraction (VAF) bin in the VCF of every sample that has a `vcf` in the samplesheet, and adds the section "Variant allele fraction" to the MultiQC report: one histogram per sample (switch with the buttons above the plot), x-axis VAF, y-axis number of variants. The VCF can come from any variant caller. Samples without `vcf` (e.g. the normal of a tumor/normal pair whose VCF is given in the tumor row) get no histogram; if no sample has a `vcf`, the pipeline stops with an `ERROR` (remove `--enable_vaf_histogram` to run without).

| Parameter                        | Default               | Meaning                                                                                                                                                             |
| -------------------------------- | --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `enable_vaf_histogram`           | `false`               | Run the VAF histogram; needs the samplesheet column `vcf`.                                                                                                          |
| `vaf_histogram_allele_fraction`  | `auto`                | VAF source: `FORMAT_AF` (FORMAT/AF), `FORMAT_AD` (FORMAT/AD[alt] / FORMAT/DP) or `auto` (`FORMAT_AF` if the VCF header declares FORMAT/AF, else `FORMAT_AD`).       |
| `vaf_histogram_bin_width`        | `0.05`                | Bin width, a fraction in (0, 1] (`0.05` = 5 % bins). Bins are [start, end); the last bin includes VAF 1.0.                                                           |
| `vaf_histogram_sample_regex`     | not set (all columns) | Count only VCF columns whose sample name matches this Python regular expression, e.g. `-T$` for the tumor column of a tumor/normal VCF. Stops if a VCF has no match. |
| `vaf_histogram_include_filtered` | `false`               | Also count records whose FILTER is not `PASS` or `.`.                                                                                                               |

The VAF is determined as in [cio-abcd/variantinterpretation](https://github.com/cio-abcd/variantinterpretation) (`--allele_fraction`): from FORMAT/AF, or as FORMAT/AD of the ALT allele divided by FORMAT/DP. With `auto` the source is chosen per VCF from its `##FORMAT` header lines and written to `vaf_histogram_summary.tsv`. There is no silent fallback: the pipeline stops with an `ERROR` naming file, line, position and column if

- the VCF header declares neither FORMAT/AF nor FORMAT/AD and FORMAT/DP (or not the fields of the forced source),
- a counted variant has no value (`.` or missing) in the field(s) of the source, FORMAT/DP is 0, or the VAF is outside [0, 1] (e.g. AD[alt] > DP),
- the VCF is malformed or cannot be read.

Each ALT allele of a multi-allelic record counts as one variant (as after `bcftools norm -m-any`); ALT `.`, `*`, `<NON_REF>` and `<*>` are no variants and are skipped. If the VCF has FORMAT/GT, only the ALT alleles in the sample's genotype are counted; alleles not in the GT (hom-ref `0/0`, no-call `./.`, e.g. the normal column of a tumor/normal VCF) are skipped and counted in `alleles_not_in_genotype`. Records with FILTER `.` (no filters applied, e.g. FreeBayes or unfiltered GATK output) are counted like `PASS`; their number is logged as a `WARNING` and written to `records_filter_missing`. With one counted VCF column the histogram is named after the samplesheet sample, with several columns `<sample>_<column>`.

Tested with synthetic VCFs in the formats of DRAGEN (tumor/normal, `GT:AD:AF:DP`), Mutect2 (tumor/normal, `GT:AD:AF:DP:F1R2:F2R1`) and FreeBayes (germline, `GT:DP:AD:RO:AO`, no FORMAT/AF) in `tests/fixtures/vaf_histogram/`, Python 3.12, MultiQC 1.35. A Strelka2-like VCF (no FORMAT/AF or AD) stops with an `ERROR`.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull core-unit-bioinformatics/exomeqc
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [core-unit-bioinformatics/exomeqc releases page](https://github.com/core-unit-bioinformatics/exomeqc/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Test profile

`-profile test` runs the pipeline on a small public dataset, downloaded on the fly, so no
other parameters are needed:

```bash
nextflow run core-unit-bioinformatics/exomeqc -profile test,<docker/singularity/conda> --outdir results_test
```

The test data are SEQC2 reference samples HCC1395 (tumor) and HCC1395BL (matched normal),
one pair each for **gene panel** (AmpliSeq), **exome** (WES) and **genome** (WGS), sliced to
`chr2:29700000-29740000` (GRCh38, ALK). With them come a FastQ Screen database (human region
and PhiX, `.tar.gz`) and a Picard haplotype map of common 1000 Genomes SNPs in that region.

| Step              | In the test profile | Expected result                                                              |
| ----------------- | ------------------- | ---------------------------------------------------------------------------- |
| FastQ Screen      | on                  | 97-100 % of reads hit `Human`, none hit `PhiX`                                |
| Fingerprint check | on                  | see below                                                                    |
| DRAGEN CNV QC     | off                 | not run; SEQC2 has no DRAGEN CNV output                                       |
| MultiQC           | always              | report with FastQ Screen and Picard fingerprint sections                     |

Fingerprint results on the test data (reference run with conda, 2026-11-17) show what the
check can and cannot tell on small regions:

- **Panel samples: `INCONCLUSIVE`** (LOD 0), even against themselves. None of the 24
  haplotype-map SNPs fall inside the AmpliSeq amplicons. This is the panel restriction
  described under [Fingerprint check](#fingerprint-check).
- **Same tissue across assays matches:** tumor exome vs tumor genome LOD +2.1, normal exome
  vs normal genome LOD +0.4.
- **Tumor vs normal is negative** (genome: LOD -41, tumor-aware -5), although both come from
  the same donor. HCC1395 carries allelic changes in this region that the tumor-aware model
  does not fully absorb. Picard labels these pairs `EXPECTED_MISMATCH` only because their
  read-group sample names differ.

The test profile checks that the steps run. It is not a benchmark of identity calling.

The data live on the `test-datasets` branch of the pipeline repository
(`params.pipelines_testdata_base_path`). How they were made, with the script to rebuild them,
is described in `seqc2/README.md` there. For offline/HPC use, download that folder, point
`--pipelines_testdata_base_path` to its parent directory, and rewrite the URLs in
`seqc2/samplesheet.csv` to local paths.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing, see [Test profile](#test-profile)
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

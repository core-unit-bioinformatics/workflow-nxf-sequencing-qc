# core-unit-bioinformatics/workflow-nxf-sequencing-qc

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/core-unit-bioinformatics/workflow-nxf-sequencing-qc)
[![GitHub Actions CI Status](https://github.com/core-unit-bioinformatics/workflow-nxf-sequencing-qc/actions/workflows/nf-test.yml/badge.svg)](https://github.com/core-unit-bioinformatics/workflow-nxf-sequencing-qc/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/core-unit-bioinformatics/workflow-nxf-sequencing-qc/actions/workflows/linting.yml/badge.svg)](https://github.com/core-unit-bioinformatics/workflow-nxf-sequencing-qc/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.2.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.2.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/core-unit-bioinformatics/workflow-nxf-sequencing-qc)

## Introduction

**core-unit-bioinformatics/workflow-nxf-sequencing-qc** collects quality control for
**short-read** sequencing runs in one MultiQC report. It does not depend on the assay type:
**gene panels, exomes and genomes** run the same way, with the same samplesheet and parameters.

The pipeline does not align reads (except optionally for the fingerprint check) or call variants.
It runs QC tools on the FASTQ, BAM and VCF files you already have, and it collects metrics from
existing pipeline output (for example DRAGEN, Picard or mosdepth) into one report.

## Features

| Step                                    | Tool                                                                                                                                                                                                    | Input                                     | Assay types            | Default |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------- | ---------------------- | ------- |
| Sample identity / swap check            | [Picard CrosscheckFingerprints](https://broadinstitute.github.io/picard/) within the run                                                                                                                | BAM                                       | panel\*, exome, genome | on      |
| Fingerprint QC and fingerprint database | [Picard ExtractFingerprint](https://broadinstitute.github.io/picard/), [samtools](https://www.htslib.org/), [bcftools](https://samtools.github.io/bcftools/), [PyVCF3](https://github.com/dridk/PyVCF3) | BAM + database of previous runs           | panel\*, exome, genome | on      |
| Contamination / species screen          | [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)                                                                                                                        | FASTQ (single- or paired-end)             | panel, exome, genome   | on      |
| CNV QC (purity, ploidy, noise, plots)   | `bin/analyze_cnv_metrics.R` + [karyoploteR](https://bioconductor.org/packages/karyoploteR/)                                                                                                             | **DRAGEN** CNV caller output              | exome, genome\*\*      | on      |
| VAF histogram                           | `bin/vaf_histogram.py`                                                                                                                                                                                  | small-variant VCF (any caller)            | panel, exome, genome   | off     |
| Report summary (CNV QC, TMB, MSI, HRD)  | `bin/exomeqc_summary.py`                                                                                                                                                                                | **DRAGEN** output, CNV and fingerprint QC | panel, exome, genome   | on      |
| Aggregated report                       | [MultiQC](https://multiqc.info/)                                                                                                                                                                        | all of the above + extra directories      | panel, exome, genome   | always  |

\* see [Restrictions](#restrictions). \*\* any assay DRAGEN somatic CNV calling supports.

MultiQC also reads any directory listed in `--multiqc_inputs` and the DRAGEN directory
given by `--dragen_dirs`. That way, DRAGEN mapping, coverage (WGS **and** target-region
coverage) and variant-calling metrics, and the output of any other
[MultiQC-supported tool](https://multiqc.info/modules/), show up in the same report.

The steps in detail:

1. Fingerprint check (`--enable_fingerprintcheck`)
   - crosscheck all BAMs of the run with each other ([`Picard CrosscheckFingerprints`](https://broadinstitute.github.io/picard/)), shown as heatmap
   - optionally align samples without BAM ([`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2), `--align_reads_for_fingerprint`)
   - extract one fingerprint per sample and run (`<sample>@<runid>`, [`Picard ExtractFingerprint`](https://broadinstitute.github.io/picard/))
   - fingerprint QC: depth and base quality of each haplotype site with adjustable thresholds ([`samtools mpileup`](https://www.htslib.org/)), optional contamination estimate (off by default, not suited for tumours), QC status per sample; failed sites are masked and status, thresholds and values are stored in the fingerprint ([`bcftools`](https://samtools.github.io/bcftools/))
   - check the fingerprint database of previous runs (`--fingerprint_database`, plain VCF, gzip or bgzip) against the haplotype map ([`PyVCF3`](https://github.com/dridk/PyVCF3)), drop the fingerprints this run replaces (same `<sample>@<runid>`) and convert it to bgzip with index ([`bcftools`](https://samtools.github.io/bcftools/))
   - crosscheck the run against the database, before the run is added, and report matches with samples of other runs, with LOD scores
   - add the run to the database: updated database `exomeqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz` for the next run
   - overview in MultiQC: haplotype map, database, replaced fingerprints and whether each sample gave a suitable fingerprint
2. Contamination screening ([`FastQ Screen`](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), `--enable_fastqscreen`)
3. DRAGEN CNV QC (`--enable_cnvqc`): QC status per sample from DRAGEN CNV metrics with adjustable thresholds, purity/coverage model plot and karyotype plots ([`karyoploteR`](https://bioconductor.org/packages/karyoploteR/))
4. VAF histogram (`--enable_vaf_histogram`): histogram of the variant allele fractions per sample from the samplesheet column `vcf`
5. Report summary (`--enable_summary`, needs `--dragen_dirs`): one row per DRAGEN case with CNV QC status, purity, ploidy, TMB, MSI and HRD on top of the report
6. Present QC ([`MultiQC`](http://multiqc.info/)), optionally named after the run (`--multiqc_report_name`)

## Restrictions

- **Short reads only.** The pipeline is built and tested for Illumina-type short reads.
  Long-read data (ONT, PacBio) is not supported.
- **CNV QC is DRAGEN-specific.** `--enable_cnvqc` only parses DRAGEN CNV caller output
  (`*.cnv_metrics.csv`, `*.cnv.vcf.gz`, and if present `*.tn.tsv.gz`, `*.tumor.baf.bedgraph.gz`,
  `*.cnv.purity.coverage.models.tsv`, `*.cnv.excluded_intervals.bed.gz`). Output from other CNV
  callers (CNVkit, GATK gCNV, ...) is not supported. Turn the step off for non-DRAGEN data
  (`enable_cnvqc: false`). The CNV thresholds have no built-in defaults and need to be calibrated
  for each assay (panel, exome or genome). The report summary also reads DRAGEN output only.
- **Fingerprinting needs BAM input.** CRAM is not supported yet; the pipeline stops with an
  error if it gets CRAM files. The haplotype map must match the reference build (`--fasta`), and
  the fingerprint database must be built with the same haplotype map. Fingerprints only work where
  the map's loci are covered: for **small gene panels**, use a haplotype map built from SNPs inside
  the panel targets, otherwise the results are `INCONCLUSIVE` and the fingerprint QC fails.
- **FastQ Screen needs its own database.** Build it in advance (`fastq_screen --get_genomes`
  or your own Bowtie2 indices). The pipeline never downloads anything.
- **Offline / HPC.** No resources are downloaded at runtime. All databases and references
  are passed in as parameters.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test_seqc2` (see [Test profile](#test-profile); `-profile test` is currently disabled) before running the workflow on actual data.

Prepare a samplesheet (one row per sample and sequencing run):

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,bam,runid,vcf
PANEL_01,PANEL_01_R1.fastq.gz,PANEL_01_R2.fastq.gz,PANEL_01.bam,20260101_LH00001_0001_AXXXXXXXX,PANEL_01.vcf.gz
WGS_01,WGS_01_R1.fastq.gz,WGS_01_R2.fastq.gz,WGS_01.bam,20260101_LH00001_0001_AXXXXXXXX,
SE_01,SE_01_R1.fastq.gz,,SE_01.bam,20260101_LH00001_0001_AXXXXXXXX,
```

- `fastq_1`/`fastq_2` are needed for every row when FastQ Screen is on. Leave `fastq_2`
  empty for single-end data.
- `bam` (BAM, not CRAM) and `runid` are needed for every row when the fingerprint check is on.
- `vcf` is optional and only used by the VAF histogram.

Put the parameters in a params file:

```yaml title="params.yaml"
input: samplesheet.csv
outdir: results
fastqscreen_database: /path/to/FastQ_Screen_Genomes
fingerprint_haplotype_map: /path/to/haplotype_map.txt
fasta: /path/to/genome.fasta
# fingerprint database of the previous run (optional)
fingerprint_database: /path/to/previous/results/fingerprintcheck/database/exomeqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz
# DRAGEN only; leave out and set enable_cnvqc: false for other callers
dragen_dirs: /path/to/dragen_output
multiqc_report_name: 20260101_LH00001_0001_AXXXXXXXX
```

Then run:

```bash
nextflow run core-unit-bioinformatics/workflow-nxf-sequencing-qc \
   -profile <docker/singularity/.../institute> \
   -params-file params.yaml
```

> [!IMPORTANT]
> Switch steps off in the params file (for example `enable_cnvqc: false`).
> `--enable_cnvqc false` on the command line is passed as a string and fails parameter validation.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For all parameters see [docs/usage.md](docs/usage.md); for the output files see [docs/output.md](docs/output.md).

## Test profile

`-profile test` and `-profile test_full` are **disabled** for now (no test dataset with direct
download links yet); the pipeline stops with an error when they are used. Use `test_seqc2`.

The `test_seqc2` profile uses SEQC2 HCC1395/HCC1395BL tumor/normal data (gene panel, exome, genome).
**The test data must be downloaded first**; Nextflow does not fetch it:

```bash
conda env create -f scripts/envs/seqc2_testdata.yaml && conda activate seqc2_testdata
scripts/download_seqc2_testdata.py --outdir testdata      # ~4 MB into ./testdata/seqc2/
nextflow run core-unit-bioinformatics/workflow-nxf-sequencing-qc -profile test_seqc2,<docker/singularity/conda> --outdir results_test
```

Run `nextflow` from the directory that contains `testdata/`, or pass
`--pipelines_testdata_base_path /path/to/testdata/`. Details and expected results:
[docs/usage.md](docs/usage.md#test-profile).

## Credits

core-unit-bioinformatics/workflow-nxf-sequencing-qc was originally written by Kai Horny
(previously named `exomeqc`).

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use core-unit-bioinformatics/workflow-nxf-sequencing-qc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

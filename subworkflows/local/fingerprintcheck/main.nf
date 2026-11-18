/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: FINGERPRINTCHECK
    1. Crosscheck the fingerprints of all samples of this run (Picard
       CrosscheckFingerprints on the BAMs, heatmap in MultiQC).
    2. Extract one fingerprint VCF per sample (Picard ExtractFingerprint) and
       check its quality: each haplotype site must reach adjustable depth and
       base quality thresholds (samtools mpileup, fingerprint_qc.py), failed
       sites are masked (flat genotype likelihoods), the contamination is
       estimated, and the sample gets a QC status. Site values go into the
       fingerprint VCF (bcftools annotate), the status with thresholds and
       values into its ##SAMPLE header line.
       If a fingerprint database of previous runs is given: check it against
       the haplotype map, drop the fingerprints this run replaces, convert it
       to bgzip with index, and crosscheck this run against it. Then add this
       run's fingerprints to it: the updated database, with a timestamp in
       the file name and one ##SAMPLE QC line per fingerprint. Report the fingerprint check and the matches with
       samples of other runs.
    Samples without a BAM are aligned with bwa-mem2 if align_reads is true.
    Fingerprint IDs are <sample>@<runid>, so that the same sample sequenced in
    several runs keeps one fingerprint per run.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX                                                         } from '../../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY                                        } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { BWAMEM2_INDEX                                                          } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                                                            } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX                                                         } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_CROSSCHECKFINGERPRINTS                                          } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'
include { PICARD_EXTRACTFINGERPRINT                                              } from '../../../modules/nf-core/picard/extractfingerprint/main'
include { CHECK_FINGERPRINT_DATABASE                                             } from '../../../modules/local/check_fingerprint_database/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FINGERPRINT_DATABASE                    } from '../../../modules/nf-core/bcftools/view/main'
include { PICARD_CROSSCHECKFINGERPRINTS as PICARD_CROSSCHECKFINGERPRINTS_DATABASE } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_FINGERPRINT_DATABASE                  } from '../../../modules/nf-core/bcftools/merge/main'
include { FINGERPRINT_DATABASE_REPORT                                            } from '../../../modules/local/fingerprint_database_report/main'
include { GAWK as GAWK_FINGERPRINT_SITES                                         } from '../../../modules/nf-core/gawk/main'
include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_FINGERPRINT                       } from '../../../modules/nf-core/samtools/mpileup/main'
include { FINGERPRINT_QC                                                         } from '../../../modules/local/fingerprint_qc/main'
include { HTSLIB_BGZIPTABIX as HTSLIB_BGZIPTABIX_FINGERPRINT_QC                  } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_FINGERPRINT_QC                  } from '../../../modules/nf-core/bcftools/annotate/main'
include { FINGERPRINT_DATABASE_HEADER                                            } from '../../../modules/local/fingerprint_database_header/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_FINGERPRINT_DATABASE            } from '../../../modules/nf-core/bcftools/reheader/main'

workflow FINGERPRINTCHECK {

    take:
    ch_samplesheet       // channel: [ val(meta), [ fastq ] ], meta has id, single_end, bam, runid
    haplotype_map        // path: Picard haplotype map
    fasta                // path: reference FASTA the BAMs are aligned to
    fasta_fai            // path or null: FASTA index, created with samtools faidx if null
    fasta_dict           // path or null: sequence dictionary, created with Picard if null
    bwamem2_index        // path or null: bwa-mem2 index directory, built from fasta if null
    align_reads          // boolean: align samples without BAM with bwa-mem2
    fingerprint_database // path or null: fingerprint database of previous runs (VCF, plain, gzip or bgzip)

    main:
    ch_versions = channel.empty()

    // sample and runid become VCF sample names and file names: CUBI file name characters only
    def name_pattern = ~/^[A-Za-z0-9._-]+$/
    def haplotype_map_file = file(haplotype_map, checkIfExists: true)
    def fasta_file = file(fasta, checkIfExists: true)
    // the updated database is named with the start time of this pipeline run
    def database_timestamp = workflow.start.format(java.time.format.DateTimeFormatter.ofPattern('yyyyMMdd-HHmmss'))

    // -------------------------------------------------------------------------
    // Validate every samplesheet row and build the fingerprint meta map:
    // [ id, runid, fingerprint_id, single_end ]
    // -------------------------------------------------------------------------
    ch_samples = ch_samplesheet
        .map { meta, fastqs ->
            // the fingerprint ID needs the RUNID, so every row must have one
            if (!meta.runid) {
                error("FINGERPRINTCHECK: samplesheet row '${meta.id}': column 'runid' is empty; it is required with --enable_fingerprintcheck (e.g. 20260101_LH00001_0001_AFIXTURE01)")
            }
            // reject sample names that would be unsafe as VCF sample or file names
            if (!(meta.id ==~ name_pattern)) {
                error("FINGERPRINTCHECK: samplesheet row '${meta.id}': sample may only contain A-Z, a-z, 0-9, '.', '_' and '-' with --enable_fingerprintcheck")
            }
            // samples without BAM can only be fingerprinted after alignment
            if (!meta.bam && !align_reads) {
                error("FINGERPRINTCHECK: samplesheet row '${meta.id}' (run ${meta.runid}): column 'bam' is empty; provide a BAM or set --align_reads_for_fingerprint to align its FASTQ files with bwa-mem2")
            }
            // Picard CrosscheckFingerprints and samtools mpileup get no reference index here
            if (meta.bam && meta.bam.toString().endsWith('.cram')) {
                error("FINGERPRINTCHECK: samplesheet row '${meta.id}' (run ${meta.runid}): CRAM input '${meta.bam}' is not supported for fingerprinting (expected BAM); convert it to BAM or set --enable_fingerprintcheck false")
            }
            // alignment takes exactly one single-end or one paired-end FASTQ set
            def expected_fastqs = meta.single_end ? 1 : 2
            if (!meta.bam && fastqs.size() != expected_fastqs) {
                error("FINGERPRINTCHECK: samplesheet row '${meta.id}' (run ${meta.runid}): ${fastqs.size()} FASTQ file(s) found, expected ${expected_fastqs} for alignment; provide one fastq_1/fastq_2 pair per sample (concatenate lanes) or a BAM")
            }
            def fingerprint_meta = [
                id            : meta.id,
                runid         : meta.runid,
                fingerprint_id: "${meta.id}@${meta.runid}".toString(),
                single_end    : meta.single_end
            ]
            [ fingerprint_meta, meta.bam, fastqs ]
        }
        .branch { meta, bam, fastqs ->
            // samples with a BAM/CRAM from the samplesheet
            bam: bam
                return [ meta, file(bam.toString(), checkIfExists: true) ]
            // samples to align from FASTQ
            align: true
                return [ meta, fastqs ]
        }

    // -------------------------------------------------------------------------
    // Reference: FASTA with index and sequence dictionary (created if missing)
    // -------------------------------------------------------------------------
    ch_fasta = channel.value([ [ id: 'reference' ], fasta_file ])

    if (fasta_fai) {
        ch_fai = channel.value(file(fasta_fai, checkIfExists: true))
    } else {
        // index the FASTA with samtools faidx
        SAMTOOLS_FAIDX(ch_fasta.map { meta, fa -> [ meta, fa, [] ] }, false)
        // value channel: the input is the value channel ch_fasta
        ch_fai = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> fai }
    }

    if (fasta_dict) {
        ch_dict = channel.value(file(fasta_dict, checkIfExists: true))
    } else {
        // create the sequence dictionary with Picard CreateSequenceDictionary
        PICARD_CREATESEQUENCEDICTIONARY(ch_fasta)
        // value channel: the input is the value channel ch_fasta
        ch_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict.map { _meta, dict -> dict }
    }

    // -------------------------------------------------------------------------
    // Alignment of samples without BAM (only with --align_reads_for_fingerprint)
    // -------------------------------------------------------------------------
    if (bwamem2_index) {
        ch_bwamem2_index = channel.value([ [ id: 'bwamem2' ], file(bwamem2_index, checkIfExists: true) ])
    } else {
        // build the index only if at least one sample has to be aligned
        BWAMEM2_INDEX(
            ch_samples.align
                .first()
                .combine(ch_fasta)
                .map { _meta, _fastqs, meta_fasta, fa -> [ meta_fasta, fa ] }
        )
        ch_bwamem2_index = BWAMEM2_INDEX.out.index.first()
    }

    // align, then sort the reads by coordinate (sort_bam = true)
    BWAMEM2_MEM(
        ch_samples.align,
        ch_bwamem2_index,
        ch_fasta,
        true
    )

    // -------------------------------------------------------------------------
    // Index all BAMs (from the samplesheet and from alignment)
    // -------------------------------------------------------------------------
    ch_bams_to_index = ch_samples.bam.mix(BWAMEM2_MEM.out.bam)

    SAMTOOLS_INDEX(ch_bams_to_index)

    // join each BAM with its index by the meta map; every BAM needs exactly one index
    ch_bam_with_index = ch_bams_to_index
        .join(SAMTOOLS_INDEX.out.index, by: 0, failOnMismatch: true, failOnDuplicate: true)
        // channel: [ meta, bam, index ]

    // -------------------------------------------------------------------------
    // 1. Crosscheck all BAMs of this run with each other (heatmap in MultiQC).
    // Build exactly one Picard tuple for the whole cohort:
    // [meta, bam_list, index_list, input2_list, input2_index_list, haplotype_map]
    //
    // collect(flat:false) is critical here; it preserves tuple boundaries so we
    // get rows shaped as [meta, bam, index].
    // -------------------------------------------------------------------------
    ch_fingerprint_input = ch_bam_with_index
        .collect(flat: false)
        .map { rows ->
            if (!rows) {
                error("FINGERPRINTCHECK: no BAM/index rows were collected")
            }

            if (rows.any { row -> row.size() != 3 }) {
                error("FINGERPRINTCHECK: expected [meta, bam, index] rows, observed: ${rows}")
            }

            def bams = rows.collect { row -> row[1] }
            def indexes = rows.collect { row -> row[2] }

            if (bams.size() != indexes.size()) {
                error("FINGERPRINTCHECK: BAM/index list size mismatch (${bams.size()} vs ${indexes.size()})")
            }

            [
                [ id: 'fingerprintcheck' ],
                bams,
                indexes,
                [],   // input2 (SECOND_INPUT) — unused
                [],   // input2_index           — unused
                haplotype_map_file
            ]
        }

    // Picard allows the reference input to be empty.
    ch_fingerprint_reference = channel
        .value([[:], []])
        .map { metaRef, fastaRef ->
            [ metaRef, fastaRef ]
        }

    PICARD_CROSSCHECKFINGERPRINTS(
        ch_fingerprint_input,
        ch_fingerprint_reference
    )

    // -------------------------------------------------------------------------
    // 2a. Extract one fingerprint VCF per sample, sample name <sample>@<runid>
    // -------------------------------------------------------------------------
    PICARD_EXTRACTFINGERPRINT(
        ch_bam_with_index,
        haplotype_map_file,
        fasta_file,
        ch_fai,
        ch_dict
    )

    // join each fingerprint VCF with its tabix index: [ meta, vcf, tbi ]
    ch_fingerprints_raw = PICARD_EXTRACTFINGERPRINT.out.vcf
        .join(PICARD_EXTRACTFINGERPRINT.out.tbi, by: 0, failOnMismatch: true, failOnDuplicate: true)

    // -------------------------------------------------------------------------
    // 2b. Fingerprint QC of each sample
    // -------------------------------------------------------------------------
    // BED of all SNPs of the haplotype map (value channel: the input is a value)
    GAWK_FINGERPRINT_SITES(
        channel.value([ [ id: 'haplotype_map_sites' ], haplotype_map_file, 'bed' ]),
        [],
        false
    )
    ch_sites_bed = GAWK_FINGERPRINT_SITES.out.output.map { _meta, bed -> bed }

    // every base at the map SNPs with base and mapping quality and read name
    ch_fasta_with_fai = ch_fasta
        .combine(ch_fai)
        .map { meta, fa, fai -> [ meta, fa, fai ] }
        .first()

    SAMTOOLS_MPILEUP_FINGERPRINT(
        ch_bam_with_index
            .combine(ch_sites_bed)
            .map { meta, bam, index, bed -> [ meta, bam, index, bed ] },
        ch_fasta_with_fai
    )

    // site QC, contamination and QC status: [ meta, fingerprint vcf, pileup ]
    FINGERPRINT_QC(
        ch_fingerprints_raw
            .join(SAMTOOLS_MPILEUP_FINGERPRINT.out.mpileup, by: 0, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, _tbi, pileup -> [ meta, vcf, pileup ] },
        haplotype_map_file
    )

    // bgzip and tabix-index the site table for bcftools annotate
    HTSLIB_BGZIPTABIX_FINGERPRINT_QC(
        FINGERPRINT_QC.out.sites.map { meta, sites -> [ meta, sites, [], [] ] },
        'compress',
        true,
        ''
    )

    // write the QC values into the fingerprint, mask failed sites (flat PL)
    // and add the ##SAMPLE QC line:
    // [ meta, vcf, tbi, annotations, annotations_index, columns, header_lines, rename_chrs ]
    BCFTOOLS_ANNOTATE_FINGERPRINT_QC(
        ch_fingerprints_raw
            .join(HTSLIB_BGZIPTABIX_FINGERPRINT_QC.out.output, by: 0, failOnMismatch: true, failOnDuplicate: true)
            .join(HTSLIB_BGZIPTABIX_FINGERPRINT_QC.out.index, by: 0, failOnMismatch: true, failOnDuplicate: true)
            .join(FINGERPRINT_QC.out.header, by: 0, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi, sites, sites_index, header -> [ meta, vcf, tbi, sites, sites_index, [], header, [] ] }
    )

    // fingerprints after QC, used for the database: [ meta, vcf, tbi ]
    ch_fingerprints = BCFTOOLS_ANNOTATE_FINGERPRINT_QC.out.vcf
        .join(BCFTOOLS_ANNOTATE_FINGERPRINT_QC.out.index, by: 0, failOnMismatch: true, failOnDuplicate: true)

    // all fingerprints of this run as two lists: [ [ vcf, ... ], [ tbi, ... ] ]
    ch_fingerprint_lists = ch_fingerprints
        .map { _meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { rows -> [ rows.collect { row -> row[0] }, rows.collect { row -> row[1] } ] }

    // this run's samples as TSV for the report; the BAM file name links each
    // sample to its rows in the crosscheck of step 1 (CROSSCHECK_BY FILE)
    ch_run_samples = ch_bam_with_index
        .map { meta, bam, _index -> "${meta.fingerprint_id}\t${meta.id}\t${meta.runid}\t${bam.name}" }
        .collectFile(name: 'run_samples.tsv', seed: "fingerprint_id\tsample\trunid\tbam", newLine: true, sort: true)

    if (fingerprint_database) {
        def database_file = file(fingerprint_database, checkIfExists: true)

        // ---------------------------------------------------------------------
        // 2c. Check the given database against the haplotype map and list the
        // fingerprints this run replaces (same <sample>@<runid>)
        // ---------------------------------------------------------------------
        // this run's fingerprint IDs are known from the samplesheet, so the
        // check runs before any alignment or fingerprint extraction
        ch_run_fingerprint_ids = ch_samples.bam
            .mix(ch_samples.align)
            .map { meta, _files -> meta.fingerprint_id }
            .collectFile(name: 'run_fingerprint_ids.txt', newLine: true, sort: true)

        CHECK_FINGERPRINT_DATABASE(
            channel.value([ [ id: 'previous_database' ], database_file ]),
            haplotype_map_file,
            ch_run_fingerprint_ids
        )

        // ---------------------------------------------------------------------
        // 2d. Keep the listed fingerprints and write them as bgzipped VCF with
        // tabix index; bcftools reads plain VCF, gzip and bgzip, and needs no
        // index here. Skipped if the database holds only this run's IDs.
        // ---------------------------------------------------------------------
        ch_database_to_keep = CHECK_FINGERPRINT_DATABASE.out.keep
            .filter { _meta, keep -> keep.size() > 0 }
            .map { meta, keep -> [ meta, database_file, keep ] }

        BCFTOOLS_VIEW_FINGERPRINT_DATABASE(
            ch_database_to_keep.map { meta, database, _keep -> [ meta, database, [] ] },
            [],
            [],
            ch_database_to_keep.map { _meta, _database, keep -> keep }
        )

        // kept fingerprints with their index: [ meta, vcf.gz, tbi ]
        ch_previous_kept = BCFTOOLS_VIEW_FINGERPRINT_DATABASE.out.vcf
            .join(BCFTOOLS_VIEW_FINGERPRINT_DATABASE.out.index, by: 0, failOnMismatch: true)

        // ---------------------------------------------------------------------
        // 2e. Crosscheck this run's fingerprints against the kept database
        // fingerprints, before this run is added to the database
        // ---------------------------------------------------------------------
        ch_crosscheck_database_input = ch_fingerprint_lists
            .combine(ch_previous_kept)
            .map { vcfs, tbis, _meta, database, database_index ->
                [ [ id: 'fingerprint_database' ], vcfs, tbis, [ database ], [ database_index ], haplotype_map_file ]
            }

        PICARD_CROSSCHECKFINGERPRINTS_DATABASE(
            ch_crosscheck_database_input,
            ch_fingerprint_reference
        )

        // inputs of the merge and the report; [] where nothing was kept
        ch_previous_for_merge = ch_previous_kept
            .map { _meta, vcf, tbi -> [ [ vcf ], [ tbi ] ] }
            .ifEmpty([ [], [] ])
        ch_database_crosscheck = PICARD_CROSSCHECKFINGERPRINTS_DATABASE.out.crosscheck_metrics
            .map { _meta, metrics -> [ metrics ] }
            .ifEmpty([ [] ])
        ch_database_check = CHECK_FINGERPRINT_DATABASE.out.check
            .join(CHECK_FINGERPRINT_DATABASE.out.samples, by: 0, failOnMismatch: true)
            .map { _meta, check, samples -> [ check, samples ] }
        ch_database_check_results = CHECK_FINGERPRINT_DATABASE.out.check
        ch_database_crosscheck_metrics = PICARD_CROSSCHECKFINGERPRINTS_DATABASE.out.crosscheck_metrics
    } else {
        // no database of previous runs: the updated database is this run only
        ch_previous_for_merge = channel.value([ [], [] ])
        ch_database_crosscheck = channel.value([ [] ])
        ch_database_check = channel.value([ [], [] ])
        ch_database_check_results = channel.empty()
        ch_database_crosscheck_metrics = channel.empty()
    }

    // -------------------------------------------------------------------------
    // 2f. Add this run's fingerprints to the kept database fingerprints:
    // updated database sequencingqc_fingerprint_database_<YYYYMMDD-hhmmss>.vcf.gz
    // -------------------------------------------------------------------------
    ch_merge_input = ch_fingerprint_lists
        .combine(ch_previous_for_merge)
        .map { vcfs, tbis, previous_vcfs, previous_tbis ->
            // bcftools/merge orders the files by name, which sets the sample order
            [ [ id: "sequencingqc_fingerprint_database_${database_timestamp}".toString() ], previous_vcfs + vcfs, previous_tbis + tbis, [] ]
        }

    BCFTOOLS_MERGE_FINGERPRINT_DATABASE(ch_merge_input, [ [:], [], [] ])

    // merged database: [ meta, vcf ]
    ch_database_merged = BCFTOOLS_MERGE_FINGERPRINT_DATABASE.out.vcf

    // one ##SAMPLE QC line per fingerprint: bcftools merge and view keep
    // outdated lines of replaced fingerprints, so the header is rewritten
    FINGERPRINT_DATABASE_HEADER(
        ch_database_merged,
        FINGERPRINT_QC.out.header.map { _meta, header -> header }.collect()
    )

    BCFTOOLS_REHEADER_FINGERPRINT_DATABASE(
        ch_database_merged
            .join(FINGERPRINT_DATABASE_HEADER.out.header, by: 0, failOnMismatch: true)
            .map { meta, vcf, header -> [ meta, vcf, header, [] ] },
        [ [:], [] ]
    )

    // updated database with its index: [ meta, vcf, tbi ]
    ch_database = BCFTOOLS_REHEADER_FINGERPRINT_DATABASE.out.vcf
        .join(BCFTOOLS_REHEADER_FINGERPRINT_DATABASE.out.index, by: 0, failOnMismatch: true)

    // -------------------------------------------------------------------------
    // 2g. Report the fingerprint QC, the fingerprint check and the matches
    // with other runs
    // -------------------------------------------------------------------------
    // all QC summaries of this run as one list
    ch_qc_summaries = FINGERPRINT_QC.out.summary
        .map { _meta, summary -> summary }
        .collect()
        .map { summaries -> [ summaries ] }

    ch_report_input = ch_fingerprint_lists
        .combine(ch_run_samples)
        .combine(PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics)
        .combine(ch_database)
        .combine(ch_database_check)
        .combine(ch_database_crosscheck)
        .combine(ch_qc_summaries)
        .map { vcfs, _tbis, run_samples, _meta_run, run_crosscheck, _meta_database, database, _database_index, database_check, database_samples, database_crosscheck, qc_summaries ->
            [ [ id: 'fingerprint_database' ], run_samples, vcfs, haplotype_map_file, run_crosscheck, database, database_check, database_samples, database_crosscheck, qc_summaries ]
        }

    FINGERPRINT_DATABASE_REPORT(ch_report_input)

    emit:
    crosscheck_metrics          = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics          // channel: [ meta, crosscheck_metrics.txt ]
    fingerprints                = ch_fingerprints                                               // channel: [ meta, vcf.gz, tbi ] (after QC)
    fingerprint_qc              = FINGERPRINT_QC.out.summary                                    // channel: [ meta, tsv ]
    database                    = ch_database                                                   // channel: [ meta, vcf.gz, tbi ]
    database_check              = ch_database_check_results                                     // channel: [ meta, tsv ] (only with a given database)
    database_crosscheck_metrics = ch_database_crosscheck_metrics                                // channel: [ meta, crosscheck_metrics.txt ] (only with kept database fingerprints)
    database_matches            = FINGERPRINT_DATABASE_REPORT.out.matches                       // channel: [ meta, tsv ]
    database_fingerprints       = FINGERPRINT_DATABASE_REPORT.out.fingerprints                  // channel: [ meta, tsv ]
    database_summary            = FINGERPRINT_DATABASE_REPORT.out.summary                       // channel: [ meta, tsv ]
    fingerprint_qc_report       = FINGERPRINT_DATABASE_REPORT.out.qc                            // channel: [ meta, tsv ] (all samples of this run)
    database_multiqc            = FINGERPRINT_DATABASE_REPORT.out.overview
        .mix(FINGERPRINT_DATABASE_REPORT.out.qc_multiqc)
        .mix(FINGERPRINT_DATABASE_REPORT.out.multiqc)                                           // channel: [ meta, _mqc.html or _mqc.tsv ]
    versions                    = ch_versions
}

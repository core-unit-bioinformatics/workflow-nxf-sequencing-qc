/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: FINGERPRINTCHECK
    1. Crosscheck the fingerprints of all samples of this run (Picard
       CrosscheckFingerprints on the BAMs, heatmap in MultiQC).
    2. Extract one fingerprint VCF per sample (Picard ExtractFingerprint),
       merge them with the fingerprint database of previous runs into an
       updated database, crosscheck this run against it and report matches
       with samples of other runs.
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
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FINGERPRINT_DATABASE   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_FINGERPRINT_DATABASE } from '../../../modules/nf-core/bcftools/merge/main'
include { PICARD_CROSSCHECKFINGERPRINTS as PICARD_CROSSCHECKFINGERPRINTS_DATABASE } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'
include { FINGERPRINT_DATABASE_REPORT                                            } from '../../../modules/local/fingerprint_database_report/main'

workflow FINGERPRINTCHECK {

    take:
    ch_samplesheet       // channel: [ val(meta), [ fastq ] ], meta has id, single_end, bam, runid
    haplotype_map        // path: Picard haplotype map
    fasta                // path: reference FASTA the BAMs are aligned to
    fasta_fai            // path or null: FASTA index, created with samtools faidx if null
    fasta_dict           // path or null: sequence dictionary, created with Picard if null
    bwamem2_index        // path or null: bwa-mem2 index directory, built from fasta if null
    align_reads          // boolean: align samples without BAM with bwa-mem2
    fingerprint_database // path or null: fingerprint database (.vcf.gz) of previous runs

    main:
    ch_versions = channel.empty()

    // sample and runid become VCF sample names and file names: CUBI file name characters only
    def name_pattern = ~/^[A-Za-z0-9._-]+$/
    def haplotype_map_file = file(haplotype_map, checkIfExists: true)
    def fasta_file = file(fasta, checkIfExists: true)

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
                error("FINGERPRINTCHECK: samplesheet row '${meta.id}' (run ${meta.runid}): column 'bam' is empty; provide a BAM/CRAM or set --align_reads_for_fingerprint to align its FASTQ files with bwa-mem2")
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
    ch_fingerprints = PICARD_EXTRACTFINGERPRINT.out.vcf
        .join(PICARD_EXTRACTFINGERPRINT.out.tbi, by: 0, failOnMismatch: true, failOnDuplicate: true)

    // all fingerprints of this run as two lists: [ [ vcf, ... ], [ tbi, ... ] ]
    ch_fingerprint_lists = ch_fingerprints
        .map { _meta, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
        .map { rows -> [ rows.collect { row -> row[0] }, rows.collect { row -> row[1] } ] }

    // this run's samples as TSV (fingerprint_id, sample, runid) for the report
    ch_run_samples = ch_fingerprints
        .map { meta, _vcf, _tbi -> "${meta.fingerprint_id}\t${meta.id}\t${meta.runid}" }
        .collectFile(name: 'run_samples.tsv', seed: "fingerprint_id\tsample\trunid", newLine: true, sort: true)

    // -------------------------------------------------------------------------
    // 2b. Previous database: drop fingerprints that this run replaces
    // -------------------------------------------------------------------------
    def previous_database_file = null
    if (fingerprint_database) {
        previous_database_file = file(fingerprint_database, checkIfExists: true)
        // bcftools merge needs a bgzipped, indexed VCF
        if (!previous_database_file.name.endsWith('.vcf.gz')) {
            error("FINGERPRINTCHECK: --fingerprint_database '${fingerprint_database}': expected a bgzipped VCF (.vcf.gz), e.g. exomeqc_fingerprint_database.vcf.gz from fingerprintcheck/database/ of a previous run")
        }
        def previous_tbi = file("${previous_database_file}.tbi")
        def previous_csi = file("${previous_database_file}.csi")
        def previous_database_index = previous_tbi.exists() ? previous_tbi : previous_csi.exists() ? previous_csi : null
        if (!previous_database_index) {
            error("FINGERPRINTCHECK: --fingerprint_database '${fingerprint_database}': no index found (expected '${previous_tbi}' or '${previous_csi}'); create one with 'tabix -p vcf ${fingerprint_database}'")
        }

        // exclude every fingerprint ID of this run from the previous database
        // (absent IDs are ignored); the report script lists the replaced ones
        ch_previous_database = ch_fingerprints
            .map { meta, _vcf, _tbi -> meta.fingerprint_id }
            .collect()
            .map { ids ->
                [ [ id: 'previous_fingerprint_database', exclude_ids: ids.sort() ], previous_database_file, previous_database_index ]
            }

        BCFTOOLS_VIEW_FINGERPRINT_DATABASE(ch_previous_database, [], [], [])

        // filtered previous database as one-element lists: [ [ vcf ], [ tbi ] ]
        ch_previous_filtered = BCFTOOLS_VIEW_FINGERPRINT_DATABASE.out.vcf
            .join(BCFTOOLS_VIEW_FINGERPRINT_DATABASE.out.index, by: 0, failOnMismatch: true)
            .map { _meta, vcf, tbi -> [ [ vcf ], [ tbi ] ] }
    } else {
        ch_previous_filtered = channel.value([ [], [] ])
    }

    // -------------------------------------------------------------------------
    // 2c. Merge this run's fingerprints with the previous database
    // -------------------------------------------------------------------------
    ch_merge_input = ch_fingerprint_lists
        .combine(ch_previous_filtered)
        .map { vcfs, tbis, previous_vcfs, previous_tbis ->
            [ [ id: 'exomeqc_fingerprint_database' ], vcfs + previous_vcfs, tbis + previous_tbis, [] ]
        }

    BCFTOOLS_MERGE_FINGERPRINT_DATABASE(ch_merge_input, [ [:], [], [] ])

    // updated database with its index: [ meta, vcf, tbi ]
    ch_database = BCFTOOLS_MERGE_FINGERPRINT_DATABASE.out.vcf
        .join(BCFTOOLS_MERGE_FINGERPRINT_DATABASE.out.index, by: 0, failOnMismatch: true)

    // -------------------------------------------------------------------------
    // 2d. Crosscheck this run's fingerprints against the updated database
    // -------------------------------------------------------------------------
    ch_crosscheck_database_input = ch_fingerprint_lists
        .combine(ch_database)
        .map { vcfs, tbis, _meta, database, database_index ->
            [ [ id: 'fingerprint_database' ], vcfs, tbis, [ database ], [ database_index ], haplotype_map_file ]
        }

    PICARD_CROSSCHECKFINGERPRINTS_DATABASE(
        ch_crosscheck_database_input,
        ch_fingerprint_reference
    )

    // -------------------------------------------------------------------------
    // 2e. Report matches with samples of other runs
    // -------------------------------------------------------------------------
    ch_report_input = PICARD_CROSSCHECKFINGERPRINTS_DATABASE.out.crosscheck_metrics
        .combine(ch_run_samples)
        .combine(ch_database.map { _meta, database, _database_index -> database })
        .map { meta, metrics, run_samples, database ->
            [ meta, metrics, run_samples, database, previous_database_file ?: [] ]
        }

    FINGERPRINT_DATABASE_REPORT(ch_report_input)

    emit:
    crosscheck_metrics          = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics          // channel: [ meta, crosscheck_metrics.txt ]
    fingerprints                = ch_fingerprints                                               // channel: [ meta, vcf.gz, tbi ]
    database                    = ch_database                                                   // channel: [ meta, vcf.gz, tbi ]
    database_crosscheck_metrics = PICARD_CROSSCHECKFINGERPRINTS_DATABASE.out.crosscheck_metrics // channel: [ meta, crosscheck_metrics.txt ]
    database_matches            = FINGERPRINT_DATABASE_REPORT.out.matches                       // channel: [ meta, tsv ]
    database_summary            = FINGERPRINT_DATABASE_REPORT.out.summary                       // channel: [ meta, tsv ]
    database_multiqc            = FINGERPRINT_DATABASE_REPORT.out.multiqc                       // channel: [ meta, _mqc.tsv or _mqc.html ]
    versions                    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { FASTQSCREEN_FASTQSCREEN } from '../modules/nf-core/fastqscreen/fastqscreen/main'
include { UNTAR                   } from '../modules/nf-core/untar/main'
include { DRAGEN_CNV_QC           } from '../modules/local/dragen_cnv_qc/main'
include { VAF_HISTOGRAM           } from '../modules/local/vaf_histogram/main'
include { EXOMEQC_SUMMARY         } from '../modules/local/exomeqc_summary/main'
include { FINGERPRINTCHECK        } from '../subworkflows/local/fingerprintcheck/main'
include { samplesheetToList       } from 'plugin/nf-schema'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_exomeqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXOMEQC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // initialize empty channels for collecting software versions and multiQC files
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def ch_summary_mqc_config = channel.empty()

    // Collect output files for multiQC
    // using the nf-schema samplesheetToList function for adding directories to multiqc via config file
    // --multiqc_inputs is optional (e.g. not set in -profile test)
    if (params.multiqc_inputs) {
        channel
            .fromList(samplesheetToList(params.multiqc_inputs, "${projectDir}/assets/multiqc_input_schema.json"))
            .set { ch_multiqc_inputs }
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_inputs)
    }

    if (params.enable_fastqscreen) {
        if (!params.fastqscreen_database) {
            error("Please provide --fastqscreen_database when --enable_fastqscreen is set")
        }

        // The database is a directory or a .tar.gz of that directory (the only form
        // that can be fetched from a URL, e.g. in -profile test). The unpacked
        // directory keeps the archive name, so relative DATABASE paths in
        // fastq_screen.conf ("<db_name>/<genome>/<index>") keep working.
        def fastqscreen_db = file(params.fastqscreen_database, checkIfExists: true)
        def ch_fastqscreen_db
        if (fastqscreen_db.name.endsWith('.tar.gz')) {
            UNTAR([ [ id: fastqscreen_db.name.replaceFirst(/\.tar\.gz$/, '') ], fastqscreen_db ])
            ch_fastqscreen_db = UNTAR.out.untar.map { _meta, db -> db }.first()
        } else {
            ch_fastqscreen_db = channel.value(fastqscreen_db)
        }

        FASTQSCREEN_FASTQSCREEN(
            ch_samplesheet,
            ch_fastqscreen_db
        )

        // MultiQC parses FastQ Screen from *_screen.txt output files.
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQSCREEN_FASTQSCREEN.out.txt.map { _meta, txt -> txt }
        )
    }

    if (params.enable_fingerprintcheck) {
        if (!params.fingerprint_haplotype_map) {
            error("Please provide --fingerprint_haplotype_map when --enable_fingerprintcheck is set")
        }
        // Picard ExtractFingerprint requires the reference the BAMs are aligned to
        if (!params.fasta) {
            error("Please provide --fasta (the reference the BAMs are aligned to) when --enable_fingerprintcheck is set")
        }
        // Rows without BAM, runid and FASTQs are validated per row in FINGERPRINTCHECK.

        FINGERPRINTCHECK(
            ch_samplesheet,
            params.fingerprint_haplotype_map,
            params.fasta,
            params.fasta_fai,
            params.fasta_dict,
            params.bwamem2_index,
            params.align_reads_for_fingerprint,
            params.fingerprint_database
        )

        // MultiQC expects plain file paths; fingerprintcheck outputs are tuple(meta, path).
        // Only the crosscheck among this run's samples goes to MultiQC's Picard module
        // (heatmap); the comparison with the database is reported as custom content.
        ch_multiqc_files       = ch_multiqc_files.mix(
            FINGERPRINTCHECK.out.crosscheck_metrics.map { _meta, crosscheck_metrics -> crosscheck_metrics }
        )
        ch_multiqc_files       = ch_multiqc_files.mix(
            FINGERPRINTCHECK.out.database_multiqc.map { _meta, database_report -> database_report }
        )
    }

    if (params.enable_cnvqc) {
        if (!params.dragen_dirs) {
            error("Please provide --dragen_dirs when --enable_cnvqc is set")
        }
        
        DRAGEN_CNV_QC(params.dragen_dirs)

        ch_multiqc_files       = ch_multiqc_files.mix(DRAGEN_CNV_QC.out.metrics)
        ch_multiqc_files       = ch_multiqc_files.mix(DRAGEN_CNV_QC.out.json)
        ch_multiqc_files       = ch_multiqc_files.mix(DRAGEN_CNV_QC.out.plots)
    }

    if (params.enable_vaf_histogram) {
        // the VCFs of the samplesheet column 'vcf', sorted by sample: [ [ samples ], [ vcfs ] ]
        ch_vaf_vcfs = ch_samplesheet
            .filter { meta, _reads -> meta.vcf }
            .map { meta, _reads -> [ meta.id, file(meta.vcf.toString(), checkIfExists: true) ] }
            .ifEmpty { error("--enable_vaf_histogram is set, but no sample in the samplesheet has a 'vcf'; add the column 'vcf' or set --enable_vaf_histogram false") }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] } ] }

        VAF_HISTOGRAM(ch_vaf_vcfs)

        ch_multiqc_files       = ch_multiqc_files.mix(VAF_HISTOGRAM.out.multiqc)
    }

    //
    // MODULE: compact summary (CNV QC, purity, ploidy, TMB, MSI, HRD) on top of the report
    //
    if (params.enable_summary) {
        if (!params.dragen_dirs) {
            log.warn("Summary skipped: --dragen_dirs is not set (set --enable_summary false to silence this)")
        } else {
            ch_summary_cnv = params.enable_cnvqc ?
                DRAGEN_CNV_QC.out.metrics.collect().ifEmpty([]) :
                channel.value([])
            ch_summary_fp  = params.enable_fingerprintcheck ?
                FINGERPRINTCHECK.out.crosscheck_metrics
                    .map { _meta, crosscheck_metrics -> crosscheck_metrics }
                    .collect()
                    .ifEmpty([]) :
                channel.value([])

            EXOMEQC_SUMMARY(
                channel.fromPath(params.dragen_dirs, type: 'dir', checkIfExists: true).collect(),
                ch_summary_cnv,
                ch_summary_fp
            )
            // passed as the last MultiQC config: sets intro_text (top of report)
            ch_summary_mqc_config = EXOMEQC_SUMMARY.out.multiqc_config
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'exomeqc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //

    // add DRAGEN output directory to the MultiQC inputs
    if (params.dragen_dirs) {
        ch_multiqc_files = ch_multiqc_files.mix(channel.fromPath(params.dragen_dirs))
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    // MultiQC configs in this order: pipeline default, extended (not replaced) by
    // --multiqc_config, then the summary config (sets intro_text, top of report)
    def ch_multiqc_configs = channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        .concat(
            multiqc_config ? channel.fromPath(multiqc_config, checkIfExists: true) : channel.empty(),
            ch_summary_mqc_config
        )
        .toList()
    MULTIQC(
        ch_multiqc_files.flatten().collect()
            .map { files -> [ [id: 'exomeqc'], files ] }
            .combine(ch_multiqc_configs.map { configs -> [ configs ] })
            .map { meta, files, configs ->
                [
                    meta,
                    files,
                    configs,
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

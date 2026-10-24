process FINGERPRINT_DATABASE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.4--py312h0fa9677_0' :
        'biocontainers/pyvcf3:1.0.4--py312h0fa9677_0' }"

    input:
    tuple val(meta), path(run_samples), path(fingerprints, stageAs: 'fingerprints/*'), path(haplotype_map), path(run_crosscheck), path(database), path(database_check), path(database_samples), path(database_crosscheck)

    output:
    tuple val(meta), path("*_matches.tsv")           , emit: matches
    tuple val(meta), path("*_fingerprints.tsv")      , emit: fingerprints
    tuple val(meta), path("*_summary.tsv")           , emit: summary
    tuple val(meta), path("*_overview_mqc.html")     , emit: overview
    tuple val(meta), path("*_matches_mqc.{tsv,html}"), emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyvcf3'), eval("python3 -c 'import vcf; print(vcf.VERSION)'"), topic: versions, emit: versions_pyvcf3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // without a given database there is no check result and no crosscheck with it;
    // with a database whose fingerprints were all replaced there is no crosscheck
    def check = database_check ? "--database-check ${database_check} --database-samples ${database_samples}" : ''
    def crosscheck = database_crosscheck ? "--crosscheck ${database_crosscheck}" : ''
    """
    # summarise the fingerprint check and classify this run's comparisons with other runs
    fingerprint_database_report.py \\
        --run-samples ${run_samples} \\
        --fingerprints ${fingerprints} \\
        --haplotype-map ${haplotype_map} \\
        --run-crosscheck ${run_crosscheck} \\
        --database ${database} \\
        ${check} \\
        ${crosscheck} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # create empty placeholders of all outputs
    touch ${prefix}_matches.tsv
    touch ${prefix}_fingerprints.tsv
    touch ${prefix}_summary.tsv
    touch ${prefix}_overview_mqc.html
    touch ${prefix}_matches_mqc.html
    """
}

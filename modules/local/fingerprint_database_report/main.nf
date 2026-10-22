process FINGERPRINT_DATABASE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(crosscheck_metrics), path(run_samples), path(database), path(previous_database, stageAs: 'previous/*')

    output:
    tuple val(meta), path("*_matches.tsv")                   , emit: matches
    tuple val(meta), path("*_summary.tsv")                   , emit: summary
    tuple val(meta), path("*_matches_mqc.{tsv,html}")        , emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the previous database is optional; without it the report says the database is new
    def previous = previous_database ? "--previous-database ${previous_database}" : ''
    """
    # classify this run's comparisons with other runs and write the report tables
    fingerprint_database_report.py \\
        --crosscheck ${crosscheck_metrics} \\
        --run-samples ${run_samples} \\
        --database ${database} \\
        ${previous} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # create empty placeholders of all outputs
    touch ${prefix}_matches.tsv
    touch ${prefix}_summary.tsv
    touch ${prefix}_matches_mqc.html
    """
}

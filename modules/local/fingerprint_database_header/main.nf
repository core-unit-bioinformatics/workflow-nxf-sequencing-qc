process FINGERPRINT_DATABASE_HEADER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.4--py312h0fa9677_0' :
        'biocontainers/pyvcf3:1.0.4--py312h0fa9677_0' }"

    input:
    tuple val(meta), path(database)
    path run_qc_headers, stageAs: 'run_qc/*'

    output:
    tuple val(meta), path("*_header.txt"), emit: header
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyvcf3'), eval("python3 -c 'import vcf; print(vcf.VERSION)'"), topic: versions, emit: versions_pyvcf3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # database header with one fingerprint QC line per fingerprint
    fingerprint_database_header.py \\
        --database ${database} \\
        --run-qc ${run_qc_headers} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # create an empty placeholder of the output
    touch ${prefix}_header.txt
    """
}

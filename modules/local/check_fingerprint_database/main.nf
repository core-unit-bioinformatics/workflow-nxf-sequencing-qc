process CHECK_FINGERPRINT_DATABASE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.4--py312h0fa9677_0' :
        'biocontainers/pyvcf3:1.0.4--py312h0fa9677_0' }"

    input:
    tuple val(meta), path(database)
    path haplotype_map
    path run_fingerprint_ids

    output:
    tuple val(meta), path("*_samples.tsv"), emit: samples
    tuple val(meta), path("*_keep.txt")   , emit: keep
    tuple val(meta), path("*_check.tsv")  , emit: check
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyvcf3'), eval("python3 -c 'import vcf; print(vcf.VERSION)'"), topic: versions, emit: versions_pyvcf3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # check the database against the haplotype map and list the fingerprints to keep
    check_fingerprint_database.py \\
        --database ${database} \\
        --haplotype-map ${haplotype_map} \\
        --run-fingerprint-ids ${run_fingerprint_ids} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # create empty placeholders of all outputs
    touch ${prefix}_samples.tsv
    touch ${prefix}_keep.txt
    touch ${prefix}_check.tsv
    """
}

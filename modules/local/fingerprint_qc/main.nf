process FINGERPRINT_QC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.4--py312h0fa9677_0' :
        'biocontainers/pyvcf3:1.0.4--py312h0fa9677_0' }"

    input:
    tuple val(meta), path(fingerprint), path(pileup)
    path haplotype_map

    output:
    tuple val(meta), path("*_sites.tsv")  , emit: sites
    tuple val(meta), path("*_header.txt") , emit: header
    tuple val(meta), path("*_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyvcf3'), eval("python3 -c 'import vcf; print(vcf.VERSION)'"), topic: versions, emit: versions_pyvcf3

    when:
    task.ext.when == null || task.ext.when

    script:
    // thresholds, QC date and pipeline version come with ext.args
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # QC of each fingerprint site, contamination and QC status of the sample
    fingerprint_qc.py \\
        --fingerprint ${fingerprint} \\
        --pileup ${pileup} \\
        --haplotype-map ${haplotype_map} \\
        --fingerprint-id '${meta.fingerprint_id}' \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # create empty placeholders of all outputs
    touch ${prefix}_sites.tsv
    touch ${prefix}_header.txt
    touch ${prefix}_summary.tsv
    """
}

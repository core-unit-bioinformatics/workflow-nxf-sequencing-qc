process VAF_HISTOGRAM {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    // one VCF per sample, in the same order; staged in numbered folders so equal file names do not collide
    tuple val(samples), path(vcfs, stageAs: 'vcf*/*')

    output:
    path("${task.ext.prefix ?: 'vaf_histogram'}.tsv")        , emit: tsv
    path("${task.ext.prefix ?: 'vaf_histogram'}_summary.tsv"), emit: summary
    path("${task.ext.prefix ?: 'vaf_histogram'}_mqc.json")   , emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'vaf_histogram'
    // pair every sample with its VCF: --input SAMPLE VCF
    def inputs = [ samples, vcfs instanceof List ? vcfs : [ vcfs ] ].transpose().collect { sample, vcf -> "--input '${sample}' '${vcf}'" }.join(' ')
    """
    # count the variants per VAF bin of every sample's VCF; \$args carries the
    # VAF source, bin width, sample regex and filter option from conf/modules.config
    vaf_histogram.py \\
        ${inputs} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: 'vaf_histogram'
    """
    # empty outputs with the real names
    touch ${prefix}.tsv
    touch ${prefix}_summary.tsv
    touch ${prefix}_mqc.json
    """
}

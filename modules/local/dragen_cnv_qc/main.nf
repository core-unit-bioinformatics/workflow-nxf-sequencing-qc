process DRAGEN_CNV_QC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-karyoploter:1.32.0--r44hdfd78af_0':
        'quay.io/biocontainers/bioconductor-karyoploter:1.32.0--r44hdfd78af_0' }"

    input:
    path(dragen_dirs)

    output:
    path("*.csv") , emit: metrics
    // the JSON needs *.cnv.purity.coverage.models.tsv, the plots *.tn.tsv.gz or
    // *.tumor.baf.bedgraph.gz; both inputs are optional in the DRAGEN output
    path("*.json"), emit: json , optional: true
    path("*.png") , emit: plots, optional: true
    path "versions.yml"
    tuple val("${task.process}"), val('karyoploter'), eval(''), topic: versions, emit: versions_karyoploter

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    
    """
    # compute the CNV QC metrics and PASS/FAIL status; \$args carries the
    # QC thresholds from conf/modules.config (omitted thresholds are not checked)
    analyze_cnv_metrics.R \\
        --dragen $dragen_dirs \\
        --prefix "${prefix}" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | grep "R version")
        karyoploteR: \$(Rscript -e "cat(as.character(utils::packageVersion('karyoploteR')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    """
    echo $args
    
    # same file names as the script writes with the default ext.args
    # (per-sample purity/coverage PNGs are off: --plot-purcov-maxmodel=false)
    touch ${prefix}cnv-metrics.csv
    touch ${prefix}purcov_maxmodel_lineplot.json
    touch ${prefix}karyoplot-stub-tn_mqc.png
    touch ${prefix}karyoplot-stub-tumor.baf_mqc.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | grep "R version")
        karyoploteR: \$(Rscript -e "cat(as.character(utils::packageVersion('karyoploteR')))")
    END_VERSIONS
    """
}

// SPDX-License-Identifier: MIT — see LICENSE in the repository root (CUBI, HHU)
process SEQUENCINGQC_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path(dragen_dirs, stageAs: 'dragen??/*')
    path(cnv_metrics, stageAs: 'cnv/*')
    path(crosscheck, stageAs: 'fingerprint/*')

    output:
    path("*sequencingqc_summary_multiqc_config.yml"), emit: multiqc_config
    path("*sequencingqc_summary.html")              , emit: html
    path("*sequencingqc_summary.tsv")               , emit: tsv
    path("*biomarker.tsv")                     , emit: biomarker
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('sequencingqc_summary.py'), eval('sequencingqc_summary.py --version'), topic: versions, emit: versions_script

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    // a single staged path is a Path (iterable over its name parts), not a List
    def cnv    = cnv_metrics ? "--cnv-metrics ${[cnv_metrics].flatten().join(' ')}" : ''
    def fp     = crosscheck  ? "--crosscheck ${[crosscheck].flatten().join(' ')}"   : ''
    """
    sequencingqc_summary.py \\
        --dragen ${[dragen_dirs].flatten().join(' ')} \\
        ${cnv} \\
        ${fp} \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: ''
    """
    echo '{}' > ${prefix}sequencingqc_summary_multiqc_config.yml
    touch ${prefix}sequencingqc_summary.html
    touch ${prefix}sequencingqc_summary.tsv
    touch ${prefix}biomarker.tsv
    """
}

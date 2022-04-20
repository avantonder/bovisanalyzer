process TBPROFILER_COLLATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::tb-profiler=3.0.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:3.0.8--pypyh5e36f6f_0' :
        'quay.io/biocontainers/tb-profiler:3.0.8--pypyh5e36f6f_0' }"

    input:
    path ('tbprofiler/*')

    output:
    tuple val(meta), path("tbprofiler.variants.txt")     , emit: variants
    tuple val(meta), path("tbprofiler.txt")              , emit: summary
    tuple val(meta), path("tbprofiler.lineage.itol.txt") , emit: lineage
    tuple val(meta), path("tbprofiler.json")             , emit: json
    tuple val(meta), path("tbprofiler.dr.itol.txt")      , emit: itol_dr
    tuple val(meta), path("tbprofiler.dr.indiv.itol.txt"), emit: itol_indiv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    tb-profiler \\
        collate --dir .\\
        $args
    """
}

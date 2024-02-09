process TBPROFILER_COLLATE {
    label 'process_low'

    conda "bioconda::tb-profiler=5.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:5.0.1--pyhdfd78af_1' :
        'quay.io/biocontainers/tb-profiler:5.0.1--pyhdfd78af_1' }"

    input:

    path json
    //path ('tbprofiler/*')
    //tuple val(meta), path(csv), path(json), path(txt)

    output:
    path("tbprofiler.variants.txt")     , emit: variants
    path("tbprofiler.txt")              , emit: summary
    path("tbprofiler.lineage.itol.txt") , emit: lineage
    path("tbprofiler.json")             , emit: json
    path("tbprofiler.dr.itol.txt")      , emit: itol_dr
    path("tbprofiler.dr.indiv.itol.txt"), emit: itol_indiv

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

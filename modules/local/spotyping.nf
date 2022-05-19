process SPOTYPING {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::spotyping=2.1-3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spotyping:2.1--3' :
        'quay.io/biocontainers/spotyping:2.1--3' }"

    input:

    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.log"), emit: xls
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "-r1 $reads" : "-r1 ${reads[0]} -r2 ${reads[1]}"
    def spotyping_version = '2.1'
    """
    SpoTyping.py \\
        $args \\
        $input_reads \\
        -o ${prefix}.txt
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SpoTyping.py: ${spotyping_version}
    END_VERSIONS      
    """
}

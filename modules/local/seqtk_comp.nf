process SEQTK_COMP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")    , emit: tsv
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk \\
        comp \\
        $args \\
        $fasta \\
        > ${prefix}.tsv 
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

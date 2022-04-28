process MASHTREE {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mashtree=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mashtree:1.2.0--pl526h516909a_0' :
        'quay.io/biocontainers/mashtree:1.2.0--pl526h516909a_0' }"

    input:
    path seqs

    output:
    path "mashtree.dnd"       , emit: tree
    path "mashtree_matrix.tsv", emit: matrix
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mashtree \\
        $args \\
        --numcpus $task.cpus \\
        --outtree mashtree.dnd \\
        $seqs
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
process REMOVE_BLOCKS {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=2.7.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:2.7.9' :
        'quay.io/biocontainers/python:2.7.9' }"

    input:
    path alignment
    path tab

    output:
    path "*.fas"       , emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def args = task.ext.args   ?: ''
    def mask_version = '0.1'
    script: 
    """
    remove_blocks_from_aln.py \\
        $args \\
        -a $alignment \\
        -t $tab \\
        -o masked_alignment.fas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        remove_blocks_from_aln.py: ${mask_version}
    END_VERSIONS
    """
}
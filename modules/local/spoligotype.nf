process SPOLIGOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vsnp=2.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsnp:2.03--hdfd78af_2' :
        'quay.io/biocontainers/vsnp:2.03--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "-r1 $reads" : "-r1 ${reads[0]} -r2 ${reads[1]}"
    def spoligotype_version = '2.03'
    script: 
    """
    vsnp_spoligotype.py \\
        $args \\
        $input_reads
        
    head -n1 spoligo.txt > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsnp_spoligotype.py: ${spoligotype_version}
    END_VERSIONS
    """
}

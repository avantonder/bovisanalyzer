process VCF2PSEUDOGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.14' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.14--h88f3f91_0' :
        'quay.io/biocontainers/bcftools:1.14--h88f3f91_0' }"

    input:
    tuple val(meta), path(vcf), path(bed)
    path reference

    output:
    tuple val(meta), path("*.fas"), emit: pseudogenome
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools index \\
        $args \\
        $vcf
    
    bcftools consensus \\
        $args \\
        -f $reference \\
        -e 'TYPE="indel"' \\
        -m $bed \\
        $vcf \\
        | sed "/^>/ s/.*/>${prefix}/" > ${prefix}.fas
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
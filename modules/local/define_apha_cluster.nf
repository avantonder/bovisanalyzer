process DEFINE_APHA_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::numpy=1.15.2 conda-forge::pandas=0.23.4 conda-forge::scipy=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' :
        'quay.io/biocontainers/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' }"

    input:
    tuple val(meta), path(csv)
    tuple val(meta), path(vcf)
    path patternsDetailsFile
    path patternsBritishBTBFile
    path patternsPinnipediiFile
    path patternsMic_PinFile
    path patternsMicrotiFile
    path patternsBTBFile

    output:
    path "*.csv",        emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parser_version = '1.0'
    """
    Stage1-test.py $readstats $vcf $prefix
    mv _stage1.csv ${prefix}_stage1.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Stage1-test.py: ${parser_version}
    END_VERSIONS 
    """
}
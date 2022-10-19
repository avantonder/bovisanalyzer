process READ_STATS {
    label 'process_low'
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }
    
    conda (params.enable_conda ? 'conda-forge::numpy=1.15.2 conda-forge::pandas=0.23.4 conda-forge::scipy=1.2.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0"
    }
    
    input:
    tuple val(meta), path(raw_json)
    tuple val(meta), path(trim_json)
    tuple val(meta), path(depth)
    tuple val(meta), path(mapreads)

    output:
    path "*.csv",          emit: csv
    path  "versions.yml",  emit: versions
    
    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parser_version = '1.0'
    """
    read_stats.py
    mv read_stats.csv ${prefix}.read_stats.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_stats.py: ${parser_version}
    END_VERSIONS 
    """
}
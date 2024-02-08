process SPOLIGOPARSE {
    label 'process_low'
    
    conda (params.enable_conda ? "conda-forge::numpy=1.15.2 conda-forge::pandas=0.23.4 conda-forge::scipy=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' :
        'quay.io/biocontainers/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' }"
    
    input:
    path txt
    path db

    output:
    path "spoligotype_summary.tsv", emit: tsv
    path  "versions.yml"          , emit: versions
    
    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def parser_version = '2.0'
    """
    spoligo_parser.py
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spoligo_parser.py: ${parser_version}
    END_VERSIONS 
    """
}
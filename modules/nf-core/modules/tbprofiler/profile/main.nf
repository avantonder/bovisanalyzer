process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_high'
    label 'error_retry'

    conda (params.enable_conda ? "bioconda::tb-profiler=4.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:4.2.0--pypyh5e36f6f_0' :
        'quay.io/biocontainers/tb-profiler:tb-profiler:4.2.0--pypyh5e36f6f_0' }"

    input:
    path tbdb_barcode
    path tbdb_bed
    path tbdb_drjson
    path tbdb_fasta
    path tbdb_gff
    path tbdb_varjson
    path tbdb_verjson
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("bam/*.bam")     , emit: bam
    tuple val(meta), path("results/*.csv") , emit: csv, optional: true
    tuple val(meta), path("results/*.json"), emit: json
    tuple val(meta), path("results/*.txt") , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")  , emit: vcf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    """
    tb-profiler \\
        profile \\
        $args \\
        --csv \\
        --txt \\
        --no_trim \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        --external_db ./tbdbnew \\
        $input_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tbprofiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}

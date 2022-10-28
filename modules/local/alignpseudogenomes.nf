process ALIGNPSEUDOGENOMES {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    path pseudogenomes
    path reference

    output:
    tuple env(NUM_ALIGNMENT_GENOMES), path("aligned_pseudogenomes.fas"), emit: aligned_pseudogenomes
    path "low_quality_pseudogenomes.tsv",                                emit: low_quality_metrics

    script: // This script is bundled with the pipeline, in nf-core/bactmap/bin/
    """
    touch low_quality_pseudogenomes.tsv
    touch aligned_pseudogenomes.fas
    for pseudogenome in ${pseudogenomes}
    do
        fraction_non_GATC_bases=\$(calculate_fraction_of_non_GATC_bases.py -f \$pseudogenome | tr -d '\\n')
        if awk 'BEGIN { exit !(\$fraction_non_GATC_bases < ${params.non_GATC_threshold}) }'; then
            cat \$pseudogenome >> aligned_pseudogenomes.fas
        else
            echo "\$pseudogenome\t\$fraction_non_GATC_bases" >> low_quality_pseudogenomes.tsv
        fi
    done
    reference_to_single_sequence.py -r ${reference} -o final_reference.fas
    cat final_reference.fas >> aligned_pseudogenomes.fas

    NUM_ALIGNMENT_GENOMES=\$(grep -c ">" aligned_pseudogenomes.fas)
    """
}

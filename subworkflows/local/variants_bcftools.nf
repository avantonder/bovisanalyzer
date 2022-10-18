/*
 * Variant calling and downstream processing for BCFTools
 */

include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/modules/bcftools/mpileup/main'
include { BCFTOOLS_FILTER  } from '../../modules/nf-core/modules/bcftools/filter/main'
include { BCFTOOLS_VIEW    } from '../../modules/nf-core/modules/bcftools/view/main'

workflow VARIANTS_BCFTOOLS {
    take:
    bam       // channel: [ val(meta), [ bam ] ]
    fasta     // channel: /path/to/genome.fasta
    tsv       // channel: /path/to/DiscrimPos.tsv
    
    main:
    /*
     * MODULE Call variants
     */
    BCFTOOLS_MPILEUP ( bam, fasta )
    /*
     * MODULE Filter variants
     */
    BCFTOOLS_FILTER ( BCFTOOLS_MPILEUP.out.vcf  )
    /*
     * MODULE Filter discriminatory variants
     */
    BCFTOOLS_VIEW ( BCFTOOLS_MPILEUP.out.vcf, tsv)

    emit:
    filtered_vcf     = BCFTOOLS_FILTER.out.vcf       // channel: [ val(meta), [ vcf ] ]
    discrim_vcf      = BCFTOOLS_VIEW.out.vcf         // channel: [ val(meta), [ vcf ] ]
    vcf              = BCFTOOLS_MPILEUP.out.vcf      // channel: [ val(meta), [ vcf ] ]
    tbi              = BCFTOOLS_MPILEUP.out.tbi      // channel: [ val(meta), [ tbi ] ]
    stats            = BCFTOOLS_MPILEUP.out.stats    // channel: [ val(meta), [ txt ] ]
    bcftools_version = BCFTOOLS_MPILEUP.out.versions // path: *.version.txt
}

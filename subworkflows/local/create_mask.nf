/*
 * Create mask bed file with bedtools
 */

include { BCFTOOLS_FILTER    } from '../../modules/nf-core/modules/bcftools/filter/main'
include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/modules/bedtools/genomecov/main'
include { BEDTOOLS_MERGE     } from '../../modules/nf-core/modules/bedtools/merge/main'

workflow CREATE_MASK {
    take:
    bam_vcf
    mask      // channel: /path/to/DataDrivenMerge20.bed
    
    main:
    
    ch_versions = Channel.empty()

    /*
     * MODULE Filter variants
     */
    bam_vcf
        .map {
            meta, bam, vcf -> [ meta, vcf ]
        }
        .set { ch_vcf }
    
    BCFTOOLS_FILTER ( ch_vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    /*
     * MODULE Calculate zero coverage
     */
    BEDTOOLS_GENOMECOV ( bam_vcf )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())
    /*
     * MODULE Extract low quality regions
     */
    BCFTOOLS_FILTER.out.vcf                  // tuple val(meta), path(vcf)
        .join ( BEDTOOLS_GENOMECOV.out.bed ) // tuple val(meta), path(bed)
        .set  { ch_vcf_bed }                 // tuple val(meta), path(vcf), path(bed)
    
    BEDTOOLS_MERGE ( 
        ch_vcf_bed, 
        mask
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    emit:
    filtered_vcf = BCFTOOLS_FILTER.out.vcf         // channel: [ val(meta), [ vcf ] ]
    zerocov_bed  = BEDTOOLS_GENOMECOV.out.bed      // channel: [ val(meta), [ bed ] ]
    quality_bed  = BEDTOOLS_MERGE.out.qual_bed     // channel: [ val(meta), [ bed ] ]
    mask_bed     = BEDTOOLS_MERGE.out.bed          // channel: [ val(meta), [ tbi ] ]
    
    versions     = ch_versions                     // channel: [ versions.yml ]
}
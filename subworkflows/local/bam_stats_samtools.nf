/*
 * Run SAMtools stats, flagstat and idxstats
 * From https://github.com/nf-core/viralrecon/blob/dev/subworkflows/nf-core/bam_stats_samtools.nf
 */

include { SAMTOOLS_STATS    } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/modules/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [bai] ]
    
    main:
    SAMTOOLS_STATS    ( bam_bai )
    SAMTOOLS_FLAGSTAT ( bam_bai )
    SAMTOOLS_IDXSTATS ( bam_bai )

    emit:
    stats    = SAMTOOLS_STATS.out.stats            // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat      // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats      // channel: [ val(meta), [ idxstats ] ]
    samtools_version  = SAMTOOLS_STATS.out.versions //    path: *.version.txt
}
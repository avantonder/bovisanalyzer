//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/modules/samtools/flagstat/main'
include { SAMTOOLS_DEPTH    } from '../../modules/nf-core/modules/samtools/depth/main'
include { SAMTOOLS_VIEW     } from '../../modules/nf-core/modules/samtools/view/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai/csi] ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_STATS (
        ch_bam_bai,
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_FLAGSTAT (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    SAMTOOLS_DEPTH (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    SAMTOOLS_VIEW (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    depth    = SAMTOOLS_DEPTH.out.tsv         // channel: [ val(meta), [ depth ] ]
    mapreads = SAMTOOLS_VIEW.out.tsv         // channel: [ val(meta), [ mapreads ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
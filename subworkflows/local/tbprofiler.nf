//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { TBPROFILER_PROFILE } from '../modules/nf-core/modules/tbprofiler/profile/main'

workflow TBPROFILER {
    take:
    reads             // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    TBPROFILER_PROFILE (
        reads
    )
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE.out.versions.first())

    TBPROFILER_PROFILE
        .out
        .csv
        .join(TBPROFILER_PROFILE.out.json)
        .join(TBPROFILER_PROFILE.out.txt)
        .map { meta, csv, json, txt -> [ meta, csv, json, txt ] }
        .set { ch_tbprofiler }

    emit:
    csv        = TBPROFILER_PROFILE.out.csv          // channel: [ val(meta), [ csv  ] ]
    json       = TBPROFILER_PROFILE.out.json         // channel: [ val(meta), [ json ] ]
    txt        = TBPROFILER_PROFILE.out.txt          // channel: [ val(meta), [ txt  ] ]

    versions   = ch_versions                         // channel: [ versions.yml ]
}
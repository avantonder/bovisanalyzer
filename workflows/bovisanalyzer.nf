/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBovisanalyzer.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.reference,  params.kraken2db, params.brackendb]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.reference) { ch_reference = file(params.reference) } else { exit 1, 'Reference fasta file not specified!' }
if (params.kraken2db) { ch_kraken2db = file(params.kraken2db) } else { exit 1, 'kraken2 database not specified!' }
if (params.brackendb) { ch_brackendb = file(params.brackendb) } else { exit 1, 'bracken database not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { KRAKENPARSE                 } from '../modules/local/krakenparse'
include { TBPROFILER_COLLATE          } from '../modules/local/tbprofiler_collate'
include { VCF2PSEUDOGENOME            } from '../modules/local/vcf2pseudogenome'
include { ALIGNPSEUDOGENOMES          } from '../modules/local/alignpseudogenomes'

include { INPUT_CHECK                 } from '../subworkflows/local/input_check'
include { FASTQC_FASTP                } from '../subworkflows/local/fastqc_fastp'
include { BAM_SORT_SAMTOOLS           } from '../subworkflows/local/bam_sort_samtools'
include { VARIANTS_BCFTOOLS           } from '../subworkflows/local/variants_bcftools'
include { SUB_SAMPLING                } from '../subworkflows/local/sub_sampling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQSCAN                                               } from '../modules/nf-core/modules/fastqscan/main'
include { KRAKEN2_KRAKEN2                                         } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { BRACKEN_BRACKEN                                         } from '../modules/nf-core/modules/bracken/bracken/main'
include { BWA_INDEX                                               } from '../modules/nf-core/modules/bwa/index/main'
include { TBPROFILER_PROFILE                                      } from '../modules/nf-core/modules/tbprofiler/profile/main'
include { BWA_MEM                                                 } from '../modules/nf-core/modules/bwa/mem/main'
include { MULTIQC                                                 } from '../modules/nf-core/modules/multiqc/main'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_READS         } from '../modules/local/multiqc_tsv_from_list'
include { CUSTOM_DUMPSOFTWAREVERSIONS                             } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []
def fail_mapped_reads = [:]

workflow BOVISANALYZER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run bwa index
    //
    BWA_INDEX (
        ch_reference
    )

    //
    // MODULE: Run fastq-scan
    //
    FASTQSCAN (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQSCAN.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_FASTP (
        INPUT_CHECK.out.reads,
        params.save_trimmed_fail,
        true
    )
    ch_variants_fastq = FASTQC_FASTP.out.reads
    ch_versions = ch_versions.mix(FASTQC_FASTP.out.versions)

    //
    // Filter empty FastQ files after adapter trimming
    //
    ch_fail_reads_multiqc = Channel.empty()
    if (!params.skip_fastp) {
        ch_variants_fastq
            .join(FASTQC_FASTP.out.trim_json)
            .map {
                meta, reads, json ->
                    pass = WorkflowBovisanalyzer.getFastpReadsAfterFiltering(json) > 0
                    [ meta, reads, json, pass ]
            }
            .set { ch_pass_fail_reads }

        ch_pass_fail_reads
            .map { meta, reads, json, pass -> if (pass) [ meta, reads ] }
            .set { ch_variants_fastq }

        ch_pass_fail_reads
            .map {
                meta, reads, json, pass ->
                if (!pass) {
                    fail_mapped_reads[meta.id] = 0
                    num_reads = WorkflowBovisanalyzer.getFastpReadsBeforeFiltering(json)
                    return [ "$meta.id\t$num_reads" ]
                }
            }
            .set { ch_pass_fail_reads }

        MULTIQC_TSV_FAIL_READS (
            ch_pass_fail_reads.collect(),
            ['Sample', 'Reads before trimming'],
            'fail_mapped_reads'
        )
        .set { ch_fail_reads_multiqc }
    }

    //
    // MODULE: Run kraken2
    //  
    ch_kraken2_multiqc = Channel.empty()
    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2 (
                ch_variants_fastq,
                ch_kraken2db
            )
        ch_kraken2_bracken       = KRAKEN2_KRAKEN2.out.txt
        ch_kraken2_krakenparse   = KRAKEN2_KRAKEN2.out.txt
        ch_kraken2_multiqc       = KRAKEN2_KRAKEN2.out.txt
        ch_versions              = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first().ifEmpty(null))
    }
    
    //
    // MODULE: Run bracken
    //
    BRACKEN_BRACKEN (
            ch_kraken2_bracken,
            ch_brackendb
        )
    ch_bracken_krakenparse = BRACKEN_BRACKEN.out.reports
    ch_versions            = ch_versions.mix(BRACKEN_BRACKEN.out.versions.first())

    //
    // MODULE: Run krakenparse
    //
    KRAKENPARSE (
            ch_kraken2_krakenparse.collect{it[1]}.ifEmpty([]),
            ch_bracken_krakenparse.collect{it[1]}.ifEmpty([])
        )
    ch_versions = ch_versions.mix(KRAKENPARSE.out.versions.first())

    //
    // SUBWORKFLOW: Subsample reads
    //
    SUB_SAMPLING(
            ch_variants_fastq
        )
    ch_variants_fastq = SUB_SAMPLING.out.reads

    //
    // SUBWORKFLOW: TBprofiler
    //
    TBPROFILER_PROFILE(
            ch_variants_fastq
        )
    ch_tbprofiler_collate = TBPROFILER_PROFILE.out.json
    ch_versions           = ch_versions.mix(TBPROFILER_PROFILE.out.versions.first())

    //
    // MODULE: Collate TB-profiler outputs
    //
    TBPROFILER_COLLATE(
            ch_tbprofiler_collate.collect{it[1]}.ifEmpty([])
        )
    
    //
    // MODULE: Map reads
    //
    BWA_MEM (
        ch_variants_fastq,
        BWA_INDEX.out.index
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // SUBWORKFLOW: Sort bam files
    //
    BAM_SORT_SAMTOOLS (
        BWA_MEM.out.bam
    )
    ch_flagstat_multiqc = BAM_SORT_SAMTOOLS.out.flagstat
    ch_versions         = ch_versions.mix(BAM_SORT_SAMTOOLS.out.samtools_version.first())

    //
    // SUBWORKFLOW: Call variants
    //
    VARIANTS_BCFTOOLS (
        BAM_SORT_SAMTOOLS.out.bam,
        ch_reference
    )
    ch_bcftools_stats_multiqc = VARIANTS_BCFTOOLS.out.stats
    ch_versions               = ch_versions.mix(VARIANTS_BCFTOOLS.out.bcftools_version.first())

    //
    // MODULE: Make pseudogenome from VCF
    //
    VCF2PSEUDOGENOME (
        VARIANTS_BCFTOOLS.out.filtered_vcf,
        ch_reference
    )

    //
    // MODULE: make pseudogenome alignment
    //
    ALIGNPSEUDOGENOMES (
        VCF2PSEUDOGENOME.out.pseudogenome.map { pseudogenome -> pseudogenome[1] }.collect(),
        ch_reference
    )
    ALIGNPSEUDOGENOMES.out.aligned_pseudogenomes
        .branch {
            aligned_pseudogenomes ->
            ALIGNMENT_NUM_PASS: aligned_pseudogenomes[0].toInteger() >= 4
            ALIGNMENT_NUM_FAIL: aligned_pseudogenomes[0].toInteger() < 4
        }
        .set { aligned_pseudogenomes_branch }

    // Don't proceeed further if two few genonmes
    aligned_pseudogenomes_branch.ALIGNMENT_NUM_FAIL.view { "Insufficient (${it[0]}) genomes after filtering to continue. Check results/pseudogenomes/low_quality_pseudogenomes.tsv for details"}

    aligned_pseudogenomes_branch.ALIGNMENT_NUM_PASS
        .map{ it[1] }
        .set { aligned_pseudogenomes }

    //
    // MODULE: Collate software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBovisanalyzer.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        ch_fail_reads_multiqc.ifEmpty([]),
        FASTQC_FASTP.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
        FASTQC_FASTP.out.trim_json.collect{it[1]}.ifEmpty([]),
        ch_kraken2_multiqc.collect{it[1]}.ifEmpty([]),
        ch_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
        ch_bcftools_stats_multiqc.collect{it[1]}.ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, fail_mapped_reads)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

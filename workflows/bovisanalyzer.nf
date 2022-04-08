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
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.reference) { ch_reference = file(params.reference) } else { exit 1, 'Reference fasta file not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { KRAKENPARSE                 } from '../modules/local/krakenparse'
include { VCF2PSEUDOGENOME            } from '../modules/local/vcf2pseudogenome'
include { ALIGNPSEUDOGENOMES          } from '../modules/local/alignpseudogenomes'

include { INPUT_CHECK                 } from '../subworkflows/local/input_check'
include { FASTQC_FASTP                } from '../subworkflows/local/fastqc_fastp'            addParams( fastqc_raw_options: modules['fastqc_raw'], fastqc_trim_options: modules['fastqc_trim'], fastp_options: fastp_options )
include { BAM_SORT_SAMTOOLS           } from '../subworkflows/local/bam_sort_samtools' addParams( samtools_sort_options: modules['samtools_sort'], samtools_index_options : modules['samtools_index'], bam_stats_options: modules['bam_stats'])
include { VARIANTS_BCFTOOLS           } from '../subworkflows/local/variants_bcftools' addParams( bcftools_mpileup_options: modules['bcftools_mpileup'], bcftools_filter_options: modules['bcftools_filter'])
include { SUB_SAMPLING                } from '../subworkflows/local/sub_sampling'      addParams( mash_sketch_options: modules['mash_sketch'], rasusa_options: modules['rasusa'])
include { TBPROFILER                  } from '../subworkflows/local/tbprofiler/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQSCAN                   } from '../modules/nf-core/modules/fastqscan/main'
include { KRAKEN2_KRAKEN2             } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { BRACKEN_BRACKEN             } from '../modules/nf-core/modules/bracken/bracken/main'
include { BWA_INDEX                   } from '../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM                     } from '../modules/nf-core/modules/bwa/mem/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_to_base: true] )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

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
        INPUT_CHECK.out.reads
    )
    ch_reads    = FASTQC_FASTP.out.reads
    ch_versions = ch_versions.mix(FASTQC_FASTP.out.fastqc_version.first())
    ch_versions = ch_versions.mix(FASTQC_FASTP.out.fastp_version.first())

    //
    // MODULE: Run kraken2
    //  
    KRAKEN2_KRAKEN2 (
            ch_reads,
            ch_kraken2db
        )
    ch_kraken2_bracken       = KRAKEN2_KRAKEN2.out.txt
    ch_kraken2_krakenparse   = KRAKEN2_KRAKEN2.out.txt
    ch_versions     = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    //
    // MODULE: Run bracken
    //
    BRACKEN_BRACKEN (
            ch_kraken2_bracken,
            ch_brackendb
        )
    ch_bracken_krakenparse = BRACKEN_BRACKEN.out.reports
    ch_versions   = ch_versions.mix(BRACKEN_BRACKEN.out.versions.first())

    //
    // MODULE: Run krakenparse
    //
    KRAKENPARSE (
            ch_kraken2_krakenparse.collect{it[1]}.ifEmpty([]),
            ch_bracken_krakenparse.collect{it[1]}.ifEmpty([])
        )
    ch_versions = ch_versions.mix(KRAKENPARSE.out.version.first())

    //
    // SUBWORKFLOW: Subsample reads
    //
    SUB_SAMPLING(
            ch_reads
        )
    ch_reads = SUB_SAMPLING.out.reads

    //
    // SUBWORKFLOW: TBprofiler
    //
    TBPROFILER(
            ch_reads
        )
    ch_versions = ch_versions.mix(TBPROFILER.out.versions.first())
    
    //
    // MODULE: Map reads
    //
    BWA_MEM (
        ch_reads,
        BWA_INDEX.out.index
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // SUBWORKFLOW: Sort bam files
    //
    BAM_SORT_SAMTOOLS (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.samtools_version.first())

    //
    // SUBWORKFLOW: Call variants
    //
    VARIANTS_BCFTOOLS (
        BAM_SORT_SAMTOOLS.out.bam,
        ch_reference
    )
    ch_versions = ch_versions.mix(VARIANTS_BCFTOOLS.out.bcftools_version.first())

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

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
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
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

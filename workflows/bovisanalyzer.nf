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
//if (params.kraken2db) { ch_kraken2db = file(params.kraken2db) } else { exit 1, 'kraken2 database not specified!' }
//if (params.brackendb) { ch_brackendb = file(params.brackendb) } else { exit 1, 'bracken database not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_tbdb_barcode           = file("$projectDir/assets/tbdb/tbdbnew.barcode.bed",                           checkIfExists: true)
ch_tbdb_bed               = file("$projectDir/assets/tbdb/tbdbnew.bed",                                   checkIfExists: true)
ch_tbdb_drjson            = file("$projectDir/assets/tbdb/tbdbnew.dr.json",                               checkIfExists: true)
ch_tbdb_fasta             = file("$projectDir/assets/tbdb/tbdbnew.fasta",                                 checkIfExists: true)
ch_tbdb_gff               = file("$projectDir/assets/tbdb/tbdbnew.gff",                                   checkIfExists: true)
ch_tbdb_varjson           = file("$projectDir/assets/tbdb/tbdbnew.variables.json",                        checkIfExists: true)
ch_tbdb_verjson           = file("$projectDir/assets/tbdb/tbdbnew.version.json",                          checkIfExists: true)
ch_spoligotype_db         = file("$projectDir/assets/spoligotype_db.tsv",                                 checkIfExists: true)
ch_discrimpos             = file("$projectDir/assets/DiscrimPos.tsv",                                     checkIfExists: true)
ch_patternsDetailsFile    = file("$projectDir/assets/Stage1_patterns/CSSnewclusters_LT708304_230119.csv", checkIfExists: true)
ch_patternsBritishBTBFile = file("$projectDir/assets/Stage1_patterns/patternsBritishBTB_LT708304.csv",    checkIfExists: true)
ch_patternsPinnipediiFile = file("$projectDir/assets/Stage1_patterns/patternsPinnipedii_LT708304.csv",    checkIfExists: true)
ch_patternsMic_PinFile    = file("$projectDir/assets/Stage1_patterns/patternsMic_Pin_LT708304.csv",       checkIfExists: true)
ch_patternsMicrotiFile    = file("$projectDir/assets/Stage1_patterns/patternsMicroti_LT708304.csv",       checkIfExists: true)
ch_patternsBTBFile        = file("$projectDir/assets/Stage1_patterns/patternsBTB_LT708304.csv",           checkIfExists: true)
ch_mask                   = file("$projectDir/assets/DataDrivenMerge20.bed",                              checkIfExists: true)
ch_multiqc_config         = file("$projectDir/assets/multiqc_config.yml",                                 checkIfExists: true)
ch_multiqc_custom_config  = params.multiqc_config ? file(params.multiqc_config) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQSCANPARSE as FASTQSCANPARSE_RAW  } from '../modules/local/fastqscanparse'
include { FASTQSCANPARSE as FASTQSCANPARSE_TRIM } from '../modules/local/fastqscanparse'
include { KRAKENPARSE                           } from '../modules/local/krakenparse'
include { TBPROFILER_COLLATE                    } from '../modules/local/tbprofiler_collate'
include { SPOTYPING                             } from '../modules/local/spotyping'
include { SPOLIGOPARSE                          } from '../modules/local/spoligoparse'
include { READ_STATS                            } from '../modules/local/read_stats'
include { READSTATS_PARSE                       } from '../modules/local/readstats_parse'
include { DEFINE_APHA_CLUSTER                   } from '../modules/local/define_apha_cluster'
include { CLUSTER_PARSE                         } from '../modules/local/cluster_parse'
include { VCF2PSEUDOGENOME                      } from '../modules/local/vcf2pseudogenome'
include { SEQTK_COMP                            } from '../modules/local/seqtk_comp'
include { SEQTK_PARSE                           } from '../modules/local/seqtk_parse'
include { METADATA_COLLATE                      } from '../modules/local/metadata_collate'
include { ALIGNPSEUDOGENOMES                    } from '../modules/local/alignpseudogenomes'

include { INPUT_CHECK                           } from '../subworkflows/local/input_check'
include { FASTQC_FASTP                          } from '../subworkflows/local/fastqc_fastp'
include { BAM_SORT_SAMTOOLS                     } from '../subworkflows/local/bam_sort_samtools'
include { BAM_MARKDUPLICATES_PICARD             } from '../subworkflows/local/bam_markduplicates_picard'
include { VARIANTS_BCFTOOLS                     } from '../subworkflows/local/variants_bcftools'
include { SUB_SAMPLING                          } from '../subworkflows/local/sub_sampling'
include { CREATE_MASK                           } from '../subworkflows/local/create_mask'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQSCAN as FASTQSCAN_RAW                      } from '../modules/nf-core/modules/fastqscan/main'
include { FASTQSCAN as FASTQSCAN_TRIM                     } from '../modules/nf-core/modules/fastqscan/main'
include { KRAKEN2_KRAKEN2                                 } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { BRACKEN_BRACKEN                                 } from '../modules/nf-core/modules/bracken/bracken/main'
include { BWA_INDEX                                       } from '../modules/nf-core/modules/bwa/index/main'
include { TBPROFILER_PROFILE                              } from '../modules/nf-core/modules/tbprofiler/profile/main'
include { BWA_MEM                                         } from '../modules/nf-core/modules/bwa/mem/main'
include { PICARD_COLLECTMULTIPLEMETRICS                   } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { SNPSITES                                        } from '../modules/nf-core/modules/snpsites/main'
include { MULTIQC                                         } from '../modules/nf-core/modules/multiqc/main'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_READS } from '../modules/local/multiqc_tsv_from_list'
include { CUSTOM_DUMPSOFTWAREVERSIONS                     } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' 

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
    FASTQSCAN_RAW (
        INPUT_CHECK.out.reads
    )
    ch_fastqscanraw_fastqscanparse = FASTQSCAN_RAW.out.json
    ch_fastqscanraw_readstats      = FASTQSCAN_RAW.out.json
    ch_versions                    = ch_versions.mix(FASTQSCAN_RAW.out.versions.first())

    //
    // MODULE: Run fastqscanparse
    //
    FASTQSCANPARSE_RAW (
        ch_fastqscanraw_fastqscanparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQSCANPARSE_RAW.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_FASTP (
        INPUT_CHECK.out.reads,
        params.save_trimmed_fail,
        false
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
    // MODULE: Run fastq-scan
    //
    FASTQSCAN_TRIM (
        ch_variants_fastq
    )
    ch_fastqscantrim_fastqscanparse = FASTQSCAN_TRIM.out.json
    ch_versions                     = ch_versions.mix(FASTQSCAN_TRIM.out.versions.first())

    //
    // MODULE: Run fastqscanparse
    //
    FASTQSCANPARSE_TRIM (
            ch_fastqscantrim_fastqscanparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQSCANPARSE_TRIM.out.versions.first())
    
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
    ch_kraken_metadata = KRAKENPARSE.out.composition
    ch_versions = ch_versions.mix(KRAKENPARSE.out.versions.first())

    //
    // MODULE: Subsample reads
    //
    SUB_SAMPLING(
        ch_variants_fastq
    )
    ch_variants_fastq = SUB_SAMPLING.out.reads
    ch_versions = ch_versions.mix(SUB_SAMPLING.out.versions.first())

    //
    // MODULE: TBprofiler
    //
    //ch_tbprofiler = Channel.empty()
    TBPROFILER_PROFILE(
        ch_tbdb_barcode,
        ch_tbdb_bed,
        ch_tbdb_drjson,
        ch_tbdb_fasta,
        ch_tbdb_gff,
        ch_tbdb_varjson,
        ch_tbdb_verjson,
        ch_variants_fastq
    )
    ch_tbprofiler_collate = TBPROFILER_PROFILE.out.json
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE.out.versions.first())
    
    //
    // MODULE: Collate TB-profiler outputs
    //
    TBPROFILER_COLLATE(
        ch_tbprofiler_collate.collect{it[1]}.ifEmpty([])
    )
    ch_tbprofiler_metadata = TBPROFILER_COLLATE.out.summary

    //
    // MODULE: Run SpoTyping
    //
    SPOTYPING (
        ch_variants_fastq
    )
    ch_spotyping_spoligoparse = SPOTYPING.out.txt
    ch_versions = ch_versions.mix(SPOTYPING.out.versions.first())

    //
    // MODULE: Run spoligoparse
    //
    SPOLIGOPARSE (
        ch_spoligotype_db,
        ch_spotyping_spoligoparse.collect{it[1]}.ifEmpty([])
    )
    ch_spoligo_metadata = SPOLIGOPARSE.out.tsv
    ch_versions = ch_versions.mix(SPOLIGOPARSE.out.versions.first())
        
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
    ch_bam                = BAM_SORT_SAMTOOLS.out.bam
    ch_flagstat_multiqc   = BAM_SORT_SAMTOOLS.out.flagstat
    ch_versions           = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions.first())

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    BAM_MARKDUPLICATES_PICARD (
        ch_bam,
        ch_reference,
        BWA_INDEX.out.index
    )
    ch_bam_mask                        = BAM_MARKDUPLICATES_PICARD.out.bam
    ch_markduplicates_flagstat_multiqc = BAM_MARKDUPLICATES_PICARD.out.flagstat
    ch_markduplicates_metrics_multiqc  = BAM_MARKDUPLICATES_PICARD.out.metrics
    ch_versions                        = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

    //
    // MODULE: Picard metrics
    //
    PICARD_COLLECTMULTIPLEMETRICS (
        ch_bam,
        ch_reference,
        []
    )
    ch_collectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    ch_versions                       = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first().ifEmpty(null))
    
    //
    // MODULE: Calculate read stats
    //
    ch_fastqscanraw_readstats                           // tuple val(meta), path(json)
        .join( FASTQSCAN_TRIM.out.json )                // tuple val(meta), path(json) 
        .join( BAM_MARKDUPLICATES_PICARD.out.depth )    // tuple val(meta), path(depth)
        .join( BAM_MARKDUPLICATES_PICARD.out.mapreads ) // tuple val(meta), path(mapreads)
        .set { ch_readstats }                           // tuple val(meta), path(json), path(json), path(depth), path(mapreads)

    READ_STATS (
        ch_readstats
    )
    ch_readstats_readstatsparse = READ_STATS.out.csv
    ch_readstats_cluster        = READ_STATS.out.csv
    ch_versions                 = ch_versions.mix(READ_STATS.out.versions.first())

    //
    // MODULE: Summarise read stats outputs
    //
    READSTATS_PARSE (
        ch_readstats_readstatsparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions           = ch_versions.mix(READSTATS_PARSE.out.versions.first())
    
    //
    // SUBWORKFLOW: Call variants
    //
    VARIANTS_BCFTOOLS (
        BAM_MARKDUPLICATES_PICARD.out.bam,
        ch_reference,
        ch_discrimpos
    )
    ch_pseudo_vcf             = VARIANTS_BCFTOOLS.out.filtered_vcf
    ch_bcftools_stats_multiqc = VARIANTS_BCFTOOLS.out.stats
    ch_versions               = ch_versions.mix(VARIANTS_BCFTOOLS.out.bcftools_version.first())

    //
    // MODULE: Define APHA cluster
    //
    ch_readstats_cluster                           // tuple val(meta), path(csv)
        .join( VARIANTS_BCFTOOLS.out.discrim_vcf ) // tuple val(meta), path(vcf)
        .set { ch_apha_cluster }                   // tuple val(meta), path(csv), path(vcf)

    DEFINE_APHA_CLUSTER (
        ch_apha_cluster,
        ch_patternsDetailsFile,
        ch_patternsBritishBTBFile,
        ch_patternsPinnipediiFile,
        ch_patternsMic_PinFile,
        ch_patternsMicrotiFile,
        ch_patternsBTBFile
    )
    ch_cluster_clusterparse = DEFINE_APHA_CLUSTER.out.csv
    ch_versions             = ch_versions.mix(DEFINE_APHA_CLUSTER.out.versions.first())
    
    //
    // MODULE: Summarise APHA cluster outputs
    //
    CLUSTER_PARSE (
        ch_cluster_clusterparse.collect{it[1]}.ifEmpty([])
    )
    ch_cluster_metadata = CLUSTER_PARSE.out.tsv
    ch_versions         = ch_versions.mix(CLUSTER_PARSE.out.versions.first())
    
    //
    // SUBWORKFLOW: Create mask bed file
    //
    ch_bam_mask                            // tuple val(meta), path(bam)
        .join( VARIANTS_BCFTOOLS.out.vcf ) // tuple val(meta), path(vcf)
        .set { ch_bam_vcf_mask }           // tuple val(meta), path(bam), path(vcf)
    
    CREATE_MASK (
        ch_bam_vcf_mask,
        ch_mask
    )
    ch_versions = ch_versions.mix(CREATE_MASK.out.versions.first())

    //
    // MODULE: Make pseudogenome from VCF
    //
    ch_pseudo_vcf                         // tuple val(meta), path(vcf)
        .join( CREATE_MASK.out.mask_bed ) // tuple val(meta), path(bed)
        .set { ch_pseudogenome }          // tuple val(meta), path(vcf), path(bed)
    
    VCF2PSEUDOGENOME (
        ch_pseudogenome,
        ch_reference
    )
    
    //
    // MODULE: Calculate number of mapped positions in pseudogenome
    //
    SEQTK_COMP (
        VCF2PSEUDOGENOME.out.pseudogenome
    )
    ch_seqtk_seqtkparse = SEQTK_COMP.out.tsv
    ch_versions = ch_versions.mix(SEQTK_COMP.out.versions.first())

    //
    // MODULE: Summarise seqtk outputs
    //
    SEQTK_PARSE (
        ch_seqtk_seqtkparse.collect{it[1]}.ifEmpty([])
    )
    ch_seqtk_metadata = SEQTK_PARSE.out.tsv
    ch_versions       = ch_versions.mix(SEQTK_PARSE.out.versions.first())

    //
    // MODULE: Collate all metadata
    //
    METADATA_COLLATE (
        ch_kraken_metadata,
        ch_tbprofiler_metadata,
        ch_spoligo_metadata,
        ch_cluster_metadata,
        ch_seqtk_metadata
    )
    ch_versions = ch_versions.mix(METADATA_COLLATE.out.versions.first())
    
    //
    // MODULE: Make pseudogenome alignment
    //
    if (!params.skip_alignment) {
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
        // MODULE: Extract SNPs from masked alignment
        //
        SNPSITES (
            aligned_pseudogenomes
        )
        ch_versions = ch_versions.mix(SNPSITES.out.versions.first())
    }
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
        ch_markduplicates_metrics_multiqc.collect{it[1]}.ifEmpty([]),
        ch_markduplicates_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
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

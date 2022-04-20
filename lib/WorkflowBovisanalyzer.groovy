//
// This file holds several functions specific to the workflow/bovisanalyzer.nf in the avantonder/bovisanalyzer pipeline
//

import groovy.json.JsonSlurper
import java.util.regex.Matcher

class WorkflowBovisanalyzer {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        if (!params.reference) {
            log.error "Reference fasta file not specified! e.g. '--reference genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.kraken2db) {
            log.error "kraken2 database not specified! e.g. '--kraken2db minikraken2_v1_8GB' or via a detectable config file."
            System.exit(1)
        }

        if (!params.brackendb) {
            log.error "bracken database not specified! e.g. '--brackendb minikraken2_v1_8GB/database100mers.kmer_distrib' or via a detectable config file."
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    // Function to extract genome size from mash stat file
    public static String find_genome_size(mash_output) {
        Matcher m = mash_output =~ /Estimated genome size: (.+)/
        String genome_size = Float.parseFloat(m[0][1]).toInteger().toString() + 'b'
        return genome_size
    }

    //
    // Function that parses fastp json output file to get total number of reads after trimming
    //
    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
    }

    //
    // Function that parses fastp json output file to get total number of reads before trimming
    //
    public static Integer getFastpReadsBeforeFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['before_filtering']['total_reads'].toInteger()
    }
}

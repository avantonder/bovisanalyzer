#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import sys
import glob
import json
import argparse

def parser_args(args=None):
    """ 
    Function for input arguments for read_stats.py
    """
    Description = 'Collect fastq-scan and samtools mapping/depth outputs and create a table for each sample'
    Epilog = """Example usage: python read_stats.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-rl", "--ref_length"  , type=int, default=4349904         , help="Length of reference file (default: 4349904).")
    parser.add_argument("-md", "--min_depth"   , type=int, default=10              , help="Minimum read depth per site  (default: 10).")
    parser.add_argument("-mp", "--min_pc"      , type=int, default=60              , help="Minimum percentage of reads mapping to reference (default: 60).")
    parser.add_argument("-mr", "--min_reads"   , type=int, default=600000          , help="Minimum number of required raw reads (default: 600000).")
    parser.add_argument("-mt", "--min_aft_trim", type=int, default=60              , help="Minimum percentage of reads after trimming  (default: 60).")
    parser.add_argument("-of", "--output_file" , type=str, default="read_stats.csv", help="read stats file (default: 'read_stats.csv').")
    return parser.parse_args(args)

def make_dir(path):
    """ 
    Function for making a directory from a provided path
    """
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

def json_to_dataframe(json_files):
    """ 
    Function to take list of json files and create a summary table
    """
    json_names = [i.replace('.json', '') for i in json_files]
    json_names_df = pd.DataFrame(json_names)
    json_names_df.columns = ['Sample']
    jsons_data = {}

    for index, file in enumerate(json_files):
        with open(file, 'r') as f:
            json_text = json.loads(f.read())
            qc = json_text['qc_stats']
            jsons_data[index] = qc

    jsons_data_df = pd.DataFrame.from_dict(jsons_data, orient = 'index')
    json_merged_df = json_names_df.join(jsons_data_df)
    json_merged_df = json_merged_df.iloc[:, [0,3]]

    return json_merged_df

def mapped_to_dataframe(mapped_files):
    """ 
    Function to take list of samtools mapped read files and create a summary table
    """
    mapped_read = [i.replace('.map.reads.tsv', '') for i in mapped_files]
    mapped_read_df = pd.DataFrame(mapped_read)
    mapped_read_df.columns = ['Sample']
    mapped_read_files_read = [open(file).read() for file in mapped_files]
    mapped_read_files_df = pd.DataFrame(mapped_read_files_read)
    mapped_read_files_df.columns = ['NumMappedReads']
    mapped_read_merged_df = mapped_read_df.join(mapped_read_files_df).replace('\n','', regex=True)

    return mapped_read_merged_df 

def depth_to_dataframe(depth_files,ref_length):
    """ 
    Function to take list of samtools depth files and create a summary table
    """
    depth_read = [i.replace('.depth.tsv', '') for i in depth_files]
    depth_read_df = pd.DataFrame(depth_read)
    depth_read_df.columns = ['Sample']
    depth_files_read = [pd.read_csv(f, sep='\t', header=None) for f in depth_files]
    
    # Calculate average mapped read depth
    average_read_depth = [i.iloc[:,[2]].mean() for i in depth_files_read]
    average_read_depth_df = pd.DataFrame(average_read_depth)
    average_read_depth_df.columns = ['MeanDepth']

    # Calculate number of bases with zero coverage
    zero_count = [(i.iloc[:,[2]] == 0).sum() for i in depth_files_read]
    zero_count_df = pd.DataFrame(zero_count)
    zero_count_df.columns = ['zeroCov']

    # Merge average and zero counts
    average_zero_df = average_read_depth_df.join(zero_count_df)

    # Calculate genome coverage
    average_zero_df['GenomeCov'] = (100 - (average_zero_df['zeroCov'] * 100/ ref_length))

    # Add sample names
    coverage_df = depth_read_df.join(average_zero_df)

    return coverage_df

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of raw reads fastq-scan json files
    raw_json_files = sorted(glob.glob('*.raw.json'))

    ## Create list of trimmed reads fastq-scan json files
    trim_json_files = sorted(glob.glob('*.trim.json'))

    ## Create list of samtools mapped reads files
    mapped_read_files = sorted(glob.glob('*.map.reads.tsv'))

    ## Create list of samtools depth files
    depth_files = sorted(glob.glob('*.depth.tsv'))

    ## Create dataframe of raw reads fastq-scan results
    raw_json_df = json_to_dataframe(raw_json_files)
    raw_json_df = raw_json_df.rename(columns = {'read_total' : 'NumRawReads'})

    ## Create dataframe of trimmed reads fastq-scan results
    trim_json_df = json_to_dataframe(trim_json_files)
    trim_json_df = trim_json_df.rename(columns = {'read_total' : 'NumTrimReads'})

    ## Merge fastq-scan dataframes
    fastqscan_merged = pd.merge(raw_json_df, trim_json_df, on = ['Sample'])
    fastqscan_merged['%afterTrim'] = fastqscan_merged['NumTrimReads'] / fastqscan_merged['NumRawReads'] * 100

    ## Create dataframe of samtools mapped read results   
    mapped_read_merged_df = mapped_to_dataframe(mapped_read_files)

    ## Merge fastq-scan and mapped read dataframes
    fastqscan_mapped_read_df = pd.merge(fastqscan_merged, mapped_read_merged_df, on = ['Sample'])

    ## Calculate % of reads that mapped
    fastqscan_mapped_read_df['%Mapped'] = fastqscan_mapped_read_df['NumMappedReads'].astype(int) / fastqscan_mapped_read_df['NumTrimReads'].astype(int) * 100

    ## Create dataframe of samtools depth results
    coverage_df = depth_to_dataframe(depth_files, args.ref_length)

    ## Merge fastqscan_mapped and depth

    merged_df = pd.merge(fastqscan_mapped_read_df, coverage_df, on = ['Sample'])

    ## Add LowQualData/Pass/Contaminated/InsufficientData/CheckRequired flag

    flags = ['LowQualData','Pass','Contaminated','InsufficientData']

    ## Set logic for applying flags based on values provided
    conditions = [
        (merged_df['%afterTrim'] < args.min_aft_trim),
        (merged_df['MeanDepth'] >= args.min_depth) & (merged_df['%Mapped'] > args.min_pc),
        (merged_df['MeanDepth'] < args.min_depth) & (merged_df['%Mapped'] < args.min_pc) & (merged_df['NumTrimReads'] > args.min_reads),
        (merged_df['MeanDepth'] < args.min_depth) & (merged_df['NumTrimReads'] < args.min_reads)
    ]

    ## Add flags to merged dataframe
    merged_df['Outcome'] = np.select(conditions, flags, default = 'CheckRequired')

    ## Write output file
    merged_df.to_csv(args.output_file, sep = ',', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())
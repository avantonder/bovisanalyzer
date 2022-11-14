#!/usr/bin/env python3

import os, argparse, sys, subprocess
import pandas as pd
from functools import reduce

def parser_args(args=None):
    """ 
    Function for input arguments for metadata_collate.py
    """
    Description = 'Collect metadata from various outputs and create a final metadata summary table'
    Epilog = """Example usage: python metadata_collate.py --spoligo_summary spoligotype_summary.tsv --tbprofiler_summary tbprofiler.txt --bracken_summary species_composition.tsv --seqtk_summary mapping_summary.tsv --cluster_summary apha_cluster_summary.tsv"""
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-sp", "--spoligo_summary"   , type=str, default="spoligotype_summary.tsv" , help="Spoligotype summary file (default: 'spoligotype_summary.tsv').")
    parser.add_argument("-tb", "--tbprofiler_summary", type=str, default="tbprofiler.txt"          , help="TBprofiler summary file (default: 'tbprofiler.txt').")
    parser.add_argument("-br", "--bracken_summary"   , type=str, help="Bracken summary file (default: None).")
    parser.add_argument("-sq", "--seqtk_summary"     , type=str, default="mapping_summary.tsv"     , help="seqtk summary file (default: 'mapping_summary.tsv').")
    parser.add_argument("-cl", "--cluster_summary"   , type=str, help="APHA clusters summary file (default: None).")
    parser.add_argument("-of", "--output_file"       , type=str, default="metadata_summary.tsv"    , help="Metadata summary file (default: 'metadata_summary.tsv').")
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

def spoligo_filter(spoligo_summary):
    """ 
    Function for reading and filtering spoligotype summary table
    """
    table = pd.read_csv(spoligo_summary, sep='\t')
    table = table.iloc[:,0:2]
    
    return table

def tbprofiler_filter(tbprofiler_summary):
    """ 
    Function for reading and filtering TBprofiler summary table
    """
    table = pd.read_csv(tbprofiler_summary, sep='\t')
    table = table.iloc[:,0:3]
    
    return table

def bracken_filter(bracken_summary):
    """ 
    Function for reading and filtering Bracken summary table
    """
    table = pd.read_csv(bracken_summary, sep='\t')
    table = table[['name','Mycobacterium tuberculosis']]
    table.rename(columns = {'name': 'sample'}, inplace=True)
    
    return table

def seqtk_filter(seqtk_summary):
    """ 
    Function for reading and filtering seqtk summary table
    """
    table = pd.read_csv(seqtk_summary, sep='\t')
    table = table.iloc[:,[0,7]]
    
    return table

def cluster_filter(cluster_summary):
    """ 
    Function for reading and filtering APHA cluster summary table
    """
    table = pd.read_csv(cluster_summary, sep='\t')
    table = table[['Sample','MeanDepth', 'pcMapped', 'Outcome', 'flag', 'group']]
    table.rename(columns = {'Sample': 'sample'}, inplace=True)
    
    return table

def metadata_merge(metadata_list, keyname):
    """ 
    Function for merging a list of dataframes using the column 'sample'
    """
    merged = reduce(lambda  left,right: pd.merge(left,right,on=[keyname], how='outer'), metadata_list)

    return merged

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Read and filter spoligotype summary
    spoligo_summary_df = spoligo_filter(args.spoligo_summary)

    ## Read and filter TBprofiler summary
    tbprofiler_summary_df = tbprofiler_filter(args.tbprofiler_summary)

    ## Read and filter Bracken summary
    if args.bracken_summary is not None:
        bracken_summary_df = bracken_filter(args.bracken_summary)

    ## Read and filter seqtk summary
    seqtk_summary_df = seqtk_filter(args.seqtk_summary)

    ## Read and filter cluster summary
    if args.cluster_summary is not None:
        cluster_summary_df = cluster_filter(args.cluster_summary)
    
    ## Create metadata list
    metadata_list = []
    if args.bracken_summary is not None and args.cluster_summary is not None:
        metadata_list = [spoligo_summary_df, tbprofiler_summary_df, bracken_summary_df, seqtk_summary_df, cluster_summary_df]
    elif args.bracken_summary is not None:
        metadata_list = [spoligo_summary_df, tbprofiler_summary_df, bracken_summary_df, seqtk_summary_df]
    elif args.cluster_summary is not None:
        metadata_list = [spoligo_summary_df, tbprofiler_summary_df, seqtk_summary_df, cluster_summary_df]
      
    ## Merge metadata dataframes
    merged_table = metadata_merge(metadata_list, 'sample')

    ## Define columns for re-ordering based on composition of metadata dataframe
    columns_to_keep = []
    if args.bracken_summary is not None and args.cluster_summary is not None:
        columns_to_keep = ['sample', 'Mycobacterium tuberculosis', 'main_lineage', 'sub_lineage', 'flag', 'group', 'SB number', 'MeanDepth', 'pcMapped', '%ref mapped', 'Outcome']
    elif args.bracken_summary is not None:
        columns_to_keep = ['sample', 'Mycobacterium tuberculosis', 'main_lineage', 'sub_lineage', 'SB number', '%ref mapped']
    elif args.cluster_summary is not None:
        columns_to_keep = ['sample', 'main_lineage', 'sub_lineage', 'flag', 'group', 'SB number', 'MeanDepth', 'pcMapped', '%ref mapped', 'Outcome']
    
    ## Re-order columns in merged metadata dataframe
    merged_table = merged_table[columns_to_keep]

    ## Write output file
    merged_table.to_csv(args.output_file, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())
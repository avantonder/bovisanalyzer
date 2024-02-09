#!/usr/bin/env python3

import os 
import argparse 
import sys 
import glob
import pandas as pd

def parser_args(args=None):
    """ 
    Function for input arguments for spoligo_parser.py
    """
    Description = 'Collect SpoTyping outputs, match against spoligotype database and create a summary table'
    Epilog = """Example usage: python spoligo_parser.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-tb", "--tbprofiler_summary", type=str, default="tbprofiler.txt"         , help="TB-profiler summary file (default: 'tbprofiler.txt').")
    parser.add_argument("-sb", "--spol_db"           , type=str, default="spoligotype_db.tsv"     , help="Spoligotype database file (default: 'spoligotype_db.tsv').")
    parser.add_argument("-of", "--output_file"       , type=str, default="spoligotype_summary.tsv", help="Spoligotype summary file (default: 'spoligotype_summary.tsv').")
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

def spol_db_parse(spol_db):
    """ 
    Function for reading and formatting spoligotype database
    """
    spoligo_db_read = pd.read_csv(spol_db, sep=' ', header=None)
    spoligo_db_df = spoligo_db_read.iloc[:, [0,1]]
    spoligo_db_df.columns = ['octal', 'SB number']

    return spoligo_db_df

def tbprofiler_parse(tbprofiler_summary):
    """ 
    Function to read TB-profiler summary file and extract octal code 
    to compare against spoligotype database
    """
    tbprofiler_read = pd.read_csv(tbprofiler_summary, sep='\t')
    tbprofiler_df = tbprofiler_read.iloc[:, [0,3]]
    tbprofiler_df.columns = ['sample', 'octal']

    return tbprofiler_df
    
def spoligo_match(spoligo_db_df, tbprofiler_df):
    """ 
    Function to match TB-profiler octal code against spoligotype database
    """
    spoligo_lookup_df = pd.merge(
        left = tbprofiler_df,
        right = spoligo_db_df,
        left_on = 'octal',
        right_on = 'octal',
        how = 'left'
    ).fillna('Not found')
    
    return spoligo_lookup_df
    
def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Parse and format TB-profiler summary file
    tbprofiler_df = tbprofiler_parse(args.tbprofiler_summary)

    ## Parse and format spoligotype database
    spoligo_db = spol_db_parse(args.spol_db)
    
    ## Merge spotyping results and database
    spoligo_lookup_df = spoligo_match(spoligo_db, tbprofiler_df)

    ## Final column order
    cols = ['sample', 'SB number', 'octal']
    spoligo_ordered_df = spoligo_lookup_df[cols]

    ## Write output file
    spoligo_ordered_df.to_csv(args.output_file, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())
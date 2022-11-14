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
    parser.add_argument("-sb", "--spol_db"    , type=str, default="spoligotype_db.tsv"     , help="Spoligotype database file (default: 'spoligotype_db.tsv').")
    parser.add_argument("-of", "--output_file", type=str, default="spoligotype_summary.tsv", help="Spoligotype summary file (default: 'spoligotype_summary.tsv').")
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
    spoligo_db_df = spoligo_db_read.iloc[:, [2,1]]
    spoligo_db_df.columns = ['binary', 'SB number']

    return spoligo_db_df

def spol_to_dataframe(file_list):
    """ 
    Function to take list of SpoTyping files and create a summary table
    """
    spoligo_names = [i.replace('.txt', '') for i in file_list]
    spoligo_names_df = pd.DataFrame(spoligo_names)
    spoligo_names_df.columns = ['sample']
    spoligo_files_read = [open(file).read() for file in file_list]
    spoligo_files_read_df = pd.DataFrame(spoligo_files_read)
    spoligo_files_read_df = spoligo_files_read_df[0].str.split('\t', expand=True)
    spoligo_files_read_df = spoligo_files_read_df.replace('\n','', regex=True)
    spoligo_files_read_df.columns = ['files', 'binary', 'octal']
    spoligo_merged_df = spoligo_names_df.join(spoligo_files_read_df)
    cols = ['sample', 'binary', 'octal']
    spoligo_ordered_df = spoligo_merged_df[cols]

    return spoligo_ordered_df

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of SpoTyping output files
    spoligo_files = sorted(glob.glob('*.txt'))

    ## Create dataframe
    spoligo_df = spol_to_dataframe(spoligo_files)

    ## Parse and format spoligotype database
    spoligo_db = spol_db_parse(args.spol_db)
    
    ## Merge spotyping results and database
    spoligo_lookup_df = pd.merge(
        left=spoligo_df,
        right=spoligo_db,
        left_on='binary',
        right_on='binary',
        how='left'
    ).fillna('Not found')

    ## Final column order
    cols = ['sample', 'SB number', 'binary', 'octal']
    spoligo_ordered_df = spoligo_lookup_df[cols]

    ## Write output file
    spoligo_ordered_df.to_csv(args.output_file, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())
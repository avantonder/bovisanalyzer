#!/usr/bin/env python

import pandas as pd
import sys
import glob
import os
import argparse

def parser_args(args=None):
    """ 
    Function for input arguments for cluster_parser.py
    """
    Description = 'Collect APHA cluster outputs and create a summary table'
    Epilog = """Example usage: python cluster_parser.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-of", "--output_file", type=str, default="apha_cluster_summary.tsv", help="Cluster summary file (default: 'apha_cluster_summary.tsv').")
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

def cluster_parse(files):
    """ 
    Function for creating a dataframe from individual csv files
    """
    file_list = [pd.read_csv(f, sep=',') for f in files]
    dataframe = pd.concat(file_list, ignore_index=True)

    return dataframe

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of cluster files
    cluster_files = sorted(glob.glob('*.csv'))

    ## Create dataframe
    cluster_df = cluster_parse(cluster_files)

    ## Write output file
    cluster_df.to_csv(args.output_file, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())
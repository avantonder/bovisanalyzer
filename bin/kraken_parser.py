#!/usr/bin/env python

import sys
import glob
import os
import argparse
import pandas as pd
import numpy as np

def parser_args(args=None):
    """ 
    Function for input arguments for kraken_parser.py
    """
    Description = 'Collect Kraken 2 and Bracken outputs and create a summary table'
    Epilog = """Example usage: python kraken_parser.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-of", "--output_file", type=str, default="species_composition.tsv", help="Species composition file (default: 'species_composition.tsv').")
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

def files_to_dataframe(kraken_files, bracken_files):
	""" 
    Function for making a dataframe from lists of Kraken 2 and Bracken report files
    """
	kraken_report_list = [pd.read_csv(f, sep='\t', header=None) for f in kraken_files]
	kraken_report_list = [i.head(1).iloc[:,[5,1]] for i in kraken_report_list]
	species_abundance_list = [pd.read_csv(f, sep='\t') for f in bracken_files]
	species_abundance_names = [i.replace('_S.tsv', '') for i in bracken_files]
	
	for a,b in zip(kraken_report_list, species_abundance_names):
		a.rename(columns = {5: 'name', 1: b}, inplace=True)

	for a,b in zip(species_abundance_list, species_abundance_names):
		a.rename(columns = {'new_est_reads': b, 'fraction_total_reads': b + '_frac'}, inplace=True)
	
	species_abundance_list = [i.iloc[:,[0,5]] for i in species_abundance_list]
	species_abundance_joined = pd.concat([i.set_index('name') for i in species_abundance_list], axis=1).reset_index()
	species_abundance_joined.rename(columns = {'index': 'name'}, inplace=True)
	kraken_report_joined = pd.concat([i.set_index('name') for i in kraken_report_list], axis=1).reset_index()
	species_abundance_unclassified = species_abundance_joined.append(kraken_report_joined, ignore_index=True)
	species_abundance_unclassified.fillna(0, inplace=True)

	for column in species_abundance_unclassified.columns[1:]:
		species_abundance_unclassified[column + '_freq'] = species_abundance_unclassified[column] / species_abundance_unclassified[column].sum() * 100

	species_abundance_unclassified_filtered = pd.concat([species_abundance_unclassified['name'], species_abundance_unclassified.filter(like='_freq')], axis=1)
	species_abundance_unclassified_filtered['Max_value'] = species_abundance_unclassified_filtered.filter(like='_freq').max(axis=1)
	species_abundance_unclassified_filtered = species_abundance_unclassified_filtered[(species_abundance_unclassified_filtered['Max_value'] >= 5)]
	total = 100 - species_abundance_unclassified_filtered.filter(like='_freq').apply(np.sum)
	total['name'] = 'other'
	final_abundance = species_abundance_unclassified_filtered.append(pd.DataFrame(total.values, index=total.keys()).T, ignore_index=True)
	final_abundance_sorted = pd.concat([final_abundance['name'], final_abundance.filter(like = '_freq')], axis = 1)
	final_abundance_sorted.columns = final_abundance_sorted.columns.str.replace('_freq', '')

	return final_abundance_sorted

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of Kraken 2 report files
    kraken_report_files = sorted(glob.glob('*.kraken2.report.txt'))

    ## Create list of Bracken species abundance files
    species_abundance_files = sorted(glob.glob('*.tsv'))
	
	## Create combined Kraken 2 and Bracken dataframe
    combined_df = files_to_dataframe(kraken_report_files,species_abundance_files)

    ## Write output file
    combined_df.T.to_csv(args.output_file, sep = '\t', header = False)

if __name__ == '__main__':
    sys.exit(main())
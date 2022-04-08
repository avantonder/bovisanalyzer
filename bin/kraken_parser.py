#!/usr/bin/env python

import os, argparse, sys, subprocess, glob
import pandas as pd
import numpy as np

kraken_report_files = sorted(glob.glob('*.kraken2.report.txt'))

species_abundance_files = sorted(glob.glob('*_output_species_abundance.txt'))

kraken_report_df = [pd.read_csv(f, sep='\t', header=None) for f in kraken_report_files]

kraken_report_df_filtered = [i.head(1).iloc[:,[5,1]] for i in kraken_report_df]

species_abundance_df = [pd.read_csv(f, sep='\t') for f in species_abundance_files]

species_abundance_names = [i.replace('_output_species_abundance.txt', '') for i in species_abundance_files]

for a,b in zip(kraken_report_df_filtered, species_abundance_names):
	a.rename(columns = {5: 'name', 1: b}, inplace=True)

for a,b in zip(species_abundance_df, species_abundance_names):
	a.rename(columns = {'new_est_reads': b, 'fraction_total_reads': b + '_frac'}, inplace=True)

species_abundance_df_filtered = [i.iloc[:,[0,5]] for i in species_abundance_df]

species_abundance_joined = pd.concat([i.set_index('name') for i in species_abundance_df_filtered], axis=1).reset_index()

species_abundance_joined.rename(columns = {'index': 'name'}, inplace=True)

kraken_report_joined = pd.concat([i.set_index('name') for i in kraken_report_df_filtered], axis=1).reset_index()

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

final_abundance_tsv_name = 'Bracken_species_composition.tsv'

final_abundance_sorted.T.to_csv(final_abundance_tsv_name, sep = '\t', header = False)
#!/usr/bin/env python3

import os, argparse, sys, subprocess, glob
import pandas as pd
import numpy as np

# Read in spoligotype database file

spoligo_db = "/../assets/spoligotype_db.txt"

spoligo_db_read = pd.read_csv(spoligo_db, sep=' ', header=None)

spoligo_db_df = spoligo_db_read.iloc[:, [2,1]]

spoligo_db_df.columns = ['binary', 'SB number']

# Read in spotyping output files

spoligo_files = sorted(glob.glob('*.txt'))

spoligo_names = [i.replace('.txt', '') for i in spoligo_files]

spoligo_names_df = pd.DataFrame(spoligo_names)

spoligo_names_df.columns = ['sample']

spoligo_files_read = [open(file).read() for file in spoligo_files]

spoligo_files_read_df = pd.DataFrame(spoligo_files_read)

spoligo_files_read_df = spoligo_files_read_df[0].str.split('\t', expand=True)

spoligo_files_read_df = spoligo_files_read_df.replace('\n','', regex=True)

spoligo_files_read_df.columns = ['files', 'binary', 'octal']

spoligo_merged_df = spoligo_names_df.join(spoligo_files_read_df)

cols = ['sample', 'binary', 'octal']

spoligo_ordered_df = spoligo_merged_df[cols]

# Merge spotyping results and database

spoligo_lookup_df = pd.merge(
    left=spoligo_ordered_df,
    right=spoligo_db_df,
    left_on='binary',
    right_on='binary',
    how='left'
).fillna('Not found')

cols = ['sample', 'SB number', 'binary', 'octal']

spoligo_ordered_df = spoligo_lookup_df[cols]

# Save to tsv file

spoligo_tsv_name = 'Spoligotype_summary.tsv'

spoligo_ordered_df.to_csv(spoligo_tsv_name, sep = '\t', header = True, index = False)
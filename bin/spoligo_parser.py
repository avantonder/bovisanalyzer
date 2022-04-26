#!/usr/bin/env python3

import os, argparse, sys, subprocess, glob
import pandas as pd
import numpy as np

spoligo_files = sorted(glob.glob('*.txt'))

spoligo_names = [i.replace('.txt', '') for i in spoligo_files]

spoligo_names_df = pd.DataFrame(spoligo_names)

spoligo_names_df.columns = ['sample']

spoligo_files_read = [open(file).read() for file in spoligo_files]

spoligo_files_read_df = pd.DataFrame(spoligo_files_read)

spoligo_files_read_df = spoligo_files_read_df[0].str.split(',', expand=True)

spoligo_files_read_df = spoligo_files_read_df.replace('\n','', regex=True)

spoligo_files_read_df.columns = ['octal', 'SB number', 'binary']

spoligo_merged_df = spoligo_names_df.join(spoligo_files_read_df)

cols = ['sample', 'SB number', 'binary', 'octal']

spoligo_ordered_df = spoligo_merged_df[cols]

spoligo_tsv_name = 'Spoligotype_summary.tsv'

spoligo_ordered_df.to_csv(spoligo_tsv_name, sep = '\t', header = True, index = False)
#!/usr/bin/env python3

import glob
import pandas as pd

# Read in seqtk comp tsv outputs

seqtk_files = sorted(glob.glob('*.tsv'))

seqtk_file_list = [pd.read_csv(f, sep='\t', header=None) for f in seqtk_files]

seqtk_df = pd.concat(seqtk_file_list, ignore_index=True)

seqtk_df = seqtk_df.iloc[:, 0:6]

seqtk_df.columns = ['sample', 'ref_length', '#A', '#C', '#G', '#T']

column_list = ['#A', '#C', '#G', '#T']

seqtk_df["mapped"] = seqtk_df[column_list].sum(axis=1)

seqtk_df["%ref mapped"] = seqtk_df['mapped'] / seqtk_df['ref_length'] * 100

seqtk_tsv_name = 'mapping_summary.tsv'

seqtk_df.to_csv(seqtk_tsv_name, sep = '\t', header = True, index = False)
#!/usr/bin/env python3

import os, argparse, sys, subprocess
import pandas as pd
import numpy as np

spoligo_summary_df = pd.read_csv('Spoligotype_summary.tsv', sep='\t')

tbprofiler_summary_df = pd.read_csv('tbprofiler.txt', sep='\t')

bracken_summary_df = pd.read_csv('Bracken_species_composition.tsv', sep='\t')

seqtk_summary_df = pd.read_csv('mapping_summary.tsv', sep='\t')

spoligo_summary_df_filter = spoligo_summary_df.iloc[:,0:2]

tbprofiler_summary_df_filter = tbprofiler_summary_df.iloc[:,0:3]

bracken_summary_df_filter = bracken_summary_df[['name','Mycobacterium tuberculosis']]

bracken_summary_df_filter = bracken_summary_df_filter.rename(columns = {'name': 'sample'})

seqtk_summary_df_filter = seqtk_summary_df.iloc[:,[0,7]]

bracken_tbprofiler_merge = pd.merge(bracken_summary_df_filter, tbprofiler_summary_df_filter, on = ['sample'])

spoligo_merge = pd.merge(bracken_tbprofiler_merge, spoligo_summary_df_filter, on = ['sample'])

all_merge = pd.merge(spoligo_merge, seqtk_summary_df_filter, on = ['sample'])

summary_file_name = 'Metadata_summary.tsv'

all_merge.to_csv(summary_file_name, sep = '\t', header = True, index = False)
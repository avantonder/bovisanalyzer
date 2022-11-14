#!/usr/bin/env python3

import os, argparse, sys, subprocess
import pandas as pd
import numpy as np

spoligo_summary_df = pd.read_csv('spoligotype_summary.tsv', sep='\t')

tbprofiler_summary_df = pd.read_csv('tbprofiler.txt', sep='\t')

bracken_summary_df = pd.read_csv('species_composition.tsv', sep='\t')

seqtk_summary_df = pd.read_csv('mapping_summary.tsv', sep='\t')

cluster_summary_df = pd.read_csv('apha_cluster_summary.tsv', sep='\t')

spoligo_summary_df_filter = spoligo_summary_df.iloc[:,0:2]

tbprofiler_summary_df_filter = tbprofiler_summary_df.iloc[:,0:3]

bracken_summary_df_filter = bracken_summary_df[['name','Mycobacterium tuberculosis']]

bracken_summary_df_filter = bracken_summary_df_filter.rename(columns = {'name': 'sample'})

seqtk_summary_df_filter = seqtk_summary_df.iloc[:,[0,7]]

cluster_summary_df_filter = cluster_summary_df[['Sample','MeanDepth', 'pcMapped', 'Outcome', 'flag', 'group']]

cluster_summary_df_filter = cluster_summary_df_filter.rename(columns = {'Sample': 'sample'})

bracken_tbprofiler_merge = pd.merge(bracken_summary_df_filter, tbprofiler_summary_df_filter, on = ['sample'])

spoligo_merge = pd.merge(bracken_tbprofiler_merge, spoligo_summary_df_filter, on = ['sample'])

spoligo_cluster_merge = pd.merge(spoligo_merge, cluster_summary_df_filter, on = ['sample'])

all_merge = pd.merge(spoligo_cluster_merge, seqtk_summary_df_filter, on = ['sample'])

all_merge = all_merge[['sample', 'Mycobacterium tuberculosis', 'main_lineage', 'sub_lineage', 'flag', 'group', 'SB number', 'MeanDepth', 'pcMapped', '%ref mapped', 'Outcome']]

summary_file_name = 'metadata_summary.tsv'

all_merge.to_csv(summary_file_name, sep = '\t', header = True, index = False)
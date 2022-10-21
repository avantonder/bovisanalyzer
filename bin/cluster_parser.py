#!/usr/bin/env python

import pandas as pd
import glob

cluster_files = sorted(glob.glob('*.csv'))

cluster_file_list = [pd.read_csv(f, sep=',') for f in cluster_files]

cluster_df = pd.concat(cluster_file_list, ignore_index=True)

cluster_tsv_name = 'apha_cluster_summary.tsv'

cluster_df.to_csv(cluster_tsv_name, sep = '\t', header = True, index = False)
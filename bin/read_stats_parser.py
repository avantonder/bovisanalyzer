#!/usr/bin/env python

import pandas as pd
import glob

readstats_files = sorted(glob.glob('*.csv'))

readstats_file_list = [pd.read_csv(f, sep=',') for f in readstats_files]

readstats_df = pd.concat(readstats_file_list, ignore_index=True)

readstats_tsv_name = 'read_stats_summary.tsv'

readstats_df.to_csv(readstats_tsv_name, sep = '\t', header = True, index = False)
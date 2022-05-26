#!/usr/bin/env python

import pandas as pd
import json
import glob

# Read in fastqscan json files

json_files = sorted(glob.glob('*.json'))

json_names = [i.replace('.json', '') for i in json_files]

json_names_df = pd.DataFrame(json_names)

json_names_df.columns = ['sample']

jsons_data = {}

for index, file in enumerate(json_files):
    with open(file, 'r') as f:
        json_text = json.loads(f.read())
        qc = json_text['qc_stats']
        jsons_data[index] = qc

jsons_data_df = pd.DataFrame.from_dict(jsons_data, orient = 'index')

json_merged_df = json_names_df.join(jsons_data_df)

# Save to tsv file

json_tsv_name = 'fastq-scan_summary.tsv'

json_merged_df.to_csv(json_tsv_name, sep = '\t', header = True, index = False)
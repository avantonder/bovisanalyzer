#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob
import json

# Reference length

ref_length = 4349904

# Define thresholds for flag assignment

mindepth = 10 # require an average of at least 10 reads per site 

minpc = 60 # require at least 60% of data maps to genome

minreads = 600000 # require at least 600,000 raw reads per sample

minafttrim = 60 # require at least 60% reads to pass quality filtering steps

# Parse raw fastq-scan results

# Read in fastqscan json files

raw_json_file = sorted(glob.glob('*.raw.json'))

raw_json_names = [i.replace('.raw.json', '') for i in raw_json_file]

raw_json_names_df = pd.DataFrame(raw_json_names)

raw_json_names_df.columns = ['Sample']

raw_jsons_data = {}

for index, file in enumerate(raw_json_file):
    with open(file, 'r') as f:
        raw_json_text = json.loads(f.read())
        qc = raw_json_text['qc_stats']
        raw_jsons_data[index] = qc

raw_jsons_data_df = pd.DataFrame.from_dict(raw_jsons_data, orient = 'index')

raw_json_merged_df = raw_json_names_df.join(raw_jsons_data_df)

fastqscan_raw_df = raw_json_merged_df.iloc[:, [0,3]]

fastqscan_raw_df = fastqscan_raw_df.rename(columns = {'read_total' : 'NumRawReads'})

# Parse trimmed fastq-scan results

trim_json_file = sorted(glob.glob('*.trim.json'))

trim_json_names = [i.replace('.trim.json', '') for i in trim_json_file]

trim_json_names_df = pd.DataFrame(trim_json_names)

trim_json_names_df.columns = ['Sample']

trim_jsons_data = {}

for index, file in enumerate(trim_json_file):
    with open(file, 'r') as f:
        trim_json_text = json.loads(f.read())
        qc = trim_json_text['qc_stats']
        trim_jsons_data[index] = qc

trim_jsons_data_df = pd.DataFrame.from_dict(trim_jsons_data, orient = 'index')

trim_json_merged_df = trim_json_names_df.join(trim_jsons_data_df)

fastqscan_trim_df = trim_json_merged_df.iloc[:, [0,3]]

fastqscan_trim_df = fastqscan_trim_df.rename(columns = {'read_total' : 'NumTrimReads'})

# Merge fastq-scan dataframes

fastqscan_merged = pd.merge(fastqscan_raw_df, fastqscan_trim_df, on = ['Sample'])

fastqscan_merged['%afterTrim'] = fastqscan_merged['NumTrimReads'] / fastqscan_merged['NumRawReads'] * 100

# Parse samtools mapped read counts

mapped_read_files = sorted(glob.glob('*.map.reads.tsv'))

mapped_read = [i.replace('.map.reads.tsv', '') for i in mapped_read_files]

mapped_read_df = pd.DataFrame(mapped_read)

mapped_read_df.columns = ['Sample']

mapped_read_files_read = [open(file).read() for file in mapped_read_files]

mapped_read_files_df = pd.DataFrame(mapped_read_files_read)

mapped_read_files_df.columns = ['NumMappedReads']

mapped_read_merged_df = mapped_read_df.join(mapped_read_files_df).replace('\n','', regex=True)

# Merge fastq-scan and mapped read dataframes

fastqscan_mapped_read_df = pd.merge(fastqscan_merged, mapped_read_merged_df, on = ['Sample'])

# Calculate % of reads that mapped

fastqscan_mapped_read_df['%Mapped'] = fastqscan_mapped_read_df['NumMappedReads'].astype(int) / fastqscan_mapped_read_df['NumTrimReads'].astype(int) * 100

# Parse samtools depth files

depth_read_files = sorted(glob.glob('*.depth.tsv'))

depth_read = [i.replace('.depth.tsv', '') for i in depth_read_files]

depth_read_df = pd.DataFrame(depth_read)

depth_read_df.columns = ['Sample']

depth_read_files_read = [pd.read_csv(f, sep='\t', header=None) for f in depth_read_files]

# Calculate average mapped read depth

average_read_depth = [i.iloc[:,[2]].mean() for i in depth_read_files_read]

average_read_depth_df = pd.DataFrame(average_read_depth)

average_read_depth_df.columns = ['MeanDepth']

# Calculate number of bases with zero coverage

zero_count = [(i.iloc[:,[2]] == 0).sum() for i in depth_read_files_read]

zero_count_df = pd.DataFrame(zero_count)

zero_count_df.columns = ['zeroCov']

# Merge average and zero counts

average_zero_df = average_read_depth_df.join(zero_count_df)

# Calculate genome coverage

average_zero_df['GenomeCov'] = (100 - (average_zero_df['zeroCov'] * 100/ ref_length))

# Add sample names

coverage_df = depth_read_df.join(average_zero_df)

# Merge fastqscan_mapped and depth

merged_df = pd.merge(fastqscan_mapped_read_df, coverage_df, on = ['Sample'])

# Add LowQualData/Pass/Contaminated/InsufficientData/CheckRequired flag

flags = ['LowQualData','Pass','Contaminated','InsufficientData']

conditions = [
    (merged_df['%afterTrim'] < minafttrim),
    (merged_df['MeanDepth'] >= mindepth) & (merged_df['%Mapped'] > minpc),
    (merged_df['MeanDepth'] < mindepth) & (merged_df['%Mapped'] < minpc) & (merged_df['NumTrimReads'] > minreads),
    (merged_df['MeanDepth'] < mindepth) & (merged_df['NumTrimReads'] < minreads)
]

merged_df['Outcome'] = np.select(conditions, flags, default = 'CheckRequired')

# Write merged dataframe to csv file

summary_file_name = 'read_stats.csv'

merged_df.to_csv(summary_file_name, sep = ',', header = True, index = False)
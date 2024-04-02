#!/usr/bin/env python3

import sys
import os
import glob
import pandas as pd
from itertools import chain

snp_dirs = [line.strip() for line in open(sys.argv[1], 'r')]
if len(snp_dirs) < 1:
    sys.exit(0)

output_directory = os.path.abspath(sys.argv[2])
isolate_data = pd.read_csv(sys.argv[3],sep="\t")
mummer_data = pd.read_csv(sys.argv[4],sep="\t")

raw_snp_distance_files = list(chain.from_iterable([glob.glob(snp_dir + '/*snp_distance_pairwise.tsv') for snp_dir in snp_dirs]))
preserved_snp_distance_files = list(chain.from_iterable([glob.glob(snp_dir + '/*snp_distance_pairwise_preserved.tsv') for snp_dir in snp_dirs]))

# Read in files
if len(raw_snp_distance_files) == 0:
    raw_snp_distance_df = pd.DataFrame()
    raw_stats = pd.DataFrame()
    raw_isolates = []
else:
    raw_snp_distance_df = pd.concat([pd.read_csv(file, sep='\t').assign(Reference_ID=os.path.basename(os.path.dirname(file))) for file in raw_snp_distance_files])
    raw_snp_distance_df['Comparison'] = raw_snp_distance_df.apply(lambda row: ';'.join(sorted([str(row['Query_1']), str(row['Query_2'])])), axis=1)
    raw_isolates = list(set(raw_snp_distance_df['Query_1'].astype(str) + raw_snp_distance_df['Query_2'].astype(str)))

if len(preserved_snp_distance_files) == 0:
    preserved_snp_distance_df = pd.DataFrame()
    preserved_stats = pd.DataFrame()
    preserved_isolates = []
else:
    preserved_snp_distance_df = pd.concat([pd.read_csv(file, sep='\t').assign(Reference_ID=os.path.basename(os.path.dirname(file))) for file in preserved_snp_distance_files])
    preserved_snp_distance_df['Comparison'] = preserved_snp_distance_df.apply(lambda row: ';'.join(sorted([str(row['Query_1']), str(row['Query_2'])])), axis=1)
    preserved_isolates = list(set(preserved_snp_distance_df['Query_1'].astype(str) + preserved_snp_distance_df['Query_2'].astype(str)))


# Summarize data
if len(raw_snp_distance_files) > 1:
    raw_snp_distance_stats = raw_snp_distance_df.groupby('Comparison')['SNP_Distance'].agg(Count=('count'), Min=('min'), Mean=('mean'), Max=('max'), Spread=(lambda x: x.max() - x.min())).astype(int).reset_index()
    raw_cocalled_stats = raw_snp_distance_df.groupby('Comparison')['SNPs_Cocalled'].agg(Count=('count'), Min=('min'), Mean=('mean'), Max=('max'), Spread=(lambda x: x.max() - x.min())).astype(int).reset_index()
    raw_stats = pd.concat([raw_snp_distance_stats.assign(Measure='SNPs'), raw_cocalled_stats.assign(Measure='Cocalled')])
    raw_stats[['Query_1', 'Query_2']] = raw_stats['Comparison'].str.split(';', expand=True)


if len(preserved_snp_distance_files) > 1:
    preserved_snp_distance_stats = preserved_snp_distance_df.groupby('Comparison')['SNP_Distance'].agg(Count=('count'), Min=('min'), Mean=('mean'), Max=('max'), Spread=(lambda x: x.max() - x.min())).astype(int).reset_index()
    preserved_cocalled_stats = preserved_snp_distance_df.groupby('Comparison')['SNPs_Cocalled'].agg(Count=('count'), Min=('min'), Mean=('mean'), Max=('max'), Spread=(lambda x: x.max() - x.min())).astype(int).reset_index()
    preserved_stats = pd.concat([preserved_snp_distance_stats.assign(Measure='SNPs'), preserved_cocalled_stats.assign(Measure='Cocalled')])
    preserved_stats[['Query_1', 'Query_2']] = preserved_stats['Comparison'].str.split(';', expand=True)

print(isolate_data)
print(mummer_data)
print(raw_snp_distance_df)
print(raw_stats)
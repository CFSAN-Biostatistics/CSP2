#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import argparse

def processHeader(header_row,file_path,trim_name):
    header_cols = [item.split(':')[0] for item in header_row]
    header_vals = [item.split(':')[1] for item in header_row]

    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    
    # Replace trim_name from Reference_ID and Query_ID
    header_data['Reference_ID'] = header_data['Reference_ID'].str.replace(trim_name,'',regex=False)
    header_data['Query_ID'] = header_data['Query_ID'].str.replace(trim_name,'',regex=False)
    
    # Add file path
    header_data['SNPDiffs_File'] = os.path.abspath(file_path)
    cols = header_data.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    header_data = header_data[cols]
    
    return header_data


parser = argparse.ArgumentParser()
parser.add_argument("--snpdiffs_file", help="Path to the SNP diffs list file")
parser.add_argument("--summary_file", help="Path to the summary file")
parser.add_argument("--isolate_file", help="Path to the isolate data file")
parser.add_argument("--trim_name", help="Trim name")
parser.add_argument("--ref_id_file", help="Path to the reference IDs file")
args = parser.parse_args()

snpdiffs_list_file = args.snpdiffs_file
summary_file = args.summary_file
isolate_data_file = args.isolate_file
trim_name = args.trim_name

if os.stat(args.ref_id_file).st_size == 0:
    ref_ids = []
else:
    ref_ids = [line.strip() for line in open(args.ref_id_file, 'r')]

# Read in all lines and ensure each file exists
snpdiffs_list = [line.strip() for line in open(snpdiffs_list_file, 'r')]
snpdiffs_list = [line for line in snpdiffs_list if line]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        sys.exit("Error: File does not exist: " + snpdiffs_file)

header_rows = []
for snpdiffs_file in snpdiffs_list:
    with open(snpdiffs_file, 'r') as file:
        top_line = file.readline().lstrip('#').strip().split('\t')
        header_rows.append(processHeader(top_line,snpdiffs_file,trim_name))

output_data = pd.concat(header_rows, ignore_index=True)
output_data.to_csv(summary_file, sep='\t', index=False)

# If ref_ids is empty, save isolate data
ref_header = ['Reference_ID','Reference_Assembly','Reference_Contig_Count','Reference_Assembly_Bases','Reference_N50','Reference_N90','Reference_L50','Reference_L90','Reference_SHA256']
query_header = ['Query_ID','Query_Assembly','Query_Contig_Count','Query_Assembly_Bases','Query_N50','Query_N90','Query_L50','Query_L90','Query_SHA256']
isolate_header = ["Isolate_ID","Assembly_Path","Contig_Count","Assembly_Bases","N50","N90","L50","L90","SHA256"]

ref_df = output_data[ref_header]
query_df = output_data[query_header]

ref_df.columns = isolate_header
query_df.columns = isolate_header

combined_df = pd.concat([ref_df,query_df])

# Set combined_df[Isolate_Type] to Reference if Isolate_ID is in ref_ids
combined_df['Isolate_Type'] = np.where(combined_df['Isolate_ID'].isin(ref_ids), 'Reference', 'Query')
combined_df = combined_df.drop_duplicates()
cols = combined_df.columns.tolist()
cols = cols[:1] + cols[-1:] + cols[1:-1]
combined_df = combined_df[cols]
combined_df.to_csv(isolate_data_file, sep='\t', index=False)


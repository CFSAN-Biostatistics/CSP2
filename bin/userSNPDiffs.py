#!/usr/bin/env python3

import os
import sys
import pandas as pd
import hashlib

def checkLineExists(file_path, sha256):
    if not os.path.exists(file_path):
        return False
    with open(file_path, 'rb') as file:
        file_sha256 = hashlib.sha256(file.read()).hexdigest()
    return file_sha256 == sha256

def processHeader(header_row,snpdiffs_path):
    header_cols = [item.split(':')[0] for item in header_row]
    header_vals = [item.split(':')[1] for item in header_row]

    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    header_data['SNPDiffs_File'] = snpdiffs_path
    return header_data

snpdiffs_list_file = sys.argv[1]
trim_name = sys.argv[2]
header_rows = []

# Read in all files, and if they exist, read in the header
try:
    snpdiffs_list = [line.strip() for line in open(snpdiffs_list_file, 'r')]
except:
    sys.exit("Error: Unable to read file: " + snpdiffs_list_file)

for snpdiffs_file in snpdiffs_list:
    try:
        snpdiffs_path = os.path.abspath(snpdiffs_file)    
        with open(snpdiffs_path, 'r') as file:
            top_line = file.readline().lstrip('#').strip().split('\t')
            header_rows.append(processHeader(top_line,snpdiffs_path))
    except:
        sys.exit("Error: Unable to read file: " + snpdiffs_file)

# Create header df
try:
    all_snpdiffs_data = pd.concat(header_rows,ignore_index=True)
    all_snpdiffs_data['Reference_ID'] = all_snpdiffs_data['Reference_ID'].str.replace(trim_name,'',regex=False)
    all_snpdiffs_data['Query_ID'] = all_snpdiffs_data['Query_ID'].str.replace(trim_name,'',regex=False)
except:
    sys.exit("Error: Unable to create header dataframe")

query_sha_counts = all_snpdiffs_data.groupby('Query_ID')['Query_SHA256'].nunique()
reference_sha_counts = all_snpdiffs_data.groupby('Reference_ID')['Reference_SHA256'].nunique()
file_counts = all_snpdiffs_data['SNPDiffs_File'].value_counts()

if ((query_sha_counts > 1).any() or (reference_sha_counts > 1).any()):
    print(all_snpdiffs_data[all_snpdiffs_data['Query_ID'].isin(query_sha_counts[query_sha_counts > 1].index)])
    print(all_snpdiffs_data[all_snpdiffs_data['Reference_ID'].isin(reference_sha_counts[reference_sha_counts > 1].index)])
    sys.exit("Multiple SHA256 values found for the same Query_ID/Reference_ID")
elif (file_counts > 1).any():
    print(all_snpdiffs_data[all_snpdiffs_data['SNPDiffs_File'].isin(file_counts[file_counts > 1].index)])
    sys.exit("The same SNPDiffs file is listed multiple times")

else:
    for index, row in all_snpdiffs_data.iterrows():
        query_assembly = os.path.abspath(row['Query_Assembly']) if checkLineExists(row['Query_Assembly'], row['Query_SHA256']) else "null"
        reference_assembly = os.path.abspath(row['Reference_Assembly']) if checkLineExists(row['Reference_Assembly'], row['Reference_SHA256']) else "null"
        
        print(",".join(row['SNPDiffs_File'],
                       row['Query_ID'], query_assembly,row['Query_Contig_Count'],row['Query_Assembly_Bases'],
                       row['Query_N50'],row['Query_N90'],row['Query_L50'],row['Query_L90'],row['Query_SHA256'],
                       row['Reference_ID'],reference_assembly,row['Reference_Contig_Count'],row['Reference_Assembly_Bases'],
                       row['Reference_N50'],row['Reference_N90'],row['Reference_L50'],row['Reference_L90'],row['Reference_SHA256']))
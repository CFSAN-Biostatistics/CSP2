#!/usr/bin/env python3

import os
import sys
import pandas as pd
import concurrent.futures

def fetchHeaders(snpdiffs_file):
    
    with open(snpdiffs_file, 'r') as file:
        top_line = file.readline().strip().split('\t')[1:]

    header_cols = [item.split(':')[0] for item in top_line]
    header_vals = [item.split(':')[1] for item in top_line]
    
    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    return header_data

# Read in list of .snpdiffs
snp_diffs_file = sys.argv[1]
if not os.path.exists(snp_diffs_file):
    sys.exit(f"{snp_diffs_file} does not exist.")

# Ensure files exist
with open(snp_diffs_file, 'r') as file:
    snp_diffs_files = file.readlines()
snp_diffs_files = [line.strip() for line in snp_diffs_files]
if not all(os.path.exists(file_path) for file_path in snp_diffs_files):
    sys.exit(f"One or more files provided by {snp_diffs_file} did not exist")

# Get output directory
output_directory = os.path.normpath(os.path.abspath(sys.argv[2]))

# Read in headers for each file
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(fetchHeaders,snp_diff_file) for snp_diff_file in snp_diffs_files]
headers = [future.result() for future in concurrent.futures.as_completed(results)]
header_df = pd.concat(headers, ignore_index=True)

# Check that each query and reference ID is associated with a single SHA256 hash
reference_sha_check = header_df.groupby('Reference_ID')['Reference_SHA256'].nunique().eq(1).all()
if not reference_sha_check:
    sys.exit("There is a reference ID associated with 2+ unique SHA hashes")

reference_sha_check = header_df.groupby('Reference_SHA256')['Reference_ID'].nunique().eq(1).all()
if not reference_sha_check:
    sys.exit("There is a reference SHA hash associated with 2+ unique reference IDs")    

query_sha_check = header_df.groupby('Query_ID')['Query_SHA256'].nunique().eq(1).all()
if not query_sha_check:
    sys.exit("There is a query ID associated with 2+ unique SHA hashes")

query_sha_check = header_df.groupby('Query_SHA256')['Query_ID'].nunique().eq(1).all()
if not query_sha_check:
    sys.exit("There is a query SHA hash associated with 2+ unique query IDs")

# Check Query x Reference
query_reference_combination_check = header_df.groupby(['Query_ID', 'Reference_ID']).size().eq(1).all()
if not query_reference_combination_check:
    sys.exit("There are multiple entries for one or more query x reference combinations")

# Check QC String
qc_string_check = header_df['QC_String'].nunique() == 1
if not qc_string_check:
    sys.exit("The snpdiffs files provided were generated using different QC parameters")

# Save isolate data
header_df.drop_duplicates(subset='Query_ID', keep='first')[['Query_ID','Query_Assembly','Query_SHA256','Query_Contig_Count','Query_Assembly_Bases']].to_csv(output_directory+"/Query_Isolates.tsv",sep="\t",index=False)
header_df.drop_duplicates(subset='Reference_ID', keep='first')[['Reference_ID','Reference_Assembly','Reference_SHA256','Reference_Contig_Count','Reference_Assembly_Bases']].to_csv(output_directory+"/Reference_Isolates.tsv",sep="\t",index=False)

# Separate queries that passed QC
pass_qc_df = header_df[header_df.Category == "PASS"]
fail_qc_df = header_df[header_df.Category != "PASS"]

pass_qc_df['SNPs'] = pass_qc_df['SNPs'].astype(int)
pass_qc_df = pass_qc_df[['Query_ID','Reference_ID','SNPs','Percent_Query_Aligned','Percent_Reference_Aligned','Median_Percent_Identity','Median_SNP_Percent_Identity','Purged_Alignment', 'Purged_N', 'Purged_Indel', 'Purged_Duplicate', 'Purged_Het', 'Purged_Density', 'Filtered_Edge']].sort_values(by='SNPs')
if pass_qc_df.shape[0] > 0:
    pass_qc_df.to_csv(output_directory+"/Screening_Results_PassQC.tsv",sep="\t",index=False)

fail_qc_df = fail_qc_df[['Query_ID','Percent_Query_Aligned','Reference_ID','Percent_Reference_Aligned','Median_Percent_Identity','Category']]
if fail_qc_df.shape[0] > 0:
    fail_qc_df.to_csv(output_directory+"/Screening_Results_FailQC.tsv",sep="\t",index=False)
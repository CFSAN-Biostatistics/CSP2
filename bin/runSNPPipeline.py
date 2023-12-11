#!/usr/bin/python3

import os
import sys
import pandas as pd
import numpy as np
from glob import glob
from itertools import combinations
from itertools import product
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import concurrent.futures
import time
from collections import Counter
import datetime
from pybedtools import BedTool
import warnings
import ast
import math
from scipy.spatial.distance import pdist, squareform

warnings.filterwarnings("ignore")

def getPairwise(compare_id):

    joined_df = base_df.loc[compare_id.split(";")[0]].merge(base_df.loc[compare_id.split(";")[1]], on='Ref_Loc', how='inner')
    if joined_df.shape[0] == 0:
        return {compare_id:(0,0)}
    else:
        return {compare_id:(joined_df.shape[0],(joined_df['Base_x'] != joined_df['Base_y']).sum())}

def getPrunedPairwise(compare_id):

    joined_df = pruned_base_df.loc[compare_id.split(";")[0]].merge(pruned_base_df.loc[compare_id.split(";")[1]], on='Ref_Loc', how='inner')
    if joined_df.shape[0] == 0:
        return {compare_id:(0,0)}
    else:
        return {compare_id:(joined_df.shape[0],(joined_df['Base_x'] != joined_df['Base_y']).sum())}
    
def makeBED(bed_df):
    bed_df = bed_df.copy() 
    if bed_df.shape[1] == 2:
        bed_df.columns = ['chrom','end']
        bed_df['chrom'] = bed_df['chrom'].astype(str) 
        bed_df['end'] = bed_df['end'].astype(int)
        bed_df['start'] = bed_df['end'] - 1
        bed_df['loc'] = bed_df['chrom'] + '/' + bed_df['end'].astype(str)
        bed_file = BedTool.from_dataframe(bed_df[['chrom','start','end','loc']]).sort()
    elif bed_df.shape[1] == 3:
        bed_df.columns = ['chrom','start','end']
        bed_df['chrom'] = bed_df['chrom'].astype(str) 
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)
        bed_file = BedTool.from_dataframe(bed_df[['chrom','start','end']]).sort()
    elif bed_df.shape[1] == 4:
        bed_df.columns = ['chrom','start','end','query_id']
        bed_df['chrom'] = bed_df['chrom'].astype(str) 
        bed_df['query_id'] = bed_df['query_id'].astype(str) 
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)
        bed_df['start'] = bed_df['start'] - 1
        bed_file = BedTool.from_dataframe(bed_df[['chrom','start','end','query_id']]).sort()
    return bed_file

def processHeader(header_row):
    header_cols = [item.split(':')[0] for item in header_row]
    header_vals = [item.split(':')[1] for item in header_row]
    
    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    return header_data

def processBED(bed_rows):
    if len(bed_rows) > 0:
        bed_df = pd.DataFrame(bed_rows, columns=['chrom','start','end'])
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)
    else:
        bed_df = pd.DataFrame(columns = ['chrom','start','end'])

    return bed_df

def processSNPs(snp_rows):
    if len(snp_rows) > 0:
        snp_df = pd.DataFrame(snp_rows, columns= ["Ref_Loc","Cat","Ref_Base","Query_Base","Query_Loc","Ref_Aligned","Perc_Iden"])
        
        het_snps = snp_df[snp_df['Ref_Aligned'] == "Multiple"]
        non_het_snps = snp_df[snp_df['Ref_Aligned'] != "Multiple"]        
        non_het_snps['Ref_Aligned'] = non_het_snps['Ref_Aligned'].astype(int)

        snp_df = non_het_snps.append(het_snps)[["Ref_Loc","Cat","Ref_Base","Query_Base","Query_Loc","Ref_Aligned","Perc_Iden"]]

    else:
        snp_df = pd.DataFrame(columns = ["Ref_Loc","Cat","Ref_Base","Query_Base","Query_Loc","Ref_Aligned","Perc_Iden"])

    return snp_df

def processSNPDiffs(snpdiffs_file):

    with open(snpdiffs_file, 'r') as file:
        lines = file.readlines()

    bed_rows = []
    snp_rows = []

    line_count = 0
    for line in lines:
        line_count = line_count + 1
        if line[0:2] == "#\t":
            header_row = line.strip().split("\t")[1:]
        elif line[0:3] == "##\t":
            bed_rows.append(line.strip().split("\t")[1:])
        else:
            snp_rows.append(line.strip().split("\t"))

    alignment_data = processHeader(header_row)
    query_id = alignment_data['Query_ID'].iloc[0]
    bed_df = processBED(bed_rows)
    snp_df = processSNPs(snp_rows)
    bed_df.insert(0, 'Query', query_id)
    snp_df.insert(0, 'Query', query_id)
    return (alignment_data,bed_df,snp_df)

def returnUncovered(query_id,bed_df,yenta_bed):
    
    query_bed = BedTool.from_dataframe(bed_df[bed_df['Query'] == query_id][['chrom','start','end']])
    not_in_ref = yenta_bed.subtract(query_bed, A=True, header=False)
    if len(not_in_ref) > 0:
        not_in_ref = not_in_ref.to_dataframe()
        uncovered_list = not_in_ref['name'].tolist()
        uncovered_count = len(uncovered_list)
        return pd.DataFrame({'Query': [query_id] * uncovered_count,'Ref_Loc': uncovered_list,'Base': ["N"] * uncovered_count,'Cat':['Purged_Uncovered']*uncovered_count})
    else:
        return pd.DataFrame(columns = ['Query','Ref_Loc','Base','Cat'])

full_start_time = time.time()

# Get reference isolate
ref_isolate = str(sys.argv[1])

# Get output directory
output_directory = os.path.normpath(os.path.abspath(sys.argv[3]))
ref_directory = f"{output_directory}/SNP_{ref_isolate}"

##### Create logfile #####
global log_file
log_file = f"{output_directory}/logs/SNP_{ref_isolate}.log"

with open(log_file,"w+") as log:
    log.write("-------------------------------------------------------\n")
    log.write("SNP Analysis\n")
    log.write(f"Reference Isolate: {ref_isolate}\n")
    log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
    log.write("-------------------------------------------------------\n\n")

##### Read in snpdiffs files #####
with open(log_file,"a+") as log:
    log.write("Step 1: Reading in .snpdiffs files...")

# Read in list of .snpdiffs
snp_diffs_file = sys.argv[2]
if not os.path.exists(snp_diffs_file):
    with open(log_file,"a+") as log:
        log.write(f"\n\t- {snp_diffs_file} does not exist.")
    sys.exit(f"{snp_diffs_file} does not exist.")

# Ensure files exist
with open(snp_diffs_file, 'r') as file:
    snp_diffs_files = file.readlines()

snp_diffs_files = [line.strip() for line in snp_diffs_files]
if not all(os.path.exists(file_path) for file_path in snp_diffs_files):
    with open(log_file,"a+") as log:
        log.write(f"\n\t- One or more files provided by {snp_diffs_file} did not exist")
    sys.exit(f"One or more files provided by {snp_diffs_file} did not exist")

with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Found {str(len(snp_diffs_files))} .snpdiffs files\n")

snpdiffs_start = time.time()
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(processSNPDiffs,snp_diff_file) for snp_diff_file in snp_diffs_files]
results = [future.result() for future in concurrent.futures.as_completed(results)]

# Unpack the results
alignment_info_dfs, bed_dfs, snp_dfs = zip(*results)

# Concatenate the DataFrames along a specific axis
alignment_info = pd.concat(alignment_info_dfs, axis=0)  # Assuming axis=0 for concatenating rows
bed_df = pd.concat(bed_dfs, axis=0)
snp_df = pd.concat(snp_dfs, axis=0)

snpdiffs_end = time.time()
snp_diffs_time = snpdiffs_end - snpdiffs_start
with open(log_file,"a+") as log:
    log.write(f"\t- Processed {str(len(snp_diffs_files))} .snpdiffs files in {snp_diffs_time:.2f}s\n")
    log.write("\n-------------------------------------------------------\n\n")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 2: Performing sanity checks...")

if not alignment_info['Reference_ID'].nunique() == 1:
    sys.exit("The snpdiffs files provided contain data from more than one reference")

if not alignment_info['Reference_SHA256'].nunique() == 1:
    sys.exit("The snpdiffs files provided contain data from more than one reference")

if not alignment_info.groupby('Query_ID')['Query_SHA256'].nunique().eq(1).all():
    sys.exit("There is a query ID associated with 2+ unique SHA hashes")

if not alignment_info.groupby('Query_SHA256')['Query_ID'].nunique().eq(1).all():
    sys.exit("There is a query SHA hash associated with 2+ unique query IDs")

if not alignment_info['Query_ID'].nunique() == len(alignment_info['Query_ID']):
    sys.exit("The snpdiffs files provided contain data from some queries more than once")

if not alignment_info['QC_String'].nunique() == 1:
    sys.exit("The snpdiffs files provided were generated using different QC parameters")

query_data = alignment_info.drop_duplicates('Query_ID')[['Query_ID','Query_Assembly','Query_Contig_Count','Query_Assembly_Bases','Query_SHA256','Category','Percent_Reference_Aligned','Percent_Query_Aligned']]
reference_data = alignment_info.drop_duplicates('Reference_ID')[['Reference_ID','Reference_Assembly','Reference_Contig_Count','Reference_Assembly_Bases','Reference_SHA256']]
reference_data['Category'] = "Reference_Isolate"
reference_data['Percent_Reference_Aligned'] = "NA"
reference_data['Percent_Query_Aligned'] = "NA"
query_data.columns = ['Isolate_ID','Assembly','Contig_Count','Assembly_Bases','SHA256','Category','Percent_Reference_Aligned','Percent_Query_Aligned']
reference_data.columns = ['Isolate_ID','Assembly','Contig_Count','Assembly_Bases','SHA256','Category','Percent_Reference_Aligned','Percent_Query_Aligned']
isolate_data = pd.concat([query_data,reference_data]).reset_index(drop=True)

with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write("\n-------------------------------------------------------\n\n")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 3: Removing alignments that failed QC...")

alignments_pass_qc = alignment_info[alignment_info['Category'] == "PASS"]
alignments_fail_qc = alignment_info[alignment_info['Category'] != "PASS"]
snp_isolates = alignments_pass_qc['Query_ID'].astype(str).tolist()
query_isolates = [isolate for isolate in snp_isolates if not isolate == ref_isolate]
queries_fail_qc = alignments_fail_qc['Query_ID'].astype(str).tolist()

with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- {str(len(snp_isolates))} alignments passed QC\n")
    if len(queries_fail_qc) > 0:
        log.write("\t- The following isolate(s) are below the QC thresholds:\n")
        for fail_query in queries_fail_qc:
            log.write("\t- "+fail_query+"\n")
    if len(snp_isolates) == 0:
        log.write("\t- No comparisons remain after QC filtering. Cannot proceed.\n")
        sys.exit(0)
    else:    
        log.write("\n-------------------------------------------------------\n\n")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 4: Collecting SNPs from all queries...")

yenta_df = snp_df[snp_df.Cat == "SNP"][['Query','Ref_Loc','Query_Base','Ref_Base','Cat']].rename(columns={'Query_Base': 'Base'})
yenta_locs = np.unique(yenta_df['Ref_Loc'].values)
yenta_count = len(yenta_locs)

if yenta_count == 0:
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write("\n\t- No SNPs present in snpdiffs files. Cannot proceed as all SNP distances are 0...\n")
    sys.exit(0)
else:
    ref_isolate_df = yenta_df[['Ref_Loc', 'Ref_Base']].drop_duplicates()
    ref_isolate_df['Query'] = ref_isolate
    ref_isolate_df = ref_isolate_df[["Query","Ref_Loc","Ref_Base"]].rename(columns={"Ref_Base": "Base"})
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write(f"\t- {str(yenta_count)} SNPs detected\n")
        log.write("\n-------------------------------------------------------\n\n")


###########################################

with open(log_file,"a+") as log:
    log.write("Step 5: Processing purged loci...")

start_time = time.time()

# Remove any purged sites that have a valid SNP
purged_df = snp_df[(snp_df.Cat.str.startswith("Purged_")) & (snp_df['Ref_Loc'].isin(yenta_locs))][['Query','Ref_Loc','Cat']]
purged_df = purged_df.merge(yenta_df[['Query', 'Ref_Loc']], on=['Query', 'Ref_Loc'], how='left', indicator=True).query('_merge == "left_only"').drop(columns=['_merge'])

if purged_df.shape[0] > 0:

    purged_df = purged_df.drop_duplicates(subset=['Query','Ref_Loc'])
    purged_df['Base'] = "N"
    purged_df = purged_df[['Query','Ref_Loc','Base','Cat']]

    end_time = time.time()
    purge_time = end_time - start_time
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write("\t- Processed " + str(purged_df.shape[0]) + f" sites that weren't covered by valid alignments or failed density filters in the first MUMmer stage ({purge_time:.2f}s)...\n")
        log.write("\n-------------------------------------------------------\n\n")

else:
    purged_df = pd.DataFrame(columns=['Query','Ref_Loc','Base','Cat'])
    with open(log_file,"a+") as log:
        log.write("Done!")
        log.write("\t- No purged sites to process...\n")
        log.write("\n-------------------------------------------------------\n\n")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 6: Creating BED file for all SNP loci...") 

try:
    # Create BED file for Yenta locs
    start_time = time.time()
    yenta_bed = makeBED(pd.DataFrame([item.split('/') for item in yenta_locs], columns=['Ref_Contig','Ref_End']))
    end_time = time.time()
    bed_time = end_time - start_time
    with open(log_file,"a+") as log:
        log.write(f"Done! ({bed_time:.2f}s) \n")
        log.write("\n-------------------------------------------------------\n\n")
except:
    with open(log_file,"a+") as log:
        log.write("\n\t- Error: Cannot create Yenta BED file...\n")
    sys.exit("Cannot create Yenta BED file")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 7: Getting SNP loci that are not covered for each query...") 

try:
    start_time = time.time()
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(returnUncovered,query_id,bed_df,yenta_bed) for query_id in query_isolates]
    results = [future.result() for future in concurrent.futures.as_completed(results)]
    uncovered_df = pd.concat(results, axis=0)

    end_time = time.time()
    coverage_time = end_time - start_time
    
    with open(log_file,"a+") as log:
        log.write(f"Done! ({coverage_time:.2f}s) \n")
        if uncovered_df.shape[0] > 0:
            log.write(f"\t- Queries lacked alignment coverage for {str(uncovered_df.shape[0])} total sites\n")
        else:
            log.write("\t- No SNPs were uncovered in any query\n")
        log.write("\n-------------------------------------------------------\n\n")
except:
    with open(log_file,"a+") as log:
        log.write("\n\t- Error: Cannot get uncovered loci...\n")
    sys.exit("Cannot get uncovered loci")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 8: Getting SNP loci that have the reference base in another query...") 

try:
    start_time = time.time()
    dummy_snp_df = pd.DataFrame(list(product(query_isolates, yenta_locs)), columns=["Query", "Ref_Loc"])
    handled_loc_df = pd.concat([uncovered_df[["Query", "Ref_Loc"]], purged_df[["Query", "Ref_Loc"]], yenta_df[["Query", "Ref_Loc"]]])
    ref_base_df = dummy_snp_df.merge(handled_loc_df[['Query', 'Ref_Loc']], on=['Query', 'Ref_Loc'], how='left', indicator=True).query('_merge == "left_only"').drop(columns=['_merge']).merge(yenta_df[['Ref_Loc', 'Ref_Base']].rename(columns={'Ref_Base': 'Base'}), on=['Ref_Loc'], how='left').drop_duplicates(subset=['Query','Ref_Loc'])
    ref_base_df['Cat'] = "Ref_Base"

    end_time = time.time()
    coverage_time = end_time - start_time
    
    with open(log_file,"a+") as log:
        log.write(f"Done! ({coverage_time:.2f}s) \n")
        if ref_base_df.shape[0] > 0:
            log.write(f"\t- Queries had the reference base for for {str(ref_base_df.shape[0])} total sites\n")
        else:
            log.write("\t- No SNPs had the reference base?\n")
        log.write("\n-------------------------------------------------------\n\n")
except:
    with open(log_file,"a+") as log:
        log.write("\n\t- Error: Cannot get loci with the ref base...\n")
    sys.exit("Cannot get loci with the ref base")


###########################################

with open(log_file,"a+") as log:
    log.write("Step 9: Compiling data and checking SNP counts...") 

start_time = time.time()

ref_base_df = ref_base_df[['Query','Ref_Loc','Base','Cat']]
yenta_df = yenta_df[['Query','Ref_Loc','Base','Cat']]
uncovered_df = uncovered_df[['Query','Ref_Loc','Base','Cat']]
purged_df = purged_df[['Query','Ref_Loc','Base','Cat']]

# Save category counts
all_snp_data = pd.concat([yenta_df,ref_base_df,purged_df,uncovered_df])
snp_categories = ["SNP","Ref_Base","Purged_Uncovered","Purged_Alignment","Purged_N","Purged_Indel","Purged_Het","Purged_Density"]
locus_category_df = all_snp_data.groupby('Ref_Loc')['Cat'].value_counts().unstack(fill_value=0).reindex(columns=snp_categories, fill_value=0).reset_index()
locus_category_df.to_csv(ref_directory+"/Locus_Categories.tsv",sep="\t",index=False)
locus_category_df['Total_Covered'] = locus_category_df.iloc[:, 1:3].sum(axis=1)
locus_category_df['Total_Missing'] = locus_category_df.iloc[:, 3:9].sum(axis=1)

to_purge = locus_category_df[locus_category_df['Total_Missing'] / (locus_category_df['Total_Missing'] + locus_category_df['Total_Covered']) >= 2/3]
prune_locs = to_purge['Ref_Loc'].to_list()
purge_count = to_purge.shape[0]

query_category_df = all_snp_data.groupby('Query')['Cat'].value_counts().unstack(fill_value=0).reindex(columns=snp_categories, fill_value=0).reset_index().rename(columns={"Query":"Isolate_ID"})
isolate_data = pd.merge(isolate_data,query_category_df,on='Isolate_ID',how="left").fillna("NA")
int_cols = ['SNP','Ref_Base','Purged_Uncovered','Purged_Alignment','Purged_N','Purged_Indel','Purged_Het','Purged_Density']
isolate_data[int_cols] = isolate_data[int_cols].applymap(lambda x: f'{x:.0f}' if isinstance(x, (float, int)) and x.is_integer() else x)
isolate_data.to_csv(ref_directory+"/Isolate_Data.tsv",sep="\t",index=False)

with open(log_file,"a+") as log:
    log.write("Done!\n") 
    log.write(f"\t- Saved locus category data to {ref_directory}/Locus_Categories.tsv\n")
    log.write(f"\t- Saved isolate data to {ref_directory}/Isolate_Data.tsv\n")
    log.write("\n-------------------------------------------------------\n\n")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 10: Processing alignment data...")

start_time = time.time()

# Add reference data
all_snp_data = pd.concat([ref_isolate_df,all_snp_data[["Query","Ref_Loc","Base"]]]).sort_values(by=['Query', 'Ref_Loc'])
pruned_snp_data = all_snp_data[~all_snp_data['Ref_Loc'].isin(prune_locs)]

global base_df
base_df = all_snp_data[np.isin(all_snp_data['Base'].values, ['A', 'C', 'T', 'G','a','c','t','g'])]
base_df = base_df.set_index('Query')

# Pivot dataframe
pivoted_df = all_snp_data.pivot(index='Query', columns='Ref_Loc', values='Base')
yenta_ordered = pivoted_df.columns

with open(ref_directory+"/Loc_List.txt","w+") as snp_file:
    snp_file.write("\n".join(yenta_ordered)+"\n")

seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in pivoted_df.iterrows()]
alignment = MultipleSeqAlignment(seq_records)
AlignIO.write(alignment, ref_directory+"/snp_alignment.fasta","fasta")

if pruned_snp_data.shape[0] > 0:
    
    pruned_pivoted_df = pruned_snp_data.pivot(index='Query', columns='Ref_Loc', values='Base')
    pruned_yenta_ordered = pruned_pivoted_df.columns

    global pruned_base_df
    pruned_base_df = pruned_snp_data[np.isin(pruned_snp_data['Base'].values, ['A', 'C', 'T', 'G','a','c','t','g'])]
    pruned_base_df = pruned_base_df.set_index('Query')

    with open(ref_directory+"/Pruned_Loc_List.txt","w+") as snp_file:
        snp_file.write("\n".join(pruned_yenta_ordered)+"\n")

    pruned_seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in pruned_pivoted_df.iterrows()]
    pruned_alignment = MultipleSeqAlignment(pruned_seq_records)
    AlignIO.write(pruned_alignment, ref_directory+"/pruned_snp_alignment.fasta","fasta")

end_time = time.time()
align_time = end_time - start_time
with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Created alignment(s) in {align_time:.2f}s\n")
    log.write(f"\t- Saved alignment to "+ref_directory+"/snp_alignment.fasta\n")
    log.write(f"\t- Saved ordered list of loci to {ref_directory}/Loc_List.txt\n")
    if pruned_snp_data.shape[0] > 0:
        log.write(f"\t- Identified {str(purge_count)} sites where more isolates had missing data than had sequence data\n")
        log.write(f"\t- Saved pruned alignment to "+ref_directory+"/pruned_snp_alignment.fasta\n")
        log.write(f"\t- Saved ordered list of pruned loci to {ref_directory}/Pruned_Loc_List.txt\n")

    log.write("\n-------------------------------------------------------\n\n")

###########################################

with open(log_file,"a+") as log:
    log.write("Step 11: Processing pairwise comparisons files...")

start_time = time.time()
pairwise_combinations = [";".join(sorted(x)) for x in list(combinations(snp_isolates, 2))]

pairwise_dict = dict()
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(getPairwise,combo) for combo in pairwise_combinations]

for result in results:
    pairwise_dict.update(result.result())
pairwise_records = [{'Comparison': key, 'Cocalled_Sites': loc[0], 'SNP_Differences':loc[1]} for key, loc in pairwise_dict.items()]


pairwise_df = pd.DataFrame(pairwise_records)
print(pairwise_df)
pairwise_df[['Query_1', 'Query_2']] = pairwise_df['Comparison'].str.split(';', 1, expand=True)
pairwise_df[['Query_1','Query_2','SNP_Differences','Cocalled_Sites']].to_csv(ref_directory+"/distance_pairwise.tsv",sep="\t",index=False)


# Create matrix
idx = sorted(set(pairwise_df['Query_1']).union(pairwise_df['Query_2']))
mirrored_distance_df = pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Differences').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
mirrored_distance_df.index.name = ''
mirrored_distance_df.to_csv(ref_directory+"/distance_matrix.tsv",sep="\t")

if pruned_snp_data.shape[0] > 0:
    pruned_pairwise_dict = dict()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(getPrunedPairwise,combo) for combo in pairwise_combinations]

    for result in results:
        pruned_pairwise_dict.update(result.result())
    pruned_pairwise_records = [{'Comparison': key, 'Cocalled_Sites': loc[0], 'SNP_Differences':loc[1]} for key, loc in pruned_pairwise_dict.items()]

    pruned_pairwise_df = pd.DataFrame(pruned_pairwise_records)
    pruned_pairwise_df[['Query_1', 'Query_2']] = pruned_pairwise_df['Comparison'].str.split(';', 1, expand=True)
    pruned_pairwise_df[['Query_1','Query_2','SNP_Differences','Cocalled_Sites']].to_csv(ref_directory+"/pruned_distance_pairwise.tsv",sep="\t",index=False)

    # Create matrix
    pruned_idx = sorted(set(pruned_pairwise_df['Query_1']).union(pruned_pairwise_df['Query_2']))
    pruned_mirrored_distance_df = pruned_pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Differences').reindex(index=pruned_idx, columns=pruned_idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
    pruned_mirrored_distance_df.index.name = ''
    pruned_mirrored_distance_df.to_csv(ref_directory+"/pruned_distance_matrix.tsv",sep="\t")

end_time = time.time()
pairwise_time = end_time - start_time

with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Calculated "+str(len(pairwise_combinations)) +f" pairwise distances in {pairwise_time:.2f}s\n")
    log.write(f"\t- Saved pairwise distance data to {ref_directory}/distance_pairwise.tsv\n")
    log.write(f"\t- Saved distance matrix to {ref_directory}/distance_matrix.tsv\n")
    if pruned_snp_data.shape[0] > 0:
        log.write(f"\t- Saved pruned pairwise distance data to {ref_directory}/pruned_distance_pairwise.tsv\n")
        log.write(f"\t- Saved pruned distance matrix to {ref_directory}/pruned_distance_matrix.tsv\n")
    log.write("\n-------------------------------------------------------\n\n")

end_time = time.time()
total_runtime = end_time - full_start_time
with open(log_file,"a+") as log:
    log.write(f"Total Runtime: {total_runtime:.2f}s\n")






























"""
    # Add reference data
    ref_isolate_df = raw_snp_df[raw_snp_df['Ref_Loc'].isin(yenta_locs)][['Ref', 'Ref_Loc', 'Ref_Base']].drop_duplicates().rename(columns={"Ref_Base": "Base", "Ref": "Query"})
    
    big_df = pd.concat([big_df, ref_isolate_df]).sort_values(by=['Query', 'Ref_Loc'])
    
    
    pivoted_df = big_df.pivot(index='Query', columns='Ref_Loc', values='Base')
    yenta_ordered = pivoted_df.columns
    
    with open(ref_directory+"/ordered_loc_list.txt", "w") as file:
        file.write("\n".join(yenta_ordered) + "\n")

    loc_strings = {column:''.join(pivoted_df[column]) for column in pivoted_df.columns}
    isolate_strings = {index:''.join(row) for index,row in pivoted_df.iterrows()}

    loc_string_df = pd.concat([get_char_tally(key, val) for key, val in loc_strings.items()], ignore_index=True)
    isolate_string_df = pd.concat([get_char_tally(key, val) for key, val in isolate_strings.items()], ignore_index=True)
    
    loc_string_df.to_csv(ref_directory+"/locus_coverage_stats.tsv",sep="\t",index=False)
    isolate_string_df.to_csv(ref_directory+"/isolate_coverage_stats.tsv",sep="\t",index=False)

    global base_df
    base_df = big_df[np.isin(big_df['Base'].values, ['A', 'C', 'T', 'G'])]
    base_df.to_csv(ref_directory+"/covered_loci.tsv",header=False,sep="\t",index=False)
    end_time = time.time()
    processing_time = end_time - start_time
    
    with open(log_file,"a+") as log:
        log.write(f"\t- Processed data into single dataframe in {processing_time:.2f}s\n")
        log.write("\t- Saved base data for covered locs to "+ref_directory+"/covered_loci.tsv\n")
        log.write("\t- Saved coverage data for loci and isolates to "+ref_directory+"/locus_coverage_stats.tsv and "+ref_directory+"/isolate_coverage_stats.tsv\n")
        log.write("\n-------------------------------------------------------\n\n")

except:
    with open(log_file,"a+") as log:
        log.write("\n\t- ERROR: Cannot load coords/snp files in "+mummer_dir+"\n")
    sys.exit("Cannot load coords/snp files in "+mummer_dir)

##### Process alignment data #####
with open(log_file,"a+") as log:
    log.write("Step 6: Processing alignment data...")

start_time = time.time()

seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in pivoted_df.iterrows()]
alignment = MultipleSeqAlignment(seq_records)
AlignIO.write(alignment, ref_directory+"/snp_alignment.fasta","fasta")

mask_dict = {'n':'N','N':'N','-':'N','?':'N','A':'A','C':'C','T':'T','G':'G'}
masked_seq_records = [SeqRecord(Seq(''.join(row).translate(str.maketrans(mask_dict))), id=query,description='') for query,row in pivoted_df.iterrows()]
masked_alignment = MultipleSeqAlignment(masked_seq_records)
AlignIO.write(masked_alignment, ref_directory+"/masked_snp_alignment.fasta","fasta")

end_time = time.time()
align_time = end_time - start_time
with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Created alignment in {align_time:.2f}s\n")
    log.write("\t- Saved alignment to "+ref_directory+"/snp_alignment.fasta\n")
    log.write("\t- Saved masked alignment to "+ref_directory+"/masked_snp_alignment.fasta\n")
    log.write("\t- Saved ordered locus list for alignment to "+ref_directory+"/ordered_loc_list.txt\n")
    log.write("\n-------------------------------------------------------\n\n")

# Save pairwise comparisons
with open(log_file,"a+") as log:
    log.write("Step 7: Processing pairwise comparisons files...")

start_time = time.time()
pairwise_combinations = [";".join(sorted(x)) for x in list(combinations(snp_isolates, 2))]
base_df = base_df.set_index('Query')

pairwise_dict = dict()
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(getPairwise,combo) for combo in pairwise_combinations]
for result in results:
    pairwise_dict.update(result.result())
pairwise_records = [{'Comparison': key, 'Cocalled_Sites': loc[0], 'SNP_Differences':loc[1]} for key, loc in pairwise_dict.items()]

pairwise_df = pd.DataFrame(pairwise_records)
pairwise_df[['Query_1', 'Query_2']] = pairwise_df['Comparison'].str.split(';', 1, expand=True)
pairwise_df[['Query_1','Query_2','Cocalled_Sites','SNP_Differences']].to_csv(ref_directory+"/distance_pairwise.tsv",sep="\t",index=False)

end_time = time.time()
pairwise_time = end_time - start_time
with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Calculated "+str(len(pairwise_combinations)) +f" pairwise distances in {pairwise_time:.2f}s\n")
    log.write("\t- Saved pairwise distance data to "+ref_directory+"/distance_pairwise.tsv\n")
    log.write("\n-------------------------------------------------------\n\n")

print(density_filter(yenta_locs,density_windows,max_snps))

end_time = time.time()
total_runtime = end_time - full_start_time
with open(log_file,"a+") as log:
    log.write(f"Total Runtime: {total_runtime:.2f}s\n")

"""
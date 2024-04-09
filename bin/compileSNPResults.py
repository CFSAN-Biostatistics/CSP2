#!/usr/bin/env python3

import sys
import os
import glob
import pandas as pd
from itertools import chain
import scipy.stats
import numpy as np
import datetime
import time

def getWarnings(df):

    df_measures = list(set(df['Measure']))
    warn_list = []
    if 'Preserved_Diff' in df_measures:
        for index,row in df.iterrows():
            if pd.isna(row['Zscore']):
                warn_list.append(np.nan)
            elif 2.5 <= row['Zscore'] < 3:
                    warn_list.append("Warning")
            elif row['Zscore'] >=3:
                warn_list.append("Failure")
            else:
                warn_list.append("Pass")
    elif 'Contig_Count' in df_measures:
        for index,row in df.iterrows():
            if pd.isna(row['Zscore']):
                warn_list.append(np.nan)
            elif row['Measure'] in ["Contig_Count","L50","L90"]:
                if 2.5 <= row['Zscore'] < 3:
                    warn_list.append("Warning")
                elif row['Zscore'] >=3:
                    warn_list.append("Failure")
                else:
                    warn_list.append("Pass")
            elif row['Measure'] == "Assembly_Bases":
                if 2.5 <= abs(row['Zscore']) < 3:
                    warn_list.append("Warning")
                elif abs(row['Zscore']) >=3:
                    warn_list.append("Failure")
                else:
                    warn_list.append("Pass")     
            elif row['Measure'] in ["N50","N90"]:
                if -3 < row['Zscore'] <= -2.5:
                    warn_list.append("Warning")
                elif row['Zscore'] <= -3:
                    warn_list.append("Failure")
                else:
                    warn_list.append("Pass")
            else:
                sys.exit(f"{row['Measure']}")
    elif ('Raw_Distance_StdDev' in df_measures) | ('Preserved_Distance_StdDev' in df_measures):
        for index,row in df.iterrows():
            if pd.isna(row['Zscore']):
                warn_list.append(np.nan)
            elif 2.5 <= row['Zscore'] < 3:
                warn_list.append("Warning")
            elif row['Zscore'] >=3:
                warn_list.append("Failure")
            else:
                warn_list.append("Pass")  
    elif 'Unique_Kmers' in df_measures:
        for index,row in df.iterrows():
            if pd.isna(row['Zscore']):
                warn_list.append(np.nan)
            elif row['Measure'] in ["Align_Percent_Diff","Unique_Kmers","gIndels","Missing_Kmers"]:
                if 2.5 <= row['Zscore'] < 3:
                    warn_list.append("Warning")
                elif row['Zscore'] >=3:
                    warn_list.append("Failure")
                else:
                    warn_list.append("Pass")   
            elif row['Measure'] in ["Compare_Aligned","Kmer_Similarity","Self_Aligned","Median_Alignment_Length"]:
                if -3 < row['Zscore'] <= -2.5:
                    warn_list.append("Warning")
                elif row['Zscore'] <= -3:
                    warn_list.append("Failure")
                else:
                    warn_list.append("Pass")
            else:
                sys.exit(f"{row['Measure']}")
    
    elif ('SNPs_Cocalled' in df_measures) | ('Raw_SNPs_Cocalled' in df_measures) | ('Preserved_SNPs_Cocalled' in df_measures):
        for index,row in df.iterrows():
            if pd.isna(row['Zscore']):
                warn_list.append(np.nan)
            elif -3 < row['Zscore'] <= -2.5:
                    warn_list.append("Warning")
            elif row['Zscore'] <= -3:
                    warn_list.append("Failure")
            else:
                warn_list.append("Pass")
    else:
        sys.exit(f"{df_measures}")
    
    return warn_list

# Read in arguments
start_time = time.time()
snp_dirs = [line.strip() for line in open(sys.argv[1], 'r')]
raw_snp_distance_files = list(chain.from_iterable([glob.glob(snp_dir + '/*snp_distance_pairwise.tsv') for snp_dir in snp_dirs]))
screening_files = list(chain.from_iterable([glob.glob(snp_dir + '/*Reference_Screening.tsv') for snp_dir in snp_dirs]))

output_directory = sys.argv[2]
log_file = f"{output_directory}/Compilation.log"
mean_isolate_file = f"{output_directory}/Mean_Assembly_Stats.tsv"
isolate_assembly_stats_file = f"{output_directory}/Isolate_Assembly_Stats.tsv"
align_stats_file =  f"{output_directory}/Isolate_Alignment_Stats.tsv"
ref_mean_summary_file = f"{output_directory}/Align_Summary_by_Reference.tsv"
snp_comparison_file = f"{output_directory}/SNP_Distance_Summary.tsv"
qc_file = f"{output_directory}/QC_Warnings_Failures.tsv"

with open(log_file,"w+") as log:
    log.write("CSP2 SNP Pipeline Compiler\n")
    log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
    log.write("-------------------------------------------------------\n\n")
    if len(raw_snp_distance_files) == 0:
        log.write("\t- CSP2 SNP Pipeline Compiler cannot detected any SNP pipeline output files\n")
        log.write("\t- Compiler stopping...\n")
        sys.exit(0)
    
isolate_data = pd.read_csv(sys.argv[3],sep="\t")
raw_isolate_count = isolate_data.shape[0]

mummer_data = pd.read_csv(sys.argv[4],sep="\t")

# Get reference IDs
reference_ids = list(set(isolate_data[isolate_data['Isolate_Type'] == "Reference"]['Isolate_ID']))
raw_ref_count = len(reference_ids)

raw_snp_distance_df = pd.concat([pd.read_csv(file, sep='\t').assign(Reference_ID=os.path.basename(os.path.dirname(file))) for file in raw_snp_distance_files])
raw_snp_distance_df['Comparison'] = raw_snp_distance_df.apply(lambda row: ';'.join(sorted([str(row['Query_1']), str(row['Query_2'])])), axis=1)

# Check for preserved data
preserved_snp_distance_files = list(chain.from_iterable([glob.glob(snp_dir + '/*snp_distance_pairwise_preserved.tsv') for snp_dir in snp_dirs]))
if len(preserved_snp_distance_files) == 0:
    has_preserved = False
else:
    has_preserved = True
    preserved_snp_distance_df = pd.concat([pd.read_csv(file, sep='\t').assign(Reference_ID=os.path.basename(os.path.dirname(file))) for file in preserved_snp_distance_files])
    preserved_snp_distance_df['Comparison'] = preserved_snp_distance_df.apply(lambda row: ';'.join(sorted([str(row['Query_1']), str(row['Query_2'])])), axis=1)
    
screening_df = pd.concat([pd.read_csv(file, sep='\t').assign(Reference_ID=os.path.basename(os.path.dirname(file))) for file in screening_files])

snp_isolates = list(set(raw_snp_distance_df['Query_1'].tolist() + raw_snp_distance_df['Query_2'].tolist()))

with open(log_file,"a+") as log:
    log.write(f"- Detected SNP distance data for {len(snp_isolates)} isolates out of {raw_isolate_count} total isolates analyzed\n")
    if len(snp_isolates) <= 2:
        log.write("\t- CSP2 SNP Pipeline Compiler cannot do much with 2 or fewer isolates\n")
        log.write("\t- Compiler stopping...\n")
        sys.exit(0)
    else:
        failed_combinations = screening_df.loc[screening_df['Screen_Category'] != "Pass"]
        failed_comparisons = []

        if failed_combinations.shape[0] > 0:
            reference_query_dict = failed_combinations.groupby('Reference_ID')['Query_ID'].apply(list).to_dict()
            for index,row in failed_combinations.iterrows():
                failed_comparisons.append(";".join(sorted([row['Query_ID'], row['Reference_ID']]))) 
            log.write("\n- The following query-reference alignments did not pass QC\n")
            for key, value in reference_query_dict.items():
                log.write(f"\nReference {key}\n{', '.join(map(str, value))}\n")

# Prune isolate data
isolate_data = isolate_data.loc[isolate_data['Isolate_ID'].isin(snp_isolates)]
reference_ids = [x for x in reference_ids if x in snp_isolates]
ref_count = len(reference_ids)
with open(log_file,"a+") as log:
    log.write(f"- Detected SNP distance data for {ref_count} reference isolates out of {raw_ref_count} total reference isolates analyzed\n")
    for ref in reference_ids:
        log.write(f"\t- {ref}\n")

# Prune MUMmer data
mummer_data['Comparison'] = mummer_data.apply(lambda row: ';'.join(sorted([str(row['Query_ID']), str(row['Reference_ID'])])), axis=1)
mummer_data = mummer_data.loc[~mummer_data['Comparison'].isin(failed_comparisons), ['SNPDiffs_File','Query_ID','Reference_ID','Comparison','Reference_Percent_Aligned','Query_Percent_Aligned','Median_Alignment_Length','Kmer_Similarity','Reference_Unique_Kmers','Query_Unique_Kmers','gIndels']]
max_align_values = np.maximum(mummer_data['Reference_Percent_Aligned'], mummer_data['Query_Percent_Aligned'])
min_align_values = np.minimum(mummer_data['Reference_Percent_Aligned'], mummer_data['Query_Percent_Aligned'])
mummer_data['Align_Percent_Diff'] = 100*((max_align_values - min_align_values)/min_align_values)

# Run basic assembly stats
isolate_stats = isolate_data.melt(id_vars=['Isolate_ID', 'Isolate_Type'], value_vars = ['Contig_Count','Assembly_Bases','N50','N90','L50','L90'],value_name='Value',var_name = "Measure")
isolate_stats['Zscore'] = isolate_stats.groupby('Measure')['Value'].transform(scipy.stats.zscore).astype('float').round(3)
isolate_stats['QC'] = getWarnings(isolate_stats)
isolate_stats['Value'] = isolate_stats['Value'].astype('int')

# Reformat for final TSV
isolate_stats['Min'] = np.nan
isolate_stats['Max'] = np.nan
isolate_stats['StdDev'] = np.nan
isolate_stats = isolate_stats[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC']].rename(columns = {'Value':'Mean'})
isolate_stats['Count'] = 1

# Get mean values
isolate_mean_df = isolate_stats.groupby(by=['Measure'])['Mean'].agg(Min = 'min',Mean = "mean",Max = 'max',StdDev = 'std',Count = 'count')
isolate_mean_df['Mean'] = isolate_mean_df['Mean'].astype("int")
isolate_mean_df['StdDev'] = isolate_mean_df['StdDev'].astype("float").round(3)
with open(log_file,"a+") as log:
    log.write("- Read in and processed isolate data\n\n")
    for index, row in isolate_mean_df.iterrows():
        log.write(f"{index}:\tMin: {row['Min']}\tMean: {row['Mean']}\tMax: {row['Max']}\tStdDev: {row['StdDev']}\n")
    log.write("\n")

# Run basic alignment stats
isolate_mummer = pd.DataFrame(columns=['Isolate_ID', 'Compare_ID', 'Self_Aligned', 'Compare_Aligned','Align_Percent_Diff','Median_Alignment_Length', 'Kmer_Similarity', 'Unique_Kmers', 'Missing_Kmers', 'gIndels', 'SNPDiffs_File'])

for isolate in snp_isolates:
    
    temp_mummer = mummer_data[(mummer_data['Query_ID'] == isolate) | (mummer_data['Reference_ID'] == isolate)].drop_duplicates(subset=['Comparison']).assign(Isolate_ID = isolate)
    
    for index, row in temp_mummer.iterrows():
        if row['Query_ID'] == isolate:
            temp_isolate_mummer = row[['Isolate_ID', 'Reference_ID', 'Query_Percent_Aligned', 'Reference_Percent_Aligned','Align_Percent_Diff', 'Median_Alignment_Length', 'Kmer_Similarity', 'Query_Unique_Kmers', 'Reference_Unique_Kmers', 'gIndels', 'SNPDiffs_File']].to_frame().T
            temp_isolate_mummer.columns = ['Isolate_ID', 'Compare_ID', 'Self_Aligned', 'Compare_Aligned', 'Align_Percent_Diff','Median_Alignment_Length', 'Kmer_Similarity', 'Unique_Kmers', 'Missing_Kmers', 'gIndels', 'SNPDiffs_File']
            isolate_mummer = pd.concat([isolate_mummer,temp_isolate_mummer])

        elif row['Reference_ID'] == isolate:
            temp_isolate_mummer = row[['Isolate_ID', 'Query_ID', 'Reference_Percent_Aligned', 'Query_Percent_Aligned', 'Align_Percent_Diff', 'Median_Alignment_Length', 'Kmer_Similarity', 'Reference_Unique_Kmers', 'Query_Unique_Kmers', 'gIndels', 'SNPDiffs_File']].to_frame().T
            temp_isolate_mummer.columns = ['Isolate_ID', 'Compare_ID', 'Self_Aligned', 'Compare_Aligned', 'Align_Percent_Diff', 'Median_Alignment_Length', 'Kmer_Similarity', 'Unique_Kmers', 'Missing_Kmers', 'gIndels', 'SNPDiffs_File']
            isolate_mummer = pd.concat([isolate_mummer,temp_isolate_mummer])

isolate_mummer['Isolate_Type'] = isolate_mummer['Isolate_ID'].apply(lambda x: 'Reference' if x in reference_ids else 'Query')
isolate_mummer_df = isolate_mummer.melt(id_vars=['Isolate_ID','Isolate_Type'], value_vars = ['Self_Aligned','Compare_Aligned', 'Align_Percent_Diff','Median_Alignment_Length','Kmer_Similarity','Unique_Kmers','Missing_Kmers','gIndels'],value_name='Value',var_name = "Measure")
isolate_mummer_df['Value'] = isolate_mummer_df['Value'].astype("float")
isolate_mummer_df = isolate_mummer_df.groupby(by=['Isolate_ID','Isolate_Type','Measure'])['Value'].agg(Count = "count",Min = "min",Value = "mean",Max = "max",StdDev = 'std').reset_index()
isolate_mummer_df['Value'] = isolate_mummer_df['Value'].astype("float").round(2)

# Get Zscores
isolate_mummer_df['Zscore'] =  isolate_mummer_df.groupby('Measure')['Value'].transform(scipy.stats.zscore).astype('float').round(3)
isolate_mummer_df['QC'] = getWarnings(isolate_mummer_df)

# Reformat for final TSV
align_stats = isolate_mummer_df[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count']].copy().rename(columns = {"Value":"Mean"})
align_stats['StdDev'] = align_stats['StdDev'].astype('float').round(3)

with open(log_file,"a+") as log:
    log.write("- Read in and processed alignment data\n")
    
# Process cocalled data
raw_cocalled_df = raw_snp_distance_df[['Comparison','Query_1','Query_2','Reference_ID','SNPs_Cocalled']]
isolate_cocalled_df = pd.DataFrame(columns = ['Isolate_ID','Count','Min','Mean','Max','StdDev'])

for isolate in snp_isolates:
    temp_cocalled = raw_cocalled_df[(raw_cocalled_df['Query_1'] == isolate) | (raw_cocalled_df['Query_2'] == isolate)].drop_duplicates(subset=['Comparison']).assign(Isolate_ID = isolate)
    temp_cocalled = temp_cocalled.groupby(['Isolate_ID'])['SNPs_Cocalled'].agg(Count = "count", Min = "min", Value = "mean", Max = "max",StdDev = 'std').reset_index()
    isolate_cocalled_df = pd.concat([isolate_cocalled_df,temp_cocalled])

isolate_cocalled_df['Measure'] = 'Raw_SNPs_Cocalled'
isolate_cocalled_df['Value'] = isolate_cocalled_df['Value'].astype('int')
isolate_cocalled_df['Zscore'] = isolate_cocalled_df['Value'].transform(scipy.stats.zscore).astype('float').round(3)
isolate_cocalled_df['QC'] = getWarnings(isolate_cocalled_df)

# Format for final TSV
isolate_cocalled_df['Isolate_Type'] = isolate_cocalled_df['Isolate_ID'].apply(lambda x: 'Reference' if x in reference_ids else 'Query')
isolate_cocalled_stats = isolate_cocalled_df[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count']].copy().rename(columns={'Value':'Mean'})

if has_preserved:
    preserved_cocalled_df = preserved_snp_distance_df[['Comparison','Query_1','Query_2','Reference_ID','SNPs_Cocalled']]
    isolate_preserved_cocalled_df = pd.DataFrame(columns = ['Isolate_ID','Count','Min','Mean','Max','StdDev'])

    for isolate in snp_isolates:
        temp_cocalled = preserved_cocalled_df[(preserved_cocalled_df['Query_1'] == isolate) | (preserved_cocalled_df['Query_2'] == isolate)].drop_duplicates(subset=['Comparison']).assign(Isolate_ID = isolate)
        temp_cocalled = temp_cocalled.groupby(['Isolate_ID'])['SNPs_Cocalled'].agg(Count = "count", Min = "min", Value = "mean", Max = "max",StdDev = 'std').reset_index()
        isolate_preserved_cocalled_df = pd.concat([isolate_preserved_cocalled_df,temp_cocalled])

    isolate_preserved_cocalled_df['Measure'] = 'Preserved_SNPs_Cocalled'
    isolate_preserved_cocalled_df['Value'] = isolate_preserved_cocalled_df['Value'].astype('int')
    isolate_preserved_cocalled_df['Zscore'] = isolate_preserved_cocalled_df['Value'].transform(scipy.stats.zscore).astype('float').round(3)
    isolate_preserved_cocalled_df['QC'] = getWarnings(isolate_preserved_cocalled_df)

    # Format for final TSV
    isolate_preserved_cocalled_df['Isolate_Type'] = isolate_preserved_cocalled_df['Isolate_ID'].apply(lambda x: 'Reference' if x in reference_ids else 'Query')
    isolate_preserved_cocalled_df['StdDev'] = isolate_preserved_cocalled_df['StdDev'].astype('float').round(3)
    isolate_cocalled_stats = pd.concat([isolate_cocalled_stats,isolate_preserved_cocalled_df[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count']].copy().rename(columns={'Value':'Mean'})])
    
with open(log_file,"a+") as log:
    log.write("- Processed cocalled SNP data\n")
    
if has_preserved:
    raw_snp_df = raw_snp_distance_df[['Comparison','Query_1','Query_2','Reference_ID','SNP_Distance']].rename(columns = {'SNP_Distance':'Raw_SNP_Distance'})
    preserved_snp_df = preserved_snp_distance_df[['Comparison','Query_1','Query_2','Reference_ID','SNP_Distance']].rename(columns = {'SNP_Distance':'Preserved_SNP_Distance'})
    snp_df = pd.merge(raw_snp_df,preserved_snp_df,how="left",on=['Comparison','Query_1','Query_2','Reference_ID'])
    snp_df['Preserved_Diff'] = abs(snp_df['Preserved_SNP_Distance'] - snp_df['Raw_SNP_Distance'])
    
    isolate_snp_df = pd.DataFrame(columns = ['Isolate_ID','Count','Min','Value','Max','StdDev'])

    for isolate in snp_isolates:
        temp_snp = snp_df[(snp_df['Query_1'] == isolate) | (snp_df['Query_2'] == isolate)].drop_duplicates(subset=['Comparison']).assign(Isolate_ID = isolate)
        temp_snp = temp_snp.groupby(['Isolate_ID'])['Preserved_Diff'].agg(Count = "count", Min = "min", Value = "mean", Max = "max",StdDev = 'std').reset_index()
        isolate_snp_df = pd.concat([isolate_snp_df,temp_snp])

    isolate_snp_df['Measure'] = 'Preserved_Diff'
    isolate_snp_df['Value'] = isolate_snp_df['Value'].astype("float")
    isolate_snp_df['Zscore'] = isolate_snp_df['Value'].transform(scipy.stats.zscore).astype('float').round(3)
    isolate_snp_df['Value'] = isolate_snp_df['Value'].astype('float').round(3)
    isolate_snp_df['QC'] = getWarnings(isolate_snp_df)

    # Format for final TSV
    isolate_snp_df['Isolate_Type'] = isolate_snp_df['Isolate_ID'].apply(lambda x: 'Reference' if x in reference_ids else 'Query')
    isolate_snp_df['StdDev'] = isolate_snp_df['StdDev'].astype('float').round(3)
    isolate_snp_stats = isolate_snp_df[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count']].copy().rename(columns={'Value':'Mean'})
    with open(log_file,"a+") as log:
        log.write("- Processed preserved SNP data\n")
else:
    with open(log_file,"a+") as log:
        log.write("- No preserved SNP data to process\n")
    isolate_snp_stats = pd.DataFrame(columns = ['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count'])

# Compare SNPs across refs
if len(reference_ids) == 1:
    isolate_stdev_stats = pd.DataFrame(columns =['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count'])
    with open(log_file,"a+") as log:
        log.write("- 1 reference provided, SNP distances have no comparisons\n")
else:
    # Get comparison stats
    raw_comparison_df = raw_snp_distance_df.groupby(by=['Comparison'])['SNP_Distance'].agg(Count = 'count', Min = 'min', Mean = 'mean', Max = 'max', StdDev = 'std').reset_index()
    raw_comparison_df['StdDev'] = raw_comparison_df['StdDev'].astype('float').round(3)
    raw_comparison_df['Mean'] = raw_comparison_df['Mean'].astype('int')
    raw_comparison_df[['Query_1', 'Query_2']] = raw_comparison_df['Comparison'].str.split(';', expand=True)
    raw_comparison_df['SNP_Spread'] = abs(raw_comparison_df['Max'] - raw_comparison_df['Min'])

    comparison_df = raw_comparison_df[['Comparison','Query_1','Query_2','Mean','StdDev','Min','Max','SNP_Spread','Count']].copy()
    
    # Get isolate stats
    isolate_stdev_df = pd.DataFrame(columns = ['Isolate_ID','Count','Min','Value','Max','StdDev'])

    for isolate in snp_isolates:
        temp_compare = raw_comparison_df[(raw_comparison_df['Query_1'] == isolate) | (raw_comparison_df['Query_2'] == isolate)].drop_duplicates(subset=['Comparison']).assign(Isolate_ID = isolate)
        temp_compare = temp_compare.groupby(by=['Isolate_ID'])['StdDev'].agg(Count = 'count', Min = 'min', Value = 'mean', Max = 'max', StdDev = 'std').reset_index()
        isolate_stdev_df = pd.concat([isolate_stdev_df,temp_compare])
    
    isolate_stdev_df['Measure'] = "Raw_Distance_StdDev"
    isolate_stdev_df['Value'] = isolate_stdev_df['Value'].astype("float")
    isolate_stdev_df['Zscore'] = isolate_stdev_df['Value'].transform(scipy.stats.zscore).astype('float').round(3)
    isolate_stdev_df['Value'] = isolate_stdev_df['Value'].astype('float').round(3)
    isolate_stdev_df['Isolate_Type'] = isolate_stdev_df['Isolate_ID'].apply(lambda x: 'Reference' if x in reference_ids else 'Query')
    isolate_stdev_df['QC'] = getWarnings(isolate_stdev_df)
    
    isolate_stdev_stats = isolate_stdev_df[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count']].copy().rename(columns={'Value':'Mean'})
        
    if has_preserved:
        comparison_df.columns = ['Comparison','Query_1','Query_2','Raw_Mean','Raw_StdDev','Raw_Min','Raw_Max','Raw_SNP_Spread','Raw_Count']

        preserved_comparison_df = preserved_snp_distance_df.groupby(by=['Comparison'])['SNP_Distance'].agg(Preserved_Count = 'count', Preserved_Min = 'min', Preserved_Mean = 'mean', Preserved_Max = 'max', Preserved_StdDev = 'std').reset_index()
        preserved_comparison_df['Preserved_StdDev'] = preserved_comparison_df['Preserved_StdDev'].astype('float').round(3)
        preserved_comparison_df['Preserved_Mean'] = preserved_comparison_df['Preserved_Mean'].astype('int')
        preserved_comparison_df[['Query_1', 'Query_2']] = preserved_comparison_df['Comparison'].str.split(';', expand=True)
        preserved_comparison_df['Preserved_SNP_Spread'] = abs(preserved_comparison_df['Preserved_Max'] - preserved_comparison_df['Preserved_Min'])
        
        comparison_df = comparison_df.merge(preserved_comparison_df,how = "left", on=['Comparison','Query_1','Query_2'])
        comparison_df['Mean_Preserved_Diff'] = abs(comparison_df['Preserved_Mean'] - comparison_df['Raw_Mean'])
        comparison_df = comparison_df[['Query_1','Query_2','Raw_Mean','Preserved_Mean','Mean_Preserved_Diff','Raw_StdDev','Preserved_StdDev','Raw_SNP_Spread','Preserved_SNP_Spread','Raw_Min','Raw_Max','Preserved_Min','Preserved_Max','Raw_Count','Preserved_Count']]
        
        isolate_stdev_df = pd.DataFrame(columns = ['Isolate_ID','Count','Min','Value','Max','StdDev'])

        for isolate in snp_isolates:
            temp_compare = preserved_comparison_df[(preserved_comparison_df['Query_1'] == isolate) | (preserved_comparison_df['Query_2'] == isolate)].drop_duplicates(subset=['Comparison']).assign(Isolate_ID = isolate)
            temp_compare = temp_compare.groupby(by=['Isolate_ID'])['Preserved_StdDev'].agg(Count = 'count', Min = 'min', Value = 'mean', Max = 'max', StdDev = 'std').reset_index()
            isolate_stdev_df = pd.concat([isolate_stdev_df,temp_compare])
    
        isolate_stdev_df['Measure'] = "Preserved_Distance_StdDev"
        isolate_stdev_df['Value'] = isolate_stdev_df['Value'].astype('float')
        isolate_stdev_df['Zscore'] = isolate_stdev_df['Value'].transform(scipy.stats.zscore).astype('float').round(3)
        isolate_stdev_df['Value'] = isolate_stdev_df['Value'].astype('float').round(3)
        isolate_stdev_df['Isolate_Type'] = isolate_stdev_df['Isolate_ID'].apply(lambda x: 'Reference' if x in reference_ids else 'Query')
        isolate_stdev_df['QC'] = getWarnings(isolate_stdev_df)
        isolate_stdev_stats = pd.concat([isolate_stdev_stats,isolate_stdev_df[['Isolate_ID','Isolate_Type','Measure','Min','Value','Max','StdDev','Zscore','QC','Count']].copy().rename(columns={'Value':'Mean'})])
        with open(log_file,"a+") as log:
            log.write("- Compared results across references\n")
    else:
        with open(log_file,"a+") as log:
            log.write("- Compared results across references\n")
# Group by ref

#### Isolate ####
ref_isolate_df = isolate_stats.loc[isolate_stats['Isolate_Type'] == "Reference"][['Isolate_ID','Measure','Mean','StdDev','Min','Max','Zscore','QC','Count']].rename(columns = {'Isolate_ID':'Reference_ID'})

#### StdDev ####
ref_stdev_df = isolate_stdev_stats.loc[isolate_stdev_stats['Isolate_Type'] == "Reference"][['Isolate_ID','Measure','Mean','StdDev','Min','Max','Zscore','QC','Count']].rename(columns = {'Isolate_ID':'Reference_ID'})

#### MUMmer ####
ref_mummer_df = pd.DataFrame(columns = ['Reference_ID','Measure','Mean','StdDev','Min','Max','Count'])
for ref in reference_ids:
    ref_mummer =  isolate_mummer[(isolate_mummer['Isolate_ID'] == ref) | (isolate_mummer['Compare_ID'] == ref)].assign(Focal_Reference = ref)
    ref_mummer['Comparison'] = ref_mummer.apply(lambda row: ';'.join(sorted([str(row['Isolate_ID']), str(row['Compare_ID'])])), axis=1)
    ref_mummer = ref_mummer.drop_duplicates(subset=['Comparison'])

    ref_mummer = ref_mummer.melt(id_vars=['Focal_Reference','Isolate_ID','Compare_ID'], value_vars = ['Align_Percent_Diff','Median_Alignment_Length','Kmer_Similarity','gIndels'],value_name='Value',var_name = "Measure")
    ref_mummer['Value'] = ref_mummer['Value'].astype("float")
    ref_mummer = ref_mummer.groupby(by=['Measure'])['Value'].agg(Count = "count",Min = "min",Mean = "mean",Max = "max",StdDev = 'std').reset_index().assign(Reference_ID = ref)
    ref_mummer = ref_mummer[['Reference_ID','Measure','Mean','StdDev','Min','Max','Count']]
    ref_mummer_df = pd.concat([ref_mummer_df,ref_mummer])
ref_mummer_df['QC'] = np.nan
ref_mummer_df['Zscore'] = np.nan

ref_mummer_summary_df = pd.concat([ref_mummer_df[['Reference_ID','Measure','Mean','StdDev','Min','Max','Zscore','QC','Count']],align_stats.loc[(align_stats['Isolate_Type'] == "Reference") & (align_stats['Measure'].isin(['Self_Aligned','Compare_Aligned','Unique_Kmers','Missing_Kmers']))][['Isolate_ID','Measure','Mean','StdDev','Min','Max','Zscore','QC','Count']].rename(columns = {'Isolate_ID':'Reference_ID'})])

#### Cocalled ####
ref_cocalled_summary_df = raw_cocalled_df.groupby(by=['Reference_ID'])['SNPs_Cocalled'].agg(Mean = "mean",StdDev = 'std',Min = "min",Max = "max",Count = 'count').reset_index()
ref_cocalled_summary_df['Measure'] = "Raw_SNPs_Cocalled"
ref_cocalled_summary_df['QC'] = np.nan
ref_cocalled_summary_df['Zscore'] = np.nan
ref_cocalled_summary_df = ref_cocalled_summary_df[['Reference_ID','Measure','Mean','StdDev','Min','Max','Zscore','QC','Count']]

if has_preserved:
    preserved_cocalled_summary = preserved_cocalled_df.groupby(by=['Reference_ID'])['SNPs_Cocalled'].agg(Mean = "mean",StdDev = 'std',Min = "min",Max = "max",Count = 'count').reset_index()
    preserved_cocalled_summary['Measure'] = "Preserved_SNPs_Cocalled"
    preserved_cocalled_summary['QC'] = np.nan
    preserved_cocalled_summary['Zscore'] = np.nan
    ref_cocalled_summary_df = pd.concat([ref_cocalled_summary_df,preserved_cocalled_summary[['Reference_ID','Measure','Mean','StdDev','Min','Max','Zscore','QC','Count']]])

#### Preserved Diff ####
if has_preserved:
    ref_summary_preserved_df = snp_df.groupby(by=['Reference_ID'])['Preserved_Diff'].agg(Mean = "mean",StdDev = 'std',Min = "min",Max = "max",Count = 'count').reset_index()
    ref_summary_preserved_df['Measure'] = "Preserved_Diff"
    ref_summary_preserved_df['QC'] = np.nan
    ref_summary_preserved_df['Zscore'] = np.nan
    ref_summary_preserved_df = ref_summary_preserved_df[['Reference_ID','Measure','Mean','StdDev','Min','Max','Count','Zscore','QC']].copy()
    
#### Compile ####
ref_summary_df = pd.concat([ref_mummer_summary_df,
                            ref_cocalled_summary_df,
                            ref_isolate_df,ref_stdev_df]).sort_values(by=['Measure'])

ref_summary_df['Mean'] = ref_summary_df['Mean'].astype("float").round(3)
ref_summary_df['Min'] = ref_summary_df['Min'].astype("float").round(3)
ref_summary_df['Max'] = ref_summary_df['Max'].astype("float").round(3)
ref_summary_df['StdDev'] = ref_summary_df['StdDev'].astype("float").round(3)

# Catch warnings and failures
all_isolate_stats = pd.concat([isolate_stats,align_stats,isolate_stdev_stats,isolate_cocalled_stats,isolate_snp_stats]).sort_values(by=['Zscore'])

warn_fail_df = all_isolate_stats.loc[all_isolate_stats['QC'].isin(['Failure','Warning'])].copy()
warn_fail_df['abs_Zscore'] = warn_fail_df['Zscore'].abs()
warn_fail_df = warn_fail_df.sort_values(by='abs_Zscore',ascending=False).drop('abs_Zscore',axis=1)

warn_fail_isolates = list(set(warn_fail_df['Isolate_ID']))
if len(warn_fail_isolates) > 0:
    with open(log_file,"a+") as log:
        log.write("\n- The following samples had QC warnings or failures:\n")
        for isolate in warn_fail_isolates:
            isolate_warn_fail = warn_fail_df.loc[warn_fail_df['Isolate_ID'] == isolate]
            
            if isolate in reference_ids:
                log.write(f"\n{isolate} (Reference):\n")
            else:
                log.write(f"\n{isolate} (Query):\n")
            
            for index,row in isolate_warn_fail.iterrows():
                log.write(f"\t- {row['Measure']} - Mean: {row['Mean']}; Zscore: {row['Zscore']}; QC: {row['QC']}\n")
else:
    with open(log_file,"a+") as log:
        log.write("-There were no QC warnings or failures\n")

# Output data

# Mean assembly stats
isolate_mean_df.reset_index().to_csv(mean_isolate_file,sep='\t',index=False)

# Isolate assembly stats
isolate_assembly_stats = isolate_stats.loc[isolate_stats['Measure'].isin(['Contig_Count','Assembly_Bases','L50','L90','N50','N90'])].drop(['Min','Max','StdDev','Count'],axis=1).rename(columns = {'Mean':'Value'})
isolate_assembly_stats.to_csv(isolate_assembly_stats_file,sep='\t',index=False)

# Isolate alignment stats
isolate_align_stats = pd.concat([align_stats,isolate_cocalled_stats,isolate_snp_stats,isolate_stdev_stats]).reset_index(drop=True)
for col in ['Min', 'Mean', 'Max', 'StdDev', 'Zscore']:
    isolate_align_stats[col] = isolate_align_stats[col].astype("float").round(3)
isolate_align_stats.to_csv(align_stats_file,sep='\t',index=False)

# Reference Assembly Stats
ref_align_summary_df = ref_summary_df.loc[(~ref_summary_df['Measure'].isin(['Contig_Count','Assembly_Bases','L50','L90','N50','N90'])) & (~pd.isna(ref_summary_df['Zscore']))]
ref_mean_summary_df = ref_summary_df.loc[(~ref_summary_df['Measure'].isin(['Contig_Count','Assembly_Bases','L50','L90','N50','N90'])) & (pd.isna(ref_summary_df['Zscore']))].drop(['Zscore','QC'],axis =1)
ref_mean_summary_df['Zscore'] = np.nan
ref_mean_summary_df['QC'] = np.nan

# Add alignment stats
if has_preserved:
    ref_mean_summary_df = pd.concat([ref_mean_summary_df,ref_summary_preserved_df])

ref_isolate_align_stats = align_stats.loc[(align_stats['Isolate_Type'] == "Reference") & (align_stats['Measure'].isin(['Self_Aligned','Compare_Aligned']))].drop(['Isolate_Type'],axis=1).rename(columns = {'Isolate_ID':'Reference_ID'})[['Reference_ID','Measure','Mean','StdDev','Min','Max','Count','Zscore','QC']]

ref_mean_summary_stats = pd.concat([ref_mean_summary_df,ref_isolate_align_stats])
ref_mean_summary_stats.to_csv(ref_mean_summary_file,sep='\t',index=False)

end_time = time.time()

with open(log_file,"a+") as log:
    log.write(f"\n- Completed compilation in {end_time - start_time:.2f} seconds\n")
    log.write(f"\t- Saved mean isolate assembly data to {mean_isolate_file}\n")
    log.write(f"\t- Saved raw isolate assembly data to {isolate_assembly_stats_file}\n")
    log.write(f"\t- Saved isolate alignment data to {align_stats_file}\n")
    log.write(f"\t- Saved reference summary data to {ref_mean_summary_file}\n")
    
    # Comparisons if multiple refs
    if len(reference_ids) > 1:
        comparison_df.to_csv(snp_comparison_file,sep="\t",index = False)
        log.write(f"\t- Saved SNP distance comparisons across references to {snp_comparison_file}\n")
    
    # Failures/warnings
    if warn_fail_df.shape[0] > 0:
        warn_fail_df.to_csv(qc_file,sep="\t",index=False)
        log.write(f"\t- Saved QC warnings/failures to {qc_file}\n")
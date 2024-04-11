#!/usr/bin/env python3

import sys
import os
import pandas as pd
import datetime
from pybedtools import BedTool,helpers
import concurrent.futures
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from itertools import combinations
import numpy as np
import uuid
import traceback

def fetchHeaders(snpdiffs_file):
    
    with open(snpdiffs_file, 'r') as file:
        top_line = file.readline().strip().split('\t')[1:]

    header_cols = [item.split(':')[0] for item in top_line]
    header_vals = [item.split(':')[1] for item in top_line]
    
    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    header_data.loc[:, 'File_Path'] = snpdiffs_file
      
    return header_data

def processBED(bed_rows,snpdiffs_orientation):
    
    bed_columns = ['Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned',
                   'Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned',
                   'Perc_Iden']
        
    reverse_columns = ['Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned',
                     'Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned',
                     'Perc_Iden']
    
    int_columns = ['Ref_Start', 'Ref_End', 'Ref_Length', 'Ref_Aligned',
                'Query_Start', 'Query_End', 'Query_Length', 'Query_Aligned']

    float_columns = ['Perc_Iden']
    
    if len(bed_rows) > 0:
        
        bed_df = pd.DataFrame(bed_rows, columns=bed_columns)
        
        # Swap columns if reversed
        if snpdiffs_orientation == -1:
            bed_df = bed_df[reverse_columns].copy()
            bed_df.columns = bed_columns

        # Gather covered loci
        covered_bed_df = bed_df[(bed_df['Ref_Start'] != ".") & (bed_df['Query_Start'] != ".")].copy()
        if covered_bed_df.shape[0] > 0:
            for col in int_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(int)
            for col in float_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float)
        
        return covered_bed_df
    else:
        return pd.DataFrame(columns=bed_columns)

def processSNPs(snp_rows,snpdiffs_orientation):
    
    snp_columns = ['Ref_Contig','Start_Ref','Ref_Pos',
                'Query_Contig','Start_Query','Query_Pos',
                'Ref_Loc','Query_Loc',
                'Ref_Start','Ref_End',
                'Query_Start','Query_End',
                'Ref_Base','Query_Base',
                'Dist_to_Ref_End','Dist_to_Query_End',
                'Ref_Aligned','Query_Aligned',
                'Query_Direction','Perc_Iden','Cat']
    
    return_columns = ['Ref_Contig','Start_Ref','Ref_Pos',
            'Query_Contig','Start_Query','Query_Pos',
            'Ref_Loc','Query_Loc',
            'Ref_Start','Ref_End',
            'Query_Start','Query_End',
            'Ref_Base','Query_Base',
            'Dist_to_Ref_End','Dist_to_Query_End',
            'Ref_Aligned','Query_Aligned',
            'Perc_Iden','Cat']
    
    reverse_columns = ['Query_Contig','Start_Query','Query_Pos',
                    'Ref_Contig','Start_Ref','Ref_Pos',
                    'Query_Loc','Ref_Loc',
                    'Query_Start','Query_End',
                    'Ref_Start','Ref_End',
                    'Query_Base','Ref_Base',
                    'Dist_to_Query_End','Dist_to_Ref_End',
                    'Query_Aligned','Ref_Aligned',
                    'Query_Direction','Perc_Iden','Cat']
    
    reverse_complement = {'A':'T','T':'A','G':'C','C':'G',
                          'a':'T','t':'A','c':'G','g':'C'}
    
    # Columns to convert to integer
    int_columns = ['Start_Ref', 'Ref_Pos', 'Start_Query', 'Query_Pos',
                'Dist_to_Ref_End', 'Dist_to_Query_End', 'Ref_Aligned', 'Query_Aligned']

    # Columns to convert to float
    float_columns = ['Perc_Iden']
    
    if len(snp_rows) > 0:
        snp_df = pd.DataFrame(snp_rows, columns= snp_columns).copy()
        
        if snpdiffs_orientation == -1:
            snp_df = snp_df[reverse_columns].copy()
            snp_df.columns = snp_columns
            
            # Replace Query_Base and Reference_Base with reverse complement if Query_Direction is -1 and base is in ['A','T','G','C','a','c','t','g']
            snp_df.loc[snp_df['Query_Direction'] == '-1','Query_Base'] = snp_df.loc[snp_df['Query_Direction'] == '-1','Query_Base'].apply(lambda x: reverse_complement[x] if x in reverse_complement else x)
            snp_df.loc[snp_df['Query_Direction'] == '-1','Ref_Base'] = snp_df.loc[snp_df['Query_Direction'] == '-1','Ref_Base'].apply(lambda x: reverse_complement[x] if x in reverse_complement else x)
            
        for col in int_columns:
            snp_df.loc[:, col] = snp_df.loc[:, col].astype(int)
        for col in float_columns:
            snp_df.loc[:, col] = snp_df.loc[:, col].astype(float)
        
    else:
        snp_df = pd.DataFrame(columns = return_columns)

    return snp_df[return_columns]

def swapHeader(header_data):
            
    raw_header_cols = [x for x in header_data.columns]
    reverse_header_cols = [item.replace('Reference', 'temp').replace('Query', 'Reference').replace('temp', 'Query') for item in raw_header_cols]
    reversed_header_data = header_data[reverse_header_cols].copy()
    reversed_header_data.columns = raw_header_cols

    return reversed_header_data
    
def parseSNPDiffs(snpdiffs_file,snpdiffs_orientation):

    bed_rows = []
    snp_rows = []
    
    with open(snpdiffs_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line[0:2] == "#\t":
            pass
        elif line[0:3] == "##\t":
            bed_rows.append(line.strip().split("\t")[1:])
        else:
            snp_rows.append(line.strip().split("\t"))

    bed_df = processBED(bed_rows,snpdiffs_orientation)
    snp_df = processSNPs(snp_rows,snpdiffs_orientation)
    return (bed_df,snp_df)

def calculate_total_length(bedtool):
    return sum(len(interval) for interval in bedtool)

def filterSNPs(raw_snp_df,bed_df,log_file, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps):

    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
        
    total_snp_count = raw_snp_df.shape[0]
    orig_cols = raw_snp_df.columns
    
    # Create BED files for query and reference
    query_df = bed_df[['Query_Contig','Query_Start','Query_End','Ref_Aligned','Perc_Iden']].copy()
    ref_df = bed_df[['Ref_Contig','Ref_Start','Ref_End','Ref_Aligned','Perc_Iden']].copy()
    
    query_df = query_df[(query_df['Ref_Aligned'] >= min_len) & (query_df['Perc_Iden'] >= min_iden)].copy()
    ref_df = ref_df[(ref_df['Ref_Aligned'] >= min_len) & (ref_df['Perc_Iden'] >= min_iden)].copy()

    unique_ref_snps = raw_snp_df['Ref_Loc'].unique()
    unique_query_snps = raw_snp_df['Query_Loc'].unique()
    
    ref_snp_bed_df = pd.DataFrame([item.split('/') for item in unique_ref_snps], columns=['Ref_Contig','Ref_End'])
    ref_snp_bed_df['Ref_Start'] = ref_snp_bed_df['Ref_End'].astype(int) - 1
    ref_snp_bed = BedTool.from_dataframe(ref_snp_bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort()
    
    query_snp_bed_df = pd.DataFrame([item.split('/') for item in unique_query_snps], columns=['Query_Contig','Query_End'])
    query_snp_bed_df['Query_Start'] = query_snp_bed_df['Query_End'].astype(int) - 1
    query_snp_bed = BedTool.from_dataframe(query_snp_bed_df[['Query_Contig','Query_Start','Query_End']]).sort()
    
    # Get the coverage counts for each SNP
    intersected_query = query_snp_bed.intersect(BedTool.from_dataframe(query_df[['Query_Contig','Query_Start','Query_End']]), c=True)
    intersected_ref = ref_snp_bed.intersect(BedTool.from_dataframe(ref_df[['Ref_Contig','Ref_Start','Ref_End']]), c=True)
    
    query_snp_dict = {f"{interval.chrom}/{interval.end}": int(interval[-1]) for interval in intersected_query}
    ref_snp_dict = {f"{interval.chrom}/{interval.end}": int(interval[-1]) for interval in intersected_ref}
    
    with open(log_file,"a+") as log:
        log.write(f"\n\t- Raw SNP + Indel Count: {total_snp_count}\n")
        log.write("\n\t- Filtering SNPs and Indels that fall on short/poor alignments...\n")
    
    # Set all sites to SNP
    raw_snp_df['Filter_Cat'] = "SNP"
    
    # Filter out SNPs based on --min_len and --min_iden
    reject_length = raw_snp_df[(raw_snp_df['Ref_Aligned'] < min_len) & (raw_snp_df['Perc_Iden'] >= min_iden)].copy()
    reject_iden = raw_snp_df[(raw_snp_df['Ref_Aligned'] >= min_len) & (raw_snp_df['Perc_Iden'] < min_iden)].copy()
    reject_lenIden = raw_snp_df[(raw_snp_df['Ref_Aligned'] < min_len) & (raw_snp_df['Perc_Iden'] < min_iden)].copy()
    
    if reject_length.shape[0] > 0:
        reject_length['Filter_Cat'] = "Purged_Length"
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Alignment Length): {reject_length.shape[0]}\n")

            
    if reject_iden.shape[0] > 0:
        reject_iden['Filter_Cat'] = "Purged_Identity"
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Alignment Identity): {reject_iden.shape[0]}\n")


    if reject_lenIden.shape[0] > 0:
        reject_lenIden['Filter_Cat'] = "Purged_LengthIdentity"
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Alignment Length + Identity): {reject_lenIden.shape[0]}\n")

    pass_filter = raw_snp_df[(raw_snp_df['Ref_Aligned'] >= min_len) & (raw_snp_df['Perc_Iden'] >= min_iden)].copy().reset_index(drop=True)
    reject_filter = pd.concat([reject_length,reject_iden,reject_lenIden]).reset_index(drop=True)
    
    # Reduce any multibase indels to a single positions
    multisnp_reject_df = reject_filter.groupby('Ref_Loc').filter(lambda x: len(x['Ref_Base'].tolist()) > 1 & all(x['Ref_Base'] == '.')).drop_duplicates(subset=['Ref_Loc'], keep='first')
    reject_filter = pd.concat([reject_filter.loc[~reject_filter['Ref_Loc'].isin(multisnp_reject_df['Ref_Loc'].tolist())],multisnp_reject_df]).reset_index(drop=True)

    # Find rows in pass_filter where either Ref_Loc or Query_Loc occur more than once
    with open(log_file,"a+") as log:
        log.write("\n\t- Filtering SNPs and Indels to remove duplicate mappings or heterozygotes...\n")
    
    ref_loc_counts = pass_filter['Ref_Loc'].value_counts().reset_index().rename(columns={'Ref_Loc':'SNP_Count','index':'Ref_Loc',})
    query_loc_counts = pass_filter['Query_Loc'].value_counts().reset_index().rename(columns={'Query_Loc':'SNP_Count','index':'Query_Loc'})

    solo_ref_locs = ref_loc_counts[ref_loc_counts['SNP_Count'] == 1]['Ref_Loc']
    dup_ref_locs = ref_loc_counts[ref_loc_counts['SNP_Count'] > 1]['Ref_Loc']

    solo_query_locs = query_loc_counts[query_loc_counts['SNP_Count'] == 1]['Query_Loc']
    dup_query_locs = query_loc_counts[query_loc_counts['SNP_Count'] > 1]['Query_Loc']

    check_dup_df = pass_filter[(~pass_filter['Ref_Loc'].isin(solo_ref_locs)) | (~pass_filter['Query_Loc'].isin(solo_query_locs))].copy()
    pass_filter = pass_filter[(pass_filter['Ref_Loc'].isin(solo_ref_locs)) & (pass_filter['Query_Loc'].isin(solo_query_locs))].copy()
    
    # For any locs that have more than one SNP, keep the SNP with the highest Ref_Aligned and Perc_Iden
    if len(dup_ref_locs) > 0:
        for ref_loc in dup_ref_locs:
            ref_loc_overlaps = ref_snp_dict[ref_loc]
            ref_loc_df = check_dup_df[check_dup_df['Ref_Loc'] == ref_loc].copy()
            
            if (ref_loc_df['Ref_Base'] == ".").all():
                if ref_loc_overlaps > 1:
                    ref_loc_df['Filter_Cat'] = "Purged_Heterozygous_Ref"
                    reject_filter = pd.concat([reject_filter,ref_loc_df]).reset_index(drop=True)
                else:
                    ref_loc_df['Filter_Cat'] = "Purged_Multibase_Indel"
                    reject_filter = pd.concat([reject_filter,ref_loc_df]).reset_index(drop=True)
            elif ref_loc_df.shape[0] != ref_loc_overlaps:
                ref_loc_df['Filter_Cat'] = "Purged_Heterozygous_Ref"
                reject_filter = pd.concat([reject_filter,ref_loc_df]).reset_index(drop=True)
            elif ref_loc_df['Query_Base'].nunique() > 1:
                ref_loc_df['Filter_Cat'] = "Purged_Heterozygous"
                reject_filter = pd.concat([reject_filter,ref_loc_df]).reset_index(drop=True)
            else:
                best_snp = ref_loc_df.sort_values(by=['Ref_Aligned', 'Perc_Iden'], ascending=[False, False]).head(1)
                pass_filter = pd.concat([pass_filter,best_snp]).reset_index(drop=True)
                
                # Add remaining rows to reject_filter with the Filter_Cat of "Purged_Duplicate"
                dup_snps = ref_loc_df[~ref_loc_df.apply(lambda x: x in best_snp, axis=1)]
                dup_snps['Filter_Cat'] = "Purged_Duplicate"
                reject_filter = pd.concat([reject_filter,dup_snps]).reset_index(drop=True)
    
    if len(dup_query_locs) > 0:
        for query_loc in dup_query_locs:
            query_loc_overlaps = query_snp_dict[query_loc]
            query_loc_df = check_dup_df[check_dup_df['Query_Loc'] == query_loc].copy()
            
            # If query_loc_df['Ref_Base'] contains more than one value, or there are more overlaps than SNPs, all are heterozygotes
            if (query_loc_df['Query_Base'] == ".").all():
                if query_loc_overlaps > 1:
                    query_loc_df['Filter_Cat'] = "Purged_Heterozygous_Ref"
                    reject_filter = pd.concat([reject_filter,query_loc_df]).reset_index(drop=True)
                else:
                    query_loc_df['Filter_Cat'] = "Purged_Multibase_Indel"
                    reject_filter = pd.concat([reject_filter,query_loc_df]).reset_index(drop=True)
            elif query_loc_df.shape[0] != query_loc_overlaps:
                query_loc_df['Filter_Cat'] = "Purged_Heterozygous_Ref"
                reject_filter = pd.concat([reject_filter,query_loc_df])
            elif query_loc_df['Ref_Base'].nunique() > 1:
                query_loc_df['Filter_Cat'] = "Purged_Heterozygous"
                reject_filter = pd.concat([reject_filter,query_loc_df])
            else:                
                best_snp = query_loc_df.sort_values(by=['Ref_Aligned', 'Perc_Iden'], ascending=[False, False]).head(1)
                pass_filter = pd.concat([pass_filter,best_snp])
                
                # Add remaining rows to reject_filter with the Filter_Cat of "Purged_Duplicate"
                dup_snps = query_loc_df[~query_loc_df.apply(lambda x: x in best_snp, axis=1)]
                dup_snps['Filter_Cat'] = "Purged_Duplicate"
                reject_filter = pd.concat([reject_filter,dup_snps])

    purged_mindel_count = reject_filter[reject_filter['Filter_Cat'] == "Purged_Multibase_Indel"].shape[0]
    purged_dup_count = reject_filter[reject_filter['Filter_Cat'] == "Purged_Duplicate"].shape[0]
    purged_het_count = reject_filter[reject_filter['Filter_Cat'].isin(["Purged_Heterozygous","Purged_Heterozygous_Ref"])].shape[0]

    # Assert that Ref_Loc and Query_Loc are unique in pass_filter
    assert pass_filter['Ref_Loc'].nunique() == pass_filter.shape[0]
    assert pass_filter['Query_Loc'].nunique() == pass_filter.shape[0]
        
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Duplicate): {purged_dup_count}\n")
        log.write(f"\t\t- Purged (Heterozygote): {purged_het_count}\n")
        log.write(f"\t\t- Purged (Multibase Indel): {purged_mindel_count}\n")
    
    # Density filtering
    density_locs = []
    ref_locs = pass_filter['Ref_Loc'].tolist()
        
    if len(density_windows) == 0:
        with open(log_file,"a+") as log:
            log.write("\n\t- Density filtering disabled...\n")
    else:
        with open(log_file,"a+") as log:
            log.write("\n\t- Filtering SNPs based on reference genome density...\n")
        
        if len(ref_locs) > 0:
            density_df = pd.DataFrame([item.split('/') for item in ref_locs], columns=['Ref_Contig','Ref_End'])
            density_df['Ref_Start'] = density_df['Ref_End'].astype(int) - 1
            density_df['Ref_Loc'] = ref_locs
            density_bed = BedTool.from_dataframe(density_df[['Ref_Contig','Ref_Start','Ref_End']])

            for i in range(0,len(density_windows)):
                window_bed = density_bed.window(density_bed,c=True, w=density_windows[i])
                window_df = window_bed.to_dataframe()
                if window_df.shape[0] > 0:
                    window_df = window_df[window_df['name'] > max_snps[i]]
                    density_locs = density_locs + ["/".join([str(x[0]),str(x[1])]) for x in list(zip(window_df.chrom, window_df.end))]
                    density_bed = BedTool.from_dataframe(density_df[~density_df.Ref_Loc.isin(density_locs)][['Ref_Contig','Ref_Start','Ref_End']])

    reject_density = pass_filter[pass_filter['Ref_Loc'].isin(density_locs)].copy()
    pass_filter = pass_filter[~pass_filter['Ref_Loc'].isin(density_locs)].copy()
    
    if reject_density.shape[0] > 0:
        reject_density['Filter_Cat'] = "Purged_Density"
        reject_filter = pd.concat([reject_filter,reject_density]).reset_index(drop=True)
    
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Density): {reject_density.shape[0]}\n")
    
    with open(log_file,"a+") as log:
        log.write("\n\t- Filtering invalid sites...\n")
    
    reject_invalid = pass_filter[pass_filter['Cat'] == "Invalid"].copy()
    
    if reject_invalid.shape[0] > 0:
        reject_invalid['Filter_Cat'] = "Purged_Invalid"
        reject_filter = pd.concat([reject_filter,reject_invalid]).reset_index(drop=True)
    
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Invalid): {reject_invalid.shape[0]}\n")
    
    pass_filter = pass_filter[pass_filter['Cat'] != "Invalid"].copy()
    
    with open(log_file,"a+") as log:
        log.write("\n\t- Filtering indel sites...\n")
    
    reject_indel = pass_filter[(pass_filter['Cat'] == "Indel")].copy()
    
    if reject_indel.shape[0] > 0:
        reject_indel['Filter_Cat'] = "Purged_SNP_Indel"
        reject_filter = pd.concat([reject_filter,reject_indel]).reset_index(drop=True)
                                                              
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (SNP Indel): {reject_indel.shape[0]}\n")
    
    pass_filter = pass_filter[pass_filter['Cat'] != "Indel"].copy()
            
    with open(log_file,"a+") as log:
        log.write("\n\t- Filtering for SNPs near query and reference edges...\n")
    
    reject_query_edge = pass_filter[(pass_filter['Dist_to_Query_End'] < query_edge) & (pass_filter['Dist_to_Ref_End'] >= ref_edge)].copy()
    reject_ref_edge = pass_filter[(pass_filter['Dist_to_Ref_End'] < ref_edge) & (pass_filter['Dist_to_Query_End'] >= query_edge)].copy()
    reject_both_edge = pass_filter[(pass_filter['Dist_to_Query_End'] < query_edge) & (pass_filter['Dist_to_Ref_End'] < ref_edge)].copy()
    
    if reject_query_edge.shape[0] > 0:
        reject_query_edge['Filter_Cat'] = "Filtered_Query_Edge"
        reject_filter = pd.concat([reject_filter,reject_query_edge]).reset_index(drop=True)
    
    if reject_ref_edge.shape[0] > 0:
        reject_ref_edge['Filter_Cat'] = "Filtered_Ref_Edge"
        reject_filter = pd.concat([reject_filter,reject_ref_edge]).reset_index(drop=True)
        
    if reject_both_edge.shape[0] > 0:
        reject_both_edge['Filter_Cat'] = "Filtered_Both_Edge"
        reject_filter = pd.concat([reject_filter,reject_both_edge]).reset_index(drop=True)
        
    with open(log_file,"a+") as log:
        log.write(f"\t\t- Purged (Ref Edge): {reject_ref_edge.shape[0]}\n")
        log.write(f"\t\t- Purged (Query Edge): {reject_query_edge.shape[0]}\n")
        log.write(f"\t\t- Purged (Both Edge): {reject_both_edge.shape[0]}\n")
        
    pass_filter = pass_filter[(pass_filter['Dist_to_Query_End'] >= query_edge) & (pass_filter['Dist_to_Ref_End'] >= ref_edge)].copy()
    
    return_df = pd.concat([pass_filter,reject_filter]).reset_index(drop=True).sort_values(by=['Ref_Loc'])
    merged_df = return_df[orig_cols].merge(raw_snp_df[orig_cols], how='outer', indicator=True)
    assert merged_df[merged_df['_merge'] == 'both'].shape[0] == return_df.shape[0]
    helpers.cleanup(verbose=False,remove_all = False)
    
    return return_df.drop(columns=['Cat']).rename({'Filter_Cat':'Cat'}, axis=1)
    
def screenSNPDiffs(snpdiffs_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,reference_id,log_directory):
    
    screen_start_time = time.time()
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)

    # Set CSP2 variables to NA
    csp2_screen_snps = purged_length = purged_identity = purged_invalid = purged_indel = purged_lengthIdentity = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = "NA"
    filtered_snp_df = pd.DataFrame()
    good_reference_bed_df = pd.DataFrame()
    
    # Ensure snpdiffs file exists
    if not os.path.exists(snpdiffs_file) or not snpdiffs_file.endswith('.snpdiffs'):
        sys.exit(f"Invalid snpdiffs file provided: {snpdiffs_file}")
        
    # Ensure header can be read in
    try:
        header_data = fetchHeaders(snpdiffs_file)
        header_query = header_data['Query_ID'][0].replace(trim_name,'')
        header_ref = header_data['Reference_ID'][0].replace(trim_name,'')
    except:       
        sys.exit(f"Error reading headers from snpdiffs file: {snpdiffs_file}")
        
    # Check snpdiffs orientation
    if (header_ref == reference_id):
        snpdiffs_orientation = 1
        query_id = header_query
    elif (header_query == reference_id):
        snpdiffs_orientation = -1
        query_id = header_ref
        header_data = swapHeader(header_data)
    else:
        sys.exit(f"Error: Reference ID not found in header of {snpdiffs_file}...")      

    # Establish log file
    log_file = f"{log_directory}/{query_id}__vs__{reference_id}.log"
    with open(log_file,"w+") as log:
        log.write("Reference Screening for SNP Pipeline Analysis\n")
        log.write(f"Query Isolate: {query_id}\n")
        log.write(f"Reference Isolate: {reference_id}\n")
        log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
        log.write("-------------------------------------------------------\n\n")
        if snpdiffs_orientation == 1:
            log.write("\t- SNPDiffs file is in the forward orientation\n")
            log.write("-------------------------------------------------------\n\n")
        else:
            log.write("\t- SNPDiffs file is in the reverse orientation\n")
            log.write("-------------------------------------------------------\n\n")

    
    # Set variables from header data
    raw_snps = int(header_data['SNPs'][0])
    kmer_similarity = float(header_data['Kmer_Similarity'][0])
    shared_kmers = int(header_data['Shared_Kmers'][0])
    query_unique_kmers = int(header_data['Query_Unique_Kmers'][0])
    reference_unique_kmers = int(header_data['Reference_Unique_Kmers'][0])
    mummer_gsnps = int(header_data['gSNPs'][0])
    mummer_gindels = int(header_data['gIndels'][0])
    
    query_bases = int(header_data['Query_Assembly_Bases'][0])
    reference_bases = int(header_data['Reference_Assembly_Bases'][0])
    
    query_contigs = int(header_data['Query_Contig_Count'][0])
    reference_contigs = int(header_data['Reference_Contig_Count'][0])
    
    raw_query_percent_aligned = float(header_data['Query_Percent_Aligned'][0])
    raw_ref_percent_aligned = float(header_data['Reference_Percent_Aligned'][0])
        
    # If the reference is not covered by at least min_cov, STOP
    if raw_ref_percent_aligned < min_cov:
        query_percent_aligned = raw_query_percent_aligned
        reference_percent_aligned = raw_ref_percent_aligned
        screen_category = "Low_Coverage"
        with open(log_file,"a+") as log:
            log.write(f"\t- Reference genome coverage: {raw_ref_percent_aligned}% \n")
            log.write(f"\t- Query covers less than --min_cov ({min_cov}%)...Screen halted...\n")
            log.write("-------------------------------------------------------\n\n")

    elif raw_snps > 10000:
        query_percent_aligned = raw_query_percent_aligned
        reference_percent_aligned = raw_ref_percent_aligned
        screen_category = "SNP_Cutoff"
        with open(log_file,"a+") as log:
            log.write(f"\t- {raw_snps} detected...\n")
            log.write("\t- > 10,000 SNPs detected by MUMmer...Screen halted...\n")
            log.write("-------------------------------------------------------\n\n")
   
    else:
    
        ##### 02: Read in BED/SNP data #####
        with open(log_file,"a+") as log:
                log.write("Step 1: Reading in snpdiffs BED/SNP data...")
        try:
            bed_df,snp_df = parseSNPDiffs(snpdiffs_file,snpdiffs_orientation)
            
            with open(log_file,"a+") as log:
                log.write("Done!\n")
                log.write("-------------------------------------------------------\n\n")
                
        except:
            with open(log_file,"a+") as log:
                log.write(f"Error reading BED/SNP data from file: {snpdiffs_file}")
            sys.exit(f"Error reading BED/SNP data from file: {snpdiffs_file}")
            
        ##### 03: Filter genome overlaps #####
        with open(log_file,"a+") as log:
            log.write("Step 2: Filtering for short overlaps and low percent identity...")
        
        good_bed_df = bed_df[(bed_df['Ref_Aligned'] >= min_len) & (bed_df['Perc_Iden'] >= min_iden)].copy()
        
        if good_bed_df.shape[0] == 0:
            screen_category = "Low_Quality_Coverage"
            with open(log_file,"a+") as log:
                log.write(f"\n\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}) , no valid alignments remain...Screen halted...\n")
                log.write("-------------------------------------------------------\n\n")

                
        else:
            # Create a BED file for alignments that pass basic QC
            good_query_bed_df = good_bed_df[['Query_Contig','Query_Start','Query_End']].copy()
            good_reference_bed_df = good_bed_df[['Ref_Contig','Ref_Start','Ref_End']].copy()
            good_reference_bed_df.loc[:, 'Query_ID'] = query_id

            good_query_aligned = calculate_total_length(BedTool.from_dataframe(good_query_bed_df).sort().merge())
            good_reference_aligned = calculate_total_length(BedTool.from_dataframe(good_reference_bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort().merge())
            
            query_percent_aligned = (good_query_aligned / query_bases) * 100
            reference_percent_aligned = (good_reference_aligned / reference_bases) * 100
            
            if reference_percent_aligned < min_cov:
                screen_category = "Low_Quality_Coverage"
                with open(log_file,"a+") as log:
                    log.write(f"\n\t- Raw reference genome coverage was {raw_ref_percent_aligned}% \n")
                    log.write(f"\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}), reference genome coverage was {reference_percent_aligned:.2f}% \n")
                    log.write(f"\t- Query covers less than --min_cov ({min_cov}%) of reference after filtering...Screen halted...\n")
                    log.write("-------------------------------------------------------\n\n")

            else:
                screen_category = "Pass"
                with open(log_file,"a+") as log:
                    log.write("Done!\n")
                    log.write(f"\t- Raw reference genome coverage was {raw_ref_percent_aligned}% \n")
                    log.write(f"\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}), reference genome coverage was {reference_percent_aligned:.2f}% \n")
                    log.write("-------------------------------------------------------\n\n")

            
                # Filter SNPs
                with open(log_file,"a+") as log:
                    log.write("Step 3: Filtering SNPs to get final SNP distances...")
                
                if raw_snps == 0:
                    csp2_screen_snps = purged_length = purged_identity = purged_lengthIdentity = purged_indel = purged_invalid = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = 0
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write("\t- No SNPs detected in MUMmer output, no filtering required\n")
                        log.write("-------------------------------------------------------\n\n")

                else:
                    filtered_snp_df = filterSNPs(snp_df,bed_df,log_file, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps)
                    csp2_screen_snps = filtered_snp_df[filtered_snp_df.Cat == "SNP"].shape[0]
                    purged_length = filtered_snp_df[filtered_snp_df.Cat == "Purged_Length"].shape[0]
                    purged_identity = filtered_snp_df[filtered_snp_df.Cat == "Purged_Identity"].shape[0]
                    purged_lengthIdentity = filtered_snp_df[filtered_snp_df.Cat == "Purged_LengthIdentity"].shape[0]
                    
                    purged_duplicate = filtered_snp_df[filtered_snp_df.Cat == "Purged_Duplicate"].shape[0]
                    purged_het = filtered_snp_df[filtered_snp_df['Cat'].isin(["Purged_Heterozygous", "Purged_Heterozygous_Ref"])].shape[0]                    
                    purged_invalid = filtered_snp_df[filtered_snp_df.Cat == "Purged_Invalid"].shape[0]
                    purged_indel = filtered_snp_df[filtered_snp_df['Cat'].isin(["SNP_Indel", "Multibase_Indel"])].shape[0]                    

                    purged_density = filtered_snp_df[filtered_snp_df.Cat == "Purged_Density"].shape[0]
                    filtered_query_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Query_Edge"].shape[0]
                    filtered_ref_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Ref_Edge"].shape[0]
                    filtered_both_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Both_Edge"].shape[0]
                    filtered_edge = filtered_query_edge + filtered_ref_edge + filtered_both_edge
                    
                    # Write filtered SNP data to file
                    snp_file = log_file.replace(".log","_SNPs.tsv")
                    filtered_snp_df.to_csv(snp_file, sep="\t", index=False)
                    
                    filtered_snp_df.loc[:, 'Query_ID'] = query_id
                    
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"\t- {csp2_screen_snps} SNPs detected between {query_id} and {reference_id} after filtering\n")
                        log.write(f"\t- {filtered_edge} SNPs were ignored because of edge proximity (Query: {filtered_query_edge}, Ref: {filtered_ref_edge}, Both: {filtered_both_edge})\n")
                        log.write(f"\t- SNP data saved to {snp_file}\n")
                        log.write("-------------------------------------------------------\n\n")
    
    screen_end_time = time.time()
    with open(log_file,"a+") as log:
        log.write(f"Screening Time: {screen_end_time - screen_start_time:.2f} seconds\n")
    
    # Clean up pybedtools temp
    helpers.cleanup(verbose=False, remove_all=False)
    
    return ([str(item) for item in [query_id,reference_id,screen_category,csp2_screen_snps,
            f"{query_percent_aligned:.2f}",f"{reference_percent_aligned:.2f}",
            query_contigs,query_bases,reference_contigs,reference_bases,
            raw_snps,purged_length,purged_identity,purged_lengthIdentity,purged_invalid,purged_indel,purged_duplicate,purged_het,purged_density,
            filtered_query_edge,filtered_ref_edge,filtered_both_edge,
            kmer_similarity,shared_kmers,query_unique_kmers,reference_unique_kmers,
            mummer_gsnps,mummer_gindels]],good_reference_bed_df,filtered_snp_df)

def assessCoverage(query_id,site_list):
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
        
    if len(site_list) == 0:
        return pd.DataFrame(columns=['Ref_Loc','Query_ID','Cat'])
    else:
        coverage_df = pass_filter_coverage_df[pass_filter_coverage_df['Query_ID'] == query_id].copy()
        
        if coverage_df.shape[0] == 0:
            uncovered_loc_df = pd.DataFrame({
                'Ref_Loc': site_list,
                'Query_ID': [query_id] * len(site_list),
                'Cat': ["Uncovered"] * len(site_list)
            })
            return uncovered_loc_df
        else:
            coverage_bed = BedTool.from_dataframe(coverage_df[['Ref_Contig','Ref_Start','Ref_End']]).sort()
            snp_bed_df = pd.DataFrame([item.split('/') for item in site_list], columns=['Ref_Contig','Ref_End'])
            snp_bed_df['Ref_Start'] = snp_bed_df['Ref_End'].astype(int) - 1
            snp_bed_df['Ref_Loc'] = site_list
            snp_bed = BedTool.from_dataframe(snp_bed_df[['Ref_Contig','Ref_Start','Ref_End','Ref_Loc']]).sort()
            
            # Ref_Locs from snp_bed that intersect with coverage_bed go into covered_locs, the rest go into uncovered_locs
            covered_locs = snp_bed.intersect(coverage_bed, wa=True)
            uncovered_locs = snp_bed.intersect(coverage_bed, v=True, wa=True)

            covered_loc_df = pd.DataFrame({
                'Ref_Loc': [snp.fields[3] for snp in covered_locs],
                'Query_ID': [query_id] * covered_locs.count(),
                'Cat': ["Ref_Base"] * covered_locs.count()
            }) if covered_locs.count() > 0 else pd.DataFrame(columns=['Ref_Loc','Query_ID','Cat'])

            uncovered_loc_df = pd.DataFrame({
                'Ref_Loc': [snp.fields[3] for snp in uncovered_locs],
                'Query_ID': [query_id] * uncovered_locs.count(),
                'Cat': ["Uncovered"] * uncovered_locs.count()
            }) if uncovered_locs.count() > 0 else pd.DataFrame(columns=['Ref_Loc','Query_ID','Cat'])
            
            # Clean up pybedtools temp
            helpers.cleanup(verbose=False, remove_all=False)
            
            return pd.concat([covered_loc_df.drop_duplicates(['Ref_Loc']),uncovered_loc_df])

def getPairwise(pair, type = "Raw"):
    
    if type == "Preserved":
        seq1 = [record.seq for record in preserved_alignment if record.id == pair[0]][0]
        seq2 = [record.seq for record in preserved_alignment if record.id == pair[1]][0]
    else:
        seq1 = [record.seq for record in alignment if record.id == pair[0]][0]
        seq2 = [record.seq for record in alignment if record.id == pair[1]][0]
    
    snps_cocalled = sum((base1 in 'ACTGactg' and base2 in 'ACTGactg') for base1, base2 in zip(seq1, seq2))
    snp_distance = sum((base1 in 'ACTGactg' and base2 in 'ACTGactg' and base1 != base2) for base1, base2 in zip(seq1, seq2)) if snps_cocalled > 0 else np.nan
    return [pair[0],pair[1],snp_distance,snps_cocalled]

# Read in arguments
start_time = time.time()
reference_id = str(sys.argv[1])

output_directory = os.path.abspath(sys.argv[2])

log_directory = os.path.abspath(sys.argv[4])
log_file = f"{output_directory}/CSP2_SNP_Pipeline.log"

# Establish log file
with open(log_file,"w+") as log:
    log.write("CSP2 SNP Pipeline Analysis\n")
    log.write(f"Reference Isolate: {reference_id}\n")
    log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
    log.write("-------------------------------------------------------\n\n")
    log.write("Reading in SNPDiffs files...")

    
# Read in all lines and ensure each file exists
snpdiffs_list = [line.strip() for line in open(sys.argv[3], 'r')]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        sys.exit("Error: File does not exist: " + snpdiffs_file)

snpdiffs_list = list(set(snpdiffs_list))

if len(snpdiffs_list) == 0:
    sys.exit("No SNPdiffs files provided...")
    
with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Read in {len(snpdiffs_list)} SNPdiffs files\n")
    log.write("-------------------------------------------------------\n\n")

                
min_cov = float(sys.argv[5])
min_len = int(sys.argv[6])
min_iden = float(sys.argv[7])

ref_edge = int(sys.argv[8])
query_edge = int(sys.argv[9])

input_density = str(sys.argv[10])
input_maxsnps = str(sys.argv[11])

if input_density == "0":
    density_windows = []
    max_snps = []
else:
    density_windows = [int(x) for x in sys.argv[10].split(",")]
    max_snps = [int(x) for x in sys.argv[11].split(",")]
assert len(density_windows) == len(max_snps)

trim_name = sys.argv[12]
max_missing = float(sys.argv[13])
random_temp_id = str(uuid.uuid4())

global temp_dir
if sys.argv[14] != "":
    random_temp_id = str(uuid.uuid4())
    temp_dir = f"{os.path.normpath(os.path.abspath(sys.argv[14]))}/{random_temp_id}"
    try:
        os.mkdir(temp_dir)
        helpers.set_tempdir(temp_dir)
    except OSError as e:
        print(f"Error: Failed to create directory '{temp_dir}': {e}")
else:
    temp_dir = ""

try:
    # Establish output files
    reference_screening_file = f"{output_directory}/Reference_Screening.tsv"
    locus_category_file = f"{output_directory}/Locus_Categories.tsv"
    query_coverage_file = f"{output_directory}/Query_Coverage.tsv"
    raw_loclist = f"{output_directory}/snplist.txt"
    raw_alignment = f"{output_directory}/snpma.fasta"
    preserved_loclist = f"{output_directory}/snplist_preserved.txt"
    preserved_alignment_file = f"{output_directory}/snpma_preserved.fasta"
    raw_pairwise = f"{output_directory}/snp_distance_pairwise.tsv"
    raw_matrix = f"{output_directory}/snp_distance_matrix.tsv"
    preserved_pairwise = f"{output_directory}/snp_distance_pairwise_preserved.tsv"
    preserved_matrix = f"{output_directory}/snp_distance_matrix_preserved.tsv"

    with open(log_file,"a+") as log:
        log.write("Screening all queries against reference...")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(screenSNPDiffs,snp_diff_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,reference_id,log_directory) for snp_diff_file in snpdiffs_list]

    # Combine results into a dataframe
    output_columns = ['Query_ID','Reference_ID','Screen_Category','CSP2_Screen_SNPs',
                'Query_Percent_Aligned','Reference_Percent_Aligned',
                'Query_Contigs','Query_Bases','Reference_Contigs','Reference_Bases',
                'Raw_SNPs','Purged_Length','Purged_Identity','Purged_LengthIdentity','Purged_Invalid','Purged_Indel','Purged_Duplicate','Purged_Het','Purged_Density',
                'Filtered_Query_Edge','Filtered_Ref_Edge','Filtered_Both_Edge',
                'Kmer_Similarity','Shared_Kmers','Query_Unique_Kmers','Reference_Unique_Kmers',
                'MUMmer_gSNPs','MUMmer_gIndels']

    # Save reference screening
    results_df = pd.DataFrame([item.result()[0] for item in results], columns = output_columns)
    results_df.to_csv(reference_screening_file, sep="\t", index=False)

    # Get reference bed dfs
    covered_df = pd.concat([item.result()[1] for item in results])

    # Get snp dfs
    filtered_snp_df = pd.concat([item.result()[2] for item in results])

    # Separate isolates that pass QC
    pass_qc_isolates = list(set(results_df[results_df['Screen_Category'] == "Pass"]['Query_ID']))
    fail_qc_isolates = list(set(results_df[results_df['Screen_Category'] != "Pass"]['Query_ID']))

    if len(pass_qc_isolates) == 0:
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Reference screening data saved to {reference_screening_file}\n")
            log.write(f"\t- Of {len(snpdiffs_list)} comparisons, no isolates passed QC. Pipeline cannot continue.\n")
            log.write(f"\t- {len(fail_qc_isolates)} comparisons failed QC\n")
            for isolate in fail_qc_isolates:
                isolate_category = results_df[results_df['Query_ID'] == isolate]['Screen_Category'].values[0]
                log.write(f"\t\t- {isolate}: {isolate_category}\n")
            log.write("-------------------------------------------------------\n\n")
        sys.exit(0)
    else:    
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Reference screening data saved to {reference_screening_file}\n")
            log.write(f"\t- Of {len(snpdiffs_list)} comparisons, {len(pass_qc_isolates)} covered at least {min_cov}% of the reference genome after removing poor alignments\n")
            if len(fail_qc_isolates) > 0:
                log.write(f"\t- {len(fail_qc_isolates)} comparisons failed QC\n")
                for isolate in fail_qc_isolates:
                    isolate_category = results_df[results_df['Query_ID'] == isolate]['Screen_Category'].values[0]
                    log.write(f"\t\t- {isolate}: {isolate_category}\n")
            log.write("-------------------------------------------------------\n\n")

    with open(log_file,"a+") as log:
        log.write(f"Compiling SNPs across {len(pass_qc_isolates)} samples...\n")

    # Remove samples that failed QC
    pass_filter_snps = filtered_snp_df[filtered_snp_df['Query_ID'].isin(pass_qc_isolates)].copy()

    global pass_filter_coverage_df
    pass_filter_coverage_df = covered_df[covered_df['Query_ID'].isin(pass_qc_isolates)].copy()

    # Note SNPs lost irrevocably to reference edge trimming
    ref_edge_df = pass_filter_snps[pass_filter_snps['Cat'].isin(["Filtered_Ref_Edge",'Filtered_Both_Edge'])].copy()
    if ref_edge_df.shape[0] > 0:
        ref_edge_list = list(set(ref_edge_df['Ref_Loc']))
        with open(log_file,"a+") as log:
            log.write(f"\t- {len(ref_edge_list)} unique SNPs were within {ref_edge}bp of a reference contig end and were not considered in any query...\n")

    # Get SNP counts
    snp_df = pass_filter_snps[pass_filter_snps['Cat'] == "SNP"].copy()

    if snp_df.shape[0] == 0:
        snp_count = 0
        with open(log_file,"a+") as log:
            log.write("\t- No SNPs detected across all samples...Skipping to output...\n")    
            log.write("-------------------------------------------------------\n\n")
        pass

    else:
        snp_list = list(set(snp_df['Ref_Loc']))
        snp_count = len(snp_list)
        
        with open(log_file,"a+") as log:
            log.write(f"\t- {snp_count} unique SNPs detected across all samples...\n")        
        
        # Create Ref_Base df
        ref_base_df = snp_df[['Ref_Loc','Ref_Base']].copy().drop_duplicates().rename(columns = {'Ref_Base':'Query_Base'})
        ref_base_df['Query_ID'] = reference_id
        ref_base_df['Cat'] = "Reference_Isolate"
        ref_base_df = ref_base_df[['Ref_Loc','Query_ID','Query_Base','Cat']]
        
        # Rescue SNPs that are near the edge if they are valid SNPs in other samples
        rescued_edge_df = pass_filter_snps[(pass_filter_snps['Cat'] == "Filtered_Query_Edge") & (pass_filter_snps['Ref_Loc'].isin(snp_list))].copy()

        if rescued_edge_df.shape[0] > 0:
            rescued_counts = rescued_edge_df.groupby('Query_ID')['Ref_Loc'].count().reset_index().rename(columns={'Ref_Loc':'Rescued_Count'})
            with open(log_file,"a+") as log:
                log.write(f"\t- {rescued_edge_df.shape[0]} total SNPs were rescued from edge proximity filtering based on being more central in other samples...\n")
                for row in rescued_counts.itertuples():
                    log.write(f"\t\t- {row.Query_ID}: {row.Rescued_Count}\n")
            
            # Remove rescued sites from pass_filter_snps
            rescue_merge = pass_filter_snps.merge(rescued_edge_df, indicator=True, how='outer')
            pass_filter_snps = rescue_merge[rescue_merge['_merge'] == 'left_only'].drop(columns=['_merge']).copy()            
            
            # Add rescued SNPs to snp_df
            rescued_edge_df['Cat'] = "Rescued_SNP"    
            snp_df = pd.concat([snp_df,rescued_edge_df])
                    
        # Get purged SNPs 
        purged_df = pass_filter_snps[pass_filter_snps['Cat'] != "SNP"].copy()
        # Get purged SNPs where no query has a valid SNP
        non_snp_df = purged_df[~purged_df['Ref_Loc'].isin(snp_list)].copy()
        
        if non_snp_df.shape[0] > 0:
            non_snp_merge = purged_df.merge(non_snp_df, indicator=True, how='outer')
            purged_df = non_snp_merge[non_snp_merge['_merge'] == 'left_only'].drop(columns=['_merge']).copy()
            with open(log_file,"a+") as log:
                log.write(f"\t- {non_snp_df.shape[0]} total SNPs were purged in all queries they were found in, and were not considered in the final dataset...\n")

        purged_snp_df = pd.DataFrame(columns=['Ref_Loc','Query_ID','Query_Base','Cat'])
        if purged_df.shape[0] > 0:
            purged_counts = purged_df.groupby('Query_ID')['Ref_Loc'].count().reset_index().rename(columns={'Ref_Loc':'Purged_Count'})
            purged_df['Query_Base'] = "N"
            purged_snp_df = purged_df[['Ref_Loc','Query_ID','Query_Base','Cat']].copy()

            with open(log_file,"a+") as log:
                log.write(f"\t- {purged_df.shape[0]} unique SNPs were purged...\n")

        # Gather base data for all valid SNPs
        snp_base_df = snp_df[['Ref_Loc','Query_ID','Query_Base','Cat']].copy()
        
        # Genomic positions that do not occur- in the SNP data are either uncovered or match the reference base
        covered_snps = pd.concat([snp_base_df,purged_snp_df]).copy()
        ref_loc_sets = covered_snps.groupby('Query_ID')['Ref_Loc'].apply(set).to_dict()

        isolates_with_missing = [isolate for isolate in pass_qc_isolates if len(set(snp_list) - ref_loc_sets.get(isolate, set())) > 0]
        missing_df = pd.DataFrame(columns=['Ref_Loc','Query_ID','Query_Base','Cat'])
        
        if len(isolates_with_missing) > 0:
            
            isolate_data = [(isolate, list(set(snp_list) - ref_loc_sets.get(isolate, set()))) for isolate in isolates_with_missing]        
            
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = [executor.submit(assessCoverage, query, sites) for query, sites in isolate_data]

            coverage_df = pd.concat([item.result() for item in results])
        
            covered_df = coverage_df[coverage_df['Cat']=='Ref_Base'].copy()
            uncovered_df = coverage_df[coverage_df['Cat']=='Uncovered'].copy()

            if uncovered_df.shape[0] > 0:
                uncovered_df['Query_Base'] = "?"
                missing_df = pd.concat([missing_df,uncovered_df[['Ref_Loc','Query_ID','Query_Base','Cat']]])

                uncovered_counts = uncovered_df.groupby('Query_ID')['Ref_Loc'].count().reset_index().rename(columns={'Ref_Loc':'Uncovered_Count'})
                with open(log_file,"a+") as log:
                    log.write(f"\t- {uncovered_df.shape[0]} SNPs were not covered by the query isolate...\n")
            
            if covered_df.shape[0] > 0:
                ref_base_snp_df = covered_df.merge(ref_base_df[['Ref_Loc','Query_Base']], on='Ref_Loc', how='left')
                missing_df = pd.concat([missing_df,ref_base_snp_df[['Ref_Loc','Query_ID','Query_Base','Cat']]])
        
        final_snp_df = pd.concat([snp_base_df,purged_snp_df,missing_df]).sort_values(by=['Ref_Loc','Query_ID']).reset_index(drop=True)
        
        # Reduce Multibase indels to one site
        multibase_indel_df = final_snp_df[final_snp_df['Cat'] == 'Purged_Multibase_Indel'].drop_duplicates(subset=['Ref_Loc','Query_ID','Cat'], keep='first')
        final_snp_df = pd.concat([multibase_indel_df, final_snp_df[final_snp_df['Cat'] != 'Purged_Multibase_Indel']])    
        snp_counts = final_snp_df.groupby('Query_ID')['Ref_Loc'].count().reset_index().rename(columns={'Ref_Loc':'SNP_Count'})

        # Assert that all snp_counts == snp_count
        assert snp_counts['SNP_Count'].nunique() == 1
        assert snp_counts['SNP_Count'].values[0] == snp_count
        final_snp_df = pd.concat([final_snp_df,ref_base_df])

        # Get locus coverage stats        
        snp_coverage_df = final_snp_df[final_snp_df['Cat'].isin(['SNP','Rescued_SNP'])].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'SNP_Count'})
        ref_base_coverage_df = final_snp_df[final_snp_df['Cat'].isin(["Ref_Base","Reference_Isolate"])].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'Ref_Base_Count'}) 
        
        if uncovered_df.shape[0] > 0:
            uncovered_count_df = final_snp_df[final_snp_df['Cat'] == "Uncovered"].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'Uncovered_Count'}).copy()
        else:
            uncovered_count_df = pd.DataFrame(columns=['Ref_Loc','Uncovered_Count'])

        if purged_snp_df.shape[0] > 0:
            purged_count_df = final_snp_df[~final_snp_df['Cat'].isin(['SNP','Rescued_SNP','Reference_Isolate','Uncovered','Ref_Base'])].groupby('Ref_Loc')['Query_ID'].count().reset_index().rename(columns={'Query_ID':'Purged_Count'}).copy()
        else:
            purged_count_df = pd.DataFrame(columns=['Ref_Loc','Purged_Count'])

        locus_coverage_df = snp_coverage_df.merge(ref_base_coverage_df, how='outer', on='Ref_Loc').merge(uncovered_count_df, how='outer', on='Ref_Loc').merge(purged_count_df, how='outer', on='Ref_Loc').fillna(0)
        locus_coverage_df.loc[:, ['SNP_Count','Ref_Base_Count','Uncovered_Count','Purged_Count']] = locus_coverage_df.loc[:, ['SNP_Count','Ref_Base_Count','Uncovered_Count','Purged_Count']].astype(int)
        locus_coverage_df['Missing_Ratio'] = ((locus_coverage_df['Uncovered_Count'] + locus_coverage_df['Purged_Count']) / (1+len(pass_qc_isolates))) * 100
        locus_coverage_df.to_csv(locus_category_file, sep="\t", index=False)
        
        # Get isolate coverage stats
        min_isolate_cols = ['Query_ID','SNP','Ref_Base','Percent_Missing','Purged','Uncovered']
        isolate_coverage_df = final_snp_df.groupby('Query_ID')['Cat'].value_counts().unstack().fillna(0).astype(int).reset_index().drop(columns=['Reference_Isolate'])
        isolate_coverage_df.loc[isolate_coverage_df['Query_ID'] == reference_id, 'Ref_Base'] = snp_count
        
        if "Rescued_SNP" not in isolate_coverage_df.columns.tolist():
            pass
        else:
            isolate_coverage_df['SNP'] = isolate_coverage_df['SNP'] + isolate_coverage_df['Rescued_SNP']    
            isolate_coverage_df = isolate_coverage_df.drop(columns=['Rescued_SNP'])
        
        if "Uncovered" not in isolate_coverage_df.columns.tolist():
            isolate_coverage_df['Uncovered'] = 0
            
        purged_cols = [col for col in isolate_coverage_df.columns.tolist() if col not in min_isolate_cols]
        if len(purged_cols) > 0:
            isolate_coverage_df['Purged'] = isolate_coverage_df[purged_cols].sum(axis=1)
        else:
            isolate_coverage_df['Purged'] = 0
        
        isolate_coverage_df['Percent_Missing'] = (isolate_coverage_df['Uncovered'] + isolate_coverage_df['Purged'])/(isolate_coverage_df['Uncovered'] + isolate_coverage_df['Purged'] + isolate_coverage_df['Ref_Base'] + isolate_coverage_df['SNP']) * 100
        isolate_coverage_df = isolate_coverage_df[min_isolate_cols + purged_cols].sort_values(by = 'Percent_Missing',ascending = False).reset_index(drop=True)
        isolate_coverage_df.to_csv(query_coverage_file, sep="\t", index=False)
        
        with open(log_file,"a+") as log:
            log.write(f"\t- SNP coverage information: {locus_category_file}\n")
            log.write(f"\t- Query coverage information: {query_coverage_file}\n")
            log.write("-------------------------------------------------------\n\n")
        
        with open(log_file,"a+") as log:
            log.write("Processing alignment data...")

        alignment_df = final_snp_df[['Query_ID','Ref_Loc','Query_Base']].copy().rename(columns={'Query_Base':'Base'}).pivot(index='Query_ID', columns='Ref_Loc', values='Base')
        csp2_ordered = alignment_df.columns

        with open(raw_loclist,"w+") as loclist:
            loclist.write("\n".join(csp2_ordered)+"\n")
        
        seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in alignment_df.iterrows()]
        
        global alignment
        alignment = MultipleSeqAlignment(seq_records)
        AlignIO.write(alignment,raw_alignment,"fasta")

        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Saved alignment of {snp_count} SNPs to {raw_alignment}\n")
            log.write(f"\t- Saved ordered loc list to {raw_loclist}\n")
        
        global preserved_alignment
        if max_missing == float(100):
            locs_pass_missing = csp2_ordered
            preserved_alignment = alignment
            
            AlignIO.write(preserved_alignment,preserved_alignment_file,"fasta")
            with open(preserved_loclist,"w+") as loclist:
                loclist.write("\n".join(csp2_ordered)+"\n")
            
            with open(log_file,"a+") as log:
                log.write("Skipping SNP preservation step...\n")
                log.write(f"\t- Saved duplicate alignment to {preserved_alignment_file}\n")
                log.write(f"\t- Saved duplicate ordered loc list to {preserved_loclist}\n")
        else:
            with open(log_file,"a+") as log:
                log.write(f"Preserving SNPs with at most {max_missing}% missing data...\n")
            
            # Parse missing data
            locs_pass_missing = list(set(locus_coverage_df[locus_coverage_df['Missing_Ratio'] <= max_missing]['Ref_Loc']))
            
            if len(locs_pass_missing) == 0:
                with open(log_file,"a+") as log:
                    log.write(f"\t- Of {snp_count} SNPs, no SNPs pass the {max_missing}% missing data threshold...\n")
                    log.write("-------------------------------------------------------\n\n")
            else:
                preserved_alignment_df = alignment_df[locs_pass_missing].copy()
                
                preserved_ordered = preserved_alignment_df.columns
                with open(preserved_loclist,"w+") as loclist:
                    loclist.write("\n".join(preserved_ordered)+"\n")
                
                seq_records = [SeqRecord(Seq(''.join(row)), id=query,description='') for query,row in preserved_alignment_df.iterrows()]
                preserved_alignment = MultipleSeqAlignment(seq_records)
                AlignIO.write(preserved_alignment,preserved_alignment_file,"fasta")
                with open(log_file,"a+") as log:
                    log.write(f"\t- Of {snp_count} SNPs, {len(locs_pass_missing)} SNPs pass the {max_missing}% missing data threshold...\n")
                    log.write(f"\t- Saved preserved alignment to {preserved_alignment_file}\n")
                    log.write(f"\t- Saved preserved ordered loc list to {preserved_loclist}\n")
                    log.write("-------------------------------------------------------\n\n")

    with open(log_file,"a+") as log:
        log.write("Processing pairwise comparisons files...")

    # Get pairwise comparisons between all pass_qc_isolates and reference_id
    pairwise_combinations = [sorted(x) for x in list(combinations([reference_id] + pass_qc_isolates, 2))]
        
    if snp_count == 0:
        pairwise_df = pd.DataFrame([(pairwise[0], pairwise[1], 0,np.nan) for pairwise in pairwise_combinations],columns = ['Query_1','Query_2','SNP_Distance','SNPs_Cocalled'])
        preserved_pairwise_df = pairwise_df.copy()
        
        pairwise_df.to_csv(raw_pairwise, sep="\t", index=False)
        preserved_pairwise_df.to_csv(preserved_pairwise, sep="\t", index=False)
        
    else:
        raw_distance_results = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {executor.submit(getPairwise, pair) for pair in pairwise_combinations}
            for future in concurrent.futures.as_completed(futures):
                raw_distance_results.append(future.result())
        raw_pairwise_df = pd.DataFrame(raw_distance_results, columns=['Query_1', 'Query_2', 'SNP_Distance', 'SNPs_Cocalled'])
        raw_pairwise_df.to_csv(raw_pairwise, sep="\t", index=False)

        if len(locs_pass_missing) == snp_count:
            preserved_pairwise_df = raw_pairwise_df.copy()
            preserved_pairwise_df.to_csv(preserved_pairwise, sep="\t", index=False)
        elif len(locs_pass_missing) == 0:
            preserved_pairwise_df = pd.DataFrame([(pairwise[0], pairwise[1], 0,np.nan) for pairwise in pairwise_combinations],columns = ['Query_1','Query_2','SNP_Distance','SNPs_Cocalled'])
            preserved_pairwise_df.to_csv(preserved_pairwise, sep="\t", index=False)
        else:
            preserved_distance_results = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = {executor.submit(getPairwise, pair,"Preserved") for pair in pairwise_combinations}
                for future in concurrent.futures.as_completed(futures):
                    preserved_distance_results.append(future.result())
            preserved_pairwise_df = pd.DataFrame(preserved_distance_results, columns=['Query_1', 'Query_2', 'SNP_Distance', 'SNPs_Cocalled'])
            preserved_pairwise_df.to_csv(preserved_pairwise, sep="\t", index=False)
    
    # Create matrix
    idx = sorted(set(raw_pairwise_df['Query_1']).union(raw_pairwise_df['Query_2']))
    mirrored_distance_df = raw_pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
    mirrored_distance_df.index.name = ''
    mirrored_distance_df.to_csv(raw_matrix,sep="\t")
    
    idx = sorted(set(preserved_pairwise_df['Query_1']).union(preserved_pairwise_df['Query_2']))
    mirrored_distance_df = preserved_pairwise_df.pivot(index='Query_1', columns='Query_2', values='SNP_Distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).applymap(lambda x: format(x, '.0f'))
    mirrored_distance_df.index.name = ''
    mirrored_distance_df.to_csv(preserved_matrix,sep="\t")

    # Clean up pybedtools temp
    helpers.cleanup(verbose=False,remove_all = False)

    end_time = time.time()
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        if snp_count == 0:
            log.write(f"\t- No SNPs detected, zeroed pairwise distance files saved to {raw_pairwise}/{preserved_pairwise}/{raw_matrix}/{preserved_matrix}\n")
        else:
            log.write(f"\t- Saved raw pairwise distances to {raw_pairwise}\n")
            log.write(f"\t- Saved raw pairwise matrix to {raw_matrix}\n")

            if max_missing == float(100):
                log.write("Skipped SNP preservation step...\n")
                log.write(f"\t- Saved duplicated preserved pairwise distances to {preserved_pairwise}\n")
                log.write(f"\t- Saved duplicated preserved pairwise matrix to {preserved_matrix}\n")
            elif len(locs_pass_missing) == 0:
                log.write(f"\t- No SNPs passed the {max_missing}% missing data threshold, zeroed pairwise distance files saved to {preserved_pairwise}/{preserved_matrix}\n")
            else:
                log.write(f"\t- Saved preserved pairwise distances to {preserved_pairwise}\n")
                log.write(f"\t- Saved preserved pairwise matrix to {preserved_matrix}\n")
        log.write(f"Total Time: {end_time - start_time:.2f} seconds\n")
        log.write("-------------------------------------------------------\n\n")

except:
    print("Exception occurred:\n", traceback.format_exc())
    helpers.cleanup(verbose=False,remove_all = False)





#!/usr/bin/env python3

import sys
import os
import pandas as pd
import datetime
from pybedtools import BedTool,helpers
import concurrent.futures
import time

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
            
        # Remove any rows where Query_Contig or Ref_Contig == "." (Unaligned)
        covered_bed_df = bed_df[(bed_df['Ref_Start'] != ".") & (bed_df['Query_Start'] != ".")].copy()
        
        if covered_bed_df.shape[0] > 0:
            for col in int_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(int)
            for col in float_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float)
            return covered_bed_df        
        else:
            return pd.DataFrame(columns=bed_columns)           
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
        
    reverse_columns = ['Query_Contig','Start_Query','Query_Pos',
                    'Ref_Contig','Start_Ref','Ref_Pos',
                    'Query_Loc','Ref_Loc',
                    'Query_Start','Query_End',
                    'Ref_Start','Ref_End',
                    'Query_Base','Ref_Base',
                    'Dist_to_Query_End','Dist_to_Ref_End',
                    'Query_Aligned','Ref_Aligned',
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
    
    return return_df.drop(columns=['Cat']).rename({'Filter_Cat':'Cat'}, axis=1)
    
def screenSNPDiffs(snpdiffs_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,ref_ids):
    
    screen_start_time = time.time()

    # Set CSP2 variables to NA
    csp2_screen_snps = purged_length = purged_identity = purged_invalid = purged_indel = purged_lengthIdentity = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = "NA"

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
    if ref_ids == []:
        snpdiffs_orientation = 1
        query_id = header_query
        reference_id = header_ref
    elif (header_query not in ref_ids) and (header_ref in ref_ids):
        snpdiffs_orientation = 1
        query_id = header_query
        reference_id = header_ref
    elif (header_query in ref_ids) and (header_ref not in ref_ids):
        snpdiffs_orientation = -1
        query_id = header_ref
        reference_id = header_query
        header_data = swapHeader(header_data)
    else:
        snpdiffs_orientation = 2
        query_id = header_query
        reference_id = header_ref        

    # Establish log file
    log_file = f"{log_dir}/{query_id}__vs__{reference_id}.log"
    with open(log_file,"w+") as log:
        log.write("Screening Analysis\n")
        log.write(f"Query Isolate: {query_id}\n")
        log.write(f"Reference Isolate: {reference_id}\n")
        log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
        log.write("-------------------------------------------------------\n\n")
        if ref_ids == []:
            log.write("\t- No explicit references set, processing in forward orientation\n")
            log.write("-------------------------------------------------------\n\n")    
        elif snpdiffs_orientation == 1:
            log.write("\t- SNPDiffs file is in the forward orientation\n")
            log.write("-------------------------------------------------------\n\n")
        elif snpdiffs_orientation == -1:
            log.write("\t- SNPDiffs file is in the reverse orientation\n")
            log.write("-------------------------------------------------------\n\n")
        else:
            snpdiffs_orientation = 1
            log.write("\t- SNPDiffs file not contain a reference and non-reference sample, processing in forward orientation\n")
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
            
            good_query_aligned = calculate_total_length(BedTool.from_dataframe(good_query_bed_df).sort().merge())
            good_reference_aligned = calculate_total_length(BedTool.from_dataframe(good_reference_bed_df).sort().merge())
            
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
                    
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"\t- {csp2_screen_snps} SNPs detected between {query_id} and {reference_id} after filtering\n")
                        log.write(f"\t- {filtered_edge} SNPs were ignored because of edge proximity (Query: {filtered_query_edge}, Ref: {filtered_ref_edge}, Both: {filtered_both_edge})\n")
                        log.write(f"\t- SNP data saved to {snp_file}\n")
                        log.write("-------------------------------------------------------\n\n")
    
    screen_end_time = time.time()
    with open(log_file,"a+") as log:
        log.write(f"Screening Time: {screen_end_time - screen_start_time:.2f} seconds\n")
        
    return [str(item) for item in [query_id,reference_id,screen_category,csp2_screen_snps,
            f"{query_percent_aligned:.2f}",f"{reference_percent_aligned:.2f}",
            query_contigs,query_bases,reference_contigs,reference_bases,
            raw_snps,purged_length,purged_identity,purged_lengthIdentity,purged_invalid,purged_indel,purged_duplicate,purged_het,purged_density,
            filtered_query_edge,filtered_ref_edge,filtered_both_edge,
            kmer_similarity,shared_kmers,query_unique_kmers,reference_unique_kmers,
            mummer_gsnps,mummer_gindels]]

# Read in arguments
# Read in all lines and ensure each file exists
snpdiffs_list = [line.strip() for line in open(sys.argv[1], 'r')]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        sys.exit("Error: File does not exist: " + snpdiffs_file)        

snpdiffs_list = list(set(snpdiffs_list))

log_dir = os.path.normpath(os.path.abspath(sys.argv[2]))

min_cov = float(sys.argv[3])
min_len = int(sys.argv[4])
min_iden = float(sys.argv[5])

ref_edge = int(sys.argv[6])
query_edge = int(sys.argv[7])

input_density = str(sys.argv[8])
input_maxsnps = str(sys.argv[9])

if input_density == "0":
    density_windows = []
    max_snps = []
else:
    density_windows = [int(x) for x in sys.argv[8].split(",")]
    max_snps = [int(x) for x in sys.argv[9].split(",")]
assert len(density_windows) == len(max_snps)

trim_name = sys.argv[10]

output_file = os.path.abspath(sys.argv[11])

if os.stat(sys.argv[12]).st_size == 0:
    ref_ids = []
else:
    ref_ids = [line.strip() for line in open(sys.argv[12], 'r')]

if sys.argv[13] != "":
    helpers.set_tempdir(sys.argv[13])
    
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(screenSNPDiffs,snp_diff_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,ref_ids) for snp_diff_file in snpdiffs_list]

# Clean up pybedtools temp
helpers.cleanup(verbose=False, remove_all=False)

# Combine results into a dataframe
output_columns = ['Query_ID','Reference_ID','Screen_Category','CSP2_Screen_SNPs',
            'Query_Percent_Aligned','Reference_Percent_Aligned',
            'Query_Contigs','Query_Bases','Reference_Contigs','Reference_Bases',
            'Raw_SNPs','Purged_Length','Purged_Identity','Purged_LengthIdentity','Purged_Invalid','Purged_Indel','Purged_Duplicate','Purged_Het','Purged_Density',
            'Filtered_Query_Edge','Filtered_Ref_Edge','Filtered_Both_Edge',
            'Kmer_Similarity','Shared_Kmers','Query_Unique_Kmers','Reference_Unique_Kmers',
            'MUMmer_gSNPs','MUMmer_gIndels']

results_df = pd.DataFrame([item.result() for item in results], columns = output_columns)
results_df.to_csv(output_file, sep="\t", index=False)






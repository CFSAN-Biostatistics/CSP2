#!/usr/bin/env python3

import sys
import os
import pandas as pd
import datetime
from pybedtools import BedTool,helpers
import concurrent.futures
import time
import uuid
import traceback
import shutil
import argparse


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
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float).astype(int)
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
            snp_df.loc[:, col] = snp_df.loc[:, col].astype(float).astype(int)
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
    
    # Grab raw data    
    total_snp_count = raw_snp_df.shape[0]
    
    # Get unique SNPs relative to the reference genome
    unique_ref_snps = raw_snp_df['Ref_Loc'].unique()
    unique_snp_count = len(unique_ref_snps)
    
    snp_tally_df = pd.DataFrame()
    
    with open(log_file,"a+") as log:
        log.write(f"\n\t- Raw SNP + indel count: {total_snp_count}\n")
        log.write(f"\n\t- Unique SNP positions in reference genome: {unique_snp_count}\n")
    
    # Set all sites to SNP
    raw_snp_df['Filter_Cat'] = "SNP"
    
    # Filter out SNPs based on --min_len and --min_iden
    reject_length = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] < min_len) & (raw_snp_df['Perc_Iden'] >= min_iden)].copy()
    if reject_length.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Length): {reject_length.shape[0]}\n")
        reject_length['Filter_Cat'] = "Purged_Length"
        snp_tally_df = pd.concat([snp_tally_df,reject_length]).reset_index(drop=True)
        
    reject_iden = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] >= min_len) & (raw_snp_df['Perc_Iden'] < min_iden)].copy()
    if reject_iden.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Identity): {reject_iden.shape[0]}\n")
        reject_iden['Filter_Cat'] = "Purged_Identity"
        snp_tally_df = pd.concat([snp_tally_df,reject_iden]).reset_index(drop=True)

    reject_lenIden = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] < min_len) & (raw_snp_df['Perc_Iden'] < min_iden)].copy()
    if reject_lenIden.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Length + Identity): {reject_lenIden.shape[0]}\n")
        reject_lenIden['Filter_Cat'] = "Purged_LengthIdentity"
        snp_tally_df = pd.concat([snp_tally_df,reject_lenIden]).reset_index(drop=True)
    
    pass_filter = raw_snp_df.loc[(raw_snp_df['Ref_Aligned'] >= min_len) & (raw_snp_df['Perc_Iden'] >= min_iden)].copy().reset_index(drop=True)
    
    # Invalid processing
    reject_invalid = pass_filter[pass_filter['Cat'] == "Invalid"].copy()
    if reject_invalid.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Invalid Base): {reject_invalid.shape[0]}\n")
        reject_invalid['Filter_Cat'] = "Purged_Invalid"
        snp_tally_df = pd.concat([snp_tally_df,reject_invalid]).reset_index(drop=True)
    pass_filter = pass_filter.loc[pass_filter['Cat'] != "Invalid"].copy()
    
    # Indel processing
    reject_indel = pass_filter[pass_filter['Cat'] == "Indel"].copy()
    if reject_indel.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Indel): {reject_indel.shape[0]}\n")
        reject_indel['Filter_Cat'] = "Purged_Indel"
        snp_tally_df = pd.concat([snp_tally_df,reject_indel]).reset_index(drop=True)
    pass_filter = pass_filter.loc[pass_filter['Cat'] != "Indel"].copy()
    
    # Check for heterozygous SNPs
    check_heterozygous = pass_filter.groupby('Ref_Loc').filter(lambda x: x['Query_Base'].nunique() > 1)
    if check_heterozygous.shape[0] > 0:      
        reject_heterozygous = pass_filter.loc[pass_filter['Ref_Loc'].isin(check_heterozygous['Ref_Loc'])].copy()
        reject_heterozygous['Filter_Cat'] = "Purged_Heterozygous"
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Heterozygotes): {reject_heterozygous.shape[0]}\n")  
        snp_tally_df = pd.concat([snp_tally_df,reject_heterozygous]).reset_index(drop=True)
        pass_filter = pass_filter.loc[~pass_filter['Ref_Loc'].isin(check_heterozygous['Ref_Loc'])].copy()
   
    # Check for duplicate SNPs and take the longest, best hit
    check_duplicates = pass_filter.groupby('Ref_Loc').filter(lambda x: x.shape[0] > 1)
    if check_duplicates.shape[0] > 0:
        reject_duplicate = pass_filter.loc[pass_filter['Ref_Loc'].isin(check_duplicates['Ref_Loc'])].copy()
        pass_filter = pass_filter.loc[~pass_filter['Ref_Loc'].isin(check_duplicates['Ref_Loc'])].copy()
                
        best_snp = reject_duplicate.groupby('Ref_Loc').apply(lambda x: x.sort_values(by=['Ref_Aligned', 'Perc_Iden'], ascending=[False, False]).head(1))
        pass_filter = pd.concat([pass_filter,best_snp]).reset_index(drop=True)

        dup_snps = reject_duplicate[~reject_duplicate.apply(lambda x: x in best_snp, axis=1)]        
        dup_snps['Filter_Cat'] = "Purged_Duplicate"
        
        snp_tally_df = pd.concat([snp_tally_df,dup_snps]).reset_index(drop=True)
                        
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Duplicates): {dup_snps.shape[0]}\n")    
    
    # Assert that Ref_Loc and Query_Loc are unique in pass_filter
    helpers.cleanup(verbose=False,remove_all = False)
    assert pass_filter['Ref_Loc'].nunique() == pass_filter.shape[0]
    assert pass_filter['Query_Loc'].nunique() == pass_filter.shape[0]
    
    # Density filtering
    density_locs = []
    ref_locs = pass_filter['Ref_Loc'].tolist()
        
    if len(density_windows) == 0:
        with open(log_file,"a+") as log:
            log.write("\n\t- Density filtering disabled...\n")
    elif len(ref_locs) > 0:
        density_df = pd.DataFrame([item.split('/') for item in ref_locs], columns=['Ref_Contig','Ref_End'])
        density_df['Ref_Start'] = density_df['Ref_End'].astype(float).astype(int) - 1
        density_bed = BedTool.from_dataframe(density_df[['Ref_Contig','Ref_Start','Ref_End']])

        # For each density window, remove all SNPs that fall in a window with > max_snps
        for i in range(0,len(density_windows)):
            window_df = density_bed.window(density_bed,c=True, w=density_windows[i]).to_dataframe()
            problematic_windows = window_df[window_df['name'] > max_snps[i]].copy()
            if not problematic_windows.empty:
                temp_locs = []            
                for _, row in problematic_windows.iterrows():
                        purge_window_df = window_df[window_df['chrom'] == row['chrom']].copy()
                        purge_window_df['Dist'] = abs(purge_window_df['end'] - row['end'])
                        window_snps = purge_window_df.sort_values(by=['Dist'],ascending=True).head(row['name'])
                        temp_locs = temp_locs + ["/".join([str(x[0]),str(x[1])]) for x in list(zip(window_snps.chrom, window_snps.end))]
                density_locs.extend(list(set(temp_locs)))
    
    density_locs = list(set(density_locs))
    reject_density = pass_filter[pass_filter['Ref_Loc'].isin(density_locs)].copy() 
           
    if reject_density.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Density): {reject_density.shape[0]}\n")
        reject_density['Filter_Cat'] = "Purged_Density"
        snp_tally_df = pd.concat([snp_tally_df,reject_density]).reset_index(drop=True)
        pass_filter = pass_filter[~pass_filter['Ref_Loc'].isin(density_locs)].copy()

    reject_query_edge = pass_filter[(pass_filter['Dist_to_Query_End'] < query_edge) & (pass_filter['Dist_to_Ref_End'] >= ref_edge)].copy()
    reject_ref_edge = pass_filter[(pass_filter['Dist_to_Ref_End'] < ref_edge) & (pass_filter['Dist_to_Query_End'] >= query_edge)].copy()
    reject_both_edge = pass_filter[(pass_filter['Dist_to_Query_End'] < query_edge) & (pass_filter['Dist_to_Ref_End'] < ref_edge)].copy()
    
    if reject_query_edge.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Query Edge): {reject_query_edge.shape[0]}\n")
        reject_query_edge['Filter_Cat'] = "Filtered_Query_Edge"
        snp_tally_df = pd.concat([snp_tally_df,reject_query_edge]).reset_index(drop=True)
            
    if reject_ref_edge.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Ref Edge): {reject_ref_edge.shape[0]}\n")
        reject_ref_edge['Filter_Cat'] = "Filtered_Ref_Edge"
        snp_tally_df = pd.concat([snp_tally_df,reject_ref_edge]).reset_index(drop=True)
            
    if reject_both_edge.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Both Edge): {reject_both_edge.shape[0]}\n")
        reject_both_edge['Filter_Cat'] = "Filtered_Both_Edge"
        snp_tally_df = pd.concat([snp_tally_df,reject_both_edge]).reset_index(drop=True)
        
    pass_filter = pass_filter[(pass_filter['Dist_to_Query_End'] >= query_edge) & (pass_filter['Dist_to_Ref_End'] >= ref_edge)].copy()
    
    helpers.cleanup(verbose=False,remove_all = False)

    assert snp_tally_df.shape[0] + pass_filter.shape[0] == total_snp_count
    return_df = pd.concat([pass_filter,snp_tally_df]).reset_index(drop=True).sort_values(by=['Ref_Loc'])
        
    return return_df.drop(columns=['Cat']).rename({'Filter_Cat':'Cat'}, axis=1)
 

def screenSNPDiffs(snpdiffs_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,ref_ids):
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)

    screen_start_time = time.time()

    # Set CSP2 variables to NA
    csp2_screen_snps = purged_length = purged_identity = purged_invalid = purged_indel = purged_lengthIdentity = purged_duplicate = purged_het = purged_density = filtered_ref_edge = filtered_query_edge = filtered_both_edge = "NA"

    # Ensure snpdiffs file exists
    if not os.path.exists(snpdiffs_file) or not snpdiffs_file.endswith('.snpdiffs'):
        run_failed = True
        sys.exit(f"Invalid snpdiffs file provided: {snpdiffs_file}")
        
    # Ensure header can be read in
    try:
        header_data = fetchHeaders(snpdiffs_file)
        header_query = header_data['Query_ID'][0].replace(trim_name,'')
        header_ref = header_data['Reference_ID'][0].replace(trim_name,'')
    except:
        run_failed = True       
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
    raw_indels = int(header_data['Indels'][0])
    raw_invalid = int(header_data['Invalid'][0])
    
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

    elif raw_snps + raw_indels + raw_invalid > 10000:
        query_percent_aligned = raw_query_percent_aligned
        reference_percent_aligned = raw_ref_percent_aligned
        screen_category = "SNP_Cutoff"
        with open(log_file,"a+") as log:
            log.write(f"\t- {raw_snps} detected...\n")
            log.write("\t- > 10,000 SNPs, indels, or invalid sites detected by MUMmer...Screen halted...\n")
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
            run_failed = True
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
                    
                    # Write filtered SNP data to file
                    snp_file = log_file.replace(".log","_SNPs.tsv")
                    filtered_snp_df.to_csv(snp_file, sep="\t", index=False)
                    
                    csp2_screen_snps = filtered_snp_df[filtered_snp_df.Cat == "SNP"].shape[0]
                    
                    purged_length = filtered_snp_df[filtered_snp_df.Cat == "Purged_Length"].shape[0]
                    purged_identity = filtered_snp_df[filtered_snp_df.Cat == "Purged_Identity"].shape[0]
                    purged_lengthIdentity = filtered_snp_df[filtered_snp_df.Cat == "Purged_LengthIdentity"].shape[0]
                    purged_invalid = filtered_snp_df[filtered_snp_df.Cat == "Purged_Invalid"].shape[0]
                    purged_indel = filtered_snp_df[filtered_snp_df.Cat == "Purged_Indel"].shape[0]
                    purged_het = filtered_snp_df[filtered_snp_df.Cat == "Purged_Heterozygous"].shape[0]                    
                    purged_duplicate = filtered_snp_df[filtered_snp_df.Cat == "Purged_Duplicate"].shape[0]
                    purged_density = filtered_snp_df[filtered_snp_df.Cat == "Purged_Density"].shape[0]
                    filtered_query_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Query_Edge"].shape[0]
                    filtered_ref_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Ref_Edge"].shape[0]
                    filtered_both_edge = filtered_snp_df[filtered_snp_df.Cat == "Filtered_Both_Edge"].shape[0]
                   
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"\t- {csp2_screen_snps} SNPs detected between {query_id} and {reference_id} after filtering\n")
                        log.write(f"\t- SNP data saved to {snp_file}\n")
                        log.write("-------------------------------------------------------\n\n")
    
    screen_end_time = time.time()
    helpers.cleanup(verbose=False, remove_all=False)

    with open(log_file,"a+") as log:
        log.write(f"Screening Time: {screen_end_time - screen_start_time:.2f} seconds\n")
    
    # Clean up pybedtools temp
    helpers.cleanup(verbose=False, remove_all=False)
    
    return [str(item) for item in [query_id,reference_id,screen_category,csp2_screen_snps,
            f"{query_percent_aligned:.2f}",f"{reference_percent_aligned:.2f}",
            query_contigs,query_bases,reference_contigs,reference_bases,
            raw_snps,purged_length,purged_identity,purged_lengthIdentity,purged_invalid,purged_indel,purged_duplicate,purged_het,purged_density,
            filtered_query_edge,filtered_ref_edge,filtered_both_edge,
            kmer_similarity,shared_kmers,query_unique_kmers,reference_unique_kmers,
            mummer_gsnps,mummer_gindels]]

# Read in arguments
global run_failed
run_failed = False

parser = argparse.ArgumentParser()
parser.add_argument("--snpdiffs_file", help="Path to the file containing SNP diffs")
parser.add_argument("--log_dir", help="Path to the log directory")
parser.add_argument("--min_cov", default=85, type=float, help="Minimum coverage")
parser.add_argument("--min_len", default=500,type=int, help="Minimum length")
parser.add_argument("--min_iden", default=99,type=float, help="Minimum identity")
parser.add_argument("--ref_edge", default=150,type=int, help="Reference edge")
parser.add_argument("--query_edge", default=150,type=int, help="Query edge")
parser.add_argument("--density_windows",default="1000,125,15", help="Density windows (comma-separated)")
parser.add_argument("--max_snps", default="3,2,1",help="Maximum SNPs (comma-separated)")
parser.add_argument('--trim_name', type=str, default="", help='trim name')
parser.add_argument("--output_file", help="Output file")
parser.add_argument("--ref_id", help="Reference IDs file")
parser.add_argument("--tmp_dir",default="", help="TMP dir")

args = parser.parse_args()

snpdiffs_list = [line.strip() for line in open(args.snpdiffs_file, 'r')]
snpdiffs_list = [line for line in snpdiffs_list if line]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        run_failed = True
        sys.exit("Error: File does not exist: " + snpdiffs_file)

snpdiffs_list = list(set(snpdiffs_list))

log_dir = os.path.normpath(os.path.abspath(args.log_dir))

min_cov = args.min_cov
min_len = args.min_len
min_iden = args.min_iden

ref_edge = args.ref_edge
query_edge = args.query_edge

input_density = args.density_windows
input_maxsnps = args.max_snps

if input_density == "0":
    density_windows = []
    max_snps = []
else:
    density_windows = [int(x) for x in args.density_windows.split(",")]
    max_snps = [int(x) for x in args.max_snps.split(",")]
assert len(density_windows) == len(max_snps)

trim_name = args.trim_name

output_file = os.path.abspath(args.output_file)

if os.stat(args.ref_id).st_size == 0:
    ref_ids = []
else:
    ref_ids = [line.strip() for line in open(args.ref_id, 'r')]

global temp_dir
if args.tmp_dir != "":
    random_temp_id = str(uuid.uuid4())
    temp_dir = f"{os.path.normpath(os.path.abspath(args.tmp_dir))}/{random_temp_id}"
    try:
        os.mkdir(temp_dir)
        helpers.set_tempdir(temp_dir)
    except OSError as e:
        run_failed = True
        print(f"Error: Failed to create directory '{temp_dir}': {e}")
else:
    temp_dir = ""

try:    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(screenSNPDiffs,snp_diff_file,trim_name, min_cov, min_len, min_iden, ref_edge, query_edge, density_windows, max_snps,ref_ids) for snp_diff_file in snpdiffs_list]

    # Clean up pybedtools temp
    helpers.cleanup(verbose=False,remove_all = False)

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
except:
    run_failed = True
    print("Exception occurred:\n", traceback.format_exc())
finally:
    helpers.cleanup(verbose=False, remove_all=False)
    if temp_dir != "":
        shutil.rmtree(temp_dir)
    if run_failed:
        sys.exit(1)






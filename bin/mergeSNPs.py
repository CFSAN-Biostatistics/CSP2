#!/usr/bin/python3

import os
import sys
import pandas as pd
import numpy as np
from glob import glob
from collections import defaultdict
from itertools import combinations
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import concurrent.futures
import time
from collections import Counter
import datetime

def getPairwise(compare,alignment):

    # Define valid bases
    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']

    # Populate temp alignment with two taxa from compare
    temp_alignment = alignment[0:0]
    for i in compare:
        temp_alignment.append(SeqRecord(Seq(str(alignment[snp_isolates.index(i)].seq)), id=str(alignment[snp_isolates.index(i)].id)))
    
    # If resulting alignment is empty, raise exception
    if int(temp_alignment.get_alignment_length()) == 0:
        sys.exit("ERROR: Alignment processed, but appears to have no bases...")
    else:
        align_length = temp_alignment.get_alignment_length()
    
    # Get columns where both species have bases, and if so, find identical sites
    cols_with_base = []
    identical_sites = []
    
    for i in range(0,align_length):
        seq_string = temp_alignment[:, i]
            
        # Get alleles
        allele_list = list(Counter(seq_string))
        allele_count = len(allele_list)
        
        if all(b in bases for b in seq_string):
            cols_with_base.append(i)
            
            if allele_count == 1:
                identical_sites.append(i)

    compare_id = ";".join(sorted([compare[0],compare[1]]))
    
    if len(cols_with_base) == 0:
        if not len(identical_sites) == 0:
            sys.exit("ERROR: Identical sites found but no co-called?")
        else:
            return(pd.DataFrame([[compare[0],compare[1],compare_id,0,0,0]],columns=['Sample_A','Sample_B','Comparison','Cocalled_Sites','Identical_Sites','Prop_Identical']))
    else:
        return(pd.DataFrame([[compare[0],compare[1],compare_id,len(cols_with_base),len(identical_sites),float(len(identical_sites))/float(len(cols_with_base))]],columns=['Sample_A','Sample_B','Comparison','Cocalled_Sites','Identical_Sites','Prop_Identical']))

def process_yenta_tuple(yenta_tuple,isolate_dict):
    isolate, loc = yenta_tuple
    to_add = [isolate + "_yenta_" + loc] + [key[2] + "_yenta_" + value for key, value in isolate_dict.items() if (key[0], key[1]) == (isolate, loc)]
    return list(set(sorted(to_add, key=lambda x: x[0])))

def dfs(node, current_list):
    visited.add(node)
    current_list.append(node)
    for neighbor in graph[node]:
        if neighbor not in visited:
            dfs(neighbor, current_list)

# Function to return alignments and coordinates for Yenta SNPs
def getIsolateData(tuple_list):
    
    start_time = time.time()
    snp_alignment = MultipleSeqAlignment([])

    # Get loc data
    loc_list = []
    for isolate in snp_isolates:
        if isolate in [item[0] for item in tuple_list]:
            loc = [item[1] for item in tuple_list if item[0] == isolate][0]
            loc_list.append((isolate,loc))
        else:
            loc_list.append((isolate,"NA"))
    
    loc_df = pd.DataFrame(columns=snp_isolates)
    loc_df.loc[len(loc_df)] = [loc[1] for loc in loc_list]
    loc_time = time.time()

    # Get isolates with data
    site_isolates = [loc[0] for loc in loc_list if not (loc[1] == "NA")]
    reference_isolate = site_isolates[0]
    missing_isolates = [loc[0] for loc in loc_list if loc[1] == "NA"]

    # Create site dict
    site_dict = {key: val for key, val in snp_dict.items() if (key[2], val[0]) in loc_list or (key[0], key[1]) in loc_list}    
    dict_list = [{'Focal': key[0], 'Focal_Loc': key[1], 'Compare': key[2], 'Compare_Loc': value[0],'Focal_Base':value[1],'Focal_Direction':value[2],'Compare_Base':value[3],'Compare_Direction':value[4]} for key, value in site_dict.items()]
    site_df = pd.DataFrame(dict_list)
    site_df['Sort'] = site_df.apply(lambda row: tuple(sorted([row['Focal'], row['Compare']])), axis=1)
    site_df = site_df[~site_df.duplicated(subset='Sort', keep='first')]
    site_df.drop('Sort', axis=1, inplace=True)
    dict_time = time.time()

    # Get base data
    base_list = []
    if len(site_isolates) < 2:
        exit("Fewer than two isolates can't have a SNP?")
    elif len(site_isolates) == 2:
        for isolate in snp_isolates:
            if isolate in site_isolates:
                if site_df['Focal'][0] == isolate:
                    base_list.append((isolate,site_df['Focal_Base'][0]))
                else:
                    base_list.append((isolate,site_df['Compare_Base'][0]))
            else:
                base_list.append((isolate,"N"))
        set_bases = base_list
    else:
        if site_df['Focal_Direction'].nunique() == 1 and site_df['Compare_Direction'].nunique() == 1:
            for isolate in snp_isolates:
                if isolate in site_isolates:
                    isolate_row = site_df.loc[(site_df['Focal'] == isolate) | (site_df['Compare'] == isolate)].iloc[0]
                    if isolate_row['Focal'] == isolate:
                        base_list.append((isolate,isolate_row['Focal_Base']))
                    else:
                        base_list.append((isolate,isolate_row['Compare_Base']))
                else:
                    base_list.append((isolate,"N"))
            set_bases = base_list
        else:
            reverse_complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G','a':'t','t':'a','c':'g','g':'c'}
            
            # Process reference sample
            reference_comparisons = site_df.loc[(site_df['Focal'] == reference_isolate) | (site_df['Compare'] == reference_isolate)]
            
            set_bases = []

            focal_isolate = reference_comparisons['Focal'][0]
            focal_direction = reference_comparisons['Focal_Direction'][0]
            compare_direction = reference_comparisons['Compare_Direction'][0]
            focal_base = reference_comparisons['Focal_Base'][0]
            compare_base = reference_comparisons['Compare_Base'][0]

            if focal_isolate == reference_isolate:
                if focal_direction == 1:
                    ref_base = (reference_isolate,focal_base)
                else:
                    ref_base = (reference_isolate,reverse_complements[focal_base])
            else:
                if compare_direction == 1:
                    ref_base = (reference_isolate,compare_base)
                else:
                    ref_base = (reference_isolate,reverse_complements[compare_base])

            set_bases.append(ref_base)

            # Set bases based on ref base
            for index, row in reference_comparisons.iterrows():
                if row['Focal'] == reference_isolate:
                    if row['Focal_Base'] == ref_base[1]:
                        set_bases.append((row['Compare'], row['Compare_Base']))
                    else:
                        set_bases.append((row['Compare'], reverse_complements[row['Compare_Base']]))
                else:
                    if row['Compare_Base'] == ref_base:
                        set_bases.append((row['Focal'], row['Focal_Base']))
                    else:
                        set_bases.append((row['Focal'], reverse_complements[row['Focal_Base']]))

            # Any isolates with data not in the reference comparisons matches the reference
            set_isolates =  [isolate[0] for isolate in set_bases]
            for isolate in [isolate for isolate in site_isolates if isolate not in set_isolates]:
                set_bases.append((isolate,ref_base[1]))

            if len(missing_isolates) > 0: 
                for isolate in missing_isolates:
                    set_bases.append((isolate,"N"))
    
    base_time = time.time()
    if len(set_bases) == len(snp_isolates):
        # Add base data to alignment
        for isolate in snp_isolates:
            iso_base = [item[1] for item in set_bases if item[0] == isolate][0]
            snp_alignment.append(SeqRecord(Seq(iso_base), id=str(isolate),description=""))
    else:
        print(site_isolates)
        print(set_bases)
        exit("Alignment error")
    

    for_loc = loc_time - start_time
    for_dict = dict_time - loc_time
    for_base = base_time - dict_time

    with open(log_file,"a+") as log:
        log.write(f"\t- Locs: {for_loc:.2f}s; Dict: {for_dict:.2f}; Base: {for_base:.2f}\n")
    return [snp_alignment,loc_df]

# Read in arguments
output_dir = os.path.abspath(sys.argv[1])
mummer_dir = os.path.abspath(str(sys.argv[2]))
snp_dir = os.path.abspath(str(sys.argv[3]))
align_cov = float(sys.argv[4])
max_perc_n = float(sys.argv[5])


# Create logfiles
global log_file
log_file = snp_dir+"/SNP_Analysis.log"
with open(log_file,"w+") as log:
    log.write("-------------------------------------------------------\n")
    log.write("Yenta SNP Analysis\n")
    log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
    log.write("-------------------------------------------------------\n\n")

# Read in SNP files

snp_files = sorted(glob(mummer_dir+"/*MUmmer_SNPs.tsv"))
if len(snp_files) <= 0:
    with open(log_file,"a+") as log:
        log.write("ERROR: No files ending in _MUmmer_SNPs.tsv in "+mummer_dir+"\n")
    exit("No files ending in _MUmmer_SNPs.tsv in "+mummer_dir)

else:
    # Read in pairwise file
    with open(log_file,"a+") as log:
        log.write("Step 1: Reading in pairwise distance file...")
    try:
        pairwise_df = pd.read_table(output_dir+"/Raw_Pairwise_Distances.tsv",sep="\t")
        full_isolate_list = sorted(list(set(pairwise_df.Query_ID) | set(pairwise_df.Reference_ID)))
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write("\t- Found pairwise file ("+output_dir+"/Raw_Pairwise_Distances.tsv), which contains data on " + str(len(full_isolate_list)) + " isolates...\n")
            log.write("\n-------------------------------------------------------\n\n")
    except:
        with open(log_file,"a+") as log:
            log.write("ERROR: Cannot find pairwise file at "+output_dir+"/Raw_Pairwise_Distances.tsv\n")
        exit("Cannot find pairwise file at "+output_dir+"/Raw_Pairwise_Distances.tsv")


    
    # Get median coverage stats
    with open(log_file,"a+") as log:
        log.write("Step 2: Removing isolates with a mean genome coverage < "+str(align_cov) + "% ...")
    median_cov_df = pd.DataFrame({'ID':  pd.concat([pairwise_df['Query_ID'], pairwise_df['Reference_ID']]), 'Coverage': pd.concat([pairwise_df['Percent_Reference_Covered'], pairwise_df['Percent_Query_Covered']])}).groupby('ID')['Coverage'].median().reset_index()

    # Extract isolates with a median percent coverage < min_median_coverage
    low_cov_ids = median_cov_df[median_cov_df['Coverage'] < align_cov]['ID'].tolist()
    global snp_isolates
    snp_isolates = median_cov_df[median_cov_df['Coverage'] >= align_cov]['ID'].tolist()
    filtered_snp_files = [filename for filename in snp_files if not any(("/" + prefix + "_") in filename or ("_" + prefix + "_") in filename for prefix in low_cov_ids)]
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        
    if len(snp_isolates) < 2:
        with open(log_file,"a+") as log:
            log.write("\t- Fewer than two isolates remain after filtering for low coverage.\n")
        sys.exit()
    elif len(snp_isolates) == 2:
        with open(log_file,"a+") as log:
            log.write("\t- Two isolates remain after filtering for low coverage, results from pairwise will not change.\n")
            sys.exit()
    else:
        if len(low_cov_ids) > 0:
            with open(log_file,"a+") as log:
                log.write("\t- The median genome coverage for the following isolate(s) are below the level set via --align_cov ("+str(align_cov) + "):\n")
                for low_cov in low_cov_ids:
                    log.write("\t- "+low_cov+"\n")
                    log.write("\n-------------------------------------------------------\n\n")
        else:
            with open(log_file,"a+") as log:
                log.write("\t- No isolates were purged due to genomic coverage.\n")
                log.write("\n-------------------------------------------------------\n\n")

    # Read in filtered snp files
    try:
        start_time = time.time()
        with open(log_file,"a+") as log:
            log.write("Step 3: Reading in SNP files...")
        raw_snp_df = pd.concat((pd.read_table(filename,sep="\t") for filename in filtered_snp_files))[['Ref','Ref_Loc','Ref_Base','Ref_Direction','Query','Query_Loc','Query_Base','Query_Direction','Cat']]
        end_time = time.time()
        snp_time = end_time - start_time
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Read in "+str(len(filtered_snp_files)) + f" SNP files in {snp_time:.2f}s\n")
            log.write("\n-------------------------------------------------------\n\n")

    except:
        with open(log_file,"a+") as log:
            log.write("\n\t- ERROR: Cannot load SNP files in "+mummer_dir+"\n")
        exit("Cannot load SNP files in "+mummer_dir)     
    
    if raw_snp_df.shape[0] == 0:
        with open(log_file,"a+") as log:
            log.write("\n\t- SNP files located, but no SNPs present. Cannot proceed...\n")
        sys.exit()
    
    # Separate purged locs
    with open(log_file,"a+") as log:
        log.write("Step 4: Removing purged sites...")
    purged_df = raw_snp_df[raw_snp_df.Cat.str.startswith("Purged_")]
    purged_tuples = set(zip(purged_df['Ref'], purged_df['Ref_Loc'])) | set(zip(purged_df['Query'], purged_df['Query_Loc']))
    if purged_df.shape[0] > 0:
        purged_counts = pd.concat([purged_df['Ref'], purged_df['Query']]).value_counts()
        pd.DataFrame({'ID': purged_counts.index, 'Count': purged_counts.values}).sort_values(by='Count', ascending=False).to_csv(snp_dir+"/Purged_SNPs_by_Isolate.tsv",index=False,sep="\t")
    
    # Process non-purged locs
    nonpurged_df = raw_snp_df[~raw_snp_df.Cat.str.startswith("Purged_")] 
    
    with open(log_file,"a+") as log:
        log.write("Done!\n")
    
    if nonpurged_df.shape[0] <= 0:
        with open(log_file,"a+") as log:
            log.write("\t- All SNPs purged...\n")

    # Get Yenta sites
    yenta_df = raw_snp_df[raw_snp_df.Cat == "Yenta_SNP"]
    yenta_tuples = set(zip(yenta_df['Ref'], yenta_df['Ref_Loc'])) | set(zip(yenta_df['Query'], yenta_df['Query_Loc']))
    yenta_tuples = [tuple_pair for tuple_pair in yenta_tuples if ~(tuple_pair in purged_tuples)]
    if yenta_df.shape[0] == 0:
        with open(log_file,"a+") as log:
            log.write("\t- SNP files located, but no Yenta SNPs present. Cannot proceed...\n")
        sys.exit()
    else:
        with open(log_file,"a+") as log:
            log.write("\n-------------------------------------------------------\n\n")

    # Fetch loc data
    with open(log_file,"a+") as log:
        log.write("Step 5: Create dictionary of locs...")
    
    start_time = time.time()
    global snp_dict
    snp_dict = defaultdict(lambda:(np.na))
    
    for isolate in snp_isolates:
                
        # Create df based on focal isolate
        isolate_df = nonpurged_df[nonpurged_df.Ref == isolate].rename(columns = {"Ref": "Focal_Isolate", "Ref_Loc": "Focal_Loc", "Ref_Base":"Focal_Base","Ref_Direction":"Focal_Direction","Query":"Compare_Isolate","Query_Loc":"Compare_Loc","Query_Base":"Compare_Base","Query_Direction":"Compare_Direction"})[["Focal_Isolate","Focal_Loc","Focal_Base","Focal_Direction","Compare_Isolate","Compare_Loc","Compare_Base","Compare_Direction"]].append(nonpurged_df[nonpurged_df.Query == isolate].rename(columns = {"Query": "Focal_Isolate", "Query_Loc": "Focal_Loc", "Query_Base":"Focal_Base","Query_Direction":"Focal_Direction","Ref":"Compare_Isolate","Ref_Loc":"Compare_Loc","Ref_Base":"Compare_Base","Ref_Direction":"Compare_Direction"})[["Focal_Isolate","Focal_Loc","Focal_Base","Focal_Direction","Compare_Isolate","Compare_Loc","Compare_Base","Compare_Direction"]])

        # Add to SNP dict unless counterpart is already there
        for index, row in isolate_df.iterrows():
            if (row['Compare_Isolate'],row['Compare_Loc'],row['Focal_Isolate']) not in snp_dict:
                snp_dict[(row['Focal_Isolate'],row['Focal_Loc'],row['Compare_Isolate'])] = (row['Compare_Loc'],row['Focal_Base'],row['Focal_Direction'],row['Compare_Base'],row['Compare_Direction'])

    end_time = time.time()
    dict_time = end_time - start_time
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write(f"\t- Created loc dictionary in {dict_time:.2f}s\n")
        log.write("\n-------------------------------------------------------\n\n")

    # Compile locs into a single list
    with open(log_file,"a+") as log:
        log.write("Step 6: Compile a single list of locs...\n")
    
    yenta_locs = []
    for isolate in snp_isolates:
        
        start_time = time.time()
        isolate_tuples = [yenta_tuple for yenta_tuple in yenta_tuples if yenta_tuple[0] == isolate]
        isolate_dict = {(k[0],k[1],k[2]):v[0] for k, v in snp_dict.items() if k[0] == isolate}
        isolate_dict.update({(k[2],v[0],k[0]):k[1] for k, v in snp_dict.items() if k[2] == isolate})
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = [executor.submit(process_yenta_tuple, yenta_tuple,isolate_dict) for yenta_tuple in isolate_tuples]
            
            for future in concurrent.futures.as_completed(results):
                yenta_locs.append(future.result())
        end_time = time.time()
        isolate_time = end_time - start_time
        with open(log_file,"a+") as log:
            log.write(f"\t- Processed locs for " + isolate + f" in {isolate_time:.2f}s\n")

    with open(log_file,"a+") as log:
        log.write("\nMerging locs...")
    
    # Create a graph using a dictionary
    start_time = time.time()
    graph = defaultdict(set)

    # Add edges based on shared elements
    for sublist in yenta_locs:
        for pair in combinations(sublist, 2):
            graph[pair[0]].add(pair[1])
            graph[pair[1]].add(pair[0])

    # Perform connected component analysis
    visited = set()
    merged_lists = []

    for element in graph:
        if element not in visited:
            current_list = []
            dfs(element, current_list)
            merged_lists.append([(part.split('_yenta_')[0],part.split('_yenta_')[1]) for part in current_list])
    
    # Remove any sites with purged locs
    merged_not_purged = [loc for loc in merged_lists if not any(any(elem in purged_tuples for elem in tpl) for tpl in loc)]
    removed_locs = [loc for loc in merged_lists if any(any(elem in purged_tuples for elem in tpl) for tpl in loc)]
    end_time = time.time()
    merge_time = end_time - start_time
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write(f"\t- Merged locs in {merge_time:.2f}s\n")
        if len(removed_locs) > 0:
            log.write("\t- Removed "+str(len(removed_locs))+" locs due to QC, " + str(len(merged_not_purged))+" remain.\n")
            log.write("\t- Data regarding which isolates contributed to purged locs can be found in "+snp_dir+"/Purged_SNPs_by_Isolate.tsv\n")
        else:
            log.write("\t- "+str(len(merged_not_purged)) + " sites identified after merging.\n")
        log.write("\n-------------------------------------------------------\n\n")

    # Process data
    full_align_list = []
    df_list = []
    time_list = []
    # Compile locs into a single list
    with open(log_file,"a+") as log:
        log.write("Step 7: Get base data for all locs...\n")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = [executor.submit(getIsolateData,merged_not_purged[i]) for i in range(len(merged_not_purged))]
    
    for future in concurrent.futures.as_completed(results):
        temp_align,temp_df = future.result()
        full_align_list.append(temp_align)
        df_list.append(temp_df)

        
    # Save loc list
    final_snp_df = pd.concat(df_list)
    final_snp_df['Yenta_ID'] = ['Yenta_Loc' + str(i+1).zfill(4) for i in range(len(merged_not_purged))]
    final_snp_df[['Yenta_ID']+snp_isolates].to_csv(snp_dir+"/Loc_List.tsv",index=False,sep="\t")

    # Create alignment
    global final_full_alignment
    final_full_alignment = MultipleSeqAlignment([SeqRecord(Seq(''.join(str(record.seq) for record in records)), id=records[0].id,description="") for records in zip(*full_align_list)])

    # Get pairwise distances
    pairwise_combinations = list(combinations(snp_isolates, 2))
    pairwise_list = [] 
    with concurrent.futures.ThreadPoolExecutor() as executor:
        pairwise_results = [executor.submit(getPairwise, compare,final_full_alignment) for compare in pairwise_combinations]
    
    for future in concurrent.futures.as_completed(pairwise_results):
        pairwise_list.append(future.result())
    
    pairwise_results = pd.concat(pairwise_list)
    pairwise_results['SNP_Difference'] = pairwise_results['Cocalled_Sites'] - pairwise_results['Identical_Sites']
    pairwise_results[["Sample_A","Sample_B","SNP_Difference","Cocalled_Sites"]].to_csv(snp_dir+"/Pairwise_SNP_Distances.tsv",sep="\t",index=False)
    AlignIO.write(final_full_alignment, snp_dir+"/SNP_Alignment.fasta","fasta")

    # Merge pairwise results with previous
    pairwise_results['Combined'] = pairwise_results.apply(lambda row: tuple(sorted([row['Sample_A'], row['Sample_B']])), axis=1)
    pairwise_df['Combined'] = pairwise_df.apply(lambda row: tuple(sorted([row['Query_ID'], row['Reference_ID']])), axis=1)
    merged_pairwise_df = pd.merge(pairwise_df[['Query_ID','Reference_ID','Combined','Yenta_SNPs']],pairwise_results[['Combined','SNP_Difference','Cocalled_Sites']],on='Combined',how='inner')
    merged_pairwise_df.rename(columns={'Yenta_SNPs': 'Raw_Yenta_SNPs', 'SNP_Difference': 'Full_SNP_Difference','Cocalled_Sites':'Full_Cocalled_Sites'}, inplace=True)

    n_data = []
    for record in final_full_alignment:
        sequence_id = record.id
        n_count = Counter(record.seq)['N']
        n_data.append((sequence_id, n_count))
    
    pd.DataFrame(n_data, columns=['Sequence_ID', 'N_Count']).sort_values('N_Count',ascending=False).to_csv(snp_dir+"/SNP_Alignment_Ns.tsv",sep="\t",index=False)
    
    # Filter alignment based on <max_perc_n> per column
    filtered_seqs = []
    filtered_ids = []
    yenta_ids = final_snp_df['Yenta_ID'].tolist()

    for column in range(final_full_alignment.get_alignment_length()):
        n_count = final_full_alignment[:, column].count('N')
        if float(n_count)/len(snp_isolates) <= max_perc_n:
            filtered_seqs.append(final_full_alignment[:, column:column+1])
            filtered_ids = yenta_ids[column]

    # Create a new alignment from the filtered sequences
    if len(filtered_seqs) == final_full_alignment.get_alignment_length():
        merged_pairwise_df.drop('Combined', axis=1).to_csv(snp_dir+"/Merged_Pairwise_Distances.tsv",sep="\t",index=False)
        with open(log_file,"a+") as log:
            log.write("\n\t- The final alignment ("+snp_dir+"/SNP_Alignment.fasta) contains "+str(final_full_alignment.get_alignment_length()) + " sites.\n")
            log.write("\n\t- No sites contained more than " + str(max_perc_n) + "% Ns, so no filtered dataset was generated.\n")
    elif len(filtered_seqs) == 0:
        merged_pairwise_df.drop('Combined', axis=1).to_csv(snp_dir+"/Merged_Pairwise_Distances.tsv",sep="\t",index=False)
        with open(log_file,"a+") as log:
            log.write("\n\t- The final alignment ("+snp_dir+"/SNP_Alignment.fasta) contains "+str(final_full_alignment.get_alignment_length()) + " sites.\n")
            log.write("\n\t- All sites contained more than " + str(max_perc_n) + "% Ns, so no filtered dataset was generated.\n")
    else:
        with open(snp_dir+"/Filtered_Locs.txt","w+") as filt:
            for id in filtered_ids:
                filt.write(id+"\n")
        filtered_alignment = MultipleSeqAlignment([SeqRecord(Seq(''.join(str(record.seq) for record in records)), id=records[0].id,description="") for records in zip(*filtered_seqs)])
        pairwise_list = [] 
        with concurrent.futures.ThreadPoolExecutor() as executor:
            pairwise_results = [executor.submit(getPairwise, compare,filtered_alignment) for compare in pairwise_combinations]
    
        for future in concurrent.futures.as_completed(pairwise_results):
            pairwise_list.append(future.result())
    
        filtered_pairwise_results = pd.concat(pairwise_list)
        filtered_pairwise_results['SNP_Difference'] = filtered_pairwise_results['Cocalled_Sites'] - filtered_pairwise_results['Identical_Sites']
        filtered_pairwise_results[["Sample_A","Sample_B","SNP_Difference","Cocalled_Sites"]].to_csv(snp_dir+"/Filtered_Pairwise_SNP_Distances.tsv",sep="\t",index=False)
        filtered_pairwise_results['Combined'] = filtered_pairwise_results.apply(lambda row: tuple(sorted([row['Sample_A'], row['Sample_B']])), axis=1)
        merged_pairwise_df = pd.merge(merged_pairwise_df,pairwise_results[['Combined','SNP_Difference','Cocalled_Sites']],on='Combined',how='inner')
        merged_pairwise_df.rename(columns={'SNP_Difference': 'Filtered_SNP_Difference','Cocalled_Sites':'Filtered_Cocalled_Sites'}, inplace=True)
        merged_pairwise_df.drop('Combined', axis=1).to_csv(snp_dir+"/Merged_Pairwise_Distances.tsv",sep="\t",index=False)

        AlignIO.write(filtered_alignment, snp_dir+"/Filtered_SNP_Alignment.fasta","fasta")
        n_data = []
        for record in final_full_alignment:
            sequence_id = record.id
            n_count = Counter(record.seq)['N']
            n_data.append((sequence_id, n_count))
    
        pd.DataFrame(n_data, columns=['Sequence_ID', 'N_Count']).sort_values('N_Count',ascending=False).to_csv(snp_dir+"/Filtered_SNP_Alignment_Ns.tsv",sep="\t",index=False)

        with open(log_file,"a+") as log:
            log.write("\n- The raw SNP alignment ("+snp_dir+"/SNP_Alignment.fasta) contains "+str(final_full_alignment.get_alignment_length()) + " sites.\n")
            log.write("\n\t- "+str(final_full_alignment.get_alignment_length() - filtered_alignment.get_alignment_length()) +" sites contained more than " + str(max_perc_n) + "% Ns, so a filtered dataset containing " + str(filtered_alignment.get_alignment_length()) +" was generated at "+snp_dir+"/Filtered_SNP_Alignment.fasta\n")

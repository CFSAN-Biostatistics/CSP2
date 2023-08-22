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

def apply_reverse_complement(row):
    reverse_complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G','a':'t','t':'a','c':'g','g':'c'}
    if row['Focal_Direction'] == -1:
        row['Focal_Base'] = reverse_complements[row['Focal_Base']]
        row['Compare_Base'] = reverse_complements[row['Compare_Base']]
    return row

def getPairwise(compare):

    # Define valid bases
    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']

    # Populate temp alignment with two taxa from compare
    temp_alignment = final_full_alignment[0:0]
    for i in compare:
        temp_alignment.append(SeqRecord(Seq(str(final_full_alignment[isolate_list.index(i)].seq)), id=str(final_full_alignment[isolate_list.index(i)].id)))
    
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
def getIsolateData(yenta_id,tuple_list):
    
    snp_alignment = MultipleSeqAlignment([])

    # Get loc data
    loc_list = []
    for isolate in isolate_list:
        if isolate in [item[0] for item in tuple_list]:
            loc = [item[1] for item in tuple_list if item[0] == isolate][0]
            loc_list.append((isolate,loc))
        else:
            loc_list.append((isolate,"NA"))
    
    loc_df = pd.DataFrame(columns=['Yenta_ID']+isolate_list)
    loc_df.loc[len(loc_df)] = [yenta_id] + [loc[1] for loc in loc_list]

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
    
    # Get base data
    base_list = []
    if len(site_isolates) < 2:
        exit("Fewer than two isolates can't have a SNP?")
    elif len(site_isolates) == 2:
        for isolate in isolate_list:
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
            for isolate in isolate_list:
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

            # For any isolates without a base, go fetch
            while not len(set_bases) == len(site_isolates):
                set_isolates =  [isolate[0] for isolate in set_bases]
                unset_isolates = [isolate for isolate in site_isolates if not isolate in set_isolates]

                print("Site: "+ str(len(site_isolates))+"; Set: "+str(len(set_isolates))+"; Unset: "+str(len(unset_isolates)),flush=True)
                for isolate in unset_isolates:
                    
                    isolate_comparisons = site_df.loc[(site_df['Focal'] == isolate) | (site_df['Compare'] == isolate)]
                    selected_row = isolate_comparisons.loc[(isolate_comparisons['Focal'].isin(set_isolates)) | (isolate_comparisons['Compare'].isin(set_isolates))].iloc[0]
                    if selected_row.shape[0] > 0:
                        selected_focal = selected_row['Focal'][0]
                        selected_focal_base = selected_row['Focal_Base'][0]
                        selected_compare = selected_row['Compare'][0]
                        selected_compare_base = selected_row['Compare_Base'][0]

                        if selected_focal == isolate:
                            compare_base = [base[1] for base in set_bases if base[0] == selected_compare]
                            if selected_row['Compare_Base'][0] == compare_base:
                                set_bases.append((isolate,selected_focal_base))
                            else:
                                set_bases.append((isolate,reverse_complements[selected_focal_base]))
                        else:
                            focal_base = [base[1] for base in set_bases if base[0] == selected_focal]
                            if selected_row['Focal_Base'][0] == focal_base:
                                set_bases.append((isolate,selected_compare_base))
                            else:
                                set_bases.append((isolate,reverse_complements[selected_compare_base]))
    
            if len(missing_isolates) > 0: 
                for isolate in missing_isolates:
                    set_bases.append((isolate,"N"))

    if len(set_bases) == len(isolate_list):
        # Add base data to alignment
        for isolate in isolate_list:
            iso_base = [item[1] for item in set_bases if item[0] == isolate][0]
            snp_alignment.append(SeqRecord(Seq(iso_base), id=str(isolate),description=""))
    else:
        print("Set: "+str(len(set_bases)),flush=True)
        print(set_bases,flush=True)
        exit("Alignment error")
    
    return [snp_alignment,loc_df]

# Read in arguments and set path
mummer_dir = os.path.abspath(str(sys.argv[1]))
snp_dir = os.path.abspath(str(sys.argv[2]))
snp_files = sorted(glob(mummer_dir+"/*MUmmer_SNPs.tsv"))

if len(snp_files) <= 0:
    exit("No files ending in _MUmmer_SNPs.tsv in "+mummer_dir)
else:
    # Read in all SNP files
    print("Reading in SNP files...")
    raw_snp_df = pd.concat((pd.read_table(filename,sep="\t") for filename in snp_files))[['Ref','Ref_Loc','Ref_Base','Ref_Direction','Query','Query_Loc','Query_Base','Query_Direction','Cat']]
    
    if raw_snp_df.shape[0] == 0:
        exit("Files found, but no SNPs detected...")
    else:
        # Get isolate list
        global isolate_list
        isolate_list = sorted(list(set(raw_snp_df.Query) | set(raw_snp_df.Ref)))
        
        # Separate purged locs
        purged_df = raw_snp_df[raw_snp_df.Cat.str.startswith("Purged_")]
        global purged_tuples
        purged_tuples = set(zip(purged_df['Ref'], purged_df['Ref_Loc'])) | set(zip(purged_df['Query'], purged_df['Query_Loc']))

        nonpurged_df = raw_snp_df[~raw_snp_df.Cat.str.startswith("Purged_")]
        if nonpurged_df.shape[0] <= 0:
            exit("No non-purged SNPs detected")
        else:
            # Get Yenta sites
            print("Finding Yenta SNPs...")
            yenta_df = raw_snp_df[raw_snp_df.Cat == "Yenta_SNP"]
            if yenta_df.shape[0] <= 0:
                exit("No Yenta SNPs detected")
            else:
                # Get Yenta loc list
                yenta_tuples = set(zip(yenta_df['Ref'], yenta_df['Ref_Loc'])) | set(zip(yenta_df['Query'], yenta_df['Query_Loc']))

                # Removed purged locs
                yenta_tuples = [tuple_pair for tuple_pair in yenta_tuples if ~(tuple_pair in purged_tuples)]

                print("Getting loc data from isolates...")
                global snp_dict
                snp_dict = defaultdict(lambda:(np.na))
                
                for isolate in isolate_list:
                    
                    print("Gathering Yenta locs for "+isolate+"...")
                    
                    # Create df based on focal isolate
                    isolate_df = nonpurged_df[nonpurged_df.Ref == isolate].rename(columns = {"Ref": "Focal_Isolate", "Ref_Loc": "Focal_Loc", "Ref_Base":"Focal_Base","Ref_Direction":"Focal_Direction","Query":"Compare_Isolate","Query_Loc":"Compare_Loc","Query_Base":"Compare_Base","Query_Direction":"Compare_Direction"})[["Focal_Isolate","Focal_Loc","Focal_Base","Focal_Direction","Compare_Isolate","Compare_Loc","Compare_Base","Compare_Direction"]].append(nonpurged_df[nonpurged_df.Query == isolate].rename(columns = {"Query": "Focal_Isolate", "Query_Loc": "Focal_Loc", "Query_Base":"Focal_Base","Query_Direction":"Focal_Direction","Ref":"Compare_Isolate","Ref_Loc":"Compare_Loc","Ref_Base":"Compare_Base","Ref_Direction":"Compare_Direction"})[["Focal_Isolate","Focal_Loc","Focal_Base","Focal_Direction","Compare_Isolate","Compare_Loc","Compare_Base","Compare_Direction"]])

                    # Add to full SNP dict
                    for index, row in isolate_df.iterrows():
                        snp_dict[(row['Focal_Isolate'],row['Focal_Loc'],row['Compare_Isolate'])] = (row['Compare_Loc'],row['Focal_Base'],row['Focal_Direction'],row['Compare_Base'],row['Compare_Direction'])

                # Compile locs into a single list
                yenta_locs = []
                for isolate in isolate_list:
                    
                    isolate_tuples = [yenta_tuple for yenta_tuple in yenta_tuples if yenta_tuple[0] == isolate]
                    isolate_dict = {(k[0],k[1],k[2]):v[0] for k, v in snp_dict.items() if k[0] == isolate}
                    
                    with concurrent.futures.ThreadPoolExecutor() as executor:
                        results = [executor.submit(process_yenta_tuple, yenta_tuple,isolate_dict) for yenta_tuple in isolate_tuples]
                        start_time = time.time()
                        
                        for future in concurrent.futures.as_completed(results):
                            yenta_locs.append(future.result())
                        
                        total_time = time.time() - start_time
                        print(f"Processed {isolate} tuples in {total_time:.2f}")


                print("Merging Locs...")
                # Create a graph using a dictionary
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
            
            print("Removed "+str(len(removed_locs))+" locs due to QC, " + str(len(merged_not_purged))+" remain...")

            # Process data
            full_align_list = []
            df_list = []

            with concurrent.futures.ThreadPoolExecutor() as executor:
                results = [executor.submit(getIsolateData, 'Yenta_Loc' + str(i+1).zfill(4),merged_not_purged[i]) for i in range(len(merged_not_purged))]
                start_time = time.time()
            
            for future in concurrent.futures.as_completed(results):
                temp_align,temp_df = future.result()
                full_align_list.append(temp_align)
                df_list.append(temp_df)

            final_snp_df = pd.concat(df_list, ignore_index=True)

            global final_full_alignment
            final_full_alignment = MultipleSeqAlignment([SeqRecord(Seq(''.join(str(record.seq) for record in records)), id=records[0].id,description="") for records in zip(*full_align_list)])

            # Get pairwise distances
            pairwise_combinations = list(combinations(isolate_list, 2))
            pairwise_list = [] 
            with concurrent.futures.ThreadPoolExecutor() as executor:
                pairwise_results = [executor.submit(getPairwise, compare) for compare in pairwise_combinations]
            
            for future in concurrent.futures.as_completed(pairwise_results):
                pairwise_list.append(future.result())
            
            pairwise_results = pd.concat(pairwise_list)
            pairwise_results['SNP_Difference'] = pairwise_results['Cocalled_Sites'] - pairwise_results['Identical_Sites']
            pairwise_results[["Sample_A","Sample_B","SNP_Difference","Cocalled_Sites"]].to_csv(snp_dir+"/Pairwise_SNP_Distances.tsv",sep="\t",index=False)
            final_snp_df.to_csv(snp_dir+"/Loc_List.tsv",index=False,sep="\t")
            AlignIO.write(final_full_alignment, snp_dir+"/SNP_Alignment.fasta","fasta")

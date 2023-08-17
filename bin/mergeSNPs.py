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
from Bio import AlignIO, SeqIO
import concurrent.futures
import time
from tqdm import tqdm

def process_yenta_tuple(yenta_tuple):
    isolate, loc = yenta_tuple
    to_add = [isolate + "_yenta_" + loc] + [key[2] + "_yenta_" + value for key, value in snp_dict.items() if (key[0], key[1]) == (isolate, loc)] + [key[0] + "_yenta_" + key[1] for key, value in snp_dict.items() if (key[2], value) == (isolate, loc)]
    return list(set(sorted(to_add, key=lambda x: x[0])))

def dfs(node, current_list):
    visited.add(node)
    current_list.append(node)
    for neighbor in graph[node]:
        if neighbor not in visited:
            dfs(neighbor, current_list)

# Function to return alignments and coordinates for Yenta SNPs
def getIsolateData(yenta_id,tuple_list):
    
    # Create alignment
    full_snp_alignment = MultipleSeqAlignment([])
    loc_list = []

    # Check if any locs are on an edge
    if any(tuple_pair in edge_tuples for tuple_pair in tuple_list):
        full_data = False
    else:
        full_data = True
    
    # Get loc + base (if present)
    for isolate in isolate_list:
        if isolate in [item[0] for item in tuple_list]:
            loc = [item[1] for item in tuple_list if item[0] == isolate][0]
        else:
            loc = "NA"

        loc_list.append(loc)
        full_snp_alignment.append(SeqRecord(Seq(base_dict[(isolate,loc)]), id=str(isolate)))
    
    loc_df = pd.DataFrame(columns=['Yenta_ID','Full_Dataset']+isolate_list)
    loc_df.loc[len(loc_df)] = [yenta_id] + [full_data] + loc_list

    if full_data:
        return [full_snp_alignment,full_snp_alignment,loc_df]
    else:
        return [full_snp_alignment,full_snp_alignment[:,0:0],loc_df]

# Read in arguments and set path
mummer_dir = os.path.abspath(str(sys.argv[1]))
snp_files = sorted(glob(mummer_dir+"/*MUmmer_SNPs.tsv"))

if len(snp_files) <= 0:
    exit("No files ending in _MUmmer_SNPs.tsv in "+mummer_dir)
else:
    # Read in all SNP files
    print("Reading in SNP files...")
    raw_snp_df = pd.concat((pd.read_table(filename,sep="\t") for filename in snp_files))[['Ref','Ref_Loc','Ref_Base','Query','Query_Loc','Query_Base','Cat']]
    
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
                # Set aside edge locs
                edge_df = raw_snp_df[raw_snp_df.Cat == "Filtered_Edge"]
                global edge_tuples
                edge_tuples = set(zip(edge_df['Ref'], edge_df['Ref_Loc'])) | set(zip(edge_df['Query'], edge_df['Query_Loc']))

                # Get Yenta loc list
                yenta_tuples = set(zip(yenta_df['Ref'], yenta_df['Ref_Loc'])) | set(zip(yenta_df['Query'], yenta_df['Query_Loc']))

                # Removed purged locs
                yenta_tuples = [tuple_pair for tuple_pair in yenta_tuples if ~(tuple_pair in purged_tuples)]
                
                global base_dict
                base_dict = defaultdict(lambda:"N")

                global snp_dict
                snp_dict = defaultdict(lambda:(np.na))
                print("Getting isolate data...")

                for isolate in isolate_list:
                    
                    # Create df based on focal isolate
                    isolate_df = nonpurged_df[nonpurged_df.Ref == isolate].rename(columns = {"Ref": "Focal_Isolate", "Ref_Loc": "Focal_Loc", "Ref_Base":"Focal_Base","Query":"Compare_Isolate","Query_Loc":"Compare_Loc","Query_Base":"Compare_Base"})[["Focal_Isolate","Focal_Loc","Focal_Base","Compare_Isolate","Compare_Loc","Compare_Base"]].append(nonpurged_df[nonpurged_df.Query == isolate].rename(columns = {"Query": "Focal_Isolate", "Query_Loc": "Focal_Loc", "Query_Base":"Focal_Base","Ref":"Compare_Isolate","Ref_Loc":"Compare_Loc","Ref_Base":"Compare_Base"})[["Focal_Isolate","Focal_Loc","Focal_Base","Compare_Isolate","Compare_Loc","Compare_Base"]])
                    
                    # Create base dict after removing duplicates
                    for index, row in isolate_df.drop_duplicates(subset = ['Focal_Isolate','Focal_Loc'],keep='first').iterrows():
                        base_dict[(row['Focal_Isolate'],row['Focal_Loc'])] = row['Focal_Base']

                    # Create SNP dict
                    for index, row in isolate_df.iterrows():
                        snp_dict[(row['Focal_Isolate'],row['Focal_Loc'],row['Compare_Isolate'])] = row['Compare_Loc']

                print("Compiling locs...")            
                yenta_locs = []
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    results = [executor.submit(process_yenta_tuple, yenta_tuple) for yenta_tuple in yenta_tuples]
                    count = 0
                    start_time = time.time()
                    
                    # Use tqdm to display the progress bar
                    for future in tqdm(concurrent.futures.as_completed(results), total=len(results)):
                        yenta_locs.append(future.result())
                        count += 1
                    
                    total_time = time.time() - start_time
                print(f"Processed a total of {count} tuples in {total_time:.2f} seconds")

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

            # Process data
            full_align_list = []
            filtered_align_list = []
            df_list = []
            for i in range(len(merged_lists)):
                temp_full,temp_filtered,temp_df = getIsolateData('Yenta_Loc' + str(i+1).zfill(4),merged_lists[i])
                full_align_list.append(temp_full)
                filtered_align_list.append(temp_filtered)
                df_list.append(temp_df)

            final_snp_df = pd.concat(df_list, ignore_index=True)
            final_full_alignment = MultipleSeqAlignment([SeqRecord(Seq(''.join(str(record.seq) for record in records)), id=records[0].id) for records in zip(*full_align_list)])
            final_filtered_alignment = MultipleSeqAlignment([SeqRecord(Seq(''.join(str(record.seq) for record in records)), id=records[0].id) for records in zip(*filtered_align_list)])

            print(final_snp_df)
            print(final_full_alignment)
            print(final_filtered_alignment)

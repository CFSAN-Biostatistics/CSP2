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
from pybedtools import BedTool
import warnings
warnings.filterwarnings("ignore")

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
        bed_df['start'] = bed_df['start'] - 1
        bed_file = BedTool.from_dataframe(bed_df[['chrom','start','end']]).sort()
    return bed_file

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

def processLoc(loc_df,ref_loc):
    snp_alignment = MultipleSeqAlignment([])
    ref_base = loc_df['Ref_Base'].iloc[0]

    for isolate in snp_isolates:
        if isolate == ref_isolate:
            snp_alignment.append(SeqRecord(Seq(ref_base), id=str(isolate),description=""))
        elif isolate in loc_df['Query'].values:
            snp_alignment.append(SeqRecord(Seq(loc_df.loc[loc_df['Query'] == isolate, 'Query_Base'].values[0]), id=str(isolate),description=""))
        elif ref_loc in coords_dict[isolate]:
            snp_alignment.append(SeqRecord(Seq(ref_base), id=str(isolate),description=""))
        else:
            snp_alignment.append(SeqRecord(Seq("N"), id=str(isolate),description=""))

    # Query coords_dict to figure out which isolates have coverage
    return [snp_alignment,ref_loc]

def parseMUmmerCoords(yenta_bed,coords_dir,query_id,ref_id,perc_iden,min_len):
    coords_file = pd.read_csv(coords_dir+"/"+query_id+"_vs_"+ref_id+".1coords",sep="\t",index_col=False,
    names=['Ref_Start','Ref_End','Query_Start','Query_End',
    'Ref_Aligned','Query_Aligned','Perc_Iden',
    'Ref_Length','Query_Length','Ref_Cov',
    'Query_Cov','Ref_Contig','Query_Contig'])

    coords_file = coords_file.loc[(coords_file['Ref_Aligned'] >= min_len) & (coords_file['Perc_Iden'] >= perc_iden)]
    ref_bed = makeBED(coords_file[['Ref_Contig','Ref_Start','Ref_End']])
    try:
        yenta_intersect = ref_bed.intersect(yenta_bed,wa=True,wb=True,header=False).to_dataframe().iloc[:, 6].tolist()
        coords_dict.update({query_id:yenta_intersect})
    except:
        sys.exit("Intersect failure")
    
# Read in arguments
output_dir = os.path.abspath(sys.argv[1])
align_cov = float(sys.argv[2])
min_perc_iden = float(sys.argv[3])
min_length = float(sys.argv[4])
max_perc_n = float(sys.argv[5])
ref_isolate = str(sys.argv[6])

# Set paths
mummer_dir = output_dir + "/MUmmer_Output"
coords_dir = mummer_dir + "/Raw/1coords"
snp_dir = output_dir + "/SNP_Analysis"
ref_directory = snp_dir+"/"+ref_isolate

# Create logfile
global log_file
log_file = ref_directory+"/SNP_Analysis_"+ref_isolate+".log"

# Create reference folder
if not os.path.exists(ref_directory):
    os.makedirs(ref_directory)
    with open(log_file,"w+") as log:
        log.write("-------------------------------------------------------\n")
        log.write("Yenta SNP Analysis\n")
        log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
        log.write("-------------------------------------------------------\n\n")
else:
    sys.exit(ref_directory+" already exists...")

# Read in SNP files
snp_files = sorted(glob(mummer_dir+"/*_vs_"+ref_isolate+"_MUmmer_SNPs.tsv"))
if len(snp_files) == 0:
    with open(log_file,"a+") as log:
        log.write("ERROR: No files ending in _vs_"+ref_isolate+"_MUmmer_SNPs.tsv in "+mummer_dir+"\n")
    sys.exit("No files ending in _vs_"+ref_isolate+"_MUmmer_SNPs.tsv in "+mummer_dir)

else:
    # Read in pairwise file
    with open(log_file,"a+") as log:
        log.write("Step 1: Reading in pairwise distance file...")
    try:
        pairwise_df = pd.read_table(output_dir+"/Raw_Pairwise_Distances.tsv",sep="\t")
        pairwise_df = pairwise_df[pairwise_df["Reference_ID"] == ref_isolate]
        full_isolate_list = sorted(list(set(pairwise_df.Query_ID) | set(pairwise_df.Reference_ID)))
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write("\t- Found pairwise file ("+output_dir+"/Raw_Pairwise_Distances.tsv), which contains data on " + str(len(full_isolate_list)) + " isolates...\n")
            log.write("\t- Pairwise file contains "+str(pairwise_df.shape[0])+" rows...\n")
            log.write("\n-------------------------------------------------------\n\n")
    except:
        with open(log_file,"a+") as log:
            log.write("ERROR: Cannot find pairwise file at "+output_dir+"/Raw_Pairwise_Distances.tsv\n")
        sys.exit("Cannot find pairwise file at "+output_dir+"/Raw_Pairwise_Distances.tsv")

    # Get median coverage stats
    with open(log_file,"a+") as log:
        log.write("Step 2: Removing isolates with a mean genome coverage < "+str(align_cov) + "% ...")
    median_cov_df = pd.DataFrame({'ID':  pd.concat([pairwise_df['Query_ID'], pairwise_df['Reference_ID']]), 'Coverage': pd.concat([pairwise_df['Percent_Reference_Covered'], pairwise_df['Percent_Query_Covered']])}).groupby('ID')['Coverage'].median().reset_index()

    # Extract isolates with a median percent coverage < min_median_coverage
    low_cov_ids = median_cov_df[median_cov_df['Coverage'] < align_cov]['ID'].tolist()
    
    global snp_isolates
    snp_isolates = median_cov_df[median_cov_df['Coverage'] >= align_cov]['ID'].tolist()
    
    if not ref_isolate in snp_isolates:
        with open(log_file,"a+") as log:
            log.write("Reference isolate "+ref_isolate+" has less than the required genomic coverage to proceed...\n")
        sys.exit("Reference isolate "+ref_isolate+" has less than the required genomic coverage to proceed...")

    query_isolates = [isolate for isolate in snp_isolates if not isolate == ref_isolate]

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
        raw_snp_df = pd.concat((pd.read_table(filename,sep="\t") for filename in filtered_snp_files))
        end_time = time.time()
        snp_time = end_time - start_time
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write(f"\t- Read in "+str(len(filtered_snp_files)) + f" SNP files in {snp_time:.2f}s\n")
            log.write("\n-------------------------------------------------------\n\n")

    except:
        with open(log_file,"a+") as log:
            log.write("\n\t- ERROR: Cannot load SNP files in "+mummer_dir+"\n")
        sys.exit("Cannot load SNP files in "+mummer_dir)     
    
    if raw_snp_df.shape[0] == 0:
        with open(log_file,"a+") as log:
            log.write("\n\t- SNP files located, but no SNPs present. Cannot proceed...\n")

    else:    
        
        # Separate purged locs
        with open(log_file,"a+") as log:
            log.write("Step 4: Removing purged sites...")
        purged_df = raw_snp_df[raw_snp_df.Cat.str.startswith("Purged_")]
        purged_locs = purged_df['Ref_Loc'].values
        if purged_df.shape[0] > 0:
            purged_counts = purged_df['Query'].value_counts()
            pd.DataFrame({'ID': purged_counts.index, 'Count': purged_counts.values}).sort_values(by='Count', ascending=False).to_csv(ref_directory+"/Purged_SNPs_by_Isolate.tsv",index=False,sep="\t")
        
        # Process non-purged locs
        nonpurged_df = raw_snp_df[~raw_snp_df.Cat.str.startswith("Purged_")] 
        
        # Check that all ref directions are positive
        if -1 in nonpurged_df['Ref_Direction'].values:
            with open(log_file,"a+") as log:
                log.write("\t- Reference sequences should all be in the positive direction...\n")
            sys.exit("Reference sequences should all be in the positive direction...")

        with open(log_file,"a+") as log:
            log.write("Done!\n")
        
        if nonpurged_df.shape[0] == 0:
            with open(log_file,"a+") as log:
                log.write("\t- All SNPs purged...\n")
        
        else:

            # Get Yenta sites
            yenta_df = raw_snp_df[raw_snp_df.Cat == "Yenta_SNP"]
            yenta_locs = np.unique(yenta_df['Ref_Loc'].values)
            removed_locs = [loc for loc in yenta_locs if loc in purged_locs]
            #yenta_locs = [loc for loc in yenta_locs if not loc in purged_locs]
                        
            if len(yenta_locs) == 0:
                with open(log_file,"a+") as log:
                    log.write("\t- SNP files located, but no Yenta SNPs present. Cannot proceed...\n")

            else:
                if len(removed_locs) > 0:
                    with open(log_file,"a+") as log:
                        log.write("\t- Removed "+str(len(removed_locs))+" locs due to QC, " + str(len(yenta_locs))+" remain.\n")
                        log.write("\t- Data regarding which isolates contributed to purged locs can be found in "+ref_directory+"/Purged_SNPs_by_Isolate.tsv\n")
                        log.write("\n-------------------------------------------------------\n\n")
                else:
                    with open(log_file,"a+") as log:
                        log.write("\t- "+str(len(yenta_locs)) + " sites identified after merging.\n")
                        log.write("\n-------------------------------------------------------\n\n")

                # Create BED file for Yenta locs
                yenta_bed = makeBED(pd.DataFrame([item.split('/') for item in yenta_locs], columns=['Ref_Contig','Ref_End']))

                # Read in coordinate files
                with open(log_file,"a+") as log:
                    log.write("Step 5: Processing coordinate files...")
                
                try:
                    start_time = time.time()

                    global coords_dict
                    coords_dict = defaultdict(lambda:[])
                    with concurrent.futures.ThreadPoolExecutor() as executor:
                        results = [executor.submit(parseMUmmerCoords,yenta_bed,coords_dir,isolate,ref_isolate,min_perc_iden,min_length) for isolate in query_isolates]
                    end_time = time.time()
                    coords_time = end_time - start_time
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"\t- Read in "+str(len(query_isolates)) + f" 1coords files in {coords_time:.2f}s\n")
                        log.write("\n-------------------------------------------------------\n\n")
                except:
                    with open(log_file,"a+") as log:
                        log.write("\n\t- ERROR: Cannot load coords files in "+coords_dir+"\n")
                    sys.exit("Cannot load coords files in "+coords_dir)

                # Process loc data
                with open(log_file,"a+") as log:
                    log.write("Step 6: Processing "+str(len(yenta_locs))+" loc(s)...")

                full_align_list = []
                loc_list = []
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    results = [executor.submit(processLoc, nonpurged_df[nonpurged_df.Ref_Loc == loc],loc) for loc in yenta_locs]

                for future in concurrent.futures.as_completed(results):
                    temp_align,loc = future.result()
                    full_align_list.append(temp_align)
                    loc_list.append(loc)
                
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
                pairwise_results[["Sample_A","Sample_B","SNP_Difference","Cocalled_Sites"]].to_csv(ref_directory+"/Pairwise_SNP_Distances.tsv",sep="\t",index=False)
                AlignIO.write(final_full_alignment, ref_directory+"/SNP_Alignment.fasta","fasta")

                n_data = []
                for record in final_full_alignment:
                    sequence_id = record.id
                    n_count = Counter(record.seq)['N']
                    n_data.append((sequence_id, n_count))
                
                pd.DataFrame(n_data, columns=['Sequence_ID', 'N_Count']).sort_values('N_Count',ascending=False).to_csv(ref_directory+"/SNP_Alignment_Ns.tsv",sep="\t",index=False)
                
                # Filter alignment based on <max_perc_n> per column
                filtered_seqs = []
                filtered_ids = []

                for column in range(final_full_alignment.get_alignment_length()):
                    n_count = final_full_alignment[:, column].count('N')
                    if float(n_count)/len(snp_isolates) <= max_perc_n:
                        filtered_seqs.append(final_full_alignment[:, column:column+1])
                        filtered_ids = loc_list[column]

                # Create a new alignment from the filtered sequences
                if len(filtered_seqs) == final_full_alignment.get_alignment_length():
                    with open(log_file,"a+") as log:
                        log.write("\n\t- The final alignment ("+ref_directory+"/SNP_Alignment.fasta) contains "+str(final_full_alignment.get_alignment_length()) + " sites.\n")
                        log.write("\n\t- No sites contained more than " + str(max_perc_n) + "% Ns, so no filtered dataset was generated.\n")
                elif len(filtered_seqs) == 0:
                    with open(log_file,"a+") as log:
                        log.write("\n\t- The final alignment ("+ref_directory+"/SNP_Alignment.fasta) contains "+str(final_full_alignment.get_alignment_length()) + " sites.\n")
                        log.write("\n\t- All sites contained more than " + str(max_perc_n) + "% Ns, so no filtered dataset was generated.\n")
                else:
                    with open(ref_directory+"/Filtered_Locs.txt","w+") as filt:
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
                    filtered_pairwise_results[["Sample_A","Sample_B","SNP_Difference","Cocalled_Sites"]].to_csv(ref_directory+"/Filtered_Pairwise_SNP_Distances.tsv",sep="\t",index=False)
                    filtered_pairwise_results['Combined'] = filtered_pairwise_results.apply(lambda row: tuple(sorted([row['Sample_A'], row['Sample_B']])), axis=1)

                    AlignIO.write(filtered_alignment, ref_directory+"/Filtered_SNP_Alignment.fasta","fasta")
                    n_data = []
                    for record in filtered_alignment:
                        sequence_id = record.id
                        n_count = Counter(record.seq)['N']
                        n_data.append((sequence_id, n_count))
                
                    pd.DataFrame(n_data, columns=['Sequence_ID', 'N_Count']).sort_values('N_Count',ascending=False).to_csv(ref_directory+"/Filtered_SNP_Alignment_Ns.tsv",sep="\t",index=False)

                    with open(log_file,"a+") as log:
                        log.write("\n- The raw SNP alignment ("+ref_directory+"/SNP_Alignment.fasta) contains "+str(final_full_alignment.get_alignment_length()) + " sites.\n")
                        log.write("\n\t- "+str(final_full_alignment.get_alignment_length() - filtered_alignment.get_alignment_length()) +" sites contained more than " + str(max_perc_n) + "% Ns, so a filtered dataset containing " + str(filtered_alignment.get_alignment_length()) +" was generated at "+ref_directory+"/Filtered_SNP_Alignment.fasta\n")

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

#def processLoc(loc_df,ref_loc):
#    print(ref_loc,flush=True)
#    snp_alignment = MultipleSeqAlignment([])
#    ref_base = loc_df['Ref_Base'].iloc[0]

#    for isolate in snp_isolates:
#        if isolate == ref_isolate:
#            snp_alignment.append(SeqRecord(Seq(ref_base), id=str(isolate),description=""))
#        elif isolate in loc_df['Query'].values:
#            snp_alignment.append(SeqRecord(Seq(loc_df.loc[loc_df['Query'] == isolate, 'Query_Base'].values[0]), id=str(isolate),description=""))
#        elif ref_loc in coords_dict[isolate]:
#            if ref_loc in purged_dict[isolate]:
#                snp_alignment.append(SeqRecord(Seq("N"), id=str(isolate),description=""))
#            else:
#                snp_alignment.append(SeqRecord(Seq(ref_base), id=str(isolate),description=""))
#        else:
#            snp_alignment.append(SeqRecord(Seq("?"), id=str(isolate),description=""))

#    # Query coords_dict to figure out which isolates have coverage
#    return [snp_alignment,ref_loc]

def processLoc(loc_df, ref_loc):
    print(ref_loc, flush=True)
    
    ref_base = loc_df['Ref_Base'].iloc[0]

    # Create a list of SeqRecord objects
    snp_records = []

    for isolate in snp_isolates:
        if isolate == ref_isolate:
            snp_records.append(SeqRecord(Seq(ref_base), id=str(isolate), description=""))
        else:
            query_match = loc_df[loc_df['Query'] == isolate]
            if not query_match.empty:
                snp_records.append(SeqRecord(Seq(query_match['Query_Base'].iloc[0]), id=str(isolate), description=""))
            elif ref_loc in purged_dict[isolate]:
                snp_records.append(SeqRecord(Seq("N"), id=str(isolate), description=""))
            elif ref_loc in coords_dict[isolate]:
                snp_records.append(SeqRecord(Seq(ref_base), id=str(isolate), description=""))
            else:
                snp_records.append(SeqRecord(Seq("?"), id=str(isolate), description=""))

    snp_alignment = MultipleSeqAlignment(snp_records)
    return [snp_alignment, ref_loc]

def parseMUmmerCoords(yenta_bed,coords_dir,query_id,ref_id,perc_iden,min_len):
    coords_file = pd.read_csv(coords_dir+"/"+query_id+"_vs_"+ref_id+".1coords",sep="\t",index_col=False,
    names=['Ref_Start','Ref_End','Query_Start','Query_End',
    'Ref_Aligned','Query_Aligned','Perc_Iden',
    'Ref_Length','Query_Length','Ref_Cov',
    'Query_Cov','Ref_Contig','Query_Contig'])

    coords_file = coords_file.loc[(coords_file['Ref_Aligned'] >= min_len) & (coords_file['Perc_Iden'] >= perc_iden) & (coords_file['Ref_Contig'].isin(yenta_contigs))]
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
    sys.exit(0)
elif len(snp_isolates) == 2:
    with open(log_file,"a+") as log:
        log.write("\t- Two isolates remain after filtering for low coverage, results from pairwise will not change.\n")
        sys.exit(0)
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
    sys.exit(0)

if -1 in raw_snp_df['Ref_Direction'].values:
    with open(log_file,"a+") as log:
        log.write("\t- Reference sequences should all be in the positive direction...\n")
    sys.exit("Reference sequences should all be in the positive direction...")
    
with open(log_file,"a+") as log:
    log.write("Step 4: Collecting Yenta SNPs...")  

# Get Yenta sites
yenta_df = raw_snp_df[raw_snp_df.Cat == "Yenta_SNP"]
yenta_locs = np.unique(yenta_df['Ref_Loc'].values)

global yenta_contigs
yenta_contigs = np.unique([item.split("/")[0] for item in yenta_locs])

yenta_count = len(yenta_locs)

if len(yenta_locs) == 0:
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write("\t- SNP files located, but no Yenta SNPs present. Cannot proceed...\n")
    sys.exit(0)

with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write("\t- "+str(yenta_count) + " SNP positions identified\n")
    log.write("\n-------------------------------------------------------\n\n")

# Separate purged locs
with open(log_file,"a+") as log:
    log.write("Step 5: Identifying purged sites...")

global purged_dict
purged_dict = defaultdict(lambda:[])

purged_df = raw_snp_df[raw_snp_df.Cat.str.startswith("Purged_") & raw_snp_df['Ref_Loc'].isin(yenta_locs)]
if purged_df.shape[0] > 0:
    for query, ref_loc in zip(purged_df['Query'], purged_df['Ref_Loc']):
        if query in purged_dict:
            purged_dict[query].append(ref_loc)
        else:
            purged_dict[query] = [ref_loc]

with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write("\n-------------------------------------------------------\n\n")

with open(log_file,"a+") as log:
    log.write("Step 6: Creating Yenta BED...") 

try:
    # Create BED file for Yenta locs
    yenta_bed = makeBED(pd.DataFrame([item.split('/') for item in yenta_locs], columns=['Ref_Contig','Ref_End']))

    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write("\n-------------------------------------------------------\n\n")
except:
        with open(log_file,"a+") as log:
            log.write("\n\t- Error: Cannot create Yenta BED file...\n")
        sys.exit("Cannot create Yenta BED file")

# Read in coordinate files
with open(log_file,"a+") as log:
    log.write("Step 7: Processing coordinate files...")
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
    log.write("Step 8: Processing "+str(yenta_count)+" loc(s)...")
start_time = time.time()
full_align_list = []
loc_list = []
nonpurged_df = raw_snp_df[~(raw_snp_df.Cat.str.startswith("Purged_") & raw_snp_df['Ref_Loc'].isin(yenta_locs))]
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = [executor.submit(processLoc, nonpurged_df[nonpurged_df.Ref_Loc == loc],loc) for loc in yenta_locs]

iteration_counter = 0
for future in concurrent.futures.as_completed(results):
    temp_align,loc = future.result()
    full_align_list.append(temp_align)
    loc_list.append(loc)
    iteration_counter += 1

    # Check if it's a multiple of 1000 and print a message
    if iteration_counter % 1000 == 0:
        iteration_counter = 0
        temp_time = time.time()
        thous_time = start_time - temp_time
        print(f"- Processed {iteration_counter} in {thous_time:.2f}s\n")

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
end_time = time.time()
align_time = end_time - start_time
with open(log_file,"a+") as log:
    log.write("Done!\n")
    log.write("\t- Generated the alignment file "+ref_directory+f"/SNP_Alignment.fasta in {align_time:.2f}s\n")
    log.write("\n-------------------------------------------------------\n\n")

# Process loc data
with open(log_file,"a+") as log:
    log.write("Step 9: Saving LocList and characterizing purged and missing data...")

n_data = []
uncovered_data = []

for record in final_full_alignment:
    sequence_id = record.id
    n_count = Counter(record.seq)['N']
    q_count = Counter(record.seq)['?']
    n_data.append((sequence_id, n_count))
    uncovered_data.append((sequence_id, q_count))

n_df = pd.DataFrame(n_data, columns=['Isolate', 'Purged_Count'])
q_df = pd.DataFrame(uncovered_data, columns=['Isolate', 'Missing_Count'])

pd.merge(n_df, q_df, on='Isolate', how='outer').sort_values('Missing_Count',ascending=False).to_csv(ref_directory+"/Missing_Purged_SNPs_by_Isolate.tsv",sep="\t",index=False)

loc_data = []
with open(ref_directory+"/Loc_List.txt","w+") as loc_file:
    for column in range(final_full_alignment.get_alignment_length()):
        n_count = final_full_alignment[:, column].count('N')
        q_count = final_full_alignment[:, column].count('?')
        loc_data.append((loc_list[column],n_count,q_count))
        loc_file.write(loc_list[column]+"\n")

site_df = pd.DataFrame(loc_data, columns=['SNP_Loc','Purged_Count','Missing_Count']).sort_values('Missing_Count',ascending=False).to_csv(ref_directory+"/Missing_Purged_SNPs_by_Loc.tsv",sep="\t",index=False)

with open(log_file,"a+") as log:
    log.write("Done!")
#!/usr/bin/env python3

import sys
import os
import pandas as pd
import re
from pybedtools import BedTool
import warnings
import ast
import numpy as np
import hashlib
from Bio import SeqIO
warnings.filterwarnings("ignore")

def makeBED(bed_df):
    if bed_df.shape[1] == 2:
        bed_df.columns = ['Ref_Contig','Ref_End']
        bed_df['Ref_End'] = bed_df['Ref_End'].astype(int)
        bed_df['Ref_Start'] = bed_df['Ref_End'] - 1
        bed_file = BedTool.from_dataframe(bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort()
    elif bed_df.shape[1] == 3:
        bed_df.columns = ['Ref_Contig','Ref_End','Metadata']
        bed_df['Ref_End'] = bed_df['Ref_End'].astype(int)
        bed_df['Ref_Start'] = bed_df['Ref_End'] - 1
        bed_file = BedTool.from_dataframe(bed_df[['Ref_Contig','Ref_Start','Ref_End','Metadata']]).sort()
    else:
        sys.exit("bed_df should have 2 or 3 columns")
    return bed_file

    bed_df.columns = ['Ref_Contig','Ref_Start','Ref_End']
    bed_df['Ref_Start'] = bed_df['Ref_Start'] - 1
    bed_file = BedTool.from_dataframe(bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort()
    return bed_file

def parseMUmmerReport(mum_report_dir,report_id):
    report_data = []
    with open(mum_report_dir+"/"+report_id+".report", 'r') as report_file:

        for line in report_file.readlines():
            split_line = re.split(r'\s+', line.strip())
            if len(split_line) == 3:
                report_data.append(
                    {
                    'Measure': split_line[0].split("(", 1)[0],
                    'Ref_Value':split_line[1].split("(", 1)[0],
                    'Query_Value':split_line[2].split("(", 1)[0]
                    })

    report_data = pd.DataFrame(report_data)
    report_data = report_data[report_data['Measure'].isin(['TotalSeqs','TotalBases','AlignedBases','TotalGSNPs','TotalGIndels'])]
    report_data = report_data.drop_duplicates(subset=['Measure']) # Drop M-to-M

    report_data = report_data.astype({"Ref_Value":"int","Query_Value":"int"})

    # Fetch assembly data
    ref_bases = report_data[report_data['Measure'] == 'TotalBases']['Ref_Value'].iloc[0]
    query_bases = report_data[report_data['Measure'] == 'TotalBases']['Query_Value'].iloc[0]

    # Fetch alignment data
    ref_aligned = report_data[report_data['Measure'] == 'AlignedBases']['Ref_Value'].iloc[0]
    query_aligned = report_data[report_data['Measure'] == 'AlignedBases']['Query_Value'].iloc[0]

    percent_ref_aligned = 100*(float(ref_aligned)/ref_bases)
    percent_query_aligned = 100*(float(query_aligned)/query_bases)
    
    return [ref_bases,percent_ref_aligned,query_bases,percent_query_aligned]

def parseMUmmerCoords(mum_coords_dir,report_id,min_iden,min_len):
    
    return_columns = ['Ref_Contig','Ref_Length','Ref_Start','Ref_End','Ref_Aligned','Query_Contig','Query_Length','Query_Start','Query_End','Query_Aligned','Perc_Iden']
    coords_file = pd.read_csv(mum_coords_dir+"/"+report_id+".1coords",sep="\t",index_col=False,
    names=['Ref_Start','Ref_End','Query_Start','Query_End',
    'Ref_Aligned','Query_Aligned','Perc_Iden',
    'Ref_Length','Query_Length','Ref_Cov',
    'Query_Cov','Ref_Contig','Query_Contig'])
    
    bad_coords_file = coords_file[(coords_file.Ref_Aligned < min_len) | (coords_file.Perc_Iden < min_iden)]
    coords_file = coords_file[(coords_file.Ref_Aligned >= min_len) & (coords_file.Perc_Iden >= min_iden)]

    # Fix start coordinates and create a BED file
    if coords_file.shape[0] > 0:
        coords_file['Ref_Start'] = (coords_file['Ref_Start'].astype(int) - 1)
        coords_file['Query_Start'] = (coords_file['Query_Start'].astype(int) - 1)
        ref_bed = BedTool.from_dataframe(coords_file[['Ref_Contig','Ref_Start','Ref_End']]).sort().merge()
    else:
        ref_bed = BedTool([])
    
    # Fix start coordinates
    if bad_coords_file.shape[0] > 0:
        bad_coords_file['Ref_Start'] = (bad_coords_file['Ref_Start'].astype(int) - 1)
        bad_coords_file['Query_Start'] = (bad_coords_file['Query_Start'].astype(int) - 1)

    return [coords_file[return_columns],bad_coords_file[return_columns],ref_bed]
    
def parseMUmmerSNPs(mum_snps_dir,report_id):

    return_columns = ['Ref_Contig','Ref_Pos','Ref_Loc','Query_Contig','Query_Pos','Query_Loc','Dist_to_Ref_End','Dist_to_Query_End','Ref_Base','Query_Base','Ref_Direction','Query_Direction']

    snp_file = pd.read_csv(mum_snps_dir+"/"+report_id+".snps",sep="\t",index_col=False,
        names=['Ref_Pos','Ref_Base','Query_Base','Query_Pos',
        'SNP_Buffer','Dist_to_End','Ref_Length','Query_Length',
        'Ref_Direction','Query_Direction','Ref_Contig','Query_Contig'])

    # Ignore SNPs where the reference has an indel
    snp_file = snp_file[snp_file['Ref_Base'] != "."]

    # If no nucleotide SNPs are present, return empty df
    if snp_file.shape[0] == 0:
        return snp_file
    else:
        # Get minimum dist_to_end for ref and query
        snp_file['Dist_to_Ref_End'] = [min([x,y]) for x,y in zip(snp_file['Ref_Pos'],snp_file['Ref_Length'] - snp_file['Ref_Pos'])]
        snp_file['Dist_to_Query_End'] = [min([x,y]) for x,y in zip(snp_file['Query_Pos'],snp_file['Query_Length'] - snp_file['Query_Pos'])]

        # Add Loc data
        snp_file['Ref_Loc'] = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(snp_file.Ref_Contig, snp_file.Ref_Pos))]
        snp_file['Query_Loc'] = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(snp_file.Query_Contig, snp_file.Query_Pos))]    

        return snp_file[return_columns]

def filterSNPs(snp_file,coords_file,bad_coords_file,density_windows,max_snps,ref_edge,query_edge):

    # Get raw SNP count
    raw_snp_count = snp_file.shape[0]

    # SNPs from bad alignments (ignore downstream)
    bad_snp_coords = pd.merge(snp_file,bad_coords_file,how='left')
    snps_fail_idenlen = bad_snp_coords[~bad_snp_coords.isnull().any(1)]
    snps_fail_idenlen = snps_fail_idenlen[snps_fail_idenlen.apply(lambda x: x['Ref_Start'] <= x['Ref_Pos'] <= x['Ref_End'], axis=1)]
    snps_fail_idenlen['Cat'] = "Purged_Alignment"

    # Process duplicate locs from bad coords
    loc_tally = snps_fail_idenlen['Ref_Loc'].value_counts()
    solo_locs = loc_tally[loc_tally == 1].index.tolist()
    
    # Grab locs covered by 1 bad alignment
    dup_fail_idenlen = snps_fail_idenlen[snps_fail_idenlen.Ref_Loc.isin(solo_locs)]
    
    # Process locs covered by 2+ alignments
    dup_query = loc_tally[loc_tally > 1]
    dup_locs = dup_query.index.tolist()

    if len(dup_locs) > 0:
        snps_check_dup = snps_fail_idenlen[snps_fail_idenlen.Ref_Loc.isin(dup_locs)]
        for loc in dup_locs:
            # Find all query overlaps with each reference locus
            dup_check_df = snps_check_dup[snps_check_dup['Ref_Loc'] == loc]
            overlap_count = dup_check_df.shape[0]
            dup_check_snps = snp_file[snp_file['Ref_Loc'] == loc]
            dup_snp_count = dup_check_snps.shape[0]

            # If the same number of SNPs and overlaps appear:
            if overlap_count == dup_snp_count:
                dup_fail_idenlen = dup_fail_idenlen.append(dup_check_df).reset_index(drop=True)
            # If there are more overlaps than SNPs (SNPs + reference base overlaps), add a redacted entry for each SNP
            else:
                dup_check_snps['Ref_Length'] = "Multiple"
                dup_check_snps['Ref_Start'] = "Multiple"
                dup_check_snps['Ref_End'] = "Multiple"
                dup_check_snps['Ref_Aligned'] = "Multiple"
                dup_check_snps['Query_Length'] = "Multiple"
                dup_check_snps['Query_Start'] = "Multiple"
                dup_check_snps['Query_End'] = "Multiple"
                dup_check_snps['Query_Aligned'] = "Multiple"
                dup_check_snps['Perc_Iden'] = "Multiple"
                dup_check_snps['Cat'] = "Purged_Alignment"
                dup_fail_idenlen = dup_fail_idenlen.append(dup_check_snps)
                
    # SNPs from good alignments
    snp_coords = pd.merge(snp_file, coords_file, how='left')
    snps_pass_idenlen = snp_coords[~snp_coords.isnull().any(1)]
    snps_pass_idenlen = snps_pass_idenlen[snps_pass_idenlen.apply(lambda x: x['Ref_Start'] <= x['Ref_Pos'] <= x['Ref_End'], axis=1)]
    snps_pass_idenlen['Cat'] = "Unchecked"

    # Create empty dataframe to hold SNP data
    processed_snps = snp_coords.iloc[0:0]
    
    # Process duplicate locs from good coords
    loc_tally = snps_pass_idenlen['Ref_Loc'].value_counts()
    solo_locs = loc_tally[loc_tally == 1].index.tolist()
    dup_query = loc_tally[loc_tally > 1]
    dup_locs = dup_query.index.tolist()
    
    # Grab locs covered by 2+ alignments
    snps_pass_dup = snps_pass_idenlen[snps_pass_idenlen.Ref_Loc.isin(solo_locs)]
    snps_check_dup = snps_pass_idenlen[snps_pass_idenlen.Ref_Loc.isin(dup_locs)]
    
    snps_fail_dup = processed_snps.iloc[0:0]
    snps_fail_het = processed_snps.iloc[0:0]

    if len(dup_locs) > 0:
        for loc in dup_locs:
            
            # Find all query overlaps with each reference locus
            dup_check_df = snps_check_dup[snps_check_dup['Ref_Loc'] == loc]
            overlap_count = dup_check_df.shape[0]
            dup_check_snps = snp_file[snp_file['Ref_Loc'] == loc]
            dup_snp_count = dup_check_snps.shape[0]

            # If the same number of SNPs and overlaps appear:
            if overlap_count == dup_snp_count:
                
                # If there is a single base:
                if dup_check_df['Query_Base'].nunique() == 1:

                    # Add the longest/best match to the list
                    sorted_df = dup_check_df.sort_values(by=['Ref_Aligned', 'Perc_Iden'], ascending=[False, False])                    
                    longest_df = sorted_df.head(1)
                    snps_pass_dup = snps_pass_dup.append(longest_df).reset_index(drop=True)
                    
                    # Add worse matches to fail list
                    lower_df = sorted_df.tail(dup_snp_count - 1)
                    lower_df['Cat'] = "Purged_Dup" 
                    snps_fail_dup = snps_fail_dup.append(purged_dup_df).reset_index(drop=True)
                
                # If there are different SNPs:
                else:
                    dup_check_df['Cat'] = "Purged_Het"
                    snps_fail_het = snps_fail_het.append(dup_check_df)
            
            # If there are more overlaps than SNPs (SNPs + reference base overlaps), add a redacted Het entry for each SNP
            else:
                dup_check_snps['Ref_Length'] = "Multiple"
                dup_check_snps['Ref_Start'] = "Multiple"
                dup_check_snps['Ref_End'] = "Multiple"
                dup_check_snps['Ref_Aligned'] = "Multiple"
                dup_check_snps['Query_Length'] = "Multiple"
                dup_check_snps['Query_Start'] = "Multiple"
                dup_check_snps['Query_End'] = "Multiple"
                dup_check_snps['Query_Aligned'] = "Multiple"
                dup_check_snps['Perc_Iden'] = "Multiple"
                dup_check_snps['Cat'] = "Purged_Het"
                snps_fail_het = snps_fail_het.append(dup_check_snps)

    # Process indels        
    snps_pass_indel = snps_pass_dup[snps_pass_dup['Query_Base'] != "."]
    snps_fail_indel = snps_pass_dup[snps_pass_dup['Query_Base'] == "."]

    # Find any SNPs where either sequence has non ACTG.
    valid_bases = ['a', 'A', 'c', 'C', 'g', 'G', 't', 'T','.']

    snps_pass_n = snps_pass_indel[(snps_pass_indel['Ref_Base'].isin(valid_bases)) & (snps_pass_indel['Query_Base'].isin(valid_bases))]        
    snps_fail_n = snps_pass_indel[(~snps_pass_indel['Ref_Base'].isin(valid_bases)) | (~snps_pass_indel['Query_Base'].isin(valid_bases))]
    
    if dup_fail_idenlen.shape[0] > 0:
        dup_fail_idenlen['Cat'] = "Purged_Alignment"

    if snps_fail_n.shape[0] > 0:
        snps_fail_n['Cat'] = "Purged_N"

    if snps_fail_indel.shape[0] > 0:
        snps_fail_indel['Cat'] = "Purged_Indel"

    if snps_fail_dup.shape[0] > 0:
        snps_fail_dup['Cat'] = "Purged_Dup"

    if snps_fail_het.shape[0] > 0:
        snps_fail_het['Cat'] = "Purged_Het"

    # Update processed SNPs
    processed_snps = processed_snps.append(dup_fail_idenlen).append(snps_fail_dup).append(snps_fail_indel).append(snps_fail_n).append(snps_fail_het).reset_index(drop=True)

    if snps_pass_n.shape[0] == 0:
        assert processed_snps.shape[0] == raw_snp_count

    else:

        if len(density_windows) > 0:
            density_locs = density_filter(snps_pass_n['Ref_Loc'],density_windows,max_snps)
            
            snps_fail_density = snps_pass_n[snps_pass_n.Ref_Loc.isin(density_locs)]
            snps_pass_n = snps_pass_n[~snps_pass_n.Ref_Loc.isin(density_locs)]

            if snps_fail_density.shape[0] > 0:
                snps_fail_density['Cat'] = "Purged_Density"
                processed_snps = processed_snps.append(snps_fail_density).reset_index(drop=True)
        
        if snps_pass_n.shape[0] == 0:
            assert processed_snps.shape[0] == raw_snp_count
            
        else:            
            snps_pass_edge = snps_pass_n[(snps_pass_n.Dist_to_Ref_End >= ref_edge) & (snps_pass_n.Dist_to_Query_End >= query_edge)]
            snps_fail_edge = snps_pass_n[(snps_pass_n.Dist_to_Ref_End < ref_edge) | (snps_pass_n.Dist_to_Query_End < query_edge)]
                
            if snps_fail_edge.shape[0] > 0:
                snps_fail_edge['Cat'] = "Filtered_Edge"
                processed_snps = processed_snps.append(snps_fail_edge).reset_index(drop=True)
            
            if snps_pass_edge.shape[0] > 0:
                snps_pass_edge['Cat'] = "SNP"
                processed_snps = processed_snps.append(snps_pass_edge).reset_index(drop=True)

    assert (processed_snps[processed_snps['Cat'] == "Unchecked"].shape[0] == 0) & (processed_snps.shape[0] == raw_snp_count)
    return processed_snps

def density_filter(locs,density_windows,max_snps):

    if len(density_windows) == 0:
        return []
    elif len(density_windows) != len(max_snps):
        sys.exit("density_windows must be the same length as max_snps")
    elif len(locs) == 0:
        return []
    else:
        density_df = pd.DataFrame([item.split('/') for item in locs], columns=['Ref_Contig','Ref_End'])
        density_df['Ref_Loc'] = locs
        
        density_locs = []
        density_bed = makeBED(density_df[['Ref_Contig','Ref_End']])
        
        for i in range(0,len(density_windows)):
            window_bed = density_bed.window(density_bed,c=True, w=density_windows[i]).to_dataframe()
            window_df = window_bed[window_bed['name'] > max_snps[i]]
            if window_df.shape[0] > 0:
                density_locs = density_locs + ["/".join([str(x[0]),str(x[1])]) for x in list(zip(window_df.chrom, window_df.end))]
                density_bed = makeBED(density_df[~density_df.Ref_Loc.isin(density_locs)][['Ref_Contig','Ref_End']])
    
    return density_locs

def fasta_info(file_path):
    records = list(SeqIO.parse(file_path, 'fasta'))
    contig_count = int(len(records))
    assembly_bases = int(sum(len(record) for record in records))
    with open(file_path, 'rb') as file:
        sha256 = hashlib.sha256(file.read()).hexdigest()
    return [file_path, contig_count, assembly_bases, sha256]

#### 01: Read in arguments ####

# Get isolate data
query = str(sys.argv[1])
query_fasta = str(sys.argv[2])
query_data = [query] + fasta_info(query_fasta)
query_string = [x+":"+str(y) for x,y in zip(['Query_ID','Query_Assembly','Query_Contig_Count','Query_Assembly_Bases','Query_SHA256'],query_data)]

reference = str(sys.argv[3])
reference_fasta = str(sys.argv[4])
reference_data = [reference] + fasta_info(reference_fasta)
reference_string = [x+":"+str(y) for x,y in zip(['Reference_ID','Reference_Assembly','Reference_Contig_Count','Reference_Assembly_Bases','Reference_SHA256'],reference_data)]

report_id = query + "__vs__" + reference

# Get directories
output_dir = os.path.normpath(os.path.abspath(sys.argv[5]))

mummer_dir = output_dir + "/MUMmer_Output"
mum_snps_dir = mummer_dir + "/snps"
mum_report_dir = mummer_dir + "/report"
mum_coords_dir = mummer_dir + "/1coords"

log_dir = output_dir + "/logs"

snpdiffs_dir = output_dir + "/snpdiffs"
snpdiffs_file = snpdiffs_dir + "/" + report_id + ".snpdiffs"

# Set filtering criteria
min_cov = float(sys.argv[6])
min_iden = float(sys.argv[7])
min_len = int(sys.argv[8])

density_windows = [int(x) for x in sys.argv[9].split(",")]
max_snps = [int(x) for x in sys.argv[10].split(",")]

qc_density =  ",".join([str(x) for x in density_windows])
qc_maxsnps =  ",".join([str(x) for x in max_snps])

assert len(density_windows) == len(max_snps)

ref_edge = int(sys.argv[11])
query_edge = int(sys.argv[12])

run_mode = str(sys.argv[13])

# Create QC String
qc_string = "_".join([str(x) for x in [min_cov,min_iden,min_len,ref_edge,query_edge,qc_density,qc_maxsnps]]) 
merged_ref_bed = BedTool([])

#### 02: Read in MUmmer report data ####
[ref_bases,percent_ref_aligned,query_bases,percent_query_aligned] = parseMUmmerReport(mum_report_dir,report_id)

if (percent_ref_aligned < min_cov) | (percent_query_aligned < min_cov):
    sample_category = "Purged_Min_Coverage"
    percent_ref_aligned_filtered = "NA"
    percent_query_aligned_filtered = "NA"
    median_percent_identity = "NA"
    final_snp_count = "NA"
    median_snp_perc_iden = "NA"
    reject_snps_alignment_count = "NA"
    reject_snps_n_count = "NA"
    reject_snps_indel_count = "NA"
    reject_snps_dup_count = "NA"
    reject_snps_het_count = "NA"
    reject_snps_density_count = "NA"
    reject_snps_edge_count = "NA"

else:
    #### 03: Process MUmmer coords file ####
    [coords_file,bad_coords_file,merged_ref_bed] = parseMUmmerCoords(mum_coords_dir,report_id,min_iden,min_len)

    # STOP if the coordinates file is empty after filtering based on <perc_iden>
    if coords_file.shape[0] == 0:
        sample_category = "Purged_Filter_Coverage"
        percent_ref_aligned_filtered = "NA"
        percent_query_aligned_filtered = "NA"
        median_percent_identity = "NA"
        final_snp_count = "NA"
        median_snp_perc_iden = "NA"
        reject_snps_alignment_count = "NA"
        reject_snps_n_count = "NA"
        reject_snps_indel_count = "NA"
        reject_snps_dup_count = "NA"
        reject_snps_het_count = "NA"
        reject_snps_density_count = "NA"
        reject_snps_edge_count = "NA"
    
    else:

        # Get information for filtered mappings
        median_percent_identity = f"{ coords_file.Perc_Iden.median():.2f}"
        filtered_ref_covered = sum(int(interval.length) for interval in merged_ref_bed)
        percent_ref_aligned_filtered = 100*(filtered_ref_covered/ref_bases)
        percent_query_aligned_filtered = 100*(filtered_ref_covered/query_bases)
        
        # STOP if the reference or query is not covered by at least <min_cov>
        if (percent_ref_aligned_filtered < min_cov) & (percent_query_aligned_filtered < min_cov):

            percent_ref_aligned_filtered = f"{percent_ref_aligned_filtered:.2f}"
            percent_query_aligned_filtered = f"{percent_query_aligned_filtered:.2f}"
            sample_category = "Purged_Filter_Coverage"
            final_snp_count = "NA"
            median_snp_perc_iden = "NA"
            reject_snps_alignment_count = "NA"
            reject_snps_n_count = "NA"
            reject_snps_indel_count = "NA"
            reject_snps_dup_count = "NA"
            reject_snps_het_count = "NA"
            reject_snps_density_count = "NA"
            reject_snps_edge_count = "NA"
        
        else:

            percent_ref_aligned_filtered = f"{percent_ref_aligned_filtered:.2f}"
            percent_query_aligned_filtered = f"{percent_query_aligned_filtered:.2f}"            

            #### 04: Process MUmmer SNPs file ####

            # Sample passed coverage QC
            sample_category = "PASS"

            # Read in SNP file
            snp_file = parseMUmmerSNPs(mum_snps_dir,report_id)
            
            # STOP if no SNPs detected
            if snp_file.shape[0] == 0:
                median_snp_perc_iden = "NA"
                final_snp_count = 0
                median_snp_perc_iden = "NA"
                reject_snps_alignment_count = "NA"
                reject_snps_n_count = "NA"
                reject_snps_indel_count = "NA"
                reject_snps_dup_count = "NA"
                reject_snps_het_count = "NA"
                reject_snps_density_count = "NA"
                reject_snps_edge_count = "NA"
            
            else:

                # Characterize and filter SNPs based on user criteria (perc_iden,query_edge/ref_edge,density)
                filtered_snps = filterSNPs(snp_file,coords_file,bad_coords_file,density_windows,max_snps,ref_edge,query_edge)
                
                final_snp_df = filtered_snps[filtered_snps.Cat == "SNP"]
                final_snp_count = final_snp_df.shape[0]

                if final_snp_count == 0:
                    median_snp_perc_iden = "NA"
                else:
                    median_snp_perc_iden = f"{final_snp_df.Perc_Iden.median():.2f}"
                
                reject_snps_alignment_count = filtered_snps[filtered_snps.Cat =="Purged_Alignment"].shape[0]
                reject_snps_n_count = filtered_snps[filtered_snps.Cat =="Purged_N"].shape[0]
                reject_snps_indel_count = filtered_snps[filtered_snps.Cat =="Purged_Indel"].shape[0]
                reject_snps_dup_count = filtered_snps[filtered_snps.Cat =="Purged_Dup"].shape[0]
                reject_snps_het_count = filtered_snps[filtered_snps.Cat =="Purged_Het"].shape[0]
                reject_snps_density_count = filtered_snps[filtered_snps.Cat =="Purged_Density"].shape[0]
                reject_snps_edge_count = filtered_snps[filtered_snps.Cat =="Filtered_Edge"].shape[0]

# Create header
snpdiffs_header=[]
snpdiffs_header.append("#\t" + "\t".join(query_string) + "\t" + "\t".join(reference_string) + "\t" + "\t".join([
"Category:"+str(sample_category),"SNPs:"+str(final_snp_count),"Median_Percent_Identity:"+str(median_percent_identity),"Median_SNP_Percent_Identity:"+str(median_snp_perc_iden),
"Percent_Query_Aligned:"+str(percent_query_aligned_filtered),"Percent_Reference_Aligned:"+str(percent_ref_aligned_filtered),"Purged_Alignment:"+str(reject_snps_alignment_count),
"Purged_N:"+str(reject_snps_n_count),"Purged_Indel:"+str(reject_snps_indel_count),"Purged_Duplicate:"+str(reject_snps_dup_count),"Purged_Het:"+str(reject_snps_het_count),
"Purged_Density:"+str(reject_snps_density_count),"Filtered_Edge:"+str(reject_snps_edge_count),"QC_String:"+qc_string]))

# Save file
with open(snpdiffs_file,"w") as file:
    file.write("\n".join(snpdiffs_header) + "\n")

if len(merged_ref_bed) > 0:
    bed_df = merged_ref_bed.to_dataframe()
    with open(snpdiffs_file,"a") as file:
        for index, row in bed_df.iterrows():
            file.write("##\t" + "\t".join([query,reference]) + "\t" + "\t".join(map(str, row))+"\n")

if str(reject_snps_alignment_count) != "NA":
        het_snps = filtered_snps[filtered_snps['Ref_Aligned'] == "Multiple"]
        non_het_snps = filtered_snps[filtered_snps['Ref_Aligned'] != "Multiple"]
        non_het_snps['Ref_Aligned'] = non_het_snps['Ref_Aligned'].apply(lambda x: f'{x:.0f}')
        filtered_snps = non_het_snps.append(het_snps)
        filtered_snps['Query'] = query
        filtered_snps['Reference'] = reference
        filtered_snps[["Query","Reference","Ref_Loc","Cat","Ref_Base","Query_Base","Query_Loc","Ref_Aligned","Perc_Iden"]].to_csv(snpdiffs_file, sep="\t",mode='a', header=False, index=False)

if run_mode == "align":
    print(",".join([str(x) for x in query_data]))
    print(",".join([str(x) for x in reference_data]))
else:
    print(",".join([query_data,reference,snpdiffs_file]))
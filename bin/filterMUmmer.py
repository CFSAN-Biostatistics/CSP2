#!/usr/bin/env python3

import sys
import os
import pandas as pd
import re
from pybedtools import BedTool
import warnings
warnings.filterwarnings("ignore")

# Return a BEDTool object from basic dataframe of (1) ref contig, (2) 1-based ref pos, (3) metadata
def makeBED(bed_df):
    if bed_df.shape[1] == 2:
        bed_df.columns = ['Ref_Contig','Ref_End']
        bed_df['Ref_Start'] = bed_df['Ref_End'] - 1
        bed_file = BedTool.from_dataframe(bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort()
    elif bed_df.shape[1] == 3:
        bed_df.columns = ['Ref_Contig','Ref_End','Metadata']
        bed_df['Ref_Start'] = bed_df['Ref_End'] - 1
        bed_file = BedTool.from_dataframe(bed_df[['Ref_Contig','Ref_Start','Ref_End','Metadata']]).sort()
    else:
        sys.exit("bed_df should have 2 or 3 columns")
    return bed_file

def parseMUmmerReport(mummer_dir,report_id):
    report_data = []
    with open(mummer_dir+"/"+report_id+".report", 'r') as report_file:

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

    return report_data

def parseMUmmerCoords(mummer_dir,report_id,perc_iden,min_len):
    
    coords_file = pd.read_csv(mummer_dir+"/"+report_id+".1coords",sep="\t",index_col=False,
    names=['Ref_Start','Ref_End','Query_Start','Query_End',
    'Ref_Aligned','Query_Aligned','Perc_Iden',
    'Ref_Length','Query_Length','Ref_Cov',
    'Query_Cov','Ref_Contig','Query_Contig'])

    coords_file = coords_file[coords_file.Ref_Aligned >= min_len]
    coords_file = coords_file[coords_file.Perc_Iden >= perc_iden]
    return coords_file[['Ref_Contig','Ref_Length','Ref_Start','Ref_End','Ref_Aligned','Query_Contig','Query_Length','Query_Start','Query_End','Query_Aligned']]
    
def parseMUmmerSNPs(mummer_dir,report_id):
    
    snp_file = pd.read_csv(mummer_dir+"/"+report_id+".snps",sep="\t",index_col=False,
        names=['Ref_Pos','Ref_Base','Query_Base','Query_Pos',
        'SNP_Buffer','Dist_to_End','Ref_Length','Query_Length',
        'Ref_Direction','Query_Direction','Ref_Contig','Query_Contig'])

    # Print indels if there are any
    indel_file = snp_file[(snp_file.Ref_Base == ".") | (snp_file.Query_Base == ".")]

    if indel_file.shape[0] > 0:
        indel_file['Ref_Loc'] = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(indel_file.Ref_Contig, indel_file.Ref_Pos))]
        indel_file['Query_Loc'] = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(indel_file.Query_Contig, indel_file.Query_Pos))]
        indel_file.to_csv(mummer_dir+"/"+report_id+"_Indels.tsv",sep="\t",index=False)

    # Remove indels
    snp_file = snp_file[(snp_file.Ref_Base != ".") & (snp_file.Query_Base != ".")]
    
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

        return snp_file[['Ref_Contig','Ref_Pos','Ref_Loc','Query_Contig','Query_Pos','Query_Loc','Dist_to_Ref_End','Dist_to_Query_End','Ref_Base','Query_Base','Ref_Direction','Query_Direction']]

def filterSNPs(snp_coords,ref_edge,query_edge):

    # Identify SNPs on contigs that fail perc_iden filter
    rejected_snps_iden = snp_coords[snp_coords.isnull().any(1)]
    rejected_snps_iden_count = rejected_snps_iden.shape[0]
    if rejected_snps_iden_count > 0:
        rejected_snps_iden['Cat'] = "Purged_Identity_Length"

    snps_pf_iden = snp_coords[~snp_coords.isnull().any(1)]
    
    # Identify duplicated SNPs (Reference positions covered by multiple mappings) and select the SNP deriving from the longest alignment
    snps_pf_iden_dup = snps_pf_iden.drop_duplicates(subset='Ref_Loc', keep=False)
    duplicated_mask = snps_pf_iden.duplicated(subset='Ref_Loc', keep=False)
    rejected_snps_dup = snps_pf_iden[duplicated_mask]

    if len(duplicated_mask) > 0:
        longest_dup_df = rejected_snps_dup.loc[rejected_snps_dup.groupby('Ref_Loc')['Ref_Aligned'].idxmax()]
        rejected_snps_dup = snps_pf_iden[~snps_pf_iden.index.isin(longest_dup_df.index) & ~snps_pf_iden.index.isin(snps_pf_iden_dup.index)]
        rejected_snps_dup['Cat'] = "Purged_Dup"
        snps_pf_iden_dup = pd.concat([snps_pf_iden_dup, longest_dup_df])

    # Look for high-density regions at the 1000bp (6 or fewer), 125bp (4 or fewer), 15bp (2 or fewer)

    # Create BED file for preserved SNPs
    preserved_bed = makeBED(snps_pf_iden_dup[['Ref_Contig','Ref_Pos']])
    density_locs = []
    
    w_1000 = preserved_bed.window(preserved_bed,c=True, w=1000)
    w_1000_df = pd.read_table(w_1000.fn, names=['Ref_Contig', 'Ref_Pos', 'Ref_End', 'Count']).query("`Count` > 6")
    w_1000_locs = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(w_1000_df.Ref_Contig, w_1000_df.Ref_End))]
    rejected_density_1000 = snps_pf_iden_dup[snps_pf_iden_dup.Ref_Loc.isin(w_1000_locs)]

    if w_1000_df.shape[0] > 0:
        density_locs = density_locs + w_1000_locs
        preserved_bed = makeBED(snps_pf_iden_dup[~snps_pf_iden_dup.Ref_Loc.isin(density_locs)][['Ref_Contig','Ref_Pos']])
        rejected_density_1000['Cat'] = "Purged_Density_1000"
    
    w_125 = preserved_bed.window(preserved_bed,c=True, w=125)
    w_125_df = pd.read_table(w_125.fn, names=['Ref_Contig', 'Ref_Pos', 'Ref_End', 'Count']).query("`Count` > 4")
    w_125_locs = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(w_125_df.Ref_Contig, w_125_df.Ref_End))]
    rejected_density_125 = snps_pf_iden_dup[snps_pf_iden_dup.Ref_Loc.isin(w_125_locs)]

    if w_125_df.shape[0] > 0:
        density_locs = density_locs + w_125_locs
        preserved_bed = makeBED(snps_pf_iden_dup[~snps_pf_iden_dup.Ref_Loc.isin(density_locs)][['Ref_Contig','Ref_Pos']])
        rejected_density_125['Cat'] = "Purged_Density_125"

    w_15 = preserved_bed.window(preserved_bed,c=True, w=15)
    w_15_df = pd.read_table(w_15.fn, names=['Ref_Contig', 'Ref_Pos', 'Ref_End', 'Count']).query("`Count` > 2")
    w_15_locs = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(w_15_df.Ref_Contig, w_15_df.Ref_End))]
    rejected_density_15 = snps_pf_iden_dup[snps_pf_iden_dup.Ref_Loc.isin(w_15_locs)]

    if w_15_df.shape[0] > 0:
        density_locs = density_locs + w_15_locs
        preserved_bed = makeBED(snps_pf_iden_dup[~snps_pf_iden_dup.Ref_Loc.isin(density_locs)][['Ref_Contig','Ref_Pos']])
        rejected_density_15['Cat'] = "Purged_Density_15"
        
    snps_pf_iden_dup_density = snps_pf_iden_dup[~snps_pf_iden_dup.Ref_Loc.isin(density_locs)]

    # Identify SNPs that are too close to the contig edges
    rejected_snps_edge = snps_pf_iden_dup_density[(snps_pf_iden_dup_density.Dist_to_Ref_End < ref_edge) | (snps_pf_iden_dup_density.Dist_to_Query_End < query_edge)]
    rejected_snps_edge_count = rejected_snps_edge.shape[0]
    if rejected_snps_edge_count > 0:
        rejected_snps_edge['Cat'] = "Filtered_Edge"

    final_snp_df = snps_pf_iden_dup_density[(snps_pf_iden_dup_density.Dist_to_Ref_End >= ref_edge) & (snps_pf_iden_dup_density.Dist_to_Query_End >= query_edge)]
    final_snp_df['Cat'] = "Yenta_SNP"

    all_snps = final_snp_df.append(rejected_snps_iden).append(rejected_snps_edge).append(rejected_snps_dup).append(rejected_density_1000).append(rejected_density_125).append(rejected_density_15)
    
    return all_snps

#### 01: Read in arguments ####

# Get isolate data
query = str(sys.argv[1])
reference = str(sys.argv[2])
report_id = str(sys.argv[3])

# Get raw MUmmer dir
mummer_dir = os.path.normpath(sys.argv[4])
mummer_parent_dir = os.path.dirname(mummer_dir)

# Set filtering criteria
align_cov = float(sys.argv[5])
perc_iden = float(sys.argv[6])
ref_edge = int(sys.argv[7])
query_edge = int(sys.argv[8])
min_len = int(sys.argv[9])

# Create dummy columns
final_columns = ['Ref_Contig','Ref_Pos','Ref_Loc','Query_Contig','Query_Pos','Query_Loc','Dist_to_Ref_End','Dist_to_Query_End','Ref_Base','Query_Base','Ref_Direction','Query_Direction','Ref_Length','Ref_Start','Ref_End','Ref_Aligned','Query_Length','Query_Start','Query_End','Query_Aligned','Cat','Ref','Query\n']

#### 02: Read in MUmmer report data ####
report_data = parseMUmmerReport(mummer_dir,report_id)

# Fetch assembly data
ref_seqs = report_data[report_data['Measure'] == 'TotalSeqs']['Ref_Value'].iloc[0]
query_seqs = report_data[report_data['Measure'] == 'TotalSeqs']['Query_Value'].iloc[0]

ref_bases = report_data[report_data['Measure'] == 'TotalBases']['Ref_Value'].iloc[0]
query_bases = report_data[report_data['Measure'] == 'TotalBases']['Query_Value'].iloc[0]

# Fetch alignment data
ref_aligned = report_data[report_data['Measure'] == 'AlignedBases']['Ref_Value'].iloc[0]
query_aligned = report_data[report_data['Measure'] == 'AlignedBases']['Query_Value'].iloc[0]

# Get gSNPs + gIndels
gsnps = report_data[report_data['Measure'] == 'TotalGSNPs']['Query_Value'].iloc[0]
gindels = report_data[report_data['Measure'] == 'TotalGIndels']['Query_Value'].iloc[0]

# STOP if the reference or query is not covered by at least <align_cov>
percent_ref_aligned = 100*(float(ref_aligned)/ref_bases)
percent_query_aligned = 100*(float(query_aligned)/query_bases)

if (percent_ref_aligned < align_cov) & (percent_query_aligned < align_cov):
    sample_category = "Purged_Filter_Coverage"
    percent_ref_aligned_filtered = "NA"
    percent_query_aligned_filtered = "NA"
    final_snp_count = "NA"
    rejected_snps_iden_count = "NA"
    rejected_snps_edge_count = "NA"
    rejected_snps_dup_count = "NA"
    rejected_snps_density1000_count = "NA"
    rejected_snps_density125_count = "NA"
    rejected_snps_density15_count = "NA"
    
    # Save empty TSV
    with open(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.tsv","w+") as file:
        file.write("\t".join(final_columns))

else:

    #### 03: Process MUmmer coords file ####
    coords_file = parseMUmmerCoords(mummer_dir,report_id,perc_iden,min_len)

    # STOP if the coordinates file is empty after filtering based on <perc_iden>
    if coords_file.shape[0] == 0:
        sample_category = "Purged_Filter_Coverage"
        percent_ref_aligned_filtered = "NA"
        percent_query_aligned_filtered = "NA"
        final_snp_count = "NA"
        rejected_snps_iden_count = "NA"
        rejected_snps_edge_count = "NA"
        rejected_snps_dup_count = "NA"
        rejected_snps_density1000_count = "NA"
        rejected_snps_density125_count = "NA"
        rejected_snps_density15_count = "NA"

        # Save empty TSV
        with open(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.tsv","w+") as file:
            file.write("\t".join(final_columns))
    else:
        # Get information for filtered mappings
        filtered_ref_bases = sum(coords_file.Ref_Aligned)
        percent_ref_aligned_filtered = 100*(float(filtered_ref_bases)/ref_bases)

        filtered_query_bases = sum(coords_file.Query_Aligned)
        percent_query_aligned_filtered = 100*(float(filtered_query_bases)/query_bases)

        # STOP if the reference or query is not covered by at least <align_cov>
        if (percent_ref_aligned_filtered < align_cov) & (percent_query_aligned_filtered < align_cov):
            sample_category = "Purged_Filter_Coverage"
            final_snp_count = "NA"
            rejected_snps_iden_count = "NA"
            rejected_snps_edge_count = "NA"
            rejected_snps_dup_count = "NA"
            rejected_snps_density1000_count = "NA"
            rejected_snps_density125_count = "NA"
            rejected_snps_density15_count = "NA"

            # Save empty TSV
            with open(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.tsv","w+") as file:
                file.write("\t".join(final_columns))

        else:

            #### 04: Process MUmmer SNPs file ####

            # Sample passed coverage QC
            sample_category = "PASS"

            # Read in SNP file
            snp_file = parseMUmmerSNPs(mummer_dir,report_id)

            # STOP if no SNPs detected
            if snp_file.shape[0] == 0:
                final_snp_count = 0
                rejected_snps_iden_count = 0
                rejected_snps_edge_count = 0
                rejected_snps_dup_count = 0
                rejected_snps_density1000_count = 0
                rejected_snps_density125_count = 0
                rejected_snps_density15_count = 0

                # Save empty TSV
                with open(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.tsv","w+") as file:
                    file.write("\t".join(final_columns))
            
            else:
                # Merge SNP data and coordinate data
                snp_coords = pd.merge(snp_file, coords_file, how='left')
                snp_coords['Cat'] = "Undefined"
    
                # Characterize and filter SNPs based on user criteria (perc_iden,query_edge/ref_edge,density)
                filtered_snps = filterSNPs(snp_coords,ref_edge,query_edge)
                filtered_snps['Ref'] = reference
                filtered_snps['Query'] = query
                final_snp_count = filtered_snps[filtered_snps.Cat == "Yenta_SNP"].shape[0]
                rejected_snps_iden_count = filtered_snps[filtered_snps.Cat =="Purged_Identity_Length"].shape[0]
                rejected_snps_edge_count = filtered_snps[filtered_snps.Cat =="Filtered_Edge"].shape[0]
                rejected_snps_dup_count = filtered_snps[filtered_snps.Cat =="Purged_Dup"].shape[0]
                rejected_snps_density1000_count = filtered_snps[filtered_snps.Cat =="Purged_Density_1000"].shape[0]
                rejected_snps_density125_count = filtered_snps[filtered_snps.Cat =="Purged_Density_125"].shape[0]
                rejected_snps_density15_count = filtered_snps[filtered_snps.Cat =="Purged_Density_15"].shape[0]

                # Save SNP data
                if filtered_snps.shape[0] > 0:
                    filtered_snps.to_csv(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.tsv",sep="\t",index=False)
                    final_bed = makeBED(filtered_snps[['Ref_Contig','Ref_Pos','Cat']])
                    pd.read_table(final_bed.fn).to_csv(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.bed",sep="\t",index=False)
                else:
                    # Save empty TSV
                    with open(mummer_parent_dir+"/"+report_id+"_MUmmer_SNPs.tsv","w+") as file:
                        file.write("\t".join(final_columns))

# Print output for Nextflow
print(",".join([
    str(query),
    str(reference),
    str(query_seqs),
    str(ref_seqs),
    str(query_bases),
    str(ref_bases),
    str(percent_query_aligned_filtered),
    str(percent_ref_aligned_filtered),
    str(sample_category),
    str(final_snp_count),
    str(gsnps),
    str(rejected_snps_iden_count),
    str(rejected_snps_edge_count),
    str(rejected_snps_dup_count),
    str(rejected_snps_density1000_count),
    str(rejected_snps_density125_count),
    str(rejected_snps_density15_count)]))
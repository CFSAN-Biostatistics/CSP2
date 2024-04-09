#!/usr/bin/env python3

import sys
import os
import pandas as pd
import re
from pybedtools import BedTool,helpers,Interval
import warnings
import numpy as np
import hashlib
from Bio import SeqIO
import subprocess

warnings.filterwarnings("ignore")

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
    report_data = report_data[report_data['Measure'].isin(['TotalSeqs','TotalBases','AlignedBases',
                                                           'TotalSNPs','TotalGSNPs','TotalIndels','TotalGIndels',
                                                           'Breakpoints','Relocations','Translocations','Inversions',
                                                           'Insertions','TandemIns'])]
    
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
    
    # Fetch SNP data
    g_snps = report_data[report_data['Measure'] == 'TotalGSNPs']['Ref_Value'].iloc[0]
    g_indels = report_data[report_data['Measure'] == 'TotalGIndels']['Ref_Value'].iloc[0]    
    
    ref_breakpoints = report_data[report_data['Measure'] == 'Breakpoints']['Ref_Value'].iloc[0]
    query_breakpoints = report_data[report_data['Measure'] == 'Breakpoints']['Query_Value'].iloc[0]
    
    ref_relocations = report_data[report_data['Measure'] == 'Relocations']['Ref_Value'].iloc[0]
    query_relocations = report_data[report_data['Measure'] == 'Relocations']['Query_Value'].iloc[0]
    
    ref_translocations = report_data[report_data['Measure'] == 'Translocations']['Ref_Value'].iloc[0]
    query_translocations = report_data[report_data['Measure'] == 'Translocations']['Query_Value'].iloc[0]

    ref_inversions = report_data[report_data['Measure'] == 'Inversions']['Ref_Value'].iloc[0]
    query_inversions = report_data[report_data['Measure'] == 'Inversions']['Query_Value'].iloc[0]                                                                                                 

    ref_insertions = report_data[report_data['Measure'] == 'Insertions']['Ref_Value'].iloc[0]
    query_insertions = report_data[report_data['Measure'] == 'Insertions']['Query_Value'].iloc[0]                                                                                                 
    
    ref_tandem = report_data[report_data['Measure'] == 'TandemIns']['Ref_Value'].iloc[0]
    query_tandem = report_data[report_data['Measure'] == 'TandemIns']['Query_Value'].iloc[0]
    
    return [ref_bases,percent_ref_aligned,
            query_bases,percent_query_aligned,
            g_snps,g_indels,
            ref_breakpoints,query_breakpoints,
            ref_relocations,query_relocations,
            ref_translocations,query_translocations,
            ref_inversions,query_inversions,
            ref_insertions,query_insertions,
            ref_tandem,query_tandem]

def parseMUmmerCoords(mum_coords_dir,report_id,reference_chr_bed,query_chr_bed):
    
    coords_file = pd.read_csv(mum_coords_dir+"/"+report_id+".1coords",sep="\t",index_col=False,
    names=['Ref_Start','Ref_End','Query_Start','Query_End',
    'Ref_Aligned','Query_Aligned','Perc_Iden',
    'Ref_Length','Query_Length','Ref_Cov',
    'Query_Cov','Ref_Contig','Query_Contig'],   
    dtype={
        'Ref_Start': int,
        'Ref_End': int,
        'Query_Start': int,
        'Query_End': int,
        'Ref_Aligned': int,
        'Query_Aligned': int,
        'Perc_Iden': float,
        'Ref_Length': int,
        'Query_Length': int,
        'Ref_Cov': float,
        'Query_Cov': float,
        'Ref_Contig': str,
        'Query_Contig': str
    })
    
    if coords_file.shape[0] > 0:
        
        coords_file['Ref_Start'] = (coords_file['Ref_Start'] - 1)
        coords_file = coords_file.apply(adjust_query_row, axis=1)
                
        merged_ref_bed = BedTool.from_dataframe(coords_file[['Ref_Contig','Ref_Start','Ref_End']]).sort().merge()
        merged_query_bed = BedTool.from_dataframe(coords_file[['Query_Contig','Query_Start','Query_End']]).sort().merge()
        
        covered_ref_regions = reference_chr_bed.intersect(merged_ref_bed)
        uncovered_ref_regions = reference_chr_bed.subtract(merged_ref_bed)
                
        covered_query_regions = query_chr_bed.intersect(merged_query_bed)
        uncovered_query_regions = query_chr_bed.subtract(merged_query_bed)
    else:
        merged_ref_bed = BedTool([])
        merged_query_bed = BedTool([])
        
        covered_ref_regions = BedTool([])
        uncovered_ref_regions = BedTool([])  
        
        covered_query_regions = BedTool([])
        uncovered_query_regions = BedTool([])  
    
    ref_genome_size = calculate_total_length(reference_chr_bed)
    ref_covered_size = calculate_total_length(covered_ref_regions)
    ref_uncovered_size = calculate_total_length(uncovered_ref_regions)

    query_genome_size = calculate_total_length(query_chr_bed)
    query_covered_size = calculate_total_length(covered_query_regions)
    query_uncovered_size = calculate_total_length(uncovered_query_regions)
    
    assert ref_covered_size + ref_uncovered_size == ref_genome_size
    assert query_covered_size + query_uncovered_size == query_genome_size
    
    # Create return_df
    return_columns = ['Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned','Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned','Perc_Iden']

    return_df = pd.DataFrame(columns=return_columns)
    
    if coords_file.shape[0] > 0:
        return_df = pd.concat([return_df,coords_file[return_columns]],ignore_index=True)

    if uncovered_ref_regions.count() > 0:
        return_uncovered_df = uncovered_ref_regions.to_dataframe().sort_values(['chrom','start']).rename(columns={'chrom':'Ref_Contig','start':'Ref_Start','end':'Ref_End'}).assign(Ref_Aligned = ".",Query_Aligned = ".",Query_Contig=".",Query_Start=".",Query_End=".",Ref_Length = ".",Query_Length=".",Perc_Iden = ".")
        return_df = pd.concat([return_df,return_uncovered_df[return_columns]],ignore_index=True)
    
    if uncovered_query_regions.count() > 0:
        return_uncovered_df = uncovered_query_regions.to_dataframe().sort_values(['chrom','start']).rename(columns={'chrom':'Query_Contig','start':'Query_Start','end':'Query_End'}).assign(Ref_Aligned = ".",Query_Aligned = ".",Ref_Contig=".",Ref_Start=".",Ref_End=".",Ref_Length = ".",Query_Length=".",Perc_Iden = ".")
        return_df = pd.concat([return_df,return_uncovered_df[return_columns]],ignore_index=True) 
    
    return return_df

def adjust_query_row(row):
    if set(['Query_Start', 'Query_End', 'Query_Contig']).issubset(row.index):
        if row['Query_Start'] > row['Query_End']:
            row['Query_Start'], row['Query_End'] = row['Query_End'] - 1, row['Query_Start']
        else:
            row['Query_Start'] = row['Query_Start'] - 1
    return row

def parseMUmmerSNPs(mum_snps_dir,report_id,coords_file):
    
    return_columns = ['Ref_Contig','Start_Ref','Ref_Pos',
                'Query_Contig','Start_Query','Query_Pos',
                'Ref_Loc','Query_Loc',
                'Ref_Start','Ref_End',
                'Query_Start','Query_End',
                'Ref_Base','Query_Base',
                'Dist_to_Ref_End','Dist_to_Query_End',
                'Ref_Aligned','Query_Aligned',
                'Query_Direction','Perc_Iden','Cat']
    
    coords_columns = ['Ref_Contig','Ref_Start','Ref_End',
                      'Query_Contig','Query_Start','Query_End',
                      'Ref_Aligned','Query_Aligned','Perc_Iden']

    snp_file = pd.read_csv(mum_snps_dir+"/"+report_id+".snps",sep="\t",index_col=False,
        names=['Ref_Pos','Ref_Base','Query_Base','Query_Pos',
        'SNP_Buffer','Dist_to_End','Ref_Length','Query_Length',
        'Ref_Direction','Query_Direction','Ref_Contig','Query_Contig'])
    
    if snp_file.shape[0] == 0:
        return pd.DataFrame(columns=return_columns)
    else:
        total_snp_count = snp_file.shape[0]

        snp_file['Start_Ref'] = snp_file['Ref_Pos'] - 1
        snp_file['Start_Query'] = snp_file['Query_Pos'] - 1
        
        snp_file['Dist_to_Ref_End'] = [min([x,y]) for x,y in zip(snp_file['Ref_Pos'],snp_file['Ref_Length'] - snp_file['Ref_Pos'])]
        snp_file['Dist_to_Query_End'] = [min([x,y]) for x,y in zip(snp_file['Query_Pos'],snp_file['Query_Length'] - snp_file['Query_Pos'])]

        snp_file['Ref_Loc'] = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(snp_file.Ref_Contig, snp_file.Ref_Pos))]
        snp_file['Query_Loc'] = ["/".join([str(x[0]),str(x[1])]) for x in list(zip(snp_file.Query_Contig, snp_file.Query_Pos))]    

        valid_bases = ['a', 'A', 'c', 'C', 'g', 'G', 't', 'T',"."]
        
        invalid_file = snp_file[(~snp_file['Ref_Base'].isin(valid_bases)) | (~snp_file['Query_Base'].isin(valid_bases))]
        snp_file = snp_file[(snp_file['Ref_Base'].isin(valid_bases)) & (snp_file['Query_Base'].isin(valid_bases))]
        indel_file = snp_file[(snp_file['Query_Base'] == ".") | (snp_file['Ref_Base'] == ".")]
        snp_file = snp_file[~((snp_file['Query_Base'] == ".") | (snp_file['Ref_Base'] == "."))]
        
        # Gather coordinate information       
        if snp_file.shape[0] == 0:
            snp_coords = pd.DataFrame(columns=return_columns)
        else:
            snp_coords = makeSNPCoords(snp_file,coords_file,snp_file.shape[0])
            snp_coords['Cat'] = "SNP"
        
        if indel_file.shape[0] == 0:
            indel_coords = pd.DataFrame(columns=return_columns)
        else:
            indel_coords = makeSNPCoords(indel_file,coords_file,indel_file.shape[0])
            indel_coords['Cat'] = "Indel"
        
        if invalid_file.shape[0] == 0:
            invalid_coords = pd.DataFrame(columns=return_columns)
        else:
            invalid_coords = makeSNPCoords(invalid_file,coords_file,invalid_file.shape[0])
            invalid_coords['Cat'] = "Invalid"
                  
        return_df = pd.concat([snp_coords,indel_coords,invalid_coords],ignore_index=True)[return_columns]
        
        assert return_df.shape[0] == total_snp_count
        return return_df
    
def makeSNPCoords(snp_df, coords_df, row_count):
    snp_coords = pd.merge(snp_df, coords_df, how='left').dropna()
    
    condition = ((snp_coords['Ref_Start'] <= snp_coords['Ref_Pos']) & 
                 (snp_coords['Ref_Pos'] <= snp_coords['Ref_End']) & 
                 (snp_coords['Query_Start'] <= snp_coords['Query_Pos']) & 
                 (snp_coords['Query_Pos'] <= snp_coords['Query_End']))
    snp_coords = snp_coords[condition]
    
    if row_count != snp_coords.shape[0]:
        snp_coords.sort_values(by=['Ref_Aligned', 'Perc_Iden'], ascending=False, inplace=True)
        snp_coords.drop_duplicates(subset=['Ref_Loc','Query_Loc'], keep='first', inplace=True)
            
    return snp_coords
    
def fasta_info(file_path):
    records = list(SeqIO.parse(file_path, 'fasta'))
    contig_count = int(len(records))
    lengths = sorted([len(record) for record in records], reverse=True)
    assembly_bases = sum(lengths)

    with open(file_path, 'rb') as file:
        sha256 = hashlib.sha256(file.read()).hexdigest()

    cumulative_length = 0
    n50 = None
    n90 = None
    l50 = None
    l90 = None
    
    for i, length in enumerate(lengths, start=1):
        cumulative_length += length
        if cumulative_length >= assembly_bases * 0.5 and n50 is None:
            n50 = length
            l50 = i
        if cumulative_length >= assembly_bases * 0.9 and n90 is None:
            n90 = length
            l90 = i
        if n50 is not None and n90 is not None:
            break

    return [file_path,contig_count,assembly_bases,n50,n90,l50,l90,sha256]

def fasta_to_bedtool(fasta_file):
    intervals = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom_name = record.id
        chrom_length = len(record.seq)
        interval = Interval(chrom_name, 0, chrom_length)
        intervals.append(interval)
    bedtool = BedTool(intervals).sort()
    return bedtool

def calculate_total_length(bedtool):
    return sum(len(interval) for interval in bedtool)

def compare_kmers(query_file,reference_file):
    ref_kmers = set()
    ref_process = subprocess.run(["kmercountexact.sh", f"in={reference_file}", "threads=1","fastadump=f", "out=stdout", "|", "cut", "-f1"], capture_output=True, text=True)
    ref_kmers.update(ref_process.stdout.strip().split('\n'))

    query_kmers = set()
    query_process = subprocess.run(["kmercountexact.sh", f"in={query_file}", "threads=1","fastadump=f", "out=stdout", "|", "cut", "-f1"], capture_output=True, text=True)
    query_kmers.update(query_process.stdout.strip().split('\n'))
    
    intersection = len(ref_kmers.intersection(query_kmers))
    similarity = 100*(intersection/(len(ref_kmers) + len(query_kmers) - intersection))
    unique_ref = len(ref_kmers.difference(query_kmers))
    unique_query = len(query_kmers.difference(ref_kmers))
    
    return [len(ref_kmers),len(query_kmers),
                unique_ref,unique_query,
                intersection,similarity]
    
#### 01: Read in arguments ####

query = str(sys.argv[1])
query_fasta = str(sys.argv[2])
reference = str(sys.argv[3])
reference_fasta = str(sys.argv[4])
mummer_dir = os.path.normpath(os.path.abspath(sys.argv[5]))
snpdiffs_dir = os.path.normpath(os.path.abspath(sys.argv[6]))

if sys.argv[7] != "":
    helpers.set_tempdir(sys.argv[7])

# Get query data
query_data = [query] + fasta_info(query_fasta)
query_string = [x+":"+str(y) for x,y in zip(['Query_ID','Query_Assembly','Query_Contig_Count','Query_Assembly_Bases',
                                             'Query_N50','Query_N90','Query_L50','Query_L90','Query_SHA256'],query_data)]

# Get reference data
reference_data = [reference] + fasta_info(reference_fasta)
reference_string = [x+":"+str(y) for x,y in zip(['Reference_ID','Reference_Assembly','Reference_Contig_Count','Reference_Assembly_Bases',
                                                 'Reference_N50','Reference_N90','Reference_L50','Reference_L90','Reference_SHA256'],reference_data)]

# Get kmer distance
[ref_kmers,query_kmers,
 unique_ref_kmers,unique_query_kmers,
 kmer_intersection,kmer_similarity] = compare_kmers(query_fasta,reference_fasta)

# Create reference BED file using fasta seq lengths
query_chr_bed = fasta_to_bedtool(query_fasta)
reference_chr_bed = fasta_to_bedtool(reference_fasta)

# Set report ID
report_id = query + "__vs__" + reference
snpdiffs_file = snpdiffs_dir + "/" + report_id + ".snpdiffs"

# Create NA variables for all downstream options
median_percent_identity = "NA"
median_alignment_length = "NA"

total_snp_count = "NA"
total_indel_count = "NA"
total_invalid_count = "NA"

#### 02: Read in MUMmer report data ####
[ref_bases,percent_ref_aligned,
            query_bases,percent_query_aligned,
            g_snps,g_indels,
            ref_breakpoints,query_breakpoints,
            ref_relocations,query_relocations,
            ref_translocations,query_translocations,
            ref_inversions,query_inversions,
            ref_insertions,query_insertions,
            ref_tandem,query_tandem] = parseMUmmerReport(mummer_dir,report_id)

if percent_ref_aligned > 0:
    
    #### 03: Process MUMmer coords file ####
    coords_file = parseMUmmerCoords(mummer_dir,report_id,reference_chr_bed,query_chr_bed)
    aligned_coords = coords_file[~((coords_file['Query_Contig'] == ".") | (coords_file['Ref_Contig'] == "."))]
    if aligned_coords.shape[0] > 0:
        aligned_coords['Perc_Iden'] = aligned_coords['Perc_Iden'].astype(float)
        aligned_coords['Ref_Aligned'] = aligned_coords['Ref_Aligned'].astype(int)

        median_percent_identity = np.median(aligned_coords['Perc_Iden'])
        median_alignment_length = np.median(aligned_coords['Ref_Aligned'])

        ##### 04: Process MUMmer SNP file ####
        processed_snps = parseMUmmerSNPs(mummer_dir,report_id,aligned_coords)
        
        total_snp_count = processed_snps[processed_snps['Cat'] == "SNP"].shape[0]
        total_indel_count = processed_snps[processed_snps['Cat'] == "Indel"].shape[0]
        total_invalid_count = processed_snps[processed_snps['Cat'] == "Invalid"].shape[0]
        
# Clean up pybedtools temp
helpers.cleanup(verbose=False, remove_all=False)

# Create header
percent_ref_aligned = f"{percent_ref_aligned:.2f}" if percent_ref_aligned != "NA" else percent_ref_aligned
percent_query_aligned = f"{percent_query_aligned:.2f}" if percent_query_aligned != "NA" else percent_query_aligned
median_percent_identity = f"{median_percent_identity:.2f}" if median_percent_identity != "NA" else median_percent_identity
median_alignment_length = f"{median_alignment_length:.2f}" if median_alignment_length != "NA" else median_alignment_length
total_snp_count = f"{total_snp_count:.0f}" if total_snp_count != "NA" else total_snp_count
total_indel_count = f"{total_indel_count:.0f}" if total_indel_count != "NA" else total_indel_count
total_invalid_count = f"{total_invalid_count:.0f}" if total_invalid_count != "NA" else total_invalid_count
ref_breakpoints = f"{ref_breakpoints:.0f}"
ref_relocations = f"{ref_relocations:.0f}"
ref_translocations = f"{ref_translocations:.0f}"
ref_inversions = f"{ref_inversions:.0f}"
ref_insertions = f"{ref_insertions:.0f}"
ref_tandem = f"{ref_tandem:.0f}"
query_breakpoints = f"{query_breakpoints:.0f}"
query_relocations = f"{query_relocations:.0f}"
query_translocations = f"{query_translocations:.0f}"
query_inversions = f"{query_inversions:.0f}"
query_insertions = f"{query_insertions:.0f}"
query_tandem = f"{query_tandem:.0f}"
g_snps = f"{g_snps:.0f}"
g_indels = f"{g_indels:.0f}"
ref_kmers = f"{ref_kmers:.0f}"
query_kmers = f"{query_kmers:.0f}"
unique_ref_kmers = f"{unique_ref_kmers:.0f}"
unique_query_kmers = f"{unique_query_kmers:.0f}"
kmer_intersection = f"{kmer_intersection:.0f}"
kmer_similarity = f"{kmer_similarity:.2f}" 

snpdiffs_header=[]
snpdiffs_header.append("#\t" +
                       "\t".join(query_string) +
                       "\t" + "\t".join(reference_string) +
                       "\t" + "\t".join([
"SNPs:"+total_snp_count,  
"Reference_Percent_Aligned:"+percent_ref_aligned,
"Query_Percent_Aligned:"+percent_query_aligned,
"Median_Percent_Identity:"+median_percent_identity,
"Median_Alignment_Length:"+median_alignment_length,
"Kmer_Similarity:"+kmer_similarity,    
"Shared_Kmers:"+kmer_intersection,
"Reference_Unique_Kmers:"+unique_ref_kmers,
"Query_Unique_Kmers:"+unique_query_kmers,
"Reference_Breakpoints:"+ref_breakpoints,
"Query_Breakpoints:"+query_breakpoints,
"Reference_Relocations:"+ref_relocations,
"Query_Relocations:"+query_relocations,
"Reference_Translocations:"+ref_translocations,
"Query_Translocations:"+query_translocations,
"Reference_Inversions:"+ref_inversions,
"Query_Inversions:"+query_inversions,
"Reference_Insertions:"+ref_insertions,
"Query_Insertions:"+query_insertions,
"Reference_Tandem:"+ref_tandem,
"Query_Tandem:"+query_tandem,
"Indels:"+total_indel_count,
"Invalid:"+total_invalid_count,
"gSNPs:"+g_snps,
"gIndels:"+g_indels]))

with open(snpdiffs_file,"w") as file:
    file.write("\n".join(snpdiffs_header) + "\n")
    for index, row in coords_file.iterrows():
        file.write("##\t" + "\t".join(map(str, row))+"\n")
    for index, row in processed_snps.iterrows():
        file.write("\t".join(map(str, row))+"\n")

print(",".join([query,reference,snpdiffs_file]))
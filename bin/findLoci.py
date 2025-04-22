#!/usr/bin/env python3

import sys
import os
import pandas as pd
import datetime
from pybedtools import BedTool,helpers
import concurrent.futures
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO,SeqIO,pairwise2
from Bio.pairwise2 import format_alignment
from itertools import combinations
import numpy as np
import uuid
import traceback
import shutil
import argparse
import textwrap

def safe_int(val):
    try:
        return int(val) if val != "NA" else 0
    except:
        return 0

def safe_float(val):
    try:
        return float(val) if val != "NA" else 0
    except:
        return 0

# SNPDiffs processing functions
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

        # Gather covered loci
        covered_bed_df = bed_df[(bed_df['Ref_Start'] != ".") & (bed_df['Query_Start'] != ".")].copy()
        if covered_bed_df.shape[0] > 0:
            for col in int_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float).astype(int)
            for col in float_columns:
                covered_bed_df.loc[:, col] = covered_bed_df.loc[:, col].astype(float)
        
        return covered_bed_df
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
    
    return_columns = ['Ref_Contig','Start_Ref','Ref_Pos',
            'Query_Contig','Start_Query','Query_Pos',
            'Ref_Loc','Query_Loc',
            'Ref_Start','Ref_End',
            'Query_Start','Query_End',
            'Ref_Base','Query_Base',
            'Dist_to_Ref_End','Dist_to_Query_End',
            'Ref_Aligned','Query_Aligned',
            'Perc_Iden','Cat']
    
    reverse_columns = ['Query_Contig','Start_Query','Query_Pos',
                    'Ref_Contig','Start_Ref','Ref_Pos',
                    'Query_Loc','Ref_Loc',
                    'Query_Start','Query_End',
                    'Ref_Start','Ref_End',
                    'Query_Base','Ref_Base',
                    'Dist_to_Query_End','Dist_to_Ref_End',
                    'Query_Aligned','Ref_Aligned',
                    'Query_Direction','Perc_Iden','Cat']
    
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

def calculate_total_length(bedtool):
    return sum(len(interval) for interval in bedtool)

# Return a dataframe of overlaps from MUMmer
def getLocusCoords(bed_df, min_len, min_iden):

    locus_columns = ['Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned',
                'Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned',
                'Perc_Iden','Overlap_Cat']
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
    
    covered_locus_df = bed_df.loc[(bed_df['Ref_Contig'] != ".") & (bed_df['Query_Contig'] != ".")].copy()
    
    # Filter out overlaps based on --min_len and --min_iden
    reject_length = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] < min_len) & (covered_locus_df['Perc_Iden'] >= min_iden)].copy()
    if reject_length.shape[0] > 0:
        reject_length['Overlap_Cat'] = "Purged_Length"
        
    reject_iden = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] >= min_len) & (covered_locus_df['Perc_Iden'] < min_iden)].copy()
    if reject_iden.shape[0] > 0:
        reject_iden['Overlap_Cat'] = "Purged_Identity"

    reject_lenIden = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] < min_len) & (covered_locus_df['Perc_Iden'] < min_iden)].copy()
    if reject_lenIden.shape[0] > 0:
        reject_lenIden['Overlap_Cat'] = "Purged_LengthIdentity"
    
    pass_filter = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] >= min_len) & (covered_locus_df['Perc_Iden'] >= min_iden)].copy().reset_index(drop=True)
    pass_filter['Overlap_Cat'] = "Pass"

    helpers.cleanup(verbose=False,remove_all = False)

    return_df = pd.concat([pass_filter,reject_length,reject_iden,reject_lenIden]).reset_index(drop=True)
    return return_df[locus_columns]

# This process assumes an in-frame protein-coding locus. Non-coding or operon functions to be added later
def processLocus(reference_id,snpdiffs_file):
    
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
    if (header_ref == reference_id):
        pass
    elif (header_query == reference_id):
        header_data = swapHeader(header_data)
    else:
        run_failed = True       
        sys.exit(f"Error: Reference ID not found in header of {snpdiffs_file}...")      

    # Get locus fasta file
    locus_fasta = header_data['Reference_Assembly'][0]
    
    # Assert locus contains a single record
    locus_assembly = SeqIO.to_dict(SeqIO.parse(locus_fasta, "fasta"))
    
    # For now, one > per FASTA, but this could be a good place to work in some redundancy
    assert len(locus_assembly) == 1, f"Locus file {locus_fasta} contains more than one record"
    locus_seq_record = next(iter(locus_assembly.values()))
    
    # Assert the locus is in-frame and protein coding (Can extend for non-coding loci later)
    locus_length = len(locus_seq_record.seq)   
    assert locus_length % 3 == 0, f"Locus file {locus_fasta} does not appear to have codons (not an even multiple of 3)"
 
    # Perform translation
    positive_translation = locus_seq_record.seq.translate(to_stop=False)
    
    # Check if positive translation is in frame and terminal, then reverse
    if positive_translation.count("*") == 1 and positive_translation.endswith("*"):
        locus_type = "terminal_forward"
        locus_record = SeqRecord(locus_seq_record.seq, id=reference_id, name=reference_id, description=f"{reference_id}_terminal_forward")
        protein_record = SeqRecord(positive_translation, id=reference_id, name=reference_id, description=f"{reference_id}_terminal_forward")

    elif positive_translation.count("*") == 0:
        locus_type = "nonterminal_forward"
        locus_record = SeqRecord(locus_seq_record.seq, id=reference_id, name=reference_id, description=f"{reference_id}_nonterminal_forward")
        protein_record = SeqRecord(positive_translation, id=reference_id, name=reference_id, description=f"{reference_id}_nonterminal_forward")

    else:
        reverse_translation = locus_seq_record.seq.reverse_complement().translate(to_stop=False)
        if reverse_translation.count("*") == 1 and reverse_translation.endswith("*"):
            locus_type = "terminal_reverse"
            locus_record = SeqRecord(locus_seq_record.seq.reverse_complement(), id=reference_id, name=reference_id, description=f"{reference_id}_terminal_reverse")
            protein_record = SeqRecord(reverse_translation, id=reference_id, name=reference_id, description=f"{reference_id}_terminal_reverse")
        elif reverse_translation.count("*") == 0:
            locus_type = "nonterminal_reverse"
            locus_record = SeqRecord(locus_seq_record.seq.reverse_complement(), id=reference_id, name=reference_id, description=f"{reference_id}_nonterminal_reverse")
            protein_record = SeqRecord(reverse_translation, id=reference_id, name=reference_id, description=f"{reference_id}_nonterminal_reverse")
        else:
            locus_type = "noncoding"
            locus_record = SeqRecord(locus_seq_record.seq, id=reference_id, name=reference_id, description=f"{reference_id}_noncoding_forward")
            protein_record = SeqRecord(Seq(""),id=reference_id, name=reference_id, description=f"{reference_id}_noncoding_forward")    

    return(locus_record,protein_record,locus_type,locus_fasta)

# This function gets the query FASTA from the snpdiffs file, and if there is a MUMmer hit, it returns the hit with 5' and  3' extention if needed/possible
def processQuery(snpdiffs_file,trim_name, min_cov, min_len, min_iden, max_contigs,reference_id,locus_record,protein_record,log_directory,output_directory):
    
    #### Initialize unprocessed data ####
    query_record = SeqRecord(Seq(""), id="", description="",name="")
    query_type = "Unprocessed"
    start_extension = 0
    end_extension = 0
    start_extension_type = "Unextended"
    end_extension_type = "Unextended"
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
    
    #### Process SNPDiffs data ####
    if not os.path.exists(snpdiffs_file) or not snpdiffs_file.endswith('.snpdiffs'):
        sys.exit(f"Invalid snpdiffs file provided: {snpdiffs_file}")
        
    try:
        header_data = fetchHeaders(snpdiffs_file)
        header_query = header_data['Query_ID'][0].replace(trim_name,'')
        header_ref = header_data['Reference_ID'][0].replace(trim_name,'')
    except:
        sys.exit(f"Error reading headers from snpdiffs file: {snpdiffs_file}")
    
    if (header_ref == reference_id):
        snpdiffs_orientation = 1
    elif (header_query == reference_id):
        snpdiffs_orientation = -1
        header_data = swapHeader(header_data)
    else:
        sys.exit(f"Error: Reference ID not found in header of {snpdiffs_file}...")     
    
    query_id = header_data['Query_ID'][0]
    
    fasta_name = f"{query_id}_{reference_id}"
    query_bases = safe_int(header_data['Query_Assembly_Bases'][0])
    query_contigs = safe_int(header_data['Query_Contig_Count'][0])
    query_assembly = header_data['Query_Assembly'][0]
    locus_assembly = header_data['Reference_Assembly'][0]
    query_coords = [query_id,reference_id,".",np.nan,np.nan]

    #### Initialize log file ####
    log_file = f"{log_directory}/{query_id}__vs__{reference_id}.log"

    if os.path.exists(query_assembly):
        with open(log_file,"w+") as log:
            log.write("CSP2 Locus Finder Analysis Log\n")
            log.write(f"Query Isolate: {query_id}\n")
            log.write(f"Query Assembly: {query_assembly}\n")
            log.write(f"Query Contigs: {query_contigs}; Query Assembly Size: {query_bases}\n\n")
            log.write(f"Reference Locus: {reference_id}\n")
            log.write(f"Locus FASTA: {locus_assembly}\n")
            log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
            log.write("-------------------------------------------------------\n\n")

    else:
        query_type = "Missing_Assembly"
        with open(log_file,"a+") as log:
            log.write(f"Query assembly {query_assembly} does not exist! Locus finding halted...\n")
            log.write("-------------------------------------------------------\n\n")
        return(query_id,query_type,query_contigs,query_bases,query_record,query_coords)
    
    # If the reference is not covered, STOP
    raw_locus_percent_aligned = safe_float(header_data['Reference_Percent_Aligned'][0])
    if raw_locus_percent_aligned == 0:
        query_type = "No_Coverage"
        with open(log_file,"a+") as log:
            log.write(f"MUMmer detected no overlaps between {query_id} and locus...Locus finding halted...\n")
            log.write("-------------------------------------------------------\n\n")
        return(query_id,query_type,query_contigs,query_bases,query_record,query_coords)

    # If the reference is covered, continue
    with open(log_file,"a+") as log:
        log.write("MUMmer overlap detected, getting coordinates from 1coords...")
    try:
        bed_df,snp_df = parseSNPDiffs(snpdiffs_file,snpdiffs_orientation)
        with open(log_file,"a+") as log:
            log.write("Done!\n")
            log.write("-------------------------------------------------------\n\n")
    except Exception as e:
        with open(log_file,"a+") as log:
            log.write(f"\nError reading BED/SNP data from file: {snpdiffs_file}\n{str(e)}")
        sys.exit(f"Error reading BED/SNP data from file: {snpdiffs_file}\n{str(e)}")
    
    # Ensure that overlaps are covered by min_len with at least min_iden
    with open(log_file,"a+") as log:
        log.write("Filtering for short overlaps and low percent identity...")
    
    good_bed_df = bed_df[(bed_df['Ref_Aligned'] >= min_len) & (bed_df['Perc_Iden'] >= min_iden)].copy()

    if good_bed_df.shape[0] == 0:
        query_type = "Low_Quality_Coverage"
        with open(log_file,"a+") as log:
            log.write(f"\n\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}) , no valid alignments remain...Locus finding halted...\n")
            log.write("-------------------------------------------------------\n\n")
        return(query_id,query_type,query_contigs,query_bases,query_record,query_coords)

    # If overlaps pass QC, extract sequence data from query assembly
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write("-------------------------------------------------------\n\n")

    with open(log_file,"a+") as log:
        log.write("Getting FASTA coordinates for locus...")    
    
    locus_coords_df = getLocusCoords(bed_df, min_len, min_iden)
    
    # Not handling multicontigs now
    if max_contigs > 1:
        with open(log_file,"a+") as log:
            log.write(f"\n-\t --max_contigs is set to a value higher than one ({max_contigs}). This feature is not available yet, sorry! Setting to 1...\n")
        max_contigs = 1
        
    # If the locus requires more than max_contigs, STOP
    if locus_coords_df['Query_Contig'].nunique() > max_contigs:
        unique_contigs = locus_coords_df['Query_Contig'].unique()
        with open(log_file,"a+") as log:
            log.write(f"\n-\t The locus {reference_id} covers multiple contigs in {query_id}\n")
            for contig in unique_contigs:
                log.write(f"\t- {contig}\n")
        query_record = SeqRecord(Seq(""), id=fasta_name, description="Multicontig",name=fasta_name)
        query_type = "Fail_Max_Contigs"
        return(query_id,query_type,query_contigs,query_bases,query_record,query_coords)
    
    # If the locus hit is contained within a single contig (for now), pull out the hit
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        log.write(" Extracting sequence data from query assembly...")

    # If the locus aligns to a single contig, process one or more hits
    assembly = SeqIO.to_dict(SeqIO.parse(query_assembly, "fasta"))    
    chrom = locus_coords_df['Query_Contig'][0]
    if chrom not in assembly:
        sys.exit(f"{chrom} not found in assembly {query_assembly}.")
        
    query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type = extendHit(query_id,reference_id,locus_coords_df,chrom,assembly[chrom].seq,locus_record,protein_record,fasta_name)
    
    # Save results to the log
    with open(log_file,"a+") as log:
        log.write("Done!\n")
        if query_type == "Overlapping":
            log.write("\t- Query has overlapping intervals, unable to extract sequence.\n")
        elif query_type == "Unknown_Frame":
            log.write("\t- Query has a weak frame signature, returning the full hit for inspection.\n")
        elif query_type == "Inversion":
            log.write("\t- Multiple orientations in hit, possible inversion, returning the full hit for inspection.\n")
        elif query_type == "Translocation":
            log.write("\t- Possible translocation, returning the full hit for inspection.\n")
        
        # Hits with analyzable FASTAs
        elif query_type == "Complete_Hit":
            log.write("\t- Query had a complete hit to the locus.\n")
        
        elif query_type == "Start_to_EndContig":
            log.write(f"\t- Query was extended {end_extension} to the end of the contig.\n")            
        elif query_type == "Start_to_ExtStop":
            log.write(f"\t- Query was extended {end_extension}bp to the first in-frame stop codon.\n")     
        elif query_type == "Start_to_Premature_Stop":
            log.write(f"\t- Query hit was reduced by {end_extension}bp due to a premature stop codon.\n")
            
        elif query_type == "ExtStop_to_Stop":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the first in-frame stop codon.\n")
        elif query_type == "ExtStop_to_ExtStop":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the first in-frame stop codon, and was extended {end_extension}bp to the first in-frame stop codon.\n")
        elif query_type == "ExtStop_to_EndContig":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the first in-frame stop codon, and was extended {end_extension}bp to the end of the contig.\n")
        elif query_type == "ExtStop_to_Premature_Stop":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the first in-frame stop codon, and was reduced {end_extension}bp to a premature stop codon.\n")
            
        elif query_type == "StartContig_to_Stop":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the end of the contig.\n")            
        elif query_type == "StartContig_to_ExtStop":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the end of the contig, and was extended {end_extension}bp to the first in-frame stop codon.\n")
        elif query_type == "StartContig_to_Premature_Stop":
            log.write(f"\t- Query was extended {start_extension}bp in the 5' direction to the end of the contig, and was reduced {end_extension}bp to a premature stop codon.\n")
        elif query_type == "Non_Terminal_Contig":
            log.write("\t- The entire query contig was in frame and no stop codon was detected, returning the in-frame contig.\n")

    return(query_id,query_type,query_contigs,query_bases,query_record,query_coords)
    
def extendHit(query_id,reference_id,locus_coords_df,contig_id,contig_seq,locus_record,protein_record,fasta_name):
    
    # Establish start and end extensions
    start_extension = 0
    start_extension_type = "Unextended"
    
    end_extension = 0
    end_extension_type = "Unextended"
        
    # Get length of query contig and locus
    contig_length = len(contig_seq)
    locus_length = len(locus_record.seq)
    query_coords = [query_id,reference_id,contig_id,0,contig_length]
    
    # Determine the min and max overlaps of the locus
    ref_start = locus_coords_df['Ref_Start'].min()
    ref_end = locus_coords_df['Ref_End'].max()   
    
    # Get boolean for coverage of the 5' and 3' end
    start_cov = True if ref_start == 0 else False
    end_cov = True if ref_end == locus_length else False
    
    # Get the lowest/highest values for start/end and grab the entire hit
    ultimate_min = min(locus_coords_df['Query_Start'].min(), locus_coords_df['Query_End'].min())
    ultimate_max = max(locus_coords_df['Query_Start'].max(), locus_coords_df['Query_End'].max())
    min_max_hit = contig_seq[ultimate_min:ultimate_max]
    
    # Add the sequence for each hit
    locus_coords_df['Extracted_Seq'] = locus_coords_df.apply(lambda row: ''.join(contig_seq[row['Query_Start']:row['Query_End']]),axis=1)

    # Get the orientation for each hit
    locus_coords_df['Orientation'] = locus_coords_df['Extracted_Seq'].apply(lambda seq: checkOrientation(seq, locus_record.seq))
    
    # Ensure the orientation for each sequence is the same
    if locus_coords_df.shape[0] == 1:
        pass
    else:
        # If the orientation varies between hits, call it an inversion and report the min/max hit
        if locus_coords_df['Orientation'].nunique() > 1:
            query_record = SeqRecord(min_max_hit, id=fasta_name, description="Inversion",name=fasta_name)
            query_type = "Inversion"
            query_coords = [query_id,reference_id,contig_id,ultimate_min,ultimate_max]
            return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)

        # Check for translocations, and return min/max hit if found
        query_order = locus_coords_df.sort_values('Ref_Start')['Query_Start'].rank(method="first").tolist()
        expected_order = list(range(1, len(query_order)+1))
        reversed_order = expected_order[::-1]

        if query_order != expected_order and query_order != reversed_order:
            query_record = SeqRecord(min_max_hit, id=fasta_name, description="Translocation", name=fasta_name)
            query_type = "Translocation"
            query_coords = [query_id,reference_id,contig_id,ultimate_min,ultimate_max]
            return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
    
    # Get the hit orientation
    hit_orientation = locus_coords_df['Orientation'].values[0]

    lowest_hit = locus_coords_df[locus_coords_df['Ref_Start'] == ref_start]    
    highest_hit = locus_coords_df[locus_coords_df['Ref_End'] == ref_end]

    # Only one hit should map to the front and back, return empty seq
    if lowest_hit.shape[0] > 1 or highest_hit.shape[0] > 1:
        query_type = "Overlapping"
        query_record = SeqRecord(Seq(""), id=fasta_name, description="Overlapping",name=f"{fasta_name}")
        return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
    
    # Reverse the contig and start/end if in the - orientation
    if hit_orientation == "+":
        low_start = lowest_hit['Query_Start'].values[0]
        low_end = lowest_hit['Query_End'].values[0]
        high_end = highest_hit['Query_End'].values[0]
    else:
        
        contig_seq = contig_seq.reverse_complement()
    
        low_start = contig_length - lowest_hit['Query_End'].values[0]
        low_end = contig_length - lowest_hit['Query_Start'].values[0]
        high_end = contig_length - highest_hit['Query_Start'].values[0]
        

       
    # Get full hit
    oriented_hit = contig_seq[low_start:high_end]
    
    # If the hit is complete, return the hit if a multiple of 3, otherwise return as Unknown_Frame
    if start_cov and end_cov:
        if len(oriented_hit) % 3 == 0:
            query_type = "Complete_Hit"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="Complete_Hit",name=fasta_name)
            if hit_orientation == "+":
                query_coords = [query_id,reference_id,contig_id,low_start,high_end]
            else:
                query_coords = [query_id,reference_id,contig_id,contig_length - high_end,contig_length - low_start]
            return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
        else:
            query_type = "Unknown_Frame"
            query_record = SeqRecord(min_max_hit, id=fasta_name, description="Unknown_Frame",name=fasta_name)
            query_coords = [query_id,reference_id,contig_id,ultimate_min,ultimate_max]
            return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
    
    # If either the start/stop is not covered, try to extend starting from the 5' end to a stop
    oriented_5_hit = contig_seq[low_start:low_end]

    # Get best best frame for the 5' hit
    scores = []
    for frame in range(3):
        nuc = oriented_5_hit[frame:]
        nuc = nuc[:len(nuc) - (len(nuc) % 3)]
        aa_seq = nuc.translate(to_stop=False)
        aln = pairwise2.align.localms(protein_record.seq, str(aa_seq),2, -1, -2, -0.5,one_alignment_only=True)[0]
        scores.append({'frame': frame, 'score': aln.score,'aa_seq':aa_seq})
    scores.sort(key=lambda x: x['score'], reverse=True)
    
    # If one frame does not appear to match the protein, set as Unknown_Frame and return
    if scores[0]['score'] < 1.3 * scores[1]['score']:
        query_type = "Unknown_Frame"
        query_record = SeqRecord(min_max_hit, id=fasta_name, description="Unknown_Frame",name=fasta_name)
        query_coords = [query_id,reference_id,contig_id,ultimate_min,ultimate_max]
        return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
    
    # If one frame performs much better, grab the frame and amino acid
    best_frame = scores[0]['frame']
    best_aa = scores[0]['aa_seq']

    # Get the in-frame start/stop coordinates for the first and last codon
    codon_start = low_start + best_frame
    codon_end = codon_start + (len(best_aa) * 3)
    
    # If the amino acid for the 5' hit contains a *, return the hit up to the stop codon
    if "*" in best_aa:
        stop_index = best_aa.find("*")        
        stop_position = codon_start + (3*(stop_index+1))
        end_extension = codon_end - stop_position
        if end_cov and abs(stop_position - high_end) < 3:
            end_extension_type = "unnecessary"
        else:
            end_extension_type = "Premature_Stop"
    else:
        framed_downstream = contig_seq[codon_end:]
        if len(framed_downstream) < 3:
            end_extension_type = "Contig"
        else:
            trimmed_downstream = framed_downstream[:len(framed_downstream) - (len(framed_downstream) % 3)]
            translated_downstream = trimmed_downstream.translate(to_stop=False)
            
            # If there is a * in the translated downstream, extend to *
            if "*" in translated_downstream:
                downstream_stop_index = translated_downstream.find("*")
                end_extension = (downstream_stop_index + 1) * 3
                if abs((codon_end + end_extension) - high_end) < 3: 
                    end_extension_type = "unnecessary"
                else:
                    end_extension_type = "Stop"
            
            # If there is no * in the translated downstream, extend to end of contig
            else:
                end_extension = len(translated_downstream) * 3
                end_extension_type = "Contig"

    #### 5' Extension #####

    # Only extend 5' if necessary
    if start_cov and best_frame == 0:
        start_extension_type = "unnecessary"
    elif codon_start < 3:
        if not start_cov:
            start_extension_type = "Contig"
        else:
            # Odd case...if the frame isn't 0, then the start isn't lined up. Return to min/max investigate.
            query_type = "Unknown_Frame"
            query_record = SeqRecord(min_max_hit, id=fasta_name, description="Unknown_Frame",name=fasta_name)
            query_coords = [query_id,reference_id,contig_id,ultimate_min,ultimate_max]
            return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
    
    else:
        # If the alignment doesn't go the 5' end, or if the frame is off, extend backwards until first stop codon (MAY OVERSHOOT)
        framed_upstream = contig_seq[:codon_start]
        trimmed_upstream = framed_upstream[len(framed_upstream) % 3:]
        codons = ''.join(textwrap.wrap(str(trimmed_upstream), 3)[::-1])
        translated_upstream = Seq(str(codons)).translate(to_stop=False)
        if "*" in translated_upstream:
            upstream_stop_index = translated_upstream.find("*")
            start_extension = 3 * upstream_stop_index
            start_extension_type = "Stop"
        else:
            start_extension_type = "Contig"
            start_extension = 3 * len(translated_upstream)
        
    new_start = codon_start - start_extension
    new_end = codon_end + end_extension
    oriented_hit = contig_seq[new_start:new_end]    
    assert (new_end - new_start) % 3 == 0, f"The extended start and stop positions for {query_id} are not a multiple of three"

    # Set query_coords
    if hit_orientation == "+":
        query_coords = [query_id,reference_id,contig_id,new_start,new_end]
    else:
        query_coords = [query_id,reference_id,contig_id,contig_length - new_end,contig_length - new_start]
        
    # Set query_type based on start/stop extension
    if start_extension_type == "unnecessary":
        
        if end_extension_type == "unnecessary":
            query_type = "Complete_Hit"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="Complete_Hit",name=fasta_name)
    
        elif end_extension_type == "Stop":
            query_type = "Start_to_ExtStop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="Start_to_ExtStop",name=fasta_name)
            
        elif end_extension_type == "Contig":
            query_type = "Start_to_EndContig"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="Start_to_EndContig",name=fasta_name)
    
        elif end_extension_type == "Premature_Stop":
            query_type = "Start_to_Premature_Stop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="Start_to_Premature_Stop",name=fasta_name)
        
    elif start_extension_type == "Stop":
        
        if end_extension_type == "unnecessary":
            query_type = "ExtStop_to_Stop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="ExtStop_to_Stop",name=fasta_name)
    
        elif end_extension_type == "Stop":
            query_type = "ExtStop_to_ExtStop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="ExtStop_to_ExtStop",name=fasta_name)
            
        elif end_extension_type == "Contig":
            query_type = "ExtStop_to_EndContig"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="ExtStop_to_EndContig",name=fasta_name)
            
        elif end_extension_type == "Premature_Stop":
            query_type = "ExtStop_to_Premature_Stop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="ExtStop_to_Premature_Stop",name=fasta_name)
        
    elif start_extension_type == "Contig":
        
        if end_extension_type == "unnecessary":
            query_type = "StartContig_to_Stop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="StartContig_to_Stop",name=fasta_name)
    
        elif end_extension_type == "Stop":
            query_type = "StartContig_to_ExtStop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="StartContig_to_ExtStop",name=f"{fasta_name}_StartContig_to_ExtStop")
            
        elif end_extension_type == "Contig":
            query_type = "Non_Terminal_Contig"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="Non_Terminal_Contig",name=f"{fasta_name}_Non_Terminal_Contig")
        
        elif end_extension_type == "Premature_Stop":
            query_type = "StartContig_to_Premature_Stop"
            query_record = SeqRecord(oriented_hit, id=fasta_name, description="StartContig_to_Premature_Stop",name=f"{fasta_name}_Premature_Stop")
    
    return(query_type,query_record,query_coords,start_extension,start_extension_type,end_extension,end_extension_type)
             
def checkOrientation(query_seq, locus_seq):
    
    fwd = pairwise2.align.localms(str(locus_seq).upper(), str(query_seq).upper(), 2, -1, -2, -0.5, one_alignment_only=True)[0].score
    rev = pairwise2.align.localms(str(locus_seq).upper(), str(Seq(query_seq).reverse_complement()).upper(), 2, -1, -2, -0.5, one_alignment_only=True)[0].score
    
    return '+' if fwd >= rev else '-'   

def saveAbberantAlignment(reference_id, locus_record, query_id, query_record, aberrant_directory):

    alignments = pairwise2.align.globalms(
        query_record.seq,
        locus_record.seq,
        match=2,           
        mismatch=-1,      
        open=-10,     
        extend=-1, 
        one_alignment_only=True)

    if alignments:
        best_alignment = format_alignment(*alignments[0])
        out_file = os.path.join(aberrant_directory, f"{query_id}_vs_{reference_id}.aln")
        with open(out_file, "w") as f:
            f.write(best_alignment)
            
#### MAIN SCRIPT ####
    
# Read in arguments
global run_failed
run_failed = False

start_time = time.time()

parser = argparse.ArgumentParser(description='CSP2 Locus Finder')
parser.add_argument('--reference_id', type=str, help='Reference Locus')
parser.add_argument('--output_directory', type=str, help='Output Directory')
parser.add_argument('--log_directory', type=str, help='Log Directory')
parser.add_argument('--snpdiffs_file', type=str, help='Path to SNPdiffs file')
parser.add_argument('--min_cov', default=85,type=float, help='Minimum coverage')
parser.add_argument('--min_len', default=500,type=int, help='Minimum length')
parser.add_argument('--min_iden', default=99,type=float, help='Minimum identity')
parser.add_argument('--max_contigs',default=1,type=int,help="Maximum contigs containing a single locus")
parser.add_argument('--trim_name', nargs='?', const="", default="", type=str, help='Trim name')
parser.add_argument('--tmp_dir',default="", type=str, help='Temporary directory')
args = parser.parse_args()

reference_id = args.reference_id
output_directory = os.path.abspath(args.output_directory)
timestamp = time.strftime("%Y%m%d_%H%M%S")

alignment_directory = f"{output_directory}/Alignments"
if os.path.exists(alignment_directory):
    alignment_directory = f"{alignment_directory}_{timestamp}"

aberrant_directory = f"{output_directory}/Aberrant_Sequences"
if os.path.exists(aberrant_directory):
    aberrant_directory = f"{aberrant_directory}_{timestamp}"

valid_directory = f"{output_directory}/Valid_Query_Sequences"
if os.path.exists(valid_directory):
    valid_directory = f"{valid_directory}_{timestamp}"
    
log_directory = os.path.abspath(args.log_directory)

reference_screening_file = f"{output_directory}/Locus_Screening.tsv"
if os.path.exists(reference_screening_file):
    reference_screening_file = f"{output_directory}/Locus_Screening_{timestamp}.tsv"

hit_coords_file = f"{output_directory}/Hit_Coordinates.tsv"
if os.path.exists(hit_coords_file):
    hit_coords_file = f"{output_directory}/Hit_Coordinates_{timestamp}.tsv"

log_file = f"{output_directory}/CSP2_SNP_Pipeline.log"
if os.path.exists(log_file):
    log_file = f"{output_directory}/CSP2_SNP_Pipeline_{timestamp}.log"
    
snpdiffs_file = args.snpdiffs_file
min_cov = args.min_cov
min_len = args.min_len
min_iden = args.min_iden
trim_name = args.trim_name
max_contigs = args.max_contigs

# Establish log file
with open(log_file,"w+") as log:
    log.write("CSP2 Locus Finder Analysis\n")
    log.write(f"Reference Locus: {reference_id}\n")
    log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
    log.write("-------------------------------------------------------\n\n")
    log.write("Reading in SNPDiffs files...")

# Read in all lines and ensure each file exists
snpdiffs_list = [line.strip() for line in open(snpdiffs_file, 'r')]
snpdiffs_list = [line for line in snpdiffs_list if line]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        run_failed = True
        sys.exit("Error: File does not exist: " + snpdiffs_file)

snpdiffs_list = list(set(snpdiffs_list))

if len(snpdiffs_list) == 0:
    run_failed = True
    sys.exit("No SNPdiffs files provided...")

with open(log_file, "a+") as log:
    log.write("Done!\n")
    log.write(f"\t- Read in {len(snpdiffs_list)} SNPdiffs files\n")
    log.write("-------------------------------------------------------\n\n")

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
    # Process the locus information
    with open(log_file,"a+") as log:
        log.write("Processing locus data...")
    
    # Get the locus data from the first SNPDiffs file
    locus_record,protein_record,locus_type,locus_fasta = processLocus(reference_id,snpdiffs_list[0])

    with open(log_file, "a+") as log:
        log.write("Done!\n")
        log.write(f"\t- Locus is determined to be {locus_type}\n")
        log.write("-------------------------------------------------------\n\n")

    if locus_type != "noncoding":
        with open(log_file,"a+") as log:
            log.write("Searching all query assemblies for loci...")
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(processQuery,snp_diff_file,trim_name, min_cov, min_len, min_iden, max_contigs,reference_id,locus_record,protein_record,log_directory,output_directory) for snp_diff_file in snpdiffs_list]
        with open(log_file,"a+") as log:
            log.write("Done!\n")
    
        results = []
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()  # Unpacks the return value of processQuery
                query_id,query_type,query_contigs,query_bases,query_record,query_coords = result
                results.append((query_id,query_type,query_contigs,query_bases,query_record,query_coords))        
            except Exception as exc:
                traceback.print_exc()
                sys.exit(1)

        df = pd.DataFrame([{
            "Reference_Locus": reference_id,
            "Query_ID": query_id,
            "Query_Type": query_type,
            "Query_Contigs": query_contigs,
            "Query_Bases": query_bases,
            "Query_Contig": query_coords[2],
            "Query_Start": query_coords[3],
            "Query_End": query_coords[4],
            "Query_Hit_Length": len(query_record.seq)} for query_id, query_type, query_contigs, query_bases, query_record, query_coords in results])
        
        # Save output
        df.to_csv(reference_screening_file, sep="\t", index=False)
        
        # Process aberrant sequences if present
        ind_aligners = ['Inversion','Translocation','Unknown_Frame']
        group_aligners = ["Complete_Hit", "ExtStop_to_EndContig", "ExtStop_to_ExtStop", "ExtStop_to_Premature_Stop", "ExtStop_to_Stop", "Non_Terminal_Contig", "StartContig_to_ExtStop", "StartContig_to_Premature_Stop", "StartContig_to_Stop", "Start_to_EndContig", "Start_to_ExtStop", "Start_to_Premature_Stop"]

        individual_alignment_df = df[df['Query_Type'].isin(ind_aligners)]     
        group_align_df = df[df['Query_Type'].isin(group_aligners)]

        # Create valid/aberrant directory if needed
        if individual_alignment_df.shape[0] > 0:
            os.makedirs(aberrant_directory)
        if group_align_df.shape[0] > 0:
            os.makedirs(alignment_directory)
            os.makedirs(valid_directory)
            
        # Save FASTA files and gather valid_records
        valid_records = []
        for query_id, query_type, query_contigs, query_bases, query_record,query_coords in results:
            if str(query_record.seq) != "":
                if query_type in ind_aligners:
                    SeqIO.write(query_record, f"{aberrant_directory}/{query_id}_{query_type}.fasta", "fasta")
                    saveAbberantAlignment(reference_id,locus_record,query_id,query_record,aberrant_directory)            
                elif query_type in group_aligners:
                    SeqIO.write(query_record, f"{valid_directory}/{query_id}.fasta", "fasta")
                    valid_records.append(query_record)
                else:
                    print(f"WHAT IS THIS {query_type}")
                    
        # Save a multifasta for alignment
        if len(valid_records) > 1:
            valid_records.append(locus_record)
            multifasta_file = f"{alignment_directory}/{reference_id}.fasta"
            with open(multifasta_file,"w") as file:
                for record in valid_records:
                    SeqIO.write(record, file, "fasta")
            print(multifasta_file,end="")
            
    # Process non-coding genes here in the future
    else:
        with open(log_file,"a+") as log:
            log.write("For now, locus mode currently supports only protein-coding loci. Locus finding halting ....")
            
except:
    run_failed = True
    sys.exit("Exception occurred:\n", traceback.format_exc())
finally:
    helpers.cleanup(verbose=False, remove_all=False)
    if temp_dir != "":
        shutil.rmtree(temp_dir)
    if run_failed:
        sys.exit(1)
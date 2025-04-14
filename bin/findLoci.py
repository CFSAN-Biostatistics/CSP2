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
from itertools import combinations
import numpy as np
import uuid
import traceback
import shutil
import argparse

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

def checkOrientation(record,locus_record):
    
    ref_seq = str(locus_record.seq).upper()
    query_seq = str(record.seq).upper()
    query_rc = str(record.seq.reverse_complement()).upper()

    fwd_alignment = pairwise2.align.globalxx(ref_seq, query_seq, one_alignment_only=True)[0]
    rev_alignment = pairwise2.align.globalxx(ref_seq, query_rc, one_alignment_only=True)[0]
    
    if rev_alignment.score > fwd_alignment.score:
        return SeqRecord(seq=record.seq.reverse_complement(),id = record.name,name=record.name,description=record.description + "_REVCOMP")
    else:
        return SeqRecord(seq=record.seq,name=record.name,id=record.name,description=record.description)
    
def getFasta(locus_coords_df, query_assembly,query_id,reference_id,locus_record):
    
    fasta_name = f"{query_id}_{reference_id}"
    # If the locus is entirely on one query contig, find the min/max and pull out the entire overlapping region
    if locus_coords_df['Query_Contig'].nunique() > 1:
        print("Multicontig",flush=True)
    elif locus_coords_df.shape[0] > 1:
        print("Multi-interval",flush=True)
    else:    
        assembly = SeqIO.to_dict(SeqIO.parse(query_assembly, "fasta"))
        
        chrom = locus_coords_df['Query_Contig'][0]
        start = int(locus_coords_df['Query_Start'][0])

        end = int(locus_coords_df['Query_End'][0])
        if chrom not in assembly:
            sys.exit(f"{chrom} not found in assembly {query_assembly}.")

        seq = assembly[chrom].seq[start:end]
        record = SeqRecord(Seq(str(seq)), name =fasta_name,id=fasta_name,description=f"{chrom}:{start}-{end}")
    
        return(checkOrientation(record,locus_record))

def getLocusCoords(bed_df,log_file, min_len, min_iden):

    locus_columns = ['Ref_Contig','Ref_Start','Ref_End','Ref_Length','Ref_Aligned',
                'Query_Contig','Query_Start','Query_End','Query_Length','Query_Aligned',
                'Perc_Iden','Overlap_Cat']
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)
    
    covered_locus_df = bed_df.loc[(bed_df['Ref_Contig'] != ".") & (bed_df['Query_Contig'] != ".")].copy()
    # Filter out overlaps based on --min_len and --min_iden
    reject_length = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] < min_len) & (covered_locus_df['Perc_Iden'] >= min_iden)].copy()
    if reject_length.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Overlap Length): {reject_length.shape[0]}\n")
        reject_length['Overlap_Cat'] = "Purged_Length"
        
    reject_iden = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] >= min_len) & (covered_locus_df['Perc_Iden'] < min_iden)].copy()
    if reject_iden.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Identity): {reject_iden.shape[0]}\n")
        reject_iden['Overlap_Cat'] = "Purged_Identity"

    reject_lenIden = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] < min_len) & (covered_locus_df['Perc_Iden'] < min_iden)].copy()
    if reject_lenIden.shape[0] > 0:
        with open(log_file,"a+") as log:
            log.write(f"\t\t- Purged (Alignment Length + Identity): {reject_lenIden.shape[0]}\n")
        reject_lenIden['Overlap_Cat'] = "Purged_LengthIdentity"
    
    pass_filter = covered_locus_df.loc[(covered_locus_df['Ref_Aligned'] >= min_len) & (covered_locus_df['Perc_Iden'] >= min_iden)].copy().reset_index(drop=True)
    pass_filter['Overlap_Cat'] = "Pass"

    helpers.cleanup(verbose=False,remove_all = False)

    return_df = pd.concat([pass_filter,reject_length,reject_iden,reject_lenIden]).reset_index(drop=True)
    return return_df[locus_columns]
    
def findLoci(snpdiffs_file,trim_name, min_cov, min_len, min_iden, max_contigs,reference_id,log_directory,output_directory):
    
    query_record = SeqRecord(Seq(""), id="", name="", description="")

    screen_start_time = time.time()
    
    if temp_dir != "":
        helpers.set_tempdir(temp_dir)

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
        snpdiffs_orientation = 1
        query_id = header_query
    elif (header_query == reference_id):
        snpdiffs_orientation = -1
        query_id = header_ref
        header_data = swapHeader(header_data)
    else:
        run_failed = True       
        sys.exit(f"Error: Reference ID not found in header of {snpdiffs_file}...")      

    # Establish log file
    log_file = f"{log_directory}/{query_id}__vs__{reference_id}.log"
    with open(log_file,"w+") as log:
        log.write("CSP2 Locus Finder Analysis Log\n")
        log.write(f"Query Isolate: {query_id}\n")
        log.write(f"Reference Locus: {reference_id}\n")
        log.write(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n")
        log.write("-------------------------------------------------------\n\n")
        if snpdiffs_orientation == 1:
            log.write("\t- SNPDiffs file is in the forward orientation\n")
            log.write("-------------------------------------------------------\n\n")
        else:
            log.write("\t- SNPDiffs file is in the reverse orientation\n")
            log.write("-------------------------------------------------------\n\n")

    
    # Set variables from header data
    query_assembly = header_data['Query_Assembly'][0]
    locus_sequence = header_data['Reference_Assembly'][0]
    
    # Assert that query assembly and locus file exist
    assert os.path.exists(query_assembly), f"Query assembly not found: {query_assembly}"
    assert os.path.exists(locus_sequence), f"Locus sequence not found: {locus_sequence}"

    # Assert locus contains a single record
    locus_assembly = SeqIO.to_dict(SeqIO.parse(locus_sequence, "fasta"))
    assert len(locus_assembly) == 1, f"Locus file {locus_sequence} contains more than one record"

    locus_seq_record = next(iter(locus_assembly.values()))
    locus_record = SeqRecord(Seq(str(locus_seq_record.seq)),name=reference_id,id=reference_id,description=f"Ref_{reference_id}")
    query_bases = safe_int(header_data['Query_Assembly_Bases'][0])
    reference_bases = safe_int(header_data['Reference_Assembly_Bases'][0])

    query_contigs = safe_int(header_data['Query_Contig_Count'][0])
    reference_contigs = safe_int(header_data['Reference_Contig_Count'][0])

    raw_query_percent_aligned = safe_float(header_data['Query_Percent_Aligned'][0])
    raw_ref_percent_aligned = safe_float(header_data['Reference_Percent_Aligned'][0])
        
    # If the reference is not covered by at least min_cov, STOP
    if raw_ref_percent_aligned < min_cov:
        query_percent_aligned = raw_query_percent_aligned
        reference_percent_aligned = raw_ref_percent_aligned
        screen_category = "Low_Coverage"
        with open(log_file,"a+") as log:
            log.write(f"\t- Reference genome coverage: {raw_ref_percent_aligned}% \n")
            log.write(f"\t- Query covers less than --min_cov ({min_cov}%)...Locus finding halted...\n")
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
                
        except Exception as e:
            run_failed = True       
            with open(log_file,"a+") as log:
                log.write(f"\nError reading BED/SNP data from file: {snpdiffs_file}\n{str(e)}")
            sys.exit(f"Error reading BED/SNP data from file: {snpdiffs_file}\n{str(e)}")
            
        ##### 03: Filter genome overlaps #####
        with open(log_file,"a+") as log:
            log.write("Step 2: Filtering for short overlaps and low percent identity...")
        
        good_bed_df = bed_df[(bed_df['Ref_Aligned'] >= min_len) & (bed_df['Perc_Iden'] >= min_iden)].copy()
        
        if good_bed_df.shape[0] == 0:
            query_percent_aligned = raw_query_percent_aligned
            reference_percent_aligned = raw_ref_percent_aligned
            screen_category = "Low_Quality_Coverage"
            with open(log_file,"a+") as log:
                log.write(f"\n\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}) , no valid alignments remain...Locus finding halted...\n")
                log.write("-------------------------------------------------------\n\n")
   
        else:
            # Create a BED file for alignments that pass basic QC
            good_query_bed_df = good_bed_df[['Query_Contig','Query_Start','Query_End']].copy()
            good_reference_bed_df = good_bed_df[['Ref_Contig','Ref_Start','Ref_End']].copy()
            good_reference_bed_df.loc[:, 'Query_ID'] = query_id

            good_query_aligned = calculate_total_length(BedTool.from_dataframe(good_query_bed_df).sort().merge())
            good_reference_aligned = calculate_total_length(BedTool.from_dataframe(good_reference_bed_df[['Ref_Contig','Ref_Start','Ref_End']]).sort().merge())
            
            query_percent_aligned = (good_query_aligned / query_bases) * 100
            reference_percent_aligned = (good_reference_aligned / reference_bases) * 100
            
            if reference_percent_aligned < min_cov:
                screen_category = "Low_Quality_Coverage"
                with open(log_file,"a+") as log:
                    log.write(f"\n\t- Raw reference locus coverage was {raw_ref_percent_aligned}% \n")
                    log.write(f"\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}), reference locus coverage was {reference_percent_aligned:.2f}% \n")
                    log.write(f"\t- Query covers less than --min_cov ({min_cov}%) of reference locus after filtering...Locus finding halted...\n")
                    log.write("-------------------------------------------------------\n\n")

            else:
                with open(log_file,"a+") as log:
                    log.write("Done!\n")
                    log.write(f"\t- Raw reference locus coverage was {raw_ref_percent_aligned}% \n")
                    log.write(f"\t- After filtering based on --min_len ({min_len}) and --min_iden ({min_iden}), reference locus coverage was {reference_percent_aligned:.2f}% \n")
                    log.write("-------------------------------------------------------\n\n")

                with open(log_file,"a+") as log:
                    log.write("Step 3: Getting FASTA coordinates for locus...")
                
                locus_coords_df = getLocusCoords(bed_df,log_file, min_len, min_iden)
                                
                with open(log_file,"a+") as log:
                    log.write("Done!\n")
                    log.write("-------------------------------------------------------\n\n")
                
                # Check contig count for query
                if locus_coords_df['Query_Contig'].nunique() > max_contigs:
                    screen_category = "Fail_Max_Contigs"
                    with open(log_file,"a+") as log:
                        log.write(f"\t- More than {max_contigs} contigs required to cover the locus \n")
                        log.write("-------------------------------------------------------\n\n")
                else:
                    screen_category = "Pass"
                    with open(log_file,"a+") as log:
                        log.write("Step 3: FASTA generation...")

                    query_record = getFasta(locus_coords_df,query_assembly,query_id,reference_id,locus_record)
                    with open(f"{output_directory}/{query_id}.fasta", "w") as output_handle:
                        SeqIO.write(query_record, output_handle, "fasta")
                    with open(log_file,"a+") as log:
                        log.write("Done!\n")
                        log.write(f"Saved FASTA to {output_directory}/{query_id}.fasta")
                        log.write("-------------------------------------------------------\n\n")


    screen_end_time = time.time()
        
    with open(log_file,"a+") as log:
        log.write(f"Locus Search Time: {screen_end_time - screen_start_time:.2f} seconds\n")
    
    return [query_id,reference_percent_aligned,screen_category]
    
    # Clean up pybedtools temp
    helpers.cleanup(verbose=False, remove_all=False)
    
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
log_directory = os.path.abspath(args.log_directory)
log_file = f"{output_directory}/CSP2_SNP_Pipeline.log"
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
    # Establish output files
    reference_screening_file = f"{output_directory}/Reference_Screening.tsv"

    with open(log_file,"a+") as log:
        log.write("Searching all queries for locus...")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(findLoci,snp_diff_file,trim_name, min_cov, min_len, min_iden, max_contigs,reference_id,log_directory,output_directory) for snp_diff_file in snpdiffs_list]

except:
    run_failed = True
    print("Exception occurred:\n", traceback.format_exc())
finally:
    helpers.cleanup(verbose=False, remove_all=False)
    if temp_dir != "":
        shutil.rmtree(temp_dir)
    if run_failed:
        sys.exit(1)
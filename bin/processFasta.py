#!/usr/bin/env python3

from Bio import SeqIO
import hashlib
import sys
import os
import pandas as pd

def fasta_info(sample_name,file_path):
    
    if not os.path.exists(file_path):
        sys.exit(f"File {file_path} does not exist.")
    elif not file_path.lower().endswith(('.fa', '.fasta', '.fna')):
        sys.exit(f"File {file_path} is not a .fa, .fasta, or .fna file.")
    
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

    print(f"{sample_name},{file_path},{contig_count},{assembly_bases},{n50},{n90},{l50},{l90},{sha256}")

fasta_tsv = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['Sample_ID','Fasta_Path'])
[fasta_info(sample_name, file_path) for sample_name, file_path in fasta_tsv.values]
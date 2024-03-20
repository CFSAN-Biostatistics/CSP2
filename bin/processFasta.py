#!/usr/bin/env python3

from Bio import SeqIO
import hashlib
import sys
import os

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

    return int(contig_count), int(assembly_bases), int(n50),int(l50), int(n90),int(l90), sha256

sample_name = str(sys.argv[1])
fasta_file = str(sys.argv[2])

if not os.path.exists(fasta_file):
    sys.exit(f"File {fasta_file} does not exist.")
elif not fasta_file.lower().endswith(('.fa', '.fasta', '.fna')):
    sys.exit(f"File {fasta_file} is not a .fa, .fasta, or .fna file.")

contig_count, assembly_bases, n50, l50, n90, l90, sha256 = fasta_info(fasta_file)
print(f"{sample_name},{fasta_file},{contig_count},{assembly_bases},{n50},{n90},{l50},{l90},{sha256}")

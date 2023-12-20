#!/usr/bin/env python3

from Bio import SeqIO
import hashlib
import sys
import os

def fasta_info(file_path):
    records = list(SeqIO.parse(file_path, 'fasta'))
    contig_count = len(records)
    assembly_bases = sum(len(record) for record in records)
    with open(file_path, 'rb') as file:
        sha256 = hashlib.sha256(file.read()).hexdigest()
    return contig_count, assembly_bases, sha256

sample_name = str(sys.argv[1])
read_type = str(sys.argv[2])
read_location = str(sys.argv[3])
fasta_file = str(sys.argv[4])

if not os.path.exists(fasta_file):
    sys.exit(f"File {fasta_file} does not exist.")
elif not fasta_file.lower().endswith(('.fa', '.fasta', '.fna')):
    sys.exit(f"File {fasta_file} is not a .fa, .fasta, or .fna file.")

contig_count, assembly_bases, sha256 = fasta_info(fasta_file)


print(f"{sample_name},{read_type},{read_location},{fasta_file},{contig_count},{assembly_bases},{sha256}")

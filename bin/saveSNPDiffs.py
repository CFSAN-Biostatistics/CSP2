#!/usr/bin/env python3

import os
import sys
import pandas as pd

def processHeader(header_row):
    header_cols = [item.split(':')[0] for item in header_row]
    header_vals = [item.split(':')[1] for item in header_row]

    header_data = pd.DataFrame(columns = header_cols)
    header_data.loc[0] = header_vals
    return header_data

snpdiffs_list_file = sys.argv[1]
summary_file = sys.argv[2]

# Read in all lines and ensure each file exists
snpdiffs_list = [line.strip() for line in open(snpdiffs_list_file, 'r')]
for snpdiffs_file in snpdiffs_list:
    if not os.path.exists(snpdiffs_file):
        sys.exit("Error: File does not exist: " + snpdiffs_file)

header_rows = []
for snpdiffs_file in snpdiffs_list:
    with open(snpdiffs_file, 'r') as file:
        top_line = file.readline().lstrip('#').strip().split('\t')
        header_rows.append(processHeader(top_line))

pd.concat(header_rows, ignore_index=True).to_csv(summary_file, sep='\t', index=False)
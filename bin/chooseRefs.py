#!/usr/bin/env python3

import numpy as np
import os
import pandas as pd
import sys
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import scipy.stats
from itertools import combinations
from Bio import SeqIO
import argparse

def getOptimalK(data, ref_count):
    
    silhouette_scores = []
    
    kmeans_1 = KMeans(n_clusters=1, random_state=0, n_init='auto').fit(data)
    kmeans_2 = KMeans(n_clusters=2, random_state=0, n_init='auto').fit(data)
    
    # Compare 1 vs. 2
    inertia_1 = kmeans_1.inertia_
    inertia_2 = kmeans_2.inertia_
    if inertia_1 > inertia_2:
        negative_movements = 1
    else:
        negative_movements = 0
    
    # Add k2 data
    labels = kmeans_2.labels_
    score = silhouette_score(data, labels)
    silhouette_scores.append(score)
    prev_score = score
    
    for k in range(3, ref_count + 3):
        kmeans = KMeans(n_clusters=k, random_state=0, n_init='auto').fit(data)
        labels = kmeans.labels_
        score = silhouette_score(data, labels)

        if score < prev_score:
            negative_movements += 1
        else:
            negative_movements = 0
        
        silhouette_scores.append(score)
        
        # Stop if two consecutive negative movements occur
        if negative_movements == 2:
            break
        
        prev_score = score

    if (inertia_1 < inertia_2) & (silhouette_scores[0] > silhouette_scores[1]):
        optimal_k = 1
    else: 
        optimal_k = np.argmax(silhouette_scores) + 2

    return optimal_k

def fasta_info(file_path):
    records = list(SeqIO.parse(file_path, 'fasta'))
    contig_count = int(len(records))
    lengths = sorted([len(record) for record in records], reverse=True)
    assembly_bases = sum(lengths)

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

    return [file_path,contig_count,assembly_bases,n50,n90,l50,l90]

parser = argparse.ArgumentParser(description='Choose reference isolates based on FASTA metrics and mean distances.')
parser.add_argument('--ref_count', type=int, default=1, help='Number of reference isolates to select')
parser.add_argument('--mash_triangle_file', type=str, help='Path to the mash triangle file')
parser.add_argument('--trim_name', nargs='?', const="", default="", type=str, help='Trim name')
args = parser.parse_args()

ref_count = args.ref_count
mash_triangle_file = os.path.abspath(args.mash_triangle_file)
trim_name = args.trim_name
ref_file = os.path.join(os.path.dirname(mash_triangle_file), 'CSP2_Ref_Selection.tsv')

# Get Sample IDs
sample_df = pd.read_csv(mash_triangle_file, sep='\t', usecols=[0], skip_blank_lines=True).dropna()
sample_df = sample_df[sample_df[sample_df.columns[0]].str.strip() != '']
sample_df.columns = ['Path']
sample_df['Isolate_ID'] = [os.path.splitext(os.path.basename(file))[0].replace(trim_name, '') for file in sample_df[sample_df.columns[0]].tolist()]
assembly_names = [os.path.splitext(os.path.basename(file))[0].replace(trim_name, '') for file in sample_df[sample_df.columns[0]].tolist()]
num_isolates = sample_df.shape[0]

# Get FASTA metrics
metrics_df = pd.DataFrame(sample_df['Path'].apply(fasta_info).tolist(), columns=['Path', 'Contigs', 'Length', 'N50','N90','L50','L90'])
metrics_df['Assembly_Bases_Zscore'] =  metrics_df['Length'].transform(scipy.stats.zscore).astype('float').round(3).fillna(0)
metrics_df['Contig_Count_Zscore'] =  metrics_df['Contigs'].transform(scipy.stats.zscore).astype('float').round(3).fillna(0)
metrics_df['N50_Zscore'] =  metrics_df['N50'].transform(scipy.stats.zscore).astype('float').round(3).fillna(0)

# Find outliers
inlier_df = metrics_df.loc[(metrics_df['N50_Zscore'] > -3) & 
                        (metrics_df['Assembly_Bases_Zscore'] < 3) & 
                        (metrics_df['Assembly_Bases_Zscore'] > -3) & 
                        (metrics_df['Contig_Count_Zscore'] < 3)]

inlier_count = inlier_df.shape[0]
inlier_isolates = [os.path.splitext(os.path.basename(file))[0].replace(trim_name, '') for file in inlier_df[inlier_df.columns[0]].tolist()]

# If not enough or just enough inliers, script is done
if ref_count > inlier_count:
    sys.exit("Error: Fewer inliers than requested references?")
elif ref_count == inlier_count:
    print(",".join(inlier_df['Path'].tolist()))
    sys.exit(0)
    
# Left join metrics_df and inlier_df
sample_df = inlier_df.merge(sample_df, on = "Path", how='left')[['Isolate_ID','Path','Contigs','Length','N50','N90','L50','L90','N50_Zscore']]

# Create distance matrix
with open(mash_triangle_file) as mash_triangle:
    a = np.zeros((num_isolates, num_isolates))
    mash_triangle.readline()
    mash_triangle.readline()
    idx = 1
    for line in mash_triangle:
        tokens = line.split()
        distances = [float(token) for token in tokens[1:]]
        a[idx, 0: len(distances)] = distances
        a[0: len(distances), idx] = distances
        idx += 1

dist_df = pd.DataFrame(a, index=assembly_names, columns=assembly_names).loc[inlier_isolates,inlier_isolates]
# Get mean distances after masking diagonal
mask = ~np.eye(dist_df.shape[0], dtype=bool)
mean_distances = dist_df.where(mask).mean().reset_index()
mean_distances.columns = ['Isolate_ID', 'Mean_Distance']

sample_df = sample_df.merge(mean_distances, on='Isolate_ID', how='left')
sample_df['Mean_Distance_Zscore'] = sample_df['Mean_Distance'].transform(scipy.stats.zscore).astype('float').round(3)
sample_df['Base_Score'] = sample_df['N50_Zscore'] - sample_df['Mean_Distance_Zscore'].fillna(0)

if ref_count == 1:
    print(",".join(sample_df.nlargest(1, 'Base_Score')['Path'].tolist()))
    sys.exit(0)

optimal_k = getOptimalK(dist_df, ref_count)

if optimal_k == 1:
    print(",".join(sample_df.nlargest(ref_count, 'Base_Score')['Path'].tolist()))
    sys.exit(0)

kmeans = KMeans(n_clusters=optimal_k, random_state=0,n_init='auto').fit(dist_df)
clusters = kmeans.labels_

cluster_df = pd.DataFrame({'Isolate_ID': dist_df.index, 'Cluster': clusters}).merge(sample_df, on='Isolate_ID',how='left')
cluster_counts = cluster_df['Cluster'].value_counts().reset_index()
cluster_counts.columns = ['Cluster', 'count']
cluster_counts['Prop'] = cluster_counts['count'] / cluster_counts['count'].sum()
cluster_df = cluster_df.merge(cluster_counts[['Cluster', 'Prop']], on='Cluster')

# Grab top ref 
final_ref_df = cluster_df.nlargest(1, 'Base_Score')
refs_chosen = final_ref_df['Isolate_ID'].tolist()

possible_refs = cluster_df.loc[~cluster_df['Isolate_ID'].isin(refs_chosen)].copy()

while len(refs_chosen) < ref_count:
    possible_refs['Mean_Ref_Distance'] = possible_refs['Isolate_ID'].apply(lambda isolate_id: np.mean(dist_df.loc[isolate_id, refs_chosen].values))
    possible_refs['Mean_Ref_Distance_Zscore'] = possible_refs['Mean_Ref_Distance'].transform(scipy.stats.zscore).astype('float').round(3)      
    possible_refs['Sort_Score'] = possible_refs.apply(lambda row: (row['Base_Score'] + row['Mean_Ref_Distance_Zscore']) if row['Mean_Ref_Distance_Zscore'] <= 0 else (row['Base_Score'] + (row['Mean_Ref_Distance_Zscore']*row['Prop'])), axis=1)
    
    final_ref_df = pd.concat([final_ref_df, possible_refs.nlargest(1, 'Sort_Score').drop(['Sort_Score','Mean_Ref_Distance','Mean_Ref_Distance_Zscore'],axis=1)])
    refs_chosen = final_ref_df['Isolate_ID'].tolist()
    possible_refs = possible_refs.loc[~possible_refs['Isolate_ID'].isin(refs_chosen)].copy()

non_ref_df = cluster_df.loc[~cluster_df['Isolate_ID'].isin(refs_chosen)].sort_values('Base_Score', ascending=False)
non_ref_df['Is_Ref'] = False
final_ref_df['Is_Ref'] = True
pd.concat([final_ref_df, non_ref_df]).reset_index(drop=True).to_csv(ref_file, index=False, sep="\t")

print(",".join(final_ref_df['Path'].tolist()))
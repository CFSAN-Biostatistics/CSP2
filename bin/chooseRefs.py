#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import scipy.stats
from itertools import combinations


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

ref_count = int(sys.argv[1])
mean_distance_df = pd.read_csv(sys.argv[2], sep="\t").rename(columns = {'Assembly':'Isolate_ID'})
isolate_count = mean_distance_df.shape[0]

mean_distance_df['Assembly_Bases_Zscore'] =  mean_distance_df['Length'].transform(scipy.stats.zscore).astype('float').round(3)
mean_distance_df['Contig_Count_Zscore'] =  mean_distance_df['Contigs'].transform(scipy.stats.zscore).astype('float').round(3)
mean_distance_df['N50_Zscore'] =  mean_distance_df['N50'].transform(scipy.stats.zscore).astype('float').round(3)

inlier_df = mean_distance_df.loc[(mean_distance_df['N50_Zscore'] > -3) & 
                        (mean_distance_df['Assembly_Bases_Zscore'] < 3) & 
                        (mean_distance_df['Assembly_Bases_Zscore'] > -3) & 
                        (mean_distance_df['Contig_Count_Zscore'] < 3)]

inlier_count = inlier_df.shape[0]

if ref_count > inlier_count:
    sys.exit("Error: Fewer inliers than requested references?")
elif ref_count == inlier_count:
    print(",".join(inlier_df['Path'].tolist()))
    sys.exit(0)

inlier_isolates = inlier_df['Isolate_ID'].tolist()

distance_matrix = pd.read_csv(sys.argv[3], sep="\t", index_col=0)
pruned_distance_matrix = distance_matrix.loc[inlier_isolates,inlier_isolates]

pruned_mean_data = [(x, pruned_distance_matrix.loc[x][pruned_distance_matrix.columns != x].mean()) for x in inlier_isolates]
pruned_mean_df = pd.DataFrame(pruned_mean_data, columns=['Isolate_ID','Pruned_Mean_Distance'])

inlier_df = inlier_df.merge(pruned_mean_df, on='Isolate_ID')
inlier_df['Mean_Distance_Zscore'] = inlier_df['Pruned_Mean_Distance'].transform(scipy.stats.zscore).astype('float').round(3)
inlier_df['Base_Score'] = inlier_df['N50_Zscore'] - inlier_df['Mean_Distance_Zscore']

if ref_count == 1:
    print(",".join(inlier_df.nlargest(1, 'Base_Score')['Path'].tolist()))
    sys.exit(0)

optimal_k = getOptimalK(pruned_distance_matrix, ref_count)

if optimal_k == 1:
    print(",".join(inlier_df.nlargest(ref_count, 'Base_Score')['Path'].tolist()))
    sys.exit(0)

kmeans = KMeans(n_clusters=optimal_k, random_state=0,n_init='auto').fit(pruned_distance_matrix)
clusters = kmeans.labels_

cluster_df = pd.DataFrame({'Isolate_ID': pruned_distance_matrix.index, 'Cluster': clusters}).merge(inlier_df, on='Isolate_ID',how='left')

cluster_size_df = cluster_df['Cluster'].value_counts().reset_index().rename(columns={'index':'Cluster','Cluster':'count'})
cluster_size_df['Prop'] = cluster_size_df['count']/cluster_size_df['count'].sum()
cluster_df = cluster_df.merge(cluster_size_df[['Cluster','Prop']], on='Cluster')

# Grab top ref 
final_ref_df = cluster_df.nlargest(1, 'Base_Score')
refs_chosen = final_ref_df['Isolate_ID'].tolist()

possible_refs = cluster_df.loc[~cluster_df['Isolate_ID'].isin(refs_chosen)].copy()

while len(refs_chosen) < ref_count:
    possible_refs['Mean_Ref_Distance'] = possible_refs['Isolate_ID'].apply(lambda isolate_id: np.mean(pruned_distance_matrix.loc[isolate_id, refs_chosen].values))
    possible_refs['Mean_Ref_Distance_Zscore'] = possible_refs['Mean_Ref_Distance'].transform(scipy.stats.zscore).astype('float').round(3)      
    possible_refs['Sort_Score'] = possible_refs.apply(lambda row: (row['Base_Score'] + row['Mean_Ref_Distance_Zscore']) if row['Mean_Ref_Distance_Zscore'] <= 0 else (row['Base_Score'] + (row['Mean_Ref_Distance_Zscore']*row['Prop'])), axis=1)
    
    final_ref_df = pd.concat([final_ref_df, possible_refs.nlargest(1, 'Sort_Score').drop(['Sort_Score','Mean_Ref_Distance','Mean_Ref_Distance_Zscore'],axis=1)])
    refs_chosen = final_ref_df['Isolate_ID'].tolist()
    possible_refs = possible_refs.loc[~possible_refs['Isolate_ID'].isin(refs_chosen)].copy()

print(cluster_df.sort_values(by='Base_Score',ascending=False).head(20))
print(final_ref_df)
print(",".join(final_ref_df['Path'].tolist()))

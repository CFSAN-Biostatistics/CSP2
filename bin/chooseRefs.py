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
    prev_score = None
    
    kmeans_1 = KMeans(n_clusters=1, random_state=0, n_init='auto').fit(data)
    kmeans_2 = KMeans(n_clusters=2, random_state=0, n_init='auto').fit(data)
    
    inertia_1 = kmeans_1.inertia_
    inertia_2 = kmeans_2.inertia_

    if inertia_1 > inertia_2:
        negative_movements = 1
    else:
        negative_movements = 0
        
    for k in range(2, ref_count + 2):
        kmeans = KMeans(n_clusters=k, random_state=0, n_init='auto').fit(data)
        labels = kmeans.labels_
        score = silhouette_score(data, labels)
        silhouette_scores.append(score)
        
        if prev_score is None:
            pass
        elif score < prev_score:
            negative_movements += 1
        else:
            negative_movements = 0
        
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
mean_distance_df = pd.read_csv(sys.argv[2], sep="\t")
isolate_count = mean_distance_df.shape[0]

if ref_count >= isolate_count:
    print(",".join(mean_distance_df['Path'].tolist()))
    sys.exit(0)

mean_distance_df['Assembly_Bases_Zscore'] =  mean_distance_df['Length'].transform(scipy.stats.zscore).astype('float').round(3)
mean_distance_df['Contig_Count_Zscore'] =  mean_distance_df['Contigs'].transform(scipy.stats.zscore).astype('float').round(3)
mean_distance_df['N50_Zscore'] =  mean_distance_df['N50'].transform(scipy.stats.zscore).astype('float').round(3)
mean_distance_df['N90_Zscore'] =  mean_distance_df['N90'].transform(scipy.stats.zscore).astype('float').round(3)

inlier_df = mean_distance_df.loc[(mean_distance_df['N50_Zscore'] > -3) & 
                        (mean_distance_df['N90_Zscore'] > -3) & 
                        (mean_distance_df['Assembly_Bases_Zscore'] < 3) & 
                        (mean_distance_df['Assembly_Bases_Zscore'] > -3) & 
                        (mean_distance_df['Contig_Count_Zscore'] < 3)]

inlier_count = inlier_df.shape[0]
inlier_isolates = inlier_df['Assembly'].tolist()
outlier_count = mean_distance_df.shape[0] - inlier_count

if ref_count > inlier_count:
    sys.exit("Error: Fewer inliers than requested references?")

distance_matrix = pd.read_csv(sys.argv[3], sep="\t", index_col=0)
pruned_distance_matrix = distance_matrix.loc[inlier_isolates,inlier_isolates]
optimal_k = getOptimalK(pruned_distance_matrix, ref_count)

if optimal_k == 1:
    cluster_df = pd.DataFrame({'Isolate_ID': pruned_distance_matrix.index, 'Cluster': [0]*inlier_count}).merge(mean_distance_df, left_on='Isolate_ID', right_on='Assembly',how='left')
    cluster_df['Prop'] = 1
else:
    kmeans = KMeans(n_clusters=optimal_k, random_state=0,n_init='auto').fit(pruned_distance_matrix)
    clusters = kmeans.labels_

    cluster_df = pd.DataFrame({'Isolate_ID': pruned_distance_matrix.index, 'Cluster': clusters}).merge(mean_distance_df, left_on='Isolate_ID', right_on='Assembly',how='left')
    cluster_size_df = cluster_df['Cluster'].value_counts().reset_index()
    cluster_size_df['Prop'] = cluster_size_df['count']/cluster_size_df['count'].sum()
    cluster_df = pd.merge(cluster_df, cluster_size_df[['Cluster','Prop']], on='Cluster')

pairwise_combinations = [sorted(x) for x in list(combinations(inlier_isolates, 2))]
pairwise_distances = [pruned_distance_matrix.loc[x,y] for x,y in pairwise_combinations]

pairwise_df = pd.DataFrame(pairwise_combinations, columns=['Isolate_1','Isolate_2'])
pairwise_df['Distance'] = pairwise_distances

pairwise_df = pd.merge(pairwise_df, cluster_df[['Isolate_ID','Cluster']], left_on='Isolate_1', right_on='Isolate_ID').rename(columns = {'Cluster':'Cluster_1'}).drop('Isolate_ID',axis=1)
pairwise_df = pd.merge(pairwise_df, cluster_df[['Isolate_ID','Cluster']], left_on='Isolate_2', right_on='Isolate_ID').rename(columns = {'Cluster':'Cluster_2'}).drop('Isolate_ID',axis=1)

pairwise_df['Dist_Type'] = np.where(pairwise_df['Cluster_1'] == pairwise_df['Cluster_2'], 'Within', 'Across')
    
isolate_mean_data = []
for isolate in inlier_isolates:
    mean_df = pairwise_df.loc[(pairwise_df['Isolate_1'] == isolate) | (pairwise_df['Isolate_2'] == isolate)].groupby('Dist_Type')['Distance'].agg(Mean_Distance = 'mean')
    if optimal_k > 1:
        isolate_mean_data.append((isolate, mean_df.loc['Within'].values[0], mean_df.loc['Across'].values[0]))
    else:
        isolate_mean_data.append((isolate, mean_df.loc['Within'].values[0]))

if optimal_k > 1:
    isolate_mean_df = pd.DataFrame(isolate_mean_data, columns=['Isolate_ID','Within_Mean','Across_Mean'])
else:
    isolate_mean_df = pd.DataFrame(isolate_mean_data, columns=['Isolate_ID','Within_Mean'])
    isolate_mean_df['Across_Mean'] = 1
    
cluster_df = pd.merge(cluster_df,isolate_mean_df,on='Isolate_ID',how='left')
cluster_df['New_Score'] = (cluster_df['Prop']*cluster_df['N50'])/(cluster_df['Within_Mean'] * cluster_df['Across_Mean'])
cluster_df.sort_values(by='New_Score',ascending=False,inplace=True)

if ref_count == 1:
    print(cluster_df['Path'].tolist()[0])
    sys.exit(0)    
elif optimal_k == 1:
    print(",".join(cluster_df.head(ref_count)['Path'].tolist()))
    sys.exit(0)

# Grab top ref 
final_ref_df = cluster_df.head(1)
refs_chosen = final_ref_df['Assembly'].tolist()

possible_refs = cluster_df.loc[~cluster_df['Isolate_ID'].isin(refs_chosen)].copy()

while len(refs_chosen) < ref_count:
    possible_refs['Mean_Ref_Distance'] = possible_refs['Isolate_ID'].apply(lambda isolate_id: np.mean(pruned_distance_matrix.loc[isolate_id, refs_chosen].values))
    possible_refs['Min_Ref_Distance'] = possible_refs['Isolate_ID'].apply(lambda isolate_id: np.min(pruned_distance_matrix.loc[isolate_id, refs_chosen].values))
    possible_refs['Sort_Score'] = possible_refs['New_Score'] * possible_refs['Mean_Ref_Distance'] * possible_refs['Min_Ref_Distance']
    
    final_ref_df = pd.concat([final_ref_df, possible_refs.nlargest(1, 'Sort_Score').drop('Sort_Score',axis=1)])
    
    refs_chosen = final_ref_df['Assembly'].tolist()
    possible_refs = possible_refs.loc[~possible_refs['Isolate_ID'].isin(refs_chosen)].copy()

print(",".join(final_ref_df['Path'].tolist()))
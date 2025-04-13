"""
File: clustering_plot.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering


def plot_clustering(
    df_struc, df_pvals, 
    df_pvals_clust, dist, 
    workdir, input_gene, 
    distances, yvalues, 
    psig_columns=[f'SUM_LFC3D_neg_05_psig', f'SUM_LFC3D_pos_05_psig'], 
    names=['Negative', 'Positive'], 
    pthr_cutoffs=['p<0.05', 'p<0.05'], 
    screen_name = 'Meta', score_type='LFC3D',  
    merge_col=['unipos', 'chain'], 
    clustering_kwargs = {"n_clusters": None, "metric": "euclidean", "linkage": "single"}, 

    horizontal=False, 
    line_subplots_kwargs={'figsize':(10,7)}, 
    dendogram_subplots_kwargs={'figsize':(15, 12)}, 
    save_type='png', 
): 
    """
    Description
        Calculates number of clusters for one clustering radius
    """
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / f'cluster_{score_type}'):
        os.mkdir(working_filedir / f'cluster_{score_type}')
    if not os.path.exists(working_filedir / f'cluster_{score_type}/plots'):
        os.mkdir(working_filedir / f'cluster_{score_type}/plots')

    # ASSERT df_struc HAS CHAIN AND POSITIONAL INFO #
    coord_columns = ["x_coord", "y_coord", "z_coord"]
    structure_columns = ["unipos", "unires", "chain"] + coord_columns
    for column in structure_columns: 
        assert column in df_struc.columns

    # ASSERT df_pvals HAS CATEGORICAL VARIABLES INDICATING WHAT TO CLUSTER #
    for column in ["unipos", "unires", "chain"] + psig_columns: 
        assert column in df_pvals.columns

    # CHECK INPUTS ARE SELF CONSISTENT #
    assert len(df_struc) == len(df_pvals)
    assert len(psig_columns) == len(names) == len(pthr_cutoffs)
    
    # SETUP DF #
    df_hits_clust = pd.concat([df_struc[structure_columns], df_pvals[psig_columns]], axis=1)
    prefix = f'{input_gene}_{screen_name}_{score_type}'
    pos_col, chain_col = merge_col[0], merge_col[1]

    # PLOT CLUSTERING DIST VS NUM OF CLUSTERS #
    clust_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_Aggr_Hits_List.tsv" 
    plot_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_cluster_distance.{save_type}"
    plot_cluster_distance(distances, yvalues, 
                          names, input_gene, 
                          clust_filename, plot_filename, 
                          line_subplots_kwargs, save_type)

    # OPEN CLUSTERING FILE #
    for name, pthr, colname in zip(names, pthr_cutoffs, psig_columns): 
        # EXTRACT ROWS ABOVE CUTOFF #
        df_pvals_temp = df_hits_clust.loc[(df_hits_clust[colname] == pthr), ].reset_index(drop=True)
        # REMOVE ROWS WITHOUT POSITION INFO FOR PDBs #
        df_pvals_temp = df_pvals_temp[~df_pvals_temp[coord_columns].isin(['-']).any(axis=1)]

        np_hits_coord = np.array(df_pvals_temp[coord_columns]).copy()
        if np_hits_coord.shape[0] < 2: 
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            raise None
        
        # RUN CLUSTERING #
        func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=dist)
        clustering = func_clustering.fit(np_hits_coord)

        dend_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_{name}_Dendogram_{str(dist)}A.{save_type}"
        title = f'{input_gene} {score_type} {name} Clusters'
        plot_dendrogram(clustering, df_pvals_temp, 
                        dist, horizontal, pos_col, chain_col, 
                        title, dend_filename, 
                        dendogram_subplots_kwargs, save_type)

        # CLUSTERS RESIDUES AND LENGTH OF EACH CLUSTER #
        df_pvals_clust_i = df_pvals_clust.loc[(df_pvals_clust[colname] == pthr), ].reset_index(drop=True)
        clust_indices = df_pvals_clust_i[f'{colname}_Clust_{str(dist)}A'].unique()

        txt_filename = working_filedir / f"cluster_{score_type}/{prefix}_{name}_Dendrogram_{str(dist)}A.txt"
        with open(txt_filename, "w") as f: 
            for c in clust_indices: 
                c_data = df_pvals_clust_i.loc[df_pvals_clust_i[f'{colname}_Clust_{str(dist)}A'] == c, ].reset_index(drop=True)

                # WRITE LENGTH, RANGE, AND ALL RESIDUES #
                if len(c_data) > 0: 
                    all_unipos = c_data[pos_col].tolist()
                    all_chains = c_data[chain_col].tolist()
                    f.write(f'Cluster {c} : Length {len(c_data)} :\n')
                    f.write('   ')
                    for unipos, chain in zip(all_unipos, all_chains): 
                        f.write(f'{chain}-{unipos} ')
                    f.write(f'\n')

    return clust_indices


def plot_cluster_distance(
        distances, yvalues, 
        names, input_gene, 
        clust_filename, plot_filename, 
        subplots_kwargs, save_type, 
): 

    dist_dict = {'clust_dist': distances}
    for n, y in zip(names, yvalues): 
        dist_dict[n] = y
    dist_stat = pd.DataFrame(dist_dict)
    dist_stat.to_csv(clust_filename, sep='\t', index=False)

    fig, ax = plt.subplots(**subplots_kwargs)
    for n in names: 
        sns.lineplot(data=dist_stat, x="clust_dist", y=n)

    plt.xlabel('Cluster Radius')
    plt.ylabel('Number of Clusters')
    plt.title(f'Positive vs Negative Clusters {input_gene}')
    plt.savefig(plot_filename, dpi=100, transparent=True, format=save_type)
    plt.close()

def plot_dendrogram(
        clustering, df_pvals_temp, 
        dist, horizontal, pos_col, chain_col, 
        title, dend_filename, 
        subplots_kwargs, save_type, 
):  
    fig, ax = plt.subplots(**subplots_kwargs)
    counts = np.zeros(clustering.children_.shape[0]) # CREATE COUNTS OF SAMPLE #
    n_samples = len(clustering.labels_)
    
    for i, merge in enumerate(clustering.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [clustering.children_, clustering.distances_, counts]).astype(float)
    xlbl_pos = list(df_pvals_temp[pos_col])
    xlbl_chain = list(df_pvals_temp[chain_col])
    xlbl = [f'{pos}-{chain}' for pos, chain in zip(xlbl_pos, xlbl_chain)]

    # PLOT CORRESPONDING DENDROGRAM #
    if horizontal: dendrogram(linkage_matrix, color_threshold=dist, labels=xlbl, orientation='right')
    else: dendrogram(linkage_matrix, color_threshold=dist, labels=xlbl, leaf_rotation=90.)

    plt.title(title)
    plt.savefig(dend_filename, dpi=100, transparent=True, format=save_type)
    plt.close()

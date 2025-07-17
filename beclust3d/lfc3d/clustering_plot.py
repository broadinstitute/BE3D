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
    dendrogram_subplots_kwargs={'figsize':(15, 12)}, 
    save_type='png', 
): 
    """
    Calculates number of clusters for one clustering radius and its associated plots.

    Parameters
    ----------
    df_str_cons : pd.DataFrame
        DataFrame containing structural data for residues. 
        Must include ['unipos', 'unires', 'chain', 'x_coord', 'y_coord', 'z_coord'].

    df_pvals : pd.DataFrame
        DataFrame containing per-residue statistical significance categories.
        Must include ['unipos', 'unires', 'chain'] plus columns listed in `psig_columns`.

    df_pvals_clust : pd.DataFrame
        DataFrame containing structure and significance information plus cluster labels assigned at each distance.

    dist : int, optional (default=25)
        Tadius (in Angstroms) to consider for clustering. 

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    distances : list of int
        List of distances (from 1 to `max_distances`) at which clustering was performed.

    yvalue_lists : list of list of int
        List of lists, containing the number of clusters found at each distance for each psig_column.

    psig_columns : list of str, optional
        List of column names in `df_pvals` indicating categorical significance labels 
        (e.g., 'p<0.05' for significant residues to cluster).

    names : list of str, optional
        List of names corresponding to `psig_columns`.

    pthr_cutoffs : list of str, optional
        List of significance thresholds corresponding to `psig_columns`.
        Only residues matching the given thresholds are included in clustering.

    screen_name : str
        Name of the screens corresponding to df_missense.

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed (e.g., 'LFC3D', 'LFC', etc.).

    max_distances : int, optional (default=25)
        Maximum radius (in Angstroms) to consider for clustering. Clustering is repeated at every integer from 1 to `max_distances`.

    merge_cols : list of str, optional (default=['unipos', 'chain'])
        Columns used to merge clustering results back into the main DataFrame.

    clustering_kwargs : dict, optional
        Dictionary of additional keyword arguments passed to `AgglomerativeClustering`.
        Must include keys like "metric" and "linkage".
        "n_clusters" should be set to None to enable distance-threshold clustering.

    Returns
    -------
    None
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
            continue
        
        # RUN CLUSTERING #
        func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=dist)
        clustering = func_clustering.fit(np_hits_coord)

        dend_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_{name}_Dendrogram_{str(int(dist))}A.{save_type}"
        title = f'{input_gene} {score_type} {name} Clusters'
        plot_dendrogram(clustering, df_pvals_temp, 
                        dist, horizontal, pos_col, chain_col, 
                        title, dend_filename, 
                        dendrogram_subplots_kwargs, save_type)

        # CLUSTERS RESIDUES AND LENGTH OF EACH CLUSTER #
        df_pvals_clust_i = df_pvals_clust.loc[(df_pvals_clust[colname] == pthr), ].reset_index(drop=True)
        clust_indices = df_pvals_clust_i[f'{colname}_Clust_{str(int(dist))}A'].unique()

        txt_filename = working_filedir / f"cluster_{score_type}/{prefix}_{name}_Dendrogram_{str(int(dist))}A.txt"
        with open(txt_filename, "w") as f: 
            for c in clust_indices: 
                c_data = df_pvals_clust_i.loc[df_pvals_clust_i[f'{colname}_Clust_{str(int(dist))}A'] == c, ].reset_index(drop=True)

                # WRITE LENGTH, RANGE, AND ALL RESIDUES #
                if len(c_data) > 0: 
                    all_unipos = c_data[pos_col].tolist()
                    all_chains = c_data[chain_col].tolist()
                    f.write(f'Cluster {c} : Length {len(c_data)} :\n')
                    f.write('   ')
                    for unipos, chain in zip(all_unipos, all_chains): 
                        f.write(f'{chain}-{unipos} ')
                    f.write(f'\n')

    return None


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
        sns.lineplot(data=dist_stat, x="clust_dist", y=n, ax=ax, label=n)

    ax.legend(title='')
    plt.xlabel('Cluster Radius')
    plt.ylabel('Number of Clusters')
    plt.title(f'Positive vs Negative Clusters {input_gene}')
    plt.savefig(plot_filename, dpi=100, transparent=False, format=save_type)
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
    plt.savefig(dend_filename, dpi=100, transparent=False, format=save_type)
    plt.close()

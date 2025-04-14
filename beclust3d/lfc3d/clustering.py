"""
File: clustering.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.cluster import AgglomerativeClustering


def clustering(
        df_struc, df_pvals, 
        workdir, input_gene, 
        psig_columns=[f'SUM_LFC3D_neg_05_psig', f'SUM_LFC3D_pos_05_psig'], # CATEGORICAL NOT QUANTITATIVE #
        pthr_cutoffs=['p<0.05', 'p<0.05'], 
        screen_name='Meta', score_type='LFC3D', 
        max_distances=25, merge_cols=['unipos', 'chain'], 
        clustering_kwargs = {"n_clusters": None, "metric": "euclidean", "linkage": "single"}, 
): 
    """
    Description
        Calculates number of clusters for a range of clustering radii
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / f'cluster_{score_type}'):
        os.mkdir(working_filedir / f'cluster_{score_type}')

    # ASSERT df_struc HAS CHAIN AND POSITIONAL INFO #
    coord_columns = ["x_coord", "y_coord", "z_coord"]
    structure_columns = ["unipos", "unires", "chain"] + coord_columns
    for column in structure_columns: 
        assert column in df_struc.columns

    # ASSERT df_pvals HAS CATEGORICAL VARIABLES INDICATING WHAT TO CLUSTER #
    for column in ["unipos", "unires", "chain"] + psig_columns: 
        assert column in df_pvals.columns, df_pvals.columns

    # CHECK INPUTS ARE SELF CONSISTENT #
    assert len(df_struc) == len(df_pvals)
    assert len(psig_columns) == len(pthr_cutoffs)
    
    # SETUP DF #
    df_hits_clust = pd.concat([df_struc[structure_columns], df_pvals[psig_columns]], axis=1)

    # CLUSTERING #
    distances = [int(i+1) for i in range(max_distances)] # CLUSTERING DISTANCE HYPERPARAM
    yvalue_lists = [[] for _ in psig_columns]
    
    for column, pthr, y_arr in zip(psig_columns, pthr_cutoffs, yvalue_lists): 
        # EXTRACT ROWS ABOVE CUTOFF #
        dict_hits = {}
        df_pvals_temp = df_hits_clust.loc[(df_hits_clust[column] == pthr), ].reset_index(drop=True)

        # REMOVE ROWS WITHOUT POSITION INFO FOR PDBs #
        df_pvals_temp = df_pvals_temp[~df_pvals_temp[coord_columns].isin(['-']).any(axis=1)]
        # POSITIONS AND CHAINS OF CLUSTERING CATEGORY #
        for col in merge_cols: 
            dict_hits[col] = list(df_pvals_temp[col])

        # EXTRACT X Y Z OF HITS ABOVE CUTOFF #
        np_hits_coord = np.array(df_pvals_temp[coord_columns].copy())
        if np_hits_coord.shape[0] < 2: # NO DATA TO CLUSTER ON #
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            y_arr.extend[[0 for _ in distances]]
            continue

        # FOR RANGE OF RADIUS, RUN CLUSTERING #
        for dist in distances: 
            func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=dist)
            clus_lbl = func_clustering.fit(np_hits_coord).labels_

            num_clusters = int(max(clus_lbl)+1) 
            # print(f'Number of clusters for {name} hits: d={dist} {num_clusters}')
            y_arr.append(num_clusters)

            dict_hits[f"{column}_Clust_{str(dist)}A"] = clus_lbl

        # CONSTRUCT A WHOLE DATAFRAME OF CLUSTERS FOR EVERY RESIDUE #
        df_hits_clust = df_hits_clust.merge(pd.DataFrame(dict_hits), how='left', on=merge_cols)

    df_hits_clust.fillna('-')

    # SAVE FILE #
    hits_filename = working_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_Aggr_Hits.tsv"
    df_hits_clust.to_csv(hits_filename, sep='\t', index=False)

    return df_hits_clust, distances, yvalue_lists

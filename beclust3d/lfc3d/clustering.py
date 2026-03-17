"""
File: clustering.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
    Performs spatial agglomerative clustering of significant residues over a range of distance thresholds.
"""

import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import cKDTree

def clustering(
    df_struc, 
    df_pvals, 
    workdir, 
    input_gene, 
    max_distances=25, 
    psig_columns=[f'SUM_LFC3D_neg_05_psig', f'SUM_LFC3D_pos_05_psig'], # CATEGORICAL NOT QUANTITATIVE #
    pthr_cutoffs=['p<0.05', 'p<0.05'], 
    screen_name='Meta', 
    score_type='LFC3D', 
    merge_cols=['unipos', 'chain'], 
    clustering_kwargs={"n_clusters": None, "metric": "euclidean", "linkage": "single"},
    atom_level=False,
): 
    """
    Performs spatial agglomerative clustering of significant residues over a range of distance thresholds.

    Parameters
    ----------
    df_struc : pd.DataFrame
        Structural feature DataFrame from sequence_structural_features(). 
        Must include columns ['unipos', 'unires', 'chain', 'x_coord', 'y_coord', 'z_coord'].

    df_pvals : pd.DataFrame
        Z-score and significance label DataFrame from znorm_score(), znorm_meta(), or prioritize_by_sequence(). 
        Must include columns ['unipos', 'unires', 'chain'] plus columns listed in `psig_columns`.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed (e.g., 'DNMT3A', 'MEN1'). 

    max_distances : int, optional (default=25)
        Maximum radius (in Angstroms) to consider for clustering. Clustering is repeated at every integer from 3 to `max_distances`.
        
    psig_columns : list of str, optional
        Column names in df_pvals containing categorical significance labels to cluster on
        (e.g., ['SUM_LFC3D_neg_05_psig', 'SUM_LFC3D_pos_05_psig']).

    pthr_cutoffs : list of str, optional
        Significance threshold values corresponding to each column in psig_columns.
        Only residues matching these values are included in clustering (e.g., ['p<0.05', 'p<0.05']).

    screen_name : str
        Screen identifier used in output filenames.

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed. Either 'LFC' or 'LFC3D'.

    merge_cols : list of str, optional (default=['unipos', 'chain'])
        Columns used to merge clustering results back into the main DataFrame.

    clustering_kwargs : dict, optional
        Keyword arguments passed to AgglomerativeClustering. 'n_clusters' should be None to enable distance-threshold clustering.
        Default: {"n_clusters": None, "metric": "euclidean", "linkage": "single"}.
        
    atom_level : bool, optional (default=False)
        If True, performs clustering at the atom level rather than residue level.

    Returns
    -------
    df_hits_clust : pd.DataFrame
        DataFrame containing structural and significance data plus cluster labels assigned at each distance threshold.

    distances : list of int
        List of distances from 1 to max_distances at which clustering was performed.

    yvalue_lists : list of list of int
        Number of clusters found at each distance for each column in psig_columns.
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
    distances = [int(i+1) for i in range(3,max_distances)] # CLUSTERING DISTANCE HYPERPARAM
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

        chain_unipos_list = list(zip(dict_hits['chain'], dict_hits['unipos'])) # this is only for the atom-level agglomerative clustering
        # EXTRACT X Y Z OF HITS ABOVE CUTOFF #
        np_hits_coord = np.array(df_pvals_temp[coord_columns].copy())
        if np_hits_coord.shape[0] < 2: # NO DATA TO CLUSTER ON #
            warnings.warn(f"Not enough data to perform agglomerative clustering")
            y_arr.extend([0 for _ in distances])
        # FOR RANGE OF RADIUS, RUN CLUSTERING #
        for dist in distances: 
            if np_hits_coord.shape[0] < 2: # NO DATA TO CLUSTER ON #
                dict_hits[f"{column}_Clust_{str(dist)}A"] = None
            else:
                clus_lbl = list()
                if atom_level:
                    clus_lbl = cluster_residues_from_pdb(pdb_file=atom_level, residue_list=chain_unipos_list, linkage='average')
                else:
                    func_clustering = AgglomerativeClustering(**clustering_kwargs, distance_threshold=dist)
                    clus_lbl = func_clustering.fit(np_hits_coord).labels_

                num_clusters = int(max(clus_lbl)+1) 
                y_arr.append(num_clusters)
                dict_hits[f"{column}_Clust_{str(dist)}A"] = clus_lbl

        # CONSTRUCT A WHOLE DATAFRAME OF CLUSTERS FOR EVERY RESIDUE #
        df_hits_clust = df_hits_clust.merge(pd.DataFrame(dict_hits), how='left', on=merge_cols)

    df_hits_clust.fillna('-')

    # SAVE FILE #
    hits_filename = working_filedir / f"cluster_{score_type}/{input_gene}_{screen_name}_Aggr_Hits.tsv"
    df_hits_clust.to_csv(hits_filename, sep='\t', index=False)

    return df_hits_clust, distances, yvalue_lists

def cluster_residues_from_pdb(pdb_file, residue_list, distance_threshold=6, linkage="single"):
    """
    Performs agglomerative clustering on residues (based on all-atom coordinates).

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.
    residue_list : list of tuples
        Each tuple is (chain_id, residue_number). Example: [("A", 42), ("A", 45)].
    distance_threshold : float, optional
        Distance threshold (Å) for clustering. Default = 3.5.
    linkage : str, optional
        Linkage method for clustering ('single', 'average', 'complete', etc.). Default = 'single'.

    Returns
    -------
    labels : np.ndarray
        Cluster labels for each residue in the provided list.
    """
    # --- Load structure ---
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("S", pdb_file)
    model = structure[0]

    # --- Collect heavy-atom coordinates per residue ---
    residues_coords = []
    for chain_id, resnum in residue_list:
        res = model[chain_id][resnum]
        if not is_aa(res, standard=True):
            continue
        coords = [atom.get_coord() for atom in res.get_atoms() if atom.element != "H"]
        residues_coords.append(np.array(coords))

    # --- Compute min interatomic distance matrix ---
    n = len(residues_coords)
    D = np.zeros((n, n))
    trees = [cKDTree(r) for r in residues_coords]
    for i in range(n):
        for j in range(i + 1, n):
            if len(residues_coords[i]) <= len(residues_coords[j]):
                dmin = trees[j].query(residues_coords[i], k=1)[0].min()
            else:
                dmin = trees[i].query(residues_coords[j], k=1)[0].min()
            D[i, j] = D[j, i] = dmin

    # --- Run clustering ---
    model = AgglomerativeClustering(
        metric="precomputed",
        linkage=linkage,
        distance_threshold=distance_threshold,
        n_clusters=None
    ).fit(D)

    return model.labels_

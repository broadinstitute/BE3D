"""
File: clustering_plot.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
    Generates line plots and dendrograms for clustering results at a specified distance threshold.
"""

import os
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl

from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from sklearn.cluster import AgglomerativeClustering
mpl.rcParams['svg.fonttype'] = 'none'

def plot_clustering(
    df_struc, 
    df_pvals, 
    df_pvals_clust, 
    dist, 
    workdir, 
    input_gene, 
    distances, 
    yvalues, 
    psig_columns=[f'SUM_LFC3D_neg_05_psig', f'SUM_LFC3D_pos_05_psig'], 
    names=['Negative', 'Positive'], 
    pthr_cutoffs=['p<0.05', 'p<0.05'], 
    screen_name = 'Meta', 
    score_type='LFC3D',  
    merge_col=['unipos', 'chain'], 
    clustering_kwargs = {"n_clusters": None, "metric": "euclidean", "linkage": "single"}, 
    horizontal=False, 
    line_subplots_kwargs={'figsize':(10,7)}, 
    dendrogram_subplots_kwargs={'figsize':(15, 12)}, 
    save_type='png', 
): 
    """
    Generates line plots and dendrograms for clustering results at a specified distance threshold.

    Parameters
    ----------
    df_struc : pd.DataFrame
        Structural feature DataFrame from sequence_structural_features(). 
        Must include columns ['unipos', 'unires', 'chain', 'x_coord', 'y_coord', 'z_coord'].

    df_pvals : pd.DataFrame
        Z-score and significance label DataFrame from znorm_score(), znorm_meta(), or prioritize_by_sequence(). 
        Must include columns ['unipos', 'unires', 'chain'] plus columns listed in `psig_columns`.

    df_pvals_clust : pd.DataFrame
        Cluster label DataFrame from clustering().

    dist : int, optional (default=25)
        Clustering radius in Angstroms to plot results for.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed (e.g., 'DNMT3A', 'MEN1'). 

    distances : list of int
        List of distances at which clustering was performed; output from clustering().

    yvalue_lists : list of list of int
        Number of clusters at each distance for each psig_column; output from clustering().

    psig_columns : list of str, optional
        Column names in df_pvals containing categorical significance labels to cluster on.
        Should match psig_columns used in clustering().

    names : list of str, optional
        Display names corresponding to each column in psig_columns
        (e.g., ['Negative', 'Positive']).

    pthr_cutoffs : list of str, optional
        Significance threshold values corresponding to each column in psig_columns.
        Should match pthr_cutoffs used in clustering().

    screen_name : str
        Screen identifier used in output filenames.

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed. Either 'LFC' or 'LFC3D'.

    merge_cols : list of str, optional (default=['unipos', 'chain'])
        Columns used to merge clustering results back into the main DataFrame.

    clustering_kwargs : dict, optional
        Dictionary of additional keyword arguments passed to `AgglomerativeClustering`.
        Must include keys like "metric" and "linkage". "n_clusters" should be set to None to enable distance-threshold clustering.
        
    horizontal : bool, optional (default=False)
        If True, renders plots with a horizontal layout.
        
    line_subplots_kwargs : dict, optional (default={'figsize': (10, 7)})
        Keyword arguments passed to the line plot figure.
        
    dendrogram_subplots_kwargs : dict, optional (default={'figsize': (15, 12)})
        Keyword arguments passed to the dendrogram figure.
        
    save_type : str, optional (default='png')
        File format for saved plots (e.g., 'png', 'pdf', 'svg').
        
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
    
    # GENERATE CONSISTENT COLORS FOR ALL CLUSTERS
    cluster_colors = generate_cluster_colors(100)
    
    # Save global color mapping to file for JavaScript reference
    color_map_file = working_filedir / f"cluster_{score_type}/cluster_color_mapping.json"
    with open(color_map_file, 'w') as f:
        json.dump({f"cluster_{i}": color for i, color in enumerate(cluster_colors)}, f, indent=2)

    # SETUP DF #
    df_hits_clust = pd.concat([df_struc[structure_columns], df_pvals[psig_columns]], axis=1)
    prefix = f'{input_gene}_{screen_name}_{score_type}'
    pos_col, chain_col = merge_col[0], merge_col[1]

    cluster_dist_list = list()
    for name, pthr, yvalue in zip(names, pthr_cutoffs, yvalues):
        cluster_dist_list.append([name,pthr,yvalue])
    cluster_dist_pd = pd.DataFrame(cluster_dist_list,columns=['name','pthr','yvalue'])
    
    for gid, gcont in cluster_dist_pd.groupby('pthr'):    
        clust_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_{gid}_Aggr_Hits_List.tsv" 
        plot_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_{gid}_cluster_distance.{save_type}"

        _names = gcont['name'].to_list()
        _yvalues = gcont['yvalue'].to_list()
        plot_cluster_distance(
            distances, _yvalues,
            _names, input_gene, 
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

        dend_filename = working_filedir / f"cluster_{score_type}/plots/{prefix}_{name}_Dendrogram_{pthr}_{str(int(dist))}A.{save_type}"
        title = f'{input_gene} {score_type} {name} Clusters'
        
        # MODIFIED: plot_dendrogram now returns detailed cluster info
        cluster_info = plot_dendrogram(
            clustering, df_pvals_temp, 
            dist, horizontal, pos_col, chain_col, 
            title, dend_filename, 
            dendrogram_subplots_kwargs, save_type,
            cluster_colors=cluster_colors)
        
        # NEW: Save detailed color mapping with active/inactive status
        analysis_color_file = working_filedir / f"cluster_{score_type}/{prefix}_{name}_{pthr}_{int(dist)}A_color_mapping.json"
        with open(analysis_color_file, 'w') as f:
            json.dump(cluster_info, f, indent=2)
        #print(f"Saved color mapping: {analysis_color_file}")

        # CLUSTERS RESIDUES AND LENGTH OF EACH CLUSTER #
        df_pvals_clust_i = df_pvals_clust.loc[(df_pvals_clust[colname] == pthr), ].reset_index(drop=True)
        clust_indices = df_pvals_clust_i[f'{colname}_Clust_{str(int(dist))}A'].unique()

        txt_filename = working_filedir / f"cluster_{score_type}/{prefix}_{name}_Dendrogram_{pthr}_{str(int(dist))}A.txt"
        with open(txt_filename, "w") as f: 
            for c in clust_indices: 
                c_data = df_pvals_clust_i.loc[df_pvals_clust_i[f'{colname}_Clust_{str(int(dist))}A'] == c, ].reset_index(drop=True)

                # WRITE LENGTH, RANGE, AND ALL RESIDUES #
                if len(c_data) > 0: 
                    all_unipos = c_data[pos_col].tolist()
                    all_chains = c_data[chain_col].tolist()
                    
                    # Include color and status information
                    cluster_key = f"cluster_{int(c)}"
                    cluster_color = cluster_info['clusters'].get(cluster_key, {}).get('color', '#CCCCCC')
                    is_active = cluster_info['clusters'].get(cluster_key, {}).get('active', False)
                    status = "ACTIVE" if is_active else "GRAYED OUT"
                    
                    f.write(f'Cluster {c} (Color: {cluster_color}, Status: {status}) : Length {len(c_data)} :\n')
                    
                    f.write('   ')
                    for unipos, chain in zip(all_unipos, all_chains): 
                        f.write(f'{chain}-{unipos} ')
                    f.write(f'\n')

    return None

def generate_cluster_colors(n_clusters=50):
    """
    Generate maximally distinct colors suitable for scientific visualization.
    Combines hand-picked colors for small sets with algorithmic generation for larger sets.
    """
    import colorsys
    
    # Hand-picked distinct colors for first 12 clusters
    base_colors = [
        '#e6194B',  # Red
        '#3cb44b',  # Green
        '#ffe119',  # Yellow
        '#4363d8',  # Blue
        '#f58231',  # Orange
        '#911eb4',  # Purple
        '#42d4f4',  # Cyan
        '#f032e6',  # Magenta
        '#bfef45',  # Lime
        '#fabed4',  # Pink
        '#469990',  # Teal
        '#dcbeff',  # Lavender
    ]
    
    if n_clusters <= len(base_colors):
        return base_colors[:n_clusters]
    
    # For more clusters, start with base and add generated ones
    colors = base_colors.copy()
    
    # Generate additional colors using golden ratio
    golden_ratio = 0.618033988749895
    h = 0.5  # Start at a different hue
    
    for i in range(n_clusters - len(base_colors)):
        h = (h + golden_ratio) % 1.0
        
        # Vary lightness and saturation more dramatically
        pattern = i % 6
        if pattern == 0:
            s, l = 0.95, 0.35
        elif pattern == 1:
            s, l = 0.60, 0.75
        elif pattern == 2:
            s, l = 0.90, 0.55
        elif pattern == 3:
            s, l = 0.70, 0.45
        elif pattern == 4:
            s, l = 0.85, 0.65
        else:
            s, l = 0.75, 0.85
        
        r, g, b = colorsys.hls_to_rgb(h, l, s)
        hex_color = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
        colors.append(hex_color)
    
    return colors

def plot_cluster_distance(
    distances, 
    yvalues, 
    names, 
    input_gene, 
    clust_filename, 
    plot_filename, 
    subplots_kwargs, 
    save_type, 
): 
    # print(distances,yvalues,names)
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
    clustering, 
    df_pvals_temp, 
    dist, 
    horizontal, 
    pos_col, 
    chain_col, 
    title, 
    dend_filename, 
    subplots_kwargs, 
    save_type,
    cluster_colors=None,
):  
    """
    Plot hierarchical clustering dendrogram with consistent cluster coloring.
    Only colors links below the distance threshold.
    """
    from scipy.cluster.hierarchy import dendrogram
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig, ax = plt.subplots(**subplots_kwargs)
    counts = np.zeros(clustering.children_.shape[0])
    n_samples = len(clustering.labels_)
    
    # Build linkage matrix from AgglomerativeClustering result
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
    
    # Create labels with cluster information
    xlbl_pos = list(df_pvals_temp[pos_col])
    xlbl_chain = list(df_pvals_temp[chain_col])
    xlbl_cluster = list(clustering.labels_)
    
    xlbl = [f'{pos}-{chain} (C{clust})' 
            for pos, chain, clust in zip(xlbl_pos, xlbl_chain, xlbl_cluster)]
    
    # Generate colors if not provided
    if cluster_colors is None or len(cluster_colors) < len(xlbl_cluster):
        cluster_colors = generate_cluster_colors(len(xlbl_cluster))
    
    # Helper function: map cluster ID to color
    def get_cluster_color(cluster_id):
        """Always return the same color for the same cluster ID"""
        return cluster_colors[int(cluster_id)]
    
    # Helper function: trace node to its cluster
    def get_cluster_for_node(node_id, linkage_mat, labels):
        """Find which cluster a node belongs to by tracing to leaves."""
        if node_id < len(labels):
            # It's a leaf node - return its cluster label
            return labels[node_id]
        else:
            # It's an internal node - trace down to left child
            merge_idx = int(node_id - len(labels))
            if merge_idx < len(linkage_mat):
                left_child = int(linkage_mat[merge_idx, 0])
                return get_cluster_for_node(left_child, linkage_mat, labels)
        return 0
    
    # Helper function: get distance for a node
    def get_node_distance(node_id, linkage_mat, n_samples):
        """Get the distance at which this node was formed."""
        if node_id < n_samples:
            # Leaf nodes have distance 0
            return 0.0
        else:
            # Internal node - look up merge distance
            merge_idx = int(node_id - n_samples)
            if merge_idx < len(linkage_mat):
                return linkage_mat[merge_idx, 2]
        return 0.0
    
    # Color function for dendrogram links
    def link_color_func(node_id):
        """
        Determine color for each link in the dendrogram.
        - If link distance <= threshold: use cluster color
        - If link distance > threshold: use gray
        """
        # Get the distance at which this node was formed
        node_distance = get_node_distance(node_id, linkage_matrix, n_samples)
        
        # If above threshold, return gray
        if node_distance > dist:
            return '#CCCCCC'
        
        # If below threshold, color by cluster
        cluster_id = get_cluster_for_node(node_id, linkage_matrix, clustering.labels_)
        return get_cluster_color(cluster_id)
    
    # Plot dendrogram with custom coloring
    dend_result = dendrogram(
        linkage_matrix,
        color_threshold=dist,
        above_threshold_color='#CCCCCC',
        labels=xlbl,
        orientation='right' if horizontal else 'top',
        leaf_rotation=90. if not horizontal else 0.,
        ax=ax,
        link_color_func=link_color_func
    )
    
    # ANALYZE THE ACTUAL COLORS USED IN THE DENDROGRAM
    # dend_result['color_list'] contains the actual colors assigned to each link
    actual_colors_used = set(dend_result.get('color_list', []))
    
    # Check which clusters were actually displayed in gray vs color
    all_cluster_ids = np.unique(clustering.labels_)
    cluster_color_status = {}
    
    for cluster_id in all_cluster_ids:
        expected_color = get_cluster_color(cluster_id).upper()
        # Check if this cluster's color appears in the dendrogram
        # Convert to uppercase for comparison
        cluster_is_colored = any(
            expected_color == actual_color.upper() 
            for actual_color in actual_colors_used
        )
        cluster_color_status[cluster_id] = cluster_is_colored
    
    # Count active vs grayed
    active_clusters = [cid for cid, is_colored in cluster_color_status.items() if is_colored]
    grayed_clusters = [cid for cid, is_colored in cluster_color_status.items() if not is_colored]
    
    # Create detailed cluster information
    cluster_info = {
        "threshold_distance": float(dist),
        "total_clusters": int(len(all_cluster_ids)),
        "active_clusters": int(len(active_clusters)),
        "grayed_clusters": int(len(grayed_clusters)),
        "clusters": {}
    }
    
    # Populate cluster information based on actual dendrogram display
    for cluster_id in all_cluster_ids:
        cluster_key = f"cluster_{int(cluster_id)}"
        is_active = cluster_color_status[cluster_id]
        assigned_color = get_cluster_color(cluster_id)
        
        cluster_info["clusters"][cluster_key] = {
            "cluster_id": int(cluster_id),
            "color": assigned_color,
            "active": is_active,
            "displayed_color": assigned_color if is_active else "#CCCCCC",
            "status": "active" if is_active else "grayed_out"
        }

    # Styling
    ax.set_title(title, fontsize=10, pad=10)
    if horizontal:
        ax.set_xlabel('Distance (Å)', fontsize=8)
        ax.tick_params(axis='y', labelsize=6)
        ax.axvline(x=dist, color='red', linestyle='--', linewidth=1, alpha=0.5, label=f'Threshold: {dist}Å')
    else:
        ax.set_ylabel('Distance (Å)', fontsize=8)
        ax.tick_params(axis='x', labelsize=6, rotation=90)
        ax.axhline(y=dist, color='red', linestyle='--', linewidth=1, alpha=0.5, label=f'Threshold: {dist}Å')
    
    ax.legend(loc='best', fontsize=6)
    
    plt.tight_layout()
    plt.savefig(dend_filename, dpi=100, transparent=False, format=save_type, bbox_inches='tight')
    plt.close()
    
    return cluster_info

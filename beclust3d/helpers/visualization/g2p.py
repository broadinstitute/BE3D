import os
import glob
import pandas as pd


def g2p_formatted_hit_cluster(results_dir,
                              gene_list: list,
                              screen_names: list,
                              lfc_pthr='05',
                              lfc3d_pthr='05',
                              meta_pthr='001',
                              dist=6,
                              meta=False,
                              conservation=False,
                              input_gene=None):
    """
    Description
        Prepare TSV for Hit Cluster Visualization on G2P

    Parameters
    ----------
    results_dir : str
        Path to the working directory where output files and results will be saved.

    gene_list : list
        Name of the gene being processed. 

    screen_names : list
        List of screens

    lfc_pthr : str (default = 05)
        P-value threshold for LFC hit prioritization, choices : [05,01,001]

    lfc3d_pthr : str (default = 05)
        P-value threshold for LFC3D hit prioritization, choices : [05,01,001]

    meta_pthr : str (default = 001)
        P-value threshold for Meta-aggregation hit prioritization, choices : [05,01,001]

    dist : int (default = 6)
        Type the distance cutoff for agglomerative clustering, which is only for title
        
    meta : bool
        True for getting Meta-aggregation analysis

    conservation : bool
        True for changing gene_list to ['Merged']
       
    Returns
    -------
    None
    """
    
    # screen_names = [f'{x}' for x in screen_names.split(',')]
    # This will gather distance-based cluster information across LFC, LFC3D and Union clustering approaches.
    result_both_directions_pd = pd.DataFrame()
    result_pos_pd = pd.DataFrame()
    result_neg_pd = pd.DataFrame()
    
    for gene_name, screen_name in zip(gene_list, screen_names):
        print(gene_name,screen_name)
        lfc_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'cluster_LFC',f'{gene_name}_{screen_name}_Aggr_Hits.tsv'))[0],sep='\t').set_index('unipos')
        lfc3d_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'cluster_LFC3D',f'{gene_name}_{screen_name}_Aggr_Hits.tsv'))[0],sep='\t').set_index('unipos')        
        union_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'cluster_union',f'{gene_name}_{screen_name}_Aggr_Hits.tsv'))[0],sep='\t').set_index('unipos')
        lfc_scores_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'LFC',f'{input_gene}_NonAggr_LFC.tsv'))[0],sep='\t').set_index('unipos')
        lfc3d_scores_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'LFC3D',f'{input_gene}_NonAggr_LFC3D.tsv'))[0],sep='\t').set_index('unipos')
        
        lfc_scores_pd = lfc_scores_pd.replace('-',None)
        lfc3d_scores_pd = lfc3d_scores_pd.replace('-',None)
        
        for direction in ['pos','neg']:
            lfc_scores_pd[f'{screen_name}_LFC_{direction}'] = lfc_scores_pd[f'{screen_name}_LFC_{direction}'].astype(float)
            lfc3d_scores_pd[f'{screen_name}_LFC3D_{direction}'] = lfc3d_scores_pd[f'{screen_name}_LFC3D_{direction}'].astype(float)            
            result_both_directions_pd = pd.concat([result_both_directions_pd,lfc_pd.filter(regex=rf'unires|{direction}_{lfc_pthr}_psig$|{direction}_{lfc_pthr}_psig_Clust_{dist}A$'),lfc3d_pd.filter(regex=rf'{direction}_{lfc3d_pthr}_psig$|{direction}_{lfc3d_pthr}_psig_Clust_{dist}A$'),union_pd.filter(regex=rf'{direction}_{lfc3d_pthr}_psig$|{direction}_{lfc3d_pthr}_psig_Clust_{dist}A$')],axis=1)  

        for direction in ['pos']:                        
            result_pos_pd = pd.concat([result_pos_pd,lfc_pd.filter(regex=rf'unires|{direction}_{lfc_pthr}_psig$|{direction}_{lfc_pthr}_psig_Clust_{dist}A$'),lfc3d_pd.filter(regex=rf'{direction}_{lfc3d_pthr}_psig$|{direction}_{lfc3d_pthr}_psig_Clust_{dist}A$'),union_pd.filter(regex=rf'{direction}_{lfc3d_pthr}_psig$|{direction}_{lfc3d_pthr}_psig_Clust_{dist}A$'),lfc_scores_pd.filter(regex=rf'.*{screen_name}_LFC_{direction}$'),lfc3d_scores_pd.filter(regex=rf'.*{screen_name}_LFC3D_{direction}$')],axis=1)
            
            rename_dict = {f'{screen_name}_LFC_{direction}_{lfc_pthr}_psig':f'{screen_name} LFC hit',
                           f'{screen_name}_LFC3D_{direction}_{lfc3d_pthr}_psig':f'{screen_name} LFC3D hit',
                           f'{screen_name}_union_{direction}_{lfc3d_pthr}_psig':f'{screen_name} union hit',
                            f'{screen_name}_LFC_{direction}_{lfc_pthr}_psig_Clust_{dist}A':f'{screen_name} LFC hit cluster',
                           f'{screen_name}_LFC3D_{direction}_{lfc3d_pthr}_psig_Clust_{dist}A':f'{screen_name} LFC3D hit cluster',
                           f'{screen_name}_union_{direction}_{lfc3d_pthr}_psig_Clust_{dist}A':f'{screen_name} union hit cluster'                           
                           }
            result_pos_pd = result_pos_pd.rename(rename_dict,axis=1)

        for direction in ['neg']:
            result_neg_pd = pd.concat([result_neg_pd,lfc_pd.filter(regex=rf'unires|{direction}_{lfc_pthr}_psig$|{direction}_{lfc_pthr}_psig_Clust_{dist}A$'),lfc3d_pd.filter(regex=rf'{direction}_{lfc3d_pthr}_psig$|{direction}_{lfc3d_pthr}_psig_Clust_{dist}A$'),union_pd.filter(regex=rf'{direction}_{lfc3d_pthr}_psig$|{direction}_{lfc3d_pthr}_psig_Clust_{dist}A$'),lfc_scores_pd.filter(regex=rf'.*{screen_name}_LFC_{direction}$'),lfc3d_scores_pd.filter(regex=rf'.*{screen_name}_LFC3D_{direction}$')],axis=1)           
            rename_dict = {f'{screen_name}_LFC_{direction}_{lfc_pthr}_psig':f'{screen_name} LFC hit',
                           f'{screen_name}_LFC3D_{direction}_{lfc3d_pthr}_psig':f'{screen_name} LFC3D hit',
                           f'{screen_name}_union_{direction}_{lfc3d_pthr}_psig':f'{screen_name} union hit',
                            f'{screen_name}_LFC_{direction}_{lfc_pthr}_psig_Clust_{dist}A':f'{screen_name} LFC hit cluster',
                           f'{screen_name}_LFC3D_{direction}_{lfc3d_pthr}_psig_Clust_{dist}A':f'{screen_name} LFC3D hit cluster',
                           f'{screen_name}_union_{direction}_{lfc3d_pthr}_psig_Clust_{dist}A':f'{screen_name} union hit cluster'                           
                           }
            result_neg_pd = result_neg_pd.rename(rename_dict,axis=1)
            
        # Create 'Pos'/'Neg' labels conditionally and combine them
        pos_mask = lfc_pd[f'{screen_name}_LFC_pos_{lfc_pthr}_psig'] == f'p<0.{lfc_pthr}'
        neg_mask = lfc_pd[f'{screen_name}_LFC_neg_{lfc_pthr}_psig'] == f'p<0.{lfc_pthr}'

        # Initialize the hits column with None
        lfc_pd[f'{screen_name}_LFC_{lfc_pthr}_hits'] = None

        # Fill in 'Pos' and 'Neg' based on masks
        lfc_pd.loc[pos_mask, f'{screen_name}_LFC_{lfc_pthr}_hits'] = 'Pos'
        lfc_pd.loc[neg_mask, f'{screen_name}_LFC_{lfc_pthr}_hits'] = 'Neg'

        # Optional: flag both if both are significant
        lfc_pd.loc[pos_mask & neg_mask, f'{screen_name}_LFC_{lfc_pthr}_hits'] = 'Pos+Neg'

        # Create 'Pos'/'Neg' labels conditionally and combine them
        pos_mask = lfc3d_pd[f'{screen_name}_LFC3D_pos_{lfc3d_pthr}_psig'] == f'p<0.{lfc3d_pthr}'
        neg_mask = lfc3d_pd[f'{screen_name}_LFC3D_neg_{lfc3d_pthr}_psig'] == f'p<0.{lfc3d_pthr}'

        # Initialize the hits column with None
        lfc3d_pd[f'{screen_name}_LFC3D_{lfc3d_pthr}_hits'] = None

        # Fill in 'Pos' and 'Neg' based on masks
        lfc3d_pd.loc[pos_mask, f'{screen_name}_LFC3D_{lfc3d_pthr}_hits'] = 'Pos'
        lfc3d_pd.loc[neg_mask, f'{screen_name}_LFC3D_{lfc3d_pthr}_hits'] = 'Neg'

        # Optional: flag both if both are significant
        lfc3d_pd.loc[pos_mask & neg_mask, f'{screen_name}_LFC3D_{lfc3d_pthr}_hits'] = 'Pos+Neg'
        result_both_directions_pd = pd.concat([result_both_directions_pd, lfc_pd[f'{screen_name}_LFC_{lfc_pthr}_hits'],lfc3d_pd[f'{screen_name}_LFC3D_{lfc3d_pthr}_hits']],axis=1)
                
    if meta:
        if conservation:
            gene_name = 'Merged'
        screen_name = 'SUM'
        lfc_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'cluster_LFC',f'{gene_name}_Meta_Aggr_Hits.tsv'))[0],sep='\t').set_index('unipos')
        lfc3d_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'cluster_LFC3D',f'{gene_name}_Meta_Aggr_Hits.tsv'))[0],sep='\t').set_index('unipos')
        union_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'cluster_union',f'{gene_name}_Meta_Aggr_Hits.tsv'))[0],sep='\t').set_index('unipos')
        lfc_scores_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'meta-aggregate',f'{gene_name}_MetaAggr_LFC.tsv'))[0],sep='\t').set_index('unipos')
        lfc3d_scores_pd = pd.read_csv(glob.glob(os.path.join(results_dir,'meta-aggregate',f'{gene_name}_MetaAggr_LFC3D.tsv'))[0],sep='\t').set_index('unipos')
        
        lfc_scores_pd = lfc_scores_pd.replace('-',None)
        lfc3d_scores_pd = lfc3d_scores_pd.replace('-',None)

        for direction in ['pos','neg']:
            lfc_scores_pd[f'{screen_name}_LFC_{direction}'] = lfc_scores_pd[f'{screen_name}_LFC_{direction}'].astype(float)
            lfc3d_scores_pd[f'{screen_name}_LFC3D_{direction}'] = lfc3d_scores_pd[f'{screen_name}_LFC3D_{direction}'].astype(float)                    
            result_both_directions_pd = pd.concat([result_both_directions_pd,lfc_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),lfc3d_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),union_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$')],axis=1)
            
        for direction in ['pos']:
            result_pos_pd = pd.concat([result_pos_pd,lfc_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),lfc3d_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),union_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),lfc_scores_pd.filter(regex=rf'.*{screen_name}_LFC_{direction}$'),lfc3d_scores_pd.filter(regex=rf'.*{screen_name}_LFC3D_{direction}$')],axis=1)
            rename_dict = {f'{screen_name}_LFC_{direction}_{meta_pthr}_psig':f'Meta LFC hit',
                           f'{screen_name}_LFC3D_{direction}_{meta_pthr}_psig':f'Meta LFC3D hit',
                           f'{screen_name}_union_{direction}_{meta_pthr}_psig':f'Meta union hit',
                           f'{screen_name}_LFC_{direction}_{meta_pthr}_psig_Clust_{dist}A':f'Meta LFC hit cluster',
                           f'{screen_name}_LFC3D_{direction}_{meta_pthr}_psig_Clust_{dist}A':f'Meta LFC3D hit cluster',
                           f'{screen_name}_union_{direction}_{meta_pthr}_psig_Clust_{dist}A':f'Meta union hit cluster',
                           f'{screen_name}_LFC_{direction}': 'Meta LFC',
                           f'{screen_name}_LFC3D_{direction}': 'Meta LFC3D'
                           }
            result_pos_pd = result_pos_pd.rename(rename_dict,axis=1)
            result_pos_pd.rename(columns={col: col.replace("_pos", "") for col in result_pos_pd.columns if "_pos" in col}, inplace=True)
            
        for direction in ['neg']:
            result_neg_pd = pd.concat([result_neg_pd,lfc_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),lfc3d_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'), union_pd.filter(regex=rf'{direction}_{meta_pthr}_psig$|{direction}_{meta_pthr}_psig_Clust_{dist}A$'),lfc_scores_pd.filter(regex=rf'.*{screen_name}_LFC_{direction}$'),lfc3d_scores_pd.filter(regex=rf'.*{screen_name}_LFC3D_{direction}$')],axis=1)        
            rename_dict = {f'{screen_name}_LFC_{direction}_{meta_pthr}_psig':f'Meta LFC hit',
                           f'{screen_name}_LFC3D_{direction}_{meta_pthr}_psig':f'Meta LFC3D hit',
                           f'{screen_name}_union_{direction}_{meta_pthr}_psig':f'Meta union hit',
                            f'{screen_name}_LFC_{direction}_{meta_pthr}_psig_Clust_{dist}A':f'Meta LFC hit cluster',
                           f'{screen_name}_LFC3D_{direction}_{meta_pthr}_psig_Clust_{dist}A':f'Meta LFC3D hit cluster',
                           f'{screen_name}_union_{direction}_{meta_pthr}_psig_Clust_{dist}A':f'Meta union hit cluster',
                           f'{screen_name}_LFC_{direction}': 'Meta LFC',
                           f'{screen_name}_LFC3D_{direction}': 'Meta LFC3D'                           
                           }
            result_neg_pd = result_neg_pd.rename(rename_dict,axis=1)       
            result_neg_pd.rename(columns={col: col.replace("_neg", "") for col in result_neg_pd.columns if "_neg" in col}, inplace=True)

        # Create 'Pos'/'Neg' labels conditionally and combine them
        pos_mask = lfc_pd[f'{screen_name}_LFC_pos_{meta_pthr}_psig'] == f'p<0.{meta_pthr}'
        neg_mask = lfc_pd[f'{screen_name}_LFC_neg_{meta_pthr}_psig'] == f'p<0.{meta_pthr}'

        # Initialize the hits column with None
        lfc_pd[f'{screen_name}_LFC_{meta_pthr}_hits'] = None

        # Fill in 'Pos' and 'Neg' based on masks
        lfc_pd.loc[pos_mask, f'{screen_name}_LFC_{meta_pthr}_hits'] = 'Pos'
        lfc_pd.loc[neg_mask, f'{screen_name}_LFC_{meta_pthr}_hits'] = 'Neg'

        # Optional: flag both if both are significant
        lfc_pd.loc[pos_mask & neg_mask, f'{screen_name}_LFC_{meta_pthr}_hits'] = 'Pos+Neg'

        # Create 'Pos'/'Neg' labels conditionally and combine them
        pos_mask = lfc3d_pd[f'{screen_name}_LFC3D_pos_{meta_pthr}_psig'] == f'p<0.{meta_pthr}'
        neg_mask = lfc3d_pd[f'{screen_name}_LFC3D_neg_{meta_pthr}_psig'] == f'p<0.{meta_pthr}'

        # Initialize the hits column with None
        lfc3d_pd[f'{screen_name}_LFC3D_{meta_pthr}_hits'] = None

        # Fill in 'Pos' and 'Neg' based on masks
        lfc3d_pd.loc[pos_mask, f'{screen_name}_LFC3D_{meta_pthr}_hits'] = 'Pos'
        lfc3d_pd.loc[neg_mask, f'{screen_name}_LFC3D_{meta_pthr}_hits'] = 'Neg'

        # Optional: flag both if both are significant
        lfc3d_pd.loc[pos_mask & neg_mask, f'{screen_name}_LFC3D_{meta_pthr}_hits'] = 'Pos+Neg'
        
        result_both_directions_pd = pd.concat([result_both_directions_pd, lfc_pd[f'{screen_name}_LFC_{meta_pthr}_hits'],lfc3d_pd[f'{screen_name}_LFC3D_{meta_pthr}_hits']],axis=1)

    result_pos_pd.rename(columns={col: col.replace("_pos", "") for col in result_pos_pd.columns if "_pos" in col}, inplace=True)           
    result_neg_pd.rename(columns={col: col.replace("_neg", "") for col in result_neg_pd.columns if "_neg" in col}, inplace=True)

    os.makedirs(os.path.join(results_dir,'g2p_visualization'),exist_ok=True)
    result_both_directions_pd.to_csv(os.path.join(results_dir,'g2p_visualization','lfc_lfc3d_union_cluster_for_g2p.tsv'),sep='\t')
    result_pos_pd.to_csv(os.path.join(results_dir,'g2p_visualization','POS_lfc_lfc3d_union_cluster_for_g2p.tsv'),sep='\t')
    result_neg_pd.to_csv(os.path.join(results_dir,'g2p_visualization','NEG_lfc_lfc3d_union_cluster_for_g2p.tsv'),sep='\t')
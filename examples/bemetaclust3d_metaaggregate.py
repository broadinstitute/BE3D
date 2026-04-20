"""
Step 3a: Meta-aggregation and meta-level clustering.
Usage example: python bemetaclust3d_metaaggregate.py ./yaml/dnmt3a_local.yaml

Requires Steps 1 and 2 to have been run first.
Only runs meta-aggregation when the number of screens > 1.
"""

import os
import sys
import pandas as pd

from be3d_helper import get_required, get_optional, load_config, find_union


def main(**kwargs):
    # REQUIRED
    input_gene       = kwargs['input_gene']
    input_uniprot    = kwargs['input_uniprot']
    screens          = kwargs['screens']
    output_dir       = kwargs['output_dir']
    conservation_run = kwargs['conservation_run']
    alt_gene_name    = kwargs['alt_gene_name']
    alt_uniprot_id   = kwargs['alt_uniprot_id']
    alt_screen_start = kwargs['alt_screen_start']

    # OPTIONAL
    user_pdb                = kwargs['user_pdb']
    function_for_meta       = kwargs['function_for_meta']
    nRandom                 = kwargs['nRandom']
    multi_screen_pthr       = kwargs['multi_screen_pthr']
    clustering_radius       = kwargs['clustering_radius']
    qa_passed_only          = kwargs['qa_passed_only']
    qa_only                 = kwargs['qa_only']
    qa_controls             = kwargs['qa_controls']
    qa_cases                = kwargs['qa_cases']
    priority_on_alternative = kwargs['priority_on_alternative']
    atom_level_naa          = kwargs['atom_level_naa']

    if user_pdb:
        structureid = f'PDB-{input_uniprot}'
    else:
        structureid = f'AF-{input_uniprot}-F1-model_v6'

    pdb_file = os.path.join(output_dir, "sequence_structure", f"{structureid}_processed.pdb")

    print(f'All results will be saved in the following directory: {output_dir}')

    if isinstance(screens, str): screens = [screen.strip() for screen in screens.split(',')]
    screen_names = [s.split('.')[0] for s in screens]

    df_residuemap = pd.DataFrame()
    conserv_dfs = []
    gene_list = []
    if conservation_run:
        df_residuemap_filename = os.path.join(output_dir, f"conservation/{input_gene}{input_uniprot}_{alt_gene_name}{alt_uniprot_id}_residuemap_conservation.tsv")
        df_residuemap = pd.read_csv(df_residuemap_filename, sep='\t')

        if priority_on_alternative:
            for screen_name in screen_names:
                conserv_dfs.append(df_residuemap)
                gene_list.append(alt_gene_name)
        else:
            for screen_name in screen_names:
                if alt_gene_name and screen_name.startswith(alt_screen_start):
                    conserv_dfs.append(df_residuemap)
                    gene_list.append(alt_gene_name)
                else:
                    conserv_dfs.append(None)
                    gene_list.append(input_gene)
    else:
        for screen_name in screen_names:
            conserv_dfs.append(None)
            gene_list.append(input_gene)

    if qa_only:
        print("qa_only=True: skipping meta-aggregate step.")
        sys.exit()

    if qa_passed_only:
        h2_ks_test_pd = pd.read_csv(f'{output_dir}/hypothesis_qc/KolmogorovSmirnov_hypothesis2.tsv', sep='\t')
        h2_ks_test_pd = h2_ks_test_pd.replace(-999, None)
        white_screen_list = h2_ks_test_pd[
            (h2_ks_test_pd[f"p_{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}"] < 0.05) &
            (h2_ks_test_pd['gene_name'].isin(gene_list))
        ]['screenid'].to_list()
        print(f'QA filter: {len(list(set(screen_names).intersection(white_screen_list)))} screens passing.')
        screen_names = list(set(screen_names).intersection(white_screen_list))
        conserv_dfs = []
        gene_list = []
        for screen_name in screen_names:
            if screen_name.startswith(alt_screen_start):
                conserv_dfs.append(df_residuemap)
                gene_list.append(alt_gene_name)
            else:
                conserv_dfs.append(None)
                gene_list.append(input_gene)

    if len(screen_names) <= 1:
        print("Only 1 screen present — meta-aggregation requires >1 screen. Exiting.")
        sys.exit()

    # Determine gene label used in Step 2 merged outputs
    gene = 'Merged' if conservation_run else input_gene

    # Load the merged LFC/LFC3D table written at the end of Step 2
    merged_df_LFC_LFC3D = pd.read_csv(f'{output_dir}/LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr.tsv.gz', sep='\t')

    # META-AGGREGATION

    # LFC3D meta-aggregate
    df_bidir_meta = average_split_meta(merged_df_LFC_LFC3D, output_dir, gene, screen_names, nRandom=nRandom, score_type='LFC3D', aggr_func_name=function_for_meta)
    df_dis, _, _ = bin_meta(df_bidir_meta, output_dir, gene, score_type='LFC3D', aggr_func_name=function_for_meta)
    znorm_meta(df_dis, output_dir, gene, screen_names, pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', aggr_func_name=function_for_meta)
    print(f'Output from beclust3d.aggregate.metaaggregate is saved in the following directory: {output_dir}meta-aggregate/')

    df_lfc3d = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC3D.tsv', sep='\t')
    average_split_bin_plots(
        df_lfc3d, workdir=output_dir, input_gene=gene,
        screen_name='', func=function_for_meta, pthr=multi_screen_pthr,
        score_type='LFC3D', aggregate_dir='meta-aggregate', save_type='svg'
    )
    print(f'Output from beclust3d.aggregate.aggregate_plot is saved in the following directory: {output_dir}meta-aggregate/')

    # LFC meta-aggregate
    df_bidir_meta = average_split_meta(merged_df_LFC_LFC3D, output_dir, gene, screen_names, nRandom=nRandom, score_type='LFC', aggr_func_name=function_for_meta)
    df_dis = bin_meta(df_bidir_meta, output_dir, gene, score_type='LFC', aggr_func_name=function_for_meta)
    znorm_meta(df_bidir_meta, output_dir, gene, screen_names,pthrs=[0.05, 0.01, 0.001], score_type='LFC', aggr_func_name=function_for_meta)
    print(f'Output from beclust3d.aggregate.metaaggregate is saved in the following directory: {output_dir}meta-aggregate/')

    df_lfc = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC.tsv', sep='\t')
    for screen_name in screen_names:
        average_split_bin_plots(
            df_lfc, workdir=output_dir, input_gene=gene,
            screen_name='', func=function_for_meta, pthr=multi_screen_pthr,
            score_type='LFC', aggregate_dir='meta-aggregate', save_type='svg'
        )
    print(f'Output from beclust3d.aggregate.aggregate_plot is saved in the following directory: {output_dir}meta-aggregate/')

    # META CLUSTERING
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    for score_type in ['LFC', 'LFC3D']:
        df_pvals = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_{score_type}.tsv', sep='\t')

        pref = f'{function_for_meta}_{score_type}'
        df_hits_clust, distances, yvalues = clustering(
            df_struc, df_pvals,
            output_dir, gene,
            psig_columns=[f'{pref}_neg_05_psig', f'{pref}_pos_05_psig', f'{pref}_neg_01_psig',
                          f'{pref}_pos_01_psig', f'{pref}_neg_001_psig', f'{pref}_pos_001_psig'],
            pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
            screen_name='Meta', score_type=score_type,
            max_distances=20, merge_cols=['unipos', 'chain'],
            atom_level=pdb_file if atom_level_naa == True else False,
        )
        plot_clustering(
            df_struc, df_pvals, df_hits_clust, clustering_radius,
            output_dir, gene, distances, yvalues,
            names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
            psig_columns=[f'{pref}_neg_05_psig', f'{pref}_pos_05_psig', f'{pref}_neg_01_psig',
                          f'{pref}_pos_01_psig', f'{pref}_neg_001_psig', f'{pref}_pos_001_psig'],
            pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
            screen_name='Meta', score_type=score_type,
            merge_col=['unipos', 'chain'],
            save_type='svg',
            dendrogram_subplots_kwargs={'figsize': (15, 3.5)},
        )
        print(f'Output from beclust3d.lfc3d.clustering is saved in the following directory: {output_dir}cluster_{score_type}/')

    # Union clustering
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
    df_pvals_LFC3D = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC3D.tsv', sep='\t')
    df_pvals_LFC   = pd.read_csv(f'{output_dir}/meta-aggregate/{gene}_MetaAggr_LFC.tsv', sep='\t')
    df_pvals = pd.concat([df_pvals_LFC3D, df_pvals_LFC.drop(columns=['unipos', 'unires', 'chain'])], axis=1)

    for each_pthr in ['05', '01', '001']:
        df_pvals[f'{function_for_meta}_union_neg_{each_pthr}_psig'] = df_pvals[
            [f'{function_for_meta}_LFC_neg_{each_pthr}_psig', f'{function_for_meta}_LFC3D_neg_{each_pthr}_psig']
        ].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)
        df_pvals[f'{function_for_meta}_union_pos_{each_pthr}_psig'] = df_pvals[
            [f'{function_for_meta}_LFC_pos_{each_pthr}_psig', f'{function_for_meta}_LFC3D_pos_{each_pthr}_psig']
        ].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)

    df_hits_clust, distances, yvalues = clustering(
        df_struc, df_pvals, output_dir, gene,
        psig_columns=[f'{function_for_meta}_union_neg_05_psig', f'{function_for_meta}_union_pos_05_psig',
                      f'{function_for_meta}_union_neg_01_psig', f'{function_for_meta}_union_pos_01_psig',
                      f'{function_for_meta}_union_neg_001_psig', f'{function_for_meta}_union_pos_001_psig'],
        pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
        screen_name='Meta', score_type='union',
        max_distances=20, merge_cols=['unipos', 'chain'],
        atom_level=pdb_file if atom_level_naa == True else False
    )
    plot_clustering(
        df_struc, df_pvals, df_hits_clust, clustering_radius,
        output_dir, gene, distances, yvalues,
        names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
        psig_columns=[f'{function_for_meta}_union_neg_05_psig', f'{function_for_meta}_union_pos_05_psig',
                      f'{function_for_meta}_union_neg_01_psig', f'{function_for_meta}_union_pos_01_psig',
                      f'{function_for_meta}_union_neg_001_psig', f'{function_for_meta}_union_pos_001_psig'],
        pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
        screen_name='Meta', score_type='union',
        merge_col=['unipos', 'chain'],
        save_type='svg',
        dendrogram_subplots_kwargs={'figsize': (15, 3.5)}
    )
    print(f'Output from beclust3d.lfc3d.clustering is saved in the following directory: {output_dir}cluster_union/')

    print("Step 3a complete: meta-aggregation and clustering done")


if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    beclust3d_path = get_optional(config, 'beclust3d_path', str, '.')
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

    from beclust3d.lfc3d.clustering import clustering
    from beclust3d.lfc3d.clustering_plot import plot_clustering
    from beclust3d.aggregate.metaaggregate import average_split_meta, bin_meta, znorm_meta
    from beclust3d.aggregate.aggregate_plot import average_split_bin_plots

    # REQUIRED
    input_gene    = get_required(config, 'input_gene',    str)
    input_uniprot = get_required(config, 'input_uniprot', str)
    screens       = get_required(config, 'screens',       (str, list))
    output_dir    = get_optional(config, 'output_dir',    (str, type(None)), './')

    conservation_run  = get_required(config, 'conservation.run',               bool)
    alt_gene_name     = get_required(config, 'conservation.alt_gene_name',     (str, type(None)))
    alt_uniprot_id    = get_required(config, 'conservation.alt_uniprot_id',    (str, type(None)))
    alt_screen_start  = get_required(config, 'conservation.alt_screen_start',  (str, type(None)))

    # OPTIONAL
    user_pdb                = get_optional(config, 'user_pdb',                (str, type(None)),  None)
    function_for_meta       = get_optional(config, 'function_for_meta',       str,                'mean')
    nRandom                 = get_optional(config, 'nRandom',                 int,                500)
    multi_screen_pthr       = get_optional(config, 'pthr.multi_screen',       (int, float),       0.05)
    clustering_radius       = get_optional(config, 'clustering_radius',       (int, float),       6.0)
    qa_passed_only          = get_optional(config, 'qa.qa_passed_only',       bool,               False)
    qa_only                 = get_optional(config, 'qa.qa_only',              bool,               False)
    qa_controls             = get_optional(config, 'qa.controls',             (list, type(None)), ['No Mutation'])
    qa_cases                = get_optional(config, 'qa.cases',                (list, type(None)), ['Nonsense', 'Splice'])
    priority_on_alternative = get_optional(config, 'priority_on_alternative', bool,               False)
    atom_level_naa          = get_optional(config, 'atom_level_naa',          bool,               False)

    main(
        input_gene=input_gene, input_uniprot=input_uniprot,
        screens=screens, output_dir=output_dir,
        user_pdb=user_pdb, nRandom=nRandom,
        multi_screen_pthr=multi_screen_pthr,
        clustering_radius=clustering_radius,
        function_for_meta=function_for_meta,
        conservation_run=conservation_run, alt_gene_name=alt_gene_name,
        alt_uniprot_id=alt_uniprot_id, alt_screen_start=alt_screen_start,
        qa_passed_only=qa_passed_only, qa_only=qa_only,
        qa_controls=qa_controls, qa_cases=qa_cases,
        priority_on_alternative=priority_on_alternative,
        atom_level_naa=atom_level_naa,
    )

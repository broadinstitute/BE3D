import os
import sys
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main(**kwargs):
    ## REQUIRED
    input_gene=kwargs['input_gene']
    input_uniprot=kwargs['input_uniprot']
    input_chain=kwargs['input_chain']
    workdir=kwargs['workdir'] 
    screen_dir=kwargs['screen_dir']
    screens=kwargs['screens']
    mut_list_col=kwargs['mut_list_col']
    mut_col=kwargs['mut_col']
    val_col=kwargs['val_col']
    gene_col=kwargs['gene_col']
    edits_col=kwargs['edits_col']
    gRNA_col=kwargs['gRNA_col']  
    output_dir=kwargs['output_dir']
    
    ## OPTIONAL
    input_fasta=kwargs['user_fasta']
    user_pdb=kwargs['user_pdb']
    user_dssp=kwargs['user_dssp']
    nRandom=kwargs['nRandom']
    pthr=kwargs['pthr']
    structure_radius=kwargs['structure_radius']    
    clustering_radius=kwargs['clustering_radius']

    structureid = f'PDB-{input_uniprot}-F1-model_v4'    
    pthr_str = str(pthr)
    pthr_str_short = str(pthr).split('.')[1]


    ## ASSIGN VALUES FOR EMPTY VARIABLES
    if output_dir == '' or not output_dir:
        output_dir = os.path.join(workdir, f'{input_gene}')

    os.makedirs(output_dir, exist_ok=True)
    print('All results will be saved in the following directory:')
    print(output_dir)

    screens = [screen.strip() for screen in screens.split(',')]

    screen_names = [s.split('.')[0] for s in screens]
    input_dfs = [pd.read_csv(os.path.join(screen_dir,s), sep='\t') for s in screens]
        
    sequence_structural_features(
        output_dir,
        input_gene, input_uniprot, structureid,
        user_fasta=user_fasta,
        user_pdb=user_pdb,
        user_dssp=user_dssp,
        chains=[input_chain],
        radius=structure_radius,
        )

    hypothesis_test(
        output_dir,
        input_dfs, screen_names,
        cases=['Nonsense'],
        controls=['Silent', 'No Mutation'],
        comp_name='Nonsense_vs_SilentNoMut',
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
    )

    parse_be_data(
        output_dir,
        input_dfs, input_gene, screen_names,
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        edits_col=edits_col,
        mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"],
        mut_delimiter = ';',
        conserv_dfs = [None for _ in screen_names],
        )

    plot_rawdata(
        output_dir,
        input_dfs,
        screen_names,
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"],
        )

    df_missense_list = [
        pd.read_csv(f'{output_dir}/screendata/{input_gene}_{sn}_Missense.tsv',
                    sep='\t') for sn in screen_names
    ]

    for df_missense, screen_name in zip(df_missense_list, screen_names):
        randomize_data(
            df_missense,
            output_dir, input_gene,
            screen_name,
            nRandom=nRandom,
            seed=True,
            )

    ## PRIOTIZE
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
    df_consrv = None
    target_res_pos, target_res = 'original_res_pos', 'unires'

    for screen_name in screen_names:

        df_control = pd.read_csv(f'{output_dir}/screendata/{input_gene}_{screen_name}_No_Mutation.tsv', sep='\t', index_col=0)
        df_dict = {}

        for mut in ['Missense', 'Silent', 'Nonsense']:
            filepath = f'{output_dir}/screendata/{input_gene}_{screen_name}_{mut}.tsv'
            if os.path.exists(filepath):
                df_dict[mut] = pd.read_csv(filepath, sep='\t', index_col=0)

        df_missense = prioritize_by_sequence(
            df_struc, df_consrv, df_control,
            output_dir,
            input_gene, screen_name,
            df_dict,
            target_res_pos=target_res_pos, target_res=target_res,
        )

        df_rand = pd.read_csv(f'{output_dir}/screendata_rand/{input_gene}_{screen_name}_Missense_rand.tsv', sep='\t')

        randomize_sequence(
            df_missense, df_rand,
            output_dir,
            input_gene, screen_name,
            nRandom=nRandom, conservation=False,
            muttype='Missense',
            function_name='max',
            target_pos='unipos', target_res=None,
            )

        plot_screendata_sequence(
            df_missense,
            output_dir,
            input_gene, screen_name, function_name='max', muttype='Missense',
        )

        # file_path = f'{output_dir}/screendata_sequence/plots/{input_gene}_{screen_name}_Missense_lfcz_scatter_by_bin_posneg.png'

    ## CLUSTERING ON PRIORITIZED
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    for screen_name in screen_names:
        df_missense = pd.read_csv(f'{output_dir}/screendata_sequence/{input_gene}_{screen_name}_protein_edits.tsv', sep='\t')
        df_missense['mean_Missense_LFC_plab_input'] = df_missense['mean_Missense_LFC_p'].apply(lambda x: f'p<{pthr_str}' if x<pthr else f'p>={pthr_str}')

        df_hits_clust, distances, yvalues = clustering(
            df_struc, df_missense,
            output_dir, input_gene,
            psig_columns=[f'mean_Missense_LFC_plab_input'],
            pthr_cutoffs=[f'p<{pthr_str}'],
            screen_name=screen_name, score_type='lfc',
            max_distances=25, merge_cols=['unipos', 'chain'],
        )

        # PLOTTING #
        plot_clustering(
            df_struc, df_missense,
            df_hits_clust, clustering_radius,
            output_dir, input_gene,
            distances, yvalues,
            psig_columns=[f'mean_Missense_LFC_plab_input'],
            names=['Bidirection'],
            pthr_cutoffs=[f'p<{pthr_str}'],
            screen_name = screen_name, score_type='lfc',
            merge_col=['unipos', 'chain'],
        )

        # file_path = f'{output_dir}/cluster_lfc/plots/{input_gene}_{screen_name}_lfc_Bidirection_Dendogram_6A.png'

    ## LFC3D AND META-AGGREGATION ON LFC/LFC3D
    df_edits_list = []
    for screen_name in screen_names:
        df = pd.read_csv(f'{output_dir}/screendata_sequence/{input_gene}_{screen_name}_protein_edits.tsv', sep='\t')
        df_edits_list.append(df)
    df_rand_list = []
    for screen_name in screen_names:
        df = pd.read_csv(f'{output_dir}/screendata_sequence_rand/{input_gene}_{screen_name}_Missense_protein_edits_rand.tsv', sep='\t')
        df_rand_list.append(df)


    df_LFC_LFC3D = calculate_lfc3d(
        df_struc, df_edits_list, df_rand_list,
        output_dir, input_gene, screen_names,
        nRandom=nRandom,  muttype='Missense',
        function_type='max', function_aggr=np.mean,
        conserved_only=False,
    )

    # LFC3D #
    df_bidir_meta = average_split_meta(
        df_LFC_LFC3D,
        output_dir, input_gene, screen_names,
        nRandom=nRandom,
        score_type='LFC3D', aggr_func=np.sum, aggr_func_name='SUM',
    )
    df_dis, df_neg_stats, df_pos_stats = bin_meta(
        df_bidir_meta,
        output_dir, input_gene,
        score_type='LFC3D', aggr_func_name='SUM',
    )
    znorm_meta(
        df_dis, df_neg_stats, df_pos_stats,
        output_dir, input_gene,
        pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', aggr_func_name='SUM',
    )

    df_lfc3d = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC3D.tsv', sep='\t')

    for screen_name in screen_names:
        average_split_bin_plots(
            df_lfc3d,
            workdir = output_dir,
            input_gene = input_gene,
            screen_name='', # BLANK FOR META #
            func='SUM', # BLANK FOR NON AGGR #
            pthr=pthr,
            score_type='LFC3D',
            aggregate_dir='meta-aggregate',
            )

    # LFC #
    df_bidir_meta = average_split_meta(
        df_LFC_LFC3D,
        output_dir, input_gene, screen_names,
        nRandom=500,
        score_type='LFC', aggr_func=np.sum, aggr_func_name='SUM',
    )
    df_dis, df_neg_stats, df_pos_stats = bin_meta(
        df_bidir_meta,
        output_dir, input_gene,
        score_type='LFC', aggr_func_name='SUM',
    )
    df_z = znorm_meta(
        df_bidir_meta, df_neg_stats, df_pos_stats,
        output_dir, input_gene,
        pthrs=[0.05, 0.01, 0.001], score_type='LFC', aggr_func_name='SUM',
    )

    df_lfc = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC.tsv', sep='\t')

    for screen_name in screen_names:
        average_split_bin_plots(
            df_lfc,
            workdir = output_dir,
            input_gene = input_gene,
            screen_name='', # BLANK FOR META #
            func='SUM', # BLANK FOR NON AGGR #
            pthr=pthr,
            score_type='LFC',
            aggregate_dir='meta-aggregate',
            )
        
    ## CLUSTERING
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    for score_type in ['LFC', 'LFC3D']:
        df_pvals = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_{score_type}.tsv', sep='\t')

        df_hits_clust, distances, yvalues = clustering(
            df_struc, df_pvals,
            output_dir, input_gene,
            psig_columns=[f'SUM_{score_type}_neg_{pthr_str_short}_psig',
                        f'SUM_{score_type}_pos_{pthr_str_short}_psig'],
            pthr_cutoffs=[f'p<{pthr_str}', f'p<{pthr_str}'],
            screen_name='Meta', score_type=score_type,
            max_distances=25, merge_cols=['unipos', 'chain'],
        )

        # PLOTTING #
        plot_clustering(
            df_struc, df_pvals,
            df_hits_clust, clustering_radius,
            output_dir, input_gene,
            distances, yvalues,
            psig_columns=[f'SUM_{score_type}_neg_{pthr_str_short}_psig',
                        f'SUM_{score_type}_pos_{pthr_str_short}_psig'],
            names=['Negative', 'Positive'],
            pthr_cutoffs=[f'p<{pthr_str}', f'p<{pthr_str}'],
            screen_name='Meta', score_type=score_type,
            merge_col=['unipos', 'chain'],
        )

        # file_path = f'{output_dir}/meta-aggregate/plots/{input_gene}_{score_type}_cutoff{pthr_str_short}_scatter_colored.png'
        # file_path = f'{output_dir}/cluster_{score_type}/plots/{input_gene}_Meta_{score_type}_Positive_Dendogram_6A.png'

    ## CLUSTERING ON UNION
    def find_union(input):
        if input[0] == f'p<{pthr_str}' or input[1] == f'p<{pthr_str}':
            return f'p<{pthr_str}'
        else:
            return f'p>={pthr_str}'

    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    # Load both LFC and lFC3D dataframes
    df_pvals_LFC3D = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC3D.tsv', sep='\t')
    df_pvals_LFC = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC.tsv', sep='\t')
    df_pvals = pd.concat([df_pvals_LFC3D, df_pvals_LFC.drop(columns=['unipos', 'unires', 'chain'])], axis=1)

    # Find union of LFC and LFC3D
    df_pvals[f'SUM_union_neg_{pthr_str_short}_psig'] = df_pvals[[f'SUM_LFC_neg_{pthr_str_short}_psig', f'SUM_LFC3D_neg_{pthr_str_short}_psig']].apply(find_union, axis=1)
    df_pvals[f'SUM_union_pos_{pthr_str_short}_psig'] = df_pvals[[f'SUM_LFC_pos_{pthr_str_short}_psig', f'SUM_LFC3D_pos_{pthr_str_short}_psig']].apply(find_union, axis=1)

    df_hits_clust, distances, yvalues = clustering(
        df_struc, df_pvals,
        output_dir, input_gene,
        psig_columns=[f'SUM_union_neg_{pthr_str_short}_psig',
                    f'SUM_union_pos_{pthr_str_short}_psig'],
        pthr_cutoffs=[f'p<{pthr_str}', f'p<{pthr_str}'],
        screen_name='Meta', score_type='union',
        max_distances=25, merge_cols=['unipos', 'chain'],
    )

    # PLOTTING #
    plot_clustering(
        df_struc, df_pvals,
        df_hits_clust, clustering_radius,
        output_dir, input_gene,
        distances, yvalues,
        psig_columns=[f'SUM_union_neg_{pthr_str_short}_psig',
                    f'SUM_union_pos_{pthr_str_short}_psig'],
        names=['Negative', 'Positive'],
        pthr_cutoffs=[f'p<{pthr_str}', f'p<{pthr_str}'],
        screen_name='Meta', score_type='union',
        merge_col=['unipos', 'chain'],
    )

    # file_path = f'{output_dir}/cluster_union/plots/{input_gene}_Meta_union_Positive_Dendogram_6A.png'
    
    ## CHARACTERIZATION
    df_domains = pd.read_csv(f'{output_dir}/sequence_structure/{input_gene}_{input_uniprot}_domains.tsv', sep='\t')
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    for cutoff in [pthr]:
        for score_type in ['LFC', 'LFC3D']:

            # ENRICHMENT TEST #
            cutoff_str = str(cutoff).split('.')[1]
            df_meta = pd.read_csv(f'{output_dir}/meta-aggregate/{input_gene}_MetaAggr_{score_type}.tsv', sep='\t')

            input_df = pd.concat([df_meta, df_domains['Domain'], df_struc['pLDDT_dis']], axis=1)

            hit_columns = [f'SUM_{score_type}_neg_{cutoff_str}_p', f'SUM_{score_type}_pos_{cutoff_str}_p']
            input_df[hit_columns] = input_df[hit_columns].replace('-', np.nan).astype(float)

            results = enrichment_test(
                input_df,
                workdir=output_dir,
                input_gene=input_gene,
                hit_columns=hit_columns,
                hit_threshold=cutoff,
                feature_column='pLDDT_dis',
                feature_values=['confident', 'low', 'very low'],
                confidence_level=0.95,
            )

            plot_enrichment_test(
                enrichment_results=results,
                workdir=output_dir,
                input_gene=input_gene,
                hit_value=cutoff,
                feature_values=['confident', 'low', 'very low'],
            )

            os.rename(f'{output_dir}/characterization/{input_gene}_enrichment_test.pickle',
                    f'{output_dir}/characterization/{input_gene}_enrichment_test_{score_type}_{cutoff_str}.pickle')
            os.rename(f'{output_dir}/characterization/plots/{input_gene}_enrichment_test.png',
                    f'{output_dir}/characterization/plots/{input_gene}_enrichment_test_{score_type}_{cutoff_str}.png')

            # BARPLOTS #
            colnames = [f'SUM_{score_type}_neg_{cutoff_str}_psig', f'SUM_{score_type}_pos_{cutoff_str}_psig']
            input_df = pd.concat([df_struc, df_meta[colnames], df_domains['Domain']], axis=1)

            hits_feature_barplot(
                input_df,
                workdir=output_dir,
                input_gene=input_gene,
                category_col='pLDDT_dis',
                values_cols=colnames, values_vals=[f'p<0.{cutoff_str}', f'p<0.{cutoff_str}'], value_names=['NEG', 'POS'],
                plot_type='Count', colors = ['darkred', 'darkblue'],
            )

            os.rename(f'{output_dir}/characterization/plots/{input_gene}_Count_pLDDT_dis_barplot.png',
                    f'{output_dir}/characterization/plots/{input_gene}_Count_pLDDT_dis_barplot_{score_type}_{cutoff_str}.png')

        # SCATTERPLOT #
        df_lfc_dis = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_LFC_dis_wght.tsv", sep='\t')
        df_lfc3d_dis = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_LFC3D_dis_wght.tsv", sep='\t')

        df_dis_input = pd.DataFrame()
        df_dis_input['unipos'] = df_lfc_dis['unipos']
        for screen_name in screen_names:
            df_dis_input[f'Meta_LFC'] = df_lfc_dis[f'SUM_LFC']
            df_dis_input[f'Meta_LFC3D'] = df_lfc3d_dis[f'SUM_LFC3D']
            df_dis_input[f'Meta_LFC3D_dis'] = df_lfc3d_dis[f'SUM_LFC3D_dis']

        df_lfc = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC.tsv", sep='\t')
        df_lfc3d = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_MetaAggr_LFC3D.tsv", sep='\t')
        df_dis_input[f'SUM_LFC3D_neg_{cutoff_str}_psig'] = df_lfc3d[f'SUM_LFC3D_neg_{cutoff_str}_psig']
        df_dis_input[f'SUM_LFC3D_pos_{cutoff_str}_psig'] = df_lfc3d[f'SUM_LFC3D_pos_{cutoff_str}_psig']

        df_dis_input = df_dis_input.rename(columns={
            f'SUM_LFC3D_neg_{cutoff_str}_psig': f"Meta_LFC3D_neg_psig",
            f'SUM_LFC3D_pos_{cutoff_str}_psig': f"Meta_LFC3D_pos_psig",
            })

        lfc_lfc3d_scatter(
            df_input=df_dis_input,
            workdir=output_dir,
            input_gene=input_gene, screen_name='Meta',
            lfc3d_hit_threshold=cutoff,
        )
        os.rename(f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter.png',
                f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter_{cutoff_str}.png')

    # SCATTERPLOT #
    df_meta_wght = pd.read_csv(f"{output_dir}/meta-aggregate/{input_gene}_{score_type}_dis_wght.tsv", sep='\t')
    df_input = pd.concat([df_struc, df_meta_wght], axis=1)

    df_input[f'{score_type}_wght'] = df_input[f'SUM_{score_type}_wght'].abs() * 100
    df_input['direction'] = np.where(
        df_meta_wght[f'SUM_{score_type}_wght'].astype(float) > 0, 'POS', np.where(df_meta_wght[f'SUM_{score_type}_wght'].astype(float) < 0, 'NEG', 'ZERO')
    )
    df_input = df_input[df_input['direction'].isin(['NEG', 'POS'])]
    df_input = df_input[~df_input['bfactor_pLDDT'].isin(['-'])]
    df_input = df_input[~df_input['RSA'].isin(['-'])]
    df_input['bfactor_pLDDT'] = df_input['bfactor_pLDDT'].astype(float)
    df_input['RSA'] = df_input['RSA'].astype(float)

    pLDDT_RSA_scatter(
        df_input,
        workdir=output_dir,
        input_gene=input_gene,
        pLDDT_col='bfactor_pLDDT', RSA_col='RSA', size_col=f'{score_type}_wght', direction_col='direction',
        color_map = {'NEG': 'darkred', 'POS': 'darkblue'}
    )
if __name__ == '__main__':
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

    from beclust3d.lfc3d.structure import sequence_structural_features
    from beclust3d.lfc3d.preprocess_data import parse_be_data
    from beclust3d.lfc3d.preprocess_data_plot import plot_rawdata
    from beclust3d.qc.hypothesis_tests import hypothesis_test
    from beclust3d.lfc3d.randomize_data import randomize_data
    from beclust3d.lfc3d.prioritize_sequence import prioritize_by_sequence
    from beclust3d.lfc3d.prioritize_sequence_plot import plot_screendata_sequence
    from beclust3d.lfc3d.randomize_sequence import randomize_sequence
    from beclust3d.lfc3d.clustering import clustering
    from beclust3d.lfc3d.clustering_plot import plot_clustering
    from beclust3d.lfc3d.calculate_lfc3d import calculate_lfc3d
    from beclust3d.aggregate.metaaggregate import average_split_meta, bin_meta, znorm_meta
    from beclust3d.aggregate.aggregate_plot import average_split_bin_plots
    from beclust3d.lfc3d.characterization import enrichment_test
    from beclust3d.lfc3d.characterization_plot import plot_enrichment_test, hits_feature_barplot, lfc_lfc3d_scatter, pLDDT_RSA_scatter

    ## REQUIRED
    input_gene = 'MEN1'
    input_uniprot = 'O00255'
    input_chain = 'A' 
    workdir = './' 
    screen_dir = '../data/' 
    screens = 'PernerNature2023-MOLM13-Screen.tsv, PernerNature2023-MV411-Screen.tsv'
    mut_list_col = None
    mut_col = 'mut_type' 
    val_col = 'delta_beta_score' 
    gene_col = 'gene' 
    edits_col = 'predicted_edit' 
    gRNA_col = None
    output_dir= None

    # OPTIONAL
    user_fasta = './men1.fasta'
    user_pdb = './men1_AF3.pdb'
    user_dssp = None
    nRandom = 500 
    pthr = 0.001
    structure_radius = 6.0 
    clustering_radius = 6.0 
    
    main(input_gene=input_gene, input_uniprot=input_uniprot, input_chain=input_chain, workdir=workdir, screen_dir=screen_dir,\
        screens=screens, mut_list_col=mut_list_col, mut_col=mut_col, val_col=val_col, gene_col=gene_col, edits_col=edits_col,\
        gRNA_col=gRNA_col, output_dir=output_dir, user_fasta=user_fasta, user_pdb=user_pdb, user_dssp=user_dssp, nRandom=nRandom,\
        pthr=pthr, structure_radius=structure_radius, clustering_radius=clustering_radius
        )

"""
Step 3b (per-screen): Characterization of individual screens.
Usage example: python beclust3d_characterization_per_screen.py ./yaml/dnmt3a_local.yaml

Requires Steps 1 and 2 to have been run first.
Runs for any number of screens (>= 1).
"""

import os
import sys
import pandas as pd
import numpy as np

from be3d_helper import get_required, get_optional, load_config


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
    function_for_lfc        = kwargs['function_for_lfc']
    function_for_lfc3d      = kwargs['function_for_lfc3d']
    single_screen_pthr      = kwargs['single_screen_pthr']
    qa_passed_only          = kwargs['qa_passed_only']
    qa_controls             = kwargs['qa_controls']
    qa_cases                = kwargs['qa_cases']
    priority_on_alternative = kwargs['priority_on_alternative']

    if user_pdb:
        structureid = f'PDB-{input_uniprot}'
    else:
        structureid = f'AF-{input_uniprot}-F1-model_v6'

    single_screen_pthr_str = str(single_screen_pthr)

    print(f'All results will be saved in the following directory: {output_dir}')

    if isinstance(screens, str): screens = [screen.strip() for screen in screens.split(',')]
    screen_names = [s.split('.')[0] for s in screens]

    df_residuemap = pd.DataFrame()
    conserv_dfs = []
    gene_list = []
    if conservation_run:
        df_residuemap_filename = os.path.join(output_dir, f"conservation/{input_gene}{input_uniprot}_{alt_gene_name}{alt_uniprot_id}_residuemap_conservation.tsv")
        df_residuemap = pd.read_csv(df_residuemap_filename)

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

    if qa_passed_only:
        h2_ks_test_pd = pd.read_csv(f'{output_dir}/hypothesis_qa/KolmogorovSmirnov_hypothesis2.tsv', sep='\t')
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

    # CHARACTERIZATION (PER-SCREEN)
    df_domains = pd.read_csv(f'{output_dir}/sequence_structure/{input_gene}_{input_uniprot}_domains.tsv', sep='\t')
    df_struc   = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    # Map score_type -> the function name used to produce that score in Step 2

    for gene, screen_name in zip(gene_list, screen_names):
        cutoff     = single_screen_pthr
        cutoff_str = single_screen_pthr_str.split('.')[1]
        for score_type in ['LFC', 'LFC3D']:

            gene_type  = 'Alternative' if (conservation_run and screen_name.startswith(alt_screen_start)) else 'Original'
            df_screen  = pd.read_csv(f'{output_dir}/{score_type}/{gene_type}_{gene}_NonAggr_{score_type}.tsv', sep='\t')

            input_df = pd.concat([df_screen, df_domains['Domain'], df_struc['pLDDT_dis']], axis=1)
            hit_columns = [f'{screen_name}_{score_type}_neg_{cutoff_str}_p', f'{screen_name}_{score_type}_pos_{cutoff_str}_p']
            input_df[hit_columns] = input_df[hit_columns].replace('-', np.nan).astype(float)

            results = enrichment_test(
                input_df, workdir=output_dir, input_gene=input_gene,
                hit_columns=hit_columns, hit_threshold=cutoff,
                feature_column='pLDDT_dis', feature_values=['confident', 'low', 'very low'],
                confidence_level=0.95,
            )
            plot_enrichment_test(
                enrichment_results=results, workdir=output_dir, input_gene=input_gene,
                hit_value=cutoff, feature_values=['confident', 'low', 'very low'],
            )
            os.rename(f'{output_dir}/characterization/{input_gene}_enrichment_test.pickle',
                      f'{output_dir}/characterization/{input_gene}_enrichment_test_{score_type}_{cutoff_str}_{screen_name}.pickle')
            os.rename(f'{output_dir}/characterization/plots/{input_gene}_enrichment_test.png',
                      f'{output_dir}/characterization/plots/{input_gene}_enrichment_test_{score_type}_{cutoff_str}_{screen_name}.png')

            # Barplots
            colnames = [f'{screen_name}_{score_type}_neg_{cutoff_str}_psig', f'{screen_name}_{score_type}_pos_{cutoff_str}_psig']
            input_df2 = pd.concat([df_screen[colnames], df_struc, df_domains['Domain']], axis=1)

            if len(input_df2['Domain'].unique()) > 1 and all([len(input_df2[x].unique()) > 2 for x in colnames]): 
                hits_feature_barplot(
                    input_df2, workdir=output_dir, input_gene=input_gene,
                    category_col='Domain',
                    values_cols=colnames, values_vals=[f'p<0.{cutoff_str}', f'p<0.{cutoff_str}'],
                    value_names=['NEG', 'POS'],
                    plot_type='Count', 
                )
                os.rename(f'{output_dir}/characterization/plots/{input_gene}_Count_Domain_barplot.png',
                          f'{output_dir}/characterization/plots/{input_gene}_Count_Domain_barplot_{score_type}_{cutoff_str}_{screen_name}.png')

            if len(input_df2['pLDDT_dis'].unique()) > 1 and all([len(input_df2[x].unique()) > 2 for x in colnames]): 
                hits_feature_barplot(
                    input_df2, workdir=output_dir, input_gene=input_gene,
                    category_col='pLDDT_dis',
                    values_cols=colnames, values_vals=[f'p<0.{cutoff_str}', f'p<0.{cutoff_str}'],
                    value_names=['NEG', 'POS'],
                    plot_type='Count', 
                )
                os.rename(f'{output_dir}/characterization/plots/{input_gene}_Count_pLDDT_dis_barplot.png',
                          f'{output_dir}/characterization/plots/{input_gene}_Count_pLDDT_dis_barplot_{score_type}_{cutoff_str}_{screen_name}.png')

        # Scatter (per-screen)
        df_lfc_dis   = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_LFC_dis_wght.tsv", sep='\t')
        df_lfc3d_dis = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_LFC3D_dis_wght.tsv", sep='\t')
        df_lfc3d     = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv", sep='\t')

        df_dis_input = pd.DataFrame()
        df_dis_input['unipos']             = df_lfc_dis['unipos']
        df_dis_input[f'{screen_name}_LFC']      = df_lfc_dis[f'{screen_name}_LFC']
        df_dis_input[f'{screen_name}_LFC3D']    = df_lfc3d_dis[f'{screen_name}_LFC3D']
        df_dis_input[f'{screen_name}_LFC3D_dis'] = df_lfc3d_dis[f'{screen_name}_LFC3D_dis']
        df_dis_input[f'{screen_name}_LFC3D_neg_{cutoff_str}_psig'] = df_lfc3d[f'{screen_name}_LFC3D_neg_{cutoff_str}_psig']
        df_dis_input[f'{screen_name}_LFC3D_pos_{cutoff_str}_psig'] = df_lfc3d[f'{screen_name}_LFC3D_pos_{cutoff_str}_psig']
        df_dis_input = df_dis_input.rename(columns={
            f'{screen_name}_LFC3D_neg_{cutoff_str}_psig': f'{screen_name}_LFC3D_neg_psig',
            f'{screen_name}_LFC3D_pos_{cutoff_str}_psig': f'{screen_name}_LFC3D_pos_psig',
        })

        lfc_lfc3d_scatter(
            df_input=df_dis_input, workdir=output_dir,
            input_gene=input_gene, screen_name=screen_name, pthr=single_screen_pthr,
        )
        os.rename(f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter.png',
                  f'{output_dir}/characterization/plots/{input_gene}_LFC_LFC3D_scatter_{cutoff_str}_{screen_name}.png')

        # pLDDT/RSA scatter
        score_type   = 'LFC3D'  # default final score
        df_wght      = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_LFC3D_dis_wght.tsv", sep='\t')

        df_input = pd.concat([df_struc, df_wght], axis=1)
        df_input = df_input[df_input[f'{screen_name}_{score_type}'] != '-']
        df_input[f'{score_type}_wght'] = df_input[f'{screen_name}_{score_type}'].astype(float).abs() * 100
        df_input['direction'] = np.where(
            df_input[f'{screen_name}_{score_type}'].astype(float) > 0, 'POS',
            np.where(df_input[f'{screen_name}_{score_type}'].astype(float) < 0, 'NEG', 'ZERO')
        )
        df_input = df_input[df_input['direction'].isin(['NEG', 'POS'])]
        df_input = df_input[~df_input['bfactor_pLDDT'].isin(['-'])]
        df_input = df_input[~df_input['RSA'].isin(['-'])]
        df_input['bfactor_pLDDT'] = df_input['bfactor_pLDDT'].astype(float)
        df_input['RSA'] = df_input['RSA'].astype(float)

        pLDDT_RSA_scatter(
            df_input, workdir=output_dir, input_gene=input_gene,
            pLDDT_col='bfactor_pLDDT', RSA_col='RSA', size_col=f'{score_type}_wght', direction_col='direction',
            color_map={'NEG': 'darkred', 'POS': 'darkblue'}, 
        )

        os.rename(f'{output_dir}/characterization/plots/{input_gene}_pLDDT_RSA_scatter.png',
                  f'{output_dir}/characterization/plots/{input_gene}_pLDDT_RSA_scatter_{screen_name}.png')

    print(f'Output from beclust3d.lfc3d.characterization is saved in the following directory: {output_dir}characterization/')

    print("Step 3b complete: per-screen characterization done")


if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    beclust3d_path = get_optional(config, 'beclust3d_path', str, '.')
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

    from beclust3d.lfc3d.characterization import enrichment_test
    from beclust3d.lfc3d.characterization_plot import plot_enrichment_test, hits_feature_barplot, lfc_lfc3d_scatter, pLDDT_RSA_scatter
    from beclust3d.helpers.visualization.g2p import g2p_formatted_hit_cluster

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
    function_for_lfc        = get_optional(config, 'function_for_lfc',        str,                'mean')
    function_for_lfc3d      = get_optional(config, 'function_for_lfc3d',      str,                'mean')
    single_screen_pthr      = get_optional(config, 'pthr.single_screen',      (int, float),       0.05)
    qa_passed_only          = get_optional(config, 'qa.qa_passed_only',       bool,               False)
    qa_controls             = get_optional(config, 'qa.controls',             (list, type(None)), ['No Mutation'])
    qa_cases                = get_optional(config, 'qa.cases',                (list, type(None)), ['Nonsense', 'Splice'])
    priority_on_alternative = get_optional(config, 'priority_on_alternative', bool,               False)

    main(
        input_gene=input_gene, input_uniprot=input_uniprot,
        screens=screens, output_dir=output_dir,
        user_pdb=user_pdb,
        single_screen_pthr=single_screen_pthr,
        function_for_lfc=function_for_lfc, function_for_lfc3d=function_for_lfc3d,
        conservation_run=conservation_run, alt_gene_name=alt_gene_name,
        alt_uniprot_id=alt_uniprot_id, alt_screen_start=alt_screen_start,
        qa_passed_only=qa_passed_only, 
        qa_controls=qa_controls, qa_cases=qa_cases,
        priority_on_alternative=priority_on_alternative,
    )

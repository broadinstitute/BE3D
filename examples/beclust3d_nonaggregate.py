"""
Step 2b: Nonaggregate scoring, per-screen clustering, and characterization.
Usage example: python be3d_aggregate_per_screen.py ./yaml/dnmt3a_local.yaml

Requires Step 2a (be3d_calculate_lfc3d.py) to have been run first.
"""

import os
import sys
import pandas as pd

from be3d_helper import get_required, get_optional, load_config


def find_union(input, pthr_str):
    if input[0] == f'p<{pthr_str}' or input[1] == f'p<{pthr_str}':
        return f'p<{pthr_str}'
    elif input[0] == f'p>={pthr_str}' or input[1] == f'p>={pthr_str}':
        return f'p>={pthr_str}'
    else:
        return '-'

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
    single_screen_pthr      = kwargs['single_screen_pthr']
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

    single_screen_pthr_str = str(single_screen_pthr)
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
        print("qa_only=True: skipping aggregate per-screen step.")
        sys.exit()

    if qa_passed_only:
        h2_ks_test_pd = pd.read_csv(f'{output_dir}/hypothesis_qc/KolmogorovSmirnov_hypothesis2.tsv', sep='\t')
        h2_ks_test_pd = h2_ks_test_pd.replace(-999, None)
        white_screen_list = h2_ks_test_pd[
            (h2_ks_test_pd[f"p_{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}"] < 0.05) &
            (h2_ks_test_pd['gene_name'].isin(gene_list))
        ]['screenid'].to_list()
        print(f'original screen size: {len(screen_names)}, screen white list size: {len(white_screen_list)}, '
              f'QA-passed screen size: {len(list(set(screen_names).intersection(white_screen_list)))}')
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

    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')

    def run_aggregate_per_species(gene, screen_names, gene_type, pdb_file):
        df_LFC_LFC3D = pd.read_csv(f'{output_dir}LFC3D/{gene_type}_{gene}_LFC_LFC3D_LFC3Dr.tsv.gz', sep='\t')

        # Per-screen clustering
        for score_type in ['LFC', 'LFC3D']:
            df_bidir = average_split_score(df_LFC_LFC3D, output_dir, gene, screen_names, score_type=score_type, gene_type=gene_type)
            df_dis, _, _ = bin_score(df_bidir, output_dir, gene, screen_names, score_type=score_type, gene_type=gene_type)
            znorm_score(df_bidir, output_dir, gene, screen_names, pthrs=[0.05, 0.01, 0.001], score_type=score_type, gene_type=gene_type)

            df_scoretype = pd.read_csv(f'{output_dir}{score_type}/{gene_type}_{gene}_NonAggr_{score_type}.tsv', sep='\t')
            for screen_name in screen_names:
                average_split_bin_plots(
                    df_scoretype,
                    workdir=output_dir,
                    input_gene=gene,
                    screen_name=screen_name,
                    func='',
                    pthr=single_screen_pthr,
                    score_type=score_type,
                    aggregate_dir=score_type,
                    save_type='svg',
                )
            print(f'Output from beclust3d.aggregate.aggregate_plot is saved in the following directory: {output_dir}{score_type}/')

            df_pvals = pd.read_csv(f'{output_dir}{score_type}/{gene_type}_{gene}_NonAggr_{score_type}.tsv', sep='\t')

            for screen_name in screen_names:
                pref = f'{screen_name}_{score_type}'
                df_hits_clust, distances, yvalues = clustering(
                    df_struc, df_pvals,
                    output_dir, gene,
                    psig_columns=[f'{pref}_neg_05_psig', f'{pref}_pos_05_psig', f'{pref}_neg_01_psig',
                                  f'{pref}_pos_01_psig', f'{pref}_neg_001_psig', f'{pref}_pos_001_psig'],
                    pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
                    screen_name=screen_name, score_type=score_type,
                    max_distances=20, merge_cols=['unipos', 'chain'],
                    atom_level=pdb_file if atom_level_naa == True else False,
                )

                plot_clustering(
                    df_struc, df_pvals,
                    df_hits_clust, clustering_radius,
                    output_dir, gene,
                    distances, yvalues,
                    names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
                    psig_columns=[f'{pref}_neg_05_psig', f'{pref}_pos_05_psig', f'{pref}_neg_01_psig',
                                  f'{pref}_pos_01_psig', f'{pref}_neg_001_psig', f'{pref}_pos_001_psig'],
                    pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
                    screen_name=screen_name, score_type=score_type,
                    merge_col=['unipos', 'chain'],
                    save_type='svg',
                    dendrogram_subplots_kwargs={'figsize': (15, 3.5)},
                )
            print(f'Output from beclust3d.lfc3d.clustering is saved in the following directory: {output_dir}cluster_{score_type}/')

        # Per-screen scatter and union clustering
        cutoff_str = str(single_screen_pthr).split('.')[1]
        df_lfc = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_NonAggr_LFC.tsv", sep='\t')
        df_lfc_dis = pd.read_csv(f"{output_dir}/LFC/{gene_type}_{gene}_LFC_dis_wght.tsv", sep='\t')
        df_lfc3d = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv", sep='\t')
        df_lfc3d_dis = pd.read_csv(f"{output_dir}/LFC3D/{gene_type}_{gene}_LFC3D_dis_wght.tsv", sep='\t')

        for screen_name in screen_names:

            df_dis_input = pd.DataFrame()
            df_dis_input['unipos'] = df_lfc_dis['unipos']
            df_dis_input = pd.concat([
                df_dis_input,
                df_lfc_dis.filter(regex=r'LFC$'),
                df_lfc3d_dis.filter(regex=r'LFC3D$'),
                df_lfc3d.filter(regex=fr'LFC3D_dis$|{cutoff_str}_psig$')
            ])
            df_dis_input = df_dis_input.rename({
                f'{screen_name}_LFC3D_neg_{cutoff_str}_psig': f'{screen_name}_LFC3D_neg_psig',
                f'{screen_name}_LFC3D_pos_{cutoff_str}_psig': f'{screen_name}_LFC3D_pos_psig',
            }, axis=1)
            lfc_lfc3d_scatter(
                df_input=df_dis_input,
                workdir=output_dir,
                input_gene=gene, screen_name=screen_name,
                pthr=single_screen_pthr,
            )
            os.rename(f'{output_dir}/characterization/plots/{gene}_LFC_LFC3D_scatter.png',
                      f'{output_dir}/characterization/plots/{gene}_LFC_LFC3D_scatter_{cutoff_str}_{screen_name}.png')

        # Union clustering
        df_pvals_LFC3D = pd.read_csv(f'{output_dir}/LFC3D/{gene_type}_{gene}_NonAggr_LFC3D.tsv', sep='\t')
        df_pvals_LFC = pd.read_csv(f'{output_dir}/LFC/{gene_type}_{gene}_NonAggr_LFC.tsv', sep='\t')
        df_pvals = pd.concat([df_pvals_LFC3D, df_pvals_LFC.drop(columns=['unipos', 'unires', 'chain'])], axis=1)

        for screen_name in screen_names:
            for each_pthr in ['05', '01', '001']:
                df_pvals[f'{screen_name}_union_neg_{each_pthr}_psig'] = df_pvals[
                    [f'{screen_name}_LFC_neg_{each_pthr}_psig', f'{screen_name}_LFC3D_neg_{each_pthr}_psig']
                ].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)
                df_pvals[f'{screen_name}_union_pos_{each_pthr}_psig'] = df_pvals[
                    [f'{screen_name}_LFC_pos_{each_pthr}_psig', f'{screen_name}_LFC3D_pos_{each_pthr}_psig']
                ].apply(lambda row: find_union(row, f'0.{each_pthr}'), axis=1)

            df_hits_clust, distances, yvalues = clustering(
                df_struc, df_pvals,
                output_dir, gene,
                psig_columns=[f'{screen_name}_union_neg_05_psig', f'{screen_name}_union_pos_05_psig',
                               f'{screen_name}_union_neg_01_psig', f'{screen_name}_union_pos_01_psig',
                               f'{screen_name}_union_neg_001_psig', f'{screen_name}_union_pos_001_psig'],
                pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
                screen_name=screen_name, score_type='union',
                max_distances=20, merge_cols=['unipos', 'chain'],
                atom_level=pdb_file if atom_level_naa == True else False, 
            )

            plot_clustering(
                df_struc, df_pvals,
                df_hits_clust, clustering_radius,
                output_dir, gene,
                distances, yvalues,
                names=['Negative', 'Positive', 'Negative', 'Positive', 'Negative', 'Positive'],
                psig_columns=[f'{screen_name}_union_neg_05_psig', f'{screen_name}_union_pos_05_psig',
                               f'{screen_name}_union_neg_01_psig', f'{screen_name}_union_pos_01_psig',
                               f'{screen_name}_union_neg_001_psig', f'{screen_name}_union_pos_001_psig'],
                pthr_cutoffs=['p<0.05', 'p<0.05', 'p<0.01', 'p<0.01', 'p<0.001', 'p<0.001'],
                screen_name=screen_name, score_type='union',
                merge_col=['unipos', 'chain'],
                save_type='svg',
                dendrogram_subplots_kwargs={'figsize': (15, 3.5)}, 
            )
        print(f'Output from beclust3d.lfc3d.clustering is saved in the following directory: {output_dir}cluster_union/')

    if conservation_run:
        run_aggregate_per_species(input_gene,   [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is None],     gene_type='Original',    pdb_file=pdb_file)
        run_aggregate_per_species(alt_gene_name, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is not None], gene_type='Alternative', pdb_file=pdb_file)
        merged_df_LFC_LFC3D = pd.read_csv(f'{output_dir}/LFC3D/Original_{input_gene}_LFC_LFC3D_LFC3Dr.tsv.gz', sep='\t').merge(
            pd.read_csv(f'{output_dir}/LFC3D/Alternative_{alt_gene_name}_LFC_LFC3D_LFC3Dr.tsv.gz', sep='\t'),
            on='unipos', how='inner', suffixes=('', f'_{alt_gene_name}')
        ).sort_index()
    else:
        run_aggregate_per_species(input_gene, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is None], gene_type='Original', pdb_file=pdb_file)

    def merge_two_tsvs(gene_name, results_dir, file_pattern, original_gene, alternative_gene, priority_on_alternative=False, compression=False):
        ext = '.tsv.gz' if compression else '.tsv'
        original_tsv = os.path.join(results_dir, f'Original_{original_gene}_{file_pattern}{ext}')
        alternative_tsv = os.path.join(results_dir, f'Alternative_{alternative_gene}_{file_pattern}{ext}')

        if os.path.exists(original_tsv) and os.path.exists(alternative_tsv):
            df1 = pd.read_csv(original_tsv, sep='\t')
            df2 = pd.read_csv(alternative_tsv, sep='\t')
            if not priority_on_alternative:
                common_cols = list(set(df1.columns) & set(df2.columns))
                merged = pd.merge(df1, df2, on=common_cols, how="outer")
            else:
                merged = pd.merge(df2, df1, how="outer")
            merged.to_csv(os.path.join(results_dir, f'{gene_name}_{file_pattern}{ext}'), sep="\t", index=False, compression='gzip' if compression else None)
        elif os.path.exists(original_tsv):
            df1 = pd.read_csv(original_tsv, sep='\t')
            df1.to_csv(os.path.join(results_dir, f'{gene_name}_{file_pattern}{ext}'), sep="\t", index=False, compression='gzip' if compression else None)
        else:
            df2 = pd.read_csv(alternative_tsv, sep='\t')
            df2.to_csv(os.path.join(results_dir, f'{gene_name}_{file_pattern}{ext}'), sep="\t", index=False, compression='gzip' if compression else None)

    merge_two_tsvs(input_gene, f'{output_dir}/LFC',  'LFC_bidirectional',  input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene, f'{output_dir}/LFC',  'LFC_dis_wght',       input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene, f'{output_dir}/LFC',  'NonAggr_LFC',        input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene, f'{output_dir}/LFC3D','LFC3D_bidirectional', input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene, f'{output_dir}/LFC3D','LFC3D_dis_wght',      input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene, f'{output_dir}/LFC3D','NonAggr_LFC3D',       input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative)
    merge_two_tsvs(input_gene, f'{output_dir}/LFC3D','LFC_LFC3D_LFC3Dr',    input_gene, alt_gene_name, priority_on_alternative=priority_on_alternative, compression=True)

    print("Step 2b complete: nonaggregate scoring, per-screen clustering and characterization")


if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    beclust3d_path = get_optional(config, 'beclust3d_path', str, '.')
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

    from beclust3d.lfc3d.clustering import clustering
    from beclust3d.lfc3d.clustering_plot import plot_clustering
    from beclust3d.lfc3d.characterization_plot import lfc_lfc3d_scatter
    from beclust3d.aggregate.nonaggregate import average_split_score, bin_score, znorm_score
    from beclust3d.aggregate.aggregate_plot import average_split_bin_plots

    # REQUIRED
    input_gene    = get_required(config, 'input_gene',    str)
    input_uniprot = get_required(config, 'input_uniprot', str)
    screens       = get_required(config, 'screens',       (str, list))
    output_dir    = get_required(config, 'output_dir',    str)

    conservation_run  = get_required(config, 'conservation.run',               bool)
    alt_gene_name     = get_required(config, 'conservation.alt_gene_name',     (str, type(None)))
    alt_uniprot_id    = get_required(config, 'conservation.alt_uniprot_id',    (str, type(None)))
    alt_screen_start  = get_required(config, 'conservation.alt_screen_start',  (str, type(None)))

    # OPTIONAL
    user_pdb                = get_optional(config, 'user_pdb',                (str, type(None)),  None)
    single_screen_pthr      = get_optional(config, 'pthr.single_screen',      (int, float),       0.05)
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
        user_pdb=user_pdb, single_screen_pthr=single_screen_pthr,
        clustering_radius=clustering_radius,
        conservation_run=conservation_run, alt_gene_name=alt_gene_name,
        alt_uniprot_id=alt_uniprot_id, alt_screen_start=alt_screen_start,
        qa_passed_only=qa_passed_only, qa_only=qa_only,
        qa_controls=qa_controls, qa_cases=qa_cases,
        priority_on_alternative=priority_on_alternative,
        atom_level_naa=atom_level_naa,
    )

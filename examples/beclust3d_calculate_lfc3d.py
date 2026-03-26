"""
Step 2a: Parse/preprocess BE data, prioritize, and LFC3D calculation.
Usage example: python be3d_calculate_lfc3d.py ./yaml/dnmt3a_local.yaml

Requires Step 1 to have been run first (sequence/structural features + hypothesis test outputs).
"""

import os
import sys
import pandas as pd

from be3d_helper import get_required, get_optional, load_config


def main(**kwargs):
    # REQUIRED
    input_gene       = kwargs['input_gene']
    input_uniprot    = kwargs['input_uniprot']
    input_chain      = kwargs['input_chain']
    screen_dir       = kwargs['screen_dir']
    screens          = kwargs['screens']
    mut_col          = kwargs['mut_col']
    val_col          = kwargs['val_col']
    gene_col         = kwargs['gene_col']
    edits_col        = kwargs['edits_col']
    output_dir       = kwargs['output_dir']
    mut_categories   = kwargs['mut_categories']
    mut_delimiter    = kwargs['mut_delimiter']
    conservation_run = kwargs['conservation_run']
    alt_gene_name    = kwargs['alt_gene_name']
    alt_uniprot_id   = kwargs['alt_uniprot_id']
    alt_screen_start = kwargs['alt_screen_start']

    # OPTIONAL
    user_pdb                = kwargs['user_pdb']
    function_for_lfc        = kwargs['function_for_lfc']
    function_for_lfc3d      = kwargs['function_for_lfc3d']
    nRandom                 = kwargs['nRandom']
    single_screen_pthr      = kwargs['single_screen_pthr']
    qa_passed_only          = kwargs['qa_passed_only']
    qa_controls             = kwargs['qa_controls']
    qa_cases                = kwargs['qa_cases']
    priority_on_alternative = kwargs['priority_on_alternative']
    ppi_chain_gene_dict     = kwargs['ppi_chain_gene_dict']
    ppi_gene_edits_dict     = kwargs['ppi_gene_edits_dict']
    v_score_threshold       = kwargs['v_score_threshold']
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
    input_dfs = [pd.read_csv(os.path.join(screen_dir, s), sep='\t') for s in screens]

    df_residuemap = pd.DataFrame()
    conserv_dfs = []
    gene_list = []
    # If considering conservation, conserv_dfs/gene_list need to be different
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

    if qa_passed_only:
        h2_ks_test_pd = pd.read_csv(f'{output_dir}/hypothesis_qc/KolmogorovSmirnov_hypothesis2.tsv', sep='\t')
        h2_ks_test_pd = h2_ks_test_pd.replace(-999, None)
        white_screen_list = h2_ks_test_pd[
            (h2_ks_test_pd[f"p_{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}"] < 0.05) &
            (h2_ks_test_pd['gene_name'].isin(gene_list))
        ]['screenid'].to_list()
        print(f'original screen size: {len(screen_names)}, screen white list size: {len(white_screen_list)}, '
              f'QA-passed screen size: {len(list(set(screen_names).intersection(white_screen_list)))}')
        input_dfs = [pd.read_csv(os.path.join(screen_dir, f'{s}.tsv'), sep='\t') for s in screen_names]
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

    # Parse BE data
    parse_be_data(
        output_dir,
        input_dfs, input_gene, screen_names,
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        edits_col=edits_col,
        mut_categories=mut_categories,
        mut_delimiter=mut_delimiter,
        conserv_dfs=conserv_dfs,
        conserv_col='alternative_res_pos',
        gene_list=gene_list,
        v_score_threshold=v_score_threshold,
    )
    print(f'Output from beclust3d.lfc3d.preprocess_data is saved in the following directory: {output_dir}screendata/')

    # Plot raw data
    plot_rawdata(
        output_dir,
        input_dfs,
        screen_names,
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        mut_categories=mut_categories,
        save_type='svg',
    )

    df_missense_list = [
        pd.read_csv(f'{output_dir}/screendata/{gene}_{screen_name}_Missense.tsv', sep='\t')
        for gene, screen_name in zip(gene_list, screen_names)
    ]
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
    sanitary_check(df_struc, df_missense_list)

    # Randomize data
    for df_missense, screen_name, gene in zip(df_missense_list, screen_names, gene_list):
        randomize_data(
            df_missense,
            output_dir, gene,
            screen_name,
            nRandom=nRandom,
            seed=True,
        )
    print(f'Output from beclust3d.lfc3d.randomize_data is saved in the following directory: {output_dir}screendata_rand/')

    # Prioritize by sequence
    for gene, screen_name, df_consrv in zip(gene_list, screen_names, conserv_dfs):
        df_control = pd.read_csv(f'{output_dir}/screendata/{gene}_{screen_name}_No_Mutation.tsv', sep='\t', index_col=0)
        df_dict = {}

        for mut in ['Missense', 'Silent', 'Nonsense']:
            filepath = f'{output_dir}/screendata/{gene}_{screen_name}_{mut}.tsv'
            if os.path.exists(filepath):
                df_dict[mut] = pd.read_csv(filepath, sep='\t', index_col=0)

        if df_consrv is not None:
            target_res_pos, target_res = 'original_res_pos', 'unires'
            alternate_res_pos, alternate_res = 'alternative_res_pos', 'alternative_res'
            df_missense = prioritize_by_sequence(
                df_dict,
                df_struc, df_consrv, df_control,
                output_dir,
                gene, screen_name,
                target_res_pos=target_res_pos,
                alt_res_pos=alternate_res_pos,
                alt_res=alternate_res,
            )
        else:
            df_missense = prioritize_by_sequence(
                df_dict,
                df_struc, df_consrv, df_control,
                output_dir,
                gene, screen_name,
            )

        df_rand = pd.read_csv(f'{output_dir}/screendata_rand/{gene}_{screen_name}_Missense_rand.tsv.gz', sep='\t')
        randomize_sequence(
            df_missense, df_rand,
            output_dir,
            gene, screen_name,
            nRandom=nRandom, conservation=False,
            muttype='Missense',
            function_name=function_for_lfc,
            target_pos='unipos', target_res=None,
        )

        plot_screendata_sequence(
            df_missense,
            output_dir,
            gene, screen_name, function_name=function_for_lfc, muttype='Missense',
        )

    print(f'Output from beclust3d.lfc3d.prioritize_sequence is saved in the following directory: {output_dir}screendata_sequence/')
    print(f'Output from beclust3d.lfc3d.randomize_sequence is saved in the following directory: {output_dir}screendata_sequence_rand/')

    # Calculate LFC3D for original gene
    df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
    df_edits_list = []
    df_rand_list = []

    for screen_name in [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is None]:
        df_missense = pd.read_csv(f'{output_dir}/screendata_sequence/{input_gene}_{screen_name}_protein_edits.tsv', sep='\t')
        df_missense[f'{function_for_lfc}_Missense_LFC_plab_input'] = df_missense[f'{function_for_lfc}_Missense_LFC_p'].apply(
            lambda x: f'p<{single_screen_pthr_str}' if x < single_screen_pthr else f'p>={single_screen_pthr_str}'
        )

        clustering(
            df_struc, df_missense,
            output_dir, input_gene,
            psig_columns=[f'{function_for_lfc}_Missense_LFC_plab_input'],
            pthr_cutoffs=[f'p<{single_screen_pthr_str}'],
            screen_name=screen_name, score_type='LFC',
            max_distances=20, merge_cols=['unipos', 'chain'],
            atom_level=pdb_file if atom_level_naa == True else False,
        )

        df_protein_edits = pd.read_csv(f'{output_dir}/screendata_sequence/{input_gene}_{screen_name}_protein_edits.tsv', sep='\t')
        df_edits_list.append(df_protein_edits)
        df_protein_edits_rand = pd.read_csv(f'{output_dir}/screendata_sequence_rand/{input_gene}_{screen_name}_Missense_protein_edits_rand.tsv.gz', sep='\t')
        df_rand_list.append(df_protein_edits_rand)

    df_LFC_LFC3D_original = calculate_lfc3d(
        df_struc, df_edits_list, df_rand_list,
        output_dir, input_gene, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is None],
        nRandom=nRandom, muttype='Missense',
        function_type_lfc=function_for_lfc,
        function_type_lfc3d=function_for_lfc3d,
        conserved_only=False, gene_type='Original',
        target_gene_chain=input_chain,
        ppi_chain_gene_dict=ppi_chain_gene_dict,
        ppi_gene_edits_dict=ppi_gene_edits_dict,
    )

    if conservation_run:
        # Calculate LFC3D for alternate gene
        df_struc = pd.read_csv(f'{output_dir}/sequence_structure/{structureid}_coord_struc_features.tsv', sep='\t')
        df_edits_list = []
        df_rand_list = []

        for screen_name in [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is not None]: 
            df_missense = pd.read_csv(f'{output_dir}/screendata_sequence/{alt_gene_name}_{screen_name}_protein_edits.tsv', sep='\t')
            df_missense[f'{function_for_lfc}_Missense_LFC_plab_input'] = df_missense[f'{function_for_lfc}_Missense_LFC_p'].apply(
                lambda x: f'p<{single_screen_pthr_str}' if x < single_screen_pthr else f'p>={single_screen_pthr_str}'
            )

            clustering(
                df_struc, df_missense,
                output_dir, alt_gene_name,
                psig_columns=[f'{function_for_lfc}_Missense_LFC_plab_input'],
                pthr_cutoffs=[f'p<{single_screen_pthr_str}'],
                screen_name=screen_name, score_type='LFC',
                max_distances=20, merge_cols=['unipos', 'chain'],
                atom_level=pdb_file if atom_level_naa == True else False,
            )

            df_protein_edits = pd.read_csv(f'{output_dir}/screendata_sequence/{alt_gene_name}_{screen_name}_protein_edits.tsv', sep='\t')
            df_edits_list.append(df_protein_edits)
            df_protein_edits_rand = pd.read_csv(f'{output_dir}/screendata_sequence_rand/{alt_gene_name}_{screen_name}_Missense_protein_edits_rand.tsv.gz', sep='\t')
            df_rand_list.append(df_protein_edits_rand)

        df_LFC_LFC3D_alt = calculate_lfc3d(
            df_struc, df_edits_list, df_rand_list,
            output_dir, alt_gene_name, [x[1] for x in zip(conserv_dfs, screen_names) if x[0] is not None], 
            nRandom=nRandom, muttype='Missense',
            function_type_lfc=function_for_lfc,
            function_type_lfc3d=function_for_lfc3d,
            conserved_only=False, gene_type='Alternative',
            target_gene_chain=input_chain,
            ppi_chain_gene_dict=ppi_chain_gene_dict,
            ppi_gene_edits_dict=ppi_gene_edits_dict,
        )

    print(f'Output from beclust3d.lfc3d.calculate_lfc3d is saved in the following directory: {output_dir}LFC3D/')

    print("Step 2a complete: LFC3D calculated")


if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    beclust3d_path = get_optional(config, 'beclust3d_path', str, '.')
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

    from beclust3d.lfc3d.preprocess_data import parse_be_data, sanitary_check
    from beclust3d.lfc3d.preprocess_data_plot import plot_rawdata
    from beclust3d.lfc3d.randomize_data import randomize_data
    from beclust3d.lfc3d.prioritize_sequence import prioritize_by_sequence
    from beclust3d.lfc3d.prioritize_sequence_plot import plot_screendata_sequence
    from beclust3d.lfc3d.randomize_sequence import randomize_sequence
    from beclust3d.lfc3d.clustering import clustering
    from beclust3d.lfc3d.calculate_lfc3d import calculate_lfc3d

    # REQUIRED
    input_gene    = get_required(config, 'input_gene',    str)
    input_uniprot = get_required(config, 'input_uniprot', str)
    input_chain   = get_required(config, 'input_chain',   str)
    screen_dir    = get_required(config, 'screen_dir',    str)
    screens       = get_required(config, 'screens',       (str, list))
    output_dir    = get_required(config, 'output_dir',    str)

    conservation_run  = get_required(config, 'conservation.run',               bool)
    v_score_threshold = get_required(config, 'conservation.v_score_threshold', (int, float))
    alt_gene_name     = get_required(config, 'conservation.alt_gene_name',     (str, type(None)))
    alt_uniprot_id    = get_required(config, 'conservation.alt_uniprot_id',    (str, type(None)))
    alt_screen_start  = get_required(config, 'conservation.alt_screen_start',  (str, type(None)))

    mut_col       = get_required(config, 'database.mut_col',       str)
    val_col       = get_required(config, 'database.val_col',       str)
    gene_col      = get_required(config, 'database.gene_col',      str)
    edits_col     = get_required(config, 'database.edits_col',     str)
    mut_delimiter = get_required(config, 'database.mut_delimiter', str)

    mut_categories = []
    mut_categories.extend(get_required(config, 'mutation_category.nonsense',    list))
    mut_categories.extend(get_required(config, 'mutation_category.splice',      list))
    mut_categories.extend(get_required(config, 'mutation_category.missense',    list))
    mut_categories.extend(get_required(config, 'mutation_category.silent',      list))
    mut_categories.extend(get_required(config, 'mutation_category.no_mutation', list))
    mut_categories.extend(get_required(config, 'mutation_category.intron',      list))

    # OPTIONAL
    user_pdb                = get_optional(config, 'user_pdb',                (str, type(None)),  None)
    function_for_lfc        = get_optional(config, 'function_for_lfc',        str,                'mean')
    function_for_lfc3d      = get_optional(config, 'function_for_lfc3d',      str,                'mean')
    nRandom                 = get_optional(config, 'nRandom',                 int,                500)
    single_screen_pthr      = get_optional(config, 'pthr.single_screen',      (int, float),       0.05)
    qa_passed_only          = get_optional(config, 'qa.qa_passed_only',       bool,               False)
    qa_controls             = get_optional(config, 'qa.controls',             (list, type(None)), ['No Mutation'])
    qa_cases                = get_optional(config, 'qa.cases',                (list, type(None)), ['Nonsense', 'Splice'])
    priority_on_alternative = get_optional(config, 'priority_on_alternative', bool,               False)
    ppi_chain_gene_dict     = get_optional(config, 'ppi_chain_gene_dict',     (dict, type(None)), None)
    ppi_gene_edits_dict     = get_optional(config, 'ppi_gene_edits_dict',     (dict, type(None)), None)
    atom_level_naa          = get_optional(config, 'atom_level_naa',          bool,               False)
    muscle_path             = get_optional(config, 'muscle_path',             (str, type(None)),  'muscle')

    main(
        input_gene=input_gene, input_uniprot=input_uniprot, input_chain=input_chain,
        screen_dir=screen_dir, screens=screens, mut_col=mut_col, val_col=val_col,
        gene_col=gene_col, edits_col=edits_col, output_dir=output_dir,
        user_pdb=user_pdb, nRandom=nRandom, single_screen_pthr=single_screen_pthr,
        function_for_lfc=function_for_lfc, function_for_lfc3d=function_for_lfc3d,
        mut_categories=mut_categories, mut_delimiter=mut_delimiter,
        conservation_run=conservation_run, alt_gene_name=alt_gene_name,
        alt_uniprot_id=alt_uniprot_id, alt_screen_start=alt_screen_start,
        v_score_threshold=v_score_threshold,
        qa_passed_only=qa_passed_only, 
        qa_controls=qa_controls, qa_cases=qa_cases,
        priority_on_alternative=priority_on_alternative,
        ppi_chain_gene_dict=ppi_chain_gene_dict,
        ppi_gene_edits_dict=ppi_gene_edits_dict,
        atom_level_naa=atom_level_naa, muscle_path=muscle_path,
    )

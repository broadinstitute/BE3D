"""
Step 1: Extract sequence/structural features and BE-QA (hypothesis testing)

Usage example: python be3d_step1_hypothesis.py ./yaml/dnmt3a_local.yaml
"""

import os
import sys
import pandas as pd
import shutil

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
    output_dir       = kwargs['output_dir']
    conservation_run = kwargs['conservation_run']
    alt_gene_name    = kwargs['alt_gene_name']
    alt_uniprot_id   = kwargs['alt_uniprot_id']
    alt_screen_start = kwargs['alt_screen_start']
    config_yaml      = kwargs['config_yaml']
    # OPTIONAL
    user_fasta              = kwargs['user_fasta']
    user_pdb                = kwargs['user_pdb']
    user_dssp               = kwargs['user_dssp']
    structure_radius        = kwargs['structure_radius']
    qa_controls             = kwargs['qa_controls']
    qa_cases                = kwargs['qa_cases']
    atom_level_naa          = kwargs['atom_level_naa']
    muscle_path             = kwargs['muscle_path']
    priority_on_alternative = kwargs['priority_on_alternative']

    # Run BE-QA
    if user_pdb:
        structureid = f'PDB-{input_uniprot}'
    else:
        structureid = f'AF-{input_uniprot}-F1-model_v6'

    os.makedirs(output_dir, exist_ok=True)
    try:
        shutil.copy2(config_yaml, os.path.join(output_dir, os.path.basename(config_yaml)))
    except shutil.SameFileError:
        print(f"{config_yaml} and {os.path.join(output_dir, os.path.basename(config_yaml))} are the same file.")
    print(f'All results will be saved in the following directory: {output_dir}')

    if isinstance(screens, str): screens = [screen.strip() for screen in screens.split(',')]
    screen_names = [s.split('.')[0] for s in screens]
    input_dfs = [pd.read_csv(os.path.join(screen_dir, s), sep='\t') for s in screens]

    conserv_dfs = []
    gene_list = []
    # If considering conservation, conserv_dfs/gene_list need to be different
    if conservation_run:
        _, df_residuemap = conservation(
            output_dir,
            input_gene, alt_gene_name,
            input_uniprot, alt_uniprot_id,
            muscle_path=muscle_path,
        )
        print(f'Output from beclust3d.lfc3d.conservation is saved in the following directory: {output_dir}conservation/')

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

    # Sequence / structural features (original gene only)
    sequence_structural_features(
        output_dir,
        input_gene, input_uniprot, structureid,
        user_fasta=user_fasta,
        user_pdb=user_pdb,
        user_dssp=user_dssp,
        target_chainid=input_chain,
        radius=structure_radius,
        atom_level_naa=atom_level_naa, 
    )
    print(f'Output from beclust3d.lfc3d.structure is saved in the following directory: {output_dir}sequence_structure/')

    # Hypothesis testing (BE-QA)
    hypothesis_test(
        output_dir,
        input_dfs, screen_names,
        cases=qa_cases,
        controls=qa_controls,
        comp_name=f"{'_'.join(qa_cases)}_vs_{'_'.join(qa_controls)}",
        mut_col=mut_col,
        val_col=val_col,
        gene_col=gene_col,
        save_type='svg', 
    )
    print(f'Output from beclust3d.qc.hypothesis_tests is saved in the following directory: {output_dir}hypothesis_qc/')

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
    
    print("BE-QA complete: sequence/structural features and hypothesis tests")


if __name__ == '__main__':
    config_yaml = sys.argv[1]
    config = load_config(config_yaml)
    beclust3d_path = config['beclust3d_path']
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), beclust3d_path)))

    from beclust3d.lfc3d.structure import sequence_structural_features
    from beclust3d.qc.hypothesis_tests import hypothesis_test
    from beclust3d.lfc3d.conservation import conservation
    from beclust3d.lfc3d.preprocess_data_plot import plot_rawdata

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

    qa_controls = get_required(config, 'qa.controls', list)
    qa_cases    = get_required(config, 'qa.cases',    list)

    # OPTIONAL
    user_fasta              = get_optional(config, 'user_fasta',              (str, type(None)), None)
    user_pdb                = get_optional(config, 'user_pdb',                (str, type(None)), None)
    user_dssp               = get_optional(config, 'user_dssp',               (str, type(None)), None)
    structure_radius        = get_optional(config, 'structure_radius',        (int, float, type(None)), 6.0)
    priority_on_alternative = get_optional(config, 'priority_on_alternative', bool, False)
    atom_level_naa          = get_optional(config, 'atom_level_naa',          bool, False)
    muscle_path             = get_optional(config, 'muscle_path',             (str, type(None)), 'muscle')

    main(
        input_gene=input_gene, input_uniprot=input_uniprot, input_chain=input_chain,
        screen_dir=screen_dir, screens=screens, mut_col=mut_col, val_col=val_col,
        gene_col=gene_col, edits_col=edits_col, output_dir=output_dir,
        user_fasta=user_fasta, user_pdb=user_pdb, user_dssp=user_dssp,
        structure_radius=structure_radius, mut_categories=mut_categories,
        mut_delimiter=mut_delimiter, conservation_run=conservation_run,
        alt_gene_name=alt_gene_name, alt_uniprot_id=alt_uniprot_id,
        alt_screen_start=alt_screen_start, v_score_threshold=v_score_threshold,
        qa_controls=qa_controls, qa_cases=qa_cases,
        priority_on_alternative=priority_on_alternative,
        atom_level_naa=atom_level_naa, muscle_path=muscle_path,
        config_yaml=config_yaml,
    )

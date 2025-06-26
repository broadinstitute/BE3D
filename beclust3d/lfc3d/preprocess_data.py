"""
File: preprocess_data.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
import warnings
from pathlib import Path

from .preprocess_data_helpers import *

def parse_be_data(
    workdir, 
    input_dfs, input_gene, screen_names, 
    mut_col='Mutation category', 
    val_col='logFC', 
    gene_col='Target Gene Symbol', 
    edits_col='Amino Acid Edits', 
    mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"], 
    mut_delimiter = ',', 
    conserv_dfs=[], conserv_col='mouse_res_pos',
    conserv_score_col='v_score', ### conserv_col
    gene_list = False
): 
    """
    Parse raw base editing screen data and create separate DataFrames for each mutation type per screen.

    Parameters
    ----------
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_dfs : list of pd.DataFrame
        List of input dataframes, one for each screen.

    input_gene : str
        Name of the gene being processed. 

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in input_dfs.

    mut_col : str, optional (default='Mutation category')
        Column name in input_dfs specifying the mutation category (e.g., 'Missense', 'Nonsense').

    val_col : str, optional (default='logFC')
        Column name in input_dfs specifying the value measurement (e.g., log fold-change).

    gene_col : str, optional (default='Target Gene Symbol')
        Column name specifying the target gene name in input_dfs.

    edits_col : str, optional (default='Amino Acid Edits')
        Column name specifying the amino acid edits or mutation information in input_dfs.

    mut_categories : list of str, optional
        List of mutation categories to extract separately. 
        Default includes ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"].

    mut_delimiter : str, optional (default=',')
        Delimiter used to separate multiple mutations within the edits_col field.

    conserv_dfs : list of pd.DataFrame, optional (default=[])
        List of conservation DataFrames, one per screen, used to optionally filter mutations based on conserved residues.

    conserv_col : str, optional (default='mouse_res_pos')
        Column name in the conservation DataFrames containing residue positions considered conserved.

    Returns
    -------
    mut_dfs : dict
        Nested dictionary where:
          - Keys are screen names (from screen_names)
          - Values are dictionaries mapping mutation types (e.g., 'Missense') to processed DataFrames
            containing parsed mutation information and LFC values.
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'screendata'):
        os.mkdir(working_filedir / 'screendata')
    if not os.path.exists(working_filedir / 'screendata/plots'):
        os.mkdir(working_filedir / 'screendata/plots')

    # CHECK INPUTS ARE SELF CONSISTENT #
    for df in input_dfs: 
        assert mut_col in df.columns, 'Check [mut_col] input'
        assert val_col in df.columns, 'Check [val_col] input'
        assert gene_col in df.columns, 'Check [gene_col] input'
        assert edits_col in df.columns, 'Check [edits_col] input'
    # for df in conserv_dfs: 
    #     if df is not None: 
    #         assert conserv_col in df.columns, 'Check [conserv_col] input'
    #         assert conserv_score_col in df.columns, 'Check [conserv_col] input'

    assert len(input_dfs) == len(screen_names) == len(conserv_dfs), 'Lengths of [input_dfs] and [screen_names] and [conservation_dfs] must match'

    mut_dfs = {}
    # OUTPUT TSV BY INDIVIDUAL SCREENS #
    for input_gene, screen_df, screen_name, conserv_df in zip(gene_list, input_dfs, screen_names, conserv_dfs): 
        print('Processing', screen_name)
        # IF WE LOOK AT CONSERVATION #
        if conserv_df is not None: 
            conserv_list = [str(x) for x in conserv_df[conserv_df[conserv_score_col]==3][conserv_col].tolist()]
        # NARROW DOWN TO INPUT_GENE #
        df_gene = screen_df.loc[screen_df[gene_col] == input_gene, ]
        mut_dfs[screen_name] = {}

        # NARROW DOWN TO EACH MUTATION TYPE #
        gene_mut_df = {}
        for mut in mut_categories: 

            # MAKE SURE MUT CATEGORY APPEARS IN DF #
            if not mut in df_gene[mut_col].unique(): 
                warnings.warn(f'{mut} not in Dataframe')
                continue

            # IF USER WANTS TO CATEGORIZE BY ONE SINGLE MUTATION PER GUIDE OR MULTIPLE MUTATIONS PER GUIDE #
            df_mut = df_gene.loc[df_gene[mut_col] == mut, ]
            df_mut = df_mut.reset_index(drop=True)
            gene_mut_df[mut] = len(df_mut)

            # ASSIGN position refAA altAA #
            df_mut[edits_col] = df_mut[edits_col].str.strip(mut_delimiter) # CLEAN
            df_mut[edits_col] = df_mut[edits_col].str.split(mut_delimiter) # STR to LIST
            df_mut[edits_col] = df_mut[edits_col].apply(lambda xs: identify_mutations(xs)) # FILTER FOR MUTATIONS ONLY #

            df_exploded = df_mut.explode(edits_col) # EACH ROW IS A MUTATION #
            df_exploded['edit_pos'] = df_exploded[edits_col].str.extract('(\d+)')
            df_exploded['refAA'] = df_exploded[edits_col].str.extract('([A-Za-z*]+)')
            df_exploded['altAA'] = df_exploded[edits_col].str.extract('[A-Za-z]+\d+([A-Za-z*]+)$')
            # IF 3 LETTER CODES ARE USED, TRANSLATE TO 1 LETTER CODE #
            df_exploded['refAA'] = df_exploded['refAA'].str.upper().apply(lambda x: aa_map.get(x, x))
            df_exploded['altAA'] = df_exploded['altAA'].str.upper().apply(lambda x: aa_map.get(x, x))

            # FILTER OUT SCORES WHERE POS DOES NOT APPEAR IN CONSERVED #
            if conserv_df is not None: 
                df_exploded = df_exploded[df_exploded['edit_pos'].isin(conserv_list) | df_exploded['edit_pos'].isna() | (df_exploded['edit_pos'] == "")]

            df_subset = df_exploded[[edits_col, 'edit_pos', 'refAA', 'altAA', val_col]]
            df_subset = df_subset.rename(columns={edits_col: 'this_edit', val_col: 'LFC'})

            # FOR PARTICULAR MUTATIONS, NEED TO SUBSET FURTHER #
            if mut == 'Missense': 
                df_subset = df_subset[(df_subset['refAA'] != df_subset['altAA']) & (df_subset['altAA'] != '*')]
            elif mut == 'Silent': # SILENT BEFORE NONSENSE (ie *248* MUTATION IS SILENT NOT NONSENSE)
                df_subset = df_subset[df_subset['refAA'] == df_subset['altAA']]
            elif mut == 'Nonsense': 
                df_subset = df_subset[df_subset['altAA'] == '*']
            else: 
                df_subset = df_subset[df_subset['LFC'] != df_subset['LFC'].shift()]
                df_subset = df_subset['LFC']

            # WRITE LIST OF MUT AND THEIR LFC VALUES #
            screen_name_nospace, mut_nospace = screen_name.replace(' ','_'), mut.replace(' ','_')
            edits_filename = f"screendata/{input_gene}_{screen_name_nospace}_{mut_nospace}.tsv"
            df_subset.to_csv(working_filedir / edits_filename, sep='\t')
            
            mut_dfs[screen_name][mut] = df_subset
        
        print(gene_mut_df) # OUTPUT COUNTS PER GENE FOR EA MUTATION #

    return mut_dfs

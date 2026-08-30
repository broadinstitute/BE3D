"""
File: preprocess_data.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
    Parses raw base editing screen data into per-mutation-type DataFrames for each screen.
"""

import os
import warnings
import pandas as pd
from pathlib import Path

from .preprocess_data_helpers import *

def parse_be_data(
    workdir, 
    input_dfs, 
    input_gene, 
    screen_names, 
    mut_col='Mutation category', 
    val_col='logFC', 
    gene_col='Target Gene Symbol', 
    edits_col='Amino Acid Edits', 
    mut_categories=["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"], 
    mut_delimiter=',', 
    conserv_dfs=[],
    conserv_col='mouse_res_pos',
    v_score_threshold=3, ### conserv_col
    gene_list=False,
    mutation_priority=None,
    invert_score=False,
):
    """
    Parses raw base editing screen data into per-mutation-type DataFrames for each screen.

    Optionally filters mutations based on conservation scores from conservation DataFrames produced by conservation().

    Parameters
    ----------
    
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_dfs : list of pd.DataFrame
        List of input dataframes, one for each screen, each containing mutation category, gene, and value columns.
    
    input_gene : str
        Name of the gene being processed (e.g., 'DNMT3A', 'MEN1'). 
        
    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in input_dfs, used in plot labels and output filenames.

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
        Column name in conserv_dfs containing residue positions to filter on.
        
    v_score_threshold : int, optional (default=3)
        Conservation score for filtering. Scores are: -1 (not conserved), 1 (weakly similar), 2 (similar), 3 (conserved).
        
    gene_list : bool, optional (default=False)
        If True, processes a list of genes rather than a single gene.

    mutation_priority : list of str or None, optional (default=None)
        Priority order (most to least deleterious) used to collapse a guide's
        mut_col value into a single category when it is a delimiter-joined
        list of per-edit categories (e.g. 'Silent;Missense;'). If None, mut_col
        values are used as-is, assuming they are already single categories.

    invert_score : bool, optional (default=False)
        If True, negate the numeric value column (val_col) of every screen
        BEFORE any per-mutation parsing/aggregation, so that the sign of the
        resulting neg/pos channels reflects the intended biology.

        BE3D assumes a DROPOUT sign convention: loss-of-function (LOF) is
        expected to produce a NEGATIVE score, so LOF hits land in the negative
        channel and downstream interpretation treats negative as deleterious.
        Two classes of screen have the OPPOSITE polarity and should set
        invert_score=True:

          - Activity-reporter screens (e.g. a methylation->fluorescence
            reporter) where LOF gives a HIGH / POSITIVE signal.
          - Drug-resistance / positive-selection ENRICHMENT screens where the
            hit is enrichment, i.e. a POSITIVE score.

        Fed raw (invert_score=False), these screens silently flip the biology:
        LOF hits fall into the positive channel. Because the built-in QA uses
        two-sided tests (Mann-Whitney U / KS), it is unaffected by the sign and
        will still "pass", hiding the error. Setting invert_score=True fixes the
        polarity so LOF hits carry a negative score as the rest of the pipeline
        expects. Default (False) leaves all values unchanged.

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
            conserv_df['v_score'] = conserv_df['v_score'].astype(int)
            conserv_list = [str(x) for x in conserv_df[conserv_df['v_score']>=v_score_threshold][conserv_col].tolist()]
        # NARROW DOWN TO INPUT_GENE #
        df_gene = screen_df.loc[screen_df[gene_col] == input_gene, ]
        # INVERT SCORE SIGN FOR ENRICHMENT / ACTIVITY SCREENS (OPT-IN) #
        # Done BEFORE any per-mutation parsing so neg/pos channels carry the
        # intended biology. QA is two-sided, so it is unaffected by the sign.
        if invert_score:
            df_gene = df_gene.copy()
            df_gene[val_col] = -df_gene[val_col]
        if mutation_priority:
            df_gene = df_gene.copy()
            df_gene[mut_col] = df_gene[mut_col].apply(
                lambda x: reduce_mutation_type(x, mut_delimiter, mutation_priority))
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

def sanitary_check(df_struc, df_missense_list, mute=True):
    """
        Check how the number of missense edits mapped to the target protein.

        Parameters
        ----------
        df_struc : pd.DataFrame
            Dataframe for target structure-sequence information.

        df_missense_list : list of pd.DataFrame
            List of missense dataframes, one for each screen.

        Returns
        -------
        """    
    struc_refAA_pos_list = (df_struc['unires']+df_struc['unipos'].astype(str)).to_list()

    for each_df_missense in df_missense_list:
        missense_refAA_pos_list = each_df_missense['this_edit'].str[:-1].to_list()
        if not mute: 
            print('-----[SANITARY CHECK]-----')
            print(f'#of missense edits:{len(set(missense_refAA_pos_list))},\
                  #of mapped missense edits:{len(set(missense_refAA_pos_list).intersection(set(struc_refAA_pos_list)))},\
                  #of not mapped missense edits:{len(set(missense_refAA_pos_list).difference(set(struc_refAA_pos_list)))},\
                  list of not mapped missense edits: {list(set(missense_refAA_pos_list).difference(set(struc_refAA_pos_list)))}')

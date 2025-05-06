"""
File: hypothesis_tests.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-23
Description: 
"""

import os
from pathlib import Path

from .hypothesis_tests_helpers import *

def hypothesis_test(
    workdir, 
    input_dfs, screen_names, 
    cases, controls, 
    comp_name='CaseVsControl', 
    mut_col='Mutation category', 
    val_col='logFC', 
    gene_col='Target Gene Symbol', 
    save_type='png', 
): 
    """
    Conduct hypothesis 1 (one screen vs control from same screens) and hypothesis 2
    (one screen vs control from all screens) on the set of input screens and genes. 

    Parameters
    ----------
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_dfs : list of pd.DataFrame
        List of input dataframes, one for each screen.

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in input_dfs.

    cases : list of str
        List of mutation categories considered as the "case" group (e.g., ['Nonsense']).

    controls : list of str
        List of mutation categories considered as the "control" group (e.g., ['No Mutation', 'Silent']).

    comp_name : str
        Comparison label used for naming plots and outputs (e.g., 'Nonsense_vs_Control').
        
    mut_col : str, optional (default='Mutation category')
        Column name in input_dfs specifying the mutation category (e.g., 'Missense', 'Nonsense').

    val_col : str, optional (default='logFC')
        Column name in input_dfs specifying the value measurement (e.g., log fold-change).

    gene_col : str, optional (default='Target Gene Symbol')
        Column name specifying the target gene name in input_dfs.

    save_type : str, optional (default='png')
        Format for saving output plots (e.g., 'png', 'pdf').

    Returns
    -------
    df_MW1_input : pd.DataFrame
        DataFrame containing Mann-Whitney U test results for Hypothesis 1 (within each screen).

    df_MW2_input : pd.DataFrame
        DataFrame containing Mann-Whitney U test results for Hypothesis 2 (across all screens).

    df_KS1_input : pd.DataFrame
        DataFrame containing Kolmogorov-Smirnov test results for Hypothesis 1 (within each screen).

    df_KS2_input : pd.DataFrame
        DataFrame containing Kolmogorov-Smirnov test results for Hypothesis 2 (across all screens).
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'hypothesis_qc'):
        os.mkdir(working_filedir / 'hypothesis_qc')

    # CHECK INPUTS ARE SELF CONSISTENT #
    for df in input_dfs: 
        assert mut_col in df.columns, 'Check [mut_col] input'
        assert val_col in df.columns, 'Check [val_col] input'
        assert gene_col in df.columns, 'Check [gene_col] input'

    unique_genes = []
    for df in input_dfs: 
        unique = df[gene_col].unique().tolist()
        unique_genes = list(set(unique_genes+unique))

    unique_mutations = []
    for df in input_dfs: 
        unique = df[mut_col].unique().tolist()
        unique_mutations = list(set(unique_mutations+unique))

    for c in cases+controls: 
        assert c in unique_mutations, f'{c} not found in mutation types'

    assert len(input_dfs) == len(screen_names), 'Lengths of [input_dfs] and [screen_names] must match'

    # AGGREGATE ACROSS SCREENS FOR HYPOTHESIS #

    # MW AND KS TESTS HYPOTHESIS 1 #
    df_MW1_input = hypothesis_one(working_filedir, input_dfs, screen_names, unique_genes, cases, controls, comp_name, 
                                  gene_col, mut_col, val_col, testtype='MannWhitney')
    df_KS1_input = hypothesis_one(working_filedir, input_dfs, screen_names, unique_genes, cases, controls, comp_name, 
                                  gene_col, mut_col, val_col, testtype='KolmogorovSmirnov')
    
    if len(unique_genes) > 1:
        hypothesis_plot(working_filedir, df_MW1_input, df_KS1_input, screen_names, 'screenid', 'gene_name', 
                        testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='1', header=comp_name, save_type=save_type)
    if len(screen_names) > 1:
        hypothesis_plot(working_filedir, df_MW1_input, df_KS1_input, unique_genes, 'gene_name', 'screenid', 
                        testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='1', header=comp_name, save_type=save_type)

    # MW AND KS TESTS HYPOTHESIS 2 #
    df_MW2_input = hypothesis_two(working_filedir, input_dfs, screen_names, unique_genes, cases, controls, comp_name, 
                                  gene_col, mut_col, val_col, testtype='MannWhitney')
    df_KS2_input = hypothesis_two(working_filedir, input_dfs, screen_names, unique_genes, cases, controls, comp_name, 
                                  gene_col, mut_col, val_col, testtype='KolmogorovSmirnov')
    
    if len(unique_genes) > 1:
        hypothesis_plot(working_filedir, df_MW2_input, df_KS2_input, screen_names, 'screenid', 'gene_name', 
                        testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='2', header=comp_name, save_type=save_type)
    if len(screen_names) > 1:
        hypothesis_plot(working_filedir, df_MW2_input, df_KS2_input, unique_genes, 'gene_name', 'screenid', 
                        testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='2', header=comp_name, save_type=save_type)
    
    return df_MW1_input, df_MW2_input, df_KS1_input, df_KS2_input

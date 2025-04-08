"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1

"""

import os
import numpy as np
from pathlib import Path

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

from .preprocess_data_helpers import *

def plot_rawdata(
    input_dfs, 
    workdir, screen_names, 
    mut_col='Mutation category', gene_col='Target Gene Symbol', val_col='logFC', 
    mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"], 
): 
    """
    Description
        Parse raw data and create plots for each input screen.
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

    assert len(screen_names) == len(input_dfs), 'Lengths of [input_dfs] and [screen_names] must match'

    # INDIVIDUAL BARPLOTS AND VIOLIN PLOTS FOR EACH SCREEN #
    for df, screen_name in zip(input_dfs, screen_names): 
        counts_by_gene(df=df, working_filedir=working_filedir, 
                        gene_col=gene_col, mut_col=mut_col, title=screen_name, 
                        mut_categories=mut_categories)
        violin_by_gene(df=df, working_filedir=working_filedir, 
                        gene_col=gene_col, mut_col=mut_col, val_col=val_col, title=screen_name, 
                        mut_categories=mut_categories)

    return None

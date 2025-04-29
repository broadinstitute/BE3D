"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1

"""

import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def plot_rawdata(
    workdir, 
    input_dfs, screen_names, 
    mut_col='Mutation category', val_col='logFC', gene_col='Target Gene Symbol', 
    mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"], 
    save_type='png', 
): 
    """
    Parse raw data and create plots for each input screen.

    Parameters
    ----------
    workdir : str
        Path to the working directory where output files and results will be saved.

    input_dfs : list of pd.DataFrame
        List of input dataframes, one for each screen.

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in input_dfs.

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
    None
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
                        mut_categories=mut_categories, save_type=save_type)
        violin_by_gene(df=df, working_filedir=working_filedir, 
                        gene_col=gene_col, mut_col=mut_col, val_col=val_col, title=screen_name, 
                        mut_categories=mut_categories, save_type=save_type)

    return None

def counts_by_gene(
    df, 
    working_filedir, 
    gene_col, mut_col, 
    title, mut_categories, 
    save_type, 
): 
    """
    Description
        Graph a bar plot of counts by category per gene
    """

    # FIND UNIQUE GENES #
    unique_genes = sorted(df[gene_col].unique().tolist())
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))

    # COMPUTE MUTATION COUNTS FOR EACH GENE AND MUT CATEGORY #
    df_mutation_counts = (
        df.groupby([gene_col, mut_col])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=mut_categories, fill_value=0)
        .reset_index()
    )
    df_mutation_counts.columns = ['Gene'] + mut_categories

    # BARPLOT #
    df_plot = df_mutation_counts.melt("Gene", var_name="Mut Type", value_name="Count")
    sns.set_style("darkgrid")
    ax = sns.catplot(data=df_plot, kind="bar", col_wrap=plot_dim, 
                     x="Mut Type", y="Count", hue="Mut Type", col="Gene", sharex=False)

    for ax in ax.axes.flat:
        for container in ax.containers:
            ax.bar_label(container)
    plt.subplots_adjust(top=0.9)
    plt.suptitle(title)

    # SAVE BARPLOT #
    plotname = f"screendata/plots/{title}_barplot_by_muttype.pdf"
    plt.savefig(working_filedir / plotname, dpi=100, transparent=True, format=save_type)
    plt.close()
    del df_mutation_counts, df_plot

def violin_by_gene(
    df, 
    working_filedir, 
    gene_col, mut_col, val_col, 
    title, 
    mut_categories, 
    save_type, 
): 
    """
    Description
        Graph a violin plot of LFC values by category per gene
    """
    
    # FIND HOW MANY PLOTS NEEDED #
    unique_genes = sorted(df[gene_col].unique())
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))

    # VIOLIN PLOT SETUP #
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(nrows=plot_dim, ncols=plot_dim, sharex=False, sharey=True, 
                             figsize=(19,17), gridspec_kw={'hspace':0.3, 'wspace':0.1})
    if plot_dim == 1: axes = [axes]
    else: axes = axes.flatten()

    legend_handles, legend_labels = None, None

    for idx, current_gene in enumerate(unique_genes): 
        df_gene = pd.DataFrame()
        # AGGREGATE OVER EVERY SCREEN #
        df_current_gene = df.loc[df[gene_col] == current_gene,]

        # ENSURE ALL MUTATION TYPES PRESENT #
        for mutcat in mut_categories: 
            if mutcat not in df_current_gene[mut_col].unique():
                df_temp = pd.DataFrame({val_col: [np.nan], mut_col: [mutcat]})
                df_current_gene = pd.concat([df_current_gene, df_temp], ignore_index=True)
        df_gene = pd.concat([df_gene, df_current_gene], ignore_index=True)
        del df_current_gene

        # CALC MEAN STD #
        filtered_df_gene = df_gene[df_gene[mut_col].isin(mut_categories)]
        Means = filtered_df_gene.groupby(mut_col)[val_col].mean()

        # PLOT VIOLIN #
        df_gene.loc[:, mut_col] = pd.Categorical(df_gene[mut_col], categories=mut_categories)
        df_gene = df_gene.sort_values(by=[mut_col]).reset_index(drop=True)
        violin = sns.violinplot(ax=axes[idx], data=df_gene, x=val_col, y=mut_col, 
                                inner=None, hue=mut_col)
        axes[idx].set_title(current_gene)

        # EXTRACT LEGEND #
        if legend_handles is None and legend_labels is None:
            legend_handles, legend_labels = axes[idx].get_legend_handles_labels()

        # REMOVE LEGEND FROM PLOTS #
        if axes[idx].get_legend() is not None: 
            axes[idx].get_legend().remove()
            axes[idx].axvline(df_gene[val_col].mean(), c="gray", linestyle="dashed")
            axes[idx].scatter(y=range(len(Means)), x=Means, c="violet", alpha=.9) # MEANS #
        del df_gene
    
    # REMOVE UNUSSED AXES
    for j in range(idx + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.suptitle(title)
    plt.subplots_adjust(top=0.9, wspace=0.1)

    # SAVE VIOLIN #
    plotname = f"screendata/plots/{title}_violinplot_by_muttype.pdf"
    plt.savefig(working_filedir / plotname, dpi=100, transparent=True, format=save_type)

    # CREATE SEPARATE LEGEND #
    if legend_handles and legend_labels:
        legend_fig = plt.figure(figsize=(4, 2))  # Adjust size as needed
        legend_ax = legend_fig.add_subplot(111)
        legend_ax.axis('off')  # Hide axes
        legend_ax.legend(legend_handles, legend_labels, loc='center')
        
        legend_path = f"screendata/plots/{title}_legend.pdf"
        legend_fig.savefig(working_filedir / legend_path, dpi=100, transparent=True, format=save_type)

    plt.close()
    
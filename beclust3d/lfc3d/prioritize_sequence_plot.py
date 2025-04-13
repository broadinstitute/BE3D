"""
File: prioritize_sequence_plot.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def plot_screendata_sequence(
    df_protein, 
    workdir, 
    input_gene, screen_name, function_name='mean', muttype='Missense', 
    save_type='png', 
): 
    """
    Description
        Parse raw data and create plots for each input screen.
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'screendata_sequence'):
        os.mkdir(working_filedir / 'screendata_sequence')
    if not os.path.exists(working_filedir / 'screendata_sequence/plots'):
        os.mkdir(working_filedir / 'screendata_sequence/plots')

    # PLOT SCATTERPLOT AND COUNTS PLOT #
    counts_by_residue(df_protein, working_filedir, 
                      input_gene=input_gene, 
                      screen_name=screen_name, 
                      muttype=muttype, 
                      save_type=save_type)

    stdev_by_residue(df_protein, working_filedir, 
                     input_gene=input_gene, 
                     screen_name=screen_name, 
                     muttype=muttype, 
                     colname=f'{function_name}_{muttype}_LFC', 
                     colname_stdev=f'{function_name}_{muttype}_LFC_stdev', 
                     centered=False, 
                     save_type=save_type)
    stdev_by_residue(df_protein, working_filedir, 
                     input_gene=input_gene, 
                     screen_name=screen_name, 
                     muttype=muttype, 
                     colname=f'{function_name}_{muttype}_LFC', 
                     colname_stdev=f'{function_name}_{muttype}_LFC_stdev', 
                     centered=True, 
                     save_type=save_type)
    
    scatterplot_by_residue(df_protein, working_filedir, 
                           input_gene=input_gene, 
                           screen_name=screen_name, 
                           muttype=function_name, 
                           colname=f'{function_name}_{muttype}_LFC', 
                           save_type=save_type)
    
    dual_scatterplot_by_residue(df_protein, working_filedir, 
                                input_gene=input_gene, 
                                screen_name=screen_name, 
                                muttype=muttype, 
                                colname=f'{function_name}_{muttype}_LFC', 
                                colname_z=f'{function_name}_{muttype}_LFC_Z', 
                                colname_plab=f'{function_name}_{muttype}_LFC_plab', 
                                save_type=save_type)
    dual_histogram_by_residue(df_protein, working_filedir, 
                              input_gene=input_gene, 
                              screen_name=screen_name, 
                              muttype=muttype, 
                              colname=f'{function_name}_{muttype}_LFC', 
                              colname_plab=f'{function_name}_{muttype}_LFC_plab', 
                              save_type=save_type)

    return None


def counts_by_residue(
    df_struc_consvr, 
    working_filedir, 
    input_gene, screen_name, muttype, 
    save_type, 
): 
    # PREP DATA #
    counts = df_struc_consvr[f'all_{muttype}_edits'].str.count(';').fillna(0).astype(int)+1
    counts[df_struc_consvr[f'all_{muttype}_edits'] == '-'] = 0

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax = sns.barplot(x=df_struc_consvr['unipos'], y=counts, 
                     color='steelblue', edgecolor='steelblue')
    ax.set_ylabel(f"Count of {muttype} Mutations")
    ax.set_xlabel(f"unipos")
    ax.set_title(f"{input_gene} Count of {muttype} Mutations {screen_name}")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    counts_filename = f"screendata_sequence/plots/{input_gene}_{screen_name}_num_{muttype}_per_res.{save_type}"
    plt.savefig(working_filedir / counts_filename, dpi=100, transparent=True, format=save_type)
    plt.close(fig)

def stdev_by_residue(
    df_struc_consvr, 
    working_filedir, 
    input_gene, screen_name, muttype, 
    colname, colname_stdev, 
    centered,
    save_type, 
): 
    # PREP DATA #
    xvals = df_struc_consvr['unipos']
    yvals = pd.to_numeric(df_struc_consvr[colname], errors='coerce').fillna(0)
    stdevs = pd.to_numeric(df_struc_consvr[colname_stdev], errors='coerce').fillna(0)
    xvals_filtered = xvals[stdevs != 0]
    if not centered: yvals_filtered = yvals[stdevs != 0]
    else:            yvals_filtered = [0]*len(xvals_filtered)
    stdevs_filtered = stdevs[stdevs != 0]

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.errorbar(x=xvals_filtered, y=yvals_filtered, yerr=stdevs_filtered, 
                color='steelblue', ls=' ', marker='o', capsize=3, capthick=1, ecolor='black', 
                )
    ax.set_ylabel(f"Standard Deviations of {muttype} Mutations")
    ax.set_xlabel(f"unipos")
    ax.set_title(f"{input_gene} Standard Dev of {muttype} Mutations {screen_name}")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    stdev_filename = f"screendata_sequence/plots/{input_gene}_{screen_name}_stdev_{muttype}_per_res.{save_type}"
    plt.savefig(working_filedir / stdev_filename, dpi=100, transparent=True, format=save_type)
    plt.close(fig)
    
def scatterplot_by_residue(
    df_struc_consvr, 
    working_filedir, 
    input_gene, screen_name, muttype, colname, 
    save_type, 
): 
    # PREP DATA #
    x_list = df_struc_consvr['unipos'].tolist()
    y_list = df_struc_consvr[colname].tolist()
    x_vals = [x for x, y in zip(x_list, y_list) if y!='-']
    y_vals = [float(y) for y in y_list if y!='-']

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.axhline(-1.0, c="red", linestyle="--")
    ax.axhline(1.0, c="blue", linestyle="--")
    ax.axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=ax, x=x_vals, y=y_vals, color='steelblue', edgecolor='steelblue')
    ax.set_ylabel(f"{muttype} LFC Score")
    ax.set_xlabel(f"unipos")
    ax.set_title(f'{input_gene} {muttype} LFC Score By Residue {screen_name}')
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    scatter_filename = f"screendata_sequence/plots/{input_gene}_{screen_name}_{colname}_score_by_res.{save_type}"
    plt.savefig(working_filedir / scatter_filename, dpi=100, transparent=True, format=save_type)
    plt.close(fig)

def dual_scatterplot_by_residue(
    df_struc_consvr, 
    working_filedir, 
    input_gene, screen_name, muttype, 
    colname, colname_z, colname_plab, 
    save_type, 
): 
    df_struc_consvr = df_struc_consvr[df_struc_consvr[colname] != '-']
    df_struc_consvr[colname] = df_struc_consvr[colname].astype(float)
    df_struc_consvr_pos = df_struc_consvr[df_struc_consvr[colname] > 0.0]
    df_struc_consvr_neg = df_struc_consvr[df_struc_consvr[colname] < 0.0]

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 6))
    for ax in axs: 
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

    axs[0].axhline(-1.0, c="red", linestyle="--")
    axs[0].axhline(1.0, c="blue", linestyle="--")
    axs[0].axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=axs[0], data=df_struc_consvr_pos, x="unipos", 
                    y=colname_z, hue=colname_plab, palette='tab10')
    axs[0].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    axs[0].set_title(f'Positive LFC Values')

    axs[1].axhline(-1.0, c="red", linestyle="--")
    axs[1].axhline(1.0, c="blue", linestyle="--")
    axs[1].axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=axs[1], data=df_struc_consvr_neg, x="unipos", 
                    y=colname_z, hue=colname_plab, palette='tab10')
    axs[1].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    axs[1].set_title(f'Negative LFC Values')

    plt.subplots_adjust(wspace=0.3)
    plt.suptitle(f'{input_gene} LFC_Z Score {screen_name}')

    scatter_filename = f"screendata_sequence/plots/{input_gene}_{screen_name}_{muttype}_lfcz_scatter_by_bin_posneg.{save_type}"
    plt.savefig(working_filedir / scatter_filename, dpi=100, transparent=True, format=save_type)
    plt.close(fig)

def dual_histogram_by_residue(
    df_struc_consvr, 
    working_filedir, 
    input_gene, screen_name, muttype, 
    colname, colname_plab, 
    save_type, 
):  
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 6))
    for ax in axs: 
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

    df_struc_consvr = df_struc_consvr[df_struc_consvr[colname] != '-']
    df_struc_consvr[colname] = df_struc_consvr[colname].astype(float)
    df_struc_consvr_pos = df_struc_consvr[df_struc_consvr[colname] > 0.0]
    df_struc_consvr_neg = df_struc_consvr[df_struc_consvr[colname] < 0.0]
    plot1 = sns.histplot(ax=axs[0], data=df_struc_consvr_pos, x=colname, hue=colname_plab, bins=80, palette='tab10')
    plot2 = sns.histplot(ax=axs[1], data=df_struc_consvr_neg, x=colname, hue=colname_plab, bins=80, palette='tab10')

    # give corresponding titles
    plot1.set_title(f'Positive LFC Counts')
    plot2.set_title(f'Negative LFC Counts')
    plt.suptitle(f'{input_gene} Mean Missense LFC Counts {screen_name}')
    plt.subplots_adjust(wspace=0.1)

    hist_filename = f"screendata_sequence/plots/{input_gene}_{screen_name}_{muttype}_lfc_hist_by_bin_posneg.{save_type}"
    plt.savefig(working_filedir / hist_filename, dpi=100, transparent=True, format=save_type)
    plt.close(fig)

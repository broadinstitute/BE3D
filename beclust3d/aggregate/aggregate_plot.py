"""
File: aggregate_plot.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-01-01
Description: 
"""

import os
import warnings
from pathlib import Path
import seaborn as sns
import matplotlib.pylab as plt

import statistics
import scipy.stats as stats
from scipy.stats import mannwhitneyu

from .aggregate_helpers import *

def average_split_bin_plots(
    df_z, workdir, input_gene, pthr=0.05, 
    screen_name='', func='SUM', score_type='LFC3D', 
    aggregate_dir='meta-aggregate', save_type='png', 
): 
    """
    Generates histograms, histplots, and scatterplots for positive and negative scores 
    with binning and significance thresholds.

    Parameters
    ----------
    df_z : pd.DataFrame
        DataFrame containing z-scores, p-values, and significance labels for scores at multiple thresholds.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    pthr : float, optional (default=0.05)
        p-value threshold for labeling statistical significance.

    screen_name : str
        Name of the screens corresponding to df_missense.

    func : str, optional (default='SUM')
        Name corresponding to 'aggr_func' (e.g., 'SUM', 'MEAN').

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed (e.g., 'LFC3D', 'LFC', etc.).

    aggregate_dir : str, optional (default='meta-aggregate')
        Subdirectory name where plots are stored (e.g., 'meta-aggregate', 'LFC3D', 'LFC').

    save_type : str, optional (default='png')
        Format for saving output plots (e.g., 'png', 'pdf').
    
    Returns
    ----------
    None
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / f'{aggregate_dir}/plots'):
        os.mkdir(working_filedir / f'{aggregate_dir}/plots')
    
    assert aggregate_dir == 'meta-aggregate' or aggregate_dir == 'LFC' or aggregate_dir == 'LFC3D'

    # SETUP PARAMS #
    pthr_str = str(pthr).split('.')[1]
    output_prefix = f"{input_gene}_{screen_name}_{score_type}_cutoff{pthr_str}"
    output_prefix = output_prefix.replace('__','_')

    neg_label = '_'.join([screen_name, func, score_type, 'neg'])
    neg_label = neg_label.replace('__', '_').strip('_')
    pos_label = '_'.join([screen_name, func, score_type, 'pos'])
    pos_label = pos_label.replace('__', '_').strip('_')

    # HISTOGRAMS #
    if screen_name == '': 
        histogram_params = [(f'{func}_{score_type}r_neg', neg_label, 'Negative'), 
                            (f'{func}_{score_type}r_pos', pos_label, 'Positive') ] # META #
    else: 
        histogram_params = [(f'{screen_name}_AVG_{score_type}r_neg', neg_label, 'Negative'), 
                            (f'{screen_name}_AVG_{score_type}r_pos', pos_label, 'Positive') ] # NON AGGR #
    
    histogram_filename = f"{output_prefix}_signal_vs_background.{save_type}"
    histogram_filename = histogram_filename.replace('__','_')
    (res_neg, res_pos) = metaaggregation_histogram(
        df_z, histogram_params, 
        working_filedir / f"{aggregate_dir}/plots/{histogram_filename}", save_type)

    # ERROR #
    if res_neg is None or res_pos is None: 
        warnings.warn('Invalid input data.')
        return None
    
    # BINNING FOR FUTURE STEPS #
    df_z = binning_lfc3d(df_z, neg_label, pos_label)

    # HISPLOTS #
    hisplots_params = [(f'{neg_label}_dis', neg_label, f'{neg_label}_{pthr_str}_psig', 'Negative P-Value'), 
                       (f'{neg_label}_dis', neg_label, f'{neg_label}_dis', 'Negative P-Value'), 
                       (f'{pos_label}_dis', pos_label, f'{pos_label}_{pthr_str}_psig', 'Positive P-Value'), 
                       (f'{pos_label}_dis', pos_label, f'{pos_label}_dis', 'Positive P-Value') ]
    
    histplot_filename = f"{output_prefix}_histplot.{save_type}"
    histplot_filename = histplot_filename.replace('__','_')
    metaaggregation_hisplot(
        df_z, hisplots_params, 
        working_filedir / f"{aggregate_dir}/plots/{histplot_filename}", save_type)

    # SCATTERPLOT #
    scatterplot_params = [(f'{neg_label}_dis', f'{neg_label}_{pthr_str}_psig', neg_label, 'Negative'), 
                          (f'{pos_label}_dis', f'{pos_label}_{pthr_str}_psig', pos_label, 'Positive')]
    
    scatterplot_filename = f"{output_prefix}_scatter_cutoff.{save_type}"
    scatterplot_filename = scatterplot_filename.replace('__','_')
    metaaggregation_scatterplot(
        df_z, scatterplot_params, input_gene, pthr, 
        working_filedir / f"{aggregate_dir}/plots/{scatterplot_filename}", save_type, colors=False)
    
    # Z SCORE SCATTERPLOT #
    scatterplot_params = [(f'{neg_label}_dis', f'{neg_label}_dis', f'{neg_label}_{pthr_str}_z', 'Negative'), 
                          (f'{pos_label}_dis', f'{pos_label}_dis', f'{pos_label}_{pthr_str}_z', 'Positive')]
    
    scatterplot_filename = f"{output_prefix}_scatter_colored.{save_type}"
    scatterplot_filename = scatterplot_filename.replace('__','_')
    metaaggregation_scatterplot(
        df_z, scatterplot_params, input_gene, pthr, 
        working_filedir / f"{aggregate_dir}/plots/{scatterplot_filename}", save_type, colors=True)


def metaaggregation_histogram(
    df_input, params, out_filename, save_type, 
): 
    """
    Description
        Helper function to plot histograms of the values along the length of the gene
    """
    
    fig, ax = plt.subplots(1, len(params), figsize=(12, 5), dpi=100)
    results_list = []

    for i, (avg, sum, out) in enumerate(params): 
        res = {}
        # PICK OUT DATA FOR PLOTTING #
        df_plot = pd.DataFrame()
        df_plot['unipos'] = df_input['unipos']
        df_plot[sum] = df_input[sum].replace('-', np.nan).astype(float)
        df_plot[avg] = df_input[avg].replace('-', np.nan).astype(float)
        df_filtered = df_plot.dropna(subset=[sum, avg])

        # MW AND PEARSON TESTS #
        U1, p = mannwhitneyu(df_plot[sum], df_plot[avg], method="asymptotic" )
        res['mannwhitneyu U1'], res['mannwhitneyu p'] = U1, p
        r, p = stats.pearsonr(df_filtered[sum].tolist(), df_filtered[avg].tolist())
        res['pearsonr r'], res['pearsonr p'] = r, p

        # SUM #
        res['sum min'], res['sum mean'] = df_plot[sum].min(), df_plot[sum].mean()
        res['sum med'], res['sum std'] = df_plot[sum].median(), df_plot[sum].std()

        if res['sum std'] == 0: 
            results_list.append(res)
            continue

        # Z AND P VALUE #
        z = statistics.NormalDist(mu=res['sum mean'], sigma=res['sum std']).zscore(-4.6)
        res['z'], res['p cdf'], res['p sf'] = z, stats.norm.cdf(z), stats.norm.sf(abs(z))

        # AVG #
        res['avg min'], res['avg mean'] = df_plot[avg].min(), df_plot[avg].mean()
        res['avg med'], res['avg std'] = df_plot[avg].median(), df_plot[avg].std()

        # PLOT #
        plot = df_plot.plot.area(x='unipos', alpha=0.5, stacked = False, ax=ax[i])
        plot.legend_.set_title(None)
        ax[i].axhline(y = res['sum mean'], color = 'r', linestyle = '-')
        ax[i].axhline(y = res['sum mean']-res['sum std'], color = 'r', linestyle = '--')
        ax[i].set_xticks(np.arange(0, len(df_input), 100))
        ax[i].set_title(out)

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

        del df_plot, df_filtered
        results_list.append(res)
    
    plt.subplots_adjust(wspace=0.15)
    plt.savefig(out_filename, dpi=100, transparent=False, format=save_type)
    plt.close()
    return results_list

def metaaggregation_hisplot(
    df_input, params, out_name, save_type, 
): 
    """
    Description
        Helper function to plot the distributions for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, len(params), figsize=(24, 5), dpi=100)

    for i, (dis, x, hue, name) in enumerate(params): 
    
        df_clean = df_input.loc[df_input[dis] != '-', ].reset_index(drop=True)
        df_clean[x] = df_clean[x].astype(float)
        plot = sns.histplot(df_clean, x=x, hue=hue, bins=50, palette='tab10', ax=ax[i])
        plot.legend_.set_title(None)
        ax[i].set_title(name)

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

    plt.subplots_adjust(wspace=0.15)
    plt.savefig(out_name, dpi=100, transparent=False, format=save_type)
    plt.close()
    return None

def metaaggregation_scatterplot(
    df_meta, params, input_gene, pthr, 
    outname, save_type, colors, 
): 
    """
    Description
        Helper function to plot scatterplots for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 2, figsize=(12, 6), dpi=300)

    for i, (dis, pval, y, out) in enumerate(params): 

        df_combined_clean = df_meta.loc[df_meta[dis] != '-', ]
        df_combined_clean[y] = df_combined_clean[y].astype(float)
        # MULTIPLE COLORS #
        if colors: 
            if 'pos' in dis: factor=1
            if 'neg' in dis: factor=-1
            ax[i].axhline(y = factor*1.65, color = 'r', linestyle = '--')
            ax[i].axhline(y = factor*1.96, color = 'r', linestyle = '--')
            ax[i].axhline(y = factor*2.58, color = 'r', linestyle = '--')

            plot = sns.scatterplot(data=df_combined_clean, x="unipos", y=y, 
                                   hue=pval, palette='tab10', ax=ax[i])

        # ABOVE AND BELOW THRESHOLD #
        else: 
            df_combined_psig = df_meta.loc[df_meta[pval] == 'p>='+str(pthr), ]
            line_list = df_combined_psig[y][df_combined_psig[y] != '-'].astype(float)
            if 'pos' in dis: line_val = line_list.max()
            if 'neg' in dis: line_val = line_list.min()
            ax[i].axhline(y = line_val, color = 'r', linestyle = '--')
            plot = sns.scatterplot(data=df_combined_clean, x="unipos", y=y, 
                                   hue=pval, palette='tab10', ax=ax[i])
            
        plot.legend_.set_title(None)
        ax[i].set_xticks(np.arange(0, len(df_meta), 100))
        ax[i].set_title(f"{input_gene} {out}")

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

    plt.subplots_adjust(wspace=0.15)
    plt.savefig(outname, dpi=100, transparent=False, format=save_type)
    plt.close()

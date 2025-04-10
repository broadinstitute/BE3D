"""
File: aggregate_plot.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-01-01
Description: 
"""

import os
from pathlib import Path
import seaborn as sns
import matplotlib.pylab as plt

import statistics
import scipy.stats as stats
from scipy.stats import mannwhitneyu

from .aggregate_helpers import *

def average_split_bin_plots(
        df_Z, workdir, input_gene, pthr=0.05, 
        name='', func='SUM', score_type='LFC3D', 
): 
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    pthr_str = str(pthr).split('.')[1]
    neg = '_'.join([name, func, score_type, 'neg']).replace('__', '_').strip('_')
    pos = '_'.join([name, func, score_type, 'pos']).replace('__', '_').strip('_')

    # HISTOGRAMS #
    if name == '': 
        histogram_params = [(f'{func}_{score_type}r_neg', neg, 'Negative'), 
                            (f'{func}_{score_type}r_pos', pos, 'Positive'), ]
    else: 
        histogram_params = [(f'{name}_AVG_{score_type}r_neg', neg, 'Negative'), 
                            (f'{name}_AVG_{score_type}r_pos', pos, 'Positive'), ]
    res_neg, res_pos = metaaggregation_histogram(df_Z, histogram_params, 
                                                 edits_filedir / f"plots/{input_gene}_{name}_signal_vs_background.png" )

    if res_neg is None or res_pos is None: return None
    df_Z = binning_lfc3d(df_Z, neg, pos)

    # HISPLOTS #
    hisplots_params = [(f'{neg}_dis', neg, f'{neg}_{pthr_str}_psig', 'Negative P-Value'), 
                       (f'{neg}_dis', neg, f'{neg}_dis', 'Negative P-Value'), 
                       (f'{pos}_dis', pos, f'{pos}_{pthr_str}_psig', 'Positive P-Value'), 
                       (f'{pos}_dis', pos, f'{pos}_dis', 'Positive P-Value'), ]
    metaaggregation_hisplot(df_Z, hisplots_params, 
                            edits_filedir / f"plots/{input_gene}_{name}_{score_type}_histogram.png" )

    # SCATTERPLOT #
    scatterplot_params = [(f'{neg}_dis', f'{neg}_{pthr_str}_psig', neg, 'Negative'), 
                          (f'{pos}_dis', f'{pos}_{pthr_str}_psig', pos, 'Positive')]
    metaaggregation_scatterplot(df_Z, scatterplot_params, input_gene, pthr, 
                                edits_filedir / f"plots/{input_gene}_{name}_{score_type}_scatter.png" )
    
    # Z SCORE SCATTERPLOT #
    scatterplot_params = [(f'{neg}_dis', f'{neg}_dis', f'{neg}_{pthr_str}_z', 'Negative'), 
                          (f'{pos}_dis', f'{pos}_dis', f'{pos}_{pthr_str}_z', 'Positive')]
    metaaggregation_scatterplot(df_Z, scatterplot_params, input_gene, pthr, 
                                edits_filedir / f"plots/{input_gene}_{name}_{score_type}_scatter_colored.png", colors=True )

def metaaggregation_histogram(
        df_meta, params, out_filename, 
): 
    """
    Description
        Helper function to plot histograms of the values along the length of the gene
    """
    
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), dpi=300)
    results_list = []

    for i, (avg, sum, out) in enumerate(params): 
        res = {}
        df_meta_plot = pd.DataFrame()
        df_meta_plot['unipos'] = df_meta['unipos']
        df_meta_plot[sum] = df_meta[sum].replace('-', np.nan).astype(float) ### can we fix the default format 250120
        df_meta_plot[avg] = df_meta[avg].replace('-', np.nan).astype(float) ### can we fix the default format 250120
        # df_meta_plot_sum = df_meta_plot[sum].dropna().tolist()
        # df_meta_plot_avg = df_meta_plot[avg].dropna().tolist()
        df_meta_filtered = df_meta_plot.dropna(subset=[sum, avg])
        df_meta_plot_sum = df_meta_filtered[sum].tolist()
        df_meta_plot_avg = df_meta_filtered[avg].tolist()

        U1, p = mannwhitneyu(df_meta_plot[sum], df_meta_plot[avg], method="asymptotic" )
        res['mannwhitneyu U1'], res['mannwhitneyu p'] = U1, p
        r, p = stats.pearsonr(df_meta_plot_sum, df_meta_plot_avg )
        res['pearsonr r'], res['pearsonr p'] = r, p
        # SUM #
        res['sum min'], res['sum mean'] = df_meta_plot[sum].min(), df_meta_plot[sum].mean()
        res['sum med'], res['sum std'] = df_meta_plot[sum].median(), df_meta_plot[sum].std()

        if res['sum std'] == 0: 
            return None
        z = statistics.NormalDist(mu=res['sum mean'], sigma=res['sum std']).zscore(-4.6)
        res['z'], res['p cdf'], res['p sf'] = z, stats.norm.cdf(z), stats.norm.sf(abs(z))
        # AVG #
        res['avg min'], res['avg mean'] = df_meta_plot[avg].min(), df_meta_plot[avg].mean()
        res['avg med'], res['avg std'] = df_meta_plot[avg].median(), df_meta_plot[avg].std()

        # PLOT #
        df_meta_plot.plot.area(x='unipos', alpha=0.55, stacked = False, ax=ax[i])
        ax[i].axhline(y = res['sum mean'], color = 'r', linestyle = '-')
        ax[i].axhline(y = res['sum mean']-res['sum std'], color = 'r', linestyle = '--')

        ax[i].legend(loc='lower left', borderaxespad=0)
        ax[i].set_xticks(np.arange(0,len(df_meta), 100))
        ax[i].set_title(out)

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

        del df_meta_plot
        results_list.append(res)
    
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_filename, dpi=300)
    return results_list[0], results_list[1]

def metaaggregation_hisplot(
        df_meta, params, out_name
): 
    """
    Description
        Helper function to plot the distributions for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 4, figsize=(36, 6), dpi=300)

    for i, (dis, x, hue, name) in enumerate(params): 
    
        df_combined_clean = df_meta.loc[df_meta[dis] != '-', ].reset_index(drop=True)
        df_combined_clean[x] = df_combined_clean[x].astype(float)
        sns.histplot(df_combined_clean, x=x, hue=hue, bins=50, palette='tab10', ax=ax[i])
        ax[i].set_title(name)

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_name, dpi=300) 

def metaaggregation_scatterplot(
        df_meta, params, input_gene, pthr, 
        outname, colors=False, 
): 
    """
    Description
        Helper function to plot scatterplots for the top 10 and bottom 10 % of points
    """
    fig, ax = plt.subplots(1, 2, figsize=(18, 6), dpi=300)

    for i, (dis, pval, y, out) in enumerate(params): 

        df_combined_clean = df_meta.loc[df_meta[dis] != '-', ]
        if colors: # MULTIPLE COLORS
            if 'pos' in dis: 
                ax[i].axhline(y = 1.65, color = 'r', linestyle = '--')
                ax[i].axhline(y = 1.96, color = 'r', linestyle = '--')
                ax[i].axhline(y = 2.58, color = 'r', linestyle = '--')
            if 'neg' in dis: 
                ax[i].axhline(y = -1.65, color = 'r', linestyle = '--')
                ax[i].axhline(y = -1.96, color = 'r', linestyle = '--')
                ax[i].axhline(y = -2.58, color = 'r', linestyle = '--')
            sns.scatterplot(data=df_combined_clean, x="unipos", y=y, hue=pval, palette='tab10', ax=ax[i])

        else: # 2 COLORS #
            df_combined_psig = df_meta.loc[df_meta[pval] == 'p>='+str(pthr), ]
            if 'pos' in dis: 
                line_val = df_combined_psig[y][df_combined_psig[y] != '-'].astype(float).max()
            if 'neg' in dis: 
                line_val = df_combined_psig[y][df_combined_psig[y] != '-'].astype(float).min()
            ax[i].axhline(y = line_val, color = 'r', linestyle = '--')
            sns.scatterplot(data=df_combined_clean, x="unipos", y=y, hue=pval, palette='tab10', ax=ax[i])

        ax[i].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
        ax[i].set_xticks(np.arange(0, len(df_meta), 100))
        ax[i].set_title(f"{input_gene} {out}")

        # SET BACKGROUND #
        ax[i].set_facecolor('#EBEBEB')
        [ax[i].spines[side].set_visible(False) for side in ax[i].spines]
        ax[i].grid(which='major', color='white', linewidth=0.5)
        ax[i].set_axisbelow(True)

    plt.subplots_adjust(wspace=0.3)
    plt.savefig(outname, dpi=300)

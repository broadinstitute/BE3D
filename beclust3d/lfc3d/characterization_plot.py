"""
File: characterization_plot.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-04
Description: 
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

# PLOT ENRICHMENT #

def plot_enrichment_test(
    enrichment_results, 
    workdir, 
    input_gene, hit_value, feature_values, 
    padding=0.5, save_type='png', 
):
    """
    Description
        Plot enrichment test results
    """
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'characterization/plots'):
        os.mkdir(working_filedir / 'characterization/plots')

    fig, ax = plt.subplots(figsize=(8, len(feature_values)))
    y_positions = [int(i)+1 for i in range(len(feature_values))]

    for i, result in enumerate(enrichment_results):
        odds_ratio = result['odds_ratio']
        ci = result['ci']
        y = y_positions[i]
        color = 'red' if i % 2 == 0 else 'blue'
        if np.isnan(odds_ratio) or np.isinf(odds_ratio):
            x_min, x_max = ax.get_xlim()
            x_mid = (x_min + x_max) / 2
            ax.plot([x_mid], [y], 'o', color="grey", markerfacecolor='none', linestyle='None')
            continue

        # DETERMINE STYLING #
        is_significant = result['p_value'] <= hit_value
        marker_style = 'o' if is_significant else 'o'
        marker_fill = color if is_significant else 'none'
        line_style = '-' if is_significant else ':'

        # Error bars for valid odds ratios
        error = [[odds_ratio - ci.low], [ci.high - odds_ratio]]
        ax.errorbar(
            x=[odds_ratio], y=[y],
            xerr=error, fmt=marker_style,
            color=color, linestyle=line_style,
            markerfacecolor=marker_fill
        )

    # CUSTOMIZE PLOT #
    ax.set_yticks(y_positions)
    ax.set_yticklabels(feature_values)
    ax.set_ylim(min(y_positions) - padding, max(y_positions) + padding)
    ax.set_xlabel('Odds Ratio')
    ax.set_title(f'{input_gene} Enrichment Test Odds Ratios')
    
    out_filename = working_filedir / f"characterization/plots/{input_gene}_enrichment_test.{save_type}"
    plt.savefig(out_filename, dpi=100, transparent=True, format=save_type)
    plt.close()

# CHARACTERIZATION PLOTS #

def lfc_lfc3d_scatter(
    df_input, 
    workdir, 
    input_gene, screen_name, 
    lfc3d_hit_threshold=0.05, save_type='png', 
): 
    """
    Description
        Generate LFC vs LFC3D scatter plot
    """
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'characterization'):
        os.mkdir(working_filedir / 'characterization')
    if not os.path.exists(working_filedir / 'characterization/plots'):
        os.mkdir(working_filedir / 'characterization/plots')

    # Load LFC and LFC3D scores and distributions
    df_input.rename(columns={
        f"{screen_name}_LFC": "LFC",
        f"{screen_name}_LFC3D": "LFC3D",
        f"{screen_name}_LFC3D_dis": "LFC3D_dis", 
        f'{screen_name}_LFC3D_neg_psig': "LFC3D_neg_psig",
        f'{screen_name}_LFC3D_pos_psig': "LFC3D_pos_psig", 
    }, inplace=True)

    # Assign p-significance label for hue coloring
    psig_dict = {'above': f'p>={lfc3d_hit_threshold}', 'below': f'p<{lfc3d_hit_threshold}'}

    def assign_psig_label(row):
        if row['LFC3D_neg_psig'] == psig_dict['above'] and row['LFC3D_pos_psig'] == psig_dict['above']:
            return 'not hit'
        elif row['LFC3D_neg_psig'] == psig_dict['above'] and row['LFC3D_pos_psig'] == psig_dict['below']:
            return 'positive hit'
        elif row['LFC3D_neg_psig'] == psig_dict['below'] and row['LFC3D_pos_psig'] == psig_dict['above']:
            return 'negative hit'
        elif row['LFC3D_neg_psig'] == psig_dict['below'] and row['LFC3D_pos_psig'] == psig_dict['below']:
            return 'pos/neg hit'
        return None
    df_input['psig_label'] = df_input.apply(assign_psig_label, axis=1)

    # Remove dashes from table and replace with 0
    df_input['LFC'] = df_input['LFC'].replace('-', 0.0).astype(float)
    df_input['LFC3D'] = df_input['LFC3D'].replace('-', 0.0).astype(float)

    y_min = df_input['LFC3D'].min()
    x_min = df_input['LFC'].min()
    df_input['LFC'] = df_input['LFC'].replace(0.0, x_min-1).astype(float)

    # Hit Type Colors
    custom_palette = {
        'not hit': 'grey', 'positive hit': 'blue',
        'negative hit': 'red', 'pos/neg hit': 'magenta'
    }
    # Scatter plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df_input, x='LFC', y='LFC3D', hue="psig_label", palette=custom_palette)
    plt.axhline(y_min, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(x_min, color="gray", linestyle="--", linewidth=0.8)
    plt.title(f"{input_gene} LFC vs LFC3D Scatter Plot")
    plt.xlabel(f"{screen_name} (LFC)")
    plt.ylabel(f"{screen_name} (LFC3D)")
    plt.grid(True, linestyle="--", alpha=0.5)

    out_filename = f'characterization/plots/{input_gene}_LFC_LFC3D_scatter.{save_type}'
    plt.savefig(working_filedir / out_filename, dpi=100, transparent=True, format=save_type)
    plt.close()

def pLDDT_RSA_scatter(
    df_filtered, 
    workdir, 
    input_gene, 
    pLDDT_col='bfactor_pLDDT', RSA_col='RSA', size_col='LFC3D_wght', direction_col='direction', 
    color_map = {'NEG': 'darkred', 'POS': 'darkblue'}, 
    save_type='png', 
):
    """
    Description
        Generate RSA vs pLDDT barplot
    """
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'characterization'):
        os.mkdir(working_filedir / 'characterization')
    if not os.path.exists(working_filedir / 'characterization/plots'):
        os.mkdir(working_filedir / 'characterization/plots')

    colors = df_filtered[direction_col].map(color_map)

    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(
        df_filtered[pLDDT_col],
        df_filtered[RSA_col],
        s=df_filtered[size_col],
        c=colors, alpha=0.7
        )

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='POS', markerfacecolor='darkred', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='NEG', markerfacecolor='darkblue', markersize=10)
    ]

    sizes = [5, 50, 95]
    for size in sizes:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', label=f'Size {size}',
                   markerfacecolor='gray', markersize=np.sqrt(size)) )

    plt.legend(handles=legend_elements, title="Legend")
    plt.xlabel('pLDDT')
    plt.ylabel('RSA')
    plt.title(f"{input_gene} RSA vs. pLDDT Scatterplot")

    out_filename = working_filedir / f'characterization/plots/{input_gene}_pLDDT_RSA_scatter.{save_type}'
    plt.savefig(out_filename, dpi=100, transparent=True, format=save_type)
    plt.close()

def hits_feature_barplot(
    df, 
    workdir, 
    input_gene, 
    category_col,
    values_cols, values_vals, value_names, 
    plot_type='Count', colors = ['darkred', 'darkblue'], 
    save_type='png', 
):
    """
    Description
        Generate hit count barplot for specified feature (ex. RSA, pLDDT, Domain, etc.)
    """
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'characterization'):
        os.mkdir(working_filedir / 'characterization')
    if not os.path.exists(working_filedir / 'characterization/plots'):
        os.mkdir(working_filedir / 'characterization/plots')
        
    assert plot_type in ['Count', 'Fraction'], "Check plot_type must be 'Count' or 'Fraction'"

    # CREATE DF WITH ORIGINAL COUNTS FOR EACH XCOL CATEGORY AND DIRECTION #
    count_data_list = []
    for col, val, name in zip(values_cols, values_vals, value_names): 
        count_data = df.groupby([category_col, col]).size().unstack(fill_value=0)
        count_data = count_data[val]
        count_data = count_data.rename(name)
        count_data_list.append(count_data)
    
    counts_df = pd.concat(count_data_list, axis=1)
    if plot_type == 'fraction': 
        counts_df = counts_df.div(counts_df.sum(axis=0), axis=1)

    # DRAW PLOT #
    if colors is None: 
        colors = plt.cm.get_cmap('tab10').colors[:len(counts_df.columns)] 
    ax = counts_df.plot(kind='bar', figsize=(6, 4), color=colors, edgecolor='black')

    # LABELS AND OUTPUT #
    plt.xlabel(category_col)
    plt.xticks(rotation=45)
    plt.ylabel(f'{plot_type} of Hits')
    plt.title(f"{input_gene} {category_col} Hit Barplot")
    plt.legend(title='Direction')

    out_filename = working_filedir / f"characterization/plots/{input_gene}_{plot_type}_{category_col}_barplot.{save_type}"
    plt.savefig(out_filename, dpi=100, transparent=True, format=save_type)
    plt.close()

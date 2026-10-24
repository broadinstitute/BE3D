"""
File: hypothesis_tests_helpers.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-23
Description: 
"""

import pandas as pd
# import matplotlib.pylab as plt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp

from matplotlib.lines import Line2D

# HYPOTHESIS 1: There is a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations and none (silent/no mutations) per screen, per gene #

def hypothesis_one(
    working_filedir, 
    df_inputs, 
    screen_names, 
    unique_genes, 
    cases, 
    controls, 
    comp_name, 
    gene_col, 
    mut_col, 
    val_col, 
    testtype, 
    col_names, 
    screen_wide_categories=None, 
): 
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('D_', 'p_')])
    col_names.extend(['num_of_cases', 'num_of_controls', 'status'])
    screen_wide = set(screen_wide_categories or [])

    df_output = pd.DataFrame(columns=col_names)
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 
            df_edits = df_input[df_input[gene_col] == current_gene]
            if df_edits.empty: # GENE NOT IN THIS SCREEN (e.g. A MOUSE GENE IN A HUMAN SCREEN), NOTHING TO TEST #
                continue
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_case = pd.DataFrame()
            for case in cases: 
                df_case = pd.concat([df_case, df_edits.loc[df_edits[mut_col]==case].reset_index(drop=True)])
            df_control = pd.DataFrame()
            for control in controls: 
                # A SCREEN-WIDE CONTROL (e.g. INTERGENIC) BELONGS TO NO GENE, SO EVERY GENE USES THE WHOLE SCREEN'S #
                df_source = df_input if control in screen_wide else df_edits
                df_control = pd.concat([df_control, df_source.loc[df_source[mut_col]==control].reset_index(drop=True)])
            new_row.extend(add_to_row(df_control, df_case, val_col, testtype))
            new_row.extend([len(df_case), len(df_control), comparison_status(df_case, df_control)])

            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            del new_row, df_case, df_control

    # SAVE FILE #
    qc_filename = f"hypothesis_qa/{testtype}_hypothesis1.tsv"
    df_output.to_csv(working_filedir / qc_filename, sep = '\t', index=False)

    return df_output

# HYPOTHESIS 2: There's a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations per gene and none (silent/no mutations) from entire screen #

def hypothesis_two(
    working_filedir, 
    df_inputs, 
    screen_names, 
    unique_genes, 
    cases, 
    controls, 
    comp_name, 
    gene_col, 
    mut_col, 
    val_col, 
    testtype, 
    col_names, 
    screen_wide_categories=None, 
): 
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('D_', 'p_')])
    col_names.extend(['num_of_cases', 'num_of_controls', 'status'])
    
    df_output = pd.DataFrame(columns=col_names)
    df_controls = list()

    # GLOBAL SILENT AND NO MUTATION #
    # EVERY CONTROL GUIDE OF THE SCREEN, WHATEVER ITS GENE -- INCLUDING ONES WITH NO GENE, LIKE INTERGENIC CONTROLS --
    # SO screen_wide_categories NEEDS NO SEPARATE HANDLING HERE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        df_control = df_input.loc[df_input[mut_col].isin(controls)].reset_index(drop=True)
        df_controls.append(df_control)
        
    # PER SCREEN PER GENE #
    for df_input, screen_name, df_control in zip(df_inputs, screen_names, df_controls): 
        for current_gene in unique_genes: 
            df_edits = df_input[df_input[gene_col] == current_gene]
            if df_edits.empty: # GENE NOT IN THIS SCREEN (e.g. A MOUSE GENE IN A HUMAN SCREEN), NOTHING TO TEST #
                continue
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_case = pd.DataFrame()
            for case in cases: 
                df_case = pd.concat([df_case, df_edits.loc[df_edits[mut_col]==case].reset_index(drop=True)])

            # df_case IS ALREADY THIS GENE'S (BUILT FROM df_edits); FILTERING IT AGAIN BY gene_col CRASHED WHEN NO CASE
            # CATEGORY EXISTED AT ALL, SINCE AN EMPTY df_case HAS NO COLUMNS #
            new_row.extend(add_to_row(df_control, df_case, val_col, testtype))
            new_row.extend([len(df_case), len(df_control), comparison_status(df_case, df_control)])
            
            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            del new_row, df_case
    del df_control

    # SAVE FILE #
    qc_filename = f"hypothesis_qa/{testtype}_hypothesis2.tsv"
    df_output.to_csv(working_filedir / qc_filename, sep = '\t', index=False)

    return df_output

def add_to_row(
    df1, 
    df2, 
    val_col, 
    function, 
): 
    if len(df1) > 0 and len(df2) > 0: 
        if function == 'KolmogorovSmirnov': 
            D, p = ks_2samp(df1[val_col], df2[val_col])
            return [D, p]
        if function == 'MannWhitney': 
            U1, p = mannwhitneyu(df1[val_col], df2[val_col], method="asymptotic")
            return [U1, p]
    # NOT TESTED: LEFT EMPTY RATHER THAN A SENTINEL, SO NO p-VALUE FILTER CAN MISTAKE IT FOR SIGNIFICANT #
    return [np.nan, np.nan]

def comparison_status(
    df_case, 
    df_control, 
): 
    """
    Why a comparison was or was not run: 'tested', 'no_cases', 'no_controls' or 'no_cases_no_controls'.
    """
    if len(df_case) and len(df_control): 
        return 'tested'
    if not len(df_case) and not len(df_control): 
        return 'no_cases_no_controls'
    return 'no_cases' if not len(df_case) else 'no_controls'

### PLOTTING FUNCTION #

def hypothesis_plot(
    working_filedir, 
    df_MW_input, 
    df_KS_input, 
    category_names, 
    cat_colname, 
    hue_colname, 
    hypothesis, 
    header, 
    save_type, 
): 

    # SETUP PLOT BY NAME (SCREEN or GENE) #
    plt.rcParams.update({'font.size': 8})
    fig, axes = plt.subplots(nrows=len(category_names), ncols=2, sharey=True, 
                             figsize=(12, 5*len(category_names)))

    # PREP DATAFRAME MW #
    df_MW_input[f"p_{header}"] = df_MW_input[f"p_{header}"].apply(negative_log_transformation)

    all_handles = []
    all_labels = []
    # DUMMY PLOT TO EXTRACT LEGEND #
    fig_dummy, ax_dummy = plt.subplots()
    scatter_dummy = sns.scatterplot(
        ax=ax_dummy, data=df_MW_input, 
        x=f"U_{header}", y=f"p_{header}", 
        hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k'
    )
    all_handles, all_labels = ax_dummy.get_legend_handles_labels()
    plt.close(fig_dummy)

    if len(category_names) == 1: axes_list = [axes[0]] # FOR ONE SCREEN/GENE #
    else: axes_list = [axes[i,0] for i in range(len(category_names))] # FOR MULTIPLE SCREEN/GENE #

    # PLOT MW #
    handles, labels = None, None
    for ax, name in zip(axes_list, category_names):
        plot1 = sns.scatterplot(
            ax=ax, 
            data=df_MW_input[df_MW_input[cat_colname]==name], 
            x=f"U_{header}", y=f"p_{header}", 
            hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05 (-log10 ≈ 1.3)')

        # REMOVE LEGEND #
        if handles is None and labels is None: handles, labels = plot1.get_legend_handles_labels()
        # TITLE AND Y AXIS #
        ax.set_ylabel(f'-log10({f"p_{header}"})')
        ax.set_title(f'Hypothesis {hypothesis}: MW {name}')

        # GRAY BACKGROUND #
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

        # LABELS #
        ax.set_xlabel('MW U-Value')
        ax.set_ylabel('MW -log(P-Value)')

    # PREP DATAFRAME KS #
    df_KS_input[f"p_{header}"] = df_KS_input[f"p_{header}"].apply(negative_log_transformation)

    if len(category_names) == 1: axes_list = [axes[1]] # FOR ONE SCREEN/GENE #
    else: axes_list = [axes[i,1] for i in range(len(category_names))] # FOR MULTIPLE SCREEN/GENE #

    for ax, name in zip(axes_list, category_names):
        plot1 = sns.scatterplot(
            ax=ax, 
            data=df_KS_input[df_KS_input[cat_colname]==name], 
            x=f"D_{header}", y=f"p_{header}", 
            hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05 (-log10 ≈ 1.3)')
        
        # TITLE AND Y AXIS #
        ax.set_title(f'Hypothesis {hypothesis}: KS {name}')

        # GRAY BACKGROUND #
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

        # LABELS #
        ax.set_xlabel('KS D-Value')
        ax.set_ylabel('KS -log(P-Value)')

    # SAVE PLOT #
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    plot_filename = f"hypothesis_qa/hypothesis{hypothesis}_scatterplot_by_{hue_colname}.{save_type}"
    plt.savefig(working_filedir / plot_filename, dpi=100, transparent=False, format=save_type)
    plt.close()
    
    # Create p=0.05 legend entry
    pval_handle = Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='p=0.05')
    all_handles.append(pval_handle)
    all_labels.append('p=0.05')

    # CREATE SEPARATE LEGEND PLOT #
    legend_fig, legend_ax = plt.subplots(figsize=(4, len(all_handles) * 0.6))
    legend_ax.axis('off')
    legend_ax.legend(all_handles, all_labels, title=hue_colname, loc='center', fontsize='medium', frameon=False)

    # SAVE LEGEND SEPARATELY #
    legend_filename = f"hypothesis_qa/hypothesis{hypothesis}_legend_by_{hue_colname}.{save_type}"
    legend_fig.savefig(working_filedir / legend_filename, dpi=100, transparent=False, format=save_type)
    plt.close()

def negative_log_transformation(
    value
):
    if pd.notna(value) and value > 0:
        return -np.log10(value)
    return value

"""
File: hypothesis_tests_helpers.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
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

# HYPOTHESIS 1: There is a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations and none (silent/no mutations) per screen, per gene #

def hypothesis_one(
    working_filedir, 
    df_inputs, screen_names, 
    unique_genes, 
    cases, controls, comp_name, 
    gene_col, mut_col, val_col, 
    testtype, 
): 
    col_names = ['screenid', 'gene_name']
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('D_', 'p_')])
    col_names.extend(['num_of_cases','num_of_controls'])

    df_output = pd.DataFrame(columns=col_names)
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 
            df_edits = df_input[df_input[gene_col] == current_gene]
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_case = pd.DataFrame()
            for case in cases: 
                df_case = pd.concat([df_case, df_edits.loc[df_edits[mut_col]==case].reset_index(drop=True)])
            df_control = pd.DataFrame()
            for control in controls: 
                df_control = pd.concat([df_control, df_edits.loc[df_edits[mut_col]==control].reset_index(drop=True)])
            new_row.extend(add_to_row(df_case, df_control, val_col, testtype))
            new_row.extend([len(df_case),len(df_control)])

            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            del new_row, df_case, df_control

    # SAVE FILE #
    qc_filename = f"hypothesis_qc/{testtype}_hypothesis1.tsv"
    df_output.to_csv(working_filedir / qc_filename, sep = '\t', index=False)

    return df_output

# HYPOTHESIS 2: There's a significant difference in the signal (LFC) #
# between knockout (nonsense/splice) mutations per gene and none (silent/no mutations) from entire screen #

def hypothesis_two(
    working_filedir, 
    df_inputs, screen_names, 
    unique_genes, 
    cases, controls, comp_name, 
    gene_col, mut_col, val_col, 
    testtype, 
): 
    col_names = ['screenid', 'gene_name']
    if testtype == 'MannWhitney': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('U_', 'p_')])
    if testtype == 'KolmogorovSmirnov': 
        col_names.extend([pref+comp for comp in [comp_name] for pref in ('D_', 'p_')])
    col_names.extend(['num_of_cases','num_of_controls'])
    
    df_output = pd.DataFrame(columns=col_names)
    df_control = pd.DataFrame()

    # GLOBAL SILENT AND NO MUTATION #
    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes:
            df_edits = df_input[df_input[gene_col] == current_gene]

            # PARSE DF FOR EACH MUT TYPE, CONCAT TO PREVIOUS GENE #
            for control in controls: 
                df_control = pd.concat([df_control, df_edits.loc[df_edits[mut_col]==control].reset_index(drop=True)])

    # PER SCREEN PER GENE #
    for df_input, screen_name in zip(df_inputs, screen_names): 
        for current_gene in unique_genes: 
            df_edits = df_input[df_input[gene_col] == current_gene]
            new_row = [screen_name, current_gene]

            # PARSE DF FOR EACH MUT TYPE #
            df_case = pd.DataFrame()
            for case in cases: 
                df_case = pd.concat([df_case, df_edits.loc[df_edits[mut_col]==case].reset_index(drop=True)])

            df_control_in = df_control[df_control[gene_col]==current_gene]
            df_case_in = df_case[df_case[gene_col]==current_gene]
            new_row.extend(add_to_row(df_control_in, df_case_in, val_col, testtype))
            new_row.extend([len(df_case),len(df_control)])
            
            # ADD NEW ROW #
            df_output.loc[len(df_output)] = new_row
            del new_row, df_case
    del df_control

    # SAVE FILE #
    qc_filename = f"hypothesis_qc/{testtype}_hypothesis2.tsv"
    df_output.to_csv(working_filedir / qc_filename, sep = '\t', index=False)

    return df_output

def add_to_row(
    df1, df2, val_col, function, 
): 
    if len(df1) > 0 and len(df2) > 0: 
        if function == 'KolmogorovSmirnov': 
            D, p = ks_2samp(df1[val_col], df2[val_col])
            return [D, p]
        if function == 'MannWhitney': 
            U1, p = mannwhitneyu(df1[val_col], df2[val_col], method="asymptotic")
            return [U1, p]
    return [-999, -999]

### PLOTTING FUNCTION #

def hypothesis_plot(
        working_filedir, 
        df_MW_input, df_KS_input, 
        category_names, cat_colname, hue_colname, 
        testtype1, testtype2, hypothesis, 
        header, save_type, 
): 

    # SETUP PLOT BY NAME (SCREEN or GENE) #
    plt.rcParams.update({'font.size': 6})
    fig, axes = plt.subplots(nrows=len(category_names), ncols=2, sharey=True, 
                             figsize=(12, 5*len(category_names)))

    # PREP DATAFRAME MW #
    df_MW_input.replace(-999, pd.NA, inplace=True)  # replace -999 with NaN
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
        plot1 = sns.scatterplot(ax=ax, 
                                data=df_MW_input[df_MW_input[cat_colname]==name], 
                                x=f"U_{header}", y=f"p_{header}", 
                                hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05 (-log10 ≈ 1.3)')
        ax.axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p=0.1 (-log10 ≈ 1.0)')

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
    df_KS_input.replace(-999, pd.NA, inplace=True)  # replace -999 with NaN
    df_KS_input[f"p_{header}"] = df_KS_input[f"p_{header}"].apply(negative_log_transformation)

    if len(category_names) == 1: axes_list = [axes[1]] # FOR ONE SCREEN/GENE #
    else: axes_list = [axes[i,1] for i in range(len(category_names))] # FOR MULTIPLE SCREEN/GENE #

    for ax, name in zip(axes_list, category_names):
        plot1 = sns.scatterplot(ax=ax, 
                                data=df_KS_input[df_KS_input[cat_colname]==name], 
                                x=f"D_{header}", y=f"p_{header}", 
                                hue=hue_colname, palette='tab20', s=100, alpha=0.7, edgecolor='k', legend=False)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05 (-log10 ≈ 1.3)')
        ax.axhline(y=-np.log10(0.1), color='blue', linestyle='--', label='p=0.1 (-log10 ≈ 1.0)')
        
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
    plot_filename = f"hypothesis_qc/hypothesis{hypothesis}_scatterplot_by_{cat_colname}.{save_type}"
    plt.savefig(working_filedir / plot_filename, dpi=100, transparent=False, format=save_type)
    plt.close()

    # CREATE SEPARATE LEGEND PLOT #
    legend_fig, legend_ax = plt.subplots(figsize=(4, len(all_handles) * 0.3))
    legend_ax.axis('off')
    legend_ax.legend(all_handles, all_labels, title=hue_colname, loc='center', fontsize='small', frameon=False)

    # SAVE LEGEND SEPARATELY #
    legend_filename = f"hypothesis_qc/hypothesis{hypothesis}_legend_by_{cat_colname}.{save_type}"
    legend_fig.savefig(working_filedir / legend_filename, dpi=100, transparent=False, format=save_type)
    plt.close()

def negative_log_transformation(value):
    if pd.notna(value) and value > 0:
        return -np.log10(value)
    return value

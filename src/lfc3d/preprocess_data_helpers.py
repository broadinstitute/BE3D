"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-08-02
Description: Translated from Notebook 3.1

"""

import re
from scipy.stats import mannwhitneyu

aa_map = {
    'ALA': 'A',  'ARG': 'R',  'ASN': 'N',  'ASP': 'D',  'CYS': 'C',  
    'GLN': 'Q',  'GLU': 'E',  'GLY': 'G',  'HIS': 'H',  'ILE': 'I',  
    'LEU': 'L',  'LYS': 'K',  'MET': 'M',  'PHE': 'F',  'PRO': 'P',  
    'SER': 'S',  'THR': 'T',  'TRP': 'W',  'TYR': 'Y',  'VAL': 'V', 
    'TER': '*', 
}

def identify_mutations(xs): 
    if not isinstance(xs, float): 
        return [x.strip() for x in xs if re.match('^[A-Z*][0-9]{1,4}[A-Z*]$', x.strip())]
    return []

# PLOT #

# def mann_whitney_test(
#     edits_filedir, screen_names, input_gene, 
# ): 
#     """
#     Description
#         Run the Mann Whitney test on 'Missense', 'Silent', 'Nonsense', 'No Mutation'
#     """

#     muts_dicts_list = [{} for _ in mut_categories_unspaced]
#     for mut, mut_dict in zip(mut_categories_unspaced, muts_dicts_list): 
#         list_mut = []
#         for screen_name in screen_names: 
#             edits_filename = f"screendata/{input_gene}_{screen_name}_{mut}.tsv"
#             try:
#                 df_temp = pd.read_csv(edits_filedir / edits_filename, sep='\t')
#                 list_mut.extend(df_temp['LFC'].tolist())
#                 del df_temp
#             except FileNotFoundError:
#                 list_mut.extend([])
#         mut_dict['LFC'] = list_mut
#         mut_dict['muttype'] = mut
#     muts_dicts = dict(map(lambda i, j : (i, j) , mut_categories_unspaced, muts_dicts_list))

#     # MANN WHITNEY TEST #
#     mannwhiteney_results = {}
#     for comp1, comp2 in comparisons: 
#         if len(muts_dicts[comp1]['LFC']) > 0 and len(muts_dicts[comp2]['LFC']) > 0: 
#             U1, p = mannwhitneyu(muts_dicts[comp1]['LFC'], muts_dicts[comp2]['LFC'], method="asymptotic")
#             mannwhiteney_results[f'{comp1} vs {comp2}'] = {'U1': U1, 'p': p}

#     df_muts = pd.concat([pd.DataFrame(d) for d in muts_dicts_list]).reset_index(drop=True)
#     return df_muts, mannwhiteney_results

# def violin_plot(
#         df_muts, edits_filedir, input_gene, screen_name, 
# ): 
#     """
#     Description
#         Graph a violin plot of LFC distribution by category
#     """

#     plt.rcParams.update({'font.size': 10})
#     fig, ax = plt.subplots(1, 2, figsize=(12,6))

#     means = df_muts.groupby('muttype')['LFC'].mean()
#     stds = df_muts.groupby('muttype')['LFC'].std()
#     medians = df_muts.groupby('muttype')['LFC'].median()
#     sorted_muttypes = sorted(means.sort_values().index.tolist())

#     sns.violinplot(ax=ax[0], data=df_muts, x="LFC", y="muttype", 
#                    order=sorted_muttypes, inner=None).set(title='LFC by Mutation Type')
#     sns.violinplot(ax=ax[1], data=df_muts, x="LFC", y="muttype", 
#                    order=sorted_muttypes, hue="LFC_direction", inner=None).set(title='LFC by Mutation Type')
#     plt.axvline(df_muts["LFC"].mean(), c="gray", linestyle="dashed")
#     plt.scatter(y=range(len(means)), x=means, c="violet", alpha=.9)
#     plt.suptitle(f'{screen_name}_{input_gene}')

#     plt.tight_layout()
#     plotname = edits_filedir / f"plots/{screen_name}_{input_gene}_LFC_dist_by_muttype.pdf"
#     plt.savefig(plotname, dpi=300, transparent=True)

#     return means, stds, medians

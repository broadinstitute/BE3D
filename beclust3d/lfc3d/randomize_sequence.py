"""
File: randomize_sequence.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: 
    Randomizes per-residue scores weighted by structural and conservation features to create a baseline distribution for significance testing.

"""

import os
import pandas as pd
from pathlib import Path

import warnings
warnings.filterwarnings('ignore')

def randomize_sequence(
    df_missense, 
    df_rand, 
    workdir, 
    input_gene, 
    screen_name, 
    nRandom=1000, 
    conservation=False, 
    muttype='Missense', 
    function_name='mean', 
    target_pos='unipos', 
    target_res=None, 
    # THERE ARE 2 MEAN FUNCTIONS, 
        # MEAN FOR CALCULATING LFC3D WHICH IS TUNABLE, 
        # AND MEAN FOR AVG RANDOMIZATIONS WHICH IS NOT TUNABLE, SO MEAN IS HARD CODED HERE #
): 
    """
    Description
    Randomizes per-residue scores weighted by structural and conservation features 
    to create a baseline distribution for significance testing.
            
    Parameters
    ----------
    df_missense : pd.DataFrame
        Per-residue LFC DataFrame from prioritize_by_sequence() to be randomized.
        
    df_rand : pd.DataFrame
        Randomized per-guide LFC DataFrame from randomize_data().

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed (e.g., 'DNMT3A', 'MEN1'). 

    screen_name : str
        Name of the screen corresponding to df_missense.

    nRandom : int, optional (default=1000)
        Number of randomizations per screen for calculating randomized LFC and LFC3D scores.
    
    conservation : bool, optional (default=False)
        If True, aggregates LFC only for residues marked as conserved. Non-conserved residues are skipped. 
        Should match whether df_consrv was provided in prioritize_by_sequence().
        
    muttype : str, optional (default='Missense')
        Mutation category to randomize (e.g., 'Missense', 'Nonsense'). 
        Should match one of the mutation types in df_dict from prioritize_by_sequence().
        
    function_name : str, optional (default='mean')
        Aggregation function name to apply per residue. 
        Should match one of the function_names used in prioritize_by_sequence().
        
    target_pos : str, optional (default='unipos')
        Column name in df_missense containing the primary sequence residue positions.
        
    target_res : str or None, optional (default=None)
        Column name in df_missense containing the primary sequence residue identities.

    Returns
    -------
    df_missense : pd.DataFrame
        DataFrame containing the randomized LFC score per residue. 
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'screendata_sequence_rand'):
        os.mkdir(working_filedir / 'screendata_sequence_rand')

    # GET RID OF NON CONSERVED FOR ALTERNATE RES ONLY #
    if target_res is not None: 
        assert target_res in df_missense.columns, 'Check [target_res]'
        res_list = df_missense[target_res].tolist()
    else: res_list = None

    # COPY VALUES FROM MISSENSE RANDOMIZED DF TO STRUCTURE SEQUENCE DF #
    position_list = df_missense[target_pos].tolist()
    rand_columns = [col for col in df_rand.columns if col.startswith('LFC')]
    new_rand_columns = [f'{function_name}_{muttype}_LFC'] + [f'{function_name}_{muttype}_LFCr{j+1}' for j in range(nRandom)]
    df_mis_positions = pd.DataFrame(columns=new_rand_columns)

    for i in range(len(df_missense)): 
        df_mis_pos = df_rand.loc[df_rand['edit_pos'] == position_list[i]]
        df_mis_pos = df_mis_pos[rand_columns]

        if df_mis_pos.shape[0] == 0 or (res_list is not None and res_list[i] == '-'): 
            # FILL WITH '-' PLACEHOLDERS IF THERE IS NOT DATA, OR IF POSITION IS NOT CONSERVED #
            res = ['-' for _ in range(df_mis_pos.shape[1])]
        else: 
            # AVERAGE ACROSS COLUMNS FOR ONE OR MORE ROWS #
            res = df_mis_pos.mean().tolist()
        df_mis_positions.loc[i] = res # ADD ROW FROM MISSENSE RANDOMIZED #

    missense_columns = ['unipos', 'unires', 'chain', 'conservation', 
                        f'{function_name}_{muttype}_LFC', f'{function_name}_{muttype}_LFC_stdev', 
                        f'{function_name}_{muttype}_LFC_Z', f'{function_name}_{muttype}_LFC_p', 
                        f'{function_name}_{muttype}_LFC_plab', f'all_{muttype}_edits', ]
    if conservation: missense_columns += ['alternate_res_pos', 'alternate_res']

    # CONSTRUCT FINAL DF FROM RELEVANT df_missense COLUMNS, AND COLLAPSED df_rand COLUMNS #
    df_missense = pd.concat([df_missense[missense_columns], df_mis_positions], axis=1)
    del df_mis_positions, df_rand

    out_filename = f"screendata_sequence_rand/{input_gene}_{screen_name}_Missense_protein_edits_rand.tsv.gz"
    df_missense.to_csv(working_filedir / out_filename, sep = '\t', index=False, compression='gzip')

    return df_missense

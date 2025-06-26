"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2

"""

import os
import numpy as np
import pandas as pd
import statistics
from pathlib import Path
from scipy.stats import norm

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def prioritize_by_sequence(
    df_dict, 
    df_struc, df_consrv, df_control, 
    workdir, 
    input_gene, screen_name, 
    pthr=0.05, 
    functions=[statistics.mean, min, max], 
    function_names=['mean', 'min', 'max'], 
    target_res_pos='human_res_pos',
    alt_res_pos='mouse_res_pos',
    alt_res='mouse_res', 
): 
    """
    Takes in results across multiple edit types for a screen, and
    aggregates the edits for each residue with sequence and conservation information. 
    
    Parameters
    ----------
    df_dict : dict of {str: pd.DataFrame}
        Dictionary mapping mutation types (e.g., 'Missense', 'Nonsense') to their respective DataFrames.
        Each DataFrame must include columns ['edit_pos', 'LFC', 'this_edit'].
            
    df_struc : pd.DataFrame
        DataFrame containing structural data for residues. Must include columns ['unipos', 'unires', 'chain'].

    df_consrv : pd.DataFrame
        DataFrame containing conservation data for residues. If None, conservation is ignored.
        Must include columns 'original_res_pos', 'alternate_res_pos', 'alternate_res', and 'conservation'.

    df_control : pd.DataFrame or None
        DataFrame of control or no-mutation LFC measurements, used to estimate background mean and std for z-scores.
        Must include 'LFC' column.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_name : str
        Name of the screens corresponding to df_missense.

    pthr : float, optional (default=0.05)
        p-value threshold for labeling statistical significance.

    functions : list of callable, optional
        List of aggregation functions to apply to LFC scores for all edits per residue (e.g., mean, min, max).

    function_names : list of str, optional
        Names corresponding to the 'functions'. Must be the same length and order as 'functions'.

    target_res_pos : str, optional (default='human_res_pos')
        Column name specifying the target residue position from df_consrv.
        
    alternate_res_pos : str, optional (default='mouse_res_pos')
        Column name specifying the alternate residue position from df_consrv.

    alternate_res : str, optional (default='mouse_res')
        Column name specifying the alternate residue information from df_consrv.

    Returns
    -------
    df_protein : pd.DataFrame
        DataFrame combining structure, conservation, aggregated LFC values, standard deviations, z-scores, p-values, and significance labels.
    """
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'screendata_sequence'):
        os.mkdir(working_filedir / 'screendata_sequence')

    # CHECK INPUTS ARE SELF CONSISTENT #
    required_columns = ['unipos', 'unires', 'chain']
    for col in required_columns: 
        assert col in df_struc.columns
    # BOTH df_consrv AND df_control SHOULD BE OPTIONAL #
    if df_consrv is not None: 
        # assert 'original_res_pos' in df_consrv.columns, 'Check [df_consrv]'
        # assert 'alternate_res_pos' in df_consrv.columns, 'Check [df_consrv]'
        # assert 'alternate_res' in df_consrv.columns, 'Check [df_consrv]'
        assert 'conservation' in df_consrv.columns, 'Check [df_consrv]'
        assert len(df_struc) == len(df_consrv)
    if df_control is not None: 
        assert 'LFC' in df_control.columns, 'Check [df_control]'
    assert len(functions) == len(function_names)

    # ADD COLUMNS FROM CONSERVATION #
    df_protein = df_struc[required_columns]
    if df_consrv is None: 
        df_protein[target_res_pos]  = df_protein['unipos']
        df_protein['conservation']      = 'None'
    else: 
        df_protein[target_res_pos]  = df_consrv[target_res_pos].to_list()
        df_protein[alt_res_pos] = df_consrv[alt_res_pos].to_list()
        df_protein[alt_res]     = df_consrv[alt_res].to_list()
        df_protein['conservation']      = df_consrv['conservation'].to_list()

    # THIS IS FOR Z-SCORE CALCULATION LATER #
    # FOR NEG POS SEPARATELY, CALC Z SCORE BASED ON THE MEAN STD PER SCREEN / GENE / DIRECTION #
    neg_mask = df_control['LFC'] < 0.0 # NEG #
    pos_mask = df_control['LFC'] > 0.0 # POS #
    df_nomut_neg = df_control.loc[neg_mask, 'LFC'] # NEG #
    df_nomut_pos = df_control.loc[pos_mask, 'LFC'] # POS #
    mu_neg, sigma_neg = df_nomut_neg.mean(), df_nomut_neg.std()
    mu_pos, sigma_pos = df_nomut_pos.mean(), df_nomut_pos.std()
    del df_nomut_neg, df_nomut_pos

    # FOR EACH EDIT TYPE, AGGREGATE LFC AND EDITS WITH CONSERVATION #
    for mut, df_edit in df_dict.items(): 

        for function, function_name in zip(functions, function_names): 
            # CALCULATE LFC SCORE PER POSITION #
            arr_LFC, arr_LFC_stdev, arr_all_edits = [], [], []
            # COLLAPSE DOWN DATA BY POSITION [edit_pos] #
            df_edit_grouped = {k: v.reset_index(drop=True) for k, v in df_edit.groupby('edit_pos')}

            # PRECOMPUTE MAPPING ORIGINAL POS AND ALT RES DICTS #
            original_res_pos_dict = df_protein[target_res_pos].to_dict()

            # FOR EACH RESIDUE #
            for i in range(len(df_protein)): 
                res_pos = int()
                structured_res = True if df_protein.iloc[i]['chain'] != '-' else False
                df_consrv_res_pos_dict,df_consrv_res_dict = dict(), dict()
                if df_consrv is not None:     
                    df_consrv_res_pos_dict = df_protein[alt_res_pos].to_dict()
                    df_consrv_res_dict = df_protein[alt_res].to_dict()
                    res_pos = df_consrv_res_pos_dict[i]
                else:
                    res_pos = original_res_pos_dict[i]
                
                # PULL DATAFRAME FOR THE CURRENT POSITION, KEEP VALUE AND EDIT #
                df_pos_edits = df_edit_grouped.get(int(res_pos), 
                                                pd.DataFrame(columns=['LFC', 'this_edit']))
                
                LFC_res, all_edits_res, stdev_res = '-', '-', '-'
                # IF df_conserv IS NOT PROVIDED, AND THERE IS A VALUE AT THAT POSITION #
                # IF df_conserv IS PROVIDED, POSITION IS CONSERVED, AND THERE IS A VALUE AT THAT POSITION #
                if structured_res:
                    if (df_consrv is None) or (df_consrv_res_dict[i] != '-'):
                        if len(df_pos_edits) > 1: 
                            score_list = df_pos_edits['LFC'].tolist()
                            LFC_res = function(score_list)
                            stdev_res = np.std(score_list)

                            pos_edits_list = df_pos_edits['this_edit'].tolist()
                            all_edits_res = ';'.join(set(pos_edits_list))
                        elif len(df_pos_edits) == 1: 
                            LFC_res = df_pos_edits.at[0, 'LFC']
                            stdev_res = 0
                            all_edits_res = df_pos_edits.at[0, 'this_edit']

                arr_LFC.append(LFC_res)
                arr_LFC_stdev.append(stdev_res)
                arr_all_edits.append(all_edits_res)
                del df_pos_edits

            df_protein[f'{function_name}_{mut}_LFC'] = arr_LFC
            df_protein[f'{function_name}_{mut}_LFC_stdev'] = arr_LFC_stdev
            df_protein[f'all_{mut}_edits'] = arr_all_edits
            del arr_LFC, arr_LFC_stdev, arr_all_edits

            # CALCULATE Z SCORE PER POSITION #
            list_z_LFC, list_p_LFC, list_plab_LFC, list_plab_thr_LFC = [], [], [], []
            LFC_raws_dict = df_protein[f'{function_name}_{mut}_LFC'].to_dict()

            for i in range(len(df_protein)):
                LFC_raw = LFC_raws_dict[i]

                if LFC_raw == '-': # or float(LFC_raw) == 0.0: 
                    z_LFC, p_LFC, plab_LFC, plab_thr = '-', 1.0, 'p=1.0', '-'
                else: 
                    LFC = float(LFC_raw)
                    if (LFC < 0.0):
                        z_LFC = statistics.NormalDist(mu=mu_neg, sigma=sigma_neg).zscore(LFC)
                        p_LFC = norm.sf(abs(z_LFC))
                        plab_LFC = get_plabel(z_LFC, direction='negative')
                        plab_thr = f'p<{str(pthr)}' if p_LFC < pthr else f'p>={str(pthr)}' ### 1 or 2 tail
                    elif (LFC > 0.0):
                        z_LFC = statistics.NormalDist(mu=mu_pos, sigma=sigma_pos).zscore(LFC)
                        p_LFC = norm.sf(abs(z_LFC))
                        plab_LFC = get_plabel(z_LFC, direction='positive')
                        plab_thr = f'p<{str(pthr)}' if p_LFC < pthr else f'p>={str(pthr)}' ### 1 or 2 tail

                list_z_LFC.append(z_LFC)
                list_p_LFC.append(p_LFC)
                list_plab_LFC.append(plab_LFC)
                list_plab_thr_LFC.append(plab_thr)

            df_protein[f'{function_name}_{mut}_LFC_Z'] = list_z_LFC
            df_protein[f'{function_name}_{mut}_LFC_p'] = list_p_LFC
            df_protein[f'{function_name}_{mut}_LFC_plab'] = list_plab_LFC
            df_protein[f'{function_name}_{mut}_LFC_plab_thr'] = list_plab_thr_LFC
            del list_z_LFC, list_p_LFC, list_plab_LFC

    strcons_edits_filename = f"screendata_sequence/{input_gene}_{screen_name}_protein_edits.tsv"
    df_protein.to_csv(working_filedir / strcons_edits_filename, sep = '\t', index=False)
    return df_protein

def get_plabel(z_LFC, direction):
    # # TWO TAIl #
    # if direction == 'negative': 
    #     thresholds = [(-3.29, '-p=0.001'), (-2.58, '-p=0.01'), 
    #                   (-1.96, '-p=0.05'), (-1.65, '-p=0.1'), (-1.0, '-p=0.3')] 
    # if direction == 'positive': 
    #     thresholds = [(3.29, '+p=0.001'), (2.58, '+p=0.01'), 
    #                   (1.96, '+p=0.05'), (1.65, '+p=0.1'), (1.0, '+p=0.3')]
    # ONE TAIL #
    if direction == 'negative': 
        thresholds = [(-3.09, '-p=0.001'), (-2.32, '-p=0.01'), 
                      (-1.64, '-p=0.05'), (-1.28, '-p=0.1'), (-0.52, '-p=0.3')] 
    if direction == 'positive': 
        thresholds = [(3.09, '+p=0.001'), (2.32, '+p=0.01'), 
                      (1.64, '+p=0.05'), (1.28, '+p=0.1'), (0.52, '+p=0.3')]

    for threshold, label in thresholds:
        if (direction == 'negative' and z_LFC < threshold) or (direction == 'positive' and z_LFC > threshold):
            return label
    return '-p=1.0' if direction=='negative' else '+p=1.0'

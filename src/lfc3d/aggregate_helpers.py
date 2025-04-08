"""
File: aggregate_helpers.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: 
"""

import numpy as np
import pandas as pd

import statistics
import scipy.stats as stats

# HELPER FUNCTIONS #

def sum_dash(values): 
    new_values = [x for x in values if x != '-']
    if len(new_values) == 0: return '-'
    else: return sum(new_values)

def filter_dash(x, mode): 
    if mode == 'neg': 
        return float(x) if x != '-' and float(x) < 0 else np.nan
    if mode == 'pos': 
        return float(x) if x != '-' and float(x) > 0 else np.nan

def calculate_stats(signal, param, pthr):
    """
    Helper function to calculate stats: z, p, plabel
    """

    # NORMALIZE ACROSS WHOLE DATASET, NOT JUST THE ROW #
    if np.isnan(signal) or param['s'] == 0: 
        return '-','-','-'
    signal_z = statistics.NormalDist(mu=param['mu'], sigma=param['s']).zscore(signal)
    signal_p = stats.norm.sf(abs(signal_z))
    signal_plabel = f'p<{str(pthr)}' if signal_p < pthr and abs(signal) > abs(param['mu']) else f'p>={str(pthr)}' ### 1 or 2 tail
    return signal_z, signal_p, signal_plabel

def binning_neg_pos(
        df_LFC_LFC3D, df_neg_stats, df_pos_stats, 
        quantile_vals, LFC3D_header
): 
    NEG_10p_v, POS_90p_v, NEG_05p_v, POS_95p_v = quantile_vals
    # BIN AND WEIGHT #
    arr_LFC3D_disc, arr_LFC3D_weight = [], []

    for i in range(0, len(df_LFC_LFC3D)): 
        LFC3D = df_LFC_LFC3D.at[i, LFC3D_header]
        if LFC3D == '-' or LFC3D == 0.0:
            LFC3D_disc, LFC3D_weight = '-', 0.0
        else: 
            LFC3Df = float(LFC3D)
            # ALIGNED FOR BETTER READABILITY #
            if                         LFC3Df <= NEG_05p_v:           LFC3D_disc, LFC3D_weight = 'NEG_05p', -0.95
            elif           NEG_05p_v < LFC3Df <= NEG_10p_v:           LFC3D_disc, LFC3D_weight = 'NEG_10p', -0.9
            elif           NEG_10p_v < LFC3Df <= df_neg_stats['25%']: LFC3D_disc, LFC3D_weight = 'NEG_25p', -0.75
            elif df_neg_stats['25%'] < LFC3Df <= df_neg_stats['50%']: LFC3D_disc, LFC3D_weight = 'NEG_50p', -0.5
            elif df_neg_stats['50%'] < LFC3Df <= df_neg_stats['75%']: LFC3D_disc, LFC3D_weight = 'NEG_75p', -0.25
            elif df_neg_stats['75%'] < LFC3Df <= df_neg_stats['max']: LFC3D_disc, LFC3D_weight = 'NEG_100p', -0.05
            
            elif df_pos_stats['25%'] > LFC3Df >= df_pos_stats['min']: LFC3D_disc, LFC3D_weight = 'POS_0p', 0.05
            elif df_pos_stats['50%'] > LFC3Df >= df_pos_stats['25%']: LFC3D_disc, LFC3D_weight = 'POS_25p', 0.25
            elif df_pos_stats['75%'] > LFC3Df >= df_pos_stats['50%']: LFC3D_disc, LFC3D_weight = 'POS_50p', 0.50
            elif           POS_90p_v > LFC3Df >= df_pos_stats['75%']: LFC3D_disc, LFC3D_weight = 'POS_75p', 0.75
            elif           POS_95p_v > LFC3Df >= POS_90p_v:           LFC3D_disc, LFC3D_weight = 'POS_90p', 0.90
            elif                       LFC3Df >= POS_95p_v:           LFC3D_disc, LFC3D_weight = 'POS_95p', 0.95
            else: LFC3D_disc, LFC3D_weight = 'NA', 0.0

        arr_LFC3D_disc.append(LFC3D_disc)
        arr_LFC3D_weight.append(LFC3D_weight)

    return arr_LFC3D_disc, arr_LFC3D_weight

def binning_lfc3d(
        df_meta, neg, pos
): 
    """
    Description
        Helper function to bin the top 10 and bottom 10 % of points
    """
    
    df_3d_list = [pd.DataFrame(), pd.DataFrame()]
    quantile_numbers = {neg: (0.1, 0.05), pos: (0.9, 0.95)}
    result = {}

    for colname, df in zip([neg, pos], df_3d_list): 
        res = {}
        df_3d_clean = df_meta.loc[df_meta[colname] != 0.0, ].reset_index(drop=True)
        df[colname] = df_3d_clean[colname].replace('-', np.nan).astype(float)

        res['dfstats'] = df[colname].describe()
        res['p1'] = df[colname].quantile(quantile_numbers[colname][0])
        res['p2'] = df[colname].quantile(quantile_numbers[colname][1])
        result[colname] = res

    df_neg_stats, df_pos_stats = result[neg]['dfstats'], result[pos]['dfstats']
    bins = [result[neg]['p1'], result[pos]['p1'], 
            result[neg]['p2'], result[pos]['p2']]

    for colname, df in zip([neg, pos], df_3d_list): 
        arr_LFC3D_disc, _ = binning_neg_pos(df_meta, df_neg_stats, df_pos_stats, bins, colname)
        df_meta[colname+'_dis'] = arr_LFC3D_disc

    return df_meta

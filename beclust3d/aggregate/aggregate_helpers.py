"""
File: aggregate_helpers.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: 
"""

import numpy as np
import pandas as pd

import statistics
import scipy.stats as stats
import glob, os

# HELPER FUNCTIONS #

def benjamini_hochberg(
    pvals
):
    """
    Benjamini-Hochberg (BH) step-up false discovery rate (FDR) q-values.

    BE3D calls per-residue significance directly on the raw p-values
    (e.g. p<0.05/0.01/0.001) with no correction for the many residues tested
    simultaneously. In practice the fraction of residues flagged at raw p<0.05
    (the "base rate") has ranged from ~7% to ~50% across targets, so the raw
    p<0.05 flag over-calls and forces users to rank hits by hand. The standard
    multiple-testing fix is to report BH/FDR q-values alongside the raw p-values
    and threshold on those instead. A q-value threshold of q<0.1 (i.e. accepting
    an expected 10% false discovery rate) is the recommended, multiplicity-aware
    way to call significant residues.

    The q-value for a p-value p_(i) (with p-values sorted ascending, rank i,
    m tested) is::

        q_(i) = min_{k >= i} ( m / k * p_(k) )

    capped at 1.0, which makes the sorted q-values monotonically non-decreasing
    and guarantees q >= p for every entry.

    Parameters
    ----------
    pvals : array-like
        Raw p-values. May contain NaN or the string '-' (BE3D's missing-value
        marker); those entries are ignored when computing the FDR correction and
        are returned as NaN in the same position. Order is preserved and ties are
        handled with a stable sort.

    Returns
    -------
    numpy.ndarray
        Float array of BH q-values, same length and order as ``pvals``, with NaN
        wherever the input was NaN / '-'.
    """
    # COERCE TO FLOAT, TREATING '-' AND None AS MISSING (NaN) #
    arr = np.array(
        [np.nan if (v is None or (isinstance(v, str) and v == '-')) else v
         for v in pvals],
        dtype=float,
    )

    q = np.full(arr.shape, np.nan, dtype=float)
    mask = ~np.isnan(arr)
    p = arr[mask]
    m = p.size
    if m == 0:
        return q

    # STEP-UP BH: SORT ASCENDING (STABLE), SCALE BY m / rank #
    order = np.argsort(p, kind='stable')
    ranked = p[order]
    ranks = np.arange(1, m + 1)
    q_sorted = ranked * m / ranks

    # ENFORCE MONOTONICITY (CUMULATIVE MIN FROM THE LARGEST P) AND CAP AT 1 #
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.minimum(q_sorted, 1.0)

    # UNSORT BACK TO ORIGINAL ORDER, THEN SCATTER INTO NON-MISSING SLOTS #
    q_valid = np.empty(m, dtype=float)
    q_valid[order] = q_sorted
    q[mask] = q_valid
    return q

def sum_dash(
    values
): 
    """
    Sum a list that contains '-' that needs to be skipped. 
    """
    new_values = [x for x in values if x != '-']
    if len(new_values) == 0: return '-'
    else: return sum(new_values)

def filter_dash(
    x, 
    mode
): 
    """
    Filter a list for only non dash and neg or non dash and pos. 
    """
    if mode == 'neg': 
        return float(x) if x != '-' and float(x) < 0 else np.nan
    if mode == 'pos': 
        return float(x) if x != '-' and float(x) > 0 else np.nan

def calculate_stats(
    signal, 
    param, 
    pthr
):
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
    df_LFC_LFC3D, 
    df_neg_stats, 
    df_pos_stats, 
    quantile_vals, 
    LFC3D_header
): 
    """
    Binning a score in df_LFC_LFC3D into NEG or POS percentiles. 
    """
    
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
            if                         LFC3Df <= NEG_05p_v:           LFC3D_disc, LFC3D_weight = 'NEG 5th', -0.95
            elif           NEG_05p_v < LFC3Df <= NEG_10p_v:           LFC3D_disc, LFC3D_weight = 'NEG 10th', -0.9
            elif           NEG_10p_v < LFC3Df <= df_neg_stats['25%']: LFC3D_disc, LFC3D_weight = 'NEG 25th', -0.75
            elif df_neg_stats['25%'] < LFC3Df <= df_neg_stats['50%']: LFC3D_disc, LFC3D_weight = 'NEG 50th', -0.5
            elif df_neg_stats['50%'] < LFC3Df <= df_neg_stats['75%']: LFC3D_disc, LFC3D_weight = 'NEG 75th', -0.25
            elif df_neg_stats['75%'] < LFC3Df <= df_neg_stats['max']: LFC3D_disc, LFC3D_weight = 'NEG 100th', -0.05
            
            elif df_pos_stats['25%'] > LFC3Df >= df_pos_stats['min']: LFC3D_disc, LFC3D_weight = 'POS 0th', 0.05
            elif df_pos_stats['50%'] > LFC3Df >= df_pos_stats['25%']: LFC3D_disc, LFC3D_weight = 'POS 25th', 0.25
            elif df_pos_stats['75%'] > LFC3Df >= df_pos_stats['50%']: LFC3D_disc, LFC3D_weight = 'POS 50th', 0.50
            elif           POS_90p_v > LFC3Df >= df_pos_stats['75%']: LFC3D_disc, LFC3D_weight = 'POS 75th', 0.75
            elif           POS_95p_v > LFC3Df >= POS_90p_v:           LFC3D_disc, LFC3D_weight = 'POS 90th', 0.90
            elif                       LFC3Df >= POS_95p_v:           LFC3D_disc, LFC3D_weight = 'POS 95th', 0.95
            else: LFC3D_disc, LFC3D_weight = 'NA', 0.0

        arr_LFC3D_disc.append(LFC3D_disc)
        arr_LFC3D_weight.append(LFC3D_weight)

    return arr_LFC3D_disc, arr_LFC3D_weight

def binning_lfc3d(
    df_meta, 
    neg, 
    pos
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

def pooled_mean_std(
    means, 
    stds, 
    ns
):
    """
    Compute pooled mean and standard deviation.
    """
    
    # Total sample size
    N = sum(ns)
    
    # Pooled mean
    pooled_mean = sum(n * mu for n, mu in zip(ns, means)) / N
    
    # Pooled variance
    pooled_var = sum(n * (s**2 + (mu - pooled_mean)**2) for n, mu, s in zip(ns, means, stds)) / N
    
    # Standard deviation
    pooled_std = np.sqrt(pooled_var)
    
    return pooled_mean, pooled_std

def mu_sigma_screens(
    workdir, 
    screen_names,
): 
    """
    """
    
    neg_stats_list = list()
    pos_stats_list = list()
    
    for screen_name in screen_names:
        neg_mean, neg_std, pos_mean, pos_std = float(), float(), float(), float()
        control_tsv = glob.glob(os.path.join(f'{workdir}/screendata/',f'*_{screen_name}_No_Mutation.tsv'))[0]
        df_control = pd.read_csv(control_tsv, sep='\t', index_col=0)        
        neg_mask = df_control['LFC'] < 0.0 # NEG #
        pos_mask = df_control['LFC'] > 0.0 # POS #
        df_nomut_neg = df_control.loc[neg_mask, 'LFC'] # NEG #
        df_nomut_pos = df_control.loc[pos_mask, 'LFC'] # POS #
        
        neg_mean, neg_std, neg_count = df_nomut_neg.mean(), df_nomut_neg.std(), df_nomut_neg.count()
        pos_mean, pos_std, pos_count = df_nomut_pos.mean(), df_nomut_pos.std(), df_nomut_pos.count()
        
        neg_stats_list.append({'mean':neg_mean,'std':neg_std,'count':neg_count})       
        pos_stats_list.append({'mean':pos_mean,'std':pos_std,'count':pos_count})
        
    return (neg_stats_list, pos_stats_list)

def mu_sigma_screens_both(
    workdir, 
    screen_names,
): 
    """
    """
    stats_list = list()
    
    for screen_name in screen_names:
        mean_all, std_all = float(), float()
        control_tsv = glob.glob(os.path.join(f'{workdir}/screendata/',f'*_{screen_name}_No_Mutation.tsv'))[0]
        df_control = pd.read_csv(control_tsv, sep='\t', index_col=0)        
        mu_all = df_control[df_control['LFC']!='-'].mean().values[0]
        sigma_all = df_control[df_control['LFC']!='-'].std().values[0]
        count_all = df_control.count()
        stats_list.append({'mean':mu_all,'std':sigma_all,'count': count_all})       
        
    return stats_list

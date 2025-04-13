"""
File: metaaggregate.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: 
"""

import os
from pathlib import Path
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

from .aggregate_helpers import *

def average_split_meta(
    df_LFC_LFC3D, 
    workdir, input_gene, screen_names, 
    nRandom=1000, 
    score_type='LFC3D', aggr_func=np.sum, aggr_func_name='SUM', 
): 
    """
    Description
        A point to meta aggregate across multiple screens or just one screen, 
        calculate signal vs background, bin this new metaaggregated signal 
        into top and bottom 10 %, and plot QC graphics
    """
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'meta-aggregate'):
        os.mkdir(working_filedir / 'meta-aggregate')

    # CHECK INPUTS ARE SELF CONSISTENT #
    assert 'unipos' in df_LFC_LFC3D.columns and 'unires' in df_LFC_LFC3D.columns and 'chain' in df_LFC_LFC3D.columns

    # INITALIZE DF #
    df_bidir_meta = df_LFC_LFC3D[['unipos', 'unires', 'chain']].copy()

    # AGGR LFC3D VALUES ACROSS SCREENS FOR EACH RESIDUE #
    # EXCLUSIVE TO META-AGGREGATE, NOT IN NON-AGGREGATE #
    # SETUP PARAMS #
    list_LFC3D_neg, list_LFC3D_pos, list_LFC3D = [], [], []
    header_scores = [f"{screen_name}_{score_type}" for screen_name in screen_names]
    for header in header_scores: 
        assert header in df_LFC_LFC3D.columns
    screen_name_dicts = [df_LFC_LFC3D[header].to_dict() for header in header_scores]

    for i in range(len(df_LFC_LFC3D)): 
        values_LFC3D_neg, values_LFC3D_pos, values_LFC3D = [], [], []

        # ADD POS AND NEG VALS SEPARATELY FOR EACH RESIDUE #
        for screen_dict in screen_name_dicts: 
            LFC3D = screen_dict[i]
            if LFC3D != '-': 
                LFC3D_value = float(LFC3D)
                if LFC3D_value < 0.0: 
                    values_LFC3D_neg.append(LFC3D_value)
                    values_LFC3D.append(LFC3D_value)
                elif LFC3D_value > 0.0: 
                    values_LFC3D_pos.append(LFC3D_value)
                    values_LFC3D.append(LFC3D_value)

        # APPLY AGGR FUNCTION FOR EVRY RESIDUE FOR ALL SCREEN #
        list_LFC3D_neg.append(aggr_func(values_LFC3D_neg) if values_LFC3D_neg else '-')
        list_LFC3D_pos.append(aggr_func(values_LFC3D_pos) if values_LFC3D_pos else '-')
        list_LFC3D.append(aggr_func(values_LFC3D) if values_LFC3D else '-')

    df_bidir_meta[f'{aggr_func_name}_{score_type}_neg'] = list_LFC3D_neg
    df_bidir_meta[f'{aggr_func_name}_{score_type}_pos'] = list_LFC3D_pos
    df_bidir_meta[f'{aggr_func_name}_{score_type}'] = list_LFC3D
    del list_LFC3D_neg, list_LFC3D_pos, list_LFC3D

    # WE RAND AND THEN SPLIT INTO POS NEG, SO THERE CAN BE NEG/POS RAND DATA FOR ROWS WITH NO NEG/POS DATA #
    # PULL RANDOMIZED DATA AND AGGR LFC3D VALUES ACROSS SCREENS FOR EACH RESIDUE #
    for n in range(nRandom): 
        new_col_neg_list, new_col_pos_list = [], []
        for screen_name in screen_names: 
            # SPLIT RAND INTO POS NEG #
            colname = f"{screen_name}_{score_type}r{str(n+1)}"
            new_col_neg = df_LFC_LFC3D[f"{colname}"].apply(lambda x: filter_dash(x, 'neg'))
            new_col_pos = df_LFC_LFC3D[f"{colname}"].apply(lambda x: filter_dash(x, 'pos'))
            new_col_neg_list.append(new_col_neg.rename(f"{colname}_neg"))
            new_col_pos_list.append(new_col_pos.rename(f"{colname}_pos"))
        df_temp = pd.concat(new_col_neg_list + new_col_pos_list, axis=1)
        
        # SUM ACROSS ALL SCREENS FOR EACH RANDOMIZATION #
        headers_neg = [f"{sn}_{score_type}r{str(n+1)}_neg" for sn in screen_names]
        headers_pos = [f"{sn}_{score_type}r{str(n+1)}_pos" for sn in screen_names]
        aggr_col_neg = df_temp[headers_neg].replace('-', np.nan).sum(axis=1)
        aggr_col_pos = df_temp[headers_pos].replace('-', np.nan).sum(axis=1)
        aggr_col_neg = aggr_col_neg.rename(f"SUM_{score_type}r{str(n+1)}_neg").replace(0.0, '-')
        aggr_col_pos = aggr_col_pos.rename(f"SUM_{score_type}r{str(n+1)}_pos").replace(0.0, '-')
        df_bidir_meta = pd.concat([df_bidir_meta, aggr_col_neg, aggr_col_pos], axis=1)
        del df_temp, new_col_neg_list, new_col_pos_list, aggr_col_neg, aggr_col_pos

    # AVG ACROSS ALL RANDOMIZATIONS #
    headers_neg = [f"SUM_{score_type}r{str(n+1)}_neg" for n in range(nRandom)]
    headers_pos = [f"SUM_{score_type}r{str(n+1)}_pos" for n in range(nRandom)]
    new_col_neg = df_bidir_meta[headers_neg].replace({'-': np.nan, 0.0: np.nan}).mean(axis=1)
    new_col_pos = df_bidir_meta[headers_pos].replace({'-': np.nan, 0.0: np.nan}).mean(axis=1)
    new_col_neg = new_col_neg.rename(f"SUM_{score_type}r_neg").replace(0.0, '-')
    new_col_pos = new_col_pos.rename(f"SUM_{score_type}r_pos").replace(0.0, '-')
    df_bidir_meta = pd.concat([df_bidir_meta, new_col_neg, new_col_pos], axis=1)
    
    # SAVE #
    out_filename_bidir = working_filedir / f"meta-aggregate/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir_meta.to_csv(out_filename_bidir, sep='\t', index=False)
    return df_bidir_meta

def bin_meta(
    df_bidir_meta, 
    workdir, input_gene, 
    score_type='LFC3D', aggr_func_name='SUM', 
): 
    """
    Description
        A point to meta aggregate across multiple screens or just one screen, 
        calculate signal vs background, bin this new metaaggregated signal 
        into top and bottom 10 %, and plot QC graphics
    """
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'meta-aggregate'):
        os.mkdir(working_filedir / 'meta-aggregate')

    # SETUP PARAMS #
    header_main = f'{aggr_func_name}_{score_type}'
    random_neg, random_pos = f'SUM_{score_type}r_neg', f'SUM_{score_type}r_pos'
    headers = [header_main, f'{header_main}_neg', f'{header_main}_pos', random_neg, random_pos]
    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}

    # CHECK INPUTS ARE SELF CONSISTENT #
    for header in headers: 
        assert header in df_bidir_meta.columns

    # INITALIZE DF #
    df_dis = df_bidir_meta[['unipos', 'unires', 'chain'] + headers].copy()

    # GENERATE THRESHOLDS FOR BINNING #
    mask_neg = df_dis[random_neg] != 0.0
    mask_pos = df_dis[random_pos] != 0.0
    df_neg_stats = df_dis[random_neg][mask_neg].describe()
    df_pos_stats = df_dis[random_pos][mask_pos].describe()

    # CALCULATE QUANTILES #
    quantile_values = {}
    for name, q in quantiles.items(): 
        df_dis_clean = df_dis[header_main].replace('-', np.nan).astype(float)
        quantile_values[name] = df_dis_clean.quantile(q)

    # CALCULATE BINS #
    arr_disc, arr_weight = binning_neg_pos(df_bidir_meta, df_neg_stats, df_pos_stats, 
                                           quantile_values.values(), header_main)
    df_dis[f"{header_main}_dis"]  = arr_disc
    df_dis[f"{header_main}_wght"] = arr_weight

    # SAVE #
    out_filename_dis = working_filedir / f"meta-aggregate/{input_gene}_{score_type}_dis_wght.tsv"
    df_dis.to_csv(out_filename_dis, sep = '\t', index=False)
    return df_dis, df_neg_stats, df_pos_stats

def znorm_meta(
    df_bidir_meta, neg_stats, pos_stats, 
    workdir, input_gene, 
    pthrs=[0.05, 0.01, 0.001], score_type='LFC3D', aggr_func_name='SUM', 
): 
    """
    Description
        A point to meta aggregate across multiple screens or just one screen, 
        calculate signal vs background, bin this new metaaggregated signal 
        into top and bottom 10 %, and plot QC graphics
    """
    # THE META-AGGREGATED RESULTS ARE Z SCORED TO THE WHOLE SET OF RANDOMIZED CONTROLS #
    # ASSUMED NEG AND POS FOR EACH #
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'meta-aggregate'):
        os.mkdir(working_filedir / 'meta-aggregate')

    # CHECK INPUTS ARE SELF CONSISTENT #
    assert f'{aggr_func_name}_{score_type}_neg' in df_bidir_meta.columns
    assert f'{aggr_func_name}_{score_type}_pos' in df_bidir_meta.columns

    assert all(isinstance(item, float) for item in pthrs), '[pthrs] must be a list of p-values'
    pthrs_str = [str(pthr).split('.')[1] for pthr in pthrs]

    # INITIALIZE DF #
    df_meta_Z = df_bidir_meta[['unipos', 'unires', 'chain']].copy()

    header_main = f'{aggr_func_name}_{score_type}'
    df_meta_Z[f'{header_main}_neg'] = df_bidir_meta[f'{header_main}_neg']
    df_meta_Z[f'{header_main}_pos'] = df_bidir_meta[f'{header_main}_pos']
    df_meta_Z[f'{header_main}r_neg'] = df_bidir_meta[f'{header_main}r_neg']
    df_meta_Z[f'{header_main}r_pos'] = df_bidir_meta[f'{header_main}r_pos']

    # SETUP PARAMS FOR CALCULATING Z SCORE #
    colnames = [f'{header_main}_neg', f'{header_main}_pos']
    params = [{'mu': neg_stats['mean'], 's': neg_stats['std']}, 
              {'mu': pos_stats['mean'], 's': pos_stats['std']}, ]

    result_data = {f'{header_main}_{sign}_{pthr_str}_{suffix}': [] 
                   for sign in ['neg', 'pos'] 
                   for suffix in ['z', 'p', 'psig'] 
                   for pthr_str in pthrs_str
                   }

    # CONVERT SIGNAL TO Z SCORE FOR META-AGGREGATED RESULTS #
    for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
        signals_dict = df_meta_Z[colname].replace('-', np.nan).to_dict()

        for pthr, pthr_str in zip(pthrs, pthrs_str): 
            for i in range(len(df_meta_Z)):
                signal = float(signals_dict[i])
                signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
                
                # APPEND RESULTS TO DICT #
                result_data[f'{header_main}_{sign}_{pthr_str}_z'].append(signal_z)
                result_data[f'{header_main}_{sign}_{pthr_str}_p'].append(signal_p)
                result_data[f'{header_main}_{sign}_{pthr_str}_psig'].append(signal_plabel)

    df_temp = pd.DataFrame(result_data).replace(0,'-')
    df_meta_Z = pd.concat([df_meta_Z, df_temp], axis=1)

    # SAVE #
    filename = working_filedir / f"meta-aggregate/{input_gene}_MetaAggr_{score_type}.tsv"
    df_meta_Z.to_csv(filename, "\t", index=False)
    return df_meta_Z

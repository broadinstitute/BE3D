"""
File: nonaggregate.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import os
from pathlib import Path
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

from .aggregate_helpers import *

def average_split_score(
        df_LFC_LFC3D, 
        workdir, input_gene, screen_names, 
        score_type='LFC3D', 
): 
    """
    Splits LFC or LFC3D scores into positive and negative components and aggregates randomized scores.

    Parameters
    ----------
    df_LFC_LFC3D : pd.DataFrame
        DataFrame containing per-residue mutation scores (e.g., LFC3D), along with randomization scores.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in df_edits_list and df_rand_list.

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed (e.g., 'LFC3D', 'LFC', etc.).

    Returns
    -------
    df_bidir : pd.DataFrame
        DataFrame containing split positive/negative scores and randomized averages for each screen.
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / score_type):
        os.mkdir(working_filedir / score_type)

    # CHECK INPUTS ARE SELF CONSISTENT #
    assert 'unipos' in df_LFC_LFC3D.columns and 'unires' in df_LFC_LFC3D.columns

    # INITALIZE DF #
    df_bidir = df_LFC_LFC3D[['unipos', 'unires', 'chain']]
    
    # SETUP PARAMS #
    header_scores = [f"{screen_name}_{score_type}" for screen_name in screen_names]

    # FOR EVERY SCREEN INDIVIDUALLY #
    for screen_name, header in zip(screen_names, header_scores): 
        df_bidir[header] = df_LFC_LFC3D[header] # LFC or LFC3D per screen
        taa_wise_LFC3D_pos, taa_wise_LFC3D_neg = [], []
        
        # PULL THE SCORE AND SPLIT POS NEG #
        taa_LFC3D_raws_dict = df_LFC_LFC3D[header].to_dict()
        for aa in range(len(df_LFC_LFC3D)): 
            taa_LFC3D_raw = taa_LFC3D_raws_dict[aa] # TARGET SCORE PER AA PER SCREEN 
            
            # SEPARATE INTO POS AND NEG #
            taa_LFC3D = float(taa_LFC3D_raw) if taa_LFC3D_raw != '-' else 0.0
            taa_wise_LFC3D_neg.append(taa_LFC3D if taa_LFC3D < 0 else '-') # EITHER THE VALUE OR 0.0
            taa_wise_LFC3D_pos.append(taa_LFC3D if taa_LFC3D > 0 else '-') # EITHER THE VALUE OR 0.0

        df_bidir[f"{header}_neg"] = taa_wise_LFC3D_neg # LFC3D_neg per SCREEN
        df_bidir[f"{header}_pos"] = taa_wise_LFC3D_pos # LFC3D_pos per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r"]     = df_LFC_LFC3D[f"{screen_name}_AVG_{score_type}r"] # AVG_LFC3Dr per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r_neg"] = df_LFC_LFC3D[f"{screen_name}_AVG_{score_type}r_neg"] # AVG_LFC3Dr_neg per SCREEN
        df_bidir[f"{screen_name}_AVG_{score_type}r_pos"] = df_LFC_LFC3D[f"{screen_name}_AVG_{score_type}r_pos"] # AVG_LFC3Dr_pos per SCREEN

    # SAVE #
    out_filename_bidir = working_filedir / f"{score_type}/{input_gene}_{score_type}_bidirectional.tsv"
    df_bidir.to_csv(out_filename_bidir, sep='\t', index=False)
    return df_bidir

def bin_score(
        df_bidir, 
        workdir, input_gene, screen_names, 
        score_type='LFC3D', 
): 
    """
    Bins positive and negative LFC or LFC3D scores into percentile thresholds.

    Parameters
    ----------
    df_bidir : pd.DataFrame
        DataFrame containing split positive/negative scores and randomized averages for each screen.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in df_edits_list and df_rand_list.

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed (e.g., 'LFC3D', 'LFC', etc.).

    Returns
    -------
    df_dis : pd.DataFrame
        DataFrame containing percentile bins and weighted scores for each residue and screen.

    df_neg_stats_list : list of pd.Series
        List containing descriptive statistics for negative scores in each screen.

    df_pos_stats_list : list of pd.Series
        List containing descriptive statistics for positive scores in each screen.
    """
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / score_type):
        os.mkdir(working_filedir / score_type)

    # SETUP PARAMS #
    quantiles = {'NEG_10p_v':0.1, 'POS_90p_v':0.9, 'NEG_05p_v':0.05, 'POS_95p_v':0.95}
    headers_LFC3D = [f"{screen_name}_{score_type}" for screen_name in screen_names]
    df_neg_stats_list, df_pos_stats_list = [], []
    
    # INITALIZE DF #
    df_dis = df_bidir[['unipos', 'unires', 'chain'] + headers_LFC3D].copy()
    
    # FOR EVERY SCREEN INDIVIDUALLY #
    for screen_name, header_LFC3D in zip(screen_names, headers_LFC3D): 

        # GENERATE THRESHOLDS FOR BINNING #
        df_temp = df_dis[header_LFC3D].replace('-', np.nan).astype(float)
        mask_neg = df_temp < 0.0
        mask_pos = df_temp > 0.0
        df_neg_stats = df_temp[mask_neg].describe()
        df_pos_stats = df_temp[mask_pos].describe()
        df_neg_stats_list.append(df_neg_stats)
        df_pos_stats_list.append(df_pos_stats)

        # CALCULATE QUANTILES #
        quantile_values = {}
        for name, q in quantiles.items(): 
            df_dis_clean = df_dis[header_LFC3D].replace('-', np.nan).astype(float)
            quantile_values[name] = df_dis_clean.quantile(q)

        # CALCULATE BINS #
        arr_disc, arr_weight = binning_neg_pos(df_bidir, df_neg_stats, df_pos_stats, 
                                               quantile_values.values(), header_LFC3D)
        df_dis[f"{screen_name}_{score_type}_dis"]  = arr_disc
        df_dis[f"{screen_name}_{score_type}_wght"] = arr_weight

    # SAVE #
    out_filename_dis = working_filedir / f"{score_type}/{input_gene}_{score_type}_dis_wght.tsv"
    df_dis.to_csv(out_filename_dis, sep = '\t', index=False)
    return df_dis, df_neg_stats_list, df_pos_stats_list

def znorm_score(
        df_bidir, workdir, input_gene, screen_names, 
        score_type='LFC3D',
        pthrs=[0.05, 0.01, 0.001], 
): 
    """
    Z-normalizes scores against randomized control distributions and assigns significance labels.

    Parameters
    ----------
    df_dis : pd.DataFrame
        DataFrame containing percentile bins and weighted scores for each residue and screen.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in df_edits_list and df_rand_list.

    score_type : str, optional (default='LFC3D')
        Label for the type of mutation score analyzed (e.g., 'LFC3D', 'LFC', etc.).

    pthrs : list of float, optional
        List of p-value thresholds used to define significance (default [0.05, 0.01, 0.001]).

    Returns
    -------
    df_z : pd.DataFrame
        DataFrame containing z-scores, p-values, and significance labels for scores at multiple thresholds.
    """
    # EACH SCREEN IS Z SCORED TO ITS OWN SET OF RANDOMIZED CONTROLS #
    # ASSUMED NEG AND POS FOR EACH #
    
    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / score_type):
        os.mkdir(working_filedir / score_type)

    # CHECK INPUTS ARE SELF CONSISTENT #
    for screen_name in screen_names: 
        assert f'{screen_name}_{score_type}_neg' in df_bidir.columns
        assert f'{screen_name}_{score_type}_pos' in df_bidir.columns
        assert f'{screen_name}_AVG_{score_type}r_neg' in df_bidir.columns
        assert f'{screen_name}_AVG_{score_type}r_pos' in df_bidir.columns
    
    assert all(isinstance(item, float) for item in pthrs), '[pthrs] must be a list of p-values'
    pthrs_str = [str(pthr).split('.')[1] for pthr in pthrs]

    # INITIALIZE DF #
    df_z = df_bidir[['unipos', 'unires', 'chain']].copy()
   
    # FOR EVERY SCREEN INDIVIDUALLY #
    for idx, screen_name in enumerate(screen_names):
        neg_mean, neg_std, pos_mean, pos_std = float(),float(),float(),float()

        # COPY COLUMNS FOR VALUES AND RAND VALUES #
        header_main = f'{screen_name}_{score_type}'
        df_z[f'{header_main}_neg'] = df_bidir[f'{header_main}_neg']
        df_z[f'{header_main}_pos'] = df_bidir[f'{header_main}_pos']
        df_z[f'{screen_name}_AVG_{score_type}r_neg'] = df_bidir[f'{screen_name}_AVG_{score_type}r_neg']
        df_z[f'{screen_name}_AVG_{score_type}r_pos'] = df_bidir[f'{screen_name}_AVG_{score_type}r_pos']
        
        if score_type == 'LFC':
            neg_stats_list, pos_stats_list = mu_sigma_screens(workdir,screen_names)            
            neg_mean = neg_stats_list[idx]['mean']
            neg_std = neg_stats_list[idx]['std']
            pos_mean = pos_stats_list[idx]['mean']
            pos_std = pos_stats_list[idx]['std']
                        
        else:
            _temp_neg_pd = df_z[df_z[f'{screen_name}_AVG_{score_type}r_neg'] != '-']
            _temp_pos_pd = df_z[df_z[f'{screen_name}_AVG_{score_type}r_pos'] != '-']
            avgr_neg_list = _temp_neg_pd[f'{screen_name}_AVG_{score_type}r_neg'].to_list()
            avgr_pos_list = _temp_pos_pd[f'{screen_name}_AVG_{score_type}r_pos'].to_list()
            neg_mean = np.mean(avgr_neg_list)
            neg_std = np.std(avgr_neg_list)
            pos_mean = np.mean(avgr_pos_list)
            pos_std = np.std(avgr_pos_list)

        # SETUP PARAMS FOR CALCULATING Z SCORE #
        colnames = [f'{header_main}_neg', f'{header_main}_pos']

        params = [{'mu': neg_mean, 's': neg_std}, 
                  {'mu': pos_mean, 's': pos_std}]

        result_data = {f'{header_main}_{sign}_{pthr_str}_{suffix}': [] 
                    for sign in ['neg', 'pos'] 
                    for suffix in ['z', 'p', 'psig'] 
                    for pthr_str in pthrs_str
                    }

        # CONVERT SIGNAL TO Z SCORE #
        for colname, param, sign in zip(colnames, params, ['neg', 'pos']): 
            signals_dict = df_z[colname].replace('-', np.nan).to_dict()

            for pthr, pthr_str in zip(pthrs, pthrs_str): 
                for i in range(len(df_z)):
                    signal = float(signals_dict[i])
                    signal_z, signal_p, signal_plabel = calculate_stats(signal, param, pthr)
                    
                    # APPEND RESULTS TO DICT #
                    result_data[f'{header_main}_{sign}_{pthr_str}_z'].append(signal_z)
                    result_data[f'{header_main}_{sign}_{pthr_str}_p'].append(signal_p)
                    result_data[f'{header_main}_{sign}_{pthr_str}_psig'].append(signal_plabel)

        df_temp = pd.DataFrame(result_data)
        df_z = pd.concat([df_z, df_temp], axis=1)

    # SAVE #
    filename = working_filedir / f"{score_type}/{input_gene}_NonAggr_{score_type}.tsv"
    df_z.to_csv(filename, "\t", index=False)
    return df_z

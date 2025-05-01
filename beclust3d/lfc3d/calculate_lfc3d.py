"""
File: calculate_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3

"""

import pandas as pd
import numpy as np
from pathlib import Path
import os
import warnings
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def calculate_lfc3d(
        df_str_cons, df_edits_list, df_rand_list, 
        workdir, input_gene, screen_names, 
        nRandom=1000, muttype='Missense', function_type='mean', function_aggr=np.mean, 
        LFC_only=False, conserved_only=False, 
        # THERE ARE 2 MEAN FUNCTIONS, MEAN FOR CALCULATING LFC3D WHICH IS TUNABLE, AND MEAN FOR AVG RANDOMIZATIONS WHICH IS NOT TUNABLE #
): 
    """
    Calculates LFC3D scores using structural data. 

    Parameters
    ----------
    df_str_cons : pd.DataFrame
        DataFrame containing structural conservation data for residues. 
        Must include columns 'unipos', 'unires', 'chain', 'Naa_pos', 'Naa_chain'.

    df_edits_list : list of pd.DataFrame
        List of mutation DataFrames for each screen. 

    df_rand_list : list of pd.DataFrame
        List of randomized mutation DataFrames for each screen.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_names : list of str
        Names of the different screens corresponding to each DataFrame in df_edits_list and df_rand_list.

    nRandom : int, optional (default=1000)
        Number of randomizations per screen for calculating randomized LFC and LFC3D scores.

    muttype : str, optional (default='Missense')
        Type of mutation to focus on (e.g., 'Missense', 'Nonsense', etc.).

    function_type : str, optional (default='mean')
        String label for the type of aggregation function used to compute LFC3D scores.

    function_aggr : function, optional (default=np.mean)
        Aggregation function used to summarize neighboring mutation effects when computing LFC3D scores.
        Function should take a list or array of values and return a scalar (e.g., np.mean, np.median).

    LFC_only : bool, optional (default=False)
        If True, skips the LFC3D computation.

    conserved_only : bool, optional (default=False)
        If True, calculates LFC3D only for residues marked as 'conserved' in the conservation data.
        Non-conserved residues will be skipped (set to NaN or '-').

    Returns
    -------
    df_struct_3d : pd.DataFrame
        DataFrame containing the structural data, LFC, LFC3D, and randomized scores. 
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'LFC3D'):
        os.mkdir(working_filedir / 'LFC3D')

    # CHECK INPUTS ARE SELF CONSISTENT #
    for str_cons_df, str_cons_rand_df in zip(df_edits_list, df_rand_list): 
        assert len(df_str_cons) == len(str_cons_df) == len(str_cons_rand_df)
    assert 'unipos' in df_str_cons.columns and 'unires' in df_str_cons.columns and 'chain' in df_str_cons.columns
    # COLUMNS FOR SMOOTHING ACROSS RESIDUES AND#
    assert 'Naa_pos' in df_str_cons.columns
    assert 'Naa_chain' in df_str_cons.columns
    structure_columns = ['Naa_pos', 'Naa_chain']
    core_columns = ['unipos', 'unires', 'chain']

    assert len(df_edits_list) == len(df_rand_list) == len(screen_names)

    df_struct_3d = df_str_cons[core_columns + structure_columns].copy()
    naa_pos_dict = df_str_cons['Naa_pos'].to_dict()
    naa_pos_chain_dict = {
        (row['unipos'], row['chain']) : (row['Naa_pos'], row['Naa_chain'])
        for _, row in df_str_cons.iterrows()
    }

    # FOR EVERY SCREEN #
    for screen_name, df_edits, df_rand in zip(screen_names, df_edits_list, df_rand_list):

        taa_conserv_dict = df_edits['conservation'].to_dict() ###
        # ADD LFC COLUMNS FROM DF #
        lfc_colname = f'{function_type}_{muttype}_LFC'
        df_struct_3d = pd.concat([df_struct_3d, 
                                  df_edits[[lfc_colname]].rename(columns={lfc_colname: f"{screen_name}_LFC"}), 
                                  df_edits[[f'{lfc_colname}_Z']].rename(columns={f'{lfc_colname}_Z': f"{screen_name}_LFC_Z"})], axis=1)
        
        # CALCULATE LFC3D, IF LFC_only SKIP OVER #
        if not LFC_only: 
            aggr_vals = []
            taa_LFC_dict = df_edits[lfc_colname].to_dict() ###

            for aa in range(len(df_edits)): # FOR EVERY RESIDUE # ###
                # RESIDUE NEEDS TO BE CONSERVED #
                if conserved_only and taa_conserv_dict[aa] != 'conserved':  ###
                    aggr_vals.append('-')
                    continue
                # CALCULATE LFC3D #
                taa_naa_LFC_vals = helper(aa, taa_LFC_dict, taa_conserv_dict, naa_pos_dict[aa], conserved_only) ###
                if len(taa_naa_LFC_vals) == 0:
                    aggr_vals.append('-')
                else: 
                    aggr_vals.append(str(function_aggr(taa_naa_LFC_vals)))
            
            df_struct_3d = pd.concat([df_struct_3d, pd.DataFrame({f"{screen_name}_LFC3D": aggr_vals})], axis=1)
            del taa_LFC_dict, aggr_vals

        # REPEAT LFC LFC3D CALCULATIONS FOR RANDOMIZED DATA #
        
        dict_temp = {}
        for r in range(nRandom):
            # ADD LFC RANDOMIZATION COLUMNS FROM DF #
            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = df_rand[f'{lfc_colname}r{str(r+1)}']

            # CALCULATE LFC3D RANDOMIZATION, IF LFC_only SKIP OVER #
            if not LFC_only: 
                aggr_vals = []
                taa_LFC_rand_dict = df_rand[f'{lfc_colname}r{str(r+1)}'].to_dict() ###

                for aa in range(len(df_rand)): # FOR EVERY RESIDUE # ###
                    # RESIDUE NEEDS TO BE CONSERVED #
                    if conserved_only and taa_conserv_dict[aa] != 'conserved': ###
                        aggr_vals.append('-')
                        continue
                    # CALCULATE LFC3D RANDOMIZATION #
                    taa_naa_LFC_vals = helper(aa, taa_LFC_rand_dict, taa_conserv_dict, naa_pos_dict[aa], conserved_only) ###
                    if len(taa_naa_LFC_vals) == 0:
                        aggr_vals.append('-')
                    else:
                        aggr_vals.append(function_aggr(taa_naa_LFC_vals))
                
                dict_temp[f"{screen_name}_LFC3Dr{str(r+1)}"] = aggr_vals
                del aggr_vals

        df_struct_3d = pd.concat((df_struct_3d, pd.DataFrame(dict_temp)), axis=1)
        # CONVERT '-' TO NAN FOR EASIER CALCULATIONS #
        df_struct_3d = df_struct_3d.replace('-', np.nan).infer_objects(copy=False)
        df_struct_3d = df_struct_3d.apply(lambda col: pd.to_numeric(col, errors='coerce'))
        del dict_temp

        # AVG OVER LFC RANDOMIZATION COLUMNS FROM DF #
        LFC_colnames = [f"{screen_name}_LFCr{str(r+1)}" for r in range(nRandom)]
        df_struct_3d[f"{screen_name}_AVG_LFCr"]     = df_struct_3d[LFC_colnames].mean(axis=1) # AVG ALL
        df_struct_3d[f"{screen_name}_AVG_LFCr_neg"] = (df_struct_3d[LFC_colnames]
                                                       .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                       .sum(axis=1) / nRandom) # AVG NEG
        df_struct_3d[f"{screen_name}_AVG_LFCr_pos"] = (df_struct_3d[LFC_colnames]
                                                       .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                       .sum(axis=1) / nRandom) # AVG POS
        
        # AVG OVER LFC3D RANDOMIZATION COLUMNS FROM DF #
        if not LFC_only: 
            LFC3D_colnames = [f"{screen_name}_LFC3Dr{str(r+1)}" for r in range(nRandom)]
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr"]     = df_struct_3d[LFC3D_colnames].mean(axis=1) # AVG ALL
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr_neg"] = (df_struct_3d[LFC3D_colnames]
                                                            .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                            .sum(axis=1) / nRandom) # AVG NEG
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr_pos"] = (df_struct_3d[LFC3D_colnames]
                                                            .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                            .sum(axis=1) / nRandom) # AVG POS
        
        # CONVERT NAN TO '-' FOR REPRESENTATION #
        df_struct_3d = df_struct_3d.fillna('-')
        print('Calculated LFC3D for', screen_name)

    df_struct_3d[core_columns + structure_columns] = df_str_cons[core_columns + structure_columns]
    out_filename = working_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr.tsv"
    df_struct_3d.to_csv(out_filename, sep = '\t', index=False)

    return df_struct_3d

def helper(
    aa, taa_LFC_dict, df_struc_edits_dict, naa_pos_str, conserved_only
): 
    # naa IS NEIGHBORING AMINO ACIDS #
    # taa IS THIS AMINO ACID #
    taa_naa_LFC_vals = []
    taa_LFC = taa_LFC_dict[aa] ###

    # VALUE FOR THIS RESIDUE #
    if taa_LFC != '-': 
        if not conserved_only or df_struc_edits_dict[aa] == 'conserved': ###
            taa_naa_LFC_vals.append(float(taa_LFC))

    # CHECK NEIGHBORING RESIDUES #
    if isinstance(naa_pos_str, str):  ###
        naa_pos_list = naa_pos_str.split(';') ###
        for naa_pos in naa_pos_list:  ###
            # VALUE FOR NEIGHBORING RESIDUE #
            if not conserved_only or df_struc_edits_dict[int(naa_pos)-1] == 'conserved': ###
                naa_LFC = taa_LFC_dict[int(naa_pos)-1] ###
                if naa_LFC != '-': 
                    taa_naa_LFC_vals.append(float(naa_LFC))

    return taa_naa_LFC_vals

"""
File: randomize_data.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: 
"""

import pandas as pd
from pathlib import Path
import numpy as np
import os

def randomize_data(
        df_missense, 
        workdir, input_gene, 
        screen_name, 
        nRandom=1000, val_colname = 'LFC', muttype='Missense', 
        seed=False, 
): 
    """
    Takes reformatted Missense dataframe and randomizes them to provide a baseline signal. 
    This function is run per gene per screen. 

    Parameters
    ----------
    df_missense : pd.DataFrame
        DataFrame from preprocess_data which contains the LFC per guide information to be randomized. 

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    screen_name : str
        Name of the screens corresponding to df_missense.

    nRandom : int, optional (default=1000)
        Number of randomizations per screen for calculating randomized LFC and LFC3D scores.

    val_colname : str, optional (default='LFC')
        Column name in df_missense specifying the value measurement.

    muttype : str, optional (default='Missense')
        Type of mutation to focus on (e.g., 'Missense', 'Nonsense', etc.).

    seed : bool, optional (default=False)
        Whether or not to randomize with a seed. 

    Returns
    -------
    df_missense : pd.DataFrame
        DataFrame containing the randomized LFC score per guide. 
    """

    # NAME VARIABLES, PATHS, CREATE DIRECTORIES #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'screendata_rand'):
        os.mkdir(working_filedir / 'screendata_rand')

    assert val_colname in df_missense.columns, 'Check [mut_col] input'

    # SHUFFLE AND ADD TO DATAFRAME #
    LFC_list = df_missense[val_colname].tolist()
    dict_temp = {}

    for i in range(nRandom): 
        if seed: 
            rng = np.random.default_rng(i)
            dict_temp[f"{val_colname}r{str(i+1)}"] = rng.permutation(LFC_list)
        else: 
            dict_temp[f"{val_colname}r{str(i+1)}"] = np.random.permutation(LFC_list)
    df_missense = pd.concat((df_missense, pd.DataFrame(dict_temp)), axis=1)

    # SAVE RESULTS #
    screen_name_nospace = screen_name.replace(' ','_')
    out_filename = working_filedir / f"screendata_rand/{input_gene}_{screen_name_nospace}_{muttype}_rand.tsv"
    df_missense.to_csv(out_filename, sep='\t', index=False)
    
    return df_missense

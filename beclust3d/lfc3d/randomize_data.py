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
    Description
        Takes reformatted Missense dataframe and randomizes them to provide a baseline signal. 
        This function is run per gene per screen. 
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

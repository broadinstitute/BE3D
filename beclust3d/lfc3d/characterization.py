"""
File: characterization.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-23
Description: 
"""

import os
import csv
import requests
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt

import scipy.stats as stats
from scipy.stats import fisher_exact
from pathlib import Path

def enrichment_test(
        df, 
        workdir, input_gene, 
        hit_columns, hit_threshold, # HITS DESCRIBED BY P-VALUE
        feature_column, feature_values, # FEATURES SUCH AS DOMAINS, pLDDT
        confidence_level=0.95, 
): 
    """
    Description
        Perform enrichment tests (e.g., Fisher's exact test) for structural features and plots the results

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing residue-level mutation data, significance scores, and features.

    workdir : str
        Path to the working directory where output files and results will be saved.

    input_gene : str
        Name of the gene being processed. 

    hit_columns : list of str
        List of column names representing hit scores (e.g., p-values, LFC scores) to test for enrichment.

    hit_threshold : float
        Threshold on the hit scores used to define "hits" (e.g., p < 0.05).

    feature_column : str
        Column name in `df` representing the structural or functional feature of interest 
        (e.g., domain assignment, pLDDT score category, secondary structure).

    feature_values : list of str
        List of specific feature values to test for enrichment (e.g., domain names, "Low pLDDT").

    confidence_level : float, optional (default=0.95)
        Confidence level for calculating confidence intervals on the odds ratios.
        
    save_type : str, optional (default='png')
        Format for saving output plots (e.g., 'png', 'pdf').

    Returns
    -------
    results : list of dict
        A list where each element is a dictionary containing:
            - 'score_type': Name of the hit column tested.
            - 'odds_ratio': Estimated odds ratio for enrichment.
            - 'ci': Confidence interval tuple (lower, upper) for the odds ratio.
            - 'p_value': P-value from the enrichment test.
    """

    # MKDIR #
    working_filedir = Path(workdir)
    if not os.path.exists(working_filedir): 
        os.mkdir(working_filedir)
    if not os.path.exists(working_filedir / 'characterization'):
        os.mkdir(working_filedir / 'characterization')

    # CHECK INPUTS ARE SELF CONSISTENT #
    assert feature_column in df.columns, 'Check [feature_column] input'

    for hit_column in hit_columns: 
        assert hit_column in df.columns, f'Check {hit_column} input'
    
    for feature_value in feature_values: 
        if feature_value not in df[feature_column].unique(): 
            warnings.warn(f'{feature_value} not found in df')
    
    # ENRICHMENT TEST #
    results = []

    for hit_col in hit_columns: 
        try: 
            test1, odds1, ci1 = category_test(data=df, 
                                              hit_col=hit_col, hit_threshold=hit_threshold, 
                                              feature_col=feature_column, feature_values=feature_values, 
                                              confidence_level=confidence_level)
            results.append({
                'score_type': hit_col, 
                'odds_ratio': test1[0], 'ci': ci1, 'p_value': test1[1]
            })
        except (ValueError, AttributeError):
            results.append({
                'score_type': hit_col, 
                'odds_ratio': np.nan, 'ci': None, 'p_value': np.nan
            })

    out_filename = working_filedir / f"characterization/{input_gene}_enrichment_test.pickle"
    with open(out_filename, 'wb') as f: 
        pickle.dump(results, f)
    return results

def category_test(data, hit_col, hit_threshold, 
                  feature_col, feature_values, 
                  confidence_level):
    
    # SEPARATE DATA INTO ABOVE AND BELOW A THRESHOLD #
    below = data[data[hit_col].astype(float) < hit_threshold]
    above = data[data[hit_col].astype(float) >= hit_threshold]

    # DEFINE IN AND OUT CATEGORIES #
    in_list = feature_values
    out_list = [x for x in data[feature_col].unique().tolist() if not x in feature_values]

    table = np.array([
        [len(below[below[feature_col].isin(in_list)]),   # BELOW THRESHOLD, IN #
         len(below[below[feature_col].isin(out_list)])], # BELOW THRESHOLD, OUT #
        [len(above[above[feature_col].isin(in_list)]),   # ABOVE THRESHOLD, IN #
         len(above[above[feature_col].isin(out_list)])]  # ABOVE THRESHOLD, OUT #
    ])

    ftest = fisher_exact(table, alternative='two-sided')
    odds_ratio = stats.contingency.odds_ratio(table)
    ci = odds_ratio.confidence_interval(confidence_level=confidence_level)
    return ftest, odds_ratio, ci

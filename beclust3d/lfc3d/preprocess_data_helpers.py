"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-08-02
Description: Translated from Notebook 3.1
"""

import re
from scipy.stats import mannwhitneyu

aa_map = {
    'ALA': 'A',  'ARG': 'R',  'ASN': 'N',  'ASP': 'D',  'CYS': 'C',  
    'GLN': 'Q',  'GLU': 'E',  'GLY': 'G',  'HIS': 'H',  'ILE': 'I',  
    'LEU': 'L',  'LYS': 'K',  'MET': 'M',  'PHE': 'F',  'PRO': 'P',  
    'SER': 'S',  'THR': 'T',  'TRP': 'W',  'TYR': 'Y',  'VAL': 'V', 
    'TER': '*', 
}

def identify_mutations(
    xs, 
): 
    if not isinstance(xs, float): 
        return [x.strip() for x in xs if re.match('^[A-Z*][0-9]{1,4}[A-Z*]$', x.strip())]
    return []

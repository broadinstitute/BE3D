"""
File: preprocess_data_helpers.py
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

def reduce_mutation_type(
    mut_type,
    mut_delimiter,
    priority_order,
):
    """
    Collapse a delimiter-joined multi-category mutation type string
    (e.g. 'Silent;Missense;', one category per edit in the guide) into a
    single category, keeping whichever category appears first in
    priority_order (most to least deleterious). Categories not listed in
    priority_order fall back to whichever token appears first in mut_type.
    Single-category values (no delimiter, or all tokens identical) are
    returned unchanged.
    """
    if not isinstance(mut_type, str):
        return mut_type
    tokens = [t.strip() for t in mut_type.split(mut_delimiter) if t.strip()]
    if not tokens:
        return mut_type
    if len(set(tokens)) == 1:
        return tokens[0]
    for category in priority_order:
        if category in tokens:
            return category
    return tokens[0]

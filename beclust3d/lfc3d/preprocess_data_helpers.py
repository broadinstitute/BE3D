"""
File: preprocess_data_helpers.py
Author: Calvin XiaoYang Hu, Yoochan Myung, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-08-02
Description: Translated from Notebook 3.1
"""

import re
import pandas as pd
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

# MOST-TO-LEAST-DELETERIOUS ORDER USED WHEN A CATEGORY LIST IS COLLAPSED WITHOUT A USER-GIVEN ONE.
# IT MATCHES THE PER-GUIDE CATEGORIES THE BUNDLED MORC2 EXAMPLE SCREENS WERE ORIGINALLY ASSIGNED #
DEFAULT_CATEGORY_PRIORITY = ['Nonsense', 'Splice-donor', 'Splice-acceptor', 'Splice Site', 'Missense',
                             'Intron', 'Silent', 'UTR', 'Flank', 'Intergenic', 'No Mutation']

def split_categories(
    value,
    delimiter,
):
    """
    The labels of a delimiter-joined category list, with blanks and trailing delimiters dropped
    ('Silent;Missense;' -> ['Silent', 'Missense']); an empty or missing value gives [].
    """
    if not isinstance(value, str):
        return []
    return [t.strip() for t in value.split(delimiter) if t.strip()]

def resolve_guide_category(
    df,
    guide_category_col=None,
    edit_categories_col=None,
    delimiter=';',
    priority=None,
    edits_col=None,
    label_roles=None,
    empty_label='No Mutation',
):
    """
    One mutation category per guide, from whichever category column(s) are named.

    Parameters
    ----------
    df : pd.DataFrame
        One screen, one row per guide.

    guide_category_col : str, optional
        Column already holding one category per guide; used as is. A value with more than one
        distinct label raises, since that column is then a per-edit list.

    edit_categories_col : str, optional
        Column holding a delimiter-joined category per edit (e.g. 'Silent;Missense;'). Alone, it is
        collapsed to one category with `priority`; next to guide_category_col, it is collapsed only
        to cross-check that column.

    delimiter : str, optional (default=';')
        Delimiter inside the category list and edits columns.

    priority : list of str, optional
        Most-to-least-deleterious order for collapsing; DEFAULT_CATEGORY_PRIORITY when None.

    edits_col : str, optional
        Amino-acid edits column. A guide with an empty category list becomes `empty_label` only
        when it also has no edits; without edits_col every empty list becomes `empty_label`.

    label_roles : dict, optional
        Label -> role (e.g. {'Splice Site': 'splice', 'Splice-donor': 'splice'}); the cross-check
        compares roles, so two labels of the same role agree. Unlisted labels compare as themselves.

    empty_label : str, optional (default='No Mutation')

    Returns
    -------
    categories : pd.Series
        One category per guide, aligned with df.

    report : dict
        'mode' ('guide', 'collapse' or 'guide+check'), 'n_guides', 'n_empty' (guides left without a
        category), 'n_collapsed' (multi-label guides collapsed), 'n_filled' (empty lists set to
        empty_label), 'unknown_labels' (labels missing from the priority order), 'n_fallback'
        (multi-label guides none of whose labels is ranked, which keep the first label of their list)
        and, for 'guide+check', 'mismatches' (a DataFrame of guides whose two columns disagree).
    """
    if not guide_category_col and not edit_categories_col:
        raise ValueError('name guide_category_col (one category per guide) or edit_categories_col (per-edit list)')
    priority = list(priority) if priority else DEFAULT_CATEGORY_PRIORITY
    report = {'n_guides': len(df), 'n_collapsed': 0, 'n_filled': 0, 'n_fallback': 0, 'unknown_labels': set()}

    def _from_guide_column():
        labels = df[guide_category_col].apply(lambda v: split_categories(v, delimiter))
        multi = labels[labels.apply(lambda ls: len(set(ls)) > 1)]
        if len(multi):
            raise ValueError(f'{guide_category_col} holds per-edit category lists (e.g. {df.loc[multi.index[0], guide_category_col]!r}), '
                             f'not one category per guide; name it as edit_categories_col instead')
        # A SINGLE LABEL MAY STILL CARRY A TRAILING DELIMITER ('Missense;') OR REPEAT ('Missense;Missense') #
        return labels.apply(lambda ls: ls[0] if ls else pd.NA)

    def _from_edit_column():
        labels = df[edit_categories_col].apply(lambda v: split_categories(v, delimiter))
        if edits_col:
            no_edits = df[edits_col].apply(lambda v: not isinstance(v, str) or not v.strip())
        else:
            no_edits = pd.Series(True, index=df.index)
        out = []
        for ls, empty_edits in zip(labels, no_edits):
            if not ls:
                if empty_edits:
                    report['n_filled'] += 1
                    out.append(empty_label)
                else:
                    out.append(pd.NA)
                continue
            report['unknown_labels'].update(l for l in ls if l not in priority)
            if len(set(ls)) > 1:
                report['n_collapsed'] += 1
                # NONE OF ITS LABELS IS RANKED, SO reduce_mutation_type KEEPS THE FIRST ONE #
                report['n_fallback'] += not any(l in priority for l in ls)
            out.append(reduce_mutation_type(delimiter.join(ls), delimiter, priority))
        return pd.Series(out, index=df.index, dtype=object)

    if guide_category_col:
        categories = _from_guide_column()
        report['mode'] = 'guide'
        if edit_categories_col:
            report['mode'] = 'guide+check'
            collapsed = _from_edit_column()
            role = lambda v: None if pd.isna(v) else (label_roles or {}).get(v, v)
            differ = pd.Series([role(x) != role(y) for x, y in zip(categories, collapsed)], index=df.index)
            report['mismatches'] = pd.DataFrame({
                guide_category_col: df.loc[differ, guide_category_col],
                edit_categories_col: df.loc[differ, edit_categories_col],
                'collapsed': collapsed[differ],
            })
    else:
        categories = _from_edit_column()
        report['mode'] = 'collapse'
    report['n_empty'] = int(categories.isna().sum())
    return categories, report

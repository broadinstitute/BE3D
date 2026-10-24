import numpy as np
import pandas as pd

from beclust3d.qa.hypothesis_tests import hypothesis_test
from beclust3d.qa.hypothesis_tests_helpers import comparison_status


def screen(gene_rows, rng):
    """Guides as (gene, category, count, mean score); cases are well separated from controls."""
    rows = [(g, cat, rng.normal(mean, 0.2)) for g, cat, n, mean in gene_rows for _ in range(n)]
    df = pd.DataFrame(rows, columns=['Gene', 'Mutation_type', 'sgRNA_score'])
    df['Gene'] = df['Gene'].where(df['Gene'] != '')  # an empty gene cell reads back as NaN
    return df


def run(tmp_path, dfs, names, cases, controls, screen_wide=None):
    return hypothesis_test(str(tmp_path), dfs, names, cases, controls, comp_name='x', mut_col='Mutation_type',
                           val_col='sgRNA_score', gene_col='Gene', save_type='png', screen_wide_categories=screen_wide)


def test_comparison_status():
    some, none = pd.DataFrame({'v': [1]}), pd.DataFrame()
    assert comparison_status(some, some) == 'tested'
    assert comparison_status(none, some) == 'no_cases'
    assert comparison_status(some, none) == 'no_controls'
    assert comparison_status(none, none) == 'no_cases_no_controls'


def test_untested_rows_are_empty_with_a_status_and_absent_genes_are_skipped(tmp_path):
    rng = np.random.default_rng(0)
    human = screen([('MORC2', 'Nonsense', 10, -1.5), ('MORC2', 'No Mutation', 20, 0), ('ATF7IP', 'No Mutation', 20, 0)], rng)
    mouse = screen([('Morc2a', 'Nonsense', 10, -1.5), ('Morc2a', 'No Mutation', 20, 0)], rng)
    mw1, mw2, ks1, ks2 = run(tmp_path, [human, mouse], ['H', 'M'], ['Nonsense'], ['No Mutation'])
    rows = ks1.set_index(['screenid', 'gene_name'])
    # a mouse gene in the human screen (and vice versa) gets no row at all
    assert sorted(rows.index) == [('H', 'ATF7IP'), ('H', 'MORC2'), ('M', 'Morc2a')]
    assert rows.loc[('H', 'MORC2'), 'status'] == 'tested' and rows.loc[('H', 'MORC2'), 'p_x'] < 1e-3
    assert rows.loc[('H', 'ATF7IP'), 'status'] == 'no_cases'
    assert np.isnan(rows.loc[('H', 'ATF7IP'), 'p_x']) and np.isnan(rows.loc[('H', 'ATF7IP'), 'D_x'])
    saved = pd.read_csv(tmp_path / 'hypothesis_qa' / 'KolmogorovSmirnov_hypothesis1.tsv', sep='\t')
    assert (saved['p_x'] != -999).all() and 'status' in saved


def test_hypothesis_two_pools_controls_of_the_whole_screen(tmp_path):
    rng = np.random.default_rng(1)
    df = screen([('G', 'Nonsense', 10, -1.5), ('G', 'No Mutation', 5, 0), ('H', 'No Mutation', 7, 0), ('', 'No Mutation', 3, 0)], rng)
    _, _, ks1, ks2 = run(tmp_path, [df], ['S'], ['Nonsense'], ['No Mutation'])
    assert ks1.set_index('gene_name').loc['G', 'num_of_controls'] == 5
    assert ks2.set_index('gene_name').loc['G', 'num_of_controls'] == 15  # including the guides with no gene
    assert sorted(ks2['gene_name']) == ['G', 'H']  # no row for the gene-less guides


def test_screen_wide_controls_reach_every_gene(tmp_path):
    rng = np.random.default_rng(2)
    df = screen([('G', 'Nonsense', 10, -1.5), ('INTERGENIC', 'Intergenic', 30, 0), ('', 'Intergenic', 30, 0)], rng)
    _, _, ks1, ks2 = run(tmp_path, [df], ['S'], ['Nonsense'], ['Intergenic'], screen_wide=['Intergenic'])
    for table in (ks1, ks2):
        assert table['gene_name'].tolist() == ['G']  # no rows for INTERGENIC or the empty gene
        assert table.loc[0, 'num_of_controls'] == 60 and table.loc[0, 'status'] == 'tested'


def test_no_case_category_anywhere_does_not_crash(tmp_path):
    rng = np.random.default_rng(3)
    df = screen([('G', 'No Mutation', 10, 0)], rng)
    _, _, ks1, ks2 = run(tmp_path, [df], ['S'], ['Splice Site'], ['No Mutation'])
    assert ks2.loc[0, 'status'] == 'no_cases'


def test_genes_are_listed_in_sorted_order(tmp_path):
    rng = np.random.default_rng(4)
    df = screen([(g, 'No Mutation', 3, 0) for g in ['ZNF1', 'ATF7IP', 'MORC2']], rng)
    _, _, ks1, _ = run(tmp_path, [df], ['S'], ['Nonsense'], ['No Mutation'])
    assert ks1['gene_name'].tolist() == ['ATF7IP', 'MORC2', 'ZNF1']

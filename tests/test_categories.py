import numpy as np
import pandas as pd
import pytest

from beclust3d.lfc3d.preprocess_data_helpers import (
    DEFAULT_CATEGORY_PRIORITY, reduce_mutation_type, resolve_guide_category, split_categories,
)

SPLICE_ROLES = {'Splice Site': 'splice', 'Splice-donor': 'splice', 'Splice-acceptor': 'splice'}


def test_split_categories_drops_blanks_and_trailing_delimiters():
    assert split_categories('Silent;Missense;', ';') == ['Silent', 'Missense']
    assert split_categories(' Missense ; ', ';') == ['Missense']
    assert split_categories(np.nan, ';') == []
    assert split_categories('', ';') == []


def test_reduce_mutation_type_ranks_then_falls_back_to_first_label():
    priority = ['Nonsense', 'Missense', 'Silent']
    assert reduce_mutation_type('Silent;Missense;', ';', priority) == 'Missense'
    assert reduce_mutation_type('Missense;Missense;', ';', priority) == 'Missense'
    assert reduce_mutation_type('UTR;Flank', ';', priority) == 'UTR'


def test_default_priority_ranks_intron_above_silent():
    assert DEFAULT_CATEGORY_PRIORITY.index('Intron') < DEFAULT_CATEGORY_PRIORITY.index('Silent')


def test_guide_column_is_used_as_is():
    df = pd.DataFrame({'cat': ['Missense', 'No Mutation', 'Missense;', 'Silent;Silent']})
    categories, report = resolve_guide_category(df, guide_category_col='cat')
    assert categories.tolist() == ['Missense', 'No Mutation', 'Missense', 'Silent']
    assert report['mode'] == 'guide' and report['n_empty'] == 0


def test_guide_column_with_per_edit_lists_is_rejected():
    df = pd.DataFrame({'cat': ['Missense', 'Silent;Missense']})
    with pytest.raises(ValueError, match='edit_categories_col'):
        resolve_guide_category(df, guide_category_col='cat')


def test_collapse_uses_priority_and_fills_guides_without_edits():
    df = pd.DataFrame({
        'types': ['Silent;Missense;', 'Silent;Intron', np.nan, np.nan, 'UTR;Flank'],
        'edits': ['S3S;L4F', 'K5K', np.nan, 'A1T', '(NC)'],
    })
    categories, report = resolve_guide_category(df, edit_categories_col='types', edits_col='edits',
                                                priority=['Missense', 'Intron', 'Silent', 'No Mutation'])
    assert categories.iloc[:3].tolist() == ['Missense', 'Intron', 'No Mutation']
    assert pd.isna(categories.iloc[3])  # an empty list with edits is left unresolved
    assert categories.iloc[4] == 'UTR'
    assert report == {**report, 'mode': 'collapse', 'n_collapsed': 3, 'n_filled': 1, 'n_fallback': 1,
                      'n_empty': 1, 'unknown_labels': {'UTR', 'Flank'}}


def test_collapse_defaults_to_the_built_in_order():
    df = pd.DataFrame({'types': ['Silent;Intron', 'Missense;Nonsense', 'Missense;Splice-donor']})
    categories, _ = resolve_guide_category(df, edit_categories_col='types')
    assert categories.tolist() == ['Intron', 'Nonsense', 'Splice-donor']


def test_both_columns_cross_check_by_role():
    df = pd.DataFrame({
        'cat': ['Missense', 'Splice Site', 'Silent'],
        'types': ['Silent;Missense', 'Missense;Splice-donor', 'Missense'],
    })
    categories, report = resolve_guide_category(df, 'cat', 'types', label_roles=SPLICE_ROLES)
    assert categories.tolist() == ['Missense', 'Splice Site', 'Silent']  # the per-guide column wins
    assert report['mode'] == 'guide+check'
    assert report['mismatches'].index.tolist() == [2]
    assert report['mismatches'].loc[2, 'collapsed'] == 'Missense'


def test_cross_check_without_roles_compares_raw_labels():
    df = pd.DataFrame({'cat': ['Splice Site'], 'types': ['Splice-donor']})
    _, report = resolve_guide_category(df, 'cat', 'types')
    assert len(report['mismatches']) == 1


def test_a_category_column_is_required():
    with pytest.raises(ValueError):
        resolve_guide_category(pd.DataFrame({'x': [1]}))


@pytest.mark.parametrize('cfg, priority, expected', [
    ({'guide_category_col': 'Mutation_Category'}, None, ('Mutation_Category', 'Mutation_Category', None)),
    ({'edit_categories_col': 'Mutation_Type'}, None, ('Mutation_Type', None, 'Mutation_Type')),
    ({'guide_category_col': 'C', 'edit_categories_col': 'T'}, ['Missense'], ('C', 'C', 'T')),
    # the older mut_col: "collapse" with a priority, "use as is" without one
    ({'mut_col': 'Mutation_type'}, ['Missense'], ('Mutation_type', None, 'Mutation_type')),
    ({'mut_col': 'Mutation_Category'}, None, ('Mutation_Category', 'Mutation_Category', None)),
])
def test_resolve_category_config(be3d_local, cfg, priority, expected):
    out = be3d_local.resolve_category_config(cfg, priority)
    assert (out['mut_col'], out['guide_category_col'], out['edit_categories_col']) == expected
    assert out['priority'] == priority


def test_resolve_category_config_partner_default_column(be3d_local):
    out = be3d_local.resolve_category_config({}, ['Missense'], legacy_default_col='Mutation_type')
    assert out['edit_categories_col'] == 'Mutation_type'


def test_resolve_category_config_needs_a_column(be3d_local):
    with pytest.raises(ValueError):
        be3d_local.resolve_category_config({}, None)


def test_cfg_get_prefers_the_new_key(be3d_local, capsys):
    assert be3d_local.cfg_get({'aa_edits_col': 'A', 'edits_col': 'B'}, 'aa_edits_col', 'edits_col') == 'A'
    assert be3d_local.cfg_get({'edits_col': 'B'}, 'aa_edits_col', 'edits_col', where='database.') == 'B'
    assert 'database.edits_col is deprecated' in capsys.readouterr().out
    assert be3d_local.cfg_get({}, 'aa_edits_col', 'edits_col', default=';') == ';'

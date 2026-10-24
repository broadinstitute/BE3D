import io

import pandas as pd
import pytest

from beclust3d.lfc3d.preprocess_data import load_controls, parse_be_data

# A REAL SCREEN ALWAYS HAS SOME EDITS, WHICH KEEPS ITS EDITS COLUMN TEXT EVEN WHERE A GUIDE HAS NONE #
EDITED_GUIDE = ('G', 'Missense', 'L4F', 0.5)


def screen(rows):
    """A screen table as parse_be_data gets it: written to and read back from a TSV, so empty cells are NaN."""
    df = pd.DataFrame(rows, columns=['Gene', 'Mutation_type', 'Mutation_list', 'sgRNA_score'])
    return pd.read_csv(io.StringIO(df.to_csv(sep='\t', index=False)), sep='\t')


def parse(tmp_path, df, categories, screen_wide=None):
    return parse_be_data(str(tmp_path), [df], 'G', ['S'], mut_col='Mutation_type', val_col='sgRNA_score',
                         gene_col='Gene', edits_col='Mutation_list', mut_categories=categories,
                         mut_delimiter=';', conserv_dfs=[None], gene_list=['G'],
                         screen_wide_categories=screen_wide)


def test_each_missense_edit_gets_its_own_row_with_the_guide_lfc(tmp_path):
    df = screen([
        ('G', 'Missense', 'E977K;E978K', 0.67),
        ('G', 'Missense', 'Q866Q;E867K;G868S', -0.14),  # the silent edit is dropped
        ('G', 'Missense', 'M443I', 0.02),
    ])
    out = parse(tmp_path, df, ['Missense'])['S']['Missense']
    assert out['this_edit'].tolist() == ['E977K', 'E978K', 'E867K', 'G868S', 'M443I']
    assert out.index.tolist() == [0, 0, 1, 1, 2]  # rows of one guide share its index
    assert out['LFC'].tolist() == [0.67, 0.67, -0.14, -0.14, 0.02]


def test_silent_and_nonsense_groups_keep_only_their_own_edits(tmp_path):
    df = screen([('G', 'Silent', 'S3S;K5K', 0.1), ('G', 'Nonsense', 'L4F;Q10*', -1.2)])
    out = parse(tmp_path, df, ['Silent', 'Nonsense'])['S']
    assert out['Silent']['this_edit'].tolist() == ['S3S', 'K5K']
    assert out['Nonsense']['this_edit'].tolist() == ['Q10*']


def test_guides_without_edits_keep_one_row_each_even_with_equal_lfcs(tmp_path):
    df = screen([('G', 'No Mutation', '', 0.1), ('G', 'No Mutation', '', 0.1), ('G', 'No Mutation', '', -0.3), EDITED_GUIDE])
    out = parse(tmp_path, df, ['No Mutation'])['S']['No Mutation']
    assert out.tolist() == [0.1, 0.1, -0.3]


def test_other_genes_are_left_out(tmp_path):
    df = screen([('G', 'No Mutation', '', 0.1), ('OTHER', 'No Mutation', '', 0.2), EDITED_GUIDE])
    assert parse(tmp_path, df, ['No Mutation'])['S']['No Mutation'].tolist() == [0.1]


def test_screen_wide_categories_come_from_the_whole_screen(tmp_path):
    df = screen([('G', 'No Mutation', '', 0.1), ('INTERGENIC', 'Intergenic', '', 0.2), ('', 'Intergenic', '', 0.3), EDITED_GUIDE])
    out = parse(tmp_path, df, ['No Mutation', 'Intergenic'], screen_wide=['Intergenic'])['S']
    assert out['Intergenic'].tolist() == [0.2, 0.3]
    assert (tmp_path / 'screendata' / 'G_S_Intergenic.tsv').exists()


def test_without_screen_wide_intergenic_guides_are_lost(tmp_path):
    df = screen([('G', 'No Mutation', '', 0.1), ('INTERGENIC', 'Intergenic', '', 0.2), EDITED_GUIDE])
    assert 'Intergenic' not in parse(tmp_path, df, ['No Mutation', 'Intergenic'])['S']


def test_load_controls_pools_categories_one_row_per_guide(tmp_path):
    df = screen([('G', 'No Mutation', '', 0.1), ('G', 'Silent', 'S3S;K5K', 0.2), ('G', 'Silent', 'A1A', 0.3)])
    parse(tmp_path, df, ['No Mutation', 'Silent'])
    assert sorted(load_controls(str(tmp_path), 'G', 'S')['LFC']) == [0.1]
    # the two-edit silent guide is counted once
    assert sorted(load_controls(str(tmp_path), 'G', 'S', ['No Mutation', 'Silent'])['LFC']) == [0.1, 0.2, 0.3]
    assert sorted(load_controls(str(tmp_path), '*', 'S', ['No Mutation', 'Silent'])['LFC']) == [0.1, 0.2, 0.3]


def test_load_controls_reports_missing_categories(tmp_path, capsys):
    parse(tmp_path, screen([('G', 'No Mutation', '', 0.1), EDITED_GUIDE]), ['No Mutation'])
    assert len(load_controls(str(tmp_path), 'G', 'S', ['No Mutation', 'Intergenic'])) == 1
    assert "no ['Intergenic'] guides" in capsys.readouterr().out
    with pytest.raises(FileNotFoundError):
        load_controls(str(tmp_path), 'G', 'S', ['Intergenic'])

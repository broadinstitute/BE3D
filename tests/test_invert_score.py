"""
Tests for the opt-in `invert_score` parameter of parse_be_data.

`invert_score=True` negates the numeric value column (val_col) of every screen
BEFORE per-mutation parsing/aggregation, so that enrichment / activity-reporter
screens (where loss-of-function is a POSITIVE score) end up with the same
dropout sign convention (LOF -> negative) the rest of the pipeline assumes.

These tests are fully self-contained and offline: they build a small synthetic
screen DataFrame and only exercise parse_be_data (no FASTA/PDB/structure step).
"""

import pandas as pd
import pytest

from beclust3d.lfc3d.preprocess_data import parse_be_data


GENE = "TESTGENE"
SCREEN = "screenA"


def _synthetic_screen():
    """A tiny screen with one guide per mutation category.

    The Nonsense guide is a POSITIVE-scored 'hit' (as in an enrichment /
    activity-reporter screen), where LOF shows up as a positive score.
    """
    return pd.DataFrame(
        {
            "Mutation category": [
                "Nonsense",
                "Missense",
                "Silent",
                "No Mutation",
            ],
            "Amino Acid Edits": [
                "Q100*",   # nonsense (LOF)
                "A50T",    # missense
                "G30G",    # silent
                "",        # no mutation
            ],
            "Target Gene Symbol": [GENE, GENE, GENE, GENE],
            "logFC": [2.5, -1.0, 0.3, 0.1],
        }
    )


def _lfc_series(subset):
    """Return the LFC values of a parse_be_data output as a plain list.

    parse_be_data returns a DataFrame (with an 'LFC' column) for most
    categories and a bare Series for 'No Mutation'.
    """
    if isinstance(subset, pd.Series):
        return list(subset.values)
    return list(subset["LFC"].values)


def _run(tmp_path, invert):
    df = _synthetic_screen()
    return parse_be_data(
        workdir=str(tmp_path),
        input_dfs=[df],
        input_gene=GENE,
        screen_names=[SCREEN],
        conserv_dfs=[None],
        gene_list=[GENE],
        invert_score=invert,
    )


def test_values_exactly_negated(tmp_path):
    """(a) Every per-mutation LFC value is exactly negated between runs."""
    out_false = _run(tmp_path / "noinv", invert=False)
    out_true = _run(tmp_path / "inv", invert=True)

    muts_false = out_false[SCREEN]
    muts_true = out_true[SCREEN]
    assert set(muts_false.keys()) == set(muts_true.keys())

    # There must be at least one non-empty category to make the test meaningful.
    seen_any = False
    for mut in muts_false:
        vals_false = _lfc_series(muts_false[mut])
        vals_true = _lfc_series(muts_true[mut])
        assert len(vals_false) == len(vals_true)
        for a, b in zip(vals_false, vals_true):
            if pd.isna(a) or pd.isna(b):
                assert pd.isna(a) and pd.isna(b)
                continue
            seen_any = True
            assert a == pytest.approx(-b)
    assert seen_any, "expected at least one numeric LFC value to compare"


def test_positive_hit_becomes_negative(tmp_path):
    """(b) A POSITIVE-scored nonsense hit (invert=False) lands in the
    LOF / negative direction when invert=True."""
    out_false = _run(tmp_path / "noinv", invert=False)
    out_true = _run(tmp_path / "inv", invert=True)

    non_false = _lfc_series(out_false[SCREEN]["Nonsense"])
    non_true = _lfc_series(out_true[SCREEN]["Nonsense"])

    assert len(non_false) == 1 and len(non_true) == 1
    assert non_false[0] > 0          # positive score without inversion
    assert non_true[0] < 0           # flipped to the LOF / neg direction
    assert non_true[0] == pytest.approx(-non_false[0])


def test_default_is_unchanged(tmp_path):
    """Default (invert_score omitted) matches invert_score=False byte-for-byte."""
    df = _synthetic_screen()
    out_default = parse_be_data(
        workdir=str(tmp_path / "default"),
        input_dfs=[df],
        input_gene=GENE,
        screen_names=[SCREEN],
        conserv_dfs=[None],
        gene_list=[GENE],
    )
    out_false = _run(tmp_path / "false", invert=False)

    for mut in out_default[SCREEN]:
        assert _lfc_series(out_default[SCREEN][mut]) == _lfc_series(
            out_false[SCREEN][mut]
        )

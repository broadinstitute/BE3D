"""
File: tests/test_categories.py
Tests for the configurable / tolerant mutation-category handling added to
beclust3d.lfc3d.preprocess_data.parse_be_data:

  (a) a screen whose control token is 'UTR' (not 'No Mutation') parses and
      writes its control file when control_category='UTR' is passed,
  (b) an extra non-standard category ('Frameshift') passed in mut_categories
      gets its own parsed bucket rather than being dropped,
  (c) a category present in the data but absent from mut_categories raises a
      warning (not a silent drop),
  (d) the default call (no new args) is unchanged.

All tests are self-contained and offline: they build tiny in-memory
DataFrames and call parse_be_data only, with user_fasta/user_pdb unused
(the structure step is never reached).
"""

import warnings

import pandas as pd
import pytest

from beclust3d.lfc3d.preprocess_data import parse_be_data


# Column names parse_be_data expects (its defaults).
MUT_COL = "Mutation category"
VAL_COL = "logFC"
GENE_COL = "Target Gene Symbol"
EDITS_COL = "Amino Acid Edits"

GENE = "GENEX"
SCREEN = "screenA"


def _make_df(rows):
    """rows: list of (mutation_category, amino_acid_edits, logFC)."""
    return pd.DataFrame(
        {
            GENE_COL: [GENE] * len(rows),
            MUT_COL: [r[0] for r in rows],
            EDITS_COL: [r[1] for r in rows],
            VAL_COL: [r[2] for r in rows],
        }
    )


def _run(tmp_path, df, mut_categories, control_category="No Mutation"):
    return parse_be_data(
        workdir=str(tmp_path),
        input_dfs=[df],
        input_gene=GENE,
        screen_names=[SCREEN],
        mut_categories=mut_categories,
        conserv_dfs=[None],
        gene_list=[GENE],
        control_category=control_category,
    )


def test_a_custom_control_token_utr(tmp_path):
    """Control tokened 'UTR' parses and writes its control file."""
    df = _make_df(
        [
            ("Missense", "G12D", -1.5),
            ("Missense", "R80K", -0.4),
            ("UTR", "", 0.02),
            ("UTR", "", -0.05),
            ("UTR", "", 0.10),
        ]
    )
    mut_dfs = _run(
        tmp_path,
        df,
        mut_categories=["Missense", "UTR"],
        control_category="UTR",
    )

    # The control bucket exists ...
    assert "UTR" in mut_dfs[SCREEN]
    # ... and its per-screen control file was written where downstream reads it.
    control_file = tmp_path / "screendata" / f"{GENE}_{SCREEN}_UTR.tsv"
    assert control_file.exists(), "control file for 'UTR' should be written"


def test_b_extra_category_frameshift_not_dropped(tmp_path):
    """A non-standard prime-editing category gets its own bucket."""
    df = _make_df(
        [
            ("Missense", "G12D", -1.5),
            ("Frameshift", "G12X", -2.0),
            ("Frameshift", "L44X", -1.8),
            ("No Mutation", "", 0.01),
        ]
    )
    mut_dfs = _run(
        tmp_path,
        df,
        mut_categories=["Missense", "Frameshift", "No Mutation"],
    )

    assert "Frameshift" in mut_dfs[SCREEN], "Frameshift should be retained, not dropped"
    frameshift_file = tmp_path / "screendata" / f"{GENE}_{SCREEN}_Frameshift.tsv"
    assert frameshift_file.exists()


def test_c_unknown_category_warns(tmp_path):
    """A category present in data but absent from mut_categories warns."""
    df = _make_df(
        [
            ("Missense", "G12D", -1.5),
            ("Intron", "", 0.0),  # present in data, NOT in mut_categories
            ("No Mutation", "", 0.01),
        ]
    )
    with pytest.warns(UserWarning, match="dropped"):
        _run(
            tmp_path,
            df,
            mut_categories=["Missense", "No Mutation"],
        )


def test_d_default_call_unchanged(tmp_path):
    """Default call (no new args) behaves as before: no 'dropped' warning."""
    df = _make_df(
        [
            ("Missense", "G12D", -1.5),
            ("Missense", "R80K", -0.4),
            ("Silent", "L44L", 0.05),
            ("No Mutation", "", 0.01),
            ("No Mutation", "", -0.02),
        ]
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mut_dfs = parse_be_data(
            workdir=str(tmp_path),
            input_dfs=[df],
            input_gene=GENE,
            screen_names=[SCREEN],
            conserv_dfs=[None],
            gene_list=[GENE],
        )

    # Expected buckets are present.
    for cat in ("Missense", "Silent", "No Mutation"):
        assert cat in mut_dfs[SCREEN]

    # No "dropped categories" warning for an all-standard screen.
    dropped_msgs = [
        str(w.message) for w in caught if "dropped" in str(w.message)
    ]
    assert dropped_msgs == [], f"unexpected dropped-category warning: {dropped_msgs}"

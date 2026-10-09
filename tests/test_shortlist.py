"""
Deterministic, offline tests for beclust3d.lfc3d.prioritize_report.validation_shortlist.

Constructs a small synthetic per-residue table spanning all three tiers and
asserts the ranking / tiering / fallback / truncation / sentinel behavior.
"""

import numpy as np
import pandas as pd

from beclust3d.lfc3d.prioritize_report import validation_shortlist


def _base_df():
    """Synthetic per-residue table with known, tier-spanning values.

    R10: significant + confident + reachable      -> Tier A
    R20: significant + LOW pLDDT + reachable       -> Tier B (low confidence)
    R30: significant + confident + UNREACHABLE     -> Tier B (unreachable)
    R40: NOT significant + confident + reachable   -> Tier C
    """
    return pd.DataFrame(
        {
            "unipos": [10, 20, 30, 40],
            "unires": ["R", "K", "D", "E"],
            "LFC3D": [-2.0, -3.0, -1.5, -0.2],
            "q_value": [0.01, 0.02, 0.005, 0.50],
            "p_value": [0.001, 0.002, 0.0005, 0.40],
            "bfactor_pLDDT": [90.0, 40.0, 95.0, 88.0],
            "reachable": [True, True, False, True],
        }
    )


def test_tier_a_confident_reachable_outranks_low_plddt():
    df = _base_df()
    out = validation_shortlist(
        df, score_col="LFC3D", q_col="q_value",
        plddt_col="bfactor_pLDDT", reachable_col="reachable",
    )
    row = {r.position: r for r in out.itertuples()}
    # (a) significant + confident + reachable residue is Tier A
    assert row[10].priority_tier == "A"
    # the significant-but-low-pLDDT residue is only Tier B
    assert row[20].priority_tier == "B"
    # and Tier A outranks Tier B
    assert row[10].priority_rank < row[20].priority_rank
    # ranks are a contiguous 1..N
    assert list(out["priority_rank"]) == list(range(1, len(out) + 1))


def test_unreachable_significant_is_not_tier_a():
    df = _base_df()
    out = validation_shortlist(
        df, score_col="LFC3D", q_col="q_value",
        plddt_col="bfactor_pLDDT", reachable_col="reachable",
    )
    row = {r.position: r for r in out.itertuples()}
    # (b) unreachable significant residue is NOT Tier A (demoted to B)
    assert row[30].priority_tier == "B"
    assert row[30].reachable == False  # noqa: E712
    # non-significant residue lands in Tier C
    assert row[40].priority_tier == "C"


def test_pvalue_fallback_when_no_q():
    df = _base_df()
    out = validation_shortlist(
        df, score_col="LFC3D", q_col=None, p_col="p_value",
        plddt_col="bfactor_pLDDT", reachable_col="reachable",
    )
    # (c) fallback to p_col works; sig_type reports 'p'
    assert (out["sig_type"] == "p").all()
    row = {r.position: r for r in out.itertuples()}
    assert row[10].priority_tier == "A"
    # p-value fallback should note the significance column used
    assert row[10].significance == 0.001


def test_missing_q_col_name_falls_back_to_p():
    df = _base_df()
    # ask for a q_col that isn't present -> should fall back to p and note it
    out = validation_shortlist(
        df, score_col="LFC3D", q_col="does_not_exist", p_col="p_value",
        reachable_col="reachable",
    )
    assert (out["sig_type"] == "p").all()
    assert any("falling back" in n for n in out.attrs["notes"])


def test_top_n_truncation_and_tsv_columns(tmp_path):
    df = _base_df()
    tsv = tmp_path / "shortlist.tsv"
    out = validation_shortlist(
        df, score_col="LFC3D", q_col="q_value",
        plddt_col="bfactor_pLDDT", reachable_col="reachable",
        top_n=2, out_tsv=str(tsv),
    )
    # (d) top_n truncation
    assert len(out) == 2
    assert list(out["priority_rank"]) == [1, 2]
    # TSV round-trips with the expected columns
    reread = pd.read_csv(tsv, sep="\t")
    assert len(reread) == 2
    expected_cols = [
        "priority_rank", "position", "residue", "effect", "significance",
        "sig_type", "pLDDT", "reachable", "priority_tier",
    ]
    assert list(reread.columns) == expected_cols


def test_nan_and_sentinel_handling():
    df = _base_df()
    # inject a '-' sentinel q (object dtype) and a NaN pLDDT
    df["q_value"] = df["q_value"].astype(object)
    df.loc[df["unipos"] == 10, "q_value"] = "-"
    df.loc[df["unipos"] == 20, "bfactor_pLDDT"] = np.nan
    out = validation_shortlist(
        df, score_col="LFC3D", q_col="q_value",
        plddt_col="bfactor_pLDDT", reachable_col="reachable",
    )
    row = {r.position: r for r in out.itertuples()}
    # (e) '-' significance -> not significant -> Tier C, and NaN preserved
    assert row[10].priority_tier == "C"
    assert np.isnan(row[10].significance)
    # NaN pLDDT -> not confident -> Tier B (still significant)
    assert row[20].priority_tier == "B"
    assert np.isnan(row[20].pLDDT)


def test_no_reachable_col_treats_all_reachable_and_notes():
    df = _base_df().drop(columns=["reachable"])
    out = validation_shortlist(
        df, score_col="LFC3D", q_col="q_value", plddt_col="bfactor_pLDDT",
    )
    assert out["reachable"].all()
    assert any("reachable" in n.lower() for n in out.attrs["notes"])
    # R30 (significant + confident) is now Tier A since reachability is assumed
    row = {r.position: r for r in out.itertuples()}
    assert row[30].priority_tier == "A"

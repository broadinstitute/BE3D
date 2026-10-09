"""
Self-contained, offline tests for beclust3d.reachability_report.

Builds a small synthetic per-residue missense edits table (residues 1..20,
missense edits present only at {3, 4, 7, 12}) and checks the reachability
accounting, target-residue reporting, and the written TSV. No network or
real structures are required.
"""

import pandas as pd

from beclust3d.lfc3d.reachability import reachability_report


AA = "ACDEFGHIKLMNPQRSTVWY"  # 20 residues, one per position 1..20


def _make_table(edited_positions):
    """Per-residue edits table in the BE3D `*_protein_edits.tsv` shape."""
    rows = []
    edited = set(edited_positions)
    for pos in range(1, 21):
        unires = AA[pos - 1]
        if pos in edited:
            # two distinct missense edits at this residue
            edits = f"{unires}{pos}A;{unires}{pos}G"
        else:
            edits = "-"  # BE3D sentinel for "no edit maps here"
        rows.append({"unipos": pos, "unires": unires,
                     "all_Missense_edits": edits})
    return pd.DataFrame(rows)


def test_reachability_report(tmp_path):
    df = _make_table({3, 4, 7, 12})
    out_tsv = tmp_path / "reachability.tsv"

    summary = reachability_report(
        df,
        total_residues=20,
        target_residues=[4, 5, 12, 18],
        out_tsv=out_tsv,
    )

    # Reachable set is exactly the edited positions
    assert set(summary["reachable_positions"]) == {3, 4, 7, 12}
    assert summary["n_reachable"] == 4
    assert summary["n_unreachable"] == 16
    assert summary["total_residues"] == 20
    assert summary["coverage"] == 4 / 20

    # Target-residue accounting
    targets = summary["targets"]
    assert targets["n_targets"] == 4
    assert targets["reachable"] == [4, 12]
    assert targets["unreachable"] == [5, 18]
    assert targets["n_reachable"] == 2
    assert targets["n_unreachable"] == 2

    # TSV written with the expected columns and per-residue reachability
    assert out_tsv.exists()
    df_out = pd.read_csv(out_tsv, sep="\t")
    assert list(df_out.columns) == [
        "unipos", "unires", "n_missense_guides", "reachable",
    ]
    assert len(df_out) == 20
    reached = set(df_out.loc[df_out["reachable"], "unipos"])
    assert reached == {3, 4, 7, 12}
    # each edited residue carries its two synthetic guides
    assert df_out.loc[df_out["unipos"] == 4, "n_missense_guides"].iloc[0] == 2
    assert df_out.loc[df_out["unipos"] == 5, "n_missense_guides"].iloc[0] == 0


def test_reachability_defaults_total_and_no_targets(tmp_path):
    df = _make_table({3, 4, 7, 12})
    summary = reachability_report(df)  # no total, no targets, no TSV
    assert summary["total_residues"] == 20  # defaults to len(df)
    assert summary["n_reachable"] == 4
    assert summary["coverage"] == 4 / 20
    assert "targets" not in summary
    assert summary["out_tsv"] is None

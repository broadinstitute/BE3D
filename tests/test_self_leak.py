"""
Tests for the LFC3D "self-leak" fix in beclust3d/lfc3d/calculate_lfc3d.py.

Background
----------
LFC3D for a residue is an aggregation (default mean) of the residue's OWN 1-D LFC together with
its structural neighbours' LFCs. Historically the code never required that any neighbour actually
contributed, so a residue with ZERO structural neighbours still received an LFC3D value equal to
its own 1-D LFC exactly (Delta == 0.0). Such self-leak residues can be spuriously reported as 3-D
hotspots.

These tests pin down:
  (a) a residue with 0 contributing neighbours has ``_n_neighbors == 0`` and, by default,
      ``LFC3D == its own 1-D LFC`` (documents the pre-existing behaviour);
  (b) with ``min_neighbors=1`` that residue's LFC3D becomes the missing sentinel ``'-'`` while a
      well-connected residue is unchanged;
  (c) the default call (no new arg) reproduces the pre-change LFC3D values exactly for connected
      residues.
Plus a focused unit test of the pure neighbour-counting helper.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

# Import the module by file path (the lfc3d/ dir has no __init__.py). #
_MOD_PATH = Path(__file__).resolve().parents[1] / "beclust3d" / "lfc3d" / "calculate_lfc3d.py"
_spec = importlib.util.spec_from_file_location("calculate_lfc3d", _MOD_PATH)
clfc3d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(clfc3d)

calculate_lfc3d = clfc3d.calculate_lfc3d
_count_contributing_neighbors = clfc3d._count_contributing_neighbors

LFC_COL = "mean_Missense_LFC"  # f'{function_type_lfc}_{muttype}_LFC'
SCREEN = "s1"
NRAND = 2


def _make_inputs():
    """
    4 residues on chain A (positions 1..4):
      pos1,2,3 are mutually connected (each has 2 neighbours with data);
      pos4 has NO neighbour -> self-leak case.
    LFC values are distinct so LFC3D can be checked exactly.
    """
    df_struc = pd.DataFrame({
        "unipos": [1, 2, 3, 4],
        "unires": ["A", "R", "N", "D"],
        "chain": ["A", "A", "A", "A"],
        # semicolon-separated 1-based neighbour positions; pos4 has none ('-') #
        "Naa_pos":   ["2;3", "1;3", "1;2", "-"],
        "Naa_chain": ["A;A", "A;A", "A;A", "-"],
    })

    lfc_vals = [1.0, 2.0, 3.0, 5.0]  # pos4's own LFC == 5.0 #
    df_edits = pd.DataFrame({
        "conservation": ["conserved"] * 4,
        LFC_COL: lfc_vals,
        f"{LFC_COL}_Z": [0.0, 0.0, 0.0, 0.0],
    })

    rand = {}
    for r in range(NRAND):
        rand[f"{LFC_COL}r{r + 1}"] = [0.1 * (r + 1)] * 4
    df_rand = pd.DataFrame(rand)

    return df_struc, [df_edits], [df_rand]


def _run(tmp_path, **kwargs):
    df_struc, df_edits_list, df_rand_list = _make_inputs()
    return calculate_lfc3d(
        df_struc=df_struc,
        df_edits_list=df_edits_list,
        df_rand_list=df_rand_list,
        workdir=str(tmp_path),
        input_gene="TEST",
        screen_names=[SCREEN],
        nRandom=NRAND,
        **kwargs,
    )


def test_count_contributing_neighbors_helper():
    """Pure helper: self (sources[0]) is never counted; '-' neighbours don't count."""
    taa = {0: 1.0, 1: 2.0, 2: "-", 3: 5.0}
    # self only -> 0 neighbours #
    assert _count_contributing_neighbors([("local", 3)], taa, {}) == 0
    # self + two real neighbours #
    assert _count_contributing_neighbors(
        [("local", 0), ("local", 1), ("local", 3)], taa, {}) == 2
    # self + one real + one missing neighbour -> 1 #
    assert _count_contributing_neighbors(
        [("local", 0), ("local", 1), ("local", 2)], taa, {}) == 1


def test_a_self_leak_default_behavior(tmp_path):
    """(a) 0-neighbour residue: n_neighbors==0 and LFC3D == its own 1-D LFC by default."""
    df = _run(tmp_path)  # default min_neighbors=0 #
    n_col = f"{SCREEN}_LFC3D_n_neighbors"
    assert n_col in df.columns, "additive _n_neighbors column must exist"

    # pos4 is row index 3 #
    assert float(df[n_col].iloc[3]) == 0.0
    # LFC3D collapses to the residue's own 1-D LFC (the self-leak) #
    assert float(df[f"{SCREEN}_LFC3D"].iloc[3]) == float(df[f"{SCREEN}_LFC"].iloc[3]) == 5.0
    # connected residues have 2 contributing neighbours #
    assert [float(x) for x in df[n_col].iloc[0:3]] == [2.0, 2.0, 2.0]


def test_b_min_neighbors_gates_self_leak(tmp_path):
    """(b) min_neighbors=1 -> self-leak residue LFC3D becomes '-', connected residue unchanged."""
    df = _run(tmp_path, min_neighbors=1)
    # excluded residue emitted as missing sentinel #
    assert df[f"{SCREEN}_LFC3D"].iloc[3] == "-"
    # its randomised LFC3D is excluded too (null stays consistent) #
    assert df[f"{SCREEN}_LFC3Dr1"].iloc[3] == "-"
    # well-connected residue (pos1) is unchanged: mean(1,2,3) == 2.0 #
    assert float(df[f"{SCREEN}_LFC3D"].iloc[0]) == 2.0
    # n_neighbors column still reports the true count #
    assert float(df[f"{SCREEN}_LFC3D_n_neighbors"].iloc[3]) == 0.0


def test_c_default_reproduces_connected_values(tmp_path):
    """(c) default call reproduces exact LFC3D for connected residues (mean of self+neighbours)."""
    df = _run(tmp_path)
    lfc3d = [float(x) for x in df[f"{SCREEN}_LFC3D"].iloc[0:3]]
    # pos1: mean(1,2,3)=2.0 ; pos2: mean(2,1,3)=2.0 ; pos3: mean(3,1,2)=2.0 #
    assert lfc3d == [2.0, 2.0, 2.0]
    # and the self-leak residue still present by default #
    assert float(df[f"{SCREEN}_LFC3D"].iloc[3]) == 5.0

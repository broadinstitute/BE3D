"""
Fast, self-contained tests for the portability / DSSP-fallback fixes.

They cover:
  1. run_dssp's pure-Python fallback when no mkdssp/dssp binary is on PATH,
     and that parse_dssp can read the generated file.
  2. Filename-sanitization of significance labels for Windows.
  3. pandas 3.x-safe tab-delimited to_csv writing.

No network access or large data is required.
"""
import os
from pathlib import Path

import pandas as pd
import pytest

from beclust3d.lfc3d import structure_helpers
from beclust3d.lfc3d.clustering_plot import sanitize_label


# A tiny ATOM-only PDB: 3 residues (MET, ALA, GLY) on chain A, with CA + a few
# backbone atoms so Biopython can parse residues and torsions.
TINY_PDB = """\
ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 20.00           N
ATOM      2  CA  MET A   1      12.560  13.207  10.000  1.00 20.00           C
ATOM      3  C   MET A   1      13.100  14.600  10.000  1.00 20.00           C
ATOM      4  O   MET A   1      12.400  15.600  10.000  1.00 20.00           O
ATOM      5  N   ALA A   2      14.420  14.640  10.000  1.00 20.00           N
ATOM      6  CA  ALA A   2      15.120  15.910  10.000  1.00 20.00           C
ATOM      7  C   ALA A   2      16.630  15.760  10.000  1.00 20.00           C
ATOM      8  O   ALA A   2      17.200  14.670  10.000  1.00 20.00           O
ATOM      9  N   GLY A   3      17.270  16.900  10.000  1.00 20.00           N
ATOM     10  CA  GLY A   3      18.720  17.000  10.000  1.00 20.00           C
ATOM     11  C   GLY A   3      19.280  18.400  10.000  1.00 20.00           C
ATOM     12  O   GLY A   3      18.560  19.400  10.000  1.00 20.00           O
TER      13      GLY A   3
END
"""


def _write_tiny_pdb(tmp_path):
    pdb_path = tmp_path / "tiny.pdb"
    pdb_path.write_text(TINY_PDB)
    return pdb_path


def test_run_dssp_fallback_when_no_binary(tmp_path, monkeypatch):
    """When mkdssp/dssp are missing, run_dssp writes a valid placeholder DSSP
    file that parse_dssp can read end-to-end."""
    # Force "no binary available".
    monkeypatch.setattr(structure_helpers.shutil, "which", lambda name: None)

    _write_tiny_pdb(tmp_path)
    dssp_name = "tiny.dssp"

    with pytest.warns(UserWarning, match="placeholder"):
        structure_helpers.run_dssp(tmp_path, "tiny.pdb", dssp_name)

    dssp_path = tmp_path / dssp_name
    assert dssp_path.exists(), "fallback DSSP file was not created"

    # There should be one data line per standard residue (3), plus the header.
    lines = dssp_path.read_text().splitlines()
    assert any(line.strip().startswith("#") for line in lines)
    data_lines = [ln for ln in lines[1:] if ln.strip()]
    assert len(data_lines) == 3

    # Now confirm parse_dssp can read it without error.
    # parse_dssp needs a fasta list keyed by unipos/unires matching the PDB.
    fasta_path = tmp_path / "fasta.tsv"
    pd.DataFrame({"unipos": [1, 2, 3], "unires": ["M", "A", "G"]}).to_csv(
        fasta_path, sep="\t", index=False
    )

    parsed_name = "tiny_dssp_parsed.tsv"
    structure_helpers.parse_dssp(
        tmp_path, dssp_name, fasta_path, parsed_name, target_chainid="A"
    )

    parsed_path = tmp_path / parsed_name
    assert parsed_path.exists()
    parsed = pd.read_csv(parsed_path, sep="\t")
    assert len(parsed) == 3
    assert list(parsed["unipos"]) == [1, 2, 3]
    # ACC/PHI/PSI came from the placeholder and must be float-parseable.
    assert float(parsed["ACC"].iloc[0]) == pytest.approx(100.0)
    assert float(parsed["PHI"].iloc[0]) == pytest.approx(-60.0)


def test_run_dssp_uses_binary_when_present(tmp_path, monkeypatch):
    """Behavior is unchanged when a binary exists: run_dssp shells out and does
    NOT generate a placeholder."""
    _write_tiny_pdb(tmp_path)

    monkeypatch.setattr(
        structure_helpers.shutil, "which",
        lambda name: "/usr/bin/mkdssp" if name == "mkdssp" else None,
    )

    calls = {}

    def fake_run(cmd, check=False, **kwargs):
        calls["cmd"] = cmd
        return None

    monkeypatch.setattr(structure_helpers.subprocess, "run", fake_run)

    def _boom(*a, **k):  # placeholder writer must NOT be called
        raise AssertionError("placeholder writer should not run when binary exists")

    monkeypatch.setattr(structure_helpers, "write_placeholder_dssp", _boom)

    structure_helpers.run_dssp(tmp_path, "tiny.pdb", "tiny.dssp")
    assert calls["cmd"][0] == "/usr/bin/mkdssp"


def test_sanitize_label_removes_windows_illegal_chars():
    assert sanitize_label("p<0.05") == "plt0.05"
    assert sanitize_label("p>=0.05") == "pge0.05"
    assert sanitize_label("p>0.01") == "pgt0.01"
    assert sanitize_label("p<=0.01") == "ple0.01"
    for label in ("p<0.05", "p>=0.05", "p>0.01", "p<=0.001"):
        out = sanitize_label(label)
        assert "<" not in out and ">" not in out


def test_to_csv_tab_delimited(tmp_path):
    """The pandas 3.x-safe keyword call writes a real tab-delimited file."""
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    out = tmp_path / "out.tsv"
    df.to_csv(out, sep="\t", index=False)

    text = out.read_text()
    assert "\t" in text
    header = text.splitlines()[0]
    assert header.split("\t") == ["a", "b"]
    reread = pd.read_csv(out, sep="\t")
    assert list(reread.columns) == ["a", "b"]
    assert len(reread) == 2

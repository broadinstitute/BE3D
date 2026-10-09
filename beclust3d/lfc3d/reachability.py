"""
File: reachability.py
Author: Amit Shenoy
Date: 2026-08-29
Description:
    Base-editing reachability report.

    A base editor (CBE / ABE) can only install a subset of amino-acid
    substitutions, so many validated functional residues are simply NOT
    REACHABLE by a tiling screen (e.g. EGFR T790M / C797S / L858R,
    BRAF V600E, PIK3CA H1047R, AKT1 E17K, DNMT3A R882H all require
    transversions or indels that no base editor makes). A residue with no
    editing guide can never be scored as a hotspot, so BE3D "misses" are
    frequently the assay's coverage limit rather than an algorithmic failure.

    This module is an additive, off-critical-path reporting utility: it does
    NOT change LFC3D scoring or clustering. Given the per-residue missense
    edits table BE3D already builds (the `*_protein_edits.tsv` produced by
    `prioritize_by_sequence`), it reports which residues actually received a
    missense edit ("reachable") and, optionally, whether a user-supplied set
    of known functional / hotspot residues was reachable at all.
"""

import os
from pathlib import Path

import pandas as pd


def _count_edits(value):
    """Count the number of distinct edits encoded in an `all_*_edits` cell.

    BE3D stores per-residue edits as a ';'-joined string (e.g. 'G12C;G12D'),
    or the sentinel '-' / empty / NaN when no edit maps to that residue.
    Returns an int (0 when no edit is present).
    """
    if value is None:
        return 0
    # NaN check (floats) without importing numpy
    if isinstance(value, float) and pd.isna(value):
        return 0
    text = str(value).strip()
    if text == '' or text == '-':
        return 0
    return len([tok for tok in text.split(';') if tok.strip() not in ('', '-')])


def reachability_report(
    df_protein,
    total_residues=None,
    target_residues=None,
    out_tsv=None,
    mutation_type='Missense',
    pos_col='unipos',
    res_col='unires',
    edits_col=None,
):
    """Compute a base-editing reachability report for a screened protein.

    A residue is "reachable" if the screen installed at least one missense
    edit at it (i.e. the per-residue edits cell is not the '-' sentinel).
    Unreachable residues can never become hotspots, so this report lets a
    user distinguish an editor-coverage gap from a true biological negative.

    Parameters
    ----------
    df_protein : pd.DataFrame
        Per-residue missense edits table as produced by
        `prioritize_by_sequence` (the `*_protein_edits.tsv`). Must contain a
        residue-position column (`pos_col`) and an edits column
        (`edits_col`, default `all_{mutation_type}_edits`). A residue name
        column (`res_col`) is used if present.
    total_residues : int, optional
        Total number of residues in the protein (its sequence/structure span).
        Coverage is computed against this. Defaults to `len(df_protein)`.
    target_residues : list of int, optional
        Positions of known functional / hotspot residues to check explicitly
        (e.g. clinical driver or resistance residues). The report states how
        many were reachable vs unreachable and lists the unreachable ones.
    out_tsv : str or Path, optional
        Path to write the per-residue TSV. If None, no file is written.
    mutation_type : str, optional (default 'Missense')
        Mutation type whose edits define reachability. Used to build the
        default `edits_col` name (`all_{mutation_type}_edits`).
    pos_col : str, optional (default 'unipos')
        Column holding the residue position.
    res_col : str, optional (default 'unires')
        Column holding the residue identity (one-letter). Optional.
    edits_col : str, optional
        Column holding the ';'-joined edits. Defaults to
        `all_{mutation_type}_edits`.

    Returns
    -------
    dict
        Summary dictionary with keys:
        `total_residues`, `n_reachable`, `n_unreachable`, `coverage`
        (fraction in [0, 1]), `reachable_positions` (sorted list),
        `unreachable_positions` (sorted list), `mutation_type`, `out_tsv`,
        and, when `target_residues` is given, a nested `targets` dict with
        `n_targets`, `n_reachable`, `n_unreachable`, `reachable`,
        `unreachable` (sorted lists of positions).
    """
    if edits_col is None:
        edits_col = f'all_{mutation_type}_edits'
    if pos_col not in df_protein.columns:
        raise KeyError(f"Position column '{pos_col}' not found in df_protein.")
    if edits_col not in df_protein.columns:
        raise KeyError(
            f"Edits column '{edits_col}' not found in df_protein. "
            f"Pass edits_col=... or set mutation_type to match your table."
        )

    rows = []
    reachable_positions = []
    for _, row in df_protein.iterrows():
        unipos = int(row[pos_col])
        unires = row[res_col] if res_col in df_protein.columns else '-'
        n_guides = _count_edits(row[edits_col])
        reachable = n_guides > 0
        if reachable:
            reachable_positions.append(unipos)
        rows.append({
            'unipos': unipos,
            'unires': unires,
            'n_missense_guides': n_guides,
            'reachable': reachable,
        })

    df_report = pd.DataFrame(rows, columns=['unipos', 'unires',
                                            'n_missense_guides', 'reachable'])

    reachable_positions = sorted(set(reachable_positions))
    reachable_set = set(reachable_positions)

    if total_residues is None:
        total_residues = len(df_protein)
    total_residues = int(total_residues)
    n_reachable = len(reachable_set)
    n_unreachable = max(total_residues - n_reachable, 0)

    # Unreachable = every position in [1..total] with no missense edit, plus
    # any positions present in the table but unreachable (kept consistent via
    # the report rows). Report unreachable positions from the full span.
    all_positions = set(range(1, total_residues + 1))
    # include out-of-span table positions in the reachable set accounting only
    all_positions |= set(int(p) for p in df_report['unipos'].tolist())
    unreachable_positions = sorted(all_positions - reachable_set)

    coverage = (n_reachable / total_residues) if total_residues > 0 else 0.0

    summary = {
        'mutation_type': mutation_type,
        'total_residues': total_residues,
        'n_reachable': n_reachable,
        'n_unreachable': n_unreachable,
        'coverage': coverage,
        'reachable_positions': reachable_positions,
        'unreachable_positions': unreachable_positions,
        'out_tsv': None,
    }

    if target_residues is not None:
        targets = [int(t) for t in target_residues]
        t_reachable = sorted([t for t in targets if t in reachable_set])
        t_unreachable = sorted([t for t in targets if t not in reachable_set])
        summary['targets'] = {
            'n_targets': len(targets),
            'n_reachable': len(t_reachable),
            'n_unreachable': len(t_unreachable),
            'reachable': t_reachable,
            'unreachable': t_unreachable,
        }

    if out_tsv is not None:
        out_path = Path(out_tsv)
        if out_path.parent and not os.path.exists(out_path.parent):
            os.makedirs(out_path.parent, exist_ok=True)
        df_report.to_csv(out_path, sep='\t', index=False)
        summary['out_tsv'] = str(out_path)

    return summary

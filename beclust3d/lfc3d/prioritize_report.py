"""
File: prioritize_report.py
Author: Amit Shenoy
Date: 2026-08-30
Description:
    Additive utility that compresses a BE3D per-residue hotspot table into a
    single ranked "what-to-validate-next" shortlist for the bench.

    BE3D over-calls hotspots (base rate of 7-50% of residues), so a scientist is
    left with a long, unranked list and no way to tell which residues are worth
    scarce wet-lab validation. This module fuses three signals BE3D already (or
    now) produces into one interpretable ranking so a lab can validate the right
    ~10 residues instead of 100:

      * statistical significance -- FDR q-value (from the FDR PR) if present,
        else falls back to the p-value;
      * structural confidence -- AlphaFold pLDDT (``bfactor_pLDDT`` in the
        structure table);
      * base-editing reachability -- was the residue editable at all (from the
        reachability PR).

    This is a read-only reporting layer. It does NOT change LFC3D, clustering,
    or any existing behavior.
"""

import numpy as np
import pandas as pd

__all__ = ["validation_shortlist"]

# Sentinel strings that BE3D uses in its tables to mark "no value".
_SENTINELS = {"-", "", "nan", "NaN", "NA", "N/A", "none", "None"}


def _to_numeric(series):
    """Coerce a column to float, mapping BE3D sentinel strings ('-', etc.) to NaN."""
    if series is None:
        return None
    s = series.copy()
    if s.dtype == object:
        s = s.map(lambda v: np.nan if (isinstance(v, str) and v.strip() in _SENTINELS) else v)
    return pd.to_numeric(s, errors="coerce")


def _to_bool_reachable(series):
    """Coerce a reachability column to a boolean flag.

    Truthy: True, 1, 'yes'/'y'/'true'/'reachable'/'edited', any nonzero number.
    Sentinels ('-', '', NaN) are treated as NOT reachable (conservative).
    """
    truthy = {"1", "true", "t", "yes", "y", "reachable", "editable", "edited"}

    def conv(v):
        if isinstance(v, str):
            t = v.strip().lower()
            if t in _SENTINELS:
                return False
            return t in truthy
        if pd.isna(v):
            return False
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        try:
            return float(v) != 0.0
        except (TypeError, ValueError):
            return False

    return series.map(conv)


def validation_shortlist(
    df,
    *,
    score_col,
    q_col=None,
    p_col=None,
    plddt_col="bfactor_pLDDT",
    reachable_col=None,
    pos_col="unipos",
    res_col="unires",
    direction_col=None,
    plddt_min=70.0,
    q_max=0.1,
    top_n=None,
    out_tsv=None,
):
    """
    Build a ranked "what-to-validate-next" shortlist from a BE3D per-residue table.

    The result is one row per candidate residue, sorted best-first, with an
    interpretable tier (A/B/C) and a numeric rank so a bench scientist can pick
    the top handful of residues to validate.

    Parameters
    ----------
    df : pandas.DataFrame
        Per-residue table (e.g. a merged hits / structure frame). Must contain
        ``score_col``, ``pos_col`` and ``res_col``. Other columns are optional.
    score_col : str
        Column holding the per-residue effect score (e.g. LFC or LFC3D). The
        ranking uses its magnitude ``|score|`` as the effect size.
    q_col : str or None, optional
        Column holding the FDR-adjusted q-value. Comes from the FDR PR (#20). If
        None (or absent from ``df``), the function falls back to ``p_col``.
    p_col : str or None, optional
        Column holding the raw p-value, used when ``q_col`` is unavailable.
    plddt_col : str, optional (default='bfactor_pLDDT')
        Column holding the AlphaFold pLDDT confidence. If absent, every residue
        is treated as low-confidence (never satisfies the confidence criterion).
    reachable_col : str or None, optional
        Column flagging whether the residue was base-editable / reachable. Comes
        from the reachability PR (#19). If None, all residues are treated as
        reachable and a note is emitted (see ``attrs['notes']``).
    pos_col, res_col : str, optional
        Position and residue-identity columns (defaults 'unipos', 'unires').
    direction_col : str or None, optional
        Optional column describing edit direction (e.g. 'pos'/'neg'); carried
        through to the output if present.
    plddt_min : float, optional (default=70.0)
        pLDDT threshold at/above which a residue counts as "confident".
    q_max : float, optional (default=0.1)
        Significance threshold on the q-value. When falling back to p-values,
        the fixed p<0.05 cutoff is used instead.
    top_n : int or None, optional
        If given, truncate the shortlist to the top ``top_n`` rows after ranking.
    out_tsv : str or None, optional
        If given, write the shortlist to this path as a tab-separated file.

    Returns
    -------
    pandas.DataFrame
        Columns: position, residue, effect (|score|), significance, sig_type
        ('q' or 'p'), pLDDT, reachable, priority_tier ('A'/'B'/'C') and
        priority_rank (1-based, best first). A ``direction`` column is added
        when ``direction_col`` is supplied. ``result.attrs['notes']`` lists any
        assumptions made (e.g. missing reachability/pLDDT columns).

    Ranking logic
    -------------
    Significant = (q < q_max) if a q-column is available, else (p < 0.05).
    Confident   = pLDDT >= plddt_min.

      * Tier A : significant AND confident AND reachable -- validate these first.
      * Tier B : significant but low-confidence OR unreachable.
      * Tier C : everything else.

    Within each tier rows are sorted by significance (ascending q/p) and then by
    effect magnitude |score| (descending). ``priority_rank`` numbers the final
    order 1..N.
    """
    if score_col not in df.columns:
        raise ValueError(f"score_col '{score_col}' not found in df.columns")
    if pos_col not in df.columns:
        raise ValueError(f"pos_col '{pos_col}' not found in df.columns")
    if res_col not in df.columns:
        raise ValueError(f"res_col '{res_col}' not found in df.columns")

    notes = []
    work = df.copy()

    # --- effect magnitude -------------------------------------------------
    effect = _to_numeric(work[score_col]).abs()

    # --- significance: prefer q, fall back to p ---------------------------
    use_q = q_col is not None and q_col in work.columns
    if q_col is not None and not use_q:
        notes.append(f"q_col '{q_col}' not present; falling back to p-value.")
    if use_q:
        significance = _to_numeric(work[q_col])
        sig_type = "q"
        sig_cutoff = q_max
    else:
        if p_col is None or p_col not in work.columns:
            notes.append(
                "No q_col or p_col available; no residue can be significant "
                "(all default to Tier C)."
            )
            significance = pd.Series(np.nan, index=work.index)
        else:
            significance = _to_numeric(work[p_col])
        sig_type = "p"
        sig_cutoff = 0.05

    significant = significance < sig_cutoff  # NaN comparisons -> False

    # --- confidence (pLDDT) ----------------------------------------------
    if plddt_col in work.columns:
        plddt = _to_numeric(work[plddt_col])
    else:
        notes.append(
            f"plddt_col '{plddt_col}' not present; treating all residues as "
            "low-confidence."
        )
        plddt = pd.Series(np.nan, index=work.index)
    confident = plddt >= plddt_min  # NaN -> False

    # --- reachability -----------------------------------------------------
    if reachable_col is not None and reachable_col in work.columns:
        reachable = _to_bool_reachable(work[reachable_col])
    else:
        if reachable_col is not None:
            notes.append(
                f"reachable_col '{reachable_col}' not present; treating all "
                "residues as reachable."
            )
        else:
            notes.append(
                "No reachable_col provided; treating all residues as reachable."
            )
        reachable = pd.Series(True, index=work.index)

    # --- tier assignment --------------------------------------------------
    tier = pd.Series("C", index=work.index, dtype=object)
    tier_b = significant & (~confident | ~reachable)
    tier[tier_b] = "B"
    tier_a = significant & confident & reachable
    tier[tier_a] = "A"

    out = pd.DataFrame(
        {
            "position": work[pos_col].values,
            "residue": work[res_col].values,
            "effect": effect.values,
            "significance": significance.values,
            "sig_type": sig_type,
            "pLDDT": plddt.values,
            "reachable": reachable.values,
            "priority_tier": tier.values,
        }
    )
    if direction_col is not None and direction_col in work.columns:
        out["direction"] = work[direction_col].values

    # --- ordering ---------------------------------------------------------
    tier_order = {"A": 0, "B": 1, "C": 2}
    out["_tier_ord"] = out["priority_tier"].map(tier_order)
    # ascending significance (best/smallest first), descending effect.
    out["_sig_sort"] = out["significance"].fillna(np.inf)
    out["_eff_sort"] = -out["effect"].fillna(-np.inf)
    out = out.sort_values(
        ["_tier_ord", "_sig_sort", "_eff_sort"], kind="mergesort"
    ).reset_index(drop=True)
    out = out.drop(columns=["_tier_ord", "_sig_sort", "_eff_sort"])
    out.insert(0, "priority_rank", np.arange(1, len(out) + 1))

    if top_n is not None:
        out = out.head(int(top_n)).reset_index(drop=True)

    out.attrs["notes"] = notes

    if out_tsv is not None:
        out.to_csv(out_tsv, sep="\t", index=False)

    return out

"""
Local spatial hotspot statistics for BE3D.

This module provides an *additive*, opt-in utility that does NOT modify the
existing LFC3D / clustering pipeline. It implements the Getis-Ord Gi\\* local
hotspot statistic (Getis & Ord 1992; Ord & Getis 1995), a standardized
local z-score that is the geospatial-statistics analogue of BE3D's LFC3D.

Motivation
----------
BE3D's LFC3D is essentially an *un-standardized* local sum of per-residue LFC
scores over a hard 6 A spatial neighborhood. The Getis-Ord Gi\\* statistic is
the textbook local hotspot statistic used across GIS, single-cell, and ecology.
For each location it compares the (optionally distance-weighted) local sum of a
variable to what would be expected under the global mean, and returns a
STANDARDIZED z-score. This makes hotspot scores:

- principled (an explicit null of spatial randomness),
- comparable across proteins (z-scores, not raw sums), and
- able to accept distance-decay / Gaussian weights instead of a hard cutoff.

Using BE3D's existing neighbor table
------------------------------------
BE3D already builds a per-residue spatial adjacency in ``structure_helpers`` and
stores it in the ``Naa`` / ``Naa_pos`` columns of the coordinate DataFrame
(``Naa_pos`` is a ``;``-joined string of neighbor ``unipos`` values per residue).
A user can convert those columns into the ``neighbor_index_lists`` argument and
run :func:`getis_ord_gi_star` directly on the per-residue LFC values as a
standardized alternative to raw LFC3D. See :func:`neighbor_lists_from_naa_pos`
for a small helper that does this conversion.
"""

from __future__ import annotations

import math
from typing import Iterable, List, Optional, Sequence

import numpy as np

__all__ = ["getis_ord_gi_star", "neighbor_lists_from_naa_pos"]


def _to_float_or_nan(v) -> float:
    """Coerce a value to float, mapping BE3D's '-'/'' sentinels and bad values to NaN."""
    if v is None:
        return float("nan")
    if isinstance(v, str):
        s = v.strip()
        if s in ("", "-", "nan", "NaN", "None"):
            return float("nan")
        try:
            return float(s)
        except ValueError:
            return float("nan")
    try:
        f = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return f


def getis_ord_gi_star(
    values: Sequence,
    neighbor_index_lists: Sequence[Iterable[int]],
    weights: Optional[Sequence[Iterable[float]]] = None,
    include_self: bool = True,
) -> np.ndarray:
    """Compute the Getis-Ord Gi\\* local hotspot z-score for each residue.

    The Gi\\* statistic (Ord & Getis 1995) for focal residue ``i`` is::

                 sum_j w_ij x_j  -  Xbar * sum_j w_ij
        Gi* = ----------------------------------------------------------
                 S * sqrt( ( n * sum_j w_ij^2 - (sum_j w_ij)^2 ) / (n - 1) )

    where ``x_j`` are the per-residue values, ``Xbar`` and ``S`` are the global
    mean and (population) standard deviation of ``x`` over all ``n`` valid
    residues, and ``j`` ranges over the neighbors of ``i`` (plus ``i`` itself
    when ``include_self=True``, which is what the star in Gi\\* denotes).

    Parameters
    ----------
    values : sequence of float
        1D array of per-residue scores (e.g. LFC by residue). NaN and BE3D's
        ``'-'`` / ``''`` sentinels are treated as missing and ignored in both
        the global statistics and the local sums.
    neighbor_index_lists : list of lists of int
        For each residue ``i``, the integer indices (into ``values``) of its
        spatial neighbors. This is exactly the adjacency BE3D already builds as
        the ``Naa`` / ``Naa_pos`` columns; use
        :func:`neighbor_lists_from_naa_pos` to convert those, or pass any
        equivalent list-of-lists.
    weights : list of lists of float, optional
        Per-neighbor weights ``w_ij`` aligned element-for-element with
        ``neighbor_index_lists`` (e.g. distance-decay or Gaussian weights).
        Defaults to binary weights of ``1.0``. When ``include_self=True`` the
        self weight is ``1.0`` (binary case) or, if weights are supplied, the
        maximum weight among the residue's neighbors (falling back to ``1.0``),
        so the focal residue is never down-weighted relative to its neighbors.
    include_self : bool, default True
        If True, include residue ``i`` itself in its own neighborhood (the
        ``*`` in Gi\\*). Standard usage is ``True``.

    Returns
    -------
    numpy.ndarray
        1D float array of Gi\\* z-scores, one per residue, same length as
        ``values``. Residues whose own value is missing, or that have no valid
        neighbors, or where the denominator is zero (e.g. degenerate/uniform
        data) are returned as ``NaN``.

    Notes
    -----
    - Positive Gi\\* indicates a hotspot (local clustering of high values);
      negative indicates a coldspot. Values near 0 indicate no local structure.
    - The global std ``S`` uses the population definition (``ddof=0``), matching
      Ord & Getis (1995).
    - This function is pure and side-effect free; it does not modify or depend
      on any BE3D pipeline state.
    """
    # --- coerce values, tracking which are valid (finite) ------------------
    x = np.array([_to_float_or_nan(v) for v in values], dtype=float)
    n_total = x.shape[0]
    if n_total == 0:
        return np.array([], dtype=float)

    valid_mask = np.isfinite(x)
    n = int(valid_mask.sum())

    out = np.full(n_total, np.nan, dtype=float)
    if n < 2:
        # Need at least 2 valid residues for the (n-1) term / a meaningful std.
        return out

    xbar = float(np.mean(x[valid_mask]))
    # population standard deviation over valid residues
    S = float(np.std(x[valid_mask], ddof=0))
    if not math.isfinite(S) or S == 0.0:
        # Uniform (or degenerate) data => no hotspots; Gi* undefined -> all NaN.
        return out

    for i in range(n_total):
        # Focal residue must itself be valid to be scored.
        if not valid_mask[i]:
            continue

        neigh = list(neighbor_index_lists[i]) if i < len(neighbor_index_lists) else []
        if weights is not None and i < len(weights):
            wlist = list(weights[i])
        else:
            wlist = None

        # Build (index, weight) pairs for valid neighbors only.
        j_idx: List[int] = []
        j_w: List[float] = []
        for k, j in enumerate(neigh):
            jj = int(j)
            if jj < 0 or jj >= n_total or not valid_mask[jj]:
                continue
            w = 1.0 if wlist is None else float(wlist[k])
            if not math.isfinite(w):
                continue
            j_idx.append(jj)
            j_w.append(w)

        if include_self:
            if i not in j_idx:
                if wlist is None:
                    self_w = 1.0
                else:
                    self_w = max(j_w) if j_w else 1.0
                j_idx.append(i)
                j_w.append(self_w)

        if not j_idx:
            continue

        w_arr = np.array(j_w, dtype=float)
        x_arr = x[np.array(j_idx, dtype=int)]

        sum_w = float(np.sum(w_arr))
        sum_w2 = float(np.sum(w_arr ** 2))
        sum_wx = float(np.sum(w_arr * x_arr))

        numerator = sum_wx - xbar * sum_w
        denom_inner = (n * sum_w2 - sum_w ** 2) / (n - 1)
        if denom_inner <= 0:
            continue
        denominator = S * math.sqrt(denom_inner)
        if denominator == 0.0 or not math.isfinite(denominator):
            continue

        out[i] = numerator / denominator

    return out


def neighbor_lists_from_naa_pos(
    unipos: Sequence,
    naa_pos: Sequence,
    sep: str = ";",
) -> List[List[int]]:
    """Convert BE3D's ``Naa_pos`` column into ``neighbor_index_lists``.

    BE3D stores neighbors as ``;``-joined ``unipos`` strings (``'-'`` when a
    residue has no neighbors). This helper maps those neighbor ``unipos`` values
    to row indices so the result can be passed straight to
    :func:`getis_ord_gi_star`.

    Parameters
    ----------
    unipos : sequence
        The residue ``unipos`` values in row order (e.g. ``df['unipos']``).
    naa_pos : sequence
        The matching ``Naa_pos`` column (``;``-joined neighbor ``unipos``
        strings, or ``'-'`` / empty for none), in the same row order.
    sep : str, default ';'
        Field separator used inside ``naa_pos`` strings.

    Returns
    -------
    list of list of int
        Per-row lists of neighbor row indices, suitable as
        ``neighbor_index_lists``.
    """
    pos_to_idx = {}
    for idx, p in enumerate(unipos):
        try:
            pos_to_idx[int(p)] = idx
        except (TypeError, ValueError):
            continue

    result: List[List[int]] = []
    for entry in naa_pos:
        neighbors: List[int] = []
        if entry is None:
            result.append(neighbors)
            continue
        s = str(entry).strip()
        if s in ("", "-", "nan", "NaN", "None"):
            result.append(neighbors)
            continue
        for tok in s.split(sep):
            tok = tok.strip()
            if not tok:
                continue
            try:
                key = int(float(tok))
            except ValueError:
                continue
            if key in pos_to_idx:
                neighbors.append(pos_to_idx[key])
        result.append(neighbors)
    return result

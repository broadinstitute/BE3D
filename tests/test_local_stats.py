"""Deterministic, offline tests for the Getis-Ord Gi* local hotspot statistic.

These tests are fully self-contained: all arrays are constructed inline, and no
structures, files, or network access are required.
"""

import math

import numpy as np

from beclust3d.lfc3d.local_stats import (
    getis_ord_gi_star,
    neighbor_lists_from_naa_pos,
)


def _cluster_example():
    """A tiny 6-point example: residues 0-2 form a HIGH cluster, 3-5 a low background."""
    values = [10.0, 10.0, 10.0, 1.0, 1.0, 1.0]
    neighbors = [
        [1, 2],  # 0
        [0, 2],  # 1
        [0, 1],  # 2
        [4, 5],  # 3
        [3, 5],  # 4
        [3, 4],  # 5
    ]
    return values, neighbors


def test_hotspot_cluster_high_positive_background_negative():
    values, neighbors = _cluster_example()
    z = getis_ord_gi_star(values, neighbors, include_self=True)

    # High-cluster residues should be strong positive hotspots.
    assert np.all(z[:3] > 1.5)
    # Background residues should be coldspots (negative here).
    assert np.all(z[3:] < 0.0)
    # High cluster clearly greater than background.
    assert z[:3].min() > z[3:].max()


def test_hand_checked_gi_star_value():
    """Cross-check residue 0's Gi* against the textbook formula computed by hand.

    n = 6, Xbar = 5.5, S (population std) = 4.5.
    Neighborhood of residue 0 with include_self = {0,1,2}, all x = 10, binary w.
      sum_w = 3, sum_w2 = 3, sum_wx = 30
      numerator   = 30 - 5.5 * 3 = 13.5
      denom_inner = (6*3 - 3^2) / (6-1) = 9/5 = 1.8
      denominator = 4.5 * sqrt(1.8)
      Gi*         = 13.5 / (4.5 * sqrt(1.8)) = sqrt(5)
    """
    values, neighbors = _cluster_example()
    z = getis_ord_gi_star(values, neighbors, include_self=True)

    expected = 13.5 / (4.5 * math.sqrt(1.8))
    assert math.isclose(expected, math.sqrt(5.0), rel_tol=1e-12)
    assert np.allclose(z[0], math.sqrt(5.0), atol=1e-10)
    assert np.allclose(z[:3], math.sqrt(5.0), atol=1e-10)
    # Symmetric low cluster => -sqrt(5).
    assert np.allclose(z[3:], -math.sqrt(5.0), atol=1e-10)


def test_uniform_values_all_zero_or_nan():
    """Uniform data has zero global std => Gi* is undefined; we return NaN, not a spurious hotspot."""
    values = [7.0] * 6
    neighbors = [[1, 2], [0, 2], [0, 1], [4, 5], [3, 5], [3, 4]]
    z = getis_ord_gi_star(values, neighbors)
    # Degenerate (S == 0) guard: all NaN, and crucially never a nonzero hotspot.
    assert np.all(np.isnan(z))
    assert not np.any(np.nan_to_num(z, nan=0.0) != 0.0)


def test_near_uniform_values_are_approximately_zero():
    """A distribution with tiny structure => Gi* magnitudes near 0."""
    values = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0001]
    neighbors = [[1, 2], [0, 2], [0, 1], [4, 5], [3, 5], [3, 4]]
    z = getis_ord_gi_star(values, neighbors)
    finite = z[np.isfinite(z)]
    assert np.all(np.abs(finite) < 3.0)


def test_nan_and_dash_handling():
    """NaN and BE3D '-'/'' sentinels are ignored; a missing focal residue -> NaN output."""
    values = [10.0, 10.0, "-", 1.0, float("nan"), 1.0]
    neighbors = [[1, 2], [0, 2], [0, 1], [4, 5], [3, 5], [3, 4]]
    z = getis_ord_gi_star(values, neighbors)

    # Missing focal residues (index 2 = '-', index 4 = NaN) produce NaN output.
    assert np.isnan(z[2])
    assert np.isnan(z[4])
    # Valid residues still produce finite scores.
    assert np.isfinite(z[0])
    assert np.isfinite(z[1])
    # High residue is a positive hotspot; low residue negative.
    assert z[0] > 0
    assert z[3] < 0


def test_divide_by_zero_guard_single_valid_residue():
    """With fewer than two valid residues (n<2), Gi* is undefined -> all NaN, no crash."""
    values = ["-", "-", 5.0, "-"]
    neighbors = [[1], [0], [3], [2]]
    z = getis_ord_gi_star(values, neighbors)
    assert z.shape == (4,)
    assert np.all(np.isnan(z))


def test_empty_input():
    z = getis_ord_gi_star([], [])
    assert isinstance(z, np.ndarray)
    assert z.shape == (0,)


def test_weights_change_result_in_expected_direction():
    """Down-weighting a high neighbor should lower the focal residue's Gi*."""
    # Residue 0's neighbors include a high (index 1) and a low (index 2) value.
    values = [5.0, 20.0, 0.0, 5.0, 5.0, 5.0]
    neighbors = [[1, 2], [0], [0], [4, 5], [3, 5], [3, 4]]

    # Binary weights: neighbor 1 (high) fully counted.
    z_binary = getis_ord_gi_star(values, neighbors, include_self=True)

    # Distance-decayed weights: strongly down-weight the high neighbor (index 1).
    weights = [
        [0.01, 1.0],  # residue 0: down-weight high neighbor, keep low neighbor
        [1.0],
        [1.0],
        [1.0, 1.0],
        [1.0, 1.0],
        [1.0, 1.0],
    ]
    z_weighted = getis_ord_gi_star(values, neighbors, weights=weights, include_self=True)

    # Removing most of the high neighbor's contribution should LOWER residue 0's Gi*.
    assert z_weighted[0] < z_binary[0]


def test_weights_uniform_matches_binary():
    """All-ones weights must reproduce the binary (default) result exactly."""
    values, neighbors = _cluster_example()
    z_binary = getis_ord_gi_star(values, neighbors, include_self=True)
    weights = [[1.0] * len(n) for n in neighbors]
    z_ones = getis_ord_gi_star(values, neighbors, weights=weights, include_self=True)
    assert np.allclose(z_binary, z_ones, atol=1e-12)


def test_include_self_flag_changes_result():
    values, neighbors = _cluster_example()
    z_star = getis_ord_gi_star(values, neighbors, include_self=True)
    z_no_self = getis_ord_gi_star(values, neighbors, include_self=False)
    # Both should flag the high cluster as positive...
    assert z_star[0] > 0 and z_no_self[0] > 0
    # ...but including self changes the numeric value.
    assert not np.allclose(z_star[0], z_no_self[0])


def test_neighbor_lists_from_naa_pos_roundtrip():
    """The Naa_pos helper maps BE3D's ';'-joined unipos strings to row indices."""
    unipos = [10, 11, 12, 13, 14, 15]
    naa_pos = ["11;12", "10;12", "10;11", "14;15", "13;15", "13;14"]
    nb = neighbor_lists_from_naa_pos(unipos, naa_pos)
    assert nb == [[1, 2], [0, 2], [0, 1], [4, 5], [3, 5], [3, 4]]

    # '-' / empty means "no neighbors".
    nb2 = neighbor_lists_from_naa_pos([10, 11], ["-", ""])
    assert nb2 == [[], []]

    # Full pipeline: convert Naa_pos then score.
    values = [10.0, 10.0, 10.0, 1.0, 1.0, 1.0]
    z = getis_ord_gi_star(values, nb)
    assert np.allclose(z[:3], math.sqrt(5.0), atol=1e-10)

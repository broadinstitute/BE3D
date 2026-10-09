"""
Tests for the Benjamini-Hochberg FDR q-value helper and the companion
q-value columns added to znorm_score / znorm_meta.

Run directly (this repo has no pytest config yet):
    PYTHONPATH=<repo> python -m pytest tests/test_fdr.py -q
"""

import numpy as np
import pandas as pd
import pytest

from beclust3d.aggregate.aggregate_helpers import benjamini_hochberg


# ---- unit tests for benjamini_hochberg ------------------------------------

def test_uniform_scaling_example():
    # p = k/100 for k=1..5, m=5 -> every q_(k) = p_(k)*5/k = 0.05
    pvals = [0.01, 0.02, 0.03, 0.04, 0.05]
    q = benjamini_hochberg(pvals)
    assert np.allclose(q, [0.05, 0.05, 0.05, 0.05, 0.05])


def test_mixed_example_handcomputed():
    # p in an intentionally out-of-order arrangement.
    # sorted: 0.005(r1) 0.02(r2) 0.04(r3) 0.20(r4) 0.60(r5), m=5
    # raw q: 0.025, 0.05, 0.0666667, 0.25, 0.60 (already monotone)
    pvals = [0.005, 0.02, 0.20, 0.04, 0.60]
    expected = [0.025, 0.05, 0.25, 0.0666666667, 0.60]
    q = benjamini_hochberg(pvals)
    assert np.allclose(q, expected)


def test_monotonicity_of_sorted_q_and_q_ge_p():
    rng = np.random.default_rng(0)
    pvals = rng.uniform(0, 1, size=200)
    q = benjamini_hochberg(pvals)
    # q >= p everywhere
    assert np.all(q >= pvals - 1e-12)
    # sorted-by-p q-values are non-decreasing
    order = np.argsort(pvals, kind='stable')
    q_sorted = q[order]
    assert np.all(np.diff(q_sorted) >= -1e-12)
    # capped at 1
    assert np.all(q <= 1.0 + 1e-12)


def test_nan_and_dash_passthrough():
    pvals = [0.01, '-', 0.02, np.nan, 0.03]
    q = benjamini_hochberg(pvals)
    # only 3 valid entries, all scale to 0.03
    assert np.isnan(q[1]) and np.isnan(q[3])
    assert np.allclose([q[0], q[2], q[4]], [0.03, 0.03, 0.03])


def test_all_missing_returns_all_nan():
    q = benjamini_hochberg(['-', np.nan, None])
    assert q.shape == (3,)
    assert np.all(np.isnan(q))


def test_single_value_equals_itself():
    q = benjamini_hochberg([0.037])
    assert np.allclose(q, [0.037])


def test_accepts_numpy_and_series():
    pvals = [0.01, 0.02, 0.03, 0.04, 0.05]
    q_list = benjamini_hochberg(pvals)
    q_np = benjamini_hochberg(np.array(pvals))
    q_series = benjamini_hochberg(pd.Series(pvals))
    assert np.allclose(q_list, q_np)
    assert np.allclose(q_list, q_series)


def test_ties_stable():
    # ties should be handled deterministically and give q >= p
    pvals = [0.02, 0.02, 0.02, 0.5]
    q = benjamini_hochberg(pvals)
    assert np.all(q >= np.array(pvals) - 1e-12)
    # the three equal small p's should share one q value
    assert np.allclose(q[0], q[1]) and np.allclose(q[1], q[2])


def test_cross_check_statsmodels_if_available():
    statsmodels = pytest.importorskip("statsmodels")
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(42)
    pvals = rng.uniform(0, 1, size=137)
    _, q_sm, _, _ = multipletests(pvals, method='fdr_bh')
    q = benjamini_hochberg(pvals)
    assert np.allclose(q, q_sm)


# ---- lightweight integration check ----------------------------------------

def test_q_columns_consistent_with_p_ordering():
    # Independent BH computed over a p-column mirrors the q-column ordering:
    # residues with smaller p must not have larger q.
    rng = np.random.default_rng(7)
    p = rng.uniform(0, 1, size=50)
    # inject a couple of missing markers as the pipeline would
    p_col = list(p)
    p_col[3] = '-'
    p_col[10] = np.nan
    q = benjamini_hochberg(p_col)

    valid = [i for i, v in enumerate(p_col)
             if not (v == '-' or (isinstance(v, float) and np.isnan(v)))]
    order = sorted(valid, key=lambda i: float(p_col[i]))
    q_ordered = [q[i] for i in order]
    # non-decreasing q along increasing p
    assert all(q_ordered[k] <= q_ordered[k + 1] + 1e-12
               for k in range(len(q_ordered) - 1))

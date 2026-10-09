"""
Tests for the opt-in inverse-variance / empirical-Bayes shrinkage
meta-aggregation option added to BE3D.

Covers:
  * inverse_variance_mean() unit behavior (primary), and
  * a light integration check that average_split_meta accepts the new
    'INVERSE_VARIANCE' option and that the default 'SUM' path is unchanged.
"""

import numpy as np
import pandas as pd
import pytest

from beclust3d.aggregate.aggregate_helpers import inverse_variance_mean
from beclust3d.aggregate.metaaggregate import average_split_meta


# ------------------------- helper unit tests (primary) ------------------------- #

def test_hand_example():
    # (2/1 + 4/4) / (1/1 + 1/4) = 3 / 1.25 = 2.4
    result = inverse_variance_mean([2.0, 4.0], [1.0, 4.0])
    assert np.isclose(result, 2.4)


def test_high_variance_screen_is_downweighted():
    values = [2.0, 4.0]
    plain_mean = np.mean(values)  # 3.0
    # second screen is very noisy -> result should be pulled toward the low-variance value 2.0
    result = inverse_variance_mean(values, [1.0, 100.0])
    assert result < plain_mean
    assert abs(result - 2.0) < abs(result - 4.0)


def test_nan_and_dash_ignored():
    # NaN and '-' entries are dropped; remaining [2.0, 4.0] w/ vars [1.0, 4.0] -> 2.4
    result_nan = inverse_variance_mean([2.0, np.nan, 4.0], [1.0, 2.0, 4.0])
    result_dash = inverse_variance_mean([2.0, '-', 4.0], [1.0, 2.0, 4.0])
    assert np.isclose(result_nan, 2.4)
    assert np.isclose(result_dash, 2.4)


def test_degenerate_zero_variance_falls_back_to_unweighted_mean():
    # all-zero (degenerate) variances -> plain mean
    result = inverse_variance_mean([2.0, 4.0, 6.0], [0.0, 0.0, 0.0])
    assert np.isclose(result, 4.0)


def test_dash_variance_falls_back_to_unweighted_mean():
    # unusable ('-') variance on a screen -> equal-weight (unweighted) mean
    result = inverse_variance_mean([2.0, 4.0], ['-', 4.0])
    assert np.isclose(result, 3.0)


def test_all_missing_returns_nan():
    assert np.isnan(inverse_variance_mean(['-', np.nan], [1.0, 2.0]))


def test_empty_returns_nan():
    assert np.isnan(inverse_variance_mean([], []))


# ------------------------- light integration test ------------------------- #

def _tiny_two_screen_frame(nRandom):
    rng = np.random.default_rng(0)
    n_res = 4
    screens = ['scA', 'scB']
    data = {
        'unipos': list(range(1, n_res + 1)),
        'unires': ['A', 'R', 'N', 'D'],
        'chain': ['A'] * n_res,
    }
    base = {
        'scA': [-1.0, 2.0, -0.5, 0.8],
        'scB': [-3.0, 1.0, 0.4, -0.2],
    }
    for s in screens:
        data[f'{s}_LFC3D'] = base[s]
        for n in range(nRandom):
            data[f'{s}_LFC3Dr{n+1}'] = rng.normal(0.0, 1.0, size=n_res).round(4)
    return pd.DataFrame(data), screens


def test_average_split_meta_inverse_variance_runs(tmp_path):
    nRandom = 5
    df, screens = _tiny_two_screen_frame(nRandom)

    out = average_split_meta(
        df.copy(), str(tmp_path), 'GENE', screens,
        score_type='LFC3D', nRandom=nRandom,
        aggr_func_name='INVERSE_VARIANCE',
    )
    assert 'INVERSE_VARIANCE_LFC3D' in out.columns
    assert 'INVERSE_VARIANCE_LFC3D_neg' in out.columns
    assert 'INVERSE_VARIANCE_LFC3D_pos' in out.columns
    assert len(out) == len(df)


def test_default_sum_path_unchanged(tmp_path):
    nRandom = 5
    df, screens = _tiny_two_screen_frame(nRandom)

    # Reference SUM output computed the old way (plain per-residue sum of screen scores).
    out_sum = average_split_meta(
        df.copy(), str(tmp_path / 'sum'), 'GENE', screens,
        score_type='LFC3D', nRandom=nRandom,
        aggr_func_name='SUM',
    )

    expected_all = []
    for i in range(len(df)):
        vals = []
        for s in screens:
            v = df.at[i, f'{s}_LFC3D']
            if v != '-' and v != 0.0:
                vals.append(float(v))
        expected_all.append(sum(vals) if vals else '-')

    got = out_sum['SUM_LFC3D'].tolist()
    for g, e in zip(got, expected_all):
        if e == '-':
            assert g == '-'
        else:
            assert np.isclose(float(g), e)


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-q']))

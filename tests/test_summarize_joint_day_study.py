"""Synthetic summary checks; no production inference or raw ancestry reads."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('summary_joint', SCRIPTS/'summarize_joint_day_study.py')
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


def toy_rows():
    rows = []
    for d in ('1d', '2d'):
        for a in s.ALPHAS:
            for h in (2, 3):
                cases = s.CASES + (['legacy_recovery', 'cartesian_mean_old_background',
                                   'frozen_mixture_old_background'] if h == 2 else [])
                for case in cases:
                    for mode in ('observed', 'matched_zero'):
                        rows.append(dict(dimension=d, alpha=a, contrasts=h, case=case,
                                         mode=mode, c95=float(h), n=2,
                                         cbg_trace=2., cbg_eigen_min=.1))
    return rows


def test_counts_duplicates_and_missing():
    rows = toy_rows()
    assert len(rows) == 672
    s.validate_counts(rows, s.ALPHAS)
    with pytest.raises(ValueError, match='duplicate'):
        s.index_rows(rows+rows[:1])
    with pytest.raises(ValueError, match='counts'):
        s.validate_counts(rows[:-1], s.ALPHAS)


def test_compare_preserves_weaker_limits_and_matches_keys():
    rows = toy_rows()[:2]
    other = [dict(r, c95=2*r['c95']) for r in reversed(rows)]
    result = s.compare_rows(rows, other)
    assert all(r['ratio'] == 2 for r in result)
    np.testing.assert_allclose([r['delta_dex'] for r in result], np.log10(2))
    with pytest.raises(ValueError, match='keys'):
        s.compare_rows(rows, other[:1])


def test_npz_standard_posterior_and_reject_axis(tmp_path):
    c = np.geomspace(1/137, 1000, 120)
    ll = -c
    density, cdf = s.normalized_posterior_cdf(c, ll)
    limit = s.posterior_upper_limit(c, ll)
    fields = dict(couplings=c, loglike=ll, density=density, cdf=cdf,
                  c95=limit, phase_loglike=np.zeros((120, 64)))
    row = dict(key='toy', c95=limit, lower_bin_mass=cdf[1], upper_bin_mass=1-cdf[-2],
               phase_count=64, likelihood_peak_at_boundary=True)
    path = tmp_path/'posterior.npz'
    np.savez(path, **fields)
    assert s.validate_posterior(path, row)['likelihood_peak_side'] == 'lower'
    fields['density'] = density[:, None]
    np.savez(path, **fields)
    with pytest.raises(ValueError, match='axis'):
        s.validate_posterior(path, row)


def test_all_four_plots_synthetic(tmp_path):
    s.make_plots(toy_rows(), tmp_path)
    assert {p.name for p in tmp_path.glob('*.png')} == {
        'c95-vs-alpha.png', 'alpha030-factorial.png',
        'background-geometry.png', 'h3-h2-ratios.png'}
    assert all(p.stat().st_size > 1000 for p in tmp_path.glob('*.png'))

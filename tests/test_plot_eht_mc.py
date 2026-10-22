"""Small geometry/branch tests; no live ensemble or FITS required."""
import importlib.util
from pathlib import Path

import numpy as np

SPEC = importlib.util.spec_from_file_location(
    'plot_eht_mc', Path(__file__).resolve().parents[1] / 'scripts/plot_eht_mc.py')
plot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plot)


def test_ticks_north_and_east():
    segments = plot.tick_segments(np.array([0., 0.]), np.array([0., 0.]),
                                  np.array([0., np.pi/2]), np.array([4., 4.]))
    np.testing.assert_allclose(segments[0], [[0, -2], [0, 2]], atol=1e-14)
    np.testing.assert_allclose(segments[1], [[-2, 0], [2, 0]], atol=1e-14)
    # Reversed FITS East axis is an axes transform, not an EVPA sign flip.
    np.testing.assert_allclose(segments.mean(axis=1), 0, atol=1e-14)


def test_circular_band_uses_residuals_and_breaks_seam():
    mean = np.deg2rad([88., 89., 91., 92.])
    samples = (mean[None, :] + np.deg2rad([-1., 0., 1.])[:, None]) % np.pi
    center, lower, upper = plot.circular_band(samples, mean, np.ones(4, bool))
    assert np.isnan(center[2]) and np.isnan(lower[2]) and np.isnan(upper[2])
    np.testing.assert_allclose(center[[0, 1, 3]], [88, 89, -88])
    np.testing.assert_allclose(upper[[0, 1, 3]]-lower[[0, 1, 3]], 1.36)


def test_circular_band_respects_invalid_bins_and_does_not_mutate():
    mean = np.array([0., .1, .2])
    samples = np.tile(mean, (4, 1))
    original = samples.copy()
    curves = plot.circular_band(samples, mean, [True, False, True])
    assert all(np.isnan(c[1]) for c in curves)
    np.testing.assert_array_equal(samples, original)


def test_display_mask_is_not_analysis_mask():
    stokes = np.zeros((2, 2, 4))
    stokes[..., 0] = [[1., .09], [.5, .5]]
    stokes[..., 1] = [[.5, .5], [.09, .2]]
    np.testing.assert_array_equal(plot.display_mask(stokes), [[True, False], [False, True]])
    np.testing.assert_array_equal(plot.display_mask(stokes, np.zeros((2, 2), bool)), False)


def test_incomplete_source_rejected_without_output(tmp_path):
    import pytest
    source = tmp_path / 'partial'
    source.mkdir()
    destination = tmp_path / 'figures'
    with pytest.raises(ValueError, match='incomplete'):
        plot.main(['--ensemble', str(source), '--output-root', str(destination)])
    assert not destination.exists()


def test_render_synthetic_products(tmp_path):
    """Exercise plotting against the declared schema, not actual MC inference."""
    source, out = tmp_path / 'source', tmp_path / 'plots'
    out.mkdir()
    e, n = np.linspace(4, -4, 8), np.linspace(-4, 4, 8)
    axes = {'east_uas': e, 'north_uas': n, 'angle_deg': np.arange(360.)}
    mean = np.linspace(.1, .4, 360)
    samples = mean[None, :] + np.linspace(-.01, .01, 10)[:, None]
    s = np.ones((8, 8, 4))
    s[..., 1] = .2
    s[..., 2] = .1
    records = {}
    for mode in plot.MODES:
        base = source / mode / plot.RESOLUTION
        for day in plot.DAYS:
            directory = base / ('day-'+day)
            directory.mkdir(parents=True)
            np.savez(directory / 'map-summary.npz', mean_stokes=s,
                     mean_rad=np.full((8, 8), .2), valid=np.ones((8, 8), bool))
            np.savez(directory / 'aperture.npz', samples_rad=samples, mean_rad=mean,
                     valid=np.ones(360, bool), angle_deg=np.arange(360.))
            columns = ['total_I', 'm_net', 'm_resolved', 'evpa_rad', 'beta2_real', 'beta2_imag']
            np.savez(directory / 'global.npz', columns=columns,
                     values=np.tile([.6, .1, .2, .3, .1, .01], (10, 1)))
        (base / 'pairs').mkdir()
        for a, b in plot.PAIRS:
            np.savez(base / 'pairs' / f'{a}-{b}.npz', difference_rad=mean,
                     covariance_rad2=np.eye(360)*.001, valid=np.ones(360, bool))
        plot.plot_maps(source, out, axes, mode, records)
    plot.plot_apertures(source, out)
    plot.plot_globals(source, out, records)
    plot.plot_differences(source, out, axes, records)
    assert len(list(out.glob('*.png'))) == 5
    assert all(p.stat().st_size > 1000 for p in out.glob('*.png'))

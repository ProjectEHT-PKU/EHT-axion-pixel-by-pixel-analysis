"""Analysis: limits on c from the EHT ensembles, with the zero-data control.

Inputs: the expected signal of scripts/analysis/expected_signal.py and the fixed-D EHT ensemble
summary of scripts/eht/summarize.py. Feature selection, cloud model, background and the
inference grid/prior come from experiment.yaml (section analysis). The EHT images
are used as reconstructed (the Paper VII Sec. 5.2 alignment is a zero shift here).
D-draw image scatter is a working covariance, not calibrated noise. Both day
pairs share one phase nuisance parameter, marginalized after multiplying likelihoods.
"""
from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import argparse
import csv
import json

from eht_axion.environment import add_data_root_option, fill_defaults
import numpy as np

from eht_axion.analysis.likelihood import MODES, LimitOptions, calculate_observed
from eht_axion.analysis.observations import (
    DAYS, branch_risk, build_observations, covariance_diagnostics,
)
from eht_axion.analysis.results import write_results
from eht_axion.analysis.expected_signal import load_expected_signal
from eht_axion.eht.ensemble import map_evpa
from eht_axion.experiment import (
    ANGULAR_BINS, DEFAULT_PATH, analysis_geometry, physical_scale, read_experiment,
)
from eht_axion.provenance import (
    describe_inputs, manifest_path, require_declared_outputs, stage_directory, verify_stage,
)
from eht_axion.raptor.images import input_cadence_hours


# Paper VII Sec. 5.2 centring = zero shift in this reproduction: no registration.
CENTERING = 'zero shift'


def options_from_config(config, background_sigma_radians):
    """Bind all physical/prior settings; background must come from verified stage."""
    sigma = background_sigma_radians
    if (isinstance(sigma, (bool, np.bool_)) or not np.isscalar(sigma)
            or not np.isreal(sigma) or not np.isfinite(sigma) or sigma < 0):
        raise ValueError('Background sigma must be finite nonnegative radians')
    settings, cloud = config['analysis']['inference'], config['analysis']['cloud']
    return LimitOptions(
        phase_count=settings['phase_count'], credibility=settings['credibility'],
        coupling_points=settings['coupling_points'], coupling_min=settings['coupling_min'],
        coupling_max=settings['coupling_max'], prior=settings['prior'],
        time_conversion=1 / physical_scale(config).tg_hours,
        exposure_hours=settings['exposure_hours'], day_gap_hours=settings['day_gap_hours'],
        pair_separation_hours=settings['pair_separation_hours'],
        n_principal=config['shared']['principal_n'],
        # Template amplitudes are per unit of the RAPTOR axion normalisation.
        reference_norm=config['raptor']['response']['reference_norm'], cloud_coefficient=cloud['coefficient'],
        cloud_pivot_alpha=cloud['pivot_alpha'], cloud_exponent=cloud['exponent'],
        grmhd_sigma_degrees=float(np.rad2deg(sigma)))


def validate_expected_signal(root, config):
    """Verify analysis:expected-signal: made from the RAPTOR images of this experiment with its geometry."""
    from scripts.analysis.expected_signal import image_inputs, stage_config

    record = verify_stage(root)
    if record['stage'] != 'analysis:expected-signal':
        raise ValueError('Unexpected expected-signal stage')
    cubes = [Path(p) for p in record['inputs'] if Path(p).name == 'stokes.npy']
    if not cubes:
        raise ValueError('Expected-signal stage declares no RAPTOR image cube')
    _, inputs, labels = image_inputs(config, cubes[0].parents[2])
    if record['inputs'] != describe_inputs(inputs):
        raise ValueError('Expected signal was not made from the RAPTOR images of this experiment')
    if record['config'] != json.loads(json.dumps(stage_config(config, labels))):
        raise ValueError('Expected signal differs from the selected geometry')
    path = root / 'expected.npz'
    require_declared_outputs(root, [path])
    products = {}
    for dim in ('1d', '2d'):
        values = load_expected_signal(path, dim)
        if not np.array_equal(values[0], config['shared']['alphas']):
            raise ValueError('Expected-signal alpha grid differs from configuration')
        shape = (ANGULAR_BINS,) if dim == '1d' else (config['shared']['image']['shape'],) * 2
        if values[3].shape != shape:
            raise ValueError('Expected-signal spatial grid differs from configuration')
        products[dim] = values
    return products, [manifest_path(root), path]


def save_results(out, products, prepared, config, sigmas, identity):
    """Run all six dimension/mode cases, each with the zero-data control (same C, d = 0)."""
    results, table = {}, []
    sections = dict(observations={}, background={}, cases={})
    for dim in ('1d', '2d'):
        alpha, amp, phase, _ = products[dim]
        selected, samples, coordinates, selection = prepared[dim]
        observed, zero_data, diagnostics = build_observations(
            samples, selected, resultant_min=identity['resultant_min'])
        positions = diagnostics['retained_initial_positions']
        sections['observations'][dim] = dict(
            diagnostics, selection=selection, selected_indices=observed.selected_indices,
            coordinates=np.asarray(coordinates)[positions], means=np.asarray(observed.means),
            standard_deviations=np.asarray(observed.standard_deviations),
            covariances=np.asarray(observed.covariances),
            pair_differences=observed.pair_differences,
            pair_covariances=np.asarray(observed.pair_covariances),
            zero_data_pair_differences=zero_data.pair_differences,
            evpa_samples_radians=np.stack(samples)[:, :, positions])
        sigmas[dim], sections['background'][dim] = selected_background(
            dim, observed.selected_indices, config, sigmas[dim])
        options = options_from_config(config, sigmas[dim])
        for mode in MODES:
            for label, data in (('observed', observed), ('zero-data', zero_data)):
                a, cuts, diagnostic = calculate_observed(
                    data, amp, phase, alphas=alpha, dimension=dim, mode=mode, options=options)
                effective_covariances = diagnostic.get('effective_pair_covariances')
                if effective_covariances is not None:
                    effective_covariances = np.asarray(effective_covariances)
                    diagnostic['working_covariance_diagnostics'] = covariance_diagnostics(
                        effective_covariances)
                    diagnostic['branch_risk'] = branch_risk(data.pair_differences, effective_covariances)
                diagnostic['observed_mean_definition'] = diagnostics['mean_definition']
                key = f'{label}-{mode}-{dim}'
                sections['cases'][key] = dict(diagnostic, log10_c95=cuts)
                results[key] = dict(alphas=a.tolist(), c95=(10**cuts).tolist(),
                                    log10_c95=cuts.tolist(), features=len(observed.selected_indices))
                table.extend(dict(data=label, mode=mode, dimension=dim, alpha=float(x),
                                  log10_c95=float(y), c95=float(10**y)) for x, y in zip(a, cuts))
                print(key, cuts, flush=True)
    with (out / 'limits.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            'data', 'mode', 'dimension', 'alpha', 'log10_c95', 'c95'])
        writer.writeheader()
        writer.writerows(table)
    plot_results(out, results)
    report = dict(identity, results=results,
                  supported_cases=[f'{mode}-{dim}' for dim in ('1d', '2d') for mode in MODES],
                  background_sigma_radians=sigmas,
                  credibility=config['analysis']['inference']['credibility'],
                  response_model_variance_included=False,
                  warnings=['D-draw scatter is working covariance, not calibrated measurement noise.',
                            'Independent days assumed; common generating D does not identify cross-day C.',
                            'Local EVPA Gaussian branch; no daily EVPA alignment or offset removal.',
                            'Finite prior/grid conditional credible bounds, not certified exclusions.',
                            'cor_fixed_background adds a fixed scalar GRMHD background variance, '
                            'not a response-model variance.'])
    write_results(out, report, sections)
    return report


def plot_results(out, results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for ax, dim in zip(axes, ('1d', '2d')):
        for mode, color in zip(MODES, ('#2679b5', '#e58f24', '#467b46')):
            for label, style in (('observed', '-'), ('zero-data', '--')):
                row = results[f'{label}-{mode}-{dim}']
                ax.plot(row['alphas'], row['log10_c95'], style, color=color, marker='o',
                        markersize=3, label=f'{label}: {mode}')
        ax.set(xlabel='Axion alpha', ylabel='log10 c95', title=dim)
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    fig.suptitle('Observed day differences and zero-data control; D-term scatter covariance')
    fig.savefig(out / 'observed-vs-zero-data.png', dpi=160)
    plt.close(fig)


def selected_background(dimension, selected, config, source):
    """Pool lagged no-axion differences on exactly the final likelihood features.

    1D integrates the unchanged theory aperture before selecting angular bins;
    masking sparse 2D pixels before integration would define a different datum.
    No sqrt(2) correction: this sigma describes a day-pair difference directly.
    """
    from eht_axion.evpa import circular_difference
    from eht_axion.analysis.expected_signal import _integrate

    maps = np.load(source['stokes'], mmap_mode='r', allow_pickle=False)
    n = config['shared']['image']['shape']
    if maps.ndim != 5 or maps.shape[1:] != (1, n, n, 4):
        raise ValueError('Background requires processed (time,1,N,N,4) no-axion data')
    lag = config['analysis']['background']['lag_frames']
    if lag >= len(maps):
        raise ValueError('Background lag exceeds the number of axion-off frames')
    operator = None
    if dimension == '1d':
        operator = theory_aperture_operator(config)
    values = []
    for frame in maps:
        frame = frame[0]
        if dimension == '1d':
            theta = _integrate(frame, operator)[selected]
        else:
            pixels = frame.reshape(-1, 4)[selected]
            if np.any((pixels[:, 1] == 0) & (pixels[:, 2] == 0)):
                raise ValueError('Undefined selected background EVPA')
            theta = .5 * np.arctan2(pixels[:, 2], pixels[:, 1])
        if not np.isfinite(theta).all():
            raise ValueError('No finite-only omission of selected background features allowed')
        values.append(theta)
    values = np.asarray(values)
    differences = circular_difference(values[lag:], values[:-lag])
    center = float(differences.mean())
    sigma = float(np.sqrt(np.mean((differences - center)**2)))
    record = dict(source, dimension=dimension, lag=lag, sigma_radians=sigma, mean_radians=center,
                  lag_hours=lag * source['cadence_hours'], selected_indices=np.asarray(selected),
                  method='pooled centered population MLE; fixed scalar pair-difference variance',
                  cross_feature_covariance_used=False, sqrt2_correction=False,
                  differences=differences, spatial_covariance=np.cov(differences, rowvar=False, ddof=0))
    return sigma, record


def theory_aperture_operator(config):
    """Exactly the RAPTOR fit/background bilinear CSR quadrature in Rg units."""
    from eht_axion.analysis.expected_signal import _operator

    settings = analysis_geometry(config)['aperturefit']
    return _operator(config['shared']['image']['shape'], settings['fov'],
                     settings['radial_max'], settings['window_degrees'], settings['angular_bins'])


def matched_aperture_samples(native_maps, operator):
    """Sector statistics of every EHT draw with exactly the operator used for the expected signal.

    The native maps [draw, y, x] are transposed to the template order [x, y] (axis 0
    west, axis 1 north); sky bin k is centred k degrees east of north. Returns, each of
    shape (draws, ANGULAR_BINS): the sector EVPA (the 1D observable), the sector
    intensity and the I-weighted mean polarized brightness (selection statistics).
    """
    from eht_axion.analysis.expected_signal import _brightness, _integrate

    evpa = np.empty((len(native_maps), ANGULAR_BINS))
    intensity, polarized = np.empty_like(evpa), np.empty_like(evpa)
    for k, frame in enumerate(native_maps):
        frame = frame.transpose(1, 0, 2)
        evpa[k] = _integrate(frame, operator)
        intensity[k], polarized[k] = _brightness(frame, operator)
    return evpa, intensity, polarized


def prepare_limit_inputs(ensemble, expected, config, products):
    """Read the EHT ensemble (scripts/eht/summarize.py) and select the likelihood features.

    Native[y,x] -> template[x,y], physical (-E,N). Sector EVPAs and sector selection
    statistics use the same quadrature operator as the expected signal. The ensemble holds the
    blurred EHT images (shared.image.blur_fwhm_uas, the same blur as the simulated images).
    """
    features = config['analysis']['features']
    stride = features['pixel_stride']
    record = verify_stage(ensemble)
    if record['stage'] != 'eht:ensemble':
        raise ValueError('Unexpected EHT ensemble stage')
    report_path = ensemble / 'report.json'
    require_declared_outputs(ensemble, [report_path])
    report = json.loads(report_path.read_text())
    draws = config['eht_imaging']['uncertainty_draws']
    if (report.get('schema') != 'eht-ensemble'
            or report.get('days') != [str(d) for d in DAYS] or report.get('draws_per_day') != draws
            or record['config'].get('draw_count') != draws * len(DAYS)
            or report.get('angle_convention') != 'EVPA radians modulo pi, east of north'
            or report.get('map_order') !=
            'day,draw,y,x,IQUV; native FITS E decreases along x, N increases along y'):
        raise ValueError('Complete audited fixed-D ensemble contract required')
    axes_path = ensemble / 'axes.npz'
    require_declared_outputs(ensemble, [axes_path])
    with np.load(axes_path, allow_pickle=False) as z:
        east_axis, north_axis = z['east_uas'], z['north_uas']
    geometry = analysis_geometry(config)
    # Native maps are [y, x] with E decreasing along x; the expected signal is [x, y] on (-E, N).
    n = config['shared']['image']['shape']
    if east_axis.shape != (n,) or north_axis.shape != (n,):
        raise ValueError('Ensemble grid size differs from shared.image.shape')
    east = np.broadcast_to(np.asarray(east_axis, dtype=float)[:, None], (n, n))
    north = np.broadcast_to(np.asarray(north_axis, dtype=float)[None, :], (n, n))
    if not (np.isfinite(east).all() and np.isfinite(north).all()):
        raise ValueError('Finite physical template-grid east/north coordinates required')
    grid = (np.arange(n) - (n - 1) / 2) * config['shared']['image']['fov_uas'] / n
    if (not np.allclose(east, -np.broadcast_to(grid[:, None], (n, n)), atol=1e-12, rtol=0)
            or not np.allclose(north, np.broadcast_to(grid[None, :], (n, n)), atol=1e-12, rtol=0)):
        raise ValueError('Ensemble sky axes differ from the shared sky grid (-E,N)')
    stokes_path = ensemble / 'stokes.npy'
    require_declared_outputs(ensemble, [stokes_path])
    inputs = [manifest_path(ensemble), report_path, axes_path, stokes_path]
    ensemble_maps = np.load(stokes_path, mmap_mode='r', allow_pickle=False)
    if ensemble_maps.shape != (len(DAYS), draws, n, n, 4):
        raise ValueError('Ensemble Stokes sample shape mismatch')
    data, masks, angular_masks = [], [], []
    operator = theory_aperture_operator(config)
    for maps in ensemble_maps:
        aperture, sector_intensity, sector_polarized = matched_aperture_samples(maps, operator)
        # Draw-mean sector statistics; both thresholds always apply
        # (polarization_fraction 0 selects by brightness only).
        intensity_1d, polarized_1d = sector_intensity.mean(axis=0), sector_polarized.mean(axis=0)
        if not (np.isfinite(intensity_1d).all() and np.isfinite(polarized_1d).all()):
            raise ValueError('Nonfinite sector brightness (a sector without positive intensity)')
        angular_masks.append((intensity_1d >= features['intensity_fraction'] * intensity_1d.max())
                             & (polarized_1d >= features['polarization_fraction'] * polarized_1d.max()))
        evpa_maps = map_evpa(np.asarray(maps))
        mean = np.mean(maps, axis=0).transpose(1, 0, 2)
        if not np.isfinite(mean).all():
            raise ValueError('Nonfinite ensemble mean Stokes')
        intensity, polarized = mean[..., 0], np.hypot(mean[..., 1], mean[..., 2])
        if intensity.max() <= 0 or polarized.max() <= 0:
            raise ValueError('Positive intensity and polarization required')
        mask = ((intensity >= features['intensity_fraction'] * intensity.max())
                & (polarized >= features['polarization_fraction'] * polarized.max()))
        mask &= np.isfinite(evpa_maps).all(axis=0).T
        masks.append(mask)
        data.append((evpa_maps, aperture))
    row, col = np.indices((n, n))
    radius = np.hypot(east, north)
    physical_mask = ((radius >= features['min_radius_rg'] * physical_scale(config).angular_rg_uas)
        & (radius <= geometry['metadata']['radial_max_uas']))
    common = np.logical_and.reduce(masks) & physical_mask
    prepared = {}
    for dim in ('1d', '2d'):
        if dim == '1d':
            finite = np.logical_and.reduce([np.isfinite(v).all(axis=0) for _, v in data])
            sector = (np.arange(ANGULAR_BINS) % features['sector_spacing_degrees']
                      == features['sector_offset_degrees'])
            selected = np.flatnonzero(products[dim][3] & finite
                & np.logical_and.reduce(angular_masks) & sector)
            values = [np.asarray(aperture[:, selected], dtype=float) for _, aperture in data]
            coords = selected[:, None].astype(float)
        else:
            valid = products[dim][3] & common & (row % stride == 0) & (col % stride == 0)
            selected = np.flatnonzero(valid)
            r, c = np.unravel_index(selected, (n, n))
            native_r, native_c = c, r  # template [x, y] -> native [y, x]
            values = []
            for evpa_maps, _ in data:
                values.append(np.asarray(evpa_maps[:, native_r, native_c], dtype=float))
            coords = np.column_stack((east.ravel()[selected], north.ravel()[selected]))
        if not len(selected):
            raise ValueError(f'No initial {dim} features')
        prepared[dim] = (selected, values, coords, dict(
            intensity_fraction=features['intensity_fraction'],
            polarization_fraction=features['polarization_fraction'],
            min_radius_rg=features['min_radius_rg'], pixel_stride=stride,
            angular_spacing=features['sector_spacing_degrees'],
            angular_offset=features['sector_offset_degrees'],
            coordinates='east/north uas; 1d template bin degrees',
            no_evpa_alignment=True, adapter=dict(centering=CENTERING, resolution='blurred',
                spatial='native[y,x] -> template[x,y]', template_axes='(-E,N)',
                angular='sky bin k centred k degrees east of north',
                sector_definition='aperture_operator quadrature on transposed per-draw Stokes, as for the '
                                  'expected signal: EVPA 0.5 arg sum W I (Q+iU); selection by draw means of '
                                  'sum W I and sum W I |P| / sum W I',
                quadrature=geometry['aperturefit'], qu_rotation='none')))
    stokes_path = expected / 'background.npy'
    require_declared_outputs(expected, [stokes_path])
    # Frame spacing from the GRMHD snapshot times recorded with the background image cube.
    (cube,) = [Path(p) for p in verify_stage(expected)['inputs']
               if Path(p).name == 'stokes.npy' and Path(p).parent.name == 'background']
    cadence = input_cadence_hours(cube, config)
    background_inputs = [stokes_path, cube.parent / 'generated.json']
    sigmas = {dim: dict(stokes=str(stokes_path), cadence_hours=cadence) for dim in ('1d', '2d')}
    return prepared, inputs, dict(stage=record['stage'], centering=CENTERING,
        background_selection='pooled scalar on final selected features, no response variance'), \
        sigmas, background_inputs


def run(args):
    config = read_experiment(args.config)
    features = config['analysis']['features']
    if args.output_root.exists():
        raise FileExistsError('New output root required')
    products, inputs = validate_expected_signal(args.expected_signal.resolve(), config)
    prepared, source_inputs, source_identity, sigmas, background_inputs = prepare_limit_inputs(
        args.ensemble.resolve(), args.expected_signal.resolve(), config, products)
    identity = dict(schema='analysis-limits',
                    scientific_profile=config,
                    resultant_min=features['resultant_min'], centering=CENTERING,
                    pixel_stride=features['pixel_stride'],
                    ensemble=source_identity, covariance='D-draw population scatter; independent days',
                    scope='credible limits conditional on grid and prior; D-term scatter as covariance')
    inputs += source_inputs + background_inputs
    with stage_directory(args.output_root.resolve(), 'analysis:limits',
                         identity, inputs) as out:
        if out is not None:
            return save_results(out, products, prepared, config, sigmas, identity)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=DEFAULT_PATH, help='default: configs/experiment.yaml')
    p.add_argument('--ensemble', type=Path, help='scripts/eht/summarize.py stage (default: eht/ensemble)')
    p.add_argument('--expected-signal', type=Path,
                   help='scripts/analysis/expected_signal.py stage (default: analysis/expected-signal)')
    p.add_argument('--output-root', type=Path, help='new directory (default: analysis/limits)')
    add_data_root_option(p)
    args = p.parse_args(argv)
    fill_defaults(p, args, {'ensemble': ('eht.ensemble', {}),
                            'expected_signal': ('analysis.expected_signal', {}),
                            'output_root': ('analysis.limits', {})})
    run(args)


if __name__ == '__main__':
    main()

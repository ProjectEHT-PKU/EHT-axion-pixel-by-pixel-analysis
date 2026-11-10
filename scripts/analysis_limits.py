"""Analysis: limits on c from the EHT ensembles, with the zero-data control.

Inputs: the templates root of analysis_templates.py and the fixed-D EHT ensemble
summary of eht_summarize.py. Feature selection, cloud model, background and the
inference grid/prior come from experiment.yaml (section analysis). The EHT images
are used as reconstructed (the Paper VII Sec. 5.2 alignment is a zero shift here).
D-draw image scatter is a working covariance, not calibrated noise. Both day
pairs share one phase nuisance parameter, marginalized after multiplying likelihoods.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
from _data_root import add_data_root_option, fill_defaults
import numpy as np

from eht_axion.analysis.likelihood import MODES, LimitOptions, calculate_observed
from eht_axion.analysis.observations import (
    DAYS, branch_risk, build_observations, covariance_diagnostics,
)
from eht_axion.analysis.templates import load_templates
from eht_axion.experiment import (
    ANGULAR_BINS, DEFAULT_PATH, analysis_geometry, physical_scale, read_experiment,
)
from eht_axion.provenance import (
    manifest_path, require_declared_outputs, stage_directory, verification_session, verify_stage,
)
from eht_axion.raptor.delivery import input_cadence_hours, validate_processed_source
from eht_axion.raptor.layout import alpha_label


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


def validate_templates(root, config, seen):
    """Verify the template stages: delivered-data identity for all alphas and the analysis geometry."""
    geometry = analysis_geometry(config)
    for alpha in config['shared']['alphas']:
        label = alpha_label(alpha)
        processed = root / 'processed' / label
        record = verify_stage(processed, recursive=True, _seen=seen)
        validate_processed_source(record, config, alpha=alpha,
                                  resolution=config['raptor']['response']['resolution'], verified=seen)
        if record['config'] != geometry['image_preprocess']:
            raise ValueError('Theory preprocessing differs from selected geometry')
        record = verify_stage(root / 'response' / label, recursive=True, _seen=seen)
        if {k: v for k, v in record['config'].items() if k != 'snapshot_labels'} != geometry['aperturefit']:
            raise ValueError('Theory aperture differs from selected geometry')
    products, inputs = {}, []
    for dim in ('1d', '2d'):
        directory = root / 'templates' / dim
        record = verify_stage(directory, recursive=True, _seen=seen)
        if record['stage'] != 'analysis:templates':
            raise ValueError('Unexpected theory template stage')
        path = directory / 'templates.npz'
        require_declared_outputs(directory, [path])
        values = load_templates(path, dim)
        if not np.array_equal(values[0], config['shared']['alphas']):
            raise ValueError('Theory alpha grid differs from configuration')
        expected = (ANGULAR_BINS,) if dim == '1d' else (config['shared']['image']['shape'],) * 2
        if values[3].shape != expected:
            raise ValueError('Template spatial grid differs from configuration')
        products[dim] = values
        inputs += [manifest_path(directory), path]
    return products, inputs


def save_results(out, products, prepared, config, sigmas, identity):
    """Run all six dimension/mode cases, each with the zero-data control (same C, d = 0)."""
    results, table = {}, []
    for dim in ('1d', '2d'):
        alpha, amp, phase, _ = products[dim]
        selected, samples, coordinates, selection = prepared[dim]
        observed, zero_data, diagnostics = build_observations(
            samples, selected, resultant_min=identity['resultant_min'])
        positions = diagnostics['retained_initial_positions']
        np.savez(out / f'observations-{dim}.npz', selected_indices=observed.selected_indices,
                 coordinates=np.asarray(coordinates)[positions], means=observed.means,
                 standard_deviations=observed.standard_deviations, covariances=observed.covariances,
                 pair_differences=observed.pair_differences, pair_covariances=observed.pair_covariances,
                 zero_data_pair_differences=zero_data.pair_differences, resultant=diagnostics['resultant'],
                 sample_counts=diagnostics['sample_counts'],
                 local_branch_mean_minus_circular_radians=
                 diagnostics['local_branch_mean_minus_circular_radians'],
                 evpa_samples_radians=np.stack(samples)[:, :, positions])
        (out / f'observations-{dim}.json').write_text(json.dumps(
            dict(diagnostics, selection=selection), indent=2, allow_nan=False) + '\n')
        sigmas[dim] = selected_background(out, dim, observed.selected_indices, config, sigmas[dim])
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
                np.savez(out / f'{key}.npz', alphas=a, log10_c95=cuts, c95=10**cuts,
                         couplings=diagnostic['couplings'], log_likelihood=diagnostic['log_likelihood'],
                         posterior_density=diagnostic['posterior_density'], cdf=diagnostic['cdf'])
                (out / f'{key}.json').write_text(json.dumps(diagnostic, indent=2, allow_nan=False)+'\n')
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
    (out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
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


def selected_background(out, dimension, selected, config, source):
    """Pool lagged no-axion differences on exactly the final likelihood features.

    1D integrates the unchanged theory aperture before selecting angular bins;
    masking sparse 2D pixels before integration would define a different datum.
    No sqrt(2) correction: this sigma describes a day-pair difference directly.
    """
    from eht_axion.evpa import circular_difference
    from eht_axion.analysis.templates import _integrate

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
    np.savez(out / f'background-selected-{dimension}.npz', selected_indices=selected,
             differences=differences, sigma_radians=sigma, mean_radians=center,
             spatial_covariance=np.cov(differences, rowvar=False, ddof=0))
    metadata = dict(source, dimension=dimension, lag=lag, sigma_radians=sigma, mean_radians=center,
                    lag_hours=lag * source['cadence_hours'], selected_indices=selected.tolist(),
                    method='pooled centered population MLE; fixed scalar pair-difference variance',
                    cross_feature_covariance_used=False, sqrt2_correction=False)
    (out / f'background-selected-{dimension}.json').write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + '\n')
    return sigma


def theory_aperture_operator(config):
    """Exactly the RAPTOR fit/background bilinear CSR quadrature in Rg units."""
    from eht_axion.analysis.templates import _operator

    settings = analysis_geometry(config)['aperturefit']
    return _operator(config['shared']['image']['shape'], settings['fov'],
                     settings['radial_max'], settings['window_degrees'], settings['angular_bins'])


def matched_aperture_samples(native_maps, operator):
    """Sector statistics of every EHT draw with exactly the operator used for the templates.

    The native maps [draw, y, x] are transposed to the template order [x, y] (axis 0
    west, axis 1 north); sky bin k is centred k degrees east of north. Returns, each of
    shape (draws, ANGULAR_BINS): the sector EVPA (the 1D observable), the sector
    intensity and the I-weighted mean polarized brightness (selection statistics).
    """
    from eht_axion.analysis.templates import _brightness, _integrate

    evpa = np.empty((len(native_maps), ANGULAR_BINS))
    intensity, polarized = np.empty_like(evpa), np.empty_like(evpa)
    for k, frame in enumerate(native_maps):
        frame = frame.transpose(1, 0, 2)
        evpa[k] = _integrate(frame, operator)
        intensity[k], polarized[k] = _brightness(frame, operator)
    return evpa, intensity, polarized


def prepare_limit_inputs(ensemble, templates, config, products, *, verified):
    """Read the EHT ensemble (eht_summarize.py) and select the likelihood features.

    Native[y,x] -> template[x,y], physical (-E,N). Sector EVPAs and sector selection
    statistics use the same quadrature operator as the templates. The ensemble holds the
    blurred EHT images (shared.image.blur_fwhm_uas, the same blur as the simulated images).
    """
    features = config['analysis']['features']
    stride = features['pixel_stride']
    record = verify_stage(ensemble, recursive=True, _seen=verified)
    if record['stage'] != 'eht:ensemble':
        raise ValueError('Unexpected EHT ensemble stage')
    report_path = ensemble / 'report.json'
    require_declared_outputs(ensemble, [report_path])
    report = json.loads(report_path.read_text())
    draws = config['eht_imaging']['uncertainty_draws']
    if (report.get('schema') != 'eht-ensemble-v2'
            or report.get('days') != [str(d) for d in DAYS] or report.get('draws_per_day') != draws
            or record['config'].get('draw_count') != draws * len(DAYS)
            or report.get('angle_convention') != 'EVPA radians modulo pi, east of north'
            or report.get('map_order') !=
            'draw,y,x,IQUV; native FITS E decreases along x, N increases along y'):
        raise ValueError('Complete audited fixed-D ensemble contract required')
    axes_path = ensemble / 'axes.npz'
    require_declared_outputs(ensemble, [axes_path])
    with np.load(axes_path, allow_pickle=False) as z:
        east_axis, north_axis = z['east_uas'], z['north_uas']
    days = {str(d): dict(stokes=f'day-{d}/stokes.npy', evpa_2d=f'day-{d}/evpa.npy') for d in DAYS}
    geometry = analysis_geometry(config)
    # Native maps are [y, x] with E decreasing along x; templates are [x, y] on (-E, N).
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
    inputs = [manifest_path(ensemble), report_path, axes_path]
    data, masks, angular_masks = [], [], []
    operator = theory_aperture_operator(config)
    for day in DAYS:
        entry = days[str(day)]
        paths = [ensemble / entry[key] for key in ('stokes', 'evpa_2d')]
        require_declared_outputs(ensemble, paths)
        inputs.extend(paths)
        maps = np.load(paths[0], mmap_mode='r', allow_pickle=False)
        if maps.ndim != 4 or maps.shape[1:] != (n, n, 4) or len(maps) < 2:
            raise ValueError('Ensemble Stokes sample shape mismatch')
        aperture, sector_intensity, sector_polarized = matched_aperture_samples(maps, operator)
        # Draw-mean sector statistics; both thresholds always apply
        # (polarization_fraction 0 selects by brightness only).
        intensity_1d, polarized_1d = sector_intensity.mean(axis=0), sector_polarized.mean(axis=0)
        if not (np.isfinite(intensity_1d).all() and np.isfinite(polarized_1d).all()):
            raise ValueError('Nonfinite sector brightness (a sector without positive intensity)')
        angular_masks.append((intensity_1d >= features['intensity_fraction'] * intensity_1d.max())
                             & (polarized_1d >= features['polarization_fraction'] * polarized_1d.max()))
        evpa_maps = np.load(paths[1], mmap_mode='r', allow_pickle=False)
        if evpa_maps.shape != maps.shape[:-1]:
            raise ValueError('Ensemble map EVPA shape mismatch')
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
                                  'templates: EVPA 0.5 arg sum W I (Q+iU); selection by draw means of '
                                  'sum W I and sum W I |P| / sum W I',
                quadrature=geometry['aperturefit'], qu_rotation='none')))
    processed = templates / 'background'
    processed_record = verify_stage(processed, recursive=True, _seen=verified)
    if processed_record['stage'] != 'analysis:background':
        raise ValueError('Unexpected background preprocessing stage')
    preprocessing = {k: v for k, v in processed_record['config'].items() if k != 'state_indices'}
    if (preprocessing != geometry['image_preprocess']
            or processed_record['config'].get('state_indices') != [0]):
        raise ValueError('Background preprocessing geometry/state mismatch')
    stokes_path = processed / 'stokes.npy'
    require_declared_outputs(processed, [stokes_path])
    # Frame spacing from the GRMHD snapshot times recorded with the delivered cubes.
    validate_processed_source(processed_record, config, alpha=None,
                              resolution=config['raptor']['background']['resolution'], verified=verified)
    delivered = [Path(e['path']) for e in processed_record['inputs'] if Path(e['path']).suffix == '.npy']
    cadence = input_cadence_hours(delivered, config)
    background_inputs = [manifest_path(processed), stokes_path,
                         *(path.parent / 'generated.json' for path in delivered)]
    sigmas = {dim: dict(stokes=str(stokes_path), cadence_hours=cadence) for dim in ('1d', '2d')}
    return prepared, inputs, dict(stage=record['stage'], centering=CENTERING,
        background_selection='pooled scalar on final selected features, no response variance'), \
        sigmas, background_inputs


def run(args):
    config = read_experiment(args.config)
    features = config['analysis']['features']
    if args.output_root.exists():
        raise FileExistsError('New output root required')
    with verification_session():
        seen = set()
        products, inputs = validate_templates(args.templates.resolve(), config, seen)
        prepared, source_inputs, source_identity, sigmas, background_inputs = prepare_limit_inputs(
            args.ensemble.resolve(), args.templates.resolve(), config, products, verified=seen)
        identity = dict(schema='analysis-limits-v2',
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
    p.add_argument('--config', type=Path, default=DEFAULT_PATH, help='default: experiment.yaml')
    p.add_argument('--ensemble', type=Path, help='eht_summarize.py stage (default: eht/ensemble)')
    p.add_argument('--templates', type=Path, help='analysis_templates.py root (default: analysis/templates)')
    p.add_argument('--output-root', type=Path, help='new directory (default: analysis/limits)')
    add_data_root_option(p)
    args = p.parse_args(argv)
    fill_defaults(p, args, {'ensemble': ('eht_ensemble', {}), 'templates': ('analysis_templates', {}),
                            'output_root': ('analysis_limits', {})})
    run(args)


if __name__ == '__main__':
    main()

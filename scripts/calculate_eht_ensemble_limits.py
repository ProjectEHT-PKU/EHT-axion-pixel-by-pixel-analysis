"""New-root EHT120 observed / matched-NULL fixed-working-covariance limits.

Not the canonical NULL pipeline. D-draw image scatter is not calibrated noise.
No coupling-squared signal-response variance is invented. Both pairs share one
phase nuisance parameter, marginalized only after multiplying likelihoods.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.limits import LimitOptions
from eht_axion.pipeline.common import (
    manifest_path, require_declared_outputs, stage_directory, verification_session, verify_stage,
)
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.constants import PhysicalScale, analysis_geometry
from eht_axion.pipeline.layout import alpha_label
from eht_axion.pipeline.products import load_templates
from eht_axion.producer.delivery import validate_processed_source
from eht_axion.research.ensemble_limits import (
    DAYS, MODES, branch_risk, build_observations, covariance_diagnostics,
)
from eht_axion.research.observed import calculate_observed


def options_from_config(config, background_sigma_radians):
    """Bind all physical/prior settings; background must come from verified stage."""
    sigma = background_sigma_radians
    if (isinstance(sigma, (bool, np.bool_)) or not np.isscalar(sigma)
            or not np.isreal(sigma) or not np.isfinite(sigma) or sigma < 0):
        raise ValueError('Background sigma must be finite nonnegative radians')
    settings, cloud = config['inference'], config['cloud']
    return LimitOptions(
        legacy=False, phase_count=settings['phase_count'],
        coupling_points=settings['coupling_points'], coupling_min=settings['coupling_min'],
        coupling_max=settings['coupling_max'], prior=settings['prior'],
        time_conversion=1 / PhysicalScale(**config['physical']).tg_hours,
        exposure_hours=settings['exposure_hours'], day_gap_hours=settings['day_gap_hours'],
        pair_separation_hours=settings['pair_separation_hours'], n_principal=settings['principal_n'],
        reference_norm=cloud['reference_norm'], cloud_coefficient=cloud['coefficient'],
        cloud_pivot_alpha=cloud['pivot_alpha'], cloud_exponent=cloud['exponent'],
        grmhd_sigma_degrees=float(np.rad2deg(sigma)), std_covariance='covariance',
        duplicate_phase_endpoint=False, combine_epochs_by_product=True,
        posterior_method='trapezoid')


def validate_templates(root, config, seen):
    """Strict producer binding including all eight alpha identities and geometry."""
    geometry = analysis_geometry(config)
    for alpha in config['data']['alphas']:
        label = alpha_label(alpha)
        processed = root / 'raptor/processed' / label
        record = verify_stage(processed, recursive=True, allow_imported=False, _seen=seen)
        validate_processed_source(record, config, alpha=alpha,
                                  resolution=config['producer']['response']['resolution'], verified=seen)
        if record['config'] != geometry['image_preprocess']:
            raise ValueError('Theory preprocessing differs from selected geometry')
        record = verify_stage(root / 'raptor/response' / label, recursive=True,
                              allow_imported=False, _seen=seen)
        if {k: v for k, v in record['config'].items() if k != 'snapshot_labels'} != geometry['aperturefit']:
            raise ValueError('Theory aperture differs from selected geometry')
    products, inputs = {}, []
    for dim in ('1d', '2d'):
        directory = root / 'raptor/templates' / dim
        record = verify_stage(directory, recursive=True, allow_imported=False, _seen=seen)
        if record['stage'] != 'raptor:summarize':
            raise ValueError('Unexpected theory template stage')
        path = directory / 'templates.npz'
        require_declared_outputs(directory, [path])
        values = load_templates(path, dim)
        if not np.array_equal(values[0], config['data']['alphas']):
            raise ValueError('Theory alpha grid differs from configuration')
        expected = (360,) if dim == '1d' else (config['image']['final_shape'],) * 2
        if values[3].shape != expected:
            raise ValueError('Template spatial grid differs from configuration')
        products[dim] = values
        inputs += [manifest_path(directory), path]
    return products, inputs


def save_results(out, products, prepared, config, sigmas, identity):
    """Run all six dimension/mode cases and their identical-C NULL controls."""
    results, table = {}, []
    for dim in ('1d', '2d'):
        alpha, amp, phase, _ = products[dim]
        selected, samples, coordinates, selection = prepared[dim]
        observed, null, diagnostics = build_observations(
            samples, selected, resultant_min=identity['resultant_min'])
        positions = diagnostics['retained_initial_positions']
        np.savez(out / f'observations-{dim}.npz', selected_indices=observed.selected_indices,
                 coordinates=np.asarray(coordinates)[positions], means=observed.means,
                 standard_deviations=observed.standard_deviations, covariances=observed.covariances,
                 pair_differences=observed.pair_differences, pair_covariances=observed.pair_covariances,
                 null_pair_differences=null.pair_differences, resultant=diagnostics['resultant'],
                 sample_counts=diagnostics['sample_counts'],
                 local_branch_mean_minus_circular_radians=
                 diagnostics['local_branch_mean_minus_circular_radians'],
                 evpa_samples_radians=np.stack(samples)[:, :, positions])
        (out / f'observations-{dim}.json').write_text(json.dumps(
            dict(diagnostics, selection=selection), indent=2, allow_nan=False) + '\n')
        if isinstance(sigmas[dim], dict):
            sigmas[dim] = selected_background(out, dim, observed.selected_indices,
                                             config, sigmas[dim])
        options = options_from_config(config, sigmas[dim])
        for mode in MODES:
            for label, data in (('observed', observed), ('null-comparison', null)):
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
                table.extend(dict(center=label, mode=mode, dimension=dim, alpha=float(x),
                                  log10_c95=float(y), c95=float(10**y)) for x, y in zip(a, cuts))
                print(key, cuts, flush=True)
    with (out / 'limits.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            'center', 'mode', 'dimension', 'alpha', 'log10_c95', 'c95'])
        writer.writeheader()
        writer.writerows(table)
    plot_results(out, results)
    report = dict(identity, results=results,
                  supported_cases=[f'{mode}-{dim}' for dim in ('1d', '2d') for mode in MODES],
                  background_sigma_radians=sigmas, confidence=.95,
                  response_model_variance_included=False, canonical_null_modified=False,
                  warnings=['D-draw scatter is working covariance, not calibrated measurement noise.',
                            'Independent days assumed; common generating D does not identify cross-day C.',
                            'Local EVPA Gaussian branch; no daily EVPA alignment or offset removal.',
                            'Finite prior/grid conditional credible bounds, not certified exclusions.',
                            'cor_fixed_background is fixed scalar GRMHD noise, NOT canonical cor_grmhd '
                            'signal variance.'])
    (out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


def plot_results(out, results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for ax, dim in zip(axes, ('1d', '2d')):
        for mode, color in zip(MODES, ('#2679b5', '#e58f24', '#467b46')):
            for label, style in (('observed', '-'), ('null-comparison', '--')):
                row = results[f'{label}-{mode}-{dim}']
                ax.plot(row['alphas'], row['log10_c95'], style, color=color, marker='o',
                        markersize=3, label=f'{label}: {mode}')
        ax.set(xlabel='Axion alpha', ylabel='log10 c95', title=dim)
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    fig.suptitle('RESEARCH: actual day differences / matched NULL; working D-scatter covariance')
    fig.savefig(out / 'observed-vs-null.png', dpi=160)
    plt.close(fig)


def selected_background(out, dimension, selected, config, source):
    """Pool lagged no-axion differences on exactly the final likelihood features.

    1D integrates the unchanged theory aperture before selecting angular bins;
    masking sparse 2D pixels before integration would define a different datum.
    No sqrt(2) correction: this sigma describes a day-pair difference directly.
    """
    from eht_axion.angles import circular_difference
    from eht_axion.pipeline.raptor import _integrate

    geometry = analysis_geometry(config)
    maps = np.load(source['stokes'], mmap_mode='r', allow_pickle=False)
    n = config['image']['final_shape']
    if maps.ndim != 5 or maps.shape[1:] != (1, n, n, 4):
        raise ValueError('Background requires processed (time,1,N,N,4) no-axion data')
    lag = source['lag']
    if (isinstance(lag, bool) or not isinstance(lag, int) or lag < 1 or lag >= len(maps)
            or lag != config['background']['lag_frames']):
        raise ValueError('Background lag differs from configuration')
    aperture = geometry['aperturefit']
    operator = None
    if dimension == '1d':
        operator = theory_aperture_operator(config)
    values = []
    for frame in maps:
        frame = frame[0]
        if dimension == '1d':
            theta = _integrate(frame, operator, aperture['weighting'])[selected]
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
    metadata = dict(source, dimension=dimension, sigma_radians=sigma, mean_radians=center,
                    lag_hours=lag * source['cadence_hours'], selected_indices=selected.tolist(),
                    method='pooled centered population MLE; fixed scalar pair-difference variance',
                    cross_feature_covariance_used=False, sqrt2_correction=False,
                    unmasked_sigma_used=False)
    (out / f'background-selected-{dimension}.json').write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + '\n')
    return sigma


def theory_aperture_operator(config):
    """Exactly the RAPTOR fit/background bilinear CSR quadrature in Rg units."""
    from eht_axion.pipeline.raptor import _operator

    settings = analysis_geometry(config)['aperturefit']
    return _operator(config['image']['final_shape'], settings['fov'], None,
                     settings['radial_max'], settings['window_degrees'], settings['angular_bins'])


def matched_aperture_samples(native_maps, config, operator):
    """Recompute native EHT draws with exactly the theory aperture observable.

    Transpose spatial axes only. RAPTOR's _integrate applies its internal -90
    bin conversion once; the returned bins are North0/East90. Native producer
    pixel-sector curves are separate diagnostics, never inference samples.
    """
    from eht_axion.pipeline.raptor import _integrate

    curves = np.empty((len(native_maps), 360))
    weighting = config['aperture']['weighting']
    for k, frame in enumerate(native_maps):
        curves[k] = _integrate(frame.transpose(1, 0, 2), operator, weighting)
    return curves


def prepare_limit_inputs(ensemble, templates, config, products, *, centering,
                         pixel_stride, mask_mode, verified):
    """Bind actual fixed-D producer schema using the independently audited axes.

    Native[y,x] -> template[x,y], physical (-E,N). Recompute 1D using theory
    CSR quadrature; native sector apertures are diagnostic only. Read display19 products,
    matching the declared additional19uas theory kernel approximation.
    """
    record = verify_stage(ensemble, recursive=True, allow_imported=False, _seen=verified)
    if record['stage'] != 'research:eht-fixed-D-ensemble':
        raise ValueError('Unexpected EHT120 ensemble stage')
    report_path = ensemble / 'report.json'
    require_declared_outputs(ensemble, [report_path])
    report = json.loads(report_path.read_text())
    if (report.get('schema') != 'eht-fixed-D-ensemble-summary-v1'
            or report.get('days') != [str(d) for d in DAYS] or report.get('draws_per_day') != 1000
            or record['config'].get('draw_count') != 4000
            or report.get('angle_convention') !=
            'EVPA radians modulo pi, E of N; aperture angles degrees E of N; no roll'
            or report.get('map_order') !=
            'draw,y,x,IQUV; native FITS E decreases along x, N increases along y'):
        raise ValueError('Complete audited fixed-D ensemble contract required')
    axes_path = ensemble / 'axes.npz'
    require_declared_outputs(ensemble, [axes_path])
    with np.load(axes_path, allow_pickle=False) as z:
        east_axis, north_axis, angles = z['east_uas'], z['north_uas'], z['angle_deg']
    if not np.array_equal(angles, np.arange(360)):
        raise ValueError('Audited North0/East90 aperture bins require identity mapping, no extra roll')
    mode = 'native' if centering == 'native' else 'icentroid'
    adapter = dict(orientation_verified=True, shape=64, fov_uas=120.,
        radial_max_uas=report['rmax_uas'], window_degrees=report['aperture_full_width_deg'],
        weighting='intensity', native_to_template_axes=[1, 0],
        template_to_native_bins=list(range(360)), east_uas=np.broadcast_to(east_axis[:, None], (64, 64)),
        north_uas=np.broadcast_to(north_axis[None, :], (64, 64)), days={
            str(d): dict(stokes=f'{mode}/display19/day-{d}/stokes.npy',
                        evpa_1d=f'{mode}/display19/day-{d}/aperture.npz',
                        evpa_2d=f'{mode}/display19/day-{d}/evpa.npy', sample_count=1000) for d in DAYS})
    geometry = analysis_geometry(config)
    if (adapter.get('shape') != config['image']['final_shape']
            or adapter.get('fov_uas') != geometry['metadata']['final_fov_uas']
            or adapter.get('radial_max_uas') != geometry['metadata']['radial_max_uas']
            or adapter.get('window_degrees') != config['aperture']['window_degrees']
            or adapter.get('weighting') != config['aperture']['weighting']):
        raise ValueError('Ensemble sky/aperture geometry differs from templates')
    axes = adapter.get('native_to_template_axes')
    bins = np.asarray(adapter.get('template_to_native_bins'))
    if axes != [1, 0]:
        raise ValueError('Audited EHT native[y,x] to template[x,y] transpose required')
    if (bins.shape != (360,) or bins.dtype.kind not in 'iu'
            or not np.array_equal(bins, np.arange(360))):
        raise ValueError('Audited North0/East90 aperture bins require identity mapping, no extra roll')
    n = config['image']['final_shape']
    east, north = (np.asarray(adapter.get(k), dtype=float) for k in ('east_uas', 'north_uas'))
    if any(v.shape != (n, n) or not np.isfinite(v).all() for v in (east, north)):
        raise ValueError('Finite physical template-grid east/north coordinates required')
    grid = (np.arange(n) - (n - 1) / 2) * adapter['fov_uas'] / n
    if (not np.allclose(east, -np.broadcast_to(grid[:, None], (n, n)), atol=1e-12, rtol=0)
            or not np.allclose(north, np.broadcast_to(grid[None, :], (n, n)), atol=1e-12, rtol=0)):
        raise ValueError('Adapter physical coordinates differ from audited template axes (-E,N)')
    inputs = [manifest_path(ensemble), report_path, axes_path]
    data, masks, angular_masks = [], [], []
    operator = theory_aperture_operator(config)
    for day in DAYS:
        entry = adapter['days'][str(day)]
        paths = [ensemble / entry[key] for key in ('stokes', 'evpa_1d', 'evpa_2d')]
        require_declared_outputs(ensemble, paths)
        inputs.extend(paths)
        maps = np.load(paths[0], mmap_mode='r', allow_pickle=False)
        aperture = matched_aperture_samples(maps, config, operator)
        with np.load(paths[1], allow_pickle=False) as z:
            # Native pixel-sector brightness defines a declared diagnostic mask,
            # but its EVPA samples are NOT the matched inference observable.
            intensity_1d = z['intensity'].mean(axis=0)
            polarized_1d = z['polarized_brightness'].mean(axis=0)
            angular_mask = intensity_1d >= .1 * intensity_1d.max()
            if mask_mode == 'intensity-polarization':
                angular_mask &= polarized_1d >= .2 * polarized_1d.max()
            angular_masks.append(angular_mask)
            if not np.array_equal(z['angle_deg'], angles):
                raise ValueError('Per-day aperture bins differ from audited axes')
        evpa_maps = np.load(paths[2], mmap_mode='r', allow_pickle=False)
        if evpa_maps.shape != maps.shape[:-1]:
            raise ValueError('Ensemble map EVPA shape mismatch')
        if (maps.ndim != 4 or maps.shape[1:] != (n, n, 4)
                or aperture.shape != (len(maps), 360) or len(maps) < 2):
            raise ValueError('Ensemble Stokes/aperture sample shape mismatch')
        if len(maps) != entry['sample_count']:
            raise ValueError('Declared ensemble sample count mismatch')
        mean = np.mean(maps, axis=0).transpose(*axes, 2)
        if not np.isfinite(mean).all():
            raise ValueError('Nonfinite ensemble mean Stokes')
        intensity, polarized = mean[..., 0], np.hypot(mean[..., 1], mean[..., 2])
        if intensity.max() <= 0 or polarized.max() <= 0:
            raise ValueError('Positive intensity and polarization required')
        mask = intensity >= .1 * intensity.max()
        if mask_mode == 'intensity-polarization':
            mask &= polarized >= .2 * polarized.max()
        mask &= np.isfinite(evpa_maps).all(axis=0).T
        masks.append(mask)
        data.append((evpa_maps, aperture))
    row, col = np.indices((n, n))
    radius = np.hypot(east, north)
    physical_mask = ((radius >= config['eht']['min_radius_rg'] * PhysicalScale(
        **config['physical']).angular_rg_uas)
        & (radius <= geometry['metadata']['radial_max_uas']))
    common = np.logical_and.reduce(masks) & physical_mask
    prepared = {}
    for dim in ('1d', '2d'):
        if dim == '1d':
            finite = np.logical_and.reduce([np.isfinite(v).all(axis=0) for _, v in data])
            selected = np.flatnonzero(products[dim][3] & finite
                & np.logical_and.reduce(angular_masks) & (np.arange(360) % 10 == 4))
            values = [np.asarray(aperture[:, bins[selected]], dtype=float) for _, aperture in data]
            coords = selected[:, None].astype(float)
        else:
            valid = products[dim][3] & common & (row % pixel_stride == 0) & (col % pixel_stride == 0)
            selected = np.flatnonzero(valid)
            r, c = np.unravel_index(selected, (n, n))
            native_r, native_c = (r, c) if axes == [0, 1] else (c, r)
            values = []
            for evpa_maps, _ in data:
                values.append(np.asarray(evpa_maps[:, native_r, native_c], dtype=float))
            coords = np.column_stack((east.ravel()[selected], north.ravel()[selected]))
        if not len(selected):
            raise ValueError(f'No initial {dim} features')
        prepared[dim] = (selected, values, coords, dict(
            mask_mode=mask_mode, intensity_fraction=.1, polarization_fraction=.2,
            min_radius_rg=config['eht']['min_radius_rg'], pixel_stride=pixel_stride,
            angular_spacing=10, angular_offset=4, coordinates='east/north uas; 1d template bin degrees',
            no_evpa_alignment=True, adapter=dict(mode=mode, resolution='display19',
                spatial='native[y,x] -> template[x,y]', template_axes='(-E,N)',
                angular='theory _integrate internal -90 once; output North0/East90',
                inference_1d='recomputed bilinear CSR theory quadrature on transposed per-draw Stokes',
                native_aperture_role='brightness selection diagnostics only; EVPA samples not used',
                quadrature=geometry['aperturefit'], qu_rotation='none')))
    processed = templates / 'raptor/background'
    processed_record = verify_stage(processed, recursive=True, allow_imported=False, _seen=verified)
    if processed_record['stage'] != 'raptor:background-preprocess':
        raise ValueError('Unexpected background preprocessing stage')
    preprocessing = {k: v for k, v in processed_record['config'].items() if k != 'state_indices'}
    if preprocessing != geometry['image_preprocess']:
        raise ValueError('Background preprocessing geometry mismatch')
    stokes_path = processed / 'stokes.npy'
    require_declared_outputs(processed, [stokes_path])
    sigmas, background_inputs = {}, [manifest_path(processed), stokes_path]
    for dim in ('1d', '2d'):
        directory = templates / 'background-unmasked' / dim
        background = verify_stage(directory, recursive=True, allow_imported=False, _seen=verified)
        if background['stage'] != 'raptor:background-unmasked':
            raise ValueError('Unexpected fixed-background source stage')
        summary = directory / 'summary.json'
        diagnostic = directory / 'evpa-difference.npz'
        require_declared_outputs(directory, [summary, diagnostic])
        metadata = json.loads(summary.read_text())
        if metadata['units'] != 'rad' or metadata['domain'] != dim:
            raise ValueError('Fixed-background domain/units mismatch')
        with np.load(diagnostic, allow_pickle=False) as z:
            sigma = float(z['sigma'])
        if sigma != metadata['sigma'] or not np.isfinite(sigma) or sigma < 0:
            raise ValueError('Invalid fixed-background sigma')
        sigmas[dim] = dict(stokes=str(stokes_path), cadence_hours=metadata['cadence_hours'],
                           lag=metadata['lag'], source_unmasked_sigma=sigma)
        background_inputs += [manifest_path(directory), summary, diagnostic]
    return prepared, inputs, dict(stage=record['stage'], centering=centering,
        background_selection='pooled scalar on final selected features, no response variance'), \
        sigmas, background_inputs


def run(args):
    config = read_config(args.config)
    if args.output_root.exists():
        raise FileExistsError('New output root required')
    with verification_session():
        seen = set()
        products, inputs = validate_templates(args.templates.resolve(), config, seen)
        prepared, source_inputs, source_identity, sigmas, background_inputs = prepare_limit_inputs(
            args.ensemble.resolve(), args.templates.resolve(), config, products,
            centering=args.centering, pixel_stride=args.pixel_stride,
            mask_mode=args.mask_mode, verified=seen)
        identity = dict(schema='eht120-observed-ensemble-limits-v1', scientific_profile=config,
                        resultant_min=args.resultant_min, centering=args.centering,
                        pixel_stride=args.pixel_stride, mask_mode=args.mask_mode,
                        ensemble=source_identity, covariance='D-draw population scatter; independent days',
                        scope='working-covariance conditional research credible limits')
        project = Path(__file__).resolve().parents[1]
        inputs += source_inputs + background_inputs + [args.config.resolve(), Path(__file__).resolve()]
        inputs += sorted((project / 'eht_axion').rglob('*.py'))
        inputs += [project / 'pixi.toml', project / 'pixi.lock']
        with stage_directory(args.output_root.resolve(), 'research:eht120-observed-limits',
                             identity, inputs) as out:
            if out is not None:
                return save_results(out, products, prepared, config, sigmas, identity)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('ensemble', 'templates', 'config', 'output-root'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--centering', choices=['native', 'centroid'], default='centroid')
    p.add_argument('--pixel-stride', type=int, default=5)
    p.add_argument('--mask-mode', choices=['intensity-polarization', 'intensity'],
                   default='intensity-polarization')
    p.add_argument('--resultant-min', type=float, default=.7)
    args = p.parse_args(argv)
    if args.pixel_stride < 1 or not np.isfinite(args.resultant_min) or not 0 < args.resultant_min <= 1:
        p.error('Invalid stride or resultant minimum')
    run(args)


if __name__ == '__main__':
    main()

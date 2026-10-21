"""Independent observed-image RESEARCH limits; original NULL pipeline is untouched.

Ensemble spread is a working covariance, NOT calibrated measurement noise.
No coupling-dependent response variance or GRMHD noise is introduced here.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.angles import circular_difference
from eht_axion.limits import LimitOptions
from eht_axion.observations import ObservationData
from eht_axion.pipeline.common import (
    manifest_path, stage_directory, verification_session, verify_stage, require_declared_outputs,
)
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.constants import PhysicalScale, analysis_geometry
from eht_axion.pipeline.products import load_templates
from eht_axion.research.observed import calculate_observed

DAYS = (95, 96, 100, 101)


def moments(values):
    """Population covariance on a recorded local pi-periodic branch."""
    center = .5 * np.angle(np.mean(np.exp(2j * values), axis=0))
    residual = circular_difference(values, center)
    mean = center + residual.mean(axis=0)
    residual -= residual.mean(axis=0)
    covariance = residual.T @ residual / len(values)
    return mean, covariance, np.abs(np.mean(np.exp(2j * values), axis=0))


def options_from_config(config):
    settings, cloud = config['inference'], config['cloud']
    return LimitOptions(
        phase_count=settings['phase_count'], coupling_points=settings['coupling_points'],
        coupling_min=settings['coupling_min'], coupling_max=settings['coupling_max'],
        prior=settings['prior'], time_conversion=1 / PhysicalScale(**config['physical']).tg_hours,
        exposure_hours=settings['exposure_hours'], day_gap_hours=settings['day_gap_hours'],
        pair_separation_hours=settings['pair_separation_hours'], n_principal=settings['principal_n'],
        reference_norm=cloud['reference_norm'], cloud_coefficient=cloud['coefficient'],
        cloud_pivot_alpha=cloud['pivot_alpha'], cloud_exponent=cloud['exponent'])


def run(args):
    config = read_config(args.config)
    study, templates = args.study.resolve(), args.templates.resolve()

    if args.output_root.exists():
        raise FileExistsError('New output root required')
    if args.variant == 'bestI':
        raise ValueError('A single fiducial image cannot estimate ensemble covariance')
    inputs = [args.config.resolve(), Path(__file__).resolve(), manifest_path(study), study / 'report.json']
    for dim in ('1d', '2d'):
        inputs += [manifest_path(templates / 'raptor/templates' / dim),
                   templates / 'raptor/templates' / dim / 'templates.npz']
    with verification_session():
        study_manifest = verify_stage(study, recursive=True, allow_imported=False)
        if study_manifest['stage'] != 'research:eht-polarization':
            raise ValueError('Unexpected study stage')
        consumed = [study / 'report.json']
        for day in DAYS:
            consumed += [study / f'{day}_{args.variant}_indices.npy',
                         study / f'{day}_{args.variant}_summary.npz',
                         study / f'{day}_candidates.npz', study / f'{day}_aligned_stokes.npy']
        require_declared_outputs(study, consumed)
        inputs += consumed
        report = json.loads((study / 'report.json').read_text())
        if report['partial']:
            raise ValueError('Cannot infer limits from a partial candidate study')
        geometry = analysis_geometry(config)
        cfg = study_manifest['config']
        if (not np.isclose(cfg['additional_blur_fwhm_uas'],
                           geometry['metadata']['gaussian_fwhm_uas'])
                or config['image']['final_shape'] != 64
                or not np.isclose(geometry['metadata']['final_fov_uas'], 128)
                or config['aperture']['weighting'] != 'intensity'
                or config['aperture']['window_degrees'] != 10
                or not np.isclose(geometry['metadata']['radial_max_uas'], 64)):
            raise ValueError('Observation and theory sky/aperture/kernel contracts differ')
        from eht_axion.producer.delivery import validate_processed_source

        seen = set()
        for alpha in config['data']['alphas']:
            from eht_axion.pipeline.layout import alpha_label

            label = alpha_label(alpha)
            processed = templates / 'raptor/processed' / label
            record = verify_stage(processed, recursive=True, allow_imported=False, _seen=seen)
            validate_processed_source(record, config, alpha=alpha, resolution=100, verified=seen)
            if record['config'] != geometry['image_preprocess']:
                raise ValueError('Theory preprocessing differs from selected research profile')
            fitted = verify_stage(templates / 'raptor/response' / label,
                                  recursive=True, allow_imported=False, _seen=seen)
            fit_config = {k: v for k, v in fitted['config'].items() if k != 'snapshot_labels'}
            if fit_config != geometry['aperturefit']:
                raise ValueError('Theory aperture differs from selected research profile')
        for dim in ('1d', '2d'):
            directory = templates / 'raptor/templates' / dim
            manifest = verify_stage(directory, recursive=True, allow_imported=False, _seen=seen)
            if manifest['stage'] != 'raptor:summarize':
                raise ValueError('Unexpected theory template stage')
            require_declared_outputs(directory, [directory / 'templates.npz'])
        identity = dict(schema='observed-image-research-v1', variant=args.variant, ridge=args.ridge,
                        resultant_min=args.resultant_min, pixel_stride=args.pixel_stride,
                        scientific_profile=config,
                        covariance='per-day hyperparameter ensemble; zero cross-day blocks assumed',
                        scope='working-covariance conditional credible limits; NOT calibrated exclusions')
        identity['alignment'] = getattr(args, 'alignment', 'aligned')
        identity['dimensions'] = getattr(args, 'dimensions', ['1d', '2d'])
        with stage_directory(args.output_root, 'research:observed-limits', identity, inputs) as out:
            results = {}
            for dim in identity['dimensions']:
                template_path = templates / 'raptor/templates' / dim / 'templates.npz'
                alpha, amp, phase, valid = load_templates(template_path, dim)
                if not np.array_equal(alpha, config['data']['alphas']):
                    raise ValueError('Theory alpha grid differs from selected research profile')
                selected, samples = prepare_samples(study, args, dim, valid)
                daily = [moments(v) for v in samples]
                # Selection uses ensemble concentration, not day-difference size.
                good = np.logical_and.reduce([v[2] >= args.resultant_min for v in daily])
                selected = selected[good]
                if not selected.size:
                    raise ValueError(f'No concentrated {dim} features remain')
                means = tuple(v[0][good] for v in daily)
                covs = []
                for _, cov, _ in daily:
                    cov = cov[np.ix_(good, good)]
                    cov = (1 - args.ridge) * cov + args.ridge * np.diag(np.diag(cov))
                    # No automatic eigenvalue clipping or hidden regularization.
                    np.linalg.cholesky(cov)
                    covs.append(cov)
                differences = tuple(circular_difference(means[b], means[a]) for a, b in ((0, 1), (2, 3)))
                obs = ObservationData(means, tuple(np.sqrt(np.diag(c)) for c in covs), tuple(covs),
                                      differences, (covs[0]+covs[1], covs[2]+covs[3]), selected)
                null = replace(obs, pair_differences=tuple(np.zeros_like(v) for v in differences))
                np.savez(out / f'observations-{dim}.npz', selected_indices=selected, means=means,
                         covariances=covs, pair_differences=differences,
                         resultant=np.stack([v[2][good] for v in daily]),
                         sample_counts=[len(v) for v in samples])
                for mode in ('std', 'cor'):
                    for label, data in (('observed', obs), ('null-comparison', null)):
                        a, cuts, diagnostics = calculate_observed(
                            data, amp, phase, alphas=alpha, dimension=dim, mode=mode,
                            options=options_from_config(config))
                        key = f'{label}-{mode}-{dim}'
                        np.savez(out / f'{key}.npz', alphas=a, log10_c95=cuts, c95=10**cuts)
                        (out / f'{key}.json').write_text(
                            json.dumps(diagnostics, indent=2, allow_nan=False)+'\n')
                        results[key] = {'alphas': a.tolist(), 'c95': (10**cuts).tolist(),
                                        'features': len(selected), 'log10_c95': cuts.tolist()}
                        print(key, cuts, flush=True)
            plot_results(out, results, args.variant)
            (out / 'report.json').write_text(json.dumps(dict(identity, results=results,
                warnings=['Hyperparameter spread is not a calibrated sampling covariance.',
                          'Independent day errors assumed; not inferred from matched identities.',
                          'EVPA local Gaussian branch, finite prior and phase/coupling grids.',
                          'Translation reference is local Apr11 best-I image, not official consensus.',
                          'No GRMHD response variance; observed/null share the same covariance.']),
                indent=2, allow_nan=False)+'\n')


def plot_results(out, results, variant):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    for ax, dim in zip(axes, ('1d', '2d')):
        if f'observed-std-{dim}' not in results:
            ax.set_visible(False)
            continue
        for mode, color in (('std', '#2679b5'), ('cor', '#e58f24')):
            for label, style in (('observed', '-'), ('null-comparison', '--')):
                row = results[f'{label}-{mode}-{dim}']
                ax.plot(row['alphas'], row['log10_c95'], style, color=color,
                        label=f'{label} {mode}', marker='o', markersize=3)
        ax.set(xlabel='Axion alpha', ylabel='log10 c95', title=f'{variant}: {dim}')
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle('RESEARCH ONLY: reconstructed day differences / working ensemble covariance')
    fig.savefig(out / 'observed-vs-null.png', dpi=160)
    plt.close(fig)


def prepare_samples(study, args, dim, template_valid):
    if dim == '1d':
        selected = np.flatnonzero(template_valid & (np.arange(360) % 10 == 4))
    else:
        masks = []
        for day in DAYS:
            with np.load(study / f'{day}_{args.variant}_summary.npz') as z:
                image = z['mean_stokes'][1]
            i, p = image[..., 0], np.hypot(image[..., 1], image[..., 2])
            masks.append(((i >= .1*i.max()) & (p >= .2*p.max())).T)
        y, x = np.indices(template_valid.shape)
        selected = np.flatnonzero(template_valid & np.logical_and.reduce(masks)
                                  & (y % args.pixel_stride == 0) & (x % args.pixel_stride == 0))
    if not len(selected):
        raise ValueError(f'No initial {dim} features')
    samples = []
    for day in DAYS:
        indices = np.load(study / f'{day}_{args.variant}_indices.npy')
        if len(indices) < 2:
            raise ValueError('At least two images required for covariance')
        if dim == '1d':
            with np.load(study / f'{day}_candidates.npz') as z:
                mode = 0 if getattr(args, 'alignment', 'aligned') == 'raw' else 1
                values = z['chi'][mode, indices][:, selected]
        else:
            maps = np.load(study / f'{day}_aligned_stokes.npy', mmap_mode='r')
            # Native(y,north;x,east-decreasing) -> existing template(axis0,axis1).
            # Constant Q/U sign reversal rotates both days90deg; differences invariant.
            rows, cols = np.unravel_index(selected, template_valid.shape)
            q = np.asarray(maps[indices[:, None], cols[None, :], rows[None, :], 1], dtype=float)
            u = np.asarray(maps[indices[:, None], cols[None, :], rows[None, :], 2], dtype=float)
            values = .5*np.arctan2(u, q)
            if np.any((q == 0) & (u == 0)):
                raise ValueError('Zero-polarization samples cannot define EVPA')
            del maps
        if not np.isfinite(values).all():
            raise ValueError('Missing/undefined EVPA samples: no silent finite-only covariance')
        samples.append(values)
    return selected, samples


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--templates', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--variant', choices=['all4800', 'legacy_evpa', 'datafit_perday', 'datafit_shared'],
                   required=True)
    p.add_argument('--alignment', choices=['raw', 'aligned'], default='aligned')
    p.add_argument('--dimensions', nargs='+', choices=['1d', '2d'], default=['1d', '2d'])
    p.add_argument('--ridge', type=float, default=0, help='Explicit covariance shrinkage fraction')
    p.add_argument('--resultant-min', type=float, default=.7, help='Declared branch-concentration selection')
    p.add_argument('--pixel-stride', type=int, default=5,
                   help='2D grid stride in pixels; research sampling choice')
    args = p.parse_args(argv)
    if not 0 <= args.ridge < 1 or not 0 < args.resultant_min <= 1 or args.pixel_stride < 1:
        p.error('Invalid ridge, concentration or stride')
    if args.alignment == 'raw' and '2d' in args.dimensions:
        p.error('Unregistered sensitivity currently supports only the saved 1D curves')
    run(args)


if __name__ == '__main__':
    main()

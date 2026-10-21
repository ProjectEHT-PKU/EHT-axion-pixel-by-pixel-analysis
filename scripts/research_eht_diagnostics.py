"""Plot measured EHT aperture distributions and selection/alignment sensitivity.

This is a single-pipeline hyperparameter study, not official D-term samples.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from eht_axion.angles import circular_difference
from eht_axion.pipeline.common import (
    manifest_path, require_declared_outputs, stage_directory, verification_session, verify_stage,
)

DAYS = (95, 96, 100, 101)
VARIANTS = ('all4800', 'legacy_evpa', 'datafit_perday', 'datafit_shared')


def quantiles(values):
    a = np.asarray(values)
    a = a[np.isfinite(a)]
    return None if not len(a) else np.quantile(a, [.05, .5, .95]).tolist()


def run(args):
    study = args.study.resolve()
    if args.output_root.exists():
        raise FileExistsError('New research diagnostic output required')
    with verification_session():
        verify_stage(study, recursive=True, allow_imported=False)
        consumed = [study / 'report.json']
        for day in DAYS:
            consumed += [study / f'{day}_candidates.npz']
            consumed += [study / f'{day}_{v}_indices.npy' for v in VARIANTS]
            consumed += [study / f'{day}_{v}_summary.npz' for v in VARIANTS]
        require_declared_outputs(study, consumed)
        report = json.loads((study / 'report.json').read_text())
        if report['partial']:
            raise ValueError('Full-ensemble diagnostics require completed full study')
        identity = dict(schema='eht-polarization-diagnostics-v1',
                        concentration_threshold=args.resultant_min,
                        scope='working ensemble distributions; not repeated observations')
        with stage_directory(args.output_root, 'research:eht-diagnostics', identity,
                             [Path(__file__).resolve(), manifest_path(study)] + consumed) as out:
            results = {}
            for variant in VARIANTS:
                arrays, means, summaries = [], [], []
                for day in DAYS:
                    indices = np.load(study / f'{day}_{variant}_indices.npy')
                    with np.load(study / f'{day}_candidates.npz') as z:
                        arrays.append(z['chi'][:, indices])
                    with np.load(study / f'{day}_{variant}_summary.npz') as z:
                        means.append(z['mean_stokes'].copy())
                        summaries.append({k: z[k].copy() for k in z.files})
                fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
                angular = np.arange(360)
                variant_result = {}
                for mode, alignment in enumerate(('raw', 'aligned')):
                    axis = axes[0, mode]
                    for k, day in enumerate(DAYS):
                        sample = arrays[k][mode]
                        center = .5*np.angle(np.mean(np.exp(2j*sample), axis=0))
                        residual = circular_difference(sample, center)
                        q = np.quantile(residual, [.16, .84], axis=0)
                        axis.plot(angular, np.rad2deg(center), label=f'Apr{day-90}')
                        axis.fill_between(angular, np.rad2deg(center+q[0]),
                                          np.rad2deg(center+q[1]), alpha=.1)
                        image = means[k][mode]
                        intensity, p = image[..., 0], np.hypot(image[..., 1], image[..., 2])
                        east, north = summaries[k]['east_uas'], summaries[k]['north_uas']
                        x, y = np.meshgrid(east, north)
                        skyphi = np.arctan2(x, y)
                        beta = np.sum((image[..., 1]+1j*image[..., 2])*np.exp(-2j*skyphi))/intensity.sum()
                        variant_result[f'{alignment}-{day}'] = dict(
                            net_polarization_fraction=float(abs(np.sum(image[..., 1]+1j*image[..., 2]))
                                                            /intensity.sum()),
                            resolved_polarization_fraction=float(p.sum()/intensity.sum()),
                            beta2_abs=float(abs(beta)), beta2_phase_degrees=float(np.rad2deg(np.angle(beta))),
                            aperture_resultant_quantiles=quantiles(abs(np.mean(np.exp(2j*sample), axis=0))))
                    axis.set(title=f'{variant}: {alignment}, ensemble directions',
                             xlabel='Sky PA (deg, E of N)', ylabel='EVPA (deg, modulo180)', ylim=(-90, 90))
                    axis.legend()
                    axis.grid(alpha=.2)
                    axis = axes[1, mode]
                    for a, b in ((0, 1), (2, 3)):
                        za = np.mean(np.exp(2j*arrays[a][mode]), axis=0)
                        zb = np.mean(np.exp(2j*arrays[b][mode]), axis=0)
                        delta = .5*np.angle(zb*za.conj())
                        good = (abs(za) >= args.resultant_min) & (abs(zb) >= args.resultant_min)
                        axis.plot(angular, np.rad2deg(delta), alpha=.25)
                        axis.plot(angular, np.where(good, np.rad2deg(delta), np.nan),
                                  label=f'Apr{DAYS[b]-90} - Apr{DAYS[a]-90}')
                        variant_result[f'{alignment}-{DAYS[a]}-{DAYS[b]}'] = dict(
                            concentrated_bins=int(good.sum()),
                            delta_degrees_quantiles=quantiles(np.rad2deg(delta[good])),
                            absolute_delta_degrees_quantiles=quantiles(abs(np.rad2deg(delta[good]))))
                    axis.set(xlabel='Sky PA (deg, E of N)', ylabel='Later minus earlier EVPA (deg)',
                             ylim=(-90, 90))
                    axis.axhline(0, color='k', lw=.6)
                    axis.grid(alpha=.2)
                    axis.legend()
                fig.suptitle('RESEARCH: hyperparameter distributions, NOT official D-term ensemble')
                fig.savefig(out / f'{variant}-aperture-diagnostics.png', dpi=150)
                plt.close(fig)
                # Fig8-like distribution plot uses counts, never claims1000officialsamples.
                fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
                for ax, k in zip(axes, (0, 3)):
                    sample = np.rad2deg(arrays[k][1])
                    edges = np.linspace(-90, 90, 181)
                    hist = np.array([np.histogram(sample[:, j], bins=edges)[0] for j in range(360)]).T
                    mesh = ax.pcolormesh(np.arange(361)-.5, edges, hist, shading='flat', cmap='magma')
                    ax.set(title=f'Apr{DAYS[k]-90}: {len(sample)} hyperparameter images',
                           xlabel='Sky PA (deg, E of N)', ylabel='EVPA (deg)')
                    fig.colorbar(mesh, ax=ax, label='Candidate count (not probability)')
                fig.savefig(out / f'{variant}-fig8-like-counts.png', dpi=150)
                plt.close(fig)
                results[variant] = variant_result
            (out / 'report.json').write_text(json.dumps(dict(identity, results=results,
                notes=['Angular bands cross the modulo180 seam; not linear confidence bands.',
                       'Registration reference is local best-I Apr11, not official consensus.',
                       'Native restored5uas plus additional19uas kernel; no known effective PSF.',
                       'Do not use Fig6 display masks as calibrated significance selections.']),
                indent=2, allow_nan=False)+'\n')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--resultant-min', type=float, default=.7)
    args = p.parse_args(argv)
    if not 0 < args.resultant_min <= 1:
        p.error('resultant-min must be positive and <=1')
    run(args)


if __name__ == '__main__':
    main()

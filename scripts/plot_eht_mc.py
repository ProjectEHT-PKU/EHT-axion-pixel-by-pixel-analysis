#!/usr/bin/env python3
"""Plot a completed conditional fixed-D MC summary, never partial producer data.

python scripts/plot_eht_mc.py --ensemble ROOT --output-root NEW_ROOT
No official FITS, thermal-noise bootstrap, or observed-limit products are inferred.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

from eht_axion.pipeline import common  # noqa: E402
from eht_axion.research.ensemble import DAYS, PAIRS, branch_difference  # noqa: E402

MODES = ('icentroid', 'native')
RESOLUTION = 'display19'
SCOPE = 'Conditional fixed-I D-calibration MC only; not formal all-noise uncertainty or official FITS'


def tick_segments(east, north, chi, length):
    """Sky-coordinate endpoints: EVPA zero is North, +90 degrees is East."""
    dx, dy = .5 * length * np.sin(chi), .5 * length * np.cos(chi)
    return np.stack((np.stack((east-dx, north-dy), axis=-1),
                     np.stack((east+dx, north+dy), axis=-1)), axis=-2)


def circular_band(samples, mean, valid):
    """16/84 local residual percentiles; break display branches at pi seams.

    The band is not the percentile of wrapped angles. Split also when either
    boundary jumps, so fill_between never bridges an artificial EVPA seam.
    """
    mean = np.asarray(mean)
    lo, hi = np.percentile(branch_difference(samples, mean), [16, 84], axis=0)
    center = branch_difference(mean, 0.)
    lower, upper = center + lo, center + hi
    keep = np.asarray(valid, dtype=bool) & np.isfinite(center + lower + upper)
    for curve in (center, lower, upper):
        jump = np.flatnonzero(np.abs(np.diff(curve)) > np.pi/2)
        keep[jump + 1] = False
    return tuple(np.where(keep, np.rad2deg(a), np.nan) for a in (center, lower, upper))


def display_mask(stokes, valid=None):
    i, q, u = (stokes[..., k] for k in range(3))
    p = np.hypot(q, u)
    mask = np.isfinite(i+p) & (i > .1*np.max(i)) & (p > .2*np.max(p))
    if valid is not None:
        mask &= valid
    return mask


def load_npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def save(fig, out, name):
    fig.savefig(out / (name + '.png'), dpi=170, bbox_inches='tight')
    plt.close(fig)


def plot_maps(source, out, axes, mode, records):
    e, n = axes['east_uas'], axes['north_uas']
    # Pixel-edge extent preserves the FITS pixel centers, including reversed E.
    extent = [e[0]-(e[1]-e[0])/2, e[-1]+(e[1]-e[0])/2,
              n[0]-(n[1]-n[0])/2, n[-1]+(n[1]-n[0])/2]
    ee, nn = np.meshgrid(e, n)
    fig, axs = plt.subplots(1, 4, figsize=(14, 4), sharex=True, sharey=True)
    lc = None
    for ax, day in zip(axs, DAYS):
        data = load_npz(source / mode / RESOLUTION / ('day-'+day) / 'map-summary.npz')
        s = data['mean_stokes']
        i, q, u = (s[..., k] for k in range(3))
        p = np.hypot(q, u)
        chi = .5*np.arctan2(u, q)  # deliberately NOT circular mean(per-draw EVPA)
        mask = display_mask(s)
        sparse = np.zeros(mask.shape, bool)
        sparse[::3, ::3] = True
        selected = mask & sparse
        ax.imshow(i / i.max(), origin='lower', extent=extent, cmap='gray', vmin=0, vmax=1)
        segments = tick_segments(ee[selected], nn[selected], chi[selected], 9*p[selected]/p.max())
        lc = LineCollection(segments, cmap='turbo', norm=plt.Normalize(0, .3), linewidths=1.4)
        lc.set_array(p[selected]/i[selected])
        ax.add_collection(lc)
        ax.set(title=f'April {int(day)-90}, 2017 (day {day})', xlabel='East offset (µas)', aspect='equal')
        records[f'{mode}/day-{day}/ticks'] = dict(eligible_pixels=int(mask.sum()),
                                                 drawn_ticks=int(selected.sum()))
    axs[0].set_ylabel('North offset (µas)')
    fig.colorbar(lc, ax=list(axs), label='|mean Q+i mean U| / mean I (color clipped at 0.30)', shrink=.75)
    fig.suptitle(f'{mode}, additional 19 µas display beam; mean I / daily peak\n'
                 'Ticks: EVPA of mean Stokes; length = 9 µas × P / daily peak P; I>0.1 peak, P>0.2 peak')
    fig.text(.5, -.03, SCOPE, ha='center', fontsize=9)
    save(fig, out, 'mean-stokes-'+mode)


def plot_apertures(source, out):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for row, mode in enumerate(MODES):
        for col, days in enumerate((DAYS[:2], DAYS[2:])):
            ax = axs[row, col]
            for day, color in zip(days, ('tab:blue', 'tab:orange')):
                d = load_npz(source / mode / RESOLUTION / ('day-'+day) / 'aperture.npz')
                c, lo, hi = circular_band(d['samples_rad'], d['mean_rad'], d['valid'])
                x = d['angle_deg']
                ax.plot(x, c, color=color, label=f'day {day}')
                ax.fill_between(x, lo, hi, color=color, alpha=.25)
            ax.set(title=mode, ylabel='EVPA (deg, local branch)', xlim=(0, 359))
            ax.legend()
            ax.grid(alpha=.2)
    for ax in axs[-1]:
        ax.set_xlabel('Aperture azimuth (deg E of N; North=0, East=90)')
    fig.suptitle('Circular mean of per-draw aperture EVPA + 16/84% local-residual interval\n'
                 '10° full-width windows, 1° steps, radius ≤60 µas; additional 19 µas display beam')
    fig.text(.5, .005, SCOPE, ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, .92))
    save(fig, out, 'aperture-evpa')


def plot_globals(source, out, records):
    fig, axs = plt.subplots(2, 5, figsize=(17, 7))
    labels = ['Net m', 'Resolved m', '|β₂|', 'arg β₂ (deg, local branch)',
              'Integrated EVPA (deg, local branch)']
    for row, mode in enumerate(MODES):
        for day in DAYS:
            d = load_npz(source / mode / RESOLUTION / ('day-'+day) / 'global.npz')
            columns = list(d['columns'])
            values = {k: d['values'][:, columns.index(k)] for k in columns}
            beta = values['beta2_real'] + 1j*values['beta2_imag']
            # Center circular histograms, not arbitrary wrapped distributions.
            phase = np.angle(beta)
            pc = np.angle(np.mean(np.exp(1j*phase)))
            chi = values['evpa_rad']
            cc = .5*np.angle(np.mean(np.exp(2j*chi)))
            arrays = [values['m_net'], values['m_resolved'], abs(beta),
                      np.rad2deg(pc + np.angle(np.exp(1j*(phase-pc)))),
                      np.rad2deg(cc + branch_difference(chi, cc))]
            records[f'{mode}/day-{day}/global'] = {
                label: dict(percentile_16_50_84=np.percentile(a, [16, 50, 84]).tolist())
                for label, a in zip(labels, arrays)}
            for ax, a, label in zip(axs[row], arrays, labels):
                ax.hist(a, bins=32, histtype='step', density=True, label=day)
                ax.set_xlabel(label)
        axs[row, 0].set_ylabel(mode+'\nProbability density')
        axs[row, 0].legend(title='Day')
    fig.suptitle('Whole-field per-draw diagnostics; β₂ uses sky azimuth E of N, no annular cut\n'
                 'additional 19 µas display beam; histogram angles locally centered per day')
    fig.text(.5, .005, SCOPE, ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, .91))
    save(fig, out, 'global-diagnostics')


def plot_differences(source, out, axes, records):
    fig, axs = plt.subplots(2, 4, figsize=(17, 8))
    e, n = axes['east_uas'], axes['north_uas']
    extent = [e[0]-(e[1]-e[0])/2, e[-1]+(e[1]-e[0])/2,
              n[0]-(n[1]-n[0])/2, n[-1]+(n[1]-n[0])/2]
    im = None
    for row, mode in enumerate(MODES):
        for col, (later, earlier) in enumerate(PAIRS):
            base = source / mode / RESOLUTION
            pair = load_npz(base / 'pairs' / f'{later}-{earlier}.npz')
            diff = np.rad2deg(pair['difference_rad'])
            sigma = np.rad2deg(np.sqrt(np.maximum(0, np.diag(pair['covariance_rad2']))))
            valid = pair['valid']
            ax = axs[row, col]
            # Neither the line nor band may bridge a branch discontinuity.
            keep = valid.copy()
            keep[1:] &= np.abs(np.diff(diff)) < 90
            ax.plot(axes['angle_deg'], np.where(keep, diff, np.nan))
            ax.fill_between(axes['angle_deg'], np.where(keep, diff-sigma, np.nan),
                            np.where(keep, diff+sigma, np.nan), alpha=.25)
            ax.set(title=f'{mode}: {later} − {earlier}', xlabel='Azimuth (deg E of N)',
                   ylabel='Δ EVPA (deg)', xlim=(0, 359))
            ax.axhline(0, color='gray', lw=.6)
            a = load_npz(base / ('day-'+later) / 'map-summary.npz')
            b = load_npz(base / ('day-'+earlier) / 'map-summary.npz')
            mask = display_mask(a['mean_stokes'], a['valid']) & display_mask(b['mean_stokes'], b['valid'])
            delta = np.rad2deg(branch_difference(a['mean_rad'], b['mean_rad']))
            ax = axs[row, col+2]
            im = ax.imshow(np.where(mask, delta, np.nan), origin='lower', extent=extent,
                           cmap='RdBu_r', vmin=-90, vmax=90)
            ax.set(title=f'{mode}: {later} − {earlier}',
                   xlabel='East offset (µas)', ylabel='North offset (µas)')
            records[f'{mode}/{later}-{earlier}'] = dict(map_display_pixels=int(mask.sum()),
                                                       aperture_valid_bins=int(valid.sum()))
    fig.colorbar(im, ax=list(axs[:, 2:].ravel()), label='Δ circular-mean EVPA (deg)', shrink=.7)
    fig.suptitle('Independent-day differences of circular means; additional 19 µas display beam\n'
                 'Left bands: Gaussian ±1σ, σ²=diag(C_later+C_earlier), NOT paired-draw percentiles; '
                 'right: brightness + resultant mask intersection')
    fig.text(.5, .005, SCOPE, ha='center', fontsize=9)
    fig.subplots_adjust(wspace=.4, hspace=.4, top=.86, bottom=.1)
    save(fig, out, 'day-differences')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ensemble', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args(argv)
    source, destination = args.ensemble.resolve(), args.output_root.resolve()
    if source.is_relative_to(destination) or destination.is_relative_to(source):
        parser.error('ensemble and output roots must be separate and non-nested')
    with common.verification_session():
        # Only summary outputs are consumed: deliberately no 97 GB FITS ancestry walk.
        common.verify_stage(source, recursive=False)
        report = json.loads((source / 'report.json').read_text())
        if (report.get('schema') != 'eht-fixed-D-ensemble-summary-v1'
                or report.get('days') != list(DAYS) or report.get('draws_per_day') != 1000):
            raise ValueError('expected complete four-day, 1000 draws/day summary')
        products = [source / 'axes.npz', source / 'report.json']
        for mode in MODES:
            for day in DAYS:
                products += [source / mode / RESOLUTION / ('day-'+day) / name
                             for name in ('map-summary.npz', 'aperture.npz', 'global.npz')]
            products += [source / mode / RESOLUTION / 'pairs' / f'{a}-{b}.npz' for a, b in PAIRS]
        common.require_declared_outputs(source, products)
        config = dict(source_manifest_sha256=common.sha256(common.manifest_path(source)),
                      source_verification='all summary output SHA256; no raw FITS ancestry rehash',
                      modes=list(MODES), resolution=RESOLUTION, scope=SCOPE,
                      display_I_peak_threshold=.1, display_P_peak_threshold=.2)
        inputs = products + [Path(__file__).resolve(), ROOT / 'eht_axion/research/ensemble.py',
                             ROOT / 'eht_axion/pipeline/common.py']
        with common.stage_directory(destination, 'research:eht-fixed-D-mc-figures', config, inputs) as out:
            if out is not None:
                axes = load_npz(source / 'axes.npz')
                records = {}
                for mode in MODES:
                    plot_maps(source, out, axes, mode, records)
                plot_apertures(source, out)
                plot_globals(source, out, records)
                plot_differences(source, out, axes, records)
                summary = dict(schema='eht-fixed-D-mc-figures-v1', source=str(source), **config,
                               frequency_hz=report['frequency_hz'], nominal_RT_frequency_hz=230e9,
                               frequency_note='Actual FITS frequency, not an exact 230 GHz match',
                               bunit=report['bunit'], registrations=report['registrations'],
                               frame=report['registration'], angle_convention=report['angle_convention'],
                               aperture_rmin=report['rmin'],
                               map_estimator='ticks: EVPA(mean IQU); differences: circular mean(draw EVPA)',
                               mask_note='Display-only I>0.1 daily peak and P>0.2 daily peak; '
                                         'difference maps also intersect existing per-day valid masks. '
                                         'No analysis selection changed.',
                               uncertainty='EVPA rad; covariance rad²; plots deg and Gaussian sigma deg. '
                                           'β2 phase is 2pi periodic, EVPA pi periodic.',
                               records=records)
                (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

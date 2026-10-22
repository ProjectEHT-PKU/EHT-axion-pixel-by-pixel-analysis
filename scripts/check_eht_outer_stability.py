#!/usr/bin/env python3
"""Compare saved outer iterates in physical maps; no optimization or registration.

Round numbers refer to ADDITIONAL rounds after the original 50-round fiducial.
Saved round-N P precedes that round's D update; final P is aligned to final D.
Diagnostics do not turn a numeric budget or a fixed-D gap into statistical errors.
"""
import argparse
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from eht_axion.pipeline import common
from eht_axion.research.polarimetry import load_stokes_fits, aperture_curves
from diagnose_eht_apr11 import map_metrics


def compare(before, after, mask):
    a, b = before[..., 1]+1j*before[..., 2], after[..., 1]+1j*after[..., 2]
    result = map_metrics(before[..., 0], a, after[..., 0], b, mask)
    result['integrated_evpa_difference_deg'] = float(np.rad2deg(.5*np.angle(b.sum()*a.sum().conj())))
    return result


def summary(image, phi):
    intensity = image[..., 0]
    p = image[..., 1]+1j*image[..., 2]
    beta = np.sum(p*np.exp(-2j*phi))/intensity.sum()
    return dict(net_m=float(abs(p.sum())/intensity.sum()),
                resolved_m=float(abs(p).sum()/intensity.sum()),
                integrated_evpa_deg=float(np.rad2deg(.5*np.angle(p.sum()))),
                beta2_abs=float(abs(beta)), beta2_phase_deg=float(np.rad2deg(np.angle(beta))))


def main(argv=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    args = parser.parse_args(argv)
    labels = ['additional50', 'additional80', 'additional100', 'final-P']
    stems = ['round-50', 'round-80', 'round-100', 'final']
    with common.verification_session():
        record = common.verify_stage(args.source_root, recursive=True, allow_imported=False)
        if record['stage'] != 'research:eht-convex:fiducial':
            raise ValueError('Expected certified-P outer calibration stage')
        files = [args.source_root/(s+'-'+suffix+'.fits') for s in stems
                 for suffix in ['raw-IQU', 'display-IQU-additional19uas']]
        common.require_declared_outputs(args.source_root, files)
        inputs = list(common.source_closure([Path(__file__).resolve(), common.manifest_path(args.source_root),
            Path(__file__).with_name('diagnose_eht_apr11.py'),
            Path(__file__).resolve().parents[1]/'eht_axion/research/polarimetry.py', *files]))
        if args.output_root.exists() or any(args.output_root.absolute().is_relative_to(p) for p in inputs):
            raise ValueError('new separate output root required')
        loaded = [load_stokes_fits(p) for p in files]
        ref = loaded[0]
        for item in loaded:
            np.testing.assert_allclose(item['east_uas'], ref['east_uas'], rtol=0, atol=1e-9)
            np.testing.assert_allclose(item['north_uas'], ref['north_uas'], rtol=0, atol=1e-9)
        raw = [item['stokes'] for item in loaded[::2]]
        display = [item['stokes'] for item in loaded[1::2]]
        # Identical I within each raw/display group; no hidden position adjustment.
        for group in (raw, display):
            for image in group[1:]:
                np.testing.assert_array_equal(image[..., 0], group[0][..., 0])
        mask = np.logical_and.reduce([(im[..., 0] >= .1*im[..., 0].max())
            & (np.hypot(im[..., 1], im[..., 2]) >= .2*np.hypot(im[..., 1], im[..., 2]).max())
            for im in display])
        east, north = ref['east_uas'], ref['north_uas']
        xx, yy = np.meshgrid(east, north)
        centroid = [float(np.sum(display[0][..., 0]*v)/display[0][..., 0].sum()) for v in (xx, yy)]
        phi = np.arctan2(xx-centroid[0], yy-centroid[1])
        angles = np.arange(360.)
        curves = [aperture_curves(im, east-centroid[0], north-centroid[1], angles, window=10)
                  for im in display]
        data = dict(day=record['config']['day'], source=str(args.source_root),
            round_interpretation='additional rounds after 50 original; round P is before its D update',
            display='same additional 19uas IQU kernel; NOT effective PSF',
            mask='four-map intersection: I>=.1peak and P>=.2peak; diagnostic not official',
            mask_pixels=int(mask.sum()), centroid_uas=centroid,
            beta2='whole-map I-centroid diagnostic; not certified official centering procedure',
            summaries={label: summary(im, phi) for label, im in zip(labels, display)}, comparisons={})
        for i, j in [(0, 1), (1, 2), (2, 3), (0, 3), (1, 3)]:
            c = compare(display[i], display[j], mask)
            c['raw_qu_relative_l2'] = float(np.linalg.norm(raw[j][..., 1:3]-raw[i][..., 1:3])
                                           /np.linalg.norm(raw[i][..., 1:3]))
            delta = .5*np.angle(np.exp(2j*(curves[j]['evpa_rad']-curves[i]['evpa_rad'])))
            bright = np.minimum(curves[i]['polarized_brightness'], curves[j]['polarized_brightness'])
            good = np.isfinite(delta) & (bright >= .2*bright.max())
            c['aperture_evpa_rms_deg'] = float(np.rad2deg(np.sqrt(np.mean(delta[good]**2))))
            c['aperture_bins'] = int(good.sum())
            data['comparisons'][labels[i]+'-to-'+labels[j]] = c
        config = dict(source=str(args.source_root), labels=labels, mask_I=.1, mask_P=.2,
                      apertures='10degree full width, 1degree steps; intensity-weighted coherent P',
                      no_registration=True, scientific_acceptance_threshold=None)
        with common.stage_directory(args.output_root, 'research:eht-outer-stability', config, inputs) as out:
            (out/'report.json').write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
            np.save(out/'common-mask.npy', mask)
            with (out/'apertures.csv').open('w') as stream:
                writer = csv.writer(stream)
                writer.writerow(['label', 'angle_deg', 'evpa_deg', 'polarized_brightness_Jy_pixel'])
                for label, c in zip(labels, curves):
                    writer.writerows(zip([label]*360, angles, np.rad2deg(c['evpa_rad']),
                                         c['polarized_brightness']))
            fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
            for label, c in zip(labels, curves):
                # Unwrap relative to final for display only; CSV preserves axial angles.
                delta = .5*np.angle(np.exp(2j*(c['evpa_rad']-curves[-1]['evpa_rad'])))
                axes[0, 0].plot(angles, np.rad2deg(delta), label=label)
                axes[0, 1].plot(angles, c['polarized_brightness'], label=label)
            axes[0, 0].set(title='Aperture EVPA difference to final', ylabel='degrees', xlabel='PA E of N')
            axes[0, 1].set(title='Intensity-weighted polarized brightness', xlabel='PA E of N')
            axes[0, 0].legend()
            vals = [data['summaries'][label] for label in labels]
            axes[1, 0].plot(labels, [100*v['net_m'] for v in vals], 'o-', label='net m')
            axes[1, 0].plot(labels, [100*v['resolved_m'] for v in vals], 'o-', label='resolved m')
            axes[1, 0].set(ylabel='percent', title='Same-display-kernel polarization fractions')
            axes[1, 0].legend()
            initial = display[0][..., 1]+1j*display[0][..., 2]
            final = display[-1][..., 1]+1j*display[-1][..., 2]
            delta = np.rad2deg(.5*np.angle(final*initial.conj()))
            panel = axes[1, 1].imshow(np.where(mask, delta, np.nan), origin='lower',
                extent=[east[0], east[-1], north[0], north[-1]], cmap='RdBu_r')
            axes[1, 1].set(title='Additional50 to final: EVPA change', xlabel='East offset (uas)',
                           ylabel='North offset (uas)')
            fig.colorbar(panel, ax=axes[1, 1], label='degrees')
            fig.savefig(out/'stability.png', dpi=160)
            plt.close(fig)
    return data


if __name__ == '__main__':
    main()

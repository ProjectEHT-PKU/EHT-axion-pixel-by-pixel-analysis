"""Render a completed isolated ehtim rebuild, preserving its execution provenance."""
from pathlib import Path
import argparse
import json

import _bootstrap  # noqa: F401
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eht_axion.pipeline.common import (
    verification_session, verify_stage, require_declared_outputs, stage_directory, manifest_path,
)
from eht_axion.research.polarimetry import load_stokes_fits


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args(argv)
    source = args.input_root.resolve()
    paths = [source / n for n in ('final-raw-IQU.fits', 'final-display-IQU-additional19uas.fits',
                                  'steps.json', 'termination.json', 'observation.json', 'plan.json')]
    with verification_session():
        verify_stage(source, recursive=True, allow_imported=False)
        require_declared_outputs(source, paths)
        termination = json.loads((source / 'termination.json').read_text())
        steps = json.loads((source / 'steps.json').read_text())
        with stage_directory(args.output_root, 'research:eht-fiducial-render',
                             {'termination': termination}, [Path(__file__).resolve(), manifest_path(source)]
                             + paths) as out:
            fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout='constrained')
            titles = ['Raw reconstruction', 'Additional19uas I/Q/U kernel']
            for ax, path, title in zip(axes, paths[:2], titles):
                loaded = load_stokes_fits(path)
                image, east, north = (loaded[k] for k in ('stokes', 'east_uas', 'north_uas'))
                area = np.prod(loaded['pixel_uas'])
                i, q, u = [image[..., j] for j in range(3)]
                p = np.hypot(q, u)
                angle = .5*np.arctan2(u, q)
                x, y = np.meshgrid(east, north)
                mesh = ax.pcolormesh(x, y, i/area*1e6, cmap='gray_r', shading='nearest')
                mask = (i >= .1*i.max()) & (p >= .2*p.max()) & (i > 0)
                mask[1::2] = False
                mask[:, 1::2] = False
                frac = np.divide(p, i, out=np.zeros_like(p), where=i > 0)
                ax.quiver(x[mask], y[mask], p[mask]/area*1e6/3*np.sin(angle[mask]),
                          p[mask]/area*1e6/3*np.cos(angle[mask]), frac[mask], cmap='turbo', clim=(0, .3),
                          angles='xy', scale_units='xy', scale=1, pivot='mid',
                          headwidth=0, headlength=0, headaxislength=0)
                ax.set(xlim=(60, -60), ylim=(-60, 60), aspect='equal', title=title,
                       xlabel='East (uas)', ylabel='North (uas)')
                fig.colorbar(mesh, ax=ax, label='I (microJy/uas²)')
            fig.suptitle('Apr11 official HOPS UVFITS + ehtim1.3.2: NOT certified official fiducial')
            fig.savefig(out / 'apr11-reconstruction.png', dpi=160)
            plt.close(fig)
            rounds = [s for s in steps if s['stage'] == 'P+D']
            fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
            for metric in ('pvis', 'm'):
                axes[0].plot([r['round'] for r in rounds], [r['chi2_after_D'][metric] for r in rounds],
                             marker='.', label=metric)
            for metric in ('relative_qu_change', 'max_delta_d'):
                axes[1].semilogy([r['round'] for r in rounds], [r[metric] for r in rounds], label=metric)
            axes[0].set(xlabel='Alternating round', ylabel='Working data-fit statistic')
            axes[1].set(xlabel='Alternating round', ylabel='Change / not a convergence certificate')
            for ax in axes:
                ax.legend()
                ax.grid(alpha=.2)
            fig.savefig(out / 'iteration-diagnostics.png', dpi=160)
            plt.close(fig)
            (out / 'report.json').write_text(json.dumps({'termination': termination,
                'final_round': rounds[-1], 'claim': 'Rebuild output, not official scientific reproduction'},
                indent=2)+'\n')


if __name__ == '__main__':
    main()

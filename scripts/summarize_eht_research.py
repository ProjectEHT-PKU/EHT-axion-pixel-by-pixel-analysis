"""Gather completed, declared observed-research results without changing source products."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from eht_axion.pipeline.common import (
    manifest_path, require_declared_outputs, stage_directory, verification_session, verify_stage,
)


def run(args):
    roots = [p.resolve() for p in args.results]
    inputs, records = [Path(__file__).resolve()], []
    with verification_session():
        seen = set()
        for root in roots:
            manifest = verify_stage(root, recursive=True, allow_imported=False, _seen=seen)
            if manifest['stage'] != 'research:observed-limits':
                raise ValueError('Not an observed-research stage')
            require_declared_outputs(root, [root / 'report.json'])
            report = json.loads((root / 'report.json').read_text())
            inputs += [manifest_path(root), root / 'report.json']
            label = report['variant'] + ('-raw-1d' if report.get('alignment') == 'raw' else '')
            for key, row in report['results'].items():
                path = root / f'{key}.json'
                require_declared_outputs(root, [path])
                inputs.append(path)
                diagnostics = json.loads(path.read_text())
                for i, alpha in enumerate(row['alphas']):
                    fit, boundary = diagnostics['fits'][i], diagnostics['boundaries'][i]
                    records.append(dict(variant=label, result=key, alpha=alpha,
                        c95=row['c95'][i], log10_c95=row['log10_c95'][i], features=row['features'],
                        data_chi2=fit['data_quadratic'], best_chi2=fit['joint_grid_chi2'],
                        delta_chi2=fit['data_quadratic']-fit['joint_grid_chi2'],
                        best_c=fit['joint_grid_coupling'], limit_in_last_bin=boundary['limit_in_last_bin']))
            for dim in report['dimensions']:
                path = root / f'observations-{dim}.npz'
                require_declared_outputs(root, [path])
                inputs.append(path)
        identity = dict(schema='observed-research-comparison-v1', source_roots=list(map(str, roots)),
                        scope='conditional hyperparameter-ensemble research, not calibrated exclusions')
        with stage_directory(args.output_root, 'research:eht-observed-comparison', identity, inputs) as out:
            with (out / 'conditional-results.csv').open('w') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
            for ax, dim in zip(axes, ('1d', '2d')):
                for root in roots:
                    report = json.loads((root / 'report.json').read_text())
                    label = report['variant'] + ('-raw-1d' if report.get('alignment') == 'raw' else '')
                    key = f'observed-cor-{dim}'
                    if key not in report['results']:
                        continue
                    row = report['results'][key]
                    ax.plot(row['alphas'], row['log10_c95'], marker='o', label=label,
                            linestyle='--' if report.get('alignment') == 'raw' else '-')
                ax.set(title=f'{dim}: correlated working covariance', xlabel='Axion alpha',
                       ylabel='log10 c95')
                ax.grid(alpha=.2)
                ax.legend(fontsize=8)
            fig.suptitle('RESEARCH: selection / alignment sensitivity; selected feature sets differ')
            fig.savefig(out / 'selection-sensitivity.png', dpi=160)
            plt.close(fig)
            chosen = roots[0]
            report = json.loads((chosen / 'report.json').read_text())
            fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout='constrained')
            for col, dim in enumerate(('1d', '2d')):
                key = f'observed-cor-{dim}'
                if key not in report['results']:
                    continue
                data = json.loads((chosen / f'{key}.json').read_text())
                index = int(np.argmin(abs(np.asarray(data['alphas'])-.3)))
                fit = data['fits'][index]
                with np.load(chosen / f'observations-{dim}.npz') as z:
                    daydiff = z['pair_differences']
                    sigma = np.sqrt(np.diagonal(z['covariances'][0]+z['covariances'][1]))
                    sigma2 = np.sqrt(np.diagonal(z['covariances'][2]+z['covariances'][3]))
                for pair, errors in enumerate((sigma, sigma2)):
                    ax = axes[pair, col]
                    x = np.arange(len(daydiff[pair]))
                    ax.errorbar(x, np.rad2deg(daydiff[pair]), yerr=np.rad2deg(errors), fmt='o',
                                label='Measured day difference / marginal ensemble spread')
                    ax.plot(x, np.rad2deg(fit['prediction_radians'][pair]), 's-',
                            label='Joint-grid best template')
                    ax.axhline(0, color='k', lw=.5)
                    ax.set(title=f'{dim}, pair {pair+1}, alpha={data["alphas"][index]}',
                           xlabel='Selected feature index (not independent beam)',
                           ylabel='EVPA difference (deg)')
                    ax.grid(alpha=.2)
                    ax.legend(fontsize=7)
            fig.suptitle('RESEARCH: common-phase fit; marginal error bars, full-covariance fit')
            fig.savefig(out / 'daydiff-fit-alpha030.png', dpi=160)
            plt.close(fig)
            (out / 'report.json').write_text(json.dumps(dict(identity, records=len(records),
                note='Different selections retain different features; curves do not isolate selection alone.',
                caution='Delta chi2 is not a calibrated detection significance.'), indent=2)+'\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output_root.exists():
        parser.error('New output root required')
    run(args)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Analysis check: split-MC covariance sensitivity at fixed features and observed center.

Each half estimates covariance, not covariance of the mean. The selected features
and real day differences stay fixed; no mask selection or covariance ridge.
"""
import argparse
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from eht_axion import provenance
from eht_axion.analysis.likelihood import MODES, calculate_observed
from eht_axion.analysis.observations import ObservationData
from eht_axion.analysis.templates import load_templates
from eht_axion.evpa import circular_difference
from eht_axion.experiment import DEFAULT_PATH, read_experiment
from analysis_limits import options_from_config


def covariance(samples, center):
    residual = circular_difference(samples, center)
    residual -= residual.mean(axis=0)
    return residual.T @ residual / len(samples)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=DEFAULT_PATH, help='default: experiment.yaml')
    for name in ('source', 'templates', 'output-root'):
        p.add_argument('--'+name, required=True, type=Path)
    args = p.parse_args(argv)
    if args.output_root.exists():
        raise FileExistsError(f'New output root required: {args.output_root}')
    with provenance.verification_session():
        provenance.verify_stage(args.source, recursive=False)
        provenance.require_declared_outputs(args.source, [args.source/'report.json'])
        source_report = json.loads((args.source/'report.json').read_text())
        config = json.loads(json.dumps(read_experiment(args.config)))
        if config != source_report['scientific_profile']:
            raise ValueError('The source limits were produced from a different experiment')
        draws = config['eht_imaging']['uncertainty_draws']
        half = draws // 2
        inputs = [Path(__file__).resolve(), args.config.resolve(), provenance.manifest_path(args.source),
                  args.source/'report.json', Path(__file__).with_name('analysis_limits.py')]
        products = {}
        for dim in ('1d', '2d'):
            root = args.templates/'templates'/dim
            provenance.verify_stage(root, recursive=False)
            provenance.require_declared_outputs(root, [root/'templates.npz'])
            products[dim] = load_templates(root/'templates.npz', dim)
            path = args.source/f'observations-{dim}.npz'
            provenance.require_declared_outputs(args.source, [path])
            inputs += [provenance.manifest_path(root), root/'templates.npz', path]
        inputs += sorted((Path(__file__).resolve().parents[1]/'eht_axion').rglob('*.py'))
        identity = dict(source=str(args.source), fixed_observed_centers=True, fixed_selection=True,
                        covariance_halves=[list(range(half)), list(range(half, draws))],
                        covariance_regularization='none', scope='finite MC sensitivity, not new error model')
        with provenance.stage_directory(args.output_root, 'analysis:check-covariance',
                                    identity, inputs, verify_lineage=False) as out:
            results, table = {}, []
            for dim in ('1d', '2d'):
                with np.load(args.source/f'observations-{dim}.npz') as z:
                    samples, means = z['evpa_samples_radians'], z['means']
                    differences, selected = z['pair_differences'], z['selected_indices']
                if samples.shape[:2] != (4, draws):
                    raise ValueError('Expected saved four-day selected samples for every draw')
                alpha, amp, phase, _ = products[dim]
                halves = [(f'first{half}', slice(0, half)), (f'last{draws - half}', slice(half, draws))]
                for name, sl in halves:
                    covs = tuple(covariance(s[sl], center) for s, center in zip(samples, means))
                    pair = (covs[0]+covs[1], covs[2]+covs[3])
                    observed = ObservationData(tuple(means), tuple(np.sqrt(np.diag(c)) for c in covs),
                        covs, differences, pair, selected)
                    options = options_from_config(config, source_report['background_sigma_radians'][dim])
                    for mode in MODES:
                        a, cuts, _ = calculate_observed(observed, amp, phase, alphas=alpha,
                            dimension=dim, mode=mode, options=options)
                        baseline = np.asarray(source_report['results'][f'observed-{mode}-{dim}']['log10_c95'])
                        key = f'{name}-{mode}-{dim}'
                        results[key] = dict(alphas=a.tolist(), log10_c95=cuts.tolist(),
                                            c95=(10**cuts).tolist(), delta_dex=(cuts-baseline).tolist())
                        table += [dict(half=name, mode=mode, dimension=dim, alpha=float(x),
                                       c95=float(10**y), delta_dex=float(delta))
                                  for x, y, delta in zip(a, cuts, cuts-baseline)]
                        print(key, 'max_abs_delta_dex', float(np.max(abs(cuts-baseline))), flush=True)
            with (out/'covariance-sensitivity.csv').open('w') as f:
                writer = csv.DictWriter(f, fieldnames=list(table[0]))
                writer.writeheader()
                writer.writerows(table)
            (out/'report.json').write_text(json.dumps(dict(identity, results=results,
                max_abs_delta_dex=max(abs(v) for r in results.values() for v in r['delta_dex'])),
                indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()

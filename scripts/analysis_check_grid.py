#!/usr/bin/env python3
"""Analysis check: re-evaluate saved selected observations at a denser grid / wider prior.

--phase-count, --coupling-points and --coupling-max override analysis.inference;
everything else must be the experiment that produced the source limits. No EHT
reselection, covariance refit or raw RAPTOR read: the source result's outputs are
verified here; its raw lineage was verified when that stage was made, so the
source manifests are recorded as inputs without re-hashing their whole lineage.
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
from eht_axion.experiment import DEFAULT_PATH, read_experiment
from analysis_limits import options_from_config, plot_results


def load_observation(path):
    with np.load(path, allow_pickle=False) as z:
        means = tuple(z['means'])
        std = tuple(z['standard_deviations'])
        cov = tuple(z['covariances'])
        differences, pairs = z['pair_differences'], tuple(z['pair_covariances'])
        selected = z['selected_indices']
    observed = ObservationData(means, std, cov, differences, pairs, selected)
    null = ObservationData(means, std, cov, np.zeros_like(differences), pairs, selected)
    return observed, null


OVERRIDES = ('phase_count', 'coupling_points', 'coupling_max')


def run(args):
    with provenance.verification_session():
        # Verify the complete RESULT outputs. This check need not rehash every
        # unchanged raw RT ancestor once more for a numerical-only repeat.
        record = provenance.verify_stage(args.source, recursive=False)
        if record['stage'] != 'research:eht120-observed-limits':
            raise ValueError('Expected completed primary observed analysis')
        source_config = record['config']['scientific_profile']
        experiment = json.loads(json.dumps(read_experiment(args.config)))
        if experiment != source_config:
            raise ValueError('The source limits were produced from a different experiment')
        config = json.loads(json.dumps(source_config))
        for key in OVERRIDES:
            if getattr(args, key) is not None:
                config['analysis']['inference'][key] = getattr(args, key)
        changed = {k: v for k, v in config['analysis']['inference'].items()
                   if v != source_config['analysis']['inference'][k]}
        if not changed:
            raise ValueError('No grid/prior override differs from the source')
        reportpath = args.source/'report.json'
        obsfiles = [args.source/f'observations-{d}.npz' for d in ('1d', '2d')]
        provenance.require_declared_outputs(args.source, [reportpath, *obsfiles])
        source_report = json.loads(reportpath.read_text())
        products, inputs = {}, [Path(__file__).resolve(), args.config.resolve(),
            provenance.manifest_path(args.source), reportpath, *obsfiles,
            Path(__file__).with_name('analysis_limits.py')]
        for dim in ('1d', '2d'):
            root = args.templates/'raptor/templates'/dim
            provenance.verify_stage(root, recursive=False)
            provenance.require_declared_outputs(root, [root/'templates.npz'])
            products[dim] = load_templates(root/'templates.npz', dim)
            inputs += [provenance.manifest_path(root), root/'templates.npz']
        inputs += sorted((Path(__file__).resolve().parents[1]/'eht_axion').rglob('*.py'))
        identity = dict(changed=changed, source=str(args.source),
            source_verified='completed result outputs + recorded producer inputs; no raw RT rehash',
            same_selected_observations=True, same_covariance=True, canonical_NULL_modified=False)
        with provenance.stage_directory(args.output_root, 'research:eht120-limit-numerics',
                                    identity, inputs, verify_lineage=False) as out:
            results, rows = {}, []
            for dim in ('1d', '2d'):
                observed, null = load_observation(args.source/f'observations-{dim}.npz')
                alpha, amp, phase, _ = products[dim]
                # Numerical check: only the declared grid/prior knobs differ; physical model fixed.
                options = options_from_config(config, source_report['background_sigma_radians'][dim])
                for mode in MODES:
                    for label, data in [('observed', observed), ('null-comparison', null)]:
                        a, cuts, diagnostic = calculate_observed(data, amp, phase, alphas=alpha,
                            dimension=dim, mode=mode, options=options)
                        key = f'{label}-{mode}-{dim}'
                        np.savez(out/(key+'.npz'), alphas=a, log10_c95=cuts, c95=10**cuts,
                            couplings=diagnostic['couplings'], cdf=diagnostic['cdf'],
                            posterior_density=diagnostic['posterior_density'])
                        (out/(key+'.json')).write_text(json.dumps(diagnostic, indent=2, allow_nan=False)+'\n')
                        baseline = np.asarray(source_report['results'][key]['log10_c95'])
                        results[key] = dict(alphas=a.tolist(), c95=(10**cuts).tolist(),
                            log10_c95=cuts.tolist(),
                            delta_log10_c95=(cuts-baseline).tolist(), features=len(data.selected_indices))
                        rows += [dict(center=label, mode=mode, dimension=dim, alpha=float(x),
                            log10_c95=float(y), c95=float(10**y), delta_log10_c95=float(delta))
                            for x, y, delta in zip(a, cuts, cuts-baseline)]
                        print(key, 'max delta dex', float(np.max(abs(cuts-baseline))), flush=True)
            with (out/'limits.csv').open('w') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            plot_results(out, results)
            (out/'report.json').write_text(json.dumps(dict(identity, results=results,
                max_delta_dex=max(abs(v) for r in results.values() for v in r['delta_log10_c95'])),
                indent=2, allow_nan=False)+'\n')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=DEFAULT_PATH, help='default: experiment.yaml')
    for key in ('source', 'templates', 'output-root'):
        p.add_argument('--'+key, required=True, type=Path)
    p.add_argument('--phase-count', type=int)
    p.add_argument('--coupling-points', type=int)
    p.add_argument('--coupling-max', type=float)
    run(p.parse_args(argv))


if __name__ == '__main__':
    main()

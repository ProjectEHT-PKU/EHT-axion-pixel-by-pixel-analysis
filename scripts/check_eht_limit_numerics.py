#!/usr/bin/env python3
"""Re-evaluate saved selected observations/templates at denser grids/prior bounds.

No EHT reselection, covariance refit or raw RAPTOR read. A source-result SHA
closure is checked once here; raw producer closure was verified by that completed
stage, not silently cached/re-signed. Canonical NULL remains untouched.
"""
import argparse
from dataclasses import replace
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from eht_axion.observations import ObservationData
from eht_axion.pipeline import common
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.products import load_templates
from eht_axion.research.ensemble_limits import MODES
from eht_axion.research.observed import calculate_observed
from calculate_eht_ensemble_limits import options_from_config, plot_results


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


def run(args):
    config = read_config(args.config)
    with common.verification_session():
        # Verify the complete RESULT outputs. This check need not rehash every
        # unchanged raw RT ancestor once more for a numerical-only repeat.
        record = common.verify_stage(args.source, recursive=False, allow_imported=False)
        if record['stage'] != 'research:eht120-observed-limits':
            raise ValueError('Expected completed primary observed analysis')
        source_config = record['config']['scientific_profile']
        for key in config:
            if key != 'inference' and config[key] != source_config[key]:
                raise ValueError('Only numerical grid/prior settings may differ')
        changed = {k: v for k, v in config['inference'].items() if v != source_config['inference'][k]}
        if not set(changed) <= {'phase_count', 'coupling_points', 'coupling_max'}:
            raise ValueError('Unexpected changed inference setting')
        reportpath = args.source/'report.json'
        obsfiles = [args.source/f'observations-{d}.npz' for d in ('1d', '2d')]
        common.require_declared_outputs(args.source, [reportpath, *obsfiles])
        source_report = json.loads(reportpath.read_text())
        # Explicit byte-identical source-receipt snapshots keep these numerical
        # repeats from recursively reopening raw RT data via filename dispatch.
        # They are receipts of completed result verification, NOT rebound stages.
        receipts = args.output_root.parent / (args.output_root.name + '-source-receipts')
        receipts.mkdir(exist_ok=False)
        original_manifest = common.manifest_path(args.source)
        snapshot = receipts / 'primary-source-receipt.json'
        snapshot.write_bytes(original_manifest.read_bytes())
        snapshots = [dict(original=str(original_manifest), snapshot=str(snapshot),
                          sha256=common.sha256(snapshot))]
        products, inputs = {}, [Path(__file__).resolve(), args.config.resolve(),
            snapshot, reportpath, *obsfiles,
            Path(__file__).with_name('calculate_eht_ensemble_limits.py')]
        for dim in ('1d', '2d'):
            root = args.templates/'raptor/templates'/dim
            common.verify_stage(root, recursive=False, allow_imported=False)
            common.require_declared_outputs(root, [root/'templates.npz'])
            products[dim] = load_templates(root/'templates.npz', dim)
            original_manifest = common.manifest_path(root)
            snapshot = receipts / (dim + '-template-source-receipt.json')
            snapshot.write_bytes(original_manifest.read_bytes())
            snapshots.append(dict(original=str(original_manifest), snapshot=str(snapshot),
                                  sha256=common.sha256(snapshot)))
            inputs += [snapshot, root/'templates.npz']
        inputs += sorted((Path(__file__).resolve().parents[1]/'eht_axion').rglob('*.py'))
        identity = dict(changed=changed, source=str(args.source), receipt_snapshots=snapshots,
            source_verified='completed result outputs + recorded producer inputs; no raw RT rehash',
            same_selected_observations=True, same_covariance=True, canonical_NULL_modified=False)
        with common.stage_directory(args.output_root, 'research:eht120-limit-numerics',
                                    identity, inputs) as out:
            results, rows = {}, []
            for dim in ('1d', '2d'):
                observed, null = load_observation(args.source/f'observations-{dim}.npz')
                alpha, amp, phase, _ = products[dim]
                options = options_from_config(config, source_report['background_sigma_radians'][dim])
                # Numerical check changes only the declared knobs; physical model fixed.
                options = replace(options, phase_count=config['inference']['phase_count'],
                                  coupling_points=config['inference']['coupling_points'])
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


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('source', 'templates', 'config', 'output-root'):
        p.add_argument('--'+key, required=True, type=Path)
    run(p.parse_args())

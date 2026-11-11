#!/usr/bin/env python3
"""Analysis check: re-evaluate saved selected observations at a denser grid / wider prior.

--phase-count, --coupling-points and --coupling-max override analysis.inference;
everything else must be the experiment that produced the source limits. No EHT
reselection, covariance refit or raw RAPTOR read; the source manifests are recorded
as inputs.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import argparse
import csv
import json

import numpy as np
from eht_axion import provenance
from eht_axion.analysis.likelihood import MODES, calculate_observed
from eht_axion.analysis.observations import ObservationData
from eht_axion.analysis.results import read_results, write_results
from eht_axion.analysis.templates import load_templates
from eht_axion.experiment import DEFAULT_PATH, read_experiment
from scripts.analysis.limits import options_from_config, plot_results
from eht_axion.environment import add_data_root_option, fill_defaults


def load_observation(record):
    """Observed and zero-data ObservationData of one observations record of a limits stage."""
    means = tuple(record['means'])
    std = tuple(record['standard_deviations'])
    cov = tuple(record['covariances'])
    differences, pairs = record['pair_differences'], tuple(record['pair_covariances'])
    selected = record['selected_indices']
    observed = ObservationData(means, std, cov, differences, pairs, selected)
    zero_data = ObservationData(means, std, cov, np.zeros_like(differences), pairs, selected)
    return observed, zero_data


OVERRIDES = ('phase_count', 'coupling_points', 'coupling_max')


def run(args):
    if args.output_root.exists():
        raise FileExistsError(f'New output root required: {args.output_root}')
    record = provenance.verify_stage(args.source)
    if record['stage'] != 'analysis:limits':
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
    files = [args.source/'report.json', args.source/'results.npz']
    provenance.require_declared_outputs(args.source, files)
    source_report, source = read_results(args.source, ('observations',))
    root = args.templates/'templates'
    provenance.verify_stage(root)
    provenance.require_declared_outputs(root, [root/'templates.npz'])
    products = {dim: load_templates(root/'templates.npz', dim) for dim in ('1d', '2d')}
    inputs = [provenance.manifest_path(args.source), *files,
              provenance.manifest_path(root), root/'templates.npz']
    identity = dict(changed=changed, source=str(args.source),
        same_selected_observations=True, same_covariance=True)
    with provenance.stage_directory(args.output_root, 'analysis:check-grid',
                                identity, inputs) as out:
        results, rows, cases = {}, [], {}
        for dim in ('1d', '2d'):
            observed, zero_data = load_observation(source['observations'][dim])
            alpha, amp, phase, _ = products[dim]
            # Numerical check: only the declared grid/prior knobs differ; physical model fixed.
            options = options_from_config(config, source_report['background_sigma_radians'][dim])
            for mode in MODES:
                for label, data in [('observed', observed), ('zero-data', zero_data)]:
                    a, cuts, diagnostic = calculate_observed(data, amp, phase, alphas=alpha,
                        dimension=dim, mode=mode, options=options)
                    key = f'{label}-{mode}-{dim}'
                    cases[key] = dict(diagnostic, log10_c95=cuts)
                    baseline = np.asarray(source_report['results'][key]['log10_c95'])
                    results[key] = dict(alphas=a.tolist(), c95=(10**cuts).tolist(),
                        log10_c95=cuts.tolist(),
                        delta_log10_c95=(cuts-baseline).tolist(), features=len(data.selected_indices))
                    rows += [dict(data=label, mode=mode, dimension=dim, alpha=float(x),
                        log10_c95=float(y), c95=float(10**y), delta_log10_c95=float(delta))
                        for x, y, delta in zip(a, cuts, cuts-baseline)]
                    print(key, 'max delta dex', float(np.max(abs(cuts-baseline))), flush=True)
        with (out/'limits.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        plot_results(out, results)
        write_results(out, dict(identity, results=results,
            max_delta_dex=max(abs(v) for r in results.values() for v in r['delta_log10_c95'])),
            dict(cases=cases))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=DEFAULT_PATH, help='default: configs/experiment.yaml')
    p.add_argument('--source', type=Path, help='scripts/analysis/limits.py stage (default: analysis/limits)')
    p.add_argument('--templates', type=Path, help='templates root (default: analysis/templates)')
    p.add_argument('--output-root', type=Path,
                   help='new directory (default: analysis/check-grid or analysis/check-prior, see --check)')
    p.add_argument('--check', choices=('grid', 'prior'),
                   help='which default output root of the data layout to use when --output-root is not given')
    p.add_argument('--phase-count', type=int)
    p.add_argument('--coupling-points', type=int)
    p.add_argument('--coupling-max', type=float)
    add_data_root_option(p)
    args = p.parse_args(argv)
    if args.output_root is None and args.check is None:
        p.error('Provide --output-root, or --check grid|prior with a data root')
    defaults = {'source': ('analysis.limits', {}), 'templates': ('analysis.templates', {})}
    if args.output_root is None:
        defaults['output_root'] = (f'analysis.check_{args.check}', {})
    fill_defaults(p, args, defaults)
    run(args)


if __name__ == '__main__':
    main()

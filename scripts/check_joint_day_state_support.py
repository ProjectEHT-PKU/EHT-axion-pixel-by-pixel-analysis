#!/usr/bin/env python3
"""Finite empirical-state sensitivity on saved four-day study inputs.

Same EHT centers/covariance, full selected feature vector and prior. Split model
states/windows, not observational data. No covariance ridge or response rescale.
"""
import argparse
import csv
import inspect
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from eht_axion.pipeline import common
from eht_axion.research import joint_days
from eht_axion.research.joint_days import calculate_joint_days, background_moments
from eht_axion.research.time_sampling import H2, H3


def subset_states(states, start, stop):
    """All eight-alpha windows are sliced together; paired B/T remain paired."""
    response = states['evolving_daily'][:, start:stop]
    background = states['background_h3_principal'][start:stop]
    w = np.full(len(background), 1/len(background))
    return response, background, w


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    a = p.parse_args(argv)
    with common.verification_session():
        common.verify_stage(a.inputs, recursive=False, allow_imported=False)
        common.verify_stage(a.reference, recursive=False, allow_imported=False)
        files = [a.inputs/'report.json'] + [a.inputs/f'{name}-{dim}.npz'
                 for name in ['states', 'observations'] for dim in ['1d', '2d']]
        common.require_declared_outputs(a.inputs, files)
        sources = files+[Path(__file__).resolve(), Path(inspect.getfile(joint_days)),
                         a.reference/'summary.csv']
        common.require_declared_outputs(a.reference, [a.reference/'summary.csv'])
        receipts = a.output_root.parent/(a.output_root.name+'-source-receipts')
        receipts.mkdir(exist_ok=False)
        for name, root in [('inputs', a.inputs), ('reference', a.reference)]:
            q = receipts/(name+'-receipt.json')
            q.write_bytes(common.manifest_path(root).read_bytes())
            sources.append(q)
        with (a.reference/'summary.csv').open() as f:
            ref = {(r['dimension'], int(r['contrasts']), r['case'], r['mode']):float(r['c95'])
                   for r in csv.DictReader(f) if float(r['alpha']) == .3}
        rows = []
        config = dict(alpha=.3, EHT_frozen=True, same_prior=True, selection_frozen=True,
                      state_subsets='first/last half windows; each subset re-estimates its own moments',
                      scope='model-state finite-window sensitivity, not independent-data error bars')
        with common.stage_directory(a.output_root, 'research:joint-day-state-support',
                                    config, sources) as out:
            for dim in ['1d', '2d']:
                with np.load(a.inputs/f'states-{dim}.npz') as z:
                    s = {k:z[k] for k in z.files}
                with np.load(a.inputs/f'observations-{dim}.npz') as z:
                    cov = z['covariances']
                idx = np.flatnonzero(np.isclose(s['alphas'], .3))[0]
                windows = len(s['background_h3_principal'])
                for label, first, last in [('first_half', 0, windows//2),
                                           ('last_half', windows//2, windows)]:
                    responses, background, weights = subset_states(s, first, last)
                    for h in [H2, H3]:
                        cols = [0, 2] if len(h) == 2 else [0, 1, 2]
                        b = background[:, cols]
                        mean, cb = background_moments(b, weights=weights)
                        centered = b-mean
                        t = np.einsum('kd,mdn->mkn', h, responses[idx])
                        cases = [
                            ('mean_background_full', dict(signal_coeff=np.mean(
                                np.einsum('kd,mdn->mkn', h, s['frozen_daily'][idx]), axis=0)[None],
                                background_cov=cb)),
                            ('evolving_mixture_full', dict(signal_coeff=t, weights=weights,
                                                          background_cov=cb)),
                            ('joint_background_evolving_mixture', dict(signal_coeff=t,
                                weights=weights, backgrounds=centered)),
                            ('joint_gaussian_moments', dict(signal_coeff=t, weights=weights,
                                                           backgrounds=centered, model='gaussian'))]
                        for case, settings in cases:
                            for target in ['observed', 'matched_zero']:
                                daily = s['observed_daily_local_unwrapped'] if target == 'observed' else (
                                    np.zeros_like(s['observed_daily_local_unwrapped']))
                                r = calculate_joint_days(daily, cov, h,
                                    couplings=np.logspace(np.log10(1/137), 3, 120),
                                    phases=np.arange(64)*2*np.pi/64, **settings)
                                baseline = ref[dim, len(h), case, target]
                                key = f'{dim}-h{len(h)}-{label}-{case}-{target}'
                                np.savez(out/(key+'.npz'), c95=r['c95'], couplings=r['couplings'],
                                         density=r['density'], cdf=r['cdf'], loglike=r['loglike'])
                                row = dict(dimension=dim, contrasts=len(h), subset=label, case=case,
                                    mode=target, states=len(b), c95=r['c95'], baseline_c95=baseline,
                                    delta_dex=float(np.log10(r['c95']/baseline)))
                                rows.append(row)
                                print(row, flush=True)
            with (out/'sensitivity.csv').open('w') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            (out/'report.json').write_text(json.dumps(dict(config, results=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()

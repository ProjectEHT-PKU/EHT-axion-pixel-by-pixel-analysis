"""Additional 24h embargo after raw-frame purge; no posterior/EHT fitting.

Tests sensitivity to temporal separation. A gap is not proof of independence.
Same fixed zero predictive mean, covariance candidates and frame coordinates.
"""
import argparse
import inspect
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from scipy.linalg import solve_triangular
from eht_axion.pipeline import common
from eht_axion.research import background_validation as validation


def embargo_training(support, train, test, times, hours):
    lo, hi = times[support[test, 0].min()], times[support[test, 1].max()]
    return train[(times[support[train, 1]] < lo-hours) |
                 (times[support[train, 0]] > hi+hours)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--embargo-hours', type=float, default=24.)
    a = p.parse_args(argv)
    if not np.isfinite(a.embargo_hours) or a.embargo_hours < 0:
        raise ValueError('nonnegative finite embargo required')
    with common.verification_session():
        common.verify_stage(a.inputs, recursive=False, allow_imported=False)
        files = [a.inputs/f'states-{dim}.npz' for dim in ['1d', '2d']]
        common.require_declared_outputs(a.inputs, files)
        receipts = a.output_root.parent/(a.output_root.name+'-source-receipts')
        receipts.mkdir(exist_ok=False)
        snapshot = receipts/'source-receipt.json'
        snapshot.write_bytes(common.manifest_path(a.inputs).read_bytes())
        sources = files+[snapshot, Path(__file__).resolve(), Path(inspect.getfile(validation))]
        records = []
        for dim, file in zip(['1d', '2d'], files):
            with np.load(file) as z:
                data, starts, times = (z['background_h3_principal'],
                    z['window_starts_hours'], z['snapshot_times_hours'])
            support = validation.support_indices(starts, times)
            for count in [5, 3]:
                folds = validation.purged_folds(starts, times, n_splits=count)
                for k, rows in [(2, [0, 2]), (3, [0, 1, 2])]:
                    b = data[:, rows]
                    for number, fold in enumerate(folds):
                        train = embargo_training(support, fold['train'], fold['test'], times,
                                                 a.embargo_hours)
                        for name in validation.MAIN_MODELS+validation.DIAGNOSTIC_MODELS:
                            c, meta = validation.model_covariance(b[train], name)
                            r = dict(dimension=dim, contrasts=k, n_splits=count, fold=number,
                                     model=name, n_train=len(train), n_test=len(fold['test']),
                                     scoreable=meta['scoreable'])
                            if meta['scoreable']:
                                t = b[fold['test']].reshape(len(fold['test']), -1)
                                chol = np.linalg.cholesky(c)
                                w = solve_triangular(chol, t.T, lower=True)
                                q = (w*w).sum(axis=0)
                                nll = .5*(t.shape[1]*np.log(2*np.pi)+
                                          2*np.log(np.diag(chol)).sum()+q)/t.shape[1]
                                r['mean_nll_per_dimension'] = float(nll.mean())
                            records.append(r)
        config = dict(embargo_hours=a.embargo_hours, EHT_used=False, hyperparameter_tuning=False,
                      note='extra separation after raw-frame purge; not a claim of independence')
        with common.stage_directory(a.output_root, 'research:background-embargo-check',
                                    config, sources) as out:
            (out/'report.json').write_text(json.dumps(dict(config, records=records),
                                                      indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()

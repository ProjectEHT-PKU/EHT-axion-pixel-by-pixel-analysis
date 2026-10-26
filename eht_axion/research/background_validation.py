"""Finite, background-only covariance validation; no EHT likelihood or I/O.

References: sklearn TimeSeriesSplit official API (``gap`` excludes training
samples before test): https://scikit-learn.org/stable/modules/generated/
sklearn.model_selection.TimeSeriesSplit.html ; Roberts et al. (2017), Ecography
40, 913–929, https://doi.org/10.1111/ecog.02881 (dependent-data blocking).
The sklearn HTML was consulted; the Roberts publisher HTML was access-blocked.
Our split is NOT TimeSeriesSplit: training may lie on either side of a held-out
contiguous block. It assesses stationary background prediction, not forecasting.
Purging uses original snapshot indices, not merely window labels. This removes
shared raw support, not all temporal dependence. No iid standard errors, p-values,
one-SE rule, or automatic model selection are supplied.

Samples must be the prepared, fixed-branch (M,K,N) backgrounds, NOT globally
centered inputs. Every training covariance uses its own centered population
moment (ddof=0), but the declared predictive mean is zero for ALL models.
Train-mean drift is reported, not fitted away at prediction time. No observation
noise, eigenvalue floors, diagonal rescue or coefficient tuning is allowed.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_triangular

MAIN_MODELS = ('iso', 'temporal_iso', 'partial_diag', 'shrink_full_0.25')
DIAGNOSTIC_MODELS = ('diagonal', 'full')


def _real(value, name):
    a = np.asarray(value)
    if a.dtype.kind not in 'fiu' or not a.size or not np.isfinite(a).all():
        raise ValueError(f'{name} must be nonempty, finite, real numeric data')
    return a.astype(float)


def support_indices(starts, snapshot_times, day_hours=(0., 24., 120., 144.),
                    width=6., *, window_centers=None):
    """Inclusive [first,last] original snapshot indices, all times in HOURS.

    Enclose min(center)-width/2 through max(center)+width/2 in floor/ceil
    original grid knots. Exact knots need no extra bracket. Interior snapshots
    are conservatively included even between daily exposures. Prepared centers
    may be supplied explicitly; otherwise centers = starts[:,None]+day_hours.
    Never truncate an exposure extending outside the snapshot range.
    """
    starts = _real(starts, 'starts')
    times = _real(snapshot_times, 'snapshot_times')
    days = _real(day_hours, 'day_hours')
    if (starts.ndim != 1 or times.ndim != 1 or days.ndim != 1
            or len(times) < 2 or np.any(np.diff(times) <= 0)
            or np.any(np.diff(starts) <= 0)):
        raise ValueError('require increasing 1D starts and snapshot_times, 1D day_hours')
    if not np.isscalar(width) or not np.isfinite(width) or width <= 0:
        raise ValueError('width must be positive finite hours')
    centers = (starts[:, None] + days if window_centers is None
               else _real(window_centers, 'window_centers'))
    if centers.shape != (len(starts), len(days)):
        raise ValueError('window_centers must have shape (windows, days)')
    lo, hi = centers.min(axis=1)-width/2, centers.max(axis=1)+width/2
    if np.any(lo < times[0]) or np.any(hi > times[-1]):
        raise ValueError('exposure outside original snapshot support')
    first = np.searchsorted(times, lo, side='right') - 1
    last = np.searchsorted(times, hi, side='left')
    return np.column_stack((first, last))


def purged_folds(starts, snapshot_times, day_hours=(0., 24., 120., 144.),
                 width=6., n_splits=5, *, window_centers=None):
    """Return contiguous array_split test blocks with support-disjoint training.

    Dictionaries expose train/test/purged window indices and inclusive raw
    support bounds. Sharing even one interpolation bracket causes purging.
    """
    support = support_indices(starts, snapshot_times, day_hours, width,
                              window_centers=window_centers)
    m = len(support)
    if (isinstance(n_splits, bool) or not isinstance(n_splits, (int, np.integer))
            or not 2 <= n_splits <= m):
        raise ValueError('n_splits must be an integer in [2, windows]')
    indices = np.arange(m)
    folds = []
    for test in np.array_split(indices, n_splits):
        first, last = support[test, 0].min(), support[test, 1].max()
        train = indices[(support[:, 1] < first) | (support[:, 0] > last)]
        purged = np.setdiff1d(indices, np.concatenate((train, test)))
        folds.append(dict(train=train, test=test, purged=purged,
                          support_indices=support.copy(),
                          test_support=np.array([first, last])))
    return folds


def model_covariance(samples, name):
    """Return (covariance, metadata) from training (M,K,N) samples only.

    temporal_iso pools Ctime[p,q] = mean_i C[(p,i),(q,i)] and returns
    Ctime ⊗ I_N (contrast-major ordering). shrink_full_0.25 is fixed
    0.75*temporal_iso + 0.25*empirical_full. partial_diag is fixed half
    empirical diagonal plus half pooled variance. Parameter counts describe
    fitted moments, not effective independent observations or free CV tuning.
    Singular PSD matrices are returned as unscoreable, never regularized.
    """
    b = _real(samples, 'samples')
    if b.ndim != 3 or len(b) < 2:
        raise ValueError('samples must have shape (M>=2,K,N)')
    if name not in MAIN_MODELS + DIAGNOSTIC_MODELS:
        raise ValueError(f'unknown covariance model: {name}')
    m, k, n = b.shape
    d = k*n
    mean = b.mean(axis=0)
    x = (b-mean).reshape(m, d)
    full = x.T @ x / m
    variances = np.diag(full)
    pooled = variances.mean()
    time = np.einsum('mpn,mqn->pq', b-mean, b-mean)/(m*n)
    temporal = np.kron(time, np.eye(n))
    counts = dict(iso=1, temporal_iso=k*(k+1)//2, partial_diag=d,
                  diagonal=d, full=d*(d+1)//2,
                  shrink_full_0_25=d*(d+1)//2)
    if name == 'iso':
        cov = pooled*np.eye(d)
    elif name == 'temporal_iso':
        cov = temporal
    elif name == 'partial_diag':
        cov = np.diag(.5*variances + .5*pooled)
    elif name == 'diagonal':
        cov = np.diag(variances)
    elif name == 'full':
        cov = full
    else:
        cov = .75*temporal + .25*full
    rank = int(np.linalg.matrix_rank(cov))
    spd = rank == d
    if spd:
        try:
            np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            spd = False
    meta = dict(model=name, n_train=m, dimension=d,
                parameter_count=counts[name.replace('.', '_')], rank=rank,
                empirical_rank_bound=min(m-1, d),
                condition_number=float(np.linalg.cond(cov)), scoreable=spd,
                failure=None if spd else 'singular_or_not_positive_definite',
                train_mean=mean, train_mean_rms=float(np.sqrt(np.mean(mean**2))),
                covariance_center='training_mean', predictive_mean='zero', ddof=0)
    return cov, meta


def evaluate_folds(samples, starts, snapshot_times, contrast_rows, n_splits=5,
                   *, day_hours=(0., 24., 120., 144.), width=6.,
                   window_centers=None):
    """Return fold/model metric records and descriptive summaries (no winner).

    contrast_rows are indices into the prepared contrast axis (H2: [0,2],
    H3: [0,1,2]). The entire four-day support is retained even for subrows.
    NLL and quadratic are per held-out window divided by vector dimension;
    quadratic is NOT a calibrated chi-square statistic. Pooled NLL weights
    folds by actual test counts; median/worst and paired differences to iso
    describe instability, not independent-replicate uncertainty. An unscoreable
    fold makes the aggregate unavailable; it is not silently dropped.
    """
    b = _real(samples, 'samples')
    rows = np.asarray(contrast_rows)
    if (b.ndim != 3 or rows.ndim != 1 or not rows.size
            or rows.dtype.kind not in 'iu' or len(np.unique(rows)) != len(rows)
            or np.any(rows < 0) or np.any(rows >= b.shape[1])):
        raise ValueError('samples must be (M,K,N), contrast_rows unique valid indices')
    folds = purged_folds(starts, snapshot_times, day_hours, width, n_splits,
                         window_centers=window_centers)
    if len(b) != len(folds[0]['support_indices']):
        raise ValueError('samples and starts must have equal window counts')
    b = b[:, rows, :]
    records = []
    for number, fold in enumerate(folds):
        for name in MAIN_MODELS + DIAGNOSTIC_MODELS:
            record = dict(fold=number, model=name, n_test=len(fold['test']),
                          n_train=len(fold['train']), n_purged=len(fold['purged']))
            if len(fold['train']) < 2:
                record.update(scoreable=False, failure='fewer_than_two_training_windows')
            else:
                cov, meta = model_covariance(b[fold['train']], name)
                record.update(meta)
                if meta['scoreable']:
                    residual = b[fold['test']].reshape(len(fold['test']), -1)
                    chol = np.linalg.cholesky(cov)
                    whitened = solve_triangular(chol, residual.T, lower=True)
                    q = np.sum(whitened**2, axis=0)
                    d = residual.shape[1]
                    nll = .5*(d*np.log(2*np.pi)
                               + 2*np.log(np.diag(chol)).sum() + q)/d
                    record.update(nll_per_dimension=nll,
                                  mean_nll_per_dimension=float(nll.mean()),
                                  mean_quadratic_per_dimension=float(q.mean()/d))
            records.append(record)
    summary = {}
    baseline = [r for r in records if r['model'] == 'iso']
    for name in MAIN_MODELS + DIAGNOSTIC_MODELS:
        rr = [r for r in records if r['model'] == name]
        valid = all(r['scoreable'] for r in rr)
        item = dict(all_folds_scoreable=valid,
                    scoreable_folds=sum(r['scoreable'] for r in rr),
                    diagnostic_only=name in DIAGNOSTIC_MODELS)
        if valid:
            scores = np.array([r['mean_nll_per_dimension'] for r in rr])
            item.update(weighted_mean_nll_per_dimension=float(np.average(
                scores, weights=[r['n_test'] for r in rr])),
                median_fold_nll_per_dimension=float(np.median(scores)),
                worst_fold_nll_per_dimension=float(scores.max()))
            if all(r['scoreable'] for r in baseline):
                delta = scores - [r['mean_nll_per_dimension'] for r in baseline]
                item.update(fold_nll_difference_to_iso=delta,
                            win_fraction_vs_iso=float(np.mean(delta < 0)))
        summary[name] = item
    return dict(folds=folds, records=records, summary=summary,
                predictive_mean='zero', validation='two_sided_purged_stationarity',
                uncertainty='descriptive_only_temporally_dependent',
                selected_model=None)

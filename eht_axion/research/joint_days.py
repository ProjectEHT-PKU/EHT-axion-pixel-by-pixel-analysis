"""Conditional local-branch four-day likelihoods with frozen EHT covariance.

All vectors flatten in (contrast, feature) order. Days are independent under
this observation model: C[(k,i),(l,j)] = sum_d H[k,d] H[l,d] C_daily[d,i,j].
The caller selects a consistent unwrapped angular branch; nothing here wraps,
reselects features, estimates EHT errors, or calibrates covariance against chi².
Templates are normalized at unit coupling, INCLUDING any physical b prefactor
supplied by the caller. The convention is c * Re(exp(i*psi) * template).
A latent index denotes ONE whole state, shared across all contrasts/features.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp

from ..statistics import normalized_posterior_cdf, posterior_upper_limit


def _array(value, name, *, complex_ok=False):
    raw = np.asarray(value)
    if raw.dtype.kind not in ("fiuc" if complex_ok else "fiu"):
        raise ValueError(f"{name} must be numeric (not bool) and real unless a template")
    result = raw.astype(complex if complex_ok else float)
    if not np.isfinite(result).all() or not result.size:
        raise ValueError(f"{name} must be nonempty and finite")
    return result


def _spd(value, size, name):
    matrix = _array(value, name)
    if matrix.shape != (size, size):
        raise ValueError(f"{name} must have shape {(size, size)}")
    if not np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-14):
        raise ValueError(f"{name} must be symmetric")
    return matrix, np.linalg.cholesky(matrix)


def contrast_observation(daily_means, daily_cov, H):
    """Return (H @ means, full contrast covariance), without angular wrapping.

    Daily means/covariances have shapes (4,N)/(4,N,N). H has K independent
    zero-sum rows (K <= 3). Shared-day off-diagonal contrast blocks are retained.
    """
    means = _array(daily_means, "daily_means")
    cov = _array(daily_cov, "daily_cov")
    h = _array(H, "H")
    if means.ndim != 2 or means.shape[0] != 4:
        raise ValueError("daily_means must have shape (4,N)")
    n = means.shape[1]
    if cov.shape != (4, n, n):
        raise ValueError("daily_cov must have shape (4,N,N)")
    if (h.ndim != 2 or h.shape[1] != 4 or not 1 <= h.shape[0] <= 3
            or not np.allclose(h.sum(axis=1), 0, rtol=0, atol=1e-12)
            or np.linalg.matrix_rank(h) != len(h)):
        raise ValueError("H must have 1..3 independent zero-sum rows and four columns")
    for block in cov:
        _spd(block, n, "daily_cov block")
    joint = np.einsum("kd,ld,dij->kilj", h, h, cov).reshape(len(h)*n, len(h)*n)
    _spd(joint, len(h)*n, "contrast covariance")
    return h @ means, joint


def background_moments(backgrounds, *, weights=None, covariance="full"):
    """Empirical predictive moments of (M,K,N) states, with population ddof=0.

    No division by state count: this is not uncertainty of an estimated mean.
    scalar pools centered feature variances; diagonal keeps feature variances;
    full includes every cross-contrast covariance. No shrinkage is applied.
    """
    b = _array(backgrounds, "backgrounds")
    if b.ndim != 3:
        raise ValueError("backgrounds must have shape (M,K,N)")
    w = _weights(weights, len(b))
    flat = b.reshape(len(b), -1)
    mean = w @ flat
    centered = flat - mean
    full = (centered.T * w) @ centered
    if covariance == "full":
        result = full
    elif covariance == "diagonal":
        result = np.diag(np.diag(full))
    elif covariance == "scalar":
        result = np.eye(flat.shape[1]) * np.trace(full) / flat.shape[1]
    elif covariance == "none":
        result = np.zeros_like(full)
    else:
        raise ValueError("covariance must be none, scalar, diagonal or full")
    return mean.reshape(b.shape[1:]), result


def _weights(weights, count):
    if weights is None:
        return np.full(count, 1 / count)
    w = _array(weights, "weights")
    if w.shape != (count,) or np.any(w <= 0):
        raise ValueError("weights must be positive with shape (M,)")
    # Scaling first also avoids overflow for large finite unnormalized weights.
    w = w / w.max()
    return w / w.sum()


def calculate_joint_days(
    daily_means, daily_cov, H, responses=None, couplings=None, phases=None, *,
    signal_coeff=None, backgrounds=None, background_cov=None,
    background_mode="none", background_mean=None, center_background=False,
    weights=None, model="empirical", prior="log-uniform",
):
    """Evaluate a full 2N/3N Gaussian observation model (also supports K=1).

    Supply exactly one of responses (M,4,N) or signal_coeff (M,K,N), complex;
    the latter permits an unchanged legacy mean pair template (M=1). A single
    phase psi and latent state m govern the ENTIRE observation vector. Uniform
    phase quadrature uses the supplied nodes as-is (no sorting or grid tuning).

    backgrounds: optional paired (M,K,N) latent offsets. Their empirical mean
    is retained by default, not treated as a calibration correction. Explicit
    center_background=True subtracts their weighted mean and reports the choice.
    background_mean: optional fixed (K,N) offset, e.g. from background_moments.
    background_mode: none/scalar/diagonal/full estimates a fixed ddof=0 Gaussian
    covariance from backgrounds. In this mode offsets are replaced by their
    mean and response-state mixing is independent of that Gaussian background.
    Alternatively background_cov accepts a fixed PSD (K*N,K*N) covariance,
    including a legacy scalar variance times identity. PSD validation allows
    roundoff-scale negative eigenvalues but never clips or repairs them.
    Empirical paired offsets
    cannot be combined with additional background covariance (double counting).

    model='empirical' averages normalized Gaussian densities across whole
    states and phases, with fixed observation covariance. model='gaussian'
    moment-matches B+c*Re(exp(i*psi)*T) at EACH c,psi, retaining linear B/T
    cross terms, full c² response covariance and parameter-dependent logdet.
    This is an explicit empirical predictive model, not independent-time data.

    Returns a dict containing loglike, density, CDF, c95, phase_loglike,
    component_min_chi2 and diagnostics. Component minimum chi² is descriptive,
    NOT a mixture goodness-of-fit statistic. Marginal deviance is -2 loglike.
    Prior support is exactly the supplied positive increasing coupling grid.
    No files are read/written and no observation/posterior draws are generated.
    """
    data, eht_cov = contrast_observation(daily_means, daily_cov, H)
    k, n = data.shape
    d = data.size
    h = _array(H, "H")
    if (responses is None) == (signal_coeff is None):
        raise ValueError("supply exactly one of responses or signal_coeff")
    if responses is not None:
        t = _array(responses, "responses", complex_ok=True)
        if t.ndim != 3 or t.shape[1:] != (4, n):
            raise ValueError("responses must have shape (M,4,N)")
        t = np.einsum("kd,mdn->mkn", h, t)
    else:
        t = _array(signal_coeff, "signal_coeff", complex_ok=True)
        if t.ndim != 3 or t.shape[1:] != (k, n):
            raise ValueError("signal_coeff must have shape (M,K,N)")
    m = len(t)
    w = _weights(weights, m)
    c = _array(couplings, "couplings")
    psi = _array(phases, "phases")
    if c.ndim != 1 or len(c) < 2 or np.any(c <= 0) or np.any(np.diff(c) <= 0):
        raise ValueError("couplings must be positive and strictly increasing, length >= 2")
    if psi.ndim != 1:
        raise ValueError("phases must be a nonempty 1-D quadrature grid")
    if model not in ("empirical", "gaussian"):
        raise ValueError("model must be empirical or gaussian")
    if background_mode not in ("none", "scalar", "diagonal", "full"):
        raise ValueError("invalid background_mode")
    if not isinstance(center_background, (bool, np.bool_)):
        raise ValueError("center_background must be bool")
    b = np.zeros((m, k, n))
    if backgrounds is not None:
        b = _array(backgrounds, "backgrounds")
        if b.shape != (m, k, n):
            raise ValueError("backgrounds must match signal_coeff state and contrast shape")
    elif center_background or background_mode != "none":
        raise ValueError("backgrounds required for centering or background_mode")
    empirical_mean = np.einsum("m,mkn->kn", w, b)
    if center_background:
        b = b - empirical_mean
    if background_cov is not None and (backgrounds is not None or background_mode != "none"):
        raise ValueError("background_cov cannot be combined with empirical backgrounds")
    fixed = eht_cov.copy()
    if background_mode != "none":
        mean, cbg = background_moments(b, weights=w, covariance=background_mode)
        fixed += cbg
        b = np.broadcast_to(mean, b.shape).copy()
    if background_cov is not None:
        cbg = _array(background_cov, "background_cov")
        if (cbg.shape != (d, d)
                or not np.allclose(cbg, cbg.T, rtol=1e-12, atol=1e-14)
                or np.linalg.eigvalsh(cbg)[0] < -1e-12 * max(1., np.linalg.norm(cbg, ord=2))):
            raise ValueError("background_cov must be symmetric positive semidefinite (K*N,K*N)")
        fixed += cbg
    if background_mean is not None:
        mean = _array(background_mean, "background_mean")
        if mean.shape != (k, n):
            raise ValueError("background_mean must have shape (K,N)")
        b += mean
    fixed, chol = _spd(fixed, d, "fixed covariance")
    signals = np.real(np.exp(1j * psi[:, None, None]) * t.reshape(m, d)[None, :, :])
    b = b.reshape(m, d)
    y = data.ravel()
    phase_like = np.empty((len(c), len(psi)))
    min_chi2 = np.empty_like(phase_like)
    norm = d * np.log(2 * np.pi) + 2 * np.log(np.diag(chol)).sum()
    if model == "empirical":
        # Only P*M*D and P*M temporaries; never a C*P*M*D residual cube.
        residual = solve_triangular(chol, (y - b).T, lower=True).T
        whitened = solve_triangular(chol, signals.reshape(-1, d).T, lower=True)
        whitened = whitened.T.reshape(signals.shape)
        q0 = np.einsum("md,md->m", residual, residual)
        cross = np.einsum("md,pmd->pm", residual, whitened)
        power = np.einsum("pmd,pmd->pm", whitened, whitened)
        for i, coupling in enumerate(c):
            q = q0[None, :] - 2 * coupling * cross + coupling**2 * power
            q = np.maximum(q, 0)  # quadratic cancellation roundoff, not covariance repair
            phase_like[i] = logsumexp(-.5 * (q + norm) + np.log(w), axis=1)
            min_chi2[i] = q.min(axis=1)
    else:
        for i, coupling in enumerate(c):
            for j, signal in enumerate(signals):
                states = b + coupling * signal
                mean = w @ states
                centered = states - mean
                covariance = fixed + (centered.T * w) @ centered
                factor = np.linalg.cholesky(covariance)
                residual = solve_triangular(factor, y - mean, lower=True)
                q = residual @ residual
                phase_like[i, j] = -.5 * (
                    q + d * np.log(2 * np.pi) + 2 * np.log(np.diag(factor)).sum())
                min_chi2[i, j] = q
    loglike = logsumexp(phase_like, axis=1) - np.log(len(psi))
    density, cdf = normalized_posterior_cdf(c, loglike, prior=prior)
    eigenvalues = np.linalg.eigvalsh(fixed)
    best = np.unravel_index(np.argmax(phase_like), phase_like.shape)
    return dict(
        couplings=c.copy(), phases=psi.copy(), loglike=loglike,
        phase_loglike=phase_like, density=density, cdf=cdf,
        c95=posterior_upper_limit(c, loglike, prior=prior),
        component_min_chi2=min_chi2, marginal_deviance=-2 * loglike,
        data=data, covariance=fixed,
        diagnostics=dict(
            eigenvalue_min=float(eigenvalues[0]), condition_number=float(eigenvalues[-1]/eigenvalues[0]),
            prior=prior, prior_bounds=(float(c[0]), float(c[-1])),
            edge_interval_mass=(float(cdf[1]), float(1 - cdf[-2])),
            state_count=m, effective_state_count=float(1 / (w @ w)),
            effective_count_is_not_independent_data=True,
            background_empirical_mean=empirical_mean,
            center_background=bool(center_background), model=model,
            best_grid_coupling=float(c[best[0]]), best_grid_phase=float(psi[best[1]]),
            minimum_component_chi2=float(min_chi2.min()),
            angular_branch="caller-supplied unwrapped local Gaussian; no wrapping",
        ),
    )

"""Exact convex fixed-I/fixed-D polarimetric (P) imaging solver: diagonal PDHG with a
floating-point duality-gap certificate (NumPy only).

J = official P objective + 2 + 100 log(2), pvis=m=ptv=1, hw=100.
Complex inner products are Re(vdot). TV has exterior-zero right/bottom rows.
Numerical primal/dual certificates are floating-point, NOT interval proofs or
EVPA/image uncertainty. Source lineage: eht-cartesian-exact.py,
eht-primaldual.py, eht-primaldual-precond.py (ported, not imported).

Schur proof: a=abs(K), c=column sums, d=row sums. tau*c=theta;
sigma*d<=theta (the two TV channels share theta/max(dx,dy)). Weighted
Cauchy-Schwarz gives ||sqrt(sigma) K sqrt(tau)|| <= theta=.99 < 1.
Entropy h has radial/tangential curvature >=1, hence H=sum(w*h) is
strongly convex with modulus min(w)>0 in x=P/I. The minimizer is unique;
this is not a scientific acceptance criterion. No sky-angle or D heuristic
is used to stop a P solve.
"""
import time
import numpy as np

THETA = .99


def D(p, shape):
    p = np.asarray(p).reshape(shape)
    dx, dy = -p.copy(), -p.copy()
    dx[:, :-1] += p[:, 1:]
    dy[:-1, :] += p[1:, :]
    return np.array([dx, dy])


def DH(r):
    out = -r[0]-r[1]
    out[:, 1:] += r[0, :, :-1]
    out[1:, :] += r[1, :-1, :]
    return out.ravel()


def project(r):
    return r / np.maximum(1., np.sqrt(np.sum(abs(r)**2, axis=0)))


def prox(v, tau, w):
    """Exact radial entropy prox via safeguarded vector Newton/bisection.

    m+tau*w*atanh(m)=|v|. nextafter is representability, not smoothing.
    """
    v = np.asarray(v, complex)
    a = abs(v)
    hi = np.minimum(a, np.nextafter(1., 0.))
    lo = np.zeros_like(a)
    m = hi/2
    for _ in range(60):
        f = m+tau*w*np.arctanh(m)-a
        hi = np.where(f > 0, m, hi)
        lo = np.where(f <= 0, m, lo)
        trial = m-f/(1+tau*w/(1-m*m))
        m = np.where((trial > lo) & (trial < hi), trial, (lo+hi)/2)
        if np.max(hi-lo) < 4e-15:
            break
    return np.divide(v, a, out=np.zeros_like(v), where=a != 0)*m


def h(m):
    m = np.asarray(m)
    if np.any(m < 0) or np.any(m > 1):
        raise ValueError('entropy requires the physical disk')
    # Continuous endpoint value, without 0*log(0).
    t = 1-m
    return .5*((1+m)*np.log1p(m)+t*np.log(np.where(t > 0, t, 1.)))


def conjugate(v, w):
    t = abs(v)/w
    return float(np.sum(w*(np.logaddexp(t, -t)-np.log(2.))))


class CartesianProblem:
    def __init__(self, intensity, shape, data, tvnorm, flux):
        self.I, self.shape = np.array(intensity, float, copy=True).ravel(), tuple(shape)
        self.flux, self.tvnorm = float(flux), float(tvnorm)
        if (self.I.size != np.prod(shape) or not np.isfinite(self.I).all()
                or np.any(self.I <= 0) or not np.isfinite([flux, tvnorm]).all()
                or min(flux, tvnorm) <= 0):
            raise ValueError('positive finite fixed I, flux and TV norm required')
        self.data = []
        for name in ('pvis', 'm'):
            d, sigma, A = map(np.asarray, data[name])
            if (d.ndim != 1 or not len(d) or sigma.shape != d.shape
                    or A.shape != (len(d), self.I.size)
                    or not all(np.isfinite(v).all() for v in (d, sigma, A))
                    or np.any(sigma <= 0)):
                raise ValueError('invalid official data tuple')
            if name == 'm':
                ai = A @ self.I
                if np.any(ai == 0):
                    raise ValueError('undefined m: A I = 0')
                A = A / ai[:, None]
            self.data.append((d, sigma**2, A, A.conj().T.copy()))

    def physical(self, P):
        P = np.asarray(P, complex).ravel()
        m = abs(P)/self.I
        if not np.isfinite(P).all() or np.any(m > 1):
            raise ValueError('P outside physical disk; no clipping')
        f, g = 0., np.zeros_like(P)
        for d, variance, A, AH in self.data:
            residual = A @ P-d
            f += np.sum(abs(residual)**2/variance)/(2*len(d))
            g += AH @ (residual/variance)/len(d)
        f += 100*np.sum(self.I*h(m))/self.flux
        unit = np.divide(P, abs(P), out=np.zeros_like(P), where=abs(P) != 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            g += 100/self.flux*np.arctanh(m)*unit
        t = D(P, self.shape)
        den = np.sqrt(np.sum(abs(t)**2, axis=0))
        u = np.divide(t, den, out=np.zeros_like(t), where=den != 0)
        return float(f+den.sum()/self.tvnorm), g+DH(u)/self.tvnorm


def exact_cartesian_problem(fixed, obs, library=None):
    """Use official cached tuples; constructing Imager never runs an optimizer."""
    if library is None:
        from .calibration import _library
        library = _library()
    if hasattr(library, '_data_tuples'):
        im = library
    else:
        initial = fixed.copy()
        initial.add_qu((.1*fixed.imvec).reshape(fixed.ydim, fixed.xdim),
                       (.1*fixed.imvec).reshape(fixed.ydim, fixed.xdim))
        initial.add_pol_image(np.zeros((fixed.ydim, fixed.xdim)), 'V')
        im = library.imager.Imager(obs, initial, prior_im=initial,
            flux=float(fixed.total_flux()), pol='P', transform=['mcv'],
            data_term={'pvis': 1, 'm': 1}, reg_term={'hw': 100, 'ptv': 1},
            norm_reg=True, norm_init=False, ttype='direct')
    if (im.dat_term_next != {'pvis': 1, 'm': 1}
            or im.reg_term_next != {'hw': 100, 'ptv': 1}
            or not im.norm_reg or im._ttype != 'direct' or not np.all(im._embed_mask)
            or im.chisq_transform or im.chisq_offset_gradient != 0):
        raise ValueError('unsupported official objective/weights/mask')
    np.testing.assert_array_equal(im.prior_next.imvec, fixed.imvec)
    flux = float(fixed.total_flux())
    if im.flux_next != flux:
        raise ValueError('official flux differs from fixed I')
    return CartesianProblem(fixed.imvec, (fixed.ydim, fixed.xdim),
                            im._data_tuples, flux*fixed.psize/im.beam_size, flux)


class Split:
    def __init__(self, problem):
        self.p, self.I = problem, problem.I
        self.w = 100*self.I/problem.flux
        self.B = np.concatenate([A*self.I[None, :]/np.sqrt(var[:, None]*len(d))
                                 for d, var, A, _ in problem.data])
        self.BH = self.B.conj().T.copy()
        self.b = np.concatenate([d/np.sqrt(var*len(d)) for d, var, _, _ in problem.data])

    def T(self, x):
        return D(self.I*x, self.p.shape)/self.p.tvnorm

    def TH(self, r):
        return self.I*DH(r)/self.p.tvnorm

    def dual(self, x):
        t = self.T(x)
        den = np.sqrt(np.sum(abs(t)**2, axis=0))
        return self.B@x-self.b, np.divide(t, den, out=np.zeros_like(t), where=den != 0)

    def cert(self, x, q=None, r=None):
        if q is None or r is None:
            q, r = self.dual(x)
        if (not np.isfinite(q).all() or not np.isfinite(r).all()
                or np.max(np.sqrt(np.sum(abs(r)**2, axis=0))) > 1+1e-14):
            raise ValueError('nonfinite or infeasible dual')
        primal = self.p.physical(self.I*x)[0]
        lower = float(-.5*np.vdot(q, q).real-np.vdot(self.b, q).real
                      -conjugate(-self.BH@q-self.TH(r), self.w))
        gap = checked_gap(primal, lower)
        return dict(primal=primal, lower=lower, gap=gap)


def checked_gap(primal, lower):
    gap = float(primal-lower)
    tol = 128*np.finfo(float).eps*max(1., abs(primal), abs(lower))
    if not np.isfinite(gap) or gap < -tol:
        raise ValueError('invalid numerical primal/dual gap')
    return gap


def steps(s):
    ab = abs(s.B)
    intensity = s.I.reshape(s.p.shape)/s.p.tvnorm
    rx, ry = intensity.copy(), intensity.copy()
    rx[:, :-1] += intensity[:, 1:]
    ry[:-1, :] += intensity[1:, :]
    col = 2*intensity.copy()
    col[:, 1:] += intensity[:, 1:]
    col[1:, :] += intensity[1:, :]
    # Zero data rows impose no restriction; choose a finite positive step.
    rows = ab.sum(axis=1)
    tau = THETA/(ab.sum(axis=0)+col.ravel())
    sd = THETA/np.where(rows > 0, rows, 1.)
    st = THETA/np.maximum(rx, ry)
    if not all(np.all(np.isfinite(t) & (t > 0)) for t in (tau, sd, st)):
        raise ValueError('invalid diagonal steps')
    return tau, sd, st


def solve_problem(problem, initialP, *, gap_atol, maxiter, check_every, time_limit, callback=None):
    """Return best feasible physical P, report, scalar-only checkpoint trace.

    Best primal and best feasible lower bound may come from different iterations;
    their pairing is valid. No full-image history or dual storage is returned.
    """
    if (not np.isfinite([gap_atol, time_limit]).all()
            or gap_atol < 0 or time_limit <= 0
            or type(maxiter) is not int or maxiter < 1
            or type(check_every) is not int or check_every < 1):
        raise ValueError('invalid numerical budget/tolerances')
    started = time.monotonic()
    s = Split(problem)
    x = np.asarray(initialP, complex).ravel()/s.I
    if x.shape != s.I.shape or not np.isfinite(x).all() or np.any(abs(x) > 1):
        raise ValueError('initial P must be feasible; no clipping')
    tau, sd, st = steps(s)
    q, r = s.dual(x)
    initial = s.cert(x, q, r)
    bestp, bestlb, bestx = initial['primal'], initial['lower'], x.copy()
    bestpi = bestli = 0
    trace = []
    xbar = x.copy()
    def checkpoint(k):
        nonlocal bestp, bestlb, bestx, bestpi, bestli
        current, algorithm = s.cert(x), s.cert(x, q, r)
        if current['primal'] < bestp:
            bestp, bestx, bestpi = current['primal'], x.copy(), k
        lower = max(current['lower'], algorithm['lower'])
        if lower > bestlb:
            bestlb, bestli = lower, k
        row = dict(iteration=k, seconds=time.monotonic()-started,
                   current=current, algorithm_dual=algorithm, best_primal=bestp,
                   best_lower=bestlb, best_gap=checked_gap(bestp, bestlb))
        trace.append(row)
        if callback is not None:
            callback(dict(row))
        return row['best_gap'] <= gap_atol
    certified = checkpoint(0)
    k = 0
    while not certified and k < maxiter and time.monotonic()-started < time_limit:
        k += 1
        q = (q+sd*(s.B@xbar-s.b))/(1+sd)
        r = project(r+st*s.T(xbar))
        xn = prox(x-tau*(s.BH@q+s.TH(r)), tau, s.w)
        xbar, x = 2*xn-x, xn
        if k % check_every == 0 or k == maxiter or time.monotonic()-started >= time_limit:
            certified = checkpoint(k)
    if trace[-1]['iteration'] != k:
        certified = checkpoint(k)
    # Re-evaluate the actual returned primal against the best historical LB.
    saved = problem.physical(s.I*bestx)[0]
    gap = checked_gap(saved, bestlb)
    report = dict(initial=initial, best_primal=saved, best_lower=bestlb, best_gap=gap,
                  best_primal_iteration=bestpi, best_lower_iteration=bestli,
                  iterations=k, seconds=time.monotonic()-started, theta=THETA,
                  gap_atol=gap_atol, gap_rtol=0., numerically_certified=certified,
                  termination='gap' if certified else ('maxiter' if k == maxiter else 'time_limit'),
                  certificate='floating-point primal/historical dual lower; not interval arithmetic')
    return s.I*bestx, report, trace

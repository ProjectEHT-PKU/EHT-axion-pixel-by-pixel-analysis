"""Two small numerical contracts; no production images or observations."""
import numpy as np
import pytest
from eht_axion.eht import convex as o
pytest_plugins = ['test_eht_calibration']  # shared tiny synthetic fixture
UNIT = dict(pvis=1., m=1., hw=100., ptv=1.)
PAPER = dict(p_data={'pvis': 1, 'm': 1}, p_reg={'hw': 100, 'ptv': 1}, norm_reg=True)


def test_exact_oracles_and_diagonal_certificate():
    rng = np.random.default_rng(202)
    shape = (2, 3)
    intensity = np.linspace(.8, 1.4, 6)
    A = rng.normal(size=(9, 6))+1j*rng.normal(size=(9, 6))
    data = {k: (rng.normal(size=9)+1j*rng.normal(size=9), np.linspace(.7, 1.3, 9), A)
            for k in ('pvis', 'm')}
    p = o.CartesianProblem(intensity, shape, data, 1.7, intensity.sum(), weights=UNIT)
    s = o.Split(p)
    known = np.arange(1, 7).reshape(shape).astype(complex)
    np.testing.assert_array_equal(o.D(known, shape),
        [[[1, 1, -3], [1, 1, -6]], [[3, 3, 3], [-4, -5, -6]]])
    x = .1*(rng.normal(size=6)+1j*rng.normal(size=6))
    r = rng.normal(size=(2, *shape))+1j*rng.normal(size=(2, *shape))
    q = rng.normal(size=18)+1j*rng.normal(size=18)
    np.testing.assert_allclose(np.vdot(s.T(x), r).real, np.vdot(x, s.TH(r)).real, atol=1e-13)
    np.testing.assert_allclose(np.vdot(s.B@x, q).real, np.vdot(x, s.BH@q).real, atol=1e-13)
    tau, sd, st = o.steps(s)
    T = np.column_stack([s.T(v).ravel() for v in np.eye(6, dtype=complex)])
    K = np.vstack([s.B, T])
    sigma = np.r_[sd, np.broadcast_to(st, (2, *shape)).ravel()]
    assert np.max(sigma*abs(K).sum(axis=1)) <= o.THETA*(1+1e-14)
    assert np.max(tau*abs(K).sum(axis=0)) <= o.THETA*(1+1e-14)
    assert np.linalg.norm(np.sqrt(sigma)[:, None]*K*np.sqrt(tau)[None, :], 2) < 1
    v = rng.normal(size=6)+1j*rng.normal(size=6)
    z = o.prox(v, tau, s.w)
    np.testing.assert_allclose(z+tau*s.w*np.arctanh(abs(z))*z/abs(z), v, atol=2e-12)
    assert np.all(o.prox(np.zeros(6, complex), tau, s.w) == 0)
    grad = s.w*np.arctanh(abs(x))*x/abs(x)
    assert abs(np.sum(s.w*o.h(abs(x)))+o.conjugate(grad, s.w)-np.vdot(x, grad).real) < 1e-12
    assert np.isfinite(o.conjugate(np.full(6, 1e5), s.w))
    own = (.5*np.linalg.norm(s.B@x-s.b)**2+np.sum(s.w*o.h(abs(x)))
           +np.sqrt(np.sum(abs(s.T(x))**2, axis=0)).sum())
    np.testing.assert_allclose(own, p.physical(intensity*x)[0], atol=1e-13)
    dx = rng.normal(size=6)+1j*rng.normal(size=6)
    eps = 1e-6
    f, g = p.physical(intensity*x)
    fd = (p.physical(intensity*x+eps*dx)[0]-p.physical(intensity*x-eps*dx)[0])/(2*eps)
    assert abs(fd-np.vdot(g, dx).real) < 1e-7
    assert s.cert(x, q, o.project(r))['gap'] >= 0
    with pytest.raises(ValueError):
        s.cert(x, q*np.nan, o.project(r))
    with pytest.raises(ValueError):
        s.cert(x, q, np.ones_like(r))
    with pytest.raises(ValueError):
        o.checked_gap(1., 1.001)
    # Known zero optimum: nonnegative quadratic, entropy and unsmoothed TV.
    zero = {k: (np.zeros(9, complex), np.ones(9), A) for k in data}
    p0 = o.CartesianProblem(intensity, shape, zero, 1.7, intensity.sum(), weights=UNIT)
    f0, g0 = p0.physical(np.zeros(6, complex))
    assert f0 == 0 and np.all(g0 == 0)
    assert o.Split(p0).cert(np.zeros(6, complex))['gap'] == 0
    P, report, trace = o.solve_problem(p0, intensity*x, gap_atol=1e-8, maxiter=6000, check_every=100,
                                       time_limit=300)
    assert report['numerically_certified'] and report['best_gap'] <= 1e-8
    assert np.all(abs(P) < intensity)
    assert all(b['best_lower'] >= a['best_lower'] for a, b in zip(trace, trace[1:]))
    assert report['best_gap'] == pytest.approx(p0.physical(P)[0]-report['best_lower'])


def test_real_132_tiny_objective_and_baseline(toy):
    import ehtim as eh
    from scipy.optimize import minimize
    from eht_axion.eht import calibration as core
    obs, fixed = toy
    obs = core.make_imaging_obs(obs)
    im = eh.imager.Imager(obs, fixed, prior_im=fixed, flux=float(fixed.total_flux()),
        pol='P', transform=['mcv'], data_term={'pvis': 1, 'm': 1},
        reg_term={'hw': 100, 'ptv': 1}, norm_reg=True, norm_init=False, ttype='direct')
    p = o.exact_cartesian_problem(fixed, obs, library=im, **PAPER)
    initial = fixed.qvec+1j*fixed.uvec
    np.testing.assert_allclose(p.physical(initial)[0], im.objfunc(im._xinit)+2+100*np.log(2), atol=2e-12)
    # Independent upstream analytic-gradient optimizer, only on this 2x2 toy.
    baseline = minimize(im.objfunc, im._xinit, jac=im.objgrad, method='L-BFGS-B',
                        options={'maxiter': 2000, 'ftol': 1e-14, 'gtol': 1e-10})
    P, report, trace = o.solve_problem(p, initial, gap_atol=1e-7, maxiter=20000, check_every=100,
                                       time_limit=300)
    assert report['numerically_certified']
    official = float(baseline.fun)+2+100*np.log(2)
    assert report['best_lower'] <= official+1e-9
    assert p.physical(P)[0] <= official+2e-7
    # Construction path has no optimization and uses the same cached objective.
    rebuilt = o.exact_cartesian_problem(fixed, obs, library=eh, **PAPER)
    np.testing.assert_allclose(rebuilt.physical(initial)[0], p.physical(initial)[0], atol=1e-13)
    im.dat_term_next = {'pvis': 2, 'm': 1}
    with pytest.raises(ValueError):
        o.exact_cartesian_problem(fixed, obs, library=im, **PAPER)
    with pytest.raises(ValueError, match='normalized'):
        o.exact_cartesian_problem(fixed, obs, library=eh, **dict(PAPER, norm_reg=False))


@pytest.mark.parametrize('weights', [dict(pvis=2.5, m=.3, hw=40., ptv=3.),
                                     dict(pvis=1., m=1., hw=100., ptv=1.)])
def test_general_weights_match_official_objective_and_certify(toy, weights):
    """The solver objective equals the ehtim objective + pvis + m + hw ln 2 for any weights."""
    import ehtim as eh
    from eht_axion.eht import calibration as core
    obs, fixed = toy
    obs = core.make_imaging_obs(obs)
    p_data, p_reg = {k: weights[k] for k in ('pvis', 'm')}, {k: weights[k] for k in ('hw', 'ptv')}
    im = eh.imager.Imager(obs, fixed, prior_im=fixed, flux=float(fixed.total_flux()),
        pol='P', transform=['mcv'], data_term=p_data, reg_term=p_reg, norm_reg=True,
        norm_init=False, ttype='direct')
    p = o.exact_cartesian_problem(fixed, obs, library=im, p_data=p_data, p_reg=p_reg, norm_reg=True)
    initial = fixed.qvec+1j*fixed.uvec
    constant = weights['pvis']+weights['m']+weights['hw']*np.log(2)
    np.testing.assert_allclose(p.physical(initial)[0], im.objfunc(im._xinit)+constant, atol=2e-12)
    P, report, _ = o.solve_problem(p, initial, gap_atol=1e-7, maxiter=20000, check_every=100,
                                   time_limit=300)
    assert report['numerically_certified'] and report['best_gap'] <= 1e-7

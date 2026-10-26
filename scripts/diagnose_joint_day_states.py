#!/usr/bin/env python3
"""Peak-point latent concentration and contrast-basis invariance diagnostics.

No new posterior fits. Values are decompositions of completed joint-state
mixture results at their saved maximum phase-likelihood nodes.
"""
import argparse
import inspect
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp
from eht_axion.pipeline import common
from eht_axion.research import joint_days
from eht_axion.research.joint_days import contrast_observation
from eht_axion.research.time_sampling import H2, H3


def state_weights(data, covariance, signals, background, coupling, phase, weights):
    chol = np.linalg.cholesky(covariance)
    residual = data.ravel()[None]-background.reshape(len(background), -1)-coupling*np.real(
        np.exp(1j*phase)*signals.reshape(len(signals), -1))
    white = solve_triangular(chol, residual.T, lower=True).T
    q = np.einsum('md,md->m', white, white)
    logparts = np.log(weights)-q/2
    normalized = np.exp(logparts-logsumexp(logparts))
    ll = logsumexp(logparts)-.5*(data.size*np.log(2*np.pi)+2*np.log(np.diag(chol)).sum())
    return normalized, q, float(ll)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'study', 'output-root']:
        p.add_argument('--'+name, required=True, type=Path)
    args = p.parse_args(argv)
    with common.verification_session():
        for root in [args.inputs, args.study]:
            common.verify_stage(root, recursive=False, allow_imported=False)
        used, results = [Path(__file__).resolve(), Path(inspect.getfile(joint_days))], []
        snapshots = args.output_root.parent/(args.output_root.name+'-source-receipts')
        snapshots.mkdir(exist_ok=False)
        for label, root in [('inputs', args.inputs), ('study', args.study)]:
            dest = snapshots/(label+'-receipt.json')
            dest.write_bytes(common.manifest_path(root).read_bytes())
            used.append(dest)
        for dim in ['1d', '2d']:
            files = [args.inputs/f'states-{dim}.npz', args.inputs/f'observations-{dim}.npz']
            common.require_declared_outputs(args.inputs, files)
            used += files
            with np.load(files[0]) as z:
                s = {k:z[k] for k in z.files}
            with np.load(files[1]) as z:
                cov = z['covariances']
            ai = np.flatnonzero(np.isclose(s['alphas'], .3))[0]
            for h in [H2, H3]:
                filename = f'{dim}-a0.300000-h{len(h)}-joint_background_evolving_mixture-observed.npz'
                path = args.study/filename
                common.require_declared_outputs(args.study, [path])
                used.append(path)
                with np.load(path) as z:
                    likelihood = z['phase_loglike']
                    i, j = np.unravel_index(np.argmax(likelihood), likelihood.shape)
                    coupling, phase, stored = float(z['couplings'][i]), float(z['phases'][j]), likelihood[i, j]
                data, c = contrast_observation(s['observed_daily_local_unwrapped'], cov, h)
                rows = [0, 2] if len(h) == 2 else [0, 1, 2]
                background = s['background_h3_principal'][:, rows]
                background = background-background.mean(axis=0)
                signals = np.einsum('kd,mdn->mkn', h, s['evolving_daily'][ai])
                weights = s['evolving_weights']
                r, q, ll = state_weights(data, c, signals, background, coupling, phase, weights)
                np.testing.assert_allclose(ll, stored, rtol=1e-12, atol=1e-9)
                # Use det(A)=1 alternative basis, without re-wrapping any edge.
                transform = np.tril(np.ones((len(h), len(h))))
                big = np.kron(transform, np.eye(data.shape[1]))
                rr, qq, lll = state_weights(transform@data, big@c@big.T,
                    np.einsum('ij,mjn->min', transform, signals),
                    np.einsum('ij,mjn->min', transform, background), coupling, phase, weights)
                np.testing.assert_allclose(q, qq, rtol=1e-8, atol=1e-6)
                np.testing.assert_allclose(r, rr, rtol=1e-7, atol=1e-8)
                np.testing.assert_allclose(ll, lll, rtol=1e-9, atol=1e-6)
                top = np.argsort(r)[-5:][::-1]
                results.append(dict(dimension=dim, contrasts=len(h), coupling=coupling, phase=phase,
                    posterior_state_ESS=float(1/(r@r)), max_weight=float(r.max()),
                    first_half_mass=float(r[:len(r)//2].sum()),
                    state_count=len(r), peak_log_likelihood=ll, minimum_component_chi2=float(q.min()),
                    top_states=[dict(index=int(k), weight=float(r[k]),
                        start_hours=float(s['window_starts_hours'][k])) for k in top],
                    basis_invariant=True, scope='one saved peak, not a posterior-wide ESS'))
        with common.stage_directory(args.output_root, 'research:joint-day-state-diagnostics',
                                    dict(alpha=.3, no_refits=True), used) as out:
            (out/'report.json').write_text(json.dumps(dict(results=results), indent=2, allow_nan=False)+'\n')
            print(json.dumps(results, indent=2), flush=True)


if __name__ == '__main__':
    main()

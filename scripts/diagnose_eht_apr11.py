#!/usr/bin/env python3
"""April11 diagnostics only: no imaging, selfcal, leakage solve, or inference.

Requires completed common stages. Checkpoint selection uses DECLARED filenames,
not directory globbing. Statistics use propagated diagonal errors, not independent
noise/DOF claims. norm_reg certification remains pending unless separately audited.
"""
import argparse
import csv
from importlib.metadata import version
import json
from pathlib import Path
import re

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.pipeline import common
from rebuild_eht_fiducial import pol_noise, pol_chi2


def baseline_labels(t1, t2):
    return np.array(['-'.join(sorted((str(a), str(b)))) for a, b in zip(t1, t2)])


def aggregate(residual, sigma, labels):
    """Marginal complex chi2=sum(|r/s|²)/(2N); no fitted-parameter subtraction."""
    residual, sigma, labels = map(np.asarray, (residual, sigma, labels))
    if residual.shape != sigma.shape or labels.shape != sigma.shape:
        raise ValueError('residual/sigma/labels shapes differ')
    if not np.isfinite(residual).all() or not np.all(np.isfinite(sigma) & (sigma > 0)):
        raise ValueError('finite residuals and positive sigmas required')
    rows = []
    for label in sorted(set(labels)):
        z = residual[labels == label] / sigma[labels == label]
        rows.append(dict(group=str(label), count=len(z), chi2=float(np.mean(abs(z)**2)/2),
                         real_chi2=float(np.mean(z.real**2)), imag_chi2=float(np.mean(z.imag**2)),
                         median_abs_sigma=float(np.median(abs(z))),
                         median_real_sigma=float(np.median(z.real)),
                         median_imag_sigma=float(np.median(z.imag))))
    return rows


def map_metrics(i, p, j, q, mask):
    """Symmetric mean-P/I weights; axial EVPA distance, no alignment/rotation."""
    i, p, j, q, mask = map(np.asarray, (i, p, j, q, mask))
    valid = mask & (i > 0) & (j > 0) & (abs(p) > 0) & (abs(q) > 0)
    if not valid.any():
        raise ValueError('empty polarized comparison mask')
    w = (abs(p[valid])/i[valid] + abs(q[valid])/j[valid])/2
    phase = np.angle(q[valid]*p[valid].conj())
    delta = phase/2
    return dict(pixels=int(valid.sum()),
                evpa_rms_deg=float(np.rad2deg(np.sqrt(np.average(delta**2, weights=w)))),
                cos2_delta=float(np.average(np.cos(phase), weights=w)),
                qu_relative_l2=float(np.linalg.norm((q-p)[valid])/np.linalg.norm(p[valid])),
                net_m_a=float(abs(p.sum())/i.sum()), net_m_b=float(abs(q.sum())/j.sum()))


def read_dterms(path):
    raw = json.loads(Path(path).read_text())
    raw = raw.get('dterms', raw)
    result = {}
    for site, pair in raw.items():
        a = np.asarray(pair, float)
        if a.shape != (2, 2) or not np.isfinite(a).all():
            raise ValueError(f'invalid D terms: {site}')
        result[site] = a[:, 0] + 1j*a[:, 1]
    return result


def artifacts(root):
    record = common.verify_stage(root, recursive=True, allow_imported=False)
    names = {r['path'] for r in record['outputs']}
    required = ['final-raw-IQU.fits'] + (['steps.json'] if 'steps.json' in names else [])
    rounds = sorted((int(m[1]), name) for name in names
                    if (m := re.fullmatch(r'round-(\d+)-raw-IQU\.fits', name)))
    if not rounds:
        raise ValueError('no declared round checkpoints')
    final_round = rounds[-1][0]
    candidates = ['final-dterms.json', f'round-{final_round:03d}-dterms.json',
                  f'round-{final_round:02d}-dterms.json', f'round-{final_round:03d}.json']
    dname = next((n for n in candidates if n in names), None)
    if dname is None:
        raise ValueError('no declared final D snapshot')
    required += [dname] + [n for _, n in rounds]
    common.require_declared_outputs(root, [root/n for n in required])
    return record, rounds, root/dname, [root/n for n in required]


def corrected(base, dterms):
    from ehtim.observing import obs_simulate
    obs = pol_noise(base, .005)  # one fresh inflation BEFORE Jones propagation
    if set(dterms) != set(map(str, obs.tarr['site'])):
        raise ValueError('D table must match station names exactly')
    for row in obs.tarr:
        row['dr'], row['dl'] = dterms[str(row['site'])]
    obs.frcal, obs.dcal = True, False
    obs.data = obs_simulate.apply_jones_inverse(obs, frcal=True, dcal=False, verbose=False)
    obs.dcal = True
    return obs


def fit_rows(obs, image, minutes):
    fields = ['t1', 't2', 'time', 'u', 'v', 'vis', 'qvis', 'uvis', 'sigma', 'qsigma', 'usigma']
    data = obs.unpack(fields, conj=True)
    vi, vq, vu, _ = image.sample_uv(np.column_stack((data['u'], data['v'])),
                                    polrep_obs='stokes', ttype='direct', verbose=False)
    p = data['qvis'] + 1j*data['uvis']
    ps = np.hypot(data['qsigma'], data['usigma'])
    m = p/data['vis']
    ms = np.sqrt(ps**2 + (data['sigma']*abs(m))**2)/abs(data['vis'])
    signals = dict(I=(data['vis']-vi, data['sigma']), Q=(data['qvis']-vq, data['qsigma']),
                   U=(data['uvis']-vu, data['usigma']), pvis=(p-vq-1j*vu, ps),
                   m=(m-(vq+1j*vu)/vi, ms))
    labels = dict(all=np.repeat('ALL', len(data)), baseline=baseline_labels(data['t1'], data['t2']),
                  time=np.array([str(int(t)) for t in np.floor(data['time']*60/minutes)]))
    rows = []
    for signal, (r, s) in signals.items():
        for kind, groups in labels.items():
            rows += [dict(signal=signal, grouping=kind, **v) for v in aggregate(r, s, groups)]
        for site in sorted(set(data['t1']) | set(data['t2'])):
            keep = (data['t1'] == site) | (data['t2'] == site)
            rows += [dict(signal=signal, grouping='station', **v)
                     for v in aggregate(r[keep], s[keep], np.repeat(str(site), keep.sum()))]
    upstream = pol_chi2(obs, image, 'direct')
    own = {v['signal']: v['chi2'] for v in rows if v['grouping'] == 'all'}
    for key in ('pvis', 'm'):
        np.testing.assert_allclose(own[key], upstream[key], rtol=1e-9, atol=1e-11)
    return rows, upstream


def vectors(image):
    return image.imvec, image.qvec + 1j*image.uvec


def physical(image):
    i, p = vectors(image)
    if not np.isfinite(i).all() or not np.isfinite(p).all() or np.any(i <= 0):
        raise ValueError('finite IQU with positive I required')
    return dict(max_p_over_i=float(np.max(abs(p)/i)), pixels_p_ge_i=int(np.sum(abs(p) >= i)))


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', nargs='+', type=Path, required=True)
    parser.add_argument('--base-input', type=Path, required=True, help='declared gaincal clean-pol UVFITS')
    parser.add_argument('--feed-file', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--i-cut', type=float, default=.1)
    parser.add_argument('--p-cut', type=float, default=.2)
    parser.add_argument('--time-bin-minutes', type=float, default=10)
    parser.add_argument('--tail-rounds', type=int, choices=[10, 20], default=20)
    parser.add_argument('--norm-reg-audit-pending', action='store_true', required=True,
                        help='explicitly acknowledge that these diagnostics do not certify norm_reg')
    args = parser.parse_args(argv)
    if not (0 < args.i_cut <= 1 and 0 < args.p_cut <= 1 and 0 < args.time_bin_minutes < float('inf')):
        parser.error('positive finite bins and thresholds in (0,1] required')
    common.reject_directory_symlinks(args.output_root)
    if args.output_root.exists():
        raise FileExistsError('new output root required')
    if version('ehtim') != '1.3.2':
        raise ValueError('ehtim 1.3.2 required')
    import ehtim as eh
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.utils import iers
    iers.conf.auto_download = False
    with common.verification_session():
        records = [artifacts(root) for root in args.runs]
        common.verify_stage(args.base_input.parent, recursive=True, allow_imported=False)
        common.require_declared_outputs(args.base_input.parent, [args.base_input])
        inputs = [Path(__file__).resolve(), Path(__file__).with_name('rebuild_eht_fiducial.py'),
                  Path(__file__).with_name('_bootstrap.py'), Path(common.__file__), args.feed_file,
                  common.manifest_path(args.base_input.parent), args.base_input]
        for root, (_, _, _, paths) in zip(args.runs, records):
            inputs += [common.manifest_path(root), *paths]
        inputs += sorted(Path(eh.__file__).resolve().parent.rglob('*.py'))
        inputs = list(common.source_closure(inputs))
        destination = args.output_root.resolve()
        if any(p == destination or p.is_relative_to(destination) for p in inputs):
            raise ValueError('output/source overlap')
        base = eh.obsdata.load_uvfits(str(args.base_input), polrep='stokes').flag_uvdist(uv_min=1e8)
        if base.mjd != 57854 or 'M87' not in base.source.upper() or not 220e9 < base.rf < 240e9:
            raise ValueError('only April11 M87 230GHz supported')
        table = {r['site']: r for r in json.loads(args.feed_file.read_text())['table']}
        for row in base.tarr:
            for key in ('fr_par', 'fr_elev', 'fr_off'):
                row[key] = table[str(row['site'])][key]
        base.tarr['dr'] = base.tarr['dl'] = 0j
        base.frcal, base.dcal = True, False
        raw = [eh.image.load_fits(str(root/'final-raw-IQU.fits')) for root in args.runs]
        def geometry(im):
            return im.xdim, im.ydim, im.psize, im.ra, im.dec, im.rf
        if any(geometry(im) != geometry(raw[0]) for im in raw):
            raise ValueError('common image grid/coordinates required; no implicit registration')
        for im in raw:
            physical(im)
            np.testing.assert_allclose(im.imvec, raw[0].imvec, rtol=1e-12, atol=1e-15)
        display = [im.blur_circ(19*eh.RADPERUAS, 19*eh.RADPERUAS) for im in raw]
        mask = np.logical_and.reduce([(im.imvec >= args.i_cut*im.imvec.max()) &
                    (abs(im.qvec+1j*im.uvec) >= args.p_cut*abs(im.qvec+1j*im.uvec).max()) for im in display])
        report = dict(settings=vars(args) | {'runs': list(map(str, args.runs)),
            'base_input': str(args.base_input),
            'feed_file': str(args.feed_file), 'output_root': str(args.output_root)}, ehtim=version('ehtim'),
            norm_reg_certified=False,
            mask='intersection across all final common19uas maps; diagnostic, not official',
            covariance='propagated diagonal approximation; conjugate rows and station groups NOT independent',
            normalization='sum(abs(residual/sigma)**2)/(2N), no fitted DOF subtraction', runs=[])
        with common.stage_directory(args.output_root, 'research:eht-apr11-diagnostics',
                                    report['settings'], inputs) as out:
            np.save(out/'common-mask.npy', mask)
            allrows = []
            fig, axes = plt.subplots(len(raw), 3, squeeze=False, figsize=(14, 4*len(raw)))
            for idx, (root, image, shown, record) in enumerate(zip(args.runs, raw, display, records)):
                manifest, rounds, dpath, _ = record
                d = read_dterms(dpath)
                rows, chi = fit_rows(corrected(base, d), image, args.time_bin_minutes)
                allrows += [dict(run=str(root), **r) for r in rows]
                tail = []
                for k, name in rounds:
                    if k < rounds[-1][0]-args.tail_rounds:
                        continue
                    checkpoint = eh.image.load_fits(str(root/name))
                    if geometry(checkpoint) != geometry(image):
                        raise ValueError('checkpoint grid changed')
                    np.testing.assert_allclose(checkpoint.imvec, image.imvec, rtol=1e-12, atol=1e-15)
                    tail.append(dict(round=k, physical=physical(checkpoint),
                        raw=map_metrics(*vectors(image), *vectors(checkpoint), np.ones_like(mask)),
                        display=map_metrics(*vectors(shown), *vectors(checkpoint.blur_circ(
                            19*eh.RADPERUAS, 19*eh.RADPERUAS)), mask)))
                report['runs'].append(dict(path=str(root), source_config=manifest['config'], chi2=chi,
                    physical=physical(image), tail_to_final=tail,
                    common19_to_first=map_metrics(*vectors(display[0]), *vectors(shown), mask)))
                plt.sca(axes[idx, 0])
                shown.display(plotp=True, pcut=args.i_cut, mcut=0, nvec=min(16, shown.xdim),
                              show=False, axis=axes[idx, 0], has_cbar=False, has_title=False)
                axes[idx, 0].set_title(root.name+'; additional19uas')
                delta = (shown.qvec+1j*shown.uvec)-(display[0].qvec+1j*display[0].uvec)
                extent = (np.array([shown.xdim, -shown.xdim, -shown.ydim, shown.ydim])
                          * shown.psize/eh.RADPERUAS/2)
                axes[idx, 1].imshow(abs(delta).reshape(shown.ydim, shown.xdim), extent=extent, origin='upper')
                axes[idx, 1].set(xlabel='East offset (uas; reversed)', ylabel='North offset (uas)',
                                 title='|Δ(Q+iU)| vs first (Jy/pixel)')
                for signal in ('Q', 'U', 'pvis', 'm'):
                    selected = [r for r in rows if r['grouping'] == 'time' and r['signal'] == signal]
                    axes[idx, 2].plot([int(r['group'])*args.time_bin_minutes/60 for r in selected],
                                      [r['chi2'] for r in selected], '.', label=signal)
                axes[idx, 2].set(xlabel='Time (hours)', ylabel='Marginal complex chi2')
                axes[idx, 2].legend()
            fig.tight_layout()
            fig.savefig(out/'diagnostics.png', dpi=150)
            plt.close(fig)
            with (out/'residual-groups.csv').open('w') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(allrows[0]))
                writer.writeheader()
                writer.writerows(allrows)
            dump(out/'diagnostics.json', report)
    return report


if __name__ == '__main__':
    main()

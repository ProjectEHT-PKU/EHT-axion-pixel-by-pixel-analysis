#!/usr/bin/env python3
"""Controlled April11 working-model P/D comparisons, NOT certified official imaging.

--plan only verifies/hashes and prints; no ehtim/numpy import or output creation.
Execution requires --audit-file: completed parent normalization/objective evidence.
The --plan audit_template supplies exact fields; only the parent may mark it
complete and working_model_approved (NOT Paper VII certification). Selection
excludes case/round budget, so all three controls share one audited objective.
No continuation is supported before the normalization audit is settled.
"""
import argparse
import hashlib
from importlib import metadata
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW, UV = 'official-I-pre-restoration.fits', 'official-I-gaincal-clean-pol.uvfits'
FROZEN = dict(schema=1, ehtim='1.3.2', npix=64, fov_uas=120., flux_jy=.6,
    p_data={'pvis': 1., 'm': 1.}, p_reg={'hw': 100., 'ptv': 1.}, maxit=100,
    stop=1e-4, fsys=.005, leakage_tol=.1, leakage_inverse=False, norm_init=False,
    chisq_transform=False, ttype='direct', initial_fraction=.05, restart_blur_uas=20.,
    solve_sites=['AZ', 'LM', 'PV'], fixed_sites=['AA', 'AP', 'JC', 'SM'])


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=ROOT/'configs/eht-apr11-validation.yaml')
    p.add_argument('--case', choices=['baseline', 'evpa45', 'restart'], required=True)
    p.add_argument('--rounds', type=int, default=100)
    p.add_argument('--norm-reg', action=argparse.BooleanOptionalAction, default=None)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--audit-file', type=Path)
    p.add_argument('--plan', action='store_true')
    return p


def common():
    import _bootstrap  # noqa: F401
    from eht_axion.pipeline import common as c
    return c


def feed_table(path):
    feed = json.loads(Path(path).read_text())
    expected = dict(AA=(1, 0, 0), AP=(1, 1, 0), AZ=(1, 1, 0),
                    LM=(1, -1, 0), PV=(1, -1, 0), SM=(1, -1, 45), JC=(1, 0, 0))
    rows = feed['table']
    table = {r['site']: tuple(r[k] for k in ('fr_par', 'fr_elev', 'fr_off')) for r in rows}
    if feed.get('status') != 'verified' or len(rows) != 7 or table != expected:
        raise ValueError('verified seven-site feed map required; fr_off remains DEGREES')
    return table


def verified_products(root, stage, names):
    c = common()
    record = c.verify_stage(root, recursive=True, allow_imported=False)
    if record['stage'] != stage:
        raise ValueError(f'expected source stage {stage}')
    paths = [root / name for name in names]
    c.require_declared_outputs(root, paths)
    return [c.manifest_path(root), *paths]


def build_plan(args):
    import yaml
    c = yaml.safe_load(args.config.read_text())
    if not isinstance(c, dict) or set(c) != set(FROZEN) | {'norm_reg', 'i_root', 'feed_file'}:
        raise ValueError('unknown/missing configuration keys')
    # JSON comparison distinguishes booleans from numbers; numeric int/float equivalent is OK.
    for k, v in FROZEN.items():
        if c[k] != v or (isinstance(c[k], bool) != isinstance(v, bool)):
            raise ValueError(f'frozen {k} requires {v!r}')
    if type(c['schema']) is not int or type(c['npix']) is not int or type(c['maxit']) is not int:
        raise ValueError('integer schema/npix/maxit required')
    if type(c['norm_reg']) is not bool or not 1 <= args.rounds <= 100:
        raise ValueError('boolean norm_reg and 1..100 rounds required')
    if args.norm_reg is not None:
        c['norm_reg'] = args.norm_reg
    for k in ('i_root', 'feed_file'):
        if not isinstance(c[k], str) or not Path(c[k]).is_absolute():
            raise ValueError(f'{k} must be absolute')
    feed_table(c['feed_file'])
    inputs = [args.config.resolve(), Path(__file__).resolve(),
              Path(__file__).with_name('rebuild_eht_fiducial.py'),
              Path(__file__).with_name('rebuild_eht_pol_warmstart.py'),
              ROOT/'scripts/_bootstrap.py', ROOT/'eht_axion/pipeline/common.py', Path(c['feed_file'])]
    inputs += verified_products(Path(c['i_root']), 'research:official-eht-I', [RAW, UV])
    inputs = list(common().source_closure(inputs))
    out = args.output_root.absolute()
    common().reject_directory_symlinks(out)
    if out.exists():
        raise FileExistsError('new output root required')
    if args.audit_file and (out == args.audit_file.resolve() or out in args.audit_file.resolve().parents):
        raise ValueError('audit/output overlap')
    source_root = Path(c['i_root']).resolve(strict=True)
    if out == source_root or out.is_relative_to(source_root) or source_root.is_relative_to(out):
        raise ValueError('output/input root overlap')
    if any(out == p or out in p.parents or p in out.parents for p in inputs):
        raise ValueError('output/input overlap')
    selection = hashlib.sha256(json.dumps(c, sort_keys=True).encode()).hexdigest()
    return dict(parameters=c, selection_sha256=selection, case=args.case, rounds=args.rounds,
        evpa_deg=45 if args.case == 'evpa45' else 0,
        restart_rounds=[2, 12, 22, 32, 42] if args.case == 'restart' else [],
        output_root=str(out), hashes={str(p): digest(p) for p in inputs},
        claim='approved current working model; NOT certified official; no appearance-based selection',
        noise='once-inflated fresh base for every P/D solve; propagated diagonal covariance approximation',
        termination='fixed round budget, no early stop; partial if rounds<100',
        audit_template=dict(complete=False, working_model_approved=False,
            selection_sha256=selection, hashes={str(p): digest(p) for p in inputs}),
        audit='parent-completed normalization/objective evidence; working-model approval NOT certification')


def projected_init(image, fixed, beam):
    import numpy as np
    init = image.blur_circ(beam, beam) if beam else image.copy()
    init.imvec = fixed.imvec.copy()
    power = np.hypot(init.qvec, init.uvec)
    scale = np.minimum(1., (1-1e-8)*init.imvec / np.maximum(power, np.finfo(float).tiny))
    init.qvec *= scale
    init.uvec *= scale
    return init, int(np.count_nonzero(scale < 1))


def objective(imgr, image):
    from ehtim.imager import make_initarr
    # Physical I,rho,phi,psi array, not the inverse-transformed optimizer vector.
    arr = make_initarr(image, imgr._embed_mask, norm_init=False, pol=True)
    chi = {k: float(v) for k, v in imgr.make_chisq_dict(arr).items()}
    reg = {k: float(v) for k, v in imgr.make_reg_dict(arr).items()}
    weighted = {k: imgr.dat_term_next[k]*(v-1) for k, v in chi.items()}
    weighted.update({k: imgr.reg_term_next[k]*v for k, v in reg.items()})
    return dict(chi2=chi, regularizers=reg, weighted_components=weighted, total=sum(weighted.values()))


def residuals(obs, image):
    """Independent direct Fourier residuals, no polchisq reduction; no conjugate duplication."""
    import numpy as np
    from ehtim.observing.obs_helpers import ftmatrix
    data = obs.switch_polrep('stokes').unpack(
        ['t1', 't2', 'u', 'v', 'vis', 'qvis', 'uvis', 'sigma', 'qsigma', 'usigma'], conj=False)
    a = ftmatrix(image.psize, image.xdim, image.ydim,
                 np.column_stack((data['u'], data['v'])), pulse=image.pulse)
    result = {}
    pairs = np.array(['-'.join(sorted((str(x), str(y)))) for x, y in zip(data['t1'], data['t2'])])
    for baseline in ['ALL', *sorted(set(pairs))]:
        mask = np.ones(len(data), bool) if baseline == 'ALL' else pairs == baseline
        row = dict(rows=int(mask.sum()))
        for name, vector, field, sigma in [('I', image.imvec, 'vis', 'sigma'),
                ('Q', image.qvec, 'qvis', 'qsigma'), ('U', image.uvec, 'uvis', 'usigma')]:
            r = (data[field] - a @ vector)[mask]
            s = data[sigma][mask]
            if not np.isfinite(r).all() or not np.all(np.isfinite(s) & (s > 0)):
                raise ValueError('invalid residual/sigma')
            row[name] = dict(rms_jy=float(np.sqrt(np.mean(abs(r)**2))),
                             reduced_complex_chi2=float(np.mean(abs(r/s)**2)/2))
        result[baseline] = row
    return result


def run(plan, eh, pol_cal, out):
    import numpy as np
    from rebuild_eht_fiducial import pol_noise, pol_chi2, optimizer_records
    c = plan['parameters']
    fixed = eh.image.load_fits(str(Path(c['i_root'])/RAW)).regrid_image(c['fov_uas']*eh.RADPERUAS, c['npix'])
    fixed.imvec *= c['flux_jy']/fixed.imvec.sum()
    if not np.all(np.isfinite(fixed.imvec) & (fixed.imvec > 0)):
        raise ValueError('strictly positive finite fixed I required')
    base = eh.obsdata.load_uvfits(str(Path(c['i_root'])/UV), polrep='stokes').flag_uvdist(uv_min=1e8)
    sites = set(base.data['t1']) | set(base.data['t2'])
    if base.mjd != 57854 or 'M87' not in base.source.upper() or not 220e9 < base.rf < 240e9:
        raise ValueError('only M87 April11 230GHz approved')
    if not set(c['solve_sites']) <= sites or not sites <= set(base.tarr['site']) or not set(base.tarr['site']) <= set(c['solve_sites']+c['fixed_sites']):
        raise ValueError('unexpected/missing stations')
    table = feed_table(c['feed_file'])
    for row in base.tarr:
        row['fr_par'], row['fr_elev'], row['fr_off'] = table[str(row['site'])]
    base.tarr['dr'] = base.tarr['dl'] = 0j
    base.frcal, base.dcal = True, False
    # Single injection only; every total-D solve receives a COPY of the same noisy base.
    noisy = pol_noise(base, c['fsys'])
    current, image = noisy.copy(), fixed.copy()
    phase = np.deg2rad(2*plan['evpa_deg'])
    p = c['initial_fraction']*fixed.imvec*np.exp(1j*phase)
    image.add_qu(p.real.reshape(c['npix'], c['npix']), p.imag.reshape(c['npix'], c['npix']))
    image.add_pol_image(np.zeros((c['npix'], c['npix'])), 'V')
    previous_d = np.zeros((len(base.tarr), 2), complex)
    for k in range(1, plan['rounds']+1):
        previous = image
        beam = c['restart_blur_uas']*eh.RADPERUAS if k in plan['restart_rounds'] else 0
        init, projected = projected_init(image, fixed, beam)
        imgr = eh.imager.Imager(current.copy(), init, init, c['flux_jy'], pol='P', transform=['mcv'],
            data_term=c['p_data'], reg_term=c['p_reg'], maxit=c['maxit'], stop=c['stop'],
            ttype=c['ttype'], norm_reg=c['norm_reg'], norm_init=c['norm_init'], chisq_transform=False)
        before = objective(imgr, init)
        with optimizer_records() as actual:
            image = imgr.make_image_P(show_updates=False)
            after = objective(imgr, image)
            chi_before = pol_chi2(current, image, 'direct')
            current = pol_cal.leakage_cal(noisy.copy(), image, sites=c['solve_sites'],
                leakage_tol=c['leakage_tol'], inverse=c['leakage_inverse'], ttype='direct', show_solution=False)
        np.testing.assert_allclose(image.imvec, fixed.imvec, rtol=1e-12, atol=1e-15)
        if not np.all(np.isfinite(image.qvec+1j*image.uvec)) or np.any(np.hypot(image.qvec, image.uvec) >= fixed.imvec):
            raise ValueError('nonphysical mcv result')
        if list(current.tarr['site']) != list(base.tarr['site']):
            raise ValueError('leakage solver changed station ordering')
        d = np.column_stack((current.tarr['dr'], current.tarr['dl']))
        if not np.isfinite(d).all() or np.any(d[np.isin(current.tarr['site'], c['fixed_sites'])] != 0):
            raise ValueError('invalid/fixed-site D terms')
        rel = float(np.linalg.norm(image.qvec+1j*image.uvec-previous.qvec-1j*previous.uvec)/max(np.linalg.norm(previous.qvec+1j*previous.uvec), 1e-30))
        row = dict(round=k, blur_restart=bool(beam), projected_pixels=projected,
            objective_initial=before, objective_final=after, optimizers=actual,
            chi2_before_D=chi_before, chi2_after_D=pol_chi2(current, image, 'direct'),
            relative_qu_change=rel, max_delta_d=float(np.max(abs(d-previous_d))),
            dterms={str(s): [[float(z.real), float(z.imag)] for z in pair] for s, pair in zip(current.tarr['site'], d)})
        dump(out/f'round-{k:03d}.json', row)
        if k % 10 == 0 or k == plan['rounds']:
            image.save_fits(str(out/f'round-{k:03d}-raw-IQU.fits'))
        print(json.dumps(row, separators=(',', ':'), allow_nan=False), flush=True)
        previous_d = d.copy()
    image.save_fits(str(out/'final-raw-IQU.fits'))
    dump(out/'final-datafit.json', dict(chi2=row['chi2_after_D'],
        raw_uncorrected_model_mismatch=residuals(noisy, image), calibrated=residuals(current, image),
        covariance='diagonal rescaling approximation; raw mismatch is NOT a forward-Jones fit'))
    dump(out/'termination.json', dict(rounds=plan['rounds'], reason='fixed_budget', scientifically_converged=False))


def main(argv=None):
    args = parser().parse_args(argv)
    with common().verification_session():
        plan = build_plan(args)
        if args.plan:
            print(json.dumps(plan, indent=2, allow_nan=False))
            return plan
        if args.audit_file is None:
            raise ValueError('complete objective audit required BEFORE any case launch (--audit-file)')
        audit_file = args.audit_file.resolve(strict=True)
        audit = json.loads(audit_file.read_text())
        if audit.get('complete') is not True or audit.get('working_model_approved') is not True or audit.get('selection_sha256') != plan['selection_sha256'] or audit.get('hashes') != plan['hashes']:
            raise ValueError('parent audit not approved for these exact source bytes/objective selection')
        plan['audit_evidence'] = dict(path=str(audit_file), sha256=digest(audit_file),
            approval_scope='current working model only, NOT certified Paper VII')
        if metadata.version('ehtim') != '1.3.2':
            raise ValueError('optional eht-rebuild environment with ehtim==1.3.2 required')
        import ehtim as eh
        from ehtim.calibrating import pol_cal
        from astropy.utils import iers
        iers.conf.auto_download = False
        code = sorted(Path(eh.__file__).resolve().parent.rglob('*.py'))
        plan['runtime_code_hashes'] = {str(p): digest(p) for p in code}
        inputs = [Path(p) for p in plan['hashes']] + [audit_file] + code
        with common().stage_directory(Path(plan['output_root']), 'research:eht-apr11-validation', plan, inputs) as out:
            if out is None:
                raise FileExistsError('existing stage forbidden')
            dump(out/'plan.json', plan)
            run(plan, eh, pol_cal, out)
    return plan


if __name__ == '__main__':
    main()

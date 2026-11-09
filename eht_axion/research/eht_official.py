"""Small Paper VII/Appendix H primitives, not a reconstruction driver.

Pinned to ehtim 1.3.2. Circular products and their separate per-real-component
sigmas remain native. Jones correction propagates diagonal marginal variances;
it does NOT retain the induced full covariance. No imaging objective weights,
production execution, thermal bootstrap, or D refits inside samples live here.
"""

from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
from importlib.metadata import version
import inspect
import json
from pathlib import Path
import re

import numpy as np


_PRODUCTS = ("rr", "ll", "rl", "lr")
_PRECORRECTED = ("AA", "AP", "SM", "JC")


def _library():
    if version("ehtim") != "1.3.2":
        raise ValueError("ehtim 1.3.2 required")
    import ehtim
    return ehtim


def _native(obs):
    if obs.polrep != "circ":
        raise ValueError("native circular Obsdata required; do not round-trip through Stokes")
    if not len(obs.data):
        raise ValueError("empty observation")
    for p in _PRODUCTS:
        if (not np.isfinite(obs.data[p + "vis"]).all()
                or not np.isfinite(obs.data[p + "sigma"]).all()
                or np.any(obs.data[p + "sigma"] <= 0)):
            raise ValueError("all retained rows require four finite products and positive sigmas")


def _metadata(obs, **updates):
    meta = deepcopy(getattr(obs, "official_metadata", {}))
    meta.update(updates)
    obs.official_metadata = meta
    return obs


def load_native(path, feed_file, *, scan_average):
    """Load released HOPS zbl-dtcal+selfcal data without a Stokes conversion.

    Stage flags are from the 2023-D01-01 release README (not the loader's
    ambiguous dcal default). This recognizes the release filename convention;
    cryptographic verification of input provenance belongs to the caller.
    Reject/flag an entire row if ANY native product/error is missing, before
    averaging. This conservative complete-record policy avoids invented hands
    and undefined Jones mixing. Averaging uses the official coherent API once;
    pass scan_average=False if the I pipeline will perform that operation.
    """
    eh = _library()
    path, feed_file = Path(path), Path(feed_file)
    if not re.fullmatch(r"SR2_M87_2017_\d{3}_(?:lo|hi)_hops_zbl-dtcal\+selfcal\.uvfits", path.name):
        raise ValueError("only the declared released HOPS zbl-dtcal+selfcal stage is supported")
    feed = json.loads(feed_file.read_text())
    if feed.get("status") != "verified" or "degrees" not in feed.get("units", {}).get("fr_off", ""):
        raise ValueError("verified feed table with fr_off in degrees required")
    table = {str(r["site"]): r for r in feed["table"]}
    obs = eh.obsdata.load_uvfits(str(path), polrep="circ", remove_nan=False)
    before = len(obs.data)
    keep = np.ones(before, dtype=bool)
    for p in _PRODUCTS:
        keep &= (np.isfinite(obs.data[p + "vis"]) & np.isfinite(obs.data[p + "sigma"])
                 & (obs.data[p + "sigma"] > 0))
    obs = obs.copy()
    obs.data = obs.data[keep].copy()
    _native(obs)
    for row in obs.tarr:
        site = str(row["site"])
        if site not in table:
            raise ValueError(f"missing feed coefficients for {site}")
        for key in ("fr_par", "fr_elev", "fr_off"):
            value = float(table[site][key])
            if not np.isfinite(value):
                raise ValueError("nonfinite feed coefficient")
            row[key] = value
    obs.tarr["dr"] = obs.tarr["dl"] = 0j
    obs.ampcal = obs.phasecal = obs.frcal = True
    obs.dcal = False  # AA/AP/SM/JC precorrected, AZ/LM/PV still unsolved.
    scan_source = "not_requested"
    if scan_average:
        scan_source = "UVFITS_NX"
        if obs.scans is None or not len(obs.scans):
            obs.add_scans()
            scan_source = "official_add_scans_default_gap"
        obs = obs.avg_coherent(0, scan_avg=True)
        _native(obs)
    meta = dict(input=str(path.resolve()), feed_file=str(feed_file.resolve()),
                calibration_source="2023-D01-01 release README: HOPS zbl-dtcal+selfcal",
                provenance_hash_verified_here=False, native_polrep="circ", fr_off_unit="degrees",
                precorrected_sites=list(_PRECORRECTED), remaining_d_sites=["AZ", "LM", "PV"],
                missing_product_policy="drop complete row if any product/error invalid",
                rows_loaded=before, rows_flagged=int((~keep).sum()), rows_output=len(obs.data),
                scan_average=bool(scan_average), scan_source=scan_source,
                flags=dict(ampcal=True, phasecal=True, frcal=True, dcal=False),
                fsys_applied=False, full_covariance_retained=False)
    _metadata(obs, **meta)
    return obs, deepcopy(meta)


def add_fsys(obs, fraction):
    """Add (fraction*abs((RR+LL)/2))**2 to RL/LR per-real variance ONCE.

    No parallel-hand inflation, visibility perturbation, or averaging here.
    An already augmented copy raises rather than silently accumulating noise.
    """
    _native(obs)
    if not np.isfinite(fraction) or fraction < 0:
        raise ValueError("finite nonnegative fsys required")
    if getattr(obs, "official_metadata", {}).get("fsys_applied", False):
        raise ValueError("fsys already applied")
    out = obs.copy()
    variance = (fraction * abs((obs.data["rrvis"] + obs.data["llvis"]) / 2)) ** 2
    for p in ("rl", "lr"):
        out.data[p + "sigma"] = np.sqrt(obs.data[p + "sigma"] ** 2 + variance)
    return _metadata(out, fsys_applied=True, fsys=float(fraction),
                     fsys_sigma_convention="per real component, independently in RL and LR")


@lru_cache(maxsize=1)
def _adapter_type():
    Obsdata = _library().obsdata.Obsdata

    class NativeImagingObs(Obsdata):
        """Local ratio fix and rectangular-tlist dtype compatibility."""

        def tlist(self, conj=False, t_gather=0., scan_gather=False):
            # Upstream np.array(groups, dtype=object) loses the structured
            # dtype if every time has the same number of baselines (NumPy 2).
            # Reconstitute only those groups; no change to grouping semantics.
            groups = super().tlist(conj=conj, t_gather=t_gather, scan_gather=scan_gather)
            return [g if g.dtype.names else np.array(list(g), dtype=self.poltype) for g in groups]

        def unpack_dat(self, data, fields, conj=False, ang_unit="deg", debias=False, timetype=False):
            out = super().unpack_dat(data, fields, conj=conj, ang_unit=ang_unit,
                                     debias=debias, timetype=timetype)
            requested = [fields] if isinstance(fields, str) else fields
            # Obsdata.unpack(conj=True) already passes data_conj(): LR* and its
            # native sigma become RL on the reversed baseline. Never reorder it.
            if self.polrep == "circ" and any(f in requested for f in ("msigma", "msnr", "msigma_phase")):
                intensity = (data["rrvis"] + data["llvis"]) / 2
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = data["rlvis"] / intensity
                    sigma = np.sqrt(data["rlsigma"] ** 2 + abs(ratio) ** 2
                                    * (data["rrsigma"] ** 2 + data["llsigma"] ** 2) / 4) / abs(intensity)
                    if "msigma" in requested:
                        out["msigma"] = sigma
                    if "msnr" in requested:
                        out["msnr"] = abs(ratio) / sigma
                    if "msigma_phase" in requested:
                        out["msigma_phase"] = sigma / abs(ratio) / (np.pi / 180 if ang_unit == "deg" else 1)
            return out

    return NativeImagingObs


def make_imaging_obs(circ):
    """Native-circular adapter for official Imager(pol='P')/pvis/m terms.

    Images must be Stokes; observations need not be. Official pvis uses RL
    and conj(LR) with distinct sigmas. Correct the 1.3.2 merr2 missing /4 in
    Var(I). This first-order diagonal approximation does not restore Jones-
    induced covariance or account for correlation between pvis and m terms.
    Re-wrap after library operations that construct a plain Obsdata.
    """
    _native(circ)
    out = circ.copy()
    out.__class__ = _adapter_type()
    return _metadata(out, ratio_error="Var(I)=(sigma_RR^2+sigma_LL^2)/4; diagonal delta method",
                     full_covariance_retained=False)


def _weight_copy(base):
    _native(base)
    if not (base.ampcal and base.phasecal and base.frcal):
        raise ValueError("gain- and field-rotation-calibrated base required")
    if base.dcal:
        raise ValueError("use the uncorrected remaining-D base, not an already corrected output")
    out = make_imaging_obs(base)
    original = {key: out.tarr[key].copy() for key in ("sefdr", "sefdl")}
    # apply_jones_inverse otherwise recomputes sigmas when all SEFDs > 0.
    # Set only the private working table to zero and restore on output.
    out.tarr["sefdr"] = out.tarr["sefdl"] = 0
    out.tarr["dr"] = out.tarr["dl"] = 0j
    return out, original


def _finish_correction(out, base, original, method):
    for key, values in original.items():
        out.tarr[key] = values
    out.dcal = True
    out.official_metadata = deepcopy(getattr(base, "official_metadata", {}))
    return _metadata(out, d_application=method,
                     jones_weight_policy="propagate input diagonal sigmas; private SEFDs zeroed/restored",
                     full_covariance_retained=False)


def _complex(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        if np.shape(value) != (2,):
            raise ValueError("complex values must be complex scalars or [real,imag]")
        value = complex(*value)
    value = complex(value)
    if not np.isfinite(value):
        raise ValueError("finite complex D required")
    return value


def _hands(value):
    if isinstance(value, dict):
        return _complex(value["R"]), _complex(value["L"])
    if len(value) != 2:
        raise ValueError("D entry requires R and L")
    return _complex(value[0]), _complex(value[1])


def apply_fixed_d(base, dterms):
    """Apply one full inverse Jones solution, no fit and no random noise.

    dterms={site: {'R': complex, 'L': complex}}; [R,L] and JSON
    [[ReR,ImR],[ReL,ImL]] entries also accepted. Unlisted sites have zero D.
    """
    _library()
    from ehtim.observing import obs_simulate
    out, original = _weight_copy(base)
    for site, values in dterms.items():
        if site not in out.tkey:
            raise ValueError(f"unknown station {site}")
        dr, dl = _hands(values)
        if abs(1 - dr * dl) < 1e-12:
            raise ValueError("singular leakage Jones matrix")
        out.tarr["dr"][out.tkey[site]] = dr
        out.tarr["dl"][out.tkey[site]] = dl
    out.data = obs_simulate.apply_jones_inverse(out, dcal=False, frcal=True, opacitycal=True, verbose=False)
    return _finish_correction(out, base, original, "fixed_D_once_no_refit")


# Exact installed 1.3.2 function source, including whitespace. Fail closed if
# upstream changes; never adapt a guessed/new solver based only on its version.
_FORWARD_SOURCE_SHA256 = "3e092df9fd6f7fca8d45adf23b9ba457765ba10b0bd7151da2c319f6cb15d500"
_FORWARD_EDITS = (
    ("        prior = np.sum((np.abs(Dpar)**2)/(leakage_tol**2))",
     "        prior = 0.0"),
    ("        priorgrad = 2*Dpar / (leakage_tol**2)",
     "        priorgrad = np.zeros_like(Dpar, dtype=np.float64)"),
)


def _adapt_forward_source(source):
    digest = sha256(source.encode("utf-8")).hexdigest()
    if digest != _FORWARD_SOURCE_SHA256:
        raise ValueError("unrecognized official forward solver source SHA256; refusing adaptation")
    adapted = source
    for old, new in _FORWARD_EDITS:
        if adapted.count(old) != 1:
            raise ValueError("official forward solver prior anchor must occur exactly once")
        adapted = adapted.replace(old, new, 1)
    return adapted


def _forward_solver():
    """Compile a private two-line adaptation; leave upstream globals untouched."""
    _library()
    from ehtim.calibrating import pol_cal_new
    source = inspect.getsource(pol_cal_new.leakage_cal_new)
    adapted = _adapt_forward_source(source)
    filename = str(Path(pol_cal_new.__file__).resolve()) + ":explicit-zero-prior"
    namespace = dict(vars(pol_cal_new))
    exec(compile(adapted, filename, "exec"), namespace)
    receipt = dict(ehtim="1.3.2", upstream_function="pol_cal_new.leakage_cal_new",
                   upstream_file=str(Path(pol_cal_new.__file__).resolve()),
                   official_source_sha256=sha256(source.encode("utf-8")).hexdigest(),
                   adapted_source_sha256=sha256(adapted.encode("utf-8")).hexdigest(),
                   changes=[dict(old=old, new=new) for old, new in _FORWARD_EDITS],
                   prior="explicit constant zero", prior_gradient="explicit zeros_like(Dpar)")
    return namespace["leakage_cal_new"], receipt


def forward_solver_identity():
    """JSON-serializable source-hash/anchor receipt for driver provenance.

    SHA256 covers inspect.getsource(leakage_cal_new), not the entire module.
    Both validation and compilation run afresh; no mutation/patch of ehtim.
    """
    return _forward_solver()[1]


def fit_forward_d(base, image, *, sites):
    """Official nonlinear forward RL/LR fit, fixed errors, V=0, NO D prior.

    Always solve against a fresh gain-calibrated base, not the preceding
    corrected output. Corrected marginals follow the official inverse Jones.
    Optimizer stopping options are official defaults; callers may record the
    optimizer via their existing instrumentation, not a second solver here.
    """
    leakage_cal_new, identity = _forward_solver()
    work, original = _weight_copy(base)
    sites = list(sites)
    present = set(work.data["t1"]) | set(work.data["t2"])
    if not sites or len(set(sites)) != len(sites) or not set(sites) <= present:
        raise ValueError("unique nonempty solve sites must occur in data")
    model = image.switch_polrep("stokes")
    if len(model.vvec) and np.any(model.vvec != 0):
        raise ValueError("forward D model requires V=0")
    model.add_pol_image(np.zeros((model.ydim, model.xdim)), "V")
    out = leakage_cal_new(work, model, sites=sites, leakage_tol=None,
                          rescale_leakage_tol=False, pol_fit=["RL", "LR"],
                          ttype="direct", use_grad=True, show_solution=False, apply_solution=True)
    _finish_correction(out, base, original, "official_full_forward_RL_LR")
    return _metadata(out, d_prior="explicit zero prior and gradient", d_solve_sites=sites,
                     forward_solver_identity=identity,
                     forward_errors="fixed native input marginal sigmas")


def d_distribution(dailyD, *, ddof):
    """Per-site/hand mean and 2x2 (Re,Im) covariance, with source IDs.

    dailyD={source_id: {site: {'R': complex, 'L': complex}}}. Require four
    distinct days to avoid silently calling another estimator Appendix H.
    Values are fractional complex D, NOT percentages. ddof is an explicit
    implementation convention (not uniquely specified by the paper).
    """
    if len(dailyD) != 4:
        raise ValueError("four daily solutions required")
    if isinstance(ddof, bool) or not isinstance(ddof, (int, np.integer)) or not 0 <= ddof < 4:
        raise ValueError("ddof must be an integer in [0,3]")
    ids = sorted(dailyD)
    sites = sorted(dailyD[ids[0]])
    if not sites or any(set(dailyD[i]) != set(sites) for i in ids):
        raise ValueError("all daily solutions must contain the same nonempty station set")
    entries = {}
    for site in sites:
        values = np.array([_hands(dailyD[i][site]) for i in ids])
        entries[site] = {}
        for j, hand in enumerate(("R", "L")):
            xy = np.column_stack((values[:, j].real, values[:, j].imag))
            entries[site][hand] = dict(mean=xy.mean(axis=0).tolist(),
                                       covariance=np.cov(xy, rowvar=False, ddof=ddof).tolist())
    return dict(source_ids=ids, ddof=int(ddof), units="fractional D", entries=entries,
                independence="between stations, hands and separate calls/days; no temporal covariance")


def sample_dterms(distribution, sites, n, seed, *, residual_component_sigma):
    """Deterministic Generator draws, no thermal resampling or D refitting.

    AA/AP/SM/JC use ZERO-mean residuals with the given per-component SD (complex
    RMS sqrt(2) times larger), not their previously applied full D. Re/Im covariance is
    preserved for other stations. For different-day independent ensembles use
    distinct seeds and record them; reusing a seed intentionally repeats draws.
    """
    sites = sorted(sites)
    if (not sites or len(set(sites)) != len(sites) or isinstance(n, bool)
            or not isinstance(n, (int, np.integer)) or n < 1):
        raise ValueError("unique nonempty sites and positive integer n required")
    if not np.isfinite(residual_component_sigma) or residual_component_sigma < 0:
        raise ValueError("nonnegative finite residual component SD required")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("explicit nonnegative integer seed required")
    rng = np.random.default_rng(seed)
    draws = [{} for _ in range(n)]
    used = {}
    for site in sites:
        used[site] = {}
        for hand in ("R", "L"):
            if site in _PRECORRECTED:
                mean, covariance = np.zeros(2), np.eye(2) * residual_component_sigma ** 2
            else:
                entry = distribution["entries"][site][hand]
                mean, covariance = np.asarray(entry["mean"]), np.asarray(entry["covariance"])
            if (mean.shape != (2,) or covariance.shape != (2, 2)
                    or not np.isfinite(mean).all() or not np.isfinite(covariance).all()
                    or not np.allclose(covariance, covariance.T, rtol=0, atol=1e-14)
                    or np.linalg.eigvalsh(covariance).min() < -1e-14):
                raise ValueError("finite mean and positive-semidefinite symmetric 2x2 covariance required")
            xy = rng.multivariate_normal(mean, covariance, size=n, check_valid="raise")
            used[site][hand] = dict(mean=mean.tolist(), covariance=covariance.tolist())
            for k in range(n):
                draws[k].setdefault(site, {})[hand] = xy[k].tolist()
    return dict(seed=int(seed), bit_generator=type(rng.bit_generator).__name__, n=int(n), sites=sites,
                source_ids=list(distribution["source_ids"]), ddof=distribution["ddof"],
                residual_component_sigma=float(residual_component_sigma),
                residual_sigma_convention="per real/imaginary component; zero mean",
                independence=distribution["independence"], distributions=used, draws=draws,
                thermal_resampled=False, refit_d_after_draw=False)

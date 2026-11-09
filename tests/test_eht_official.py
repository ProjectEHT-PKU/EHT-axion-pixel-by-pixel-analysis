"""Tiny synthetic API/noise/Jones oracles only; no released-data reconstruction."""

import json
import numpy as np
import pytest

from eht_axion.research import eht_official as core


@pytest.fixture
def toy():
    eh = pytest.importorskip("ehtim")
    from astropy.utils import iers
    iers.conf.auto_download = False
    from ehtim import const_def as c
    sites = ["AZ", "LM", "PV"]
    tarr = np.zeros(3, dtype=c.DTARR)
    tarr["site"] = sites
    # Geocentric coordinates: avoid degenerate feed rotation geometry.
    xyz = [[-1828796, -5054407, 3427865], [-768716, -5988542, 2063276], [5088968, -301682, 3825016]]
    for k, (x, y, z) in enumerate(xyz):
        tarr[k]["x"], tarr[k]["y"], tarr[k]["z"] = x, y, z
    tarr["fr_par"] = 1
    tarr["fr_elev"] = [1, -1, -1]
    tarr["fr_off"] = [0, 15, 45]
    tarr["sefdr"], tarr["sefdl"] = 1000, 2000  # must NOT replace input errors
    rows = []
    rng = np.random.default_rng(71)
    for time in np.linspace(0, 10, 12):
        for a, b in ((0, 1), (0, 2), (1, 2)):
            row = np.zeros((), dtype=c.DTPOL_CIRC)
            row["time"], row["tint"] = time, 10
            row["t1"], row["t2"] = sites[a], sites[b]
            row["u"], row["v"] = rng.uniform(1e8, 2e9, 2)
            for p, sigma in zip(core._PRODUCTS, (.011, .023, .037, .059)):
                row[p + "sigma"] = sigma
            rows.append(row)
    obs = eh.obsdata.Obsdata(12.5, 12.4, 230e9, 2e9, np.array(rows), tarr,
                            polrep="circ", ampcal=True, phasecal=True, frcal=True, dcal=False,
                            mjd=57854, scantable=np.array([[-.1, 10.1]]))
    image = eh.image.Image(np.array([[.07, .13], [.23, .17]]), 15 * eh.RADPERUAS, 12.5, 12.4,
                           rf=230e9, mjd=57854)
    image.add_qu(np.array([[.005, .015], [-.01, .008]]), np.array([[.006, -.007], [.012, .009]]))
    image.add_pol_image(np.zeros((2, 2)), "V")
    simulated = image.observe_same_nonoise(obs, ttype="direct", zero_empty_pol=True, verbose=False)
    simulated.dcal = False
    return simulated, image


def angles(obs):
    fields = obs.unpack(["el1", "el2", "par_ang1", "par_ang2"], ang_unit="rad")
    result = []
    for leg in (1, 2):
        table = obs.tarr[[obs.tkey[s] for s in obs.data[f"t{leg}"]]]
        result.append(table["fr_elev"] * fields[f"el{leg}"] + table["fr_par"] * fields[f"par_ang{leg}"]
                      + np.deg2rad(table["fr_off"]))
    return result


def matrices(obs, dterms):
    phi1, phi2 = angles(obs)
    result = []
    for leg, phi in ((1, phi1), (2, phi2)):
        values = []
        for site, angle in zip(obs.data[f"t{leg}"], phi):
            dr, dl = core._hands(dterms.get(site, [0j, 0j]))
            values.append([[1, dr * np.exp(2j * angle)], [dl * np.exp(-2j * angle), 1]])
        result.append(np.array(values))
    return result


def coherency(obs):
    return np.array([[obs.data["rrvis"], obs.data["rlvis"]],
                     [obs.data["lrvis"], obs.data["llvis"]]]).transpose(2, 0, 1)


def set_coherency(obs, value):
    for p, i, j in (("rr", 0, 0), ("rl", 0, 1), ("lr", 1, 0), ("ll", 1, 1)):
        obs.data[p + "vis"] = value[:, i, j]


def corrupt(obs, dterms):
    j1, j2 = matrices(obs, dterms)
    out = obs.copy()
    set_coherency(out, j1 @ coherency(obs) @ j2.conj().transpose(0, 2, 1))
    return out


def test_fsys_exact_once_and_copies(toy):
    obs, _ = toy
    before = obs.data.copy()
    out = core.add_fsys(obs, .005)
    added = (.005 * abs((before["rrvis"] + before["llvis"]) / 2)) ** 2
    for p in core._PRODUCTS:
        np.testing.assert_array_equal(out.data[p + "vis"], before[p + "vis"])
        np.testing.assert_allclose(out.data[p + "sigma"] ** 2,
                                   before[p + "sigma"] ** 2 + (added if p in ("rl", "lr") else 0))
    np.testing.assert_array_equal(obs.data, before)
    with pytest.raises(ValueError, match="already"):
        core.add_fsys(out.copy(), .005)
    assert not np.shares_memory(out.data, obs.data)


def test_adapter_real_library_pvis_ratio_and_imager(toy):
    eh = pytest.importorskip("ehtim")
    from ehtim.imaging import pol_imager_utils as pol
    from ehtim.observing import obs_helpers
    obs, image = toy
    original_merr = obs_helpers.merr2
    adapted = core.make_imaging_obs(obs)
    dat = obs.data_conj()
    pvis, psigma, _ = pol.chisqdata_pvis(adapted, image, [])
    np.testing.assert_array_equal(pvis, dat["rlvis"])
    np.testing.assert_array_equal(psigma, dat["rlsigma"])
    assert len(np.unique(psigma)) == 2  # no RL/LR symmetrization
    m, sigma, _ = pol.chisqdata_m(adapted, image, [])
    intensity = (dat["rrvis"] + dat["llvis"]) / 2
    expected = np.sqrt(dat["rlsigma"] ** 2 + abs(m) ** 2
                       * (dat["rrsigma"] ** 2 + dat["llsigma"] ** 2) / 4) / abs(intensity)
    np.testing.assert_allclose(sigma, expected)
    assert np.any(abs(obs.unpack("msigma", conj=True)["msigma"] - expected) > 1e-8)
    assert obs_helpers.merr2 is original_merr
    np.testing.assert_array_equal(adapted.copy().unpack("msigma")["msigma"],
                                   adapted.unpack("msigma")["msigma"])
    imgr = eh.imager.Imager(adapted, image, image, .6, pol="P", transform=["mcv"],
                           data_term={"pvis": 1, "m": 1}, reg_term={"hw": 1},
                           ttype="direct", norm_reg=True, norm_init=False, maxit=1)
    result = imgr.make_image_P(show_updates=False)
    np.testing.assert_array_equal(result.imvec, image.imvec)


def test_fixed_d_matrix_values_and_noise_oracle(toy):
    obs, _ = toy
    dterms = {"AZ": [.04 + .07j, -.03 + .02j], "LM": [-.02 + .04j, .01 - .05j],
              "PV": [-.13 - .01j, .15 + .006j]}
    measured = corrupt(obs, dterms)
    saved = measured.copy()
    result = core.apply_fixed_d(measured, dterms)
    np.testing.assert_allclose(coherency(result), coherency(obs), rtol=1e-11, atol=1e-13)
    j1, j2 = matrices(measured, dterms)
    a, b = np.linalg.inv(j1), np.linalg.inv(j2)
    # Independent vec(C) linear-map marginal oracle, with unequal native errors.
    for k in range(len(measured.data)):
        input_sigmas = [measured.data[p + "sigma"][k] for p in ("rr", "rl", "lr", "ll")]
        linear = np.kron(a[k], b[k].conj())
        covariance = linear @ np.diag(np.square(input_sigmas)) @ linear.conj().T
        actual = [result.data[p + "sigma"][k] for p in ("rr", "rl", "lr", "ll")]
        np.testing.assert_allclose(np.square(actual), covariance.diagonal().real, rtol=1e-12)
    for p in core._PRODUCTS:
        np.testing.assert_array_equal(measured.data[p + "sigma"], saved.data[p + "sigma"])
    np.testing.assert_array_equal(result.tarr["sefdr"], measured.tarr["sefdr"])
    assert np.all(result.tarr["sefdr"] > 0)
    assert result.dcal and not measured.dcal
    np.testing.assert_array_equal(measured.tarr, saved.tarr)
    assert result.official_metadata["full_covariance_retained"] is False
    with pytest.raises(ValueError, match="already corrected"):
        core.apply_fixed_d(result, dterms)
    zero = core.apply_fixed_d(obs, {})
    for p in core._PRODUCTS:
        np.testing.assert_allclose(zero.data[p + "sigma"], obs.data[p + "sigma"])
    np.testing.assert_allclose(coherency(zero), coherency(obs))


def test_forward_objective_gradient_and_injected_recovery(toy, monkeypatch):
    from ehtim.calibrating import pol_cal_new
    obs, image = toy
    true_d = {"AZ": [.035 + .06j, -.04 + .025j], "LM": [-.02 + .035j, .01 - .04j],
              "PV": [-.12 - .015j, .14 + .008j]}
    measured = corrupt(obs, true_d)
    original_minimize = pol_cal_new.opt.minimize
    original_solver = pol_cal_new.leakage_cal_new
    records = []

    def capture(fun, x0, **kwargs):
        rng = np.random.default_rng(91)
        x = rng.normal(0, .04, len(x0))
        d = {s: x[4*k:4*k+4].copy().view(complex) for k, s in enumerate(("AZ", "LM", "PV"))}
        predicted = corrupt(obs, d)
        oracle = sum(np.sum(abs(measured.data[p + "vis"] - predicted.data[p + "vis"]) ** 2
                                   / measured.data[p + "sigma"] ** 2)
                     for p in ("rl", "lr")) / (4 * len(obs.data))
        np.testing.assert_allclose(fun(x), oracle, rtol=1e-12)
        numeric = []
        for j in range(len(x)):
            step = np.zeros(len(x))
            step[j] = 1e-6
            numeric.append((fun(x + step) - fun(x - step)) / 2e-6)
        np.testing.assert_allclose(kwargs["jac"](x), numeric, rtol=2e-7, atol=2e-8)
        # Exercise the actual adapted closures: the prior does not inspect
        # Dpar at all, even if a line-search trial contains inf/nan. Nonfinite
        # forward-model arithmetic is outside this prior-only guarantee.
        closure = dict(zip(fun.__code__.co_freevars, [c.cell_contents for c in fun.__closure__]))
        grad = kwargs["jac"]
        assert "leakage_tol" not in closure and not closure["rescale_leakage_tol"]
        for bad in (np.inf, np.nan):
            sentinel = np.full_like(x0, bad)
            # Replace only closure cell contents, not any module function.
            cell = fun.__closure__[fun.__code__.co_freevars.index("chisq_total")]
            grad_cell = grad.__closure__[grad.__code__.co_freevars.index("chisq_total_grad")]
            saved, saved_grad = cell.cell_contents, grad_cell.cell_contents
            try:
                cell.cell_contents = lambda _: 7.0
                grad_cell.cell_contents = lambda _: np.full_like(x0, 3.0)
                with np.errstate(all="raise"):
                    assert fun(sentinel) == 7.0
                    np.testing.assert_array_equal(grad(sentinel), np.full_like(x0, 3.0))
            finally:
                cell.cell_contents, grad_cell.cell_contents = saved, saved_grad
        # This is the actual official optimizer, not a replacement solver.
        res = original_minimize(fun, x0, **kwargs)
        records.append((fun(x0), fun(res.x)))
        return res

    monkeypatch.setattr(pol_cal_new.opt, "minimize", capture)
    before = measured.copy()
    fitted = core.fit_forward_d(measured, image, sites=("AZ", "LM", "PV"), optimizer_options=None)
    assert pol_cal_new.leakage_cal_new is original_solver
    assert fitted.official_metadata["d_prior"] == "explicit zero prior and gradient"
    assert fitted.official_metadata["forward_solver_identity"] == core.forward_solver_identity()
    assert records[0][1] < records[0][0] * 1e-8
    for site, target in true_d.items():
        row = fitted.tarr[fitted.tkey[site]]
        np.testing.assert_allclose([row["dr"], row["dl"]], target, atol=2e-5)
    for p in core._PRODUCTS:
        np.testing.assert_array_equal(measured.data[p + "sigma"], before.data[p + "sigma"])
    np.testing.assert_array_equal(measured.tarr, before.tarr)
    np.testing.assert_allclose(coherency(fitted), coherency(obs), atol=2e-5)


def test_forward_source_receipt_and_fail_closed():
    pytest.importorskip("ehtim")
    from hashlib import sha256
    import inspect
    from ehtim.calibrating import pol_cal_new
    source = inspect.getsource(pol_cal_new.leakage_cal_new)
    adapted = core._adapt_forward_source(source)
    identity = core.forward_solver_identity()
    assert identity["official_source_sha256"] == sha256(source.encode()).hexdigest()
    assert identity["adapted_source_sha256"] == sha256(adapted.encode()).hexdigest()
    assert identity["official_source_sha256"] != identity["adapted_source_sha256"]
    assert len(identity["changes"]) == 2
    restored = adapted
    for change in identity["changes"]:
        assert source.count(change["old"]) == 1
        restored = restored.replace(change["new"], change["old"], 1)
    assert restored == source
    json.dumps(identity, allow_nan=False)
    with pytest.raises(ValueError, match="SHA256"):
        core._adapt_forward_source(source + "\n")


def test_forward_zero_leakage_recovery(toy):
    obs, image = toy
    fitted = core.fit_forward_d(obs, image, sites=("AZ", "LM", "PV"), optimizer_options=None)
    np.testing.assert_allclose(fitted.tarr["dr"], 0, atol=1e-10)
    np.testing.assert_allclose(fitted.tarr["dl"], 0, atol=1e-10)
    np.testing.assert_allclose(coherency(fitted), coherency(obs), atol=1e-12)


def test_load_native_flags_feed_and_official_scan_api(toy, monkeypatch, tmp_path):
    eh = pytest.importorskip("ehtim")
    obs, _ = toy
    raw = obs.copy()
    raw.data["rlsigma"][0] = np.nan
    raw.data["lrvis"][1] = np.nan
    raw.data["rrsigma"][2] = 0
    raw.dcal = True  # loader default must NOT be trusted
    table = dict(status="verified", units={"fr_off": "degrees"},
                 table=[{key: (str(row[key]) if key == "site" else float(row[key]))
                         for key in ("site", "fr_par", "fr_elev", "fr_off")} for row in raw.tarr])
    feed = tmp_path / "feed.json"
    feed.write_text(json.dumps(table))
    calls = []

    def loader(path, **kwargs):
        calls.append(kwargs)
        return raw.copy()

    monkeypatch.setattr(eh.obsdata, "load_uvfits", loader)
    path = tmp_path / "SR2_M87_2017_101_lo_hops_zbl-dtcal+selfcal.uvfits"
    loaded, meta = core.load_native(path, feed, scan_average=False)
    assert calls == [dict(polrep="circ", remove_nan=False)]
    assert meta["rows_flagged"] == 3 and len(loaded.data) == len(raw.data) - 3
    assert loaded.polrep == "circ" and loaded.ampcal and loaded.phasecal and loaded.frcal
    assert not loaded.dcal
    np.testing.assert_array_equal(loaded.tarr["fr_off"], raw.tarr["fr_off"])
    np.testing.assert_array_equal(loaded.tarr["dr"], 0)
    averaged, receipt = core.load_native(path, feed, scan_average=True)
    assert receipt["scan_average"] is True
    expected = loaded.avg_coherent(0, scan_avg=True)
    for p in core._PRODUCTS:
        np.testing.assert_array_equal(averaged.data[p + "sigma"], expected.data[p + "sigma"])
    assert averaged.polrep == "circ"
    with pytest.raises(ValueError, match="declared"):
        core.load_native(tmp_path / "unknown.uvfits", feed, scan_average=True)


def daily_solutions():
    return {f"day-{k}": {"AZ": {"R": complex(.01*k, .02*k), "L": complex(-.02*k, .01*k)},
                         "PV": {"R": complex(-.1 + .01*k, .02 - .01*k),
                                "L": complex(.14 - .005*k, .003*k)}} for k in range(4)}


def test_distributions_covariance_and_sampling():
    daily = daily_solutions()
    distribution = core.d_distribution(daily, ddof=1)
    xy = np.array([[.01*k, .02*k] for k in range(4)])
    np.testing.assert_allclose(distribution["entries"]["AZ"]["R"]["mean"], xy.mean(axis=0))
    np.testing.assert_allclose(distribution["entries"]["AZ"]["R"]["covariance"],
                               np.cov(xy, rowvar=False, ddof=1))
    population = core.d_distribution(daily, ddof=0)
    np.testing.assert_allclose(population["entries"]["AZ"]["R"]["covariance"],
                               np.array(distribution["entries"]["AZ"]["R"]["covariance"]) * .75)
    sites = ["AZ", "PV", "AA", "AP", "SM", "JC"]
    state = np.random.get_state()
    samples = core.sample_dterms(distribution, sites, 2000, 2021, residual_component_sigma=.01)
    assert samples == core.sample_dterms(distribution, reversed(sites), 2000, 2021,
                                         residual_component_sigma=.01)
    assert samples != core.sample_dterms(distribution, sites, 2000, 2022, residual_component_sigma=.01)
    np.testing.assert_array_equal(np.random.get_state()[1], state[1])
    az = np.array([d["AZ"]["R"] for d in samples["draws"]])
    # Rank-one covariance gives a stringent direct test of within-D correlation.
    np.testing.assert_allclose(az[:, 1], 2 * az[:, 0], atol=1e-8)
    aa = np.array([d["AA"]["R"] for d in samples["draws"]])
    np.testing.assert_allclose(aa.mean(axis=0), 0, atol=.0007)
    np.testing.assert_allclose(aa.std(axis=0), .01, rtol=.07)
    for site in core._PRECORRECTED:
        assert samples["distributions"][site]["R"]["mean"] == [0, 0]
        assert samples["distributions"][site]["L"]["covariance"] == [[.0001, 0], [0, .0001]]
    assert samples["source_ids"] == sorted(daily)
    assert not samples["thermal_resampled"] and not samples["refit_d_after_draw"]


@pytest.mark.parametrize("ddof", [-1, 4, 1.5, True])
def test_distribution_bad_ddof(ddof):
    with pytest.raises(ValueError, match="ddof"):
        core.d_distribution(daily_solutions(), ddof=ddof)


def test_validation(toy):
    obs, image = toy
    with pytest.raises(ValueError, match="circular"):
        core.add_fsys(obs.switch_polrep("stokes"), .005)
    invalid = obs.copy()
    invalid.data["lrsigma"][0] = np.inf
    with pytest.raises(ValueError, match="finite"):
        core.make_imaging_obs(invalid)
    with pytest.raises(ValueError, match="unknown station"):
        core.apply_fixed_d(obs, {"OTHER": [0j, 0j]})
    with pytest.raises(ValueError, match="singular"):
        core.apply_fixed_d(obs, {"AZ": [1, 1]})
    model = image.copy()
    model.add_pol_image(np.ones((2, 2)), "V")
    with pytest.raises(ValueError, match="V=0"):
        core.fit_forward_d(obs, model, sites=("AZ", "LM", "PV"), optimizer_options=None)
    with pytest.raises(ValueError, match="four daily"):
        core.d_distribution({"day": {}}, ddof=1)

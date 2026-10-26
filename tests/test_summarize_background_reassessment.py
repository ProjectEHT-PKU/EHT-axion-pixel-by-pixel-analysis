"""Synthetic postprocessing tests; do not run production summary or science."""

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "background_summary", SCRIPTS / "summarize_background_reassessment.py"
)
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


def toy_rows():
    import itertools

    rows = []
    for d, a, h, b, response, mode in itertools.product(
        ("1d", "2d"), s.ALPHAS, (2, 3), s.MODELS, s.RESPONSES, s.MODES
    ):
        rows.append(
            dict(
                key=f"{d}-a{a:.6f}-h{h}-{b}-{response}-{mode}",
                dimension=d,
                alpha=a,
                contrasts=h,
                background=b,
                response=response,
                mode=mode,
                c95=2.0 + a + h,
                state_count=141 if response == "frozen_mixture" else 1,
                fixed_cov_condition=2.0,
            )
        )
    return rows


def toy_cv():
    return [
        dict(
            dimension_label=d,
            contrasts=h,
            n_splits=n,
            fold=f,
            model=b,
            n_test=24 + f,
            n_train=65,
            dimension=h * {"1d": 23, "2d": 20}[d],
            scoreable=True,
            mean_nll_per_dimension=-2.0 + 0.1 * f - 0.05 * s.MODELS.index(b) if b in s.MODELS else 100.0 + f,
            condition_number=1e6 if b == "full" else 2.0,
        )
        for d in ("1d", "2d")
        for h in (2, 3)
        for n in (5, 3)
        for f in range(n)
        for b in s.MODELS + ("diagonal", "full")
    ]


def test_all_512_keys_and_missing_duplicate_rejected():
    rows = toy_rows()
    assert len(s.index_rows(rows, s.ALPHAS)) == 512
    with pytest.raises(ValueError, match="duplicate"):
        s.index_rows(rows + rows[:1], s.ALPHAS)
    with pytest.raises(ValueError, match="counts"):
        s.index_rows(rows[:-1], s.ALPHAS)
    with pytest.raises(ValueError, match="filename"):
        s.index_rows([dict(rows[0], key="../escape")])
    with pytest.raises(ValueError, match="state count"):
        s.index_rows([dict(rows[0], state_count=3)])


def test_posterior_normal_monotone_sha_and_axes(tmp_path):
    c = np.geomspace(1 / 137, 1000, 120)
    ll = -c
    density, cdf = s.normalized_posterior_cdf(c, ll, prior="log-uniform")
    limit = s.posterior_upper_limit(c, ll, prior="log-uniform")
    fields = dict(
        couplings=c, loglike=ll, density=density, cdf=cdf, c95=limit, phase_loglike=np.zeros((120, 64))
    )
    row = dict(key="toy", c95=limit)
    meta = dict(coupling_points=120, phase_count=64, prior_bounds=[1 / 137, 1000])
    p = tmp_path / "posterior.npz"
    np.savez(p, **fields)
    result = s.validate_posterior(p, row, meta)
    assert result["sha256"] == s.sha256(p)
    assert result["phase_axis_checked"]
    assert result["normalized_mass"] == pytest.approx(1)
    np.savez(p, **dict(fields, cdf=np.zeros(120)))
    with pytest.raises(ValueError, match="CDF"):
        s.validate_posterior(p, row, meta)
    np.savez(p, **dict(fields, density=density[:, None]))
    with pytest.raises(ValueError, match="axis"):
        s.validate_posterior(p, row, meta)
    np.savez(p, **dict(fields, phase_loglike=np.zeros((64, 120))))
    with pytest.raises(ValueError, match="phase axis"):
        s.validate_posterior(p, row, meta)
    np.savez(p, **fields)
    with pytest.raises(ValueError, match="c95 mismatch"):
        s.validate_posterior(p, dict(row, c95=2 * limit), meta)
    del fields["phase_loglike"]
    np.savez(p, **fields)
    assert not s.validate_posterior(p, row, meta, require_phase=False)["phase_axis_checked"]


def test_no_partial_fold_average_and_weighted_score():
    records = toy_cv()
    aggregates = s.aggregate_cv(records)
    assert len(records) == 192
    assert len(aggregates) == 48
    item = next(
        x
        for x in aggregates
        if (x["dimension"], x["contrasts"], x["n_splits"], x["model"]) == ("1d", 2, 5, "iso")
    )
    assert item["weighted_nll"] == pytest.approx(
        np.average([-2 + 0.1 * f for f in range(5)], weights=range(24, 29))
    )
    for r in records:
        if (r["dimension_label"], r["contrasts"], r["model"], r["fold"]) == ("1d", 3, "full", 0):
            r["scoreable"] = False
            r["mean_nll_per_dimension"] = None
    bad = [
        r
        for r in s.aggregate_cv(records)
        if (r["dimension"], r["contrasts"], r["model"]) == ("1d", 3, "full")
    ]
    assert all(r["weighted_nll"] is None and r["delta_to_iso"] is None for r in bad)
    assert sorted(r["scoreable_folds"] for r in bad) == [2, 4]
    with pytest.raises(ValueError, match="CV fold"):
        s.aggregate_cv(records[:-1])


def test_all_four_plots_include_singular_and_negative_nll(tmp_path):
    records = toy_cv()
    for r in records:
        if r["model"] == "full" and r["contrasts"] == 3:
            r.update(scoreable=False, mean_nll_per_dimension=None, condition_number=None)
    s.make_plots(toy_rows(), records, s.aggregate_cv(records), tmp_path)
    assert {p.name for p in tmp_path.glob("*.png")} == set(s.PLOTS)
    assert all((tmp_path / p).stat().st_size > 1000 for p in s.PLOTS)


def test_output_collision_rejected_without_touching_inputs(tmp_path):
    from argparse import Namespace

    with pytest.raises(FileExistsError, match="new noncanonical"):
        s.run(Namespace(output_root=tmp_path, inputs=tmp_path / "missing"))


def test_receipt_names_do_not_expand_manifest_ancestry(tmp_path):
    from eht_axion.pipeline.common import describe_inputs

    receipt = tmp_path / "source-receipt-00.json"
    receipt.write_text('{"complete":true,"inputs":[{"path":"/nonexistent/raw"}]}')
    records = describe_inputs([receipt])
    assert len(records) == 1
    assert records[0]["sha256"] == s.sha256(receipt)
    assert Path(records[0]["path"]).name == "source-receipt-00.json"

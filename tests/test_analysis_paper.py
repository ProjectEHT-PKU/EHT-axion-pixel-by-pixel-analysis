"""scripts/analysis/paper.py: number formatting, tables, figures, one end-to-end run on synthetic data."""

import csv
import json
from pathlib import Path
import re

import numpy as np
import pytest

from eht_axion.analysis import paper_figures as figures
from scripts.analysis import paper
from eht_axion import provenance
from eht_axion.experiment import read_experiment

MODES = ("std", "cor", "cor_fixed_background")


def test_formatting_helpers():
    assert paper.word(10) == "Ten" and paper.word(44) == "FortyFour" and paper.word(5) == "Five"
    assert paper.alpha_word(0.35) == "AlphaThirtyFive"
    assert paper.sig(0.5, 3) == "0.500" and paper.sig(132.3, 3) == "132" and paper.sig(0.0992, 3) == "0.0992"
    assert paper.sig(104.0, 2) == "100" and paper.approx(104.0) == "104" and paper.approx(0.2916) == "0.29"
    assert paper.sci(7.39e7, 2) == r"7.4\times10^{7}" and paper.sci(9.96e5, 2) == r"1.0\times10^{6}"
    assert paper.order(11579.0) == "10^{4}"
    assert paper.percent_change(np.log10(1.5)) == pytest.approx(50.0)
    assert paper.angle_runs([74, 84, 94, 104, 144, 154, 334], 10) == (
        r"$74^\circ$--$104^\circ$, $144^\circ$--$154^\circ$ and $334^\circ$"
    )
    assert paper.angle_runs([4, 14], 10) == r"$4^\circ$--$14^\circ$"


def test_numbers_reject_bad_names_and_write_sources():
    numbers = paper.Numbers()
    numbers.add("NumExample", "1.5", 1.5, "file.json field", "an example")
    with pytest.raises(ValueError):
        numbers.add("NumExample", "1", 1, "s", "duplicate")
    with pytest.raises(ValueError):
        numbers.add("Num2D", "1", 1, "s", "digits are not allowed in macro names")
    text = numbers.tex()
    assert "% an example; source: file.json field\n\\newcommand{\\NumExample}{1.5}" in text
    assert numbers.json()["numbers"]["NumExample"]["value"] == 1.5


def test_window_constants(tmp_path):
    source = tmp_path / "axion.c"
    source.write_text("double w;\n    window = exp(-pow(r / 25.0, 15.0));\n")
    assert paper.window_constants(source) == (25.0, 15.0)
    source.write_text("window = exp(-r);\n")
    with pytest.raises(ValueError):
        paper.window_constants(source)


def test_two_emitter_formula_matches_direct_evaluation():
    for ratio, chi12, phi12 in paper.TWO_EMITTER["check_points"]:
        k = figures.kappa(ratio, chi12, phi12)
        assert figures.kappa_direct(ratio, chi12, phi12, 1e-5, 256) == pytest.approx(k, rel=1e-6)
    assert figures.kappa(0.5, 0.3, 0.0) == pytest.approx(1.0)


def test_tables_are_complete_tabular_environments():
    alphas = [0.1, 0.3]
    results = {
        f"observed-{m}-{d}": dict(alphas=alphas, c95=[123.4, 0.0991]) for d in ("1d", "2d") for m in MODES
    }
    table = paper.table_limits(results, {"1d": 24, "2d": 20})
    assert table.startswith(r"\begin{tabular}{c ccc ccc}") and table.rstrip().endswith(r"\end{tabular}")
    assert (
        r"\multicolumn{3}{c}{24 sectors}" in table and "0.30 & 0.0991 &" in table and "0.10 & 123 &" in table
    )
    rows = {(d, m): (48 if d == "1d" else 40, 30.04, 25.66) for d in ("1d", "2d") for m in MODES}
    chi2 = paper.table_chi2(rows)
    assert "Sectors & 48 & 30.0 & 30.0 & 30.0 & 25.7 & 25.7 & 25.7 \\\\" in chi2
    assert chi2.rstrip().endswith(r"\end{tabular}")


# ----------------------------------------------------------------------------- synthetic data root
def stage(path, name, files):
    """Publish a completed stage with the given {relative path: writer(path)}."""
    with provenance.stage_directory(path, name, {"synthetic": name}, []) as out:
        for relative, writer in files.items():
            target = out / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            writer(target)


def text(value):
    return lambda path: path.write_text(value if isinstance(value, str) else json.dumps(value))


def array(value):
    return lambda path: np.save(path, value)


def arrays(**values):
    return lambda path: np.savez(path, **values)


def stokes(rng, shape):
    n = shape[-2]
    x, y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n), indexing="ij")
    ring = np.exp(-(((np.hypot(x, y) - 0.5) / 0.2) ** 2)) + 0.01
    data = np.empty(shape + (4,))
    data[..., 0] = ring
    angle = np.arctan2(y, x) + 0.1 * rng.standard_normal(shape)
    data[..., 1] = 0.3 * ring * np.cos(2 * angle)
    data[..., 2] = 0.3 * ring * np.sin(2 * angle)
    data[..., 3] = 0
    return data


@pytest.fixture
def data_root(tmp_path):
    rng = np.random.default_rng(1)
    config = json.loads(json.dumps(read_experiment()))
    alphas = config["shared"]["alphas"]
    root = tmp_path / "root"
    n = config["shared"]["image"]["shape"]
    selected = {
        "1d": np.array([74, 84, 94, 144, 154]),
        "2d": np.array([5 * 64 + 10, 30 * 64 + 30, 30 * 64 + 31]),
    }
    limits, results = {}, {}
    for d, idx in selected.items():
        k = len(idx)
        base = 0.01 * (np.eye(k) + 0.3)
        coordinates = (
            np.column_stack((idx, np.zeros(k)))
            if d == "1d"
            else np.column_stack((idx // 64, idx % 64)) * 1.875
        )
        limits[f"observations-{d}.npz"] = arrays(
            selected_indices=idx,
            coordinates=coordinates,
            pair_covariances=np.stack((base, 2 * base)),
            pair_differences=0.05 * rng.standard_normal((2, k)),
        )
        limits[f"observations-{d}.json"] = text({"pair_condition_numbers": [7.4e7, 1.1e6]})
        differences = 0.03 * rng.standard_normal((138, k))
        limits[f"background-selected-{d}.npz"] = arrays(
            selected_indices=idx,
            differences=differences,
            sigma_radians=differences.std(),
            mean_radians=differences.mean(),
        )
        limits[f"background-selected-{d}.json"] = text({"lag_hours": 26.68})
        for m in MODES:
            fits = [
                dict(
                    data_quadratic=30.0 + j,
                    joint_grid_chi2=25.0 + j,
                    prediction_radians=np.zeros((2, k)).tolist(),
                )
                for j in range(len(alphas))
            ]
            limits[f"observed-{m}-{d}.json"] = text(
                dict(
                    alphas=alphas,
                    fits=fits,
                    selected_indices=idx.tolist(),
                    effective_pair_covariances=np.stack((base, 2 * base)).tolist(),
                    branch_risk=[
                        {"distance_to_seam_in_sigma": [9.0, 5.5]},
                        {"distance_to_seam_in_sigma": [7.0]},
                    ],
                )
            )
            for data in ("observed", "zero-data"):
                c95 = list(np.geomspace(100, 0.1, len(alphas)) * (1.5 if data == "observed" else 1))
                results[f"{data}-{m}-{d}"] = dict(alphas=alphas, c95=c95, log10_c95=list(np.log10(c95)))
    limits["report.json"] = text(dict(scientific_profile=config, results=results))
    stage(root / "analysis/limits", "analysis:limits", limits)
    shifted = {
        key: dict(value, log10_c95=list(np.asarray(value["log10_c95"]) + 1e-3))
        for key, value in results.items()
    }
    stage(
        root / "analysis/check-grid",
        "analysis:check-grid",
        {"report.json": text(dict(results=shifted, max_delta_dex=1e-3))},
    )
    stage(
        root / "analysis/check-prior",
        "analysis:check-grid",
        {"report.json": text(dict(results=results, max_delta_dex=0.0))},
    )

    def covariance_csv(path):
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=["half", "mode", "dimension", "alpha", "c95", "delta_dex"]
            )
            writer.writeheader()
            for d in ("1d", "2d"):
                for m in MODES:
                    writer.writerow(
                        dict(half="first500", mode=m, dimension=d, alpha=0.1, c95=1.0, delta_dex=-0.05)
                    )

    stage(
        root / "analysis/check-covariance",
        "analysis:check-covariance",
        {"covariance-sensitivity.csv": covariance_csv, "report.json": text({})},
    )
    templates = root / "analysis/templates"
    stage(
        templates / "processed/alpha-0.30",
        "analysis:preprocess",
        {"stokes.npy": array(stokes(rng, (1, 1, n, n)))},
    )
    stage(
        templates / "response/alpha-0.30",
        "analysis:fit",
        {
            "response-2d.npz": arrays(
                amplitude=0.01 + 0.01 * rng.random((1, n, n)),
                phase=rng.uniform(-np.pi, np.pi, (1, n, n)),
                snapshot_labels=np.array([2300]),
            )
        },
    )
    for d, shape in (("1d", (360,)), ("2d", (n, n))):
        stage(
            templates / f"templates/{d}",
            "analysis:templates",
            {
                "templates.npz": arrays(
                    alphas=np.asarray(alphas),
                    amplitude_mean=0.04 + 0.01 * rng.random((len(alphas),) + shape),
                    amplitude_std=0.005 * rng.random((len(alphas),) + shape) + 1e-4,
                )
            },
        )
    stage(templates / "background", "analysis:background", {"stokes.npy": array(stokes(rng, (6, 1, n, n)))})
    days = {}
    for day in (95, 96, 100, 101):
        values = np.column_stack(
            (
                np.full(10, 0.6),
                np.full(10, 0.02),
                np.full(10, 0.065),
                np.zeros(10),
                np.full(10, 0.03),
                np.full(10, 0.01),
            )
        )
        days[f"day-{day}/global.npz"] = arrays(
            values=values,
            columns=np.array(["total_I", "m_net", "m_resolved", "evpa_rad", "beta2_real", "beta2_imag"]),
        )
    coordinates = (np.arange(n) - (n - 1) / 2) * 120 / n
    stage(
        root / "eht/ensemble",
        "eht:ensemble",
        dict(
            days,
            **{
                "axes.npz": arrays(east_uas=-coordinates, north_uas=coordinates),
                "report.json": text(dict(bunit="JY/PIXEL", frequency_hz=227.0707e9)),
            },
        ),
    )
    shard = root / "raptor/delivered/data/raptor/res-400/background/t1600-t3000"
    stage(
        shard,
        "raptor:delivery",
        {
            "stokes-t1600-t3000.npy": array(np.ones((141, 1, 2, 2, 4)) * 1e9),
            "generated.json": text({"times_tg": list(np.arange(15010.0, 15151.0))}),
        },
    )
    stage(
        root / "raptor/runs/main/build",
        "raptor:build",
        {"source/src/axion.c": text("    window = exp(-pow(r / 25.0, 15.0));\n")},
    )
    tables = {}
    for key, spin, members in (
        ("low", 0.9375, [a for a in alphas if a <= 0.3]),
        ("high", 0.99, [a for a in alphas if a > 0.3]),
    ):
        r = np.linspace(1, 240, 200)
        rows = [(a, x, np.exp(-(((x - 20 / a) / (10 / a)) ** 2))) for a in members for x in r]
        path = tmp_path / f"radial-{key}.csv"
        np.savetxt(path, rows, delimiter=",")
        tables[key] = dict(path=str(path), spin=spin, sha256=provenance.sha256(path))
    (root / "raptor/plan.json").write_text(json.dumps({"radial_tables": tables}))
    return root


def test_run_writes_one_complete_paper_directory(data_root, monkeypatch):
    monkeypatch.delenv("EHT_RADIAL_LOW", raising=False)
    monkeypatch.delenv("EHT_RADIAL_HIGH", raising=False)
    paper.main(["--data-root", str(data_root), "--no-env-file"])
    out = data_root / "analysis/paper"
    manifest = provenance.verify_stage(out)
    assert manifest["stage"] == "analysis:paper"
    assert sorted(p.name for p in (out / "figures").iterdir()) == sorted(f + ".pdf" for f in figures.FIGURES)
    tex = (out / "numbers.tex").read_text()
    macros = re.findall(r"\\newcommand\{\\(Num[A-Za-z]+)\}\{(.*)\}", tex)
    values = json.loads((out / "numbers.json").read_text())["numbers"]
    assert [name for name, _ in macros] == list(values)
    assert all(values[name]["tex"] == value and values[name]["source"] for name, value in macros)
    lookup = dict(macros)
    assert lookup["NumSectors"] == "5" and lookup["NumPixelComponents"] == "6"
    assert lookup["NumDtermCorrClosePixelsMin"] == "0.2"
    assert lookup["NumSectorAngles"] == r"$74^\circ$--$94^\circ$ and $144^\circ$--$154^\circ$"
    assert lookup["NumMinBranchDistanceSigma"] == "5.5"
    assert lookup["NumConditionSectorsPairOne"] == r"7.4\times10^{7}"
    assert lookup["NumObsFrequencyGHz"] == "227.07" and lookup["NumReconPolMaxPercent"] == "6.5"
    assert lookup["NumRatioCorBgSectorsMin"] == "1.5" and lookup["NumChiGainMax"] == "5.0"
    assert lookup["NumCovCheckCorBgMaxPercentApprox"] == "12" and lookup["NumWindowAtTwentyRg"] == "0.97"
    assert lookup["NumSnapshots"] == "141" and lookup["NumSpanTg"] == "140"
    assert lookup["NumCloudBAlphaThirty"] == "0.69"
    assert (out / "table-limits.tex").read_text().startswith(r"\begin{tabular}")
    assert r"\multicolumn{3}{c}{5 sectors}" in (out / "table-limits.tex").read_text()
    with pytest.raises(FileExistsError):
        paper.main(["--data-root", str(data_root), "--no-env-file"])


def test_products_from_another_experiment_are_refused(data_root, tmp_path):
    config = read_experiment()
    config["analysis"]["inference"]["credibility"] = 0.9
    import yaml

    other = tmp_path / "other.yaml"
    other.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="different experiment"):
        paper.main(
            [
                "--data-root",
                str(data_root),
                "--no-env-file",
                "--config",
                str(other),
                "--output-root",
                str(tmp_path / "new"),
            ]
        )


def test_radial_table_must_match_the_plan(data_root, tmp_path):
    wrong = tmp_path / "wrong.csv"
    wrong.write_text("0.1,1,1\n")
    with pytest.raises(ValueError, match="differs from the one recorded"):
        paper.main(
            [
                "--data-root",
                str(data_root),
                "--no-env-file",
                "--radial-low",
                str(wrong),
                "--output-root",
                str(tmp_path / "new"),
            ]
        )
    assert not (tmp_path / "new").exists()
    assert Path(data_root / "raptor/plan.json").exists()

"""Synthetic wiring regressions: real workflow/EHT/inference, no science data.

Only aperture quadrature, posterior evaluation and generated-delivery identity
are substituted. Stage transactions, lineage checks and compact products remain
real; the delivery substitute admits a tiny footer fixture, not production data.
"""

from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from test_pipeline_end_to_end import make_inputs
from test_pipeline_eht import tiny_operator
from eht_axion import polarization
from eht_axion.pipeline import eht, inference, raptor, workflow
from eht_axion.pipeline.common import verify_stage
from eht_axion.pipeline.config import read_config
from eht_axion.producer import delivery

ROOT = Path(__file__).resolve().parents[1]


def independent_scales(config):
    """SI expressions, independent of analysis_geometry/PhysicalScale helpers."""
    physical = config["physical"]
    gm = 1.3271244e20 * physical["mass_solar"]
    angular_rg = (gm / 299792458.0**2
                  / (physical["distance_mpc"] * 1e6 * 3.085677581491367e16)
                  * 180 / np.pi * 3.6e9)
    return angular_rg, gm / 299792458.0**3 / 3600


@pytest.mark.parametrize("mass_factor,distance_factor,fwhm,fraction,footer_step", [
    (1.0, 1.0, 19.503299817566703, 1.0, 2.5),
    (1.25, 1.0, 19.503299817566703, 1.0, 2.5),
    (1.0, 0.8, 19.503299817566703, 1.0, 2.5),
    (1.0, 1.0, 12.0, 0.6, 5.0),
])
def test_input_parameters_reach_real_workflow_and_inference(
    tmp_path, monkeypatch, mass_factor, distance_factor, fwhm, fraction, footer_step,
):
    config = read_config(ROOT / "configs/analysis.yaml")
    # Non-default template normalization is an analysis-only contract: the
    # producer's first-scope validator deliberately fixes its own norm to 0.1.
    del config["producer"]
    config["physical"]["mass_solar"] *= mass_factor
    config["physical"]["distance_mpc"] *= distance_factor
    config["data"]["alphas"] = [0.1, 0.3]
    config["image"].update(final_shape=8, gaussian_fwhm_uas=fwhm, rotation_degrees=37.0)
    config["aperture"]["radial_max_fraction"] = fraction
    config["eht"].update(topset_size=6, min_radius_rg=0.0, historical_mask=False,
                         covariance_ridge=0.001)
    config["eht_input"].update(
        shape=8, pixel_size_uas=16.0, angle_roll=3, exclude_columns=1,
        exclusion_rectangles=[[2, 3, 2, 4]],
        rank_offsets={str(day): [0, 0] for day in eht.DAYS},
        final_offsets={str(day): [0, 0] for day in eht.DAYS},
    )
    config["background"].update(source_resolution=100, lag_frames=3)
    config["cloud"].update(coefficient=0.047, pivot_alpha=0.12, exponent=1.3, reference_norm=0.27)
    config["inference"].update(phase_count=4, coupling_points=5,
                               exposure_hours=2.5, day_gap_hours=23.0, pair_separation_hours=97.0)
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(config))
    assert read_config(path) == config
    before = deepcopy(config)

    raw, out = tmp_path / "input", tmp_path / "output"
    make_inputs(raw)
    labels = list(range(1600, 1690, 10))
    for alpha in config["data"]["alphas"]:
        # Separate directories avoid sharing footer metadata across alpha files.
        original = raw / "raptor-stokes/res-100" / f"alpha-{alpha}.npy"
        directory = raw / "raptor/res-100" / f"alpha-{alpha:.2f}"
        directory.mkdir(parents=True)
        original.rename(directory / "stokes.npy")
        (directory / "generated.json").write_text(json.dumps({
            "snapshot_labels": labels,
            "times_tg": [15010 + footer_step * i for i in range(len(labels))],
        }))
    contract = config["eht_input"]
    evidence = raw / "eht-input-geometry.json"
    evidence.write_text(json.dumps({
        "input_geometry": {key: contract[key] for key in ("shape", "fov_uas", "pixel_size_uas")},
        "fits_samples": [],
    }))

    def tiny_generated_inputs(paths, profile, *, alpha, resolution, **kwargs):
        # Do not mock cadence: workflow and inference must read actual footer
        # times. **kwargs also permits the shared optional verified/seen cache.
        assert profile == config
        assert alpha in config["data"]["alphas"] and resolution == 100
        inputs = list(map(Path, paths))
        metadata = [p.parent / "generated.json" for p in inputs]
        return inputs + metadata, labels

    monkeypatch.setattr(delivery, "generated_inputs", tiny_generated_inputs)
    calls = {name: [] for name in ("preprocess", "fit", "eht", "background", "aperture", "limits")}

    def spy(module, name, key):
        original = getattr(module, name)

        def wrapped(*args, **kwargs):
            calls[key].append(deepcopy(kwargs))
            return original(*args, **kwargs)

        monkeypatch.setattr(module, name, wrapped)

    spy(raptor, "preprocess_cube", "preprocess")
    spy(raptor, "fit_processed", "fit")
    spy(raptor, "estimate_background", "background")
    spy(eht, "run_eht_stage", "eht")

    def cheap_aperture(*, shape, **kwargs):
        calls["aperture"].append(dict(shape=shape, **kwargs))
        return tiny_operator(shape, bins=kwargs["angular_bins"])

    def capture_limits(observations, amplitude, phase, *, alphas, dimension, mode, options):
        assert set(amplitude) == set(phase) == set(alphas)
        assert all(amplitude[a].shape == phase[a].shape for a in alphas)
        assert len(observations.selected_indices) > 0
        calls["limits"].append((dimension, mode, options))
        # Engineering sentinel, NOT a computed scientific limit.
        return np.asarray(alphas), np.full(len(alphas), -0.25), {"synthetic_wiring_only": True}

    monkeypatch.setattr(polarization, "aperture_operator", cheap_aperture)
    monkeypatch.setattr(inference, "calculate_from_products", capture_limits)
    # Do not use stage='all': that also produces figures, which are irrelevant
    # to parameter wiring. Every selected stage below executes its real body.
    for stage in ("raptor", "eht", "background", "observations", "inference"):
        assert workflow.run_stage(path, out, input_root=raw, stage=stage, trust_pickle=True) == out

    theta, tg_hours = independent_scales(config)
    fov = 128.0 / theta
    radius = fraction * fov / 2
    sigma = fwhm / theta / np.sqrt(8 * np.log(2))
    expected_preprocess = dict(
        sigma=sigma, init_fov=config["image"]["initial_fov_rg"], final_fov=fov,
        final_shape=8, rotation_angle=np.deg2rad(37.0),
    )
    assert len(calls["preprocess"]) == 3
    for actual in calls["preprocess"][:2]:
        assert actual == pytest.approx(expected_preprocess)
    background_preprocess = dict(calls["preprocess"][2])
    assert background_preprocess.pop("state_indices") == [0]
    assert background_preprocess == pytest.approx(expected_preprocess)
    assert len(calls["fit"]) == 2
    for actual in calls["fit"]:
        assert actual == dict(
            fov=pytest.approx(fov), radial_max=pytest.approx(radius), snapshot_labels=labels,
            weighting=config["aperture"]["weighting"], window_degrees=10, angular_bins=360,
        )
    assert len(calls["aperture"]) == 3  # two RAPTOR fits and the actual EHT stage
    for actual in calls["aperture"]:
        assert actual == dict(shape=8, fov=pytest.approx(fov), radial_max=pytest.approx(radius),
                              window_degrees=10, angular_bins=360)
    assert len(calls["eht"]) == 1
    actual_eht = calls["eht"][0]
    assert actual_eht["fov"] == pytest.approx(fov)
    assert actual_eht["radial_max"] == pytest.approx(radius)
    assert actual_eht["expected_shape"] == 8
    assert actual_eht["angular_fov_uas"] == config["image"]["final_fov_uas"] == 128.0
    assert actual_eht["pixel_size_uas"] == 16.0
    assert actual_eht["angle_roll"] == 3
    for key in ("rank_offsets", "final_offsets"):
        assert actual_eht[key] == {int(day): value for day, value in contract[key].items()}
    assert actual_eht["exclude_columns"] == contract["exclude_columns"]
    assert actual_eht["exclusion_rectangles"] == contract["exclusion_rectangles"]
    eht_report = json.loads((out / "eht/report.json").read_text())
    assert eht_report["geometry"] == dict(
        mode="declared", expected_shape=8, angular_fov_uas=128.0, pixel_size_uas=16.0,
        fov_rg=pytest.approx(fov),
    )

    assert len(calls["background"]) == 1
    assert calls["background"][0]["cadence_hours"] == pytest.approx(footer_step * tg_hours)
    assert calls["background"][0]["lag"] == config["background"]["lag_frames"]
    summary = json.loads((out / "background/summary.json").read_text())
    assert summary["cadence_hours"] == pytest.approx(footer_step * tg_hours)
    assert summary["cadence_hours"] != pytest.approx(0.112411)
    assert {(dim, mode) for dim, mode, _ in calls["limits"]} == {
        (dim, mode) for dim in ("1d", "2d") for mode in ("std", "cor", "cor_grmhd")
    }
    assert len(calls["limits"]) == 6
    for _, _, options in calls["limits"]:
        assert options.time_conversion == pytest.approx(1 / tg_hours)
        assert options.time_conversion != pytest.approx(0.112411)
        assert options.reference_norm == 0.27 != 0.1
        assert options.cloud_coefficient == config["cloud"]["coefficient"]
        assert options.cloud_pivot_alpha == config["cloud"]["pivot_alpha"]
        assert options.cloud_exponent == config["cloud"]["exponent"]
        assert options.exposure_hours == 2.5
        assert options.day_gap_hours == 23.0
        assert options.pair_separation_hours == 97.0
        assert options.n_principal == config["inference"]["principal_n"]
        assert options.grmhd_sigma_degrees == pytest.approx(np.rad2deg(summary["sigma"]))
    manifest = verify_stage(out / "inference", recursive=True)
    assert manifest["config"]["options"]["reference_norm"] == 0.27
    assert manifest["config"]["options"]["time_conversion"] == pytest.approx(1 / tg_hours)
    assert manifest["config"]["cloud"] == config["cloud"]
    assert manifest["config"]["derived_geometry"]["image_preprocess"] == pytest.approx(expected_preprocess)
    assert not (out / "figures").exists()
    assert config == before

    # The real inference entrypoint must reject changed upstream geometry,
    # rather than silently calculate with products made under a different M/D.
    stale = deepcopy(config)
    stale["physical"]["mass_solar"] *= 1.1
    with pytest.raises(ValueError, match="differs from completed EHT"):
        inference.run_inference(out, stale, rebuild=True)
    assert len(calls["limits"]) == 6

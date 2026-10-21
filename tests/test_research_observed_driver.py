"""Tiny synthetic driver checks: no raw images, producer jobs or scientific inference."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from eht_axion.pipeline.common import SCHEMA_VERSION, verify_stage
from eht_axion.pipeline.constants import analysis_geometry
from eht_axion.pipeline.layout import alpha_label
import eht_axion.producer.delivery as delivery

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/research_observed_limits.py"
SPEC = importlib.util.spec_from_file_location("research_observed_driver", SCRIPT)
driver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(driver)


@pytest.fixture
def synthetic(tmp_path, monkeypatch):
    # Read the actual full producer profile, including all eight alpha identities.
    # No config patching, raw producer files, geometry kernels or large maps.
    args = argparse.Namespace(config=ROOT / "configs/eht-observed-research.yaml",
                              study=tmp_path / "study", templates=tmp_path / "theory",
                              output_root=tmp_path / "output", variant="all4800",
                              ridge=.25, resultant_min=.7, pixel_stride=5)
    config = driver.read_config(args.config)
    geometry = analysis_geometry(config)
    records, verified, validated, staged = {}, [], [], []

    def manifest(directory, stage, cfg, files):
        directory.mkdir(parents=True, exist_ok=True)
        outputs = []
        for name, content in files.items():
            path = directory / name
            path.write_bytes(content)
            outputs.append(dict(path=name, size=len(content),
                                sha256=hashlib.sha256(content).hexdigest()))
        record = dict(schema=SCHEMA_VERSION, complete=True, stage=stage,
                      config=deepcopy(cfg), inputs=[], outputs=outputs)
        records[directory] = record
        (directory / "manifest.json").write_text(json.dumps(record))

    consumed = [args.study / "report.json"]
    for day in driver.DAYS:
        consumed.extend(args.study / name for name in (
            f"{day}_all4800_indices.npy", f"{day}_all4800_summary.npz",
            f"{day}_candidates.npz", f"{day}_aligned_stokes.npy"))
    manifest(args.study, "research:eht-polarization", {"additional_blur_fwhm_uas": 19},
             {path.name: b'{"partial": false}' if path.name == "report.json" else b"synthetic"
              for path in consumed})
    for alpha in config["data"]["alphas"]:
        label = alpha_label(alpha)
        manifest(args.templates / "raptor/processed" / label, "raptor:preprocess",
                 geometry["image_preprocess"], {"placeholder": b"synthetic"})
        manifest(args.templates / "raptor/response" / label, "raptor:aperturefit",
                 dict(geometry["aperturefit"], snapshot_labels=[1600, 1610]),
                 {"placeholder": b"synthetic"})
    for dim, shape in (("1d", (360,)), ("2d", (5, 5))):
        directory = args.templates / "raptor/templates" / dim
        directory.mkdir(parents=True)
        values = np.ones((len(config["data"]["alphas"]),) + shape)
        np.savez(directory / "templates.npz", alphas=config["data"]["alphas"],
                 amplitude_mean=values, amplitude_std=values * .01,
                 phase_mean=values * .1, phase_std=values * .01)
        manifest(directory, "raptor:summarize", {},
                 {"templates.npz": (directory / "templates.npz").read_bytes()})

    def verify(directory, **kwargs):
        verified.append((directory, kwargs))
        return records[directory]

    def validate(record, cfg, *, alpha, resolution, verified):
        assert record == records[args.templates / "raptor/processed" / alpha_label(alpha)]
        assert cfg == config and resolution == 100
        assert isinstance(verified, set)
        validated.append(alpha)

    @contextmanager
    def stage(directory, name, identity, inputs):
        assert directory == args.output_root
        staged.append((name, identity, list(inputs)))
        directory.mkdir()
        yield directory

    def unexpected(*args, **kwargs):
        pytest.fail("Rejected input reached sample preparation or inference")

    monkeypatch.setattr(driver, "verify_stage", verify)
    monkeypatch.setattr(delivery, "validate_processed_source", validate)
    monkeypatch.setattr(driver, "stage_directory", stage)
    monkeypatch.setattr(driver, "prepare_samples", unexpected)
    monkeypatch.setattr(driver, "calculate_observed", unexpected)
    return SimpleNamespace(args=args, config=config, geometry=geometry, records=records,
                           consumed=consumed, verified=verified, validated=validated, staged=staged)


def test_moments_local_pi_branch_and_population_covariance():
    mean = np.array([np.pi / 2 - .02, -.3])
    residual = np.array([[-.1, -.2], [0, .1], [.1, .1]])
    values = mean + residual
    values += np.pi * np.array([[0, 1], [-1, 0], [2, -2]])
    actual, covariance, resultant = driver.moments(values)
    np.testing.assert_allclose(driver.circular_difference(actual, mean), 0, atol=1e-14)
    np.testing.assert_allclose(covariance, [[.02 / 3, .01], [.01, .02]], atol=1e-14)
    expected_resultant = [(1 + 2 * np.cos(.2)) / 3,
                          abs(np.exp(-.4j) + 2 * np.exp(.2j)) / 3]
    np.testing.assert_allclose(resultant, expected_resultant, atol=1e-14)


@pytest.mark.parametrize("problem", ["extra-undeclared-study", "missing-study", "undeclared-template"])
def test_consumed_files_must_be_present_and_declared(synthetic, monkeypatch, problem):
    s = synthetic
    if problem == "undeclared-template":
        directory = s.args.templates / "raptor/templates/2d"
        target = directory / "templates.npz"
        s.records[directory]["outputs"] = [dict(path="placeholder", size=1,
                                                   sha256=hashlib.sha256(b"x").hexdigest())]
        (directory / "placeholder").write_bytes(b"x")
    else:
        directory = s.args.study
        target = directory / "95_candidates.npz"
        if problem == "missing-study":
            target.unlink()
        else:
            # The file exists, but verification of other outputs must not bless it.
            s.records[directory]["outputs"] = [e for e in s.records[directory]["outputs"]
                                                if e["path"] != target.name]
    (directory / "manifest.json").write_text(json.dumps(s.records[directory]))
    stub = driver.verify_stage

    def verify_with_real_integrity(path, **kwargs):
        stub(path, **kwargs)
        return verify_stage(path, **kwargs)

    monkeypatch.setattr(driver, "verify_stage", verify_with_real_integrity)
    message = "Unsafe/missing stage output" if problem == "missing-study" else "not a declared stage output"
    with pytest.raises(ValueError, match=message):
        driver.run(s.args)
    assert not s.args.output_root.exists()
    assert not s.staged
    if problem != "missing-study":
        assert target.is_file()


@pytest.mark.parametrize("problem, message", [
    ("study-stage", "Unexpected study stage"),
    ("template-stage", "Unexpected theory template stage"),
    ("partial", "partial candidate study"),
    ("kernel", "sky/aperture/kernel"),
    ("preprocess-profile", "Theory preprocessing differs"),
    ("aperture-profile", "Theory aperture differs"),
])
def test_rejects_wrong_identity_or_geometry_before_publication(synthetic, problem, message):
    s = synthetic
    label = alpha_label(s.config["data"]["alphas"][0])
    if problem == "study-stage":
        s.records[s.args.study]["stage"] = "raptor:summarize"
    elif problem == "template-stage":
        s.records[s.args.templates / "raptor/templates/2d"]["stage"] = "raptor:aperturefit"
    elif problem == "partial":
        (s.args.study / "report.json").write_text('{"partial": true}')
    elif problem == "kernel":
        s.records[s.args.study]["config"]["additional_blur_fwhm_uas"] = 20
    elif problem == "preprocess-profile":
        s.records[s.args.templates / "raptor/processed" / label]["config"]["rotation_angle"] += .1
    else:
        s.records[s.args.templates / "raptor/response" / label]["config"]["window_degrees"] = 20
    with pytest.raises(ValueError, match=message):
        driver.run(s.args)
    assert not s.args.output_root.exists()
    assert not s.staged


def test_run_passes_actual_differences_and_same_covariance_to_null(synthetic, monkeypatch):
    s = synthetic
    residual = np.array([[-.1, -.2], [0, .1], [.1, .1]])
    centers = np.array([[1.54, -.2], [1.64, -.1], [1.4, .1], [1.2, .4]])
    scales = np.array([1, 2, 1.5, .5])
    samples = [np.column_stack((center + scale * residual, [0, np.pi / 3, -np.pi / 3]))
               for center, scale in zip(centers, scales)]
    prepared, calls = [], []

    def prepare(study, args, dim, valid):
        assert study == s.args.study and args is s.args
        assert valid.shape == ((360,) if dim == "1d" else (5, 5))
        assert valid.all()
        prepared.append(dim)
        return np.array([4, 14, 24]), samples

    def calculate(obs, amp, phase, *, alphas, dimension, mode, options):
        assert list(amp) == list(phase) == s.config["data"]["alphas"]
        np.testing.assert_array_equal(alphas, s.config["data"]["alphas"])
        assert options.phase_count == s.config["inference"]["phase_count"]
        assert options.reference_norm == s.config["cloud"]["reference_norm"]
        calls.append((dimension, mode, obs))
        return alphas, np.zeros(len(alphas)), {"synthetic": True}

    monkeypatch.setattr(driver, "prepare_samples", prepare)
    monkeypatch.setattr(driver, "calculate_observed", calculate)
    driver.run(s.args)
    assert prepared == ["1d", "2d"]
    assert len(calls) == 8
    raw_covariance = np.array([[.02 / 3, .01], [.01, .02]])
    shrunk = .75 * raw_covariance + .25 * np.diag(np.diag(raw_covariance))
    expected_covariances = scales[:, None, None] ** 2 * shrunk
    expected_differences = np.array([[.1, .1], [-.2, .3]])
    for i in range(0, len(calls), 2):
        dim, mode, obs = calls[i]
        null_dim, null_mode, null = calls[i + 1]
        assert (dim, mode) == (null_dim, null_mode) == (("1d", "2d")[i // 4],
                                                       ("std", "cor")[(i // 2) % 2])
        np.testing.assert_array_equal(obs.selected_indices, [4, 14])
        np.testing.assert_allclose(driver.circular_difference(np.array(obs.means), centers), 0, atol=1e-14)
        np.testing.assert_allclose(obs.covariances, expected_covariances, atol=1e-14)
        np.testing.assert_allclose(obs.standard_deviations,
                                   np.sqrt(np.diagonal(expected_covariances, axis1=1, axis2=2)))
        np.testing.assert_allclose(obs.pair_differences, expected_differences, atol=1e-14)
        np.testing.assert_allclose(obs.pair_covariances,
                                   [expected_covariances[0] + expected_covariances[1],
                                    expected_covariances[2] + expected_covariances[3]])
        np.testing.assert_array_equal(null.pair_differences, np.zeros((2, 2)))
        assert null.covariances is obs.covariances
        assert null.pair_covariances is obs.pair_covariances
        assert null.means is obs.means

    expected_verified = [s.args.study]
    for alpha in s.config["data"]["alphas"]:
        expected_verified.extend(s.args.templates / f"raptor/{kind}" / alpha_label(alpha)
                                 for kind in ("processed", "response"))
    expected_verified.extend(s.args.templates / "raptor/templates" / dim for dim in ("1d", "2d"))
    assert [path for path, _ in s.verified] == expected_verified
    assert all(kw["recursive"] and kw["allow_imported"] is False for _, kw in s.verified)
    shared_seen = s.verified[1][1]["_seen"]
    assert all(kw["_seen"] is shared_seen for _, kw in s.verified[1:])
    assert s.validated == s.config["data"]["alphas"]
    assert len(s.staged) == 1
    name, identity, inputs = s.staged[0]
    assert name == "research:observed-limits"
    assert identity["scientific_profile"] == s.config
    assert identity["variant"] == s.args.variant and identity["ridge"] == s.args.ridge
    expected_inputs = {s.args.config.resolve(), SCRIPT.resolve(), s.args.study / "manifest.json",
                       *s.consumed}
    for dim in ("1d", "2d"):
        directory = s.args.templates / "raptor/templates" / dim
        expected_inputs.update([directory / "manifest.json", directory / "templates.npz"])
        with np.load(s.args.output_root / f"observations-{dim}.npz", allow_pickle=False) as saved:
            np.testing.assert_array_equal(saved["selected_indices"], [4, 14])
            np.testing.assert_allclose(saved["pair_differences"], expected_differences, atol=1e-14)
            np.testing.assert_array_equal(saved["sample_counts"], [3, 3, 3, 3])
    assert set(inputs) == expected_inputs
    report = json.loads((s.args.output_root / "report.json").read_text())
    assert set(report["results"]) == {f"{label}-{mode}-{dim}"
                                      for label in ("observed", "null-comparison")
                                      for mode in ("std", "cor") for dim in ("1d", "2d")}

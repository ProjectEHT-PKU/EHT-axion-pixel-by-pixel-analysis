"""Whole tiny-chain lineage regressions; no historical/real science inputs."""

from contextlib import contextmanager
from pathlib import Path
import json
import numpy as np
import pytest

from test_pipeline_end_to_end import make_inputs
from eht_axion.pipeline.workflow import run_stage, read_config
from eht_axion.pipeline.inference import run_inference
from eht_axion.pipeline.convergence import check_convergence
from eht_axion.pipeline.common import verify_stage, describe_inputs
from eht_axion.pipeline.figures import plot_results


@contextmanager
def restored(path):
    before = path.read_bytes()
    try:
        yield
    finally:
        path.write_bytes(before)


@pytest.fixture(scope="module")
def tiny_chain(tmp_path_factory):
    root = tmp_path_factory.mktemp("integrity")
    raw, out = root / "input", root / "output"
    make_inputs(raw)
    source = Path(__file__).resolve().parents[1] / "configs/research.yaml"
    text = source.read_text()
    replacements = {
        "alphas: [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]": "alphas: [0.1, 0.3]",
        "final_shape: 64": "final_shape: 8", "topset_size: 1572": "topset_size: 40",
        "min_radius_rg: 2.0": "min_radius_rg: 0.0", "historical_mask: true": "historical_mask: false",
        "covariance_ridge: 0.0": "covariance_ridge: 0.001",
        "source_resolution: 400": "source_resolution: 100", "phase_count: 64": "phase_count: 16",
        "coupling_points: 120": "coupling_points: 30",
        # Keep the pre-N2 numerical oracle explicit; production defaults are N=2.
        "principal_n: 2": "principal_n: 4",
    }
    for before, after in replacements.items():
        assert before in text
        text = text.replace(before, after)
    config = root / "config.yaml"
    config.write_text(text)
    run_stage(config, out, input_root=raw, trust_pickle=True)
    return root, raw, out, config


def test_valid_chain_preserves_verified_numerics(tiny_chain):
    _, _, out, _ = tiny_chain
    # Historical N=4 engineering oracle, explicitly selected by tiny_chain.
    # The 211 default correction has separate default/analytic regression coverage.
    expected = {
        "1d": [[-0.0001886100856455447, -1.1073034273258424],
               [0.13089591098059092, -0.9823464577669513],
               [0.26385984920594635, -0.85723284116235]],
        "2d": [[0.3086294792552724, -0.8159039273929453],
               [0.1917774745621371, -0.9296676922787799],
               [0.2547205604436642, -0.8663517220227757]],
    }
    for dim, cuts in expected.items():
        with np.load(out / f"inference/limits-{dim}.npz") as z:
            np.testing.assert_allclose(z["cuts"], cuts, rtol=0, atol=1e-12)
    verify_stage(out / "inference", recursive=True, allow_imported=False)
    for dim in ("1d", "2d"):
        m = verify_stage(out / "raptor/templates" / dim)
        assert sum(Path(x["path"]).name == "manifest.json" for x in m["inputs"]) == 2


def test_stale_raw_refused_by_actual_inference(tiny_chain):
    _, raw, out, config = tiny_chain
    path = raw / "raptor-stokes/res-100/alpha-0.1.npy"
    stamp = (out / "inference/manifest.json").read_bytes()
    with restored(path):
        x = np.load(path)
        x[..., 1] *= .91
        np.save(path, x)
        with pytest.raises(ValueError, match="Stale upstream lineage"):
            run_inference(out, read_config(config), rebuild=True)
    assert (out / "inference/manifest.json").read_bytes() == stamp


def test_imported_response_cannot_evade_template_chain(tiny_chain):
    _, _, out, config = tiny_chain
    path = out / "raptor/response/alpha-0.10/manifest.json"
    with restored(path):
        value = json.loads(path.read_text())
        value["imported_from"] = {"path": "archive"}
        path.write_text(json.dumps(value))
        with pytest.raises(ValueError, match="Stale upstream lineage|reuse-imported"):
            run_inference(out, read_config(config), rebuild=True)


def test_convergence_refuses_stale_template(tiny_chain):
    root, _, out, _ = tiny_chain
    path = out / "raptor/templates/1d/templates.npz"
    target = root / "bad-convergence.json"
    with restored(path):
        with np.load(path) as z:
            value = {k: z[k] for k in z.files}
        value["amplitude_mean"] *= 2
        np.savez(path, **value)
        with pytest.raises(ValueError, match="Stale upstream lineage|checksum"):
            check_convergence(out, target, dimensions=("1d",), modes=("std",))
    assert not target.exists()


def test_convergence_rejects_undeclared_run_config(tiny_chain):
    root, _, out, _ = tiny_chain
    path = out / "inference/manifest.json"
    with restored(path):
        value = json.loads(path.read_text())
        value["outputs"] = [x for x in value["outputs"] if x["path"] != "run.json"]
        path.write_text(json.dumps(value))
        with pytest.raises(ValueError, match="not a declared stage output"):
            check_convergence(out, root / "undeclared-run-convergence.json",
                              dimensions=("1d",), modes=("std",))
    assert not (root / "undeclared-run-convergence.json").exists()


def test_convergence_matches_baseline_first(tiny_chain):
    root, _, out, _ = tiny_chain
    report = check_convergence(out, root / "convergence.json", dimensions=("1d",), modes=("std",))
    assert report["checks"]["std_1d"]["baseline_max_delta_dex"] < 1e-12


@pytest.mark.parametrize("tool", ["inference", "figures", "convergence"])
def test_old_missing_template_edges_fail_closed(tiny_chain, tool):
    root, raw, out, config = tiny_chain
    manifest_paths = [p for p in out.rglob("manifest.json") if p.parent != out]
    saved = {p: p.read_bytes() for p in manifest_paths}
    raw_path = raw / "raptor-stokes/res-100/alpha-0.1.npy"
    raw_before = raw_path.read_bytes()
    try:
        # Simulate historical manifests with bare response files only, then
        # rehash declared downstream edges: ordinary recursion cannot see raw.
        for dim in ("1d", "2d"):
            path = out / "raptor/templates" / dim / "manifest.json"
            value = json.loads(path.read_text())
            value["inputs"] = [x for x in value["inputs"] if Path(x["path"]).suffix == ".npz"]
            path.write_text(json.dumps(value))
        for stage in (out / "observations/1d", out / "observations/2d", out / "inference"):
            path = stage / "manifest.json"
            value = json.loads(path.read_text())
            value["inputs"] = describe_inputs([x["path"] for x in value["inputs"]])
            path.write_text(json.dumps(value))
        x = np.load(raw_path)
        x[..., 1] *= .9
        np.save(raw_path, x)
        with pytest.raises(ValueError, match="Stale upstream lineage|different run/root"):
            if tool == "inference":
                run_inference(out, read_config(config), rebuild=True)
            elif tool == "figures":
                plot_results(out / "inference", root / "old-edge-figures")
            else:
                check_convergence(out, root / "old-edge-convergence.json",
                                  dimensions=("1d",), modes=("std",))
    finally:
        raw_path.write_bytes(raw_before)
        for path, before in saved.items():
            path.write_bytes(before)


def test_figures_cannot_attach_external_inference_to_local_run(tiny_chain):
    root, _, out, _ = tiny_chain
    external = root / "unrelated.txt"
    external.write_text("not this run")
    path = out / "inference/manifest.json"
    with restored(path):
        value = json.loads(path.read_text())
        value["inputs"] = describe_inputs([external])
        path.write_text(json.dumps(value))
        with pytest.raises(ValueError, match="different run/root"):
            plot_results(out / "inference", root / "foreign-figures")
    assert not (root / "foreign-figures").exists()


def test_output_root_link_never_rebuilt(tiny_chain):
    root, _, out, config = tiny_chain
    alias = root / "run-alias"
    alias.symlink_to(out, target_is_directory=True)
    original = (out / "figures/manifest.json").read_bytes()
    with pytest.raises(ValueError, match="shared directory symlink"):
        run_stage(config, alias, stage="figures", rebuild=True)
    assert (out / "figures/manifest.json").read_bytes() == original

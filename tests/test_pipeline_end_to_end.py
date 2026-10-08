"""Tiny oracle-like chain: no historical files, notebooks, hardcoded data locations."""

import pickle
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from eht_axion.pipeline.workflow import run_stage
from eht_axion.pipeline.common import verify_stage
from eht_axion.pipeline.figures import plot_results
from eht_axion.pipeline.products import load_observations


def run(args):
    return run_stage(
        args.config,
        args.output_root,
        input_root=args.input_root,
        stage=args.stage,
        trust_pickle=args.trust_pickle,
        rebuild=args.rebuild,
        reuse_imported=getattr(args, "reuse_imported", False),
    )


def make_inputs(root):
    rng = np.random.default_rng(128)
    shape, times, candidates = 8, 9, 60
    fields = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
    metadata = [{field: (i if field == "zbl" else 1) for field in fields} for i in range(candidates)]
    for day in (95, 96, 100, 101):
        path = root / "eht-image-ensemble" / f"230_{day}.pkl"
        path.parent.mkdir(parents=True, exist_ok=True)
        angle = 0.15 + rng.normal(0, 0.035, (shape, shape, candidates))
        intensity = np.ones_like(angle)
        order = rng.permutation(candidates)
        data = dict(
            I_data=intensity[..., order],
            Q_data=np.cos(2 * angle)[..., order],
            U_data=np.sin(2 * angle)[..., order],
            meta_data=[metadata[k] for k in order],
        )
        with path.open("wb") as stream:
            pickle.dump(data, stream)
    phase = np.arange(16) * 2 * np.pi / 16
    for alpha in (0.1, 0.3):
        image = np.zeros((times, 17, shape, shape, 4))
        for t in range(times):
            base = 0.2 + rng.normal(0, 0.02, (shape, shape))
            angle = base + 0.015 * np.cos(phase[:, None, None] + 0.3)
            image[t, 0, ..., 1] = np.cos(2 * base)
            image[t, 0, ..., 2] = np.sin(2 * base)
            image[t, 1:, ..., 1] = np.cos(2 * angle)
            image[t, 1:, ..., 2] = np.sin(2 * angle)
            image[t, ..., 0] = 1
        path = root / "raptor-stokes" / "res-100" / f"alpha-{alpha}.npy"
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path, image)


def test_script_only_chain_and_cache(tmp_path):
    raw, out = tmp_path / "input", tmp_path / "output"
    make_inputs(raw)
    template = Path(__file__).resolve().parents[1] / "configs" / "research.yaml"
    text = template.read_text()
    for before, after in [
        ("alphas: [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]", "alphas: [0.1, 0.3]"),
        ("final_shape: 64", "final_shape: 8"),
        ("topset_size: 1572", "topset_size: 40"),
        ("min_radius_rg: 2.0", "min_radius_rg: 0.0"),
        ("historical_mask: true", "historical_mask: false"),
        ("covariance_ridge: 0.0", "covariance_ridge: 0.001"),
        ("source_resolution: 400", "source_resolution: 100"),
        ("phase_count: 64", "phase_count: 16"),
        ("coupling_points: 120", "coupling_points: 30"),
    ]:
        assert before in text
        text = text.replace(before, after)
    config = tmp_path / "config.yaml"
    config.write_text(text)
    args = Namespace(
        config=config,
        input_root=raw,
        output_root=out,
        stage="all",
        trust_pickle=True,
        rebuild=False,
    )
    # Exercise the actual direct script from outside the repository, not only a helper.
    import os
    import subprocess
    import sys

    script = Path(__file__).resolve().parents[1] / "scripts/run_all.py"
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    for name in list(environment):
        if name.startswith("EHT_"):
            environment.pop(name)
    environment.update(
        EHT_CONFIG=str(config),
        EHT_INPUT_ROOT=str(raw),
        EHT_OUTPUT_ROOT=str(tmp_path / "must-not-be-used"),
        EHT_TRUST_PICKLE="true",
    )
    result = subprocess.run(
        [sys.executable, str(script), "--output-root", str(out)],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert not (tmp_path / "must-not-be-used").exists()
    curves = list((out / "inference").glob("limits-*.npz"))
    assert len(curves) == 2
    with np.load(out / "inference/limits-1d.npz", allow_pickle=False) as z:
        assert z["cuts"].shape == (3, 2)
        assert z["modes"].tolist() == ["std", "cor", "cor_grmhd"]
    assert (out / "figures/limits.png").stat().st_size > 10_000
    assert (out / "config.yaml").read_bytes() == config.read_bytes()
    assert not (out / "config.toml").exists()
    assert verify_stage(out / "figures")["complete"]
    assert (out / "eht/2017-04-05/stokes.npy").is_file()
    assert (out / "raptor/processed/alpha-0.10/stokes.npy").is_file()
    obs = load_observations(out / "observations/2d/observations.npz")
    for sigma, covariance in zip(obs.standard_deviations, obs.covariances):
        np.testing.assert_allclose(sigma**2, np.diag(covariance))
    stamp = (out / "inference/limits-1d.npz").stat().st_mtime_ns
    run(args)
    assert (out / "inference/limits-1d.npz").stat().st_mtime_ns == stamp
    root_identity = (out / "manifest.json").read_bytes()
    config.write_text(text.replace("prior: log-uniform", "prior: uniform"))
    args.stage = "figures"
    run(args)
    assert (out / "manifest.json").read_bytes() == root_identity
    assert len(list((out / "invocations").glob("*.json"))) == 3
    config.write_text(text)
    args.stage = "all"
    config.write_text(text.replace("weighting: intensity", "weighting: flux"))
    args.stage = "inference"
    args.rebuild = True
    with pytest.raises(ValueError, match="differs from completed EHT"):
        run(args)
    config.write_text(text.replace("source_resolution: 100", "source_resolution: 400"))
    with pytest.raises(ValueError, match="background source differs"):
        run(args)
    config.write_text(text.replace("alphas: [0.1, 0.3]", "alphas: [0.1]"))
    args.rebuild = False
    with pytest.raises(ValueError, match="alpha grid differs|lineage belongs"):
        run(args)
    args.rebuild = True
    config.write_text(text)
    with np.load(out / "inference/limits-2d.npz") as z:
        values = {k: z[k] for k in z.files}
    values["prior"] = np.asarray("uniform")
    np.savez(out / "inference/limits-2d.npz", **values)
    with pytest.raises(ValueError, match="Mixed"):
        plot_results(out / "inference", tmp_path / "bad-figures")
    with pytest.raises(ValueError, match="checksum"):
        verify_stage(out / "inference")

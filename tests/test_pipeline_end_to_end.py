"""Tiny oracle-like chain: no historical files, notebooks, hardcoded data locations."""

import pickle
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from eht_axion.pipeline.cli import run
from eht_axion.pipeline.common import verify_stage
from eht_axion.pipeline.figures import plot_results
from eht_axion.pipeline.products import load_observations


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
    template = Path(__file__).resolve().parents[1] / "configs" / "research.toml"
    text = template.read_text()
    for before, after in [
        ("alphas = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]", "alphas = [0.1, 0.3]"),
        ("final_shape = 64", "final_shape = 8"),
        ("topset_size = 1572", "topset_size = 40"),
        ("min_radius_rg = 2.0", "min_radius_rg = 0.0"),
        ("historical_mask = true", "historical_mask = false"),
        ("covariance_ridge = 0.0", "covariance_ridge = 0.001"),
        ("source_resolution = 400", "source_resolution = 100"),
        ("phase_count = 64", "phase_count = 16"),
        ("coupling_points = 120", "coupling_points = 30"),
    ]:
        assert before in text
        text = text.replace(before, after)
    config = tmp_path / "config.toml"
    config.write_text(text)
    args = Namespace(
        config=config,
        input_root=raw,
        output_root=out,
        stage="all",
        trust_pickle=True,
        rebuild=False,
    )
    run(args)
    curves = list((out / "07_inference").glob("*.npz"))
    assert len(curves) == 6
    assert (out / "08_figures" / "limits.png").stat().st_size > 10_000
    assert verify_stage(out / "08_figures")["complete"]
    obs = load_observations(out / "06_observations" / "2d" / "observations.npz")
    for sigma, covariance in zip(obs.standard_deviations, obs.covariances):
        np.testing.assert_allclose(sigma**2, np.diag(covariance))
    stamp = (out / "07_inference" / "std_1d.npz").stat().st_mtime_ns
    run(args)
    assert (out / "07_inference" / "std_1d.npz").stat().st_mtime_ns == stamp
    config.write_text(text.replace('weighting = "intensity"', 'weighting = "flux"'))
    args.stage = "inference"
    args.rebuild = True
    with pytest.raises(ValueError, match="differs from completed EHT"):
        run(args)
    config.write_text(text.replace("source_resolution = 100", "source_resolution = 400"))
    with pytest.raises(ValueError, match="background source differs"):
        run(args)
    config.write_text(text.replace("alphas = [0.1, 0.3]", "alphas = [0.1]"))
    args.rebuild = False
    with pytest.raises(ValueError, match="alpha grid differs|lineage belongs"):
        run(args)
    args.rebuild = True
    config.write_text(text)
    with np.load(out / "07_inference" / "cor_2d.npz") as z:
        values = {k: z[k] for k in z.files}
    values["prior"] = np.asarray("uniform")
    np.savez(out / "07_inference" / "cor_2d.npz", **values)
    with pytest.raises(ValueError, match="Mixed"):
        plot_results(out / "07_inference", tmp_path / "bad-figures")
    with pytest.raises(ValueError, match="checksum"):
        verify_stage(out / "07_inference")

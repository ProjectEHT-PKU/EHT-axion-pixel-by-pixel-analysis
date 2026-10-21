"""Synthetic-only checks; never opens production FITS or historical pickles."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import pickle
import subprocess
import sys

import numpy as np
import pytest
import yaml

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "research_eht_polarization.py"
SPEC = importlib.util.spec_from_file_location("research_eht_driver", SCRIPT)
driver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(driver)


def test_help_without_optional_imports():
    completed = subprocess.run([sys.executable, "-S", str(SCRIPT), "--help"],
                               check=True, capture_output=True, text=True)
    assert "--config" in completed.stdout
    assert "--output-root" in completed.stdout
    assert "--max-candidates" in completed.stdout


def test_canonical_filename_and_metadata(tmp_path):
    key = (.4, 40, .01, 10, 10, 100, 1000)
    assert driver.candidate_path(tmp_path, key).name == "0.4_40.0_0.01_10.0_10.0_100.0_1000.0.fits"
    row = dict(zip(driver.FIELDS, key), **dict.fromkeys(driver.METRICS, 1.0), chi2=1.0)
    path = tmp_path / "meta.pkl"
    path.write_bytes(pickle.dumps({"meta_data": [row], "I_data": np.zeros((2, 2, 1))}))
    assert list(driver.read_metadata(path)) == [key]
    row["chi2"] = 4
    path.write_bytes(pickle.dumps({"meta_data": [row]}))
    with pytest.raises(ValueError, match="chi2"):
        driver.read_metadata(path)


def test_threshold_intersection_empty_not_relaxed():
    keys = [(float(i),)*7 for i in range(3)]
    metadata = {d: {k: (np.full(5, 1.0 if i == j else 3.0), float(i+1))
                    for i, k in enumerate(keys)} for j, d in enumerate(driver.DAYS)}
    selected = driver.selections(keys, metadata, {keys[1]}, dict.fromkeys(driver.METRICS, 2))
    for j, day in enumerate(driver.DAYS):
        assert selected[day]["datafit_shared"].size == 0
        assert selected[day]["legacy_evpa"].tolist() == [1]
        assert selected[day]["datafit_perday"].tolist() == ([j] if j < 3 else [])


def test_bounded_translation_zero_padding():
    ref = np.zeros((20, 20, 4))
    ref[7:10, 8:11, 0] = 1
    moved = driver.translated(ref, (2, -3))
    offset, correlation = driver.translation(ref[..., 0], moved[..., 0], 4)
    assert offset.tolist() == [-2, 3]
    assert correlation == pytest.approx(1)
    edge = np.zeros((8, 8, 4))
    edge[-1, -1] = 1
    assert driver.translated(edge, (1, 1)).sum() == 0


def test_sky_curves_reuses_existing_operator():
    from eht_axion.polarization import aperture_operator, integrate_stokes

    image = np.zeros((8, 8, 4))
    image[..., 0] = 1
    sky_angle = np.deg2rad(32)
    image[..., 1], image[..., 2] = np.cos(2*sky_angle), np.sin(2*sky_angle)
    operator = aperture_operator(shape=8, fov=128, radial_max=64,
                                 radial_resolution=4, angular_resolution_degrees=2)
    actual, brightness = driver.sky_curves(image, operator, 256)
    np.testing.assert_allclose(actual, sky_angle, atol=1e-12)
    assert np.all(brightness > 0)
    old = image.transpose(1, 0, 2).copy()
    old[..., 1:3] *= -1
    expected = driver.wrap(np.roll(integrate_stokes(old, operator, "intensity"), -90)-np.pi/2)
    np.testing.assert_allclose(actual, expected)


def test_bright_mask_and_actual_day_difference():
    image = np.zeros((2, 2, 4))
    image[..., 0] = [[1, .09], [.5, .5]]
    image[..., 1] = [[1, 1], [.19, .2]]
    assert driver.bright_mask(image).tolist() == [[True, False], [False, True]]
    assert np.rad2deg(driver.wrap(np.deg2rad(-178))) == pytest.approx(2)


@pytest.fixture
def synthetic_config(tmp_path):
    fits = pytest.importorskip("astropy.io.fits")
    keys = [(.4, 40, .01, 10, 10, 100, 1000), (.5, 40, .01, 10, 10, 100, 1000)]
    cfg = {"expected_candidates": 2, "legacy_selection": str(tmp_path / "selection.npz"),
           "thresholds": dict.fromkeys(driver.METRICS, 2), "days": {},
           "additional_blur_fwhm_uas": 19, "max_shift_pixels": 3}
    np.savez(cfg["legacy_selection"], selected_identity=np.asarray([json.dumps(keys[1])]))
    yy, xx = np.mgrid[:64, :64]
    intensity = np.exp(-((xx-31.5)**2+(yy-31.5)**2)/90)
    for j, day in enumerate(driver.DAYS):
        folder = tmp_path / str(day)
        folder.mkdir()
        metadata = []
        for i, key in enumerate(keys):
            row = dict(zip(driver.FIELDS, key), **dict.fromkeys(driver.METRICS, 1+i), chi2=1+i)
            metadata.append(row)
            angle = np.deg2rad(j*5)
            image = np.stack([intensity, .2*intensity*np.cos(2*angle),
                              .2*intensity*np.sin(2*angle), intensity*0], axis=-1)
            header = fits.Header({"CDELT1": -2/3.6e9, "CDELT2": 2/3.6e9,
                                  "CRPIX1": 32.5, "CRPIX2": 32.5,
                                  "CRVAL1": 187.7, "CRVAL2": 12.4, "FREQ": 230e9,
                                  "CTYPE1": "RA---SIN", "CTYPE2": "DEC--SIN", "BUNIT": "JY/PIXEL"})
            hdus = [fits.PrimaryHDU(image[..., 0], header)]
            hdus += [fits.ImageHDU(image[..., k], name=name) for k, name in enumerate("IQUV") if k]
            fits.HDUList(hdus).writeto(driver.candidate_path(folder, key))
        meta = tmp_path / f"{day}.pkl"
        meta.write_bytes(pickle.dumps({"meta_data": metadata, "I_data": np.zeros((2, 2, 2))}))
        cfg["days"][day] = {"fits_dir": str(folder), "metadata": str(meta)}
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))
    return config


def test_config_safe_and_explicit(synthetic_config, tmp_path):
    cfg = driver.load_config(synthetic_config)
    assert set(cfg["days"]) == set(driver.DAYS)
    cfg["days"][95]["metadata"] = "relative.pkl"
    synthetic_config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="absolute"):
        driver.load_config(synthetic_config)
    bad = tmp_path / "bad.yaml"
    bad.write_text("!!python/object/apply:os.system ['false']")
    with pytest.raises(yaml.YAMLError):
        driver.load_config(bad)


def test_unknown_fits_fails_instead_of_skipping(synthetic_config, tmp_path):
    pytest.importorskip("eht_axion.research.polarimetry")
    (tmp_path / "95" / "unrecognized.fits").write_bytes(b"not a FITS")
    args = argparse.Namespace(config=synthetic_config, output_root=tmp_path / "bad-study",
                              trust_pickle=True, max_candidates=1, study_id="bad")
    with pytest.raises(ValueError, match="filename mismatch.*unexpected"):
        driver.run(args)
    assert not args.output_root.exists()


def test_synthetic_publication(synthetic_config, tmp_path, monkeypatch):
    pytest.importorskip("eht_axion.research.polarimetry")
    # Exact production kernel is checked above; avoid its expensive construction in fixture runs.
    import eht_axion.polarization as polarization
    original = polarization.aperture_operator
    monkeypatch.setattr(polarization, "aperture_operator", lambda **kw: original(
        **kw, radial_resolution=2, angular_resolution_degrees=2))
    args = argparse.Namespace(config=synthetic_config, output_root=tmp_path / "study",
                              trust_pickle=True, max_candidates=1, study_id="synthetic")
    report = driver.run(args)
    assert report["partial"] and report["research_only"]
    assert not report["official_reproduction"]
    out = args.output_root
    assert (out / "manifest.json").is_file()
    maps = np.load(out / "95_aligned_stokes.npy", mmap_mode="r")
    assert maps.shape == (1, 64, 64, 4) and maps.dtype == np.float64
    assert np.load(out / "95_aperture.npy", allow_pickle=False).shape == (1, 360)
    assert np.load(out / "95_legacy_evpa_indices.npy").size == 0
    with np.load(out / "95_all4800_summary.npz", allow_pickle=False) as summary:
        assert summary["mean_stokes"].shape == (2, 64, 64, 4)
        assert summary["mean_stokes_curve_chi"].shape == (2, 360)
    with np.load(out / "all4800_aligned_95_96_daydiff.npz", allow_pickle=False) as data:
        np.testing.assert_allclose(np.rad2deg(data["curve_delta_chi"]), 5, atol=1e-6)
        assert data["mask"].dtype == bool
    with pytest.raises(FileExistsError):
        driver.run(args)
    args.output_root = tmp_path / "untrusted"
    args.trust_pickle = False
    with pytest.raises(ValueError, match="trust-pickle"):
        driver.run(args)
    assert not args.output_root.exists()

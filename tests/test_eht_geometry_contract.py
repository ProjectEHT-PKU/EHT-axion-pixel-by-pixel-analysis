"""Small synthetic input-contract tests; never read historical ensembles."""

import json
import pickle

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from eht_axion.pipeline import eht


@pytest.mark.parametrize("declaration", [
    {"expected_shape": 64},
    {"angular_fov_uas": 128.0},
    {"pixel_size_uas": 2.0},
    {"expected_shape": 64, "angular_fov_uas": 128.0},
    {"expected_shape": 64, "pixel_size_uas": 2.0},
    {"angular_fov_uas": 128.0, "pixel_size_uas": 2.0},
])
def test_geometry_declaration_is_all_or_none_before_io(tmp_path, monkeypatch, declaration):
    def forbidden(*args, **kwargs):
        pytest.fail("invalid geometry reached input IO")

    monkeypatch.setattr(eht, "input_eht_file", forbidden)
    monkeypatch.setattr(eht, "_load_day", forbidden)
    with pytest.raises(ValueError, match="all supplied or all None"):
        eht.run_eht_stage(tmp_path / "missing", tmp_path / "out", trust_pickle=True, **declaration)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("field,value", [
    *(('expected_shape', value) for value in (True, False, 1, 0, -2, 64.0, "64", np.nan, np.inf)),
    *((field, value) for field in ("angular_fov_uas", "pixel_size_uas")
      for value in (True, False, 0, -1, "2", np.nan, np.inf, -np.inf, 10 ** 400)),
    ("expected_shape", 10 ** 400),  # product overflows, not an accepted infinite extent
    ("angular_fov_uas", 120),
])
def test_bad_geometry_rejected_before_data_read(tmp_path, monkeypatch, field, value):
    declaration = dict(expected_shape=64, angular_fov_uas=128.0, pixel_size_uas=2.0)
    declaration[field] = value

    def forbidden(*args, **kwargs):
        pytest.fail("invalid geometry reached input IO")

    monkeypatch.setattr(eht, "_load_day", forbidden)
    monkeypatch.setattr(eht, "input_eht_file", forbidden)
    with pytest.raises(ValueError):
        eht.run_eht_stage(tmp_path / "missing", tmp_path / "out", trust_pickle=True, **declaration)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("field", ["topset_size", "batch_size", "pixel_batch_size", "angular_bins"])
@pytest.mark.parametrize("value", [0, -1, True, 2.5, np.nan, np.inf])
def test_bad_sizes_rejected_before_data_read(tmp_path, monkeypatch, field, value):
    def forbidden(*args, **kwargs):
        pytest.fail("invalid sizes reached input IO")

    monkeypatch.setattr(eht, "_load_day", forbidden)
    with pytest.raises(ValueError, match="sizes must be positive integers"):
        eht.run_eht_stage(tmp_path / "missing", tmp_path / "out", trust_pickle=True, **{field: value})
    assert not (tmp_path / "out").exists()


def test_validate_contract_needs_no_mass_distance_and_never_rounds_geometry():
    assert eht.validate_geometry_contract() is None
    assert eht.validate_geometry_contract(64, 128.0, 2.0) is None
    assert eht.validate_geometry_contract(np.int64(64), np.float64(128), np.float64(2)) is None
    # Relative tolerance admits FITS decimal roundoff, not a different FOV.
    assert eht.validate_geometry_contract(64, 127.99999999999872, 2) is None
    with pytest.raises(ValueError, match="must equal"):
        eht.validate_geometry_contract(64, 128.000001, 2)
    with pytest.raises(ValueError, match="must equal"):
        eht.validate_geometry_contract(2, 1e-200, 1e-200)


@pytest.mark.parametrize("field", ["I_data", "Q_data", "U_data", "V_data"])
def test_raw_shape_mismatch_precedes_calculation_and_products(tmp_path, monkeypatch, field):
    raw = {key: np.ones((8, 8, 3)) for key in ("I_data", "Q_data", "U_data", "V_data")}
    raw["meta_data"] = []
    raw[field] = np.ones((6, 6, 3))
    monkeypatch.setattr(eht, "input_eht_file", lambda root, day: root / f"{day}.pkl")
    monkeypatch.setattr(eht, "_load_day", lambda path: raw)

    def forbidden(*args, **kwargs):
        pytest.fail("shape mismatch reached calculation")

    for name in ("align_metadata", "transform_stokes", "open_memmap", "summarize_evpa"):
        monkeypatch.setattr(eht, name, forbidden)
    with pytest.raises(ValueError, match=rf"{field} raw shape .*expected_shape"):
        eht.run_eht_stage(tmp_path / "in", tmp_path / "out", trust_pickle=True,
                          expected_shape=8, angular_fov_uas=16, pixel_size_uas=2)
    assert not list((tmp_path / "out").iterdir())


@pytest.fixture
def tiny_inputs(tmp_path):
    root = tmp_path / "input"
    root.mkdir()
    rng = np.random.default_rng(19)
    for day in eht.DAYS:
        raw = {"meta_data": [dict(zip(eht.FIELDS, (i, 1, 0, 1, 1, 1, 1))) for i in range(5)]}
        theta = rng.normal(0.1, 0.05, (8, 8, 5))
        raw.update(I_data=np.ones(theta.shape), Q_data=np.cos(2 * theta), U_data=np.sin(2 * theta))
        with (root / f"230_{day}.pkl").open("wb") as stream:
            pickle.dump(raw, stream)
    return root


def test_declared_report_keeps_rg_separate_and_legacy_outputs_identical(tmp_path, tiny_inputs):
    offsets = {day: (0, 0) for day in eht.DAYS}
    options = dict(trust_pickle=True, topset_size=3, fov=35.13846177,
                   historical_mask=True, exclude_columns=1, exclusion_rectangles=((2, 3, 2, 4),),
                   min_radius=0, angular_bins=4, angle_roll=1, rank_offsets=offsets,
                   final_offsets=offsets, operator=csr_matrix(np.eye(4, 64)), batch_size=2)
    legacy = eht.run_eht_stage(tiny_inputs, tmp_path / "legacy", **options)
    report = eht.run_eht_stage(tiny_inputs, tmp_path / "declared", expected_shape=8,
                               angular_fov_uas=16, pixel_size_uas=2, **options)
    assert legacy["geometry"]["mode"] == "legacy"
    assert legacy["expected_shape"] is None
    assert report["geometry"] == dict(mode="declared", expected_shape=8, angular_fov_uas=16.0,
                                       pixel_size_uas=2.0, fov_rg=options["fov"])
    assert report["fov"] == options["fov"] != report["angular_fov_uas"]
    assert report["expected_shape"] == report["shape"] == 8
    assert report["pixel_size_uas"] == 2
    assert report["rank_offsets"] == report["final_offsets"] == offsets
    assert report["angle_roll"] == 1
    persisted = json.loads((tmp_path / "declared" / "report.json").read_text())
    assert persisted["geometry"] == report["geometry"]
    assert persisted["mask_options"]["exclusion_rectangles"] == [[2, 3, 2, 4]]
    mask = np.load(tmp_path / "declared" / "mask.npy")
    assert not mask[:, :1].any() and not mask[2:3, 2:4].any()
    for path in (tmp_path / "legacy").rglob("*"):
        other = tmp_path / "declared" / path.relative_to(tmp_path / "legacy")
        if path.suffix == ".npy":
            np.testing.assert_array_equal(np.load(path), np.load(other))
        elif path.suffix == ".npz":
            with np.load(path) as old, np.load(other) as new:
                assert old.files == new.files
                for key in old.files:
                    np.testing.assert_array_equal(old[key], new[key])


def test_final_pass_rechecks_declared_shape(tmp_path, tiny_inputs, monkeypatch):
    load_day = eht._load_day
    calls = 0

    def changed_input(path):
        nonlocal calls
        raw = load_day(path)
        calls += 1
        if calls == 5:
            for field in ("I_data", "Q_data", "U_data"):
                raw[field] = raw[field][:6, :6]
        return raw

    monkeypatch.setattr(eht, "_load_day", changed_input)
    offsets = {day: (0, 0) for day in eht.DAYS}
    with pytest.raises(ValueError, match="raw shape .*expected_shape"):
        eht.run_eht_stage(tiny_inputs, tmp_path / "out", trust_pickle=True, topset_size=3,
                          expected_shape=8, angular_fov_uas=16, pixel_size_uas=2,
                          historical_mask=False, min_radius=0, rank_offsets=offsets,
                          final_offsets=offsets, angular_bins=4, operator=csr_matrix(np.eye(4, 64)))
    assert calls == 5
    assert not (tmp_path / "out" / "report.json").exists()
    assert not list((tmp_path / "out").glob("*/stokes.npy"))
    assert not list((tmp_path / "out").glob(".eht-work-*"))

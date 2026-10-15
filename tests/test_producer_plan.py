"""Synthetic BHAC footers and CSVs only: no ray tracing or historical outputs."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
import pytest

from eht_axion.model import omega
from eht_axion.producer import plan_generation, read_production_config, read_snapshot_footer
from eht_axion.producer.plan import validate_radial_table

ROOT = Path(__file__).resolve().parents[1]


def write_snapshot(path, *, time=15010.0, spin=0.9375, gamma=4 / 3, nx=(2, 2, 2),
                   ndim=3, neq=13, truncate=False):
    nleafs, nw, nws = 1, 12, 0
    metadata = [gamma, 0.0025, 1.0, spin] + [0.0] * (neq - 4)
    # Include a minimal payload and forest; light preflight only inspects the tail.
    payload = bytes(nleafs * np.prod(nx) * nw * 8) + struct.pack("<i", 1)
    if truncate:
        payload = payload[:8]
    path.write_bytes(payload + struct.pack("<3i", *nx) + struct.pack(f"<{neq}d", *metadata)
                     + struct.pack("<8id", nleafs, 3, ndim, 3, nw, nws, neq, 12, time))


def write_table(path, alphas=(0.1, 0.3, 0.44)):
    np.savetxt(path, [[alpha, r, alpha * r] for alpha in alphas for r in (1.0, 10.0, 50.0)],
               delimiter=",")


@pytest.fixture
def inputs(tmp_path):
    config = read_production_config(ROOT / "configs/production.yaml")
    config["producer"]["snapshots"] = {"start": 1600, "stop": 1610, "step": 10}
    fluid = tmp_path / "fluid"
    fluid.mkdir()
    write_snapshot(fluid / "data001600.dat")
    write_snapshot(fluid / "data001610.dat", time=15011.0)
    tables = {"low": tmp_path / "low.csv", "high": tmp_path / "high.csv"}
    for path in tables.values():
        write_table(path)
    return config, fluid, tables


def test_task_contract_and_physical_time(inputs, tmp_path):
    config, fluid, tables = inputs
    build = tmp_path / "future-build"
    plan = plan_generation(config, fluid, tables, build)
    assert not build.exists()
    assert plan["schema_version"] == 1
    assert plan["verification"] == "unverified"
    assert plan["hash_inputs"] is False
    assert plan["counts"] == {"snapshots": 2, "baseline": 4, "response": 256,
                              "background": 2, "total": 260}
    assert [s["label"] for s in plan["snapshots"]] == [1600, 1610]
    assert [s["time_tg"] for s in plan["snapshots"]] == [15010.0, 15011.0]
    assert plan["requirements"]["label_time_offset_tg"] == 14850.0
    assert len(set(t["id"] for t in plan["tasks"])) == 260
    expected_keys = {
        "id", "label", "snapshot_name", "time_tg", "source_path", "source_sha256", "kind",
        "resolution", "alpha", "omega", "norm", "phase_index", "output_slot", "table_id",
        "table_path", "table_spin", "table_sha256", "parameters",
    }
    for task in plan["tasks"]:
        assert set(task) == expected_keys
        assert task["source_sha256"] is None
        assert task["source_path"] == str(fluid / task["snapshot_name"])
        assert task["time_tg"] == (15010 if task["label"] == 1600 else 15011)
        params = task["parameters"]
        assert len(params) == 19
        assert all(set(entry) == {"unit", "value"} for entry in params.values())
        assert params["MBH"] == {"unit": "(Msun)", "value": 6.5e9}
        assert params["DISTANCE"] == {"unit": "(kpc)", "value": 16800.0}
        assert params["CAM_SIZE_X"] == {"unit": "(Rg)", "value": 46.2261}
        assert params["IMG_WIDTH"]["value"] == task["resolution"]
        assert params["FREQ_MIN"]["value"] == 230e9
        assert params["MAX_LEVEL"]["value"] == 1
        if task["kind"] == "baseline":
            assert task["alpha"] is None
            assert task["omega"] == task["norm"] == task["phase_index"] == task["output_slot"] == 0
            assert task["table_id"] is task["table_path"] is task["table_sha256"] is None
            assert params["AXION_ALPHA"]["value"] == 0
            assert params["AXION_TABLE"]["value"] == "disabled.csv"
        else:
            assert task["omega"] == omega(task["alpha"], 2)
            assert task["omega"] != task["alpha"]
            assert task["norm"] == 0.1
            assert task["output_slot"] == task["phase_index"] + 1
            assert task["table_id"] == ("low" if task["alpha"] <= 0.3 else "high")
            assert params["AXION_TABLE"]["value"] == str(tables[task["table_id"]])
    response = [t for t in plan["tasks"] if t["kind"] == "response" and t["label"] == 1600]
    assert {t["phase_index"] for t in response} == set(range(16))
    assert {t["output_slot"] for t in response} == set(range(1, 17))
    assert json.loads(json.dumps(plan, allow_nan=False)) == plan


def test_default_plan_never_hashes_or_writes(inputs, monkeypatch):
    config, fluid, tables = inputs
    before = {p: p.stat().st_mtime_ns for p in fluid.parent.rglob("*")}
    def forbidden(*args, **kwargs):
        pytest.fail("light plan must not hash full inputs")
    monkeypatch.setattr("eht_axion.producer.plan.file_sha256", forbidden)
    plan_generation(config, fluid, tables)
    assert {p: p.stat().st_mtime_ns for p in fluid.parent.rglob("*")} == before


def test_explicit_hashes_once_per_unique_file(inputs, monkeypatch):
    from eht_axion.producer.plan import file_sha256
    config, fluid, tables = inputs
    calls = []
    def tracked(path):
        calls.append(path)
        return file_sha256(path)
    monkeypatch.setattr("eht_axion.producer.plan.file_sha256", tracked)
    light = plan_generation(config, fluid, tables)
    plan = plan_generation(config, fluid, tables, hash_inputs=True)
    assert len(calls) == len(set(calls)) == 4
    assert plan["verification"] == "input_hashes_only"
    assert [t["id"] for t in plan["tasks"]] == [t["id"] for t in light["tasks"]]
    for snapshot in plan["snapshots"]:
        expected = hashlib.sha256(Path(snapshot["source_path"]).read_bytes()).hexdigest()
        assert snapshot["source_sha256"] == expected
    assert all(t["source_sha256"] for t in plan["tasks"])
    assert all(t["table_sha256"] for t in plan["tasks"] if t["kind"] == "response")


def test_first_scope_rejects_incompatible_resolution(inputs):
    config, fluid, tables = inputs
    config["producer"]["response"]["resolution"] = 400
    with pytest.raises(ValueError, match="response.resolution=100"):
        plan_generation(config, fluid, tables)


def test_stable_ids_and_baseline_semantic_dependencies(inputs):
    config, fluid, tables = inputs
    first = plan_generation(config, fluid, tables)
    assert first == plan_generation(config, fluid, tables)
    changed = deepcopy(config)
    changed["inference"]["coupling_points"] = 121
    second = plan_generation(changed, fluid, tables)

    def baselines(plan):
        return [t["id"] for t in plan["tasks"] if t["kind"] == "baseline"]

    def responses(plan):
        return [t["id"] for t in plan["tasks"] if t["kind"] == "response"]
    assert first["config_sha256"] != second["config_sha256"]
    assert baselines(first) == baselines(second)
    assert responses(first) != responses(second)
    changed["image"]["initial_fov_rg"] = 46.3
    assert baselines(first) != baselines(plan_generation(changed, fluid, tables))


def test_table_change_keeps_baseline_identity(inputs):
    config, fluid, tables = inputs
    first = plan_generation(config, fluid, tables, hash_inputs=True)
    write_table(tables["low"], (0.1, 0.2, 0.44))
    second = plan_generation(config, fluid, tables, hash_inputs=True)
    assert first["radial_tables"]["low"]["sha256"] != second["radial_tables"]["low"]["sha256"]
    assert [t["id"] for t in first["tasks"] if t["kind"] == "baseline"] == [
        t["id"] for t in second["tasks"] if t["kind"] == "baseline"]


def test_no_lexicographic_slice_or_extra_file_reads(inputs):
    config, fluid, tables = inputs
    (fluid / "data001590.dat").write_text("bad old input, must not read")
    (fluid / "data001605.dat").write_text("off cadence, must not read")
    (fluid / "data003010.dat").write_text("outside inclusive range, must not read")
    assert plan_generation(config, fluid, tables)["counts"]["snapshots"] == 2


@pytest.mark.parametrize("mutation,match", [
    ({"time": 15010.0}, "cadence"), ({"time": 15012.0}, "cadence"),
    ({"time": float("nan")}, "Nonfinite"), ({"spin": 0.99}, "actual spin"),
    ({"gamma": 1.0}, "gamma"), ({"gamma": 1.5}, "gamma changes"),
    ({"nx": (5, 2, 2)}, "divide"), ({"ndim": 2}, "layout"),
    ({"truncate": True}, "Truncated BHAC payload"),
])
def test_bad_footer_or_cadence_fails(inputs, mutation, match):
    config, fluid, tables = inputs
    args = {"time": 15011.0, **mutation}
    write_snapshot(fluid / "data001610.dat", **args)
    with pytest.raises(ValueError, match=match):
        plan_generation(config, fluid, tables)


def test_truncated_and_missing_footer(inputs):
    config, fluid, tables = inputs
    (fluid / "data001610.dat").write_bytes(b"short")
    with pytest.raises(ValueError, match="Truncated BHAC footer"):
        plan_generation(config, fluid, tables)
    (fluid / "data001610.dat").unlink()
    with pytest.raises(FileNotFoundError):
        plan_generation(config, fluid, tables)


@pytest.mark.parametrize("content,match", [
    ("alpha,r,R\n0.1,1,1\n", "headerless"), ("# comment\n0.1,1,1\n", "headerless"),
    ("0.1,1,1,4\n0.1,2,2,4\n", "three columns"), ("0.1,1,1\n", "two rows"),
    ("0.1,1,nan\n0.1,2,2\n", "finite"), ("0.1,1,-1\n0.1,2,2\n", "R >= 0"),
    ("0.1,2,1\n0.1,1,2\n", "increasing"),
    ("0.1,1,1\n0.1,1,2\n", "increasing"),
    ("0.3,1,1\n0.3,2,2\n0.1,1,1\n0.1,2,2\n", "rectangular"),
    ("0.1,1,1\n0.1,2,2\n0.3,1,1\n0.3,3,2\n", "identical"),
    ("0.1,1,1\n0.1,2,2\n0.3,1,1\n", "rectangular"),
])
def test_invalid_radial_tables(tmp_path, content, match):
    path = tmp_path / "table.csv"
    path.write_text(content)
    with pytest.raises(ValueError, match=match):
        validate_radial_table(path, [])


def test_alpha_coverage(inputs):
    config, fluid, tables = inputs
    write_table(tables["high"], (0.1, 0.3))
    with pytest.raises(ValueError, match="cover required alphas"):
        plan_generation(config, fluid, tables)


def test_table_validation_once_per_unique_file_and_no_renormalization(inputs, monkeypatch):
    config, fluid, tables = inputs
    tables["high"] = tables["low"]
    calls = []
    def tracked(path, required):
        calls.append(path)
        return validate_radial_table(path, required)
    monkeypatch.setattr("eht_axion.producer.plan.validate_radial_table", tracked)
    before = tables["low"].read_bytes()
    plan = plan_generation(config, fluid, tables)
    assert len(calls) == 1
    assert before == tables["low"].read_bytes()
    assert plan["radial_tables"]["low"]["metadata"]["normalization"] == "unchanged"
    assert plan["radial_tables"]["high"]["metadata"]["provenance_verified"] is False


def test_paths_and_operational_arguments(inputs, tmp_path):
    config, fluid, tables = inputs
    with pytest.raises(ValueError, match="exactly low and high"):
        plan_generation(config, fluid, {"low": tables["low"]})
    with pytest.raises(ValueError, match="boolean"):
        plan_generation(config, fluid, tables, hash_inputs=1)
    with pytest.raises(ValueError, match="single C path token"):
        plan_generation(config, fluid, tables, tmp_path / "has spaces")
    with pytest.raises(ValueError, match="255"):
        plan_generation(config, fluid, tables, tmp_path / ("é" * 130))
    no_producer = {k: v for k, v in config.items() if k != "producer"}
    with pytest.raises(ValueError, match="requires producer"):
        plan_generation(no_producer, fluid, tables)


def test_resource_lower_bounds_and_official_setup(inputs):
    config, fluid, tables = inputs
    plan = plan_generation(config, fluid, tables)
    req = plan["requirements"]
    assert req["output_stokes_lower_bound_bytes"] == (258 * 100**2 + 2 * 400**2) * 4 * 8
    assert req["per_process_memory_lower_bound_bytes"] > 0
    assert req["live_rehash_required_at_execution"] is True
    assert req["free_space_check_owner"] == "caller"
    assert "--sf=sfc" in req["official_setup_args"]
    assert "--sigma-cut=1.0" in req["official_setup_args"]
    assert "--uniform=0" in req["official_setup_args"]
    metadata = read_snapshot_footer(fluid / "data001600.dat", config["producer"]["grid"])
    assert metadata["gamma"] == 4 / 3
    assert metadata["spin"] == 0.9375
    assert metadata["body_verified"] is False

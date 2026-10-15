"""Synthetic HDF delivery and downstream identity wiring, no real ray tracing."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from eht_axion.pipeline.common import stage_directory, verify_stage
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.layout import input_raptor_files
from eht_axion.producer.delivery import deliver_stokes, generated_inputs
from test_producer_outputs import write_products

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def delivered(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("generated-delivery")
    config = read_config(ROOT / "configs/production.yaml")
    config["data"]["alphas"] = [0.2]
    config["producer"]["alpha_table"] = {"0.20": "low"}
    config["producer"]["snapshots"] = {"start": 1600, "stop": 1610, "step": 10}
    run = tmp_path / "run"
    from test_producer_plan import write_snapshot, write_table
    from eht_axion.producer.plan import plan_generation
    from eht_axion.producer.jobs import _prepare_task
    from eht_axion.producer.outputs import validate_job_output, CORE_CONSTANTS

    fluid = tmp_path / "fluid"
    fluid.mkdir()
    for label in (1600, 1610):
        write_snapshot(fluid / f"data{label:06d}.dat", time=15010+(label-1600)/10)
    tables = {key: tmp_path / f"{key}.csv" for key in ("low", "high")}
    for table in tables.values():
        write_table(table)
    plan = plan_generation(config, fluid, tables, run / "build")
    plan["expected_core_commit"] = "synthetic-commit"
    build = run / "build"
    with stage_directory(build, "producer-build", {
        "source_commit": "synthetic-commit", "settings": config["producer"]["build"]
    }) as temp:
        (temp / "RAPTOR").write_text("synthetic binary, never executed")
        (temp / "RAPTOR").chmod(0o755)
        grid = config["producer"]["grid"]
        fields = {f"nxlone{i+1}": v for i, v in enumerate(grid["nxlone"])}
        fields.update({f"xprobminx{i+1}": v for i, v in enumerate(grid["xprobmin"])})
        fields.update({f"xprobmax{i+1}": v for i, v in enumerate(grid["xprobmax"])})
        fields["hslope"] = grid["hslope"]
        (temp / "grid_mks.in").write_text("".join(f"{k} = {v}\n" for k, v in fields.items()))
        (temp / "definitions.h").write_text(
            "".join(f"#define {key} ({value})\n" for key,value in CORE_CONSTANTS.items()))
        macros = {"POL":"1", "num_frequencies":"1", "num_pixels_1d":"10", "NSPIN":"3",
                  "SMR":"1", "AMR":"0", "UNIF":"0", "metric":"MKSBHAC", "TH":"1",
                  "DF":"TH", "IMGFILE":"1", "SPECFILE":"1", "int_method":"RK4",
                  "SFC":"1", "BPOL":"PLUS", "SIGMA_CUT":"1.0",
                  "THETAE_MIN":"0.001", "THETAE_MAX":"100.0"}
        (temp / "macros.txt").write_text(
            # Real preprocessor output intersperses empty include-guard macros.
            "".join(f"#define GUARD_{key}\n#define {key} ({value})\n" for key,value in macros.items()))
        (temp / "receipt.json").write_text(json.dumps({"source_commit":"synthetic-commit"}))
    base_inputs = [build / "manifest.json", build / "RAPTOR", build / "grid_mks.in", build / "definitions.h"]
    for task in plan["tasks"]:
        prepared, identity, inputs = _prepare_task(task)
        prepared["core_constants"] = CORE_CONSTANTS
        prepared["core_constants_source"] = str(build / "definitions.h")
        with stage_directory(run / "jobs" / task["id"], "producer-job", identity,
                             base_inputs + inputs) as temp:
            write_products(temp, prepared["parameters"], prepared["output_slot"])
            (temp / "success.json").write_text(json.dumps({
                "schema": 1, "returncode": 0, "timed_out": False, "output_validated": True}))
            (temp / "output.json").write_text(json.dumps(validate_job_output(temp, prepared)))
    products = deliver_stokes(plan, run)
    return config, run, products


def test_delivery_shapes_axes_labels_and_ancestry(delivered):
    config, run, products = delivered
    assert len(products) == 2
    for alpha, resolution, shape in ((.2, 100, (2, 17, 100, 100, 4)),
                                     (.3, 400, (2, 1, 400, 400, 4))):
        paths = input_raptor_files(run / "data", alpha, resolution)
        assert len(paths) == 1
        array = np.load(paths[0], mmap_mode="r")
        assert array.shape == shape
        expected = (1 + np.arange(resolution)[:, None]*2 + np.arange(resolution)[None, :]*3)*1e20
        np.testing.assert_array_equal(array[0, 0, ..., 0], expected)
        inputs, labels = generated_inputs(paths, config, alpha=alpha, resolution=resolution)
        assert labels == [1600, 1610]
        assert paths[0].parent / "manifest.json" in inputs
        record = verify_stage(paths[0].parent, recursive=True, allow_imported=False)
        assert all(".building-" not in e["path"] for e in record["inputs"])
        assert len([e for e in record["inputs"] if e["path"].endswith("manifest.json")]) == shape[0]*shape[1]


@pytest.mark.parametrize("key,value", [("mass_solar", 6.2e9), ("distance_mpc", 16.9)])
def test_generated_physics_cannot_be_silently_changed(delivered, key, value):
    config, run, _ = delivered
    changed = deepcopy(config)
    changed["physical"][key] = value
    paths = input_raptor_files(run / "data", .2, 100)
    with pytest.raises(ValueError, match="science identity"):
        generated_inputs(paths, changed, alpha=.2, resolution=100)


def test_generated_metadata_and_job_ancestors_live_verified(delivered):
    config, run, products = delivered
    path = products[0]
    metadata = path.parent / "generated.json"
    contents = metadata.read_text()
    changed = json.loads(contents)
    changed["omega"] = .2
    metadata.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="checksum"):
        generated_inputs([path], config, alpha=.2, resolution=100)
    metadata.write_text(contents)
    hdf = next((run / "jobs").glob("*/output/*.h5"))
    original = hdf.read_bytes()
    try:
        hdf.write_bytes(b"changed")
        with pytest.raises(ValueError, match="checksum|lineage"):
            generated_inputs([path], config, alpha=.2, resolution=100)
    finally:
        hdf.write_bytes(original)


def test_generated_data_reaches_real_preprocess_fit_and_background(delivered, tmp_path):
    import yaml
    from eht_axion.pipeline.workflow import run_stage
    from eht_axion.producer.delivery import validate_processed_source

    config, run, _ = delivered
    profile = deepcopy(config)
    profile["image"]["final_shape"] = 8
    profile["background"]["lag_frames"] = 1
    profile["background"]["method"] = "mle"
    path = tmp_path / "production.yaml"
    path.write_text(yaml.safe_dump(profile))
    output = tmp_path / "analysis"
    run_stage(path, output, raptor_input_root=run / "data", stage="raptor")
    fitted = output / "raptor/response/alpha-0.20"
    with np.load(fitted / "response-2d.npz") as z:
        np.testing.assert_array_equal(z["snapshot_labels"], [1600, 1610])
    record = verify_stage(output / "raptor/processed/alpha-0.20", recursive=True)
    validate_processed_source(record, profile, alpha=.2, resolution=100)
    # The actual background numerical path is tested without inventing EHT data:
    # it only needs a declared matching mask from the separate EHT source stage.
    with stage_directory(output / "eht", "eht:ensemble", {}) as temporary:
        np.save(temporary / "mask.npy", np.ones((8, 8), dtype=bool))
    run_stage(path, output, raptor_input_root=run / "data", stage="background")
    bg = verify_stage(output / "raptor/background", recursive=True)
    validate_processed_source(bg, profile, alpha=.3, resolution=400)
    assert np.load(output / "raptor/background/stokes.npy", mmap_mode="r").shape == (2, 1, 8, 8, 4)


def test_missing_metadata_cannot_downgrade_to_legacy(delivered):
    config, run, products = delivered
    path = products[0]
    metadata = path.parent / "generated.json"
    content = metadata.read_bytes()
    metadata.unlink()
    try:
        legacy = {k: v for k, v in config.items() if k != "producer"}
        with pytest.raises(ValueError, match="legacy downgrade"):
            generated_inputs([path], legacy, alpha=.2, resolution=100)
    finally:
        metadata.write_bytes(content)


def test_convert_only_rejects_wrong_build_even_with_valid_jobs(delivered):
    config, run, _ = delivered
    from eht_axion.producer.plan import plan_generation

    source = run.parent / "fluid"
    tables = {key: run.parent / f"{key}.csv" for key in ("low", "high")}
    plan = plan_generation(config, source, tables, run / "build")
    plan["expected_core_commit"] = "other-commit"
    with pytest.raises(ValueError, match="Build source commit"):
        deliver_stokes(plan, run)


def test_production_never_falls_back_to_legacy(tmp_path):
    config = read_config(ROOT / "configs/production.yaml")
    cube = tmp_path / "stokes.npy"
    np.save(cube, np.zeros((2, 17, 100, 100, 4)))
    with pytest.raises(ValueError, match="requires verified generated"):
        generated_inputs([cube], config, alpha=.2, resolution=100)
    legacy = {k: v for k, v in config.items() if k != "producer"}
    assert generated_inputs([cube], legacy, alpha=.2, resolution=100) == ([cube], None)

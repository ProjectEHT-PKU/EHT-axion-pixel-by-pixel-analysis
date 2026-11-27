"""Synthetic HDF image cubes and downstream identity wiring, no real ray tracing."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from eht_axion.experiment import read_experiment
from eht_axion.provenance import stage_directory, verify_stage
from eht_axion.raptor.layout import image_cube
from eht_axion.raptor.images import write_images, generated_inputs, input_cadence_hours
from test_raptor_outputs import write_products

REPOSITORY = "github.com/example/raptor"

@pytest.fixture(scope="module")
def raptor_images(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("generated-images")
    config = read_experiment()
    config["shared"]["alphas"] = [0.2]
    config["raptor"]["alpha_table"] = {"0.20": "low"}
    config["raptor"]["snapshots"] = {"start": 1600, "stop": 1610, "step": 10}
    run = tmp_path / "run"
    from test_raptor_plan import write_snapshot, write_table
    from eht_axion.raptor.plan import plan_generation
    from eht_axion.raptor.jobs import _prepare_task
    from eht_axion.raptor.outputs import validate_job_output, CORE_CONSTANTS

    fluid = tmp_path / "fluid"
    fluid.mkdir()
    for label in (1600, 1610):
        write_snapshot(fluid / f"data{label:06d}.dat", time=15010+(label-1600)/10)
    tables = {key: tmp_path / f"{key}.csv" for key in ("low", "high")}
    for table in tables.values():
        write_table(table)
    plan = plan_generation(config, fluid, tables, run / "build")
    plan["raptor_source"] = {"repository": REPOSITORY, "commit": "synthetic-commit"}
    build = run / "build"
    with stage_directory(build, "raptor:build", {
        "source_repository": REPOSITORY, "source_commit": "synthetic-commit",
        "settings": config["raptor"]["build"]
    }) as temp:
        (temp / "RAPTOR").write_text("synthetic binary, never executed")
        (temp / "RAPTOR").chmod(0o755)
        grid = config["raptor"]["grid"]
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
        (temp / "receipt.json").write_text(json.dumps({"source_repository": REPOSITORY,
                                                       "source_commit": "synthetic-commit"}))
    base_inputs = [build / "manifest.json", build / "RAPTOR", build / "grid_mks.in", build / "definitions.h"]
    for task in plan["tasks"]:
        prepared, identity, inputs = _prepare_task(task)
        prepared["core_constants"] = CORE_CONSTANTS
        prepared["core_constants_source"] = str(build / "definitions.h")
        with stage_directory(run / "jobs" / task["id"], "raptor:job", identity,
                             base_inputs + inputs) as temp:
            write_products(temp, prepared["parameters"], prepared["output_slot"])
            validate_job_output(temp, prepared)
    products = write_images(plan, run)
    return config, run, products


def test_image_shapes_axes_labels_and_ancestry(raptor_images):
    config, run, products = raptor_images
    assert len(products) == 2
    for alpha, resolution, shape in ((.2, 100, (2, 17, 100, 100, 4)),
                                     (None, 400, (2, 1, 400, 400, 4))):
        path = image_cube(run / "images", alpha, resolution)
        array = np.load(path, mmap_mode="r")
        assert array.shape == shape
        expected = (1 + np.arange(resolution)[:, None]*2 + np.arange(resolution)[None, :]*3)*1e20
        np.testing.assert_array_equal(array[0, 0, ..., 0], expected)
        inputs, labels = generated_inputs(path, config, alpha=alpha, resolution=resolution)
        assert labels == [1600, 1610]
        assert path.parent / "manifest.json" in inputs
        record = verify_stage(path.parent)
        assert all(".building-" not in e for e in record["inputs"])
        assert len([e for e in record["inputs"] if e.endswith("manifest.json")]) == shape[0]*shape[1]


@pytest.mark.parametrize("key,value", [("mass_solar", 6.2e9), ("distance_mpc", 16.9)])
def test_generated_physics_cannot_be_silently_changed(raptor_images, key, value):
    config, run, _ = raptor_images
    changed = deepcopy(config)
    changed["shared"][key] = value
    paths = image_cube(run / "images", .2, 100)
    with pytest.raises(ValueError, match="identity differs"):
        generated_inputs(paths, changed, alpha=.2, resolution=100)


def test_analysis_ignores_sky_grid_but_not_ray_tracing_changes(raptor_images):
    config, run, _ = raptor_images
    changed = deepcopy(config)
    changed["shared"]["image"]["fov_uas"] = 100.0
    changed["analysis"]["inference"]["phase_count"] = 32
    paths = image_cube(run / "images", .2, 100)
    assert generated_inputs(paths, changed, alpha=.2, resolution=100)[1] == [1600, 1610]
    changed["raptor"]["plasma"]["mass_unit_g"] = 4e25
    with pytest.raises(ValueError, match="identity differs"):
        generated_inputs(paths, changed, alpha=.2, resolution=100)


def test_generated_data_reaches_the_expected_signal_and_background(raptor_images, tmp_path):
    import yaml

    from scripts.analysis.expected_signal import run as run_expected
    from scripts.analysis.constraints import validate_expected_signal

    config, run, _ = raptor_images
    profile = deepcopy(config)
    profile["shared"]["image"]["shape"] = 8
    profile["analysis"]["background"]["lag_frames"] = 1
    path = tmp_path / "experiment.yaml"
    path.write_text(yaml.safe_dump(profile))
    output = tmp_path / "expected-signal"
    run_expected(path, run / "images", output)
    record = verify_stage(output)
    assert record["stage"] == "analysis:expected-signal"
    assert record["config"]["snapshot_labels"] == [1600, 1610]
    assert sorted(p.name for p in output.iterdir()) == ["background.npy", "expected.npz", "manifest.json"]
    products, inputs = validate_expected_signal(output, profile)
    assert inputs == [output / "manifest.json", output / "expected.npz"]
    assert products["2d"][3].shape == (8, 8) and products["1d"][3].shape == (360,)
    assert np.load(output / "background.npy", mmap_mode="r").shape == (2, 1, 8, 8, 4)
    (cube,) = [Path(e) for e in record["inputs"] if e.endswith("background/stokes.npy")]
    assert input_cadence_hours(cube, profile) > 0
    with pytest.raises(FileExistsError):
        run_expected(path, run / "images", output)
    changed = deepcopy(profile)
    changed["analysis"]["rotation_degrees"] += 1
    with pytest.raises(ValueError, match="geometry"):
        validate_expected_signal(output, changed)


def test_missing_metadata_is_rejected(raptor_images):
    config, run, products = raptor_images
    path = products[0]
    metadata = path.parent / "generated.json"
    content = metadata.read_bytes()
    metadata.unlink()
    try:
        with pytest.raises(ValueError, match="Missing stage output"):
            generated_inputs(path, config, alpha=.2, resolution=100)
    finally:
        metadata.write_bytes(content)


def test_convert_only_rejects_wrong_build_even_with_valid_jobs(raptor_images):
    config, run, _ = raptor_images
    from eht_axion.raptor.plan import plan_generation

    source = run.parent / "fluid"
    tables = {key: run.parent / f"{key}.csv" for key in ("low", "high")}
    plan = plan_generation(config, source, tables, run / "build")
    plan["raptor_source"] = {"repository": REPOSITORY, "commit": "other-commit"}
    with pytest.raises(ValueError, match="Build RAPTOR source"):
        write_images(plan, run)


def test_unimage_cube_is_rejected(tmp_path):
    cube = tmp_path / "stokes.npy"
    np.save(cube, np.zeros((2, 17, 100, 100, 4)))
    with pytest.raises(ValueError, match="no manifest"):
        generated_inputs(cube, read_experiment(), alpha=.2, resolution=100)

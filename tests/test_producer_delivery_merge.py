"""Synthetic multi-run delivery: real stage receipts, no compiler or ray tracing."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from eht_axion.experiment import read_experiment
from eht_axion.pipeline import common
from eht_axion.pipeline.layout import input_raptor_files
from eht_axion.producer.delivery import deliver_stokes, generated_inputs, validate_processed_source
from eht_axion.producer.jobs import _prepare_task, serialize_parameters
from eht_axion.producer.outputs import CORE_CONSTANTS, validate_job_output
from eht_axion.producer.plan import plan_generation
from test_producer_outputs import write_products
from test_producer_plan import write_snapshot, write_table

def _write_build(root, config, compiler, *, settings=None):
    """Two scientifically equivalent builds may have different compiler/binary bytes."""
    build = root / "build"
    identity = {"source_commit": "fake", "settings": config["raptor"]["build"]
                if settings is None else settings}
    with common.stage_directory(build, "producer-build", identity) as temporary:
        # Executable permission is part of the build contract; this is never launched.
        (temporary / "RAPTOR").write_text(f"synthetic RAPTOR built by {compiler}; never execute\n")
        (temporary / "RAPTOR").chmod(0o755)
        grid = config["raptor"]["grid"]
        fields = {f"nxlone{i+1}": value for i, value in enumerate(grid["nxlone"])}
        fields.update({f"xprobminx{i+1}": value for i, value in enumerate(grid["xprobmin"])})
        fields.update({f"xprobmax{i+1}": value for i, value in enumerate(grid["xprobmax"])})
        fields["hslope"] = grid["hslope"]
        (temporary / "grid_mks.in").write_text("".join(f"{key} = {value}\n" for key, value in fields.items()))
        (temporary / "definitions.h").write_text(
            "".join(f"#define {key} ({value})\n" for key, value in CORE_CONSTANTS.items()))
        macros = {"POL": "1", "num_frequencies": "1", "num_pixels_1d": "10", "NSPIN": "3",
                  "SMR": "1", "AMR": "0", "UNIF": "0", "metric": "MKSBHAC", "TH": "1",
                  "DF": "TH", "IMGFILE": "1", "SPECFILE": "1", "int_method": "RK4",
                  "SFC": "1", "BPOL": "PLUS", "SIGMA_CUT": "1.0",
                  "THETAE_MIN": "0.001", "THETAE_MAX": "100.0"}
        (temporary / "macros.txt").write_text(
            "".join(f"#define GUARD_{key}\n#define {key} ({value})\n" for key, value in macros.items()))
        (temporary / "receipt.json").write_text(json.dumps({
            "source_commit": "fake", "toolchain": {"compiler": compiler}}))
    return [build / name for name in ("manifest.json", "RAPTOR", "grid_mks.in", "definitions.h")]


def _write_job(root, task, build_inputs, *, extra_inputs=(), omit_binary=False):
    prepared, identity, raw_inputs = _prepare_task(task)
    prepared.update(core_constants=CORE_CONSTANTS, core_constants_source=str(build_inputs[-1]))
    inputs = [path for path in build_inputs if not (omit_binary and path.name == "RAPTOR")]
    with common.stage_directory(root / "jobs" / task["id"], "producer-job", identity,
                                inputs + raw_inputs + list(extra_inputs)) as temporary:
        write_products(temporary, prepared["parameters"], prepared["output_slot"])
        (temporary / "task.json").write_text(json.dumps(identity))
        (temporary / "model.in").write_text(serialize_parameters(prepared["parameters"]))
        shutil.copyfile(build_inputs[-2], temporary / "grid_mks.in")
        (temporary / "success.json").write_text(json.dumps({
            "schema": 1, "returncode": 0, "timed_out": False, "output_validated": True}))
        (temporary / "output.json").write_text(json.dumps(validate_job_output(temporary, prepared)))


def _unchanged_records(roots, raw_inputs):
    """Keep exact bytes and mtimes, so an equivalent re-sign is not mistaken for preservation."""
    paths = set(raw_inputs)
    for root in roots:
        for name in ("manifest.json", "receipt.json", "success.json", "output.json", "RAPTOR"):
            paths.update(root.rglob(name))
    return {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


@pytest.fixture(scope="module")
def merged(tmp_path_factory):
    # Like the build unit fixture, isolate identity from concurrent checkout edits.
    # Do not test signature-cache reuse while delivery.py is under development.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(common, "revision", lambda stage=None: {
            "git_commit": "merge-fixture", "python_source_sha256": "merge-fixture-implementation"})
        directory = tmp_path_factory.mktemp("merged-delivery")
        config = read_experiment()
        config["shared"]["alphas"] = [0.2]
        config["raptor"]["alpha_table"] = {"0.20": "low"}
        config["raptor"]["snapshots"] = {"start": 1600, "stop": 1610, "step": 10}
        fluid = directory / "fluid"
        fluid.mkdir()
        for label in (1600, 1610):
            write_snapshot(fluid / f"data{label:06d}.dat", time=15010 + (label - 1600) / 10)
        table_root = directory / "tables"
        table_root.mkdir()
        tables = {key: table_root / f"{key}.csv" for key in ("low", "high")}
        for path in tables.values():
            write_table(path)
        roots = [directory / "source-one", directory / "source-two"]
        plan = plan_generation(config, fluid, tables, roots[0] / "build")
        plan["expected_core_commit"] = "fake"
        assert plan["counts"] == {"snapshots": 2, "baseline": 4, "response": 32,
                                  "background": 2, "total": 36}
        builds = [_write_build(root, config, compiler)
                  for root, compiler in zip(roots, ("compiler-one", "compiler-two"))]
        owners = {}
        for task in plan["tasks"]:
            index = int(task["label"] == 1610)
            owners[task["id"]] = roots[index]
            _write_job(roots[index], task, builds[index])
        raw_inputs = list(fluid.iterdir()) + list(tables.values())
        before = _unchanged_records(roots, raw_inputs)
        output = directory / "new-delivery"
        # run_root deliberately has no build/jobs: only job_roots own the sources.
        run = directory / "unused-run"
        products = deliver_stokes(plan, run, job_roots=roots, output_root=output)
        fixture = SimpleNamespace(config=config, plan=plan, roots=roots, builds=builds,
                                  owners=owners, raw_inputs=raw_inputs, before=before,
                                  run=run, output=output, products=products, fluid=fluid, tables=table_root)
        assert _unchanged_records(roots, raw_inputs) == before
        yield fixture
        assert _unchanged_records(roots, raw_inputs) == before


def test_merged_shapes_canonical_layout_and_labels(merged):
    m = merged
    expected_paths = [
        m.output / "data/raptor/res-100/alpha-0.20/stokes.npy",
        m.output / "data/raptor/res-400/alpha-0.30/t1600-t1610/stokes-t1600-t1610.npy",
    ]
    assert m.products == expected_paths
    assert not m.run.exists()
    for root in m.roots:
        assert not (root / "data").exists()
    for alpha, resolution, states, expected_path in zip((.2, .3), (100, 400), (17, 1), expected_paths):
        paths = input_raptor_files(m.output / "data", alpha, resolution)
        assert paths == [expected_path]
        cube = np.load(paths[0], mmap_mode="r", allow_pickle=False)
        assert cube.shape == (2, states, resolution, resolution, 4)
        expected = (1 + np.arange(resolution)[:, None] * 2 + np.arange(resolution)[None, :] * 3) * 1e20
        for time in range(2):
            for state in range(states):
                for channel, factor in enumerate((1, .2, -.3, .01)):
                    np.testing.assert_allclose(cube[time, state, ..., channel], expected * factor)
        inputs, labels = generated_inputs(paths, m.config, alpha=alpha, resolution=resolution)
        assert labels == [1600, 1610]
        assert paths[0].parent / "manifest.json" in inputs
        metadata = json.loads((paths[0].parent / "generated.json").read_text())
        assert metadata["times_tg"] == [15010., 15011.]
        assert metadata["axes"] == ["time", "state", "x", "y", "stokes"]


def test_merged_sources_and_recursive_ancestry_preserve_both_actual_builds(merged):
    m = merged
    assert m.builds[0][1].read_bytes() != m.builds[1][1].read_bytes()
    for name in ("grid_mks.in", "definitions.h", "macros.txt"):
        assert (m.roots[0] / "build" / name).read_bytes() == (m.roots[1] / "build" / name).read_bytes()
    build_records = [common.verify_stage(root / "build", recursive=True)
                     for root in m.roots]
    assert build_records[0]["config"] == build_records[1]["config"]
    for path, resolution, count in zip(m.products, (100, 400), (34, 2)):
        record = common.verify_stage(path.parent, recursive=True)
        sources = json.loads((path.parent / "sources.json").read_text())
        selected = [task for task in m.plan["tasks"] if task["resolution"] == resolution]
        assert sources == [{"task_id": task["id"], "root": str(m.owners[task["id"]])}
                           for task in selected]
        assert len(sources) == count
        assert {item["root"] for item in sources} == {str(root) for root in m.roots}
        assert "sources.json" in {entry["path"] for entry in record["outputs"]}
        expected_inputs = set()
        for task in selected:
            job = m.owners[task["id"]] / "jobs" / task["id"]
            expected_inputs.update((job / "manifest.json",
                                    job / "output" / f"img_data_{task['output_slot']}.h5"))
            job_record = common.verify_stage(job, recursive=True)
            _, _, raw = _prepare_task(task)
            build = m.builds[m.roots.index(m.owners[task["id"]])]
            assert {Path(entry["path"]) for entry in job_record["inputs"]} == set(build + raw)
        assert {Path(entry["path"]) for entry in record["inputs"]} == expected_inputs
        closure = set(common.source_closure([path.parent / "manifest.json"]))
        assert {p for p in closure if p.name == "RAPTOR"} == {build[1] for build in m.builds}
        assert {build[0] for build in m.builds} <= closure
        assert all(".building-" not in str(p) for p in closure)
        assert ".building-" not in (path.parent / "sources.json").read_text()
        assert ".building-" not in (path.parent / "manifest.json").read_text()
    assert _unchanged_records(m.roots, m.raw_inputs) == m.before


def test_merged_background_reaches_real_downstream_path(merged, tmp_path):
    from analysis_templates import run as run_templates

    profile = deepcopy(merged.config)
    profile["shared"]["image"]["shape"] = 8
    profile["analysis"]["background"]["lag_frames"] = 1
    config_path = tmp_path / "experiment.yaml"
    config_path.write_text(yaml.safe_dump(profile))
    analysis = tmp_path / "analysis"
    run_templates(config_path, merged.output / "data", analysis)
    processed = analysis / "raptor/background"
    record = common.verify_stage(processed, recursive=True)
    validate_processed_source(record, profile, alpha=.3, resolution=400)
    assert np.load(processed / "stokes.npy", mmap_mode="r").shape == (2, 1, 8, 8, 4)
    closure = set(common.source_closure([processed / "manifest.json"]))
    assert {p for p in closure if p.name == "RAPTOR"} == {build[1] for build in merged.builds}


def _copy_root(merged, index, tmp_path):
    """Copies keep original signed edges; never repair/re-sign a source record."""
    root = tmp_path / f"source-copy-{index}"
    shutil.copytree(merged.roots[index], root)
    roots = list(merged.roots)
    roots[index] = root
    return root, roots


@pytest.mark.parametrize("case", ["ambiguous", "missing"])
def test_merge_requires_exactly_one_task_directory(merged, tmp_path, case):
    copied, roots = _copy_root(merged, 1, tmp_path)
    task = next(task for task in merged.plan["tasks"] if task["label"] == 1600)
    if case == "ambiguous":
        shutil.copytree(merged.roots[0] / "jobs" / task["id"], copied / "jobs" / task["id"])
        count = 2
    else:
        task = next(task for task in merged.plan["tasks"] if task["label"] == 1610)
        shutil.rmtree(copied / "jobs" / task["id"])
        count = 0
    target = tmp_path / "rejected"
    with pytest.raises(ValueError, match=f"exactly one source for task {task['id']}, found {count}"):
        deliver_stokes(merged.plan, merged.run, job_roots=roots, output_root=target)
    assert not target.exists()


@pytest.mark.parametrize("index", [0, 1])
def test_each_source_build_must_match_plan_settings(merged, tmp_path, index):
    copied, roots = _copy_root(merged, index, tmp_path)
    shutil.rmtree(copied / "build")
    settings = deepcopy(merged.config["raptor"]["build"])
    settings["sigma_cut"] = 2.0
    _write_build(copied, merged.config, "wrong-settings", settings=settings)
    common.verify_stage(copied / "build", recursive=True)
    target = tmp_path / "rejected"
    with pytest.raises(ValueError, match="Build settings differ from production plan"):
        deliver_stokes(merged.plan, merged.run, job_roots=roots, output_root=target)
    assert not target.exists()


@pytest.mark.parametrize("case", ["copied-build", "other-build", "extra-input", "missing-binary"])
def test_job_inputs_must_name_exact_actual_source_build(merged, tmp_path, case):
    copied, roots = _copy_root(merged, 1, tmp_path)
    if case != "copied-build":
        task = next(task for task in merged.plan["tasks"] if task["label"] == 1610)
        shutil.rmtree(copied / "jobs" / task["id"])
        build = merged.builds[0] if case == "other-build" else [
            copied / "build" / path.name for path in merged.builds[1]]
        extras = [merged.tables / "high.csv"] if case == "extra-input" else []
        _write_job(copied, task, build, extra_inputs=extras, omit_binary=case == "missing-binary")
    # These are complete valid stages, not checksum-corruption failures.
    for task in merged.plan["tasks"]:
        if task["label"] == 1610:
            common.verify_stage(copied / "jobs" / task["id"], recursive=True)
    target = tmp_path / "rejected"
    with pytest.raises(ValueError, match="Job lineage does not declare the selected build and source inputs"):
        deliver_stokes(merged.plan, merged.run, job_roots=roots, output_root=target)
    assert not target.exists()


@pytest.mark.parametrize("case", ["same-root", "inside-root", "ancestor", "fluid", "table"])
def test_delivery_output_must_not_overlap_sources(merged, tmp_path, case):
    targets = {"same-root": merged.roots[0], "inside-root": merged.roots[1] / "new-delivery",
               "ancestor": merged.roots[0].parent, "fluid": merged.fluid / "new-delivery",
               "table": merged.tables / "new-delivery"}
    target = targets[case]
    existed = target.exists()
    with pytest.raises(ValueError, match="disjoint from source job roots|overlaps a consumed source"):
        deliver_stokes(merged.plan, merged.run, job_roots=merged.roots, output_root=target)
    assert target.exists() == existed
    assert _unchanged_records(merged.roots, merged.raw_inputs) == merged.before

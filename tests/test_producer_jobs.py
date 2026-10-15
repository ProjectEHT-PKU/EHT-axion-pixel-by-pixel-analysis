"""Fake Python executable exercises orchestration, not scientific simulation."""
from copy import deepcopy
import inspect
import json
import os
from pathlib import Path
import sys

import pytest

from eht_axion.pipeline import common
from eht_axion.producer.jobs import PARAMETER_UNITS, run_jobs, serialize_parameters
from eht_axion.producer.outputs import CORE_CONSTANTS
from test_producer_outputs import make_task, write_products


@pytest.fixture
def setup(tmp_path):
    source = tmp_path / "fluid" / "fluid.dat"
    source.parent.mkdir()
    source.write_bytes(b"synthetic snapshot bytes, not a real fluid")
    table = tmp_path / "tables" / "table.csv"
    table.parent.mkdir()
    table.write_text(".2,1,1\n.2,100,0\n")
    mode = tmp_path / "mode"
    mode.write_text("success")
    script = f"#!{sys.executable}\n" + inspect.getsource(write_products)
    script += f'''
import json, os, sys, time
from pathlib import Path
mode = Path({str(mode)!r}).read_text().strip()
tokens = Path(sys.argv[1]).read_text().split()
p = {{}}
for i in range(0, len(tokens), 3):
    name, unit, value = tokens[i:i+3]
    p[name] = value if name == "AXION_TABLE" else (int(value) if name in
        {{"IMG_WIDTH", "IMG_HEIGHT", "AXION_PHASE", "MAX_LEVEL"}} else float(value))
print("fake stdout", flush=True)
print("fake stderr", file=sys.stderr, flush=True)
Path("environment.json").write_text(json.dumps({{k: os.environ[k] for k in
    ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]}}))
assert Path("grid_mks.in").read_text() == "private grid\\n"
assert Path(sys.argv[2]).read_bytes()
if mode == "timeout":
    import subprocess
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    Path({str(tmp_path / 'child.pid')!r}).write_text(str(child.pid))
    time.sleep(60)
write_products(Path.cwd(), p, int(sys.argv[3]))
if mode == "invalid":
    (Path("output") / f"img_data_{{sys.argv[3]}}.h5").write_bytes(b"invalid")
if mode == "mutate":
    Path(sys.argv[2]).write_bytes(b"changed during run")
if mode == "failure":
    sys.exit(7)
'''
    build = tmp_path / "build"
    with common.stage_directory(build, "producer-build", {"fake": True}) as temporary:
        binary = temporary / "RAPTOR"
        binary.write_text(script)
        binary.chmod(0o755)
        (temporary / "grid_mks.in").write_text("private grid\n")
        (temporary / "definitions.h").write_text(
            "".join(f"#define {k} ({v})\n" for k, v in CORE_CONSTANTS.items()))
    task = make_task(source)
    run = tmp_path / "run"
    return dict(source=source, table=table, mode=mode, build=build, task=task, run=run)


def invoke(setup, tasks=None, **kwargs):
    return run_jobs({"tasks": tasks or [setup["task"]]}, setup["build"], setup["run"], **kwargs)


def test_success_resume_stable_lineage_and_private_cwd(setup):
    result = invoke(setup, threads=3)
    assert result["counts"] == {"complete": 1, "cached": 0, "failed": 0}
    job = setup["run"] / "jobs/baseline"
    manifest = common.verify_stage(job, recursive=True)
    assert all(".building-" not in e["path"] for e in manifest["inputs"])
    assert str(setup["source"]) in {e["path"] for e in manifest["inputs"]}
    assert json.loads((job / "environment.json").read_text()) == {
        "OMP_NUM_THREADS": "3", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    setup["mode"].write_text("failure")
    again = invoke(setup, threads=3)
    assert again["counts"]["cached"] == 1
    assert len(list((setup["run"] / "attempts/baseline").iterdir())) == 1
    assert (job / "stdout.log").read_text().strip() == "fake stdout"
    assert "fake stderr" in (job / "stderr.log").read_text()


@pytest.mark.parametrize("mode,returncode", [("failure", 7), ("invalid", 0), ("timeout", -9)])
def test_failed_attempt_retained_and_resume_reruns(setup, mode, returncode):
    setup["mode"].write_text(mode)
    result = invoke(setup, timeout_seconds=.4 if mode == "timeout" else 10)
    assert result["counts"]["failed"] == 1
    assert not (setup["run"] / "jobs/baseline").exists()
    attempt = Path(result["jobs"][0]["attempt"])
    record = json.loads((attempt / "result.json").read_text())
    assert record["returncode"] == returncode
    assert record["complete"] is False
    assert (attempt / "work/model.in").is_file()
    assert (attempt / "stderr.log").is_file()
    setup["mode"].write_text("success")
    result = invoke(setup)
    assert result["counts"]["complete"] == 1
    assert len(list((setup["run"] / "attempts/baseline").iterdir())) == 2


def test_max_jobs_only_selects_first_n(setup):
    second = deepcopy(setup["task"])
    second["id"] = "second"
    result = invoke(setup, [setup["task"], second], max_jobs=1)
    assert result["planned_tasks"] == 2
    assert result["selected_tasks"] == 1
    assert not (setup["run"] / "jobs/second").exists()


def test_plan_explicit_units_and_response_outputslot(setup):
    task = make_task(setup["source"], phase=0)
    task["table_path"] = str(setup["table"])
    task["parameters"]["AXION_TABLE"] = str(setup["table"])
    task["parameters"] = {name: {"unit": PARAMETER_UNITS[name], "value": value}
                          for name, value in task["parameters"].items()}
    result = invoke(setup, [task])
    assert result["counts"]["complete"] == 1
    directory = setup["run"] / "jobs/response-0"
    assert (directory / "output/img_data_1.h5").is_file()
    metadata = json.loads((directory / "output.json").read_text())
    assert metadata["state_slot"] == 1
    assert metadata["physical_phase_index"] == 0
    assert str(setup["table"]) in {e["path"] for e in common.verify_stage(directory)["inputs"]}


def test_baseline_does_not_depend_on_tables(setup):
    task = deepcopy(setup["task"])
    task.update(table_id="low", table_path=str(setup["table"]), table_spin=.99)
    result = invoke(setup, [task])
    assert result["counts"]["complete"] == 1
    setup["table"].unlink()
    task.update(table_id="high", table_path="/does/not/exist", table_spin=.1)
    result = invoke(setup, [task])
    assert result["counts"]["cached"] == 1
    manifest = common.verify_stage(setup["run"] / "jobs/baseline")
    assert not any("table.csv" in e["path"] for e in manifest["inputs"])


def test_source_mutation_rejects_complete(setup):
    setup["mode"].write_text("mutate")
    with pytest.raises(ValueError, match="changed"):
        invoke(setup)
    assert not (setup["run"] / "jobs/baseline").exists()
    receipt = next((setup["run"] / "attempts/baseline").glob("*/result.json"))
    assert json.loads(receipt.read_text())["complete"] is False


def test_cross_session_source_bytes_verified(setup):
    invoke(setup)
    before = setup["source"].stat()
    setup["source"].write_bytes(b"x" * before.st_size)
    os.utime(setup["source"], ns=(before.st_atime_ns, before.st_mtime_ns))
    result = invoke(setup)
    assert result["counts"]["failed"] == 1
    assert "Stale" in result["jobs"][0]["error"]


def test_shared_source_read_once_per_invocation(setup, monkeypatch):
    original = Path.open
    count = 0

    def counted(path, *args, **kwargs):
        nonlocal count
        if path == setup["source"] and args and args[0] == "rb":
            count += 1
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", counted)
    tasks = [deepcopy(setup["task"]) for _ in range(3)]
    for i, task in enumerate(tasks):
        task["id"] = f"baseline-{i}"
    result = invoke(setup, tasks)
    assert result["counts"]["complete"] == 3
    assert count == 1
    invoke(setup, tasks)
    assert count == 2


def test_sixteen_response_phases_share_one_source_hash(setup, monkeypatch):
    original = Path.open
    reads = {setup["source"]: 0, setup["table"]: 0}

    def counted(path, *args, **kwargs):
        if path in reads and args and args[0] == "rb":
            reads[path] += 1
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", counted)
    tasks = [make_task(setup["source"], phase=phase) for phase in range(16)]
    for task in tasks:
        task["table_path"] = str(setup["table"])
        task["parameters"]["AXION_TABLE"] = str(setup["table"])
    result = invoke(setup, tasks)
    assert result["counts"]["complete"] == 16
    assert reads == {setup["source"]: 1, setup["table"]: 1}
    result = invoke(setup, tasks)
    assert result["counts"]["cached"] == 16
    assert reads == {setup["source"]: 2, setup["table"]: 2}


def test_task_lock_prevents_launch(setup):
    jobroot = setup["run"] / "jobs"
    jobroot.mkdir(parents=True)
    (jobroot / "baseline.lock").write_text("owner")
    result = invoke(setup)
    assert result["counts"]["failed"] == 1
    assert "locked" in result["jobs"][0]["error"]
    assert not (setup["run"] / "attempts").exists()


@pytest.mark.parametrize("root", ["source-parent", "build", "linked"])
def test_unsafe_roots_fail_before_launch(setup, root):
    if root == "source-parent":
        setup["run"] = setup["source"].parent
    elif root == "build":
        setup["run"] = setup["build"] / "inside"
    else:
        target = setup["source"].parent / "target"
        target.mkdir()
        setup["run"].symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError):
        invoke(setup)


def test_failure_does_not_rollback_other_jobs(setup, monkeypatch):
    from eht_axion.producer import jobs
    original = jobs._execute
    calls = 0

    def execute(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            setup["mode"].write_text("failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(jobs, "_execute", execute)
    second = deepcopy(setup["task"])
    second["id"] = "second"
    result = invoke(setup, [setup["task"], second])
    assert result["counts"] == {"complete": 1, "cached": 0, "failed": 1}
    common.verify_stage(setup["run"] / "jobs/baseline", recursive=True)
    assert not (setup["run"] / "jobs/second").exists()


def test_build_child_of_run_container_is_allowed(setup):
    setup["run"].mkdir()
    relocated = setup["run"] / "build"
    setup["build"].rename(relocated)
    setup["build"] = relocated
    assert invoke(setup)["counts"]["complete"] == 1


def test_cancel_kills_process_group_and_keeps_attempt(setup, monkeypatch):
    from eht_axion.producer import jobs
    called = []
    original_kill = jobs._kill_group

    def kill(process):
        called.append(process.pid)
        original_kill(process)

    monkeypatch.setattr(jobs, "_kill_group", kill)
    monkeypatch.setattr(jobs.subprocess.Popen, "wait", _cancel_wait(jobs.subprocess.Popen.wait))
    with pytest.raises(KeyboardInterrupt):
        invoke(setup)
    assert called
    receipt = next((setup["run"] / "attempts/baseline").glob("*/result.json"))
    assert json.loads(receipt.read_text())["cancelled"] is True
    assert not (setup["run"] / "jobs/baseline").exists()


def _cancel_wait(original):
    interrupted = False

    def wait(process, *args, **kwargs):
        nonlocal interrupted
        if process.args and str(process.args[0]).endswith("RAPTOR") and not interrupted:
            interrupted = True
            raise KeyboardInterrupt
        return original(process, *args, **kwargs)

    return wait


@pytest.fixture
def production(setup):
    from eht_axion.producer.config import read_production_config
    from eht_axion.producer.plan import plan_generation
    from test_producer_plan import write_snapshot, write_table

    config = read_production_config(Path(__file__).resolve().parents[1] / "configs/production.yaml")
    config["producer"]["snapshots"] = {"start": 1600, "stop": 1600, "step": 10}
    source = setup["source"].with_name("data001600.dat")
    write_snapshot(source)
    write_table(setup["table"])
    tables = {"low": setup["table"], "high": setup["table"]}
    plan = plan_generation(config, source.parent, tables)
    plan["expected_core_commit"] = "synthetic-core-commit"
    production_build = setup["build"].with_name("production-build")
    settings = {"settings": config["producer"]["build"], "source_commit": "synthetic-core-commit"}
    with common.stage_directory(production_build, "producer-build", settings) as work:
        for name in ("RAPTOR", "definitions.h"):
            import shutil
            shutil.copy2(setup["build"] / name, work / name)
        grid = config["producer"]["grid"]
        text = "".join(f"nxlone{i+1} = {v}\n" for i, v in enumerate(grid["nxlone"]))
        text += "".join(f"xprobminx{i+1} = {v}\n" for i, v in enumerate(grid["xprobmin"]))
        text += "".join(f"xprobmax{i+1} = {v}\n" for i, v in enumerate(grid["xprobmax"]))
        (work / "grid_mks.in").write_text(text + f"hslope = {grid['hslope']}\n")
        (work / "macros.txt").write_text("\n".join(f"#define {k} {v}" for k, v in {
            "POL": "1", "num_frequencies": "1", "num_pixels_1d": "10", "NSPIN": "3",
            "SMR": "1", "AMR": "0", "UNIF": "0", "metric": "MKSBHAC", "DF": "TH",
            "IMGFILE": "1", "SPECFILE": "1", "int_method": "RK4", "SFC": "1", "BPOL": "PLUS",
            "SIGMA_CUT": str(config["producer"]["build"]["sigma_cut"]),
            "THETAE_MIN": str(config["producer"]["build"]["thetae_min"]),
            "THETAE_MAX": str(config["producer"]["build"]["thetae_max"]),
        }.items()))
        (work / "receipt.json").write_text(json.dumps({"source_commit": "synthetic-core-commit"}))
    setup["build"] = production_build
    return setup, plan


@pytest.mark.parametrize("product,content", [("grid_mks.in", "hslope = 99\n"),
                                             ("macros.txt", "#define POL 0\n"),
                                             ("receipt.json", '{"source_commit": "wrong"}')])
def test_production_rejects_verified_but_wrong_build_products(production, product, content):
    setup, plan = production
    build = setup["build"]
    (build / product).write_text(content)
    # Another legitimately published build can have internally verified bytes
    # without matching this plan. Rehash the fixture to simulate that situation.
    path = build / "manifest.json"
    record = json.loads(path.read_text())
    for entry in record["outputs"]:
        if entry["path"] == product:
            entry.update(size=(build / product).stat().st_size, sha256=common.sha256(build / product))
    path.write_text(json.dumps(record))
    common.verify_stage(build, recursive=True)
    with pytest.raises(ValueError):
        run_jobs(plan, build, setup["run"], max_jobs=1)
    assert not setup["run"].exists()


def test_production_gate_accepts_bound_plan(production, monkeypatch):
    from eht_axion.producer import jobs
    setup, plan = production
    executed = []

    def execute(command, cwd, env, stdout, stderr, timeout):
        executed.append(command)
        tokens = (cwd / "model.in").read_text().split()
        p = {}
        for i in range(0, len(tokens), 3):
            name, _, value = tokens[i:i+3]
            p[name] = value if name == "AXION_TABLE" else float(value)
        for name in ("IMG_WIDTH", "IMG_HEIGHT", "AXION_PHASE", "MAX_LEVEL"):
            p[name] = int(p[name])
        write_products(cwd, p, int(command[3]))
        stdout.write_text("synthetic\n")
        stderr.write_text("")
        return 0, False

    monkeypatch.setattr(jobs, "_execute", execute)
    result = run_jobs(plan, setup["build"], setup["run"], max_jobs=1)
    assert result["counts"]["complete"] == 1
    assert len(executed) == 1


@pytest.mark.parametrize("change", ["commit", "missing-commit", "build", "time", "mass", "omega",
                                    "table-spin", "source-hash"])
def test_production_gate_rejects_semantic_drift_before_launch(production, change):
    setup, plan = production
    task = plan["tasks"][0]
    if change == "commit":
        plan["expected_core_commit"] = "wrong"
    elif change == "missing-commit":
        del plan["expected_core_commit"]
    elif change == "build":
        plan["config"]["producer"]["build"]["sigma_cut"] = 2.
    elif change == "time":
        task["time_tg"] += 1
    elif change == "mass":
        task["parameters"]["MBH"]["value"] *= 2
    elif change == "source-hash":
        task["source_sha256"] = "f" * 64
    else:
        task = next(t for t in plan["tasks"] if t["kind"] == "response")
        plan["tasks"] = [task]
        if change == "omega":
            task["omega"] *= 2
            task["parameters"]["AXION_OMEGA"]["value"] = task["omega"]
        else:
            task["table_spin"] = .1
    with pytest.raises(ValueError):
        run_jobs(plan, setup["build"], setup["run"], max_jobs=1)
    assert not setup["run"].exists()


@pytest.mark.parametrize("mutation", ["extra", "missing", "unit", "phase", "integer", "path"])
def test_exact_parameter_contract(setup, mutation):
    p = deepcopy(setup["task"]["parameters"])
    if mutation == "extra":
        p["UNKNOWN"] = 0
    elif mutation == "missing":
        del p["AXION_ALPHA"]
    elif mutation == "unit":
        p["DISTANCE"] = {"unit": "(Mpc)", "value": 16.8}
    elif mutation == "phase":
        p["AXION_PHASE"] = 16
    elif mutation == "integer":
        p["IMG_WIDTH"] = 20.0
    else:
        p["AXION_TABLE"] = "has whitespace"
    with pytest.raises(ValueError):
        serialize_parameters(p)

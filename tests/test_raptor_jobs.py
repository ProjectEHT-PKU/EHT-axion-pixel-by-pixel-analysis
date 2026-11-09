"""Fake Python executable exercises orchestration, not scientific simulation."""
from copy import deepcopy
import inspect
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

import pytest

from eht_axion import provenance
from eht_axion.raptor.jobs import PARAMETER_UNITS, run_jobs, serialize_parameters
from eht_axion.raptor.outputs import CORE_CONSTANTS
from test_raptor_outputs import make_task, write_products


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
task_id = json.loads(Path("task.json").read_text())["id"]
Path("binding-environment.json").write_text(json.dumps({{k: os.environ.get(k) for k in
    ["OMP_PROC_BIND", "OMP_PLACES", "GOMP_CPU_AFFINITY"]}}))
Path("actual-affinity.json").write_text(json.dumps(sorted(os.sched_getaffinity(0))))
if mode == "overlap":
    time.sleep(.2)
if mode == "thread-binding":
    import threading
    cpus = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {{cpus[0]}})
    def bound_thread(cpu):
        os.sched_setaffinity(0, {{cpu}})
        # Keep daemon affinity present until this short fake process exits;
        # host contention must not make it disappear before the last sample.
        threading.Event().wait(30)
    for cpu in cpus[1:]:
        threading.Thread(target=bound_thread, args=(cpu,), daemon=True).start()
    time.sleep(.3)
    Path("leader-affinity.json").write_text(json.dumps(sorted(os.sched_getaffinity(0))))
if mode == "timeout":
    import subprocess
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    Path({str(tmp_path / 'child.pid')!r}).write_text(str(child.pid))
    (Path({str(tmp_path)!r}) / f"{{task_id}}.child.pid").write_text(str(child.pid))
    time.sleep(60)
write_products(Path.cwd(), p, int(sys.argv[3]))
if mode == "invalid":
    (Path("output") / f"img_data_{{sys.argv[3]}}.h5").write_bytes(b"invalid")
if mode == "mutate":
    Path(sys.argv[2]).write_bytes(b"changed during run")
if mode == "failure" or (mode == "mixed" and task_id == "baseline-1"):
    sys.exit(7)
'''
    build = tmp_path / "build"
    with provenance.stage_directory(build, "raptor:build", {"fake": True}) as temporary:
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
    manifest = provenance.verify_stage(job, recursive=True)
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


def test_any_ancestor_manifest_is_refused(setup):
    setup["run"] = setup["run"].parent / "output-container/runs/main"
    container = setup["run"].parents[1]
    container.mkdir()
    (container / "manifest.json").write_text(json.dumps({"layout": "container"}))
    with pytest.raises(ValueError, match="existing output stage"):
        invoke(setup)


def test_actual_stage_ancestor_still_refused(setup):
    parent = setup["run"].parent / "completed"
    with provenance.stage_directory(parent, "scientific-stage", {}) as temp:
        (temp / "result.txt").write_text("immutable")
    setup["run"] = parent / "runs/main"
    with pytest.raises(ValueError, match="existing output stage"):
        invoke(setup)


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
    assert str(setup["table"]) in {e["path"] for e in provenance.verify_stage(directory)["inputs"]}


def test_baseline_does_not_depend_on_tables(setup):
    task = deepcopy(setup["task"])
    task.update(table_id="low", table_path=str(setup["table"]), table_spin=.99)
    result = invoke(setup, [task])
    assert result["counts"]["complete"] == 1
    setup["table"].unlink()
    task.update(table_id="high", table_path="/does/not/exist", table_spin=.1)
    result = invoke(setup, [task])
    assert result["counts"]["cached"] == 1
    manifest = provenance.verify_stage(setup["run"] / "jobs/baseline")
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
    result = invoke(setup, tasks, workers=3)
    assert result["counts"]["complete"] == 3
    assert count == 1
    invoke(setup, tasks, workers=2)
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
    result = invoke(setup, tasks, workers=4)
    assert result["counts"]["complete"] == 16
    assert reads == {setup["source"]: 1, setup["table"]: 1}
    result = invoke(setup, tasks, workers=3)
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
    from eht_axion.raptor import jobs
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
    provenance.verify_stage(setup["run"] / "jobs/baseline", recursive=True)
    assert not (setup["run"] / "jobs/second").exists()


def test_build_child_of_run_container_is_allowed(setup):
    setup["run"].mkdir()
    relocated = setup["run"] / "build"
    setup["build"].rename(relocated)
    setup["build"] = relocated
    assert invoke(setup)["counts"]["complete"] == 1


def test_cancel_kills_process_group_and_keeps_attempt(setup, monkeypatch):
    from eht_axion.raptor import jobs
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
    from eht_axion.raptor.plan import plan_generation
    from test_raptor_plan import write_snapshot, write_table

    from eht_axion.experiment import read_experiment

    config = read_experiment()
    config["raptor"]["snapshots"] = {"start": 1600, "stop": 1600, "step": 10}
    source = setup["source"].with_name("data001600.dat")
    write_snapshot(source)
    write_table(setup["table"])
    tables = {"low": setup["table"], "high": setup["table"]}
    plan = plan_generation(config, source.parent, tables)
    plan["expected_core_commit"] = "synthetic-core-commit"
    production_build = setup["build"].with_name("production-build")
    settings = {"settings": config["raptor"]["build"], "source_commit": "synthetic-core-commit"}
    with provenance.stage_directory(production_build, "raptor:build", settings) as work:
        for name in ("RAPTOR", "definitions.h"):
            import shutil
            shutil.copy2(setup["build"] / name, work / name)
        grid = config["raptor"]["grid"]
        text = "".join(f"nxlone{i+1} = {v}\n" for i, v in enumerate(grid["nxlone"]))
        text += "".join(f"xprobminx{i+1} = {v}\n" for i, v in enumerate(grid["xprobmin"]))
        text += "".join(f"xprobmax{i+1} = {v}\n" for i, v in enumerate(grid["xprobmax"]))
        (work / "grid_mks.in").write_text(text + f"hslope = {grid['hslope']}\n")
        (work / "macros.txt").write_text("\n".join(f"#define {k} {v}" for k, v in {
            "POL": "1", "num_frequencies": "1", "num_pixels_1d": "10", "NSPIN": "3",
            "SMR": "1", "AMR": "0", "UNIF": "0", "metric": "MKSBHAC", "DF": "TH",
            "IMGFILE": "1", "SPECFILE": "1", "int_method": "RK4", "SFC": "1", "BPOL": "PLUS",
            "SIGMA_CUT": str(config["raptor"]["build"]["sigma_cut"]),
            "THETAE_MIN": str(config["raptor"]["build"]["thetae_min"]),
            "THETAE_MAX": str(config["raptor"]["build"]["thetae_max"]),
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
            entry.update(size=(build / product).stat().st_size, sha256=provenance.sha256(build / product))
    path.write_text(json.dumps(record))
    provenance.verify_stage(build, recursive=True)
    with pytest.raises(ValueError):
        run_jobs(plan, build, setup["run"], max_jobs=1)
    assert not setup["run"].exists()


def test_production_gate_accepts_bound_plan(production, monkeypatch):
    from eht_axion.raptor import jobs
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
        plan["config"]["raptor"]["build"]["sigma_cut"] = 2.
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


def baseline_tasks(setup, count):
    tasks = [deepcopy(setup["task"]) for _ in range(count)]
    for i, task in enumerate(tasks):
        task["id"] = f"baseline-{i}"
    return tasks


def test_bounded_rolling_overlap_private_states_and_worker_count_resume(setup, monkeypatch):
    from eht_axion.raptor import jobs
    original = jobs._execute
    lock = threading.Lock()
    barrier = threading.Barrier(2)
    active = peak = calls = 0
    contexts = []

    def execute(command, cwd, env, stdout, stderr, timeout):
        nonlocal active, peak, calls
        with lock:
            calls += 1
            call = calls
            active += 1
            peak = max(peak, active)
            contexts.append((provenance._VERIFICATION.get(), provenance._STAGE_TEMPORARIES.get()))
        try:
            if call <= 2:
                barrier.wait(timeout=10)
            return original(command, cwd, env, stdout, stderr, timeout)
        finally:
            with lock:
                active -= 1

    monkeypatch.setattr(jobs, "_execute", execute)
    setup["mode"].write_text("overlap")
    tasks = baseline_tasks(setup, 6)
    result = invoke(setup, tasks, workers=2, threads=1, max_jobs=5)
    assert peak == result["scheduling"]["peak_in_flight"] == 2
    assert calls == 5
    assert result["counts"] == {"complete": 5, "cached": 0, "failed": 0}
    assert [r["id"] for r in result["jobs"]] == [t["id"] for t in tasks[:5]]
    assert len({id(session) for session, _ in contexts}) == 1
    assert all(len(temporaries) == 1 for _, temporaries in contexts)
    assert len({temporaries[0] for _, temporaries in contexts}) == 5
    assert not provenance._STAGE_TEMPORARIES.get()
    assert provenance._VERIFICATION.get() is None
    for entry in result["jobs"]:
        destination = Path(entry["directory"])
        manifest = provenance.verify_stage(destination, recursive=True)
        assert manifest["config"]["id"] == entry["id"]
        assert "workers" not in manifest["config"]
        assert json.loads((destination / "task.json").read_text()) == manifest["config"]
        receipt = json.loads((Path(entry["attempt"]) / "result.json").read_text())
        assert receipt["task_id"] == entry["id"]
        assert receipt["pid"] == receipt["process_group"] == entry["pid"]
        assert receipt["wall_seconds"] >= receipt["process_wall_seconds"] > 0
        assert receipt["actual_cpu_set"] == sorted(os.sched_getaffinity(0))
        assert entry["wall_seconds"] > 0
        assert not any(".building-" in e["path"] for e in manifest["inputs"])
    assert not (setup["run"] / "jobs/baseline-5").exists()
    resumed = invoke(setup, tasks, workers=3, threads=1, max_jobs=5)
    assert resumed["counts"]["cached"] == 5
    assert calls == 5


def test_concurrent_failure_is_independent_and_later_jobs_launch(setup):
    setup["mode"].write_text("mixed")
    result = invoke(setup, baseline_tasks(setup, 5), workers=2, threads=1)
    assert result["counts"] == {"complete": 4, "cached": 0, "failed": 1}
    failed = result["jobs"][1]
    assert failed["id"] == "baseline-1"
    assert failed["returncode"] == 7
    assert not Path(failed["directory"]).exists()
    assert (Path(failed["attempt"]) / "work/model.in").is_file()
    for entry in result["jobs"]:
        if entry["status"] == "complete":
            provenance.verify_stage(entry["directory"], recursive=True)


@pytest.mark.parametrize("options", [
    {"workers": 0}, {"workers": -1}, {"workers": True}, {"workers": 1.5},
    {"threads": 0}, {"threads": True}, {"cpu_affinity": "bogus"},
    {"cpu_list": "0"}, {"workers": 10**6},
])
def test_invalid_allocation_fails_before_launch(setup, options):
    with pytest.raises(ValueError):
        invoke(setup, **options)
    assert not setup["run"].exists()


def test_affinity_uses_parent_policy_and_sanitizes_only_child_env(setup, monkeypatch):
    from eht_axion.raptor import affinity, jobs
    cpus = sorted(os.sched_getaffinity(0))[:2]
    if len(cpus) < 2:
        pytest.skip("needs two allowed CPUs")
    allocations = [[cpus[0]], [cpus[1]]]
    allocate_calls, prefix_calls, commands = [], [], []
    original_prefix = affinity.launch_prefix
    original_execute = jobs._execute

    def allocate(workers, threads, mode="none", cpu_list=None):
        allocate_calls.append((workers, threads, mode, cpu_list))
        return allocations

    def prefix(cpu_set):
        prefix_calls.append(cpu_set)
        return original_prefix(cpu_set)

    def execute(command, cwd, env, stdout, stderr, timeout):
        assert command[0] == original_prefix(allocations[0])[0]
        assert Path(command[0]).is_absolute()
        assert command[1] == "--cpu-list"
        assert command[2] in {str(cpu) for cpu in cpus}
        assert command[3] == str(setup["build"] / "RAPTOR")
        assert env["OMP_PROC_BIND"] == "close"
        assert env["OMP_PLACES"] == "cores"
        assert "GOMP_CPU_AFFINITY" not in env
        assert not any(k.startswith("_PARENT_RAPTOR") for k in env)
        commands.append(command)
        return original_execute(command, cwd, env, stdout, stderr, timeout)

    monkeypatch.setattr(affinity, "allocate_cpu_sets", allocate)
    monkeypatch.setattr(affinity, "launch_prefix", prefix)
    monkeypatch.setattr(jobs, "_execute", execute)
    monkeypatch.setenv("OMP_PROC_BIND", "spread")
    monkeypatch.setenv("OMP_PLACES", "threads")
    monkeypatch.setenv("GOMP_CPU_AFFINITY", "old-invalid-mask")
    setup["mode"].write_text("overlap")
    result = invoke(setup, baseline_tasks(setup, 4), workers=2, threads=1,
                    cpu_affinity="physical", cpu_list=",".join(map(str, cpus)))
    assert result["counts"]["complete"] == 4
    assert allocate_calls == [(2, 1, "physical", ",".join(map(str, cpus)))]
    assert prefix_calls == allocations
    assert result["cpu_sets"] == allocations
    for entry in result["jobs"]:
        mask = allocations[entry["worker_slot"]]
        destination = Path(entry["directory"])
        assert entry["cpu_set"] == entry["actual_cpu_set"] == mask
        assert json.loads((destination / "actual-affinity.json").read_text()) == mask
        assert json.loads((destination / "binding-environment.json").read_text()) == {
            "OMP_PROC_BIND": "close", "OMP_PLACES": "cores", "GOMP_CPU_AFFINITY": None}
        receipt = json.loads((Path(entry["attempt"]) / "result.json").read_text())
        assert receipt["cpu_set"] == receipt["actual_cpu_set"] == mask
        assert receipt["omp_proc_bind"] == "close"
        assert receipt["omp_places"] == "cores"
        assert provenance.verify_stage(destination)["config"]["cpu_set"] == mask
    assert os.environ["GOMP_CPU_AFFINITY"] == "old-invalid-mask"
    # Cache hits complete much faster than original subprocesses. Slot identity
    # must stay deterministic instead of following completion-order refill.
    setup["mode"].write_text("failure")
    again = invoke(setup, baseline_tasks(setup, 4), workers=2, threads=1,
                   cpu_affinity="physical", cpu_list=",".join(map(str, cpus)))
    assert again["counts"] == {"complete": 0, "cached": 4, "failed": 0}
    assert [e["cpu_set"] for e in again["jobs"]] == [allocations[i % 2] for i in range(4)]
    assert os.environ["OMP_PROC_BIND"] == "spread"
    assert os.environ["OMP_PLACES"] == "threads"


@pytest.mark.parametrize("allocation", [[[0], [0]], [[0, 1], [2]], [None, [1]], [[0]]])
def test_malformed_policy_allocation_rejected(setup, monkeypatch, allocation):
    from eht_axion.raptor import affinity
    monkeypatch.setattr(affinity, "allocate_cpu_sets", lambda *a, **kw: allocation)
    with pytest.raises(ValueError, match="CPU|affinity"):
        invoke(setup, workers=2, threads=1, cpu_affinity="logical")
    assert not setup["run"].exists()


def _assert_process_stopped(pid):
    # A killed grandchild may remain a zombie until its external reaper runs.
    for _ in range(100):
        try:
            state = Path(f"/proc/{pid}/stat").read_text().split(") ", 1)[1].split()[0]
        except FileNotFoundError:
            return
        if state == "Z":
            return
        time.sleep(.02)
    pytest.fail(f"process {pid} is still running")


def test_parallel_timeout_kills_descendants_but_continues_selected_prefix(setup):
    setup["mode"].write_text("timeout")
    result = invoke(setup, baseline_tasks(setup, 3), workers=2, threads=1, timeout_seconds=1.5)
    assert result["counts"] == {"complete": 0, "cached": 0, "failed": 3}
    for entry in result["jobs"]:
        receipt = json.loads((Path(entry["attempt"]) / "result.json").read_text())
        assert receipt["timed_out"] is True
        assert receipt["returncode"] == -9
        _assert_process_stopped(receipt["pid"])
        child = setup["mode"].parent / f"{entry['id']}.child.pid"
        assert child.is_file()
        _assert_process_stopped(int(child.read_text()))


@pytest.mark.parametrize("exception_type", [KeyboardInterrupt, SystemExit])
def test_main_cancellation_kills_all_groups_joins_and_stops_launches(setup, monkeypatch, exception_type):
    from eht_axion.raptor import jobs
    setup["mode"].write_text("timeout")
    real_wait = jobs.wait
    called = False

    def interrupt(fs, **kwargs):
        nonlocal called
        if called:
            return real_wait(fs, **kwargs)
        called = True
        deadline = time.monotonic() + 10
        while len(list(setup["mode"].parent.glob("baseline-*.child.pid"))) < 2:
            if time.monotonic() > deadline:
                pytest.fail("fake subprocesses did not start")
            time.sleep(.02)
        raise exception_type("test cancellation")

    monkeypatch.setattr(jobs, "wait", interrupt)
    with pytest.raises(exception_type, match="test cancellation"):
        invoke(setup, baseline_tasks(setup, 5), workers=2, threads=1)
    receipts = list((setup["run"] / "attempts").glob("*/*/result.json"))
    assert len(receipts) == 2
    for path in receipts:
        receipt = json.loads(path.read_text())
        assert receipt["cancelled"] is True
        assert receipt["complete"] is False
        assert receipt["returncode"] == -9
        _assert_process_stopped(receipt["pid"])
        child = setup["mode"].parent / f"{receipt['task_id']}.child.pid"
        _assert_process_stopped(int(child.read_text()))
        assert (path.parent / "work/model.in").is_file()
    assert not list((setup["run"] / "jobs").glob("*.building-*"))
    assert not list((setup["run"] / "jobs").glob("*.lock"))
    assert not any(t.name.startswith("raptor-job") for t in threading.enumerate())
    assert provenance._VERIFICATION.get() is None


def test_real_main_thread_sigint_cancels_eight_fake_groups_and_joins(setup):
    if len(os.sched_getaffinity(0)) < 8:
        pytest.skip("eight workers require eight allowed CPUs")
    setup["mode"].write_text("timeout")
    payload = json.dumps({"tasks": baseline_tasks(setup, 10), "build": str(setup["build"]),
                          "run": str(setup["run"]), "marker": str(setup["mode"].parent / "joined.json")})
    harness = '''
import json, sys, threading
from pathlib import Path
from eht_axion import provenance
from eht_axion.raptor.jobs import run_jobs
p = json.loads(sys.argv[1])
try:
    run_jobs({"tasks": p["tasks"]}, Path(p["build"]), Path(p["run"]), workers=8, threads=1)
except KeyboardInterrupt:
    workers = [t.name for t in threading.enumerate() if t.name.startswith("raptor-job")]
    Path(p["marker"]).write_text(json.dumps({"workers": workers,
        "session_closed": provenance._VERIFICATION.get() is None,
        "temporaries_closed": not provenance._STAGE_TEMPORARIES.get()}))
    sys.exit(130)
else:
    sys.exit("expected SIGINT cancellation")
'''
    process = subprocess.Popen([sys.executable, "-c", harness, payload],
                               cwd=Path(__file__).resolve().parents[1],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, start_new_session=True)
    try:
        deadline = time.monotonic() + 20
        while len(list(setup["mode"].parent.glob("baseline-*.child.pid"))) < 8:
            if process.poll() is not None:
                out, err = process.communicate()
                pytest.fail(f"runner exited before readiness: {out}\n{err}")
            if time.monotonic() > deadline:
                pytest.fail("eight fake groups did not start")
            time.sleep(.02)
        os.kill(process.pid, signal.SIGINT)
        out, err = process.communicate(timeout=15)
        assert process.returncode == 130, (out, err)
        marker = json.loads((setup["mode"].parent / "joined.json").read_text())
        assert marker == {"workers": [], "session_closed": True, "temporaries_closed": True}
        receipts = list((setup["run"] / "attempts").glob("*/*/result.json"))
        assert len(receipts) == 8
        assert len(list((setup["run"] / "attempts").iterdir())) == 8
        for path in receipts:
            receipt = json.loads(path.read_text())
            assert receipt["cancelled"] is True
            assert receipt["complete"] is False
            assert receipt["returncode"] == -9
            assert receipt["pid"] == receipt["process_group"]
            _assert_process_stopped(receipt["pid"])
            child = setup["mode"].parent / f"{receipt['task_id']}.child.pid"
            _assert_process_stopped(int(child.read_text()))
            assert (path.parent / "work/model.in").is_file()
        assert not list((setup["run"] / "jobs").iterdir())
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
        process.communicate(timeout=10)
        # Keep failures in this test from leaking fake processes, whose groups
        # intentionally differ from the runner's session.
        for child in setup["mode"].parent.glob("baseline-*.child.pid"):
            try:
                os.killpg(os.getpgid(int(child.read_text())), signal.SIGKILL)
            except ProcessLookupError:
                pass


def test_cpu_observation_unions_thread_masks_not_only_openmp_leader(tmp_path, monkeypatch):
    from eht_axion.raptor import jobs
    task_root = tmp_path / "123/task"
    task_root.mkdir(parents=True)
    for tid in (123, 124, 125):
        (task_root / str(tid)).mkdir()
    original_glob = Path.glob

    def glob(path, pattern):
        if path == Path("/proc/123/task"):
            return original_glob(task_root, pattern)
        return original_glob(path, pattern)

    def mask(tid):
        if tid == 125:
            raise ProcessLookupError()
        return {0} if tid == 123 else {2, 4}

    monkeypatch.setattr(Path, "glob", glob)
    monkeypatch.setattr(jobs.os, "sched_getaffinity", mask)
    assert jobs._observed_cpu_set(123) == [0, 2, 4]


def test_receipt_observes_all_thread_masks_with_single_core_leader(setup):
    cpus = sorted(os.sched_getaffinity(0))[:2]
    if len(cpus) < 2:
        pytest.skip("needs two allowed CPUs")
    setup["mode"].write_text("thread-binding")
    result = invoke(setup, workers=1, threads=2, cpu_affinity="logical",
                    cpu_list=",".join(map(str, cpus)))
    assert result["counts"]["complete"] == 1
    entry = result["jobs"][0]
    assert entry["actual_cpu_set"] == cpus
    assert json.loads((Path(entry["directory"]) / "leader-affinity.json").read_text()) == cpus[:1]
    receipt = json.loads((Path(entry["attempt"]) / "result.json").read_text())
    assert receipt["actual_cpu_set"] == receipt["cpu_set"] == cpus


def test_production_memory_guard_uses_actual_selected_inflight(production, monkeypatch):
    from eht_axion.raptor import affinity
    setup, plan = production
    calls = []

    def reject(workers):
        calls.append(workers)
        raise ValueError("synthetic memory headroom failure")

    monkeypatch.setattr(affinity, "check_memory", reject)
    with pytest.raises(ValueError, match="headroom"):
        run_jobs(plan, setup["build"], setup["run"], workers=4, threads=1, max_jobs=2)
    assert calls == [2]
    assert not setup["run"].exists()
    # An empty authorized prefix has no in-flight memory demand.
    result = run_jobs(plan, setup["build"], setup["run"], workers=4, threads=1, max_jobs=0)
    assert result["selected_tasks"] == 0
    assert calls == [2]


def test_tasks_only_skips_production_memory_guard(setup, monkeypatch):
    from eht_axion.raptor import affinity

    def reject(*args, **kwargs):
        pytest.fail("synthetic tasks-only API must not call memory guard")

    monkeypatch.setattr(affinity, "check_memory", reject)
    assert invoke(setup)["counts"]["complete"] == 1


def test_none_clears_unsigned_inherited_binding_and_keeps_binary_command(setup, monkeypatch):
    from eht_axion.raptor import jobs
    original = jobs._execute
    inherited = {"OMP_PROC_BIND": "spread", "OMP_PLACES": "threads", "GOMP_CPU_AFFINITY": "0"}
    for key, value in inherited.items():
        monkeypatch.setenv(key, value)

    def execute(command, cwd, env, stdout, stderr, timeout):
        assert command[0] == str(setup["build"] / "RAPTOR")
        assert all(key not in env for key in inherited)
        return original(command, cwd, env, stdout, stderr, timeout)

    monkeypatch.setattr(jobs, "_execute", execute)
    result = invoke(setup, workers=1, threads=1)
    assert result["counts"]["complete"] == 1
    job = setup["run"] / "jobs/baseline"
    assert json.loads((job / "binding-environment.json").read_text()) == dict.fromkeys(inherited)
    for key, value in inherited.items():
        assert os.environ[key] == value
    monkeypatch.setenv("OMP_PROC_BIND", "true")
    assert invoke(setup, threads=1)["counts"]["cached"] == 1


def test_per_job_source_checks_do_not_rescan_all_selected_sources(setup, monkeypatch):
    from eht_axion.raptor import jobs
    tasks = baseline_tasks(setup, 5)
    sources = []
    for i, task in enumerate(tasks):
        path = setup["source"].with_name(f"fluid-{i}.dat")
        path.write_bytes(b"synthetic")
        sources.append(path)
        task["source_path"] = str(path)
    original = jobs._assert_sources_unchanged
    scans = []

    def check(paths):
        scans.append(set(paths).intersection(sources))
        original(paths)

    monkeypatch.setattr(jobs, "_assert_sources_unchanged", check)
    result = invoke(setup, tasks, workers=3, threads=1)
    assert result["counts"]["complete"] == 5
    assert scans[0] == set(sources)  # one main-thread root-guard closure
    assert len(scans) == 1 + 2 * len(tasks)
    assert all(len(scan) == 1 for scan in scans[1:])

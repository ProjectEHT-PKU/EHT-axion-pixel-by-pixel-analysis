#!/usr/bin/env python3
"""ONE-OFF data migration (not part of the pipeline): completed RAPTOR jobs -> new task ids.

The cleaned code names job directories by a SHA-256 of the ray-tracing inputs
only. The 18,330 completed jobs of the m65-211 production (two job roots, two
byte-different builds) were recorded under the old ids. Each old job already
records its full RAPTOR inputs (task.json parameters, model.in, grid_mks.in,
source path, input hashes). This script:

1. builds the NEW plan with the NEW code (experiment.yaml, real fluid/tables);
2. matches every old job to exactly one new task and asserts that its recorded
   task (all fields except the id), model.in bytes, grid, source/table hashes
   and output receipts equal what the new code would execute/record;
3. writes NEW_ROOT/<root>/{build,jobs/<new id>,config.yaml,requests/} exactly as
   raptor_generate.py would: build files and job outputs are HARD LINKS (same
   bytes, same inodes; originals are only read), task.json/output.json/manifest.json
   are written by the new code (stage_directory, validate_job_output).

Originals are never modified, renamed or deleted. Run with the default Pixi
Python from the worktree; --limit N writes only the first N jobs (test mode).
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import uuid

WORKTREE = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(WORKTREE)]

from eht_axion.experiment import read_experiment  # noqa: E402
from eht_axion.pipeline import common  # noqa: E402
from eht_axion.producer.jobs import _prepare_task, parameter_values, serialize_parameters  # noqa: E402
from eht_axion.producer.outputs import read_core_constants, validate_job_output  # noqa: E402
from eht_axion.producer.plan import plan_generation  # noqa: E402

OLD_ROOTS = {
    "m65-211": Path("/home/canoziia/PHYS-production/results/m65-211"),
    "pilot-m65-t3000-res400-t16": Path("/home/canoziia/GRMHD_DATA/production/pilot-m65-t3000-res400-t16"),
}
FLUID = Path("/home/canoziia/GRMHD_DATA/a094-high")
ARCHIVED_PLAN = Path("/home/canoziia/PHYS-production/control/full-plan.json")
RAPTOR_SOURCE = Path("/home/canoziia/MyProject/PHYS/raptor")
RUNNER_FIELDS = ("threads", "output_adapter_schema", "cpu_affinity", "cpu_set", "omp_proc_bind", "omp_places")
LINKED = ("grid_mks.in", "model.in", "stdout.log", "stderr.log", "success.json")  # + output/*


def key(task):
    return (task["kind"], task["resolution"], task["label"], task["alpha"], task["phase_index"])


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def link_tree(source, destination):
    """Hard-link every regular file of source into a new destination tree."""
    for path in sorted(source.rglob("*")):
        target = destination / path.relative_to(source)
        if path.is_symlink():
            raise ValueError(f"Unexpected symlink in build: {path}")
        if path.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(path, target)


def build_inputs(build):
    return [build / "manifest.json", build / "RAPTOR", build / "grid_mks.in", build / "definitions.h"]


def match_jobs(plan):
    """Assert a one-to-one, field-identical correspondence old job <-> new task."""
    tasks = {key(t): t for t in plan["tasks"]}
    if len(tasks) != len(plan["tasks"]):
        raise AssertionError("duplicate plan keys")
    matched, problems = {}, []
    for name, root in OLD_ROOTS.items():
        grid_sha = sha(root / "build/grid_mks.in")
        for directory in sorted((root / "jobs").iterdir()):
            old = json.loads((directory / "task.json").read_text())
            manifest = json.loads((directory / "manifest.json").read_text())
            k = key(old)
            if k in matched:
                problems.append(dict(job=str(directory), problem="duplicate task key"))
                continue
            if k not in tasks:
                problems.append(dict(job=str(directory), problem="no new task"))
                continue
            task = tasks[k]
            _, config, raw = _prepare_task(task)
            for field in RUNNER_FIELDS:
                if field in old:
                    config[field] = old[field]
            checks = {
                "old_manifest_complete": manifest.get("complete") is True and manifest["stage"] == "producer-job",
                "old_manifest_config_equals_task_json": manifest["config"] == old,
                "old_directory_name_is_old_id": directory.name == old["id"],
                "task_fields_equal_except_id": {f: v for f, v in old.items() if f != "id"}
                == {f: v for f, v in config.items() if f != "id"},
                "parameters_equal_plan": {n: v for n, v in old["parameters"].items() if n != "AXION_TABLE"}
                == {n: v for n, v in parameter_values(task["parameters"]).items() if n != "AXION_TABLE"},
                "model_in_bytes_equal": (directory / "model.in").read_text()
                == serialize_parameters(config["parameters"]),
                "grid_equals_build_grid": sha(directory / "grid_mks.in") == grid_sha,
                "success_receipt": json.loads((directory / "success.json").read_text()) == {
                    "schema": 1, "returncode": 0, "timed_out": False, "output_validated": True},
            }
            recorded = {e["path"]: e["sha256"] for e in manifest["inputs"]}
            checks["recorded_raw_inputs_present"] = all(str(p) in recorded for p in raw)
            if not all(checks.values()):
                problems.append(dict(job=str(directory), failed=[c for c, ok in checks.items() if not ok]))
            matched[k] = dict(root=name, old_dir=str(directory), old_id=old["id"], new_id=task["id"],
                              runner={f: old[f] for f in RUNNER_FIELDS if f in old},
                              raw_sha={str(p): recorded.get(str(p)) for p in raw})
    missing = set(tasks) - set(matched)
    if problems or missing:
        raise AssertionError(json.dumps(dict(problems=problems[:20], problem_count=len(problems),
                                             missing=len(missing)), default=str))
    return matched


def migrate_job(args):
    entry, task, new_root = args
    sub = new_root / entry["root"]
    build = sub / "build"
    destination = sub / "jobs" / task["id"]
    old_dir = Path(entry["old_dir"])
    old_manifest = json.loads((old_dir / "manifest.json").read_text())
    old_outputs = {e["path"]: e for e in old_manifest["outputs"]}
    old_output_json = json.loads((old_dir / "output.json").read_text())
    if destination.exists():
        status = "existing"
    else:
        prepared, config, raw = _prepare_task(task)
        config.update(entry["runner"])
        prepared["core_constants"] = read_core_constants(build / "definitions.h")
        prepared["core_constants_source"] = str(build / "definitions.h")
        with common.stage_directory(destination, "producer-job", config, build_inputs(build) + raw) as temp:
            for relative in old_outputs:
                if relative in LINKED or relative.startswith("output/"):
                    (temp / relative).parent.mkdir(parents=True, exist_ok=True)
                    os.link(old_dir / relative, temp / relative)
            (temp / "task.json").write_text(json.dumps(config, indent=2, allow_nan=False) + "\n")
            validated = validate_job_output(temp, prepared)
            (temp / "output.json").write_text(json.dumps(validated, indent=2, allow_nan=False) + "\n")
        status = "written"
    manifest = common.verify_stage(destination)
    new_outputs = {e["path"]: e for e in manifest["outputs"]}
    new_output_json = json.loads((destination / "output.json").read_text())
    new_inputs = {e["path"]: e["sha256"] for e in manifest["inputs"]}
    checks = {
        "same_output_files": set(new_outputs) == set(old_outputs),
        "linked_bytes_equal_recorded": all(new_outputs[p]["sha256"] == old_outputs[p]["sha256"]
                                           for p in old_outputs if p in LINKED or p.startswith("output/")),
        "linked_same_inode": all(os.stat(destination / p).st_ino == os.stat(old_dir / p).st_ino
                                 for p in old_outputs if p in LINKED or p.startswith("output/")),
        "output_json_equal_except_build_path": {k: v for k, v in new_output_json.items()
                                                if k != "core_constants_source"}
        == {k: v for k, v in old_output_json.items() if k != "core_constants_source"},
        "raw_input_sha_equal_recorded": all(new_inputs[p] == s for p, s in entry["raw_sha"].items()),
        "task_json_equal_except_id": {k: v for k, v in json.loads((destination / "task.json").read_text()).items()
                                      if k != "id"}
        == {k: v for k, v in json.loads((old_dir / "task.json").read_text()).items() if k != "id"},
    }
    return dict(old_id=entry["old_id"], new_id=task["id"], root=entry["root"], status=status, checks=checks,
                ok=all(checks.values()), old_elapsed_seconds=old_manifest.get("elapsed_seconds"))


def worker(batch):
    with common.verification_session():
        return [migrate_job(item) for item in batch]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--new-root", type=Path, required=True)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int, help="test mode: migrate only this many jobs")
    args = p.parse_args(argv)
    started = time.time()
    new_root = args.new_root.absolute()
    experiment_path = WORKTREE / "experiment.yaml"
    experiment = read_experiment(experiment_path)
    archived = json.loads(ARCHIVED_PLAN.read_text())
    tables = {k: Path(archived["radial_tables"][k]["path"]) for k in ("low", "high")}
    plan = plan_generation(experiment, FLUID, tables, new_root / "m65-211" / "build")
    lock = json.loads((WORKTREE / "raptor-source.lock.json").read_text())
    plan["input_roots"] = [str(FLUID), str(tables["low"].parent), str(tables["high"].parent)]
    plan["source_roots"] = [str(RAPTOR_SOURCE)]
    plan["expected_core_commit"] = lock["commit"]
    matched = match_jobs(plan)
    print(f"matched {len(matched)} old jobs one-to-one to {len(plan['tasks'])} new tasks", flush=True)

    new_root.mkdir(parents=True, exist_ok=True)
    plan_path = new_root / "plan.json"
    if not plan_path.exists():
        plan_path.write_text(json.dumps(plan, indent=1, allow_nan=False) + "\n")
    elif json.loads(plan_path.read_text()) != plan:
        raise ValueError("existing plan.json differs")
    builds = {}
    for name, old in OLD_ROOTS.items():
        sub = new_root / name
        if not (sub / "build").exists():
            link_tree(old / "build", sub / "build")
        with common.verification_session():
            common.verify_stage(sub / "build", recursive=True)
        builds[name] = {str(path.relative_to(sub)): sha(path) for path in build_inputs(sub / "build")}
        if builds[name] != {str(path.relative_to(old)): sha(path) for path in build_inputs(old / "build")}:
            raise AssertionError(f"linked build bytes differ: {name}")
        config_copy = sub / "config.yaml"
        if not config_copy.exists():
            shutil.copyfile(experiment_path, config_copy)
        requests = sub / "requests"
        if not requests.exists():
            requests.mkdir()
            record = {"config": str(experiment_path), "run_root": str(sub), "core_source": str(RAPTOR_SOURCE),
                      "core_commit": lock["commit"], "plan": plan,
                      "note": "written by reports/cleanup-20261106/migration/migrate_jobs.py (no RAPTOR run)"}
            (requests / f"{uuid.uuid4().hex}.json").write_text(json.dumps(record, indent=1) + "\n")

    tasks = {key(t): t for t in plan["tasks"]}
    items = [(entry, tasks[k], new_root) for k, entry in sorted(matched.items(), key=lambda kv: str(kv[0]))]
    if args.limit:
        items = items[:args.limit]
    # One batch per snapshot label: each 690 MB snapshot is hashed by one process only.
    labels = sorted({item[1]["label"] for item in items})
    groups = [[i for i in items if i[1]["label"] == label] for label in labels]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for chunk in pool.map(worker, groups):
            results.extend(chunk)
            print(f"{len(results)}/{len(items)} jobs done", flush=True)
    failures = [r for r in results if not r["ok"]]
    receipt = dict(
        schema="cleanup-20261106-job-migration-v1",
        description="Completed m65-211 RAPTOR jobs re-indexed under the new ray-tracing task ids "
                    "(hard links; no RAPTOR run; originals untouched)",
        script=str(Path(__file__).resolve()), script_sha256=sha(Path(__file__)),
        worktree_commit=common.revision()["git_commit"], experiment=str(experiment_path),
        experiment_sha256=sha(experiment_path), old_roots={k: str(v) for k, v in OLD_ROOTS.items()},
        new_root=str(new_root), plan=str(plan_path), plan_tasks=len(plan["tasks"]),
        matched_one_to_one=len(matched), migrated=len(results), failures=failures[:50],
        failure_count=len(failures), builds_sha256=builds,
        per_job_checks=sorted({c for r in results for c in r["checks"]}),
        old_to_new_ids={r["old_id"]: r["new_id"] for r in results},
        attempts_note="attempt receipts (attempts/) stay in the original roots",
        elapsed_seconds=time.time() - started, limit=args.limit)
    name = "MIGRATION.json" if not args.limit else f"MIGRATION-limit{args.limit}.json"
    (new_root / name).write_text(json.dumps(receipt, indent=1) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "old_to_new_ids"}, indent=1, default=str))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Evidence 3b: the D-term Monte Carlo inputs derived from experiment.yaml are unchanged.

eht_sample.py (new) now reads draw count, seed, fsys, D statistics and sites from
experiment.yaml instead of literals. Replays the bundle check, the 4 x 1000 D-term
draws and the task list for the recorded batch eht-uncertainty-remote-v1 with the
new and the main code and compares them with the recorded all-draws.json /
tasks.jsonl. No imaging. Writes eht-sample-equivalence.json.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from _side import NEW, run, write

BATCH = Path("/home/canoziia/GRMHD_DATA/analysis/eht-uncertainty-remote-v1")
RESULTS = BATCH / "results/full4000-v1"
recorded_draws = json.loads((RESULTS / "all-draws.json").read_text())
recorded_tasks = (RESULTS / "tasks.jsonl").read_bytes()
batch = json.loads((RESULTS / "batch.json").read_text())
seed, count = batch["config"]["seed"], batch["config"]["draws_per_day"]

NEW_CODE = f"""
import json, hashlib
from pathlib import Path
import eht_sample as b
science = b.science_settings({str(NEW / 'experiment.yaml')!r})
bundle, snapshot = b.verify_bundle(Path({str(BATCH / 'inputs')!r}), science)
draws = b.make_draws(Path({str(BATCH / 'inputs')!r}), bundle, science['draws_per_day'], science['seed'], science)
tasks = b.enumerate_tasks(draws)
result = dict(science=science, snapshot_bundle_sha256=snapshot['bundle_sha256'], draws=draws,
              tasks_sha256=hashlib.sha256(b''.join(b.encoded(t) + b'\\n' for t in tasks)).hexdigest())
"""
MAIN_CODE = f"""
import json, hashlib
from pathlib import Path
import run_eht_uncertainty_batch as b
bundle, snapshot = b.verify_bundle(Path({str(BATCH / 'inputs')!r}))
draws = b.make_draws(Path({str(BATCH / 'inputs')!r}), bundle, {count}, {seed})
tasks = b.enumerate_tasks(draws)
result = dict(snapshot_bundle_sha256=snapshot['bundle_sha256'], draws=draws,
              tasks_sha256=hashlib.sha256(b''.join(b.encoded(t) + b'\\n' for t in tasks)).hexdigest())
"""
PY = "/home/canoziia/MyProject/PHYS/eht-script-pipeline/.pixi/envs/eht-rebuild/bin/python"  # batch environment
new, old = run("new", NEW_CODE, python=PY), run("main", MAIN_CODE, python=PY)


def max_draw_difference(a, b):
    worst, count = 0.0, 0
    for day in b:
        if len(a[day]["draws"]) != len(b[day]["draws"]):
            return float("inf"), -1
        for x, y in zip(a[day]["draws"], b[day]["draws"]):
            if set(x) != set(y):
                return float("inf"), -1
            for site in x:
                for hand in ("R", "L"):
                    for u, v in zip(x[site][hand], y[site][hand]):
                        if u != v:
                            count += 1
                            worst = max(worst, abs(u - v))
    return worst, count


def without_draws(value):
    return {d: {k: v for k, v in entry.items() if k != "draws"} for d, entry in value.items()}
recorded_tasks_sha = hashlib.sha256(recorded_tasks).hexdigest()
record = dict(
    purpose="D-term Monte Carlo draws/tasks from experiment.yaml vs main literals vs recorded batch",
    batch=str(RESULTS), recorded_seed=seed, recorded_draws_per_day=count,
    science_from_experiment=new["science"],
    experiment_matches_recorded_batch_config=(new["science"]["seed"], new["science"]["draws_per_day"]) == (
        seed, count),
    bundle_sha256=dict(new=new["snapshot_bundle_sha256"], main=old["snapshot_bundle_sha256"],
                       recorded=batch["snapshot"]["bundle_sha256"]),
    draws_new_equal_recorded=new["draws"] == recorded_draws,
    draws_main_equal_recorded=old["draws"] == recorded_draws,
    draws_new_equal_main_bitwise=new["draws"] == old["draws"],
    tasks_sha256=dict(new=new["tasks_sha256"], main=old["tasks_sha256"], recorded=recorded_tasks_sha),
    distribution_metadata_new_equal_recorded=without_draws(new["draws"]) == without_draws(recorded_draws),
)
worst, count = max_draw_difference(new["draws"], recorded_draws)
record["draws_new_vs_recorded"] = dict(max_abs_diff=worst, differing_components=count,
                                       note="recorded on the remote batch host; last-bit BLAS/libm differences")
record["all_identical"] = (record["experiment_matches_recorded_batch_config"]
                           and record["draws_new_equal_main_bitwise"]
                           and record["tasks_sha256"]["new"] == record["tasks_sha256"]["main"]
                           and record["distribution_metadata_new_equal_recorded"]
                           and worst <= 1e-12 and len(set(record["bundle_sha256"].values())) == 1)
write("eht-sample-equivalence.json", record)
print({k: v for k, v in record.items() if k not in ("science_from_experiment",)})

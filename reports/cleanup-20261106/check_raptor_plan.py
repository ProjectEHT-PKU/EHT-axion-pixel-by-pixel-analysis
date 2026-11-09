#!/usr/bin/env python3
"""Evidence 2: RAPTOR task plan of the new code == main's plan == archived frozen plan.

Real inputs: BHAC fluid /home/canoziia/GRMHD_DATA/a094-high and the two radial
tables recorded in /home/canoziia/PHYS-production/control/full-plan.json; plans
are built with hash_inputs=False (footers and tables only). Only task ids (new
ray-tracing digest) and table hashes (now always recorded) may differ.
Writes raptor-plan-equivalence.json.
"""
from __future__ import annotations

import json
from pathlib import Path
import tempfile

from _side import MAIN_COMMIT, NEW, run, write

FLUID = "/home/canoziia/GRMHD_DATA/a094-high"
ARCHIVED = Path("/home/canoziia/PHYS-production/control/full-plan.json")
archived = json.loads(ARCHIVED.read_text())
TABLES = {key: archived["radial_tables"][key]["path"] for key in ("low", "high")}
FIELDS = ("label", "snapshot_name", "time_tg", "source_path", "kind", "resolution", "alpha", "omega", "norm",
          "phase_index", "output_slot", "table_id", "table_path", "table_spin", "parameters")
tmp = Path(tempfile.mkdtemp(prefix="cleanup-plan-"))

main_plan = tmp / "main-plan.json"
run("main", f"""
import json, subprocess, tempfile
from pathlib import Path
from eht_axion.producer import plan_generation, read_production_config
text = subprocess.run(['git', '-C', {str(NEW)!r}, 'show', '{MAIN_COMMIT}:configs/production.yaml'],
                      check=True, capture_output=True, text=True).stdout
path = Path(tempfile.mkdtemp()) / 'production.yaml'
path.write_text(text)
plan = plan_generation(read_production_config(path), {FLUID!r}, {{k: Path(v) for k, v in {TABLES!r}.items()}})
Path({str(main_plan)!r}).write_text(json.dumps(plan))
result = len(plan['tasks'])
""")
new_plan = tmp / "new-plan.json"
run("new", f"""
import json
from pathlib import Path
from eht_axion.experiment import read_experiment
from eht_axion.producer.plan import plan_generation
plan = plan_generation(read_experiment(), {FLUID!r}, {{k: Path(v) for k, v in {TABLES!r}.items()}})
Path({str(new_plan)!r}).write_text(json.dumps(plan))
result = len(plan['tasks'])
""")
old, new = json.loads(main_plan.read_text()), json.loads(new_plan.read_text())


def key(task):
    return (task["kind"], task["resolution"], task["label"], task["alpha"], task["phase_index"])


def compare_tasks(reference, candidate):
    ref = {key(t): t for t in reference["tasks"]}
    cand = {key(t): t for t in candidate["tasks"]}
    mismatches = []
    for k in sorted(set(ref) | set(cand), key=str):
        if k not in ref or k not in cand:
            mismatches.append(dict(task=list(map(str, k)), problem="missing"))
            continue
        for field in FIELDS:
            if ref[k][field] != cand[k][field]:
                mismatches.append(dict(task=list(map(str, k)), field=field, old=ref[k][field],
                                       new=cand[k][field]))
    return dict(reference_tasks=len(ref), candidate_tasks=len(cand), unique_keys=len(ref) == len(
        reference["tasks"]) and len(cand) == len(candidate["tasks"]), mismatches=mismatches[:50],
        mismatch_count=len(mismatches), identical=not mismatches and len(ref) == len(cand))


def strip_tables(tables):
    return {k: {f: v for f, v in t.items() if f != "sha256"} for k, t in tables.items()}


requirements_old = {k: v for k, v in old["requirements"].items() if k != "official_setup_args"}
requirements_new = {k: v for k, v in new["requirements"].items() if k != "official_setup_args"}
record = dict(
    purpose="RAPTOR plan: new code vs main code vs archived frozen plan (per-task ray-tracing inputs)",
    inputs=dict(fluid_root=FLUID, radial_tables=TABLES, hash_inputs=False, archived_plan=str(ARCHIVED)),
    counts=dict(main=old["counts"], new=new["counts"], archived=archived["counts"],
                identical=old["counts"] == new["counts"] == archived["counts"]),
    new_vs_main_tasks=compare_tasks(old, new),
    new_vs_archived_tasks=compare_tasks(archived, new),
    task_order_identical=[key(t) for t in old["tasks"]] == [key(t) for t in new["tasks"]],
    snapshots_identical=old["snapshots"] == new["snapshots"],
    radial_tables_identical_except_hash=strip_tables(old["radial_tables"]) == strip_tables(new["radial_tables"]),
    new_table_sha256={k: t["sha256"] for k, t in new["radial_tables"].items()},
    archived_table_sha256={k: t["sha256"] for k, t in archived["radial_tables"].items()},
    build_settings_identical=old["config"]["producer"]["build"] == new["config"]["raptor"]["build"],
    grid_identical=old["config"]["producer"]["grid"] == new["config"]["raptor"]["grid"],
    setup_arguments=dict(main=old["requirements"]["official_setup_args"],
                         new=new["requirements"]["official_setup_args"],
                         same_set=sorted(old["requirements"]["official_setup_args"])
                         == sorted(new["requirements"]["official_setup_args"]),
                         note="new plan lists build._arguments order = the command the build really runs"),
    other_requirements_identical=requirements_old == requirements_new,
    ids=dict(new_unique=len({t["id"] for t in new["tasks"]}) == len(new["tasks"]),
             unchanged_ids=len({t["id"] for t in old["tasks"]} & {t["id"] for t in new["tasks"]}),
             note="all ids change by design: new digest of ray-tracing inputs only"),
    plan_config=dict(main_keys=sorted(old["config"]), new_keys=sorted(new["config"]),
                     main_config_sha256=old["config_sha256"], new_config_sha256=new["config_sha256"]),
)
record["all_ray_tracing_inputs_identical"] = all((
    record["counts"]["identical"], record["new_vs_main_tasks"]["identical"],
    record["new_vs_archived_tasks"]["identical"], record["task_order_identical"], record["snapshots_identical"],
    record["radial_tables_identical_except_hash"], record["build_settings_identical"], record["grid_identical"],
    record["setup_arguments"]["same_set"], record["other_requirements_identical"], record["ids"]["new_unique"]))
write("raptor-plan-equivalence.json", record)
print(json.dumps({k: v for k, v in record.items() if k in ("counts", "all_ray_tracing_inputs_identical")},
                 indent=1))
print("main mismatches", record["new_vs_main_tasks"]["mismatch_count"], "archived mismatches",
      record["new_vs_archived_tasks"]["mismatch_count"])

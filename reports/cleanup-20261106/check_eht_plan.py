#!/usr/bin/env python3
"""Evidence 3: what the eht-imaging driver derives from configuration (dry-run --plan).

For each recorded Paper VII stage (eht-paper-{I,PD}-*-v1, MC pilot) the same
CLI is replayed in --plan mode with main's reconstruct_eht_official.py (configs/
eht-official.yaml) and with the new eht_reconstruct.py (experiment.yaml + the
four locations as CLI options). No imaging is executed. Compared: every plan
field except `cli` and the code/config entries of `inputs`; data, upstream-stage
and ehtim inputs must match by path and sha256. Writes eht-plan-equivalence.json.
Run with the eht-rebuild environment (ehtim 1.3.2 metadata is enumerated).
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile

from _side import MAIN, MAIN_COMMIT, NEW, diff, run, write

ANALYSIS = Path("/home/canoziia/GRMHD_DATA/analysis")
PY = "/home/canoziia/MyProject/PHYS/eht-script-pipeline/.pixi/envs/eht-rebuild/bin/python"
OLD_YAML = Path(tempfile.mkdtemp(prefix="cleanup-eht-")) / "eht-official.yaml"
OLD_YAML.write_text(subprocess.run(["git", "-C", str(NEW), "show", f"{MAIN_COMMIT}:configs/eht-official.yaml"],
                                   check=True, capture_output=True, text=True).stdout)
LOCATIONS = ("input_root", "input_receipt", "feed_rotation_file", "official_script")


def replay(side, recorded, output):
    cli = recorded["cli"]
    argv = ["--stage", cli["stage"], "--day", str(cli["day"]), "--output-root", str(output), "--plan"]
    if cli.get("i_root"):
        argv += ["--i-root", cli["i_root"]]
    if cli.get("d_roots"):
        argv += ["--d-roots", *cli["d_roots"]]
    if cli.get("published_dterms"):
        argv += ["--published-dterms", cli["published_dterms"]]
    if cli.get("draws") is not None:
        argv += ["--draws", str(cli["draws"])]
    if cli.get("review_file") and Path(cli["review_file"]).is_file():
        argv += ["--review-file", cli["review_file"]]
    if side == "main":
        argv += ["--config", str(OLD_YAML)]
        module = "reconstruct_eht_official"
    else:
        for key in LOCATIONS:
            argv += ["--" + key.replace("_", "-"), recorded["parameters"][key]]
        argv += ["--no-env-file"]
        module = "eht_reconstruct"
    code = f"""
import {module} as driver
try:
    plan = driver.build_plan(driver.parser().parse_args({argv!r}))
    result = dict(ok=True, plan=plan)
except Exception as exc:
    result = dict(ok=False, error=f'{{type(exc).__name__}}: {{exc}}')
"""
    return run(side, code, python=PY)


def split_inputs(plan, root):
    """(code/config paths, data + upstream-stage inputs, installed-library inputs)."""
    code, other, library = [], {}, {}
    for entry in plan["inputs"]:
        path = entry["path"]
        if "/site-packages/" in path:
            library[path] = entry["sha256"]
        elif (path.startswith(str(root) + "/") or path.startswith(str(NEW) + "/") or path.startswith(str(MAIN) + "/")
                or path == str(OLD_YAML) or path.endswith("/eht-script-pipeline/configs/eht-official.yaml")):
            code.append(path)
        else:
            other[path] = entry["sha256"]
    return code, other, library


def comparable(plan):
    return {k: v for k, v in plan.items() if k not in ("cli", "inputs", "output_root")}


cases = []
for directory in sorted(ANALYSIS.glob("eht-paper-*")):
    recorded = json.loads((directory / "plan.json").read_text())
    output = Path(tempfile.mkdtemp(prefix="cleanup-eht-out-")) / "never-created"
    sides = {side: replay(side, recorded, output) for side in ("main", "new")}
    case = dict(stage_root=str(directory), stage=recorded["stage"], day=recorded["day"],
                main_ok=sides["main"]["ok"], new_ok=sides["new"]["ok"],
                errors={s: v.get("error") for s, v in sides.items() if not v["ok"]})
    if sides["new"]["ok"]:
        new = sides["new"]["plan"]
        new_code, new_other, new_lib = split_inputs(new, NEW)
        rec_code, rec_other, rec_lib = split_inputs(recorded, Path("/home/canoziia/MyProject/PHYS/eht-script-pipeline"))
        case["new_vs_recorded"] = dict(
            fields_identical=comparable(new) == comparable(json.loads(json.dumps(recorded))),
            field_differences=diff(comparable(recorded), comparable(new)),
            data_inputs_identical=new_other == rec_other,
            data_input_differences=sorted(set(new_other.items()) ^ set(rec_other.items()))[:20],
            library_paths_identical=set(new_lib) == set(rec_lib),
            library_files_with_changed_bytes_since_recording=sorted(
                p for p in set(new_lib) & set(rec_lib) if new_lib[p] != rec_lib[p]),
            code_inputs=dict(recorded=len(rec_code), new=len(new_code)))
        if sides["main"]["ok"]:
            main = sides["main"]["plan"]
            main_code, main_other, main_lib = split_inputs(main, MAIN)
            case["new_vs_main"] = dict(
                fields_identical=comparable(new) == comparable(main),
                field_differences=diff(comparable(main), comparable(new)),
                data_and_library_inputs_identical=new_other == main_other and new_lib == main_lib,
                input_differences=sorted(set(new_other.items()) ^ set(main_other.items()))[:20])
    else:
        # The parameter set alone is still derivable without replaying upstream verification.
        params = run("new", f"""
from eht_axion.experiment import read_experiment, eht_imaging_parameters
import json
recorded = json.load(open({str(directory / 'plan.json')!r}))
loc = {{k: recorded['parameters'][k] for k in {LOCATIONS!r}}}
result = eht_imaging_parameters(read_experiment(), loc)
""")
        case["parameters_vs_recorded"] = dict(
            identical=json.loads(json.dumps(params)) == recorded["parameters"],
            differences=diff(recorded["parameters"], json.loads(json.dumps(params))))
    cases.append(case)
    print(directory.name, case["main_ok"], case["new_ok"], case.get("new_vs_recorded", {}).get("fields_identical"),
          case.get("new_vs_main", {}).get("fields_identical"), case["errors"])

record = dict(
    purpose="eht-imaging plan (parameters, raw input, draws, label, scope) from config: new vs main vs recorded",
    compared_fields="all plan fields except cli/output_root; inputs compared except code/config files",
    cases=cases,
    all_replayed_new_plans_identical_to_recorded=all(
        c.get("new_vs_recorded", {}).get("fields_identical", False)
        and c.get("new_vs_recorded", {}).get("data_inputs_identical", False)
        for c in cases if c["new_ok"]),
    all_replayed_new_plans_identical_to_main=all(
        c.get("new_vs_main", {}).get("fields_identical", False)
        and c.get("new_vs_main", {}).get("data_and_library_inputs_identical", False)
        for c in cases if c["new_ok"]),
    all_parameters_identical=all(
        c["new_vs_recorded"]["fields_identical"] if c["new_ok"] else c["parameters_vs_recorded"]["identical"]
        for c in cases),
    notes=["cli differs by design: config path (experiment.yaml) and the four location options",
           "code inputs differ by design (renamed/removed modules)",
           "eht-paper-PD-* cannot be replayed by main or new code: their I-stage manifests list main-checkout "
           "source files that changed after recording (stale lineage); their parameter sets are compared instead",
           "ehtim files whose bytes changed in the shared .pixi environment after recording are listed per case; "
           "main and new replays see the same library bytes"])
write("eht-plan-equivalence.json", record)
print("replayed identical:", record["all_replayed_new_plans_identical_to_recorded"],
      "vs main:", record["all_replayed_new_plans_identical_to_main"],
      "parameters identical:", record["all_parameters_identical"])

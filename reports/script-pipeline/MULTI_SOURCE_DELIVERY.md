# Multi-source RAPTOR delivery

User authorized merging the completed remote and local400 jobs and converting
them to canonical Stokes NPY. `deliver_stokes()` now accepts explicit `job_roots`
and an optional independent `output_root`. The default single-run API remains.
Direct entrypoint: `scripts/merge_raptor_runs.py` / Pixi `merge-raptor` task.

Each planned task must appear in exactly one source run. Every source's actual
build, commit, grid/macros/science settings, successful task receipt, raw/table
inputs, HDF schema/coordinates/flux and declared-output lineage are checked.
No job/build is copied into another run or re-signed. Conversion products retain
real job/HDF ancestors; sources.json records per-task roots and inventory/index.json
lists the complete published delivery. Original HDF and inputs remain live.

Tests: **1075 passed in196.41s**, Ruff passed. New16-case synthetic suite spans
two compilers/binaries with equal science/core, bothresponse/background across
roots, exact17/1 state layouts, downstreambackground, livebothbuild ancestry,
missing/ambiguoustask, wrongbuildsettings/lineage and output/source containment.

The initial parallel test run overlapped newly written fixture files and produced
an intermediate-grid-name mismatch; rerunning frozen final files passed. A
pre-existing short affinity fakeworker could expire during severe host load;
the fake daemon now remains until its short parent process exits, without
changing production behavior or reducing receipt assertions.

Real inventory preflight:
- /home/canoziia/PHYS-production/results/m65-211:18234jobs
  (18189res100 +45res400).
- /home/canoziia/GRMHD_DATA/production/pilot-m65-t3000-res400-t16:96res400.
- Total18330tasks,0missing/duplicate.

Real conversion was started in tmux `eht-merge-convert`, targeting:
`/home/canoziia/GRMHD_DATA/production/m65-211-merged/data/`.
Expected8fullres100NPYs,15res400time-shardNPYs,total6858242944bytes(6.858GB)
before metadata. Eachres100shape(141,17,100,100,4),backgroundtotal(141,1,400,400,4).

Conversion first live-verifies all consumed rawbytes, including97GBsnapshotinputs;
no existing SHA receipt substitutes for current content validation. Hash session
reads shared inputs once per invocation. A separate temporary verification script
will compare all deliveredframe/channelvalues exactly againstoriginalHDF and
then validate the full inventory lineage; no normalization/rotation is added.

At this report creation the real job is **in progress, not declared complete**.
Completion requires `inventory/manifest.json` and root `ARRAY_VERIFIED.json` with
complete=true/all_IQUV_exact=true, plus exit0logs. Monitor report/logs will supply
actual completion status. No postprocess/inference or newRAPTORrun was started.

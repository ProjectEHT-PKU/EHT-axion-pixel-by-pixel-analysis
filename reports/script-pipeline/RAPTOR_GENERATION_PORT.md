# RAPTOR run generation port — implementation and verification

## Delivered

Python entrypoints `scripts/build_raptor.py` and `scripts/generate_raptor.py` are
available through the optional Linux `raptor-build` Pixi environment. They
reuse the pinned core's **official setup.sh and makefile**, not a duplicate
header generator or source list. The default analysis pipeline never launches
ray tracing. Detailed usage: `docs/raptor-generation.md`.

New modules in `eht_axion/producer/`: strict production configuration, lightweight
footer/table preflight, explicit job planning, official build, private job
execution/resume, HDF validation and atomic canonical Stokes delivery.

The config contains an explicit new migration benchmark, not silently revised
historical science: 6.2e9 Msun / 16.8 Mpc / camera46.2261, thermal/MUNIT3e25/R1/1,
cut1, low-five cloud spin.9375 / high-three .99, background spin.9375, 16 phases,
reference norm.1. Other analysis geometry values are copied, not certified as
new calibration. No original inputs or outputs were changed.

Alpha→omega delegates to existing `model.omega(alpha,2)`; no Dolan solver or
frequency digitization dependency is added. C alpha/omega inputs remain separate.
Changing the source of omega later must change this common model and recorded
identities; no orchestration rewrite is required. Existing cached data may not
be relabeled as recomputed at the new frequency.

## Core changes are engineering-only and separate commits

RAPTOR source remains a separate Git clone with retained history, local branch
`fix/axion-source`. `raptor-source.lock.json` pins:

`20056006820b6536c96de790dd5fdaad5ada4eae`.

- `9c6049e`: official setup option validation/failure exit, reliable template
  copying, grid case, explicit sigma/temperature-window/UNIF options; make
  compiler/include flag separation.
- `2005600`: compile recipe puts h5cc wrapper options before compiler flags,
  required by the resolved HDF5 2.2 wrapper.

Neither commit updates original plasma solvers/EOS/Faraday science. No push
was performed. Production requires the local clean locked checkout; no remote
clone/pull is assumed.

## Safety and integration

- Official setup receives explicit options in a private source/build snapshot.
  Final macros, source bytes, toolchain, flags and binary are recorded.
- Jobs use independent cwd/model/grid, explicit output slot, raw footer time,
  checked subprocess exit, timeout/cancel process groups, independent failure
  attempts and per-job atomic completion. No existence-only resume.
- Baseline identity excludes unconsumed cloud/table parameters; compute once
  per snapshot/resolution. Enabled jobs bind alpha/omega/table/spin.
- Hash verification context is invocation-only: shared raw sources are read
  once per batch, stat identity checked on reuse/exit; each new invocation fully
  rehashes. Active transaction roots are explicitly registered, not detected by
  arbitrary '.building-' substrings in input names. Inputs are treated as
  immutable; this is not a filesystem snapshot against adversarial mutation.
- HDF adapter checks all IQUV/tau/tauF/coordinate datasets, exact frequency,
  finite values, x-major block/pixel order, geometry, and Jy/sr area integral
  against spectrum Jy. No flips/transpose/U sign/area rescaling of delivered NPY.
- Response delivery is (T,17,100,100,4); background is (T,1,400,400,4).
  High-resolution phase jobs are not generated. Conversion uses existing
  convert_hdf5 and mmap, with raw HDF and receipts retained as ancestors.
- convert-only verifies the same production/build/parameter bindings as run,
  exact job build/raw/table dependencies and successful output receipts.
- Generated metadata cannot be removed to downgrade data into legacy mode;
  inference checks producer settings and source metadata independently. Real
  labels flow through fit; separate EHT/RAPTOR input roots have no silent fallback.
- .env is never executed by Python. CLI overrides exported EHT variables. The
  user's existing .env was not changed; additions are in .env.example.

## Validation performed

| Check | Result |
|---|---|
| Full default Pixi synthetic suite | **1001 passed in142.89s** |
| Ruff | passed |
| git diff --check | passed |
| pixi lock --check | passed |
| Default environment package records vs pre-port lock | exactly unchanged |
| Optional Linux compiler/HDF5/GSL/make/rsync environment | resolved and installed |
| Actual official setup→full C compile/link through Python/Pixi entrypoint | passed |
| Actual build reuse | verified CACHE HIT |
| Actual source/grid/macros/first baseline+response plan binding, no C launch | passed |
| Read-only real preflight:141 copied snapshot footers +2 existing radial CSV | passed |
| Full real plan | 18,048 responses+282 baselines=18,330 tasks; **not launched** |

Actual build (no scientific invocation) lives at:

`../eht-pipeline-runs/production/engineering-validation/build/`.

A separate build used for intermediate checking is retained under
`/home/canoziia/GRMHD_DATA/producer-engineering-validation/build/`.
Neither is promoted as a scientific dataset. Logs are archived alongside this
report in `raptor-generation-validation/`.

New tests exercise valid and malformed YAML/footer/table/model input, selected
source identity, setup failures/macro drift, compiler wrapper behavior, source
revision cache rejection, shared-input hash count, mutation, paths/symlinks,
fake process success/failure/partial/timeout/cancellation/resume/locks, HDF
coordinates/channels/flux, source/build task binding, missing generated metadata,
wrong conversion commit, producer delivery and actual analysis preprocess/fit/
background entrypoints. No test reads scientific GRMHD bodies or runs old binary.

An independent GPT-6.1 Sol read-only review identified convert-only binding,
metadata downgrade, temporary-name hash exclusion and active-macro parsing
holes; they were closed and regression-covered before final validation.

## Not done

No real RAPTOR snapshot, pilot or full transfer run. No new frequency/radial
solve, flux calibration, beam/crop/WCS recalibration or final limit rerun.
No old run deletion/re-signing, no output cleanup, no Git push.

Adding a Pixi environment changes the **whole lock hash**, and shared-source
changes also invalidate old stage identities, despite unchanged default package
records. This is intentional fail-closed behavior; old caches are not rewritten
or silently rerun. Use new production run IDs, and separately authorize a pilot
before deciding compute resources or scientific production.

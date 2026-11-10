# RAPTOR simulation

## Purpose

This page describes how the RAPTOR input is produced: building the pinned RAPTOR core,
planning the ray-tracing tasks, running them as content-addressed jobs, validating each
job's output, converting the outputs to Stokes cubes and delivering (or merging) those
cubes for the analysis. RAPTOR simulation does not read any EHT imaging product.

## Inputs

- A local checkout of the RAPTOR repository (`--source-root` or `EHT_RAPTOR_SOURCE_ROOT`)
  that is clean, is at exactly `raptor.source.commit`, and has a git remote pointing to
  `raptor.source.repository`.
- BHAC GRMHD snapshots `data<label>.dat` for the labels of `raptor.snapshots`
  (`--fluid-root` or `EHT_GRMHD_ROOT`).
- Two headerless CSV radial tables with columns $\alpha, r, R$ (`--radial-low`,
  `--radial-high` or `EHT_RADIAL_LOW`, `EHT_RADIAL_HIGH`).
- `experiment.yaml`: `shared` physics keys and the `raptor` section
  ([Experiment configuration](experiment-config.md)).

## Outputs

All outputs are below one run root (`--run-root` or `EHT_RAPTOR_RUN_ROOT`; the default
layout uses `raptor/runs/<name>/`, see [Data layout](data-layout.md#default-layout)):

- `build/`: the compiled RAPTOR and its receipts (stage `raptor:build`).
- `jobs/<task id>/`: one stage per task (stage `raptor:job`); failed attempts are kept in
  `attempts/<task id>/`.
- `data/raptor/res-100/alpha-<α>/stokes.npy`: one delivered cube per $\alpha$, and
  `data/raptor/res-400/background/t<A>-t<B>/stokes-t<A>-t<B>.npy`: axion-off time shards
  (stage `raptor:delivery`).
- `config.yaml` (a copy of the experiment file) and `requests/<uuid>.json` (one record per
  invocation).

`raptor_merge.py` writes the same `data/raptor/...` tree into a new output root, plus
`inventory/` (stage `raptor:merge`). File details are on
[Data layout](data-layout.md).

## How to run

Environment `raptor-build` (Linux only). Paths are taken from the command line, the
exported environment, `.env.raptor`, then `.env` ([Development](development.md)). With a data root (`M87_AXION_DATA_ROOT` or `--data-root`), every root option below has a
default from the [default layout](data-layout.md#default-layout); the commands give the roots
explicitly.

```bash
pixi run -e raptor-build raptor-build                       # build into RUN_ROOT/build
pixi run -e raptor-build raptor-generate --plan --plan-output plan.json   # list tasks only
pixi run -e raptor-build raptor-generate                    # run jobs, then deliver cubes
pixi run -e raptor-build raptor-generate --convert-only     # deliver already complete jobs
pixi run -e raptor-build raptor-merge --plan plan.json \
    --job-root RUN_A --job-root RUN_B --output-root MERGED  # merge runs of several hosts
```

Options of `raptor_generate.py`: `--max-jobs N` (run only the first N tasks; delivery waits
for the rest), `--threads` (OpenMP threads per RAPTOR process), `--workers` (concurrent
processes), `--cpu-affinity none|physical|logical`, `--cpu-list`, `--timeout-seconds`. Both
build and generate accept `--config`, `--source-root`, `--run-root`, `--env-file` and
`--no-env-file`.

## Details

### Build

`build_raptor` (`eht_axion/raptor/build.py`) checks that the source checkout is clean
(including untracked files), that HEAD is exactly `raptor.source.commit`, and that at least
one git remote URL equals `raptor.source.repository` after normalization (scheme, user,
port and a trailing `.git` removed, `git@host:` and `https://host/` treated alike, lower
case). It then copies the tracked source bytes into `build/source/`, and runs the
official `setup.sh` with arguments derived from `raptor.build` (`--code`, `--metric`,
`--int`, `--grid`, `--rad`, `--sf`, `--bflip`, `--sigma-cut`, `--uniform`, `--thetae-min`,
`--thetae-max`). It records the compiler command lines and the preprocessor macros of the
build (`macros.txt`) and verifies that the macros match the YAML settings. The build
record contains the settings, the RAPTOR repository and commit and the toolchain (the commit
identifies the source; no per-file source hashes) and the SHA-256 of the built binary, which every job
records as an input; radial tables and snapshots are not build inputs.

### Task plan

`plan_generation` (`eht_axion/raptor/plan.py`) reads only small inputs:

- For each label it reads the BHAC footer of `data<label>.dat` (little-endian block
  layout, time in $t_g$, adiabatic index, spin) and checks that the spin equals
  `raptor.fluid_spin` (0.9375),
  that the block size divides `raptor.grid.nxlone`, that the file is at least as long as
  the footer implies, and that consecutive footer times differ by
  `label_step_tg` × label difference.
- Each radial table is loaded, checked to be a rectangular increasing $(\alpha, r)$ grid
  covering the $\alpha$ values that use it, and hashed. Amplitudes are not rescaled.

For every snapshot the plan contains:

- one axion-off ("baseline") task at `response.resolution` (100) and one at
  `background.resolution` (400);
- for every $\alpha$ and every phase index $k = 0,\dots,15$ one axion-on ("response") task
  at 100 pixels with `AXION_ALPHA` $=\alpha$,
  `AXION_OMEGA` $=\Omega(\alpha)=\alpha(1-\alpha^2/(2n^2))$, `AXION_NORM` $=$
  `reference_norm`, `AXION_PHASE` $=k$ (the core applies the phase $k\pi/8$) and the radial
  table of `raptor.alpha_table`.

With the current values this gives $141 \times 130 = 18\,330$ tasks. Each task carries the
19 `model.in` values (`MBH` in $M_\odot$, `DISTANCE` in kpc, `M_UNIT` in g, `R_LOW`,
`R_HIGH`, `INCLINATION` in degrees, image size, camera size in $r_g$, one frequency,
`STEPSIZE`, the four axion values, `AXION_TABLE`, `MAX_LEVEL` = 1) and its output slot
(0 for baseline, $k+1$ for phase $k$).

### Task identifier

The task identifier is the SHA-256 of a canonical JSON document (`task_identity`) with:

- `schema_version` (2);
- the snapshot file name and its footer metadata;
- the 19 `model.in` values except the table path;
- for axion-on tasks, the table spin and the SHA-256 of the table content;
- `raptor.grid` and `raptor.build`.

Sky grid, EHT imaging and analysis settings never enter a task identifier, so changing
them never invalidates ray tracing.

### Jobs

`run_jobs` (`eht_axion/raptor/jobs.py`) runs the tasks with a fixed number of worker slots.
Before launching it verifies the build stage (including its lineage), checks that the
build's repository and commit equal the plan's `raptor_source` (which `raptor_generate.py`
copies from `raptor.source`), re-reads every
snapshot footer and radial table and compares them with the plan (radial-table hashes
always; snapshot hashes when the plan contains them), compares the build grid file and
macros with the plan, and checks memory headroom. The SHA-256 of every snapshot and table
is recorded in the manifest of each job that reads it. Each task runs in its own
stage directory `jobs/<task id>/` with inputs build manifest, binary, `grid_mks.in`,
`definitions.h`, snapshot and table. The run root may not lie inside a stage directory
(no ancestor may contain a `manifest.json`). RAPTOR is called as
`RAPTOR model.in <snapshot> <output slot>` in a temporary directory; if it exits
successfully, its output is validated (next section) and the stage is published. A failed
attempt is preserved under `attempts/<task id>/attempt-*/` (logs, receipt and the partial
working directory) and never published. Re-running resumes: a published job is reused only
if its manifest, success receipt, validated output metadata, `model.in`, task receipt and
grid file all match. Cancellation kills all process groups and keeps their attempts.

### Job output validation

`validate_job_output` (`eht_axion/raptor/outputs.py`) checks one job:

- the HDF5 file `output/img_data_<slot>.h5` contains exactly the datasets I, Q, U, V, tau,
  tauF, alpha, beta at the single frequency, each of shape
  $((N/10)^2, 100)$ (10×10-pixel blocks), with hard links and no external storage;
- every block's camera coordinates `alpha`, `beta` equal the expected pixel centres
  (x-major block and pixel order);
- the spectrum file `output/spectrum_<slot>_<inclination>.dat` has one row
  (frequency, I, Q, U, V in Jy) whose fluxes equal the HDF5 intensities (Jy/sr) times the
  pixel solid angle $\Delta\Omega = (\Delta x\,\Delta y)\,(r_g/D)^2$, using the constants of
  the build's `definitions.h`, to a relative tolerance of $10^{-10}$.

The result (`output.json`) also records counts of pixels with $I<0$ and with
$\sqrt{Q^2+U^2+V^2} > I$ as diagnostics; it makes no scientific acceptance.

### Conversion

`convert_hdf5` (`eht_axion/raptor/convert.py`) reads a rectangular time × state grid of
job HDF5 files and writes one float64 array of shape $(T, S, N, N, 4)$ (Stokes I, Q, U, V
in Jy/sr; axis 2 is camera x, axis 3 is camera y). Every 10×10 block is placed at its
recorded origin; the blocks must tile the image exactly.

### Delivery and merge

`deliver_stokes` (`eht_axion/raptor/delivery.py`) assigns every planned task to exactly one
job root, verifies each job stage and its lineage, checks that each job declares the build
and source inputs of its task, and checks the job's success receipt. The job HDF5 is bound
to the validated `output.json` by the SHA-256 values in the job manifest and is not
re-validated. Then it writes:

- for each $\alpha$: `res-100/alpha-<α>/stokes.npy` of shape (141, 17, 100, 100, 4); state 0
  is the axion-off baseline and state $k+1$ has axion phase $k\pi/8$;
- for the axion-off background: `res-400/background/t<A>-t<B>/stokes-t<A>-t<B>.npy`, shards
  of at most 10 snapshots, shape $(T, 1, 400, 400, 4)$.

Each product directory contains `generated.json`: the ray-tracing configuration, $\alpha$
(null for the background), resolution, snapshot labels, footer times in $t_g$, state count,
$\Omega$, axis order and units. The analysis accepts a product only if this record equals
what the current `experiment.yaml` would deliver (`generated_inputs`), and it takes the
frame cadence from the recorded footer times (`input_cadence_hours`).

`raptor_merge.py` calls the same delivery with several `--job-root` arguments and a new
`--output-root`, verifies the frozen plan's ray-tracing configuration, and writes an
inventory stage listing the products. Original job directories are not modified.

### CPU affinity

`allocate_cpu_sets` (`eht_axion/raptor/affinity.py`) gives each worker a disjoint set of
`threads` CPUs. `physical` uses one logical CPU per physical core (from
`/sys/devices/system/cpu/*/topology`), `logical` uses any allowed CPU, `none` does not pin.
Pinning uses `taskset`; `--cpu-list` restricts the allowed CPUs. With pinning, the CPU set
and OpenMP binding are part of the job identity. `check_memory` requires memory headroom
per worker (at least 2.5 GiB, and in `raptor_generate.py` twice the plan's per-process
lower bound), including cgroup-v2 limits.

### .env.raptor

`raptor_build.py` and `raptor_generate.py` read `.env.raptor` and `.env` (copy
`.env.raptor.example`). Keys: `EHT_RAPTOR_SOURCE_ROOT`, `EHT_GRMHD_ROOT`, `EHT_RADIAL_LOW`,
`EHT_RADIAL_HIGH`, `EHT_RAPTOR_RUN_ROOT`, `EHT_RAPTOR_THREADS` (default 2), `EHT_RAPTOR_WORKERS` (default 1),
`EHT_RAPTOR_CPU_AFFINITY` (default `none`), `EHT_RAPTOR_CPU_LIST`,
`EHT_RAPTOR_TIMEOUT_SECONDS`. The source root and the run root have no default. The run
root must not overlap the source, the fluid
directory or the table directories, and a run root keeps the experiment file it was
started with (`config.yaml`); a different file requires a new run.

## Related pages

- [Documentation index](index.md)
- [Overview](overview.md)
- [Experiment configuration](experiment-config.md)
- [Provenance](provenance.md)
- [Data layout](data-layout.md)
- [Development](development.md)

# RAPTOR Simulation

[← Documentation index](index.md)

The `RAPTOR` simulation pipeline compiles a pinned version of the polarized general-relativistic ray tracer `RAPTOR`, validates BHAC GRMHD snapshots and axion radial-profile tables, executes $18\,330$ ray-tracing jobs with bounded concurrency and optional CPU pinning, validates every HDF5/spectrum pair, and delivers float64 Stokes cubes for the downstream analysis. It does not depend on any EHT imaging product.

## Prerequisites and inputs

Running the `RAPTOR` scripts requires Linux (`linux-64`) and the `raptor-build` Pixi environment (`c-compiler`, `make`, `gsl`, `hdf5`, `pkg-config`, `rsync`). Four external inputs are required (located through the three roots `EHT_RAPTOR_SOURCE_ROOT`, `EHT_INPUT_ROOT`, and `EHT_DATA_ROOT` in `.env`, or explicit CLI flags; see [Configuration](configuration.md#machine-environment-variables-env)):

- **`RAPTOR` source checkout** (`--source-root` / `EHT_RAPTOR_SOURCE_ROOT`): a clean git checkout (no modified or untracked files) whose `HEAD` equals `raptor.source.commit` (`20056006820b6536c96de790dd5fdaad5ada4eae`) and that has at least one git remote normalizing to `raptor.source.repository` (`github.com/projecteht-pku/raptor`).
- **BHAC GRMHD snapshots** (`--fluid-root`, default `$EHT_INPUT_ROOT/grmhd`): binary files `data<label:06d>.dat` (`data001600.dat` through `data003000.dat` in steps of `10`, 141 snapshots). Resolved file paths must be single tokens of at most 255 UTF-8 bytes without whitespace.
- **Radial-profile tables** (`--radial-low`, `--radial-high`, default `$EHT_INPUT_ROOT/radial/result_a=0.9375_l1m1.txt` and `$EHT_INPUT_ROOT/radial/result_a=0.99_l1m1.txt`): two headerless three-column CSV files $(\alpha, r, R)$ for `low` (`spin = 0.9375`) and `high` (`spin = 0.99`), also at paths of at most 255 UTF-8 bytes without whitespace.
- **Experiment parameters** (`--config`, default `configs/experiment.yaml`): `shared` physics scalars and the `raptor` section.

All `RAPTOR` products live in `raptor/` below the data root (`plan.json`, `build/`, `jobs/`, `delivered/`); it must not overlap the repository root, the `RAPTOR` source checkout, the snapshot directory, or the radial-table directories, and cannot lie inside an existing stage directory.

## Running the simulation

Copy `.env.example` to `.env` and set the three roots `EHT_RAPTOR_SOURCE_ROOT`, `EHT_INPUT_ROOT`, and `EHT_DATA_ROOT`. The full workflow is:

```bash
# 1. Compile the pinned RAPTOR binary into $EHT_DATA_ROOT/raptor/build
pixi run -e raptor-build raptor-build

# 2. Validate inputs and list the tasks without launching C processes or writing files
pixi run -e raptor-build raptor-generate --plan

# 3. Write raptor/plan.json, execute the ray-tracing jobs (resumable) and deliver the Stokes cubes
pixi run -e raptor-build raptor-generate

# Deliver the cubes of already completed jobs without launching RAPTOR:
pixi run -e raptor-build raptor-generate --convert-only
```

**Execution options on `scripts/raptor/generate.py`:**
- `--max-jobs N`: run only the first `N` planned tasks (skips cube delivery while tasks remain incomplete).
- `--workers W` (`EHT_RAPTOR_WORKERS`, default `1`) and `--threads T` (`EHT_RAPTOR_THREADS`, default `2`): concurrent `RAPTOR` processes and OpenMP threads per process.
- `--cpu-affinity {none,physical,logical}` (`EHT_RAPTOR_CPU_AFFINITY`, fallback `none`) and `--cpu-list RANGE` (`EHT_RAPTOR_CPU_LIST`, e.g. `0-47`): worker CPU pinning via `taskset`.
- `--timeout-seconds SEC` (`EHT_RAPTOR_TIMEOUT_SECONDS`): per-process wall-time limit.

`scripts/raptor/build.py` needs only the `RAPTOR` checkout (`--source-root`) and the data root (`--data-root` / `EHT_DATA_ROOT`); it neither requires nor resolves the input root. The first `scripts/raptor/generate.py` run writes `raptor/plan.json` (hashing every snapshot once); every later run must plan the same tasks or is refused.

## Building the pinned RAPTOR core

`build_raptor` (`eht_axion/raptor/build.py`) publishes stage `raptor:build` at `raptor/build/`:

1. Verifies that `--source-root` is a clean git repository at `raptor.source.commit` with a normalized remote matching `raptor.source.repository`.
2. Extracts tracked git blobs (`git ls-tree` / `git cat-file`) into `build/source/`, sanitizes ambient `make`/include variables (`VPATH`, `CPATH`, `MAKEFLAGS`, `RAPTOR`, etc.), and resolves `h5cc` (`HDF5_CC`).
3. Runs `setup.sh` with `--code`, `--metric`, `--int`, `--grid`, `--rad`, `--sf`, `--bflip`, `--sigma-cut`, `--uniform`, `--thetae-min`, and `--thetae-max` derived from `raptor.build`.
4. Captures compiler command lines via `make -nB all` (`compiler-flags.txt`) and active preprocessor definitions via `h5cc ... -dM -E` (`macros.txt`), asserting that `metric`, `int_method`, `AMR`, `SMR`, `POL`, `SFC`, `BPOL`, `SIGMA_CUT`, `UNIF`, `THETAE_MIN`, and `THETAE_MAX` match `raptor.build`.
5. Writes `receipt.json` recording the toolchain, macro values, and the path of the `RAPTOR` executable, which downstream jobs declare as an input.

## Task planning and task IDs

`plan_generation` (`eht_axion/raptor/plan.py`) inspects only lightweight metadata and performs no ray tracing:

- **Snapshot footers** (`read_snapshot_footer`): reads the trailing 40-byte little-endian BHAC footer plus `nx` and `neqpar` records of each `data<label:06d>.dat`. Verifies 3D block layout, block dimensions `nx` dividing `raptor.grid.nxlone`, black-hole spin equal to `raptor.fluid_spin` (`0.9375` within $10^{-12}$), file size at least the footer-implied payload lower bound, constant adiabatic index $\gamma$, and uniform footer time steps equal to `shared.label_step_tg * step` ($1.0\,t_g$ per 10 labels).
- **Radial tables** (`validate_radial_table`): loads each headerless CSV, requires finite $(\alpha, r, R)$ with $R \ge 0$ arranged as a rectangular grid with strictly increasing $\alpha$ blocks sharing an identical strictly increasing $r$ grid ($\ge 2$ points), verifies $\alpha$ coverage for `raptor.alpha_table`, and hashes the file bytes (amplitudes are never rescaled).

For each of the 141 snapshots, the plan emits 130 tasks ($18\,330$ total):
- **2 `baseline` tasks** (axion off: `AXION_ALPHA = 0`, `AXION_OMEGA = 0`, `AXION_NORM = 0`, `AXION_PHASE = 0`, `AXION_TABLE = disabled`, `output_slot = 0`) at `response.resolution` (`100`) and `background.resolution` (`400`);
- **128 `response` tasks** at `100` pixels for all 8 values of $\alpha$ and 16 phases $k \in \{0, \dots, 15\}$ (`AXION_ALPHA` $= \alpha$, `AXION_OMEGA` $= \Omega(\alpha) = \alpha(1 - \alpha^2/(2n^2))$, `AXION_NORM` $= 0.1$, `AXION_PHASE` $= k$ corresponding to phase $k\pi/8$, `output_slot` $= k + 1$, and the radial CSV assigned by `raptor.alpha_table`).

Each task specifies the 19 `model.in` parameters (`MBH` in $M_\odot$, `DISTANCE` in $\mathrm{kpc}$, `M_UNIT` in $\mathrm{g}$, `R_LOW`, `R_HIGH`, `INCLINATION` in $\mathrm{deg}$, `IMG_WIDTH`, `IMG_HEIGHT`, `CAM_SIZE_X`, `CAM_SIZE_Y` in $r_g$, `FREQS_PER_DEC = 1`, `FREQ_MIN` in $\mathrm{Hz}$, `STEPSIZE`, `AXION_ALPHA`, `AXION_OMEGA`, `AXION_NORM`, `AXION_PHASE`, `AXION_TABLE`, `MAX_LEVEL = 1`).

**Task identifier (`task_id`):** the relative job directory named by the task parameters: `res-<n>/alpha-<α>/t<label>/phase-<k>` for an axion response (for example `res-100/alpha-0.30/t1600/phase-07`) and `res-<n>/baseline/t<label>` for an axion-off image (for example `res-400/baseline/t1600`). A job stage records its full `model.in` configuration; an existing job directory whose recorded configuration differs from the plan is refused, never silently reused.

## Job execution, CPU affinity, and output validation

Before launching jobs, `scripts/raptor/generate.py` and `run_jobs` (`eht_axion/raptor/jobs.py`):
- check free disk space ($\ge 4 \times \text{output\_stokes\_lower\_bound\_bytes} + 2\,\mathrm{GiB}$) and available RAM ($\ge \text{workers} \times \max(2.5\,\mathrm{GiB},\,2 \times \text{per\_process\_memory\_lower\_bound\_bytes}) + 2\,\mathrm{GiB}$, honoring Linux `MemAvailable` and cgroup-v2 `memory.max - memory.current`);
- verify the `raptor:build` stage, match `raptor_source`, `grid_mks.in`, and required C macros (`POL=1`, `num_frequencies=1`, `num_pixels_1d=10`, `NSPIN=3`, `SMR=1`, `AMR=0`, `UNIF=0`, `metric=MKSBHAC`, `DF=TH`, `IMGFILE=1`, `SPECFILE=1`), and recheck live snapshot footers, the radial-table SHA-256 recorded in the plan, and (when the plan records them) the snapshot SHA-256.

**Scheduling and CPU affinity (`eht_axion/raptor/affinity.py`):** task index `i` is deterministically assigned to worker slot `i % workers`. `allocate_cpu_sets` assigns each worker `threads` disjoint CPUs (`physical` selects one logical CPU per `(socket, core)` from `/sys/devices/system/cpu/*/topology`; `logical` selects any allowed CPU; `none` does not pin). When affinity is active, workers launch under `taskset --cpu-list ...` with `OMP_PROC_BIND=close` and `OMP_PLACES=cores`, and those binding settings are recorded in the `raptor:job` config so pinned cache hits remain deterministic on resume.

**Execution, resume, and failure preservation:** each new execution creates `raptor/attempts/<task_id>/attempt-*/` (`input.json`, `stdout.log`, `stderr.log`, `result.json`) and invokes `RAPTOR model.in <snapshot> <output_slot>` in a temporary stage directory. If `RAPTOR` fails, times out, or is interrupted (`SIGINT`/`SIGTERM` kills all active process groups), the partial temporary directory is copied to `attempts/<task_id>/attempt-*/work/` and no stage is published. On re-run, an existing `jobs/<task_id>/` stage is reused only if its manifest records the same configuration and inputs and its `success.json`, `output.json`, `model.in`, `task.json`, and `grid_mks.in` agree with the task.

**Per-job output validation (`validate_job_output` in `eht_axion/raptor/outputs.py`):**
- checks `output/img_data_<slot>.h5` has exclusively internal hard-linked datasets `{I,Q,U,V,tau,tauF}<freq:e>`, `alpha`, and `beta`, each finite of shape `((N/10)^2, 100)` ($10 \times 10$ x-major blocks);
- verifies block camera coordinates `alpha` and `beta` match expected pixel centres in $r_g$;
- verifies the single-row `output/spectrum_<slot>_<inclination:.2f>.dat` (`frequency`, $I, Q, U, V$ in $\mathrm{Jy}$) matches $\sum \text{HDF5}_{IQUV}\,\Delta\Omega$ within relative tolerance $10^{-10}$, where $\Delta\Omega = (\Delta x\,\Delta y)(r_g / D)^2$ uses the literal C constants `GGRAV`, `MSUN`, `SPEED_OF_LIGHT`, and `KPCTOCM` parsed from `build/definitions.h`;
- records pixel counts with $I < 0$ or $\sqrt{Q^2 + U^2 + V^2} > I$ in `output.json` as diagnostics without filtering data.

## Stokes cube conversion and delivery

`convert_hdf5` (`eht_axion/raptor/convert.py`) assembles a rectangular `time × state` grid of validated job HDF5 files into a float64 `(T, S, N, N, 4)` array (`[time, state, x, y, IQUV]` in $\mathrm{Jy\,sr^{-1}}$) and writes `convert.json`.

`deliver_stokes` (`eht_axion/raptor/delivery.py`) maps every planned task to its completed job stage in `raptor/jobs/`, checks each job's recorded build and inputs and its `success.json`, and publishes `raptor:delivery` stages under `raptor/delivered/`:
- `res-100/alpha-<alpha:.2f>/stokes.npy` of shape `(141, 17, 100, 100, 4)` (8 stages: state `0` axion-off baseline, state `k+1` phase $k\pi/8$);
- `res-400/background/stokes.npy` of shape `(141, 1, 400, 400, 4)` (axion-off only).

Each delivery stage also writes `generated.json` (schema `raptor-generated`, storing `ray_tracing_config`, `alpha`, `resolution`, `snapshot_labels`, BHAC footer `times_tg`, `state_count`, `omega`, axes, and units). Downstream analysis stages call `generated_inputs` to enforce exact agreement between `generated.json` and `configs/experiment.yaml`.


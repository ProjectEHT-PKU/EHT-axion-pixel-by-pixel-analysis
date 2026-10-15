# M87* EHT axion analysis — script pipeline

Research-specific Python scripts from RAPTOR Stokes images and EHT image ensembles
to six axion–photon coupling limit curves, with reusable angle/HDF5/plotting/numerical
modules. This workspace is **not an installable distribution or a generic CLI product**.
No Notebook, Nix/NixOS, system Python, LaTeX or old local intermediates are required.

**Start here: [Pixi 中文运行指南](docs/PIXI_QUICKSTART.md)**

**New data generation:** [RAPTOR生成端指南](docs/raptor-generation.md) — Python/Pixi
entrypoints reuse the official RAPTOR setup/make workflow. Generation is explicit;
`pixi run pipeline` never launches ray tracing. The core remains a separately
pinned repository; no duplicate header templates are maintained here.

## Run

Input data are already available on this machine through `data/inputs/`.
Outputs go directly into `../eht-pipeline-runs/` as configured in `.env.example`;
no extra run-name subdirectory is required. Verified caches are reused, and stale
configurations are rejected rather than silently overwritten.

```bash
pixi install --locked
pixi run check
# First use only; do not overwrite an existing .env.
cp .env.example .env
# Review .env; enable EHT_TRUST_PICKLE=true only for the trusted supplied inputs.
pixi run pipeline
```

Analysis-stage scripts automatically read the repository root's ignored `.env` on
each invocation; no manual shell loading is needed. Precedence is **CLI > exported
shell environment > env file**. Parsing uses `python-dotenv` with `interpolate=False`:
no shell execution or variable expansion. File-supplied relative paths use the
repository root; CLI and externally exported relative paths use the current working
directory. There is no current-directory env-file search.

Use `--env-file PATH` to replace the default file selection, or `--no-env-file` to
disable file loading without clearing the external environment (provide required
inputs explicitly via CLI for a self-contained run). The env-file selector itself
is relative to the current working directory.
The tracked `.env.example` contains analysis/shared settings only; producer scripts
read root `.env` then `.env.raptor`, with the latter overriding file values. For the
one-time producer template setup, see the [generation guide](docs/raptor-generation.md).
Use the unchanged `pixi run` tasks; no interactive shell activation is required.

For another machine, obtain the provided scientific dataset first; see the
[download appendix](docs/PIXI_QUICKSTART.md). Local data links and `.pixi/` are not
committed to Git. This checkout does not contain a second copy of the data archive.

## Code and data flow

```text
configs/          explicit physical scales and analysis choices

eht_axion/
  pipeline/       input conversion, RAPTOR/EHT processing, cache, inference, plots
  angles.py       circular statistics
  polarization.py aperture integration
  model.py        exposure and day-difference signal templates
  observations.py observation type and compatibility loaders
  limits.py       likelihood and posterior limits
  statistics.py   numerical probability utilities

tests/            synthetic physics, schemas, cache/safety and end-to-end tests
scripts/          direct processing/inference/plotting scripts and regression checks
reports/          small, versioned validation evidence; no raw arrays
```

```text
RAPTOR IQUV → processed Stokes → snapshot response → template distribution ─┐
EHT ensembles → matched topset → relocated Stokes/EVPA → covariance ──────┤
No-axion time series → background EVPA-difference scatter ────────────────┤
                                                                        ▼
                                  1D/2D × std/cor/cor_grmhd → figures/tables
```

Physical intermediate data are cached and checksummed. Large Stokes/EVPA arrays
remain mmap-capable NPY; response/template/limit products use named numeric NPZ
fields without pickle. Imported older products require explicit `--reuse-imported`.

| Mode | Observation noise | GRMHD background | Template uncertainty |
|---|---|---|---|
| `std` | diagonal of the same covariance used by cor | no | no |
| `cor` | full covariance | no | no |
| `cor_grmhd` | full covariance | yes | yes |

Scientific settings are YAML, separate from machine-local `.env` paths/options.
`configs/research.yaml` is the only scientific profile. Observed EVPA day differences
are always fixed at zero; observation uncertainties and covariance are retained.
There is no measured-day-difference mode or `observation` selector. This null scenario
is not automatically a median expected limit.
The profile uses intensity-weighted coherent 1D EVPA and ensemble-based statistics.
It is a documented research profile, not a frozen publication result.

## Documentation

- [顺序运行与常见问题](docs/PIXI_QUICKSTART.md)
- [Pipeline design and scientific assumptions](docs/SCRIPT_PIPELINE.md)
- [可复用模块与独立脚本](docs/REUSABLE_MODULES.md)
- [Canonical filenames, array schema and migration](docs/DATA_LAYOUT_PLAN.md)
- [RAPTOR stages and HDF5 conversion](docs/raptor-stage.md)
- [EHT processing](docs/eht-stage.md)
- [Transport / Figure-4 diagnostic branch](docs/transport-stage.md)
- [M87* constants and unresolved calibration](docs/M87_CONSTANTS_AUDIT.md)
- [Real-data validation evidence](reports/script-pipeline/README.md)

## Original project safety

The original project is **separate**:
`/home/canoziia/MyProject/PHYS/EHT_ZC`, branch `data/analysis-inputs-draft`.
Its notebooks, source code, simulations, data archives and old figures remain there.
This cleaned worktree deliberately omits those historical copies. The two worktrees
share Git object history but have separate files and indexes; do not remove `.git`.

The compatibility loaders and regression tests still used by this pipeline are
retained. The historical utility library and its optional statsmodels environment
are not needed here. All current tests run in the single default Pixi environment.

## Reusable modules without packaging

`eht_axion/` is a local source namespace, not a separately installed library. Scripts
import its angle, HDF5/image, plotting, statistics and cache functions directly.
`workflow.py` holds our data-specific run order, while argument parsing is in scripts/.
There is no `eht-pipeline` command, editable installation or package build metadata.
Dependencies and test/lint settings are all in `pixi.toml`; there is no
`pyproject.toml`, `pytest.ini` or `ruff.toml`. Use `pixi run test` / `pixi run lint`
so the explicit search paths and rules are applied. Existing task names stay familiar;
for example, `pixi run raptor` executes
`python scripts/process_raptor.py`.

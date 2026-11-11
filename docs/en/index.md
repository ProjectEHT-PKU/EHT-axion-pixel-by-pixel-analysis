# M87* EHT Axion Pipeline Documentation

English | [简体中文](../zh/index.md)

This repository implements an end-to-end search for axion-induced electric-vector position angle (EVPA) oscillations in Event Horizon Telescope (EHT) 2017 polarimetric observations of M87\*. It combines polarized general-relativistic ray-tracing simulations (`RAPTOR`), EHT Paper VII polarimetric imaging with D-term Monte Carlo, and phase-marginalized Bayesian inference on the dimensionless coupling $c$. The codebase is organized as direct CLI scripts managed by [Pixi](https://pixi.prefix.dev/) rather than an installable Python package.

## Quick start

### Requirements

- [Pixi](https://pixi.prefix.dev/) (provisions Python `3.12.*` and all locked scientific and build dependencies from `pixi.lock`) and `git`.
- **Platform support**: `default`, `eht-summary`, and `eht-imaging` support `linux-64`, `linux-aarch64`, `osx-64`, `osx-arm64`, and `win-64`. Compiling and running the C `RAPTOR` ray tracer (`raptor-build` environment) requires `linux-64`.
- **Clean git working tree**: writing any stage output requires zero uncommitted changes to tracked repository files (`git status --untracked-files=no`), because the repository's `HEAD` commit serves as the code provenance identity.

### Install environments and configure paths

Clone the repository and install the locked Pixi environments you need:

```bash
git clone <repository-url> eht-script-pipeline
cd eht-script-pipeline

pixi install --locked                  # default environment: analysis, verification, lint, core tests
pixi install --locked -e eht-imaging   # EHT imaging + full test suite (adds ehtim 1.3.2 and astropy)
pixi install --locked -e raptor-build  # Linux only: RAPTOR C toolchain and job orchestration
```

Copy `.env.example` to `.env` and set the three path roots of your machine:

```bash
cp .env.example .env
# Edit .env and set:
# EHT_RAPTOR_SOURCE_ROOT=/path/to/raptor
# EHT_INPUT_ROOT=/path/to/m87-axion-inputs
# EHT_DATA_ROOT=/path/to/m87-axion
```

With `EHT_DATA_ROOT` set in `.env` (or exported in the shell, or passed via `--data-root`), every script automatically resolves its input and output directories under the [default data layout](data-layout.md#default-layout), and every external input under the read-only [input layout](data-layout.md#input-layout) of `EHT_INPUT_ROOT`. All scientific parameters live in [`configs/experiment.yaml`](../../configs/experiment.yaml); machine paths live in the single `.env` (see [Configuration](configuration.md)).

### Run the tests and linter

All unit and integration tests run on synthetic data in temporary directories and require neither observational inputs nor a compiled `RAPTOR` binary:

```bash
pixi run --locked check                                       # default-env tests + Ruff linter
OPENBLAS_NUM_THREADS=2 pixi run --locked -e eht-imaging test  # full test suite (incl. ehtim & astropy)
pixi run --locked lint                                        # Ruff check only
```

### Where input data come from

`RAPTOR` simulation and EHT imaging are independent upstream pipelines whose finished products feed the analysis:

1. **Delivered `RAPTOR` Stokes cubes** (`$EHT_DATA_ROOT/raptor/delivered/data/`): either staged from an existing simulation archive or produced by running `raptor-build`, `raptor-generate`, and `raptor-merge` on BHAC GRMHD snapshots (`data001600.dat`–`data003000.dat`), two headerless $(\alpha, r, R)$ radial-profile CSVs, and a clean checkout of `github.com/projecteht-pku/raptor` at commit `20056006820b6536c96de790dd5fdaad5ada4eae` (located by `EHT_INPUT_ROOT` and `EHT_RAPTOR_SOURCE_ROOT` in `.env`; see [RAPTOR simulation](raptor-simulation.md)).
2. **Summarized EHT polarimetric ensemble** (`$EHT_DATA_ROOT/eht/ensemble/`): either staged from an existing imaging archive or produced by running `eht-reconstruct`, `eht-refine`, `eht-export-bundle`, `eht-sample`, and `eht-summarize` on the four official 2017 M87\* low-band HOPS `zbl-dtcal+selfcal` `uvfits` files (observing days `95`, `96`, `100`, `101`), feed-rotation JSON, and official EHT imaging script (located in `$EHT_INPUT_ROOT/eht/`; see [EHT imaging](eht-imaging.md)).

### Run the analysis chain on existing products

Once `$EHT_DATA_ROOT/raptor/delivered/data` and `$EHT_DATA_ROOT/eht/ensemble` exist, run the complete analysis and numerical sensitivity checks in the `default` environment:

```bash
pixi run analysis-templates
pixi run analysis-limits
pixi run analysis-check-grid --check grid --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --check prior --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance
pixi run analysis-summarize
```

Each command publishes a new stage directory under `$EHT_DATA_ROOT/analysis/` (`templates/`, `limits/`, `check-grid/`, `check-prior/`, `check-covariance/`, and `summary/`, including `limits/limits.csv` and `summary/RESULTS.md`) and refuses to overwrite an existing output directory.

To generate the publication figures, LaTeX macros (`numbers.tex`, `numbers.json`), and tables (`table-limits.tex`, `table-chi2.tex`) in `$EHT_DATA_ROOT/analysis/paper/` — which additionally reads `$EHT_DATA_ROOT/raptor/plan.json`, `src/axion.c` of the `RAPTOR` checkout at `$EHT_RAPTOR_SOURCE_ROOT`, and the radial-profile CSVs in `$EHT_INPUT_ROOT/radial/` — run:

```bash
pixi run analysis-paper
```

### Verify results

Every published stage contains a `manifest.json` recording parameters, the repository `git_commit`, and SHA-256 checksums of all declared inputs and outputs. Verify a stage (and with `--recursive`, its entire upstream lineage) without recomputing any science products:

```bash
pixi run verify $EHT_DATA_ROOT/analysis/limits --recursive
```

## Documentation

| Guide | Contents |
|---|---|
| [Overview and concepts](overview.md) | Pipeline architecture, Mermaid data-flow diagram, physical scales, polarimetric notation, and core design rules |
| [Configuration](configuration.md) | Complete key-by-key reference for `configs/experiment.yaml`, validation and derived views, and the three path roots and execution settings of `.env` |
| [RAPTOR simulation](raptor-simulation.md) | Pinned `RAPTOR` build, BHAC snapshot and radial-table validation, content-addressed tasks, job runner, CPU affinity, delivery, and multi-host merge |
| [EHT imaging](eht-imaging.md) | EHT Paper VII polarimetric imaging reproduction, convex PDHG refinement, D-term Monte Carlo, ensemble summary, and deviations from Paper VII |
| [Analysis](analysis.md) | Image preprocessing, sector aperture operator, harmonic templates, feature selection, covariance modes, Bayesian likelihood/posterior, checks, and paper outputs |
| [Data layout and provenance](data-layout.md) | Default `$EHT_DATA_ROOT` tree, read-only `$EHT_INPUT_ROOT` input layout, per-stage file/array specifications, atomic `stage_directory` manifests, non-manifest receipts, and verification |
| [Development](development.md) | Pixi environments and tasks, test suite and linter, module/script map, architectural invariants, and guides for extending parameters or stages |

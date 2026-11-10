# Development

## Purpose

This page describes the software environment of the repository: the Pixi environments and
tasks, how to run the tests and the linter, how machine-local settings are given through
`.env` files, and the rules that new code must follow. The repository is not an installable
package; scripts are run directly or through `pixi run`.

## Inputs

- `pixi.toml` and `pixi.lock`: environments and tasks.
- `.env` and `.env.raptor`: machine-local paths and execution settings, copied from
  `.env.example` and `.env.raptor.example` (both files are git-ignored).

## Outputs

None; tests write only to temporary directories.

## How to run

```bash
pixi install --locked
pixi run check                       # tests and lint in the default environment
pixi run -e eht-imaging test         # all tests, including those that need ehtim and astropy
pixi run lint
pixi run verify PATH/TO/STAGE --recursive
```

## Details

### Environments

| Environment | Adds | Used for |
|---|---|---|
| `default` | Python 3.12, numpy, scipy, matplotlib, h5py, PyYAML, pytest, ruff, python-dotenv | analysis scripts, tests, lint, `verify` |
| `raptor-build` | C compiler, make, GSL, HDF5, pkg-config, rsync (linux-64 only) | `raptor-build`, `raptor-generate`, `raptor-merge` |
| `eht-summary` | astropy | `eht-summarize` |
| `eht-imaging` | astropy and ehtim 1.3.2 | `eht-reconstruct`, `eht-refine`, `eht-export-bundle`, `eht-sample` |

All environments set `MPLBACKEND=Agg`, `OPENBLAS_NUM_THREADS=2`, `OMP_NUM_THREADS=2` and
`PYTHONUNBUFFERED=1`.

### Tasks

| Task | Command |
|---|---|
| `raptor-build`, `raptor-generate`, `raptor-merge` | `scripts/raptor/*.py` |
| `eht-reconstruct`, `eht-refine`, `eht-export-bundle`, `eht-sample`, `eht-summarize` | `scripts/eht/*.py` |
| `analysis-templates`, `analysis-limits`, `analysis-check-grid`, `analysis-check-covariance`, `analysis-summarize`, `analysis-paper` | `scripts/analysis/*.py` |
| `verify` | `scripts/verify_products.py` |
| `test` | `pytest -q -ra` with `testpaths=tests` and `pythonpath=.` |
| `lint` | `ruff check --isolated --target-version py312 --line-length 110 --select E,F,W eht_axion tests scripts` |
| `check` | `test` and `lint` |

### Tests

The tests use synthetic inputs only; they never read observational data or run RAPTOR. In
the default environment the tests that need ehtim or astropy are skipped; the
`eht-imaging` environment runs all of them. Some tests are not unit tests of a function but
checks of repository rules:

- `test_experiment.py` checks validation of every section of `configs/experiment.yaml` and that no
  public function in `eht_axion` has a numeric default value (a short allow-list covers
  execution settings and array-layout arguments).
- `test_docs_links.py` checks both documentation trees, `docs/en/` and `docs/zh/`: every
  page is reachable from `README.md` through relative links, every relative link and anchor
  resolves, every page has one H1 title and links back to the `index.md` of its language,
  and the two trees have the same file names and the same heading levels per page.

### .env files

`eht_axion/environment.py` reads dotenv files as data: no shell expansion, no variable
interpolation, no search in the caller's working directory. Relative paths are taken
relative to the repository root. Unknown keys and duplicate keys are errors. The order of
precedence is: command-line option, exported environment variable, `.env.raptor` (only for
`scripts/raptor/build.py` and `scripts/raptor/generate.py`), `.env`. `--env-file PATH` replaces the default
files and `--no-env-file` disables them.

| Key | Read by | Meaning |
|---|---|---|
| `EHT_DATA_ROOT` | every script | data root of the default layout ([Data layout](data-layout.md#default-layout)); overridden by `--data-root` |
| `EHT_IMAGING_INPUT_ROOT` | `scripts/eht/reconstruct.py` | directory of the HOPS uvfits files |
| `EHT_IMAGING_INPUT_RECEIPT` | `scripts/eht/reconstruct.py` | input receipt JSON |
| `EHT_IMAGING_FEED_ROTATION` | `scripts/eht/reconstruct.py` | feed-rotation table JSON |
| `EHT_IMAGING_OFFICIAL_SCRIPT` | `scripts/eht/reconstruct.py` | official EHT imaging script |
| `EHT_RAPTOR_SOURCE_ROOT`, `EHT_GRMHD_ROOT`, `EHT_RADIAL_LOW`, `EHT_RADIAL_HIGH`, `EHT_RAPTOR_RUN_ROOT`, `EHT_RAPTOR_THREADS`, `EHT_RAPTOR_WORKERS`, `EHT_RAPTOR_CPU_AFFINITY`, `EHT_RAPTOR_CPU_LIST`, `EHT_RAPTOR_TIMEOUT_SECONDS` | RAPTOR scripts | see [RAPTOR simulation](raptor-simulation.md) |

The analysis scripts, `scripts/eht/refine.py`, `scripts/eht/export_bundle.py`, `scripts/eht/sample.py` and
`scripts/eht/summarize.py` read only `EHT_DATA_ROOT` from `.env`; their other paths are
command-line options or defaults below the data root.

### Rules for code changes

- A scientific value belongs in `configs/experiment.yaml`, once, with validation in
  `eht_axion/experiment.py`; library functions receive it as an explicit argument.
- Machine paths and execution settings belong on the command line or in `.env`.
- A script that writes scientific products uses `stage_directory` and declares every data
  file it reads from outside the repository as an input; repository files are identified
  by the commit and are never inputs ([Provenance](provenance.md)). Commit before running:
  results are refused while tracked files have uncommitted changes.
- A product is read only in the format the current code writes; a changed product needs a
  new output root. Scripts refuse existing output roots except where
  resuming is part of the design ([Provenance](provenance.md#output-roots)).

## Related pages

- [Documentation index](index.md)
- [Package map](package-map.md)
- [Provenance](provenance.md)
- [Experiment configuration](experiment-config.md)
- [RAPTOR simulation](raptor-simulation.md)

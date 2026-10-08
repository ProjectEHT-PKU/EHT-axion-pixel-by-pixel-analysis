# EHT axion pixel-by-pixel analysis

## Script-first development pipeline

This branch provides a complete **Python-only** image-output → inference → figure
pipeline. No notebook execution or author-local intermediate paths are required.

**Start here: [Pixi quick start (中文)](docs/PIXI_QUICKSTART.md) · [Script pipeline guide](docs/SCRIPT_PIPELINE.md)**

Install [Pixi](https://pixi.prefix.dev/latest/installation/) first. Nix/NixOS and a
preinstalled Python are not required. Dependencies and platforms are locked in `pixi.lock`.

```bash
pixi install --locked
pixi run check
pixi run pipeline --config configs/research.toml \
  --input-root data/inputs \
  --output-root /path/to/large-disk/research-run --trust-pickle
```

Create `data/inputs` from the published archive without copying large arrays:

```bash
pixi run migrate inputs --source /path/to/shared_data/analysis-inputs/draft --destination data/inputs
```

Dates use `2017-04-05`, alpha identities use `alpha-0.30`, response/limit products
use named NPZ fields. `draft` is publication status only, not a runtime requirement.
See [data naming and migration](docs/DATA_LAYOUT_PLAN.md).

Physical cache layers: processed Stokes → snapshot response fits → template moments;
EHT identity-matched topset → observation covariance; measured GRMHD background →
six limit curves → PNG/PDF/SVG/CSV. Completed stages are content-verified and reusable.

The default profile restores explicit **intensity-weighted** 1D apertures, uses
ensemble-centered covariance and fixes the five-day phase advance. This is a new
research profile, **not** a one-bug replacement for previous figures. Pending
physical choices, literature constants and provenance are documented in
[M87_CONSTANTS_AUDIT.md](M87_CONSTANTS_AUDIT.md).

The remainder below documents the previous runner/data-release workflow retained
for reference and regression. Prefer the script pipeline above for new work.


This repository evaluates limits on an axion-induced, time-dependent rotation of
M87* EVPA using EHT image ensembles and GRMHD/RAPTOR signal templates.  It
contains two analyses:

- **2-D / pixel-by-pixel:** retain spatial EVPA pixels and their covariance;
- **1-D / angular:** use coherent aperture EVPA and spaced angular bins (spacing alone
  does not establish statistical independence). The previous rebuilt profile used
  direct Q/U flux integration; the new script profile defaults to intensity weighting.

Each analysis reports three uncertainty models:

| mode | observation uncertainty | GRMHD background | axion-template uncertainty |
|---|---|---|---|
| `std` | diagonal standard deviations | no | no |
| `cor` | full image-ensemble covariance | no | no |
| `cor_grmhd` | full covariance | yes | yes |

## Installation

```bash
python -m pip install -e '.[dev]'
```

## Reproducible final-limit runner

The historical workflow lives in `GRMHD/main.ipynb`,
`EHT-Imaging/main.ipynb`, and `result.ipynb`.  The tested command-line path is:

```bash
export PYTHONPATH=$PWD
python scripts/inspect_intermediates.py
python scripts/run_final_limits.py --legacy --dimension 2d --mode cor_grmhd

# Validate/rebuild corrected EHT products before the corrected final run.
python scripts/validate_topset_ranking.py
python scripts/rebuild_corrected_topset.py
python scripts/rebuild_1d_products.py
python scripts/fit_rebuilt_1d_products.py
python scripts/run_final_limits.py --dimension 2d --mode cor_grmhd \
  --observation-dir audit/rebuilt_2d
python scripts/run_final_limits.py --dimension 1d --mode cor_grmhd \
  --observation-dir audit/rebuilt_1d \
  --fit-parameter-dir audit/rebuilt_1d
```

Run all six combinations by selecting `--dimension {1d,2d}` and
`--mode {std,cor,cor_grmhd}`.  `--legacy` deliberately reproduces the historical
Notebook formulas; without it the implementation uses:

- 16 unique nuisance phases (no duplicated 0/2π endpoint),
- products of independent epoch-pair likelihoods in log space,
- the Gaussian log determinant when covariance depends on coupling,
- explicit trapezoidal posterior integration and an explicit prior,
- circular interpolation of fitted phases,
- the appropriate model-error phase at both epochs.

The default output is under `audit/results/`.  Generate the six-curve result with:

```bash
python scripts/plot_final_limits.py
```

See [`audit/REPORT.md`](audit/REPORT.md) for reproduced values, verified fixes,
runtime measurements, and remaining scientific choices.

## Shared final-inference data

Collaborator-facing binary inputs live only under
`shared_data/final-inference/draft/` and are managed by Git LFS. The `draft`
revision is intentionally mutable while the analysis is under active review;
immutable `v1` will be created only when the paper result is frozen.

After `git lfs pull`, run directly from the shared dataset:

```bash
python scripts/run_final_limits.py --dimension 2d --mode cor_grmhd \
  --observation-dir shared_data/final-inference/draft/2d \
  --fit-parameter-dir shared_data/final-inference/draft/2d
python scripts/run_final_limits.py --dimension 1d --mode cor_grmhd \
  --observation-dir shared_data/final-inference/draft/1d \
  --fit-parameter-dir shared_data/final-inference/draft/1d
```

See [`shared_data/README.md`](shared_data/README.md) for update, freeze, and
partial-download policy.

A separate, larger development dataset is available under
`shared_data/analysis-inputs/draft/` (about 20.9 GB via Git LFS). It contains the
RAPTOR Stokes cubes and four EHT-Imaging 4800-image ensembles required to
reconstruct the full downstream analysis. It starts *after* RAPTOR and
EHT-Imaging have run; it is not raw GRMHD/visibility input. See
[`shared_data/analysis-inputs/README.md`](shared_data/analysis-inputs/README.md).

## Intermediate data

Raw and generated arrays are intentionally not stored in ordinary Git.  The
final calculation expects the following local products:

```text
GRMHD/data_out/EVPA_cos_fit_result_alpha={alpha}.npy
GRMHD/data_out/EVPA_cos_fit_result_1d_10deg_alpha={alpha}.npy
GRMHD/data_out/{A_params,phi0_params}{,_1d}_dict.npy
EHT-Imaging/data_out/{95,96,100,101}_all_data.npy
EHT-Imaging/data_out/{95,96,100,101}.npz
EHT-Imaging/data_out/{95,96,100,101}_EVPA_1d.npy
EHT-Imaging/data_out/{95,96,100,101}_1d.npz
EHT-Imaging/data_out/condition.npy
```

The EHT raw image ensembles are read from:

```text
EHT-Imaging/data_in/230_{95,96,100,101}.pkl
```

Large RAPTOR ray tables `GRMHD/data_in/GRMHD_cut*` are raw CSV trajectories.
They are reduced to small `r_coord*.npy` and `axion_field_coord*.npy` products and
must remain outside Git.

## Tests

```bash
export PYTHONPATH=$PWD
pytest -q
```

Tests cover EVPA circular arithmetic, phase interpolation, Stokes integration,
Gaussian likelihood normalization, independent likelihood combination, posterior
integration, and repaired legacy utility behavior.

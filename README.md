# EHT axion pixel-by-pixel analysis

This repository evaluates limits on an axion-induced, time-dependent rotation of
M87* EVPA using EHT image ensembles and GRMHD/RAPTOR signal templates.  It
contains two analyses:

- **2-D / pixel-by-pixel:** retain spatial EVPA pixels and their covariance;
- **1-D / angular:** radially integrate Stokes Q/U and use independent angular bins.

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

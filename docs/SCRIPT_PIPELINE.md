# Script-only research pipeline

> Development branch: `refactor/script-pipeline`. New runs do not execute or parse notebooks.
> Previous notebooks remain as historical provenance, not runtime dependencies.
> Scientific choices pending author confirmation are labeled, not silently declared fixed.

## Structure

```text
configs/research.toml         physical scales + preserved image settings + analysis choices

eht_axion/
  pipeline/
    workflow.py              research-stage orchestration (no argument parser/entrypoint)
    common.py                SHA-256 identity, locks, transactional/cache manifests
    constants.py             M/D → time, angular scale, axion mass conversion
    raptor.py                HDF5 conversion, image processing, response fits, background
    eht.py                   candidate identity matching, ranking, relocation, summaries
    products.py              named NPZ schema, centered observation covariance
    inference.py             six-mode batch, effective settings, numerical diagnostics
    figures.py               figures/tables, no fitting or hidden data paths
    diagnostics.py           intermediate EHT/mask and GRMHD-lag visual checks
    convergence.py           denser coupling/phase grids, reuse all physical data
    transport.py             Figure-4 transport decomposition diagnostic branch
  polarization.py            coherent aperture operators and integration
  model.py                   exposure/day-difference templates
  limits.py                  likelihood and posterior upper limits
  statistics.py              Gaussian and quadrature utilities

scripts/                     direct scripts; Pixi tasks execute these files
  process_raptor.py           RAPTOR processing stage
  process_eht.py              EHT processing stage
  estimate_background.py     background stage
  build_observations.py      observation products
  calculate_limits.py        six limit curves
  plot_limits.py             plot cached run results
  run_all.py                 run these stages in sequence
  raptor_tools.py            explicit HDF5/single-stage tools

tests/test_pipeline_*.py      synthetic scientific and failure-mode tests
```

Old notebooks, utility libraries, audit outputs and data archives are not copied into
this cleaned worktree. They remain intact in the separate original `EHT_ZC` worktree
(branch `data/analysis-inputs-draft`, baseline `74ebae1`) and its Git history.
Only active comparison/regression entrypoints and tested compatibility loaders remain.

## Run with Pixi (recommended; no Nix or system Python required)

Install Pixi using the [official instructions](https://pixi.prefix.dev/latest/installation/),
then open a terminal in the repository root:

```bash
pixi install --locked
pixi run help
pixi run check
cp .env.example .env
# Review paths and enable EHT_TRUST_PICKLE only for trusted inputs.
set -a
source .env
set +a
pixi run pipeline
```

`pixi.toml` defines Python, NumPy, SciPy, Matplotlib, HDF5 and test/lint dependencies.
`pixi.lock` fixes the resolved environment. Pixi does **not** install this project.
Tasks execute `python scripts/*.py`; scripts load the adjacent reusable source modules
through a small repository-relative bootstrap. No pip install, editable distribution,
console entrypoint or manual PYTHONPATH is needed. Dependencies and explicit pytest/Ruff
options are unified in `pixi.toml`; use `pixi run test` and `pixi run lint` rather than
bare tool commands. No separate packaging/test/lint config files remain.
Platforms are Linux x86-64/ARM64,
macOS Intel/Apple Silicon and Windows x86-64. The actual execution validation reported
here is Linux x86-64; a solved lockfile for another platform is not an execution test.

No copying of raw cubes into the repository is necessary. Inputs are opened read-only.
One default Pixi environment runs all current tests. The removed historical utility
library needed statsmodels; current processing and regression do not.

Pixi activation sets Agg rendering and two BLAS/OpenMP threads. Run long computations
in tmux on Linux/macOS or a suitable persistent terminal/job scheduler. Windows users
can run individual `pixi run` commands directly without Bash.

Pickle is executable: only enable `--trust-pickle` for inputs whose source you trust.
SHA-256 detects corruption/changes but does not establish trust.

Main research-stage scripts use command-line arguments over exported `EHT_*`
environment values. Python does not auto-load `.env`; source a trusted local file
in the shell as shown above. Destructive/reuse options default false and can be
explicitly disabled with `--no-rebuild`/`--no-reuse-imported`. Numerical choices
remain in the TOML science config, not environment variables.

Individual stages use the same exported configuration:

```bash
pixi run raptor
pixi run eht
pixi run background
pixi run observations
pixi run inference
pixi run figures
pixi run verify /path/to/run-name/eht
```

Use explicit config/input/output arguments only for temporary overrides.
Expert standalone tools (conversion, migration, direct plotting) keep their own explicit
arguments; `.env` defaults apply to the six main stages and `run_all.py`.
The whole-chain order is RAPTOR → EHT → background → observations → inference → figures.
RAPTOR and EHT may be computed independently, in different stage directories.

## Naming and migration

The canonical input root is `data/inputs`, not a publication-status path:
`raptor/res-100/alpha-0.30/stokes.npy`, `eht/2017-04-05/image-ensemble.pkl`.
Original shared `analysis-inputs/draft` is supported as a legacy read adapter.
See [DATA_LAYOUT_PLAN.md](DATA_LAYOUT_PLAN.md) for the complete schema.

```bash
pixi run migrate inputs --source /path/to/analysis-inputs/draft --destination data/inputs
pixi run migrate run --source /path/to/old-run --destination /path/to/canonical-run
```

Large arrays are linked, not duplicated. Imported manifests preserve original code
identity and mark scientific_recomputation=false. Explicit `--reuse-imported` is required
to reuse these stages. Config/input changes still fail validation. Do not delete source
archives that a migrated symlink view depends on. New science runs use the canonical
names directly and do not require migration.

## Physical cache layers

| Stage directory | Scientific meaning | Main outputs | Recompute when |
|---|---|---|---|
| `raptor/processed/alpha-*` | simulated telescope-plane Stokes after blur/rotation/crop | `stokes.npy` `(T,17,64,64,4)` | image geometry, blur or input changes |
| `raptor/response/alpha-*` | response to axion phase for every GRMHD snapshot | `response-2d.npz`, `response-1d.npz` with named amplitude/phase/baseline/snapshot_labels | phase fitting or 1D observable changes |
| `raptor/templates/{1d,2d}` | temporal distribution of responses | `templates.npz` mean/std and valid counts | response distribution model changes |
| `eht/` (ISO date subdirectories) | identity-matched topset and observed EVPA ensembles | `selection.npz`, `YYYY-MM-DD/{stokes.npy,evpa-1d.npy,summary.npz}`, mask | ranking, registration, apertures or inputs change |
| `raptor/background` | high-resolution-origin processed Stokes | `stokes.npy` | background image preprocessing changes |
| `background` | pooled no-axion EVPA difference scatter at recorded lag | `summary.json`, `evpa-difference.npz` | lag, region or fit changes |
| `observations/{1d,2d}` | selected day differences and centered covariance | `observations.npz`, diagnostics JSON | mask, selected angular offset or covariance model changes |
| `inference` | posterior constraints and integration diagnostics | `limits-{1d,2d}.npz` with mode axis, `posterior-{mode}-{dim}.json`, `run.json` | timing, coupling/phase grid, prior or noise model changes |
| `figures` | views of the saved numerical results | PNG/PDF/SVG, CSV, figure JSON | plot styling changes |

Only raw EHT inputs use pickle. New science products are ordinary numeric NPZ/NPY,
loadable with `allow_pickle=False`. Never infer stage completeness from the existence
of a single array: each stage's `manifest.json` is written last. Old `stage.json`
is accepted for read-only compatibility. Root run manifests are descriptive, not completion markers.

The research profile currently requires 360 one-degree bins; arbitrary angular grids
are rejected rather than silently misaligning RAPTOR and EHT angle origins.

The stage signature includes input bytes, settings and relevant Python source bytes.
A cache hit verifies output checksums. Changed inputs/config/source produce a stale-cache
error, not silent reuse. Use a new output root or `--rebuild`. Rebuild creates a temporary
sibling stage and replaces only after success; an ordinary failure leaves old output intact.
A process crash between renames can leave `.previous`/`.building-*`; inspect and recover,
never automatically treat these as completed products. Locks record PID for manual recovery.

## Explicit methodological changes vs previous code

These changes are separate, not a claim that all historical definitions were bugs:

- **Formula fix:** second day-pair advances by omega × 120h, not omega × 60h.
- **Numerics:** centered 10° apertures, exact radial endpoint and trapezoidal endpoint weights.
- **Definition restored:** default `intensity` uses integral I(Q+iU), consistent with the old
  intensity-weighted statistic and EHT-VII. `flux` is an explicit alternative, not a silent correction.
- **Estimator choice:** new profile uses local circular lifting + Gaussian population MLE,
  not histogram-fit main-peak estimates. Concentration/finite diagnostics are retained;
  Gaussianity is not guaranteed merely by finite values. Topset membership can change.
- **Coherent covariance:** observed center and covariance refer to the same ensemble lift;
  std uses diag(cor). A std/cor comparison now isolates spatial association, not different marginals.
- **Background choice:** re-estimated from state0, zeros retained; lag3 is about26.68h, not exact24h.
  Gaussian MLE and historical histogram fitting are distinct explicit options in the background CLI.
- **Statistics:** independent day-pair likelihoods multiply; parameter-dependent logdet retained;
  prior/bounds explicit; phase endpoint unique; posterior segment inversion matches quadrature.
- **Performance:** fixed-covariance likelihood quadratic coefficients are cached; physical
  intermediate arrays are streamed rather than loading the res400 cube into memory.

## Choices intentionally preserved, not externally certified

`46.2262/32.6868 Rg`, sigma `2.16871 Rg`, rotation108°, fixed two-stage EHT offsets,
manual exclusions, 1D angular offset4 and no inheritance of the 2D mask are retained
and recorded. Original RAPTOR calibration/FITS headers are needed to settle angular
FOV and Stokes basis. Mass6.5e9/distance16.8 are a declared working convention, not
proof that the original transfer simulation used exactly that calibration.

Principal quantum number4 matches the historical frequency formula but is not asserted
to be the correct interpretation of the draft's211 mode. Axion template shape is held
fixed when analytic temporal frequency is evaluated; changing the physical state requires
checking the original simulated templates, not just editing a frequency.

In this implementation the plotted `c` is the coefficient multiplying the unit
response **before** `b(alpha)`; following the current code convention it is interpreted
as `2*pi*g_agamma*f_a`. Predicted signal scales as `c*b(alpha)*A/AXION_NORM` times
the time-difference response. The old draft sometimes defined c including b; do not
compare those vertical axes without converting the normalization.

The full noise mode assumes spatially diagonal model errors and independent model
realizations between day-pairs. Cross-day observation noise is also assumed independent.
These are different assumptions. GRMHD background is scalar diagonal and pooled over
space/time. Shared latent templates would require cross-pair covariance not implemented
by this profile. Prior/phase/coupling convergence and these assumptions need paper review.

The figure is a reproducible **research-profile result**, not a frozen publication claim.
Do not compare its shifts with old results as a one-bug effect: observable, estimator,
selection and time-model changes are all listed above.

## HDF5 start and missing upstream inputs

`pixi run convert-hdf5 --help` documents explicit named Stokes
datasets, snapshot/state association and block layout. No dependence on dictionary key
order or absolute labels used as array indices. HDF5 attributes should be archived as
provenance. Conversion is tested synthetically; the current local dataset begins at
NumPy cubes, so actual historical HDF5 conversion is not claimed verified against source.

Raw GRMHD fluid snapshots, original RAPTOR execution configuration and EHT visibility
reconstruction configuration remain outside this repository. “End-to-end” here means
RAPTOR image outputs and EHT reconstructed ensembles → final numerical/graphical results.

## Additional script-only checks and diagnostics

```bash
pixi run test
pixi run diagnostics --run-root /path/to/run-name \
  --output-dir /path/to/run-name/diagnostics
pixi run convergence --run-root /path/to/run-name \
  --output /path/to/convergence.json --modes std cor cor_grmhd
pixi run transport --help
```

The new observation product preserves the EHT summary's saved circular lift
(`center` and `unwrapped_mean`); it never wraps a broad ensemble again around its
arithmetic mean. Selected finite template pixels are validated after masking,
so excluded undefined phases cannot invalidate useful pixels. Stage replacements
reject input files inside the output directory, including resolved symlink aliases.
Inference checks upstream lineage and effective configuration, not just array readability.

## Constant/provenance references

- [M87* constant and calibration audit](M87_CONSTANTS_AUDIT.md)
- Original code and earlier classification/run-chain records: separate `EHT_ZC` worktree,
  not runtime files in this cleaned checkout.
- EHT-I https://arxiv.org/abs/1906.11238
- EHT-VII https://arxiv.org/abs/2105.01169
- EHT-VIII https://arxiv.org/abs/2105.01173

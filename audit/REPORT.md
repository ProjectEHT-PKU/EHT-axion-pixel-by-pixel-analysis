# Reproduction and bug-fix audit

Date: 2026-09-14/15

Branch: `audit/reproduce-fix-final-results`

## Scope and resource policy

The final six curves were rerun personally from validated intermediates in detached
`tmux` sessions.  The 22 GB RAPTOR/GRMHD raw stage was not rerun because its
products would take much longer and the machine has 31 GiB RAM and no swap.
Instead, the inputs immediately upstream of inference were validated, and the
EHT 4800→1572 topset selection plus corrected 1-D Stokes integration were rerun.
No historical `.npy`/`.npz` products were overwritten.

## What was reproduced

The historical result consists of three uncertainty modes for both 1-D and 2-D:

1. `std`: diagonal observation standard deviations;
2. `cor`: full observation covariance;
3. `cor_grmhd`: observation covariance + scalar GRMHD evolution variance +
   propagated axion-template variance.

`scripts/run_final_limits.py --legacy` reproduced every value in
`Results/data/{std,cor,cor_grmhd}{,_1d}.npz`.  Across all 48 values, the maximum
absolute difference was below `3.4e-15`; see `legacy_reproduction.json`.

## Intermediate validation

- Four saved 1572-image EHT Stokes stacks reproduce their saved I and P means:
  max I difference `0`, max P difference `5.5e-20`.
- The historical 2-D condition selects 770 pixels.
- GRMHD Gaussian products were independently refit at 256 sampled locations per
  dimension.  Maximum differences were:
  - 1-D amplitude `1.27e-11`, phase `9.57e-8` rad;
  - 2-D amplitude `1.72e-10`, phase `9.76e-9` rad.
- The historical 1572-image selection was inferred directly from saved images
  for all four days. Every day agrees 1572/1572 with the independently
  recomputed legacy ranking, with exact Stokes signatures; all four use the
  same hyperparameter index set.

## Verified bugs fixed

### EHT topset coordinate mismatch

The first ranking pass translated fitted means but not per-image EVPAs before
forming residuals.  Correct alignment changes 198 of 1572 selected images;
rank correlation is 0.9725 and maximum individual rank displacement is 1380.
The corrected second-stage topset has 757 valid pixels, 756 overlapping the
historical 770-pixel condition.

### Incorrect combination of independent epochs

The two epoch-pair likelihoods were added.  Independent data must multiply,
therefore log likelihoods now add.  In a controlled 2-D `cor_grmhd` ablation,
this changes `log10(c95)` by up to 0.177 dex.

### Missing parameter-dependent Gaussian normalization

`M_axion` scales with coupling, so `log(det(M))` is not a removable constant.
The corrected Gaussian log likelihood includes it.  Maximum ablation impact:
0.271 dex.

### Accidental posterior quadrature and hidden prior

The old equal-weight cumulative sum on a log grid only implicitly selected a
log-uniform prior and used rectangular point weights.  Prior choice is now
explicit and trapezoidal integration is performed in coupling.  Maximum impact:
0.060 dex.

### Duplicated nuisance-phase endpoint

`linspace(0, 2*pi, 16)` sampled 0 and 2π twice.  Corrected sampling uses 16 unique
phases.  Maximum impact: 0.0035 dex.

### Linear interpolation across a phase branch cut

Fitted phase means are now interpolated on the unit circle rather than linearly
across ±π.

### Wrong uncertainty phase at the second epoch pair

The second signal is advanced by five days, but its propagated model uncertainty
previously reused the first epoch phase.  It now uses its own phase.  Current
numerical impact is small (`1.1e-5` dex) but the former expression was wrong.

### 1-D Stokes integration

EHT pickle diagnostics show Q/U have the same flux scale as I and P/I has median
0.13–0.15: Q/U are Stokes fluxes, not fractional q/u.  The old integrator
multiplied Q and U by I a second time.  It also called a three-degree sum a
10-degree aperture.  The corrected path integrates Q/U directly over an actual
10-degree circular window.

A cached sparse linear operator reproduces the corrected spline implementation
to `7.7e-14` degree maximum error.  It rebuilds all four EHT 1-D stacks and all
eight GRMHD alpha stacks in 11.6 s.  Correcting both the topset and Stokes
integration changes median EHT 1-D EVPAs by 9.5–21.6 degrees depending on day.

### Deterministic code defects

- fixed scalar and non-symmetric interval behavior in `CircleCalc`;
- `ImageS(blur=False)` now honors `rotation_angle`;
- stopped `ImageS` from corrupting its own `data` attribute;
- removed incomplete `plt.imshow(` syntax from `result.ipynb`;
- expensive Notebook regeneration is opt-in rather than immediately reset true;
- removed a deliberate `KeyboardInterrupt`;
- renamed two colliding `Fit_EVPA_different_omega` functions;
- comparator now returns zero for equal EHT metadata;
- original EHT Notebook recomputes per-image EVPA after translation.

## Corrected expected/null limits

The historical final code explicitly multiplied measured changes by zero.  The
main figure therefore retains that interpretation as **expected/null limits**.
Values below are `log10(c95)`.

| alpha | std 1-D | std 2-D | cor 1-D | cor 2-D | full 1-D | full 2-D |
|---:|---:|---:|---:|---:|---:|---:|
|0.10|1.952|1.329|1.969|1.785|1.988|1.630|
|0.15|1.206|0.576|1.219|1.036|1.238|0.898|
|0.20|0.709|0.079|0.720|0.532|0.738|0.435|
|0.25|0.344|-0.295|0.357|0.154|0.368|0.067|
|0.30|0.046|-0.604|0.060|-0.163|0.058|-0.309|
|0.35|-0.190|-0.867|-0.179|-0.428|-0.196|-0.633|
|0.40|-0.387|-1.089|-0.377|-0.649|-0.403|-0.886|
|0.44|-0.497|-1.207|-0.491|-0.773|-0.515|-1.022|

The single six-curve deliverable is:

- `figures/final_limits_corrected.pdf`
- `figures/final_limits_corrected.png`

`figures/legacy_vs_corrected.*` shows before/after differences.  Machine-readable
values are in `figures/final_limits.csv` and `figures/final_limits.json`.

## Observed-path validation

The residual is now explicit: default null data, or `--observed` for measured
April 5→6 and April 10→11 differences.  All six observed configurations were
also executed successfully and saved under
`audit/results/rebuilt_corrected_observed/`.  They are not mixed into the null
figure because that would change the scientific interpretation of the original
requested result.

## Performance improvements

- Build only selected 770/757 covariance blocks instead of four 4096² matrices:
  roughly 28× fewer covariance entries and ~500 MiB avoided.
- Cache Cholesky factors for covariance independent of alpha/coupling/phase.
- Work in log-likelihood space (`logsumexp`) to avoid underflow.
- Replace per-image spline construction with a cached 4.4%-dense Stokes operator.
- Replace histogram `curve_fit` where appropriate with unbinned Gaussian MLE;
  corrected 1-D fits fell from 6.3 s to 0.20 s and no longer reject finite bins.
- The corrected six inference runs complete in about 63 s total on this machine;
  corrected EHT topset reconstruction took 43.7 s.
- Doubling the coupling grid from 60 to 120 points changes the full-model limits
  by at most 0.0094 dex (`grid_convergence.json`).

## Remaining scientific choices (not silently changed)

1. Cross-day EHT reconstruction covariance is assumed zero.  Shared hyperparameter
   systematics may require `C22 + C11 - C12 - C21` rather than `C22 + C11`.
2. GRMHD evolution uncertainty is represented by one scalar diagonal variance,
   not a spatial covariance.
3. The main plot is a null/expected limit.  An observed limit is materially
   weaker and must be labeled separately.
4. The coupling-grid bounds and log-uniform prior remain analysis choices and
   should be stated in any paper caption/method section.
5. Raw data distribution/checksums still need a durable archive (Zenodo/OSF);
   generated multi-GB tables must not be committed to ordinary Git.

## Verification commands

```bash
export PYTHONPATH=$PWD
pytest -q
ruff check eht_axion scripts tests
python -m compileall -q eht_axion scripts tests utils
python scripts/verify_legacy_results.py
python scripts/plot_final_limits.py
```

At audit completion: 24 tests pass; new/changed Python code passes Ruff and
compile checks; all three main Notebooks pass per-cell Python syntax parsing.

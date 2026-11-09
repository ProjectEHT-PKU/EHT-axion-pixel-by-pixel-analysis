# EHT imaging

## Purpose

This page describes the EHT imaging input: a reproduction of the EHT 2017 M87\*
polarimetric imaging method of EHT Paper VII (with the total-intensity workflow of Paper IV)
using ehtim 1.3.2, followed by a Monte Carlo over the instrumental leakage (D-term)
solution that gives 1000 polarimetric images per day. The per-day ensembles of these
images are the EHT input of the [Analysis](analysis.md). EHT imaging does not read any
RAPTOR product.

## Inputs

- The four HOPS low-band uvfits files `SR2_M87_2017_<day>_lo_hops_zbl-dtcal+selfcal.uvfits`
  for days 95, 96, 100 and 101 (2017 April 5, 6, 10, 11), in the directory
  `EHT_IMAGING_INPUT_ROOT`. Their names and SHA-256 are `eht_imaging.inputs`.
- An input receipt (`EHT_IMAGING_INPUT_RECEIPT`): JSON with `complete: true` and a `files`
  list (`path`, `sha256`, `size`) that must agree with `eht_imaging.inputs`.
- The feed-rotation table (`EHT_IMAGING_FEED_ROTATION`): JSON with `status: verified` and,
  per station, `fr_par`, `fr_elev`, `fr_off` (degrees).
- The official EHT imaging script (`EHT_IMAGING_OFFICIAL_SCRIPT`).
- `experiment.yaml`: `shared.image` and the `eht_imaging` section
  ([Experiment configuration](experiment-config.md)).

The four locations are machine paths, given on the command line or in `.env`.

## Outputs

- I stage (one per day): fixed total-intensity image `final-fixedI.fits` and gain-calibrated
  native circular data `clean-gaincal-pol.uvfits`.
- Fiducial stage (one per day): P images and D terms of every round, final D terms.
- Refined stage (one per day): exact P images, refined D terms (`final-dterms.json`),
  certificates.
- Bundle: the four I products, the four refined D-term files and the feed table, with
  `bundle.json`.
- Monte Carlo batch: one certified P image per draw (4 × 1000), with receipts.
- Ensemble summary: per day `stokes.npy`, `evpa.npy`, `aperture.npz`, `global.npz`, plus
  `axes.npz` and `report.json`.

File details are on [Data layout](data-layout.md).

## How to run

Environment `eht-rebuild` (ehtim 1.3.2, astropy); `eht_summarize.py` needs only
`eht-research`. Every command writes one new root.

```bash
# one root per day and stage; repeat for 95, 96, 100, 101
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --plan
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-rebuild eht-reconstruct --stage fiducial --day 95 --i-root I95 \
    --output-root PD95 --review-file REVIEW
pixi run -e eht-rebuild eht-refine --day 95 --source-root PD95 --i-root I95 --output-root D95

# after the four days
pixi run -e eht-rebuild eht-export-bundle --i-root I95 --i-root I96 --i-root I100 --i-root I101 \
    --d-root D95 --d-root D96 --d-root D100 --d-root D101 --output-root BUNDLE
pixi run -e eht-rebuild eht-sample --bundle-root BUNDLE --output-root MC --workers 64 --threads 2
pixi run -e eht-research eht-summarize --dataset-root DATASET --output-root ENSEMBLE
```

`eht_reconstruct.py --plan` only verifies and hashes the inputs and prints the plan.
Execution requires `--review-file`, a nonempty file recording an external review decision;
its bytes become a stage input, and its presence is not an automatic approval.

`eht_sample.py` options: `--workers` (default 64), `--threads` per worker (default 2),
`--max-active-memory-gib` (default 100), `--resource-wait` (default 300 s), `--task-limit`,
`--resume`. A run with the same output root resumes; it requires an identical batch
identity (bundle, science values and code).

`eht_summarize.py` reads `DATASET/COPY_VERIFIED.json` and the Monte Carlo results in
`DATASET/results/full4000-v1/` (see [Ensemble summary](#ensemble-summary)).

## Details

### Data loading

`load_native` (`eht_axion/eht/calibration.py`) loads a HOPS file in native circular
polarization (RR, LL, RL, LR) without conversion to Stokes. A whole row is dropped if any
of the four products or errors is missing or nonpositive. The feed-rotation coefficients of
the table are written into the station table, D terms are set to zero, and the
calibration flags are `ampcal = phasecal = frcal = True`, `dcal = False`: the stations
AA, AP, SM, JC are already leakage-corrected in the HOPS products, and AZ, LM, PV are not.
Data are scan-averaged coherently once.

### I stage

`eht_reconstruct.py --stage I`:

- runs the official EHT imaging script on the raw file with exactly two source
  substitutions: the Fourier transform type `nfft` becomes `direct`, and the field of view
  128 µas becomes `shared.image.fov_uas` (120 µas). The final image `im_out` of the script
  (after its own restoring blur) is the fixed I image; it must already be positive, finite
  and on the 64 × 64, 120 µas grid;
- solves **shared station gains**: on the native scan-averaged data with baselines shorter
  than `uv_min` removed, an ehtim self-calibration of all stations (including PV) against
  the fixed I image with `method='both'`, `gain_tol=[0.02, 0.2]`, and equal right and left
  gains;
- applies these gains to the native circular data and writes
  `clean-gaincal-pol.uvfits`, `shared-gains.npz` and `polarization-input.json`.

The total-intensity image is fixed from here on; all later steps solve only for linear
polarization $P = Q + iU$ (with $V = 0$).

### Fiducial P/D alternation

`eht_reconstruct.py --stage fiducial` reads the I stage of the same day. The base data are
the gain-calibrated data with `uv_min` applied. Starting from $P = m_0\, I\, e^{2i\chi_0}$
(`initial_m`, `initial_evpa_deg`), each of `rounds` (50) rounds:

- restarts P at the rounds in `restart_rounds` from the I-weighted fractional
  polarization blurred with a `restart_blur_uas` (20 µas) Gaussian, scaled by
  $1-\epsilon$ (`restart_domain_epsilon`) to stay inside $\lvert P/I\rvert < 1$;
- adds the systematic error `fsys` to the RL/LR errors (`add_fsys`);
- applies the current D terms by one inverse Jones correction (`apply_fixed_d`);
- images P with the ehtim `Imager` (`pol='P'`, transform `mcv`, data terms `p_data`:
  `pvis` and `m`; regularizers `p_reg`: `hw` entropy and `ptv` total variation;
  `norm_reg`; `maxit`; `stop`; direct Fourier transform; flux of the fixed I image);
  with `shift_objective_constant` a constant is added to the objective, which does not
  change the minimizer;
- refits the D terms of `solve_sites` (AZ, LM, PV) against the same uncorrected noisy base
  and the new P image with the official forward RL/LR leakage fit of ehtim
  (`fit_forward_d`), in which the D prior and its gradient are replaced by zero.

A final P-only solve aligns the delivered image with the last D terms. Each round's D
terms and raw and display images are written; `termination.json` records that the stop is
the fixed round budget.

### Refinement of the fiducial D terms

`eht_refine.py` refines the fiducial D terms of one day (the "tightD" stage) with an exact
convex P solver (`eht_axion/eht/convex.py`). For fixed I and fixed D, the Paper VII P
objective with `pvis` = `m` = `ptv` = 1, `hw` = 100 and normalized regularizers is convex
in $P$:

$$
J(P) = \sum_{t\in\{\mathrm{pvis},m\}} \frac{1}{2N_t}\sum_j \frac{\lvert (A_t P)_j - d_{t,j}\rvert^2}{\sigma_{t,j}^2}
 + \frac{100}{F}\sum_x I_x\, h\!\left(\frac{\lvert P_x\rvert}{I_x}\right)
 + \frac{\mathrm{TV}(P)}{\tau_\mathrm{TV}},
$$

with $h(m) = \tfrac12[(1+m)\ln(1+m) + (1-m)\ln(1-m)]$, $F$ the total flux and
$\tau_\mathrm{TV}$ the ehtim TV normalization. The data tuples are those of the official
ehtim imager. `solve_problem` runs a diagonally preconditioned primal-dual hybrid gradient
(PDHG) method and computes, every `check_every` iterations, a feasible dual lower bound;
the returned P is certified when the best primal value minus the best lower bound is at
most `gap_atol`. This is a floating-point certificate, not interval arithmetic.

Starting from the fiducial final P and D, each of `refine.outer_rounds` (100) rounds solves
P to the gap tolerance (`gap_atol`, `maxiter`, `check_every`, `time_limit` of
`eht_imaging.refine`) and refits D with `fit_forward_d` using the L-BFGS-B options
`eht_imaging.refine.d_optimizer`; every optimizer call is recorded in
`D-optimizer-records.json`. A final P solve uses the last D. `termination.json` records
`all_P_numerically_certified`; joint P/D optimality is not claimed.

### Bundle export

`eht_export_bundle.py` takes the four I roots and the four refined roots, verifies each
(full lineage), checks that their imaging parameters equal `eht_imaging_parameters` of the
current `experiment.yaml`, that each refined stage descends from the I stage of the same
day, that every P solve of each refined stage is certified, and that all four I stages used
one feed table. It copies, read-only:

- per day `fixedI.fits`, `gaincal-pol.uvfits` and `fiducial-dterms.json`;
- lineage files (`lineage/<day>/I-manifest.json`, `D-manifest.json`, `D-plan.json`,
  `D-termination.json`, `D-steps.json`);
- `feed-rotation.json` and the script itself.

`bundle.json` lists every file with SHA-256 and size and records the science contract:
`fsys`, `covariance_ddof`, `residual_component_sigma` (`residual_d_sigma`), the sky grid and
the P objective weights, and `source_kind = localrebuiltfourdayfiducial`. The equivalence
of this script's output with the bundle recorded for the existing Monte Carlo run has not
yet been checked.

### D-term Monte Carlo

`eht_sample.py` verifies the bundle (every file hash, the science contract against
`experiment.yaml`) and then:

- forms the D-term distribution (`d_distribution`): for each station in `solve_sites` and
  each hand R, L, the mean and the 2 × 2 covariance of $(\mathrm{Re}\,D, \mathrm{Im}\,D)$ over
  the four refined daily solutions, with `covariance_ddof` = 1;
- draws `uncertainty_draws` = 1000 D sets per day (`sample_dterms`) with seed
  `SeedSequence([seed, day])`: bivariate normal for AZ, LM, PV; zero-mean normal with
  standard deviation `residual_d_sigma` = 0.01 per real and imaginary component for the
  pre-corrected stations AA, AP, SM, JC. Thermal noise is not resampled;
- for each draw: adds `fsys`, applies the drawn D terms once, and solves the exact convex P
  problem from $P = 0.05\,I$ (`initial_p_fraction`) with `eht_imaging.monte_carlo`
  (`gap_atol`, `maxiter` 50000, `check_every`, `time_limit` 600 s). An uncertified draw gets
  one second attempt with `retry_maxiter` 200000 and `retry_time_limit` 1800 s; a task with
  two uncertified attempts or two errors is reported as exhausted;
- writes, per certified draw, the raw and the display image (19 µas blur,
  `best-display-IQU-additional19uas.fits`) and a certificate.

The batch identity is the SHA-256 of the plan (bundle snapshot, science values, source and
package file hashes, execution threads); `finish.json` is written only when all 4000 draws
are certified.

### Ensemble summary

`eht_summarize.py` consumes a copy of the Monte Carlo output: it expects
`DATASET/COPY_VERIFIED.json` (schema `eht-copyback-verified-v1`, with the batch identity,
the SHA-256 of `batch.json` and `finish.json`, and flags confirming that the copied files
were hash-checked) and the batch in `DATASET/results/full4000-v1/`. It checks the receipts,
identities, certificates and the SHA-256 of every consumed file, then reads the 1000
display19 images per day (no further blur) and writes, per day under
`native/display19/day-<day>/`:

- `stokes.npy`: shape (1000, 64, 64, 4), axes draw, y, x, IQUV in the native FITS
  orientation (east offset decreases along x, north increases along y);
- `evpa.npy`: EVPA $\tfrac12\arg(Q+iU)$ modulo $\pi$ in radians, NaN where $I \le 0$ or
  $Q+iU=0$;
- `aperture.npz`: native sector curves at 360 one-degree angles east of north, each sector
  of full width `analysis.aperture.window_degrees` and radius up to 60 µas: `intensity`
  $=\sum I$ and `polarized_brightness` $=\sum I\lvert P\rvert/\sum I$. The analysis uses
  these only for 1D feature selection;
- `global.npz`: whole-field diagnostics per draw: total $I$,
  $m_\mathrm{net} = \lvert\sum P\rvert/\sum I$,
  $m_\mathrm{res} = \sum\lvert P\rvert/\sum I$, EVPA of $\sum P$, and
  $\beta_2 = \sum P\,e^{-2i\varphi}/\sum I$ with $\varphi$ the azimuth east of north about
  the grid centre (no annulus cut). These are for comparison with Paper VII and are not
  used by the analysis.

`axes.npz` holds `east_uas`, `north_uas` and `angle_deg`; `report.json` records the
conventions.

### Image centre

EHT Paper VII Sec. 5.2 aligns images by cross-correlation. In this reproduction that
alignment is a zero shift: the images are used as reconstructed on the common grid, and
the analysis has no re-centring option.

### Known deviations from EHT Paper VII

- This is a method reproduction, not a bitwise reproduction of the published images. The
  official I script runs with a direct Fourier transform and a 120 µas field of view.
- The locally solved D terms differ from Paper VII Table 5. The largest difference noted
  during validation is for PV on 2017 April 10 (day 100), where $D_R$ differs by about 19
  percentage points. The code does not compare with Table 5; this is a recorded validation
  result.
- `restart_rounds`, `maxit` and the refinement budget are numerical choices, not values
  from the paper.

## Related pages

- [Documentation index](index.md)
- [Overview](overview.md)
- [Experiment configuration](experiment-config.md)
- [Analysis](analysis.md)
- [Provenance](provenance.md)
- [Data layout](data-layout.md)

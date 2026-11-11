# EHT Imaging

[← Documentation index](index.md)

The EHT imaging pipeline reproduces the EHT 2017 M87\* polarimetric reconstruction method of EHT Paper VII (with total-intensity imaging following Paper IV) on `ehtim == 1.3.2`, refines the daily instrumental leakage (D-term) solutions with an exact convex polarimetric solver, and runs a 1000-draw-per-day D-term Monte Carlo whose summarized ensembles feed the [Analysis](analysis.md). It does not read any `RAPTOR` simulation product.

## Prerequisites and inputs

`scripts/eht/reconstruct.py`, `refine.py`, `export_bundle.py`, and `sample.py` run in the `eht-imaging` Pixi environment (`ehtim 1.3.2`, `astropy`); `scripts/eht/summarize.py` requires only `eht-summary` (`astropy`). `scripts/eht/reconstruct.py` reads three machine-local locations from the `eht/` directory of the read-only input root (`EHT_INPUT_ROOT` / `--input-root`, see the [input layout](data-layout.md#input-layout)); each can be overridden by an explicit option, and relative paths are resolved against the working directory:

- **HOPS `uvfits` directory** (`--uvfits-dir`, default `$EHT_INPUT_ROOT/eht`): contains the four released low-band files `SR2_M87_2017_<day>_lo_hops_zbl-dtcal+selfcal.uvfits` for days `95`, `96`, `100`, `101` (2017 April 5, 6, 10, 11), whose filenames and SHA-256 digests must match `eht_imaging.inputs` in `configs/experiment.yaml`.
- **Feed-rotation table** (`--feed-rotation-file`, default `$EHT_INPUT_ROOT/eht/feed-rotation.json`): JSON with `"status": "verified"`, `"units": {"fr_off": "...degrees..."}`, and per-station `fr_par`, `fr_elev`, `fr_off`.
- **Official EHT imaging script** (`--official-script`, default `$EHT_INPUT_ROOT/eht/eht-imaging_pipeline.py`): official total-intensity imaging pipeline script.

The SHA-256 of the selected day's `uvfits` file is checked against `eht_imaging.inputs` before any imaging; the `uvfits` file, feed-rotation table, and official script are hashed into the stage inputs, and their paths are recorded as the `input_root`, `feed_rotation_file`, and `official_script` imaging parameters (`EHT_LOCATIONS`).

## Running the EHT pipeline

With `EHT_INPUT_ROOT` and `EHT_DATA_ROOT` set in `.env`, run the per-day reconstruction and refinement for all four observing days (`95`, `96`, `100`, `101`), then export the bundle, draw the Monte Carlo ensemble, and summarize:

```bash
# Optional dry-run plan verification (no ehtim import or disk writes)
pixi run -e eht-imaging eht-reconstruct --stage I --day 95 --plan

# 1–3. Per-day fixed-I, fiducial P/D, and convex D refinement
for day in 95 96 100 101; do
  pixi run -e eht-imaging eht-reconstruct --stage I --day $day
  pixi run -e eht-imaging eht-reconstruct --stage fiducial --day $day
  pixi run -e eht-imaging eht-refine --day $day
done

# 4. Bundle the four fixed-I and refined-D stages for Monte Carlo sampling
pixi run -e eht-imaging eht-export-bundle

# 5. Run the 4 x 1000 D-term Monte Carlo (resumable in $EHT_DATA_ROOT/eht/sample)
pixi run -e eht-imaging eht-sample --workers 64 --threads 2

# 6. Summarize certified blurred FITS draws into $EHT_DATA_ROOT/eht/ensemble
pixi run -e eht-summary eht-summarize
```

Every output directory except `eht/sample` must be new. `scripts/eht/sample.py` supports `--workers` (default `64`), `--threads` (default `2`), `--max-active-memory-gib` (default `100`), `--resource-wait` (default `300` s; exits with code `3` if memory admission remains blocked, safe to resume), `--task-limit`, and `--resume`. Re-running `eht-sample` on the same output directory resumes automatically if the batch plan SHA-256 (`batch.json` `identity`) is identical.

## Reconstruction method

### Native circular data loading (`eht_axion/eht/calibration.py`)

`load_native` loads each HOPS `uvfits` file in native circular polarization (`polrep="circ"`, products `RR`, `LL`, `RL`, `LR`) without round-tripping through Stokes:
- Any visibility record where any of the four correlations or sigmas is non-finite or has $\sigma \le 0$ is dropped completely before averaging.
- Feed-rotation angles (`fr_par`, `fr_elev`, `fr_off`) are populated from `feed-rotation.json`, station D terms are initialized to `0j`, and calibration flags are set to `ampcal = phasecal = frcal = True`, `dcal = False` (ALMA/APEX/SMA/JCMT in `fixed_sites` are pre-corrected in HOPS; `AZ`, `LM`, `PV` in `solve_sites` remain uncorrected).
- Visibilities are coherently scan-averaged once (`avg_coherent(0, scan_avg=True)`).
- `make_imaging_obs` wraps `Obsdata` to preserve structured dtypes on rectangular `tlist` groups under NumPy 2 and fixes the missing factor of $1/4$ in $\operatorname{Var}(I) = (\sigma_{RR}^2 + \sigma_{LL}^2)/4$ in `ehtim` 1.3.2's `msigma` error propagation.

### Stage `I`: fixed total intensity and shared gains

`scripts/eht/reconstruct.py --stage I` publishes stage `eht:i`:
1. Copies the official EHT script and executes an adapted copy with **only two exact string substitutions**: `"ttype     = 'nfft'"` $\to$ `"ttype     = 'direct'"` and `'fov       = 128*eh.RADPERUAS'` $\to$ `'fov       = 120*eh.RADPERUAS'`. All other settings — including the compact zero-baseline flux `zbl = 0.60` Jy and the $5\,\mu\mathrm{as}$ restoring blur producing `im_out` — remain unchanged. `im_out` is verified to be positive, finite, and on the $64 \times 64$, $120\,\mu\mathrm{as}$ grid, and saved as `final-fixedI.fits`.
2. Excludes intra-site baselines shorter than `uv_min` ($10^8\,\lambda = 100\,\mathrm{M}\lambda$) from the scan-averaged native circular data, then solves **shared $R = L$ station gains** across all present stations (including `PV`) against `final-fixedI.fits` via `eh.selfcal(..., pol='I', method='both', ttype='direct', gain_tol=[0.02, 0.2])`.
3. Applies the shared gains (`applycal`) to the native circular visibilities and writes `clean-gaincal-pol.uvfits`, `shared-gains.npz`, and `polarization-input.json`.

All subsequent stages hold $I$ fixed to `final-fixedI.fits` and solve exclusively for linear polarization $P = Q + iU$ ($V = 0$).

### Stage `fiducial`: alternating $P$ imaging and forward $D$ fits

`scripts/eht/reconstruct.py --stage fiducial` publishes stage `eht:fiducial`. Starting from uniform fractional polarization $P = m_0 I e^{2i\chi_0}$ (`initial_m = 0.1`, `initial_evpa_deg = 0.0`) and `dterms = {}` on the `uv_min`-cut gain-calibrated data, each round $r = 1, \dots, 50$ (`rounds`):
1. At $r \in \{2, 12, 22, 32, 42\}$ (`restart_rounds`), re-initializes $P$ from the $I$-weighted blurred complex polarization fraction,
   $$P_{\mathrm{init}} = (1 - \epsilon)\,I\,\frac{K_{20} * P}{K_{20} * I},$$
   with $20\,\mu\mathrm{as}$ Gaussian `restart_blur_uas` and $\epsilon = 10^{-8}$ (`restart_domain_epsilon`) to remain strictly inside $\lvert P/I\rvert < 1$.
2. Adds fractional systematic noise `fsys = 0.005` in quadrature to `RL` and `LR` variances (`add_fsys`: $\sigma_{RL,LR}^2 \to \sigma_{RL,LR}^2 + \bigl(f_{\mathrm{sys}}\lvert(RR + LL)/2\rvert\bigr)^2$) and applies the current D terms once by inverse Jones calibration (`apply_fixed_d`, temporarily zeroing table SEFDs so input marginal sigmas are preserved).
3. Solves for $P$ with `ehtim` `Imager` (`pol='P'`, `transform=['mcv']`, `data_term={'pvis': 1.0, 'm': 1.0}`, `reg_term={'hw': 100.0, 'ptv': 1.0}`, `norm_reg=True`, `ttype='direct'`, `maxit=500`, `stop=1e-8`; `shift_objective_constant=True` adds the constant $w_{\mathrm{pvis}} + w_m + h_w\ln 2$ without changing the minimizer).
4. Refits D terms on `solve_sites` (`AZ`, `LM`, `PV`) against the **uncorrected** `fsys`-augmented base visibilities and the new $P$ image using `fit_forward_d`: a private two-line adaptation of `ehtim.calibrating.pol_cal_new.leakage_cal_new` (`pol_fit=['RL', 'LR']`, `ttype='direct'`) that sets the D-term Gaussian prior and its gradient identically to zero.

After round 50, one final `final-P-only` solve aligns `final-raw-IQU.fits` and `final-blurred-IQU.fits` ($19\,\mu\mathrm{as}$ blur) with `final-dterms.json`.

### Convex refinement (`scripts/eht/refine.py` and `eht_axion/eht/convex.py`)

For fixed $I > 0$ and fixed D terms, the normalized-regularizer `ehtim` polarimetric objective is strictly convex in $P \in \mathbb{C}^{64 \times 64}$ on the physical disk $\lvert P_x\rvert \le I_x$:

$$
J(P) = \sum_{t \in \{\mathrm{pvis}, m\}} \frac{w_t}{2 N_t} \sum_{j=1}^{N_t} \frac{\lvert (A_t P)_j - d_{t,j}\rvert^2}{\sigma_{t,j}^2} + \frac{h_w}{F} \sum_x I_x\,h\!\left(\frac{\lvert P_x\rvert}{I_x}\right) + p_{\mathrm{tv}}\,\frac{\operatorname{TV}(P)}{\tau_{\mathrm{TV}}},
$$

where $h(m) = \frac{1}{2}\bigl[(1+m)\ln(1+m) + (1-m)\ln(1-m)\bigr]$, $F = \sum_x I_x$, $\tau_{\mathrm{TV}}$ is the `ehtim` beam-scaled TV normalization, and $(d_t, \sigma_t, A_t)$ are extracted directly from `ehtim`'s `Imager._data_tuples`. $J(P)$ equals the `ehtim` objective plus $w_{\mathrm{pvis}} + w_m + h_w\ln 2$.

`solve_problem` minimizes $J(P)$ in reduced variables $x = P / I$ via diagonally preconditioned primal–dual hybrid gradient (PDHG, step-size safety factor $\theta = 0.99 < 1$, safeguarded Newton–bisection radial entropy prox) and evaluates a feasible Fenchel–Rockafellar dual lower bound every `check_every` (`500`) iterations. A solution is **numerically certified** when $\text{best\_primal} - \text{best\_lower} \le \text{gap\_atol}$ (`1e-6`; floating-point duality gap, not interval arithmetic).

`scripts/eht/refine.py` (stage `eht:refine`) starts from the fiducial `final-raw-IQU.fits` and `final-dterms.json`, runs `outer_rounds = 100` alternations of certified PDHG $P$ solves (`maxiter = 20000`, `time_limit = 300` s) and forward `L-BFGS-B` D refits (`eht_imaging.refine.d_optimizer`, logged to `D-optimizer-records.json`), followed by a 101st `final-P` solve aligned with the last D terms.

### Bundle export, D-term Monte Carlo, and ensemble summary

1. **Bundle export (`scripts/eht/export_bundle.py`)**: verifies full lineage of all four `eht:i` and four `eht:refine` stages, checks `all_P_numerically_certified == True` on every day, and copies read-only (`0444`) `<day>/fixedI.fits`, `<day>/gaincal-pol.uvfits`, `<day>/fiducial-dterms.json`, `lineage/<day>/*`, and `feed-rotation.json` with `bundle.json` (schema `eht-bundle`).
2. **D-term Monte Carlo (`scripts/eht/sample.py`)**:
   - Computes the 4-day mean and $2 \times 2$ $(\operatorname{Re} D, \operatorname{Im} D)$ sample covariance (`covariance_ddof = 1`) per hand (`R`, `L`) and station in `solve_sites` (`d_distribution`).
   - Draws `uncertainty_draws = 1000` D-term sets per day using `SeedSequence([20210405, day])`: bivariate normal for `AZ`, `LM`, `PV`, and independent zero-mean normal with per-component SD `residual_d_sigma = 0.01` for `AA`, `AP`, `SM`, `JC` (`sample_dterms`). Thermal visibility noise is not resampled and D terms are not refitted per draw.
   - Solves the exact convex $P$ problem cold-started from $P = 0.05\,I$ (`initial_p_fraction`) with `gap_atol = 1e-6`, `maxiter = 50000`, `time_limit = 600` s; uncertified draws receive one retry (`retry_maxiter = 200000`, `retry_time_limit = 1800` s). Writes `best-raw-IQU.fits`, `best-blurred-IQU.fits` ($19\,\mu\mathrm{as}$ blur), `certificate.json`, and `taskreceipt.json` per draw, and `finish.json` once all $4000$ draws are certified.
3. **Ensemble summary (`scripts/eht/summarize.py`)**: verifies `batch.json`, `finish.json`, all $4000$ `taskreceipt.json` and `certificate.json` files, and SHA-256 digests of every `best-blurred-IQU.fits`, then publishes stage `eht:ensemble` containing `axes.npz` (`east_uas`, `north_uas`), `report.json`, and per-day `day-<day>/{stokes.npy, evpa.npy, global.npz}`:
   - `stokes.npy`: float64 `(1000, 64, 64, 4)` (`[draw, y, x, IQUV]` in native FITS orientation);
   - `evpa.npy`: float64 `(1000, 64, 64)` EVPA in $[0, \pi)\,\mathrm{rad}$ (`NaN` where $I \le 0$ or $Q = U = 0$);
   - `global.npz`: Paper VII Sec. 5.3 whole-field diagnostics (`total_I`, `m_net`, `m_resolved`, `evpa_rad`, `beta2_real`, `beta2_imag`) computed over the $19\,\mu\mathrm{as}$-blurred images for comparison with Paper VII (used by `paper.py`, not by `limits.py`).

## Known deviations from EHT Paper VII

- **Method reproduction, not bitwise image reproduction**: the official total-intensity script is executed with two adaptations (`ttype = 'direct'` instead of `'nfft'`, and `fov = 120` $\mu\mathrm{as}$ instead of `128` $\mu\mathrm{as}$).
- **Local D-term solutions vs. Paper VII Table 5**: D terms solved locally by `reconstruct.py` and `refine.py` differ from the values tabulated in EHT Paper VII Table 5; the largest discrepancy is on 2017 April 10 (day `100`) at station `PV`, where $D_R$ differs by roughly 19 percentage points. The pipeline does not compare against or substitute Table 5 values.
- **Zero-shift image centring**: EHT Paper VII Sec. 5.2/5.3 aligns reconstructed images by maximizing the Stokes-$I$ cross-correlation with the Paper IV consensus image. In this reproduction that alignment is a zero shift: all $\beta_2$ diagnostics and all 1D/2D analysis features are evaluated about the native grid centre, with no re-centring option.
- **Numerical schedules and solver budgets**: `restart_rounds`, `restart_blur_uas`, `maxit`, and the PDHG / `L-BFGS-B` refinement and Monte Carlo budgets in `eht_imaging.refine` and `eht_imaging.monte_carlo` are numerical implementation choices rather than paper-mandated parameters.

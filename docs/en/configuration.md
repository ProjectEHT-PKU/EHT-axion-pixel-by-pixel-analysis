# Configuration

[← Documentation index](index.md)

Configuration is split strictly between scientific parameters in [`configs/experiment.yaml`](../../configs/experiment.yaml) and machine-local paths or execution settings in `.env` and CLI options (three path roots: `EHT_RAPTOR_SOURCE_ROOT`, `EHT_INPUT_ROOT`, `EHT_DATA_ROOT`).

## Parameter file structure and validation

`configs/experiment.yaml` is the single scientific parameter file of the repository (the default of every `--config` flag). Every script loads it through `read_experiment()` in `eht_axion/experiment.py`:

```bash
pixi run python -c "from eht_axion.experiment import read_experiment; read_experiment()"
```

An alternative YAML file passed via `--config PATH` must contain the exact same four top-level sections and nested key sets:

- `shared`: black-hole scales, $\alpha$ grid, bound-state quantum number, and the common $64 \times 64$ sky grid.
- `raptor`: `RAPTOR` source repository/commit, camera geometry, compile-time switches, plasma model, snapshot range, resolutions, radial-table assignments, and BHAC grid coordinates.
- `eht_imaging`: released HOPS `uvfits` identities, calibration tolerances, imaging weights, fiducial schedule, convex refinement budget, and D-term Monte Carlo settings.
- `analysis`: sky rotation, sector aperture geometry, feature selection thresholds, superradiant cloud scaling, GRMHD background lag, and Bayesian inference grid/prior settings.

**Validation rules (`eht_axion/experiment.py`):**
- Uses a custom `SafeLoader` (`_UniqueKeyLoader`) that rejects duplicate YAML keys (including merged keys).
- Enforces exact key sets in every mapping (`SCHEMA`), strict Python types (`bool` is rejected where `int` or `float` is expected; integers must be `int`), finiteness, and domain bounds; it fills zero defaults.
- Checks two cross-section geometric invariants: `shared.image.fov_uas` converted to $r_g$ ($31.42\,r_g$) must not exceed `raptor.camera_fov_rg` ($46.2261\,r_g$), and `analysis.features.min_radius_rg` ($2.0\,r_g$) must be strictly smaller than the aperture radius (`radial_max_fraction * final_fov / 2`, $15.71\,r_g$).

### Derived views and module constants

`eht_axion/experiment.py` exposes four derived views recorded in stage manifests and plans, plus two fixed constants:

| Symbol | Result | Consumers |
|---|---|---|
| `physical_scale(config)` | `PhysicalScale` (`tg_seconds` $= 32\,015.7\,\mathrm{s}$, `tg_hours` $= 8.893\,\mathrm{h}$, `angular_rg_uas` $= 3.819\,\mu\mathrm{as}$, unit conversions) | Analysis scripts, `RAPTOR` frame cadence (`input_cadence_hours`) |
| `ray_tracing_config(config)` | `shared` physics scalars (`mass_solar`, `distance_mpc`, `label_step_tg`, `alphas`, `principal_n`) plus `raptor` without `raptor.source` | `RAPTOR` task plans, `generated.json` delivery identity |
| `analysis_geometry(config)` | `image_preprocess` (`sigma`, `init_fov`, `final_fov`, `final_shape`, `rotation_angle` in $r_g$ and $\mathrm{rad}$), `aperturefit` (`fov`, `radial_max`, `window_degrees`, `angular_bins`), and `metadata` | `scripts/analysis/templates.py`, `scripts/analysis/limits.py` |
| `eht_imaging_parameters(config, locations)` | Flat dictionary combining `EHT_LOCATIONS`, `shared.image` (`npix`, `fov_uas`, `blur_fwhm_uas`), and `eht_imaging` excluding `refine` and `monte_carlo` | `scripts/eht/reconstruct.py`, `scripts/eht/refine.py`, `scripts/eht/export_bundle.py` |
| `EHT_DAYS = (95, 96, 100, 101)` | Observing days of year 2017 (`eht_imaging.inputs` must match this exact set) | All EHT imaging, summary, and analysis modules |
| `ANGULAR_BINS = 360` | Number of $1^\circ$ sector bins in the 1D aperture observable | `analysis_geometry`, `scripts/analysis/limits.py` |

## Parameter reference: `configs/experiment.yaml`

### `shared`

| Key | Value | Unit | Meaning and validation | Read by |
|---|---|---|---|---|
| `mass_solar` | `6500000000.0` | $M_\odot$ | M87\* mass ($> 0$) | `physical_scale`; `RAPTOR` `MBH` |
| `distance_mpc` | `16.8` | $\mathrm{Mpc}$ | Distance to M87\* ($> 0$; multiplied by $1000$ for `RAPTOR` `DISTANCE` in $\mathrm{kpc}$) | `physical_scale`; `RAPTOR` `DISTANCE` |
| `label_step_tg` | `0.1` | $t_g$ / label | Expected GRMHD footer time step per snapshot label unit ($> 0$) | `plan_generation` footer cadence check; `PhysicalScale` |
| `alphas` | `[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]` | – | Strictly increasing positive two-decimal $\alpha$ grid (8 values) | `RAPTOR` planning/delivery, templates, likelihood, paper |
| `principal_n` | `2` | – | Principal quantum number $n \ge 1$ in $\Omega(\alpha) = \alpha(1 - \alpha^2/(2n^2))$ | `RAPTOR` `AXION_OMEGA`; `analysis/model.py` |
| `image.shape` | `64` | pixels | Square sky-grid side length ($\ge 2$) | EHT `npix`, ensemble summary, analysis resampling |
| `image.fov_uas` | `120.0` | $\mu\mathrm{as}$ | Sky-grid field of view ($1.875\,\mu\mathrm{as}$/pixel; replaces `128` $\mu\mathrm{as}$ in the official EHT script) | EHT `fov_uas`, analysis resampling and aperture |
| `image.blur_fwhm_uas` | `19.0` | $\mu\mathrm{as}$ | Circular Gaussian matching kernel FWHM ($\ge 0$) applied to EHT blurred images and to simulated images (not an effective PSF) | EHT `save_image`; analysis preprocessing $\sigma$ |

### `raptor`

| Key | Value | Unit | Meaning and validation | Read by |
|---|---|---|---|---|
| `source.repository` | `github.com/projecteht-pku/raptor` | – | Normalized lowercase `host/owner/name`; at least one git remote of the local checkout must match | `build_raptor`, `generate.py` plan `raptor_source`, `run_jobs` |
| `source.commit` | `20056006820b6536c96de790dd5fdaad5ada4eae` | – | Exact 40-hex commit required at `HEAD` of a clean local `RAPTOR` checkout | `build_raptor`, `run_jobs`, `deliver_stokes`, `paper.py` (`src/axion.c`) |
| `camera_fov_rg` | `46.2261` | $r_g$ | Square camera width ($176.5\,\mu\mathrm{as}$) | `RAPTOR` `CAM_SIZE_X/Y`; `analysis_geometry` `init_fov` |
| `fluid_spin` | `0.9375` | – | Dimensionless black-hole spin $a \in [0, 1]$ checked against every BHAC snapshot footer | `read_snapshot_footer`, `paper.py` |
| `build.model` | `bhac` | – | `setup.sh --code` (`bhac`, `harm3d`, `hamr`; `run_jobs` requires `MKSBHAC`) | `build_raptor`, `task_identity` |
| `build.metric` | `mks` | – | `setup.sh --metric` (`mks`, `cks`) | `build_raptor`, `task_identity` |
| `build.integrator` | `rk4` | – | `setup.sh --int` (`rk4`, `rk2`, `rk45`, `ver`) | `build_raptor`, `task_identity` |
| `build.grid` | `SMR` | – | `setup.sh --grid` (`SMR`, `AMR`) | `build_raptor`, `task_identity` |
| `build.polarized` | `true` | – | `setup.sh --rad=pol` (`true` required by `run_jobs`) | `build_raptor`, `task_identity` |
| `build.sfc` | `true` | – | `setup.sh --sf=sfc` | `build_raptor`, `task_identity` |
| `build.magnetic_polarity` | `plus` | – | `setup.sh --bflip` (`plus`, `minus`) | `build_raptor`, `task_identity` |
| `build.sigma_cut` | `1.0` | – | Magnetization cut `setup.sh --sigma-cut` ($> 0$) | `build_raptor`, `task_identity` |
| `build.thetae_min` | `0.001` | – | Minimum dimensionless electron temperature `--thetae-min` ($> 0$, $< \text{thetae\_max}$) | `build_raptor`, `task_identity` |
| `build.thetae_max` | `100.0` | – | Maximum dimensionless electron temperature `--thetae-max` | `build_raptor`, `task_identity` |
| `build.uniform_output` | `false` | – | `setup.sh --uniform=0` (`false` required by `run_jobs`) | `build_raptor`, `task_identity` |
| `plasma.mass_unit_g` | `3.0e+25` | $\mathrm{g}$ | Accretion-flow mass scale `M_UNIT` ($> 0$) | `RAPTOR` `model.in` |
| `plasma.r_low` | `1.0` | – | Ion-to-electron temperature ratio parameter `R_LOW` ($> 0$) | `RAPTOR` `model.in` |
| `plasma.r_high` | `1.0` | – | Ion-to-electron temperature ratio parameter `R_HIGH` ($> 0$) | `RAPTOR` `model.in` |
| `plasma.inclination_degrees` | `163.0` | $\mathrm{deg}$ | Observer inclination `INCLINATION` ($\in [0, 180]$) | `RAPTOR` `model.in`, spectrum file name |
| `plasma.frequency_hz` | `230.0e+9` | $\mathrm{Hz}$ | Single observing frequency `FREQ_MIN` ($> 0$) | `RAPTOR` `model.in`, HDF5 dataset names |
| `plasma.step_size` | `0.005` | – | Geodesic integrator step `STEPSIZE` ($> 0$) | `RAPTOR` `model.in` |
| `snapshots.{start, stop, step}` | `1600, 3000, 10` | label | Inclusive snapshot label grid ($141$ snapshots, `0 <= start <= stop <= 999999`) | `plan_generation`, `deliver_stokes`, `generated_inputs` |
| `response.resolution` | `100` | pixels | Square side of axion-on and 100-px baseline images (multiple of 10) | `plan_generation`, `deliver_stokes`, `templates.py` |
| `response.reference_norm` | `0.1` | – | `AXION_NORM` reference amplitude ($> 0$) | `RAPTOR` `model.in`, `limits.py`, `paper.py` |
| `background.resolution` | `400` | pixels | Square side of axion-off background images (multiple of 10) | `plan_generation`, `deliver_stokes`, `templates.py`, `paper.py` |
| `radial_tables.low.spin` | `0.9375` | – | Declared spin $a \in [0, 1]$ of the `low` radial CSV | `plan_generation`, `task_identity`, `paper.py` |
| `radial_tables.high.spin` | `0.99` | – | Declared spin $a \in [0, 1]$ of the `high` radial CSV (intentionally distinct from `fluid_spin`) | `plan_generation`, `task_identity`, `paper.py` |
| `alpha_table` | `0.10`–`0.30`: `low`; `0.35`–`0.44`: `high` | – | Maps every two-decimal string key `'0.10'`, …, `'0.44'` to `low` or `high` | `plan_generation`, `run_jobs`, `paper.py` |
| `grid.nxlone` | `[96, 48, 48]` | cells | Base BHAC `grid_mks.in` domain cells (snapshot block `nx` must divide `nxlone`) | `read_snapshot_footer`, `run_jobs`, `task_identity` |
| `grid.xprobmin` | `[0.17395330712343798, 0.0, 0.0]` | MKS coords | `xprobminx1..3` in `grid_mks.in` (raw BHAC angle units before $2\pi$ scaling) | `run_jobs`, `task_identity` |
| `grid.xprobmax` | `[7.824046010856292, 0.5, 1.0]` | MKS coords | `xprobmax1..3` in `grid_mks.in` | `run_jobs`, `task_identity` |
| `grid.hslope` | `0.0` | – | `hslope` in `grid_mks.in` | `run_jobs`, `task_identity` |

### `eht_imaging`

Note that the compact zero-baseline flux density (`zbl = 0.60` Jy) of the total-intensity reconstruction is fixed inside the official EHT imaging script rather than `configs/experiment.yaml`; polarimetric stages use `fixed.total_flux()` from the reconstructed Stage `I` image.

| Key | Value | Unit | Meaning and validation | Read by |
|---|---|---|---|---|
| `inputs` | `{95: [...], 96: [...], 100: [...], 101: [...]}` | – | `[filename, sha256]` of each day's HOPS `SR2_M87_2017_<day>_lo_hops_zbl-dtcal+selfcal.uvfits` | `reconstruct.py` (SHA-256 of the file bytes in `$EHT_INPUT_ROOT/eht/`) |
| `gain_tol` | `[0.02, 0.2]` | – | `[lower, upper]` fractional gain tolerances for shared-gain `selfcal` | `reconstruct.py` (`--stage I`) |
| `fsys` | `0.005` | – | Fractional systematic noise added in quadrature to `RL`/`LR` sigmas ($> 0$) | `add_fsys` in `reconstruct.py`, `refine.py`, `sample.py`; `bundle.json` |
| `solve_sites` | `['AZ', 'LM', 'PV']` | – | Stations whose D terms are solved (disjoint from `fixed_sites`) | `reconstruct.py`, `refine.py`, `sample.py` |
| `fixed_sites` | `['AA', 'AP', 'SM', 'JC']` | – | Stations pre-corrected for leakage in released HOPS `uvfits` | `load_native` metadata, `sample_dterms` |
| `p_data` | `{pvis: 1.0, m: 1.0}` | – | Data-term weights in the $P$ objective ($> 0$) | `reconstruct.py`, `convex.py`, `bundle.json` |
| `p_reg` | `{hw: 100.0, ptv: 1.0}` | – | Regularizer weights: hold-weight polarization entropy `hw` and polarimetric TV `ptv` ($> 0$) | `reconstruct.py`, `convex.py`, `bundle.json` |
| `norm_reg` | `true` | – | Normalize `ehtim` regularizers (`true` required by `convex.py`) | `reconstruct.py`, `convex.py`, `bundle.json` |
| `initial_m` | `0.1` | – | Initial fractional polarization $m_0 \in (0, 1)$ of fiducial round 1 | `reconstruct.py` |
| `initial_evpa_deg` | `0.0` | $\mathrm{deg}$ | Initial uniform EVPA $\chi_0$ of fiducial round 1 | `reconstruct.py` |
| `rounds` | `50` | – | Fiducial $P$/$D$ alternation rounds ($\ge 1$) | `reconstruct.py` |
| `restart_rounds` | `[2, 12, 22, 32, 42]` | – | Distinct rounds in `[2, rounds]` that re-initialize $P$ from blurred fractional polarization | `reconstruct.py` |
| `restart_blur_uas` | `20.0` | $\mu\mathrm{as}$ | Circular Gaussian FWHM used in `restart_rounds` ($> 0$) | `reconstruct.py` |
| `restart_domain_epsilon` | `1.0e-08` | – | Contraction $\epsilon \in (0, 1)$ scaling restarted $P$ by $1 - \epsilon$ so $\lvert P/I\rvert < 1$ | `reconstruct.py` |
| `maxit` | `500` | iterations | Per-call `ehtim` optimizer iteration cap in `fiducial` | `reconstruct.py` |
| `stop` | `1.0e-08` | – | `ehtim` optimizer convergence tolerance in `fiducial` | `reconstruct.py` |
| `shift_objective_constant` | `true` | – | Add constant $w_{\mathrm{pvis}} + w_m + h_w\ln 2$ to `ehtim` $P$ objective (minimizer unchanged) | `reconstruct.py` |
| `uv_min` | `100000000.0` | $\lambda$ | Minimum baseline length ($100\,\mathrm{M}\lambda$) excluding intra-site zero baselines before shared-gain calibration and in all $P$ imaging | `reconstruct.py`, `refine.py`, `sample.py` |
| `uncertainty_draws` | `1000` | draws/day | D-term Monte Carlo realizations per day ($4000$ total) | `sample.py`, `eht/summarize.py`, `limits.py`, `check_covariance.py` |
| `seed` | `20210405` | – | Base integer RNG seed ($\ge 0$); day $d$ uses `SeedSequence([seed, day])` | `sample.py` |
| `covariance_ddof` | `1` | – | Delta degrees of freedom ($\in [0, 3]$) for the 4-day $(\operatorname{Re} D, \operatorname{Im} D)$ covariance | `d_distribution`, `export_bundle.py`, `sample.py` |
| `residual_d_sigma` | `0.01` | – | Per-component ($\operatorname{Re}$, $\operatorname{Im}$) zero-mean Gaussian SD for `fixed_sites` | `sample_dterms`, `export_bundle.py`, `sample.py` |
| `refine.outer_rounds` | `100` | rounds | Alternating exact-$P$ / forward-$D$ rounds before the final $P$ solve | `refine.py` |
| `refine.gap_atol` | `1.0e-06` | objective | Primal–dual gap certificate tolerance per PDHG $P$ solve | `refine.py` |
| `refine.maxiter` | `20000` | iterations | PDHG iteration cap per $P$ solve in `refine.py` | `refine.py` |
| `refine.check_every` | `500` | iterations | Iteration interval between dual-gap evaluations | `refine.py` |
| `refine.time_limit` | `300.0` | $\mathrm{s}$ | Wall-time cap per $P$ solve in `refine.py` | `refine.py` |
| `refine.d_optimizer.{ftol, gtol, maxiter, maxcor, maxls, maxfun}` | `1e-13, 1e-9, 10000, 50, 50, 100000` | – | `scipy.optimize.minimize` `L-BFGS-B` options for each forward $D$ refit | `refine.py` → `fit_forward_d` |
| `monte_carlo.initial_p_fraction` | `0.05` | – | Cold-start initialization $P = 0.05\,I$ (real, $\in (0, 1)$) for every Monte Carlo draw | `sample.py` |
| `monte_carlo.gap_atol` | `1.0e-06` | objective | Primal–dual gap tolerance for certifying each draw | `sample.py` |
| `monte_carlo.maxiter` | `50000` | iterations | First-attempt PDHG iteration cap | `sample.py` |
| `monte_carlo.check_every` | `500` | iterations | Iteration interval between dual-gap checks | `sample.py` |
| `monte_carlo.time_limit` | `600.0` | $\mathrm{s}$ | First-attempt wall-time cap | `sample.py` |
| `monte_carlo.retry_maxiter` | `200000` | iterations | Second-attempt iteration cap if attempt 1 is uncertified | `sample.py` |
| `monte_carlo.retry_time_limit` | `1800.0` | $\mathrm{s}$ | Second-attempt wall-time cap | `sample.py` |

### `analysis`

| Key | Value | Unit | Meaning and validation | Read by |
|---|---|---|---|---|
| `rotation_degrees` | `108.0` | $\mathrm{deg}$ | Active rotation angle $\theta$ mapping `RAPTOR` camera coordinates onto the sky (rotates image and multiplies $Q + iU$ by $e^{2i\theta}$) | `analysis_geometry` |
| `aperture.window_degrees` | `10` | $\mathrm{deg}$ | Angular full width $w \in (0, 360]$ of each sector aperture | `analysis_geometry` → `aperture_operator` |
| `aperture.radial_max_fraction` | `1.0` | – | Fraction $\in (0, 1]$ of the half-FOV ($60\,\mu\mathrm{as}$) defining the outer aperture radius $r_{\max}$ | `analysis_geometry` → `aperture_operator`, 2D radius mask |
| `features.intensity_fraction` | `0.1` | – | Minimum daily mean brightness relative to the daily maximum ($\in [0, 1)$) | `limits.py` |
| `features.polarization_fraction` | `0.2` | – | Minimum daily mean polarized brightness relative to the daily maximum ($\in [0, 1)$; `0` selects by intensity only) | `limits.py` |
| `features.min_radius_rg` | `2.0` | $r_g$ | Inner radius cut ($7.64\,\mu\mathrm{as}$) for 2D pixel features | `limits.py` |
| `features.pixel_stride` | `5` | pixels | Subsampling stride along both spatial axes for 2D pixel features | `limits.py` |
| `features.sector_spacing_degrees` | `10` | $\mathrm{deg}$ | Angular spacing between candidate 1D sectors ($\in [1, 360]$) | `limits.py`, `paper.py` |
| `features.sector_offset_degrees` | `4` | $\mathrm{deg}$ | Starting sector index east of north ($4^\circ, 14^\circ, \dots, 354^\circ$; $\in [0, \text{spacing}-1]$) | `limits.py` |
| `features.resultant_min` | `0.7` | – | Minimum circular resultant length $R = \lvert\langle e^{2i\chi}\rangle\rvert \in (0, 1]$ required on all four days | `build_observations` |
| `cloud.coefficient` | `0.03` | – | Normalization in $b(\alpha) = \text{coefficient} \times (\alpha / \text{pivot\_alpha})^{\text{exponent}}$ ($> 0$) | `cloud_amplitude_scaling` |
| `cloud.pivot_alpha` | `0.05` | – | Pivot $\alpha$ in $b(\alpha)$ ($> 0$) | `cloud_amplitude_scaling` |
| `cloud.exponent` | `1.75` | – | Power-law exponent in $b(\alpha)$ | `cloud_amplitude_scaling` |
| `background.lag_frames` | `3` | frames | Frame lag $L \ge 1$ for axion-off EVPA differences ($3\,t_g = 26.68\,\mathrm{h}$) | `selected_background`, `paper.py` |
| `inference.phase_count` | `64` | grid points | Uniform discrete grid size $N_\phi \ge 2$ on $[0, 2\pi)$ for marginalizing $\phi_0$ | `evpa_templates` (`check_grid.py` `--phase-count` override) |
| `inference.coupling_points` | `120` | grid points | Log-spaced grid size $\ge 2$ on `[coupling_min, coupling_max]` | `calculate_observed` (`check_grid.py` `--coupling-points` override) |
| `inference.coupling_min` | `0.0072992700729927005` | – | Lower prior bound $c_{\min} = 1/137 > 0$ | `calculate_observed` |
| `inference.coupling_max` | `1000.0` | – | Upper prior bound $c_{\max} > c_{\min}$ | `calculate_observed` (`check_grid.py` `--coupling-max` override) |
| `inference.prior` | `log-uniform` | – | Prior density on $c$: `log-uniform` ($\pi(c) \propto 1/c$) or `uniform` | `normalized_posterior_cdf` |
| `inference.credibility` | `0.95` | – | Credible quantile $\in (0, 1)$ defining $c_{95}$ | `posterior_upper_limit`, `RESULTS.md` |
| `inference.exposure_hours` | `6.0` | $\mathrm{h}$ | Rectangular exposure duration $T \ge 0$ in $\operatorname{sinc}(\omega T / 2)$ | `evpa_templates` |
| `inference.day_gap_hours` | `24.0` | $\mathrm{h}$ | Intra-pair day separation $\tau > 0$ ($96 - 95$ and $101 - 100$) | `evpa_templates`, `paper.py` |
| `inference.pair_separation_hours` | `120.0` | $\mathrm{h}$ | Separation $t_2 - t_1 \ge 0$ between the two pair midpoints | `evpa_templates` |

## Machine environment variables (`.env`)

Every machine-local location derives from exactly three path roots, set in one git-ignored `.env` at the repository root (copy from [`.env.example`](../../.env.example)), exported in the shell, or passed on the command line:

1. `EHT_RAPTOR_SOURCE_ROOT` (`--source-root` / `--raptor-source-root`): the local `RAPTOR` git checkout pinned by `raptor.source`.
2. `EHT_INPUT_ROOT` (`--input-root`): read-only external inputs in the fixed [input layout](data-layout.md#input-layout) (`INPUT_LAYOUT`): GRMHD snapshots, radial tables, the official EHT `uvfits` files, feed-rotation table, and imaging script.
3. `EHT_DATA_ROOT` (`--data-root`): every generated product in the [default layout](data-layout.md#default-layout) (`LAYOUT`).

`eht_axion/environment.py` parses `.env` strictly as key–value data:

- **No shell execution or interpolation**: `$VAR` and `$(...)` are not expanded; syntax errors, unknown variable names, missing `=value`, or duplicate keys raise `ValueError`.
- **No working-directory search**: the default dotenv file is always `ROOT / ".env"` in the repository root. `--env-file PATH` replaces it with another file; `--no-env-file` disables dotenv loading.
- **Precedence**: `CLI option > exported environment (os.environ) > .env`. An explicit per-file option (`--run-root`, `--fluid-root`, `--radial-low`, `--radial-high`, `--uvfits-dir`, `--feed-rotation-file`, `--official-script`, or any `--...-root` of a stage) overrides the default derived from its root.
- **Path anchoring**: relative path values inside `.env` are anchored to the repository root. Relative CLI paths (and exported values) are resolved against the caller's working directory.

| Variable | Kind | Fallback / default | Read by | Description |
|---|---|---|---|---|
| `EHT_RAPTOR_SOURCE_ROOT` | Path | required for `build.py`, `generate.py`, `paper.py` | `scripts/raptor/build.py`, `generate.py`, `scripts/analysis/paper.py` | Local checkout of `github.com/projecteht-pku/raptor` at `raptor.source.commit` (`--source-root`; `--raptor-source-root` in `paper.py`) |
| `EHT_INPUT_ROOT` | Path | required unless every input file is given explicitly | `scripts/raptor/generate.py`, `scripts/eht/reconstruct.py`, `scripts/analysis/paper.py` | Root of the read-only [input layout](data-layout.md#input-layout): `grmhd/` (`--fluid-root`), `radial/` tables (`--radial-low`, `--radial-high`), and `eht/` (`--uvfits-dir`, `--feed-rotation-file`, `--official-script`); overridden by `--input-root` |
| `EHT_DATA_ROOT` | Path | required unless explicit `--...-root` flags are given | Every pipeline script | Root of the [default directory layout](data-layout.md#default-layout), including the run root `raptor/runs/main` of `build.py` / `generate.py`; overridden by `--data-root` |
| `EHT_RAPTOR_THREADS` | Value (`int`) | `2` | `scripts/raptor/generate.py` | OpenMP threads per `RAPTOR` worker process (`--threads`) |
| `EHT_RAPTOR_WORKERS` | Value (`int`) | `1` | `scripts/raptor/generate.py` | Concurrent `RAPTOR` worker processes (`--workers`) |
| `EHT_RAPTOR_CPU_AFFINITY` | Value (`str`) | `none` (`.env.example` uses `physical`) | `scripts/raptor/generate.py` | Worker CPU pinning policy: `none`, `physical`, or `logical` (`--cpu-affinity`) |
| `EHT_RAPTOR_CPU_LIST` | Value (`str`) | unset (all allowed CPUs) | `scripts/raptor/generate.py` | Allowed CPU IDs/ranges such as `0-47` when affinity is `physical` or `logical` (`--cpu-list`) |
| `EHT_RAPTOR_TIMEOUT_SECONDS` | Value (`float`) | unset (no per-job timeout) | `scripts/raptor/generate.py` | Wall-time limit in seconds per `RAPTOR` job (`--timeout-seconds`) |

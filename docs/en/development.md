# Development

[← Documentation index](index.md)

The repository is structured as a library folder (`eht_axion/`) paired with direct CLI scripts (`scripts/`) and a synthetic test suite (`tests/`), managed reproducibly by [Pixi](https://pixi.prefix.dev/) (`pixi.toml` and `pixi.lock`). It is not packaged as a wheel or installed into `site-packages`.

## Environments and Pixi tasks

All four Pixi environments set `MPLBACKEND=Agg`, `OPENBLAS_NUM_THREADS=2`, `OMP_NUM_THREADS=2`, and `PYTHONUNBUFFERED=1` upon activation.

| Environment | Platforms | Dependencies added beyond `default` | Purpose |
|---|---|---|---|
| `default` | `linux-64`, `linux-aarch64`, `osx-64`, `osx-arm64`, `win-64` | Python `3.12.*`, `numpy`, `scipy`, `matplotlib`, `h5py`, `pyyaml`, `pytest`, `ruff`, `python-dotenv` | Analysis scripts, `lint`, and non-`ehtim`/`astropy` tests |
| `raptor-build` | `linux-64` | `c-compiler`, `make`, `gsl`, `hdf5`, `pkg-config`, `rsync` | Compiling `RAPTOR`, running jobs, and writing the image cubes |
| `eht-summary` | All 5 platforms | `astropy` | Summarizing D-term Monte Carlo FITS draws into `eht:ensemble` |
| `eht-imaging` | All 5 platforms | `astropy`, `ehtim == 1.3.2` (PyPI) | EHT reconstruction, convex refinement, Monte Carlo, and full test suite |

| Task | Environment | Command / target |
|---|---|---|
| `raptor-build`, `raptor-generate` | `-e raptor-build` | `scripts/raptor/{build,generate}.py` |
| `eht-reconstruct`, `eht-refine`, `eht-sample` | `-e eht-imaging` | `scripts/eht/{reconstruct,refine,sample}.py` |
| `eht-summarize` | `-e eht-summary` (or `-e eht-imaging`) | `scripts/eht/summarize.py` |
| `analysis-expected-signal`, `analysis-limits`, `analysis-summarize`, `analysis-paper` | `default` | `scripts/analysis/{expected_signal,limits,summarize,paper}.py` |
| `test` | `default` or `-e eht-imaging` | `python -m pytest -q -ra -o "testpaths=tests" -o "pythonpath=."` |
| `lint` | `default` | `ruff check --isolated --target-version py312 --line-length 110 --select E,F,W eht_axion tests scripts` |
| `check` | `default` | Runs `test` and `lint` |

## Testing and linting

Run the full test suite (including `ehtim` and `astropy` tests) and the isolated Ruff linter before reviewing or committing changes:

```bash
OPENBLAS_NUM_THREADS=2 pixi run --locked -e eht-imaging test
pixi run --locked lint
```

In the `default` environment (`pixi run --locked check`), tests that import `ehtim` or `astropy` skip cleanly via `pytest.importorskip`. Every test uses synthetic fixtures in `pytest`'s `tmp_path` — no observational `uvfits` files, snapshots, or `RAPTOR` binaries are read.

Several test modules enforce repository-wide invariants rather than single-function behavior:
- `tests/test_experiment.py`: checks `configs/experiment.yaml` schema/range validation and scans the AST of `eht_axion/` (`test_library_signatures_have_no_numeric_default_values`) to forbid numeric default parameter values on public functions (outside a five-entry execution/layout allow-list).
- `tests/test_data_root.py`: verifies that every entry of `eht_axion.environment.LAYOUT` appears in the `<data root>/` tree and every entry of `INPUT_LAYOUT` in the `<input root>/` tree documented in `docs/en/data-layout.md` and `docs/zh/data-layout.md`, and checks CLI/dotenv root precedence across all scripts.
- `tests/test_docs_links.py`: verifies that `docs/en/` and `docs/zh/` have identical page sets and heading-level sequences, that every page is reachable from `README.md` and links back to its language `index.md`, and that all relative Markdown links and `#anchors` resolve.
- `tests/test_stage_safety.py` and `tests/test_provenance.py`: test symlink rejection, source/output overlap guards, dirty-tree rejection, and atomic rollback on failure.

## Codebase structure

### Top-level library modules (`eht_axion/`)

| Module | Responsibility | Key public API |
|---|---|---|
| `eht_axion/experiment.py` | Strict loader, schema/range validator, and derived views of `configs/experiment.yaml` | `read_experiment`, `validate`, `validate_ray_tracing`, `physical_scale`, `ray_tracing_config`, `analysis_geometry`, `eht_imaging_parameters`, `ANGULAR_BINS`, `EHT_DAYS` |
| `eht_axion/environment.py` | Safe reader of the single `.env`, resolution of the three roots (`EHT_RAPTOR_SOURCE_ROOT`, `EHT_INPUT_ROOT`, `EHT_DATA_ROOT`), and the input and product layouts | `resolve_environment`, `add_environment_options`, `add_data_root_option`, `add_input_root_option`, `data_root`, `input_root`, `raptor_source_root`, `layout_path`, `input_path`, `fill_defaults`, `fill_input_defaults`, `LAYOUT`, `INPUT_LAYOUT` |
| `eht_axion/provenance.py` | Transactional stage directories, `manifest.json`, git commit identity, and input-file SHA-256 | `stage_directory`, `verify_stage`, `require_declared_outputs`, `sha256`, `manifest_path`, `source_closure`, `contained_path`, `code_identity` |
| `eht_axion/physics.py` | SI/IAU constants, M87\* gravitational time/angle scales, and 211-state frequency | `PhysicalScale`, `omega` |
| `eht_axion/evpa.py` | Linear-polarization EVPA and $\pi$-periodic circular differences | `finite_evpa`, `circular_difference` |

### Subpackages

| Module | Responsibility | Key public API |
|---|---|---|
| `eht_axion/raptor/plan.py` | BHAC footer and radial-CSV validation, task list, and parameter-named task IDs | `plan_generation`, `task_id`, `read_snapshot_footer`, `validate_radial_table` |
| `eht_axion/raptor/build.py` | Clean-checkout verification and `setup.sh` compilation of pinned `RAPTOR` | `build_raptor` |
| `eht_axion/raptor/jobs.py` | Bounded concurrent `RAPTOR` process runner with deterministic slots and resume | `run_jobs`, `serialize_parameters`, `parameter_values` |
| `eht_axion/raptor/affinity.py` | Disjoint `taskset` CPU allocation and Linux/`cgroup-v2` RAM headroom checks | `allocate_cpu_sets`, `check_memory`, `launch_prefix` |
| `eht_axion/raptor/outputs.py` | Single-job HDF5 block/coordinate and integrated-spectrum validation | `validate_job_output`, `read_core_constants`, `hdf_mapping` |
| `eht_axion/raptor/convert.py` | Tiled job HDF5 assembly into float64 `(T, S, H, W, 4)` `stokes.npy` | `convert_hdf5` |
| `eht_axion/raptor/images.py` | Stokes image cubes from the job HDF5 and downstream experiment-identity checks | `write_images`, `generated_inputs`, `input_cadence_hours` |
| `eht_axion/raptor/layout.py` | Image-cube directory naming (`alpha-<α>/`, `background/`) | `alpha_label`, `product_directory`, `image_cube` |
| `eht_axion/raptor/options.py` | Shared CLI/dotenv parser of the `RAPTOR` scripts | `options` |
| `eht_axion/eht/calibration.py` | Native circular HOPS loading, `fsys`, inverse Jones D application, zero-prior forward D fit, and 4-day D sampling | `load_native`, `add_fsys`, `make_imaging_obs`, `apply_fixed_d`, `fit_forward_d`, `d_distribution`, `sample_dterms` |
| `eht_axion/eht/convex.py` | Exact convex fixed-$I$/fixed-$D$ PDHG $P$ solver with duality-gap certification | `exact_cartesian_problem`, `solve_problem`, `CartesianProblem` |
| `eht_axion/eht/ensemble.py` | Strict 4-HDU Stokes FITS reader, Paper VII global diagnostics, and per-pixel EVPA maps | `load_stokes_fits`, `global_diagnostics`, `global_values`, `map_evpa` |
| `eht_axion/analysis/results.py` | Limit-stage result files: arrays in `results.npz`, everything else in `report.json` | `write_results`, `read_results` |
| `eht_axion/analysis/expected_signal.py` | Gaussian blur + sky rotation resampling, 16-phase DFT fit per snapshot, MLE over the snapshots, and reader | `preprocess_frame`, `fit_phase_dft`, `snapshot_response`, `alpha_responses`, `expected_signal`, `processed_background`, `load_expected_signal` |
| `eht_axion/analysis/aperture.py` | Sparse CSR bilinear quadrature operator for 1D sector EVPA and brightness | `aperture_operator`, `integrate_stokes`, `sector_brightness` |
| `eht_axion/analysis/model.py` | Cloud amplitude $b(\alpha)$, circular $\alpha$ interpolation, and two-pair EVPA signal model | `cloud_amplitude_scaling`, `interpolate_fit_parameters`, `evpa_templates` |
| `eht_axion/analysis/observations.py` | Daily circular moments, population D-scatter covariances, pair differences, and seam/eigenvalue diagnostics | `ObservationData`, `moments`, `build_observations`, `branch_risk`, `covariance_diagnostics` |
| `eht_axion/analysis/likelihood.py` | Covariance modes (`std`, `cor`, `cor_fixed_background`), Gaussian $\chi^2$, and phase marginalization | `LimitOptions`, `MODES`, `calculate_observed` |
| `eht_axion/analysis/posterior.py` | Trapezoidal posterior normalization and analytical piecewise-linear CDF upper-limit inversion | `normalized_posterior_cdf`, `posterior_upper_limit` |
| `eht_axion/analysis/paper_figures.py` | Deterministic vector-PDF plotting routines for the eight paper figures | `save`, `radial_function`, `two_photon`, `signal_maps`, `template_maps`, `background`, `observational_uncertainty`, `data_vs_model`, `limits` |

### Entry-point scripts (`scripts/`)

Every script starts with a two-line `sys.path` bootstrap locating the repository root (the nearest ancestor containing `eht_axion/`), so `python scripts/<part>/<name>.py` works from any working directory. Test files mirror module and script names (`tests/test_<part>_<name>.py` for `scripts/<part>/<name>.py`).

## Extending the pipeline

### Architectural invariants

1. **Define scientific parameters once in `configs/experiment.yaml`** and validate them in `eht_axion/experiment.py`; pass them explicitly into `eht_axion/` functions without numeric defaults.
2. **Keep machine paths and execution knobs out of YAML**: use CLI flags, `add_data_root_option` / `fill_defaults`, or declared keys in `eht_axion/environment.py` (`PATH_KEYS`, `VALUE_KEYS`).
3. **Consume only verified declared stage outputs**: call `verify_stage` and `require_declared_outputs` on upstream stages, and wrap writes in `stage_directory(...)`, listing upstream `manifest_path(dir)` and consumed files in `input_paths` (never repository source files).
4. **Keep lines within 110 characters** and free of Ruff `E`, `F`, and `W` violations (`py312`).

### Adding a scientific parameter

1. Add the key and value under the appropriate section (`shared`, `raptor`, `eht_imaging`, or `analysis`) in `configs/experiment.yaml`.
2. Add the key name to `SCHEMA` in `eht_axion/experiment.py` and add type/domain checks in `_validate_shared`, `_validate_raptor`, `_validate_eht_imaging`, or `_validate_analysis`.
3. If the parameter affects `RAPTOR` ray tracing, ensure `ray_tracing_config` (and the task parameters if per-task) includes it; if it affects EHT imaging stages, ensure `eht_imaging_parameters` (`_EHT_PARAMETER_ORDER`) or solver plan dictionaries include it.
4. Thread the value explicitly from the script into the `eht_axion/` function signature (with no numeric default) and add positive/negative validation cases in `tests/test_experiment.py`.
5. Document the key in `docs/en/configuration.md` and `docs/zh/configuration.md`.

### Adding a new analysis stage

1. Add a layout key (for example `"analysis.my_check": "analysis/my-check"`) to `LAYOUT` in `eht_axion/environment.py`, and add the matching line under `analysis/` in the `<data root>/` tree of both `docs/en/data-layout.md` and `docs/zh/data-layout.md` (checked by `tests/test_data_root.py`).
2. Implement pure numerical logic in `eht_axion/analysis/` without file I/O or default scientific values.
3. Create `scripts/analysis/my_check.py` with the standard `sys.path` preamble, `add_data_root_option`, and `fill_defaults`. Verify upstream stages via `verify_stage` + `require_declared_outputs`, refuse an existing `--output-root`, and publish outputs inside `stage_directory(destination, "analysis:my-check", identity, inputs)`.
4. Register an `analysis-my-check` task in `pixi.toml`, add synthetic tests in `tests/test_analysis_my_check.py` and `tests/test_data_root.py`, and document the stage in `docs/en/` and `docs/zh/`.

# RAPTOR stages without notebooks

Module: `eht_axion.pipeline.raptor`. No notebook is executed or imported. Runtime
numerics use NumPy/SciPy and the parent `eht_axion.polarization` API, not ImageS,
matplotlib, or statsmodels. Only HDF5 conversion additionally requires `h5py`.

## Integration API

All stages return `dict[str, pathlib.Path]`. Output directories need not exist;
they are created. Existing product files are **not overwritten**, including
aliases/symlinks/hardlinks to inputs. These are low-level stage writers, **not
transactional/cache wrappers**: the orchestrator should invoke them inside its
`common.stage_directory` temporary directory, publishing only on success. A
failed low-level call can leave partial products in its output directory.

```python
from eht_axion.pipeline.raptor import (
    convert_hdf5, preprocess_cube, fit_processed, summarize_fits,
    estimate_background,
)

# Raw cube or explicitly ordered time shards, never a directory/glob guess.
preprocess_cube(
    sources, output_dir,
    sigma=2.16871, init_fov=46.2262, final_fov=32.6868,
    rotation_angle=0.6 * np.pi, final_shape=64, dtype="float64",
    state_indices=None,
)
fit_processed(
    input_dir, output_dir, fov=32.6868, operator=None,
    weighting="intensity", radial_max=None, window_degrees=10,
    angular_bins=360, dtype="float64",
)
summarize_fits(
    {0.1: "alpha-0.1/fits_2d.npy", 0.3: "alpha-0.3/fits_2d.npy"},
    output_dir, minimum_samples=2,
)
```

`preprocess_stokes(inputs, ...)` is equivalent to `preprocess_cube(sources, ...)`.
`fit_stokes(inputs, ...)` accepts processed arrays/paths/shards directly;
`fit_processed(input_dir, ...)` reads `input_dir/stokes.npy`. Both preserve the
contracts above. `preprocess_frame(stokes, ...)` and
`fit_phase_dft(evpa, phase_axis=0)` are independently callable numerical kernels.

### Input and memory contracts

- Stokes channel order is **I,Q,U,V**, axis0 of each image is **x**, axis1 is **y**.
- Raw and processed cubes: `(T,S,N,N,4)`. Preprocessing also accepts a noaxion
  time series `(T,N,N,4)`.
- Inputs may be a `.npy` path, explicitly ordered list of `.npy` time shards, or
  a NumPy array. File inputs use `np.load(mmap_mode="r", allow_pickle=False)`.
  They are never written, copied wholesale, or concatenated. Shards must have
  the same dtype and non-time shape; repeated paths are rejected. **Caller owns
  chronological order and continuity**; it is not inferred from filenames.
- Preprocessing holds one image's working arrays, fitting one time/state stack.
  Summary holds spatial accumulators, not a time cube. Background iterates
  time pairs for moments and again for histogram counts, not a copied cube.
- Real integer/floating inputs are accepted; complex/object/bool arrays are
  rejected. Output cubes/fits are explicitly float32 or float64 (default float64);
  arithmetic and summary arrays use float64. No pickled outputs.
- `state_indices=[0]` processes **only** state0 and retains axis1, yielding
  `(T,1,N,N,4)`. This avoids processing 16 unnecessary phase states for the
  res400 background. Selected states retain caller order; duplicate/out-of-range
  selections are rejected. Metadata records the original state indices.

## HDF5 conversion

```python
convert_hdf5(
    sources, output_dir,
    datasets={"I": "images/I", "Q": "images/Q", "U": "images/U", "V": "images/V"},
    image_shape=(100, 100),
    block_origins=origins, pixel_shape=(10, 10), dtype="float64",
)
```

Dataset names are mandatory, distinct and validated; HDF5 key ordering is
**never** used to infer Stokes labels. `block_origins[b]` gives source block b's
integer `(x,y)` origin **in output pixels**, not block-grid coordinates. Each
source dataset is `(blocks, pixels)` or `(blocks, pixel_x, pixel_y)`; flattened
pixels are reshaped in **C order**. Tiles must cover the entire image with no
holes or overlap. For already-dense `(H,W)` datasets omit both block parameters.

A single source path outputs `(H,W,4)`; a list outputs `(T,H,W,4)`; a rectangular
nested time/state path list outputs `(T,S,H,W,4)`. Named datasets and headers
are validated before allocating the output; tiled data are read block by block.
Products: `stokes.npy`, `convert.json`.

## Preprocessing and scientific conventions

Products: `stokes.npy`, `preprocess.json`.

1. Blur each IQUV channel with a Gaussian, nearest boundary extension.
2. Inverse spatial rotation and bilinear resampling onto centered output pixel
   coordinates `(index-(N-1)/2)*final_fov/N`. The output FOV implements cropping.
3. Active spin-2 rotation:
   `Q' = Q*cos(2θ)-U*sin(2θ)`, `U' = Q*sin(2θ)+U*cos(2θ)`; I and V are unchanged
   by the polarization rotation.

**sigma is a standard deviation in FOV units, not pixels and not FWHM**;
`sigma_pixels = sigma * input_N / init_fov`. `rotation_angle` is in radians.
No integrated-flux normalization is applied during resampling. Out-of-domain
samples clamp to the image boundary, matching the old bilinear spline behavior.

The legacy choices `sigma=2.16871`, `init_fov=46.2262`, `final_fov=32.6868`,
`rotation_angle=0.6*pi` must be supplied **explicitly**, rather than silently
inherited. `sigma=0` still executes rotation/resampling/cropping; it does **not**
enter ImageS's different historical `blur=False` code path. Constant images and
the legacy SciPy bilinear blur path are tested. No `calc_1D` work is performed
in preprocessing.

## Phase fitting and aperture alignment

Products:

| File | Shape | Columns / units |
|---|---|---|
| `fits_2d.npy` | `(T,H,W,3)` | amplitude, phi, mean EVPA; radians |
| `fits_1d.npy` | `(T,360,3)` | amplitude, phi, mean EVPA; radians |
| `fit.json` | JSON | geometry, dtype, weighting, phase and alignment conventions |

The fit accepts S=16 phase states or S=17 with noaxion state0. Only the **last
16 phase states** enter fitting, including branch selection; state0 is entirely
excluded, even when invalid. Phase samples are `ψ_k=2*pi*k/16`, with **no duplicate
endpoint**. EVPA is `0.5*atan2(U,Q)`, period pi. Zero Q and U is undefined (NaN),
not an EVPA of zero.

EVPAs are lifted about their doubled-angle circular mean and centered by the
arithmetic mean on that branch. The first DFT coefficient gives
`A = 2*abs(sum(delta*exp(-iψ)))/16` and `phi=arg(coefficient)`, corresponding to
**`mean_EVPA + A*cos(ψ + phi)`**, not `cos(ψ-phi)`. phi is in `[-pi,pi)`, mean EVPA
in `[-pi/2,pi/2)`. This branch choice is independent of state0 (unlike the old
notebook). It assumes the EVPA modulation fits on a meaningful pi-periodic
branch: very large/winding modulations are not uniquely recovered by this
first-harmonic model. A missing phase or degenerate circular mean invalidates
that location. Numerically zero amplitude (<=64 float64 eps rad) is reported as
0 with phi NaN, avoiding an arbitrary noise phase.

The 1D path calls parent APIs:

```python
operator = aperture_operator(shape=N, fov=fov, radial_max=radial_max,
                             window_degrees=10, angular_bins=360)
evpa = integrate_stokes(stokes, operator, weighting="intensity")
evpa = np.roll(evpa, -90, axis=-1)  # Before phase DFT: required EHT alignment!
```

The aperture zero direction is RAPTOR +x, whereas EHT zero corresponds to RAPTOR
+y: `output_bin[k] = aperture_bin[(k+90)%360]`. No image-axis swap is added. This
roll is also applied to 1D background EVPA. Nonstandard grids must have bin
counts divisible by four and use `-angular_bins//4`; default products have 360
bins. The fit metadata records this shift.

Default **intensity** weighting coherently integrates `I*(Q+iU)` and matches the
intended legacy-paper statistic. Optional **flux** weighting coherently integrates
`Q+iU`; it is scientifically different when I varies. The parent's explicit
centered-aperture operator is used rather than silently falling back to the
older flux-only radial helper. Callers may inject a compatible operator
`(angular_bins,N*N)`; they are responsible for its geometry/provenance.

## Time distribution / templates

`summarize_fits` processes **one spatial dimensionality per invocation**. Use
separate output directories for 2D and 1D. Alpha keys are sorted numerically;
values may be fit paths, explicitly ordered fit shards, or arrays.

`templates.npz` contains named, non-pickled arrays:

- `alphas`: `(A,)`, float64;
- `amplitude_mean`, `amplitude_std`, `phase_mean`, `phase_std`:
  `(A,H,W)` for 2D or `(A,360)` for 1D, float64;
- `amplitude_count`, `phase_count`: same spatial shapes, int64 diagnostics.

Amplitude uses Gaussian maximum-likelihood moments, population std (`ddof=0`).
Phase is centered on its unit-circle mean, lifted onto that branch, then fitted
with the same Gaussian moments; it is **not** the alternative circular standard
deviation `sqrt(-2*log(R))`. Means are wrapped to `[-pi,pi)`; fewer than
`minimum_samples` finite values or a degenerate phase mean produces NaN. Zero
amplitude samples contribute to amplitude statistics but their undefined phase
does not contribute to phase statistics. `summarize.json` records the method.

**Scientific change from older utilities:** no Lilliefors/statsmodels normality
acceptance gate or histogram-based rejection is applied to template summaries.
The transparent MLE summaries should not be interpreted as proof of Gaussianity.

## Noaxion background

```python
preprocess_cube(res400_sources, background_processed_dir,
                sigma=2.16871, init_fov=46.2262, final_fov=32.6868,
                rotation_angle=0.6*np.pi, final_shape=64, state_indices=[0])
estimate_background(background_processed_dir / "stokes.npy", output_dir,
                    lag=3, cadence_hours=26.68/3, noaxion_state=0,
                    mask=None, difference="angular", method="mle",
                    bins=100, minimum_samples=2, domain="2d")
```

`estimate_background` takes **processed Stokes**, not fit products. For an already
selected 4D `(T,H,W,4)` noaxion series explicitly set `noaxion_state=None`.
Otherwise the state index is required (default 0). A mask is a static boolean
`(H,W)` array or `.npy` path with True = included; it is not inferred from zero
values. In `domain="2d"` it selects pixels. In optional `domain="1d"` it zeroes
excluded Stokes before aperture integration; `fov` is then required and the same
weighting/aperture options as the fit are available.

Differences are **later minus earlier**, `EVPA(t+lag)-EVPA(t)`. Default `angular`
wraps each difference into `[-pi/2,pi/2)`; optional `linear` is the unwrapped
numeric subtraction of the principal EVPAs. A true zero difference is kept.
Only nonfinite/undefined values are removed. Masking is not implemented by
replacing zeros with NaN.

Default lag3 corresponds to **26.68 hours, not exactly 24 hours / one day**.
Metadata records lag, cadence and actual lag hours. The result is the scatter
of the pooled time differences, with no division by sqrt(2), spatial-mean
subtraction, or normalization. Correlated pixel/time samples are pooled; this
is not an independent-sample uncertainty estimate. Angular differences near
opposite branch boundaries can remain non-Gaussian even after wrapping.

Default `method="mle"` estimates a free global mean and Gaussian MLE sigma
(`ddof=0`) across all finite differences. This intentionally replaces the old
histogram-fit/rejection behavior. Products:

- `background.json`: numeric `sigma`, `mean`, `sigma_mle`, `mean_mle`, units `rad`,
  finite sample count, method, lag hours, mask and scientific conventions.
- `background_diagnostics.npz`: `hist_counts`, `hist_edges` (full finite range),
  scalar `count`, `mean`, `sigma`, `mean_mle`, `sigma_mle`.

Optional `method="legacy_hist"` reproduces the old normalized Gaussian PDF fit
with no free amplitude: density histogram over MLE mean +/-5 sigma, free mean
and positive sigma, and sigma covariance rejection controlled by
`legacy_threshold` (default 0.2). It requires >=40 finite samples unless scatter
is exactly zero; zero-scatter returns sigma0 without attempting a singular fit.
There is no Lilliefors gate. Diagnostics additionally include fit histogram
counts/edges and sigma error. Failure raises rather than silently emitting NaN.

## Standalone CLI

Every subcommand is available without the parent orchestrator:

```sh
python -m eht_axion.pipeline.raptor --help
python -m eht_axion.pipeline.raptor convert img0.h5 img1.h5 \
  --mapping mapping.json --states-per-time 1 --output-dir converted
python -m eht_axion.pipeline.raptor preprocess cube.npy --output-dir processed \
  --sigma 2.16871 --init-fov 46.2262 --final-fov 32.6868 \
  --rotation-angle 1.8849555921538759 --final-shape 64
python -m eht_axion.pipeline.raptor fit processed/stokes.npy --output-dir fitted \
  --fov 32.6868 --weighting intensity
python -m eht_axion.pipeline.raptor summarize \
  --fit 0.1=alpha01/fits_1d.npy --fit 0.3=alpha03/fits_1d.npy --output-dir templates1d
python -m eht_axion.pipeline.raptor preprocess shard1.npy shard2.npy \
  --output-dir bg_processed --state-indices 0 --sigma 2.16871 \
  --init-fov 46.2262 --final-fov 32.6868 --rotation-angle 1.8849555921538759
python -m eht_axion.pipeline.raptor background bg_processed/stokes.npy \
  --output-dir background --lag 3 --cadence-hours 8.893333333333333
```

`mapping.json` has `datasets`, `image_shape`, and optional `block_origins` and
`pixel_shape` keys matching conversion kwargs. CLI `--states-per-time` groups
HDF5 inputs in explicit time-major/state-minor order. Repeated `--fit` entries
for the same alpha denote ordered time shards. Each command prints a JSON
product-path mapping on success. These CLI examples illustrate usage, **not
an instruction to run large data**.

## Validation performed

`tests/test_pipeline_raptor.py`: 32 synthetic-only tests passed, including
spin-2 constant images, rotated coordinate ramps, legacy bilinear numerical
comparison, wrapped phase/sign oracle with state0 excluded, exact -90 bin
alignment, state0-only processing count, read-only mmap and source hashes,
whole-cube/shard equality, explicit HDF5 mapping, named/dtyped products,
MLE/circular summaries, zero background differences, masks/lag metadata,
legacy histogram option, alias/dtype checks, and all five CLI subcommands.
No large cube was read and no production pipeline was run by this worker.

# Script-only transport / Figure-4 bypass

`eht_axion.pipeline.transport` ports source logic from `GRMHD/main.ipynb`
(zero-based cells 7, 15, 31, 36, 43–64). Production does **not** read, parse,
execute, or require any notebook; it does not import `utils.libanaly.ImageS`.
This is an independent diagnostic stage, not an extra time-series inference stage.
The four input variants are calculated separately, not synthesized from one curve.

## Parent API and production CLI

```python
from eht_axion.pipeline.transport import run_transport_diagnostics

products = run_transport_diagnostics(
    input_dir, output_dir,
    radii=(4, 5, 6, 7, 8, 9, 10), r_width=0.5,
    fov=46.2262, sigma=2.16871, rotation_angle=0,
    angular_bins=360, mask_order="mask-before-blur",
    amplitude_method="legacy-rms",
    angular_interpolation="spin2-periodic", crossing="first-outward",
    spin=0.9375, denominator_floor=1e-12,
)
# dict[str, pathlib.Path]; no transaction/cache ownership assumed
```

Run from the project environment with NumPy, SciPy and Matplotlib installed:

```bash
cd /home/canoziia/eht-script-pipeline
python -m eht_axion.pipeline.transport \
  shared_data/analysis-inputs/draft/raptor-transport-diagnostics/snapshot-2110-res-400 \
  --output-dir /tmp/eht-transport-snapshot2110 \
  --radii 4 5 6 7 8 9 10 --r-width 0.5 \
  --fov 46.2262 --sigma 2.16871 --rotation-angle 0 \
  --angular-bins 360 --mask-order mask-before-blur \
  --amplitude-method legacy-rms \
  --angular-interpolation spin2-periodic --crossing first-outward \
  --spin 0.9375 --denominator-floor 1e-12
```

Use a fresh output directory: existing products are never overwritten. CLI prints
only the product-path JSON on success. It also works outside the checkout when
this package is installed or the checkout is on `PYTHONPATH`. No Jupyter, TeX,
external fonts, or interactive display is needed.

Inputs must be **materialized** NPY files, not Git LFS pointer text. This checkout
had pointer files during development; this stage deliberately does not fetch or
modify them. The parent/data owner supplies the actual data and runs the real
command above. No real-data result or producer configuration is claimed verified
by the synthetic tests.

## Input contract

All paths below are relative to `input_dir`, loaded with `mmap_mode="r"`,
`allow_pickle=False`; arrays are never mutated.

| File | Contract |
| --- | --- |
| `common.npy` | `(17,N,N,4)`, state0 noaxion + 16 phase states |
| `jonly.npy`, `eq_jonly.npy`, `fake.npy` | `(16 or 17,N,N,4)` |
| `common_noaxion.npy` | `(N,N,4)` or singleton `(1,N,N,4)` |
| `r_coord_400.npy` | `(N,N)`, source Boyer–Lindquist radius map |
| `axion_field_coord_400_16.npy` | `(16,N,N)`, additive EVPA rotation in radians |

N is inferred from shapes, not filenames; synthetic tests use much smaller grids.
Channel order is I,Q,U,V; axis0=x, axis1=y. Pixel centers are
`(arange(N)-(N-1)/2)*fov/N`. All angular results are radians.

The names suggest reference transport (`common`), emissivity-only (`jonly`),
equatorial emissivity-only (`eq_jonly`) and a screen/synthetic control (`fake`).
These switch meanings are **inferred, not verified** RAPTOR configuration.
Figures use filenames rather than claiming physical effects (1), (2), (3).

## Ported behavior and explicit method choices

### Radius intersections and phase conventions

`get_360_x_y(r_coord, r_s, ...)` uses the old ray orientation
`(x,y)=(-rho*sin(varphi), rho*cos(varphi))`: zero at +y, counterclockwise.
Bilinear `RegularGridInterpolator` and `brentq` find the radius crossing.
Search is restricted to the pixel-center square and `rho <= 2*r_s`.
Nonfinite **radius geometry only** is replaced by zero, as in the old source.
Missing intersections remain NaN.

Default `first-outward` scans 2049 radial positions, reports sign-change/exact-root
counts, and chooses the innermost root. Multiple crossings are explicitly saved.
The scan does not guarantee detection of tangencies or arbitrarily close roots.
`--crossing legacy-bracket` instead tests just the bracket endpoints (old behavior,
except safely clipping to the grid). NPZ saves coordinates, validity, crossing
counts and interpolated radius residuals.

DFT uses `mean + A*cos(2*pi*k/16 + phi)` with k=0..15. State0 is discarded
before the fit, including circular centering. Positive A and phi in [-pi,pi)
are supplied by the parent's `raptor.fit_phase_dft`. The producer's absolute phase
origin has **not** been independently verified. The synthetic sign/phase convention,
orientation and interpolation residuals have been tested.

### Mask and blur order (corrected, not bitwise reproduction)

Old ImageS masked only I, leaving sampled Q/U unaffected. **That path is not used.**
This implementation explicitly applies the same boolean source annulus
`finite(r_coord) & (abs(r_coord-r_s) <= r_width/2)` to **all IQUV**.
Both orders are always calculated and retained:

1. `mask-before-blur` (default primary): isolate annular emission, then blur IQUV.
2. `blur-before-mask`: blur full IQUV first, then select the source annulus.

These answer different questions and can differ substantially. `--mask-order`
selects the main figures only; both orders remain in NPZ and in the comparison
figure. Sigma=0 and the positive configured sigma are both calculated.
`2.16871` is the old ImageS `None` default, interpreted as Gaussian **standard
 deviation in native FOV units**, not pixels or FWHM; it is not a verified physical
beam configuration. Boundary extension is `nearest`; channels are never blurred
into each other. Sampling outside the map is NaN, not extrapolation.

### Native method-1 sampling and rotation

`get_raptor_evpa(...)` samples native-resolution processed Q and U at radius
coordinates with bilinear interpolation, then computes EVPA. This closely follows
`EVPA_1d_from_coord` / `oneDmethod=1`; the old method-1 actually bypassed the
64x64 rotated/regridded output. Accordingly, final_shape/final_fov and aperture
integration are not silently introduced. No extra -90-bin roll is applied because
the coordinates already start at +y.

Angular rotation remaps `varphi -> varphi-theta`, then adds theta to EVPA
(active spin-2 convention). Default `spin2-periodic` interpolates exp(2i*EVPA)
periodically to avoid branch-cut and wrap-seam artifacts; analytical amplitude
interpolation is also periodic in this mode. `--angular-interpolation legacy-angle`
reproduces the old nonperiodic scalar `np.interp` of EVPA, including its branch
artifacts. Zero polarization is undefined/NaN rather than a fabricated zero angle.

The original cell-44 analytical helper hardcodes theta=0.6*pi while cells 47/49
call transport with theta=0. The production stage uses the **same explicit
rotation** for both analytical and transport curves (default zero).
`get_analytical_A_1d(...)` retains the helper's theta=0.6*pi and nonperiodic
interpolation defaults for direct port comparisons; the stage overrides these.

### Amplitude definitions and analytical screen

The old `_fit_cos_evpa_A_only` is **RMS**, not a first harmonic fit:
`sqrt(sum(delta_EVPA**2)/8)`. It includes higher harmonics; for a 16-state array
its centering uses states[1:] but its RMS uses all 16. `legacy_rms_amplitude`
retains this unusual behavior exactly for auditability. `--amplitude-method
legacy-rms` selects this as primary (default). DFT amplitude, phase and mean are
also saved regardless of the primary selection. `--amplitude-method dft` selects
the first-harmonic amplitude for figures and ratios instead.

`get_analytical_A_1d(...)` fits the field map, then interpolates its amplitude at
radius crossings, following cell 44. This unblurred map amplitude is the ratio
denominator even for blurred transport. Denominators <=1e-12 rad by default,
nonfinite values, and missing crossings produce NaN, never infinity.

`analytical_stokes(...)` / `get_analytical_evpa(...)` separately implement cell 36:
use common state0's P and EVPA, add each field[k] directly to EVPA, reconstruct
Q/U, copy I from each input state, and set V=0 (as in the notebook model).
This analytical screen is processed by exactly the same mask/blur sampler and
saved as the fifth comparison curve, not substituted for any real variant.
`common_noaxion` supplies the source-image background; the full pixelwise
`common_state0_minus_noaxion` is saved so their equality is not assumed.

## Products and axes

```text
output_dir/
  transport.npz
  transport.json
  figures/
    amplitudes.png
    analytic_ratios.png
    mask_order_comparison.png
    source_rings.png
```

The first two figures show BL r=9 and r=4, all four real variants plus the
analytical screen, dashed/faint unblurred versus solid blurred. The order comparison
shows both processing orders for all four real variants. Values are not forcibly
clipped to the notebook's ratio limits.

The native source color map uses `log10(common_noaxion I / Imax)` plus colored
BL r=9 and photon-orbit reference curves. Color is **unblurred A(jonly)/A(common)**
in native angular labels, including when the other figures use nonzero rotation.
The photon-reference geometry retains the old formula
`2*(1+cos(2/3*acos(spin)))` mapped through `r_coord`. It is **not verified as an
observed photon ring**. As in old cell 64 its color uses the **BL r=4 ratio proxy**,
not the ratio evaluated at the photon-reference BL radius; this is labeled both
in the figure and JSON. No unverified coupling, physical switch legend, or
notebook's enlarged-polarization quiver is presented as calibrated evidence.

NPZ contains no pickle/object arrays. Main axes:

- `variants`: common, jonly, eq_jonly, fake, analytical_screen_common_state0.
- `radii`: caller order; must include 4 and 9 for source-map compatibility.
- `mask_orders`: mask-before-blur, blur-before-mask.
- `blur_sigmas`: 0, configured positive sigma.
- `amplitude_rms`, `amplitude`, `analytic_ratio`:
  `(variant, radius, mask_order, blur, angle)`.
- `fits_dft`: same axes plus `(amplitude, phase, mean_EVPA)` last.
- `analytical_map_rms`, `analytical_map_dft`: `(radius, angle)`.
- `jonly_common_ratio`: `(radius, mask_order, blur, angle)`.
- `source_color_ratio`: `(radius, angle)`, native rotation=0.
- `mask_order_amplitude_delta`: mask-before minus blur-before,
  `(variant, radius, blur, angle)`.
- `ring_xy`: `(radius, angle, xy)`; `ring_valid`, `crossing_count`,
  `crossing_residual`: `(radius, angle)`.
- Photon-reference coordinates/counts/residuals, noaxion intensity and
  common-state0-minus-noaxion are also retained.

JSON records input paths/shapes/dtypes/bytes, primary and compared methods,
finite ratio counts, missing/multiple crossings, units, phase and angular
conventions, all inferred/unverified caveats and product paths. It is a provenance
summary, not a claim of verified producer switches or a content-hash cache.

## Validation performed

```bash
python -m pytest -q -p no:cacheprovider tests/test_pipeline_transport.py
ruff check --no-cache eht_axion/pipeline/transport.py tests/test_pipeline_transport.py
```

Ten synthetic tests cover bilinear crossing residuals and cardinal orientation,
multiple/missing crossings, all-channel mask and order noncommutation, exact old
RMS versus harmonic DFT, analytical Stokes reconstruction, wrapped phase sign and
state0 exclusion, native Q/U sampling and active rotation, periodic EVPA branch
handling, readonly input hashes/NPZ/PNG/JSON, outside-checkout CLI and actionable
LFS-pointer failure. No shared-data computation or notebook execution is part of
these tests. The real production command is for the parent to execute.

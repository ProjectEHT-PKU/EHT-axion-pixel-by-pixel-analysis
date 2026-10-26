# Four-day physical-time model: research handoff

## Scope and bounded external check

Only new research code/tests are introduced; no production workflow, frozen EHT
selection or existing outputs are changed. Online check before implementation:
Wikipedia requests returned HTTP 403; the following pages returned HTTP 200 and
were read for the Gaussian transform and polarimetry context:

- https://statproofbook.github.io/P/mvn-ltt.html : explicitly proves
  `AX+b ~ N(A mu+b, A Sigma A^T)`.
- https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol : NRAO
  polarimetry guide. Spin-two conventions used here are explicitly
  `P=Q+iU=|P| exp(2i chi)`, `chi=arg(P)/2`, periodic modulo pi.

These support the mathematical construction, not the adequacy of a six-hour
GRMHD interpolant or a claim that evolving four-day predictions are validated.

## Frozen inputs and inspected convention

Primary analysis root: `/home/canoziia/GRMHD_DATA/analysis/`.

- `eht-120uas-observed-limits-v1/observations-{1d,2d}.npz`: `means` (4,N),
  `covariances` (4,N,N), `evpa_samples_radians` (4,1000,N), and frozen
  `selected_indices`: N=23 (1D), N=20 (2D). Keep these selections unchanged.
- `eht-120uas-observed-templates-v1/raptor/background/stokes.npy`:
  float64 (141,1,64,64,4), all I,Q,U,V channels, singleton state0 axis.
- Response is under **raptor/response**, not root-level response:
  `alpha-0.30/response-{1d,2d}.npz`. Keys `amplitude`, `phase`, `baseline`,
  `snapshot_labels`; shapes (141,360) and (141,64,64), labels 1600:10:3000.
- Pipeline `PhysicalScale`: label step 0.1 tg, tg=8.893247544352288 h
  at the inspected default M=6.5e9 solar masses. Labels are not hours:
  `t=(labels-labels[0])*label_step_tg*tg_hours`. Parent must verify against
  stage configuration/manifest rather than silently assuming defaults.
- Native snapshot spacing is 8.893247544 h. 141 snapshots span **140** gaps:
  1245.054656209 h, not 141 gaps (~1254 h).

## Four days and covariance

Day centers relative to Apr 5: `[0,24,120,144]` h. Consecutive contrasts:

```
H3 = [[-1, 1, 0, 0], [0,-1, 1, 0], [0, 0,-1, 1]]
H2 = H3[[0,2]]
```

The middle contrast spans **96 h**, not 24 h. The old pair midpoints are
12 and 132 h, separated by 120 h, consistent with existing `model.py`.
For fixed conditional-independent daily errors, flatten day then selected
spatial column: `C_H = (H kron I) blockdiag(C5,C6,C10,C11) (H kron I)^T`.
H3 adjacent crossblocks are **-C6** and **-C10**; first/last crossblock is zero.
Never replace this by three independent pair covariance blocks.

H3 has rank three and removes one common temporal baseline **per spatial
column**. All six day pairs are redundant (rank three); treating them as six
independent measurements overcounts. Invertible changes of three-contrast basis
preserve the Mahalanobis statistic; Gaussian density normalization changes by
the fixed Jacobian. This equivalence is local-linear, not a guarantee for
independently principal-wrapped contrasts. A single scalar baseline shared
across all spatial columns is a different nuisance model.

## Exposure and model states

`eht-script-pipeline/eht_axion/research/time_sampling.py` provides:

- `integrate_stokes_at(times, stokes, centers, width=6)` -> centers.shape plus
  untouched Stokes trailing axes. Piecewise-linear interpolation in physical
  hours, exact rectangular average, including all interior knots.
- `four_day_stokes(times, stokes, width=6, starts=None)` yields `(t0, days)`;
  each days shape is (4,1,64,64,4) for this cube. Default first center is
  first snapshot+3 h, subsequent starts follow native cadence. Both ends of
  every exposure must be in range; explicit invalid requests raise, not clip.
- `average_modulated_response(centers, coefftimes, complexcoeff, omega,
  width=6)` -> average of `z(t) exp(i omega t)`, where `z=A exp(i phase)`.
  Linearly interpolate real/imaginary coefficients, not phase angles. Uses
  analytic integration of each linear segment with stable sinc/derivative.
  Constant envelope yields `z exp(i omega center) sinc(omega width/2)`.
- `contrast_covariance(day_covariances,h=H3)` -> (K*N,K*N).
- `evpa_contrasts(chi_days)` accepts (...,4,N), returning raw H3,
  principal H3, integer winding and six direct-edge closure winding.

Stokes averages precede EVPA and aperture evaluation. For 1D use the existing
CSR `_operator`/`_integrate` with exactly the established geometry, intensity
weighting and 90-degree alignment, then frozen indices. For 2D compute selected
pixel EVPA from averaged Q,U. No source re-render or source-cube hashing is
performed by this module; provenance is stage manifest plus new-output hashes
at parent orchestration time. Parent builds Bdays (nwindow,4,N), selected
responses and synchronized background/response windows. This module does not
claim response coefficients fitted per snapshot equal re-fitting an
exposure-averaged nonlinear Stokes response; that approximation needs a separate
state0/phase-Stokes comparison.

Exposure 6 h is **underresolved** by 8.893 h cadence: exact integration of the
interpolant does not reveal missing sub-cadence variability. Matched exposure
and lag are improvements in bookkeeping, not proof of physical accuracy.
Windows overlap heavily and may be temporally correlated. Equal empirical
weights represent model-state marginalization, **not independent measurements**;
no sqrt(number of windows) error reduction or iid bootstrap is justified.

## Phase, branches, and parent responsibilities

Use `omega_hour=omega(alpha,n_principal=2)/tg_hours` and common realization
`Re[response*exp(i psi)]`, plus parent fit normalization/cloud scaling. The old
pair formula absorbs a common cosine-to-sine phase rotation and sets phase zero
at first midpoint. Uniform phase marginalization is unchanged by a common
phase shift/sign, but phase-index recovery tests must map this convention.

Do not chronologically `unwrap(2*chi)/2` across 24/96 h gaps and claim the
chosen winding is physical. Retain raw daily angles. Principal per-contrast
background deltas are usable only as an explicitly local periodic model;
report winding frequencies and six-edge closure failures. Closure winding zero
checks compatibility, not uniqueness of a physical branch. No automatic
observation/sample exclusions, selection changes or unbounded temporal drift
covariance are introduced. Parent must inspect circular concentration and
branch sensitivity, preserving day-error local branches and distinguishing
observational C from correlated model-state scatter. Background-response
crossdependence is not manufactured here; only aligned parent computations can
quantify it.

## Executed checks

`cd eht-script-pipeline && .pixi/envs/default/bin/python -m pytest -q tests/test_time_sampling.py`

**13 passed**: linear Stokes and interior-knot average, sinc exposure at zero,
small/positive/negative frequency, complex linear-envelope comparison to
independent adaptive quadrature, bounds/no extrapolation, zero-width endpoints,
window count, shared-day covariance, basis-invariant Mahalanobis, spin-two
winding/closure diagnostics.

Actual stage cube was mmap-read and averaged one window at a time (no RT):
124 valid windows, t0 first=3 h, last=1096.869447955 h. Each has full four-day
coverage. All windows completed; first-window Stokes min/max were
-4.766006691e18 / 1.660656617e19 in unchanged stage units. These extrema are
only a processing sanity check, not a physical-unit calibration or likelihood
result. No Bdays/response science artifacts or constraints are claimed here.

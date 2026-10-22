# EHT P-optimization: mathematical/runtime audit

## Bottom line / immediate handoff to budget-driver implementer

**The recorded iteration caps are real, but “just increase maxit” is not a complete convergence diagnosis. No gross analytic-gradient error was found in the tested real late-state directions. The evidence instead supports a combination of capped work, strongly nonuniform mcv coordinates, and sharp unsmoothed TV geometry.** April 11 also demonstrates that `success=True` can mean ftol termination with a substantial residual gradient.

**Next step: validate an actual same-objective, no-blur continuation on a small representative case before committing a larger production budget. Do not change solver/objective now.** This audit did not perform that continuation and therefore does not establish the needed budget or eventual convergence rate. Cartesian Q/U convexity is a useful future diagnostic/design fact, not a recommendation to replace the current solver before that test.

All new artifacts are in `/tmp/eht-optimization-audit/`. No repository or installed-package code was changed. No imaging method, D fit, or `minimize` was run. The probe replaces `scipy.optimize.minimize` with an exception, constructs actual Imagers, and evaluates objectives/analytic gradients only. No new scientific thresholds are introduced.

## 1. Completed-stage evidence

Source: existing `/home/canoziia/GRMHD_DATA/analysis/eht-paper-PD-{day}-v1/steps.json`. Counts below concern P calls, not D optimizers.

| Day | P calls | Total P iterations | Calls reporting success | Final nit / reason | Final shifted objective | Final transformed g_inf |
|---|---:|---:|---:|---|---:|---:|
| 95 | 51 | 23,720 | 5 | 500 / iteration cap | 4.281823342613123 | 8.09486e-4 |
| 96 | 51 | 24,621 | 4 | 500 / iteration cap | 4.368754136082956 | 3.73888e-4 |
| 100 (April 10) | 51 | 25,500 | 0 | 500 / iteration cap | 3.284292257772549 | 1.85596e-4 |
| 101 (April 11) | 51 | 21,915 | 10 | 34 / relative function reduction | 3.762322183375474 | 8.42958e-4 |

The successful P calls in these files terminated on relative reduction, not the 1e-8 gradient condition. April 10 final call still decreased its objective by 0.00427469374; April 11 final call decreased it by 0.0000287858892. Earlier rounds use different D values and some apply blur restarts: **do not interpret the sum of their iteration counts as a single continuous fixed-objective L-BFGS run** or their objectives as one monotonic sequence.

Existing MC `/home/canoziia/GRMHD_DATA/analysis/eht-paper-MC-reference-101-pilot10-v2/statistics.json`: 10 draws, 60 calls, 10,000 iterations, **zero successful calls**, all hitting their iteration limits. Each draw has five 100-iteration passes (blur before passes 2–5), then an unblurred 500-iteration pass. The last pass alone reduces the objective by 0.09238–0.26420; its final g_inf remains 5.89434e-4–2.03388e-3. Blur restarts raise the input objective from roughly 4–6 to roughly 23–29. These receipts establish substantial unfinished numerical work; they do not establish a necessary new science threshold or guarantee a specific added iteration count will suffice. The MC images were not gradient-probed in this audit.

Machine-readable summary: `existing-stop-summary.json`. No full source-closure hashing or rerun of the four science stages was performed.

## 2. Exact stopping/options semantics

Repo executor: `eht-script-pipeline/scripts/reconstruct_eht_official.py:image_p` constructs `ShiftedImager(pol='P', transform=['mcv'], norm_reg=True, norm_init=False, ttype='direct')`, with configured data weights `{pvis:1,m:1}`, regularizers `{hw:100,ptv:1}`, `maxit=500`, `stop=1e-8`.

Installed source actually used:

- `eht-script-pipeline/.pixi/envs/eht-rebuild/lib/python3.12/site-packages/ehtim/imager.py`
- sibling `ehtim/imaging/pol_imager_utils.py`
- ehtim distribution version **1.3.2**, SciPy **1.18.1**, Python **3.12.15**.

`imager.py:42–43,604–616` explicitly supplies L-BFGS-B with:

```
maxiter = maxit_next
ftol = stop_next = 1e-8
gtol = stop_next = 1e-8
maxcor = NHIST = 50
maxls = MAXLS = 40
jac = objgrad
bounds = None
```

Thus maxcor/maxls are **not** SciPy defaults (10/20). `maxfun` is not supplied: installed `_lbfgsb_py.py:272–274` defaults to **15000** evaluations. A long continuous solve can therefore hit maxfun independently of maxiter. Record `nfev/njev` and the actual termination message, not just nit.

SciPy upstream documentation was fetched before the mathematical/runtime probe: https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html, snapshot `scipy-lbfgsb.html`. It defines:

- ftol: `(f_k-f_{k+1}) / max(|f_k|,|f_{k+1}|,1) <= ftol`;
- gtol: maximum absolute projected gradient component <= gtol;
- these are alternative stopping conditions, not a requirement that both hold.

There are no explicit optimizer bounds here; the reported gradient is an ordinary gradient in mcv solver coordinates, not a physical Q/U stationarity certificate. The added constant `sum(data weights)+100*ln(2)` leaves gradients unchanged but **does change ftol's relative denominator**. Retain the same shifted objective in continuation comparisons. `chisq_transform=False`, `chisq_offset_gradient=0` were verified in both probes.

## 3. Real late-state central directional derivative checks

Artifacts: `probe.py`, successful final `probe.log`, `probe-100.json`, `probe-101.json`. Reproduction command, from `/home/canoziia/MyProject/PHYS`:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  eht-script-pipeline/.pixi/envs/eht-rebuild/bin/python \
  /tmp/eht-optimization-audit/probe.py
```

For each day, read its I-stage `final-fixedI.fits`, `clean-gaincal-pol.uvfits`, PD `final-raw-IQU.fits`, and `final-dterms.json`; restore the I-stage feed/flag contract with `load_gaincal_base`; apply fsys once and fixed D once using the existing core adapter. No D optimization. The late FITS I was asserted exactly equal to fixed I.

Crucially: construct an actual `Imager`, obtaining its actual **inverse-transformed** `_xarr` of shape `(4,4096)` and packed `_xinit` of length 8192, with solve mask `[0,1,1,0]`. This is **not** a physical `(I,rho,phi,psi)` tuple falsely passed as optimizer coordinates. Forward-transforming those coordinates reproduces FITS complex P (asserted). Recomputed objectives match the stored final objectives exactly at printed precision; g_inf matches to rounding. Repeated analytic gradient evaluations were bitwise identical and finite.

For every direction, record its raw norm and normalized L2 norm, analytic `g·d`, both objective evaluations, central derivative, absolute error, and relative error `|FD-g·d|/max(|FD|,|g·d|)` at **h=1e-3,1e-4,1e-5,1e-6** in JSON. All tested directions have norm 1. These are absolute solver-coordinate displacements along a direction, not per-component relative steps. Random seeds are fixed. Six directions per state: joint random, mprime-only random, phi-only random, normalized gradient, phi on lowest-I half, phi on highest-I half.

### Numerical evidence

| Day / direction | Analytic directional derivative | Relative errors at h = 1e-3, 1e-4, 1e-5, 1e-6 |
|---|---:|---|
| 100 joint random | -1.77143e-6 | 4.18e-4, 2.71e-5, 2.71e-5, 2.77e-3 |
| 100 mprime random | 1.28469e-6 | 1.96e-3, 1.32e-4, 5.75e-4, 4.45e-3 |
| 100 phi random | -2.41033e-6 | 2.27e-4, 1.20e-5, 7.10e-5, 6.61e-4 |
| 100 gradient | 6.76529e-4 | 9.51e-3, 7.80e-5, 1.02e-6, **2.66e-8** |
| 100 low-I phi | 1.04897e-5 | 5.13e-4, 1.03e-5, 3.26e-6, 4.77e-4 |
| 100 high-I phi | 7.09419e-6 | 3.95e-4, 1.14e-6, 8.13e-5, 1.42e-3 |
| 101 joint random | 1.63501e-6 | 5.94e-1, 2.20e-2, 4.66e-4, 3.87e-3 |
| 101 mprime random | 1.04468e-6 | 5.52e-3, 1.70e-4, 5.10e-4, 1.33e-2 |
| 101 phi random | 9.12938e-6 | 1.26e-4, 1.21e-5, 4.24e-5, 6.58e-4 |
| 101 gradient | 1.32811e-3 | 5.70e-1, 5.14e-1, 5.24e-2, **4.41e-4** |
| 101 low-I phi | -2.43275e-6 | 5.27e-2, 6.43e-4, 5.91e-5, 1.11e-3 |
| 101 high-I phi | -8.44287e-6 | 9.69e-1, 1.31e-1, 1.45e-3, 1.92e-4 |

The small random directional derivatives make excessively small differences cancellation-limited: some objective differences are around 1e-12. The larger errors at h=1e-3 on day 101 are **not sufficient evidence of a wrong gradient**. Term-resolved evaluation of the gradient-aligned direction localizes the changing derivative to TV:

| Day 101 weighted term | Analytic | Central h=1e-3 | Central h=1e-6 |
|---|---:|---:|---:|
| m | -0.001130960563 | -0.001130960698 | -0.001130959992 |
| pvis | -0.000631637935 | -0.000631638083 | -0.000631637187 |
| 100*hw | 0.001516765090 | 0.001516765114 | 0.001516764492 |
| ptv | 0.001573940297 | 0.000817287747 | 0.001573359526 |

TV agreement improves dramatically with smaller steps. This supports sharp local TV geometry rather than a large smooth-term chain-rule defect. **These are directional tests at two real states, not an exhaustive componentwise proof, not a verification of every polarization mode, and not a convergence certificate.** The current installed four-slot mcv/physical-gradient-slot handling was used as-is.

## 4. TV denominators and mcv conditioning

Installed TV is unsmoothed isotropic complex-P TV with zero-padded forward neighbors. `stv_pol_grad` divides directly by three norm denominators; it has no denominator epsilon or general zero-subgradient handler. At an exact flat stencil these divisions can produce 0/0. The objective itself remains defined and convex there. Back-neighbor boundary numerators are masked in current installed code; this does not resolve all possible zero denominators.

At the two actual tested images, **none of the three denominator arrays has an exact zero**. Their minima are 7.52314e-10 Jy (100) and 4.13188e-10 Jy (101), versus median forward norms 1.74691e-6 / 1.85229e-6 Jy. The gradients are finite, but proximity to flat stencils permits sharp curvature/nonsmooth limiting behavior. Do not add TV smoothing merely to make L-BFGS happier: that changes the objective.

For V=0, the actual map is

```
m = 1/2 + atan(mprime / TANWIDTH_M)/pi
P = I*m*exp(i*phi),  phi = 2*EVPA, TANWIDTH_M = 0.5
```

The local Cartesian-P Jacobian columns are orthogonal, with singular magnitudes `I*dm/dmprime` and `I*m=|P|`. This gives direct coordinate-scaling evidence:

| Quantity | Day 100 | Day 101 |
|---|---:|---:|
| m min / max | 0.000669 / 0.7130 | 0.000772 / 0.7887 |
| max absolute mprime | 237.858 | 206.160 |
| min dm/dmprime | 2.81309e-6 | 3.74464e-6 |
| radial Jacobian magnitude range (Jy) | 1.476e-10–2.261e-4 | 5.647e-11–1.836e-4 |
| angular Jacobian magnitude range (Jy) | 1.339e-8–2.090e-4 | 1.164e-8–1.855e-4 |
| largest per-pixel radial/angular anisotropy | 237.858 | 206.160 |

These are Jacobian/scaling diagnostics, **not measured Hessian condition numbers**. At these states there is no near-m=1 physical-boundary pileup. The more evident mcv saturation is **near m=0**, where mprime runs toward minus infinity and radial derivatives collapse; EVPA is also weakly identified as |P| shrinks. Thus “physical mcv concentration” should not be described as proven high-polarization-boundary saturation.

Fixed I is strictly positive, not zero: minimum 3.009e-6 / 3.846e-6 Jy per pixel, maximum 1.355e-3 / 1.405e-3. The lowest-I half of pixels carries only 4.49% / 4.80% of total I but 9.34% / 9.45% of total |P|; its median fractional m is 0.204 / 0.189. Faint fractional-polarization/EVPA “ghost” structure is therefore numerically plausible and can occupy many coordinates without much flux. This is not proof those features are artifacts, and not permission to clip/mask/raise I. For day 101 the maximum residual gradient is not in the lowest-I half, so faint pixels alone do not explain every residual.

## 5. Why fixed-I, fixed-D is convex in Cartesian Q/U

This follows from the installed expressions, not an assumption about the mcv objective:

1. `chisq_p`: weighted squared complex residual of `A P`. For `P=Q+iU`, this is a positive-semidefinite real quadratic.
2. `chisq_m`: weighted squared residual of `(A P)/(A I)`. With fixed I, each nonzero denominator is constant. It is again a linear complex model in P with a positive-semidefinite quadratic loss. Fixed observation sigmas/Jones correction are essential here.
3. The minimized HW contribution, after the harmless constant shift, is `100 * sum I*h(|P|/I) / flux`, where `h(m)=((1+m)/2)ln((1+m)/2)+((1-m)/2)ln((1-m)/2)+ln(2)`. On 0<=m<1, h is convex and nondecreasing (`h'=atanh(m)`, `h''=1/(1-m^2)>0`), so its radial composition is convex.
4. Minimized P-TV is a sum of Euclidean norms of linear neighbor differences in Q/U, including zero-padding; convex, nonsmooth at zero stencils.
5. Per-pixel physical constraints `Q_i^2+U_i^2 <= I_i^2` are convex disks. The mcv map covers the disk interior apart from exact zero amplitude at finite mprime; it is nonlinear and poorly scaled near zero. Convexity in Cartesian variables does **not** imply convexity or good conditioning in these polar/mcv coordinates.

The mathematical rules used here are standard: convex quadratic losses, norms composed with affine maps, convex nondecreasing radial composition, and convex disk constraints. Reference: S. Boyd and L. Vandenberghe, *Convex Optimization* (2004), Chapter 3, especially §3.1 and §3.2.4 (composition); official freely available text: https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf. This is an analytic application of those rules to the inspected installed formulas, not a tested new solver.

There are 2N real Q/U unknowns, not a free ratio-denominator optimization. This argument assumes the actual `chisq_transform=False`, fixed I, fixed D and weights, and finite nonzero sampled AI denominators. It does not establish joint convexity in D, validate the full alternating scientific procedure, or license replacing the objective. A future Cartesian/proximal/primal-dual approach could exploit this structure while preserving exact TV, but **no such implementation or solver switch is proposed for immediate execution**.

## 6. Budget/continuation and runtime recommendations

1. **Do a tiny representative actual continuation first**, using final D, identical fixed I/weights/objective, the delivered late image, and no blur. Compare extra work, achieved **same-D objective** reduction, gradient trend, finite values, and termination messages on day 100 versus 101. Compare image deltas both in raw Q/U/complex P and after the **same additional 19 uas kernel** on both images; do not confuse that common display kernel with an effective PSF. Use already-reviewed delta definitions/criteria, not thresholds invented by this audit. Do not require J to be monotonic **across D updates**, since those change the objective. This audit has not done that optimization; it cannot select a scientifically justified extra budget.
2. **Separate numerical budget exhaustion, ftol stop, gtol stop, and line-search/evaluation failure in the driver.** Report all, without declaring science convergence. A plain `success` flag conflates two distinct stopping tests. Do not invent a new acceptance threshold from these two probes.
3. **Prefer monitoring/checkpoint chunks within one continuous L-BFGS-B call** when the objective is unchanged: retain the 50-pair history and line-search state. Merely breaking into repeated `minimize` calls with the previous image/x is a **warm restart**, not continuation of the same L-BFGS state. `hess_inv` returned by SciPy is not an input that restores full L-BFGS-B state. If implementation uses actual short restarts, label and measure the cost honestly; do not claim state retention. The public SciPy minimize API does not offer generic resume serialization.
4. Curvature reuse is not automatically valid across new D values; those change the objective. Blur changes the iterate and intentionally sacrifices fit. Do not schedule repeated blur solely because a numerical budget chunk ended. Preserve any already-reviewed scientific restart schedule distinctly from budget management.
5. Beware hidden **maxfun=15000** in extended calls. Capture nfev, njev, wall time, objective, g_inf, status/message and actual configured options. Keep the original **maxcor=50** for the immediate continuous, no-blur budget test; do not automatically lower it to 20. Increasing maxcor is not an obvious first fix either: Imager already uses 50, not 10.
6. Runtime evidence is limited but concrete: one-thread direct Imager construction measured about **0.10 s / 0.39 s**, single objective **0.0034 s / 0.0067 s**, single gradient **0.0127 s / 0.0276 s** (100 / 101). There are 65 / 150 retained native observation rows. Both probes plus all finite differences took seconds, not an optimization run. These are single-point timings, not throughput benchmarks or production forecasts. For 8192 variables, 50 pairs alone require about 6.55 MB (decimal) for s/y vectors; matrix rebuild and lost curvature are separate restart costs. On long 500-iteration calls matrix setup may be small relative to evaluation work; for short chunks its fraction rises. **No measured claim about restart iteration overhead** can be made without continuation trials.
7. Keep science checks in physical/flux-weighted observables separate from transformed-gradient bookkeeping, using definitions/thresholds already reviewed elsewhere. Spatial EVPA in low-P pixels is intrinsically unstable; integrated complex P and raw Q/U are better diagnostics than interpreting unweighted faint-pixel EVPA changes as scientific nonconvergence.

## Artifacts and limitations

- `probe.py`: read-only, fail-closed against minimize, two real late states only.
- `probe.log`: both cases completed successfully; an initial harmless `eh.__version__` reporting AttributeError was fixed to `importlib.metadata.version('ehtim')` before the successful rerun. It never triggered imaging/optimization.
- `probe-100.json`, `probe-101.json`: actual objectives, gradients, direction norms, complete h sweeps, term decomposition, denominator statistics, conditioning diagnostics and timings.
- `existing-stop-summary.json`: concise extraction from existing completed-stage optimizer receipts.
- `scipy-lbfgsb.html`: upstream documentation snapshot; installed source is authoritative for actual defaults.

No large hashes, no all-four-day science reruns, no MC reruns, no push, no repo edits, no new scientific convergence standard. The central unresolved empirical question is **how much a true no-blur continuation buys for a measured runtime**, not whether this audit has already demonstrated a replacement solver or a sufficient final iteration count.

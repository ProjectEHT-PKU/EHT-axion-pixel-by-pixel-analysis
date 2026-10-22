# EHT optimizer budget and convergence resolution

## What has been solved

The fixed-I/fixed-D polarized image optimization is now solved with an explicit
floating-point primal-dual objective-suboptimality bound, instead of interpreting
L-BFGS-B maxiter or success flags as convergence. The objective is unchanged:
pvis1, m1, normalized HW100, exact unsmoothed isotropic P-TV1, same native data,
weights, fixed I/D and zero-exterior boundaries. No beam, data, prior, mask or
regularizer-weight tuning was used to obtain convergence.

Independent Gemini review approved this numerical solver substitution; report
`optimization-evidence/GEMINI_PDHG_REVIEW.md`. Certificates are conditional on
this fixed data/model and floating-point evaluation, not interval arithmetic,
EVPA confidence intervals or a joint P/D global-optimality proof.

## Diagnosis and actual experiments

Read-only gradient probes on real Apr10/11 final states reproduced stored
objectives/gradients and did not reveal a gross analytic-gradient error. Strong
TV curvature near flat stencils and poorly scaled low-m polar coordinates were
measured. Both P and D solvers use function-reduction stopping independently of
gradient stopping; success=True is not evidence that gtol was reached.

Apr10 SAME archived start, frozen D, no blur restarts:

|solver/budget|seconds|shifted objective|outcome|
|---|---:|---:|---|
|archived P|—|3.284292257773|500step budget result|
|mcv L-BFGS-B extra10000|133.6|3.278450921311|maxiter, gradient2.35e-4|
|mcv L-BFGS-B extra30000|392.7|3.276105873690|maxiter, gradient1.27e-4|
|exact Cartesian diagonal PDHG1500|7.61|3.270397168408|gap8.13e-7|

The runs share the same physical start/target but unequal work/budgets; these
numbers are observed work/results, not a general cross-platform speed claim.
The diagnostic Cartesian z prototype initially had an isolated-zero polar-TV
recovery defect; it was blocked, replaced with direct Cartesian gradients and
independently checked. No defective prototype was used for delivered stages.

## Method

For x=P/I, split weighted data Bx-b, exact TV T x, and radial HW entropy
`g(x)=sum(100*I/F*h(abs(x)))`. A diagonal preconditioned Chambolle-Pock/PDHG
iteration uses exact radial entropy proximal solves and unit-ball TV dual
projections. Equal diagonal steps within each TV pair preserve the projection.
Frozen absolute row/column sums prove a Schur contraction<=.99<1; no guessed
power norm, TV epsilon, phase branch choice, or variable clipping is required.

Feasible primal J and feasible dual lower bound LB give J-LB>=0. Best historical
LB remains valid for subsequent primals of the SAME fixed-D problem. Every D
update rebuilds B/b, diagonal steps, duals and certificate history from scratch.
The numeric absolute gap tolerance1e-6 is a solver accuracy target, not an EHT
scientific threshold. No claim that objective accuracy equals EVPA accuracy.

Modules and entrypoints:
- `eht_axion/research/eht_optimizer.py`: objective, adjoints, entropy prox,
  dual bounds, diagonal steps and fixed-P solver.
- `scripts/refine_eht_convex.py`: frozen-D polish, MC draw polish, coupled tests.
- `scripts/refine_eht_optimization.py`: original mcv continuation diagnostic.
- `scripts/refine_eht_d_tolerance.py`: scoped D-options/receipts diagnostic.

## All four archived D sets: P optimization complete

All at unchanged archived D, one certified P solve per day:

|day|iterations|solver seconds|final gap|
|---|---:|---:|---:|
|95|5500|55.08|8.3062e-7|
|96|9500|86.00|7.7013e-7|
|100|1500|7.61|8.1292e-7|
|101|7500|54.03|8.9655e-7|

Outputs `/home/canoziia/GRMHD_DATA/analysis/eht-convex-fixedD-<day>-v1/`.

## All ten reference MC draws: P optimization complete

Every previous Apr11 pilot draw was polished with its SAME sampled D and I.
No leakage refit or resampling occurred. All10 met gap1e-6, max gap9.7767e-7.
Solver total599.45s, mean59.95s, range11.98–115.54s in this execution, with some
concurrent jobs. Do not treat these timings as a full4000-image forecast.
Each image and numeric certificate is in
`.../eht-convex-MC101-draw0-v1/` through `.../eht-convex-MC101-draw9-v1/`.
This resolves the previous 'all pilot optimizers exhausted maxit' P-solver block.
It does not by itself validate a statistical EVPA ensemble or replace the pending
aperture/beta2 summaries.

## Separate outer D issue and current continuation

Ten extra Apr10 P/D rounds with certified P solves showed D max component deltas
falling .0019087->.0015269. Rechecking actual D fits revealed official D optimizer
explicitly uses ftol=gtol1e-6; a5step success had gradient9.14e-4, not1e-6.
Actual-state all12component finite differences agreed within1.74e-10 absolute.

A SAME-model numeric test tightened D ftol1e-13/gtol1e-9, retaining its objective,
prior-free forward Jones model and zero starting vector. D gradients improved
into~3e-6 after8iterations, but outer drift stayed essentially unchanged
(.0015267 at10rounds). Thus loose D numerical stopping is not sufficient to
explain/remove the coupling drift. Original solver/model was not rewritten.

P and D optimize different procedural targets; they have not been established
as restrictions of one common joint objective. Do not demand P-J monotonicity
across D updates or invent a new joint prior/m data term to force it. Fixed-D P
uniqueness/optimality does not establish a unique D fixed point.

An explicit100outer-round tighter-D test is running at
`.../eht-convex-tightD-100-100round-v1/`, log
`/tmp/eht-convex-validation/tightD100-100round.log`, tmux
eht-calibration-stability. A background monitor checks it every120s and will
report budget completion, P certificates, D drift and any failure. It does not
automatically declare science convergence or launch full MC. No four-dayI rerun,
raw RAPTOR recomputation, full4000 ensemble or axion inference has been started.

## Validation and provenance

Full optional regression1532passed211.12s before the final D-hook test addition.
Final focused13tests passed; new files all Ruff E/F/W pass. Tests include exact
objective comparison with real ehtim1.3.2 synthetic inputs, independent boundary
adjoints, prox/Fenchel identities, explicit Schur norm bound, valid and invalid
dual-gap checks, actual-zero optimum and scoped D-option behavior.

The first formal continuation invocation failed on undeclared threadpoolctl;
removed the unnecessary import before publishing any stage. The first convex
invocation rejected a conda-added SciPy Python file absent from wheel RECORD;
source enumeration was corrected before stage creation. Both failures are
retained as logs, not hidden reruns of published results.

All old executed sources and results are preserved. New numerical versions have
new run IDs; no old manifest was re-signed or rebound. No canonical NULL changes.

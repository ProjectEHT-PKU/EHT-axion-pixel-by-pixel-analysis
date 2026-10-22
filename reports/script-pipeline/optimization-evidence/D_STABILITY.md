# Focused D stability audit — day 100

## Scope and bottom line
Read existing coupled receipts and source; perform exactly ONE official D refit on the saved final 64×64 source, plus objective/gradient evaluations. No P solve, MC, production continuation, repository edits, or objective changes. Artifacts: `d_probe.py`, `d-probe.json`, `d-probe.log` in this directory.

**The outer coupling is still drifting. The tested D refit succeeds by function reduction, NOT by satisfying gradient tolerance. No large analytic-gradient mismatch is found at the actual saved/refitted D.** This does not establish D global optimality or instability of the physical solution.

## Existing full outer trace
Source: `/home/canoziia/GRMHD_DATA/analysis/eht-convex-coupled-100-10round-v1`.
All 11 fixed-D P calls are numerically certified, gaps 6.37e-7–9.58e-7.
`D_max_delta` is maximum absolute complex station/hand change, in fractional D.

| Call | D max change | P objective before D | Same P objective after D |
|---:|---:|---:|---:|
| 1 | .001908710 | 3.270397168 | 3.283807910 |
| 2 | .001869887 | 3.280529764 | 3.293591147 |
| 3 | .001829538 | 3.290552094 | 3.303243692 |
| 4 | .001787985 | 3.300420605 | 3.312737591 |
| 5 | .001745512 | 3.310112274 | 3.322054400 |
| 6 | .001702367 | 3.319612714 | 3.331184452 |
| 7 | .001658789 | 3.328913636 | 3.340120311 |
| 8 | .001614899 | 3.338008541 | 3.348856563 |
| 9 | .001570930 | 3.346893632 | 3.357389083 |
| 10 | .001526933 | 3.355565179 | 3.365714909 |

Call 11 is final P alignment with unchanged last D, not an eleventh D update. Raw QU change is 1.00805e-5 Jy L2 there. Net m changes .00946338→.01031409 and integrated EVPA .507308→.447337 rad across P endpoints 1→11. These are drift diagnostics, not science acceptance criteria. Every subsequent P call lowers its own initial fixed-D target; comparing different-D endpoints as a monotone descent sequence is invalid.

## Actual native D closure and one refit
Loaded the original gain-calibrated day-100 base using the driver contract, applied the existing uv cut (65 rows) and fsys once, and passed it to `core.fit_forward_d`. Native fixed input marginal errors were unchanged by this audit. Used final raw IQU, V=0, AZ/LM/PV, official direct Fourier model, explicit zero prior/gradient. Intercepted SciPy only to capture/evaluate; forwarded original arguments unchanged and restored it afterward.

**Actual options are NOT SciPy defaults:** L-BFGS-B, maxiter=10000, ftol=1e-6, gtol=1e-6, maxcor=50, maxls=40. The installed official solver explicitly supplies both tolerances. The guessed 1e-9/1e-5 defaults do not apply.

| Quantity | Value |
|---|---:|
| Initialization | all 12 real D coordinates zero |
| Objective at zero | 21.874990477655 |
| Objective at saved D, final source fixed | 1.081546926256 |
| Gradient infinity norm at saved D | .794110098336 |
| Refit objective | 1.080746779817 |
| Refit gradient infinity norm | .000914161088 |
| Success/status | true / 0 |
| Stop | RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| nit / nfev / njev | 5 / 7 / 7 |
| Optimizer-only elapsed | .00552 s |
| Max complex refit-minus-saved D | .001720739645 |

Thus **gtol is not met**: refit g_inf is about 914 times the configured tolerance. `success=True` is not a stationarity certificate. Saved D was fitted to round-10 P, whereas this refit uses the changed final-P source: its .794 gradient and D displacement mix coupling drift with optimizer stopping. They are NOT a same-image repeatability failure. No original round-10 optimizer result was saved, so its stationarity cannot be reconstructed from these receipts alone.

All 12 component central finite differences were checked at saved and refitted D. For h=1e-4, 1e-5, 1e-6, maximum absolute analytic/numeric discrepancies respectively:
- saved D: 1.31e-12, 1.94e-11, 1.64e-10;
- refitted D: 1.56e-12, 1.55e-11, 1.73e-10.
These actual-state checks strongly disfavor an f/g mismatch here; they are not a global proof. D forward RL/LR includes station-product terms, so the joint D fit is not established convex. Zero resets do not prove uniqueness or rule out other local solutions. No multistart or tightened-tolerance solve was performed.

## What alternation does—and does not—optimize
D minimizes native forward RL/LR data residuals only: sum of two normalized cross-hand chi-squares divided by two, i.e. total weighted squared residual/(4 Nvis), with explicit prior zero. P minimizes corrected-observation pvis + m + 100 HW + TV at fixed I/D. These are **not the two restrictions of a demonstrated common joint objective**. The inspected implementation realizes procedural alternating calibration/imaging, not proven mathematical block-coordinate descent. This audit does not independently establish the paper authors' mathematical intent. Do not invent a joint prior or add m to D to manufacture monotonicity. The changing corrected errors/observations also belong to each declared P problem.

Positive fixed I and strongly convex HW yield a unique exact Cartesian P minimizer for each fixed D under this model; finite-gap endpoints are approximations. This says nothing about uniqueness of D or an outer fixed point. No shape-invariance explanation is required by these facts.

## Actionable next steps (not executed)
1. Persist original D optimizer result per call: options, x0, final x/f/g_inf, success/status/message, nit/nfev/njev and source-image identity. Distinguish ftol termination from gtol; present receipts do not support a D stationarity claim.
2. If strict D stationarity is needed, first authorize a narrow same-source tighter-ftol numerical test while preserving data, weights, prior and initialization. More maxiter alone will not address an ftol stop after five iterations. Do not silently change official options.
3. More outer rounds (e.g. an authorized 50→100 budget continuation) are a numerical budget experiment, not a physics change or convergence guarantee. Current ten-round trace gives no defensible required total. Parent may continue the approved budget, preserving checkpoints; this audit starts none.
4. Record full complex D-vector changes, raw QU changes, integrated observables, and comparisons after the SAME 19 uas display kernel on a common explicitly defined bright-map support. A detector is a monitoring convention, not heuristic scientific truth or a new uncertainty threshold. No display-map detector was evaluated here.
5. Keep the existing D uncertainty prescription unchanged (including residual per-component sigma=.01). Optimizer gradients and outer deltas are not statistical D error bars. No posterior/covariance recalibration follows from this test.

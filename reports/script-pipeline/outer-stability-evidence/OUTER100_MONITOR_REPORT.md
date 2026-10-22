# OUTER100 monitor report

## Completion and scope
- Monitored existing tmux `eht-calibration-stability` at 120-second intervals until log `EXIT=0`; no restart/stop, source edits, new science, recomputation, full hashes or FITS reads.
- Case: `/home/canoziia/GRMHD_DATA/analysis/eht-convex-tightD-100-100round-v1` (initially `.building-35ibf3gp`).
- Log: `/tmp/eht-convex-validation/tightD100-100round.log`.
- Completed **100 outer P/D rounds + final P = 101 P solves, 100 D optimizer records**.
- No uncertified P, stale-log >=15-minute condition, or incomplete worker disappearance observed; normal exit. Previously completed 10 MC draws were not recomputed (completion/gap status supplied by parent, not independently re-audited here).

## Termination (verbatim)
```json
{
  "all_P_numerically_certified": true,
  "joint_D_optimality": "not verified",
  "D_refitted": true,
  "outer_rounds": 100,
  "reason": "explicit outer budget; final P aligned with last D",
  "blur_restarts": false
}
```
**Outer run ended by explicit 100-round budget, not a Ddelta scientific acceptance criterion.** All-P numerical certification does not establish joint D optimality or statistically true D/uncertainty.

## D delta comparison
| Round | recorded D_max_delta | max complex-component delta | L2 real-parameter-vector delta | D g_inf |
|---|---:|---:|---:|---:|
| 1 | 0.00190874464693 | prior D not reconstructed | prior D not reconstructed | 3.0363016e-06 |
| 10 | 0.00152667662937 | 0.0017838750053352764 | 0.0017947693611597708 | 3.26866768e-06 |
| 20 | 0.00111017951938 | 0.0012243254475358383 | 0.0012299453700383004 | 2.42688289e-07 |
| 50 | 0.000348442572571 | 0.00035991108573667815 | 0.00036307361998761447 | 2.43655631e-07 |
| 100 | 4.14979045882e-05 | 4.192522890616682e-05 | 4.243313171016846e-05 | 2.42007115e-07 |

Max complex-component change and real-parameter-vector L2 norm are different norms. The latter includes both real and imaginary coordinates, not just real parts of complex D. Values reconstructed only from consecutive small optimizer records; round 1 prior D was not assumed zero merely because optimizer initial_coordinates are zero.

## P certificates
- numerically_certified=true: 101/101 (100%).
- Termination counts: {'gap': 101}.
- Best-gap range: 6.3712209819e-07 .. 9.9788641883e-07; final P gap 8.82544068048e-07.

These are floating-point primal/historical-dual certificates, not interval arithmetic or statistical truths.

## D optimizer / gradient phase
- Options: {'maxiter': 10000, 'ftol': 1e-13, 'gtol': 1e-09, 'maxcor': 50, 'maxls': 50, 'maxfun': 100000}.
- success=true: 100/100.
- Termination messages: {'CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH': 100}.
- g_inf range: 2.42007115284e-07 .. 3.3690292254e-06; final 2.42007115284e-07.
- g_inf <= gtol: 0/100.
The optimizer success flag must not be confused with reaching gradient tolerance: relative-function-reduction termination can precede gradient convergence. Records contain endpoint gradient norms, not a full within-iteration gradient trajectory. No statistical D uncertainty certificate is established; joint D optimality remains not verified.

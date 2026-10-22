# EHT120 observed limits — completed result review

All four requested runs completed with EXIT=0 and complete manifests. Primary/native stage: `research:eht120-observed-limits`; grid/prior stage: `research:eht120-limit-numerics`.
Earlier grid startup TERM143 was parent-stopped before science output; current successful check uses source-receipt snapshots and output SHA verification, not another raw-lineage scan.
No production edits, science starts, restarts or raw rehashes were performed by this reviewer.

## Completion and identity
- primary: `eht-120uas-observed-limits-v1`, complete=True, 96 CSV rows (12 cases × 8 alpha).
- native: `eht-120uas-observed-limits-native-v1`, complete=True, 96 CSV rows (12 cases × 8 alpha).
- gridcheck: `eht-120uas-observed-limits-gridcheck-v1`, complete=True, 96 CSV rows (12 cases × 8 alpha).
- priorcheck: `eht-120uas-observed-limits-priorcheck-v1`, complete=True, 96 CSV rows (12 cases × 8 alpha).

## α=0.3: conditional 95% upper bounds in coupling c
No explicit g-GeV mapping is used. NULL controls use identical selection/covariance, zero pair differences.
| Case | Centroid observed | Matched NULL | Native observed | Native NULL |
|---|---:|---:|---:|---:|
| std-1d | 1.408437 | 0.4399731 | 1.037302 | 0.4440919 |
| cor-1d | 1.115501 | 0.4153914 | 1.132305 | 0.4087758 |
| cor_fixed_background-1d | 1.797924 | 0.6336564 | 1.618459 | 0.6332724 |
| std-2d | 0.9009636 | 0.4353386 | 0.5110435 | 0.4374047 |
| cor-2d | 1.49221 | 0.2646022 | 1.374724 | 0.2574695 |
| cor_fixed_background-2d | 1.382772 | 0.5149715 | 0.8894947 | 0.50525 |

## Independent numeric comparison (all 8 alpha, observed + NULL)
- gridcheck: maximum |Δlog10 c95|=0.0044350424 dex; settings {'phase_count': 128, 'coupling_points': 240}.
  1d: observed-std=0.00215999, null-comparison-std=0.00119341, observed-cor=0.00244087, null-comparison-cor=0.00123486, observed-cor_fixed_background=0.00139747, null-comparison-cor_fixed_background=0.00115542
  2d: observed-std=0.00130771, null-comparison-std=0.0011816, observed-cor=0.00443504, null-comparison-cor=0.00125553, observed-cor_fixed_background=0.00179959, null-comparison-cor_fixed_background=0.00118605
- priorcheck: maximum |Δlog10 c95|=0.0042112456 dex; settings {'phase_count': 128, 'coupling_points': 260, 'coupling_max': 10000.0}.
  1d: observed-std=0.00199901, null-comparison-std=0.00111285, observed-cor=0.00228522, null-comparison-cor=0.00115496, observed-cor_fixed_background=0.0012989, null-comparison-cor_fixed_background=0.00107251
  2d: observed-std=0.00121708, null-comparison-std=0.00110474, observed-cor=0.00421125, null-comparison-cor=0.00116789, observed-cor_fixed_background=0.00167154, null-comparison-cor_fixed_background=0.00110783

## Selected data and local-branch diagnostics
- primary 1d: N=23, sample counts=[1000, 1000, 1000, 1000], coordinates shape=(23, 1); selected-background σ=2.15354°. Condition numbers=[66220892.21955174, 962585.7960951008].
  Pair 1: max |Δχ|=18.855°; σ° min/median/max={'min': 1.914722904171573, 'median': 3.47333658048718, 'max': 18.59286747358356}; minimum seam distance=71.145°, minimum seam/σ=4.68518; 1σ/3σ seam counts=0/0.
  Pair 2: max |Δχ|=17.2071°; σ° min/median/max={'min': 2.880720431866061, 'median': 5.042115065653181, 'max': 14.817696200965853}; minimum seam distance=72.7929°, minimum seam/σ=5.97157; 1σ/3σ seam counts=0/0.
- primary 2d: N=20, sample counts=[1000, 1000, 1000, 1000], coordinates shape=(20, 2); selected-background σ=2.10327°. Condition numbers=[9824.560616796523, 9638.744930363753].
  Pair 1: max |Δχ|=16.4087°; σ° min/median/max={'min': 1.787217102404994, 'median': 3.220944153280359, 'max': 9.352650239386875}; minimum seam distance=73.5913°, minimum seam/σ=8.43497; 1σ/3σ seam counts=0/0.
  Pair 2: max |Δχ|=11.3396°; σ° min/median/max={'min': 2.3437256873733063, 'median': 4.880214233649347, 'max': 12.74529985893257}; minimum seam distance=78.6604°, minimum seam/σ=6.71435; 1σ/3σ seam counts=0/0.
- native 1d: N=23, sample counts=[1000, 1000, 1000, 1000], coordinates shape=(23, 1); selected-background σ=2.15354°. Condition numbers=[68873801.16533391, 1043993.681622147].
  Pair 1: max |Δχ|=16.0086°; σ° min/median/max={'min': 1.9467044829059892, 'median': 3.7760187834967534, 'max': 18.457202217250053}; minimum seam distance=73.9914°, minimum seam/σ=4.84664; 1σ/3σ seam counts=0/0.
  Pair 2: max |Δχ|=18.3375°; σ° min/median/max={'min': 2.8869665210330884, 'median': 5.1557988842111735, 'max': 14.408840836360246}; minimum seam distance=71.6625°, minimum seam/σ=6.02627; 1σ/3σ seam counts=0/0.
- native 2d: N=20, sample counts=[1000, 1000, 1000, 1000], coordinates shape=(20, 2); selected-background σ=2.10327°. Condition numbers=[11579.570892202397, 8609.245461808692].
  Pair 1: max |Δχ|=14.234°; σ° min/median/max={'min': 1.74627174770874, 'median': 3.2827776016791193, 'max': 8.899478390480501}; minimum seam distance=75.766°, minimum seam/σ=9.27272; 1σ/3σ seam counts=0/0.
  Pair 2: max |Δχ|=11.605°; σ° min/median/max={'min': 2.421280393852915, 'median': 4.800991091557028, 'max': 11.972913291110737}; minimum seam distance=78.395°, minimum seam/σ=7.15155; 1σ/3σ seam counts=0/0.
- 1D is a 360-bin theory-matched aperture, selected before likelihood; not a fabricated 350-feature mask. Final indices and actual coordinate shapes are in JSON.
- Native and centroid both have N1D=23/N2D=20 here; equal counts alone do not establish identical data or attribute all limit differences to sky shift. No axis-overlap test is invented.

## Actual fit diagnostic at α=0.3
| Case | χ²(c=0) | Best joint-grid χ² | data count |
|---|---:|---:|---:|
| std-1d | 39.29402 | 30.64711 | 46 |
| cor-1d | 158.8174 | 152.9526 | 46 |
| cor_fixed_background-1d | 25.04433 | 19.687 | 46 |
| std-2d | 37.75712 | 34.10532 | 40 |
| cor-2d | 196.0129 | 171.4305 | 40 |
| cor_fixed_background-2d | 44.89452 | 38.05914 | 40 |
The pure correlated working model has poor absolute fit, notably cor-2D; a finite posterior upper bound does not validate this Gaussian covariance model. Fixed background materially changes the fit.

## Posterior vectors (independently read NPZ)
- primary: all 12×8 posterior norms range 1–1; CDF ends at 1; monotone=True; maximum normalized last-bin mass=9.9920072e-16.
- native: all 12×8 posterior norms range 1–1; CDF ends at 1; monotone=True; maximum normalized last-bin mass=7.7715612e-16.
- gridcheck: all 12×8 posterior norms range 1–1; CDF ends at 1; monotone=True; maximum normalized last-bin mass=7.7715612e-16.
- priorcheck: all 12×8 posterior norms range 1–1; CDF ends at 1; monotone=True; maximum normalized last-bin mass=8.8817842e-16.

## Finite-MC covariance sensitivity and scope
- Auxiliary covcheck complete: first/last 500 draws, fixed centers/selection and same pooled background; 96 observed rows. Saved report maximum 0.1061015654 dex is substantially larger than grid error; full saved diagnostics retained in summary JSON. No eigenvalue repair.
- Cmax 1000→10000 (260 coupling points, 128 phases) and negligible uppermost-bin mass do not indicate upper-prior domination here. Log-uniform lower prior remains part of the conditional answer.
- Numerical grid convergence is not a 1% total-accuracy claim: finite-MC covariance, model misspecification and branch assumptions remain.
- D-model covariance omits cross-day calibration covariance; days and the two pairs are assumed independent. Scatter is not divided by √1000 and is not calibrated all-noise uncertainty.
- All six primary effective-covariance cases have zero observed features within 1σ/3σ of the branch seam. Local unwrapped Gaussian residuals still need caution for broad controls; seam diagnostics are not formal all-noise exclusions.
- Fixed-background scalar σ ignores cross-feature background covariance and is not canonical coupling-dependent signal covariance.
- Original 128-μas templates/canonical NULL are retained and unmodified. No new all-systematic exclusion or noise-calibration claim.

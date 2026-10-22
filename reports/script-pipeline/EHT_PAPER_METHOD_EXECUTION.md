# Paper VII reconstruction execution after independent Gemini review

## Review and scope

User requested Gemini review first, then execution. The independent review at
`/home/canoziia/MyProject/PHYS/reports/eht-official-audit/GEMINI_PLAN_REVIEW.md`
returned APPROVE WITH MANDATORY CONDITIONS. Required changes were implemented
before new science was launched. This is the published-method reconstruction,
not a claim of bit-identical historical pixels or five-method consensus.

New code:
- `eht_axion/research/eht_official.py`: native circular loading, error handling,
  ratio adapter, forward D fit, fixed-D correction and D-distribution sampling.
- `scripts/reconstruct_eht_official.py`: single I/fiducial/uncertainty entrypoint.
- `scripts/reconstruct_eht_fixed_d.py`: fixed-D pilot with explicit absent-station
  handling and multi-pass optimization. No leakage fit is called by this path.
- `configs/eht-official.yaml`: explicit physical/method/numerical conventions.

All previous executed sources and results remain. No RAPTOR changes, NULL root
replacement or new axion inference in this execution.

## Changes implemented

1. Execute the entire released I script with only two unique source changes:
   128->120uas from the start, and NFFT->direct Fourier backend. Keep64² and the
   .6Jy compact-flux constraint; no post-hoc I regrid or renormalization.
   Use its final5uas-restored I. Perform one shared R=L gain solve against that
   same I on original untapered native circular data, excluding uv<1e8 first.
   Gain application retains separate RL/LR errors. No I2% error inflation is
   carried into the P data. Station tables matched by name, not position.
2. Reuse the official full nonlinear forward D model and analytic gradient with
   two exact source adaptations: explicit zero prior and zero prior gradient.
   Upstream function source SHA and unique replacement anchors are validated;
   source adaptation receipt is available via `forward_solver_identity()`.
   Fixed native measurement errors, V=0, no hidden leakage_tol=.1 penalty.
3. Native circular values/weights throughout P/D: no irreversible Stokes
   roundtrip. Cross-hand sigma² adds(.005|VI|)² exactly once. Parallel weights
   remain gain-calibrated thermal values. Local ratio adapter uses
   Var(I)=(sigmaRR²+sigmaLL²)/4. Jones propagation retains marginal errors, not
   complete covariance; this limitation remains explicit.
4. P objective keeps normalized HW100/PTV1/Pvis1/m1. Numerical stop1e-8, explicit
   constant shift removes irrelevant negative entropy/data offsets without
   changing gradients. maxit500 is a numerical work ceiling, not a paper
   acceptance threshold. All iterations/gradient norms/failures recorded.
   Fiducial P/D budget50.20uas restarts have a declared numerical schedule;
   fixed-I restart uses an I-weighted fractional complex-P interpolation to stay
   in m<1, not clipping fitted outputs or changing the objective.
5. Separate Appendix H sampling path. Published Table5 reference and rebuilt
   four-day D are distinct source tracks, never mixed. Default ddof1; residual
   AA/AP/SM/JC complex Gaussian per-real-component SD.01, matching library
   convention. Same-day absent station draws are recorded, not applied to a
   nonexistent station. Fresh base -> fsys once -> sampledD correction -> P-only.
   No thermal bootstrap, extra hyperparameter sweep or D refit within draws.

Gemini linked Paper IV normalization factors to the current norm_reg=True
implementation. This resolves the earlier 'unknown normalization' blanket block;
keep exact formulas/units recorded rather than change beta to fit a desired plot.

## Preflight

-Full optional-environment suite:1520passed195.07s,EXIT0.
-35 targeted new core/driver/old-wrapper tests passed; new sources Ruff passed.
-31 focused tests including absent-station/multi-pass fixed-D wrapper passed.
-Core includes direct matrix Jones/noise oracle, analytic vs finite-difference
  gradient, known injected-D recovery, zero-D recovery and source identity tests.
-Input inspection: all4 released HOPS files have finite native four-product data;
  the new complete-row validity policy dropped no rows before scan averaging.
-Four-file original input receipt:
  `/home/canoziia/GRMHD_DATA/EHT_DATA/official-2023-D01-01/four-day-inputs-verified.json`.
-Original source identities checked against completed I stage; all match current
  source bytes. No post-hoc source rewriting or old-manifest resigning.

## Actual execution status

All four I stages completed EXIT0:
`/home/canoziia/GRMHD_DATA/analysis/eht-paper-I-{101,95,96,100}-v1/`.

Fiducial controller runs sequentially101,95,96,100:
`.../eht-paper-PD-<day>-v1/`. Apr11(101) completed EXIT0, remaining days queued/running.
Each stage is published only after completion.
Logs `/tmp/eht-official-paper/PD-<day>.log`, tmux `eht-paper-pd`.
Mid-Apr11 D values already differ materially from the old penalized experiments;
not a final scientific comparison and not evidence one prior alone caused changes.

Fixed-D reference pilot:
- v1 rejected an unused JC draw before P optimization (Apr11 has no JC).
  Failure log preserved; it did not modify any completed stage.
- v2 wrapper filters only station draws absent from that day's observation and
  records them in each optimizer row. Same full7-station draws remain saved.
- v2 executes10 samples with5 P passes(maxit100) and final unblurred polish500.
  Root `.../eht-paper-MC-reference-101-pilot10-v2/`.
  Log `/tmp/eht-official-paper/MC-reference-101-pilot10-v2.log`.
- v2 completed10draws EXIT0. All10 final optimizers reached their iteration
  ceiling; final transformed-gradient infinity norms range .000589 to .002034.
  The receipt records every unsuccessful draw. This is NOT an accepted final EVPA
  uncertainty ensemble. Full1000/day is NOT queued until optimization/schedule
  checks pass. Program EXIT0 does not imply optimization convergence.
- The inherited sampling-receipt text still says single P call; actual v2 schedule
  is declared in plan.fixed_D_schedule and per-optimizer p_only_pass1..6 records.
  This metadata inconsistency is documented, not patched in executed source or
  re-signed. Use actual trace and declared wrapper plan for execution identity.

The second Gemini code review initially hit a retryable Cloudflare520; after
backoff/retry it COMPLETED and wrote GEMINI_IMPLEMENTATION_REVIEW.md. Verdict:
I/fiducial P/D APPROVED with no remaining execution blocker; fixed-D pilot wrapper
valid but full1000 ensemble GATED pending optimizer-quality checks and summaries.
The report is copied into eht-paper-evidence/; no scientific source was changed.

## Remaining outputs to implement/validate

Sampling/core and fixed-D imaging exist. The current new executor reports whole
image net/resolved polarization and integrated EVPA, and saves every raw/display
IQU map. Paper-style360-step10degree aperture profiles and beta2 statistics are
explicitly deferred; existing research code is available but not yet bound into
this new source chain. Do not label currentpilot as the complete published figure
or final daily covariance. Full samples require optimization acceptance and then
these periodic-angle/spatial summaries. No inference before validation.

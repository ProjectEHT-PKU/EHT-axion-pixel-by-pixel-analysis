# Independent Gemini Review: Paper VII Polarimetric & Appendix H Reconstruction Plan

**Verdict: APPROVE WITH MANDATORY CONDITIONS**
`IMPLEMENTATION_PLAN_ZH.md` and `CODE_AUDIT.md` accurately diagnose the 5 methodological divergences in the experimental scripts. Main may proceed to implementation and scientific execution **immediately upon enforcing the concrete fixes and locked conventions below**. Do not stall on unrecoverable historical script minutiae.

---

## 1. Prioritized Concrete Fixes (Gating Implementation)

### Fix 1: Stokes $I$ Grid (`120 uas` from Start), Taper Domain & Single Post-$I$ Self-Cal
* **Evidence**: Paper VII p.22 C.2 specifies fixing $I$ FOV at $120\,\mu\text{as}$ ($64\times 64$, $F=0.6\,\text{Jy}$) during the Paper IV `eht-imaging` script, then "(re)self-calibrating the station amplitude and phase gains (assuming $G_R=G_L$) to our final Stokes $\mathcal{I}$ image", and dropping zero baselines ($\text{uv\_min}=10^8\lambda$) for polarimetric imaging.
* **Domain & Gain-Lineage Resolution**:
  1. In `official-I-original.py`, replace `fov = 128*eh.RADPERUAS` with `fov = 120*eh.RADPERUAS` from the start. Eliminate post-hoc `regrid_image(120, 64)` and flux renormalization.
  2. In `official-I-original.py`, `im_out = imgr.out_last().blur_circ(5.0*eh.RADPERUAS)` is expressly constructed to match the **original untapered** data (`obs_orig`), whereas `imgr.out_last()` lives in the $5\,\mu\text{as}$ reverse-tapered domain (`obs_sc_init`). Pairing `out_last()` with untapered `obs_orig` was a genuine domain mismatch.
  3. **Do NOT compound intermediate selfcal tables**. The released HOPS UVFITS (`zbl-dtcal+selfcal`) already contain Paper IV SMILI gains (VII p.22). Inside the Paper IV script, the Round 4 `caltab` is solved directly against `obs_sc_init` (from scratch, not incrementally). However, `obs_sc_init` has zero-baseline rescaling (`uv < 0.1e9`), $5\,\mu\text{as}$ reverse taper, and $2\%$ fractional $I$ noise.
  4. **Exact Implementation**: Take the restored $120\,\mu\text{as}$ final image `im_out`. Prepare clean untapered scan-averaged data in `polrep='circ'`. To prevent short baselines ($<0.1\,\text{G}\lambda$, where extended jet flux $\sim 0.6\,\text{Jy}$ is missing from the compact $0.6\,\text{Jy}$ model `im_out`) from corrupting station amplitude gains, either:
     - Perform the single post-$I$ `eh.selfcal(..., method='both', gain_tol=[0.02, 0.2], sites=['SM','JC','AA','AP','LM','SP','AZ'], caltable=True)` on the untapered data with short baselines rescaled (`rescale_zerobaseline` on $uv < 10^8\lambda$, or after flagging $uv < 10^8\lambda$), and apply that `caltab` (`rscale == lscale`) to the un-inflated, untapered `obs_orig` (`polrep='circ'`), then flag `uv_min=1e8`; **OR**
     - Apply the Paper IV Round 4 `caltab` (aligned by station name) to `obs_orig` (`polrep='circ'`), flag `uv_min=1e8`, and run one final `selfcal(..., method='both', gain_tol=[0.02, 0.2])` against restored `im_out` (check gain shift is near unity). Lock `fixed_I = im_out` ($120\,\mu\text{as}$ restored) for all subsequent $P$ and $D$ steps.

### Fix 2: Full Nonlinear Forward-Jones $D$-Term Solver (`pol_cal_new`) Without Prior
* **Evidence**: Paper VII p.23 C.2: *"we iteratively solve for the remaining D-terms by minimizing the $\chi^2$ between the real (gain-calibrated) data and sampled data from the current image reconstruction corrupted with Jones matrices (Equation (4)). We do not use any linearized approximations... throughout we assume the model image has no circular polarization ($\mathcal{V}=0$)."*
* **Mathematical & Code Verification of `pol_cal_new.leakage_cal_new`**:
  1. Legacy `pol_cal.leakage_cal(inverse=False)` inverse-corrects data and recomputes $\sigma$ at each trial $D$ (in `obs_simulate.py:1115`, `tarr['sefdr'] == 0` triggers fallback to data sigmas, making the denominator $\sigma(D)$ depend on $D$ and biasing $D\to 0$ alongside the `abs(D/0.1)^2` penalty). Legacy `inverse=True` ignores the image when `const_fpol=False`.
  2. `pol_cal_new.leakage_cal_new` implements the exact nonlinear forward RIME for `pol_fit=['RL', 'LR']` with fixed `sigma_RL`, `sigma_LR` and $\mathcal{V}=0$ (`im.vvec = 0` gives `ft_RR == ft_LL == ft_I`).
  3. **Analytic Gradient Verified**: In `pol_cal_new.py:216-315`, for $f = \frac{1}{2 N_{\text{vis}} N_{\text{pol}}} \sum |V - V_{\text{leak}}|^2 / \sigma^2$, the Wirtinger derivatives with respect to $\text{Re}(D_{R1}), \text{Im}(D_{R1})$ (where $V_{\text{leak}}$ is linear in $D_1$ and $D_2^*$) have signs `[-Re(resid * dV), +Im(resid * dV)]` for site 1 and `[-Re(resid * dV), -Im(resid * dV)]` for site 2 (conjugate), divided by `Nvis * len(pol_fit)` (the factor of 2 cancels). This is mathematically exact.
  4. **How to Disable the Unsupported Prior**: Passing `leakage_tol=np.inf` in IEEE-754 float64 evaluates `Dpar / (np.inf**2) == 0.0` cleanly whenever `Dpar` is finite, **but** a cleaner, zero-risk wrapper or local solver function calling `chisq_total` and `chisq_total_grad` directly (or setting `prior = 0.0; priorgrad = 0.0` when `leakage_tol is None or np.isinf(leakage_tol)`) is **mandatory** so line-search steps cannot hit `inf/inf = nan`.
  5. **Required Unit Checks before 4-Day Run**: (a) Finite-difference vs analytic gradient check (`rtol < 1e-6`); (b) Zero-leakage synthetic recovery ($D_{\text{fit}} \approx 0$); (c) Known injected $D$ recovery on synthetic data with `V=0`.
  6. **Note on `apply_jones_inverse` after $D$-fit**: When `apply_solution=True` produces the $D$-corrected observation for the next `make_image_P` round, `apply_jones_inverse` mixes `rrsigma, llsigma, rlsigma, lrsigma` by $|J^{-1}|^2$. Because `leakage_cal_new` is always called on a fresh copy of `noisy_base` (before $D$-correction), this does not affect the $D$-fit, only the next $P$ step's `psigma`/`msigma`.

### Fix 3: Native Circular (`polrep='circ'`) Noise Path & $f_{\text{sys}}=0.005$ Convention
* **Evidence**: Paper VII p.22 Eq.(C3) & p.31 Appendix G.1: *"we add a term equal to $f_{\text{sys}} \times |\tilde{\mathcal{I}}|$ in quadrature to our baseline thermal noise estimates, where $f_{\text{sys}}$ is an overall multiplicative factor"* ($f_{\text{sys}}=0.005$, Table 6).
* **Circular vs Stokes Pitfall Discovered in `ehtim 1.3.2`**:
  1. Loading/converting `circ -> stokes -> circ` (`obsdata.py:302-336`) replaces separate `rlsigma`, `lrsigma` with $\sqrt{(\sigma_{RL}^2+\sigma_{LR}^2)/2}$, erasing unequal cross-hand weights and single-pol flags (e.g. JCMT).
  2. **CRITICAL `ehtim` Quirk**: If `obs` is kept in `polrep='circ'` when passed to `Imager(..., pol='P')`, `obs.unpack(['pvis','psigma'], conj=True)` (`obsdata.py:863-876`) unpacks `rlvis` and `conj(lrvis)` with their native `rlsigma` and `lrsigma`! And `obs.unpack(['m','msigma'], conj=True)` calls `obsh.merr2(rlsigma, rrsigma, llsigma, I, m)` (`obs_helpers.py:701`), which computes:
     $$\sigma_m = \sqrt{\frac{\sigma_{RL}^2 + (\sigma_{RR}^2 + \sigma_{LL}^2)|m|^2}{|I|^2}}$$
     Notice `merr2` is missing the factor of $1/4$ on $(\sigma_{RR}^2 + \sigma_{LL}^2)$ because $\sigma_I^2 = \frac{1}{4}(\sigma_{RR}^2 + \sigma_{LL}^2)$ (compare `merr` in Stokes which uses `sigma**2 = 0.25*(rrsigma**2 + llsigma**2)`). Since $|m|\sim 0.1$, this difference is tiny ($\sim 0.3\%$), **provided** `rrsigma` and `llsigma` are NOT inflated by $f_{\text{sys}}$.
  3. **Locked Noise Rule**: Keep `base` in `polrep='circ'` through UVFITS load (`polrep='circ'`), `avg_coherent`, and `caltab.applycal`. Inject $f_{\text{sys}}=0.005$ **once** on cross-hands only:
     $$\tilde{I} = \tfrac{1}{2}(V_{RR} + V_{LL}), \quad \sigma_{RL,\text{new}} = \sqrt{\sigma_{RL}^2 + (f_{\text{sys}}|\tilde{I}|)^2}, \quad \sigma_{LR,\text{new}} = \sqrt{\sigma_{LR}^2 + (f_{\text{sys}}|\tilde{I}|)^2}$$
     Keep `rrsigma`, `llsigma` at their gain-calibrated thermal values. Pass this `polrep='circ'` observation directly to both `leakage_cal_new` and `Imager(..., pol='P')` (note `Imager.make_image` only switches `prior_next` and `init_next` to Stokes, leaving `obs` in `circ`, so `chisqdata_pvis` and `chisqdata_m` will use native `rlsigma`/`lrsigma` via `data_conj()`!).

### Fix 4: Regularizer Normalization (`norm_reg=True`) & Optimizer Stopping Tolerance
* **Published Regularizer Definitions Verified**:
  1. Chael et al. (2016) Eqs.(13), (15), (C4), (C7) define unnormalized $S_{\text{HW}}$ ($m_{\max}=1$) and $S_{\text{TV}}$.
  2. Paper IV p.33 Appendix A Eqs.(27)–(35) explicitly defines the `eht-imaging` library normalization factors $\zeta$: for entropy terms $\zeta = F$, and for TV terms $\zeta = F \frac{\Delta\theta}{\Theta_{\text{beam}}}$ (`norm = flux * psize / beam_size` in `imager_utils.py:2291` and `pol_imager_utils.py:986`, where `NORM_REGULARIZER = True`).
  3. Therefore, Paper VII Table 6 weights `hw: 100, ptv: 1, pvis: 1, m: 1` were surveyed **with `ehtim`'s native `norm_reg=True` convention**. **Lock `norm_reg=True`**; do not toggle `norm_reg=False` or invent custom $\beta$ rescalings, and remove `--norm-reg-audit-pending` as a blocker.
* **Stopping Criterion Fix (`stop`)**:
  1. With `norm_reg=True`, $-100 S_{\text{HW}}/F \approx -100\ln 2 \approx -69.31$, so `J_pol` sits near $-64$. In `validate_eht_apr11.py`, `stop=1e-4` sets L-BFGS-B `ftol=1e-4`, causing `(|J_k - J_{k+1}|)/|J| <= 1e-4` to trigger after **`nit=1`** whenever $|\Delta J| < 6.4\times 10^{-3}$ (seen in rounds 50–100 of all three validation runs).
  2. Shifting `shw` by $+\ln 2$ (so $S_{\text{HW}}(m=0)=0$) or setting `stop=1e-8` (or `1e-10`) inside `Imager` prevents artificial 1-step relative-`ftol` exits while `maxit=100` bounds per-round work. **Action**: Set `stop=1e-8` (or subtract the constant offset in `objfunc`) for `Imager(..., pol='P')`, lock outer budget at Paper VII's `n_iter = 50` (`maxit = 100`), and record actual `nit` and gradient norms. Do **not** add self-invented early-stop rules (`map_rtol`/`dterm_atol`).

### Fix 5: Appendix H.1 Monte Carlo Executor (Fixed-$D$ Only)
* **Evidence**: Paper VII p.33–34 Appendix H.1 & Fig. 24:
  - Draw 1000 random $D$-term sets per observing day.
  - For AZ (SMT), LM (LMT), PV: draw $D_R, D_L$ from 2D real-imag normal distributions with mean and $2\times 2$ covariance computed across the 4 observing days (`eht-imaging` includes $\text{Cov}(\text{Re}\,D, \text{Im}\,D)$, independent across stations and $R/L$ hands).
  - For AA, AP, SM, JC: draw residual $D_R, D_L$ from zero-mean circular complex Gaussians with $1\%$ standard deviation.
  - Apply the drawn 7-station $D$-terms once via `apply_jones_inverse` to the day's gain-calibrated base observation (with $f_{\text{sys}}=0.005$ added once).
  - Reconstruct $Q, U$ keeping $I$ fixed ($\mathcal{V}=0$) **with NO $D$-term calibration loop**.

---

## 2. Locked Implementation Conventions (No Further Literature Fishing)

1. **Two Distinct Sampler Tracks (Never Mix)**:
   - **Track A (`published_table5_reference`)**: Mean and $2\times 2$ covariance of $(\text{Re}\,D, \text{Im}\,D)$ across the 4 days from `/home/canoziia/GRMHD_DATA/EHT_DATA/official-2023-D01-01/reference/paperVII-table5-ehtim.json`.
   - **Track B (`rebuilt_fiducial_selfconsistent`)**: Mean and $2\times 2$ covariance across the 4 days from the local 4-day fiducial rebuild using Fixes 1–4.
   - *Rule*: Never pair Track A covariances with Track B means (or vice versa).
2. **4-Day Sample Covariance `ddof`**:
   - **Default**: `ddof=1` (`np.cov(..., ddof=1)`, unbiased 4-epoch sample covariance). Log `ddof=0` ($3/4$ scale) as a secondary metadata diagnostic on the drawn $D$ cloud only; do not double the 4000-image MC runs.
3. **1% Circular Complex Gaussian Convention (AA, AP, SM, JC)**:
   - **Default**: Match `ehtim.observing.obs_helpers.cerror(sigma=0.01)` (`obs_helpers.py:709`):
     $$\text{Re}(D) \sim \mathcal{N}(0, 0.01^2), \quad \text{Im}(D) \sim \mathcal{N}(0, 0.01^2) \quad (\text{per-component SD} = 1\%)$$
     Document that total complex modulus RMS $\sqrt{\mathbb{E}|D|^2} = 0.01\sqrt{2}$ under this `ehtim` convention (whereas complex RMS $= 1\%$ would divide per-component SD by $\sqrt{2}$).
4. **Fixed-$D$ MC Imaging Schedule**:
   - In each MC realization, $D$ is fixed, so the $P$ objective landscape is static. Run fixed-$I$ `make_image_P` (`pvis: 1, m: 1, hw: 100, ptv: 1, norm_reg: True, transform: ['mcv'], ttype: 'direct'`) initialized at $m=0.05, \chi=0, V=0$. Because $D$ does not update between rounds, 50 outer rounds without blur restarts is simply 1 continuous L-BFGS-B optimization (or a short sequence of $20\,\mu\text{as}$ blur-restart polish steps, e.g., 5 cycles of `maxit=100` as in Chael 2016 Fig. 1, or 50 rounds matching the fiducial script with $D$-solve omitted). Verify on a 10-sample pilot that the chosen fixed-$D$ schedule converges cleanly (`nit`, objective decrease) before launching $4\times 1000$ samples.
5. **Post-Processing & Paper VII Metrics (p.10, 13–14 Eqs. 10–14)**:
   - Convolve $I, Q, U$ with the common $19\,\mu\text{as}$ FWHM circular Gaussian kernel (`blur_circ(19*eh.RADPERUAS, 19*eh.RADPERUAS)`).
   - Compute $|m|_{\text{net}}$ (Eq. 12), $\langle |m|\rangle$ (Eq. 13), $\beta_2$ amplitude and phase $\angle\beta_2$ (Eq. 14, centered on Stokes $I$), and azimuthal $\langle\mathcal{P}\rangle(\varphi)$, $\langle\chi\rangle(\varphi)$ in $10^\circ$ windows stepped by $1^\circ$ (Eqs. 10–11).

---

## 3. Action Sequence Approved for Main

1. Implement the 4-day $120\,\mu\text{as}$ Stokes $I$ + consistent post-$I$ gain-cal pipeline and run all 4 days (`095, 096, 100, 101`).
2. Implement the prior-free forward-Jones $D$ solver (`pol_cal_new` without prior) + native `polrep='circ'` noise path + `stop=1e-8`, run unit verification tests, and execute the 4-day fiducial $P/D$ reconstruction (50 rounds). Compare recovered $D$-terms against Table 5.
3. Run a 10-sample fixed-$D$ pilot on Apr 11, then execute the Appendix H.1 Monte Carlo ensembles into new stage directories (never overwriting existing `analysis/` or `limits` outputs).

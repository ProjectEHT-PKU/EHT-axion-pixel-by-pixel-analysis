# Numerical Review: Exact Fixed-I/Fixed-D Diagonal PDHG Solver

**Decision**: **APPROVE**

## 1. Executive Summary & Priority Actions
- **Status**: **APPROVED** as a drop-in numerical solver replacement for `mcv` L-BFGS-B on the fixed-$I$, fixed-$D$ linear polarization subproblem.
- **Classification**: Pure numerical optimization change on the exact same objective; no change to paper physics, regularizers, weights, or data model.
- **Approved Production Scope**:
  1. **4-day fixed-$D$ polish + MC10**: Approved for immediate port and execution.
  2. **Coupled outer $D$-update loop**: Approved, **provided** that after every outer $D$ update, the operator $B(D)$, data vector $b(D)$, diagonal step sizes $(\tau, \sigma_{\text{data}}, \sigma_{\text{TV}})$, and primal-dual certificate are recomputed from scratch for the newly defined fixed-$D$ subproblem.

## 2. Objective & Mathematical Equivalence Audit
1. **Exact Objective Match**:
   - Matches the official fixed-$I$, fixed-$D$ objective up to the exact constant shift $+2 + 100\ln 2$: data terms `pvis=1`, `m=1` (with exact $A_m = \operatorname{diag}(A I)^{-1} A$ and variance/length normalization), regularizers `hw=100`, `ptv=1`, `norm_reg=True`, `flux = sum(I)`, `tvnorm = flux * psize / beam_size`, and no polar reference prior.
   - Forward-difference PTV operator $D(P)$ preserves the upstream zero-exterior boundary on the right column and bottom row (`dx[:, :-1] += p[:, 1:]`, `dy[:-1, :] += p[1:, :]`).
2. **Direct Cartesian Gradient (No Polar Recovery)**:
   - For $F_{\text{HW}}(P) = \frac{100}{\text{flux}}\sum_i I_i h(|P_i|/I_i)$ with $h(m) = \frac{1}{2}((1+m)\ln(1+m) + (1-m)\ln(1-m))$, the chain rule $\frac{\partial m_i}{\partial |P_i|} = \frac{1}{I_i}$ cancels $I_i$, giving $\nabla_P F_{\text{HW}}(P)_i = \frac{100}{\text{flux}}\operatorname{atanh}(m_i)\frac{P_i}{|P_i|}$.
   - In normalized coordinates $x_i = P_i / I_i$, $g(x) = \sum_i w_i h(|x_i|)$ with $w_i = \frac{100 I_i}{\text{flux}}$ has gradient $w_i \operatorname{atanh}(|x_i|)\frac{x_i}{|x_i|} = I_i (\nabla_P F_{\text{HW}})_i$. Verified analytically and by finite differences (`6.95e-10`). No polar gradient transformation is used.
3. **Inner Products, Adjoints, Conjugate, and Proximal Maps**:
   - Real inner product $\langle u, v \rangle_{\mathbb{R}} = \operatorname{Re}(\operatorname{vdot}(u, v))$ is used consistently across primal and dual complex spaces; adjoints $B^H$ and $T^H(r) = I D^H(r)/\text{tvnorm}$ match to $<10^{-13}$.
   - Fenchel conjugate $g^*(v) = \sum_i w_i \ln\cosh(|v_i|/w_i)$ is evaluated overflow-safely via `logaddexp(t, -t) - log(2)` and satisfies Fenchel-Young equality to `0.0` (`<1e-12`).
   - Dual lower bound $\text{LB}(q, r) = -\frac{1}{2}\|q\|_2^2 - \operatorname{Re}\langle b, q\rangle - g^*(-B^H q - T^H r)$ for $\|r_{\text{pixel}}\|_2 \le 1$ has exact signs and dual feasibility checks.
   - Safeguarded Newton-bisection solves the radial prox equation $m + \tau_i w_i \operatorname{atanh}(m) = |v_i|$ on $[0, 1)$ to `1.90e-15` max FOC residual without TV or entropy smoothing.
4. **Open vs. Closed Domain Equivalence**:
   - While `physical(P)` evaluates on the open interior $|P_i| < I_i$ ($|x_i| < 1$), $h(m)$ is continuous on $[0, 1]$ with $h'(m) \to +\infty$ as $m \to 1^-$. Thus $\inf_{|x|<1} J(Ix) = \min_{|x|\le 1} J(Ix)$, and every proximal iterate remains strictly inside $|x_i| < 1$ (`max_m = 0.7113` on Apr 10).

## 3. Diagonal Preconditioning & Contraction Proof
- Column and row $\ell_1$ accumulators $c_i = \sum_j |K_{ji}|$ and $d_j = \sum_i |K_{ji}|$ accurately include interior and boundary stencils (`col = 2*I` plus interior left/up neighbors; `rx, ry = I` plus interior right/down neighbors).
- Using the **same** $\sigma_{\text{TV},\text{pixel}} = \theta / \max(d_{x,\text{pixel}}, d_{y,\text{pixel}})$ on both $(dx, dy)$ components at each pixel ensures:
  1. The $\Sigma^{-1}$-metric dual proximal step on the coupled isotropic TV ball remains the standard Euclidean 2-channel unit-ball projection (`project(r + st * s.T(xbar))`).
  2. The complex Schur bound $\|\Sigma^{1/2} K \mathrm{T}^{1/2}\|_2 \le \sqrt{\max_j(\sigma_j d_j)\max_i(\tau_i c_i)} \le \theta = 0.99 < 1$ holds strictly.
- On Apr 10 fixed-$I$/fixed-$D$ data, diagonal PDHG reaches $J = 3.27039786797$, $\text{LB} = 3.27039590460$, gap $= 1.963 \times 10^{-6} < 10^{-5}$ in 1,000 steps (`6.43s`), outperforming 30k-step `mcv` L-BFGS-B (`3.276106`, `393s`).

## 4. Required Caveats & Production Guardrails
1. **Conditional Floating-Point Certificate**: The duality gap is a valid **conditional numerical certificate** in IEEE-754 double precision (subject to standard floating-point roundoff), not rigorous interval arithmetic.
2. **Subproblem Scope**: Convexity and the primal-dual gap certificate apply strictly to the **fixed-$I$, fixed-$D$ subproblem**, not global joint optimality over $(P, D)$.
3. **Outer $D$-Loop Invalidation Rule**: Updating $D$ changes the measurement operator $B(D)$ and data $b(D)$, defining a new optimization problem. Any prior lower bound LB is invalid once $D$ changes; production code must rebuild the operator, diagonal steps, and certificate at each outer $D$ step.
4. **Numerical Criterion Only**: The target gap $< 10^{-5}$ is strictly a numerical solver stopping/certification tolerance and must not be claimed as a scientific or image-space uncertainty threshold.

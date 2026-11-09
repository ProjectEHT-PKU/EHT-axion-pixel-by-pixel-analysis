# Glossary

## Purpose

This page defines the symbols and terms used in the documentation and in the code, with
their units and the YAML key or code name where one exists, so that the other pages can
use them without repeating definitions.

## Details

### Physical quantities

| Term | Definition |
|---|---|
| $t_g$ | gravitational time $GM/c^3$; 8.893 h for $M = 6.5\times10^9\,M_\odot$ (`PhysicalScale.tg_hours`) |
| $r_g$ | gravitational radius $GM/c^2$; 3.819 µas at 16.8 Mpc (`PhysicalScale.angular_rg_uas`) |
| $\alpha$ | gravitational fine-structure constant $GM\mu/(\hbar c)$ of an axion of mass $\mu$ (`shared.alphas`); $\mu c^2 = \alpha \times 2.06\times10^{-20}$ eV for the M87\* mass |
| 211 state | axion cloud bound state with $n = 2$, $l = 1$, $m = 1$ |
| $\Omega(\alpha)$ | dimensionless frequency $M\omega = \alpha(1-\alpha^2/(2n^2))$ with $n$ = `shared.principal_n` (`physics.omega`); $\omega = \Omega/t_g$ |
| $c$ (coupling) | dimensionless model coupling that scales the axion-induced EVPA signal; called `g` and `couplings` in `likelihood.py`. The code makes no conversion to an axion-photon coupling in GeV$^{-1}$ |
| $c_{95}$ | posterior upper limit of $c$ at `analysis.inference.credibility` (0.95) |
| $b(\alpha)$ | cloud amplitude scaling `analysis.cloud.coefficient` $\times(\alpha/$`pivot_alpha`$)^{\text{exponent}}$ |
| $N_\mathrm{ref}$ | RAPTOR axion normalization `raptor.response.reference_norm`; templates are per this norm |
| axion phase | RAPTOR phase index $k$ with phase $k\pi/8$, $k = 0,\dots,15$; in the likelihood, the unknown phase $\phi_0$ |

### Polarimetry

| Term | Definition |
|---|---|
| Stokes I, Q, U, V | intensity and polarization; RAPTOR cubes in Jy/sr, EHT FITS in the units of their `BUNIT` |
| $P$ | linear polarization $Q + iU$; $\lvert P\rvert = \sqrt{Q^2+U^2}$ |
| EVPA | electric-vector position angle $\chi = \tfrac12\arg(Q+iU)$, in radians, defined modulo $\pi$ (`evpa.finite_evpa`) |
| $\mathrm{wrap}_\pi$ | difference wrapped to $[-\pi/2, \pi/2)$ (`evpa.circular_difference`) |
| resultant length | $R = \lvert\langle e^{2i\chi}\rangle\rvert$ over draws; 1 for identical EVPAs |
| $m_\mathrm{net}$, $m_\mathrm{res}$ | net and resolved fractional polarization $\lvert\sum P\rvert/\sum I$ and $\sum\lvert P\rvert/\sum I$ |
| $\beta_2$ | $\sum P\,e^{-2i\varphi}/\sum I$, $\varphi$ the azimuth east of north (whole field) |

### EHT imaging

| Term | Definition |
|---|---|
| day 95, 96, 100, 101 | day of year 2017: April 5, 6, 10, 11 |
| day pair | (95, 96) and (100, 101); the data are the differences 96−95 and 101−100 |
| D terms | complex instrumental polarization leakage $D_R$, $D_L$ per station (fractional) |
| solve sites, fixed sites | AZ, LM, PV (D terms solved here) and AA, AP, SM, JC (pre-corrected in the HOPS data) |
| fsys | fractional systematic error added to RL and LR errors (`eht_imaging.fsys`) |
| fixed I | the total-intensity image of the I stage; all later steps solve only for $P$ |
| fiducial | the P/D alternation result of `eht_reconstruct.py --stage fiducial` |
| tightD | the refined fiducial D terms of `eht_refine.py` |
| display19 | images blurred with the 19 µas FWHM Gaussian `shared.image.blur_fwhm_uas` |
| draw | one D-term sample and its certified P image; 1000 per day |
| gap certificate | best primal objective minus best dual lower bound of the exact P solver, at most `gap_atol` |

### Analysis

| Term | Definition |
|---|---|
| state | index along axis 1 of a RAPTOR cube: 0 axion-off, $k+1$ axion phase $k\pi/8$ |
| snapshot | one GRMHD time; 141 snapshots, 1 $t_g$ apart |
| response | per-snapshot first-harmonic fit $\bar\chi + A\cos(\psi+\phi)$ over the 16 phases |
| template | mean over snapshots of the response amplitude $A$ and phase $\phi$ |
| 1D, 2D | aperture-sector observable (360 one-degree sectors) and pixel observable (64 × 64 grid) |
| feature | one selected sector (1D) or pixel (2D) used in the likelihood |
| GRMHD background | $\sigma_\mathrm{bg}$, EVPA scatter of lagged axion-off frames on the selected features |
| zero-data control | the same likelihood with day differences set to zero (`null-comparison`) |

## Related pages

- [Documentation index](index.md)
- [Overview](overview.md)
- [Analysis](analysis.md)
- [EHT imaging](eht-imaging.md)

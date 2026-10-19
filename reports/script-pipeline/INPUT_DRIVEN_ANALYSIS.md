# Input-driven analysis audit and rerun

## Scope
User requested a fresh analysis of the new merged RAPTOR products and EHT
reconstruction ensembles, explicitly checking that derived scales no longer
silently use historical numbers. This is NOT a new RAPTOR or EHT-imaging run.
Fundamental physical constants remain versioned/documented; independent model
choices are explicit inputs, not incorrectly described as derivable from M/D.

## New profile
`configs/analysis.yaml` keeps the completed producer identity unchanged, including
M6.5e9/D16.8, camera46.2261Rg, referenceNorm0.1, mixed cloudspins and principalN2.
New sky-coordinate input is128uas across64pixels, verified from original four-day
FITS headers and trusted pickle shapes. Actual observation frequency in those
headers is226191789062.5Hz;230GHz is the band label, NOT exact frequency equality
to the RAPTOR230GHz. No frequency correction or reconstruction is implied.

`analysis_geometry()` derives:
- theta_g=3.818993252351109uas/Rg;
- tg=8.893247544352288hours, inverse=0.11244486280324628/hour;
- target128uas =>33.51668660875444Rg, pixel2uas;
- aperture outerradius=halfFOV fraction1 =>16.75834330437722Rg;
- matchingkernelFWHM19.503299817566703uas =>sigma2.16871Rg.

Kernel value explicitly preserves the prior sigma at current M/D; it is NOT an
unapproved switch to20uas. This angular FWHM is now an independent input and
converts if M/D changes. It is not an effective-PSF calibration. The previous
32.6868Rg target was124.83066864uas, inconsistent with128uasEHT; this rerun fixes
that sky-grid mismatch. Old research/production profiles and outputs are intact.

Background cadence is now derived from generated GRMHD footer times rather
than fixed10*label_step. Actualspacing1tg yields8.893247544352288h andlag3
26.679742633056865h, not an exact24h; no interpolation/sqrt2 change introduced.

Template referenceNorm explicitly flows fromproducer profile+verifiedgenerated
identity toLimitOptions/evpa_templates; signal and propagated sigma scale by its
inverse. Cloud prescription coefficient0.03/pivot0.05/exponent1.75 is now explicit
config/provenance; it is a retained model assumption, not computed from M/D.

EHT shape/pixel/FOV contract checked against raw arrays before computations.
Rank/final per-day shifts, angle_roll, maskrectangles/columns, TopSet1572,
intensityweight, aperture10deg, lag3, ridge0, NULLcentres and other existing
choices retained and explicitly recorded; no covert TopSet/statistics revision.
Four originalFITS samples and input-geometry evidence are declared EHT ancestors.

## Compatibility boundary
TIME_CONVERSION0.112411, AXION_NORM0.1 and oldexpertFOV/cadence defaults remain
only as legacy publicAPI defaults, preserving historicalregressions. New angular
workflow requires explicit cloud+EHT contracts and supplies all derived/science
parameters. Thus 'no constants' means no hidden historicalscale fallback on the
new main execution path, NOT deletion ofc/pi/GM or all supportedlegacydefaults.

RAPTOR effectivebuild uses inherited roundedCGS constants. Analysis uses modern
IAU nominal solarGM/exactc/parsec. ProducerGM/Rg/tg differ about0.02518%, netangular
scale differs about0.01483%. These origins remain distinct; no existingbuild/data
was re-signed or relabeled to pretend the constants were identical.

## Execution (pending at document creation)
Inputs:
- /home/canoziia/GRMHD_DATA/production/m65-211-merged/data
- /home/canoziia/GRMHD_DATA/EHT_DATA/analysis-inputs (symlinkview to backed-uporiginal
  fourdaypickles plus headergeometry evidence anddeclaredoriginalFITS)
New outputroot planned:
`/home/canoziia/GRMHD_DATA/analysis/m65-211-input-derived-v1`.

Tests and exact launch/completion evidence will be recorded separately. Do not
interpret this document as a completed science result. No prior scienceoutputs
are overwritten and no backgroundraw-generationqueue is relaunched.

Remaining interpretation caveats (not silently fixed): TopSetdefinition/coverage,
A/phiindependence,scalarbackgroundfor1D,discretelag/exposuretimecovariance,
realradialenvelope/spectrumapproximation,fluxcalibration and screenconventions.

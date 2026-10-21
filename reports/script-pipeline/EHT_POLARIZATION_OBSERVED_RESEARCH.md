# EHT official-method audit, raw-FITS diagnostics, and observed-data research

## Status and what was actually done
User explicitly authorized a separate observed-data study and a persistent goal.
No NULL-only inference API was changed and no historical result was overwritten.
No RAPTOR C generation or full UVFITS imaging was rerun. Local backed-up original
19200 FITS (fourdays x4800) were processed and verified against declared metadata.
All six observed sensitivity cases completed successfully, `failures=[]`, EXIT0.
The matched19uas RAPTOR postprocessing/template run also completed EXIT0.

This study reproduces documented *types of processing and diagnostics*, not the
EHT Collaboration's five-pipeline imaging/calibration survey or official images.
We have one historical pipeline with a4800-hyperparameter grid, not1000calibration
realizations and not five independently validated fiducials. Resulting c95 values
are **conditional working-covariance research limits**, not calibrated exclusions.

## Primary sources and corrected interpretation of Top Set
See `eht-polarization-evidence/OFFICIAL_METHOD.md` and `LOCAL_PROVENANCE.md`.
PaperVII: https://arxiv.org/abs/2105.01169 DOI10.3847/2041-8213/abe71d.
PaperIV: https://arxiv.org/abs/1906.11241 DOI10.3847/2041-8213/ab0e85.
PDFs were rendered and pages read as images, not OCR/pdf2text. Key renderedpages
are additionally preserved under PHYS/reports/eht-polarization-research/paper-pages.

- VII Fig6: four days by five methods; Fig7: average of five equally blurred
  fiducials **on each day**, averaging Stokes separately before P/EVPA.
- VII Fig8:1000Dterm-perturbed reconstructions or posterior samples, not a
  hyperparametergrid and not EVPA-residual-center selection.
- VII§5.1:eht-imaging comparisonkernel19uasFWHM, not automatically a measuredPSF.
- VII§5.2:I-only correlation centering and coherent10degree angular Stokes
  averages at1degree steps; lowP/north andApr10coverage especially uncertain.
- IV TopSet1572 arose from37500StokesIparametercombinations, accepted by real
  closurefit AND synthetic recovery.1572 is not a polarized sample-count rule.
- VIIpolarizedfiducial selection uses six polarizedsynthetics andDterm recovery;
  it is not IV'sTopSet copied into this grid.

Our old1572 is selected by EVPA centrality. It overlaps only680/1572 with a
fourday mean-I datafit ranking of the same size. This establishes a selection
rule mismatch, not a complete causal explanation of every daydifference.
Originalpol_sorted is a different output branch, NOT merely reorderedrawFITS;
day100metadata has6742rows/1942duplicates. It was not mixed into this study.
The complete rawgrid, seven-parameteridentities, fivefinitefitmetrics are used.

## New reproducible components
Optional Pixi environment `eht-research` adds Astropy; default workflow unchanged.
- `eht_axion/research/polarimetry.py`: nativeFITS/WCS/Stokes validation, skyEVPA,
  coherentmapquantities, convolution andI-onlyrelative registration utilities.
- `scripts/research_eht_polarization.py`: streamoriginalFITS tofloat64mmap;
  selection/blur/alignment/meanmaps/realdaydiff artifacts.
- `scripts/research_eht_diagnostics.py`: fourday/aperture distributions,
  directionconcentration andsky-integratedpolarization diagnostics.
- `eht_axion/research/observed.py`: separateactualdaydifference likelihood,
  fixedobs covariance, commonphasejointproduct; originalNULL APIuntouched.
- `scripts/research_observed_limits.py`: verifiesdeclaredstudyandtemplateinputs,
  generatedphysicalprofile,preprocess/aperture/normalization/alphagrid, produces
  observed andmatchedzero-center comparisons.
- `scripts/summarize_eht_research.py`: declaredresultaggregation/plots/CSV.
- Configs `eht-polarization-study.yaml` and `eht-observed-research.yaml`.
Docs: `docs/eht-polarization-research.md`.

Tests:1456passed215.45s inoptionalAstropyenvironment,Ruffpassed; following
rawalignmententrysupport regression116passed5.27s+Ruffpassed. Independentreview
foundandfixed3provenancegaps: actualtemplate/configidentity,declaredfilemembership,
observedscriptbytesasinput. Coordinate/dayorderreview foundno blockingtranspose,
spin2constantoffsetorlatersignerror. This is implementation evidence, notphysics
orlikelihoodcalibrationcertification.

## Exact processing contract and caveats
NativeFITS(y,x,IQUV), eastdecreasing/northincreasing,64²/128uas/2uas-pixel,
frequency226191789062.5Hz, Jy/pixel; noQ/U signflips in savedskymaps. EVPAEofN.
FITS often useOBSRA/OBSDEC instead ofCRVAL; loader supports this actualwriter.
HDU3isV, NOT EVPA despite oldpicklevariable `theta_data`.
Additional19uas convolutionappliedsameIQUV; original5uasrestoration retained.
This is not a certifiedeffectivebeam or a full reproduction of officialconsensus
matching. Croppeddisplay120uas, absolute10uas ticklength=30microJy/uas²P,
colorfraction0..0.3. DisplayI>=.1peak AND P>=.2peak isexplicitstricterchoice;
papercaption's proseAND ambiguity not falselyresolvedasofficialplotcode.

RegisteredtoonefixedlocalApr11bestmean-Iχ²image, integertranslationmax8pixels,
notofficialconsensusastrometry. AllStokestranslatedtogether,zero-paddinglosses
recorded. Rawvsalignedoutputsavailable. Somecandidates hitsearchboundary:shared
Apr11=15.49%,Apr6=3.97%,Apr5=.05%,Apr10=0%. Registration is not solved merely
because a correlation wasmaximized.
NativearraytransposeplusconstantQ/Uminus maps tooldtemplatebasis; fixed90deg
EVPAoffset cancels in differences/covariance. Aperturequadrature androll-90
match existingRAPTORfitspatialconvention, not a newly chosen sign toalignresults.

Artifact metadata note: full-study `artifact_schema.aligned_stokes` text retains
an erroneousfloat32 label fromdraft; actualNPYheaders/manifest schemas/arraysare
float64. This descriptive typo is documented, **not re-signed** after publication.
Its executing source bytes are preserved so liveancestry remains valid.

## Selection comparisons
All conditions below are local diagnostics, NOT officialTopSet:
|Variant|Apr5|Apr6|Apr10|Apr11|
|---|---:|---:|---:|---:|
|all4800|4800|4800|4800|4800|
|legacyEVPAcenter|1572|1572|1572|1572|
|each storedvis/cphase/logcamp/pvis/m <2, perday|2608|2124|2148|2228|
|sameidentitypassesallfourdays|1989|1989|1989|1989|

Storedfitmetrics are optimizer-stage values BEFORE final5uasrestoration; this is
not a fresh final-FITS/UVFITSfit anddoesnot certifyDterms/selfcal. Exacthistorical
customehtimlibrary/gains are missing, so equalweights don't guarantee sameobjective.

## Daydifference diagnostics: drift is not wholly a TopSet artifact
Forfourdayshared1989ensemble, registeredmeanStokesmaps give median Δχ:
Apr6-Apr5=-4.43deg;Apr11-Apr10=+6.97deg (brightintersection499/538pixels).
Changingselection doesn't removeit:earlymedians-4.97/-4.60/-4.74/-4.43deg for
all/legacy/perday/shared;latemedians+5.34/+7.54/+6.74/+6.97deg.
These are differentmasks, so they are not an exactfixed-pixelcausalcontrast.

1DconcentrationR>=.7 latepairretains97/202/127/137of360angularbins for
all/legacy/perday/shared. ±90degverticalplotjumps areEVPAmod180seams, notphysical
jumps. LowP/northbroadness,selectionandregistrationallmatter. Brightsouthwest
PA~180..270deg retains differences ofseveral to~10deg acrossvariants;not a rigid
rotation ofthewholering andnot anaxiondetection.

RegisteredsharedmeanStokes Apr5/6/10/11:
netfraction2.40/2.07/2.73/2.66%;resolvedfraction8.84/7.30/6.05/6.53%;
|beta2|.0368/.0410/.0158/.0346;argbeta2=-169.66/-161.27/-152.51/-160.96deg.
Thesearelocalensemblemeanquantities, notofficialDterm/posteriorintervals;
Apr10differsespecially. Noagreementthresholdusedtotunethesevalues.

## Observed posterior study and exact comparison boundary
Two actual paircenters,Apr6-Apr5 andApr11-Apr10, areused, NOTsettozero. Same
covariance/mask/parameters thenusedforamatchedzero-centercomparison. Perday
hyperparameterpopulationcovariances onlocalcircularbranches are summedassuming
zerocrossdaycovariance. Noautomaticridge/clipping;concentration>=.7;1Dnonoverlap
10degspacingoffset4;2D5pixelstrideplusbrightnessintersection. Featurecounts vary
withselection and are deliberatelyreported. RMLimages are NOT independentnoise
realizations;thislikelihood isaconditionalworkingsensitivityonly.

Atalpha.30, correlatedobservationalworkingcovariance (noresponsemodelvariance):
|Selection|1Dfeatures|observedc95|samecovNULLc95|2Dfeatures|observedc95|samecovNULLc95|
|---|---:|---:|---:|---:|---:|---:|
|shareddatafit|10|4.1123|2.1537|13|4.9746|2.5366|
|perdaydatafit|9|4.0859|2.2720|13|5.0108|2.5389|
|all4800|7|4.8471|3.2724|12|4.1581|3.5653|
|legacyEVPA1572|16|3.8460|1.0599|16|4.4456|0.7353|

Raw1D(noextraIregistration)shared:8features,c95obs4.9168/null2.7705;
rawlegacy:16features,c95obs4.4371/null1.1487. Thiscomparisonalsochangesfeature
selection via concentration, so don'tattribute allvariationto translation alone.
Sharedactualdataα.3jointchi²NULL=5.2294,best=2.3439(1D20datacomponents),
NULL=6.4060,best=3.4729(2D26components). Improvements~2.9arenotcalibrateddetection
significances. Legacygiveslarger apparentimprovement butdifferentfeatures and
centrality-selectedcovariance;cannotclaimdiscoveryorselectitbecausemorestringent.
Observed95%boundsarenotinlastgridbin;lowalpha~.1canremainpriorboundsensitive.
No formalphase/coupling-grid/branch/priorcoverageconvergenceisclaimed.

Do NOTcomparethesenumbersdirectlytoearlierNULLsixcurvesasifonlycenterschanged:
newobsstudyalsodiffersinkernel,registration,selection,featuremaskandcovariance.
Only dashedNULLcurves **inside each researchcase** arecontrolledcentercomparisons.
NoGRMHDresponsevariance/logdetdebateisresolvedbythisfixedcovarianceexperiment.

## Locations and final interpretation
14TBoutputsunder `/home/canoziia/GRMHD_DATA/analysis/`:
- `eht-polarization-study-v1`:104declaredoutputs,IQUVmmap/means/differences;
- `eht-polarization-diagnostics-v1`:9declaredoutputs,aperture/histogramfigures;
- `eht-observed-research-templates-v1`:19uasmatchingRAPTORpostprocess/templates;
- `eht-observed-{datafit_shared,datafit_perday,all4800,legacy_evpa}-v1`;
- `eht-observed-{datafit_shared,legacy_evpa}-raw-1d-v1`;
- `eht-observed-research-summary-v1`:aggregationstage (completionseeitsmanifest).

**Answer:** we can reconstruct the localFITSintoofficial-inspiredfourday
polarizationmaps and quantify selection/alignmenteffects withoutnewimaging.
Our oldTopSetdefinitionisnot theofficialdefinitionandcanartificially narrow
spread. Butrealdaydifferencesremainafterchangingselection;theirphysicalorigin
cannotbeassignedfromtheseimagesalone. Underexplicitworkingcovariance,actual
centersweakenlimitsrelativetomatchedNULL. Thisisusefuldiagnosticevidence,NOTa
newcertifiedaxionexclusion. ReliableobservedconstraintsstillneedfinalFITSdatafit,
calibration/Dtermrealizations orsyntheticclosure/noisevalidation,representative
image/modelensemble andstablebeam/registration/systematicuncertainties.

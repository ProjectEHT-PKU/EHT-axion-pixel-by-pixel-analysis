# Completed: verified4000-image copyback → EHT errors → real-day-difference limits

## Delivered and verified

The user-authorized post-run Goal has executed, not merely produced a plan.
Remote4×1000 fixed-D images completed EXIT0. Entire result/control/input/source
snapshot copied locally; all4015attempts retained (4000certified,15uncertified
attempts recovered byretry). Verified count **24,164 files**, total
1,241,843,797verifiedbytes. Output/input/runtimecode SHA checks passed; daily
membership exactly1000each. Remoteoriginals retained. Sourcebinaryenvironment
excluded as reproducible; dependency/version/sourcehash receipts retained.

Localdataset:
`/home/canoziia/GRMHD_DATA/analysis/eht-uncertainty-remote-v1/`
with `COPY_VERIFIED.json`. Batch identity
`06a8f63baa95b7ab959c951e31ea45e6b6933a3db06455e54a65fe8c7a328af0`.

## Completed new stages (all EXIT0)

All paths below relative to `/home/canoziia/GRMHD_DATA/analysis/`:

|Stage|Root|Products|
|---|---|---|
|1000draw/day summaries|eht-mc-summary-v1|92declaredoutputs:raw/display19 native/I-centroid maps,perdrawEVPAs,sectors,covariances,means,globalmetrics|
|MC diagnostic figures|eht-mc-figures-v1|5figures+summary|
|120uasRAPTORprocessing|eht-120uas-observed-templates-v1|8αresponsepreprocess/fit,1D/2Dtemplates,backgroundstate0andunmaskeddiagnostics|
|Observed+matchedNULLprimary|eht-120uas-observed-limits-v1|96records=6cases×8α×2centers,posteriors,CDFs,selectedsamples,Cobs/Cbgdiagnostics|
|Nativecenteringcontrol|eht-120uas-observed-limits-native-v1|same96records,noI-centroidrelocation|
|Densegrids|eht-120uas-observed-limits-gridcheck-v1|phase128/coupling240 vs64/120|
|Upperpriorsensitivity|eht-120uas-observed-limits-priorcheck-v1|cmax10000 vs1000,phase128/coupling260|
|FiniteMCcovariancecheck|eht-120uas-observed-covcheck-v1|first/last500drawcovariances,fixedcenters/selection/background|
|Aggregate|eht-120uas-results-summary-v1|limits/fitsplots,tables,sourceSHAreceipt snapshots,diagnostics|

No rawRT regeneration. Original128uasproducts andcanonicalNULLcode/results were
not edited or rebound. Newobservedresearchlikelihoodusesactual96−95,101−100
differences, commonphasejointproduct thenphaseaverage, noexp-externalsumrule.

## Geometry, observables and covariance

64²/120uas=1.875uas/pixel. Actualobservingfrequency227.070703125GHz vsnominalRT230GHz
recorded, not relabeled. RAPTORoriginalFOV46.2261Rg androtation108 retained,
reprocessed to120uas. Commonadditional19uasIQUkernel is an explicitlydeclared
matching approximation, not certifiedeffectivePSF. FixedI restoration details
and original reconstruction provenance retained.

NativeFITS(y=Nincreasing,x=Edecreasing) is transposedspatially to template(x,y)
withinternalaxes(-E,N). NoadditionalRAPTORrotation orper-dayEVPAalignment applied.
I-onlycentroidshifts<1pixel; allStokeschannels/draws/resolutions shareeachshift;
edgeI losses~.33–.40% fordisplay19 products. Native sensitivity kept separately.

Inference1D recomputesall1000draws using the SAME bilinearCSR aperture operator
asRAPTOR, includingtheexistinginternal−90 conversion exactlyonce. Nativepixel
sectorcurves remainpaper-style diagnostics, not likelihoodsamples.2Dusesdirect
spatialtranspose. Brightness/resultantmasks selectedwithoutlookingatdaydifference
magnitudes;23angularfeatures and20pixelsinprimary. No covariance ridge, eigenvalue
repair, phase recentering or successful-subset averaging. Actualcircularmeans,
centeredpopulationlocalbranchcovariance(ddof0), NOT covarianceofmean or÷1000.
Independentconditionaldays assumed: pairC=Cdaylate+Cdayearly, no fake matching
ofindependentlygenerateddrawIDs. Crossdaycalibrationcovariance is not known.

SelectedGRMHD pair-differencescalarσ=2.15354deg1D/2.10327deg2D, estimatedafterfinal
featureselection. This is fixedbackgroundvariance, notcanonicalc-dependentresponse
uncertainty. Actualsimulationlag26.6797h vsobservational24h isrecorded benchmark
approximation; no silent temporalrescale. A/phase response uncertainty andits
c²variance are excluded inthis observed engine, not silently claimed included.

## Actual α=.30 results: model coupling c95

|Mode|1Dobserved|1DmatchedNULL|2Dobserved|2DmatchedNULL|
|---|---:|---:|---:|---:|
|std(diagonalDcov)|1.40844|.439973|.900964|.435339|
|cor(fullDcov)|1.11550|.415391|1.49221|.264602|
|cor+fixedbackground|1.79792|.633656|1.38277|.514972|

Theseare95%conditionalBayesianuppercrediblelimitsonmodelc, NOTGeVconversion,
frequentistexclusion,detectionormedianexpectedlimit. Curves/tablescontainall8α
(.10,.15,.20,.25,.30,.35,.40,.44), andidenticalCzero-centerNULLcontrols.

## Checks that change interpretation

- Densegridmaxchange.004435dex≈1.03%. Priorextensioncomparedtodensegridmaxchange
  .000316647dex≈.073%; savedCDFsmonotoneandnormalized,uppermostbinmass~1e-15.
  Lowerloguniformpriorremainsconditionalassumption.
- SplittingMCcovarianceintofirst/last500draws atfixedobservationcenters/features
  changesboundbyupto.10610dex≈27.7%, especially ill-conditionedcor1D. Gridaccuracy
  isnot1%totalphysicalprecision. Raw1Dpairconditionnumbers~6.62e7/9.63e5,
  2D~9.82e3/9.64e3. Covarianceestimationuncertainty not marginalized.
- Nativevscentroidmaxlimitdifference.569847dex (up to~3.71ratio). HereNcounts
  happenboth23/20; equalcountsdoNOTproveidenticalselection/data. Treatasjoint
  registration/selection/model sensitivity, not isolatedcausalshiftmeasurement.
- Observedfeaturesare>3marginalσfrom±90degseams; reportingseamdistanceisnotproof
  allpriorpredictiveGaussianresidualsvalid. No data-dependentseamexclusions.
- Atα.30,corD-onlyχ² is158.8→153.0for46components1D;196.0→171.4for40components2D.
  LargeabsoluteresidualsmeanD-onlycovariancemodelisstrained, notaxiondetection.
  Withfixedbackground:25.0→19.7(1D),44.9→38.1(2D). Improvedχ² doesnotbyitself
  independentlycalibratethebenchmarknoise. Noformalphysicalexclusionclaim.

Whole-map/sectorsMCdiagnostics arein `eht-120uas-evidence/MC_DIAGNOSTICS.md`.
Apr10netm~.995%withbroadintegratedEVPAstd26.55deg; brightSWsectors stillhave
conditionalσmedian~3.37deg andlatepairseveraldegreevariation. Netvectorcancellation
isnotwhole-imageunconstrainedness. Otherdayssectorσ andactualdaydifferencesare
reportedquantitativelywithoutconvertingmaximumfeaturesintodetections.

## Testing, source integrity, remaining limits

Finalfulloptionalregression1573passed236.08s after optionaldependency/path testfixes;
newnumericsreuse/aggregation/covariancetests6passed. Targetedensemble/inference/figures33pass.
FinalSHAverificationof10completedresult/templateoutputstagesalsocompletedEXIT0;
receipt `eht-120uas-final-output-verification.json` copiedintoevidence.
FullnewsourceRuffpassesexceptknownone111-charE501instartrunorchestrator,whosebytes
werepreservedafteractualexecution ratherthanrewritingitsliveancestry. Theexecuted
covariancecheckhasone112-charE501likewisepreserved. Noglobalsuppression.

Initialregressionfailureswereoptional-ehtimdependencytestwithoutskip andastale
syntheticbackgroundpath;tests corrected, notsciencealgorithmchanged. Gridcheckfirst
invocationwasstoppedbeforeoutputuponunnecessaryrawrehash; successfulnewversionuses
byte-identicalsource-receiptsnapshotsunderexplicitnames(notfakecanonicalmanifests)
andchecksconsumedresultoutputSHAs. Mainandnativestagesverifiedfullrawlineageseparately.
Finalaggregateverifiesallcompletedconsumedresultoutputs; nosourcesresigned.

Sourcescripts/data/FITS/publicationreportsretained; noautomaticpush. Figurelayout
infirstdaydifferenceplot hasoverlappingcolorbar, documenteddisplayissue; summary
limits/fitfiguresareseparatecompletedproducts. Notusedtoalteranyinference.

**Delivered:** reproducibleconditionalpolarizationerrorandreal-day-differenceconstraint
pipelinewithverifiedcopy,matching120uasresponses,backgrounds,controlsandchecks.
**Notasserted:** all-noisecalibrationcoverage,officialfive-pipelinemean,formalaxion
exclusion,completejointcrossdayDposteriororexacteffectivePSF. Thosearemodel/science
limits, not unexecutedparts ofthis authorizedworkinganalysis.

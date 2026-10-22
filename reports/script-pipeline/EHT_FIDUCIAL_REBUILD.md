# Apr11 UVFITS reconstruction with current official ehtim

## Actual completion and version
PyPI and GitHub latest formal release verified as ehtim1.3.2 (2026-06-16),
not historical1.2.7. Official1.3.2 is pinned in optional `eht-rebuild` Pixi env;
legacy1.2.7/customfiles remain untouched. RuntimePython3.12, numpy2.5.4,
scipy1.18.1, astropy7.2.0. All reconstruction described below was actually run,
not merely planned. No oldFITS/NULLproducts were overwritten.

OfficialHOPSApr11low-bandUVFITSdownloadedfrompinnedEHT2023-D01-01commit:
`b71399714fdfe13579c741075d757b4ce65fd3d1`.
SHA2565cf936dce926b5b13a0f39065096da1f836ab52ca95fdfffef581057029e8d5b,
524160bytes, actualrf227070703125Hz,5877rowsbefore22scanaveraging.
AA/AP/JC/SMzero-baselineDtermsalreadycorrected;AZ/LM/PVremain tosolve.
JCnotpresentthisday. AllabsoluteEVPA/parallacticcorrections arealreadyapplied.
`dcal=True`UVFITSmetadata isNOTproofallstationsarecorrected.

OfficialpolarizedimageFITSorcompletePaperVIIexecutablepipelinewasnotretrieved.
Onlyaone-methodofficial-inspiredreconstructionisclaimed,notbitwiseofficialor
five-methodconsensus. Accessresearchrecordandfeedtableevidencearepreserved
underreports/script-pipeline/eht-fiducial-evidence/.

## Real staged attempts, not interchangeable successes
1. Engineeringpilotv1failedinitialTVgradient(nonfinite);retainedlog, nofalseimage.
   AddedofficialI-script1e-3offcenterGaussianinitialization andfluxrenormalization,
   preflightfiniteobjective/gradient gate;noepsilonTV orupstreamcodepatch.
2. Pilotv2twoP/Drounds completedEXIT0, partialnotconverged.
3. SimplifiedI+50P/DroundscompletedEXIT0,butIampchi²~18 unacceptable,notafiducial.
4. SwitchedtoactualofficialPaperIVreleasedIscript,withoutrewritingitsstages.
   runpyexecutesexactsourcewithonly nfft->directDTFT backendadaptation;IERSoffline.
   Releasedfullselfcal/convergence/regularization/short-baselinescaling retained.
   Outputactualfinalselfcalchi²:vis.43043,amp.38329,cphase.80662,logcamp.55924.
5. FirstP/DfromIcheckpointcarriedI2%noiseinflation;preservedasdiagnostic,NOT
   thermal+0.5%polarimetricerror. Thiswasexplicitlyidentifiedandnotpresentedfinal.
6. CleanPinput:applyofficialIfinalgainstabletooriginalscanaveragedobsBEFORE
   reverse-taper/I2%noise. Stationtablemustreorderbyidentitybeforeapplycal.
   Firstexportfailedtableordermismatch;fixednewrunv2EXIT0,neverresignedoldrun.
7. CleanP/D50roundsEXIT0:χ²pvis2.2653,m2.3291;every-roundblurrestartprevented
   smallmapchanges;noengineeringorscientificconvergenceclaimed.
8. Warmstartcomparisonusesblurrestartsonlyrounds2,12,22,32,42 ('severalrestarts'
   ratherthaneveryDupdate),allotherstepsfixed,50roundsEXIT0. Finalpvis2.31535,
   m2.37255,finalQUrelativechange.0030653,Dchange.00019065. Stillaboveexplicit
   .001/.0001engineeringthresholds;scientifically_converged=false.

## Exact implemented P/D choices and limitations
FixedIderivedfromofficialreleasedIcheckpoint128uas,regridded120uas64²andcompact
fluxrenormalized.6Jy;thisadaptationisrecorded,notclaimedidenticalPaperVIIIimage.
Pvis/mweights1/1,HW100,PTV1,fsys.005,100iterationsperPoptimizer,
AZ/LM/PVfullJonesleakagesolvewithtol.1;zero-correctedsitesheldzero.
EachroundDsolveusesfreshgain-correctedoriginalbase,notpreviouscorrecteddata,
soDtermsarenotcumulativelyapplied. I isheldfixedduringPstep.

FeedrotationcoefficientsfromPaperVIITable1andpinnedupstreamEHT2017_template.txt,
sha256a5fcb2ae7e991d2fdd259187d37bf0253a5ff2212f1f71084afeff1d9e0c85e7.
fr_offusesdegrees(SMA45),notradians. Correctedfrcaldata stillneedsfeedcoeffs
forleakage'sexp(±2iphi). UVFITSloaderdefaultzerosarenotusedasphysicalmounts.

Public1.3.2Imagerdoesnotapplysystematic_noisetoPvis/m. Wrapperadds
(.005*absVI)^2topsigma²bysplittingqsigma/usigmavarianceequally,onceperfreshbase,
forbothP/D. Thisisexplicitworkingmarginalerrorpropagation;Stokes/circularchanges
andJonesstepsdon'tcarrycompletecross-handcovariance. Noassertionoffullnoise
covariancecorrectness. Dprior.1isretainedimplementationchoice,notPaperVIItraining
score'sdtol.05. Don'tconfusetheseparameters.

BlurQ/UthenrestorefixedIcanbreakP<I;restartprojectsviolatingpixelsinside
mcvphysicaldomainandrecordscounts. Thatprojectionisininitialization,notdata
clipping. Eachoptimizeractualnit/statusstored;I/Pmax100,Dupstreammax1000,
selfcalupstreammax10000. EXIT0doesnotimplyeachoptimizercertifiedconverged.

Finalrawimageandadditional19uasI/Q/UdisplayFITSsavedexplicitlywithbothblur
arguments;nosilentI-onlyblur. FinalwarmstartimageIflux.6Jy,netfraction1.1123%,
resolvedfractionraw10.9994%,afterdisplaykernel7.1908%;notunedtobepaperinterval.
No1000Dtermrealizations,syntheticselectiontrainingornewconstraintsrun.

## Outputs (all on14TB)
Base `/home/canoziia/GRMHD_DATA/analysis/`:
- `eht-fiducial-apr11-v132-pilot-v2` engineeringpilot;
- `eht-fiducial-apr11-v132-full-v1` simplifiedI diagnostic;
- `eht-official-I-apr11-v132-v1` releasedIstage;
- `eht-official-PD-apr11-v132-v1` inheritedI-noise diagnostic;
- `eht-official-I-cleanpol-apr11-v132-v2` releasedI+cleanpolarizationinput;
- `eht-official-PD-cleanpol-apr11-v132-v2` cleaninputevery-roundblur;
- `eht-official-PD-warmstart-apr11-v132-v1` fewerblurrestartcomparison;
- `eht-official-PD-warmstart-apr11-v132-render-v1` finalimageanditerationplots.
Eachcompletedstagehasdeclaredmanifest,sourceandpackagebytes,lineageandreceipts.

Tests:1473passed192.56s including12newplanning/APItestsand5officialwrappertests.
RuffhasoneE501(119-characterprovenancedescriptioninexecutedwarmstartscript).
Notgloballyignored,norchangedafterexecutiontoruinmanifestsourceidentity.
Allothercheckednewsourcepassed;thisisexplicitstyletodo,notnumericalfailure.

## What this establishes
We CAN reconstructfromofficialpolarizedUVFITSusingcurrentehtim,notjustreplotold
images. OfficialIreleasedstagesworkwith1.3.2viaexplicitdirectbackendadaptation.
P/Drealfitimprovesandproducesringpolarization,butthisfirstApr11runhasnotmet
convergence/fitcertification;itmustNOTberepresentedasofficialFigure6reproduction.
Beforefour-dayproduction/observedconstraints:examinefixedI/calibrationcovariance,
properPaperVIIsystematicdefinitionandobjectiveweightnormalization,currentvs
historicalreleasebehavior,restart/optimizerconvergenceandDtermstability,then
compareagainstverifiableofficialimage/metricdatawhenobtained. Noautomatic
fullhyperparameter/Dtermensemblebatchwaslaunched.

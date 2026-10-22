# Apr11 controlled validation: objective audit, initialization and restart

## Completed task
User authorized continuing the Apr11 baseline validation, not four-dayproduction
or newaxionconstraints. Newentrypoint `scripts/validate_eht_apr11.py` replaces
further branching experiments; oldexecutedscripts andmanifests remain untouched.
Explicitconfig `configs/eht-apr11-validation.yaml` and exact-byte audit approval
coveronlythe currentworkingmodel,NOTPaperVIIcertification. No token/turnbudgets
were applied to the final independent review; no child Goal was used.

Three controlled100-roundcases all completed EXIT0:
- baseline: uniform5%fraction,EVPA0°,no blur restarts;
- evpa45: same except initialEVPA45°;
- restart: EVPA0°,20uasblurrestartsonrounds2,12,22,32,42.
SameIcheckpoint,UVFITS,feedrotationmap,objectiveweights,noise,Dprior,100Psteps,
upstreamcalibrationbudgetsandgrid;allcasesactually100Drounds.
No case selected by visual agreement with paper.

## Objective audit result
Fullreportcopiedas `apr11-validation-evidence/OBJECTIVE_AUDIT.md`.
ReadPaperVII22/23/31renderedpages andinstalledofficialehtim1.3.2sources; bounded
upstreambytecomparison. Exacthistoricalβnormalizationisnotdocumented sufficiently
to certify aone-to-onemapping. Keepnorm_reg=True;do nottoggleittomakefitlookgood.
CurrentPobjective (CP,Cm are meancomplexχ² withdenominator2N):
`J=CP+Cm-2-100*Hraw/F+Traw*b/(F*pixel)`.
F=.6Jy,pixel1.875uas,b25.010878uas;HWnorm1.6666667×raw,PTVnorm22.2318918×raw.
Theseareexactcurrentfunctions,nothistoricalweightsconversion. Prior=init does
notcreateaper-roundPimageKLprior:HWdependscurrentm/fixedI,notpriorQ/U.

fsyshelperadds(.005|VI|)²topsigma²andRL/LRsigma²afterStokesconversion;nofactor2
bugestablished. Covariancepropagationisdiagonalapproximation,notfullcrosshand
noise. Currentinverse-JonesDsolveusesD-dependenterrorsandzerocenteredDpenalty
sum|D/.1|²;paperdescribesforwardJonesmodelanddoesnotidentifythisprior/scale.
Foronefrozenearlysnapshot,Dobjective4.99329+.639985=5.63328,exactreproduced;
priorisnotnegligible. pol_cal_new forwardAPIisfuturecontrolledcomparison,not
silentlysubstitutedinthisexperiment. Noformalofficialobjectivecertification.

## Actual final results
Inputs: officialreleasedIpipelineclean-polcheckpoint
`/home/canoziia/GRMHD_DATA/analysis/eht-official-I-cleanpol-apr11-v132-v2/`,
plusfixedofficialfeedtable. Same64²120uasgrid/fixedI.6Jy,positiveI;
pol HW100/PTV1/Pvis1/m1,norm_regTrue,fsys.005once,freshbaseDsolve,
AZ/LM/PVremainingDterms,otherpresentstationsheldzero.

|case|finalCP|finalCm|lastQUrelativechange|lastDmaxchange|
|---|---:|---:|---:|---:|
|baseline|2.655588|2.677915|.00326581|.00023147|
|initialEVPA45|2.551239|2.588424|.00233418|.00000852|
|scheduledrestart|2.317903|2.376780|.00197173|.00003806|

NumbersCP/Cmrecomputedindependentlyfromdirectmodelvisibilitiesandmatchedto
upstreamreducersatrtol1e-9. χ²normalizationis2Ncomplexsampleconvention,no fitted
Ndofsubtraction. Theseareworkingdata-fitstatistics,notqualitypass/failthresholds.
FinalPobjectivevalues−63.89069/−63.99389/−64.10883arepre-Datfixedcurrentdata;
notacompletecommonjointlikelihoodmodelselectionstatistic.

## Stabilization within each run vs between-run sensitivity
Commoncomparisonmask:568pixels,intersectionofI>=10%peakandP>=20%peakafter
sameadditional19uasIQUkernel;diagnosticnotclaimedofficialplotmask/effectivePSF.
PixelEVPAweightisdeclaredsymmetricfractionalP/I;nodataalignment/anglecorrection.

|case|19uasEVPARMS80→100|19uasEVPARMS90→100|QUrelativeL2(80→100)|
|---|---:|---:|---:|
|baseline|.04750deg|.02325deg|.0021983|
|initial45|.04074deg|.01564deg|.0026571|
|restart|.01413deg|.00967deg|.0006879|

Betweenfinalruns on EXACT SAME568pixels:
- initial45vsbaseline:EVPARMS2.49377deg,QUrelativeL2.1041527(10.42%);
- restartvsbaseline:EVPARMS4.92819deg,QUrelativeL2.2807175(28.07%).
Netpolarizationfractionsfullmap:baseline1.71485%,initial45.31145%,restart1.11150%.
Thesearecoherentsumsandcanamplifycancellation;notindependentuncertaintybars.
Thussmalllatechangesatdisplayresolutiondonotimplyinitialization-independentor
uniqueimage. Changeslargerbetweenrunsarenotexplainedawaybyquotingastopping
criterion. Donotchoosearestartbecauseitsdatametricissmallestorlooksmoreofficial.

## Independent residual diagnostics
`diagnose_eht_apr11.py` readsfinalDbydeclaredcheckpoint,correctsfreshsamebaseonce,
usesofficialdirectsample_uv prediction,thenindependentlyaggregatesI/Q/U/P/m
residualsbybaseline,stationand10minutetimebins,withreal/imag contributions.
Lastcheckpoints80/90/100arecomparedwithoutfittinganewtranslationorEVPArotation.
Noextraimaging/calibrationoptimizationorconstraintsinthisdiagnostic.

BaselinecaseworstPvisbaselineχ²:AA-PV6.76984(30conjugateincludedrows),
AA-AZ4.06604(32),LM-PV3.87101(16). Restart:AA-PV4.58870,AA-AZ3.56424,
AA-LM2.90649. Earlytimebinscanbeverylarge;someonly2or4rows. Stationgroups
andconjugatesoverlapandarenotindependenttrials;cannotdeclarefaultyatelescope,
sourcevariabilityornewphysicsfromthesegroupings.

## Validation and lineage
Finalfulloptional-environmentregression:1489passed200.37s,EXIT0.
Planning/APItests12passed;diagnosticpurehelpers4passed;16combinedpassed.
Independentstaticreviewconfirmedphysicalarray4slots,stationmapping,
actualclean-polinputs,noncumulativeerrors/sharedobjectiveandcaseisolation.
Planningtestidentifiednestedoutputrootloophole,fixedbeforeallsciencecases.
Allcompletedstagesdeclareinputs/packagebytes/helperfiles/checkpoints/outputs.
NooldNULLfunctionorresultwasreplaced,nofourdaybatchoraxioninferenceinitiated.

Styleboundary:oldexecutedwarmstartscriptretainsoneE501;newexecutedvalidation
runnerretainslonglinesdocumentedbylint. Sourcebyteswerefrozenintheaudit and
allthreecases,andwerenotchangedpost-hoc merelytosatisfystylechecks. No global
Ruff suppression added. Latercodeformattingrequiresanewversion/sourceidentity,
notresigningtheseexistingruns. Diagnosticnewfilelintpasses.

## Output locations and interpretation
Allunder `/home/canoziia/GRMHD_DATA/analysis/`:
- `eht-apr11-validation-{baseline,evpa45,restart}-v1`;
- `eht-apr11-validation-diagnostics-v1`:diagnostics.png,diagnostics.json,
  residual-groups.csv,common-mask.npy,manifest.json;
- `eht-apr11-baseline-diagnostic-smoke-v1`:earlierone-caseAPIvalidation.

**Conclusion:**currentworkingApr11reconstructionislocallystableatsmoothedbright
featuresoverlateiterations,butnotinitialization/restart-independent;unmodeled
residualsremainandexactPaperVIIobjectivenormalization/Dmodelisuncertified.
Beforefourdayproduction:resolveDobjective/normalizationwithrecoveredoriginal
codeorauthorinputs,orcarefullylabelalternativeworkingmodels;performcontrolled
forward-vs-inverseDlikelihoodcomparisonseparatelyfromprior/weightchanges.
Moreiterationsofexactlythesamesetupisnotbyitselfasolution.Nonewconstraints.

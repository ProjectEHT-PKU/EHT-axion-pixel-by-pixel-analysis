# Remote full EHT fixed-D uncertainty ensemble: 64 workers × 2 threads

## Authorized and verified resources

User authorized starting64×2 if account memory allowance is128GiB. SSH inspection
confirmed epyc-jwl:96physicalcores/192hardwarethreads,2NUMAnodes, account ancestor
`/user.slice/user-1000.slice/memory.max=137438953472` (128GiB),
`memory.high=128849018880` (120GiB). The per-SSH-session limit is unlimited but
is NOT the relevant account allowance. Hardwarefree memory is not permission.
Existing account usage was~25GiB. Controller uses100GiB account admission ceiling,
full-pool1GiB/worker reserve adjusted for materialized childRSS, and physical
headroom. It pauses admission, never kills science workers, and exits resumably
if idle resource-blocked for300s. No unbounded128worker launch.

## Actual code and data

New entrypoint `scripts/run_eht_uncertainty_batch.py` is input-only and independent
of old absolute live producer roots. Runtime sources and distribution metadata
are hashed into batch identity. No canonical common stage is forged or re-signed.
SourceSHA256 at launch:
`612fe695a8404e8e08db302a83a53ee798c2459b61ad061a2ae13929220d52e0`.

Verified export of fixedI, gain-calibrated native UVFITS and finalD from the FOUR
completed certified-P/tighter-D100additionalround local stages:
`/home/canoziia/GRMHD_DATA/EHT_DATA/official-uncertainty-inputs-v1/`.
34files,10,270,721bytes; originals and producer manifests retained intact.
Snapshotbundle.json lists originalpaths/SHA/bytes, four sourceI/D roots,
objective1/1/100/1normTrue, fixed64²120uas, residualcomponentSD.01, ddof1,
fsys.005 and source_kind=localrebuiltfourdayfiducial. Source roots are explicitly
export lineage, not assumed available or imported as live stages remotely.
Remote exported input bytes verified and mode made read-only.

All4000tasks are fixed-D solves from initialP=.05I real, NOT warm starts from
old pilot images. D distributions use this local-four-day mean/covariance,
not publishedTable5reference means or mixedsource covariances.1000draws/day,
masterseed20210405, day-specific deterministicseeds, independenthands/stations,
withinD real/imag covariance retained. No thermal noise bootstrap or D refit.

No objectiveweights,TVsmoothing,beam,calibrationrecipe orRAPTORparameterschanged.
ExactdiagonalPDHG gap1e-6 numericalcertificatesrequired. Defaultcap50k/600s;
uncertified tasks retainattempts and receive one200k/1800s numerical retry.
Errors, partial or uncertified draws never count as complete. Every attempt
stores raw/displayIQU, solver trace/stats/RSS, outputSHA andtaskreceipt.
Customfinish.json is written only after4000 numericalcertifications and explicitly
states scientifically_certified=false. This is not validation of allnoise/systematics.

## Installation and smoke

Newisolatedremoteproject copiedfromGit82a241b plusnew batchentrypoint:
`~/PHYS-production/eht-uncertainty-v1/{code,inputs,control,logs,results}/`.
LockedPixienvironment eht-rebuild, ehtim1.3.2. Directconda.anaconda.org downloads
were unreachable; only thisproject's `.pixi/config.toml` sets aTsinghuamirror.
Lock/packagehashesremainunchanged. No globalconfig or existingproduction modified.
InstallcompletedEXIT0. Originalunreachableattemptlogretained.

Two-tasksmokefrominitialPcompletedbothcertified:
-95-0000:6000iterations,gap9.56096e-7,solver16.65s,total43.34s,
  worker-lifetimepeakRSS.486GiB;
-96-0000:5500iterations,gap6.04392e-7,solver19.00s,total45.74s,
  peakRSS.594GiB.
Separateoutputroot:results/smoke-v1, no misleadingfullfinishreceipt.
Smokeexecutedsourcecopyretainedbeforefullresourceguardfix; fullusesnewidentity.
Tests3passed, newdriver/testsRuffpassed. Independentreviewtwo resourceguards
(in-flightreserve, boundedidleblock)fixedandunit-testedbeforefullstart.

## Full actual launch

Remote tmux:eht-uncertainty-64x2;64processworkers,eachBLAS2threads.
Rollinginterleavedfourdaytasks,not4000simultaneousprocesses.
Results `~/PHYS-production/eht-uncertainty-v1/results/full4000-v1/`.
Batchidentity:
`06a8f63baa95b7ab959c951e31ea45e6b6933a3db06455e54a65fe8c7a328af0`.

Firstchecks64active, noexhaustedtasks;workertotalRSS~31.5GiB,
accountmemory~49–52GiB, peak52.24GiB;projectedreservepeak~90.17GiB.
Launchsnapshot retainedunderremote-uncertainty-evidence:71/4000certified,
64active, noalert. Thisislaunchprogress,NOTcompletion.

Read-onlyremotemonitorrunsin eht-uncertainty-monitor every60s:
logs/monitor-status.json,monitor-history.jsonl,MONITOR_FINISHED.json.
Independentbackgroundmonitor(noGoal/token/turnbudget)pollsremote every120s and
reportsactualcompletion/exception/resourcefailure. Noautomaticworker kill,
scientificsettingchanges,oldoutputdeletionor fullretryloopfromscratch.
Resume verifies identical snapshot, executed sourceidentity, deterministicdraws
andcompletedartifactSHAs before skippingthem. Workers/threadcount canbeoperational
butthreadcountisrecordedinthisrunidentity; don'tsilentlychangeandrebindoldreceipts.

Fullfinishdoesnotyetprovideofficial-styleaperture/beta2summary,crossdayEVPA
covarianceoraxionlimits; these require a later clearlyidentified statisticsstage.
Original localinputs, old results andremote formerRAPTORproductionremainretained.

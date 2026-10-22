# Remote uncertainty run — completed

## Source and completion
- SSH source: `canoziia@epyc-jwl:/home/canoziia/PHYS-production/eht-uncertainty-v1`
- Controller tmux: `eht-uncertainty-64x2`; existing read-only monitor: `eht-uncertainty-monitor`.
- Finish receipt: `results/full4000-v1/finish.json`.
- Finish/run identity: `06a8f63baa95b7ab959c951e31ea45e6b6933a3db06455e54a65fe8c7a328af0`.
- Source snapshot bundle identity (batch metadata, not rehashed): `9b46f39dca06dc3087ac32f619ad58fc0cf5ca793ca71a7263a07ff13549756e`.
- Actual controller exit: **0**, independently present in `logs/full4000.log` and final monitor status.
- **4000 numerically certified**, 4000 finish receipt keys; days 95, 96, 100, 101 each have exactly 1000 successful receipts. Active/pending/exhausted: 0/0/0.

## Timing and resources
- First attempt metadata mtime (run-start proxy): 2026-10-22 17:01:32.280316 UTC.
- Finish receipt mtime: 2026-10-22 19:04:37.790735 UTC.
- Elapsed: 7385.51 seconds (2 h 3 min 5.51 s); throughput: 1949.76 certified tasks/hour.
- Existing monitor history: 123 samples; observed account cgroup memory peak including recorded running peak: 59,119,169,536 bytes (55.06 GiB). Peak projected memory: 96,819,171,328 bytes (90.17 GiB).
- Limits: memory.high 120 GiB; hard maximum 128 GiB. No monitor alerts recorded. Sampled monitoring does not prove absence of unsampled instantaneous excursions.

## Attempts and bounded integrity checks
- Successful-attempt histogram: attempt 1 = 3985; attempt 2 = 15.
- All preserved task-receipt statuses: certified = 4000; uncertified = 15. Thus **15 failed/uncertified attempts**, recovered by retry; no task-error status receipts, no exhausted tasks.
- Read all 4000 successful task receipts: identity/status/task-key checks passed; every output named by these receipts exists. No bad receipts or missing referenced outputs.
- Read four certificate samples, one per day: all numerically certified, best gaps 6.04e-7 to 9.56e-7, below 1e-6.
- No bulk FITS reads, no full output hashing, no large output copy. Copy-stage hash validation remains the parent's responsibility.
- Machine-readable bounded-check summary: `/tmp/eht-remote-uncertainty/completion-summary.json`.

## Numerical certification is not uncertainty-coverage validation
The finish explicitly states `scientifically_certified: false`. Floating-point primal/historical dual-gap certification is not interval arithmetic and does not establish scientific acceptance, polarization uncertainty coverage, or true inter-day restrictions. This run uses the local rebuilt 4D fiducial (not published), 1000 draws/day, ddof=1, residual SD 0.01/component, fsys=0.005, native circular basis, initial P=0.05 I real, exact PDHG gap tolerance 1e-6, two threads/worker. First budget: 50k iterations/600 seconds; retry: 200k/1800 seconds. Canonical NULL must remain as the downstream control.

## Parent handoff — not performed here
Remote complete. Copy-back source is `canoziia@epyc-jwl:/home/canoziia/PHYS-production/eht-uncertainty-v1/`, especially the entire `results/full4000-v1/` tree including preserved attempts and provenance/metadata/logs needed for reproduction. Local bulk destination has not been selected or created by this monitoring agent; `/tmp/eht-remote-uncertainty/` contains only this report and small summary, not the results.

Per the user's latest authorization, only after receiving this completion handoff should the parent open persistentGoal, copy all results back and validate hashes, then undertake the requested complete 120 uas observed-polarization-error plus genuine day-difference restriction analysis, retaining canonical NULL as control. This agent has not opened Goal, changed science/code/locks/outputs, restarted/killed jobs, or initiated extra science.

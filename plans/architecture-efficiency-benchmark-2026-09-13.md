# Architecture efficiency audit: before/after benchmark

Date: 2026-09-13. Branch `playground`, `/data/work/EMC2_dev_oo`.

Compares the package before the audit plan (`62801947`, the parent of C1) with
the package after it (`ad491d90`: C1-C17, plus two fixes this benchmark
prompted). The later liveness fix `4f6f4ae1` changes only how dying worker
processes are detected; it was not re-benchmarked.

## Bottom line

- End-to-end hierarchical fits changed by -2.4% to +3.4%, with identical
  posteriors in all 17 cases. The largest realistic case (140 subjects, P = 64)
  went from 999.7 s to 975.5 s (-2.4%). Eleven of the twelve P = 64 cases were
  0.3-2.4% faster; the exception is 8 subjects at +3.4%. The P = 8 cases
  ranged from -1.1% to +1.8%.
- The large gains are in the likelihood call for wide designs. The 256-cell
  cliff is gone (width 257: 113.6 -> 45.0 ms, -60%). Width 256 is 18% faster,
  and 64-level designs are 7-8% faster per call for BAwL, RDM and DDM.
- Narrow designs, the single particle step and the worker pool are level
  with the pre-plan code.
- There is no consistent change in peak memory. Single-run peaks are too
  noisy to support one.
- The benchmark found two defects, both now fixed and committed: a 12%
  worker-pool slowdown introduced by the plan, and a pre-existing silent
  chain stall in direct `run_stages()` calls.

Realistic fits barely moved because the likelihood kernel is about 85% of a
particle step, and the part of the call the plan optimised was already small:
the parameter prologue is 1.8 ms of a 37-40 ms narrow-design call. The
mapper's big win applies past 256 design cells, and none of the realistic
fixtures reach that. The P = 250 `stress` fixtures that would were left out as
unrepresentative of real fits.

## Method

- **Builds.** Clean `git archive` exports of each commit, built with
  `R CMD INSTALL` under the same configure profile (legacy, g++-13,
  `-O3 -march=native -ffast-math`). A second independent install of
  `62801947` provided an A/A check; its `EMC2.so` was byte-identical to the
  first install's.
- **Environment.** 32-vCPU EPYC Genoa VM, 62 GB, OpenBLAS pthreads with
  `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`. This matches the
  RStudio session's `~/.Renviron`. Nothing else was running on the machine.
- **Ordering.** Pre and post runs alternated, and the order flipped from case
  to case, so drift affected both builds equally. The likelihood benchmarks
  ran twice in each order.
- **End-to-end fits.** `WorkingTests/bench_hierarchical_regression.R` with
  fixed iteration counts: preburn, burn and adapt 100 each, sample 300,
  step 100. Wall time covers the whole fit including preparation. The only
  change to the script was to take its PSS memory reader from a standalone
  copy of HEAD's `.emc_profile_memory`, so both builds were measured by the
  same code; the pre-plan namespace does not have that function.
- **Equivalence.** Every comparison checks the posterior: the alpha and
  subject log-likelihood checksum for fits, and the final sampler state for
  the pool benchmark. All were identical between pre and post.

### Noise floor (A/A: two installs of the pre-plan commit)

| Case | pre | pre2 | Difference | Peak memory pre / pre2 |
|---|---|---|---|---|
| smoke | 26.9 s | 26.7 s | -0.4% | 1.32 / 1.28 GB |
| n32_p8 | 47.1 s | 47.3 s | +0.3% | 2.10 / 2.11 GB |
| cores4 | 225.3 s | 225.3 s | -0.0% | 3.20 / 3.13 GB |

Wall-time noise between installs is below 0.5% on fits. On the shorter pool
benchmark it reached about 4% in one earlier round, which is why every pool
figure below is a median of repeated, interleaved runs.

## Defects found by this benchmark

### 1. The plan slowed the worker pool by 12% (fixed in `c9ad349e`)

The pool preburn benchmark (LBA, 24 subjects, 200 trials, 4 workers,
100 particles), run three times interleaved with `search_width` corrected,
took 2.05 -> 2.30 s at 20 iterations (+12%) and 6.72 -> 7.09 s at 200 (+5.5%),
with identical results. Median master CPU at 200 iterations rose from 1.07 to
1.56 s. Wrapping the pool functions with timers in both builds found four
causes:

| Cause | Commit | Measured cost |
|---|---|---|
| `.emc_wpool_done_take` asked `readBin()` for 64 KiB on every poll of the wait loop | C4 | GC 0.17 -> 0.41 s over 200 iterations; fixing it together with the connection count (next row) cut master CPU from 1.82 to 1.11 s (pre-plan: 1.00-1.17 s) |
| `.emc_wpool_feasible_workers` counted connections with `showConnections(all = TRUE)` | C5 | ~0.2 s per pool start |
| `.emc_kernel_totals` built a data frame twice per worker request, even with measurement off | C9 | ~200 us per call, against ~1 us for the flag check |
| ST_START records made the caller re-enter the wait loop from scratch, restarting the spin | C4 | 1600 wait calls per 200 iterations instead of 800; ~14 us per empty poll |

The pool-function timings in this table were taken before the `search_width`
correction (defect 2), so from iteration 26 the sampler was repeating states.
The pool I/O costs they measure don't depend on that; the kernel-totals cost
is a standalone microbenchmark.

After the fixes the pool benchmark is level with the pre-plan code (table
below). What remains, roughly 0.5 ms per iteration, is the cost of bounded
non-blocking reads and the larger completion record.

### 2. `run_stages()` without `search_width` stalled every chain at iteration 26 (fixed in `ad491d90`)

This predates the plan. `run_stages()` defaults to `search_width = NULL`, and
`set_p_accept()` computes `0.02 * (1/NULL)`, which is `numeric(0)` rather than
an error. At the first adaptation (iteration 26) every subject's epsilon
became empty, so the second proposal component's covariance was NA and about
half of each proposal cloud was NaN. Every weight was then NA, and
`safe_new_particle()` repeated the previous state. `fit()` and `run_emc()` pass
`search_width = 1`, so ordinary fits were never affected. Direct callers
were, including `WorkingTests/bench_large_fit_sampling.R`, the plan's own
reproduction benchmark: its 200-iteration timings had been measuring caught
errors from iteration 26 onward. C6's failure classification is what exposed
it ("4200 updates were answered by repeating the previous state").

### 3. Worker liveness could flip from dead back to alive (fixed in `4f6f4ae1`)

Found while fixing the intermittent failure in "workers whose reaper has died
are terminated by their owner". `.emc_wpool_pid_alive()` reported a process
as alive if it was reaped between its `kill(pid, 0)` probe and its `/proc`
read, or if it was in state X. Polling 400 exiting processes gave 230 flips
from dead back to alive; after the fix there were none in 6.6 million polls.
The test's own premise, that workers outlive a killed template, had also
stopped holding once workers armed the parent-death signal: they now die
36-42 ms after the template.

## Results

### Likelihood call: one subject, 2000 trials, 100 particles

`plans/architecture-efficiency-bench.R`. Figures are ms per call, each the
mean of two runs (pre-post order and post-pre order); "prologue" is the
parameter mapping, transform and bounds step.

| Fixture | Pre | Post | Change | Prologue pre | Prologue post |
|---|---|---|---|---|---|
| BAwL, narrow | 40.7 | 40.3 | -0.8% | 1.8 | 1.8 |
| BAwL, 64 levels | 45.0 | 41.8 | -7.0% | 6.0 | 3.2 |
| RDM, narrow | 37.2 | 37.2 | +0.0% | 1.8 | 1.8 |
| RDM, 64 levels | 41.6 | 38.6 | -7.3% | 6.1 | 3.1 |
| DDM, narrow | 17.5 | 17.4 | -0.2% | 0.9 | 0.9 |
| DDM, 64 levels | 19.4 | 17.8 | -8.0% | 2.7 | 1.3 |
| RDM, covariate with 255 cells | 37.5 | 37.2 | -0.8% | 1.9 | 1.9 |
| RDM, 256 cells | 37.4 | 37.1 | -0.8% | 1.9 | 1.9 |
| RDM, 257 cells | 39.9 | 37.2 | -6.7% | 4.0 | 1.9 |
| RDM, 512 cells | 39.5 | 37.4 | -5.5% | 4.0 | 2.1 |

Small particle batches through `calc_ll_manager` (1, 2, 5, 25 and 100
particles) are unchanged, within -3.4% to +0.4%.

### Wide designs: the 256-cell cliff

`AUDIT_MODE=wide`, RDM with `v ~ 0 + condition`.

| Width | Pre | Post | Change |
|---|---|---|---|
| 256 | 54.6 ms | 44.8 ms | -17.9% |
| 257 | 113.6 ms | 45.0 ms | -60.4% |

The same fit with one unused 257th design row appended took 53 -> 113 ms
before the plan (2.12x) and 45 -> 45 ms after it (1.00x), with bit-identical
log-likelihoods.

### LBA prologue benchmark (CPU ms)

`WorkingTests/bench_likelihood_prologue.R`.

| Fixture | Whole call pre -> post | Prologue pre -> post |
|---|---|---|
| All-finite RTs | 36.40 -> 37.05 (+1.8%) | 1.60 -> 1.58 |
| 15% omissions, UC = 2.0 | 33.42 -> 33.23 (-0.6%) | 1.27 -> 1.40 |

Cost against particle count (1 to 100 particles) is unchanged, within -2.6% to +1.5%.

### One subject's particle step: 500 steps, 200 trials

`AUDIT_MODE=particle`. With A/A control, BAwL, RDM and DDM were all
within +/-1% across repeated rounds; this run gave BAwL 2.38 -> 2.37 s,
RDM 2.14 -> 2.13 s and DDM 1.12 -> 1.14 s. A 2-4% slowdown in an earlier
single run was install noise.

### Preburn through the worker pool

`WorkingTests/bench_large_fit_sampling.R`: LBA, 24 subjects, 200 trials,
4 workers, 100 particles, spawn backend. Final sampler states were identical.

| Iterations | Pre | Post (`ad491d90`) | Change |
|---|---|---|---|
| 20 | 2.08 s | 2.08 s | -0.1% |
| 200 | 6.73 s | 6.85 s | +1.7% |

Each figure is the mean of two runs, one in each order. Before the pool fixes,
a separate round of three interleaved runs gave 2.05 -> 2.30 s (+12%) and
6.72 -> 7.09 s (+5.5%); see defect 1.

### End-to-end hierarchical fits

RDM; preburn/burn/adapt 100 each, then sample 300. Wall time covers the whole
fit. "Stage change" gives preburn/burn/adapt/sample in percent. Peak memory is
the PSS of the whole process tree. The posterior was the same in every case.

| Case | N | P | Chains | Cores/chain | Particle factor | Pre | Post | Change | Stage change | Peak memory, pre -> post |
|---|---|---|---|---|---|---|---|---|---|---|
| smoke | 8 | 8 | 2 | 2 | 50 | 26.9 s | 27.2 s | +1.4% | +1/+1/+5/+0 | 1.32 -> 1.34 GB |
| n8_p8 | 8 | 8 | 2 | 4 | 50 | 22.2 s | 22.0 s | -1.1% | -1/+8/-5/-3 | 1.58 -> 1.65 GB |
| n8_p64 | 8 | 64 | 2 | 4 | 50 | 81.3 s | 84.1 s | +3.4% | -4/+3/+6/+5 | 4.67 -> 4.72 GB |
| n32_p8 | 32 | 8 | 2 | 4 | 50 | 47.1 s | 47.8 s | +1.5% | +3/+1/+1/+1 | 2.10 -> 1.94 GB |
| n32_p64 | 32 | 64 | 2 | 4 | 50 | 250.7 s | 247.3 s | -1.4% | -3/-2/-2/-1 | 6.84 -> 7.35 GB |
| n140_p8 | 140 | 8 | 2 | 4 | 50 | 179.0 s | 181.0 s | +1.1% | +1/+1/+2/+1 | 3.15 -> 3.15 GB |
| n140_p64 | 140 | 64 | 2 | 4 | 50 | 999.7 s | 975.5 s | -2.4% | -4/-2/-3/-2 | 15.16 -> 14.79 GB |
| cores1 | 32 | 64 | 1 | 1 | 50 | 724.7 s | 711.5 s | -1.8% | -5/-5/+1/-1 | 1.17 -> 0.89 GB |
| cores2 | 32 | 64 | 1 | 2 | 50 | 403.1 s | 393.8 s | -2.3% | -2/-4/-2/-2 | 1.88 -> 1.88 GB |
| cores4 | 32 | 64 | 1 | 4 | 50 | 225.3 s | 224.4 s | -0.4% | -5/-2/+3/-0 | 3.20 -> 2.99 GB |
| cores8 | 32 | 64 | 1 | 8 | 50 | 149.0 s | 147.3 s | -1.2% | -3/+2/-3/-1 | 5.59 -> 5.24 GB |
| cores16 | 32 | 64 | 1 | 16 | 50 | 94.5 s | 92.7 s | -1.9% | +1/-5/+0/-3 | 10.11 -> 9.14 GB |
| chains2 | 32 | 64 | 2 | 4 | 50 | 247.6 s | 246.0 s | -0.7% | -1/-2/+1/-1 | 6.63 -> 6.71 GB |
| chains4 | 32 | 64 | 4 | 4 | 50 | 300.0 s | 299.1 s | -0.3% | -2/-1/-1/+0 | 15.31 -> 13.75 GB |
| particles25 | 32 | 64 | 1 | 4 | 25 | 168.5 s | 166.1 s | -1.4% | -3/-9/-1/-0 | 3.25 -> 2.99 GB |
| particles100 | 32 | 64 | 1 | 4 | 100 | 336.8 s | 332.6 s | -1.2% | -4/-5/-2/+1 | 3.12 -> 3.11 GB |
| n32_p8, step 10 | 32 | 8 | 2 | 4 | 50 | 175.5 s | 178.6 s | +1.8% | +2/+4/+2/+1 | 1.68 -> 2.07 GB |

`cores4`, `particles50` and `chains1` specify the same fit, so it was run once,
as `cores4`.

### Memory

There is no consistent change. Single-run peak PSS depends on when the
0.25 s memory poll lands relative to worker recycles. The step-10 case read
1.68 -> 2.07 GB in the main run and 2.07 -> 1.95 GB (pre -> post) when
repeated. That repeat also counted live pool processes every 0.1 s: neither
build ever exceeded the expected 10 (2 chains x (1 template + 4 workers)),
so retired workers never overlapped the next pool. The larger peaks at high
core counts were lower after the plan (cores16 -9.6%, chains4 -10.2%), but
these are single runs.

## Not measured

- The P = 250 `stress` fits, which are the only end-to-end cases above 256
  design cells. The likelihood-level wide benchmark covers the cliff instead.
- C17's converged particle-count comparison (several chains, full length).
- Mixture quality, blocking and parameterisation.
- Any platform other than this Linux machine.

## Reproducing

From the repository root, with `lib_pre` and `lib_post` built from clean
`git archive` exports of the two commits:

```sh
# likelihood level (no EMC_LIB support in this script: use R_LIBS)
R_LIBS=lib_pre  AUDIT_MODE=likelihood Rscript plans/architecture-efficiency-bench.R
R_LIBS=lib_post AUDIT_MODE=wide AUDIT_REPS=6 Rscript plans/architecture-efficiency-bench.R
R_LIBS=lib_post AUDIT_MODE=particle Rscript plans/architecture-efficiency-bench.R
EMC_LIB=lib_post Rscript WorkingTests/bench_likelihood_prologue.R

# pool preburn (the script now passes search_width = 1)
EMC_LIB=lib_post EMC_N=24 EMC_TRIALS=200 EMC_ITER=200 EMC_WORKERS=4 EMC_PARTICLES=100 \
  Rscript WorkingTests/bench_large_fit_sampling.R

# end-to-end fits, one case at a time, alternating libraries
EMC_CASES=n140_p64 EMC_LIB=lib_post EMC_LABEL=post EMC_OUT=n140_p64_post.rds \
  Rscript WorkingTests/bench_hierarchical_regression.R
```

Put `EMC_LIB=` before `Rscript`. Placed after, it becomes a script argument
and the stale system EMC2 loads instead. Check the namespace path each script
prints.

The raw outputs (per-case RDS and logs, the driver script and the comparison
script) were written to the session scratchpad under `fullbench/`, which is
not part of the repository.

# Large-fit sampling architecture benchmark

Date: 2026-08-06

Baseline: `playground` at `138ffe16`

Candidate: `perf/large-fit-sampling`, created from the same commit

The benchmark driver is `WorkingTests/bench_large_fit_sampling.R`. Both
revisions were installed into separate libraries and run in fresh R processes
with matched L'Ecuyer seeds. The installed-package runs used the clean spawned
worker backend. This host exposes one logical CPU, so these results measure
allocation, process-management, and IPC costs reliably, but they cannot establish
multi-core speedup.

## Result

Gemini correctly identified two useful targets: repeated history-array copying
and over-frequent worker recycling. Sending every subject's current state to
every worker was also unnecessary, but eliminating that broadcast was a much
smaller win than claimed in the tested fits. Persistent subject ownership was
not implemented: static ownership would discard the existing measured-cost LPT
balancing, complicate failure recovery, and was not justified by the residual
dispatch cost after slicing.

The candidate now:

- sends each worker only its assigned `alpha`, subject population means,
  likelihoods, adaptation state, and RNG streams;
- broadcasts only the population covariance and precomputed Cholesky cache;
- writes sample-history slices and copies extension prefixes with small native
  routines, retaining the public R matrix/array representation;
- explicitly detaches the compact initial sample store because `make_emc()`
  initially shares it between replicated chains;
- recycles clean spawned workers every 50 iterations by default, while retaining
  the conservative 10-iteration default for inherited/fork workers;
- offers opt-in per-iteration phase, IPC-size, and worker timing through
  `options(emc2.sampler_profile = TRUE)`.

## Matched measurements

| Benchmark | Current | Candidate | Change |
|---|---:|---:|---:|
| History writes, P=40, N=140, 1,000 saved + 50 writes | 5.139 s | 0.012 s | 428x faster |
| LBA fit, P=5, N=140, 200 trials, 80 iterations, default recycling | 10.717 s | 8.562 s | 20.1% faster |
| Same LBA fit, recycle=50 on both revisions | 9.150 s | 9.189 s | candidate 0.4% slower |
| High-P LBA fit, P=43, N=140, 80 trials, 20 iterations, recycle=50 (healthy pools only) | 5.935 s | 5.793 s | 2.4% faster |

Low-P and history times are medians of three interleaved fresh-process runs.
Final `alpha` and subject log-likelihood arrays were identical between revisions
for every matched fit. The low-P default result shows that its end-to-end
improvement comes from less recycling at this short history length. At equal
recycling, sliced IPC is neutral for low P.

The high-P protocol stress exposed a reliability difference: 7 of 9 baseline
runs lost the worker pool with SIGPIPE and correctly completed through the
serial fallback, versus 0 of 3 candidate runs. The high-P timing table excludes
degraded runs (two clean baseline runs and three clean candidate runs), so it is
useful directional evidence rather than an equally replicated performance
estimate. Including degraded runs would conflate transport reliability with
steady-state speed.

With P=43 and four workers, the candidate profiler measured a 44,944-byte shared
payload per active worker and 157,088 private bytes in total, approximately
336,864 request bytes per iteration. The subject-scaled matrices now cross the
worker pipes once in aggregate rather than once per worker, and the population
covariance is reused from the factor cache instead of appearing twice in the
shared blob. The remaining shared payload is dominated by O(P^2) factor data;
stateful subject ownership would not remove it.

## Recycling and memory

In an 80-iteration, N=140, 200-trial run, recycling every 10 iterations reduced
peak proportional-set size from 625,574 KiB to 611,124 KiB (2.3%) but increased
elapsed time from 8.427 s to 12.441 s (47.6%). In the shorter memory probe, the
PSS difference was only 240 KiB while recycling cost 12.8% elapsed time. These
measurements support a less aggressive default for the clean spawned backend,
not disabling recycling universally.

## Assessment of the five claims

1. **Stateful workers:** directionally valid for IPC, but the suggested static
   ownership conflicts with dynamic load balancing and fault recovery. Slicing
   removes the N-by-worker broadcast without moving authoritative chain state
   out of the master. The remaining measured gain did not warrant the added
   state protocol.
2. **The 64 KiB cliff:** large blocking pipe writes can impose backpressure and
   serialize dispatch, but exceeding pipe capacity is not by itself a deadlock.
   Slicing is still worthwhile and is now implemented. At P=43 the shared blob
   is 44,944 bytes, while a complete per-worker request averages about 84 KiB
   after its private subject state; the blocking protocol must therefore remain
   correct above the nominal buffer size.
3. **Worker recycling:** correct for the clean spawned backend. Recycling every
   10 iterations bought little PSS reduction in the measured large fit and had
   substantial process-creation cost. Fork/inherited workers retain the old
   conservative default because their copy-on-write exposure differs.
4. **History-array copying:** correct and the largest isolated optimization.
   Native writes require an ownership invariant; without explicitly detaching
   replicated initial stores they can mutate sibling chains. The implementation
   and tests cover that case.
5. **Two-stage estimation:** not an optimization of the same sampler. It changes
   the inferential method and should be evaluated as a separate modeling choice,
   not substituted to fix PMwG implementation overhead.

## Reproduction

Install each revision into a separate R library, then run the driver in a fresh
process, for example:

```sh
EMC_LIB=/path/to/library EMC_LABEL=current EMC_BENCH_MODE=fit \
EMC_N=140 EMC_TRIALS=200 EMC_ITER=80 EMC_WORKERS=4 \
EMC_PARTICLES=25 EMC_RECYCLE=default EMC_OUT=/tmp/current.rds \
Rscript WorkingTests/bench_large_fit_sampling.R
```

Use `EMC_BENCH_MODE=history`, `EMC_P=40`, `EMC_SAVED=1000`, and
`EMC_WRITES=50` for the history test. Set `EMC_PROFILE=true` for phase and
wire-size diagnostics; profiling is off by default because measuring serialized
message sizes itself performs extra object walks.

# Overhead baseline and optimization gates, 2026-11-11

Baseline `7e161f75`; warming candidate `0c4d5112`. Frozen isolated R
installations, R 4.6.1, GCC 13, CPU affinity 0–26, OMP/OpenBLAS/MKL thread
counts 1, machine benchmark lock. Build logs, session information, raw fitted
objects and profiles are retained in the ignored `results/` directory. Data
are generated with seed 123; sampler seeds are 101, 202 and 303. See
`PROTOCOL.md` for exact stopping criteria and the proposed later-step gates.

| Fixture | Baseline elapsed (s) | Warmed elapsed (s) | Comparison |
|---|---:|---:|---|
| EXG default fit, 3 seeds | 21.797, 22.554, 23.058 | 17.025, 17.029, 17.042 | Samples, RNG and compared continuation fields identical |
| LNR default fit, 3 seeds | 104.929, 105.700, 107.337 | 101.444, 102.699, 101.955 | Samples, RNG and compared continuation fields identical |
| EXG fixed transitions, 3 seeds | 20.556, 20.538, 20.802 | 16.774, 16.539, 16.782 | Samples, chain/caller RNG and compared continuation fields identical |
| LNR fixed transitions, 3 seeds | 134.651, 133.838, 135.893 | 131.864, 132.119, 132.900 | Samples, chain/caller RNG and compared continuation fields identical |
| EXG precision target, 3 seeds | 21.313, 21.070, 21.064 | 16.785, 17.031, 16.785 | Target reached at 400 draws/chain for every run |
| LNR precision target, 3 seeds | 148.668, 134.675, 134.609 | 145.943, 131.603, 131.334 | Target reached at 500, 400, 400 draws/chain in both arms |
| Single-subject likelihood pool | 60.795 | 56.989 | Samples, RNG and compared continuation fields identical |
| Hierarchical forced rebuild, partial concurrency | 81.975 | 82.126 | Draws and RNG identical; cost state differs |
| Pooled, ordinary allocator | 91.577 | 88.083 | Draws, RNG and allocator state differ |
| Pooled, fixed allocation | 126.176 | 120.488 | Draws and RNG identical; cost state differs |

The fixed-transition runs use 50 preburn, 200 burn, 200 adapt, and 400 sample
iterations. The precision target requires minimum alpha ESS above 800 and at
least 400 sample iterations. Paired precision runs stopped at the same draw
counts with the same ESS. Two additional seed-101 fixed-transition timing pairs, run candidate first,
took 20.041 and 20.287 s baseline versus 16.305 and 15.784 s warmed for EXG,
and 132.555 and 133.071 s versus 130.566 and 130.830 s for LNR.
Individual quantiles, quantile MCSE, minimum bulk and
tail ESS per draw, and ESS per second are in `results/comparison.txt`; the default-fit equivalents
are in `results/default_comparison.txt`. The EXG baseline default fit reproduces
the plan's approximately 22.5 s starting measurement.

A pooled baseline-versus-baseline repeat changed all chains' samples and RNG
streams, while all 75 pooled parameter quantiles remained within 2.53 combined
Monte Carlo errors. Holding allocation fixed restored exact baseline/candidate
draw and RNG identity. The hierarchical baseline repeat changed allocator and
HMC cost fields while draws and RNG remained identical. These controls show
that whole-fit bit identity is not a useful correctness gate for this ordinary
allocator: timing affects its later schedule even without the optimization.
The candidate is restored in `50cb57e7`; the pooled posterior comparison below
assesses the resulting schedule variation directly.

### Ordinary pooled allocation posterior check

Fourteen paired seeds used the same pooled LBA data (hash
`d3650406e44a8b6307adab451250714a`) and 400 production draws in each of three
chains (42 chains per arm). The corrected analysis covers 37 quantities:
subject parameters, group means, group variances, derived log group SDs, and
the pooled parameter. Per-seed checks found six quantile exceedances out of
1,554 comparisons, all from seeds 404 and 606; the largest was 5.99 combined
MCSE. Combining all 42 chains, all 111 parameter quantiles were within four
combined MCSEs (largest 1.85).

This is a fixed-work comparison, not a comparison of default stopping rules.
The fixture explicitly requests preburn 50, burn 200, adapt 200 with
`min_unique = 1`, and sample 400. That bypasses both packages' default burn and
sample convergence stops. For standard hierarchical fits, ordinary `fit()` in
both packages defaults to R-hat checks on subject parameters (`alpha`) and group
means (`mu`). The additional scale criteria the user flagged belong to DMC2's
separate `hardfit()` defaults: at least 1,000 sample draws, the 99th-percentile
alpha and maximum `mu`/`sigma2` R-hats below 1.1, and minimum bulk ESS of 100
over group means and log group SDs. `hardfit()` also checks likelihood gaps and
subject spread; its hierarchical guidance recommends six chains. EMC2 `fit()`
does not apply those `hardfit()` scale/ESS checks by default. The 14-seed
fixture has three chains and does not implement the extra gap/spread checks, so
its logged scale diagnostics are partial post-hoc checks, not a claim that
`hardfit()` passed. Every 400-draw fit fell short of its 100 group-mean/log-SD
bulk ESS floor.

Mean elapsed time was 93.924 s baseline and 90.226 s candidate, a 3.93% paired
reduction. The paired geometric speedup was 1.041 [1.029, 1.053] (95% t
interval). Across-seed minimum group-mean/log-SD ESS ratios were wide: 1.531
[0.803, 2.919] for bulk and 1.422 [0.892, 2.268] for tail ESS. Across the
eight group-mean/log-SD coordinates, median bulk and tail ESS ratios were
1.001 and 0.969; the 10th percentiles were 0.776 and 0.707. The 0.90
non-inferiority claim was unresolved at 400 draws. The summaries are in
`results/pooled14_posterior.txt` and `results/pooled14_aggregate.txt`.

### Pooled 1,200-draw scale diagnostics and same-build control

Four paired seeds (1506, 1607, 1708 and 1809) extended the same fixed-work
comparison to 1,200 sample draws per chain. Mean elapsed time was 699.399 s
baseline versus 651.615 s warmed, a 6.82% reduction; geometric speedup was
1.073 [1.043, 1.105]. All 111 aggregate parameter quantiles were within four
combined MCSEs (maximum 3.52). The partial DMC2 `hardfit()` R-hat-plus-ESS
subset passed on 3/4 baseline fits and 1/4 warmed fits. This omits the
likelihood-gap and subject-spread checks and uses three chains.

Across the eight DMC2 group-mean/log-SD coordinates and four paired seeds, the
candidate/baseline ESS ratios had median 0.572 (bulk) and 0.654 (tail); 29/32
bulk and 22/32 tail ratios were below 0.90. The geometric mean of each run's
minimum group/log-SD ESS ratio was 0.338 [0.111, 1.027] for bulk and 0.538
[0.275, 1.052] for tail. These fixed-work results raise an ESS concern despite
the runtime gain.

Two more paired 1,200-draw fixed-work runs used the frozen baseline library in
both runner arms. The partial `hardfit()` subset passed on 0/2 fits in each
arm. Same-build per-coordinate ESS ratios had median 1.150 (bulk) and 1.123
(tail); 5/16 coordinates in each measure were below 0.90, with minima 0.395
and 0.155. This confirms substantial run-to-run schedule variation.

An additional pinned, machine-locked run measured time to a shared scale-aware
target patterned after DMC2 `hardfit()`. Both EMC2 builds were checked against
this same external target; it is not either package's default `fit()` stop or
an exact replay of DMC2 `hardfit()`. DMC2's rule can recheck after dropping the
oldest third of sample draws when its R-hat check fails. This continuation gate
checks each full sample prefix and does not apply that discard-and-recheck path.
It also omits likelihood-gap, subject-spread and culling rules and uses three
chains. It used the original 1,200-draw fits as starting states,
checked at 100-draw intervals, and continued misses up to 4,000 draws. The
table reports the first passing check and elapsed time; a miss at 4,000 is
censored at the cap.

| Seed | Baseline draws / time (s) | Warmed draws / time (s) |
|---:|---:|---:|
| 1506 | 1,500 / 207.4 | Missed at 4,000 / 427.4 |
| 1607 | 1,000 / 158.6 | 1,200 / 163.2 |
| 1708 | 1,000 / 154.8 | 1,900 / 230.7 |
| 1809 | 1,000 / 157.0 | 2,000 / 238.2 |

Across the three pairs where both builds passed, warmed/baseline geometric
time-to-target speedup was 0.755 [0.438, 1.301]. The two-seed same-build
time-to-target control produced ratios of 1.311 and 1.211 between its runner
arms, geometric mean 1.260 [0.764, 2.080]. The paired result points toward
slower convergence after warming, and one warmed fit missed the cap, while the
same-build control and wide intervals show that this four-seed gate cannot
separate the effect from allocator variation. Keep the PSIS warm change under
evaluation; do not claim a time-to-DMC2-hardfit advantage from these data.

Raw summaries are in
`results/pooled_long_posterior.txt`,
`results/pooled_long_repeat_posterior.txt`,
`results/pooled_stop_summary.txt`,
`results/pooled_scale_control_summary.txt`, and the corresponding ignored RDS
files under `results/pooled_long/`, `results/pooled_long_repeat/`,
`results/pooled_stop/` and `results/pooled_scale_control/`.

### Step 4a: shared-tuple REXG truncation cache

The Step 4a cache stores natural and log truncation terms for up to 16 exact
`(mu, sigma, tau)` tuples per raw REXG batch. Batches with more than 16 tuples
retain the prior row-by-row route. On the 704-row, 256-particle fixture, all
256 outputs were bit-identical to the frozen baseline; median call time over
ten repetitions was 0.028 s before and 0.014 s after (2.00x). A second fixture
with four parameter cells and 19,492 rows also produced bit-identical outputs,
with median time 0.934 s before and 0.623 s after (1.50x). A full EXG
fixed-transition fit using the warmed pre-cache build versus the cache build
took 16.529 s and 11.255 s; every chain's samples, RNG, and compared
continuation fields were identical. Related EXG truncation, tail, stop-signal
and raw-path tests passed 113 expectations; one optional development-validation
test was skipped. See `results/rexg_cache_comparison.txt`,
`results/rexg_cells_comparison.txt`, `results/rexg_kernel_fit_comparison.txt`,
and `results/rexg_kernel_tests.log`.

The configurable-bound stop-success routes in SSEXG and SSRDEX now cache the
stop EXG normalization terms per integral; SSEXG also caches each go EXG
survivor normalization. Both use the model's actual `exg_lb` and `exgS_lb`
values. With two distinct go parameter tuples and nonzero lower bounds, 500
128-node integrations produced identical results before and after. Median
SSEXG time was 0.031 s versus 0.014 s (2.21x); SSRDEX was 0.021 s versus
0.015 s (1.40x). A full raw SSEXG likelihood fixture (542 rows, 271 unique
trials, 32 particles, 64-node stop integration) had the same data checksum
(`81f548b49a7896402c189b5bcf21e98a`) and bit-identical likelihood outputs;
median `calc_ll_oo` time was 0.066 s versus 0.033 s (2.00x). The focused
stop-signal tests passed 91 expectations. See `results/ss_cache_base.log`,
`results/ss_cache_candidate.log`, `results/ss_likelihood_base.log`,
`results/ss_likelihood_candidate.log`, and `results/ss_cache_tests.log`.

The shared-tuple cache now also covers `rexg_logS_at_t()`, which evaluates
truncated race-survival endpoints. A finite-UT two-choice REXG fixture with
12,960 expanded rows and 512 repeated particles produced bit-identical native
likelihoods. Across 15 calls, median time fell from 194 ms to 191 ms (1.016x).
Three paired fixed-transition fits used the same 8x replicated one-subject
Forstmann data, a finite `UT = 2`, and 10 preburn, 10 burn, 10 adapt and 40
sample iterations. Fit objects compared equal, and samples, chain RNG, and
caller RNG states were identical within every pair. Baseline/candidate elapsed
times were 9.399/9.363 s,
10.261/10.154 s, and 9.924/9.840 s. This is a small fit-level gain, about
0.4–1.1% per seed, with no transition or posterior change. The untruncated
single-accumulator control did not enter this endpoint callback and is not used
to estimate its effect. Raw runs and exact outputs are in
`results/step4a_logS_truncated_race/` and
`results/step4a_logS_truncated_fit/`; runners are
`run_rexg_logS_cache.sh` and `run_rexg_logS_fit.sh`.

This completes the measured configurable-bound stop-success paths in
`model_SS_EXG.h`, `model_SS_RDEX.h` and the live raw adapter.

The scalar REXG pdf/cdf callbacks now reuse an exact-key truncation context
across quadrature calls, using the supplied race context when available and a
thread-local fallback otherwise. The focused 8-racer, 512-particle native
likelihood probe was bit-identical and tied the direct-callback control at a
0.587 s median in a reverse-order pair. The fit-level fixture exercised the
finite, nonuniform `UT` scalar route on two racers and eight replicated
subjects. Across seeds 101, 202 and 303, the direct-callback control took
51.836, 55.701 and 55.378 s; the cached version took 38.387, 38.810 and
39.769 s (1.35–1.44x). Every paired fit object and final RNG state was
identical. A reverse-order seed-101 pair confirmed the result: 51.709 s
without the cache and 38.225 s with it, again with identical fit and RNG.
Raw calls and fit records are in
`results/step4a_logS_scalar_cache_vs_noop_race8/`,
`results/step4a_logS_scalar_cache_vs_noop_race8_reverse/`,
`results/step4a_logS_scalar_fit_cache_vs_noop_2race/` and
`results/step4a_logS_scalar_fit_cache_vs_noop_2race_reverse/`.

### Step 4b: stable EXG tails and prepared kernel parameters

The REXG raw likelihood now evaluates the density through a scaled normal-tail
factor, computes the lower CDF from a scaled-tail difference, and forms the
upper CDF from the direct survivor terms. This removes the earlier `z < -8`
Mills-series approximation and avoids cancellation in the zero-truncation
normalizer at large `sigma / tau`. Each exact cached parameter tuple also
retains inverse sigma, the sigma/tau ratio, and log scales, avoiding repeated
division and logarithms for every likelihood row.

The focused EXG and stop-success tests passed 134 expectations with no failures
or warnings. The numerical checks cover the former tail-approximation boundary,
deep left and right tails, a 100:1 sigma/tau ratio, and scales down to 1e-6.
Against the frozen pre-Step-4b build, the 704-row by 256-particle native
likelihood call used identical data and differed by at most 2.27e-12 in log
likelihood; median time fell from 0.014 s to 0.012 s (1.17x).

Three paired default-host EXG fits (seeds 101, 202 and 303) each reached 1,000
draws per chain. Mean elapsed time fell from 12.102 s to 11.603 s (4.1%). All
5%, 50% and 95% quantiles were within four combined MCSEs; the largest
standardized difference was 3.14e-8. Minimum bulk and tail ESS ratios were
exactly 1.0 for each seed, with identical caller and chain RNG states. The
alpha arrays differ only at floating-point scale (maximum absolute difference
3.30e-9). Raw timings, fits, and the posterior-gate report are in the ignored
`results/step4b_final/` directory.

The common one-tuple REXG raw path now reuses prepared density parameters and
skips per-row cache-key matching. Multi-tuple and invalid-row cases retain the
general path. On the same 704-row by 256-particle likelihood fixture, the
likelihood-vector difference from the frozen Step 4b build was at most
4.55e-13; median call time fell from 12.5 ms to 11.0 ms (1.136x). Six paired
default-host fits, with 400 sample draws per chain and seeds 101, 202, 303, 404,
505 and 606, took 10.247/8.767, 9.498/8.995, 9.751/8.755, 9.748/9.284,
9.748/8.996 and 9.747/8.759 seconds (baseline/candidate). The geometric mean
speedup was 1.097x, with a paired-seed 95% t interval of 1.052–1.144x.
Chain and caller RNG states matched. Sample arrays were not bit-identical, but
all paired 5%, 50% and 95% quantiles differed by less than 5.6e-7 combined
MCSEs overall and 4.5e-7 MCSEs within any chain. Across-seed minimum bulk
ESS-per-draw ratio was 0.9998 [0.9993, 1.000]; minimum tail ESS-per-draw was
1.000 in all pairs. Mean user CPU fell from 22.05 to 19.80 seconds; system
CPU, minor faults and peak RSS were effectively unchanged. Raw calls, fits,
resource records and the reproducible runner/analyzer are in
`results/step4b_simd_live/`, `results/step4b_simd_fit_live/`,
`run_step4b_simd_fit_gate.sh` and `analyse_step4b_simd_fit_gate.R`.

A follow-up gathered the common density rows into a separate compute loop and
used GCC 13's vectorization report to verify actual SIMD generation (64-byte
and 32-byte vectors). Kernel outputs differed by at most 4.55e-13, but median
time rose from 11 ms to 12 ms on the 704-row by 256-particle fixture. The
six-seed fit gate also regressed: baseline/candidate elapsed times were
8.993/9.004, 8.997/9.020, 8.743/9.257, 8.750/9.248, 8.746/9.248 and
8.754/9.249 seconds. Geometric speedup was 0.963x [0.935, 0.991]. Quantile
differences remained negligible (largest 2.83e-6 combined MCSE), and bulk and
tail ESS ratios were approximately 1. The gather/compute/scatter candidate was
therefore reverted; the prepared one-tuple path above remains, since its
separate six-seed gate showed a 1.097x fit speedup. In the current source, the
SIMD pragma on that prepared loop does not make GCC vectorize the density or
survivor loop. Raw records and compiler diagnostics are under
`results/step4b_true_simd_{kernel,direct_kernel,fit}/`.

An exploratory prepared log-survivor formula (`log S(t) - log S(0)`) was also
checked against the current stable `ptexg` path. It differed by at most
2.27e-13 in the raw likelihood vector, tied in the 704-row kernel call (11 ms
in both arms), and gave a 1.001x fit speedup across three seeds with a 95% t
interval of 0.997–1.005x. It did not demonstrate a cost reduction; the existing
stable survivor path is retained. The probe records are in
`results/step4b_logsurv_live/` and `results/step4b_logsurv_fit_live/`.

### Step 3: adaptive rebuild dispatch

The rebuild queue now times keyed Laplace refresh/check jobs and estimates a
serial total against an LPT parallel makespan plus observed startup overhead.
Serial execution is selected only when its prediction is at least 10% lower;
otherwise the existing parallel route remains in use. Focused dispatcher tests
cover selection, order, RNG preservation, error propagation and timing
calibration.

Three paired forced-rebuild replays used the same saved three-chain, 40-subject
starting state (fixture RDS MD5 `5311a12528551a748ec0bb76e8564390`), with two
sample blocks and two active chains. Candidate/baseline in-process elapsed
times were 91.920/93.199 s, 93.405/92.692 s and 93.877/94.012 s. The median
paired ratio was 0.999; this is effectively tied. Chain samples and subject
RNG states were identical in all pairs. In the two instrumented, order-reversed
replays, all four recorded job batches stayed parallel: they contained 40, 65
or 120 jobs, and predicted serial work was 2.7–3.9 times the parallel estimate.
The whole-fit fixture therefore validates the no-regression path but does not
exercise serial selection. Candidate wall time, CPU and minor faults were also
within run-to-run variation (raw records in `results/step3_dispatch/`).

A targeted cached EXG refresh now exercises the serial route. After one real
refresh populated the fit/check caches and eight parallel refreshes calibrated
the fork cost, ten paired measurements per seed compared the same three-chain,
one-subject state on each route. Median cached-refresh time was 13.5–16 ms
through forks and 1.5–2 ms serially (6.75–10.67x); the predicted serial cost
was about 5 ms against a 14 ms parallel makespan. Samples, RNG and all compared
continuation fields were identical for seeds 101, 202 and 303. This is a
focused dispatch win of roughly 12–14 ms per cached refresh, not a claim of a
large whole-fit speedup. The large hierarchical fixture remains on the
parallel route, where the measured serial estimate was 2.7–3.9 times higher.
Step 3 passes its no-regression and targeted-benefit gates. Raw timing pairs
are in `results/step3_dispatch_cached/`.

The controlled EXG profile had seven rebuild blocks around 2.43–2.70 s before
warming and 1.78–2.11 s after; its no-rebuild sample block took 1.49 and 1.55 s.
The separate instrumented EXG runs each had 913 likelihood batch calls (12,786
rows), 21 optimizations, 773 objective and 115 gradient evaluations. An
`strace` run counted 87 `clone` calls in each arm; traced startup timing is not
comparable to ordinary timings. `/usr/bin/time` recorded 2.65–2.66 million minor
faults in the three baseline EXG fixed runs and 1.33–1.37 million warmed.
Sampled process-tree peaks were about 977 MiB PSS / 833 MiB private before and
879 MiB / 696 MiB after; 0.1 s sampling misses short-lived processes and bins
are approximate. These patterns are consistent with less child startup and
copy-on-write work, not proof of a causal breakdown.

The analytic two-mode target has exact right-mode probability 0.5. A short
custom-model fit produced identical paired draws and RNG, with 0.463 of draws
in the right mode (indicator MCSE 0.024, bulk ESS 416, rank R-hat 1.008);
this is a coverage probe, not evidence of accurate mode weights. The target is
an analytic test distribution, not a generative model for the supplied data. The custom-model initialization fallback and the one-subject pool also
matched in the paired checks. Hierarchical profiles showed worker growth above
the initial allocation, but the sampled profile is not a proof of peak total
core usage or reclamation at a subsequent block.

A separate expensive-profile run of the hierarchical worker pool covered 60
chain iterations. It recorded 11,748,842 serialized bytes in total (maximum
196,272 in one recorded iteration), including 11,683,080 shared bytes and
43,053 private-message bytes; remaining bytes are shared-file path traffic.
Twelve sampled private-memory readings peaked at 805,711,872 bytes. These
counts describe this short fixture and its profiler's measured transport, not
all serialization in the process tree.

In a separate short hierarchical role-instrumented replay, direct Linux
`smaps_rollup` snapshots identified one master, six block-chain processes,
three templates and sixteen workers. The largest observed RSS / PSS / private
MiB were respectively 229.6 / 224.0 / 223.4 for the master, 402.3 / 252.0 /
188.0 for chains, 229.7 / 214.3 / 210.5 for templates, and 214.1 / 76.4 / 1.6
for workers. Templates and workers were sampled at launch; these are not
worker peak or steady-state measurements. Chain snapshots were at block start
and end. The external tree sampler alone missed reparented workers and
therefore must not be used as a complete role breakdown.

The isolated native likelihood probes reproduced the EXG arithmetic gap at
approximately 150 ns per row-particle for EMC2 versus 33 ns for the archived
reference. The LNR probes measured approximately 56 ns versus 54 ns. These
are performance comparisons: the fixtures are generated by their respective
packages and compress to different row counts, and the EXG distributions
differ in truncation. The archived reference was rebuilt with GCC 13 and the
recorded optimization flags; its OpenMP flags still differ from EMC2's.

The focused worker-pool, core-budget and rebuild-overhead test selection passed
1,402 expectations, with two clean-process-backend tests skipped in the local
backend and one warning emitted by the intentional worker-recycling failure
test. The earlier focused rebuild test passed 68 expectations.

### Step 5: bounded chain-state message prototype

The first continuation-state schema retains `data`, `model` and `prior` in the
chain process and sends the remaining sampler fields and attributes at a block
boundary. On a saved EXG chain, the serialized message was 26,612 bytes versus
845,898 bytes for the full sampler (3.15%). Applying a serialized message to
the retained setup reproduced the complete sampler under `all.equal`. A local
100-round probe took 0.008 s to serialize and 0.015 s to serialize plus
deserialize the message. These are codec feasibility measurements only; chain
process persistence is now available as the opt-in
`options(emc2.persistent_chains = TRUE)` for Unix fits with
`cores_per_chain = 1`, and for higher `cores_per_chain` when all chains stay
active together. A synthetic three-chain, 30-block replay was 3.12x faster
at the median than forking on every block, including process-pool startup and
shutdown, with identical final state.

The end-to-end gate resumed the same saved EXG state for 60 transitions in
three 20-transition blocks, with three active chains. Across seeds 101, 202
and 303, fork-per-block times were 1.108, 1.133 and 1.071 s; persistent-chain
times were 0.755, 0.737 and 0.724 s (mean speedup 1.50x). All samples, RNG
streams, caller RNG and compared continuation fields were identical in every
pair. On the seed-101 resource run, whole-tree minor faults fell from 218,884
to 151,087 (31%), system CPU from 0.93 to 0.60 s (35%), and peak RSS from
262,128 to 256,888 KiB. The focused worker-pool tests passed 277 expectations;
two backend-specific tests were skipped and the existing intentional
recycling-failure test emitted its warning. The feature remains opt-in until
multi-core worker-retention, allocation and failure-recovery gates are
complete; this initial one-core result did not close them. Raw fits and
resource records are under
`results/step5_chain/`; `chain_persistence.R` runs each fixed-transition arm.

The full-concurrency path also supports `cores_per_chain = 2`. On the same
saved EXG fixture, three paired runs with three chains and 60 transitions took
3.136, 3.232 and 3.055 s with block forks, versus 2.891, 3.022 and 2.908 s
with persistent chains (1.05–1.09x per-seed speedup). Samples, chain and caller
RNG states, and compared continuation fields matched exactly. The current
worker-pool test selection passed 287 expectations; two backend-specific tests
were skipped and the intentional worker-recycling failure test emitted its
expected warning. This smaller gain is consistent with subject likelihood
workers still being restarted at block boundaries. At this point, retained
subject-worker, partial-concurrency, allocation and failure-recovery gates
were still outstanding. Raw pairs are in
`results/step5_chain/{base,candidate}_cpc2_{101,202,303}.rds`.

An initial 100-block run without process recycling showed continued private
memory growth: chain PSS rose from about 76–78 MiB initially to 261 MiB for
one chain by block 90, while RSS reached 293 MiB. The process now restarts
after each 25 acknowledged blocks, reconstructing from the returned sampler
state. In the 100-block, 1-transition-per-block check, each chain used four
processes; within each 25-block segment PSS rose from roughly 82–100 MiB to
126–133 MiB and then reset. The final block also triggered the scheduled
proposal refresh, briefly reaching 228 MiB PSS / 261 MiB RSS before the
acknowledged restart. A paired fork-per-block run took 10.408 s versus 7.264 s
with persistent chains and 25-block restarts (1.43x); all three chains had
identical samples, RNG streams and compared continuation state. This tests the
chain process memory bound only; retained subject workers and full allocation
reclamation are still open. The 100-block trace and fits are in
`results/step5_chain/memory_100_restart25.csv` and the sibling RDS/log files.

### Step 5: retain sampling workers across blocks

Persistent chains now retain their forked subject and one-subject likelihood
workers across block boundaries, refresh mutable context with acknowledged
epochs, grow pools when cores are donated, and reclaim excess workers before
the next block. A focused integration test confirms worker PID stability and
context advancement; another exercises context refresh, growth from two to
four workers, and reclamation back to two. The worker-pool selection passed
308 expectations, with two backend-specific skips and the expected warning
from the deliberate recycling-failure test.

Three paired two-core-per-chain runs (60 transitions in three blocks) took a
mean 3.007 s inside `run_emc` with block workers and 1.501 s with retained
workers (2.00x). Samples, chain RNG, caller RNG, and compared continuation
fields matched exactly for all three seeds. Likelihood worker PIDs stayed
constant and their context epochs advanced 0, 1, 2. Whole-process peak RSS was
about 256 MiB baseline and 246 MiB candidate; system CPU rose from a mean 1.19
to 2.06 s and minor faults from 249k to 314k. The elapsed-time win therefore
comes with more system activity on this short two-core fixture. Forcing a
full collection before these two-worker forks increased faults and peak RSS
further, so that extra collection was removed; the existing collection gate
for larger fork batches remains unchanged.

A 60-block, one-transition-per-block paired run with the normal 25-block chain
restart policy took 12.631 s baseline and 8.396 s candidate. Fits remained
identical. Candidate peak RSS was 250 MiB versus 267 MiB, system CPU 10.92 s
versus 13.20 s, and minor faults 969k versus 1,991k. Chain and worker PSS rose
within each 25-block process lifetime and reset when the chain restarted; the
largest chain PSS was 217,196 KiB. With chain restarts disabled for 60 blocks,
candidate peak RSS reached 350,328 KiB, so the restart policy is necessary to
bound this retained-worker path.

Persistent execution now also supports partial multi-core concurrency. Three
paired runs used two active chain slots out of three with two cores per chain;
mean in-`run_emc` time was 9.953 s with block forks and 2.113 s with persistent
chains (4.71x). Samples, RNG, caller RNG, and continuation fields matched
exactly. In the trace, the later chain received four likelihood workers after
the earlier chains completed, then returned to its two-worker baseline at the
next block. This validates the queued-chain donation/reclamation path on this
fixture. The final uneven-cost active-chain donation and live-interrupt gates
are recorded below; worker/task retry, chain failure/kill cleanup, and
retained-worker rebuilds are covered in the following sections.
Raw runs are in
`results/step5_chain/worker_retention_final_pairs/`,
`worker_retention_60blocks/`, `worker_retention_60blocks_restart25/`, and
`worker_retention_partial_cpc2_cfc2/`.

### Step 5: rebuilds on retained subject workers

Persistent chain processes now dispatch Laplace refresh and conditional-check
tasks to their retained subject worker pools. A cache-miss regression was
caught during the first benchmark attempt: the synthetic sampler's alpha store
was a matrix although the mode-start code indexes it in three dimensions. The
worker now constructs the expected `parameter × subject × draw` array, and a
focused regression test exercises a real cache miss. A worker transport
failure test also verifies that the original subject tasks and seeds are
retried once, without committing a partial result.

The fixed-schedule gate resumed the same saved three-chain, 40-subject state
(fixture MD5 `5311a12528551a748ec0bb76e8564390`) for two 100-transition sample
blocks, with two active chains and two cores per chain. Across seeds 101, 202
and 303, baseline elapsed times were 240.744, 242.293 and 242.194 s; retained
rebuild times were 241.535, 241.373 and 241.557 s. The means were 241.744 and
241.488 s (1.001x), an effective tie. Every pair had identical samples, chain
and caller RNG, continuation state after removing measured timing fields, and
allocation/HMC schedule. All six rebuild batches per candidate used retained
subject pools, worker PIDs stayed stable within each chain, and there were no
mode failures. Treat this as a correctness pass with no demonstrated runtime
gain.

At rebuild entry, the seed-303 trace measured combined chain-plus-worker PSS
across two rebuilds per chain. PSS increased by 2.3–8.3 MiB and private memory
by 2.2–9.4 MiB; each chain retained the same worker PIDs. One chain had four
workers after receiving queued-chain capacity. This is a short repeated-rebuild
check, not a long-run memory-growth bound. The final uneven-cost donation and
live-interrupt gates are recorded below. Focused worker-pool tests cover task
retry, chain failure and kill cleanup. Raw fits, timings and the memory trace
are in `results/step5_retained_rebuild_fixed/`.

### Step 5: uneven donation and live interrupt cleanup

The uneven-block integration test starts two persistent chains with two
retained subject workers each. One chain completes first and releases its
allocation; the still-running chain observes the new budget, grows from two to
four workers, then shrinks to its two-worker base before releasing the next
block. Both donated worker PIDs are gone at release, the original workers stay
alive, and both chains enter the following block with two workers. The focused
test passed all nine expectations. The paired partial-concurrency fit above
provides the end-to-end throughput check for this allocation path: persistent
execution with dynamic worker growth was 4.71x faster than block forks on the
three-chain fixture, with identical samples and continuation state. The
uneven-block test itself is a correctness and lifecycle fixture, not a
separate likelihood-throughput estimate.

A live SIGINT was sent while a persistent chain was inside a block with two
nested retained workers. The test used the same on-exit pool cleanup registered
by `run_emc()`; it confirmed the interrupted R process caught the signal and
both the chain and subject-worker PIDs were reaped. The focused test passed all
six expectations. Together with the failure/restart tests and the bounded
restart memory run above, this closes the Step 5 interrupt and allocation
lifecycle gates for the tested fork backend. These checks do not exercise the
separate clean-spawn backend on this Unix configuration.

### Step 2: conditional subject-fit cache

The durable per-chain fit/check cache passed 258 focused expectations with no
failures or warnings; one optional development-validation test was skipped.
This includes exact-conditioning reuse, changed-prior refits, checkpoint/RDS
durability, duplicate stripping, and a configurable fresh-fit cadence.

Full mode reuse was compared against the pre-cache build on seven EXG and
fourteen LNR seeds, with 1,000 draws in each of three chains. Raw data hashes
matched within every pair. The cache intentionally changes chain RNG and
draws; the caller RNG remained identical. No seed had a 5%, 50% or 95%
quantile outside four combined MCSEs. In the combined chains, the largest
standardized quantile difference was 1.02 for EXG and 1.78 for LNR.

For EXG, mean fit elapsed fell from 13.56 s to 11.73 s (13.5%). Across-seed
minimum ESS ratios were 1.025 [0.988, 1.064] for bulk and 1.053 [0.981, 1.131]
for tail ESS; both intervals clear the 0.90 margin. Combined minimum bulk
ESS/draw was 0.931 baseline and 0.920 candidate; tail ESS/draw was 0.919 and
0.921. The EXG cache is enabled by default.

For LNR, mean fit elapsed fell from 103.60 s to 95.20 s (8.1%), but combined
minimum bulk ESS/draw fell from 0.772 to 0.689 and tail ESS/draw from 0.593 to
0.512. Across-seed minimum ESS ratios were 0.919 [0.832, 1.014] for bulk and
0.867 [0.683, 1.102] for tail; neither interval establishes the required
non-inferiority margin. The posterior quantiles agree, but full LNR mode reuse
does not pass the efficiency gate. A seven-seed policy that refit after two
cache hits saved 4.0% elapsed time, with aggregate bulk ESS/draw 0.750 versus
0.722 and tail ESS/draw 0.595 versus 0.514; periodic refitting did not resolve
the tail-ESS loss. A diagnostic-only LNR replay had exactly the same samples
and RNG as full reuse, pointing to cached mode reuse as the cause. Therefore,
default fit reuse is limited to the measured REXG model; other models can opt
in with `options(emc2.subject_fit_cache_max_reuses = Inf)` while their gates
are evaluated. The production-default LNR replay matched baseline samples and
all RNG states exactly (103.98 s baseline, 105.02 s candidate).

A supplementary same-build LNR comparison used the current package with the
model-default fresh-fit policy versus `max_reuses = Inf`, isolating cache
reuse from the surrounding code revision. On three additional seeds with
1,000 draws per chain, mean elapsed time fell from 98.63 s to 89.43 s (9.3%).
All 5%, 50% and 95% quantiles were within four combined MCSEs (largest
standardized difference 1.92). Per-seed minimum bulk ESS/draw was
0.747/0.793, 0.657/0.668 and 0.461/0.765; minimum tail ESS/draw was
0.580/0.638, 0.464/0.541 and 0.237/0.583. The three-seed ESS ratios are
highly variable and do not establish the 0.90 non-inferiority margin. Chain
RNG changes are expected; caller RNG matched in two pairs and differed in one.
The same seeds also reached the alpha ESS 800 precision target at the 400-draw
floor in both arms. Mean time was 126.79 s fresh and 120.00 s with reuse
(5.4% faster); minimum ESS/s ratios were 1.043, 1.051 and 0.928. This pilot
supports keeping LNR reuse available as an opt-in, while leaving the default
unchanged until its across-seed efficiency gate is resolved. Fits and logs are
in `results/step2_lnr_cache_fixed_same_build/` and
`results/step2_lnr_cache_precision_same_build/`.

A separate 14-seed precision comparison tested exact-key optimizer warm starts
while still refitting the Laplace approximation on every refresh
(`emc2.subject_fit_cache_warm_start = TRUE`, with mode reuse disabled in both
arms). Every run reached the alpha ESS 800 target; all 5%, 50% and 95%
quantiles were within four combined MCSEs, with a largest standardized
difference of 2.90. The candidate took a mean 125.72 s versus 128.35 s for
fresh optimizer starts (2.1% lower elapsed time); the paired geometric
speedup was 1.020 [0.998, 1.043], which includes a tie. Minimum ESS-per-draw
ratios were 0.993 [0.892, 1.106] for bulk and 1.047 [0.880, 1.247] for tail.
The lower bounds narrowly miss the 0.90 non-inferiority gate, so this is not
enough evidence to enable the option by default. Keep it available as an
opt-in for further evaluation. This policy intentionally changes chain RNG
streams; caller RNG differed in two of the fourteen pairs. One pair ran in
reverse arm order as a timing-order check. Raw fits are in
`results/step2_lnr_warm_start_precision/`; run
`analyse_lnr_warm_start.R` with the same seed list to reproduce the summary.
Minimum ESS-per-second ratios were 0.982 [0.888, 1.085] for bulk and
1.035 [0.874, 1.226] for tail; these also include no change and miss the 0.90
lower-bound gate, so they do not support promoting the option.

The analytic bimodal cache gate was extended from 400 to 5,000 sample draws per
chain and from one to eighteen paired seeds. Both arms used the same current
package build; `max_reuses = 0` disabled mode reuse and `Inf` enabled it. Each
arm therefore contributes 54 chains and 270,000 draws. The target's exact
right-mode mass is 0.5. Fresh fits estimated 0.4991 (MCSE 0.00168), and cached
fits estimated 0.5009 (MCSE 0.00172); the difference is 0.75 combined MCSE.
Both estimates are within 0.6 MCSE of the exact mass. The 5% and 95% quantile
shifts were 0.85 and 0.25 combined MCSE. The 50% quantile is not informative
for an exactly balanced bimodal target because its median is non-unique; mode
mass is the direct coverage measure here.

Across-seed minimum mode-indicator ESS ratios were 0.956 [0.911, 1.002], theta
bulk ESS ratios 0.979 [0.942, 1.017], and theta tail ESS ratios 0.992
[0.963, 1.022]. The cache/fresh geometric speedup was 1.013x [1.004, 1.021],
about a 1.2% reduction on this cheap custom target. This adversarial fixture
supports opt-in reuse for a pure custom likelihood, but arbitrary custom R
likelihoods can capture mutable state that the cache key cannot describe, so
they remain opt-in. LNR remains off by default because its separate tail-ESS
gate does not pass. Raw fits and timing files are in
`results/multimodal_precision_live/`; `run_multimodal_precision.sh` and
`analyse_multimodal_precision.R` reproduce this gate.

## Frozen reference on common raw data

For seeds 101, 202 and 303, both libraries read the same exported raw data;
all six paired data hashes matched. The archived reference requires validation
attributes absent from this input, so the wrapper supplies them with a 1e-8
boundary tolerance and omits an unused lower-censor attribute. Data columns
are unchanged. EMC2 / reference compression produced 704 / 705 unique EXG
rows and 6,213 / 6,217 unique LNR trials. EXG also has different truncation
semantics. Consequently the following timings and ESS are exploratory
performance references, not equivalent-posterior or non-inferiority tests.

| Fixture, seed order 101/202/303 | EMC2 | Archived reference |
|---|---|---|
| EXG elapsed seconds | 22.298, 22.051, 22.048 | 1.567, 1.596, 1.615 |
| EXG minimum bulk ESS/draw | 0.933, 0.928, 0.907 | 0.364, 0.429, 0.283 |
| EXG minimum bulk ESS/second | 125.6, 126.3, 123.5 | 696.7, 807.2, 525.3 |
| LNR elapsed seconds | 105.448, 106.483, 105.982 | 54.569, 54.586, 55.311 |
| LNR minimum bulk ESS/draw | 0.794, 0.816, 0.803 | 0.464, 0.390, 0.293 |
| LNR minimum bulk ESS/second | 22.58, 22.99, 22.74 | 25.51, 21.45, 15.88 |

Per-parameter quantiles, Monte Carlo errors, bulk and tail ESS, R-hat, and
minimum tail ESS per draw and second are in
`results_optimized/shared_comparison.txt`. Three seeds are descriptive and do
not establish uncertainty intervals for an efficiency or non-inferiority
claim. Frozen reference build logs, saved fits, timing records and shared raw
inputs are retained in the ignored `results_optimized/` directory. See the
protocol for compiler and legacy-validation qualifications.

## Step 6a: block-scoped likelihood cache validation

`run_stage()` now prepares the likelihood partition once at block entry and
passes a validated-context flag through particle, marginal, pooled, ensemble,
HMC and retained-worker likelihood routes. Public likelihood calls still
validate on entry, and custom R likelihoods keep that path inside the manager.
All-finite data are rebuilt at each block entry because their fast cache form
cannot detect response-column mutation; censored/non-finite data validate the
cached partition once at that boundary. Focused tests cover stale partitions,
list attributes, registered/custom dispatch and worker context propagation.

The controlled seed-101 profiles (`HOST=0`, sampler profiler enabled) took
9.795 s baseline and 9.112 s candidate for EXG, and 135.238 s and 126.310 s
for LNR. Both pairs had identical samples and RNG states. The end-to-end gate
used the frozen Step 5 install at `a10dc96d`, default fit host, fixed
transitions (50 preburn, 200 burn, 200 adapt and 400 sample iterations), data
generation seed 123, sampler seeds 101/202/303, CPUs 0–26, and single-threaded
OMP/OpenBLAS/MKL. All six pairs had identical data checksums, samples, chain
and caller RNG, and continuation fields.

| Model | Baseline elapsed (s), seeds 101/202/303 | Candidate elapsed (s) | Mean speedup |
|---|---|---|---:|
| EXG | 9.504, 10.009, 9.250 | 9.250, 9.248, 8.999 | 1.046x |
| LNR | 129.880, 130.651, 130.175 | 128.382, 125.370, 125.845 | 1.029x |

A same-seed reverse-order pair also favored the candidate: EXG took 9.003 s
candidate versus 9.502 s baseline, and LNR took 127.111 s versus 129.857 s.
Those reverse-order fits had identical samples and RNG. The default-host
speedups are smaller than in the controlled profiles, so the measured fit
gains are 4.6% for EXG and 2.9% for LNR. Mean system CPU, minor faults and peak
RSS were 3.75 s / 905k / 319 MiB baseline versus 3.92 s / 950k / 319 MiB for
EXG, and 9.50 s / 2.72m / 364 MiB versus 8.12 s / 2.32m / 356 MiB for LNR.
The EXG resource differences are small and mixed; LNR improved on these three
measures. Raw profiles, fits and `/usr/bin/time -v` records are in
`results/step6a_profiles/`, `results/step6a_end_to_end/` and
`results/step6a_reverse_order/`.

An Rprof follow-up found that the Laplace optimizer and its diagnostic probes
still validated the likelihood partition on every `.emc_burn_ll_batch()` call.
Mode starts now prepare data once before their parallel fits; refresh and
rebuild contexts pass that validated state through optimization and diagnostic
calls. Standalone batches default to validation, and custom R likelihoods
continue to validate inside `calc_ll_manager()`. The focused cache and mode
start tests passed 92 expectations. Fixed-transition EXG and LNR profiles
used seed 101, the same generated data, three serial chains, `HOST=0`, and one
BLAS/OpenMP thread. Samples and chain RNG states matched exactly in all six
chains. EXG took 17.364 s before and 17.253 s after. LNR took 359.606 s before
and 361.029 s after; this single timing pair shows no end-to-end speedup.
However, the prior LNR profile attributed 0.94 s self time to
`.is_valid_ll_cache()` (1.93 s total) during repeated optimizer calls; the
candidate profile measured 0.01 s self time (0.04 s total), with remaining
samples at data-preparation boundaries. That removes about 1.9 s of sampled
validation work, but not enough to move this likelihood-dominated whole-fit
timing. Raw profiles, fits and resource records are in `results/step6c_rprof/` and
`results/step6c_laplace_validation/`.

## Step 6b: likelihood-context cost gate

The prologue probe now materializes the same natural-scale parameter buffer
as `calc_ll_oo` before returning, so its timing includes mapping, transforms,
bounds and materialization. On the paper EXG/LNR data, 256 jittered particles
gave bit-identical likelihoods through `calc_ll_manager()` and direct
`calc_ll_oo()`, and every row passed the bound check. The manager's block-
validated path represents Step 6a's internal sampler call; the default manager
path represents public validation.

| Model | DADM rows | Public manager (ms) | Block-validated manager (ms) | Native call (ms) | Full prologue (ms) | Prologue share |
|---|---:|---:|---:|---:|---:|---:|
| EXG | 704 | 12.4 | 12.2 | 12.4 | 0.42 | 3.4% |
| LNR | 18,639 | 274.8 | 265.0 | 262.6 | 13.47 | 5.1% |

Scaling the particle count from 1 to 256 separated the per-particle work from
the per-call cost. Linear fits to the prologue timings had intercepts of
0.067 ms for EXG and 1.083 ms for LNR, with slopes of 0.00137 and 0.04855 ms
per particle. The intercept is an upper bound on what a retained native plan
could save: creating a fresh mapper for the first particle and some per-call
setup also contribute. At 256 particles that upper bound is about 0.5% of the
LNR native likelihood call; most of its 13.47 ms prologue scales with particle
count. The block-validated R manager added 2.4 ms over the direct native call
in LNR, so a context may still save roughly another percent if it avoids that
dispatch work. This makes a narrowly scoped context prototype worth measuring,
but does not support assuming that all prologue work can be removed.

The run used the frozen Step 6b probe install, seed 123 data, seed 101 particle
clouds, CPU affinity 0–26 and one OMP/OpenBLAS/MKL thread. Raw logs, RDS timing
vectors, and `/usr/bin/time -v` records are in
`results/step6b_prologue_scaling/`; the reproducible runner is
`run_step6b_prologue.sh`.

The retained-mapper prototype was then tested as an explicit, process-local
context. It caches the parameter-table layout and transform plan; particle
values, data, and likelihood scratch stay call-local. The manager only selects
it for validated internal calls to registered native models. Its cache key
covers trial count, model type, proposal columns, parameter types, constant
values, designs, bounds and transforms. A PID change rebuilds the context, so
forked workers and checkpoint restarts do not reuse native pointers. Custom
R likelihoods and trend models use the existing path. The focused tests passed
12 expectations, including bit-identical likelihoods after changing proposal
values and rejection of stale proposal columns, trial counts, and constants.

In a repeated prologue probe, the 256-particle EXG native call was 12.0 ms both
with and without the context. LNR was 265.8 ms without it and 263.8 ms with it
(1.008x); at 128 particles the medians both rounded to 133.0 ms. Likelihoods
were bit-identical. The manager medians were 12.2/12.0 ms (public/validated)
for EXG and 268.2/264.2 ms for LNR. This points to a small possible LNR gain,
not a meaningful EXG gain.

The fixed-transition fit gate used seed 101, three serial chains, `HOST=0`,
identical generated data within each pair, and one BLAS/OpenMP thread. Context
off/on elapsed times were 7.755/7.873 s for EXG and 123.485/123.622 s for LNR.
All six chains had identical samples and RNG states, and caller RNG states
matched. One timing pair per model is descriptive; it does not establish a
small speed difference. The context remains available with
`options(emc2.native_pt_context = TRUE)`, but is off by default because these
fit timings show no end-to-end gain. Raw prologue and fit records are in
`results/step6b_context_prologue_live/` and
`results/step6b_context_fit_live/`; runners and the fit-pair analyzer are
`run_step6b_prologue.sh`, `run_step6b_context_fit_gate.sh` and
`analyse_step6b_context_fit_gate.R`.

## Step 6c: particle-step materiality gate

A fresh fixed-transition EXG run on the current code recorded 2,550 iteration
rows. Median particle-phase time was 6 ms in burn/adapt/sample and 11 ms in
preburn; these phase timings include likelihood work. The latest available
worker-side function profiles for this sampler path (seed 101, before the
context prototype) report `new_particle` self time of 0.52 s in EXG and 0.99 s
in LNR, about 0.20 and 0.39 ms per iteration. Its inclusive time (14.97 s
EXG, 332.51 s LNR) is dominated by likelihood evaluation (`calc_ll_oo`:
13.13 s and 349.18 s). The context prototype does not change `new_particle`.
The plan's materiality marker is about 1.5 ms of R overhead for 15 ms of
likelihood work. The fresh iteration records are in
`results/step6c_particle_gate_live/`; the function profiles are in
`results/step6c_rprof/`.

A corrected one-chain in-process profile, which does not lose forked children,
recorded 14.49 s inclusive in `new_particle` over 2,550 transitions, including
13.08 s in the likelihood manager. This leaves about 1.41 s (0.55 ms per
transition) in the enclosing particle phase beyond likelihood evaluation,
roughly 8% of the full fit. That is below the plan's 1.5 ms materiality marker,
but close enough to test a narrow R-to-C++ weight calculation before ruling
out a broader port.

The first opt-in C++ prototype covers only the common non-marginalized,
single-component weight calculation. Proposal generation, likelihood calls,
selection RNG, adaptation and EMC2's independent/local composition remain in
R. Twelve focused seeded transitions matched proposal state and RNG exactly;
86 particle, adaptation and transition expectations passed. On three EXG
seeds, full fits took 7.865/7.709, 7.710/7.873 and 7.893/7.801 seconds
(R/C++), for a 1.004x geometric mean. One LNR seed took 124.687/123.596
seconds (1.009x). All draw arrays, stages, likelihood samples and RNG states
matched exactly; the analyzer ignores only `pm_settings` metadata, which
records the different option value. This slice does not show a useful fit
speedup. It remains opt-in as a measured prototype while the rest of the
particle-step overhead is considered; these results do not decide the value
of a complete transition-preserving port. Raw fits and timings are in
`results/step6c_particle_fit_live/`, with
`run_step6c_particle_weights_gate.sh` and
`analyse_step6c_particle_weights_gate.R`.

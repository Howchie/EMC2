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

Seven paired seeds used the same pooled LBA data (hash
`d3650406e44a8b6307adab451250714a`) and 400 production draws in each of three
chains. Five pairs had every parameter's 5%, 50% and 95% quantiles within four
combined MCSEs; seed 404 had one exceedance and seed 606 had five, with the
largest at 5.99. Combining the 21 independent chains per arm, all 75 quantile
comparisons were within four combined MCSEs (largest 2.07). This supports
posterior agreement at the aggregate level, though some pooled coordinates
remain weakly mixed (maximum rank R-hat 1.14 baseline, 1.18 candidate).

Mean elapsed time was 94.28 s baseline and 91.19 s candidate, a 3.28% paired
reduction. Aggregate minimum bulk ESS per draw was 0.01198 baseline and 0.00905
candidate; minimum tail ESS per draw was 0.03317 and 0.01964. The corresponding
minimum bulk ESS/s was 0.152 and 0.119, and minimum tail ESS/s was 0.422 and
0.259. Across-seed
candidate/baseline minimum ESS ratios were highly variable: the geometric mean
and 95% t interval were 1.11 [0.48, 2.58] for bulk and 1.19 [0.40, 3.52] for
tail ESS. The interval is too wide to establish the proposed 10% ESS
non-inferiority margin, and the pooled run is too short to make a precise
minimum-ESS comparison. Treat the time gain as measured, posterior quantiles
as consistent in aggregate, and pooled efficiency as inconclusive. Full
per-seed summaries and aggregate diagnostics are in
`results/pooled_posterior_comparison.txt` and
`results/pooled_aggregate_comparison.txt`.

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

This completes the measured configurable-bound stop-success paths in
`model_SS_EXG.h`, `model_SS_RDEX.h` and the live raw adapter. Scalar contexts
that supply a stable parameter tuple remain to be measured before claiming the
full plan's Step 4a.

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
within run-to-run variation (raw records in
`results/step3_dispatch/`). A focused small/uneven likelihood workload is still
needed to establish a benefit before calling Step 3 complete.

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

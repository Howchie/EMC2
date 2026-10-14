# Response to the WS2 independent review

Reviewed against `ws2/rebuild-overhead` at `74dccb9b` and its source in
`/data/work/EMC2_ws2`. The review is useful, particularly in identifying the
work still needed to approach the original timing target. Several proposed
shortcuts are not redundant in the implementation as written, however, and
would change proposal decisions, random-number use, or the fitted atlas.

## Scope and measured result

The review is right that the under-10-second rebuild target was not reached,
and that a complete burn-plus-adapt wall-time comparison is absent. The report
does not claim either result. The original plan's 27-core design was superseded
for further measurements by the user's explicit instruction to use nine cores
(3 chains × 3 cores) and accept relative comparisons. There is no 27-core
measurement or extrapolation in the report. The LBA saved-state rebuild was
79.154 to 65.147 seconds at nine cores; this is not a whole-fit saving.

The small-model requirement also deserves a direct qualification: the
Forstmann saved-history sample-next route was about 1% slower in three short
profiles, while its burn and sample-entry routes were faster. The sample-next
profile advances a saved fit to its next rebuild gate; it does not simulate
the transitions before that gate. This is a small observed regression on that
route, not evidence that all Forstmann runs regress. Its cause is not isolated.

## Findings

| Review finding | Assessment |
|---|---|
| Duplicate Pareto importance-sampling checks | There are two similar checks, but they are not an exact duplicate that can be skipped without changing behavior. |
| Laplace refresh is obsolete in sample stage | Incorrect as stated; the refreshed likelihood approximation is still consumed by sample proposals. A different refresh cadence would change proposals. |
| Reuse the persistent sampling pool | Plausible architectural work, not a drop-in change. The existing pool serves a specific sampling protocol. The claimed memory cause is unproven. |
| Orbit scheduling and inner subject work | Credible follow-up to measure. Cost ordering is absent at this call site, though the claim queue already assigns jobs dynamically. The proposed savings are not yet evidenced. |
| Prune nearby atlas nodes | Research idea that changes the regression data and therefore the fitted proposal surface; no exactness guarantee is supplied. |
| Forward differences or analytic gradients | Forward differences change numerical optimization accuracy and results. Analytic gradients need model-specific support and a fallback for custom likelihoods. |
| Single-pass C++ stopping diagnostics | Worth profiling as a later implementation, but the proposed replacement must preserve the package's existing diagnostic versions and inputs. The sub-50-ms estimate is unsupported. |
| Adaptive Laplace-refresh cadence | Changes when fits and proposals are refreshed. The workstream explicitly said not to implement changes to adaptation timing or gates. |

### 1. Importance-sampling checks

`.emc_laplace_refresh()` checks the optimized mode with
`.emc_start_cov(result)`, then records the summary in `burn_mode_stats`.
`.emc_check_likelihood_fits()` rebuilds a conditional Gaussian from the
canonical likelihood quadratic and current group prior, checks that fit, and
stores the result in `lik_laplace[[s]]$check`. The latter check gates whether
`.emc_attach_lik()` supplies that approximation to the efficient proposal.

For an exact quadratic with an unregularized Hessian and a fully converged
mode, the two Gaussian inputs can coincide. The code does not guarantee those
conditions: the canonical Hessian is projected to positive semidefinite,
regularization and optimizer tolerances can differ, and the checks draw
independent 256-particle samples. The first check also advances the subject RNG
stream before the second. Reusing the first result would therefore change the
second check's Monte Carlo estimate, its possible proposal gate, and the RNG
trajectory. The review's possible CPU saving is worth investigating only with
an explicitly accepted change in those semantics; it is not a zero-risk
deduplication under WS2's identity requirement.

### 2. Sample-stage refresh

The review correctly observes that sample-stage `create_chain_proposals()`
does not use `burn_mode_var` for its empirical chain covariance. That does not
make the full refresh unused. The refresh updates `lik_laplace`; sample
proposals attach that fit and use `.emc_laplace_conditional()` to form the
efficient conditional proposal when its check passes. The atlas also consumes
the likelihood fits. Skipping refresh or retaining an older fit could be
studied as a change to the proposal adaptation, but it does not preserve the
current proposal path by inspection.

### 3. Worker-pool reuse and memory

The work queue forks only on Unix when a job batch is larger than twice its
core width; otherwise it takes the existing apply path. The persistent pool
created by `run_emc()` is used by the sampling protocol, with chain slots
assigned for sampling after proposal rebuilding. Reusing it for arbitrary
mode optimizations and atlas jobs would require extending that protocol and
handling model data, errors, RNG streams, and worker state across blocks. It
may reduce process startup, but it could also retain data or closures in the
long-lived workers, the lifecycle problem WS2 was designed to avoid.

The measured master RSS increase during the profiled 11-block run is real, but
the run does not establish fork/join heap fragmentation as its cause. Profiling
allocations, ordinary R heap retention, and copy-on-write accounting are other
possibilities. The stable template RSS establishes worker stability for that
run; it does not establish bounded master RSS. A pool migration should be
supported by a separate memory and lifetime experiment.

### 4. Orbit work

This is the clearest follow-up candidate. Orbit-width jobs go through the
dynamic claim queue, but `.emc_shear_orbit_steps()` does not pass a predicted
cost order. Each job's width search can evaluate the orbit up to four times,
and `.emc_shear_orbit_terms()` loops through subjects serially inside that
outer job. A cost order may shorten the final worker tail. However, the queue
already hands out new jobs as workers finish, and the reported idle-core
seconds also include contention and unmeasured overhead; neither fact alone
shows that 30.6 core-seconds are recoverable by sorting. Parallelizing the
inner loop while retaining outer workers also risks oversubscribing the core
budget. A nine-core per-job profile should first establish variation in job
cost and where time is spent. A batched likelihood API would need to handle
the per-subject data/model interface rather than assume all subjects form one
existing `calc_ll_manager()` call.

### 5–8. Proposal-surface and numerical changes

The atlas already selects subjects using orbit contribution spread. Filtering
nearby nodes is a different operation: those nodes become training rows for a
curve fit, so dropping them can alter the fit and proposal. It needs a
criterion that preserves the fitted surface, or must be treated as a sampler
change.

Central finite differences currently use symmetric evaluations. Switching to
forward differences reduces evaluations but changes truncation error and can
change optimizer stopping points and Hessians. Analytic gradients could be a
strong optimization for supported models, with a general fallback, but the
review's CPU estimate requires implementation and measurement.

The stopping code supports both traditional split-chain R-hat with its
transform and multivariate options, and rank-normalized/folded R-hat; ESS also
has old and new paths. A C++ kernel over raw `alpha` alone would not cover
`mu`, `sigma2`, parameter selections, transforms, and both versions with the
same outputs. The present parallel implementation deliberately keeps each
existing diagnostic group and function intact. A compiled replacement is a
reasonable later project after specifying and testing parity for each path.

Finally, a group-displacement refresh gate changes adaptation timing. It may
be a worthwhile algorithmic design, but it needs explicit conditions for
stale-fit reuse and production-stage adaptation; it was correctly left as a
future option under the workstream's boundary.

## Follow-up priority

1. Profile individual orbit-width jobs on the authorized nine-core setup; try
   a deterministic predicted-cost ordering only if runtimes are uneven, then
   compare wall time and exact output/RNG state.
2. Measure a small-workload serial/parallel cutoff to address the observed
   Forstmann sample-next slowdown without penalizing larger workloads.
3. Treat check reuse, sample-stage refresh cadence, atlas node pruning,
   gradient changes, and worker-pool migration as separate proposals with
   explicit sampler-semantic and memory evidence.
4. Keep the target status as missed until a comparable nine-core result shows
   otherwise. Do not infer a 27-core result from the nine-core profiles.

The source-level assessments above do not change the implementation or the
reported measurements. No new benchmark or test was run for this response.

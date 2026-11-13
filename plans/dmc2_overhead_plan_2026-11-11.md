# Closing the fixed-overhead gap to DMC2 (2026-11-11)

Keep EMC2's sampler architecture and model definitions while removing repeated
process startup, likelihood setup and rebuild work. DMC2 (`/data/work/DMC2`)
provides references for kernel arithmetic and persistent likelihood setup.
EMC2 retains its proposal schedule, adaptation, group moves, marginalisation,
shear, HMC, collapsed moves and allocator.

Source baseline: EMC2 HEAD `7e161f75`, DMC2 HEAD `fe8341a5`. Recheck the named
functions against the implementation baseline when starting each step.

## Measurements and starting assumptions

Fixture: `CensoringTruncation/paper.R` ex-Gaussian, 10% upper censoring,
`type = "single"`, 704 unique rows, 3 chains, `cores_per_chain = 1`.
Reported fit time: EMC2 22.5 s (user 49 s, sys 14–16 s), DMC2 1.8 s (sys ~0).
The two test fits reported EMC2 ESS per draw at approximately 1.45 times
DMC2's; establish reproducible EMC2 baselines before using this as a target.

| Item | Starting measurement |
|---|---|
| Blocks in the fit | 15; 9 rebuild |
| Block without rebuild | 0.82 s; the chain takes 0.53 s in-process for the same 100 iterations |
| Block with rebuild | 1.8–2.6 s; `.emc_laplace_refresh` ~0.49 s and `add_proposals` ~0.47 s, including `.emc_check_likelihood_fits` |
| Rebuild work in-process | Mode and Hessian 0.037 s per chain; importance check ~0.03 s for 256 candidate draws |
| PSIS startup | First `loo::psis` call ~0.28 s; subsequent calls ~0.4 ms |
| Forks | 107 `clone` calls; R heap ~130 MB, RSS 235 MB; full GC 0.16–0.19 s |
| REXG likelihood | EMC2 155 ns per row-particle; DMC2 39 ns |
| LNR likelihood | EMC2 59 ns per row-particle; DMC2 66 ns |
| In-chain R overhead | ~1.2–2 ms per subject step, including ~0.4 ms for likelihood-data cache handling |

These fit and kernel measurements need reproduction in Step 0. The original
fixture does not exercise hierarchical fits, worker pools or core reallocation.
The system-time pattern is consistent with copy-on-write faults during child
GC; confirm that attribution with minor-fault and process-memory measurements.

A separate startup probe reproduced cold forked PSIS calls of 0.287–0.298 s
and 0.001 s after warming the parent, with no change to `.Random.seed`.

EMC2's REXG distribution truncates at zero and renormalises. DMC2's race
kernel uses the untruncated distribution. Use DMC2 for performance references;
compare before/after EMC2 posteriors under the same EMC2 model.

## Step 0. Establish the benchmark harness and acceptance gates

Fixtures:

- The ex-Gaussian fit above and paper.R LNR `datTC`.
- An uneven-cost hierarchical fit with `cores_per_chain > 1` and enough subjects
  to make donation of idle cores observable.
- A pooled-parameter fit and a difficult or multimodal fit for cache validation.

Exercise full and partial chain concurrency, the single-subject likelihood
pool, and a custom-model fallback. Use small focused fixtures where a full
fit is unnecessary to test transport, invalidation or cleanup.

Record wall time, user/system CPU across the process tree, block profiles,
likelihood and optimiser evaluation counts, fork counts, minor faults,
serialized bytes, and memory by chain/template/worker. Pin package versions,
compiler flags, BLAS threads, data and seeds. Run controlled profiles with
`options(emc2.sampler_profile = TRUE, emc2.fit_host = FALSE)` and end-to-end
benchmarks with the default fit host enabled.

Measure both fixed-transition throughput and time to a common precision target.
Use multiple seeds and repeated timing runs; define Monte Carlo error tolerances
and an ESS-per-draw non-inferiority margin before comparing variants. Record
posterior quantiles, their Monte Carlo uncertainty, minimum ESS per draw and
ESS per second. Estimate uncertainty across runs rather than requiring each
noisy minimum ESS estimate to increase.

Acceptance gates:

- Steps 1, 3, 4a, 5 and setup-only changes in 6 preserve draws, RNG states and
  continuation-relevant sampler state against their immediate baseline.
- Step 2 is gated at the posterior level. Skipping diagnostic draws changes
  subsequent subject streams; that change is accepted. Check posterior agreement
  within the declared Monte Carlo tolerances and ESS-per-draw non-inferiority.
- Step 4b meets explicit numerical error limits against independent references,
  followed by the posterior and efficiency gates.
- Every step demonstrates its intended cost reduction without a material
  regression on the other relevant fixtures. Record remaining bottlenecks
  before proceeding.

Implement in this order: 0, 1, 3, 2, 4a, 4b, 5, then whichever parts of 6 the
new profiles justify. Step numbers identify the work areas; Step 3's serial
route provides a cheap rebuild baseline before adding the cache.

## Step 1. Warm PSIS in the fork parent

Run one deterministic, nonconstant finite-weight `loo::psis` call once in each
process that will fork rebuild jobs. At present this is the caller executing
`run_emc`, including the fit host. When Step 5 introduces rebuild jobs on
workers, warm any independent template that forks those workers as well.
Warming an interactive caller does not warm a separately launched Rscript.

Use the existing `loo::psis`, preserving RNG state and diagnostic behavior.
Verify that children inherit the warmed code and that the startup cost is paid
once per parent. Measure the reduction in rebuild-block time before estimating
further gains.

Entry points: `.emc_laplace_is_start` and `run_emc` in `R/mode_start.R` and
`R/fitting.R`; process launch in `R/fit_host.R` and `R/chain_pool.R`.

## Step 2. Cache successful conditional fits and reuse matching checks

Add a durable, chain-local subject-fit cache. Its conditioning key includes
subject prior mean/covariance, pooled values after `.emc_pooled_sync`, the
identity or version of likelihood data/model configuration, and the
proposal-prior scale. Construct/version immutable configuration at setup
boundaries rather than hashing full data on every likelihood call.

Store the conditional mode (`result$alpha`), widened proposal centre
(`result$proposal_mode`), proposal and conditional covariances, likelihood
quadratic, diagnostic inputs/results and cache key. Populate the existing
one-use `burn_mode_mu`/`burn_mode_var` fields from this cache at refresh time.
Those fields are consumed and cleared by `create_chain_proposals`; they are
not durable cache storage. Include the cache in checkpoints and define a
cache-miss path for older saved fits.

For `type = "single"` with unchanged likelihood inputs, reuse a successful fit.
For hierarchical or pooled fits, refit whenever the conditioning key changes.
Do not cache failures indefinitely: retry at the next scheduled refresh, and
invalidate entries when their inputs change.

Caching a successful mode deliberately replaces repeated optimisation from
the current draw. An unchanged objective can have several local optima, so
validate mode coverage and efficiency on the difficult fixture. If first-fit
reuse fails those gates, retain periodic fresh-start searches for that class
of fit and measure the resulting cost. For necessary refits, evaluate warm
starts from the cached conditional mode; keep them only if they save work
without losing coverage or efficiency.

Reuse an importance check only when the actual assessed Gaussian mean and
covariance and its target conditioning state match a stored check. Compare
these after curvature projection, radius limiting, covariance shrinking,
eigenvalue flooring and regularisation. Use exact equality initially; retain
separate checks when equivalence cannot be established. Refresh checks when
their inputs change. Propagate reused diagnostics consistently to proposal
fallback, shear and collapsed planning.

The checks advance each subject's stream. Reusing or skipping one changes
later draws and the noise in diagnostic-based decisions. Apply the posterior
and efficiency gates from Step 0, including checkpoint/resume consistency for
the new implementation.

Entry points: `.emc_burn_subject_prior`, `.emc_burn_mode_start`,
`.emc_laplace_refresh`, `.emc_check_likelihood_fits`, `.emc_lik_canonical`
and `.emc_hessian_covariance` in `R/mode_start.R`; proposal construction in
`R/fitting.R`.

## Step 3. Dispatch cheap rebuild jobs in-process

Extend `.emc_rebuild_lapply` (`R/rebuild_queue.R`) to choose serial execution
when it is cheaper than fresh forks. Estimate serial elapsed cost from recent
likelihood timings and observed evaluation counts, including optimiser,
Hessian and diagnostic work. Compare it with estimated parallel makespan
plus fork/startup overhead, accounting for unequal subject costs.

Calibrate fork overhead once per fit and update job-cost estimates as work is
observed. Include calibration cost in end-to-end timing and retain a simple
fallback when estimates are unavailable. Recalibrate the dispatch decision
after PSIS warming and cache reuse change the workload.

Preserve job ordering, explicit subject streams, error propagation and the
caller's RNG state across either route. Test serial/fork identity and recovery
from a failed job. Reuse of persistent sampling workers is implemented in
Step 5, when those workers actually survive the block boundary.

## Step 4. Reduce REXG kernel cost while preserving the distribution

### 4a. Cache repeated normalization calculations

In `src/model_EXG.cpp`, `drexg_raw`, `prexg_raw`, `rexg_logS_at_t` and scalar
pdf/cdf routes repeatedly evaluate parameter-only truncation quantities.
Compute them once per reusable parameter tuple within a particle/evaluation.

Expose the required `ParamTable` constant/cell metadata through the model
context or an explicit preparation interface; raw kernels currently receive
column pointers rather than the table. Key reuse on the joint mu/sigma/tau
and truncation-bound tuple. Respect masks, reset particle-dependent caches,
and fall back to per-row evaluation when metadata is unavailable, including
trend paths where constant tracking is disabled.

Preserve each existing density/CDF calculation's arithmetic and branch choices.
Cache natural and log quantities separately where needed for draw identity.
Cover `model_SS_EXG.h` and `model_SS_RDEX.h` using their actual configurable
lower bounds. Reuse scalar preparation where its caller supplies a stable
parameter context.

Verify bit identity of kernel outputs and fit draws before measuring the gain.

### 4b. Fuse the common kernel path with stable tail fallbacks

Use DMC2's `dexg_pexg_fast_impl` in `src/model_exgaussian.cpp` as the reference
for gather, SIMD compute and scatter, with log output selected at compile time.
Adapt `ExGaussian::weighted_tail`, `raw_pdf`, `raw_prob` and their
`pnorm_utils.h` helpers. Keep the 1e-12 scale clamp and disable optional normal
or shifted-exponential distribution approximations. EMC2 already builds with
`-fopenmp-simd`.

For positive effective time t, preserve EMC2's normalized density f(t)/S(0),
survivor S(t)/S(0), and CDF (F(t)-F(0))/S(0), together with the existing support,
mask, invalid-parameter and likelihood-floor rules. The raw race API returns
log values; include all conversions and normalization in the benchmark.

Retain stable log evaluation for extreme densities, upper-tail survivors,
cancellation near the truncation bound, and small S(0). An unnormalized
natural value may underflow while its normalized value remains representable;
use log ratios in that case. The existing `z < -8` density branch is one
fallback, not the full tail policy.

Validate density normalization, CDF/survivor agreement, support boundaries,
tiny sigma/tau, extreme scale ratios and far tails. Compare against the old
kernel and independent high-accuracy references; account for the old
far-left density approximation when setting tolerances. Then run the fit-level
gates. Measure how much of the kernel gap remains after EMC2's normalization
and log-output requirements are included.

## Step 5. Keep chain processes and their worker pools across blocks

Introduce one persistent process per chain, driven by the length-prefixed FIFO
transport in `R/chain_pool.R`. Start with a chain-state round-trip prototype,
then add persistent worker ownership and allocation messages. Activate at most
`cores_for_chains` chains at a time, even when more chain processes exist.

The current `.emc_wpool_session_*` machinery retains clean templates; sampling
workers stop at `run_stage` exit. Extend that lifecycle explicitly so workers
can serve both sampling and rebuild phases.

### Block state and memory

The master owns accumulated history, proposal rebuilding and stopping logic.
A chain process holds immutable setup, current mutable sampler state, the
block's proposal settings, and at most one block of draws. At completion it
returns that block plus every mutable continuation field needed by the master:
RNG streams, adaptation/allocator state, counters and diagnostics. After
acknowledgement, release the block and retain only current state.

Send proposal summaries/plans and current-state updates back to the chain;
never send accumulated history. Preserve any bounded working summaries needed
by the transition. Define the request/result schema explicitly and verify all
continuation fields across stage changes, thinning, trimming and save/resume.

DMC2's `SamplerContexts` and `sample_chains` in `src/sampler/runner.cpp`
illustrate persistent likelihood setup with per-call chain state. Its
`finish()` returns updated sampler state, and its store clones the supplied
sample object. Bound the input and retained history in EMC2's own protocol.

Keep `emc2.worker_recycle` for subject/likelihood workers. Add an explicit
chain-process restart policy, restoring from the acknowledged continuation
state, if recycling is needed to bound chain memory. Over a fixed-size
50-block run, measure chain/template/worker RSS and private memory after
warmup and require no sustained growth attributable to retained history/jobs.

### Core allocation and reclamation

Replace the fork-inherited allocation descriptor with control messages carrying
a block/allocation epoch, chain/job identity and current budget. Chains report
completion and capacity release; the controller donates idle capacity to
unfinished chains. Deliver updates at safe iteration boundaries, preserving
the current dynamic reallocation behavior of `.emc_core_ctl`/`.emc_cores_now`.
Reject stale replies or releases from earlier epochs and keep control messages
serviceable while a block runs.

Add pool shrinking/release as well as growth. Reclaim donated workers before
the next block activates all eligible chains. Count subject workers,
single-subject likelihood workers and inner BLAS/likelihood threads within
the budget, including temporary capacity during recycling. Preserve partial
chain concurrency and handle failed growth or worker replacement without
blocking progress.

The acceptance fixture must show donation during an uneven block, reclamation
at the next block, no budget excess, and no material throughput regression.
Compare draws and continuation state across worker counts and allocation
changes. Include worker/chain failure, restart, interrupt and parent-exit cleanup.
These checks gate rollout of persistent chains.

### Rebuilds on retained workers

Add a rebuild message handler carrying subject identity, conditioning prior,
pooled values, optimiser start, cache/version information and subject seed.
Return the fit/check and updated seed in deterministic subject order. Construct
closures only within the handler and release job-local objects after replying.
Keep immutable data/setup and the intentional bounded caches.

Schedule rebuild and sampling work within the same core budget. On a transport
failure, retry from the original input seed and commit each logical result
once. Verify result/state identity against serial dispatch and measure memory
across repeated rebuilds. Warm PSIS in the appropriate worker template.

## Step 6. Reduce remaining in-chain overhead where profiles justify it

### 6a. Validate immutable likelihood setup at block boundaries

Move repeated `.cache_ll_data_attrs`/`.is_valid_ll_cache` work in `R/utils.R`
to block/context preparation for internal sampler calls. Retain validation
for public likelihood entry points. Define invalidation for data, pooled
values, designs, trends and model-configuration updates; refresh mutable
inputs without silently reusing stale partitions or setup.

### 6b. Retain a native likelihood context

Build designs, bounds, transforms and column maps once per chain/context and
pass C++ a managed context pointer. Use DMC2's `Likelihoods::Part` in
`src/sampler/likelihood.cpp` and `src/LikelihoodContext.h` as references.
Give the context explicit ownership, invalidation and restart behavior;
reconstruct native pointers when loading a checkpoint. Keep custom-model
fallbacks and verify identical likelihoods and draws. Measure the remaining
`calc_ll_oo` prologue before committing to this work.

### 6c. Move the particle step into C++ if it remains material

Port EMC2's `new_particle` from `R/sampling.R`, using DMC2's `particle.cpp`
and `proposals.cpp` for implementation patterns. Preserve EMC2's composition
of independent and auxiliary-centre local steps, weighting, retained particle,
marginalisation, RNG order and adaptation. DMC2 already has an auxiliary-centre
local step, but selects between local and independent kernels; use its
arithmetic without adopting that different transition schedule.

Keep group moves, shear, HMC, collapsed planning and the allocator in R.
Require transition/state parity on supported paths and posterior/efficiency
checks. Proceed only if cheap-likelihood fits still spend a material share
of runtime here: ~1.5 ms of overhead is under 10% when the corresponding
likelihood work costs 15 ms.

## Performance targets and scope

The immediate targets are to pay PSIS startup once per fork parent, remove
repeat fixed-conditioning fits/checks, avoid forks for cheap rebuild work,
and eliminate repeated parameter-only kernel calculations. Step 5 targets
chain-fork and worker-startup cost while preserving idle-core gains.

Re-estimate total fit time after every accepted step. Savings overlap: cache
reuse changes dispatch costs, warming changes startup cost, and persistent
processes change both. Establish targets from the residual profile rather
than treating cumulative wall-time forecasts as acceptance criteria. Repeat
the LNR benchmark, whose kernel is already competitive, to isolate the
benefit from setup and process changes.

Retain EMC2's runner architecture. Threaded execution is a separate future
project requiring per-kernel thread-safety checks, including static scratch
buffers, memoised solvers and R-callback custom models. DMC2's `thread_safe()`
gating is a reference if process overhead remains material after Step 5.

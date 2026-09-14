# Ordered execution plan for the remaining efficiency work

Date: 2026-09-14
Parent documents:

- `plans/architecture-efficiency-audit.md` — the original architecture plan;
- `plans/architecture-efficiency-review-2026-09-14.md` — what actually shipped
  and why the end-to-end benchmark did not improve.

This is an implementation plan, not a request to run another multi-hour fit.
Every stage below has a cheap gate. A later stage starts only if the earlier
stage is correct, numerically safe, and measurable on a short fixture.

## Operating rules

1. Keep the current worker-lifecycle fixes. Do not trade orphan prevention,
   generation checks, bounded reads or deterministic streams for a timing win.
2. Change one computational mechanism at a time. Keep a clean pre-change build
   and a clean post-change build; alternate run order and report median and
   spread over repeated short runs.
3. Use the native C++ path as the numerical reference. The obsolete R
   likelihoods are not an oracle. For storage, scheduling and routing changes,
   require bit-identical likelihood vectors and sampler state. For arithmetic
   hoists, state the permitted ulp/absolute error before implementation and
   preserve support, floors, branch decisions and reduction order.
4. Keep timing tests out of `testthat` assertions. Unit tests check behaviour;
   `WorkingTests/` scripts report time, CPU, memory and process counts. A noisy
   timing result is a reason to repeat, not to weaken a numerical test.
5. Do not start an hours-long hierarchical fit until a targeted test predicts a
   material gain. A final confirmation fit is one short/medium case per model,
   followed by one representative large case only when the predicted saving is
   at least a few percent.

## Recommended first four commits

This is the shortest safe path from the current tree to evidence:

1. **E0: baseline harness** — add the PCOUNTER fixture and the short direct-call
   benchmark; make the baseline matrix reproducible.
2. **E1: spare-core safety gate** — stop automatic multi-subject spare-core
   donation and add route/process-budget tests.
3. **E2: PCOUNTER probe** — add branch/reuse counters only; do not change
   arithmetic or scheduling.
4. **E3: one kernel experiment** — choose PCOUNTER only if E2 clears its gate;
   otherwise implement the RDM preparation experiment and leave PCOUNTER as a
   documented negative result.

Do not combine E0--E3 in one commit. Each commit must leave the package
installable and the targeted tests runnable on its own.

## Stage 0 — freeze a cheap baseline and expose the decisions

### Objective

Make each subsequent result attributable to one change, and ensure that the
benchmark is measuring useful likelihood work rather than an invalid floor or a
worker failure.

### Work

- Extend the existing audit fixtures (`tests/testthat/helper-audit-harness.R`)
  with a PCOUNTER fixture and a common direct-call timing helper. Keep the
  existing BAwL/RDM/DDM fixtures unchanged.
- Add a short script, preferably
  `WorkingTests/bench_efficiency_targets.R`, that runs only direct native
  likelihood calls and one particle step. It should accept `MODEL`, `TRIALS`,
  `CELLS`, `PARTICLES`, `BRANCH`, `REPS` and `PROFILE` through the environment.
- For each case: warm up three calls, measure at least seven calls, report
  median, MAD, CPU time, likelihood checksum, number of native calls, kernel
  rows/cells/seconds, and whether any result equals the invalid floor.
- Run with profiling off for speed numbers; run the same short cases with
  profiling on only to explain a result. Keep compiler flags and BLAS settings
  identical between builds.

### Required baseline matrix

| Workload | Trials | Particles | Design shape |
|---|---:|---:|---|
| direct mapper/kernel smoke | 64, 600 | 1, 5, 25 | one cell, 4 cells, 64 cells |
| mapper cliff | 600 | 5 | 255, 256, 257, 512 cells |
| short sampler step | 200 | 25, 100 | 1 and 8 subjects |
| high-core scheduler | 200 | 25, 100 | 2, 4, 8 subjects; 8/16/32 cores |

The high-core cases are deliberately short. They are intended to expose nested
forks, queue tails and process leaks, not to estimate posterior quality.

### Gate

The baseline must produce finite, informative likelihoods, identical checksums
on repeated runs, and no live workers after the call. If it does not, fix the
harness or the pre-existing defect before doing optimisation work.

## Stage 1 — remove the unsafe spare-core path before measuring new speed

### Objective

Prevent the newly introduced low-N/high-core route from masking later gains.

### Work

- Until a persistent combined executor exists, do not widen `ctx$r_cores` for
  multi-subject pools. In practice this means keeping `allow_spare = FALSE` for
  that route, or gating it behind an explicitly disabled development option.
- Preserve the existing behaviour when the user explicitly requests
  `r_cores > 1`; this stage only removes automatic donation of otherwise spare
  cores.
- Add a route/profile field saying whether a call used serial, persistent pool,
  or nested per-call splitting. A result that is faster because it silently
  changed route is not a valid comparison.

### Targeted tests

- Extend `tests/testthat/test-particle-core-budget.R` with the cases
  `N = 1, 2, 8`, total cores `1, 4, 8, 16, 32`, explicit `r_cores = 1` and
  `r_cores = 2`. Assert that the product of active subject workers, likelihood
  workers and BLAS threads never exceeds the requested budget.
- Add a short process-tree test around `WorkingTests/bench_worker_pool.R`:
  after 20 iterations, no retired worker, template or nested child remains.
- Exercise the no-FIFO/Windows fallback with the pool constructor mocked or
  unavailable. It must not silently claim spare-core use.

### Gate

No nested per-call executor is used for an ordinary multi-subject fit. The
short high-core matrix must be no slower than the current baseline by more than
measurement noise; otherwise inspect route selection before proceeding.

## Stage 2 — PCOUNTER feasibility probe (no optimisation yet)

PCOUNTER is a good low-risk laboratory because it has a self-contained native
closed-form implementation (`src/model_PCOUNTER.cpp`), an explicit raw batch
adapter, and extensive existing correctness tests (`tests/testthat/test-pcounter.R`).
It is not a promise that PCOUNTER is the largest package-level win; the probe
answers whether parameter-cell preparation can pay off in a model whose inner
formula is more expensive than the mapper.

### Workload design

Add PCOUNTER to the audit fixture with formulas that independently exercise:

- all parameters constant across trials;
- `nu` varying over 4 and 64 condition cells;
- `k` varying across cells (changes the integer threshold and recurrence range);
- continuous `nu` or `t0` (should produce no useful cell reuse);
- `sv = 0`, `gamma = 0`, `omega = 0` degenerate branches;
- non-zero `sv`, `gamma`, `omega`, including truncation/omission data;
- `k = 0`, moderate `k`, and the existing large-threshold case (`k = 1023`).

Use 64, 600 and 2,000 trials and 1, 5, 25 and 100 particles. Keep every call
informative and use the same proposal cloud for pre/post builds.

### Measurements to add

Instrument only the PCOUNTER raw path, behind the existing off-by-default
profile flag:

- calls, rows and elapsed time for `dpcounter_raw`, `ppcounter_raw` and
  `pcounter_logS_at_t`;
- branch counts for `sv_zero`, `gamma_zero`, `omega_zero`, fallback-tail use and
  invalid/support exits;
- distribution of `K = 2 + floor(k + 0.5)` and the number of Stirling/rising-factor
  terms evaluated;
- time spent preparing parameter-only values versus time spent in the
  RT-dependent sums.

The instrumentation must not allocate an R object or change a floating-point
operation when disabled. Add a declaration for PCOUNTER to the kernel-reuse
metadata only if the declaration describes an actually reusable subexpression;
do not label all six parameter columns reusable because `t0` and the time
dependent terms still vary per row.

### Decision gate

Proceed to an optimisation only if at least one realistic PCOUNTER case shows:

- a repeatable preparation/recomputation component of at least 10% of raw
  kernel time, and
- a cell pattern that occurs in ordinary designs (not only an artificial
  all-constant fixture).

If the result is below that threshold, record PCOUNTER as a negative result and
move the production experiment to RDM. Do not force a PCOUNTER special case for
the sake of completing the plan.

## Stage 3 — first production kernel: the smallest validated preparation

Choose exactly one of the following after Stage 2. The default choice is RDM;
choose PCOUNTER only if its Stage 2 counters meet the gate.

### Option A: PCOUNTER

Implement only parameter-cell preparation that has been measured. Likely
candidates are the exact `K`, `shape`, `rate`, logarithms and branch flags that
are recomputed inside `pcounter_log_h()`/`pcounter_log_pi_phi()`. Keep RT- and
`t0`-dependent quantities in the row loop. Treat Stirling-row storage as a
separate experiment: the current rolling implementation is O(K) memory, so do
not replace it with an unbounded K-squared table.

Pass a small immutable prepared context from the raw adapter to the evaluator,
with a bounded cell cache. Fall back to the current scalar path when cells are
not repeated, K is large, memory is insufficient, or a branch is unsupported.

### Option B: RDM (default production candidate)

Prepare the exact shared scale/geometry terms (`1/s`, `(B + A/2)/s`, `v/s`,
`A/(2s)`) once for the common refinement of the columns that feed them. Feed
the same prepared values to the winner-density and survivor paths. Preserve the
existing operand order deliberately; do not rely on algebraic equivalence for
bit identity.

### Numerical tests before timing

- Add a test-specific reference route that runs the current scalar formula, or
  compare two clean package builds from before and after the change.
- Compare per-trial and total log likelihoods on ordinary, boundary and
  adversarial grids. Include zero/near-zero scale, invalid parameters,
  `rt <= t0`, `Inf`, omissions and truncation.
- Compare all race responses and all particle rows, not only the summed total.
- Check that changing particle order, chunk size, cell order and worker count
  does not change results.
- Run the existing model-specific tests and the audit differential/numerical
  suites.

### Timing gate

Require a repeatable improvement of at least 10% in the targeted direct kernel
case, no more than 2% regression in the narrow/continuous-covariate cases, and
no increase in peak memory beyond the bounded cache budget. If it fails, revert
the specialisation and retain the instrumentation.

## Stage 4 — move the validated kernel pattern across models

Do not implement BAwL, DDM and PCOUNTER/RDM together. Repeat the same loop for
each model:

1. declare the exact parameter dependency set;
2. measure cell reuse and preparation time;
3. implement one immutable prepared context with a scalar fallback;
4. run numerical and branch tests;
5. run short direct-call timings;
6. promote only if the timing gate passes.

Recommended order:

1. RDM — few reusable scale terms and a common raw/adaptor structure;
2. PCOUNTER — only if Stage 2 demonstrates a material repeated preparation;
3. BAwL — normaliser and kill/guess branches, with exact denominator/floor
   handling;
4. DDM — scale/log/sine terms, preserving series selection, reflection and
   summation order.

For each model, include narrow, 4/64-cell, 257-cell, continuous-covariate,
omission/mixed and truncation fixtures. The 257-cell case verifies that the
kernel change does not accidentally restore the mapper cliff.

## Stage 5 — finish the mapper only where route counters justify it

### Scalar representation

Add an explicit design-plan classification for a genuinely scalar reader. Do
not infer it from a repeated first value. Auto-added constant parameters may be
represented as constants, while a general `~ 1` design that a transform or
self-reference can read at arbitrary rows must remain materialised.

Tests:

- extend `test-param-table-prologue.R` for automatic constants, explicit
  intercepts, self-intercepts, split transforms, bounds and mixed cell/scalar
  columns;
- compare mapped natural parameters and bound masks to the independent wrapper
  bit-for-bit;
- run assertion builds with `EMC2_PT_CONST_CHECK`;
- report the fraction of columns and rows taking scalar, cell and row routes.

### Cell admission

Use short measurements to choose a minimum reuse ratio or a cost model. Test
255/256/257/512 cells, one-cell designs, continuous covariates and small trial
counts. Never make the cutoff a magic number without a benchmark showing that
the gather/scratch overhead is paid back.

Gate: no route may change mapped values, bounds or rejection counts; a new
cutoff must improve the direct prologue on the cases it admits and not regress
narrow/continuous cases.

## Stage 6 — immutable native context and mixed-data path

This stage is a prerequisite for a safe persistent combined executor, not a
stand-alone promise of a large speedup.

### Work

- Separate data/design-dependent immutable state from particle values.
- Reuse `PtMapper` plans, transform specifications, data masks and mixed-path
  partition views within a worker or likelihood context.
- Remove remaining per-particle vector copies in the mixed path only when the
  owner lifetime is explicit.
- Invalidate/rebuild on data, factor levels, designs, constants, transforms,
  bounds, trend and guess-window changes; rebuild after deserialization and
  never serialize an external pointer as if it were portable.

### Tests

Use `tests/testthat/test-native-context.R`, `test-ll-data-cache.R`,
`test-ss-raw-path.R`, `test-nan-censoring.R` and the differential harness. Add
repeat-call, reordered-particle, mixed-omission, truncation, trend and
deserialization cases. Compare context-enabled and context-disabled routes
bit-for-bit.

Timing is limited to 64/600/2,000 trial direct calls and a 20-iteration particle
step. Do not run a full fit until context construction is below the saved work.

## Stage 7 — replace nested spare cores with one bounded executor

Only start this stage after Stage 1 and Stage 6 are complete.

### Design constraint

There must be one explicit budget for subject tasks, particle tasks, native/BLAS
threads and memory. A subject worker must not start a fresh process pool for
each likelihood call. Deterministic likelihood tasks carry particle indices,
draw no random numbers and write ordered result slots; subject-local RNG streams
remain in the owning sampler process.

### Incremental implementation

1. Keep static LPT as the control route.
2. Add a bounded whole-subject queue using the existing framed protocol, but
   only when measured tail imbalance exceeds queue overhead.
3. Add particle sub-tasks only for the low-subject/high-core case, using the
   same persistent executor and a fixed memory cap.
4. Preserve serial, fork-pool and Windows/macOS fallback routes. If a platform
   cannot provide the persistent executor, use the safe static route rather
   than per-call spawning.

### Cheap tests and benchmarks

- `test-particle-core-budget.R`: budget, explicit-width and reproducibility
  assertions;
- `test-audit-differential.R`: worker count, partition, recycle and fallback
  equivalence;
- `test-audit-faults.R` and the liveness tests: kill workers/templates at send,
  receive and recycle boundaries;
- `bench_worker_pool.R`: 20 and 50 iterations, 2/4/8 subjects, 8/16/32 cores,
  25/100 particles, with process-tree snapshots;
- one 200-iteration pool run only after the short matrix is clean.

Compare static LPT, the new queue and the serial fallback. Require no orphaned
processes, no budget oversubscription, identical final sampler state, and a
measurable reduction in worker tail/queue time before enabling the queue by
default.

## Stage 8 — blocked proposal density and proposal work

Benchmark this separately from the ordinary one-component path.

- Add short blocked/marginal fixtures with two and three proposal components,
  `P = 8, 32, 64`, 20--50 iterations and 25/100 particles.
- Verify C15's full prior factor cache is hit only for the same group draw and
  never replaces the correlated prior with independent component priors.
- Profile matrix slicing, centring, proposal-density products and mixture
  allocation before fusing anything into C++.
- If a native proposal-density batch is justified, preserve R-side proposal
  generation and RNG order first; move only deterministic density evaluation.

Gate on ESS/sec and numerical equivalence, not raw proposal-call time alone.
Skip this stage if blocked proposals are not a material fraction of real fits.

## Stage 9 — history and worker retention, only if block profiling pays

Use the C16 block profile to decide whether this is worth implementation.

- If history joins dominate, implement owned growable capacity with separate
  allocated and committed lengths. Keep public arrays, names, checkpoints,
  resume and copy isolation unchanged.
- If startup dominates, retain a clean worker template/context within a fit and
  explicitly version stage, proposal and transform changes. Keep hard recycle
  limits based on private memory and a safe fallback.
- Test `test-history-join.R`, `test-sample-store-inplace.R`, checkpoint/resume,
  multiple stage blocks and deserialization with short arrays.

Use synthetic wide/long stores to measure joins; do not use a full posterior fit
as the first test. A full fit is justified only if block profiling predicts a
material wall-time or memory benefit.

## Stage 10 — particle count adaptation experiment

Do not change the rule from C17 based on particle-weight ESS alone.

1. Run short fixed-particle chains at factors 25, 50 and 100 on small BAwL,
   RDM and DDM fixtures, recording wall time, rejection rates and particle
   weight ESS.
2. Run enough iterations for chain ESS/R-hat to be meaningful on a reduced
   parameter set; report posterior-chain ESS/sec for population parameters,
   covariance parameters and at least one condition contrast.
3. Only then test a bounded adaptation rule in a short multi-chain run, with
   frozen-particle controls and posterior equivalence/coverage checks.

This stage is statistical validation, not a free compute optimisation. Stop if
the saved particle work does not translate to equal-or-better posterior ESS/sec.

## Stage 11 — final promotion and one representative fit

After all preceding gates pass:

- run the full test suite and the assertion build;
- repeat the direct-call matrix with profiling disabled and enabled;
- run one short/medium hierarchical fit for each of BAwL, RDM and DDM (for
  example 8--32 subjects, 50--100 iterations per stage, 25--100 particles);
- verify posterior checksums, sampler state, worker counts, peak process-tree
  memory and no orphaned workers;
- run the large `N = 140, P = 64` case only if the short results predict at
  least a 3% end-to-end gain;
- run a wide `P >= 257` case to demonstrate the mapper/kernel interaction;
- document negative results and disabled routes, not just successful changes.

Do not start a wholesale PMwG rewrite or a GPU/native sampler migration unless
this staged work demonstrates a remaining bottleneck that cannot be addressed
without changing the public inference semantics.

## Definition of done

The remaining work is complete when each enabled optimisation has:

- a measured affected workload and an unchanged-workload control;
- unit, differential, edge-case and cross-platform fallback tests;
- explicit numerical equivalence/error bounds;
- a bounded memory/process budget and no orphan regression;
- a repeatable direct-call benefit that survives into a short particle step; and
- a documented reason for enabling it in the default route.

An optimisation that does not pass those gates remains an instrumented,
opt-in experiment or is removed. That is preferable to another broad commit
whose only observable effect is extra orchestration overhead.
 
## Execution record — 2026-09-14

The staged implementation and promotion gates produced the following final
routes and decisions:

- Stage 0 is enabled. The direct-call harness covers BAwL, RDM, DDM and
  PCOUNTER, including default/varying branches, 64/600 trials, 1/4/64 cells
  and 1/5 particles. The split profiling-disabled matrix completed with
  finite likelihoods, reproducible checksums and no invalid-floor results.
- Stage 1 is enabled. Automatic spare-core donation is disabled for the
  default single-core route; explicit `r_cores > 1` remains available.
  Route and budget diagnostics distinguish serial, persistent-pool and nested
  fallback execution. Core-budget, process-liveness, fault-injection and
  worker smoke checks passed.
- Stage 2 is instrumented but negative for a PCOUNTER kernel rewrite.
  Realistic probes measured preparation fractions from approximately
  0.0035% to 0.11%, below the 10% gate. Branch, support, threshold,
  Stirling/rising-term and preparation/RT counters remain available behind
  the off-by-default profile flag.
- The RDM prepared-geometry experiment was reverted after bit-identical
  numerical checks passed but direct timings ranged from a 2.8% regression
  to no improvement. The scalar RDM route is retained.
- Stage 5 is enabled. The mapper now uses structural scalar-reader
  classification, bounded cell admission, route counters and explicit
  scalar-marker retirement on full-width writers. Normal and
  `EMC2_PT_CONST_CHECK` builds passed the prologue and equivalence checks.
- Stage 6 retains explicit per-call native ownership. Invalidation-key,
  mixed-data, truncation/omission, trend, deserialization and differential
  checks passed; cross-call context caching was not promoted because its
  fixed construction cost was small and a mutable cache would add risk.
- Stage 7 retains static LPT and the persistent worker pool. Two- and
  eight-subject worker runs were balanced and orphan-free, but measured
  queue delay did not provide a tail benefit, so a combined queue and
  particle-subtask executor were not enabled.
- Stage 8 was not promoted. Short blocked proposal-density runs failed
  safely before producing a usable workload (insufficient iterations and a
  downstream indexing error); the existing proposal-cache and density
  correctness tests pass, so no speculative native batch was added.
- Stage 9 was not promoted. Synthetic history-join timings were noisy and
  did not establish a consistent dominant cost; the existing bounded worker
  retention and sample-store paths remain unchanged.
- Stage 10 was not promoted. The representative RDM hierarchical regression
  measured nearly unchanged ESS/sec across particle factors 25, 50 and 100,
  so no adaptive default was introduced.
- Final validation included normal and assertion installs, a complete
  installed-package `testthat` suite with no failures (35 development or
  availability skips and six warnings, including two worker-pool
  fault-recovery warnings), focused testthat suites, the split direct matrix,
  and a representative RDM hierarchical fit. A full `R CMD check` was
  attempted but timed out during the executable-file/FIFO check; it is not
  recorded as a passing `R CMD check` result.

Only routes with a measured safety or correctness case are enabled by
default. The negative experiments remain documented here rather than being
left as dormant production special cases.

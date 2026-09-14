# Why the architecture-efficiency implementation did not speed up sampling

Date: 2026-09-14
Compared against the audit plan in `plans/architecture-efficiency-audit.md` and
the measurements in `plans/architecture-efficiency-benchmark-2026-09-13.md`.

The ordered follow-on implementation and test sequence is in
`plans/architecture-efficiency-execution-plan-2026-09-14.md`.

## Executive diagnosis

The result is consistent with the code and the benchmark. The implementation
removed one real mapper cliff, but most of the proposed sampling wins were
either deliberately deferred, implemented only as instrumentation/safety work,
or activated only outside the benchmark's ordinary workload. At the same time,
the worker protocol now has a small per-iteration cost. Consequently the
end-to-end result is near the noise floor and can be slower on short or
under-parallelised fits.

The most important conclusions are:

1. The large measured mapper win is real, but is a small fraction of an
   ordinary particle step and is absent from the benchmark's end-to-end cases:
   the `P = 250` stress cases were not run, while the normal fixtures have
   `P <= 64`.
2. C9 did not implement model-specific kernel hoisting. It records how many
   cell evaluations a future specialised kernel could use; the BAwL/RDM/DDM
   kernels do not consume that plan yet.
3. C13 did not implement dynamic whole-subject scheduling. Static LPT remains
   the scheduler, and the new spare-core path can create short-lived nested
   `mclapply` processes inside persistent subject workers. That is precisely the
   kind of process churn that can make low-subject/high-core fits slower.
4. C14, C15, C16 and C17 do not change the dominant path for the benchmark:
   history joining is still quadratic, the prior cache is gated to blocked
   proposals, worker/context retention was not implemented, and the particle
   adaptation rule was not changed.
5. C4-C6/C9 and the lifecycle fixes are valuable correctness and observability
   work, not speedups. The benchmark found and fixed an initial 12% pool
   regression; approximately 0.5 ms per iteration of bounded protocol overhead
   remains.

This is therefore not evidence that PMwG or the native likelihood architecture
needs to be replaced wholesale. It is evidence that the implementation was
mostly a prerequisite tranche, while the expected gains were attributed to it
as if the later optimisation stages had shipped.

## What the benchmark actually proves

The benchmark is internally coherent and its posterior-equivalence checks are
important. Its results should be read as follows:

| Measurement | Interpretation |
|---|---|
| BAwL/RDM/DDM, 64 design levels: about 7--8% faster per native call | The row-constant mapper path works. This is a prologue saving, not a kernel saving. |
| 256 cells: about 18% faster; 257 cells: 113.6 -> 45.0 ms | C7's dynamic scratch removes the fixed 256-cell cliff. This is the clearest delivered speedup. |
| Narrow calls and one-subject particle steps: unchanged | The native race kernel still dominates; the changed mapper/prologue is too small to move the total. |
| Pool after fixes: about level, with roughly 0.5 ms/iteration remaining | Protocol/lifecycle code is now safer, but it is not free. |
| Full fits: -2.4% to +3.4% | Most differences are within the stated fit noise, except that the small `N = 8, P = 64` case exposes overhead. |

The end-to-end ceiling was low even under ideal integration. The benchmark
reports roughly 1.8 ms of prologue in a 37--40 ms narrow call. Optimising that
slice cannot produce a large fit-level gain unless the design has many cells or
the native kernel itself is changed. The only end-to-end cases where the
mapper's wide-design behaviour could dominate (`P = 250`) were explicitly
omitted.

## Plan-to-code audit

### P0: worker lifecycle, framing and telemetry (C1--C6, later fixes)

This work largely landed and fixed real reliability defects. It was not an
expected source of positive throughput. C4/C5/C9 initially added enough
polling, allocation and connection overhead to slow the pool by 12% at 20
iterations and 5.5% at 200 iterations. The follow-up fixes (`c9ad349e`,
`ad491d90`, `4f6f4ae1`) brought the measured pool back to approximately the
baseline, but they cannot make it faster than the old protocol. Treat this
tranche as a prerequisite and a correctness win, not as a performance result.

The current per-iteration path still performs framed completion reads, deadline
and generation checks, and R-level request/reply bookkeeping. That cost is
visible when particle work is short. It is a reasonable trade for eliminating
orphaned workers and false liveness, but it must be included in any speed
model.

### P1 mapper work (C7, C8, C10)

C7 is the one complete, demonstrated optimisation. `ParamTable` now admits
design-cell scratch according to a byte budget instead of a fixed 256-cell
constant (`src/ParamTable.h:464-493`). The benchmark confirms the cliff is gone.

C8's `JointCells` abstraction exists, but it is not a consuming kernel API.
The runtime uses it in `emc_reuse_cells()` to count potential reuse for
instrumentation (`src/particle_ll.cpp:67-86`, `:848-850`, `:1251-1253`) and in a
diagnostic interface (`src/param_table_interface.cpp:380-400`). No BAwL, RDM or
DDM hot loop evaluates a reusable term once per joint cell. Thus the apparent
"joint cell partition" feature should not be counted as a kernel speedup.

The deferred scalar-fill part of C8 is also mostly inactive in the standard
user workflow. `put_scalar()` defers only when `out_to_entry[b] < 0`
(`src/ParamTable.h:1215-1233`). `make_design()` adds a `parameter ~ 1` design
for every unspecified model parameter (`R/design.R:360-368`), so ordinary
constant parameters still have an output design entry. The row-constant path
can still help, which explains the 64-level prologue results, but the intended
compact scalar representation is not broadly engaged. The safe fix is to add a
design-plan flag for genuinely scalar readers (or represent auto-constant
parameters as constants), not to remove the guard globally.

C10's explicit update mask is semantically correct, but the C++ side verifies
every column marked unchanged by scanning all particles
(`src/particle_ll.cpp:570-608`). This is an O(particles x parameters) check on
each native call. It is a useful correctness gate, but it is recurring work and
has no cached proof for the common proposal structure; it can cancel a small
mapper saving for short batches.

The cell admission code intentionally has no minimum reuse ratio
(`src/ParamTable.h:464-480`). A design with only one fewer cell than trials can
therefore pay cell bookkeeping and gathers without enough reuse to win. This is
not the cause of the reported broad slowdown, but it is a sensible follow-up
once real per-design route counters exist.

### P1 kernel preparation (C9)

The plan called for measured, model-specific preparation of repeated BAwL,
RDM and DDM terms. The implementation stops at measurement. `kernel_stats.cpp`
declares the reusable columns, and `particle_ll.cpp` records `rows`, hypothetical
`cells` and elapsed time; there is no cached BAwL normaliser, RDM geometry or DDM
scale/log/sine table passed to the kernels. This is the largest unshipped item
relative to the plan and the main reason the likelihood-dominated portion of a
narrow fit did not improve.

The next implementation should only proceed if the counters show a useful
reuse ratio on actual BAwL/RDM/DDM designs. It must preserve operand order,
support/floor decisions and the existing finite/mixed path split; algebraically
equivalent hoists are not automatically bit-identical.

### P1 spare cores and scheduling (C13)

The source comments explicitly retain static LPT and defer dynamic whole-subject
assignment pending measurement (`R/sampling.R:328-345`). The pool re-partitions
between iterations, but a worker cannot take a subject from another worker's
partition after it finishes (`R/chain_pool.R:1831-1940`,
`:2181-2200`). The dynamic queue in `chain_pool.R:1623-1718` is a queue for a
single-subject likelihood pool, not the multi-subject sampler's tail scheduler.

More seriously, the spare-core budget is wired into `ctx$r_cores`
(`R/sampling.R:690-700`) and each subject worker calls `new_particle()` with
that width (`R/chain_pool.R:1513-1525`). For a typical `N = 8`, `n_cores = 32`,
`r_cores = 1` fit, the budget is eight subject workers and four likelihood
cores per worker. There is no persistent inner executor for those calls;
`calc_ll_pooled()` falls back to `calc_ll_manager()`, which repeatedly launches
`auto_mclapply()` for the proposal rows (`R/sampling.R:1867-1876`,
`:1927-1974`). On Unix this is nested fork/process creation inside persistent
workers; if that width is enabled on Windows, the same abstraction creates a
new PSOCK cluster per call. This is a credible direct cause of the `N = 8,
P = 64` slowdown and is not exercised by
the published core grid (`N = 32`, at most 16 cores).

The non-pool fallback omits `allow_spare = TRUE` (`R/sampling.R:931-947`), so
Windows and other systems without FIFOs do not receive the spare-capacity
allocation at all. That is safe but means the claimed package-wide utilisation
improvement is platform-dependent.

Do not widen `r_cores` this way until there is one bounded, persistent executor
for the combined subject/particle work, or until the spare path is disabled by
default. A nested per-call fork is not a speed optimisation.

### P1 immutable native context and mixed path (C11)

C11 fixed ownership and reduced some mixed-path vector copies, but the central
cache is explicitly not reused across calls. `PtMapper` is constructed per
`calc_ll_oo` call (`src/particle_ll.cpp:516-570`, `:785-787`), and its comments
state that cross-call reuse is intentionally absent. This is a prerequisite
for native threading, not a delivered speedup.

### P2 history and block setup (C14, C16)

C14 replaces `abind` with a contiguous copy for the compatible array case, but
its own comment states that total copying remains quadratic
(`R/objects.R:540-590`). The benchmark's four large stage blocks do not make
this a dominant cost. The planned appendable capacity/committed-length design
was deferred.

C16 records block-level timings only (`R/fitting.R:120-170`,
`R/profile_schema.R:508-540`). `run_stage()` still starts a pool and registers
`on.exit(.emc_wpool_stop(...))` for each block (`R/sampling.R:684-715`); workers,
templates and native context are not retained across blocks. Therefore the
planned startup amortisation is not present. Repeated pool startup can matter
for short stages, but this commit cannot improve it.

### P2 proposal density (C15)

C15 caches the full group factor only when there is more than one unique
proposal component (`R/sampling.R:1221-1248`, `:1424-1436`). The standard
benchmark uses the default one-block/one-component path, where the existing
proposal factor is already sufficient and the new cache is not used. The cache
is a valid targeted optimisation for blocked fits, but blocked/marginal fits
were not benchmarked, so it cannot explain the headline result.

### P2 particle adaptation (C17)

C17 records particle-weight ESS but leaves the adaptation rule unchanged
(`R/sampling.R:1607-1650`). A fixed particle-factor benchmark cannot show a
speedup from this commit. Particle-weight ESS is also not posterior-chain ESS;
changing the rule requires converged, multi-chain ESS-per-second experiments.

## Why the average can be slower

The observed small regressions have three plausible, code-supported sources:

- residual framed-pool and R orchestration overhead on every iteration;
- repeated pool/template startup at stage blocks because C16 retention was not
  implemented; and
- the new spare-core allocation creating nested, short-lived inner process
  pools for small-N/high-core fits.

The first two are measurable fixed costs. The third is a path-selection problem
that should be treated as a correctness-of-budget issue even if it is not the
dominant path for large-N fits. The benchmark did not include the combinations
that select that path, so it cannot rule it in or out from timing alone.

## Recommended corrective sequence

1. **Stop the nested spare path immediately.** Keep `allow_spare = FALSE` for
   multi-subject pools unless a persistent combined executor is available.
   Add a benchmark matrix with `N = 2, 4, 8`, cores 8--32, and both pool and
   PSOCK/fallback backends. Record process counts, per-worker CPU, queue delay,
   and likelihood route; require no worker growth beyond the configured budget.
2. **Measure before adding more scheduler machinery.** Use the corrected
   `worker_cpu_max` and queue telemetry to compare static LPT against a bounded
   whole-subject queue. Only switch when tail idle time exceeds queue/IPC cost.
   Preserve subject-local RNG streams and ordered result slots.
3. **Finish C9 or remove its speed claim.** Implement one model at a time (RDM
   is the cleanest first candidate), using the measured joint partition and an
   immutable per-call/per-worker preparation object. Benchmark narrow, 64-cell,
   mixed/omission and continuous-covariate cases with exact numerical gates.
4. **Make scalar compact storage match real designs.** Distinguish a true
   intercept-only/scalar design from a general output design, then defer or
   materialise safely. Add route counters so a benchmark reports how many
   columns/cells actually took each path.
5. **Only then address retention/history.** If block profiling shows material
   cost, retain worker templates/context within a fit and replace quadratic
   joins with owned capacity plus committed length. Otherwise leave this as a
   large-fit memory project rather than a sampler-speed promise.
6. **Run a representative benchmark before any wholesale rewrite.** Include
   BAwL, RDM and DDM; `P = 8, 64, 257, 512, 250`; omitted/truncated data;
   blocked proposals; `n_blocks > 1`; small-N/high-core fits; Windows/macOS
   fallback; repeated quiet-machine runs; and converged ESS/sec, not only fixed
   wall time. Keep the existing posterior equivalence checks.

## Bottom line

The implementation did not fail mysteriously: it delivered a robust worker
protocol and a real wide-mapper fix, while the main sampling optimisations in
the plan remain conditional or unimplemented. The benchmark's near-zero
end-to-end result is therefore the expected outcome. The highest-confidence
action is to remove or redesign the nested spare-core route, then finish one
measured native kernel specialisation and validate it on the workloads where
its reuse is actually present.

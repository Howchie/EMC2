# Postmortem: why the architecture work did not speed up sampling

Date: 2026-09-14

This document records the result after implementing the staged work in
`plans/architecture-efficiency-execution-plan-2026-09-14.md`. The result is
not a mysterious failure of the sampler. The work was mostly a safety,
measurement and routing tranche; the few changes that reached production did
not remove enough work from the user-facing sampling path.

## Executive result

The package has a real direct-call improvement for wide design matrices, but
there is no accepted optimisation of the dominant ordinary sampling path.
The benchmark therefore correctly reports an end-to-end result near the noise
floor (and occasional regressions on short or lightly parallel fits).

The previous benchmark measured approximately 1.8 ms of mapper/prologue work
inside a 37--40 ms narrow native likelihood call. Even eliminating that whole
slice could improve a narrow call by only about 4--5%, before accounting for
the sampler, proposal work, R/IPC and worker overhead. The native race kernel
and sampler orchestration were not changed. This puts a hard ceiling on the
benefit available from the shipped mapper changes.

## What actually shipped

| Area | Final state | Why it did not improve ordinary sampling |
|---|---|---|
| Worker lifecycle/framing | Enabled and safer | Liveness, bounded reads and generation checks add a small fixed protocol cost; this work was not a throughput optimisation. |
| Spare-core handling | Unsafe automatic donation disabled; explicit width retained | No additional parallel capacity is being used by default. |
| Mapper/cell scratch | Enabled | The 256-cell cliff is removed, but normal fixtures do not reach it; cell values are still scattered to all trial rows. |
| Scalar deferral | Enabled in structurally safe cases | The deferred column is widened by `ensure_full()` before the ordinary row-oriented kernels consume it, so much of the write is shifted rather than removed. |
| Model kernel reuse | Not implemented | The reuse metadata and counters are observational; BAwL/RDM/DDM kernels do not consume prepared terms. |
| PCOUNTER rewrite | Not promoted | The negative gate was based on an invalid preparation timer (see below). |
| Native context reuse | Not promoted | Contexts and mapper plans are still rebuilt per likelihood call. |
| Dynamic scheduler/combined executor | Not promoted | Static LPT remains; measured queue tail did not pay for extra IPC. |
| Blocked proposal fusion | Not promoted | The short blocked workload did not reach a usable gate. |
| History capacity/startup retention | Not promoted | Synthetic timings did not identify a dominant cost. |
| Particle adaptation | Not changed | Recording particle-weight ESS is not a speedup and did not establish posterior ESS/sec benefit. |

Thus “the whole plan was implemented” should not be read as “the whole plan
was promoted into hot production paths”. The promotion gates rejected every
large computational change, as they were designed to do when evidence was
insufficient.

## Important measurement defect: PCOUNTER preparation was not measured

The PCOUNTER negative result is not currently a valid feasibility result.
In `src/model_PCOUNTER.cpp`, all three raw adapters set `prepared` immediately
after extracting column pointers:

```text
dpcounter_raw       : 315--327
ppcounter_raw       : 359--371
pcounter_logS_at_t  : 410--423
```

The row loop, including every call to `pcounter_log_eval()`, starts after that
timestamp. Parameter-dependent work such as branch setup, `K`, `shape`,
`rate`, logarithms, rising terms and Stirling rows is therefore charged to
`rt_sum_seconds`, not `preparation_seconds`. The reported 0.0035--0.11%
“preparation fraction” is only adapter pointer/entry setup. It cannot justify
rejecting a PCOUNTER prepared-context experiment.

This does not imply that PCOUNTER will win. It means the experiment must be
corrected before that conclusion is made. The corrected split must time or
count the parameter-only work inside the evaluator, while keeping RT- and
`t0`-dependent work in the row loop. Timing instrumentation must remain
off-by-default and must not be used as a speed number.

## Two implementation effects that can hide small gains

### 1. Diagnostic route counters are still on every mapper call

`ParamTable.h` calls `emc::mapper_count_map()` unconditionally on scalar, cell
and row routes. The definitions in `src/emc_scratch.cpp` are out-of-line and
perform an atomic flag load even when counters are disabled. Without LTO this
is a function call plus an atomic operation in the mapper hot path. That can
erase a small mapper saving, particularly for narrow designs.

The counter API should be changed to a genuinely zero-cost disabled route
(for example an inline branch through a cached non-atomic flag, or a compile-
time no-op for production builds). Then run an A/A direct matrix with all
telemetry disabled and enabled separately. Do not claim a mapper result until
the disabled path is no slower than the pre-instrumentation control.

### 2. Cell and scalar representations do not reach the kernel

The cell route computes one value per cell but then executes
`outc[r] = post_val[ex[r]]` for every trial row. The ordinary race path then
materialises full-width parameter columns. Scalar deferral similarly calls
`ensure_full()` from `materialize_into()` before the kernel receives the
matrix. The current API is therefore still a trial-row matrix API; it cannot
capture the larger saving suggested by cell-level reuse.

To obtain a material benefit, a model kernel must consume an immutable
cell/scalar representation directly, with a bounded fallback to the existing
row matrix. Optimising the mapper while retaining a mandatory full-width
materialisation is only a prologue optimisation.

## Evidence classification

Strong evidence:

- posterior and sampler-state equivalence checks passed for the reported
  comparisons;
- the wide-design mapper cliff disappeared;
- worker fault/liveness and budget tests passed;
- ordinary narrow particle steps and full fits remained approximately level.

Insufficient or misleading evidence:

- PCOUNTER preparation fractions, because the timer boundary is misplaced;
- any claim that `JointCells` is a kernel optimisation, because no production
  kernel consumes it;
- any claim of dynamic tail balancing, because static LPT remains the default;
- any claim of context/startup amortisation, because cross-call ownership was
  not promoted.

The correct conclusion is not “all optimisation ideas failed”. It is “the
implemented changes did not alter the dominant work, and one negative probe
needs a corrected measurement.”

## Correct next sequence

Do not start another broad optimisation pass or a wholesale PMwG rewrite yet.
Use this short sequence:

1. **Restore measurement integrity.** Fix the PCOUNTER timing split and make
   disabled mapper counters zero-cost. Rebuild and run the direct A/A matrix
   (64/600 trials, 1/4/64 cells, 1/5/25 particles) plus the existing numerical
   and assertion tests.
2. **Profile one complete particle iteration.** Separate proposal generation,
   proposal-density evaluation, native likelihood, group Gibbs, R/IPC wait,
   history/update and pool startup. Report medians and A/A spread; do not infer
   a target from a direct kernel benchmark alone.
3. **Choose one component that is at least 20% of the iteration.** For a
   likelihood target, use RDM first and implement one immutable prepared
   context consumed by the kernel. For a sampler target, move only
   deterministic proposal-density work to native code while preserving R-side
   RNG order. Keep the scalar/row fallback.
4. **Gate cheaply.** Require bit-identical likelihoods and sampler state,
   at least 10% improvement in the targeted direct/iteration fixture, no more
   than 2% regression in continuous/narrow controls, and no process or memory
   regression. Revert immediately if a gate fails.
5. **Run short end-to-end fits only after the iteration gate.** Use BAwL,
   RDM and DDM, 8 subjects, 64/600 trials, 25/100 particles and 20--50
   iterations. Require at least 5% wall-time and no relevant ESS/sec loss.
   A long fit is justified only if these short runs predict a gain.

If no single component is at least 20% of a complete iteration, the remaining
work is dominated by fixed PMwG/R orchestration. At that point, further mapper
micro-optimisation is not a good use of time; the choices are a bounded native
proposal/Gibbs implementation or a carefully benchmarked sampler redesign,
not more route counters.

## Tests required for the next experiment

- `test-pcounter.R` plus a new counter test proving that parameter-only work
  is separated from RT-dependent work;
- `test-param-table-prologue.R`, `test-audit-mapper-equivalence.R` and
  `EMC2_PT_CONST_CHECK` for scalar/cell/row equivalence;
- direct A/A timing with profiling disabled and enabled;
- `test-particle-core-budget.R`, worker fault/liveness tests and process-tree
  cleanup checks;
- fixed particle-step checksums across serial, explicit multi-core and
  fallback routes;
- only after the above pass, short BAwL/RDM/DDM fits with wall time and ESS/sec.

No hours-long fit is needed to establish any of these gates.

## Validation run — 2026-09-14

The current tree was installed into a clean temporary library and the planned
cheap checks were run.

- Focused tests passed: PCOUNTER instrumentation (23), parameter-table
  prologue (48), mapper equivalence (625), differential routes (114, with
  three CRAN skips), particle-core budget (640), audit faults (164, with eight
  CRAN skips), and worker pool (205, with four CRAN skips and two expected
  fault-recovery warnings).
- The full source test directory exited successfully with no test failures.
  It reported the package's normal availability/development skips and five
  warnings, including the intentionally injected worker-recycle warning.
- The 60-case direct matrix (BAwL, RDM, DDM, PCOUNTER and LBA; 64/600 trials;
  1/4/64 cells; 1/5 particles) produced finite, reproducible checksums with
  no invalid-floor results.
- The profiled BAwL/RDM/DDM matrix likewise stayed finite and reproducible.
- A clean `EMC2_PT_CONST_CHECK` assertion build installed successfully; its
  48 prologue and 625 mapper-equivalence checks passed.
- Short smoke fits completed without numerical failure: RDM 10.80 s, BAwL
  11.01 s, DDM 10.44 s and PCOUNTER 14.39 s at eight subjects, 64 trials,
  particle factor 5 and 10/10/10/30 stage iterations.
- A higher-resolution PCOUNTER truncation probe (2,000 trials, 25 particles)
  reported 0.000032 s “preparation” versus 1.367 s row time, while enabling
  the profile increased elapsed time from 0.144 s to 0.157 s. This confirms
  that the current split is not a valid preparation measurement and that the
  profile path is unsuitable for speed comparisons.

The first direct `test_file()` invocation failed only because testthat did not
attach the package's lazy data (`forstmann`) in that invocation. Rerunning with
the installed package attached passed all affected tests; this was a test
runner setup issue, not a package failure.

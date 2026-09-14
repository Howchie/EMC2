# All-model sampling-efficiency benchmark protocol

Date: 2026-09-14

## Purpose

Measure whether the architecture work improves the speed users actually see
when sampling. PCOUNTER is a feasibility probe only. It is included in every
comparison, but no PCOUNTER-only production route is accepted: a measured
kernel pattern must either be implemented through the same prepared-context
contract for every applicable model or be reverted.

## Comparison discipline

Run a clean pre-change build and a clean post-change build with the same
compiler, BLAS, CPU affinity, R version, seed, data and environment. Alternate
pre/post order for each case. First run five A/A repeats of the same build to
estimate the machine and installation noise; then run at least seven
interleaved pre/post repeats for cases that survive the short gate.

Use one process per benchmark run and keep the machine otherwise idle. Pin
BLAS/native threads explicitly (`OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`) before R starts. Record the
library path, commit, compiler flags, BLAS vendor, detected cores and git
status in the output file.

Every matched run must report:

- wall time and CPU time for preparation, each sampler stage, cleanup and the
  complete fit;
- per-iteration particle time, worker CPU/elapsed time, queue/response wait,
  startup, route (`serial`, `persistent_pool`, `nested_per_call`) and peak PSS
  for the entire process tree;
- worst relevant bulk and tail ESS per total wall second, with sample-stage
  ESS/sec beside it;
- likelihood and posterior checksums, finiteness/floor checks, native call
  counts and worker/process cleanup status.

Do not compare iteration counts alone. A change that makes iterations cheaper
but changes the route, sampler state, or ESS is not a speed improvement.

## Models and matched fixtures

Run the same matrix for `RDM`, `LBA`, `BAwL`, `DDM` and `PCOUNTER`. The current
hierarchical harness accepts `EMC_MODEL=PCOUNTER`; use `condition` only as a
design factor for PCOUNTER and keep its omission/guess terms fixed.

### 1. Direct native likelihood attribution

Use `WorkingTests/bench_efficiency_targets.R` with profiling off for speed and
on for explanation:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for model in RDM LBA BAwL DDM; do
  MODEL="$model" TRIALS=64,600,2000 CELLS=1,4,64 \
    PARTICLES=1,5,25,100 REPS=9 BRANCH=default \
    Rscript WorkingTests/bench_efficiency_targets.R
done
MODEL=PCOUNTER TRIALS=64,600,2000 CELLS=1,4,64 \
  PARTICLES=1,5,25,100 REPS=9 \
  BRANCH=constant,varying,k_varying,continuous,mixed,truncation,omission \
  Rscript WorkingTests/bench_efficiency_targets.R
```

Run `k_zero` and `k_1023` separately for PCOUNTER. Run the 255/256/257/512
cell cliff separately for each model whose mapper supports that design. The
direct gate is attribution only: it says which kernel or mapper is worth
testing, not that a user-facing fit became faster.

Required direct outputs are median/MAD elapsed and CPU time, informative
finite likelihoods, repeated checksums, rows/cells, route-independent native
call counts, and (for PCOUNTER) branch, K, Stirling/rising-term,
preparation/RT counters.

### 2. Fixed particle-step throughput

Use the same prepared sampler state and proposal cloud for every build. Measure
200 calls to the particle-step function after three warm-up calls, with no
adaptation and no convergence checks. Run:

| Subjects | Trials/subject | Particles | Core budgets |
|---:|---:|---:|---:|
| 1, 8 | 200 | 25, 100 | 1, 4, 8, 16, 32 |
| 2, 4, 8 | 200 | 25, 100 | 1, 4, 8, 16, 32 |

For each model, compare serial, explicit `r_cores=2`, and the default
`r_cores=1` routes. Require bit-identical particle results and record the
active worker/BLAS budget. This isolates proposal, likelihood, scheduling and
serialization costs before a full fit adds adaptation and diagnostics.

### 3. Short end-to-end sampling

Use `WorkingTests/bench_hierarchical_regression.R`, with fixed stage lengths
and `search_width=1`. The new overrides keep the PCOUNTER proof-of-concept
cheap while preserving matched work:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export EMC_CASES=smoke,cores,particles EMC_TRIALS=64 EMC_PARTICLE_FACTOR=5
export EMC_ITER_SCALE=0.1 EMC_STEP=10 EMC_REPS=5
mkdir -p bench-results
for model in RDM LBA BAwL DDM PCOUNTER; do
  EMC_MODEL="$model" EMC_LABEL="post-$model" \
    EMC_OUT="bench-results/post-$model.rds" \
    Rscript WorkingTests/bench_hierarchical_regression.R
done
```

Repeat the same loop against the pre-change build with `post-` changed to
`pre-`. To print paired comparisons while producing the post-change result,
also set `EMC_BASELINE=bench-results/pre-$model.rds` for that model.

For the models that pass the short gate, repeat with the ordinary workload:

- 8 subjects, approximately 8 and 64 parameters, 200 and 600 trials;
- particle factors 25, 50 and 100;
- 1, 4, 8, 16 and 32 cores per chain;
- one and two chains with the per-chain budget held fixed.

Use the same case names and output files for pre/post builds so the comparison
is paired. PCOUNTER must appear in the same summary table as the other models;
do not report it as a separate success criterion.

### 4. Mapper and route stress cases

Only after the short end-to-end run is clean, run 600-trial cases at 255, 256,
257 and 512 cells, plus a continuous covariate with one cell per trial. These
cases are diagnostic for the mapper and scalar-reader work. They must report
the same likelihood checksum, sampler checksum and route identity as the
narrow controls. A speedup that comes from falling back to a different route
is not credited.

## Decision rules

1. **Noise gate.** A result must exceed the A/A median noise estimate and be
   repeatable in both run orders. Treat changes smaller than 2% as noise unless
   the A/A spread is lower and the result is reproduced in at least seven
   repeats.
2. **Direct-kernel gate.** A model-specific arithmetic change needs at least
   10% improvement in a realistic direct case, no more than 2% regression in
   its narrow/continuous control, and unchanged per-trial likelihoods.
3. **Sampling gate.** Promote a change only if complete-fit wall time improves
   by at least 5% and worst relevant ESS/sec does not regress, with no route,
   memory, liveness or reproducibility failure.
4. **Generalisation gate.** A PCOUNTER result is a proof of concept. If its
   measured preparation is not reusable across RDM/LBA/BAwL/DDM under the same
   contract, remove the PCOUNTER-specific production optimisation. If the
   pattern is reusable, implement and benchmark it model by model; do not
   leave PCOUNTER on a privileged path.
5. **Reversion gate.** If a model fails either the direct or sampling gate,
   revert its specialisation and retain only shared instrumentation that is
   demonstrably off-by-default and numerically inert.

## Deliverables

For each build, retain the raw `.rds`/text output, commit and environment
metadata, the A/A noise summary, paired pre/post tables, and a short decision
record per model. The final claim should be phrased as total ESS/sec and total
wall time by model, not as a PCOUNTER-only kernel result.

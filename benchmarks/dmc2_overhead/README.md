# DMC2 overhead baseline: fixture runner

`run.R` is a small fresh-process runner for building the Step 0 baseline in
`plans/dmc2_overhead_plan_2026-11-11.md`. A fixture is an R script defining
`run_fixture(seed)`; it must construct its data and fit and return the fit
and any fixture-specific measurements. Run each seed and repetition in a
separate process, with a distinct output path:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript benchmarks/dmc2_overhead/run.R fixture.R run-1.rds 101
```

The output records the fixture checksum, seed and RNG states, elapsed and R
process CPU times, session/package versions, thread settings and fit. Fixture
scripts should explicitly set `emc2.sampler_profile` and `emc2.fit_host` for
controlled profiles; use the default fit host for end-to-end runs. Use the
repository's benchmark lock for timed runs (see `benchmarks/allocator/README.md`).
R's `system.time` here does **not** account reliably for the whole process tree;
collect tree CPU, forks, faults and memory separately. This runner alone does
not satisfy Step 0: the specified fixtures, repeated measurements, precision
and posterior gates, and process-level instrumentation still need to be built
and recorded before accepting optimizations.

`psis_probe.R cold|warm` measures three concurrent child PSIS calls, launched
from a fresh parent. Run each observation under the benchmark lock with BLAS
threads fixed. On 2026-11-11, three cold process runs took 0.306–0.341 s wall;
three parent-warmed runs took 0.013–0.014 s. Each child returned the same
Pareto k (-0.2286573658) and the parent RNG state was unchanged. This isolates
startup; it is not a rebuild-block or full-fit benchmark.

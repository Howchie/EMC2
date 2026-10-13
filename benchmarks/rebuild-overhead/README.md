# Rebuild overhead measurements

Use an isolated baseline checkout and private R libraries. Use nine cores
(3 chains × 3 cores) for further runs, with the same CPU affinity, input state
and seed for each pair. The user requested relative
comparisons alongside other work, so these runs do not acquire the shared lock. Scripts print the loaded package path.

```bash
EMC_LIB=/path/to/lib R_LIBS=/path/to/lib \
EMC_IN=/path/to/saved-fit.rds STAGE=burn PROF_OUT=/path/to/profile.rds \
CORES=9 taskset -c 23-31 \
Rscript benchmarks/rebuild-overhead/profile.R
```

`profile.R` refreshes a saved state's Laplace fits and proposals. Set `MODE=run`,
`ITER=22`, and `STEP=2` for eleven sampling blocks. Its RDS output contains:

- `profile`: exclusive component wall, master CPU and worker CPU; nested
  inclusive wall is also retained. Exclusive components can be added without
  counting their children's work again.
- `jobs`: actual worker counts, requested width, worker CPU and maximum RSS
  for each fork dispatch. A component without a dispatch runs in the caller.
- `blocks`: the package's ordinary block records, with complete proposal
  accounting on the updated build.
- `memory`: per-process RSS at each block boundary, including pool templates
  adopted by a subreaper.
- `state` and `seed`: output sampler state and caller RNG state. These can
  contain private data; keep the output on the fixture's machine.

Worker CPU travels with each fork's result. `proc.time()` child CPU alone is
unsuitable for phase attribution: reaping can occur after `mclapply` returns.
Nested fork work is included once. Sampling workers that use the persistent
pool do not return through this wrapper; sampling CPU in this profiler is
partial. The non-sampling components use the measured fork CPU. Instrumentation
adds a small amount of overhead and may affect the allocator's timing probes.

`busy_cores` is CPU seconds divided by exclusive wall seconds, including the
master; it measures average useful CPU occupancy. `idle_core_seconds` is
`max(0, requested_cores * wall - master_cpu - worker_cpu)`. Both depend on CPU
contention and on the configured capacity. Worker counts in `jobs` describe
actual processes; they are not a claim that those processes worked for the
whole component.

`screen.R` sources candidate R implementations over the baseline namespace,
then compares each changed component's complete result and RNG state using
`identical()`. It runs paired baseline/candidate measurements twice. Run it
from a frozen copy of the candidate R sources and pass saved fits as arguments.
Raw data fixtures need an initialized checkpoint first. `SCREEN_OUT` selects
its output file.

`identity.R` runs two blocks and saves draws, RNG streams and discrete plans.
Run it separately with the baseline and candidate libraries, then compare its
outputs with `identical()`. `HOLD_ALLOCATION=1` uses the existing sampler tests'
fixed move allocation, allowing a draw comparison without wall-time probes
selecting different move counts. This is validation scaffolding only.

`monitor.py out.csv <command...>` polls descendant RSS every 0.2 seconds.
The block-boundary records supplement it for pool processes reparented to a
subreaper. RSS includes shared pages; summing RSS across forks overstates unique
physical memory. Compare a worker's readings across blocks.

For saved-state sample rebuild profiles, `SAMPLE_REBUILD=first` clears the next
rebuild gate and `SAMPLE_REBUILD=next` advances the transition counter to that
gate. These are component profiles on the supplied history, not a run of the
intervening sample transitions.

`diagnostics.R` isolates the array-slicing change. `diagnostics-installed.R`
runs complete stopping checks against one installed library, saving the full
results and RNG for a comparison with another library. Set `DIAG_OUT` to the
output path. Both use `mu`, `sigma2` and `alpha`, supported by the standard and
diagonal-gamma fixtures. `REPEATS=12 KEEP_RESULTS=0 MEMORY=1` runs repeated
checks with collection between checks and records worker RSS and master file
descriptor counts. The master RNG and complete check result must repeat exactly.
Raw result states remain private.

Large diagnostic groups use Unix forks when there are at least eight groups
and a million scalar draws. Each existing MCMC group stays intact, including
its multivariate R-hat calculation. Public diagnostic summaries default to one
core; sampler stopping checks supply the current core budget. Custom summary
functions retain their serial evaluation and argument forwarding.

`decisions.R` compares complete stopping results and RNG with `identical()` on
a shared input. Run from the candidate source directory against the baseline
library, pass fixture paths, and set `DECISIONS_OUT`. This avoids comparing
recreated model/prior function environments across separately deserialized
files. `summarize.R /path/to/raw-results /path/to/exports` exports timing and
process-memory tables without sampler states.

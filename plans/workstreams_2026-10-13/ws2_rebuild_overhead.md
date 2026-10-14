# WS2: Cut the time between sampling blocks

Read `README.md` in this directory first. It covers the shared rules (worktree, private library, machine lock, comment style, exactness) and the current state of the code.

## Goal

Cut the wall time `run_emc` spends outside sampling (between blocks) without changing what the sampler does.
- On the user's design, the target is to bring non-sampling time per burn block from about 27-35 s to under 10 s.
- Small models must not get slower in absolute terms.

The default is changes that leave every draw identical for a given seed. Anything that changes the draws must keep the transition statistically exact, and must be argued, tested and reported separately.

## Why

The block profile of a burn on the user's design (`lba7r`: LBA, 106 subjects, p = 20, 3 chains x 9 cores) shows 27-35 s of every 115-145 s block outside sampling. Most of it runs serially in the master R process while all 27 workers sit idle. That is about 20% of burn wall time. In the sample stage, rebuilds follow a doubling schedule and are rare, so the gain there is smaller but real for short fits.

Measured on `lba7r_after_adapt.rds` with `harness/time_prop.R`, one burn rebuild (Rlib built from `0a5773dd`):

| Component | Time | Where |
|---|---|---|
| `.emc_laplace_refresh` | about 9-10 s | `block_prepare`, `R/fitting.R` `run_emc` |
| `.emc_alloc_rebuild` | about 5 s | allocator; **owned by WS1**, measure but do not edit |
| `.emc_shear_orbit_steps` | about 5.9 s | `R/shear_move.R`; one batched `auto_mclapply` over (chain, coordinate) jobs |
| rest of `.emc_shear_set_plans` | about 3 s | `R/shear_move.R`; frames, curve gates, atlas, validation |
| `.emc_check_likelihood_fits` | about 1.3 s | |
| `.emc_collapsed_set_plans` | about 1.2 s | |
| `create_chain_proposals` | about 1 s | |
| `create_eff_proposals` | about 0.2-0.8 s per call, several calls | |
| `check_progress` | about 2.5 s | R-hat and ESS for the stop rule |

## Work

### 1. Measure properly first

Build a profiler that attributes all non-sampling wall time in `run_emc`'s block loop to named components, and records how many cores each component kept busy. Report wall, master CPU and worker CPU, so serial phases are obvious.

The existing block profile's `block_proposals` covers only the chain-plan part of `add_proposals`, not the allocator or the efficient proposals. Fix that field as part of this work.

Profile the burn, adapt, and first two sample-stage rebuilds on:
- `lba7r` (all three saved stages);
- `fixture_nback31.rds`;
- T30 (`benchmarks/subject-hmc/fixtures/fixture_T30.rds`);
- Forstmann (package data; small model, where fixed overheads dominate).

`harness/prof_run3.R` and `harness/time_prop.R` are starting points.

### 2. Rank the opportunities by idle core-seconds, not seconds

A 3 s phase that idles 27 cores costs more than a 5 s phase that uses them all. For each component, classify it as:
- inherently serial;
- parallelisable across chains (three chains' work is independent);
- parallelisable across subjects or coordinates;
- redundant (recomputed when an input has not changed, or computed when its result is unused).

### 3. Implement, in order of payoff

Candidates to examine. These are not prescriptions; verify each against the code.
- **Laplace refresh:** check whether it already uses the worker pool, and how it balances load (subjects differ widely in cost).
- **Per-chain work:** much of `add_proposals` and plan building is per chain, but runs as a loop in the master. Some inputs pool across chains (variances over every chain's window), so separate the pooled steps from the chain-local ones before parallelising.
- **Orbit steps and curve validation:** these evaluate many likelihood rows; they could use the claim-counter subject queue (`project_emc2_shear_queue.md`, `project_emc2_worker_pool.md`) instead of forked `auto_mclapply` jobs.
- **`check_progress`:** R-hat and ESS over growing windows; check for recomputation and needless copying.
- **Redundant work:** fits checked again when unchanged, proposals rebuilt when the window has not moved, and so on.
- **Overlapping rebuild with sampling:** proposals built from draws up to one block back, while the next block runs. This changes the draws. It is valid adaptive MCMC only under conditions you must state (adaptation a function of the past only, diminishing adaptation in production). Do not build it: write it up as an option with the expected saving and the argument, for the user to decide.

### 3b. Not yours: serial work inside iterations

The ensemble update's in-chain cost, the Gibbs step and the other per-iteration serial steps belong to WS5 (`ws5_ensemble_cost.md`). If your profiler shows in-iteration costs, report them for WS5.

### 4. Verify each change

- **Identity:** for a fixed seed and state, the draws after one or two blocks are bit-identical to the baseline build. Where parallel order changes floating-point summation, show the differences are rounding-level and that no discrete decision (plan, gate, curve selection) changes.
- **Memory:** worker RSS must not grow across blocks. Run 10+ blocks and record per-worker RSS. A nested model-definition leak once crashed the machine, so anything new sent to workers must not accumulate.
- **Timing:** paired, under `bench_lock.sh`, on an idle machine. Same seed, same state, baseline and new build back to back, repeated twice.

## Deliverables

1. The profiler (with the corrected `block_proposals` field) and per-fixture profiles before and after, in `ws2_report.md`.
2. Commits on `ws2/<topic>`, one per optimisation, each with its identity and memory evidence and its measured saving.
3. A table: component x fixture, before and after, wall and idle core-seconds; plus whole-fit wall saved for a burn plus adapt on `lba7r`.
4. A written note on the overlap option: expected saving, validity conditions, risks. Not built.

## Boundaries

- Do not edit `R/move_allocation.R` (WS1). Report allocator costs, and pass any idea for them to WS1 via your report.
- Do not change the particle step (WS3) or sampler behaviour.
- If an optimisation would alter adaptation timing, plans or gates, stop and write it up instead.
- All timing runs go through `bench_lock.sh`. Stop at once if the user starts a fit.

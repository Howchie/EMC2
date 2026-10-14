# Sampler wrap-up: four parallel workstreams

Five independent agents, each with its own prompt in this directory:

| # | File | Owns | Depends on |
|---|---|---|---|
| WS1 | `ws1_allocator.md` | `R/move_allocation.R`, `src/alloc_need.cpp` | WS4's mock geometries, once available (not a blocker) |
| WS2 | `ws2_rebuild_overhead.md` | Speed of everything between sampling blocks, *except* `R/move_allocation.R` | none |
| WS3 | `ws3_particle_step.md` | The particle step's count and adaptation: `new_particle()`, `update_pm_settings()`, `update_local_steps()`, `particle_factor`, `tune$*` in `R/sampling.R` and `R/fitting.R` | WS4's mock geometries, once available (not a blocker) |
| WS5 | `ws5_ensemble_cost.md` | The ensemble update's speed, then its pool size and rounds (`R/ensemble_move.R`); the Gibbs group step's speed; the iteration profile fields | WS4's geometries for Phase 2 (not a blocker) |
| WS4 | `ws4_test_geometries.md` | New code under `benchmarks/geometry/`, plus the new mock model structures (R model definitions and C++ kernels, REXG-style) | none; the others consume it |

Paste the matching prompt file (in full) as the agent's task, then this README.

## Why these and not an ablation matrix now

The sampler now has many mechanisms that EMC2 never had:
- group ensemble update;
- shear (interweaving) rounds with curved, state-indexed frames and orbit-measured steps;
- the collapsed group move;
- subject HMC;
- the minimax move allocator;
- an optimiser-seeded burn start.

Each was built against a bottleneck found on one fixture. Fixing it usually exposed the next bottleneck, and no fixture isolates one geometry. So an ablation matrix today would be expensive and uninterpretable. These four workstreams make the pieces auditable first:
- **WS4** builds test problems whose geometry is *known*.
- **WS1** checks that the allocator's predictions are calibrated.
- **WS2** removes avoidable serial overhead.
- **WS3** puts the particle count on a principled, general footing.

The ablation comes after, on WS4's set.

## State of the code (2026-10-13)

- Repository `/data/work/EMC2_dev_oo`, branch `playground`, HEAD `0a5773dd`.
- The working tree has unrelated tracked deletions under `tools/` and untracked `plans/` files. These belong to the user: never commit, revert or delete them.
- Legacy baseline: `/data/work/EMC2_cens_trunc2` (no group moves, still has the bugs in `plans/sampler_bugs.md`). Its `make_emc` floors rts to the 1/60 s grid with no midpoint correction. Compare posteriors with it only on pre-floored data and `rt_resolution = NULL` (see the Forstmann facts below).

**Forstmann (LNR, package data), seed 1, 3 chains x 9 cores, same build:**
- Minimum bulk ESS is over 3,000 sample draws.
- Lag-k and tail behaviour come from `fit_arm.R`'s diagnostics.

| Arm | Wall | CPU | mu | lsd | alpha bulk / tail |
|---|---|---|---|---|---|
| full sampler | 90 s | 819 s | 2,947 | 2,693 | 2,308 / 2,176 |
| group moves off (`nogroup`) | 50 s | 520 s | 1,627 | 611 | 478 / 163 |
| also no mode-start burn (`bare`) | 46 s | 505 s | 2,041 | 718 | 989 / 667 |

Single seed: on an easy model the mechanisms cost about 1.7x wall and CPU and buy 2.3-4x worst-coordinate ESS. On identical data the legacy and current likelihoods agree exactly, and all group posteriors agree within |z| < 2.

Where the extra cost goes. These are profiled seed-1 runs (`PROFILE=1`), sample stage, chain-seconds summed over the 3 chains x 1,000 iterations:

| Arm | Total | Particle | Gibbs | Shear | HMC | Other (in-chain R) | Worst ESS (lsd) |
|---|---|---|---|---|---|---|---|
| no group moves | 63 | 54 | 3.6 | - | - | 5.5 | 611 |
| ensemble only | 108 | 63 | 7.7 | - | - | 37.7 | 2,764 |
| full | 124 | 64 | 8.6 | 6.7 | 3.0 | 41.8 | 2,693 |

- On Forstmann the ensemble update accounts for about three-quarters of the extra cost, and for nearly all of the ESS gain.
- Its own step is about 11 ms of serial R per iteration in each chain process, while that chain's workers wait. That cost does not scale with likelihood cost, so cheap-likelihood models pay the most, proportionally.
- It runs every iteration (`emc2.ensemble_rounds = 10`, `emc2.ensemble_pool = 32`) and is not priced by the allocator. It is the core group mechanism and always runs. Whether it runs is not allocated or tuned; its pool size and rounds are constants with no recorded derivation, and WS5 makes it cheap and then sizes them on principle.

**User's real design (lba7r_gL: LBA, 106 subjects, p = 20), one burn rebuild:**
- `add_proposals` takes 18.7 s:
  - allocator about 5 s;
  - orbit steps 5.9 s;
  - other shear planning 3 s;
  - likelihood-fit checks 1.3 s;
  - collapsed plans 1.2 s;
  - chain proposals 1 s.
- On top of that: the Laplace refresh is about 9-10 s (in `block_prepare`) and `check_progress` about 2.5 s.
- Together that is 27-35 s of non-sampling time per ~120-140 s burn block, mostly serial in the master while 27 workers idle.

## Shared resources

- `/data/work/EMC2_workstreams/fixtures/`:
  - `lba7r_after_preburn.rds`, `lba7r_after_adapt.rds`, `lba7r_sample_stage.rds`: the user's design at three stages. Real data; keep it on this machine.
  - `lbaE_fit.rds`: an LBA with a trend, a completed fit.
  - `forstmann_full_fit.rds`: a completed full-sampler Forstmann fit (LNR).
  - `fixture_nback31.rds`, `fixture_nback.rds`: real N-back data with its prior and type (diagonal-gamma).
  - Synthetic fixtures in the repo: `benchmarks/subject-hmc/fixtures/` (T30 Wald RDM with correlated group means, lbaA, lbaT, lbaE).
  - Caution: lbaT leaves A unidentified (multi-optima), so it is not a test of mixing.
- `/data/work/EMC2_workstreams/harness/`:
  - `time_prop.R` times each `add_proposals` component on a saved state by wrapping namespace functions.
  - `prof_alloc.R` profiles the allocator rebuild.
  - `prof_run3.R` runs stages with the block profiler on.
  - `fit_arm2.R` is `benchmarks/subject-hmc/fit_arm.R` with an `ARM_R=<file.R>` hook for options.
  - `compare.R` compares ESS and posterior agreement.
  - `prof_ens.R` profiles the ensemble move over a 200-iteration sample block.
  - `runtests.R` runs `EMC_LIB=... Rscript runtests.R test-x.R ...` from the repo root.
- Saved fits have their duplicates stripped. Call `EMC2:::restore_duplicates()`, then `restore_custom_kernel_pointers()`, before calling internals. Without that, `n_pars` is NULL.

## Rules every workstream follows

**Isolation**
- Work in your own git worktree on your own branch:
  ```
  git -C /data/work/EMC2_dev_oo worktree add /data/work/EMC2_ws<N> -b ws<N>/<topic> playground
  ```
- Install into your own private library, with both variables set:
  ```
  R CMD INSTALL --no-docs --no-test-load -l <lib> .
  EMC_LIB=<lib> R_LIBS=<lib> Rscript ...
  ```
- `library(EMC2)` without `lib.loc` loads a stale system install. Every script prints `find.package("EMC2")`, and you check it.
- Run `rm src/*.o` after editing any header: the makefiles ignore header dependencies.
- Never commit to `playground`. Deliver a branch plus a report; the user merges.

**The machine (32 cores; one heavy job at a time)**
- Any job using more than 4 cores, and *any* job whose wall or CPU time you will report, runs under the lock:
  ```
  /data/work/EMC2_workstreams/bench_lock.sh <command...>
  ```
- Fits use 3 chains x 9 cores, pinned with `taskset -c 0-26`. Never run several smaller fits side by side: the allocator prices moves by wall time, so concurrent load changes its plans, not just the timings.
- If the user starts a fit (an R process you did not start, using many cores), stop your heavy jobs at once and wait. Load from someone else faked a 2x regression before.
- Never edit a script, or reinstall a library, that a running job reads. Rscript reads files incrementally. Copy the script, and give each launched fit its own frozen library copy.
- Workers must not accumulate state between tasks: a nested model-definition leak once exhausted RAM and crashed the machine. Anything you send to the worker pool goes through the existing broadcast file and claim-counter paths.
- Short screens before long runs: mechanism-level checks on saved states first, then at most one short continuation (about 400 iterations), and stop-rule fits only once a change shows an obvious effect cheaply. Never queue hours of fits to explore.

**Code**
- One default sampler (fullc: ensemble + shear rounds + subject HMC + collapsed move, allocated). No options framework.
- `options(emc2.*)` on/off switches are scaffolding for testing only.
- Never remove or disable a mechanism because a short, noisy run looked flat. A mechanism that never builds or engages is a bug to find, not a verdict. Zero acceptance is acceptable off its target, but every mechanism must be shown to accept, with counts, on a case it was designed for.
- No user-facing warnings as a fix for poor mixing; fix the sampler.
- Comments say what the code does and why. No fixture names, benchmark numbers, dates or investigation history in code comments.
- Statistical exactness is non-negotiable:
  - every adaptive quantity is fixed within a sampling block and changes only between blocks, as a function of the past;
  - any new move gets an invariance test, including a deliberately broken control that the test must catch.
- Only the C++ likelihood path is authoritative; the R likelihood path is outdated.
- A resumed checkpoint applies its recorded `production_options` over the session's options. `options()` set before a continuation can therefore be silently ignored. Check the chain's `production_options`.
- Tests: `NOT_CRAN=true`, plus `EMC2_TEST_LEVEL=full` for the full level, run from the repository root with the shared `runtests.R`. The suites that touch the sampler are `test-move-allocation.R`, `test-shear-move.R`, `test-shear-standard-invariance.R`, `test-collapsed-move.R` and `test-subject-hmc.R`. Run the sampling and fitting tests relevant to your change too.

**Reporting**
- Every claim states its evidence: fixture, seed count, build, and the measured numbers.
- Single seeds are "indicative", never "established".
- Say plainly what did not work or was not tested.
- Write your report to `plans/workstreams_2026-10-13/ws<N>_report.md` in your worktree, and commit it on your branch.
- Background reading: the memory files in `/home/ubuntu/.claude/projects/-data-CensoringTest/memory/`, indexed in `MEMORY.md`. They reflect what was true when written; the code is the source of truth, so verify a memory against the code before relying on it.

# WS1: Is the move allocator making the right decisions?

Read `README.md` in this directory first. It covers the shared rules (worktree, private library, machine lock, comment style, exactness) and the current state of the code.

## Goal

Establish, with evidence, whether the allocator (`R/move_allocation.R`) buys the right moves:
- its predictions of how much each move speeds up each coordinate are calibrated;
- its plans are stable against measurement noise;
- it engages where a move is needed and abstains where it is not.

Fix what the evidence shows to be wrong, and make what remains simpler where that loses nothing. This is an audit with fixes, not a redesign. Propose a redesign only if calibration shows the model is wrong in a way local fixes cannot cure, and stop to report before building it.

## Why

The allocator decides how many shear rounds, collapsed moves, subject HMC transitions and particle steps each chain runs. It has carried most of the sampler's recent bugs, each of which silently degraded mixing rather than failing:
- a regime holding every particle step (305d47fc);
- immovable units flooring the target so a chain bought nothing (55b29f4d);
- phantom credit for the collapsed move (d8357bc0);
- measurements carried across kernels, which killed the collapsed move for a whole fit;
- plan flips from noisy 100-draw windows;
- 13-20 s of serial R per rebuild (0a5773dd).

It has many constants. Nobody has checked its *predictions* against what the next block actually does. Until that is done, nobody can say whether an observed result reflects a mechanism or the allocator's misjudgement of it.

## What the allocator is (verify each point in the code; this is a summary)

- **Measurement, at each rebuild, per chain** (`.emc_alloc_measure`):
  - per coordinate, a lag-one correlation `rho = 1 - ESJD / (2 var)` of "the rest of the iteration";
  - `var` is pooled over every chain's window draws (all production draws in the sample stage);
  - an upper confidence bound on ESJD makes `rho` a lower bound;
  - coordinates are group mu along the shear directions, group log sds, and whitened subject alphas.
- **Move factors:**
  - shear round factor `g_k`;
  - collapsed move factor;
  - subject HMC factor `h_i`;
  - particle step factor `p_i`;
  - all measured from the moves' own jumps;
  - blended with pooled or prior references by evidence weight (`.EMC_ALLOC_ROUND_EVIDENCE = 20`);
  - decayed by 0.5 per production rebuild.
- **Composition:** k shear or collapsed moves give a factor `1 - c(1 - a^k)`. The saturation `a` is measured per coordinate (`.emc_alloc_saturation`), geometric until identified.
- **Subjects:** particle rate `q` on a grid of eighths, HMC rate `r`, jointly priced under coupled selection (`.emc_alloc_dev_factor`). There is a floor: no coordinate may end up slower than under the particle step alone.
- **Objective:** minimise `max over coordinates of tau * (C_rest + costs)` with `tau = (1 + rho) / (1 - rho)` (`.emc_alloc_solve`, `.emc_alloc_solve_units`).
- **Costs:**
  - measured wall time on the chain's reserved cores;
  - marginal per move, plus a fixed cost per family (`.EMC_ALLOC_FAMILIES`);
  - activation probabilities (`.emc_alloc_activation`).
- **Other rules:**
  - `.EMC_ALLOC_MIN_GAIN`: a plan must beat doing nothing by this margin;
  - exploration floors for under-measured coordinates;
  - HMC regimes (`emc2.alloc_hmc_regimes`, default on);
  - particle rates shared across chains as the maximum (`plan_chain(j, q_fixed)`, settle).

Background in memory: `project_emc2_graded_allocator.md`, `project_emc2_alloc_flip.md`, `project_emc2_alloc_floor_veto.md`, `project_emc2_particle_rate_alloc.md`, `project_emc2_regime_particle_bug.md`, `project_emc2_alloc_need_kernel.md`, `project_emc2_subject_hmc.md`. Also `plans/handoff_graded_allocator_2026-10-12.md`.

## Work

### 1. Write the model down from the code

Produce a short specification:
- every equation the allocator uses;
- every constant, with its value and where it comes from: derived, measured, or chosen;
- every default-on branch;
- which paths a default `fit()` actually takes.

Note any branch that no default fit can reach (dead code or scaffolding). This document is the baseline for everything else.

### 2. A replay harness

Given a saved fit state, run the measurement and solve exactly as a rebuild would, deterministically and without sampling. Output, per coordinate:
- the plan;
- the measured inputs;
- the predicted lag-one `rho` (and `tau`) under the plan;
- the predicted `rho` with nothing bought.

There was once a `solver_input` record for this (see `project_emc2_alloc_flip.md`); reuse the idea. Behaviour-preserving refactors must give identical plans on replay.

### 3. Calibration: predicted against realised

For each fixture, take a stationary saved state, force a rebuild, and run one block (100-200 iterations) under the plan. Then compare, per coordinate:
- the predicted `rho` under the plan;
- the realised lag-one autocorrelation in that block, using the same estimator (ESJD over pooled variance) *and* a standard estimator from the draws;
- the realised integrated autocorrelation time. Lag one may understate slow modes: the saturation story says moves at fixed residuals stop helping, so lag one can look fine while the IACT stays long.

Also run the counterfactuals that matter. For a coordinate where the allocator bought nothing, force a modest allocation and check that its prediction of "no gain" holds. Do the reverse where it bought a lot.

Fixtures:
- Forstmann (easy; little should be bought);
- T30 (correlated group means; shear and ensemble designed-for);
- `lbaE_fit.rds` (log-sd funnel and curved frames);
- `fixture_nback31.rds` (many subjects, p = 21);
- `lba7r_sample_stage.rds` (the user's design).

When WS4 publishes Gaussian mock geometries (see `ws4_test_geometries.md`; it announces them in `/data/work/EMC2_workstreams/geometry/STATUS.md`), use them. Their true autocorrelations under simple kernels are computable, so they turn calibration into an exact check.

Summarise calibration per coordinate kind: group mean directions, log sds, subjects. Report bias and spread, and say which kinds the allocator systematically over- or under-credits.

### 4. Decision stability

Bootstrap or block-resample the measurement window and re-solve. How often does the plan change materially, and which inputs drive the flips? A plan that flips on noise wastes the moves it buys. Quantify this before proposing any smoothing.

### 5. Known-answer decision tests

Extend `tests/testthat/test-move-allocation.R` with unit sets whose optimal allocation is known by brute force over a grid oracle. One exists for simple cases; extend it to the joint particle-rate/HMC pricing, the collapsed joint grid, floors, regimes and `MIN_GAIN`. These pin behaviour for the fixes that follow.

### 6. Fix and simplify

Fix defects the calibration or decision tests demonstrate, each with a test.

Remove a constant or branch only when:
- the specification shows it is unreachable or redundant, *or*
- calibration on several fixtures shows it changes no decision.

State which in the commit message. Do not tune constants to one fixture.

### 6b. The ensemble is out of scope

The group ensemble update always runs, every iteration. Whether it runs is not for the allocator to decide. Its speed, pool size and rounds belong to WS5. If your calibration shows that the allocator misjudges the other moves *because* the ensemble's contribution is folded into "the rest of the iteration", report that.

### 7. Cost

The allocator rebuild is now about 5 s on lba7r (`harness/time_prop.R`, `harness/prof_alloc.R`). It is about 24k calls of `.emc_alloc_unit_need` at roughly 80 µs of R overhead each, plus activation and settle.
- Pricing whole units (all subjects x all targets) in one C++ call is the obvious next step.
- Keep results equal up to the package's fast-math rounding, and verify on every call of a real rebuild, as was done for 0a5773dd: run the old R implementation alongside the new one and compare.

## Deliverables

1. `plans/workstreams_2026-10-13/ws1_allocator_spec.md`: the specification from step 1.
2. The replay and calibration harness under `benchmarks/allocator/`, with a README on how to run it.
3. Calibration tables per fixture and coordinate kind, plus the stability results, in `ws1_report.md`.
4. Fixes and simplifications as separate commits on `ws1/<topic>`, each with tests. All sampler suites must pass at the default and full levels.
5. A one-page verdict:
   - is the allocator's model adequate;
   - what it gets wrong;
   - what was fixed;
   - what remains;
   - what WS4's test set must contain to judge it properly.

## Boundaries

- You own `R/move_allocation.R` and `src/alloc_need.cpp`.
- The particle *count* and its adaptation belong to WS3. The particle *rate* `q` is yours. If WS3 changes how the particle step performs, your measurement of `p_i` follows automatically; note any interaction you find in your report.
- Speed changes outside `R/move_allocation.R` belong to WS2.
- Heavy runs go through `bench_lock.sh`; calibration continuations are short (one block or a few). Do not run stop-rule fits.

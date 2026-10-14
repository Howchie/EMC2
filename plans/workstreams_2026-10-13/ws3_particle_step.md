# WS3: A principled, general particle count for the subject step

Read `README.md` in this directory first. It covers the shared rules (worktree, private library, machine lock, comment style, exactness) and the current state of the code.

## Goal

Replace the particle step's ad hoc sizing with a rule derived from:
- what the step achieves per particle;
- what each particle costs.

The rule must:
- adapt from measured quantities;
- carry no constant chosen to suit a dataset;
- behave sensibly across parameter counts, data informativeness, proposal quality and likelihood cost.

First find out where the particle step's cost goes, and remove avoidable overhead. A cheaper step changes the optimum.

## Why

The particle step (a conditional importance-sampling update of each subject's alpha) is most of the sampler's work on many models. On N-back it is about 80% of the iteration.

Earlier screens (N-back 31, three seeds per arm) found this:
- N = 60 particles gave about 1.7x the median-direction ESS per CPU-second of the default (about 229).
- The slowest direction overlapped, and log-sd minima improved 7x.
- See `project_emc2_benchmark_protocol.md` in memory.

That says the default count is far from optimal on at least one real model. But nobody knows why the default is what it is, so no one change can be trusted to generalise.

## The current scheme (verify in the code)

**Initial count:** `round(particle_factor * sqrt(n_pars))`, with `particle_factor = 50` by default (`R/fitting.R`, `run_stages`). No derivation of 50 or of the square root is recorded; treat both as legacy heuristics.

**Adaptation** (`update_pm_settings()`, `R/sampling.R`):
- The count moves toward a target weight ESS of the independence component: `tune$target_ESS = 2.5 * sqrt(n_pars)`.
- It moves by `(target / ess)^tune$ESS_scale` per iteration, with `ESS_scale = 0.05`.
- It is bounded to `[25, 1.2 x initial]` (`tune$max_particles`).
- It adapts only when `gd_good` and the mixture has more than three components.
- Tuning stops in the sample stage (`plans/sampler_bugs.md` item 3).

**Mixture:** independence components (group-conditional prior, chains/efficient proposals) and state-centred random-walk components.
- Each random walk has its own conditional-IS step, steered to a share target of Phi(-2.38 / sqrt(2)) (`update_local_steps()`).
- Mixture weights follow each component's expected squared jump per particle (`plans/sampler_bugs.md` items 1, 9, 10).

**Particle rate:** the allocator (WS1) decides how often each subject's particle step runs (rate `q`), and prices it from measured worker CPU. The allocator keeps a floor: no coordinate may end up slower than under the particle step alone.

**Cost:** on T30 (Wald RDM, 40 subjects) a particle step cost about 5.5 ms of worker CPU per subject, against about 0.5 ms of likelihood evaluation. Most of it was overhead outside the likelihood. That figure is old (`project_emc2_subject_hmc.md`); re-measure it.

## Work

### 1. Where does a particle step's time go?

Profile one subject's particle step on Forstmann (LNR), T30, `fixture_nback31.rds` and `lba7r`. Break it down into:
- proposal draws (mixture sampling, multivariate normal draws);
- likelihood rows;
- density evaluations of the proposals, including the auxiliary-centre terms;
- weights and resampling;
- R overhead and data movement.

Express each as fixed cost per step and cost per particle (`c0 + c1 N`). If the step's overhead dwarfs its likelihood work, cut it first. Typical targets are moving the density and weight computations into C++, avoiding repeated matrix factorisations and avoiding copies. Verify such changes as draw-identical for a fixed seed.

### 2. Theory: optimal N for this step

Derive the expected movement per unit cost as a function of N and of measurable properties of the target and proposals. Start from the literature on iterated conditional importance sampling, multiple-try Metropolis and i-CSMC: Andrieu, Lee and Vihola, and related results.

For a mixture of independence proposals, the probability of leaving the reference particle rises with N roughly as `(N-1) e / (1 + (N-1) e)`, where `e` is the per-particle relative efficiency (weight ESS / N). Expected squared jump per cost is then maximised at an interior N set by `c0 / c1` and `e`.

Work out the corresponding result for the random-walk components, whose share is pinned by their step controller, and for the mixture as a whole. State the assumptions and where they fail: heavy-tailed weights, curved posteriors, parameter count. Check the predictions against measurements in step 3 before building on them.

### 3. Measure, then build

On saved stationary states, sweep N for single subjects without running fits:
- the subject's ESJD in whitened coordinates;
- the move probability;
- the weight ESS;
- the CPU, against the theory's predictions.

Use subjects that differ in informativeness and curvature: the T30, lba7r and N-back fixtures, and WS4's mock geometries when they are published (`/data/work/EMC2_workstreams/geometry/STATUS.md`). Then design the adaptation:
- N per subject (or per chain), chosen from measured `e`, `c0` and `c1`;
- bounded only where the theory says the measurements are unreliable;
- no dataset-tuned constant.

Decide whether `particle_factor` survives:
- If it does, give it a meaning the user can reason about.
- If the rule replaces it, keep the argument accepted for backward compatibility, and document what it now does.

### 4. Group-level effect

Group mixing depends on how far subjects move each iteration, which the allocator does not fully model. A cheaper, weaker particle step can lower subject movement and slow the group even if per-subject efficiency rises.

So check the group effect on short continuations, about 400 iterations from saved states:
- lag-k decorrelation of group means and log sds, against the checkpoint variance;
- CPU.

Compare the current rule with the new rule on Forstmann, T30, `fixture_nback31.rds` and `lba7r_sample_stage.rds`. One seed each first; repeat only where a result is close. Report the allocator's chosen particle rates and HMC rates in both arms, because they will shift.

### 5. Exactness and tests

- The count must not change within a sampling block, nor in the sample stage except at rebuilds.
- Any change to the step's proposal densities needs the existing particle-step invariance tests to pass, or new ones (with a deliberately broken control).
- Run `test-subject-hmc.R`, `test-move-allocation.R`, the sampling tests and the shear identity tests, which force particle rates.

## Deliverables

1. The cost breakdown (step 1) and any draw-identical speedups, as separate commits on `ws3/<topic>`.
2. A note with the derivation and its empirical check (steps 2-3), in `ws3_report.md`.
3. The new adaptation rule with tests, as its own commit.
4. Continuation results (step 4): subject and group level, CPU, and allocator rates for both arms, per fixture. Single seeds are labelled indicative.
5. A recommendation for the user on `particle_factor` and the default.

## Boundaries

- You own the particle count, its adaptation and the step's own code in `R/sampling.R` (`new_particle`, `update_pm_settings`, `update_local_steps`, and the `tune` defaults), plus the `particle_factor` plumbing in `R/fitting.R`.
- The particle *rate* `q` and the allocator belong to WS1. If you find the allocator's floor or pricing misjudges a changed particle step, document it for WS1; do not edit `R/move_allocation.R`.
- Speedups to code between blocks belong to WS2.
- The ensemble update's pool rows (M - 1 per subject) ride on the same worker round as the particle step. WS5 sizes the pool by the same importance-efficiency argument you are developing. Share your derivation and your cost measurements of that round through your report, and read WS5's. Do not edit `R/ensemble_move.R`.
- Heavy runs go through `bench_lock.sh`. N sweeps on single subjects are light, but if you report their timings, they still go through the lock.

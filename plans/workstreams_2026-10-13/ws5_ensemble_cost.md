# WS5: The ensemble update: make it cheap, then size it on principle

Read `README.md` in this directory first. It covers the shared rules (worktree, private library, machine lock, comment style, exactness) and the current state of the code.

## Goal

Make the group ensemble update cost as little as possible, *without changing what it does*:
- same rounds and pool size;
- same proposals and the same acceptance;
- draws identical for a given seed, or equal up to floating-point rounding with no discrete decision changed.

Then apply the same treatment to the other serial in-iteration steps of the chain process, the Gibbs group update first.

The ensemble is the sampler's core group mechanism. It **always runs**, every iteration. Whether it runs is not up for allocation or tuning.

Its size is a different matter. The pool size `M` (`emc2.ensemble_pool = 32`) and the number of rounds `R` (`emc2.ensemble_rounds = 10`) are fixed constants with no recorded derivation. They do not change with:
- the number of parameters p;
- the number of subjects;
- how informative the data are;
- how well the pool proposal matches each subject's posterior.

**Phase 2** replaces them with a principled, general rule. Phase 1 must finish first, because the cost model Phase 2 optimises changes once the rounds are cheap. In Phase 1, every change preserves the draws.

## Why

On Forstmann (LNR, 19 subjects, p = 4), profiled:
- The ensemble update is about three-quarters of the full sampler's extra cost over a sampler with no group moves.
- It is also nearly all of the full sampler's ESS gain on that model. See the README table.

Its own step (`.emc_ensemble_move`, `R/ensemble_move.R`, called from `R/sampling.R`) costs:
- about 7.3 ms per iteration for one chain on an idle machine;
- about 11 ms with three chains running;
- about 45% of the chain process's sampled time.

All of that is serial R in the chain process while the chain's 9 workers wait. Its likelihood candidates are evaluated separately, as extra pool rows in the particle round, so they are not part of this cost.

The arithmetic involved is tiny. The time is R overhead spread across hundreds of small matrix calls per iteration. An `Rprof` of 200 sample iterations showed:

| Function | Time |
|---|---|
| `.emc_ensemble_draw_pools` | 0.32 s |
| `.emc_ensemble_log_g` | 0.19 s |
| `.emc_ensemble_log_q` / `.emc_ensemble_log_q_m2` | 0.14 s |
| `.emc_ensemble_draw_q` | 0.13 s |
| `.emc_ensemble_select_all` | 0.11 s |

The rest is in the move's body: `as.matrix`, `matrix`, `t`, `forwardsolve`, `solve`, `chol`, `colSums`.

This overhead does not scale with likelihood cost, so cheap-likelihood models pay the most, proportionally. On models with costly likelihoods it is a smaller share, but it is still idle worker time.

The profiling script is `/data/work/EMC2_workstreams/harness/prof_ens.R` (`EMC_IN=/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds`). It runs `fit(stage = "sample", iter = 200)` on one chain with `.emc_ensemble_move` wrapped for timing, under `Rprof`.

## Work

### 1. Profile properly

Measure per-iteration time of the ensemble move and its parts on:
- Forstmann (p = 4);
- T30 (p of about 10, 40 subjects);
- `fixture_nback31.rds` (p = 21, about 31 subjects, diagonal-gamma);
- `lba7r_sample_stage.rds` (p = 20, 106 subjects).

Cost will scale differently with p, subjects, pool size and rounds. Find out how, so the fix targets the regime that matters. Include the pool draws that feed the particle round (`.emc_ensemble_draw_pools`) and the score refresh (`.emc_ensemble_refresh_scores`).

### 2. Read the move and write down its arithmetic

Write down, per round, what is computed:
- the proposal distribution;
- the candidate densities;
- the subject-prior terms under each candidate group state;
- the selection.

Identify:
- what can be computed once per iteration rather than per round;
- what is recomputed redundantly (for example, factorisations of matrices that have not changed);
- what is independent across candidates and can be vectorised.

Do this for every supported group type (standard, diagonal, diagonal-gamma, and any others `.emc_ensemble_supported` accepts).

### 3. Implement

The most likely shape is one C++ function per iteration that runs all rounds: candidate densities, subject-prior terms, selection. It takes the pool and current state, and returns the chosen states plus the bookkeeping the sampler records. R keeps the orchestration.
- Random draws must come from the same streams, in the same order, as now, so a fixed seed gives the same chain. The sampler's RNG streams live in the `emc` object, so check how they are consumed.
- If that is impossible for one piece, state why and show the alternative is exact.
- Keep each change separately verifiable: port one part, verify identity, then the next.
- The package builds with `-ffast-math`, so expect rounding-level differences from R. Compare every call of real runs, running the old R implementation side by side, and show that no selection changes. This was done for the allocator in commit `0a5773dd`.

### 4. Then the Gibbs step and other in-chain work

`gibbs_step_standard` (and its `rwish`, factorisations and so on) showed similar overhead: roughly 1-3 ms per iteration on Forstmann. Apply the same approach to it and to any other serial in-chain step the iteration profile shows is significant. The profile has no fields of its own for ensemble or collapsed time; add them (`R/profile_schema.R` and the recording sites).

### 5. Verify (Phase 1)

- **Identity:** fixed seed and state, 200 iterations, old against new build: identical draws, or rounding-level differences with every selection identical.
- **Invariance tests:** the ensemble's existing tests must pass. Add a test that the C++ path and an R reference agree on random states for every supported group type.
- **Speed:** paired timings under `bench_lock.sh` on the four fixtures above, per-iteration ensemble time and whole sample-stage wall, before and after.
- Suites: `test-ensemble*` (whatever exists), `test-shear-move.R`, `test-collapsed-move.R`, `test-subject-hmc.R`, `test-move-allocation.R`, and the sampling tests. Default and full levels.

## Phase 2: size the pool and the rounds on principle

Start once Phase 1's speed-ups are in, so the costs you optimise are the real ones. Use WS4's test geometries when they are published (`/data/work/EMC2_workstreams/geometry/STATUS.md`). Until then, use the four fixtures above, plus T30 for strongly correlated group means.

### 6. What M and R buy

The move is an ensemble (Neal) / data-augmentation update: the group state moves as if each subject could be any of its M pool members.

- **Pool size M:**
  - costs M - 1 likelihood rows per subject per iteration (on the workers, riding on the particle round);
  - controls how far the group can move;
  - its value depends on the importance efficiency of the pool proposal `q_s` against each subject's posterior. That efficiency falls with p and with mismatch between `q_s` and the posterior, so a fixed M buys less on larger or worse-matched models.
- **Rounds R:**
  - cost hierarchy arithmetic only (cheap after Phase 1);
  - each round is a Gibbs sweep on the pool-augmented target, so the rounds saturate once the group has explored what the fixed pools allow. Extra rounds then buy nothing.

**Measure, on saved stationary states, without fits:**
- the group's movement per iteration (squared jump of the group means and log sds, standardised) as a function of M and R;
- how often subjects' selections change;
- the pool importance efficiency (weight ESS / M) per subject;
- the CPU, split into likelihood rows and arithmetic.

Sweep M over roughly 8-128 and R over 1-20.

### 7. Derive the rule

From step 6 and the ensemble-MCMC literature (Neal 2011; Shestopaloff and Neal 2013), and from the particle-step analysis WS3 is doing (the same importance-efficiency argument applies), derive:
- **M**, from measured pool efficiency and per-row cost. For example: the M at which the expected group movement per unit cost peaks, given each subject's efficiency.
- **R**, from measured saturation, for example the round after which an extra round's expected gain falls below its cost.

The rule must:
- adapt from quantities the sampler already measures, or can measure cheaply at rebuilds;
- change M and R only between blocks, never within one, so the kernel stays exact;
- contain no constant chosen to suit one dataset.

If a fixed default is genuinely as good across the test set, say so with the evidence and keep it simple. That is an acceptable outcome. Coordinate with WS3 through your reports: the pool rows share the particle round with the particle step's own rows, so both counts set that round's cost.

### 8. Check the rule

Check it on short continuations (about 400 iterations from saved states, one seed first, a repeat only where close). Measure:
- group lag-k decorrelation against the checkpoint variance;
- worst-coordinate ESS per CPU-second and per wall-second;
- the chosen M and R per fixture;
- the allocator's plans in both arms (they will shift, because the ensemble is part of "the rest of the iteration" the allocator measures).

The existing ensemble invariance tests must pass for any M and R the rule can choose. Add a test that the move is invariant when M and R differ between blocks.

## Deliverables

1. Phase 1: commits on `ws5/<topic>`, one per part, each with identity evidence and its measured saving.
2. `ws5_report.md`, containing:
   - the cost model (how the move's cost scales with p, subjects, pool and rounds), before and after;
   - per-fixture timings;
   - what remains, and why.
3. The new profile fields for ensemble and collapsed time.
4. Phase 2: the sweep results (step 6), the derivation (step 7), the rule as its own commit with tests, and the continuation results (step 8), in `ws5_report.md`. Single seeds are labelled indicative.

## Boundaries

- You own `R/ensemble_move.R`, the ensemble's call sites in `R/sampling.R`, the Gibbs group-update functions, any new C++ for them, and the iteration profile fields.
- Do not change the particle step (WS3) or the allocator (WS1).
- Do not change between-block code (WS2) except the profile schema, and tell WS2 in your report if you touch it.
- The move always runs every iteration. Phase 1 changes no behaviour. Phase 2 changes only M and R, by the derived rule, between blocks.

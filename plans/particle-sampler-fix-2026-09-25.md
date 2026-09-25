# Particle sampler bias: summary and fix

2026-09-25. Commits on `playground`: `1fd2e30f`, `3447e715`, and the
`n_blocks`/`par_groups` fix that follows them.

## Symptom

Two fits of mathematically equivalent models (RDMSWTN, and RDMSWTN_UT with
u = 0) had the same data, priors and likelihood. The likelihoods matched
exactly at 26,880 parameter vectors. The fits agreed on group means but not on
between-subject SDs (e.g. 0.74 vs 0.89), and pWAIC differed by about 90. Each
fit looked converged on its own (sigma2 R-hat ≤ 1.02).

## Cause

The subject-level particle update in `new_particle()` (`R/sampling.R`) did not
leave the posterior invariant. Two defects did this:

1. **A proposal component centred on the current state.** Every stage has a
   random-walk component around the subject's current parameters; burn has two.
   Conditional importance sampling weights each particle by
   posterior / proposal density. For the current (reference) particle, that
   density came from a component centred on the particle itself, so it sat at
   the component's peak. Staying put was therefore always penalised. Chains
   were pushed off their state, which spread subjects out and inflated the
   group variance. How much depended on each run's adapted mix, epsilon and
   proposals, so every run settled on a different, wrong answer.
2. **A forced minimum of one particle per component.** `numbers_from_proportion()`
   used `pmax(1, rmultinom(...))`, but the weights assume particles drawn in the
   mix proportions. Rare components were over-sampled and over-weighted.

The means were mostly unaffected. The damage showed up in variances, so
pWAIC-driven model comparisons from affected fits are unreliable.

## Fix

- The random-walk components are now centred on an auxiliary point
  z ~ N(current, S) instead of the current point. The weights include the
  matching N(x | z, S) factor (the extended-target construction), so the
  proposal no longer depends on the reference particle.
- Particles are allocated by a plain multinomial draw.
- Each block runs **two conditional-IS steps**: the independent proposals
  (group, chains, eff) with ordinary weights, then the auxiliary-centred random
  walk around the first step's winner. A composition of invariant kernels is
  invariant. In a single combined step, the auxiliary factor penalised every
  independent particle that jumped far. Splitting the steps recovers
  +28% min group-variance ESS, +42% mu ESS and +19% ESS per second. The
  `marginalise` path keeps the single combined step (correct, but slower).
- `n_blocks > 1` now works. It used to crash from preburn onward because
  per-block settings were never created, and it was a parallel (Jacobi) update
  in any case. Blocks now update sequentially (proper Gibbs), and the stored
  likelihood is the final-state value per likelihood component. A
  single-parameter block also crashed the sample stage's efficient proposals
  (`drop = FALSE`). Validated against exact Gibbs at 2 and 3 blocks. Each block
  uses the full particle count, so an iteration costs more, but ESS per second
  was slightly better than with one block in the test model.
- `par_groups` labels no longer need to be 1..k. Previously `c(2, 2, 5, 5)`
  silently treated the second block as diagonal. The same applied to the prior
  sampler.

## Evidence

- **Exact kernel test (conjugate Gaussian):** the old code was off by up to
  z = 23; the new code has |z| ≤ 2.5 across 16 stage/mix configurations.
  This is now `tests/testthat/test-particle-invariance.R`
  (`EMC2_TEST_LEVEL=full`, 15 s).
- **Full `fit()` against exact Gibbs.** The likelihood was mocked as a
  conjugate Gaussian; everything else (stages, adaptation, worker pool, Gibbs
  step) ran for real.
  - Old code: group variance inflated by up to 15% (z up to 7.3), with the
    affected parameters changing from seed to seed. This happened under
    diagonal-gamma and blocked `par_groups`.
  - New code: within Monte Carlo error for diagonal-gamma, blocked and
    standard (full inverse-Wishart).
- **Ruled out:** the likelihood and truncation, stopping rules, RNG streams,
  particle-count adaptation, and the group-level Gibbs updates. The Gibbs
  conditionals were checked for every covariance type.

## What to do

- **Refit anything whose conclusions rest on group variances, pWAIC/WAIC/LOO
  or subject-level spread.** That means any fit run before `1fd2e30f`.
- ESS per iteration is lower than the old sampler reported. Part of the old
  ESS was an artefact of the same bias (chains were forced to move), so budget
  somewhat longer sample stages.
- Tuning experiments that were measured and rejected: auxiliary scale ≠ 1,
  and a lower random-walk acceptance target/epsilon floor. The synthetic-grid
  gains did not survive full adaptive fits, so validate any future tuning in
  real `fit()` runs.

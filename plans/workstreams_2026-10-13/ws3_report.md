# WS3 particle step report

## Work completed

`R/sampling.R::new_particle()` now builds each candidate cloud in its destination matrix. It no longer creates a list of proposal matrices, row-binds that list, and then makes a second full-state matrix from a repeated vector. The active parameter columns are materialized once per cloud and reused for the Gaussian density calls.

The row maximum in the mixture log density now reduces the columns with `pmax()` instead of calling `apply(..., 1, max)`. It selects the same maxima in the same component order and leaves the following log-sum-exp calculation unchanged. Neither change alters proposal draws, particle allocation, proposal densities, or adaptation.

## Draw identity and tests

I compared the pre-change function from `playground` with the changed source on the saved Forstmann state. For each of preburn, burn, and sample, I checked both a single parameter block and two blocks. The sample cases also used a mixture with zero-weight components. In all six comparisons, the full result and final `.Random.seed` were identical.

The private-library test run passed `test-particle-invariance.R`, `test-particle-adaptation.R`, and `test-sampling-rejection.R`: 22 test blocks, no failures or errors, one skip. `git diff --check` and an R parse check also passed.

## Cost profile

The following profiles repeatedly called one subject's particle step at the last saved sample state. They are mechanism screens, not fit-level ESS benchmarks.

| State | Particles | Before | After | Profile evidence |
|---|---:|---:|---:|---|
| Forstmann LNR, first subject, p = 4 | 100 | 1.812 s / 1,000 calls | 1.624 s / 1,000 calls | Likelihood was 56.5% of the old profile. `apply()` used 0.10 s; it no longer appears in the changed profile. |
| lba7r sample state, first subject, p = 20 | 224 | 46.31 ms / call | 47.06 ms / call | Likelihood was 96.0% before and 96.8% after. R allocation changes are below the timing noise at this cost. The runs used 1,000 and 200 calls, respectively. |

Each profile is one run, so the Forstmann time difference is indicative; the speedup is not established. The LBA result shows why the same R optimization matters much less when each likelihood row is expensive.

## Particle-count screen

I swept total particles N = 25, 50, 100, 200, and 400 at each saved state, with 40 transitions per N from the same fixed subject and group state. Proposal settings were reset for each transition. Weight efficiency is `weight ESS / (independence particles + reference)`. ESJD uses the group covariance to whiten the final subject displacement. CPU cost is per transition.

| Model and N | CPU ms | Relative weight ESS | Move probability | Independence move probability | ESJD / CPU ms |
|---|---:|---:|---:|---:|---:|
| Forstmann, 25 | 1.05 | 0.545 | 1.000 | 1.000 | 1.221 |
| Forstmann, 50 | 1.15 | 0.545 | 1.000 | 0.975 | 1.321 |
| Forstmann, 100 | 1.80 | 0.510 | 1.000 | 1.000 | 0.786 |
| Forstmann, 200 | 2.60 | 0.525 | 1.000 | 1.000 | 0.614 |
| Forstmann, 400 | 4.45 | 0.494 | 1.000 | 1.000 | 0.366 |
| lba7r, 25 | 7.55 | 0.0584 | 0.075 | 0.025 | 0.167 |
| lba7r, 50 | 16.05 | 0.0310 | 0.125 | 0.050 | 0.116 |
| lba7r, 100 | 21.58 | 0.0170 | 0.325 | 0.100 | 0.213 |
| lba7r, 200 | 40.58 | 0.0099 | 0.400 | 0.125 | 0.131 |
| lba7r, 400 | 79.35 | 0.0062 | 0.450 | 0.250 | 0.097 |

The least-squares cost fits over these five counts were about `0.79 ms + 0.0092 ms × N` for Forstmann and `4.00 ms + 0.187 ms × N` for lba7r. Forstmann's highest measured ESJD per CPU was at N = 50; lba7r's was at N = 100. These are single-state, single-seed screens with 40 transitions per count, not stable optima. For lba7r, relative weight efficiency falls sharply with N, which violates the constant-efficiency approximation that would make a simple closed-form optimum predictive.

For the plan's approximate independence-step model, let the probability of selecting away from the reference be

`P_move(N) = e (N - 1) / (1 + e (N - 1))`,

where e is relative importance efficiency. If the mean squared jump conditional on moving is J and the step cost is `c0 + c1 N`, maximizing `J P_move(N) / (c0 + c1 N)` over continuous N gives

`N* = 1 + sqrt((1 + c0/c1) / e)`.

This derivation assumes e and J stay constant as N changes, and treats one independence step in isolation. The screen already shows those assumptions are weak for lba7r, and the production kernel composes an independence step with random-walk steps. The result is therefore a hypothesis to check, not a replacement rule.

My reading of the cited i-cSMC results is that they establish convergence conditions and particle-count bounds for conditional SMC and particle Gibbs, but do not directly optimize this sampler's ESJD per CPU. The multiple-try literature similarly supplies valid multi-candidate kernels, not this workstream's proposed move-probability approximation. See [Andrieu, Lee & Vihola (2018)](https://doi.org/10.3150/15-BEJ785) and [Liu, Liang & Wong (2000)](https://doi.org/10.1080/01621459.2000.10473908). The approximation above is the workstream's empirical model, with the optimum derived here under its stated assumptions.

## Count rule and remaining scope

I left `particle_factor`, `target_ESS`, and the adaptation rule unchanged. Two subject-state screens are not enough to set a general default, and their apparent best counts differ. The available N-back and T30 resources are data/design fixtures rather than stationary fit checkpoints, so they were not included in the saved-state sweep. I did not run short group continuations or measure allocator rates; those checks are needed before changing the production count rule.

The changes in this report are limited to `new_particle()`. The particle-count derivation still needs repeated seeds and additional states, including WS4 geometries, before it can support an adaptation change.

# Joint group–subject ensemble move: initial implementation

2026-09-27. Design plan only; no sampler code or new fits have been run.

## Recommendation

Replace the experimental group transport with an **ensemble group–subject move**. In each selected sweep, the sampler keeps several candidate values for every subject. It proposes a new group state after summing over those subject choices, then selects one subject value from each candidate set. The group proposal uses the EAM likelihoods already computed for subject particles, so it needs no additional model likelihood evaluations. This is the main advantage over adding more transported AM or RAM candidates.

The current transport in `R/sampling.R` shifts population intercepts and rescales subject deviations while keeping population correlations fixed. AM learns covariance from group draws, although this move targets a different, conditional distribution. The RAM worktree path evaluates 32 candidates per subject during adaptation. Existing real-RDM comparisons do not show a reliable gain over the corrected sampler with group moves disabled. Those findings motivate changing the conditioning of the group update, not further tuning the same transport.

This method is **exact for a finite candidate set** when the steps below are followed. Its benefit and overhead in EMC2 still need measurement. No sampler can guarantee faster convergence for every design or make an unidentified parameter data-identified.

## Method

Let `h` denote the free group parameters, `a_s` the full parameter vector for subject `s`, `L_s(a_s)` the subject's EAM likelihood, and `g_s(a_s | h)` its hierarchy density. The target is

\[
\pi(h,a\mid y)\propto p(h)\prod_s L_s(a_s)g_s(a_s\mid h).
\]

For each subject, choose a normalized proposal density `q_s(a)` with full support. Its fitted location, scale and mixture weights are **frozen before retained sampling** and do not depend on the current `h`. A practical first proposal combines EMC2's frozen empirical subject proposal with a broad Student-t component. On a scheduled sweep, place the current `a_s` at a uniformly chosen position in a pool `A_s` of `M_s` candidates; draw the remaining candidates from `q_s`. Cache each candidate's EAM log likelihood and `log q_s`.

For fixed pools, use a Metropolis step on the collapsed group target

\[
\rho(h\mid A)\propto p(h)\prod_s
\left[\frac{1}{M_s}\sum_{j=1}^{M_s}
\frac{L_s(a_{sj})g_s(a_{sj}\mid h)}{q_s(a_{sj})}\right].
\]

Evaluate the sums with log-sum-exp. Include the group prior, any coordinate Jacobian and the forward/reverse proposal ratio. After accepting or rejecting the group proposal, select each subject's new state with probability proportional to `L_s(a_sj) g_s(a_sj | h) / q_s(a_sj)`, then discard the unused candidates. Keeping the current subject in each pool and selecting its position uniformly gives the original joint posterior as the marginal of this augmented transition. The general ensemble construction is described by [Neal](https://arxiv.org/abs/1101.0387) and its parameter-update use by [Shestopaloff and Neal](https://arxiv.org/html/1305.0320#S5).

This cannot be implemented by blindly reusing all of `new_particle()`'s current candidates: its group-centred and auxiliary-centred components have different generation laws. Start with a separate, group-independent pool step that uses part of the **existing subject-particle budget**. Keep the corrected local subject update as a separate step after selecting from the pool. Retain the ordinary group Gibbs update. Do not run a Gibbs update conditioned on the old selected subjects between the collapsed group move and subject selection.

## Initial code work

1. Add a small hierarchy interface that evaluates `p(h)` and batched `g_s(a | h)`, maps unconstrained group coordinates to native parameters, and rebuilds the accepted group state. For the first implementation, cover standard, blocked standard and diagonal-gamma hierarchies, including group-design slopes. Use a Cholesky parameterization so a proposal can change **all** free covariance correlations. Hold auxiliary hyperparameters fixed during this move and update them through their existing valid kernel.
2. Add the pool builder and collapsed group transition in `R/sampling.R` or a focused `R/ensemble_move.R`. Use one group candidate per selected sweep. Adapt its proposal during discarded warmup from this move's own acceptance feedback; freeze it for production. Start with a fixed sweep cadence. Measure hierarchy arithmetic, likelihood calls, memory and worker communication before changing the cadence.
3. Carry the candidate values, likelihoods and proposal densities through the existing subject-worker path. Cache the current group target, factor each proposed covariance once and update stored subject likelihoods after selection. Failed evaluations must retain the current joint state and expose their cause. Checkpoint only after a complete sweep, with proposal settings and RNG state saved.
4. Freeze the subject mixture, epsilon, particle count and efficient proposal at the production boundary as well. `update_pm_settings()` currently continues adjusting these during sampling; a frozen group proposal alone does not freeze the full sampler.

The core is model-independent when the observation likelihood depends on `h` only through `a_s`. Other hierarchy representations can use the same move once their native group and subject densities are supplied. The initial implementation does not claim package-wide coverage until those adapters exist. A model with a direct group effect in its observation likelihood needs that effect included in the group-target evaluator.

## Minimal tests

1. **Exact transition on real EAM likelihoods.** Generate a tiny identified RDM dataset with `make_random_effects()` and `make_data()`, fixing `s` with `constants = c(s = log(1))`. With two subjects and two or three candidates each, enumerate every subject-index combination. Compare that exact sum with the implemented product of sums, check the Cholesky Jacobian and forward/reverse Metropolis probabilities, and check the selection probabilities after both acceptance and rejection. Use the actual RDM likelihood; no Gaussian observation surrogate is needed.
2. **Small correlated RDM fit.** Generate one identified, mixed-sign correlated dataset with `make_random_effects()` and the compiled RDM simulator, using roughly 10–12 subjects and a modest trial count. Fit the corrected baseline with group moves disabled and the new move, using equal likelihood-evaluation budgets, three dispersed chains and two fitting seeds. Report group means, SDs, correlations, selected joint projections, cross-seed agreement, maximum R-hat and elapsed time. Keep the R-hat target at **1.1**. If the baseline fails to explore the posterior, do not use agreement with it as a correctness certificate; the exact-transition check and agreement among dispersed new chains remain necessary. The newly fixed `WorkingTests/correlated_group_rdm_fixture.rds` is an optional larger follow-up, not the minimum test.
3. **Easy-model overhead check.** Run one small, identified LBA fit generated the same way, with diagonal or mildly correlated random effects. Compare total time, model likelihood evaluations and the subject/group posterior summaries with the corrected baseline. The new move should not require extra EAM likelihood evaluations or materially delay this easy fit.

Check one serial and one worker run for identical candidate/subject mapping and correct saved likelihoods. Do not run a large SBC campaign for this first decision. Use ESS only to describe efficiency after examining R-hat and posterior agreement; it is not a convergence threshold.

If the finite-pool identity fails, the posterior shape differs across dispersed runs, or pool reweighting costs more than it saves on the small EAM fits, stop before extending the engine. A successful first test justifies implementing the remaining hierarchy adapters and broader validation; it does not justify changing the package default by itself.

# Proposal: joint group–subject translation move

2026-09-25. Status: proposal, nothing implemented.

## Problem

The sampler draws the group mean from its conditional given the subject parameters (`gibbs_step_diag_gamma`, `R/variant_diag_gamma.R`). Given alpha, mu has sd sigma/sqrt(n) per parameter. When subjects are weakly informed about a parameter, the marginal posterior of mu is much wider than that. mu can then only move as far as all subjects drift together, so it mixes like a centred hierarchy in a funnel.

Applied witness: the post-fix UT0 refit, `Fits/Match_Control_Exp1_rdm7s_mass_noA_t0R_UT0_vs2_mt_d42_pfix.RData` in the N-back project. It ran 3 chains of 1,500 draws, used up max_tries, and reached max mu R-hat 1.22. The ratio of the conditional sd sigma/sqrt(n) to the marginal sd of mu predicts the mu ESS ordering across all 16 parameters (log-scale correlation 0.92).

| parameter | mu ESS | marginal sd | sigma/sqrt(n) |
|---|---:|---:|---:|
| v.M | 43 | 0.193 | 0.038 |
| v.d_match_LoadLow | 56 | 0.169 | 0.047 |
| B.w_prac | 2302 | 0.031 | 0.029 |

The particle update can't fix this: it moves each subject given mu. The data pin a combination of mu and all subjects moving together, and no single-block Gibbs step makes that move.

## Move

After the particle step of each iteration (so alpha and `subj_ll` are current), and before `fill_samples()`:

1. Draw delta ~ N(0, lambda^2 V), where V is a proposal covariance for mu over the non-nuisance, non-marginalised parameters.
2. Propose mu' = mu + delta and alpha'_s = alpha_s + delta for every subject s. sigma is unchanged.
3. Accept with log ratio sum_s [ll_s(alpha'_s) - ll_s(alpha_s)] + log p(mu') - log p(mu).

Why this ratio is correct:
- p(alpha | mu, Sigma) is invariant under a common translation, and the Jacobian is 1.
- The proposal is symmetric, so no Hastings correction is needed.
- It is exactly the Metropolis update of mu in the non-centred parameterisation eta_s = alpha_s - mu. Alternating it with the existing centred Gibbs draw is an ASIS-style interweaving (Yu & Meng 2011).

Implementation notes:
- **Cost:** one likelihood evaluation per subject, run in parallel on the existing worker pool, i.e. one particle's worth. With tens to hundreds of particles per subject that is at most a few percent per iteration, so K > 1 moves per iteration are affordable.
- **V:** the across-chain covariance of mu over recent draws, built at the same point as `chains_var` in `create_chain_proposals()` (`R/fitting.R`). Before adaptation, use sigma^2 / n. A diagonal V, or per-`par_groups` blocks, is a safe first version.
- **lambda:** adapt towards about 0.25–0.3 acceptance during burn/adapt and freeze it in sample. Mix adaptation currently continues in sample; don't copy that.
- **Covariance types:** location invariance holds for diagonal, diagonal-gamma, blocked and standard. For regression group models (`group_design`), translate the intercept columns only, so that X_s beta shifts by delta for every subject.
- **Exclusions:** nuisance and marginalised parameters stay fixed, and subject-specific constants are untouched.
- **Scale move (optional second step):** set alpha'_s = mu + e^kappa (alpha_s - mu) and Sigma' = e^(2 kappa) Sigma. It needs the Jacobian e^(n p kappa) and the Sigma prior ratio. It targets slow group variances, and is worth adding only if the translation alone leaves sigma2 ESS poor.

## Validation

1. Extend `tests/testthat/test-particle-invariance.R` (conjugate Gaussian, exact Gibbs reference) with the move on. Moments must match within MCSE for every covariance type.
2. Mocked-likelihood `fit()` against exact Gibbs, as in `plans/particle-sampler-fix-2026-09-25.md`, with a deliberately weakly informed subject level (subject posterior sd >= between-subject sd) so the old kernel is slow.
3. Applied: refit the UT0 witness. Success means v.M and v.d_match_LoadLow mu ESS rise several-fold per wall-clock second (the ceiling is roughly the marginal-to-conditional sd ratio squared, up to 25x for v.M), and the fit meets the default stopping rule.

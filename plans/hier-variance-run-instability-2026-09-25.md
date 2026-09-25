# Investigation: hierarchical variance posterior differs between runs of the same model

> **Resolved 2026-09-25.** The cause was a bias in the subject-level particle update, not mixing, RNG or the likelihood. See `plans/particle-sampler-fix-2026-09-25.md` (commits `1fd2e30f`, `3447e715`, `a994a8a7`). Refit affected models before comparing group variances or pWAIC.

Date: 2026-09-25

Status: open. Saved-fit witness: `WorkingTests/repro_hier_variance_saved_fits.R`.
The synthetic search script `WorkingTests/repro_hier_variance_run_instability.R`
uses the default standard prior, whereas the applied fits use `diagonal-gamma`;
it is not yet a reproduction of the applied variance-prior case.

## Sigma2 stopping-rule test (2026-09-25)

Running EMC2's actual `check_progress()` on each saved fit with
`selection = c("alpha", "mu", "sigma2")`, `max_gd = 1.1`, and `iter = 1000`
returns `done = TRUE` for both. The within-fit maximum sigma2 R-hats are 1.020
and 1.011. Thus adding sigma2 to the sample-stage R-hat selection at the
current threshold would have stopped both fits at the same 1000 draws and
cannot by itself fix this discrepancy. The rank-normalized R-hats are smaller
(1.015 and 1.009), so switching R-hat versions alone does not flag it either.

The saved-fit witness verifies all 42 subjects' trial columns and every prior
value are equal, and tests the two current likelihood implementations at
posterior parameter vectors (maximum difference zero). The six chains pooled
across runs have maximum SD R-hat 1.249, versus within-fit maxima near 1.01.
For `v.d_match_LoadHigh`, the two posterior SD medians are 0.737 and 0.888.
The diagonal-gamma conditional variance expectations agree with the saved
variance draws within 0.8% for all 16 parameters; the discrepancy enters
through the sampled subject-parameter dispersion, not the gamma variance draw.

Short same-seed RDMSWTN fits reproduce exactly under fixed parallel settings,
with distinct chains. A self-contained four-parameter simulated RDMSWTN with
the diagonal-gamma prior agreed across three different seeds (12 subjects,
30 trials per design cell, 500 sample draws). This negative control does not
reproduce the applied failure. Continuation of the two applied fits made their
new-draw SD medians agree within 0.0081 across all 16 parameters; the initially
high-variance fit's maximum sigma2 R-hat rose from 1.011 to 1.200 while it
moved. This supports an initially misleading convergence diagnostic, although
the exact sampler mechanism remains open. See
`plans/hier-variance-handoff-2026-09-25.md` for the completed experiment and
next steps.

## Observation

In the N-back PM project, two fits of an identical model gave different hierarchical variance posteriors:

- **What was identical:** data, priors, bounds and transforms. The summed log-likelihood at the same parameter vectors matched to machine precision: RDMSWTN vs RDMSWTN_UT with u fixed at 0, -3705.431846 both, max per-subject difference 0.
- **What still differed:**
  - Group means agreed, but every between-subject SD was larger in one run (e.g. 0.74 vs 0.89 and 0.28 vs 0.35), far beyond Monte Carlo error.
  - pWAIC was 468 vs 556, and higher in all 42 subjects (about +2.2 each).
  - lpd differed by only about 4.
- **Each run looked converged:** max R-hat on sigma2 <= 1.02 and min ESS about 500. Both used the default Mersenne-Twister RNG, 3 chains in parallel, 9 cores per chain, and default `fit()` stopping rules.

Consequence: information-criterion differences driven by pWAIC (up to about 90 LOOIC) are not reproducible, so model comparison currently relies on lpd and posterior-predictive checks.

The applied fits are `Fits/Match_Control_Exp1_rdm7s_mass_noA_t0R_vs2_mt_d42.RData` and `..._UT0_...` in `/data/work/PM/NirvanaHons_Nback`. Pointwise values are in `Fits/pw_*.rds` there.

## Hypotheses to test, roughly in order

1. **Chain initialisation / RNG streams.** Do the 3 chains in one run start from the same or correlated states? If so, R-hat measures agreement between near-copies, not convergence. Relevant history: the earlier L'Ecuyer stream-reuse bug (`plans/rng-stream-overlap-2026-09-22.md`). Check whether a similar reuse survives under Mersenne-Twister with `mclapply`, e.g. forked workers inheriting the same `.Random.seed`, or the same seed being set per stage or per block.
2. **Seed handling across `fit()` stages.** Is `set.seed()` before `make_emc()`/`fit()` enough to make a serial run reproducible? If not, something reseeds or draws from an unseeded source (C++ RNG without `GetRNGstate`/`PutRNGstate`, system time, a per-worker seed).
3. **Slow mixing of the group variances.** Variance parameters mix slowly under PMwG, so the stopping rule (R-hat on mu) ends sampling while sigma2 is still drifting, and different runs stop in different places. Test: longer sampling stages or `stop_criteria` including `sigma2`, then check whether between-run agreement improves.
4. **Stage-specific differences** (preburn/burn/adapt proposal construction) that let early randomness fix which region of the variance posterior each run ends up in.

## What the reproduction reports

- Posterior median between-subject SD per run, against the simulated true sample SD.
- Pairwise between-run differences in MCSE units. |z| > 4 means the difference is not Monte Carlo noise.
- Whether same-seed runs are bit-identical, with parallel chains (A vs B) and with serial chains (S1 vs S2).
- The first retained sigma draw per chain. Identical or near-identical columns suggest shared initialisation or streams.
- lpd and pWAIC per run.

Expected with a healthy sampler:
- Serial same-seed runs are identical.
- Parallel same-seed runs with the same core settings should be bit-identical
  under the current per-chain L'Ecuyer stream implementation, even when the
  caller uses Mersenne-Twister.
- All runs agree on SD medians within a few MCSE.
- pWAIC agrees within a few units.

Any |z| well above 4 between C, D, A and E reproduces the problem.

## Done when

Either a root cause is found (with a fix and a regression test), or it is shown that the default stopping rule is too loose for sigma2, with a recommended `fit()` configuration that makes between-run SD medians and pWAIC agree within MCSE.

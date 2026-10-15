# Allocator follow-up: decisions, joint predictors, and the kappa floor

The allocation performance problem remains open. This evaluation found a correctness error in the untested kappa code, tested its correction, ruled out several simple explanations for the spending regression, and tested a joint predictor that improves PM forecasts but fails the cross-fixture check. None of the experimental scoring policies below is installed in the sampler.

**Concurrent edits:** around 07:00–07:02 UTC another session removed kappa and its tests from the working tree and added other allocation changes. Those edits were retained. The corrected kappa implementation and the version used for the offline comparisons are preserved in `evaluated_allocation.R`, extracted from the tested private package. The kappa correction is consequently not a current working-tree runtime change. The full test pass below applies to the private build and test version used during this evaluation, not to the subsequently edited allocator. The standalone `kappa_floor_regression.R` checks the preserved correction independently.

## Recommendation

Keep the point estimate and the held-HMC solve correction. Do not treat a scalar kappa multiplier, removal of subject credit alone, or a joint covariance tau estimate as a demonstrated performance fix. In particular, do not replace the current objective with the joint predictor on the strength of its PM burn average: its first three decisions and its continuation and T30 results are substantially weaker.

The next allocation experiment should change what is purchased: evaluate group movement over the burn horizon, with subject purchases constrained by their own conditional progress, and compare the resulting policy with fixed subject-HMC rates zero and one. A larger estimate of group tau is insufficient when the model still buys moves using conditional bounce. Those two end-to-end controls and three paired seeds for a successful candidate remain to be run; the current offline results do not establish a speedup.

## Correctness correction tested in the private build

`.emc_alloc_eff_inv()` now returns an unreachable target below the kappa floor. With nonnegative raw rho,

```
tau_effective = kappa * (1 + rho) / (1 - rho) >= kappa.
```

Previously, inverting an effective target below that floor clamped raw rho to zero. A shear unit with raw rho zero and kappa four was consequently priced at zero cost for an effective rho target 0.2, even though its effective rho was 0.6. The inverse now returns -1 for this target, which both existing pricing paths reject. The boundary at the floor remains reachable.

The regression test covers both scalar and batched pricing at targets below, at, and above the floor. The inverse round-trip assertion now applies to reachable targets. The stage-tau simulation test pins its random generator and uses 20,000 draws rather than 2,000: the original realization estimated tau 25.25 for an AR(1) whose population tau is 19, failing its population-value assertion despite agreement with `posterior::ess_basic()`.

This fix establishes arithmetic correctness; it does not validate kappa as a spending policy. The external removal of the experiment supersedes the runtime fix; it was not reintroduced into the subsequently edited files.

## Replays of allocation decisions

`replay.R` measures each saved state once, evaluates candidate policies against those measurements, shares the HMC choice across chains, takes the largest chosen particle rate per subject, and re-solves at held rates. No likelihood transition or new fit is run. The baseline allocator source is pinned to `bf677e14`; the point policies include the held-HMC correction.

Costs below are **predicted reserved CPU seconds per iteration**, averaged over chains. They are not realized continuation timings. Objectives from different tau models are not comparable.

| Saved state | Point cost | Without subject credit | Cost change | Kappa cost | Point shear rounds | Kappa shear rounds |
|---|---:|---:|---:|---:|---:|---:|
| pair_new_2 | 9.68 | 10.00 | +3.3% | 9.07 | 2.99 | 2.14 |
| pair_new_3 | 9.98 | 9.98 | approximately zero | 9.72 | 3.73 | 3.38 |
| pair_point_2 | 13.08 | 12.63 | -3.4% | 12.88 | 3.74 | 3.49 |
| pair_point_3 | 13.27 | 13.28 | approximately zero | 12.95 | 4.26 | 3.78 |
| arm_pointsave_1 | 15.20 | 14.46 | -4.8% | 13.63 | 6.71 | 4.78 |

The computed values, including every chain, are in `replay.csv` and [RESULTS.md](RESULTS.md). Removing credit does affect some decisions, but these replays do not support it as a sufficient explanation or fix for the previous 50% spending regression. Scaling group tau by kappa can reduce cost partly by cutting shear, rather than redirecting the budget toward it.

Two additional decompositions were evaluated in `subject_objective.R`. Optimizing each subject's own tau times its allocated share of iteration cost bought HMC at average rates approximately 0.8–1.0; it recreates the overspending tendency. Buying only the cheapest subject kernels that preserve the particle-only coordinate targets produced much smaller changes, often similar to removing group credit. That policy also removes the incentive to improve a genuinely slow subject beyond its baseline, so it is not a general replacement.

`conditional_subjects.R` is a further offline prototype: it computes subject positions in the shear's changing Gaussian conditional frames, then approximates their particle and HMC factors using the reference-frame jump evidence. At `arm_pointsave_1` the predicted cost changes from 15.20 to 14.26; at `pair_point_2`, from 13.08 to 12.44. It changes spending modestly and does not resolve group tau. The approximation does not reproduce the wall quantile maps or remeasure each jump in its instantaneous frame; it must not be used as runtime measurement code.

### Replay of the concurrently edited policy

`latest_allocation.R` captures the other session's source at approximately 07:06 UTC. That version removes group credit for subject moves by default, adds deviation coordinates to subject refresh requirements, and prefers cheaper HMC regimes when their objective is within the existing gain threshold. It was evaluated read-only using the same measurements; it was not installed into the private test library.

| Saved state | Point cost | Captured latest cost | Cost change |
|---|---:|---:|---:|
| pair_new_2 | 9.68 | 10.11 | +4.5% |
| pair_new_3 | 9.98 | 10.08 | +1.0% |
| pair_point_2 | 13.08 | 12.74 | -2.6% |
| pair_point_3 | 13.27 | 13.37 | +0.7% |
| arm_pointsave_1 | 15.20 | 14.39 | -5.3% |

The original PM checkpoint was also replayed as a sample-stage rebuild. Mean predicted cost changes from 10.65 for point to 10.47 for the captured latest policy, about -1.7%. Its shared particle rate rises from 0.611 to 0.619. This is a saved-measurement rebuild, not a recreation of initial calibration or a continuation. These results do not establish that the new policy closes the observed spending gap; realized cost and ESS require paired continuations. They also do not validate later edits made after the snapshot was taken. See `replay_latest.csv` and `replay_checkpoint.csv`.

## Joint covariance prediction

`prediction.R` fits on 100 draws per chain. It predicts the next block's lag-10 autocorrelation and compares tau with a Geyer/ESS reference from the following 300 draws where available. Training and evaluation draws do not overlap. The eight slowest coordinates are selected from the evaluation draws to describe slow-tail error, not to fit the predictor.

The joint model uses covariance C and the full one-step jump covariance S. In the reversible linear approximation, its coordinate tau is

```
tau_j = 4 * (C %*% solve(S) %*% C)[j, j] / C[j, j] - 1.
```

In one dimension this is the current ESJD/AR(1) score. Keeping the cross terms permits slow combinations of coordinates even when their separate lag-one correlations look moderate. The implementation works in covariance-whitened coordinates, discards covariance eigenvalues below 1e-8 of the largest, and clips the inferred retention eigenvalues into [0, 1). It is a projected surrogate, not an exact tau identity for the complete sampler. The motivation is related to the [variational approach to slow processes](https://arxiv.org/abs/1211.7103); nonreversible chains require additional care, as in [Huang and Mao's asymptotic variance formulas](https://arxiv.org/abs/2012.13895).

The augmented version adds the subjects' mean parameters. A quadratic variant adds 15 products of the leading five principal components. Shrinking the jump covariance and fitting at lags 5, 10, and 25 were also evaluated. Longer-lag roots inflate tau on the fast Forstmann control and do not provide a reliable correction.

| Dataset | Scalar: lag-10 error / slow tau ratio | Joint group: lag-10 error / slow tau ratio | Augmented joint: lag-10 error / slow tau ratio |
|---|---:|---:|---:|
| PM burn, fixed shear directions | 0.146 / 0.139 | 0.083 / 0.672 | 0.098 / 0.826 |
| PM new seed 2, same directions | 0.156 / 0.148 | 0.087 / 0.376 | unavailable |
| PM new seed 3, same directions | 0.217 / 0.118 | 0.123 / 0.288 | unavailable |
| PM point seed 2, same directions | 0.141 / 0.102 | 0.102 / 0.206 | unavailable |
| PM point seed 3, same directions | 0.144 / 0.234 | 0.114 / 0.474 | unavailable |
| T30 production | 0.086 / 0.246 | 0.099 / 0.428 | 0.107 / 0.459 |
| Forstmann production | 0.044 / 1.008 | 0.044 / 1.089 | 0.044 / 1.113 |

The tau ratio is predicted / measured, with one the desired value. The lag-10 column is median absolute error. PM's fixed direction basis is taken from the checkpoint, as in the original `predtest.R`; these transformed forecasts do not recreate the historical directions available at every early burn iteration. Raw-coordinate results are also saved.

The attractive PM burn average hides the timing problem. For its first three 100-draw training blocks, augmented joint predictions are only 0.49, 0.66, and 0.46 of the subsequent slow-coordinate tau in the fixed direction basis. The quadratic predictor raises these to 0.56, 0.74, and 0.58, while worsening lag-10 prediction relative to the simpler joint group predictor.

T30 burn fails in the other direction. On raw group coordinates, the joint group prediction is approximately 1.8 and 3.1 times the next block's slow-coordinate tau for the two available training blocks; the augmented predictor is approximately 2.0 and 3.4 times. Changing burn geometry makes the fitted stationary linear model a poor forecast. Its improved PM score is insufficient evidence for promotion.

## Per-family operators

`family_operator.R` estimates a positive semidefinite jump covariance for each family, selects the slow joint direction, and compares family shares with Bartlett persistent transport. It also estimates directions on one half of a block and measures transport on the other half. This comparison includes subject mean coordinates, so particle and HMC transport are not forced to zero by choosing a pure group coordinate.

For PM, median whole-block shares are approximately 70% shear, 20% group, 5% HMC, and 5% particle. Measured persistent transport on the same directions gives approximately 66%, 21%, 6%, and 8%. This is useful directional evidence, but it is not a validation of the absolute tau estimate or a changed allocation. The PM paths span two 100-row recorded windows plus a final 200-row window; the script joins and splits them by actual rows rather than assuming every saved path has 100 rows. T30 and Forstmann paths are checked against stored state increments before use.

`kernel_composition.R` also estimates each family's reversible linear action and composes them in sampler order. Family multipliers alter only that family's action. This honors more of the move structure than a scalar average of subject correlations, but its forecasts still fail: PM's median slow-tail tau ratio is about 0.37, T30 production about 0.22, and T30 burn is badly overpredicted. Some PM family eigenvalues require clipping. The ordered linear-kernel approximation is therefore also rejected as a production predictor.

## Why a scalar tail correction is insufficient

Even an exact lag-one estimate does not determine tau. For example, an autocorrelation mixture

```
r(k) = 0.9 * a^k + 0.1 * b^k
a = (0.8 - 0.1 * b) / 0.9
```

has r(1) = 0.8 for every b. With b = 0.8 its tau is 9. As b approaches one its tau grows without bound, while the first 100 lags approach a finite curve. Short-lag shape can reveal a persistent component without precisely determining its relaxation time. This is not a proof that allocation cannot work: it identifies what a successful model must add, namely reliable kernel structure or an objective tied to the observable burn horizon.

## Reproduction and validation

The scripts use the handoff's saved states under its scratch directory and the PM checkpoint path. They load a private package for offline numerical helpers and the preserved allocator definitions in `evaluated_allocation.R`; the old checkpoint's likelihood is never evaluated. No new continuation fit was launched and no comparative ESS/CPU improvement is claimed. The constraints on serial work inside chain processes do not apply to these standalone offline experiments; no new runtime predictor was added.

From the repository root:

```bash
mkdir -p .Rlib_alloc_review
MAKEFLAGS=-j24 nice -n 19 R CMD INSTALL --no-test-load --library=.Rlib_alloc_review .

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/prediction.R")' "$PWD"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/kernel_composition.R")' "$PWD"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/replay.R")' "$PWD" \
  pair_point_2 pair_new_2 pair_new_3 pair_point_3 arm_pointsave_1
VARIANTS=own_subject,baseline_subject REPLAY_CSV=replay_subject.csv Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/replay.R")' "$PWD" \
  pair_new_2 pair_point_2 pair_point_3 arm_pointsave_1
VARIANTS=conditional REPLAY_CSV=replay_conditional.csv Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/replay.R")' "$PWD" \
  pair_point_2 pair_new_2 arm_pointsave_1
VARIANTS=latest REPLAY_CSV=replay_latest.csv Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/replay.R")' "$PWD" \
  pair_point_2 pair_new_2 pair_new_3 pair_point_3 arm_pointsave_1
VARIANTS=baseline,point,point_no_credit,latest REPLAY_CSV=replay_checkpoint.csv Rscript -e \
  'source("plans/alloc_prediction_2026-10-15/eval_operator/replay.R")' "$PWD" checkpoint
python3 plans/alloc_prediction_2026-10-15/eval_operator/summarize.py
Rscript plans/alloc_prediction_2026-10-15/eval_operator/kappa_floor_regression.R

NOT_CRAN=true EMC2_TEST_LEVEL=full R_LIBS="$PWD/.Rlib_alloc_review" \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript -e \
  'library(EMC2); testthat::test_dir("tests/testthat", filter="^(move-allocation|collapsed-move|shear-move|subject-hmc)$", package="EMC2", reporter="summary", stop_on_failure=TRUE)'
```

The initial test run identified the two kappa failures described above. After correction, all four requested files passed with `NOT_CRAN=true EMC2_TEST_LEVEL=full`: move allocation, collapsed moves, shear moves, and subject HMC (`tests_final.log`, exit status zero). This run used the private corrected build before the concurrent removal of kappa. The full test command above tests whichever checkout is currently present; it does not recreate the previous test sources. `kappa_floor_regression.R` explicitly tests the preserved snapshot. `git diff --check` passed and all evaluation R scripts parsed. The package was installed successfully into the private library. No changes were committed.

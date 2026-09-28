# Release-gate status

The holdout below was rerun with the package's actual proposal counts: four
for legacy and one for AM. Results report the worst R-hat and mean across seeds
of each run's minimum directional bulk ESS. ESS per 100k likelihood evaluations
is averaged over the three seed-level minimum ESS values.

Earlier exploratory tables used stale proposal counts and are superseded by
the corrected results below.

The retained package methods are pooled legacy and regularized AM. The
analytic holdout's `disabled` arm is the original exact Gibbs sampler with no
added group proposal. The oracle is a benchmark comparator, not a package
method.

## Locked Gaussian holdout

The holdout used 500 preburn, 6,000 adaptation, 500 settling, and 4,000
production transitions per chain, with four chains and three seeds. Legacy used
1,408,000 subject likelihood evaluations per run; AM and oracle each used
352,000. The disabled control made no group-move likelihood evaluations.

| Dimension | Method | Worst R-hat | Mean minimum bulk ESS | Mean min ESS/sec | ESS / 100k evaluations |
|---:|---|---:|---:|---:|---:|
| 2 | Original sampler | 1.148 | 23 | 15.3 | — |
| 2 | Legacy K=4 | 1.001 | 3,051 | 81.7 | 217 |
| 2 | AM K=1 | 1.005 | 1,745 | 62.2 | 496 |
| 2 | Oracle | 1.003 | 2,126 | 77.6 | 604 |
| 8 | Original sampler | 2.579 | 5 | 3.1 | — |
| 8 | Legacy K=4 | 2.244 | 9 | 0.22 | 0.64 |
| 8 | AM K=1 | 1.022 | 466 | 15.3 | 132 |
| 8 | Oracle | 1.020 | 574 | 19.1 | 163 |

Against the original sampler, AM increases the mean minimum ESS about 76-fold
in 2D and 93-fold in 8D; minimum ESS per second rises about 4-fold and 5-fold.
The original arm performs the exact Gaussian group-mean and subject Gibbs
updates each iteration, then skips only the added group proposal. Against
legacy K=4, AM uses one quarter the likelihood evaluations. Legacy has higher
fixed-transition ESS and ESS/sec in 2D, while AM yields about 2.3 times the
minimum ESS per likelihood evaluation. In 8D, AM recovers the difficult
directions where legacy remains stuck and reaches about 81% of oracle minimum
ESS. The analytic runner also resamples subject effects each transition, so
residuals move while the group transport preserves them. Detailed directional
results and draws are in `results/holdout/`.

## Non-Gaussian stress suite

The corrected three-seed stress suite used 500 preburn, 6,000 adaptation,
500 settling, and 4,000 production transitions per chain. Legacy used 704,000
subject likelihood evaluations per run and AM used 176,000. Each row reports
the worst R-hat over all coordinates and seeds, the mean across seeds of
minimum bulk ESS, and mean minimum ESS per 100k likelihood evaluations.

| Target | Method | Worst R-hat | Mean minimum bulk ESS | ESS / 100k evaluations |
|---|---|---:|---:|---:|
| Rotated Student-t, 8D | Legacy | 1.055 | 93 | 13.3 |
| Rotated Student-t, 8D | AM | 1.019 | 414 | 235 |
| Neal funnel, 8D | Legacy | 2.567 | 5.5 | 0.78 |
| Neal funnel, 8D | AM | 1.403 | 15.8 | 8.99 |
| Separated two-mode, 2D | Legacy | 1.735 | 6.1 | 0.87 |
| Separated two-mode, 2D | AM | 1.737 | 6.1 | 3.47 |

AM improved Student-t and funnel exploration at one quarter the likelihood
evaluation count, though the funnel remained poorly mixed. In this isolated
translation harness, `disabled` is a fixed-state control rather than the
original sampler. Neither proposal crossed modes in the two-mode target; all
four chains retained their starting mode. Per-seed metrics and draws are in
`results/stress/`.

## Real-model probe

`real-models.R` ran two chains for 1,000 retained draws on 50 trials from each
of two subjects, with 250 AM shape/scale updates and 100 settling updates.
Both methods included all five scale coordinates. Legacy used four group
candidates per transition; AM used one. The minimum ESS values are very low
and R-hat remains high, so this is an integration/performance probe rather
than convergence evidence.

| Model | Method | Worst R-hat | Minimum bulk ESS | Minimum ESS / second |
|---|---|---:|---:|---:|
| RDM | Legacy | 1.605 | 1.67 | 0.087 |
| RDM | AM | 1.621 | 1.65 | 0.107 |
| LBA | Legacy | 1.209 | 3.38 | 0.206 |
| LBA | AM | 1.308 | 2.52 | 0.165 |

This one-seed probe does not show a consistent AM advantage: AM improved
RDM ESS/sec but reduced LBA ESS/sec by about 20%. The per-run metrics are in
`results/real-models/metrics.csv`.

## Gate status

| Gate | Status |
|---|---|
| Prior signature survives reconstructed design closures | Focused regression passes, including serialization normalization |
| AM excludes preburn and burn moments; applies scale feedback immediately | Focused regression passes |
| Failed candidate likelihood is recorded and rejected | Focused and worker-retry regressions pass |
| Gaussian holdout against disabled, legacy, and oracle | Complete at 2D and 8D; results above |
| Hierarchical subject residuals move jointly | Exercised by the analytic Gibbs hierarchy holdout |
| Scale-coordinate posterior invariance | Full targeted unknown-variance regression passes for standard and diagonal-gamma |
| Real-model AM integration | Short two-subject RDM and LBA fits reached production with scale coordinates enabled; tiny burn-in emitted max-tries warnings, so these are mechanical checks only |
| Installed spawned-worker parity | AM and disabled runs produced identical draws and RNG states at one and two workers |
| Non-Gaussian stress suite | Complete at 2D and 8D; AM remains weak on the funnel and neither method crosses the separated modes |
| 32D/128D stress suite | Not run |
| Short real-model RDM/LBA efficiency probe | Complete; chains were far from converged and AM did not consistently outperform legacy |
| Longer real-model posterior-accuracy comparison | Not run |
| Package-wide regression | Still failing. The latest `devtools::test(reporter = "summary")` run reached testthat's 10-failure reporting cap and noted 47 additional failures. Reported failures include sampling snapshots in `test-fit.R`, `test-group.R`, `test-joint.R`, and `test-model_functions.R`, a `par_group` type snapshot in `test-make_emc.R`, and an omission-mass assertion in `test-likelihoods.R`. The focused group-move engine and particle-invariance tests pass, with development-only skips. Reconcile the package-wide failures before changing the default. |

Legacy remains the default. The Gaussian holdout supports AM over legacy for
correlated 8D targets, but full package regression and longer real-model
efficiency comparisons are still required before changing the default.

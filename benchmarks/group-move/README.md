> ## USER OVERRIDE (2026-09-27) -- supersedes everything below
>
> **Stop all group-move benchmarking. Do not change sampler defaults.**
> Verdict from the N-back witness and the analytic holdout: the rigid
> interweaving move cannot deliver a material gain at realistic dimension,
> whatever the adaptation scheme. Even the ideal proposal (exact curvature,
> optimal 2.38/sqrt(d) scaling) moves the ridge parameters no further per
> proposal than the Gibbs step already does (v.M: 0.23 vs 0.20 posterior SD,
> at ~23% acceptance). A controlled on/off continuation found no gain on the
> ridge. The group move is now OFF by default (`emc2.group_move = FALSE`),
> K defaults to 1, and RAM (Vihola 2012, batched) is available as
> `emc2.group_move_method = "ram"` for experiments only. Scripts:
> `mvt.R` (synthetic RDM), `applied.R` (N-back witness). Await the user.

# Group-move release gates

This harness compares the original sampler with no added group proposal, the
existing pooled legacy move, and regularized adaptive Metropolis (AM). Legacy
remains the package default until the gates below pass. Candidate algorithms
without release evidence have been removed from the code, tests, and runner.

## Analytic hierarchy

`run.R` composes exact Gaussian subject and group Gibbs updates with the real
group-move transition. Its `disabled` arm is the original Gibbs-only sampler;
it compares that baseline with legacy, AM, and an oracle
translation proposal using the exact conditional covariance
`(P0^-1 + sum(Omega_s^-1))^-1`, scaled by `2.38 / sqrt(d)`. The oracle is a
benchmark reference, not a package method. The marginal posterior is known
exactly, so standardized posterior eigen-directions test both exploration and
uncertainty calibration. Subject effects are resampled each iteration, so the
hierarchical residuals move while the group transport keeps them fixed.
Development targets use two dimensions; the locked holdout includes an 8D
rotated covariance. Each method uses its package default proposal count (four
for legacy and one for AM); reports include likelihood evaluations and wall
time as well as ESS. These runs assess translation geometry, not the particle
sampler.

Run the short mechanical check or the analytic suites from the repository root:

```sh
Rscript benchmarks/group-move/run.R smoke
Rscript benchmarks/group-move/run.R development
Rscript benchmarks/group-move/run.R holdout
Rscript benchmarks/group-move/real-models.R
```

## Target-family and dimension checks

`stress.R smoke` checks the harness. `stress.R` compares disabled, legacy, and
AM on a rotated Student-t target, a funnel, and a separated two-mode target.
In this separate synthetic stress suite, subject residuals remain fixed at
zero, so it isolates group translation. Its disabled arm is a no-move control,
not the original sampler. `stress.R dimensions` compares those methods with an
exact-covariance oracle on 32D and 64D rotated Gaussians. The oracle is a
measurement reference, not a package method. Results from this suite do not
establish end-to-end fit performance.

## Remaining release gates

The locked Gaussian comparison and full non-Gaussian stress suite are complete;
see [RESULTS.md](RESULTS.md). Targeted tests and real-model runs exercise scale
coordinates, and the analytic hierarchy moves subject residuals. The
one-seed `real-models.R` probe compares the original sampler, legacy, and AM on
short RDM/LBA fits; the chains do not converge, so it measures integration and
efficiency only. The 32D/64D stress suite and a clean package-wide regression
remain before a default change.

Each runner fingerprints its source and settings before reusing result files.
Do not reuse results from before the AM lifecycle and scale-feedback fixes;
those caches were removed when the release suite was reset.

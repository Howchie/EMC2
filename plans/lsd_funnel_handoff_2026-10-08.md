# Handoff: slow group log-sd coordinates (the "log-sd funnel")

Branch `playground`, worktree `/data/work/EMC2_dev_oo`. Starting point: commit
`e43fdc30` (shear rounds, subject HMC, shared move allocation). Background and
results: `benchmarks/subject-hmc/README.md`.

## Problem

On hierarchical fits with a full group covariance, some group log-sd
coordinates mix slowly, and the subject random effects tied to them do too.
That is what keeps the `fit()` stop rule (mu and alpha R-hat < 1.1) open on T30
(10-19 sample tries, where 10 is the minimum) and on lbaA (2 of 4 `full` runs
hit 20+ tries). mu is fine.

T30 seed 1 (`benchmarks/subject-hmc/results_funnel/T30_seed1_emc.rds`, 3 chains
x 1800 production draws):

| coordinate | log sd (typical) | marginal sd of log sd | lag-1 rho | ESS bulk |
|---|---|---|---|---|
| B_Eneutral | -3.0 to -3.3 | 0.94 | 0.93 | 24 |
| v_Eneutral | -1.5 to -2.1 | 0.36 | 0.85 | 93 |
| t0 | -2.2 to -2.4 | 0.19 | 0.56 | 288 |

The slowest coordinate is a near-zero group sd (sd ~0.04), i.e. the classic
hierarchical funnel: when the group sd is small, every subject's deviation on
that parameter must be small, and vice versa.

## What is known

- Shear log-sd rounds (`R/shear_move.R`) sample one coordinate along a
  one-dimensional orbit with subjects transported non-centred. The orbit's
  width is 40-60% of the marginal sd for B_Eneutral (15-35% of variance),
  narrower than even a linear regression on all other group coordinates
  gives (sd ratio 0.70; `funnel/lsd_condsd.R`), so the dependence is
  non-linear. Residual-sd orbits and orbits along the posterior's own
  principal axes of whitened log V were no wider (`funnel/orbit_dirs.R`).
- Rounds help, then stop: mean-only rounds leave min log-sd ESS ~9-21, a few
  log-sd rounds per iteration give ~3x, and up to ~16 rounds per iteration
  gave no further gain at +60% CPU per iteration.
- A naive joint non-centred Gaussian proposal for (mu, chol V) from the
  subjects' Laplace fits failed: in prior-dominated directions most
  proposals had negative Cholesky diagonals (`funnel/nc_group.R`).
- Varying how many rounds the allocator buys moves T30's tries but not
  total CPU (5 seeds per arm), so tuning the allocator is not the lever.

## Constraints

- One default sampler; temporary on/off switches while testing are fine,
  an options framework is not. Don't discard a component that is merely
  tuned badly.
- No warnings as a fix for poor mixing.
- Judge by R-hat < 1.1, then ESS per CPU second. CPU from the cgroup
  (`run_arm.sh` appends `cpu.stat`), not elapsed time.
- Every move must leave the posterior exactly invariant; adaptation only at
  proposal rebuilds (AirMCMC schedule). Must work across group types
  (standard, diagonal-gamma, blocked, factor/SEM) or degrade cleanly.
- Comments: what/why only, no benchmark numbers, dates or history.
- Another session may be editing `R/utils.R` and `src/` in the same worktree.
  Build test libraries from a clean tree (`git archive HEAD` plus your
  changed files, R CMD INSTALL into a temp lib) and never touch their files.
- Don't benchmark while other fits run. Don't commit unless asked.

## Assets

- Fixtures: `benchmarks/subject-hmc/fixtures/fixture_T30.rds`,
  `fixture_lbaA.rds`; forstmann (package data) and corr15 as no-regression
  controls (see the README manifest).
- Saved fits (git-ignored): `benchmarks/subject-hmc/results_funnel/`
  (`T30_seed1_emc.rds`, `T30_seed1_panel5_emc.rds`, `lbaA_seed2_emc.rds`).
- Diagnostics: `benchmarks/subject-hmc/funnel/` (`worst.R` slowest
  coordinates with block means, `lsd_condsd.R`, `orbit_dirs.R`,
  `nc_group.R`). Scripts that call internals read `EMC_LIB`.
- Driver: `benchmarks/subject-hmc/run_arm.sh <lib> <case> full <seed> <dir>`.
- Reusable internals: `.emc_shear_fits()` (subject likelihood fits),
  `.emc_shear_v_spec()` / `.emc_shear_log_prior_v()` (the V prior), the
  orbit transport in `R/shear_move.R`, and the allocation solver
  (`R/move_allocation.R`), into which a new move can enter as a unit with
  its own measured factor and cost.

## Success criterion

T30 and lbaA: min log-sd and alpha ESS per CPU second up materially and stop
rule met in fewer tries across at least 3 seeds, with forstmann and corr15
no worse than `e43fdc30`.

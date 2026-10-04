# Mode-start benchmarks (2026-09-29)

Scripts behind the burn mode-start work in `R/mode_start.R`: the batched Hessian
and the Laplace importance-sampling diagnostic. The optimiser comparisons (BFGS
vs nlminb, restarted BFGS, Pathfinder) ran through an internal switch that has
since been removed; their results are recorded below. Run from the
repo root. `EMC_LIB` selects the EMC2 library (default: first `.libPaths()`);
outputs go to `MODE_START_OUT` (default `benchmarks/mode-start/results/`, git-ignored).

| Script | What it measures |
|---|---|
| `setup.R` | 150-iteration preburn for a case, saved as `preburn_<case>.rds`. Cases: `lba71` (forstmann LBA + dummy factor, p = 71), `joint` (forstmann LBA + LNR, p = 74), the RDM fixtures `corr`, `freeA`, `wide` from `benchmarks/group-move/fixture_fits.R` (generate those first), and the simulated `corr15`, `freeA15`, `wide15`, `lba71sim`, `jointsim`. |
| `laplace_k.R` | Pareto k and share outside the support of the Laplace importance weights, from the real burn mode start, across all cases. |
| `laplace_k_by_n.R` | Pareto k at 64 to 4096 importance draws. |
| `hessian_timing.R` | `stats::optimHess` vs the batched gradient-difference Hessian, p = 23 to 101. |
| `burn_e2e.R` | Mode start plus 100 burn iterations on the p = 71 hierarchical fit. |
| `simulate_fixtures.R` | The simulated fixtures corr15, freeA15, wide15, lba71sim and jointsim. |
| `burn_arms.R` | Burn from a saved preburn with or without the burn efficient proposals; R-hat, stuck subjects and step sizes. `STOP=fixed` runs 500 iterations. |

Findings (see the comments in `R/mode_start.R`):

- The batched Hessian is 9-14x faster than `optimHess` (2p^2 + 2p rows instead of 4p^2 serial calls).
- BFGS stays. `nlminb` only won on the smooth LBA, and lost by 10-230 nats on the RDM and joint cases near the t0 edge.
  Stationarity by case, BFGS vs nlminb: corr 97/100%, freeA 90/93%, wide 20/33%, lba71 58/100%, joint 0/11%.
- Laplace Pareto k > 0.7 for about 15% of subject-chains at p = 13, and for about 84% at p = 71-74. No draws fell outside the support, and k does not fall with more draws (it is stable from 256 up).

Example:

    CASE=wide Rscript benchmarks/mode-start/setup.R
    Rscript benchmarks/mode-start/laplace_k.R

## Findings and open items (2026-09-29)

### Pathfinder: tried and removed

Codex's `R/pathfinder_mode.R` used Pathfinder only to find the mode, which
reduces it to a multi-start L-BFGS. On wide it stopped short on 19 of 30 subjects
(up to 561 nats) and produced a covariance for 20% of them. On lba71 it reached
the mode but with 7.5x the likelihood rows.

Reworked to use what Pathfinder actually outputs (the best-ELBO Gaussian along
the path, as the start centre and covariance), it fit the conditional posterior
worse than the Laplace Gaussian on every fixture:

| case | Laplace median k | Pathfinder median k (history 5 / 20 / p) |
|---|---|---|
| corr | 0.46 | 1.06 / 0.81 / 0.88 |
| wide | 0.71 | 1.78 / 1.44 / 1.42 |

It was better on only 2-4 subjects per fixture, at 3-4x the cost. Pathfinder's
advantage is avoiding the Hessian, which pays off only at p in the thousands.
At p <= 100 the batched Hessian is exact and cheap. Removed 2026-09-29.

### How well BFGS seeds the modes

Stationarity (scaled gradient < 0.1) is a harsh proxy. What matters for a start is
the shortfall below the best mode found by any method (BFGS, nlminb, Pathfinder),
against the typical-set scale of about p/2 nats:

| case | p | stationary | within 1 nat | 1 nat to p/2 | beyond p/2 |
|---|---|---|---|---|---|
| corr | 13 | 29/30 | 30 | 0 | 0 |
| freeA | 14 | 27/30 | 30 | 0 | 0 |
| wide | 24 | 6/30 | 27 | 1 | 2 (worst 49.5) |
| lba71 | 71 | 11/19 | 18 | 1 (26.8) | 0 |
| joint | 74 | 0/19 | 18 | 1 (9.2) | 0 |

Restarting BFGS from its end point gains at most 2e-4 nats. So its stops are real
local optima or support edges, not stalls. The few badly seeded subjects sit in
other local basins, and only a multi-start would reach them.

### Burn-stage proposals: what was kept

Three changes, all in the burn stage only:

- **Efficient proposals from burn iteration 200** (`create_eff_proposals(from_burn = TRUE)`):
  always on. On every fixture they gave lower subject R-hat than burn without
  them. Without them, lba71 and joint ended with a chain that never moved (infinite
  R-hat). The `emc2.burn_eff` switch is removed.
- **Importance-resampled mode starts: removed.** No consistent gain (better on two
  fixtures, worse on two). The Pareto k diagnostic stays; the gate below uses it.
- **Pareto k gate.** A subject whose Laplace fit is poor (k > 0.7 or over half its
  mass outside the support) takes the covariance of its own burn draws once there is
  a usable one.
- **Local step controller** (`plans/sampler_bugs.md`, item 9). The random-walk step
  adapts to the share Phi(-2.38 / sqrt(2)) = 0.046, the share at the random-walk
  scale 2.38 / sqrt(d), from each iteration's share, with a gain matched to that
  target. Each state-centred component runs as its own importance step, and restarts
  at 2.38 / sqrt(d) when its covariance changes.

The old controller targeted 0.2 (about 1.19 / sqrt(d)), and its stage floor of 0.6
hid that at d <= 15 but froze every chain at d ~ 70. A step fixed at
min(0.6, 2.38 / sqrt(d)) fixed the freezing. Setting the old controller's target to
0.046 mixed worse than that fixed step. The target was not the cause:

- Two random walks shared one importance step, so their auxiliary densities
  multiplied and the tighter one pinned every particle. Both steps collapsed to the
  floor in the first 200 burn iterations.
- The step slot changed component at the efficient-proposal switch and inherited the
  floored value.
- A running-mean share and an 8x gain made it overshoot to 4-8 / sqrt(d) for ~150
  iterations.

With those fixed, the controller against the fixed step (pairs R-hat > 1.1, seeds
1 / 2 / 3; median step x sqrt(d) in brackets):

| fixture (d) | fixed min(0.6, 2.38 / sqrt(d)) | adaptive |
|---|---|---|
| freeA15 (14) | 13 / 46 / 62 (2.24) | 38 / 57 / 56 (1.8-2.1) |
| wide15 (24) | 1 / 4 / 5 (2.38) | 0 / 2 / 4 (2.28) |
| lba71sim (71) | 483 / 595 / 608 (2.38) | 378 / 392 / 379 (2.83) |
| jointsim (41 + 33) | 314 / 249 / 300 (2.38) | 177 / 180 / 219 (2.2-2.4) |

freeA15's group level is unconverged in both arms (mu R-hat 1.5-2.7), and it varies
about sixfold between seeds. On frozen-group subject posteriors (wide15), efficiency
is flat for step x sqrt(d) from 1.2 to 3.0, and the controller lands inside that range.

Earlier arms, before the controller fix: the 0.6 floor against the dimension floor.
Fixed 500 burn iterations from the same preburn, R-hat over the second half,
"stuck" = a chain moved on under 30% of the last 100 iterations
(`burn_arms.R`, `STOP=fixed`):

| fixture (d) | seed | pairs R-hat > 1.1: 0.6 floor / dimension floor | stuck subjects | max subject R-hat |
|---|---|---|---|---|
| corr15 (13) | 1 | 0 / 0 | 0 / 0 | 1.06 / 1.06 |
| corr15 (13) | 2 | 0 / 0 | 0 / 0 | 1.04 / 1.04 |
| freeA15 (14) | 1 | 73 / 31 | 0 / 0 | 2.56 / 1.83 |
| freeA15 (14) | 2 | 12 / 23 | 0 / 0 | 1.93 / 1.99 |
| wide15 (24) | 1 | 3 / 0 | 2 / 0 | 1.25 / 1.21 |
| wide15 (24) | 2 | 1 / 3 | 1 / 0 | 1.27 / 1.25 |
| lba71sim (71) | 1 | 1101 / 497 | 20 / 0 | 7.83 / 2.53 |
| lba71sim (71) | 2 | 1126 / 649 | 20 / 0 | 6.40 / 2.29 |
| jointsim (41 + 33) | 1 | 676 / 286 | 10 / 0 | 4.71 / 2.32 |
| jointsim (41 + 33) | 2 | 733 / 424 | 10 / 0 | 4.55 / 1.91 |

Seed 1, the two components separately (pairs R-hat > 1.1 / stuck):

| fixture | 0.6 floor | 0.1 floor only | gate + 0.1 floor | gate + dimension floor |
|---|---|---|---|---|
| freeA15 | 73 / 0 | 58 / 0 | 32 / 0 | 31 / 0 |
| wide15 | 3 / 2 | 42 / 0 | 61 / 0 | 0 / 0 |
| lba71sim | 1101 / 20 | 963 / 0 | 623 / 0 | 497 / 0 |
| jointsim | 676 / 10 | 895 / 0 | 583 / 0 | 286 / 0 |

freeA15 varies about sixfold between seeds at the same build, so its single-seed
differences are noise. lba71sim and jointsim are not converged at 500 iterations
under either build, but every chain moves.

Measurement traps:

- The chains draw from random streams stored in the preburn object, so
  `set.seed()` alone does not change a replicate. `burn_arms.R` clears them for
  seeds other than 1.
- The default burn stopping rule averages R-hat over all subject parameters and
  group means. It stopped freeA15 at 200-215 iterations with group-mean R-hat
  still 2.5, so burn length under that rule is noisy and was not used to compare
  builds.

### Cross-chain mode sharing at the burn mode start

Each chain's BFGS converges to the mode it is already in. A chain caught in a
local mode during preburn widens its own group variance around it, so under its
own group state the trap still wins. At the mode start, the chains' modes for
each subject are therefore scored by their Laplace log mass averaged over all
chains' group states, and a chain re-optimises from another chain's mode when that
is better by more than 20 nats. On the old freeA fixture this moved 2 of 90
chain modes (one subject's chains had split on `v`); subject-parameter pairs with
R-hat > 2 after 400 iterations fell from 17 to 4.

### New fixtures (`simulate_fixtures.R`)

The earlier fixtures were misleading for convergence:

- corr/freeA/wide inflated group SDs (x3, and SD 0.5 on all 24 wide terms), so the
  tail subjects had thresholds near zero (natural B of 0.01-0.02 in a cell). Their
  RTs are then t0 plus noise and the posterior is a ridge bounded only by the
  group prior.
- lba71/joint fit a random 10/5-level factor to about 14 trials per cell of the
  Forstmann data (24% of cells without errors), and joint fit the same data twice.

corr15/freeA15 use group SDs x1.5 with 400 trials per cell, and wide15 uses SD 0.3
with 600, keeping each group SD at the same multiple of the subject SE as before.
lba71sim/jointsim simulate the same designs from set truths (F effects: SD 0.3
around zero), 20 subjects x 40 trials per cell. `setup.R` builds their preburns.

### Bugs listed in `plans/convergence/HANDOFF.md` (all resolved)

- **"`type = "diagonal"` gives one full covariance block": wrong as stated.** `par_group`
  is `1:p` and every Gibbs sweep is diagonal. The full covariance came from two
  other bugs, both fixed 2026-09-29:
  - The starting `theta_var` was a full inverse-Wishart draw (`get_startpoints_standard`).
    Only the init slice was affected.
  - Prior sampling read a nonexistent `sampler$par_groups` (`get_type_objects.R`).
- **`predict()` with a free parameter constant across draws** (e.g. integer
  sampling): fixed. `predict()` now keeps rows without variation across draws
  (`remove_constants = FALSE`). It gives identical output to before on models with
  design constants (LNR, RDM with `s` and `A`, LBA with `sv`), because design constants
  come from the design, not the sampled matrix.
- **Joint models:** `get_data()`, `sampled_pars()` and `predict()` overflowed the C
  stack on every joint `emc` (`get_joint_names` -> `sampled_pars` -> `get_data` ->
  `get_joint_names`). Fixed.

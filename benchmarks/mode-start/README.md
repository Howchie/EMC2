# Mode-start benchmarks (2026-09-29)

Scripts behind the burn mode-start work in `R/mode_start.R`: the batched Hessian
and the Laplace importance-sampling diagnostic. The optimiser comparisons (BFGS
vs nlminb, restarted BFGS, Pathfinder) ran through an internal switch that has
since been removed; their results are recorded below. Run from the
repo root. `EMC_LIB` selects the EMC2 library (default: first `.libPaths()`);
outputs go to `MODE_START_OUT` (default `benchmarks/mode-start/results/`, git-ignored).

| Script | What it measures |
|---|---|
| `setup.R` | 150-iteration preburn for a case, saved as `preburn_<case>.rds`. Cases: `lba71` (forstmann LBA + dummy factor, p = 71), `joint` (forstmann LBA + LNR, p = 74), and the RDM fixtures `corr`, `freeA`, `wide` from `benchmarks/group-move/fixture_fits.R` (generate those first). |
| `laplace_k.R` | Pareto k and share outside the support of the Laplace importance weights, from the real burn mode start, across all cases. |
| `laplace_k_by_n.R` | Pareto k at 64 to 4096 importance draws. |
| `hessian_timing.R` | `stats::optimHess` vs the batched gradient-difference Hessian, p = 23 to 101. |
| `burn_e2e.R` | Mode start plus 100 burn iterations on the p = 71 hierarchical fit. |

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

### Burn-stage efficient proposals (Codex, `R/fitting.R` / `R/sampling.R`)

- From 200 burn iterations, every burn block rebuilds `create_eff_proposals(from_burn = TRUE)`,
  using the last 250 warm-up draws per chain. The subject kernel then runs the
  adapt-style four-component mixture (group, RW, chains, eff). `reset_pm_settings()`
  clears `mix` each burn block, so the default `c(0.1, 0.25, 0.5, 0.15)` does apply.
- Evidence so far: paired reruns from a 200-iteration checkpoint of the N-back LBA
  (`plans/convergence/repro_burn_wide_eff_from_200.R`). With eff proposals and the 3x
  widening (the kept configuration), burn exited at mean R-hat 1.033 after 4 tries;
  without them it hit max tries at 1.117.
- Still needed: more seeds and models, and a test covering the `burn_eff_ready` path
  (none exists).
- The first rebuild's window still contains the post-mode-start transient.
- The 3x widening of the burn `ensemble_q` covariance was reverted, so it is still 3x.

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

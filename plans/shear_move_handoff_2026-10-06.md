# Shear group move: handoff (2026-10-06)

## Goal
Fix non-convergence on lba71sim (LBA, 71 pars, 20 subjects × 2400 trials) with
one default, mathematically exact sampler. Chains stay independent. The shear
move (`options(emc2.group_move = "shear")`) is the candidate to replace the
ensemble group update.

## User constraints
- One default sampler, with no stacked or competing proposals and no
  arbitrary or unvalidated heuristics.
- Change only what was asked for. An earlier rewrite changed several things at
  once and confounded the comparison.
- Comments say what and why only, with no history, dates or benchmark numbers.
- Never edit a script while a run is using it. Don't benchmark while the
  user's fits are running. Commit only when asked.
- Use `posterior::ess_bulk` and `rhat`, never summed per-chain coda ESS.
  Verify a diagnosis before stating it.

## Current design (`R/shear_move.R`, untracked file)
Target: p(mu) p(V) Π_s L_s(a_s) N(a_s | mu, V). Map for a step u = (u_mu, u_sd):
- mu → mu + u_mu and V → DVD, with D = diag(exp(u_sd)).
- Each subject maps as a_s → m_s' + R_s'⁻¹ R_s (a_s − m_s), where
  (m_s, R_sᵀR_s) is the subject's conditional N(m, Q⁻¹) under its Gaussian
  likelihood fit `lik_laplace` (Q = H_s + V⁻¹). Without a fit, H = 0.
- Jacobian: Σ_i (g_i+1) u_sd,i + Σ_s log|R_s| − log|R_s'|, where g_i is the
  size of the covariance block holding parameter i.
- V prior: diagonal-gamma uses gamma(shape, rate) on each precision. Standard
  uses an IW(v+d−1, 2v diag(1/a_half)) per blocked group and gamma(v/2,
  v/a_half) for unblocked parameters.

There are two kinds of step, each a separate Metropolis–Hastings step
(`.emc_shear_mh`):
1. **Mean steps (restored to the old working version).** u_sd = 0. u_mu runs
   along the top 5 slowest lag-10 directions of the chain's own mu window
   (`.emc_shear_basis`, old code). It has its own Robbins–Monro scale
   targeting 0.3 and takes `emc2.shear_repeats` = K steps per iteration. With
   u_sd = 0 the map equals the old a_s += (H+V⁻¹)⁻¹V⁻¹u.
2. **Scale steps (new).** One per parameter i every iteration, with u_sd =
   e_i·step_i·z. Each parameter has its own Robbins–Monro step targeting
   0.44, initialised at 1/sqrt(2n). Steps adapt only in the discarded stages.

**Re-centring of the likelihood fits (new, always on when the shear is
active).** `.emc_shear_recentre_lik` runs at each proposal rebuild
(`R/fitting.R`, `create_chain_proposals`):
- it computes lin += (mean Q⁻¹)⁻¹ · mean(a_s − m_s) over the chain's window,
  thinned to 50 draws;
- H is kept;
- the particle Laplace-conditional proposal reads the same `lik_laplace`, so
  it changes too.

## Results on lba71sim, seed 2, strict burn rule, independent chains, 27 cores
`posterior` estimators:

| arm | mu_A ESS / R-hat | logsd_A ESS | worst alpha | CPU (s) |
|---|---|---|---|---|
| ensemble | 9 / 1.25 | 33 | 8 / 1.31 | 1188 |
| old shear (mean only) K=4 | 65 / 1.04 | 142 | 17 / 1.15 | 1798 |
| v2: joint 10-direction basis, shared scale (bad redesign, reverted) K=4 | 21 / 1.11 | 132 | 22 / 1.10 | 1898 |
| **v4 (current) K=4** | 19 / 1.11 | 90 | 5 / 1.54 | 14192 |
| v4 K=1 | 10 / 1.21 | 61 | 5 / 1.57 | 13066 |

Excluding chain 1 (chains 2–3 only), K=4:

| | mu_A ESS | logsd_A ESS | median logsd ESS | median mu ESS |
|---|---|---|---|---|
| old | 46 | 94 | 449 | 100 |
| v4 | 47 | 121 | 901 | 133 |

The mean part matches the old move. The scale steps double the typical
logsd ESS, but the funnel on A is only modestly better.

## Open problems
1. **Subject 11, chain 1, is in a separate region for all of production**
   (inherited from the shared warm-up; entered during adapt).
   - Chain means of A / B_Eacc / t0_Eacc: −0.29 / −0.03 / 0.21. Chains 2–3:
     −0.7 / 0.19 / 0.04.
   - Truth: −0.84 / 0.37 / −0.14.
   - Log-likelihood is about 2–3 lower. Draws have sd(A) 0.04–0.09, while the
     Laplace conditional sd is 0.17.
   - Stored log-likelihoods are exact (checked), and the chain isn't at an RT
     boundary (min RT 0.445).
   - Unverified: is it a genuine second mode or a stuck chain? Re-centring now
     centres chain 1's fit (and so its particle proposal) on that region and
     may keep it there.
   - Suggested check: maximise subject 11's conditional posterior from each
     chain's position and compare.
2. **Cost.**
   - The scan is about 1,400 serial single-subject likelihood calls per
     iteration on each chain's master: 2.5–2.7 s per move against 0.77 s per
     iteration before.
   - The likelihood is 89% of the move's time (Rprof).
   - Fix: dispatch each scale step's 20 subject evaluations to the chain's
     worker pool (see `.emc_wpool_iter` and how `R/ensemble_move.R`
     dispatches). Not started.
3. A's funnel (corr(logsd_A, subject spread) ≈ 0.72; subject A's
   likelihood-only sd 0.23–0.30 against a true group sd of 0.20) is still the
   slowest group quantity.

## Files changed for this work (all uncommitted)
- `R/shear_move.R`: the whole move, basis, recentring, priors, Jacobian.
- `R/fitting.R` (`create_chain_proposals`): calls
  `.emc_shear_set_basis(chain, mu_window)` and then
  `.emc_shear_recentre_lik(...)` when the shear is active.
- `tests/testthat/test-shear-move.R`, which covers:
  - the basis;
  - the map's inverse and a numerical Jacobian (IW block and diagonal);
  - the V prior against an independent density;
  - exact invariance with V free, for K = 1 and 4 and correct and wrong fits,
    covering both step types;
  - recovery of a known location by re-centring.
  - All pass. Removing a Jacobian term fails both the Jacobian test and the
    invariance test.
- `NEWS.md`: the shear entry.

Other uncommitted work on this branch (`playground`) predates this and isn't
part of this problem: GD rules, the Laplace-conditional proposal, the
independent-chains fix, `emc2.sample_rebuilds`, and the benchmarks.

## How to run
- Test: `Rscript <scratch>/runtest.R tests/testthat/test-shear-move.R`.
- Scratch: `.Rtmp/claude-1000/-data-work-EMC2-dev-oo/4630b319-1e4b-441f-aec7-589268bf5557/scratchpad`
- Benchmark queue: `<scratch>/shear_v4.sh`. It installs into
  `<scratch>/Rlib_padapt`, then runs `benchmarks/production-adaptation/warmup.R`
  and `production.R` with
  `EMC_LIB=... CORES_PER_CHAIN=9 CASE=lba71sim SEED=2 INDEP=1 GROUP_MOVE=shear BURN_RULE=strict`,
  then `SHEAR_REPEATS=k ARM=doubling TAG=shearK$k SAVE_EMC=1`.
  - Copy the script to a new name and a new `OUT`. Don't edit it while it
    runs. Launch with `setsid nohup`, and check by name that it has really
    stopped after any kill.
- Results:
  - `benchmarks/production-adaptation/results_shear/{ensemble,shear,shearK,shear_v2,shear_v4}/`
    hold `prod_*.rds` (`per_var`, `stop`) and the saved `emc_*.rds`.
  - The fixture with true values is `benchmarks/mode-start/results/fixture_lba71sim.rds`.
- Comparison scripts in scratch: `cmp.R` (headline table), `chk11.R`
  (subject 11 by chain), `time_scan.R` (cost of one move), `probe_scale.R`.
- The benchmark's `prod_cpu` misses the worker pool's CPU. Compare cost on
  elapsed time for sequential runs on the same machine.

## Resolution (2026-10-07)

The shear is now the second half of one group update: the ensemble sweep
(centred) followed by the shear (non-centred), i.e. interweaving. Changes
from v4, all in `R/shear_move.R`:

- The slowest-direction basis from lagged draws is replaced by the closed-form
  lag-1 operator of the centred Gibbs step under the likelihood fits; mean and
  scale steps are each taken with that predicted autocorrelation / prior share,
  so well-informed hierarchies get few steps.
- Re-centring the likelihood fits on chain draws is removed (positive feedback:
  a stuck chain re-centred its fit on itself).
- Each direction and scale adapts its own step size.
- Delayed acceptance on the fits when every subject has a checked fit.
- The Huang-Wand `a_half` is integrated out of each step (closed form) and
  redrawn after it; holding it fixed pinned weakly informed group sds.
- Subject 11 on lba71sim was a tail-stuck chain, not a second mode. The
  remaining slow mode on free-A RDM fits is the group A-B-t0 ridge: the
  Gaussian frames transport subjects linearly along curved subject ridges.

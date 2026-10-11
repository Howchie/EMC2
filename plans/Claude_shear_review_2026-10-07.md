# Shear group update: overhead, the seed-2 failure, and the fix

Independent review, Claude, 2026-10-07. Branch `playground` at 11993e52; `1aec215d`
(plan built after the fit checks) landed during the review, and the "fixed" runs below
use an equivalent change. Nothing in the working tree was edited. Scripts, logs and
libraries are in the session scratchpad (see the last section).

How things were measured:

- ESS is `posterior::ess_bulk` and R-hat is `posterior::rhat`.
- CPU is the cgroup total (`systemd-run --user --scope`, then `cpu.stat`), unless stated otherwise.
- The README's CPU figure (`proc.time()` user+sys including children) does not see the spawned
  likelihood workers. It reports about 1/5 to 1/8 of the true cost. For example, the fixed-plan
  freeA15 production reported 779 s, but the cgroup measured 6301 s. Ratios between arms happened
  to come out similar, but the absolute ESS per CPU-second in `cmp.R` is overstated about 5x.

## Summary

1. **Most of the README's overhead and the seed-2 result come from one bug.** Every shear plan
   built during warm-up was the plan for "no usable likelihood fit". That gives λ ≈ 1 and w = 1
   for all 14 parameters, hence 28 steps per iteration. The per-direction step sizes were
   adapted against that plan and then frozen in production.
   Seed 1 (the "success") ran an older overlay that used unchecked fits. So seed 1 vs seed 2
   compared different code, not different seeds.
2. **With the plan fixed, the shear still costs more and still doesn't fix freeA15.** On corr15
   it costs 1.33x the ensemble for the same ESS. On freeA15 it costs 1.31x and still fails
   (R-hat 1.43/1.47).
3. **The freeA15 slowness is weak identification, and R-hat can't see it.** Each subject's
   likelihood is flat in log A below about −1. The ensemble's seed-2 "convergence" depends on
   where the chains start: started in the tail, its chains stay there for 1000 iterations.
4. **The README's curved transport is the wrong fix.** It widens the exact path along the slow
   direction by only 1.27–1.36x. The loss is per subject: each subject's B and t0 must
   follow its own A along its own ridge.
   A per-subject **ridge transport** widens the path 1.81x (an oracle reaches 2.3x). It passes
   an exact-draw invariance test. On freeA15 it gives mu_A ESS 26–33 and R-hat 1.09–1.13 from
   three start sets, against 4–7 and 1.41–1.88 for the ensemble, at 1.19x the ensemble's CPU.
5. **Gating the shear makes it free on easy problems.** Taking only directions with λ ≥ 0.9 and
   scales with w ≥ 0.8 brings corr15 from 1.33x to 1.06x the ensemble's CPU; the residual is
   noise in the particle step.

## 1. Overhead: where the extra cost goes

Profile setup: one chain at a time with 5 workers, from the seed-2 shear warm-up state,
300 iterations × 3 chains, current code (so the plan was still the buggy one). Worker CPU is
counted per message. Arms: shear 0.894 s/iter vs ensemble 0.507 s/iter (1.76x).

| component of the extra 0.387 s/iter | s/iter | share |
|---|---|---|
| worker likelihood rounds (13.9 rounds/iter × 30 subjects, about 20 ms of worker CPU each) | 0.273 | 71% |
| master side of those rounds (serialisation, spin-wait polling, 2.8 ms per round) | 0.039 | 10% |
| Cholesky frames for each proposed state (30 × chol/chol2inv, 1.4 ms per step) | 0.022 | 6% |
| subject map, hierarchy, screen, context, a_half | 0.007 | 2% |
| planning (0.012 s in total) | ~0 | 0% |
| difference in particle-step worker CPU (noise or state) | 0.044 | 11% |

Rounds are the cost, and two things inflate them:

- **The plan bug:** 15.4 steps per iteration here, against 4.4–8.9 under a correct plan.
- **The cheap screen rarely rejects:** it passed 91% of steps. Along the slow direction, the
  Gaussian (Laplace) path is about 4x wider than the exact one (median width 0.92 vs 0.215
  posterior sd, 12 states). So the screen accepts steps that the exact stage then rejects.

Clean comparison, same fixed-plan warm-up start, 200 iterations × 3 chains, cgroup CPU:

| case | ensemble | shear, fixed plan | shear, gated + ridge |
|---|---|---|---|
| corr15 | 0.266 s/iter | 0.353 (1.33x; 2.9 rounds/iter) | 0.282 (1.06x; 0 steps) |
| freeA15 | 0.750 | 0.985 (1.31x; 8.2 rounds/iter) | 0.889 (1.19x; 4.6 rounds/iter, ridge tables 21.7 s) |

What can be cut without changing the target distribution:

1. **The plan/check order.** Already done in 1aec215d.
2. **Gating.** Step a direction only if λ ≥ 0.9 and a scale only if w ≥ 0.8, with at most two
   directions. The choice depends only on the frozen plan, so the move stays valid. Effect:
   corr15 1.33x → 1.06x. Without the cap, one chain gated in many unadapted directions and ran
   at 33 rounds/iter with 5% acceptance.
3. **Mean steps without Cholesky.** A mean step leaves V, and so every R_s, unchanged; only the
   centres move. This gives the same transitions to 5e-13 over 200 moves. It saves about half of
   the 6% frame cost.
4. **The cheap screen.** It costs 0.03 ms per step but saves little on hard problems, because the
   Gaussian path is far too wide. Keep it, but don't count on it.
5. **Polling.** The master spin-waits for replies (10% of the shear's cost, and the same pattern
   in the particle step). A blocking wait, or a longer backoff floor, would save it. This is
   pool-level, not shear-specific.

## 2. Why seed 2 failed

- **The plan bug (confirmed numerically).** `create_chain_proposals()` set the plan before
  `add_proposals()` checked the Laplace fits. Every `.emc_laplace_refresh` and burn mode start
  stores fits without `$check`, so `.emc_shear_fits()` saw no usable fit.
  The stored seed-2 plans (freeA15 and corr15) equal the no-fit plan: |Δw| < 1e-15 and
  |Δλ| < 7e-4. The plan from the checked fits would have had λ = 0.94, 0.72, 0.57, … and 6.5–8.9
  steps per iteration, not 28.
- **Frozen step sizes made the mean steps do nothing.** The steps adapted against that plan were
  0.009–0.075 posterior sd. They were frozen in production after the first rebuild produced a
  correct plan. The mean steps were then accepted 70% of the time, moved about 1% of a posterior
  sd, and still cost 14–18 rounds per iteration.
  I reran the seed-2 production from the saved warm-up with `SAVE_EMC=1`; it is bit-identical
  (R-hat 1.318, ESS 10.66).
- **Seed 1 used different code.** Its saved fit carries no `$check` on any of its 90 fits, yet its
  stored plan is the checked-fit kind (λ from 0.99 down to 0.002). So that overlay used unchecked
  fits, and seed 1's success says nothing about the packaged code.
- **Fixing the plan doesn't rescue seed 2.** Fixed warm-up plus production: R-hat 1.43/1.47,
  min mu ESS 9. Step adaptation now works (acceptance 0.44–0.48). But chain 3 ends warm-up in
  the tail (mu_A −1.89, log sd_A −0.8) and needs about 600 iterations to come back.
- **The real geometry: A is weakly identified.** For each subject's profile log-likelihood in
  log A (other 13 parameters re-optimised):
  - From log A = −5 to −1, the median drop from the maximum is 0.05–0.33 nats.
  - Above that it falls steeply: about 3.9 nats at log A = 0 and 22.6 at 0.5.
  - The prior share of A is w_A ≈ 0.9–0.98, and (mu_A, sd_A) form a funnel.
- **The ensemble's seed-2 convergence depends on the start.** I started the default ensemble
  from tail states:
  - From −2.0: one chain stays at −2.12 to −2.19 for 1000 iterations, with log sd_A down to
    −4.2. R-hat 1.65.
  - Dispersed starts, seed 91: R-hat 1.41.
  - Dispersed starts, seed 92: one chain sits between −2.0 and −3.0 the whole time. R-hat 1.88.

  A sampler that does mix (below) puts mu_A at mean ≈ −1.2 to −1.3, sd 0.2–0.3, with under 2% of
  draws below −1.8. So the ensemble's seed-2 answer (−1.17, sd 0.18) was about right. But it got
  there because of where warm-up ended, not because it mixes. "Ensemble converged, shear
  failed" in the README table reflects start states plus the bug, not the quality of the moves.
- **Ruled out:** a second mode, step-size adaptation once the bug is fixed, and a few bad fits.
  The surrogate does not systematically point the wrong way; it is too wide (see §1).

A profile-likelihood estimate of the (mu_A, sd_A) posterior (`a_reference.R`) put mu_A at
−1.54 ± 0.45. It ignores the coupling of B and t0 through the group, and it is contradicted by
three converged start sets. Don't use it as a reference.

## 3. The fix

Exact log target along the slowest plan direction, from posterior states of the seed-2 fit.
Width is the sd, in posterior-sd units, of the normalised exp(path); a wider path allows bigger
steps.

| transport of a mean step | median width | vs linear |
|---|---|---|
| Gaussian surrogate (what the screen sees) | 0.92 | — |
| linear Laplace frames (current) | 0.215 | 1 |
| README curved transport (κ regressed on the chain's draws) | 0.264 | 1.36 |
| exact conditional-mode κ along the direction (Newton per subject) | 0.249 | 1.27 |
| frames from each subject's empirical covariance | 0.186 | 0.85 |
| **ridge transport (below)** | **0.41** | **1.81** |
| oracle: each subject's non-A coordinates re-optimised (not a valid move) | 0.511 | 2.3 |

The first four rows use 12 states; the last three and their linear baseline (0.218) use 6.

**Why a curved centre fails.** At a one-sd step the exact path loses 4–29 nats in total. That
loss is spread across most of the 30 subjects; the worst single subject loses 2–7 nats. Moving
A alone, without B/t0 compensation, costs 30–180 nats. So the compensation is what matters, and
it differs by subject: each subject's B and t0 must follow its own A along its own curved ridge.
A group-level κ shifts every subject the same way and can't represent that. Each subject's small
mismatch then adds up over 30 subjects.

**Ridge transport.** It is general: it needs only the plan and the likelihood.

- For subject s, let e_s be the unit vector along its Laplace response to the slowest plan
  direction, (H_s + V⁻¹)⁻¹ V⁻¹ d.
- At each proposal rebuild, tabulate κ_s(ξ) on a grid of ξ = e_s′a (9 points, mean ± 4 sd of
  the window). κ_s(ξ) is the part of the subject's conditional mode on the slice e_s′a = ξ that the
  Gaussian surrogate misses, projected off e_s. The mode comes from Newton steps with the Laplace
  precision as the metric and central finite-difference gradients, about 1,100 likelihood
  evaluations per subject.
- Let T(a) = a − κ_s(e_s′a). The move is T⁻¹ ∘ L_u ∘ T, where L_u is the existing linear frame
  map.
- The cheap screen uses the Gaussian surrogate on T(a).

**Why it is valid.** Because e_s′κ_s = 0, T is a shear with unit determinant, and
T⁻¹(b) = b + κ_s(e_s′b). Since L₋ᵤ = L_u⁻¹, the map for −u is the inverse of the map for u, and
the Jacobian is L_u's. The tables are frozen between rebuilds, like every other proposal, and the
screen may use any deterministic surrogate.

**Toy check with exact draws.** Setup:

- Subjects (a, b) with y ~ N(log(eᵃ + eᵇ), 0.05²): a one-sided plateau in a, with b
  compensating along a curved ridge, as in RDM.
- Huang-Wand prior on V. Group, subjects and data are drawn from the joint, so each state is an
  exact posterior draw.
- 3000 replicates, 3 moves each, paired z-tests on about 50 functionals.

Results:

- Every proposed arm is consistent (max |z| ≤ 2.9): the current policy, every step forced,
  gated with 3 repeats, the fast mean path, ridge, and gated + ridge + fast.
- Both broken controls are detected: no Jacobian gives |z| = 28; a ridge shear that isn't
  volume-preserving gives |z| = 9.7.

**freeA15.** Three chains start from dispersed states (mu_A −2.1, −1.6, −1.0), run 100
re-adapting iterations, then 1000 production iterations. Two arms ran at a time, CPU from the
cgroup.

| arm | R-hat mu_A (all / 2nd half) | ESS mu_A | ESS lsd_A | ESS mu_B | ESS mu_t0 | CPU |
|---|---|---|---|---|---|---|
| ensemble (seed 91) | 1.41 / 1.11 | 6.7 | 48.9 | 8.1 | 15.3 | 2370 s |
| ensemble (seed 92) | 1.88 / 2.13 | 4.3 | 22.9 | 4.6 | 7.1 | 2417 s |
| shear, fixed plan | 1.66 / 1.75 | 4.9 | 14.0 | 6.7 | 13.0 | 3212 s |
| gated with 4 repeats, no ridge | 1.27 / 1.51 | 8.7 | 23.0 | 12.4 | 18.7 | 4491 s (uncapped) |
| **gated + ridge (seed 91)** | **1.11 / 1.06** | **29.1** | **101.5** | **49.6** | **158** | 3107 s |
| **gated + ridge (seed 92)** | **1.09 / 1.16** | **32.9** | 23.2 | **57.6** | **137** | 2878 s |
| gated + ridge, deep starts (−2.3, −2.2, −2.05) | 1.13 / 1.09 | 26.2 | 23.2 | 40.2 | 174 | — |

Gated with and without the ridge differ only by the ridge, which isolates it as the cause.
Per CPU-second, against the ensemble, mu_A mixes 3–6x better, mu_B 5–10x and mu_t0 7–17x.
log sd_A is still the weakest (ESS 23–100), and R-hat is not yet below 1.1 in every run.

**Other candidates:**

- *A subject-level move along the ridge* won't help. Per the README's frozen-group test, subjects
  already mix well given the group (ESS/iter 0.25–0.8 vs 0.02–0.04 jointly); the bottleneck is
  the coupling between levels.
- *Reparameterising b = A + B* is model-specific.
- *Joint HMC in non-centred coordinates* with finite-difference gradients would follow curved
  ridges, but costs about 2–4x per iteration. Keep it as a later, opt-in option.

## 4. Recommendations, ranked, each with a cheap test

1. **Keep 1aec215d. Rebaseline the evidence table with cgroup CPU and dispersed starts.** The
   current table compares buggy code, measures about 1/5 of the CPU, and trusts R-hat from
   warm-ups that all end near the same mode.
   *Test:* `plan_check.R` on a new warm-up state shows the stored plan equals the checked-fit
   plan. `tail_start.R` / `mix_start.R` from tail states, with the arm wrapped in
   `systemd-run --user --scope`.
2. **Gate the shear and add the fast mean path.** Directions with λ ≥ 0.9 (at most 2, optionally
   repeated) and scales with w ≥ 0.8. On easy problems it costs essentially nothing.
   *Test:* cgroup CPU ratio against the ensemble ≤ 1.05 on corr15 (measured 1.06, within noise)
   and on forstmann LNR, plus the toy invariance test.
3. **Replace the README's curved transport with the ridge transport for gated directions.** Build
   its tables only when some direction is gated in (21.7 s of CPU over 600 freeA15
   chain-iterations here; zero on corr15). Consider parallelising the build across subjects.
   *Test:* `path_ridge.R` on a second hard fixture (lba71sim, or a free-t0 GNG) should give a
   width ratio ≥ 1.5. `mix_start.R` on freeA15 seeds 1 and 3 should give R-hat ≤ 1.1 and
   mu_A ESS ≥ 3x the ensemble's per CPU.
4. **Make step sizes robust to rebuilds.** Steps are indexed by eigen-order, so a rebuild that
   rotates a direction inherits the wrong step, and production freezes it. Reset adaptation for a
   direction whose new vector's cosine with the old one is below 0.9. Never freeze a step adapted
   against a different plan.
   *Test:* per-direction production acceptance stays between 0.3 and 0.6 across rebuilds.
5. **Warn on weak identification.** When w_i ≥ 0.9 for most subjects and the plan has λ ≥ 0.95,
   warn that parameter i is prior-dominated and that mode-started chains under-disperse. Start
   chains spread along that direction.
   *Test:* on freeA15 the warning names A; with the spread starts, the ensemble's R-hat
   reflects its poor mixing instead of passing.
6. **Small fixes.**
   - `options(emc2.sampler_profile = TRUE)` aborts with the ensemble on: "unknown profile
     field(s): reject_ensemble_*". *Test:* any 5-iteration fit with profiling on.
   - Pool spin-wait (§1). *Test:* master CPU per round before and after.

## 5. Not established

- The ridge transport has been tried on one fixture (freeA15: two seeds plus deep starts), on the
  single slowest direction, and only for mean steps against linear frames. Its screen pass rate is
  still about 93%.
- Neither log sd_A mixing nor the funnel is solved. Scale steps keep their linear frames.
- Gating thresholds (0.9/0.8), the cap of 2, the 4 repeats and the 9-point grid were not tuned.
- Mixing-arm CPU was measured with two arms running concurrently. The clean per-iteration costs
  are in §1.

## Files

Session scratchpad:
`/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/011d1589-4aea-4fb2-9488-3fc543953fdd/scratchpad/`

- `wt/`: a detached git worktree with the overlay. It holds profile timers (`R/zz_sprof.R`,
  gated by `emc2.shear_prof`), the check-before-plan swap, and the options
  `emc2.shear_policy="gated"`, `_lam0`, `_w0`, `_reps`, `_max_dirs`, `_fast_mean` and
  `_ridge`. The diff is in `shear_overlay.patch`.
  Libraries: `lib_base` (HEAD), `lib_fix` (fix + options), `lib_fix3` (adds the ridge gate and
  the direction cap). Remove the worktree with `git worktree remove <scratchpad>/wt`.
- `prof/` (profile), `cost/` (clean cgroup costs), `ana/` (plan_check, a_profile, a_reference,
  curve_mode_path, path_decomp, path_oracle, path_ridge), `tail/` and `mix/` (dispersed and tail
  starts), `toy/` (toy_invariance.R, fast_eq2.R).
- `repro/shear/emc_freeA15_doubling_shear_2.rds`: the reproduced seed-2 run with all draws.
- `runs/fix_shear*/`: fixed-plan warm-up and production for freeA15 and corr15 seed 2, with
  `cgroup_cpu_*.txt`.
- `arm.sh`: warm-up plus production inside a cgroup scope.

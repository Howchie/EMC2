# RDMSWTN Gaussian sv > 0: sample stage is 3-5x slower than the gamma fit (handoff, 2026-10-20)

## Symptom

Two fits of the same design, data, cores and sampler build differ only in the drift law:

- `RDMSWTN(posdrift = TRUE, correlated = FALSE)` (Gaussian sv)
- `RDMSWTN(..., drift_dist = "gamma")`

The Gaussian fit's sample stage takes 6.4-18.9 min per 100 iterations, against 3.6-3.7 min for the gamma fit. That's 3-5x slower, and the block times vary erratically. It used all 30 sample tries (1845 iterations, ~6.5 h of sampling, ~7.5 h in total). At the end alpha max Rhat was 1.100 and mu max Rhat 1.060. This is not normal for this model size; the gamma fit is the reference.

Neither the per-call likelihood cost nor the total CPU explains it (see the evidence below). The time goes into wall-clock stragglers inside the shear rounds: one worker evaluates something slow while the others idle.

## Fit details

- Project: `/data/work/PM/NirvanaHons_Nback`.
- Script snapshots:
  - `Fits/run_scripts/rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gaussian_20261019.R`
  - `..._gamma_20261019.R`

  They are identical apart from `drift_dist`.
- Run settings: `RDM_VARIANT=noA,noS,sv,sepA,sepV,sepL,gL,t0L,blk`, 3 chains x 9 cores, `max_tries = 30`, `rhat_version = "new"`.
- Model: A = 0 (constant), s ~ 1, sv ~ 0 + lM (sv_lMTRUE, sv_lMFALSE), t0 ~ 0 + Load, posdrift = TRUE. 103 subjects, 20353 trials, LT 0.25 / UT 3.5.
- Package: playground HEAD 2369cb2f plus the uncommitted `.EMC_MODE_MAX_EXTRA_SD` guard in `R/mode_start.R` (see "Related fix" below).
- Gaussian output:
  - `Fits/Match_Control_Exp1_rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gaussian.RData`
  - log `Fits/rdm7r_..._gaussian.log`
- Gamma output: `Samples/samples_control_Exp1_rdm7r_..._gamma.RData`, in progress.

## Evidence

### Block timings (from the logs)

| Stage | Gaussian, min per block | Gamma, min per block |
|---|---|---|
| preburn | 4.7 | 5.1 |
| burn | 2.9-8.2 (first block 15.4) | 3.1-3.4 (first block 8.4) |
| sample | 6.4-18.9 | 3.6-3.7 |

### Sampler bookkeeping (`emc[[c]]$subject_hmc`)

| | Gaussian chain 1 / 2 / 3 | Gamma chains 1-3 |
|---|---|---|
| `alloc$c_round`, cost of one shear round | 8.92 / 2.01 / 14.31 | 1.43-1.48 |
| `alloc$hmc_fixed` | 2.71 / 0.60 / 2.23 | 1.39-1.41 |
| HMC CPU per iteration (`totals$cpu` / iterations) | ~8.6 s | ~10 s |
| HMC CPU per likelihood row | 0.88 ms | 1.37 ms |
| HMC acceptance | 0.69-0.71 | 0.65-0.72 |

- **Shear rounds are 1.5-10x more expensive in the Gaussian fit,** and the cost differs a lot between chains.
- **HMC is not the cause.** Its CPU per iteration and per row is lower in the Gaussian fit, and its acceptance is the same.
- **Per-subject particle cost is normal:** `c_part` median 0.19 s (Gaussian) vs 0.26 s (gamma).
- So total CPU is similar and the Gaussian wall time is lost to imbalance: some calls inside shear rounds (and probably collapsed rounds) are very slow, and the round waits on them.

### The likelihood near the posterior is fast

- Gaussian, at the final sample draws: about 1.0-1.4 ms per subject call (10-draw average, 4 subjects).
- Forcing sv_lMFALSE to 1e-2 ... 1e-6, or s = 0.1 with sv_lMFALSE = 1e-4, also gives ~1.0-1.3 ms.
- So it is not the near-zero mismatch sv by itself. The Gaussian posterior does put sv_lMFALSE near 0 (median 0.011, min 2e-4; the gamma fit has 0.32), which makes the lM = FALSE accumulator nearly deterministic in drift. The shear and collapsed moves then propose combinations far from that posterior.
- Scripts used: `plans/rdmswtn_gaussian_sv_slow_2026-10-20/svF_time.R` and `cost_cmp.R`. Run them from the project directory.

## Most likely mechanism (unconfirmed)

It is the same class of bug as the mode-search one fixed on 2026-10-19 (below). A move evaluates the likelihood at parameter vectors far from the posterior. There the Gaussian-sv survivor is finite but very slow:

- `rdmswtn_spv::log_posterior_mills` (called via `seed_log_T` / `post_log_TU` in `src/rdmswtn_spv.h`) took ~30 s per subject call at s ~ 3e-7, sv ~ 2e5, v ~ e^22.
- Its panel loop (`side()`) doubles panel widths when `sl^2 -> 1` and can cover hundreds of units.
- The Newton window search (`window_end`) and the 200-step mode loop also have no cost bound.

The gamma law has closed forms (`src/rdmswtn_gamma.h`) and does not slow down like this. That would explain why only the Gaussian fit is slow.

Candidate culprits, in order:

1. **Shear moves** (`R/shear_move.R`): proposals along group directions shift every subject at once. With a ridge in s, sv_lMTRUE, B and v, a long shear step can push some subjects to extreme s/sv. `c_round` is the cost that blows up.
2. **Collapsed moves** (`R/collapsed_move.R`): the same kind of joint proposal, with its cost folded into `C_rest`.
3. **HMC leapfrog positions** with finite-difference probes (h = 1e-4) are local, so they are unlikely; their per-row cost is fine.

## How to confirm

1. **Profile a live run.** Restart the Gaussian fit from the snapshot script with `RDM_TAG_SUFFIX=_gaussian_prof`, so the finished file is not overwritten. When CPU idles with one busy R process, run `perf record -F 300 -g -p <pid> -- sleep 10` and then `perf report --no-children --sort sym`. perf works on this machine; gdb isn't installed and ptrace is restricted. On 2026-10-19, `log_posterior_mills`'s lambda #2 plus `seed_log_T` were 95% of the straggler's time.
2. **Capture the parameters.** Wrap the likelihood entry the move uses (as `.emc_burn_ll_batch` was wrapped for the mode-search diagnosis) with `assignInNamespace` in a forked test. Save any call over 1 s together with its rows, then label the parameters. The pattern from 2026-10-19 is in `plans/rdmswtn_gaussian_sv_slow_2026-10-20/`: `replay_modes.R` (it replays every subject's mode search from a saved sampler state in `mclapply` and logs slow calls; it needs env `SP` set to an output directory) and `trace_mode.R`.
3. **Replay shear and collapsed rounds** from `Fits/Match_..._gaussian.RData` in the same way and time each subject's likelihood batch per round.

## Fix options

- **Likelihood side** (needed whatever the caller): bound the work in `log_posterior_mills`.
  - Cap the panel count and the window length, scaled with kPostDrop.
  - Fall back to the closed form, or a coarse bound, when the window is absurd.
  - Or reject (-Inf) parameter combinations whose drift and noise scales lie outside any plausible range (v/s or sv/s beyond ~1e6) before the quadrature.

  The cost of one call must be bounded. A 30 s call violates the pool's reply deadline and stalls a whole chain.
- **Move side:** apply the same prior-radius guard as `.EMC_MODE_MAX_EXTRA_SD` to shear and collapsed proposals. Reject them without evaluating the likelihood when any subject lands more than ~20 group-prior SDs from its current position. This keeps detailed balance if the rule depends only on the (symmetric) proposal geometry; check that for each move, or make the guard part of the target (treat those points as outside the support).
- **Regression test:** add one for each fix (the mode-search guard has one in `tests/testthat/test-burn-mode-start.R`).

## Related fix already made (2026-10-19, uncommitted)

`R/mode_start.R` `.emc_optimize_subject`: BFGS's first unit step along the raw gradient threw steep subjects 40-700 prior SDs away, into ~30 s likelihood calls. Burn then ran on one core. The fix returns 1e10 without evaluating the likelihood beyond `.EMC_MODE_MAX_EXTRA_SD` (20) Mahalanobis SDs past the start. Slow subjects went from 11-32 s to ~2 s with identical modes. The fix is installed to the main R library and covered by a test in `test-burn-mode-start.R`, but is **not committed**. It cured burn, not the sample-stage slowness described here.

## Also worth checking

- The Gaussian fit stopped at `max_tries` (alpha Rhat 1.100). Check whether the erratic shear cost made the allocator (`R/move_allocation.R`) drop or ration shear/collapsed rounds, which would hurt mixing as well as speed. Chain 3 has `tau` 10 and `cost` 95.5, against `tau` ~15-24 and `cost` ~20 elsewhere.
- 31 observations have Pareto k > 0.7 in the Gaussian fit (WAIC 19209.0, pWAIC 929.5). That may be a model issue rather than a sampler one, but compare it with the gamma fit when it finishes.

## Resolution (2026-10-20)

Confirmed by replaying the saved Gaussian fit: forced shear rounds from chain 3's sample draws (single core), and the rebuild's orbit-width probes. Every slow call had sv_lMFALSE transported to log sv = 14-60. There were two causes, and they compound.

1. **The Laplace fits extrapolate along flat directions** (`R/mode_start.R` `.emc_lik_canonical`). The posterior puts sv_lMFALSE near 0, and there the likelihood is flat in log sv. The likelihood Hessian has ~0 (often exactly 0 after the PSD projection) curvature along it, while `lin` carries the gradient P(x − mu). So the approximation's own mode sat at infinity: implied mode17 = 135 for one subject, and 71/103 subjects had singular precisions. The shear frame centre is (H + V⁻¹)⁻¹(lin + V⁻¹mu), and it slides toward that phantom mode as a move widens V. Log-sd steps on coordinate 35 (lsd sv_lMFALSE) and the orbit-width probes therefore transported subjects to sv ~ e^60.
   - Fix: raise the curvature in each eigendirection just enough to keep the approximation's mode within `.EMC_LIK_MODE_RADIUS_SD` (10) group-prior SDs of x. The gradient at x is kept, so the conditional mode under the fitting prior is unchanged.
   - Effect on the orbit probe for coordinate 35: 89.6 s → 2.7 s, worst call 87 s → 0.02 s, max |alpha| 167 → 41. The measured width was unchanged (2.11 vs 2.13).
   - Effect on forced shear rounds for chain 3: iterations with a call of 39 s or more went from 4 of 6 to 1 of 6.
   - Test: `test-laplace-conditional.R`.
2. **The truncation normaliser chased rounding noise** (`src/race_integrands.cpp` `get_trunc_normaliser_rowmajor_cpp`). A wide log-sd candidate still legitimately reaches sv/s ~ 1e12. There the sv > 0 CDF is noise at the 1e-12 level and is sometimes wrong outright (F(0.25) = 6e-5, F(3.5) = 1.7e-5, where both should be ≈ 1 − 1e-12). S(UT) then exceeds S(LT), the closed form is NaN, and the GSL fallback ran QNG/QAG at 200 and then 1000 subdivisions, natural and then log-shifted, per unique trial and winner. That was ~750k CDF calls and 75 s for one 8-row batch.
   - Fix: when S(LT) ≤ eps or S(UT) > S(LT), return the defective atom (or −Inf) without integrating. That is what the existing final check concludes anyway.
   - Effect on the 17 captured slow batches: 1.3-125 s → ≤ 0.07 s. All rows with log sv_lMFALSE < 5 are bit-identical. Rows at log sv ≥ 14 change only downward (more trials at min_ll), and those were already more than 700 nats worse than the sane candidates.

### Accuracy follow-up (2026-10-20, uncommitted)

Checked against an mpmath reference: the drift integral at 30 digits, sv 0.3-1e12, t 0.05-3.4, A = 0 and 0.5. There were three faults.

1. **`src/rdmswtn_spv.h`.** The posterior-Mills terms formed 1 − rho² and the posterior offset r − post_mean by cancellation. Through them, F broke down from sv/s ~ 1e4: error 5e-8 at 1e4, 6.5e-4 at 1e6, and outright wrong values by 1e8.
   - Fix: `Law::om = s²/Q`, formed directly, is carried into `log_posterior_mills` (curvature and recursion) and `post_log_TU`.
   - Now: log F is within 5e-15 up to sv = 1e12, log S within 1.2e-9 (an existing t = 0.05 point), and A > 0 within 6.4e-11 (unchanged).
2. **Truncation normaliser.** It took the scalar natural 1 − F path, because RDMSWTN had not opted into batched defective truncation. That path loses every survivor below the resolution of F.
   - Fix: opt in (`race_dispatch.cpp`). `rdmswtn_logS_at_t` now reports a NaN survivor as NaN (scalar fallback) instead of silently skipping it.
   - It also routes t = +Inf on the sv ≈ 0 branch to `wald_k0_log_surv`. The natural CDF is NaN there, which had counted the never-finish atom as 1.
3. **`src/particle_ll.cpp`.** Truncated trials were floored per accumulator: any component below min_ll set the whole trial to min_ll, even when Z was just as small and f/Z was moderate.
   - Fix: with truncation active, the kernels now return raw logs, and the floor applies after dividing by Z. Untruncated trials re-apply the per-component floor, so their values are unchanged.

Effect on the 136 captured rows:
- the 73 sane rows are unchanged (≤ 1e-10);
- rows at log sv 14-22 went from about −900…−4800 to −330…−650;
- rows at log sv 23-38 went from a floor of about −4600 to −320…−620, apart from rows where t0 > LT makes trials genuinely impossible.

The cross-variant check (posdrift × Erlang type × sv × A × UT) is identical to 9e-16.

Not changed: a normaliser below double epsilon (log Z < −36, reached at log sv ≳ 45) is still refused and sends the trial to min_ll. Most models' survivors, and RDMSWTN's kill/guess combination branches, are not resolved that far, and accepting such Z values could inflate the likelihood.

Cancellation test: `trunc_window_resolved` (`src/race_contract.h`) accepts the closed form S(LT) − S(UT) + atom only if the two survivors differ by more than sqrt(eps) relatively, or a retained atom carries more than sqrt(eps) of S(LT). Otherwise the trial integrates its density over the window. It is applied in the batched, scalar and correlated-node normalisers. Posterior draws are bit-identical, with no measurable cost. A 1e-14-wide window is now exact to 6 digits (before: 1.3e-2 off).

Test: `test-rdmswtn-extreme-sv.R`.

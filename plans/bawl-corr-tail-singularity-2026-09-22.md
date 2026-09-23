# Correlated BAwL returns huge positive log-likelihoods in a parameter corner

Date: 2026-09-22

Status: fixed and verified, uncommitted (see Resolution at the end). Present
in the build installed 2026-09-22 08:41 (after the RNG fix).

## What happened

Refitting the N-back trend7 LBA reference with independent RNG streams
(`/data/work/PM/NirvanaHons_Nback`, `BAwL(correlated = TRUE, posdrift = FALSE)`,
`filter_defective = TRUE`, LT = 0.25, UT = 4, shape-2 variance prior) produced
max R-hat of 4.7 → 6.6 during sampling. The group-level parameters were fine:
R-hat ≤ 1.07 for mu, ≤ 1.05 for sigma².

One subject, 52331, drove all of it. Its summed log-likelihood was about
**+10,900 in every chain** (10,938 / 10,767 / 10,998), i.e. about +53 per
trial. Every other subject's chains agreed within 1 unit. In the earlier
(healthy) saved fit, this subject's summed LL was −88.

The earlier fits under `RNGkind("L'Ecuyer-CMRG")` never reached this region,
which is plausibly because their chains shared recycled streams
(`plans/rng-stream-overlap-2026-09-22.md`). Once the sampler reaches it, the
region is an absorbing attractor.

## Evidence that the value is wrong

At the offending draw, ordinary Low-load Novel correct trials (RT 0.50–0.60 s,
no duplicate-RT oddities) score **+280 to +344** each. A conditional density
cannot reach e^343 at many different RTs; it must integrate to 1 over the
RT window.

Mapped parameters for those trials:

| accumulator | v | sv | B | A | b | t0 | rho (loading) |
|---|---:|---:|---:|---:|---:|---:|---:|
| matching (winner) | 3.553 | 0.096 | 0.005 | 0.471 | 0.475 | 0.052 | +0.92 |
| mismatching | −1.105 | 1.000 | 0.001 | 0.471 | 0.471 | 0.052 | −0.92 |

The pair correlation is −0.92 (sampled probit −1.751). With sv ≈ 0.1, the
winner essentially always finishes by about 0.19 s, so an RT of 0.5 s is deep in
the tail (z ≈ −26).

| worst trial (RT 0.5) | log-lik |
|---|---:|
| correlated route, as fitted | **+343.7** |
| identical parameters, rho = 0 (independent route) | −23.0 (the log 1e−10 floor) |
| analytic independent race density at RT 0.5, unnormalised | −340.5 |
| `filter_defective` removed, correlated route | +343.7 (unchanged) |

The correlated value is almost exactly **minus the true log density**. That
suggests a sign flip or an inverted ratio in a log-space fallback, for example
one used when the density or a normaliser underflows. `filter_defective`
does not matter here, so the finite-response normaliser is not the culprit. The
finite-UT truncation normaliser has not been ruled out.

## Bisection

Starting from the singular draw, each parameter was reset in turn to the
subject's value from the healthy fit. The subject total was 11005 at the
singular draw and −83 at the healthy draw.

- **Resetting any one of these removed the spike** (total then −940 to −3528):
  `v_LoadLow`, `v.e_novel`, `B`, `B_lR_target`, `B_LoadLow`, `A`, `t0`,
  `t0_LoadLow`, `sv_lMTRUE`, `rho`.
- **Resetting these did not** (total stayed +5,000 to +11,000):
  `v_LoadHigh`, `strength_*`, `v.phi`, `v.fam`, `B.w_prac`,
  `t0_repeatTrialRepeat`, `rho_coupled_repeats`.

So the trigger is a conjunction on Low-load Novel trials: tiny B, winner
sv ≈ 0.1, t0 at 0.052 (near its 0.05 lower bound), strong negative rho and a
high winner drift.

## Minimal reproduction attempt (does not trigger)

One subject, eight trials (RT 0.2–0.6), `BAwL(correlated = TRUE,
posdrift = FALSE)`, `v ~ lM, sv ~ lM, B ~ 1, A ~ 1, t0 ~ 1, rho ~ 1`, at the
same natural values (rho ±0.92), with and without LT = 0.25 / UT = 4, and with
and without `filter_defective`. It returned sane values (floored at −23) in
every combination. The full model therefore adds something needed for the
fault. Candidates:
- per-accumulator thresholds (B 0.005 vs 0.001);
- the custom trend kernels (the drift trend on `v`, the practice trend on `B`);
- the `rho ~ coupled_repeats` structure;
- compressed / unique-trial evaluation (`calc_ll_oo_pw`).

## Reproduction with the real object

From `/data/work/PM/NirvanaHons_Nback`:

```
Rscript Fits/debug/reproduce_bawl_corr_singularity.R
# subject total ll: correlated 11005.4 | rho = 0 -2598.8
# worst trial 42: correlated 343.7 | rho = 0 -23.0
```

The script uses the samples in
`Fits/debug/lba_trend7_vs2_mt_singular_samples.RData`. That file is kept
outside `Samples/` so it cannot be mistaken for a result.

## Suggested next steps

1. Evaluate the correlated pair density and survivor for the two-row trial at
   the table's values directly in the C++ entry point, splitting natural-space
   and log-space paths. Look for a branch that returns `-log(x)` or divides by
   an underflowed quantity.
2. Add a guard or test: for any parameter vector, a trial's conditional
   log-density must not exceed `-log(bin width)` for observed-RT resolution,
   or at least must not flip sign relative to the rho = 0 route in the far tail.
3. Separately (modelling, not a package bug), t0 sitting at the 0.05 bound and
   winner sv → 0.1 is a degenerate LBA corner. A prior or bound on sv may be
   worth considering for this model even after the numeric fix.

## Resolution (2026-09-22, uncommitted)

Status: fixed and verified. Subject 52331 at the singular draw now totals
-2657.0 (was +11005.4; rho = 0 gives -2598.8), and its largest trial is +0.2.

### Cause

Two garbage values, one divided by the other. The exact pair kernel
(`src/bawl_corr_exact.h`) builds every survivor and cause density by
subtracting bivariate normal CDF corners. Those corners have ABSOLUTE
accuracy only: tvpack is good to ~3e-16 absolute and the Drezner branch of
`norm_cdf_2d_hybrid` to 3.8e-7. At x ~ -26, rho = -0.92, tvpack cancels
Phi(x)Phi(y) ~ e^-343 against an integral of the same size and leaves noise
around 1% of that. So, for trial 42 (true values from a 256-bit reference):

| quantity | exact route | truth |
|---|---:|---:|
| log cause density at RT .5 | -348.0 | -1904.4 |
| log S(LT = .25) | -691.8 | -280.45 |
| trial ll | +343.7 | -1624 (floors at -23.03) |

"Almost exactly minus the true log density" was a coincidence. The
normaliser cross-check that should have caught the bad Z compared it with
the numeric pair route, but that route only scanned +/-12 marginal SD. The
mass here sits 12-26 SD out, so it returned 0, and "keep the larger valid
estimate" kept the garbage. The minimal reproduction needs per-accumulator
B (with a common B the residual is an exact 0, which floors):
`test-bawl-correlated.R`, "cancellation residuals cannot turn an impossible
trial into a spike".

### Fix

- Every exact pair value carries an absolute error bound
  (`BAwLCorrPairResult::noise`): per-corner bounds (tvpack 1e-15; Drezner
  banded by |rho|, 1e-12 to 1e-6), propagated through the rectangle moments
  and the weighted terms, and accumulated over the numerator and Z.
- A trial is trusted only if both its numerator and its Z clear their bound
  by 1e3. Otherwise it is retried on tvpack corners, then on the numeric
  route; routes are never mixed within a trial. An unresolved numerator
  whose upper bound is already below the floor is floored without a retry.
  This replaces the "keep the larger Z" cross-check.
- The numeric route no longer uses a fixed window. Every pair integrand is
  log-concave in the first racer's drift, so the route scans for the mode,
  finds where each flank falls 36 nats, and integrates with adaptive G7/K15.
  It matches a 256-bit reference to ~1e-8 even at e^-696.
- `pnorm_std(x, lower = FALSE)` returned 1 - Phi(x) under USE_FAST_PNORM:
  no relative precision past x ~ 5 and zero past 8.3. It now returns
  Phi(-x). This also fixes the natural-space helpers in `wald_functions.h`
  (the point-start Wald CDF was -3% off at an early RT).

### Verification

- 700 random pair configurations (posdrift on/off, |rho| up to 1 - 1e-4,
  sv .05-2, A 1e-3-2, deep-tail RTs) against the 256-bit reference: every
  value the exact route accepts is within 1.9e-6 (log); the numeric route is
  within ~3e-6 wherever the truth is representable and not below e^-690.
- Healthy lba_trend7_vs2 posterior (40 subjects x 10 draws): 134 of 81,370
  trial values move, 99% of them by < 7e-5; the largest move (2.3e-3, a
  3.68 s RT) is the old value being wrong: the new one matches the reference
  to 1e-9.
- CPU cost (min of 3 interleaved runs, loaded box): +2% at posterior draws,
  +6% at particles spread 3x the posterior SD. Numeric fallbacks fall from
  3,240 to ~330 per 2.4M trial evaluations, because the tvpack tier absorbs
  what the old unstable guard sent to the numeric route.
- Full suite at EMC2_TEST_LEVEL=full: no new failures. The remaining 31
  failures fail identically without this fix. 29 are stale snapshots from the
  RNG-stream fix (committed HEAD plus only R/sampling.R, R/fitting.R and
  R/chain_pool.R reproduces all of them). 2 fail on HEAD too: seed-dependent
  Monte Carlo checks that change once another file has set L'Ecuyer-CMRG.

Library with the fix: /data/tmp/Rlib_bctail_final. It is not yet installed
into the default library, because fits were running against that one.

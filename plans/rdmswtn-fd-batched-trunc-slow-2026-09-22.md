# filter_defective makes RDMSWTN about 1.7x slower at a finite UT

Date: 2026-09-22

Status: fixed and verified via closed-form wald_k0_log_surv_closed in src/wald_functions.h.
Performance only: the log-likelihood values are identical with and without it.

## Symptom

The N-back trend7 natural-scale RDM (`RDMSWTN(posdrift = FALSE)`,
`TC = list(filter_defective = TRUE)`, LT = 0.25, UT = 4) samples at 56–66 s
per 100 iterations. The posdrift = TRUE fit of the same design, with no
filter_defective, takes 36–44 s. The difference is not the natural scale
or `_IO`. It comes from `filter_defective` combined with a finite UT, and the
positive-drift model slows down in the same way when the flag is added.

## Measurements (Before Fix)

Setup: `EMC2:::calc_ll_manager` (the sampler's route), 40 subjects, 60
posterior alpha draws each, median of 3 runs, single core. Tested on the
saved natural-scale and positive-drift fits in `/data/work/PM/NirvanaHons_Nback`.

| model | flag | UT | time | summed ll |
|---|---|---:|---:|---:|
| posdrift = TRUE | none | 4 | 0.69 s | −4101.5 |
| posdrift = TRUE | fd attr added | 4 | **1.18 s** | −4101.5 |
| posdrift = TRUE | none | Inf | 0.56 s | −4135.4 |
| posdrift = FALSE | fd | 4 | **1.21 s** | −4055.4 |
| posdrift = FALSE | fd dropped | 4 | 0.67 s | −4055.4 |
| posdrift = FALSE | fd | 2 | 0.69 s | −3503.8 |
| posdrift = FALSE | fd dropped | 2 | 0.67 s | −3503.8 |
| posdrift = FALSE | fd | 1.2 | 0.70 s | −788.9 |
| posdrift = FALSE | fd | Inf | 0.54 s | −4092.9 |

The pointwise route (`calc_ll_pw`) shows no difference. The likelihood values
agree to 0.1, because with one positive-drift accumulator the never-finish
atom is about 0.

## Measurements (After Closed-Form Fix)

Tested with `bench_mgr3.R` on `samples_control_Exp1_rdm_trend7_s_natural.RData`:

| model | flag | UT | time | summed ll |
|---|---|---:|---:|---:|
| posdrift = FALSE | fd | 4 | **0.69 s** | −4092.2 |
| posdrift = FALSE | fd dropped | 4 | 0.67 s | −4092.2 |
| posdrift = FALSE | fd | 2 | 0.68 s | −3532.9 |
| posdrift = FALSE | fd dropped | 2 | 0.68 s | −3532.9 |
| posdrift = FALSE | fd | 1.2 | 0.67 s | −788.9 |
| posdrift = FALSE | fd dropped | 1.2 | 0.64 s | −788.9 |

The 1.8x slowdown at UT = 4 s is completely eliminated (0.69 s vs 0.67 s).

## Cause

`particle_ll.cpp` (around lines 2505–2520) decides which normaliser route to take:
- **Without the flag,** RDMSWTN is `defective_upper_tail` with a finite UT and
  has no `supports_batched_defective_truncation`. It therefore takes the
  per-trial scalar route (`get_trunc_normaliser_rowmajor_cpp`, which uses the
  natural `cdf1`).
- **With the flag,** `defective_finite_UT` is false, so the batched route
  `rdmswtn_logS_at_t` is used.

In the batched route, `rdmswtn_k0_logsurv` (`model_RDMSWTN.cpp:154`) checks the
natural CDF. When it is at least `1 - EMC2_CDF_SAT_MARGIN` (1e-8), it calls
`wald_k0_log_surv`, which ran `wald_k0_log_avg_over_d`: 4 panels × 20-node
Gauss–Legendre, i.e. 80 log-space point-survivor evaluations plus
`log_sum_exp` for each accumulator.

At UT = 4 s it ran on almost every trial: a winning accumulator with a strong
drift has F(4 − t0) ≈ 1.

## Fix Implemented

Added `wald_k0_log_surv_closed` and `wald_k0_surv_exp_antideriv` in `src/wald_functions.h`.
The start-point average decomposes into:
  A * S_A(t) = J1 - I2
where J1 is evaluated via `log_normal_phi_integral()`, and I2 = (T(d_hi) - T(d_lo)) / (2 mu).
When z1 > 0, T(d) = V(d) - 1 where V(d) = phi(v) * (R(x1) + R(x2)), cancelling the
-1 identically upon differencing and eliminating catastrophic cancellation.
`wald_k0_log_surv` uses the closed form with `wald_k0_log_avg_over_d` retained as a
safe fallback if outer differencing ever loses precision.

## Reproduction

The benchmark scripts are in the session scratchpad (`bench_mgr3.R`). In short:
load both fits, take posterior alpha draws, and time
`EMC2:::calc_ll_manager(P, dadm, model)` while toggling
`attr(dadm, "emc2_filter_defective")` and `dadm$UT`.

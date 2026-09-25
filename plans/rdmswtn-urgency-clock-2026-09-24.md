# Plan: RDMSWTN under an accelerating (urgency) clock

Date: 2026-09-24

Status: linear and exponential clocks implemented in EMC2. Refactored the same day (see "Implementation as built") so RDMSWTN_TT and RDMSWTN_UT share one clock layer over the ordinary RDMSWTN kernels.

## Motivation

In the N-back PM project (`/data/work/PM/NirvanaHons_Nback`) the working model is an RDMSWTN race with variance-additive noise, sigma^2 = s^2 (1 + rho v). This is the constant-theta case of Smith's two-source diffusion. It matches median RTs and choice proportions in every Load x trialType cell. It overshoots the slow tail (q90) on lure trials, and the cause is known:

- The baseline variance term adds relative spread in proportion to 1/v.
- Lure accumulators have the weakest correct-response drift, so they get the most. In the High-load lure cell the baseline supplies 56% of the noise variance and CV is 0.70, against 0.51-0.57 elsewhere.
- Removing the baseline (pure Poisson, constant CV) fixes q90 but destroys accuracy and medians.

So the remaining lever is the within-trial time course, not the noise law. A gain that rises during the trial ("urgency") shortens long decisions most and barely touches fast ones.

`RDMSWTN_TT` already implements an exact time change of the RDMSWTN process, but its clock *decelerates* (exhaustion: gain 1 - x/tau falls to zero and produces omission mass). Urgency needs the mirror image: an accelerating clock with no plateau.

## Mathematical kernel

Let x = t - t0 be decision time and theta(x) > 0 an operational clock rate. Operational time is q(x) = integral_0^x theta(z) dz. If Y is an ordinary RDMSWTN accumulator (constant v, b, A, s, sv), then X(x) = Y(q(x)) exactly. Conditional on its drift draw, the physical-time SDE is `dX = v theta(x) dx + s sqrt(theta(x)) dW_x`, so infinitesimal variance is multiplied by theta and diffusion standard deviation by sqrt(theta), not by theta. The same rule applies to every additive diffusion noise source. Thus

    f(t) = f_RDMSWTN(q(x)) * theta(x)
    S(t) = S_RDMSWTN(q(x))
    F(t) = F_RDMSWTN(q(x))

These are the same identities `drdmswtn_tt` / `prdmswtn_tt_log_surv` use. Only q, theta and q^{-1} change.

Conditions:
- The clock is global within an accumulator. It scales drift and all infinitesimal noise variances by the clock rate. The time change is exact for any s, including trend-driven s (e.g. sigma^2 = s^2(1 + rho v)).
- If the baseline noise is meant *not* to be gain-modulated (Smith's premature-sampling reading), the construction is no longer exact and a numerical first-passage solver would be needed. That is out of scope. Document that the urgency gain multiplies the whole accumulator.
- Start-point variability A and between-trial drift variability sv compose unchanged, because they act on Y.
- With posdrift = TRUE and q unbounded, there is no new omission mass. Under posdrift = FALSE the intrinsic defective mass is unchanged: S(inf) = S_Y(inf).
- The clock is per accumulator (per row), so the design can let it vary by accumulator. The principled default is one gain per trial, shared by all accumulators (`u ~ 1` or `u ~ Load`).
- The clock does not replace a rate-dependent noise trend on s. The trend is a law across accumulators at a fixed time (sigma^2 = s^2(1 + rho v)); the clock scales every accumulator equally over time. The combined process is dX = v theta dt + s sqrt((1 + rho v) theta) dW. Keep the s trend feeding the `s` column; `u` only adds the clock. Without the trend, an accelerating clock would reduce right skew further in the high-drift cells that constant noise already makes too peaked.

This time change is distinct from multiplicative state urgency `Z(x) = U(x)Y(x)`. The latter has a pathwise equivalent bound `b/U(x)` for the ungated process. Smith and Ratcliff (2022, doi:10.1037/rev0000301) explicitly distinguish this state transformation from a time transformation: when urgency also transforms the time coordinate, the collapsing-bound equivalence does not follow. The linear and exponential clocks below therefore should not be described as particular collapsing-bound shapes.

### Not equivalent to collapsing bounds

Two different operations both get called "urgency". This model implements only the first.

- **Rate gain (this model):** theta(t) multiplies the input rates, dX = v theta dt + s sqrt(theta) dW. This is an exact time change, so choice probabilities and the quantile CAF are unchanged; only the RT axis is warped.
- **State gain (urgency-gating):** the decision variable is u(t) X(t) against a fixed criterion b. This is exactly a collapsing boundary b/u(t) on X (Smith & Ratcliff, 2022, Psychological Review, integral-equation treatment of urgency-gating and collapsing boundaries). It changes choices: late responses become less accurate. It has no closed-form RDM density except in special cases.
- **Linear collapse is already covered:** for a single-boundary Wald accumulator, a linearly collapsing boundary b - c t equals adding c to the drift. In a race, linear collapse is a common additive drift boost that the existing drift parameters absorb.
- **Diagnostic:** the quantile CAF separates the two. In the N-back `_d42` data, Novel and Lure CAFs match the stationary model with no extra late accuracy fall. That favours rate gain over global collapse. Target trials show an extra late fall (High load: -12 vs -3 points from quintile 3 to 5), which may need a Target-specific account instead.

### Primary clock: linear urgency (implement first)

    theta(x) = 1 + u x,  u >= 0 (units 1/s)
    q(x)     = x + u x^2 / 2
    q^{-1}(y) = 2y / (1 + sqrt(1 + 2 u y))     (stable form, as in rdmswtn_tt_qinv)
    log theta = log1p(u x)

- u = 0 is exactly RDMSWTN. One parameter, cleanly nested.
- Exact inverse, so the simulator stays closed form.
- Interpretation: gain doubles after 1/u seconds of decision time.

## Nonlinear clock options

All of these are viable. Any theta(x) > 0 with a closed-form integral gives an exact density. The simulator needs q^{-1}: closed form where available, otherwise a guarded monotone Newton or bisection solve (q is strictly increasing and convex for rising gain, so Newton from y is safe).

| Clock | theta(x) | q(x) | q^{-1} | Params | Nests RDMSWTN at | Notes |
|---|---|---|---|---|---|---|
| Linear (primary) | 1 + u x | x + u x^2/2 | closed | u | u = 0 | Simplest; moderate tail compression. |
| Exponential | e^{kappa x} | (e^{kappa x} - 1)/kappa | log(1 + kappa y)/kappa | kappa | kappa = 0 (use expm1/log1p) | Mirror of FRQfade's A_kappa with the sign flipped. Very strong truncation of long tails; may over-shorten lures. |
| Delayed linear (hinge) | 1 + u (x - d)_+ | x + u (x - d)_+^2 / 2 | closed, piecewise | u, d | u = 0 | Urgency only after d s. Leaves the whole fast range untouched. Kink in the density's derivative at x = d; the density itself is continuous. |
| Saturating (hyperbolic) | 1 + u x/(x + h) | x + u (x - h log1p(x/h)) | Newton | u, h | u = 0 | Gain rises then plateaus at 1 + u. Bounded urgency; avoids runaway gain for very slow trials. |
| Power increment | 1 + u k x^{k-1}, k >= 1 | x + u x^k | closed only for k = 1, 2; Newton otherwise | u, k | u = 0 | k = 2 is the linear clock. Only fit k if linear/exponential both fail; u and k trade off. |

Implementation scope: `clock = "linear"` and `clock = "exponential"`. The delayed, saturating, and power shapes are excluded from this model's current API.

Pure power clocks q = x^gamma (theta = gamma x^{gamma-1}) are **not** recommended: theta(0) is 0 or infinite for gamma != 1, which changes the leading edge as well as the tail and confounds t0.

## Model contract

New sibling model; `RDMSWTN_TT` stays unchanged so existing fits remain reproducible.

- Constructor: `RDMSWTN_UT(posdrift = TRUE, correlated = FALSE, correlate = c("times", "drifts"), clock = c("linear", "exponential"))`, mirroring `RDMSWTN_TT()` in `R/model_RDM.R`.
- C++ names: `RDMSWTN_UT`, `RDMSWTN_UT_IO`, `RDMSWTN_UT_CORR`, `RDMSWTN_UT_CORRD`; exponential variants add `_EXP` before the correlation suffix.
- Parameters: v, B, A, t0, s, sv, u.
- u: transform `exp`, default `log(0)` so a design that omits it is exactly RDMSWTN, `minmax` `c(0, Inf)`, and listed in `exception` like A/sv (0 allowed).
- `p_types_canonical = c("v","B","A","t0","s","sv","u")`.
- Column registry: new `namespace rdmswtn_ut { enum : int { v = 0, B, A, t0, s, sv, u, N_REQ }; }` in `src/col_registry.h`. u is a **required** column (see the positional-column caveat in the FRQfade plan).

## Implementation as built

The first version copied the RDMSWTN_TT adapters, simulators and constructor. The review found that the package already had a generic operational-time warp (`src/time_warp.{h,cpp}`, installed over an adapter's five entry points for the ballistic `eta`). Both TT and UT are now that warp over the ordinary RDMSWTN kernels:

- **`src/time_clock.h`** (new, no Rcpp dependency) holds the whole clock family in `emc2tw::`: the ballistic power clock (unchanged, bitwise), `CLOCK_LINEAR`, `CLOCK_EXPONENTIAL` and `CLOCK_EXHAUSTION` (TT). Each has an identity value, a validity rule, `clock_fwd`, `clock_log_jac`, `clock_inv` and `clock_budget` (finite only for exhaustion).
- **`TimeWarpPlan`** (`src/race_contract.h`) gained `clock`, `par_name` and `required`. `configure_time_warp_context()` resolves the clock column by name. A required clock (TT `tau`, UT `u`) must be present.
- **Wrappers** (`src/time_warp.cpp`) are clock-generic. Two new cases:
  - an invalid clock parameter zeroes the row's density and survivor;
  - `logS_at_t(+Inf)` on the exhaustion clock is evaluated row-wise through the raw survivor at `t0 + tau/2`, because the common-t trick cannot represent per-row frozen times.
- **Dispatch** (`src/race_dispatch.cpp`): one RDMSWTN branch. TT/UT use the RDMSWTN adapters, register no timer columns, and set the clock. The per-model TT/UT adapters and the `ut_exponential` context flag are gone.
- **R-facing scalars** (`src/model_RDM.cpp`): `drdmswtn_clock` / `prdmswtn_clock`, the vectorised `dRDMSWTN_clock_cpp` / `pRDMSWTN_clock_cpp`, and `rdmswtn_clock_qinv`.
- **Simulators** (`src/model_rng.cpp`): `rrdmswtn_clock_cpp`, `rrdmswtn_clock_drift_corr_cpp` and `rrdmswtn_clock_corr_cpp` take the clock code. `draw_pair_copula_uniforms()` and `validate_rdmswtn_copula_rows()` are shared with the plain RDMSWTN copula. LNR's copula draws in a different RNG order and was left alone.
- **R** (`R/model_RDM.R`): `.rdmswtn_clocks` (code, parameter, defaults, bounds) and `.rdmswtn_clock_model()` build both constructors. `dRDMSWTN_clock` / `pRDMSWTN_clock`, `.rRDMSWTN_clock` / `.rRDMSWTN_clock_corr`, `.rdmswtn_pair_copula_uniforms()` and `.resolve_rdmswtn_race()` serve TT, UT and the RDMSWTN copula. `.rfun_RDMSWTN_clock()` is in `R/model_rng.R`.

Consequences, measured on the likelihood:
- TT and UT inherit RDMSWTN's closed-form sv = 0 kernels. At sv = 0, UT went from 1.49x to 1.15x the cost of plain RDMSWTN, and TT improved the same way.
- The survivor route is now RDMSWTN's own for both models.
- Totals agree with the pre-refactor implementation to 6e-13 over 81 TT/UT configurations (finite, omission, truncated; independent, copula, correlated drifts). The ballistic `eta` warp is bit-identical.

## Tests (`tests/testthat/test-rdmswtn-ut.R`, modelled on test-rdmswtn-tt.R)

1. The constructor exposes the contract: parameter names, defaults, transforms, the u = 0 exception.
2. u = 0 reproduces RDMSWTN exactly (density, CDF, log-survivor, race likelihood) for posdrift TRUE/FALSE, with and without A and sv.
3. Operational-time identities:
   - numeric d/dt of F equals f;
   - integral of f over (t0, inf) equals 1 (posdrift) or S_Y defect (IO);
   - F(t) = F_RDMSWTN(q(x)) on a grid.
4. Clock inversion is stable at small and large u and y: q(qinv(y)) = y to 1e-12.
5. Compiled and R reference simulators agree with the analytic CDF (KS test / quantiles), independent and correlated.
6. `make_data` dispatches through the compiled simulator.
7. The compiled likelihood handles finite, truncated (LT/UT) and `filter_defective` trials, matching a direct R computation.
8. Correlated versions nest rho = 0 and reject malformed designs, as in the TT tests.
9. Monotonicity: for fixed parameters, increasing u shortens every RT quantile, most at high quantiles.

## Validation in the N-back project (after the package passes its tests)

- **Search-script flag:** e.g. `uG` in `PMDC_control_Exp1_rdm_trend7_search.R` sets `model = RDMSWTN_UT(...)` with `u ~ 1` (then `u ~ Load`), on top of `t0R, vs2, noA` and the variance-additive noise trend.
- **Prior:** log u ~ N(log(0.5), 1). A gain doubling after ~2 s is moderate relative to 0.5-2 s decision times.
- **Parameter recovery first:** simulate from the fitted `rdm7s_noA_t0R_vs2_mt_d42` parameters plus a known u, refit, and check u is identified jointly with B, s and rho. Urgency and threshold both shape tails, so identifiability is the main risk.
- **Compare against `rdm7s_noA_t0R_vs2_mt_d42`** on:
  - ordinary and equal-cell-weight LOOIC (`cell_weighted_elpd.R`);
  - the lure q90 and the medians/accuracy in the subject-predictives overlay (`overlay_predictives.R`).
- **Success criterion:** lure q90 overshoot removed without degrading medians (currently within 15 ms in all cells) or accuracy.

## Risks

- **Identifiability:** u against B and the noise law. Mitigate with recovery simulations before interpreting u.
- **Boundary behaviour:** if u -> 0 in the posterior, the log-scale sampling puts mass near the boundary. That is fine and interpretable ("no urgency"), but R-hat may be slow.
- **Numerics:** q(x) grows quadratically, so for large u and slow trials the RDMSWTN kernel is evaluated at large operational times. The existing kernels already handle long times, but add a large-q test.

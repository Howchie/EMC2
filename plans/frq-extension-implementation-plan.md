# FRQ extension implementation plan

> **Superseded parameterisation note (2026-09-20).** The implemented FRQ
> contract now fits the generative coordinates `(alpha, beta, p, lambda, t0,
> delta, cv_u)`. `h` is derived for reporting and the former `tau` conditional
> median chart is no longer exposed. The historical `(h, tau)` inversion notes
> below describe the pre-refactor design and should not be used for new model
> formulas.

Date: 2026-09-20

Scope: extend the existing FRQ race with analytic threshold variability, unit-level
rate frailty, and an explicit proper (`p = 1`) boundary. The implementation now
follows this plan; the document records the model decisions and validation contract.

## Recommendation

Use one nested FRQ family

\[
F(u) = H_\delta\{I_{pG_{c}(u)}(\alpha,\beta)\},
\]

where `H_delta` is the existing continuous criterion-variability map and

\[
G_c(u) = 1 - (1 + c^2\lambda u)^{-1/c^2}.
\]

`c` is the coefficient of variation of an individual unit's registration rate.
The limit `c = 0` is the current exponential-registration FRQ. Equivalently,
each unit has an independent rate

\[
R_i \sim \operatorname{Gamma}(a, a/\lambda), \qquad a = 1/c^2.
\]

This preserves `lambda` as the mean unit-registration rate, retains arbitrary
continuous `alpha` and `beta`, and stays closed form. It is not the same latent
mechanism as PCOUNTER's trialwise common rate draw, which is intentional.

Do not add PCOUNTER's integer geometric threshold (`omega`) to FRQ. The existing
`delta` construction is the better continuous replacement: it is already
implemented analytically, has a closed inverse, and does not require integer
thresholds. Do not add `phi`/`kappa`; they would be a different coordinate chart
and, especially when shared across racers, impose a different cross-racer model.

## What the repository already has

The current FRQ implementation already supplies the important pieces of the
threshold extension:

- `delta` is a uniform log-odds criterion spread with analytic `H`, `H^{-1}`, and
  `H'` in `src/model_FRQ.h`.
- `h` and `tau` are inverted through `H^{-1}` and two beta quantiles, so they
  retain their meanings when `delta > 0` (`R/model_FRQ.R`).
- The compiled density, CDF, survivor, truncation callback, and C++ simulator
  all carry `delta` (`src/model_FRQ.cpp`, `src/model_rng.cpp`).
- The FRQ cache is shared by the race callbacks; its key will need one more
  extension parameter when frailty is added.

The gaps are that unit frailty is absent, `h = 1` is rejected by the model
bounds even though the C++ kernel can evaluate it, and several comments still
describe FRQ as unconditionally defective.

## Parameterisation decisions

### Keep the current fitted chart

The public sampled chart remains

`(alpha, beta, h, tau, t0, delta, cv_u)`

with `cv_u = 0` as the exact no-frailty member. `tau` remains the conditional
median decision time. The existing derived `lambda` remains the mean unit rate
reported by `Ttransform`; it is not replaced by `phi`, and a new free `lambda`
column is not added because that would change the current `h`/`tau` contract and
create a scale migration.

If “preserve lambda” instead means that users must design directly on a sampled
`lambda`, that is a separate public reparameterisation (`(p, lambda)` or
`(h, lambda)`) and should be designed and versioned separately rather than
silently mixed into this extension.

### Unit-frailty parameter

Expose `cv_u`, the CV of the latent rate of one evidence unit. Internally use
`rho = cv_u^2` and `a = 1/rho`.

- transform: `exp`;
- default: `log(0)`;
- exact exception: `cv_u = 0`;
- keep it out of `p_types_canonical`, like `delta`, so existing designs remain
  unchanged and the extension is normally fixed/shared initially.

This is preferable to fitting `a` directly: the ordinary exponential model is
an exact, reachable boundary instead of the unattainable limit `a = Inf`.

### Proper versus defective

Allow `h = 1` as a bound exception while retaining the free-parameter upper
bound just below one. With the existing inversion, `h = 1` gives `p = 1`, so
availability is fixed to one and the accumulator is proper. This is a fixed
parameterisation, not a value the sampler should approach through a probit
intercept. Keep the default at `h = 0.95`, so the current defective model is
unchanged unless the user explicitly fixes `h` to one.

The race context may continue to advertise a possible defective upper tail:
when `h = 1`, `log S(+Inf) = -Inf`, so the existing atom/truncation code becomes
a no-op for that row. Update the misleading “always defective” comments.

## Analytic core to implement

### Registration CDF and density

For finite `c > 0`, evaluate on stable log scales:

\[
\log S_R(u) = -\log1p(c^2\lambda u)/c^2,
\quad G_c(u) = -\operatorname{expm1}(\log S_R(u)),
\]

\[
\log g_c(u) = \log\lambda -(1/c^2+1)\log1p(c^2\lambda u).
\]

Branch exactly at `c = 0` to the current exponential formulas, preserving the
existing zero-frailty results and fast path. At `u = Inf`, use `G = 1` and
`S_R = 0` explicitly.

### `(h, tau)` inversion

Let

`z_h = H_delta^{-1}(h)`, `z_tau = H_delta^{-1}(h/2)`,
`p = qbeta(z_h, alpha, beta)`, and `u_tau = qbeta(z_tau, alpha, beta)`.

For `c = 0`, retain the current

`lambda = -log1p(-u_tau / p) / tau`.

For `c > 0`, use

`lambda = expm1(-c^2 * log1p(-u_tau / p)) / (c^2 * tau)`.

This keeps `tau` as the conditional median and `lambda` as the mean unit rate.
Reject non-finite or collapsed quantiles exactly as the current `frq_derive`
does; additionally reject non-finite `c^2`.

### Density, CDF, and survivor

Generalise `FrqState` from the exponential `q = p(1-exp(-lambda*u))` to
`q = p*G_c(u)` and its cancellation-safe complement
`1-q = (1-p) + p*S_R(u)`. Replace only the registration-chain part of the
log-density with `log(p) + log(g_c)`. The incomplete-beta and `H_delta` paths
remain unchanged.

The survivor still uses the beta reflection and the symmetric `H_delta` map;
there is no subtraction from one. This preserves deep-tail accuracy and the
defective atom `1-h`.

## File-level work

1. `R/model_FRQ.R`
   - Add `cv_u` after `delta` in `p_types`, transforms, bounds, and the model
     documentation.
   - Default and exception must be exact zero; recommend fitting it shared or
     fixed before allowing racer/condition-specific values.
   - Pass `cv_u` by name through `dFRQ`, `pFRQ`, `sFRQ`, `rFRQ`, and `Ttransform`.
   - Report `lambda`, `cv_u`, and `a_u = 1/cv_u^2` (`Inf` at zero), while
     retaining the existing `p`, `N`, `d`, and `sQ` reports.
   - Document that integer `alpha`/`beta` interpretations are optional shape
     interpretations, not an implementation target for PCOUNTER equivalence.

2. `src/col_registry.h`
   - Preserve the current prefix order `alpha, beta, h, tau, t0, delta`.
   - Append `cv_u` and make it part of the required FRQ prefix used by the
     compiled likelihood. Update the local contract comment.

3. `src/model_FRQ.h`
   - Extend `FrqPars`, `frq_derive`, `FrqState`, and the cache key with `cv_u`.
   - Add one shared registration-state helper for `G_c`, `S_R`, and `log g_c`;
     do not duplicate the formula in the density and simulator.
   - Preserve the exact `cv_u == 0` branch and all existing `delta == 0`
     branches.
   - Keep the C++ checks for defective `p < 1`; permit `p == 1` only when
     `h == 1`.

4. `src/model_FRQ.cpp` and generated `R/RcppExports.R`
   - Extend direct `dfrq`, `pfrq`, and `frq_rate` entry points with a trailing
     `cv_u = 0` argument so existing positional calls keep their meaning.
   - Update scalar/raw/truncation callbacks to read the new column and include it
     in memoisation. Regenerate bindings with `Rcpp::compileAttributes()`; do
     not hand-edit the generated file.

5. `src/model_rng.cpp` and `R/model_rng.R`
   - Resolve `cv_u` by name, defaulting to zero only for the direct/reference
     simulator path where the column is absent.
   - Use the same shared inverse-registration helper as the likelihood.
   - Keep the existing direct `rbeta` draw when both extensions are inactive;
     use the `H^{-1}` path only when `delta > 0`.
   - Guard `u == 1`/`p == 1` endpoint roundoff so a proper model cannot acquire
     artificial omissions from a numerically infinite inverse.

6. `src/race_dispatch.cpp`
   - Keep the FRQ adapter and batched defective-truncation capability.
   - Update comments from “always defective” to “supports a defective atom, with
     a proper `h = 1` boundary”; no new dispatcher variant is needed.

## Validation sequence

### Kernel and mapping tests

- `cv_u = 0` is bit-identical to the current density, CDF, survivor, rate
  inversion, and simulator paths (including the existing `delta = 0` default).
- Independent R references for `G_c`, `g_c`, and the full FRQ density/CDF agree
  over small, moderate, large, and heavy-tail grids.
- Finite-difference density equals the CDF derivative; integrating the density
  gives `h` for defective cases and one for `h = 1`.
- `pfrq(tau)/h = 1/2` for all tested `delta`/`cv_u` combinations, and
  `pfrq(Inf) = h`, `sFRQ(Inf) = 1-h`.
- `cv_u -> 0` converges smoothly to the exponential member without relying on
  cancellation-prone `1/cv_u^2` arithmetic.

### Simulator and race-contract tests

- Simulate independent unit rates `R_i ~ Gamma(a, a/lambda)` and compare
  quantiles and omission frequency with the compiled FRQ simulator for integer
  and non-integer shapes.
- Exercise the combined `delta > 0`, `cv_u > 0` path, C++ and R fallback paths,
  loser survivors, omissions, censoring, and finite upper truncation.
- Add a proper-boundary fixture with `h = 1` fixed through the design/constants
  path. Assert `p = 1`, no simulated `Inf` responses, `s(+Inf) = 0`, and an
  omission likelihood of `-Inf` rather than a finite accidental mass.
- Pin the new column order and verify that a missing `cv_u` is accepted only by
  named direct/reference simulator inputs, never by the positional compiled
  likelihood contract.

### Recovery and identifiability gate

Run short simulations for: base FRQ, threshold variability only, unit frailty
only, both extensions, defective `h < 1`, and proper `h = 1`. Initially keep
`alpha`, `beta`, and `cv_u` shared across racers/conditions; let `tau` carry
the intended speed manipulation. Check posterior predictive q10/q50/q90,
omissions, and the `h`/`tau` anchors before permitting richer design formulas.

Do not add a trialwise common gamma rate (`sv`) in this model: with continuous
`alpha`, the exact PCOUNTER gamma-rate mixture loses its finite polynomial
closure. If common-mode rate variability is later required, make it a separate
quadrature or finite-mixture model with its own validation and race-dependence
contract.

## Correlation semantics and the time coordinate

`cv_u` is a marginal unit-rate heterogeneity parameter, not a correlation
parameter. It may sensibly be shared across racers and conditions through the
design (and should be shared initially), but that only says that the racers have
the same frailty distribution; it does not make their latent registration
processes or response times correlated.

The existing `PCOUNTERcorr` `rho` mechanism is specifically a shared gamma-rate
factor. It cannot be transplanted by relabelling `rho` as `cv_u`: doing so would
change the FRQ joint distribution and would invalidate the current
winner-density-times-loser-survivor factorisation. If correlated racers become a
requirement, add a separate `rho_u`-type shared-frailty model with an explicit
joint likelihood (probably one-dimensional quadrature or a finite mixture).
That is a later model family, not part of the first closed-form FRQ extension.

The time coordinate should remain `tau`. It is already the conditional median
completion time,

`P(T <= t0 + tau | T < Inf) = 1/2`,

not an approximation. For the frailty extension, `lambda` remains the mean
unit-registration rate and is solved so that this same median anchor holds. If
`r = u_tau / p`, the inversion is

`lambda = ((1-r)^(-cv_u^2) - 1) / (cv_u^2 * tau)`,

with the exponential expression recovered exactly at `cv_u = 0`. Thus `tau`
is the stable user-facing speed parameter, while `lambda` remains the
generative/reporting quantity. This is preferable to switching the fitted chart
to a free `lambda`, especially because finite frailty can produce heavy tails
and make means unstable or poorly identified. The main identifiability guard is
to keep `cv_u` shared/fixed initially and not fit it freely against `alpha`,
`beta`, `delta`, and `tau` in the same design.

## Explicit non-goals

- No literal integer PCOUNTER implementation, integer threshold sampler, or
  geometric threshold `omega` in FRQ.
- No `phi = gamma/(nu+gamma)`/`kappa` coordinate migration.
- No claim that independent unit frailty is equivalent to PCOUNTER's shared
  trialwise rate variability; it is an FRQ-native shape mechanism.
- No removal of the existing PCOUNTER model or its independent tests.

## Completion criteria

The extension is ready only when the base branch is unchanged, the combined
kernel/simulator passes the mapping and normalization tests, `h = 1` works
through the actual design and race likelihood, and the recovery study shows
that `cv_u` is not merely absorbing `alpha`, `beta`, or `delta` under the
intended default design.

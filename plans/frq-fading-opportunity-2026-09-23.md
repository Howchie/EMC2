# Plan: FRQ with fading evidence opportunity

Date: 2026-09-23

Status: implemented (2026-09-23).

## Goal

Add a sibling FRQ model for which omissions arise from the same evidence-registration process that generates response times. Preserve the existing `FRQ()` parameterization and behavior unchanged. The alternate model replaces fitted `h` with a finite-opportunity decay parameter `kappa`; it is not an exact reparameterization of current FRQ.

Retain the `delta` threshold-variability and `cv_u` unit-rate-frailty extensions, which a design can switch off by setting them to zero. Keep `pContaminant` as the existing optional external omission mixture for lapse-like omissions outside the evidence process.

## Model contract

Use a distinct constructor and C++ model identifier. `FRQfade()` / `FRQ_FADE` are provisional names; settle naming before implementation. The fitted parameters are:

| Parameter | Role |
|---|---|
| `alpha`, `beta` | Continuous quorum and spare-capacity shapes, with the same fitted support as FRQ: each at least 1. |
| `lambda` | Mean latent registration rate at decision time zero. |
| `kappa` | Exponential decay rate of evidence-registration opportunity; exact zero is a supported boundary. |
| `t0` | Non-decision time. |
| `cv_u` | CV of independent Gamma-distributed unit rates; exact zero recovers fixed unit rates. |
| `delta` | Half-width of the existing uniform log-odds shift in effective quorum percentile; exact zero recovers fixed threshold. |

`cv_u` and `delta` are optional in the design sense only: they default to zero and the model works with them switched off. In the C++ column contract they are **required** columns, as in the existing FRQ contract, because optional trailing columns resolve by position and a design that lacks one would silently read the wrong slot. Provisional kernel order: `alpha, beta, lambda, kappa, t0, delta, cv_u`.

Do not fit `h` in this variant. Report derived `h` and `q_inf` in mapped output. Unlike current FRQ, `delta` will generally change derived `h`; do not compensate with an availability inversion, since that would reintroduce a separate completion coordinate. Document that `cv_u` also affects `h` whenever `kappa > 0`.

## Mathematical kernel

For decision time `u = rt - t0`, define

\[
A_\kappa(u)=\begin{cases}
  (1-e^{-\kappa u})/\kappa,&\kappa>0,\\
  u,&\kappa=0.
\end{cases}
\]

For `cv_u > 0`, let `c2 = cv_u^2` and

\[
q(u)=1-[1+\lambda c2 A_\kappa(u)]^{-1/c2},\qquad
q'(u)=\lambda e^{-\kappa u}[1+\lambda c2 A_\kappa(u)]^{-1/c2-1}.
\]

At `cv_u = 0`, use the limits

\[
q(u)=1-e^{-\lambda A_\kappa(u)},\qquad
q'(u)=\lambda e^{-\kappa u}e^{-\lambda A_\kappa(u)}.
\]

With `delta = 0`, the CDF is `I_q(alpha, beta)`. With `delta > 0`, compose the existing FRQ threshold map:

\[
F(u)=H_\delta(I_{q(u)}(\alpha,\beta)),
\]

\[
f(u)=H_\delta'(I_q)\,
\frac{q^{\alpha-1}(1-q)^{\beta-1}}{B(\alpha,\beta)}q'(u).
\]

Use the reflected beta identity and symmetry of `H_delta` for the survivor:

\[
S(u)=H_\delta(I_{1-q(u)}(\beta,\alpha)).
\]

For `kappa = 0`, `q_inf = 1` and `h = 1`, including when `delta` and `cv_u` are active. For `kappa > 0`,

\[
q_\infty=1-[1+\lambda c_u^2/\kappa]^{-1/c_u^2},\qquad
h=H_\delta(I_{q_\infty}(\alpha,\beta)),
\]

with `q_inf = 1 - exp(-lambda/kappa)` at `cv_u = 0`. The survivor at infinity must be `1 - h`.

Equivalent representation (useful for the interpretation and the simulator): with `cv_u = 0`, `q(u)` is the proper (`h = 1`) FRQ curve evaluated on the internal clock `s = A_kappa(u)`, which stops at `s = 1/kappa`. A response occurs if and only if the undecayed quorum time `T*` falls below `1/kappa`, and the observed decision time is then `u = -log1p(-kappa T*)/kappa`. For fixed shapes, `q_inf` depends only on `m = lambda/kappa` and `kappa` is a pure time scale, so the conditional RT law is a scale family whose shape is set by the completion probability. This is where the omission-RT coupling comes from, and it does not vanish at `alpha = beta = 1` (unlike current FRQ, where the conditional RT law there is exactly `h`-free).

Implement cancellation-safe branches: use `expm1` for `A_kappa`, `log1p`/`expm1` for the frailty survival and `q`, explicit zero-parameter limits, and a dedicated infinite-time endpoint. Avoid computing complements by subtracting nearly equal probabilities where the existing reflected FRQ path can be reused.

## Implementation steps

1. **Add the public model variant.** Put its constructor in a separate R model file or a clearly isolated section of `R/model_FRQ.R`. Define its own canonical parameter order, defaults, bounds, exceptions, mapped output, `dfun`, `pfun`, simulator hook, and optional nuisance parameters. Keep `FRQ()` unchanged.

2. **Add a distinct C++ column contract and dispatch.** Register a required positional contract for the alternate parameters in `src/col_registry.h`. Add an explicit `FRQ_FADE` dispatch branch before the existing substring match for `FRQ` in `src/race_dispatch.cpp`; otherwise the alternate name would select the current `alpha, beta, h, ...` contract. Update the comment on the existing `FRQ` branch, which currently claims that no other c_name contains `"FRQ"` and that the branch position is therefore free. (The global prefix/suffix probes — `_E2`, `_EMIX`, `_GLOBAL`, `_IO`, etc. — were checked and none match `FRQ_FADE`.) Mark the adapter as defective-tail capable and provide its `logS_at_t` callback so batched finite-UT truncation can use the analytic endpoint.

3. **Add the rate/opportunity state.** Implement a separate state builder for `q`, `1-q`, and `log(q')`, with `lambda`, `kappa`, and `cv_u` validation and stable limits. The current `FrqPars` availability half derives `p` from fitted `h` through an iterative `qbeta` solve cached in `FrqAvail`/`FrqMemo`; none of that applies here. The variant has no availability inversion at all (so it should be cheaper per call than current FRQ), and only `lbeta(alpha, beta)` is worth caching. Give it its own small cache rather than reusing `FrqCache`.

4. **Share the quorum and threshold calculations.** Reuse the incomplete-beta evaluator and `H_delta` value, inverse, and derivative routines from `src/model_FRQ.h`. Structure the alternate evaluator around its own state and the shared beta-quorum operations. Keep the existing FRQ path numerically and behaviorally stable; any common evaluator refactor must preserve its outputs.

5. **Implement density, CDF, survivor, and simulator entry points.** The density uses `log(q') - lbeta(alpha,beta)` rather than current FRQ's `log(p) + log(G') - lbeta`. For simulation, reuse the latent quorum-percentile construction (including the `delta` transform); the inverse of `q(u)` is closed form: `A = -log1p(-U)/lambda` at `cv_u = 0`, or `A = expm1(-c2 log1p(-U))/(lambda c2)` with frailty; the draw is an omission (`+Inf`) when `kappa A >= 1`, otherwise `u = -log1p(-kappa A)/kappa` (and `u = A` at `kappa = 0`). Decide omission with the `kappa A >= 1` test in the same arithmetic used for `u`, not by comparing `U` with a separately computed `q_inf`, so borderline draws cannot produce NaN or negative times. Add the compiled C++ simulator as the default path and retain an R reference simulator if required by the model conventions.

6. **Document the interpretation and mappings.** State that `kappa` limits total integrated opportunity, rather than imposing a hard finite cutoff in physical time; equivalently, a hard horizon `1/kappa` on the internal clock `A_kappa(u)` (see the representation above). Explain that `h` is derived, `delta` can alter it, and `cv_u` affects it under finite opportunity. Report `h`, `q_inf`, and existing FRQ-style shape summaries as mapped quantities where meaningful.

## Acceptance checks for implementation

- At `kappa = 0`, density, CDF, survivor, and simulation agree with current FRQ fixed at `h = 1` under matching `alpha`, `beta`, `lambda`, `cv_u`, `delta`, and `t0`.
- At positive `kappa`, `F(Inf) = h`, `S(Inf) = 1-h`, and the density integrates to `h`.
- Check the `cv_u = 0` and `delta = 0` branches, exact zero `kappa`, small positive `kappa`, and valid parameter-bound behavior.
- Simulated completion frequency and finite response times agree with the analytic CDF, including the `delta` branch.
- Race likelihoods preserve the model-implied omission probability and work with the existing defective-tail/truncation contract.
- Existing `FRQ()` density, CDF, survivor, `logS_at_t`, and seeded simulation outputs are bit-identical before and after the change (guards the dispatch reordering and any shared-evaluator refactor).
- The omission-RT coupling the variant exists for is actually present and recoverable: a validation-level check (not a package test) that the conditional RT law changes with completion probability at fixed shapes, including `alpha = beta = 1`, and a parameter-recovery run on a multi-condition design where omission rate and RT move together. Reference magnitudes from the review (`cv_u = delta = 0`, `kappa` refit freely): KL of the conditional RT law for `h` 0.95 -> 0.60 is about 1.1e-2 nats/response at `(1,1)` and 9.7e-3 at `(2,3)`, and about 7.7e-3 even with `alpha, beta` free in the second condition.

The repository requires a compiled default simulator and a regression test for a new model. These checks are implementation requirements; none were added or run while writing this plan.

## Scope limits

- Do not change or remove current `FRQ()` or reinterpret its fitted `h`.
- Do not fit both `h` and `kappa` in this variant unless a separate persistent-availability mechanism is explicitly introduced.
- Do not change PCOUNTER or combine fading opportunity with PCOUNTER self-excitation as part of this work.

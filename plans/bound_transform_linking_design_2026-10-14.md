# Bound/transform linking: design note

Companion to `bound_transform_linking_2026-10-14.md`. Implementation branch:
`bound-transform-link` (worktree `/data/work/EMC2_bound_link`).

## The linking rule

Applied in `fill_transform`, which now receives the filled bound. Per model
parameter (`p_types` name):

1. **Explicit sides win.** A `transform$lower`/`upper` declared by the model
   constructor or the user is used as given. This keeps FRQ/FRQfade
   (`alpha`/`beta`, lower 1, already equal to the bound) and RLF `alpha`
   (explicit `pnorm` window [1, 2], guard bound [1.01, 1.99], exception 2)
   exactly as they are.
2. **Unspecified sides are filled from the bound.** The link's default side
   (`exp` lower 0, `pnorm` lower 0 / upper 1) moves to the bound when the
   bound is tighter: `lower = max(default, bound_lo)`,
   `upper = min(default, bound_hi)`. `identity`, and the `exp` upper side,
   cannot absorb anything. The family is never changed.
3. **Parameters with a constant keep the plain link.** If any design column of
   a parameter is fixed by `constants` (including parameters left out of the
   formula, which become constants at their `p_types` default), the default
   sides are used. Constants keep their natural meaning (`sv = log(1)` is
   still 1, `t0 = log(0)` is still 0), and every exception value (`t0 = 0`,
   `A = 0`, `pContaminant = 0`, `sv = 0`, `pi = 0/1`, ...) is still reached
   exactly. The cost is that such a parameter keeps its strict-test wall on
   its sampled coefficients.
4. **A sampled ±Inf keeps its plain meaning.** A linked side records its
   default in `transform$lower_inf` / `upper_inf`, and the C++ transform
   maps an exact −Inf / +Inf there (0, or 1 for a `pnorm` upper side). So
   `log(0)` and `qnorm(0)` in a sampled vector (for example `pGuess =
   qnorm(0)` in `make_data`, or a test switching a mechanism off) still give
   0 and reach the exception. A finite value that rounds onto the linked side
   (`0.05 + exp(-50)`) is not the exception and is still rejected by the strict
   test. Samplers never propose ±Inf, so fits are unaffected.
5. **Bound overrides move the link.** `design(bound = ...)` is merged by
   `fill_bound` before the transform is filled, so the rule applies to the
   user's bound. A user transform on the same side still wins (rule 1).
6. **Validation.** It is an error if a parameter's bound and its transform
   window do not overlap (nothing could ever be in bounds). A transform window
   narrower than the bound is allowed. A wider one is allowed too: the strict
   test keeps guarding that side, as now.

## What the strict test still does

`do_bound` / `c_do_bound` are unchanged and remain the backstop for:

- sides the family cannot absorb (`exp` upper guards such as sv ≤ 10, a ≤ 10;
  `identity` guards such as v in ±20);
- parameters on the plain link because they have a constant;
- natural-scale paths after the transform (`posttransform` trends,
  `Ttransform`). `premap` and `pretransform` trends act before the transform,
  so the link covers them;
- floating-point saturation on an absorbed side: `0.05 + exp(x)` rounds to
  0.05 below x ≈ −39.5, and `pnorm` saturates past |x| ≈ 8.3. Under any
  sensible group distribution these sit many sds out, so they are left
  to the strict test rather than relaxed.

Data-dependent limits (t0 above the fastest RT) are kernel zeros, not bound
failures. They are unchanged.

## Downstream code that already follows the transform

- C++ transform specs (`src/transform_utils.cpp`) read `lower`/`upper` per
  parameter; no change.
- Marginalised t0 (`src/marginal_likelihood.cpp`) inverts the transform at the
  natural bound, so `log(0.05 - 0.05) = -Inf` removes the edge automatically.
- Shear/collapsed walls (`.emc_shear_support`) already skip a bound at the
  transform's own floor; the t0 wall simply disappears. No edit to
  `R/shear_move.R`.

## Old objects

`design()` stores the filled transform in the model closure, and nothing
re-derives it on load. `update2version` passes the old filled transform back
as an explicit transform, so rule 1 keeps it. Saved fits, priors and designs
map exactly as before.

## User-facing scale (phase 2, after phase 1 lands)

Phase 1 changes what sampled-scale numbers mean for linked parameters (only
t0-class bounds and a few others move noticeably; 1e-4 floors are
invisible). Phase 1 also prints each linked parameter's link in
`summary(design)` (for example `t0 = 0.05 + exp(x)`), so the sampled scale is
labelled where users first see it. Phase 2 options:

- natural-scale cell values (`mapped_pars`, `get_pars(map = TRUE)`) as the
  default in summaries and credints, with the sampled scale opt-in and labelled;
- prior and `make_data` location given on the natural scale through a helper
  that inverts the link (spreads stay on the sampled scale);
- default priors whose natural-scale median matches today's.

`make_random_effects` reads the `exp` cell spread as a lognormal CV of the
natural value. For a linked parameter it is a CV of (value − lower). That is
documented rather than changed.

## Risk found during implementation: the boundary tail

A bijection from the real line onto (L, ∞) puts L at x = −∞. Where the
likelihood is flat in natural t0 near the bound (the data want t0 at or below
0.05), the sampled-scale posterior of that subject has a left tail held up
only by the group distribution. The group variance of log(t0 − 0.05) then
inflates to cover it, which weakens shrinkage for every subject. The old
model's cliff truncated the tail at x = log(0.05).

This was seen in a 2-subject, 25-iteration LNR preburn (`test-variant_funs`):
one subject went to t0 = 0.050 and theta_var(t0) rose to 30-100 (sampled
scale), against 0.01-1 on the plain link. Forstmann's subjects sit well above
the floor (lowest 2.5% quantile 0.116), so it is not affected. The N-back fit
(17 subjects on the floor) is the deciding check for whether the link is the
right model, before merging.

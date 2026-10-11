# Binned-rt likelihood: open issues after the af8d830e review

Two accuracy issues found while checking the binned response-time likelihood
(`src/rt_bins.h`) on 2026-10-08. Neither is fixed. The two bugs the review did
fix (DDM bins built from CDF differences, and the SS raw path's separate
integrator) are in the "Binned rt" commit.

## 1. Midpoint rule near onset

**What happens.** A trial's bin probability is integrated exactly only when
the onset (t0, or the earliest onset in the trial) lies inside the bin or
less than one bin width below it. `rt_bin_quadrature_weight` blends from the
integral to the midpoint over that one bin with a C1 smoothstep. Further
above onset the bin is scored as `h * f(midpoint)`.

**Why it matters.** The midpoint rule's error grows with the density's
curvature over the bin, and the leading edge of an RT density is where the
curvature is largest. On SSRDEX (B = 1, rt_resolution = 1/60), fast
outlier trials one to two bins above t0, with log densities between -10 and
-18, scored 1-3 nats away from a 1e-10-tolerance reference integral. The
error is per trial and lands on the trials that inform t0 most, so it can
move t0 and the variability parameters. On a convex leading edge the
midpoint underestimates the bin. I measured the size of the error but did
not tabulate its sign.

**What is arbitrary.** The one-bin blend window is a fixed constant. It is
too short for steep leading edges and longer than needed in the body of the
distribution.

**Proposed fix.** Make the switch depend on an error estimate:

- Per trial, estimate the midpoint error from values the rule already
  computes or can get cheaply. For example, compare the midpoint with a
  two-point Gauss rule in z = sqrt(t - onset).
- Set the quadrature weight to a smooth function of log(error / tolerance),
  instead of distance to onset in bins.
- The error estimate is continuous in the parameters, so the likelihood stays
  continuous. That property is what the current smoothstep exists to
  guarantee.
- Cost is limited to the trials whose bins actually need quadrature,
  wherever they are relative to onset.

## 2. DDM density jumps with st0 > 0

**Status (2026-10-11): fixed, uncommitted.** `src/ddm_variability.h` replaces
`ddiff` with fixed-node rules (closed-form start-point integral at small
u/a^2, fixed rule over the large-time series otherwise, log-sqrt map over the
st0 window). st0 > 0 bins use the midpoint window plus a per-side curvature
correction, with a side integrated exactly where it lies on the leading edge
(within ~sqrt(h) d_min of u = 0, at onset and at t0 + st0) or where a bin spans
more than the density's decay length. The plain midpoint biased st0 by up to
~3.6 SE on simulated 600-trial data sets (bias from bins near t0 + st0, not
onset); the corrected rule keeps ML shifts below 0.1 SE at ~1.4x the shipped
binned cost (exact-everywhere was ~2.6-2.9x). The CDF (`pdiff`) moved to the
same fixed rules. The absolute tolerance also made the shipped density
inaccurate (q99 0.4 nats), not only discontinuous.

**What happens.** With between-trial variability (st0, and the same code
path for sv and SZ), the DDM log density has small discontinuities in the
parameters. Example: v = 1, log a ≈ -0.965, sv = 0.5, SZ = 0.1, st0 = 0.05,
rt - t0 = one bin (1/60 s). At one value of `a` the log density jumps by
5.7e-3 nats; the median step between neighbouring grid points is 3e-5.
With st0 = 0.02 the jump is 4.9e-3 at a nearby `a`. The jump appears with
and without rt binning, so it predates af8d830e.

**Cause (confirmed by varying the tolerance).** `ddiff`/`pdiff` integrate
over the variability distributions with `hcubature` to an absolute
tolerance of 0.9 eps on the density, with `eps = 5e-3` and at most 6000
evaluations (`src/model_DDM.h`, `src/fncs_seven.cpp`).

- The jump is the cubature accepting its first-pass estimate. Over a band of
  `a` that estimate is about 5e-3 nats off, but its error estimate collapses
  toward zero, so it passes the tolerance.
- On the right of the jump the shipped values equal the first pass to seven
  digits (cap the evaluations at 1 or 17 and you get the same numbers).
- Raising the evaluation cap changes nothing. Tightening `eps` only narrows
  the band where the bad estimate passes. On a 2,501-point grid over log a
  in [-0.99, -0.94], the points more than 1e-3 off a tight reference
  (eps 1e-8) were:

  | eps  | 5e-3 | 5e-4 | 5e-5 | 5e-6 | 5e-7 |
  |------|------|------|------|------|------|
  | points | 833 | 83 | 8 | 1 | 0 |

  Every accepted point is still off by the full ~5e-3.
- Per-evaluation cost: 0.11 ms at 5e-3, 0.13 ms at 5e-6, 0.47 ms at the
  1e-8 reference.

The same integrator and tolerance give the Wiener CDF, which is why that
CDF is accurate only to about 1e-3 near onset (the reason the binned DDM now
integrates the density instead).

**Why it matters.**

- Each jump is small, but every affected trial carries one, at its own
  parameter values.
- Finite-difference gradients and Hessians with steps near 1e-4 (the
  mode-start Hessian, Laplace fits used by the shear and collapsed moves) see
  slopes of several hundred at a jump.
- Subject HMC sees a discontinuous energy surface there.

**Proposed fix.**

- Integrate over st0 (uniform) and SZ (uniform) with fixed-node
  Gauss-Legendre rules. Their nodes do not move with the parameters, so the
  result is smooth.
- Choose the node count once by accuracy against a tight reference rather
  than adaptively per call.
- If the adaptive integrator stays: never accept the first pass. Force at
  least one subdivision so the error estimate compares two refinements
  rather than two rules on one region.
- Tightening `eps` alone shrinks the region where the bad estimate passes,
  but every point that still passes keeps the full error. It also makes
  where the jumps fall depend on the tolerance.

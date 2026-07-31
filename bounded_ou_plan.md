# Bounded OU: the DDM with leak, on the Fokker-Planck solver

## What this is

Smith & Ratcliff (2004, *Psychol Rev* 111:333-367, Appendix pp. 41-45) describe a
two-choice OU diffusion:

```
dX = ( xi - beta*(X - z) ) dt + s dW,    X(0) = z,    absorbing at 0 and a
```

with across-trial variability in drift (eta), start point (s_z) and non-decision
time (s_t). One process, two absorbing barriers; which barrier is hit encodes the
response. It is the DDM with a leak term, and at `beta = 0` it *is* the Wiener
diffusion.

The paper computes it by the Buonocore/Fortet integral-equation recursion. We
compute it with the existing Fokker-Planck solver (`src/fpe_solver.h`) instead, by
adding a second absorbing boundary. That solver previously supported exactly one
absorbing boundary (at `xi = 1`), with `xi = 0` a no-flux far-field face placed
6 diffusive sd away by `fpe_x_lo_*`.

**Design rule:** at `beta = 0` the model must reproduce the package's own
Navarro-Fuss DDM (`dDDM`/`pDDM`) exactly, while travelling a completely different
code path (PDE march vs series expansion). That makes the existing DDM a free and
decisive oracle. It is the first test in `tests/testthat/test-fpe-bounded-ou.R`.

## Conventions

Per-response quantities follow the package's DDM (`R/model_DDM.R:70,83`) exactly:

- `dfun(rt, R, pars)` / `pfun(rt, R, pars)` return the **defective** density and
  cdf for that response. They sum to 1 across responses, not each on their own.
- `lower` is the first level of `R`, `upper` the second.
- `FPE_Result::pdf`/`cdf` are the **upper** barrier; `pdf_lower`/`cdf_lower` the
  lower. `surv` is shared (total unabsorbed mass) and keeps its existing meaning.

## Decisions taken

**Decay anchor.** S&R's drift decays toward the *starting point* `z`, and their
footnote 2 argues for that over the alternative. Note that "decay toward zero" is
not available in this parameterisation at all -- zero is one of the response
boundaries -- so the only other sensible anchor is the midpoint `a/2`.
`FPE_ModelBoundedOU::anchor` is a field defaulting to `z`, which costs one
addition in `atil_affine` and leaves the choice open. It also matters for cost:
because the anchor is `z`, **`z` enters the operator and not just the initial
condition**, so across-trial start-point variability is a genuine second
quadrature dimension (`n_eta x n_z` solves per condition). An anchor independent
of `z` would collapse the whole uniform start into one solve.

**Likelihood path.** The bounded OU is DDM-shaped, not race-shaped. Rather than
reimplement truncation/censoring, generalise `c_log_likelihood_DDM_pt`
(`particle_ll.cpp:1927`) behind a context object holding the two primitives it
actually needs -- `d_DDM_Wien_raw` and `p_DDM_Wien_raw` (`model_DDM.h:53,133`)
share the signature `(rts, Rs, cols, n_rows, mask, is_ok, out, min_ll)`. Default
to the Wiener pair so DDM stays byte-identical.

---

## Stage 1 -- solver core (DONE, validated)

`src/fpe_solver.h`, `src/fpe_models.h`, `src/fpe_diffusion.cpp`,
`tests/testthat/test-fpe-bounded-ou.R`.

- `FPE_Mesh::build(M, grade, symmetric)`: symmetric `tanh` stretch clustering at
  both ends, `c = acosh(sqrt(grade))` so `grade` keeps meaning "widest cell over
  boundary cell" for both maps. `cL_A`/`cL_B` at `xi = 0` from the same one-sided
  quadratic fit that produces `cA`/`cB`, so on a symmetric mesh they coincide.
- `build_op`: `if constexpr (Model::lower_absorbing)` on the `i == 0` face.
  Signs mirror the `xi = 1` row: the outward flux is
  `J_0 = D*(cL_A*q[0] - cL_B*q[1]) > 0`, the signed face flux is `F_0 = -J_0`.
- `flux_out_lower()` alongside `flux_out()`, with the same negative-value floor
  (the one-sided difference goes to ~-3e-16 in the far tail and `log()` of that
  would poison a likelihood).
- `FPE_Result` gains `pdf_lower`, `cdf_lower`, filled only under `if constexpr`.
- `fpe_solve` gains a `cdf0_lower` argument and splits the absorbed mass:
  `cdf_lower` from the flux trapezoid, held monotone and inside `[prev, cdf_mass]`;
  `cdf_upper = cdf_mass - cdf_lower`. **The exact mass identity anchors the pair,
  so only the split carries quadrature error** -- and because `cdf_mass` never
  decreases and previously equalled the sum, both defective cdfs come out
  monotone with no third clamp that could break the identity.
  `flux_mass_mismatch` now compares against the SUM of both flux integrals.
- `fpe_seed` gains the image reflected about the lower barrier, a `t_seed` cap of
  `0.25*min(a - z, z - x_lo)`, and an `absorbed_lower` out-parameter. That last
  is not a nicety: the cap leaves each barrier 4 sd away, so ~3e-5 of the mass is
  absorbed before `t_seed` and would otherwise be charged to the wrong response.
  Both barriers being >= 4 sd away is also what licenses truncating the infinite
  image series after this pair (the next images sit >= 8 sd out, < 1e-14).
- `FPE_ModelBoundedOU`: `drift(x) = v + beta*(anchor - x)`, domain exactly
  `[0, a]` with no `fpe_x_lo_*` padding, `a0 = (v + beta*(anchor - xlo))/L`,
  `a1 = (-beta*L - L')/L`. At `beta = 0` these are `FPE_ModelBM`'s exactly.
- `fpe_bou_fht_pdf_cdf_vec()` validation entry point in `fpe_diffusion.cpp`.

### Measured, and it contradicted the original plan twice

**Grading never wins here -- the default is uniform.** The plan assumed the fix
was symmetric grading. Measured (`beta = 0` vs the DDM oracle, nx=384, nt=3000,
worst `|dlog f|` over the central mass):

| grade | 1 | 4 | 8 |
|---|---|---|---|
| sigma=1.00 | 1.22e-4 | 1.24e-4 | 1.91e-4 |
| sigma=0.50 | 5.52e-3 | 8.18e-3 | 1.08e-2 |
| sigma=0.25 | 1.69e-1 | 2.36e-1 | 3.08e-1 |

Uniform wins everywhere and the gap widens as the domain gets wide in units of
sigma. There is no far field to economise on: the domain is exactly `[0,a]` and
the density is O(1) across it, so graded cells are taken from where the solution
varies. Hence `FPE_GRADE_BOUNDED = 1.0`. The symmetric map is still worth having
-- it is what makes the symmetry test pass at `grade = 8` -- but it is not the
default.

**The error is time-dominated, not space-dominated.** Clean O(dt^2) from
Crank-Nicolson; space plateaus by nx ~ 1024:

```
nx    128     256     512     1024    2048       nt    1500    3000    6000    12000
err  1.28e-2 4.32e-3 1.68e-3 9.16e-4 7.17e-4    err  1.07e-2 2.70e-3 7.17e-4 2.19e-4
```

So the accuracy budget belongs on `nt`/`tgrade`, inverting the ROU grid's
priorities (`fpe_race.h:116`: "nx is the binding constraint and dt is nearly
free"). Stage 2 must not simply inherit `FPE_Grid`'s defaults.

**Measurement caveat.** An early sweep reported 9.5e-2 log-error; it came from
screening on `ref > 1e-8`, which picks up the extreme rising flank (t=0.07,
density 2.9e-8, 3e-10 of the mass below it). On the central 99% of mass the same
solve is 1.7e-3. The tests use a mass-based window for densities and unfiltered
comparison for cdfs.

### Validation status

129 assertions passing; no regression (`test-fpe-fht.R` 39/39, `test-rou.R` 62/62).

| Test | Result |
|---|---|
| `beta = 0` vs Navarro-Fuss DDM, 24 param combos, both boundaries | cdf < 2e-3, log-pdf < 2e-2 |
| Symmetry (`v=0, Z=0.5`), `beta` in {0,2,4,8}, `grade` in {1,8} | 1e-13 |
| `cdf_up + cdf_lo + surv - 1` | 2.2e-16 |
| Leak actually changes the model (no-op guard) | P(upper) 0.72 -> 0.46 |

The symmetry test is the load-bearing one: at `grade = 8` it is what confirms the
symmetric mesh resolves both boundaries alike, and it is what would have failed
loudly under the original one-sided grading.

---

## Stage 2 -- production path (DONE)

`src/fpe_bou.h` + `src/bou_diffusion.cpp`.

Rather than template the race lane kernel (`fpe_race.h:597`), which carries
machinery this model has no use for -- loser survivors, sparse queries, per-lane
schedules and lane retirement -- the bounded model got its own cache and
lane-batched march. Its parallelism is across QUADRATURE NODES, which share a
mesh, a horizon and a clock exactly, so a common-clock batch is not an
approximation here (it is precisely the objection that makes
`fpe_solve_batch_ou_common_clock` deprecated for racing accumulators).

- `Key` is bit-exact and quadrature nodes are folded in **before** the key is
  built, so a node is just another key and the cache dedupes across nodes free.
- `SolveCache::t_horizon` is set ONCE per batch from the largest decision time.
  Growing it per row re-solved every node as a sorted RT vector was walked; that
  is what made the first version unusably slow, and it is now pinned by a test
  asserting the solve count is exactly `n_sv * n_sz`.
- `fpe_seed` gained `t_seed_force` so every lane starts from one clock: the
  batch adopts the smallest `t_seed` any lane would have chosen, so no lane is
  seeded later, and therefore less accurately, than it would have been alone.
- The layout is lane-interleaved and left to auto-vectorisation rather than
  hand-written intrinsics -- every lane runs an identical instruction sequence,
  which is the case the compiler handles well, for a fraction of the code.
  Measured (118 rows, sv+SZ+st0, 49 nodes, nx=512):
  scalar 0.855 s, 4 lanes 0.458 s, 8 lanes 0.365 s.
- Defaults are NOT `fperace::FPE_Grid`'s: `dt_target = 5e-4` with
  `grade = FPE_GRADE_BOUNDED` (uniform), because this model is time-dominated.

## Stage 3 -- across-trial variability (DONE)

Gauss-Hermite over `sv`, Gauss-Legendre over `SZ`, both folded into the cache
key; `st0` is an outer Gauss-Legendre loop over the query time, so it triggers
no solves at all. Cost is `n_sv * n_sz` marches per parameter set.

Two conventions had to be established by measurement, not assumed, and both are
now pinned by tests:

- **`st0` is `U(t0, t0+st0)`** -- the lower-edge form. Reconstructing an
  `st0 > 0` density from point-`t0` densities, the lower-edge average matches to
  2e-3 where the centred form is off by 3.4e-1.
- **`SZ` arrives RAW** at the C++ entry points and is widened internally by
  `2*SZ*min(Z,1-Z)`, exactly as `d_DDM_Wien_raw` does (`model_DDM.h:122`).
  Ttransform is an R-side step and does NOT run on this path. The R-level
  `dDDM`/`pDDM` are the opposite contract -- they take the already-widened value,
  because by the time they are called Ttransform has run. Comparing the two
  without accounting for this shows a spurious 0.6% disagreement at Z = 0.45.

Validated against the DDM oracle at `beta = 0`: density and cdf agree to ~3e-5
on both responses for every combination of `sv`, `SZ` and `st0`.

## Stage 4 -- R wrappers and simulator (DONE)

`R/model_BOU.R` -- wrappers only, no R-side likelihood (everything is C++).
DDM's `p_types`, `transform`, `bound` and `Ttransform` reused verbatim plus
`beta`, so a DDM design converts by adding `beta~1`. `log_likelihood` is the
shared `log_likelihood_ddm`.

The simulator (`rbou_cpp`) steps the exact OU transition -- `beta = 0` regular
through `expm1`, not a special case -- and applies a Brownian-bridge correction
at BOTH barriers for paths that cross and return within a step. Without it RTs
are biased up and error rates down, which would make it useless as the
independent check it exists to be. It agrees with the solver on P(upper) within
Monte Carlo error and on quantiles to 1-4 ms.

## Stage 5 -- collapsing bounds (DONE)

Both barriers close **symmetrically on the midpoint** `a/2`, so the separation
runs from `a` at `t = 0` down to `aInf`. That is the only symmetric choice
available: there is no fixed level to collapse toward, because zero and `a` are
the two response boundaries. It also preserves the invariance under the response
relabelling `(v, Z) -> (-v, 1-Z)` that the symmetry test checks.

The internal form is `(midpoint, half-width)`: `FPE_ModelBoundedOU` holds `mid`
and hands the existing `FPE_Boundary` module the HALF-separation `h(t)`, with
barriers at `mid +/- h(t)`. All four collapse forms (Weibull, exponential, the
two linear ones) therefore come for free and behave exactly as they do for ROU.
`set_separation(a, kind, aInf, tau, pw)` is the single place the user-facing
`(a, aInf)` parameterisation meets the internal one.

`0.5*a` is exact in binary and `2*(0.5*a) == a`, so a fixed bound reproduces the
old `(0, a)` geometry bit for bit -- which is what lets the degenerate case
(`aInf == a`, where `set_kind` reports `fixed = true`) be asserted with
`expect_identical` against the Navarro-Fuss-validated fixed-bound model.

### The frame term

`x_lo` was promoted from a constant to `x_lo(t)` across the model concept (all
three structs; BM/OU ignore the argument). With `x = x_lo(t) + xi*L(t)` the
transformed drift picks up `-(x_lo' + xi*L')/L`, so a moving lower edge adds
`-x_lo'/L` to `a0`. Here `x_lo' = -h' = -L'/2` exactly, so it enters as
`+L'/2` -- an outward drift on the lower side, which is the correct sign:
mass standing still in the lab frame is moving toward `xi = 0` relative to a
rising floor. There is deliberately no `x_lo_prime()` accessor; the quantity is
already exact from the `Lp` the operator builder passes in, and a second route
to it would be a second thing to keep in step.

`fpe_seed` reads `x_hi(t)`/`x_lo(t)` rather than `bnd.a(t)`/`xlo`, and caps
`t_seed` at `0.1*tau` so the two-sided gap it measures at `t = 0` is still valid
at the time it seeds. In practice that never binds -- `t_seed` is ~1 ms.

### The separation floor

`h(t)` is floored at `FPE_BOU_MIN_SEP = 2%` of `a`, so `aInf = 0` means "collapse
essentially all the way" rather than exactly so. Nothing in the SPACE
discretisation needs this -- the mesh lives on the normalised `xi` and follows
the barriers down -- but `D(t) = (B/L)^2/2` grows like `1/L^2`, so the time step
stops resolving the flux long before the width does. At 2% of `a` the remaining
survivor crosses in `(0.02a/s)^2` seconds, so the floor moves where that last
sliver is absorbed by well under a millisecond.

### A latent bug the kink exposed

`fpe_solve` used to clamp the lower defective cdf to `[cdf_lo_prev, cdf_mass]`
and hand the upper boundary the remainder, on the argument that the remainder is
then automatically monotone. **That argument is wrong**, and the comment
asserting it has been corrected. `cdf_mass` being non-decreasing bounds the
lower cdf's LEVEL but not its INCREMENT, so a lower flux trapezoid that outruns
the mass route over a step pushes the remainder below `cdf_up_prev`. It never
fired while every bound was fixed; it fires on the kinked `linear_additive`
collapse, where `b'(t)` jumps at `t = tau`.

The window is now two-sided, `[cdf_lo_prev, cdf_mass - cdf_up_prev]`, which is
never empty (`cdf_mass >= cdf_lo_prev + cdf_up_prev` because `cdf_mass` does not
decrease) and preserves the exact-sum identity. Fixed in both `fpe_solve` and
the lane kernel.

### Cost

`static_op()` goes false, which forfeits the one-time factorisation: the operator
is rebuilt and refactorised every step instead of once per schedule block. To
keep that as cheap as possible the lane kernel got `bou_build_op_lanes`, which
writes the tridiagonal operator for all lanes STRAIGHT into the interleaved
layout -- same stencil, same Scharfetter-Gummel weights, same uniform-mesh
Peclet recurrence as `fpe::build_op`, with the lane loop innermost, so no
per-lane `FPE_Op` is materialised and scattered. The batch carries a second
operator buffer because Crank-Nicolson needs the old-time operator for the
right-hand side while the left-hand side needs the new one.

Measured (200 query times, `sv`+`SZ`+`st0` at 7x7 nodes, nx=384):

    fixed bounds       0.204 s
    exponential        0.502 s     <- 2.5x

That is the whole price of a collapse: a per-step build plus a per-step
factorisation on top of the two Thomas sweeps that were already there.

Convergence is unchanged in character -- still time-dominated, not
space-dominated (relative pdf error against a refined reference):

    nx=192  nt=1500  9.8e-4      nx=384  nt=1500  9.9e-4
    nx=192  nt=3000  2.8e-4      nx=384  nt=3000  2.0e-4
    nx=192  nt=6000  2.9e-4      nx=384  nt=6000  5.2e-5

### Validation

No closed form exists for a diffusion between two converging boundaries, so the
evidence is of three independent kinds:

1. **Degenerate reduction.** `aInf == a` is `expect_identical` to the
   fixed-bound model, for all four forms.
2. **Two independent implementations of the moving march.** `fpe_solve` rebuilds
   an `FPE_Op` and refactorises an `FPE_Tri`; the batch runs the interleaved
   builder and its own Thomas elimination. With fixed bounds both build the
   operator once, so this comparison only becomes load-bearing for the collapse.
3. **The simulator**, which steps the exact OU transition against linearly
   interpolated barriers (the Brownian-bridge probability against a linear
   boundary is the same expression with the barrier taken at each end of the
   step) and shares none of the PDE's numerics. Used as a check, never as a
   calibration target.

Plus the mass identity, monotone defective cdfs, and symmetry under collapse.

Measured: degenerate reduction exact to 0; mass identity 2e-16 and symmetry
~1e-11 for all four forms; lane-vs-scalar 5e-5 for the three smooth forms and
1.9e-3 for the kinked `linear_additive` (the kink costs an order locally, as
`FPE_Boundary`'s own comment predicts); simulator P(upper) within 0.5-1.3 Monte
Carlo standard errors and quantiles within 4.5 ms.

### R side

`BOU(boundary_collapse = c("fixed", "exponential", "linear_additive",
"linear_multiplicative", "weibull"))`, mirroring `ROU()`. Non-fixed forms add
`aInf` and `tau` (plus `pw` for Weibull) AFTER the nine required columns, so the
positional column contract is untouched and a fixed-bound BOU is unchanged by
the collapse existing. The form reaches C++ through the `c_name` suffix
(`BOU_BEXP`, `BOU_BWEIB`, `BOU_BLIN_ADD`, `BOU_BLIN_MULT`), read by
`bou_bnd_kind_from_type`; `bou_row` gates the optional columns on it.

Three fixes made along the way, all pre-existing:

* `init_ddm_shared_state` resolved only the REQUIRED columns into `cols`, so a
  DDM-family model with optional trailing columns read past the end of the
  vector -- a segfault, found the first time a collapsing BOU likelihood ran. It
  now walks `keep_names` positionally with padding, the way the race path
  already did.
* the simulator drew `t0` from the CENTRED range while the solver's quadrature
  uses the lower-edge `Uniform(t0, t0+st0)`, putting the two half an `st0` apart;
* `calc_ll_oo_pw` validated the column prefix against `emc2col::ddm::spec()`
  rather than the adapter's, so BOU's `beta` column went unchecked there.

Still hardcoded to `emc2col::ddm::spec()`: the materialized fallback in
`c_log_likelihood_DDM` (`particle_ll.cpp:2385`). BOU never reaches it -- dispatch
happens in `calc_ll_oo`/`calc_ll_oo_pw` -- but a future DDM-family model that
does would silently get a Wiener likelihood there.

`BOU()` also now declares `compress_ok = FALSE`, as `ROU()` and `RLF()` do: the
solve is cached per parameter tuple and `rt` only picks readout points off a
march already paid for, so binning saves nothing -- and it keeps `make_emc`'s
default `rt_resolution = 1/60` from evaluating the density at floored RTs.

## Verification beyond Stage 1

- `beta > 0` against the paper's own method: the Buonocore/Fortet recursion
  (S&R A5a-A7; Smith 2000 eqs 47a/47b) at `delta = 10 ms`, as a test-only
  reference. **Do not calibrate the PDE against Monte Carlo** -- MC is a sanity
  check only.
- `nx` refinement on a pinned domain: clean second-order decay on both fluxes. A
  first-order tail would indicate the boundary-layer singularity and would
  motivate the Richardson pairing used in the RLF solver
  (`model_RLF.h:1546-1631`) -- note that is **not** present in the FPE path
  today, so if added it must extrapolate both fluxes.
- No race regression: `test-fpe-fht.R`, `test-rou.R`, plus a timing check on an
  ROU fit to confirm the `if constexpr` branches compile out.
- End to end: `design()` -> `make_emc()` -> short `fit()` with recovery at known
  `(v, a, Z, beta, t0)`. Use `rt_resolution = NULL` in recovery studies, or the
  density is evaluated at floored RTs and the recovery is biased.

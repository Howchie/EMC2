# Plan: Fokker–Planck (Crank–Nicolson) first-passage solver for leaky accumulators

Date: 2026-07-29
Target branch: `bawl-correlated-race`

## Context

The ported Volterra solver (`src/utils_reducible_diffusion.h`, `src/model_OU_Volterra.h`,
`src/model_BM_Volterra.h`) is now mathematically correct — five kernel/CDF bugs were
fixed on 2026-07-29 and it validates against closed forms and Monte Carlo. But timing
shows it cannot ever be a likelihood:

```
OU chunked, 200 RTs, fixed bound:   max_t=1.0    172 ms
                                    max_t=2.0    662 ms
                                    max_t=4.0   2613 ms    (clean O(max_t^2))
OU chunked, moving bound 1 -> 0.5,  max_t=2.0   1697 ms
```

Per accumulator per design cell per parameter set. A modest race design is ~8 cells
=> ~5 s per likelihood evaluation; EMC2 needs ~1 ms. That is 3–4 orders of magnitude,
and it is structural: the Volterra history sum is O(N_t^2).

The Fokker–Planck route is **O(N_x · N_t)** — linear in `max_rt` instead of quadratic —
and yields pdf and cdf from one solve. It also subsumes the "one solver, four models"
goal more cleanly than the Volterra (Psi, Phi) route, because it never requires
reducibility to Brownian motion.

The two open Volterra defects (OU chunked mid-RT ~1e-3; `lambda*t` tail degeneracy from
`v_max = 1 - exp(-lambda*t) -> 1`) are both artefacts of the Volterra grid layout and
simply do not exist in this formulation. They are deliberately not fixed.

**Outcome**: a new, independent FPE solver reachable from R, validated against the
existing (now-trusted) Volterra solver, the closed forms, and the simulators. The
Volterra code is retained untouched as the accuracy oracle.

---

## Mathematical specification

### Model class

One Itô diffusion, one absorbing upper boundary `a(t)`:

```
dX = A(x,t) dt + B(x,t) dW
```

| model | A(x,t) | B | note |
|---|---|---|---|
| BM | `mu` | `sigma` | |
| OU | `-lambda * (x - theta)` | `sigma` | |
| GBM | — | — | solve `Y = log X` => **BM** with `A = mu - sigma^2/2`, `B = sigma` |
| Gompertz | — | — | solve `Y = log X` => **OU**; reuse `transform_gompertz_to_ou()` (`utils_reducible_diffusion.h:119`) |

So the PDE core needs **only two coefficient functions**; GBM/Gompertz are a state
log-transform applied before the solve and to the boundary. Critically, **B is constant
in x for all four models** after transform, which removes the `d(D)/dxi` term below.

Boundary `a(t)`: reuse the existing collapsing-bound family verbatim —
`exp_decay_scalar()` (`utils_reducible_diffusion.h:56`), `default_boundary_decay()` (:396),
`fixed_boundary_decay()` (:400), `evaluate_boundary_decay()` (:404). Same `(b0, binf, tau, pow)`
parameterisation as the Volterra path, so parameters transfer 1:1.

### Governing equation

Sub-density `p(x,t)` of paths not yet absorbed, on `x in (x_lo, a(t))`:

```
dp/dt = -d/dx [ A p ] + (1/2) d2/dx2 [ B^2 p ]
```

with `p(a(t), t) = 0` (absorbing) and no-flux at `x_lo`.

### Moving-boundary normalisation

Map to a fixed domain with `L(t) = a(t) - x_lo`, `xi = (x - x_lo)/L(t) in [0,1]`, and
`q(xi,t) = L(t) * p(x,t)` (so `integral q dxi = integral p dx`):

```
dq/dt = -d/dxi [ Atil q ] + (1/2) d2/dxi2 [ Btil^2 q ]

Atil(xi,t) = ( A(x_lo + xi*L, t) - xi * L'(t) ) / L(t)
Btil(t)    = B / L(t)                        (no xi dependence)
D(t)       = 0.5 * Btil^2                    (constant in xi)
```

Because `D` is xi-independent, the flux is exactly `F = Atil*q - D*dq/dxi` with **no
effective-drift correction** — state this in a comment, it is the reason the scheme
stays simple.

Code-level:

```cpp
// BM  (also GBM after log-transform)
Atil = (pars.mu - xi * Lp) / L;
// OU  (also Gompertz after log-transform)
Atil = (-pars.lambda * (x_lo + xi * L - pars.theta) - xi * Lp) / L;
// both
D    = 0.5 * (pars.B / L) * (pars.B / L);
```

`L'(t)` by central difference on `evaluate_boundary_decay`, or exactly 0 when `fixed_b`.

### FPT density = boundary flux

At an absorbing boundary `p = 0`, so the moving-frame correction `-a'(t) p` vanishes and

```
g(t) = -D * dq/dxi |_{xi=1}
```

This is the single strongest reason to use a **finite-volume** discretisation: `g(t)`
*is* the numerical face flux `F_M`, computed to the scheme's own order, with no separate
one-sided derivative hack.

### Discretisation

Cell-centred FV on a **graded mesh** of `M` cells with faces `xf_0 = 0 < ... < xf_M = 1`,
widths `dx_i`, centres `xc_i = (xf_i + xf_{i+1})/2`:

```
dq_i/dt = -( F_{i+1} - F_i ) / dx_i
```

The balance is exact for any cell widths, which is why finite volume is the right
framework here; only the stencil coefficients become cell-dependent.

**Grading.** Cells are uniform in a stretched coordinate `eta` and clustered at the
absorbing barrier:

```
xi(eta) = 1 - sinh(c (1 - eta)) / sinh(c),   c = acosh(grade)
```

`grade` is the ratio of far-field to barrier cell width; `grade = 1` recovers the uniform
mesh exactly. **Default `grade = 8`** — a flat optimum across BM, OU, Gompertz, collapsing
bounds and uniform starts (see the log). This matters because essentially all of the
accuracy lives in the last few percent of the domain: the barrier is at `xi = 1` and, once
the domain is sized as below, the start point sits just under it, while the lower 80–90%
holds only the far tail. All mesh geometry is precomputed once per solve in `FPE_Mesh`;
the time march only ever touches `D(t)` and the two affine drift coefficients.

Interior faces `j = 1..M-1` at `xi = xf_j` — **Scharfetter–Gummel** (exponentially
fitted) flux, which is central at low Péclet and upwind at high Péclet, so small
`sigma` does not oscillate:

```cpp
inline double bern(double z) {                 // z/(exp(z)-1), Bern(0)=1
    if (std::abs(z) < 1e-8) return 1.0 - 0.5 * z;
    if (z >  700.0) return 0.0;
    if (z < -700.0) return -z;
    return z / std::expm1(z);
}
// dc_j = xc_j - xc_{j-1};   P = Atil(xf_j) * dc_j / D
F_j = (D / dc_j) * ( bern(-P) * q[j-1] - bern(P) * q[j] );
```

Use the identity `Bern(-z) = Bern(z) + z`, so the two face weights cost one transcendental
between them rather than two.

Boundary faces:

```cpp
F_0 = 0.0;                                     // no-flux far field (exact mass identity)
F_M = D * (cA * q[M-1] - cB * q[M-2]);         // 2nd-order one-sided, q(1) = 0
```

`F_M` fits `u(s) = c1 s + c2 s^2` through `(0, 0)`, `(s1, q_{M-1})`, `(s2, q_{M-2})` with
`s = 1 - xi`, `s1 = 1 - xc_{M-1}`, `s2 = 1 - xc_{M-2}`, and takes `F_M = D u'(0)`:

```
cA = s2 / (s1 (s2 - s1)),   cB = s1 / (s2 (s2 - s1))
```

which for a uniform mesh collapses to the familiar `(9 q_{M-1} - q_{M-2}) / (3h)`.

Do **not** reach for a higher-order face stencil — tested and rejected, see the log.

Because `F_0 = 0`, mass loss is exactly the absorbed probability:

```
CDF(t_n) = 1 - sum_i dx_i q_i^n
```

so the cdf is free. **Consistency check**: compare this against the trapezoidal integral
of `g(t)` over the same time grid — the two agree only if the scheme is converged, so
their difference is a genuine per-solve error estimate (the mass identity alone is a
tautology and is *not* a check).

### Time stepping

Crank–Nicolson on the tridiagonal operator, `Thomas` algorithm, `O(M)` per step:

```
(I - dt/2 * Lop^{n+1}) q^{n+1} = (I + dt/2 * Lop^{n}) q^{n}
```

**Rannacher start is mandatory** — CN rings on the near-singular initial condition. Take
the first 2 full steps as backward Euler at `dt/2`, then switch to CN.

**Fixed-boundary fast path.** When `b(t)` is constant, `L`, `L'`, `D` and the affine drift
coefficients are all constant, so the spatial operator never changes. Better still, a
Rannacher half-step (`c = step = dt/2`) and a CN full step (`c = step/2 = dt/2`) use the
*same* `c`, so the left-hand side matrix is constant for the entire march: build it once
and factorise it once (`FPE_Tri::factor` / `FPE_Tri::solve`). This is worth ~2.75x and is
what pays for the grading.

### Initial condition

- **Uniform start-point variability** (`z0` = upper limit of `Uniform(0, z0)`): exact
  cell-overlap averaging, so the discontinuous edges do not cost an order of accuracy.
  No delta. This replaces the whole `averaged_image_term` / `calculate_f_term` apparatus.

  **The start is uniform in the PHYSICAL state.** For the geometric models, which solve
  `Y = log X`, its density in the solver coordinate is `e^y/(Zhi - Zlo)`, not flat. Take
  the cell mass in the physical variable and both cases are handled exactly.

- **Point start** (`z0 = 0` in the existing convention): do not seed a delta. Seed at a
  short `t_seed` from the analytic method-of-images absorbed Gaussian in a
  frozen-coefficient local frame (`A(z,0)` constant, `B` constant):

  ```
  p(x, t_seed) = phi(x; z + A t, B^2 t)
                 - exp(2 A (a - z)/B^2) * phi(x; 2a - z + A t, B^2 t)
  ```

  The image mean is `2a - z + A t` — reflect the START point about the barrier, *then*
  drift. Integrate exact cell masses (via `erfc`), not midpoint samples; that is what
  licenses a small `t_seed`.

  Size `t_seed` from `s = B sqrt(t_seed) ~= 4 * min(dx at z, 1/M) * L`, capped so the
  barrier is at least `4s` away and at `0.25 * t_max`. The `min(..., 1/M)` cap matters:
  `t_seed` is where the frozen-coefficient approximation substitutes for the solver, and a
  start point landing in the coarse far field must not be made to pay for the grading.

### Domain lower edge

```cpp
// BM/GBM
x_lo = z_min + std::min(0.0, A * t_max) - 6.0 * sigma * std::sqrt(t_max);
// OU/Gompertz
x_lo = std::min(z_min, theta)
       - 6.0 * sigma * std::sqrt((1.0 - std::exp(-2.0*lambda*t_max)) / (2.0*lambda));
```

Sized as (lowest reachable **mean**) minus 6 diffusive sd, *not* as a blanket multiple of
`|A| t_max + sigma sqrt(t_max)`: every unit of slack at the bottom is resolution stolen
from the barrier. Drift may only *lower* the floor. For BM(mu=1, sigma=1, t_max=2) this
gives `L = 9.5`, against `L = 28.3` for the naive rule and an 8x worse cdf error.

Default sizing `M = 256`, `N_t = 512`, `grade = 8`; all exposed. Measured cost for a
fixed bound: **0.6 ms** for 200 RTs at `max_t = 2`.

---

## Implementation (as built)

### `src/fpe_solver.h` — model-agnostic core

No Rcpp, no model knowledge. The model is a **template parameter, never a
`std::function`**: the Volterra path pays ~330 ns per kernel evaluation almost entirely
because `KernelFn` (`utils_reducible_diffusion.h:72`) is an un-inlinable indirect call
inside its innermost loop.

- `bern`, `bern_from_exp`, `bern_neg` — Bernoulli function and its `Bern(-z) = Bern(z)+z`
  companion.
- `FPE_Tri` — tridiagonal system split into `factor()` and `solve()`, so the fixed-boundary
  fast path can eliminate once for the whole march.
- `FPE_Mesh` — graded mesh geometry, built once per solve: faces, centres, widths and
  reciprocals, face spacings, the Péclet basis `pu/pv`, and the absorbing-face weights
  `cA/cB`. `local_dx(xi)` sizes the seed.
- `FPE_Op` — the tridiagonal spatial operator plus `D(t)`; grown but never re-zeroed,
  since `build_op` overwrites every entry and runs once per step.
- `build_op` — SG fluxes plus the two boundary faces. On a uniform mesh `Atil` is affine
  in `xi`, so the face Péclet numbers are an arithmetic sequence and `exp(P_j)` follows by
  one multiply, removing both `expm1` calls per face. A graded mesh breaks the progression
  and pays one `exp` per face.
- `flux_out`, `apply_shifted`, `grid_lookup`.
- `fpe_solve` — Rannacher start + CN loop, operator pointers swapped rather than copied;
  returns `FPE_Result { t, pdf, cdf, flux_mass_mismatch }`.

### `src/fpe_models.h` — the four models

`FPE_Boundary` (`b(t)`, `b_prime(t)`, and their log-state forms `a`, `a_prime`) plus two
model structs, `FPE_ModelBM` and `FPE_ModelOU`, each supplying `x_lo`, `B`, `drift`,
`length`, `length_prime`, `static_op` and `atil_affine`. `drift` survives only for
seeding; the march uses `atil_affine`.

GBM and Gompertz are not separate models — they are the same two structs with a log-state
boundary and transformed parameters, set up at the call site in `fpe_diffusion.cpp`. The
Gompertz map matches `transform_gompertz_to_ou()` (`utils_reducible_diffusion.h:119`),
restated here so this TU need not include the Volterra headers.

Also `fpe_x_lo_bm` / `fpe_x_lo_ou`, `fpe_norm_mass`, `fpe_seed`, and `fpe_run`
(build mesh → seed → march).

### `src/fpe_diffusion.cpp` — the only TU

**Does not include the Volterra headers** — those define non-inline free functions owned
by `volterra_diffusion.cpp`. Cross-validation happens in R, not C++.

```
fpe_bm_fht_pdf_cdf_vec(t, mu, sigma, z0, b0, binf, tau, pow, nx, nt, grade)
fpe_ou_fht_pdf_cdf_vec(t, lambda, theta, sigma, z0, b0, binf, tau, pow, nx, nt, grade)
fpe_gbm_fht_pdf_cdf_vec(t, mu, sigma, z0, b0, binf, tau, pow, nx, nt, start_floor, grade)
fpe_gompertz_fht_pdf_cdf_vec(t, alpha, beta, z0, k0, kinf, tau, pow, nx, nt,
                             start_floor, grade)
```

NB `grade` comes **after** `start_floor` for the geometric models. Each returns
`list(pdf, cdf, mismatch, t_grid, pdf_grid)`.

**Gotcha: `shlib.mk` has no header dependency tracking — `rm -f src/fpe_diffusion.o
src/EMC2.so` before every rebuild.** Verify with `R CMD INSTALL` into a temp lib and print
`getNamespaceInfo("EMC2","path")`; never trust a bare `library(EMC2)`. Rcpp cannot parse a
namespaced C++ default, so `grade` defaults to a literal `8.0` in the export signatures.

### Validation ladder (`.Rtmp/fpe_validate.R`)

Run in order; each must pass before the next matters.

1. **BM fixed bound** vs analytic Wald.
   `wald_pdf <- function(t,mu,b) b/sqrt(2*pi*t^3)*exp(-(b-mu*t)^2/(2*t))`
   `wald_cdf <- function(t,mu,b) pnorm((mu*t-b)/sqrt(t)) + exp(2*mu*b)*pnorm(-(b+mu*t)/sqrt(t))`
2. **OU, `b0 == theta`** vs `ou_fht_pdf_vec_closed_form` / `ou_fht_cdf_vec_closed_form`.
   *Note this case is nearly vacuous for Volterra (it zeroes the kernel) but is a real
   test for FPE* — it does not degenerate here.
3. **OU, `b0 != theta`** (`theta=2, b0=1`) vs `ou_fht_pdf_vec` (Volterra, **single-solve**)
   and vs `simulate_ou_hit_times_bb`. Only the single-solve Volterra is a reference; the
   chunked path is approximate and known to have mid-RT error, so it is not an oracle.
4. **Collapsing bounds** vs `simulate_bm_hit_times_bb` / `simulate_ou_hit_times_bb`.
5. **Start-point variability** `z0 = 0.5` vs `simulate_bm_hit_times_bb`.
6. **GBM / Gompertz** vs `simulate_gbm_hit_times_bb` / `simulate_gompertz_hit_times_bb`.
7. **`flux_mass_mismatch`** on every case above.
8. **High `lambda*t`**: `lambda=5, t=2` — the Volterra failure case. Must converge
   monotonically in `nx`/`nt` (Volterra oscillated 0.13/0.48/1.0).
9. **Convergence order**: halve `h` and `dt`, confirm ~4x error reduction (2nd order).

**MC trap**: `simulate_*_hit_times_bb` returns `NA` for paths that never hit within
`t_max`. Divide the empirical cdf by the full `n`, not by the count of finite hits —
renormalising over the hitters invents a ~4% bias against the solver.

### Bake-off (`.Rtmp/fpe_bakeoff.R`)

Accuracy-at-fixed-wall-clock against the Volterra chunked path (`max_t in {1,2,4}`, fixed
and moving bound, 200 RTs).

**Timing trap**: in R a promise **memoises** on first force, so
`for (i in 1:n) force(expr)` times one call and divides by `n` — 5x-optimistic. Use
`eval(substitute(expr), parent.frame())`.

**Decision gate:**
- FPE reaches `<= 1e-4` rel on the mid-RT band in `<= 5 ms` at `max_t = 2` => adopt.
- FPE lands `1e-4 .. 1e-3` in `<= 5 ms` => adopt, and record the accuracy ceiling.
- FPE cannot beat ~50 ms => stop and reconsider; only then revisit the Volterra
  micro-optimisations (`std::function` -> template, `std::map` -> indexed vector),
  which are otherwise wasted effort on a path we are retiring.

### `tests/testthat/test-fpe-fht.R`

Deterministic assertions only — the MC arms of the ladder stay in `.Rtmp/fpe_validate.R`,
being slow and non-deterministic.

Harness gotcha: running `test_dir` against an *installed* namespace throws bogus
"could not find function" errors for unexported helpers — pass
`env = new.env(parent = asNamespace("EMC2"))` and `TESTTHAT_PARALLEL=false`.

---

## Explicitly out of scope

- Wiring into model dispatch / `design()`. The FPE solver stays reachable from R only,
  exactly as the Volterra code is today.
- Fixing OU chunked mid-RT error, and the `lambda*t` tail degeneracy. Both are Volterra
  grid artefacts that this formulation does not have.
- Volterra performance work. The bake-off gate settled this: Volterra is retired as a
  likelihood and kept as the accuracy oracle, where its speed does not matter.
- A BNR/Smith (Buonocore–Nobile–Ricciardi) regular-kernel Volterra solver. Was held in
  reserve as a third oracle in case FPE and Volterra disagreed with neither clearly right;
  never needed, since FPE self-converges where they differ.

---

# Implementation log

Branch `bawl-correlated-race`. Commits `9ce2af3e` (core + 4 models + exports),
`80c62bac` (physical-vs-log uniform start; `expm1` removed from the inner loop),
`9e56b7b4` (tests + bake-off), plus the graded mesh and fixed-boundary fast path.

Everything above this line is the **as-built** specification, corrected in place; this
section records results and the reasoning that is not visible from the code.

## Status

Complete. Reachable from R, **not** wired into dispatch / `design()` — same scope as the
Volterra path. Full suite: **830 pass, 0 fail, 0 error**, 103 `On CRAN` skips, 2
pre-existing warnings.

## Validation ladder — all items pass

| item | result |
|---|---|
| 1. BM vs Wald | cdf 2.6e-5 … 2.3e-4 abs (nx=256, four parameter sets) |
| 2. OU `b0 == theta` vs closed form | pdf 1.9e-4 … 1.3e-3 at nx=512, all lambda |
| 3. OU `b0 != theta` vs single-solve Volterra + MC | agrees with MC to ≤8e-4 (MC noise); see below |
| 4. Collapsing bounds vs MC | BM 5.4e-4, OU 9.8e-4 (MC se ~7.9e-4) |
| 5. Start-point variability vs MC | BM 1.5e-3, OU 1.4e-3 |
| 6. GBM / Gompertz vs MC | 1.5e-3 / 1.0e-3 |
| 7. `flux_mass_mismatch` | 1e-13 … 1e-14 once resolved |
| 8. high `lambda*t` (Volterra failure case) | converges monotonically 0.9085→0.9096 (Volterra oscillated 0.13/0.48/1.0) |
| 9. convergence order | cdf error ratios 3.69 / 3.85 / 4.10 |

**Item 3.** FPE and Volterra differ ~0.3% at the OU pdf peak (`lambda=2, theta=2, b0=1`).
FPE is self-converged to 1e-4 across nx 256→1024 and Volterra is not; MC confirms the cdf
for both to 8e-4 but cannot resolve a localised 0.3% pdf error. The single-solve Volterra
is the only Volterra reference worth citing — the chunked path is approximate and carries
known mid-RT error. The planned BNR/Smith third oracle was therefore never needed.

## Bake-off — decision gate: ADOPT

Accuracy at matched wall clock, BM vs Wald, mid-RT band `t in [0.15, 1.5]`, `max_t = 2`:

| solver | config | wall clock | mid-RT rel |
|---|---|---|---|
| FPE | nx=128 nt=256 | 0.20 ms | 4.2e-3 |
| FPE | nx=256 nt=512 (default) | 0.60 ms | 1.3e-3 |
| FPE | nx=384 nt=768 | 1.40 ms | 5.6e-4 |
| FPE | nx=512 nt=1024 | 2.40 ms | 3.3e-4 |
| Volterra | base_panels=150 | 1320 ms | 2.9e-7 |

Timing vs `max_t` (OU, 200 RTs, fixed bound): FPE 0.8 / 0.6 / 0.8 ms against Volterra
407 / 1110 / 3456 ms at `max_t` 1 / 2 / 4 — **509x / 1850x / 4319x**. The speedup grows
with `max_t` because Volterra's history sum is O(N_t^2) while the march is O(M·N_t).

**It is a trade, not a clean win.** Per solve Volterra is ~3 orders more accurate and ~3
orders slower. FPE is the only one of the two that can ever be called inside a sampler;
Volterra stays as the accuracy oracle.

**Moving-bound asymmetry.** A collapsing bound costs FPE 5.3x (0.60 → 3.20 ms) against
Volterra's 2.6x — but only because the fixed case got so much faster. The operator is
time-varying, so neither the constant-LHS factorisation nor the arithmetic-Péclet
recurrence applies, and every step pays an `exp` per face. In absolute terms it is still
~900x faster than Volterra's 2920 ms.

## The graded mesh

At converged `nt`, `grade = 8` buys a clean **5.2x** on both cdf and pdf at every `nx`,
and the scheme stays exactly 2nd order (ratios 3.7–4.1). The grade sweep is a flat optimum
at 6–12 across BM, OU (three lambdas), Gompertz, collapsing bounds and uniform starts;
beyond ~16 the far field starves and the error turns back up.

Two things are easy to get wrong here:

- **The gain does not show up at low `nt`.** A first sweep at `nx=256, nt=512` saw only
  1.6x and appeared to plateau past grade 6 — that plateau was the *time* error taking
  over once the space error had been cut. Sweep grading against a converged `nt`, then
  rebalance.
- **Grading initially made the point-start seed worse.** Sizing `t_seed` by the local cell
  width means a start point in the coarse far field gets a *larger* `t_seed` and leans
  harder on the frozen-coefficient formula. Gompertz (start range 500x wide in the
  physical state, so the start sits low in the domain) regressed from 2.3e-3 to 5.8e-3
  against the point-start mixture oracle. Capping at the uniform width fixed it and turned
  the regression into a 4x gain (2.3e-3 → 5.7e-4).

Net effect on the frontier: the old uniform-mesh 1.2e-3 @ 4.2 ms became 1.3e-3 @ 0.6 ms
and 3.3e-4 @ 2.4 ms — roughly an order of magnitude in accuracy per millisecond, of which
the fixed-boundary factorisation is ~2.75x and the grading the rest.

## Dead ends — tested and rejected, do not retry

- **Higher-order absorbing face.** Replacing the 2nd-order one-sided flux with the
  3rd-order 4-point stencil `D*(3.75 q[M-1] - (5/6) q[M-2] + 0.15 q[M-3])/h` made the pdf
  **worse** at every resolution (4.8e-3 → 6.5e-3 at nx=256; 1.19e-3 → 1.63e-3 at nx=512).
  `q` itself is only 2nd order, so differentiating it harder amplifies solution error.
- **Seed width.** Sweeping `FPE_SEED_CELLS` over 1.0 / 1.5 / 2.0 / 3.0 / 4.0 moved the pdf
  error only 4.0e-3 → 4.8e-3. Left at 4.0.

## Remaining headroom

The residual is the scheme's global `O(h^2 + dt^2)`, uniform in `t`, which shows up as a
large *relative* error only in the early tail (`t = 0.05`, density ~1e-3 of peak). With the
mesh now graded and the fixed-bound operator factorised once, the obvious next levers are
the moving-bound path (which pays full price every step) and the point-start seed (whose
frozen-coefficient error is what the mixture oracle in the test suite measures at 5.7e-4).
Neither is currently the binding constraint on using this as a likelihood.

## Verification

```bash
cd /data/work/EMC2_dev_oo
rm -f src/fpe_diffusion.o src/EMC2.so
Rscript -e 'Rcpp::compileAttributes()'
R CMD INSTALL --library=.Rtmp/Rlib_fpe .
Rscript .Rtmp/fpe_validate.R      # ladder, must be clean through item 9
Rscript .Rtmp/fpe_bakeoff.R       # decision gate
TESTTHAT_PARALLEL=false Rscript -e '.libPaths(c(".Rtmp/Rlib_fpe",.libPaths())); \
  library(EMC2); testthat::test_dir("tests/testthat", env=new.env(parent=asNamespace("EMC2")))'
```

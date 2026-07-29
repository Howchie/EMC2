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

Cell-centred FV, `M` cells, `h = 1/M`, centres `xi_i = (i + 1/2) h`, `i = 0..M-1`:

```
dq_i/dt = -( F_{i+1} - F_i ) / h
```

Interior faces `j = 1..M-1` at `xi = j*h` — **Scharfetter–Gummel** (exponentially
fitted) flux, which is central at low Péclet and upwind at high Péclet, so small
`sigma` does not oscillate:

```cpp
inline double bern(double z) {                 // z/(exp(z)-1), Bern(0)=1
    if (std::abs(z) < 1e-8) return 1.0 - 0.5 * z;
    return z / std::expm1(z);
}
// P = Atil_j * h / D
F_j = (D / h) * ( bern(-P) * q[j-1] - bern(P) * q[j] );
```

Boundary faces:

```cpp
F_0 = 0.0;                                        // no-flux far field (exact mass id.)
F_M = D * (9.0 * q[M-1] - q[M-2]) / (3.0 * h);    // 2nd-order one-sided, q(1)=0
```

`F_M` is the 2nd-order one-sided derivative through the three points
`(s=0, q=0), (s=h/2, q_{M-1}), (s=3h/2, q_{M-2})` with `s = 1 - xi`; its coefficients
are `(-8*u0 + 9*u1 - u2)/(3h)`.

Because `F_0 = 0`, mass loss is exactly the absorbed probability:

```
CDF(t_n) = 1 - h * sum_i q_i^n
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
the first 4 steps as backward Euler at `dt/2`, then switch to CN.

### Initial condition

- **Uniform start-point variability** (`z0` = upper limit of `Uniform(0, z0)`): directly
  `q = 1/z0` on `[0, z0]`, zero elsewhere. Smooth; no delta. This replaces the whole
  `averaged_image_term` / `calculate_f_term` apparatus.
- **Point start** (`z0 = 0` in the existing convention): do not seed a delta. Seed at a
  short `t_seed` from the analytic method-of-images absorbed Gaussian in a frozen-coefficient
  local frame (`A(z,0)` constant, `B` constant), using `gaussian_pdf`/`gaussian_cdf` from
  `src/gaussian.h`:

  ```
  p(x, t_seed) = phi(x; m, s^2) - exp(2*A*(a - z)/B^2) * phi(x; 2a - m, s^2)
  m = z + A*t_seed,   s = B*sqrt(t_seed)
  ```

  Pick `t_seed` so `s ~= 4*h*L` (Gaussian resolved by the mesh); shrink it if the
  boundary is closer than `4s`. Initialise `CDF(t_seed)` to the matching Wald cdf.

### Domain lower edge

```cpp
// BM/GBM
x_lo = std::min(0.0, z_min) - 8.0 * (std::abs(mu) * t_max + sigma * std::sqrt(t_max));
// OU/Gompertz
x_lo = std::min(z_min, theta)
       - 8.0 * sigma * std::sqrt((1.0 - std::exp(-2.0*lambda*t_max)) / (2.0*lambda));
```

Default sizing `M = 256`, `N_t = ceil(t_max / dt)` with `dt = t_max/512`; both exposed.
Estimated cost ~1.3e5 tridiagonal ops => **~0.1–0.3 ms**.

---

## Implementation order

### Step 1 — `src/fpe_solver.h` (model-agnostic core)

No Rcpp, no model knowledge. Contains:

- `struct FPE_Coeffs { double D; std::function-free drift via templated functor; }` —
  **pass the drift as a template parameter, not `std::function`.** The existing
  `KernelFn = std::function<...>` in `utils_reducible_diffusion.h:72` is invoked 3x per
  panel inside an O(N^2) loop and is a measured ~330 ns/eval; do not repeat that mistake.
- `bern(double)`, `thomas_solve(a, b, c, d, out)`.
- `struct FPE_Grid { int M; double h; std::vector<double> xi_centre; }`.
- `build_tridiagonal(...)` assembling `Lop` rows from the SG fluxes + the two boundary faces.
- `fpe_step(...)` — one BE or CN step.
- `fpe_solve(...)` — Rannacher start + CN loop; returns
  `struct FPE_Result { std::vector<double> t, pdf, cdf; double flux_mass_mismatch; }`.

### Step 2 — `src/fpe_models.h` (parameter prep + boundary)

- `struct FPE_Params { double A_const, lambda, theta, sigma, z_lo, z_hi, b0, binf, tau, pow; bool fixed_b, sp_var, log_state; std::string model; };`
- `prepare_fpe_bm(...)`, `prepare_fpe_ou(...)`, `prepare_fpe_gbm(...)`, `prepare_fpe_gompertz(...)`
  — mirroring the argument order of the existing `prepare_bm_params` (`model_BM_Volterra.h:23`)
  and `prepare_ou_params` (`model_OU_Volterra.h:33`) so R-level calls are drop-in.
  - GBM: `log_state = true`, `A_const = mu - sigma^2/2`, boundary `log(b0)`, start `log(z)`.
  - Gompertz: `log_state = true`, then `transform_gompertz_to_ou` for `lambda/theta/sigma`.
- `fpe_boundary(t, pars)` / `fpe_boundary_prime(t, pars)` wrapping `evaluate_boundary_decay`
  (log-transformed when `log_state`).
- Two drift functors `BMDrift`, `OUDrift` satisfying the Step 1 template interface.

### Step 3 — `src/fpe_diffusion.cpp` (the only TU; mirrors `volterra_diffusion.cpp`)

**Must not include the volterra headers** — those define non-inline free functions owned
by `volterra_diffusion.cpp`. Validation against them happens in R, not C++.

Exports:

```
fpe_bm_fht_pdf_cdf_vec(t, mu, sigma, z0, b0, binf, tau, pow, nx, nt)
fpe_ou_fht_pdf_cdf_vec(t, lambda, theta, sigma, z0, b0, binf, tau, pow, nx, nt)
fpe_gbm_fht_pdf_cdf_vec(...)
fpe_gompertz_fht_pdf_cdf_vec(...)
```

each returning `list(pdf, cdf, mismatch)`. Interpolate the internal time grid onto the
requested `t` with the existing `rd_lookup_grid_value` pattern
(`utils_reducible_diffusion.h:2146`) — reimplement locally (3 lines) rather than include.

Then `Rcpp::compileAttributes()` and rebuild. **Gotcha: `shlib.mk` has no header
dependency tracking — `rm -f src/fpe_diffusion.o src/volterra_diffusion.o src/EMC2.so`
before every rebuild.** Verify with `R CMD INSTALL` into a temp lib and print
`getNamespaceInfo("EMC2","path")`; never trust a bare `library(EMC2)`.

### Step 4 — validation ladder (`.Rtmp/fpe_validate.R`)

Run in order; each must pass before the next matters.

1. **BM fixed bound** vs analytic Wald — target `< 1e-6` rel on pdf and cdf.
   `wald_pdf <- function(t,mu,b) b/sqrt(2*pi*t^3)*exp(-(b-mu*t)^2/(2*t))`
   `wald_cdf <- function(t,mu,b) pnorm((mu*t-b)/sqrt(t)) + exp(2*mu*b)*pnorm(-(b+mu*t)/sqrt(t))`
2. **OU, `b0 == theta`** vs `ou_fht_pdf_vec_closed_form` / `ou_fht_cdf_vec_closed_form`.
   *Note this case is nearly vacuous for Volterra (it zeroes the kernel) but is a real
   test for FPE* — it does not degenerate here.
3. **OU, `b0 != theta`** (`theta=2, b0=1`) vs `ou_fht_pdf_vec` (Volterra, single-solve)
   **and** vs `simulate_ou_hit_times_bb`. Volterra is the tighter oracle (~1e-6 self-
   converged); MC floors at ~5e-4.
4. **Collapsing bounds** vs `simulate_bm_hit_times_bb` / `simulate_ou_hit_times_bb`.
5. **Start-point variability** `z0 = 0.5` vs `simulate_bm_hit_times_bb`.
6. **GBM / Gompertz** vs `simulate_gbm_hit_times_bb` / `simulate_gompertz_hit_times_bb`.
7. **`flux_mass_mismatch`** < 1e-8 on every case above.
8. **High `lambda*t`**: `lambda=5, t=2` — the Volterra failure case. Must converge
   monotonically in `nx`/`nt` (Volterra oscillated 0.13/0.48/1.0).
9. **Convergence order**: halve `h` and `dt`, confirm ~4x error reduction (2nd order).

### Step 5 — bake-off (`.Rtmp/fpe_bakeoff.R`)

Accuracy-at-fixed-wall-clock against the Volterra chunked path, on the timing cases
already measured (`max_t in {1,2,4}`, fixed and moving bound, 200 RTs).

**Decision gate:**
- FPE reaches `<= 1e-4` rel on the mid-RT band in `<= 5 ms` at `max_t = 2` => adopt.
- FPE lands `1e-4 .. 1e-3` in `<= 5 ms` => adopt, and record the accuracy ceiling.
- FPE cannot beat ~50 ms => stop and reconsider; only then revisit the Volterra
  micro-optimisations (`std::function` -> template, `std::map` -> indexed vector),
  which are otherwise wasted effort on a path we are retiring.

### Step 6 — `tests/testthat/test-fpe-fht.R`

Encode ladder items 1, 2, 3, 7, 9 with loose-but-real tolerances (no MC in the suite —
too slow and non-deterministic). Follow the existing style in
`tests/testthat/test-lba-logspace.R`.

Harness gotcha: running `test_dir` against an *installed* namespace throws bogus
"could not find function" errors for unexported helpers — pass
`env = new.env(parent = asNamespace("EMC2"))` and `TESTTHAT_PARALLEL=false`.

---

## Explicitly out of scope

- Wiring into model dispatch / `design()`. The FPE solver stays reachable from R only,
  exactly as the Volterra code is today.
- Fixing OU chunked mid-RT error, and the `lambda*t` tail degeneracy. Both are Volterra
  grid artefacts that this formulation does not have.
- Volterra performance work — conditional on the Step 5 gate (see above).
- A BNR/Smith (Buonocore–Nobile–Ricciardi) regular-kernel Volterra solver. Held in
  reserve as a third independent oracle, to be built **only if** FPE and Volterra
  disagree in Step 4.3/4.4 and neither is clearly right.

## Verification summary

```bash
cd /data/work/EMC2_dev_oo
rm -f src/fpe_diffusion.o src/volterra_diffusion.o src/EMC2.so
Rscript -e 'Rcpp::compileAttributes()'
R CMD INSTALL --library=.Rtmp/Rlib_fpe .
Rscript .Rtmp/fpe_validate.R      # ladder, must be clean through item 9
Rscript .Rtmp/fpe_bakeoff.R       # decision gate
Rscript -e '.libPaths(c(".Rtmp/Rlib_fpe",.libPaths())); testthat::test_dir("tests/testthat", env=new.env(parent=asNamespace("EMC2")))'
```

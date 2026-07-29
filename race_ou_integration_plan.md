# Integrating the new FHT solvers as EMC2 race models

**Scope.** Systematic evaluation of the two new solver families (reducible
diffusions via Fokker–Planck; Race Lévy Flight) and a concrete pathway to wire
them into the package's race-model machinery. The proof of concept is the
**racing Ornstein–Uhlenbeck (leaky accumulator)**; the architecture generalises
to BM/GBM/Gompertz and RLF with no further structural work.

Written 2026-07-29 against `36d2e521`.

---

## 0. Executive summary

* The FPE solver (`src/fpe_solver.h`, `src/fpe_models.h`) is **ready to be a
  likelihood backend**. Measured cost is 0.2–0.6 ms per solve at usable
  accuracy, and one solve serves *every* RT that shares a parameter vector.
* The RLF solver (`src/model_RLF.h`) is **not ready** — 25–210 ms per solve,
  with an internal refinement loop that `stop()`s on non-convergence. It needs
  a fixed-cost mode before it can sit under MCMC. Defer it.
* The **only real architectural mismatch** is that EMC2's race kernels are
  pointwise (`RaceRawFun`, one row = one accumulator × one trial), while the
  PDE solvers are grid-valued (one solve = a whole `t` grid). This is resolved
  by a **group-solve adapter**: partition rows by their parameter tuple, solve
  once per group, interpolate. That is a new `*_raw` kernel and nothing else —
  no change to `c_log_likelihood_race`, the design system, or the samplers.
* Four concrete defects must be fixed before wiring up (§2.3). The most
  important: the `flux_mass_mismatch` diagnostic is **not** a usable accuracy
  proxy, and the solver returns small **negative** densities in the far tail.
* Parameterisation: estimate **(v, k)** — accumulation rate and leak — not
  (lambda, theta). This makes the racing OU **nest the RDM exactly at k = 0**,
  which gives a free correctness oracle, a free fast path, and a sane default.

---

## 1. What exists, and what "integration" actually requires

### 1.1 The two new code paths

| | reducible diffusions | Race Lévy Flight |
|---|---|---|
| numerical core | `src/fpe_solver.h` (474 ln) | `src/model_RLF.h` (764 ln) |
| models | `src/fpe_models.h` — BM, OU, GBM, Gompertz | RLF only |
| R entry points | `src/fpe_diffusion.cpp` (4 exports) | `src/rlf_diffusion.cpp` (2 exports) |
| tests | `tests/testthat/test-fpe-fht.R` (27 assertions, pass) | `tests/testthat/test-rlf-fht.R` (59 assertions, pass) |
| simulators | `simulate_{ou,bm,gbm,gompertz}_hit_times_bb` (Volterra TU) | `simulate_rlf_hit_times_cpp` |
| superseded | `utils_reducible_diffusion.h` + `model_{OU,BM}_Volterra.h` | — |

Both are deliberately disconnected: `src/fpe_diffusion.cpp:10-12` states the
entry points are "reachable from R for development and validation only". The
Volterra path stays as an independent cross-validation oracle — it should
**not** be deleted, and should **not** become the likelihood path (509–4319×
slower per the FPE plan log).

### 1.2 What a race model must supply

The contract is four function pointers plus a column spec, assembled in
`resolve_race_model_adapter()` (`src/particle_ll.cpp:485-624`) and consumed by
`c_log_likelihood_race()` (`src/particle_ll.cpp:7130`):

```c
// src/utils.h:29-40
typedef double (*RacePdf1Fun)(double rt, const double* par, void* ctx);   // GSL integrand
typedef double (*RaceCdf1Fun)(double rt, const double* par, void* ctx);   // GSL integrand
typedef void (*RaceRawFun)(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok,
                           double* out, double min_ll, void* ctx);        // batch, log scale
typedef void (*RaceLogSAtTFun)(double t, const double* const* cols, int n_rows_total,
                               int n_lR, int n_par, const int* trunc_mask,
                               int n_unique_trials, const int* isok_all,
                               void* ctx, double* logS_out);              // truncation norms
```

Answering the question in the brief directly: **yes**, the "scalar" functions
are just single-RT versions of the same math (compare `drdm_scalar`,
`src/utils.h:202`, against `drdm_raw`, `src/utils.h:298`). They exist only
because GSL's `gsl_function` needs a `double -> double` callback for the
censoring integrals. They are *not* free for a PDE model — see §3.4.

`model_pfun_raw` returns the **log-survivor**, not the log-CDF. `logS_at_t` is
the same quantity summed over a trial's accumulators at one scalar `t`.

### 1.3 What a model definition file must supply

From `RDM()` (`R/model_RDM.R:213-241`) the required list elements are:
`type` (`"RACE"`), `c_name`, `p_types` (named defaults, **order must equal the
`col_registry.h` enum**), `p_types_canonical`, `transform$func`,
`bound$minmax` + `bound$exception`, `Ttransform`, `rfun`, `dfun`, `pfun`,
`log_likelihood`. Plus a `col_registry.h` namespace, an `@export` roxygen block
with the parameter table, and a `NAMESPACE` entry.

---

## 2. Evaluation

### 2.1 FPE solver — cost

`fpe_ou_fht_pdf_cdf_vec`, OU with `lambda=2, theta=2, sigma=1, b0=1`, 200 output
times, elapsed per solve (this machine):

| M \ nt | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|
| 64  | 0.05 | 0.10 | 0.15 | 0.30 |
| 128 | 0.10 | 0.20 | 0.30 | 0.60 |
| 256 | 0.20 | 0.35 | 0.60 | 1.15 |
| 512 | 0.35 | 0.65 | 1.15 | 2.20 |

(ms). Cost is linear in `M*nt` as designed, and independent of how many RTs are
requested. **M=256, nt=512 → 0.6 ms** is the working point (§2.2).

Cost projection for a real fit: a `forstmann`-style design
(`v~lM, s~lM, B~E+lR, A~1, t0~1`) has **12 unique (accumulator × design-cell)
parameter rows** covering all 524 compressed dadm rows for one subject. So
12 × 0.6 ms ≈ **7 ms per particle per subject**. That is roughly three orders of
magnitude above the RDM's analytic kernels, but it is a *fixed* cost set by the
design, not by the data size — and it is unambiguously tractable.

The hard constraint that follows: **cost scales with the number of distinct
parameter rows, not trials.** A trend/covariate model (`R/trend.R`) makes every
trial's parameters unique and turns 12 solves into 4433. Trends must be
rejected at `design()` time for these models until a batched-solve strategy
exists.

### 2.2 FPE solver — accuracy

Reference = self-converged solve at M=1024, nt=4096. Leaky-accumulator regime
`k=4` equivalent (`lambda=4, theta=2, sigma=1, z0=0, b0=1`), 60 times over
[0.05, 2.5]:

| M \ nt | 256 | 512 | 1024 |
|---|---|---|---|
| 128 | 6.7e-3 | 2.1e-3 | 1.3e-3 |
| 256 | 5.1e-3 | 1.4e-3 | 4.1e-4 |
| 512 | 5.5e-3 | 1.4e-3 | 3.2e-4 |

(max abs CDF error). Second-order in `dt`, saturating in `M` at about 256.
Cross-checked against 4×10^5 Brownian-bridge Euler paths: max |empirical −
reference| = 8.5e-4 at `dt=1e-4` and 5.1e-4 at `dt=2e-5`, i.e. the MC path is
the less accurate of the two.

**The log-survivor holds relative accuracy far into the tail** — this is the
single most important result for race integration, because a race likelihood
multiplies loser survivors and needs `log S` to stay accurate where `S` is
tiny. Measured `|Δ log S|` at M=256/nt=512 vs the reference:

| S | 1.0 | 0.18 | 3.7e-3 | 1.0e-6 | 2.7e-10 | 7.5e-14 |
|---|---|---|---|---|---|---|
| Δ log S | 0.001 | 0.003 | 0.001 | 0.004 | 0.009 | 0.014 |

This works because the mass route computes `S = sum(dx*q)` directly as a sum of
positive quantities, never as `1 - CDF`. Below `S ≈ 1e-14` the sub-density
underflows and `S` becomes exactly 0 → `log S = -Inf`; must be floored.

Density relative error in the useful range is 1e-3 to 3e-3 at M=256/nt=512.

### 2.3 FPE solver — the four blocking defects (all now FIXED)

**(a) `flux_mass_mismatch` is not an accuracy proxy.** `src/fpe_solver.h` claimed
the two CDF routes "agree only when the scheme has resolved the solution".
Measured: at M=256/nt=512 the mismatch is **4.3e-8** while the true CDF error is
**1.4e-3**. Both routes are functionals of the same `q`, so a discretisation
error common to both cancels out of their difference. It detects an unstable
march or a mis-built operator — worth having — but says nothing about
resolution. *Fixed:* the header now states this explicitly and forbids its use
as a runtime accept/reject gate.

**(b) Negative densities in the far tail.** At `t=2.251`, where the true pdf is
1.5e-17, the solver returned **−2.99e-16**; `log()` of that is `NaN`, which
would silently poison a likelihood. `1 - sum(dx*q)` was likewise not reliably
inside [0,1] nor monotone out there. *Fixed:* the flux is floored at 0 at the
source, and the reported CDF is clamped to [0,1] and forced non-decreasing. The
**raw** unclamped value still feeds the mismatch diagnostic, so clamping cannot
mask an unstable march.

**(c) Small-`t` error — misdiagnosed in the first draft, and the correction
matters.** It is *not* the frozen-coefficient seed. It is the **uniform time
grid**, and it is a function of `dt` alone: at M=256, relative pdf error at
t=0.05 is 83% / 26% / 8% / 2.6% / 1.3% for dt = 9.8 / 4.9 / 2.4 / 1.2 / 0.61 ms,
*identically* at `t_max`=0.25 and `t_max`=2.5. A uniform `dt = t_max/nt` spends
its steps in proportion to elapsed time, but the density climbs from ~0 to its
mode inside the first ~100 ms and then decays smoothly for a second or more.
With `t_max` set by the largest RT in a data set, the rising flank — where `t0`
is identified — is starved. The seed residual, once `dt` converges, is only ~6%
at t=0.03 (where the density is 1e-3 against a mode of 7.2) and 0.8% at t=0.05.

*Fixed:* **graded time stepping** (`FPE_TimeSchedule`, `tgrade`), the time-axis
analogue of the existing graded spatial mesh. `dt` is held constant within
blocks and doubles between them, so `tgrade` is the ratio of the last step to
the first. Piecewise-constant `dt` is what preserves the one-time O(M)
factorisation of the fixed-boundary operator: it is redone once per block,
i.e. `log2(tgrade)` times for the whole march, rather than once per step. CN is
second order for any step sequence. Measured at M=256, nt=512, `t_max`=2.5 —
**same step budget, same 0.65 ms**:

| tgrade | rel. err @ t=0.05 | max CDF err | max ⎮Δlog pdf⎮ (pdf>1e-4) | max ⎮Δlog S⎮ (S>1e-3) |
|---|---|---|---|---|
| 1 (uniform) | 26% | 2.2e-3 | 2.35 | 3.4e-3 |
| 8 | 2.8% | 2.7e-4 | 0.45 | — |
| **32 (default)** | **1.2%** | **1.8e-4** | **0.107** | **1.5e-3** |
| 64 | 1.0% | 2.1e-4 | 0.082 | — |

`tgrade=32` improves the CDF 12×, the winner log-density 22×, *and* the
practically-relevant log-survivor. The cost is confined to the extreme tail:
⎮Δlog S⎮ for `S < 1e-6` degrades from 1.4e-2 to 8.1e-2. That is the right trade
— an 0.08 error on a `log S` of −30 is 0.3%, whereas 2.35 nats on a winner
density is fatal. `FPE_TGRADE = 32`; `tgrade <= 1` reproduces the uniform grid
exactly, so every pre-existing test still measures what it says it does.

**(d) Sub-`t_seed` lookup returned a constant.** `grid_lookup` clamped to
`vg.front()` below the grid, handing back the flux *at* `t_seed` for every
earlier time — a positive constant where the density is ~0. *Fixed:*
`zero_below_grid` for the PDF; the CDF keeps the clamp, since `cdf[1]` is mass
genuinely absorbed before `t_seed`.

All four are committed with regression tests (`test-fpe-fht.R`, now 39
assertions; the 27 pre-existing ones pass unchanged, and `test-rlf-fht.R`'s 59
are untouched).

### 2.4 RLF solver — evaluation and the speed problem

Correct and well-instrumented (59 assertions covering the α→2 Brownian limit,
mass/flux consistency, monotone CDF, positivity, CMS-simulator agreement, and
explicit rejection of out-of-budget regimes). Three things block it:

* **25–210 ms per solve** (`v=2, sigma=1, alpha=1.7, b0=1`), 40–350× the FPE.
* **Cost is parameter-dependent and discontinuous.** `nx=200` triggers domain
  expansion to `nx_used=660` and costs 140 ms; `nx=400` satisfies the check
  immediately (`refinement_skipped=TRUE`) and costs 25 ms. A likelihood whose
  cost jumps 5× across a parameter boundary is hostile to adaptive MCMC.
* **It `stop()`s when refinement fails**, aborting the sampler rather than
  returning `min_ll` for a bad proposal.

Directions, in expected order of payoff:

1. **Exploit the Toeplitz structure.** The Grünwald–Letnikov weights `g_k`
   depend only on `|i-j|`, so the nonlocal operator is Toeplitz. That makes a
   matvec O(nx log nx) by FFT instead of O(nx²), and lets the CN solve run as
   preconditioned GMRES with a circulant preconditioner — retiring the
   O(nx³) dense LU (`RLF_DenseLU`) entirely. At `nx=660` that is ~2.9e8 flops of
   factorisation and ~3.5e8 of back-substitution both replaced by ~1e7.
2. **Replace the refinement loop with an a priori `(nx, nt)` rule**, calibrated
   offline over the `(v, sigma, alpha, b0)` grid against the current adaptive
   solver. Removes both the 5× discontinuity and the `stop()`.
3. **Tighten the lower closure.** The conservative no-flux censoring at the
   lower edge is what forces the wide domain. Matching to the known `x^{-alpha}`
   power-law tail instead would let `nx` shrink substantially.

Only after (1) and (2) is RLF a candidate; the §3 architecture then applies to
it unchanged. Not scheduled here.

---

## 3. Architecture

### 3.1 The mismatch, and the fix

Analytic race models evaluate `f(rt_i | par_i)` independently per row. A PDE
solver amortises one march over a whole time grid; evaluating it per row would
be ~500× more expensive than necessary.

**Group-solve, backed by a solve cache.** All four entry points in the race
contract — `dfun_raw`, `pfun_raw`, `logS_at_t`, and the scalar `pdf1`/`cdf1` —
are served by one mechanism:

```
solve_cache : parameter tuple  ->  { t_grid, log_pdf_grid, log_S_grid, t_max }
```

* **Lookup** with a required `t_max`: if the cached solve reaches at least that
  far, interpolate; otherwise re-solve to the new horizon and replace.
* **Batch kernels** (`dfun_raw`, `pfun_raw`) make one O(n_rows × n_par) pass to
  hash rows into parameter groups, take `t_max = max(rt - t0)` per group, then
  fill every row by interpolation.
* **Scalar functions** (`pdf1`, `cdf1`) hash their single parameter row and hit
  the same cache — so a `cdf1` call is an interpolation, not a solve.
* **Invalidation** by an epoch counter on `ContextForRaceModels`, bumped once
  per particle at the ~3 particle-loop sites in `calc_ll_oo`.

This is entirely contained in the model's kernels: `c_log_likelihood_race`
never learns that the model is grid-based. Note the cache is not an
optimisation to be deferred — §3.4 shows it is what makes censoring viable at
all, so it is phase-1 core.

Grouping is data-fixed in practice (parameters are constant within a design
cell), so it could be precomputed into `RaceSharedState`
(`src/particle_ll.cpp:686`, built at `:2845`) — but hashing per call is ~50 µs
against milliseconds of solving, so do the stateless version and only precompute
if profiling asks for it.

### 3.2 Choosing the grid per group

`t_max` varies per group, so a fixed `nt` gives a parameter-dependent `dt`.
Scale it: `M = 256`, `grade = 8`, `tgrade = 32`, and
`nt = clamp(ceil(t_max / 2e-3), 256, 4096)`. Graded time means the *effective*
step on the rising flank is ~`dt/32`, which is what §2.3(c) measured as
sufficient. Expose `M`, `dt_target`, `grade` and `tgrade` as package options
(`emc2.fpe_nx` / `emc2.fpe_dt` / `emc2.fpe_grade` / `emc2.fpe_tgrade`) so
accuracy can be swept against a real fit without a recompile.

### 3.3 Trends and covariates: slow, not forbidden

A trend (`R/trend.R`) makes every trial's parameters unique, so the 12 solves of
a `forstmann`-style design become one per trial. That is slow — hundreds of
times slower — but there is nothing *wrong* with it, and the group-solve
mechanism degenerates to one-solve-per-row gracefully and correctly. So:
**no error.** Emit a one-time `message()` at `make_emc()` reporting the number
of distinct parameter rows and the implied per-particle cost, so the user is
choosing the cost knowingly. This is the same courtesy as the existing
"Likelihood speedup factor: 3.6, 4433 unique trials" report.

### 3.4 Censoring and truncation are in scope for phase 1

Reading the actual code paths (`src/particle_ll.cpp:7483` and the
"other trials" loop at `:8003`) confirms these are much cheaper than the first
draft assumed:

* **Unknown winner** — the overwhelmingly common case — is
  `log_surv_cm(t, ...)`, which is `sum_k log(1 - cdf1(t, par_k))`: one scalar
  CDF evaluation per accumulator at `LT` / `LC` / `UC` / `UT`. Upper censoring
  is `log_diff_exp(log S(UC), log S(UT))`, lower censoring
  `log_diff_exp(log S(LT), log S(LC))`, and the truncation normaliser is the
  same quantity at the truncation bounds. **No integration at all.**
* **Known winner** (`R_j_idx != NA`) and the active-nogo branch enter
  `integrate_interval` → GSL over the defective win density. This is the Go/NoGo
  case and is rare in practice.
* GSL is also the *fallback* when `log_diff_exp` underflows.

So censoring costs a handful of extra `cdf1` calls at a handful of distinct `t`
values, on parameter rows that are already in the cache. It is nearly free. The
GSL path is the only expensive one, and even there the memo turns each integrand
evaluation into an interpolation, provided the cached solve reaches the
integration upper limit — so `integrate_interval` must warm the cache to `upp`
before handing GSL the callback.

`logS_at_t` is one scalar `t` across all trials, served by the same cache.

### 3.5 R-side dfun / pfun / rfun

* `dRACEOU(rt, pars)` / `pRACEOU(rt, pars)` — R mirrors used by `make_data()`,
  `predict()` and the plotting path. Same grouping in R: `split()` on the
  parameter rows, one `fpe_ou_fht_pdf_cdf_vec()` call per group. `pfun` returns
  the CDF (R convention) even though the C++ `pfun_raw` returns the log-survivor.
* `rfun` — `simulate_ou_hit_times_bb()` already exists and is Brownian-bridge
  corrected. Wrap it as `rRACEOU(lR, pars, ok)` following `rRDM`
  (`R/model_RDM.R:101`): scalar-parameter simulator, so loop over parameter
  groups, assemble the `n_acc × n_trials` finishing-time matrix, take
  `max.col(-t(dt))`, add `t0`. Register `.rfun_RACEOU` in `R/model_rng.R`
  alongside `.rfun_RDM` (`R/model_rng.R:42`). A C++ kernel in `model_rng.cpp`
  is optional and not needed for phase 1.

---

## 4. Parameterisation of the racing OU

### 4.1 The model

Solver form: `dX = -lambda (X - theta) dt + sigma dW`, absorbed at `b(t)`,
started `X_0 ~ U(0, A)`.

Leaky-accumulator form (Usher & McClelland 2001; Smith & Ratcliff 2004,
appendix — `Math/Smith & Ratcliff (2004).pdf`): `dX = (v - k X) dt + sigma dW`,
`v` the stimulus input and `k` the leak. Same process, with

```
lambda = k,    theta = v / k
```

### 4.2 Estimate (v, k), not (lambda, theta)

1. **The k → 0 limit is the Wiener race.** With `theta = v/k` the limit is a
   coordinate singularity (`theta → ∞`); with `(v, k)` the drift `A(x) = v - kx`
   is perfectly regular and at `k = 0` it is constant, i.e. Brownian motion with
   drift. So `k` has a meaningful zero and the racing OU **nests the RDM**.
2. `v` is what the design system wants contrasts on. `v ~ lM` means "input
   differs by stimulus match"; `theta ~ lM` would confound input with leak.
3. `v` and `k` are close to orthogonal in the likelihood; `theta` is a ratio of
   the two.

**Required solver change** (small): give `FPE_ModelOU` a `(v, k)` form.
`atil_affine` becomes

```c
a0 = (v - lambda * xlo) / L;      // was -lambda*(xlo - theta)/L
a1 = (-lambda * L - Lp) / L;      // unchanged
```

and `fpe_x_lo_ou` takes the lowest reachable mean as `min(z_min, v/k)` when
`k > eps` and `z_min + min(0, v*t_max)` otherwise — matching `fpe_x_lo_bm`.
`theta` becomes a derived quantity used only for domain sizing. The existing
`theta`-based R entry points stay as they are, for the Volterra cross-checks.

### 4.3 No code fallback at k = 0

The nesting is **mathematical, not a dispatch branch**. At `k = 0` the affine
drift coefficients above collapse to exactly `FPE_ModelBM`'s, so the OU model
solves the Brownian problem through its own PDE code — no branch, no special
case. It must **not** call the analytic Wald kernel, for two reasons:

* it would couple two unrelated code paths, with the usual consequences; and
* more importantly it would **destroy the oracle**. The value of `k = 0` is that
  it lets an independent implementation be checked against a trusted one. A
  fallback makes that check vacuous.

This follows the precedent already set by RDMSWTN, which likewise does not nest
back into RDM in code and is validated against it instead. (It is the *opposite*
of the BAwL/LBA arrangement, where LBA genuinely is the `k=0` member of the same
kernel — that choice was made to avoid duplicating one formula, and it is not
the right choice when the whole point is cross-validation.)

### 4.4 Two clarifications on the geometry

**Start-point variability is always positive.** `z_lo = 0`, `z_hi = max(0, A)`
(`src/fpe_diffusion.cpp:67`), so `X_0 ~ U(0, A)` exactly as in the RDM and LBA.
There is no negative start point and `A` is bounded below by 1e-4 as elsewhere.

**There is no reflecting lower barrier, and this is a real difference from the
LCA.** Usher & McClelland's leaky *competing* accumulator rectifies activation at
zero — a genuine reflecting boundary that is part of the model. Here the process
is free to go negative, and only positive boundary hits are counted. The
solver's domain does have a no-flux face, but it sits at `x_lo`, roughly
`6 sigma sqrt(var)` below the lowest reachable mean, purely as a truncation of
an infinite domain — it is placed where essentially no mass reaches (~1e-9), not
at a psychologically meaningful zero.

The consequence is that this model is the **leaky accumulator** of Smith &
Ratcliff's appendix (an OU race), not the LCA: no rectification and no
lateral inhibition between accumulators. Worth stating plainly in the roxygen,
because "leaky accumulator" is used loosely in the literature for both. Adding
rectification later is a *reflecting* condition at `x = 0` — a one-line change
to the `i = 0` face in `build_op` — but it changes the model, so it should be a
separate named variant, not a parameter.

### 4.5 Collapsing boundaries: a pluggable boundary module

Collapsing bounds are a primary use case, not an extension, and they are the
thing these solvers can do that the analytic race models cannot. They are in
scope. The cost is real — a moving boundary makes `static_op()` false, forfeiting
the one-time factorisation — but it buys something unavailable elsewhere.

The architecture should follow the kill/guess-clock precedent: a
**`boundary_collapse` argument to the model constructor, defaulting to
`"fixed"`**, selecting a boundary *form*, with each form contributing its own
named parameters to `p_types` as optional columns after `N_REQ` (exactly how
`bawl` handles `mG`/`mK`/`omega`, `src/col_registry.h:72`).

The key design point is the one raised in review: **the boundary functions must
be independent of both the solver and the simulator.** They answer two questions
and nothing else:

```c
struct FPE_Boundary {
  int kind;                        // FIXED | WEIBULL | EXPONENTIAL | LINEAR | ...
  double b0, binf, p1, p2;         // form-specific shape parameters
  double b(double t) const;        // boundary at time t
  double b_prime(double t) const;  // its derivative
  bool fixed() const;              // enables the static-operator fast path
};
```

`build_op` calls `length()`/`length_prime()` once per *time step* (not per cell),
so a `switch` on `kind` inside `b()` is free. Neither the solver nor the
simulator learns any parameter names; the **model definition** owns the mapping
from named parameters to `(kind, b0, binf, p1, p2)`, in `Ttransform`. This is
what lets the *same* boundary module serve `simulate_ou_hit_times_bb`, which
today hard-codes the Weibull form — that shared use is the main argument for
factoring it out rather than growing the existing struct in place.

Forms to define now: `"fixed"`. Forms to leave stubbed but architecturally
reachable: `"weibull"` (`binf, tau, pw` — the current hard-coded form),
`"exponential"` (`binf, tau`), `"linear"` (`binf, tau`). Only `"fixed"` needs to
work in phase 1; what must exist in phase 1 is the *seam*.

### 4.6 Proposed `p_types`

| Parameter | Transform | Natural scale | Default | Mapping | Interpretation |
|---|---|---|---|---|---|
| `v`  | log | [0, ∞) | `log(1)` | | Accumulation rate (stimulus input) |
| `k`  | log | [0, ∞) | `log(0)` | | Leak rate (1/s); `k = 0` is the Wiener race |
| `B`  | log | [0, ∞) | `log(1)` | `b = B + A` | Threshold above the start-point range |
| `A`  | log | [0, ∞) | `log(0)` | `z0 = A` | Start-point range, `X_0 ~ U(0, A)` |
| `t0` | log | [0, ∞) | `log(0)` | | Non-decision time |
| `s`  | log | [0, ∞) | `log(1)` | | Within-trial diffusion sd (fix to 1) |
| `pContaminant` | probit | [0, 1] | `qnorm(0)` | | Generic nuisance infrastructure |

Optional columns, gated on `boundary_collapse` (§4.5): `binf`, `tau`, `pw`.

* Deliberately identical to `RDM()` apart from `k`, so the two are directly
  comparable and an RDM design converts by adding `k~1`.
* `s = 1` fixed for scale identification. `k` has units of 1/time so it is *not*
  absorbed by that scaling — with `s` fixed, `(v, k, B, A)` are identified.
* `bound$minmax`: `v = c(1e-3, Inf)`, `k = c(0, Inf)`, `B = c(0, Inf)`,
  `A = c(1e-4, Inf)`, `t0 = c(0.05, Inf)`, `s = c(0, Inf)`;
  `bound$exception = c(A = 0, v = 0, k = 0)`.
* **`v > 0` enforced.** Negative input is meaningful for a mismatching
  accumulator, but it puts `theta = v/k < 0` far below the barrier, forcing
  `x_lo` down and costing most of the mesh resolution. Revisit only if the
  graded mesh is retuned for that regime.

### 4.7 Column registry

```c
// R/model_RACEOU.R — racing Ornstein-Uhlenbeck (leaky accumulator).
// binf/tau/pw exist only for boundary_collapse != "fixed" and must be gated on
// the context flag, per the convention for optional columns.
namespace raceou {
  enum : int { v = 0, k, B, A, t0, s, N_REQ, binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"v", "k", "B", "A", "t0", "s"};
    return {n, N_REQ, "RACEOU"};
  }
}
```

`c_name = "RACEOU"`. Its branch goes **first** in
`resolve_race_model_adapter()` — dispatch is by substring
(`src/particle_ll.cpp:517`), and while `"RACEOU"` collides with no existing key
today, ordering it first with a comment is what keeps that true.

---

## 5. Implementation steps

~~Step 0~~ **DONE:** the four §2.3 defects are fixed and regression-tested.

**Steps 1–9 are DONE** (model name `ROU`, not `RACEOU`). What was built, and how
each differs from the plan above, is recorded in §7. The remainder:

1. ~~**`(v, k)` reparameterisation** of `FPE_ModelOU` and `fpe_x_lo_ou`.~~ DONE.
2. ~~**Extract the boundary module** (§4.5).~~ DONE, with one carve-out: the
   Weibull form is now behind `FPE_BoundaryKind` and `set_kind()`, but
   `simulate_ou_hit_times_bb` was **not** repointed at it. That simulator lives
   in the Volterra TU, which is the independent cross-validation oracle for the
   solver; making the two share a boundary module would couple the oracle to the
   thing it checks. Repoint it only alongside step 12, when a collapse form is
   actually needed in both.
3. ~~**The solve cache** (§3.1).~~ DONE as `fperace::SolveCache`. No epoch
   counter was needed: the adapter outlives the particle loop, so the cache is
   simply cleared at the top of each particle iteration. Keys are compared
   bit-exactly, so the clear bounds memory rather than protecting correctness.
4. ~~**The kernels.**~~ DONE, in `src/model_ROU.h` (`fpe_race.h` ended up holding
   the model-agnostic cache and the solver adapter, which is the cleaner split).
5. ~~**Registry + dispatch.**~~ DONE. One correction to the plan:
   `defective_upper_tail` is **true**, not false — a leaky accumulator whose
   asymptote `v/k` sits below threshold has a genuine never-finish probability.
6. ~~**`R/model_ROU.R`.**~~ DONE.
7. ~~**Validation against the reference simulator.**~~ DONE, but against a *new*
   simulator (`rou_hit_times_vec`) rather than `simulate_ou_hit_times_bb`, for
   the reason in §7.
8. ~~**`k = 0` equivalence at the R level.**~~ DONE — see §7 for the measured
   numbers. The plan's "~1e-6" was wrong: the target is discretisation error,
   not machine precision, so the test asserts *convergence under refinement*
   instead of a fixed tolerance.
9. ~~**Censoring / truncation tests.**~~ DONE for unknown-winner censoring
   against `RDM` at `k = 0`. The known-winner GNG path through GSL is **not**
   yet covered — see §7.

10. ~~**Grid-resolution sweep.**~~ DONE — §8.1. The shipped default moved to
    `nx = 384`.

11. ~~**Parameter recovery.**~~ DONE — §8.2.

12. ~~**Then the collapsing forms.**~~ DONE — §8.3, with one design change on the
    user's instruction: the collapse asymptote is measured from **zero**, not
    from the top of the start-point range.

**Deferred, explicitly:** RLF (needs §2.4 items 1–2 first); BM/GBM/Gompertz
races (same architecture, different `fpe::` model struct and column spec —
mechanical once OU lands); rectified/LCA variant (§4.4); negative `v`.

---

## 6. Risks

| Risk | Assessment |
|---|---|
| ~7 ms/particle/subject is too slow for real fits | Real but bounded, and measured. Solve-cache reuse across `dfun`/`pfun` halves it; parallelism over subjects already exists. Confirm on a full fit at step 11 before committing to more models. |
| Trend models are hundreds of times slower again | Accepted, not forbidden (§3.3). Report the cost at `make_emc()` so it is a knowing choice. |
| Grid error roughens the LL surface | Step 10 measures it directly. Pointwise accuracy is not sufficient evidence — see the `rel_tol` precedent. |
| Graded time degrades the extreme survivor tail | Measured: ⎮Δlog S⎮ goes 1.4e-2 → 8.1e-2 for `S < 1e-6`, while the winner density improves 22×. Right trade, but re-check if a design produces many accumulators with `S < 1e-6` at observed RTs. |
| `k` and `B` trade off | Expected — a leaky accumulator's threshold and asymptote are related. Step 11 will show it; if severe the answer is a prior on `k`, not a reparameterisation. |
| Losing the static-operator fast path for collapsing bounds | **Measured 5.3×, not the 2–3× guessed here, and the reasoning was wrong.** The graded time schedule amortises refactorisation across *`dt` changes*, which happen `log2(tgrade)` times; a moving boundary changes the operator every step, which the block structure cannot help with. `Binf = b0` degenerates to `fixed` and recovers 1.03×. |
| `S` underflows to exactly 0 past ~1e-14 | Floored at `min_ll`. Only bites for RTs far beyond any accumulator's plausible range. |

---

## 7. What was actually built (2026-07-29, steps 1–9)

### 7.1 Files

| File | Role |
|---|---|
| `src/fpe_race.h` | **new.** Model-agnostic solve cache (`Key`, `Entry`, `SolveCache`), the `(v,k)` → solver adapter `rou_solve`, log-space interpolation, and the reference OU simulator `rou_hit_time`. |
| `src/model_ROU.h` | **new.** The five race-contract kernels: `drou_raw`, `prou_raw`, `rou_logS_at_t`, `drou_scalar`, `prou_scalar`, plus `rou_configure_grid`. |
| `src/rou_diffusion.cpp` | **new.** `rou_pdf_cdf_vec`, `rou_hit_times_vec` — the R-facing exports behind `dROU`/`pROU`/`rROU`. |
| `R/model_ROU.R` | **new.** `ROU()`, `dROU`, `pROU`, `rROU`, `.rou_grid()`. |
| `tests/testthat/test-rou.R` | **new.** 27 assertions. |
| `src/fpe_solver.h` | `FPE_Result::surv` added (see §7.3). |
| `src/fpe_models.h` | `FPE_ModelOU` reparameterised to `(v, lambda)`; `FPE_Boundary` given a `kind` seam; `fpe_x_lo_ou` takes `(z_min, v, lambda, sigma, t_max)`. |
| `src/fpe_diffusion.cpp` | Two OU call sites moved to `set_lambda_theta()`; `surv` added to the returned list. |
| `src/col_registry.h`, `src/particle_ll.cpp`, `R/model_rng.R`, `NAMESPACE` | Registration and dispatch. |

### 7.2 Measured behaviour

**k = 0 against the analytic Wald**, through the PDE, no branch (`v=1.5, B=1, A=0`):

| | max abs |
|---|---|
| pdf | 6.1e-4 |
| cdf | 6.8e-5 |
| log pdf | 3.7e-2 (at t = 0.05, the flank) |
| log S | 3.2e-4 |

With `A = 0.5`, against the package's own `dRDM`/`pRDM`: 4.8e-4 / 7.9e-5.

**`s` is scaled out bit-exactly**: `(v,B,A,s) = (1.5,1,0.5,1)` and `(3,2,1,2)`
return `identical()` densities at `k = 2`, confirming `k` is correctly left
un-rescaled (it has units of 1/time).

**Simulator vs solver**, 2×10^5 draws, `v=1.5, B=1, A=0.5`: max |ΔCDF| = 1.9e-3
at `k = 0` and 1.2e-3 at `k = 2`, against an MC standard error of ~1e-3. At
`k = 2` the simulator reports 1.8% of paths never finishing, which the solver
reproduces as defective upper mass.

**Full likelihood path**, forstmann subject 1, `v~lM, B~E+lR, A~1, t0~1`,
524 compressed rows, 810 trials, ROU with `k` fixed at 0 against RDM:

| `emc2.fpe_nx` | `emc2.fpe_dt` | ROU ll | ll − RDM | s/call |
|---|---|---|---|---|
| 128 | 4.0e-3 | −929.283 | +9.010 | 0.004 |
| **256** | **2.0e-3** | **−935.704** | **+2.590** | **0.007** |
| 384 | 1.0e-3 | −937.151 | +1.143 | 0.020 |
| 512 | 5.0e-4 | −937.710 | +0.584 | 0.052 |
| 768 | 2.5e-4 | −938.041 | +0.253 | 0.154 |

(RDM analytic: −938.2939.) The residual falls monotonically toward zero under
refinement, which is what makes it discretisation error rather than a different
model. **At the shipped default the offset is 2.59 nats over 810 trials, i.e.
0.0032/trial.** Cost is 6.9 ms/particle/subject, matching the §2.1 projection of
~7 ms almost exactly.

### 7.3 Design decisions that departed from the plan

**The survivor is returned by the solver, not derived from the CDF.**
`FPE_Result` gained a `surv` vector carrying `sum(dx*q)` clamped and forced
non-increasing. §2.2 established that this quantity keeps relative accuracy to
`S ≈ 1e-14` because it is a sum of positive cell masses; computing `1 - cdf`
instead would throw away exactly the accuracy a race likelihood depends on. The
race kernels return `log(surv)`, never `log1p(-cdf)`.

**Grids are stored and interpolated in log space.** Both the density and the
survivor are close to exponential over a grid interval, so linear interpolation
of the log is far better than of the value — and it is the log that the race
likelihood consumes. A floored neighbour (`LOG_FLOOR = -700`) is never averaged
into a live value.

**A new simulator rather than `simulate_ou_hit_times_bb`.** That function is
parameterised by `(lambda, theta)` and rejects `lambda <= 0`, so it cannot reach
`k = 0` at all — the one case with an independent analytic answer. `rou_hit_time`
steps the exact OU transition with the `k → 0` limits taken through `expm1`, plus
the Brownian-bridge crossing correction, so `k = 0` is an ordinary value there
too. The Volterra simulator remains untouched as the oracle for the solver.

**`emc2.fpe_*` options are honoured by the sampled likelihood, not only by
`dROU`.** This was missing in the first cut and is not cosmetic: without it there
is no way to demonstrate convergence in `(nx, dt)` on a real design, and the
§7.2 table could not have been produced. `rou_configure_grid()` reads the options
once per likelihood call, and `nt_max` follows `dt_target` so that asking for a
finer step is not silently ignored on the long-horizon groups.

**`defective_upper_tail = true`.** §4.6 had this as false. It is wrong: with
`v/k` below threshold a leaky accumulator has real never-finish mass (1.8%
measured at `k = 2` above).

### 7.4 Not yet done

* **Known-winner / GNG censoring through GSL.** The unknown-winner path is
  tested against RDM at `k = 0`; the `integrate_interval` path is reachable but
  uncovered. The concern is specifically that `integrate_interval` should warm
  the cache to its upper limit before handing GSL the callback, so that the
  sweep pays for one solve rather than one per integrand evaluation
  (`rou_scalar_horizon` currently makes this likely but not guaranteed).
* **The §3.3 trend cost message** at `make_emc()`. Trends work — the group-solve
  degenerates to one solve per row correctly — but the user is not yet told what
  that costs.
* **Steps 10–12**: the grid-resolution *surface* sweep (§5.10; the pointwise
  numbers in §7.2 are not sufficient evidence — see the `rel_tol` precedent),
  parameter recovery (§5.11), and the collapsing boundary forms (§5.12).
* **A resolution recommendation.** The default `(256, 2e-3)` costs 0.0032
  nats/trial of bias against a known-exact reference. That is almost certainly
  fine for parameter estimation and almost certainly *not* fine for model
  comparison by information criteria, where a systematic per-trial offset does
  not cancel. Step 10 should settle this.

---

## 8. Steps 10–12 (2026-07-29)

### 8.1 Step 10 — the resolution sweep, and the shipped default

The §5.10 worry, taken from the `gsl_ctl.rel_tol` precedent, was that pointwise
accuracy could be adequate while the log-likelihood *surface* stayed rough
enough to break a sampler. **It does not happen here.** Five parameters were
scanned over 61 points each (a step of ~0.005 in the sampling scale, comparable
to an adaptive DE-MCMC proposal), at twelve `(nx, dt, tgrade)` settings, on
forstmann subject 1 with `v~lM, k~1, B~E+lR, A~1, t0~1`.

Roughness is measured as the RMS residual of a 7-point local quadratic fit,
applied to `LL_setting − LL_reference` rather than to `LL_setting`. That
differencing matters: the raw roughness of a scan is dominated by the *true*
curvature of the likelihood, which is not an error and does not fall under
refinement — measured on the raw scans, `t0` looked *worse* at the reference
resolution than at the coarsest one. Differencing cancels it.

| setting | t0 | k | B | v | v_lMd |
|---|---|---|---|---|---|
| nx=128 dt=4e-3 | 1.5e-1 | 1.1e-1 | 9.3e-2 | 1.7e-8 | 4.7e-9 |
| nx=256 dt=2e-3 | 8.8e-2 | 8.7e-2 | 8.4e-2 | 4.4e-9 | 1.2e-9 |
| **nx=384 dt=2e-3** | **4.7e-2** | **3.1e-2** | **4.0e-2** | **2.0e-9** | **5.9e-10** |
| nx=512 dt=1e-3 | 1.8e-2 | 2.5e-2 | 1.9e-2 | 7.6e-10 | 2.2e-10 |

Against a median |ΔLL| **per scan step** of 17 (t0), 14 (k), 39 (B), 12 (v) and
6.7 (v_lMd) nats. So the jitter is two to three orders of magnitude below the
signal one proposal step carries, at *every* setting tested. The drift
parameters are numerically exact to ~1e-9. There is nothing here for a sampler
to trip over.

**`nx` is the binding constraint, not `dt`.** This corrects the framing of
§2.3(c), which diagnosed the accuracy problem as a `dt` problem — true of a
*uniform* time grid, but `tgrade = 32` already fixed the time axis. Bias against
a converged reference, in nats/trial:

| | dt=4e-3 | dt=2e-3 | dt=1e-3 | dt=5e-4 |
|---|---|---|---|---|
| nx=128 | 1.0e-2 | | | |
| nx=256 | | 2.6e-3 | 2.4e-3 | |
| nx=384 | | 1.1e-3 | 9.5e-4 | |
| nx=512 | | | 4.4e-4 | 3.8e-4 |

Halving `dt` buys 7–16%; raising `nx` by half buys a factor of two. `tgrade` is
saturated at 32: at nx=256 the bias is 2.9e-3 / 2.6e-3 / 2.6e-3 for tgrade =
8 / 32 / 64, so 64 buys nothing and 8 is measurably worse.

**Fit quality is insensitive to all of it.** The MLE was found at each setting
and then *scored on the reference likelihood*. Every one lands within **0.22
nats over 810 trials (2.7e-4 nats/trial)** of the best, including nx=128. Note
the reference setting's own MLE is not the best-scoring one — Nelder-Mead on
nine parameters does not converge tightly enough for that — which is also why
the per-setting MLEs themselves are *not* a usable measure of resolution error.
An earlier reading of this sweep took a 20–30% spread in the fitted `k` as
resolution bias; it is optimiser noise on a correlated ridge. Profiling `k` at
the reference MLE shows it is sharply identified conditionally (−31 nats at
0.5×, −133 nats at 2×), so the spread is a joint-direction effect, not flatness.

**Shipped default: `nx = 384`, `dt = 2e-3`, `grade = 8`, `tgrade = 32`.**
Cost 10.4 ms/particle/subject, 1.5× the old default and 46× the RDM's 0.225 ms
on the same design. The justification is not estimation accuracy — nx=128 is
already adequate for that — but the absolute likelihood value, which is what
information criteria consume and where a systematic per-trial offset does not
cancel. At 384 the offset is 1.1e-3 nats/trial (0.9 nats over 810 trials); at
the old default it was 2.6e-3 (2.1 nats), enough to move an AIC comparison by
the equivalent of a parameter.

Two caveats worth stating in the docs rather than burying:

* Comparing **two ROU variants** to each other is far safer than comparing ROU
  to an analytic model. The bias is nearly constant over a parameter scan — over
  the `k` scan its sd is 1.2e-4 against a mean of 2.6e-3, i.e. 20× smaller — so
  it substantially cancels in a ROU-to-ROU difference.
* For an IC comparison against an analytic model, use `nx = 512` or above.

| resolution | ms/particle | ×RDM | bias (nats/trial) |
|---|---|---|---|
| 128 / 4e-3 | 2.97 | 13 | 1.0e-2 |
| 256 / 2e-3 | 7.03 | 31 | 2.6e-3 |
| **384 / 2e-3** | **10.43** | **46** | **1.1e-3** |
| 512 / 1e-3 | 26.33 | 117 | 4.4e-4 |
| 768 / 2.5e-4 | 153.07 | 680 | reference |

### 8.3 Step 12 — collapsing boundaries

All three forms are implemented in the §4.5 module and reachable from
`ROU(boundary_collapse = )`, which appends `Binf`/`tau` (and `pw` for Weibull)
as optional columns before `pContaminant`, following BAwL's `mG`/`mK`/`omega`
precedent, and suffixes `c_name` with `_BEXP` / `_BLIN` / `_BWEIB`.

**The asymptote is measured from zero, not from `A`.** The first cut anchored it
at `Binf + A` so that the bound could never enter the start-point range; on the
user's instruction it is now `Binf` outright. The bound is therefore free to
descend into `[0, A]` and, in the limit, to meet the start point — a forced
response, which is a state the model should be able to express. `Binf` is
bounded away from zero (1e-3) only so the solver's domain cannot collapse.

`b(0) = B + A` in every form, so `Binf = B + A` is the degenerate case;
`set_kind` detects it and reports `fixed = true`, which restores the one-time
factorisation. That is tested by `expect_identical` against the fixed-bound
model, for all three forms — not an aesthetic check, since without it a sampler
wandering to a non-collapsing boundary would silently pay the full moving-bound
cost.

Cost of losing the static-operator fast path, measured at nx=256/dt=2e-3:

| | ms/solve | × fixed |
|---|---|---|
| fixed | 1.17 | 1.00 |
| exponential | 6.20 | 5.31 |
| linear | 6.17 | 5.29 |
| weibull | 6.40 | 5.49 |
| exponential, `Binf = b0` | 1.20 | 1.03 |

**5.3×, not the 2–3× guessed in §6.** The guess assumed the graded time
schedule's block structure would amortise the refactorisation; it does not. That
structure amortises refactorisation across *`dt` changes*, of which there are
`log2(tgrade)` in the whole march, whereas a moving boundary changes the
operator on *every step*.

The reference simulator was extended rather than replaced: `rou_hit_time_bnd`
evaluates `b` at both ends of each step and applies the Brownian-bridge crossing
probability against the *linear* barrier over the step, which keeps the
correction exact to the order the step already is. Comparing against it at
10^5 draws, max |ΔCDF| is within the 7e-3 tolerance (≈2× the MC standard error)
for all three forms. `simulate_ou_hit_times_bb` in the Volterra TU is still
**not** repointed at the shared boundary module, for the §5.2 reason: it is the
independent oracle for the solver, and sharing code would couple the oracle to
what it checks.

`test-rou.R` is now 54 assertions, all passing.

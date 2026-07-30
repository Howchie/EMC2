# BAwD: thresholded transient-drive ballistic accumulator — integration plan

Source spec: `Math/ballistic_decay.txt`. Target: a new RACE model alongside
LBA / BAwL, sharing the ballistic-accumulator numerical machinery.

The model is `BAwD(drift_distribution = c("lognormal", "normal"))`. Both launch
distributions are built in one pass: the process geometry is identical for both
and only the two integrals over the launch distribution differ, so §3 factors
the kernel into a shared geometry layer plus a launch policy and §4 carries the
selection through `c_name` into the race context.

Every formula below has been checked against Monte Carlo and, where relevant, a
high-precision reference. The scripts are committed and are the starting point
for the test file:

| script | checks |
|---|---|
| `Math/check_bd.R` | truncated-normal launch: pdf, cdf, `P(T=∞)` vs 400k-draw MC |
| `Math/check_lnorm.R` | lognormal launch: same, plus `integrate(f) == F` |
| `Math/check_stoploss.R` | `log C(v)` variants vs a 400-bit `Rmpfr` reference |

**Verification discipline:** `library(EMC2)` loads a stale install. Validate
every step by `R CMD INSTALL --library=<temp lib>` and check the harness path
banner before trusting a number.

---

## 0. Relation to existing models

BAwD is related to BAwL but distinct; the two coincide in exactly one corner of
the parameter space.

| | BAwL | BAwD |
|---|---|---|
| trajectory | `X(t) = z e^{-kt} + (V/k)(1-e^{-kt})` | `X(t) = z + (V/k)(1-e^{-kt}) - ℓ t` |
| what decays | **accumulated evidence** (leak on the state, start point included) | **the drive** `U(t)=V e^{-kt}`; the state never leaks |
| opposing force | proportional to `X` | constant `ℓ` (clearance / sampling floor) |
| asymptote | `V/k` (monotone, approached from below) | rises to a peak at `log(V/ℓ)/k`, then falls |
| never-finish condition | `V ≤ k·b` | `V < V_c`, `V_c = ℓ·r_c`, `r_c - 1 - log r_c = k(b-z)/ℓ` |
| support of RTs | `(0, ∞)`, defective | **`(0, T_max]`, defective — hard right endpoint** |
| density near the top | `f → 0` as `t → ∞` | `f ∝ (T_max - t)` at `A = 0`, `∝ (T_max - t)²` for `A > 0` |

Limits, all confirmed numerically:

* `ell = 0, A = 0` **is** BAwL — reproduces `EMC2:::pleakyba()` to the last
  digit (`k` ↔ leak rate; both give the `V > k·b` condition).
* `ell = 0, A > 0` is a *static-start-point* BAwL, materially different because
  BAwL decays the start point and BAwD does not (`F(0.2)` = 0.0359 vs 0.0100 at
  `v=2, sv=1, b=1, A=0.5, k=1.5`).
* `k → 0` is exactly **LBA with drift `v - ell`**.
* `ell > 0` is the regime with no existing counterpart: a finite response
  window, defective mass with a hard endpoint, and an effectively quadratic
  collapsing bound `B(t) = k t (b + ℓ t)/(1 - e^{-kt})`.

The behavioural distinction from BAwL is the RT support: BAwL turns weak drive
into an arbitrarily slow response, BAwD converts it into an omission. That is a
testable difference in the right tail plus the omission rate.

---

## 1. Parameterisation

### 1a. Truncated-normal launch (`drift_distribution = "normal"`)

BAwL's parameter vector plus one parameter:

| par | transform | default | role |
|---|---|---|---|
| `v` | identity | `1` | mean launch strength (drive at `t=0`) |
| `sv` | log | `log(1)` | between-trial SD of launch strength; fix to 1 for scale |
| `B` | log | `log(1)` | `b = B + A` |
| `A` | log | `log(0)` | start-point range, `z ~ U(0,A)` (**static**, not decayed) |
| `t0` | log | `log(0)` | non-decision time |
| `k` | log | `log(0)` | drive-decay rate (`k = 0` → LBA limit) |
| `ell` | log | `log(0)` | tonic clearance / minimum effective sampling rate |
| `pContaminant` | probit | `qnorm(0)` | as elsewhere |

`V ~ N(v, sv²)`, truncated positive when `posdrift = TRUE`. The decay rate is
named `k` rather than `gamma` so that an existing BAwL design converts to BAwD
by adding `ell ~ 1`; `k` plays the same dimensional role (1/time curvature) in
both models. The mechanistic difference (drive decay vs state leak) belongs in
the roxygen block, using the §0 table.

### 1b. Lognormal launch (`drift_distribution = "lognormal"`, the default)

`log V ~ N(mu, sigma²)`. Same `B`, `A`, `t0`, `k`, `ell`; `v`/`sv` replaced by
`mu` (identity) and `sigma` (log). `V > 0` holds by construction, so the whole
`posdrift` / `IO` / normalizer-floor branch disappears — no truncation constant,
no `LBA_DENOM_FLOOR`, no `pnorm(v/sv)` denominator anywhere. This is the default
because its likelihood is analytic on the ordinary path (§2b).

Scale identification differs from 1a and must be documented prominently. The
evidence axis is defined only up to a scale: `V, ell, b, A → cV, cell, cb, cA`
leaves every crossing time unchanged. In 1a, `sv = 1` fixes it. Here `sigma` is
**dimensionless** and cannot, so one of `ell`, `B`, `A` must be fixed instead.
**Fix `ell = 1`**: it is the natural unit (the processing floor) and it makes the
spec's dimensionless quantities literal, `r = V` and `c = k·b`. Omitting this
produces a ridge, not an error.

### 1c. Rejected alternative, and an identifiability warning

The spec's `(t0, γ, c, π, σ_r)` form — `log r ~ N(µ_r, σ_r²)` with `µ_r`
reparameterised through the response probability `π` — is better identified for a
single accumulator but is not used here. `π` is a *per-accumulator marginal*
quantity, whereas in a race the omission rate is a joint property of all
accumulators, so a parameter named `π` would not mean what it says; and `v`/`b`
are what the `design()` formulas and `lM`/`lR` contrasts are built to carry. The
docs should still note that `r = V/ell` and `c = k·b/ell` are the dimensionless
quantities controlling shape.

Identifiability warning for the roxygen block: at `k = 0` only `v - ell` is
identified, so `ell` is identified purely through curvature and the finite
window. Recommend `ell ~ 1` (or a coarse condition factor) and *not* crossing
`ell` with the same factors as `v`.

---

## 2. Shared geometry and the truncated-normal likelihood

Let `q(u) = (1 - e^{-ku})/k` (`expm1`-safe; `q = u` at `k = 0`). With start
point `z`, crossing at `u` requires launch strength

    V*(u,z) = (b - z + ℓu) / q(u)

which is **U-shaped** in `u`, unlike BAwL's monotone `V*`. Its minimum is the
critical drift; the decreasing branch is the first-passage branch. Define

    c(u) = (v - (b + ℓu)/q(u)) / sv,     m(u) = 1 / (sv q(u)),     E = e^{-ku}
    denom = Φ(v/sv) if posdrift else 1

so that `P(V ≥ V*(u,z)) = Φ(c + m z)/denom` — exactly BAwL's kernel form with
new `c` and `m`, which is what allows the existing numerics to be reused whole.

**Saturation.** With `y = ku`, the accumulator started at `z` peaks at
`T_max(z)` solving `e^y - 1 - y = c_z ≡ k(b-z)/ℓ`, at critical drift
`V_c(z) = ℓ e^y`; for `u > T_max(z)` the hit probability is frozen. Inverting
that relation gives the closed-form boundary start point

    z*(u) = b - (ℓ/k)(e^{ku} - 1 - ku),        Z(u) = min(max(z*(u), 0), A)

so `z < Z` is still live and `z > Z` is frozen. Only two Lambert-W / Newton
solves are needed per row, never per quadrature node:

    y_A solves e^y - 1 - y = k(b-A)/ℓ    (earliest saturation, T_max(A) = y_A/k)
    y_0 solves e^y - 1 - y = k·b/ℓ       (full saturation, T_max = T_max(0) = y_0/k)

Use Newton on the convex `h(y) = e^y - 1 - y - c` (init `sqrt(2c)` for small
`c`, `log(c + 1 + log(c+1))` for large); 4–6 iterations, with no branch to
select. This is `-W_{-1}(-e^{-(1+c)})` computed robustly.

**Density** (`hi = c + m Z`, `ΔΦ = Φ(hi) - Φ(c)`, `Δφ = φ(hi) - φ(c)`):

    f(u) = [ (v·E - ℓ)·ΔΦ  +  sv·E·Δφ ] / (A · denom)          for Z > 0
    f(u) = 0                                                    for Z ≤ 0 (u ≥ T_max)
    f(u) = φ(c) · ((b + ℓu)E - ℓ q) / (sv q² · denom)           A = 0 limit

At `k = 0` this is `[(v-ℓ)ΔΦ + sv Δφ]/(A denom)`, the LBA density with drift
`v - ℓ`, i.e. the existing kernel with one substitution.

**CDF.** With `H(w) = w Φ(w) + φ(w)`:

    F(u) = [ (H(c + mZ) - H(c))/m  +  Ψ(y_A, min(ku, y_0)) ] / (A · denom)
    Ψ(y1,y2) = (ℓ/k) ∫_{y1}^{y2} Φ((v - ℓ e^y)/sv) (e^y - 1) dy

The first term is BAwL's `∫Φ` verbatim with the upper limit `A` replaced by
`Z(u)`. `Ψ` is the frozen contribution of already-saturated start points,
written in `y` rather than `z` so that no root-finding appears inside the
integral: the integrand is smooth, monotone and bounded — 20–32 fixed
Gauss–Legendre nodes via `src/gl_quad.h`. Three regimes:

* `u ≤ T_max(A)`: `Z = A`, `Ψ = 0` — pure closed form, identical in shape to BAwL.
* `T_max(A) < u < T_max`: both terms; `f` is still pure closed form.
* `u ≥ T_max`: `F = F_max = Ψ(y_A, y_0)/(A·denom)`, constant; `f = 0`.

**Defective mass** `P(T = ∞) = 1 - F_max`, available in the same closed/quad
form and exact at `t = Inf` — no `qagiu`, unlike BAwL's clock paths. At `A = 0`
there is no quadrature at all and `F_max = Φ((v - ℓ e^{y_0})/sv)/denom`.

`Math/check_bd.R` compares all of the above against 400k-draw Monte Carlo
(root-found first passages) over four parameter sets including `A = 0`, `k → 0`
and a heavily defective set: CDF agrees to 3–4 decimals throughout, density
within MC noise, `P(T=∞)` to 4 decimals (0.71199 vs 0.71203; 0.93149 vs
0.93190).

### The flat-CDF contract

This is load-bearing for censoring and truncation and must hold exactly:

    for t >= T_max:   F(t) = F_max  (exactly constant),  S(t) = 1 - F_max,  f(t) = 0
    at t = Inf:       F(Inf) = F_max, NOT 1

Nothing about the endpoint is invalid or needs bypassing. An upper censoring
bound or timeout `UC > T_max` is well posed: the response mass in `(UC, ∞)` is
`F_max - F(UC) = 0` and the omission mass is `S(UC) = 1 - F_max`, which is
exactly the intrinsic never-finish mass. Same for `UT`. The existing
defective-tail plumbing (`ctx.defective_upper_tail`,
`log_surv_cm(R_PosInf, ...)` at `particle_ll.cpp:8292–8380`) already expects
this shape from BAwL and needs no change — it only has to receive `F_max` rather
than a clamped `1` or an extrapolated `NaN` from an evaluator that assumed
monotone-to-one. Tests §6.5–6.6 pin the invariant at both regime seams.

One consequence, which `pContaminant` does not address: a **finite** observed RT
above `t0 + T_max` has density exactly zero, so that trial contributes `min_ll`.
This is a constraint on the posterior's support — the supported window must cover
the observed RTs — handled by initialisation and priors, not by a mixture.
`pContaminant` is a Bernoulli omission mixture (`1 - pContaminant` multiplying
finite observations, `pContaminant` the `rt = Inf` mass) and must not be used to
patch late responses. Late contaminant *responses*, if they turn out to be
needed, are a separate uniform-RT mixture component and out of scope here.

### Special cases to branch on explicitly

* `ell ≤ eps`: `c_z → ∞` and the `y` parameterisation blows up. No saturation
  ever occurs: `Z = A`, `Ψ = 0`, `F(∞) = Φ((v - k b)/sv)/denom` (BAwL's
  defective mass with a static start point). Guard `ell` before dividing.
* `k ≤ eps`: `q = u`; `T_max = sqrt(2b/(ℓk)) → ∞`. Route to the LBA-with-`v-ℓ`
  formulas rather than evaluating `e^{ku}-1-ku`, which cancels catastrophically.
* `b = A` and `z = A`: `b - z = 0`, instantaneous crossing — the same edge case
  LBA already has.

---

## 2b. Lognormal launch: analytic on the ordinary path

Under lognormal `V` the frozen integral also has an elementary antiderivative,
so every regime is closed form in `Φ`/`φ`. Let `log V ~ N(mu, sigma²)`,
`M = e^{mu + sigma²/2}` (= `E[V]`),

    Gbar(w) = P(V >= w) = Φ((mu - log w)/sigma)
    d1(w) = (mu + sigma² - log w)/sigma,     d2(w) = (mu - log w)/sigma

**Live part.** Substituting `w = V*(u,z)` (so `dz = -q dw`), with
`w_hi = (b + ℓu)/q(u)` and `w_lo = (b + ℓu - Z(u))/q(u)`:

    I_live = q(u) * [ J(w_hi) - J(w_lo) ],
    J(w)   = w Φ(d2(w)) - M Φ(d1(w)),        J'(w) = Gbar(w)

`-J(w)` is `E[(V-w)_+]`, an undiscounted Black call price; §2c layer 1 covers
its stable evaluation.

**Frozen part.** A single closed antiderivative `K(s)` exists (by parts on the
`e^s` term after completing the square), but **do not implement it that way.**
Factor it instead over the critical-launch interval
`[v_a, v_b] = [ℓ e^{s_A}, ℓ e^{min(k u, s_0)}]`:

    I_frozen = (1/k) * [ I_0 - ℓ·I_{-1} ]
    I_0    = ∫ Gbar(v) dv       over [v_a, v_b]
    I_{-1} = ∫ Gbar(v)/v dv     over [v_a, v_b]

This is the same quantity — substituting `v = ℓ e^s`, `dv = v ds` gives
`(ℓ/k)∫(e^s - 1) Gbar(ℓ e^s) ds`, identical to `K` — but every piece now lands
on a primitive that already exists in `wald_functions.h`:

    I_0    -> C(v) = E[(V-v)_+], the stop-loss / Black call    [NEW, §2c layer 1]
    I_{-1} -> sigma * ∫ Q(x) dx = sigma * [D(x_a) - D(x_b)]
              == log_normal_q_interval(x_a, x_b)               [ALREADY EXISTS]

`D(x) = φ(x) - x Q(x)` is verbatim `log_normal_q_antiderivative_abs()`
(`wald_functions.h:145`), including the `1 - x R(x)` regrouping for `x > 0` and
the `φ(x)/x²` asymptotic branch, so `log I_{-1}` is
`log(sigma) + log_normal_q_interval(x_a, x_b)` and needs no new code at all.
`s_A`, `s_0` are the same two Newton solves; the geometry is unchanged.

    F(u) = (I_live + I_frozen) / A          (no denom: no truncation constant)
    F(Inf) = F_max = I_frozen(s_A -> s_0) / A

**Density** — simpler than the truncated-normal case, because `∫ w g(w) dw` is
the lognormal partial expectation:

    f(u) = [ E·M·(Φ(d1(w_lo)) - Φ(d1(w_hi))) - ℓ·(Φ(d2(w_lo)) - Φ(d2(w_hi))) ] / A
    f(u) = g(w_hi) (w_hi·E - ℓ) / q(u)      for A = 0
    f(u) = 0                                for u >= T_max

`Math/check_lnorm.R` covers three parameter sets including `A = 0` and both
partial and full saturation: CDF vs 400k MC to 3–4 decimals, `P(T=∞)` 0.68670 vs
0.68751 / 0.66143 vs 0.66163 / 0.95965 vs 0.95998, and `integrate(f) == F` to
**9 significant digits** — a sharper internal consistency check than MC, and it
pins the live/frozen seam.

**Validation constraint.** This variant has no exact oracle: `ell = 0, A = 0` is
not BAwL and `k → 0` is not LBA once `V` is lognormal. It inherits the geometry
certified by §6.1–6.3, but its own two integrals are validated only against MC,
`integrate(f) == F`, and direct quadrature. Hence the ordering in §7.

---

## 2c. Numerics: three cancellation layers

Three separate cancellations must be handled: inside `J` itself; in the endpoint
differences `J(a)-J(b)` / `K(a)-K(b)`; and in the frozen combination
`I_0 - ℓ I_{-1}`. All three have stable alternatives, and the pattern already
exists in this tree — `wald_k0_log_cdf_closed()` and `wald_pt_log_surv()` in
`wald_functions.h` do exactly this for the Wald: Mills-ratio regrouping, a
retained-fraction guard at the shared constant `EMC2_LOG_CANCEL_MIN`
(`= log(1e-6)`), and a `NA_REAL` return that hands off to a short quadrature.
Follow that contract rather than introducing a parallel one.

### Already present — do not rewrite

| need | existing |
|---|---|
| `R(x) = Q(x)/φ(x)`, continued fraction past `x = 7.07` | `mills_ratio_std()` |
| `log φ(x)` without underflow | `log_phi_std()` |
| `D(x) = φ(x) - x Q(x)`, incl. `1-xR(x)` and `φ/x²` branches | `log_normal_q_antiderivative_abs()` |
| `log ∫Q dx` = `log[D(lo) - D(hi)]` | `log_normal_q_interval()` |
| `log ∫Φ dz` | `log_normal_phi_integral()` |
| `log(Φ(hi) - Φ(lo))` via the correct tail | `log_normal_interval()` |
| `log1mexp`, `logdiffexp`, signed log arithmetic | `log1m_exp()`, `log_diff_exp()`, `signed_log_*` |
| lognormal CDF/PDF/log-survivor | `plnorm_std()`, `dlnorm_std()`, `lnorm_log_surv_std()` |
| retained-fraction constants | `EMC2_LOG_CANCEL_MIN`, `EMC2_NAT_REL_CANCEL`, `EMC2_CDF_SAT_MARGIN` |

There is no `erfcx` in the tree and none is needed: `mills_ratio_std` covers
`x >= 0` by continued fraction and `x < 0` is well conditioned directly.

### Layer 1 — the stop-loss function (the only new primitive)

    C(v) = E[(V-v)_+] = M·Q(x-sigma) - v·Q(x),     x = (log v - mu)/sigma
    upper tail:  C(v) = v·φ(x)·[R(x-sigma) - R(x)]
    log form:    log C = log v + log_phi_std(x) + log(R(x-sigma) - R(x))

The upper-tail identity follows from `M φ(x-sigma) = v φ(x)`. Add
`log_lognormal_stoploss(v, mu, sigma)` next to `mills_ratio_std`, shaped like
`wald_pt_log_surv`: naive `log_diff_exp` of the two log terms while
well conditioned, Mills-difference regrouping past that, guard at
`EMC2_LOG_CANCEL_MIN`. **Build the `log φ` term analytically via `log_phi_std`,
never via `dnormP(x, log=TRUE)` routing through a natural `dnorm`** — that is
the difference between working and returning `NaN` in the deep tail.

Measured in `Math/check_stoploss.R` against a 400-bit `Rmpfr` reference: the
retained fraction of the subtraction is `≈ sigma/x`, and the naive **log**
difference still holds ~1e-11 relative on `log C` at a retained fraction of
`1e-7` (`sigma = 1e-4, x = 1e3`), degrading gracefully rather than failing. The
Mills form recovers full accuracy there (2e-15). Both forms stay finite
arbitrarily deep provided the `log φ` and `log Φ` terms are analytic; the naive
form's failure mode is lost digits, not `NaN`.

Consequence for effort allocation: the ordinary guarded ladder covers the
realistic range at layer 1 and the Mills branch is cheap insurance for the
small-`sigma` deep tail. Layer 3 is the one with a genuinely ill-conditioned
region that no rearrangement can fix.

### Layer 2 — endpoint differences

Never subtract two raw `C` or `D` values. Both intervals are needed only as
differences, and both are positive:

    log[C(a) - C(b)]  = log C(a) + log1m_exp(log C(b) - log C(a))     (a < b)
    log[D(x_a) - D(x_b)] = log_normal_q_interval(x_a, x_b)            (exists)

`log1m_exp` is the existing stable `log(1 - e^d)`. Same guard: if the retained
fraction falls below `1e-6`, hand off to layer 3's fallback.

### Layer 3 — the frozen combination `I_0 - ℓ·I_{-1}`

    log(I_0 - ℓ I_{-1}) = L_0 + log1m_exp(L_1 - L_0),  L_1 = log ℓ + log I_{-1}

The combination is non-negative because the critical launch values satisfy
`v ≥ ℓ`, so `I_0 - ℓ I_{-1} = ∫ Gbar(v)(1 - ℓ/v) dv ≥ 0`. It is
ill-conditioned exactly when the critical-launch interval concentrates near
`v = ℓ`, and no rearrangement of two separately evaluated terms recovers
relative precision there. Fall back to the **combined non-negative integral**,
which has no cancellation anywhere:

    I_frozen = (ℓ/k) ∫_{s_a}^{s_b} expm1(s) · Q(x_ell + s/sigma) ds,
    x_ell = (log ℓ - mu)/sigma

`expm1` keeps it accurate near `s = 0`; the integrand is smooth, positive,
bounded, and contains no Lambert-W. 8–12 GL nodes; trigger on
`retained = -expm1(L_1 - L_0) < 1e-6`.

Prefer this positive-integrand fallback to an asymptotic Mills series wherever
both apply: measured against the `Rmpfr` reference, the leading-order
`R(x-sigma) - R(x) ≈ sigma/(x(x-sigma))` is only 0.3–8% accurate over
`x ∈ [9, 30]`, whereas the quadrature is exact to node count for a comparable
handful of instructions.

Two structural notes. The **live** part and the **density** take the same
treatment: `C(w_lo) - C(w_hi) = ∫ Gbar` and `f·A = ∫ g(w)(w·E - ℓ) dw` are both
positive integrands (the density's vanishes at `w_lo` exactly when partially
saturated), so one positive-integrand GL routine serves the live, frozen and
density fallbacks. And the truncated-normal variant uses that same GL
scaffolding as its **ordinary** path, so the quadrature is written once and the
two variants differ only in when they reach for it.

---

## 3. C++ likelihood — `src/model_LBA.h`

Add after the shared BAwL core, which stays untouched. The new code is
structurally a copy of `bawl_threshold_terms` / `ba_natural_cdf` /
`ba_natural_pdf` / `log_ba_cdf` / `log_ba_pdf` with the new `c, m` and the `Z(u)`
upper limit, so all four existing acceptance modes
(`BA_ACCEPT_STRICT/RAW/CLAMP`) and the whole signed-log fallback ladder carry
over unchanged.

Factor it in two layers, since the launch distribution is the only difference
between 1a and 1b:

* *geometry* (shared, distribution-free): `q, E, w_hi, w_lo, Z, z*, y_A/s_A,
  y_0/s_0`, regime selection, `T_max`.
* *launch policy* (two implementations): the live integral, the frozen integral,
  and their two derivative counterparts.

Use a compile-time policy struct (`BawdNormalLaunch`, `BawdLognormalLaunch`)
with static methods: zero-cost, and it avoids a second copy of the regime logic,
which is where the subtle bugs would live. Each policy supplies exactly four
methods — `log_live`, `log_frozen`, `log_density_bracket`, and
`log_positive_fallback` (§2c layer 3) — so a third launch distribution later is
a new struct rather than a new branch in twelve places. The single
`bawd_launch` runtime branch sits at the top of each adapter, not inside the
per-row loops.

1. `bawd_geometry(u, b, A, k, ell, ...)` → the shared quantities above, plus
   flags `saturated`, `partially_saturated`, `ell_zero`, `k_zero`. Cache
   `y_A`/`y_0` in a small POD struct: they depend only on `(k, ell, b, A)`, not
   on `u`, so the truncation/censoring paths that call the CDF at several `t`
   for one row reuse them.
2. `bawd_newton_y(c)` → the `W_{-1}` solve of §2.
3. Frozen integral. **Normal launch:** `bawd_psi_natural(y1, y2, v, sv, k, ell)`
   by GL quadrature (`src/gl_quad.h`), with a log-space twin using
   `pnorm_log_direct` + log-sum-exp over nodes. **Lognormal launch:** the
   `I_0`/`I_{-1}` factoring of §2b plus `J(w)` for the live part, both in the
   cancellation-safe forms of §2c.
4. `ba_natural_cdf_bawd` / `ba_natural_pdf_bawd`, mirroring the existing
   signatures and acceptance semantics, including the `A ≤ BAWL_A_EPS` and
   `span < BAWL_NATURAL_MIN_SPAN` midpoint limits — the midpoint weight here is
   `-V*'(u, Z/2)/sv` in place of BAwL's `(b - A/2)`.
5. `log_bawd_cdf` / `log_bawd_pdf` — authoritative log path via
   `log_normal_phi_integral`, `log_normal_interval`, and `signed_log_*`. The
   density numerator is `(vE - ℓ)ΔΦ + sv·E·Δφ`, whose lead term is signed: it
   *does* go negative near `T_max`, unlike BAwL's `(bm + c)`, so the
   `BAWL_LOG_BRACKET_MIN` cancellation guard plus midpoint fallback is required,
   not optional.
6. `bawd_cdf_norm` / `bawd_pdf_norm` natural-then-log wrappers,
   `bawd_cdf_scalar_natural` / `bawd_pdf_scalar_natural` CLAMP variants for
   truncation normalisers, and exported vectorised `dbawd` / `pbawd` (plus
   scalar `dbawd_norm` / `pbawd_norm`) for the R `dfun`/`pfun`. Include
   `log_out` from the start; `pleakyba` lacks it and that is a wart worth not
   repeating.
7. `t == R_PosInf` must return `F_max`, not 1, and `F` must be *bit-identical*
   for every `t >= T_max` — evaluate the frozen branch, never extrapolate the
   live formula. This is the flat-CDF contract of §2, relied on by the omission
   and censoring paths (`particle_ll.cpp:8316`, `:8324`,
   `log_surv_cm(R_PosInf, ...)`).

No guess/kill clocks in v1: the model already generates omissions and a bounded
window mechanistically, and the clock cross product is where most of BAwL's
complexity lives. Leave the column spec extensible (`N_REQ` then optional
`mG`/`mK`) so they can be added later.

---

## 4. Wiring

### 4a. Selecting the launch distribution

A `c_name` suffix mapped to a `ctx` field reaches every consumer that needs it:

* `resolve_race_model_adapter()` has exactly **two** callers, `calc_ll_oo`
  (`particle_ll.cpp:3959`) and `calc_ll_oo_pw` (`:4398`); both then use
  `adapter.ctx` for the whole evaluation, so setting the field once in the
  adapter covers both entry points.
* The `ContextForRaceModels*` is threaded to every kernel signature: `RaceRawFun`
  / `RaceLogSAtTFun` (`void* ctx_`), `RacePdf1Fun` / `RaceCdf1Fun` (`void* ctx`),
  `gsl_race_params_scalar::ctx`, and the `model_ctx` /
  `model_context_for_funcs` pointers used by the truncation and censoring
  branches (`:4663`, `:4700`, `:5440`, `:5482`, `:5981`, `:7040`, `:7596`) —
  roughly 100 references, with no evaluation path lacking it.
* Existing precedent: `pcounter_integer_K` from `"INTK"`, `bawl_correlated` from
  `"_CORR"`, `kill_shape` from `"_E2"` / `"_EMIX"`, ROU's `bnd_kind` from
  `"_BWEIB"` / `"_BEXP"` / `"_BLIN"`.

Add `int bawd_launch = 0;  // 0 = truncated normal, 1 = lognormal` to
`ContextForRaceModels` and set it from the suffix. Three constraints:

1. **Suffix collision.** Dispatch is `std::string::find()`, so the suffix must
   not be a substring of any other key. Use **`_LOGN`**. Do *not* use `_LNR`,
   which collides with the LNR branch, and do not rely on branch ordering to
   avoid it; `_LN` is safe today but `"LNR"` shows how close the namespace is.
2. **The R-callable exports bypass `ctx` entirely.** `dbawd` / `pbawd` (reached
   from R via `dBAwD`/`pBAwD`) and `rbawd_cpp` are plain `[[Rcpp::export]]`
   functions with no context argument, so they need an explicit `int launch = 0`
   parameter, the way `posdrift`, `log_out` and `kill_shape` are passed today.
   The R constructor closure must derive it from the *same*
   `drift_distribution` value it puts in `c_name`, in **one** place, or
   `dfun`/`pfun` will silently disagree with the sampled likelihood. This is the
   failure `.rou_cols()` was written to prevent; assert agreement in a test
   rather than trusting it.
3. **Two ColSpecs, because `p_types` differ** (`v, sv` vs `mu, sigma`) and
   `validate_col_prefix()` checks names positionally. Keep the *positions*
   identical so all kernel indexing is shared, and select the spec in the same
   adapter branch.

### 4b. File-by-file

**`src/col_registry.h`** — two namespaces next to `bawl`, identical layout:

```cpp
namespace bawd {                       // BAwD / BAwDIO — truncated-normal launch
  enum : int { v = 0, sv, B, A, t0, k, ell, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v","sv","B","A","t0","k","ell"};
    return {n, N_REQ, "BAwD"};
  }
}
namespace bawd_logn {                  // BAwD_LOGN — lognormal launch
  enum : int { mu = 0, sigma, B, A, t0, k, ell, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu","sigma","B","A","t0","k","ell"};
    return {n, N_REQ, "BAwD_LOGN"};
  }
}
```
The enum order **must** equal `p_types` order in the constructor;
`validate_col_prefix()` enforces it per call.

**`src/utils.h`** — scalar/raw/batch adapters modelled on `dbawl_*`:
`dbawd_scalar`, `pbawd_scalar` (`RacePdf1Fun`/`RaceCdf1Fun`, natural scale,
CLAMP evaluators), `dbawd_raw`, `pbawd_raw` (`RaceRawFun`, writing log-density /
**log-survivor**, floored via `raw_log_value`/`raw_log_zero`), and
`bawd_logS_at_t` (`RaceLogSAtTFun`). Kernel contract to respect: `pfun_raw`
returns the log-survivor, not the CDF. The `bawl_k_fixed_zero` /
`bawl_clocks_fixed_off` fast-path pattern is a useful model but those flags must
not be reused — BAwD needs no per-row optional-column gating.

**`src/particle_ll.cpp`, `resolve_race_model_adapter()`** — new branch
immediately before the `BAwL` branch. Dispatch is by substring and `"BAwD"` is a
substring of nothing and contains neither `"BAwL"` nor `"LBA"`, so placement is
safe either way; keep it adjacent for readability and say so in a comment.

```cpp
} else if (type_std.find("BAwD") != std::string::npos) {
  const bool logn = type_std.find("_LOGN") != std::string::npos;
  out.pdf1_ptr = &dbawd_scalar;  out.cdf1_ptr = &pbawd_scalar;
  out.model_dfun_raw = &dbawd_raw; out.model_pfun_raw = &pbawd_raw;
  out.logS_at_t_ptr = &bawd_logS_at_t;
  out.col_spec = logn ? emc2col::bawd_logn::spec() : emc2col::bawd::spec();
  out.ctx.t0_index = emc2col::bawd::t0;   // same position in both layouts
  out.ctx.bawd_launch = logn ? 1 : 0;
  out.ctx.defective_upper_tail = true;    // ALWAYS, not only under IO
  // posdrift is meaningless for the lognormal launch (V > 0 by construction).
  if (!logn && type_std.find("IO") != std::string::npos)
    out.ctx.use_posdrift = false;
}
```
The kernels then branch once on `ctx->bawd_launch` to select the launch policy;
the geometry call above it is shared. `kill_active` already resolves to `false`
for an unsuffixed name; leave `mean_g_index`/`mean_k_index` at `-1`.

**`R/model_LBA.R`** — `dBAwD`, `pBAwD`, `rBAwD` (R reference simulator), and the
exported constructor

```r
BAwD(drift_distribution = c("lognormal", "normal"), posdrift = TRUE)
```

with `type = "RACE"`, `c_name` `"BAwD_LOGN"` / `"BAwD"` / `"BAwDIO"`,
`p_types_canonical` matching the chosen launch form, `Ttransform` adding
`b = B + A`, and a roxygen block covering the parameter table, the BAwL/LBA
nesting for `drift_distribution = "normal"`, the flat-CDF and bounded-support
behaviour, and the scale-fixing rule of §1a/§1b. `bound$minmax`:
`ell = c(1e-4, Inf)` with `exception = c(A = 0, k = 0, ell = 0)` so `ell = 0`
(the static-start BAwL limit) stays exactly reachable. `posdrift` applies only
to the normal launch; raise an error if it is passed as `FALSE` alongside the
lognormal form rather than silently ignoring it.

`log_likelihood` is a `stop()` pointing at the compiled path, as
`LogicalRulesLBA` and `BAwLcorr` already do. The R likelihood is not a supported
route in this project, so `R/likelihood.R` needs no change and
`log_likelihood_race_missing`'s `is_defective` handling is not in scope.

**`R/make_data.R:825`** — `posdrift <- !grepl("IO", model()$c_name)` works for
`BAwDIO` unchanged; confirm by test rather than by inspection. `dfun`/`pfun`
*are* live (posterior prediction, `make_data`, diagnostics) and must call the
same compiled kernels as the likelihood, following the ROU pattern.

**`NAMESPACE` / `man/`** — `export(BAwD)` via roxygen, regenerate docs.
Regenerate `src/RcppExports.*` for the new `[[Rcpp::export]]` entries.

---

## 5. Simulation

Given a launch strength `V` (from whichever distribution is configured) and
`z ~ U(0,A)`:

1. If `V ≤ ell` → `Inf`. Else peak `u_p = log(V/ell)/k`,
   `X_max = V q(u_p) - ell·u_p`; if `X_max < b - z` → `Inf`.
2. Else Newton/bisect `V q(u) - ell·u - (b - z) = 0` on `(0, u_p]`. The function
   is concave increasing there, so Newton from `u_p` is monotone and safe; 5–8
   iterations to 1e-12. The `W_{-1}` closed form is equivalent; Newton avoids a
   second branch-selection site.
3. `k ≤ eps` → `(b - z)/(V - ell)`; `ell ≤ eps` → BAwL's inversion
   `-(1/k) log((V - k(b-z))/V)` when `V > k(b-z)`, else `Inf`.
4. Add `t0`, race by row-min per trial. An all-`Inf` column ⇒ `R = NA`,
   `rt = Inf`, the package's existing omission convention (`bad_col` in `rBAwL`,
   `resolve_race`/`pack_result` in `model_rng.cpp`, and the GNG withheld
   convention). Nothing downstream needs changing; `make_data()` already handles
   `rt = Inf` with `R = NA`.

Ship the R reference `rBAwD` **and** `rbawd_cpp` in `src/model_rng.cpp` with
`.rfun_BAwD` in `R/model_rng.R` following `.rfun_BAwL` (C++ under
`.use_cpp_rfun()`, R fallback otherwise). The pair is what makes the R-vs-C++
equivalence test possible, and posterior prediction on a defective model draws a
lot of samples.

---

## 6. Tests — `tests/testthat/test-bawd.R`

Organised like `test-rou.R`, around known oracles first.

1. **`ell = 0, A = 0` ≡ BAwL** (normal launch): `dBAwD`/`pBAwD` vs
   `dleakyba`/`pleakyba`, exact to 1e-12. This is the anchor test; it validates
   the shared geometry layer for both variants.
2. **`ell = 0, A > 0` ≠ BAwL**: assert a material difference, with a comment
   explaining why (static vs decayed start point) so a later "simplification"
   cannot quietly collapse the two.
3. **`k → 0` ≡ LBA with drift `v - ell`**: vs `dlba`/`plba`, 1e-8.
4. **Monte Carlo** (`Math/check_bd.R` and `Math/check_lnorm.R` promoted):
   pdf/cdf/`P(T=∞)` for **both** launch forms over ≥3 parameter sets each,
   spanning unsaturated / partially saturated / fully saturated `u`, `A = 0` and
   `A > 0`, `posdrift` both ways.
5. **The flat-CDF contract** of §2: `f(u) = 0` exactly for `u >= T_max`; `F`
   bit-identical at `T_max`, `2·T_max`, `1e6` and `Inf`; `S(Inf) = 1 - F_max`
   equal to the simulator's omission rate; `f ∝ (T_max - u)` at `A = 0` and
   `∝ (T_max - u)²` for `A > 0` (the live interval collapses linearly *and* the
   integrand vanishes at its lower end) just
   below the endpoint.
6. **`d`/`p` consistency across the seams**: `integrate(f) == F` — good to ~9
   digits for the lognormal form, so use a tight tolerance there — and `pBAwD`
   monotone non-decreasing through both seams (`T_max(A)` and `T_max`), where a
   live/frozen mismatch would surface as a kink.
7. **Compiled likelihood**: the C++ race path against an independent R reference
   written in the test file (the `test-rou.R` `wald_pdf` pattern) on a small
   `design()` — omission trials, `LT/UT` truncation, `LC/UC` censoring with `UC`
   deliberately placed **beyond `T_max`**, and `pContaminant`. The
   defective-tail branches at `particle_ll.cpp:8292–8380` are what needs
   coverage; the beyond-`T_max` censoring case is what would break if an
   evaluator clamped `F` to 1.
8. **rfun**: `.rfun_BAwD` C++ vs R reference under a fixed seed; recovery smoke
   test that `make_emc`/`fit` runs a few iterations.
9. **Launch-form agreement**: `dBAwD`/`pBAwD` (R, via the exported `launch`
   argument) vs the compiled likelihood for *both* forms on the same design.
   This is the §4a constraint 2 failure mode and it is silent if untested.
10. **Cancellation layers** (§2c), the part most likely to rot:
    * `log_lognormal_stoploss` vs an `Rmpfr` reference across
      `sigma ∈ {0.5, 1e-2, 1e-4} × x ∈ {10, 1e3, 1e5}`, the grid in
      `Math/check_stoploss.R`; assert finite, no `NaN`, and agreement to the
      tolerance each branch claims. Note that `Phi(-1e5)` is outside `Rmpfr`'s
      default exponent range, so those rows bound the reference, not the code.
    * layer 3 forced into its fallback (critical-launch interval concentrated at
      `v ≈ ℓ`): the closed and positive-integrand routes must agree, and the
      trigger must actually fire — assert on a branch counter, not just the
      value, or a later refactor can quietly make the fallback dead code.
11. **Support constraint**: a finite RT above `t0 + T_max` must reach the
    `min_ll` floor cleanly, with no raw `-Inf` or `NaN` escaping.

Log-space coverage: extend the `test-lba-logspace.R` style of extreme-parameter
sweep (tiny `sv`/`sigma`, large `|v|`/`|mu|`, `u` within 1e-9 of `T_max`, `w`
deep in the lognormal upper tail where `C` cancels), asserting that the natural
and log branches agree where both are valid and that neither returns `NaN`.

---

## 7. Order of work

One pass covering both distributions. The validation order within it is fixed:
the truncated-normal launch is the only variant with 1e-12 oracles, so it
certifies the shared geometry before the lognormal policy is built on top.

1. Shared geometry (`bawd_geometry`, `bawd_newton_y`, regime selection), the
   positive-integrand GL routine of §2c layer 3, and the **normal-launch**
   policy; exported `dbawd`/`pbawd` with the `launch` argument. Validate against
   `Math/check_bd.R` and the exact BAwL/LBA corners. Everything after this
   stands on this geometry, so do not proceed while it is red.
2. `log_lognormal_stoploss` in `wald_functions.h` next to `mills_ratio_std`,
   with its own `Rmpfr`-referenced test. This is the only new numerical
   primitive in the change; `D` / `I_{-1}` reuse `log_normal_q_interval`
   unchanged.
3. **Lognormal-launch** policy on the certified geometry: the `I_0`/`I_{-1}`
   factoring, the partial-expectation density, and the layer-2 and layer-3
   guards. Validate against `Math/check_lnorm.R`, `integrate(f) == F` to ~1e-9,
   and direct quadrature of `Gbar` for live and frozen **separately** —
   separately, because a compensating error in the pair passes the combined
   check.
4. `col_registry.h` (both specs), `utils.h` adapters, the single
   `particle_ll.cpp` branch, and `ctx.bawd_launch`.
5. `R/model_LBA.R` constructor with `drift_distribution`, `dBAwD`/`pBAwD`
   threading `launch` from one place, roxygen.
6. `rBAwD` (R, both forms), then `rbawd_cpp` and `.rfun_BAwD`.
7. `test-bawd.R` in the order of §6.
8. Regenerate `RcppExports`, `NAMESPACE`, `man/`; full `testthat` run.

---

## 8. As built

Implemented on branch `bawl-correlated-race`. Where the delivered code departs
from the plan above, this is the record.

**File layout.** The numerics are in a new `src/model_BAwD.h`, included at the
end of `src/model_LBA.h` (so every guard constant it uses is already defined,
and it remains confined to the one translation unit that may hold
`[[Rcpp::export]]` definitions). The precedent is `src/bawl_geometry.h`, which
does the same. The R side is `R/model_BAwD.R` rather than an addition to
`R/model_LBA.R`, for the same reason: 700 lines of a different process do not
belong in the LBA file.

**One evaluation path, not two.** §3 anticipated BAwL's structure — a guarded
natural-space evaluator plus an authoritative log-space fallback. Only the log
path was built; `bawd_cdf_norm` / `bawd_pdf_norm` and the
`*_scalar_natural` evaluators exponentiate at the boundary. The reason is that
BAwD has three time regimes and two launch distributions, and a natural-space
duplicate of that dispatch is where the subtle bugs would live. It costs a few
`log` calls per evaluation against a model whose default variant is closed form.
The log primitives used are the ones BAwL's own log branch uses
(`log_normal_phi_integral`, `log_normal_interval`, `log_normal_q_interval`),
plus the one new `log_lognormal_stoploss`.

**`log_lognormal_stoploss` never returns `NA`.** §2c layer 1 proposed a
retained-fraction guard handing off to quadrature. Measurement showed the Mills
form good to 2e-15 well past the proposed trigger, so the guard sits at a
retained fraction of `1e-12` and below it the function uses the *symbolically
differenced* asymptotic
`R(y) - R(x) = sigma/(xy) - sigma(x²+xy+y²)/(x³y³) + ...`, whose relative error
is `O(1/x²)` and therefore negligible wherever that branch is reachable. Layer
3's positive-integrand quadrature remains, for the frozen combination
`I_0 - ell·I_{-1}`, which does have a genuinely ill-conditioned region.

**`ell` defaults to `log(1)`, not `log(0)`.** Leaving `ell` out of the `design()`
formula then fixes the evidence scale, which is exactly the convention §1b
requires for the lognormal launch. `ell = 0` remains reachable as a bound
exception.

**Endpoint order.** §0 and §6 originally claimed `f ∝ (T_max - t)`. That holds
only for a point start. With `A > 0` the live interval `[w_lo, w_hi]` collapses
linearly *and* the integrand vanishes at its lower end, so `f ∝ (T_max - t)²`.
Both orders are asserted in `test-bawd.R`.

**Saturation test.** Saturation is `z*(u) <= 0`, **not** the clamped
`Z(u) <= 0`. With a point start `Z` is identically zero, so the latter reports
every time as saturated and the CDF collapses to the constant `F_max`. This was
a real bug in the first build, caught by the `A = 0` reference comparison.

**Extra exports for testing.** `bawd_tmax()` (the right endpoint) and
`lognormal_stoploss_log()`, because both the flat-CDF contract and the
cancellation layers are stated in terms of quantities a test would otherwise
have to recompute — and a test that recomputes them is not testing them.

**Simulator.** `bawd_hit_time_r()` in `src/model_rng.h` inverts the trajectory by
Newton from the LBA-limit time `d/(V - ell)`, which lies below the root because
`q(u) <= u`; the trajectory is concave and increasing up to its peak, so the
iteration rises monotonically and cannot overshoot. `R/model_BAwD.R` mirrors it
for the pure-R reference path. The `W_{-1}` closed form is equivalent but needs
a branch choice at each call site.

**`e^{-ku}` is carried in log space.** `E = e^{-ku}` underflows to zero for
`ku > 745` while the density it multiplies is still representable, so every term
with a bare factor of `E` is formed as `log(coefficient) + log_E`. For the
truncated-normal launch that means regrouping the numerator as
`E(v ΔΦ + sv Δφ) - ℓ ΔΦ` rather than `(vE - ℓ)ΔΦ + sv E Δφ` — algebraically the
same, but with the single factor of `E` outside the signed-log sum. `BawdAtU`
therefore carries `log_E` beside `E`.

**Deep-tail oracle.** `Rmpfr`'s `pnorm` underflows to zero past `x` of a few
thousand, so the definition-based reference `M Q(x-sigma) - v Q(x)` cannot check
the deepest rows. The test adds a second 400-bit reference built from a
200-level Mills-ratio continued fraction, which stays `O(1/x)` and so reaches
arbitrarily deep — and is algorithmically independent of the kernel's own
asymptotic branch.

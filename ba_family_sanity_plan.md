# Ballistic Accumulator Family: Consistency and Correctness Plan

## 1. Context

A sanity pass over the six ballistic accumulator models (LBA, BAwL, BAwD/BAwDp,
BAwF, BAwR, BTAwL) found that BAwD/BAwF/BAwR share a clean, well-factored
template, and that closed-form coverage is already correct: the places that fall
back to quadrature (the normal-launch frozen term everywhere, and the lognormal
one when the critical-launch interval collapses) genuinely have no elementary
closed form. `src/col_registry.h` is consistent across the family, with identical
column *positions* and only names differing between launch laws.

Three things need fixing.

### 1.1 No BA model has an accurate log-survivor

All of them obtain `log(1-F)` by complementing the CDF, so the survivor saturates
near `log(eps) = -36` and then returns `-Inf` or non-monotone noise. Measured
against an exact, cancellation-free reference (integrating the launch *CDF* over
start points, so no subtraction occurs anywhere):

```
plain LBA           BAwL (k=2)           BAwF (lognormal)
v=6   ref  -62.1    v=8   ref  -97.4     mu=5   ref  -64.2
      got  -36.0          got   -Inf            got   -Inf
v=14  ref -699.5    v=10  ref -209.9     mu=12  ref -591.5
      got  -36.0          got  -36.0            got  -22.4
```

BAwD and BAwR show the same non-monotone junk by `mu = 5`. On the real likelihood
path (30-trial LBA, sweeping the winner's drift) the log-likelihood falls
correctly to -1200 and then goes flat and non-monotone at `60 x min_ll`:

```
v_lMTRUE:   1      2      3      4      6       8      12      20      40      80
logLik:  +3.58  -15.8  -88.7 -217.4 -645.6 -1200.7 -1380.9 -1381.0 -1381.4 -1381.2
```

This is the survivor of every losing accumulator, so it fires whenever one
accumulator is confidently faster than the observed RT: exactly where the sampler
needs a wall rather than a plateau.

### 1.2 BTAwL was written against a different template

No log-space path (it is `log()` of a `[0,1]`-clamped natural value), no
`BA_ACCEPT_*` contract, four duplicated helpers, and two avoidable root solves in
the hot path. Measured per call, same harness, each model across its own support:

```
BTAwL cdf   9.80 us     BAwF cdf 0.40    BAwR cdf 0.25    BAwD cdf 0.62
BTAwL pdf   2.45 us     BAwF pdf 0.27    BAwL cdf 0.28    LBA  cdf 0.18
BTAwLmix  111.50 us
```

Broken down: `btawl_tmax` alone is 2.48 us (versus 0.14 us for `bawf_tmax_vec`)
and runs on *every* evaluation; the frozen branch adds ~14 us in the late-`t`
region.

### 1.3 `T_max` is a derived diagnostic where BAwD now samples it

For BAwR the endpoint mixes `b`, `kappa` and `p`, so the rate chart puts the
well-determined coordinate across an axis rather than along one.

**Outcome:** all six models numerically sound in both tails, BTAwL structurally
identical to its siblings and roughly an order of magnitude faster, and the
endpoint on a sampled axis wherever the model has one.

---

## 2. Working-tree prerequisite

The tree is **not clean**. Unrelated in-flight ROUp work touches
`R/model_ROUp.R`, `R/model_rng.R`, `src/model_ROUp.h`, `src/rou_diffusion.cpp`,
`src/fpe_race.h`, `tests/testthat/test-roup.R`, `man/ROUp.Rd`, and carries
`NAMESPACE` / `R/RcppExports.R` / `src/RcppExports.cpp` / `src/col_registry.h` /
`src/particle_ll.cpp` edits with it. Those last five are exactly the files
Phases 2 and 5 must modify, and both phases regenerate RcppExports.

**Settle this before Phase 1**: commit the ROUp work, stash it, or branch from
before it. Then branch for this work.

---

## 3. Phase 0 — Port the test battery to BTAwL

**Do this first, before touching any BTAwL code.** `tests/testthat/test-btawl.R`
is 103 lines / 7 tests against `test-bawf.R`'s 544 / 21 and `test-bawr.R`'s
557 / 21. Every later phase is unguarded until this exists, and it costs nothing
numerically.

Mirror `tests/testthat/test-bawf.R`, reusing its helper shape (`bawf_dat`,
`bawf_mk`, `bawf_ll`, `ref_F_bawf`, `ref_f_bawf`, `ref_race_ll_bawf`). Add for
BTAwL:

- CDF and PDF against independent start-point quadrature
- `integrate(f) == F` across the live/frozen seam
- Monte Carlo agreement for both launch distributions
- log output equals log of natural output
- CDF monotone, density non-negative across both seams
- `F(Inf)` is the finishing mass and the tail is defective
- endpoint shutdown quadratic for `A > 0`, linear for `A = 0`
- compiled race likelihood vs the R reference
- `p_types` reordering rejected by the column contract
- omissions past `t0 + T_max` well posed
- R simulator agrees with the CDF, and with `make_data` round-tripped through a fit
- `T_max = 2*tau` at `k*tau = 1` (an exact identity, and a good chart check)

The family's "k = 0 is exactly the LBA" test deliberately does **not** apply:
BTAwL's transient shape has no off switch, so `k = 0` is a finite-total-drive
model with no hard endpoint, not an LBA.

Also add a **stiff-window** case. At `k*tau = 20` the frozen range is extreme —
at `k=20, tau=1`: `u1 = 1.05258`, `T_max = 1.05263` (a window 5e-5 wide) with
`max Gamma = 1.3e6`. That regime must be in the suite before Phase 4b changes the
quadrature.

---

## 4. Phase 1 — Minor consistency fixes

Cheap, low risk, no intended numerical change.

- **`src/model_BTAwL.h:123-125`** — the comment claims the tangency start point is
  "monotone decreasing from b at onset to zero at t_max". It is not. `Gamma`
  starts at 1, rises to a maximum at exactly `u = tau`, then falls to 0 at
  `T_max` (verified across `k` in [0.3, 20] and `tau` in [0.05, 1]; at
  `k=2, tau=0.2` the peak is 1.124). The code is correct **only because** it
  clamps to `[0, b]` at `:132`, and because
  `d(u) = E*eta(u)*[z - b*Gamma(u)]` still has exactly one sign change for
  `z < b`. Rewrite the comment and say the clamp is load-bearing — a future
  Newton solve seeded from an assumed monotonicity would be wrong.
  `Math/BallisticAccumulators.tex:2076` already describes this correctly.
- **`R/model_BTAwL.R:301-304, 372-375`** — `log_likelihood` routes to
  `log_likelihood_race_missing`; BAwD/BAwF/BAwR all `stop()` because the R
  likelihood route is outdated. Match them.
- **`R/model_BTAwL.R:285, 356`** — hardcoded `p_types_canonical`; use
  `setdiff(names(p_types), .nuisance_par_names)` as everywhere else.
- **`R/model_BTAwL.R:274, 340-341`** — `k = c(0, Inf)`; use `c(1e-4, Inf)`. The
  `exception = c(k = 0)` already carries the limit, as in BAwF/BAwR.
- **`src/model_BTAwL.h:22-26`** — alias `BTAWL_LAUNCH_*`, `BTAWL_K_EPS`,
  `BTAWL_A_EPS` to the `BAWL_*` definitions, the way `src/model_BAwF.h:54-64`
  aliases BAwD's. Removes a fourth independent definition of the 0/1 launch pair.
- **`src/utils.h:136`** — `ctx.btawl_launch` defaults to LOGNORMAL while BTAwL's R
  default is `"normal"`; `bawl_launch` (`utils.h:150`) explicitly defaults to
  match its model. Flip it, and align the `dbtawl`/`pbtawl`/`dbtawlmix`/
  `pbtawlmix` Rcpp defaults (`launch = 1`) with the R wrappers' `launch = 0L`.
  Re-run `Rcpp::compileAttributes()`.
- **`src/model_BAwR.h`** — 35 `Bawp*` struct/type names left from the BAwP rename;
  rename to `Bawr*`. Mechanical, confined to that header.
- **`NAMESPACE`** — `BTAwLMix` and `BTAwL_mixed` are both exported and
  `BTAwLTransient` is defined but not. Checked: no test references the
  `BTAwLMix` *constructor* (`test-btawl.R` calls `EMC2:::pBTAwLMix`/`dBTAwLMix`,
  which are different functions), so keep `BTAwL` + `BTAwL_mixed` exported and
  demote the two aliases.

### 4.1 Trap: the denominator floor

`btawl_normal_denom` (`src/model_BTAwL.h:152-156`) floors at
`BTAWL_DENOM_FLOOR = BAWL_DENOM_FLOOR = 1e-300`. Both
`log_positive_normalizer` (`src/model_LBA.h:62-63`) and `natural_normalizer`
(`:76-78`) default to `LBA_DENOM_FLOOR = 1e-10`. `src/model_LBA.h:36-42` states
in so many words that this difference is a **model contract**, not an
implementation detail, and that changing it moves BAwL values in the
`z < about -37` regime.

So when deleting `btawl_normal_denom`, pass `BAWL_DENOM_FLOOR` **explicitly at
every call site**. Add a regression test pinning `v/sv` near -40 with
`posdrift = TRUE`.

For the record, `btawl_surv` at `w <= 0` (`:160-164`) is **correct** and should
be preserved when it is folded into the shared idiom: untruncated it returns
`P(V >= 0) = Phi(v/sv)`, truncated it returns `Phi(v/sv)/Phi(v/sv) = 1`. The
apparent asymmetry is the truncation normaliser doing its job.

**Verify:** `R CMD INSTALL` into a temp lib, full `testthat` run, and confirm
likelihood values are bit-identical to the pre-change build on a saved fit.
**Commit.**

---

## 5. Phase 2 — BAwR endpoint chart

Replace `kappa` with `Tmax` as the sampled coordinate, mirroring BAwD exactly.
`T_max = (b(p+1)/(kappa*p))^(1/(p+1))` was verified against `bawr_tmax_vec` at 16
parameter combinations with relative error 0, and the inverse

```
kappa = b (p+1) / (p * Tmax^(p+1))
```

is closed-form and an **unconstrained bijection** onto `(0, Inf)` for every
`b > 0, p > 0`. Simpler than BAwD in one respect: `bawr_geometry` takes `kappa`
directly, so there is no circularity to break — BAwD needed `bawd_shape_flags`
(`src/model_BAwD.h:189-220`) only because its inverse had to evaluate the shape
function without the `ell` it was computing.

- **`src/model_BAwR.h`** — add `bawr_kappa_from_tmax(b, pw, Tmax)` next to
  `bawr_sat_time`. `Tmax = Inf` returns `0.0` (the inert `kappa = 0` LBA limit,
  matching `bawd_ell_from_tmax`); NaN inputs return `R_NaN` so a bad value
  propagates rather than being swallowed. Export `bawr_kappa_vec(Tmax, b, pw)`
  alongside `bawr_tmax_vec` / `bawr_vcrit_vec`.
- **`src/utils.h`** — one `bawr_clear_to_kappa(clear, b, pw)` helper mirroring
  `bawd_clear_to_ell` (`src/utils.h:1557-1564`), used by *every* column read in
  `dbawr_scalar`, `pbawr_scalar`, `dbawr_raw`, `pbawr_raw`, `bawr_logS_at_t`, so
  no kernel can forget the back-solve. BAwR has a single chart, so unlike BAwD no
  `uses_tmax` predicate is needed.
- **`src/col_registry.h:166-180`** — rename slot 5 from `kappa` to `Tmax` in both
  `bawr::spec()` and `bawr_logn::spec()`, and rename the enum member to `clear`
  (BAwD's convention). No second spec, no c_name suffix. `validate_col_prefix()`
  fails loudly on any mismatch.
- **`R/model_BAwR.R:257-267`** — `p_types` gains `Tmax = Inf`, drops `kappa`;
  `transform` exp; `minmax = c(1e-4, Inf)`; `exception = c(A = 0, Tmax = Inf)`.
- **`R/model_BAwR.R:287-297` `Ttransform`** — derive `kappa` via `bawr_kappa_vec`
  and `cbind` it on, so `rBAwR` / `.bawr_hit_time` / `dBAwR` / `pBAwR` /
  `mapped_pars()` keep seeing the mechanistic rate and need no change at all
  (this is exactly why BAwD's simulator was untouched). Then recompute `Tmax` via
  `bawr_tmax_vec` and **overwrite** the sampled column so the reported endpoint
  is the realised one, as `R/model_BAwD.R:530-536` does. Keep `Vcrit`, `b`,
  `rt_max`.
- **`man/BAwR.Rd`** — restate BAwD's gauge caveat: `kappa` is now derived and can
  no longer be held at a constant to fix the evidence scale, so pin an intercept
  of `mu` (or `v`), `B` or `A`.

**Commit.**

---

## 6. Phase 3 — Accurate log-survivors across the family

The adapter needs **no** surgery: `model_pfun_raw`'s contract *is already*
`log(1-F)` (`src/utils.h:1542-1545`). The defect is purely that every BA
implementation produces it as `log1m_exp(log_cdf)`.

### The fix

`S(u) = (1/A) * integral_0^A G(v_req(z)) dz` is the same integral as the CDF with
the launch survivor replaced by the launch **CDF**. Every normal factor becomes a
*lower* tail, which `pnorm_log_direct(x, true)` evaluates accurately past -1700.
Add alongside each `log_<model>_cdf_*`:

- `log_<model>_surv_normal(u, g, v, sv, posdrift, denom_floor)`
- `log_<model>_surv_logn(u, g, mu, sigma)`

reusing `log_normal_phi_integral` (its argument order reversed already gives the
complementary `integral Phi`), the `signed_log` helpers, and `bawd_log_gl_split`
for the frozen term with `pnorm_log_direct(..., true)` substituted for
`(..., false)`. Write the lognormal "put" primitives **directly** rather than
deriving them from `log_lognormal_*_stoploss` by put-call parity, which would
reintroduce the very cancellation being fixed.

### Wiring

In `src/utils.h`, change `pbawl_raw`, `pbawd_raw`, `pbawdp_raw`, `pbawf_raw`,
`pbawr_raw` and their five `*_logS_at_t` partners to keep the existing
natural-space fast path in front — it is accurate while `cdf < 1 - 1e-8`, i.e.
`logS > -18.4`, which covers the overwhelming majority of evaluations and
preserves current values there — and route to the new evaluator *instead of*
`log1m_exp(log_cdf)` when it is rejected.

Scope is the BA family only. BTAwL's survivor arrives with its log path in
Phase 4. The other 13 race kernels are out of scope for this pass.

### Deferred to 3b — decide after measuring

The *scalar* path (`log_survivor_rowmajor`, `src/particle_ll.cpp:5624-5658`) goes
through `cdf1_ptr`, which returns a natural CDF, and applies `safe_log1m_race`.
Fixing it needs a new `logS1_ptr` adapter slot with a nullptr default and a
fallback, plus assignments in the six BA branches of
`resolve_race_model_adapter`. That path feeds censoring/omission (`log_surv_cm`,
~8727), the truncation-normaliser fallback
(`get_trunc_normaliser_rowmajor_cpp`, 5875-5991), LogicalRules, and the
`defective_finite_UT` route. Measure how far it actually binds first.

**Verify:** the exact-survivor reference becomes a testthat helper applied to all
five models; the `v_lMTRUE` likelihood sweep becomes a monotonicity regression
test. **Commit.**

---

## 7. Phase 4 — BTAwL numerical overhaul

Conform BTAwL to the existing `XxxGeom` / `XxxAtU` / `BA_ACCEPT_*` template. Do
**not** extract a shared header yet — see §11.

### 4a. `btawl_tmax` Newton — no intended value change

`src/model_BTAwL.h:76-90` does 80 doubling probes plus 100 bisections, called by
`btawl_geometry` on every evaluation: 2.48 us, versus 0.14 us for
`bawf_tmax_vec`. Replace with a bracketed Newton on `phi(u) = eta(u) - k*K(u)`,
seeded at `tau*(1 + 1/(k*tau))` and bracketed **below by `tau`**, since
`T_max > tau` is proved in the .tex. Assert 1e-12 agreement with the old
bisection across a `(k, tau)` grid before deleting it.

Do this early: it is isolated, it is a large win, and it makes every later test
run faster.

### 4b. Frozen reparametrisation z -> s — changes values at the ~1e-7 qags level

Replaces `btawl_frozen_cdf` (`:259-283`) and kills `btawl_tangent_time`
(`:138-150`) and `btawl_frozen_integrand` (`:252`). Reparametrise to tangency
time `s`, where the Jacobian is closed form:

```
|dz/ds| = b * k * K(s) * |eta'(s)| / (exp(-k*s) * eta(s)^2)
```

verified against numerical differentiation to 1e-9, and the resulting frozen mass
verified against the z-coordinate reference. Integrate with `bawd_log_gl_split`
in log space, exactly as `bawf_log_frozen_quad` and `bawr_log_frozen_quad` do.
No inner root solve.

The lower limit is the **descending**-branch root of `b*Gamma(s) = A`; nothing
saturates before `u1`, where `Gamma` falls back through 1.

> **Checked:** the `|eta'(s)|` factor vanishes at `s = tau`, but that point is
> *outside* the integration range. `argmax Gamma = tau` exactly, and `u1 > tau`
> at every `(k, tau)` combination tested, so `(s_A, T_max)` lies entirely on the
> descending branch. There is no interior kink.

Do **not** assume BAwF's node count and single split transfer. Verify against the
retained qags reference, especially in the stiff large-`k*tau` window from
Phase 0. Retain the old `btawl_frozen_cdf` as a test-only `[[Rcpp::export]]` for
exactly one commit, assert ~1e-8 relative agreement, then delete.

**Must land alone.**

### 4c-4d. Log-space path, `BA_ACCEPT_*`, and `BtawlAtU` — must land together

The `BA_ACCEPT_*` contract's entire content is "return `false` so the caller uses
the log path" (`src/model_LBA.h:51-60`), so shipping the guards without a log
path just converts accepted values into zeros. Introducing `BtawlAtU` forces the
live CDF and live PDF to move together.

Deliverables, mirroring BAwR: `btawl_at_u`, `log_btawl_cdf_{normal,logn}`,
`log_btawl_pdf_{normal,logn}`, `log_btawl_surv_{normal,logn}`,
`btawl_natural_{cdf,pdf}_{normal,logn}`, `ba_natural_{cdf,pdf}_btawl`,
`btawl_log_{cdf,pdf}`, `btawl_{cdf,pdf}_norm`, `btawl_{cdf,pdf}_scalar_natural`.
Delete `btawl_cdf_log` / `btawl_pdf_log` (`:325-334`) and the
`fmin(fmax(out,0),1)` clamp at `:305` that currently destroys both tails.

Delete the duplicated helpers: `btawl_normal_denom` (subject to §4.1),
`btawl_J` -> `log_normal_phi_integral` / `normal_phi_integral_nat`,
`btawl_surv` / `btawl_pdf_v` -> the inline `log_gbar` lambda idiom. Note
`btawl_surv` / `btawl_pdf_v` are **also used by `BTAwL_mixed`**, so either keep
them alive until Phase 6 or sweep the mixture in the same commit.

Two structural details to carry over deliberately:

- **`w_lo` must be *set*, not recomputed, at an unclipped upper limit.**
  `src/model_BAwR.h:190-195` is explicit that this identity is what makes the
  endpoint shutdown clean and the live/frozen junction continuous. BTAwL
  currently recomputes it by formula (`:227`). Without the identity, the density
  will not integrate to the CDF near `t_max`.
- **`btawl_cdf` clamps `u` to `t_max` and evaluates there** (`:288-291`); the
  template instead sets `saturated = true` and lets the frozen branch carry it
  (`src/model_BAwR.h:168-176`). Equivalent in exact arithmetic, different code
  path — the `F(Inf)` defective-mass test from Phase 0 is what catches a mistake.

### 4e. Switch the `utils.h` adapters — own commit

`dbtawl_raw` (`utils.h:1707`), `pbtawl_raw` (`:1738`), `btawl_logS_at_t`
(`:1764`), `dbtawl_scalar` (`:1683`), `pbtawl_scalar` (`:1695`) onto
`btawl_log_pdf` / `btawl_log_cdf` / `log_btawl_surv_*`, matching `pbawf_raw`
(`utils.h:1978-2003`).

This is where fitted likelihoods move, and where the current damage is worst:
`btawl_logS_at_t:1788-1791` does `if (!(lc < 0.0)) { bad = true; break; }`, so
one CDF that the `[0,1]` clamp rounded to exactly 1 sends the **entire trial's**
survivor to `-Inf`. Fixing it will change truncated-likelihood values
substantially. Keep it separate so "the likelihood moved" can be bisected from
"the model got fast".

**Commit boundaries:** 4a; 4b; 4c-4d; 4e.

---

## 8. Phase 5 — BTAwL `(k, Ttrans)` chart

Adopting `(k, T_T)` over `(tau, T_max)`: at the sustained limit `pi = 1` the
`(tau, T_max)` pair leaves `k` encoded jointly through both coordinates, whereas
`(k, T_T)` keeps `k` directly identified and lets `T_T` become the irrelevant
one, mirroring the model structure. Two supporting facts, both verified:

- `r -> r*psi(r)` (where `psi(r) = T_max/tau`, a universal function of
  `r = k*tau` only — confirmed by scale invariance) is **monotone from 0 to Inf**,
  so `(k, T_T)` is an unconstrained bijection needing no joint bound.
  `(tau, T_max)` by contrast carries the hard constraint `T_max > tau`.
- `tau` back-solves from a **single** scalar equation
  `eta(T_T; tau) = k*K(T_T; k, tau)` — no nested `T_max` solve — monotone in
  `tau` with a **guaranteed bracket `(0, T_T)`**. Round-trips to machine
  precision across `k` in [0.3, 20] and `tau` in [0.05, 3].

Naming: **`Ttrans`, not `Tmax`**, since for `pi > 0` the mixture has no hard
endpoint and `T_T` is the transient channel's intrinsic window. `Ttransform`
reports `Tmax = Ttrans` for the pure transient model only.

- `chart = c("endpoint", "rate")` on both `BTAwL()` and `BTAwL_mixed()`,
  defaulting to `"endpoint"`, so `pi = 0` nesting stays exact in either chart.
  The rate chart exists because `k = 0` (finite total drive, no hard endpoint —
  a member the .tex discusses) has no endpoint to sample. **The endpoint chart
  therefore drops the `k = 0` exception and bounds `k` to `[1e-4, Inf)`**, which
  is what stops a chain wandering into the singularity mid-run; the rate chart
  keeps today's `exception = c(A = 0, k = 0)`. This mirrors BAwD, where
  `gamma = 1` keeps the rate chart for exactly this reason
  (`R/model_BAwD.R:470-476`) — except that `gamma` is a constructor constant
  whereas `k` is sampled, which is why the bound rather than a predicate does the
  work here.
- c_name suffix `_RATE` for the rate chart, endpoint unsuffixed. Parsed inside
  the existing `BTAwL_MIX` / `BTAwL` branches (`src/particle_ll.cpp:788, 803`),
  so the substring-ordering hazard (`BTAwL_MIX` must precede `BTAwL`) is
  unaffected.
- `src/col_registry.h:185-212` — `spec()` names slot 6 `Ttrans` and `spec_rate()`
  names it `tau`; `btawl_mix` likewise for slot 7 (`tau_t`). Enum member renamed
  `clear`. A single predicate `btawl_uses_ttrans(ctx)` drives **both** the spec
  choice and every column read, as `col_registry.h:107-114` documents for BAwD.
- `src/model_BTAwL.h` — `btawl_tau_from_ttrans(k, Ttrans)` (bracketed Newton),
  and `src/utils.h` gets one `btawl_clear_to_tau(ctx, clear, k)` funnel mirroring
  `bawd_clear_to_ell`.

### 8.1 One solver, and store the result

Unlike `bawd_ell_from_tmax`, this back-solve is a root find, so two rules:

1. **Store the solved `tau` and `t_max = Ttrans` on the geometry**, so
   `btawl_geometry` never runs a forward `T_max` solve at all. Landing Phase 5
   without this re-introduces a per-evaluation solve and undoes Phase 4a.
2. **Export exactly one solver.** `Ttransform` needs `tau` for reporting and for
   the simulator; call the exported `btawl_tau_vec()` from R, exactly as
   `R/model_BAwD.R:522` calls `bawd_ell_vec`. Do **not** write a second solver in
   R — two solvers with different tolerances is the classic "R and C++ disagree
   at 1e-9" test failure.

`Ttransform` reports the derived `tau` so `rBTAwL` (which reads
`pars[j, "tau"]`, `R/model_BTAwL.R:209`) is untouched, plus `Tmax`, `rt_max` and
`Vcrit`. Add the missing `btawl_vcrit_vec` export for parity with
`bawf_vcrit_vec` / `bawr_vcrit_vec`.

### 8.2 The R simulator will drift

`R/model_BTAwL.R:134-165` carries its own `.btawl_tmax` (doubling + `uniroot`)
and `.btawl_hit_time`. There is no C++ BTAwL sampler, so R is the only sampler
and the `make_data` -> refit tests depend on it agreeing with the kernel. Both
Phase 4a and Phase 5 invalidate it. Point `.btawl_tmax` at the exported C++
solver.

### 8.3 The mixture nesting invariant

`test-btawl.R:24` asserts `pi = 0` is exactly transient BTAwL. If BTAwL moves to
`(k, Ttrans)` while the mixture stays on `tau_t`, that test either breaks or must
compare through the converted `tau`. Move both charts together.

Document that `Ttrans` is unidentified as `pi -> 1` and should be held constant
there.

**Commit.**

---

## 9. Phase 6 — BTAwL_mixed closed-form live branch

Currently adaptive `qags` over `z`, with a 96-point logarithmic grid scan plus
bisections inside each node, and the pdf integrand repeating the entire scan for
its running-minimum test: 111 us/call, ~280x the family. Almost all of it is
avoidable, because the mixture is as closed-form as the pure model:

- `V*_mix(u,z) = (b - z*exp(-k*u))/H_mix(u)` is **affine in z**, so
  `btawl_live_cdf` / `btawl_live_pdf` apply verbatim with `H -> H_mix`.
- `btawl_mix_d(u,z)` is also affine in `z`, so the live/frozen cut is the same
  closed-form `zcut = b*H'_mix(u) / (exp(-k*u)*input(u))`.
- Even the temporary-freeze case stays an interval: the running minimum is a min
  of functions affine in `z`, hence concave, so `V* - m` is convex and `{<= 0}`
  is an interval.

So the mixture needs at most two 1-D root solves per `u`, not a 2-D numerical
scheme. Largest payoff, most design work; land last, behind a test pinning the
current values to MC/quadrature tolerance first. Note that once BTAwL is fast the
mixture becomes the only cost anyone notices.

---

## 10. Verification

1. **Always install into a temp lib** and check the harness path banner —
   `library(EMC2)` loads a stale install.
2. **Exact-survivor harness.** `S = (1/A) * integral_0^A G(v_req(z)) dz`
   integrates the launch CDF, so it has no cancellation and is the reference for
   all six models. Assert 1e-8 relative agreement down to `log S = -600`.
3. **Likelihood monotonicity.** The `v_lMTRUE` sweep must decrease monotonically
   well past the current plateau at `60 x min_ll`.
4. **Timing thresholds** from the measured baselines: BTAwL cdf 9.8 -> <1 us,
   BTAwL pdf 2.45 -> <0.6 us, `btawl_tmax` 2.48 -> <0.15 us, mixture 111 -> <2 us.
   (BTAwL will not reach BAwR's 0.25 us: `btawl_h` needs two `exp()` calls plus a
   series branch where BAwR has a power law.) BAwF 0.40 / BAwR 0.25 / BAwD 0.62 /
   BAwL 0.28 / LBA 0.18 must not regress.
5. **Bit-comparison gates.** Phases 1, 2 and 4a must not change likelihood values
   beyond 1e-12 on a saved fit. Phase 4b is ~1e-7 against the retained qags
   reference. Phases 3, 4c-4e change tails by design and use the exact reference
   instead.
6. **Full `testthat` run** against the ~5,286-test baseline at each commit.
7. **Profile `p` under both BAwR charts.** `p -> Inf` is a known BAwR pathology;
   under the rate chart it collapses `T_max` to 1, under the endpoint chart it
   sends `kappa` to 0 or Inf. Confirm the chart change did not make it worse.

---

## 11. Explicitly deferred

### 11.1 Shared-header extraction

After name normalisation, `src/model_BAwF.h:342-800` and
`src/model_BAwR.h:277-734` differ **only** in four comment blocks and two
function signatures — the four `log_*` and four `*_natural_*` bodies are
identical. BAwD's CDF half is identical too; its PDF half is the affine-in-`w`
generalisation `wgt = E_g*(E_rel*w - ell)` of BAwF/BAwR's `wgt = w - c`.

Notably, **BTAwL's PDF is already in that generalised form**
(`src/model_BTAwL.h:236-247`, `out = (c0*mass + c1*first)/(A*q)`), and its `q` is
`h/E` rather than the other three models' `u`. That is the argument for
conforming BTAwL **first**: the interface cannot be designed correctly against
three clients that agree on both points, and extraction done afterwards is a
provable bit-identical no-op across four clients, which extraction done first is
not.

When it happens, the shape should be `template <class Policy>` free functions
over a struct of `static` members (matching the existing `bawd_log_gl` /
`bawd_log_gl_split` idiom at `src/model_BAwD.h:423, 447`) — never virtual
dispatch. **Do not** parameterise over the launch distribution: the normal path
uses `log_normal_phi_integral` while the lognormal path uses the closed
`log_lognormal_*_stoploss` primitives, and hiding that behind a common survivor
policy would force the closed forms back into quadrature. Keep the `launch` int
dispatch in the thin shims where it already lives.

Two interface fixes to make at that point: take `const Geom&` rather than raw
parameters (today `bawf_cdf_norm` builds the geometry twice per evaluation on the
fallback path — tolerable at 0.14 us, not at BTAwL's cost), and stage the merge
CDF-first, PDF-second, because BAwD's `gamma_zero`/`ell_zero` log-accuracy
branches (`src/model_BAwD.h:726-734`) are the most likely place to silently
perturb BAwD.

This is the down-payment on `cpp_architecture_refactor_plan.md:163`. Note the
include graph is already circular (`src/model_BAwD.h:26` includes
`model_LBA.h`; `src/model_LBA.h:1096` includes `model_BAwD.h`, working only
because LBA is always entered first), so the new header must be a **leaf** and
should absorb `BA_ACCEPT_*`, `BAWL_NATURAL_*`, the `*_DENOM_FLOOR` constants,
`BAWL_LAUNCH_*`, `log_positive_normalizer`, `natural_normalizer`,
`natural_cdf_safe`, `natural_normal_interval_safe` and the `bawd_log_gl*`
templates — which breaks the cycle as a side effect and resolves the complaint
already recorded at `src/model_LBA.h:45-48`.

### 11.2 Other

- The `logS1_ptr` scalar-path slot (§6, Phase 3b).
- Log-survivor accuracy for the 13 non-BA race kernels.
- The `gsl_set_error_handler_off()` save/restore pattern mutates a process-global
  handler. It appears at 13 sites across 8 files (`model_LBA.h`, `model_RDM.h`,
  `model_SS_EXG.h`, `ss_raw.h`, `tools.cpp`, `particle_ll.cpp`, `model_BTAwL.h`),
  so it is the house pattern, and it is safe under the current fork-based worker
  pool. It would become a race under the OpenMP threading contemplated in
  `cpp_architecture_refactor_plan.md`; flag it there, not here.

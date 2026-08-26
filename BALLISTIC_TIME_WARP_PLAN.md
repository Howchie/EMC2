# Implementation plan: the single-parameter operational-time warp for the ballistic family

Source specification: `Math/ballistic-time.md`.
Target branch: `playground` (worktree `/data/work/EMC2_dev_oo`).
Status: **not started**. Everything below is design; no code has been written.

This document is written so that an implementer who has not read the codebase
can execute it end to end. Read sections 1–4 before touching anything: they
contain the decisions that make the rest a mechanical exercise.

---

## 1. What is being added, in one paragraph

Every ballistic accumulator in EMC2 is written as a deterministic trajectory in
*accumulation time* `u = rt - t0`. The warp inserts one strictly increasing
reparameterisation of that clock,

```
s = c_eta(u) = ((1 + u)^omega - 1) / omega ,   omega = exp(eta)
```

between physical accumulation time `u` and the *operational* time `s` that the
existing kernels already consume. Because every parent model's CDF, density,
survivor, endpoint, and defect are already functions of `s`, the entire
extension is

```
F_eta(u) = F_0(c_eta(u))
log f_eta(u) = log f_0(c_eta(u)) + (omega - 1) * log1p(u)
u = c_eta^{-1}(s) = (1 + omega s)^{1/omega} - 1        (simulation)
```

and nothing inside any kernel changes. `eta = 0` gives `omega = 1`,
`c_0(u) = u`, and zero log-Jacobian, i.e. the exact parent model.

**One new sampled parameter, `eta`, per ballistic model. No new kernels, no new
quadrature, no change to any closed form.**

---

## 2. Scope

### 2.1 In scope (models that gain `eta`)

| Model | Constructor(s) | c_name prefix | Dispatch branch (`src/race_dispatch.cpp`) |
|---|---|---|---|
| LBA | `LBA()` | `LBA`, `LBAIO` | line 514 |
| BAwL (clock-free only, see 2.2) | `BAwL(erlang_type = "none", correlated = FALSE)` | `BAwL*` | line 471 |
| BAwD | `BAwD()` | `BAwD*` | line 423 |
| BAwDp | `BAwDp()` | `BAwDp*` | line 398 |
| BAwF | `BAwF()` | `BAwF*` | line 327 |
| BAwR | `BAwR()` | `BAwR*` | line 365 |
| BTAwL (full local race) | `BTAwL()` | `BTAwL*` | line 296 |
| BTAwL transient | `BTAwLTransient()` | `BTAwL_TRANSIENT*` | line 265 |
| BTAwL sustained | `BTAwLSustained()` | `BTAwL_SUSTAINED*` | line 241 |

Line numbers are from the state of `src/race_dispatch.cpp` at the time of
writing; re-derive them with
`grep -n '} else if (type_std.find(' src/race_dispatch.cpp` before editing.
Dispatch is by substring and the branch order is load-bearing (`BAwDp` must
precede `BAwD`, `BTAwL_SUSTAINED`/`BTAwL_TRANSIENT` must precede `BTAwL`), so
add lines inside existing branches — never reorder them.

All launch distributions (`normal`, `lognormal`, `splitlognormal`, `weibull`),
all `IO` (`posdrift = FALSE`) variants, all fixed-kernel suffixes
(`_RHO1/2/4`, `_GAM12/23/34/100`), and both BTAwL charts (`rate`, `endpoint`)
are covered automatically, because the warp is installed on the resolved
adapter's function pointers rather than inside any kernel.

### 2.2 Out of scope for this change (and why)

These get **no `eta` in `p_types`**, so `design(formula = list(eta ~ 1))`
fails naturally with "eta ... not found in model p_types".

1. **BAwL with an Erlang kill/guess clock** (`erlang_type != "none"`).
   `dkilledleakyba_norm` / `pkilledleakyba_norm` (`src/model_LBA.h:1026,1108`)
   run the Erlang clocks on **raw** `t` (from stimulus onset), not on
   accumulation time. `Math/ballistic-time.md` says nothing about them, and
   "does the operational clock also govern an external Poisson timer?" is a
   modelling decision, not an implementation detail. Composing the warp there
   would either warp the clocks (changing `mG`/`mK`'s units silently) or leave
   them unwarped (requiring a rewrite of the killed kernels). Refuse instead.
   *Extension path:* decide the semantics, then either warp the raw-time API
   per §15 of the spec (`T* = t0 + c_eta(T - t0)`, clocks then run on
   operational time measured from `t0`) or add an unwarped clock argument.

2. **BAwLcorr** (`BAwL(correlated = TRUE)`). The correlated-drift path has its
   own inlined kernels and fused fast routes
   (`src/correlated_likelihood.cpp:300-360`, `src/bawl_corr_exact.h`,
   `src/bawl_corr_counters.h`) that bypass `adapter.pdf1_ptr` /
   `model_dfun_raw`, so installing the wrapper would silently miss them.
   *Extension path:* warp `tau` inside `bawl_corr_row_geometry`,
   `bawl_corr_single_log_density` (add the Jacobian) and
   `bawl_corr_single_log_survival`, plus the fused raw route at
   `src/correlated_likelihood.cpp:1277`.

3. **`LBA_LogicalRules`** (`R/model_LBA.R:265`). Builds its own `p_types` and
   runs through `src/logicalrules_likelihood.cpp`, which does its own `t0`
   bookkeeping and never calls the race adapter's raw kernels.

4. Non-ballistic race models (RDM, RDMSWTN, ROU, ROUp, RLF, GOM, LNR, REXG,
   PCOUNTER, FRQ) and the two-boundary models (DDM, BOU). They are not in the
   spec.

### 2.3 Explicitly preserved

* Every closed form. The warp adds no integral (spec §13).
* The intrinsic defect: `F_eta(Inf) = F_0(Inf)` for every `eta` (spec §3).
* Critical launches `V_crit`, saturation launches, and the *operational*
  endpoint `T_max`. Only the *physical* endpoint `rt_max` moves.
* Bit-identical parent behaviour at `eta = 0`. This is a hard requirement, not
  an aspiration — see 4.1.

---

## 3. Six decisions, with rationale

Do not re-litigate these while implementing; if one turns out to be wrong,
change it deliberately and update this file.

### D1. `eta` is a trailing, always-present `p_type`, not a `c_name` variant

`R/utils.R:70-73` already documents the mechanism:

> Both are *trailing* p_types: `design()` turns any p_type absent from the
> formula into a constant at its default, and `emc2col::validate_col_prefix`
> only checks the canonical prefix of the column order, so appending here is
> free.

So appending `eta = 0` to each ballistic model's `p_types`:

* changes nothing for any existing design or saved fit (it becomes a constant
  at 0, and 0 is the exact parent);
* needs **zero** changes to `src/col_registry.h` (no new `ColSpec`, no new
  enum, no new namespace — a `_TW` suffix would have needed ~30 of them, one
  per model × launch distribution);
* keeps one `c_name` per model, so `compare()` labels and saved fits are
  unaffected;
* is turned on by the user simply writing `eta ~ 1` in the formula.

### D2. `eta` is resolved **by name**, not positionally

`src/col_registry.h`'s own comment on `frq::delta` explains why an *optional*
trailing column resolved *positionally* is dangerous. We sidestep that entirely
by copying the mechanism already used for `rho`
(`configure_corr_drift_context`, `src/race_dispatch.cpp:576`): scan
`keep_names` (which is exactly `names(model$p_types)`; see
`src/particle_ll.cpp:418` and `R/sampling.R:1635`) for `"eta"` and store the
index in the context. Absent ⇒ index `-1` ⇒ the warp is never installed.

### D3. The warp is a **generic outer layer**, installed on the adapter's
function pointers — not edited into each model's kernels

This is spec §20 ("Do **not** duplicate the existing ... mathematical kernels")
taken literally, and it is what guarantees the user's requirement that `eta`
means *exactly* the same thing in every model. There are 7 model families ×
5 adapter entry points ≈ 35 places where a hand-written warp could drift
(BAwL's `dbawl_raw` alone has three branches). One wrapper cannot drift.

Concretely, `configure_time_warp_context()` moves the resolved adapter's five
function pointers into `ctx.tw.base_*` and replaces them with five generic
wrappers. Every kernel keeps reading `tt = rt - t0` and receives `s` there
instead of `u`.

### D4. `eta` is **not** in `p_types_canonical`

`R/design.R:358-370` prints
`"Parameter(s) X not specified in formula and assumed constant."` for every
*canonical* p_type missing from the formula. Making `eta` canonical would emit
that message for every pre-existing ballistic design. Treat it exactly like
`pContaminant`/`pGuess`: a real parameter that is off by default and silent
until asked for.

### D5. `eta` is on the **accumulation clock**, measured from `t0`, in **seconds**

`u = rt - t0` for every model, so `eta` has one meaning family-wide. The
spec hard-codes the 1-second reference (`1 + u`), which means **`eta` is not
scale-free**: fitting the same data in milliseconds would give a different
`eta`. EMC2 is a seconds package throughout (`t0` lower bound `0.05`,
`rt_resolution` default `1/60`), so this is consistent, but it must be said
out loud in the docs.

`c'(0) = 1`: operational and physical clocks agree at accumulation onset, so
`eta` is a *curvature*, not a rate rescaling. `eta < 0` slows operational time
progressively (convex `c^{-1}`, more right skew); `eta > 0` speeds it up
(concave `c^{-1}`, less right skew).

### D6. Warp per accumulator row is *permitted* but *discouraged*

The kernel reads `eta` from the accumulator's own parameter row, so
`eta ~ lM` is expressible and the likelihood stays mathematically valid
(independent racers, each with its own warped marginal). It is **not** the
mechanism the spec describes (§12: separate warps change winner ordering), so
the documentation must recommend `eta ~ 1`. Do not add a runtime check — EMC2
does not police design formulas elsewhere.

---

## 4. The canonical numerics

### 4.1 Identity at `eta == 0` must be *exact*, not approximate

`expm1(log1p(u))` is **not** bit-identical to `u`. Every helper must
short-circuit:

```cpp
if (eta == 0.0) return u;        // fwd
if (eta == 0.0) return 0.0;      // log_jac
if (eta == 0.0) return s;        // inv
```

and the wrappers must additionally tail-call the base function pointer when no
row in the batch has a non-zero `eta`. That way an `eta`-carrying model with
`eta` constant at 0 produces byte-identical log-likelihoods and byte-identical
simulator draws (the RNG stream is untouched because no extra variate is drawn).
This is what makes D1 safe.

Note `eta == 0.0` is false for `NaN`, so a `NaN` eta correctly propagates `NaN`
into the density, where the existing `emc2_isfinite` guards floor it.

### 4.2 Reference implementation (C++)

```cpp
namespace emc2tw {

// s = c_eta(u) = ((1+u)^omega - 1)/omega,  omega = exp(eta)
inline double fwd(double u, double eta) {
  if (eta == 0.0) return u;                 // exact parent
  if (!(u > 0.0)) return u;                 // u <= 0 and NaN pass through
  if (!R_FINITE(u)) return u;               // c(+Inf) = +Inf for every omega > 0
  const double l1p   = std::log1p(u);
  const double omega = std::exp(eta);
  if (omega < 1e-12) return l1p;            // omega -> 0 limit: c(u) -> log1p(u)
  const double x = omega * l1p;
  if (x > 709.0) return R_PosInf;           // overflow guard
  return std::expm1(x) / omega;
}

// log c'_eta(u) = (omega - 1) log(1 + u)
inline double log_jac(double u, double eta) {
  if (eta == 0.0) return 0.0;
  if (!(u > 0.0) || !R_FINITE(u)) return 0.0;
  return (std::exp(eta) - 1.0) * std::log1p(u);
}

// u = c_eta^{-1}(s) = (1 + omega s)^{1/omega} - 1
inline double inv(double s, double eta) {
  if (eta == 0.0) return s;
  if (!(s > 0.0)) return s;
  if (!R_FINITE(s)) return s;
  const double omega = std::exp(eta);
  if (omega < 1e-12) return std::expm1(s);  // omega -> 0 limit: c^{-1}(s) -> e^s - 1
  const double y = std::log1p(omega * s) / omega;
  if (y > 709.0) return R_PosInf;
  return std::expm1(y);
}

}  // namespace emc2tw
```

Never write `pow(1 + u, omega)` or `pow(1 + omega*s, 1/omega)` anywhere. The
`log1p`/`expm1` forms above are the only permitted spellings (spec §14).

### 4.3 Bounds on `eta`

`minmax = c(-5, 5)`, i.e. `omega` in `(0.0067, 148)`. Rationale: `omega = 148`
already makes `(1+u)^(omega-1)` ~1e150 at `u = 1 s`, and `eta` beyond ±5 is not
identifiable from RT data — the bound exists to keep the overflow guards in 4.2
unreachable in practice, not to express a belief. No `exception` entry; `0` is
interior. Default prior is EMC2's generic `N(0, 1)` on the identity scale,
which puts `omega` in roughly `(0.37, 2.7)` at ±1 sd — a sensible default.

---

## 5. Work plan

Do the stages in order. Stages 1–2 are independently testable; do not start
stage 3 before stage 1 compiles.

### Stage 0 — Baseline (do this first, it takes 10 minutes and saves hours)

```bash
cd /data/work/EMC2_dev_oo
export R_LIBS_TW=$(mktemp -d)
rm -f src/*.o src/EMC2.so
R CMD INSTALL --preclean --no-multiarch --library="$R_LIBS_TW" . 2>&1 | tail -20
R_LIBS_USER="$R_LIBS_TW" Rscript -e 'library(EMC2); cat(find.package("EMC2"), "\n")'
R_LIBS_USER="$R_LIBS_TW" Rscript -e 'testthat::test_local(reporter="summary")' \
  2>&1 | tail -60 > /tmp/tw_baseline.txt
```

**Record which tests already fail on `playground` before your change.** This
branch is known to carry pre-existing failures; you must not be blamed for
them, and you must not "fix" them. Diff against `/tmp/tw_baseline.txt` at the
end.

> Trap (previously burned this repo): plain `library(EMC2)` can load a **stale
> installed copy** rather than your build. Always print `find.package("EMC2")`
> and confirm it is under `$R_LIBS_TW`.

> Trap: `src/Makevars` has **no header dependency tracking**. Adding a field to
> `ContextForRaceModels` in `src/race_contract.h` changes the struct layout for
> ~20 translation units. A partial rebuild produces an ODR/ABI mismatch that
> segfaults at run time with no compile error. **`rm -f src/*.o src/EMC2.so`
> before every rebuild in stages 1–2.**

---

### Stage 1 — C++ likelihood core

#### 1a. `src/race_contract.h` — add the plan struct and one context member

Insert immediately after the `RaceLogSAtTFun` typedef (currently ~line 80,
after the recent `EndpointQueryPlan` addition):

```cpp
// Operational-time warp (Math/ballistic-time.md).  The warp is installed as a
// generic outer layer over a ballistic model's five adapter entry points:
// configure_time_warp_context() moves the model's own pointers into base_*
// here and replaces the adapter's with the tw_* wrappers in time_warp.h.  The
// kernels themselves are untouched -- they keep reading tt = rt - t0 and are
// simply handed s = c_eta(tt) instead of tt.
struct TimeWarpPlan {
  int eta_index = -1;                       // keep_names position of `eta`; -1 = absent
  bool supported = false;                   // set true only by ballistic dispatch branches
  RacePdf1Fun     base_pdf1  = nullptr;
  RaceCdf1Fun     base_cdf1  = nullptr;
  RaceRawFun      base_d_raw = nullptr;
  RaceRawFun      base_p_raw = nullptr;
  RaceLogSAtTFun  base_logS  = nullptr;
  bool installed() const { return eta_index >= 0 && base_pdf1 != nullptr; }
};
```

Add to `ContextForRaceModels`, next to `t0_index`:

```cpp
    TimeWarpPlan tw;
```

#### 1b. New file `src/time_warp.h`

```cpp
#ifndef EMC2_TIME_WARP_H
#define EMC2_TIME_WARP_H

#include <Rcpp.h>
#include <cmath>
#include <string>

#include "race_contract.h"
#include "race_dispatch.h"

namespace emc2tw { /* fwd, log_jac, inv exactly as in section 4.2 */ }

// Generic warp wrappers.  Signatures match RacePdf1Fun / RaceCdf1Fun /
// RaceRawFun / RaceLogSAtTFun so they can be dropped into RaceModelAdapter.
double tw_pdf1(double t, const double* par, void* ctx_);
double tw_cdf1(double t, const double* par, void* ctx_);
void   tw_d_raw(const double* rt, const double* const* cols, int n_rows,
                const int* mask, const int* isok, double* out,
                double min_ll, void* ctx_);
void   tw_p_raw(const double* rt, const double* const* cols, int n_rows,
                const int* mask, const int* isok, double* out,
                double min_ll, void* ctx_);
void   tw_logS_at_t(double t, const double* const* cols, int n_rows_total,
                    int n_lR, int n_par, const int* trunc_mask,
                    int n_unique_trials, const int* isok_all, void* ctx_,
                    double* logS_out);

// Resolve the `eta` column by name and install the wrappers.  No-op when the
// design has no `eta`; hard error when it has one but the model does not
// support the warp.  Mirrors configure_corr_drift_context().
void configure_time_warp_context(RaceModelAdapter& adapter,
                                 const Rcpp::CharacterVector& keep_names,
                                 const std::string& caller);

#endif  // EMC2_TIME_WARP_H
```

#### 1c. New file `src/time_warp.cpp`

`configure_time_warp_context`:

```cpp
void configure_time_warp_context(RaceModelAdapter& adapter,
                                 const Rcpp::CharacterVector& keep_names,
                                 const std::string& caller) {
  int idx = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == "eta") { idx = j; break; }
  }
  if (idx < 0) return;                       // no warp column: leave the adapter alone
  if (!adapter.ctx.tw.supported) {
    Rcpp::stop("%s: the design supplies an 'eta' column but this model does not "
               "implement the operational-time warp (Math/ballistic-time.md).",
               caller.c_str());
  }
  if (adapter.ctx.t0_index < 0) {
    Rcpp::stop("%s: the operational-time warp requires a purely additive t0 "
               "column.", caller.c_str());
  }
  adapter.ctx.tw.eta_index  = idx;
  adapter.ctx.tw.base_pdf1  = adapter.pdf1_ptr;
  adapter.ctx.tw.base_cdf1  = adapter.cdf1_ptr;
  adapter.ctx.tw.base_d_raw = adapter.model_dfun_raw;
  adapter.ctx.tw.base_p_raw = adapter.model_pfun_raw;
  adapter.ctx.tw.base_logS  = adapter.logS_at_t_ptr;

  if (adapter.pdf1_ptr       != nullptr) adapter.pdf1_ptr       = &tw_pdf1;
  if (adapter.cdf1_ptr       != nullptr) adapter.cdf1_ptr       = &tw_cdf1;
  if (adapter.model_dfun_raw != nullptr) adapter.model_dfun_raw = &tw_d_raw;
  if (adapter.model_pfun_raw != nullptr) adapter.model_pfun_raw = &tw_p_raw;
  if (adapter.logS_at_t_ptr  != nullptr) adapter.logS_at_t_ptr  = &tw_logS_at_t;
}
```

Scalar wrappers. These are the GSL censoring/truncation integrand path, so they
must stay allocation-free:

```cpp
double tw_pdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double t0  = par[ctx->t0_index];
  const double eta = par[tw.eta_index];
  const double u   = t - t0;
  if (eta == 0.0 || !(u > 0.0) || !R_FINITE(u))
    return tw.base_pdf1(t, par, ctx_);        // includes t == +Inf: base returns 0
  const double p = tw.base_pdf1(t0 + emc2tw::fwd(u, eta), par, ctx_);
  if (!(p > 0.0) || !R_FINITE(p)) return 0.0;
  return std::exp(std::log(p) + emc2tw::log_jac(u, eta));
}

double tw_cdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double t0  = par[ctx->t0_index];
  const double eta = par[tw.eta_index];
  const double u   = t - t0;
  if (eta == 0.0 || !(u > 0.0) || !R_FINITE(u))
    return tw.base_cdf1(t, par, ctx_);        // t == +Inf must reach the base:
                                              // cdf1(Inf) is F_max, not 1, and
                                              // defect invariance depends on it
  return tw.base_cdf1(t0 + emc2tw::fwd(u, eta), par, ctx_);
}
```

Batch density. The **only** subtle part of the whole change is the interaction
with `min_ll` flooring: the Jacobian must be added *before* the floor, or a
floored zero density picks up the Jacobian and lands slightly off the floor.
`ContextForRaceModels::floor_raw_log_lik` is already a per-call switch that
`particle_ll.cpp` toggles (search `dense_floor_raw_log_lik_prev`), so borrow it:

```cpp
void tw_d_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double* t0_  = cols[ctx->t0_index];
  const double* eta_ = cols[tw.eta_index];

  bool any = false;
  if (eta_ != nullptr && t0_ != nullptr) {
    for (int i = 0; i < n_rows && !any; ++i) if (mask[i] && eta_[i] != 0.0) any = true;
  }
  if (!any) {                                  // bit-identical parent path
    tw.base_d_raw(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
    return;
  }

  static thread_local std::vector<double> rt_buf, lj_buf;
  rt_buf.assign(rt, rt + n_rows);
  lj_buf.assign(static_cast<size_t>(n_rows), 0.0);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double eta = eta_[i];
    if (eta == 0.0) continue;
    const double u = rt[i] - t0_[i];
    if (!(u > 0.0) || !R_FINITE(u)) continue;  // guards fire identically on u and s
    rt_buf[i] = t0_[i] + emc2tw::fwd(u, eta);
    lj_buf[i] = emc2tw::log_jac(u, eta);
  }

  const bool floor_raw = ctx->floor_raw_log_lik;
  ctx->floor_raw_log_lik = false;              // unfloored, so the Jacobian is exact
  tw.base_d_raw(rt_buf.data(), cols, n_rows, mask, isok, out, min_ll, ctx_);
  ctx->floor_raw_log_lik = floor_raw;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    out[i] = raw_log_value(out[i] + lj_buf[i], min_ll, floor_raw);
  }
}
```

`raw_log_value` (declared in `race_contract.h`) maps a non-finite input to
`raw_log_zero(min_ll, floor_raw)`, so `-Inf + lj = -Inf` correctly becomes
`min_ll` again when flooring is on. **Do not** toggle `floor_raw_log_lik` in
`tw_p_raw` — it writes log-survivors, has no Jacobian, and uses `0.0` (not
`min_ll`) for its skip branches:

```cpp
void tw_p_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* ctx_) {
  /* identical prologue building rt_buf; then simply: */
  tw.base_p_raw(rt_buf.data(), cols, n_rows, mask, isok, out, min_ll, ctx_);
}
```

Batch survivor at a scalar `t` (the truncation/censoring normaliser). Here `t`
is shared across rows but `t0` and `eta` are per row, so the warp cannot be
folded into `t`. Instead **substitute the `t0` column**, chosen so that the base
kernel's own `tt = t - t0` computes `s`:

```cpp
void tw_logS_at_t(double t, const double* const* cols, int n_rows_total,
                  int n_lR, int n_par, const int* trunc_mask,
                  int n_unique_trials, const int* isok_all, void* ctx_,
                  double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double* t0_  = cols[ctx->t0_index];
  const double* eta_ = cols[tw.eta_index];

  bool any = false;
  if (eta_ != nullptr && t0_ != nullptr && R_FINITE(t)) {
    for (int r = 0; r < n_rows_total && !any; ++r) if (eta_[r] != 0.0) any = true;
  }
  if (!any) {
    tw.base_logS(t, cols, n_rows_total, n_lR, n_par, trunc_mask,
                 n_unique_trials, isok_all, ctx_, logS_out);
    return;
  }

  static thread_local std::vector<double> t0_buf;
  static thread_local std::vector<const double*> cols_buf;
  t0_buf.assign(t0_, t0_ + n_rows_total);
  for (int r = 0; r < n_rows_total; ++r) {
    const double eta = eta_[r];
    if (eta == 0.0) continue;
    const double u = t - t0_[r];
    if (!(u > 0.0)) continue;                  // leave t0 alone so the base's
                                               // "not started" branch is unchanged
    t0_buf[r] = t - emc2tw::fwd(u, eta);       // => t - t0_buf[r] == c_eta(u)
  }

  // Mirror particle_ll.cpp's padding contract: kernels may FETCH (never
  // dereference) pointers for optional trailing columns this variant lacks.
  const int n_slots = std::max(n_par, 16);
  cols_buf.assign(static_cast<size_t>(n_slots), nullptr);
  for (int j = 0; j < n_par; ++j) cols_buf[j] = cols[j];
  cols_buf[ctx->t0_index] = t0_buf.data();

  tw.base_logS(t, cols_buf.data(), n_rows_total, n_lR, n_par, trunc_mask,
               n_unique_trials, isok_all, ctx_, logS_out);
}
```

Only copy `cols[j]` for `j < n_par`; reading past `n_par` in the source is
undefined (`race_endpoint_prepare_groups` sizes `compact_col_ptrs` at exactly
`n_par`, `src/race_integrands.cpp:99`).

#### 1d. `src/race_dispatch.cpp` — mark the nine supported branches

Add exactly one line inside each of the branches listed in 2.1:

```cpp
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
```

For the `BAwL` branch (line ~471) make it conditional, matching 2.2:

```cpp
    // The Erlang kill/guess clocks run on RAW time and the correlated path has
    // its own inlined kernels, so neither composes with the operational-time
    // warp yet.  The R constructor already withholds `eta` from those variants;
    // this is the defensive half of the same contract.
    out.ctx.tw.supported = !out.ctx.kill_active && !out.ctx.corr_drift_active;
```

Place it *after* `kill_active` and `corr_drift_active` are set. The `LBA` branch
(line ~514) is unconditional (`bawl_clocks_fixed_zero`/`_off` are always true
there).

#### 1e. `src/race_dispatch.h` — nothing

`configure_time_warp_context` is declared in `time_warp.h`, which includes
`race_dispatch.h`. Keep it that way so `race_dispatch.cpp` does not need to know
about the wrappers.

#### 1f. `src/particle_ll.cpp` — two call sites

After each existing `configure_rdmswtn_corr_context(...)` call, add:

```cpp
    configure_time_warp_context(adapter, keep_names, "calc_ll_oo");
```

at `src/particle_ll.cpp:788` (in `calc_ll_oo`) and

```cpp
    configure_time_warp_context(adapter, keep_names, "calc_ll_oo_pw");
```

at `src/particle_ll.cpp:1264` (in `calc_ll_oo_pw`). Add
`#include "time_warp.h"` to the include block.

**Order matters**: it must run *after* the dispatch branch has set `t0_index`
and `tw.supported`, and after `validate_col_prefix`.

#### 1g. `src/col_registry.h` — a comment only

Add, near the top with the other contract notes:

```cpp
// `eta` (the operational-time warp of Math/ballistic-time.md) is deliberately
// NOT in any enum here: it is resolved BY NAME in configure_time_warp_context()
// (src/time_warp.cpp), like `rho`, so that adding it to a model's p_types needs
// no new ColSpec and cannot shift any positional index.
```

#### 1h. Build and prove stage 1

```bash
rm -f src/*.o src/EMC2.so
R CMD INSTALL --preclean --no-multiarch --library="$R_LIBS_TW" .
R_LIBS_USER="$R_LIBS_TW" Rscript -e 'testthat::test_local(reporter="summary")' | tail -40
```

At this point **no R model exposes `eta`**, so the wrappers are never installed
and the full suite must match `/tmp/tw_baseline.txt` exactly. If it does not,
you have an ABI/rebuild problem (see the Stage 0 trap), not a logic problem.

---

### Stage 2 — C++ simulators (`src/model_rng.cpp`)

Spec §16: simulate the parent's operational finishing time `S`, then return
`U = c_eta^{-1}(S)`. Do this **per accumulator, before the race is resolved and
before `t0` is added**, so that the winner is decided on physical times.

Add `#include "time_warp.h"` and, inside each simulator, right after
`const auto ci = col_index_map(pars);`:

```cpp
  const int ieta = ci.count("eta") ? ci.at("eta") : -1;    // -1: no warp column
```

then a local helper:

```cpp
  auto warp_u = [&](double u, int r) -> double {
    return (ieta < 0) ? u : emc2tw::inv(u, pars(r, ieta));
  };
```

The seven simulators fall into **two conventions** — get this right or `t0`
gets warped too:

| Function | line | current | change to |
|---|---|---|---|
| `rlba_cpp` | 359 | `dt[r] = (d < 0.0) ? R_PosInf : d;` (t0 added by `resolve_race`) | `dt[r] = (d < 0.0) ? R_PosInf : warp_u(d, r);` |
| `rbawl_cpp_impl` | ~610 | `dt[r] = d + pars(r, it0);` (t0 added here) | `dt[r] = warp_u(d, r) + pars(r, it0);` |
| `rbta_wl_cpp` | 704 | `dt[r] = R_FINITE(hit) ? hit + t0 : R_PosInf;` (×3 branches: mode 0, 1, full) | `dt[r] = R_FINITE(hit) ? warp_u(hit, r) + t0 : R_PosInf;` |
| `rbawd_cpp` | 868 | `dt[r] = (R_FINITE(u) && u >= 0.0) ? u : R_PosInf;` | `... ? warp_u(u, r) : R_PosInf;` |
| `rbawf_cpp` | 933 | same | same |
| `rbawr_cpp` | 994 | same | same |
| `rbawdp_cpp` | 1054 | same | same |

`rbta_wl_cpp` has **three** assignment sites (transient-only, sustained-only,
and the `fmin(t_T, t_S)` full race). Warp all three. Because the warp is
monotone and shared by both channels, warping the `fmin` is identical to
warping each channel and then taking the `fmin` — the spec's §11 requirement
that the two channels share one clock is satisfied either way, but warping
after the `fmin` is one line instead of two.

Leave `rbawl_corr_cpp` (line 840) untouched — out of scope (2.2).

`emc2tw::inv(+Inf, eta) == +Inf`, so intrinsic omissions stay omissions; that is
spec §16 step 2 and it needs no extra branch.

---

### Stage 3 — R parameter machinery

#### 3a. `R/utils.R` — shared helpers

Add next to `add_nuisance_pars` (~line 88):

```r
# Name of the operational-time warp parameter (Math/ballistic-time.md).  Like
# the nuisance parameters this is a *trailing* p_type kept out of
# p_types_canonical: eta = 0 is the exact parent model, so a design that never
# mentions it is unchanged and should not be warned about.
.time_warp_par_name <- "eta"

# Append the operational-time warp parameter to a ballistic model's parameter
# machinery.  Call this immediately BEFORE add_nuisance_pars() so the column
# order stays  <model parameters>, eta, pContaminant, pGuess.  The kernels
# resolve eta by name, so the position is cosmetic.
add_time_warp_par <- function(p_types, transform, minmax, exception = NULL) {
  nm <- .time_warp_par_name
  if (!(nm %in% names(p_types))) {
    p_types[[nm]] <- 0
    transform[[nm]] <- "identity"
    # +-5 keeps omega = exp(eta) inside (0.007, 148): wide enough that the bound
    # never binds on real data, tight enough that (1+u)^(omega-1) cannot overflow.
    minmax <- cbind(minmax, c(-5, 5))
    colnames(minmax)[ncol(minmax)] <- nm
  }
  list(p_types = p_types, transform = transform, minmax = minmax,
       exception = exception)
}

# --- the warp itself, vectorised; the exact mirror of emc2tw:: in src/time_warp.h
# eta == 0 short-circuits to the identity so that a design carrying a constant
# eta = 0 reproduces the parent model bit for bit.
.tw_active <- function(eta) !is.na(eta) & eta != 0

.tw_fwd <- function(u, eta) {                       # s = c_eta(u)
  eta <- rep_len(eta, length(u))
  out <- u
  act <- .tw_active(eta) & !is.na(u) & is.finite(u) & u > 0
  if (!any(act)) return(out)
  om <- exp(eta[act]); l1p <- log1p(u[act]); x <- om * l1p
  s <- ifelse(om < 1e-12, l1p, expm1(x) / om)
  s[x > 709] <- Inf
  out[act] <- s
  out
}

.tw_log_jac <- function(u, eta) {                   # log c'_eta(u)
  eta <- rep_len(eta, length(u))
  out <- numeric(length(u))
  act <- .tw_active(eta) & !is.na(u) & is.finite(u) & u > 0
  out[act] <- (exp(eta[act]) - 1) * log1p(u[act])
  out
}

.tw_jac <- function(u, eta) exp(.tw_log_jac(u, eta))

.tw_inv <- function(s, eta) {                       # u = c_eta^{-1}(s)
  eta <- rep_len(eta, length(s))
  out <- s
  act <- .tw_active(eta) & !is.na(s) & is.finite(s) & s > 0
  if (!any(act)) return(out)
  om <- exp(eta[act])
  y <- log1p(om * s[act]) / om
  u <- ifelse(om < 1e-12, expm1(s[act]), expm1(y))
  u[y > 709] <- Inf
  out[act] <- u
  out
}

# eta column of a mapped-parameter matrix; zeros when the model has no warp.
.tw_eta <- function(pars) {
  if (!is.null(colnames(pars)) && .time_warp_par_name %in% colnames(pars))
    pars[, .time_warp_par_name]
  else rep(0, NROW(pars))
}
```

#### 3b. Model constructors — add the parameter

Anchors (verify with
`grep -n 'add_nuisance_pars\|p_types_canonical' R/model_*.R`):

| Constructor | file | `add_nuisance_pars` call | `p_types_canonical` |
|---|---|---|---|
| `BAwD()`   | `R/model_BAwD.R`  | 527 | 545 |
| `BAwDp()`  | `R/model_BAwD.R`  | 661 | 674 |
| `BAwR()`   | `R/model_BAwR.R`  | 288 | 301 |
| `BAwF()`   | `R/model_BAwF.R`  | 309 | 324 |
| BTAwL transient | `R/model_BTAwL.R` | 435 | 441 |
| BTAwL sustained | `R/model_BTAwL.R` | 488 | 494 |
| BTAwL full race | `R/model_BTAwL.R` | 534 | 545 |
| `BAwL()`   | `R/model_LBA.R`   | 799 | 826 (explicit vector, no change) |
| `LBA()`    | `R/model_LBA.R`   | n/a (literal `p_types`, line 158) | 159 (no change) |

For each of `BAwD()`, `BAwDp()`, `BAwF()`, `BAwR()`, and
`.btawl_constructor()`'s three branches, insert immediately
before the existing `.nuis <- add_nuisance_pars(...)` line:

```r
  # Operational-time warp; eta = 0 is exactly this model (Math/ballistic-time.md).
  .tw <- add_time_warp_par(p_types, transform, minmax, exception)
  p_types <- .tw$p_types; transform <- .tw$transform
  minmax <- .tw$minmax; exception <- .tw$exception
```

Then change the canonical list from

```r
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
```

to

```r
    p_types_canonical = setdiff(names(p_types),
                                c(.time_warp_par_name, .nuisance_par_names)),
```

`BAwL()` (`R/model_LBA.R:717`) is conditional:

```r
  # The Erlang clocks run on raw time and the correlated path has its own
  # kernels, so neither composes with the warp yet; withholding eta from
  # p_types makes `eta ~ 1` fail in design() rather than silently do nothing.
  if (erlang_type == "none" && !correlated) {
    .tw <- add_time_warp_par(p_types, transform, minmax, exception)
    p_types <- .tw$p_types; transform <- .tw$transform
    minmax <- .tw$minmax; exception <- .tw$exception
  }
```

`BAwL`'s `p_types_canonical` is an explicit vector
(`c(launch_pars, "B", "A", "t0", "k")`), so it needs no change.

`LBA()` (`R/model_LBA.R:148`) hard-codes everything in one literal. Edit the
three literals so that `eta` lands before `pContaminant`:

```r
    p_types = c("v" = 1, "sv" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0),
                "eta" = 0, "pContaminant" = qnorm(0), "pGuess" = qnorm(0)),
    p_types_canonical = c("v", "sv", "B", "A", "t0"),   # unchanged: eta stays optional
    transform = list(func = c(v = "identity", sv = "exp", B = "exp", A = "exp",
                              t0 = "exp", eta = "identity",
                              pContaminant = "pnorm", pGuess = "pnorm")),
    bound = list(minmax = cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf),
                                A = c(1e-4, Inf), B = c(1e-4, Inf),
                                t0 = c(0.05, Inf), eta = c(-5, 5),
                                pContaminant = c(0.001, 0.999),
                                pGuess = c(0.001, 0.999)),
                 exception = c(A = 0, pContaminant = 0, pGuess = 0)),
```

Do **not** touch `LBA_LogicalRules` (`R/model_LBA.R:265`).

#### 3c. `Ttransform` — warp the reported physical endpoint

Every finite-endpoint model reports `Tmax` (operational) and
`rt_max` (physical). Under the warp only the second moves. In
`R/model_BAwD.R` (BAwD ~line 548, BAwDp ~line 677), `R/model_BAwF.R` (~327),
`R/model_BAwR.R` (~310), and the three `R/model_BTAwL.R` branches
(~444, ~497, ~548), change

```r
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax, ...)
```

to

```r
      # Tmax is the OPERATIONAL endpoint and is invariant under the warp
      # (spec §17); only the physical ceiling moves.
      cbind(pars, b = b, Tmax = Tmax,
            rt_max = pars[, "t0"] + .tw_inv(Tmax, .tw_eta(pars)), ...)
```

`.tw_inv(Inf, eta) == Inf`, so the no-endpoint cases (BAwD `gamma = 1`, BAwR
`kappa = 0`, BTAwL with a sustained channel, BAwDp `lambda = 0`) are unchanged.

**BTAwL `chart = "endpoint"` note.** That chart *samples* `Ttrans` and
back-solves `tau`, so `Ttrans` is by construction the **operational** transient
endpoint once `eta` is free, and `rt_max` is the physical one. The existing
`out[, "Ttrans"] <- Tmax` assignment stays as is. This must be documented (5e)
because it changes the interpretability claim that chart was built for.

---

### Stage 4 — R reference paths

These are not the sampled likelihood (the compiled race path is authoritative;
`log_likelihood()` for BAwD/BAwF/BAwR/BTAwL deliberately `stop()`s), but they
are the documented single-accumulator wrappers and the `emc2.cpp_rfun = FALSE`
reference simulators, and tests use both.

#### 4a. `dXXX` / `pXXX` wrappers

Sixteen functions, all with the same shape:

| function | file:line | | function | file:line |
|---|---|---|---|---|
| `.lba_dfun` / `.lba_pfun` | `R/model_LBA.R:1,16` | | `dBAwD` / `pBAwD` | `R/model_BAwD.R:111,131` |
| `dBAwL` / `pBAwL` | `R/model_LBA.R:331,357` | | `dBAwDp` / `pBAwDp` | `R/model_BAwD.R:162,178` |
| `dBAwF` / `pBAwF` | `R/model_BAwF.R:67,84` | | `dBAwR` / `pBAwR` | `R/model_BAwR.R:48,64` |
| `dBTAwLTransient` / `pBTAwLTransient` | `R/model_BTAwL.R:30,47` | | `dBTAwL` / `pBTAwL` | `R/model_BTAwL.R:81,99` |
| `dBTAwLSustained` / `pBTAwLSustained` | `R/model_BTAwL.R:125,145` | | | |

The pattern is always

```r
  dt  <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ...
  s <- .tw_fwd(dt[ok], eta[ok])
  out[ok] <- dbawr(t = s, ...) * .tw_jac(dt[ok], eta[ok])   # density
  out[ok] <- pbawr(t = s, ...)                              # CDF: no Jacobian
```

Do **not** add an `eta` argument to the `[[Rcpp::export]]`ed `dbawr`/`pbawr`/…
kernels. Keeping the warp in R means `Rcpp::compileAttributes()` never has to
run for this change, and `EMC2:::dbawr()` keeps its meaning as "the parent
kernel at operational time `s`", which the tests rely on as an oracle.

#### 4b. R fallback simulators

`rBAwD`, `rBAwDp` (`R/model_BAwD.R:290,369`), `rBAwF` (`R/model_BAwF.R:140`),
`rBAwR` (`R/model_BAwR.R:118`), `rBAwL` (`R/model_LBA.R:383`), `.lba_rfun`
(`R/model_LBA.R:36`), `.rBTAwLTransient_R`, `.rBTAwL_R`, `.rBTAwLSustained_R`
(`R/model_BTAwL.R:212,265,345`).

Each computes a vector of accumulation-time hits before assembling `dt`. Insert
one line right after that vector exists, e.g. in `rBAwR`:

```r
    hit <- mapply(.bawr_hit_time, V, p[, "b"] - z, p[, "kappa"], p[, "p"])
    hit <- .tw_inv(hit, .tw_eta(p))            # operational -> physical time
    hit[!is.finite(hit) | hit < 0] <- Inf
```

Same idea in the others; make sure it lands **before** `t0` is added.

---

### Stage 5 — Documentation

#### 5a. Parameter tables

Add one row to the roxygen parameter table of every constructor that gained
`eta` (`LBA`, `BAwL`, `BAwD`, `BAwDp`, `BAwF`, `BAwR`, `BTAwL`,
`BTAwLTransient`, `BTAwLSustained`):

```
#' | *eta* | identity | \[-5, 5\] | 0 | | Operational-time warp; `0` is ordinary time. |
```

#### 5b. Shared `@details` paragraph (paste into each of the nine)

```
#' Every ballistic model in EMC2 accepts an optional operational-time warp
#' `eta`, which reparameterises accumulation time as
#' \verb{s = ((1 + u)^omega - 1)/omega} with \verb{omega = exp(eta)} and
#' \verb{u = rt - t0}. The trajectory, the required launch strength, the
#' critical launches, and the intrinsic omission probability are exactly the
#' parent model's, evaluated at `s`; only the mapping from internal to physical
#' time changes, and the density picks up the Jacobian
#' \verb{(1 + u)^(omega - 1)}. `eta = 0` is the parent model exactly, and is the
#' default, so a formula that does not mention `eta` fits the unchanged model.
#'
#' `eta < 0` makes operational time run progressively slower, stretching the
#' right tail (an exponential physical tail becomes a stretched exponential);
#' `eta > 0` speeds it up and reduces right skew. The clock starts at ordinary
#' speed (\verb{c'(0) = 1}), so `eta` is a curvature rather than a rescaling of
#' time, and the reference time is a hard-coded **one second**: `eta` is only
#' interpretable when RTs are in seconds.
#'
#' `eta` is intended as a subject-level parameter shared by every accumulator
#' and condition (`eta ~ 1`). Per-accumulator warps are expressible but change
#' which accumulator wins a given latent race, which is not the mechanism this
#' parameter represents. See `Math/ballistic-time.md`.
```

#### 5c. `@param` note where the warp is withheld

In `BAwL()`'s `@param erlang_type` and `@param correlated`, add: "Clock and
correlated variants do not support the operational-time warp `eta`."

#### 5d. `NEWS.md`

One bullet under the current development heading.

#### 5e. `Math/ballistic-time.md`

Append an "Implementation status" section recording: what was implemented,
the three carve-outs of 2.2, the operational-vs-physical reading of `Tmax`,
`Ttrans` and `rt_max`, and the fact that BAwDp's Jacobian composes as
`m'(c(u)) c'(u)`. Also fix the two typos in the source spec while you are
there: §11 writes `c_\eta'((u)` (doubled paren), and §11's `\eta_j(s)` channel
gains collide notationally with the warp parameter `\eta` — rename the gains to
`g_j(s)` in that section.

#### 5f. Regenerate `man/`

```bash
Rscript -e 'roxygen2::roxygenise(".")'
git status --short man/
```

---

### Stage 6 — Tests

New file `tests/testthat/test-time-warp.R`. Start it with
`skip_model_validation()` only for the slow Monte-Carlo blocks — the fast
invariants (1, 2, 8, 9) should run unconditionally, and one cheap contract
assertion belongs in `tests/testthat/test-model-family-contracts.R`.

Build a table once and loop it, so every model is held to the same standard —
that is the point of the whole design:

```r
tw_models <- list(
  LBA   = list(model = function() LBA(),
               formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, eta ~ 1)),
  BAwL  = list(model = function() BAwL(),  ...),
  BAwD  = list(model = function() BAwD(),  ...),
  BAwDp = list(model = function() BAwDp(), ...),
  BAwF  = list(model = function() BAwF(),  ...),
  BAwR  = list(model = function() BAwR(),  ...),
  BTAwL_T = list(model = function() BTAwLTransient(), ...),
  BTAwL_S = list(model = function() BTAwLSustained(), ...),
  BTAwL   = list(model = function() BTAwL(), ...)
)
```

Required tests, keyed to spec §19:

1. **Parent recovery (the critical one).** For every model in `tw_models`,
   build two designs on the same small data set — one whose model has `eta`
   constant at 0, one built from a constructor call with `eta` absent (e.g.
   compare `BAwR()` against `BAwR()` with `eta` in `constants`) — and assert
   the compiled log-likelihoods are **identical** (`expect_identical` on the
   numeric vector, not `expect_equal`). Also assert `make_data()` with a fixed
   seed produces identical `R`/`rt`.
2. **Warp algebra.** `.tw_inv(.tw_fwd(u, eta), eta) == u` to 1e-12 over a grid
   `u in {1e-6, .01, .1, .5, 1, 3, 10}`, `eta in {-3,-1,-.3,0,.3,1,3}`;
   `.tw_fwd(u, 0) identical u`; `.tw_log_jac(u, 0) identical 0`; the
   `omega -> 0` limits `c(u) -> log1p(u)` and `c^{-1}(s) -> expm1(s)`;
   `.tw_fwd(Inf, eta) == Inf`.
3. **C++/R agreement.** For each model, compare the compiled likelihood
   (`EMC2:::calc_ll_oo`) at a non-zero `eta` against a hand-written R reference
   built only from the exported kernels:
   `log dbawr(s) + log_jac + sum_j log(1 - pbawr_j(s))`. This is the test that
   proves the wrapper composes correctly and is worth more than all the others.
4. **CDF/density consistency.** `numDeriv::grad(F_eta, u) == f_eta(u)` for a
   handful of `u` per model, both signs of `eta`.
5. **Defect invariance.** `pXXX(Inf, pars_eta)` equals `pXXX(Inf, pars_0)` for
   `eta in {-2, -0.5, 0.5, 2}`, every model with `defective_upper_tail`.
6. **Endpoint mapping.** For BAwD (`gamma < 1`), BAwDp (`0 < lambda < 1`),
   BAwF, BAwR, BTAwL-transient: `F_eta` is flat and `f_eta == 0` for
   `u > .tw_inv(Tmax, eta)`, and *not* flat just below it. Also check that
   `Ttransform`'s `rt_max` equals `t0 + .tw_inv(Tmax, eta)`.
7. **Simulation recovery.** Simulate ~2e4 trials from each model with
   `eta = -0.7` and `eta = +0.7` (default C++ simulator, per `AGENTS.md`), then
   KS-compare the empirical accumulation-time CDF against `F_0(c_eta(u))`.
   Repeat under `options(emc2.cpp_rfun = FALSE)` for at least one model to hold
   the R fallback to the same standard.
8. **Race ordering.** With a common `eta`, simulate under a fixed seed with
   `eta = 0` and with `eta = 1` and assert the response identity `R` is
   trial-by-trial identical while `rt` is not (spec §19, "Race ordering").
9. **Censoring / truncation.** Fit-free check: with `UC` finite, the compiled
   censored log-likelihood equals `log S_0(c_eta(UC - t0))` summed over
   accumulators (spec §3), and the truncated likelihood at `LT/UT` matches a
   numerically integrated reference. This is the only test that exercises
   `tw_logS_at_t` and `tw_cdf1`, so do not skip it.
10. **Refusals.** `design(..., model = function() BAwL(erlang_type =
    "local_kill"), formula = list(..., eta ~ 1))` errors; ditto
    `BAwL(correlated = TRUE)`.
11. **Recovery study** (not a unit test) in
    `WorkingTests/time_warp_recovery.R`: simulate a two-subject hierarchical
    BAwR and LBA at `eta = -0.6`, fit, and report `eta` recovery plus its
    posterior correlation with `t0`, `B`, and the launch spread. Expect `eta`
    to trade off with the launch-spread parameter (`sigma`/`sv`) — both control
    right skew. **Use `rt_resolution = NULL` in `make_emc()`**: `make_emc`
    floors RTs to `1/60 s` by default, which biases any parameter that reads
    the density near small `u`, and `eta` reads exactly there.

Finally, extend `tests/testthat/test-model-family-contracts.R`'s first
`test_that` with a loop asserting that every ballistic constructor exposes
`"eta"` in `names(m$p_types)`, that its default is `0`, and that it is *absent*
from `p_types_canonical`.

---

## 6. Semantics contract — what `eta` means in each model

This table is the answer to "consistency between each model in what the time
parameter means". In every row the *only* change is `u -> s = c_eta(u)` plus the
density Jacobian; the third column is the mechanism-specific consequence.

| Model | Trajectory under the warp | What `eta < 0` (slower operational time) does |
|---|---|---|
| LBA | `X = z + V·c(u)` | Stretches the whole RT distribution's right tail; no endpoint to move. |
| BAwL | `X = z e^{-k c(u)} + (V/k)(1 - e^{-k c(u)})` | Exponential physical tail becomes a **stretched exponential** `u^{omega-1} exp(-k u^omega / omega)`; `V_crit = kb` and hence the omission rate are unchanged. |
| BAwD | `X = z + V Q_rho(c(u)) - ell R_{rho,gamma}(c(u))` | Pushes the finite wall `S_max` out in physical time; live/frozen split, `z*`, `Z`, and the frozen quadrature all unchanged. |
| BAwDp | `X = z + V q_lambda(c(u))` | Freeze time `S* = log(1/lambda)/k` is unchanged; the physical response-free gap after `U* = c^{-1}(S*)` moves. |
| BAwF | `X = h_rho(c(u))[z + V c(u)]` | Same: `x_max` and `S_max = x_max/k` fixed; `U_max = c^{-1}(S_max)` moves. The b-free-endpoint property is preserved in operational time. |
| BAwR | `X = z + V c(u) - kappa c(u)^{p+1}/(p+1)` | `S_sat(z)` and `S_max` (which *do* depend on `b`) are unchanged; `U_max = c^{-1}(S_max)`. |
| BTAwL | `X_j = z_j e^{-k c(u)} + V_j H_j(c(u))`, `j in {T, S}` | **One warp shared by both channels** — this is mandatory (spec §11/§12); the Jacobian appears once in the local-race density, never twice. |

Invariant across all of them, and worth asserting in tests:
`F_eta(Inf) = F_0(Inf)`, `V_crit` unchanged, `Tmax` (operational) unchanged,
`rt_max` (physical) `= t0 + c^{-1}(Tmax)`.

---

## 7. Traps

1. **Rebuild everything.** `src/Makevars` tracks no header dependencies. Adding
   `TimeWarpPlan tw` to `ContextForRaceModels` changes the struct layout for
   every TU that includes `race_contract.h`. `rm -f src/*.o src/EMC2.so` or you
   get a silent ODR mismatch and a segfault.
2. **`eta == 0` must short-circuit.** `expm1(log1p(u)) != u` bitwise. Without
   the short-circuit, every existing test that compares log-likelihoods to
   stored values starts failing by ~1e-16 per trial, which compounds.
3. **Jacobian before flooring.** See `tw_d_raw` in 1c. Adding `log_jac` to an
   already-floored `min_ll` is wrong; toggle `floor_raw_log_lik` instead.
4. **`t0` conventions differ between simulators.** `rlba_cpp`/`rbawd_cpp`/
   `rbawf_cpp`/`rbawr_cpp`/`rbawdp_cpp` store bare accumulation time in `dt`
   and let `resolve_race` add `t0`; `rbawl_cpp_impl` and `rbta_wl_cpp` add `t0`
   themselves. Warping a value that already contains `t0` warps the
   non-decision time too and is silently wrong (it will pass a KS test at small
   `t0` and fail at large `t0`).
5. **`cdf1(+Inf)` must reach the base unwarped.** Several kernels document that
   `t == Inf` returns `F_max`, not `1`. `t0 + fwd(Inf, eta)` is `Inf`, so this
   works either way, but the explicit `!R_FINITE(u)` early return makes it
   obvious and avoids `Inf - Inf` when `t0` is also infinite.
6. **`tw_logS_at_t` may be handed compacted columns.**
   `race_endpoint_prepare_groups` (`src/race_integrands.cpp:44`) dedupes trials
   by their full parameter block and passes `compact_col_ptrs` (length exactly
   `n_par`) with a different `n_rows_total`. Build the warped `t0` buffer from
   whatever `cols`/`n_rows_total` you were given, never from an outer variable.
   `eta` is part of the dedup key automatically, so grouping stays correct.
7. **Do not copy pointers past `n_par`.** Padding slots must be `nullptr`, not
   copies of out-of-range source entries.
8. **BTAwL `Ttrans` changes meaning.** With `eta` free it is the *operational*
   transient endpoint. Say so in the docs; a user who chose that chart for
   interpretability deserves to know.
9. **`Tmax` vs `rt_max`.** Only `rt_max` warps. Warping `Tmax` too would
   double-apply the transform wherever `Tmax` feeds a kernel (BTAwL's endpoint
   chart back-solves `tau` from it).
10. **`eta` collides with a sampler-internal name.** `variant_factor.R`,
    `variant_SEM.R`, and `variant_infnt_factor.R` store factor scores in
    `sampler$samples$eta`. That is a different list slot from `samples$alpha`,
    so there is no functional collision and `"eta"` is not a valid
    `selection = ` value — but note it in a code comment so the next reader does
    not spend an hour on it.
11. **Recovery studies need `rt_resolution = NULL`.** `make_emc()` floors RTs to
    `1/60 s` by default. `eta` is largely identified by the shape of the density
    at short `u`, exactly where the flooring bites.
12. **`eta` trades off with launch spread.** `sv`/`sigma` and `eta` both control
    right skew. Expect a strong posterior correlation and weak identification
    when both are free with few trials; this is a property of the model, not a
    bug, but say it in the docs so nobody files it as one.

---

## 8. Acceptance criteria

- [ ] `rm -f src/*.o && R CMD INSTALL` succeeds with no new warnings.
- [ ] Full `testthat::test_local()` output is identical to
      `/tmp/tw_baseline.txt` except for the new `test-time-warp.R` block.
- [ ] For all nine models: log-likelihood with `eta` constant `0` is
      `identical()` to the pre-change build's log-likelihood on the same data
      and parameters.
- [ ] For all nine models: `make_data()` under a fixed seed with `eta = 0` is
      `identical()` to the pre-change build's output.
- [ ] `test-time-warp.R` items 1–10 pass.
- [ ] `WorkingTests/time_warp_recovery.R` recovers `eta` within its 95% CI for
      both models.
- [ ] `roxygen2::roxygenise()` run and `man/` committed.
- [ ] `Math/ballistic-time.md` gained its implementation-status section.

---

## 9. Files touched (checklist)

**New**
- `src/time_warp.h`
- `src/time_warp.cpp`
- `tests/testthat/test-time-warp.R`
- `WorkingTests/time_warp_recovery.R`

**C++ modified**
- `src/race_contract.h` — `TimeWarpPlan`, `ContextForRaceModels::tw`
- `src/race_dispatch.cpp` — 9 `tw.supported` lines
- `src/particle_ll.cpp` — 1 include, 2 `configure_time_warp_context` calls
- `src/model_rng.cpp` — 1 include, 7 simulators (9 assignment sites)
- `src/col_registry.h` — comment only

**R modified**
- `R/utils.R` — `.time_warp_par_name`, `add_time_warp_par`, `.tw_*` helpers
- `R/model_LBA.R` — `LBA()`, `BAwL()`, `.lba_dfun`/`.lba_pfun`,
  `dBAwL`/`pBAwL`, `.lba_rfun`, `rBAwL`
- `R/model_BAwD.R` — `BAwD()`, `BAwDp()`, `dBAwD`/`pBAwD`,
  `dBAwDp`/`pBAwDp`, `rBAwD`, `rBAwDp`, both `Ttransform`s
- `R/model_BAwF.R` — `BAwF()`, `dBAwF`/`pBAwF`, `rBAwF`, `Ttransform`
- `R/model_BAwR.R` — `BAwR()`, `dBAwR`/`pBAwR`, `rBAwR`, `Ttransform`
- `R/model_BTAwL.R` — `.btawl_constructor()` × 3 branches, the three
  `dBTAwL*`/`pBTAwL*` wrappers, the three `.rBTAwL*_R` fallbacks, three
  `Ttransform`s
- `tests/testthat/test-model-family-contracts.R` — contract assertion

**Docs**
- `man/*.Rd` (regenerated), `NEWS.md`, `Math/ballistic-time.md`

**Not touched, deliberately**
- `src/model_BAwD.cpp`, `src/model_BAwF.cpp`, `src/model_BAwR.cpp`,
  `src/model_BAwL.cpp`, `src/model_BTAwL.cpp`, `src/model_LBA.h`,
  `src/correlated_likelihood.cpp`, `src/logicalrules_likelihood.cpp`,
  `src/marginal_likelihood.cpp`, `src/RcppExports.cpp`, `NAMESPACE`.
  If you find yourself editing a `model_BAw*.cpp`, stop: the design says the
  kernels never learn about the warp, and something has gone wrong upstream.

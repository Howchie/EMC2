# Fix: Erlang Process t0 Timing Bug (rdmswtn, rdmgbm, bawl)

## Context

Race model likelihoods use the trick `tt = rt - t0` in adapter wrappers, passing `tt` to core functions. This correctly excludes the zero-density `rt < t0` region for the EAM accumulation process. However, when erlang guess/kill processes are added (local mixture variants), the same `tt` is passed to `erlang_log_surv`/`erlang_log_pdf` inside the core functions — incorrectly making the erlang clock start at `t0` instead of `t=0`. This biases erlang rate estimates.

A second related bug: when `tt <= 0` (`rt < t0`), adapters unconditionally return `min_ll`/0, even when erlang processes are active. But erlang responses can legitimately occur before `t0` (e.g., fast guess responses), so the `rt < t0` region must not be discarded when erlang is active.

## Fix Strategy

**Least-destructive approach**: Add `double t0 = 0.0` as a new trailing default parameter to each core likelihood function. When `t0 = 0` (all existing internal callers), behavior is identical. Adapters in `utils.h` then pass raw `rt[i]` as first argument and `t0_[i]` explicitly, so the core functions split:
- EAM parts use `tt = t - t0` internally  
- Erlang parts use raw `t` directly

## Critical Files

- `src/model_RDM.h` — core RDM/Wald/GBM likelihoods
- `src/model_LBA.h` — core BAwL likelihoods  
- `src/utils.h` — adapter wrappers (scalar, batch-raw, logS_at_t) for all 3 models

## Implementation Steps

### 1. Core functions in `src/model_RDM.h`

Add `double t0 = 0.0` parameter to each, compute `const double t_eam = t - t0;` at the top, then split usage:

| Function | Lines (approx) | EAM uses | Erlang uses |
|----------|---------------|----------|-------------|
| `dwald` | ~314 | `t_eam` for `var`, `mu_new`, `log(t_eam)` | `t` in `erlang_log_surv`, `erlang_log_pdf` |
| `pwald` | ~386 | `t_eam` for `st`, `c1`, `c2`, barrier calcs | `t` in `erlang_log_surv` |
| `dgbm` | ~593 | `t_eam` for `var`, `mu_new` | `t` in `erlang_log_surv`, `erlang_log_pdf` |
| `pgbm` | ~650 | `t_eam` for `st`, `c1`, `c2` | `t` in `erlang_log_surv` |
| `dswtn_core` | ~987 | `t_eam` for `tv`, `den_common` etc. | `t` in `erlang_log_surv`, `erlang_log_pdf` |
| `drdmswtn` | ~1282 | `t_eam` guard + passes to `dwald`/`dswtn_core` with `t0` | — |
| `prdmswtn` | ~1323 | same | — |
| `drdmswtn_local_combo` | ~932 | `t_eam` for EAM density call | `t` passed to `local_combo_response_pdf` |
| `prdmswtn_local_combo` | ~945 | `t_eam` lower bound for CDF | `t` in erlang kernel |
| `dgbm_local_combo` | ~963 | same pattern | same |
| `pgbm_local_combo` | ~972 | same pattern | same |

**Early exit restructuring in each core function:**
```cpp
// OLD top guard:
if (var <= FPM_EPSILON || t <= FPM_EPSILON ...) return 0;

// NEW:
const double t_eam = t - t0;
if (t <= 0.0) return log_out ? R_NegInf : 0.0;  // physically impossible raw rt
if (t_eam <= 0.0) {
  // EAM contribution is 0; erlang guess can still contribute
  if (!guess || k <= 0.0) return log_out ? R_NegInf : 0.0;
  // S_R(0) = 1 (EAM not started), so:
  const double log_fG = erlang_log_pdf(t, k, kill_shape);
  return log_out ? log_fG : std::exp(log_fG);
}
// ... rest of function unchanged, using t_eam for EAM, t for erlang
```

Recursive internal calls (e.g., `pwald` called from `dwald`'s guess branch) must also forward `t0`.

### 2. Core functions in `src/model_LBA.h`

Add `double t0 = 0.0` to `dkilledleakyba_norm` (~line 262) and `pkilledleakyba_norm` (~line 289).

- `dleakyba_norm(t, ...)` → `dleakyba_norm(t_eam, ...)`
- `pleakyba_norm(t, ...)` → `pleakyba_norm(t_eam, ...)`
- `erlang_log_surv(t, lambda_k, ks)` and `erlang_log_pdf(t, lambda_g, ks)` — keep raw `t`

For the `pkilledleakyba_norm` GL20 combined kill+guess branch, the integration variable `u` runs from `0` to raw `t`, and the integrand uses `u_eam = max(0.0, u - t0)` for EAM parts:
```cpp
const double u  = 0.5 * t * (nodes[j] + 1.0);   // raw time
const double u_eam = std::max(0.0, u - t0);       // EAM time
const double fR = dleakyba_norm(u_eam, ...);
const double SR = 1.0 - pleakyba_norm(u_eam, ...);
// SG, SK, fG still use raw u
```

When `t_eam <= 0` (entire range is pre-EAM):
- `use_kill` only: return 0 (no EAM hits, kill only suppresses)
- `use_guess` only: return `finish(log1p(-exp(erlang_log_surv(t, lambda_g, ks))))`
- Both: return `finish(log1p(-exp(erlang_log_surv(t, lg, ks) + erlang_log_surv(t, lk, ks))))`

### 3. Adapters in `src/utils.h`

**6 scalar functions** (`dbawl_scalar`, `pbawl_scalar`, `drdmswtn_scalar`, `prdmswtn_scalar`, `drdmgbm_scalar`, `prdmgbm_scalar`):
```cpp
// OLD:
const double tt = t - par[3];  // or par[4] for bawl
if (tt <= 0.0) return 0.0;
core_fn(tt, ...);

// NEW:
const double t0_val = par[3];  // (or par[4])
const double tt = t - t0_val;
const bool erl = (lg > 1e-12 || lk > 1e-12);  // must read lg/lk before guard
if (tt <= 0.0 && !erl) return 0.0;
if (t <= 0.0) return 0.0;
core_fn(t, ..., /*t0=*/t0_val);  // pass raw t and t0 separately
```

**6 batch-raw functions** (`dbawl_raw`, `pbawl_raw`, `drdmswtn_raw`, `prdmswtn_raw`, `drdmgbm_raw`, `prdmgbm_raw`):

Move `lg`/`lk` reads before the early-exit guard. Same restructuring:
```cpp
const double t0_i = t0_[i];
const double tt   = rt[i] - t0_i;
const double lg   = ...;
const double lk   = ...;
const bool erl    = (lg > 1e-12 || lk > 1e-12);
if (tt <= 0.0 && !erl) { out[i] = min_ll; continue; }
if (rt[i] <= 0.0)      { out[i] = min_ll; continue; }
// ... call core_fn(rt[i], ..., t0_i)
```

Note: the `k_use <= 0` fast path (`dwald_k0`, `pwald_k0`) only fires when no erlang, so its internal `if (tt <= 0)` guard is still correct there.

**3 logS_at_t functions** (`bawl_logS_at_t`, `rdmswtn_logS_at_t`, `rdmgbm_logS_at_t`):

These compute log survivor at a scalar truncation time `t`. When `t < t0[r]` and erlang is active, the accumulator hasn't started but erlang has been running:
```cpp
const double t0_r = t0_[r];
const double tt   = t - t0_r;
const double lg   = ...;
const double lk   = ...;
const bool erl    = (lg > 1e-12 || lk > 1e-12);

if (tt <= 0.0) {
  if (!erl) continue;  // EAM not started, no erlang → log-survivor contribution = 0
  // EAM not started → S_R = 1; erlang processes running since t=0
  // log-survivor contribution = erlang_log_surv(t, lambda_total, ks)
  if (is_local_kill_guess && lg > 1e-12 && lk > 1e-12) {
    logS += erlang_log_surv(t, lg, ks) + erlang_log_surv(t, lk, ks);
  } else {
    const double lam = (is_guess_type ? lg : 0.0) + lk;
    if (lam > 1e-12) logS += erlang_log_surv(t, lam, ks);
  }
  continue;
}
// tt > 0: call pcore_fn(t, ..., t0_r) — pass raw t and t0_r
```

## Correctness Invariants

| Condition | Expected PDF output |
|-----------|-------------------|
| `rt <= 0` | `min_ll` (impossible) |
| `0 < rt < t0`, no erlang | `min_ll` (EAM zero) |
| `0 < rt < t0`, guess only | `erlang_pdf(rt, lambda_g)` (S_R=1) |
| `0 < rt < t0`, kill only | `min_ll` (kill only suppresses, no density source) |
| `0 < rt < t0`, both | `erlang_pdf(rt, lg) * erlang_surv(rt, lk)` |
| `rt >= t0`, no erlang | unchanged (existing behavior) |
| `rt >= t0`, with erlang | erlang uses raw `rt`, EAM uses `rt - t0` |

## Backward Compatibility

When `lambda_g = 0` AND `lambda_k = 0`: the `erl` flag is false, so the adapter early-exit fires identically to before. The core functions receive `t0 = 0.0` by default from all non-adapter callers, making `t_eam = t`, preserving all existing behavior exactly.

## Verification

1. Run full existing test suite — all tests must pass unchanged (no erlang or `t0=0` cases dominate)
2. Add `tests/testthat/test-erlang-t0.R` with:
   - Identity test: `dwald(rt, ..., k=0.5, t0=0)` == `dwald(rt, ..., k=0.5)` 
   - Guess density at `rt < t0`: verify equals `erlang_pdf(rt, lg)` analytically
   - PDF/CDF numerical derivative consistency at `rt < t0` and `rt > t0`
   - Mass normalization: `prdmswtn(Inf, ..., erlang)` ≤ 1
   - logS_at_t with `t < t0` and erlang: verify equals `erlang_log_surv(t, lambda, n)`
3. Build with `devtools::document()` + `devtools::install()` or `R CMD INSTALL`

## Implementation Order

1. `dwald` + `pwald` in `model_RDM.h` → test with `dWald`/`pWald` R exports
2. `dgbm` + `pgbm` → test with `dGBMspv`/`pGBMspv`
3. `dswtn_core`, `drdmswtn`, `prdmswtn`, `*_local_combo` → test with `dSWTNspv`/`pSWTNspv`
4. `dkilledleakyba_norm` + `pkilledleakyba_norm` in `model_LBA.h` → test with `dBAwL`
5. All adapters in `utils.h` — scalar, batch-raw, logS_at_t for all 3 models
6. Run full test suite + add new targeted tests

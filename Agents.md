## 1. EMC2 Philosophy

This document captures the optimisation and architectural philosophy that guides ongoing development of EMC2. It is intended for contributors (human or AI) who are adding new models, refactoring existing code, or reviewing changes.

The overarching constraints are absolute:
1. **Speed can never change results**: Vectorising a call or switching to a faster kernel that produces an identical result is good; changing user-facing semantics or model behaviour is not acceptable.
2. **Leanness**: When changing the behavior of a function or model, do not preserve backward compatibility shims by default unless specifically instructed to do so. Dead code, unused solvers, and obsolete interfaces should be removed decisively.
3. **Modular separation**: Model implementations must be strictly decoupled from likelihood orchestration engines. Analytic cores live in `.cpp` files; model geometry and public declarations live in model headers; kernel column contracts live in `col_registry.h`; shared contracts live in `race_contract.h`.
4. **Minimalist comments**: Comments should describe specific functionality, never documenting development history, alternative past attempts, or rationale logs. This is particularly true for user-facing roxygen comments — never use comments as a train of thought or changelog under any circumstances.

# Repository requirements

- Every model constructor in this package (including race, diffusion,
  discrete-choice, signal-detection, and MRI models) MUST use a compiled C++
  simulator by default. A pure-R generator may exist only as an explicit
  opt-out/reference path (`options(emc2.cpp_rfun = FALSE)`).
- Adding or changing a model requires an exported C++ simulator entry point,
  its Rcpp binding, and a regression test that exercises the default path.
- Keep comments short and local; do not duplicate implementation details in
  prose when the code and tests already make the behavior clear.
- All models must document their parameter table in the Details section in a consistent fashion;
  when adding new models ensure the documentation is updated and made consistent with existing models.

---

## 2. Compiler-Level Optimisation (C++ / Rcpp)

### 2.1 Build flags

`src/Makevars` (and `Makevars.win` / `Makevars.ucrt`) set the following flags, which must be kept:

```makefile
override CXXFLAGS += -O3 $(MARCH_CXXFLAGS) -ffast-math -fno-math-errno -Igsl_bundled -DUSE_FAST_PNORM -DR_NO_REMAP
```

- **`-O3`**: full inlining, loop unrolling, and auto-vectorisation. R's default `-O2` is insufficient for tight likelihood loops.
- **`$(MARCH_CXXFLAGS)`**: the best `-march` target for the target machine, selected by `./configure` probing (Linux x86_64). It emits `-march=<level>` when the compiler's `-march=native` resolution is a recognised but stale level (e.g. GCC 11 mapping AMD EPYC Genoa to znver3), and otherwise falls back to `-march=native`. Column-major data layouts and `#pragma omp simd` loops are written to exploit SIMD.
- **Parallel compilation**: `src/Makevars` self-injects a safe `-jN` into `MAKEFLAGS` unless the caller already supplied parallelism, so `R CMD INSTALL .` gets parallel compilation without manually exporting `MAKEFLAGS=-jN`.
- **`-ffast-math`**: allows the compiler to reorder and contract floating-point operations, enabling auto-vectorisation of FP loops. This is safe here because MCMC log-likelihood differences matter; per-sample ULP precision does not.
- **`-DUSE_FAST_PNORM`**: activates the Hart rational-approximation normal CDF in `pnorm_std()` (see §3).
- **`-DR_NO_REMAP`**: required for clean builds, ensures R API macros do not pollute standard C++ symbols.
- **`-Igsl_bundled`**: ensures bundled GSL headers are used across Linux and Windows toolchains to eliminate external build dependencies.

### 2.2 `-ffast-math` incompatibility guard

`-ffast-math` enables `-ffinite-math-only`, which breaks `std::isfinite`, `std::isinf`, and `std::isnan` because the compiler assumes all values are finite. **Never use the standard library versions directly.** Instead use the macros defined in `src/utility_functions.h`:

```cpp
emc2_isfinite(x)   // → R_FINITE(x)
emc2_isinf(x)      // → !R_FINITE(x) && !ISNAN(x)
emc2_isnan(x)      // → ISNAN(x)
```

These delegate to Rinternals macros that are unconditionally correct under `-ffast-math`. Apply this rule to any new C++ code in the package.

### 2.3 Header hygiene

- Prefer `<Rcpp.h>` over `<RcppArmadillo.h>` whenever Armadillo matrix types (`arma::mat`, etc.) are not directly required. Gratuitous Armadillo includes inflate compilation times and memory footprint.
- Never place `using namespace Rcpp;` or `using namespace std;` in header files.

---

## 3. Fast Mathematical Approximations & Reusable Numerics

### 3.1 `pnorm_std` & Fast Distribution Wrappers

`src/utility_functions.h` (and `src/wald_functions.h`) provide fast mathematical approximations:
- `pnorm_std(x, lower, log_p)`: Hart 7-term rational polynomial approximation (when `USE_FAST_PNORM` is defined), 3–8× faster than R's implementation and accurate to machine precision for practical ranges.
- `dnorm_std(x, log_p)`: fast normal PDF without R's wrapper overhead.
- `qnorm_std(p, lower, log_p)`: fast normal quantile approximation.
- `plnorm_std(x, meanlog, sdlog, lower, log_p)` and `dlnorm_std(x, meanlog, sdlog, log_p)`: fast log-normal CDF/PDF routing through `pnorm_std`.

**Policy**: All model headers and likelihood helpers must use `pnorm_std` / `dnorm_std` / `plnorm_std` / `dlnorm_std`. Direct calls to `R::pnorm`, `R::dnorm`, `R::plnorm`, `R::dlnorm` in inner loops are a red flag in code review.

### 3.2 Cancellation-Safe Arithmetic (`src/composite_functions.h`)

Use dedicated numerically stable functions when handling probabilities and logs:
- `std::log1p(-cdf)` instead of `std::log(1.0 - cdf)` for log-survivors when `cdf ≈ 1`.
- `std::expm1(x)` instead of `std::exp(x) - 1.0` when `x ≈ 0`.
- `log_sum_exp(a, b)`: stable log-domain addition $\log(e^a + e^b)$.
- `log1m(p)`: stable computation of $\log(1 - p)$.
- `log1p_exp(x)`: stable computation of $\log(1 + e^x)$.

Compute likelihoods in probability space when safe for speed, and divert to log-space where there is genuine risk of numerical underflow/cancellation.

### 3.3 Shared Callable Quadrature Templates (`src/quad_templates.h`)

Do not write ad-hoc numerical quadrature loops in model files. Use the shared templates:
- `integrate_positive_drift_quad(mu, sv, kernel_fn, n_nodes)`: Gauss-Legendre quadrature for truncated normal positive drift distributions in natural scale.
- `integrate_positive_drift_quad_log(mu, sv, log_kernel_fn, n_nodes)`: log-domain Gauss-Legendre quadrature with `log_sum_exp` accumulation for severe tails.
- `integrate_density_gl20_finite(t_upper, density_fn)`: fixed 20-node Gauss-Legendre finite-interval integral.
- `integrate_density_gl20_infinite(rate_scale, density_fn)`: fixed 20-node mapped infinite-interval integral.
- `integrate_density_adaptive_finite(t_upper, density_fn)` and `integrate_density_adaptive_infinite(rate_scale, density_fn)`: adaptive GSL quadrature fallback.
- `bawd_log_gl(log_f, a, b, n)` and `bawd_log_gl_split(log_f, a, b, mid, n)`: split-panel Gauss-Legendre quadrature for ballistic decay kernels.

### 3.4 Contaminant Mixtures (`src/contaminant_mixture.h`)

EMC2 standardizes contaminant arithmetic into a nested stick-breaking mixture:
- `pContaminant` ($p_C$): omission mixture contributing probability mass only at $t = +\infty$.
- `pGuess` ($p_G$): uniform outlier mixture contributing flat density over the valid response window $[L_G, U_G]$.

Use the unified mixture helpers exclusively:
```cpp
mix_contaminants(ll_proc, pC, pG, log_g, is_omission)
mix_contaminants_rt(ll_proc, pC, pG, guess_kernel, rt, R_known)
```
Do not duplicate contaminant logic inside individual model likelihoods.

### 3.5 Timer Helpers & Clock Conversions (`src/timer_helpers.h`)

For models with Erlang guess/kill timers (BAwL, RDM, etc.), the C++ likelihood receives raw sampled timer means ($m_G, m_K$). Use `src/timer_helpers.h` to convert means to rate parameters:
- `erlang_lambda_from_mean(mean, kill_shape)`: shape-dependent rate conversion (Erlang-1: $1/\mu$; Erlang-2: $2/\mu$; EMIX: $1/\mu$).
- `erlang_omega_for_shape(kill_shape, par, omega_index)`: mixture weight resolution.
- `timed_lambda_dispatch(ctx, lambda_g, lambda_k)`: handles local guess vs. local kill vs. combination dispatch.

---

## 4. Modular Architecture & Contract System

The C++ backend is organized into decoupled modules with explicit contracts:

```
                  ┌─────────────────────────────────────┐
                  │          calc_ll_oo Engine          │
                  │        (src/particle_ll.cpp)        │
                  └──────┬───────────────────────┬──────┘
                         │                       │
         ┌───────────────▼──────────────┐ ┌──────▼────────────────────────┐
         │     Race Likelihood Engine   │ │     DDM / BOU Engine          │
         │   (c_log_likelihood_race)    │ │ (c_log_likelihood_DDM_pt)     │
         └───────────────┬──────────────┘ └──────┬────────────────────────┘
                         │                       │
         ┌───────────────▼──────────────┐ ┌──────▼────────────────────────┐
         │      RaceModelAdapter        │ │         DDMAdapter            │
         │   (src/race_dispatch.cpp)    │ │   (src/likelihood_ddm.cpp)    │
         └───────────────┬──────────────┘ └──────┬────────────────────────┘
                         │                       │
         ┌───────────────▼───────────────────────▼────────────────────────┐
         │              Column Registry (src/col_registry.h)              │
         │         Single source of truth for parameter layouts           │
         └────────────────────────────────────────────────────────────────┘
```

### 4.1 Parameter Column Registry (`src/col_registry.h`)

`src/col_registry.h` is the **single source of truth** for parameter column ordering expected by all C++ batch and scalar kernels.

Every model family must define:
1. An `enum : int` specifying the exact 0-indexed column order.
2. Required columns listed before `N_REQ`; optional columns (e.g. `mG`, `mK`, `omega`, `Binf`, `tau`, `pw`) listed after `N_REQ`.
3. A `ColSpec spec()` returning `{names_array, N_REQ, "LABEL"}`.

```cpp
namespace emc2col {
  namespace lba {
    enum : int { v = 0, sv, B, A, t0, N_REQ };
    inline ColSpec spec() {
      static const char* n[] = {"v", "sv", "B", "A", "t0"};
      return {n, N_REQ, "LBA"};
    }
  }
}
```

At adapter configuration time, `validate_col_prefix(col_names, spec)` is called once per likelihood execution. If R's `p_types` order does not match C++'s expected column layout, it errors immediately, preventing silent parameter misalignment.

**Rule**: Kernels must gate every dereference of optional columns ($j \ge \text{N\_REQ}$) on the corresponding context flag (`ctx->kill_active`, `ctx->kill_shape == 3`, etc.).

### 4.2 Race Adapter Contract (`src/race_contract.h` & `src/race_dispatch.h`)

Race accumulator models interact with the central race engine via `RaceModelAdapter`:

```cpp
struct RaceModelAdapter {
  RacePdf1Fun pdf1_ptr = nullptr;         // scalar PDF for GSL / censoring
  RaceCdf1Fun cdf1_ptr = nullptr;         // scalar CDF for GSL / censoring
  RaceRawFun model_dfun_raw = nullptr;    // fast-path batch log-density
  RaceRawFun model_pfun_raw = nullptr;    // fast-path batch log-survivor
  RaceLogSAtTFun logS_at_t_ptr = nullptr; // batch log-survivor at scalar t
  emc2col::ColSpec col_spec{};            // column layout contract
  ContextForRaceModels ctx;               // model metadata and solve caches
};
```

#### Function Pointer Signatures

1. **`RaceRawFun`** (Fast-path batch raw kernel over all trials):
   ```cpp
   typedef void (*RaceRawFun)(const double* rt,
                              const double* const* cols, // array of column pointers
                              int n_rows,
                              const int* mask,           // 1 = compute, 0 = skip
                              const int* isok,           // row parameter validity
                              double* out,               // write log-density / log-survivor
                              double min_ll,
                              void* ctx);
   ```
2. **`RacePdf1Fun` / `RaceCdf1Fun`** (Scalar PDF/CDF for a single accumulator):
   ```cpp
   typedef double (*RacePdf1Fun)(double rt, const double* par, void* ctx);
   typedef double (*RaceCdf1Fun)(double rt, const double* par, void* ctx);
   ```
3. **`RaceLogSAtTFun`** (Batch log-survivor at scalar $t$ for truncation normalization):
   ```cpp
   typedef void (*RaceLogSAtTFun)(double t, const double* const* cols,
                                  int n_rows_total, int n_lR, int n_par,
                                  const int* trunc_mask, int n_unique_trials,
                                  const int* isok_all, void* ctx, double* logS_out);
   ```

#### Model Registration

All race model variants are registered in `src/race_dispatch.cpp` inside `resolve_race_model_adapter(type_std, caller)`. When adding a new race model:
1. Implement the analytic core and wrappers in `src/model_<NAME>.h` and `src/model_<NAME>.cpp`.
2. Register column definitions in `src/col_registry.h`.
3. Wire up the adapter function pointers, `ColSpec`, and context flags in `src/race_dispatch.cpp`.

#### Raw-Log Helpers (`src/race_contract.h`)

Batch raw kernels must use the shared helpers for consistent log-flooring:
- `raw_floor_log_lik(ctx)`: returns true if raw values should be floored at `min_ll`.
- `raw_log_zero(min_ll, floor_raw)`: returns `min_ll` when floored, `R_NegInf` otherwise.
- `raw_log_value(log_x, min_ll, floor_raw)`: clamps valid finite `log_x` to $\ge \text{min\_ll}$ when floored.

### 4.3 Two-Boundary (DDM / BOU) Contract (`src/likelihood_ddm.h`)

Two-boundary models (Wiener DDM, Bounded Ornstein-Uhlenbeck) share the decoupled engine in `src/likelihood_ddm.cpp`:
- `DDMAdapter`: holds `d_raw`, `p_raw`, `col_spec`, `endpoint_cdf_cache`, and `ContextForDDMModels ctx`.
- `DDMRawFun`:
  ```cpp
  using DDMRawFun = void (*)(const double* rts, const int* Rs,
                             const double* const* cols, int n_rows,
                             const int* mask, const int* is_ok,
                             double* out, double min_ll,
                             ContextForDDMModels* ctx);
  ```
- `c_log_likelihood_DDM_pt`: handles all truncation, censoring, go/no-go, and contaminant mixing without heap allocations, delegating exclusively to `d_raw` / `p_raw`.

### 4.4 Specialized Likelihood Modules

- **`src/marginal_likelihood.h` / `.cpp`**: Integrates $t_0$ out of the complete subject likelihood via Gauss-Legendre quadrature on the log-$t_0$ axis (`calc_ll_oo_marginal`). Treats the underlying likelihood evaluator as a black box.
- **`src/correlated_likelihood.h` / `.cpp`**: Shared single-factor drift correlation (`drift_factor.h`) and Gaussian finishing-time copula likelihoods.
- **`src/logicalrules_likelihood.h` / `.cpp`**: Logical rules and capacity-counter racing architectures.
- **`src/model_SS_adapters.h` / `src/ss_raw.h`**: Stop-signal race likelihoods (SSEXG, SSRDEX).
- **`src/race_integrands.h` / `.cpp`**: GSL scalar race integrands and `RaceEndpointGroupCache`.

---

## 5. Memory Model, Pointers, and Data Layout

### 5.1 Column-Major Indexing & Pointer Array Conventions

`ParamTable` maintains an underlying column-major matrix `base(n_trials, p)`.

#### Batch Kernels: Array of Column Pointers (`const double* const* cols`)
In all batch raw kernels (`RaceRawFun`, `DDMRawFun`, `RaceLogSAtTFun`):
- `cols` is passed as a `const double* const*` array of column pointers.
- Column $j$ is accessed via `cols[j]`, which is a contiguous pointer (`const double*`) to the $N$ rows of that parameter.
- Row $i$ of parameter $j$ is accessed as `cols[j][i]`.
- Parameter column indices are referenced via `emc2col::<model>::<param>` enums:

```cpp
void dmyaccumulator_raw(const double* rt, const double* const* cols, int n_rows,
                        const int* mask, const int* isok, double* out,
                        double min_ll, void* ctx_) {
  const double* v  = cols[emc2col::myaccumulator::v];
  const double* B  = cols[emc2col::myaccumulator::B];
  const double* A  = cols[emc2col::myaccumulator::A];
  const double* t0 = cols[emc2col::myaccumulator::t0];

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // Sequential cache-hot access: v[i], B[i], A[i], t0[i]
  }
}
```

#### Scalar Callbacks: Single-Row Pointer (`const double* par`)
In scalar callbacks (`RacePdf1Fun`, `RaceCdf1Fun`):
- `par` is a `const double*` pointing to the parameter array for a single accumulator.
- Parameter $j$ is accessed as `par[emc2col::<model>::<param>]`.

#### GSL Integration Parameters (`gsl_race_params_scalar`)
- In `gsl_race_params_scalar`, `pars` is a contiguous row-major flattened array of shape `(n_lR, n_par)`.
- Accumulator $a \in [0, n\_lR)$ begins at `pars + a * n_par`.

### 5.2 Planned ParamTable Evaluation

`ParamTable` separates metadata resolution from numerical evaluation:
- Invariant design masks and transform specifications are resolved once per likelihood call.
- `fill_from_particle_row_planned()` fills only modified parameter columns directly via `std::fill` and leaves invariant natural-scale columns untouched.
- `update_pt_only()` uses a planned fast lane that avoids re-allocating string sets or converting `STRSXP` during particle iterations.

### 5.3 Zero Heap Allocation Inside Particle Loops

`ModelSharedState` (and `RaceSharedState`, `DDMSharedState`) is constructed **once per dataset** outside the particle loop. It pre-allocates:
- Pre-read censoring/truncation bound vectors (`LT_vec`, `UT_vec`, `LC_vec`, `UC_vec`).
- Pre-computed finite trial masks (`finite_mask`, `finite_mask_int`).
- Pre-partitioned index vectors for finite vs. non-finite unique trials (`finite_unique_idx`, `other_unique_idx`).
- Pre-allocated scratch buffers (`res_buf`, `idx_win`, `idx_loss`, `ok_int_buf`, `alt_res_buf`, `lF_LC_1_buf`, etc.).
- Pre-resolved contaminant and guess column indices (`pc_col`, `pg_col`, `GuessKernel`).

**Policy**: No `Rcpp::NumericVector`, `std::vector` allocations, or dataframe attribute queries may occur inside per-particle likelihood loops. All scratch memory must come from `ModelSharedState`.

### 5.4 Bulk Memory & SIMD Primitives

- Use `std::copy`, `std::fill`, and `std::fill_n` for bulk array operations (compilers emit vectorised `memcpy`/`memset`).
- Use `std::unordered_map` for setup-time name lookups (`O(1)` average); never perform linear string scans in inner loops.
- Annotate pure reduction loops over trials with `#pragma omp simd reduction(+:total_ll)` when iterations are data-independent and alias-free.

---

## 6. Caching Strategies & State Lifecycle

| Cache Level | Object | Scope | Lifetime / Reset Policy |
|---|---|---|---|
| **Per-Dataset** | `ModelSharedState` | Trial bounds, finite masks, scratch buffers, guess windows | Built once per dataset outside particle loop. |
| **Per-Dataset** | `TrendRuntime` / `ParamTable` | Design plans, filtered transform specs, fill plans | Resolved once at likelihood setup. |
| **Per-Particle** | `fperace::SolveCache`, `rlf::SolveCache`, `fpebou::SolveCache` | Fokker-Planck / PDE finite-difference marches | `std::shared_ptr` in context; cleared once per particle via `cache->clear()`. Amortizes PDE marches across trials sharing parameter tuples. |
| **Per-Particle** | `DDMEndpointCache` | Wiener numerical endpoint CDFs | Cleared once per particle (`values.clear()`), retaining map bucket allocations. |
| **Per-Particle** | `RaceEndpointGroupCache` | Scalar GSL race integration across shared parameter tuples | Sized to unique tuples; re-evaluated per particle. |

---

## 7. DDM & Boundary-Crossing Numerical Hardening

The Wiener diffusion likelihood (`src/ddm_functions_inline.h`, `src/model_DDM.h`, `src/likelihood_ddm.cpp`) enforces dedicated numerical hardening:
- Detect whether all truncated trials have `LT == 0` before calling `p_DDM_Wien`; if so, skip those CDF evaluations entirely (analytically `R_NegInf`).
- Detect whether any truncated trial has finite `UT`; if none do, skip `UT` CDF evaluations.
- Separate inlined Wiener PDF/CDF from generic routines to ensure full translation-unit inlining.
- Short-circuit boundary checks before entering numerical evaluation loops.

---

## 8. R-Side Vectorisation Standards

R loops and `apply` statements over large matrices must use vectorised replacements:

| Old Pattern | Replacement | Notes |
|---|---|---|
| `apply(mat, 2, sum)` | `colSums(mat)` | C-level, no per-column R dispatch |
| `apply(mat, 1, sum)` | `rowSums(mat)` | C-level, no per-row R dispatch |
| `apply(mat, 2, mean)` | `colMeans(mat)` | |
| `apply(mat, 1, which.max)` | `max.col(mat, ties.method="first")` | Fully vectorised |
| `apply(mat, 1, paste, collapse=" ")` | `do.call(paste, c(unname(as.data.frame(mat)), sep=" "))` | Avoids per-row R loop & matrix coercion |
| `rbind` inside a loop | Accumulate in a list, then `do.call(rbind, list)` | Avoids quadratic copy growth |
| `sapply(x, function(v) v[1])` | `vapply(x, `[[`, 1L, FUN.VALUE=numeric(1))` or direct indexing | Type-stable, no dispatch overhead |

---

## 9. Proactive Review Checklist

When writing or reviewing C++ likelihood code, verify:

- [ ] Parameter column layout is registered in `src/col_registry.h` with an enum and `ColSpec`.
- [ ] `validate_col_prefix()` is invoked during adapter configuration.
- [ ] All finite-check calls use `R_FINITE` / `emc2_isfinite`, `emc2_isinf`, `emc2_isnan` (never raw `std::isfinite`).
- [ ] Normal/lognormal distribution calls use `pnorm_std` / `dnorm_std` / `plnorm_std` / `dlnorm_std`.
- [ ] Log-survivors use `std::log1p(-cdf)` and log-additions use `log_sum_exp()`.
- [ ] Quadrature integrals reuse templates from `src/quad_templates.h`.
- [ ] Contaminant mixing uses `mix_contaminants` / `mix_contaminants_rt` from `src/contaminant_mixture.h`.
- [ ] Erlang timer conversions use `erlang_lambda_from_mean` from `src/timer_helpers.h`.
- [ ] Batch raw kernels accept `const double* const* cols` and index parameters via `emc2col::<model>::<param>`.
- [ ] Scalar callbacks accept `const double* par` and index parameters via `emc2col::<model>::<param>`.
- [ ] Model adapter is wired up in `src/race_dispatch.cpp` via `resolve_race_model_adapter()`.
- [ ] Optional trailing columns ($j \ge \text{N\_REQ}$) are guarded by context flags before dereferencing.
- [ ] Quantities depending only on data are pre-computed in `ModelSharedState`.
- [ ] Zero Rcpp/heap allocations inside the per-particle likelihood loop.
- [ ] Pure reduction loops over trials use `#pragma omp simd reduction(+:...)`.
- [ ] No `#include <RcppArmadillo.h>` unless Armadillo types are genuinely required.
- [ ] Comments are minimalist, factual, and omit changelogs or developmental reasoning.

When writing or reviewing R code, verify:

- [ ] No `apply(mat, 1/2, sum/mean)` — use `rowSums`/`colSums`/`rowMeans`/`colMeans`.
- [ ] No `apply(mat, 1, which.max)` — use `max.col`.
- [ ] No `rbind` in a loop — accumulate to list and `do.call(rbind, ...)`.
- [ ] No `apply(df, 1, paste)` — use `do.call(paste, ...)`.
- [ ] Row-wise string operations are kept off the MCMC critical path.

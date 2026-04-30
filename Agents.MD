# EMC2 Performance Philosophy

This document captures the optimisation philosophy that guides ongoing development of EMC2. It is intended for contributors (human or AI) who are adding new code, refactoring existing code, or reviewing changes. The overarching constraint is absolute: **speed can never change results**. Vectorising a call that produces an identical result is good; changing user-facing semantics or model behaviour is not acceptable.

---

## 1. The Golden Rule

Every optimisation must be semantics-preserving. Before shipping any performance change:

- Confirm numerical output is identical (or within documented floating-point tolerance) against the previous path.
- Run the full test suite. Snapshot tests in `tests/testthat/_snaps/` are the ground truth for likelihood values — regenerate them only when a deliberate numerical change is made, and document why.
- Never alter the public API of a model or a likelihood function in the name of performance. Inputs, outputs, and parameter names seen by users must remain stable.

---

## 2. Compiler-Level Optimisation (C++ / Rcpp)

### 2.1 Build flags

`src/Makevars` (and `Makevars.win`) set the following flags, which must be kept:

```makefile
override CXXFLAGS += -O3 -march=native -ffast-math -fno-math-errno -Igsl_bundled -DUSE_FAST_PNORM -DR_NO_REMAP
```

- **`-O3`**: full inlining, loop unrolling, and auto-vectorisation. R's default `-O2` is insufficient for tight likelihood loops.
- **`-march=native`**: enables AVX2/SSE4 SIMD on the build machine. Column-major data layouts and `#pragma omp simd` loops are written to exploit this.
- **`-ffast-math`**: allows the compiler to reorder and contract floating-point operations, enabling auto-vectorisation of FP loops. This is safe here because MCMC log-likelihood differences matter; per-sample ULP precision does not.
- **`-DUSE_FAST_PNORM`**: activates the Hart rational-approximation normal CDF in `pnorm_std()` (see §3).
- **`-DR_NO_REMAP`**: required for build on some machines, ensures R plays nicely with C++.
- **`-Igsl_bundled`**: Ensures the bundled gsl functions are used to remove external dependency (Rtools installs gsl for windows users but linux users were encountering missing dependency issues)

### 2.2 `-ffast-math` incompatibility guard

`-ffast-math` enables `-ffinite-math-only`, which breaks `std::isfinite`, `std::isinf`, and `std::isnan` because the compiler assumes all values are finite. **Never use the standard versions directly.** Instead use the wrappers defined in `src/utility_functions.h`:

```cpp
emc2_isfinite(x)   // → R_FINITE(x)
emc2_isinf(x)      // → !R_FINITE(x) && !ISNAN(x)
emc2_isnan(x)      // → ISNAN(x)
```

These delegate to Rinternals macros that are unconditionally correct. Apply this rule to any new C++ code in the package.

---

## 3. Fast Mathematical Approximations

### 3.1 `pnorm_std` — Hart rational approximation

`src/wald_functions.h` provides `pnorm_std(x, lower, log_p)` as a drop-in replacement for `R::pnorm(x, 0, 1, lower, log_p)`. When `USE_FAST_PNORM` is defined it uses a 7-term rational polynomial (Hart algorithm) that is 3–8× faster than R's implementation. The approximation is accurate to machine precision for all practical argument ranges; it is used by LBA, RDM, Wald, LNR, exGaussian, and transform utilities.

The same header also provides `plnorm_std` and `dlnorm_std` as fast wrappers around the log-normal CDF/PDF, routing through `pnorm_std`.

**Policy**: all new model headers and likelihood helpers must use `pnorm_std` / `plnorm_std` / `dlnorm_std`. Direct calls to `R::pnorm`, `R::plnorm`, `R::dlnorm` are a red flag in code review.

### 3.2 `log1p` / `expm1`

Use `std::log1p` and `std::expm1` wherever the argument is near zero. In particular, log-survivor probabilities `log(1 - CDF)` must be written as `std::log1p(-cdf)` to avoid catastrophic cancellation when `cdf ≈ 1`.

### 3.3 Log-space arithmetic

Compute likelihoods in log-space end-to-end. Avoid materialising intermediate probabilities only to immediately take their log. Where possible, use functions that return the log directly (e.g. `dlnorm_std(..., true)` for the log-PDF) instead of computing the linear value and then calling `std::log`.

---

## 4. C++ Memory and Data Layout

### 4.1 Column-major parameter buffers

Rcpp `NumericMatrix` is column-major. Raw-buffer likelihood functions (see §4.2) accept a `const double* pars_cm` pointer where column `j` starts at `pars_cm + j * n_rows`. This layout gives sequential cache lines when iterating over all rows of a single parameter — the inner loop over trials is cache-hot.

New raw-buffer functions must follow this convention. Never pass a row-major copy of the parameter matrix into a hot loop.

### 4.2 Raw-buffer fast paths

The central race-model likelihood (`c_log_likelihood_race`) supports two function-pointer interfaces:

- **`RaceCdfFun` / `RacePdf1Fun`**: the original Rcpp-vector interface, used for censored and truncated trials where the per-trial logic is complex.
- **`RaceRawFun`**: a raw-buffer interface introduced for the all-finite, no-truncation fast path. Signature:

```cpp
typedef void (*RaceRawFun)(const double* rt,
                           const double* pars_cm,   // column-major
                           int n_rows,
                           const int* mask,          // 1 = compute, 0 = skip
                           const int* isok,          // parameter validity
                           double* out,              // write log-density here
                           double min_ll,
                           void* ctx);
```

Every accumulator model (LBA, RDM, LNR, Wald, ...) must implement both a `_raw` density and a `_raw` survivor variant. The fast path eliminates all per-particle Rcpp heap allocations for the common case (finite RTs, no truncation), which is the dominant computational path in MCMC.

When adding a new accumulator model, implement `_raw` variants alongside the standard CDF/PDF wrappers. Do not leave a model with only the Rcpp-vector path.

### 4.3 Pre-allocated scratch buffers via `ModelSharedState`

`ModelSharedState` (and its subtypes `RaceSharedState`, `DDMSharedState`, etc.) are built **once per data set**, outside the particle loop, and reused across all particles. They hold:

- Pre-read censoring/truncation bound vectors (`LT`, `UT`, `LC`, `UC`).
- Pre-computed finite-trial masks (as both `LogicalVector` and `std::vector<int>` for SIMD loops).
- Pre-partitioned index vectors for finite vs. non-finite unique trials.
- Pre-allocated `double` scratch buffers sized to `n_trials`.

**Policy**: Any quantity that depends only on the data (not the particle parameters) must be computed in the shared-state builder, not inside the per-particle likelihood call. Moving work outside the inner loop is the highest-yield optimisation pattern in this codebase.

### 4.4 `std::copy` / `std::fill` over manual loops

Use `std::copy` and `std::fill` (or `std::fill_n`) for bulk array operations. The compiler recognises these and emits `memcpy`/`memset` or vectorised equivalents. Manual element loops over plain arrays are acceptable only when per-element conditional logic prevents vectorisation.

### 4.5 `std::unordered_map` for lookups

Linear scans (`O(N·M)`) over parameter name vectors are replaced with `std::unordered_map<std::string, double>` lookups (`O(1)` average). This applies to parameter extraction in `ParamTable`, transform-spec construction, and any other place where a name is looked up inside a loop over trials or particles.

### 4.6 `#pragma omp simd`

Pure reduction loops (e.g. summing log-likelihoods across trials) are annotated with `#pragma omp simd reduction(+:total_ll)`. This instructs the compiler to emit SIMD reduction instructions without requiring OpenMP threading. It is safe in contexts where each iteration is independent and the only shared state is the accumulator.

Use this pragma on any new loop that:
1. Iterates over `n_trials` or `n_particles` without data-dependent branching.
2. Accumulates a scalar result.
3. Has no pointer aliasing that the compiler cannot prove away.

Do not use it on loops that call Rcpp functions or that branch on per-element data.

---

## 5. R-Side Vectorisation

R loops and `apply` statements over many rows are slow. The following replacements must be used by default:

| Old pattern | Replacement | Notes |
|---|---|---|
| `apply(mat, 2, sum)` | `colSums(mat)` | C-level, no R evaluation per column |
| `apply(mat, 1, sum)` | `rowSums(mat)` | |
| `apply(mat, 2, mean)` | `colMeans(mat)` | |
| `apply(mat, 1, which.max)` | `max.col(mat, ties.method="first")` | Fully vectorised |
| `apply(mat, 1, paste, collapse=" ")` | `do.call(paste, c(unname(as.data.frame(mat)), sep=" "))` | Avoids per-row R loop and matrix coercion |
| `rbind` inside a loop | Accumulate as a list, then `do.call(rbind, list)` | Pre-allocation pattern; avoid quadratic copy growth |
| `sapply(x, function(v) v[1])` | `vapply(x, `[[`, 1L, FUN.VALUE=numeric(1))` or direct indexing | Type-stable, no dispatch overhead |

The `do.call(rbind, list)` pattern in particular must replace any loop that grows a matrix by successive `rbind` calls. Quadratic copy cost from `rbind` in a loop is a recurring performance anti-pattern.

---

## 6. Caching Computed Quantities

### 6.1 TrendRuntime filter caching

`TrendRuntime` caches filtered `TransformSpec` vectors so that the spec-selection step is paid once at setup rather than once per particle. When adding new transform types, the filter result must be captured in the runtime object, not recomputed in the likelihood body.

### 6.2 Column-index caching

Attribute and column reads on `DataFrame`/`NumericMatrix` inside per-particle loops are replaced with a one-time index lookup stored in `ModelSharedState`. For example, the `pContaminant` column is searched once and stored as an integer index (`-1` if absent); subsequent per-particle calls read the index directly.

---

## 7. DDM-Specific Numerical Hardening

The Wiener diffusion likelihood (`src/ddm_functions_inline.h`, `src/model_DDM.h`) has dedicated hardening beyond the general rules:

- Detect whether all truncated trials have `LT == 0` before calling `p_DDM_Wien`; if so, skip those CDF evaluations entirely (the result is analytically `R_NegInf`).
- Detect whether any truncated trial has finite `UT`; if none do, skip the `UT` CDF calls.
- The inlined Wiener PDF/CDF is separated from the generic `pdf_fncs.cpp`/`cdf_fncs.cpp` path to keep the critical DDM path in a single translation unit for inlining.

New CDF-based models should apply the same boundary-condition short-circuits before entering the evaluation loop.

---

## 8. Proactive Review Checklist

When reviewing or writing any C++ likelihood code, check:

- [ ] All finite-check calls use `R_FINITE` / `emc2_isfinite`, not `std::isfinite`.
- [ ] Normal CDF/PDF calls use `pnorm_std` / `dlnorm_std`, not `R::pnorm` / `R::dlnorm`.
- [ ] Log-survivor values use `std::log1p(-cdf)`, not `std::log(1 - cdf)`.
- [ ] Parameter matrix access inside the trial loop is column-major.
- [ ] Quantities that depend only on data are moved into a `SharedState` builder.
- [ ] Scratch buffers are pre-allocated; no `Rcpp::NumericVector` constructed inside the particle loop.
- [ ] Pure reduction loops over trials have `#pragma omp simd reduction`.
- [ ] New accumulator models supply `_raw` buffer variants for the all-finite fast path.

When reviewing or writing any R code, check:

- [ ] No `apply(mat, 1/2, sum/mean)` — use `rowSums`/`colSums`/`rowMeans`/`colMeans`.
- [ ] No `apply(mat, 1, which.max)` — use `max.col`.
- [ ] No `rbind` in a loop — accumulate to list and `do.call(rbind, ...)` after.
- [ ] No `apply(df, 1, paste)` — use `do.call(paste, ...)`.
- [ ] Row-wise string keys built in a loop are a known hot-spot; verify they are not on the critical MCMC path.

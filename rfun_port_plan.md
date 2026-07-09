# C++ simulation kernels for posterior prediction (`rfun` port) — implementation sketch

Status: **design only** (plan.md commit 14). Nothing here is implemented; this
document is the blueprint for a future PR. Scope was agreed as: sketch in
detail now, implement later; validation via fast simulate-and-recover runs
(KillTests.R-style), **not** SBC.

## 1. Problem and scope

`make_data()` / `predict()` route every simulated dataset through the R
`rfun`s (`rLBA`, `rBAwL`, `rRDM`, `rRDMSWTN`, plus the SS samplers). These are
vectorized but allocate heavily (per-call `matrix(Inf, nr, n_trials)`,
`data.frame` staging, `msm::rtnorm`), and posterior prediction calls them once
per posterior draw. The port moves the four race families onto trial-vectorized
C++ kernels sharing the likelihood's column-major parameter layout:

| family  | R sampler   | notes |
|---------|-------------|-------|
| LBA / LBAIO | `rLBA` (model_LBA.R:135) | closed form; `posdrift` toggles truncated drift |
| BAwL    | `rBAwL` (model_LBA.R:347) | leak hit time analytic (same C1/C2 algebra as `pleakyba_norm`, model_LBA.h:143–200) |
| RDM (+pContaminant) | `rRDM` → `rWald` (model_RDM.R:30,101) | Michael–Schucany–Haas IG sampler |
| RDMSWTN | `rRDMSWTN` → `rSWTN` (model_RDM.R:495,872) | trial drift draw + Wald + Erlang guess/kill clocks |

Out of scope for 14a/14b: DDM (`rtdists`-quality DDM simulation is a separate
project), SS families (staircase logic is inherently sequential and already
cheap relative to fitting), LNR/EXG/GAMMA (already trivial draws), MRI.

## 2. Architecture

### 2.1 New header `src/model_rng.h`

Pure functions, no Rcpp types in the inner loops; all draws go through R's RNG
(`R::rnorm`, `R::runif`, `R::rgamma`, `R::qnorm`) so `set.seed()` reproducibility
and single-threaded semantics are preserved. `GetRNGstate/PutRNGstate` is
handled by Rcpp attributes at the export boundary. **No OpenMP** in these
kernels — R's RNG is not thread-safe and prediction is IO-bound on data-frame
assembly anyway.

Core primitives (all scalar, inlined; `_r` suffix = consumes RNG):

```cpp
// One-sided truncated normal N(mu, sd) on [lo, Inf). Uses inverse-CDF for
// modest truncation (matching rSWTN's qnorm(lo + u*(1-lo)) exactly in
// distribution) and Robert (1995) exponential rejection when
// (lo - mu)/sd > ~4 to avoid qnorm tail cancellation. msm::rtnorm equivalent.
double rtnorm_lower_r(double mu, double sd, double lo);

// Wald/IG first passage: criterion k, rate l, diffusion s.
// Exactly the R rwaldt algorithm (model_RDM.R:33-63):
//   l <= tiny  -> Levy: c / qnorm(1 - u/2)^2 with c = (k/s)^2
//   else       -> MSH transform: y = rnorm()^2; mu = k/l; lam = (k/s)^2;
//                 x0 = mu + mu^2 y/(2 lam) - mu/(2 lam) sqrt(4 mu lam y + mu^2 y^2);
//                 accept x0 if u <= mu/(mu+x0) else mu^2/x0.
double rwald_fpt_r(double k, double l, double s);

// Full rWald semantics incl. start-point draw and defective negative drift:
//   bs = B + U(0, A);
//   v > 0 (or v >= 0 when !posdrift) -> rwald_fpt_r(bs, v, s)
//   v < 0 && !posdrift -> hit w.p. exp(2 v bs / s^2); if hit,
//                         rwald_fpt_r(bs, -v, s); else Inf.
//   v <= 0 && posdrift -> Inf.
double rwald_acc_r(double B, double v, double A, double s, bool posdrift);

// Erlang clock draw shared by RDMSWTN/BAwL guess & kill timers
// (.rdmswtn_erlang_omega semantics):
//   shape 1/2: rgamma(shape, rate = lambda)
//   shape 3 (mixed): with prob omega use shape 1 rate lambda,
//                    else shape 2 rate 2*lambda.
// lambda <= 0 or NA -> +Inf (clock disabled).
double rerlang_clock_r(double lambda, int shape, double omega);
```

Per-family trial kernels (column-major `pars` in the model's canonical
Ttransformed layout — same column meanings the R rfuns index by name):

```cpp
// LBA: dt = (b - A*u) / rtnorm_lower_r(v, sv, posdrift ? 0 : -Inf);
// dt < 0 -> Inf.   pars cols: v, sv, b, A, t0.
void rlba_trials_r(const double* pars_cm, int n_rows, int n_acc,
                   const int* ok, bool posdrift,
                   double* dt_out /* n_rows, decision time incl. nothing */);

// BAwL: draw drift ~ rtnorm_lower_r(v, sv, posdrift?0:-Inf), start u*A;
//   k < eps           -> LBA formula
//   drift > k*b       -> dt = -(1/k) log((drift - k b)/(drift - k A u)),
//                        ratio clamped to [double_xmin, 1-1e-15]
//   else              -> Inf (never reaches leaky asymptote)
// Matches rBAwL and inverts x(t) = b under the same E/G leak terms used by
// pleakyba_norm (model_LBA.h:143-200): x(t) = (v/k)(1-e^{-kt}) + a e^{-kt}.
void rbawl_trials_r(...);

// RDM: dt = rwald_acc_r(B, v, A, s, posdrift); pars cols v, B, A, t0 (+ s).
void rrdm_trials_r(...);

// RDMSWTN: v_draw = (sv > 1e-12) ? (posdrift ? rtnorm_lower_r(v, sv, 0)
//                                            : R::rnorm(v, sv)) : v;
//          dt = rwald_acc_r(b - A, v_draw, A, s, posdrift)   [rSWTN semantics]
void rrdmswtn_trials_r(...);
```

Race resolution + clocks, shared epilogue (mirrors the tail of `rRDMSWTN` /
`rBAwL` exactly — this is where the semantics live, keep it one function):

```cpp
struct RaceSimResult { std::vector<int> R;      // 1-based winner code, 0 = NA
                       std::vector<double> rt;  // NA_REAL when R == 0
                       std::vector<int> isTime; // only filled when "time" level present
                     };

RaceSimResult resolve_race_r(
    const double* dt,        // n_acc x n_trials finish times, raw axis (t0 added)
    const double* t0_col, int n_acc, int n_trials,
    const int* lR_codes,     // 1..n_acc per row
    int time_code, int nogo_code,           // -1 when absent
    const ErlangSpec* guess, const ErlangSpec* kill, bool global_kill, ...);
```

Ordering constraints that must be preserved from R (all determined by reading
`rRDMSWTN` model_RDM.R:872–980 and `rBAwL` model_LBA.R:347–470):

1. **Raw-time axis**: `dt` gets `t0` added *before* Erlang clocks compete
   (kill/guess clocks run from stimulus onset — same convention the
   likelihood uses since 564af55a).
2. **global_kill**: one timer per trial, `lambda_k` must be constant across
   accumulators (validate, error otherwise); trial killed ⇒ all
   accumulators' finishes after `tk` become losers; if nothing finishes
   before `tk` the trial is an omission (R = NA, rt = NA) unless a guess
   clock fires.
3. **local kill**: per-accumulator timer; `tk <= dt` ⇒ that accumulator's
   finish → Inf.
4. **guess**: per-accumulator guess timer on non-`nogo` accumulators; the
   effective finish is `min(dt, tg)`; if the winner is a guess-timer win, R
   is that accumulator but flagged (matches R where guess time replaces dt).
5. **Winner**: argmin over effective finishes; all-Inf column ⇒ R = NA.
6. **`time` accumulator**: if winner is the `time` level, resample R
   uniformly from non-time/non-nogo levels (`.apply_timed_guess_winner`,
   utils.R:53–72) and set `isTime = TRUE`. Do this **inside** the kernel so
   the RACE wrapper does not need a second pass.
7. **ok mask**: rows with `!ok[row]` never finish (dt = Inf) and trials whose
   first-accumulator row is !ok return R = NA (see `rLBA`'s
   `ok <- matrix(ok, nrow)[1,]` reduction).

### 2.2 Exports (RcppExports)

One entry point per family, thin Rcpp wrappers:

```cpp
// [[Rcpp::export]]
Rcpp::List rlba_cpp(Rcpp::NumericMatrix pars, Rcpp::IntegerVector lR_codes,
                    Rcpp::CharacterVector lR_levels, Rcpp::LogicalVector ok,
                    bool posdrift);
// rrdm_cpp, rrdmswtn_cpp(..., int erlang_shape, std::string erlang_type),
// rbawl_cpp(..., bool guess, bool global)
```

Returning `list(R = integer 1-based, rt = numeric, isTime = logical|NULL)`;
the R wrapper builds the factor + data.frame (cheap, once).

### 2.3 R dispatch (14b)

In each model file the `rfun` closure grows a C++ fast path, keeping the R
implementation as the reference fallback (same pattern as the likelihood
ports, toggleable for tests):

```r
rfun = function(data, pars) {
  if (isTRUE(getOption("EMC2.cpp_rfun", TRUE)) &&
      is.null(attr(pars, "staircase")))
    return(.rfun_cpp_race("LBA", data$lR, pars, attr(pars, "ok"), posdrift))
  rLBA(data$lR, pars, ok = attr(pars, "ok"), posdrift = posdrift)
}
```

- `RACE_rfun` (make_data.R:574) and `LogicalRules_rfun` keep their R loop
  over RACE subsets / rule accumulators but each subset call lands on the
  C++ kernel — no change needed there beyond the closure above.
- Staircase SS simulation is untouched (R loop stays).
- `pars` arrive **already Ttransformed and bound-checked** (make_data.R:497–
  511), so kernels never re-check bounds; they only honor `ok`.

## 3. RNG parity and testing consequences

The C++ kernels consume R's RNG stream in a *different order* than the R
rfuns (e.g. `msm::rtnorm` uses rejection internally; `rwaldt` draws vectors
`y`, `z` per subset rather than interleaved per trial). Therefore:

- **Seeded simulated datasets will differ** from the R implementation at the
  same seed. Every test that `set.seed(); make_data()` and asserts on exact
  values (snapshot tests, `_snaps/`) must either pin `options(EMC2.cpp_rfun
  = FALSE)` in setup or have snapshots regenerated once, in the same commit,
  with the flag's default flipped. Audit: `grep -l make_data tests/testthat`.
- Distributional equivalence is the correctness criterion, not stream
  equality.

## 4. Validation plan (no SBC)

1. **Kernel-level distribution tests** (new `test-rfun-cpp.R`):
   for each family × edge case (posdrift on/off, A = 0, sv = 0, negative v
   with posdrift=FALSE, k→0 for BAwL, each erlang_type for RDMSWTN), draw
   n = 2e4 trials from R rfun and C++ rfun at fixed params;
   compare (a) response proportions within binomial 99.9% CI, (b) two-sample
   KS statistic on winner RTs < 0.02, (c) omission rates. Fast (<2 s/family).
2. **Likelihood cross-check**: empirical CDF of C++-simulated data against
   the model's own `pfun` at the generating params (this catches parameter-
   column mix-ups that two matched-but-wrong samplers would miss).
3. **Simulate-and-recover** (KillTests.R / demo_cens_trunc.R style, ~20 s
   fits, single condition): RDMSWTN local_kill + global_kill, RDM +
   pContaminant, LBA censored+truncated, BAwL local_kill. Recovery criterion:
   true values inside 95% CIs for ≥ 90% of parameters, posterior means within
   0.15 of truth on the sampled scale — same informal bar the existing
   WorkingTests scripts use.
4. **Performance gate**: `predict()` on a fitted RDMSWTN emc object
   (100 posterior draws × 1000 trials): expect ≥ 5x wall-clock reduction;
   record numbers in the commit message.

## 5. Commit slicing

- **14a** `src/model_rng.h` + exports + `test-rfun-cpp.R` (kernels validated
  against R rfuns, dispatch still off by default).
- **14b** flip `EMC2.cpp_rfun` default to TRUE, wire the four model files'
  `rfun` closures, regenerate affected snapshots, run the simulate-and-
  recover suite, benchmark `predict()`.

## 6. Risks / open questions

- `rtnorm` tail behaviour: `qnorm`-inversion loses precision for
  `pnorm(0, v, sv)` ≈ 1 (strongly negative v with posdrift). Robert's
  rejection sampler avoids this; verify with KS tests at v/sv = -6.
- BAwL `guess=TRUE` + `erlang=3` mixed-omega interaction has no existing
  test coverage in R either — add the R-vs-C++ comparison case regardless.
- `-ffast-math`: kernels must use `emc2_isfinite`/`R_FINITE` (never
  `std::isfinite`) and avoid NaN-based control flow, consistent with the
  likelihood code.
- `rgamma(shape, rate)` in C++ is `R::rgamma(shape, scale)` — **scale, not
  rate**; invert explicitly (classic porting bug).

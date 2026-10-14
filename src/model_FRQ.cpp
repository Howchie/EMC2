// R-callable entry points for FRQ.
//
// These live in their own translation unit rather than in model_FRQ.h because
// the header is also included by model_rng.cpp for the simulator, and a header
// carrying [[Rcpp::export]] definitions can only be seen by one TU.  They
// bypass ContextForRaceModels entirely and so take the parameters directly.

#include "utility_functions.h"
#include "model_FRQ.h"
#include "race_contract.h"

// ---------------------------------------------------------------------------
// R-callable entry points.  All take the DECISION time (rt - t0); R/model_FRQ.R
// subtracts t0 before calling, exactly as it does for the other race models.
// Arguments recycle against `t`, as in dbawd/pbawd: each must be length 1 or
// the full length, and anything else is an error (see frq_check_len).
// ---------------------------------------------------------------------------

// Only scalars and full-length vectors recycle.  A partial-length argument is
// an error rather than R-style wrap-around: the loops index directly, so
// silently accepting one would read past the end of the vector.
static inline void frq_check_len(const NumericVector& x, int n, const char* nm) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("FRQ: `%s` must be length 1 or length %d, not %d.",
               nm, n, static_cast<int>(x.size()));
}

static inline double frq_pick(const NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

static inline FrqCache* frq_cache_from_context(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr && ctx->frq_cache
    ? ctx->frq_cache.get() : nullptr;
}

// [[Rcpp::export]]
NumericVector dfrq(NumericVector t, NumericVector K, NumericVector v,
                   NumericVector sv, NumericVector h,
                   NumericVector delta = NumericVector::create(0.0),
                   bool log_out = false,
                   NumericVector clock = NumericVector::create(1.0)) {
  const int n = t.size();
  frq_check_len(K, n, "K");         frq_check_len(v, n, "v");
  frq_check_len(sv, n, "sv");       frq_check_len(h, n, "h");
  frq_check_len(delta, n, "delta");
  frq_check_len(clock, n, "clock");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(K, i), pick(v, i), pick(sv, i),
                                pick(h, i), pick(delta, i), pick(clock, i));
    const double lp = frq_log_pdf_dt(t[i], s);
    out[i] = log_out ? lp : ((lp > R_NegInf) ? std::exp(lp) : 0.0);
  }
  return out;
}

// `lower_tail = false` returns the SURVIVOR 1 - F, which tends to the
// defective mass 1 - h rather than to zero.  It is a separate branch rather
// than 1 - pfrq(...) because the whole point of the reflection identity is to
// avoid that subtraction.
// [[Rcpp::export]]
NumericVector pfrq(NumericVector t, NumericVector K, NumericVector v,
                   NumericVector sv, NumericVector h,
                   NumericVector delta = NumericVector::create(0.0),
                   bool lower_tail = true, bool log_out = false,
                   NumericVector clock = NumericVector::create(1.0)) {
  const int n = t.size();
  frq_check_len(K, n, "K");         frq_check_len(v, n, "v");
  frq_check_len(sv, n, "sv");       frq_check_len(h, n, "h");
  frq_check_len(delta, n, "delta");
  frq_check_len(clock, n, "clock");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(K, i), pick(v, i), pick(sv, i),
                                pick(h, i), pick(delta, i), pick(clock, i));
    double lp;
    if (lower_tail) {
      lp = frq_log_cdf_dt(t[i], s);
    } else {
      lp = frq_log_surv_dt(t[i], s);
      if (ISNAN(lp)) lp = R_NaN;
    }
    if (ISNAN(lp)) { out[i] = R_NaN; continue; }
    out[i] = log_out ? lp : ((lp > R_NegInf) ? std::exp(lp) : 0.0);
  }
  return out;
}

// Map the fitted coordinates to the generative ones: the Beta shape
// beta = (v/sv)^2, the availability p, the registration rate lambda and the
// expected total evidence E (Inf when h = 1), and log(1 - p).  Used by
// Ttransform and the R reference simulator.
// [[Rcpp::export]]
NumericMatrix frq_rate(NumericVector K, NumericVector v, NumericVector sv,
                       NumericVector h,
                       NumericVector delta = NumericVector::create(0.0),
                       NumericVector clock = NumericVector::create(1.0)) {
  // Length is the longest argument, not K's: a scalar K against a vector v has
  // to give one row per v.
  const int n = static_cast<int>(std::max(std::max(std::max(K.size(), v.size()),
                                                   std::max(sv.size(), h.size())),
                                          std::max(delta.size(), clock.size())));
  frq_check_len(K, n, "K");         frq_check_len(v, n, "v");
  frq_check_len(sv, n, "sv");       frq_check_len(h, n, "h");
  frq_check_len(delta, n, "delta");
  frq_check_len(clock, n, "clock");
  NumericMatrix out(n, 5);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(K, i), pick(v, i), pick(sv, i),
                                pick(h, i), pick(delta, i), pick(clock, i));
    out(i, 0) = s.ok ? s.beta : NA_REAL;
    out(i, 1) = s.ok ? s.p : NA_REAL;
    out(i, 2) = s.ok ? ((s.poisson && s.h == 1.0) ? NA_REAL : s.lambda) : NA_REAL;
    out(i, 3) = s.ok ? s.E : NA_REAL;
    out(i, 4) = s.ok ? s.log_om_p : NA_REAL;
  }
  colnames(out) = CharacterVector::create("beta", "p", "lambda", "Etot",
                                          "log1mp");
  return out;
}

// H^{-1} for the threshold-variability generator. The R simulator pulls
// uniforms back through the SAME map the likelihood uses, and a second
// implementation in R could drift from it.
// Recycles `delta` against `y`.
// [[Rcpp::export]]
NumericVector frq_h_inv_r(NumericVector y, NumericVector delta) {
  const int n = y.size();
  frq_check_len(delta, n, "delta");
  NumericVector out(n);
  double last = R_NaN;
  FrqH H;
  for (int i = 0; i < n; ++i) {
    const double d = frq_pick(delta, i);
    if (d != last) { H = frq_h_make(d); last = d; }
    out[i] = H.ok ? frq_h_inv(y[i], H) : NA_REAL;
  }
  return out;
}

// ============================================================
// FRQ (Finite Reservoir Quorum) adapters
// Column layout: K=0, v=1, sv=2, h=3, t0=4, delta=5, clock=6.
// delta is the threshold-variability half-width (zero by default) and clock the
// registration clock shape (one by default, a constant-rate counter).
//
// The upper tail may be defective: an accumulator terminates with probability
// h = I_p(alpha, beta), and the leftover mass 1 - h sits at t = +Inf.  At
// h = 1 the survivor reaches zero and the same paths become proper.
//
// Every entry point resolves its row through FrqMemo, which splits the state
// on (K, sv/v, h, delta) -- the beta-quantile half, cached per particle --
// and (v, clock), which is recomputed per row for the cost of one log.
// ============================================================

// The scalar entry points are called once per evaluation point, often several
// times on one row (bin edges and midpoint), so they share a per-thread memo.
// Its keys are the parameter values themselves, so a hit is valid across
// particles; only a miss consults the per-particle cache.
static inline const FrqPars& frq_scalar_row(const double* par, void* ctx_) {
  static thread_local FrqMemo memo;
  memo.shared = frq_cache_from_context(ctx_);
  return memo.get(par[emc2col::frq::K], par[emc2col::frq::v],
                  par[emc2col::frq::sv], par[emc2col::frq::h],
                  par[emc2col::frq::delta], par[emc2col::frq::clock]);
}

double dfrq_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq::K]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  return frq_pdf_natural_dt(tt, frq_scalar_row(par, ctx_));
}
double pfrq_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq::K]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf deliberately reaches the kernel: the CDF there is h, not one
  // when h < 1, and is one at the proper h = 1 boundary.
  return frq_cdf_natural_dt(tt, frq_scalar_row(par, ctx_));
}

void dfrq_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* k_  = cols[emc2col::frq::K];
  const double* v_  = cols[emc2col::frq::v];
  const double* sv_ = cols[emc2col::frq::sv];
  const double* h_  = cols[emc2col::frq::h];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  const double* cu_ = cols[emc2col::frq::clock];
  FrqMemo memo(frq_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(k_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const FrqPars& s = memo.get(k_[i], v_[i], sv_[i], h_[i], de_[i], cu_[i]);
    const double log_pdf = frq_log_pdf_dt(tt, s);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

void pfrq_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* k_  = cols[emc2col::frq::K];
  const double* v_  = cols[emc2col::frq::v];
  const double* sv_ = cols[emc2col::frq::sv];
  const double* h_  = cols[emc2col::frq::h];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  const double* cu_ = cols[emc2col::frq::clock];
  FrqMemo memo(frq_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one, matching
    // every other race adapter: the trial is failed by its winner's density.
    if (R_IsNA(k_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = 0.0; continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const FrqPars& s = memo.get(k_[i], v_[i], sv_[i], h_[i], de_[i], cu_[i]);
    const double log_surv = frq_log_surv_dt(tt, s);
    if (ISNAN(log_surv)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (log_surv > R_NegInf) ? log_surv
                                   : raw_log_zero(min_ll, floor_raw);
  }
}

void frq_logS_at_t(double t, const double* const* cols,
                          int /*n_rows_total*/, int n_lR, int /*n_par*/,
                          const int* trunc_mask, int n_unique_trials,
                          const int* isok_all, void* ctx_, double* logS_out) {
  const double* k_  = cols[emc2col::frq::K];
  const double* v_  = cols[emc2col::frq::v];
  const double* sv_ = cols[emc2col::frq::sv];
  const double* h_  = cols[emc2col::frq::h];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  const double* cu_ = cols[emc2col::frq::clock];
  FrqMemo memo(frq_cache_from_context(ctx_));
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(k_[r]) || !R_FINITE(t0_[r])) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      const FrqPars& s = memo.get(k_[r], v_[r], sv_[r], h_[r], de_[r], cu_[r]);
      const double log_surv = !R_FINITE(t)
        ? (s.ok ? frq_log_surv_dt(R_PosInf, s) : R_NaN)
        : frq_log_surv_dt(tt, s);
      if (ISNAN(log_surv)) { bad = true; break; }
      logS += log_surv;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

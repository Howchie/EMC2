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

// [[Rcpp::export]]
NumericVector dfrq(NumericVector t, NumericVector alpha, NumericVector beta,
                   NumericVector h, NumericVector tau,
                   NumericVector delta = NumericVector::create(0.0),
                   bool log_out = false) {
  const int n = t.size();
  frq_check_len(alpha, n, "alpha"); frq_check_len(beta, n, "beta");
  frq_check_len(h, n, "h");         frq_check_len(tau, n, "tau");
  frq_check_len(delta, n, "delta");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(alpha, i), pick(beta, i), pick(h, i),
                                pick(tau, i), pick(delta, i));
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
NumericVector pfrq(NumericVector t, NumericVector alpha, NumericVector beta,
                   NumericVector h, NumericVector tau,
                   NumericVector delta = NumericVector::create(0.0),
                   bool lower_tail = true, bool log_out = false) {
  const int n = t.size();
  frq_check_len(alpha, n, "alpha"); frq_check_len(beta, n, "beta");
  frq_check_len(h, n, "h");         frq_check_len(tau, n, "tau");
  frq_check_len(delta, n, "delta");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(alpha, i), pick(beta, i), pick(h, i),
                                pick(tau, i), pick(delta, i));
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

// (alpha, beta, h, tau) -> (p, lambda).  Exported so that the R Ttransform
// reports the generative coordinates through the SAME inversion the likelihood
// uses; recomputing it in R would let the two drift apart silently.  Returns a
// two-column matrix, NA in both columns for a rejected parameter row.
// [[Rcpp::export]]
NumericMatrix frq_rate(NumericVector alpha, NumericVector beta,
                       NumericVector h, NumericVector tau,
                       NumericVector delta = NumericVector::create(0.0)) {
  // Length is the longest argument, not alpha's: a scalar alpha against a
  // vector beta has to give one row per beta.
  const int n = static_cast<int>(std::max(std::max(alpha.size(), beta.size()),
                                 std::max(std::max(h.size(), tau.size()),
                                          delta.size())));
  frq_check_len(alpha, n, "alpha"); frq_check_len(beta, n, "beta");
  frq_check_len(h, n, "h");         frq_check_len(tau, n, "tau");
  frq_check_len(delta, n, "delta");
  NumericMatrix out(n, 2);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(alpha, i), pick(beta, i), pick(h, i),
                                pick(tau, i), pick(delta, i));
    out(i, 0) = s.ok ? s.p : NA_REAL;
    out(i, 1) = s.ok ? s.lambda : NA_REAL;
  }
  colnames(out) = CharacterVector::create("p", "lambda");
  return out;
}

// The hard-coded conditional-quantile level, exposed so the tests assert
// against the compiled constant rather than a copy of it.
// [[Rcpp::export]]
double frq_quantile_level() { return FRQ_QUANTILE; }

// H^{-1} for the threshold-variability generator, exported for the same reason
// frq_rate() is: the R simulator has to pull uniforms back through the SAME map
// the likelihood uses, and a second implementation in R could drift from it.
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
// Column layout: alpha=0, beta=1, h=2, tau=3, t0=4, delta=5.  No optional
// columns and no context flags -- the model has a single variant, and the
// conditional quantile defining tau is the compile-time constant FRQ_QUANTILE.
// delta is the threshold-variability half-width; it is zero by default, and
// frq_derive() then builds an inactive generator that costs nothing.
//
// The upper tail is ALWAYS defective: an accumulator terminates only with
// probability h = I_p(alpha, beta) < 1, and the leftover mass 1 - h sits at
// t = +Inf.  pfun_raw therefore returns the log-SURVIVOR (the race kernel
// contract), which saturates at log(1 - h) rather than falling to -Inf, and
// that is exactly what makes an omission trial score sum_j log(1 - h_j).
//
// Every entry point derives (p, lambda) through frq_derive(), memoised on the
// exact parameter bits because the two qbeta inversions dominate the cost and
// consecutive compressed rows usually repeat.
// ============================================================

double dfrq_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[emc2col::frq::alpha]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  const FrqPars s = frq_derive(par[emc2col::frq::alpha], par[emc2col::frq::beta],
                               par[emc2col::frq::h], par[emc2col::frq::tau],
                               par[emc2col::frq::delta]);
  return frq_pdf_natural_dt(tt, s);
}
double pfrq_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[emc2col::frq::alpha]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf deliberately reaches the kernel: the CDF there is h, not one.
  const FrqPars s = frq_derive(par[emc2col::frq::alpha], par[emc2col::frq::beta],
                               par[emc2col::frq::h], par[emc2col::frq::tau],
                               par[emc2col::frq::delta]);
  return frq_cdf_natural_dt(tt, s);
}

void dfrq_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* al_ = cols[emc2col::frq::alpha];
  const double* be_ = cols[emc2col::frq::beta];
  const double* h_  = cols[emc2col::frq::h];
  const double* ta_ = cols[emc2col::frq::tau];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  FrqMemo memo;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(al_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const FrqPars& s = memo.get(al_[i], be_[i], h_[i], ta_[i], de_[i]);
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
  const double* al_ = cols[emc2col::frq::alpha];
  const double* be_ = cols[emc2col::frq::beta];
  const double* h_  = cols[emc2col::frq::h];
  const double* ta_ = cols[emc2col::frq::tau];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  FrqMemo memo;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one, matching
    // every other race adapter: the trial is failed by its winner's density.
    if (R_IsNA(al_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = 0.0; continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const FrqPars& s = memo.get(al_[i], be_[i], h_[i], ta_[i], de_[i]);
    const double log_surv = frq_log_surv_dt(tt, s);
    if (ISNAN(log_surv)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (log_surv > R_NegInf) ? log_surv
                                   : raw_log_zero(min_ll, floor_raw);
  }
}

void frq_logS_at_t(double t, const double* const* cols,
                          int /*n_rows_total*/, int n_lR, int /*n_par*/,
                          const int* trunc_mask, int n_unique_trials,
                          const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* al_ = cols[emc2col::frq::alpha];
  const double* be_ = cols[emc2col::frq::beta];
  const double* h_  = cols[emc2col::frq::h];
  const double* ta_ = cols[emc2col::frq::tau];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  FrqMemo memo;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(al_[r]) || !R_FINITE(t0_[r])) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      const FrqPars& s = memo.get(al_[r], be_[r], h_[r], ta_[r], de_[r]);
      const double log_surv = frq_log_surv_dt(tt, s);
      if (ISNAN(log_surv)) { bad = true; break; }
      logS += log_surv;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

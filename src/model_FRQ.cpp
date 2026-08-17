// R-callable entry points for FRQ (Math/FRQ.tex).
//
// These live in their own translation unit rather than in model_FRQ.h because
// the header is also included by model_rng.cpp for the simulator, and a header
// carrying [[Rcpp::export]] definitions can only be seen by one TU.  They
// bypass ContextForRaceModels entirely and so take the parameters directly.

#include "model_FRQ.h"

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
                   NumericVector h, NumericVector tau, bool log_out = false) {
  const int n = t.size();
  frq_check_len(alpha, n, "alpha"); frq_check_len(beta, n, "beta");
  frq_check_len(h, n, "h");         frq_check_len(tau, n, "tau");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(alpha, i), pick(beta, i), pick(h, i),
                                pick(tau, i));
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
                   NumericVector h, NumericVector tau, bool lower_tail = true,
                   bool log_out = false) {
  const int n = t.size();
  frq_check_len(alpha, n, "alpha"); frq_check_len(beta, n, "beta");
  frq_check_len(h, n, "h");         frq_check_len(tau, n, "tau");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(alpha, i), pick(beta, i), pick(h, i),
                                pick(tau, i));
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
                       NumericVector h, NumericVector tau) {
  // Length is the longest argument, not alpha's: a scalar alpha against a
  // vector beta has to give one row per beta.
  const int n = static_cast<int>(std::max(std::max(alpha.size(), beta.size()),
                                          std::max(h.size(), tau.size())));
  frq_check_len(alpha, n, "alpha"); frq_check_len(beta, n, "beta");
  frq_check_len(h, n, "h");         frq_check_len(tau, n, "tau");
  NumericMatrix out(n, 2);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(alpha, i), pick(beta, i), pick(h, i),
                                pick(tau, i));
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


// Racing subordinator adapters and R-callable kernels.  The math lives in
// model_RSUB.h, which model_rng.cpp also includes for the simulator.

#include "utility_functions.h"
#include "model_RSUB.h"
#include "race_contract.h"
#include "col_registry.h"

// ---------------------------------------------------------------------------
// R-callable entry points on the DECISION time (rt - t0).  Arguments other
// than `t` must be length 1 or the full length.
// ---------------------------------------------------------------------------

// The evidence law: dense (inverse-Gaussian) or packets (compound Poisson).
static inline double rsub_log_pdf(bool packets, double u, double v, double B,
                                  double A, double s) {
  return packets ? rsub::log_pdf<rsub::cp::Kernel>(u, v, B, A, s)
                 : rsub::log_pdf<rsub::ig::Kernel>(u, v, B, A, s);
}
static inline double rsub_log_cdf(bool packets, double u, double v, double B,
                                  double A, double s) {
  return packets ? rsub::log_cdf<rsub::cp::Kernel>(u, v, B, A, s)
                 : rsub::log_cdf<rsub::ig::Kernel>(u, v, B, A, s);
}
static inline double rsub_log_surv(bool packets, double u, double v, double B,
                                   double A, double s) {
  return packets ? rsub::log_surv<rsub::cp::Kernel>(u, v, B, A, s)
                 : rsub::log_surv<rsub::ig::Kernel>(u, v, B, A, s);
}

static inline void rsub_check_len(const NumericVector& x, int n, const char* nm) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("RSUB: `%s` must be length 1 or length %d, not %d.",
               nm, n, static_cast<int>(x.size()));
}

static inline double rsub_pick(const NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

// [[Rcpp::export]]
NumericVector drsub(NumericVector t, NumericVector v, NumericVector B,
                    NumericVector A, NumericVector s, bool packets = false,
                    bool log_out = false) {
  const int n = t.size();
  rsub_check_len(v, n, "v");
  rsub_check_len(B, n, "B");
  rsub_check_len(A, n, "A");
  rsub_check_len(s, n, "s");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double vi = rsub_pick(v, i), Bi = rsub_pick(B, i),
                 Ai = rsub_pick(A, i), si = rsub_pick(s, i);
    double lp;
    if (ISNAN(t[i]) || !rsub::pars_ok(vi, Bi, Ai, si)) lp = R_NaN;
    else lp = rsub_log_pdf(packets, t[i], vi, Bi, Ai, si);
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// `lower_tail = false` returns the survivor, computed directly rather than
// as 1 - F.
// [[Rcpp::export]]
NumericVector prsub(NumericVector t, NumericVector v, NumericVector B,
                    NumericVector A, NumericVector s, bool packets = false,
                    bool lower_tail = true, bool log_out = false) {
  const int n = t.size();
  rsub_check_len(v, n, "v");
  rsub_check_len(B, n, "B");
  rsub_check_len(A, n, "A");
  rsub_check_len(s, n, "s");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double vi = rsub_pick(v, i), Bi = rsub_pick(B, i),
                 Ai = rsub_pick(A, i), si = rsub_pick(s, i);
    double lp;
    if (ISNAN(t[i]) || !rsub::pars_ok(vi, Bi, Ai, si)) lp = R_NaN;
    else lp = lower_tail ? rsub_log_cdf(packets, t[i], vi, Bi, Ai, si)
                         : rsub_log_surv(packets, t[i], vi, Bi, Ai, si);
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// ---------------------------------------------------------------------------
// Race adapters.  Column layout (col_registry.h, emc2col::rsub): v, B, A, t0, s;
// ctx->rsub_packets selects the evidence law.
// ---------------------------------------------------------------------------

static inline bool rsub_packets(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr && ctx->rsub_packets != 0;
}

static inline bool rsub_row_ok(const double* par) {
  return !R_IsNA(par[emc2col::rsub::v]) && R_FINITE(par[emc2col::rsub::t0]) &&
         rsub::pars_ok(par[emc2col::rsub::v], par[emc2col::rsub::B],
                       par[emc2col::rsub::A], par[emc2col::rsub::s]);
}

double drsub_scalar(double t, const double* par, void* ctx_) {
  if (!rsub_row_ok(par)) return 0.0;
  const double tt = t - par[emc2col::rsub::t0];
  if (!(tt > 0.0) || !R_FINITE(tt)) return 0.0;
  const double lp = rsub_log_pdf(rsub_packets(ctx_), tt, par[emc2col::rsub::v], par[emc2col::rsub::B],
                                  par[emc2col::rsub::A], par[emc2col::rsub::s]);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

double prsub_scalar(double t, const double* par, void* ctx_) {
  if (!rsub_row_ok(par)) return 0.0;
  const double tt = t - par[emc2col::rsub::t0];
  if (!(tt > 0.0)) return 0.0;
  const double lp = rsub_log_cdf(rsub_packets(ctx_), tt, par[emc2col::rsub::v], par[emc2col::rsub::B],
                                  par[emc2col::rsub::A], par[emc2col::rsub::s]);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

void drsub_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = cols[emc2col::rsub::v];
  const double* B_  = cols[emc2col::rsub::B];
  const double* A_  = cols[emc2col::rsub::A];
  const double* t0_ = cols[emc2col::rsub::t0];
  const double* s_  = cols[emc2col::rsub::s];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double tt = rt[i] - t0_[i];
    if (!isok[i] || R_IsNA(v_[i]) || !R_FINITE(t0_[i]) ||
        !rsub::pars_ok(v_[i], B_[i], A_[i], s_[i]) ||
        !(tt > 0.0) || !R_FINITE(tt)) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    out[i] = raw_log_value(rsub_log_pdf(rsub_packets(ctx_), tt, v_[i], B_[i], A_[i], s_[i]),
                           min_ll, floor_raw);
  }
}

void prsub_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = cols[emc2col::rsub::v];
  const double* B_  = cols[emc2col::rsub::B];
  const double* A_  = cols[emc2col::rsub::A];
  const double* t0_ = cols[emc2col::rsub::t0];
  const double* s_  = cols[emc2col::rsub::s];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one; the
    // trial is failed by its winner's density.
    if (!isok[i] || R_IsNA(v_[i]) || !R_FINITE(t0_[i]) ||
        !rsub::pars_ok(v_[i], B_[i], A_[i], s_[i])) {
      out[i] = 0.0;
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0)) { out[i] = 0.0; continue; }
    const double log_surv = rsub_log_surv(rsub_packets(ctx_), tt, v_[i], B_[i], A_[i], s_[i]);
    if (ISNAN(log_surv)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (log_surv > R_NegInf) ? std::fmin(log_surv, 0.0)
                                   : raw_log_zero(min_ll, floor_raw);
  }
}

void rsub_logS_at_t(double t, const double* const* cols,
                    int /*n_rows_total*/, int n_lR, int /*n_par*/,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out) {
  const double* v_  = cols[emc2col::rsub::v];
  const double* B_  = cols[emc2col::rsub::B];
  const double* A_  = cols[emc2col::rsub::A];
  const double* t0_ = cols[emc2col::rsub::t0];
  const double* s_  = cols[emc2col::rsub::s];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(v_[r]) || !R_FINITE(t0_[r]) ||
          !rsub::pars_ok(v_[r], B_[r], A_[r], s_[r])) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;  // not started: survivor one
      const double log_surv = rsub_log_surv(rsub_packets(ctx_), tt, v_[r], B_[r], A_[r], s_[r]);
      if (ISNAN(log_surv)) { bad = true; break; }
      logS += log_surv;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

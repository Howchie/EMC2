// Opponent Poisson race adapters and R-callable kernels.  The math lives in
// model_OPR.h, which model_rng.cpp also includes for the simulator.

#include <initializer_list>
#include "utility_functions.h"
#include "model_OPR.h"
#include "race_contract.h"
#include "col_registry.h"

// ---------------------------------------------------------------------------
// R-callable entry points on the DECISION time (rt - t0) and the event rates
// lambda+ = a, lambda- = b; R/model_OPR.R maps the sampled chart to the rates.
// Every argument is length 1 or the length of the longest.
// ---------------------------------------------------------------------------

static inline void opr_check_len(const NumericVector& x, int n, const char* nm) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("OPR: `%s` must be length 1 or length %d, not %d.",
               nm, n, static_cast<int>(x.size()));
}

static inline double opr_pick(const NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

static inline int opr_len(std::initializer_list<int> sizes) {
  int n = 0;
  for (int k : sizes) if (k == 0) return 0; else if (k > n) n = k;
  return n;
}

// [[Rcpp::export]]
NumericVector dopr(NumericVector t, NumericVector K, NumericVector lambda_p,
                   NumericVector lambda_m, bool log_out = false) {
  const int n = opr_len({static_cast<int>(t.size()), static_cast<int>(K.size()), static_cast<int>(lambda_p.size()), static_cast<int>(lambda_m.size())});
  opr_check_len(t, n, "t");
  opr_check_len(K, n, "K");
  opr_check_len(lambda_p, n, "lambda_p");
  opr_check_len(lambda_m, n, "lambda_m");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = opr_pick(t, i);
    const double k = opr_pick(K, i), a = opr_pick(lambda_p, i),
                 b = opr_pick(lambda_m, i);
    double lp;
    if (ISNAN(ti) || ISNAN(k) || ISNAN(a) || ISNAN(b)) lp = R_NaN;
    else lp = opr::pars_ok(k, a, b) ? opr::log_pdf(ti, k, a, b) : R_NaN;
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// `lower_tail = false` returns the survivor, which tends to the never-finish
// mass 1 - h; it is computed directly rather than as 1 - F.
// [[Rcpp::export]]
NumericVector popr(NumericVector t, NumericVector K, NumericVector lambda_p,
                   NumericVector lambda_m, bool lower_tail = true,
                   bool log_out = false) {
  const int n = opr_len({static_cast<int>(t.size()), static_cast<int>(K.size()), static_cast<int>(lambda_p.size()), static_cast<int>(lambda_m.size())});
  opr_check_len(t, n, "t");
  opr_check_len(K, n, "K");
  opr_check_len(lambda_p, n, "lambda_p");
  opr_check_len(lambda_m, n, "lambda_m");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = opr_pick(t, i);
    const double k = opr_pick(K, i), a = opr_pick(lambda_p, i),
                 b = opr_pick(lambda_m, i);
    double lp;
    if (ISNAN(ti) || !opr::pars_ok(k, a, b)) lp = R_NaN;
    else lp = lower_tail ? opr::log_cdf(ti, k, a, b)
                         : opr::log_surv(ti, k, a, b);
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// The Gamma rate mixture on the ratio chart: net rate v (the mean of the
// Gamma), opposing ratio rho and rate CV cv > 0.
// [[Rcpp::export]]
NumericVector dopr_gamma(NumericVector t, NumericVector K, NumericVector v,
                         NumericVector rho, NumericVector cv,
                         bool log_out = false) {
  const int n = opr_len({static_cast<int>(t.size()), static_cast<int>(K.size()), static_cast<int>(v.size()), static_cast<int>(rho.size()), static_cast<int>(cv.size())});
  opr_check_len(t, n, "t");
  opr_check_len(K, n, "K");
  opr_check_len(v, n, "v");
  opr_check_len(rho, n, "rho");
  opr_check_len(cv, n, "cv");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = opr_pick(t, i);
    const opr::gm::GState g = opr::gm::make_gstate(opr_pick(K, i), opr_pick(v, i),
                                                   opr_pick(rho, i), opr_pick(cv, i));
    const double lp = (ISNAN(ti) || !g.ok) ? R_NaN : opr::gm::log_pdf(ti, g);
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector popr_gamma(NumericVector t, NumericVector K, NumericVector v,
                         NumericVector rho, NumericVector cv,
                         bool lower_tail = true, bool log_out = false) {
  const int n = opr_len({static_cast<int>(t.size()), static_cast<int>(K.size()), static_cast<int>(v.size()), static_cast<int>(rho.size()), static_cast<int>(cv.size())});
  opr_check_len(t, n, "t");
  opr_check_len(K, n, "K");
  opr_check_len(v, n, "v");
  opr_check_len(rho, n, "rho");
  opr_check_len(cv, n, "cv");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = opr_pick(t, i);
    const opr::gm::GState g = opr::gm::make_gstate(opr_pick(K, i), opr_pick(v, i),
                                                   opr_pick(rho, i), opr_pick(cv, i));
    double lp;
    if (ISNAN(ti) || !g.ok) lp = R_NaN;
    else lp = lower_tail ? opr::gm::log_cdf(ti, g) : opr::gm::log_surv(ti, g);
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// ---------------------------------------------------------------------------
// Race adapters.  Column layout (col_registry.h, emc2col::opr): K, the two
// rate coordinates, t0, and the rate CV under the ratio chart.  ctx->opr_chart selects (v, rho), (lambda_p,
// lambda_m) or (v, lambda0).  Every accumulator with lambda- > lambda+ is defective.
// ---------------------------------------------------------------------------

static inline int opr_chart_of(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr ? ctx->opr_chart : 0;
}

// The rate CV column exists under the ratio chart only; elsewhere it is zero.
static inline const double* opr_cv_col(const double* const* cols, int chart) {
  return chart == 0 ? cols[emc2col::opr::cv] : nullptr;
}

// State for one row; null for an invalid row.
static inline const opr::Acc* opr_state(opr::StateMemo& memo, double K,
                                        double c1, double c2, double t0,
                                        double cv) {
  if (R_IsNA(K) || R_IsNA(c1) || R_IsNA(c2) || !R_FINITE(t0) || ISNAN(cv))
    return nullptr;
  const opr::Acc& s = memo.get(K, c1, c2, cv);
  return s.ok() ? &s : nullptr;
}

double dopr_scalar(double t, const double* par, void* ctx_) {
  const int chart = opr_chart_of(ctx_);
  opr::StateMemo memo(chart);
  const double t0 = par[emc2col::opr::t0];
  const double cv = chart == 0 ? par[emc2col::opr::cv] : 0.0;
  const opr::Acc* s = opr_state(memo, par[emc2col::opr::K],
                                par[emc2col::opr::r1], par[emc2col::opr::r2], t0, cv);
  if (s == nullptr) return 0.0;
  const double tt = t - t0;
  if (!(tt > 0.0)) return 0.0;
  const double lp = opr::log_pdf(tt, *s);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

double popr_scalar(double t, const double* par, void* ctx_) {
  const int chart = opr_chart_of(ctx_);
  opr::StateMemo memo(chart);
  const double t0 = par[emc2col::opr::t0];
  const double cv = chart == 0 ? par[emc2col::opr::cv] : 0.0;
  const opr::Acc* s = opr_state(memo, par[emc2col::opr::K],
                                par[emc2col::opr::r1], par[emc2col::opr::r2], t0, cv);
  if (s == nullptr) return 0.0;
  const double tt = t - t0;
  if (!(tt > 0.0)) return 0.0;
  // tt == Inf reaches the kernel: the CDF there is h.
  const double lp = opr::log_cdf(tt, *s);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

void dopr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int chart = opr_chart_of(ctx_);
  opr::StateMemo memo(chart);
  const double* cv_ = opr_cv_col(cols, chart);
  const double* k_  = cols[emc2col::opr::K];
  const double* c1_ = cols[emc2col::opr::r1];
  const double* c2_ = cols[emc2col::opr::r2];
  const double* t0_ = cols[emc2col::opr::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const opr::Acc* s =
        isok[i] ? opr_state(memo, k_[i], c1_[i], c2_[i], t0_[i],
                            cv_ ? cv_[i] : 0.0) : nullptr;
    if (s == nullptr) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0) || !R_FINITE(tt)) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    out[i] = raw_log_value(opr::log_pdf(tt, *s), min_ll, floor_raw);
  }
}

void popr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int chart = opr_chart_of(ctx_);
  opr::StateMemo memo(chart);
  const double* cv_ = opr_cv_col(cols, chart);
  const double* k_  = cols[emc2col::opr::K];
  const double* c1_ = cols[emc2col::opr::r1];
  const double* c2_ = cols[emc2col::opr::r2];
  const double* t0_ = cols[emc2col::opr::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one; the
    // trial is failed by its winner's density.
    const opr::Acc* s =
        isok[i] ? opr_state(memo, k_[i], c1_[i], c2_[i], t0_[i],
                            cv_ ? cv_[i] : 0.0) : nullptr;
    if (s == nullptr) {
      out[i] = 0.0;
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0)) { out[i] = 0.0; continue; }
    const double log_surv = opr::log_surv(tt, *s);
    if (ISNAN(log_surv)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (log_surv > R_NegInf) ? std::fmin(log_surv, 0.0)
                                   : raw_log_zero(min_ll, floor_raw);
  }
}

// At t = +Inf this returns log S(+Inf), the never-finish atom, which the
// batched defective truncation normaliser adds to the finite-window mass.
void opr_logS_at_t(double t, const double* const* cols,
                   int /*n_rows_total*/, int n_lR, int /*n_par*/,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* ctx_, double* logS_out) {
  const int chart = opr_chart_of(ctx_);
  opr::StateMemo memo(chart);
  const double* cv_ = opr_cv_col(cols, chart);
  const double* k_  = cols[emc2col::opr::K];
  const double* c1_ = cols[emc2col::opr::r1];
  const double* c2_ = cols[emc2col::opr::r2];
  const double* t0_ = cols[emc2col::opr::t0];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = start + kk;
      const opr::Acc* s =
          isok_all[r] ? opr_state(memo, k_[r], c1_[r], c2_[r], t0_[r],
                                  cv_ ? cv_[r] : 0.0) : nullptr;
      if (s == nullptr) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;  // not started: survivor one
      const double log_surv = opr::log_surv(tt, *s);
      if (ISNAN(log_surv)) { bad = true; break; }
      logS += log_surv;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

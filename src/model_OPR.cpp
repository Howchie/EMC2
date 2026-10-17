// Opponent Poisson race adapters and R-callable kernels.  The math lives in
// model_OPR.h, which model_rng.cpp also includes for the simulator.

#include "utility_functions.h"
#include "model_OPR.h"
#include "race_contract.h"
#include "col_registry.h"

// ---------------------------------------------------------------------------
// R-callable entry points on the DECISION time (rt - t0) and the event rates
// lambda+ = a, lambda- = b; R/model_OPR.R maps the sampled chart to the rates.
// Arguments other than `t` must be length 1 or the full length.
// ---------------------------------------------------------------------------

static inline void opr_check_len(const NumericVector& x, int n, const char* nm) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("OPR: `%s` must be length 1 or length %d, not %d.",
               nm, n, static_cast<int>(x.size()));
}

static inline double opr_pick(const NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

// [[Rcpp::export]]
NumericVector dopr(NumericVector t, NumericVector K, NumericVector lambda_p,
                   NumericVector lambda_m, bool log_out = false) {
  const int n = t.size();
  opr_check_len(K, n, "K");
  opr_check_len(lambda_p, n, "lambda_p");
  opr_check_len(lambda_m, n, "lambda_m");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double k = opr_pick(K, i), a = opr_pick(lambda_p, i),
                 b = opr_pick(lambda_m, i);
    double lp;
    if (ISNAN(t[i]) || ISNAN(k) || ISNAN(a) || ISNAN(b)) lp = R_NaN;
    else lp = opr::pars_ok(k, a, b) ? opr::log_pdf(t[i], k, a, b) : R_NaN;
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
  const int n = t.size();
  opr_check_len(K, n, "K");
  opr_check_len(lambda_p, n, "lambda_p");
  opr_check_len(lambda_m, n, "lambda_m");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double k = opr_pick(K, i), a = opr_pick(lambda_p, i),
                 b = opr_pick(lambda_m, i);
    double lp;
    if (ISNAN(t[i]) || !opr::pars_ok(k, a, b)) lp = R_NaN;
    else lp = lower_tail ? opr::log_cdf(t[i], k, a, b)
                         : opr::log_surv(t[i], k, a, b);
    out[i] = log_out ? lp : (ISNAN(lp) ? R_NaN : std::exp(lp));
  }
  return out;
}

// ---------------------------------------------------------------------------
// Race adapters.  Column layout (col_registry.h, emc2col::opr): K, the two
// rate coordinates, t0.  ctx->opr_chart selects (v, rho), (lambda_p,
// lambda_m) or (v, lambda0).  Every accumulator with lambda- > lambda+ is defective.
// ---------------------------------------------------------------------------

static inline int opr_chart_of(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr ? ctx->opr_chart : 0;
}

// State for one row; null for an invalid row.
static inline const opr::State* opr_state(opr::StateMemo& memo, double K,
                                          double c1, double c2, double t0) {
  if (R_IsNA(K) || R_IsNA(c1) || R_IsNA(c2) || !R_FINITE(t0)) return nullptr;
  const opr::State& s = memo.get(K, c1, c2);
  return s.ok ? &s : nullptr;
}

double dopr_scalar(double t, const double* par, void* ctx_) {
  opr::StateMemo memo(opr_chart_of(ctx_));
  const double t0 = par[emc2col::opr::t0];
  const opr::State* s = opr_state(memo, par[emc2col::opr::K],
                                  par[emc2col::opr::r1], par[emc2col::opr::r2], t0);
  if (s == nullptr) return 0.0;
  const double tt = t - t0;
  if (!(tt > 0.0)) return 0.0;
  const double lp = opr::log_pdf(tt, *s);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

double popr_scalar(double t, const double* par, void* ctx_) {
  opr::StateMemo memo(opr_chart_of(ctx_));
  const double t0 = par[emc2col::opr::t0];
  const opr::State* s = opr_state(memo, par[emc2col::opr::K],
                                  par[emc2col::opr::r1], par[emc2col::opr::r2], t0);
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
  opr::StateMemo memo(opr_chart_of(ctx_));
  const double* k_  = cols[emc2col::opr::K];
  const double* c1_ = cols[emc2col::opr::r1];
  const double* c2_ = cols[emc2col::opr::r2];
  const double* t0_ = cols[emc2col::opr::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const opr::State* s =
        isok[i] ? opr_state(memo, k_[i], c1_[i], c2_[i], t0_[i]) : nullptr;
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
  opr::StateMemo memo(opr_chart_of(ctx_));
  const double* k_  = cols[emc2col::opr::K];
  const double* c1_ = cols[emc2col::opr::r1];
  const double* c2_ = cols[emc2col::opr::r2];
  const double* t0_ = cols[emc2col::opr::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one; the
    // trial is failed by its winner's density.
    const opr::State* s =
        isok[i] ? opr_state(memo, k_[i], c1_[i], c2_[i], t0_[i]) : nullptr;
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
  opr::StateMemo memo(opr_chart_of(ctx_));
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
      const opr::State* s =
          isok_all[r] ? opr_state(memo, k_[r], c1_[r], c2_[r], t0_[r]) : nullptr;
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

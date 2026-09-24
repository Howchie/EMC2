#include "model_FRQfade.h"
#include "race_contract.h"
#include "utility_functions.h"

#include <algorithm>

namespace {

static inline void frqfade_check_len(const Rcpp::NumericVector& x, int n,
                                     const char* name) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("FRQfade: `%s` must be length 1 or length %d, not %d.",
               name, n, static_cast<int>(x.size()));
}

static inline double frqfade_pick(const Rcpp::NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

static inline FrqFadeCache* frqfade_cache_from_context(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr && ctx->frq_fade_cache
    ? ctx->frq_fade_cache.get() : nullptr;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector dfrqfade(Rcpp::NumericVector t,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector lambda,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector delta,
                             Rcpp::NumericVector cv_u) {
  const int n = t.size();
  frqfade_check_len(alpha, n, "alpha"); frqfade_check_len(beta, n, "beta");
  frqfade_check_len(lambda, n, "lambda"); frqfade_check_len(kappa, n, "kappa");
  frqfade_check_len(delta, n, "delta"); frqfade_check_len(cv_u, n, "cv_u");
  Rcpp::NumericVector out(n);
  FrqFadeMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqFadePars& s = memo.get(frqfade_pick(alpha, i), frqfade_pick(beta, i),
      frqfade_pick(lambda, i), frqfade_pick(kappa, i), frqfade_pick(delta, i),
      frqfade_pick(cv_u, i));
    out[i] = frqfade_pdf_natural_dt(t[i], s);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector pfrqfade(Rcpp::NumericVector t,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector lambda,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector delta,
                             Rcpp::NumericVector cv_u,
                             bool lower_tail = true) {
  const int n = t.size();
  frqfade_check_len(alpha, n, "alpha"); frqfade_check_len(beta, n, "beta");
  frqfade_check_len(lambda, n, "lambda"); frqfade_check_len(kappa, n, "kappa");
  frqfade_check_len(delta, n, "delta"); frqfade_check_len(cv_u, n, "cv_u");
  Rcpp::NumericVector out(n);
  FrqFadeMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqFadePars& s = memo.get(frqfade_pick(alpha, i), frqfade_pick(beta, i),
      frqfade_pick(lambda, i), frqfade_pick(kappa, i), frqfade_pick(delta, i),
      frqfade_pick(cv_u, i));
    const double lp = lower_tail ? frqfade_log_cdf_dt(t[i], s)
                                 : frqfade_log_surv_dt(t[i], s);
    out[i] = (lp > R_NegInf) ? std::exp(lp) : 0.0;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix frq_fade_summary(Rcpp::NumericVector alpha,
                                     Rcpp::NumericVector beta,
                                     Rcpp::NumericVector lambda,
                                     Rcpp::NumericVector kappa,
                                     Rcpp::NumericVector delta,
                                     Rcpp::NumericVector cv_u) {
  const int n = std::max({alpha.size(), beta.size(), lambda.size(),
                          kappa.size(), delta.size(), cv_u.size()});
  frqfade_check_len(alpha, n, "alpha"); frqfade_check_len(beta, n, "beta");
  frqfade_check_len(lambda, n, "lambda"); frqfade_check_len(kappa, n, "kappa");
  frqfade_check_len(delta, n, "delta"); frqfade_check_len(cv_u, n, "cv_u");
  Rcpp::NumericMatrix out(n, 2);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("q_inf", "h");
  FrqFadeMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqFadePars& s = memo.get(frqfade_pick(alpha, i), frqfade_pick(beta, i),
      frqfade_pick(lambda, i), frqfade_pick(kappa, i), frqfade_pick(delta, i),
      frqfade_pick(cv_u, i));
    frqfade_ensure_endpoint(s);
    out(i, 0) = s.ok ? s.q_inf : NA_REAL;
    out(i, 1) = s.ok ? std::exp(s.log_h) : NA_REAL;
  }
  return out;
}

double dfrqfade_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq_fade::alpha]) ||
      !R_FINITE(par[emc2col::frq_fade::t0])) return 0.0;
  const double u = t - par[emc2col::frq_fade::t0];
  if (t <= 0.0 || u <= 0.0) return 0.0;
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  const FrqFadePars& s = memo.get(
    par[emc2col::frq_fade::alpha], par[emc2col::frq_fade::beta],
    par[emc2col::frq_fade::lambda], par[emc2col::frq_fade::kappa],
    par[emc2col::frq_fade::delta], par[emc2col::frq_fade::cv_u]);
  return frqfade_pdf_natural_dt(u, s);
}

double pfrqfade_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq_fade::alpha]) ||
      !R_FINITE(par[emc2col::frq_fade::t0])) return 0.0;
  const double u = t - par[emc2col::frq_fade::t0];
  if (t <= 0.0 || u <= 0.0) return 0.0;
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  const FrqFadePars& s = memo.get(
    par[emc2col::frq_fade::alpha], par[emc2col::frq_fade::beta],
    par[emc2col::frq_fade::lambda], par[emc2col::frq_fade::kappa],
    par[emc2col::frq_fade::delta], par[emc2col::frq_fade::cv_u]);
  return frqfade_cdf_natural_dt(u, s);
}

void dfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* alpha = cols[emc2col::frq_fade::alpha];
  const double* beta = cols[emc2col::frq_fade::beta];
  const double* lambda = cols[emc2col::frq_fade::lambda];
  const double* kappa = cols[emc2col::frq_fade::kappa];
  const double* t0 = cols[emc2col::frq_fade::t0];
  const double* delta = cols[emc2col::frq_fade::delta];
  const double* cv_u = cols[emc2col::frq_fade::cv_u];
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(alpha[i]) || !isok[i] || !R_FINITE(t0[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double u = rt[i] - t0[i];
    if (rt[i] <= 0.0 || u <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const FrqFadePars& s = memo.get(alpha[i], beta[i], lambda[i], kappa[i],
                                    delta[i], cv_u[i]);
    const double lp = frqfade_log_pdf_dt(u, s);
    out[i] = lp > R_NegInf && emc2_isfinite(lp)
      ? raw_log_value(lp, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

void pfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* alpha = cols[emc2col::frq_fade::alpha];
  const double* beta = cols[emc2col::frq_fade::beta];
  const double* lambda = cols[emc2col::frq_fade::lambda];
  const double* kappa = cols[emc2col::frq_fade::kappa];
  const double* t0 = cols[emc2col::frq_fade::t0];
  const double* delta = cols[emc2col::frq_fade::delta];
  const double* cv_u = cols[emc2col::frq_fade::cv_u];
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(alpha[i]) || !isok[i] || !R_FINITE(t0[i])) {
      out[i] = 0.0;
      continue;
    }
    const double u = rt[i] - t0[i];
    if (rt[i] <= 0.0 || u <= 0.0) {
      out[i] = 0.0;
      continue;
    }
    const FrqFadePars& s = memo.get(alpha[i], beta[i], lambda[i], kappa[i],
                                    delta[i], cv_u[i]);
    const double ls = frqfade_log_surv_dt(u, s);
    if (ISNAN(ls)) out[i] = raw_log_zero(min_ll, floor_raw);
    else out[i] = ls > R_NegInf ? ls : raw_log_zero(min_ll, floor_raw);
  }
}

void frqfade_logS_at_t(double t, const double* const* cols,
                       int /*n_rows_total*/, int n_lR, int /*n_par*/,
                       const int* trunc_mask, int n_unique_trials,
                       const int* isok_all, void* ctx_, double* logS_out) {
  const double* alpha = cols[emc2col::frq_fade::alpha];
  const double* beta = cols[emc2col::frq_fade::beta];
  const double* lambda = cols[emc2col::frq_fade::lambda];
  const double* kappa = cols[emc2col::frq_fade::kappa];
  const double* t0 = cols[emc2col::frq_fade::t0];
  const double* delta = cols[emc2col::frq_fade::delta];
  const double* cv_u = cols[emc2col::frq_fade::cv_u];
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double total = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(alpha[r]) || !R_FINITE(t0[r])) {
        bad = true;
        break;
      }
      const double u = t - t0[r];
      if (u <= 0.0) continue;
      const FrqFadePars& s = memo.get(alpha[r], beta[r], lambda[r], kappa[r],
                                     delta[r], cv_u[r]);
      const double ls = frqfade_log_surv_dt(u, s);
      if (ISNAN(ls)) { bad = true; break; }
      total += ls;
    }
    logS_out[j] = bad ? R_NegInf : total;
  }
}

#include "model_FRQfade.h"
#include "race_contract.h"
#include "utility_functions.h"

#include <algorithm>

// ============================================================
// FRQfade adapters
// Column layout: K=0, v=1, sv=2, kappa=3, t0=4, plo=5, clock=6.
// The base FRQ state (h = 1) goes through the shared per-particle FrqCache.
// ============================================================

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

static inline FrqCache* frqfade_cache_from_context(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr && ctx->frq_cache ? ctx->frq_cache.get() : nullptr;
}

static inline const FrqFadePars& frqfade_row(FrqFadeMemo& memo,
                                             const double* par) {
  return memo.get(par[emc2col::frq_fade::K], par[emc2col::frq_fade::v],
                  par[emc2col::frq_fade::sv], par[emc2col::frq_fade::kappa],
                  par[emc2col::frq_fade::plo], par[emc2col::frq_fade::clock]);
}

// Per-thread memo for the scalar entry points, as frq_scalar_row in
// model_FRQ.cpp: keyed on parameter values, so hits are valid across particles.
static inline const FrqFadePars& frqfade_scalar_row(const double* par,
                                                    void* ctx_) {
  static thread_local FrqFadeMemo memo;
  memo.base.shared = frqfade_cache_from_context(ctx_);
  return frqfade_row(memo, par);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector dfrqfade(Rcpp::NumericVector t, Rcpp::NumericVector K,
                             Rcpp::NumericVector v, Rcpp::NumericVector sv,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector plo = Rcpp::NumericVector::create(0.0),
                             Rcpp::NumericVector clock = Rcpp::NumericVector::create(1.0)) {
  const int n = t.size();
  frqfade_check_len(K, n, "K"); frqfade_check_len(v, n, "v");
  frqfade_check_len(sv, n, "sv"); frqfade_check_len(kappa, n, "kappa");
  frqfade_check_len(plo, n, "plo"); frqfade_check_len(clock, n, "clock");
  Rcpp::NumericVector out(n);
  FrqFadeMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqFadePars& s = memo.get(frqfade_pick(K, i), frqfade_pick(v, i),
      frqfade_pick(sv, i), frqfade_pick(kappa, i), frqfade_pick(plo, i),
      frqfade_pick(clock, i));
    out[i] = frqfade_pdf_natural_dt(t[i], s);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector pfrqfade(Rcpp::NumericVector t, Rcpp::NumericVector K,
                             Rcpp::NumericVector v, Rcpp::NumericVector sv,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector plo = Rcpp::NumericVector::create(0.0),
                             Rcpp::NumericVector clock = Rcpp::NumericVector::create(1.0),
                             bool lower_tail = true) {
  const int n = t.size();
  frqfade_check_len(K, n, "K"); frqfade_check_len(v, n, "v");
  frqfade_check_len(sv, n, "sv"); frqfade_check_len(kappa, n, "kappa");
  frqfade_check_len(plo, n, "plo"); frqfade_check_len(clock, n, "clock");
  Rcpp::NumericVector out(n);
  FrqFadeMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqFadePars& s = memo.get(frqfade_pick(K, i), frqfade_pick(v, i),
      frqfade_pick(sv, i), frqfade_pick(kappa, i), frqfade_pick(plo, i),
      frqfade_pick(clock, i));
    const double lp = lower_tail ? frqfade_log_cdf_dt(t[i], s)
                                 : frqfade_log_surv_dt(t[i], s);
    out[i] = (lp > R_NegInf) ? std::exp(lp) : 0.0;
  }
  return out;
}

// Derived completion probability h = F_1(1/kappa) and expected total
// evidence; NA for invalid rows.
// [[Rcpp::export]]
Rcpp::NumericMatrix frq_fade_summary(Rcpp::NumericVector K,
                                     Rcpp::NumericVector v,
                                     Rcpp::NumericVector sv,
                                     Rcpp::NumericVector kappa,
                                     Rcpp::NumericVector plo = Rcpp::NumericVector::create(0.0),
                                     Rcpp::NumericVector clock = Rcpp::NumericVector::create(1.0)) {
  const int n = static_cast<int>(std::max({K.size(), v.size(), sv.size(),
                                           kappa.size(), plo.size(),
                                           clock.size()}));
  frqfade_check_len(K, n, "K"); frqfade_check_len(v, n, "v");
  frqfade_check_len(sv, n, "sv"); frqfade_check_len(kappa, n, "kappa");
  frqfade_check_len(plo, n, "plo"); frqfade_check_len(clock, n, "clock");
  Rcpp::NumericMatrix out(n, 2);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("h", "Etot");
  FrqFadeMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqFadePars& s = memo.get(frqfade_pick(K, i), frqfade_pick(v, i),
      frqfade_pick(sv, i), frqfade_pick(kappa, i), frqfade_pick(plo, i),
      frqfade_pick(clock, i));
    out(i, 0) = s.ok ? std::exp(frqfade_log_cdf_dt(R_PosInf, s)) : NA_REAL;
    out(i, 1) = s.ok ? frqfade_total_evidence(s) : NA_REAL;
  }
  return out;
}

double dfrqfade_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq_fade::K]) ||
      !R_FINITE(par[emc2col::frq_fade::t0])) return 0.0;
  const double u = t - par[emc2col::frq_fade::t0];
  if (t <= 0.0 || u <= 0.0) return 0.0;
  return frqfade_pdf_natural_dt(u, frqfade_scalar_row(par, ctx_));
}

double pfrqfade_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq_fade::K]) ||
      !R_FINITE(par[emc2col::frq_fade::t0])) return 0.0;
  const double u = t - par[emc2col::frq_fade::t0];
  if (t <= 0.0 || u <= 0.0) return 0.0;
  // u == Inf reaches the kernel: the CDF there is the derived h.
  return frqfade_cdf_natural_dt(u, frqfade_scalar_row(par, ctx_));
}

void dfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* k_  = cols[emc2col::frq_fade::K];
  const double* v_  = cols[emc2col::frq_fade::v];
  const double* sv_ = cols[emc2col::frq_fade::sv];
  const double* ka_ = cols[emc2col::frq_fade::kappa];
  const double* t0_ = cols[emc2col::frq_fade::t0];
  const double* pl_ = cols[emc2col::frq_fade::plo];
  const double* cu_ = cols[emc2col::frq_fade::clock];
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(k_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double u = rt[i] - t0_[i];
    if (rt[i] <= 0.0 || u <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const FrqFadePars& s = memo.get(k_[i], v_[i], sv_[i], ka_[i], pl_[i], cu_[i]);
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
  const double* k_  = cols[emc2col::frq_fade::K];
  const double* v_  = cols[emc2col::frq_fade::v];
  const double* sv_ = cols[emc2col::frq_fade::sv];
  const double* ka_ = cols[emc2col::frq_fade::kappa];
  const double* t0_ = cols[emc2col::frq_fade::t0];
  const double* pl_ = cols[emc2col::frq_fade::plo];
  const double* cu_ = cols[emc2col::frq_fade::clock];
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(k_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = 0.0;
      continue;
    }
    const double u = rt[i] - t0_[i];
    if (rt[i] <= 0.0 || u <= 0.0) {
      out[i] = 0.0;
      continue;
    }
    const FrqFadePars& s = memo.get(k_[i], v_[i], sv_[i], ka_[i], pl_[i], cu_[i]);
    const double ls = frqfade_log_surv_dt(u, s);
    if (ISNAN(ls)) out[i] = raw_log_zero(min_ll, floor_raw);
    else out[i] = ls > R_NegInf ? ls : raw_log_zero(min_ll, floor_raw);
  }
}

void frqfade_logS_at_t(double t, const double* const* cols,
                       int /*n_rows_total*/, int n_lR, int /*n_par*/,
                       const int* trunc_mask, int n_unique_trials,
                       const int* isok_all, void* ctx_, double* logS_out) {
  const double* k_  = cols[emc2col::frq_fade::K];
  const double* v_  = cols[emc2col::frq_fade::v];
  const double* sv_ = cols[emc2col::frq_fade::sv];
  const double* ka_ = cols[emc2col::frq_fade::kappa];
  const double* t0_ = cols[emc2col::frq_fade::t0];
  const double* pl_ = cols[emc2col::frq_fade::plo];
  const double* cu_ = cols[emc2col::frq_fade::clock];
  FrqFadeMemo memo(frqfade_cache_from_context(ctx_));
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double total = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(k_[r]) || !R_FINITE(t0_[r])) {
        bad = true;
        break;
      }
      const double u = t - t0_[r];
      if (u <= 0.0) continue;
      const FrqFadePars& s = memo.get(k_[r], v_[r], sv_[r], ka_[r], pl_[r], cu_[r]);
      const double ls = frqfade_log_surv_dt(u, s);
      if (ISNAN(ls)) { bad = true; break; }
      total += ls;
    }
    logS_out[j] = bad ? R_NegInf : total;
  }
}

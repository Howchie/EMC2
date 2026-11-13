#include <Rcpp.h>
#include "model_EXG.h"
#include "race_contract.h"
#include "col_registry.h"
#include "utility_functions.h"
#include "exgaussian_functions.h"
#include "simd_math.h"
#include <array>
#include <cmath>
#include <vector>

namespace {

struct RexgTruncationCacheSet {
  static constexpr size_t max_entries = 16;
  std::array<ExgTruncationCache, max_entries> entries;
  size_t size = 0;

  const ExgTruncationCache* find(double mu, double sigma, double tau) const {
    for (size_t i = 0; i < size; ++i) {
      const ExgTruncationCache& entry = entries[i];
      if (entry.matches(mu, sigma, tau, 0.0, R_PosInf)) return &entry;
    }
    return nullptr;
  }
};

bool rexg_truncation_caches(const double* mu, const double* sigma,
                            const double* tau, const double* t0,
                            int n_rows, const int* mask, const int* isok,
                            RexgTruncationCacheSet& caches) {
  caches.size = 0;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i] || !isok[i] || std::isnan(mu[i]) || std::isnan(sigma[i]) ||
        std::isnan(tau[i]) || std::isnan(t0[i]) || sigma[i] <= 0.0 || tau[i] <= 0.0 ||
        !R_FINITE(t0[i]) || t0[i] < 0.0) continue;
    if (caches.find(mu[i], sigma[i], tau[i]) != nullptr) continue;
    if (caches.size == RexgTruncationCacheSet::max_entries) {
      caches.size = 0;
      return false;
    }
    caches.entries[caches.size++] = make_exg_truncation_cache(
        mu[i], sigma[i], tau[i], 0.0, R_PosInf);
  }
  return caches.size > 0;
}

bool rexg_logS_truncation_caches(const double* mu, const double* sigma,
                                 const double* tau, const double* t0,
                                 int n_lR, int n_unique_trials,
                                 const int* trunc_mask, const int* isok,
                                 RexgTruncationCacheSet& caches) {
  caches.size = 0;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    for (int k = 0; k < n_lR; ++k) {
      const int r = start + k;
      if (!isok[r] || std::isnan(mu[r]) || std::isnan(sigma[r]) || std::isnan(tau[r]) ||
          std::isnan(t0[r]) || sigma[r] <= 0.0 || tau[r] <= 0.0 ||
          !R_FINITE(t0[r]) || t0[r] < 0.0) continue;
      if (caches.find(mu[r], sigma[r], tau[r]) != nullptr) continue;
      if (caches.size == RexgTruncationCacheSet::max_entries) {
        caches.size = 0;
        return false;
      }
      caches.entries[caches.size++] = make_exg_truncation_cache(
          mu[r], sigma[r], tau[r], 0.0, R_PosInf);
    }
  }
  return caches.size > 0;
}

const ExgTruncationCache* rexg_scalar_truncation_cache(
    double mu, double sigma, double tau, double lower, double upper,
    void* ctx_) {
  if (ctx_ == nullptr)
    return exg_scalar_truncation_memo().get(mu, sigma, tau, lower, upper);
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (!ctx->exg_truncation_cache)
    ctx->exg_truncation_cache = std::make_shared<ExgTruncationCache>();
  return exg_truncation_cache_get(*ctx->exg_truncation_cache,
                                  mu, sigma, tau, lower, upper);
}

// Dense per-branch workspace for the shared-parameter density pass.
struct RexgDenseScratch {
  std::vector<double> neg_z, neg_e, neg_base, pos_arg, pos_base;
  std::vector<int> neg_row, pos_row;
};

RexgDenseScratch& rexg_dense_scratch() {
  static thread_local RexgDenseScratch scratch;
  return scratch;
}

// Log density for every row when all valid rows share one parameter tuple.
// Rows are gathered by branch of exg_log_weighted_tail, the exp/log1p/log
// calls run over each dense group, and the results are scattered back; rows
// outside the two common branches take dexg_prepared as before.
void drexg_shared_dense(const double* rt, const double* mu_, const double* sigma_,
                        const double* tau_, const double* t0_, int n_rows,
                        const int* mask, const int* isok,
                        const ExgTruncationCache& cache,
                        double* out, double min_ll, bool floor_raw) {
  RexgDenseScratch& ws = rexg_dense_scratch();
  const size_t cap = static_cast<size_t>(n_rows);
  if (ws.neg_z.size() < cap) {
    ws.neg_z.resize(cap); ws.neg_base.resize(cap); ws.neg_row.resize(cap);
    ws.pos_arg.resize(cap); ws.pos_base.resize(cap); ws.pos_row.resize(cap);
  }
  int n_neg = 0, n_pos = 0;
  const double ratio = cache.sigma_tau_ratio;
  const bool ratio_finite = R_FINITE(ratio);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (std::isnan(mu_[i]) || std::isnan(sigma_[i]) || std::isnan(tau_[i]) || std::isnan(t0_[i]) ||
        !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0 ||
        !R_FINITE(t0_[i]) || t0_[i] < 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0) || !R_FINITE(tt)) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double y = (tt - cache.mu) * cache.inverse_sigma;
    const double w = ratio - y;
    if (ratio_finite && R_FINITE(y)) {
      if (w < 0.0) {
        // log Phi(-w) on its bulk side: log1p(-Q(-w)).
        ws.neg_z[n_neg] = -w;
        ws.neg_base[n_neg] = ratio * (0.5 * ratio - y);
        ws.neg_row[n_neg] = i;
        ++n_neg;
        continue;
      }
      // log H(w): the central branch logs its rational factor alone; the
      // other branches are rational, so only the log remains.
      double add = 0.0, arg;
      if (w <= 0.67448975) {
        const double u = w * w;
        const double numerator = (((0.065682337918207449113 * u +
          2.2352520354606839287) * u + 161.02823106855587881) * u +
          1067.6894854603709582) * u + 18154.981253343561249;
        const double denominator = (((u + 47.20258190468824187) * u +
          976.09855173777669322) * u + 10260.932208618978205) * u +
          45507.789335026729956;
        add = 0.5 * u;
        arg = 0.5 - w * numerator / denominator;
      } else {
        arg = exg_normal_tail_factor(w);
      }
      ws.pos_arg[n_pos] = arg;
      ws.pos_base[n_pos] = -0.5 * y * y + add;
      ws.pos_row[n_pos] = i;
      ++n_pos;
      continue;
    }
    const double log_density = dexg_prepared(
        tt, cache.mu, cache.inverse_sigma, cache.tau_prepared,
        cache.sigma_tau_ratio, cache.log_sigma_tau_ratio, cache.log_tau,
        true) - cache.log_normalizer;
    out[i] = raw_log_value(log_density, min_ll, floor_raw);
  }

  // w < 0: Q(z) for z = -w > 0 as in fast_norm_phi, then log1p(-Q).
  if (n_neg > 0) {
    const double* z = ws.neg_z.data();
    ws.neg_e.resize(static_cast<size_t>(n_neg));
    double* e = ws.neg_e.data();
    for (int k = 0; k < n_neg; ++k) e[k] = z[k] <= 37.0 ? -z[k] * z[k] / 2.0 : 0.0;
    emc2_simd_exp(e, n_neg);
#pragma omp simd
    for (int k = 0; k < n_neg; ++k) {
      const double zk = z[k];
      const double n = (((((FAST_NORM_N6 * zk + FAST_NORM_N5) * zk + FAST_NORM_N4) * zk + FAST_NORM_N3) * zk + FAST_NORM_N2) * zk + FAST_NORM_N1) * zk + FAST_NORM_N0;
      const double d = ((((((FAST_NORM_M7 * zk + FAST_NORM_M6) * zk + FAST_NORM_M5) * zk + FAST_NORM_M4) * zk + FAST_NORM_M3) * zk + FAST_NORM_M2) * zk + FAST_NORM_M1) * zk + FAST_NORM_M0;
      const double f = zk + 1.0 / (zk + 2.0 / (zk + 3.0 / (zk + 4.0 / (zk + 13.0 / 20.0))));
      const double q = zk < FAST_NORM_SPLIT ? e[k] * n / d : e[k] / (FAST_NORM_RT2PI * f);
      e[k] = zk <= 37.0 ? -q : 0.0;
    }
    emc2_simd_log1p(e, n_neg);
    const double* base = ws.neg_base.data();
    const int* row = ws.neg_row.data();
    for (int k = 0; k < n_neg; ++k) {
      const double log_density = (base[k] + e[k]) - cache.log_tau - cache.log_normalizer;
      out[row[k]] = raw_log_value(log_density, min_ll, floor_raw);
    }
  }

  if (n_pos > 0) {
    double* arg = ws.pos_arg.data();
    emc2_simd_log(arg, n_pos);
    const double* base = ws.pos_base.data();
    const int* row = ws.pos_row.data();
    for (int k = 0; k < n_pos; ++k) {
      const double log_density = (base[k] + arg[k]) - cache.log_tau - cache.log_normalizer;
      out[row[k]] = raw_log_value(log_density, min_ll, floor_raw);
    }
  }
}

}  // namespace

void drexg_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* mu_    = cols[emc2col::rexg::mu];
  const double* sigma_ = cols[emc2col::rexg::sigma];
  const double* tau_   = cols[emc2col::rexg::tau];
  const double* t0_    = cols[emc2col::rexg::t0];
  RexgTruncationCacheSet norm_caches;
  const bool have_norm_caches = rexg_truncation_caches(
      mu_, sigma_, tau_, t0_, n_rows, mask, isok, norm_caches);
  if (have_norm_caches && norm_caches.size == 1) {
    drexg_shared_dense(rt, mu_, sigma_, tau_, t0_, n_rows, mask, isok,
                       norm_caches.entries[0], out, min_ll, floor_raw);
    return;
  }
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (std::isnan(mu_[i]) || std::isnan(sigma_[i]) || std::isnan(tau_[i]) || std::isnan(t0_[i]) ||
        !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0 ||
        !R_FINITE(t0_[i]) || t0_[i] < 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0) || !R_FINITE(tt)) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const ExgTruncationCache* norm_cache =
      have_norm_caches ? norm_caches.find(mu_[i], sigma_[i], tau_[i]) : nullptr;
    out[i] = raw_log_value(dtexg(tt, mu_[i], sigma_[i], tau_[i], 0.0, R_PosInf,
                                 true, norm_cache),
                           min_ll, floor_raw);
  }
}

void prexg_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* mu_    = cols[emc2col::rexg::mu];
  const double* sigma_ = cols[emc2col::rexg::sigma];
  const double* tau_   = cols[emc2col::rexg::tau];
  const double* t0_    = cols[emc2col::rexg::t0];
  RexgTruncationCacheSet norm_caches;
  const bool have_norm_caches = rexg_truncation_caches(
      mu_, sigma_, tau_, t0_, n_rows, mask, isok, norm_caches);
  if (have_norm_caches && norm_caches.size == 1) {
    const ExgTruncationCache& cache = norm_caches.entries[0];
#pragma omp simd
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (std::isnan(mu_[i]) || std::isnan(sigma_[i]) || std::isnan(tau_[i]) || std::isnan(t0_[i]) ||
          !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0 ||
          !R_FINITE(t0_[i]) || t0_[i] < 0.0) {
        out[i] = 0.0;
        continue;
      }
      const double tt = rt[i] - t0_[i];
      const double log_surv = tt <= 0.0
        ? 0.0
        : ptexg(tt, cache.mu, cache.sigma, cache.tau, 0.0, R_PosInf,
                false, true, &cache);
      out[i] = R_FINITE(log_surv) ? std::fmin(log_surv, 0.0)
                                  : (log_surv == R_NegInf
                                       ? raw_log_zero(min_ll, floor_raw) : 0.0);
    }
    return;
  }
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (std::isnan(mu_[i]) || std::isnan(sigma_[i]) || std::isnan(tau_[i]) || std::isnan(t0_[i]) ||
        !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0 ||
        !R_FINITE(t0_[i]) || t0_[i] < 0.0) {
      out[i] = 0.0;
      continue;
    }
    const ExgTruncationCache* norm_cache =
      have_norm_caches ? norm_caches.find(mu_[i], sigma_[i], tau_[i]) : nullptr;
    const double log_surv = (rt[i] - t0_[i] <= 0.0)
      ? 0.0 : ptexg(rt[i] - t0_[i], mu_[i], sigma_[i], tau_[i], 0.0,
                    R_PosInf, false, true, norm_cache);
    out[i] = R_FINITE(log_surv) ? std::fmin(log_surv, 0.0)
                                : (log_surv == R_NegInf ? raw_log_zero(min_ll, floor_raw) : 0.0);
  }
}

void rexg_logS_at_t(double t, const double* const* cols,
                    int /*n_rows_total*/, int n_lR, int /*n_par*/,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* mu_    = cols[emc2col::rexg::mu];
  const double* sigma_ = cols[emc2col::rexg::sigma];
  const double* tau_   = cols[emc2col::rexg::tau];
  const double* t0_    = cols[emc2col::rexg::t0];
  RexgTruncationCacheSet norm_caches;
  const bool have_norm_caches = rexg_logS_truncation_caches(
      mu_, sigma_, tau_, t0_, n_lR, n_unique_trials, trunc_mask, isok_all,
      norm_caches);
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || std::isnan(mu_[r]) || std::isnan(sigma_[r]) || std::isnan(tau_[r]) ||
          std::isnan(t0_[r]) || !R_FINITE(t0_[r]) ||
          sigma_[r] <= 0.0 || tau_[r] <= 0.0 || t0_[r] < 0.0) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;
      const ExgTruncationCache* norm_cache = have_norm_caches
        ? norm_caches.find(mu_[r], sigma_[r], tau_[r]) : nullptr;
      const double log_surv = ptexg(tt, mu_[r], sigma_[r], tau_[r], 0.0,
                                    R_PosInf, false, true, norm_cache);
      if (!R_FINITE(log_surv)) {
        if (log_surv == R_NegInf) { bad = true; break; }
        continue;
      }
      logS += std::fmin(log_surv, 0.0);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

double dexg_scalar(double x, const double* par, void* ctx) {
  const double t0 = par[3];
  if (!R_FINITE(t0) || t0 < 0.0) return 0.0;
  const double tt = x - t0;
  if (tt <= 0.0) return 0.0;
  const ExgTruncationCache* cache = rexg_scalar_truncation_cache(
      par[0], par[1], par[2], 0.0, R_PosInf, ctx);
  return dtexg(tt, par[0], par[1], par[2], 0.0, R_PosInf, false, cache);
}

double pexg_scalar(double x, const double* par, void* ctx) {
  const double t0 = par[3];
  if (!R_FINITE(t0) || t0 < 0.0) return 0.0;
  const double tt = x - t0;
  if (tt <= 0.0) return 0.0;
  const ExgTruncationCache* cache = rexg_scalar_truncation_cache(
      par[0], par[1], par[2], 0.0, R_PosInf, ctx);
  return ptexg(tt, par[0], par[1], par[2], 0.0, R_PosInf, true, false, cache);
}

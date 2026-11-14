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

// Whether every masked row is valid and carries one parameter tuple; if so
// `first` is one such row.
bool rexg_rows_uniform(const double* mu, const double* sigma, const double* tau,
                       const double* t0, int n_rows, const int* mask,
                       const int* isok, int& first) {
  first = 0;
  while (first < n_rows && !mask[first]) ++first;
  if (first == n_rows) return false;
  const double m = mu[first], sg = sigma[first], ta = tau[first], on = t0[first];
  if (!isok[first] || std::isnan(m) || !(sg > 0.0) || !(ta > 0.0) ||
      !emc2_isfinite(on) || on < 0.0) return false;
  int differ = 0;
#pragma omp simd reduction(+:differ)
  for (int i = first; i < n_rows; ++i)
    differ += (mask[i] != 0) & ((isok[i] == 0) | (mu[i] != m) | (sigma[i] != sg) |
                                (tau[i] != ta) | (t0[i] != on));
  return differ == 0;
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
  std::vector<double> w[4], y[4], e;
  std::vector<int> row[4];
};

RexgDenseScratch& rexg_dense_scratch() {
  static thread_local RexgDenseScratch scratch;
  return scratch;
}

// Log density for every row when all valid rows share one parameter tuple.
// Rows are gathered by branch of exg_log_weighted_tail (w < 0, then the
// three ranges of exg_normal_tail_factor), each branch's rational and its
// exp/log1p/log run over a dense group, and the results are scattered back;
// rows outside those branches take dexg_prepared. `rows_valid` says every
// masked row already passed the parameter checks.
void drexg_shared_dense(const double* rt, const double* mu_, const double* sigma_,
                        const double* tau_, const double* t0_, int n_rows,
                        const int* mask, const int* isok,
                        const ExgTruncationCache& cache,
                        double* out, double min_ll, bool floor_raw,
                        bool rows_valid = false) {
  RexgDenseScratch& ws = rexg_dense_scratch();
  const size_t cap = static_cast<size_t>(n_rows);
  if (ws.e.size() < cap) {
    ws.e.resize(cap);
    for (int c = 0; c < 4; ++c) {
      ws.w[c].resize(cap); ws.y[c].resize(cap); ws.row[c].resize(cap);
    }
  }
  double* const wv[4] = {ws.w[0].data(), ws.w[1].data(), ws.w[2].data(), ws.w[3].data()};
  double* const yv[4] = {ws.y[0].data(), ws.y[1].data(), ws.y[2].data(), ws.y[3].data()};
  int* const rv[4] = {ws.row[0].data(), ws.row[1].data(), ws.row[2].data(), ws.row[3].data()};
  int n[4] = {0, 0, 0, 0};
  const double ratio = cache.sigma_tau_ratio;
  const bool ratio_finite = R_FINITE(ratio);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!rows_valid &&
        (std::isnan(mu_[i]) || std::isnan(sigma_[i]) || std::isnan(tau_[i]) || std::isnan(t0_[i]) ||
         !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0 ||
         !R_FINITE(t0_[i]) || t0_[i] < 0.0)) {
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
      const int c = (w >= 0.0) + (w > 0.67448975) + (w > 5.6568542494923801952);
      const int k = n[c]++;
      wv[c][k] = w; yv[c][k] = y; rv[c][k] = i;
      continue;
    }
    const double log_density = dexg_prepared(
        tt, cache.mu, cache.inverse_sigma, cache.tau_prepared,
        cache.sigma_tau_ratio, cache.log_sigma_tau_ratio, cache.log_tau,
        true) - cache.log_normalizer;
    out[i] = raw_log_value(log_density, min_ll, floor_raw);
  }

  // w < 0: log Phi(-w) on its bulk side, log1p(-Q(z)) for z = -w > 0 with Q
  // as in fast_norm_phi.
  if (n[0] > 0) {
    const int m = n[0];
    double* z = wv[0];
    double* base = yv[0];
    double* e = ws.e.data();
#pragma omp simd
    for (int k = 0; k < m; ++k) {
      z[k] = -z[k];
      base[k] = ratio * (0.5 * ratio - base[k]);
      e[k] = z[k] <= 37.0 ? -z[k] * z[k] / 2.0 : 0.0;
    }
    emc2_simd_exp(e, m);
#pragma omp simd
    for (int k = 0; k < m; ++k) {
      const double zk = z[k];
      const double nn = (((((FAST_NORM_N6 * zk + FAST_NORM_N5) * zk + FAST_NORM_N4) * zk + FAST_NORM_N3) * zk + FAST_NORM_N2) * zk + FAST_NORM_N1) * zk + FAST_NORM_N0;
      const double d = ((((((FAST_NORM_M7 * zk + FAST_NORM_M6) * zk + FAST_NORM_M5) * zk + FAST_NORM_M4) * zk + FAST_NORM_M3) * zk + FAST_NORM_M2) * zk + FAST_NORM_M1) * zk + FAST_NORM_M0;
      const double f = zk + 1.0 / (zk + 2.0 / (zk + 3.0 / (zk + 4.0 / (zk + 13.0 / 20.0))));
      const double q = zk < FAST_NORM_SPLIT ? e[k] * nn / d : e[k] / (FAST_NORM_RT2PI * f);
      e[k] = zk <= 37.0 ? -q : 0.0;
    }
    emc2_simd_log1p(e, m);
    const int* row = rv[0];
    for (int k = 0; k < m; ++k) {
      const double log_density = (base[k] + e[k]) - cache.log_tau - cache.log_normalizer;
      out[row[k]] = raw_log_value(log_density, min_ll, floor_raw);
    }
  }

  // w >= 0: log H(w) = -y^2/2 + log of the scaled tail's rational factor
  // (the central range keeps its exp(w^2/2) outside the log).
  if (n[1] > 0) {
    double* arg = wv[1];
    double* base = yv[1];
#pragma omp simd
    for (int k = 0; k < n[1]; ++k) {
      const double x = arg[k];
      const double u = x * x;
      const double numerator = (((0.065682337918207449113 * u +
        2.2352520354606839287) * u + 161.02823106855587881) * u +
        1067.6894854603709582) * u + 18154.981253343561249;
      const double denominator = (((u + 47.20258190468824187) * u +
        976.09855173777669322) * u + 10260.932208618978205) * u +
        45507.789335026729956;
      base[k] = -0.5 * base[k] * base[k] + 0.5 * u;
      arg[k] = 0.5 - x * numerator / denominator;
    }
  }
  if (n[2] > 0) {
    double* arg = wv[2];
    double* base = yv[2];
#pragma omp simd
    for (int k = 0; k < n[2]; ++k) {
      const double x = arg[k];
      const double numerator = (((((((1.0765576773720192317e-8 * x +
        0.39894151208813466764) * x + 8.8831497943883759412) * x +
        93.506656132177855979) * x + 597.27027639480026226) * x +
        2494.5375852903726711) * x + 6848.1904505362823326) * x +
        11602.651437647350124) * x + 9842.7148383839780218;
      const double denominator = (((((((x + 22.266688044328115691) * x +
        235.38790178262499861) * x + 1519.377599407554805) * x +
        6485.558298266760755) * x + 18615.571640885098091) * x +
        34900.952721145977266) * x + 38912.003286093271411) * x +
        19685.429676859990727;
      base[k] = -0.5 * base[k] * base[k];
      arg[k] = numerator / denominator;
    }
  }
  if (n[3] > 0) {
    double* arg = wv[3];
    double* base = yv[3];
#pragma omp simd
    for (int k = 0; k < n[3]; ++k) {
      const double inv = 1.0 / arg[k];
      const double u = inv * inv;
      const double numerator = ((((0.02307344176494017303 * u +
        0.21589853405795699) * u + 0.1274011611602473639) * u +
        0.022235277870649807) * u + 0.001421619193227893466) * u +
        2.9112874951168792e-5;
      const double denominator = ((((u + 1.28426009614491121) * u +
        0.468238212480865118) * u + 0.0659881378689285515) * u +
        0.00378239633202758244) * u + 7.29751555083966205e-5;
      base[k] = -0.5 * base[k] * base[k];
      arg[k] = (0.39894228040143267794 - u * (numerator / denominator)) * inv;
    }
  }
  for (int c = 1; c < 4; ++c) {
    if (n[c] == 0) continue;
    double* arg = wv[c];
    emc2_simd_log(arg, n[c]);
    const double* base = yv[c];
    const int* row = rv[c];
    for (int k = 0; k < n[c]; ++k) {
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
  int first = 0;
  if (rexg_rows_uniform(mu_, sigma_, tau_, t0_, n_rows, mask, isok, first)) {
    const ExgTruncationCache cache = make_exg_truncation_cache(
        mu_[first], sigma_[first], tau_[first], 0.0, R_PosInf);
    drexg_shared_dense(rt, mu_, sigma_, tau_, t0_, n_rows, mask, isok,
                       cache, out, min_ll, floor_raw, true);
    return;
  }
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

#include "utility_functions.h"
#include "race_contract.h"
#include "col_registry.h"
#include "model_RLF.h"
#include "model_RLF_kernels.h"

void drlf_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* context) {
  const bool floor_raw = raw_floor_log_lik(context);
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) {
    for (int i = 0; i < n_rows; ++i) {
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    }
    return;
  }
  auto* ctx = static_cast<ContextForRaceModels*>(context);
  rlf_prepare_rows(*cache, rt, cols, n_rows, isok,
                   ctx != nullptr ? ctx->endpoint_queries : nullptr);
  const double* t0 = cols[emc2col::rlf::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int group = cache->row_group[i];
    if (group < 0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    out[i] = raw_log_value(
      rlf::entry_log_pdf(cache->entries[group], rt[i] - t0[i]),
      min_ll, floor_raw);
  }
}

void prlf_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* context) {
  const bool floor_raw = raw_floor_log_lik(context);
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  auto* ctx = static_cast<ContextForRaceModels*>(context);
  rlf_prepare_rows(*cache, rt, cols, n_rows, isok,
                   ctx != nullptr ? ctx->endpoint_queries : nullptr);
  const double* t0 = cols[emc2col::rlf::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int group = cache->row_group[i];
    if (group < 0) {
      out[i] = 0.0;
      continue;
    }
    out[i] = raw_log_value(
      rlf::entry_log_S(cache->entries[group], rt[i] - t0[i]),
      min_ll, floor_raw);
  }
}

double drlf_scalar(double time, const double* par, void* context) {
  if (R_IsNA(par[emc2col::rlf::v])) return 0.0;
  const double tt = time - par[emc2col::rlf::t0];
  if (!(tt > 0.0)) return 0.0;
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) return 0.0;
  const double horizon = rlf_scalar_horizon(tt);
  rlf::Key key;
  if (!rlf_key_from_par(par, horizon, cache->grid, key)) return 0.0;
  const int group = rlf::cache_get(*cache, key, horizon, tt);
  if (group < 0) return 0.0;  // no usable solve; treated as zero density
  const double log_pdf =
    rlf::entry_log_pdf(cache->entries[group], tt);
  return log_pdf <= rlf::RLF_LOG_FLOOR ? 0.0 : std::exp(log_pdf);
}

double prlf_scalar(double time, const double* par, void* context) {
  if (R_IsNA(par[emc2col::rlf::v])) return 0.0;
  const double tt = time - par[emc2col::rlf::t0];
  if (!(tt > 0.0)) return 0.0;
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) return 0.0;
  const double horizon = rlf_scalar_horizon(tt);
  rlf::Key key;
  if (!rlf_key_from_par(par, horizon, cache->grid, key)) return 0.0;
  const int group = rlf::cache_get(*cache, key, horizon, tt);
  if (group < 0) return 0.0;  // no usable solve; treated as zero CDF
  const double log_survivor =
    rlf::entry_log_S(cache->entries[group], tt);
  if (log_survivor >= 0.0) return 0.0;
  if (log_survivor <= rlf::RLF_LOG_FLOOR) return 1.0;
  return -std::expm1(log_survivor);
}

void rlf_logS_at_t(double time, const double* const* cols,
                   int n_rows_total, int n_lR, int /*n_par*/,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* context,
                   double* logS_out) {
  rlf::SolveCache* cache = rlf_cache(context);
  const double* v = cols[emc2col::rlf::v];
  const double* B = cols[emc2col::rlf::B];
  const double* A = cols[emc2col::rlf::A];
  const double* t0 = cols[emc2col::rlf::t0];
  const double* s = cols[emc2col::rlf::s];
  const double* alpha = cols[emc2col::rlf::alpha];

  for (int trial = 0; trial < n_unique_trials; ++trial) {
    if (!trunc_mask[trial]) continue;
    const int start = trial * n_lR;
    double total = 0.0;
    bool bad = false;
    for (int accumulator = 0; accumulator < n_lR && !bad; ++accumulator) {
      const int row = start + accumulator;
      if (row >= n_rows_total) break;
      if (!isok_all[row] || R_IsNA(v[row])) {
        bad = true;
        break;
      }
      const double tt = time - t0[row];
      if (!(tt > 0.0)) continue;
      rlf::Key key;
      if (!rlf::rlf_key(v[row], s[row], alpha[row], B[row], A[row], key)) {
        bad = true;
        break;
      }
      if (cache->grid.horizon_split) {
        key.bucket = rlf::rlf_horizon_bucket(key, tt);
      }
      const int group = rlf::cache_get(*cache, key, tt, tt);
      if (group < 0) {
        // No usable solve for this accumulator: the truncation normaliser for
        // the trial is undefined, so flag it the same way a bad key is.
        bad = true;
        break;
      }
      const double log_survivor =
        rlf::entry_log_S(cache->entries[group], tt);
      if (log_survivor <= rlf::RLF_LOG_FLOOR) {
        bad = true;
      } else if (log_survivor < 0.0) {
        total += log_survivor;
      }
    }
    logS_out[trial] = bad ? R_NegInf : total;
  }
}

#include "utility_functions.h"
#include "race_contract.h"
#include "col_registry.h"
#include "model_ROU.h"
#include "model_GOM.h"

void dgomp_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok, double* out,
               double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = gomp_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i)
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    return;
  }
  gomp_prepare_rows(*C, rt, cols, n_rows, isok);
  const double* t0_ = cols[emc2col::gompertz::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = raw_log_value(fperace::entry_log_pdf(C->e[g], rt[i] - t0_[i]),
                           min_ll, floor_raw);
  }
}

void pgomp_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok, double* out,
               double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = gomp_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  gomp_prepare_rows(*C, rt, cols, n_rows, isok);
  const double* t0_ = cols[emc2col::gompertz::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    out[i] = (g < 0) ? 0.0 : raw_log_value(
      fperace::entry_log_S(C->e[g], rt[i] - t0_[i]), min_ll, floor_raw);
  }
}

double dgomp_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = gomp_cache(ctx_);
  if (C == nullptr || R_IsNA(par[emc2col::gompertz::alpha])) return 0.0;
  const double tt = t - par[emc2col::gompertz::t0];
  if (!(tt > 0.0)) return 0.0;
  fperace::Key p;
  const GompColIdx ci = gomp_col_idx();
  if (!gomp_key(par[ci.alpha], par[ci.beta], par[ci.K], par[ci.B], par[ci.A],
                gomp_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, (tt > 1.0) ? 2.0 * tt : tt + 1.0);
  const double lp = fperace::entry_log_pdf(C->e[g], tt);
  return lp <= fperace::LOG_FLOOR ? 0.0 : std::exp(lp);
}

double pgomp_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = gomp_cache(ctx_);
  if (C == nullptr || R_IsNA(par[emc2col::gompertz::alpha])) return 0.0;
  const double tt = t - par[emc2col::gompertz::t0];
  if (!(tt > 0.0)) return 0.0;
  fperace::Key p;
  const GompColIdx ci = gomp_col_idx();
  if (!gomp_key(par[ci.alpha], par[ci.beta], par[ci.K], par[ci.B], par[ci.A],
                gomp_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, (tt > 1.0) ? 2.0 * tt : tt + 1.0);
  const double lS = fperace::entry_log_S(C->e[g], tt);
  if (lS >= 0.0) return 0.0;
  if (lS <= fperace::LOG_FLOOR) return 1.0;
  return -std::expm1(lS);
}

void gomp_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int /*n_par*/,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out) {
  fperace::SolveCache* C = gomp_cache(ctx_);
  if (C == nullptr) return;
  const GompColIdx ci = gomp_col_idx();
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double total = 0.0;
    bool bad = false;
    for (int m = 0; m < n_lR && !bad; ++m) {
      const int r = j * n_lR + m;
      if (r >= n_rows_total) break;
      if (!isok_all[r]) { bad = true; break; }
      const double tt = t - cols[ci.t0][r];
      if (!(tt > 0.0)) continue;
      fperace::Key p;
      if (!gomp_key(cols[ci.alpha][r], cols[ci.beta][r], cols[ci.K][r],
                    cols[ci.B][r], cols[ci.A][r],
                    gomp_bnd_row(C->bnd_kind, cols, r), p)) {
        bad = true; break;
      }
      const int g = fperace::cache_get(*C, p, tt);
      const double lS = fperace::entry_log_S(C->e[g], tt);
      if (lS <= fperace::LOG_FLOOR) { bad = true; break; }
      if (lS < 0.0) total += lS;
    }
    logS_out[j] = bad ? R_NegInf : total;
  }
}

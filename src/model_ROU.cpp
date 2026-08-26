#include "utility_functions.h"
#include "race_contract.h"
#include "col_registry.h"
#include "model_ROU.h"

void drou_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) {
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    }
    return;
  }
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  rou_prepare_rows(*C, rt, cols, n_rows, isok,
                   ctx != nullptr ? ctx->endpoint_queries : nullptr);

  // The solve is in the rescaled state Y = X/s, whose first-passage time is the
  // SAME random variable, so the density in t needs no Jacobian -- exactly as
  // the RDM divides its arguments by s and applies none.
  const double* t0_ = cols[rou_col_idx(C->par_kind).t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double lp = fperace::entry_log_pdf(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(lp, min_ll, floor_raw);
  }
}

void prou_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  rou_prepare_rows(*C, rt, cols, n_rows, isok,
                   ctx != nullptr ? ctx->endpoint_queries : nullptr);

  const double* t0_ = cols[rou_col_idx(C->par_kind).t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) {
      // Before t0, or an unusable row.  An accumulator that cannot have
      // finished has survived, exactly as the analytic kernels report.
      out[i] = 0.0;
      continue;
    }
    const double ls = fperace::entry_log_S(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(ls, min_ll, floor_raw);
  }
}

void rou_logS_at_t(double t, const double* const* cols,
                   int n_rows_total, int n_lR, int /*n_par*/,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* ctx_, double* logS_out) {
  fperace::SolveCache* C = rou_cache(ctx_);
  const RouColIdx ci = rou_col_idx(C->par_kind);
  const double* p1_ = cols[ci.p1];
  const double* p2_ = cols[ci.p2];
  const double* p3_ = cols[ci.p3];
  const double* B_  = cols[ci.B];
  const double* A_  = cols[ci.A];
  const double* t0_ = cols[ci.t0];

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int m = 0; m < n_lR && !bad; ++m) {
      const int r = start + m;
      if (r >= n_rows_total) break;
      if (!isok_all[r] || R_IsNA(p1_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;                       // cannot have finished
      fperace::Key p;
      if (!fperace::rou_key_par(C->par_kind, p1_[r], p2_[r], p3_[r], B_[r],
                                A_[r], rou_bnd_row(C->bnd_kind, cols, r), p)) {
        bad = true; break;
      }
      const int g = fperace::cache_get(*C, p, tt, tt);
      const double lS = fperace::entry_log_S(C->e[g], tt);
      if (lS <= fperace::LOG_FLOOR) { bad = true; break; }
      if (lS < 0.0) logS += lS;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

double drou_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) return 0.0;
  const RouColIdx ci = rou_col_idx(C->par_kind);
  if (R_IsNA(par[ci.p1])) return 0.0;
  const double tt = t - par[ci.t0];
  if (!(tt > 0.0)) return 0.0;
  fperace::Key p;
  if (!fperace::rou_key_par(C->par_kind, par[ci.p1], par[ci.p2], par[ci.p3],
                            par[ci.B], par[ci.A],
                            rou_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, rou_scalar_horizon(tt), tt);
  const double lp = fperace::entry_log_pdf(C->e[g], tt);
  return (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
}

double prou_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) return 0.0;
  const RouColIdx ci = rou_col_idx(C->par_kind);
  if (R_IsNA(par[ci.p1])) return 0.0;
  const double tt = t - par[ci.t0];
  if (!(tt > 0.0)) return 0.0;
  fperace::Key p;
  if (!fperace::rou_key_par(C->par_kind, par[ci.p1], par[ci.p2], par[ci.p3],
                            par[ci.B], par[ci.A],
                            rou_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, rou_scalar_horizon(tt), tt);
  const double lS = fperace::entry_log_S(C->e[g], tt);
  if (lS >= 0.0) return 0.0;
  if (lS <= fperace::LOG_FLOOR) return 1.0;
  return -std::expm1(lS);            // 1 - S, accurate for S near 1
}

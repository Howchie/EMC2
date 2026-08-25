#include "utility_functions.h"
#include "race_contract.h"
#include "col_registry.h"
#include "model_ROUp.h"

void droup_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) {
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    }
    return;
  }
  roup_prepare_rows(*C, rt, cols, n_rows, isok);

  const double* t0_ = cols[roup_col_idx(C->par_kind).t0];
  if (C->roup_local) {
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      const int gs = C->row_group_s[i], gt = C->row_group_t[i];
      if (gs < 0 && gt < 0) {
        out[i] = raw_log_zero(min_ll, floor_raw);
        continue;
      }
      double lp = 0.0, ls = 0.0;
      roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0),
                      gs >= 0, gt >= 0, rt[i] - t0_[i], lp, ls);
      out[i] = raw_log_value(lp, min_ll, floor_raw);
    }
    return;
  }
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double lp = fperace::entry_log_pdf(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(lp, min_ll, floor_raw);
  }
}
void proup_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  roup_prepare_rows(*C, rt, cols, n_rows, isok);

  const double* t0_ = cols[roup_col_idx(C->par_kind).t0];
  if (C->roup_local) {
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      const int gs = C->row_group_s[i], gt = C->row_group_t[i];
      if (gs < 0 && gt < 0) {
        out[i] = 0.0;
        continue;
      }
      double lp = 0.0, ls = 0.0;
      roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0),
                      gs >= 0, gt >= 0, rt[i] - t0_[i], lp, ls);
      out[i] = raw_log_value(ls, min_ll, floor_raw);
    }
    return;
  }
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) {
      out[i] = 0.0;
      continue;
    }
    const double ls = fperace::entry_log_S(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(ls, min_ll, floor_raw);
  }
}

void roup_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int /*n_par*/,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out) {
  fperace::SolveCache* C = roup_cache(ctx_);
  const RoupColIdx ci = roup_col_idx(C->par_kind);
  const double* v_S_   = cols[ci.v_S];
  const double* transient_ = cols[ci.transient];
  const double* tau_S_ = cols[ci.tau_S];
  const double* tau_T_ = cols[ci.tau_T];
  const double* k_     = cols[ci.k];
  const double* B_     = cols[ci.B];
  const double* A_     = cols[ci.A];
  const double* t0_    = cols[ci.t0];
  const double* s_     = cols[ci.s];

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int m = 0; m < n_lR && !bad; ++m) {
      const int r = start + m;
      if (r >= n_rows_total) break;
      if (!isok_all[r] || R_IsNA(v_S_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;
      if (C->roup_local) {
        fperace::Key ks, kt; bool hs = false, ht = false;
        if (!roup_local_keys(C->par_kind, v_S_[r], transient_[r],
                             tau_S_[r], tau_T_[r], k_[r], B_[r], A_[r], s_[r],
                             roup_bnd_row(C->bnd_kind, cols, r),
                             ks, kt, hs, ht)) {
          bad = true;
          break;
        }
        const int gs = hs ? fperace::cache_get(*C, ks, tt) : -1;
        const int gt = ht ? fperace::cache_get(*C, kt, tt) : -1;
        double lp = 0.0, lS = 0.0;
        roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0),
                        hs, ht, tt, lp, lS);
        if (lS <= fperace::LOG_FLOOR) {
          bad = true;
          break;
        }
        if (lS < 0.0) logS += lS;
        continue;
      }
      fperace::Key p;
      if (!fperace::roup_key_par(C->par_kind, v_S_[r], transient_[r], tau_S_[r],
                                 tau_T_[r], k_[r], B_[r], A_[r], s_[r],
                                 roup_bnd_row(C->bnd_kind, cols, r), p)) {
        bad = true; break;
      }
      const int g = fperace::cache_get(*C, p, tt);
      const double lS = fperace::entry_log_S(C->e[g], tt);
      if (lS <= fperace::LOG_FLOOR) { bad = true; break; }
      if (lS < 0.0) logS += lS;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

double droup_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) return 0.0;
  const RoupColIdx ci = roup_col_idx(C->par_kind);
  if (R_IsNA(par[ci.v_S])) return 0.0;
  const double tt = t - par[ci.t0];
  if (!(tt > 0.0)) return 0.0;
  if (C->roup_local) {
    fperace::Key ks, kt; bool hs = false, ht = false;
    if (!roup_local_keys(C->par_kind, par[ci.v_S], par[ci.transient],
                         par[ci.tau_S], par[ci.tau_T], par[ci.k],
                         par[ci.B], par[ci.A], par[ci.s],
                         roup_bnd_par(C->bnd_kind, par),
                         ks, kt, hs, ht)) return 0.0;
    const int gs = hs ? fperace::cache_get(*C, ks, roup_scalar_horizon(tt)) : -1;
    const int gt = ht ? fperace::cache_get(*C, kt, roup_scalar_horizon(tt)) : -1;
    double lp = 0.0, lS = 0.0;
    roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0),
                    hs, ht, tt, lp, lS);
    return (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
  }
  fperace::Key p;
  if (!fperace::roup_key_par(C->par_kind, par[ci.v_S], par[ci.transient],
                             par[ci.tau_S], par[ci.tau_T], par[ci.k],
                             par[ci.B], par[ci.A], par[ci.s],
                             roup_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, roup_scalar_horizon(tt));
  const double lp = fperace::entry_log_pdf(C->e[g], tt);
  return (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
}

double proup_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) return 0.0;
  const RoupColIdx ci = roup_col_idx(C->par_kind);
  if (R_IsNA(par[ci.v_S])) return 0.0;
  const double tt = t - par[ci.t0];
  if (!(tt > 0.0)) return 0.0;
  if (C->roup_local) {
    fperace::Key ks, kt; bool hs = false, ht = false;
    if (!roup_local_keys(C->par_kind, par[ci.v_S], par[ci.transient],
                         par[ci.tau_S], par[ci.tau_T], par[ci.k],
                         par[ci.B], par[ci.A], par[ci.s],
                         roup_bnd_par(C->bnd_kind, par),
                         ks, kt, hs, ht)) return 0.0;
    const int gs = hs ? fperace::cache_get(*C, ks, roup_scalar_horizon(tt)) : -1;
    const int gt = ht ? fperace::cache_get(*C, kt, roup_scalar_horizon(tt)) : -1;
    double lp = 0.0, lS = 0.0;
    roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0),
                    hs, ht, tt, lp, lS);
    if (lS >= 0.0) return 0.0;
    if (lS <= fperace::LOG_FLOOR) return 1.0;
    return -std::expm1(lS);
  }
  fperace::Key p;
  if (!fperace::roup_key_par(C->par_kind, par[ci.v_S], par[ci.transient],
                             par[ci.tau_S], par[ci.tau_T], par[ci.k],
                             par[ci.B], par[ci.A], par[ci.s],
                             roup_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, roup_scalar_horizon(tt));
  const double lS = fperace::entry_log_S(C->e[g], tt);
  if (lS >= 0.0) return 0.0;
  if (lS <= fperace::LOG_FLOOR) return 1.0;
  return -std::expm1(lS);
}

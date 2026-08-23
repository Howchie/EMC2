#include <RcppArmadillo.h>
#include "model_lnr.h"
#include "race_contract.h"
#include "col_registry.h"
#include "wald_functions.h"

// LNR race-model kernels.  These bodies are moved verbatim from utils.h; they
// keep the same NA/Inf/floor behaviour.  utils.h is intentionally NOT included
// here so this translation unit only depends on the narrow race contract and
// the lognormal helpers in wald_functions.h.

double dlnr_scalar(double t, const double* par, void* /*ctx_*/) {
  const double m  = par[0];
  const double s  = par[1];
  const double t0 = par[2];

  if (R_IsNA(m)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return dlnorm_std(tt, m, s, false);
}

double plnr_scalar(double t, const double* par, void* /*ctx_*/) {
  const double m  = par[0];
  const double s  = par[1];
  const double t0 = par[2];

  if (R_IsNA(m)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return plnorm_std(tt, m, s, true, false);
}

void dlnr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* m_  = cols[emc2col::lnr::m];
  const double* s_  = cols[emc2col::lnr::s];
  const double* t0_ = cols[emc2col::lnr::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(m_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    // Direct log density: no natural round trip, so far-tail log values stay
    // finite instead of collapsing to min_ll once exp() underflows.
    out[i] = raw_log_value(dlnorm_std(tt, m_[i], s_[i], true), min_ll, floor_raw);
  }
}

void plnr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* /*ctx_*/) {
  const double* m_  = cols[emc2col::lnr::m];
  const double* s_  = cols[emc2col::lnr::s];
  const double* t0_ = cols[emc2col::lnr::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(m_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double logS = lnorm_log_surv_std(tt, m_[i], s_[i]);
    if (!R_FINITE(logS)) { out[i] = 0.0; continue; }
    out[i] = logS;
  }
}

void lnr_logS_at_t(double t, const double* const* cols,
                   int n_rows_total, int n_lR, int /*n_par*/,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* m_  = cols[emc2col::lnr::m];
  const double* s_  = cols[emc2col::lnr::s];
  const double* t0_ = cols[emc2col::lnr::t0];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(m_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;
      const double logSk = lnorm_log_surv_std(tt, m_[r], s_[r]);
      if (!R_FINITE(logSk)) { bad = true; break; }
      logS += logSk;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

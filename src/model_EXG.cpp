#include <RcppArmadillo.h>
#include "model_EXG.h"
#include "race_contract.h"
#include "col_registry.h"
#include "utility_functions.h"
#include "exgaussian_functions.h"
#include <cmath>

void drexg_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* mu_    = cols[emc2col::rexg::mu];
  const double* sigma_ = cols[emc2col::rexg::sigma];
  const double* tau_   = cols[emc2col::rexg::tau];
  const double* t0_    = cols[emc2col::rexg::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(mu_[i]) || R_IsNA(sigma_[i]) || R_IsNA(tau_[i]) || R_IsNA(t0_[i]) ||
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
    out[i] = raw_log_value(dtexg(tt, mu_[i], sigma_[i], tau_[i], 0.0, R_PosInf, true),
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(mu_[i]) || R_IsNA(sigma_[i]) || R_IsNA(tau_[i]) || R_IsNA(t0_[i]) ||
        !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0 ||
        !R_FINITE(t0_[i]) || t0_[i] < 0.0) {
      out[i] = 0.0;
      continue;
    }
    const double log_surv = (rt[i] - t0_[i] <= 0.0)
      ? 0.0 : ptexg(rt[i] - t0_[i], mu_[i], sigma_[i], tau_[i], 0.0, R_PosInf, false, true);
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
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(mu_[r]) || R_IsNA(sigma_[r]) || R_IsNA(tau_[r]) ||
          R_IsNA(t0_[r]) || !R_FINITE(t0_[r]) ||
          sigma_[r] <= 0.0 || tau_[r] <= 0.0 || t0_[r] < 0.0) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;
      const double log_surv = ptexg(tt, mu_[r], sigma_[r], tau_[r], 0.0, R_PosInf, false, true);
      if (!R_FINITE(log_surv)) {
        if (log_surv == R_NegInf) { bad = true; break; }
        continue;
      }
      logS += std::fmin(log_surv, 0.0);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

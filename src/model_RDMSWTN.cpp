#include <RcppArmadillo.h>
#include "model_RDMSWTN.h"
#include "model_RDM.h"
#include "col_registry.h"

// ============================================================
// RDMSWTN_TT adapters
// Column layout: v=0, B=1, A=2, t0=3, s=4, sv=5, tau=6
// ============================================================

double drdmswtn_tt_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::rdmswtn_tt::v])) return 0.0;
  return drdmswtn_tt(
      t, par[emc2col::rdmswtn_tt::v],
      par[emc2col::rdmswtn_tt::B] + par[emc2col::rdmswtn_tt::A],
      par[emc2col::rdmswtn_tt::A], par[emc2col::rdmswtn_tt::s],
      par[emc2col::rdmswtn_tt::t0], par[emc2col::rdmswtn_tt::sv],
      par[emc2col::rdmswtn_tt::tau], false,
      ctx ? ctx->use_posdrift : true);
}

double prdmswtn_tt_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::rdmswtn_tt::v])) return 0.0;
  return prdmswtn_tt(
      t, par[emc2col::rdmswtn_tt::v],
      par[emc2col::rdmswtn_tt::B] + par[emc2col::rdmswtn_tt::A],
      par[emc2col::rdmswtn_tt::A], par[emc2col::rdmswtn_tt::s],
      par[emc2col::rdmswtn_tt::t0], par[emc2col::rdmswtn_tt::sv],
      par[emc2col::rdmswtn_tt::tau], false,
      ctx ? ctx->use_posdrift : true);
}

void drdmswtn_tt_raw(
    const double* rt, const double* const* cols, int n_rows,
    const int* mask, const int* isok, double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_ = cols[emc2col::rdmswtn_tt::v];
  const double* B_ = cols[emc2col::rdmswtn_tt::B];
  const double* A_ = cols[emc2col::rdmswtn_tt::A];
  const double* t0_ = cols[emc2col::rdmswtn_tt::t0];
  const double* s_ = cols[emc2col::rdmswtn_tt::s];
  const double* sv_ = cols[emc2col::rdmswtn_tt::sv];
  const double* tau_ = cols[emc2col::rdmswtn_tt::tau];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(v_[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double log_pdf = drdmswtn_tt(
        rt[i], v_[i], B_[i] + A_[i], A_[i], s_[i], t0_[i], sv_[i],
        tau_[i], true, pd);
    out[i] = raw_log_value(log_pdf, min_ll, floor_raw);
  }
}

void prdmswtn_tt_raw(
    const double* rt, const double* const* cols, int n_rows,
    const int* mask, const int* isok, double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_ = cols[emc2col::rdmswtn_tt::v];
  const double* B_ = cols[emc2col::rdmswtn_tt::B];
  const double* A_ = cols[emc2col::rdmswtn_tt::A];
  const double* t0_ = cols[emc2col::rdmswtn_tt::t0];
  const double* s_ = cols[emc2col::rdmswtn_tt::s];
  const double* sv_ = cols[emc2col::rdmswtn_tt::sv];
  const double* tau_ = cols[emc2col::rdmswtn_tt::tau];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(v_[i])) {
      out[i] = 0.0;
      continue;
    }
    const double log_cdf = prdmswtn_tt(
        rt[i], v_[i], B_[i] + A_[i], A_[i], s_[i], t0_[i], sv_[i],
        tau_[i], true, pd);
    if (log_cdf == R_NegInf) {
      out[i] = 0.0;
    } else if (ISNAN(log_cdf) || log_cdf >= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
    } else {
      out[i] = log1m_exp(log_cdf);
    }
  }
}

void rdmswtn_tt_logS_at_t(
    double t, const double* const* cols, int /*n_rows_total*/, int n_lR,
    int /*n_par*/, const int* trunc_mask, int n_unique_trials,
    const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_ = cols[emc2col::rdmswtn_tt::v];
  const double* B_ = cols[emc2col::rdmswtn_tt::B];
  const double* A_ = cols[emc2col::rdmswtn_tt::A];
  const double* t0_ = cols[emc2col::rdmswtn_tt::t0];
  const double* s_ = cols[emc2col::rdmswtn_tt::s];
  const double* sv_ = cols[emc2col::rdmswtn_tt::sv];
  const double* tau_ = cols[emc2col::rdmswtn_tt::tau];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(v_[r])) {
        bad = true;
        break;
      }
      const double log_cdf = prdmswtn_tt(
          t, v_[r], B_[r] + A_[r], A_[r], s_[r], t0_[r], sv_[r],
          tau_[r], true, pd);
      if (ISNAN(log_cdf) || log_cdf >= 0.0) {
        bad = true;
        break;
      }
      if (R_FINITE(log_cdf)) logS += log1m_exp(log_cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

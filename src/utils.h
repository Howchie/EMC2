#ifndef utils_h
#define utils_h

#include <RcppArmadillo.h>
#include "exgaussian_functions.h"
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <limits>
#include "col_registry.h"
#include "fpe_race.h"
#include "fpe_bou.h"
#include "model_RLF.h"
#include "utility_functions.h"
#include "model_RDM.h"
#include "model_RDM_adapters.h"
#include "model_RDMSWTN.h"
#include "model_GBM.h"
#include "model_LBA.h"
#include "model_BAwL.h"
#include "model_BAwD.h"
#include "model_BAwF.h"
#include "model_BAwR.h"
#include "model_BTAwL.h"
#include "race_contract.h"

using namespace Rcpp;

#include "timer_helpers.h"

// Included here, not at the top: the ROU kernels need ContextForRaceModels and
// the raw_log_* helpers above.  (The cache TYPE they store in the context comes
// from fpe_race.h, which has no such dependency and is included at the top.)
#include "model_ROU.h"
#include "model_ROUp.h"
#include "model_GOM.h"
#include "model_RLF_kernels.h"
// Likewise: the BOU primitives need ContextForDDMModels defined above.
#include "model_BOU.h"
// FRQ needs nothing from this file -- it is pure Rmath -- but it carries
// [[Rcpp::export]] entry points, so it must be seen by exactly one translation
// unit and this header is included by exactly one (particle_ll.cpp).
#include "model_FRQ.h"
#include "model_PCOUNTER.h"


// BAwF shares BAwD's launch and rho context fields: the two models never
// coexist in one adapter, and duplicating the fields would create a second
// place for the "0 = normal, 1 = lognormal" contract to drift.
inline int bawf_launch_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_launch : BAWF_LAUNCH_LOGNORMAL;
}
inline double bawf_rho_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_rho : R_PosInf;
}

inline double dbawf_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawf::v])) return 0.0;
  const double tt = t - par[emc2col::bawf::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  return bawf_pdf_scalar_natural(
    tt, par[emc2col::bawf::A],
    par[emc2col::bawf::B] + par[emc2col::bawf::A],
    par[emc2col::bawf::v], par[emc2col::bawf::sv],
    par[emc2col::bawf::k],
    bawf_launch_of(ctx), ctx ? ctx->use_posdrift : true, bawf_rho_of(ctx));
}

inline double pbawf_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawf::v])) return 0.0;
  const double tt = t - par[emc2col::bawf::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf is handled inside the kernel and returns F_max, not 1.
  return bawf_cdf_scalar_natural(
    tt, par[emc2col::bawf::A],
    par[emc2col::bawf::B] + par[emc2col::bawf::A],
    par[emc2col::bawf::v], par[emc2col::bawf::sv],
    par[emc2col::bawf::k],
    bawf_launch_of(ctx), ctx ? ctx->use_posdrift : true, bawf_rho_of(ctx));
}

inline void dbawf_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawf_launch_of(ctx);
  const double rho = bawf_rho_of(ctx);
  const double* p1_ = cols[emc2col::bawf::v];
  const double* p2_ = cols[emc2col::bawf::sv];
  const double* B_  = cols[emc2col::bawf::B];
  const double* A_  = cols[emc2col::bawf::A];
  const double* t0_ = cols[emc2col::bawf::t0];
  const double* k_  = cols[emc2col::bawf::k];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double log_pdf = bawf_log_pdf(tt, A_[i], B_[i] + A_[i], p1_[i],
                                        p2_[i], k_[i], launch, pd, rho);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbawf_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawf_launch_of(ctx);
  const double rho = bawf_rho_of(ctx);
  const double* p1_ = cols[emc2col::bawf::v];
  const double* p2_ = cols[emc2col::bawf::sv];
  const double* B_  = cols[emc2col::bawf::B];
  const double* A_  = cols[emc2col::bawf::A];
  const double* t0_ = cols[emc2col::bawf::t0];
  const double* k_  = cols[emc2col::bawf::k];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    double cdf = 0.0;
    if (ba_natural_cdf_bawf(tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i], k_[i],
                            launch, pd, rho, BAWF_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = bawf_log_surv(tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i],
                                      k_[i], launch, pd, rho, BAWF_DENOM_FLOOR);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

inline void bawf_logS_at_t(double t, const double* const* cols,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawf_launch_of(ctx);
  const double rho = bawf_rho_of(ctx);
  const double* p1_ = cols[emc2col::bawf::v];
  const double* p2_ = cols[emc2col::bawf::sv];
  const double* B_  = cols[emc2col::bawf::B];
  const double* A_  = cols[emc2col::bawf::A];
  const double* t0_ = cols[emc2col::bawf::t0];
  const double* k_  = cols[emc2col::bawf::k];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(p1_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      double cdf = 0.0;
      if (ba_natural_cdf_bawf(tt, A_[r], B_[r] + A_[r], p1_[r], p2_[r], k_[r],
                              launch, pd, rho, BAWF_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
        if (cdf > 0.0) logS += std::log1p(-cdf);
      } else {
        const double ls = bawf_log_surv(tt, A_[r], B_[r] + A_[r], p1_[r], p2_[r],
                                        k_[r], launch, pd, rho, BAWF_DENOM_FLOOR);
        if (!(ls > R_NegInf) || ISNAN(ls)) { bad = true; break; }
        logS += ls;
      }
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// BAwR shares BAwD's launch context field, for the same reason BAwF does: the
// models never coexist in one adapter and the "0 = normal, 1 = lognormal"
// contract should have exactly one definition.  BAwR needs no rho: its decay
// shape is the sampled exponent `p`, not a fixed kernel index.
inline int bawr_launch_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_launch : BAWR_LAUNCH_LOGNORMAL;
}

// BAwR always uses the endpoint chart: the sampled clearance slot is Tmax,
// while the mechanistic kernel continues to receive kappa.  Funnel every read
// through this inverse so scalar, raw, and truncation paths cannot disagree.
inline double bawr_clear_to_kappa(double clear, double b, double pw) {
  return bawr_kappa_from_tmax(b, pw, clear);
}

inline double dbawr_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawr::v])) return 0.0;
  const double tt = t - par[emc2col::bawr::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  return bawr_pdf_scalar_natural(
    tt, par[emc2col::bawr::A],
    par[emc2col::bawr::B] + par[emc2col::bawr::A],
    par[emc2col::bawr::v], par[emc2col::bawr::sv],
    bawr_clear_to_kappa(par[emc2col::bawr::clear],
                        par[emc2col::bawr::B] + par[emc2col::bawr::A],
                        par[emc2col::bawr::p]), par[emc2col::bawr::p],
    bawr_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline double pbawr_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawr::v])) return 0.0;
  const double tt = t - par[emc2col::bawr::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf is handled inside the kernel and returns F_max, not 1.
  return bawr_cdf_scalar_natural(
    tt, par[emc2col::bawr::A],
    par[emc2col::bawr::B] + par[emc2col::bawr::A],
    par[emc2col::bawr::v], par[emc2col::bawr::sv],
    bawr_clear_to_kappa(par[emc2col::bawr::clear],
                        par[emc2col::bawr::B] + par[emc2col::bawr::A],
                        par[emc2col::bawr::p]), par[emc2col::bawr::p],
    bawr_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline void dbawr_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawr_launch_of(ctx);
  const double* p1_ = cols[emc2col::bawr::v];
  const double* p2_ = cols[emc2col::bawr::sv];
  const double* B_  = cols[emc2col::bawr::B];
  const double* A_  = cols[emc2col::bawr::A];
  const double* t0_ = cols[emc2col::bawr::t0];
  const double* clear_ = cols[emc2col::bawr::clear];
  const double* pw_ = cols[emc2col::bawr::p];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double log_pdf = bawr_log_pdf(tt, A_[i], B_[i] + A_[i], p1_[i],
                                        p2_[i], bawr_clear_to_kappa(clear_[i], B_[i] + A_[i], pw_[i]), pw_[i], launch, pd);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbawr_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawr_launch_of(ctx);
  const double* p1_ = cols[emc2col::bawr::v];
  const double* p2_ = cols[emc2col::bawr::sv];
  const double* B_  = cols[emc2col::bawr::B];
  const double* A_  = cols[emc2col::bawr::A];
  const double* t0_ = cols[emc2col::bawr::t0];
  const double* clear_ = cols[emc2col::bawr::clear];
  const double* pw_ = cols[emc2col::bawr::p];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const double kap = bawr_clear_to_kappa(clear_[i], B_[i] + A_[i], pw_[i]);
    double cdf = 0.0;
    if (ba_natural_cdf_bawr(tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i], kap, pw_[i],
                            launch, pd, BAWR_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = bawr_log_surv(tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i],
                                      kap, pw_[i], launch, pd, BAWR_DENOM_FLOOR);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

inline void bawr_logS_at_t(double t, const double* const* cols,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawr_launch_of(ctx);
  const double* p1_ = cols[emc2col::bawr::v];
  const double* p2_ = cols[emc2col::bawr::sv];
  const double* B_  = cols[emc2col::bawr::B];
  const double* A_  = cols[emc2col::bawr::A];
  const double* t0_ = cols[emc2col::bawr::t0];
  const double* clear_ = cols[emc2col::bawr::clear];
  const double* pw_ = cols[emc2col::bawr::p];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(p1_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      const double kap = bawr_clear_to_kappa(clear_[r], B_[r] + A_[r], pw_[r]);
      double cdf = 0.0;
      if (ba_natural_cdf_bawr(tt, A_[r], B_[r] + A_[r], p1_[r], p2_[r], kap, pw_[r],
                              launch, pd, BAWR_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
        if (cdf > 0.0) logS += std::log1p(-cdf);
      } else {
        const double ls = bawr_log_surv(tt, A_[r], B_[r] + A_[r], p1_[r], p2_[r],
                                        kap, pw_[r], launch, pd, BAWR_DENOM_FLOOR);
        if (!(ls > R_NegInf) || ISNAN(ls)) { bad = true; break; }
        logS += ls;
      }
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
  (void)n_rows_total;
}




// Helper to safely get a column from a DataFrame with a default value if missing
// Also fills NA values with the default for backward compatibility.
inline Rcpp::NumericVector get_col_with_default(const Rcpp::DataFrame& df, const std::string& name, double default_val) {
  if (df.containsElementNamed(name.c_str())) {
    Rcpp::NumericVector col = df[name];
    // Check for NAs and replace if necessary
    bool has_na = false;
    for (int i = 0; i < col.size(); ++i) {
      if (Rcpp::NumericVector::is_na(col[i])) {
        has_na = true;
        break;
      }
    }
    if (has_na) {
      Rcpp::NumericVector res = Rcpp::clone(col);
      for (int i = 0; i < res.size(); ++i) {
        if (Rcpp::NumericVector::is_na(res[i])) {
          res[i] = default_val;
        }
      }
      return res;
    }
    return col;
  }
  return Rcpp::NumericVector(df.nrow(), default_val);
}

#endif

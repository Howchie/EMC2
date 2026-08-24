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

// ============================================================
// BAwDp (proportional-clearance drive clock)
// ============================================================

struct BawDpClock {
  bool ok = false;
  bool frozen = false;
  double m = 0.0;
  double dm = 0.0;
};

// Internal clock m(u) = (1-exp(-ku))/k - lambda*u.  For 0 < lambda < 1,
// the clock freezes at u* = log(1/lambda)/k and m_max is finite.  This is an
// elementary change of variable; the downstream LBA evaluation is always at
// k = 0, so no numerical integration is introduced.
inline BawDpClock bawdp_clock(double u, double k, double lambda) {
  BawDpClock out;
  if (!(u > 0.0) || !(k >= 0.0) || !(lambda >= 0.0) || !(lambda < 1.0) ||
      !emc2_isfinite(k) || !emc2_isfinite(lambda)) return out;

  if (k <= BAWL_K_EPS) {
    out.ok = true;
    out.frozen = false;
    out.m = (u == R_PosInf) ? R_PosInf : (1.0 - lambda) * u;
    out.dm = 1.0 - lambda;
    return out;
  }

  // The clock ceiling: 1/k at lambda = 0, approached only as u -> Inf, and
  // m(u*) for lambda > 0, attained at the universal freeze time.
  const double log_lambda = (lambda > 0.0) ? std::log(lambda) : R_NegInf;
  const double u_star = (lambda > 0.0) ? -log_lambda / k : R_PosInf;
  const double m_max = (lambda > 0.0)
    ? (-std::expm1(log_lambda) + lambda * log_lambda) / k
    : 1.0 / k;
  if (!emc2_isfinite(m_max) || !(m_max > 0.0)) return out;
  auto freeze = [&]() {
    out.ok = true; out.frozen = true; out.m = m_max; out.dm = 0.0;
    return out;
  };
  if (u == R_PosInf || u >= u_star) return freeze();

  const double ku = k * u;
  const double e = std::exp(-ku);
  const double m = -std::expm1(-ku) / k - lambda * u;
  const double dm = e - lambda;
  // exp(-k u) underflows to zero for k u > ~745.  The clock has then reached
  // m_max to machine precision and the density below it is under 1e-320, so
  // freezing is the correct limit.  Reporting failure instead would return a
  // zero CDF where the defective ceiling belongs -- a survivor of one for
  // every loser in the race.
  if (!(dm > 0.0)) return freeze();
  if (!emc2_isfinite(m) || !(m > 0.0)) return out;
  out.ok = true;
  out.m = m;
  out.dm = dm;
  return out;
}

inline double dbawdp_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = ctx ? ctx->bawl_launch : BAWL_LAUNCH_LOGNORMAL;
  const double tt = t - par[emc2col::bawdp::t0];
  if (t <= 0.0 || !(tt > 0.0)) return 0.0;
  const BawDpClock g = bawdp_clock(tt, par[emc2col::bawdp::k],
                                   par[emc2col::bawdp::lambda]);
  if (!g.ok || g.frozen) return 0.0;
  const double lp = log_ba_pdf_launch(g.m, par[emc2col::bawdp::A],
                                      par[emc2col::bawdp::B] + par[emc2col::bawdp::A],
                                      par[emc2col::bawdp::v], par[emc2col::bawdp::sv],
                                      0.0, ctx ? ctx->use_posdrift : true,
                                      BAWL_DENOM_FLOOR, launch);
  return (lp > R_NegInf && emc2_isfinite(lp)) ? std::exp(lp + std::log(g.dm)) : 0.0;
}

inline double pbawdp_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = ctx ? ctx->bawl_launch : BAWL_LAUNCH_LOGNORMAL;
  const double tt = t - par[emc2col::bawdp::t0];
  if (t <= 0.0 || !(tt > 0.0)) return 0.0;
  const BawDpClock g = bawdp_clock(tt, par[emc2col::bawdp::k],
                                   par[emc2col::bawdp::lambda]);
  if (!g.ok) return 0.0;
  const double lp = log_ba_cdf_launch(g.m, par[emc2col::bawdp::A],
                                      par[emc2col::bawdp::B] + par[emc2col::bawdp::A],
                                      par[emc2col::bawdp::v], par[emc2col::bawdp::sv],
                                      0.0, ctx ? ctx->use_posdrift : true,
                                      BAWL_DENOM_FLOOR, launch);
  return (lp > R_NegInf) ? std::fmin(std::exp(lp), 1.0) : 0.0;
}

inline void dbawdp_raw(const double* rt, const double* const* cols, int n_rows,
                       const int* mask, const int* isok, double* out,
                       double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = ctx ? ctx->bawl_launch : BAWL_LAUNCH_LOGNORMAL;
  const double* p1 = cols[emc2col::bawdp::v];
  const double* p2 = cols[emc2col::bawdp::sv];
  const double* B = cols[emc2col::bawdp::B];
  const double* A = cols[emc2col::bawdp::A];
  const double* t0 = cols[emc2col::bawdp::t0];
  const double* k = cols[emc2col::bawdp::k];
  const double* lam = cols[emc2col::bawdp::lambda];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double tt = rt[i] - t0[i];
    const BawDpClock g = bawdp_clock(tt, k[i], lam[i]);
    if (!isok[i] || rt[i] <= 0.0 || !(tt > 0.0) || !g.ok || g.frozen) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double lp = log_ba_pdf_launch(g.m, A[i], B[i] + A[i], p1[i], p2[i],
                                        0.0, pd, BAWL_DENOM_FLOOR, launch);
    const double out_lp = (lp > R_NegInf && emc2_isfinite(lp))
      ? lp + std::log(g.dm) : R_NegInf;
    out[i] = (out_lp > R_NegInf && emc2_isfinite(out_lp))
      ? raw_log_value(out_lp, min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbawdp_raw(const double* rt, const double* const* cols, int n_rows,
                       const int* mask, const int* isok, double* out,
                       double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = ctx ? ctx->bawl_launch : BAWL_LAUNCH_LOGNORMAL;
  const double* p1 = cols[emc2col::bawdp::v];
  const double* p2 = cols[emc2col::bawdp::sv];
  const double* B = cols[emc2col::bawdp::B];
  const double* A = cols[emc2col::bawdp::A];
  const double* t0 = cols[emc2col::bawdp::t0];
  const double* k = cols[emc2col::bawdp::k];
  const double* lam = cols[emc2col::bawdp::lambda];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double tt = rt[i] - t0[i];
    if (!isok[i] || rt[i] <= 0.0 || !(tt > 0.0)) { out[i] = 0.0; continue; }
    const BawDpClock g = bawdp_clock(tt, k[i], lam[i]);
    if (!g.ok) { out[i] = 0.0; continue; }
    const double b_i = B[i] + A[i];
    double cdf = 0.0;
    if (ba_natural_cdf_launch(g.m, A[i], b_i, p1[i], p2[i], 0.0, pd,
                              BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, launch)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = log_ba_surv_launch(g.m, A[i], b_i, p1[i], p2[i],
                                           0.0, pd, BAWL_DENOM_FLOOR, launch);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

inline void bawdp_logS_at_t(double t, const double* const* cols,
                            int /*n_rows_total*/, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = ctx ? ctx->bawl_launch : BAWL_LAUNCH_LOGNORMAL;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double logS = 0.0; bool bad = false;
    for (int a = 0; a < n_lR && !bad; ++a) {
      const int r = j * n_lR + a;
      if (!isok_all[r]) { bad = true; break; }
      const double tt = t - cols[emc2col::bawdp::t0][r];
      if (!(tt > 0.0)) continue;
      const BawDpClock g = bawdp_clock(tt, cols[emc2col::bawdp::k][r],
                                       cols[emc2col::bawdp::lambda][r]);
      if (!g.ok) { bad = true; break; }
      const double b_r = cols[emc2col::bawdp::B][r] + cols[emc2col::bawdp::A][r];
      double cdf = 0.0;
      if (ba_natural_cdf_launch(g.m, cols[emc2col::bawdp::A][r], b_r,
                                cols[emc2col::bawdp::v][r], cols[emc2col::bawdp::sv][r],
                                0.0, pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, launch)) {
        if (cdf > 0.0) logS += std::log1p(-cdf);
      } else {
        const double ls = log_ba_surv_launch(
          g.m, cols[emc2col::bawdp::A][r], b_r,
          cols[emc2col::bawdp::v][r], cols[emc2col::bawdp::sv][r], 0.0, pd,
          BAWL_DENOM_FLOOR, launch);
        if (!(ls > R_NegInf)) { bad = true; break; }
        logS += ls;
      }
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

inline double bawdp_pdf_norm(double t, double A, double b, double p1, double p2,
                             double k, double lambda, int launch, bool posdrift,
                             bool log_out) {
  const BawDpClock g = bawdp_clock(t, k, lambda);
  if (!g.ok || g.frozen) return log_out ? R_NegInf : 0.0;
  const double lp = log_ba_pdf_launch(g.m, A, b, p1, p2, 0.0, posdrift,
                                      BAWL_DENOM_FLOOR, launch);
  if (!(lp > R_NegInf) || !emc2_isfinite(lp)) return log_out ? R_NegInf : 0.0;
  const double out = lp + std::log(g.dm);
  return log_out ? out : std::exp(out);
}

inline double bawdp_cdf_norm(double t, double A, double b, double p1, double p2,
                             double k, double lambda, int launch, bool posdrift,
                             bool log_out) {
  const BawDpClock g = bawdp_clock(t, k, lambda);
  if (!g.ok) return log_out ? R_NegInf : 0.0;
  const double lp = log_ba_cdf_launch(g.m, A, b, p1, p2, 0.0, posdrift,
                                      BAWL_DENOM_FLOOR, launch);
  if (log_out) return lp;
  return (lp > R_NegInf) ? std::fmin(std::exp(lp), 1.0) : 0.0;
}

// [[Rcpp::export]]
NumericVector dbawdp(NumericVector t, NumericVector A, NumericVector b,
                     NumericVector p1, NumericVector p2, NumericVector k,
                     NumericVector lambda, int launch = 1, bool posdrift = true,
                     bool log_out = false) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawdp_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                            pick(p2, i), pick(k, i), pick(lambda, i),
                            launch, posdrift, log_out);
  return out;
}

// [[Rcpp::export]]
NumericVector pbawdp(NumericVector t, NumericVector A, NumericVector b,
                     NumericVector p1, NumericVector p2, NumericVector k,
                     NumericVector lambda, int launch = 1, bool posdrift = true,
                     bool log_out = false) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawdp_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                            pick(p2, i), pick(k, i), pick(lambda, i),
                            launch, posdrift, log_out);
  return out;
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

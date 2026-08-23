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
#include "model_LBA.h"
#include "race_contract.h"

using namespace Rcpp;

inline bool raw_floor_log_lik(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx == nullptr || ctx->floor_raw_log_lik;
}

inline double erlang_omega_for_shape(int kill_shape, const double* par = nullptr,
                                     int omega_index = -1) {
  if (kill_shape <= 1) return 1.0;
  if (kill_shape == 2) return 0.0;
  if (par == nullptr || omega_index < 0) return 1.0;
  return std::fmax(0.0, std::fmin(1.0, par[omega_index]));
}

// Convert a timer mean to the Erlang rate parameter used by erlang_log_surv /
// erlang_log_pdf.  The C++ particle likelihood receives raw timer means (mG,
// mK) directly from the sampled parameter space; Ttransform is NOT applied on
// the C++ path.  This function performs the shape-dependent conversion:
//   Erlang-1 (exponential): rate = 1 / mean
//   Erlang-2:               rate = 2 / mean   (so that E[T] = 2/rate = mean)
//   EMIX (shape 3):         rate = 1 / mean   (each component is rescaled
//                                               inside erlang_log_surv for n=3)
inline double erlang_lambda_from_mean(double mean, int kill_shape) {
  if (!(mean > 0.0) || !emc2_isfinite(mean)) return 0.0;
  return ((kill_shape == 2) ? 2.0 : 1.0) / mean;
}

inline double raw_log_zero(double min_ll, bool floor_raw) {
  return floor_raw ? min_ll : R_NegInf;
}

inline double raw_log_value(double log_x, double min_ll, bool floor_raw) {
  if (!R_FINITE(log_x)) return raw_log_zero(min_ll, floor_raw);
  return floor_raw ? ((log_x > min_ll) ? log_x : min_ll) : log_x;
}

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

struct TimedLambdaDispatch {
  double lambda_g;
  double lambda_k;
  bool guess;
  bool use_combo;
};

inline TimedLambdaDispatch timed_lambda_dispatch(const ContextForRaceModels* ctx,
                                                 double lambda_g,
                                                 double lambda_k) {
  constexpr double kLamEps = 1e-12;
  const bool local_guess_only = ctx && ctx->is_local_guess;
  const bool local_kill_guess = ctx && ctx->is_local_kill_guess;
  const bool has_guess = lambda_g > kLamEps;
  const bool has_kill  = lambda_k > kLamEps;

  TimedLambdaDispatch out{0.0, 0.0, false, local_kill_guess && has_guess && has_kill};
  if (out.use_combo) return out;

  if (local_guess_only || (local_kill_guess && has_guess)) {
    out.lambda_g = lambda_g;
    out.guess = true;
  } else if (has_kill) {
    out.lambda_k = lambda_k;
  }
  return out;
}

// Column layout: v=0, B=1, A=2, t0=3, s=4
inline double drdm_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[3];
  if (tt <= 0.0) return 0.0;
  const double inv_s = 1.0 / par[4];
  return digt_impl(tt, par[1] * inv_s + 0.5 * par[2] * inv_s,
                       par[0] * inv_s, 0.5 * par[2] * inv_s);
}

inline double prdm_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[3];
  if (tt <= 0.0) return 0.0;
  const double inv_s = 1.0 / par[4];
  return pigt_impl(tt, par[1] * inv_s + 0.5 * par[2] * inv_s,
                       par[0] * inv_s, 0.5 * par[2] * inv_s);
}

// GBM: column layout v=0, B=1, A=2, t0=3, s=4, mG=5 (guess-clock mean), mK=6 (kill-clock mean)
inline double drdmgbm_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double lg = (ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[5], ks) : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? erlang_lambda_from_mean(par[6], ks) : 0.0;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
  if (dispatch.use_combo) {
    return dgbm_local_combo(t, par[0], 1.0 + par[1] + par[2],
                            par[2], par[4], t0_val, lg, lk, false, ks, omega);
  }
  return dgbm(t,
              par[0],
              1.0 + par[1] + par[2],
              par[2],
              par[4],
              t0_val,
              dispatch.lambda_g, dispatch.lambda_k,
              false, ks, dispatch.guess, omega);
}

inline double prdmgbm_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double lg = (ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[5], ks) : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? erlang_lambda_from_mean(par[6], ks) : 0.0;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
  if (dispatch.use_combo) {
    return pgbm_local_combo(t, par[0], 1.0 + par[1] + par[2],
                            par[2], par[4], t0_val, lg, lk, false, ks, omega);
  }
  return pgbm(t,
              par[0],
              1.0 + par[1] + par[2],
              par[2],
              par[4],
              t0_val,
              dispatch.lambda_g, dispatch.lambda_k,
              false, ks, dispatch.guess, omega);
}

inline double dlnr_scalar(double t, const double* par, void* /*ctx_*/) {
  const double m  = par[0];
  const double s  = par[1];
  const double t0 = par[2];

  if (R_IsNA(m)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return dlnorm_std(tt, m, s, false);
}

inline double plnr_scalar(double t, const double* par, void* /*ctx_*/) {
  const double m  = par[0];
  const double s  = par[1];
  const double t0 = par[2];

  if (R_IsNA(m)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return plnorm_std(tt, m, s, true, false);
}

// Column order per src/col_registry.h.
inline void drdm_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = cols[emc2col::rdm::v];
  const double* B_  = cols[emc2col::rdm::B];
  const double* A_  = cols[emc2col::rdm::A];
  const double* t0_ = cols[emc2col::rdm::t0];
  const double* s_  = cols[emc2col::rdm::s];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double inv_s = 1.0 / s_[i];
    const double pdf = digt_impl(tt, (B_[i] + 0.5 * A_[i]) * inv_s,
                                     v_[i] * inv_s, 0.5 * A_[i] * inv_s);
    out[i] = (pdf > 0.0 && emc2_isfinite(pdf)) ? raw_log_value(std::log(pdf), min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void prdm_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = cols[emc2col::rdm::v];
  const double* B_  = cols[emc2col::rdm::B];
  const double* A_  = cols[emc2col::rdm::A];
  const double* t0_ = cols[emc2col::rdm::t0];
  const double* s_  = cols[emc2col::rdm::s];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double inv_s = 1.0 / s_[i];
    const double cdf = pigt_impl(tt, (B_[i] + 0.5 * A_[i]) * inv_s,
                                     v_[i] * inv_s, 0.5 * A_[i] * inv_s);
    if (cdf >= 1.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

// Truncation survivor helpers (log-survivor at a scalar T, used for normalization)
inline void rdm_logS_at_t(double t, const double* const* cols,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* v_  = cols[emc2col::rdm::v];
  const double* B_  = cols[emc2col::rdm::B];
  const double* A_  = cols[emc2col::rdm::A];
  const double* t0_ = cols[emc2col::rdm::t0];
  const double* s_  = cols[emc2col::rdm::s];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;
      const double inv_s = 1.0 / s_[r];
      const double cdf = pigt_impl(tt, (B_[r] + 0.5 * A_[r]) * inv_s,
                                       v_[r] * inv_s, 0.5 * A_[r] * inv_s);
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

inline void drdmgbm_raw(const double* rt, const double* const* cols, int n_rows,
                        const int* mask, const int* isok,
                        double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool global_kill = ctx ? ctx->has_global_kill() : false;
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const double* v_        = cols[emc2col::rdmgbm::v];
  const double* B_        = cols[emc2col::rdmgbm::B];
  const double* A_        = cols[emc2col::rdmgbm::A];
  const double* t0_       = cols[emc2col::rdmgbm::t0];
  const double* s_        = cols[emc2col::rdmgbm::s];
  const double* lambda_g_ = cols[emc2col::rdmgbm::mG];
  const double* lambda_k_ = cols[emc2col::rdmgbm::mK];
  const double* omega_    = (kill_shape == 3) ? cols[emc2col::rdmgbm::omega] : nullptr;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double lg = (!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lambda_g_[i], kill_shape);
    const double lk = (global_kill || !ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lambda_k_[i], kill_shape);
    const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(kill_shape);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (rt[i] <= 0.0)      { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    double log_pdf;
    const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
    if (dispatch.use_combo) {
      log_pdf = dgbm_local_combo(rt[i], v_[i], 1.0 + B_[i] + A_[i],
                                 A_[i], s_[i], t0_i, lg, lk, true, kill_shape, omega);
    } else {
      log_pdf = dgbm(rt[i],
                     v_[i],
                     1.0 + B_[i] + A_[i],
                     A_[i],
                     s_[i],
                     t0_i,
                     dispatch.lambda_g, dispatch.lambda_k,
                     true, kill_shape, dispatch.guess, omega);
    }
    out[i] = raw_log_value(log_pdf, min_ll, floor_raw);
  }
}

inline void prdmgbm_raw(const double* rt, const double* const* cols, int n_rows,
                        const int* mask, const int* isok,
                        double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool global_kill = ctx ? ctx->has_global_kill() : false;
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const double* v_        = cols[emc2col::rdmgbm::v];
  const double* B_        = cols[emc2col::rdmgbm::B];
  const double* A_        = cols[emc2col::rdmgbm::A];
  const double* t0_       = cols[emc2col::rdmgbm::t0];
  const double* s_        = cols[emc2col::rdmgbm::s];
  const double* lambda_g_ = cols[emc2col::rdmgbm::mG];
  const double* lambda_k_ = cols[emc2col::rdmgbm::mK];
  const double* omega_    = (kill_shape == 3) ? cols[emc2col::rdmgbm::omega] : nullptr;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double lg = (!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lambda_g_[i], kill_shape);
    const double lk = (global_kill || !ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lambda_k_[i], kill_shape);
    const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(kill_shape);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = 0.0; continue; }
    if (rt[i] <= 0.0)      { out[i] = 0.0; continue; }
    const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
    if (dispatch.use_combo) {
      const double log_cdf = pgbm_local_combo(rt[i], v_[i], 1.0 + B_[i] + A_[i],
                                              A_[i], s_[i], t0_i, lg, lk, true, kill_shape, omega);
      if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
      if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
      out[i] = log1m_exp(log_cdf);
      continue;
    }
    const double log_cdf = pgbm(rt[i],
                                v_[i],
                                1.0 + B_[i] + A_[i],
                                A_[i],
                                s_[i],
                                t0_i,
                                dispatch.lambda_g, dispatch.lambda_k,
                                true, kill_shape, dispatch.guess, omega);
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

// GBM: column layout v=0, B=1, A=2, t0=3, s=4, mG=5 (guess-clock mean), mK=6 (kill-clock mean)
inline void rdmgbm_logS_at_t(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool global_kill = ctx ? ctx->has_global_kill() : false;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* v_        = cols[emc2col::rdmgbm::v];
  const double* B_        = cols[emc2col::rdmgbm::B];
  const double* A_        = cols[emc2col::rdmgbm::A];
  const double* t0_       = cols[emc2col::rdmgbm::t0];
  const double* s_        = cols[emc2col::rdmgbm::s];
  const double* lambda_g_ = cols[emc2col::rdmgbm::mG];
  const double* lambda_k_ = cols[emc2col::rdmgbm::mK];
  const double* omega_    = (ks == 3) ? cols[emc2col::rdmgbm::omega] : nullptr;

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double t0_r = t0_[r];
      const double tt = t - t0_r;
      const double lg = (!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lambda_g_[r], ks);
      const double lk = (global_kill || !ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lambda_k_[r], ks);
      const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[r])) :
                           erlang_omega_for_shape(ks);
      const bool erl = (lg > 1e-12 || lk > 1e-12);
      const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
      if (tt <= 0.0) {
        if (!erl) continue;
        // EAM not started; erlang processes have been running since t=0
        if (dispatch.use_combo)
          logS += erlang_log_surv(t, lg, ks, omega) + erlang_log_surv(t, lk, ks, omega);
        else {
          const double lam = dispatch.guess ? dispatch.lambda_g : dispatch.lambda_k;
          if (lam > 1e-12) logS += erlang_log_surv(t, lam, ks, omega);
        }
        continue;
      }
      if (dispatch.use_combo) {
        const double log_cdf = pgbm_local_combo(t, v_[r], 1.0 + B_[r] + A_[r],
                                                A_[r], s_[r], t0_r, lg, lk, true, ks, omega);
        if (!R_FINITE(log_cdf) || log_cdf >= 0.0) { bad = true; break; }
        logS += log1m_exp(log_cdf);
        continue;
      }
      const double log_cdf = pgbm(t,
                                  v_[r],
                                  1.0 + B_[r] + A_[r],
                                  A_[r],
                                  s_[r],
                                  t0_r,
                                  dispatch.lambda_g, dispatch.lambda_k,
                                  true, ks, dispatch.guess, omega);
      if (!R_FINITE(log_cdf)) { bad = true; break; }
      if (log_cdf >= 0.0) { bad = true; break; }
      logS += log1m_exp(log_cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// Column order per src/col_registry.h.
inline void dlnr_raw(const double* rt, const double* const* cols, int n_rows,
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

inline void plnr_raw(const double* rt, const double* const* cols, int n_rows,
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

inline void lnr_logS_at_t(double t, const double* const* cols,
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

// Refactored PCOUNTER kernels.  These implement the closed forms:
// gamma-distributed trialwise input, pure-birth self-excitation,
// and a geometric threshold excess.  The raw `k` parameter is a non-negative
// offset; pcounter_k_int converts it to the canonical event threshold
// K = 2 + floor(k + 0.5) shared by every C++ entry point.
static constexpr double PC_EPS = 1e-8;
// Stirling numbers are generated one row at a time by pcounter_log_eval.
// There is deliberately no model-level upper bound on the finite threshold;
// pcounter_k_supported only rejects values that cannot be represented by the
// integer loop indices used by the exact finite-K formulas.

inline bool pcounter_k_supported(double k) {
  if (!R_FINITE(k) || k < 0.0) return false;
  const double z = std::floor(k + 0.5);
  // Tail formulas request rows through K+65; keep every int loop index
  // representable without imposing a fitted-model threshold cap.
  return z <= static_cast<double>(std::numeric_limits<int>::max() - 67);
}

inline bool pcounter_needs_stirling(double sv, double gamma) {
  return R_FINITE(sv) && R_FINITE(gamma) &&
         sv >= PC_EPS && gamma >= PC_EPS;
}

inline double pcounter_logsumexp(const std::vector<double>& x) {
  if (x.empty()) return R_NegInf;
  double m = R_NegInf;
  for (double z : x) if (z > m) m = z;
  if (!R_FINITE(m)) return R_NegInf;
  double s = 0.0;
  for (double z : x) if (R_FINITE(z)) s += std::exp(z - m);
  return m + std::log(s);
}

inline double pcounter_logdiffexp(double a, double b) {
  if (!R_FINITE(a)) return R_NegInf;
  if (!R_FINITE(b)) return a;
  if (b >= a) return R_NegInf;
  return a + std::log1p(-std::exp(b - a));
}

inline int pcounter_k_int(double k) {
  if (!pcounter_k_supported(k)) return 0;
  const double z = std::floor(k + 0.5);
  return static_cast<int>(z) + 2;
}

// A single rolling row of log Stirling numbers of the second kind.  The
// recurrence S(n,j) = S(n-1,j-1) + (n-1)S(n-1,j) is applied in descending j
// order, so the storage is O(n), rather than the former O(n^2) table.
class PcounterStirlingRows {
 public:
  PcounterStirlingRows() : n_(0), row_(1, 0.0) {}

  double at(int n, int j) {
    if (n < 0 || j < 0 || j > n) return R_NegInf;
    if (n < n_) {
      n_ = 0;
      row_.assign(1, 0.0);
    }
    while (n_ < n) {
      const int next = n_ + 1;
      row_.push_back(R_NegInf);
      for (int col = next; col >= 1; --col) {
        const double a = row_[static_cast<size_t>(col - 1)];
        const double log_b = (col < next && R_FINITE(row_[static_cast<size_t>(col)]))
          ? std::log(static_cast<double>(next - 1)) +
                row_[static_cast<size_t>(col)]
          : R_NegInf;
        row_[static_cast<size_t>(col)] =
          (a == R_NegInf) ? log_b :
          (log_b == R_NegInf) ? a :
          (a > log_b ? a + std::log1p(std::exp(log_b - a))
                      : log_b + std::log1p(std::exp(a - log_b)));
      }
      row_[0] = R_NegInf;
      n_ = next;
    }
    return row_[static_cast<size_t>(j)];
  }

 private:
  int n_;
  std::vector<double> row_;
};



inline double pcounter_log_rising(double a, int n) {
  double out = 0.0;
  for (int j = 0; j < n; ++j) out += std::log(a + static_cast<double>(j));
  return out;
}

inline void pcounter_log_lm(double a, double nu, double sv,
                            double& logL, double& logM) {
  if (sv <= PC_EPS) {
    logL = -nu * a;
    logM = std::log(nu) + logL;
    return;
  }
  const double shape = nu * nu / (sv * sv);
  const double rate = nu / (sv * sv);
  logL = -shape * std::log1p(a / rate);
  logM = std::log(shape) - std::log(rate + a) + logL;
}

inline double pcounter_log_h(int n, double t, double nu, double sv, double gamma,
                             PcounterStirlingRows& rows) {
  double logL, logM;
  pcounter_log_lm(t, nu, sv, logL, logM);
  if (sv <= PC_EPS)
    return logL + pcounter_log_rising(nu / gamma, n);
  const double shape = nu * nu / (sv * sv);
  const double rate = nu / (sv * sv);
  std::vector<double> terms;
  std::vector<double> rising(static_cast<size_t>(n + 1), 0.0);
  for (int j = 1; j <= n; ++j)
    rising[static_cast<size_t>(j)] =
      rising[static_cast<size_t>(j - 1)] +
      std::log(shape + static_cast<double>(j - 1));
  terms.reserve(static_cast<size_t>(n + 1));
  for (int j = 0; j <= n; ++j) {
    const double c = rows.at(n, j);
    if (!R_FINITE(c)) { terms.push_back(R_NegInf); continue; }
    terms.push_back(c + rising[static_cast<size_t>(j)] - j * std::log(gamma) -
                       j * std::log(rate + t));
  }
  return logL + pcounter_logsumexp(terms);
}

inline void pcounter_log_pi_phi(int n, double t, double nu, double sv, double gamma,
                                bool gamma_zero,
                                PcounterStirlingRows& rows,
                                double& logpi, double& logphi) {
  if (gamma_zero) {
    if (sv <= PC_EPS) {
      logpi = -nu * t;
      if (n > 0) logpi += n * std::log(nu * t) - R::lgammafn(n + 1.0);
      logphi = std::log(nu) + logpi;
      return;
    }
    const double shape = nu * nu / (sv * sv);
    const double rate = nu / (sv * sv);
    double logL, logM;
    pcounter_log_lm(t, nu, sv, logL, logM);
    logpi = logL + pcounter_log_rising(shape, n) - R::lgammafn(n + 1.0);
    if (n > 0) logpi += n * std::log(t);
    logpi -= n * std::log(rate + t);
    logphi = logpi + std::log(shape + n) - std::log(rate + t);
    return;
  }
  const double q = -std::expm1(-gamma * t);
  const double logq = q > 0.0 ? std::log(q) : R_NegInf;
  logpi = n * logq - R::lgammafn(n + 1.0) +
          pcounter_log_h(n, t, nu, sv, gamma, rows);
  logphi = std::log(gamma) + n * logq - R::lgammafn(n + 1.0) +
           pcounter_log_h(n + 1, t, nu, sv, gamma, rows);
}

inline double pcounter_log_tail(int start, double t, double nu, double sv, double gamma,
                                double logr, bool gamma_zero,
                                PcounterStirlingRows& rows, bool phi) {
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, rows, lp, lf);
    terms.push_back((phi ? lf : lp) + n * logr);
  }
  return pcounter_logsumexp(terms);
}

inline double pcounter_log_fixed_cdf(int start, double t, double nu, double sv,
                                     double gamma, bool gamma_zero,
                                     PcounterStirlingRows& rows) {
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, rows, lp, lf);
    terms.push_back(lp);
  }
  return pcounter_logsumexp(terms);
}

inline double pcounter_log_geom_cdf(int start, double t, double nu, double sv,
                                    double gamma, double omega, bool gamma_zero,
                                    PcounterStirlingRows& rows) {
  const double logr = std::log(omega) - std::log1p(omega);
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    if (n == start) { terms.push_back(R_NegInf); continue; }
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, rows, lp, lf);
    terms.push_back(lp + std::log(-std::expm1((n - start) * logr)));
  }
  return pcounter_logsumexp(terms);
}

inline void pcounter_log_eval(double t, double nu, double sv, double gamma,
                              double k, double omega,
                              double& logf, double& logS, double& logF) {
  logf = R_NegInf; logS = 0.0; logF = R_NegInf;
  if (ISNAN(t) || !R_FINITE(nu) || !R_FINITE(sv) || !R_FINITE(gamma) ||
      !pcounter_k_supported(k) || !R_FINITE(omega) || nu <= 0.0 ||
      sv < 0.0 || gamma < 0.0 || omega < 0.0) return;
  if (R_PosInf == t) { logS = R_NegInf; logF = 0.0; return; }
  if (!R_FINITE(t) || t <= 0.0) return;

  const int kk = pcounter_k_int(k);
  const bool gamma_zero = gamma < PC_EPS;
  const bool sv_zero = sv < PC_EPS;
  const bool omega_zero = omega < PC_EPS;
  PcounterStirlingRows rows;
  if (omega_zero) {
    std::vector<double> probs;
    probs.reserve(static_cast<size_t>(kk));
    for (int n = 0; n < kk; ++n) {
      double lp, lf;
      pcounter_log_pi_phi(n, t, nu, sv_zero ? 0.0 : sv,
                          gamma_zero ? 0.0 : gamma, gamma_zero, rows, lp, lf);
      probs.push_back(lp);
    }
    logS = std::min(0.0, pcounter_logsumexp(probs));
    double lp, lf;
    pcounter_log_pi_phi(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                        gamma_zero ? 0.0 : gamma, gamma_zero, rows, lp, lf);
    logf = lf;
    logF = logS > -1e-7 ? pcounter_log_fixed_cdf(kk, t, nu,
      sv_zero ? 0.0 : sv, gamma_zero ? 0.0 : gamma, gamma_zero, rows)
      : (logS < 0.0 ? std::log(-std::expm1(logS)) : R_NegInf);
    return;
  }

  const double logp = -std::log1p(omega);
  const double logr = std::log(omega) + logp;
  const double p = std::exp(logp);
  double ar, logd;
  if (gamma_zero) {
    ar = p * t; logd = 0.0;
  } else {
    const double q = -std::expm1(-gamma * t);
    logd = std::log1p(-std::exp(logr) * q);
    ar = t + logd / gamma;
  }
  double logLar, logMar;
  pcounter_log_lm(ar, nu, sv_zero ? 0.0 : sv, logLar, logMar);
  const int nlow = kk - 2;
  std::vector<double> low_pi, low_pi_r, low_phi_r;
  if (nlow >= 0) {
    low_pi.reserve(static_cast<size_t>(nlow + 1));
    low_pi_r.reserve(static_cast<size_t>(nlow + 1));
    low_phi_r.reserve(static_cast<size_t>(nlow + 1));
    for (int n = 0; n <= nlow; ++n) {
      double lp, lf;
      pcounter_log_pi_phi(n, t, nu, sv_zero ? 0.0 : sv,
                          gamma_zero ? 0.0 : gamma, gamma_zero, rows, lp, lf);
      low_pi.push_back(lp);
      low_pi_r.push_back(lp + n * logr);
      low_phi_r.push_back(lf + n * logr);
    }
  }
  const double loglow = pcounter_logsumexp(low_pi);
  const double loglowr = pcounter_logsumexp(low_pi_r);
  double logtail = pcounter_logdiffexp(logLar, loglowr);
  if (nlow >= 0 && (loglowr > logLar - 1e-7 || !R_FINITE(logtail)))
    logtail = pcounter_log_tail(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                gamma_zero ? 0.0 : gamma, logr, gamma_zero,
                                rows, false);
  logS = std::min(0.0, pcounter_logsumexp(
      std::vector<double>{loglow, (1.0 - kk) * logr + logtail}));
  const double logA = logMar - (gamma_zero ? 0.0 : logd);
  const double loglowf = pcounter_logsumexp(low_phi_r);
  double logtail_flux = pcounter_logdiffexp(logA, loglowf);
  logf = logp + (1.0 - kk) * logr + logtail_flux;
  if (nlow >= 0 && (loglowf > logA - 1e-7 || !R_FINITE(logf))) {
    logtail_flux = pcounter_log_tail(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                     gamma_zero ? 0.0 : gamma, logr, gamma_zero,
                                     rows, true);
    logf = logp + (1.0 - kk) * logr + logtail_flux;
  }
  if (logS > -1e-7)
    logF = pcounter_log_geom_cdf(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                 gamma_zero ? 0.0 : gamma, omega, gamma_zero,
                                 rows);
  else logF = logS < 0.0 ? std::log(-std::expm1(logS)) : R_NegInf;
}

inline double dpcounter_scalar(double t, const double* par, void* /*ctx_*/) {
  for (int j = 0; j < emc2col::pcounter::N_REQ; ++j)
    if (!R_FINITE(par[j])) return 0.0;
  if (par[emc2col::pcounter::nu] <= 0.0 ||
      par[emc2col::pcounter::sv] < 0.0 ||
      par[emc2col::pcounter::gamma] < 0.0 ||
      par[emc2col::pcounter::k] < 0.0 ||
      par[emc2col::pcounter::omega] < 0.0 ||
      !pcounter_k_supported(par[emc2col::pcounter::k])) return 0.0;
  double lf, ls, lF;
  pcounter_log_eval(t - par[emc2col::pcounter::t0], par[0], par[1], par[2],
                    par[3], par[4], lf, ls, lF);
  return R_FINITE(lf) ? std::exp(lf) : 0.0;
}

inline double ppcounter_scalar(double t, const double* par, void* /*ctx_*/) {
  for (int j = 0; j < emc2col::pcounter::N_REQ; ++j)
    if (!R_FINITE(par[j])) return 0.0;
  if (par[0] <= 0.0 || par[1] < 0.0 || par[2] < 0.0 || par[3] < 0.0 ||
      par[4] < 0.0 || !pcounter_k_supported(par[3])) return 0.0;
  double lf, ls, lF;
  pcounter_log_eval(t - par[5], par[0], par[1], par[2], par[3], par[4],
                    lf, ls, lF);
  return lF == 0.0 ? 1.0 : (R_FINITE(lF) ? std::exp(lF) : 0.0);
}

inline void dpcounter_raw(const double* rt, const double* const* cols, int n_rows,
                          const int* mask, const int* isok, double* out,
                          double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* nu = cols[emc2col::pcounter::nu];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* kk = cols[emc2col::pcounter::k];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || !R_FINITE(nu[i]) || !R_FINITE(sv[i]) ||
        !R_FINITE(ga[i]) || !R_FINITE(kk[i]) || !R_FINITE(om[i]) ||
        !R_FINITE(t0[i]) || nu[i] <= 0.0 || sv[i] < 0.0 ||
        ga[i] < 0.0 || om[i] < 0.0 || !pcounter_k_supported(kk[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    double lf, ls, lF;
    pcounter_log_eval(rt[i] - t0[i], nu[i], sv[i], ga[i], kk[i], om[i],
                      lf, ls, lF);
    out[i] = raw_log_value(lf, min_ll, floor_raw);
  }
}
inline void ppcounter_raw(const double* rt, const double* const* cols, int n_rows,
                          const int* mask, const int* isok, double* out,
                          double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* nu = cols[emc2col::pcounter::nu];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* kk = cols[emc2col::pcounter::k];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || !R_FINITE(nu[i]) || !R_FINITE(sv[i]) ||
        !R_FINITE(ga[i]) || !R_FINITE(kk[i]) || !R_FINITE(om[i]) ||
        !R_FINITE(t0[i]) || nu[i] <= 0.0 || sv[i] < 0.0 ||
        ga[i] < 0.0 || om[i] < 0.0 || !pcounter_k_supported(kk[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    double lf, ls, lF;
    pcounter_log_eval(rt[i] - t0[i], nu[i], sv[i], ga[i], kk[i], om[i],
                      lf, ls, lF);
    // Log-survivors are NOT floored at min_ll: a loser accumulator legitimately
    // carries large negative log-survival, and clamping it flattens the
    // likelihood surface out in the tails (and would disagree with
    // pcounter_logS_at_t, which sums the unfloored value for truncation).
    // Only an exactly-zero survivor falls back to raw_log_zero, as in
    // prdm_raw / plnr_raw / prexg_raw.
    out[i] = R_FINITE(ls) ? std::fmin(ls, 0.0)
                          : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pcounter_logS_at_t(double t, const double* const* cols,
                               int /*n_rows_total*/, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  (void)ctx_;
  const double* kk = cols[emc2col::pcounter::k];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  const double* nu = cols[emc2col::pcounter::nu];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double sum = 0.0;
    for (int k = 0; k < n_lR; ++k) {
      const int r = j * n_lR + k;
      if (!isok_all[r] || !R_FINITE(nu[r]) || !R_FINITE(sv[r]) ||
          !R_FINITE(ga[r]) || !R_FINITE(kk[r]) || !R_FINITE(om[r]) ||
          !R_FINITE(t0[r]) || nu[r] <= 0.0 || sv[r] < 0.0 ||
          ga[r] < 0.0 || om[r] < 0.0 || !pcounter_k_supported(kk[r])) {
        sum = R_NegInf; break;
      }
      double lf, ls, lF;
      pcounter_log_eval(t - t0[r], nu[r], sv[r], ga[r], kk[r], om[r],
                        lf, ls, lF);
      sum += ls;
    }
    logS_out[j] = sum;
  }
}

inline void drexg_raw(const double* rt, const double* const* cols, int n_rows,
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

inline void prexg_raw(const double* rt, const double* const* cols, int n_rows,
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

inline void rexg_logS_at_t(double t, const double* const* cols,
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

// ============================================================
// BAwL (Ballistic Accumulator with Leak + killing/guessing) adapters
// Column layout: v=0, sv=1, B=2, A=3, t0=4, k=5, mG=6 (guess-clock mean), mK=7 (kill-clock mean)
// ============================================================
inline double dbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[4];
  const double tt = t - t0_val;
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    if (tt <= 0.0 || t <= 0.0) return 0.0;
    // Scalar consumers (truncation normalisers, GSL integrands) clamp to
    // [0, 1], so use the natural-scale evaluator that admits saturation.
    return bawl_pdf_scalar_natural(
      tt, par[3], par[2] + par[3], par[0], par[1], 0.0,
      ctx->use_posdrift, LBA_DENOM_FLOOR
    );
  }
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const int ks = ctx ? ctx->kill_shape : 1;
  const double k_val = (ctx && ctx->bawl_k_fixed_zero) ? 0.0 : par[5];
  const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[6], ks) : 0.0);
  const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? erlang_lambda_from_mean(par[7], ks) : 0.0);
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  // Pass raw t and t0_val; core function splits EAM (t - t0) from erlang (t).
  return dkilledleakyba_norm(
    t, par[0], par[2] + par[3], par[3], par[1], t0_val, k_val, lg, lk,
    ctx->use_posdrift, false, ks, local_guess, omega, ctx->bawl_launch
  );
}

inline double pbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[4];
  const double tt = t - t0_val;
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    if (tt <= 0.0 || t <= 0.0) return 0.0;
    if (tt == R_PosInf)
      return ctx->use_posdrift ? 1.0 : pnorm_std(par[0] / par[1], true, false);
    return bawl_cdf_scalar_natural(
      tt, par[3], par[2] + par[3], par[0], par[1], 0.0,
      ctx->use_posdrift, LBA_DENOM_FLOOR
    );
  }
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const int ks = ctx ? ctx->kill_shape : 1;
  const double k_val = (ctx && ctx->bawl_k_fixed_zero) ? 0.0 : par[5];
  const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[6], ks) : 0.0);
  const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? erlang_lambda_from_mean(par[7], ks) : 0.0);
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  return pkilledleakyba_norm(
    t, par[0], par[2] + par[3], par[3], par[1], t0_val, k_val, lg, lk,
    ctx->use_posdrift, false, ks, local_guess, omega, ctx->bawl_launch
  );
}

inline void dbawl_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
  // Lognormal launch strengths occupy the v/sv slots as (mu, sigma); posdrift
  // and the normalizer floor are then inert inside the kernels.
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const double* v_  = cols[emc2col::bawl::v];
  const double* sv_  = cols[emc2col::bawl::sv];
  const double* B_  = cols[emc2col::bawl::B];
  const double* A_  = cols[emc2col::bawl::A];
  const double* t0_ = cols[emc2col::bawl::t0];
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    // Exact LBA member: retain the shared numerical kernel, but skip the
    // killed-clock wrapper and all optional-column bookkeeping per row.
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) {
        out[i] = raw_log_zero(min_ll, floor_raw);
        continue;
      }
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0 || rt[i] <= 0.0) {
        out[i] = raw_log_zero(min_ll, floor_raw);
        continue;
      }
      const double b_i = B_[i] + A_[i];
      double pdf = 0.0;
      if (ba_natural_pdf(tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd,
                         LBA_DENOM_FLOOR, BA_ACCEPT_RAW, pdf)) {
        out[i] = (pdf > 0.0)
          ? raw_log_value(std::log(pdf), min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      } else {
        const double log_pdf = log_ba_pdf(tt, A_[i], b_i, v_[i], sv_[i],
                                          0.0, pd, LBA_DENOM_FLOOR);
        out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
          ? raw_log_value(log_pdf, min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      }
    }
    return;
  }
  const bool local_guess = ctx->is_local_guess || ctx->is_local_kill_guess;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* k_  = (ctx && ctx->bawl_k_fixed_zero) ? nullptr : cols[emc2col::bawl::k];
  const double* lg_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[emc2col::bawl::mG];
  const double* lk_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[emc2col::bawl::mK];
  const double* omega_ = (ks == 3) ? cols[emc2col::bawl::omega] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double kval = (ctx && ctx->bawl_k_fixed_zero) ? 0.0 : k_[i];
    const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[i], ks));
    const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[i], ks));
    const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(ks);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (rt[i] <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (!erl) {
      // The no-clock BAwL member is evaluated directly in natural space in
      // the scalar adapter.  Keep that fast path in the raw callback too;
      // routing every conditional GH node through dkilledleakyba_norm(...,
      // log_out=true) needlessly enters the strict log wrapper and is costly
      // in the correlated model.  The log evaluator remains the tail fallback.
      double pdf = 0.0;
      if (ba_natural_pdf_launch(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval,
                                pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, pdf, lau)) {
        out[i] = (pdf > 0.0)
          ? raw_log_value(std::log(pdf), min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      } else {
        const double log_pdf = log_ba_pdf_launch(tt, A_[i], B_[i] + A_[i], v_[i],
                                                 sv_[i], kval, pd,
                                                 BAWL_DENOM_FLOOR, lau);
        out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
          ? raw_log_value(log_pdf, min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      }
      continue;
    }
    // Pass raw rt and t0; core function uses t0 to split EAM vs erlang time.
    const double log_pdf = dkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, kval, lg, lk,
      pd, true, ks, local_guess, omega, lau
    );
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbawl_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
  // Lognormal launch strengths occupy the v/sv slots as (mu, sigma); posdrift
  // and the normalizer floor are then inert inside the kernels.
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const double* v_  = cols[emc2col::bawl::v];
  const double* sv_  = cols[emc2col::bawl::sv];
  const double* B_  = cols[emc2col::bawl::B];
  const double* A_  = cols[emc2col::bawl::A];
  const double* t0_ = cols[emc2col::bawl::t0];
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    // Exact LBA member; pfun writes the log-survivor expected by the raw
    // likelihood path, while log_ba_cdf supplies a stable log CDF.
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) {
        out[i] = 0.0;
        continue;
      }
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0 || rt[i] <= 0.0) {
        out[i] = 0.0;
        continue;
      }
      if (tt == R_PosInf) {
        const double ls = log_ba_surv_launch(tt, A_[i], B_[i] + A_[i],
                                             v_[i], sv_[i], 0.0, pd,
                                             LBA_DENOM_FLOOR, lau);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
        continue;
      }
      const double b_i = B_[i] + A_[i];
      double cdf = 0.0;
      if (ba_natural_cdf(tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd,
                         LBA_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
        out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
      } else {
        const double ls = log_ba_surv_launch(tt, A_[i], b_i, v_[i], sv_[i],
                                             0.0, pd, LBA_DENOM_FLOOR, lau);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
      }
    }
    return;
  }
  const bool local_guess = ctx->is_local_guess || ctx->is_local_kill_guess;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* k_  = (ctx && ctx->bawl_k_fixed_zero) ? nullptr : cols[emc2col::bawl::k];
  const double* lg_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[emc2col::bawl::mG];
  const double* lk_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[emc2col::bawl::mK];
  const double* omega_ = (ks == 3) ? cols[emc2col::bawl::omega] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double kval = (ctx && ctx->bawl_k_fixed_zero) ? 0.0 : k_[i];
    const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[i], ks));
    const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[i], ks));
    const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(ks);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = 0.0; continue; }
    if (rt[i] <= 0.0) { out[i] = 0.0; continue; }
    if (!erl) {
      double cdf = 0.0;
      if (ba_natural_cdf_launch(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval,
                                pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, lau)) {
        out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
      } else {
        const double ls = log_ba_surv_launch(tt, A_[i], B_[i] + A_[i], v_[i],
                                             sv_[i], kval, pd,
                                             BAWL_DENOM_FLOOR, lau);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
      }
      continue;
    }
    const double log_cdf = pkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, kval, lg, lk,
      pd, true, ks, local_guess, omega, lau
    );
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

// BAwL: column layout v=0, sv=1, B=2, A=3, t0=4, k=5, mG=6 (guess-clock mean), mK=7 (kill-clock mean)
inline void bawl_logS_at_t(double t, const double* const* cols,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd          = ctx->use_posdrift;
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* v_  = cols[emc2col::bawl::v];
  const double* sv_ = cols[emc2col::bawl::sv];
  const double* B_  = cols[emc2col::bawl::B];
  const double* A_  = cols[emc2col::bawl::A];
  const double* t0_ = cols[emc2col::bawl::t0];
  const double* k_  = (ctx && ctx->bawl_k_fixed_zero) ? nullptr : cols[emc2col::bawl::k];
  const double* lg_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[emc2col::bawl::mG];
  const double* lk_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[emc2col::bawl::mK];
  const double* omega_ = (ks == 3) ? cols[emc2col::bawl::omega] : nullptr;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double t0_r = t0_[r];
      const double tt = t - t0_r;
      if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
        if (tt <= 0.0) continue;
        if (tt == R_PosInf) {
          const double ls = log_ba_surv_launch(tt, A_[r], B_[r] + A_[r],
                                               v_[r], sv_[r], 0.0, pd,
                                               LBA_DENOM_FLOOR, lau);
          if (!(ls > R_NegInf)) { bad = true; break; }
          logS += ls;
          continue;
        }
        double cdf = 0.0;
        if (ba_natural_cdf(tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], 0.0, pd,
                           LBA_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
          if (cdf > 0.0) logS += std::log1p(-cdf);
        } else {
          const double ls = log_ba_surv_launch(tt, A_[r], B_[r] + A_[r], v_[r],
                                               sv_[r], 0.0, pd, LBA_DENOM_FLOOR, lau);
          if (!(ls > R_NegInf)) { bad = true; break; }
          logS += ls;
        }
        continue;
      }
      const double kval = (ctx && ctx->bawl_k_fixed_zero) ? 0.0 : k_[r];
      const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
        ((!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[r], ks));
      const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
        ((!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[r], ks));
      const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[r])) :
                           erlang_omega_for_shape(ks);
      const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
      const bool erl = (lg > 1e-12 || lk > 1e-12);
      if (tt <= 0.0) {
        if (!erl) continue;  // EAM not started, no erlang → log-survivor += 0
        // Before EAM onset, pure kill produces no observed response mass, while
        // guess paths can.  Use the same CDF logic as the scalar path so kill
        // before t0 is not treated as an observed hit removed by lower truncation.
        const double log_cdf = pkilledleakyba_norm(
          t, v_[r], B_[r] + A_[r], A_[r], sv_[r], t0_r, kval, lg, lk,
          pd, true, ks, local_guess, omega, lau
        );
        if (!R_FINITE(log_cdf)) continue;
        if (log_cdf >= 0.0) { bad = true; break; }
        logS += log1m_exp(log_cdf);
        continue;
      }
      if (!erl) {
        double cdf = 0.0;
        if (ba_natural_cdf_launch(tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], kval,
                                  pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, lau)) {
          if (cdf > 0.0) logS += std::log1p(-cdf);
        } else {
          const double ls = log_ba_surv_launch(tt, A_[r], B_[r] + A_[r], v_[r],
                                               sv_[r], kval, pd,
                                               BAWL_DENOM_FLOOR, lau);
          if (!(ls > R_NegInf) || ISNAN(ls)) { bad = true; break; }
          logS += ls;
        }
        continue;
      }
      // Both EAM and erlang contribute; pass raw t and t0_r.
      const double log_cdf = pkilledleakyba_norm(
        t, v_[r], B_[r] + A_[r], A_[r], sv_[r], t0_r, kval, lg, lk,
        pd, true, ks, local_guess, omega, lau
      );
      if (log_cdf >= 0.0) { bad = true; break; }
      if (R_FINITE(log_cdf)) logS += log1m_exp(log_cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// ============================================================
// BAwD (Ballistic Accumulator with drive Decay) adapters
// Column layout: p1=0 (v | mu), p2=1 (sv | sigma), B=2, A=3, t0=4, k=5, ell=6.
// There are no clocks and no optional columns, so every row reads the same
// seven values; ctx->bawd_launch selects the launch distribution.
//
// The upper tail is ALWAYS defective (never-finish mass exists whenever
// ell > 0 or the drift can fall short), and the CDF is exactly flat past
// t0 + T_max.  pfun_raw returns the log-SURVIVOR, per the race kernel
// contract, so the flat tail arrives as a constant log(1 - F_max) rather than
// as a clamped zero.
// ============================================================

inline int bawd_launch_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_launch : BAWD_LAUNCH_LOGNORMAL;
}
inline double bawd_gamma_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_gamma : 0.0;
}
inline double bawd_rho_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_rho : R_PosInf;
}

// Slot `clear` holds T_max under the endpoint chart and ell at gamma = 1; the
// kernels below all take ell, so every column read goes through this.
inline double bawd_clear_to_ell(const ContextForRaceModels* ctx, double clear,
                                double b, double k) {
  const double gamma = bawd_gamma_of(ctx);
  if (!bawd_uses_tmax(gamma)) return clear;
  return bawd_ell_from_tmax(b, k, clear, gamma, bawd_rho_of(ctx));
}

inline double dbawd_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawd::v])) return 0.0;
  const double tt = t - par[emc2col::bawd::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  const double b = par[emc2col::bawd::B] + par[emc2col::bawd::A];
  return bawd_pdf_scalar_natural(
    tt, par[emc2col::bawd::A], b,
    par[emc2col::bawd::v], par[emc2col::bawd::sv],
    par[emc2col::bawd::k],
    bawd_clear_to_ell(ctx, par[emc2col::bawd::clear], b,
                      par[emc2col::bawd::k]),
    bawd_launch_of(ctx), ctx ? ctx->use_posdrift : true,
    bawd_gamma_of(ctx), bawd_rho_of(ctx));
}

inline double pbawd_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawd::v])) return 0.0;
  const double tt = t - par[emc2col::bawd::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf is handled inside the kernel and returns F_max, not 1.
  const double b = par[emc2col::bawd::B] + par[emc2col::bawd::A];
  return bawd_cdf_scalar_natural(
    tt, par[emc2col::bawd::A], b,
    par[emc2col::bawd::v], par[emc2col::bawd::sv],
    par[emc2col::bawd::k],
    bawd_clear_to_ell(ctx, par[emc2col::bawd::clear], b,
                      par[emc2col::bawd::k]),
    bawd_launch_of(ctx), ctx ? ctx->use_posdrift : true,
    bawd_gamma_of(ctx), bawd_rho_of(ctx));
}

inline void dbawd_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawd_launch_of(ctx);
  const double gamma = bawd_gamma_of(ctx);
  const double rho = bawd_rho_of(ctx);
  const double* p1_ = cols[emc2col::bawd::v];
  const double* p2_ = cols[emc2col::bawd::sv];
  const double* B_  = cols[emc2col::bawd::B];
  const double* A_  = cols[emc2col::bawd::A];
  const double* t0_ = cols[emc2col::bawd::t0];
  const double* k_  = cols[emc2col::bawd::k];
  const double* clear_ = cols[emc2col::bawd::clear];
  const bool tmax_chart = bawd_uses_tmax(gamma);
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
    const double b_i = B_[i] + A_[i];
    const double ell_i = tmax_chart
      ? bawd_ell_from_tmax(b_i, k_[i], clear_[i], gamma, rho) : clear_[i];
    const double log_pdf = bawd_log_pdf(tt, A_[i], b_i, p1_[i],
                                        p2_[i], k_[i], ell_i, launch, pd,
                                        gamma, rho);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbawd_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawd_launch_of(ctx);
  const double gamma = bawd_gamma_of(ctx);
  const double rho = bawd_rho_of(ctx);
  const double* p1_ = cols[emc2col::bawd::v];
  const double* p2_ = cols[emc2col::bawd::sv];
  const double* B_  = cols[emc2col::bawd::B];
  const double* A_  = cols[emc2col::bawd::A];
  const double* t0_ = cols[emc2col::bawd::t0];
  const double* k_  = cols[emc2col::bawd::k];
  const double* clear_ = cols[emc2col::bawd::clear];
  const bool tmax_chart = bawd_uses_tmax(gamma);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const double b_i = B_[i] + A_[i];
    const double ell_i = tmax_chart
      ? bawd_ell_from_tmax(b_i, k_[i], clear_[i], gamma, rho) : clear_[i];
    double cdf = 0.0;
    if (ba_natural_cdf_bawd(tt, A_[i], b_i, p1_[i], p2_[i], k_[i], ell_i,
                            launch, pd, gamma, rho, BAWD_DENOM_FLOOR,
                            BA_ACCEPT_RAW, cdf)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double log_s = bawd_log_surv(tt, A_[i], b_i, p1_[i], p2_[i],
                                         k_[i], ell_i, launch, pd, gamma, rho,
                                         BAWD_DENOM_FLOOR);
      out[i] = (log_s > R_NegInf && emc2_isfinite(log_s))
        ? log_s : raw_log_zero(min_ll, floor_raw);
    }
  }
}

// ============================================================
// BTAwL transient-member adapters. The full local race has a separate
// nine-column adapter below; both use p1=0 (v | mu), p2=1 (sv | sigma).
// ============================================================

inline int btawl_launch_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->btawl_launch : BTAWL_LAUNCH_NORMAL;
}

inline bool btawl_uses_ttrans(const ContextForRaceModels* ctx) {
  return ctx && ctx->btawl_ttrans_chart;
}

inline double btawl_tau_of(const ContextForRaceModels* ctx, double clear,
                           double k) {
  if (!btawl_uses_ttrans(ctx)) return clear;
  if (ctx != nullptr) {
    for (int j = 0; j < ContextForRaceModels::btawl_tau_cache_size; ++j) {
      if (ctx->btawl_tau_cache_k[j] == k &&
          ctx->btawl_tau_cache_clear[j] == clear)
        return ctx->btawl_tau_cache_value[j];
    }
  }
  const double tau = btawl_tau_from_ttrans(k, clear);
  if (ctx != nullptr) {
    const int j = ctx->btawl_tau_cache_next;
    ctx->btawl_tau_cache_k[j] = k;
    ctx->btawl_tau_cache_clear[j] = clear;
    ctx->btawl_tau_cache_value[j] = tau;
    ctx->btawl_tau_cache_next =
      (j + 1) % ContextForRaceModels::btawl_tau_cache_size;
  }
  return tau;
}

inline double dbtawl_transient_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::btawl_transient::v])) return 0.0;
  const double tt = t - par[emc2col::btawl_transient::t0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[emc2col::btawl_transient::B] + par[emc2col::btawl_transient::A];
  if (btawl_uses_ttrans(ctx))
    return btawl_pdf_chart(tt, par[emc2col::btawl_transient::A], b,
                           par[emc2col::btawl_transient::v], par[emc2col::btawl_transient::sv],
                           par[emc2col::btawl_transient::k], par[emc2col::btawl_transient::clear],
                           btawl_launch_of(ctx), ctx ? ctx->use_posdrift : true, true,
                           btawl_tau_of(ctx, par[emc2col::btawl_transient::clear],
                                        par[emc2col::btawl_transient::k]));
  return btawl_pdf(tt, par[emc2col::btawl_transient::A], b,
                   par[emc2col::btawl_transient::v], par[emc2col::btawl_transient::sv],
                   par[emc2col::btawl_transient::k],
                   btawl_tau_of(ctx, par[emc2col::btawl_transient::clear], par[emc2col::btawl_transient::k]),
                   btawl_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline double pbtawl_transient_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::btawl_transient::v])) return 0.0;
  const double tt = t - par[emc2col::btawl_transient::t0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[emc2col::btawl_transient::B] + par[emc2col::btawl_transient::A];
  if (btawl_uses_ttrans(ctx))
    return btawl_cdf_chart(tt, par[emc2col::btawl_transient::A], b,
                           par[emc2col::btawl_transient::v], par[emc2col::btawl_transient::sv],
                           par[emc2col::btawl_transient::k], par[emc2col::btawl_transient::clear],
                           btawl_launch_of(ctx), ctx ? ctx->use_posdrift : true, true,
                           btawl_tau_of(ctx, par[emc2col::btawl_transient::clear],
                                        par[emc2col::btawl_transient::k]));
  return btawl_cdf(tt, par[emc2col::btawl_transient::A], b,
                   par[emc2col::btawl_transient::v], par[emc2col::btawl_transient::sv],
                   par[emc2col::btawl_transient::k],
                   btawl_tau_of(ctx, par[emc2col::btawl_transient::clear], par[emc2col::btawl_transient::k]),
                   btawl_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline void dbtawl_transient_raw(const double* rt, const double* const* cols, int n_rows,
                       const int* mask, const int* isok, double* out,
                       double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_transient::v];
  const double* p2 = cols[emc2col::btawl_transient::sv];
  const double* B = cols[emc2col::btawl_transient::B];
  const double* A = cols[emc2col::btawl_transient::A];
  const double* t0 = cols[emc2col::btawl_transient::t0];
  const double* k = cols[emc2col::btawl_transient::k];
  const double* clear = cols[emc2col::btawl_transient::clear];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double tt = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(tt > 0.0)) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double lp = btawl_uses_ttrans(ctx)
      ? btawl_log_pdf_chart(tt, A[i], B[i] + A[i], p1[i], p2[i], k[i], clear[i],
                            launch, pd, true, btawl_tau_of(ctx, clear[i], k[i]))
      : btawl_pdf_log(tt, A[i], B[i] + A[i], p1[i], p2[i], k[i], clear[i], launch, pd);
    out[i] = (lp > R_NegInf && emc2_isfinite(lp))
      ? raw_log_value(lp, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbtawl_transient_raw(const double* rt, const double* const* cols, int n_rows,
                       const int* mask, const int* isok, double* out,
                       double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_transient::v];
  const double* p2 = cols[emc2col::btawl_transient::sv];
  const double* B = cols[emc2col::btawl_transient::B];
  const double* A = cols[emc2col::btawl_transient::A];
  const double* t0 = cols[emc2col::btawl_transient::t0];
  const double* k = cols[emc2col::btawl_transient::k];
  const double* clear = cols[emc2col::btawl_transient::clear];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(tt > 0.0)) { out[i] = 0.0; continue; }
    const double tau = btawl_tau_of(ctx, clear[i], k[i]);
    const BtawlGeom g = btawl_uses_ttrans(ctx)
      ? btawl_geometry_from_ttrans(A[i], B[i] + A[i], k[i], clear[i], tau)
      : btawl_geometry(A[i], B[i] + A[i], k[i], tau);
    double cdf = 0.0;
    if (btawl_natural_cdf_from_geom(tt, g, p1[i], p2[i], launch, pd, cdf)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = launch == BTAWL_LAUNCH_LOGNORMAL
        ? log_btawl_surv_logn(tt, g, p1[i], p2[i])
        : log_btawl_surv_normal(tt, g, p1[i], p2[i], pd);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

inline void btawl_transient_logS_at_t(double t, const double* const* cols,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_transient::v];
  const double* p2 = cols[emc2col::btawl_transient::sv];
  const double* B = cols[emc2col::btawl_transient::B];
  const double* A = cols[emc2col::btawl_transient::A];
  const double* t0 = cols[emc2col::btawl_transient::t0];
  const double* k = cols[emc2col::btawl_transient::k];
  const double* clear = cols[emc2col::btawl_transient::clear];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double ls = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = j * n_lR + kk;
      if (r >= n_rows_total || !isok_all[r] || R_IsNA(p1[r])) { bad = true; break; }
      const double tt = t - t0[r];
      if (!(tt > 0.0)) continue;
      const double lsr = btawl_uses_ttrans(ctx)
        ? btawl_log_surv_chart(tt, A[r], B[r] + A[r], p1[r], p2[r], k[r], clear[r],
                              launch, pd, true, btawl_tau_of(ctx, clear[r], k[r]))
        : btawl_log_surv(tt, A[r], B[r] + A[r], p1[r], p2[r], k[r], clear[r], launch, pd);
      if (!(lsr > R_NegInf) || ISNAN(lsr)) { bad = true; break; }
      ls += lsr;
    }
    logS_out[j] = bad ? R_NegInf : ls;
  }
}

// BTAwL sustained/transient local race.
inline double dbtawl_local_race_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::btawl_local_race::v])) return 0.0;
  const double tt = t - par[emc2col::btawl_local_race::t0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[emc2col::btawl_local_race::B] + par[emc2col::btawl_local_race::A];
  return btawl_local_race_pdf(tt, par[emc2col::btawl_local_race::A], b,
                       par[emc2col::btawl_local_race::v], par[emc2col::btawl_local_race::sv],
                       par[emc2col::btawl_local_race::k], par[emc2col::btawl_local_race::tau_s],
                       btawl_tau_of(ctx, par[emc2col::btawl_local_race::clear], par[emc2col::btawl_local_race::k]), par[emc2col::btawl_local_race::pi],
                       btawl_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline double pbtawl_local_race_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::btawl_local_race::v])) return 0.0;
  const double tt = t - par[emc2col::btawl_local_race::t0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[emc2col::btawl_local_race::B] + par[emc2col::btawl_local_race::A];
  return btawl_local_race_cdf(tt, par[emc2col::btawl_local_race::A], b,
                       par[emc2col::btawl_local_race::v], par[emc2col::btawl_local_race::sv],
                       par[emc2col::btawl_local_race::k], par[emc2col::btawl_local_race::tau_s],
                       btawl_tau_of(ctx, par[emc2col::btawl_local_race::clear], par[emc2col::btawl_local_race::k]), par[emc2col::btawl_local_race::pi],
                       btawl_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline void dbtawl_local_race_raw(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok, double* out,
                           double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx); const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_local_race::v];
  const double* p2 = cols[emc2col::btawl_local_race::sv];
  const double* B = cols[emc2col::btawl_local_race::B];
  const double* A = cols[emc2col::btawl_local_race::A];
  const double* t0 = cols[emc2col::btawl_local_race::t0];
  const double* k = cols[emc2col::btawl_local_race::k];
  const double* ts = cols[emc2col::btawl_local_race::tau_s];
  const double* tt_clear = cols[emc2col::btawl_local_race::clear];
  const double* pi = cols[emc2col::btawl_local_race::pi];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double lp = btawl_local_race_pdf_log(u, A[i], B[i] + A[i], p1[i], p2[i], k[i],
                                        ts[i], btawl_tau_of(ctx, tt_clear[i], k[i]), pi[i], launch, pd);
    out[i] = (lp > R_NegInf && emc2_isfinite(lp))
      ? raw_log_value(lp, min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbtawl_local_race_raw(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok, double* out,
                           double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx); const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_local_race::v];
  const double* p2 = cols[emc2col::btawl_local_race::sv];
  const double* B = cols[emc2col::btawl_local_race::B];
  const double* A = cols[emc2col::btawl_local_race::A];
  const double* t0 = cols[emc2col::btawl_local_race::t0];
  const double* k = cols[emc2col::btawl_local_race::k];
  const double* ts = cols[emc2col::btawl_local_race::tau_s];
  const double* tt_clear = cols[emc2col::btawl_local_race::clear];
  const double* pi = cols[emc2col::btawl_local_race::pi];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = 0.0; continue; }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) { out[i] = 0.0; continue; }
    const double tau_t = btawl_tau_of(ctx, tt_clear[i], k[i]);
    const double cdf = btawl_local_race_cdf(u, A[i], B[i] + A[i], p1[i], p2[i], k[i],
                                     ts[i], tau_t, pi[i], launch, pd);
    if (emc2_isfinite(cdf) && cdf >= 0.0 && cdf < 1.0 - 1e-8) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = btawl_local_race_log_surv(u, A[i], B[i] + A[i], p1[i], p2[i], k[i],
                                           ts[i], tau_t, pi[i], launch, pd);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

inline void btawl_local_race_logS_at_t(double t, const double* const* cols,
                                int n_rows_total, int n_lR, int /*n_par*/,
                                const int* trunc_mask, int n_unique_trials,
                                const int* isok_all, void* ctx_, double* out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx); const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_local_race::v];
  const double* p2 = cols[emc2col::btawl_local_race::sv];
  const double* B = cols[emc2col::btawl_local_race::B];
  const double* A = cols[emc2col::btawl_local_race::A];
  const double* t0 = cols[emc2col::btawl_local_race::t0];
  const double* k = cols[emc2col::btawl_local_race::k];
  const double* ts = cols[emc2col::btawl_local_race::tau_s];
  const double* tt_clear = cols[emc2col::btawl_local_race::clear];
  const double* pi = cols[emc2col::btawl_local_race::pi];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double ls = 0.0; bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = j * n_lR + kk;
      if (r >= n_rows_total || !isok_all[r] || R_IsNA(p1[r])) { bad = true; break; }
      const double u = t - t0[r]; if (!(u > 0.0)) continue;
      const double lsr = btawl_local_race_log_surv(u, A[r], B[r] + A[r], p1[r], p2[r], k[r],
                                            ts[r], btawl_tau_of(ctx, tt_clear[r], k[r]),
                                            pi[r], launch, pd);
      if (!(lsr > R_NegInf) || ISNAN(lsr)) { bad = true; break; }
      ls += lsr;
    }
    out[j] = bad ? R_NegInf : ls;
  }
}

inline double dbtawl_sustained_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const double tt = t - par[emc2col::btawl_sustained::t0];
  if (R_IsNA(par[emc2col::btawl_sustained::v]) || !(t > 0.0) || !(tt > 0.0))
    return 0.0;
  return btawl_sustained_pdf(tt, par[emc2col::btawl_sustained::A],
                             par[emc2col::btawl_sustained::B] +
                               par[emc2col::btawl_sustained::A],
                             par[emc2col::btawl_sustained::v],
                             par[emc2col::btawl_sustained::sv],
                             par[emc2col::btawl_sustained::k],
                             par[emc2col::btawl_sustained::tau_s], launch,
                             ctx ? ctx->use_posdrift : true);
}

inline double pbtawl_sustained_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const double tt = t - par[emc2col::btawl_sustained::t0];
  if (R_IsNA(par[emc2col::btawl_sustained::v]) || !(t > 0.0) || !(tt > 0.0))
    return 0.0;
  return btawl_sustained_cdf(tt, par[emc2col::btawl_sustained::A],
                             par[emc2col::btawl_sustained::B] +
                               par[emc2col::btawl_sustained::A],
                             par[emc2col::btawl_sustained::v],
                             par[emc2col::btawl_sustained::sv],
                             par[emc2col::btawl_sustained::k],
                             par[emc2col::btawl_sustained::tau_s], launch,
                             ctx ? ctx->use_posdrift : true);
}

inline void dbtawl_sustained_raw(const double* rt, const double* const* cols,
                                 int n_rows, const int* mask, const int* isok,
                                 double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* p1 = cols[emc2col::btawl_sustained::v];
  const double* p2 = cols[emc2col::btawl_sustained::sv];
  const double* B = cols[emc2col::btawl_sustained::B];
  const double* A = cols[emc2col::btawl_sustained::A];
  const double* t0 = cols[emc2col::btawl_sustained::t0];
  const double* k = cols[emc2col::btawl_sustained::k];
  const double* ts = cols[emc2col::btawl_sustained::tau_s];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double p = btawl_sustained_pdf(u, A[i], B[i] + A[i], p1[i], p2[i],
                                         k[i], ts[i], launch, pd);
    out[i] = (p > 0.0 && emc2_isfinite(p))
      ? raw_log_value(std::log(p), min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbtawl_sustained_raw(const double* rt, const double* const* cols,
                                 int n_rows, const int* mask, const int* isok,
                                 double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* p1 = cols[emc2col::btawl_sustained::v];
  const double* p2 = cols[emc2col::btawl_sustained::sv];
  const double* B = cols[emc2col::btawl_sustained::B];
  const double* A = cols[emc2col::btawl_sustained::A];
  const double* t0 = cols[emc2col::btawl_sustained::t0];
  const double* k = cols[emc2col::btawl_sustained::k];
  const double* ts = cols[emc2col::btawl_sustained::tau_s];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = 0.0; continue; }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) { out[i] = 0.0; continue; }
    const double cdf = btawl_sustained_cdf(u, A[i], B[i] + A[i], p1[i], p2[i],
                                           k[i], ts[i], launch, pd);
    if (emc2_isfinite(cdf) && cdf >= 0.0 && cdf < 1.0 - 1e-8) {
      out[i] = cdf > 0.0 ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = btawl_sustained_log_surv(u, A[i], B[i] + A[i], p1[i], p2[i],
                                                  k[i], ts[i], launch, pd);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

inline void btawl_sustained_logS_at_t(double t, const double* const* cols,
                                      int n_rows_total, int n_lR, int /*n_par*/,
                                      const int* trunc_mask, int n_unique_trials,
                                      const int* isok_all, void* ctx_, double* out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* p1 = cols[emc2col::btawl_sustained::v];
  const double* p2 = cols[emc2col::btawl_sustained::sv];
  const double* B = cols[emc2col::btawl_sustained::B];
  const double* A = cols[emc2col::btawl_sustained::A];
  const double* t0 = cols[emc2col::btawl_sustained::t0];
  const double* k = cols[emc2col::btawl_sustained::k];
  const double* ts = cols[emc2col::btawl_sustained::tau_s];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double ls = 0.0; bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = j * n_lR + kk;
      if (r >= n_rows_total || !isok_all[r] || R_IsNA(p1[r])) { bad = true; break; }
      const double u = t - t0[r]; if (!(u > 0.0)) continue;
      const double lsr = btawl_sustained_log_surv(u, A[r], B[r] + A[r], p1[r], p2[r],
                                                  k[r], ts[r], launch, pd);
      if (!(lsr > R_NegInf) || ISNAN(lsr)) { bad = true; break; }
      ls += lsr;
    }
    out[j] = bad ? R_NegInf : ls;
  }
}

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

inline void bawd_logS_at_t(double t, const double* const* cols,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawd_launch_of(ctx);
  const double gamma = bawd_gamma_of(ctx);
  const double rho = bawd_rho_of(ctx);
  const double* p1_ = cols[emc2col::bawd::v];
  const double* p2_ = cols[emc2col::bawd::sv];
  const double* B_  = cols[emc2col::bawd::B];
  const double* A_  = cols[emc2col::bawd::A];
  const double* t0_ = cols[emc2col::bawd::t0];
  const double* k_  = cols[emc2col::bawd::k];
  const double* clear_ = cols[emc2col::bawd::clear];
  const bool tmax_chart = bawd_uses_tmax(gamma);
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
      const double b_r = B_[r] + A_[r];
      const double ell_r = tmax_chart
        ? bawd_ell_from_tmax(b_r, k_[r], clear_[r], gamma, rho) : clear_[r];
      double cdf = 0.0;
      if (ba_natural_cdf_bawd(tt, A_[r], b_r, p1_[r], p2_[r], k_[r], ell_r,
                              launch, pd, gamma, rho, BAWD_DENOM_FLOOR,
                              BA_ACCEPT_RAW, cdf)) {
        if (cdf > 0.0) logS += std::log1p(-cdf);
      } else {
        const double log_s = bawd_log_surv(tt, A_[r], b_r, p1_[r], p2_[r],
                                           k_[r], ell_r, launch, pd, gamma, rho,
                                           BAWD_DENOM_FLOOR);
        if (!(log_s > R_NegInf) || ISNAN(log_s)) { bad = true; break; }
        logS += log_s;
      }
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
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

// ============================================================
// FRQ (Finite Reservoir Quorum) adapters
// Column layout: alpha=0, beta=1, h=2, tau=3, t0=4, delta=5.  No optional
// columns and no context flags -- the model has a single variant, and the
// conditional quantile defining tau is the compile-time constant FRQ_QUANTILE.
// delta is the threshold-variability half-width; it is zero by default, and
// frq_derive() then builds an inactive generator that costs nothing.
//
// The upper tail is ALWAYS defective: an accumulator terminates only with
// probability h = I_p(alpha, beta) < 1, and the leftover mass 1 - h sits at
// t = +Inf.  pfun_raw therefore returns the log-SURVIVOR (the race kernel
// contract), which saturates at log(1 - h) rather than falling to -Inf, and
// that is exactly what makes an omission trial score sum_j log(1 - h_j).
//
// Every entry point derives (p, lambda) through frq_derive(), memoised on the
// exact parameter bits because the two qbeta inversions dominate the cost and
// consecutive compressed rows usually repeat.
// ============================================================

inline double dfrq_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[emc2col::frq::alpha]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  const FrqPars s = frq_derive(par[emc2col::frq::alpha], par[emc2col::frq::beta],
                               par[emc2col::frq::h], par[emc2col::frq::tau],
                               par[emc2col::frq::delta]);
  return frq_pdf_natural_dt(tt, s);
}
inline double pfrq_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[emc2col::frq::alpha]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf deliberately reaches the kernel: the CDF there is h, not one.
  const FrqPars s = frq_derive(par[emc2col::frq::alpha], par[emc2col::frq::beta],
                               par[emc2col::frq::h], par[emc2col::frq::tau],
                               par[emc2col::frq::delta]);
  return frq_cdf_natural_dt(tt, s);
}

inline void dfrq_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* al_ = cols[emc2col::frq::alpha];
  const double* be_ = cols[emc2col::frq::beta];
  const double* h_  = cols[emc2col::frq::h];
  const double* ta_ = cols[emc2col::frq::tau];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  FrqMemo memo;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(al_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const FrqPars& s = memo.get(al_[i], be_[i], h_[i], ta_[i], de_[i]);
    const double log_pdf = frq_log_pdf_dt(tt, s);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pfrq_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* al_ = cols[emc2col::frq::alpha];
  const double* be_ = cols[emc2col::frq::beta];
  const double* h_  = cols[emc2col::frq::h];
  const double* ta_ = cols[emc2col::frq::tau];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  FrqMemo memo;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one, matching
    // every other race adapter: the trial is failed by its winner's density.
    if (R_IsNA(al_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = 0.0; continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const FrqPars& s = memo.get(al_[i], be_[i], h_[i], ta_[i], de_[i]);
    const double log_surv = frq_log_surv_dt(tt, s);
    if (ISNAN(log_surv)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (log_surv > R_NegInf) ? log_surv
                                   : raw_log_zero(min_ll, floor_raw);
  }
}

inline void frq_logS_at_t(double t, const double* const* cols,
                          int /*n_rows_total*/, int n_lR, int /*n_par*/,
                          const int* trunc_mask, int n_unique_trials,
                          const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* al_ = cols[emc2col::frq::alpha];
  const double* be_ = cols[emc2col::frq::beta];
  const double* h_  = cols[emc2col::frq::h];
  const double* ta_ = cols[emc2col::frq::tau];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* de_ = cols[emc2col::frq::delta];
  FrqMemo memo;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(al_[r]) || !R_FINITE(t0_[r])) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      const FrqPars& s = memo.get(al_[r], be_[r], h_[r], ta_[r], de_[r]);
      const double log_surv = frq_log_surv_dt(tt, s);
      if (ISNAN(log_surv)) { bad = true; break; }
      logS += log_surv;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// ============================================================
// RDMSWTN adapters
// Column layout: v=0, B=1, A=2, t0=3, s=4, sv=5, mG=6 (guess-clock mean), mK=7 (kill-clock mean)
// ============================================================

// Guarded natural evaluation with a direct log fallback: the natural k = 0
// Wald density/survivor are kept for ordinary cases, and the stable log
// primitives in wald_functions.h take over on underflow, cancellation, or
// CDF saturation.  No natural round trips survive in the raw race kernels.
inline double rdmswtn_k0_logpdf(double tt, double mu, double b, double A, bool posdrift) {
  if (tt <= 0.0) return R_NegInf;
  if (posdrift && mu <= 0.0) return R_NegInf;
  double pdf;
  if (dwald_k0_natural(tt, b, mu, A, pdf))
    return (pdf > 0.0) ? std::log(pdf) : R_NegInf;
  return dwald_k0_log(tt, b, mu, A);
}

inline double rdmswtn_k0_logpdf(double tt, double mu, double b, double A, double s, bool posdrift) {
  if (!(s > 0.0)) return R_NegInf;
  const double inv_s = 1.0 / s;
  return rdmswtn_k0_logpdf(tt, mu * inv_s, b * inv_s, A * inv_s, posdrift);
}

inline double rdmswtn_k0_logsurv(double tt, double mu, double b, double A, bool posdrift) {
  if (tt <= 0.0) return 0.0;
  if (posdrift && mu <= 0.0) return 0.0;
  const double cdf = pwald_k0(tt, b, mu, A);
  const double cl = std::max(0.0, std::min(1.0, cdf));
  if (cl <= 0.0) return 0.0;  // survivor ~ 1: -F below natural resolution
  if (cl >= 1.0 - EMC2_CDF_SAT_MARGIN) {
    // Saturated natural CDF: evaluate the log survivor directly.
    return wald_k0_log_surv(tt, b, mu, A);
  }
  return std::log1p(-cl);
}

inline double rdmswtn_k0_logsurv(double tt, double mu, double b, double A, double s, bool posdrift) {
  if (!(s > 0.0)) return 0.0;
  const double inv_s = 1.0 / s;
  return rdmswtn_k0_logsurv(tt, mu * inv_s, b * inv_s, A * inv_s, posdrift);
}

inline double drdmswtn_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const double lg = (ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[6], ks) : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? erlang_lambda_from_mean(par[7], ks) : 0.0;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  // Pass raw t and t0_val; core functions split EAM (t - t0) from erlang (t).
  const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
  if (dispatch.use_combo) {
    return drdmswtn_local_combo(t, par[0], par[1] + par[2],
                                par[2], par[4], t0_val, par[5],
                                lg, lk, 20, false, ks, pd, omega);
  }
  return drdmswtn(t,
                  par[0],
                  par[1] + par[2],
                  par[2],
                  par[4], t0_val, par[5],
                  dispatch.lambda_g, dispatch.lambda_k,
                  20, false, ks, dispatch.guess, pd, omega);
}

inline double prdmswtn_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const double lg = (ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[6], ks) : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? erlang_lambda_from_mean(par[7], ks) : 0.0;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
  if (dispatch.use_combo) {
    return prdmswtn_local_combo(t, par[0], par[1] + par[2],
                                par[2], par[4], t0_val, par[5],
                                lg, lk, 20, false, ks, pd, omega);
  }
  return prdmswtn(t,
                  par[0],
                  par[1] + par[2],
                  par[2],
                  par[4], t0_val, par[5],
                  dispatch.lambda_g, dispatch.lambda_k,
                  20, false, ks, dispatch.guess, pd, omega);
}

inline void drdmswtn_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_       = cols[emc2col::rdmswtn::v];
  const double* B_       = cols[emc2col::rdmswtn::B];
  const double* A_       = cols[emc2col::rdmswtn::A];
  const double* t0_      = cols[emc2col::rdmswtn::t0];
  const double* s_       = cols[emc2col::rdmswtn::s];
  const double* sv_      = cols[emc2col::rdmswtn::sv];
  const double* lg_ = cols[emc2col::rdmswtn::mG];
  const double* lk_ = cols[emc2col::rdmswtn::mK];
  const double* omega_ = (kill_shape == 3) ? cols[emc2col::rdmswtn::omega] : nullptr;
  const double sv_eps = 1e-10;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i  = t0_[i];
    const double tt    = rt[i] - t0_i;
    const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(kill_shape);
    const double lg = (!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[i], kill_shape);
    const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[i], kill_shape);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (rt[i] <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    double log_pdf;
    const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
    if (dispatch.use_combo) {
      // Pass raw rt and t0; combo function splits EAM vs erlang time.
      log_pdf = drdmswtn_local_combo(rt[i], v_[i], B_[i] + A_[i],
                                     A_[i], s_[i], t0_i, sv_[i],
                                     lg, lk, 20, true, kill_shape, pd, omega);
    } else if (!emc2_isfinite(sv_[i]) || std::fabs(sv_[i]) <= sv_eps) {
      if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0) {
        // No kill, no sv, no erlang: use the closed-form k=0 Wald directly.
        if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
        log_pdf = rdmswtn_k0_logpdf(tt, v_[i], B_[i] + A_[i],
                                    A_[i], s_[i], pd);
      } else {
        // Pass raw rt and t0 to dwald.
        log_pdf = dwald(rt[i], v_[i], B_[i] + A_[i], A_[i],
                        s_[i], t0_i, dispatch.lambda_g, dispatch.lambda_k,
                        true, kill_shape, dispatch.guess, pd, omega);
      }
    } else {
      // Pass raw rt and t0; drdmswtn splits EAM vs erlang time.
      log_pdf = drdmswtn(rt[i], v_[i], B_[i] + A_[i],
                         A_[i], s_[i], t0_i, sv_[i],
                         dispatch.lambda_g, dispatch.lambda_k,
                         20, true, kill_shape, dispatch.guess, pd, omega);
    }
    out[i] = raw_log_value(log_pdf, min_ll, floor_raw);
  }
}

inline void prdmswtn_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_       = cols[emc2col::rdmswtn::v];
  const double* B_       = cols[emc2col::rdmswtn::B];
  const double* A_       = cols[emc2col::rdmswtn::A];
  const double* t0_      = cols[emc2col::rdmswtn::t0];
  const double* s_       = cols[emc2col::rdmswtn::s];
  const double* sv_      = cols[emc2col::rdmswtn::sv];
  const double* lg_ = cols[emc2col::rdmswtn::mG];
  const double* lk_ = cols[emc2col::rdmswtn::mK];
  const double* omega_ = (kill_shape == 3) ? cols[emc2col::rdmswtn::omega] : nullptr;
  const double sv_eps = 1e-10;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i  = t0_[i];
    const double tt    = rt[i] - t0_i;
    const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(kill_shape);
    const double lg = (!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[i], kill_shape);
    const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[i], kill_shape);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = 0.0; continue; }
    if (rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
    if (dispatch.use_combo) {
      const double log_cdf = prdmswtn_local_combo(rt[i], v_[i], B_[i] + A_[i],
                                                  A_[i], s_[i], t0_i, sv_[i],
                                                  lg, lk, 20, true, kill_shape, pd, omega);
      if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
      if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
      out[i] = log1m_exp(log_cdf);
      continue;
    }
    double log_cdf;
    if (!emc2_isfinite(sv_[i]) || std::fabs(sv_[i]) <= sv_eps) {
      if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0) {
        // No kill, no sv, no erlang: use the closed-form k=0 Wald directly.
        if (tt <= 0.0) { out[i] = 0.0; continue; }
        out[i] = rdmswtn_k0_logsurv(tt, v_[i], B_[i] + A_[i],
                                    A_[i], s_[i], pd);
        continue;
      }
      log_cdf = pwald(rt[i], v_[i], B_[i] + A_[i], A_[i],
                      s_[i], t0_i, dispatch.lambda_g, dispatch.lambda_k,
                      true, kill_shape, dispatch.guess, pd, omega);
    } else {
      log_cdf = prdmswtn(rt[i], v_[i], B_[i] + A_[i],
                         A_[i], s_[i], t0_i, sv_[i],
                         dispatch.lambda_g, dispatch.lambda_k,
                         20, true, kill_shape, dispatch.guess, pd, omega);
    }
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

// RDMSWTN: column layout v=0, B=1, A=2, t0=3, s=4, sv=5
inline void rdmswtn_logS_at_t(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const int mode_hint = ctx ? ctx->mode_hint : 0;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_      = cols[emc2col::rdmswtn::v];
  const double* B_      = cols[emc2col::rdmswtn::B];
  const double* A_      = cols[emc2col::rdmswtn::A];
  const double* t0_     = cols[emc2col::rdmswtn::t0];
  const double* s_      = cols[emc2col::rdmswtn::s];
  const double* sv_     = cols[emc2col::rdmswtn::sv];
  const double* lg_ = cols[emc2col::rdmswtn::mG];
  const double* lk_ = cols[emc2col::rdmswtn::mK];
  const double* omega_ = (kill_shape == 3) ? cols[emc2col::rdmswtn::omega] : nullptr;
  const double sv_eps = 1e-10;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double t0_r = t0_[r];
      const double tt = t - t0_r;
      const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[r])) :
                           erlang_omega_for_shape(kill_shape);
      const double lg = (!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[r], kill_shape);
      const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[r], kill_shape);
      const bool erl = (lg > 1e-12 || lk > 1e-12);
      if (tt <= 0.0) {
        if (!erl) continue;  // EAM not started, no erlang → logS += 0
        // EAM not started but erlang running: log-survivor = erlang_log_surv(t, ...)
        const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
        if (dispatch.use_combo) {
          const double log_cdf = prdmswtn_local_combo(t, v_[r], B_[r] + A_[r],
                                                      A_[r], s_[r], t0_r, sv_[r],
                                                      lg, lk, 20, true, kill_shape, pd, omega);
          if (!R_FINITE(log_cdf) || log_cdf >= 0.0) { bad = true; break; }
          logS += log1m_exp(log_cdf);
        } else {
          double log_cdf;
          if (mode_hint == 1 || (!emc2_isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps)) {
            log_cdf = pwald(t,
                            v_[r],
                            B_[r] + A_[r],
                            A_[r],
                            s_[r],
                            t0_r,
                            dispatch.lambda_g, dispatch.lambda_k,
                            true, kill_shape, dispatch.guess, pd, omega);
          } else {
            log_cdf = prdmswtn(t,
                               v_[r],
                               B_[r] + A_[r],
                               A_[r],
                               s_[r],
                               t0_r,
                               sv_[r],
                               dispatch.lambda_g, dispatch.lambda_k,
                               20, true, kill_shape, dispatch.guess, pd, omega);
          }
          if (!R_FINITE(log_cdf) || log_cdf >= 0.0) { bad = true; break; }
          logS += log1m_exp(log_cdf);
        }
        continue;
      }
      double log_cdf;
      const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
      if (dispatch.use_combo) {
        log_cdf = prdmswtn_local_combo(t, v_[r], B_[r] + A_[r],
                                       A_[r], s_[r], t0_r, sv_[r],
                                       lg, lk, 20, true, kill_shape, pd, omega);
        if (!R_FINITE(log_cdf) || log_cdf >= 0.0) { bad = true; break; }
        logS += log1m_exp(log_cdf);
        continue;
      }
      if (mode_hint == 1) {
        if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0 &&
            (!emc2_isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps)) {
          logS += rdmswtn_k0_logsurv(tt, v_[r], B_[r] + A_[r],
                                     A_[r], s_[r], pd);
          continue;
        } else {
          log_cdf = pwald(t,
                          v_[r],
                          B_[r] + A_[r],
                          A_[r],
                          s_[r],
                          t0_r,
                          dispatch.lambda_g, dispatch.lambda_k,
                          true, kill_shape, dispatch.guess, pd, omega);
        }
      } else if (mode_hint == 2) {
        log_cdf = prdmswtn(t,
                           v_[r],
                           B_[r] + A_[r],
                           A_[r],
                           s_[r],
                           t0_r,
                           sv_[r],
                           dispatch.lambda_g, dispatch.lambda_k,
                           20, true, kill_shape, dispatch.guess, pd, omega);
      } else if (!emc2_isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps) {
        if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0) {
          logS += rdmswtn_k0_logsurv(tt, v_[r], B_[r] + A_[r],
                                     A_[r], s_[r], pd);
          continue;
        } else {
          log_cdf = pwald(t,
                          v_[r],
                          B_[r] + A_[r],
                          A_[r],
                          s_[r],
                          t0_r,
                          dispatch.lambda_g, dispatch.lambda_k,
                          true, kill_shape, dispatch.guess, pd, omega);
        }
      } else {
        log_cdf = prdmswtn(t,
                           v_[r],
                           B_[r] + A_[r],
                           A_[r],
                           s_[r],
                           t0_r,
                           sv_[r],
                           dispatch.lambda_g, dispatch.lambda_k,
                           20, true, kill_shape, dispatch.guess, pd, omega);
      }
      if (log_cdf >= 0.0) { bad = true; break; }
      if (R_FINITE(log_cdf)) logS += log1m_exp(log_cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// ============================================================
// RDMSWTN_TT adapters
// Column layout: v=0, B=1, A=2, t0=3, s=4, sv=5, tau=6
// ============================================================

inline double drdmswtn_tt_scalar(double t, const double* par, void* ctx_) {
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

inline double prdmswtn_tt_scalar(double t, const double* par, void* ctx_) {
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

inline void drdmswtn_tt_raw(
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

inline void prdmswtn_tt_raw(
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

inline void rdmswtn_tt_logS_at_t(
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

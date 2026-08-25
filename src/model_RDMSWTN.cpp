#include <Rcpp.h>
#include "model_RDMSWTN.h"
#include "model_RDM.h"
#include "wald_functions.h"
#include "timer_helpers.h"
#include "utility_functions.h"
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
// ============================================================
// RDMSWTN adapters
// Column layout: v=0, B=1, A=2, t0=3, s=4, sv=5, mG=6 (guess-clock mean), mK=7 (kill-clock mean)
// ============================================================

// Guarded natural evaluation with a direct log fallback: the natural k = 0
// Wald density/survivor are kept for ordinary cases, and the stable log
// primitives in wald_functions.h take over on underflow, cancellation, or
// CDF saturation.  No natural round trips survive in the raw race kernels.
static inline double rdmswtn_k0_logpdf(double tt, double mu, double b, double A, bool posdrift) {
  if (tt <= 0.0) return R_NegInf;
  if (posdrift && mu <= 0.0) return R_NegInf;
  double pdf;
  if (dwald_k0_natural(tt, b, mu, A, pdf))
    return (pdf > 0.0) ? std::log(pdf) : R_NegInf;
  return dwald_k0_log(tt, b, mu, A);
}

static inline double rdmswtn_k0_logpdf(double tt, double mu, double b, double A, double s, bool posdrift) {
  if (!(s > 0.0)) return R_NegInf;
  const double inv_s = 1.0 / s;
  return rdmswtn_k0_logpdf(tt, mu * inv_s, b * inv_s, A * inv_s, posdrift);
}

static inline double rdmswtn_k0_logsurv(double tt, double mu, double b, double A, bool posdrift) {
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

static inline double rdmswtn_k0_logsurv(double tt, double mu, double b, double A, double s, bool posdrift) {
  if (!(s > 0.0)) return 0.0;
  const double inv_s = 1.0 / s;
  return rdmswtn_k0_logsurv(tt, mu * inv_s, b * inv_s, A * inv_s, posdrift);
}

double drdmswtn_scalar(double t, const double* par, void* ctx_) {
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

double prdmswtn_scalar(double t, const double* par, void* ctx_) {
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

void drdmswtn_raw(const double* rt, const double* const* cols, int n_rows,
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

void prdmswtn_raw(const double* rt, const double* const* cols, int n_rows,
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
void rdmswtn_logS_at_t(double t, const double* const* cols,
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

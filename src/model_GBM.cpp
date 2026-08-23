#include <RcppArmadillo.h>
#include "model_GBM.h"
#include "model_RDM.h"
#include "timer_helpers.h"
#include "col_registry.h"

// GBM: column layout v=0, B=1, A=2, t0=3, s=4, mG=5 (guess-clock mean), mK=6 (kill-clock mean)
double drdmgbm_scalar(double t, const double* par, void* ctx_) {
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

double prdmgbm_scalar(double t, const double* par, void* ctx_) {
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

void drdmgbm_raw(const double* rt, const double* const* cols, int n_rows,
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

void prdmgbm_raw(const double* rt, const double* const* cols, int n_rows,
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
void rdmgbm_logS_at_t(double t, const double* const* cols,
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

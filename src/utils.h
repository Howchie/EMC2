#ifndef utils_h
#define utils_h

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <unordered_map>
#include "utility_functions.h"
#include "model_RDM.h"
#include "model_LBA.h"

using namespace Rcpp;

// Function pointer types for race model PDF/CDF adapters
typedef Rcpp::NumericVector (*RacePdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
                                          Rcpp::LogicalVector winner,
                                          Rcpp::LogicalVector is_ok,
                                          void* model_specific_context);
typedef Rcpp::NumericVector (*RaceCdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
                                          Rcpp::LogicalVector winner,
                                          Rcpp::LogicalVector is_ok,
                                          void* model_specific_context);

// Scalar (single-RT / single-accumulator) helpers for GSL integration
typedef double (*RacePdf1Fun)(double rt, const double* par, void* model_specific_context);
typedef double (*RaceCdf1Fun)(double rt, const double* par, void* model_specific_context);

// Fast-path batch raw kernel: writes log-density or log-survivor to pre-allocated buffer
typedef void (*RaceRawFun)(const double* rt, const double* pars_cm, int n_rows,
                           const int* mask, const int* isok,
                           double* out, double min_ll, void* ctx_);

// Batch log-survivor at scalar t for truncation normalisation
typedef void (*RaceLogSAtTFun)(double t, const double* pars_cm,
                               int n_rows_total, int n_lR, int n_par,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out);

struct gsl_race_params {
  const Rcpp::NumericMatrix* p_trial;
  Rcpp::LogicalVector winner;
  Rcpp::LogicalVector isok;
  RacePdfFun model_dfun;
  RaceCdfFun model_pfun;
  int n_lR;
  void* model_specific_context;
};

struct gsl_race_params_scalar {
  const double* pars;   // row-major, length n_lR * n_par
  int n_lR;
  int n_par;
  int winner_idx0;      // 0-based
  const int* isok;      // length n_lR, 0/1
  RacePdf1Fun pdf1;
  RaceCdf1Fun cdf1;
  void* ctx;
};

double gsl_f_race_scalar(double t, void* p);

// --------------------------------------------------------------------------
// Context object for race models, holding metadata and switches.
// --------------------------------------------------------------------------
struct ContextForRaceModels {
    double min_lik_for_pdf = 1e-10;
    bool use_posdrift = true;
    bool gng = false;
    int t0_index = -1;
    int lambda_g_index = -1;
    int lambda_k_index = -1;
    // For models with infinite tails or defective upper mass (like LBA with sv).
    bool defective_upper_tail = false;
    // Optional per-particle fast-kernel hint:
    // 0 = auto detect in kernel; 1 = zero-variability branch; 2 = nonzero branch.
    int mode_hint = 0;
    // Kill-process shape: 1 = exponential, 2 = Erlang-2, 3 = omega mixture.
    int kill_shape = 1;
    int erlang_omega_index = -1;

    // 2x2 Model Flags
    bool is_local_guess = false;
    bool is_global_kill = false;
    bool is_local_kill = false;
    bool is_local_kill_guess = false;

    // Per-particle switch: disable kill bookkeeping when lambda is zero
    bool kill_active = true;

    // Raw race kernels historically floor log probabilities to min_ll.  Timed
    // race mixtures need unfloored components before the final log_sum_exp.
    bool floor_raw_log_lik = true;

    // When false, adapters ignore lambda_k (for global kill races).
    // Set to false in global-kill models, and true in single-accumulator analytic paths.
    bool apply_lk_to_racers = true;

    // Tri-state caches for optional accumulator levels:
    // -2 = unresolved (detect from data once), -1 = absent, >0 = factor code.
    int time_code = -2;
    int nogo_code = -2;

    bool has_global_kill() const {
      return is_global_kill && kill_active;
    }
};

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

// Scalar adapters (single-RT, single-parameter-row) used by GSL integration.
inline double dlba_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[4];
  if (tt <= 0.0) return 0.0;
  // par: v=0, sv=1, B=2, A=3, t0=4 (LBA uses fixed scale s=1)
  return dlba_norm(tt, par[3], par[2] + par[3], par[0], par[1], ctx->use_posdrift);
}

inline double plba_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[4];
  if (tt <= 0.0) return 0.0;
  return plba_norm(tt, par[3], par[2] + par[3], par[0], par[1], ctx->use_posdrift);
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

// GBM: column layout v=0, B=1, A=2, t0=3, s=4, lambda_g=5, lambda_k=6
inline double drdmgbm_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const double lg = (ctx && ctx->kill_active) ? par[5] : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? par[6] : 0.0;
  const int ks = ctx ? ctx->kill_shape : 1;
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
  const double lg = (ctx && ctx->kill_active) ? par[5] : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? par[6] : 0.0;
  const int ks = ctx ? ctx->kill_shape : 1;
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

// RDM: column layout v=0, B=1, A=2, t0=3, s=4
inline void drdm_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = pars_cm + 0 * n_rows;
  const double* B_  = pars_cm + 1 * n_rows;
  const double* A_  = pars_cm + 2 * n_rows;
  const double* t0_ = pars_cm + 3 * n_rows;
  const double* s_  = pars_cm + 4 * n_rows;
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

inline void prdm_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = pars_cm + 0 * n_rows;
  const double* B_  = pars_cm + 1 * n_rows;
  const double* A_  = pars_cm + 2 * n_rows;
  const double* t0_ = pars_cm + 3 * n_rows;
  const double* s_  = pars_cm + 4 * n_rows;
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
inline void rdm_logS_at_t(double t, const double* pars_cm,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* v_  = pars_cm + 0 * n_rows_total;
  const double* B_  = pars_cm + 1 * n_rows_total;
  const double* A_  = pars_cm + 2 * n_rows_total;
  const double* t0_ = pars_cm + 3 * n_rows_total;
  const double* s_  = pars_cm + 4 * n_rows_total;
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

inline void drdmgbm_raw(const double* rt, const double* pars_cm, int n_rows,
                        const int* mask, const int* isok,
                        double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool global_kill = ctx ? ctx->has_global_kill() : false;
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const double* v_        = pars_cm + 0 * n_rows;
  const double* B_        = pars_cm + 1 * n_rows;
  const double* A_        = pars_cm + 2 * n_rows;
  const double* t0_       = pars_cm + 3 * n_rows;
  const double* s_        = pars_cm + 4 * n_rows;
  const double* lambda_g_ = pars_cm + 5 * n_rows;
  const double* lambda_k_ = pars_cm + 6 * n_rows;
  const double* omega_    = (kill_shape == 3) ? pars_cm + 7 * n_rows : nullptr;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double lg = (!ctx->kill_active) ? 0.0 : lambda_g_[i];
    const double lk = (global_kill || !ctx->kill_active) ? 0.0 : lambda_k_[i];
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

inline void prdmgbm_raw(const double* rt, const double* pars_cm, int n_rows,
                        const int* mask, const int* isok,
                        double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool global_kill = ctx ? ctx->has_global_kill() : false;
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const double* v_        = pars_cm + 0 * n_rows;
  const double* B_        = pars_cm + 1 * n_rows;
  const double* A_        = pars_cm + 2 * n_rows;
  const double* t0_       = pars_cm + 3 * n_rows;
  const double* s_        = pars_cm + 4 * n_rows;
  const double* lambda_g_ = pars_cm + 5 * n_rows;
  const double* lambda_k_ = pars_cm + 6 * n_rows;
  const double* omega_    = (kill_shape == 3) ? pars_cm + 7 * n_rows : nullptr;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double lg = (!ctx->kill_active) ? 0.0 : lambda_g_[i];
    const double lk = (global_kill || !ctx->kill_active) ? 0.0 : lambda_k_[i];
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

// GBM: column layout v=0, B=1, A=2, t0=3, s=4, lambda_g=5, lambda_k=6
inline void rdmgbm_logS_at_t(double t, const double* pars_cm,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool global_kill = ctx ? ctx->has_global_kill() : false;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* v_        = pars_cm + 0 * n_rows_total;
  const double* B_        = pars_cm + 1 * n_rows_total;
  const double* A_        = pars_cm + 2 * n_rows_total;
  const double* t0_       = pars_cm + 3 * n_rows_total;
  const double* s_        = pars_cm + 4 * n_rows_total;
  const double* lambda_g_ = pars_cm + 5 * n_rows_total;
  const double* lambda_k_ = pars_cm + 6 * n_rows_total;
  const double* omega_    = (ks == 3) ? pars_cm + 7 * n_rows_total : nullptr;

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
      const double lg = (!ctx->kill_active) ? 0.0 : lambda_g_[r];
      const double lk = (global_kill || !ctx->kill_active) ? 0.0 : lambda_k_[r];
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

// LNR: column layout m=0, s=1, t0=2
inline void dlnr_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* m_  = pars_cm + 0 * n_rows;
  const double* s_  = pars_cm + 1 * n_rows;
  const double* t0_ = pars_cm + 2 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(m_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double pdf = dlnorm_std(tt, m_[i], s_[i], false);
    out[i] = (pdf > 0.0 && emc2_isfinite(pdf)) ? raw_log_value(std::log(pdf), min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void plnr_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* /*ctx_*/) {
  const double* m_  = pars_cm + 0 * n_rows;
  const double* s_  = pars_cm + 1 * n_rows;
  const double* t0_ = pars_cm + 2 * n_rows;
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

inline void lnr_logS_at_t(double t, const double* pars_cm,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* m_  = pars_cm + 0 * n_rows_total;
  const double* s_  = pars_cm + 1 * n_rows_total;
  const double* t0_ = pars_cm + 2 * n_rows_total;
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

// RGAMMA: column layout lambda=0, shape=1, shift=2
inline double drgamma_scalar(double t, const double* par, void* ctx_) {
  (void)ctx_;
  if (R_IsNA(par[0]) || R_IsNA(par[1]) || R_IsNA(par[2])) return 0.0;
  if (par[0] <= 0.0 || par[1] <= 0.0) return 0.0;
  const double tt = t - par[2];
  if (tt <= 0.0) return 0.0;
  return R::dgamma(tt, par[1], 1.0 / par[0], false);
}

inline double prgamma_scalar(double t, const double* par, void* ctx_) {
  (void)ctx_;
  if (R_IsNA(par[0]) || R_IsNA(par[1]) || R_IsNA(par[2])) return 0.0;
  if (par[0] <= 0.0 || par[1] <= 0.0) return 0.0;
  const double tt = t - par[2];
  if (tt <= 0.0) return 0.0;
  return R::pgamma(tt, par[1], 1.0 / par[0], true, false);
}

inline void drgamma_raw(const double* rt, const double* pars_cm, int n_rows,
                        const int* mask, const int* isok,
                        double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* lambda_ = pars_cm + 0 * n_rows;
  const double* shape_  = pars_cm + 1 * n_rows;
  const double* shift_  = pars_cm + 2 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(lambda_[i]) || R_IsNA(shape_[i]) || R_IsNA(shift_[i]) ||
        !isok[i] || lambda_[i] <= 0.0 || shape_[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - shift_[i];
    if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double pdf = R::dgamma(tt, shape_[i], 1.0 / lambda_[i], false);
    out[i] = (pdf > 0.0 && emc2_isfinite(pdf)) ? raw_log_value(std::log(pdf), min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void prgamma_raw(const double* rt, const double* pars_cm, int n_rows,
                        const int* mask, const int* isok,
                        double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* lambda_ = pars_cm + 0 * n_rows;
  const double* shape_  = pars_cm + 1 * n_rows;
  const double* shift_  = pars_cm + 2 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(lambda_[i]) || R_IsNA(shape_[i]) || R_IsNA(shift_[i]) ||
        !isok[i] || lambda_[i] <= 0.0 || shape_[i] <= 0.0) {
      out[i] = 0.0;
      continue;
    }
    const double tt = rt[i] - shift_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double cdf = R::pgamma(tt, shape_[i], 1.0 / lambda_[i], true, false);
    if (cdf >= 1.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

inline void rgamma_logS_at_t(double t, const double* pars_cm,
                             int n_rows_total, int n_lR, int /*n_par*/,
                             const int* trunc_mask, int n_unique_trials,
                             const int* isok_all, void* ctx_, double* logS_out) {
  (void)ctx_;
  const double* lambda_ = pars_cm + 0 * n_rows_total;
  const double* shape_  = pars_cm + 1 * n_rows_total;
  const double* shift_  = pars_cm + 2 * n_rows_total;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(lambda_[r]) || R_IsNA(shape_[r]) || R_IsNA(shift_[r]) ||
          lambda_[r] <= 0.0 || shape_[r] <= 0.0) { bad = true; break; }
      const double tt = t - shift_[r];
      if (tt <= 0.0) continue;
      const double cdf = R::pgamma(tt, shape_[r], 1.0 / lambda_[r], true, false);
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// LBA: column layout v=0, sv=1, B=2, A=3, t0=4 (fixed scale s=1)
inline void dlba_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_ = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double pdf = dlba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], pd);
    out[i] = (pdf > 0.0 && emc2_isfinite(pdf)) ? raw_log_value(std::log(pdf), min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void plba_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_ = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double cdf = plba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], pd);
    if (cdf >= 1.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

inline void lba_logS_at_t(double t, const double* pars_cm,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_  = pars_cm + 0 * n_rows_total;
  const double* sv_ = pars_cm + 1 * n_rows_total;
  const double* B_  = pars_cm + 2 * n_rows_total;
  const double* A_  = pars_cm + 3 * n_rows_total;
  const double* t0_ = pars_cm + 4 * n_rows_total;
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
      const double cdf = plba_norm(tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], pd);
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// ============================================================
// BAwL (Ballistic Accumulator with Leak + killing/guessing) adapters
// Column layout: v=0, sv=1, B=2, A=3, t0=4, k=5, lambda_g=6, lambda_k=7
// ============================================================

inline double dbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[4];
  const double tt = t - t0_val;
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const double lg = (ctx && ctx->kill_active) ? par[6] : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? par[7] : 0.0;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  // Pass raw t and t0_val; core function splits EAM (t - t0) from erlang (t).
  return dkilledleakyba_norm(
    t, par[0], par[2] + par[3], par[3], par[1], t0_val, par[5], lg, lk,
    ctx->use_posdrift, false, ks, local_guess, omega
  );
}

inline double pbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[4];
  const double tt = t - t0_val;
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const double lg = (ctx && ctx->kill_active) ? par[6] : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? par[7] : 0.0;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  return pkilledleakyba_norm(
    t, par[0], par[2] + par[3], par[3], par[1], t0_val, par[5], lg, lk,
    ctx->use_posdrift, false, ks, local_guess, omega
  );
}

inline void dbawl_raw(const double* rt, const double* pars_cm, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
  const bool local_guess = ctx->is_local_guess || ctx->is_local_kill_guess;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_  = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  const double* k_  = pars_cm + 5 * n_rows;
  const double* lg_ = pars_cm + 6 * n_rows;
  const double* lk_ = pars_cm + 7 * n_rows;
  const double* omega_ = (ks == 3) ? pars_cm + 8 * n_rows : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double lg = (!ctx->kill_active) ? 0.0 : lg_[i];
    const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : lk_[i];
    const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(ks);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (rt[i] <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    // Pass raw rt and t0; core function uses t0 to split EAM vs erlang time.
    const double pdf = dkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, k_[i], lg, lk,
      pd, false, ks, local_guess, omega
    );
    out[i] = (pdf > 0.0 && emc2_isfinite(pdf)) ? raw_log_value(std::log(pdf), min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

inline void pbawl_raw(const double* rt, const double* pars_cm, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
  const bool local_guess = ctx->is_local_guess || ctx->is_local_kill_guess;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_  = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  const double* k_  = pars_cm + 5 * n_rows;
  const double* lg_ = pars_cm + 6 * n_rows;
  const double* lk_ = pars_cm + 7 * n_rows;
  const double* omega_ = (ks == 3) ? pars_cm + 8 * n_rows : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double lg = (!ctx->kill_active) ? 0.0 : lg_[i];
    const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : lk_[i];
    const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(ks);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = 0.0; continue; }
    if (rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const double log_cdf = pkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, k_[i], lg, lk,
      pd, true, ks, local_guess, omega
    );
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

// BAwL: column layout v=0, sv=1, B=2, A=3, t0=4, k=5, lambda_g=6, lambda_k=7
inline void bawl_logS_at_t(double t, const double* pars_cm,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd          = ctx->use_posdrift;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* v_  = pars_cm + 0 * n_rows_total;
  const double* sv_ = pars_cm + 1 * n_rows_total;
  const double* B_  = pars_cm + 2 * n_rows_total;
  const double* A_  = pars_cm + 3 * n_rows_total;
  const double* t0_ = pars_cm + 4 * n_rows_total;
  const double* k_  = pars_cm + 5 * n_rows_total;
  const double* lg_ = pars_cm + 6 * n_rows_total;
  const double* lk_ = pars_cm + 7 * n_rows_total;
  const double* omega_ = (ks == 3) ? pars_cm + 8 * n_rows_total : nullptr;
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
      const double lg = (!ctx->kill_active) ? 0.0 : lg_[r];
      const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : lk_[r];
      const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[r])) :
                           erlang_omega_for_shape(ks);
      const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
      const bool erl = (lg > 1e-12 || lk > 1e-12);
      if (tt <= 0.0) {
        if (!erl) continue;  // EAM not started, no erlang → log-survivor += 0
        // EAM not started but erlang running: log-survivor = erlang_log_surv(t, ...)
        if (local_guess && lg > 1e-12 && lk > 1e-12) {
          logS += erlang_log_surv(t, lg, ks, omega) + erlang_log_surv(t, lk, ks, omega);
        } else {
          const bool is_g = ctx->is_local_guess || ctx->is_local_kill_guess;
          const double lam = (is_g ? lg : 0.0) + lk;
          if (lam > 1e-12) logS += erlang_log_surv(t, lam, ks, omega);
        }
        continue;
      }
      // Both EAM and erlang contribute; pass raw t and t0_r.
      const double cdf = pkilledleakyba_norm(
        t, v_[r], B_[r] + A_[r], A_[r], sv_[r], t0_r, k_[r], lg, lk,
        pd, false, ks, local_guess, omega
      );
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// ============================================================
// RDMSWTN adapters
// Column layout: v=0, B=1, A=2, t0=3, s=4, sv=5, lambda_g=6, lambda_k=7
// ============================================================

inline double rdmswtn_k0_logpdf(double tt, double mu, double b, double A, bool posdrift) {
  if (tt <= 0.0) return R_NegInf;
  if (posdrift && mu <= 0.0) return R_NegInf;
  const double pdf = dwald_k0(tt, b, mu, A);
  if (!(pdf > 0.0) || !emc2_isfinite(pdf)) return R_NegInf;
  return std::log(pdf);
}

inline double rdmswtn_k0_logsurv(double tt, double mu, double b, double A, bool posdrift) {
  if (tt <= 0.0) return 0.0;
  if (posdrift && mu <= 0.0) return 0.0;
  const double cdf = pwald_k0(tt, b, mu, A);
  const double cl = std::max(0.0, std::min(1.0, cdf));
  if (cl <= 0.0) return 0.0;
  if (cl >= 1.0) return R_NegInf;
  return std::log1p(-cl);
}

inline double drdmswtn_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const double inv_s = 1.0 / par[4];
  const int ks = ctx ? ctx->kill_shape : 1;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const double lg = (ctx && ctx->kill_active) ? par[6] : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? par[7] : 0.0;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  // Pass raw t and t0_val; core functions split EAM (t - t0) from erlang (t).
  const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
  if (dispatch.use_combo) {
    return drdmswtn_local_combo(t, par[0] * inv_s, (par[1] + par[2]) * inv_s,
                                par[2] * inv_s, 1.0, t0_val, par[5] * inv_s,
                                lg, lk, 20, false, ks, pd, omega);
  }
  return drdmswtn(t,
                  par[0] * inv_s,
                  (par[1] + par[2]) * inv_s,
                  par[2] * inv_s,
                  1.0, t0_val, par[5] * inv_s,
                  dispatch.lambda_g, dispatch.lambda_k,
                  20, false, ks, dispatch.guess, pd, omega);
}

inline double prdmswtn_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double t0_val = par[3];
  const double tt = t - t0_val;
  const double inv_s = 1.0 / par[4];
  const int ks = ctx ? ctx->kill_shape : 1;
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const double lg = (ctx && ctx->kill_active) ? par[6] : 0.0;
  const double lk = (ctx && ctx->kill_active && ctx->apply_lk_to_racers) ? par[7] : 0.0;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
  if (dispatch.use_combo) {
    return prdmswtn_local_combo(t, par[0] * inv_s, (par[1] + par[2]) * inv_s,
                                par[2] * inv_s, 1.0, t0_val, par[5] * inv_s,
                                lg, lk, 20, false, ks, pd, omega);
  }
  return prdmswtn(t,
                  par[0] * inv_s,
                  (par[1] + par[2]) * inv_s,
                  par[2] * inv_s,
                  1.0, t0_val, par[5] * inv_s,
                  dispatch.lambda_g, dispatch.lambda_k,
                  20, false, ks, dispatch.guess, pd, omega);
}

inline void drdmswtn_raw(const double* rt, const double* pars_cm, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_       = pars_cm + 0 * n_rows;
  const double* B_       = pars_cm + 1 * n_rows;
  const double* A_       = pars_cm + 2 * n_rows;
  const double* t0_      = pars_cm + 3 * n_rows;
  const double* s_       = pars_cm + 4 * n_rows;
  const double* sv_      = pars_cm + 5 * n_rows;
  const double* lg_ = pars_cm + 6 * n_rows;
  const double* lk_ = pars_cm + 7 * n_rows;
  const double* omega_ = (kill_shape == 3) ? pars_cm + 8 * n_rows : nullptr;
  const double sv_eps = 1e-10;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i  = t0_[i];
    const double tt    = rt[i] - t0_i;
    const double inv_s = 1.0 / s_[i];
    const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(kill_shape);
    const double lg = (!ctx->kill_active) ? 0.0 : lg_[i];
    const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : lk_[i];
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (rt[i] <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    double log_pdf;
    const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
    if (dispatch.use_combo) {
      // Pass raw rt and t0; combo function splits EAM vs erlang time.
      log_pdf = drdmswtn_local_combo(rt[i], v_[i] * inv_s, (B_[i] + A_[i]) * inv_s,
                                     A_[i] * inv_s, 1.0, t0_i, sv_[i] * inv_s,
                                     lg, lk, 20, true, kill_shape, pd, omega);
    } else if (!emc2_isfinite(sv_[i]) || std::fabs(sv_[i]) <= sv_eps) {
      if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0) {
        // No kill, no sv, no erlang: use the closed-form k=0 Wald directly.
        if (tt <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
        log_pdf = rdmswtn_k0_logpdf(tt, v_[i] * inv_s, (B_[i] + A_[i]) * inv_s,
                                    A_[i] * inv_s, pd);
      } else {
        // Pass raw rt and t0 to dwald.
        log_pdf = dwald(rt[i], v_[i] * inv_s, (B_[i] + A_[i]) * inv_s, A_[i] * inv_s,
                        1.0, t0_i, dispatch.lambda_g, dispatch.lambda_k,
                        true, kill_shape, dispatch.guess, pd, omega);
      }
    } else {
      // Pass raw rt and t0; drdmswtn splits EAM vs erlang time.
      log_pdf = drdmswtn(rt[i], v_[i] * inv_s, (B_[i] + A_[i]) * inv_s,
                         A_[i] * inv_s, 1.0, t0_i, sv_[i] * inv_s,
                         dispatch.lambda_g, dispatch.lambda_k,
                         20, true, kill_shape, dispatch.guess, pd, omega);
    }
    out[i] = raw_log_value(log_pdf, min_ll, floor_raw);
  }
}

inline void prdmswtn_raw(const double* rt, const double* pars_cm, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_       = pars_cm + 0 * n_rows;
  const double* B_       = pars_cm + 1 * n_rows;
  const double* A_       = pars_cm + 2 * n_rows;
  const double* t0_      = pars_cm + 3 * n_rows;
  const double* s_       = pars_cm + 4 * n_rows;
  const double* sv_      = pars_cm + 5 * n_rows;
  const double* lg_ = pars_cm + 6 * n_rows;
  const double* lk_ = pars_cm + 7 * n_rows;
  const double* omega_ = (kill_shape == 3) ? pars_cm + 8 * n_rows : nullptr;
  const double sv_eps = 1e-10;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i  = t0_[i];
    const double tt    = rt[i] - t0_i;
    const double inv_s = 1.0 / s_[i];
    const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(kill_shape);
    const double lg = (!ctx->kill_active) ? 0.0 : lg_[i];
    const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : lk_[i];
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = 0.0; continue; }
    if (rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
    if (dispatch.use_combo) {
      const double log_cdf = prdmswtn_local_combo(rt[i], v_[i] * inv_s, (B_[i] + A_[i]) * inv_s,
                                                  A_[i] * inv_s, 1.0, t0_i, sv_[i] * inv_s,
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
        out[i] = rdmswtn_k0_logsurv(tt, v_[i] * inv_s, (B_[i] + A_[i]) * inv_s,
                                    A_[i] * inv_s, pd);
        continue;
      }
      log_cdf = pwald(rt[i], v_[i] * inv_s, (B_[i] + A_[i]) * inv_s, A_[i] * inv_s,
                      1.0, t0_i, dispatch.lambda_g, dispatch.lambda_k,
                      true, kill_shape, dispatch.guess, pd, omega);
    } else {
      log_cdf = prdmswtn(rt[i], v_[i] * inv_s, (B_[i] + A_[i]) * inv_s,
                         A_[i] * inv_s, 1.0, t0_i, sv_[i] * inv_s,
                         dispatch.lambda_g, dispatch.lambda_k,
                         20, true, kill_shape, dispatch.guess, pd, omega);
    }
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

// RDMSWTN: column layout v=0, B=1, A=2, t0=3, s=4, sv=5
inline void rdmswtn_logS_at_t(double t, const double* pars_cm,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int kill_shape   = ctx ? ctx->kill_shape : 1;
  const int mode_hint = ctx ? ctx->mode_hint : 0;
  const bool pd = ctx ? ctx->use_posdrift : true;
  const double* v_      = pars_cm + 0 * n_rows_total;
  const double* B_      = pars_cm + 1 * n_rows_total;
  const double* A_      = pars_cm + 2 * n_rows_total;
  const double* t0_     = pars_cm + 3 * n_rows_total;
  const double* s_      = pars_cm + 4 * n_rows_total;
  const double* sv_     = pars_cm + 5 * n_rows_total;
  const double* lg_ = pars_cm + 6 * n_rows_total;
  const double* lk_ = pars_cm + 7 * n_rows_total;
  const double* omega_ = (kill_shape == 3) ? pars_cm + 8 * n_rows_total : nullptr;
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
      const double inv_s = 1.0 / s_[r];
      const double omega = (kill_shape == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[r])) :
                           erlang_omega_for_shape(kill_shape);
      const double lg = (!ctx->kill_active) ? 0.0 : lg_[r];
      const double lk = (!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : lk_[r];
      const bool erl = (lg > 1e-12 || lk > 1e-12);
      if (tt <= 0.0) {
        if (!erl) continue;  // EAM not started, no erlang → logS += 0
        // EAM not started but erlang running: log-survivor = erlang_log_surv(t, ...)
        const TimedLambdaDispatch dispatch = timed_lambda_dispatch(ctx, lg, lk);
        if (dispatch.use_combo) {
          const double log_cdf = prdmswtn_local_combo(t, v_[r] * inv_s, (B_[r] + A_[r]) * inv_s,
                                                      A_[r] * inv_s, 1.0, t0_r, sv_[r] * inv_s,
                                                      lg, lk, 20, true, kill_shape, pd, omega);
          if (!R_FINITE(log_cdf) || log_cdf >= 0.0) { bad = true; break; }
          logS += log1m_exp(log_cdf);
        } else {
          double log_cdf;
          if (mode_hint == 1 || (!emc2_isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps)) {
            log_cdf = pwald(t,
                            v_[r] * inv_s,
                            (B_[r] + A_[r]) * inv_s,
                            A_[r] * inv_s,
                            1.0,
                            t0_r,
                            dispatch.lambda_g, dispatch.lambda_k,
                            true, kill_shape, dispatch.guess, pd, omega);
          } else {
            log_cdf = prdmswtn(t,
                               v_[r]  * inv_s,
                               (B_[r] + A_[r]) * inv_s,
                               A_[r]  * inv_s,
                               1.0,
                               t0_r,
                               sv_[r] * inv_s,
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
        log_cdf = prdmswtn_local_combo(t, v_[r] * inv_s, (B_[r] + A_[r]) * inv_s,
                                       A_[r] * inv_s, 1.0, t0_r, sv_[r] * inv_s,
                                       lg, lk, 20, true, kill_shape, pd, omega);
        if (!R_FINITE(log_cdf) || log_cdf >= 0.0) { bad = true; break; }
        logS += log1m_exp(log_cdf);
        continue;
      }
      if (mode_hint == 1) {
        if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0 &&
            (!emc2_isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps)) {
          logS += rdmswtn_k0_logsurv(tt, v_[r] * inv_s, (B_[r] + A_[r]) * inv_s,
                                     A_[r] * inv_s, pd);
          continue;
        } else {
          log_cdf = pwald(t,
                          v_[r] * inv_s,
                          (B_[r] + A_[r]) * inv_s,
                          A_[r] * inv_s,
                          1.0,
                          t0_r,
                          dispatch.lambda_g, dispatch.lambda_k,
                          true, kill_shape, dispatch.guess, pd, omega);
        }
      } else if (mode_hint == 2) {
        log_cdf = prdmswtn(t,
                           v_[r]  * inv_s,
                           (B_[r] + A_[r]) * inv_s,
                           A_[r]  * inv_s,
                           1.0,
                           t0_r,
                           sv_[r] * inv_s,
                           dispatch.lambda_g, dispatch.lambda_k,
                           20, true, kill_shape, dispatch.guess, pd, omega);
      } else if (!emc2_isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps) {
        if (dispatch.lambda_g <= 0.0 && dispatch.lambda_k <= 0.0) {
          logS += rdmswtn_k0_logsurv(tt, v_[r] * inv_s, (B_[r] + A_[r]) * inv_s,
                                     A_[r] * inv_s, pd);
          continue;
        } else {
          log_cdf = pwald(t,
                          v_[r] * inv_s,
                          (B_[r] + A_[r]) * inv_s,
                          A_[r] * inv_s,
                          1.0,
                          t0_r,
                          dispatch.lambda_g, dispatch.lambda_k,
                          true, kill_shape, dispatch.guess, pd, omega);
        }
      } else {
        log_cdf = prdmswtn(t,
                           v_[r]  * inv_s,
                           (B_[r] + A_[r]) * inv_s,
                           A_[r]  * inv_s,
                           1.0,
                           t0_r,
                           sv_[r] * inv_s,
                           dispatch.lambda_g, dispatch.lambda_k,
                           20, true, kill_shape, dispatch.guess, pd, omega);
      }
      if (log_cdf >= 0.0) { bad = true; break; }
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

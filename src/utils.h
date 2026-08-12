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
typedef void (*RaceRawFun)(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok,
                           double* out, double min_ll, void* ctx_);

// Batch log-survivor at scalar t for truncation normalisation
typedef void (*RaceLogSAtTFun)(double t, const double* const* cols,
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
  // Shift for the log-space fallback integrand: the integrand is evaluated
  // as exp(log integrand - log_scale) so GSL still sees a natural-scale
  // function while the winner-pdf x loser-survivor product is formed in log
  // space.  Only used by gsl_f_race_scalar_logshift.
  double log_scale = 0.0;
};

double gsl_f_race_scalar(double t, void* p);
double gsl_f_race_scalar_logshift(double t, void* p);
double log_race_integrand_scalar(double t, const gsl_race_params_scalar* P);

// --------------------------------------------------------------------------
// Context object for race models, holding metadata and switches.
// --------------------------------------------------------------------------
struct ContextForRaceModels {
    double min_lik_for_pdf = 1e-10;
    bool use_posdrift = true;
    bool gng = false;
    int t0_index = -1;
    // Column indices for the Erlang timer mean parameters (mG, mK) in the raw
    // particle parameter matrix that calc_ll_oo passes to the C++ likelihood.
    // These store the user-visible means, NOT rates.  Call erlang_lambda_from_mean()
    // before passing values from these columns to erlang_log_surv/erlang_log_pdf.
    int mean_g_index = -1;   // column of mG (guess-clock mean)
    int mean_k_index = -1;   // column of mK (kill-clock mean)
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

    // LBA is represented by the shared BAwL kernels with k=0 and both timer
    // means fixed at their natural-scale off value (zero).  In this mode the
    // first five columns retain the LBA layout and optional BAwL columns must
    // not be read.
    bool bawl_k_fixed_zero = false;
    bool bawl_clocks_fixed_off = false;

    // BAwD launch-strength distribution: 0 = truncated normal (v, sv),
    // 1 = lognormal (mu, sigma).  Set from the c_name suffix "_LOGN" in
    // resolve_race_model_adapter(); the two layouts share column positions, so
    // this flag is the only thing that distinguishes them inside the kernels.
    // The values must match BAWD_LAUNCH_* in model_BAwD.h and the `launch`
    // argument of the exported dbawd/pbawd.
    int bawd_launch = BAWD_LAUNCH_LOGNORMAL;

    // Correlated *drift draws* through one shared standard-normal factor
    // (drift_factor.h).  The low-level row rho determines the accumulator's
    // signed factor variance share; the model's R Ttransform maps a
    // cell-level correlation into these row values before the kernel is
    // called, and rho == 0 leaves that racer independent of the shared draw.
    // Any model whose drift is N(v, sv^2) can use this route: the column
    // positions below tell the shared driver where v and sv live, and
    // corr_drift_generic_only suppresses the BAwL-specific exact and fused
    // fast routes for models that have no affine-in-drift survivor.
    bool corr_drift_active = false;
    int corr_drift_rho_index = -1;
    int corr_drift_v_col = -1;
    int corr_drift_sv_col = -1;
    bool corr_drift_generic_only = false;

    // Correlated RDMSWTN couples one directly specified pair of complete
    // finishing-time marginals with a Gaussian copula. This metadata is
    // intentionally distinct from the shared drift-factor path above; the two
    // are the "times" and "drifts" settings of RDMSWTNcorr(correlate =).
    bool rdmswtn_correlated = false;
    int rdmswtn_rho_index = -1;

    // Tri-state caches for optional accumulator levels:
    // -2 = unresolved (detect from data once), -1 = absent, >0 = factor code.
    int time_code = -2;
    int nogo_code = -2;

    // PDE-backed race models amortise one Fokker-Planck march over every
    // row that shares a parameter tuple, and over the scalar CDF calls the
    // censoring path makes at the truncation/censoring bounds.  Held by shared
    // pointer because the adapter is copied around, and allocated only by the
    // models that need it -- every analytic model leaves this null and pays
    // nothing.  Cleared once per particle; see fperace::SolveCache.
    std::shared_ptr<fperace::SolveCache> fpe_cache;
    std::shared_ptr<rlf::SolveCache> rlf_cache;

    bool has_global_kill() const {
      return is_global_kill && kill_active;
    }
};

// ---------------------------------------------------------------------------
// Two-boundary (DDM-shaped) models.
//
// c_log_likelihood_DDM_pt owns all of the truncation, censoring and go/no-go
// bookkeeping, and needs exactly two things from the model: the log DEFECTIVE
// density and the log DEFECTIVE cdf of a response, vectorised over trials.
// Everything else it does -- the interval masses, the log-sum-exp combinations,
// the all-finite fast path -- is model independent.  Routing those two through
// function pointers is what lets a second two-boundary model (the bounded OU of
// Smith & Ratcliff 2004, which is the DDM with leak) reuse the whole kernel
// rather than restate it, and it leaves the Wiener DDM on exactly the code it
// had before: its entries are thin shims that the compiler inlines away.
//
// This mirrors ContextForRaceModels above, which plays the same role for the
// race evaluator.
// ---------------------------------------------------------------------------
struct ContextForDDMModels {
  // PDE-backed two-boundary models amortise one Fokker-Planck march over every
  // trial sharing a parameter tuple.  Null for the analytic Wiener DDM, which
  // allocates nothing.  Cleared once per particle; see fpebou::SolveCache.
  std::shared_ptr<fpebou::SolveCache> bou_cache;
  int bnd_kind = 0;                 // FPE_BoundaryKind for a collapsing bound
  bool floor_raw_log_lik = false;
};

// (rts, Rs, cols, n_rows, mask, is_ok, out, floor, ctx)
//   Rs    1 = lower, 2 = upper -- the DDM's response coding
//   mask  rows with 0 are skipped and `out` is left untouched there
//   out   LOG defective density (d_raw) or LOG defective cdf (p_raw); the two
//         responses' values sum to 1 in the natural scale, not each alone
using DDMRawFun = void (*)(const double*, const int*, const double* const*, int,
                           const int*, const int*, double*, double,
                           ContextForDDMModels*);

struct DDMAdapter {
  DDMRawFun d_raw = nullptr;
  DDMRawFun p_raw = nullptr;
  emc2col::ColSpec col_spec{};
  // The analytic Wiener DDM can memoize numerical endpoint CDFs by exact
  // parameter key.  PDE-backed two-boundary adapters have their own solve
  // caches and a different parameter contract, so they leave this disabled.
  bool endpoint_cdf_cache = false;
  ContextForDDMModels ctx;
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
#include "model_GOM.h"
#include "model_RLF_kernels.h"
// Likewise: the BOU primitives need ContextForDDMModels defined above.
#include "model_BOU.h"

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
// and a geometric threshold excess.
static constexpr double PC_EPS = 1e-8;

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
  if (!R_FINITE(k)) return 0;
  const double z = std::floor(k + 0.5);
  if (z < 1.0) return 1;
  if (z > static_cast<double>(std::numeric_limits<int>::max() - 2))
    return std::numeric_limits<int>::max() - 2;
  return static_cast<int>(z);
}

inline std::vector<std::vector<double>> pcounter_stirling_logs(int nmax) {
  std::vector<std::vector<double>> st(static_cast<size_t>(nmax + 1),
                                      std::vector<double>(static_cast<size_t>(nmax + 1), R_NegInf));
  st[0][0] = 0.0;
  for (int n = 1; n <= nmax; ++n) {
    for (int j = 1; j <= n; ++j) {
      const double a = st[n - 1][j - 1];
      const double b = (j < n && n > 1) ? std::log(static_cast<double>(n - 1)) + st[n - 1][j]
                                         : R_NegInf;
      st[n][j] = (a == R_NegInf) ? b :
                 (b == R_NegInf) ? a :
                 (a > b ? a + std::log1p(std::exp(b - a))
                          : b + std::log1p(std::exp(a - b)));
    }
  }
  return st;
}

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
                             const std::vector<std::vector<double>>& st) {
  double logL, logM;
  pcounter_log_lm(t, nu, sv, logL, logM);
  if (sv <= PC_EPS)
    return logL + pcounter_log_rising(nu / gamma, n);
  const double shape = nu * nu / (sv * sv);
  const double rate = nu / (sv * sv);
  std::vector<double> terms;
  terms.reserve(static_cast<size_t>(n + 1));
  for (int j = 0; j <= n; ++j) {
    const double c = st[n][j];
    if (!R_FINITE(c)) { terms.push_back(R_NegInf); continue; }
    terms.push_back(c + pcounter_log_rising(shape, j) - j * std::log(gamma) -
                       j * std::log(rate + t));
  }
  return logL + pcounter_logsumexp(terms);
}

inline void pcounter_log_pi_phi(int n, double t, double nu, double sv, double gamma,
                                bool gamma_zero,
                                const std::vector<std::vector<double>>& st,
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
  logpi = n * logq - R::lgammafn(n + 1.0) + pcounter_log_h(n, t, nu, sv, gamma, st);
  logphi = std::log(gamma) + n * logq - R::lgammafn(n + 1.0) +
           pcounter_log_h(n + 1, t, nu, sv, gamma, st);
}

inline double pcounter_log_tail(int start, double t, double nu, double sv, double gamma,
                                double logr, bool gamma_zero,
                                const std::vector<std::vector<double>>& st, bool phi) {
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, st, lp, lf);
    terms.push_back((phi ? lf : lp) + n * logr);
  }
  return pcounter_logsumexp(terms);
}

inline double pcounter_log_fixed_cdf(int start, double t, double nu, double sv,
                                     double gamma, bool gamma_zero,
                                     const std::vector<std::vector<double>>& st) {
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, st, lp, lf);
    terms.push_back(lp);
  }
  return pcounter_logsumexp(terms);
}

inline double pcounter_log_geom_cdf(int start, double t, double nu, double sv,
                                    double gamma, double omega, bool gamma_zero,
                                    const std::vector<std::vector<double>>& st) {
  const double logr = std::log(omega) - std::log1p(omega);
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    if (n == start) { terms.push_back(R_NegInf); continue; }
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, st, lp, lf);
    terms.push_back(lp + std::log(-std::expm1((n - start) * logr)));
  }
  return pcounter_logsumexp(terms);
}

inline void pcounter_log_eval(double t, double nu, double sv, double gamma,
                              double k, double omega,
                              const std::vector<std::vector<double>>& st,
                              double& logf, double& logS, double& logF) {
  logf = R_NegInf; logS = 0.0; logF = R_NegInf;
  if (ISNAN(t) || !R_FINITE(nu) || !R_FINITE(sv) || !R_FINITE(gamma) ||
      !R_FINITE(k) || !R_FINITE(omega) || nu <= 0.0 || sv < 0.0 || gamma < 0.0 ||
      k <= 0.0 || omega < 0.0) return;
  if (R_PosInf == t) { logS = R_NegInf; logF = 0.0; return; }
  if (!R_FINITE(t) || t <= 0.0) return;

  const int kk = pcounter_k_int(k);
  const bool gamma_zero = gamma < PC_EPS;
  const bool sv_zero = sv < PC_EPS;
  const bool omega_zero = omega < PC_EPS;
  if (omega_zero) {
    std::vector<double> probs;
    probs.reserve(static_cast<size_t>(kk));
    for (int n = 0; n < kk; ++n) {
      double lp, lf;
      pcounter_log_pi_phi(n, t, nu, sv_zero ? 0.0 : sv,
                          gamma_zero ? 0.0 : gamma, gamma_zero, st, lp, lf);
      probs.push_back(lp);
    }
    logS = std::min(0.0, pcounter_logsumexp(probs));
    double lp, lf;
    pcounter_log_pi_phi(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                        gamma_zero ? 0.0 : gamma, gamma_zero, st, lp, lf);
    logf = lf;
    logF = logS > -1e-7 ? pcounter_log_fixed_cdf(kk, t, nu,
      sv_zero ? 0.0 : sv, gamma_zero ? 0.0 : gamma, gamma_zero, st)
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
                          gamma_zero ? 0.0 : gamma, gamma_zero, st, lp, lf);
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
                                gamma_zero ? 0.0 : gamma, logr, gamma_zero, st, false);
  logS = std::min(0.0, pcounter_logsumexp(
      std::vector<double>{loglow, (1.0 - kk) * logr + logtail}));
  const double logA = logMar - (gamma_zero ? 0.0 : logd);
  const double loglowf = pcounter_logsumexp(low_phi_r);
  double logtail_flux = pcounter_logdiffexp(logA, loglowf);
  logf = logp + (1.0 - kk) * logr + logtail_flux;
  if (nlow >= 0 && (loglowf > logA - 1e-7 || !R_FINITE(logf))) {
    logtail_flux = pcounter_log_tail(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                     gamma_zero ? 0.0 : gamma, logr, gamma_zero, st, true);
    logf = logp + (1.0 - kk) * logr + logtail_flux;
  }
  if (kk == 1) logF = std::log(-std::expm1(logLar));
  else if (logS > -1e-7)
    logF = pcounter_log_geom_cdf(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                 gamma_zero ? 0.0 : gamma, omega, gamma_zero, st);
  else logF = logS < 0.0 ? std::log(-std::expm1(logS)) : R_NegInf;
}

inline double dpcounter_scalar(double t, const double* par, void* /*ctx_*/) {
  for (int j = 0; j < emc2col::pcounter::N_REQ; ++j) if (R_IsNA(par[j])) return 0.0;
  if (par[emc2col::pcounter::nu] <= 0.0 || par[emc2col::pcounter::sv] < 0.0 ||
      par[emc2col::pcounter::gamma] < 0.0 || par[emc2col::pcounter::k] <= 0.0 ||
      par[emc2col::pcounter::omega] < 0.0) return 0.0;
  const int kk = pcounter_k_int(par[emc2col::pcounter::k]);
  const auto st = pcounter_stirling_logs(std::max(65, kk + 65));
  double lf, ls, lF;
  pcounter_log_eval(t - par[emc2col::pcounter::t0], par[0], par[1], par[2], par[3], par[4], st, lf, ls, lF);
  return R_FINITE(lf) ? std::exp(lf) : 0.0;
}

inline double ppcounter_scalar(double t, const double* par, void* /*ctx_*/) {
  for (int j = 0; j < emc2col::pcounter::N_REQ; ++j) if (R_IsNA(par[j])) return 0.0;
  if (par[0] <= 0.0 || par[1] < 0.0 || par[2] < 0.0 || par[3] <= 0.0 || par[4] < 0.0) return 0.0;
  const auto st = pcounter_stirling_logs(std::max(65, pcounter_k_int(par[3]) + 65));
  double lf, ls, lF;
  pcounter_log_eval(t - par[5], par[0], par[1], par[2], par[3], par[4], st, lf, ls, lF);
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
  int maxk = 1;
  for (int i = 0; i < n_rows; ++i) if (mask[i] && isok[i] && R_FINITE(kk[i])) maxk = std::max(maxk, pcounter_k_int(kk[i]));
  const auto st = pcounter_stirling_logs(std::max(65, maxk + 65));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(nu[i]) || R_IsNA(sv[i]) || R_IsNA(ga[i]) ||
        R_IsNA(kk[i]) || R_IsNA(om[i]) || R_IsNA(t0[i]) || nu[i] <= 0.0 ||
        sv[i] < 0.0 || ga[i] < 0.0 || kk[i] <= 0.0 || om[i] < 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    double lf, ls, lF;
    pcounter_log_eval(rt[i] - t0[i], nu[i], sv[i], ga[i], kk[i], om[i], st, lf, ls, lF);
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
  int maxk = 1;
  for (int i = 0; i < n_rows; ++i) if (mask[i] && isok[i] && R_FINITE(kk[i])) maxk = std::max(maxk, pcounter_k_int(kk[i]));
  const auto st = pcounter_stirling_logs(std::max(65, maxk + 65));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(nu[i]) || R_IsNA(sv[i]) || R_IsNA(ga[i]) ||
        R_IsNA(kk[i]) || R_IsNA(om[i]) || R_IsNA(t0[i]) || nu[i] <= 0.0 ||
        sv[i] < 0.0 || ga[i] < 0.0 || kk[i] <= 0.0 || om[i] < 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    double lf, ls, lF;
    pcounter_log_eval(rt[i] - t0[i], nu[i], sv[i], ga[i], kk[i], om[i], st, lf, ls, lF);
    out[i] = raw_log_value(ls, min_ll, floor_raw);
  }
}

inline void pcounter_logS_at_t(double t, const double* const* cols,
                               int /*n_rows_total*/, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  (void)ctx_;
  const double* kk = cols[emc2col::pcounter::k];
  int maxk = 1;
  for (int i = 0; i < n_unique_trials * n_lR; ++i) if (isok_all[i] && R_FINITE(kk[i])) maxk = std::max(maxk, pcounter_k_int(kk[i]));
  const auto st = pcounter_stirling_logs(std::max(65, maxk + 65));
  const double* nu = cols[emc2col::pcounter::nu];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double sum = 0.0;
    for (int k = 0; k < n_lR; ++k) {
      const int r = j * n_lR + k;
      if (!isok_all[r] || R_IsNA(nu[r]) || R_IsNA(sv[r]) || R_IsNA(ga[r]) ||
          R_IsNA(kk[r]) || R_IsNA(om[r]) || R_IsNA(t0[r]) || nu[r] <= 0.0 ||
          sv[r] < 0.0 || ga[r] < 0.0 || kk[r] <= 0.0 || om[r] < 0.0) {
        sum = R_NegInf; break;
      }
      double lf, ls, lF;
      pcounter_log_eval(t - t0[r], nu[r], sv[r], ga[r], kk[r], om[r], st, lf, ls, lF);
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(mu_[i]) || R_IsNA(sigma_[i]) || R_IsNA(tau_[i]) ||
        !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    // dexg derives the density in log terms internally; take that value
    // directly instead of exponentiating and re-logging it.
    out[i] = raw_log_value(dexg(rt[i], mu_[i], sigma_[i], tau_[i], true),
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(mu_[i]) || R_IsNA(sigma_[i]) || R_IsNA(tau_[i]) ||
        !isok[i] || sigma_[i] <= 0.0 || tau_[i] <= 0.0) {
      out[i] = 0.0;
      continue;
    }
    // Direct upper-tail log probability (protected inside pexg): stays
    // finite after the natural lower-tail CDF saturates to one.
    const double log_surv = pexg(rt[i], mu_[i], sigma_[i], tau_[i], false, true);
    if (!R_FINITE(log_surv)) {
      out[i] = (log_surv == R_NegInf) ? raw_log_zero(min_ll, floor_raw) : 0.0;
      continue;
    }
    out[i] = std::fmin(log_surv, 0.0);
  }
}

inline void rexg_logS_at_t(double t, const double* const* cols,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* mu_    = cols[emc2col::rexg::mu];
  const double* sigma_ = cols[emc2col::rexg::sigma];
  const double* tau_   = cols[emc2col::rexg::tau];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR && !bad; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(mu_[r]) || R_IsNA(sigma_[r]) || R_IsNA(tau_[r]) ||
          sigma_[r] <= 0.0 || tau_[r] <= 0.0) { bad = true; break; }
      const double log_surv = pexg(t, mu_[r], sigma_[r], tau_[r], false, true);
      if (!R_FINITE(log_surv)) {
        if (log_surv == R_NegInf) { bad = true; break; }
        continue;  // NA parameters already screened; treat as no contribution
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
    ctx->use_posdrift, false, ks, local_guess, omega
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
    ctx->use_posdrift, false, ks, local_guess, omega
  );
}

inline void dbawl_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
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
      if (ba_natural_pdf(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval,
                         pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, pdf)) {
        out[i] = (pdf > 0.0)
          ? raw_log_value(std::log(pdf), min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      } else {
        const double log_pdf = log_ba_pdf(tt, A_[i], B_[i] + A_[i], v_[i],
                                          sv_[i], kval, pd, BAWL_DENOM_FLOOR);
        out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
          ? raw_log_value(log_pdf, min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      }
      continue;
    }
    // Pass raw rt and t0; core function uses t0 to split EAM vs erlang time.
    const double log_pdf = dkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, kval, lg, lk,
      pd, true, ks, local_guess, omega
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
        const double cdf_inf = pd ? 1.0 : pnorm_std(v_[i] / sv_[i], true, false);
        if (cdf_inf >= 1.0) out[i] = raw_log_zero(min_ll, floor_raw);
        else out[i] = R_FINITE(cdf_inf) ? std::log1p(-cdf_inf) : 0.0;
        continue;
      }
      const double b_i = B_[i] + A_[i];
      double cdf = 0.0;
      if (ba_natural_cdf(tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd,
                         LBA_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
        // Acceptance guarantees cdf in [0, 1 - 1e-8).
        out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
      } else {
        const double log_cdf = log_ba_cdf(tt, A_[i], b_i, v_[i], sv_[i],
                                          0.0, pd, LBA_DENOM_FLOOR);
        if (log_cdf >= 0.0) {
          out[i] = raw_log_zero(min_ll, floor_raw);
        } else {
          out[i] = R_FINITE(log_cdf) ? log1m_exp(log_cdf) : 0.0;
        }
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
      if (ba_natural_cdf(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval,
                         pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
        out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
      } else {
        const double log_cdf = log_ba_cdf(tt, A_[i], B_[i] + A_[i], v_[i],
                                          sv_[i], kval, pd, BAWL_DENOM_FLOOR);
        if (log_cdf >= 0.0) out[i] = raw_log_zero(min_ll, floor_raw);
        else out[i] = R_FINITE(log_cdf) ? log1m_exp(log_cdf) : 0.0;
      }
      continue;
    }
    const double log_cdf = pkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, kval, lg, lk,
      pd, true, ks, local_guess, omega
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
          const double cdf_inf = pd ? 1.0 : pnorm_std(v_[r] / sv_[r], true, false);
          if (cdf_inf >= 1.0) { bad = true; break; }
          logS += std::log1p(-cdf_inf);
          continue;
        }
        double cdf = 0.0;
        if (ba_natural_cdf(tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], 0.0, pd,
                           LBA_DENOM_FLOOR, BA_ACCEPT_RAW, cdf)) {
          if (cdf > 0.0) logS += std::log1p(-cdf);
          continue;
        }
        const double log_cdf = log_ba_cdf(tt, A_[r], B_[r] + A_[r], v_[r],
                                          sv_[r], 0.0, pd, LBA_DENOM_FLOOR);
        if (log_cdf >= 0.0) { bad = true; break; }
        if (R_FINITE(log_cdf)) logS += log1m_exp(log_cdf);
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
          pd, true, ks, local_guess, omega
        );
        if (!R_FINITE(log_cdf)) continue;
        if (log_cdf >= 0.0) { bad = true; break; }
        logS += log1m_exp(log_cdf);
        continue;
      }
      // Both EAM and erlang contribute; pass raw t and t0_r.
      const double log_cdf = pkilledleakyba_norm(
        t, v_[r], B_[r] + A_[r], A_[r], sv_[r], t0_r, kval, lg, lk,
        pd, true, ks, local_guess, omega
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

inline double dbawd_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawd::v])) return 0.0;
  const double tt = t - par[emc2col::bawd::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  return bawd_pdf_scalar_natural(
    tt, par[emc2col::bawd::A],
    par[emc2col::bawd::B] + par[emc2col::bawd::A],
    par[emc2col::bawd::v], par[emc2col::bawd::sv],
    par[emc2col::bawd::k], par[emc2col::bawd::ell],
    bawd_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline double pbawd_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[emc2col::bawd::v])) return 0.0;
  const double tt = t - par[emc2col::bawd::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf is handled inside the kernel and returns F_max, not 1.
  return bawd_cdf_scalar_natural(
    tt, par[emc2col::bawd::A],
    par[emc2col::bawd::B] + par[emc2col::bawd::A],
    par[emc2col::bawd::v], par[emc2col::bawd::sv],
    par[emc2col::bawd::k], par[emc2col::bawd::ell],
    bawd_launch_of(ctx), ctx ? ctx->use_posdrift : true);
}

inline void dbawd_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawd_launch_of(ctx);
  const double* p1_ = cols[emc2col::bawd::v];
  const double* p2_ = cols[emc2col::bawd::sv];
  const double* B_  = cols[emc2col::bawd::B];
  const double* A_  = cols[emc2col::bawd::A];
  const double* t0_ = cols[emc2col::bawd::t0];
  const double* k_  = cols[emc2col::bawd::k];
  const double* ell_ = cols[emc2col::bawd::ell];
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
    const double log_pdf = bawd_log_pdf(tt, A_[i], B_[i] + A_[i], p1_[i],
                                        p2_[i], k_[i], ell_[i], launch, pd);
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
  const double* p1_ = cols[emc2col::bawd::v];
  const double* p2_ = cols[emc2col::bawd::sv];
  const double* B_  = cols[emc2col::bawd::B];
  const double* A_  = cols[emc2col::bawd::A];
  const double* t0_ = cols[emc2col::bawd::t0];
  const double* k_  = cols[emc2col::bawd::k];
  const double* ell_ = cols[emc2col::bawd::ell];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const double log_cdf = bawd_log_cdf(tt, A_[i], B_[i] + A_[i], p1_[i],
                                        p2_[i], k_[i], ell_[i], launch, pd);
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

inline void bawd_logS_at_t(double t, const double* const* cols,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawd_launch_of(ctx);
  const double* p1_ = cols[emc2col::bawd::v];
  const double* p2_ = cols[emc2col::bawd::sv];
  const double* B_  = cols[emc2col::bawd::B];
  const double* A_  = cols[emc2col::bawd::A];
  const double* t0_ = cols[emc2col::bawd::t0];
  const double* k_  = cols[emc2col::bawd::k];
  const double* ell_ = cols[emc2col::bawd::ell];
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
      const double log_cdf = bawd_log_cdf(tt, A_[r], B_[r] + A_[r], p1_[r],
                                          p2_[r], k_[r], ell_[r], launch, pd);
      if (log_cdf >= 0.0) { bad = true; break; }
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

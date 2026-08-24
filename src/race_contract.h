#ifndef EMC2_RACE_CONTRACT_H
#define EMC2_RACE_CONTRACT_H

#include <Rcpp.h>
#include <memory>

#include "col_registry.h"

namespace fperace {
struct SolveCache;
}
namespace fpebou {
struct SolveCache;
}
namespace rlf {
struct SolveCache;
}

// Launch selectors mirror the model-specific BAW* constants without making
// this contract header depend on model implementation headers.
namespace emc2race_contract {
constexpr int kBawdLaunchLognormal = 1;  // BAWD_LAUNCH_LOGNORMAL
constexpr int kBtawlLaunchNormal = 0;   // BTAWL_LAUNCH_NORMAL
constexpr int kBawlLaunchNormal = 0;    // BAWL_LAUNCH_NORMAL
}


// Scalar function pointer types for race model PDF/CDF adapters used by GSL
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

    // Timed race mixtures need unfloored components before the final log_sum_exp.
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
    int bawd_launch = emc2race_contract::kBawdLaunchLognormal;

    // BTAwL launch-strength distribution.  This is kept separate from the
    // BAwD field because both models can be selected by one process and their
    // adapters may coexist in diagnostics/tests.
    int btawl_launch = emc2race_contract::kBtawlLaunchNormal;
    // BTAwL clearance chart: false samples the intrinsic time constant tau;
    // true samples the observable transient endpoint Ttrans and back-solves
    // tau before entering the shared geometry.
    bool btawl_ttrans_chart = false;
    static constexpr int btawl_tau_cache_size = 4;
    mutable double btawl_tau_cache_k[btawl_tau_cache_size] =
      {R_NaN, R_NaN, R_NaN, R_NaN};
    mutable double btawl_tau_cache_clear[btawl_tau_cache_size] =
      {R_NaN, R_NaN, R_NaN, R_NaN};
    mutable double btawl_tau_cache_value[btawl_tau_cache_size] =
      {R_NaN, R_NaN, R_NaN, R_NaN};
    mutable int btawl_tau_cache_next = 0;

    // Fixed BAwD clearance-fade exponent parsed from the c_name suffix.
    // Allowed values are mirrored in R/model_BAwD.R.
    double bawd_gamma = 0.0;
    double bawd_rho = R_PosInf;

    // BAwL launch-strength distribution, on the same convention: 0 = normal
    // (v, sv), 1 = lognormal (mu, sigma).  Set from the "_LOGN" c_name suffix.
    // The default is Gaussian BAwL; LBA is always normal.
    int bawl_launch = emc2race_contract::kBawlLaunchNormal;

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

// Shared raw-log helpers used by the batch raw kernels across race models.
// Kept here (next to ContextForRaceModels) so both the analytic kernels in
// utils.h and the LNR kernels in model_lnr.cpp can share one definition.
inline bool raw_floor_log_lik(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx == nullptr || ctx->floor_raw_log_lik;
}

inline double raw_log_zero(double min_ll, bool floor_raw) {
  return floor_raw ? min_ll : R_NegInf;
}

inline double raw_log_value(double log_x, double min_ll, bool floor_raw) {
  if (!R_FINITE(log_x)) return raw_log_zero(min_ll, floor_raw);
  return floor_raw ? ((log_x > min_ll) ? log_x : min_ll) : log_x;
}

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

#endif // EMC2_RACE_CONTRACT_H

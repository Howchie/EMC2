#ifndef PARTICLE_LL_UTILS_H
#define PARTICLE_LL_UTILS_H
#include <RcppArmadillo.h>
#include "utility_functions.h" // For any existing utilities we might use
#include "lnr_functions.h"
#include "wald_functions.h"
#include <gsl/gsl_integration.h>

// Scalar (single-RT / single-accumulator) model helpers for fast GSL integration.
typedef double (*RacePdf1Fun)(double rt, const double* par, void* model_specific_context);
typedef double (*RaceCdf1Fun)(double rt, const double* par, void* model_specific_context);

// Raw-buffer variants for the all-finite-no-truncation fast path.
// Iterates all n_rows; where mask[i]==1 writes log-density (dfun) or log(1-cdf) (pfun)
// into out[i].  Where mask[i]==0 out[i] is left unchanged.
// pars_cm is column-major: column j starts at pars_cm + j * n_rows.
typedef void (*RaceRawFun)(const double* rt,
                           const double* pars_cm,
                           int n_rows,
                           const int* mask,
                           const int* isok,
                           double* out,
                           double min_ll,
                           void* ctx);

// Batch log-survivor at a scalar time point t.
// For each unique trial j where trunc_mask[j]==1:
//   logS_out[j] = sum_{k=0}^{n_lR-1} log(1 - F_k(t, pars[j*n_lR + k]))
// pars_cm is column-major: column c starts at pars_cm + c * n_rows_total.
typedef void (*RaceLogSAtTFun)(
    double t,
    const double* pars_cm,
    int n_rows_total,
    int n_lR,
    int n_par,
    const int* trunc_mask,
    int n_unique_trials,
    const int* isok_all,
    void* ctx,
    double* logS_out
);

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

struct RaceSharedState;  // defined in particle_ll.cpp

double c_log_likelihood_race(Rcpp::NumericMatrix pars, Rcpp::DataFrame dadm,
                                        RacePdf1Fun pdf1,      // Scalar PDF for integration
                                        RaceCdf1Fun cdf1,      // Scalar CDF for integration
                                        const int n_trials,
                                        const Rcpp::LogicalVector winner,
                                        const Rcpp::IntegerVector expand,
                                        double min_ll,
                                        const Rcpp::LogicalVector ok_params, // Parameter validity for each row in 'pars' matrix
                                        int n_lR, // Number of accumulators in the race
                                        void *model_context_for_funcs,
                                        bool all_finite_trials,
                                        RaceRawFun model_dfun_raw = nullptr,
                                        RaceRawFun model_pfun_raw = nullptr,
                                        RaceLogSAtTFun logS_at_t = nullptr,
                                        RaceSharedState* shared = nullptr,
                                        Rcpp::NumericVector* trial_ll_out = nullptr
);

// Exported wrappers
Rcpp::NumericVector calc_ll(Rcpp::NumericMatrix p_matrix, Rcpp::DataFrame data,
                            Rcpp::NumericVector constants, Rcpp::List designs,
                            Rcpp::String type, Rcpp::List bounds,
                            Rcpp::List transforms, Rcpp::List pretransforms,
                            Rcpp::CharacterVector p_types, double min_ll,
                            Rcpp::List trend);
							
// Structures
struct ContextForRaceModels {
    double min_lik_for_pdf;
    bool use_posdrift=true;
    // Whether the race can place non-zero probability mass at +Inf.
    // This controls omission/right-tail bookkeeping in likelihood code.
    bool defective_upper_tail=false;
	bool log_out=false;
	bool gng=false;
	std::string LogicalRule;
    // Optional per-particle fast-kernel hint:
    // 0 = auto detect in kernel; 1 = zero-variability branch; 2 = nonzero branch.
    int mode_hint = 0;
};

// Scalar adapters (single-RT, single-parameter-row) used by GSL integration.
inline double dlba_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const double v  = par[0];
  const double sv = par[1];
  const double b  = par[2] + par[3];
  const double A  = par[3];
  const double t0 = par[4];

  if (R_IsNA(v)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return dlba_norm(tt, A, b, v, sv, ctx->use_posdrift);
}

inline double plba_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const double v  = par[0];
  const double sv = par[1];
  const double b  = par[2] + par[3];
  const double A  = par[3];
  const double t0 = par[4];

  if (R_IsNA(v)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return plba_norm(tt, A, b, v, sv, ctx->use_posdrift);
}

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

// ---- Raw-buffer implementations (log-density and log-survivor) ----
// These avoid per-particle Rcpp vector allocation in the all-finite-no-truncation fast path.
// Each function writes into out[i] for every i where mask[i]==1; other slots are untouched.

// LBA: column layout v=0, sv=1, B=2, A=3, t0=4
inline void dlba_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx->use_posdrift;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_ = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  #pragma omp simd
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double pdf = dlba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], pd);
    out[i] = (pdf > 0.0 && emc2_isfinite(pdf)) ? std::log(pdf) : min_ll;
  }
}

inline void plba_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx->use_posdrift;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_ = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  #pragma omp simd
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double cdf = plba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], pd);
    // log(1 - cdf); plba_norm guarantees cdf in [0,1]
    if (cdf >= 1.0) { out[i] = min_ll; continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

// RDM: column layout v=0, B=1, A=2, t0=3, s=4
inline void drdm_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* /*ctx_*/) {
  const double* v_  = pars_cm + 0 * n_rows;
  const double* B_  = pars_cm + 1 * n_rows;
  const double* A_  = pars_cm + 2 * n_rows;
  const double* t0_ = pars_cm + 3 * n_rows;
  const double* s_  = pars_cm + 4 * n_rows;
  #pragma omp simd
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double inv_s = 1.0 / s_[i];  // one division, three multiplications below
    const double pdf = digt_impl(tt, B_[i] * inv_s + 0.5 * A_[i] * inv_s,
                            v_[i] * inv_s, 0.5 * A_[i] * inv_s);
    out[i] = (pdf > 0.0 && R_FINITE(pdf)) ? std::log(pdf) : min_ll;
  }
}

inline void prdm_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* /*ctx_*/) {
  const double* v_  = pars_cm + 0 * n_rows;
  const double* B_  = pars_cm + 1 * n_rows;
  const double* A_  = pars_cm + 2 * n_rows;
  const double* t0_ = pars_cm + 3 * n_rows;
  const double* s_  = pars_cm + 4 * n_rows;
  #pragma omp simd
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double inv_s = 1.0 / s_[i];  // one division, three multiplications below
    const double cdf = pigt_impl(tt, B_[i] * inv_s + 0.5 * A_[i] * inv_s,
                            v_[i] * inv_s, 0.5 * A_[i] * inv_s);
    if (cdf >= 1.0) { out[i] = min_ll; continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

// LNR: column layout m=0, s=1, t0=2
inline void dlnr_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* /*ctx_*/) {
  const double* m_  = pars_cm + 0 * n_rows;
  const double* s_  = pars_cm + 1 * n_rows;
  const double* t0_ = pars_cm + 2 * n_rows;
  #pragma omp simd
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(m_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double log_pdf = dlnorm_std(tt, m_[i], s_[i], true);
    out[i] = (R_FINITE(log_pdf) && log_pdf > min_ll) ? log_pdf : min_ll;
  }
}

inline void plnr_raw(const double* rt, const double* pars_cm, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* /*ctx_*/) {
  const double* m_  = pars_cm + 0 * n_rows;
  const double* s_  = pars_cm + 1 * n_rows;
  const double* t0_ = pars_cm + 2 * n_rows;
  #pragma omp simd
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(m_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double logS = lnorm_log_surv_std(tt, m_[i], s_[i]);
    out[i] = R_FINITE(logS) ? logS : min_ll;
  }
}

// LBA: column layout v=0, sv=1, B=2, A=3, t0=4
inline void lba_logS_at_t(double t, const double* pars_cm,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx->use_posdrift;
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
      if (tt <= 0.0) continue;            // CDF(t)=0 → log(1-0)=0, no contribution
      const double cdf = plba_norm(tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], pd);
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// RDM: column layout v=0, B=1, A=2, t0=3, s=4
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
      const double cdf = pigt_impl(tt, B_[r] * inv_s + 0.5 * A_[r] * inv_s,
                              v_[r] * inv_s, 0.5 * A_[r] * inv_s);
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// LNR: column layout m=0, s=1, t0=2
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

// Adapter prototypes
Rcpp::NumericVector lba_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector winner,
                                            Rcpp::LogicalVector is_ok,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    // Pass use_posdrift from context to dlba_c
    return dlba_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok, ctx->use_posdrift);
}

Rcpp::NumericVector lba_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector winner,
                                            Rcpp::LogicalVector is_ok,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    // Pass use_posdrift from context to plba_c
    return plba_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok, ctx->use_posdrift);
}

// Static adapter for RDM dfun
Rcpp::NumericVector rdm_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector winner,
                                            Rcpp::LogicalVector is_ok,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    return drdm_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for RDM pfun
Rcpp::NumericVector rdm_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector winner,
                                            Rcpp::LogicalVector is_ok,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    return prdm_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LNR dfun
Rcpp::NumericVector lnr_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector winner,
                                            Rcpp::LogicalVector is_ok,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context); // Reusing same context struct
    return dlnr_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LNR pfun
Rcpp::NumericVector lnr_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector winner,
                                            Rcpp::LogicalVector is_ok,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context); // Reusing same context struct
    return plnr_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// ============================================================
// BAwL (Ballistic Accumulator with Leak) adapters
// Column layout: v=0, sv=1, B=2, A=3, t0=4, k=5
// ============================================================

inline double dbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[4];
  if (tt <= 0.0) return 0.0;
  // Route through the leaky kernel in all cases; its internal k->0 branch
  // handles the exact LBA limit.
  return dleakyba_norm(tt, par[3], par[2] + par[3], par[0], par[1], par[5],
                       ctx->use_posdrift);
}

inline double pbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[4];
  if (tt <= 0.0) return 0.0;
  // Route through the leaky kernel in all cases; its internal k->0 branch
  // handles the exact LBA limit.
  return pleakyba_norm(tt, par[3], par[2] + par[3], par[0], par[1], par[5],
                       ctx->use_posdrift);
}

inline void dbawl_raw(const double* rt, const double* pars_cm, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd    = ctx->use_posdrift;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_ = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  const double* k_  = pars_cm + 5 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double pdf = dleakyba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], k_[i], pd);
    out[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : min_ll;
  }
}

inline void pbawl_raw(const double* rt, const double* pars_cm, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd    = ctx->use_posdrift;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* sv_ = pars_cm + 1 * n_rows;
  const double* B_  = pars_cm + 2 * n_rows;
  const double* A_  = pars_cm + 3 * n_rows;
  const double* t0_ = pars_cm + 4 * n_rows;
  const double* k_  = pars_cm + 5 * n_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double cdf = pleakyba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], k_[i], pd);
    if (cdf >= 1.0) { out[i] = min_ll; continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

// BAwL: column layout v=0, sv=1, B=2, A=3, t0=4, k=5
inline void bawl_logS_at_t(double t, const double* pars_cm,
                            int n_rows_total, int n_lR, int /*n_par*/,
                            const int* trunc_mask, int n_unique_trials,
                            const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd    = ctx->use_posdrift;
  const double* v_  = pars_cm + 0 * n_rows_total;
  const double* sv_ = pars_cm + 1 * n_rows_total;
  const double* B_  = pars_cm + 2 * n_rows_total;
  const double* A_  = pars_cm + 3 * n_rows_total;
  const double* t0_ = pars_cm + 4 * n_rows_total;
  const double* k_  = pars_cm + 5 * n_rows_total;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;
      const double cdf = pleakyba_norm(tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], k_[r], pd);
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

Rcpp::NumericVector bawl_dfun_adapter(Rcpp::NumericVector rt,
                                      Rcpp::NumericMatrix pars,
                                      Rcpp::LogicalVector winner,
                                      Rcpp::LogicalVector is_ok,
                                      void* context) {
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
  return dleakyba_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok, ctx->use_posdrift);
}

Rcpp::NumericVector bawl_pfun_adapter(Rcpp::NumericVector rt,
                                      Rcpp::NumericMatrix pars,
                                      Rcpp::LogicalVector winner,
                                      Rcpp::LogicalVector is_ok,
                                      void* context) {
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
  return pleakyba_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok, ctx->use_posdrift);
}

// ============================================================
// RDMSWTN adapters
// Column layout: v=0, B=1, A=2, t0=3, s=4, sv=5
// ============================================================

inline double drdmswtn_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[3];
  if (tt <= 0.0) return 0.0;
  const double inv_s = 1.0 / par[4];
  return drdmswtn(tt,
                  (par[1] + par[2]) * inv_s,  // b = (B+A)/s
                  par[0] * inv_s,              // v/s
                  par[2] * inv_s,              // A/s
                  par[5] * inv_s,              // sv/s
                  1.0, 0.0, 20, false);
}

inline double prdmswtn_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[3];
  if (tt <= 0.0) return 0.0;
  const double inv_s = 1.0 / par[4];
  return prdmswtn(tt,
                  (par[1] + par[2]) * inv_s,
                  par[0] * inv_s,
                  par[2] * inv_s,
                  par[5] * inv_s,
                  1.0, 0.0, 20, false);
}

inline void drdmswtn_raw(const double* rt, const double* pars_cm, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int mode_hint = ctx ? ctx->mode_hint : 0;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* B_  = pars_cm + 1 * n_rows;
  const double* A_  = pars_cm + 2 * n_rows;
  const double* t0_ = pars_cm + 3 * n_rows;
  const double* s_  = pars_cm + 4 * n_rows;
  const double* sv_ = pars_cm + 5 * n_rows;
  const double sv_eps = 1e-10;

  bool all_sv_zero = (mode_hint == 1);
  if (mode_hint == 0) {
    all_sv_zero = true;
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) continue;
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0) continue;
      const double sv_i = sv_[i];
      if (std::isfinite(sv_i) && std::fabs(sv_i) > sv_eps) {
        all_sv_zero = false;
        break;
      }
    }
  }

  if (all_sv_zero) {
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0) { out[i] = min_ll; continue; }
      const double inv_s = 1.0 / s_[i];
      const double pdf = dwald(tt,
                               (B_[i] + A_[i]) * inv_s,
                               v_[i] * inv_s,
                               1.0,
                               A_[i] * inv_s,
                               false);
      out[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : min_ll;
    }
    return;
  }

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    double pdf;
    if (!std::isfinite(sv_[i]) || std::fabs(sv_[i]) <= sv_eps) {
      const double inv_s = 1.0 / s_[i];
      pdf = dwald(tt,
                  (B_[i] + A_[i]) * inv_s,
                  v_[i] * inv_s,
                  1.0,
                  A_[i] * inv_s,
                  false);
    } else {
      const double inv_s = 1.0 / s_[i];
      pdf = drdmswtn(tt,
                     (B_[i] + A_[i]) * inv_s,
                     v_[i]  * inv_s,
                     A_[i]  * inv_s,
                     sv_[i] * inv_s,
                     1.0, 0.0, 20, false);
    }
    out[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : min_ll;
  }
}

inline void prdmswtn_raw(const double* rt, const double* pars_cm, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int mode_hint = ctx ? ctx->mode_hint : 0;
  const double* v_  = pars_cm + 0 * n_rows;
  const double* B_  = pars_cm + 1 * n_rows;
  const double* A_  = pars_cm + 2 * n_rows;
  const double* t0_ = pars_cm + 3 * n_rows;
  const double* s_  = pars_cm + 4 * n_rows;
  const double* sv_ = pars_cm + 5 * n_rows;
  const double sv_eps = 1e-10;

  bool all_sv_zero = (mode_hint == 1);
  if (mode_hint == 0) {
    all_sv_zero = true;
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) continue;
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0) continue;
      const double sv_i = sv_[i];
      if (std::isfinite(sv_i) && std::fabs(sv_i) > sv_eps) {
        all_sv_zero = false;
        break;
      }
    }
  }

  if (all_sv_zero) {
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0) { out[i] = 0.0; continue; }
      const double inv_s = 1.0 / s_[i];
      const double cdf = pwald(tt,
                               (B_[i] + A_[i]) * inv_s,
                               v_[i] * inv_s,
                               1.0,
                               A_[i] * inv_s,
                               false);
      if (cdf >= 1.0) { out[i] = min_ll; continue; }
      out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
    }
    return;
  }

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    double cdf;
    if (!std::isfinite(sv_[i]) || std::fabs(sv_[i]) <= sv_eps) {
      const double inv_s = 1.0 / s_[i];
      cdf = pwald(tt,
                  (B_[i] + A_[i]) * inv_s,
                  v_[i] * inv_s,
                  1.0,
                  A_[i] * inv_s,
                  false);
    } else {
      const double inv_s = 1.0 / s_[i];
      cdf = prdmswtn(tt,
                     (B_[i] + A_[i]) * inv_s,
                     v_[i]  * inv_s,
                     A_[i]  * inv_s,
                     sv_[i] * inv_s,
                     1.0, 0.0, 20, false);
    }
    if (cdf >= 1.0) { out[i] = min_ll; continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
  }
}

// RDMSWTN: column layout v=0, B=1, A=2, t0=3, s=4, sv=5
inline void rdmswtn_logS_at_t(double t, const double* pars_cm,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int mode_hint = ctx ? ctx->mode_hint : 0;
  const double* v_  = pars_cm + 0 * n_rows_total;
  const double* B_  = pars_cm + 1 * n_rows_total;
  const double* A_  = pars_cm + 2 * n_rows_total;
  const double* t0_ = pars_cm + 3 * n_rows_total;
  const double* s_  = pars_cm + 4 * n_rows_total;
  const double* sv_ = pars_cm + 5 * n_rows_total;
  const double sv_eps = 1e-10;
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
      double cdf;
      if (mode_hint == 1) {
        const double inv_s = 1.0 / s_[r];
        cdf = pwald(tt,
                    (B_[r] + A_[r]) * inv_s,
                    v_[r] * inv_s,
                    1.0,
                    A_[r] * inv_s,
                    false);
      } else if (mode_hint == 2) {
        const double inv_s = 1.0 / s_[r];
        cdf = prdmswtn(tt,
                       (B_[r] + A_[r]) * inv_s,
                       v_[r]  * inv_s,
                       A_[r]  * inv_s,
                       sv_[r] * inv_s,
                       1.0, 0.0, 20, false);
      } else if (!std::isfinite(sv_[r]) || std::fabs(sv_[r]) <= sv_eps) {
        const double inv_s = 1.0 / s_[r];
        cdf = pwald(tt,
                    (B_[r] + A_[r]) * inv_s,
                    v_[r] * inv_s,
                    1.0,
                    A_[r] * inv_s,
                    false);
      } else {
        const double inv_s = 1.0 / s_[r];
        cdf = prdmswtn(tt,
                       (B_[r] + A_[r]) * inv_s,
                       v_[r]  * inv_s,
                       A_[r]  * inv_s,
                       sv_[r] * inv_s,
                       1.0, 0.0, 20, false);
      }
      if (cdf >= 1.0) { bad = true; break; }
      if (cdf > 0.0) logS += std::log1p(-cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

Rcpp::NumericVector rdmswtn_dfun_adapter(Rcpp::NumericVector rt,
                                         Rcpp::NumericMatrix pars,
                                         Rcpp::LogicalVector winner,
                                         Rcpp::LogicalVector is_ok,
                                         void* context) {
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
  return drdmswtn_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

Rcpp::NumericVector rdmswtn_pfun_adapter(Rcpp::NumericVector rt,
                                         Rcpp::NumericMatrix pars,
                                         Rcpp::LogicalVector winner,
                                         Rcpp::LogicalVector is_ok,
                                         void* context) {
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
  return prdmswtn_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
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

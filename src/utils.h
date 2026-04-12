#ifndef PARTICLE_LL_UTILS_H
#define PARTICLE_LL_UTILS_H
#include <RcppArmadillo.h>
#include "utility_functions.h" // For any existing utilities we might use
#include <gsl/gsl_integration.h>

// Function pointer types for model's dfun and pfun
// Added LogicalVector is_ok for bound checking and void* context for extra args
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

// function pointers for integration helpers (truncation and censoring)
struct gsl_race_params {const Rcpp::NumericMatrix* p_trial;
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

double c_log_likelihood_race(Rcpp::NumericMatrix pars, Rcpp::DataFrame dadm,
                                        RacePdfFun model_dfun, // Pointer to the model's PDF adapter function
                                        RaceCdfFun model_pfun, // Pointer to the model's CDF adapter function
                                        RacePdf1Fun pdf1,      // Scalar PDF for integration
                                        RaceCdf1Fun cdf1,      // Scalar CDF for integration
                                        const int n_trials,
                                        const Rcpp::LogicalVector winner,
                                        const Rcpp::IntegerVector
                                        expand,   
                                        double min_ll, 
                                        const Rcpp::LogicalVector
                                        ok_params, // Parameter validity for each row in 'pars' matrix
                                        int n_lR, // Number of accumulators in the race
                                        void *model_context_for_funcs,
                                        bool all_finite_trials,
                                        RaceRawFun model_dfun_raw = nullptr,
                                        RaceRawFun model_pfun_raw = nullptr
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
	bool log_out=false;
	bool gng=false;
	std::string LogicalRule;
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
  const double v = par[0];
  const double B = par[1];
  const double A = par[2];
  const double s = par[4];
  return digt(tt, B / s + 0.5 * A / s, v / s, 0.5 * A / s);
}

inline double prdm_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[0])) return 0.0;
  const double tt = t - par[3];
  if (tt <= 0.0) return 0.0;
  const double v = par[0];
  const double B = par[1];
  const double A = par[2];
  const double s = par[4];
  return pigt(tt, B / s + 0.5 * A / s, v / s, 0.5 * A / s);
}

inline double dlnr_scalar(double t, const double* par, void* /*ctx_*/) {
  const double m  = par[0];
  const double s  = par[1];
  const double t0 = par[2];

  if (R_IsNA(m)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return R::dlnorm(tt, m, s, FALSE);
}

inline double plnr_scalar(double t, const double* par, void* /*ctx_*/) {
  const double m  = par[0];
  const double s  = par[1];
  const double t0 = par[2];

  if (R_IsNA(m)) return 0.0;
  const double tt = t - t0;
  if (tt <= 0.0) return 0.0;
  return R::plnorm(tt, m, s, TRUE, FALSE);
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double pdf = dlba_norm(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], pd);
    out[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : min_ll;
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double pdf = digt(tt, B_[i] / s_[i] + 0.5 * A_[i] / s_[i],
                            v_[i] / s_[i], 0.5 * A_[i] / s_[i]);
    out[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : min_ll;
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = 0.0; continue; }
    const double cdf = pigt(tt, B_[i] / s_[i] + 0.5 * A_[i] / s_[i],
                            v_[i] / s_[i], 0.5 * A_[i] / s_[i]);
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
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(m_[i]) || !isok[i]) { out[i] = min_ll; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0) { out[i] = min_ll; continue; }
    const double pdf = R::dlnorm(tt, m_[i], s_[i], FALSE);
    out[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : min_ll;
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
    const double cdf = R::plnorm(tt, m_[i], s_[i], TRUE, FALSE);
    if (cdf >= 1.0) { out[i] = min_ll; continue; }
    out[i] = (cdf <= 0.0) ? 0.0 : std::log1p(-cdf);
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

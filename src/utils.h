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
                                        Rcpp::NumericVector* pw_out = nullptr
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

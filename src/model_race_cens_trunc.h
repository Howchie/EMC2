#ifndef MODEL_RACE_CENS_TRUNC_H
#define MODEL_RACE_CENS_TRUNC_H

#include <RcppArmadillo.h>
#include "utility_functions.h" // For any existing utilities we might use

// Function pointer types for model's dfun and pfun
// Added LogicalVector is_ok for bound checking and void* context for extra args
typedef Rcpp::NumericVector (*RacePdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
                                          Rcpp::LogicalVector is_ok,
                                          bool log_p,
                                          void* model_specific_context);
typedef Rcpp::NumericVector (*RaceCdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
                                          Rcpp::LogicalVector is_ok,
                                          bool log_p,
                                          void* model_specific_context);

// Forward declaration for the main exported function
// Returns a single double: the total log-likelihood for one particle.

// --- Internal Helper Function Declarations (not directly exported to R via Rcpp::export) ---

// Integrand function equivalent to R's f_race_integrand
// p_trial_this_winner_first: (n_accumulators x n_model_params)
// model_dfun, model_pfun: function pointers to the model's density and CDF
double f_race_integrand_cpp(
    double t,
    const Rcpp::NumericMatrix& p_trial_this_winner_first, // Parameters for the winner (1st row) and losers
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    void* model_specific_context // Added context
);

Rcpp::NumericVector f_race_integrand_batch_cpp(
    const Rcpp::NumericVector& rts_batch,
    const Rcpp::List& pars_all_trials_ordered, // List of NumericMatrix, each n_acc x n_params
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    void* model_specific_context // Added context
);

double c_log_likelihood_race_cens_trunc(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, covering all dadm rows for that particle
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions, structured for all accumulators
    RacePdfFun model_dfun,                  // Pointer to the model's PDF adapter function
    RaceCdfFun model_pfun,                  // Pointer to the model's CDF adapter function
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector& ok_params,   // Parameter validity for each row in 'pars' matrix
    int n_acc,                              // Number of accumulators in the race (must be > 0 if data exists)
    const Rcpp::IntegerVector& expand_vec,  // Vector for expanding unique LLs to full trial count
    void* model_context_for_funcs           // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
);
// Numerical integration helper
// k_winner_idx: 1-based index of the accumulator treated as winner
// p_all_acc: parameters for all accumulators for the current trial condition
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    double low,
    double upp,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    double epsilon,
    void* model_specific_context
);

// Truncation correction factor helper
// Returns the factor itself (not log)
double get_trunc_corr_factor_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    double LT, double UT,
    int n_acc,
    double epsilon,
    void* model_specific_context
);

// Helper to order parameters (winner first)
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc, // n_acc x n_params
    int k_idx, // 1-based index of the winner
    int n_acc
);


#endif // MODEL_RACE_CENS_TRUNC_H
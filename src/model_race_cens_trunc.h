#ifndef MODEL_RACE_CENS_TRUNC_H
#define MODEL_RACE_CENS_TRUNC_H

#include <RcppArmadillo.h>
#include "utility_functions.h" // For any existing utilities we might use

// Function pointer types for model's dfun and pfun
typedef Rcpp::NumericVector (*RacePdfFun)(Rcpp::NumericVector rt, Rcpp::NumericMatrix pars, bool log_p);
typedef Rcpp::NumericVector (*RaceCdfFun)(Rcpp::NumericVector rt, Rcpp::NumericMatrix pars, bool log_p);

// Forward declaration for the main exported function
// Returns a single double: the total log-likelihood for one particle.
double c_log_likelihood_race_cens_trunc(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, all unique trial conditions
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions
    RacePdfFun model_dfun,                  // Pointer to the model's C++ PDF
    RaceCdfFun model_pfun,                  // Pointer to the model's C++ CDF
    double min_ll,
    const Rcpp::LogicalVector& ok_params,   // Parameter validity for this particle's pars matrix
    int n_acc,                              // Number of accumulators in the race
    const Rcpp::IntegerVector& expand_vec   // Vector for expanding unique LLs to full trial count
);


// --- Internal Helper Function Declarations (not directly exported to R via Rcpp::export) ---
// Signatures of helpers remain the same as they operate on a per-unique-trial basis.

// Integrand function equivalent to R's f_race_integrand
// p_trial_this_winner_first: (n_accumulators x n_model_params)
// model_dfun, model_pfun: function pointers to the model's density and CDF
double f_race_integrand_cpp(
    double t,
    const Rcpp::NumericMatrix& p_trial_this_winner_first,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc // Added n_acc
);

// Numerical integration helper
// k_winner_idx: 1-based index of the accumulator treated as winner
// p_all_acc: parameters for all accumulators for the current trial condition
double integrate_for_kth_winner_cpp(
    int k_winner_idx,
    const Rcpp::NumericMatrix& p_all_acc,
    double low,
    double upp,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    double epsilon // Added epsilon for integration tolerance
);

// Truncation correction factor helper
// Returns the factor itself (not log)
double get_trunc_corr_factor_for_kth_winner_cpp(
    int k_winner_idx,
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    double LT, double UT,
    int n_acc,
    double epsilon // Added epsilon
);

// Helper to order parameters (winner first)
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    int k_idx,
    int n_acc
);


#endif // MODEL_RACE_CENS_TRUNC_H

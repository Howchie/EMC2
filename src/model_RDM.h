#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include "utility_functions.h"
#include "composite_functions.h"
#include "gaussian.h"
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
double dwald(double t, double b, double mu, double sigma=1, double A=0, bool log_out=false) {
  double x_lo = 0;
  double x_hi = x_lo+A;
  
  const double var = sigma * sigma * t;
  if (var <= FPM_EPSILON) return 0.0;
  
  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;
  
  // Standard Wald pdf for a point start to a fixed boundary, given a drift and  
  if (span <= FPM_EPSILON) { // Fallback to original point-start formula
    const double d = b - x_hi;
    if (d <= 0) return (t==0.0) ? std::numeric_limits<double>::infinity() : 0.0;
    const double delta = d - mu * t;
    const double pdf_val = d/t * Gstar(var,delta);
    return log_out ? std::log(pdf_val) : pdf_val;
  }
  
  // Integrate the PDF across the uniform start range
  // The integral of (b-x) * N(x | b-mu*t, var) dx simplifies to a form
  // involving the Normal PDF and CDF.
  const double mu_new = b - mu * t;
  
  // Evaluate the Normal PDF N(mu_new, var) at the boundaries
  const double pdf_hi = Gstar(var, x_hi - mu_new);
  const double pdf_lo = Gstar(var, x_lo - mu_new);
  
  // Evaluate the Normal CDF Integral ∫ N(mu_new, var) dx over the boundaries
  const double cdf_integral = Gstar_Integral(var, mu_new, x_lo, x_hi);
  
  // The integral simplifies to μt * ∫N(x)dx + var * [N(x)]_lo^hi
  const double term1 = mu * t * cdf_integral;
  const double term2 = var * (pdf_hi - pdf_lo);
  
  const double pdf_val = ((1.0 / t) * (term1 + term2))/span;
  
  return log_out ? std::log(pdf_val) : pdf_val;
}

// Analytic solution for Geometric Brownian Motion -- Inverse Gaussian in log-space, with uniform start-point variability in PHYSICAL space (becomes exponential in log space)
//TODO figure out best way to handle mu-sigma^2 < .5 (return Inf but can we bound the parameters in the sampler?)
// [[Rcpp::export]]
double dgbm(double t, double b, double mu, double sigma, double A, bool log_out=false) {
  double x_lo = 0; // this is equivalent to starting at X0=1 in physical space (log(1) 
  double x_hi = std::log1p(A); // in physical space 1+A
  double log_b = std::log(b);
  double log_mu = mu - 0.5 * sigma * sigma;
  const double var = sigma * sigma * t; // Diffusion variance at time t
  if (var <= FPM_EPSILON) return 0.0;
  if (x_hi < x_lo) std::swap(x_lo, x_hi); // shouldn't be possible
  const double span = x_hi - x_lo;
  
  if (span <= FPM_EPSILON) { // Fallback to original point-start formula, expressed in terms of a heat kernel
    const double d = log_b - x_hi; // distance to boundary from start
    if (d <= 0) return (t==0.0) ? std::numeric_limits<double>::infinity() : 0.0;
    const double delta = d - log_mu * t; // Drift-adjusted mean of the heat kernel (distance form)
    const double pdf_val = d/t * Gstar(var,delta);
    return log_out ? std::log(pdf_val) : pdf_val;
  }
  
  const double norm_const = std::exp(x_hi) - std::exp(x_lo);
  if (norm_const <= FPM_EPSILON) return 0.0;
  
  // We complete the square, resulting in a new Gaussian
  const double mu_new_exp = log_b - log_mu * t + var;
  const double exp_factor = std::exp(log_b - log_mu * t + var / 2.0);
  
  // Evaluate the NEW Normal PDF N(mu_new_exp, var) at the boundaries
  const double pdf_hi = Gstar(var, x_hi - mu_new_exp);
  const double pdf_lo = Gstar(var, x_lo - mu_new_exp);
  
  // Evaluate the NEW Normal CDF Integral
  const double cdf_integral = Gstar_Integral(var, mu_new_exp, x_lo, x_hi);
  
  // The structure of the integral solution is the same, just with the new mean
  const double integral_term1 = (log_mu * t - var) * cdf_integral;
  const double integral_term2 = var * (pdf_hi - pdf_lo);
  const double integral_result = integral_term1 + integral_term2;
  
  const double pdf_val = (1.0 / (norm_const * t)) * exp_factor * integral_result;
  return log_out ? std::log(pdf_val) : pdf_val;
}

/**
 * @brief CDF for Wald (ABM) with uniform start-point variability.
 * Corresponds to 'dwald'.
 */
// [[Rcpp::export]]
double pwald(double t, double b, double mu, double sigma, double A, bool log_out=false) {
  if (t <= FPM_EPSILON) return log_out ? -std::numeric_limits<double>::infinity() : 0.0;
  
  double x_lo = 0;
  double x_hi = x_lo + A;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;
  
  // Fallback to standard Wald CDF
  if (span <= FPM_EPSILON) { // point-start fallback
    const double sqrt_var = sigma * std::sqrt(t);
    const double dist = b - x_hi;
    if (sqrt_var <= FPM_EPSILON) return (mu > 0) ? 1.0 : 0.0;
    
    const double term1_arg = (mu * t - dist) / sqrt_var;
    const double term2_arg = (-mu * t - dist) / sqrt_var;
    
    // Use  Gstar_CDF for the Standard Normal CDF (Φ)
    const double cdf1 = gaussian_cdf(term1_arg, 0.0, 1.0); 
    const double cdf2 = gaussian_cdf(term2_arg, 0.0, 1.0);
    
    const double exp_term = std::exp(2.0 * mu * dist / (sigma * sigma));
    double cdf_val =0.0;
    cdf_val = cdf1 + exp_term * cdf2;
    if (!std::isfinite(exp_term)) cdf_val = 1.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }
  
  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) { // t is effectively zero
    double cdf_val = (mu > 0 && b <= x_hi) ? 1.0 : 0.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }
  
  const double k = 2.0 * mu / (sigma * sigma);
  
  // First term integral: ∫ Φ((x + (μt - b))/st) dx
  const double a1 = 1.0 / st;
  const double c1 = (mu * t - b) / st;
  const double term1_integral = integrate_gaussian_cdf(a1, c1, x_lo, x_hi);
  
  // Second term integral: ∫ exp(k*(b-x)) * Φ((x - (μt+b))/st) dx
  //  = exp(kb) * ∫ exp(-kx) * Φ(a2*x + c2) dx
  const double a2 = 1.0 / st;
  const double c2 = (-mu * t - b) / st;
  const double term2_integral = std::exp(k * b) * integrate_exp_times_normal_cdf(-k, a2, c2, x_lo, x_hi);
  
  double cdf_val = (term1_integral + term2_integral) / span;
  
  if (cdf_val < 0.0) cdf_val = 0.0;
  if (cdf_val > 1.0) cdf_val = 1.0;
  
  return log_out ? std::log(cdf_val) : cdf_val;
}


/**
 * @brief CDF for Geometric Brownian Motion with physical-space uniform start.
 */
// [[Rcpp::export]]
double pgbm(double t, double b, double mu, double sigma, double A, bool log_out=false) {
  if (t <= FPM_EPSILON) return log_out ? -std::numeric_limits<double>::infinity() : 0.0;
  
  // --- Setup MUST match dgbm ---
  const double x_lo = 0.0;         // log(1)
  const double x_hi = std::log1p(A); // log(1+A)
  const double log_b = std::log(b);
  const double log_mu = mu - 0.5 * sigma * sigma; // GBM drift in log-space
  
  const double norm_const = std::exp(x_hi) - std::exp(x_lo); // This is just A
  
  if (norm_const <= FPM_EPSILON) { // point-start fallback
    const double sqrt_var = sigma * std::sqrt(t);
    const double dist = log_b - x_hi;
    if (sqrt_var <= FPM_EPSILON) return (log_mu > 0) ? 1.0 : 0.0;
    
    const double term1_arg = (log_mu * t - dist) / sqrt_var;
    const double term2_arg = (-log_mu * t - dist) / sqrt_var;
    
    // Use  Gstar_CDF for the Standard Normal CDF (Φ)
    const double cdf1 = gaussian_cdf(term1_arg, 0.0, 1.0); 
    const double cdf2 = gaussian_cdf(term2_arg, 0.0, 1.0);
    
    const double exp_term = std::exp(2.0 * log_mu * dist / (sigma * sigma));
    double cdf_val =0.0;
    cdf_val = cdf1 + exp_term * cdf2;
    if (!std::isfinite(exp_term)) cdf_val = 1.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }
  
  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) {
    double cdf_val = (log_mu > 0 && log_b <= x_hi) ? 1.0 : 0.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }
  
  const double k_log = 2.0 * log_mu / (sigma * sigma);
  
  // We need to compute (1/norm_const) * ∫[x_lo, x_hi] F(t, x0) * exp(x0) dx0
  // F(t, x0) is the point-start CDF with dist = log_b - x0, drift = log_mu
  
  // Term 1: ∫ exp(x) * Φ((x + (log_mu*t - log_b))/st) dx
  const double k1 = 1.0;
  const double a1 = 1.0 / st;
  const double c1 = (log_mu * t - log_b) / st;
  const double term1_integral = integrate_exp_times_normal_cdf(k1, a1, c1, x_lo, x_hi);
  
  // Term 2: ∫ exp(x) * exp(k_log*(log_b-x)) * Φ((x - (log_mu*t+log_b))/st) dx
  //  = exp(k_log*log_b) * ∫ exp((1-k_log)x) * Φ(a2*x + c2) dx
  const double k2 = 1.0 - k_log;
  const double a2 = 1.0 / st;
  const double c2 = (-log_mu * t - log_b) / st;
  const double exp_factor = std::exp(k_log * log_b);
  
  const double term2_integral = exp_factor * integrate_exp_times_normal_cdf(k2, a2, c2, x_lo, x_hi);
  
  double cdf_val = (term1_integral + term2_integral) / norm_const;
  
  if (cdf_val < 0.0) cdf_val = 0.0;
  if (cdf_val > 1.0) cdf_val = 1.0;
  
  return log_out ? std::log(cdf_val) : cdf_val;
}

/**
 * Calculates the log-PDF for the Shifted Wald Truncated Normal (SWTN) model.
 *
 * This is the full generalization that includes a free diffusion
 * coefficient (s) where Steingroever's published results held s = 1.
 * NB: Steingroever et al.'s JAGS code included s through a precision parameterisation but fixed it at 1 for implementation
 *
 * @param t_adj The adjusted time (t - theta).
 * @param threshold The boundary (threshold > 0).
 * @param mu_drift The mean of the drift-rate distribution.
 * @param sv The standard deviation of the drift-rate distribution.
 * @param s The diffusion coefficient
 * @param log_out Whether to return the log-PDF (true) or PDF (false).
 * @return The log-PDF or PDF value.
 */
// [[Rcpp::export]]
double dswtn(double t_adj, double threshold, double mu_drift, double sv, double s = 1.0, double c = 0.0, bool log_out=false) {
  
  // --- 1. Input Validation ---
  if (t_adj <= 1e-10) return R_NegInf; // log(0)
  if (threshold <= 1e-10) return R_NegInf; // No boundary
  if (sv < 0) return R_NaN;           // Drift variance < 0
  if (s <= 1e-10) return R_NegInf; // Diffusion variance must be > 0
  
  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
    return dwald(t_adj, threshold, mu_drift, s, 0.0, log_out);
  }
  
  double v = sv * sv; // Variance of the drift variability distribution
  double s_sq = s * s; // Diffusion variance
  double tv = t_adj * v;    // t' * drift_var
  
  // This is the common denominator term from the derivation: (t'*v + s)
  double den_common = tv + s_sq;
  double log_den_common = std::log(den_common);
  
  // Normalization for P(xi > c)
  // This is for the *original* drift distribution N(mu_drift, v)
  // log( 1 / P(xi > c) ) = -log(P(xi > c))
  double log_prob_mu_gt_c = R::pnorm(c, mu_drift, sv, false, true);
  if (std::isinf(log_prob_mu_gt_c) && log_prob_mu_gt_c < 0) return R_NegInf;
  double term_log_norm = -log_prob_mu_gt_c;
  
  //The Mixture Core
  double term_log_threshold = std::log(threshold);
  double term_log_denom = -0.5 * (std::log(M_PI) + std::log(2.0) + 3.0 * std::log(t_adj) + log_den_common);
  
  //The Exponential Term ---
  double term_log_exp = -(std::pow(threshold - mu_drift * t_adj, 2.0)) / (2.0 * t_adj * den_common);
  
  // The Truncation Integral ---
  // This is log( P(xi_new > c) ) for the *new* Gaussian
  double mu_new = (threshold * v + mu_drift * s_sq) / den_common;
  double sigma_new = std::sqrt(s_sq * v / den_common);
  
  // log(P(xi_new > c))
  double term_log_int = R::pnorm(c, mu_new, sigma_new, false, true);
  
  double log_pdf_val = term_log_threshold + term_log_denom + term_log_norm + term_log_exp + term_log_int;
  
  if (std::isnan(log_pdf_val)) return R_NegInf;
  return log_out ? log_pdf_val : std::exp(log_pdf_val);
}

/**
 * Calculates the CDF for the Shifted Wald Truncated Normal (SWTN) model.
 *
 * This is the full, correct generalization that includes a free diffusion
 * coefficient (s, i.e., sigma).
 *
 * @param t_adj The adjusted time (t - theta).
 * @param threshold The boundary (threshold > 0).
 * @param mu_drift The mean of the drift-rate distribution.
 * @param sv The standard deviation of the drift-rate distribution.
 * @param s The diffusion coefficient (sigma).
 * @param c The lower truncation bound for the drift rate (default 0).
 * @param log_out Whether to return the log-CDF (true) or CDF (false).
 */
// [[Rcpp::export]]
double pswtn(double t_adj, double threshold, double mu_drift, double sv, double s = 1.0, double c = 0.0, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (sv < 0) return R_NaN;
  
  // Handle sv == 0 case (becomes standard Wald with sigma)
  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    // Call the standard Wald CDF (which must accept 's')
    return pwald(t_adj, threshold, mu_drift, s, 0.0, log_out);
  }
  
  double v = sv * sv;     // v = drift rate variance (sigma_xi^2)
  double s_sq = s * s;  // s_sq = diffusion variance (sigma^2)
  double tv = t_adj * v;  // t' * v
  double threshold_2 = threshold * threshold;
  
  // --- Normalization Term (for original drift dist) ---
  // P(xi > c) for xi ~ N(mu_drift, v)
  double log_prob_mu_gt_c = R::pnorm(c, mu_drift, sv, false, true); // log.p=TRUE
  if (std::isinf(log_prob_mu_gt_c) && log_prob_mu_gt_c < 0) return R_NegInf;
  
  // --- Term 1 (Standard) ---
  // This is from integrating the first part of the Wald CDF
  double s_denom_new = std::sqrt(t_adj * (s_sq + tv));
  double rho_new = (std::sqrt(t_adj) * sv) / std::sqrt(s_sq + tv);
  
  double h1_new = (mu_drift * t_adj - threshold) / s_denom_new;
  double k1_new = (mu_drift - c) / sv; // <-- 'c' is incorporated here
  
  double term1;
  if (std::fabs(rho_new) > 0.9999 || std::fabs(h1_new) > 8.0 || std::fabs(k1_new) > 8.0) {
    term1 = norm_cdf_2d(h1_new, k1_new, rho_new);
  } else {
    term1 = norm_cdf_2d_fast(h1_new, k1_new, rho_new);
  }
  double log_A = (term1 > 1e-300) ? std::log(term1) : R_NegInf;
  
  // --- Term 2 (Reflected) ---
  // This is from integrating the second (exp() * Phi()) part
  
  // New exponential term: log( exp( (2*a*mu/s^2) + (2*a^2*v/s^4) ) )
  double log_exp_term = (2.0 * threshold * mu_drift / s_sq) + (2.0 * threshold_2 * v / (s_sq * s_sq));
  
  // New "reflected" mean: mu' = mu + 2*threshold*(v/s^2)
  double mu_p_new = mu_drift + 2.0 * threshold * v / s_sq;
  
  double h2_new = (-mu_p_new * t_adj - threshold) / s_denom_new; // denom is the same
  double k2_new = (mu_p_new - c) / sv; // <-- 'c' is incorporated here
  
  double term2;
  if (std::fabs(rho_new) > 0.9999 || std::fabs(h2_new) > 8.0 || std::fabs(k2_new) > 8.0) {
    term2 = norm_cdf_2d(h2_new, k2_new, -rho_new); // rho is negative
  } else {
    term2 = norm_cdf_2d_fast(h2_new, k2_new, -rho_new);
  }
  double log_B = (term2 > 1e-300) ? (log_exp_term + std::log(term2)) : R_NegInf;
  
  // --- Combine Terms ---
  // log( (A + B) / P(xi > c) ) = log(A + B) - log(P(xi > c))
  double log_numerator = log_sum_exp(log_A, log_B);
  double log_cdf_val = log_numerator - log_prob_mu_gt_c;
  
  // --- Error Checking and Return ---
  if (std::isnan(log_cdf_val)) {
    return log_out ? R_NegInf : 0.0;
  }
  
  double cdf_val = std::exp(log_cdf_val);
  if (cdf_val < 0.0) return log_out ? R_NegInf : 0.0;
  if (cdf_val > 1.0) return log_out ? 0.0 : 1.0;
  
  return log_out ? log_cdf_val : cdf_val;
}

// Top-level PDF for RDM_SWTN model
// Parameters b, A, mu_drift, sv are assumed to be already scaled by s if applicable. This integrates the dswtn function across b->b+A using a 20-point Gauss-Legendre approximation.
// [[Rcpp::export]]
double drdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s=1.0, double c=0.0, int n_gauss_nodes  = 20, bool log_out=false) {
  
  if (t_adj <= 1e-10) return 0.0;
  if (b <= 1e-7) return R_NegInf; // Max threshold non-positive
  bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
  bool no_drift_var = (sv < 1e-7);
  
  if (no_A_var && no_drift_var) {
    // Case 1: Simple Wald. Threshold is b. Drift is mu_drift.
    if (mu_drift <= 1e-10 || b <= 1e-10) return R_NegInf; // No positive drift to positive boundary
    return dwald(t_adj, b, mu_drift,s,0.0, log_out);
  } else if (!no_A_var && no_drift_var) {
    // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
    return dwald(t_adj, b, mu_drift, s, A, log_out);
  } else if (no_A_var && !no_drift_var) {
    // Case 3: SWTN with fixed threshold b (like original dswtn). Use Steingrover et al (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability.
    return dswtn(t_adj, b, mu_drift, sv,s,c,log_out);
  } else {
    // Case 4: Full model - SWTN with RDM-style SPV.
    // Integrate dswtn(t_adj, actual_k, mu_drift, sv)
    // where actual_k ~ U(b-A, b).
    // Note: 1e-7 has already confirmed A is not zero for this 'else' block.
    // So, no_A_var is false here.
    // If max threshold b is non-positive, result is -Inf. Shouldn't be possible.
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl;
    if (n_nodes == 20) {
      gl = gl20;
    } else {
      gl = Rcpp::as<Rcpp::List>(gauss_quad(n_nodes, "legendre"));
    }
    const Rcpp::NumericVector gl_nodes = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    
    // Map gl_nodes to [b, b+A]
    Rcpp::NumericVector k_nodes = b - 0.5 * A + 0.5*A*gl_nodes;
    // 3. Evaluate CDF at each (t, k) pair and integrate
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      double pdf_val = dswtn(t_adj, k_nodes[j], mu_drift, sv, s, c, false); // integral cannot take log_pdf
      integral      += gl_weights[j] * pdf_val;
    }
    // Scale for Gauss-Legendre over [b, b+A] and divide by A (uniform pdf)
    // This reduces to integral * 0.5 but it's handy to see it written out properly
    double out = (integral * (A / 2.0)) / A;
    
    if (out < 0.0 || std::isnan(out)) {
      return log_out ? R_NegInf : 0;
    }
    return log_out ? std::log(out) : out;
  }
}

// Top-level CDF for RDM_SWTN model
// Parameters b, A, mu_drift, sv are assumed to be already scaled by s if applicable. This integrates the pswtn function across b->b+A using a 20-point Gauss-Legendre approximation.
// [[Rcpp::export]]
double prdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv,double s=1.0, double c=0.0, int n_gauss_nodes  = 20, bool log_out=false)
{
  if (t_adj <= 1e-10) return 0.0;
  if ((b) <= 1e-7) return 1.0;
  bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
  bool no_drift_var = (sv < 1e-7);
  
  if (no_A_var && no_drift_var) {
    // Case 1 standard Wald (neither variability parameter)
    if (mu_drift <= 1e-7 && b > 0) return 0.0;
    return pwald(t_adj, b, mu_drift, s, 0.0, log_out);
  } else if (!no_A_var && no_drift_var) {
    // Case 2: standard RDM with start-point variability only.
    return pwald(t_adj, b, mu_drift, s, A, log_out);
  } else if (no_A_var && !no_drift_var) {
    // Case 3: SWTN with fixed threshold b (like original dswtn). Starts with Steingrover et al's (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability - but they derived only the pdf with a fixed mu_drift so here we're integrating that numerically with mu_drift =0->Inf
    if (b <= 1e-9 && mu_drift>1e-9) return 1.0;
    return pswtn(t_adj, b, mu_drift, sv, s, c, log_out);
  } else {
    
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl;
    if (n_nodes == 20) {
      gl = gl20;
    } else {
      gl = Rcpp::as<Rcpp::List>(gauss_quad(n_nodes, "legendre"));
    }
    const Rcpp::NumericVector gl_nodes = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    
    // Map gl_nodes to [b, b+A]
    Rcpp::NumericVector k_nodes = b - 0.5 * A + 0.5*A*gl_nodes;
    // 3. Evaluate CDF at each (t, k) pair and integrate
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      double cdf_val = pswtn(t_adj, k_nodes[j], mu_drift, sv, s, c, false); // integral cannot take log_cdf
      integral      += gl_weights[j] * cdf_val;
    }
    // Scale for Gauss-Legendre over [b, b+A] and divide by A (uniform pdf)
    // This reduces to integral * 0.5 but it's handy to see it written out properly
    double out = (integral * (A / 2.0)) / A;
    
    if (out < 0.0 || std::isnan(out)) {
      return log_out ? R_NegInf : 0;
    }
    if (out > 1.0) {
      return log_out ? 0 : 1;
    }
    return log_out ? std::log(out) : out;
  }
  
}
/*
#' Calculates the probability density function (PDF) for the
#' Wald distribution with trial-varying normally distributed drift rates (SWTN).
#' This distribution is sometimes referred to as an IG_TN (Inverse Gaussian
#' mixed with Truncated Normal drift rates).
#'
#' @param t Vector of reaction times.
#' @param b Vector of base threshold parameters (lower bound of SPV uniform range).
#' @param A Vector of start-point variability range parameters (width of SPV uniform range, actual_k ~ U(b, b+A)). Must be non-negative.
#' @param v Vector of mean drift rate parameters.
#' @param sv Vector of standard-deviation parameters for the drift rate distribution. Must be non-negative.
#' @param t0 Vector of non-decision time parameters.
#' @return A numeric vector of log-density values.
#' @details This function models RTs based on a Wald process where:
#'   1. The start-point (threshold `k`) varies from trial to trial, uniformly distributed over `[b/s, (b+A)/s]`.
#'   2. The drift rate (`xi`) varies from trial to trial, drawn from a normal
#'      distribution N(v/s, sv), truncated at 0.
#'   Handles vectorization of parameters. */
 // [[Rcpp::export]]
 NumericVector dSWTNspv(NumericVector t, NumericVector v, NumericVector b, NumericVector A, NumericVector t0, NumericVector sv, NumericVector s=1.0, NumericVector c=0.0, int n_gauss_nodes=20, bool log_out=false) { 
   const int n = t.size();
   NumericVector pdf(n);
   
   auto pick = [](const NumericVector &vec, int i) -> double {
     return vec.size() == 1 ? vec[0] : vec[i];
   };
   
   for (int i = 0; i < n; ++i) {
     const double shifted_t = t[i] - pick(t0, i);
     if (shifted_t <= 0.0) {
       pdf[i] = log_out ? R_NegInf : 0.0;
       continue;
     }
     pdf[i] = drdmswtn(shifted_t,
                       pick(b, i),
                       pick(v, i),
                       pick(A, i),
                       pick(sv, i),
                       pick(s, i),
                       pick(c, i),
                       n_gauss_nodes,
                       log_out);
   }
   return pdf;
 }

/*
#' Cumulative distribution function for the RDM_SWTN model
#'
#' Calculates the CDF for the RDM_SWTN model, which incorporates RDM-style
#' start-point variability and SWTN-style trial-varying drift rates.
#'
#' @param t Vector of reaction times.
#' @param b Vector of base threshold parameters.
#' @param A Vector of start-point variability range parameters.
#' @param v Vector of mean drift rate parameters.
#' @param sv Vector of variance parameters for drift rate.
#' @param t0 Vector of non-decision time parameters.
#' @param spv_abs_err,spv_rel_err,spv_max_eval Control integration for start-point variability.
#' @param dr_abs_err,dr_rel_err,dr_max_eval Control integration for drift-rate variability (inner integral).
#' @return A numeric vector of CDF values.
#' @details Parameters are handled similarly to `dRDM_SWTN`, including scaling by `s`. */
 // [[Rcpp::export]]
 NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector b, NumericVector A, NumericVector t0, NumericVector sv,  NumericVector s=1.0, NumericVector c=0.0, int n_gauss_nodes=20, bool log_out=false) { 
   const int n = t.size();
   NumericVector cdf(n);
   
   auto pick = [](const NumericVector &vec, int i) -> double {
     return vec.size() == 1 ? vec[0] : vec[i];
   };
   
   for (int i = 0; i < n; ++i) {
     const double shifted_t = t[i] - pick(t0, i);
     if (shifted_t <= 0.0) {
       cdf[i] = log_out ? R_NegInf : 0.0;
       continue;
     }
     cdf[i] = prdmswtn(shifted_t,
                       pick(b, i),
                       pick(v, i),
                       pick(A, i),
                       pick(sv, i),
                       pick(s, i),
                       pick(c, i),
                       n_gauss_nodes,
                       log_out);
   }
   return cdf;
 }

// [[Rcpp::export]]
NumericVector drdmswtn_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4, sv = 5
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        double b = (pars(i,1) + pars(i,2)) / pars(i,4); // convert B,A to scaled b
        out[k] = drdmswtn(rts[i] - pars(i,3),
                          b,
                          pars(i,0)/pars(i,4),
                          pars(i,2)/pars(i,4),
                          pars(i,5)/pars(i,4),
                          1.0,      // s already absorbed above
                          0.0,      // c
                          20,       // n_gauss_nodes
                          log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector prdmswtn_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4, sv = 5
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        double b = (pars(i,1) + pars(i,2)) / pars(i,4); // convert B,A to scaled b
        out[k] = prdmswtn(rts[i] - pars(i,3),
                          b,
                          pars(i,0)/pars(i,4),
                          pars(i,2)/pars(i,4),
                          pars(i,5)/pars(i,4),
                          1.0,      // s already absorbed above
                          0.0,      // c
                          20,       // n_gauss_nodes
                          log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector dWALDspv(NumericVector t, NumericVector b, NumericVector v, NumericVector sigma, NumericVector A, NumericVector t0, bool log_out=false) { 
  int n = t.size();
  NumericVector pdf(n);
  
  auto pick = [](const NumericVector &vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double shifted_t = t[i] - pick(t0, i);
    if (shifted_t <= 0.0) {
      pdf[i] = log_out ? R_NegInf : 0.0;
      continue;
    }
    pdf[i] = dwald(shifted_t,
                   pick(b, i),
                   pick(v, i),
                   pick(sigma, i),
                   pick(A, i),
                   log_out);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dGBMspv(NumericVector t, NumericVector b, NumericVector v, NumericVector sigma, NumericVector A, NumericVector t0, bool log_out=false) { 
  const int n = t.size();
  NumericVector pdf(n);
  
  auto pick = [](const NumericVector &vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  
  for (int i = 0; i < n; ++i) {
    const double shifted_t = t[i] - pick(t0, i);
    if (shifted_t <= 0.0) {
      pdf[i] = log_out ? R_NegInf : 0.0;
      continue;
    }
    pdf[i] = dgbm(shifted_t,
                  pick(b, i),
                  pick(v, i),
                  pick(sigma, i),
                  pick(A, i),
                  log_out);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pWALDspv(NumericVector t, NumericVector b, NumericVector v, NumericVector sigma, NumericVector A, NumericVector t0, bool log_out=false) { 
  int n = t.size();
  NumericVector cdf(n);
  
  auto pick = [](const NumericVector &vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double shifted_t = t[i] - pick(t0, i);
    if (shifted_t <= 0.0) {
      cdf[i] = log_out ? R_NegInf : 0.0;
      continue;
    }
    cdf[i] = pwald(shifted_t,
                   pick(b, i),
                   pick(v, i),
                   pick(sigma, i),
                   pick(A, i),
                   log_out);
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector pGBMspv(NumericVector t, NumericVector b, NumericVector v, NumericVector sigma, NumericVector A, NumericVector t0, bool log_out=false) { 
  const int n = t.size();
  NumericVector cdf(n);
  
  auto pick = [](const NumericVector &vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  
  for (int i = 0; i < n; ++i) {
    const double shifted_t = t[i] - pick(t0, i);
    if (shifted_t <= 0.0) {
      cdf[i] = log_out ? R_NegInf : 0.0;
      continue;
    }
    cdf[i] = pgbm(shifted_t,
                  pick(b, i),
                  pick(v, i),
                  pick(sigma, i),
                  pick(A, i),
                  log_out);
  }
  return cdf;
}

#endif

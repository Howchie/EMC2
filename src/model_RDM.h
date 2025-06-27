#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include <random>
#include <algorithm>
#include "gauss.h" // For one_d struct and gauss_kronrod constants (if gauss_kronrod is used directly)
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include "utility_functions.h"
#include "truncated_normal.hpp"
#include <string>
using namespace Rcpp;

const double L_PI = std::log(M_PI);

// RDM
// [[Rcpp::export]]
double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. + t/mu), 0., 1., true, false);
  double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. - t/mu), 0., 1., true, false);

  return std::exp(2.0 * lambda / mu + std::log(p1)) + p2;
}

// Standard Wald PDF: f(t | drift_rate_xi, boundary_alpha) equivalent to pigt0 but with threshold/drift directly parameterized
// t must be > 0
// [[Rcpp::export]]
double dwald(double t, double boundary_alpha, double drift_rate_xi) {
    if (t <= 0) return 0.0;
    if (boundary_alpha <= 0) return (drift_rate_xi > 0) ? 1.0 : 0.0; // if boundary is 0 or less, instant absorption if drift is positive
    if (drift_rate_xi <= 1e-10) { // Effectively zero or negative drift away from boundary
        return 0.0; // Cannot reach positive boundary alpha
    }

    return boundary_alpha*std::exp(-std::pow((boundary_alpha-drift_rate_xi*t),2)/(2*t))/std::sqrt(2*M_PI*std::pow(t,3));
}

// Standard Wald CDF: F(t | drift_rate_xi, boundary_alpha) equivalent to digt0 but with threshold/drift directly parameterized
// t must be > 0
// [[Rcpp::export]]
double pwald(double t, double boundary_alpha, double drift_rate_xi) {
    if (t <= 0) return 0.0;
    if (boundary_alpha <= 0) return (drift_rate_xi > 0) ? 1.0 : 0.0; // if boundary is 0 or less, instant absorption if drift is positive
    if (drift_rate_xi <= 1e-10) { // Effectively zero or negative drift away from boundary
        return 0.0; // Cannot reach positive boundary alpha
    }

    double t_sqrt = std::sqrt(t);
    double term1_arg = (drift_rate_xi * t - boundary_alpha) / t_sqrt;
    double term2_arg = -(drift_rate_xi * t + boundary_alpha) / t_sqrt;

    double cdf_val = R::pnorm(term1_arg, 0.0, 1.0, true, false) +
                     std::exp(2.0 * boundary_alpha * drift_rate_xi) * R::pnorm(term2_arg, 0.0, 1.0, true, false);

    if (std::isnan(cdf_val) || cdf_val < 0.0) return 0.0;
    if (cdf_val > 1.0) return 1.0;
    return cdf_val;
}

// [[Rcpp::export]]
double digt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.) {
  //  return 0.;
  //}
  double lambda = k * k;
  double e;
  if (l == 0.) {
    e = -.5 * lambda / t;
  } else {
    double mu = k / l;
    e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
  }
  return std::exp(e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI));
}

// [[Rcpp::export]]
double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return pigt0(t, k, l); // uses the standard wald pdf (threshold, drift) for ease of interpretation
  }

  double sqt = std::sqrt(t);
  double lgt = std::log(t);
  double cdf;

  if (l < threshold){
    double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
    double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1;

    double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

    cdf = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
  } else {
    double t1a = std::exp(- .5 * std::pow(k - a - t * l, 2) / t);
    double t1b = std::exp(- .5 * std::pow(a + k - t * l, 2) / t);
    double t1 = std::exp(.5* (lgt - M_LN2 - L_PI)) * (t1a - t1b);

    double t2a = std::exp(2. * l * (k - a) + R::pnorm(- (k - a + t * l) / sqt, 0., 1., true, true));
    double t2b = std::exp(2. * l * (k + a) + R::pnorm(- (k + a + t * l) / sqt, 0., 1., true, true));
    double t2 = a + (t2b - t2a) / (2. * l);

    double t4a = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t4b = 2. * R::pnorm((k - a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;

    cdf = .5 * (t4 + t2 + t1) / a;
  }
  if (cdf < 0. || std::isnan(cdf)) {
    return 0.;
  }
  return cdf;
}

// [[Rcpp::export]]
double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10){
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return digt0(t, k, l);
  }
  double pdf;
  if (l < threshold){
    double term = std::exp(- (k - a) * (k - a) / (2. * t)) - std::exp(- (k + a) * (k + a) / (2. * t));
    pdf = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(term) - M_LN2 - std::log(a));
  } else {
    double sqt = std::sqrt(t);

    double t1a = - std::pow(a - k + t * l, 2) / (2. * t);
    double t1b = - std::pow(a + k - t * l, 2) / (2. * t);
    double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

    double t2a = 2. * R::pnorm((- k + a) / sqt + sqt * l, 0., 1., true, false) - 1.;
    double t2b = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);

    pdf = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
  }
  if (pdf < 0. || std::isnan(pdf)) {
    return 0.;
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector drdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = digt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

// [[Rcpp::export]]
NumericVector prdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = pigt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}


// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return cdf;
}

// START OF SWTN IMPLEMENTATION
// Helper struct for CDF integration parameters
struct pswtn_Params {
    double t_adj;        // Time (adjusted for non-decision time)
    double alpha;        // Threshold
    double mu_drift;     // Mean of the drift rate distribution
    double sigma_drift; // Variance of the drift rate distribution
};

// PDF of a normal distribution N(mu, sigma) truncated at lower_bound (typically 0)
// xi is the value at which to evaluate the PDF
/*
 [[Rcpp::export]]
double truncated_normal_pdf(double xi, double mu, double sigma, double lower_bound=0) {
    if (xi < lower_bound) return 0.0;
    if (sigma <= 1e-10) { // Effectively zero variance
        if (xi >= lower_bound && std::fabs(xi - mu) < 1e-7) return R_PosInf; // Spike at mu if mu is valid
        return 0.0;
    }
    double sigma = std::sqrt(sigma);

    // Prob that X > lower_bound, where X ~ N(mu, sigma)
    double prob_above_lower = R::pnorm(lower_bound, mu, sigma, false, false); // P(X > lower_bound)
                                                                           // same as 1.0 - R::pnorm(lower_bound, mu, sigma, true, false)
                                                                           // or R::pnorm((mu - lower_bound)/sigma, 0, 1, true, false)

    if (prob_above_lower < 1e-10) { // If virtually no mass above lower_bound
      // if mu is very far below lower_bound, dnorm will be tiny, division by tiny prob_above_lower could be huge
      // Check if xi is close to mu and mu is close to lower_bound in such a case.
      // This case is tricky. If mu < lower_bound and sigma is small, density is practically 0.
      return 0.0;
    }

    double pdf_val = R::dnorm(xi, mu, sigma, false) / prob_above_lower;
    return pdf_val;
}
 */

// Integrand for the SWTN CDF: dwald * truncated_normal_a_pdf (from library)
// This is the GSL-style function that will be integrated over xi.
double pswtn_integrand(double current_xi, void* p) {
    pswtn_Params* params = static_cast<pswtn_Params*>(p);

    // current_xi comes from GSL, lower bound is handled by qagiu's 'a' parameter.
    // No need to check current_xi < 0 here if 'a' is 0.
    // However, explicit check for safety or specific logic can remain if needed.
    if (current_xi < 0) { // Should ideally not be hit if gsl_integration_qagiu is called with lower_bound_xi = 0
        return 0.0;
    }

    double wc = pwald(params->t_adj, params->alpha, current_xi);
    double tn = truncated_normal_a_pdf(current_xi, params->mu_drift, params->sigma_drift, 0.0); // lower_bound = 0 for drift rate

    return wc * tn;
}

// Log-PDF for SWTN
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation
// [[Rcpp::export]]
double dswtn(double t_adj, double alpha, double mu_drift, double sigma_drift) {
    if (t_adj <= 1e-10) return R_NegInf; // log(0)
    if (alpha <= 1e-10) return R_NegInf; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative

    // Handle sigma_drift_sq == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return dwald(t_adj, alpha, mu_drift);
    }

    double lambda = 1.0; // Matches SWTN.cc
    double d = mu_drift;
    double v = std::pow(sigma_drift,2); // parameterize the RDM with std but formula uses variance

    // Log of the normalization constant for the drift rate xi ~ TN(d, v, lower=0)
    // log(1 / P(xi > 0)) = -log(P(xi > 0))
    // P(xi > 0) = pnorm(d/sqrt(v), 0, 1, true, false) if d > 0 is assumed for TN's mode.
    // Or more generally, P(xi > 0) for xi ~ N(d,v) is 1 - pnorm(0, d, sqrt(v), true, false)
    // = pnorm(0, d, sqrt(v), false, false)
    // = pnorm( (0-d)/sqrt(v), 0, 1, false, false) = pnorm(-d/sqrt(v), 0, 1, false, false)
    // = pnorm(d/sqrt(v), 0, 1, true, false)
    //double log_prob_xi_gt_0 = R::pnorm(d / std::sqrt(v), 0.0, 1.0, true, true); // log(Phi(d/sqrt(v)))
    //if (std::isinf(log_prob_xi_gt_0) && log_prob_xi_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero

    double term1 = std::log(alpha);
    double term2 = 0.5 * (std::log(lambda) - std::log(2.0) - L_PI - 3.0 * std::log(t_adj) - std::log(lambda * t_adj * v + 1.0));
    double term3 = -R::pnorm(d / std::sqrt(v), 0.0, 1.0, true, true); // This is -log(P(xi>0)) = log (1/P(xi>0))
    double term4 = -(lambda * std::pow(d * t_adj - alpha, 2.0)) / (2.0 * t_adj * (lambda * t_adj * v + 1.0));
    double term5 = R::pnorm((lambda * alpha * v + d) / std::sqrt(lambda * t_adj * std::pow(v,2) + v), 0.0, 1.0, true, true);
	double log_pdf_val = term1 + term2 + term3 + term4 + term5;
	
	/* Alternate un-logged version to match the Steingrover (2021) paper - untested, and their JAGS code used the log version so sticking with that)
	double term1 = alpha;
    double term2 = 0.5 * (std::log(lambda) - std::log(2.0) - L_PI - 3.0 * std::log(t_adj) - std::log(lambda * t_adj * v + 1.0));
    double term3 = 1 / pnorm(d / std::sqrt(v), 0, 1, 1, 0) ; // This is -log(P(xi>0)) = log (1/P(xi>0))
    double term4 = std::exp( - (lambda * std::pow(d * t - alpha, 2)) / (2 * t * (lambda * t * v + 1)) );
    double term5 = R::pnorm( (lambda * alpha * v + d) /(std::sqrt(lambda * t * std::pow(v, 2) + v) ), 0.0, 1.0, true, false);
	double pdf_val = std::log(term1 * term2 * term3 * term4 * term5); */

    if (std::isnan(log_pdf_val)) return R_NegInf; // Should be caught by specific parameter checks earlier
    return std::exp(log_pdf_val);
}

// CDF for SWTN (numerical integration of the pdf)
// t_adj is time already adjusted for non-decision time (t - theta)
// [[Rcpp::export]]
double pswtn(double t_adj, double alpha, double mu_drift, double sigma_drift,
                             double abs_err = 1e-6, double rel_err = 1e-6, size_t max_eval = 1000) {
    if (t_adj <= 0) return 0.0;
    if (alpha <= 0) return 1.0; // Hit boundary immediately if alpha is at or below 0
    if (sigma_drift < 0) return R_NaN;

    // Handle sigma_drift == 0 case (becomes standard Wald CDF)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return 0.0; // No positive drift
        // Use standard Wald CDF: pwald(t_adj, alpha, mu_drift)
        return pwald(t_adj, alpha, mu_drift);
    }

    pswtn_Params params_struct;
    params_struct.t_adj = t_adj;
    params_struct.alpha = alpha;
    params_struct.mu_drift = mu_drift;
    params_struct.sigma_drift = sigma_drift;

    double lower_int_bound_xi = 0.0; // For gsl_integration_qags, this is the lower limit. Upper is +Inf.

    double integral_val;
    double integral_err;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_eval); // max_eval can guide workspace size
    gsl_function F;
    F.function = &pswtn_integrand;
    F.params = &params_struct;

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status = gsl_integration_qagiu(&F, lower_int_bound_xi,
                                       abs_err, rel_err, max_eval,
                                       w, &integral_val, &integral_err);
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
      Rcpp::warning("SWTN CDF GSL integration (qagiu) did not converge; result may be inaccurate. GSL_ERROR: %s", gsl_strerror(status));
      // Fallback or error handling might be needed. For now, return the possibly inaccurate value or 0.
      // Depending on the error, integral_val might be usable or not.
      // For now, let's assume if it fails, it's 0, but this might need refinement.
      if (status == GSL_EROUND) { // Result is probably ok, just couldn't reach full accuracy
          // Potentially accept integral_val
      } else {
          // For other errors, maybe return 0 or NaN
          Rcpp::Rcout << "GSL pswtn Error: " << gsl_strerror(status) << std::endl;
          // integral_val = 0.0; // Or some other indicator of failure
      }
    }

    if (std::isnan(integral_val) || integral_val < 0.0) return 0.0;
    if (integral_val > 1.0) return 1.0;

    return integral_val;
}


// Random Number Generator for Truncated Normal N(mu, sigma) lower_bound=0
// Uses Inverse Transform Sampling
// [[Rcpp::export]]
double rtnorm(double mu, double sigma) {
	Rcpp::RNGScope scope;
    if (sigma <= 1e-10) { // Zero variance
        return (mu > 0) ? mu : 0.0; // If mu <= 0 and no variance, cannot be > 0. This needs careful thought for TN.
                                    // If sigma is 0, it's a point mass at mu. If mu > 0, sample is mu. If mu <=0, sample is undefined for lower_bound=0.
                                    // For safety, if mu<=0 and sigma is tiny, this should be handled by caller or return error.
                                    // Let's assume if sigma is tiny, rtnorm is effectively mu if mu > 0, else problematic.
                                    // For now, if mu <=0 and sigma is 0, this will struggle.
                                    // Let's return mu if mu > 0, else a very small positive if mu <=0 (hackish). A robust TN might be better.
                                    // For now, if mu <= 0 and sigma is ~0, this implies no positive drift can be sampled.
        if (mu > 0) return mu;
        else return 1e-9; // Avoid division by zero later, but this is not ideal.
    }
    double lower_bound_std = -mu / sigma; // Standardized lower bound (0-mu)/sigma

    double p_lower = R::pnorm(lower_bound_std, 0.0, 1.0, true, false);
    if (p_lower > 1.0 - 1e-7) { // Virtually no mass above lower bound
        // This means mu is very far below 0. Sample a value very close to 0.
        // This can happen if mu_drift is very negative.
        return 1e-9; // A small positive number to avoid issues in alpha/xi
    }

    double u = R::runif(0,1);
    double p_sample = u * (1.0 - p_lower) + p_lower;
    if (p_sample > 1.0 - 1e-10) p_sample = 1.0 - 1e-10; // Cap at just under 1 to avoid qnorm(1) = Inf
    if (p_sample < 1e-10) p_sample = 1e-10;// Cap at just above 0

    double x_std = R::qnorm(p_sample, 0.0, 1.0, true, false);
    return x_std * sigma + mu;
}

// Random Number Generator for Inverse Gaussian (Wald)
// Based on Michael, Schucany, Haas (1976) algorithm
// mu_ig: mean of IG (alpha/xi for us)
// lambda_ig: shape of IG (alpha^2 for us)
// [[Rcpp::export]]
double rinvgauss_rng(double mu_ig, double lambda_ig) {
	Rcpp::RNGScope scope;
    if (mu_ig <= 0 || lambda_ig <= 0) return R_PosInf; // Or handle error appropriately

    double nu = R::rnorm(0,1); // Z ~ N(0,1)
    double y = nu * nu;
    double x = mu_ig + (mu_ig * mu_ig * y) / (2.0 * lambda_ig) -
               (mu_ig / (2.0 * lambda_ig)) * std::sqrt(4.0 * mu_ig * lambda_ig * y + mu_ig * mu_ig * y * y);

    double u = R::runif(0,1);
    if (u <= mu_ig / (mu_ig + x)) {
        return x;
    } else {
        return mu_ig * mu_ig / x;
    }
}

// Random Number Generator for SWTN
// Returns one sample
// [[Rcpp::export]]
double rswtn(double alpha, double mu_drift, double sigma_drift, double theta) {
	Rcpp::RNGScope scope;
    if (alpha <= 0) return theta; // Hits boundary immediately

    // 1. Sample drift rate xi from TruncatedNormal(mu_drift, sigma_drift, lower=0)
    double xi = rtnorm(mu_drift, sigma_drift);
    if (xi <= 1e-9) { // If sampled drift is effectively zero or negative (due to numerical limits or extreme params)
        return R_PosInf; // Cannot reach boundary with this drift
    }

    // 2. Sample RT from InverseGaussian(mean = alpha/xi, shape = alpha^2)
    // Standard IG parameterization: mean `m`, shape `l`. PDF: sqrt(l/(2*pi*t^3)) * exp(-l*(t-m)^2/(2*m^2*t))
    // For Wald N(v,a,s=1) -> mean RT = a/v. This is our mu_IG.
    // The shape parameter `lambda_ig` for `rinvgauss` in R is `alpha^2` if `mean = alpha/xi`.
    // Let's verify this common IG parameterization for rinvgauss_rng:
    // mu_ig = alpha/xi
    // lambda_ig = alpha^2 (this is a common one, also called 'shape')
    double rt_sample = rinvgauss_rng(alpha / xi, alpha * alpha);

    return rt_sample + theta;
}


// --- Structs for SPV integration ---
struct SWTN_SPV_Integrand_Params_PDF {
    double t_adj;
    // k_center and A_spv are not needed here as current_k is the integration variable
    double mu_drift;
    double sigma_drift;
    // A_spv is needed for normalization (1/A_spv) if done outside integrand
};

struct SWTN_SPV_Integrand_Params_CDF {
    double t_adj;
    double mu_drift;
    double sigma_drift;
    // A_spv for normalization
};


// --- Integrands for Start Point Variability ---
// GSL style integrand for swtn_spv_logpdf (integrating over k)
double gsl_swtn_spv_integrand(double current_k, void* p) {
    SWTN_SPV_Integrand_Params_PDF* params = static_cast<SWTN_SPV_Integrand_Params_PDF*>(p);

    if (current_k <= 1e-9) { // dswtn expects alpha (threshold) > 0.
        return 0.0;
    }
    // Calculate PDF for current_k (alpha for dswtn)
    double pdf_val = dswtn(params->t_adj, current_k, params->mu_drift, params->sigma_drift);

    // Check for -Inf or very small pdf_val
    if (!std::isfinite(pdf_val) || pdf_val < -1e-300) {
        return 0.0;
    } else {
        return pdf_val;
    }
}

// GSL style integrand for swtn_spv_cdf (integrating over k)
double gsl_pswtn_spv_integrand(double current_k, void* p) {
    SWTN_SPV_Integrand_Params_CDF* params = static_cast<SWTN_SPV_Integrand_Params_CDF*>(p);

    // pswtn handles current_k <= 0 internally by returning 1.0 if drift is positive.
    // Since QAGIU integrates from a lower bound (e.g., 0), current_k will be non-negative.
    // pswtn itself ensures its alpha (current_k here) is handled appropriately.
    return pswtn(params->t_adj, current_k, params->mu_drift, params->sigma_drift);
}

/* Pretty sure these are deprecated (used the k_center form, not RDM form)
// --- Core functions with Start Point Variability (SPV) ---
// Log-PDF for SWTN with Start Point Variability (uniform over k_center +/- A_spv/2)
// k_center: central threshold, A_spv: range of start point variability
// [[Rcpp::export]]
double dswtn_spv(double t_adj, double k_center, double A_spv,
                                        double mu_drift, double sigma_drift,
                                        double abs_err = 1e-6, double rel_err = 1e-6, size_t max_eval = 1000) {
    if (t_adj <= 1e-10) return R_NegInf;
    if (A_spv < 0) return R_NaN; // Invalid A_spv (negative width)

    // If no start point variability (A_spv is tiny), use the non-SPV version
    if (A_spv < 1e-7) {
        if (k_center <= 1e-9) return R_NegInf; // Central k must be positive if no SPV
        return dswtn(t_adj, k_center, mu_drift, sigma_drift);
    }

    // Max threshold k_center + A_spv/2.0. If this is <=0, then R_NegInf, as integrand will be 0 for positive k.
    if ((k_center + A_spv/2.0) <= 1e-9) return R_NegInf;

    SWTN_SPV_Integrand_Params_PDF int_params;
    int_params.t_adj = t_adj;
    int_params.mu_drift = mu_drift;
    int_params.sigma_drift = sigma_drift;

    double lower_k_bound = k_center - A_spv / 2.0;
    double upper_k_bound = k_center + A_spv / 2.0;

    // Clip lower_k_bound if it's non-positive, as dswtn expects positive threshold.
    // The integrand gsl_dswtn_spv_integrand also has a check for current_k > 1e-9.
    if (lower_k_bound < 1e-9) lower_k_bound = 1e-9;

    // If after clipping, the integration range is invalid or zero.
    if (lower_k_bound >= upper_k_bound) return R_NegInf;

    double integral_val;
    double integral_err;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_eval);
    gsl_function F;
    F.function = &gsl_swtn_spv_integrand;
    F.params = &int_params;

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status = gsl_integration_qags(&F, lower_k_bound, upper_k_bound, abs_err, rel_err, max_eval,
                                       w, &integral_val, &integral_err);
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
      // Rcpp::Rcout << "Warning: SWTN+SPV PDF GSL integration (qags) for k failed. GSL Error: " << gsl_strerror(status) << std::endl;
    }

    if (integral_val <= 0) return R_NegInf; // Integral itself should be positive for a PDF component before normalization

    // Restore normalization by A_spv
    double final_pdf = integral_val / A_spv;
    if (final_pdf <= 1e-300) return R_NegInf; // Avoid log(0) after normalization

    return final_pdf;
}


// CDF for SWTN with Start Point Variability (uniform over k_center +/- A_spv/2)
// [[Rcpp::export]]
double swtn_spv_cdf(double t_adj, double k_center, double A_spv,
                                     double mu_drift, double sigma_drift,
                                     double abs_err = 1e-6, double rel_err = 1e-6, size_t max_eval = 1000) {
    if (t_adj <= 0) return 0.0;
    if (A_spv < 0) return R_NaN; // Invalid A_spv

    // If no start point variability (A_spv is tiny)
    if (A_spv < 1e-7) {
        // pswtn handles k_center <= 0 by returning 1.0 if drift positive, or 0 if not.
        return pswtn(t_adj, k_center, mu_drift, sigma_drift, abs_err, rel_err, max_eval);
    }

    // If max possible threshold (k_center + A_spv/2.0) is effectively zero or less:
    // The integrand gsl_pswtn_spv_integrand calls pswtn.
    // pswtn(alpha<=0) returns 1.0 if drift >0.
    // So, if the entire range [lower_k_bound, upper_k_bound] is <=0,
    // the integral of pswtn(k) over this range would be integral of ~1.0.
    // The result after normalization by A_spv would be ~1.0.
    if ((k_center + A_spv / 2.0) <= 1e-9) {
         // Check if there's any positive drift, as pswtn for alpha<=0 depends on it.
         // A quick check: if mu_drift (mean of drift dist) is positive, assume some drift can be positive.
         // Or rely on pswtn to evaluate correctly for k<=0.
         // A more robust way: if the entire integration range is non-positive,
         // call pswtn with a representative non-positive k (e.g., upper_k_bound).
         // If pswtn(non_positive_k) is 1, then the result here is 1. If 0, then 0.
         return pswtn(t_adj, k_center + A_spv / 2.0, mu_drift, sigma_drift, abs_err, rel_err, max_eval);
    }

    SWTN_SPV_Integrand_Params_CDF int_params;
    int_params.t_adj = t_adj;
    int_params.mu_drift = mu_drift;
    int_params.sigma_drift = sigma_drift;

    double lower_k_bound = k_center - A_spv / 2.0;
    double upper_k_bound = k_center + A_spv / 2.0;

    // Unlike PDF, for CDF, k can be <=0 for pswtn, which handles it by returning 1 or 0.
    // So, no strict need to clip lower_k_bound to be > 0 for the integration itself,
    // as long as gsl_pswtn_spv_integrand correctly passes k to pswtn.
    // However, if the entire range is effectively non-positive, the check above handles it.
    // If lower_k_bound >= upper_k_bound (e.g. A_spv is tiny but not caught by <1e-7, or k_center makes it so),
    // the integral will be zero, and normalized result zero, which is fine.

    double integral_val;
    double integral_err;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_eval);
    gsl_function F;
    F.function = &gsl_pswtn_spv_integrand;
    F.params = &int_params;

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status = gsl_integration_qags(&F, lower_k_bound, upper_k_bound, abs_err, rel_err, max_eval,
                                       w, &integral_val, &integral_err);
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
      // Rcpp::Rcout << "Warning: SWTN+SPV CDF GSL integration (qags) for k failed. GSL Error: " << gsl_strerror(status) << std::endl;
    }

    // Restore normalization by A_spv
    double final_cdf = integral_val / A_spv;

    if (std::isnan(final_cdf) || final_cdf < 0.0) return 0.0;
    if (final_cdf > 1.0) return 1.0;

    return final_cdf;
}
*/


// --- RDM_SWTN: Combines RDM-style SPV with SWTN drift variability ---
// These are the new top-level core functions.
// They will use dswtn and pswtn as the "inner" functions
// when drift variability is present.
// Parameters for the integrand when integrating over RDM-style SPV U(B, B+A)
struct RDM_SWTN_SPV_Integrand_Params {
    double t_adj;
    double mu_drift;
    double sigma_drift;
    // B and A define the integration range, not passed in struct here.
};

// GSL-style Integrand for PDF in rdm_swtn: dswtn(t_adj, current_actual_k, mu_drift, sigma_drift)
// current_actual_k is the integration variable (threshold), integrated from 0 to Inf.
double gsl_rdm_swtn_spv_integrand(double current_actual_k, void* p) {
    RDM_SWTN_SPV_Integrand_Params* params = static_cast<RDM_SWTN_SPV_Integrand_Params*>(p);

    if (current_actual_k <= 1e-9) { // Threshold for dswtn must be positive
        return 0.0;
    }
    double pdf_val = dswtn(params->t_adj, current_actual_k, params->mu_drift, params->sigma_drift);
    if (!std::isfinite(pdf_val) || pdf_val < 1e-300) { // exp(-700) is ~0
        return 0.0;
    } else {
        return pdf_val;
    }
}

// GSL-style Integrand for CDF in rdm_swtn: pswtn(t_adj, current_actual_k, mu_drift, sigma_drift)
// current_actual_k is the integration variable (threshold), integrated from 0 to Inf.
double gsl_rdm_pswtn_spv_integrand(double current_actual_k, void* p) {
    RDM_SWTN_SPV_Integrand_Params* params = static_cast<RDM_SWTN_SPV_Integrand_Params*>(p);

    // pswtn handles current_actual_k <= 0 appropriately.
    // Integration from 0 to Inf, so current_actual_k >= 0.
    return pswtn(params->t_adj, current_actual_k, params->mu_drift, params->sigma_drift,
                 1e-7, 1e-7, 1000); // Using tighter fixed defaults for inner pswtn integration
}


// Top-level Log-PDF for RDM_SWTN model
// Parameters B, A, mu_drift, are assumed to be already scaled by s if applicable.
// [[Rcpp::export]]
double drdmswtn(double t_adj, double B, double mu_drift, double A,
                                    double sigma_drift,
                                    double spv_abs_err = 1e-6, double spv_rel_err = 1e-6, size_t spv_max_eval = 1000) {

    if (t_adj <= 1e-10) return R_NegInf;

    bool no_A_var = (A < 1e-7);
    bool no_drift_var = (sigma_drift < 1e-10);

    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald (like dwald). Threshold is B. Drift is mu_drift.
        if (B <= 1e-9) return R_NegInf; // Threshold must be positive
        if (mu_drift <= 1e-9 && B > 0) return R_NegInf; // No positive drift to positive boundary
        return dwald(t_adj, B, mu_drift); // dwald(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
        // dWald calls digt(t_adj, B + 0.5*A, mu_drift, 0.5*A)
        // So, k_for_digt = B + 0.5*A, a_for_digt = 0.5*A
        if ((B + A) <= 1e-9) return R_NegInf; // Max threshold non-positive
        double k_digt = B + 0.5 * A;
        double a_digt = 0.5 * A;
        if (k_digt - a_digt < 0 && A > 1e-7) { // if B is negative
             // digt handles negative B if A is large enough to make parts of U(B,B+A) positive.
             // However, the k for digt (center) should ideally be positive if a_digt is not dominant.
             // For safety, if B (the minimum threshold) < 1e-9, but B+A > 1e-9,
             // digt needs to be robust. The original digt seems to handle this.
        }
         // Ensure a_digt is not zero if A was meant to be non-zero.
        if (A < 1e-7) a_digt = 0.0; // force if A is tiny

        double pdf_val = digt(t_adj, k_digt, mu_drift, a_digt);
        if (pdf_val <= 1e-300) return R_NegInf;
        return pdf_val;
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold B (like original dswtn).
        if (B <= 1e-9) return R_NegInf;
        return dswtn(t_adj, B, mu_drift, sigma_drift);
    } else {
        // Case 4: Full model - SWTN with RDM-style SPV.
        // Integrate exp(dswtn(t_adj, actual_k, mu_drift, sigma_drift))
        // where actual_k ~ U(B, B+A).
        // Note: 1e-7 has already confirmed A is not zero for this 'else' block.
        // So, no_A_var is false here.

        // If max threshold B+A is non-positive, result is -Inf
        if ((B + A) <= 1e-9) return R_NegInf;

        RDM_SWTN_SPV_Integrand_Params int_params_pdf;
        int_params_pdf.t_adj = t_adj;
        int_params_pdf.mu_drift = mu_drift;
        int_params_pdf.sigma_drift = sigma_drift;

        double lower_k_bound_pdf = B;
        double upper_k_bound_pdf = B + A;

        // Clip lower bound if it's non-positive, as dswtn expects positive threshold.
        // The integrand gsl_rdm_swtn_spv_integrand also checks current_actual_k > 1e-9.
        if (lower_k_bound_pdf < 1e-9) lower_k_bound_pdf = 1e-9;

        // If after clipping, the integration range is invalid or zero.
        if (lower_k_bound_pdf >= upper_k_bound_pdf) return R_NegInf;

        double integral_val_pdf;
        double integral_err_pdf;

        gsl_integration_workspace* w_pdf = gsl_integration_workspace_alloc(spv_max_eval);
        gsl_function F_pdf;
        F_pdf.function = &gsl_rdm_swtn_spv_integrand;
        F_pdf.params = &int_params_pdf;

        gsl_error_handler_t* old_handler_pdf = gsl_set_error_handler_off();
        int status_pdf = gsl_integration_qags(&F_pdf, lower_k_bound_pdf, upper_k_bound_pdf, spv_abs_err, spv_rel_err, spv_max_eval,
                                           w_pdf, &integral_val_pdf, &integral_err_pdf);
        gsl_set_error_handler(old_handler_pdf);
        gsl_integration_workspace_free(w_pdf);

        if (status_pdf != GSL_SUCCESS) {
            // Rcpp::Rcout << "Warning: RDM_SWTN PDF GSL integration (qags) for actual_k failed. GSL Error: " << gsl_strerror(status_pdf) << std::endl;
        }

        if (integral_val_pdf <= 0) return R_NegInf; // Integral should be positive

        // Restore normalization by A (width of the uniform distribution U(B, B+A))
        double final_pdf_val = integral_val_pdf / A;
        if (final_pdf_val <= 1e-300) return R_NegInf; // Avoid log(0)
        return final_pdf_val;
    }
}

// Top-level CDF for RDM_SWTN model
// [[Rcpp::export]]
double prdmswtn(double t_adj, double B, double mu_drift, double A,
                                 double sigma_drift,
                                 double spv_abs_err = 1e-6, double spv_rel_err = 1e-6, size_t spv_max_eval = 1000) {
    if (t_adj <= 0) return 0.0;

    bool no_A_var = (A < 1e-7);
    bool no_drift_var = (sigma_drift < 1e-10);

    if (no_A_var && no_drift_var) {
        // Case 1
        if (B <= 1e-9) return 1.0;
        if (mu_drift <= 1e-9 && B > 0) return 0.0;
        return pwald(t_adj, B, mu_drift);
    } else if (!no_A_var && no_drift_var) {
        // Case 2
        if ((B + A) <= 1e-9) return 1.0;
        double k_digt = B + 0.5 * A;
        double a_digt = 0.5 * A;
        // if (A < 1e-7) a_digt = 0.0; // A is not zero here due to !no_A_var
        return pigt(t_adj, k_digt, mu_drift, a_digt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3
        if (B <= 1e-9) return 1.0;
        return pswtn(t_adj, B, mu_drift, sigma_drift, spv_abs_err, spv_rel_err, spv_max_eval);
    } else {
        // Case 4: Full model - SWTN CDF with RDM-style SPV (uniform over U(B, B+A)).
        // A is not zero here.

        // If max threshold B+A is non-positive, effectively an immediate hit.
        if ((B + A) <= 1e-9) {
             return pswtn(t_adj, B + A, mu_drift, sigma_drift, spv_abs_err, spv_rel_err, spv_max_eval);
        }

        RDM_SWTN_SPV_Integrand_Params int_params_cdf;
        int_params_cdf.t_adj = t_adj;
        int_params_cdf.mu_drift = mu_drift;
        int_params_cdf.sigma_drift = sigma_drift;

        double lower_k_bound_cdf = B;
        double upper_k_bound_cdf = B + A;

        // For CDF, pswtn handles k <= 0. No need to clip lower_k_bound_cdf if B is negative,
        // as long as the range [B, B+A] is valid.
        // If lower_k_bound_cdf >= upper_k_bound_cdf (e.g., if A was negative, though caught earlier, or tiny),
        // integral will be zero, which is fine.

        double integral_val_cdf;
        double integral_err_cdf;

        gsl_integration_workspace* w_cdf = gsl_integration_workspace_alloc(spv_max_eval);
        gsl_function F_cdf;
        F_cdf.function = &gsl_rdm_pswtn_spv_integrand;
        F_cdf.params = &int_params_cdf;

        gsl_error_handler_t* old_handler_cdf = gsl_set_error_handler_off();
        int status_cdf = gsl_integration_qags(&F_cdf, lower_k_bound_cdf, upper_k_bound_cdf, spv_abs_err, spv_rel_err, spv_max_eval,
                                           w_cdf, &integral_val_cdf, &integral_err_cdf);
        gsl_set_error_handler(old_handler_cdf);
        gsl_integration_workspace_free(w_cdf);

        if (status_cdf != GSL_SUCCESS) {
            // Rcpp::Rcout << "Warning: RDM_SWTN CDF GSL integration (qags) for actual_k failed. GSL Error: " << gsl_strerror(status_cdf) << std::endl;
        }

        // Restore normalization by A
        double final_cdf_val = integral_val_cdf / A;

        if (std::isnan(final_cdf_val) || final_cdf_val < 0.0) return 0.0;
        if (final_cdf_val > 1.0) return 1.0;
        return final_cdf_val;
    }
}

// --- Rcpp Exported Functions ---

/* Log-density of the RDM_SWTN distribution

#' Calculates the log probability density function (PDF) for the
#' Wald distribution with trial-varying normally distributed drift rates (SWTN).
#' This distribution is sometimes referred to as an IG_TN (Inverse Gaussian
#' mixed with Truncated Normal drift rates).
#'
#' @param t Vector of reaction times.
#' @param B Vector of base threshold parameters (lower bound of SPV uniform range).
#' @param A Vector of start-point variability range parameters (width of SPV uniform range, actual_k ~ U(B, B+A)). Must be non-negative.
#' @param mu_drift Vector of mean drift rate parameters.
#' @param sigma_drift Vector of standard-deviation parameters for the drift rate distribution. Must be non-negative.
#' @param t0 Vector of non-decision time parameters.
#' @param s Vector of overall scaling parameters (typically 1.0). Scales B, A, mu_drift, but not sigma_drift.
#' @return A numeric vector of log-density values.
#' @details This function models RTs based on a Wald process where:
#'   1. The start-point (threshold `k`) varies from trial to trial, uniformly distributed over `[B/s, (B+A)/s]`.
#'   2. The drift rate (`xi`) varies from trial to trial, drawn from a normal
#'      distribution N(mu_drift/s, sigma_drift), truncated at 0.
#'   Handles vectorization of parameters. */
// [[Rcpp::export]]
NumericVector dRDM_SWTN(NumericVector t, NumericVector B, NumericVector mu_drift, NumericVector A, NumericVector t0, NumericVector s, NumericVector sigma_drift) {
    int n = t.size();
    // Vector recycling
    if (B.size() == 1) B = rep(B, n);
    if (A.size() == 1) A = rep(A, n);
    if (mu_drift.size() == 1) mu_drift = rep(mu_drift, n);
    if (sigma_drift.size() == 1) sigma_drift = rep(sigma_drift, n);
    if (t0.size() == 1) t0 = rep(t0, n);
    if (s.size() == 1) s = rep(s, n);

    NumericVector pdf(n);
    for (int i = 0; i < n; ++i) {
        double s_val = s[i];
        if (s_val <= 1e-10) { // Avoid division by zero if s is too small
            pdf[i] = R_NegInf;
            continue;
        }
        double B_scaled = B[i] / s_val;
        double A_scaled = A[i] / s_val;
        double mu_drift_scaled = mu_drift[i] / s_val;

        double t_adj = t[i] - t0[i];
        pdf[i] = drdmswtn(t_adj, B_scaled, mu_drift_scaled, A_scaled, sigma_drift[i]);
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
#' @param B Vector of base threshold parameters.
#' @param A Vector of start-point variability range parameters.
#' @param mu_drift Vector of mean drift rate parameters.
#' @param sigma_drift Vector of variance parameters for drift rate.
#' @param t0 Vector of non-decision time parameters.
#' @param s Vector of overall scaling parameters.
#' @param spv_abs_err,spv_rel_err,spv_max_eval Control integration for start-point variability.
#' @param dr_abs_err,dr_rel_err,dr_max_eval Control integration for drift-rate variability (inner integral).
#' @return A numeric vector of CDF values.
#' @details Parameters are handled similarly to `dRDM_SWTN`, including scaling by `s`. */
// [[Rcpp::export]]
NumericVector pRDM_SWTN(NumericVector t, NumericVector B, NumericVector mu_drift, NumericVector A, NumericVector t0, NumericVector s, NumericVector sigma_drift, 
                         double spv_abs_err = 1e-6, double spv_rel_err = 1e-6, int spv_max_eval = 1000
                         ) { 
                         // To add finer control for inner drift integration, params for pswtn could be exposed
                         // For now, prdmswtn passes spv_abs_err etc. to pswtn if it's Case 3 (no SPV)
                         // or uses fixed defaults if SPV is active (Case 4 integrand call to pswtn).
                         // This can be refined if needed.
    int n = t.size();
    // Vector recycling
    if (B.size() == 1) B = rep(B, n);
    if (A.size() == 1) A = rep(A, n);
    if (mu_drift.size() == 1) mu_drift = rep(mu_drift, n);
    if (sigma_drift.size() == 1) sigma_drift = rep(sigma_drift, n);
    if (t0.size() == 1) t0 = rep(t0, n);
    if (s.size() == 1) s = rep(s, n);

    NumericVector cdf_val(n);
    for (int i = 0; i < n; ++i) {
        double s_val = s[i];
        if (s_val <= 1e-10) {
            cdf_val[i] = ( (t[i] - t0[i]) > 0 ) ? 1.0 : 0.0; // Effectively infinite drift if s=0
            continue;
        }
        double B_scaled = B[i] / s_val;
        double A_scaled = A[i] / s_val;
        double mu_drift_scaled = mu_drift[i] / s_val;

        double t_adj = t[i] - t0[i];
        cdf_val[i] = prdmswtn(t_adj, B_scaled, mu_drift_scaled, A_scaled, sigma_drift[i],
                                        spv_abs_err, spv_rel_err, static_cast<size_t>(spv_max_eval));
    }
    return cdf_val;
}

/*
#' Random number generation for the RDM_SWTN model
#'
#' Generates random samples from the RDM_SWTN model.
#'
#' @param n_samples The number of samples to generate.
#' @param B Base threshold parameter.
#' @param A Start-point variability range.
#' @param mu_drift Mean drift rate.
#' @param sigma_drift Variance of drift rate.
#' @param t0 Non-decision time.
#' @param s Overall scaling parameter (default 1.0).
#' @return A numeric vector of random samples.
#' @details Scaled parameters `B/s`, `A/s`, `mu_drift/s`, `sigma_drift/(s*s)` are used.
#'   A start point `actual_k` is sampled from `U(B_scaled, B_scaled + A_scaled)`.
#'   If `sigma_drift` is zero, samples from standard RDM with SPV.
#'   If `A` is zero, samples from SWTN with fixed threshold `B_scaled`.
#'   If both are zero, samples from simple Wald. */
// [[Rcpp::export]]
NumericVector rRDM_SWTN(int n_samples, double B, double mu_drift, double A, double t0, double s = 1.0, double sigma_drift=0.0) {
    NumericVector samples(n_samples);
    Rcpp::RNGScope scope;

    if (s <= 1e-10) { // effectively infinite drift or zero threshold
        for(int i=0; i<n_samples; ++i) samples[i] = t0; // or R_PosInf if B > 0
        return samples;
    }

    double B_scaled = B / s;
    double A_scaled = A / s;
    double mu_drift_scaled = mu_drift / s;

    bool no_A_var = (A_scaled < 1e-7); // Use scaled A for this check
    bool no_drift_var = (sigma_drift < 1e-10);


    for (int i = 0; i < n_samples; ++i) {
        double current_k_for_wald; // This is the threshold for the (potentially drift-varying) Wald process

        if (no_A_var) {
            current_k_for_wald = B_scaled;
        } else {
            current_k_for_wald = R::runif(B_scaled, B_scaled + A_scaled);
        }

        if (current_k_for_wald <= 1e-9) { // Threshold is effectively non-positive
             samples[i] = t0; // Instant hit at non-decision time
             continue;
        }

        if (no_drift_var) {
            // Case 1 (A=0, sigmasq=0) or Case 2 (A>0, sigmasq=0)
            // Both reduce to sampling from a simple Wald distribution with threshold = current_k_for_wald
            // and drift = mu_drift_scaled.
            // Need a simple Wald RNG here. rinvgauss_rng(mu_ig, lambda_ig, rng)
            // mu_ig = current_k_for_wald / mu_drift_scaled
            // lambda_ig = current_k_for_wald^2
            if (mu_drift_scaled <= 1e-9) { // No positive drift
                 samples[i] = R_PosInf; // Cannot reach positive boundary
                 continue;
            }
            samples[i] = rinvgauss_rng(current_k_for_wald / mu_drift_scaled, std::pow(current_k_for_wald,2)) + t0;
        } else {
            // Case 3 (A=0, sigmasq > 0) or Case 4 (A>0, sigmasq > 0)
            // Both use rswtn with threshold = current_k_for_wald
            samples[i] = rswtn(current_k_for_wald, mu_drift_scaled, sigma_drift, t0);
        }
    }
    return samples;
}

// [[Rcpp::export]]
NumericVector drdmswtn_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4, sv=5
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = digt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4), pars(i,5));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

// [[Rcpp::export]]
NumericVector prdmswtn_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = pigt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4), pars(i,5));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

// It would be good to adapt drdm_c and prdm_c or create new versions
// that can select between the original RDM and this SWTN variant.
// This would require adding a model type parameter and dispatching.
// For now, these are standalone functions for the SWTN.

// --- RDM Race Model using SWTN ---
// Calculates log-likelihood for multiple trials for an N-accumulator race model
// where each accumulator is a SWTN process.
// Parameters can vary per trial and per accumulator.
/*
#' Log-likelihood for an N-accumulator RDM_SWTN race model
#'
#' Calculates the log-likelihood for multiple trials of an N-accumulator race model,
#' where each accumulator follows the RDM_SWTN model (RDM-style start-point
#' variability combined with SWTN-style trial-varying drift rates).
#'
#' @param rts Vector of observed reaction times for each trial.
#' @param choices Vector of observed choices for each trial, as integers from 1 to N_acc.
#' @param params_B Matrix (N_trials x N_acc) of base threshold (B) parameters.
#' @param params_A Matrix (N_trials x N_acc) of start-point variability range (A) parameters.
#' @param params_mu_drift Matrix (N_trials x N_acc) of mean drift rate parameters.
#' @param params_sigma_drift Matrix (N_trials x N_acc) of drift rate variance parameters.
#' @param params_t0 Matrix (N_trials x N_acc) of non-decision time (t0) parameters.
#' @param params_s Matrix (N_trials x N_acc) of overall scaling (s) parameters.
#' @param min_log_lik Minimum log-likelihood value to return for invalid inputs or numerical issues.
#' @param cdf_abs_err Desired absolute error for the internal CDF calculations.
#' @param cdf_rel_err Desired relative error for the internal CDF calculations.
#' @param cdf_max_eval Maximum evaluations for internal CDF calculations.
#' @return A numeric vector of log-likelihood values for each trial.
#' @details This function computes the likelihood of observing each RT and choice pair,
#'   given the parameters for N SWTN accumulators racing independently. Parameters can
#'   vary per trial and per accumulator. 
// [[Rcpp::export]]
NumericVector loglik_RDM_SWTN_race(
    NumericVector rts,
    IntegerVector choices,
    NumericMatrix params_B,
    NumericMatrix params_A,
    NumericMatrix params_mu_drift,
    NumericMatrix params_sigma_drift,
    NumericMatrix params_t0,
    NumericMatrix params_s,
    double min_log_lik = -1e10,
    double spv_abs_err = 1e-6,
    double spv_rel_err = 1e-6,
    int    spv_max_eval = 1000
) { 
    // Note: Inner drift variability integration errors are currently fixed in rdm_pswtn_spv_integrand's call to pswtn
    // or inherited from spv_abs_err in rdm_swtn/cdf Case 3. This could be further parameterized if needed.
    int n_trials = rts.size();
    if (n_trials == 0) return NumericVector(0);

    // Dimension checks
    if (choices.size() != n_trials || params_B.nrow() != n_trials ||
        params_A.nrow() != n_trials || params_mu_drift.nrow() != n_trials ||
        params_sigma_drift.nrow() != n_trials || params_t0.nrow() != n_trials ||
        params_s.nrow() != n_trials) {
        Rcpp::stop("Input vector/matrix dimensions do not match number of trials.");
    }

    int n_acc = params_B.ncol();
    if (n_acc == 0) {
        if (n_trials > 0) Rcpp::stop("Number of accumulators (ncol of param matrices) is zero.");
        return NumericVector(n_trials);
    }
    if (params_A.ncol() != n_acc || params_mu_drift.ncol() != n_acc ||
        params_sigma_drift.ncol() != n_acc || params_t0.ncol() != n_acc ||
        params_s.ncol() != n_acc) {
        Rcpp::stop("Parameter matrices must have the same number of columns (accumulators).");
    }

    NumericVector trial_log_likelihoods(n_trials);

    for (int i = 0; i < n_trials; ++i) {
        double rt = rts[i];
        int winning_acc_idx = choices[i] - 1;

        if (winning_acc_idx < 0 || winning_acc_idx >= n_acc) {
            trial_log_likelihoods[i] = min_log_lik;
            continue;
        }

        // Winner parameters
        double s_winner = params_s(i, winning_acc_idx);
        if (s_winner <= 1e-10) { trial_log_likelihoods[i] = min_log_lik; continue; }
        double B_winner_scaled = params_B(i, winning_acc_idx) / s_winner;
        double A_winner_scaled = params_A(i, winning_acc_idx) / s_winner;
        double mu_drift_winner_scaled = params_mu_drift(i, winning_acc_idx) / s_winner;
        double t0_winner = params_t0(i, winning_acc_idx);
        double t_adj_winner = rt - t0_winner;

        if (t_adj_winner <= 1e-10) {
            trial_log_likelihoods[i] = min_log_lik;
            continue;
        }

        double current_log_lik = drdmswtn(
            t_adj_winner, B_winner_scaled, A_winner_scaled,
            mu_drift_winner_scaled, sigma_drift_winner,
            spv_abs_err, spv_rel_err, static_cast<size_t>(spv_max_eval)
        );

        if (!std::isfinite(current_log_lik) || current_log_lik < min_log_lik -100) {
             trial_log_likelihoods[i] = min_log_lik;
             continue;
        }

        // Loser parameters and survivor functions
        for (int j = 0; j < n_acc; ++j) {
            if (j == winning_acc_idx) continue;

            double s_loser = params_s(i, j);
            if (s_loser <= 1e-10) { current_log_lik = min_log_lik; break; } // If s_loser is 0, implies problem
            double B_loser_scaled = params_B(i, j) / s_loser;
            double A_loser_scaled = params_A(i, j) / s_loser;
            double mu_drift_loser_scaled = params_mu_drift(i, j) / s_loser;
            double sigma_drift_sq_loser_scaled = params_sigma_drift_sq(i, j) / (s_loser * s_loser);
            double t0_loser = params_t0(i, j);
            double t_adj_loser = rt - t0_loser;

            if (t_adj_loser <= 1e-10) {
                // CDF is 0, Survivor is 1, log(Survivor) is 0. No change.
                continue;
            }

            double cdf_loser = prdmswtn(
                t_adj_loser, B_loser_scaled, A_loser_scaled,
                mu_drift_loser_scaled, sigma_drift_sq_loser_scaled,
                spv_abs_err, spv_rel_err, static_cast<size_t>(spv_max_eval)
            );

            double survivor_loser = 1.0 - cdf_loser;

            if (survivor_loser <= 1e-10) {
                current_log_lik = min_log_lik;
                break;
            }
            current_log_lik += std::log(survivor_loser);
             if (!std::isfinite(current_log_lik) || current_log_lik < min_log_lik -100) {
                 current_log_lik = min_log_lik;
                 break;
            }
        }
        trial_log_likelihoods[i] = std::max(current_log_lik, min_log_lik);
    }
    return trial_log_likelihoods;
}*/
#endif
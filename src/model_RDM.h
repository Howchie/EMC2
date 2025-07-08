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
#include "bivnorm.h"
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

// Standard Wald PDF: f(t | drift_rate_xi, boundary_alpha) equivalent to pigt0 but with threshold/drift directly parameterized
// t must be > 0
// [[Rcpp::export]]
double dwald_classic(double t, double boundary_alpha, double drift_rate_xi) {
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
double pwald_classic(double t, double boundary_alpha, double drift_rate_xi) {
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
double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return pigt0(t, k, l);
  }
  double sqt = std::sqrt(t);
  double lgt = std::log(t);
  double cdf;

  if (l < threshold){
    double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
    double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1; // ZH query is -k-a correct here or is it supposed to k-a like everywhere else?

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
	//Rcout << "--- Debug: digt called. drift: " << l<<" threshold: "<< k<<"spv: "<< a<< std::endl;
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
NumericVector dRDM_c(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    double t_adj = t[i] = t[i] - t0[i];
    if (t_adj <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t_adj, B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pRDM_c(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    double t_adj = t[i] = t[i] - t0[i];
    if (t_adj <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t_adj, B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return cdf;
}


// PDF for SWTN -- expressed in natural scale terms and tested to ensure it gives the same output as Steingrover 2021 log_pdf JAGS implementation
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation
// [[Rcpp::export]]
double dswtn(double t_adj, double alpha, double mu_drift, double sigma_drift) {
    if (t_adj <= 1e-10) return R_NegInf; // log(0)
    if (alpha <= 1e-10) return R_NegInf; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative
    // Handle sigma_drift == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return dwald_classic(t_adj, alpha, mu_drift);
    }

    double d = mu_drift;
    double v = std::pow(sigma_drift,2); // parameterize the RDM with std but formula uses variance
	double prob_d_gt_0 = R::pnorm(d / std::sqrt(v), 0.0, 1.0, true, false);
	if (std::isinf(prob_d_gt_0) || prob_d_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero
	
	double term1 = alpha;
	double term2 = std::sqrt((2 * M_PI * std::pow(t_adj,3))*((t_adj*v)+1));
	double term3 = 1 / prob_d_gt_0; // This is 1/P(xi>0)
    double term4 = std::exp( - (std::pow(d * t_adj - alpha, 2)) / ((2 * t_adj) * (t_adj * v + 1)) );
    double term5 = R::pnorm( (alpha * v + d) /(std::sqrt(t_adj * std::pow(v, 2) + v) ), 0.0, 1.0, true, false);
	double pdf_val = (term1/term2) * term3 * term4 * term5;

    if (std::isnan(pdf_val)) return R_NegInf; // Should be caught by specific parameter checks earlier
    return pdf_val;
}

// CDF for SWTN -- derived from the above PDF using a substitutuion with bivariate normals
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation
// [[Rcpp::export]]
double pswtn(double t_adj, double alpha, double mu_drift, double sigma_drift) {
    if (t_adj <= 1e-10) return 0.0; // log(0)
    if (alpha <= 1e-10) return 1.0; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative
    // Handle sigma_drift == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return pwald_classic(t_adj, alpha, mu_drift);
    }

    double v = std::pow(sigma_drift,2); // parameterize the RDM with std but formula uses variance
	double prob_d_gt_0 = R::pnorm(mu_drift / sigma_drift, 0.0, 1.0, true, false);
	if (std::isinf(prob_d_gt_0) || prob_d_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero
	
	double denom = std::sqrt(t_adj*(1+t_adj*v));
	double rho = (sigma_drift*t_adj)/denom;
	// first Φ2
	double h1 =  (mu_drift*t_adj - alpha)/denom;
	double k1 =  mu_drift / sigma_drift;
	/*NumericVector upper1 = NumericVector::create(h1, k1);
    NumericMatrix corr1(2, 2);
    corr1(0, 0) = 1.0; corr1(0, 1) = rho;
    corr1(1, 0) = rho; corr1(1, 1) = 1.0;
    NumericVector term1v = pmvnorm_cpp(upper1, corr1);
	double term1 = term1v[0];*/
	//double term1 = pbivnorm_fast(h1,k1,rho);
	double term1;
    //if (std::fabs(rho) < 0.97) {
    //    term1 = norm_cdf_2d_vfast(h1, k1, rho);
    //} else {
        term1 = norm_cdf_2d_fast(h1, k1, rho);
    //}
	// second Φ2 (reflected drift)
	double mu_p = mu_drift + 2*alpha*v;
	double h2 = (-mu_p*t_adj - alpha)/denom;
	double k2 = mu_p/sigma_drift;
	/*NumericVector upper2 = NumericVector::create(h2, k2);
    NumericMatrix corr2(2, 2);
    corr2(0, 0) = 1.0; corr2(0, 1) = -rho;
    corr2(1, 0) = -rho; corr2(1, 1) = 1.0;
    NumericVector term2v = pmvnorm_cpp(upper2, corr2);
    double term2 = term2v[0];*/
	//double term2 = pbivnorm_fast(h2,k2,-rho);
	double term2;
    //if (std::fabs(rho) < 0.97) {
    //    term2 = norm_cdf_2d_vfast(h2, k2, -rho);
    //} else {
        term2 = norm_cdf_2d_fast(h2, k2, -rho);
    //}
	double cdf_val  = (term1 + std::exp(2*alpha*mu_drift + 2*std::pow(alpha,2)*v)*term2) / prob_d_gt_0;
	
	if (std::isnan(cdf_val) || cdf_val < 0.0) return 0.0;
    if (cdf_val > 1.0) return 1.0;
	return cdf_val;
}

// Top-level PDF for RDM_SWTN model
// Parameters B, A, mu_drift, sv are assumed to be already scaled by s if applicable. This integrates the dswtn function across B->B+A using a 20-point Gauss-Legendre approximation.
// [[Rcpp::export]]
double drdmswtn(double t_adj, double B, double mu_drift, double A,
                                    double sigma_drift,int n_gauss_nodes  = 20) {

    if (t_adj <= 1e-10) return 0.0;
	bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
    bool no_drift_var = (sigma_drift < 1e-7);
	
    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald (like dwald_classic). Threshold is B. Drift is mu_drift.
        if (mu_drift <= 1e-10 || B <= 1e-10) return R_NegInf; // No positive drift to positive boundary
        return dwald_classic(t_adj, B, mu_drift); // dwald_classic(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
        // dWald calls digt(t_adj, B + 0.5*A, mu_drift, 0.5*A) i.e. k-center +- half-width of A
        // So, k_for_digt = B + 0.5*A, a_for_digt = 0.5*A
        if ((B + A) <= 1e-9) return R_NegInf; // Max threshold non-positive
        double k_digt = B + 0.5 * A;
        double a_digt = 0.5 * A;
        return digt(t_adj, k_digt, mu_drift, a_digt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold B (like original dswtn). Use Steingrover et al (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability.
        if (B <= 1e-9) return R_NegInf;
        return dswtn(t_adj, B, mu_drift, sigma_drift);
    } else {
        // Case 4: Full model - SWTN with RDM-style SPV.
        // Integrate dswtn(t_adj, actual_k, mu_drift, sigma_drift)
        // where actual_k ~ U(B, B+A).
        // Note: 1e-7 has already confirmed A is not zero for this 'else' block.
        // So, no_A_var is false here.

        // If max threshold B+A is non-positive, result is -Inf. Shouldn't be possible.
        if ((B + A) <= 1e-9) return R_NegInf;

        //Get Gauss-Legendre 20-point nodes/weights from statmod (static so it only calls once)
		const Rcpp::NumericVector& nodes   = gl["nodes"];
		const Rcpp::NumericVector& weights = gl["weights"];
	
		// Map nodes to [B, B+A]
		Rcpp::NumericVector k_nodes = B + 0.5 * A * (nodes + 1.0);
	
		// 3. Evaluate CDF at each (t, k) pair and integrate
	
		double integral = 0.0;

		for (int j = 0; j < n_gauss_nodes; ++j) {
			double pdf_val = dswtn(t_adj, k_nodes[j], mu_drift, sigma_drift);
			integral      += weights[j] * pdf_val;
		}
		// Scale for Gauss-Legendre over [B, B+A] and divide by A (uniform pdf)
		// This reduces to integral * 0.5 but it's handy to see it writtent out properly
		double out = (integral * (A / 2.0)) / A;

		// Normalization by A (width of the uniform distribution U(B, B+A))
        if (std::isnan(out) || out < 0) return R_NegInf; // Integral should be positive;
		if (out > 1.0) return 1.0;
		return out;
    }
}

// Top-level CDF for RDM_SWTN model
// Parameters B, A, mu_drift, sv are assumed to be already scaled by s if applicable. This integrates the pswtn function across B->B+A using a 20-point Gauss-Legendre approximation.
// [[Rcpp::export]]
double prdmswtn(double t_adj, double B, double mu_drift, double A,
                                    double sigma_drift,int n_gauss_nodes  = 20)
{
    if (t_adj <= 1e-10) return 0.0;
	bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
    bool no_drift_var = (sigma_drift < 1e-7);

    if (no_A_var && no_drift_var) {
        // Case 1 standard Wald (neither variability parameter)
        if (B <= 1e-7) return 1.0;
        if (mu_drift <= 1e-7 && B > 0) return 0.0;
        return pwald_classic(t_adj, B, mu_drift);
    } else if (!no_A_var && no_drift_var) {
        // Case 2 if (A < 1e-7) a_digt = 0.0; // A is not zero here due to !no_A_var
        if ((B + A) <= 1e-7) return 1.0;
        double k_pigt = B + 0.5 * A;
        double a_pigt = 0.5 * A;
        return pigt(t_adj, k_pigt, mu_drift, a_pigt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold B (like original dswtn). Starts with Steingrover et al's (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability - but they derived only the pdf with a fixed mu_drift so here we're integrating that numerically with mu_drift =0->Inf
        if (B <= 1e-9 && mu_drift>1e-9) return 1.0;
        return pswtn(t_adj, B, mu_drift, sigma_drift);
    } else {

		//Get Gauss-Legendre 20-point nodes/weights from statmod (static so it only calls once)
		const Rcpp::NumericVector& nodes   = gl["nodes"];
		const Rcpp::NumericVector& weights = gl["weights"];
	
		// Map nodes to [B, B+A]
		Rcpp::NumericVector k_nodes = B + 0.5 * A * (nodes + 1.0);
	
		// 3. Evaluate CDF at each (t, k) pair and integrate
	
		double integral = 0.0;

		for (int j = 0; j < n_gauss_nodes; ++j) {
			double cdf_val = pswtn(t_adj, k_nodes[j], mu_drift, sigma_drift);
			integral      += weights[j] * cdf_val;
		}
		// Scale for Gauss-Legendre over [B, B+A] and divide by A (uniform pdf)
		// This reduces to integral * 0.5 but it's handy to see it writtent out properly
		double out = (integral * (A / 2.0)) / A;
				// Defensive clipping
		// Normalization by A (width of the uniform distribution U(B, B+A))
        if (std::isnan(out) || out < 0) return 0.0;
		if (out > 1.0) return 1.0;
		return out;
    }
    
}
/*
#' Calculates the probability density function (PDF) for the
#' Wald distribution with trial-varying normally distributed drift rates (SWTN).
#' This distribution is sometimes referred to as an IG_TN (Inverse Gaussian
#' mixed with Truncated Normal drift rates).
#'
#' @param t Vector of reaction times.
#' @param B Vector of base threshold parameters (lower bound of SPV uniform range).
#' @param A Vector of start-point variability range parameters (width of SPV uniform range, actual_k ~ U(B, B+A)). Must be non-negative.
#' @param v Vector of mean drift rate parameters.
#' @param sv Vector of standard-deviation parameters for the drift rate distribution. Must be non-negative.
#' @param t0 Vector of non-decision time parameters.
#' @return A numeric vector of log-density values.
#' @details This function models RTs based on a Wald process where:
#'   1. The start-point (threshold `k`) varies from trial to trial, uniformly distributed over `[B/s, (B+A)/s]`.
#'   2. The drift rate (`xi`) varies from trial to trial, drawn from a normal
#'      distribution N(v/s, sigma_drift), truncated at 0.
#'   Handles vectorization of parameters. */
// [[Rcpp::export]]
NumericVector dSWTNspv(NumericVector t, NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector sv) {
    int n = t.size();
    NumericVector pdf(n);
    for (int i = 0; i < n; ++i) {
		t[i] = t[i] - t0[i];
		if (t[i] <= 0){
			pdf[i] = 0.;
		} else { 
			pdf[i] = drdmswtn(t[i], B[i], v[i], A[i], sv[i]);
		}
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
#' @param v Vector of mean drift rate parameters.
#' @param sv Vector of variance parameters for drift rate.
#' @param t0 Vector of non-decision time parameters.
#' @param spv_abs_err,spv_rel_err,spv_max_eval Control integration for start-point variability.
#' @param dr_abs_err,dr_rel_err,dr_max_eval Control integration for drift-rate variability (inner integral).
#' @return A numeric vector of CDF values.
#' @details Parameters are handled similarly to `dRDM_SWTN`, including scaling by `s`. */
// [[Rcpp::export]]
NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector sv) { 
    int n = t.size();
    NumericVector cdf(n);
    for (int i = 0; i < n; i++){
		t[i] = t[i] - t0[i];
		if (t[i] <= 0){
		cdf[i] = 0.;
		} else {
			cdf[i] = prdmswtn(t[i], B[i], v[i], A[i], sv[i]);
		}
	}
    return cdf;
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
		//double sv = pars(i,0)*pars(i,5); // convert coefficient of variation to standard deviation
        out[k] = drdmswtn(rts[i] - pars(i,3), pars(i,1)/pars(i,4), pars(i,0)/pars(i,4), pars(i,2)/pars(i,4), pars(i,5));
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
  //v = 0, B = 1, A = 2, t0 = 3, s = 4, sv=5
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
		//double sv = pars(i,0)*pars(i,5); // convert coefficient of variation to standard deviation
        out[k] = prdmswtn(rts[i] - pars(i,3), pars(i,1)/pars(i,4), pars(i,0)/pars(i,4), pars(i,2)/pars(i,4), pars(i,5));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}


// Shifted Wald with Truncated-Normal Drift Variability code adapted from https://github.com/HelenSteingroever/jags-wald
// PDF for SWTN -- this one uses the log-pdf form from Steingrover's JAGS code (exponentiated to give probability)
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation
/*
// [[Rcpp::export]]
double dswtn_log(double t_adj, double alpha, double mu_drift, double sigma_drift) {
    if (t_adj <= 1e-10) return R_NegInf; // log(0)
    if (alpha <= 1e-10) return R_NegInf; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative
    // Handle sigma_drift_sq == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return dwald_classic(t_adj, alpha, mu_drift);
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

    if (std::isnan(log_pdf_val)) return R_NegInf; // Should be caught by specific parameter checks earlier
    return std::exp(log_pdf_val);
}

Random generation functions are here but standard seems to be R implementation
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
*/

/*
// START OF SWTN IMPLEMENTATION
// Helper struct for CDF integration parameters
struct pswtn_Params {
    double t_adj;        // Time (adjusted for non-decision time)
    double alpha;        // Threshold
    double mu_drift;     // Mean of the drift rate distribution
    double sigma_drift; // SD of the drift rate distribution
};

// Integrand for the SWTN CDF: dwald_classic * truncated_normal_a_pdf (from library)
// This is the GSL-style function that will be integrated over xi.
double pswtn_integrand(double current_xi, void* p) {
    pswtn_Params* params = static_cast<pswtn_Params*>(p);

    // current_xi comes from GSL, lower bound is handled by qagiu's 'a' parameter.
    // However, explicit check for safety or specific logic can remain if needed.
    if (current_xi < 1e-10) { // Should ideally not be hit if gsl_integration_qagiu is called with lower_bound_xi = 0
        return 0.0;
    }

    double wc = pwald_classic(params->t_adj, params->alpha, current_xi);
    double tn = truncated_normal_a_pdf(current_xi, params->mu_drift, params->sigma_drift, 0.0); // lower_bound = 0 for drift rate


	if (std::isnan(wc) || std::isnan(tn) || !std::isfinite(wc) || !std::isfinite(tn) || wc < 0 || tn < 0) {
        Rcpp::Rcout << "pswtn_integrand warning: NaN/Inf/Negative produced. xi=" << current_xi 
                   << ", t_adj=" << params->t_adj << ", alpha=" << params->alpha
                   << ", mu_d=" << params->mu_drift << ", sig_d=" << params->sigma_drift
                   << ", pwald_cdf=" << wc << ", tnorm_pdf=" << tn << std::endl;
        return 0.0; 
    }
	
    return wc * tn;
}

// CDF for SWTN (numerical integration of the pdf)
// t_adj is time already adjusted for non-decision time (t - theta)
// [[Rcpp::export]]
double pswtn_numeric_integral(double t_adj, double alpha, double mu_drift, double sigma_drift,
                             double abs_err = 1e-8, double rel_err = 1e-8, size_t max_eval = 10000) {
    if (t_adj <= 1e-10) return 0.0; // Never hits boundary
    if (alpha <= 1e-10) return 1.0; // Hit boundary immediately if alpha is at or below 0
    if (sigma_drift < 0) return R_NaN; // Variance can't be negative (should be bounded in EMC anyway)
    // Handle sigma_drift == 0 case (becomes standard Wald CDF)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return 0.0; // No positive drift
        // Use standard Wald CDF: pwald_classic(t_adj, alpha, mu_drift)
        return pwald_classic(t_adj, alpha, mu_drift);
    }

    pswtn_Params params_struct;
    params_struct.t_adj = t_adj;
    params_struct.alpha = alpha;
    params_struct.mu_drift = mu_drift;
    params_struct.sigma_drift = sigma_drift;

    double lower_int_bound_xi = 0.0; // For gsl_integration_qags, this is the lower limit. Upper is +Inf.

    double integral_val = 0.0;
    double integral_err = 0.0;

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
		Rcpp::Rcerr << "pswtn GSL_ERROR: " << gsl_strerror(status) 
                  << " (code " << status << ")\n"
                  << "  Inputs: t_adj=" << t_adj << ", alpha=" << alpha 
                  << ", mu_drift=" << mu_drift << ", sigma_drift=" << sigma_drift << "\n"
                  << "  Requested abs_err=" << abs_err << ", rel_err=" << rel_err 
                  << ", max_eval=" << max_eval << "\n"
                  << "  Returned integral_val=" << integral_val << ", integral_err=" << integral_err 
                  << std::endl;
      // If divergence or other critical error, returning NaN might be better to propagate error
      if (status == GSL_EDIVERGE || status == GSL_EMAXITER || status == GSL_ESING || status == GSL_EBADFUNC) {
          return R_NaN; // Indicate failure more strongly
      }
    }

    if (std::isnan(integral_val) || integral_val < 0.0) return 0.0;
    if (integral_val > 1.0) return 1.0;
	//Rcout << "--- Debug: pswtn called. adj_t= "<<t_adj<<" drift: " << mu_drift<<" threshold: "<< alpha<<" sv: "<< sigma_drift<< "integral: "<<integral_val<< std::endl;
    return integral_val;
}

// --- RDMSWTN: Combines RDM-style SPV with SWTN drift variability ---
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

// --- Integrands for Start Point Variability ---
// GSL style integrand for swtn_spv_logpdf (integrating over k, equivalent to integrating out spv)
double gsl_dswtn_spv_integrand(double current_k, void* p) {
    RDM_SWTN_SPV_Integrand_Params* params = static_cast<RDM_SWTN_SPV_Integrand_Params*>(p);

    if (current_k <= 1e-10) { // dswtn expects alpha (threshold) > 0.
        return 0.0;
    }
    // Calculate PDF for current_k (alpha for dswtn)
    double pdf_val = dswtn(params->t_adj, current_k, params->mu_drift, params->sigma_drift);

    // Check for -Inf or very small pdf_val
	return (std::isfinite(pdf_val) && pdf_val >= 0) ? pdf_val : 0.0;
}

// GSL style integrand for swtn_spv_cdf (integrating over k, equivalent to integrating out spv)
double gsl_pswtn_spv_integrand(double current_k, void* p) {
    RDM_SWTN_SPV_Integrand_Params* params = static_cast<RDM_SWTN_SPV_Integrand_Params*>(p);

    // pswtn handles current_k <= 0 internally by returning 1.0 if drift is positive.
    // Since QAGIU integrates from a lower bound (e.g., 0), current_k will be non-negative.
    // pswtn itself ensures its alpha (current_k here) is handled appropriately.
    return pswtn(params->t_adj, current_k, params->mu_drift, params->sigma_drift); // Using tighter fixed defaults for inner pswtn integration
}

// Top-level PDF for RDM_SWTN model
// Parameters B, A, mu_drift, sv are assumed to be already scaled by s if applicable.
// [[Rcpp::export]]
double drdmswtn_numeric_integral(double t_adj, double B, double mu_drift, double A,
                                    double sigma_drift,
                                    double spv_abs_err = 1e-8, double spv_rel_err = 1e-8, size_t spv_max_eval = 10000) {

    if (t_adj <= 1e-10) return 0.0;
	bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
    bool no_drift_var = (sigma_drift < 1e-7);
	
    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald (like dwald_classic). Threshold is B. Drift is mu_drift.
        if (mu_drift <= 1e-10 || B <= 1e-10) return R_NegInf; // No positive drift to positive boundary
        return dwald_classic(t_adj, B, mu_drift); // dwald_classic(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
        // dWald calls digt(t_adj, B + 0.5*A, mu_drift, 0.5*A) i.e. k-center +- half-width of A
        // So, k_for_digt = B + 0.5*A, a_for_digt = 0.5*A
        if ((B + A) <= 1e-9) return R_NegInf; // Max threshold non-positive
        double k_digt = B + 0.5 * A;
        double a_digt = 0.5 * A;
        return digt(t_adj, k_digt, mu_drift, a_digt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold B (like original dswtn). Use Steingrover et al (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability.
        if (B <= 1e-9) return R_NegInf;
        return dswtn(t_adj, B, mu_drift, sigma_drift);
    } else {
        // Case 4: Full model - SWTN with RDM-style SPV.
        // Integrate dswtn(t_adj, actual_k, mu_drift, sigma_drift)
        // where actual_k ~ U(B, B+A).
        // Note: 1e-7 has already confirmed A is not zero for this 'else' block.
        // So, no_A_var is false here.

        // If max threshold B+A is non-positive, result is -Inf. Shouldn't be possible.
        if ((B + A) <= 1e-9) return R_NegInf;

        RDM_SWTN_SPV_Integrand_Params int_params_pdf;
        int_params_pdf.t_adj = t_adj;
        int_params_pdf.mu_drift = mu_drift;
        int_params_pdf.sigma_drift = sigma_drift;

        double lower_k_bound_pdf = B;
        double upper_k_bound_pdf = B + A;

        double integral_val_pdf = 0.0;
        double integral_err_pdf = 0.0;

        gsl_integration_workspace* w_pdf = gsl_integration_workspace_alloc(spv_max_eval);
        gsl_function F_pdf;
        F_pdf.function = &gsl_dswtn_spv_integrand;
        F_pdf.params = &int_params_pdf;

        gsl_error_handler_t* old_handler_pdf = gsl_set_error_handler_off();
        int status_pdf = gsl_integration_qags(&F_pdf, lower_k_bound_pdf, upper_k_bound_pdf, spv_abs_err, spv_rel_err, spv_max_eval,
                                           w_pdf, &integral_val_pdf, &integral_err_pdf);
        gsl_set_error_handler(old_handler_pdf);
        gsl_integration_workspace_free(w_pdf);

        if (status_pdf != GSL_SUCCESS) {
            Rcpp::Rcerr << "drdmswtn (SPV) GSL_ERROR: " << gsl_strerror(status_pdf) << " (code " << status_pdf << ")\n"
                        << "  Inputs: t_adj=" << t_adj << ", B=" << B << ", A=" << A 
                        << ", mu_drift=" << mu_drift << ", sigma_drift=" << sigma_drift << std::endl;
            return 0.0; // 0? Or NaN to indicate error
        }

        if (integral_val_pdf <= 0) return R_NegInf; // Integral should be positive

        // Normalization by A (width of the uniform distribution U(B, B+A))
        if (std::isnan(integral_val_pdf)) return 0.0;
        return integral_val_pdf / A;
    }
}


// OLD NUMERIC INTEGRATION STUFF
// Top-level CDF for RDM_SWTN model
// [[Rcpp::export]]
double prdmswtn_numeric_integral(double t_adj, double B, double mu_drift, double A,
                                 double sigma_drift,
                                 double spv_abs_err = 1e-8, double spv_rel_err = 1e-8, size_t spv_max_eval = 10000) {
    if (t_adj <= 0) return 0.0;
	//Rcout << "--- Debug: prdmswtn called. drift: " << mu_drift<<" threshold: "<< B<<" spv: "<< A<<" sv: "<< sigma_drift<< std::endl;
    bool no_A_var = (A < 1e-7);
    bool no_drift_var = (sigma_drift < 1e-7);

    if (no_A_var && no_drift_var) {
        // Case 1 standard Wald (neither variability parameter)
        if (B <= 1e-7) return 1.0;
        if (mu_drift <= 1e-7 && B > 0) return 0.0;
        return pwald_classic(t_adj, B, mu_drift);
    } else if (!no_A_var && no_drift_var) {
        // Case 2 if (A < 1e-7) a_digt = 0.0; // A is not zero here due to !no_A_var
        if ((B + A) <= 1e-7) return 1.0;
        double k_pigt = B + 0.5 * A;
        double a_pigt = 0.5 * A;
        return pigt(t_adj, k_pigt, mu_drift, a_pigt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold B (like original dswtn). Starts with Steingrover et al's (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability - but they derived only the pdf with a fixed mu_drift so here we're integrating that numerically with mu_drift =0->Inf
        if (B <= 1e-9 && mu_drift>1e-9) return 1.0;
        return pswtn(t_adj, B, mu_drift, sigma_drift);
    } else {
        // Case 4: Full model - SWTN CDF with RDM-style SPV (uniform over U(B, B+A)).
        // If max threshold B+A is non-positive, effectively an immediate hit. (check this..)
        if ((B + A) <= 1e-9) {
             return 1.0;
        }

        RDM_SWTN_SPV_Integrand_Params int_params_cdf;
        int_params_cdf.t_adj = t_adj;
        int_params_cdf.mu_drift = mu_drift;
        int_params_cdf.sigma_drift = sigma_drift;

        double lower_k_bound_cdf = B;
        double upper_k_bound_cdf = B + A;

        double integral_val_cdf = 0.0;
        double integral_err_cdf = 0.0;

        gsl_integration_workspace* w_cdf = gsl_integration_workspace_alloc(spv_max_eval);
        gsl_function F_cdf;
        F_cdf.function = &gsl_pswtn_spv_integrand;
        F_cdf.params = &int_params_cdf;

        gsl_error_handler_t* old_handler_cdf = gsl_set_error_handler_off();
        int status_cdf = gsl_integration_qags(&F_cdf, lower_k_bound_cdf, upper_k_bound_cdf, spv_abs_err, spv_rel_err, spv_max_eval,
                                           w_cdf, &integral_val_cdf, &integral_err_cdf);
        gsl_set_error_handler(old_handler_cdf);
        gsl_integration_workspace_free(w_cdf);

        if (status_cdf != GSL_SUCCESS) {
             Rcpp::Rcerr << "prdmswtn (SPV) GSL_ERROR: " << gsl_strerror(status_cdf) << " (code " << status_cdf << ")\n"
                        << "  Inputs: t_adj=" << t_adj << ", B=" << B << ", A=" << A 
                        << ", mu_drift=" << mu_drift << ", sigma_drift=" << sigma_drift << std::endl;
            return 0.0; // Indicate failure
        }

		// Normalization by A (width of the uniform distribution U(B, B+A))
        if (std::isnan(integral_val_cdf) || integral_val_cdf < 0) return 0.0;
        return integral_val_cdf / A;
    }
}
*/

#endif
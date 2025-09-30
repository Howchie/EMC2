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
#include "truncated_normal.h"
#include <string>

using namespace Rcpp;

const double L_PI = std::log(M_PI);

// RDM
// [[Rcpp::export]]
double pigt0(double t, double k = 1.0, double l = 1.0, bool log_out=false){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1.0 + t/mu), 0.0, 1.0, true, false);
  double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1.0 - t/mu), 0.0, 1.0, true, false);
  double cdf_val = std::exp(2.0 * lambda / mu + std::log(p1)) + p2;
  if (cdf_val < 0.0 || std::isnan(cdf_val)) {
  	  return log_out ? R_NegInf : 0;
  }
  if (cdf_val > 1.0) {
	  return log_out ? 0 : 1;
  }
  return log_out ? std::log(cdf_val) : cdf_val;
}

// [[Rcpp::export]]
double digt0(double t, double k = 1.0, double l = 1.0, bool log_out=false){
  //if (t <= 0.) {
  //  return 0.;
  //}
  double lambda = k * k;
  double e;
  if (l == 0.0) {
    e = -.5 * lambda / t;
  } else {
    double mu = k / l;
    e = - (lambda / (2.0 * t)) * ((t * t) / (mu * mu) - 2.0 * t / mu + 1.0);
  }
  double pdf_val = std::exp(e + .5 * std::log(lambda) - .5 * std::log(2.0 * t * t * t * M_PI));
	if (pdf_val < 0.0 || std::isnan(pdf_val)) {
		return log_out ? R_NegInf : 0;
	}
  return log_out ? std::log(pdf_val) : pdf_val;
}

// Standard Wald PDF: f(t | drift_rate_xi, boundary_alpha) equivalent to pigt0 but with threshold/drift directly parameterized
// t must be > 0
// [[Rcpp::export]]
double dwald_classic(double t, double boundary_alpha, double drift_rate_xi, bool log_out=false) {
    if (t <= 0) return 0.0;
    if (boundary_alpha <= 0) return (drift_rate_xi > 0) ? 1.0 : 0.0; // if boundary is 0 or less, instant absorption if drift is positive
    if (drift_rate_xi <= 1e-10) { // Effectively zero or negative drift away from boundary
        return 0.0; // Cannot reach positive boundary alpha
    }

    double pdf_val = boundary_alpha*std::exp(-std::pow((boundary_alpha-drift_rate_xi*t),2)/(2*t))/std::sqrt(2*M_PI*std::pow(t,3));
	if (pdf_val < 0.0 || std::isnan(pdf_val)) {
		return log_out ? R_NegInf : 0;
	}
	return log_out ? std::log(pdf_val) : pdf_val;
}

// Standard Wald CDF: F(t | drift_rate_xi, boundary_alpha) equivalent to digt0 but with threshold/drift directly parameterized
// t must be > 0
// [[Rcpp::export]]
double pwald_classic(double t, double boundary_alpha, double drift_rate_xi, bool log_out=false) {
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
	if (cdf_val < 0.0 || std::isnan(cdf_val)) {
		return log_out ? R_NegInf : 0;
	}
	if (cdf_val > 1.0) {
		return log_out ? 0 : 1;
	}
    return log_out ? std::log(cdf_val) : cdf_val;
}

// [[Rcpp::export]]
double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10, bool log_out=false){
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return pigt0(t, k, l);
  }
  double sqt = std::sqrt(t);
  double lgt = std::log(t);
  double cdf_val;

  if (l < threshold){
    double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
    double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1; // ZH query is -k-a correct here or is it supposed to k-a like everywhere else?

    double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

    cdf_val = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
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

    cdf_val = .5 * (t4 + t2 + t1) / a;
  }
	if (cdf_val < 0.0 || std::isnan(cdf_val)) {
		return log_out ? R_NegInf : 0;
	}
	if (cdf_val > 1.0) {
		return log_out ? 0 : 1;
	}
  return log_out ? std::log(cdf_val) : cdf_val;
}

// [[Rcpp::export]]
double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10, bool log_out=false){
	//Rcout << "--- Debug: digt called. drift: " << l<<" threshold: "<< k<<"spv: "<< a<< std::endl;
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return digt0(t, k, l);
  }
  double pdf_val;
  if (l < threshold){
    double term = std::exp(- (k - a) * (k - a) / (2. * t)) - std::exp(- (k + a) * (k + a) / (2. * t));
    pdf_val = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(term) - M_LN2 - std::log(a));
  } else {
    double sqt = std::sqrt(t);

    double t1a = - std::pow(a - k + t * l, 2) / (2. * t);
    double t1b = - std::pow(a + k - t * l, 2) / (2. * t);
    double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

    double t2a = 2. * R::pnorm((- k + a) / sqt + sqt * l, 0., 1., true, false) - 1.;
    double t2b = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);

    pdf_val = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
  }
  if (pdf_val < 0.0 || std::isnan(pdf_val)) {
    return log_out ? R_NegInf : 0;
  }
  return log_out ? std::log(pdf_val) : pdf_val;
}

// [[Rcpp::export]]
NumericVector drdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = digt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4), log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }

  return(out);
}

// [[Rcpp::export]]
NumericVector prdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = pigt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4), log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }

  return(out);
}


// [[Rcpp::export]]
NumericVector dRDM_c(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0, bool log_out=false){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    double t_adj = t[i] = t[i] - t0[i];
    if (t_adj <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t_adj, B[i] + .5 * A[i], v[i], .5 * A[i], log_out);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pRDM_c(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0, bool log_out=false){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    double t_adj = t[i] = t[i] - t0[i];
    if (t_adj <= 0){
      cdf[i] = 0.0;
    } else {
      cdf[i] = pigt(t_adj, B[i] + .5 * A[i], v[i], .5 * A[i], log_out);
    }
  }
  return cdf;
}

// Shifted Wald with Truncated-Normal Drift Variability code adapted from https://github.com/HelenSteingroever/jags-wald
// PDF for SWTN -- this one uses the log-pdf form from Steingrover's JAGS code (exponentiated to give probability)
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation

// [[Rcpp::export]]
double dswtn(double t_adj, double alpha, double mu_drift, double sigma_drift, bool log_out=false) {
    if (t_adj <= 1e-10) return R_NegInf; // log(0)
    if (alpha <= 1e-10) return R_NegInf; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative
    // Handle sigma_drift_sq == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return dwald_classic(t_adj, alpha, mu_drift,log_out);
    }

    //double lambda = 1.0; // Matches SWTN.cc
    double v = std::pow(sigma_drift,2); // parameterize the RDM with std but formula uses variance
    // Pre-compute for numerical stability
    const double tv  = t_adj * v;            // tsigma^2
    const double log1p_tv = std::log1p(tv);  // log(1 + tsigma^2)
    // Log of the normalization constant for the drift rate xi ~ TN(mu_drift, v, lower=0)
    // log(1 / P(xi > 0)) = -log(P(xi > 0))
    // P(xi > 0) = pnorm(mu_drift/sqrt(v), 0, 1, true, false) if mu_drift > 0 is assumed for TN's mode.
    // Or more generally, P(xi > 0) for xi ~ N(mu_drift,v) is 1 - pnorm(0, mu_drift, sqrt(v), true, false)
    // = pnorm(0, mu_drift, sqrt(v), false, false)
    // = pnorm( (0-mu_drift)/sqrt(v), 0, 1, false, false) = pnorm(-mu_drift/sqrt(v), 0, 1, false, false)
    // = pnorm(mu_drift/sqrt(v), 0, 1, true, false)
    double log_prob_xi_gt_0 = R::pnorm(mu_drift / sigma_drift, 0.0, 1.0, true, true); // log(Phi(mu_drift/sqrt(v)))
    if (std::isinf(log_prob_xi_gt_0) && log_prob_xi_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero

    double term1 = std::log(alpha);
    double term2 = 0.5 * (- std::log(2.0) - L_PI - 3.0 * std::log(t_adj) - log1p_tv);
    double term3 = -log_prob_xi_gt_0; // This is -log(P(xi>0)) = log (1/P(xi>0))
    double term4 = -(std::pow(mu_drift * t_adj - alpha, 2.0)) / (2.0 * t_adj * (tv + 1.0));
    //double term5 = R::pnorm((alpha * v + mu_drift) / std::sqrt(lambda * t_adj * std::pow(v,2) + v), 0.0, 1.0, true, true); // ZH edit: Factor out v^2 (which is sigma^4) here for numerical stability
	double z     = (alpha * sigma_drift + (mu_drift/sigma_drift)) / std::sqrt(tv + 1.0);
    double term5 = R::pnorm(z, 0.0, 1.0,true, true);
	double log_pdf_val = term1 + term2 + term3 + term4 + term5;

    if (std::isnan(log_pdf_val)) return R_NegInf; // Should be caught by specific parameter checks earlier
    return log_out ? log_pdf_val : std::exp(log_pdf_val);
}

// CDF for SWTN -- derived from the above PDF using a substitutuion with bivariate normals
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation
// [[Rcpp::export]]
double pswtn(double t_adj, double alpha, double mu_drift, double sigma_drift, bool log_out=false) {
    if (t_adj <= 1e-10) return 0.0; // log(0)
    if (alpha <= 1e-10) return 1.0; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative
    // Handle sigma_drift == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return pwald_classic(t_adj, alpha, mu_drift, log_out);
    }

    double sigma_2 = std::pow(sigma_drift,2); 
	double alpha_2 = std::pow(alpha,2); 
	double prob_d_gt_0 = R::pnorm(mu_drift / sigma_drift, 0.0, 1.0, true, false);
	if (std::isinf(prob_d_gt_0) || prob_d_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero
	double log_prob_d_gt_0 = R::pnorm(mu_drift / sigma_drift, 0.0, 1.0, true, true); // log.p = TRUE
    if (std::isinf(log_prob_d_gt_0)) return R_NegInf;
	double denom = std::sqrt(t_adj*(1+t_adj*sigma_2));
	double rho = (sigma_drift*t_adj)/denom;
	// first phi2
	double h1 =  (mu_drift*t_adj - alpha)/denom;
	double k1 =  mu_drift / sigma_drift;
	double term1;
    if (std::fabs(rho) > 0.97 || std::fabs(h1)>8.0 || std::fabs(k1)>8.0) {
        term1 = norm_cdf_2d(h1, k1, rho); // Genz when in the tails
    } else {
        term1 = norm_cdf_2d_fast(h1, k1, rho); // Faster Drezner/West algorithm otherwise
    }
	double log_A = (term1 > 1e-300) ? std::log(term1) : R_NegInf; // Avoid log(0)
	// second phi2 (reflected drift)
	double mu_p = mu_drift + 2*alpha*sigma_2;
	double h2 = (-mu_p*t_adj - alpha)/denom;
	double k2 = mu_p/sigma_drift;

	double term2;
    if (std::fabs(rho) > 0.97 || std::fabs(h2)>8.0 || std::fabs(k2)>8.0) {
        term2 = norm_cdf_2d(h2, k2, -rho); // Genz when in the tails
    } else {
        term2 = norm_cdf_2d_fast(h2, k2, -rho); // Faster Drezner/West algorithm otherwise
    }
	double log_term2 = (term2 > 1e-300) ? std::log(term2) : R_NegInf; // Avoid log(0)
	double log_exp_term = 2 * alpha * mu_drift + 2 * alpha_2 * sigma_2;
	double log_B = log_exp_term + log_term2;
    // --- Combine using log-sum-exp ---
    double log_numerator = log_sum_exp(log_A, log_B);
    
    // Final log CDF
    double log_cdf_val = log_numerator - log_prob_d_gt_0;
	//double cdf_val  = (term1 + std::exp(2*alpha*mu_drift + 2*alpha_2*sigma_2)*term2) / prob_d_gt_0;
	double cdf_val = std::exp(log_cdf_val);
	if (cdf_val < 0.0 || std::isnan(cdf_val)) {
		return log_out ? R_NegInf : 0;
	}
	if (cdf_val > 1.0) {
		return log_out ? 0 : 1;
	}
	return log_out ? std::log(cdf_val) : cdf_val;
}

// Top-level PDF for RDM_SWTN model
// Parameters b, A, mu_drift, sv are assumed to be already scaled by s if applicable. This integrates the dswtn function across b->b+A using a 20-point Gauss-Legendre approximation.
// [[Rcpp::export]]
double drdmswtn(double t_adj, double b, double mu_drift, double A,
                                    double sigma_drift,int n_gauss_nodes  = 20, bool log_out=false) {

    if (t_adj <= 1e-10) return 0.0;
	if (b <= 1e-7) return R_NegInf; // Max threshold non-positive
	bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
    bool no_drift_var = (sigma_drift < 1e-7);
	
    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald (like dwald_classic). Threshold is b. Drift is mu_drift.
        if (mu_drift <= 1e-10 || b <= 1e-10) return R_NegInf; // No positive drift to positive boundary
        return dwald_classic(t_adj, b, mu_drift,log_out); // dwald_classic(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
        // dWald calls digt(t_adj, B + 0.5*A, mu_drift, 0.5*A) i.e. k-center +- half-width of A where B is the DIFFERENCE between A and threshold
		// We've use b (i.e., B+A) here, so we back-transform and use the digt function.
        double k_digt = b - 0.5 * A;
        double a_digt = 0.5 * A;
        return digt(t_adj, k_digt, mu_drift, a_digt,log_out);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold b (like original dswtn). Use Steingrover et al (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability.
        return dswtn(t_adj, b, mu_drift, sigma_drift,log_out);
    } else {
        // Case 4: Full model - SWTN with RDM-style SPV.
        // Integrate dswtn(t_adj, actual_k, mu_drift, sigma_drift)
        // where actual_k ~ U(b-A, b).
        // Note: 1e-7 has already confirmed A is not zero for this 'else' block.
        // So, no_A_var is false here.
        // If max threshold b is non-positive, result is -Inf. Shouldn't be possible.
        //Get Gauss-Legendre 20-point gl_nodes/gl_weights from statmod (static so it only calls once)
		const Rcpp::NumericVector& gl_nodes   = gl["nodes"];
		const Rcpp::NumericVector& gl_weights = gl["weights"];
	
		// Map gl_nodes to [b, b+A]
		Rcpp::NumericVector k_nodes = b - 0.5 * A + 0.5*A*gl_nodes;
		// 3. Evaluate CDF at each (t, k) pair and integrate
		double integral = 0.0;
		for (int j = 0; j < n_gauss_nodes; ++j) {
			double pdf_val = dswtn(t_adj, k_nodes[j], mu_drift, sigma_drift,false); // integral cannot take log_pdf
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
                                    double sigma_drift,int n_gauss_nodes  = 20, bool log_out=false)
{
    if (t_adj <= 1e-10) return 0.0;
	if ((b) <= 1e-7) return 1.0;
	bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
    bool no_drift_var = (sigma_drift < 1e-7);

    if (no_A_var && no_drift_var) {
        // Case 1 standard Wald (neither variability parameter)
        if (mu_drift <= 1e-7 && b > 0) return 0.0;
        return pwald_classic(t_adj, b, mu_drift,log_out);
    } else if (!no_A_var && no_drift_var) {
        // Case 2 if (A < 1e-7) a_digt = 0.0; // A is not zero here due to !no_A_var
		// pWald calls pigt(t_adj, B + 0.5*A, mu_drift, 0.5*A) i.e. k-center +- half-width of A where B is the DIFFERENCE between A and threshold
		// We've use b (i.e., B+A) here, so we back-transform and use the pigt function.
        double k_pigt = b - 0.5 * A;
        double a_pigt = 0.5 * A;
        return pigt(t_adj, k_pigt, mu_drift, a_pigt,log_out);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold b (like original dswtn). Starts with Steingrover et al's (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability - but they derived only the pdf with a fixed mu_drift so here we're integrating that numerically with mu_drift =0->Inf
        if (b <= 1e-9 && mu_drift>1e-9) return 1.0;
        return pswtn(t_adj, b, mu_drift, sigma_drift,log_out);
    } else {

		//Get Gauss-Legendre 20-point gl_nodes/gl_weights from statmod (static so it only calls once)
		const Rcpp::NumericVector& gl_nodes   = gl["nodes"];
		const Rcpp::NumericVector& gl_weights = gl["gl_weights"];
	
		// Map gl_nodes to [b, b+A]
		Rcpp::NumericVector k_nodes = b - 0.5 * A + 0.5*A*gl_nodes;
		// 3. Evaluate CDF at each (t, k) pair and integrate
		double integral = 0.0;
		for (int j = 0; j < n_gauss_nodes; ++j) {
			double cdf_val = pswtn(t_adj, k_nodes[j], mu_drift, sigma_drift,false); // integral cannot take log_cdf
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
#'      distribution N(v/s, sigma_drift), truncated at 0.
#'   Handles vectorization of parameters. */
// [[Rcpp::export]]
NumericVector dSWTNspv(NumericVector t, NumericVector v, NumericVector b, NumericVector A, NumericVector t0, NumericVector sv, bool log_out=false) {
    int n = t.size();
    NumericVector pdf(n);
    for (int i = 0; i < n; ++i) {
		t[i] = t[i] - t0[i];
		if (t[i] <= 0){
			pdf[i] = 0.;
		} else { 
			pdf[i] = drdmswtn(t[i], b[i], v[i], A[i], sv[i],log_out);
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
NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector b, NumericVector A, NumericVector t0, NumericVector sv, bool log_out=false) { 
    int n = t.size();
    NumericVector cdf(n);
    for (int i = 0; i < n; i++){
		t[i] = t[i] - t0[i];
		if (t[i] <= 0){
		cdf[i] = 0.;
		} else {
			cdf[i] = prdmswtn(t[i], b[i], v[i], A[i], sv[i],log_out);
		}
	}
    return cdf;
}

// [[Rcpp::export]]
NumericVector drdmswtn_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  //v = 0, b = 1, zA = 2, t0 = 3, s = 4, cv=5
	NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
		double sv = (pars(i,0)/pars(i,4))*pars(i,5); // convert coefficient of variation to standard deviation
        double A = (pars(i,1)/pars(i,4)) * pars(i,2); // convert zA to A
		out[k] = drdmswtn(rts[i] - pars(i,3), pars(i,1)/pars(i,4), pars(i,0)/pars(i,4), A, sv,log_out);
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
  //v = 0, b = 1, zA = 2, t0 = 3, s = 4, cv=5
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
		double sv = (pars(i,0)/pars(i,4))*pars(i,5); // convert coefficient of variation to standard deviation
        double A = (pars(i,1)/pars(i,4)) * pars(i,2); // convert zA to A
		out[k] = prdmswtn(rts[i] - pars(i,3), pars(i,1)/pars(i,4), pars(i,0)/pars(i,4), A, sv,log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }

  return(out);
}


// Numerical Integration SWTN Functions (DEPRECATED)
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
// Parameters for the integrand when integrating over RDM-style SPV U(b, b+A)
struct RDM_SWTN_SPV_Integrand_Params {
    double t_adj;
    double mu_drift;
    double sigma_drift;
    // b and A define the integration range, not passed in struct here.
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
// Parameters b, A, mu_drift, sv are assumed to be already scaled by s if applicable.
// [[Rcpp::export]]
double drdmswtn_numeric_integral(double t_adj, double b, double mu_drift, double A,
                                    double sigma_drift,
                                    double spv_abs_err = 1e-8, double spv_rel_err = 1e-8, size_t spv_max_eval = 10000) {

    if (t_adj <= 1e-10) return 0.0;
	bool no_A_var = (A < 1e-7); // setting them quite low so the reduction logic only triggers if the user has genuinely fixed the value to zero (the lower bound in EMC2 is ~1e-4 so this should never be triggered during sampling)
    bool no_drift_var = (sigma_drift < 1e-7);
	
    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald (like dwald_classic). Threshold is b. Drift is mu_drift.
        if (mu_drift <= 1e-10 || b <= 1e-10) return R_NegInf; // No positive drift to positive boundary
        return dwald_classic(t_adj, b, mu_drift); // dwald_classic(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
        // dWald calls digt(t_adj, b + 0.5*A, mu_drift, 0.5*A) i.e. k-center +- half-width of A
        // So, k_for_digt = b + 0.5*A, a_for_digt = 0.5*A
        if ((b + A) <= 1e-9) return R_NegInf; // Max threshold non-positive
        double k_digt = b + 0.5 * A;
        double a_digt = 0.5 * A;
        return digt(t_adj, k_digt, mu_drift, a_digt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold b (like original dswtn). Use Steingrover et al (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability.
        if (b <= 1e-9) return R_NegInf;
        return dswtn(t_adj, b, mu_drift, sigma_drift);
    } else {
        // Case 4: Full model - SWTN with RDM-style SPV.
        // Integrate dswtn(t_adj, actual_k, mu_drift, sigma_drift)
        // where actual_k ~ U(b, b+A).
        // Note: 1e-7 has already confirmed A is not zero for this 'else' block.
        // So, no_A_var is false here.

        // If max threshold b+A is non-positive, result is -Inf. Shouldn't be possible.
        if ((b + A) <= 1e-9) return R_NegInf;

        RDM_SWTN_SPV_Integrand_Params int_params_pdf;
        int_params_pdf.t_adj = t_adj;
        int_params_pdf.mu_drift = mu_drift;
        int_params_pdf.sigma_drift = sigma_drift;

        double lower_k_bound_pdf = b;
        double upper_k_bound_pdf = b + A;

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
                        << "  Inputs: t_adj=" << t_adj << ", b=" << b << ", A=" << A 
                        << ", mu_drift=" << mu_drift << ", sigma_drift=" << sigma_drift << std::endl;
            return 0.0; // 0? Or NaN to indicate error
        }

        if (integral_val_pdf <= 0) return R_NegInf; // Integral should be positive

        // Normalization by A (width of the uniform distribution U(b, b+A))
        if (std::isnan(integral_val_pdf)) return 0.0;
        return integral_val_pdf / A;
    }
}


// OLD NUMERIC INTEGRATION STUFF
// Top-level CDF for RDM_SWTN model
// [[Rcpp::export]]
double prdmswtn_numeric_integral(double t_adj, double b, double mu_drift, double A,
                                 double sigma_drift,
                                 double spv_abs_err = 1e-8, double spv_rel_err = 1e-8, size_t spv_max_eval = 10000) {
    if (t_adj <= 0) return 0.0;
	//Rcout << "--- Debug: prdmswtn called. drift: " << mu_drift<<" threshold: "<< b<<" spv: "<< A<<" sv: "<< sigma_drift<< std::endl;
    bool no_A_var = (A < 1e-7);
    bool no_drift_var = (sigma_drift < 1e-7);

    if (no_A_var && no_drift_var) {
        // Case 1 standard Wald (neither variability parameter)
        if (b <= 1e-7) return 1.0;
        if (mu_drift <= 1e-7 && b > 0) return 0.0;
        return pwald_classic(t_adj, b, mu_drift);
    } else if (!no_A_var && no_drift_var) {
        // Case 2 if (A < 1e-7) a_digt = 0.0; // A is not zero here due to !no_A_var
        if ((b + A) <= 1e-7) return 1.0;
        double k_pigt = b + 0.5 * A;
        double a_pigt = 0.5 * A;
        return pigt(t_adj, k_pigt, mu_drift, a_pigt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: SWTN with fixed threshold b (like original dswtn). Starts with Steingrover et al's (2021) derivation of the Shifted-Wald with Truncated Normal drift-variability - but they derived only the pdf with a fixed mu_drift so here we're integrating that numerically with mu_drift =0->Inf
        if (b <= 1e-9 && mu_drift>1e-9) return 1.0;
        return pswtn(t_adj, b, mu_drift, sigma_drift);
    } else {
        // Case 4: Full model - SWTN CDF with RDM-style SPV (uniform over U(b, b+A)).
        // If max threshold b+A is non-positive, effectively an immediate hit. (check this..)
        if ((b + A) <= 1e-9) {
             return 1.0;
        }

        RDM_SWTN_SPV_Integrand_Params int_params_cdf;
        int_params_cdf.t_adj = t_adj;
        int_params_cdf.mu_drift = mu_drift;
        int_params_cdf.sigma_drift = sigma_drift;

        double lower_k_bound_cdf = b;
        double upper_k_bound_cdf = b + A;

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
                        << "  Inputs: t_adj=" << t_adj << ", b=" << b << ", A=" << A 
                        << ", mu_drift=" << mu_drift << ", sigma_drift=" << sigma_drift << std::endl;
            return 0.0; // Indicate failure
        }

		// Normalization by A (width of the uniform distribution U(b, b+A))
        if (std::isnan(integral_val_cdf) || integral_val_cdf < 0) return 0.0;
        return integral_val_cdf / A;
    }
}


/* DEPRECATED for numerical stability purposes
// PDF for SWTN -- expressed in natural scale terms and tested to ensure it gives the same output as Steingrover 2021 log_pdf JAGS implementation
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift is drift standard deviation
// [[Rcpp::export]]
double dswtn(double t_adj, double alpha, double mu_drift, double sigma_drift, bool log_out=false) {
    if (t_adj <= 1e-10) return R_NegInf; // log(0)
    if (alpha <= 1e-10) return R_NegInf; // No boundary to hit, or ill-defined
    if (sigma_drift < 0) return R_NaN; // standard deviation cannot be negative
    // Handle sigma_drift == 0 case (becomes standard Wald)
    if (sigma_drift <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        return dwald_classic(t_adj, alpha, mu_drift, log_out);
    }

    double sigma_2 = std::pow(sigma_drift,2); // parameterize the RDM with std but formula uses variance
	double prob_d_gt_0 = R::pnorm(mu_drift / sigma_drift, 0.0, 1.0, true, false);
	if (std::isinf(prob_d_gt_0) || prob_d_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero
	
	double term1 = alpha;
	double term2 = std::sqrt((2 * M_PI * std::pow(t_adj,3))*((t_adj*sigma_2)+1));
	double term3 = 1 / prob_d_gt_0; // This is 1/P(xi>0)
    double term4 = std::exp( - (std::pow(mu_drift * t_adj - alpha, 2)) / ((2 * t_adj) * (t_adj * sigma_2 + 1)) );
    double term5 = R::pnorm( (alpha * sigma_2 + mu_drift) /(std::sqrt(t_adj * std::pow(sigma_2, 2) + sigma_2) ), 0.0, 1.0, true, false);
	double pdf_val = (term1/term2) * term3 * term4 * term5;

    if (pdf_val < 0.0 || std::isnan(pdf_val)) {
		return log_out ? R_NegInf : 0;
	}
    return log_out ? std::log(pdf_val) : pdf_val;
}
*/

#endif

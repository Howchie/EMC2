#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)

// RDM
double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. + t/mu), 0., 1., true, false);
  double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. - t/mu), 0., 1., true, false);

  return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

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

#endif

// START OF DSWTN IMPLEMENTATION

// Helper struct for CDF integration parameters
struct DSWTN_CDF_Params {
    double t_adj;        // Time (adjusted for non-decision time)
    double alpha;        // Threshold
    double mu_drift;     // Mean of the drift rate distribution
    double sigma_drift_sq; // Variance of the drift rate distribution
};

// Standard Wald CDF: F(t | drift_rate_xi, boundary_alpha)
// t must be > 0
inline double wald_cdf_classic(double t, double drift_rate_xi, double boundary_alpha) {
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

// PDF of a normal distribution N(mu, sigma_sq) truncated at lower_bound (typically 0)
// xi is the value at which to evaluate the PDF
inline double truncated_normal_pdf(double xi, double mu, double sigma_sq, double lower_bound) {
    if (xi < lower_bound) return 0.0;
    if (sigma_sq <= 1e-10) { // Effectively zero variance
        if (xi >= lower_bound && std::fabs(xi - mu) < 1e-7) return R_PosInf; // Spike at mu if mu is valid
        return 0.0;
    }
    double sigma = std::sqrt(sigma_sq);

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

// Integrand for the DSWTN CDF: wald_cdf_classic * truncated_normal_pdf
// This is the function that hcubature will integrate over xi.
int dswtn_cdf_integrand(unsigned dim, const double* xi_val, void* p, unsigned fdim, double* retval) {
    if (dim != 1 || fdim != 1) {
      // Should not happen for our 1D integral
      return 1; // Error
    }
    DSWTN_CDF_Params* params = static_cast<DSWTN_CDF_Params*>(p);
    double current_xi = xi_val[0];

    if (current_xi < 0) { // Drift rate must be positive
        retval[0] = 0.0;
        return 0;
    }

    double wc = wald_cdf_classic(params->t_adj, current_xi, params->alpha);
    double tn = truncated_normal_pdf(current_xi, params->mu_drift, params->sigma_drift_sq, 0.0); // lower_bound = 0 for drift rate

    retval[0] = wc * tn;
    return 0; // Success
}


// Log-PDF for DSWTN
// t_adj is time already adjusted for non-decision time (t - theta)
// alpha is threshold, mu_drift is mean drift, sigma_drift_sq is drift variance
inline double dswtn_logpdf_core(double t_adj, double alpha, double mu_drift, double sigma_drift_sq) {
    if (t_adj <= 1e-10) return R_NegInf; // log(0)
    if (alpha <= 1e-10) return R_NegInf; // No boundary to hit, or ill-defined
    if (sigma_drift_sq < 0) return R_NaN; // variance cannot be negative

    // Handle sigma_drift_sq == 0 case (becomes standard Wald)
    if (sigma_drift_sq <= 1e-10) {
        if (mu_drift <= 1e-10) return R_NegInf; // No positive drift
        // Standard Wald log-PDF: log(alpha) - 0.5 * log(2*pi*t_adj^3) - (alpha - mu_drift*t_adj)^2 / (2*t_adj)
        // This is for s=1. The parameters for digt0 are k=alpha, l=mu_drift
        // log(digt0(t_adj, alpha, mu_drift))
        // digt0 uses lambda_wald = k*k = alpha*alpha and mu_wald = k/l = alpha/mu_drift
        // log( (alpha / sqrt(2*PI*t_adj^3)) * exp(-(alpha - mu_drift*t_adj)^2 / (2*mu_drift^2*t_adj)) ) is not quite right.
        // Standard Wald (Johnson, Kotz, Balakrishnan Vol 1, p. 261, eq. 16.1 for a=alpha, v=mu_drift, x0=0)
        // f(t) = (alpha / sqrt(2*pi*t^3)) * exp( -(alpha - mu_drift*t)^2 / (2*t) )
        // This is if s=1 (diffusion coefficient).
        // Let's use the existing digt0 as it's tested. Parameters for digt0 are k and l.
        // k ~ boundary, l ~ drift. So k=alpha, l=mu_drift.
        return std::log(digt0(t_adj, alpha, mu_drift));
    }

    double lambda_fixed = 1.0; // Matches DSWTN.cc
    double d = mu_drift;
    double v = sigma_drift_sq;

    // Log of the normalization constant for the drift rate xi ~ TN(d, v, lower=0)
    // log(1 / P(xi > 0)) = -log(P(xi > 0))
    // P(xi > 0) = pnorm(d/sqrt(v), 0, 1, true, false) if d > 0 is assumed for TN's mode.
    // Or more generally, P(xi > 0) for xi ~ N(d,v) is 1 - pnorm(0, d, sqrt(v), true, false)
    // = pnorm(0, d, sqrt(v), false, false)
    // = pnorm( (0-d)/sqrt(v), 0, 1, false, false) = pnorm(-d/sqrt(v), 0, 1, false, false)
    // = pnorm(d/sqrt(v), 0, 1, true, false)
    double log_prob_xi_gt_0 = R::pnorm(d / std::sqrt(v), 0.0, 1.0, true, true); // log(Phi(d/sqrt(v)))
    if (std::isinf(log_prob_xi_gt_0) && log_prob_xi_gt_0 < 0) return R_NegInf; // if P(xi > 0) is zero

    double term1 = std::log(alpha);
    double term2 = 0.5 * (std::log(lambda_fixed) - std::log(2.0) - M_LN_PI - 3.0 * std::log(t_adj) - std::log(lambda_fixed * t_adj * v + 1.0));
    double term3 = -log_prob_xi_gt_0; // This is -log(P(xi>0)) = log (1/P(xi>0))
    double term4 = -(lambda_fixed * std::pow(d * t_adj - alpha, 2.0)) / (2.0 * t_adj * (lambda_fixed * t_adj * v + 1.0));
    double term5 = R::pnorm((lambda_fixed * alpha * v + d) / std::sqrt(lambda_fixed * t_adj * v * v + v), 0.0, 1.0, true, true);

    double log_pdf_val = term1 + term2 + term3 + term4 + term5;

    if (std::isnan(log_pdf_val)) return R_NegInf; // Should be caught by specific parameter checks earlier
    return log_pdf_val;
}


// Forward declaration for hcubature/gauss_kronrod if not in a shared header for model_RDM.h
// Assuming gauss_kronrod is available from hcubature.cpp linkage
// struct one_d { double result; double err; int kdivide; }; // Already in gauss.h potentially
// void gauss_kronrod(double a, double b, one_d& out, void* pars,
//                    int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval));
// --- Actually, hcubature.cpp includes gauss.h which defines one_d.
// --- We need to include "gauss.h" or ensure hcubature function is declared.
// --- For simplicity, let's assume hcubature is callable. We might need to include its header if it has one.
// --- Or, if hcubature.cpp is compiled and linked, its functions are available.
// --- Let's add #include "hcubature.cpp" if it's header-only style or ensure linkage.
// --- For now, assume it's available. The Rcpp build process should handle linkage.

// Need to declare hcubature if its definition is in hcubature.cpp
// int hcubature(int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval),
//               void* pars, unsigned n, const double* a, const double* b,
//               size_t maxEval, double reqAbsError, double reqRelError, double* val, double* err);
// This is defined in hcubature.cpp.
// The file gauss.h contains struct one_d, xd7, wd7, gwd7.
// It seems gauss.h should be included.
#include "gauss.h" // For one_d struct and gauss_kronrod constants (if gauss_kronrod is used directly)
                   // And for hcubature's own needs if it uses these directly.
                   // The hcubature.cpp has its own gauss_kronrod.

// Need to declare the hcubature function if it's not in a header.
// It's better to have a hcubature.h. For now, I'll copy the signature.
extern "C" { // If hcubature.c was compiled as C
    int hcubature(int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval),
              void* pars, unsigned n, const double* min, const double* max,
              size_t maxEval, double reqAbsError, double reqRelError, double* val, double* err);
}


// CDF for DSWTN
// t_adj is time already adjusted for non-decision time (t - theta)
inline double dswtn_cdf_core(double t_adj, double alpha, double mu_drift, double sigma_drift_sq,
                             double abs_err = 1e-6, double rel_err = 1e-6, size_t max_eval = 1000) {
    if (t_adj <= 0) return 0.0;
    if (alpha <= 0) return 1.0; // Hit boundary immediately if alpha is at or below 0
    if (sigma_drift_sq < 0) return R_NaN;

    // Handle sigma_drift_sq == 0 case (becomes standard Wald CDF)
    if (sigma_drift_sq <= 1e-10) {
        if (mu_drift <= 1e-10) return 0.0; // No positive drift
        // Use standard Wald CDF: pigt0(t_adj, alpha, mu_drift)
        return pigt0(t_adj, alpha, mu_drift);
    }

    DSWTN_CDF_Params params;
    params.t_adj = t_adj;
    params.alpha = alpha;
    params.mu_drift = mu_drift;
    params.sigma_drift_sq = sigma_drift_sq;

    double lower_int_bound_xi = 0.0;
    // Determine a practical upper integration bound for xi
    // e.g., mu_drift + some number of std devs. Max drift rate to consider.
    double upper_int_bound_xi = params.mu_drift + 8.0 * std::sqrt(params.sigma_drift_sq);
    if (upper_int_bound_xi < 1e-5) upper_int_bound_xi = 1.0; // Ensure a positive range if mu_drift is very small or negative
    if (params.mu_drift < 0 && upper_int_bound_xi < 2.0 * std::sqrt(params.sigma_drift_sq)) { // if mu is negative, ensure range covers some positive drifts
        upper_int_bound_xi = 2.0 * std::sqrt(params.sigma_drift_sq);
    }


    double integral_val;
    double integral_err;

    // hcubature expects arrays for bounds
    double a[1] = {lower_int_bound_xi};
    double b[1] = {upper_int_bound_xi};

    // Using hcubature (which might call gauss_kronrod for 1D)
    int success = hcubature(dswtn_cdf_integrand, &params, 1, a, b,
                            max_eval, abs_err, rel_err, &integral_val, &integral_err);

    if (success != 0) {
      // Rcpp::warning("DSWTN CDF integration did not converge; result may be inaccurate.");
      // Fallback or error handling might be needed. For now, return the possibly inaccurate value.
    }

    if (std::isnan(integral_val) || integral_val < 0.0) return 0.0;
    if (integral_val > 1.0) return 1.0;

    return integral_val;
}


// Random Number Generator for Truncated Normal N(mu, sigma_sq) lower_bound=0
// Uses Inverse Transform Sampling
inline double rtnorm_rng_positive(double mu, double sigma_sq, RNG& rng) {
    if (sigma_sq <= 1e-10) { // Zero variance
        return (mu > 0) ? mu : 0.0; // If mu <= 0 and no variance, cannot be > 0. This needs careful thought for TN.
                                    // If sigma_sq is 0, it's a point mass at mu. If mu > 0, sample is mu. If mu <=0, sample is undefined for lower_bound=0.
                                    // For safety, if mu<=0 and sigma_sq is tiny, this should be handled by caller or return error.
                                    // Let's assume if sigma_sq is tiny, rtnorm is effectively mu if mu > 0, else problematic.
                                    // For now, if mu <=0 and sigma_sq is 0, this will struggle.
                                    // Let's return mu if mu > 0, else a very small positive if mu <=0 (hackish). A robust TN might be better.
                                    // For now, if mu <= 0 and sigma_sq is ~0, this implies no positive drift can be sampled.
        if (mu > 0) return mu;
        else return 1e-9; // Avoid division by zero later, but this is not ideal.
    }
    double sigma = std::sqrt(sigma_sq);
    double lower_bound_std = -mu / sigma; // Standardized lower bound (0-mu)/sigma

    double p_lower = R::pnorm(lower_bound_std, 0.0, 1.0, true, false);
    if (p_lower > 1.0 - 1e-7) { // Virtually no mass above lower bound
        // This means mu is very far below 0. Sample a value very close to 0.
        // This can happen if mu_drift is very negative.
        return 1e-9; // A small positive number to avoid issues in alpha/xi
    }

    double u = rng.uniform(); // Rcpp's RNG gives U(0,1)
    double p_sample = u * (1.0 - p_lower) + p_lower;
    if (p_sample > 1.0 - 1e-10) p_sample = 1.0 - 1e-10; // Cap at just under 1 to avoid qnorm(1) = Inf
    if (p_sample < 1e-10) p_sample = 1e-10;         // Cap at just above 0

    double x_std = R::qnorm(p_sample, 0.0, 1.0, true, false);
    return x_std * sigma + mu;
}

// Random Number Generator for Inverse Gaussian (Wald)
// Based on Michael, Schucany, Haas (1976) algorithm
// mu_ig: mean of IG (alpha/xi for us)
// lambda_ig: shape of IG (alpha^2 for us)
inline double rinvgauss_rng(double mu_ig, double lambda_ig, RNG& rng) {
    if (mu_ig <= 0 || lambda_ig <= 0) return R_PosInf; // Or handle error appropriately

    double nu = rng.norm(); // Z ~ N(0,1)
    double y = nu * nu;
    double x = mu_ig + (mu_ig * mu_ig * y) / (2.0 * lambda_ig) -
               (mu_ig / (2.0 * lambda_ig)) * std::sqrt(4.0 * mu_ig * lambda_ig * y + mu_ig * mu_ig * y * y);

    double u = rng.uniform();
    if (u <= mu_ig / (mu_ig + x)) {
        return x;
    } else {
        return mu_ig * mu_ig / x;
    }
}

// Random Number Generator for DSWTN
// Returns one sample
inline double rswtn_core(double alpha, double mu_drift, double sigma_drift_sq, double theta, RNG& rng) {
    if (alpha <= 0) return theta; // Hits boundary immediately

    // 1. Sample drift rate xi from TruncatedNormal(mu_drift, sigma_drift_sq, lower=0)
    double xi = rtnorm_rng_positive(mu_drift, sigma_drift_sq, rng);
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
    double rt_sample = rinvgauss_rng(alpha / xi, alpha * alpha, rng);

    return rt_sample + theta;
}


// --- Rcpp Exported Functions ---

#' Log-density of the DSWTN distribution
#'
#' Calculates the log probability density function (PDF) for the
#' Wald distribution with trial-varying normally distributed drift rates (DSWTN).
#' This distribution is sometimes referred to as an IG_TN (Inverse Gaussian
#' mixed with Truncated Normal drift rates).
#'
#' @param t Vector of reaction times.
#' @param B Vector of base threshold parameters (lower bound of SPV uniform range).
#' @param A Vector of start-point variability range parameters (width of SPV uniform range, actual_k ~ U(B, B+A)). Must be non-negative.
#' @param mu_drift Vector of mean drift rate parameters.
#' @param sigma_drift_sq Vector of variance parameters for the drift rate distribution. Must be non-negative.
#' @param t0 Vector of non-decision time parameters.
#' @param s Vector of overall scaling parameters (typically 1.0). Scales B, A, mu_drift, and sqrt(sigma_drift_sq).
#' @return A numeric vector of log-density values.
#' @details This function models RTs based on a Wald process where:
#'   1. The start-point (threshold `k`) varies from trial to trial, uniformly distributed over `[B/s, (B+A)/s]`.
#'   2. The drift rate (`xi`) varies from trial to trial, drawn from a normal
#'      distribution N(mu_drift/s, sigma_drift_sq/(s*s)), truncated at 0.
#'   Handles vectorization of parameters.
// [[Rcpp::export]]
NumericVector dRDM_DSWTN_log(NumericVector t, NumericVector B, NumericVector A, NumericVector mu_drift, NumericVector sigma_drift_sq, NumericVector t0, NumericVector s) {
    int n = t.size();
    // Vector recycling
    if (B.size() == 1) B = rep(B, n);
    if (A.size() == 1) A = rep(A, n);
    if (mu_drift.size() == 1) mu_drift = rep(mu_drift, n);
    if (sigma_drift_sq.size() == 1) sigma_drift_sq = rep(sigma_drift_sq, n);
    if (t0.size() == 1) t0 = rep(t0, n);
    if (s.size() == 1) s = rep(s, n);

    NumericVector log_pdf(n);
    for (int i = 0; i < n; ++i) {
        double s_val = s[i];
        if (s_val <= 1e-10) { // Avoid division by zero if s is too small
            log_pdf[i] = R_NegInf;
            continue;
        }
        double B_scaled = B[i] / s_val;
        double A_scaled = A[i] / s_val;
        double mu_drift_scaled = mu_drift[i] / s_val;
        double sigma_drift_sq_scaled = sigma_drift_sq[i] / (s_val * s_val);

        double t_adj = t[i] - t0[i];
        log_pdf[i] = rdm_dswtn_logpdf_core(t_adj, B_scaled, A_scaled, mu_drift_scaled, sigma_drift_sq_scaled);
    }
    return log_pdf;
}

#' Cumulative distribution function for the RDM_DSWTN model
#'
#' Calculates the CDF for the RDM_DSWTN model, which incorporates RDM-style
#' start-point variability and DSWTN-style trial-varying drift rates.
#'
#' @param t Vector of reaction times.
#' @param B Vector of base threshold parameters.
#' @param A Vector of start-point variability range parameters.
#' @param mu_drift Vector of mean drift rate parameters.
#' @param sigma_drift_sq Vector of variance parameters for drift rate.
#' @param t0 Vector of non-decision time parameters.
#' @param s Vector of overall scaling parameters.
#' @param spv_abs_err,spv_rel_err,spv_max_eval Control integration for start-point variability.
#' @param dr_abs_err,dr_rel_err,dr_max_eval Control integration for drift-rate variability (inner integral).
#' @return A numeric vector of CDF values.
#' @details Parameters are handled similarly to `dRDM_DSWTN_log`, including scaling by `s`.
// [[Rcpp::export]]
NumericVector pRDM_DSWTN(NumericVector t, NumericVector B, NumericVector A, NumericVector mu_drift, NumericVector sigma_drift_sq, NumericVector t0, NumericVector s,
                         double spv_abs_err = 1e-6, double spv_rel_err = 1e-6, int spv_max_eval = 1000
                         // To add finer control for inner drift integration, params for dswtn_cdf_core could be exposed
                         // For now, rdm_dswtn_cdf_core passes spv_abs_err etc. to dswtn_cdf_core if it's Case 3 (no SPV)
                         // or uses fixed defaults if SPV is active (Case 4 integrand call to dswtn_cdf_core).
                         // This can be refined if needed.
                         ) {
    int n = t.size();
    // Vector recycling
    if (B.size() == 1) B = rep(B, n);
    if (A.size() == 1) A = rep(A, n);
    if (mu_drift.size() == 1) mu_drift = rep(mu_drift, n);
    if (sigma_drift_sq.size() == 1) sigma_drift_sq = rep(sigma_drift_sq, n);
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
        double sigma_drift_sq_scaled = sigma_drift_sq[i] / (s_val * s_val);

        double t_adj = t[i] - t0[i];
        cdf_val[i] = rdm_dswtn_cdf_core(t_adj, B_scaled, A_scaled, mu_drift_scaled, sigma_drift_sq_scaled,
                                        spv_abs_err, spv_rel_err, static_cast<size_t>(spv_max_eval));
    }
    return cdf_val;
}

#' Random number generation for the RDM_DSWTN model
#'
#' Generates random samples from the RDM_DSWTN model.
#'
#' @param n_samples The number of samples to generate.
#' @param B Base threshold parameter.
#' @param A Start-point variability range.
#' @param mu_drift Mean drift rate.
#' @param sigma_drift_sq Variance of drift rate.
#' @param t0 Non-decision time.
#' @param s Overall scaling parameter (default 1.0).
#' @return A numeric vector of random samples.
#' @details Scaled parameters `B/s`, `A/s`, `mu_drift/s`, `sigma_drift_sq/(s*s)` are used.
#'   A start point `actual_k` is sampled from `U(B_scaled, B_scaled + A_scaled)`.
#'   If `sigma_drift_sq` is zero, samples from standard RDM with SPV.
#'   If `A` is zero, samples from DSWTN with fixed threshold `B_scaled`.
#'   If both are zero, samples from simple Wald.
// [[Rcpp::export]]
NumericVector rRDM_DSWTN(int n_samples, double B, double A, double mu_drift, double sigma_drift_sq, double t0, double s = 1.0) {
    NumericVector samples(n_samples);
    RNG rng;

    if (s <= 1e-10) { // effectively infinite drift or zero threshold
        for(int i=0; i<n_samples; ++i) samples[i] = t0; // or R_PosInf if B > 0
        return samples;
    }

    double B_scaled = B / s;
    double A_scaled = A / s;
    double mu_drift_scaled = mu_drift / s;
    double sigma_drift_sq_scaled = sigma_drift_sq / (s * s);

    const double A_is_zero_threshold = 1e-7;
    const double sigmasq_is_zero_threshold = 1e-10;

    bool no_A_var = (A_scaled < A_is_zero_threshold); // Use scaled A for this check
    bool no_drift_var = (sigma_drift_sq_scaled < sigmasq_is_zero_threshold);


    for (int i = 0; i < n_samples; ++i) {
        double current_k_for_wald; // This is the threshold for the (potentially drift-varying) Wald process

        if (no_A_var) {
            current_k_for_wald = B_scaled;
        } else {
            current_k_for_wald = rng.uniform(B_scaled, B_scaled + A_scaled);
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
            samples[i] = rinvgauss_rng(current_k_for_wald / mu_drift_scaled, current_k_for_wald * current_k_for_wald, rng) + t0;
        } else {
            // Case 3 (A=0, sigmasq > 0) or Case 4 (A>0, sigmasq > 0)
            // Both use rswtn_core with threshold = current_k_for_wald
            samples[i] = rswtn_core(current_k_for_wald, mu_drift_scaled, sigma_drift_sq_scaled, t0, rng);
        }
    }
    return samples;
}


// It would be good to adapt drdm_c and prdm_c or create new versions
// that can select between the original RDM and this DSWTN variant.
// This would require adding a model type parameter and dispatching.
// For now, these are standalone functions for the DSWTN.


// --- RDM Race Model using DSWTN ---

// Calculates log-likelihood for multiple trials for an N-accumulator race model
// where each accumulator is a DSWTN process.
// Parameters can vary per trial and per accumulator.

#' Log-likelihood for an N-accumulator RDM_DSWTN race model
#'
#' Calculates the log-likelihood for multiple trials of an N-accumulator race model,
#' where each accumulator follows the RDM_DSWTN model (RDM-style start-point
#' variability combined with DSWTN-style trial-varying drift rates).
#'
#' @param rts Vector of observed reaction times for each trial.
#' @param choices Vector of observed choices for each trial, as integers from 1 to N_acc.
#' @param params_B Matrix (N_trials x N_acc) of base threshold (B) parameters.
#' @param params_A Matrix (N_trials x N_acc) of start-point variability range (A) parameters.
#' @param params_mu_drift Matrix (N_trials x N_acc) of mean drift rate parameters.
#' @param params_sigma_drift_sq Matrix (N_trials x N_acc) of drift rate variance parameters.
#' @param params_t0 Matrix (N_trials x N_acc) of non-decision time (t0) parameters.
#' @param params_s Matrix (N_trials x N_acc) of overall scaling (s) parameters.
#' @param min_log_lik Minimum log-likelihood value to return for invalid inputs or numerical issues.
#' @param cdf_abs_err Desired absolute error for the internal CDF calculations.
#' @param cdf_rel_err Desired relative error for the internal CDF calculations.
#' @param cdf_max_eval Maximum evaluations for internal CDF calculations.
#' @return A numeric vector of log-likelihood values for each trial.
#' @details This function computes the likelihood of observing each RT and choice pair,
#'   given the parameters for N DSWTN accumulators racing independently. Parameters can
#'   vary per trial and per accumulator.
// [[Rcpp::export]]
NumericVector loglik_RDM_DSWTN_race(
    NumericVector rts,
    IntegerVector choices,
    NumericMatrix params_B,
    NumericMatrix params_A,
    NumericMatrix params_mu_drift,
    NumericMatrix params_sigma_drift_sq,
    NumericMatrix params_t0,
    NumericMatrix params_s,
    double min_log_lik = -1e10,
    double spv_abs_err = 1e-6,
    double spv_rel_err = 1e-6,
    int    spv_max_eval = 1000
    // Note: Inner drift variability integration errors are currently fixed in rdm_dswtn_cdf_spv_integrand's call to dswtn_cdf_core
    // or inherited from spv_abs_err in rdm_dswtn_logpdf/cdf_core Case 3. This could be further parameterized if needed.
) {
    int n_trials = rts.size();
    if (n_trials == 0) return NumericVector(0);

    // Dimension checks
    if (choices.size() != n_trials || params_B.nrow() != n_trials ||
        params_A.nrow() != n_trials || params_mu_drift.nrow() != n_trials ||
        params_sigma_drift_sq.nrow() != n_trials || params_t0.nrow() != n_trials ||
        params_s.nrow() != n_trials) {
        Rcpp::stop("Input vector/matrix dimensions do not match number of trials.");
    }

    int n_acc = params_B.ncol();
    if (n_acc == 0) {
        if (n_trials > 0) Rcpp::stop("Number of accumulators (ncol of param matrices) is zero.");
        return NumericVector(n_trials);
    }
    if (params_A.ncol() != n_acc || params_mu_drift.ncol() != n_acc ||
        params_sigma_drift_sq.ncol() != n_acc || params_t0.ncol() != n_acc ||
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
        double sigma_drift_sq_winner_scaled = params_sigma_drift_sq(i, winning_acc_idx) / (s_winner * s_winner);
        double t0_winner = params_t0(i, winning_acc_idx);
        double t_adj_winner = rt - t0_winner;

        if (t_adj_winner <= 1e-10) {
            trial_log_likelihoods[i] = min_log_lik;
            continue;
        }

        double current_log_lik = rdm_dswtn_logpdf_core(
            t_adj_winner, B_winner_scaled, A_winner_scaled,
            mu_drift_winner_scaled, sigma_drift_sq_winner_scaled,
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

            double cdf_loser = rdm_dswtn_cdf_core(
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
}


// --- Structs for SPV integration ---
struct DSWTN_SPV_Integrand_Params_PDF {
    double t_adj;
    // k_center and A_spv are not needed here as current_k is the integration variable
    double mu_drift;
    double sigma_drift_sq;
    // A_spv is needed for normalization (1/A_spv) if done outside integrand
};

struct DSWTN_SPV_Integrand_Params_CDF {
    double t_adj;
    double mu_drift;
    double sigma_drift_sq;
    // A_spv for normalization
};


// --- Integrands for Start Point Variability ---

// Integrand for PDF with SPV: exp(dswtn_logpdf_core(t_adj, current_k, ...))
// Note: hcubature integrates the function directly. If A_spv is non-zero, the result
// of hcubature then needs to be divided by A_spv (or multiplied by 1/A_spv in integrand).
// For logpdf, we integrate pdf then take log.
int dswtn_pdf_spv_integrand(unsigned dim, const double* current_k_val, void* p, unsigned fdim, double* retval) {
    if (dim != 1 || fdim != 1) return 1; // Error
    DSWTN_SPV_Integrand_Params_PDF* params = static_cast<DSWTN_SPV_Integrand_Params_PDF*>(p);
    double current_k = current_k_val[0];

    if (current_k <= 0) { // Threshold must be positive
        retval[0] = 0.0;
        return 0;
    }
    // Calculate PDF for current_k (alpha for dswtn_logpdf_core)
    double log_pdf_val = dswtn_logpdf_core(params->t_adj, current_k, params->mu_drift, params->sigma_drift_sq);
    if (!std::isfinite(log_pdf_val)) {
        retval[0] = 0.0;
    } else {
        retval[0] = std::exp(log_pdf_val);
    }
    return 0; // Success
}

// Integrand for CDF with SPV: dswtn_cdf_core(t_adj, current_k, ...)
int dswtn_cdf_spv_integrand(unsigned dim, const double* current_k_val, void* p, unsigned fdim, double* retval) {
    if (dim != 1 || fdim != 1) return 1; // Error
    DSWTN_SPV_Integrand_Params_CDF* params = static_cast<DSWTN_SPV_Integrand_Params_CDF*>(p);
    double current_k = current_k_val[0];

    if (current_k <= 0) { // Threshold must be positive
        retval[0] = 0.0; // Or some other appropriate value if k can be 0 and means instant absorption for CDF context
                       // If current_k is boundary, and it's <=0, CDF should be 1 if t_adj > 0.
                       // Wald_cdf_classic handles alpha<=0 returning 1 if drift >0.
                       // dswtn_cdf_core calls wald_cdf_classic.
    }

    // Calculate CDF for current_k (alpha for dswtn_cdf_core)
    // Default integration params for inner CDF calculation for now
    retval[0] = dswtn_cdf_core(params->t_adj, current_k, params->mu_drift, params->sigma_drift_sq);
    return 0; // Success
}


// --- Core functions with Start Point Variability (SPV) ---

// Log-PDF for DSWTN with Start Point Variability (uniform over k_center +/- A_spv/2)
// k_center: central threshold, A_spv: range of start point variability
inline double dswtn_with_spv_logpdf_core(double t_adj, double k_center, double A_spv,
                                        double mu_drift, double sigma_drift_sq,
                                        double abs_err = 1e-6, double rel_err = 1e-6, size_t max_eval = 1000) {
    if (t_adj <= 1e-10) return R_NegInf;
    if (k_center <= 0 && (k_center + A_spv/2.0) <=0) return R_NegInf; // Entire SPV range is non-positive
    if (A_spv < 0) return R_NaN; // Invalid A_spv

    // If no start point variability, use the non-SPV version
    if (A_spv < 1e-7) { // Treat as no SPV
        if (k_center <= 0) return R_NegInf; // Central k must be positive if no SPV
        return dswtn_logpdf_core(t_adj, k_center, mu_drift, sigma_drift_sq);
    }

    DSWTN_SPV_Integrand_Params_PDF int_params;
    int_params.t_adj = t_adj;
    int_params.mu_drift = mu_drift;
    int_params.sigma_drift_sq = sigma_drift_sq;

    double lower_k = k_center - A_spv / 2.0;
    double upper_k = k_center + A_spv / 2.0;

    // Ensure lower_k is not less than a very small positive, as threshold k must be > 0.
    // If k_center - A_spv/2 is <=0, the integration range needs to start from a small positive value.
    // However, the integrand dswtn_pdf_spv_integrand itself handles current_k <= 0 returning 0.
    // So, we can allow lower_k to be its calculated value.
    // If lower_k < 0 and upper_k < 0, then this should result in 0 PDF.
    if (upper_k <= 1e-9) return R_NegInf; // Whole range non-positive
    if (lower_k < 1e-9) lower_k = 1e-9; // Clip integration start if it goes non-positive, to avoid issues with k in wald.
                                        // This means we are integrating over U(max(epsilon, k_center-A_spv/2), k_center+A_spv/2) effectively.
                                        // And then normalizing by full A_spv.
                                        // This is slightly different from true U(k_center-A_spv/2, k_center+A_spv/2) if range crosses 0.
                                        // A more proper way would be to let integrand handle k<=0 and integrate over full theoretical range.
                                        // The current integrand returns 0 if current_k <=0.

    double integral_val;
    double integral_err;
    double a_k[1] = {lower_k};
    double b_k[1] = {upper_k};

    int success = hcubature(dswtn_pdf_spv_integrand, &int_params, 1, a_k, b_k,
                            max_eval, abs_err, rel_err, &integral_val, &integral_err);

    if (success != 0) {
      // Rcpp::Rcout << "Warning: DSWTN+SPV PDF integration did not converge." << std::endl;
    }

    if (integral_val <= 1e-300) return R_NegInf; // Avoid log(0) from very small positive

    // Normalize by the width of the uniform distribution A_spv
    double final_pdf = integral_val / A_spv;
    if (final_pdf <= 1e-300) return R_NegInf;

    return std::log(final_pdf);
}


// CDF for DSWTN with Start Point Variability
inline double dswtn_with_spv_cdf_core(double t_adj, double k_center, double A_spv,
                                     double mu_drift, double sigma_drift_sq,
                                     double abs_err = 1e-6, double rel_err = 1e-6, size_t max_eval = 1000) {
    if (t_adj <= 0) return 0.0;
    // If k_center + A_spv/2 (max possible threshold) is <=0, effectively always hit, so CDF is 1.
    if ((k_center + A_spv / 2.0) <= 1e-9) return 1.0;
    if (A_spv < 0) return R_NaN;

    if (A_spv < 1e-7) { // No SPV
        if (k_center <= 0) return 1.0; // if k_center is 0 or less, and no variability, it's hit.
        return dswtn_cdf_core(t_adj, k_center, mu_drift, sigma_drift_sq, abs_err, rel_err, max_eval);
    }

    DSWTN_SPV_Integrand_Params_CDF int_params;
    int_params.t_adj = t_adj;
    int_params.mu_drift = mu_drift;
    int_params.sigma_drift_sq = sigma_drift_sq;

    double lower_k = k_center - A_spv / 2.0;
    double upper_k = k_center + A_spv / 2.0;

    // Integrand dswtn_cdf_spv_integrand calls dswtn_cdf_core, which calls wald_cdf_classic.
    // wald_cdf_classic handles current_k (as boundary_alpha) <= 0 by returning 1.0 if drift > 0.
    // So, integrating from a non-positive lower_k is fine.
    // If upper_k is also non-positive, then the integral should yield A_spv * 1.0, and result is 1.0.

    double integral_val;
    double integral_err;
    double a_k[1] = {lower_k};
    double b_k[1] = {upper_k};

    int success = hcubature(dswtn_cdf_spv_integrand, &int_params, 1, a_k, b_k,
                            max_eval, abs_err, rel_err, &integral_val, &integral_err);

    if (success != 0) {
      // Rcpp::Rcout << "Warning: DSWTN+SPV CDF integration did not converge." << std::endl;
    }

    // Normalize by the width of the uniform distribution A_spv
    double final_cdf = integral_val / A_spv;

    if (std::isnan(final_cdf) || final_cdf < 0.0) return 0.0;
    if (final_cdf > 1.0) return 1.0;

    return final_cdf;
}


// END OF DSWTN IMPLEMENTATION with SPV
// Original DSWTN functions (dswtn_logpdf_core, dswtn_cdf_core, rswtn_core) remain as they are called by these SPV versions.
// Rcpp exported functions will need to be updated to call these _with_spv_ versions and take appropriate parameters.


// --- RDM_DSWTN: Combines RDM-style SPV with DSWTN drift variability ---

// These are the new top-level core functions.
// They will use dswtn_logpdf_core and dswtn_cdf_core as the "inner" functions
// when drift variability is present.

// Parameters for the integrand when integrating over RDM-style SPV U(B, B+A)
struct RDM_DSWTN_SPV_Integrand_Params {
    double t_adj;
    double mu_drift;
    double sigma_drift_sq;
    // B and A define the integration range, not passed in struct here.
};

// Integrand for PDF: exp(dswtn_logpdf_core(t_adj, current_actual_k, mu_drift, sigma_drift_sq))
// current_actual_k is the integration variable, drawn from U(B, B+A)
int rdm_dswtn_pdf_spv_integrand(unsigned dim, const double* current_actual_k_val, void* p, unsigned fdim, double* retval) {
    if (dim != 1 || fdim != 1) return 1;
    RDM_DSWTN_SPV_Integrand_Params* params = static_cast<RDM_DSWTN_SPV_Integrand_Params*>(p);
    double current_k = current_actual_k_val[0];

    if (current_k <= 1e-9) { // Threshold for dswtn_logpdf_core must be positive
        retval[0] = 0.0;
        return 0;
    }
    double log_pdf_val = dswtn_logpdf_core(params->t_adj, current_k, params->mu_drift, params->sigma_drift_sq);
    if (!std::isfinite(log_pdf_val)) {
        retval[0] = 0.0;
    } else {
        retval[0] = std::exp(log_pdf_val);
    }
    return 0;
}

// Integrand for CDF: dswtn_cdf_core(t_adj, current_actual_k, mu_drift, sigma_drift_sq)
// current_actual_k is the integration variable, drawn from U(B, B+A)
int rdm_dswtn_cdf_spv_integrand(unsigned dim, const double* current_actual_k_val, void* p, unsigned fdim, double* retval) {
    if (dim != 1 || fdim != 1) return 1;
    RDM_DSWTN_SPV_Integrand_Params* params = static_cast<RDM_DSWTN_SPV_Integrand_Params*>(p);
    double current_k = current_actual_k_val[0];

    // dswtn_cdf_core internally handles current_k <= 0 by calling wald_cdf_classic,
    // which returns 1.0 if boundary <=0 and drift > 0.
    // So, no specific check for current_k <=0 needed here for the integrand value itself,
    // unless we want to restrict integration domain explicitly.
    // For now, rely on dswtn_cdf_core's handling.

    retval[0] = dswtn_cdf_core(params->t_adj, current_k, params->mu_drift, params->sigma_drift_sq,
                               1e-7, 1e-7, 1000); // Using tighter fixed defaults for inner integration
    return 0;
}


// Top-level Log-PDF for RDM_DSWTN model
// Parameters B, A, mu_drift, sigma_drift_sq are assumed to be already scaled by s if applicable.
inline double rdm_dswtn_logpdf_core(double t_adj, double B, double A,
                                    double mu_drift, double sigma_drift_sq,
                                    double spv_abs_err = 1e-6, double spv_rel_err = 1e-6, size_t spv_max_eval = 1000) {

    if (t_adj <= 1e-10) return R_NegInf;

    const double A_is_zero_threshold = 1e-7;
    const double sigmasq_is_zero_threshold = 1e-10; // as used in dswtn_logpdf_core

    bool no_A_var = (A < A_is_zero_threshold);
    bool no_drift_var = (sigma_drift_sq < sigmasq_is_zero_threshold);

    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald (like digt0). Threshold is B. Drift is mu_drift.
        if (B <= 1e-9) return R_NegInf; // Threshold must be positive
        if (mu_drift <= 1e-9 && B > 0) return R_NegInf; // No positive drift to positive boundary
        return std::log(digt0(t_adj, B, mu_drift)); // digt0(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV (like digt). Drift is mu_drift.
        // dWald calls digt(t_adj, B + 0.5*A, mu_drift, 0.5*A)
        // So, k_for_digt = B + 0.5*A, a_for_digt = 0.5*A
        if ((B + A) <= 1e-9) return R_NegInf; // Max threshold non-positive
        double k_digt = B + 0.5 * A;
        double a_digt = 0.5 * A;
        if (k_digt - a_digt < 0 && A > A_is_zero_threshold) { // if B is negative
             // digt handles negative B if A is large enough to make parts of U(B,B+A) positive.
             // However, the k for digt (center) should ideally be positive if a_digt is not dominant.
             // For safety, if B (the minimum threshold) < 1e-9, but B+A > 1e-9,
             // digt needs to be robust. The original digt seems to handle this.
        }
         // Ensure a_digt is not zero if A was meant to be non-zero.
        if (A < A_is_zero_threshold) a_digt = 0.0; // force if A is tiny

        double pdf_val = digt(t_adj, k_digt, mu_drift, a_digt);
        if (pdf_val <= 1e-300) return R_NegInf;
        return std::log(pdf_val);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: DSWTN with fixed threshold B (like original dswtn_logpdf_core).
        if (B <= 1e-9) return R_NegInf;
        return dswtn_logpdf_core(t_adj, B, mu_drift, sigma_drift_sq);
    } else {
        // Case 4: Full model - DSWTN with RDM-style SPV.
        // Integrate exp(dswtn_logpdf_core(t_adj, actual_k, mu_drift, sigma_drift_sq))
        // where actual_k ~ U(B, B+A).
        if (A <= 1e-9) { // Should have been caught by no_A_var, but defensive.
            if (B <= 1e-9) return R_NegInf;
            return dswtn_logpdf_core(t_adj, B, mu_drift, sigma_drift_sq);
        }
        if ((B + A) <= 1e-9) return R_NegInf; // Max threshold non-positive

        RDM_DSWTN_SPV_Integrand_Params int_params;
        int_params.t_adj = t_adj;
        int_params.mu_drift = mu_drift;
        int_params.sigma_drift_sq = sigma_drift_sq;

        double lower_actual_k = B;
        double upper_actual_k = B + A;

        // Integrand returns 0 if current_k <= 0.
        // If B itself is negative, need to ensure integration range is sensible or relies on integrand.
        // If B < 1e-9 and B+A > 1e-9, we integrate from B. Integrand handles non-positive k.
        // If B+A < 1e-9 (whole range non-positive), this will result in 0.

        double integral_val;
        double integral_err;
        double k_limits_a[1] = {lower_actual_k};
        double k_limits_b[1] = {upper_actual_k};

        int success = hcubature(rdm_dswtn_pdf_spv_integrand, &int_params, 1, k_limits_a, k_limits_b,
                                spv_max_eval, spv_abs_err, spv_rel_err, &integral_val, &integral_err);
        if (success != 0) {
            // Rcpp::Rcout << "Warning: RDM_DSWTN PDF integration did not converge." << std::endl;
        }

        if (integral_val <= 1e-300) return R_NegInf;

        double final_pdf = integral_val / A; // Normalize by width of U(B, B+A)
        if (final_pdf <= 1e-300) return R_NegInf;
        return std::log(final_pdf);
    }
}

// Top-level CDF for RDM_DSWTN model
inline double rdm_dswtn_cdf_core(double t_adj, double B, double A,
                                 double mu_drift, double sigma_drift_sq,
                                 double spv_abs_err = 1e-6, double spv_rel_err = 1e-6, size_t spv_max_eval = 1000) {
    if (t_adj <= 0) return 0.0;

    const double A_is_zero_threshold = 1e-7;
    const double sigmasq_is_zero_threshold = 1e-10;

    bool no_A_var = (A < A_is_zero_threshold);
    bool no_drift_var = (sigma_drift_sq < sigmasq_is_zero_threshold);

    if (no_A_var && no_drift_var) {
        // Case 1: Simple Wald CDF (like pigt0). Threshold B.
        if (B <= 1e-9) return 1.0; // Hit non-positive threshold immediately
        if (mu_drift <= 1e-9 && B > 0) return 0.0; // Cannot reach positive boundary
        return pigt0(t_adj, B, mu_drift); // pigt0(t, k, l)
    } else if (!no_A_var && no_drift_var) {
        // Case 2: Standard RDM with SPV CDF (like pigt).
        if ((B + A) <= 1e-9) return 1.0; // Max threshold non-positive
        double k_digt = B + 0.5 * A;
        double a_digt = 0.5 * A;
        if (A < A_is_zero_threshold) a_digt = 0.0;
        return pigt(t_adj, k_digt, mu_drift, a_digt);
    } else if (no_A_var && !no_drift_var) {
        // Case 3: DSWTN CDF with fixed threshold B.
        if (B <= 1e-9) return 1.0;
        return dswtn_cdf_core(t_adj, B, mu_drift, sigma_drift_sq, spv_abs_err, spv_rel_err, spv_max_eval); // Pass SPV errors as they are for the only integration here
    } else {
        // Case 4: Full model - DSWTN CDF with RDM-style SPV.
        // Integrate dswtn_cdf_core(t_adj, actual_k, mu_drift, sigma_drift_sq)
        // where actual_k ~ U(B, B+A).
        if (A <= 1e-9) { // Defensive
             if (B <= 1e-9) return 1.0;
             return dswtn_cdf_core(t_adj, B, mu_drift, sigma_drift_sq, spv_abs_err, spv_rel_err, spv_max_eval);
        }
        if ((B+A) <= 1e-9) return 1.0;


        RDM_DSWTN_SPV_Integrand_Params int_params;
        int_params.t_adj = t_adj;
        int_params.mu_drift = mu_drift;
        int_params.sigma_drift_sq = sigma_drift_sq;

        double lower_actual_k = B;
        double upper_actual_k = B + A;

        double integral_val;
        double integral_err;
        double k_limits_a[1] = {lower_actual_k};
        double k_limits_b[1] = {upper_actual_k};

        // dswtn_cdf_spv_integrand calls dswtn_cdf_core. dswtn_cdf_core handles k<=0 for its 'alpha' correctly.
        // So, we can integrate from B even if B is non-positive.

        int success = hcubature(rdm_dswtn_cdf_spv_integrand, &int_params, 1, k_limits_a, k_limits_b,
                                spv_max_eval, spv_abs_err, spv_rel_err, &integral_val, &integral_err);
        if (success != 0) {
            // Rcpp::Rcout << "Warning: RDM_DSWTN CDF integration did not converge." << std::endl;
        }

        double final_cdf = integral_val / A; // Normalize by width of U(B, B+A)

        if (std::isnan(final_cdf) || final_cdf < 0.0) return 0.0;
        if (final_cdf > 1.0) return 1.0;
        return final_cdf;
    }
}

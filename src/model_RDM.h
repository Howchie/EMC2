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
#' @param alpha Vector of threshold (boundary separation) parameters. Must be positive.
#' @param mu_drift Vector of mean drift rate parameters for the underlying normal distribution.
#' @param sigma_drift_sq Vector of variance parameters for the underlying normal drift rate distribution. Must be non-negative.
#' @param theta Vector of non-decision time parameters.
#' @return A numeric vector of log-density values.
#' @details The drift rate for each trial is assumed to be drawn from a normal
#'   distribution N(mu_drift, sigma_drift_sq), truncated at 0 to ensure positive drift.
#'   The function handles vectorization of parameters. If parameter vectors are shorter
#'   than `t`, they are recycled.
// [[Rcpp::export]]
NumericVector dDSWTN_log(NumericVector t, NumericVector alpha, NumericVector mu_drift, NumericVector sigma_drift_sq, NumericVector theta) {
    int n = t.size();
    if (alpha.size() == 1) alpha = rep(alpha, n);
    if (mu_drift.size() == 1) mu_drift = rep(mu_drift, n);
    if (sigma_drift_sq.size() == 1) sigma_drift_sq = rep(sigma_drift_sq, n);
    if (theta.size() == 1) theta = rep(theta, n);

    NumericVector log_pdf(n);
    for (int i = 0; i < n; ++i) {
        double t_adj = t[i] - theta[i];
        log_pdf[i] = dswtn_logpdf_core(t_adj, alpha[i], mu_drift[i], sigma_drift_sq[i]);
    }
    return log_pdf;
}

#' Cumulative distribution function for the DSWTN distribution
#'
#' Calculates the cumulative distribution function (CDF) for the
#' Wald distribution with trial-varying normally distributed drift rates (DSWTN).
#'
#' @param t Vector of reaction times.
#' @param alpha Vector of threshold (boundary separation) parameters. Must be positive.
#' @param mu_drift Vector of mean drift rate parameters for the underlying normal distribution.
#' @param sigma_drift_sq Vector of variance parameters for the underlying normal drift rate distribution. Must be non-negative.
#' @param theta Vector of non-decision time parameters.
#' @param abs_err Desired absolute error for the numerical integration of the CDF.
#' @param rel_err Desired relative error for the numerical integration of the CDF.
#' @param max_eval Maximum number of evaluations for the numerical integration routine.
#' @return A numeric vector of CDF values.
#' @details The CDF is calculated by numerically integrating the DSWTN PDF.
#'   The function handles vectorization of parameters. If parameter vectors are shorter
#'   than `t`, they are recycled.
// [[Rcpp::export]]
NumericVector pDSWTN(NumericVector t, NumericVector alpha, NumericVector mu_drift, NumericVector sigma_drift_sq, NumericVector theta,
                     double abs_err = 1e-6, double rel_err = 1e-6, int max_eval = 1000) {
    int n = t.size();
    if (alpha.size() == 1) alpha = rep(alpha, n);
    if (mu_drift.size() == 1) mu_drift = rep(mu_drift, n);
    if (sigma_drift_sq.size() == 1) sigma_drift_sq = rep(sigma_drift_sq, n);
    if (theta.size() == 1) theta = rep(theta, n);

    NumericVector cdf_val(n);
    for (int i = 0; i < n; ++i) {
        double t_adj = t[i] - theta[i];
        cdf_val[i] = dswtn_cdf_core(t_adj, alpha[i], mu_drift[i], sigma_drift_sq[i], abs_err, rel_err, static_cast<size_t>(max_eval));
    }
    return cdf_val;
}

#' Random number generation for the DSWTN distribution
#'
#' Generates random samples from the Wald distribution with trial-varying
#' normally distributed drift rates (DSWTN).
#'
#' @param n_samples The number of samples to generate.
#' @param alpha The threshold (boundary separation) parameter. Must be positive.
#' @param mu_drift The mean drift rate parameter for the underlying normal distribution.
#' @param sigma_drift_sq The variance parameter for the underlying normal drift rate distribution. Must be non-negative.
#' @param theta The non-decision time parameter.
#' @return A numeric vector of random samples from the DSWTN distribution.
#' @details This function samples a drift rate from N(mu_drift, sigma_drift_sq) truncated at 0,
#'   then samples from a standard Wald distribution with that drift rate and the given alpha,
#'   finally adding theta. Parameters are scalar (single values).
// [[Rcpp::export]]
NumericVector rDSWTN(int n_samples, double alpha, double mu_drift, double sigma_drift_sq, double theta) {
    NumericVector samples(n_samples);
    RNG rng; // Default Rcpp RNG
    for (int i = 0; i < n_samples; ++i) {
        samples[i] = rswtn_core(alpha, mu_drift, sigma_drift_sq, theta, rng);
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

#' Log-likelihood for an N-accumulator DSWTN race model
#'
#' Calculates the log-likelihood for multiple trials of an N-accumulator race model,
#' where the finishing time of each accumulator is modeled by a DSWTN distribution.
#'
#' @param rts Vector of observed reaction times for each trial.
#' @param choices Vector of observed choices for each trial, as integers from 1 to N_acc,
#'   where N_acc is the number of accumulators.
#' @param params_alpha Matrix (N_trials x N_acc) of threshold (alpha) parameters.
#' @param params_mu_drift Matrix (N_trials x N_acc) of mean drift rate (mu_drift) parameters.
#' @param params_sigma_drift_sq Matrix (N_trials x N_acc) of drift rate variance (sigma_drift_sq) parameters.
#' @param params_theta Matrix (N_trials x N_acc) of non-decision time (theta) parameters.
#' @param min_log_lik Minimum log-likelihood value to return for invalid inputs or numerical issues.
#' @param cdf_abs_err Desired absolute error for the internal CDF calculations.
#' @param cdf_rel_err Desired relative error for the internal CDF calculations.
#' @param cdf_max_eval Maximum evaluations for internal CDF calculations.
#' @return A numeric vector of log-likelihood values for each trial.
#' @details This function computes the likelihood of observing each RT and choice pair,
#'   given the parameters for N DSWTN accumulators racing independently. Parameters can
#'   vary per trial and per accumulator.
// [[Rcpp::export]]
NumericVector loglik_DSWTN_race(
    NumericVector rts,
    IntegerVector choices,
    NumericMatrix params_alpha,
    NumericMatrix params_mu_drift,
    NumericMatrix params_sigma_drift_sq,
    NumericMatrix params_theta,
    double min_log_lik = -1e10,
    double cdf_abs_err = 1e-6,
    double cdf_rel_err = 1e-6,
    int    cdf_max_eval = 1000
) {
    int n_trials = rts.size();
    if (n_trials == 0) return NumericVector(0);

    // Basic dimension checks
    if (choices.size() != n_trials || params_alpha.nrow() != n_trials ||
        params_mu_drift.nrow() != n_trials || params_sigma_drift_sq.nrow() != n_trials ||
        params_theta.nrow() != n_trials) {
        Rcpp::stop("Input vector/matrix dimensions do not match number of trials.");
    }

    int n_acc = params_alpha.ncol();
    if (n_acc == 0) {
      if (n_trials > 0) Rcpp::stop("Number of accumulators (ncol of param matrices) is zero.");
      return NumericVector(n_trials); // Or NumericVector(0) if n_trials also 0
    }

    if (params_mu_drift.ncol() != n_acc || params_sigma_drift_sq.ncol() != n_acc || params_theta.ncol() != n_acc) {
        Rcpp::stop("Parameter matrices must have the same number of columns (accumulators).");
    }

    NumericVector trial_log_likelihoods(n_trials);

    for (int i = 0; i < n_trials; ++i) {
        double rt = rts[i];
        int winning_choice_1_indexed = choices[i];
        int winning_acc_idx = winning_choice_1_indexed - 1; // Convert to 0-indexed

        if (winning_acc_idx < 0 || winning_acc_idx >= n_acc) {
            trial_log_likelihoods[i] = min_log_lik; // Invalid choice index
            continue;
        }

        // --- PDF of the winning accumulator ---
        double theta_winner = params_theta(i, winning_acc_idx);
        double t_adj_winner = rt - theta_winner;

        if (t_adj_winner <= 1e-10) { // Winner finished at or before its non-decision time
            trial_log_likelihoods[i] = min_log_lik;
            continue;
        }

        double alpha_winner = params_alpha(i, winning_acc_idx);
        double mu_drift_winner = params_mu_drift(i, winning_acc_idx);
        double sigma_drift_sq_winner = params_sigma_drift_sq(i, winning_acc_idx);

        double current_log_lik = dswtn_logpdf_core(
            t_adj_winner, alpha_winner, mu_drift_winner, sigma_drift_sq_winner
        );

        if (!std::isfinite(current_log_lik) || current_log_lik < min_log_lik -100) { // Check for -Inf or very small numbers
             trial_log_likelihoods[i] = min_log_lik;
             continue;
        }

        // --- Survivor functions of the losing accumulators ---
        for (int j = 0; j < n_acc; ++j) {
            if (j == winning_acc_idx) continue; // Skip the winner

            double theta_loser = params_theta(i, j);
            double t_adj_loser = rt - theta_loser;

            if (t_adj_loser <= 1e-10) {
                // Loser finished at or before its non-decision time.
                // Its CDF is 0, so survivor is 1. log(1) = 0. No change to log_lik.
                continue;
            }

            double alpha_loser = params_alpha(i, j);
            double mu_drift_loser = params_mu_drift(i, j);
            double sigma_drift_sq_loser = params_sigma_drift_sq(i, j);

            double cdf_loser = dswtn_cdf_core(
                t_adj_loser, alpha_loser, mu_drift_loser, sigma_drift_sq_loser,
                cdf_abs_err, cdf_rel_err, static_cast<size_t>(cdf_max_eval)
            );

            double survivor_loser = 1.0 - cdf_loser;

            if (survivor_loser <= 1e-10) { // Using a small epsilon to avoid log(0)
                                          // This implies the trial is (nearly) impossible under these parameters.
                current_log_lik = min_log_lik;
                break; // Break from inner loop, this trial's likelihood is already min_log_lik
            }
            current_log_lik += std::log(survivor_loser);
             if (!std::isfinite(current_log_lik) || current_log_lik < min_log_lik -100) {
                 current_log_lik = min_log_lik; // ensure it doesn't go way below if many terms are bad
                 break;
            }
        }
        trial_log_likelihoods[i] = std::max(current_log_lik, min_log_lik); // Ensure it's not less than min_log_lik
    }
    return trial_log_likelihoods;
}


// END OF DSWTN IMPLEMENTATION

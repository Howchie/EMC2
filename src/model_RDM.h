#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <RcppArmadillo.h>
#include "wald_functions.h"
#include "composite_functions.h"
#include "gaussian.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    double ti = t[i] - t0[i];
    if (ti <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt_impl(ti, B[i] + .5 * A[i], v[i], .5 * A[i]);
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
    double ti = t[i] - t0[i];
    if (ti <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt_impl(ti, B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return cdf;
}

// ==========================================================================
// RDMSWTN: Racing Diffusion Model with Shifted Wald / Truncated-Normal drift
// ==========================================================================
//
// Canonical Wald form: absorbing Brownian motion with drift mu, diffusion
// sigma, starting uniformly on [b-A, b] (so b is the upper threshold).
//
// Parameter hierarchy:
//   dwald / pwald   — Wald PDF/CDF with uniform start-point range A on [0,b].
//   dswtn / pswtn   — Wald + between-trial drift variability sv (SWTN), fixed threshold.
//   drdmswtn / prdmswtn — full model: SWTN + start-point variability A.
//
// All functions expect t_adj = t - t0 (non-decision time already removed).
// The `s` (diffusion) parameter is absorbed by the caller: pass s=1 once
// parameters have been pre-scaled by s.
// ==========================================================================

// --------------------------------------------------------------------------
// Wald (inverse-Gaussian) PDF with uniform start-point variability.
// b    = upper threshold (start-point ~ Unif(b-A, b), i.e. range A below b).
// mu   = drift rate.
// sigma = diffusion coefficient.
// A    = start-point variability range (0 = point start at b).
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double dwald(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
             bool log_out = false) {
  double x_lo = 0.0;
  double x_hi = x_lo + A;

  const double var = sigma * sigma * t;
  if (var <= FPM_EPSILON) return 0.0;

  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;

  // Point-start Wald at fixed boundary distance d = b - x.
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;
    if (d <= 0.0) return (t == 0.0) ? std::numeric_limits<double>::infinity() : 0.0;
    const double delta = d - mu * t;
    const double pdf_val = d / t * Gstar(var, delta);
    return log_out ? std::log(pdf_val) : pdf_val;
  }

  // Canonical SPV form: integrate over start-point x in [0, A].
  const double mu_new = b - mu * t;
  const double pdf_hi = Gstar(var, x_hi - mu_new);
  const double pdf_lo = Gstar(var, x_lo - mu_new);
  const double cdf_integral = Gstar_Integral(var, mu_new, x_lo, x_hi);

  const double term1 = mu * t * cdf_integral;
  const double term2 = var * (pdf_hi - pdf_lo);
  const double pdf_val = ((1.0 / t) * (term1 + term2)) / span;

  return log_out ? std::log(pdf_val) : pdf_val;
}

// [[Rcpp::export]]
double pwald(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
             bool log_out = false) {
  if (t <= FPM_EPSILON) return log_out ? -std::numeric_limits<double>::infinity() : 0.0;

  double x_lo = 0.0;
  double x_hi = x_lo + A;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;

  // Point-start fallback.
  if (span <= FPM_EPSILON) {
    const double sqrt_var = sigma * std::sqrt(t);
    const double dist = b - x_hi;
    if (sqrt_var <= FPM_EPSILON) return (mu > 0) ? 1.0 : 0.0;

    const double term1_arg = (mu * t - dist) / sqrt_var;
    const double term2_arg = (-mu * t - dist) / sqrt_var;
    const double cdf1 = gaussian_cdf(term1_arg, 0.0, 1.0);
    const double cdf2 = gaussian_cdf(term2_arg, 0.0, 1.0);
    const double exp_term = std::exp(2.0 * mu * dist / (sigma * sigma));

    double cdf_val = cdf1 + exp_term * cdf2;
    if (!std::isfinite(exp_term)) cdf_val = 1.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }

  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) {
    const double cdf_val = (mu > 0 && b <= x_hi) ? 1.0 : 0.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }

  const double k = 2.0 * mu / (sigma * sigma);

  // First term: ∫ Φ((x + (mu*t - b))/st) dx.
  const double a1 = 1.0 / st;
  const double c1 = (mu * t - b) / st;
  const double term1_integral = integrate_gaussian_cdf(a1, c1, x_lo, x_hi);

  // Second term: ∫ exp(k*(b-x)) * Φ((x - (mu*t+b))/st) dx.
  const double a2 = 1.0 / st;
  const double c2 = (-mu * t - b) / st;
  const double term2_integral =
    std::exp(k * b) * integrate_exp_times_normal_cdf(-k, a2, c2, x_lo, x_hi);

  double cdf_val = (term1_integral + term2_integral) / span;
  if (cdf_val < 0.0) cdf_val = 0.0;
  if (cdf_val > 1.0) cdf_val = 1.0;

  return log_out ? std::log(cdf_val) : cdf_val;
}

// --------------------------------------------------------------------------
// SWTN: Shifted Wald with Truncated-Normal drift variability.
// No start-point variability (fixed threshold = threshold param).
// sv  = SD of between-trial drift distribution N(mu_drift, sv^2) [truncated at 0].
// s   = within-trial diffusion coefficient.
// c   = lower truncation for drift (default 0 for positive-drift assumption).
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double dswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s = 1.0, double c = 0.0, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (sv < 0.0) return R_NaN;
  if (s <= 1e-10) return log_out ? R_NegInf : 0.0;

  // sv == 0: reduces to standard Wald
  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return dwald(t_adj, threshold, mu_drift, s, 0.0, log_out);
  }

  const double v  = sv * sv;
  const double s2 = s * s;
  const double tv = t_adj * v;
  const double den_common     = tv + s2;
  const double log_den_common = std::log(den_common);

  // Normalisation P(xi > c) for xi ~ N(mu_drift, sv^2)
  const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
  if (std::isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;

  const double term_log_threshold = std::log(threshold);
  const double term_log_denom     = -0.5 * (std::log(M_PI) + M_LN2 +
                                             3.0 * std::log(t_adj) + log_den_common);
  const double term_log_exp       = -(std::pow(threshold - mu_drift * t_adj, 2.0)) /
                                      (2.0 * t_adj * den_common);

  // Inner truncation: P(xi_new > c) for xi_new ~ N(mu_new, sigma_new^2)
  const double mu_new    = (threshold * v + mu_drift * s2) / den_common;
  const double sigma_new = std::sqrt(s2 * v / den_common);
  const double term_log_int = pnorm_std((c - mu_new) / sigma_new, false, true);

  const double log_pdf = term_log_threshold + term_log_denom +
                         (-log_prob_gt_c) + term_log_exp + term_log_int;
  if (std::isnan(log_pdf)) return log_out ? R_NegInf : 0.0;
  return log_out ? log_pdf : std::exp(log_pdf);
}

// --------------------------------------------------------------------------
// SWTN CDF.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double pswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s = 1.0, double c = 0.0, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (sv < 0.0) return R_NaN;

  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return pwald(t_adj, threshold, mu_drift, s, 0.0, log_out);
  }

  const double v   = sv * sv;
  const double s2  = s * s;
  const double tv  = t_adj * v;
  const double th2 = threshold * threshold;

  // P(xi > c) for xi ~ N(mu_drift, sv^2)
  const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
  if (std::isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;

  // Shared denominator term
  const double s_denom_new = std::sqrt(t_adj * (s2 + tv));
  const double rho_new     = (std::sqrt(t_adj) * sv) / std::sqrt(s2 + tv);

  // Term 1: bivariate normal
  const double h1 = (mu_drift * t_adj - threshold) / s_denom_new;
  const double k1 = (mu_drift - c) / sv;
  double term1;
  if (std::fabs(rho_new) > 0.9999 || std::fabs(h1) > 8.0 || std::fabs(k1) > 8.0)
    term1 = norm_cdf_2d(h1, k1, rho_new);
  else
    term1 = norm_cdf_2d_fast(h1, k1, rho_new);
  const double log_A = (term1 > 1e-300) ? std::log(term1) : R_NegInf;

  // Term 2: reflected bivariate normal
  const double log_exp_term = (2.0 * threshold * mu_drift / s2) +
                               (2.0 * th2 * v / (s2 * s2));
  const double mu_p  = mu_drift + 2.0 * threshold * v / s2;
  const double h2    = (-mu_p * t_adj - threshold) / s_denom_new;
  const double k2    = (mu_p - c) / sv;
  double term2;
  if (std::fabs(rho_new) > 0.9999 || std::fabs(h2) > 8.0 || std::fabs(k2) > 8.0)
    term2 = norm_cdf_2d(h2, k2, -rho_new);
  else
    term2 = norm_cdf_2d_fast(h2, k2, -rho_new);
  const double log_B = (term2 > 1e-300) ? (log_exp_term + std::log(term2)) : R_NegInf;

  const double log_numerator = log_sum_exp(log_A, log_B);
  const double log_cdf       = log_numerator - log_prob_gt_c;

  if (std::isnan(log_cdf)) return log_out ? R_NegInf : 0.0;
  double cdf_val = std::exp(log_cdf);
  cdf_val = std::fmax(0.0, std::fmin(1.0, cdf_val));
  return log_out ? std::log(cdf_val) : cdf_val;
}

// --------------------------------------------------------------------------
// Full RDMSWTN: SWTN + start-point variability.
// b  = upper threshold; threshold ~ Unif(b-A, b).
// Uses pre-cached GL20 nodes for the quadrature over [b-A, b].
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double drdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s = 1.0, double c = 0.0,
                int n_gauss_nodes = 20, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? R_NegInf : 0.0;

  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv < 1e-7);

  if (no_sv) {
    // sv=0: standard Wald with start-point variability A (dwald handles A=0 internally)
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return dwald(t_adj, b, mu_drift, s, A, log_out);
  } else if (no_A && !no_sv) {
    // SWTN with fixed threshold b
    return dswtn(t_adj, b, mu_drift, sv, s, c, log_out);
  } else {
    // Full model: integrate dswtn over threshold ~ Unif(b-A, b)
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl = (n_nodes == 20) ? gl20
                                    : Rcpp::as<Rcpp::List>(gauss_quad(n_nodes, "legendre"));
    const Rcpp::NumericVector gl_nodes   = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    // Map nodes from [-1,1] to [b-A, b]
    const double center     = b - 0.5 * A;
    const double half_width = 0.5 * A;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double thresh_j = center + half_width * gl_nodes[j];
      integral += gl_weights[j] * dswtn(t_adj, thresh_j, mu_drift, sv, s, c, false);
    }
    // Scale: GL on [-1,1]→[b-A,b] gives factor half_width; divide by A (uniform density)
    const double out_val = integral * half_width / A;   // = integral * 0.5
    if (out_val < 0.0 || std::isnan(out_val)) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(out_val) : out_val;
  }
}

// [[Rcpp::export]]
double prdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s = 1.0, double c = 0.0,
                int n_gauss_nodes = 20, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? 0.0 : 1.0;

  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv < 1e-7);

  if (no_sv) {
    // sv=0: standard Wald with start-point variability A (pwald handles A=0 internally)
    if (mu_drift <= 1e-7) return log_out ? R_NegInf : 0.0;
    return pwald(t_adj, b, mu_drift, s, A, log_out);
  } else if (no_A && !no_sv) {
    return pswtn(t_adj, b, mu_drift, sv, s, c, log_out);
  } else {
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl = (n_nodes == 20) ? gl20
                                    : Rcpp::as<Rcpp::List>(gauss_quad(n_nodes, "legendre"));
    const Rcpp::NumericVector gl_nodes   = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    const double center     = b - 0.5 * A;
    const double half_width = 0.5 * A;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double thresh_j = center + half_width * gl_nodes[j];
      integral += gl_weights[j] * pswtn(t_adj, thresh_j, mu_drift, sv, s, c, false);
    }
    double out_val = integral * half_width / A;
    out_val = std::fmax(0.0, std::fmin(1.0, out_val));
    return log_out ? std::log(out_val) : out_val;
  }
}

// --------------------------------------------------------------------------
// Vectorised R-callable exports (match zachdev naming so R code is portable).
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector t0, NumericVector sv,
                       NumericVector s = 1.0, NumericVector c = 0.0,
                       int n_gauss_nodes = 20, bool log_out = false) {
  const int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { pdf[i] = log_out ? R_NegInf : 0.0; continue; }
    pdf[i] = drdmswtn(t_adj, pick(b, i), pick(v, i), pick(A, i),
                      pick(sv, i), pick(s, i), pick(c, i), n_gauss_nodes, log_out);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector t0, NumericVector sv,
                       NumericVector s = 1.0, NumericVector c = 0.0,
                       int n_gauss_nodes = 20, bool log_out = false) {
  const int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { cdf[i] = log_out ? R_NegInf : 0.0; continue; }
    cdf[i] = prdmswtn(t_adj, pick(b, i), pick(v, i), pick(A, i),
                      pick(sv, i), pick(s, i), pick(c, i), n_gauss_nodes, log_out);
  }
  return cdf;
}

// --------------------------------------------------------------------------
// Vectorised adapters for the EMC2 race-likelihood machinery.
// Column layout (from Ttransform + p_types reorder):
//   v=0, B=1, A=2, t0=3, s=4, sv=5  [pContaminant and b follow, ignored here]
// Scaling convention: absorb s inside kernel; pass s=1 to drdmswtn/prdmswtn.
// --------------------------------------------------------------------------

#endif

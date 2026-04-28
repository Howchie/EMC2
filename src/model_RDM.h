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
                    NumericVector B, NumericVector A, NumericVector t0,
                    bool log_out = false){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    double ti = t[i] - t0[i];
    if (ti <= 0){
      pdf[i] = log_out? R_NegInf : 0.0;
    } else {
      pdf[i] = digt_impl(ti, B[i] + .5 * A[i], v[i], .5 * A[i], log_out);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0,
                    bool log_out = false){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    double ti = t[i] - t0[i];
    if (ti <= 0){
      cdf[i] = log_out? 0.0 : R_NegInf;
    } else {
      cdf[i] = pigt_impl(ti, B[i] + .5 * A[i], v[i], .5 * A[i], log_out);
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
double dwald_old(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
             bool log_out = false, double k = 0.0) {
  double x_lo = 0.0; double x_hi = x_lo + A; 
  const double var = sigma * sigma * t;
  if (var <= FPM_EPSILON || t <= FPM_EPSILON || sigma <= 0.0 || k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }

  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;

  // Point-start Wald at fixed boundary distance d = b - x.
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;
    if (d <= 0.0) {
      if (t == 0.0) {
        return log_out ? R_PosInf : std::numeric_limits<double>::infinity();
      }
      return log_out ? R_NegInf : 0.0;
    }
    const double delta = d - mu * t;
    const double log_pdf = std::log(d / t * Gstar(var, delta)) - k*t; 
    return log_out ? log_pdf : std::exp(log_pdf);
  }
    const double mu_new = b - mu * t;
    const double pdf_hi = Gstar(var, x_hi - mu_new);
    const double pdf_lo = Gstar(var, x_lo - mu_new);
    const double cdf_integral = Gstar_Integral(var, mu_new, x_lo, x_hi);
    const double term1 = mu * t * cdf_integral; 
    const double term2 = var * (pdf_hi - pdf_lo); 
    // todo is there a more efficient handling of the log here? It seems important for numerical stability because k could be very small 
    const double log_pdf = std::log(((1.0 / t) * (term1 + term2)) / span) - k*t;
    return log_out ? log_pdf : std::exp(log_pdf);
}

// [[Rcpp::export]]
double dwald(double t, double b, double mu, double sigma = 1.0, double A = 0.0, double k = 0.0,
             bool log_out = false) {
  double x_lo = 0.0;
  double x_hi = x_lo + A;

  const double var = sigma * sigma * t;
  if (var <= FPM_EPSILON || t <= FPM_EPSILON || sigma <= 0.0 || k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }

  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;
  double log_pdf = R_NegInf;
  // Point-start Wald at fixed boundary distance d = b - x.
  const double logt = std::log(t);
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;
    if (d <= 0.0) {
      if (t == 0.0) {
        return log_out ? R_PosInf : std::numeric_limits<double>::infinity();
      }
      return log_out ? R_NegInf : 0.0;
    }
    const double delta = d - mu * t;
    log_pdf = std::log(d) - logt + Gstar(var, delta, true) - k * t;
    return log_out ? log_pdf : std::exp(log_pdf);
  }

  // Canonical SPV form: integrate over start-point x in [0, A].
  const double mu_new = b - mu * t;
  const double log_pdf_hi = Gstar(var, x_hi - mu_new, true);
  const double log_pdf_lo = Gstar(var, x_lo - mu_new, true);
  const double log_cdf_integral = Gstar_Integral(var, mu_new, x_lo, x_hi, true);

  // term1 = mu * t * cdf_integral
  signed_log term1 = make_signed_log(R_NegInf, 0);
  if (std::abs(mu) > FPM_EPSILON && log_cdf_integral != R_NegInf) {
    term1 = make_signed_log(
      std::log(std::abs(mu)) + logt + log_cdf_integral,
      mu > 0.0 ? 1 : -1
    );
  }
  // term2 = var * (pdf_hi - pdf_lo)
  signed_log pdf_hi = make_signed_log(log_pdf_hi, 1);
  signed_log pdf_lo = make_signed_log(log_pdf_lo, 1);
  signed_log pdf_diff = signed_log_sub(pdf_hi, pdf_lo);

  signed_log term2 = make_signed_log(R_NegInf, 0);
  if (pdf_diff.sign != 0) {
    term2 = make_signed_log(std::log(var) + pdf_diff.log_abs, pdf_diff.sign);
  }
  signed_log total = signed_log_add(term1, term2);

  if (total.sign <= 0 || total.log_abs == R_NegInf) {
    return log_out ? R_NegInf : 0.0;
  }

  log_pdf =
    total.log_abs -
    logt -
    std::log(span) -
    k * t;

  return log_out ? log_pdf : std::exp(log_pdf);
}

// --------------------------------------------------------------------------
// Wald / killed Wald CDF with uniform start-point variability.
// Returns sub-CDF:
//   P(T_wald <= t AND T_kill > T_wald)
// k = 0 gives standard Wald, nested.
// t = Inf returns total finite-response probability.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double pwald(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
             double k = 0.0, bool log_out = false) {
  auto finish = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;

    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;

    return log_out ? log_p : std::exp(log_p);
  };

  if (sigma <= 0.0 || k < 0.0) {
    return NA_REAL;
  }

  if (t <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  double x_lo = 0.0;
  double x_hi = x_lo + A;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);

  const double span = x_hi - x_lo;
  const double sig2 = sigma * sigma;

  const double nu = std::sqrt(mu * mu + 2.0 * sig2 * k);
  const double eta1 = (mu - nu) / sig2;

  auto log_p_inf = [&]() {
    if (span <= FPM_EPSILON) {
      const double d = b - x_hi;
      if (d <= 0.0) return 0.0;
      return eta1 * d;
    }

    if (std::abs(eta1) < FPM_EPSILON) {
      return 0.0;
    }

    // ∫ exp(eta1 * (b - x)) dx / span
    // eta1 is <= 0 for k >= 0.
    const double q = -eta1;

    if (q <= FPM_EPSILON) {
      return 0.0;
    }

    return
      eta1 * b +
      log_diff_exp(q * x_lo, q * x_hi) -
      std::log(q) -
      std::log(span);
  };

  if (!std::isfinite(t)) {
    return finish(log_p_inf());
  }

  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  // Point-start case.
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;

    if (d <= 0.0) {
      return log_out ? 0.0 : 1.0;
    }

    const double log_prefactor = d * (mu - nu) / sig2;

    const double term1_arg = (nu * t - d) / st;
    const double term2_arg = (-nu * t - d) / st;

    const double log_cdf1 = pnorm_std(term1_arg, true, true);
    const double log_cdf2 = pnorm_std(term2_arg, true, true);

    const double log_exp_term = 2.0 * nu * d / sig2;

    const double log_cdf_nu =
      log_sum_exp(log_cdf1, log_exp_term + log_cdf2);

    return finish(log_prefactor + log_cdf_nu);
  }

  // SPV case.
  const double a = 1.0 / st;

  // First term:
  // exp[(mu - nu)(b-x)/sig2] *
  // Phi((x + nu*t - b) / st)
  const double c1 = (nu * t - b) / st;

  const double log_term1 =
    eta1 * b +
    log_integrate_exp_times_normal_cdf(-eta1, a, c1, x_lo, x_hi);

  // Second term:
  // exp[(mu + nu)(b-x)/sig2] *
  // Phi((x - nu*t - b) / st)
  const double eta2 = (mu + nu) / sig2;
  const double c2 = (-nu * t - b) / st;

  const double log_term2 =
    eta2 * b +
    log_integrate_exp_times_normal_cdf(-eta2, a, c2, x_lo, x_hi);

  const double log_cdf_val =
    log_sum_exp(log_term1, log_term2) - std::log(span);

  return finish(log_cdf_val);
}

// [[Rcpp::export]]
double pwald_old(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
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
    const double cdf1 = pnorm_std(term1_arg, true, false);
    const double cdf2 = pnorm_std(term2_arg, true, false);
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
// Killed / defective Wald CDF with uniform start-point variability.
// Returns sub-CDF:
//   P(T_wald <= t AND T_kill > T_wald)
// where T_kill ~ Exponential(k).
//
// k = 0 gives the ordinary Wald CDF.
// t = Inf returns total finite-response probability.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double pwald_old_k(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
             double k = 0.0, bool log_out = false) {

  if (sigma <= 0.0 || k < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double x_lo = 0.0;
  double x_hi = x_lo + A;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);

  const double span = x_hi - x_lo;
  const double sig2 = sigma * sigma;

  // Tilted drift induced by exponential killing.
  const double nu = std::sqrt(mu * mu + 2.0 * sig2 * k);

  // Total finite-response probability:
  // F_k(Inf) = E[exp(-k * T_wald)]
  auto p_inf = [&]() {
    if (span <= FPM_EPSILON) {
      const double d = b - x_hi;
      if (d <= 0.0) return 1.0;

      const double eta = (mu - nu) / sig2;
      return std::exp(eta * d);
    }

    const double eta = (mu - nu) / sig2;

    if (std::abs(eta) < FPM_EPSILON) {
      return 1.0;
    }

    // ∫ exp(eta * (b - x)) dx over [x_lo, x_hi]
    const double integral =
      std::exp(eta * b) *
      (std::exp(-eta * x_lo) - std::exp(-eta * x_hi)) / eta;

    return integral / span;
  };

  if (!std::isfinite(t)) {
    return log_out? R_NegInf : 0.0;
  }

  if (t <= FPM_EPSILON) {
    return log_out? R_NegInf : 0.0;
  }

  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) {
    return log_out? R_NegInf : 0.0;
  }

  // ------------------------------------------------------------------------
  // Point-start case.
  // F_k(t; mu) =
  //   exp(d * (mu - nu) / sigma^2) * F_wald(t; nu)
  // ------------------------------------------------------------------------
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;

    if (d <= 0.0) {
      return log_out? 0.0 : 1.0;
    }

    const double prefactor = std::exp(d * (mu - nu) / sig2);

    const double term1_arg = (nu * t - d) / st;
    const double term2_arg = (-nu * t - d) / st;

    const double cdf1 = pnorm_std(term1_arg, true, false);
    const double cdf2 = pnorm_std(term2_arg, true, false);

    const double exp_arg = 2.0 * nu * d / sig2;
    double cdf_nu;

    if (exp_arg > 700.0) {
      // Same spirit as your existing safeguard.
      // In this regime the Wald CDF with positive drift nu is effectively 1
      // for practical MCMC purposes.
      cdf_nu = 1.0;
    } else {
      cdf_nu = cdf1 + std::exp(exp_arg) * cdf2;
    }

    return log_out? std::log(prefactor) + std::log(cdf_nu) : prefactor * cdf_nu;
  }

  // ------------------------------------------------------------------------
  // Uniform start-point variability case.
  //
  // For d = b - x:
  //
  // F_k(t)
  // = 1/A ∫ exp[d(mu - nu)/sigma^2] F_wald(t; d, nu, sigma) dx
  //
  // This becomes two integrals of exp(kx) * Phi(ax+c), using your helper.
  // ------------------------------------------------------------------------

  const double a = 1.0 / st;

  // First Wald-CDF term:
  // exp[(mu - nu)(b-x)/sig2] *
  // Phi((x + nu*t - b) / st)
  const double eta1 = (mu - nu) / sig2;
  const double c1 = (nu * t - b) / st;

  const double term1 =
    std::exp(eta1 * b) *
    integrate_exp_times_normal_cdf(-eta1, a, c1, x_lo, x_hi);

  // Second Wald-CDF term:
  // exp[(mu - nu)(b-x)/sig2] *
  // exp[2*nu*(b-x)/sig2] *
  // Phi((x - nu*t - b) / st)
  //
  // Combined exponent:
  // exp[(mu + nu)(b-x)/sig2]
  const double eta2 = (mu + nu) / sig2;
  const double c2 = (-nu * t - b) / st;

  const double term2 =
    std::exp(eta2 * b) *
    integrate_exp_times_normal_cdf(-eta2, a, c2, x_lo, x_hi);

  const double cdf = (term1 + term2) / span;

  if (cdf < 0.0) {
    return log_out ? R_NegInf : 0.0;
  } else if (cdf > 1.0){
    return log_out ? 0.0 : 1.0;
  }
  return log_out ? std::log(cdf) : cdf;
}

// --------------------------------------------------------------------------
// GBM first-passage with physical-space start-point variability.
// Process: dX_t = mu X_t dt + sigma X_t dW_t, start X0 ~ Unif(1, 1 + A).
// Boundary: absorb at X_t = b (b > 0).
// In log-space this is a Wald kernel with drift (mu - 0.5*sigma^2), diffusion
// sigma, and start range [log(1), log(1 + A)] weighted by exp(x0).
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double dgbm(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
            bool log_out = false) {
  if (!(t > FPM_EPSILON) || !(sigma > FPM_EPSILON) || !(b > FPM_EPSILON) || !(A >= 0.0)) {
    return log_out ? R_NegInf : 0.0;
  }

  const double x_lo = 0.0;
  const double x_hi = std::log1p(A);
  const double log_b = std::log(b);
  const double log_mu = mu - 0.5 * sigma * sigma;
  const double var = sigma * sigma * t;
  const double span = x_hi - x_lo;

  // Point-start fallback: A = 0 -> X0 = 1.
  if (span <= FPM_EPSILON) {
    const double d = log_b - x_hi;
    if (d <= 0.0) return (t == 0.0) ? std::numeric_limits<double>::infinity() : (log_out ? R_NegInf : 0.0);
    const double delta = d - log_mu * t;
    const double pdf_val = d / t * Gstar(var, delta);
    if (!(pdf_val > 0.0) || !std::isfinite(pdf_val)) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(pdf_val) : pdf_val;
  }

  // Normalization for start density on [1, 1 + A] mapped to log-space.
  const double norm_const = A;
  if (!(norm_const > FPM_EPSILON)) return log_out ? R_NegInf : 0.0;

  // Complete-the-square form of the weighted Gaussian integral.
  const double mu_new_exp = log_b - log_mu * t + var;
  const double exp_factor = std::exp(log_b - log_mu * t + 0.5 * var);
  const double pdf_hi = Gstar(var, x_hi - mu_new_exp);
  const double pdf_lo = Gstar(var, x_lo - mu_new_exp);
  const double cdf_integral = Gstar_Integral(var, mu_new_exp, x_lo, x_hi);

  const double integral_term1 = (log_mu * t - var) * cdf_integral;
  const double integral_term2 = var * (pdf_hi - pdf_lo);
  const double integral_result = integral_term1 + integral_term2;
  const double pdf_val = (exp_factor * integral_result) / (norm_const * t);

  if (!(pdf_val > 0.0) || !std::isfinite(pdf_val)) return log_out ? R_NegInf : 0.0;
  return log_out ? std::log(pdf_val) : pdf_val;
}

// [[Rcpp::export]]
double pgbm(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
            bool log_out = false) {
  if (!(b > FPM_EPSILON) || !(sigma > FPM_EPSILON) || !(A >= 0.0)) {
    return log_out ? R_NegInf : 0.0;
  }
  if (t <= FPM_EPSILON) return log_out ? R_NegInf : 0.0;

  const double x_lo = 0.0;
  const double x_hi = std::log1p(A);
  const double log_b = std::log(b);
  const double log_mu = mu - 0.5 * sigma * sigma;
  const double st = sigma * std::sqrt(t);
  const double norm_const = A;

  // Point-start fallback: A = 0 -> X0 = 1.
  if (!(norm_const > FPM_EPSILON)) {
    const double dist = log_b - x_hi;
    if (st <= FPM_EPSILON) {
      const double cdf_val = (log_mu > 0.0) ? 1.0 : 0.0;
      return log_out ? std::log(cdf_val) : cdf_val;
    }
    const double term1_arg = (log_mu * t - dist) / st;
    const double term2_arg = (-log_mu * t - dist) / st;
    const double cdf1 = pnorm_std(term1_arg, true, false);
    const double cdf2 = pnorm_std(term2_arg, true, false);
    const double exp_term = std::exp(2.0 * log_mu * dist / (sigma * sigma));

    double cdf_val = cdf1 + exp_term * cdf2;
    if (!std::isfinite(exp_term)) cdf_val = 1.0;
    cdf_val = std::fmax(0.0, std::fmin(1.0, cdf_val));
    if (cdf_val <= 0.0) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }

  if (st <= FPM_EPSILON) {
    const double cdf_val = (log_mu > 0.0) ? 1.0 : 0.0;
    return log_out ? std::log(cdf_val) : cdf_val;
  }

  const double k_log = 2.0 * log_mu / (sigma * sigma);

  // (1/A) * ∫_{log(1)}^{log(1+A)} F(t | x0) * exp(x0) dx0
  const double term1_integral =
    integrate_exp_times_normal_cdf(1.0, 1.0 / st, (log_mu * t - log_b) / st, x_lo, x_hi);
  const double term2_integral =
    std::exp(k_log * log_b) *
    integrate_exp_times_normal_cdf(1.0 - k_log, 1.0 / st, (-log_mu * t - log_b) / st, x_lo, x_hi);

  double cdf_val = (term1_integral + term2_integral) / norm_const;
  cdf_val = std::fmax(0.0, std::fmin(1.0, cdf_val));
  if (cdf_val <= 0.0) return log_out ? R_NegInf : 0.0;
  return log_out ? std::log(cdf_val) : cdf_val;
}

// --------------------------------------------------------------------------
// SWTN: Shifted Wald with Truncated-Normal drift variability.
// No start-point variability (fixed threshold = threshold param).
// sv  = SD of between-trial drift distribution N(mu_drift, sv^2) [truncated at 0].
// s   = within-trial diffusion coefficient.
// c   = lower truncation for drift (default 0 for positive-drift assumption).
// --------------------------------------------------------------------------
inline double dswtn_core(double t_adj, double threshold, double mu_drift,
                         double sv, double s, double lambda, double c,
                         double log_prob_gt_c, bool log_out);

inline Rcpp::List get_gl_nodes_weights(int n_gauss_nodes);

// [[Rcpp::export]]
double dswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s = 1.0, double lambda=0.0, double c = 0.0, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (sv < 0.0) return R_NaN;
  if (s <= 1e-10) return log_out ? R_NegInf : 0.0;

  // sv == 0: reduces to standard Wald
  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return dwald(t_adj, threshold, mu_drift, s, 0.0, lambda, log_out);
  }

  // Normalisation P(xi > c) for xi ~ N(mu_drift, sv^2)
  const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
  if (emc2_isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;

  return dswtn_core(t_adj, threshold, mu_drift, sv, s, lambda, c, log_prob_gt_c, log_out);
}

// --------------------------------------------------------------------------
// SWTN CDF.
// --------------------------------------------------------------------------
inline Rcpp::List get_gl_nodes_weights(int n_gauss_nodes) {
  const int n_nodes = std::max(1, n_gauss_nodes);
  return (n_nodes == 20) ? get_gl20()
                         : Rcpp::as<Rcpp::List>(gauss_quad(n_nodes, "legendre"));
}

inline double dswtn_core(double t_adj, double threshold, double mu_drift,
                         double sv, double s, double lambda, double c,
                         double log_prob_gt_c, bool log_out) {
  const double v  = sv * sv;
  const double s2 = s * s;
  const double tv = t_adj * v;
  const double den_common     = tv + s2;
  const double log_den_common = std::log(den_common);

  const double term_log_threshold = std::log(threshold);
  const double term_log_denom     = -0.5 * (std::log(M_PI) + M_LN2 +
                                             3.0 * std::log(t_adj) + log_den_common);
  const double term_log_exp       = -(std::pow(threshold - mu_drift * t_adj, 2.0)) /
                                      (2.0 * t_adj * den_common);

  const double mu_new    = (threshold * v + mu_drift * s2) / den_common;
  const double sigma_new = std::sqrt(s2 * v / den_common);
  const double term_log_int = pnorm_std((c - mu_new) / sigma_new, false, true);

  const double log_pdf = term_log_threshold + term_log_denom +
                         (-log_prob_gt_c) + term_log_exp + term_log_int - lambda * t_adj;
  if (ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double pswtn_killed_quad(double t_adj, double threshold, double mu_drift,
                                double sv, double s, double c, double lambda,
                                int n_gauss_nodes = 20) {
  const int n_nodes = std::max(1, n_gauss_nodes);
  Rcpp::List gl = get_gl_nodes_weights(n_nodes);
  const Rcpp::NumericVector gl_nodes   = gl["nodes"];
  const Rcpp::NumericVector gl_weights = gl["weights"];

  const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
  if (!(alpha < 1.0)) return 0.0;

  double integral = 0.0;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);  // map [-1,1] -> [0,1]
    const double p = alpha + (1.0 - alpha) * u;
    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    integral += gl_weights[j] * pwald(t_adj, threshold, drift_j, s, 0.0, lambda, false);
  }

  return 0.5 * integral;
}

inline double pswtn_killed_inf_quad(double threshold, double mu_drift,
                                    double sv, double s, double c, double lambda,
                                    int n_gauss_nodes = 20,
                                    bool log_out = false) {
  auto finish_log = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;
    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;
    return log_out ? log_p : std::exp(log_p);
  };

  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (s <= 1e-10 || sv < 0.0 || lambda < 0.0) return NA_REAL;

  const double s2 = s * s;

  if (sv <= 1e-10) {
    const double nu = std::sqrt(mu_drift * mu_drift + 2.0 * s2 * lambda);
    return finish_log(threshold * (mu_drift - nu) / s2);
  }

  const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
  if (!(alpha < 1.0)) return log_out ? R_NegInf : 0.0;

  const int n_nodes = std::max(1, n_gauss_nodes);
  Rcpp::List gl = get_gl_nodes_weights(n_nodes);
  const Rcpp::NumericVector gl_nodes   = gl["nodes"];
  const Rcpp::NumericVector gl_weights = gl["weights"];

  double log_integral = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);
    double p = alpha + (1.0 - alpha) * u;
    p = std::fmax(0.0, std::fmin(std::nextafter(1.0, 0.0), p));

    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    const double nu_j = std::sqrt(drift_j * drift_j + 2.0 * s2 * lambda);
    const double log_hit_j = threshold * (drift_j - nu_j) / s2;

    if (gl_weights[j] > 0.0) {
      log_integral = log_sum_exp(log_integral, std::log(gl_weights[j]) + log_hit_j);
    }
  }

  return finish_log(log_integral - M_LN2);
}

inline double log_uniform_threshold_hit(double eta, double b, double A) {
  if (A <= 1e-10) {
    const double d = std::fmax(b, 0.0);
    return eta * d;
  }

  const double d_lo_raw = b - A;
  const double d_hi_raw = b;

  if (d_hi_raw <= 0.0) return 0.0;

  const double zero_width = std::fmax(0.0, std::fmin(d_hi_raw, 0.0) - d_lo_raw);
  const double d_lo = std::fmax(d_lo_raw, 0.0);
  const double d_hi = d_hi_raw;

  double log_pos_integral = R_NegInf;
  if (d_hi > d_lo) {
    if (std::abs(eta) < 1e-12) {
      log_pos_integral = std::log(d_hi - d_lo);
    } else if (eta < 0.0) {
      log_pos_integral = log_diff_exp(eta * d_lo, eta * d_hi) - std::log(-eta);
    } else {
      log_pos_integral = log_diff_exp(eta * d_hi, eta * d_lo) - std::log(eta);
    }
  }

  double log_total = log_pos_integral;
  if (zero_width > 0.0) {
    log_total = log_sum_exp(log_total, std::log(zero_width));
  }

  return log_total - std::log(A);
}

inline double prdmswtn_killed_inf_quad(double b, double mu_drift, double A,
                                       double sv, double s, double c, double lambda,
                                       int n_gauss_nodes = 20,
                                       bool log_out = false) {
  auto finish_log = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;
    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;
    return log_out ? log_p : std::exp(log_p);
  };

  if (b <= 1e-10) return log_out ? 0.0 : 1.0;
  if (s <= 1e-10 || sv < 0.0 || A < 0.0 || lambda < 0.0) return NA_REAL;

  const double s2 = s * s;

  if (sv <= 1e-10) {
    const double nu = std::sqrt(mu_drift * mu_drift + 2.0 * s2 * lambda);
    const double eta = (mu_drift - nu) / s2;
    return finish_log(log_uniform_threshold_hit(eta, b, A));
  }

  const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
  if (!(alpha < 1.0)) return log_out ? R_NegInf : 0.0;

  const int n_nodes = std::max(1, n_gauss_nodes);
  Rcpp::List gl = get_gl_nodes_weights(n_nodes);
  const Rcpp::NumericVector gl_nodes   = gl["nodes"];
  const Rcpp::NumericVector gl_weights = gl["weights"];

  double log_integral = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);
    double p = alpha + (1.0 - alpha) * u;
    p = std::fmax(0.0, std::fmin(std::nextafter(1.0, 0.0), p));

    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    const double nu_j = std::sqrt(drift_j * drift_j + 2.0 * s2 * lambda);
    const double eta_j = (drift_j - nu_j) / s2;
    const double log_hit_j = log_uniform_threshold_hit(eta_j, b, A);

    if (gl_weights[j] > 0.0) {
      log_integral = log_sum_exp(log_integral, std::log(gl_weights[j]) + log_hit_j);
    }
  }

  return finish_log(log_integral - M_LN2);
}

// [[Rcpp::export]]
double pswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s = 1.0, double c = 0.0, double lambda = 0.0,bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (sv < 0.0) return R_NaN;

  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return pwald(t_adj, threshold, mu_drift, s, 0.0, lambda, log_out);
  }

  if (lambda > 1e-10) {
    if (!emc2_isfinite(t_adj)) {
      return pswtn_killed_inf_quad(threshold, mu_drift, sv, s, c, lambda, 20, log_out);
    }
    const double cdf = pswtn_killed_quad(t_adj, threshold, mu_drift, sv, s, c, lambda, 20);
    if (!(cdf > 0.0)) return log_out ? R_NegInf : 0.0;
    if (cdf >= 1.0) return log_out ? 0.0 : 1.0;
    return log_out ? std::log(cdf) : cdf;
  }

  const double v   = sv * sv;
  const double s2  = s * s;
  const double tv  = t_adj * v;
  const double th2 = threshold * threshold;

  // P(xi > c) for xi ~ N(mu_drift, sv^2)
  const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
  if (emc2_isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;

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

  if (ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  double cdf = std::exp(log_cdf);
  if (cdf < 0.0) {
    return log_out ? R_NegInf : 0.0;
  } else if (cdf > 1.0){
    return log_out ? 0.0 : 1.0;
  }
  return log_out ? std::log(cdf) : cdf;
}

// --------------------------------------------------------------------------
// Full RDMSWTN: SWTN + start-point variability.
// b  = upper threshold; threshold ~ Unif(b-A, b).
// Uses pre-cached GL20 nodes for the quadrature over [b-A, b].
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double drdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s = 1.0, double lambda = 0.0, double c = 0.0,
                int n_gauss_nodes = 20, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? R_NegInf : 0.0;

  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv < 1e-7);

  if (no_sv) {
    // sv=0: standard Wald with start-point variability A (dwald handles A=0 internally)
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return dwald(t_adj, b, mu_drift, s, A, lambda, log_out);
  } else if (no_A && !no_sv) {
    // SWTN with fixed threshold b
    return dswtn(t_adj, b, mu_drift, sv, s, lambda, c, log_out);
  } else {
    // Full model: integrate dswtn over threshold ~ Unif(b-A, b)
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl = get_gl_nodes_weights(n_nodes);
    const Rcpp::NumericVector gl_nodes   = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
    if (emc2_isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;
    // Map nodes from [-1,1] to [b-A, b]
    const double center     = b - 0.5 * A;
    const double half_width = 0.5 * A;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double thresh_j = center + half_width * gl_nodes[j];
      integral += gl_weights[j] * dswtn_core(
        t_adj, thresh_j, mu_drift, sv, s, 0.0, c, log_prob_gt_c, false
      );
    }
    // Scale: GL on [-1,1]→[b-A,b] gives factor half_width; divide by A (uniform density)
    const double out_val_unkilled = integral * half_width / A;   // = integral * 0.5
    const double out_val = out_val_unkilled * std::exp(-lambda * t_adj);
    if (out_val < 0.0 || ISNAN(out_val)) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(out_val_unkilled) - lambda * t_adj : out_val;
  }
}

// [[Rcpp::export]]
double prdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s = 1.0, double c = 0.0, double lambda = 0.0,
                int n_gauss_nodes = 20, bool log_out = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? 0.0 : 1.0;

  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv < 1e-7);

  if (!emc2_isfinite(t_adj) && lambda > 1e-10) {
    return prdmswtn_killed_inf_quad(b, mu_drift, A, sv, s, c, lambda, n_gauss_nodes, log_out);
  }

  if (no_sv) {
    // sv=0: standard Wald with start-point variability A (pwald handles A=0 internally)
    if (mu_drift <= 1e-7) return log_out ? R_NegInf : 0.0;
    return pwald(t_adj, b, mu_drift, s, A, lambda, log_out);
  } else if (no_A && !no_sv) {
    return pswtn(t_adj, b, mu_drift, sv, s, c, lambda, log_out);
  } else {
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl = get_gl_nodes_weights(n_nodes);
    const Rcpp::NumericVector gl_nodes   = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
    if (!(alpha < 1.0)) return log_out ? R_NegInf : 0.0;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double u = 0.5 * (gl_nodes[j] + 1.0);
      const double p = alpha + (1.0 - alpha) * u;
      const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
      integral += gl_weights[j] * pwald(t_adj, b, drift_j, s, A, lambda, false);
    }
    double out_val = 0.5 * integral;
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
                       NumericVector s = 1.0, NumericVector c = 0.0, NumericVector lambda = 0.0,
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
                      pick(sv, i), pick(s, i), pick(lambda, i), pick(c, i), n_gauss_nodes, log_out);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector t0, NumericVector sv,
                       NumericVector s = 1.0, NumericVector c = 0.0, NumericVector lambda = 0.0,
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
                      pick(sv, i), pick(s, i), pick(c, i), pick(lambda, i),n_gauss_nodes, log_out);
  }
  return cdf;
}

// --------------------------------------------------------------------------
// GBM vectorized wrappers.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector dGBMspv(NumericVector t, NumericVector v, NumericVector b,
                      NumericVector A, NumericVector t0, NumericVector s = 1.0,
                      bool log_out = false) {
  const int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { pdf[i] = log_out ? R_NegInf : 0.0; continue; }
    pdf[i] = dgbm(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i), log_out);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pGBMspv(NumericVector t, NumericVector v, NumericVector b,
                      NumericVector A, NumericVector t0, NumericVector s = 1.0,
                      bool log_out = false) {
  const int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { cdf[i] = log_out ? R_NegInf : 0.0; continue; }
    cdf[i] = pgbm(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i), log_out);
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

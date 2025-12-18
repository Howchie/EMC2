#ifndef MODEL_SS_H
#define MODEL_SS_H

#include <Rcpp.h>
#include <cmath>

inline double ss_clamp_prob01(double p) {
  if (!std::isfinite(p)) return p;
  if (p < 0.0) return 0.0;
  if (p > 1.0) return 1.0;
  return p;
}

inline double ss_exgauss_cdf(double x, double mu, double sigma, double tau) {
  // Ex-Gaussian CDF with tau as the mean of the exponential component.
  if (!(tau > 0.0) || !(sigma > 0.0)) return NA_REAL;

  if (!R_finite(x)) {
    return (x < 0.0) ? 0.0 : 1.0;
  }

  // Shifted exponential approximation when sigma is tiny.
  if (sigma < 1e-4) {
    return R::pexp(x - mu, 1.0 / tau, /*lower_tail=*/1, /*log_p=*/0);
  }

  // Normal approximation when tau is small relative to sigma.
  if (tau < 0.05 * sigma) {
    return R::pnorm(x, mu, sigma, /*lower_tail=*/1, /*log_p=*/0);
  }

  const double s2 = sigma * sigma;
  const double z = x - mu - (s2 / tau);

  const double term1 = R::pnorm((x - mu) / sigma, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
  const double logPhi = R::pnorm(z / sigma, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);

  const double a = mu + (s2 / tau);
  const double expo = (a * a - mu * mu - 2.0 * x * (s2 / tau)) / (2.0 * s2);
  const double term2 = std::exp(logPhi + expo);

  return ss_clamp_prob01(term1 - term2);
}

inline double ss_exgauss_pdf(double x, double mu, double sigma, double tau) {
  // Ex-Gaussian PDF with tau as the mean of the exponential component.
  if (!(tau > 0.0) || !(sigma > 0.0)) return NA_REAL;

  if (!R_finite(x)) return 0.0;

  if (sigma < 1e-4) {
    return R::dexp(x - mu, 1.0 / tau, /*give_log=*/0);
  }

  if (tau < 0.05 * sigma) {
    return R::dnorm(x, mu, sigma, /*give_log=*/0);
  }

  const double s2 = sigma * sigma;
  const double z = x - mu - (s2 / tau);

  const double logPhi = R::pnorm(z / sigma, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
  const double logpdf = logPhi - std::log(tau) - (z + (s2 / (2.0 * tau))) / tau;

  const double pdf = std::exp(logpdf);
  if (!(pdf >= 0.0) || !std::isfinite(pdf)) return 0.0;
  return pdf;
}

inline double ss_exgauss_surv(double x, double mu, double sigma, double tau) {
  const double cdf = ss_exgauss_cdf(x, mu, sigma, tau);
  if (!std::isfinite(cdf)) return NA_REAL;
  return ss_clamp_prob01(1.0 - cdf);
}

inline double ss_stop_cdf_abs(double t_abs, double SSD, double muS, double sigmaS, double tauS) {
  return ss_exgauss_cdf(t_abs - SSD, muS, sigmaS, tauS);
}

inline double ss_stop_surv_abs(double t_abs, double SSD, double muS, double sigmaS, double tauS) {
  return ss_exgauss_surv(t_abs - SSD, muS, sigmaS, tauS);
}

// Scalar adapters for using ex-Gaussian go racers with existing integration helpers.
// Expects par = {mu, sigma, tau}.
inline double ss_exg_go_pdf1(double t, const double* par, void* /*ctx*/) {
  return ss_exgauss_pdf(t, par[0], par[1], par[2]);
}

inline double ss_exg_go_cdf1(double t, const double* par, void* /*ctx*/) {
  return ss_exgauss_cdf(t, par[0], par[1], par[2]);
}

#endif


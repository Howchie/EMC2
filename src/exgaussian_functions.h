#ifndef exgaussian_functions_h
#define exgaussian_functions_h

#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include "composite_functions.h"
#include "wald_functions.h"
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

const double SIG_TAU_EPS = 1e-12;

// probability density function of ex-Gaussian distribution
double dexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {

  // minimal parameter validation - assuming all inputs are finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // Numerically stable branch for extreme tails where z = (x-mu)/sigma - sigma/tau is very negative.
  {
    double y = (x - mu) / sig_p;
    double a = sig_p / tau_p;
    double z = y - a;
    if (z < -8.0) {
      // Asymptotic log-density via the Mills-ratio expansion of log Phi(z) as
      // z -> -inf. The tau-dependent pieces of the exact log-density cancel
      // exactly, leaving the Gaussian core plus the log-Mills series:
      //   log Phi(z) = log phi(z) - log(-z) + log(1 - 1/z^2 + 3/z^4 - ...)
      //              = -0.5 log(2pi) - 0.5 z^2 - log(-z) - 1/z^2 + 2.5/z^4 + ...
      // so  log f ~ -log(tau) - 0.5 log(2pi) - 0.5 y^2 - log(-z) - 1/z^2 + 2.5/z^4.
      // Keeping the -1/z^2 + 2.5/z^4 corrections (was: leading term only) makes
      // this branch match the exact density to ~4e-6 at the z=-8 switch (was
      // ~1.5%), so the integrate/GL routes agree with the exact analytic
      // stop-success form and the branch discontinuity is negligible for the
      // sampler. (Cheap; no extra special functions.)
      double zi2 = 1.0 / (z * z);
      double log_out_stable = -std::log(tau_p) - 0.5 * std::log(2.0 * M_PI)
                              - 0.5 * y * y - std::log(-z)
                              - zi2 + 2.5 * zi2 * zi2;
      return log_d ? log_out_stable : std::exp(log_out_stable);
    }
  }

  // compute Phi term
  double z = (x - mu) / sig_p - sig_p / tau_p;
  double log_phi = pnorm_std(z, true, true);
  if (emc2_isnan(log_phi)) {
    return NA_REAL;
  }
  if (emc2_isinf(log_phi)) {
    return log_d ? log_phi : std::exp(log_phi);
  }

  // compute exp term
  double log_exp = (mu - x) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // final output: log density of ex-Gaussian
  double log_out = -std::log(tau_p) + log_exp + log_phi;

  return log_d ? log_out : std::exp(log_out);
}

// Direct upper-tail log probability of the ex-Gaussian:
//   S(q) = Q(u) + exp(k) * Phi(u - sigma/tau),
//   u = (q - mu)/sigma, k = (mu - q)/tau + sigma^2/(2 tau^2).
// Both terms are non-negative, so log S is a plain log_sum_exp of two log
// terms and never round-trips through the lower-tail CDF (no exp, no
// rounding to one, no log1m).
inline double pexg_log_upper(double q, double mu, double sigma, double tau) {
  const double tau_p = std::max(tau, SIG_TAU_EPS);
  const double sig_p = std::max(sigma, SIG_TAU_EPS);
  const double u = (q - mu) / sig_p;
  const double log_q_term = pnorm_log_direct(u, false);
  const double log_exp_term = (mu - q) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);
  const double log_phi2 = pnorm_log_direct(u - sig_p / tau_p, true);
  return std::fmin(0.0, log_sum_exp(log_q_term, log_exp_term + log_phi2));
}

// cumulative distribution function of ex-Gaussian distribution
double pexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {

  // minimal parameter validation - assuming mu, sigma, and tau are all finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // handle infinite q
  if (emc2_isinf(q)) {
    double out = (q < 0.) ? 0. : 1.;
    if (!lower_tail) out = 1. - out;
    return log_p ? std::log(out) : out;
  }

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // compute the two Phi terms
  double log_phi_1 = pnorm_std((q - mu) / sig_p, true, true);
  double log_phi_2 = pnorm_std(
    (q - mu) / sig_p - sig_p / tau_p, true, true
  );

  // compute the exp term in log space
  double log_exp_term = (mu - q) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // combined second term
  double log_second_term = log_exp_term + log_phi_2;

  // now obtain ex-Gaussian log CDF
  double log_cdf_lower;
  if (log_phi_1 > log_second_term) {
    log_cdf_lower = log_diff_exp(log_phi_1, log_second_term);
  } else {
    log_cdf_lower = R_NegInf;
  }

  double out;
  if (lower_tail) {
    out = log_p ? log_cdf_lower : std::exp(log_cdf_lower);
  } else {
    if (log_cdf_lower == R_NegInf) {
      out = log_p ? 0. : 1.;
    } else {
      double cdf_lower = std::exp(log_cdf_lower);
      if (cdf_lower >= 1. - 1e-8) {
        // Near saturation 1 - cdf has lost its digits; switch to the direct
        // upper-tail expression instead of rounding the survivor to zero.
        const double log_surv = pexg_log_upper(q, mu, sig_p, tau_p);
        out = log_p ? log_surv : std::exp(log_surv);
      } else {
        out = log_p ? log1m(cdf_lower) : -std::expm1(log_cdf_lower);
      }
    }
  }

  return(out);
}

// probability density function of truncated ex-Gaussian distribution
double dtexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = R_NegInf,
    const double upper = R_PosInf,
    const bool log_d = false
) {

  if (lower == R_NegInf && upper == R_PosInf) {
    return dexg(x, mu, sigma, tau, log_d);
  }
  if (sigma <= 0. || tau <= 0.) return NA_REAL;
  if (lower >= upper) return NA_REAL;
  if (x <= lower || x >= upper) return log_d ? R_NegInf : 0.;

  double x_ld = dexg(x, mu, sigma, tau, true);
  if (x_ld == R_NegInf) return log_d ? R_NegInf : 0.;

  double lower_lcdf, upper_lcdf;
  if (lower == R_NegInf) {
    lower_lcdf = R_NegInf;
  } else {
    lower_lcdf = pexg(lower, mu, sigma, tau, true, true);
  }
  if (upper == R_PosInf) {
    upper_lcdf = 0.;
  } else {
    upper_lcdf = pexg(upper, mu, sigma, tau, true, true);
  }

  double log_normaliser;
  if (lower_lcdf == R_NegInf) {
    log_normaliser = upper_lcdf;
  } else if (lower_lcdf > -M_LN2) {
    // Both bounds sit in the upper half: the lower-tail log CDFs are nearly
    // equal (both ~0) and their log difference cancels.  The same normaliser
    // as a survivor difference S(lower) - S(upper) keeps its digits.
    const double ls_lo = pexg_log_upper(lower, mu, sigma, tau);
    const double ls_up = (upper == R_PosInf) ? R_NegInf
                                             : pexg_log_upper(upper, mu, sigma, tau);
    log_normaliser = log_diff_exp(ls_lo, ls_up);
  } else {
    log_normaliser = log_diff_exp(upper_lcdf, lower_lcdf);
  }
  if (!(log_normaliser > R_NegInf) || ISNAN(log_normaliser)) {
    return log_d ? R_NegInf : 0.;
  }

  double log_out = x_ld - log_normaliser;
  return log_d ? log_out : std::exp(log_out);
}

// cumulative distribution function of truncated ex-Gaussian distribution
double ptexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = R_NegInf,
    const double upper = R_PosInf,
    const bool lower_tail = true,
    const bool log_p = false
) {

  if (lower == R_NegInf && upper == R_PosInf) {
    return pexg(q, mu, sigma, tau, lower_tail, log_p);
  }
  if (sigma <= 0. || tau <= 0.) return NA_REAL;
  if (lower >= upper) return NA_REAL;
  if (q <= lower) {
    double out = lower_tail ? 0. : 1.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }
  if (q >= upper) {
    double out = lower_tail ? 1. : 0.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }

  double q_cdf = pexg(q, mu, sigma, tau);

  double lower_cdf, upper_cdf;
  if (lower == R_NegInf) {
    lower_cdf = 0.;
  } else {
    lower_cdf = pexg(lower, mu, sigma, tau);
  }
  if (upper == R_PosInf) {
    upper_cdf = 1.;
  } else {
    upper_cdf = pexg(upper, mu, sigma, tau);
  }

  // Natural fast path: ordinary bounds where the subtractions keep their
  // digits.  Otherwise fall back to matched log differences below.
  const double normaliser = upper_cdf - lower_cdf;
  if (normaliser > 1e-10 * std::max(upper_cdf, 1e-300)) {
    double out;
    if (lower_tail) {
      out = (q_cdf - lower_cdf) / normaliser;
    } else {
      out = (upper_cdf - q_cdf) / normaliser;
    }
    out = std::max(0., std::min(1., out));
    // A zero here can only be numerator cancellation (q is strictly inside
    // the bounds), which matters when the caller wants a log tail.
    if (out > 0. || !log_p) {
      return log_p ? std::log(out) : out;
    }
  }

  // Log-space fallback: form numerator and normaliser as log differences on
  // the better-conditioned side (lower-tail logs on the left, survivor logs
  // on the right of the median).
  auto lcdf = [&](double x) {
    return pexg(x, mu, sigma, tau, true, true);
  };
  auto lsurv = [&](double x) {
    return (x == R_PosInf) ? R_NegInf : pexg_log_upper(x, mu, sigma, tau);
  };
  const double lower_lcdf = (lower == R_NegInf) ? R_NegInf : lcdf(lower);

  double log_norm;
  if (lower_lcdf == R_NegInf) {
    log_norm = (upper == R_PosInf) ? 0.0 : lcdf(upper);
  } else if (lower_lcdf > -M_LN2) {
    log_norm = log_diff_exp(lsurv(lower), lsurv(upper));
  } else {
    log_norm = log_diff_exp((upper == R_PosInf) ? 0.0 : lcdf(upper), lower_lcdf);
  }
  if (!(log_norm > R_NegInf) || ISNAN(log_norm)) return NA_REAL;

  double log_num;
  if (lower_tail) {
    // F(q) - F(lower) = S(lower) - S(q)
    if (lower_lcdf == R_NegInf) {
      log_num = lcdf(q);
    } else if (lower_lcdf > -M_LN2) {
      log_num = log_diff_exp(lsurv(lower), lsurv(q));
    } else {
      log_num = log_diff_exp(lcdf(q), lower_lcdf);
    }
  } else {
    // F(upper) - F(q) = S(q) - S(upper)
    const double q_lcdf = lcdf(q);
    if (q_lcdf > -M_LN2) {
      log_num = log_diff_exp(lsurv(q), lsurv(upper));
    } else {
      log_num = log_diff_exp((upper == R_PosInf) ? 0.0 : lcdf(upper), q_lcdf);
    }
  }
  if (ISNAN(log_num)) log_num = R_NegInf;
  double log_out = log_num - log_norm;
  if (ISNAN(log_out) || log_out == R_NegInf) return log_p ? R_NegInf : 0.;
  if (log_out > 0.) log_out = 0.;
  return log_p ? log_out : std::exp(log_out);
}

// Scalar versions for GSL integration
inline double dtexg_scalar(double x, const double* par, void* /*ctx*/) {
  // par: mu, sigma, tau, ..., lb=8
  return dtexg(x, par[0], par[1], par[2], par[8], R_PosInf, false);
}

inline double ptexg_scalar(double x, const double* par, void* /*ctx*/) {
  // par: mu, sigma, tau, ..., lb=8
  return ptexg(x, par[0], par[1], par[2], par[8], R_PosInf, true, false);
}

inline double dexg_scalar(double x, const double* par, void* /*ctx*/) {
  const double t0 = par[3];
  if (!R_FINITE(t0) || t0 < 0.0) return 0.0;
  const double tt = x - t0;
  if (tt <= 0.0) return 0.0;
  return dtexg(tt, par[0], par[1], par[2], 0.0, R_PosInf, false);
}

inline double pexg_scalar(double x, const double* par, void* /*ctx*/) {
  const double t0 = par[3];
  if (!R_FINITE(t0) || t0 < 0.0) return 0.0;
  const double tt = x - t0;
  if (tt <= 0.0) return 0.0;
  return ptexg(tt, par[0], par[1], par[2], 0.0, R_PosInf, true, false);
}

#endif

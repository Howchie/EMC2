#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include "utility_functions.h"
#include "composite_functions.h"

using namespace Rcpp;

constexpr double L_PI = 1.1447298858494001741434; 
constexpr double LOG_2PI = 1.83787706640934548356; 
constexpr double FAST_NORM_RT2PI = 2.506628274631000502415765284811;
constexpr double FAST_NORM_SPLIT = 7.07106781186547;
constexpr double FAST_NORM_N0 = 220.206867912376;
constexpr double FAST_NORM_N1 = 221.213596169931;
constexpr double FAST_NORM_N2 = 112.079291497871;
constexpr double FAST_NORM_N3 = 33.912866078383;
constexpr double FAST_NORM_N4 = 6.37396220353165;
constexpr double FAST_NORM_N5 = 0.700383064443688;
constexpr double FAST_NORM_N6 = 3.52624965998911e-02;
constexpr double FAST_NORM_M0 = 440.413735824752;
constexpr double FAST_NORM_M1 = 793.826512519948;
constexpr double FAST_NORM_M2 = 637.333633378831;
constexpr double FAST_NORM_M3 = 296.564248779674;
constexpr double FAST_NORM_M4 = 86.7807322029461;
constexpr double FAST_NORM_M5 = 16.064177579207;
constexpr double FAST_NORM_M6 = 1.75566716318264;
constexpr double FAST_NORM_M7 = 8.83883476483184e-02;
constexpr double LOG_SQRT_2PI = 0.91893853320467274178;

inline double fast_norm_phi(double x) {
  const double z = std::fabs(x);
  double c = 0.0;

  if (z <= 37.0) {
    const double e = std::exp(-z * z / 2.0);
    if (z < FAST_NORM_SPLIT) {
      const double n = (((((FAST_NORM_N6 * z + FAST_NORM_N5) * z + FAST_NORM_N4) * z + FAST_NORM_N3) * z + FAST_NORM_N2) * z + FAST_NORM_N1) * z + FAST_NORM_N0;
      const double d = ((((((FAST_NORM_M7 * z + FAST_NORM_M6) * z + FAST_NORM_M5) * z + FAST_NORM_M4) * z + FAST_NORM_M3) * z + FAST_NORM_M2) * z + FAST_NORM_M1) * z + FAST_NORM_M0;
      c = e * n / d;
    } else {
      const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
      c = e / (FAST_NORM_RT2PI * f);
    }
  }

  return x <= 0.0 ? c : 1.0 - c;
}

// log(P(Z > z)) for z >= 0, computed directly without materialising the probability.
inline double fast_log_upper_tail(double z) {
  if (z > 37.0) return R_NegInf;
  const double z2 = z * z;
  if (z >= FAST_NORM_SPLIT) {
    const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
    return -0.5 * z2 - LOG_SQRT_2PI - std::log(f);
  }
  const double n = (((((FAST_NORM_N6 * z + FAST_NORM_N5) * z + FAST_NORM_N4) * z + FAST_NORM_N3) * z + FAST_NORM_N2) * z + FAST_NORM_N1) * z + FAST_NORM_N0;
  const double d = ((((((FAST_NORM_M7 * z + FAST_NORM_M6) * z + FAST_NORM_M5) * z + FAST_NORM_M4) * z + FAST_NORM_M3) * z + FAST_NORM_M2) * z + FAST_NORM_M1) * z + FAST_NORM_M0;
  return -0.5 * z2 + std::log(n) - std::log(d);
}

inline double pnorm_std(double x, bool lower = true, bool log_p = false) {
#ifdef USE_FAST_PNORM
  if (log_p) {
    const bool want_tail = (x >= 0.0) == (!lower);
    if (want_tail) {
      double p = fast_norm_phi(x);
      if (!lower) p = 1.0 - p;
      if (p > 1e-10) return std::log(p);
      return lower ? fast_log_upper_tail(-x) : fast_log_upper_tail(x);
    }
    // Bulk side: fast_norm_phi(-|x|) gives the small lower-tail probability (≤ 0.5),
    // so log1p is stable — no cancellation.
    return std::log1p(-fast_norm_phi(-std::fabs(x)));
  }
  double p = fast_norm_phi(x);
  if (!lower) p = 1.0 - p;
  return p;
#else
  return R::pnorm(x, 0.0, 1.0, lower, log_p);
#endif
}

// Fast log-normal CDF/PDF wrappers.
// For sdlog <= 0 or non-finite sdlog, fall back to R's implementation to preserve semantics.
inline double plnorm_std(double x, double meanlog, double sdlog,
                         bool lower_tail = true, bool log_p = false) {
  if (x <= 0.0) {
    if (log_p) return lower_tail ? R_NegInf : 0.0;
    return lower_tail ? 0.0 : 1.0;
  }
  if (!emc2_isfinite(sdlog) || !(sdlog > 0.0)) {
    return R::plnorm(x, meanlog, sdlog, lower_tail, log_p);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return pnorm_std(z, lower_tail, log_p);
}

inline double dlnorm_std(double x, double meanlog, double sdlog,
                         bool log_p = false) {
  if (x <= 0.0) return log_p ? R_NegInf : 0.0;
  if (!emc2_isfinite(sdlog) || !(sdlog > 0.0)) {
    return R::dlnorm(x, meanlog, sdlog, log_p);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  const double log_pdf = -std::log(x) - std::log(sdlog) - 0.5 * z * z - LOG_SQRT_2PI;
  return log_p ? log_pdf : std::exp(log_pdf);
}

inline double lnorm_log_surv_std(double x, double meanlog, double sdlog) {
  if (x <= 0.0) return 0.0;
  if (!emc2_isfinite(sdlog) || !(sdlog > 0.0)) {
    const double cdf = R::plnorm(x, meanlog, sdlog, true, false);
    if (cdf >= 1.0) return R_NegInf;
    if (cdf <= 0.0) return 0.0;
    return std::log1p(-cdf);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return pnorm_std(z, false, true);
}

inline double pigt0(double t, double k, double l){
  if (l == 0.0) return 2.0 * pnorm_std(-k / std::sqrt(t), true, false);
  double mu = k / l;
  double lambda = k * k;

  double p1 = pnorm_std(std::sqrt(lambda/t) * (1. + t/mu), false, true);
  double p2 = pnorm_std(std::sqrt(lambda/t) * (1. - t/mu), false, false);

  return std::exp(2.0 * k * l + p1) + p2;
}

inline double digt0(double t, double k, double l){
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

inline double pigt_impl(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
  if (t <= 0.) return 0.;
  if (a < threshold) return pigt0(t, k, l);

  const double sqt = std::sqrt(t);

  if (std::abs(l) < threshold) {
    const double p1 = pnorm_std((k + a) / sqt, true, false);
    const double p2 = pnorm_std((k - a) / sqt, true, false);
    const double diff_sq = (k + a) * (k + a) - (k - a) * (k - a);
    const double t1 = sqt * std::exp(-0.5 * (k - a) * (k - a) / t) * (-std::expm1(-0.5 * diff_sq / t)) / FAST_NORM_RT2PI;
    return (t1 + (k + a) * (1.0 - p1) - (k - a) * (1.0 - p2)) / a;
  }

  const double t1a = std::exp(- .5 * (k - a - t * l) * (k - a - t * l) / t);
  const double t1b = std::exp(- .5 * (a + k - t * l) * (a + k - t * l) / t);
  const double t1 = sqt * (t1a - t1b) / FAST_NORM_RT2PI;

  const double t2a = std::exp(2. * l * (k - a) + pnorm_std(- (k - a + t * l) / sqt, true, true));
  const double t2b = std::exp(2. * l * (k + a) + pnorm_std(- (k + a + t * l) / sqt, true, true));
  const double t2 = a + (t2b - t2a) / (2. * l);

  const double t4a = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
  const double t4b = 2. * pnorm_std((k - a) / sqt - sqt * l, true, false) - 1.;
  const double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;

  return .5 * (t4 + t2 + t1) / a;
}

inline double digt_impl(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10){
  if (t <= 0.) return 0.;
  if (a < threshold) return digt0(t, k, l);

  if (std::abs(l) < threshold) {
    const double diff_sq = (k + a) * (k + a) - (k - a) * (k - a);
    const double log_term = -0.5 * (k - a) * (k - a) / t + std::log(-std::expm1(-0.5 * diff_sq / t));
    return std::exp(-0.5 * (LOG_2PI + std::log(t)) + log_term - std::log(2.0 * a));
  }

  const double sqt = std::sqrt(t);

  const double t1a = - (a - k + t * l) * (a - k + t * l) / (2. * t);
  const double t1b = - (a + k - t * l) * (a + k - t * l) / (2. * t);
  const double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

  const double t2a = 2. * pnorm_std((- k + a) / sqt + sqt * l, true, false) - 1.;
  const double t2b = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
  const double t2 = 0.5 * l * (t2a + t2b);

  return (t1 + t2) / (2.0 * a);
}

// Global symbols for exported functions (defined with [[Rcpp::export]] in wald_functions.cpp)
double pigt(double t, double k, double l, double a, double threshold);
double digt(double t, double k, double l, double a, double threshold);

// --------------------------------------------------------------------------
// Wald SPV CDF/PDF without kill, sigma = 1 (caller pre-scales).
// Parameterised in canonical Wald form:
//   start x ~ U[0, A], threshold b, so distance d = b - x ~ U[b-A, b].
// The CDF is the closed-form integral of the point-Wald CDF over d.
// --------------------------------------------------------------------------
inline double pwald_k0(double t, double b, double mu, double A) {
  if (t <= 0.0 || b <= 0.0) return 0.0;
  if (A <= 1e-12) return pigt0(t, b, mu);
  if (A <= 0.0) return 0.0;

  const double d_raw_lo = b - A;
  const double d_lo = std::fmax(0.0, d_raw_lo);
  const double d_hi = b;
  const double immediate = (d_raw_lo < 0.0) ? std::fmin(-d_raw_lo, A) : 0.0;
  if (d_hi <= d_lo) return std::fmax(0.0, std::fmin(1.0, immediate / A));

  const double sqt = std::sqrt(t);
  auto phi_std = [](double z) {
    return std::exp(-0.5 * z * z) / FAST_NORM_RT2PI;
  };

  if (std::abs(mu) <= 1e-12) {
    auto antideriv_mu0 = [&](double d) {
      const double z = -d / sqt;
      return 2.0 * (d * pnorm_std(z, true, false) - sqt * phi_std(z));
    };
    const double p = (immediate + antideriv_mu0(d_hi) - antideriv_mu0(d_lo)) / A;
    return std::fmax(0.0, std::fmin(1.0, p));
  }

  auto antideriv = [&](double d) {
    const double z1 = (mu * t - d) / sqt;
    const double z2 = -(d + mu * t) / sqt;
    const double term1 = (d - mu * t) * pnorm_std(z1, true, false) - sqt * phi_std(z1);
    const double term2 = (std::exp(2.0 * mu * d + pnorm_std(z2, true, true)) -
                          pnorm_std(z1, true, false)) / (2.0 * mu);
    return term1 + term2;
  };

  const double p = (immediate + antideriv(d_hi) - antideriv(d_lo)) / A;
  return std::fmax(0.0, std::fmin(1.0, p));
}

inline double dwald_k0(double t, double b, double mu, double A) {
  if (t <= 0.0 || b <= 0.0) return 0.0;
  const double b_lo  = b - A;
  if (A <= 1e-12) return digt0(t, b, mu);
  const double sqt   = std::sqrt(t);
  const double inv_A = 1.0 / A;

  if (std::abs(mu) <= 1e-12) {
    const double t1 = M_SQRT1_2 * (std::exp(-0.5 * b_lo * b_lo / t)
                                  -std::exp(-0.5 * b    * b    / t))
                    / (std::sqrt(M_PI) * sqt);
    return t1 * inv_A;
  }

  const double t1a = -(b_lo - t * mu) * (b_lo - t * mu) / (2.0 * t);
  const double t1b = -(b    - t * mu) * (b    - t * mu) / (2.0 * t);
  const double t1  = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

  const double t2a = 2.0 * pnorm_std(-b_lo / sqt + sqt * mu, true, false) - 1.0;
  const double t2b = 2.0 * pnorm_std( b    / sqt - sqt * mu, true, false) - 1.0;
  const double t2  = 0.5 * mu * (t2a + t2b);

  return (t1 + t2) * inv_A;
}

#endif

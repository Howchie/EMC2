#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)

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

inline double pnorm_std(double x, bool lower = true, bool log_p = false) {
#ifdef USE_FAST_PNORM
  double p = fast_norm_phi(x);
  if (!lower) p = 1.0 - p;
  if (log_p) return (p <= 0.0) ? R_NegInf : std::log(p);
  return p;
#else
  return R::pnorm(x, 0.0, 1.0, lower, log_p);
#endif
}

double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double p1 = pnorm_std(std::sqrt(lambda/t) * (1. + t/mu), false, false);
  double p2 = pnorm_std(std::sqrt(lambda/t) * (1. - t/mu), false, false);

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
    double t5a = 2. * pnorm_std((k + a) / sqt, true, false) - 1;
    double t5b = 2. * pnorm_std((- k - a) / sqt, true, false) - 1;

    double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

    cdf = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
  } else {
    double t1a = std::exp(- .5 * (k - a - t * l) * (k - a - t * l) / t);
    double t1b = std::exp(- .5 * (a + k - t * l) * (a + k - t * l) / t);
    double t1 = std::exp(.5* (lgt - M_LN2 - L_PI)) * (t1a - t1b);

    double t2a = std::exp(2. * l * (k - a) + pnorm_std(- (k - a + t * l) / sqt, true, true));
    double t2b = std::exp(2. * l * (k + a) + pnorm_std(- (k + a + t * l) / sqt, true, true));
    double t2 = a + (t2b - t2a) / (2. * l);

    double t4a = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
    double t4b = 2. * pnorm_std((k - a) / sqt - sqt * l, true, false) - 1.;
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

    double t1a = - (a - k + t * l) * (a - k + t * l) / (2. * t);
    double t1b = - (a + k - t * l) * (a + k - t * l) / (2. * t);
    double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

    double t2a = 2. * pnorm_std((- k + a) / sqt + sqt * l, true, false) - 1.;
    double t2b = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
    double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);

    pdf = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
  }
  if (pdf < 0. || std::isnan(pdf)) {
    return 0.;
  }
  return pdf;
}

#endif

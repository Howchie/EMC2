#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include "utility_functions.h"
#include "composite_functions.h"

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

inline double pnorm_std(double x, bool lower = true, bool log_p = false) {
#ifdef USE_FAST_PNORM
  // In log-space, use fast approximation in the central region and defer to
  // R's stable tail implementation for extreme |x| values where cancellation
  // can dominate.
  if (log_p) {
    if (x <= -8.0 || x >= 8.0) return R::pnorm(x, 0.0, 1.0, lower, true);
    double p = fast_norm_phi(x);
    if (!lower) p = 1.0 - p;
    return (p <= 0.0) ? R_NegInf : std::log(p);
  }
  double p = fast_norm_phi(x);
  if (!lower) p = 1.0 - p;
  return p;
#else
  return R::pnorm(x, 0.0, 1.0, lower, log_p);
#endif
}

// Fast log-normal CDF/PDF wrappers.
struct SignedLogTerm {
  int sign;       // -1, 0, +1
  double logabs;  // log(|x|), or -Inf for zero
};

inline SignedLogTerm signed_log_zero() {
  return {0, R_NegInf};
}

inline SignedLogTerm signed_log_from_value(double x) {
  if (!R_FINITE(x) || x == 0.0) return signed_log_zero();
  return {x < 0.0 ? -1 : 1, std::log(std::fabs(x))};
}

inline SignedLogTerm signed_log_from_logabs(double logabs, int sign = 1) {
  if (sign == 0 || logabs == R_NegInf) return signed_log_zero();
  return {sign < 0 ? -1 : 1, logabs};
}

inline SignedLogTerm signed_log_mul(SignedLogTerm a, SignedLogTerm b) {
  if (a.sign == 0 || b.sign == 0) return signed_log_zero();
  if (!R_FINITE(a.logabs) || !R_FINITE(b.logabs)) return signed_log_zero();
  return {a.sign * b.sign, a.logabs + b.logabs};
}

inline SignedLogTerm signed_log_scale(SignedLogTerm a, double scale) {
  if (a.sign == 0 || scale == 0.0) return signed_log_zero();
  if (!R_FINITE(scale)) return signed_log_zero();
  return {a.sign * (scale < 0.0 ? -1 : 1), a.logabs + std::log(std::fabs(scale))};
}

inline SignedLogTerm signed_log_scale_logabs(SignedLogTerm a, double log_scale, int scale_sign = 1) {
  if (a.sign == 0 || scale_sign == 0 || !R_FINITE(log_scale)) return signed_log_zero();
  return {a.sign * (scale_sign < 0 ? -1 : 1), a.logabs + log_scale};
}

inline SignedLogTerm signed_log_add(SignedLogTerm a, SignedLogTerm b) {
  if (a.sign == 0) return b;
  if (b.sign == 0) return a;

  if (a.sign == b.sign) {
    return {a.sign, log_sum_exp(a.logabs, b.logabs)};
  }

  if (a.logabs == b.logabs) return signed_log_zero();
  if (a.logabs > b.logabs) {
    return {a.sign, log_diff_exp(a.logabs, b.logabs)};
  }
  return {b.sign, log_diff_exp(b.logabs, a.logabs)};
}

inline SignedLogTerm signed_log_diff_exp(double log_a, double log_b) {
  if (!R_FINITE(log_a) || !R_FINITE(log_b)) return signed_log_zero();
  if (log_a == log_b) return signed_log_zero();
  if (log_a > log_b) return {1, log_diff_exp(log_a, log_b)};
  return {-1, log_diff_exp(log_b, log_a)};
}

inline SignedLogTerm signed_log_2phi_minus1(double z) {
  const double log2 = std::log(2.0);
  if (z >= 0.0) {
    double log_q = pnorm_std(z, false, true);
    if (!R_FINITE(log_q)) return signed_log_zero();
    double arg = log2 + log_q;
    if (arg > 0.0) arg = 0.0;
    return signed_log_from_logabs(log1m_exp(arg), 1);
  }

  double log_p = pnorm_std(z, true, true);
  if (!R_FINITE(log_p)) return signed_log_zero();
  double arg = log2 + log_p;
  if (arg > 0.0) arg = 0.0;
  return signed_log_from_logabs(log1m_exp(arg), -1);
}

// Cancellation threshold for the hybrid direct/log-space switch.
static constexpr double LOGSPACE_REL_THRESHOLD = 1e-15;

inline double pigt0(double t, double k, double l){
  double mu = k / l;
  double lambda = k * k;

  double p1 = pnorm_std(std::sqrt(lambda/t) * (1. + t/mu), false, false);
  double p2 = pnorm_std(std::sqrt(lambda/t) * (1. - t/mu), false, false);

  return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

inline double pigt0_log(double t, double k, double l){
  if (t <= 0.) {
    return R_NegInf;
  }
  const double mu = k / l;
  const double lambda = k * k;
  const double log_p1 = pnorm_std(std::sqrt(lambda/t) * (1. + t/mu), false, true);
  const double log_p2 = pnorm_std(std::sqrt(lambda/t) * (1. - t/mu), false, true);
  const double log_term1 = std::exp(std::log(2. * lambda) - std::log(mu)) + log_p1;
  return log_sum_exp(log_term1, log_p2);
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

inline double digt0_log(double t, double k, double l){
  if (t <= 0.) {
    return R_NegInf;
  }
  const double lambda = k * k;
  double e;
  if (l == 0.) {
    e = -.5 * lambda / t;
  } else {
    const double mu = k / l;
    e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
  }
  return e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI);
}

// Returns the CDF value in [0, 1] if direct computation is numerically accurate,
// or -1.0 as a sentinel meaning "fall back to log-space".
// Valid CDFs are always >= 0, so -1.0 is unambiguous.
inline double pigt_try_direct(double t, double k, double l, double a, double threshold) {
  if (t <= 0.) return 0.0;
  if (a < threshold) {
    const double v = pigt0(t, k, l);
    if (v < 0. || v > 1.0 || !R_FINITE(v)) return -1.0;
    return v;
  }

  const double sqt = std::sqrt(t);
  const double lgt = std::log(t);

  if (l < threshold) {
    const double t5a = 2. * pnorm_std((k + a) / sqt, true, false) - 1;
    const double t5b = 2. * pnorm_std((- k - a) / sqt, true, false) - 1;
    const double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    const double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    const double e6a = std::exp(t6a);
    const double e6b = std::exp(t6b);
    const double cross = ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
    const double raw = 1.0 + e6a - e6b + cross;
    const double abs_sum = 1.0 + e6a + e6b + std::fabs(cross);
    if (!R_FINITE(raw) || !R_FINITE(abs_sum) || abs_sum <= 0.0 ||
        std::fabs(raw) <= LOGSPACE_REL_THRESHOLD * abs_sum) return -1.0;
    const double cdf = 0.5 * raw / a;
    if (cdf < 0. || cdf > 1.0) return -1.0;
    return cdf;
  }

  const double t1a = std::exp(- .5 * (k - a - t * l) * (k - a - t * l) / t);
  const double t1b = std::exp(- .5 * (a + k - t * l) * (a + k - t * l) / t);
  const double t1 = std::exp(.5 * (lgt - M_LN2 - L_PI)) * (t1a - t1b);

  const double t2a = std::exp(2. * l * (k - a) + pnorm_std(- (k - a + t * l) / sqt, true, true));
  const double t2b = std::exp(2. * l * (k + a) + pnorm_std(- (k + a + t * l) / sqt, true, true));
  const double t2 = a + (t2b - t2a) / (2. * l);

  const double t4a = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
  const double t4b = 2. * pnorm_std((k - a) / sqt - sqt * l, true, false) - 1.;
  const double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;

  const double raw = t4 + t2 + t1;
  const double abs_sum = std::fabs(t4) + std::fabs(t2) + std::fabs(t1);
  if (!R_FINITE(raw) || !R_FINITE(abs_sum) || abs_sum <= 0.0 ||
      std::fabs(raw) <= LOGSPACE_REL_THRESHOLD * abs_sum) return -1.0;
  const double cdf = .5 * raw / a;
  if (cdf < 0. || cdf > 1.0) return -1.0;
  return cdf;
}

// Returns the PDF value (>= 0) if direct computation is numerically accurate,
// or -1.0 as a sentinel meaning "fall back to log-space".
// Valid densities are always >= 0, so -1.0 is unambiguous.
inline double digt_try_direct(double t, double k, double l, double a, double threshold) {
  if (t <= 0.) return 0.0;
  if (a < threshold) {
    const double v = digt0(t, k, l);
    if (v < 0. || !R_FINITE(v)) return -1.0;
    return v;
  }

  if (l < threshold) {
    const double term1 = std::exp(- (k - a) * (k - a) / (2. * t));
    const double term2 = std::exp(- (k + a) * (k + a) / (2. * t));
    const double raw = term1 - term2;
    const double abs_sum = term1 + term2;
    if (!R_FINITE(raw) || !R_FINITE(abs_sum) || abs_sum <= 0.0 ||
        std::fabs(raw) <= LOGSPACE_REL_THRESHOLD * abs_sum) return -1.0;
    const double value = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(raw) - M_LN2 - std::log(a));
    if (value < 0. || !R_FINITE(value)) return -1.0;
    return value;
  }

  const double sqt = std::sqrt(t);

  const double t1a = - (a - k + t * l) * (a - k + t * l) / (2. * t);
  const double t1b = - (a + k - t * l) * (a + k - t * l) / (2. * t);
  const double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

  const double t2a = 2. * pnorm_std((- k + a) / sqt + sqt * l, true, false) - 1.;
  const double t2b = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
  const double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);

  const double raw = t1 + t2;
  const double abs_sum = std::fabs(t1) + std::fabs(t2);
  if (!R_FINITE(raw) || !R_FINITE(abs_sum) || abs_sum <= 0.0 ||
      std::fabs(raw) <= LOGSPACE_REL_THRESHOLD * abs_sum) return -1.0;
  const double value = std::exp(std::log(raw) - M_LN2 - std::log(a));
  if (value < 0. || !R_FINITE(value)) return -1.0;
  return value;
}

inline double pigt_log_internal(double t, double k, double l, double a, double threshold){
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return std::exp(pigt0_log(t, k, l));
  }

  const double sqt = std::sqrt(t);
  const double lgt = std::log(t);

  if (l < threshold){
    const SignedLogTerm t5a = signed_log_2phi_minus1((k + a) / sqt);
    const SignedLogTerm t5b = signed_log_2phi_minus1((- k - a) / sqt);

    const double log_a = std::log(a);
    const double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - log_a;
    const double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - log_a;

    SignedLogTerm out = signed_log_from_value(1.0);
    out = signed_log_add(out, signed_log_from_logabs(t6a, 1));
    out = signed_log_add(out, signed_log_from_logabs(t6b, -1));

    const SignedLogTerm c1 = signed_log_from_value((-k + a) / (2.0 * a));
    const SignedLogTerm c2 = signed_log_from_value(-(k - a) / (2.0 * a));
    out = signed_log_add(out, signed_log_mul(c1, t5a));
    out = signed_log_add(out, signed_log_mul(c2, t5b));

    if (out.sign <= 0 || !R_FINITE(out.logabs)) {
      return 0.;
    }
    const double log_cdf = out.logabs - log_a - M_LN2;
    if (log_cdf >= 0.0) return 1.0;
    const double cdf = std::exp(log_cdf);
    return (cdf < 0. || ISNAN(cdf)) ? 0. : cdf;
  }

  const double lt1a = - .5 * (k - a - t * l) * (k - a - t * l) / t;
  const double lt1b = - .5 * (a + k - t * l) * (a + k - t * l) / t;
  const SignedLogTerm t1 = signed_log_scale_logabs(
    signed_log_diff_exp(lt1a, lt1b),
    .5 * (lgt - M_LN2 - L_PI)
  );

  const double lt2a = 2. * l * (k - a) + pnorm_std(- (k - a + t * l) / sqt, true, true);
  const double lt2b = 2. * l * (k + a) + pnorm_std(- (k + a + t * l) / sqt, true, true);
  SignedLogTerm t2 = signed_log_from_value(a);
  t2 = signed_log_add(
    t2,
    signed_log_scale_logabs(
      signed_log_diff_exp(lt2b, lt2a),
      -std::log(2.0 * std::fabs(l)),
      l < 0.0 ? -1 : 1
    )
  );

  const SignedLogTerm t4a = signed_log_2phi_minus1((k + a) / sqt - sqt * l);
  const SignedLogTerm t4b = signed_log_2phi_minus1((k - a) / sqt - sqt * l);
  const SignedLogTerm c1 = signed_log_from_value(.5 * (t * l - a - k + .5 / l));
  const SignedLogTerm c2 = signed_log_from_value(.5 * (k - a - t * l - .5 / l));
  const SignedLogTerm t4 = signed_log_add(signed_log_mul(c1, t4a), signed_log_mul(c2, t4b));

  SignedLogTerm sum = signed_log_add(signed_log_add(t4, t2), t1);
  if (sum.sign <= 0 || !R_FINITE(sum.logabs)) return 0.;

  const double log_cdf = sum.logabs - M_LN2 - std::log(a);
  if (!R_FINITE(log_cdf)) return 0.;
  if (log_cdf >= 0.) return 1.;
  const double cdf = std::exp(log_cdf);
  return (cdf < 0. || ISNAN(cdf)) ? 0. : cdf;
}

inline double digt_log_internal(double t, double k, double l, double a, double threshold){
  if (t <= 0.){
    return 0.;
  }
  if (a < threshold){
    return std::exp(digt0_log(t, k, l));
  }

  const double sqt = std::sqrt(t);

  if (l < threshold){
    const SignedLogTerm term = signed_log_diff_exp(- (k - a) * (k - a) / (2. * t), - (k + a) * (k + a) / (2. * t));
    const double log_pdf = - .5 * (M_LN2 + L_PI + std::log(t)) + term.logabs - M_LN2 - std::log(a);
    if (!R_FINITE(log_pdf)) return 0.;
    return std::exp(log_pdf);
  }

  const double lt1a = - (a - k + t * l) * (a - k + t * l) / (2. * t);
  const double lt1b = - (a + k - t * l) * (a + k - t * l) / (2. * t);
  const SignedLogTerm t1 = signed_log_scale_logabs(
    signed_log_diff_exp(lt1a, lt1b),
    -.5 * (M_LN2 + L_PI + std::log(t))
  );

  const SignedLogTerm t2a = signed_log_2phi_minus1((- k + a) / sqt + sqt * l);
  const SignedLogTerm t2b = signed_log_2phi_minus1((k + a) / sqt - sqt * l);
  const SignedLogTerm t2 = signed_log_scale(signed_log_add(t2a, t2b), 0.5 * l);

  SignedLogTerm sum = signed_log_add(t1, t2);
  if (sum.sign <= 0 || !R_FINITE(sum.logabs)) return 0.;

  const double log_pdf = sum.logabs - M_LN2 - std::log(a);
  if (!R_FINITE(log_pdf)) return 0.;
  const double pdf = std::exp(log_pdf);
  return (pdf < 0. || ISNAN(pdf)) ? 0. : pdf;
}

inline double pigt_impl(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
  const double direct = pigt_try_direct(t, k, l, a, threshold);
  if (direct < 0.) return pigt_log_internal(t, k, l, a, threshold);
  return direct;
}

inline double digt_impl(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10){
  const double direct = digt_try_direct(t, k, l, a, threshold);
  if (direct < 0.) return digt_log_internal(t, k, l, a, threshold);
  return direct;
}

// Global symbols for exported functions (to be specified with [[Rcpp::export]] in .cpp)
double pigt_old(double t, double k, double l, double a, double threshold);
double digt_old(double t, double k, double l, double a, double threshold);
double pigt_log(double t, double k, double l, double a, double threshold);
double digt_log(double t, double k, double l, double a, double threshold);
double pigt(double t, double k, double l, double a, double threshold);
double digt(double t, double k, double l, double a, double threshold);

#endif

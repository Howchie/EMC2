#ifndef EMC2_BAWD_KERNEL_H
#define EMC2_BAWD_KERNEL_H

#include <cmath>
#include <cfloat>
#include <R_ext/Arith.h>

// Finite-rho scalar arithmetic for BAwD.  The exponential kernel deliberately
// does not use these helpers: keeping that path separate preserves its exact
// historical arithmetic.
inline double bawd_expm1_div(double y) {
  if (std::fabs(y) < 1e-8)
    return 1.0 + y * (0.5 + y / 6.0);
  return std::expm1(y) / y;
}

inline double bawd_pow_diff(double c, double L) {
  return L * bawd_expm1_div(c * L);
}

inline double bawd_pk_log_tau(double rho, double x) {
  return std::log1p(x / rho);
}

inline double bawd_pk_kq(double rho, double L) {
  return rho * bawd_pow_diff(1.0 - rho, L);
}

inline double bawd_pk_kr(double rho, double gamma, double L) {
  if (gamma >= 1.0 - 1e-12)
    return bawd_pk_kq(rho, L);
  if (gamma <= 1e-12)
    return rho * std::expm1(L);
  return rho * bawd_pow_diff(1.0 - rho * gamma, L);
}

inline double bawd_pk_log_wrel(double rho, double gamma, double L) {
  return rho * (1.0 - gamma) * L;
}

inline double bawd_pk_psi(double rho, double gamma, double L) {
  const double omega = 1.0 - gamma;
  const double m = rho * omega;
  const double M = m * L;
  const double A = (1.0 - rho * gamma) * L;
  if (std::fmax(std::fabs(M), std::fabs(A)) < 1e-4) {
    const double M2 = M * M;
    const double MA = M * A;
    const double A2 = A * A;
    const double M3 = M2 * M;
    const double M2A = M2 * A;
    const double MA2 = M * A2;
    const double A3 = A2 * A;
    return rho * L * M * (
      0.5 + (M + A) / 6.0 + (M2 + MA + A2) / 24.0 +
      (M3 + M2A + MA2 + A3) / 120.0);
  }
  return rho * (std::exp(A) * bawd_pow_diff(rho - 1.0, L) -
                bawd_pow_diff(1.0 - rho * gamma, L));
}

inline double bawd_pk_log_psi_prime(double rho, double gamma, double L) {
  const double omega = 1.0 - gamma;
  return std::log(rho * omega) - rho * gamma * L +
    std::log(bawd_pow_diff(rho - 1.0, L));
}

inline double bawd_pk_peak_x(double rho, double gamma, double log_ratio) {
  return rho * std::expm1(log_ratio / (rho * (1.0 - gamma)));
}

inline double bawd_pk_newton_x(double rho, double gamma, double c) {
  if (!(c > 0.0)) return 0.0;
  if (ISNAN(c) || !R_FINITE(c)) return R_PosInf;
  const double omega = 1.0 - gamma;
  if (!(omega > 0.0)) return R_PosInf;
  double x = std::sqrt(2.0 * c / omega);
  if (!(x > 0.0) || !R_FINITE(x)) x = 1.0;
  double lo = 0.0;
  double hi = x;
  while (bawd_pk_psi(rho, gamma, std::log1p(hi / rho)) < c) {
    lo = hi;
    if (hi > 1e300) return R_PosInf;
    hi *= 2.0;
  }
  if (!(hi > lo)) return hi;
  x = std::fmin(std::fmax(x, lo), hi);
  for (int it = 0; it < 100; ++it) {
    const double L = std::log1p(x / rho);
    const double f = bawd_pk_psi(rho, gamma, L) - c;
    if (f > 0.0) hi = x; else lo = x;
    const double log_deriv = bawd_pk_log_psi_prime(rho, gamma, L);
    double xn = (R_FINITE(log_deriv) && log_deriv < std::log(DBL_MAX))
      ? x - f / std::exp(log_deriv) : 0.5 * (lo + hi);
    if (!(xn > lo) || !(xn < hi) || !R_FINITE(xn)) xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) <= 1e-14 * std::fmax(1.0, xn)) return xn;
    x = xn;
  }
  return x;
}

#endif  // EMC2_BAWD_KERNEL_H

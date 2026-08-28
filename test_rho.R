library(EMC2)

src_code <- '
#include <Rcpp.h>
using namespace Rcpp;

// copy of bawd_kernel.h functions
inline double bawd_expm1_div(double y) {
  if (std::fabs(y) < 1e-4) {
    const double y2 = y * y;
    return 1.0 + y / 2.0 + y2 / 6.0 + y2 * y / 24.0 + y2 * y2 / 120.0;
  }
  return std::expm1(y) / y;
}
inline double bawd_pow_diff(double c, double L) {
  return L * bawd_expm1_div(c * L);
}
inline double bawd_pk_log_tau(double rho, double x) {
  return std::log1p(x / rho);
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

// [[Rcpp::export]]
double test_psi(double rho, double gamma, double x) {
  return bawd_pk_psi(rho, gamma, bawd_pk_log_tau(rho, x));
}
'

Rcpp::sourceCpp(code = src_code)

print(test_psi(0.5, 0.0, 1.0))

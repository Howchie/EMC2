#include "gaussian.h"

// Out-of-line definitions for the R-facing bivariate-normal CDF wrappers.
//
// The implementations (norm_cdf_2d_vfast / norm_cdf_2d / norm_cdf_2d_fast) are
// inline in gaussian.h so that header can be included from more than one
// translation unit.  These three wrappers are [[Rcpp::export]]'d, so they need a
// single out-of-line definition with external linkage for RcppExports.cpp to
// link against.

// [[Rcpp::export]]
double pbvn_tsay(double h, double k, double rho) {
  return norm_cdf_2d_vfast(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_tvpack(double h, double k, double rho) {
  return norm_cdf_2d(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_drezner(double h, double k, double rho) {
  return norm_cdf_2d_fast(h, k, rho);
}

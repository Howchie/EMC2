#include "wald_functions.h"

// [[Rcpp::export]]
double pigt_old(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10) {
    const double v = pigt_try_direct(t, k, l, a, threshold);
    return v < 0. ? 0. : v;
}

// [[Rcpp::export]]
double digt_old(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10) {
    const double v = digt_try_direct(t, k, l, a, threshold);
    return v < 0. ? 0. : v;
}

// [[Rcpp::export]]
double pigt_log(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10) {
    return pigt_log_internal(t, k, l, a, threshold);
}

// [[Rcpp::export]]
double digt_log(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10) {
    return digt_log_internal(t, k, l, a, threshold);
}

// [[Rcpp::export]]
double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10) {
    return pigt_impl(t, k, l, a, threshold);
}

// [[Rcpp::export]]
double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10) {
    return digt_impl(t, k, l, a, threshold);
}

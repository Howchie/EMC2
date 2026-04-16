#include "wald_functions.h"

// [[Rcpp::export]]
double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10) {
    return pigt_impl(t, k, l, a, threshold);
}

// [[Rcpp::export]]
double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10) {
    return digt_impl(t, k, l, a, threshold);
}

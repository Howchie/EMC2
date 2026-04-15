#include "wald_functions.h"

// [[Rcpp::export]]
double pigt_old(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10) {
    const PigtEval direct = pigt_direct_eval(t, k, l, a, threshold);
    if (direct.value < 0. || !std::isfinite(direct.value)) {
        return 0.;
    }
    return direct.value;
}

// [[Rcpp::export]]
double digt_old(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10) {
    const DigtEval direct = digt_direct_eval(t, k, l, a, threshold);
    if (direct.value < 0. || !std::isfinite(direct.value)) {
        return 0.;
    }
    return direct.value;
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

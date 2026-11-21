#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "composite_functions.h"

// Return NA when the nearest-positive-definite projection may change W or S.
// [[Rcpp::export(rng = false)]]
double bridge_diwish(const arma::mat& W, double v, const arma::mat& S) {
  const arma::uword p = W.n_rows;
  if (!p || W.n_cols != p || S.n_rows != p || S.n_cols != p ||
      !R_FINITE(v) || v <= p - 1)
    return NA_REAL;

  double min_scale = R_PosInf, max_scale = 0, half_logdet_S = 0;
  for (arma::uword j = 0; j < p; ++j) {
    const double scale = S(j, j);
    if (!R_FINITE(scale) || scale <= 0) return NA_REAL;
    min_scale = std::min(min_scale, scale);
    max_scale = std::max(max_scale, scale);
    half_logdet_S += std::log(scale);
    for (arma::uword i = 0; i < p; ++i) {
      if (!R_FINITE(W(i, j)) || (i != j && (!R_FINITE(S(i, j)) || S(i, j) != 0))) return NA_REAL;
    }
  }
  if (min_scale <= 1e-5 * max_scale ||
      arma::norm(W - W.t(), "inf") > 1e-14 * arma::norm(W, "inf"))
    return NA_REAL;

  arma::mat factor, inverse;
  if (!arma::chol(factor, W) || !arma::inv_sympd(inverse, W)) return NA_REAL;
  // The Frobenius norms bound the ratio of the largest and smallest eigenvalues.
  if (!(arma::norm(W, "fro") * arma::norm(inverse, "fro") < 1e5))
    return NA_REAL;

  double half_logdet_W = 0, trace = 0, gamma_part = 0;
  for (arma::uword j = 0; j < p; ++j) {
    half_logdet_W += std::log(factor(j, j));
    trace += S(j, j) * inverse(j, j);
    gamma_part += R::lgammafn((v - j) / 2);
  }
  const double log_denom = gamma_part + 0.5 * v * p * M_LN2 +
    0.25 * p * (p - 1) * std::log(M_PI);
  const double log_numer = v * (0.5 * half_logdet_S) -
    (v + p + 1) * half_logdet_W - 0.5 * trace;
  const double density = std::exp(log_numer - log_denom);
  return R_FINITE(density) ? density : 1e-100;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector bridge_log_add(const Rcpp::NumericVector& a,
                                  const Rcpp::NumericVector& b) {
  if (a.size() != b.size()) Rcpp::stop("Vectors must have the same length");
  Rcpp::NumericVector out(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i)
    out[i] = ISNAN(a[i]) || ISNAN(b[i]) ? NA_REAL : log_sum_exp(a[i], b[i]);
  return out;
}

// [[Rcpp::export(rng = false)]]
double bridge_update(const Rcpp::NumericVector& l1, const Rcpp::NumericVector& l2,
                     double log_s1, double log_s2_r) {
  if (ISNAN(log_s1) || ISNAN(log_s2_r)) return NA_REAL;
  double numerator = R_NegInf, denominator = R_NegInf;
  for (double x : l2) {
    if (ISNAN(x)) return NA_REAL;
    numerator = log_sum_exp(numerator, -log_sum_exp(log_s1, log_s2_r - x));
  }
  for (double x : l1) {
    if (ISNAN(x)) return NA_REAL;
    denominator = log_sum_exp(denominator, -log_sum_exp(log_s1 + x, log_s2_r));
  }
  return (numerator - std::log(static_cast<double>(l2.size()))) -
    (denominator - std::log(static_cast<double>(l1.size())));
}

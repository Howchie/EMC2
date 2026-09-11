// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utility_functions.h"

// [[Rcpp::export]]
arma::vec fast_dmvnorm(const arma::mat& x, const arma::rowvec& mean, const arma::mat& sigma) {
  int n = x.n_rows;
  int p = x.n_cols;
  if (n == 0 || p == 0) return arma::vec();
  arma::vec out(n);
  
  arma::mat R;
  if (!arma::chol(R, sigma)) {
    out.fill(R_NegInf);
    return out;
  }
  arma::mat rooti = arma::inv(arma::trimatu(R));
  double log_const = arma::sum(arma::log(rooti.diag())) - 0.5 * p * std::log(2.0 * M_PI);
  
  arma::mat x_centered = x;
  x_centered.each_row() -= mean;
  arma::mat z = x_centered * rooti;
  out = -0.5 * arma::sum(arma::square(z), 1) + log_const;
  
  return out;
}

// The factorisation `fast_dmvnorm` performs internally, hoisted so a caller that
// evaluates the same covariance many times can do it once.
//
// C15.  With more than one block, the prior density is evaluated per subject per
// component per iteration against the SAME group covariance, and each call
// factorised it again: N x blocks Cholesky decompositions of one matrix per
// iteration.
//
// It must be this function and not R's `chol` + `backsolve`, because the point
// is that `fast_dmvnorm_rooti()` with this factor reproduces `fast_dmvnorm()`
// bit for bit -- a different triangular inverse would agree to about a ULP,
// which is not the same thing.
//
// Singular handling is preserved rather than improved: `fast_dmvnorm` returns
// -Inf for every row when the Cholesky fails, so `ok = FALSE` here means the
// caller must do the same. The regularisation ladder in `.chol_factor()` is
// PROPOSAL regularisation and deliberately does not appear here; applying it
// would quietly change the prior.
// [[Rcpp::export]]
Rcpp::List fast_dmvnorm_factor(const arma::mat& sigma) {
  const int p = sigma.n_cols;
  arma::mat R;
  if (p == 0 || !arma::chol(R, sigma)) {
    return Rcpp::List::create(Rcpp::_["ok"] = false,
                              Rcpp::_["rooti"] = R_NilValue,
                              Rcpp::_["log_const"] = R_NegInf);
  }
  arma::mat rooti = arma::inv(arma::trimatu(R));
  const double log_const =
    arma::sum(arma::log(rooti.diag())) - 0.5 * p * std::log(2.0 * M_PI);
  return Rcpp::List::create(Rcpp::_["ok"] = true,
                            Rcpp::_["rooti"] = rooti,
                            Rcpp::_["log_const"] = log_const);
}

// [[Rcpp::export]]
arma::vec fast_dmvnorm_rooti(const arma::mat& x, const arma::rowvec& mean, const arma::mat& rooti, double log_const) {
  int n = x.n_rows;
  if (n == 0) return arma::vec();
  arma::mat x_centered = x;
  x_centered.each_row() -= mean;
  arma::mat z = x_centered * rooti;
  return -0.5 * arma::sum(arma::square(z), 1) + log_const;
}


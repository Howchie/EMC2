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

// [[Rcpp::export]]
arma::vec fast_dmvnorm_rooti(const arma::mat& x, const arma::rowvec& mean, const arma::mat& rooti, double log_const) {
  int n = x.n_rows;
  if (n == 0) return arma::vec();
  arma::mat x_centered = x;
  x_centered.each_row() -= mean;
  arma::mat z = x_centered * rooti;
  return -0.5 * arma::sum(arma::square(z), 1) + log_const;
}


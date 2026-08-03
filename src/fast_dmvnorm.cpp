// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utility_functions.h"

// [[Rcpp::export]]
arma::vec fast_dmvnorm(const arma::mat& x, const arma::rowvec& mean, const arma::mat& sigma) {
  int n = x.n_rows;
  int p = x.n_cols;
  arma::vec out(n);
  
  arma::mat R;
  if (!arma::chol(R, sigma)) {
    out.fill(R_NegInf);
    return out;
  }
  arma::mat rooti = arma::inv(arma::trimatu(R));
  double log_const = arma::sum(arma::log(rooti.diag())) - 0.5 * p * std::log(2.0 * M_PI);
  
  for (int i = 0; i < n; i++) {
    arma::vec z = rooti.t() * arma::trans(x.row(i) - mean);
    out[i] = -0.5 * arma::dot(z, z) + log_const;
  }
  return out;
}

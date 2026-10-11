// Subject conditional frames for the shear move (R/shear_move.R): for each
// subject the upper Cholesky factor R_s of H_s + V^-1, its centre
// (H_s + V^-1)^-1 (lin_s + V^-1 mu) and log|R_s|, from the likelihood
// approximations stacked as a cube (a zero slice is a subject without one).
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List emc_shear_frames_cpp(const arma::cube& prec, const arma::mat& lin,
                                const arma::mat& Vi, const arma::vec& Vi_mu) {
  const arma::uword p = Vi.n_rows, n = lin.n_cols;
  arma::mat m(p, n);
  arma::cube R(p, p, n);
  arma::vec logdet(n);
  arma::mat Q, Rs;
  arma::vec rhs;
  for (arma::uword s = 0; s < n; ++s) {
    Q = prec.slice(s) + Vi;
    if (!arma::chol(Rs, Q)) return Rcpp::List::create(Rcpp::Named("ok") = false);
    R.slice(s) = Rs;
    rhs = lin.col(s) + Vi_mu;
    m.col(s) = arma::solve(arma::trimatu(Rs), arma::solve(arma::trimatl(Rs.t()), rhs));
    logdet(s) = arma::sum(arma::log(Rs.diag()));
  }
  return Rcpp::List::create(Rcpp::Named("ok") = true, Rcpp::Named("m") = m,
                            Rcpp::Named("R") = R, Rcpp::Named("logdet") = logdet);
}

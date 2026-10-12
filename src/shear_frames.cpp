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

// Triangular solves by plain substitution, as R's backsolve/forwardsolve: no
// condition check and no approximate fallback, so a near-singular frame gives
// the same (possibly non-finite) answer the R map does.
static arma::mat upper_solve(const arma::mat& R, const arma::mat& B) {
  const arma::uword p = R.n_rows;
  arma::mat X(B);
  for (arma::uword c = 0; c < X.n_cols; ++c)
    for (arma::uword ii = p; ii-- > 0;) {
      double v = X(ii, c);
      for (arma::uword k = ii + 1; k < p; ++k) v -= R(ii, k) * X(k, c);
      X(ii, c) = v / R(ii, ii);
    }
  return X;
}

// Solves R' X = B for upper triangular R (L = R' lower).
static arma::mat upper_t_solve(const arma::mat& R, const arma::mat& B) {
  const arma::uword p = R.n_rows;
  arma::mat X(B);
  for (arma::uword c = 0; c < X.n_cols; ++c)
    for (arma::uword ii = 0; ii < p; ++ii) {
      double v = X(ii, c);
      for (arma::uword k = 0; k < ii; ++k) v -= R(k, ii) * X(k, c);
      X(ii, c) = v / R(ii, ii);
    }
  return X;
}

static arma::mat lower_solve(const arma::mat& L, const arma::mat& B) {
  const arma::uword p = L.n_rows;
  arma::mat X(B);
  for (arma::uword c = 0; c < X.n_cols; ++c)
    for (arma::uword ii = 0; ii < p; ++ii) {
      double v = X(ii, c);
      for (arma::uword k = 0; k < ii; ++k) v -= L(ii, k) * X(k, c);
      X(ii, c) = v / L(ii, ii);
    }
  return X;
}

// The quantile map between standard normals truncated below at c0 and at c1
// (R/shear_move.R .emc_shear_qmap): z keeps its conditional upper-tail
// probability, and below the wall the mirror map applies. Returns false when
// the moved point is not finite.
static bool shear_qmap(double z, double c0, double c1, double& z1, double& log_jac) {
  if (std::isfinite(z) && z < c0) {
    if (!shear_qmap(-z, -c0, -c1, z1, log_jac)) return false;
    z1 = -z1;
    return true;
  }
  const double ls0 = R::pnorm5(c0, 0.0, 1.0, 0, 1);
  const double ls1 = R::pnorm5(c1, 0.0, 1.0, 0, 1);
  const double lf = R::pnorm5(z, 0.0, 1.0, 0, 1) - ls0;
  z1 = R::qnorm5(lf + ls1, 0.0, 1.0, 0, 1);
  log_jac = R::dnorm4(z, 0.0, 1.0, 1) - R::dnorm4(z1, 0.0, 1.0, 1) - ls0 + ls1;
  return true;
}

// Transport of every subject from frames (m0, R0) to frames (m1, R1), the
// shear and collapsed moves' map (R/shear_move.R .emc_shear_propose and
// .emc_shear_wall_map). A subject without a wall takes the affine frame map
// a1 = m1 + R1^-1 R0 (a0 - m0); one with a wall (unit normal U[, s], offset
// b[s], complement basis Uperp[, , s]) keeps its half-space: the complement
// moves by the affine map between the frames' marginals and u'a by the
// quantile map between the conditionals given it. Returns the moved
// subjects; per subject the wall map's log Jacobian less the plain one's
// (zero without a wall), whether the map was computed, and |R1 (a1 - m1)|^2.
// [[Rcpp::export]]
Rcpp::List emc_shear_transport_cpp(const arma::mat& alpha,
                                   const arma::mat& m0, const arma::cube& R0,
                                   const arma::vec& logdet0,
                                   const arma::mat& m1, const arma::cube& R1,
                                   const arma::vec& logdet1,
                                   const Rcpp::LogicalVector& has,
                                   const arma::mat& U, const arma::vec& b,
                                   const arma::cube& Uperp) {
  const arma::uword p = alpha.n_rows, n = alpha.n_cols;
  const bool walls = has.size() == (R_xlen_t) n;
  arma::mat out(p, n);
  arma::vec jac(n, arma::fill::zeros), r1sq(n);
  Rcpp::LogicalVector ok(n, true);
  arma::vec a, y0, y1, Qu0, Qu1, u, d;
  arma::mat W0, W1, S0, S1, L0, L1, Up;
  for (arma::uword s = 0; s < n; ++s) {
    a = alpha.col(s);
    const arma::mat& Ra = R0.slice(s);
    const arma::mat& Rb = R1.slice(s);
    if (!walls || has[s] != TRUE) {
      out.col(s) = m1.col(s) + upper_solve(Rb, Ra * (a - m0.col(s)));
    } else {
      u = U.col(s);
      Up = Uperp.slice(s);
      // Marginal covariances of the complement coordinates: Uperp' Q^-1 Uperp
      // with Q^-1 = R^-1 R^-T.
      W0 = upper_t_solve(Ra, Up);
      W1 = upper_t_solve(Rb, Up);
      S0 = W0.t() * W0;
      S1 = W1.t() * W1;
      if (!arma::chol(L0, S0, "lower") || !arma::chol(L1, S1, "lower")) {
        ok[s] = false; out.col(s).fill(NA_REAL); continue;
      }
      y0 = Up.t() * (a - m0.col(s));
      y1 = L1 * lower_solve(L0, y0);
      Qu0 = Ra.t() * (Ra * u);
      Qu1 = Rb.t() * (Rb * u);
      const double q0 = arma::dot(u, Qu0), q1 = arma::dot(u, Qu1);
      const double sc0 = 1.0 / std::sqrt(q0), sc1 = 1.0 / std::sqrt(q1);
      const double mc0 = arma::dot(u, m0.col(s)) - arma::dot(Qu0, Up * y0) / q0;
      const double mc1 = arma::dot(u, m1.col(s)) - arma::dot(Qu1, Up * y1) / q1;
      const double z0 = (arma::dot(u, a) - mc0) / sc0;
      const double w0 = (b(s) - mc0) / sc0;
      const double w1 = (b(s) - mc1) / sc1;
      double z1 = NA_REAL, lj = NA_REAL;
      if (!std::isfinite(z0) || !shear_qmap(z0, w0, w1, z1, lj)) {
        ok[s] = false; out.col(s).fill(NA_REAL); continue;
      }
      const double t1 = mc1 + sc1 * z1;
      if (!std::isfinite(t1)) { ok[s] = false; out.col(s).fill(NA_REAL); continue; }
      out.col(s) = m1.col(s) + Up * y1 + (t1 - arma::dot(u, m1.col(s))) * u;
      const double wall_jac = arma::sum(arma::log(L1.diag())) - arma::sum(arma::log(L0.diag())) +
        lj - std::log(sc0) + std::log(sc1);
      jac(s) = wall_jac - (logdet0(s) - logdet1(s));
    }
    d = Rb * (out.col(s) - m1.col(s));
    r1sq(s) = arma::dot(d, d);
    if (!out.col(s).is_finite()) ok[s] = false;
  }
  return Rcpp::List::create(Rcpp::Named("alpha") = out, Rcpp::Named("jac") = jac,
                            Rcpp::Named("ok") = ok, Rcpp::Named("r1sq") = r1sq);
}

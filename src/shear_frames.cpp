// Subject conditional frames for the shear move (R/shear_move.R): for each
// subject the upper Cholesky factor R_s of H_s + V^-1, its centre
// (H_s + V^-1)^-1 (lin_s + V^-1 mu) and log|R_s|, from the likelihood
// approximations stacked as a cube (a zero slice is a subject without one).
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
#include <RcppArmadillo.h>
#include <vector>

// [[Rcpp::export]]
Rcpp::List emc_shear_frames_cpp(const arma::cube& prec, const arma::mat& lin,
                                const arma::mat& Vi, const arma::vec& Vi_mu) {
  const arma::uword p = Vi.n_rows, n = lin.n_cols;
  arma::mat m(p, n);
  arma::cube R(p, p, n);
  arma::vec logdet(n);
  arma::mat Q, Rs;
  arma::vec rhs, z, centre;
  if (!Vi.is_finite() || !Vi_mu.is_finite())
    return Rcpp::List::create(Rcpp::Named("ok") = false);
  for (arma::uword s = 0; s < n; ++s) {
    Q = prec.slice(s) + Vi;
    if (!Q.is_finite() || !arma::chol(Rs, Q))
      return Rcpp::List::create(Rcpp::Named("ok") = false);
    R.slice(s) = Rs;
    rhs = lin.col(s) + Vi_mu;
    // A pseudo-inverse can discard directions in an extreme log-sd
    // probe yet return finite centres. The transport requires the actual
    // conditional frame: report failure instead of accepting that fallback.
    if (!rhs.is_finite() ||
        !arma::solve(z, arma::trimatl(Rs.t()), rhs, arma::solve_opts::no_approx) ||
        !arma::solve(centre, arma::trimatu(Rs), z, arma::solve_opts::no_approx) ||
        !centre.is_finite())
      return Rcpp::List::create(Rcpp::Named("ok") = false);
    m.col(s) = centre;
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
struct ShearTransport {
  arma::mat alpha;
  arma::vec jac, r1sq;
  std::vector<bool> ok;
};

static void shear_transport(const arma::mat& alpha,
                            const arma::mat& m0, const arma::cube& R0, const arma::vec& logdet0,
                            const arma::mat& m1, const arma::cube& R1, const arma::vec& logdet1,
                            const Rcpp::LogicalVector& has, const arma::mat& U,
                            const arma::vec& b, const arma::cube& Uperp, ShearTransport& res) {
  const arma::uword p = alpha.n_rows, n = alpha.n_cols;
  const bool walls = has.size() == (R_xlen_t) n;
  res.alpha.set_size(p, n);
  res.jac.zeros(n);
  res.r1sq.set_size(n);
  res.ok.assign(n, true);
  arma::mat& out = res.alpha;
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
        res.ok[s] = false; out.col(s).fill(NA_REAL); continue;
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
        res.ok[s] = false; out.col(s).fill(NA_REAL); continue;
      }
      const double t1 = mc1 + sc1 * z1;
      if (!std::isfinite(t1)) { res.ok[s] = false; out.col(s).fill(NA_REAL); continue; }
      out.col(s) = m1.col(s) + Up * y1 + (t1 - arma::dot(u, m1.col(s))) * u;
      const double wall_jac = arma::sum(arma::log(L1.diag())) - arma::sum(arma::log(L0.diag())) +
        lj - std::log(sc0) + std::log(sc1);
      res.jac(s) = wall_jac - (logdet0(s) - logdet1(s));
    }
    d = Rb * (out.col(s) - m1.col(s));
    res.r1sq(s) = arma::dot(d, d);
    if (!out.col(s).is_finite()) res.ok[s] = false;
  }
}

// [[Rcpp::export]]
Rcpp::List emc_shear_transport_cpp(const arma::mat& alpha,
                                   const arma::mat& m0, const arma::cube& R0,
                                   const arma::vec& logdet0,
                                   const arma::mat& m1, const arma::cube& R1,
                                   const arma::vec& logdet1,
                                   const Rcpp::LogicalVector& has,
                                   const arma::mat& U, const arma::vec& b,
                                   const arma::cube& Uperp) {
  ShearTransport res;
  shear_transport(alpha, m0, R0, logdet0, m1, R1, logdet1, has, U, b, Uperp, res);
  return Rcpp::List::create(Rcpp::Named("alpha") = res.alpha, Rcpp::Named("jac") = res.jac,
                            Rcpp::Named("ok") = Rcpp::wrap(res.ok), Rcpp::Named("r1sq") = res.r1sq);
}

// The covariance prior's log density at V, as .emc_shear_log_prior_v()
// computes it: the same factorisations, and sums accumulated as R's sum().
// NA when a block has no Cholesky factor.
#if defined(__GNUC__) && !defined(__clang__)
#define EMC2_STRICT_MATH __attribute__((optimize("no-fast-math")))
#else
#define EMC2_STRICT_MATH
#endif

// Strict IEEE arithmetic: a vectorised log differs from R's in the last bit.
EMC2_STRICT_MATH
static double shear_log_prior_v(const arma::mat& V, const Rcpp::List& spec) {
  double lp = 0.0;
  const Rcpp::IntegerVector g = spec["gamma_idx"];
  const bool uses_a_half = Rcpp::as<bool>(spec["uses_a_half"]);
  if (g.size() > 0) {
    long double acc = 0.0L;
    if (uses_a_half) {
      const double v = Rcpp::as<double>(spec["v"]);
      const Rcpp::NumericVector A = spec["A"];
      for (R_xlen_t i = 0; i < g.size(); ++i) {
        const int j = g[i] - 1;
        const double tau = 1 / V(j, j);
        const double term = (v / 2 + 1) * std::log(tau) -
          (v + 1) / 2 * std::log(v * tau + 1 / (A[j] * A[j]));
        acc += term;
      }
    } else {
      const Rcpp::NumericVector shape = spec["shape"], rate = spec["rate"];
      for (R_xlen_t i = 0; i < g.size(); ++i) {
        const int j = g[i] - 1;
        const double tau = 1 / V(j, j);
        const double term = (shape[j] + 1) * std::log(tau) - rate[j] * tau;
        acc += term;
      }
    }
    lp = lp + static_cast<double>(acc);
  }
  const Rcpp::List iw = spec["iw"];
  if (iw.size() == 0) return lp;
  const double v = Rcpp::as<double>(spec["v"]);
  const Rcpp::NumericVector A = spec["A"];
  for (R_xlen_t k = 0; k < iw.size(); ++k) {
    const Rcpp::IntegerVector b = iw[k];
    const int d = b.size();
    std::vector<double> R(static_cast<size_t>(d) * d, 0.0);
    for (int c = 0; c < d; ++c)
      for (int r = 0; r <= c; ++r) R[r + static_cast<size_t>(c) * d] = V(b[r] - 1, b[c] - 1);
    int info = 0;
    F77_CALL(dpotrf)("U", &d, R.data(), &d, &info FCONE);
    if (info != 0) return NA_REAL;
    long double sum_log = 0.0L;
    for (int j = 0; j < d; ++j) sum_log += std::log(R[j + static_cast<size_t>(j) * d]);
    std::vector<double> inv(R);
    F77_CALL(dpotri)("U", &d, inv.data(), &d, &info FCONE);
    if (info != 0) return NA_REAL;
    long double sum_inv = 0.0L;
    for (int j = 0; j < d; ++j) {
      const double a = A[b[j] - 1];
      sum_inv += std::log(v * inv[j + static_cast<size_t>(j) * d] + 1 / (a * a));
    }
    const double df = v + d - 1;
    lp = lp - (df + d + 1) * static_cast<double>(sum_log) -
      (v + d) / 2 * static_cast<double>(sum_inv);
  }
  return lp;
}

// The candidates of one shear round (R/shear_move.R .emc_shear_candidates):
// for each offset t[k] the group state moves by u = t[k] v (mean part, then
// log-sd part), the frames are rebuilt there from every subject's
// likelihood approximation and the subjects transported from the current
// frames (m0, R0).
// Per candidate: whether it is valid, the moved subjects, the new frames
// (m, R, logdet), the inverse and upper Cholesky factor of the moved V, the
// summed log N(a_s | mu, V) of the moved subjects and the transport's log
// Jacobian (frames plus walls, without the V-prior term), with the walls'
// part of it (wall maps less the plain maps they replace) separately.
// [[Rcpp::export]]
Rcpp::List emc_shear_candidates_cpp(const arma::vec& t, const arma::vec& v,
                                    const arma::mat& alpha, const arma::vec& mu,
                                    const arma::mat& V,
                                    const arma::mat& m0, const arma::cube& R0,
                                    const arma::vec& logdet0,
                                    const arma::cube& prec, const arma::mat& lin,
                                    const Rcpp::LogicalVector& has,
                                    const arma::mat& U, const arma::vec& b,
                                    const arma::cube& Uperp,
                                    Rcpp::Nullable<Rcpp::List> spec = R_NilValue) {
  const arma::uword p = alpha.n_rows, n = alpha.n_cols, K = t.n_elem;
  const double logdet0_sum = arma::accu(logdet0);
  const bool with_prior = spec.isNotNull();
  Rcpp::List spec_list;
  if (with_prior) spec_list = Rcpp::List(spec.get());
  Rcpp::List out(K);
  arma::mat m1(p, n), Rs, Q, cholV, Ri, Vi;
  arma::cube R1(p, p, n);
  arma::vec logdet1(n), rhs;
  arma::mat Z;
  ShearTransport tr;
  for (arma::uword k = 0; k < K; ++k) {
    const arma::vec u = t(k) * v;
    const arma::vec mu1 = mu + u.head(p);
    const arma::vec d = arma::exp(u.tail(p));
    const arma::mat V1 = V % (d * d.t());
    bool ok = arma::chol(cholV, V1);
    if (ok) {
      Ri = upper_solve(cholV, arma::eye(p, p));
      Vi = Ri * Ri.t();
      const arma::vec Vi_mu = Vi * mu1;
      for (arma::uword s = 0; s < n && ok; ++s) {
        Q = prec.slice(s) + Vi;
        if (!arma::chol(Rs, Q)) { ok = false; break; }
        R1.slice(s) = Rs;
        rhs = lin.col(s) + Vi_mu;
        m1.col(s) = upper_solve(Rs, upper_t_solve(Rs, rhs));
        logdet1(s) = arma::sum(arma::log(Rs.diag()));
      }
    }
    if (!ok) {
      out[k] = Rcpp::List::create(Rcpp::Named("ok") = false);
      continue;
    }
    shear_transport(alpha, m0, R0, logdet0, m1, R1, logdet1, has, U, b, Uperp, tr);
    bool all_ok = tr.alpha.is_finite();
    for (arma::uword s = 0; s < n && all_ok; ++s) all_ok = tr.ok[s];
    if (!all_ok) {
      out[k] = Rcpp::List::create(Rcpp::Named("ok") = false);
      continue;
    }
    // sum_s log N(a_s | mu1, V1) without the 2 pi constant.
    Z = upper_t_solve(cholV, tr.alpha.each_col() - mu1);
    const double hier = -0.5 * arma::accu(Z % Z) - (double) n * arma::sum(arma::log(cholV.diag()));
    const double jac = logdet0_sum - arma::accu(logdet1) + arma::accu(tr.jac);
    Rcpp::List cand = Rcpp::List::create(
      Rcpp::Named("ok") = true, Rcpp::Named("alpha") = tr.alpha,
      Rcpp::Named("m") = m1, Rcpp::Named("R") = R1,
      Rcpp::Named("logdet") = logdet1, Rcpp::Named("Vi") = Vi,
      Rcpp::Named("chol_V") = cholV, Rcpp::Named("hier") = hier,
      Rcpp::Named("jac") = jac, Rcpp::Named("wall_jac") = arma::accu(tr.jac));
    if (with_prior) {
      cand["V"] = V1;
      cand["lp_v"] = shear_log_prior_v(V1, spec_list);
    }
    out[k] = cand;
  }
  return out;
}

// A group state's factors and frames (R/shear_move.R .emc_shear_state): the
// upper Cholesky factor of V and its inverse as R's chol() and chol2inv()
// return them, the subjects' frames at V^-1 and V^-1 mu, and the covariance
// prior's density. ok = FALSE when mu or V is not finite, V has no factor or
// a frame fails.
// [[Rcpp::export(rng = false)]]
Rcpp::List emc_shear_state_cpp(Rcpp::NumericVector mu, Rcpp::NumericMatrix V,
                               const arma::cube& prec, const arma::mat& lin,
                               Rcpp::List spec) {
  const int p = V.nrow();
  const Rcpp::List failed = Rcpp::List::create(Rcpp::Named("ok") = false);
  if (V.ncol() != p || mu.size() != p) return failed;
  for (int i = 0; i < p; ++i) if (!R_FINITE(mu[i])) return failed;
  for (R_xlen_t i = 0; i < V.size(); ++i) if (!R_FINITE(V[i])) return failed;
  Rcpp::NumericMatrix cholV(p, p);
  for (int c = 0; c < p; ++c)
    for (int r = 0; r <= c; ++r) cholV(r, c) = V(r, c);
  int info = 0;
  F77_CALL(dpotrf)("U", &p, cholV.begin(), &p, &info FCONE);
  if (info != 0) return failed;
  Rcpp::NumericMatrix Vi(p, p);
  for (int c = 0; c < p; ++c)
    for (int r = 0; r <= c; ++r) Vi(r, c) = cholV(r, c);
  F77_CALL(dpotri)("U", &p, Vi.begin(), &p, &info FCONE);
  if (info != 0) return failed;
  for (int c = 0; c < p; ++c)
    for (int r = c + 1; r < p; ++r) Vi(r, c) = Vi(c, r);
  for (R_xlen_t i = 0; i < Vi.size(); ++i) if (!R_FINITE(Vi[i])) return failed;
  Rcpp::NumericVector Vi_mu(p);
  const double one = 1.0, zero = 0.0;
  const int inc = 1;
  F77_CALL(dgemv)("N", &p, &p, &one, Vi.begin(), &p, mu.begin(), &inc, &zero,
                  Vi_mu.begin(), &inc FCONE);
  const Rcpp::List frames = emc_shear_frames_cpp(
    prec, lin, Rcpp::as<arma::mat>(Vi), Rcpp::as<arma::vec>(Vi_mu));
  if (!Rcpp::as<bool>(frames["ok"])) return failed;
  const SEXP dn = Rf_getAttrib(V, R_DimNamesSymbol);
  if (!Rf_isNull(dn)) Rf_setAttrib(cholV, R_DimNamesSymbol, dn);
  return Rcpp::List::create(
    Rcpp::Named("ok") = true, Rcpp::Named("Vi") = Vi, Rcpp::Named("chol_V") = cholV,
    Rcpp::Named("m") = frames["m"], Rcpp::Named("R") = frames["R"],
    Rcpp::Named("logdet") = frames["logdet"],
    Rcpp::Named("lp_v") = shear_log_prior_v(Rcpp::as<arma::mat>(V), spec));
}

// log Lhat_s of each subject's column of alpha: lin_s' a - a' H_s a / 2 for
// the subjects listed in `has` (1-based, the slices of H in order), 0 for the
// rest; the products and sums are R's.
// [[Rcpp::export(rng = false)]]
EMC2_STRICT_MATH
Rcpp::NumericVector emc_collapsed_log_lhat_cpp(Rcpp::NumericMatrix alpha,
                                               Rcpp::NumericVector H,
                                               Rcpp::NumericMatrix lin,
                                               Rcpp::IntegerVector has) {
  const int p = alpha.nrow(), n = alpha.ncol();
  Rcpp::NumericVector out(n);
  std::vector<double> Ha(p);
  const double one = 1.0, zero = 0.0;
  const int inc = 1;
  for (R_xlen_t k = 0; k < has.size(); ++k) {
    const int s = has[k] - 1;
    if (s < 0 || s >= n) Rcpp::stop("log lhat: subject out of range");
    const double* a = &alpha[static_cast<R_xlen_t>(s) * p];
    const double* Hs = &H[static_cast<R_xlen_t>(k) * p * p];
    F77_CALL(dgemv)("N", &p, &p, &one, Hs, &p, a, &inc, &zero, Ha.data(), &inc FCONE);
    long double quad = 0.0L, linear = 0.0L;
    for (int i = 0; i < p; ++i) {
      volatile double term = a[i] * Ha[i];
      quad += term;
      volatile double lterm = lin(i, k) * a[i];
      linear += lterm;
    }
    volatile double half = -0.5 * static_cast<double>(quad);
    out[s] = half + static_cast<double>(linear);
  }
  return out;
}

// Each candidate's quadratic likelihood surrogate summed over subjects:
// sum_s lin_s' a_s - a_s' P_s a_s / 2, with P_s the slices of `cube` (zero for
// a subject without a fit). Candidates are parameters x subjects.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector emc_shear_surrogate_cpp(Rcpp::List cands, Rcpp::NumericVector cube,
                                            Rcpp::NumericMatrix lin) {
  const int p = lin.nrow(), n = lin.ncol();
  if (cube.size() != static_cast<R_xlen_t>(p) * p * n)
    Rcpp::stop("emc_shear_surrogate_cpp: cube and lin disagree");
  const double* P = REAL(cube);
  const R_xlen_t K = cands.size();
  Rcpp::NumericVector out(K);
  for (R_xlen_t k = 0; k < K; ++k) {
    Rcpp::NumericMatrix a(cands[k]);
    if (a.nrow() != p || a.ncol() != n)
      Rcpp::stop("emc_shear_surrogate_cpp: candidate dimensions");
    double total = 0.0;
    for (int s = 0; s < n; ++s) {
      const double* as = &a(0, s);
      const double* Ps = P + static_cast<R_xlen_t>(s) * p * p;
      double lin_term = 0.0, quad = 0.0;
      for (int j = 0; j < p; ++j) {
        lin_term += lin(j, s) * as[j];
        double row = 0.0;
        for (int i = 0; i < p; ++i) row += Ps[i + static_cast<R_xlen_t>(j) * p] * as[i];
        quad += row * as[j];
      }
      total += lin_term - 0.5 * quad;
    }
    out[k] = total;
  }
  return out;
}

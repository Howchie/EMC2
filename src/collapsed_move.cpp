// Inner chain of the collapsed group move (R/collapsed_move.R): the
// Laplace-collapsed group posterior phat(phi) with its gradient, and
// Hamiltonian trajectories on it. phi = (mu, the log-Cholesky factor of every
// inverse-Wishart block of V, the log sd of every gamma-prior parameter).
//
// The model list carries the chart (0-based positions), the subjects'
// quadratic approximations (prec as a p x p x n cube, lin as p x n) and the
// priors, as packed by .emc_collapsed_cpp_model().

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

namespace {

struct Block {
  arma::uvec idx;  // parameters of the block
  arma::uvec at;   // positions in phi, column-major lower triangle
};

struct CollapsedModel {
  int p = 0;
  int dim = 0;
  std::vector<Block> blocks;
  arma::uvec gamma, gamma_at;
  arma::cube H;
  arma::mat lin;
  arma::vec m0;
  arma::mat P0;
  bool uses_a_half = true;
  double v = 2.0;
  arma::vec A, shape, rate;
};

CollapsedModel unpack(const Rcpp::List& m) {
  CollapsedModel out;
  out.p = Rcpp::as<int>(m["p"]);
  out.dim = Rcpp::as<int>(m["dim"]);
  Rcpp::List bi = m["block_idx"], ba = m["block_at"];
  for (R_xlen_t k = 0; k < bi.size(); ++k) {
    Block b;
    b.idx = Rcpp::as<arma::uvec>(bi[k]);
    b.at = Rcpp::as<arma::uvec>(ba[k]);
    out.blocks.push_back(b);
  }
  out.gamma = Rcpp::as<arma::uvec>(m["gamma"]);
  out.gamma_at = Rcpp::as<arma::uvec>(m["gamma_at"]);
  Rcpp::NumericVector Hv = m["H"];
  Rcpp::IntegerVector hd = Hv.attr("dim");
  out.H = arma::cube(Hv.begin(), hd[0], hd[1], hd[2], false, true);
  out.lin = Rcpp::as<arma::mat>(m["lin"]);
  out.m0 = Rcpp::as<arma::vec>(m["m0"]);
  out.P0 = Rcpp::as<arma::mat>(m["P0"]);
  out.uses_a_half = Rcpp::as<bool>(m["uses_a_half"]);
  out.v = Rcpp::as<double>(m["v"]);
  out.A = Rcpp::as<arma::vec>(m["A"]);
  out.shape = Rcpp::as<arma::vec>(m["shape"]);
  out.rate = Rcpp::as<arma::vec>(m["rate"]);
  return out;
}

// log phat(phi) and its gradient; false when phi gives no usable state.
bool surrogate(const arma::vec& phi, const CollapsedModel& M, double& lp,
               arma::vec& grad) {
  const int p = M.p;
  arma::vec mu = phi.head(p);
  arma::mat V(p, p, arma::fill::zeros);
  std::vector<arma::mat> Ls(M.blocks.size());
  double log_jac = 0.0;
  for (size_t k = 0; k < M.blocks.size(); ++k) {
    const Block& b = M.blocks[k];
    const arma::uword d = b.idx.n_elem;
    arma::mat L(d, d, arma::fill::zeros);
    arma::uword q = 0;
    for (arma::uword j = 0; j < d; ++j)
      for (arma::uword i = j; i < d; ++i) L(i, j) = phi(b.at(q++));
    for (arma::uword j = 0; j < d; ++j) {
      log_jac += (double)(d - j + 1) * L(j, j);
      L(j, j) = std::exp(L(j, j));
    }
    V.submat(b.idx, b.idx) = L * L.t();
    Ls[k] = L;
  }
  for (arma::uword g = 0; g < M.gamma.n_elem; ++g) {
    const double lg = phi(M.gamma_at(g));
    V(M.gamma(g), M.gamma(g)) = std::exp(2.0 * lg);
    log_jac += 2.0 * lg;
  }
  if (!V.is_finite()) return false;
  arma::mat cV;
  if (!arma::chol(cV, V)) return false;
  arma::mat cVi;
  if (!arma::inv(cVi, arma::trimatu(cV))) return false;
  arma::mat Vi = cVi * cVi.t();

  // Covariance prior (.emc_shear_log_prior_v) and its gradient in V.
  double lp_v = 0.0;
  arma::mat G_V(p, p, arma::fill::zeros);
  for (arma::uword g = 0; g < M.gamma.n_elem; ++g) {
    const arma::uword k = M.gamma(g);
    const double tau = 1.0 / V(k, k);
    double dtau;
    if (M.uses_a_half) {
      const double c = 1.0 / (M.A(k) * M.A(k));
      lp_v += (M.v / 2.0 + 1.0) * std::log(tau) - (M.v + 1.0) / 2.0 * std::log(M.v * tau + c);
      dtau = (M.v / 2.0 + 1.0) / tau - (M.v + 1.0) / 2.0 * M.v / (M.v * tau + c);
    } else {
      lp_v += (M.shape(k) + 1.0) * std::log(tau) - M.rate(k) * tau;
      dtau = (M.shape(k) + 1.0) / tau - M.rate(k);
    }
    G_V(k, k) += -dtau * tau * tau;
  }
  for (const Block& b : M.blocks) {
    const double d = (double)b.idx.n_elem;
    arma::mat Vb = V.submat(b.idx, b.idx);
    arma::mat Rb;
    if (!arma::chol(Rb, Vb)) return false;
    arma::mat Rbi;
    if (!arma::inv(Rbi, arma::trimatu(Rb))) return false;
    arma::mat Lb = Rbi * Rbi.t();
    const double df = M.v + d - 1.0;
    lp_v -= (df + d + 1.0) * arma::accu(arma::log(Rb.diag()));
    arma::vec wk(b.idx.n_elem);
    for (arma::uword k = 0; k < b.idx.n_elem; ++k) {
      const double c = 1.0 / (M.A(b.idx(k)) * M.A(b.idx(k)));
      lp_v -= (M.v + d) / 2.0 * std::log(M.v * Lb(k, k) + c);
      wk(k) = M.v / (M.v * Lb(k, k) + c);
    }
    G_V.submat(b.idx, b.idx) += -(df + d + 1.0) / 2.0 * Lb +
      (M.v + d) / 2.0 * (Lb * arma::diagmat(wk) * Lb);
  }
  if (!std::isfinite(lp_v)) return false;
  arma::vec dm = mu - M.m0;
  lp = -0.5 * arma::dot(dm, M.P0 * dm) + lp_v + log_jac;

  // Subjects: c_s = log Int Lhat_s(a) N(a | mu, V) da up to a constant.
  const double half_ld_vi = -arma::accu(arma::log(cV.diag()));
  arma::vec Vi_mu = Vi * mu;
  const double mu_Vi_mu = arma::dot(mu, Vi_mu);
  arma::mat S(p, p, arma::fill::zeros);
  arma::vec dsum(p, arma::fill::zeros);
  const arma::uword n_fit = M.H.n_slices;
  arma::mat R, Ri;
  for (arma::uword s = 0; s < n_fit; ++s) {
    arma::mat Q = M.H.slice(s) + Vi;
    if (!arma::chol(R, Q)) return false;
    arma::vec b = M.lin.col(s) + Vi_mu;
    arma::vec z, m;
    if (!arma::solve(z, arma::trimatl(R.t()), b, arma::solve_opts::no_approx)) return false;
    lp += half_ld_vi - arma::accu(arma::log(R.diag())) + 0.5 * arma::dot(z, z) - 0.5 * mu_Vi_mu;
    if (!arma::solve(m, arma::trimatu(R), z, arma::solve_opts::no_approx)) return false;
    arma::vec d = m - mu;
    if (!arma::inv(Ri, arma::trimatu(R))) return false;
    S += Ri * Ri.t() + d * d.t();
    dsum += d;
  }
  if (!std::isfinite(lp)) return false;
  // d/dLambda of sum_s c_s is (n V - S) / 2, and dLambda = -Lambda dV Lambda.
  arma::mat G_lam = 0.5 * ((double)n_fit * V - S);
  G_V += -Vi * G_lam * Vi;
  G_V = 0.5 * (G_V + G_V.t());

  grad.zeros(M.dim);
  grad.head(p) = Vi * dsum - M.P0 * dm;
  for (size_t k = 0; k < M.blocks.size(); ++k) {
    const Block& b = M.blocks[k];
    const arma::uword d = b.idx.n_elem;
    const arma::mat& L = Ls[k];
    arma::mat dL = 2.0 * G_V.submat(b.idx, b.idx) * L;
    arma::uword q = 0;
    for (arma::uword j = 0; j < d; ++j)
      for (arma::uword i = j; i < d; ++i) {
        double g = dL(i, j);
        if (i == j) g = g * L(j, j) + (double)(d - j + 1);
        grad(b.at(q++)) = g;
      }
  }
  for (arma::uword g = 0; g < M.gamma.n_elem; ++g) {
    const arma::uword k = M.gamma(g);
    grad(M.gamma_at(g)) = 2.0 * G_V(k, k) * V(k, k) + 2.0;
  }
  return grad.is_finite();
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List collapsed_surrogate_cpp(const arma::vec& phi, const Rcpp::List& model) {
  CollapsedModel M = unpack(model);
  double lp = R_NegInf;
  arma::vec grad;
  if (!surrogate(phi, M, lp, grad))
    return Rcpp::List::create(Rcpp::Named("lp") = R_NegInf, Rcpp::Named("grad") = R_NilValue);
  return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// n_traj Hamiltonian trajectories on phat in whitened coordinates phi = C x,
// each accepted or rejected against phat. Each trajectory jitters eps by a
// uniform factor in [0.8, 1.2] and takes ceil(len / eps) steps. Random draws,
// in order per trajectory: the jitter, the momentum, the acceptance.
// [[Rcpp::export]]
Rcpp::List collapsed_trajectories_cpp(const arma::vec& phi0, const arma::mat& C,
                                      double eps, int n_traj, double len,
                                      const Rcpp::List& model) {
  CollapsedModel M = unpack(model);
  Rcpp::NumericVector acc(n_traj);
  arma::vec phi = phi0, grad;
  double lp;
  if (!surrogate(phi, M, lp, grad))
    return Rcpp::List::create(Rcpp::Named("phi") = phi0, Rcpp::Named("lp") = R_NegInf,
                              Rcpp::Named("acc") = acc);
  const arma::uword dim = phi.n_elem;
  arma::vec x, v, v0, g, g_new;
  for (int t = 0; t < n_traj; ++t) {
    const double e = eps * (0.8 + 0.4 * R::unif_rand());
    const int L = std::max(1, (int)std::ceil(len / e));
    v0.set_size(dim);
    for (arma::uword i = 0; i < dim; ++i) v0(i) = R::norm_rand();
    x = phi;
    g = C.t() * grad;
    v = v0 + 0.5 * e * g;
    double lp_new = lp;
    arma::vec grad_new;
    bool ok = true;
    for (int l = 1; l <= L; ++l) {
      x += e * (C * v);
      if (!surrogate(x, M, lp_new, grad_new)) { ok = false; break; }
      g_new = C.t() * grad_new;
      v += (l < L ? e : 0.5 * e) * g_new;
    }
    if (!ok) { acc[t] = 0.0; continue; }
    const double la = lp_new - lp - 0.5 * (arma::dot(v, v) - arma::dot(v0, v0));
    const double a = std::isfinite(la) ? std::min(1.0, std::exp(la)) : 0.0;
    acc[t] = a;
    if (R::unif_rand() < a) {
      phi = x;
      lp = lp_new;
      grad = grad_new;
    }
  }
  return Rcpp::List::create(Rcpp::Named("phi") = phi, Rcpp::Named("lp") = lp,
                            Rcpp::Named("acc") = acc);
}

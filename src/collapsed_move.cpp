// Inner chain of the collapsed group move (R/collapsed_move.R): the
// Laplace-collapsed group posterior phat(phi) with its gradient, and
// Hamiltonian trajectories on it. phi = (mu, the log-Cholesky factor of every
// inverse-Wishart block of V, the log sd of every gamma-prior parameter).
//
// The model list carries the chart (0-based positions), the subjects'
// quadratic approximations (prec as a p x p x n cube, lin as p x n) and the
// priors, as packed by .emc_collapsed_cpp_model().

#include <RcppArmadillo.h>
#include <algorithm>
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
  // The approximations again, subject-innermost so the per-subject
  // factorisations below run in SIMD lanes: element (i, j), i >= j, of
  // subject s at Hl[(i * p + j) * n + s], and lin[i, s] at ll[i * n + s].
  std::vector<double> Hl, ll;
};

// Scratch for the subject terms, subject-innermost like Hl.
struct Work {
  std::vector<double> L, Li, z, m, inv_d, ld, sum, n1, ninf;
  explicit Work(const CollapsedModel& M) {
    const size_t p = M.p, n = M.H.n_slices;
    L.resize(p * p * n);
    Li.resize(p * p * n);
    z.resize(p * n);
    m.resize(p * n);
    inv_d.resize(p * n);
    ld.resize(n);
    sum.resize(n);
    n1.resize(2 * n);
    ninf.resize(2 * n);
  }
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
  const arma::uword p = out.H.n_rows, n = out.H.n_slices;
  out.Hl.assign(p * p * n, 0.0);
  out.ll.assign(p * n, 0.0);
  for (arma::uword s = 0; s < n; ++s)
    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword j = 0; j <= i; ++j) out.Hl[(i * p + j) * n + s] = out.H(i, j, s);
      out.ll[i * n + s] = out.lin(i, s);
    }
  return out;
}

// Sum over subjects of c_s + log det(Vi)^(-1/2), where c_s is
// log Int Lhat_s(a) N(a | mu, V) da up to a constant. With Q_s = H_s + Vi =
// L_s L_s' and L_s w_s = lin_s - H_s mu,
//   c_s = log det(Vi)^(1/2) - log det(L_s) + w_s'w_s / 2 + log Lhat_s(mu),
// which never forms Vi mu: the textbook form b'Q^(-1)b - mu'Vi mu with
// b = lin + Vi mu cancels catastrophically when Vi is large. The posterior
// means give d_s = m_s - mu = L_s'^(-1) w_s, and it accumulates
// S = sum_s Q_s^(-1) + d_s d_s' and dsum = sum_s d_s.
// All subjects are factorised together, one SIMD lane per subject. False when
// some Q_s is not numerically positive definite.
bool subject_terms(const CollapsedModel& M, Work& W, const arma::mat& Vi,
                   const arma::vec& mu, double& lp, arma::mat& S, arma::vec& dsum) {
  const arma::uword p = M.p, n = M.H.n_slices;
  const double* H = M.Hl.data();
  double* L = W.L.data();
  double* Li = W.Li.data();
  double* z = W.z.data();
  double* m = W.m.data();
  double* id = W.inv_d.data();
  double* ld = W.ld.data();
  bool bad = false;
  for (arma::uword s = 0; s < n; ++s) ld[s] = 0.0;
  // z = lin - H mu, and the sum of log Lhat_s(mu) = mu'lin_s - mu'H_s mu / 2.
  double lhat = 0.0;
  for (arma::uword i = 0; i < p; ++i) {
    double* zi = z + i * n;
    const double* li = M.ll.data() + i * n;
#pragma omp simd
    for (arma::uword s = 0; s < n; ++s) zi[s] = 0.0;
    for (arma::uword j = 0; j < p; ++j) {
      const double* Hij = H + (i >= j ? i * p + j : j * p + i) * n;
      const double muj = mu(j);
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) zi[s] += Hij[s] * muj;
    }
    const double mui = mu(i);
    double a = 0.0;
#pragma omp simd reduction(+ : a)
    for (arma::uword s = 0; s < n; ++s) {
      a += mui * (li[s] - 0.5 * zi[s]);
      zi[s] = li[s] - zi[s];
    }
    lhat += a;
  }
  // Cholesky factors and forward solves, column by column.
  for (arma::uword j = 0; j < p; ++j) {
    double* Ljj = L + (j * p + j) * n;
    const double* Hjj = H + (j * p + j) * n;
    const double vjj = Vi(j, j);
#pragma omp simd
    for (arma::uword s = 0; s < n; ++s) Ljj[s] = Hjj[s] + vjj;
    for (arma::uword k = 0; k < j; ++k) {
      const double* Ljk = L + (j * p + k) * n;
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) Ljj[s] -= Ljk[s] * Ljk[s];
    }
    double* idj = id + j * n;
#pragma omp simd reduction(|| : bad)
    for (arma::uword s = 0; s < n; ++s) {
      const double piv = Ljj[s];
      bad = bad || !(piv > 0.0);
      ld[s] += std::log(piv);
      Ljj[s] = std::sqrt(piv);
      idj[s] = 1.0 / Ljj[s];
    }
    if (bad) return false;
    for (arma::uword i = j + 1; i < p; ++i) {
      double* Lij = L + (i * p + j) * n;
      const double* Hij = H + (i * p + j) * n;
      const double vij = Vi(i, j);
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) Lij[s] = Hij[s] + vij;
      for (arma::uword k = 0; k < j; ++k) {
        const double* Lik = L + (i * p + k) * n;
        const double* Ljk = L + (j * p + k) * n;
#pragma omp simd
        for (arma::uword s = 0; s < n; ++s) Lij[s] -= Lik[s] * Ljk[s];
      }
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) Lij[s] *= idj[s];
    }
    double* zj = z + j * n;
    for (arma::uword k = 0; k < j; ++k) {
      const double* Ljk = L + (j * p + k) * n;
      const double* zk = z + k * n;
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) zj[s] -= Ljk[s] * zk[s];
    }
#pragma omp simd
    for (arma::uword s = 0; s < n; ++s) zj[s] *= idj[s];
  }
  double acc = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    const double* zj = z + j * n;
#pragma omp simd reduction(+ : acc)
    for (arma::uword s = 0; s < n; ++s) acc += zj[s] * zj[s];
  }
  double sld = 0.0;
#pragma omp simd reduction(+ : sld)
  for (arma::uword s = 0; s < n; ++s) sld += ld[s];
  lp = 0.5 * (acc - sld) + lhat;
  // Back solves L_s' d_s = w_s, d_s kept in m.
  for (arma::uword ii = p; ii-- > 0;) {
    double* mi = m + ii * n;
    const double* zi = z + ii * n;
#pragma omp simd
    for (arma::uword s = 0; s < n; ++s) mi[s] = zi[s];
    for (arma::uword k = ii + 1; k < p; ++k) {
      const double* Lki = L + (k * p + ii) * n;
      const double* mk = m + k * n;
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) mi[s] -= Lki[s] * mk[s];
    }
    const double* idi = id + ii * n;
#pragma omp simd
    for (arma::uword s = 0; s < n; ++s) mi[s] *= idi[s];
  }
  for (arma::uword i = 0; i < p; ++i) {
    const double* mi = m + i * n;
    double ds = 0.0;
#pragma omp simd reduction(+ : ds)
    for (arma::uword s = 0; s < n; ++s) ds += mi[s];
    dsum(i) = ds;
  }
  // Li = L^(-1), lower triangular, column by column.
  for (arma::uword j = 0; j < p; ++j) {
    const double* idj = id + j * n;
    double* Lijj = Li + (j * p + j) * n;
#pragma omp simd
    for (arma::uword s = 0; s < n; ++s) Lijj[s] = idj[s];
    for (arma::uword i = j + 1; i < p; ++i) {
      double* Liij = Li + (i * p + j) * n;
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) Liij[s] = 0.0;
      for (arma::uword k = j; k < i; ++k) {
        const double* Lik = L + (i * p + k) * n;
        const double* Likj = Li + (k * p + j) * n;
#pragma omp simd
        for (arma::uword s = 0; s < n; ++s) Liij[s] -= Lik[s] * Likj[s];
      }
      const double* idi = id + i * n;
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) Liij[s] *= idi[s];
    }
  }
  // Reject numerically singular factors, as a triangular solve does: the
  // reciprocal condition number of L_s in the 1- or infinity-norm below
  // machine epsilon.
  double* sum = W.sum.data();
  for (int f = 0; f < 2; ++f) {
    const double* T = f == 0 ? L : Li;
    double* n1 = W.n1.data() + f * n;
    double* ninf = W.ninf.data() + f * n;
    for (arma::uword s = 0; s < n; ++s) n1[s] = ninf[s] = 0.0;
    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword s = 0; s < n; ++s) sum[s] = 0.0;
      for (arma::uword i = j; i < p; ++i) {
        const double* Tij = T + (i * p + j) * n;
#pragma omp simd
        for (arma::uword s = 0; s < n; ++s) sum[s] += std::abs(Tij[s]);
      }
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) n1[s] = std::max(n1[s], sum[s]);
    }
    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword s = 0; s < n; ++s) sum[s] = 0.0;
      for (arma::uword j = 0; j <= i; ++j) {
        const double* Tij = T + (i * p + j) * n;
#pragma omp simd
        for (arma::uword s = 0; s < n; ++s) sum[s] += std::abs(Tij[s]);
      }
#pragma omp simd
      for (arma::uword s = 0; s < n; ++s) ninf[s] = std::max(ninf[s], sum[s]);
    }
  }
  const double eps = std::numeric_limits<double>::epsilon();
  const double* n1 = W.n1.data();
  const double* ninf = W.ninf.data();
#pragma omp simd reduction(|| : bad)
  for (arma::uword s = 0; s < n; ++s)
    bad = bad || !(n1[s] * n1[n + s] * eps < 1.0) || !(ninf[s] * ninf[n + s] * eps < 1.0);
  if (bad) return false;
  // S(i, j) = sum_s [sum_{k >= i} Li(k, i) Li(k, j) + d_i d_j], i >= j.
  for (arma::uword i = 0; i < p; ++i)
    for (arma::uword j = 0; j <= i; ++j) {
      const double* mi = m + i * n;
      const double* mj = m + j * n;
      double a = 0.0;
#pragma omp simd reduction(+ : a)
      for (arma::uword s = 0; s < n; ++s) a += mi[s] * mj[s];
      for (arma::uword k = i; k < p; ++k) {
        const double* Lki = Li + (k * p + i) * n;
        const double* Lkj = Li + (k * p + j) * n;
#pragma omp simd reduction(+ : a)
        for (arma::uword s = 0; s < n; ++s) a += Lki[s] * Lkj[s];
      }
      S(i, j) = a;
      S(j, i) = a;
    }
  return std::isfinite(lp) && S.is_finite();
}

// log phat(phi) and its gradient; false when phi gives no usable state.
bool surrogate(const arma::vec& phi, const CollapsedModel& M, Work& W, double& lp,
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
  arma::mat S(p, p);
  arma::vec dsum(p);
  const arma::uword n_fit = M.H.n_slices;
  double lp_s;
  if (!subject_terms(M, W, Vi, mu, lp_s, S, dsum)) return false;
  lp += lp_s + (double)n_fit * half_ld_vi;
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
  Work W(M);
  double lp = R_NegInf;
  arma::vec grad;
  if (!surrogate(phi, M, W, lp, grad))
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
  Work W(M);
  Rcpp::NumericVector acc(n_traj);
  arma::vec phi = phi0, grad;
  double lp;
  if (!surrogate(phi, M, W, lp, grad))
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
      if (!surrogate(x, M, W, lp_new, grad_new)) { ok = false; break; }
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

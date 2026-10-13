// Ensemble group update (R/ensemble_move.R) in C++: the pool draws from the
// frozen proposals, the stacked pools with their proposal densities, and the
// Gibbs rounds on the pool-augmented target for the standard and
// diagonal-gamma hierarchies. Random numbers come from R's generator in the
// order the R rounds (.emc_ensemble_rounds_r) consume them, so a seed gives
// the same chain up to rounding.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

const double LOG_2PI = std::log(2.0 * M_PI);

// A plain R numeric vector (no dim attribute) from an Armadillo vector.
inline Rcpp::NumericVector as_vector(const arma::vec& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

// Pool proposal density from a squared Mahalanobis distance: the mixture of
// N(mu, V) and a Student-t with scale matrix c V, both sharing V's factor.
// The terms that depend only on the dimension and the mixture are fixed per
// call.
struct LogQ {
  double a0, b0, c_nu, half_nu_p;
  LogQ(int p, double c, double tail_weight, double tail_df) {
    const double nu = tail_df;
    a0 = std::log1p(-tail_weight) - 0.5 * p * LOG_2PI;
    b0 = std::log(tail_weight) + R::lgammafn((nu + p) / 2.0) - R::lgammafn(nu / 2.0) -
      0.5 * (p * std::log(nu * M_PI) + p * std::log(c));
    c_nu = c * nu;
    half_nu_p = (nu + p) / 2.0;
  }
  double operator()(double m2, double logdet) const {
    const double a = a0 - 0.5 * (logdet + m2);
    const double b = b0 - 0.5 * logdet - half_nu_p * std::log1p(m2 / c_nu);
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
  }
};

// Covariance blocks of the standard hierarchy: the parameters sharing a group
// label, in the order the labels first appear, each block's parameters in
// ascending order. Unblocked parameters take a gamma update on their own.
struct BlockPlan {
  std::vector<std::vector<arma::uword>> blocks;
  std::vector<arma::uword> unblocked;
  arma::vec block_dim;
};

BlockPlan block_plan(const Rcpp::LogicalVector& is_blocked,
                     const Rcpp::IntegerVector& par_group) {
  const int p = is_blocked.size();
  BlockPlan plan;
  plan.block_dim.set_size(p);
  std::vector<int> seen;
  for (int i = 0; i < p; ++i) {
    if (!is_blocked[i]) {
      plan.unblocked.push_back(i);
      plan.block_dim(i) = 1.0;
      continue;
    }
    const int g = par_group[i];
    if (std::find(seen.begin(), seen.end(), g) != seen.end()) continue;
    seen.push_back(g);
    std::vector<arma::uword> idx;
    for (int k = i; k < p; ++k)
      if (is_blocked[k] && par_group[k] == g) idx.push_back(k);
    for (arma::uword k : idx) plan.block_dim(k) = static_cast<double>(idx.size());
    plan.blocks.push_back(idx);
  }
  return plan;
}

// One inverse-Wishart draw with scale S and df degrees of freedom: a Wishart
// draw W with scale S^-1 through its Bartlett factor (diagonal chi-squares,
// then the upper triangle's normals column by column), inverted. W = U'U with
// U = Z C upper triangular (both factors are), so U is the draw's precision
// factor and W its precision.
struct IWishDraw {
  arma::mat Sigma, W, U;
};

IWishDraw riwish(double df, arma::mat S, const char* context) {
  const arma::uword d = S.n_rows;
  S = 0.5 * (S + S.t());
  if (!S.is_finite())
    Rcpp::stop("%s: inverse-Wishart scale contains non-finite values", context);
  arma::mat S_inv;
  if (!arma::inv_sympd(S_inv, S))
    Rcpp::stop("%s: inverse-Wishart scale is not numerically positive definite", context);
  arma::mat CC;
  if (!arma::chol(CC, S_inv))
    Rcpp::stop("%s: inverse-Wishart scale could not be factorised (chol)", context);
  arma::mat Z(d, d, arma::fill::zeros);
  for (arma::uword i = 0; i < d; ++i) Z(i, i) = std::sqrt(R::rchisq(df - static_cast<double>(i)));
  for (arma::uword j = 1; j < d; ++j)
    for (arma::uword i = 0; i < j; ++i) Z(i, j) = R::norm_rand();
  IWishDraw out;
  out.U = arma::trimatu(Z * CC);
  out.W = out.U.t() * out.U;
  if (!arma::inv_sympd(out.Sigma, out.W))
    Rcpp::stop("%s: inverse-Wishart draw is singular", context);
  return out;
}

// Conjugate group update of the standard hierarchy without group designs:
// the joint normal draw of the means, an inverse-Wishart per covariance
// block, a gamma per unblocked variance, then the half-t scales. When `U` is
// given it receives each block's precision factor, in plan order.
void gibbs_standard(const arma::mat& alpha, const arma::mat& P0, const arma::vec& m0,
                    double v, const arma::vec& A_inv2, const BlockPlan& plan,
                    arma::vec& tmu, arma::mat& tvar, arma::mat& tvinv,
                    arma::vec& a_half, std::vector<arma::mat>* U = nullptr) {
  const arma::uword p = alpha.n_rows, n = alpha.n_cols;
  const arma::mat prec_post = P0 + static_cast<double>(n) * tvinv;
  const arma::vec b_post = P0 * m0 + tvinv * arma::sum(alpha, 1);
  arma::mat R_post;
  if (!arma::chol(R_post, prec_post))
    Rcpp::stop("group mean precision is not positive definite (chol)");
  const arma::vec w = arma::solve(arma::trimatl(R_post.t()), b_post);
  arma::vec e(p);
  for (arma::uword i = 0; i < p; ++i) e(i) = R::norm_rand();
  tmu = arma::solve(arma::trimatu(R_post), w + e);
  const arma::mat resid = alpha.each_col() - tmu;
  arma::mat tvar_new(p, p, arma::fill::zeros), tvinv_new(p, p, arma::fill::zeros);
  for (const std::vector<arma::uword>& b : plan.blocks) {
    const arma::uvec idx(b);
    const arma::uword d = idx.n_elem;
    const arma::mat rb = resid.rows(idx);
    arma::mat S = rb * rb.t();
    for (arma::uword k = 0; k < d; ++k) S(k, k) += 2.0 * v / a_half(idx(k));
    const double df = v + static_cast<double>(d) - 1.0 + static_cast<double>(n);
    IWishDraw draw = riwish(df, S, "standard Gibbs covariance block");
    tvar_new(idx, idx) = draw.Sigma;
    tvinv_new(idx, idx) = draw.W;
    if (U) U->push_back(std::move(draw.U));
  }
  const double shape_u = v / 2.0 + static_cast<double>(n) / 2.0;
  for (arma::uword i : plan.unblocked) {
    const double sum_sq = arma::dot(resid.row(i), resid.row(i));
    const double rate = v / a_half(i) + 0.5 * sum_sq;
    const double g = R::rgamma(shape_u, 1.0 / rate);
    tvinv_new(i, i) = g;
    tvar_new(i, i) = 1.0 / g;
  }
  for (arma::uword i = 0; i < p; ++i) {
    const double shape = (v + plan.block_dim(i)) / 2.0;
    const double rate = v * tvinv_new(i, i) + A_inv2(i);
    a_half(i) = 1.0 / R::rgamma(shape, 1.0 / rate);
  }
  tvar = std::move(tvar_new);
  tvinv = std::move(tvinv_new);
}

// Conjugate group update of the diagonal-gamma hierarchy: independent normal
// means when the prior precision is diagonal (a joint draw otherwise), then a
// gamma per precision.
void gibbs_diag_gamma(const arma::mat& alpha, const arma::mat& P0, const arma::vec& m0,
                      const arma::vec& shape0, const arma::vec& rate0,
                      arma::vec& tmu, arma::vec& tvinv) {
  const arma::uword p = alpha.n_rows, n = alpha.n_cols;
  const arma::vec P0d = P0.diag();
  double off = 0.0;
  for (arma::uword i = 0; i < p; ++i)
    for (arma::uword j = 0; j < p; ++j)
      if (i != j) off = std::max(off, std::abs(P0(i, j)));
  const arma::vec rs = arma::sum(alpha, 1);
  const double nn = static_cast<double>(n);
  if (off <= 1e-12 * arma::abs(P0d).max()) {
    for (arma::uword i = 0; i < p; ++i) {
      const double var_mu = 1.0 / (nn * tvinv(i) + P0d(i));
      const double mean_mu = var_mu * (rs(i) * tvinv(i) + m0(i) * P0d(i));
      const double sd = std::sqrt(var_mu);
      tmu(i) = sd == 0.0 ? mean_mu : mean_mu + sd * R::norm_rand();
    }
  } else {
    arma::mat prec = P0;
    prec.diag() += nn * tvinv;
    const arma::vec b = rs % tvinv + P0 * m0;
    arma::mat R_post;
    if (!arma::chol(R_post, prec))
      Rcpp::stop("group mean precision is not positive definite (chol)");
    arma::vec w = arma::solve(arma::trimatl(R_post.t()), b);
    for (arma::uword i = 0; i < p; ++i) w(i) += R::norm_rand();
    tmu = arma::solve(arma::trimatu(R_post), w);
  }
  const arma::mat resid = alpha.each_col() - tmu;
  for (arma::uword i = 0; i < p; ++i) {
    const double shape = shape0(i) + nn / 2.0;
    const double rate = rate0(i) + arma::dot(resid.row(i), resid.row(i)) / 2.0;
    tvinv(i) = R::rgamma(shape, 1.0 / rate);
  }
}

// The standard hierarchy's prior: theta_mu_invar, theta_mu_mean, v and A
// (one scale per parameter, or one for all).
struct StandardPrior {
  arma::mat P0;
  arma::vec m0, A_inv2;
  double v = 0.0;
};

StandardPrior standard_prior(const Rcpp::List& prior, arma::uword p) {
  StandardPrior pr;
  pr.P0 = Rcpp::as<arma::mat>(prior["theta_mu_invar"]);
  pr.m0 = Rcpp::as<arma::vec>(prior["theta_mu_mean"]);
  pr.v = Rcpp::as<double>(prior["v"]);
  const arma::vec A = Rcpp::as<arma::vec>(prior["A"]);
  if (pr.P0.n_rows != p || pr.P0.n_cols != p || pr.m0.n_elem != p ||
      (A.n_elem != 1 && A.n_elem != p))
    Rcpp::stop("standard Gibbs: prior dimension does not match the parameters");
  pr.A_inv2.set_size(p);
  for (arma::uword i = 0; i < p; ++i) {
    const double a = A(A.n_elem == 1 ? 0 : i);
    pr.A_inv2(i) = 1.0 / (a * a);
  }
  return pr;
}

}  // namespace

// The standard hierarchy's group update without group designs
// (gibbs_step_standard): means, covariance blocks, unblocked variances and
// half-t scales from the subject states and the previous precision and
// scales.
// [[Rcpp::export]]
Rcpp::List emc_gibbs_standard_cpp(const arma::mat& alpha, const arma::mat& tvinv0,
                                  const arma::vec& a_half0, const Rcpp::List& prior,
                                  const Rcpp::LogicalVector& is_blocked,
                                  const Rcpp::IntegerVector& par_group) {
  const arma::uword p = alpha.n_rows;
  if (tvinv0.n_rows != p || tvinv0.n_cols != p || a_half0.n_elem != p ||
      static_cast<arma::uword>(is_blocked.size()) != p ||
      static_cast<arma::uword>(par_group.size()) != p)
    Rcpp::stop("standard Gibbs: group state does not match the parameters");
  const StandardPrior pr = standard_prior(prior, p);
  const BlockPlan plan = block_plan(is_blocked, par_group);
  arma::vec tmu(p), a_half = a_half0;
  arma::mat tvar(p, p), tvinv = tvinv0;
  gibbs_standard(alpha, pr.P0, pr.m0, pr.v, pr.A_inv2, plan, tmu, tvar, tvinv, a_half);
  return Rcpp::List::create(Rcpp::Named("tmu") = as_vector(tmu), Rcpp::Named("tvar") = tvar,
                            Rcpp::Named("tvinv") = tvinv,
                            Rcpp::Named("a_half") = as_vector(a_half));
}

// The diagonal-gamma hierarchy's group update (gibbs_step_diag_gamma): means
// and precisions from the subject states and the previous precisions.
// [[Rcpp::export]]
Rcpp::List emc_gibbs_diag_gamma_cpp(const arma::mat& alpha, const arma::vec& tvinv0,
                                    const Rcpp::List& prior) {
  const arma::uword p = alpha.n_rows;
  const arma::mat P0 = Rcpp::as<arma::mat>(prior["theta_mu_invar"]);
  const arma::vec m0 = Rcpp::as<arma::vec>(prior["theta_mu_mean"]);
  const arma::vec shape0 = Rcpp::as<arma::vec>(prior["shape"]);
  const arma::vec rate0 = Rcpp::as<arma::vec>(prior["rate"]);
  if (tvinv0.n_elem != p || P0.n_rows != p || P0.n_cols != p || m0.n_elem != p ||
      shape0.n_elem != p || rate0.n_elem != p)
    Rcpp::stop("diagonal-gamma Gibbs: prior or state does not match the parameters");
  arma::vec tmu(p), tvinv = tvinv0;
  gibbs_diag_gamma(alpha, P0, m0, shape0, rate0, tmu, tvinv);
  return Rcpp::List::create(Rcpp::Named("tmu") = as_vector(tmu),
                            Rcpp::Named("tvinv") = as_vector(tvinv));
}

// Factors of the frozen pool proposals: per subject the lower Cholesky
// factor of V and log det V. The proposals change only between blocks, so a
// block computes these once.
// [[Rcpp::export]]
Rcpp::List emc_ensemble_factors_cpp(const Rcpp::List& var) {
  const int n = var.size();
  if (n < 1) Rcpp::stop("ensemble pools need at least one subject");
  const int p = Rcpp::as<arma::mat>(var[0]).n_rows;
  arma::cube L(p, p, n);
  arma::vec logdet(n);
  for (int s = 0; s < n; ++s) {
    const arma::mat V = Rcpp::as<arma::mat>(var[s]);
    if (static_cast<int>(V.n_rows) != p || static_cast<int>(V.n_cols) != p)
      Rcpp::stop("ensemble proposal for subject %d has the wrong dimension", s + 1);
    arma::mat Ls;
    if (!arma::chol(Ls, V, "lower"))
      Rcpp::stop("ensemble proposal covariance for subject %d is not positive definite (chol)", s + 1);
    L.slice(s) = Ls;
    logdet(s) = 2.0 * arma::sum(arma::log(Ls.diag()));
  }
  return Rcpp::List::create(Rcpp::Named("L") = L, Rcpp::Named("logdet") = as_vector(logdet));
}

// Fresh pool members for every subject from its frozen proposal
// (1 - w) N(mu, V) + w t_df(mu, c V), V = L L': per subject the normals of
// the M1 draws, then their component uniforms, then a chi-square per tail
// draw. Returns the draws (rows per subject), their log proposal densities,
// and the factors the current state's density needs later.
// [[Rcpp::export]]
Rcpp::List emc_ensemble_draw_pools_cpp(const Rcpp::List& mu, const arma::cube& L,
                                       const arma::vec& logdet, int M1, double tail_weight,
                                       double tail_df, double tail_excess,
                                       const Rcpp::CharacterVector& par_names) {
  const int n = mu.size();
  if (n < 1) Rcpp::stop("ensemble pools need at least one subject");
  const int p = L.n_rows;
  if (static_cast<int>(L.n_slices) != n || static_cast<int>(logdet.n_elem) != n)
    Rcpp::stop("ensemble pools: proposal factors do not match the subjects");
  const double c = 1.0 + tail_excess / std::sqrt(static_cast<double>(p));
  const LogQ log_q(p, c, tail_weight, tail_df);
  Rcpp::List rows(n), lq(n);
  arma::mat mus(p, n);
  arma::mat z(p, M1), X(M1, p);
  arma::vec scale(M1), lqs(M1);
  for (int s = 0; s < n; ++s) {
    const arma::vec m = Rcpp::as<arma::vec>(mu[s]);
    if (static_cast<int>(m.n_elem) != p)
      Rcpp::stop("ensemble proposal for subject %d has the wrong dimension", s + 1);
    mus.col(s) = m;
    for (int j = 0; j < M1; ++j)
      for (int i = 0; i < p; ++i) z(i, j) = R::norm_rand();
    for (int j = 0; j < M1; ++j) scale(j) = R::runif(0.0, 1.0) < tail_weight ? -1.0 : 1.0;
    for (int j = 0; j < M1; ++j)
      if (scale(j) < 0) scale(j) = std::sqrt(c * tail_df / R::rchisq(tail_df));
    const arma::mat Lz = L.slice(s) * z;
    for (int j = 0; j < M1; ++j) {
      double m2 = 0.0;
      for (int i = 0; i < p; ++i) {
        X(j, i) = m(i) + Lz(i, j) * scale(j);
        m2 += z(i, j) * z(i, j);
      }
      lqs(j) = log_q(scale(j) * scale(j) * m2, logdet(s));
    }
    Rcpp::NumericMatrix Xr = Rcpp::wrap(X);
    Rcpp::colnames(Xr) = par_names;
    rows[s] = Xr;
    lq[s] = as_vector(lqs);
  }
  return Rcpp::List::create(Rcpp::Named("rows") = rows, Rcpp::Named("lq") = lq,
                            Rcpp::Named("mu") = mus, Rcpp::Named("c") = c);
}

// The stacked pools: subject s occupies rows (s - 1) M + 1 .. s M, its current
// state at position k0[s] among the M - 1 fresh draws. Returns the stack,
// the pool log-likelihoods and the group-independent part of the log weights,
// ll - log q, with the current state's log q computed here.
// [[Rcpp::export]]
Rcpp::List emc_ensemble_stack_cpp(const Rcpp::List& rows, const Rcpp::List& lq,
                                  const arma::mat& ll_new, const arma::mat& proposals,
                                  const Rcpp::IntegerVector& k0, const arma::cube& L,
                                  const arma::mat& mu, const arma::vec& logdet, double c,
                                  double tail_weight, double tail_df) {
  const int n = rows.size();
  const int M1 = ll_new.n_rows, M = M1 + 1;
  const int p = mu.n_rows;
  if (static_cast<int>(proposals.n_rows) != p + 1 || static_cast<int>(proposals.n_cols) != n ||
      static_cast<int>(ll_new.n_cols) != n || k0.size() != n)
    Rcpp::stop("ensemble stack: inconsistent dimensions");
  const LogQ log_q(p, c, tail_weight, tail_df);
  arma::mat stack(n * M, p);
  arma::vec ll(n * M), base_w(n * M);
  for (int s = 0; s < n; ++s) {
    const arma::mat fresh = Rcpp::as<arma::mat>(rows[s]);
    const arma::vec lq_s = Rcpp::as<arma::vec>(lq[s]);
    if (static_cast<int>(fresh.n_rows) != M1 || static_cast<int>(fresh.n_cols) != p ||
        static_cast<int>(lq_s.n_elem) != M1)
      Rcpp::stop("ensemble stack: pool of subject %d has the wrong shape", s + 1);
    const int kc = k0[s] - 1;
    if (kc < 0 || kc >= M) Rcpp::stop("ensemble stack: position out of range");
    const int off = s * M;
    const arma::vec x = proposals.col(s).head(p);
    const arma::vec zc = arma::solve(arma::trimatl(L.slice(s)), x - mu.col(s));
    const double lq_cur = log_q(arma::dot(zc, zc), logdet(s));
    stack.row(off + kc) = x.t();
    ll(off + kc) = proposals(p, s);
    base_w(off + kc) = ll(off + kc) - lq_cur;
    for (int j = 0, r = 0; j < M; ++j) {
      if (j == kc) continue;
      stack.row(off + j) = fresh.row(r);
      ll(off + j) = ll_new(r, s);
      base_w(off + j) = ll_new(r, s) - lq_s(r);
      ++r;
    }
  }
  return Rcpp::List::create(Rcpp::Named("pool_stack") = stack,
                            Rcpp::Named("ll_stack") = as_vector(ll),
                            Rcpp::Named("base_w") = as_vector(base_w));
}

// The Gibbs rounds on pi~(h, k | A): the group state given the selected
// subject states (the hierarchy's conjugate update), then every subject's
// pool position given the group state (a categorical draw by the Gumbel-max
// trick, one uniform per pool row in stack order). type 0 is the standard
// hierarchy (state: tvinv, a_half; prior: theta_mu_invar, theta_mu_mean, v,
// A), type 1 the diagonal-gamma hierarchy (state: tvinv as a vector; prior:
// theta_mu_invar, theta_mu_mean, shape, rate). Returns the final selection,
// the group state and the pools' weight ESS under it.
// [[Rcpp::export]]
Rcpp::List emc_ensemble_rounds_cpp(const arma::mat& pool_stack, const arma::vec& base_w,
                                   const Rcpp::IntegerVector& k0, int M, int rounds,
                                   int type, const Rcpp::List& state, const Rcpp::List& prior,
                                   const Rcpp::LogicalVector& is_blocked,
                                   const Rcpp::IntegerVector& par_group) {
  const int n = k0.size();
  const arma::uword p = pool_stack.n_cols;
  const arma::uword nM = pool_stack.n_rows;
  if (nM != static_cast<arma::uword>(n) * M || base_w.n_elem != nM)
    Rcpp::stop("ensemble rounds: inconsistent dimensions");
  const arma::mat P0 = Rcpp::as<arma::mat>(prior["theta_mu_invar"]);
  const arma::vec m0 = Rcpp::as<arma::vec>(prior["theta_mu_mean"]);
  if (P0.n_rows != p || P0.n_cols != p || m0.n_elem != p)
    Rcpp::stop("ensemble rounds: prior dimension does not match the parameters");
  arma::vec tmu(p, arma::fill::zeros), a_half, tvinv_d, shape0, rate0, A_inv2;
  arma::mat tvar, tvinv;
  double v = 0.0;
  BlockPlan plan;
  if (type == 0) {
    tvinv = Rcpp::as<arma::mat>(state["tvinv"]);
    a_half = Rcpp::as<arma::vec>(state["a_half"]);
    const StandardPrior pr = standard_prior(prior, p);
    v = pr.v;
    A_inv2 = pr.A_inv2;
    if (tvinv.n_rows != p || tvinv.n_cols != p || a_half.n_elem != p ||
        static_cast<arma::uword>(is_blocked.size()) != p ||
        static_cast<arma::uword>(par_group.size()) != p)
      Rcpp::stop("ensemble rounds: group state does not match the parameters");
    plan = block_plan(is_blocked, par_group);
    tvar.set_size(p, p);
  } else if (type == 1) {
    tvinv_d = Rcpp::as<arma::vec>(state["tvinv"]);
    shape0 = Rcpp::as<arma::vec>(prior["shape"]);
    rate0 = Rcpp::as<arma::vec>(prior["rate"]);
    if (tvinv_d.n_elem != p || shape0.n_elem != p || rate0.n_elem != p)
      Rcpp::stop("ensemble rounds: group state does not match the parameters");
  } else {
    Rcpp::stop("ensemble rounds: unsupported hierarchy");
  }
  std::vector<int> selected(n);
  for (int s = 0; s < n; ++s) selected[s] = k0[s] - 1;
  arma::mat alpha(p, n);
  arma::vec lw(nM), wess(n, arma::fill::zeros);
  const arma::mat stack_t = pool_stack.t();
  arma::mat D(p, nM);
  arma::vec m2(nM);
  std::vector<arma::mat> U;
  for (int r = 1; r <= rounds; ++r) {
    for (int s = 0; s < n; ++s) alpha.col(s) = stack_t.col(s * M + selected[s]);
    if (type == 0) {
      U.clear();
      gibbs_standard(alpha, P0, m0, v, A_inv2, plan, tmu, tvar, tvinv, a_half, &U);
      // The group density through the block precision factors: block b adds
      // |U_b (x_b - mu_b)|^2 and -2 sum log diag U_b, an unblocked parameter
      // its precision times the squared deviation and -log precision.
      D = stack_t;
      D.each_col() -= tmu;
      m2.zeros();
      double logdet_var = 0.0;
      for (std::size_t b = 0; b < plan.blocks.size(); ++b) {
        const arma::uvec idx(plan.blocks[b]);
        const arma::mat Zb = U[b] * D.rows(idx);
        m2 += arma::sum(arma::square(Zb), 0).t();
        logdet_var -= 2.0 * arma::sum(arma::log(arma::abs(U[b].diag())));
      }
      for (arma::uword i : plan.unblocked) {
        m2 += tvinv(i, i) * arma::square(D.row(i)).t();
        logdet_var -= std::log(tvinv(i, i));
      }
      const double norm = p * LOG_2PI + logdet_var;
      for (arma::uword i = 0; i < nM; ++i) lw(i) = base_w(i) - 0.5 * (norm + m2(i));
    } else {
      gibbs_diag_gamma(alpha, P0, m0, shape0, rate0, tmu, tvinv_d);
      const arma::vec dvar = 1.0 / tvinv_d;
      const arma::vec sd = arma::sqrt(dvar);
      const double norm = p * LOG_2PI + arma::sum(arma::log(dvar));
      for (arma::uword i = 0; i < nM; ++i) {
        double q = 0.0;
        for (arma::uword k = 0; k < p; ++k) {
          const double zk = (stack_t(k, i) - tmu(k)) / sd(k);
          q += zk * zk;
        }
        lw(i) = base_w(i) - 0.5 * (norm + q);
      }
    }
    for (int s = 0; s < n; ++s) {
      const int off = s * M;
      int best = 0;
      double best_v = -std::numeric_limits<double>::infinity();
      for (int j = 0; j < M; ++j) {
        const double g = -std::log(-std::log(R::runif(0.0, 1.0)));
        const double val = lw(off + j) + g;
        if (j == 0 || val > best_v) { best_v = val; best = j; }
      }
      selected[s] = best;
    }
    if (r == rounds) {
      for (int s = 0; s < n; ++s) {
        const int off = s * M;
        const double mx = lw.subvec(off, off + M - 1).max();
        double sw = 0.0, sw2 = 0.0;
        for (int j = 0; j < M; ++j) {
          const double w = std::exp(lw(off + j) - mx);
          sw += w;
          sw2 += w * w;
        }
        wess(s) = sw * sw / sw2;
      }
    }
  }
  Rcpp::IntegerVector sel(n);
  for (int s = 0; s < n; ++s) sel[s] = selected[s] + 1;
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("selected") = sel,
                                      Rcpp::Named("tmu") = as_vector(tmu),
                                      Rcpp::Named("wess") = as_vector(wess));
  if (type == 0) {
    out["tvar"] = tvar;
    out["tvinv"] = tvinv;
    out["a_half"] = as_vector(a_half);
  } else {
    out["tvinv"] = as_vector(tvinv_d);
  }
  return out;
}

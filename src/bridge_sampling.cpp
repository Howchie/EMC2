#include <RcppArmadillo.h>
#include <vector>
#include "utility_functions.h"
#include "composite_functions.h"

namespace {

// log of 2^(v p / 2) Gamma_p(v / 2), the inverse-Wishart normaliser without
// its scale determinant.
double log_iwish_norm(double v, arma::uword p) {
  double out = 0.5 * v * p * M_LN2 + 0.25 * p * (p - 1) * std::log(M_PI);
  for (arma::uword j = 0; j < p; ++j) out += R::lgammafn(0.5 * (v - j));
  return out;
}

// Assembles log IW(W | v, S) from 0.5 log|S|, 0.5 log|W| and tr(S W^-1).
// A term that is not finite means W or S lies off the positive-definite
// cone, where the density is zero.
double log_iwish_parts(double v, arma::uword p, double half_logdet_S,
                       double half_logdet_W, double trace) {
  const double out = v * half_logdet_S - (v + p + 1) * half_logdet_W -
    0.5 * trace - log_iwish_norm(v, p);
  return R_FINITE(out) ? out : R_NegInf;
}

// Lower Cholesky factor of the lower triangle of A; false when A is not
// numerically positive definite.
bool lower_chol(const arma::mat& A, arma::mat& L) {
  const arma::uword p = A.n_rows;
  L.zeros(p, p);
  for (arma::uword j = 0; j < p; ++j) {
    double pivot = A(j, j);
    for (arma::uword k = 0; k < j; ++k) pivot -= L(j, k) * L(j, k);
    if (!R_FINITE(pivot) || pivot <= 0) return false;
    const double diag = std::sqrt(pivot);
    L(j, j) = diag;
    for (arma::uword i = j + 1; i < p; ++i) {
      double sum = A(i, j);
      for (arma::uword k = 0; k < j; ++k) sum -= L(i, k) * L(j, k);
      L(i, j) = sum / diag;
    }
  }
  return true;
}

// Fills the lower factor L from a column-major lower triangle whose diagonal
// is on the log scale and returns the sum of the log diagonal.
double unpack_log_chol(const double* theta, arma::uword p, arma::mat& L) {
  L.zeros(p, p);
  double log_diag = 0;
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = j; i < p; ++i) {
      const double x = *theta++;
      if (i == j) {
        log_diag += x;
        L(i, j) = std::exp(x);
      } else {
        L(i, j) = x;
      }
    }
  }
  return log_diag;
}

// sum_j scale[j] * |column j of L^-1|^2, which is tr(S (L L')^-1) for
// S = diag(scale).
double trace_scale_inverse(const arma::mat& L, const double* scale,
                           std::vector<double>& column) {
  const arma::uword p = L.n_rows;
  column.resize(p);
  double trace = 0;
  for (arma::uword j = 0; j < p; ++j) {
    double norm2 = 0;
    for (arma::uword i = j; i < p; ++i) {
      double sum = i == j ? 1.0 : 0.0;
      for (arma::uword k = j; k < i; ++k) sum -= L(i, k) * column[k];
      column[i] = sum / L(i, i);
      norm2 += column[i] * column[i];
    }
    trace += scale[j] * norm2;
  }
  return trace;
}

// log IW(W | v, diag(exp(log_scale))) with W = L L' (covariance = true) or
// W = (L L')^-1 (covariance = false), for a lower factor L whose log
// diagonal sums to log_diag.
double log_iwish_factor(const arma::mat& L, double log_diag, double v,
                        const double* log_scale, bool covariance,
                        std::vector<double>& scale, std::vector<double>& column) {
  const arma::uword p = L.n_rows;
  scale.resize(p);
  double half_logdet_S = 0;
  for (arma::uword j = 0; j < p; ++j) {
    half_logdet_S += 0.5 * log_scale[j];
    scale[j] = std::exp(log_scale[j]);
  }
  double trace = 0;
  if (covariance) {
    trace = trace_scale_inverse(L, scale.data(), column);
  } else {
    for (arma::uword j = 0; j < p; ++j)
      for (arma::uword i = j; i < p; ++i) trace += scale[i] * L(i, j) * L(i, j);
  }
  return log_iwish_parts(v, p, half_logdet_S, covariance ? log_diag : -log_diag, trace);
}

}  // namespace

// Log inverse-Wishart density of W with v degrees of freedom and scale S,
// from the Cholesky factors of both. -Inf when W or S is not positive
// definite.
// [[Rcpp::export(rng = false)]]
double log_diwish(const arma::mat& W, double v, const arma::mat& S) {
  const arma::uword p = W.n_rows;
  if (!p || W.n_cols != p || S.n_rows != p || S.n_cols != p)
    Rcpp::stop("W and S must be square matrices of the same size");
  if (!R_FINITE(v) || v <= p - 1)
    Rcpp::stop("The degrees of freedom must exceed the dimension minus one");
  arma::mat LW, LS;
  if (!lower_chol(W, LW) || !lower_chol(S, LS)) return R_NegInf;
  double half_logdet_W = 0, half_logdet_S = 0, trace = 0;
  for (arma::uword j = 0; j < p; ++j) {
    half_logdet_W += std::log(LW(j, j));
    half_logdet_S += std::log(LS(j, j));
  }
  // tr(S W^-1) is the squared Frobenius norm of LW^-1 LS.
  std::vector<double> column(p);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = j; i < p; ++i) {
      double sum = LS(i, j);
      for (arma::uword k = j; k < i; ++k) sum -= LW(i, k) * column[k];
      column[i] = sum / LW(i, i);
      trace += column[i] * column[i];
    }
  }
  return log_iwish_parts(v, p, half_logdet_S, half_logdet_W, trace);
}

// Log inverse-Wishart density with scale diag(exp(log_scale)) of the matrix
// whose log-Cholesky vector is theta (column-major lower triangle of L, log
// diagonal). The matrix is L L' when covariance is true and (L L')^-1 when it
// is false; no factorisation or inverse of the matrix itself is formed.
// [[Rcpp::export(rng = false)]]
double log_diwish_log_chol(const Rcpp::NumericVector& theta, double v,
                           const Rcpp::NumericVector& log_scale,
                           bool covariance = true) {
  const arma::uword p = log_scale.size();
  if (!p || static_cast<arma::uword>(theta.size()) != p * (p + 1) / 2)
    Rcpp::stop("theta must hold the lower triangle of a matrix matching log_scale");
  if (!R_FINITE(v) || v <= p - 1)
    Rcpp::stop("The degrees of freedom must exceed the dimension minus one");
  arma::mat L;
  std::vector<double> scale, column;
  const double log_diag = unpack_log_chol(theta.begin(), p, L);
  return log_iwish_factor(L, log_diag, v, log_scale.begin(), covariance, scale, column);
}

// Group likelihood, group prior and Jacobian of the standard hierarchy for
// each proposal row. `alpha` holds the subject parameters, one subject after
// another. `group` holds theta_mu, log a_half, the log variances of the
// parameters in `unblocked`, then one log-Cholesky vector per entry of
// `blocks`; both index parameters from zero. The covariance factors are read
// off the proposal, so nothing is factorised per row.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector bridge_group_standard(
    const Rcpp::NumericMatrix& alpha, const Rcpp::NumericMatrix& group,
    const Rcpp::List& group_designs, const Rcpp::IntegerVector& unblocked,
    const Rcpp::List& blocks, const arma::vec& mu_mean, const arma::mat& mu_var,
    double v, const Rcpp::NumericVector& A) {
  const arma::uword p = A.size();
  const R_xlen_t n = alpha.nrow();
  if (!p || alpha.ncol() % p || group.nrow() != n ||
      static_cast<arma::uword>(group_designs.size()) != p)
    Rcpp::stop("Proposals do not match the number of parameters");
  const arma::uword n_subj = alpha.ncol() / p;

  std::vector<arma::mat> designs(p);
  std::vector<arma::uword> first(p);
  arma::uword n_mu = 0;
  for (arma::uword k = 0; k < p; ++k) {
    designs[k] = Rcpp::as<arma::mat>(group_designs[k]);
    if (designs[k].n_rows != n_subj) Rcpp::stop("A group design does not match the subjects");
    first[k] = n_mu;
    n_mu += designs[k].n_cols;
  }
  std::vector<std::vector<arma::uword>> block_pars(blocks.size());
  arma::uword n_cols = n_mu + p + unblocked.size(), max_block = 1;
  for (R_xlen_t b = 0; b < blocks.size(); ++b) {
    block_pars[b] = Rcpp::as<std::vector<arma::uword>>(blocks[b]);
    const arma::uword d = block_pars[b].size();
    for (arma::uword par : block_pars[b])
      if (par >= p) Rcpp::stop("A covariance block indexes a missing parameter");
    n_cols += d * (d + 1) / 2;
    max_block = std::max(max_block, d);
  }
  for (int par : unblocked)
    if (par < 0 || static_cast<arma::uword>(par) >= p)
      Rcpp::stop("An unblocked index is out of range");
  if (static_cast<arma::uword>(group.ncol()) != n_cols || mu_mean.n_elem != n_mu ||
      mu_var.n_rows != n_mu || mu_var.n_cols != n_mu)
    Rcpp::stop("The group proposals do not match the hierarchy");

  arma::mat mu_chol;
  const bool mu_proper = lower_chol(mu_var, mu_chol);
  double mu_half_logdet = 0;
  if (mu_proper)
    for (arma::uword j = 0; j < n_mu; ++j) mu_half_logdet += std::log(mu_chol(j, j));
  const double log_2pi = std::log(2 * M_PI), log_v = std::log(v);
  const double lgamma_half_v = R::lgammafn(0.5 * v), lgamma_half = R::lgammafn(0.5);
  const double constant = -0.5 * log_2pi * (n_mu + static_cast<double>(n_subj) * p) -
    mu_half_logdet;

  arma::mat resid(p, n_subj), L;
  std::vector<double> row(n_cols), z(std::max(n_mu, max_block)), log_scale(max_block),
    scale, column;
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < n_cols; ++j) row[j] = group(i, j);
    const double* theta_mu = row.data();
    const double* theta_a = theta_mu + n_mu;
    const double* theta_var = theta_a + p;
    const double* theta_chol = theta_var + unblocked.size();
    double total = constant;

    // Normal prior on the group mean coefficients.
    if (!mu_proper) total = R_NegInf;
    for (arma::uword t = 0; mu_proper && t < n_mu; ++t) {
      double sum = theta_mu[t] - mu_mean[t];
      for (arma::uword k = 0; k < t; ++k) sum -= mu_chol(t, k) * z[k];
      z[t] = sum / mu_chol(t, t);
      total -= 0.5 * z[t] * z[t];
    }

    for (arma::uword k = 0; k < p; ++k) {
      const arma::mat& X = designs[k];
      for (arma::uword s = 0; s < n_subj; ++s) {
        double mean = 0;
        for (arma::uword j = 0; j < X.n_cols; ++j) mean += X(s, j) * theta_mu[first[k] + j];
        resid(k, s) = alpha(i, s * p + k) - mean;
      }
      // a_half is inverse gamma with shape 1/2 and rate 1 / A^2; its log is sampled.
      total += -std::log(A[k]) - lgamma_half - 0.5 * theta_a[k] -
        std::exp(-theta_a[k]) / (A[k] * A[k]);
    }

    for (R_xlen_t u = 0; u < unblocked.size(); ++u) {
      const arma::uword k = unblocked[u];
      const double log_var = theta_var[u];
      double ss = 0;
      for (arma::uword s = 0; s < n_subj; ++s) ss += resid(k, s) * resid(k, s);
      total -= 0.5 * (n_subj * log_var + std::exp(-log_var) * ss);
      // The variance is inverse gamma with shape v / 2 and rate v / a_half;
      // its log is sampled.
      total += 0.5 * v * (log_v - theta_a[k]) - lgamma_half_v - 0.5 * v * log_var -
        v * std::exp(-theta_a[k] - log_var);
    }

    for (const std::vector<arma::uword>& pars : block_pars) {
      const arma::uword d = pars.size();
      const double log_diag = unpack_log_chol(theta_chol, d, L);
      double ss = 0;
      for (arma::uword s = 0; s < n_subj; ++s) {
        for (arma::uword t = 0; t < d; ++t) {
          double sum = resid(pars[t], s);
          for (arma::uword k = 0; k < t; ++k) sum -= L(t, k) * z[k];
          z[t] = sum / L(t, t);
          ss += z[t] * z[t];
        }
      }
      total -= n_subj * log_diag + 0.5 * ss;
      // Inverse-Wishart prior with scale 2 v diag(1 / a_half), and the
      // Jacobian of the covariance with respect to its log-Cholesky vector.
      total += d * M_LN2;
      for (arma::uword t = 0; t < d; ++t) {
        log_scale[t] = M_LN2 + log_v - theta_a[pars[t]];
        total += (d - t + 1) * *theta_chol;
        theta_chol += d - t;
      }
      total += log_iwish_factor(L, log_diag, v + d - 1, log_scale.data(), true,
                                scale, column);
    }
    out[i] = R_FINITE(total) ? total : R_NegInf;
  }
  return out;
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

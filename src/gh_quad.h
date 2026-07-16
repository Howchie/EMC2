#ifndef EMC2_GH_QUAD_H
#define EMC2_GH_QUAD_H

// Runtime Gauss-Hermite rules for integrals against a standard normal latent
// factor.  The Hermite rule generated here integrates exp(-x^2); callers use
// z = sqrt(2) * x and divide weights by sqrt(pi).
//
// BAwLcorr integrates its shared latent factor with two batched passes of
// these rules: a fixed scan pass locates each trial's integrand mass, and a
// second pass re-centers a rule per trial (agh_center_from_scan /
// agh_log_weight below).  See c_log_likelihood_bawl_correlated.
//
// The repository carries only the adaptive subset of bundled GSL, so its
// gsl_integration_fixed Hermite implementation is not linkable here.  The
// equivalent Golub-Welsch construction uses the Hermite Jacobi matrix and
// caches the resulting rules by node count.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

struct GHRule {
  std::vector<double> x;
  std::vector<double> w;
};

inline GHRule make_gh_rule(int n) {
  if (n < 2 || n > 256)
    throw std::runtime_error("BAwLcorr GH node count must be between 2 and 256.");

  arma::mat jacobi(static_cast<arma::uword>(n), static_cast<arma::uword>(n),
                   arma::fill::zeros);
  for (int i = 0; i < n - 1; ++i) {
    const double offdiag = std::sqrt((static_cast<double>(i) + 1.0) / 2.0);
    jacobi(i, i + 1) = offdiag;
    jacobi(i + 1, i) = offdiag;
  }

  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, jacobi))
    throw std::runtime_error("BAwLcorr failed to construct a Gauss-Hermite rule.");

  GHRule rule;
  rule.x.assign(eigval.begin(), eigval.end());
  rule.w.resize(static_cast<size_t>(n));
  const double sqrt_pi = std::sqrt(std::acos(-1.0));
  for (int i = 0; i < n; ++i) {
    // The first eigenvector component squared is the normalized weight.
    rule.w[static_cast<size_t>(i)] = sqrt_pi * eigvec(0, i) * eigvec(0, i);
  }
  return rule;
}

inline const GHRule& gh_rule(int n) {
  static std::unordered_map<int, GHRule> cache;
  static std::mutex cache_mutex;
  std::lock_guard<std::mutex> lock(cache_mutex);
  auto it = cache.find(n);
  if (it == cache.end()) it = cache.emplace(n, make_gh_rule(n)).first;
  return it->second;
}

inline double gh_standard_normal_weight(const GHRule& rule, int i) {
  return rule.w[static_cast<size_t>(i)] / std::sqrt(std::acos(-1.0));
}

inline int bawl_corr_quad_nodes(const char* env_name, int default_n) {
  const char* raw = std::getenv(env_name);
  if (raw == nullptr || *raw == '\0') return default_n;
  char* end = nullptr;
  const long parsed = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || parsed < 2 || parsed > 256)
    throw std::runtime_error("BAwLcorr quadrature node count must be between 2 and 256.");
  return static_cast<int>(parsed);
}

// Per-trial center/scale for adaptive (moment-matched) Gauss-Hermite.
// Estimated from log node masses g (weights included) at latent values z; a
// second pass re-centers its rule at N(mu, sigma^2).  sigma is inflated and
// clamped so a peak mislocated by up to the scan-node spacing stays inside
// the refined rule's span.
struct AGHCenter {
  double mu = 0.0;
  double sigma = 1.0;
};

inline AGHCenter agh_center_from_scan(const double* g, const double* z, int n,
                                      double sigma_min = 0.3,
                                      double sigma_max = 1.25,
                                      double inflate = 1.3) {
  AGHCenter out;
  double m = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) m = std::max(m, g[i]);
  if (!std::isfinite(m)) return out;  // nothing seen: keep the prior rule
  double s0 = 0.0, s1 = 0.0, s2 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(g[i])) continue;
    const double p = std::exp(g[i] - m);
    s0 += p;
    s1 += p * z[i];
    s2 += p * z[i] * z[i];
  }
  if (!(s0 > 0.0)) return out;
  const double mu = s1 / s0;
  const double var = std::fmax(s2 / s0 - mu * mu, 0.0);
  if (!std::isfinite(mu)) return out;
  out.mu = mu;
  out.sigma = std::fmin(std::fmax(inflate * std::sqrt(var), sigma_min), sigma_max);
  return out;
}

// Log quadrature weight for node i of a rule adapted to N(mu, sigma^2):
// int f(z) phi(z) dz ~= sum_i sigma sqrt(2) w_i e^{x_i^2} phi(z_i) f(z_i)
// with z_i = mu + sigma sqrt(2) x_i.  `lgw` is log(w_i) + x_i^2.  At
// (mu, sigma) = (0, 1) this reduces exactly to gh_standard_normal_weight.
inline double agh_log_weight(double lgw, double sigma, double z) {
  const double log_sqrt2 = 0.5 * std::log(2.0);
  const double log_sqrt_2pi = 0.5 * std::log(2.0 * std::acos(-1.0));
  return lgw + std::log(sigma) + log_sqrt2 - log_sqrt_2pi - 0.5 * z * z;
}

#endif

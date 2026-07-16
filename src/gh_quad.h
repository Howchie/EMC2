#ifndef EMC2_GH_QUAD_H
#define EMC2_GH_QUAD_H

// Runtime Gauss-Hermite rules for integrals against a standard normal latent
// factor.  The Hermite rule generated here integrates exp(-x^2); callers use
// z = sqrt(2) * x and divide weights by sqrt(pi).
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

inline int gh_node_count(double max_abs_rho = 0.0) {
  const char* raw = std::getenv("EMC2_BAWLCORR_GH_N");
  // The convergence sweep on positive-drift trials shows that 40 nodes is
  // accurate through |rho| <= .8, while the residual-SD collapse near .95
  // needs a denser rule.  Keep the override for cheap studies and debugging;
  // production calls adapt only at the high-correlation boundary.
  if (raw == nullptr || *raw == '\0') {
    if (max_abs_rho > 0.9) return 200;
    if (max_abs_rho > 0.8) return 80;
    return 40;
  }
  char* end = nullptr;
  const long parsed = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || parsed < 2 || parsed > 256)
    throw std::runtime_error("EMC2_BAWLCORR_GH_N must be an integer between 2 and 256.");
  return static_cast<int>(parsed);
}

inline double gh_standard_normal_weight(const GHRule& rule, int i) {
  return rule.w[static_cast<size_t>(i)] / std::sqrt(std::acos(-1.0));
}

#endif

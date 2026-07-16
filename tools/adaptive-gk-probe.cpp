// Small, isolated Rcpp probe for src/adaptive_gk.h.
//
// This file is intentionally outside src/: it exercises the candidate
// quadrature without changing any production likelihood path.

// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <vector>

#include "../src/adaptive_gk.h"

namespace {

Rcpp::List pack_result(const AdaptiveLogVectorResult& result,
                       std::size_t callback_points) {
  return Rcpp::List::create(
      Rcpp::Named("log_integral") =
        Rcpp::NumericVector(result.log_integral.begin(), result.log_integral.end()),
      Rcpp::Named("log_error") =
        Rcpp::NumericVector(result.log_error.begin(), result.log_error.end()),
      Rcpp::Named("intervals") = static_cast<double>(result.intervals),
      Rcpp::Named("eval_points") = static_cast<double>(result.eval_points),
      Rcpp::Named("callback_points") = static_cast<double>(callback_points),
      Rcpp::Named("converged") = result.converged);
}

void check_same_length(const Rcpp::NumericVector& a,
                       const Rcpp::NumericVector& b,
                       const char* message) {
  if (a.size() != b.size()) Rcpp::stop(message);
}

} // namespace

// Integrate exp(amplitude - .5 * ((x - center) / width)^2) over [0, 1].
// [[Rcpp::export]]
Rcpp::List adaptive_gk_gaussian(Rcpp::NumericVector centers,
                                Rcpp::NumericVector widths,
                                Rcpp::NumericVector amplitudes,
                                double abs_tol,
                                double rel_tol,
                                int limit,
                                Rcpp::NumericVector initial_breaks) {
  check_same_length(centers, widths, "centers and widths must have the same length");
  check_same_length(centers, amplitudes, "centers and amplitudes must have the same length");

  const std::size_t n = static_cast<std::size_t>(centers.size());
  std::vector<double> breaks(initial_breaks.begin(), initial_breaks.end());
  std::size_t callback_points = 0;
  auto evaluator = [&](double x, std::vector<double>& out) {
    ++callback_points;
    for (std::size_t j = 0; j < n; ++j) {
      const double width = widths[static_cast<R_xlen_t>(j)];
      if (!(width > 0.0)) {
        out[j] = -std::numeric_limits<double>::infinity();
        continue;
      }
      const double delta = (x - centers[static_cast<R_xlen_t>(j)]) / width;
      out[j] = amplitudes[static_cast<R_xlen_t>(j)] - 0.5 * delta * delta;
    }
  };

  const AdaptiveLogVectorResult result = adaptive_gk21_log_vector(
      n, evaluator, abs_tol, rel_tol,
      static_cast<std::size_t>(std::max(limit, 0)), breaks);
  return pack_result(result, callback_points);
}

// Integrate a conditional normal observation over a standard-normal latent
// factor.  The caller supplies x = Phi(z), so the normal prior density is
// already represented by dx and must not be included in the log integrand.
// [[Rcpp::export]]
Rcpp::List adaptive_gk_latent_normal(Rcpp::NumericVector means,
                                     Rcpp::NumericVector loadings,
                                     Rcpp::NumericVector residual_sd,
                                     Rcpp::NumericVector observations,
                                     double abs_tol,
                                     double rel_tol,
                                     int limit,
                                     Rcpp::NumericVector initial_breaks) {
  check_same_length(means, loadings, "means and loadings must have the same length");
  check_same_length(means, residual_sd, "means and residual_sd must have the same length");
  check_same_length(means, observations, "means and observations must have the same length");

  const std::size_t n = static_cast<std::size_t>(means.size());
  std::vector<double> breaks(initial_breaks.begin(), initial_breaks.end());
  std::size_t callback_points = 0;
  const double log_sqrt_2pi = 0.5 * std::log(2.0 * std::acos(-1.0));
  auto evaluator = [&](double x, std::vector<double>& out) {
    ++callback_points;
    const double z = R::qnorm5(x, 0.0, 1.0, true, false);
    for (std::size_t j = 0; j < n; ++j) {
      const double sd = residual_sd[static_cast<R_xlen_t>(j)];
      if (!(sd > 0.0) || !R_FINITE(z)) {
        out[j] = -std::numeric_limits<double>::infinity();
        continue;
      }
      const double delta = (observations[static_cast<R_xlen_t>(j)] -
                            means[static_cast<R_xlen_t>(j)] -
                            loadings[static_cast<R_xlen_t>(j)] * z) / sd;
      out[j] = -log_sqrt_2pi - std::log(sd) - 0.5 * delta * delta;
    }
  };

  const AdaptiveLogVectorResult result = adaptive_gk21_log_vector(
      n, evaluator, abs_tol, rel_tol,
      static_cast<std::size_t>(std::max(limit, 0)), breaks);
  return pack_result(result, callback_points);
}

// Fixed Gauss-Hermite reference/cost comparator for the same latent-normal
// integrand. `nodes` and `weights` are the statmod/GSL convention for
// integral exp(-x^2) f(x) dx.
// [[Rcpp::export]]
Rcpp::NumericVector gh_latent_normal(Rcpp::NumericVector means,
                                     Rcpp::NumericVector loadings,
                                     Rcpp::NumericVector residual_sd,
                                     Rcpp::NumericVector observations,
                                     Rcpp::NumericVector nodes,
                                     Rcpp::NumericVector weights) {
  check_same_length(means, loadings, "means and loadings must have the same length");
  check_same_length(means, residual_sd, "means and residual_sd must have the same length");
  check_same_length(means, observations, "means and observations must have the same length");
  if (nodes.size() != weights.size())
    Rcpp::stop("nodes and weights must have the same length");

  const std::size_t n = static_cast<std::size_t>(means.size());
  const double log_sqrt_pi = 0.5 * std::log(std::acos(-1.0));
  const double log_sqrt_2pi = 0.5 * std::log(2.0 * std::acos(-1.0));
  std::vector<double> out(n, -std::numeric_limits<double>::infinity());
  for (R_xlen_t q = 0; q < nodes.size(); ++q) {
    const double z = std::sqrt(2.0) * nodes[q];
    const double log_w = std::log(weights[q]) - log_sqrt_pi;
    for (std::size_t j = 0; j < n; ++j) {
      const double sd = residual_sd[static_cast<R_xlen_t>(j)];
      if (!(sd > 0.0)) continue;
      const double delta = (observations[static_cast<R_xlen_t>(j)] -
                            means[static_cast<R_xlen_t>(j)] -
                            loadings[static_cast<R_xlen_t>(j)] * z) / sd;
      const double log_f = -log_sqrt_2pi - std::log(sd) - 0.5 * delta * delta;
      out[j] = adaptive_log_add_exp(out[j], log_w + log_f);
    }
  }
  return Rcpp::NumericVector(out.begin(), out.end());
}

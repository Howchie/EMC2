// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// The weights reproduce R's arithmetic operation for operation.
#pragma GCC optimize ("no-fast-math")

arma::vec fast_dmvnorm_rooti(const arma::mat& x, const arma::rowvec& mean,
                             const arma::mat& rooti, double log_const);

// Assemble the conditional-IS weights for one non-marginalized parameter
// component. Proposal generation, likelihood evaluation, random selection,
// diagnostics and adaptation remain in R so the transition schedule and RNG
// call order are unchanged.
// [[Rcpp::export]]
Rcpp::NumericVector emc2_particle_log_weights(
    const arma::mat& proposals,
    Rcpp::NumericVector log_likelihood,
    Rcpp::NumericVector group_log_density,
    Rcpp::NumericVector mix,
    Rcpp::IntegerVector components,
    Rcpp::List means,
    Rcpp::List rootis,
    Rcpp::NumericVector log_constants,
    Rcpp::IntegerVector auxiliary_components,
    double reference_log_likelihood) {
  const R_xlen_t n = proposals.n_rows;
  if (n == 0 || log_likelihood.size() != n || group_log_density.size() != n)
    Rcpp::stop("particle weight inputs have inconsistent row counts");

  std::vector<arma::vec> densities(rootis.size());
  for (int j = 0; j < components.size(); ++j) {
    const int k = components[j] - 1;
    if (k < 0 || k >= rootis.size() || k >= means.size() ||
        k >= log_constants.size())
      Rcpp::stop("particle weight component index is out of range");
    if (k == 0) {
      densities[k] = Rcpp::as<arma::vec>(group_log_density);
    } else if (densities[k].n_elem == 0) {
      if (Rf_isNull(rootis[k]))
        Rcpp::stop("particle weight component is missing its Gaussian factor");
      densities[k] = fast_dmvnorm_rooti(
        proposals, Rcpp::as<arma::rowvec>(means[k]),
        Rcpp::as<arma::mat>(rootis[k]), log_constants[k]);
    }
  }

  std::vector<arma::vec> log_mix_terms;
  log_mix_terms.reserve(components.size());
  for (int j = 0; j < components.size(); ++j) {
    const int k = components[j] - 1;
    const double log_weight = std::log(mix[k]);
    arma::vec term = densities[k];
    term += log_weight;
    log_mix_terms.push_back(std::move(term));
  }

  arma::vec log_mixture(n);
  for (arma::uword i = 0; i < n; ++i) {
    double largest = log_mix_terms[0][i];
    for (R_xlen_t j = 1; j < static_cast<R_xlen_t>(log_mix_terms.size()); ++j) {
      const double value = log_mix_terms[j][i];
      if (std::isnan(largest) || std::isnan(value)) {
        largest = std::numeric_limits<double>::quiet_NaN();
        break;
      }
      if (value > largest) largest = value;
    }
    // rowSums() accumulates in extended precision.
    long double sum = 0.0L;
    for (const auto& term : log_mix_terms) sum += std::exp(term[i] - largest);
    log_mixture[i] = largest + std::log(static_cast<double>(sum));
  }

  double smallest_finite = std::numeric_limits<double>::infinity();
  for (double value : log_mixture)
    if (std::isfinite(value)) smallest_finite = std::min(smallest_finite, value);
  const double replacement = std::isfinite(smallest_finite) ? smallest_finite : -1e10;
  for (double& value : log_mixture)
    if (!std::isfinite(value)) value = replacement;

  arma::vec auxiliary_density(n, arma::fill::zeros);
  for (int j = 0; j < auxiliary_components.size(); ++j) {
    const int k = auxiliary_components[j] - 1;
    if (k < 0 || k >= static_cast<int>(densities.size()) ||
        densities[k].n_elem != n)
      Rcpp::stop("particle auxiliary component is absent from the mixture");
    auxiliary_density += densities[k];
  }

  arma::vec log_weights(n);
  const double reference_ll = log_likelihood[0];
  for (arma::uword i = 0; i < n; ++i) {
    const double ll_total = (log_likelihood[i] + reference_log_likelihood) - reference_ll;
    log_weights[i] = ((ll_total + group_log_density[i]) + auxiliary_density[i]) -
      log_mixture[i];
  }

  bool has_nan = false;
  double largest_weight = -std::numeric_limits<double>::infinity();
  for (double value : log_weights) {
    if (std::isnan(value)) has_nan = true;
    else if (value > largest_weight) largest_weight = value;
  }
  Rcpp::NumericVector weights(n);
  if (has_nan) {
    std::fill(weights.begin(), weights.end(), NA_REAL);
  } else {
    for (arma::uword i = 0; i < n; ++i)
      weights[i] = std::exp(log_weights[i] - largest_weight);
  }
  return weights;
}

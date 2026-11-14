// The group draws of the factor hierarchies (src/gibbs_factor.cpp).
#ifndef EMC2_GIBBS_FACTOR_H
#define EMC2_GIBBS_FACTOR_H

#include <cstddef>
#include <vector>

namespace emc_factor {

// A column-major matrix.
struct Mat {
  int r, c;
  std::vector<double> v;
  Mat(int rows = 0, int cols = 0)
      : r(rows), c(cols), v(static_cast<std::size_t>(rows) * cols, 0.0) {}
  double& operator()(int i, int j) { return v[i + static_cast<std::size_t>(j) * r]; }
  double operator()(int i, int j) const { return v[i + static_cast<std::size_t>(j) * r]; }
};

// The hierarchy's state: scores (subjects x factors), working loadings
// (parameters x factors), means, error precisions and factor precisions; the
// draw adds the implied covariance and the loadings on unit-variance factors.
struct State {
  Mat eta, lambda, var, lambda_scaled;
  std::vector<double> mu, sig_err_inv, psi_inv;
};

// The prior: the means' precision and mean, the gamma shapes (with the
// subjects' half count added) and rates of the error and factor precisions,
// and each parameter's loading precision.
struct Prior {
  Mat mu_invar;
  std::vector<double> mu_mean, shape_err, rate_err, lambda_invar, shape_psi, rate_psi;
};

// The draw given the subjects' states (subjects x parameters) and which
// loadings are free (parameters x factors, column-major). FALSE when it is
// left to the R form; the state and the random stream are then spent.
bool gibbs(State& state, const Mat& alpha, const Prior& prior,
           const std::vector<int>& free_entry);

// The shrinkage-prior hierarchy's state: scores, loadings, the loadings'
// local shrinkage (parameters x factors), means, error precisions and the
// factors' global shrinkage; the draw adds the implied covariance.
struct InfiniteState {
  Mat eta, lambda, psi, var;
  std::vector<double> mu, epsilon_inv, delta;
};

// Its prior: the means' precisions and means, the local shrinkage's degrees
// of freedom, the global shrinkage's shapes and rates (the second pair only
// when `has_ad2`), and the error precisions' shapes (with the subjects' half
// count added) and rates.
struct InfinitePrior {
  std::vector<double> mu_invar, mu_mean, shape_err, rate_err;
  double df, ad1, bd1, ad2, bd2;
  bool has_ad2;
};

// The draw given the subjects' states; FALSE as for gibbs().
bool gibbs_infinite(InfiniteState& state, const Mat& alpha, const InfinitePrior& prior);

// The structural equation model's state: scores, loadings, the factors'
// structural coefficients B, the covariates' effects on parameters (K) and
// factors (G), the factors' disturbance precision, means and error
// precisions.
struct SemState {
  Mat eta, lambda, B, K, G, delta_inv;
  std::vector<double> mu, epsilon_inv;
};

// Its prior: the means' precisions and means, the error precisions' shapes
// (with the subjects' half count added) and rates, the prior precisions of a
// parameter's coefficients (covariates then factors) and of a factor's
// (covariates then factors), and which coefficients are free (column-major).
struct SemPrior {
  std::vector<double> mu_prec, mu_mean, shape_err, rate_err, load_prec, struct_prec;
  std::vector<int> free_lambda, free_K, free_B, free_G;
};

// The draw of everything but the disturbance precision, given the subjects'
// states and covariates (subjects x covariates); FALSE as for gibbs().
bool gibbs_sem(SemState& state, const Mat& y, const Mat& cov, const SemPrior& prior);

// What the disturbance precision's draw needs: each factor's block (numbered
// from zero in order of first appearance), the prior's shape and rate, and
// the covariates' means and covariance.
struct SemBlocks {
  std::vector<int> block;
  double a_d, b_d;
  std::vector<double> x_mean;
  Mat x_var;
};

// The subjects' moments with the scores integrated out: each subject's mean
// (parameters x subjects) and their common covariance, and both averaged
// over the covariates.
struct SemMoments {
  Mat subj_mu, subj_var, pop_var;
  std::vector<double> pop_mean;
};

// The disturbance precision given the rest of the draw, block by block, and
// the moments the completed draw implies; FALSE as for gibbs().
bool finish_sem(SemState& state, const Mat& cov, const SemBlocks& blocks, SemMoments& moments);

}  // namespace emc_factor

#endif

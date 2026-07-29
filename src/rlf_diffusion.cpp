// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Race Lévy Flight (RLF) model entry points: Rcpp interface.
//
// Solves the Space-Fractional Fokker-Planck Equation (SFFPE) for the RLF model
// using the fixed-boundary pre-factored dense LU Crank-Nicolson solver (Option A).
// Also provides the Chambers-Mallows-Stuck (CMS) path simulator.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <Rcpp.h>
#include "model_RLF.h"

using namespace Rcpp;

namespace {

Rcpp::List rlf_package(const rlf::RLF_Result& r, const NumericVector& t) {
  const int n = t.size();
  NumericVector pdf(n), cdf(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti) || ti <= 0.0) { pdf[i] = 0.0; cdf[i] = 0.0; continue; }
    pdf[i] = rlf::grid_lookup(r.t, r.pdf, ti);
    cdf[i] = rlf::grid_lookup(r.t, r.cdf, ti);
  }
  return Rcpp::List::create(
    _["pdf"] = pdf,
    _["cdf"] = cdf,
    _["t_grid"] = Rcpp::wrap(r.t),
    _["pdf_grid"] = Rcpp::wrap(r.pdf),
    _["cdf_grid"] = Rcpp::wrap(r.cdf));
}

double rlf_t_max(const NumericVector& t) {
  double m = 0.0;
  for (int i = 0; i < t.size(); ++i) if (R_finite(t[i]) && t[i] > m) m = t[i];
  return (m > 0.0) ? m : 1.0;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// RLF PDE First-Hitting-Time Solver
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rlf_fht_pdf_cdf_vec(NumericVector t, double v, double sigma,
                               double alpha, double b0, double z0 = 0.0,
                               int nx = 200, int nt = 400) {
  if (sigma <= 0.0) stop("rlf_fht_pdf_cdf_vec: sigma must be positive.");
  if (alpha <= 1.0 || alpha > 2.0) stop("rlf_fht_pdf_cdf_vec: alpha must be in (1, 2].");
  if (b0 <= 0.0) stop("rlf_fht_pdf_cdf_vec: boundary b0 must be positive.");

  const double t_max = rlf_t_max(t);

  rlf::RLF_Model m;
  m.v = v;
  m.sigma = sigma;
  m.alpha = alpha;
  m.b0 = b0;
  m.z0 = std::max(0.0, z0);

  return rlf_package(rlf::rlf_solve(m, t_max, nx, nt), t);
}

// ---------------------------------------------------------------------------
// RLF Stochastic Path Simulator
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector simulate_rlf_hit_times_cpp(int n_sims, double v, double sigma,
                                         double alpha, double b0, double z0 = 0.0,
                                         double t_max = 5.0, double dt = 0.001,
                                         unsigned int seed = 42) {
  if (n_sims <= 0) stop("simulate_rlf_hit_times_cpp: n_sims must be positive.");
  if (sigma <= 0.0) stop("simulate_rlf_hit_times_cpp: sigma must be positive.");
  if (alpha <= 0.0 || alpha > 2.0) stop("simulate_rlf_hit_times_cpp: alpha must be in (0, 2].");
  if (b0 <= 0.0) stop("simulate_rlf_hit_times_cpp: b0 must be positive.");

  const auto res = rlf::simulate_rlf_hit_times(n_sims, v, sigma, alpha, b0, z0, t_max, dt, seed);
  return Rcpp::wrap(res);
}

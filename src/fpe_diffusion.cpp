// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Fokker-Planck first-passage solvers for leaky accumulators: R entry points.
//
// This is the only translation unit that includes fpe_solver.h / fpe_models.h.
// It deliberately does NOT include the Volterra headers
// (utils_reducible_diffusion.h, model_OU_Volterra.h, model_BM_Volterra.h) --
// those define non-inline free functions owned by volterra_diffusion.cpp, and
// cross-validation between the two solvers happens in R, not in C++.
//
// None of these entry points are referenced by the model dispatch / design
// machinery: like the Volterra path, they are reachable from R for development
// and validation only.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <Rcpp.h>
#include "fpe_models.h"

using namespace Rcpp;

namespace {

// Interpolate a solve onto the requested RTs and package for R.
Rcpp::List fpe_package(const fpe::FPE_Result& r, const NumericVector& t) {
  const int n = t.size();
  NumericVector pdf(n), cdf(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti) || ti <= 0.0) { pdf[i] = 0.0; cdf[i] = 0.0; continue; }
    pdf[i] = fpe::grid_lookup(r.t, r.pdf, ti);
    cdf[i] = fpe::grid_lookup(r.t, r.cdf, ti);
  }
  return Rcpp::List::create(
    _["pdf"] = pdf,
    _["cdf"] = cdf,
    _["mismatch"] = r.flux_mass_mismatch,
    _["t_grid"] = Rcpp::wrap(r.t),
    _["pdf_grid"] = Rcpp::wrap(r.pdf));
}

double fpe_t_max(const NumericVector& t) {
  double m = 0.0;
  for (int i = 0; i < t.size(); ++i) if (R_finite(t[i]) && t[i] > m) m = t[i];
  return (m > 0.0) ? m : 1.0;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Brownian motion with drift.  Start ~ Uniform(0, z0); z0 == 0 is a point start
// at 0.  Boundary b(t) = binf + (|b0|-|binf|) * exp(-(t/tau)^pow).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List fpe_bm_fht_pdf_cdf_vec(NumericVector t, double mu, double sigma,
                                  double z0, double b0, double binf,
                                  double tau = 1.0, double pow = 1.0,
                                  int nx = 256, int nt = 512) {
  if (sigma <= 0.0) stop("fpe_bm_fht_pdf_cdf_vec: sigma must be positive.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelBM m;
  m.A = mu;
  m.sigma = sigma;
  m.bnd.set(b0, binf, tau, pow, false);
  const double z_lo = 0.0, z_hi = std::max(0.0, z0);
  m.xlo = fpe::fpe_x_lo_bm(std::min(0.0, z_lo), mu, sigma, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_bm_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt), t);
}

// ---------------------------------------------------------------------------
// Ornstein-Uhlenbeck.  dX = -lambda (X - theta) dt + sigma dW.
// Start ~ Uniform(0, z0).  Crossings require theta to pull the process across
// b0, i.e. the usual leaky-accumulator setup has theta > b0 > 0.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List fpe_ou_fht_pdf_cdf_vec(NumericVector t, double lambda, double theta,
                                  double sigma, double z0, double b0, double binf,
                                  double tau = 1.0, double pow = 1.0,
                                  int nx = 256, int nt = 512) {
  if (sigma <= 0.0) stop("fpe_ou_fht_pdf_cdf_vec: sigma must be positive.");
  if (lambda < 0.0) stop("fpe_ou_fht_pdf_cdf_vec: lambda must be non-negative.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelOU m;
  m.lambda = lambda;
  m.theta = theta;
  m.sigma = sigma;
  m.bnd.set(b0, binf, tau, pow, false);
  const double z_lo = 0.0, z_hi = std::max(0.0, z0);
  m.xlo = fpe::fpe_x_lo_ou(std::min(0.0, z_lo), lambda, theta, sigma, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_ou_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt), t);
}

// ---------------------------------------------------------------------------
// Geometric Brownian motion.  dX = mu X dt + sigma X dW, solved as Y = log X,
// which is BM with drift mu - sigma^2/2 against the barrier log b(t).
// Start ~ Uniform(start_floor, z0 + start_floor), matching
// simulate_gbm_hit_times_bb (which does `z0 += start_floor`).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List fpe_gbm_fht_pdf_cdf_vec(NumericVector t, double mu, double sigma,
                                   double z0, double b0, double binf,
                                   double tau = 1.0, double pow = 1.0,
                                   int nx = 256, int nt = 512,
                                   double start_floor = 1.0) {
  if (sigma <= 0.0) stop("fpe_gbm_fht_pdf_cdf_vec: sigma must be positive.");
  if (start_floor <= 0.0) stop("fpe_gbm_fht_pdf_cdf_vec: start_floor must be positive.");
  if (b0 <= 0.0 || binf <= 0.0) stop("fpe_gbm_fht_pdf_cdf_vec: boundaries must stay positive.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelBM m;
  m.A = mu - 0.5 * sigma * sigma;      // Ito correction for Y = log X
  m.sigma = sigma;
  m.bnd.set(b0, binf, tau, pow, true); // barrier is log b(t)

  const double z_lo = std::log(start_floor);
  const double z_hi = std::log(z0 + start_floor);
  m.xlo = fpe::fpe_x_lo_bm(z_lo, m.A, sigma, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_gbm_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt), t);
}

// ---------------------------------------------------------------------------
// Gompertz growth.  Y = log X is OU with lambda = alpha, sigma = beta and
// theta = log(kinf) - beta^2/(2 alpha) -- the same reduction the reference
// simulator uses (simulate_gompertz_hit_times_bb), where the carrying capacity
// kinf sets the reversion level and k0 -> kinf is the collapsing barrier.
// Start ~ Uniform(start_floor, z0).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List fpe_gompertz_fht_pdf_cdf_vec(NumericVector t, double alpha, double beta,
                                        double z0, double k0, double kinf,
                                        double tau = 1.0, double pow = 1.0,
                                        int nx = 256, int nt = 512,
                                        double start_floor = 1e-3) {
  if (beta <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: beta must be positive.");
  if (alpha <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: alpha must be positive.");
  if (k0 <= 0.0 || kinf <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: capacities must stay positive.");
  if (start_floor <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: start_floor must be positive.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelOU m;
  m.lambda = alpha;
  m.sigma = beta;
  m.theta = std::log(kinf) - (beta * beta) / (2.0 * alpha);
  m.bnd.set(k0, kinf, tau, pow, true);

  const double z_lo = std::log(start_floor);
  const double z_hi = std::log(std::max(z0, start_floor));
  m.xlo = fpe::fpe_x_lo_ou(z_lo, alpha, m.theta, beta, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_gompertz_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt), t);
}

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
  NumericVector pdf(n), cdf(n), surv(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti) || ti <= 0.0) { pdf[i] = 0.0; cdf[i] = 0.0; surv[i] = 1.0; continue; }
    // PDF is zero before the seed time; CDF is not (see fpe::grid_lookup).
    pdf[i] = fpe::grid_lookup(r.t, r.pdf, ti, true);
    cdf[i] = fpe::grid_lookup(r.t, r.cdf, ti);
    surv[i] = fpe::grid_lookup(r.t, r.surv, ti);
  }
  return Rcpp::List::create(
    _["pdf"] = pdf,
    _["cdf"] = cdf,
    _["surv"] = surv,
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
                                  int nx = 256, int nt = 512,
                                  double grade = 8.0,     // fpe::FPE_GRADE
                                  double tgrade = 1.0) {
  if (sigma <= 0.0) stop("fpe_bm_fht_pdf_cdf_vec: sigma must be positive.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelBM m;
  m.A = mu;
  m.sigma = sigma;
  m.bnd.set(b0, binf, tau, pow, false);
  const double z_lo = 0.0, z_hi = std::max(0.0, z0);
  m.xlo = fpe::fpe_x_lo_bm(std::min(0.0, z_lo), mu, sigma, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_bm_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt, grade, tgrade), t);
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
                                  int nx = 256, int nt = 512,
                                  double grade = 8.0,     // fpe::FPE_GRADE
                                  double tgrade = 1.0) {
  if (sigma <= 0.0) stop("fpe_ou_fht_pdf_cdf_vec: sigma must be positive.");
  if (lambda < 0.0) stop("fpe_ou_fht_pdf_cdf_vec: lambda must be non-negative.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelOU m;
  m.set_lambda_theta(lambda, theta);
  m.sigma = sigma;
  m.bnd.set(b0, binf, tau, pow, false);
  const double z_lo = 0.0, z_hi = std::max(0.0, z0);
  m.xlo = fpe::fpe_x_lo_ou(std::min(0.0, z_lo), m.v, lambda, sigma, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_ou_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt, grade, tgrade), t);
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
                                   double start_floor = 1.0,
                                   double grade = 8.0, double tgrade = 1.0) {
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
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt, grade, tgrade), t);
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
                                        double start_floor = 1e-3,
                                        double grade = 8.0, double tgrade = 1.0) {
  if (beta <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: beta must be positive.");
  if (alpha <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: alpha must be positive.");
  if (k0 <= 0.0 || kinf <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: capacities must stay positive.");
  if (start_floor <= 0.0) stop("fpe_gompertz_fht_pdf_cdf_vec: start_floor must be positive.");
  const double t_max = fpe_t_max(t);

  fpe::FPE_ModelOU m;
  m.sigma = beta;
  m.set_lambda_theta(alpha, std::log(kinf) - (beta * beta) / (2.0 * alpha));
  m.bnd.set(k0, kinf, tau, pow, true);

  const double z_lo = std::log(start_floor);
  const double z_hi = std::log(std::max(z0, start_floor));
  m.xlo = fpe::fpe_x_lo_ou(z_lo, m.v, alpha, beta, t_max);

  if (m.bnd.a(t_max) <= m.xlo) stop("fpe_gompertz_fht_pdf_cdf_vec: boundary collapses below the domain.");
  return fpe_package(fpe::fpe_run(m, z_lo, z_hi, t_max, nx, nt, grade, tgrade), t);
}

// ---------------------------------------------------------------------------
// Bounded OU (Smith & Ratcliff 2004): absorbing at BOTH 0 and a, start at
// z = Z*a, decay `beta` toward `anchor` (default the start point).  beta = 0
// reduces to the Wiener diffusion, which is the oracle this is validated
// against.
//
// Returns the DEFECTIVE per-response density and cdf under the same convention
// as dDDM/pDDM (R/model_DDM.R): `upper`/`lower` each integrate to that
// response's probability, and the two sum to 1.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List fpe_bou_fht_pdf_cdf_vec(NumericVector t, double v, double beta,
                                   double a, double Z, double sigma,
                                   double SZ = 0.0,
                                   double anchor = NA_REAL,
                                   int bkind = 0, double aInf = 0.0,
                                   double tau = 0.0, double pw = 0.0,
                                   int nx = 256, int nt = 512,
                                   double grade = 8.0,
                                   double tgrade = 1.0) {
  if (sigma <= 0.0) stop("fpe_bou_fht_pdf_cdf_vec: sigma must be positive.");
  if (a <= 0.0) stop("fpe_bou_fht_pdf_cdf_vec: a must be positive.");
  if (Z <= 0.0 || Z >= 1.0) stop("fpe_bou_fht_pdf_cdf_vec: Z must be in (0,1).");
  const double t_max = fpe_t_max(t);

  const double z = Z * a;
  // Same mapping as the DDM's Ttransform (R/model_DDM.R:155): SZ is a
  // proportion, widened to the largest symmetric range that stays inside [0,a].
  const double sz = 2.0 * SZ * std::min(z, a - z);

  fpe::FPE_ModelBoundedOU m;
  m.v = v;
  m.beta = beta;
  m.anchor = R_finite(anchor) ? anchor : z;
  m.sigma = sigma;
  m.set_separation(a, bkind, aInf, tau, pw);

  const double z_lo = z - 0.5 * sz;
  const double z_hi = z + 0.5 * sz;

  fpe::FPE_Mesh g;
  g.build(nx, grade, fpe::FPE_ModelBoundedOU::symmetric_mesh);
  std::vector<double> q0;
  double absorbed_lower = 0.0;
  const double t0 = fpe::fpe_seed(m, z_lo, z_hi, g, t_max, q0, &absorbed_lower);
  const fpe::FPE_Result r =
    fpe::fpe_solve(m, q0, t0, t_max, g, std::max(1, nt), tgrade, absorbed_lower);

  const int n = t.size();
  NumericVector pdf_up(n), cdf_up(n), pdf_lo(n), cdf_lo(n), surv(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti) || ti <= 0.0) { surv[i] = 1.0; continue; }
    pdf_up[i] = fpe::grid_lookup(r.t, r.pdf, ti, true);
    cdf_up[i] = fpe::grid_lookup(r.t, r.cdf, ti);
    pdf_lo[i] = fpe::grid_lookup(r.t, r.pdf_lower, ti, true);
    cdf_lo[i] = fpe::grid_lookup(r.t, r.cdf_lower, ti);
    surv[i]   = fpe::grid_lookup(r.t, r.surv, ti);
  }
  return Rcpp::List::create(
    _["pdf_upper"] = pdf_up, _["cdf_upper"] = cdf_up,
    _["pdf_lower"] = pdf_lo, _["cdf_lower"] = cdf_lo,
    _["surv"] = surv,
    _["mismatch"] = r.flux_mass_mismatch,
    _["t_grid"] = Rcpp::wrap(r.t),
    _["pdf_upper_grid"] = Rcpp::wrap(r.pdf),
    _["pdf_lower_grid"] = Rcpp::wrap(r.pdf_lower));
}

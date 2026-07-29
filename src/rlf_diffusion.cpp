// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Race Lévy Flight (RLF) model entry points: Rcpp interface.
//
// Solves the nonlocal Fokker-Planck equation for the RLF model using a
// fixed-boundary, pre-factorised Rannacher/Crank-Nicolson backend.  Also
// provides the Chambers-Mallows-Stuck (CMS) path simulator.
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
    _["mismatch"] = r.flux_mass_mismatch,
    _["lower_boundary_pressure"] = r.lower_boundary_pressure,
    _["operator_conservation_error"] = r.operator_conservation_error,
    _["min_density"] = r.min_density,
    _["domain_pdf_error"] = r.domain_pdf_error,
    _["domain_cdf_error"] = r.domain_cdf_error,
    _["spatial_pdf_error"] = r.spatial_pdf_error,
    _["spatial_cdf_error"] = r.spatial_cdf_error,
    _["x_lo"] = r.x_lo,
    _["dx"] = r.dx,
    _["nx_used"] = r.nx_used,
    _["nt_used"] = r.nt_used,
    _["domain_refinements"] = r.domain_refinements,
    _["spatial_refinements"] = r.spatial_refinements,
    _["refinement_checked"] = r.refinement_checked,
    _["refinement_skipped"] = r.refinement_skipped,
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
  if (!R_FINITE(v) || v <= 0.0)
    stop("rlf_fht_pdf_cdf_vec: v must be finite and positive.");
  if (!R_FINITE(sigma) || sigma <= 0.0)
    stop("rlf_fht_pdf_cdf_vec: sigma must be finite and positive.");
  if (!R_FINITE(alpha) || alpha <= 1.0 || alpha > 2.0)
    stop("rlf_fht_pdf_cdf_vec: alpha must be finite and in (1, 2].");
  if (!R_FINITE(b0) || b0 <= 0.0)
    stop("rlf_fht_pdf_cdf_vec: boundary b0 must be finite and positive.");
  if (!R_FINITE(z0) || z0 < 0.0 || z0 >= b0)
    stop("rlf_fht_pdf_cdf_vec: z0 must be finite and in [0, b0).");
  if (nx < 30) stop("rlf_fht_pdf_cdf_vec: nx must be at least 30.");
  if (nt < 50) stop("rlf_fht_pdf_cdf_vec: nt must be at least 50.");

  const double t_max = rlf_t_max(t);

  rlf::RLF_Model m;
  m.v = v;
  m.sigma = sigma;
  m.alpha = alpha;
  m.b0 = b0;
  m.z0 = z0;

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
  if (!R_FINITE(v) || v <= 0.0)
    stop("simulate_rlf_hit_times_cpp: v must be finite and positive.");
  if (!R_FINITE(sigma) || sigma <= 0.0)
    stop("simulate_rlf_hit_times_cpp: sigma must be finite and positive.");
  if (!R_FINITE(alpha) || alpha <= 1.0 || alpha > 2.0)
    stop("simulate_rlf_hit_times_cpp: alpha must be finite and in (1, 2].");
  if (!R_FINITE(b0) || b0 <= 0.0)
    stop("simulate_rlf_hit_times_cpp: b0 must be finite and positive.");
  if (!R_FINITE(z0) || z0 < 0.0 || z0 >= b0)
    stop("simulate_rlf_hit_times_cpp: z0 must be finite and in [0, b0).");
  if (!R_FINITE(t_max) || t_max <= 0.0)
    stop("simulate_rlf_hit_times_cpp: t_max must be finite and positive.");
  if (!R_FINITE(dt) || dt <= 0.0)
    stop("simulate_rlf_hit_times_cpp: dt must be finite and positive.");

  const auto res = rlf::simulate_rlf_hit_times(n_sims, v, sigma, alpha, b0, z0, t_max, dt, seed);
  return Rcpp::wrap(res);
}

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
// Grouped single-accumulator PDF/CDF queries used by dRLF/pRLF.  Parameter rows
// sharing a sigma-scaled key pay for one solve, and t0 is only a time shift.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rlf_pdf_cdf_vec(
    NumericVector rt, NumericVector v, NumericVector B, NumericVector A,
    NumericVector t0, NumericVector s, NumericVector alpha, int nx = 200,
    double dt_target = 5e-3, double tgrade = 1.0, bool adaptive = false,
    bool explicit_inverse = true, bool sparse_output = true,
    bool simd_batch = true) {
  const int n = rt.size();
  if (v.size() != n || B.size() != n || A.size() != n ||
      t0.size() != n || s.size() != n || alpha.size() != n) {
    stop("rlf_pdf_cdf_vec: all parameter vectors must match length(rt).");
  }
  if (nx < 30 || !(dt_target > 0.0) || !(tgrade >= 1.0)) {
    stop("rlf_pdf_cdf_vec: invalid grid configuration.");
  }

  NumericVector pdf(n, 0.0), cdf(n, 0.0);
  rlf::SolveCache cache;
  cache.grid.nx = nx;
  cache.grid.dt_target = dt_target;
  cache.grid.tgrade = tgrade;
  cache.grid.adaptive = adaptive;
  cache.grid.explicit_inverse = explicit_inverse;
  cache.grid.sparse_output = sparse_output;
  cache.grid.simd_batch = simd_batch;

  std::vector<rlf::Key> keys;
  std::vector<double> horizons;
  std::vector<std::vector<double>> query_times;
  std::vector<int> row_group(n, -1);
  std::unordered_map<rlf::Key, int, rlf::KeyHash> groups;
  for (int i = 0; i < n; ++i) {
    const double tt = rt[i] - t0[i];
    if (!R_finite(tt) || !(tt > 0.0)) continue;
    rlf::Key key;
    if (!rlf::rlf_key(v[i], s[i], alpha[i], B[i], A[i], key)) continue;
    const auto found = groups.find(key);
    int group = -1;
    if (found == groups.end()) {
      group = static_cast<int>(keys.size());
      groups.emplace(key, group);
      keys.push_back(key);
      horizons.push_back(tt);
      query_times.push_back(std::vector<double>(1, tt));
    } else {
      group = found->second;
      horizons[group] = std::max(horizons[group], tt);
      query_times[group].push_back(tt);
    }
    row_group[i] = group;
  }
  for (auto& times : query_times) {
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }

  std::vector<int> cache_index;
  rlf::cache_get_batch(
    cache, keys, horizons, cache_index,
    sparse_output ? &query_times : nullptr);
  for (int i = 0; i < n; ++i) {
    if (row_group[i] < 0) continue;
    const rlf::Entry& entry = cache.entries[cache_index[row_group[i]]];
    const double tt = rt[i] - t0[i];
    const double log_pdf = rlf::entry_log_pdf(entry, tt);
    pdf[i] = log_pdf <= rlf::RLF_LOG_FLOOR ? 0.0 : std::exp(log_pdf);
    const double log_survivor = rlf::entry_log_S(entry, tt);
    cdf[i] = log_survivor >= 0.0 ? 0.0 :
      (log_survivor <= rlf::RLF_LOG_FLOOR
        ? 1.0 : -std::expm1(log_survivor));
  }

  return Rcpp::List::create(
    _["pdf"] = pdf, _["cdf"] = cdf,
    _["n_solves"] = static_cast<int>(cache.solve_count),
    _["n_keys"] = static_cast<int>(keys.size()));
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

// [[Rcpp::export]]
NumericVector rlf_hit_times_vec(
    NumericVector v, NumericVector B, NumericVector A, NumericVector s,
    NumericVector alpha, double dt = 1e-3, double t_max = 30.0) {
  const int n = v.size();
  if (B.size() != n || A.size() != n || s.size() != n ||
      alpha.size() != n) {
    stop("rlf_hit_times_vec: all parameter vectors must have equal length.");
  }
  if (!R_FINITE(dt) || !R_FINITE(t_max) ||
      !(dt > 0.0) || !(t_max > 0.0)) {
    stop("rlf_hit_times_vec: dt and t_max must be finite and positive.");
  }
  const double required_steps = std::ceil(t_max / dt);
  if (!(required_steps <=
        static_cast<double>(std::numeric_limits<int>::max()))) {
    stop("rlf_hit_times_vec: t_max/dt requires too many steps.");
  }

  NumericVector out(n, R_PosInf);
  Rcpp::RNGScope scope;
  const int n_steps = static_cast<int>(required_steps);
  for (int row = 0; row < n; ++row) {
    rlf::Key key;
    if (!rlf::rlf_key(v[row], s[row], alpha[row], B[row], A[row], key)) {
      continue;
    }
    double x = A[row] > 0.0 ? A[row] * ::unif_rand() : 0.0;
    double time = 0.0;
    for (int step = 0; step < n_steps; ++step) {
      const double step_dt =
        step + 1 == n_steps ? t_max - step * dt : dt;
      const double scale =
        s[row] * std::pow(0.5 * step_dt, 1.0 / alpha[row]);
      double stable = 0.0;
      if (alpha[row] == 2.0) {
        stable = std::sqrt(2.0) * ::norm_rand();
      } else {
        const double u01 =
          std::min(1.0 - 1e-15, std::max(1e-15, ::unif_rand()));
        const double U = (u01 - 0.5) * M_PI;
        const double w01 =
          std::min(1.0 - 1e-15, std::max(1e-15, ::unif_rand()));
        const double W = -std::log(w01);
        const double aU = alpha[row] * U;
        const double num = std::sin(aU);
        const double den = std::pow(std::cos(U), 1.0 / alpha[row]);
        const double tail = std::pow(
          std::cos((1.0 - alpha[row]) * U) / W,
          (1.0 - alpha[row]) / alpha[row]);
        stable = (num / den) * tail;
      }
      x += v[row] * step_dt + scale * stable;
      time += step_dt;
      if (x >= B[row] + A[row]) {
        out[row] = time;
        break;
      }
    }
  }
  return out;
}

#ifndef model_RLF_h
#define model_RLF_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Race Lévy Flight (RLF) model and nonlocal Fokker-Planck solver.
//
// Process convention
// ------------------
//   dX_t = v dt + sigma dL_t^(alpha)
//
// where the symmetric stable increment over dt has characteristic function
//
//   E exp(i k sigma dL) = exp(-0.5 * sigma^alpha * dt * |k|^alpha).
//
// This convention is deliberately continuous at alpha = 2: dL^(2) is an
// ordinary Wiener increment and sigma is the Brownian standard deviation used
// by the rest of EMC2.
//
// The nonlocal generator cannot use the tridiagonal local-diffusion operator in
// fpe_solver.h.  This backend nevertheless follows the same architecture:
//   * no Rcpp/model-dispatch knowledge in the numerical core;
//   * a time-invariant operator factorised once for a fixed boundary;
//   * Rannacher backward-Euler startup before Crank-Nicolson;
//   * PDF from upper-boundary probability flux and CDF from surviving mass;
//   * an independent flux/mass mismatch convergence diagnostic.
//
// Spatial truncation
// ------------------
// The physical domain is (-infinity, b0).  A finite lower edge is unavoidable
// numerically, but killing paths there would incorrectly count lower exits as
// upper first passages.  Attempts to jump below the lower edge are therefore
// censored (a conservative/no-flux closure).  The solver automatically expands
// this lower domain at fixed spatial resolution, then refines the spatial grid,
// until independent PDF/CDF comparisons stabilize.  Refinement is bounded and
// a non-converged solve is rejected rather than silently returned.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rlf {

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------------------------
// Dense LU factorisation with partial pivoting.
// ---------------------------------------------------------------------------
struct RLF_DenseLU {
  int n = 0;
  std::vector<double> LU;
  std::vector<int> pvt;

  bool factor(int n_, const std::vector<double>& A) {
    n = n_;
    LU = A;
    pvt.resize(n);
    for (int i = 0; i < n; ++i) pvt[i] = i;

    for (int k = 0; k < n; ++k) {
      int max_i = k;
      double max_val = std::abs(LU[k * n + k]);
      for (int i = k + 1; i < n; ++i) {
        const double val = std::abs(LU[i * n + k]);
        if (val > max_val) {
          max_val = val;
          max_i = i;
        }
      }
      if (!(max_val > 1e-15)) return false;

      if (max_i != k) {
        std::swap(pvt[k], pvt[max_i]);
        for (int j = 0; j < n; ++j) {
          std::swap(LU[k * n + j], LU[max_i * n + j]);
        }
      }

      const double inv_pivot = 1.0 / LU[k * n + k];
      for (int i = k + 1; i < n; ++i) {
        const double mult = LU[i * n + k] * inv_pivot;
        LU[i * n + k] = mult;
        const double* __restrict lu_k = &LU[k * n + k + 1];
        double* __restrict lu_i = &LU[i * n + k + 1];
        const int len = n - (k + 1);
        for (int j = 0; j < len; ++j) {
          lu_i[j] -= mult * lu_k[j];
        }
      }
    }
    return true;
  }

  void solve(const double* __restrict b, double* __restrict x,
             double* __restrict y) const {
    for (int i = 0; i < n; ++i) {
      double sum = b[pvt[i]];
      const double* __restrict lu_i = &LU[i * n];
      for (int j = 0; j < i; ++j) sum -= lu_i[j] * y[j];
      y[i] = sum;
    }
    for (int i = n - 1; i >= 0; --i) {
      double sum = y[i];
      const double* __restrict lu_i = &LU[i * n];
      for (int j = i + 1; j < n; ++j) sum -= lu_i[j] * x[j];
      x[i] = sum / lu_i[i];
    }
  }
};

// g_k = (-1)^k binom(alpha, k).
inline std::vector<double> compute_gl_weights(double alpha, int max_k) {
  std::vector<double> g(max_k + 1, 0.0);
  g[0] = 1.0;
  for (int k = 1; k <= max_k; ++k) {
    g[k] = (1.0 - (alpha + 1.0) / static_cast<double>(k)) * g[k - 1];
  }
  return g;
}

struct RLF_Model {
  double v = 1.0;
  double sigma = 1.0;
  double alpha = 1.7;
  double b0 = 1.0;
  double z0 = 0.0;
};

struct RLF_Result {
  std::vector<double> t, pdf, cdf;
  double flux_mass_mismatch = 0.0;
  double lower_boundary_pressure = 0.0;
  double operator_conservation_error = 0.0;
  double min_density = 0.0;
  double domain_pdf_error = std::numeric_limits<double>::quiet_NaN();
  double domain_cdf_error = std::numeric_limits<double>::quiet_NaN();
  double spatial_pdf_error = std::numeric_limits<double>::quiet_NaN();
  double spatial_cdf_error = std::numeric_limits<double>::quiet_NaN();
  double x_lo = 0.0;
  double dx = 0.0;
  int nx_used = 0;
  int nt_used = 0;
  int domain_refinements = 0;
  int spatial_refinements = 0;
  bool refinement_checked = false;
  bool refinement_skipped = false;
};

inline double grid_lookup(const std::vector<double>& tg,
                          const std::vector<double>& vg, double t) {
  const int n = static_cast<int>(tg.size());
  if (n == 0) return 0.0;
  if (t <= tg.front()) return vg.front();
  if (t >= tg.back()) return vg.back();
  const int j = static_cast<int>(
    std::lower_bound(tg.begin(), tg.end(), t) - tg.begin());
  const double w = (t - tg[j - 1]) / (tg[j] - tg[j - 1]);
  return vg[j - 1] + w * (vg[j] - vg[j - 1]);
}

struct RLF_Operator {
  int n = 0;
  std::vector<double> L;
  std::vector<double> upper_kill;
  std::vector<double> lower_censor;
  double max_exit_rate = 0.0;
  double conservation_error = 0.0;
};

// Shifted symmetric Grünwald generator on the uniform grid.  For a Fourier
// mode, the two shifted sums tend to
//
//   2 cos(pi*alpha/2) |k|^alpha.
//
// The coefficient below therefore gives the model convention
// -0.5*sigma^alpha*|k|^alpha and tends continuously to
// -0.5*sigma^2*k^2 at alpha=2.
inline RLF_Operator build_rlf_operator(const RLF_Model& m, int n, double h) {
  RLF_Operator op;
  op.n = n;
  op.L.assign(static_cast<size_t>(n) * n, 0.0);
  op.upper_kill.assign(n, 0.0);
  op.lower_censor.assign(n, 0.0);

  const auto g = compute_gl_weights(m.alpha, n + 1);
  const double cosine = std::cos(M_PI * m.alpha / 2.0);
  const double C = -(0.5 * std::pow(m.sigma, m.alpha)) /
                   (2.0 * cosine * std::pow(h, m.alpha));

  // q[d] is the transition coefficient between grid points d cells apart.
  // All q[d>0] are non-negative for 1 < alpha <= 2.
  std::vector<double> q(n, 0.0);
  q[0] = 2.0 * C * g[1];
  if (n > 1) q[1] = C * (g[0] + g[2]);
  for (int d = 2; d < n; ++d) q[d] = C * g[d + 1];

  // tail[D] is the one-sided jump rate to targets at least D cells away.
  // For D=1 it is C*alpha.  For D>=2,
  //   sum_{d=D}^infinity q[d] = C * sum_{k=D+1}^infinity g[k]
  //                            = -C * sum_{k=0}^D g[k].
  std::vector<double> tail(n + 1, 0.0);
  tail[1] = C * m.alpha;
  double partial = g[0];
  for (int k = 1; k <= n; ++k) {
    partial += g[k];
    if (k >= 2) tail[k] = std::max(0.0, -C * partial);
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      op.L[static_cast<size_t>(i) * n + j] = q[std::abs(i - j)];
    }
  }

  // Column j is a source state.  Omitted targets above the barrier are killed
  // and contribute to the FPT flux.  Omitted targets below x_lo are censored:
  // adding that rate back to the diagonal suppresses those artificial exits.
  for (int j = 0; j < n; ++j) {
    const double lower = tail[j + 1];
    const double upper = tail[n - j];
    op.lower_censor[j] = lower;
    op.upper_kill[j] = upper;
    op.L[static_cast<size_t>(j) * n + j] += lower;
  }

  // Positive-drift upwind generator.  The final transition exits through b0.
  const double adv = m.v / h;
  for (int j = 0; j < n; ++j) {
    op.L[static_cast<size_t>(j) * n + j] -= adv;
    if (j + 1 < n) {
      op.L[static_cast<size_t>(j + 1) * n + j] += adv;
    } else {
      op.upper_kill[j] += adv;
    }
  }

  for (int j = 0; j < n; ++j) {
    const double exit = -op.L[static_cast<size_t>(j) * n + j];
    op.max_exit_rate = std::max(op.max_exit_rate, exit);

    double col_sum = op.upper_kill[j];
    for (int i = 0; i < n; ++i) {
      col_sum += op.L[static_cast<size_t>(i) * n + j];
    }
    op.conservation_error =
      std::max(op.conservation_error, std::abs(col_sum));
  }
  return op;
}

inline void apply_shifted(const RLF_Operator& op, double c,
                          const std::vector<double>& p,
                          std::vector<double>& out) {
  const int n = op.n;
  out.resize(n);
  for (int i = 0; i < n; ++i) {
    double value = p[i];
    const double* __restrict row = &op.L[static_cast<size_t>(i) * n];
    for (int j = 0; j < n; ++j) value += c * row[j] * p[j];
    out[i] = value;
  }
}

inline double weighted_rate(const std::vector<double>& rate,
                            const std::vector<double>& p, double h) {
  double out = 0.0;
  const int n = static_cast<int>(p.size());
  for (int j = 0; j < n; ++j) out += rate[j] * p[j];
  return h * out;
}

inline double density_mass(const std::vector<double>& p, double h) {
  double mass = 0.0;
  for (double value : p) mass += value;
  return h * mass;
}

inline double density_min(const std::vector<double>& p) {
  return *std::min_element(p.begin(), p.end());
}

// One fixed-domain/fixed-grid solve.  `M` is the number of intervals and the
// M-1 interior nodes are x_lo+h, ..., b0-h.  The public solver below owns all
// validation and convergence refinement.
inline RLF_Result rlf_solve_fixed_grid(const RLF_Model& m, double t_max,
                                       int M, int Nt,
                                       double lower_extent) {
  const double x_lo = -lower_extent;
  const double h = (m.b0 - x_lo) / M;
  const int n = M - 1;

  const RLF_Operator op = build_rlf_operator(m, n, h);

  // Crank-Nicolson is positivity preserving for this M-matrix generator when
  // I + dt/2 L is entrywise non-negative.  Refine nt automatically so its
  // diagonal is at least 0.1; this also controls high-Peclet drift cases.
  const double safe_cn = 1.8;
  const double nt_required =
    std::ceil(t_max * op.max_exit_rate / safe_cn);
  if (!(nt_required <=
        static_cast<double>(std::numeric_limits<int>::max()))) {
    throw std::invalid_argument(
      "rlf_solve: parameters require too many stable time steps.");
  }
  const int nt_positive = static_cast<int>(nt_required);
  const int n_steps = std::max(Nt, std::max(1, nt_positive));
  const double dt = t_max / n_steps;
  const double c = 0.5 * dt;

  std::vector<double> lhs(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      lhs[static_cast<size_t>(i) * n + j] =
        (i == j ? 1.0 : 0.0) - c * op.L[static_cast<size_t>(i) * n + j];
    }
  }

  RLF_DenseLU lu;
  if (!lu.factor(n, lhs)) {
    throw std::runtime_error(
      "rlf_solve: factorization of the nonlocal operator failed.");
  }

  std::vector<double> p(n, 0.0);
  if (m.z0 > 0.0) {
    // Exact overlap with node-centred control intervals, followed by a small
    // normalization correction for the two omitted Dirichlet endpoints.
    for (int i = 0; i < n; ++i) {
      const double x = x_lo + (i + 1) * h;
      const double xa = x - 0.5 * h;
      const double xb = x + 0.5 * h;
      const double overlap =
        std::min(xb, m.z0) - std::max(xa, 0.0);
      if (overlap > 0.0) p[i] = (overlap / m.z0) / h;
    }
    const double initial_mass = density_mass(p, h);
    if (!(initial_mass > 0.0)) {
      throw std::runtime_error("rlf_solve: start distribution missed the grid.");
    }
    for (double& value : p) value /= initial_mass;
  } else {
    // Preserve the start location exactly rather than snapping it to whichever
    // node happens to be closest for this domain/grid combination.  The two
    // weights have unit total mass and first moment zero.
    const double grid_position = -x_lo / h;
    const int left_node = static_cast<int>(std::floor(grid_position));
    const double weight_right = grid_position - left_node;
    const int left_i = left_node - 1;
    if (left_i >= 0 && left_i + 1 < n) {
      p[left_i] = (1.0 - weight_right) / h;
      p[left_i + 1] = weight_right / h;
    } else {
      const int nearest =
        std::max(0, std::min(n - 1,
          static_cast<int>(std::round(grid_position)) - 1));
      p[nearest] = 1.0 / h;
    }
  }

  RLF_Result res;
  res.x_lo = x_lo;
  res.dx = h;
  res.nx_used = M;
  res.nt_used = n_steps;
  res.operator_conservation_error = op.conservation_error;
  res.min_density = density_min(p);

  const int n_rann = std::min(2, n_steps);
  const int n_levels = 1 + 2 * n_rann + (n_steps - n_rann);
  res.t.reserve(n_levels);
  res.pdf.reserve(n_levels);
  res.cdf.reserve(n_levels);

  double time = 0.0;
  double flux_prev = weighted_rate(op.upper_kill, p, h);
  double lower_prev = weighted_rate(op.lower_censor, p, h);
  double cdf_flux = 0.0;
  double lower_pressure = 0.0;

  res.t.push_back(time);
  res.pdf.push_back(flux_prev);
  res.cdf.push_back(0.0);

  std::vector<double> rhs(n, 0.0), next(n, 0.0), work(n, 0.0);

  auto record_step = [&](double step, bool backward_euler) {
    if (backward_euler) {
      rhs = p;
    } else {
      apply_shifted(op, 0.5 * step, p, rhs);
    }

    lu.solve(rhs.data(), next.data(), work.data());
    p.swap(next);
    time += step;

    const double min_p = density_min(p);
    res.min_density = std::min(res.min_density, min_p);
    if (min_p < -1e-8 / h) {
      throw std::runtime_error(
        "rlf_solve: material negative density; increase nt or nx.");
    }

    const double flux = weighted_rate(op.upper_kill, p, h);
    const double lower = weighted_rate(op.lower_censor, p, h);
    if (backward_euler) {
      cdf_flux += step * flux;
      lower_pressure += step * lower;
    } else {
      cdf_flux += 0.5 * step * (flux_prev + flux);
      lower_pressure += 0.5 * step * (lower_prev + lower);
    }
    flux_prev = flux;
    lower_prev = lower;

    const double cdf_mass = 1.0 - density_mass(p, h);
    res.flux_mass_mismatch =
      std::max(res.flux_mass_mismatch, std::abs(cdf_mass - cdf_flux));
    res.t.push_back(time);
    res.pdf.push_back(flux);
    res.cdf.push_back(cdf_mass);
  };

  // The BE half-step and CN full-step LHS are both I - dt/2 L, so one dense
  // factorization serves the entire march.
  for (int k = 0; k < n_rann; ++k) {
    record_step(0.5 * dt, true);
    record_step(0.5 * dt, true);
  }
  for (int k = n_rann; k < n_steps; ++k) {
    record_step(dt, false);
  }

  res.lower_boundary_pressure = lower_pressure;
  return res;
}

struct RLF_Comparison {
  double pdf = 0.0;
  double cdf = 0.0;
};

// Compare complete hitting-time curves independently of either solver's time
// grid.  PDF error is normalized by the larger peak (with 1/t_max as a floor),
// making the tolerance dimensionless and invariant to the time unit.
inline RLF_Comparison compare_rlf_results(const RLF_Result& coarse,
                                          const RLF_Result& fine,
                                          double t_max) {
  constexpr int RLF_COMPARE_POINTS = 128;
  double pdf_scale = 1.0 / t_max;
  for (int k = 1; k <= RLF_COMPARE_POINTS; ++k) {
    const double t =
      t_max * static_cast<double>(k) / RLF_COMPARE_POINTS;
    pdf_scale = std::max(
      pdf_scale, std::abs(grid_lookup(coarse.t, coarse.pdf, t)));
    pdf_scale = std::max(
      pdf_scale, std::abs(grid_lookup(fine.t, fine.pdf, t)));
  }

  RLF_Comparison error;
  for (int k = 1; k <= RLF_COMPARE_POINTS; ++k) {
    const double t =
      t_max * static_cast<double>(k) / RLF_COMPARE_POINTS;
    error.pdf = std::max(
      error.pdf,
      std::abs(grid_lookup(coarse.t, coarse.pdf, t) -
               grid_lookup(fine.t, fine.pdf, t)) / pdf_scale);
    error.cdf = std::max(
      error.cdf,
      std::abs(grid_lookup(coarse.t, coarse.cdf, t) -
               grid_lookup(fine.t, fine.cdf, t)));
  }
  return error;
}

inline int rlf_intervals_for_spacing(double lower_extent, double b0,
                                     double target_h) {
  const double required = std::ceil((b0 + lower_extent) / target_h);
  if (!(required <=
        static_cast<double>(std::numeric_limits<int>::max()))) {
    throw std::runtime_error(
      "rlf_solve: automatic domain expansion exceeds the grid-size limit.");
  }
  return std::max(30, static_cast<int>(required));
}

// Cheap fast-path gate.  The initial solve already measures the lower-tail
// pressure and enforces positivity/conservation.  Combine those diagnostics
// with the stable spread over the relevant first-passage time; only cases near
// a truncation or spatial-resolution boundary need expensive second solves.
inline bool rlf_fast_path_safe(const RLF_Model& m, double t_max,
                               const RLF_Result& base, int M) {
  constexpr double RLF_FAST_PRESSURE_TOL = 0.002;
  constexpr double RLF_FAST_DX_SCALE_TOL = 0.10;
  constexpr int RLF_FAST_MIN_INTERVALS = 160;
  if (M < RLF_FAST_MIN_INTERVALS) return false;
  if (!(base.lower_boundary_pressure <= RLF_FAST_PRESSURE_TOL)) return false;
  if (base.min_density < -1e-10 / base.dx) return false;
  if (base.flux_mass_mismatch > 1e-10) return false;
  if (base.operator_conservation_error > 1e-10) return false;

  const double characteristic_time =
    std::min(t_max, m.b0 / m.v);
  const double stable_spread =
    m.sigma * std::pow(0.5 * characteristic_time, 1.0 / m.alpha);
  if (!(stable_spread > 0.0)) return false;
  return base.dx / stable_spread <= RLF_FAST_DX_SCALE_TOL;
}

// Automatically converged fixed-boundary RLF solve.  M and Nt are minimum
// requested resolutions: Nt may grow for CN positivity, while M may grow first
// to verify lower-domain truncation and then to verify spatial discretization.
inline RLF_Result rlf_solve(const RLF_Model& m, double t_max,
                            int M = 200, int Nt = 400) {
  if (!(m.v > 0.0)) {
    throw std::invalid_argument("rlf_solve: v must be positive.");
  }
  if (!(m.alpha > 1.0) || m.alpha > 2.0) {
    throw std::invalid_argument("rlf_solve: alpha must be in (1, 2].");
  }
  if (!(m.sigma > 0.0)) {
    throw std::invalid_argument("rlf_solve: sigma must be positive.");
  }
  if (!(m.b0 > 0.0)) {
    throw std::invalid_argument("rlf_solve: boundary b0 must be positive.");
  }
  if (m.z0 < 0.0 || m.z0 >= m.b0) {
    throw std::invalid_argument("rlf_solve: z0 must be in [0, b0).");
  }
  if (!(t_max > 0.0)) {
    throw std::invalid_argument("rlf_solve: t_max must be positive.");
  }

  M = std::max(M, 30);
  Nt = std::max(Nt, 50);

  constexpr double RLF_LOWER_SCALES = 12.0;
  constexpr double RLF_DOMAIN_FACTOR = 1.5;
  constexpr double RLF_GRID_FACTOR = 1.5;
  constexpr double RLF_DOMAIN_PDF_TOL = 0.01;
  constexpr double RLF_DOMAIN_CDF_TOL = 0.002;
  constexpr double RLF_SPATIAL_PDF_TOL = 0.05;
  constexpr double RLF_SPATIAL_CDF_TOL = 0.0125;
  constexpr int RLF_MAX_DOMAIN_REFINEMENTS = 3;
  constexpr int RLF_MAX_SPATIAL_REFINEMENTS = 3;

  const double stable_scale =
    m.sigma * std::pow(0.5 * t_max, 1.0 / m.alpha);
  double lower_extent =
    std::max(m.b0, RLF_LOWER_SCALES * stable_scale);
  int intervals = M;
  const int max_intervals =
    (M > std::numeric_limits<int>::max() / 4)
      ? std::numeric_limits<int>::max()
      : std::max(1200, 4 * M);

  RLF_Result accepted =
    rlf_solve_fixed_grid(m, t_max, intervals, Nt, lower_extent);
  const double initial_h = accepted.dx;

  // Most routine parameter points are comfortably inside the initial
  // 12-scale domain and have adequate resolution.  Avoid mandatory dense-LU
  // comparison solves in that regime; the returned flag makes this decision
  // explicit to callers.
  if (rlf_fast_path_safe(m, t_max, accepted, intervals)) {
    accepted.refinement_skipped = true;
    return accepted;
  }

  // First isolate finite-domain error: expand x_lo while preserving dx.  This
  // prevents ordinary grid refinement from masquerading as domain convergence.
  bool domain_converged = false;
  for (int refinement = 1;
       refinement <= RLF_MAX_DOMAIN_REFINEMENTS; ++refinement) {
    const double expanded_extent = lower_extent * RLF_DOMAIN_FACTOR;
    const int expanded_intervals =
      rlf_intervals_for_spacing(expanded_extent, m.b0, initial_h);
    if (expanded_intervals > max_intervals) break;

    RLF_Result candidate = rlf_solve_fixed_grid(
      m, t_max, expanded_intervals, Nt, expanded_extent);
    const RLF_Comparison error =
      compare_rlf_results(accepted, candidate, t_max);
    candidate.domain_pdf_error = error.pdf;
    candidate.domain_cdf_error = error.cdf;
    candidate.domain_refinements = refinement;
    accepted = std::move(candidate);
    lower_extent = expanded_extent;
    intervals = expanded_intervals;

    if (error.pdf <= RLF_DOMAIN_PDF_TOL &&
        error.cdf <= RLF_DOMAIN_CDF_TOL) {
      domain_converged = true;
      break;
    }
  }
  if (!domain_converged) {
    throw std::runtime_error(
      "rlf_solve: automatic lower-domain expansion failed to converge; "
      "increase nx or reduce the requested time/parameter range.");
  }

  // Then refine dx on the accepted domain.  Each candidate recomputes the
  // automatic time-step requirement, so spatial and temporal stability remain
  // coupled exactly as they are for the initial solve.
  bool spatial_converged = false;
  for (int refinement = 1;
       refinement <= RLF_MAX_SPATIAL_REFINEMENTS; ++refinement) {
    const double proposed =
      std::ceil(RLF_GRID_FACTOR * static_cast<double>(intervals));
    if (!(proposed <= static_cast<double>(max_intervals))) break;
    const int refined_intervals = static_cast<int>(proposed);

    RLF_Result candidate = rlf_solve_fixed_grid(
      m, t_max, refined_intervals, Nt, lower_extent);
    const RLF_Comparison error =
      compare_rlf_results(accepted, candidate, t_max);
    candidate.domain_pdf_error = accepted.domain_pdf_error;
    candidate.domain_cdf_error = accepted.domain_cdf_error;
    candidate.domain_refinements = accepted.domain_refinements;
    candidate.spatial_pdf_error = error.pdf;
    candidate.spatial_cdf_error = error.cdf;
    candidate.spatial_refinements = refinement;
    accepted = std::move(candidate);
    intervals = refined_intervals;

    if (error.pdf <= RLF_SPATIAL_PDF_TOL &&
        error.cdf <= RLF_SPATIAL_CDF_TOL) {
      spatial_converged = true;
      break;
    }
  }
  if (!spatial_converged) {
    throw std::runtime_error(
      "rlf_solve: automatic spatial refinement failed to converge; "
      "increase nx or use a solver suited to this high-resolution regime.");
  }

  accepted.refinement_checked = true;
  return accepted;
}

// ---------------------------------------------------------------------------
// Chambers-Mallows-Stuck symmetric alpha-stable random variable.
// Standard CMS has characteristic function exp(-|k|^alpha); the caller applies
// the 2^(-1/alpha) factor required by the EMC2 sigma convention.
// ---------------------------------------------------------------------------
inline double sample_stable_cms(double alpha, std::mt19937& rng) {
  if (alpha == 2.0) {
    std::normal_distribution<double> norm(0.0, std::sqrt(2.0));
    return norm(rng);
  }
  std::uniform_real_distribution<double> unif(
    -M_PI_2 + 1e-12, M_PI_2 - 1e-12);
  std::exponential_distribution<double> exp1(1.0);

  const double U = unif(rng);
  const double W = exp1(rng);
  const double aU = alpha * U;
  const double om_a_U = (1.0 - alpha) * U;
  const double num = std::sin(aU);
  const double den = std::pow(std::cos(U), 1.0 / alpha);
  const double term2 =
    std::pow(std::cos(om_a_U) / W, (1.0 - alpha) / alpha);
  return (num / den) * term2;
}

inline std::vector<double> simulate_rlf_hit_times(
    int n_sims, double v, double sigma, double alpha, double b0, double z0,
    double t_max = 5.0, double dt = 0.001, unsigned int seed = 42) {
  if (n_sims <= 0) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: n_sims must be positive.");
  }
  if (!(v > 0.0)) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: v must be positive.");
  }
  if (!(alpha > 1.0) || alpha > 2.0) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: alpha must be in (1, 2].");
  }
  if (!(sigma > 0.0)) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: sigma must be positive.");
  }
  if (!(b0 > 0.0)) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: b0 must be positive.");
  }
  if (z0 < 0.0 || z0 >= b0) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: z0 must be in [0, b0).");
  }
  if (!(t_max > 0.0)) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: t_max must be positive.");
  }
  if (!(dt > 0.0)) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: dt must be positive.");
  }

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif_start(0.0, z0);
  std::uniform_real_distribution<double> unif01(0.0, 1.0);

  const double n_steps_required = std::ceil(t_max / dt);
  if (!(n_steps_required <=
        static_cast<double>(std::numeric_limits<int>::max()))) {
    throw std::invalid_argument(
      "simulate_rlf_hit_times: t_max/dt requires too many steps.");
  }
  const int n_steps = static_cast<int>(n_steps_required);
  const double full_scale =
    sigma * std::pow(0.5 * dt, 1.0 / alpha);
  const double last_dt = t_max - (n_steps - 1) * dt;
  const double last_scale =
    sigma * std::pow(0.5 * last_dt, 1.0 / alpha);

  std::vector<double> hit_times(
    n_sims, std::numeric_limits<double>::quiet_NaN());

  for (int s = 0; s < n_sims; ++s) {
    double x = (z0 > 0.0) ? unif_start(rng) : 0.0;
    double time = 0.0;
    for (int step = 0; step < n_steps; ++step) {
      const bool final_step = step + 1 == n_steps;
      const double step_dt = final_step ? last_dt : dt;
      const double scale = final_step ? last_scale : full_scale;
      const double x_old = x;
      x += v * step_dt + scale * sample_stable_cms(alpha, rng);
      time += step_dt;

      bool crossed = x >= b0;
      if (!crossed && alpha == 2.0) {
        // Conditional Brownian-bridge crossing probability.  This removes the
        // O(sqrt(dt)) discrete-monitoring bias in the alpha=2 validation limit.
        const double exponent =
          -2.0 * (b0 - x_old) * (b0 - x) /
          (sigma * sigma * step_dt);
        crossed = unif01(rng) < std::exp(exponent);
      }
      if (crossed) {
        hit_times[s] = time;
        break;
      }
    }
  }
  return hit_times;
}

} // namespace rlf

#endif // model_RLF_h

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
//   * a time-invariant operator inverted once for each distinct time step;
//   * PDF from upper-boundary probability flux and CDF from surviving mass;
//   * an independent flux/mass mismatch convergence diagnostic.
// Likelihood evaluation adds dimensionless parameter caching, sparse query-time
// output, and 4/8-lane SIMD batches around this numerical core.
//
// Discretisation
// --------------
// Space is second order: the jump part uses fractional centred differences
// (not the first-order shifted Grünwald-Letnikov stencil) and the drift uses
// centred finite-volume faces (not donor-cell upwinding, whose v*h/2 numerical
// diffusion and v*p_n boundary flux dominated every other error source).
// Time is TR-BDF2, which is second order and L-stable, so the step size is set
// by accuracy rather than by a positivity limit dt ~ h^alpha; without that
// coupling, refining the grid costs O(n^2) per step instead of O(n^4) overall.
// A backward-Euler startup phase -- free, because backward Euler over
// tau = (gamma/2) dt reuses the TR-BDF2 factorisation -- damps the point-mass
// initial condition.
//
// Spatial truncation
// ------------------
// The physical domain is (-infinity, b0).  A finite lower edge is unavoidable
// numerically, but killing paths there would incorrectly count lower exits as
// upper first passages.  Attempts to jump below the lower edge are therefore
// censored (a conservative/no-flux closure).  The edge is placed at the depth
// whose first-passage probability matches the grid resolution, so domain width
// and discretisation error are balanced instead of the domain being fixed at a
// number of stable scales.  The solver can then automatically expand this
// lower domain at fixed spatial resolution, then refine the spatial grid,
// until independent PDF/CDF comparisons stabilize.  Refinement is bounded and
// a non-converged solve is rejected rather than silently returned.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <R_ext/RS.h>

#include "rlf_toeplitz.h"
#if defined(__x86_64__) || defined(_M_X64) || defined(__AVX2__)
#include <immintrin.h>
#endif

// R_ext/Lapack.h declares the entire LAPACK surface and clashes noisily with
// RcppArmadillo's declarations when this header is used by particle_ll.cpp.
// The RLF inverse only needs these two routines.
extern "C" {
void F77_NAME(dgetrf)(const int*, const int*, double*, const int*, int*, int*);
void F77_NAME(dgetri)(const int*, double*, const int*, const int*, double*,
                      const int*, int*);
void F77_NAME(dgemv)(const char*, const int*, const int*, const double*,
                     const double*, const int*, const double*, const int*,
                     const double*, double*, const int*, size_t);
void F77_NAME(dgemm)(const char*, const char*, const int*, const int*,
                     const int*, const double*, const double*, const int*,
                     const double*, const int*, const double*, double*,
                     const int*, size_t, size_t);
}

namespace rlf {

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Dense propagators are stored row-major.  BLAS sees the same bytes as the
// column-major transpose, so DGEMV with trans='T' computes the desired
// matrix-vector product without copying.  Keeping this in one helper also
// gives vendor BLAS implementations the whole operation instead of relying on
// the compiler to rediscover a tuned GEMV inside every time step.
inline void rlf_dense_matvec(const std::vector<double>& matrix,
                             const std::vector<double>& x,
                             std::vector<double>& out, int n) {
  out.resize(n);
  const char transpose = 'T';
  const int increment = 1;
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemv)(
    &transpose, &n, &n, &one, matrix.data(), &n, x.data(), &increment,
    &zero, out.data(), &increment, 1);
}

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

// Explicit inverse of the stationary Crank-Nicolson left-hand side.  A dense
// triangular solve and a dense matrix-vector product have the same O(n^2)
// operation count, but the latter has no loop-carried dependency and therefore
// vectorises well.  The inverse is formed once per distinct time-step block.
struct RLF_DenseInverse {
  int n = 0;
  std::vector<double> value;

  bool build(int n_, const std::vector<double>& A) {
    n = n_;
    const int nn = n;
    const int lda = nn;
    int info = 0;
    std::vector<double> column_major(static_cast<size_t>(n) * n);
    for (int row = 0; row < n; ++row) {
      for (int col = 0; col < n; ++col) {
        column_major[static_cast<size_t>(col) * n + row] =
          A[static_cast<size_t>(row) * n + col];
      }
    }
    std::vector<int> pivot(n);
    F77_CALL(dgetrf)(
      &nn, &nn, column_major.data(), &lda, pivot.data(), &info);
    if (info != 0) return false;

    double workspace_query = 0.0;
    const int query_size = -1;
    F77_CALL(dgetri)(
      &nn, column_major.data(), &lda, pivot.data(),
      &workspace_query, &query_size, &info);
    if (info != 0) return false;
    const int workspace_size = std::max(
      nn, static_cast<int>(workspace_query));
    std::vector<double> workspace(workspace_size);
    F77_CALL(dgetri)(
      &nn, column_major.data(), &lda, pivot.data(),
      workspace.data(), &workspace_size, &info);
    if (info != 0) return false;

    value.resize(static_cast<size_t>(n) * n);
    for (int row = 0; row < n; ++row) {
      for (int col = 0; col < n; ++col) {
        value[static_cast<size_t>(row) * n + col] =
          column_major[static_cast<size_t>(col) * n + row];
      }
    }
    return true;
  }

  void multiply(const std::vector<double>& x,
                std::vector<double>& out) const {
    rlf_dense_matvec(value, x, out, n);
  }
};

struct RLF_TimeSchedule {
  std::vector<double> dt;
  std::vector<int> steps;

  int n_steps() const {
    int out = 0;
    for (int n : steps) out += n;
    return out;
  }
};

// Piecewise-constant time grading, matching the schedule used by the local FPE
// backend.  Keeping dt constant within a block amortises the O(n^3) inverse
// setup, while doubling it between blocks concentrates resolution near t = 0.
inline RLF_TimeSchedule rlf_time_schedule(double T, int nt, double tgrade) {
  RLF_TimeSchedule s;
  const int n_steps = std::max(nt, 1);
  if (!(tgrade > 1.0 + 1e-9) || n_steps < 4) {
    s.dt.push_back(T / n_steps);
    s.steps.push_back(n_steps);
    return s;
  }

  int n_blk = 1 + static_cast<int>(
    std::floor(std::log(tgrade) / std::log(2.0) + 1e-9));
  n_blk = std::max(1, std::min(n_blk, n_steps / 2));
  const int per = std::max(2, n_steps / n_blk);
  const double dt0 =
    T / (per * (std::ldexp(1.0, n_blk) - 1.0));
  s.dt.reserve(n_blk);
  s.steps.reserve(n_blk);
  for (int b = 0; b < n_blk; ++b) {
    s.dt.push_back(std::ldexp(dt0, b));
    s.steps.push_back(per);
  }
  return s;
}

// Ceiling on the stability-driven step budget.  The budget below is set by the
// operator's exit rate, which scales as h^-alpha, and h is parameter dependent
// (rlf_lower_extent sizes the domain per parameter set) even though nx is not.
// A drift large enough to make the lower tail unreachable collapses the domain
// onto b0, and with it h, so the budget diverges: at alpha within 1e-8 of 2 and
// v/b ~ 6e6 a single solve wants ~1e9 steps.  INT_MAX is no protection -- that
// is hours of marching inside one likelihood evaluation, which presents as a
// silently wedged sampler rather than as an error.  Converged solves at the
// shipped defaults use 95 to ~2000 steps, so this leaves two orders of
// magnitude of headroom and still bails in well under a second.  Throwing hands
// the decision to the caller: from a likelihood, safe_new_particle rejects the
// particle, which is the right answer for a region this degenerate anyway.
inline int rlf_max_time_steps = 100000;

inline RLF_TimeSchedule rlf_stable_time_schedule(
    double t_max, int Nt, double max_exit_rate, double safe_cn,
    double tgrade) {
  const double step_ceiling =
    static_cast<double>(std::max(1, rlf_max_time_steps));
  const double nt_required =
    std::ceil(t_max * max_exit_rate / safe_cn);
  if (!(nt_required <= step_ceiling)) {
    throw std::invalid_argument(
      "rlf_solve: parameters require too many stable time steps.");
  }
  int step_budget = std::max(
    Nt, std::max(1, static_cast<int>(nt_required)));
  RLF_TimeSchedule schedule =
    rlf_time_schedule(t_max, step_budget, tgrade);
  for (int guard = 0; guard < 8; ++guard) {
    const double max_dt =
      *std::max_element(schedule.dt.begin(), schedule.dt.end());
    const double ratio = max_dt * max_exit_rate / safe_cn;
    if (!(ratio > 1.0 + 1e-12)) break;
    const double grown = std::ceil(step_budget * ratio * 1.01);
    if (!(grown <= step_ceiling)) {
      throw std::invalid_argument(
        "rlf_solve: graded schedule requires too many stable time steps.");
    }
    step_budget = std::max(step_budget + 1, static_cast<int>(grown));
    schedule = rlf_time_schedule(t_max, step_budget, tgrade);
  }
  return schedule;
}


// Fractional centred-difference weights,
//
//   c_k = (-1)^k Gamma(alpha+1) / [Gamma(alpha/2-k+1) Gamma(alpha/2+k+1)],
//
// for which -h^-alpha sum_k c_k u(x-kh) approximates the Riesz derivative
// d^alpha u / d|x|^alpha to second order in h.  The shifted Grünwald-Letnikov
// stencil this replaces is only first order, which dominated the solver error
// at every usable grid size.  For 1 < alpha <= 2 the weights satisfy c_0 > 0,
// c_{k != 0} < 0 and sum_k c_k = 0, so the generator below stays an M-matrix.
// At alpha = 2 they collapse to (-1, 2, -1), the standard Laplacian.
inline std::vector<double> compute_centred_weights(double alpha, int max_k) {
  std::vector<double> c(max_k + 1, 0.0);
  const double half = 0.5 * alpha;
  c[0] = std::exp(std::lgamma(alpha + 1.0) - 2.0 * std::lgamma(half + 1.0));
  for (int k = 0; k < max_k; ++k) {
    c[k + 1] = -c[k] * (half - k) / (half + k + 1.0);
  }
  return c;
}

struct RLF_Model {
  double v = 1.0;
  double sigma = 1.0;
  double alpha = 1.7;
  double b0 = 1.0;
  double z0 = 0.0;
};

struct RLF_Result {
  std::vector<double> t, pdf, cdf, surv;
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

inline double positive_log_interpolate(double x0, double x1, double w) {
  if (x0 > 0.0 && x1 > 0.0) {
    return std::exp(std::log(x0) + w * (std::log(x1) - std::log(x0)));
  }
  return x0 + w * (x1 - x0);
}

struct RLF_Operator {
  int n = 0;
  std::vector<double> L;
  std::vector<double> upper_kill;
  std::vector<double> lower_censor;
  double max_exit_rate = 0.0;
  double conservation_error = 0.0;
  bool upwind_drift = false;
};

// Symmetric centred-difference generator on the uniform grid.  The jump part
// has Fourier symbol -0.5*sigma^alpha*|k|^alpha up to O(h^2), matching the
// model convention and tending continuously to -0.5*sigma^2*k^2 at alpha = 2.
//
// The drift is discretised in flux form with centred face values
// F_{i+1/2} = v (p_i + p_{i+1}) / 2, which is second order and, crucially,
// consistent at the absorbing edge: the Dirichlet ghost p(b0) = 0 makes the
// The advective boundary flux is v*p_n/2.  Using the upwind value v*p_n adds
// numerical diffusion v*h/2 to sigma^2/2 and biases the exit flux.
//
// The Peclet fallback below is not a rare safety net.  The condition reduces to
// sigma^alpha |c_1(alpha)| h^(1-alpha) >= v, whose h exponent vanishes as
// alpha -> 1, leaving the grid-independent limit 0.424 sigma / v.  Below about
// alpha = 1.25 (about 1.5 when v = 2 sigma) it therefore fails at every nx, and
// refinement cannot restore second order.  Setting this flag keeps the centred
// faces regardless; positivity is then left to TR-BDF2's L-stability and the
// backward-Euler startup, and is watched by the min_density guards.
inline bool rlf_force_centred_drift = false;

// Diagnostic multiplier on the domain-width cap in rlf_lower_extent; 1.0 is the
// calibrated setting.  Lets the mesh and the truncation depth be varied
// independently, which nx alone cannot do because the cap grows with it.
inline double rlf_width_scale = 1.0;

inline RLF_Operator build_rlf_operator(const RLF_Model& m, int n, double h) {
  RLF_Operator op;
  op.n = n;
  op.L.assign(static_cast<size_t>(n) * n, 0.0);
  op.upper_kill.assign(n, 0.0);
  op.lower_censor.assign(n, 0.0);

  const auto c = compute_centred_weights(m.alpha, n + 1);
  const double scale =
    0.5 * std::pow(m.sigma, m.alpha) / std::pow(h, m.alpha);

  // q[d] is the transition rate between grid points d cells apart.
  // All q[d>0] are non-negative for 1 < alpha <= 2.
  std::vector<double> q(n, 0.0);
  for (int d = 0; d < n; ++d) q[d] = -scale * c[d];

  // tail[D] = sum_{d>=D} q[d] is the one-sided jump rate to targets at least D
  // cells away.  sum_{k=-inf}^{inf} c_k = 0 gives sum_{d>=1} c_d = -c_0/2, so
  //   sum_{d>=D} q[d] = scale * (c_0/2 + sum_{d=1}^{D-1} c_d),
  // which keeps jumps past the end of the grid exactly accounted for.
  std::vector<double> tail(n + 1, 0.0);
  double partial = 0.5 * c[0];
  for (int k = 1; k <= n; ++k) {
    tail[k] = std::max(0.0, scale * partial);
    partial += c[k];
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

  // Positive-drift transport, by exponential fitting of the nearest-neighbour
  // pair (Il'in / Scharfetter-Gummel).  Only |i-j| = 1 changes, so the operator
  // stays Toeplitz-plus-diagonal and the drift is uniform in j.
  //
  // A nearest-neighbour rate q[1] on spacing h is the diffusion D = q[1] h^2,
  // giving the cell Peclet number Pe = v h / D = v / (h q[1]).  Fitting the
  // constant-coefficient flux F = v p - D p' exactly on the grid gives
  //
  //   a+ = (v/h) / (1 - exp(-Pe)),   a- = (v/h) / (exp(Pe) - 1),
  //
  // with a+ - a- = v/h exactly, so the drift is transported without loss.  Both
  // rates are non-negative at every Pe, so L stays an M-matrix unconditionally.
  // The scheme interpolates the two branches this replaces: as Pe -> 0 it is the
  // centred pair q[1] +- v/(2h), and as Pe -> infinity it is donor-cell upwind.
  //
  // The gate it replaces (centred while q[1] >= v/(2h), i.e. Pe <= 2, else
  // upwind) was not a rare fallback.  Written out, the condition is
  // sigma^alpha |c_1(alpha)| h^(1-alpha) >= v, whose h exponent vanishes as
  // alpha -> 1, leaving the grid-independent limit 0.424 sigma / v.  Below about
  // alpha = 1.25 (about 1.5 when v = 2 sigma) it failed at every nx, so the
  // solver was silently first order there with no way to refine out of it: the
  // numerical diffusion v h / 2 made the process too diffusive, the CDF ran ~5%
  // above Monte Carlo at alpha = 1.1, and the resulting too-light tail biased
  // alpha-hat downwards.  Forcing the centred pair instead is not an option --
  // it produces material negative density at alpha = 1.1 for any nt -- whereas
  // fitting adds only the artificial diffusion positivity actually requires,
  // D (Pe/2 coth(Pe/2) - 1) against upwind's D Pe/2.
  const double adv = m.v / h;
  const double half_adv = 0.5 * adv;
  double delta_up = half_adv;    // a+ - q[1]
  double delta_down = -half_adv; // a- - q[1]
  op.upwind_drift = false;
  if (!(n > 1) || !(q[1] > 0.0)) {
    // No usable nearest-neighbour diffusion to fit against.
    op.upwind_drift = true;
    delta_up = adv;
    delta_down = 0.0;
  } else if (!rlf_force_centred_drift) {
    const double pe = adv / q[1];
    if (pe > 700.0) {              // exp(Pe) overflows; a- is zero to machine
      op.upwind_drift = true;
      delta_up = adv - q[1];
      delta_down = -q[1];
    } else if (pe > 1e-8) {        // below this the centred pair is exact
      delta_up = adv / -std::expm1(-pe) - q[1];
      delta_down = adv / std::expm1(pe) - q[1];
    }
  }

  for (int j = 0; j < n; ++j) {
    double outflow = delta_up;
    if (j + 1 < n) {
      op.L[static_cast<size_t>(j + 1) * n + j] += delta_up;
    } else {
      // The cell above the last one is past the barrier: killed, not moved.
      op.upper_kill[j] += delta_up;
    }
    if (j > 0) {
      op.L[static_cast<size_t>(j - 1) * n + j] += delta_down;
      outflow += delta_down;
    }
    // At j = 0 the downward face is the truncated lower edge, whose rate is
    // censored back out anyway, so only the upward change leaves the diagonal.
    op.L[static_cast<size_t>(j) * n + j] -= outflow;
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

// ---------------------------------------------------------------------------
// TR-BDF2 time integration.
//
// gamma = 2 - sqrt(2) is the classical choice for which the trapezoidal and
// BDF2 stages share the implicit matrix M = I - (gamma/2) dt L, so one
// factorisation serves both.  With u* the trapezoidal stage value,
//
//   u*   = (2 M^-1 - I) u_n,
//   u_n1 = M^-1 (a u* - b u_n),   a = 1/(gamma(2-gamma)), b = a - 1,
//
// hence u_n1 = T u_n with T = 2a M^-1 M^-1 - (a+b) M^-1.  The scheme is second
// order and L-stable, so stiff modes are damped instead of ringing.
// Crank-Nicolson (which this replaces) is only A-stable and needed both a
// Rannacher startup and dt <= 1.8/max_exit_rate ~ h^alpha to stay positive,
// which coupled the step count to the grid and made refinement cost O(n^4).
constexpr double RLF_TRBDF2_GAMMA = 0.5857864376269049;      // 2 - sqrt(2)
constexpr double RLF_TRBDF2_IMPLICIT = 0.2928932188134524;   // gamma/2
constexpr double RLF_TRBDF2_A = 1.2071067811865475;          // 1/(g(2-g))
constexpr double RLF_TRBDF2_B = 0.2071067811865475;          // a - 1

// I + (gamma/2) dt L stays entrywise non-negative below dt*max_exit_rate =
// 2/gamma, which is what makes the trapezoidal sub-stage positivity
// preserving.  Enforcing that outright would put dt ~ h^alpha and drive the
// refinement cost back to O(n^4), so the schedule instead damps the singular
// initial data with the backward-Euler startup below and keeps only a loose
// guard against genuinely pathological steps.
constexpr double RLF_TRBDF2_SAFE = 32.0 / RLF_TRBDF2_GAMMA;

// Backward Euler over tau = (gamma/2) dt uses exactly the TR-BDF2 matrix
// I - (gamma/2) dt L, so a Rannacher-style damping phase costs no extra
// factorisation.  Backward Euler is unconditionally positive on this
// M-matrix generator, which is what tames the point-mass initial condition
// before the (much cheaper) full steps take over.
constexpr int RLF_STARTUP_STEPS = 7;

// Rescale a schedule so that RLF_STARTUP_STEPS startup steps of
// tau = (gamma/2) dt[0] plus the scheduled full steps land exactly on t_max.
inline double rlf_apply_startup(RLF_TimeSchedule& schedule, double t_max) {
  if (schedule.dt.empty()) return 0.0;
  double total = 0.0;
  for (size_t b = 0; b < schedule.dt.size(); ++b) {
    total += schedule.dt[b] * schedule.steps[b];
  }
  const double startup_share =
    RLF_STARTUP_STEPS * RLF_TRBDF2_IMPLICIT * schedule.dt[0];
  const double factor = t_max / (total + startup_share);
  for (double& dt : schedule.dt) dt *= factor;
  return RLF_TRBDF2_IMPLICIT * schedule.dt[0];
}

inline std::vector<double> rlf_trbdf2_lhs(const RLF_Operator& op, double dt) {
  const int n = op.n;
  const double c = RLF_TRBDF2_IMPLICIT * dt;
  std::vector<double> lhs(static_cast<size_t>(n) * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const size_t ij = static_cast<size_t>(i) * n + j;
      lhs[ij] = (i == j ? 1.0 : 0.0) - c * op.L[ij];
    }
  }
  return lhs;
}

// Collapse both stages into one dense matrix so a step is a single matvec.
//
// The buffers are row-major, so asking BLAS (which is column-major) for
// minv * minv on the same storage returns (minv * minv)^T in column-major
// order, i.e. exactly minv * minv when read back row-major.
inline std::vector<double> rlf_trbdf2_step_matrix(
    const std::vector<double>& minv, int n) {
  std::vector<double> step(static_cast<size_t>(n) * n, 0.0);
  const char no_transpose = 'N';
  const double alpha = 2.0 * RLF_TRBDF2_A;
  const double beta = 0.0;
  F77_CALL(dgemm)(&no_transpose, &no_transpose, &n, &n, &n, &alpha,
                  minv.data(), &n, minv.data(), &n, &beta, step.data(), &n,
                  1, 1);
  const double linear = RLF_TRBDF2_A + RLF_TRBDF2_B;
  const size_t total = static_cast<size_t>(n) * n;
  for (size_t ij = 0; ij < total; ++ij) step[ij] -= linear * minv[ij];
  return step;
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

inline void rlf_initial_density(const RLF_Model& m, double x_lo, double h,
                                int n, std::vector<double>& p) {
  p.assign(n, 0.0);
  if (m.z0 > 0.0) {
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
    return;
  }

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

// Probability that the process ever drops below -L before t_max.
//
// The Brownian part is the exact drifted first-passage formula; the jump part
// is the single-big-jump estimate t * nu((-inf, -L]) with the Levy density
// nu(z) = c_alpha |z|^(-1-alpha), c_alpha = 0.5 sigma^alpha Gamma(1+alpha)
// sin(pi alpha/2)/pi, which vanishes continuously at alpha = 2.
inline double rlf_visit_probability(const RLF_Model& m, double t_max,
                                    double L) {
  const double norm_cdf = [](double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
  }(-(L + m.v * t_max) / (m.sigma * std::sqrt(t_max)));
  const double reflect = [](double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
  }((m.v * t_max - L) / (m.sigma * std::sqrt(t_max)));
  const double exponent = -2.0 * m.v * L / (m.sigma * m.sigma);
  const double brownian =
    norm_cdf + (exponent > -700.0 ? std::exp(exponent) * reflect : 0.0);

  const double c_alpha = 0.5 * std::pow(m.sigma, m.alpha) *
    std::tgamma(1.0 + m.alpha) * std::sin(M_PI * m.alpha / 2.0) / M_PI;
  const double jump = c_alpha > 0.0
    ? t_max * c_alpha * std::pow(L, -m.alpha) / m.alpha
    : 0.0;
  return brownian + jump;
}

// Initial lower truncation of the (-infinity, b0) domain.
//
// The censored closure at x_lo only misrepresents paths that actually reach
// it, so the extent is chosen as the depth whose visit probability equals a
// target tolerance -- not as a fixed number of stable scales, which ignored
// the drift and was roughly three times too wide in the Brownian limit.
//
// The tolerance is tied to the grid: every unit of extent is paid for in
// resolution (h = (b0 + L)/nx) and the discretisation error is O(h^2), so
// truncating at tol ~ 1/nx^2 balances the two error sources instead of
// letting an over-wide domain starve the boundary layer of grid points.
// Empirically this is a 4-5x accuracy gain at fixed cost.
inline double rlf_lower_extent(const RLF_Model& m, double t_max, int nx) {
  const double resolution = std::max(30.0, static_cast<double>(nx));
  const double tol =
    std::min(1e-2, std::max(1e-6, 100.0 / (resolution * resolution)));

  double hi = std::max(m.b0, m.sigma * std::pow(0.5 * t_max, 1.0 / m.alpha));
  for (int guard = 0; guard < 60 && rlf_visit_probability(m, t_max, hi) > tol;
       ++guard) {
    hi *= 2.0;
  }
  double lo = 0.0;
  for (int it = 0; it < 60; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (rlf_visit_probability(m, t_max, mid) > tol) lo = mid; else hi = mid;
  }

  // Heavy tails cannot be covered: the depth needed grows like tol^(-1/alpha),
  // so for small alpha the criterion above asks for a domain hundreds of units
  // wide, and every unit of it is paid for out of the boundary layer.  The cap
  // is where that trade is settled, so it is set by measurement rather than by
  // argument: at a fixed grid budget the CDF error against a wide converged
  // reference is minimised near six half-widths, and the optimum drifts wider
  // as the budget grows (alpha = 1.3, nx = 128: 3.0e-2 at two half-widths
  // versus 8.5e-3 at six, for exactly the same cost).  It never binds near
  // alpha = 2, where the visit-probability criterion is the tighter of the two.
  const double budget_widths =
    6.0 * std::sqrt(std::max(30.0, static_cast<double>(nx)) / 128.0);
  const double max_widths =
    rlf_width_scale * std::min(24.0, std::max(4.0, budget_widths));
  const double spread = m.sigma * std::pow(0.5 * t_max, 1.0 / m.alpha);
  const double widest = max_widths * std::max(m.b0, spread);
  return std::max(m.b0, std::min(hi, widest));
}

// Grid size past which the stages are solved matrix-free rather than by
// forming the explicit inverse.
//
// Set by measurement, not by flop count: the dense path is O(n^3) against
// O(iters * n log n), but it spends those flops inside OpenBLAS at tens of
// GFLOPS while the transforms here are hand-rolled and memory bound, which
// pushes the crossover far above where the asymptotics alone would put it.
// Measured against the dense path on identical grids (alpha = 1.5, 200
// queries): 0.12x at nx = 512, 0.28x at 1024, 0.81x at 2048, 1.32x at 3072
// and 1.74x at 4096, with the two paths agreeing to ~1e-11.  So this only
// pays for the refined grids that heavy tails need, and would be a large
// pessimisation at the shipped defaults.
constexpr int RLF_MATRIX_FREE_MIN_N = 2560;
// Tight enough that the Krylov residual is far below the discretisation error
// it is embedded in, so the march is not polluted by solver tolerance.
constexpr double RLF_MATRIX_FREE_TOL = 1e-12;

// One fixed-domain/fixed-grid solve.  `M` is the number of intervals and the
// M-1 interior nodes are x_lo+h, ..., b0-h.  The public solver below owns all
// validation and convergence refinement.
inline RLF_Result rlf_solve_fixed_grid(
    const RLF_Model& m, double t_max, int M, int Nt, double lower_extent,
    double tgrade = 1.0, bool explicit_inverse = true,
    const std::vector<double>* query_times = nullptr) {
  const double x_lo = -lower_extent;
  const double h = (m.b0 - x_lo) / M;
  const int n = M - 1;

  const RLF_Operator op = build_rlf_operator(m, n, h);

  // TR-BDF2 keeps its trapezoidal sub-stage positivity preserving while
  // I + (gamma/2) dt L is entrywise non-negative.  This is the only remaining
  // constraint linking dt to the grid, and it is nearly twice as permissive as
  // the Crank-Nicolson one it replaces.
  RLF_TimeSchedule schedule = rlf_stable_time_schedule(
    t_max, Nt, op.max_exit_rate, RLF_TRBDF2_SAFE, tgrade);
  const double startup_step = rlf_apply_startup(schedule, t_max);
  const int n_steps = schedule.n_steps() + RLF_STARTUP_STEPS;
  RLF_DenseLU lu;
  RLF_DenseInverse inverse;
  std::vector<double> step_matrix;
  bool use_collapsed_step = true;

  // Forming the explicit inverse and the TR-BDF2 step matrix costs O(n^3) and
  // dominates everything else once the grid is refined, which is exactly what
  // heavy tails ask for.  Past the crossover the stages are solved
  // matrix-free instead: O(n log n) per BiCGSTAB iteration, and the Strang
  // circulant preconditioner holds the iteration count at four to six
  // independently of n and alpha, so the cubic term disappears.
  const bool matrix_free = n >= RLF_MATRIX_FREE_MIN_N;
  rlf::RLF_ToeplitzSystem mf;
  if (matrix_free) explicit_inverse = false;

  std::vector<double> p;
  rlf_initial_density(m, x_lo, h, n, p);

  RLF_Result res;
  res.x_lo = x_lo;
  res.dx = h;
  res.nx_used = M;
  res.nt_used = n_steps;
  res.operator_conservation_error = op.conservation_error;
  res.min_density = density_min(p);

  const int n_levels = 1 + n_steps;
  res.t.reserve(n_levels);
  res.pdf.reserve(n_levels);
  res.cdf.reserve(n_levels);
  res.surv.reserve(n_levels);

  double time = 0.0;
  double flux_prev = weighted_rate(op.upper_kill, p, h);
  double lower_prev = weighted_rate(op.lower_censor, p, h);
  double cdf_flux = 0.0;
  double lower_pressure = 0.0;
  double survivor_prev = 1.0;
  double cdf_prev = 0.0;

  const bool sparse = query_times != nullptr;
  size_t query_pos = 0;
  std::vector<double> sorted_queries;
  if (sparse) {
    sorted_queries = *query_times;
    std::sort(sorted_queries.begin(), sorted_queries.end());
    sorted_queries.erase(
      std::unique(sorted_queries.begin(), sorted_queries.end()),
      sorted_queries.end());
    while (query_pos < sorted_queries.size() &&
           sorted_queries[query_pos] <= 0.0) {
      res.t.push_back(sorted_queries[query_pos]);
      res.pdf.push_back(0.0);
      res.cdf.push_back(0.0);
      res.surv.push_back(1.0);
      ++query_pos;
    }
  } else {
    res.t.push_back(time);
    res.pdf.push_back(flux_prev);
    res.cdf.push_back(0.0);
    res.surv.push_back(1.0);
  }

  std::vector<double> rhs(n, 0.0), next(n, 0.0), work(n, 0.0),
    stage(n, 0.0);

  auto record_step = [&](double step, bool backward_euler) {
    const double old_time = time;
    const double old_flux = flux_prev;
    const double old_survivor = survivor_prev;
    if (explicit_inverse) {
      if (backward_euler) {
        rlf_dense_matvec(inverse.value, p, next, n);
      } else if (use_collapsed_step) {
        rlf_dense_matvec(step_matrix, p, next, n);
      } else {
        // The collapsed TR--BDF2 propagator is
        //   2a M^-2 - (a+b) M^-1.
        // For a short block, two GEMVs are cheaper than forming M^-2 with a
        // GEMM.  This is algebraically identical to the precomputed path.
        rlf_dense_matvec(inverse.value, p, stage, n);
        rlf_dense_matvec(inverse.value, stage, next, n);
        for (int i = 0; i < n; ++i) {
          next[i] =
            2.0 * RLF_TRBDF2_A * next[i] -
            (RLF_TRBDF2_A + RLF_TRBDF2_B) * stage[i];
        }
      }
    } else if (matrix_free) {
      // Warm-starting from the current density costs nothing and saves an
      // iteration or two, since consecutive steps differ by O(dt).
      auto solve = [&](const double* b, std::vector<double>& out) {
        out = p;
        if (rlf::rlf_toeplitz_solve(mf, b, out.data(),
                                    RLF_MATRIX_FREE_TOL) < 0) {
          throw std::runtime_error(
            "rlf_solve: matrix-free stage solve did not converge.");
        }
      };
      if (backward_euler) {
        solve(p.data(), next);
      } else {
        apply_shifted(op, RLF_TRBDF2_IMPLICIT * step, p, rhs);
        solve(rhs.data(), stage);
        for (int i = 0; i < n; ++i) {
          rhs[i] = RLF_TRBDF2_A * stage[i] - RLF_TRBDF2_B * p[i];
        }
        solve(rhs.data(), next);
      }
    } else if (backward_euler) {
      lu.solve(p.data(), next.data(), work.data());
    } else {
      // Trapezoidal stage over gamma*step, then the BDF2 stage; both share the
      // factorisation of I - (gamma/2) step L.
      apply_shifted(op, RLF_TRBDF2_IMPLICIT * step, p, rhs);
      lu.solve(rhs.data(), stage.data(), work.data());
      for (int i = 0; i < n; ++i) {
        rhs[i] = RLF_TRBDF2_A * stage[i] - RLF_TRBDF2_B * p[i];
      }
      lu.solve(rhs.data(), next.data(), work.data());
    }

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
    cdf_flux += 0.5 * step * (flux_prev + flux);
    lower_pressure += 0.5 * step * (lower_prev + lower);
    flux_prev = flux;
    lower_prev = lower;

    const double mass_raw = density_mass(p, h);
    double survivor = std::min(1.0, std::max(0.0, mass_raw));
    if (survivor > survivor_prev) survivor = survivor_prev;
    survivor_prev = survivor;
    double cdf_mass = std::min(1.0, std::max(0.0, 1.0 - mass_raw));
    if (cdf_mass < cdf_prev) cdf_mass = cdf_prev;
    cdf_prev = cdf_mass;
    // Trapezoidal integration of the exit flux and the surviving mass agree
    // only to the order of the time discretisation: TR-BDF2 advances mass with
    // its own stage quadrature (weights 0.354, 0.354, 0.293 on t_n, the stage
    // time and t_n+1), which the endpoint trapezoid does not reproduce.  The
    // residual is therefore an O(dt^2) convergence diagnostic rather than the
    // exact identity the Crank-Nicolson march produced.
    res.flux_mass_mismatch =
      std::max(res.flux_mass_mismatch,
               std::abs((1.0 - mass_raw) - cdf_flux));

    if (sparse) {
      while (query_pos < sorted_queries.size() &&
             sorted_queries[query_pos] <= time) {
        const double tq = sorted_queries[query_pos];
        const double w = (tq - old_time) / (time - old_time);
        const double qpdf =
          positive_log_interpolate(old_flux, flux, w);
        const double qsurv =
          positive_log_interpolate(old_survivor, survivor, w);
        res.t.push_back(tq);
        res.pdf.push_back(std::max(0.0, qpdf));
        res.surv.push_back(std::min(1.0, std::max(0.0, qsurv)));
        res.cdf.push_back(1.0 - res.surv.back());
        ++query_pos;
      }
    } else {
      res.t.push_back(time);
      res.pdf.push_back(flux);
      res.cdf.push_back(cdf_mass);
      res.surv.push_back(survivor);
    }
  };

  // Build the left-hand side once per graded-time block.  Both TR-BDF2 stages
  // use I - (gamma/2) dt L, and with the explicit inverse available the two
  // stages collapse into one dense matvec per step.
  auto update_lhs_and_factor = [&](double step, int block_steps) {
    if (matrix_free) {
      mf = rlf::rlf_build_toeplitz_system(op.L, n,
                                          RLF_TRBDF2_IMPLICIT * step);
      return;
    }
    std::vector<double> lhs = rlf_trbdf2_lhs(op, step);
    if (explicit_inverse) {
      if (!inverse.build(n, lhs)) {
        throw std::runtime_error(
          "rlf_solve: inversion of the nonlocal operator failed.");
      }
      // One GEMM plus one GEMV per step crosses two GEMVs per step at roughly
      // half the matrix dimension on the supported BLAS backends.
      use_collapsed_step = 2 * block_steps >= n;
      if (use_collapsed_step) {
        step_matrix = rlf_trbdf2_step_matrix(inverse.value, n);
      } else {
        step_matrix.clear();
      }
    } else if (!lu.factor(n, lhs)) {
      throw std::runtime_error(
        "rlf_solve: factorization of the nonlocal operator failed.");
    }
  };

  for (size_t b = 0; b < schedule.dt.size(); ++b) {
    const double dt_b = schedule.dt[b];
    update_lhs_and_factor(dt_b, schedule.steps[b]);
    if (b == 0) {
      for (int k = 0; k < RLF_STARTUP_STEPS; ++k) {
        record_step(startup_step, true);
      }
    }
    for (int k = 0; k < schedule.steps[b]; ++k) record_step(dt_b, false);
  }

  if (sparse && !sorted_queries.empty()) {
    // The schedule reaches t_max analytically; floating accumulation can finish
    // a few ulps below it.  Match the full-grid endpoint clamp for such queries.
    while (query_pos < sorted_queries.size()) {
      res.t.push_back(sorted_queries[query_pos]);
      res.pdf.push_back(std::max(0.0, flux_prev));
      res.surv.push_back(survivor_prev);
      res.cdf.push_back(1.0 - survivor_prev);
      ++query_pos;
    }
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
  // The initial domain is now placed for boundary resolution rather than to
  // drive the censored mass to zero, so a few tenths of a percent of pressure
  // is the designed operating point, not a warning sign.
  constexpr double RLF_FAST_PRESSURE_TOL = 0.01;
  constexpr double RLF_FAST_DX_SCALE_TOL = 0.10;
  constexpr int RLF_FAST_MIN_INTERVALS = 160;
  if (M < RLF_FAST_MIN_INTERVALS) return false;
  if (!(base.lower_boundary_pressure <= RLF_FAST_PRESSURE_TOL)) return false;
  if (base.min_density < -1e-10 / base.dx) return false;
  if (base.flux_mass_mismatch > 1e-3) return false;
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
                            int M = 200, int Nt = 400, bool adaptive = true,
                            double tgrade = 1.0,
                            bool explicit_inverse = true,
                            const std::vector<double>* query_times = nullptr,
                            double lower_extent_override = -1.0) {
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

  constexpr double RLF_DOMAIN_FACTOR = 1.5;
  constexpr double RLF_GRID_FACTOR = 1.5;
  constexpr double RLF_DOMAIN_PDF_TOL = 0.01;
  constexpr double RLF_DOMAIN_CDF_TOL = 0.002;
  constexpr double RLF_SPATIAL_PDF_TOL = 0.05;
  constexpr double RLF_SPATIAL_CDF_TOL = 0.0125;
  // Expansion starts from the resolution-matched extent, which is deliberately
  // narrow for heavy tails, so the verifier needs more doublings to escape it.
  constexpr int RLF_MAX_DOMAIN_REFINEMENTS = 7;
  constexpr int RLF_MAX_SPATIAL_REFINEMENTS = 3;

  // An override pins the truncation depth, which is what lets two solves of
  // different M share a domain -- the only way to make their spacings differ
  // by exactly the ratio of their M, as Richardson extrapolation requires.
  double lower_extent = lower_extent_override > 0.0
    ? lower_extent_override
    : rlf_lower_extent(m, t_max, M);
  int intervals = M;
  const int max_intervals =
    (M > std::numeric_limits<int>::max() / 4)
      ? std::numeric_limits<int>::max()
      : std::max(1200, 4 * M);

  // A step that produces materially negative density is recoverable: the
  // schedule is refined and the march retried.  Throwing straight out of a
  // likelihood evaluation would abort a whole sampling run.
  auto solve_with_retry = [&](int M_use, double extent_use,
                              const std::vector<double>* queries) {
    int steps = Nt;
    for (int attempt = 0; attempt < 4; ++attempt) {
      try {
        return rlf_solve_fixed_grid(m, t_max, M_use, steps, extent_use,
                                    tgrade, explicit_inverse, queries);
      } catch (const std::runtime_error&) {
        // Refining past the stability ceiling would spend the budget the
        // ceiling exists to bound, so stop retrying once it is reached.
        const int ceiling = std::max(1, rlf_max_time_steps);
        if (attempt == 3 || steps >= ceiling) throw;
        steps = std::min(steps * 4, ceiling);
      }
    }
    throw std::runtime_error("rlf_solve: unreachable retry state.");
  };

  RLF_Result accepted =
    solve_with_retry(intervals, lower_extent,
                     adaptive ? nullptr : query_times);
  const double initial_h = accepted.dx;

  // Most routine parameter points are comfortably inside the initial
  // 12-scale domain and have adequate resolution.  Avoid mandatory dense-LU
  // comparison solves in that regime; the returned flag makes this decision
  // explicit to callers.
  if (!adaptive || rlf_fast_path_safe(m, t_max, accepted, intervals)) {
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

    RLF_Result candidate =
      solve_with_retry(expanded_intervals, expanded_extent, nullptr);
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

    RLF_Result candidate =
      solve_with_retry(refined_intervals, lower_extent, nullptr);
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
// Race-likelihood solve cache.
//
// sigma is absent from Key by construction.  Dividing the state by sigma gives
// the equivalent process dY = (v/sigma)dt + dL, with b and the uniform
// start-point range scaled the same way.  t0 is also absent: it only shifts the
// query time and never changes the PDE march.
// ---------------------------------------------------------------------------
constexpr double RLF_LOG_FLOOR = -700.0;

inline double rlf_safe_log(double x) {
  return (x > 0.0 && std::isfinite(x))
    ? std::max(std::log(x), RLF_LOG_FLOOR)
    : RLF_LOG_FLOOR;
}

struct Key {
  double v = 0.0;
  double alpha = 0.0;
  double b = 0.0;
  double A = 0.0;
  // Horizon bucket; see rlf_horizon_bucket.  Rows that need a long march are
  // solved separately from the bulk so they cannot coarsen it.
  int bucket = 0;

  bool operator==(const Key& other) const {
    return v == other.v && alpha == other.alpha &&
           b == other.b && A == other.A && bucket == other.bucket;
  }
};

// A key's rows all share one march, sized to the longest first-passage time
// among them, and rlf_lower_extent sizes the lower domain from that horizon as
// L = widths * max(b, (0.5 t_max)^(1/alpha)).  While the stable spread stays
// under b the domain is pinned at its floor and, at fixed nx, the absorbing
// boundary keeps a fixed share of the grid.  Past
//
//   t_crit = 2 * b^alpha
//
// the spread takes over, the domain grows like t^(1/alpha) -- nearly linearly
// for alpha near 1 -- and the boundary layer loses resolution in proportion.
// One slow trial in a heavy-tailed data set is therefore enough to unresolve
// the likelihood of every other trial sharing its parameters: at alpha = 1.1,
// b = 2, a single 16.9 s observation drops the cells between the start point
// and the boundary from 18 to 6, and the resulting profile likelihood peaks at
// alpha ~ 1.4 and does not converge under grid refinement.
//
// Splitting the rows by horizon costs one extra march per occupied bucket and
// keeps the bulk of the data at the floor domain.  Buckets double, so bucket k
// widens the domain by 2^(k/alpha) rather than letting a single outlier set it
// for everyone.
//
// Re-measured after Richardson extrapolation landed, in case correcting the
// leading h error made the split redundant.  It does not: extrapolation assumes
// an asymptotic error expansion, and a ballooned domain leaves too few cells in
// the boundary layer for one to hold.  Turning the split off on 2000-trial data
// (profile peak, converged solver in brackets):
//
//   alpha 1.10 [1.111]   raw 128  1.010 -> 1.330   pair 128+192  1.092 -> 1.245
//   alpha 1.50 [1.581]   raw 128  1.508 -> 1.517   pair  64+96   1.619 -> 1.568
//
// and the profile log-likelihood at alpha = 1.3 goes from 0.054 to 0.478 nats
// of roughness about a local quadratic.  The cost is self-limiting: the split
// only engages once t_max passes t_crit, so at alpha = 1.7 with 5 s of data no
// bucket is occupied and it is a no-op, and the 1.6-1.7x it costs at alpha 1.1
// to 1.5 is only paid where it is also worth 0.15 to 0.32 of alpha.
//
// Loosening the threshold to t_crit = 2 (c b)^alpha, so the bulk group absorbs
// a factor c of widening before anything splits off, was measured at c = 1.3,
// 1.6 and 2.0.  It saves 5-20% and is erratic in both directions (alpha = 1.3
// peak error -0.011 -> -0.002 at c = 1.3 but -0.121 at c = 1.6), so the
// threshold stays where the domain first starts to grow.
constexpr int RLF_MAX_HORIZON_BUCKET = 6;

inline int rlf_horizon_bucket(const Key& key, double t_max) {
  const double t_crit = 2.0 * std::pow(key.b, key.alpha);
  if (!(t_crit > 0.0) || !(t_max > t_crit)) return 0;
  const int bucket =
    1 + static_cast<int>(std::log2(t_max / t_crit));
  return std::min(std::max(bucket, 1), RLF_MAX_HORIZON_BUCKET);
}

struct KeyHash {
  size_t operator()(const Key& key) const noexcept {
    size_t h = std::hash<double>{}(key.v);
    auto mix = [&](double x) {
      const size_t hx = std::hash<double>{}(x);
      h ^= hx + static_cast<size_t>(0x9e3779b9U) + (h << 6) + (h >> 2);
    };
    mix(key.alpha);
    mix(key.b);
    mix(key.A);
    mix(static_cast<double>(key.bucket));
    return h;
  }
};

inline bool rlf_key(double v, double sigma, double alpha, double B, double A,
                    Key& out) {
  if (!(v > 0.0) || !std::isfinite(v) ||
      !(sigma > 0.0) || !std::isfinite(sigma) ||
      !(alpha > 1.0) || alpha > 2.0 || !std::isfinite(alpha) ||
      !(B > 0.0) || !std::isfinite(B) ||
      !(A >= 0.0) || !std::isfinite(A)) {
    return false;
  }
  const double inv_sigma = 1.0 / sigma;
  out.v = v * inv_sigma;
  out.alpha = alpha;
  out.A = A * inv_sigma;
  out.b = (B + A) * inv_sigma;
  return std::isfinite(out.v) && std::isfinite(out.b) &&
         std::isfinite(out.A) && out.b > out.A;
}

// Keep these defaults synchronized with R/model_RLF.R.  `nx` is a fit-level
// option rather than a function of alpha: varying resolution with alpha adds
// parameter-dependent discretization bias and roughness to the likelihood.
struct Grid {
  int nx = 160;
  double dt_target = 1.6e-2;
  int nt_min = 50;
  int nt_max = 20000;
  double tgrade = 1.0;
  bool adaptive = false;
  bool explicit_inverse = true;
  bool sparse_output = true;
  // Solve long-horizon rows separately from the bulk; see rlf_horizon_bucket.
  bool horizon_split = true;
  // Pair each solve with a finer one and Richardson-extrapolate; see
  // rlf_cache_solve.
  bool richardson = true;
  double richardson_ratio = 1.25;
  // SIMD batching is optional; keep it disabled unless explicitly requested.
  bool simd_batch = false;

  int nt_for(double t_max) const {
    const double wanted =
      std::ceil(t_max / std::max(dt_target, 1e-8));
    const int bounded = wanted < static_cast<double>(nt_max)
      ? static_cast<int>(wanted) : nt_max;
    return std::max(nt_min, bounded);
  }
};

struct Entry {
  Key key;
  double t_max = 0.0;
  std::vector<double> t;
  std::vector<double> log_pdf;
  std::vector<double> log_S;
  bool complete_grid = true;
};

struct SolveCache {
  Grid grid;
  std::vector<Entry> entries;
  std::unordered_map<Key, int, KeyHash> index;
  std::vector<int> row_group;
  size_t solve_count = 0;

  void new_particle() {
    entries.clear();
    index.clear();
    row_group.clear();
    solve_count = 0;
  }
};

#if defined(__AVX512F__) && defined(__FMA__)
constexpr size_t RLF_BATCH_LANES = 8;
template <size_t LANES> struct RLFSimd;
template <> struct RLFSimd<8> {
  using Vec = __m512d;
  static inline Vec zero() { return _mm512_setzero_pd(); }
  static inline Vec load(const double* x) { return _mm512_loadu_pd(x); }
  static inline void store(double* x, Vec v) { _mm512_storeu_pd(x, v); }
  static inline Vec fmadd(Vec a, Vec b, Vec c) {
    return _mm512_fmadd_pd(a, b, c);
  }
};
template <> struct RLFSimd<4> {
  using Vec = __m256d;
  static inline Vec zero() { return _mm256_setzero_pd(); }
  static inline Vec load(const double* x) { return _mm256_loadu_pd(x); }
  static inline void store(double* x, Vec v) { _mm256_storeu_pd(x, v); }
  static inline Vec fmadd(Vec a, Vec b, Vec c) {
    return _mm256_fmadd_pd(a, b, c);
  }
};
#elif defined(__AVX2__) && defined(__FMA__)
constexpr size_t RLF_BATCH_LANES = 4;
template <size_t LANES> struct RLFSimd;
template <> struct RLFSimd<4> {
  using Vec = __m256d;
  static inline Vec zero() { return _mm256_setzero_pd(); }
  static inline Vec load(const double* x) { return _mm256_loadu_pd(x); }
  static inline void store(double* x, Vec v) { _mm256_storeu_pd(x, v); }
  static inline Vec fmadd(Vec a, Vec b, Vec c) {
    return _mm256_fmadd_pd(a, b, c);
  }
};
#else
constexpr size_t RLF_BATCH_LANES = 4;
#endif

inline void rlf_store_result(const Key& key, double t_max,
                             const RLF_Result& result, bool complete_grid,
                             Entry& out) {
  out.key = key;
  out.t_max = t_max;
  out.complete_grid = complete_grid;
  out.t = result.t;
  const size_t n = result.t.size();
  out.log_pdf.resize(n);
  out.log_S.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.log_pdf[i] = rlf_safe_log(result.pdf[i]);
    const double survivor =
      i < result.surv.size() ? result.surv[i] : 1.0 - result.cdf[i];
    out.log_S[i] = rlf_safe_log(survivor);
  }
}

// Richardson extrapolation in h.
//
// With the fitted drift in place the remaining error is the absorbing
// boundary: the exact density behaves like (b0 - x)^(alpha/2) there, and
// uniform-grid centred fractional weights cannot resolve a half-integer power
// at second order.  The result is a clean first-order error for every
// alpha < 2 -- measured order 0.85 to 1.4 across alpha in [1.1, 1.7], v in
// [0.5, 3], b0 in [1, 2] -- which is exactly the situation extrapolation is
// for.  (At alpha = 2 the exponent is 1, the singularity disappears and the
// measured order is 2.01; p = 1 then under-corrects rather than over-corrects,
// so it still helps there, just less than a p = 2 combination would.)
//
// Two solves at spacings h_c > h_f give u = (r u_f - u_c) / (r - 1) with
// r = h_c / h_f.  Both share one truncation depth, which matters twice over:
// it makes r exactly the ratio of the interval counts, and it holds the
// finite-domain error common between the two so the combination cannot amplify
// it.  Choosing the depth by nx instead would be actively wrong -- the width
// cap grows like sqrt(nx), so a fine solve can land on a *larger* h than the
// coarse one (nx = 256 gives 0.0741 against nx = 128's 0.0727) and the
// extrapolation then diverges.  The time step is likewise shared, so the
// O(dt^2) part of the error is common-mode and passes through unchanged.
//
// Measured on mean |d log pdf| per trial over the central 96% of each
// distribution, against an nx = 1024 reference on the pinned domain: 0.0220
// raw at nx = 128 against 0.0048 for the (128, 192) pair, a 4.6x reduction
// for 3.6x the time, and it improves every cell of the grid.  Refining
// instead is far worse value -- raw nx = 288 costs 2.3x the pair and is still
// 1.7x less accurate.
inline bool rlf_extrapolate(const RLF_Result& coarse, const RLF_Result& fine,
                            double ratio, RLF_Result& out) {
  const size_t n = fine.t.size();
  if (n == 0 || coarse.t.size() != n || !(ratio > 1.0)) return false;
  for (size_t i = 0; i < n; ++i) {
    // The two solves were handed the same query vector, so any drift here
    // means they are not on a common time base and cannot be combined.
    if (!(std::abs(coarse.t[i] - fine.t[i]) <=
          1e-9 * (1.0 + std::abs(fine.t[i])))) {
      return false;
    }
  }

  out = fine;
  if (out.surv.size() != n) {          // sparse solves may omit the survivor
    out.surv.resize(n);
    for (size_t i = 0; i < n; ++i) out.surv[i] = 1.0 - out.cdf[i];
  }
  const double w = ratio / (ratio - 1.0);
  auto blend = [&](double c, double f) { return w * f - (w - 1.0) * c; };
  for (size_t i = 0; i < n; ++i) {
    const double pdf = blend(coarse.pdf[i], fine.pdf[i]);
    const double c_surv =
      i < coarse.surv.size() ? coarse.surv[i] : 1.0 - coarse.cdf[i];
    const double f_surv =
      i < fine.surv.size() ? fine.surv[i] : 1.0 - fine.cdf[i];
    const double surv = blend(c_surv, f_surv);
    // Extrapolation is an unconstrained linear combination and can in
    // principle overshoot into a negative density or survivor.  It did not do
    // so anywhere in the parameter sweep, but a single negative value would
    // become a -Inf log-likelihood, so keep the fine solve pointwise wherever
    // the combination is not a usable density.
    if (!(pdf > 0.0) || !std::isfinite(pdf) ||
        !(surv > 0.0) || !std::isfinite(surv) || surv > 1.0) {
      continue;
    }
    out.pdf[i] = pdf;
    out.surv[i] = surv;
    out.cdf[i] = 1.0 - surv;
  }
  return true;
}

inline void rlf_cache_solve(const Key& key, double t_max, const Grid& grid,
                            const std::vector<double>* query_times,
                            Entry& out) {
  RLF_Model model;
  model.v = key.v;
  model.sigma = 1.0;
  model.alpha = key.alpha;
  model.b0 = key.b;
  model.z0 = key.A;

  const int nx = std::max(grid.nx, 30);
  const int nt = grid.nt_for(t_max);
  const int nx_fine =
    static_cast<int>(std::lround(nx * grid.richardson_ratio));

  if (grid.richardson && !grid.adaptive && nx_fine > nx) {
    const double extent = rlf_lower_extent(model, t_max, nx);
    const RLF_Result fine = rlf_solve(
      model, t_max, nx_fine, nt, false, grid.tgrade, grid.explicit_inverse,
      query_times, extent);
    // The two solves march their own time schedules -- the stable step
    // depends on the operator's exit rate, which changes with h -- so on the
    // complete-grid path they would land on different time bases and could
    // not be combined.  Handing the coarse solve the fine grid's own times
    // puts them back on a common base, so sparse and complete output stay
    // exactly equivalent rather than differing by whether extrapolation ran.
    const RLF_Result coarse = rlf_solve(
      model, t_max, nx, nt, false, grid.tgrade, grid.explicit_inverse,
      query_times != nullptr ? query_times : &fine.t, extent);
    RLF_Result blended;
    const bool ok = rlf_extrapolate(
      coarse, fine, static_cast<double>(nx_fine) / nx, blended);
    rlf_store_result(key, t_max, ok ? blended : fine,
                     query_times == nullptr, out);
    return;
  }

  const RLF_Result result = rlf_solve(
    model, t_max, nx, nt, grid.adaptive,
    grid.tgrade, grid.explicit_inverse, query_times);
  rlf_store_result(key, t_max, result, query_times == nullptr, out);
}

struct BatchAction {
  double step = 0.0;
  size_t block = 0;
  bool backward_euler = false;
};

struct BatchLane {
  RLF_Model model;
  RLF_Operator op;
  double h = 0.0;
  std::vector<double> initial;
  std::vector<std::vector<double>> step_matrix;
  std::vector<double> startup_matrix;
  std::vector<BatchAction> actions;
  std::vector<double> queries;
  RLF_Result result;
};

inline BatchLane rlf_build_batch_lane(
    const Key& key, double t_max, const Grid& grid,
    const std::vector<double>* query_times) {
  BatchLane lane;
  lane.model.v = key.v;
  lane.model.sigma = 1.0;
  lane.model.alpha = key.alpha;
  lane.model.b0 = key.b;
  lane.model.z0 = key.A;

  const double lower_extent =
    rlf_lower_extent(lane.model, t_max, std::max(grid.nx, 30));
  const int M = std::max(grid.nx, 30);
  const int n = M - 1;
  const double x_lo = -lower_extent;
  lane.h = (lane.model.b0 - x_lo) / M;
  lane.op = build_rlf_operator(lane.model, n, lane.h);
  rlf_initial_density(lane.model, x_lo, lane.h, n, lane.initial);

  RLF_TimeSchedule schedule = rlf_stable_time_schedule(
    t_max, grid.nt_for(t_max), lane.op.max_exit_rate, RLF_TRBDF2_SAFE,
    grid.tgrade);
  const double startup_step = rlf_apply_startup(schedule, t_max);
  lane.step_matrix.resize(schedule.dt.size());
  for (size_t block = 0; block < schedule.dt.size(); ++block) {
    std::vector<double> lhs = rlf_trbdf2_lhs(lane.op, schedule.dt[block]);
    RLF_DenseInverse inverse;
    if (!inverse.build(n, lhs)) {
      throw std::runtime_error(
        "rlf_solve: batched inversion of the nonlocal operator failed.");
    }
    if (block == 0) lane.startup_matrix = inverse.value;
    lane.step_matrix[block] = rlf_trbdf2_step_matrix(inverse.value, n);
  }

  for (int step = 0; step < RLF_STARTUP_STEPS; ++step) {
    lane.actions.push_back({startup_step, 0, true});
  }
  for (size_t block = 0; block < schedule.dt.size(); ++block) {
    for (int step = 0; step < schedule.steps[block]; ++step) {
      lane.actions.push_back({schedule.dt[block], block, false});
    }
  }

  if (query_times != nullptr) {
    lane.queries = *query_times;
    std::sort(lane.queries.begin(), lane.queries.end());
    lane.queries.erase(
      std::unique(lane.queries.begin(), lane.queries.end()),
      lane.queries.end());
    lane.result.t.reserve(lane.queries.size());
    lane.result.pdf.reserve(lane.queries.size());
    lane.result.cdf.reserve(lane.queries.size());
    lane.result.surv.reserve(lane.queries.size());
  } else {
    lane.result.t.reserve(lane.actions.size() + 1);
    lane.result.pdf.reserve(lane.actions.size() + 1);
    lane.result.cdf.reserve(lane.actions.size() + 1);
    lane.result.surv.reserve(lane.actions.size() + 1);
  }
  return lane;
}

template <size_t LANES>
inline std::vector<RLF_Result> rlf_solve_batch_lanes(
    const std::vector<Key>& keys, const std::vector<double>& horizons,
    const Grid& grid,
    const std::vector<std::vector<double>>* query_times = nullptr) {
  const size_t num_lanes = keys.size();
  std::vector<BatchLane> lanes;
  lanes.reserve(num_lanes);
  size_t max_actions = 0;
  for (size_t lane = 0; lane < num_lanes; ++lane) {
    lanes.push_back(rlf_build_batch_lane(
      keys[lane], horizons[lane], grid,
      query_times == nullptr ? nullptr : &(*query_times)[lane]));
    max_actions = std::max(max_actions, lanes.back().actions.size());
  }

  const int n = std::max(grid.nx, 30) - 1;
  std::vector<double> state(static_cast<size_t>(n) * LANES, 0.0);
  std::vector<double> next(static_cast<size_t>(n) * LANES, 0.0);
  std::vector<double> matrix(
    static_cast<size_t>(n) * n * LANES, 0.0);
  const std::vector<double>* selected[LANES] = {};
  const std::vector<double>* previous_selected[LANES] = {};
  double time[LANES] = {};
  double pdf[LANES] = {};
  double survivor[LANES] = {};
  size_t query_pos[LANES] = {};

  for (size_t lane = 0; lane < LANES; ++lane) {
    survivor[lane] = 1.0;
    if (lane >= num_lanes) continue;
    for (int i = 0; i < n; ++i) {
      state[static_cast<size_t>(i) * LANES + lane] =
        lanes[lane].initial[i];
    }
    pdf[lane] = weighted_rate(
      lanes[lane].op.upper_kill, lanes[lane].initial, lanes[lane].h);
    if (query_times == nullptr) {
      lanes[lane].result.t.push_back(0.0);
      lanes[lane].result.pdf.push_back(pdf[lane]);
      lanes[lane].result.cdf.push_back(0.0);
      lanes[lane].result.surv.push_back(1.0);
    } else {
      while (query_pos[lane] < lanes[lane].queries.size() &&
             lanes[lane].queries[query_pos[lane]] <= 0.0) {
        lanes[lane].result.t.push_back(
          lanes[lane].queries[query_pos[lane]]);
        lanes[lane].result.pdf.push_back(0.0);
        lanes[lane].result.cdf.push_back(0.0);
        lanes[lane].result.surv.push_back(1.0);
        ++query_pos[lane];
      }
    }
  }

  for (size_t action_index = 0;
       action_index < max_actions; ++action_index) {
    bool matrix_changed = false;
    bool active[LANES] = {};
    for (size_t lane = 0; lane < LANES; ++lane) {
      active[lane] =
        lane < num_lanes && action_index < lanes[lane].actions.size();
      if (active[lane]) {
        const BatchAction& action = lanes[lane].actions[action_index];
        selected[lane] = action.backward_euler
          ? &lanes[lane].startup_matrix
          : &lanes[lane].step_matrix[action.block];
      } else {
        selected[lane] = nullptr;
      }
      if (selected[lane] != previous_selected[lane]) matrix_changed = true;
    }

    if (matrix_changed) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          const size_t ij = static_cast<size_t>(i) * n + j;
          const size_t base = ij * LANES;
          for (size_t lane = 0; lane < LANES; ++lane) {
            matrix[base + lane] = selected[lane] == nullptr
              ? (i == j ? 1.0 : 0.0)
              : (*selected[lane])[ij];
          }
        }
      }
      for (size_t lane = 0; lane < LANES; ++lane) {
        previous_selected[lane] = selected[lane];
      }
    }

    for (int i = 0; i < n; ++i) {
#if (defined(__AVX512F__) || defined(__AVX2__)) && defined(__FMA__)
      using Simd = RLFSimd<LANES>;
      using Vec = typename Simd::Vec;
      Vec sum = Simd::zero();
      for (int j = 0; j < n; ++j) {
        const Vec coefficient = Simd::load(
          &matrix[(static_cast<size_t>(i) * n + j) * LANES]);
        const Vec value = Simd::load(
          &state[static_cast<size_t>(j) * LANES]);
        sum = Simd::fmadd(coefficient, value, sum);
      }
      Simd::store(&next[static_cast<size_t>(i) * LANES], sum);
#else
      for (size_t lane = 0; lane < LANES; ++lane) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
          sum += matrix[
            (static_cast<size_t>(i) * n + j) * LANES + lane] *
            state[static_cast<size_t>(j) * LANES + lane];
        }
        next[static_cast<size_t>(i) * LANES + lane] = sum;
      }
#endif
    }
    state.swap(next);

    for (size_t lane = 0; lane < num_lanes; ++lane) {
      if (!active[lane]) continue;
      const BatchAction& action = lanes[lane].actions[action_index];
      const double old_time = time[lane];
      const double old_pdf = pdf[lane];
      const double old_survivor = survivor[lane];
      time[lane] += action.step;

      double mass = 0.0;
      double flux = 0.0;
      double min_density = std::numeric_limits<double>::infinity();
      for (int i = 0; i < n; ++i) {
        const double value =
          state[static_cast<size_t>(i) * LANES + lane];
        mass += value;
        flux += lanes[lane].op.upper_kill[i] * value;
        min_density = std::min(min_density, value);
      }
      if (min_density < -1e-8 / lanes[lane].h) {
        throw std::runtime_error(
          "rlf_solve: material negative density in SIMD batch.");
      }
      mass *= lanes[lane].h;
      flux *= lanes[lane].h;
      if (!(flux > 0.0)) flux = 0.0;
      survivor[lane] =
        std::min(old_survivor, std::min(1.0, std::max(0.0, mass)));
      pdf[lane] = flux;

      if (query_times == nullptr) {
        lanes[lane].result.t.push_back(time[lane]);
        lanes[lane].result.pdf.push_back(pdf[lane]);
        lanes[lane].result.surv.push_back(survivor[lane]);
        lanes[lane].result.cdf.push_back(1.0 - survivor[lane]);
      } else {
        while (query_pos[lane] < lanes[lane].queries.size() &&
               lanes[lane].queries[query_pos[lane]] <= time[lane]) {
          const double query = lanes[lane].queries[query_pos[lane]];
          const double w =
            (query - old_time) / (time[lane] - old_time);
          const double query_pdf =
            positive_log_interpolate(old_pdf, pdf[lane], w);
          const double query_survivor = positive_log_interpolate(
            old_survivor, survivor[lane], w);
          lanes[lane].result.t.push_back(query);
          lanes[lane].result.pdf.push_back(std::max(0.0, query_pdf));
          lanes[lane].result.surv.push_back(
            std::min(1.0, std::max(0.0, query_survivor)));
          lanes[lane].result.cdf.push_back(
            1.0 - lanes[lane].result.surv.back());
          ++query_pos[lane];
        }
      }
    }
  }

  std::vector<RLF_Result> results(num_lanes);
  for (size_t lane = 0; lane < num_lanes; ++lane) {
    while (query_times != nullptr &&
           query_pos[lane] < lanes[lane].queries.size()) {
      lanes[lane].result.t.push_back(
        lanes[lane].queries[query_pos[lane]]);
      lanes[lane].result.pdf.push_back(pdf[lane]);
      lanes[lane].result.surv.push_back(survivor[lane]);
      lanes[lane].result.cdf.push_back(1.0 - survivor[lane]);
      ++query_pos[lane];
    }
    results[lane] = std::move(lanes[lane].result);
  }
  return results;
}

inline std::vector<RLF_Result> rlf_solve_batch(
    const std::vector<Key>& keys, const std::vector<double>& horizons,
    const Grid& grid,
    const std::vector<std::vector<double>>* query_times = nullptr) {
#if defined(__AVX512F__) && defined(__FMA__)
  if (keys.size() <= 4) {
    return rlf_solve_batch_lanes<4>(
      keys, horizons, grid, query_times);
  }
  return rlf_solve_batch_lanes<8>(
    keys, horizons, grid, query_times);
#else
  return rlf_solve_batch_lanes<4>(
    keys, horizons, grid, query_times);
#endif
}

inline bool entry_has_queries(const Entry& entry,
                              const std::vector<double>& query_times) {
  if (entry.complete_grid) return true;
  size_t j = 0;
  for (double time : query_times) {
    while (j < entry.t.size() && entry.t[j] < time) ++j;
    if (j >= entry.t.size() || entry.t[j] != time) return false;
  }
  return true;
}

// Returns -1 when this parameter point has no usable solve.  rlf_solve()
// already refines the time grid up to the stability ceiling before giving up,
// so an exception escaping it means the point is genuinely unresolvable (the
// alpha -> 2 / large-drift corner).  Propagating that out of a likelihood
// evaluation would kill the whole chain -- and at start-point selection, where
// particles are drawn from the prior, a single such draw among a thousand was
// enough to abort the fit.  Callers floor the affected rows instead, which
// rejects the particle.  Solving into a temporary keeps a failed refinement
// from leaving a half-written entry behind.
inline int cache_get(SolveCache& cache, const Key& key, double t_need) {
  if (!(t_need > 0.0)) t_need = 1e-3;
  const auto found = cache.index.find(key);
  if (found != cache.index.end()) {
    Entry& entry = cache.entries[found->second];
    if (entry.complete_grid && entry.t_max >= t_need) return found->second;
    Entry replacement;
    try {
      rlf_cache_solve(key, std::max(t_need, entry.t_max), cache.grid, nullptr,
                      replacement);
    } catch (const std::exception&) {
      return -1;
    }
    entry = std::move(replacement);
    ++cache.solve_count;
    return found->second;
  }

  Entry fresh;
  try {
    rlf_cache_solve(key, t_need, cache.grid, nullptr, fresh);
  } catch (const std::exception&) {
    return -1;
  }
  const int idx = static_cast<int>(cache.entries.size());
  cache.entries.push_back(std::move(fresh));
  cache.index.emplace(key, idx);
  ++cache.solve_count;
  return idx;
}

inline void cache_get_batch(
    SolveCache& cache, const std::vector<Key>& keys,
    const std::vector<double>& horizons, std::vector<int>& out_indices,
    const std::vector<std::vector<double>>* query_times = nullptr) {
  const size_t n = keys.size();
  out_indices.assign(n, -1);
  std::vector<size_t> missing;
  for (size_t i = 0; i < n; ++i) {
    double t_need = horizons[i];
    if (!(t_need > 0.0)) t_need = 1e-3;
    const auto found = cache.index.find(keys[i]);
    if (found != cache.index.end()) {
      const Entry& entry = cache.entries[found->second];
      if (entry.t_max >= t_need &&
          (entry.complete_grid ||
           (query_times != nullptr &&
            entry_has_queries(entry, (*query_times)[i])))) {
        out_indices[i] = found->second;
        continue;
      }
    }
    missing.push_back(i);
  }

  std::stable_sort(missing.begin(), missing.end(),
    [&](size_t a, size_t b) { return horizons[a] < horizons[b]; });
  for (size_t start = 0; start < missing.size();
       start += RLF_BATCH_LANES) {
    const size_t count =
      std::min(RLF_BATCH_LANES, missing.size() - start);
    std::vector<Key> chunk_keys(count);
    std::vector<double> chunk_horizons(count);
    std::vector<std::vector<double>> chunk_queries;
    if (query_times != nullptr) chunk_queries.resize(count);
    for (size_t lane = 0; lane < count; ++lane) {
      const size_t source = missing[start + lane];
      chunk_keys[lane] = keys[source];
      chunk_horizons[lane] = horizons[source];
      if (query_times != nullptr) {
        chunk_queries[lane] = (*query_times)[source];
      }
    }

    std::vector<RLF_Result> batch_results;
    // The lane march builds its own operator and domain and so cannot be
    // paired with a second resolution; extrapolation wins the conflict, since
    // it changes the answer while the lane path only changes the cost.
    const bool use_batch =
      count > 1 && !cache.grid.adaptive && !cache.grid.richardson &&
      cache.grid.explicit_inverse && cache.grid.simd_batch &&
      cache.grid.nx <= 128;
    bool batched = use_batch;
    if (batched) {
      try {
        batch_results = rlf_solve_batch(
          chunk_keys, chunk_horizons, cache.grid,
          query_times == nullptr ? nullptr : &chunk_queries);
      } catch (const std::runtime_error&) {
        // One bad lane fails the whole SIMD march; retry those keys through
        // the scalar path, which can refine its own time grid.
        batched = false;
      }
    }

    for (size_t lane = 0; lane < count; ++lane) {
      const size_t source = missing[start + lane];
      Entry replacement;
      if (batched) {
        rlf_store_result(
          keys[source], horizons[source], batch_results[lane],
          query_times == nullptr, replacement);
      } else {
        const std::vector<double>* queries =
          query_times == nullptr ? nullptr : &(*query_times)[source];
        try {
          rlf_cache_solve(
            keys[source], horizons[source], cache.grid, queries, replacement);
        } catch (const std::exception&) {
          // Unresolvable parameter point: leave it uncached and unindexed and
          // report -1, so only this key's rows are floored.  See cache_get().
          out_indices[source] = -1;
          continue;
        }
      }
      ++cache.solve_count;

      const auto found = cache.index.find(keys[source]);
      if (found != cache.index.end()) {
        cache.entries[found->second] = std::move(replacement);
        out_indices[source] = found->second;
      } else {
        const int idx = static_cast<int>(cache.entries.size());
        cache.entries.push_back(std::move(replacement));
        cache.index.emplace(keys[source], idx);
        out_indices[source] = idx;
      }
    }
  }
}

inline double entry_log_value(const Entry& entry,
                              const std::vector<double>& values,
                              double time, bool pdf) {
  const size_t n = entry.t.size();
  if (n == 0) return RLF_LOG_FLOOR;
  if (time <= entry.t.front()) {
    if (pdf && time < entry.t.front()) return RLF_LOG_FLOOR;
    return values.front();
  }
  if (time >= entry.t.back()) return values.back();
  const auto it = std::lower_bound(entry.t.begin(), entry.t.end(), time);
  const size_t hi = static_cast<size_t>(it - entry.t.begin());
  if (*it == time) return values[hi];
  if (!entry.complete_grid) return RLF_LOG_FLOOR;
  const size_t lo = hi - 1;
  const double w =
    (time - entry.t[lo]) / (entry.t[hi] - entry.t[lo]);
  return values[lo] + w * (values[hi] - values[lo]);
}

inline double entry_log_pdf(const Entry& entry, double time) {
  return entry_log_value(entry, entry.log_pdf, time, true);
}

inline double entry_log_S(const Entry& entry, double time) {
  return entry_log_value(entry, entry.log_S, time, false);
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

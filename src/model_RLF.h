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
// Likelihood evaluation adds dimensionless parameter caching, sparse
// query-time output and a shared solve cache around this numerical core.
//
// Discretisation
// --------------
// Space is second order: the jump part uses fractional centred differences
// (not the first-order shifted Grünwald-Letnikov stencil) and the drift uses
// centred finite-volume faces (not donor-cell upwinding, whose v*h/2 numerical
// diffusion and v*p_n boundary flux dominated every other error source).
//
// Time is not discretised at all.  The semi-discrete problem dp/dt = L p is
// autonomous and only three linear functionals of p are ever needed (the exit
// flux, the surviving mass and the censored lower-edge rate), so the solution
// is obtained as exp(t L) p0 by shift-invert Arnoldi: a Krylov basis is built
// on (I - gamma L)^-1, the projected generator is eigendecomposed, and every
// output quantity comes out as a closed-form sum of m complex exponentials,
// evaluable at any t in O(m), with m chosen per solve by convergence.  This replaced a TR-BDF2 march;
// see the note above rlf_build_modes for why, and for the measurements.
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

// R_ext/Lapack.h declares the entire LAPACK surface and clashes noisily with
// RcppArmadillo's declarations when this header is used by particle_ll.cpp.
// The propagator only needs these five routines.
extern "C" {
void F77_NAME(dgetrf)(const int*, const int*, double*, const int*, int*, int*);
void F77_NAME(dgetrs)(const char*, const int*, const int*, const double*,
                      const int*, const int*, double*, const int*, int*,
                      size_t);
void F77_NAME(dgemv)(const char*, const int*, const int*, const double*,
                     const double*, const int*, const double*, const int*,
                     const double*, double*, const int*, size_t);
void F77_NAME(dhseqr)(const char*, const char*, const int*, const int*,
                      const int*, double*, const int*, double*, double*,
                      double*, const int*, double*, const int*, int*,
                      size_t, size_t);
void F77_NAME(dgecon)(const char*, const int*, const double*, const int*,
                      const double*, double*, double*, int*, int*, size_t);
void F77_NAME(dtrevc)(const char*, const char*, int*, const int*,
                      const double*, const int*, double*, const int*, double*,
                      const int*, const int*, int*, double*, int*,
                      size_t, size_t);
}

namespace rlf {

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// The shifted operator is assembled row-major, which BLAS reads as the
// column-major transpose.  Factoring those bytes therefore gives the LU of
// A^T, and solving with trans = 'T' recovers A x = b.  Only the factorisation
// is formed: the explicit inverse costs 4/3 n^3 more than DGETRF alone and
// only pays for itself past roughly n solves, while shift-invert Arnoldi needs
// only around 24 of them (measured at n = 199: DGETRF 0.215 ms against
// DGETRF+DGETRI 0.870 ms, and 24 DGETRS 0.187 ms against 24 DGEMV 0.084 ms).
struct RLF_ShiftedLU {
  int n = 0;
  std::vector<double> lu;
  std::vector<int> pivot;

  bool factor(int n_, std::vector<double>& A) {
    n = n_;
    lu.swap(A);
    pivot.resize(n);
    int info = 0;
    F77_CALL(dgetrf)(&n, &n, lu.data(), &n, pivot.data(), &info);
    return info == 0;
  }

  // In-place solve of A x = b.
  bool solve(double* b) const {
    const char transpose = 'T';
    const int one = 1;
    int info = 0;
    F77_CALL(dgetrs)(
      &transpose, &n, &one, lu.data(), &n, pivot.data(), b, &n, &info, 1);
    return info == 0;
  }
};


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
  // Deepest negative excursion of the density, as a fraction of its peak.
  double min_density = 0.0;
  double domain_pdf_error = std::numeric_limits<double>::quiet_NaN();
  double domain_cdf_error = std::numeric_limits<double>::quiet_NaN();
  double spatial_pdf_error = std::numeric_limits<double>::quiet_NaN();
  double spatial_cdf_error = std::numeric_limits<double>::quiet_NaN();
  double x_lo = 0.0;
  double dx = 0.0;
  int nx_used = 0;
  int krylov_dim = 0;
  // 1-norm condition estimate of the reduced eigenbasis at the accepted
  // dimension, and the largest shift-and-invert residual over the probe times.
  // Diagnostics: see rlf_reduce_modes.
  double eigen_condition = 1.0;
  double krylov_residual = 0.0;
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
  // The generator is Toeplitz apart from its diagonal: every state-dependent
  // term (the censoring rate added back at the lower edge and the sealed
  // advective face at x_lo) lives there.  Only this compact form is kept; the
  // dense n x n generator is never assembled.
  std::vector<double> toeplitz_col;
  std::vector<double> toeplitz_row;
  std::vector<double> toeplitz_diag;
  std::vector<double> upper_kill;
  std::vector<double> lower_censor;
  // Graded form.  When `dense` is non-empty the mesh is nonuniform, the
  // generator has no Toeplitz structure left to exploit, and `width` carries
  // the per-cell quadrature weights that a single h stood in for.
  std::vector<double> dense;
  std::vector<double> width;
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

inline RLF_Operator build_rlf_operator(
    const RLF_Model& m, int n, double h) {
  RLF_Operator op;
  op.n = n;
  op.toeplitz_col.assign(n, 0.0);
  op.toeplitz_row.assign(n, 0.0);
  op.toeplitz_diag.assign(n, 0.0);
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

  for (int d = 0; d < n; ++d) {
    op.toeplitz_col[d] = q[d];
    op.toeplitz_row[d] = q[d];
  }

  // Column j is a source state.  Omitted targets above the barrier are killed
  // and contribute to the FPT flux.  Omitted targets below x_lo are censored:
  // adding that rate back to the diagonal suppresses those artificial exits.
  for (int j = 0; j < n; ++j) {
    const double lower = tail[j + 1];
    const double upper = tail[n - j];
    op.lower_censor[j] = lower;
    op.upper_kill[j] = upper;
    op.toeplitz_diag[j] = q[0] + lower;
  }
  // Positive-drift transport, by exponential fitting of the nearest-neighbour
  // pair (Il'in / Scharfetter-Gummel).  Only |i-j| = 1 changes, so the
  // operator stays Toeplitz-plus-diagonal and the drift is uniform in j.
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
      // The fitted pair is translation invariant, so the compact form only
      // needs the |i-j| = 1 entries once.
      if (j == 0) op.toeplitz_col[1] += delta_up;
    } else {
      // The cell above the last one is past the barrier: killed, not moved.
      op.upper_kill[j] += delta_up;
    }
    if (j > 0) {
      if (j == 1) op.toeplitz_row[1] += delta_down;
      outflow += delta_down;
    }
    // At j = 0 the downward face is the truncated lower edge, whose rate is
    // censored back out anyway, so only the upward change leaves the diagonal.
    op.toeplitz_diag[j] -= outflow;
  }

  // The jump part's column sum is q summed over the reachable off-diagonal
  // band plus the diagonal, which prefix sums of q give in O(1) per column,
  // and the drift terms cancel in every column -- only the boundary kill at
  // the last cell survives.  So the conservation diagnostic costs O(n)
  // without L ever being formed.
  std::vector<double> q_prefix(n, 0.0);
  for (int d = 1; d < n; ++d) q_prefix[d] = q_prefix[d - 1] + q[d];
  for (int j = 0; j < n; ++j) {
    const double exit = -op.toeplitz_diag[j];
    op.max_exit_rate = std::max(op.max_exit_rate, exit);

    double col_sum = op.upper_kill[j] + op.toeplitz_diag[j] +
      q_prefix[j] + q_prefix[n - 1 - j];
    if (j + 1 < n) col_sum += delta_up;
    if (j > 0) col_sum += delta_down;
    op.conservation_error =
      std::max(op.conservation_error, std::abs(col_sum));
  }
  return op;
}


// ---------------------------------------------------------------------------
// Graded mesh.
//
// The uniform grid ties two unrelated requirements to one number.  The mesh is
// h = (b0 + L)/nx, so the depth L that keeps the lower tail on the domain is
// bought out of the resolution at the absorbing edge, and rlf_lower_extent has
// to settle that trade.  It settles it on absolute CDF error, which is the
// wrong currency for the survivor: log S is relative, and at long horizons S
// is e^-10 or smaller, so a truncation that is invisible in the CDF is worth
// ten nats in the likelihood.  Measured against a wide converged reference the
// shipped default is out by 0.68 nats in log S from truncation alone at a mesh
// where halving h is worth 0.014 -- the depth, not the mesh, is what limits
// the solver.
//
// Grading breaks the tie.  The stable tail only has to be represented, not
// resolved: it is smooth and slowly varying, so cells that grow geometrically
// downward cover a hundred units in twenty-five cells, while the core keeps
// the fine mesh the boundary layer needs.
//
// Faces run from x_lo up to b0.  The core is uniform over the top
// `core_width`; below it widths grow by `ratio` per cell until x_lo is
// reached.  The last coarse cell is stretched rather than split so that the
// bottom face lands exactly on x_lo.
struct RLF_Mesh {
  std::vector<double> face;   // n + 1 ascending faces, face[0] = x_lo
  std::vector<double> node;   // n cell centres
  std::vector<double> width;  // n cell widths
  int n = 0;
  // Cells [uniform_from, n) all have width `uniform_h`.  Face separations
  // inside that block depend only on the index difference, which is what keeps
  // the kernel table's largest block Toeplitz.
  int uniform_from = 0;
  double uniform_h = 0.0;
};

inline RLF_Mesh rlf_build_mesh(double x_lo, double b0, double core_width,
                               int n_core, double ratio,
                               int n_edge = 0, double edge_ratio = 1.0) {
  RLF_Mesh mesh;
  const double span = b0 - x_lo;
  const double core = std::min(std::max(core_width, 0.0), span);
  const int nc = std::max(1, n_core);
  const double h = core > 0.0 ? core / nc : span / nc;

  std::vector<double> widths;
  widths.reserve(nc + 64);
  // Coarse cells first, from x_lo upward, so the geometric run is built from
  // its wide end and the stretch lands on the cell nearest x_lo.
  double remaining = span - core;
  if (remaining > 1e-12 * span) {
    std::vector<double> coarse;
    double w = h * ratio;
    double covered = 0.0;
    for (int guard = 0; guard < 4096 && covered + w < remaining; ++guard) {
      coarse.push_back(w);
      covered += w;
      w *= ratio;
    }
    // Absorb what is left into one final cell rather than leaving a sliver.
    const double last = remaining - covered;
    if (last > 0.0) {
      if (!coarse.empty() && last < 0.5 * coarse.back()) {
        coarse.back() += last;
      } else {
        coarse.push_back(last);
      }
    }
    for (auto it = coarse.rbegin(); it != coarse.rend(); ++it) {
      widths.push_back(*it);
    }
  }
  for (int i = 0; i < nc; ++i) widths.push_back(h);

  // Geometric refinement into the absorbing edge.  The killed density vanishes
  // like (b0 - x)^(alpha/2), an algebraic singularity that no uniform mesh
  // resolves at any affordable width, so the top `n_edge` core cells are
  // replaced by the same number of cells covering the same span but shrinking
  // toward b0.  Geometric refinement against an algebraic singularity buys
  // accuracy exponentially in the cell count, which is why a handful of cells
  // is enough.
  if (n_edge > 1 && edge_ratio > 1.0 &&
      static_cast<int>(widths.size()) >= n_edge) {
    double span = 0.0;
    for (int k = 0; k < n_edge; ++k) span += widths[widths.size() - 1 - k];
    widths.resize(widths.size() - n_edge);
    // Widths w, w/r, w/r^2, ... summing to span, appended largest first.
    double total = 0.0, term = 1.0;
    for (int k = 0; k < n_edge; ++k) { total += term; term /= edge_ratio; }
    double w = span / total;
    for (int k = 0; k < n_edge; ++k) { widths.push_back(w); w /= edge_ratio; }
  }

  mesh.n = static_cast<int>(widths.size());
  mesh.uniform_from = mesh.n;
  mesh.uniform_h = h;
  while (mesh.uniform_from > 0 &&
         std::abs(widths[mesh.uniform_from - 1] - h) <= 1e-12 * h) {
    --mesh.uniform_from;
  }
  mesh.width = widths;
  mesh.face.assign(mesh.n + 1, 0.0);
  mesh.face[0] = x_lo;
  for (int i = 0; i < mesh.n; ++i) mesh.face[i + 1] = mesh.face[i] + widths[i];
  mesh.face[mesh.n] = b0;  // exact, against accumulated rounding
  mesh.node.assign(mesh.n, 0.0);
  for (int i = 0; i < mesh.n; ++i) {
    mesh.node[i] = 0.5 * (mesh.face[i] + mesh.face[i + 1]);
  }
  return mesh;
}

// Bernoulli function t / (e^t - 1), the exponential-fitting weight.
inline double rlf_bernoulli(double t) {
  if (t > 700.0) return 0.0;
  if (t < -700.0) return -t;
  if (std::abs(t) < 1e-8) return 1.0 - 0.5 * t;
  return t / std::expm1(t);
}

// Graded generator.
//
// The Levy measure is nu(z) = c_alpha |z|^(-1-alpha).  Writing Phi for its
// second antiderivative, Phi(u) = |u|^(1-alpha) / (alpha (alpha-1)), the total
// jump rate between two disjoint intervals is a second difference of Phi at
// the four face separations, so one table of Phi over face pairs supplies
// every entry with three additions.  That table is the only transcendental
// work in the assembly and it is symmetric, so it costs n^2/2 powers.
//
// Cells one apart cannot be treated this way -- the double integral diverges
// for alpha > 1, which is the whole supported range, and is exactly why a
// naive jump-chain discretisation is unavailable on a nonuniform mesh.  What
// diverges is the short-jump part, and short jumps are a diffusion: the net
// rate at which a linear density is carried across the shared face by jumps
// between the two cells is
//
//   D = c_alpha \int_{C_i} \int_{C_j} |x - y|^(-alpha),
//
// in closed form again, and this is used as the diffusivity of a
// Scharfetter-Gummel face flux.  There is no cutoff parameter anywhere: the
// split is by cell adjacency, and each side is integrated exactly over its own
// region.  At alpha = 2 the far kernel vanishes identically (c_alpha carries
// sin(pi alpha / 2)) and D tends to sigma^2 / 2, so the diffusion limit is
// reached continuously rather than by a special case.
//
// Both parts are exactly conservative -- the jump matrix by the symmetry of
// Omega, the face fluxes by telescoping -- and both have non-negative
// off-diagonals, so the generator stays an M-matrix and exp(tL) stays
// positive.
inline RLF_Operator build_rlf_operator_graded(const RLF_Model& m,
                                              const RLF_Mesh& mesh,
                                              std::vector<double>* g_scratch = nullptr) {
  const int n = mesh.n;
  RLF_Operator op;
  op.n = n;
  op.width = mesh.width;
  op.dense.assign(static_cast<size_t>(n) * n, 0.0);
  op.upper_kill.assign(n, 0.0);
  op.lower_censor.assign(n, 0.0);

  const double alpha = m.alpha;
  const double eps = 2.0 - alpha;
  const double c_alpha = 0.5 * std::pow(m.sigma, alpha) *
    std::tgamma(1.0 + alpha) * std::sin(M_PI * alpha * 0.5) / M_PI;
  // c_alpha / (2 - alpha), taken through the removable zero at alpha = 2 so
  // that the near-field diffusivity reaches sigma^2 / 2 rather than 0/0.
  const double sinc = std::abs(eps) < 1e-9
    ? 1.0 - (M_PI * eps * 0.5) * (M_PI * eps * 0.5) / 6.0
    : std::sin(M_PI * eps * 0.5) / (M_PI * eps * 0.5);
  const double c_over_eps =
    0.5 * std::pow(m.sigma, alpha) * std::tgamma(1.0 + alpha) * 0.5 * sinc;

  const double phi_scale = 1.0 / (alpha * (alpha - 1.0));
  auto phi = [&](double u) {
    return phi_scale * std::pow(std::abs(u), 1.0 - alpha);
  };
  // Third antiderivative, odd, carrying its (2 - alpha) pole outside so that
  // the pairing with c_alpha stays finite in the diffusion limit.  It is the
  // second one multiplied by the separation, so it costs no extra power.
  auto chi = [&](double u) { return phi(u) * u; };

  // Phi and Chi over every face pair, built once.  These tables are the only
  // transcendental work in the assembly; everything downstream is additions.
  //
  // The powers are the expensive part, and most of them are avoidable: the
  // core of the mesh is uniform, so inside that block a face separation
  // depends only on the index difference and one vector of n_core values fills
  // an n_core^2 block.  Only the graded tail, which is short by construction,
  // needs a power per pair.
  const int nf = n + 1;
  const int u0 = mesh.uniform_from;
  std::vector<double> band(nf, 0.0);
  for (int k = 1; k < nf - u0; ++k) band[k] = phi(k * mesh.uniform_h);
  std::vector<double> g_local;
  std::vector<double>& G = g_scratch != nullptr ? *g_scratch : g_local;
  G.assign(static_cast<size_t>(nf) * nf, 0.0);
  for (int p = 0; p < nf; ++p) {
    for (int q = p + 1; q < nf; ++q) {
      const double d = mesh.face[q] - mesh.face[p];
      const double gv = (p >= u0) ? band[q - p] : phi(d);
      G[static_cast<size_t>(p) * nf + q] = gv;
      G[static_cast<size_t>(q) * nf + p] = gv;
    }
  }
  auto Gv = [&](int p, int q) { return G[static_cast<size_t>(p) * nf + q]; };
  auto Xv = [&](int p, int q) {
    return Gv(p, q) * (mesh.face[p] - mesh.face[q]);
  };

  // Piecewise-linear reconstruction inside each cell.
  //
  // A cell-constant density makes every kernel integral exact for the
  // representation but leaves the scheme first order, and measurably so: the
  // error runs as (2 - alpha) h, vanishing in the diffusion limit where the
  // far kernel vanishes with c_alpha and growing as the tail gets heavier.
  // What it misses is the first moment of the density inside the source cell,
  // which the kernel weights asymmetrically because it is falling steeply
  // across the cell.  Carrying a slope fixes the order, and the moment
  // integral it needs is the same closed form one antiderivative further on.
  //
  // s_j = sl_lo[j] p_{j-1} + sl_mid[j] p_j + sl_hi[j] p_{j+1}, one-sided at
  // the two ends.
  std::vector<double> sl_lo(n, 0.0), sl_mid(n, 0.0), sl_hi(n, 0.0);
  for (int j = 0; j < n; ++j) {
    if (n < 2) break;
    if (j == 0) {
      const double g = mesh.node[1] - mesh.node[0];
      sl_mid[0] = -1.0 / g; sl_hi[0] = 1.0 / g;
    } else if (j == n - 1) {
      const double g = mesh.node[n - 1] - mesh.node[n - 2];
      sl_lo[j] = -1.0 / g; sl_mid[j] = 1.0 / g;
    } else {
      const double g = mesh.node[j + 1] - mesh.node[j - 1];
      sl_lo[j] = -1.0 / g; sl_hi[j] = 1.0 / g;
    }
  }

  double* __restrict L = op.dense.data();
  auto at = [&](int i, int j) -> double& {
    return L[static_cast<size_t>(i) * n + j];
  };

  // --- far field: every pair of cells at least one cell apart
  //
  // Omega is the total jump rate between the two cells, a second difference of
  // Phi at the four face separations.  Mu is the same integral weighted by the
  // signed offset of the source point from its cell centre, which is what the
  // slope multiplies; it is one antiderivative further on and reuses the same
  // two tables.  Both are exact, so the only approximation left is that the
  // density is linear across a cell.
  auto scatter = [&](int i, int j, double coef) {
    if (coef == 0.0) return;
    if (sl_lo[j] != 0.0) at(i, j - 1) += coef * sl_lo[j];
    if (sl_mid[j] != 0.0) at(i, j) += coef * sl_mid[j];
    if (sl_hi[j] != 0.0) at(i, j + 1) += coef * sl_hi[j];
  };
  // The loss term is the same moment with the two cells' roles swapped, so it
  // needs no second evaluation: summing mu over targets for a fixed source is
  // exactly the quantity each row subtracts from itself.
  std::vector<double> source_moment(n, 0.0);
  for (int i = 0; i < n; ++i) {
    const double inv_w = 1.0 / mesh.width[i];
    double outflow = 0.0;
    for (int j = 0; j < n; ++j) {
      if (std::abs(i - j) < 2) continue;
      const double g00 = Gv(i, j), g01 = Gv(i, j + 1);
      const double g10 = Gv(i + 1, j), g11 = Gv(i + 1, j + 1);
      double omega = c_alpha * ((g10 + g01) - (g11 + g00));
      if (!(omega > 0.0)) omega = 0.0;  // rounding in the far tail only
      at(i, j) += omega * inv_w;
      outflow += omega * inv_w;

      const double mu =
        -0.5 * mesh.width[j] * c_alpha * ((g10 - g00) + (g11 - g01)) +
        c_over_eps *
          ((Xv(i + 1, j) + Xv(i, j + 1)) - (Xv(i + 1, j + 1) + Xv(i, j)));
      scatter(i, j, mu * inv_w);
      source_moment[j] += mu;
    }
    at(i, i) -= outflow;
  }
  for (int i = 0; i < n; ++i) {
    scatter(i, i, -source_moment[i] / mesh.width[i]);
  }

  // --- adjacent faces: near-field diffusivity plus drift, by exponential
  //     fitting.  Face k sits between cells k-1 and k.
  auto near_diffusivity = [&](double wa, double wb) {
    // c_alpha \int\int |x-y|^{-alpha} over the two adjacent cells, written so
    // that the (2-alpha) zero cancels analytically.
    const double bracket =
      std::pow(wa + wb, eps) - std::pow(wa, eps) - std::pow(wb, eps);
    return c_over_eps * bracket / (1.0 - alpha);
  };

  for (int k = 1; k < n; ++k) {
    const int i = k - 1, j = k;
    const double D = near_diffusivity(mesh.width[i], mesh.width[j]);
    const double g = mesh.node[j] - mesh.node[i];
    const double T = D / g;
    double up, down;
    if (T > 0.0) {
      const double pe = m.v * g / D;
      up = T * rlf_bernoulli(-pe);    // i -> j
      down = T * rlf_bernoulli(pe);   // j -> i
    } else {
      // No near-field diffusion to fit against; pure upwind transport.
      up = m.v / g;
      down = 0.0;
    }
    at(i, i) -= up / mesh.width[i];
    at(i, j) += down / mesh.width[i];
    at(j, j) -= down / mesh.width[j];
    at(j, i) += up / mesh.width[j];
  }

  // --- absorbing edge at b0
  //
  // The cell below the barrier loses mass two ways: across the face itself,
  // and by jumping clear over the first ghost cell.  Giving the ghost a real
  // width is what keeps the second term finite -- integrated from b0 it would
  // diverge, because a source arbitrarily close to the barrier is killed by
  // arbitrarily short jumps.  Every lower cell is at least two cells from the
  // killing region, so its jump kill is the plain closed form.
  {
    const int top = n - 1;
    const double w_top = mesh.width[top];
    const double D = near_diffusivity(w_top, w_top);
    // The barrier is the face, not the ghost cell's centre.  The linear
    // reconstruction that vanishes at b0 has gradient p_top / (w_top / 2), so
    // the conductance distance is half a cell -- using the centre-to-centre
    // distance instead puts the absorbing boundary at b0 + w_top / 2 and costs
    // a clean order of accuracy at every alpha, the diffusion limit included.
    const double g = 0.5 * w_top;
    double out_rate;
    if (D > 0.0) {
      const double pe = m.v * g / D;
      out_rate = (D / g) * rlf_bernoulli(-pe);
    } else {
      out_rate = m.v / g;
    }
    const double face_kill = out_rate / w_top;
    at(top, top) -= face_kill;
    op.upper_kill[top] += face_kill;

    const double b = mesh.face[n];
    const double ghost = b + w_top;
    for (int j = 0; j < n; ++j) {
      const double lo = mesh.face[j], hi = mesh.face[j + 1];
      const double target = (j == top) ? ghost : b;
      const double inv_w = 1.0 / mesh.width[j];
      double omega = c_alpha * (phi(target - hi) - phi(target - lo));
      if (!(omega > 0.0)) omega = 0.0;
      at(j, j) -= omega * inv_w;
      op.upper_kill[j] += omega * inv_w;
      // First moment over the source cell, the killing counterpart of mu.
      const double mu = 0.5 * mesh.width[j] * c_alpha *
          (phi(target - lo) + phi(target - hi)) +
        c_over_eps * (chi(target - hi) - chi(target - lo));
      if (sl_lo[j] != 0.0) {
        at(j, j - 1) -= mu * sl_lo[j] * inv_w;
        op.upper_kill[j - 1] += mu * sl_lo[j] / mesh.width[j - 1];
      }
      if (sl_mid[j] != 0.0) {
        at(j, j) -= mu * sl_mid[j] * inv_w;
        op.upper_kill[j] += mu * sl_mid[j] * inv_w;
      }
      if (sl_hi[j] != 0.0) {
        at(j, j + 1) -= mu * sl_hi[j] * inv_w;
        op.upper_kill[j + 1] += mu * sl_hi[j] / mesh.width[j + 1];
      }
    }
  }

  // --- censored edge at x_lo
  //
  // Mass that would leave the bottom of the domain is not removed: those paths
  // are misrepresented by the truncation either way, and suppressing the exit
  // is the closure that keeps the survivor from decaying for a reason the
  // model does not have.  So nothing is subtracted here and the rate is only
  // recorded, as the pressure diagnostic that says whether the domain is deep
  // enough.
  {
    const double a = mesh.face[0];
    const double ghost = a - mesh.width[0];
    for (int j = 0; j < n; ++j) {
      const double lo = mesh.face[j], hi = mesh.face[j + 1];
      const double target = (j == 0) ? ghost : a;
      double omega = c_alpha * (phi(lo - target) - phi(hi - target));
      if (!(omega > 0.0)) omega = 0.0;
      op.lower_censor[j] = omega / mesh.width[j];
    }
  }

  for (int j = 0; j < n; ++j) {
    op.max_exit_rate = std::max(op.max_exit_rate, -at(j, j));
  }
  // Conservation: mass leaves only through the barrier, so the width-weighted
  // column sums must equal minus the kill rate.
  for (int j = 0; j < n; ++j) {
    double col = 0.0;
    for (int i = 0; i < n; ++i) col += mesh.width[i] * at(i, j);
    op.conservation_error = std::max(
      op.conservation_error,
      std::abs(col / mesh.width[j] + op.upper_kill[j]));
  }
  return op;
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

// Mesh for a parameter point, sized by one budget: the number of core cells.
//
// Grading decouples the two numbers rlf_lower_extent had to trade off, so both
// can now be set on their own terms.  The depth is generous -- sixty stable
// widths, where the uniform grid could afford six -- because reaching it costs
// a logarithmic number of cells rather than a proportional one, and because
// log S at a long horizon is a relative quantity that a shallow domain gets
// wrong by nats.  The core covers the barrier and the few widths below it that
// carry the mass; below that the cells grow by 45% each, which measurement puts
// at the flat part of the accuracy-cost curve (1.15 costs cells for nothing,
// 1.7 starts to bite).
inline double rlf_mesh_ratio = 1.45;
inline double rlf_mesh_extent_widths = 60.0;
inline double rlf_mesh_core_widths = 2.0;

inline RLF_Mesh rlf_auto_mesh(const RLF_Model& m, double t_max, int n_core) {
  const double spread =
    std::max(m.b0, m.sigma * std::pow(0.5 * t_max, 1.0 / m.alpha));
  const double extent = rlf_mesh_extent_widths * spread;
  const double core =
    std::min(extent, m.b0 + rlf_mesh_core_widths * spread);
  return rlf_build_mesh(-extent, m.b0, core, std::max(n_core, 20),
                        rlf_mesh_ratio);
}

// ---------------------------------------------------------------------------
// Shift-invert Krylov propagator.
//
// dp/dt = L p is autonomous, L is a fixed matrix over the whole horizon, and
// the likelihood only ever reads three linear functionals of p.  A time march
// is therefore doing far more work than the problem needs: it propagates all n
// state components through hundreds of steps to recover three scalars per
// query time, and it pays an O(n^3) setup per distinct step size to make each
// of those steps cheap.
//
// Arnoldi on B = (I - gamma L)^-1 instead gives B V_m = V_m H_m + h v e_m^T,
// from which the projected generator is L_m = (I - H_m^-1)/gamma.  Its
// eigendecomposition turns every functional into a closed-form sum of m
// complex exponentials, so the whole curve is available at once and any query
// time costs O(m).  The shift is what makes this work at this size: it maps
// the eigenvalues nearest zero -- exactly the modes that survive to the times
// the likelihood asks about -- to the outside of B's spectrum, so convergence
// is governed by the resolved dynamics rather than by ||L|| t.  A polynomial
// Krylov space is not competitive here for the same reason a time march is
// not: it needs m ~ ||L|| t.
//
// Measured against a converged (nt = 20000) TR-BDF2 march on the identical
// grid, as max relative PDF error over the central mass:
//
//                          m = 16    m = 24    m = 32
//   alpha 1.5, nx 160      1.4e-3    2.2e-4    6.6e-6
//   alpha 1.1, nx 160      5.7e-3    1.1e-3    3.6e-6
//   alpha 1.8, nx 160      2.4e-3    9.1e-5    2.0e-5
//   alpha 2.0, nx 160      1.3e-3    2.1e-4    2.0e-4     (||L|| t = 9021)
//   alpha 1.3, nx 320      2.7e-2    4.9e-3    1.1e-4
//
// Convergence is flat in stiffness and in the horizon, which a march is not:
// at m = 24 the error is two orders of magnitude below the ~5e-3 spatial
// discretisation floor the extrapolated pair leaves behind.  End to end
// against the march it replaced this is 2.4x to 4.9x faster, with the largest
// gains at long horizons and refined grids, and it is also more accurate,
// because the O(dt^2) time error is gone rather than merely small.
// ---------------------------------------------------------------------------

// Arnoldi dimension: where the search starts, how it grows, where it stops,
// and the log-scale movement below which it is called converged.  Ordinary
// points finish at the minimum; the maximum is only reached in the heavy-tail,
// fast-drift, fine-grid corner.  See rlf_build_modes for the calibration.
inline int RLF_KRYLOV_MIN = 20;
inline int RLF_KRYLOV_STEP = 6;
inline int RLF_KRYLOV_MAX = 64;
inline double RLF_KRYLOV_TOL = 1e-5;
// Density, as a fraction of the peak, below which the flux convergence test
// stops being relative.  Four decades covers the range that carries data.
inline double RLF_KRYLOV_FLOOR = 1e-4;
// Shift-and-invert residual, and condition estimate of the reduced
// eigenbasis, past which a solve is not allowed to skip the refinement checks.
// Both are diagnostics rather than stopping rules; see rlf_reduce_modes for
// what each one measures and rlf_build_modes for why neither replaces the
// observable convergence test.
inline double RLF_KRYLOV_RESID = 1e-4;
inline double RLF_KRYLOV_CONDITION = 1e10;
// Shift, as a fraction of the horizon.  The Krylov space resolves eigenvalues
// of size ~1/gamma, so this trades the fast modes that carry the leading edge
// against the slow ones that carry the tail.  It is calibrated rather than
// argued: sweeping it against a converged reference over alpha, v, b, the
// horizon and the grid puts the optimum in a broad basin around 0.02, and the
// conventional t_max/10 sits far enough up the fast-mode side of that basin to
// cost two orders of magnitude of accuracy at no saving.
inline double RLF_KRYLOV_SHIFT = 0.02;
// Dimension the last solve settled on.  Diagnostic only; the Richardson pair
// leaves the second member's value here.
inline int rlf_last_krylov_dim = 0;

// A curve reconstructed from the projected eigendecomposition.
//
//   f(t) = sum_j exp(re_j t) [a_j cos(im_j t) - b_j sin(im_j t)]
//
// DGEEV returns conjugate pairs adjacently with the eigenvector's real part in
// the first column and its imaginary part in the second.  The first member of
// a pair carries im > 0 and holds the whole (already real) contribution; the
// second carries im < 0 and is skipped.  Real modes have im == 0 and b == 0.
struct RLF_Curve {
  std::vector<double> a, b;

  double at(const std::vector<double>& re, const std::vector<double>& im,
            double t) const {
    double out = 0.0;
    const size_t m = a.size();
    for (size_t j = 0; j < m; ++j) {
      if (im[j] < 0.0) continue;
      const double growth = std::exp(re[j] * t);
      if (im[j] == 0.0) {
        out += a[j] * growth;
      } else {
        const double phase = im[j] * t;
        out += growth *
          (a[j] * std::cos(phase) - b[j] * std::sin(phase));
      }
    }
    return out;
  }
};

struct RLF_Modes {
  int m = 0;
  std::vector<double> re, im;
  RLF_Curve flux;      // probability flux through the absorbing boundary
  RLF_Curve surv;      // surviving mass
  RLF_Curve censor;    // rate at which the censored lower edge is pressed
  // Modal form of  beta (h_{m+1,m}/gamma) e_m^T H_m^{-1} exp(t L_m) e_1, whose
  // magnitude is the norm of the shift-and-invert residual.
  RLF_Curve defect;
  // 1-norm condition estimate of the (column-normalised) reduced eigenvector
  // matrix.  The projected generator is nonnormal, so the modal expansion can
  // in principle be a cancelling sum of huge terms; this is what would say so.
  double eigen_condition = 1.0;
  // beta * V and the modal coefficients of the reduced state, kept so the
  // density itself can be reconstructed for the positivity diagnostic.
  std::vector<double> basis;       // n x m, column major
  std::vector<double> state_a, state_b;   // m x m, column major
  int n = 0;

  // Every output time needs the flux, the survivor and the integrated flux,
  // and each of those is the same sum over the same modes with the same
  // exp/cos/sin envelope.  Building the envelope once per time and reducing
  // each curve against it turns four transcendental sweeps into one, which is
  // what the many-distinct-response-time case is made of.
  void envelope_at(double t, std::vector<double>& weight) const {
    weight.assign(static_cast<size_t>(2 * m), 0.0);
    for (int j = 0; j < m; ++j) {
      if (im[j] < 0.0) continue;
      const double growth = std::exp(re[j] * t);
      if (im[j] == 0.0) {
        weight[2 * j] = growth;
      } else {
        const double phase = im[j] * t;
        weight[2 * j] = growth * std::cos(phase);
        weight[2 * j + 1] = -growth * std::sin(phase);
      }
    }
  }

  // Reduce one curve against an envelope built by envelope_at.  Matches
  // RLF_Curve::at, which forms a * cos - b * sin.
  static double reduce(const RLF_Curve& curve,
                       const std::vector<double>& weight) {
    double out = 0.0;
    const size_t count = curve.a.size();
    for (size_t j = 0; j < count; ++j) {
      out += curve.a[j] * weight[2 * j] + curve.b[j] * weight[2 * j + 1];
    }
    return out;
  }

  double flux_at(double t) const { return flux.at(re, im, t); }
  double surv_at(double t) const { return surv.at(re, im, t); }
  double censor_at(double t) const { return censor.at(re, im, t); }
  // 2-norm of the shift-and-invert residual  p_m'(t) - L p_m(t)  at time t.
  // See rlf_reduce_modes for how the curve is built.
  double residual_at(double t) const {
    return std::abs(defect.at(re, im, t));
  }

  // Exact antiderivative of a modal curve, as a curve plus the constant that
  // makes it vanish at t = 0:  int_0^t sum_j c_j e^{lam_j s} ds
  //   = sum_j (c_j / lam_j) e^{lam_j t} - sum_j c_j / lam_j.
  RLF_Curve integrate(const RLF_Curve& curve, double& constant) const {
    RLF_Curve out;
    out.a.assign(m, 0.0);
    out.b.assign(m, 0.0);
    constant = 0.0;
    for (int j = 0; j < m; ++j) {
      if (im[j] < 0.0) continue;
      if (im[j] == 0.0) {
        if (!(std::abs(re[j]) > 0.0)) continue;
        out.a[j] = curve.a[j] / re[j];
      } else {
        // (a + i b) / (re + i im)
        const double d = re[j] * re[j] + im[j] * im[j];
        if (!(d > 0.0)) continue;
        out.a[j] = (curve.a[j] * re[j] + curve.b[j] * im[j]) / d;
        out.b[j] = (curve.b[j] * re[j] - curve.a[j] * im[j]) / d;
      }
      constant -= out.a[j];
    }
    return out;
  }

  // Smallest and largest cell of the reconstructed density at time t.  exp(t L)
  // is entrywise non-negative because L is an M-matrix generator, so anything
  // materially below zero here is Krylov truncation rather than a genuine
  // property of the discretisation -- which is exactly what makes it a useful
  // check.  The peak is returned alongside so the undershoot can be reported
  // as a fraction of it: the density scales like 1/h, so an absolute floor on
  // the trough would tighten itself every time the grid is refined.
  void density_extremes_at(double t, double& lowest, double& highest) const {
    if (basis.empty()) return;
    std::vector<double> y(m, 0.0);
    for (int j = 0; j < m; ++j) {
      if (im[j] < 0.0) continue;
      const double growth = std::exp(re[j] * t);
      const double* first = &state_a[static_cast<size_t>(j) * m];
      if (im[j] == 0.0) {
        for (int i = 0; i < m; ++i) y[i] += growth * first[i];
      } else {
        const double phase = im[j] * t;
        const double c = std::cos(phase), s = std::sin(phase);
        const double* second = &state_b[static_cast<size_t>(j) * m];
        for (int i = 0; i < m; ++i) {
          y[i] += growth * (first[i] * c - second[i] * s);
        }
      }
    }
    for (int i = 0; i < n; ++i) {
      double value = 0.0;
      for (int j = 0; j < m; ++j) {
        value += basis[static_cast<size_t>(j) * n + i] * y[j];
      }
      lowest = std::min(lowest, value);
      highest = std::max(highest, value);
    }
  }
};

// Reusable scratch for one solve.  The Richardson pair changes n between its
// two members, so these only avoid reallocation within a member, but the
// buffers are small next to the factorisation they feed.
struct RLF_KrylovWork {
  std::vector<double> shifted, basis, hessenberg, vector, projection;
  std::vector<double> reduced, eigvec, coefficient, lapack;
  // Reusable face-pair kernel table for graded operator assembly.
  std::vector<double> operator_G;
  // H_m is consumed by the Schur factorisation, but the residual estimate
  // needs it afterwards, so it is kept.
  std::vector<double> hess_copy, residual_row;
  std::vector<int> pivot, iwork;
};

// Assemble I - gamma L directly from the compact operator.  L itself is never
// formed: it is Toeplitz apart from its diagonal, and the only consumer is
// this matrix.
inline void rlf_shifted_operator(const RLF_Operator& op, double gamma,
                                 std::vector<double>& out) {
  const int n = op.n;
  out.resize(static_cast<size_t>(n) * n);
  if (!op.dense.empty()) {
    const double* __restrict src = op.dense.data();
    for (int i = 0; i < n; ++i) {
      double* __restrict row = &out[static_cast<size_t>(i) * n];
      const double* __restrict from = &src[static_cast<size_t>(i) * n];
      for (int j = 0; j < n; ++j) row[j] = -gamma * from[j];
      row[i] += 1.0;
    }
    return;
  }
  for (int i = 0; i < n; ++i) {
    double* __restrict row = &out[static_cast<size_t>(i) * n];
    for (int j = 0; j < i; ++j) row[j] = -gamma * op.toeplitz_col[i - j];
    row[i] = 1.0 - gamma * op.toeplitz_diag[i];
    for (int j = i + 1; j < n; ++j) row[j] = -gamma * op.toeplitz_row[j - i];
  }
}

// Project one rate vector onto the Krylov basis.
inline void rlf_project_rate(const std::vector<double>& rate,
                             const std::vector<double>& basis, int n, int m,
                             double scale, std::vector<double>& out) {
  out.assign(m, 0.0);
  for (int k = 0; k < m; ++k) {
    const double* column = &basis[static_cast<size_t>(k) * n];
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += rate[i] * column[i];
    out[k] = sum * scale;
  }
}

// Extend the Arnoldi basis on B = (I - gamma L)^-1 from `from` vectors to
// `to`, with one reorthogonalisation pass.  Returns the number of vectors
// actually built, which is smaller than `to` only on a lucky breakdown (an
// invariant subspace, where the remaining modes are already exact).
inline int rlf_arnoldi_extend(const RLF_ShiftedLU& lu, int n, int stride,
                              int from, int to, RLF_KrylovWork& work) {
  const char transpose = 'T', no_transpose = 'N';
  const int one = 1;
  const double d_one = 1.0, d_minus = -1.0, d_zero = 0.0;
  for (int k = from; k < to; ++k) {
    const double* source = &work.basis[static_cast<size_t>(k) * n];
    std::copy(source, source + n, work.vector.begin());
    if (!lu.solve(work.vector.data())) {
      throw std::runtime_error("rlf_solve: shifted solve failed.");
    }
    const int done = k + 1;
    for (int pass = 0; pass < 2; ++pass) {
      F77_CALL(dgemv)(
        &transpose, &n, &done, &d_one, work.basis.data(), &n,
        work.vector.data(), &one, &d_zero, work.projection.data(), &one, 1);
      F77_CALL(dgemv)(
        &no_transpose, &n, &done, &d_minus, work.basis.data(), &n,
        work.projection.data(), &one, &d_one, work.vector.data(), &one, 1);
      for (int i = 0; i < done; ++i) {
        work.hessenberg[static_cast<size_t>(k) * stride + i] +=
          work.projection[i];
      }
    }
    double norm = 0.0;
    for (double value : work.vector) norm += value * value;
    norm = std::sqrt(norm);
    work.hessenberg[static_cast<size_t>(k) * stride + k + 1] = norm;
    if (!(norm > 1e-13)) return k + 1;
    double* target = &work.basis[static_cast<size_t>(k + 1) * n];
    for (int i = 0; i < n; ++i) target[i] = work.vector[i] / norm;
  }
  return to;
}

// Turn the leading m x m block of the Hessenberg matrix into modal curves for
// the three output functionals.  Everything here is O(m^3) or smaller, so it
// can be repeated at several m against one factorisation and one basis.
inline void rlf_reduce_modes(const RLF_Operator& op, double h, double gamma,
                             double beta, int n, int m, int stride,
                             RLF_KrylovWork& work, RLF_Modes& modes) {
  const int one = 1;
  // --- eigendecomposition, taken on H_m rather than on the reduced generator
  //
  // L_m = (I - H_m^-1) / gamma is a function of H_m, so the two have the same
  // eigenvectors and their eigenvalues are related by lam = (1 - 1/mu)/gamma.
  // Working on H_m directly is worth doing twice over: it drops the explicit
  // inversion entirely, and H_m is already upper Hessenberg, so DHSEQR and
  // DTREVC can be called on it without the Hessenberg reduction that DGEEV
  // would spend most of its time on.  This is the largest single cost in the
  // solve after the shifted factorisation, and the search visits it once per
  // rung.
  work.reduced.assign(static_cast<size_t>(m) * m, 0.0);
  for (int j = 0; j < m; ++j) {
    const int rows = std::min(j + 2, m);
    for (int i = 0; i < rows; ++i) {
      work.reduced[static_cast<size_t>(j) * m + i] =
        work.hessenberg[static_cast<size_t>(j) * stride + i];
    }
  }
  work.hess_copy = work.reduced;
  const double h_next = work.hessenberg[static_cast<size_t>(m - 1) * stride + m];

  modes.m = m;
  modes.n = n;
  modes.re.assign(m, 0.0);
  modes.im.assign(m, 0.0);
  work.eigvec.resize(static_cast<size_t>(m) * m);
  int info = 0;
  {
    const char schur = 'S', identity = 'I';
    int lwork = -1;
    double wanted = 0.0;
    F77_CALL(dhseqr)(
      &schur, &identity, &m, &one, &m, work.reduced.data(), &m,
      modes.re.data(), modes.im.data(), work.eigvec.data(), &m, &wanted,
      &lwork, &info, 1, 1);
    lwork = info == 0 ? static_cast<int>(wanted) : 0;
    if (lwork < m) lwork = std::max(m, 1);
    if (static_cast<size_t>(lwork) > work.lapack.size()) {
      work.lapack.resize(lwork);
    }
    F77_CALL(dhseqr)(
      &schur, &identity, &m, &one, &m, work.reduced.data(), &m,
      modes.re.data(), modes.im.data(), work.eigvec.data(), &m,
      work.lapack.data(), &lwork, &info, 1, 1);
    if (info != 0) {
      throw std::runtime_error("rlf_solve: Schur factorisation failed.");
    }
  }
  {
    // HOWMNY = 'B' back-transforms with the Schur vectors already in eigvec,
    // giving eigenvectors of H_m -- and so of the reduced generator.  DTREVC
    // does not normalise the way DGEEV does, which does not matter here: the
    // modal coefficients come from solving against this same basis, so any
    // per-column scaling cancels.
    const char right = 'R', back = 'B';
    int found = 0;
    work.lapack.resize(std::max<size_t>(work.lapack.size(), 3 * m));
    F77_CALL(dtrevc)(
      &right, &back, nullptr, &m, work.reduced.data(), &m, nullptr, &one,
      work.eigvec.data(), &m, &m, &found, work.lapack.data(), &info, 1, 1);
    if (info != 0) {
      throw std::runtime_error("rlf_solve: eigenvector computation failed.");
    }
  }
  // DTREVC leaves the columns unscaled.  Scaling each mode to unit length --
  // jointly over a conjugate pair, so the pair convention survives -- costs
  // nothing, cancels out of every downstream formula, and is within a factor
  // sqrt(m) of the diagonal scaling that minimises the condition number of the
  // eigenbasis.  Without it the condition estimate below measures DTREVC's
  // arbitrary normalisation rather than the conditioning of the expansion.
  for (int j = 0; j < m; ) {
    const int width = (modes.im[j] > 0.0) ? 2 : 1;
    double norm = 0.0;
    for (int c = 0; c < width; ++c) {
      const double* column = &work.eigvec[static_cast<size_t>(j + c) * m];
      for (int i = 0; i < m; ++i) norm += column[i] * column[i];
    }
    norm = std::sqrt(norm);
    if (norm > 0.0) {
      const double scale = 1.0 / norm;
      for (int c = 0; c < width; ++c) {
        double* column = &work.eigvec[static_cast<size_t>(j + c) * m];
        for (int i = 0; i < m; ++i) column[i] *= scale;
      }
    }
    j += width;
  }
  // mu -> lam = (1 - 1/mu) / gamma.  1/mu = conj(mu)/|mu|^2, so the imaginary
  // part keeps its sign and DHSEQR's convention of putting the positive member
  // of a conjugate pair first carries over unchanged.
  for (int j = 0; j < m; ++j) {
    const double a = modes.re[j], b = modes.im[j];
    const double d = a * a + b * b;
    if (!(d > 0.0)) {
      throw std::runtime_error("rlf_solve: projected operator is singular.");
    }
    modes.re[j] = (1.0 - a / d) / gamma;
    modes.im[j] = (b / d) / gamma;
  }

  // Modal coefficients of the reduced state: solve VR c = e_1.  For a
  // conjugate pair DGEEV's real eigenvector columns are [Re v, Im v], and the
  // complex coefficient of v is (c_j - i c_{j+1}); the pair then contributes
  // one real part rather than twice one.
  work.coefficient.assign(m, 0.0);
  work.coefficient[0] = 1.0;
  work.pivot.resize(m);
  std::vector<double> factored = work.eigvec;
  double eigen_norm = 0.0;
  for (int j = 0; j < m; ++j) {
    double column_sum = 0.0;
    for (int i = 0; i < m; ++i) {
      column_sum += std::abs(work.eigvec[static_cast<size_t>(j) * m + i]);
    }
    eigen_norm = std::max(eigen_norm, column_sum);
  }
  F77_CALL(dgetrf)(&m, &m, factored.data(), &m, work.pivot.data(), &info);
  if (info != 0) {
    throw std::runtime_error("rlf_solve: reduced eigenbasis is singular.");
  }
  {
    // O(m^2) on a factorisation that is already paid for.  A blown-up estimate
    // means the modal expansion is a cancelling sum, which is a different
    // failure from an under-resolved Krylov space and is not fixed by more
    // vectors -- see rlf_build_modes.
    const char one_norm = '1';
    double reciprocal = 0.0;
    work.lapack.resize(std::max<size_t>(work.lapack.size(), 4 * m));
    work.iwork.resize(m);
    int cond_info = 0;
    F77_CALL(dgecon)(
      &one_norm, &m, factored.data(), &m, &eigen_norm, &reciprocal,
      work.lapack.data(), work.iwork.data(), &cond_info, 1);
    modes.eigen_condition =
      (cond_info == 0 && reciprocal > 0.0) ? 1.0 / reciprocal : HUGE_VAL;
  }
  {
    const char none = 'N';
    F77_CALL(dgetrs)(
      &none, &m, &one, factored.data(), &m, work.pivot.data(),
      work.coefficient.data(), &m, &info, 1);
  }

  // --- turn each rate functional into a modal curve
  auto modal_curve = [&](const std::vector<double>& projected,
                         RLF_Curve& curve) {
    curve.a.assign(m, 0.0);
    curve.b.assign(m, 0.0);
    for (int j = 0; j < m; ++j) {
      if (modes.im[j] < 0.0) continue;
      const double* column = &work.eigvec[static_cast<size_t>(j) * m];
      double first = 0.0;
      for (int k = 0; k < m; ++k) first += projected[k] * column[k];
      if (modes.im[j] == 0.0) {
        curve.a[j] = first * work.coefficient[j];
        continue;
      }
      const double* next = &work.eigvec[static_cast<size_t>(j + 1) * m];
      double second = 0.0;
      for (int k = 0; k < m; ++k) second += projected[k] * next[k];
      // (first + i second) * (c_j - i c_{j+1})
      curve.a[j] = first * work.coefficient[j] +
                   second * work.coefficient[j + 1];
      curve.b[j] = second * work.coefficient[j] -
                   first * work.coefficient[j + 1];
    }
  };
  auto make_curve = [&](const std::vector<double>& rate, double scale,
                        RLF_Curve& curve) {
    std::vector<double> projected;
    rlf_project_rate(rate, work.basis, n, m, scale * beta, projected);
    modal_curve(projected, curve);
  };
  // Cell quadrature weights.  On a uniform mesh every cell is h wide and the
  // two functionals below reduce to the old scalar scaling.
  std::vector<double> cell_weight;
  if (op.width.empty()) cell_weight.assign(n, h); else cell_weight = op.width;
  std::vector<double> weighted(n);
  for (int i = 0; i < n; ++i) weighted[i] = op.upper_kill[i] * cell_weight[i];
  make_curve(weighted, 1.0, modes.flux);
  make_curve(cell_weight, 1.0, modes.surv);

  // --- shift-and-invert residual, as a modal curve
  //
  // For the shift-and-invert Arnoldi approximation p_m(t) = beta V_m y_m(t),
  // y_m(t) = exp(t L_m) e_1, the defect against the equation being solved is
  //
  //   p_m'(t) - L p_m(t)
  //     = -beta (h_{m+1,m} / gamma) v_{m+1} (e_m^T H_m^{-1} y_m(t)),
  //
  // and since v_{m+1} has unit length the norm of that is the scalar in
  // brackets.  e_m^T H_m^{-1} is one transposed solve against the Hessenberg
  // matrix, and w^T y_m(t) is then a sum over the same modes as every other
  // output, so the whole thing costs O(m^2) once and O(m) per evaluation.
  //
  // This is what makes a single rung self-testing: without it the search can
  // only detect convergence by finding that two successive rungs agree, which
  // means the cheapest possible outcome is two eigensolves rather than one.
  {
    work.residual_row.assign(m, 0.0);
    work.residual_row[m - 1] = 1.0;
    work.pivot.resize(m);
    F77_CALL(dgetrf)(&m, &m, work.hess_copy.data(), &m, work.pivot.data(),
                     &info);
    if (info == 0) {
      const char transposed = 'T';
      F77_CALL(dgetrs)(
        &transposed, &m, &one, work.hess_copy.data(), &m, work.pivot.data(),
        work.residual_row.data(), &m, &info, 1);
    }
    if (info != 0) {
      // No usable estimate; leave the curve empty so residual_at reads zero
      // and the observable two-rung test governs on its own.
      modes.defect.a.clear();
      modes.defect.b.clear();
    } else {
      const double scale = beta * std::abs(h_next) / gamma;
      for (double& value : work.residual_row) value *= scale;
      modal_curve(work.residual_row, modes.defect);
    }
  }
}

// Everything that is read once per solve rather than once per rung: the
// censored-edge curve and the modal form of the state.  Runs off the
// eigendecomposition the last rlf_reduce_modes left in `work`, so settling on
// a dimension costs no second eigensolve.
inline void rlf_materialize_modes(const RLF_Operator& op, double h,
                                  double beta, int n, int m,
                                  RLF_KrylovWork& work, RLF_Modes& modes) {
  {
    std::vector<double> projected;
    std::vector<double> weighted(n);
    for (int i = 0; i < n; ++i) {
      weighted[i] = op.lower_censor[i] *
        (op.width.empty() ? h : op.width[i]);
    }
    rlf_project_rate(weighted, work.basis, n, m, beta, projected);
    modes.censor.a.assign(m, 0.0);
    modes.censor.b.assign(m, 0.0);
    for (int j = 0; j < m; ++j) {
      if (modes.im[j] < 0.0) continue;
      const double* column = &work.eigvec[static_cast<size_t>(j) * m];
      double first = 0.0;
      for (int k = 0; k < m; ++k) first += projected[k] * column[k];
      if (modes.im[j] == 0.0) {
        modes.censor.a[j] = first * work.coefficient[j];
        continue;
      }
      const double* next = &work.eigvec[static_cast<size_t>(j + 1) * m];
      double second = 0.0;
      for (int k = 0; k < m; ++k) second += projected[k] * next[k];
      modes.censor.a[j] = first * work.coefficient[j] +
                          second * work.coefficient[j + 1];
      modes.censor.b[j] = second * work.coefficient[j] -
                          first * work.coefficient[j + 1];
    }
  }

  // --- keep enough to reconstruct the density for the positivity diagnostic
  modes.basis.assign(static_cast<size_t>(m) * n, 0.0);
  for (int j = 0; j < m; ++j) {
    const double* column = &work.basis[static_cast<size_t>(j) * n];
    double* target = &modes.basis[static_cast<size_t>(j) * n];
    for (int i = 0; i < n; ++i) target[i] = beta * column[i];
  }
  modes.state_a.assign(static_cast<size_t>(m) * m, 0.0);
  modes.state_b.assign(static_cast<size_t>(m) * m, 0.0);
  for (int j = 0; j < m; ++j) {
    if (modes.im[j] < 0.0) continue;
    const double* column = &work.eigvec[static_cast<size_t>(j) * m];
    double* first = &modes.state_a[static_cast<size_t>(j) * m];
    if (modes.im[j] == 0.0) {
      for (int i = 0; i < m; ++i) first[i] = column[i] * work.coefficient[j];
      continue;
    }
    const double* next = &work.eigvec[static_cast<size_t>(j + 1) * m];
    double* second = &modes.state_b[static_cast<size_t>(j) * m];
    const double cr = work.coefficient[j], ci = -work.coefficient[j + 1];
    for (int i = 0; i < m; ++i) {
      first[i] = column[i] * cr - next[i] * ci;
      second[i] = column[i] * ci + next[i] * cr;
    }
  }
}

// Build the propagator for one grid, growing the Krylov space until the two
// functionals the likelihood actually reads stop moving.
//
// A fixed dimension cannot serve this model.  The Krylov space has to resolve
// the decay rates that matter over [0, t_max], and how many of those there are
// depends on the parameters: a heavy tail with a fast drift on a fine grid
// needs roughly twice the dimension of an ordinary point, and using the larger
// value everywhere would pay for the worst case on every solve.  Measured
// against a converged reference, a fixed 28 leaves errors of 0.3 nats at
// alpha = 1.13, v = 4.2, nx = 320 -- in the body of the distribution, not the
// tail -- while ordinary points are converged by 24.
//
// The convergence test is on log f(t) and log S(t) rather than on a residual
// bound, because those are the quantities that reach the likelihood, and it is
// deliberately taken in logs so that a Krylov-truncation sign error (which
// exp(tL) cannot produce, L being an M-matrix generator) reads as a large
// discrepancy and buys more vectors instead of being floored away.
//
// The shift-and-invert residual (rlf_reduce_modes) is a genuine single-rung
// error certificate and was tried as the stopping rule, on the reasoning that
// a rung able to retire itself would save the extra rung that "agrees with its
// predecessor" costs.  Measured over the parameter net it does not pay.  It is
// a norm on the density defect, so it is blind to the currency the survivor is
// judged in: |d log S| = |dS|/S, and S reaches e^-30 inside the horizon.  Set
// loose enough to fire early it costs three orders of magnitude of accuracy in
// log S (9.9e-2 against 1.2e-6); set tight enough to be safe it fires at a
// dimension the observable test would have reached anyway, and mean m moves
// from 33.3 to 33.2.  Coarser ladders leaning on it are slower still, since
// the search has to climb through the intermediate rungs regardless and the
// only saving available was ever the last comparison.  It is kept as a
// diagnostic, where it does separate a truncated subspace from a spatial
// problem, and not as a stopping rule.
//
// The tolerance is tight relative to the discretisation error it sits under
// because rlf_cache_solve extrapolates two grids as 5 f - 4 c: the projection
// errors of the two members are independent, so the combination amplifies them
// by around two orders of magnitude, where it holds the spatial error common
// and cancels it.
//
// `m_fixed`, when positive, skips the search entirely and takes a single rung
// at that dimension.  It is used for the second member of a Richardson pair:
// the first member has already established, by the observable test, that the
// dimension resolves the horizon, and the two grids differ by 25% in n on the
// same domain, so they see the same spectrum.  The pair is solved fine member
// first for exactly this reason -- the finer grid is the one that can need
// more vectors, so the reused dimension is an upper bound rather than a
// guess.
inline RLF_Modes rlf_build_modes(const RLF_Operator& op, double h,
                                 const std::vector<double>& p0, double t_max,
                                 RLF_KrylovWork& work, int m_fixed = 0) {
  const int n = op.n;
  const int m_cap = std::min(RLF_KRYLOV_MAX, n);
  if (!(m_cap > 1)) {
    throw std::runtime_error("rlf_solve: grid too small for the propagator.");
  }
  const double gamma = RLF_KRYLOV_SHIFT * t_max;

  rlf_shifted_operator(op, gamma, work.shifted);
  RLF_ShiftedLU lu;
  if (!lu.factor(n, work.shifted)) {
    throw std::runtime_error(
      "rlf_solve: factorisation of the shifted nonlocal operator failed.");
  }

  const int stride = m_cap + 1;
  work.basis.assign(static_cast<size_t>(n) * stride, 0.0);
  work.hessenberg.assign(static_cast<size_t>(stride) * stride, 0.0);
  double beta = 0.0;
  for (double value : p0) beta += value * value;
  beta = std::sqrt(beta);
  if (!(beta > 0.0)) {
    throw std::runtime_error("rlf_solve: start distribution missed the grid.");
  }
  for (int i = 0; i < n; ++i) work.basis[i] = p0[i] / beta;
  work.vector.resize(n);
  work.projection.resize(stride);

  // Probe times: geometric, so the fast transient and the tail both count.
  constexpr int RLF_KRYLOV_PROBES = 12;
  double probe[RLF_KRYLOV_PROBES];
  for (int k = 0; k < RLF_KRYLOV_PROBES; ++k) {
    probe[k] = t_max * std::pow(2.0, (k - (RLF_KRYLOV_PROBES - 1)) * 0.5);
  }
  double previous[2 * RLF_KRYLOV_PROBES];
  bool have_previous = false;

  RLF_Modes modes;
  int built = 0;
  const bool single_rung = m_fixed > 0;
  int m = single_rung
    ? std::min(m_fixed, m_cap)
    : std::min(RLF_KRYLOV_MIN, m_cap);
  while (true) {
    const int grown = rlf_arnoldi_extend(lu, n, stride, built, m, work);
    const bool invariant = grown < m;
    built = grown;
    m = grown;
    rlf_reduce_modes(op, h, gamma, beta, n, m, stride, work, modes);

    // The two functionals are judged in different currencies, because they
    // reach the likelihood in different ways.
    //
    // The exit flux is judged relatively, but against a denominator floored at
    // a fixed fraction of the peak.  A pure relative test on log f chases the
    // deep tail to full precision, and the deep tail is where the modal
    // reconstruction is cancellation noise -- at alpha = 2 the earliest probe
    // can carry a flux nine orders below the peak, and comparing logs there
    // never settles, so the search runs to the cap on a curve that converged
    // twenty vectors earlier.  A pure absolute test has the opposite fault: it
    // under-resolves the leading edge, which carries only a per cent or so of
    // the peak density but is what fixes t0 and the drift.  The floor keeps
    // the test relative across the four decades of density that carry the
    // data, and absolute below that.
    //
    // The survivor is judged relatively, since log S is what the omission and
    // truncation paths read, and unlike the flux it decays smoothly and never
    // collapses into noise.
    double flux_value[RLF_KRYLOV_PROBES], surv_value[RLF_KRYLOV_PROBES];
    double flux_peak = 0.0;
    for (int k = 0; k < RLF_KRYLOV_PROBES; ++k) {
      flux_value[k] = modes.flux_at(probe[k]);
      surv_value[k] = modes.surv_at(probe[k]);
      flux_peak = std::max(flux_peak, std::abs(flux_value[k]));
    }
    const double flux_floor = std::max(flux_peak * RLF_KRYLOV_FLOOR, 1e-300);
    double moved = 0.0;
    for (int k = 0; k < RLF_KRYLOV_PROBES; ++k) {
      const double log_surv = std::log(std::max(surv_value[k], 1e-12));
      if (have_previous) {
        const double denominator =
          std::max(std::abs(flux_value[k]), flux_floor);
        moved = std::max(
          moved, std::abs(flux_value[k] - previous[2 * k]) / denominator);
        moved = std::max(moved, std::abs(log_surv - previous[2 * k + 1]));
      }
      previous[2 * k] = flux_value[k];
      previous[2 * k + 1] = log_surv;
    }

    // An invariant subspace is exact, and the cap is the point past which more
    // vectors cost more than the spatial error they are chasing.  The tolerance is
// on the movement between successive rungs, which overstates the error at the
// rung finally accepted: measured over a random parameter net, 1e-5 here
// leaves 9e-6 after extrapolation, against a spatial error of about 1e-2.
    if (invariant || m >= m_cap || single_rung) break;
    if (have_previous && moved < RLF_KRYLOV_TOL) break;
    have_previous = true;
    m = std::min(m + RLF_KRYLOV_STEP, m_cap);
  }
  rlf_materialize_modes(op, h, beta, n, m, work, modes);
  rlf_last_krylov_dim = m;
  return modes;
}

// One fixed-domain/fixed-grid solve.  `M` is the number of intervals and the
// M-1 interior nodes are x_lo+h, ..., b0-h.  `n_out` sets the resolution of
// the returned time grid when no query times are supplied; it no longer has
// anything to do with how the equation is solved.  The public solver below
// owns all validation and convergence refinement.
// Everything downstream of the discretisation: propagate, then read the three
// functionals the likelihood wants.  Shared by the uniform and graded meshes,
// which differ only in how `op` and `p` were built.
inline RLF_Result rlf_propagate(
    const RLF_Operator& op, double h, const std::vector<double>& p,
    double t_max, int n_out, const std::vector<double>* query_times,
    RLF_KrylovWork* scratch, int m_fixed) {
  const int n = op.n;
  RLF_KrylovWork local;
  RLF_Modes modes = rlf_build_modes(
    op, h, p, t_max, scratch != nullptr ? *scratch : local, m_fixed);

  RLF_Result res;
  res.dx = h;
  res.krylov_dim = modes.m;
  res.eigen_condition = modes.eigen_condition;
  // Largest shift-and-invert residual over a geometric sweep of the horizon.
  // Unlike the mass mismatch, which only sees the component of the Krylov
  // error that breaks the flux/survivor relation, this sees the whole defect,
  // so the two together separate a truncated subspace from a spatial problem.
  {
    constexpr int RLF_RESIDUAL_PROBES = 10;
    for (int k = 1; k <= RLF_RESIDUAL_PROBES; ++k) {
      res.krylov_residual = std::max(
        res.krylov_residual,
        modes.residual_at(t_max * std::pow(2.0, k - RLF_RESIDUAL_PROBES)));
    }
  }
  res.operator_conservation_error = op.conservation_error;

  // Exact antiderivative of the exit flux, for the flux/mass consistency
  // check and for the censored lower-edge pressure.
  double flux_offset = 0.0, censor_offset = 0.0;
  const RLF_Curve cumulative_flux = modes.integrate(modes.flux, flux_offset);
  const RLF_Curve cumulative_censor =
    modes.integrate(modes.censor, censor_offset);
  res.lower_boundary_pressure =
    cumulative_censor.at(modes.re, modes.im, t_max) + censor_offset;

  // Output times: the requested queries, or a uniform grid over [0, t_max].
  std::vector<double> times;
  if (query_times != nullptr) {
    times = *query_times;
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  } else {
    const int levels = std::max(n_out, 1);
    times.reserve(levels + 1);
    for (int k = 0; k <= levels; ++k) {
      times.push_back(t_max * static_cast<double>(k) / levels);
    }
  }

  const size_t count = times.size();
  res.t = times;
  res.pdf.resize(count);
  res.cdf.resize(count);
  res.surv.resize(count);
  double survivor_prev = 1.0;
  std::vector<double> weight;
  for (size_t i = 0; i < count; ++i) {
    const double time = times[i];
    if (!(time > 0.0)) {
      res.pdf[i] = time < 0.0 ? 0.0 : std::max(0.0, modes.flux_at(0.0));
      res.surv[i] = 1.0;
      res.cdf[i] = 0.0;
      continue;
    }
    modes.envelope_at(time, weight);
    res.pdf[i] = std::max(0.0, RLF_Modes::reduce(modes.flux, weight));
    const double raw_survivor = RLF_Modes::reduce(modes.surv, weight);
    // The modal survivor is monotone to Krylov accuracy; clamping keeps a
    // downstream log/interpolation from ever seeing it drift the other way.
    double survivor =
      std::min(survivor_prev, std::min(1.0, std::max(0.0, raw_survivor)));
    survivor_prev = survivor;
    res.surv[i] = survivor;
    res.cdf[i] = 1.0 - survivor;
    const double integrated =
      RLF_Modes::reduce(cumulative_flux, weight) + flux_offset;
    res.flux_mass_mismatch = std::max(
      res.flux_mass_mismatch,
      std::abs((1.0 - raw_survivor) - integrated));
  }

  // Positivity probe, reported as the deepest trough as a fraction of the
  // tallest peak.  exp(t L) is entrywise non-negative because L is an M-matrix
  // generator, so this measures Krylov truncation rather than a property of
  // the discretisation; a handful of geometrically spaced times covers the
  // transient and the tail for a few tens of microseconds.
  constexpr int RLF_DENSITY_PROBES = 8;
  double lowest = 0.0, highest = 0.0;
  for (int k = 1; k <= RLF_DENSITY_PROBES; ++k) {
    const double time =
      t_max * std::pow(2.0, k - RLF_DENSITY_PROBES);
    modes.density_extremes_at(time, lowest, highest);
  }
  res.min_density = highest > 0.0 ? lowest / highest : 0.0;
  return res;
}

inline RLF_Result rlf_solve_fixed_grid(
    const RLF_Model& m, double t_max, int M, int n_out, double lower_extent,
    const std::vector<double>* query_times = nullptr,
    RLF_KrylovWork* scratch = nullptr, int m_fixed = 0) {
  const double x_lo = -lower_extent;
  const double h = (m.b0 - x_lo) / M;
  const int n = M - 1;
  const RLF_Operator op = build_rlf_operator(m, n, h);
  std::vector<double> p;
  rlf_initial_density(m, x_lo, h, n, p);
  RLF_Result res = rlf_propagate(
    op, h, p, t_max, n_out, query_times, scratch, m_fixed);
  res.x_lo = x_lo;
  res.nx_used = M;
  return res;
}

// Start distribution on a graded mesh: the uniform-[0, A] mass, or the point
// mass at zero split linearly between the two cells straddling it so that the
// start point is not silently snapped to a cell centre.
inline void rlf_initial_density_mesh(const RLF_Model& m, const RLF_Mesh& mesh,
                                     std::vector<double>& p) {
  const int n = mesh.n;
  p.assign(n, 0.0);
  if (m.z0 > 0.0) {
    double mass = 0.0;
    for (int i = 0; i < n; ++i) {
      const double overlap =
        std::min(mesh.face[i + 1], m.z0) - std::max(mesh.face[i], 0.0);
      if (overlap > 0.0) {
        p[i] = overlap / (m.z0 * mesh.width[i]);
        mass += overlap / m.z0;
      }
    }
    if (!(mass > 0.0)) {
      throw std::runtime_error("rlf_solve: start distribution missed the grid.");
    }
    for (int i = 0; i < n; ++i) p[i] /= mass;
    return;
  }
  int i = static_cast<int>(
    std::lower_bound(mesh.node.begin(), mesh.node.end(), 0.0) -
    mesh.node.begin());
  if (i <= 0) {
    p[0] = 1.0 / mesh.width[0];
  } else if (i >= n) {
    p[n - 1] = 1.0 / mesh.width[n - 1];
  } else {
    const double w =
      (0.0 - mesh.node[i - 1]) / (mesh.node[i] - mesh.node[i - 1]);
    p[i - 1] = (1.0 - w) / mesh.width[i - 1];
    p[i] = w / mesh.width[i];
  }
}

inline RLF_Result rlf_solve_graded(
    const RLF_Model& m, double t_max, const RLF_Mesh& mesh, int n_out,
    const std::vector<double>* query_times = nullptr,
    RLF_KrylovWork* scratch = nullptr, int m_fixed = 0) {
  const RLF_Operator op = build_rlf_operator_graded(
    m, mesh, scratch != nullptr ? &scratch->operator_G : nullptr);
  std::vector<double> p;
  rlf_initial_density_mesh(m, mesh, p);
  RLF_Result res = rlf_propagate(
    op, mesh.width.back(), p, t_max, n_out, query_times, scratch, m_fixed);
  res.x_lo = mesh.face.front();
  res.nx_used = mesh.n;
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
  // Relative, because min_density now is: the Krylov reconstruction leaves a
  // trough a few parts in 10^7 of the peak, which is not a resolution problem
  // and must not send an ordinary solve down the refinement path.
  constexpr double RLF_FAST_UNDERSHOOT_TOL = -1e-5;
  constexpr int RLF_FAST_MIN_INTERVALS = 160;
  if (M < RLF_FAST_MIN_INTERVALS) return false;
  if (!(base.lower_boundary_pressure <= RLF_FAST_PRESSURE_TOL)) return false;
  if (base.min_density < RLF_FAST_UNDERSHOOT_TOL) return false;
  // Two independent statements about the propagator itself: the defect against
  // the equation, and whether the modal expansion that evaluates it is a
  // cancelling sum.  Neither is a stopping rule for the Krylov search -- see
  // rlf_build_modes -- but either being out of range is a reason to make a
  // solve prove itself by refinement rather than waving it through.
  if (base.krylov_residual > RLF_KRYLOV_RESID) return false;
  if (!(base.eigen_condition < RLF_KRYLOV_CONDITION)) return false;
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
                            int M = 200, int n_out = 400, bool adaptive = true,
                            const std::vector<double>* query_times = nullptr,
                            double lower_extent_override = -1.0,
                            RLF_KrylovWork* scratch = nullptr,
                            int m_fixed = 0) {
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
  n_out = std::max(n_out, 50);

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

  auto solve = [&](int M_use, double extent_use,
                   const std::vector<double>* queries) {
    return rlf_solve_fixed_grid(
      m, t_max, M_use, n_out, extent_use, queries, scratch, m_fixed);
  };

  RLF_Result accepted =
    solve(intervals, lower_extent, adaptive ? nullptr : query_times);
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
      solve(expanded_intervals, expanded_extent, nullptr);
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
      solve(refined_intervals, lower_extent, nullptr);
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
  // Core cells of the graded mesh; the tail below the core is sized
  // logarithmically on top of this.  See rlf_auto_mesh.
  int nx = 70;
  // Spacing of the returned time grid on the complete-grid path, where later
  // arbitrary queries are answered by interpolation.  It is not a time step:
  // the propagator carries no time discretisation, so this only trades
  // interpolation error
  // against the size of a cache entry.  The sparse path ignores it and
  // evaluates the modal form at the requested times directly.
  double dt_target = 1.6e-2;
  int samples_min = 50;
  int samples_max = 20000;
  bool adaptive = false;
  bool sparse_output = true;
  // Solve long-horizon rows separately from the bulk; see rlf_horizon_bucket.
  bool horizon_split = true;
  // Pair each solve with a finer one and Richardson-extrapolate; see
  // rlf_cache_solve.
  bool richardson = true;
  double richardson_ratio = 1.4;

  int samples_for(double t_max) const {
    const double wanted =
      std::ceil(t_max / std::max(dt_target, 1e-8));
    const int bounded = wanted < static_cast<double>(samples_max)
      ? static_cast<int>(wanted) : samples_max;
    return std::max(samples_min, bounded);
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
  // Krylov scratch, held here so the O(n^2) shifted operator and the Arnoldi
  // basis are allocated once per cache rather than once per solve.
  RLF_KrylovWork scratch;
  std::vector<Entry> entries;
  std::unordered_map<Key, int, KeyHash> index;
  std::vector<int> row_group;
  size_t solve_count = 0;

  // Reuse guard for rlf_prepare_rows.  Repeated raw calls inside one
  // particle (d, p, then time-only variants at new GL nodes) rebuild the
  // identical grouping; when this state matches the incoming rows exactly the
  // regrouping and cache lookup are skipped.  Pointer identity alone is not
  // sufficient -- callers reuse the same column arrays and mutate their
  // contents -- so an exact value snapshot of every grouping input is kept
  // alongside the pointers (about 7 doubles plus one flag integer per row,
  // which is small next to one cached march).
  bool rows_prepared = false;
  int prepared_n_rows = -1;
  const double* prepared_rt = nullptr;
  const int* prepared_isok = nullptr;
  std::vector<const double*> prepared_cols;
  std::vector<double> prepared_values;
  std::vector<int> prepared_isok_values;

  void new_particle() {
    entries.clear();
    index.clear();
    row_group.clear();
    solve_count = 0;
    rows_prepared = false;
    prepared_n_rows = -1;
    prepared_rt = nullptr;
    prepared_isok = nullptr;
    prepared_cols.clear();
    prepared_values.clear();
    prepared_isok_values.clear();
  }
};

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
// Convergence order at the absorbing edge, 1 + alpha/2 rounded to a single
// number over the supported range.  Using the exact alpha-dependent exponent
// was measured and is not better: it improves the median and costs the tail.
inline double RLF_EDGE_ORDER = 1.5;

// `ratio` is the effective error ratio between the two members, r^q for a
// scheme of order q -- not the mesh ratio itself.
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
                            Entry& out, RLF_KrylovWork* scratch = nullptr) {
  RLF_Model model;
  model.v = key.v;
  model.sigma = 1.0;
  model.alpha = key.alpha;
  model.b0 = key.b;
  model.z0 = key.A;

  const int nx = std::max(grid.nx, 30);
  const int samples = grid.samples_for(t_max);
  const int nx_fine =
    static_cast<int>(std::lround(nx * grid.richardson_ratio));

  if (grid.adaptive) {
    // Validation path: the uniform grid with its self-checking domain and
    // resolution refinement, kept so the graded solver has something
    // independent to be measured against.
    const RLF_Result result = rlf_solve(
      model, t_max, nx, samples, true, query_times, -1.0, scratch);
    rlf_store_result(key, t_max, result, query_times == nullptr, out);
    return;
  }

  if (grid.richardson && nx_fine > nx) {
    // What is left after grading is the absorbing edge.  The killed density
    // vanishes like (b0 - x)^(alpha/2), an algebraic singularity, and the
    // measured order of the scheme is 1 + alpha/2 rather than 2 because of it.
    // Refining into the edge does not help -- the profile is scale free, so a
    // geometric mesh meets the same relative closure error at every level, and
    // it costs accuracy elsewhere by making the widths vary.  Extrapolating
    // does help, and the exponent is known, so the pair is combined at
    // r^1.5 rather than at r.  Measured over the parameter net that is worth
    // roughly a factor of three in the ninetieth-percentile error for a factor
    // of 1.5 in cost, which is the better half of the frontier.
    //
    // The two members share the domain and the tail grading and differ only in
    // the core cell count, so their meshes are nested in the sense the
    // combination needs, and the fine member is solved first so that the
    // coarse one can reuse its Krylov dimension without re-climbing.
    const RLF_Mesh fine_mesh = rlf_auto_mesh(model, t_max, nx_fine);
    const RLF_Mesh coarse_mesh = rlf_auto_mesh(model, t_max, nx);
    const RLF_Result fine = rlf_solve_graded(
      model, t_max, fine_mesh, samples, query_times, scratch);
    const RLF_Result coarse = rlf_solve_graded(
      model, t_max, coarse_mesh, samples, query_times, scratch,
      fine.krylov_dim);
    RLF_Result blended;
    const double ratio = static_cast<double>(fine_mesh.n) / coarse_mesh.n;
    const bool ok = rlf_extrapolate(
      coarse, fine, std::pow(ratio, RLF_EDGE_ORDER), blended);
    rlf_store_result(key, t_max, ok ? blended : fine,
                     query_times == nullptr, out);
    return;
  }

  const RLF_Mesh mesh = rlf_auto_mesh(model, t_max, nx);
  const RLF_Result result =
    rlf_solve_graded(model, t_max, mesh, samples, query_times, scratch);
  rlf_store_result(key, t_max, result, query_times == nullptr, out);
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
inline bool entry_has_query(const Entry& entry, double query_time) {
  if (!std::isfinite(query_time) || !(query_time > 0.0)) return false;
  if (entry.complete_grid) return entry.t_max >= query_time;
  return std::binary_search(entry.t.begin(), entry.t.end(), query_time);
}

inline bool same_parameter_key(const Key& a, const Key& b) {
  return a.v == b.v && a.alpha == b.alpha && a.b == b.b && a.A == b.A;
}

inline int cache_get(SolveCache& cache, const Key& key, double t_need,
                     double query_time = R_NaN) {
  if (!(t_need > 0.0)) t_need = 1e-3;
  if (std::isfinite(query_time) && query_time > 0.0) {
    for (size_t i = 0; i < cache.entries.size(); ++i) {
      if (same_parameter_key(cache.entries[i].key, key) &&
          entry_has_query(cache.entries[i], query_time)) {
        return static_cast<int>(i);
      }
    }
  }
  const auto found = cache.index.find(key);
  if (found != cache.index.end()) {
    Entry& entry = cache.entries[found->second];
    if (entry_has_query(entry, query_time)) return found->second;
    if (entry.complete_grid && entry.t_max >= t_need) return found->second;
    // Sparse entries answer exact prepared queries. Otherwise upgrade once to a
    // complete grid covering max(t_need, entry.t_max).
    Entry replacement;
    try {
      rlf_cache_solve(key, std::max(t_need, entry.t_max), cache.grid, nullptr,
                      replacement, &cache.scratch);
    } catch (const std::exception&) {
      return -1;
    }
    entry = std::move(replacement);
    ++cache.solve_count;
    return found->second;
  }

  Entry fresh;
  try {
    rlf_cache_solve(key, t_need, cache.grid, nullptr, fresh, &cache.scratch);
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
  // Two kinds of work are separated here.  A key with no cached entry is
  // solved directly.  A key whose entry exists but cannot answer this
  // request -- sparse output missing one of the requested query times, or a
  // horizon beyond the stored solve -- is upgraded once to a complete grid
  // covering max(old t_max, t_need); replacing it with another sparse solve
  // would re-solve on every later query set instead of converging after one.
  std::vector<size_t> upgrades;
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
      upgrades.push_back(i);
      continue;
    }
    missing.push_back(i);
  }

  // Upgrades run with nullptr queries so the replacement is a complete grid;
  // the old entry is kept when unresolvable, and only this request's rows are
  // floored.
  for (size_t i : upgrades) {
    double t_need = horizons[i];
    if (!(t_need > 0.0)) t_need = 1e-3;
    const int slot = cache.index[keys[i]];
    Entry replacement;
    try {
      rlf_cache_solve(keys[i], std::max(t_need, cache.entries[slot].t_max),
                      cache.grid, nullptr, replacement, &cache.scratch);
    } catch (const std::exception&) {
      continue;  // keep the previous entry; rows stay floored at -1
    }
    cache.entries[slot] = std::move(replacement);
    ++cache.solve_count;
    out_indices[i] = slot;
  }

  for (size_t source : missing) {
    const std::vector<double>* queries =
      query_times == nullptr ? nullptr : &(*query_times)[source];
    Entry replacement;
    try {
      rlf_cache_solve(keys[source], horizons[source], cache.grid, queries,
                      replacement, &cache.scratch);
    } catch (const std::exception&) {
      // Unresolvable parameter point: leave it uncached and unindexed and
      // report -1, so only this key's rows are floored.  See cache_get().
      out_indices[source] = -1;
      continue;
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

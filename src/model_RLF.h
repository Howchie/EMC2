#ifndef model_RLF_h
#define model_RLF_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Race Lévy Flight (RLF) model and Space-Fractional Fokker-Planck (SFFPE) solver.
//
// The Race Lévy Flight process models decision making with heavy-tailed evidence
// accumulation noise (alpha-stable distribution with 1 < alpha <= 2).
//
// Governing SFFPE:
//   dp/dt = -v * dp/dx + K_alpha * d^alpha p / d|x|^alpha
//
// Boundary conditions:
//   p(a, t) = 0   (fixed absorbing threshold at a > 0)
//   p(-L, t) = 0  (absorbing lower domain limit at -L)
//
// Option A solver:
//   Fixed boundary => spatial matrix L_op is time-invariant.
//   The Crank-Nicolson matrix (I - 0.5 * dt * L_op) is factorised ONCE via dense
//   LU decomposition with row pivoting (O(M^3)).
//   Time marching solves (I - 0.5 * dt * L_op) q^{n+1} = (I + 0.5 * dt * L_op) q^n
//   in O(M^2) per step, yielding solve times of ~1-5 ms.
//
// Simulator:
//   Paths are simulated using the Chambers-Mallows-Stuck (CMS) algorithm for
//   symmetric alpha-stable random variables.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <limits>
#include <stdexcept>

namespace rlf {

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

// ---------------------------------------------------------------------------
// Dense LU Factorization with Partial Pivoting
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
        double val = std::abs(LU[i * n + k]);
        if (val > max_val) {
          max_val = val;
          max_i = i;
        }
      }
      if (max_val < 1e-15) return false; // Singular matrix

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

  void solve(const double* __restrict b, double* __restrict x, double* __restrict y) const {
    // Forward substitution L y = P b
    for (int i = 0; i < n; ++i) {
      double sum = b[pvt[i]];
      const double* __restrict lu_i = &LU[i * n];
      for (int j = 0; j < i; ++j) {
        sum -= lu_i[j] * y[j];
      }
      y[i] = sum;
    }
    // Back substitution U x = y
    for (int i = n - 1; i >= 0; --i) {
      double sum = y[i];
      const double* __restrict lu_i = &LU[i * n];
      for (int j = i + 1; j < n; ++j) {
        sum -= lu_i[j] * x[j];
      }
      x[i] = sum / lu_i[i];
    }
  }
};

// ---------------------------------------------------------------------------
// Compute Grünwald-Letnikov (GL) Fractional Weights g_k^{(\alpha)}
//   g_0 = 1,  g_k = (1 - (alpha + 1)/k) * g_{k-1}
// ---------------------------------------------------------------------------
inline std::vector<double> compute_gl_weights(double alpha, int max_k) {
  std::vector<double> g(max_k + 1, 0.0);
  g[0] = 1.0;
  for (int k = 1; k <= max_k; ++k) {
    g[k] = (1.0 - (alpha + 1.0) / static_cast<double>(k)) * g[k - 1];
  }
  return g;
}

// ---------------------------------------------------------------------------
// RLF PDE Model Parameters & Solver Output
// ---------------------------------------------------------------------------
struct RLF_Model {
  double v = 1.0;       // drift rate
  double sigma = 1.0;   // noise scale parameter
  double alpha = 1.7;   // stability index (1 < alpha <= 2)
  double b0 = 1.0;      // decision threshold / barrier (fixed)
  double z0 = 0.0;      // start variability upper bound Uniform(0, z0)
};

struct RLF_Result {
  std::vector<double> t, pdf, cdf;
};

// Linear lookup interpolation
inline double grid_lookup(const std::vector<double>& tg,
                          const std::vector<double>& vg, double t) {
  const int n = static_cast<int>(tg.size());
  if (n == 0) return 0.0;
  if (t <= tg.front()) return vg.front();
  if (t >= tg.back())  return vg.back();
  const int j = static_cast<int>(
    std::lower_bound(tg.begin(), tg.end(), t) - tg.begin());
  const double w = (t - tg[j - 1]) / (tg[j] - tg[j - 1]);
  return vg[j - 1] + w * (vg[j] - vg[j - 1]);
}

// ---------------------------------------------------------------------------
// Option A PDE Solver for Race Lévy Flight
// ---------------------------------------------------------------------------
inline RLF_Result rlf_solve(const RLF_Model& m, double t_max, int M = 200, int Nt = 400) {
  RLF_Result res;
  if (m.alpha <= 1.0 || m.alpha > 2.0) {
    throw std::invalid_argument("rlf_solve: alpha must be in (1, 2].");
  }
  if (m.sigma <= 0.0) {
    throw std::invalid_argument("rlf_solve: sigma must be positive.");
  }
  if (m.b0 <= 0.0) {
    throw std::invalid_argument("rlf_solve: boundary b0 must be positive.");
  }
  if (t_max <= 0.0) t_max = 1.0;

  M = std::max(M, 30);
  Nt = std::max(Nt, 50);

  // Extend lower boundary -L based on diffusion/dispersion scale
  const double scale_t = std::pow(t_max, 1.0 / m.alpha);
  const double L = std::max(m.b0, 8.0 * m.sigma * scale_t);
  const double h = (m.b0 + L) / M;
  const double dt = t_max / Nt;

  const int n_int = M - 1; // interior grid nodes x_1 .. x_{M-1}
  std::vector<double> L_op(n_int * n_int, 0.0);

  const auto g_gl = compute_gl_weights(m.alpha, M + 5);

  double C_alpha;
  if (std::abs(m.alpha - 2.0) < 1e-6) {
    C_alpha = 0.25 * m.sigma * m.sigma / (h * h);
  } else {
    C_alpha = - (0.25 * std::pow(m.sigma, m.alpha)) /
                (2.0 * std::cos(M_PI * m.alpha / 2.0) * std::pow(h, m.alpha));
  }

  for (int i = 0; i < n_int; ++i) {
    // Upwind advection: -v * (p_i - p_{i-1}) / h
    L_op[i * n_int + i] -= m.v / h;
    if (i > 0) {
      L_op[i * n_int + (i - 1)] += m.v / h;
    }

    // Shifted GL Riesz fractional derivative
    const int node_i = i + 1;
    for (int j = 0; j < n_int; ++j) {
      const int node_j = j + 1;
      const int k_left = node_i - node_j + 1;
      if (k_left >= 0 && k_left < static_cast<int>(g_gl.size())) {
        L_op[i * n_int + j] += C_alpha * g_gl[k_left];
      }
      const int k_right = node_j - node_i + 1;
      if (k_right >= 0 && k_right < static_cast<int>(g_gl.size())) {
        L_op[i * n_int + j] += C_alpha * g_gl[k_right];
      }
    }
  }

  // Crank-Nicolson matrices:
  //   LHS = I - 0.5 * dt * L_op
  //   RHS = I + 0.5 * dt * L_op
  std::vector<double> LHS(n_int * n_int, 0.0);
  std::vector<double> RHS_mat(n_int * n_int, 0.0);
  for (int i = 0; i < n_int; ++i) {
    for (int j = 0; j < n_int; ++j) {
      const double val = L_op[i * n_int + j];
      LHS[i * n_int + j]     = (i == j ? 1.0 : 0.0) - 0.5 * dt * val;
      RHS_mat[i * n_int + j] = (i == j ? 1.0 : 0.0) + 0.5 * dt * val;
    }
  }

  // Pre-factorize LHS once
  RLF_DenseLU lu;
  if (!lu.factor(n_int, LHS)) {
    throw std::runtime_error("rlf_solve: LU factorization of spatial operator failed.");
  }

  // Initial condition: Uniform(0, z0) or point start at 0
  std::vector<double> p(n_int, 0.0);
  if (m.z0 > 1e-10) {
    const double Zhi = std::min(m.z0, m.b0);
    for (int i = 0; i < n_int; ++i) {
      const double xi_a = -L + i * h;
      const double xi_b = -L + (i + 1) * h;
      const double ov = std::min(xi_b, Zhi) - std::max(xi_a, 0.0);
      if (ov > 0.0) p[i] = (ov / Zhi) / h;
    }
  } else {
    // Point start at x = 0
    const int i_start = std::max(0, std::min(n_int - 1, static_cast<int>(std::round(L / h)) - 1));
    p[i_start] = 1.0 / h;
  }

  res.t.resize(Nt + 1);
  res.pdf.resize(Nt + 1);
  res.cdf.resize(Nt + 1);

  double mass = 0.0;
  for (double val : p) mass += val * h;
  res.t[0] = 0.0;
  res.cdf[0] = std::max(0.0, 1.0 - mass);
  res.pdf[0] = 0.0;

  std::vector<double> rhs(n_int, 0.0), p_next(n_int, 0.0), y_tmp(n_int, 0.0);

  for (int step = 1; step <= Nt; ++step) {
    // rhs = RHS_mat * p
    for (int i = 0; i < n_int; ++i) {
      double sum = 0.0;
      const double* __restrict r_row = &RHS_mat[i * n_int];
      for (int j = 0; j < n_int; ++j) {
        sum += r_row[j] * p[j];
      }
      rhs[i] = sum;
    }

    lu.solve(rhs.data(), p_next.data(), y_tmp.data());
    p = p_next;

    mass = 0.0;
    for (double val : p) mass += val * h;
    const double current_cdf = std::max(0.0, 1.0 - mass);
    const double current_t = step * dt;

    res.t[step] = current_t;
    res.cdf[step] = current_cdf;
    res.pdf[step] = std::max(0.0, (current_cdf - res.cdf[step - 1]) / dt);
  }

  return res;
}

// ---------------------------------------------------------------------------
// Chambers-Mallows-Stuck (CMS) Symmetric Alpha-Stable Random Variable Generator
// ---------------------------------------------------------------------------
inline double sample_stable_cms(double alpha, std::mt19937& rng) {
  if (alpha == 2.0) {
    std::normal_distribution<double> norm(0.0, std::sqrt(2.0));
    return norm(rng);
  }
  std::uniform_real_distribution<double> unif(-M_PI_2 + 1e-12, M_PI_2 - 1e-12);
  std::exponential_distribution<double> exp1(1.0);

  const double U = unif(rng);
  const double W = exp1(rng);

  if (std::abs(alpha - 1.0) < 1e-8) {
    return std::tan(U);
  }

  const double aU = alpha * U;
  const double om_a_U = (1.0 - alpha) * U;
  const double num = std::sin(aU);
  const double den = std::pow(std::cos(U), 1.0 / alpha);
  const double term2 = std::pow(std::cos(om_a_U) / W, (1.0 - alpha) / alpha);
  return (num / den) * term2;
}

// ---------------------------------------------------------------------------
// Efficient C++ Path Simulator for Race Lévy Flight
// ---------------------------------------------------------------------------
inline std::vector<double> simulate_rlf_hit_times(
    int n_sims, double v, double sigma, double alpha, double b0, double z0,
    double t_max = 5.0, double dt = 0.001, unsigned int seed = 42) {

  if (alpha <= 0.0 || alpha > 2.0) {
    throw std::invalid_argument("simulate_rlf_hit_times: alpha must be in (0, 2].");
  }
  if (sigma <= 0.0) {
    throw std::invalid_argument("simulate_rlf_hit_times: sigma must be positive.");
  }
  if (b0 <= 0.0) {
    throw std::invalid_argument("simulate_rlf_hit_times: b0 must be positive.");
  }
  if (dt <= 0.0) dt = 0.001;

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif_start(0.0, std::max(0.0, z0));

  const int n_steps = static_cast<int>(std::ceil(t_max / dt));
  const double scale = sigma * std::pow(dt, 1.0 / alpha);

  std::vector<double> hit_times(n_sims, std::numeric_limits<double>::quiet_NaN());

  for (int s = 0; s < n_sims; ++s) {
    double x = (z0 > 0.0) ? unif_start(rng) : 0.0;
    if (x >= b0) {
      hit_times[s] = 0.0;
      continue;
    }
    double t = 0.0;
    for (int step = 0; step < n_steps; ++step) {
      const double z = sample_stable_cms(alpha, rng);
      x += v * dt + scale * z;
      t += dt;
      if (x >= b0) {
        hit_times[s] = t;
        break;
      }
    }
  }
  return hit_times;
}

} // namespace rlf

#endif // model_RLF_h

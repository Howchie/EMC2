#ifndef utils_reducible_diffusion_h
#define utils_reducible_diffusion_h

#define _USE_MATH_DEFINES
#include "utility_functions.h"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <functional>
#include <limits>
#include <cstddef>
#include <utility>
#include <deque>
#include <Rcpp.h>
#include <quadmath.h>
#include <stdexcept>
#include <iostream> // For error messages
using namespace Rcpp;

struct RD_Params;

using KernelFn = std::function<double(double,double,const RD_Params&)>;
using ForcingFn = std::function<double(double,const RD_Params&)>;
using AbelFn = std::function<double(double,const RD_Params&)>;
using BoundaryDecayFn = std::function<double(double, const RD_Params &)>;

// Prepare and scale parameters for the OU hitting time problem.
struct RD_Params {
  double b0, t_scaled, zL_scaled, zU_scaled, b_scaled, binf, c, lambda, theta, tau, pow, omega, beta_prime;
  double mu, sigma, z0, v_max, theta_max;
  bool fixed_b, sp_var;
  BoundaryDecayFn boundary_fn;
  std::vector<double> boundary_params;
};

inline double default_boundary_decay(double t, const RD_Params& pars) {
  return exp_decay_scalar(t, pars.b0, pars.binf, pars.tau, pars.pow);
}

inline double fixed_boundary_decay(double, const RD_Params& pars) {
  return pars.b0;
}

inline double evaluate_boundary_decay(double t, const RD_Params& pars) {
  if (pars.boundary_fn) {
    return pars.boundary_fn(t, pars);
  }
  if (pars.fixed_b) {
    return pars.b0;
  }
  return default_boundary_decay(t, pars);
}

// Threshold for switching to Abel approximation, based on the scaled time.
// The paper notes solutions are "visually indistinguishable" up to t=0.02.
const double SMALL_T_SCALED_THRESHOLD = 0.02;
const double FPM_EPSILON = 1e-12;

struct BoundaryDecayCache {
    double h = 0.0;
    double inv_h = 0.0;
    double tol = 0.0;
    std::vector<double> beta_nodes;
    std::vector<double> beta_midpoints;
    std::vector<double> beta_prime_nodes;
    std::vector<double> beta_prime_midpoints;

    bool empty() const {
        return beta_nodes.empty() && beta_midpoints.empty();
    }

    bool has_derivatives() const {
        return !beta_prime_nodes.empty() || !beta_prime_midpoints.empty();
    }

    double lookup(double v) const {
        if (inv_h <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double scaled = v * inv_h;
        const double idx = std::round(scaled);
        if (!beta_nodes.empty() && std::abs(scaled - idx) <= tol) {
            const double idx_clamped = std::clamp(idx, 0.0, static_cast<double>(beta_nodes.size() - 1));
            return beta_nodes[static_cast<std::size_t>(idx_clamped)];
        }
        if (!beta_midpoints.empty()) {
            const double base = std::floor(scaled);
            const double mid_center = base + 0.5;
            if (std::abs(scaled - mid_center) <= tol) {
                const double mid_idx = std::clamp(base, 0.0, static_cast<double>(beta_midpoints.size() - 1));
                return beta_midpoints[static_cast<std::size_t>(mid_idx)];
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    double lookup_prime(double v) const {
        if (!has_derivatives() || inv_h <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double scaled = v * inv_h;
        const double idx = std::round(scaled);
        if (!beta_prime_nodes.empty() && std::abs(scaled - idx) <= tol) {
            const double idx_clamped = std::clamp(idx, 0.0, static_cast<double>(beta_prime_nodes.size() - 1));
            return beta_prime_nodes[static_cast<std::size_t>(idx_clamped)];
        }
        if (!beta_prime_midpoints.empty()) {
            const double base = std::floor(scaled);
            const double mid_center = base + 0.5;
            if (std::abs(scaled - mid_center) <= tol) {
                const double mid_idx = std::clamp(base, 0.0, static_cast<double>(beta_prime_midpoints.size() - 1));
                return beta_prime_midpoints[static_cast<std::size_t>(mid_idx)];
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    }
};

// Right-end quadratic fit (last 3 points) => nu'(tau), nu''(tau), and nu(mid) via Newton form
struct RightQuad {
    double nu_prime;   // nu'(tau)
    double nu_second;  // nu''(tau)
    double nu_at_mid;  // quadratic nu at midpoint of last panel
};
RightQuad right_end_quadratic(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     double mid) {
    const int n = (int)x.size();
    RightQuad out{0.0, 0.0, y.back()};
    if (n < 3) return out;
    const double x0 = x[n-3], x1 = x[n-2], x2 = x[n-1];
    const double y0 = y[n-3], y1 = y[n-2], y2 = y[n-1];

    // Divided differences (Newton) on uneven grid
    const double f01  = (y1 - y0) / (x1 - x0);
    const double f12  = (y2 - y1) / (x2 - x1);
    const double f012 = (f12 - f01) / (x2 - x0);

    // nu'(x2) and nu''(x2) for the quadratic that interpolates the last 3 points
    out.nu_second = 2.0 * f012;
    out.nu_prime  = f01 + f012 * ((x2 - x0) + (x2 - x1));  // p'(x2) = f[x0,x1] + f[x0,x1,x2]((x2-x0)+(x2-x1))

    // nu(mid) from the same quadratic
    const double mid_m_x0 = mid - x0;
    const double mid_m_x1 = mid - x1;
    out.nu_at_mid = y0 + f01 * (mid_m_x0) + f012 * (mid_m_x0) * (mid_m_x1);

    return out;
}

double local_quad_derivative(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  int j) {
    const int n = static_cast<int>(x.size());
    if (n < 3) return 0.0;
    int c = std::max(2, std::min(j, n - 1)); // ensure we can pick (a,b,c)
    int b = c - 1;
    int a = c - 2;
    const double x0 = x[a], x1 = x[b], x2 = x[c];
    const double y0 = y[a], y1 = y[b], y2 = y[c];
    const double f01  = (y1 - y0) / (x1 - x0);
    const double f12  = (y2 - y1) / (x2 - x1);
    const double f012 = (f12 - f01) / (x2 - x0);
    return f01 + f012 * ((x2 - x0) + (x2 - x1));
  };

double quad_interp(
    const std::vector<double>& grid,
    const std::vector<double>& nu,
    int i,            // panel index for [i, i+1]
    double t_q    // query inside [theta[i], theta[i+1]]
) {
    const int N = static_cast<int>(grid.size()) - 1;
    // Choose a local triple around the panel:
    // interior panels: (i-1,i,i+1); near left: (0,1,2); near right: (N-2,N-1,N)
    int a, b, c;
    if (i == 0) { a = 0; b = 1; c = 2; }
    else if (i == N-1) { a = N-2; b = N-1; c = N; }
    else { a = i-1; b = i; c = i+1; }

    const double xa = grid[a], xb = grid[b], xc = grid[c];
    const double ya = nu[a],    yb = nu[b],    yc = nu[c];

    // Lagrange basis at theta_q
    const double Lab = (t_q - xb) / (xa - xb);
    const double Lac = (t_q - xc) / (xa - xc);
    const double Lba = (t_q - xa) / (xb - xa);
    const double Lbc = (t_q - xc) / (xb - xc);
    const double Lca = (t_q - xa) / (xc - xa);
    const double Lcb = (t_q - xb) / (xc - xb);

    return ya * Lab * Lac + yb * Lba * Lbc + yc * Lca * Lcb;
}

inline BoundaryDecayCache make_boundary_decay_cache(double t_max, int num_steps, const RD_Params& pars,
                                                    const BoundaryDecayFn& transform) {
    BoundaryDecayCache cache;
    if (pars.fixed_b || num_steps <= 0 || t_max <= 0.0 || !transform) {
        return cache;
    }
    const double h = t_max / num_steps;
    if (h <= 0.0) {
        return cache;
    }
    cache.h = h;
    cache.inv_h = 1.0 / h;
    cache.tol = std::max(1e-10, 32.0 * std::numeric_limits<double>::epsilon() *
                                     std::max(1.0, static_cast<double>(num_steps)));
    cache.beta_nodes.resize(num_steps + 1);
    for (int k = 0; k <= num_steps; ++k) {
        const double t_k = h * k;
        cache.beta_nodes[k] = transform(t_k, pars);
    }
    cache.beta_midpoints.resize(num_steps);
    for (int k = 0; k < num_steps; ++k) {
        const double t_mid = h * (k + 0.5);
        cache.beta_midpoints[k] = transform(t_mid, pars);
    }

    cache.beta_prime_nodes.assign(num_steps + 1, 0.0);
    if (num_steps == 1) {
        const double slope = (cache.beta_nodes[1] - cache.beta_nodes[0]) / h;
        cache.beta_prime_nodes[0] = slope;
        cache.beta_prime_nodes[1] = slope;
    } else if (num_steps >= 2) {
        const double inv_two_h = 1.0 / (2.0 * h);
        cache.beta_prime_nodes[0] = (-3.0 * cache.beta_nodes[0] + 4.0 * cache.beta_nodes[1] - cache.beta_nodes[2]) * inv_two_h;
        for (int k = 1; k < num_steps; ++k) {
            cache.beta_prime_nodes[k] = (cache.beta_nodes[k + 1] - cache.beta_nodes[k - 1]) * inv_two_h;
        }
        cache.beta_prime_nodes[num_steps] = (3.0 * cache.beta_nodes[num_steps] - 4.0 * cache.beta_nodes[num_steps - 1] + cache.beta_nodes[num_steps - 2]) * inv_two_h;
    }

    cache.beta_prime_midpoints.resize(num_steps);
    for (int k = 0; k < num_steps; ++k) {
        cache.beta_prime_midpoints[k] = (cache.beta_nodes[k + 1] - cache.beta_nodes[k]) / h;
    }

    return cache;
}

// Indefinite integral of P(s)/sqrt(A-s) ds, where P(s) is a quadratic.
// P(s) = c2*s^2 + c1*s + c0
// Returns the evaluated definite integral from s=0 to s=2.
double eval_quad_sqrt_integral(double A, double c2, double c1, double c0) {
    auto term_at = [A](double s, double c2, double c1, double c0) {
        if (A <= s) return 0.0;
        double sqrt_A_minus_s = std::sqrt(A - s);
		// antiderivatives for s^2, s^1, s^0 terms (pre-derived)
        double term_s2 = c2 * (-2.0/15.0) * sqrt_A_minus_s * (3.0*s*s + 4.0*A*s + 8.0*A*A);
        double term_s1 = c1 * (-2.0/3.0)  * sqrt_A_minus_s * (s + 2.0*A);
        double term_s0 = c0 * (-2.0)      * sqrt_A_minus_s;
        return term_s2 + term_s1 + term_s0;
    };
    return term_at(2.0, c2, c1, c0) - term_at(0.0, c2, c1, c0);
}

// Analytical weight coefficients from Linz (1985) / Paper page 13.
// dt is the time difference (t_n - t_start_of_block)
// alpha = (z/2) integral (1-s)(2-s)/sqrt (A - s) ds,  beta = z integral s(2-s)/sqrt (A - s) ds,  gamma = (z/2) integral s(s-1)/sqrt (A - s) ds
// Stieltjes-trapezoid history sum for node n + helper functions
struct Panel { double a, mid, b, d; };
Panel make_panel(int i, double h) {
    const double a = i*h, b = (i+1)*h;
    return {a, 0.5*(a+b), b, b-a};
}

struct StieltjesMoments { double J0, J1, J2; };
StieltjesMoments stieltjes_moments(double tn, double a, double b) {
    const double ra = sqrt_pos(tn - a), rb = sqrt_pos(tn - b);
    const double J0 = 2.0*(ra - rb);
    const double A  = tn - a;
    const double ra2 = ra*ra, rb2 = rb*rb, ra3 = ra2*ra, rb3 = rb2*rb, ra5 = ra3*ra2, rb5 = rb3*rb2;
    const double I1 = 2.0 * ((A*ra - ra3/3.0) - (A*rb - rb3/3.0));
    const double I2 = 2.0 * ((A*A*ra - (2.0/3.0)*A*ra3 + (1.0/5.0)*ra5)
                           - (A*A*rb - (2.0/3.0)*A*rb3 + (1.0/5.0)*rb5));
    const double J1 = a*J0 + I1;
    const double J2 = a*a*J0 + 2.0*a*I1 + I2;
    return {J0,J1,J2};
}

struct LagrangeWeights { double w0, w1, w2; };
LagrangeWeights stieltjes_lagrange_weights(const Panel& p, const StieltjesMoments& M) {
    const double d2 = p.d*p.d;
    const double s0=p.a, s1=p.mid, s2=p.b;
    const double w0 = ( 2.0/d2)*(M.J2 - (s1+s2)*M.J1 + s1*s2*M.J0);
    const double w1 = (-4.0/d2)*(M.J2 - (s0+s2)*M.J1 + s0*s2*M.J0);
    const double w2 = ( 2.0/d2)*(M.J2 - (s0+s1)*M.J1 + s0*s1*M.J0);
    return {w0,w1,w2};
}
// Coefficients for quadratic product integration on [0,dt] with singularity at dt
void block_coeffs(double dt, double z_eff, double& a, double& b, double& g) {
    if (dt<=0) { a=b=g=0; return; }
    const double A = dt / z_eff;
    a = 0.5*std::sqrt(z_eff) * eval_quad_sqrt_integral(A,  1.0,-3.0, 2.0); // (1-s)(2-s)
    b =       std::sqrt(z_eff) * eval_quad_sqrt_integral(A, -1.0, 2.0, 0.0); // s(2-s)
    g = 0.5*std::sqrt(z_eff) * eval_quad_sqrt_integral(A,  1.0,-1.0, 0.0); // s(s-1)
}


// Stieltjes-trapezoid history sum for node n (exclude the last interval)
double stieltjes_trap_history(
    int n,
    int last_panel_max,
    const std::vector<double>& F,
    double h,
    const RD_Params& pars,
    const KernelFn& K)
{
    if (n <= 1) return 0.0;
    const int Nf = static_cast<int>(F.size());
    if (Nf <= 1) return 0.0;

    const double tn = n * h;
    const int max_known = std::min(Nf - 1, n - 1);
    int last_panel = std::min(last_panel_max, max_known - 1);
    if (last_panel < 0) return 0.0;

    double S = 0.0;

    for (int i = 0; i <= last_panel; ++i) {
        const Panel p = make_panel(i, h);
        if (!(p.d > 0.0)) continue;

        const double Fa = F[i];
        const double Fb = F[i + 1];

        bool have_mid = false;
        double Fm = 0.5 * (Fa + Fb);
        if (i > 0 && (i + 1) <= max_known) {
            Fm = (-0.125) * F[i - 1] + 0.75 * F[i] + 0.375 * F[i + 1];
            have_mid = true;
        } else if (i == 0 && (i + 2) <= max_known) {
            Fm = 0.375 * F[i] + 0.75 * F[i + 1] - 0.125 * F[i + 2];
            have_mid = true;
        }

        const StieltjesMoments M = stieltjes_moments(tn, p.a, p.b);
        if (have_mid) {
            const LagrangeWeights W = stieltjes_lagrange_weights(p, M);

            const double Qa = K(tn, p.a, pars) * Fa;
            const double Qm = K(tn, p.mid, pars) * Fm;
            const double Qb = K(tn, p.b, pars) * Fb;

            if (std::isfinite(W.w0) && std::isfinite(W.w1) && std::isfinite(W.w2)) {
                S += W.w0 * Qa + W.w1 * Qm + W.w2 * Qb;
                continue;
            }
        }

        const double Wa = (p.b * M.J0 - M.J1) / p.d;
        const double Wb = (M.J1 - p.a * M.J0) / p.d;

        const double Qa = K(tn, p.a, pars) * Fa;
        const double Qb = K(tn, p.b, pars) * Fb;

        S += Wa * Qa + Wb * Qb;
    }
    return S;
}

// Seed first nodes of nu^f on a theta-grid using the Abel approximation.
// Returns the number of nodes filled (>=1).
struct QuadraticPanel {
    double a{0.0};
    double mid{0.0};
    double b{0.0};
    double Fa{0.0};
    double Fm{0.0};
    double Fb{0.0};
    // Metadata for controlled coarsening
    int chunk_id{-1};
    double base_w{0.0};

    double width() const { return b - a; }
};

inline double lagrange_quadratic_value(double x,
                                       double x0, double y0,
                                       double x1, double y1,
                                       double x2, double y2) {
    const double denom0 = (x0 - x1) * (x0 - x2);
    const double denom1 = (x1 - x0) * (x1 - x2);
    const double denom2 = (x2 - x0) * (x2 - x1);

    if (std::abs(denom0) < FPM_EPSILON ||
        std::abs(denom1) < FPM_EPSILON ||
        std::abs(denom2) < FPM_EPSILON) {
        return 0.5 * (y0 + y2);
    }

    const double l0 = ((x - x1) * (x - x2)) / denom0;
    const double l1 = ((x - x0) * (x - x2)) / denom1;
    const double l2 = ((x - x0) * (x - x1)) / denom2;
    return l0 * y0 + l1 * y1 + l2 * y2;
}

inline QuadraticPanel merge_quadratic_panels(const QuadraticPanel& left,
                                             const QuadraticPanel& right) {
    QuadraticPanel merged;
    merged.a = left.a;
    merged.b = right.b;
    merged.mid = 0.5 * (merged.a + merged.b);
    merged.Fa = left.Fa;
    merged.Fb = right.Fb;
    merged.chunk_id = left.chunk_id;
    merged.base_w = (left.base_w > 0.0 ? left.base_w : right.base_w);

    const double shared_node = left.b; // equals right.a
    const double shared_value = 0.5 * (left.Fb + right.Fa);
    merged.Fm = lagrange_quadratic_value(
        merged.mid,
        left.a, left.Fa,
        shared_node, shared_value,
        right.b, right.Fb);
    return merged;
}

inline double compute_midpoint_value(int left_idx,
                                     int right_idx,
                                     const std::vector<double>& grid,
                                     const std::vector<double>& F) {
    const int N = static_cast<int>(grid.size());
    if (left_idx < 0 || right_idx >= N || left_idx >= right_idx) {
        return 0.0;
    }

    const double mid = 0.5 * (grid[left_idx] + grid[right_idx]);

    int aux = (left_idx > 0)
                  ? (left_idx - 1)
                  : ((right_idx + 1 < N) ? (right_idx + 1) : left_idx);
    if (aux == left_idx || aux == right_idx) {
        if (right_idx + 1 < N) {
            aux = right_idx + 1;
        } else if (left_idx > 0) {
            aux = left_idx - 1;
        }
    }

    if (aux == left_idx || aux == right_idx) {
        return 0.5 * (F[left_idx] + F[right_idx]);
    }

    return lagrange_quadratic_value(
        mid,
        grid[aux], F[aux],
        grid[left_idx], F[left_idx],
        grid[right_idx], F[right_idx]);
}

inline QuadraticPanel panel_from_indices(int left_idx,
                                         int right_idx,
                                         const std::vector<double>& grid,
                                         const std::vector<double>& F) {
    QuadraticPanel panel;
    if (left_idx < 0 || right_idx >= static_cast<int>(grid.size()) || left_idx >= right_idx) {
        return panel;
    }

    panel.a = grid[left_idx];
    panel.b = grid[right_idx];
    panel.mid = 0.5 * (panel.a + panel.b);
    panel.Fa = F[left_idx];
    panel.Fb = F[right_idx];
    panel.Fm = compute_midpoint_value(left_idx, right_idx, grid, F);
    panel.base_w = panel.b - panel.a;
    return panel;
}

inline QuadraticPanel panel_from_indices_chunk(int left_idx,
                                               int right_idx,
                                               const std::vector<double>& grid,
                                               const std::vector<double>& F,
                                               int chunk_id,
                                               double chunk_h_hint = 0.0) {
    QuadraticPanel p = panel_from_indices(left_idx, right_idx, grid, F);
    p.chunk_id = chunk_id;
    if (chunk_h_hint > 0.0) p.base_w = chunk_h_hint;
    return p;
}

inline double stieltjes_panel_quadratic(const QuadraticPanel& panel,
                                        double tn,
                                        const RD_Params& pars,
                                        const KernelFn& K) {
    if (!(panel.b > panel.a) || tn <= panel.a + FPM_EPSILON) {
        return 0.0;
    }
    const Panel p{panel.a, panel.mid, panel.b, panel.b - panel.a};
    const StieltjesMoments M = stieltjes_moments(tn, p.a, p.b);
    const LagrangeWeights W = stieltjes_lagrange_weights(p, M);

    const double Qa = K(tn, panel.a, pars) * panel.Fa;
    const double Qm = K(tn, panel.mid, pars) * panel.Fm;
    const double Qb = K(tn, panel.b, pars) * panel.Fb;

    return W.w0 * Qa + W.w1 * Qm + W.w2 * Qb;
}

// Near the evaluation point tn, large coarsened panels can introduce a small bias
// in the history sum. To stabilize the endpoint (e.g., at t_max), refine any
// panel whose right edge is too close to tn by splitting it once into two
// equal subpanels and applying the same quadratic product-integration on each.
inline double stieltjes_panel_quadratic_refined(const QuadraticPanel& panel,
                                                double tn,
                                                const RD_Params& pars,
                                                const KernelFn& K,
                                                double refine_factor = 2.0) {
    if (!(panel.b > panel.a)) return 0.0;
    const double w = panel.width();
    if (tn <= panel.a + FPM_EPSILON) return 0.0;

    // If far from tn relative to panel width, the single-panel evaluation is fine.
    if ((tn - panel.b) >= refine_factor * w) {
        return stieltjes_panel_quadratic(panel, tn, pars, K);
    }

    // Split [a,b] into [a,mid] and [mid,b], reconstruct midpoints by the
    // quadratic defined by (a,Fa), (mid,Fm), (b,Fb).
    QuadraticPanel left, right;
    left.a = panel.a; left.b = panel.mid; left.mid = 0.5 * (left.a + left.b);
    right.a = panel.mid; right.b = panel.b; right.mid = 0.5 * (right.a + right.b);

    left.Fa = panel.Fa; left.Fb = panel.Fm;
    right.Fa = panel.Fm; right.Fb = panel.Fb;

    auto qval = [&](double x){
        return lagrange_quadratic_value(x,
                                        panel.a,  panel.Fa,
                                        panel.mid, panel.Fm,
                                        panel.b,  panel.Fb);
    };
    left.Fm = qval(left.mid);
    right.Fm = qval(right.mid);

    return stieltjes_panel_quadratic(left, tn, pars, K)
         + stieltjes_panel_quadratic(right, tn, pars, K);
}
int seed_nu_on_grid(
    const std::vector<double>& grid,
    const RD_Params& pars,
    std::vector<double>& F,
    const AbelFn& abel,
    double t_cut = SMALL_T_SCALED_THRESHOLD)
{
    const int N = static_cast<int>(grid.size());
    int J = 0;
    for (; J < N; ++J) {
        const double v = grid[J];
        if (v > t_cut) break;
        F[J] = abel(v, pars);
    }

    if (N >= 2 && J < 2) {
        F[1] = abel(grid[1], pars);
        J = std::max(J, 2);
    }

    if (J < 1) J = 1;
    if (J & 1) --J;
    return J;
}

struct SeedInfo {
    int first_unseeded = 0;
    int last_seeded = -1;
};

inline SeedInfo seed_nu_on_nonuniform_grid(const std::vector<double>& grid,
                                           const RD_Params& pars,
                                           std::vector<double>& F,
                                           const AbelFn& abel,
                                           double t_cut = SMALL_T_SCALED_THRESHOLD) {
    SeedInfo info{};
    const int N = static_cast<int>(grid.size());
    if (N == 0) {
        info.first_unseeded = 0;
        info.last_seeded = -1;
        return info;
    }

    int last = -1;
    for (int j = 0; j < N; ++j) {
        const double t = grid[j];
        if (t > t_cut) break;
        F[j] = abel(t, pars);
        last = j;
    }

    if (last < 0) {
        F[0] = abel(grid[0], pars);
        last = 0;
    }
    if (last < 1 && N >= 2) {
        F[1] = abel(grid[1], pars);
        last = std::max(last, 1);
    }

    info.last_seeded = std::min(N - 1, std::max(0, last));

    int first_unseeded = info.last_seeded + 1;
    if (first_unseeded < 1) {
        first_unseeded = std::min(1, N - 1);
    }
    if (first_unseeded % 2 != 0) {
        ++first_unseeded;
    }
    if (first_unseeded >= N) {
        first_unseeded = N - 1;
        if (first_unseeded % 2 != 0 && first_unseeded > 0) {
            --first_unseeded;
        }
    }
    info.first_unseeded = std::max(0, first_unseeded);
    return info;
}
// Quadratic block step on a uniform grid t_j = j*h, j = 0..N, N even
std::vector<double> solve_nu_block_by_block_impl(
    double t_max,
    const RD_Params& pars,
    int num_steps,
    std::vector<double>& grid,
    std::vector<double> F,
    int k0_seeded,
    const KernelFn& K,
    const ForcingFn& G)
{
    if (num_steps < 2) num_steps = 2;
    if (num_steps % 2) ++num_steps;
    const int N = num_steps;
    const double h = t_max / N;

    grid.resize(N + 1);
    for (int i = 0; i <= N; ++i) {
        grid[i] = i * h;
    }

    if (static_cast<int>(F.size()) < N + 1) {
        F.resize(N + 1, 0.0);
    }

    if (std::abs(F[0]) < 1e-300) {
        F[0] = G(grid[0], pars);
    }

    const int M = N / 2;
    const int first_unseeded_block = std::max(0, k0_seeded / 2);
    const int m0 = std::min(M, first_unseeded_block);
    for (int m = m0; m < M; ++m) {
        const int j0 = 2 * m;
        const int j1 = 2 * m + 1;
        const int j2 = 2 * m + 2;

        const double t0 = j0 * h;
        const double t1 = j1 * h;
        const double t2 = j2 * h;
        const double tmid = t0 + 0.5 * h;

        const double S1 = stieltjes_trap_history(j1, j0 - 1, F, h, pars, K);
        const double S2 = stieltjes_trap_history(j2, j0 - 1, F, h, pars, K);

        double a1, b1, g1, a2, b2, g2;
        block_coeffs(h, 0.5 * h, a1, b1, g1);
        block_coeffs(2 * h, h, a2, b2, g2);

        const double K1_0 = K(t1, t0, pars);
        const double K1_mid = K(t1, tmid, pars);
        const double K1_1 = K(t1, t1, pars);

        const double K2_0 = K(t2, t0, pars);
        const double K2_1 = K(t2, t1, pars);
        const double K2_2 = K(t2, t2, pars);

        const double G1 = G(t1, pars);
        const double G2 = G(t2, pars);

        double F1_pred;
        if (j0 >= 1) {
            const double den1_pred = 1.0 - g1 * K1_1 - (3.0 / 8.0) * b1 * K1_mid;
            const double rhs1_pred = G1 + S1
                + a1 * K1_0 * F[j0]
                + b1 * K1_mid * ((-0.125) * F[j0 - 1] + 0.75 * F[j0]);
            F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
        } else {
            const double den1_pred = 1.0 - g1 * K1_1 - 0.5 * b1 * K1_mid;
            const double rhs1_pred = G1 + S1 + a1 * K1_0 * F[j0] + 0.5 * b1 * K1_mid * F[j0];
            F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
        }

        const double den2_pred = 1.0 - g2 * K2_2;
        const double rhs2_pred = G2 + S2 + a2 * K2_0 * F[j0] + b2 * K2_1 * F1_pred;
        const double F2_pred = rhs2_pred / ((std::abs(den2_pred) < FPM_EPSILON) ? FPM_EPSILON : den2_pred);

        const double A = 1.0 - (3.0 / 4.0) * b1 * K1_mid - g1 * K1_1;
        const double B = (1.0 / 8.0) * b1 * K1_mid;
        const double C = G1 + S1 + a1 * K1_0 * F[j0]
                        + (3.0 / 8.0) * b1 * K1_mid * F[j0];

        const double D = -b2 * K2_1;
        const double E = 1.0 - g2 * K2_2;
        const double R = G2 + S2 + a2 * K2_0 * F[j0];

        const double det = A * E - B * D;
        const double inv_det = (std::abs(det) < FPM_EPSILON) ? (1.0 / FPM_EPSILON) : (1.0 / det);

        const double F1_corr = (E * C - B * R) * inv_det;
        const double F2_corr = (-D * C + A * R) * inv_det;

        const double blend = 1.0;
        F[j1] = blend * F1_corr + (1.0 - blend) * F1_pred;
        F[j2] = blend * F2_corr + (1.0 - blend) * F2_pred;
    }
    return F;
}

// Dense non-uniform history sum: sum quadratic product-integration over all
// native panels strictly before tn, up to (but excluding) panel index cut_index.
// This restores the exact O(N^2) history on a non-uniform grid.
inline double stieltjes_history_dense_nonuniform(
    int cut_index_exclusive,
    double tn,
    const std::vector<double>& grid,
    const std::vector<double>& F,
    const RD_Params& pars,
    const KernelFn& K)
{
    const int total_nodes = static_cast<int>(grid.size());
    if (total_nodes < 2) return 0.0;
    int last_panel = std::min(std::max(cut_index_exclusive - 1, 0), total_nodes - 2);
    if (last_panel < 0) return 0.0;

    double S = 0.0;
    for (int i = 0; i <= last_panel; ++i) {
        // Panels are strictly before tn
        if (grid[i + 1] >= tn - FPM_EPSILON) break;
        QuadraticPanel p = panel_from_indices(i, i + 1, grid, F);
        if (p.b <= p.a) continue;
        S += stieltjes_panel_quadratic(p, tn, pars, K);
    }
    return S;
}

std::vector<double> solve_nu_block_with_abel(
    double t_max,
    const RD_Params& pars,
    int num_steps,
    std::vector<double>& grid,
    const KernelFn& K,
    const ForcingFn& G,
    const AbelFn& abel,
    double t_cut = SMALL_T_SCALED_THRESHOLD)
{
    std::vector<double> F(num_steps + 1, 0.0);
    grid.resize(num_steps + 1);
    for (int i = 0; i <= num_steps; ++i) {
        grid[i] = (t_max / num_steps) * i;
    }

    const int k0 = seed_nu_on_grid(grid, pars, F, abel, t_cut);
    return solve_nu_block_by_block_impl(t_max, pars, num_steps, grid, F, k0, K, G);
}

struct ChunkSpec {
    double start = 0.0;
    double end = 0.0;
    double h = 0.0;
    int num_steps = 0;
    int start_index = 0;

    int end_index() const { return start_index + num_steps; }
    double span() const { return end - start; }
};

struct ChunkingOptions {
    // Geometric growth ratio of step size: h_{m+1} = ratio * h_m.
    // Larger => coarser later chunks (fewer total panels). Smaller (>=1.5) => finer.
    double ratio = 2.0;
    // Target number of panels in the first chunk near zero (even is enforced).
    // Larger => finer resolution near zero.
    int base_panels = 128;
    // Maximum number of chunks to create. Larger => more chunk transitions and
    // potentially more total panels (if remaining span forces extra chunks).
    int max_chunks = 12;
    // Number of most-recent panels kept at native resolution in history.
    // Larger => more exact history near tn (slower, more memory), smaller => faster.
    int history_recent = 64;
    // Tolerance for treating adjacent panels as contiguous when merging tiers.
    // Larger => more aggressive merging.
    double width_tolerance = 1e-8;
    // Near-tn refinement threshold used during history evaluation:
    // refine if (tn - panel.b) < refine_ratio * panel.width(). Larger => less refinement.
    double refine_ratio = 8.0;
    // Minimum width factor to stop recursive splitting by width.
    // We stop splitting when panel.width() <= refine_min_width_factor * base_width.
    // Larger => less splitting (coarser); smaller (>=1.0) => more splitting (finer).
    double refine_min_width_factor = 1.01;
    // Disable history compression entirely. When true, keep all panels exact in
    // the recent window. Increases cost toward O(N^2) but maximizes accuracy.
    bool disable_compression = false;
    // Keep the last K chunks exact in the history (i.e., do not compress any
    // panels that belong to these last chunks). This links the history dropoff
    // to chunk boundaries instead of a raw panel count and typically improves
    // stability near tn. Larger K => more exact history, slower.
    int keep_last_chunks = 2;
    // Additional accuracy guards for compressed history evaluation:
    // If the kernel varies too much across a coarsened panel at evaluation t,
    // force a split regardless of distance. Relative curvature tolerance using
    // three-point second-difference of K(t, s) across [a, mid, b]. Smaller =>
    // more splitting (more accurate, slower).
    double kernel_curvature_tol = 0.02;
    // Also consider first-difference (relative) variation of K across a panel.
    // If max(|K(a)-K(mid)|, |K(b)-K(mid)|)/|K(mid)| exceeds this tol, split.
    double kernel_variation_tol = 0.05;
    // Hard cap on how wide a coarsened panel may be, as a multiple of the
    // base (uniform) width. Panels wider than max_coarsen_factor * base_width
    // are recursively split until under the cap, even if far from t.
    double max_coarsen_factor = 16.0;
    // Print per-chunk/tier history statistics for diagnostics
    bool print_history_stats = true;
};

inline std::vector<ChunkSpec> build_chunked_theta_layout(
    double theta_max,
    int uniform_num_steps,
    const ChunkingOptions& opts)
{
    std::vector<ChunkSpec> chunks;
    if (!(theta_max > 0.0) || uniform_num_steps <= 0) {
        return chunks;
    }

    if (uniform_num_steps % 2 != 0) {
        ++uniform_num_steps;
    }

    // Allow near-uniform chunking by permitting ratio >= 1.0. Values close to 1.0
    // mean little to no coarsening across chunks.
    const double ratio = std::max(1.0, opts.ratio);
    int panels_per_chunk = std::max(4, opts.base_panels);
    if (panels_per_chunk % 2 != 0) {
        ++panels_per_chunk;
    }
    const int max_chunks = std::max(1, opts.max_chunks);
    const double base_h = theta_max / uniform_num_steps;

    double start = 0.0;
    int start_index = 0;
    double current_h = base_h;

    for (int chunk_id = 0; chunk_id < max_chunks && start < theta_max - 1e-12; ++chunk_id) {
        int panels = panels_per_chunk;
        const double span = panels * current_h;
        double end = start + span;

        const bool is_last = (chunk_id == max_chunks - 1) || (end >= theta_max - 1e-12);
        if (is_last) {
            const double remaining = theta_max - start;
            if (remaining <= 0.0) {
                break;
            }
            panels = std::max(2, static_cast<int>(std::ceil(remaining / current_h)));
            if (panels % 2 != 0) {
                ++panels;
            }
            const double adjusted_h = remaining / panels;

            ChunkSpec chunk;
            chunk.start = start;
            chunk.end = theta_max;
            chunk.h = adjusted_h;
            chunk.num_steps = panels;
            chunk.start_index = start_index;
            chunks.push_back(chunk);
            start_index += panels;
            start = theta_max;
            break;
        }

        ChunkSpec chunk;
        chunk.start = start;
        chunk.end = end;
        chunk.h = current_h;
        chunk.num_steps = panels;
        chunk.start_index = start_index;
        chunks.push_back(chunk);

        start = end;
        start_index += panels;
        current_h *= ratio;
    }

    if (start < theta_max - 1e-12) {
        const double remaining = theta_max - start;
        double h_guess = chunks.empty() ? base_h : chunks.back().h * opts.ratio;
        if (!(h_guess > 0.0)) {
            h_guess = base_h;
        }
        int panels = std::max(2, static_cast<int>(std::ceil(remaining / h_guess)));
        if (panels % 2 != 0) {
            ++panels;
        }
        const double adjusted_h = remaining / panels;

        ChunkSpec chunk;
        chunk.start = start;
        chunk.end = theta_max;
        chunk.h = adjusted_h;
        chunk.num_steps = panels;
        chunk.start_index = start_index;
        chunks.push_back(chunk);
    }

    return chunks;
}

// Builds the full theta grid that the chunked solver will use, without running the solver.
// Useful for preparing boundary/spline data (e.g., beta and beta') that must be known
// before calling the kernel solve.
inline std::vector<double> build_chunked_theta_grid(
    double theta_max,
    int uniform_num_steps,
    const ChunkingOptions& opts)
{
    std::vector<double> grid;
    if (!(theta_max > 0.0) || uniform_num_steps <= 0) {
        grid.assign(1, 0.0);
        return grid;
    }

    if (uniform_num_steps % 2 != 0) {
        ++uniform_num_steps;
    }

    const auto chunks = build_chunked_theta_layout(theta_max, uniform_num_steps, opts);
    const int total_panels = std::accumulate(chunks.begin(), chunks.end(), 0,
        [](int acc, const ChunkSpec& c){ return acc + std::max(0, c.num_steps); });

    grid.reserve(std::max(1, total_panels + 1));
    grid.clear();
    double current = 0.0;
    grid.push_back(current);
    for (const auto& chunk : chunks) {
        current = grid.back();
        for (int j = 0; j < chunk.num_steps; ++j) {
            current += chunk.h;
            if (j == chunk.num_steps - 1) {
                current = chunk.end;
            }
            grid.push_back(current);
        }
    }
    if (!grid.empty()) {
        grid.back() = theta_max;
    }
    return grid;
}

class TieredQuadraticHistory {
public:
    explicit TieredQuadraticHistory(double base_width,
                                    int recent_capacity,
                                    int max_tiers,
                                    double width_tol,
                                    double refine_ratio = 8.0,
                                    double refine_min_width_factor = 1.01,
                                    double kernel_curvature_tol = 0.02,
                                    double kernel_variation_tol = 0.05,
                                    double max_coarsen_factor = 16.0)
        : base_width_(std::max(base_width, FPM_EPSILON)),
          recent_capacity_(std::max(1, recent_capacity)),
          max_tiers_(std::max(1, max_tiers)),
          width_tol_(std::max(width_tol, 1e-12)),
          refine_ratio_(std::max(1.0, refine_ratio)),
          refine_min_width_factor_(std::max(1.0, refine_min_width_factor)),
          kernel_curvature_tol_(std::max(0.0, kernel_curvature_tol)),
          kernel_variation_tol_(std::max(0.0, kernel_variation_tol)),
          max_coarsen_factor_(std::max(2.0, max_coarsen_factor)),
          recent_(),
          tiers_(std::max(1, max_tiers)) {}

    void clear(double base_width) {
        base_width_ = std::max(base_width, FPM_EPSILON);
        recent_.clear();
        for (auto& v : tiers_) v.clear();
    }

    void push_panel(const QuadraticPanel& panel) {
        if (!(panel.b > panel.a)) return;
        // Keep the most recent L panels at native resolution
        recent_.push_back(panel);
        while (static_cast<int>(recent_.size()) > recent_capacity_) {
            QuadraticPanel p_old = recent_.front();
            recent_.pop_front();
            add_to_tier(0, p_old);
        }
    }

    double evaluate(double t, const RD_Params& pars, const KernelFn& K) const {
        double total = 0.0;
        auto accumulate_panel = [&](const QuadraticPanel& p) {
            if (p.b >= t - FPM_EPSILON) return;
            total += eval_panel_adaptive(p, t, pars, K);
        };
        for (const auto& p : recent_) accumulate_panel(p);
        for (int lvl = 0; lvl < max_tiers_; ++lvl) {
            const auto& v = tiers_[lvl];
            for (const auto& p : v) accumulate_panel(p);
        }
        return total;
    }

    // Evaluate history contributions at evaluation time t, but only from panels
    // whose right edge is strictly less than cut_b. This mirrors the legacy
    // uniform-grid logic where, for the second node in a block, the last two
    // panels are excluded from the history and handled analytically by the
    // block update. If a compressed panel straddles cut_b, we clip it to
    // [a, cut_b] by reconstructing values from the panel's quadratic and then
    // evaluate adaptively on the clipped sub-panel.
    double evaluate_with_cut(double t, double cut_b,
                             const RD_Params& pars, const KernelFn& K) const {
        double total = 0.0;

        auto qval_from_panel = [](const QuadraticPanel& panel, double x) {
            return lagrange_quadratic_value(x,
                                            panel.a,  panel.Fa,
                                            panel.mid, panel.Fm,
                                            panel.b,  panel.Fb);
        };

        auto accumulate_with_cut = [&](const QuadraticPanel& p) {
            // Nothing to add if panel starts at/after cut or panel ends at/after t
            if (p.a >= cut_b - FPM_EPSILON) return;
            if (p.b >= t - FPM_EPSILON) return; // exclude the last panel by t

            if (p.b <= cut_b - FPM_EPSILON) {
                // Entire panel is strictly before cut_b: take full contribution
                total += eval_panel_adaptive(p, t, pars, K);
                return;
            }

            // Panel crosses the cut: clip to [a, cut_b]
            QuadraticPanel clipped;
            clipped.a = p.a;
            clipped.b = cut_b;
            clipped.mid = 0.5 * (clipped.a + clipped.b);
            clipped.Fa = p.Fa;
            clipped.Fb = qval_from_panel(p, clipped.b);
            clipped.Fm = qval_from_panel(p, clipped.mid);
            clipped.chunk_id = p.chunk_id;
            clipped.base_w = (p.base_w > 0.0) ? p.base_w : base_width_;
            if (clipped.b > clipped.a + FPM_EPSILON) {
                total += eval_panel_adaptive(clipped, t, pars, K);
            }
        };

        for (const auto& p : recent_) accumulate_with_cut(p);
        for (int lvl = 0; lvl < max_tiers_; ++lvl) {
            const auto& v = tiers_[lvl];
            for (const auto& p : v) accumulate_with_cut(p);
        }
        return total;
    }

    void dump_diagnostics_per_chunk(const std::vector<ChunkSpec>& chunks) const {
        const int C = static_cast<int>(chunks.size());
        std::vector<int> recent_counts(C, 0);
        int recent_unknown = 0;
        for (const auto& p : recent_) {
            if (p.chunk_id >= 0 && p.chunk_id < C) recent_counts[p.chunk_id]++;
            else recent_unknown++;
        }
        auto print_counts = [&](const char* label, const std::vector<int>& cnt, int unknown, int total){
            Rcout << label << ": total=" << total << ", per-chunk=[";
            for (int i = 0; i < C; ++i) {
                if (i) Rcout << ", ";
                Rcout << cnt[i];
            }
            Rcout << "]";
            if (unknown) Rcout << ", unknown=" << unknown;
            Rcout << std::endl;
        };
        int recent_total = static_cast<int>(recent_.size());
        Rcout << "history stats:" << std::endl;
        print_counts("  recent", recent_counts, recent_unknown, recent_total);

        // Tiers
        for (int lvl = 0; lvl < max_tiers_; ++lvl) {
            const auto& v = tiers_[lvl];
            std::vector<int> counts(C, 0);
            int unknown = 0;
            double sum_ratio = 0.0;
            int ratio_n = 0;
            for (const auto& p : v) {
                if (p.chunk_id >= 0 && p.chunk_id < C) counts[p.chunk_id]++;
                else unknown++;
                const double bw = (p.base_w > 0.0) ? p.base_w : base_width_;
                if (bw > 0.0) { sum_ratio += (p.width() / bw); ratio_n++; }
            }
            const int total = static_cast<int>(v.size());
            Rcout << "  tier " << lvl << ": total=" << total
                  << ", mean_width_factor=" << (ratio_n ? (sum_ratio / ratio_n) : 0.0)
                  << ", per-chunk=[";
            for (int i = 0; i < C; ++i) {
                if (i) Rcout << ", ";
                Rcout << counts[i];
            }
            Rcout << "]";
            if (unknown) Rcout << ", unknown=" << unknown;
            Rcout << std::endl;
        }
    }

private:
    static inline double max3(double a, double b, double c) {
        return std::max(a, std::max(b, c));
    }

    bool are_contiguous(const QuadraticPanel& left, const QuadraticPanel& right) const {
        // Accept exact continuity or small numerical gaps/overlaps
        if (left.chunk_id != right.chunk_id) return false;
        const double scale = max3(1.0, left.width(), right.width());
        return std::abs(left.b - right.a) <= width_tol_ * scale;
    }

    // Recursively refine panels that are too close to t relative to their
    // width, to control the error from coarsened tiers near the endpoint.
    double eval_panel_adaptive(const QuadraticPanel& panel,
                               double t,
                               const RD_Params& pars,
                               const KernelFn& K) const {
        if (!(panel.b > panel.a) || t <= panel.a + FPM_EPSILON) return 0.0;

        const double w = panel.width();
        const double dist = t - panel.b; // distance from panel's right edge to t

        // Heuristic: decide if we can safely evaluate without further splits.
        // 1) never allow panels wider than a hard multiple of the base width
        // 2) if the kernel varies too much across [a,mid,b] at time t, split
        // 3) otherwise, if sufficiently far (dist >= R*w) or already fine by width,
        //    accept a single evaluation.
        const double R = refine_ratio_;

        // Hard width cap check (per-panel base width if available)
        const double local_base = (panel.base_w > 0.0) ? panel.base_w : base_width_;
        const bool too_wide = (w > max_coarsen_factor_ * local_base);

        // Kernel curvature (relative) across panel
        double Ka = 0.0, Km = 0.0, Kb = 0.0;
        Ka = K(t, panel.a, pars);
        Km = K(t, panel.mid, pars);
        Kb = K(t, panel.b, pars);
        const double denom = std::abs(Km) + 1e-14;
        const double k_curv = std::abs(Ka - 2.0 * Km + Kb) / denom;
        const double k_var = std::max(std::abs(Ka - Km), std::abs(Kb - Km)) / denom;

        const bool wide_ok = !too_wide;
        const bool kernel_ok = (k_curv <= kernel_curvature_tol_) && (k_var <= kernel_variation_tol_);
        const bool far_enough = (dist >= R * w);
        const bool fine_enough = (w <= refine_min_width_factor_ * local_base);

        if (wide_ok && kernel_ok && (far_enough || fine_enough)) {
            return stieltjes_panel_quadratic(panel, t, pars, K);
        }

        // Split panel into two halves using the quadratic defined by
        // (a,Fa), (mid,Fm), (b,Fb) to reconstruct child midpoints.
        QuadraticPanel left, right;
        left.a = panel.a; left.b = panel.mid; left.mid = 0.5 * (left.a + left.b);
        right.a = panel.mid; right.b = panel.b; right.mid = 0.5 * (right.a + right.b);

        left.Fa = panel.Fa; left.Fb = panel.Fm;
        right.Fa = panel.Fm; right.Fb = panel.Fb;
        left.chunk_id = panel.chunk_id; right.chunk_id = panel.chunk_id;
        left.base_w = local_base; right.base_w = local_base;

        auto qval = [&](double x){
            return lagrange_quadratic_value(x,
                                            panel.a,  panel.Fa,
                                            panel.mid, panel.Fm,
                                            panel.b,  panel.Fb);
        };
        left.Fm = qval(left.mid);
        right.Fm = qval(right.mid);

        return eval_panel_adaptive(left,  t, pars, K)
             + eval_panel_adaptive(right, t, pars, K);
    }

    void add_to_tier(int level, const QuadraticPanel& panel) {
        if (level >= max_tiers_) {
            // Append to last tier without merging if we ran out of tiers
            tiers_[max_tiers_ - 1].push_back(panel);
            return;
        }

        auto& v = tiers_[level];
        if (!v.empty() && are_contiguous(v.back(), panel)) {
            QuadraticPanel merged = merge_quadratic_panels(v.back(), panel);
            v.pop_back();
            add_to_tier(level + 1, merged);
        } else {
            v.push_back(panel);
        }
    }

    double base_width_;
    int recent_capacity_;
    int max_tiers_;
    double width_tol_;
    double refine_ratio_;
    double refine_min_width_factor_;
    double kernel_curvature_tol_;
    double kernel_variation_tol_;
    double max_coarsen_factor_;
    std::deque<QuadraticPanel> recent_;
    std::vector<std::vector<QuadraticPanel>> tiers_;
};

// Chunked theta-grid Volterra solver with tiered history compression.
// Retains the quadratic integration accuracy while reducing asymptotic cost.
std::vector<double> solve_nu_chunked_with_abel(
    double t_max,
    const RD_Params& pars,
    int uniform_num_steps,
    std::vector<double>& grid,
    const KernelFn& K,
    const ForcingFn& G,
    const AbelFn& abel,
    const ChunkingOptions& opts,
    double t_cut = SMALL_T_SCALED_THRESHOLD)
{    // Single-chunk fallback: use legacy uniform-grid solver to guarantee
    // numerical identity with the baseline when no chunking is present.
    // This avoids any differences from the generalized history path.
    {
        auto __chunks_probe = build_chunked_theta_layout(t_max, uniform_num_steps, opts);
        if (__chunks_probe.size() == 1) {
            // Delegate entirely to the legacy implementation
            return solve_nu_block_with_abel(t_max, pars, uniform_num_steps, grid, K, G, abel, t_cut);
        }
    }
    if (t_max <= 0.0) {
        grid.assign(1, 0.0);
        return std::vector<double>{G(0.0, pars)};
    }

    if (uniform_num_steps < 2) {
        uniform_num_steps = 2;
    }
    if (uniform_num_steps % 2 != 0) {
        ++uniform_num_steps;
    }

    auto chunks = build_chunked_theta_layout(t_max, uniform_num_steps, opts);
    if (chunks.empty()) {
        ChunkSpec chunk;
        chunk.start = 0.0;
        chunk.end = t_max;
        chunk.h = t_max / uniform_num_steps;
        chunk.num_steps = uniform_num_steps;
        chunk.start_index = 0;
        chunks.push_back(chunk);
    }

    // Build the theta grid using the same layout logic as all other components.
    // This single-sources the grid construction and avoids drift.
    grid = build_chunked_theta_grid(t_max, uniform_num_steps, opts);

    const int total_nodes = static_cast<int>(grid.size());

    std::vector<double> F(total_nodes, 0.0);
    F[0] = G(grid[0], pars);

    SeedInfo seed = seed_nu_on_nonuniform_grid(grid, pars, F, abel, t_cut);
    int first_unseeded = std::min(seed.first_unseeded, total_nodes - 1);
    int last_seeded = std::min(std::max(seed.last_seeded, 0), total_nodes - 1);

    if (first_unseeded >= total_nodes - 1) {
        return F;
    }

    const double base_h = t_max / static_cast<double>(uniform_num_steps);

    for (size_t ci = 0; ci < chunks.size(); ++ci) {
        const auto& chunk = chunks[ci];
        const int chunk_start = chunk.start_index;
        const int chunk_end = chunk.end_index();
        if (chunk.num_steps < 2) {
            continue;
        }

        int first_block = std::max(chunk_start / 2, first_unseeded / 2);
        const int block_end = chunk_end / 2;

        // Coefficients for the quadratic product integration depend only on h within a chunk.
        const double h = chunk.h;
        double a1, b1, g1, a2, b2, g2;
        block_coeffs(h, 0.5 * h, a1, b1, g1);
        block_coeffs(2.0 * h, h, a2, b2, g2);

        for (int m = first_block; m < block_end; ++m) {
            const int j0 = 2 * m;
            const int j1 = j0 + 1;
            const int j2 = j0 + 2;

            if (j2 >= total_nodes) {
                break;
            }
            if (j0 < chunk_start || j2 > chunk_end) {
                continue;
            }
            if (j1 <= last_seeded && j2 <= last_seeded) {
                continue;
            }

            const double t0 = grid[j0];
            const double t1 = grid[j1];
            const double t2 = grid[j2];
            const double tmid = t0 + 0.5 * h;

            // Dense non-uniform history (no compression):
            // Cut at t0 => exclude [j0,j1] and [j1,j2] from S2, and [j0,j1] from S1
            const double S1 = stieltjes_history_dense_nonuniform(j0, t1, grid, F, pars, K);
            const double S2 = stieltjes_history_dense_nonuniform(j0, t2, grid, F, pars, K);

            const double K1_0 = K(t1, t0, pars);
            const double K1_mid = K(t1, tmid, pars);
            const double K1_1 = K(t1, t1, pars);

            const double K2_0 = K(t2, t0, pars);
            const double K2_1 = K(t2, t1, pars);
            const double K2_2 = K(t2, t2, pars);

            const double G1 = G(t1, pars);
            const double G2 = G(t2, pars);

            double F1_pred = 0.0;
            const int local_j0 = j0 - chunk_start;
            if (local_j0 >= 1) {
                const double Fmid = (-0.125) * F[j0 - 1] + 0.75 * F[j0] + 0.375 * F[j0 + 1];
                const double den1_pred = 1.0 - g1 * K1_1 - (3.0 / 8.0) * b1 * K1_mid;
                const double rhs1_pred = G1 + S1
                    + a1 * K1_0 * F[j0]
                    + b1 * K1_mid * Fmid;
                F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
            } else {
                const double den1_pred = 1.0 - g1 * K1_1 - 0.5 * b1 * K1_mid;
                const double rhs1_pred = G1 + S1 + a1 * K1_0 * F[j0] + 0.5 * b1 * K1_mid * F[j0];
                F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
            }

            const double den2_pred = 1.0 - g2 * K2_2;
            const double rhs2_pred = G2 + S2 + a2 * K2_0 * F[j0] + b2 * K2_1 * F1_pred;
            const double F2_pred = rhs2_pred / ((std::abs(den2_pred) < FPM_EPSILON) ? FPM_EPSILON : den2_pred);

            const double A = 1.0 - (3.0 / 4.0) * b1 * K1_mid - g1 * K1_1;
            const double B = (1.0 / 8.0) * b1 * K1_mid;
            const double C = G1 + S1 + a1 * K1_0 * F[j0]
                + (3.0 / 8.0) * b1 * K1_mid * F[j0];

            const double D = -b2 * K2_1;
            const double E = 1.0 - g2 * K2_2;
            const double R = G2 + S2 + a2 * K2_0 * F[j0];

            const double det = A * E - B * D;
            const double inv_det = (std::abs(det) < FPM_EPSILON) ? (1.0 / FPM_EPSILON) : (1.0 / det);

            const double F1_corr = (E * C - B * R) * inv_det;
            const double F2_corr = (-D * C + A * R) * inv_det;

            const double blend = 1.0;
            F[j1] = blend * F1_corr + (1.0 - blend) * F1_pred;
            F[j2] = blend * F2_corr + (1.0 - blend) * F2_pred;

            // No history compression: dense history recomputed each step

            last_seeded = std::max(last_seeded, j2);
        }
    }

    return F;
}

inline std::vector<double> solve_nu_block_with_abel_chunked(
    double t_max,
    const RD_Params& pars,
    int num_steps,
    std::vector<double>& grid,
    const KernelFn& K,
    const ForcingFn& G,
    const AbelFn& abel,
    const ChunkingOptions& opts = ChunkingOptions(),
    double t_cut = SMALL_T_SCALED_THRESHOLD)
{
    return solve_nu_chunked_with_abel(t_max, pars, num_steps, grid, K, G, abel, opts, t_cut);
}

// Linear-product weights for the integral of (Linear(H) / (x-s)^{3/2}) ds on [a,b]
void weights_panel_32(double x, double a, double b, double& Wa, double& Wb) {
    const double ra = std::sqrt(std::max(0.0, x - a));
    const double rb = std::sqrt(std::max(0.0, x - b));
    const double h  = b - a;
    if (h <= 0.0) { Wa = Wb = 0.0; return; }
    if (rb <= FPM_EPSILON) { Wa = 2.0 * ra / h; Wb = 0.0; return; } // last panel
    Wb = (2.0 / h) * ((ra - rb) * (ra - rb)) / rb;
    Wa = 2.0 * (1.0 / rb - 1.0 / ra) - Wb;
}

void weights_panel_12(double x, double a, double b, double& Wa, double& Wb) {
    const double ra = std::sqrt(std::max(0.0, x - a));
    const double rb = std::sqrt(std::max(0.0, x - b));
    const double h  = b - a;

    if (h <= 0.0) { Wa = Wb = 0.0; return; }

    // Last panel (or numerically rb ~ 0): exact limit
    if (rb <= FPM_EPSILON) {Wa = (ra > 0.0) ? (2.0 / ra) : 0.0; Wb = 0.0; return; }
    const double Wb_term = (2.0 * (x - a) / rb) + (2.0 * rb) - (4.0 * ra);
    Wb = Wb_term / h;
    Wa = 4.0 / (ra + rb) - 2.0 / ra;
}

double integrate_density(const std::vector<double>& grid,
                                const std::vector<double>& density) {
  const std::size_t n = density.size();
  if (n < 2) {
    return 0.0;
  }
  double total = 0.0;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    const double h = grid[i + 1] - grid[i];
    total += 0.5 * h * (density[i] + density[i + 1]);
  }
  return total;
}

// Heat kernel G*(t, dbeta) = exp(-dbeta^2/(2t)) / sqrt(2 pi t)
inline double Gstar(double t, double delta) {
  if (t <= FPM_EPSILON) return 0.0;
  const double denom = std::sqrt(2.0 * M_PI * t);
  return safe_exp(-(delta * delta) / (2.0 * t)) / denom;
}

// Computes the uniform-avergage image weight for the BM case
inline double averaged_image_term(double t, double beta_t, const RD_Params& pars) {
    if (t <= FPM_EPSILON) return 0.0;

    // Handle the point-start case (no averaging needed)
    if (!pars.sp_var) {
        const double num = pars.zU_scaled - beta_t; // (z - b) in scaled units
        return (num / t) * Gstar(t, num);
    }

    double z_lo = pars.zL_scaled;
    double z_hi = pars.zU_scaled;
    if (z_hi < z_lo) std::swap(z_lo, z_hi);
    const double span = z_hi - z_lo;
    if (span <= FPM_EPSILON) {
        return 0.0;
    }

    // Elegant formulation using the heat kernel (Gstar)
    const double G_at_lower_bound = Gstar(t, z_lo - beta_t);
    const double G_at_upper_bound = Gstar(t, z_hi - beta_t);

    return (1.0 / span) * (G_at_lower_bound - G_at_upper_bound);
}

/**
 * @brief Determines the number of integration steps for a given t_max.
 *
 * This function calculates the number of steps based on a desired
 * "fineness" (step size), enforcing a minimum number of steps,
 * and ensuring the final count is even (which is often required
 * for algorithms like Simpson's rule).
 *
 * @param t_max The maximum time (or theta) for the integration.
 * @param fineness The desired approximate step size (e.g., 0.01).
 * @param min_steps The minimum number of steps to use (e.g., 100).
 * @return The calculated (even) number of steps.
 */
int calculate_num_steps(double t_max, double fineness = 0.01, int min_steps = 100) {
    // 1. Calculate ideal steps based on fineness
    // We use std::ceil to ensure we don't under-shoot.
    // e.g., t_max=1.005 / 0.01 = 100.5, which needs 101 steps.
    double ideal_steps_double = std::ceil(t_max / fineness);

    // 2. Enforce the minimum number of steps
    int num_steps = static_cast<int>(
        std::max(static_cast<double>(min_steps), ideal_steps_double)
    );

    // 3. Ensure the number of steps is even
    if (num_steps % 2 != 0) {
        num_steps++;
    }

    num_steps = (num_steps < 2) ? 2 : num_steps;

    return num_steps;
}

namespace util {
/**
 * @class AkimaSpline
 * @brief Implements the Akima (or "Akima C1") piecewise cubic Hermite interpolator.
 *
 * This interpolator is C1 continuous (it has a continuous first derivative).
 * Its main advantage over a standard cubic spline is that it's "local,"
 * meaning the polynomial for a given interval is determined only by nearby points,
 * which prevents the wild oscillations (Runge's phenomenon) that can
 * plague global interpolators.
 *
 * It's particularly good for non-monotonic data or data with
 * abrupt changes in slope, as it tends to produce a more "natural"
 * and less "wiggly" fit.
 */
class AkimaSpline {
public:
    /**
     * @brief Constructs the Akima spline interpolator.
     *
     * The interpolator is built and coefficients are pre-calculated
     * immediately upon construction.
     *
     * @param t The vector of independent variable values (e.g., time).
     * Must be sorted in ascending order.
     * @param nu The vector of dependent variable values.
     * @throws std::runtime_error if t and nu sizes don't match,
     * if t is not sorted, or if there are fewer than 5 points
     * (Akima needs at least 5 points to compute its tangents correctly).
     */
    AkimaSpline(const std::vector<double>& t, const std::vector<double>& nu) {
        if (t.size() != nu.size()) {
            throw std::runtime_error("AkimaSpline: t and nu vectors must have the same size.");
        }
        if (t.size() < 5) {
            throw std::runtime_error("AkimaSpline: requires at least 5 data points for full algorithm.");
            // Note: Could implement fallbacks (linear, quadratic) for < 5 points,
            // but for this solver, it's better to enforce a minimum grid.
        }

        n_ = t.size();
        t_ = t;
        nu_ = nu;

        // Verify that t is sorted
        for (size_t i = 0; i < n_ - 1; ++i) {
            if (t_[i] >= t_[i + 1]) {
                throw std::runtime_error("AkimaSpline: t vector must be strictly increasing.");
            }
        }

        // Pre-calculate the derivatives (tangents) at each point `t_i`
        // This is the core of the Akima algorithm.
        calculate_tangents();
    }

    /**
     * @brief Interpolates the value at a new point t_val.
     * @param t_val The point at which to interpolate.
     * @return The interpolated value nu(t_val).
     */
    double interpolate(double t_val) const {
        // Find the correct interval [t_i, t_{i+1}] for t_val
        // std::upper_bound finds the first element > t_val.
        // So, `it` points to t_{i+1}, and `i` will be its index.
        auto it = std::upper_bound(t_.begin(), t_.end(), t_val);

        // Handle edge cases (extrapolation)
        if (it == t_.begin()) {
            // t_val is before the first point
            return nu_.front(); // Clamping
            // Or could do linear extrapolation:
            // return nu_[0] + (t_val - t_[0]) * d_[0];
        }
        if (it == t_.end()) {
            // t_val is after the last point
            return nu_.back(); // Clamping
            // Or could do linear extrapolation:
            // return nu_[n_-1] + (t_val - t_[n_-1]) * d_[n_-1];
        }

        // `it` points to t_[i], so we want the interval [i-1, i]
        size_t i = std::distance(t_.begin(), it) - 1;

        double h = t_[i + 1] - t_[i];
        if (h == 0.0) {
            // Should be caught by the sorted check, but good to have
            return nu_[i];
        }

        // Normalize t_val to s in [0, 1]
        double s = (t_val - t_[i]) / h;

        // Apply the standard C1 Hermite cubic polynomial
        double s2 = s * s;
        double s3 = s * s2;

        double h00 = 2 * s3 - 3 * s2 + 1;
        double h10 = s3 - 2 * s2 + s;
        double h01 = -2 * s3 + 3 * s2;
        double h11 = s3 - s2;

        return h00 * nu_[i] + h10 * h * d_[i] + h01 * nu_[i + 1] + h11 * h * d_[i + 1];
    }

    /**
     * @brief Computes the first derivative at a new point t_val.
     * @param t_val The point at which to compute the derivative.
     * @return The derivative d(nu)/d(t) at t_val.
     */
    double derivative(double t_val) const {
        // Find the correct interval [t_i, t_{i+1}]
        auto it = std::upper_bound(t_.begin(), t_.end(), t_val);

        // Handle edge cases
        if (it == t_.begin()) {
            return d_.front(); // Constant derivative extrapolation
        }
        if (it == t_.end()) {
            return d_.back(); // Constant derivative extrapolation
        }

        size_t i = std::distance(t_.begin(), it) - 1;

        double h = t_[i + 1] - t_[i];
        if (h == 0.0) {
            // This case is tricky. The derivative is technically infinite
            // or undefined. We can return the average of the tangents
            // or just the left tangent.
            return d_[i];
        }

        // Normalize t_val to s in [0, 1]
        double s = (t_val - t_[i]) / h;

        // We need the derivative of the Hermite polynomial *with respect to t_val*.
        double s2 = s * s;

        double dh00_ds = 6 * s2 - 6 * s;
        double dh10_ds = 3 * s2 - 4 * s + 1;
        double dh01_ds = -6 * s2 + 6 * s;
        double dh11_ds = 3 * s2 - 2 * s;

        double dp_ds = dh00_ds * nu_[i] + dh10_ds * h * d_[i] +
                       dh01_ds * nu_[i + 1] + dh11_ds * h * d_[i + 1];

        return dp_ds / h;
    }


private:
    size_t n_;
    std::vector<double> t_;   // x-values
    std::vector<double> nu_;  // y-values
    std::vector<double> d_;   // derivatives (tangents) at each point

    /**
     * @brief Pre-calculates the tangents (derivatives) at each data point.
     *
     * This is the core logic of the Akima spline. It uses a weighted
     * average of slopes from adjacent intervals to find a "natural"
     * looking tangent, avoiding the oscillations of standard splines.
     */
    void calculate_tangents() {
        d_.resize(n_);
        std::vector<double> m(n_ - 1); // Slopes of intervals [t_i, t_{i+1}]

        // 1. Calculate the slopes of the n-1 intervals
        for (size_t i = 0; i < n_ - 1; ++i) {
            m[i] = (nu_[i + 1] - nu_[i]) / (t_[i + 1] - t_[i]);
        }

        // 2. We need "ghost" slopes at the ends to handle boundaries.
        // We add two on each side: m[-2], m[-1], m[0], ..., m[n-2], m[n-1], m[n]
        // (total n+3 slopes for n points)
        // m[-1] = 2*m[0] - m[1]
        // m[-2] = 2*m[-1] - m[0] = 3*m[0] - 2*m[1]
        // m[n-1] = 2*m[n-2] - m[n-3]
        // m[n] = 2*m[n-1] - m[n-2] = 3*m[n-2] - 2*m[n-3]
        
        std::vector<double> m_ext(n_ + 3);
        for (size_t i = 0; i < n_ - 1; ++i) {
            m_ext[i + 2] = m[i];
        }

        m_ext[1] = 2.0 * m_ext[2] - m_ext[3];
        m_ext[0] = 2.0 * m_ext[1] - m_ext[2];
        m_ext[n_ + 1] = 2.0 * m_ext[n_] - m_ext[n_ - 1];
        m_ext[n_ + 2] = 2.0 * m_ext[n_ + 1] - m_ext[n_];


        // 3. Calculate the weighted average for the tangents
        for (size_t i = 0; i < n_; ++i) {
            // We need slopes from m_ext[i], m_ext[i+1], m_ext[i+2], m_ext[i+3]
            // which correspond to original slopes m[i-2], m[i-1], m[i], m[i+1]
            
            // Get weights w1 = |m_{i+1} - m_i| and w2 = |m_{i-1} - m_{i-2}|
            // (using the extended m_ext indices)
            double w1 = std::abs(m_ext[i + 3] - m_ext[i + 2]);
            double w2 = std::abs(m_ext[i + 1] - m_ext[i]);

            if (w1 + w2 == 0.0) {
                // Special case: all four slopes are equal (linear segment)
                // or w1=0 and w2=0 (e.g. at a local extremum)
                // Use the average of the two middle slopes.
                d_[i] = (m_ext[i + 1] + m_ext[i + 2]) / 2.0;
            } else {
                // Standard weighted average
                // d_i = (w1 * m_{i-1} + w2 * m_i) / (w1 + w2)
                d_[i] = (w1 * m_ext[i + 1] + w2 * m_ext[i + 2]) / (w1 + w2);
            }
        }
    }
};

} // namespace util


#endif // utils_reducible_diffusion_h


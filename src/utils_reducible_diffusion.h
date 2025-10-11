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
#include <Rcpp.h>
#include <quadmath.h>
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

// Right-end quadratic fit (last 3 points) => ν′(τ), ν″(τ), and ν(mid) via Newton form
struct RightQuad {
    double nu_prime;   // ν′(τ)
    double nu_second;  // ν″(τ)
    double nu_at_mid;  // quadratic ν at midpoint of last panel
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

    // ν′(x2) and ν″(x2) for the quadratic that interpolates the last 3 points
    out.nu_second = 2.0 * f012;
    out.nu_prime  = f01 + f012 * ((x2 - x0) + (x2 - x1));  // p'(x2) = f[x0,x1] + f[x0,x1,x2]((x2-x0)+(x2-x1))

    // ν(mid) from the same quadratic
    const double mid_m_x0 = mid - x0;
    const double mid_m_x1 = mid - x1;
    out.nu_at_mid = y0 + f01 * (mid_m_x0) + f012 * (mid_m_x0) * (mid_m_x1);

    return out;
}

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

    // Lagrange basis at θq
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

// Seed first nodes of ν^f on a θ-grid using the Abel approximation.
// Returns the number of nodes filled (≥1).
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

    // Last panel (or numerically rb≈0): exact limit
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

#endif // utils_reducible_diffusion_h


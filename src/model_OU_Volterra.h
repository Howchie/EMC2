#ifndef OU_HITTING_TIME_H
#define OU_HITTING_TIME_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <Rcpp.h>
#include <quadmath.h>
using namespace Rcpp;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Section 1: Core Constants and Helper Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Threshold for switching to Abel approximation, based on the scaled time.
// The paper notes solutions are "visually indistinguishable" up to t=0.02.
const double SMALL_T_SCALED_THRESHOLD = 0.02;
const double FPM_EPSILON = 1e-12;

// Standard normal CDF
double normal_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

inline double sqrt_pos(double x) {
  return std::sqrt(x > 0.0 ? x : 0.0);
}

inline double safe_exp(double x) {
  // avoid denormals/overflow
  if (x < -745.0) return 0.0;
  if (x >  709.0) return std::exp(709.0);
  return std::exp(x);
}
inline double safe_div(double num, double den) {
  if (!std::isfinite(num) || !std::isfinite(den) || den <= 0.0) return 0.0;
  return num / den;
}

// Prepare and scale parameters for the OU hitting time problem.
struct ScaledParams {
    bool is_trivial;
    double pdf, cdf;
    double t_scaled, z_scaled, b_scaled;
};

// Centralized function for parameter validation, scaling, and handling edge cases.
inline ScaledParams prepare_ou_params(double t, double lambda, double theta, double sigma, double z0, double b_level) {
    ScaledParams p;
    p.is_trivial = true;

    if (t < 0.0 || lambda <= 0.0 || sigma <= 0.0) {
        p.pdf = 0.0; p.cdf = 0.0; // Or throw error
        return p;
    }
    if (t < FPM_EPSILON) {
        p.pdf = 0.0; p.cdf = 0.0;
        return p;
    }

    const double c = std::sqrt(lambda) / sigma;
    p.t_scaled = lambda * t;
    p.z_scaled = c * (z0 - theta);
    p.b_scaled = c * (b_level - theta);

    if (p.z_scaled < p.b_scaled) { // Hitting from below -> flip problem
        p.z_scaled = -p.z_scaled;
        p.b_scaled = -p.b_scaled;
    }

    if (p.z_scaled <= p.b_scaled + FPM_EPSILON) { // Already hit
        p.pdf = 0.0; p.cdf = 1.0;
        return p;
    }

    if (std::abs(p.b_scaled) < FPM_EPSILON) { // Analytical case for b_scaled = 0
        double sinh_t = std::sinh(p.t_scaled);
        if (sinh_t < FPM_EPSILON) {
            p.pdf = 0.0; p.cdf = 0.0;
        } else {
            double common_arg = -p.z_scaled * std::exp(-0.5 * p.t_scaled) / std::sqrt(sinh_t);
            p.cdf = 2.0 * normal_cdf(common_arg);
            const double exp_term = std::exp(-0.5 * common_arg * common_arg + 0.5 * p.t_scaled);
            const double g_scaled = (p.z_scaled * exp_term) / (std::sqrt(2.0 * M_PI) * std::pow(sinh_t, 1.5));
            p.pdf = lambda * g_scaled;
        }
        return p;
    }

    p.is_trivial = false;
    return p;
}

inline void make_theta_and_tau_grids(double t_scaled, int N, double tau_max,
                                     std::vector<double>& theta_grid,
                                     std::vector<double>& tau_grid_theta_mapped,
                                     std::vector<double>& tau_grid_uniform,
                                     std::vector<double>& tau_grid_right_clustered,
                                     double& theta_max,
                                     double alpha = 1.5)
{
    theta_max = std::exp(t_scaled) - 1.0;
    const double h_theta = theta_max / N;

    theta_grid.resize(N + 1);
    tau_grid_theta_mapped.resize(N + 1);
    tau_grid_uniform.resize(N + 1);
    tau_grid_right_clustered.resize(N + 1);

    for (int j = 0; j <= N; ++j) {
        const double theta_j = j * h_theta;
        theta_grid[j] = theta_j;

        const double s = 1.0 + theta_j;
        tau_grid_theta_mapped[j] = 0.5 * (s * s - 1.0);

        tau_grid_uniform[j] = (tau_max * j) / N;

        const double u = static_cast<double>(j) / N;
        const double clustered = tau_max * (1.0 - std::pow(1.0 - u, alpha));
        tau_grid_right_clustered[j] = (j == N) ? tau_max : clustered;
    }

    tau_grid_theta_mapped.back() = tau_max;
    tau_grid_uniform.back() = tau_max;
    tau_grid_right_clustered.back() = tau_max;
}

// Convert theta values into tau or vice-versa
inline double tau_of_theta(double th) {
        // tau(theta) = 0.5 * ((1+theta)^2 - 1) = theta + 0.5 theta^2
        return th + 0.5 * th * th;
}
inline double theta_of_tau(double tau){ 
	return std::sqrt(2.0*tau + 1.0) - 1.0;
}

inline double t_from_tau_theta(double x, bool use_theta) {
    // Primary: tau = 0.5 * (exp(2 t) - 1)  ⇒  t = 0.5 * log(1 + 2 tau)
    if (use_theta) {
        // Fallback: theta = exp(t) - 1 ⇒ t = log(1 + theta)
        if (x > -1.0) {
            const double arg = std::max(0.0, 1.0 + x);
            if (arg > 0.0) {
                return std::log1p(x);
            }
        }
    } else {
        if (x > -0.5) {
            const double arg = std::max(0.0, 1.0 + 2.0 * x); // guard tiny negatives
            if (arg > 0.0) {
                return 0.5 * std::log1p(2.0 * x);
            }
        }
    }

    return 0.0;
}

enum class LastPanelStrategy { None, USplit, QuadraticClosed, Richardson };

inline void effective_positions(double v, double vp, double z, double b,
                                double& Delta, double& one_minus_vp) {
  const double z_eff = z * (1.0 - v);
  const double b_eff = b * (1.0 - vp);
  Delta = z_eff-b_eff;           
  one_minus_vp = 1.0 - vp;
}

// Linear-product weights for the integral of (Linear(H) / (x-s)^{3/2}) ds on [a,b]
inline void weights_panel_32(double x, double a, double b, double& Wa, double& Wb) {
    const double ra = std::sqrt(std::max(0.0, x - a));
    const double rb = std::sqrt(std::max(0.0, x - b));
    const double h  = b - a;
    if (h <= 0.0) { Wa = Wb = 0.0; return; }
    if (rb <= FPM_EPSILON) { Wa = 2.0 * ra / h; Wb = 0.0; return; } // last panel
    Wb = (2.0 / h) * ((ra - rb) * (ra - rb)) / rb;
    Wa = 2.0 * (1.0 / rb - 1.0 / ra) - Wb;
}

inline void weights_panel_12(double x, double a, double b, double& Wa, double& Wb) {
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

// Right-end quadratic fit (last 3 points) => ν′(τ), ν″(τ), and ν(mid) via Newton form
struct RightQuad {
    double nu_prime;   // ν′(τ)
    double nu_second;  // ν″(τ)
    double nu_at_mid;  // quadratic ν at midpoint of last panel
};

inline RightQuad right_end_quadratic(const std::vector<double>& x,
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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Section 2: Backward Method Components (for CDF)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// The non-singular part of the Volterra kernel from page 11 (backward variables).
// This is phi^b_B(v, v') in the paper's notation.
// [[Rcpp::export]]
double kernel_backward(double v, double vp, double b) {
    if (std::abs(b) < FPM_EPSILON) {
        return 0.0;
    }
	double dv = v - vp;
    double s = 2.0 - v - vp;
	if (dv < 0.0 || s <= FPM_EPSILON) return 0.0;
    double expo = safe_exp(-b*b*dv/s); // clamp extreme values
    double term1 = (2.0 * b / std::sqrt(M_PI));
    return term1 * (((1.0 - vp) * expo) / (s * std::sqrt(s)));
           
}

// Solve nu_b(v) on an arbitrary increasing grid v_grid[0..N], v_grid[0]=0, v_grid[N]=v_max, using the Stieltjes-trapezoid rule.
// This is essential for the change-of-variable strategy, avoiding interpolation.
// [[Rcpp::export]]
std::vector<double> solve_nu_b_on_grid(const std::vector<double>& v_grid, double b) {
	const int N = static_cast<int>(v_grid.size()) - 1;
    if (N < 1) return {1.0};
    
    std::vector<double> F(N + 1, 0.0);
    F[0] = 1.0; // Boundary condition F(0) = 1

    for (int k = 1; k <= N; ++k) {
        const double v_k = v_grid[k];
        double sum_term = 0.0;

        // Sum contributions from previous intervals [v_{i-1}, v_i] for i=1..k-1
        for (int i = 1; i < k; ++i) {
            const double v_i   = v_grid[i];
            const double v_im1 = v_grid[i-1];
            
            const double K_k_i   = kernel_backward(v_k, v_i,   b);
            const double K_k_im1 = kernel_backward(v_k, v_im1, b);
            
            // Stieltjes weight for this interval
            const double w = 2.0 * (std::sqrt(v_k - v_im1) - std::sqrt(v_k - v_i));
            sum_term += 0.5 * (K_k_i * F[i] + K_k_im1 * F[i-1]) * w;
        }

        // Final interval [v_{k-1}, v_k] contains the unknown F[k]
        const double v_km1 = v_grid[k-1];
        const double w_last = 2.0 * std::sqrt(v_k - v_km1);
        
        const double K_k_k   = kernel_backward(v_k, v_k,   b);
        const double K_k_km1 = kernel_backward(v_k, v_km1, b);
        
        // Equation to solve for F[k]:
        // F[k] = 1.0 + sum_term + 0.5 * (K_k_k * F[k] + K_k_km1 * F[k-1]) * w_last
        double numerator = 1.0 + sum_term + 0.5 * K_k_km1 * F[k-1] * w_last;
        double denominator = 1.0 - 0.5 * K_k_k * w_last;
        
        F[k] = numerator / ((std::abs(denominator) < FPM_EPSILON) ? FPM_EPSILON : denominator);
    }
    return F;
}

// Indefinite integral of P(s)/sqrt(A-s) ds, where P(s) is a quadratic.
// P(s) = c2*s^2 + c1*s + c0
// Returns the evaluated definite integral from s=0 to s=2.
inline double eval_quad_sqrt_integral(double A, double c2, double c1, double c0) {
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
inline double alpha_coeff(double dt, double z_eff) {
    if (dt <= 0) return 0.0;
    const double A = dt / z_eff;                // expect A = 2 in both rows
    const double I = eval_quad_sqrt_integral(A, 1.0, -3.0,  2.0); // (1-s)(2-s)
    return 0.5 * std::sqrt(z_eff) * I;
}
inline double beta_coeff(double dt, double z_eff) {
    if (dt <= 0) return 0.0;
    const double A = dt / z_eff;
    const double I = eval_quad_sqrt_integral(A, -1.0, 2.0,  0.0); // s(2-s)
    return       std::sqrt(z_eff) * I;
}
inline double gamma_coeff(double dt, double z_eff) {
    if (dt <= 0) return 0.0;
    const double A = dt / z_eff;
    const double I = eval_quad_sqrt_integral(A,  1.0,-1.0,  0.0); // s(s-1)
    return 0.5 * std::sqrt(z_eff) * I;
}

// Analytical solution for v^b(v) for small v (Abel approximation from Eq. 26)
// Here, v corresponds to the paper's theta.
// [[Rcpp::export]]
double abel_approx_nu_b(double v, double b) {
    if (std::abs(b) < FPM_EPSILON) {
        return 1.0;
    }
    return 2.0 * std::exp( (b * b * v) / 2.0) * normal_cdf(b * std::sqrt(v));
}

/**
 * @brief Solves the Volterra equation for v^b(v) using the block-by-block method.
 *
 * Implements the quadratic interpolation method from Section 5.1.2 of Lipton & Kaushansky.
 * This provides higher accuracy than the trapezoidal rule, especially for singular kernels.
 *
 * @param v_max The final value of v (time) to solve up to.
 * @param b The scaled boundary level.
 * @param num_steps The number of steps for the grid (MUST BE AN EVEN NUMBER).
 * @return A vector containing the solution v^b at each grid point.
 */

// Stieltjes-trapezoid history sum for node n (exclude the last interval)
inline double stieltjes_trap_history_backward_old(int n,
                                     const std::vector<double>& F,
                                     double h,
                                     double b)
{
    if (n <= 1) return 0.0;
    const double tn = n * h;
    double S = 0.0;
    for (int i = 0; i <= n - 2; ++i) {
        const double ti   = i     * h;
        const double tip1 = (i+1) * h;
        // weight: integral_{ti}^{ti+1} 1/sqrt (tn - s) ds = 2(sqrt (tn - ti) - sqrt (tn - tip1))
        const double w = 2.0 * (std::sqrt(tn - ti) - std::sqrt(tn - tip1));
        const double Ki   = kernel_backward(tn, ti,   b);
        const double Kip1 = kernel_backward(tn, tip1, b);
        S += 0.5 * (Ki * F[i] + Kip1 * F[i+1]) * w;
    }
    return S;
}

// Stieltjes-trapezoid history sum for node n (exclude the last interval)
inline double stieltjes_trap_history_backward(
    int n, int last_panel_max,const std::vector<double>& F, double h, double b, bool force_simple=false)
{
    if (n <= 1) return 0.0;
    if (F.empty()) return 0.0;
    const int size_as_int = static_cast<int>(F.size());
    if (size_as_int <= 1) return 0.0;
    const auto K  = kernel_backward;
    const double tn = n * h;
    const int max_known = std::min(size_as_int - 1, n - 1);
    int last_panel = std::min(last_panel_max, max_known - 1);
    if (last_panel < 0) return 0.0;

    if (max_known < 1) return 0.0;
    double S = 0.0;
    if (last_panel < 0) return 0.0;
    for (int i = 0; i <= last_panel; ++i) {
        const double a   = i * h;
        const double bnd = (i + 1) * h;
        const double mid = 0.5 * (a + bnd);
        const double d   = bnd - a;
        if (!(d > 0.0)) continue;
        const double Fa = F[i];
        const double Fb = F[i + 1];
        bool have_quad = false;
        double Fm = 0.5 * (Fa + Fb);
        if (i > 0 && i + 1 <= max_known && !force_simple) {
            const double F_prev = F[i - 1];
            const double F_next = F[i + 1];
            Fm = (-0.125) * F_prev + 0.75 * Fa + 0.375 * F_next;
            have_quad = true;
        } else if (i == 0 && i + 2 <= max_known && !force_simple) {
            const double F2 = F[i + 2];
            Fm = 0.375 * Fa + 0.75 * Fb - 0.125 * F2;
            have_quad = true;
        }
        const double ra = sqrt_pos(tn - a);
        const double rb = sqrt_pos(tn - bnd);
        if (!have_quad) {
            double Wb;
            if (rb <= 0.0) {
                Wb = (2.0 / d) * ((tn - a) * ra - (ra * ra * ra) / 3.0);
            } else {
                Wb = (2.0 / d) * (((tn - a) * (ra - rb)) - (ra * ra * ra - rb * rb * rb) / 3.0);
            }
            const double Wa = 2.0 * (ra - rb) - Wb;
            const double Ka = K(tn, a, b);
            const double Kb = K(tn, bnd, b);
            S += Wa * (Ka * Fa) + Wb * (Kb * Fb);
            continue;
        }
        const double Qa = K(tn, a, b)   * Fa;
        const double Qm = K(tn, mid, b) * Fm;
        const double Qb = K(tn, bnd, b) * Fb;
        const double A  = tn - a;
        const double A2 = A * A;
        const double ra2 = ra * ra;
        const double rb2 = rb * rb;
        const double ra3 = ra2 * ra;
        const double rb3 = rb2 * rb;
        const double ra5 = ra3 * ra2;
        const double rb5 = rb3 * rb2;
        const double J0 = 2.0 * (ra - rb);
        const double term1_a = A * ra - ra3 / 3.0;
        const double term1_b = A * rb - rb3 / 3.0;
        const double I1 = 2.0 * (term1_a - term1_b);
        const double term2_a = A2 * ra - (2.0 / 3.0) * A * ra3 + (1.0 / 5.0) * ra5;
        const double term2_b = A2 * rb - (2.0 / 3.0) * A * rb3 + (1.0 / 5.0) * rb5;
        const double I2 = 2.0 * (term2_a - term2_b);
        const double J1 = a * J0 + I1;
        const double J2 = a * a * J0 + 2.0 * a * I1 + I2;
        const double d2 = d * d;
        if (!(d2 > 0.0)) continue;
        const double s0 = a;
        const double s1 = mid;
        const double s2 = bnd;
        const double w0 = (2.0 / d2) * (J2 - (s1 + s2) * J1 + s1 * s2 * J0);
        const double w1 = (-4.0 / d2) * (J2 - (s0 + s2) * J1 + s0 * s2 * J0);
        const double w2 = (2.0 / d2) * (J2 - (s0 + s1) * J1 + s0 * s1 * J0);
        if (!std::isfinite(w0) || !std::isfinite(w1) || !std::isfinite(w2)) {
            double Wb;
            if (rb <= 0.0) {
                Wb = (2.0 / d) * ((tn - a) * ra - (ra * ra * ra) / 3.0);
            } else {
                Wb = (2.0 / d) * (((tn - a) * (ra - rb)) - (ra * ra * ra - rb * rb * rb) / 3.0);
            }
            const double Wa = 2.0 * (ra - rb) - Wb;
            S += Wa * Qa + Wb * Qb;
            continue;
        }
        S += w0 * Qa + w1 * Qm + w2 * Qb;
    }
    return S;
}

// Seed first nodes of ν^f on a θ-grid using the Abel approximation.
// Returns the number of nodes filled (≥1).
inline int seed_nu_b_on_grid(const std::vector<double>& grid,
                                        double b,
                                        std::vector<double>& F,
                                        double t_cut = 0.02) 
{
    const int N = grid.size();
    int J = 0;
    for (;J < N; ++J) {
        const double v = grid[J];
        const double t = -std::log(1-v);
        if (t > t_cut) break;
        F[J] = abel_approx_nu_b(v, b);
    }
    // Ensure we don't leave the first step under-determined.
    if (J < 2 && N >= 2) {
        const double v = grid[J];
        F[1] = abel_approx_nu_b(v, b);
        J = std::max(J, 2);
    }
    if (J < 1) J = 1;
    if (J & 1) J -= 1;  // largest even ≤ J
    return J;           // all indices 0..J are "seeded/locked"
}

// Quadratic block step on a uniform grid t_j = j*h, j = 0..N, N even
// [[Rcpp::export]]
std::vector<double> solve_nu_b_block_by_block(double v_max, double b,
                                              int num_steps, 
                                              std::vector<double> &v_grid,
                                              std::vector<double> F,
                                              int k0_seeded)
{
    // enforce even N >= 2
    if (num_steps < 2) num_steps = 2;
    if (num_steps % 2) ++num_steps;
    const auto K = kernel_backward;
    const int N = num_steps;
    const double h = v_max / N;
	v_grid.resize(N+1);
    for (int i = 0; i <= N; ++i) {
        v_grid[i] = i * h;
    }

    if (std::abs(F[0]) < 1e-300) F[0] = 1.0;

    const int M = N / 2;
    const int first_unseeded_block = std::max(0, k0_seeded / 2); // integer div
    const int m0 = std::min(M, first_unseeded_block);  // block start index in blocks
    for (int m = m0; m < M; ++m) {
        const int j0 = 2*m;       // 2m
        const int j1 = 2*m + 1;   // 2m+1
        const int j2 = 2*m + 2;   // 2m+2
		

        const double t0 = j0 * h;
        const double t1 = j1 * h;
        const double t2 = j2 * h;
		const double tmid = t0 + 0.5*h;
        // history up to v0 (i.e., up to j0) ----
        const double S1 = stieltjes_trap_history_backward(j1, j0-1, F, h, b);
        const double S2 = stieltjes_trap_history_backward(j2, j0-1, F, h, b);

        // coefficients for product integration with quadratic F(s)
        const double a1 = alpha_coeff(h,  0.5*h);
        const double b1 = beta_coeff (h,  0.5*h);
        const double g1 = gamma_coeff(h,  0.5*h);

        const double a2 = alpha_coeff(2*h, h);
        const double b2 = beta_coeff (2*h, h);
        const double g2 = gamma_coeff(2*h, h);

        // kernel evaluations
        const double K1_0 = K(t1, t0, b);
        const double K1_mid = K(t1, tmid, b);
        const double K1_1 = K(t1, t1, b);
        // K(t1, t2) is zero due to causality, so we don't need it.

        const double K2_0 = K(t2, t0, b);
        const double K2_1 = K(t2, t1, b);
        const double K2_2 = K(t2, t2, b);
		
        const double G1 = 1;
        const double G2 = 1;
        // ---- 4) Assemble and solve ----
        double F1_pred, F2_pred;

        if (j0 >= 1) {
            // Interior: solve (1 - g1*K11 - (3/8)*b1*K1mid) * F1 = RHS
            const double den1_pred = 1.0 - g1*K1_1 - (3.0/8.0)*b1*K1_mid;
            const double rhs1_pred = G1 + S1
                + a1*K1_0*F[j0]
                + b1*K1_mid*( (-0.125)*F[j0-1] + 0.75*F[j0] );
            F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
        } else {
            // First block: linear midpoint Fmid ≈ 0.5(F0 + F1)
            const double den1_pred = 1.0 - g1*K1_1 - 0.5*b1*K1_mid;
            const double rhs1_pred = G1 + S1 + a1*K1_0*F[j0] + 0.5*b1*K1_mid*F[j0];
            F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
        }

        // Even row predictor is the usual 1×1 (uses predicted F1)
        {
            const double den2_pred = 1.0 - g2*K2_2;
            const double rhs2_pred = G2 + S2 + a2*K2_0*F[j0] + b2*K2_1*F1_pred;
            F2_pred = rhs2_pred / ((std::abs(den2_pred) < FPM_EPSILON) ? FPM_EPSILON : den2_pred);
        }

        // Odd equation:  A*F1 + B*F2 = C
        const double A = 1.0 - (3.0/4.0) * b1 * K1_mid - g1 * K1_1;
        const double B = (1.0/8.0) * b1 * K1_mid;
        const double C = G1 + S1 + a1 * K1_0 * F[j0]
                            + (3.0/8.0) * b1 * K1_mid * F[j0];

        // Even equation: D*F1 + E*F2 = R
        const double D = - b2 * K2_1;
        const double E = 1.0 - g2 * K2_2;
        const double R = G2 + S2 + a2 * K2_0 * F[j0];

        // Solve the 2x2 once. 
        const double det = A*E - B*D;
        const double inv_det = (std::abs(det) < FPM_EPSILON) ? (1.0/FPM_EPSILON) : (1.0/det);

        double F1_corr = ( E*C - B*R) * inv_det;
        double F2_corr = (-D*C + A*R) * inv_det;

        const double theta = 1; // set to e.g., 0.7 if you want under-relaxation
        F[j1] = theta*F1_corr + (1.0 - theta)*F1_pred;
        F[j2] = theta*F2_corr + (1.0 - theta)*F2_pred;

    }

    return F;
}

std::vector<double> solve_nu_b_block_with_abel(double v_max, double z, double b,
                                               int num_steps,
                                               std::vector<double>& v_grid,
                                               double t_cut = SMALL_T_SCALED_THRESHOLD)
{
    std::vector<double> F(num_steps+1, 0.0);
    // seeding uses the *same* grid we’re about to build
    v_grid.resize(num_steps+1);
    for (int i=0;i<=num_steps;++i) v_grid[i] = (v_max/num_steps) * i;
    const int k0 = seed_nu_b_on_grid(v_grid, b, F, t_cut);
    return solve_nu_b_block_by_block(v_max, b, num_steps, v_grid, F, k0);
}


// Regularized integrands: factor out 1/sqrt(v - vp)^(3/2) so we can
// integrate with Stieltjes weights on the same uniform v'-grid.
inline double G_integrand_smooth(double v, double vp, double z, double b, double nu) {
    const double dv = v - vp;
    const double s = 2.0 - v - vp;
	double Delta, one_minus_vp;
	const double s_eff = (s > 1e-14) ? s : 1e-14;
	const double dv_eff = (dv > 1e-14) ? dv : 1e-14;
	effective_positions(v, vp, z, b, Delta, one_minus_vp);
    const double expo = safe_exp(-Delta * Delta / (dv_eff * s_eff));

    const double denom = std::sqrt(M_PI) * (s_eff * std::sqrt(s_eff)); 
	const double num = Delta * expo * one_minus_vp * nu;
    if (denom < FPM_EPSILON) return 0.0;
    return num / denom;
}

// High-order integrator for the final CDF calculation.
inline double integrate_cdf_backward_trapezoid(
    double v, const std::vector<double>& v_grid, const std::vector<double>& nu_vals,
    double z_scaled, double b_scaled) {
    
    double integral_sum = 0.0;
    const int N = static_cast<int>(v_grid.size()) - 1;

    for (int i = 0; i < N; ++i) { // Loop over N intervals
        const double a = v_grid[i], b = v_grid[i+1];
        const double delta_v = b - a;
        if (delta_v < FPM_EPSILON) continue;

        const double H_a = G_integrand_smooth(v, a, z_scaled, b_scaled, nu_vals[i]);
		double H_b = 0.0; 
        const double sqrt_v_minus_a = std::sqrt(std::max(0.0, v - a));
 
        // Final interval b==v and Hb=0
        double W_a, W_b;
        if (i == N - 1) {  // Check if we are in the final interval where b ~ v
			// The smooth part H_b at v'=v is zero due to the exponential term.
			 weights_panel_12(v, a, b, W_a, W_b);
        } else {
           // Use the standard product trapezoidal formula for all other intervals.
            H_b = G_integrand_smooth(v, b, z_scaled, b_scaled, nu_vals[i+1]);
			weights_panel_12(v, a, b, W_a, W_b);
        }

        integral_sum += W_a * H_a + W_b * H_b;
    }
    
    return integral_sum;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Section 3: Forward Method Components (for PDF)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// [[Rcpp::export]]
double abel_approx_nu_f(double theta, double z, double b) {
    if (theta < FPM_EPSILON) {
        return 0.0; // nu_f(0) = g(0) = 0
    }
    const double num = b * theta + z - b;

    const double exp_arg1 = 0.5 * b * b * theta + b * (z - b);
    const double cdf_arg = -num / std::sqrt(theta);
    const double term1 = b * safe_exp(exp_arg1) * normal_cdf(cdf_arg);

    const double exp_arg2 = -((b-z) * (b-z))/(2*theta);
    const double term2 = -safe_exp(exp_arg2) / std::sqrt(2.0 * M_PI * theta);

    return term1 + term2;
}

/**
 * @brief Non-singular part of the kernel for the forward Volterra equation (Eq. 7).
 * The full kernel is K(tau, tau') / sqrt(tau - tau').
 */
// [[Rcpp::export]] 
double kernel_forward_tau(double tau, double tau_p, double b) {
    if (std::abs(b) < FPM_EPSILON) return 0.0;
	
	// Pre-calculate the sqrt terms for efficiency and readability
	const double s = std::sqrt(2.0 * tau + 1.0);
    const double sp = std::sqrt(2.0 * tau_p + 1.0);
    const double sum = s + sp;
	
    // Handle diagonal explicitly ratio term cancels (difference of 0) 
    if (std::abs(tau - tau_p) < 1e-15) {
        
        if (s < FPM_EPSILON) return 0.0;
        return -(std::sqrt(2.0/M_PI) * b) / (2*s); 
    }
	
	if (tau <= tau_p + FPM_EPSILON || std::abs(b) < FPM_EPSILON) {
        return 0.0;
    }
    
    if (sum < FPM_EPSILON) return 0.0; // Avoid division by zero
    
    const double ratio = (s - sp) / sum;
    const double expo = safe_exp(-b * b * ratio);
    
    // Eq. 7 Lipton et al (2018), -sign when moving to RHS of equation
    return -(std::sqrt(2.0/M_PI) * b) * expo / sum;
}

// [[Rcpp::export]]
double kernel_forward_theta(double theta, double theta_p, double b) {
    if (std::abs(b) < FPM_EPSILON) return 0.0;

    const double s  = 1.0 + theta;
    const double sp = 1.0 + theta_p;
    const double sum = 2.0 + theta + theta_p;
	// Handle diagonal explicitly
    if (std::abs(s - sp) < 1e-15) {
        // K(tau,tau) = - b / (sqrt (2pi) sqrt s)
		if (s < FPM_EPSILON) return 0.0;
        return -(b / std::sqrt(2.0 * M_PI)) / sqrt_pos(s);
    }
    if (theta <= theta_p + FPM_EPSILON) return 0.0;

    const double ratio = (theta - theta_p) / sum;
    const double expo  = safe_exp(-b*b*ratio);

    // Eq. 10 Lipton et al (2018), -sign when moving to RHS of equation
    return -(2.0*b/std::sqrt(M_PI)) * (expo * sp) / (sum * std::sqrt(sum));
}

/**
 * @brief The non-integral forcing term g(tau) in the forward Volterra equation (Eq. 7).
 */
// [[Rcpp::export]]
double g_term_forward_tau(double tau, double z, double b) {
    if (tau < FPM_EPSILON) {
        return 0.0;
    }
    
    const double s = std::sqrt(2.0 * tau + 1.0);
    const double num = s * b - z;
    const double expo = safe_exp(-num * num / (2.0 * tau));
    const double denom = sqrt_pos(2.0 * M_PI * tau);
	if (denom <= FPM_EPSILON) return 0.0;
    return -expo / denom;
}

// theta-form
// [[Rcpp::export]]
double g_term_forward_theta(double theta, double z, double b) {
    const double s = 1.0 + theta;
	const double term1 = s*s-1.0;
    if (theta <= FPM_EPSILON) return 0.0;
    const double num = s*b - z;
	const double expo = safe_exp(-(num*num)/term1);
	const double denom = sqrt_pos(M_PI * term1);
	if (denom <= FPM_EPSILON) return 0.0;
    return -expo / denom;
}

// Stieltjes-trapezoid history sum for node n (exclude the last interval)
inline double stieltjes_trap_history_forward(
    int n, int last_panel_max,const std::vector<double>& F, double h, double b, bool use_theta, bool force_simple=false)
{
    if (n <= 1) return 0.0;
    if (F.empty()) return 0.0;
    const int size_as_int = static_cast<int>(F.size());
    if (size_as_int <= 1) return 0.0;
    const auto K  = use_theta ? kernel_forward_theta : kernel_forward_tau;
    const double tn = n * h;
    const int max_known = std::min(size_as_int - 1, n - 1);
    int last_panel = std::min(last_panel_max, max_known - 1);
    if (last_panel < 0) return 0.0;

    if (max_known < 1) return 0.0;
    double S = 0.0;
    if (last_panel < 0) return 0.0;
    for (int i = 0; i <= last_panel; ++i) {
        const double a   = i * h;
        const double bnd = (i + 1) * h;
        const double mid = 0.5 * (a + bnd);
        const double d   = bnd - a;
        if (!(d > 0.0)) continue;
        const double Fa = F[i];
        const double Fb = F[i + 1];
        bool have_quad = false;
        double Fm = 0.5 * (Fa + Fb);
        if (i > 0 && i + 1 <= max_known && !force_simple) {
            const double F_prev = F[i - 1];
            const double F_next = F[i + 1];
            Fm = (-0.125) * F_prev + 0.75 * Fa + 0.375 * F_next;
            have_quad = true;
        } else if (i == 0 && i + 2 <= max_known && !force_simple) {
            const double F2 = F[i + 2];
            Fm = 0.375 * Fa + 0.75 * Fb - 0.125 * F2;
            have_quad = true;
        }
        const double ra = sqrt_pos(tn - a);
        const double rb = sqrt_pos(tn - bnd);
        if (!have_quad) {
            double Wb;
            if (rb <= 0.0) {
                Wb = (2.0 / d) * ((tn - a) * ra - (ra * ra * ra) / 3.0);
            } else {
                Wb = (2.0 / d) * (((tn - a) * (ra - rb)) - (ra * ra * ra - rb * rb * rb) / 3.0);
            }
            const double Wa = 2.0 * (ra - rb) - Wb;
            const double Ka = K(tn, a, b);
            const double Kb = K(tn, bnd, b);
            S += Wa * (Ka * Fa) + Wb * (Kb * Fb);
            continue;
        }
        const double Qa = K(tn, a, b)   * Fa;
        const double Qm = K(tn, mid, b) * Fm;
        const double Qb = K(tn, bnd, b) * Fb;
        const double A  = tn - a;
        const double A2 = A * A;
        const double ra2 = ra * ra;
        const double rb2 = rb * rb;
        const double ra3 = ra2 * ra;
        const double rb3 = rb2 * rb;
        const double ra5 = ra3 * ra2;
        const double rb5 = rb3 * rb2;
        const double J0 = 2.0 * (ra - rb);
        const double term1_a = A * ra - ra3 / 3.0;
        const double term1_b = A * rb - rb3 / 3.0;
        const double I1 = 2.0 * (term1_a - term1_b);
        const double term2_a = A2 * ra - (2.0 / 3.0) * A * ra3 + (1.0 / 5.0) * ra5;
        const double term2_b = A2 * rb - (2.0 / 3.0) * A * rb3 + (1.0 / 5.0) * rb5;
        const double I2 = 2.0 * (term2_a - term2_b);
        const double J1 = a * J0 + I1;
        const double J2 = a * a * J0 + 2.0 * a * I1 + I2;
        const double d2 = d * d;
        if (!(d2 > 0.0)) continue;
        const double s0 = a;
        const double s1 = mid;
        const double s2 = bnd;
        const double w0 = (2.0 / d2) * (J2 - (s1 + s2) * J1 + s1 * s2 * J0);
        const double w1 = (-4.0 / d2) * (J2 - (s0 + s2) * J1 + s0 * s2 * J0);
        const double w2 = (2.0 / d2) * (J2 - (s0 + s1) * J1 + s0 * s1 * J0);
        if (!std::isfinite(w0) || !std::isfinite(w1) || !std::isfinite(w2)) {
            double Wb;
            if (rb <= 0.0) {
                Wb = (2.0 / d) * ((tn - a) * ra - (ra * ra * ra) / 3.0);
            } else {
                Wb = (2.0 / d) * (((tn - a) * (ra - rb)) - (ra * ra * ra - rb * rb * rb) / 3.0);
            }
            const double Wa = 2.0 * (ra - rb) - Wb;
            S += Wa * Qa + Wb * Qb;
            continue;
        }
        S += w0 * Qa + w1 * Qm + w2 * Qb;
    }
    return S;
}

// Seed first nodes of ν^f on a θ-grid using the Abel approximation.
// Returns the number of nodes filled (≥1).
inline int seed_nu_f_on_grid(const std::vector<double>& grid,bool use_theta,
                                        double z, double b,
                                        std::vector<double>& F,
                                        double t_cut = SMALL_T_SCALED_THRESHOLD) 
{
    const int N = grid.size();
    int J = 0;
    for (;J < N; ++J) {
        const double th = use_theta ? grid[J] : theta_of_tau(grid[J]);
        const double t = t_from_tau_theta(th,use_theta=true);
        if (t > t_cut) break;
        F[J] = abel_approx_nu_f(th, z, b);
    }
    // Ensure we don't leave the first step under-determined.
    if (J < 2 && N >= 2) {
        const double th = use_theta ? grid[J] : theta_of_tau(grid[J]);
        F[1] = abel_approx_nu_f(th, z, b);
        J = std::max(J, 2);
    }
    if (J < 1) J = 1;
    if (J & 1) J -= 1;  // largest even ≤ J
    return J;           // all indices 0..J are "seeded/locked"
}


/**
 * @brief Solves the Volterra equation for the forward function v^f(tau) on a grid.
 * This adapts the stable Stieltjes-trapezoid solver for a non-constant g(tau) term.
 */
// the boolean "theta transform" denotes which solver to use
std::vector<double> solve_nu_f_grid(const std::vector<double> &grid, double z,
                                    double b, bool use_theta,
                                    std::vector<double> F, int k0_seeded) {
    const int N = (int)grid.size()-1;
	    // Choose the right pair (kernel, g) for the coordinate
    const auto K = use_theta ? kernel_forward_theta : kernel_forward_tau;
    const auto G = use_theta ? g_term_forward_theta : g_term_forward_tau;
    if (k0_seeded < 0) k0_seeded = 0;
    if (k0_seeded == 0 && std::abs(F[0]) < 1e-300) F[0] = G(grid[0], z, b);
    for (int k = std::max(1, k0_seeded); k<=N; ++k) {
        const double xk = grid[k];
        double sum = 0.0;
        for (int i=1; i<k; ++i) {
            const double w = 2.0*( sqrt_pos(xk - grid[i-1])
                                  -sqrt_pos(xk - grid[i]) );
			const double Ki   = K(xk, grid[i],   b);
			const double Kim1 = K(xk, grid[i-1], b);


            sum += 0.5*(Ki*F[i] + Kim1*F[i-1]) * w;
			
        }
        const double w_last = 2.0*sqrt_pos(xk - grid[k-1]);
		const double Kkk    = K(xk, xk, b);
		const double Kkkm1  = K(xk, grid[k-1], b);
		const double gk     = G(xk, z, b);
        const double num = gk + sum + 0.5*Kkkm1*F[k-1]*w_last;
        const double den = 1.0 - 0.5*Kkk*w_last;
        F[k] = num / ((std::abs(den) < FPM_EPSILON) ? FPM_EPSILON : den);
    }
    return F;
}

std::vector<double> solve_nu_f_grid_with_abel(const std::vector<double>& grid,
                                              double z, double b, bool use_theta,
                                              double theta_cut = SMALL_T_SCALED_THRESHOLD)
{
    std::vector<double> F(grid.size(), 0.0);
    const int k0 = seed_nu_f_on_grid(grid, use_theta, z, b, F, theta_cut);
    return solve_nu_f_grid(grid, z, b, use_theta, F, k0);
}


// Quadratic block step on a uniform grid t_j = j*h, j = 0..N, N even
// [[Rcpp::export]]
std::vector<double> solve_nu_f_block_by_block(double t_max, double z, double b,
                                              int num_steps, bool use_theta,
                                              std::vector<double> &t_grid,
                                              std::vector<double> F,
                                              int k0_seeded)
{
    // enforce even N >= 2
    if (num_steps < 2) num_steps = 2;
    if (num_steps % 2) ++num_steps;
	const auto K = use_theta ? kernel_forward_theta : kernel_forward_tau;
    const auto G = use_theta ? g_term_forward_theta : g_term_forward_tau;
    const int N = num_steps;
    const double h = t_max / N;
	t_grid.resize(N+1);
    for (int i = 0; i <= N; ++i) {
        t_grid[i] = i * h;
    }

    if (std::abs(F[0]) < 1e-300) F[0] = G(t_grid[0], z, b);

    const int M = N / 2;
    const int first_unseeded_block = std::max(0, k0_seeded / 2); // integer div
    const int m0 = std::min(M, first_unseeded_block);  // block start index in blocks
    for (int m = m0; m < M; ++m) {
            const int j0 = 2*m;       // 2m
    const int j1 = j0 + 1;    // 2m+1
    const int j2 = j0 + 2;    // 2m+2

    const double t0 = j0 * h;
    const double t1 = j1 * h;
    const double t2 = j2 * h;
    const double tmid = t0 + 0.5*h;

    // --- 1) history strictly up to (and including) panel j0-1 ---
    // IMPORTANT: exclude [j0, j1] from both S1 and S2
    const double S1 = stieltjes_trap_history_forward(j1, j0-1, F, h, b, use_theta);
    const double S2 = stieltjes_trap_history_forward(j2, j0-1, F, h, b, use_theta);

    // --- 2) Abel weights (A = dt/z_eff = 2 in both rows) ---
    const double a1 = alpha_coeff(h,  0.5*h);
    const double b1 = beta_coeff (h,  0.5*h);
    const double g1 = gamma_coeff(h,  0.5*h);

    const double a2 = alpha_coeff(2*h, h);
    const double b2 = beta_coeff (2*h, h);
    const double g2 = gamma_coeff(2*h, h);

    // --- 3) kernel samples ---
    const double K1_0   = K(t1, t0,   b);
    const double K1_mid = K(t1, tmid, b);
    const double K1_1   = K(t1, t1,   b);

    const double K2_0 = K(t2, t0, b);
    const double K2_1 = K(t2, t1, b);
    const double K2_2 = K(t2, t2, b);

    const double G1 = G(t1, z, b);
    const double G2 = G(t2, z, b);

    double F1_pred, F2_pred;

    if (j0 >= 1) {
        // Interior: solve (1 - g1*K11 - (3/8)*b1*K1mid) * F1 = RHS
        const double den1_pred = 1.0 - g1*K1_1 - (3.0/8.0)*b1*K1_mid;
        const double rhs1_pred = G1 + S1
            + a1*K1_0*F[j0]
            + b1*K1_mid*( (-0.125)*F[j0-1] + 0.75*F[j0] );
        F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
    } else {
        // First block: linear midpoint Fmid ≈ 0.5(F0 + F1)
        const double den1_pred = 1.0 - g1*K1_1 - 0.5*b1*K1_mid;
        const double rhs1_pred = G1 + S1 + a1*K1_0*F[j0] + 0.5*b1*K1_mid*F[j0];
        F1_pred = rhs1_pred / ((std::abs(den1_pred) < FPM_EPSILON) ? FPM_EPSILON : den1_pred);
    }

    // Even row predictor is the usual 1×1 (uses predicted F1)
    {
        const double den2_pred = 1.0 - g2*K2_2;
        const double rhs2_pred = G2 + S2 + a2*K2_0*F[j0] + b2*K2_1*F1_pred;
        F2_pred = rhs2_pred / ((std::abs(den2_pred) < FPM_EPSILON) ? FPM_EPSILON : den2_pred);
    }

    // Odd equation:  A*F1 + B*F2 = C
    const double A = 1.0 - (3.0/4.0) * b1 * K1_mid - g1 * K1_1;
    const double B = (1.0/8.0) * b1 * K1_mid;
    const double C = G1 + S1 + a1 * K1_0 * F[j0]
                           + (3.0/8.0) * b1 * K1_mid * F[j0];

    // Even equation: D*F1 + E*F2 = R
    const double D = - b2 * K2_1;
    const double E = 1.0 - g2 * K2_2;
    const double R = G2 + S2 + a2 * K2_0 * F[j0];

    // Solve the 2x2 once. 
    const double det = A*E - B*D;
    const double inv_det = (std::abs(det) < FPM_EPSILON) ? (1.0/FPM_EPSILON) : (1.0/det);

    double F1_corr = ( E*C - B*R) * inv_det;
    double F2_corr = (-D*C + A*R) * inv_det;

    const double theta = 1; // set to e.g., 0.7 if you want under-relaxation
    F[j1] = theta*F1_corr + (1.0 - theta)*F1_pred;
    F[j2] = theta*F2_corr + (1.0 - theta)*F2_pred;

    }

    return F;
}

std::vector<double> solve_nu_f_block_with_abel(double t_max, double z, double b,
                                               int num_steps, bool use_theta,
                                               std::vector<double>& t_grid,
                                               double t_cut = 0.01)
{
    std::vector<double> F(num_steps+1, 0.0);
    // seeding uses the *same* grid we’re about to build
    t_grid.resize(num_steps+1);
    for (int i=0;i<=num_steps;++i) t_grid[i] = (t_max/num_steps) * i;
    const int k0 = seed_nu_f_on_grid(t_grid, use_theta, z, b, F, t_cut);
    return solve_nu_f_block_by_block(t_max, z, b, num_steps, use_theta, t_grid, F, 0);
}

// Helper functions for the forward pdf terms
inline double rho_tau(double tau, double tau_p) {
    const double s  = std::sqrt(2.0 * tau    + 1.0);
    const double sp = std::sqrt(2.0 * tau_p  + 1.0);
    const double denom = s + sp;
    if (denom <= FPM_EPSILON) return 0.0;
    return (s - sp) / denom;
}

inline double C_tau_smooth(double tau, double tau_p, double b_scaled) {
    const double r = rho_tau(tau, tau_p);
    const double b2 = b_scaled * b_scaled;
    return (1.0 - 2.0 * b2 * r) * safe_exp(-b2 * r);
}

// integral_{last panel} 2 * C(tau, tau-u^2) * (nu_lin(tau-u^2) - nu(tau)) / u^2 du  by trapezoid in u
inline double last_panel_u_linear(double tau, double a, double b_scaled,
                                  double nu_at_a, double nu_at_tau,
                                  int K_sub = 256) {
    const double u0 = std::sqrt(std::max(0.0, tau - a)); // upper u for the panel
    if (u0 <= 0.0 || K_sub <= 0) return 0.0;
    const double du = u0 / K_sub;

    auto nu_lin = [&](double s){ // linear on [a, tau]
        if (tau <= a) return nu_at_a;
        const double w = (s - a) / (tau - a);
        return (1.0 - w) * nu_at_a + w * nu_at_tau;
    };

    auto f = [&](double u)->double {
        if (u <= 0.0) return 0.0;   // finite limit at 0, measure-zero endpoint
        const double s = tau - u*u;
        const double C = C_tau_smooth(tau, s, b_scaled);
        return 2.0 * C * (nu_lin(s) - nu_at_tau) / (u*u);
    };

    double acc = 0.5 * f(0.0);
    for (int k = 1; k < K_sub; ++k) acc += f(k * du);
    acc += 0.5 * f(u0);
    return acc * du;
}

// Closed-form last panel using nu′(tau) and nu''(tau) + linearized C
// H(x) ≈ [-A x + (0.5 B + alpha A) x^2],  alpha = (3/2) b^2 / S^2,  x = tau - s
inline double last_panel_quadratic_closed(double tau, double a, double b_scaled,
                                          const std::vector<double> &tau_grid,
                                          const std::vector<double> &nu_f_vals,
                                          double nu_prime, double nu_second) {
    const double h = tau - a;
    if (h <= 0.0){
      return 0.0;
    }
    // a) Analytically integrate the dominant singular part: integral[-nu'(tau)/sqrt(tau-s)]ds
    const double S = std::sqrt(2.0 * tau + 1.0);
    const double alpha = 1.5 * b_scaled * b_scaled / (S*S);
    const double beta = 0.5*nu_second + alpha * nu_prime;
        // a) Analytically integrate the dominant parts: integral[ H'(tau)(s-tau) + H''(tau)/2 * (s-tau)^2 ] / (tau-s)^1.5 ds
    // This simplifies to: integral[ -nu'(tau)/sqrt(tau-s) + H''(tau)/2 * sqrt(tau-s) ] ds
    const double analytic_part_1 = -2.0 * nu_prime * std::sqrt(h);
    const double analytic_part_2 = (2.0/3.0) * beta * std::pow(h, 1.5);
    
   return(analytic_part_1 + analytic_part_2);

}

// Richardson on the last panel: I*(h) ≈ I(h/2) + (I(h/2)-I(h)) / (1/sqrt(2) - 1)
inline double last_panel_richardson(double tau, double a, double b_scaled,
                                    double nu_at_a, double nu_at_tau,
                                    // nu at midpoint from a quadratic fit (or pass any estimate)
                                    double nu_at_mid_from_quad,
                                    int K_sub = 128) {
    // I(h): one panel
    const double I_full = last_panel_u_linear(tau, a, b_scaled, nu_at_a, nu_at_tau, K_sub);

    // split at midpoint (in τ), but use ν(mid) from quadratic fit
    const double mid = 0.5 * (a + tau);
    const double I_half =
        last_panel_u_linear(tau, mid, b_scaled, nu_at_mid_from_quad, nu_at_tau, K_sub) +
        last_panel_u_linear(tau, a,   b_scaled, nu_at_a,              nu_at_mid_from_quad, K_sub);

    const double denom = (1.0 / std::sqrt(2.0)) - 1.0; // negative
    return I_half + (I_half - I_full) / denom;
}

inline double quad_interp_theta(
    const std::vector<double>& theta,
    const std::vector<double>& nu,
    int i,            // panel index for [i, i+1]
    double theta_q    // query inside [theta[i], theta[i+1]]
) {
    const int N = static_cast<int>(theta.size()) - 1;
    // Choose a local triple around the panel:
    // interior panels: (i-1,i,i+1); near left: (0,1,2); near right: (N-2,N-1,N)
    int a, b, c;
    if (i == 0) { a = 0; b = 1; c = 2; }
    else if (i == N-1) { a = N-2; b = N-1; c = N; }
    else { a = i-1; b = i; c = i+1; }

    const double xa = theta[a], xb = theta[b], xc = theta[c];
    const double ya = nu[a],    yb = nu[b],    yc = nu[c];

    // Lagrange basis at θq
    const double Lab = (theta_q - xb) / (xa - xb);
    const double Lac = (theta_q - xc) / (xa - xc);
    const double Lba = (theta_q - xa) / (xb - xa);
    const double Lbc = (theta_q - xc) / (xb - xc);
    const double Lca = (theta_q - xa) / (xc - xa);
    const double Lcb = (theta_q - xb) / (xc - xb);

    return ya * Lab * Lac + yb * Lba * Lbc + yc * Lca * Lcb;
}

/**
 * @brief The smooth integrand of the final PDF integral after the u-transform.
 *        f(u) = 2 * H_smooth(θ, θ-u²) * (v(θ-u²) - v(θ)) * (1+θ-u²) / u²
 *        This function computes the value of f(u) given theta and theta'.
 *
 * @param theta The final theta value (endpoint of integration).
 * @param theta_p The point theta' < theta.
 * @param nu_f_theta The value of v^f at theta.
 * @param nu_f_p The value of v^f at theta'.
 * @param b_scaled The scaled boundary level.
 * @return The value of the regularized integrand.
 */
inline double regularized_pdf_integrand_theta(
    double theta, double theta_p,
    double nu_f_theta, double nu_f_p,
    double b_scaled)
{
    const double u2 = theta - theta_p;
    if (u2 <= 0) return 0.0;

    // Kernel's smooth part H_smooth(theta, theta_p)
    const double sum_part = 2.0 + theta + theta_p;
    const double ratio = u2 / sum_part;
    const double b2 = b_scaled * b_scaled;
    const double H_smooth = (1.0 - 2.0 * b2 * ratio) * safe_exp(-b2 * ratio);

    const double nu_diff = nu_f_p - nu_f_theta;
    const double jacobian_part = 1.0 + theta_p;

    const double denom = u2*sum_part * sqrt_pos(sum_part);
    if (denom < FPM_EPSILON) return 0.0;
    // NB 2* is from the jacobian of the u-transform
    return 2.0* (H_smooth * nu_f_p - nu_f_theta) * jacobian_part / denom;
}

/**
 * @brief function to integrate the PDF's singular term using the u-transform
 *        directly on the solver's native uniform theta-grid.
 */
// [[Rcpp::export]]
double integrate_pdf_forward_theta_u(
    double theta_max,
    const std::vector<double>& theta_grid,
    const std::vector<double>& nu_f_vals,
    double b_scaled)
{
    const int N = theta_grid.size() - 1;
    if (N < 1) return 0.0;

    const double nu_f_at_max = nu_f_vals.back();

    std::vector<double> u_grid(N + 1);
    std::vector<double> f_vals(N + 1);

    for (int i = 0; i <= N; ++i) {
        u_grid[i] = sqrt_pos(theta_max - theta_grid[i]);
    }

    // Evaluate f(u_i) for i = 0 to N-1
    for (int i = 0; i < N; ++i) {
        f_vals[i] = regularized_pdf_integrand_theta(
            theta_max, theta_grid[i],
            nu_f_at_max, nu_f_vals[i],
            b_scaled
        );
    }

    // Evaluate f(u_N) = f(0) using the CORRECTED derivative limit
    const double mid = 0.5 * (theta_grid[N-1] + theta_max);
	const RightQuad rq = right_end_quadratic(theta_grid, nu_f_vals, mid);
    const double nu_prime_theta = rq.nu_prime;
    const double s_max = 1.0 + theta_max;
    const double f_right_limit = -nu_prime_theta / std::sqrt(2.0 * s_max);

    double integral_sum = 0.0;
    // Panel-wise Simpson with quadratic ν(θ′) at θ_mid
    // Panels are [i, i+1] in θ (equivalently [u_{i+1}, u_i] in u)
    for (int i = N - 1; i >= 0; --i) {
        const double uL = u_grid[i];
        const double uR = u_grid[i+1];
        const double h  = uL - uR; // panel width in u

        const double thetaL = theta_grid[i];
        const double thetaR = theta_grid[i+1];
        const double uM     = 0.5 * (uL + uR);
        const double thetaU = theta_max - uM * uM;

        // ν values
        const double nuL = nu_f_vals[i];
        const double nuR = nu_f_vals[i+1];
        const double nuU = quad_interp_theta(theta_grid, nu_f_vals, i, thetaU);

        // f at left endpoint
        const double fL = regularized_pdf_integrand_theta(theta_max, thetaL,
                                                      nu_f_at_max, nuL, b_scaled);

        // f at right endpoint: exact limit on the final panel, standard elsewhere
        const double fR = (i == N - 1) ? f_right_limit :  regularized_pdf_integrand_theta(theta_max, thetaR,
                                          nu_f_at_max, nuR, b_scaled);

        // f at midpoint with quadratic ν and midpoint θ′ in the smooth factors
        const double fM = regularized_pdf_integrand_theta(theta_max, thetaU,
                                                      nu_f_at_max, nuU, b_scaled);

        // Simpson on this u-panel
        integral_sum += (h / 6.0) * (fL + 4.0 * fM + fR);
    }

    return integral_sum;
}

// Integrate integral_0^tau H(s) / (tau - s)^{3/2} ds by linear product integration.
// On each panel [a,b], replace H by its linear interpolant and integrate
// analytically to obtain weights W_a, W_b. The last panel uses H(b)=0 
// W_a = 2sqrt (tau-a)/h, W_b = 0, which is the exact limit.
// [[Rcpp::export]]
double integrate_pdf_forward(
    double tau, const std::vector<double>& tau_grid,
    const std::vector<double>& nu_f_vals,
    double b_scaled,
    std::string lp,
    int K_sub = 256
	) {

    double integral_sum = 0.0;
    const int N = static_cast<int>(tau_grid.size()) - 1;
	double nu_f_at_tau = nu_f_vals.back();
    for (int i = 0; i < N; ++i) {
        const double a = tau_grid[i];
        const double b = (i == N-1) ? tau : tau_grid[i+1];
        const double h = b - a;
        if (b <= a) continue;

        // Evaluate the smooth numerator function H at the interval endpoints
        double H_a = C_tau_smooth(tau, a, b_scaled) * (nu_f_vals[i])-nu_f_at_tau;
        double H_b = (i == N-1) ? 0.0
                                      : C_tau_smooth(tau, b, b_scaled)*(nu_f_vals[i+1])-nu_f_at_tau;
		const double y_a = tau - a;
		const double y_b = tau - b;
		const double r_a = std::sqrt(std::max(y_a, 0.0));
		const double r_b = std::sqrt(std::max(y_b, 0.0));
        // Analytically derived weights for integral[a,b] (LinearInterp(H)) / (tau-s)^1.5 ds
        double W_a = 0.0, W_b = 0.0;
        if (i == N - 1) { // Final interval where b ~ tau, so H(b) ~ 0
			if (lp == "None") {
				// plain product rule with Hb=0 and last-panel weights
				weights_panel_32(tau, a, tau, W_a, W_b); // Wb=0 here
				integral_sum += W_a * H_a;
			}
			if (lp == "USplit") {
				integral_sum += last_panel_u_linear(tau, a, b_scaled, nu_f_vals[i], nu_f_at_tau, K_sub);
			}
			
			if (lp == "QuadraticClosed" || lp == "Richardson") {
				const double mid = 0.5 * (a + tau);
				const RightQuad rq = right_end_quadratic(tau_grid, nu_f_vals, mid);
			
				if (lp == "QuadraticClosed") {
					integral_sum += last_panel_quadratic_closed(tau, a, b_scaled, tau_grid, nu_f_vals,rq.nu_prime, rq.nu_second);
				} else { // Richardson
					integral_sum += last_panel_richardson(tau, a, b_scaled, nu_f_vals[i], nu_f_at_tau, rq.nu_at_mid, K_sub);
				}
			}
        } else {
			weights_panel_32(tau, a, b, W_a, W_b);
			integral_sum += W_a * H_a + W_b * H_b;
        }      
    }
    
    return integral_sum;
}

// [[Rcpp::export]]
inline double integrate_pdf_forward_u_robust(double tau,
    const std::vector<double>& tau_grid,
	const std::vector<double>& nu_f_vals,
    double b_scaled)
{

    const int N = tau_grid.size()-1;
    if (N < 2) return 0.0;
    const double mid = 0.5 * (tau_grid[N-1] + tau);
	const RightQuad rq = right_end_quadratic(tau_grid, nu_f_vals, mid);
    const double nu_f_at_max = nu_f_vals.back();
    std::vector<double> u_grid(N + 1);
    std::vector<double> f_vals(N + 1);
    for (int i = 0; i <= N; ++i) {
        u_grid[i] = sqrt_pos(tau - tau_grid[i]);
    }

    // Evaluate f(u_i) for i = 0 to N-1
    for (int i = 0; i < N; ++i) {
      double Hsmooth = C_tau_smooth(tau, tau_grid[i], b_scaled);
      double u2 = (u_grid[i] * u_grid[i]);
      f_vals[i] = 2.0 * (Hsmooth * nu_f_vals[i] - nu_f_at_max) / u2; 
    }
    f_vals[N] = -2.0 * rq.nu_prime;
                       
    // Perform the non-uniform trapezoidal rule sum
    double integral_sum = 0.0;
    for (int i = 0; i <N-1; ++i) {
        const double du = u_grid[i] - u_grid[i+1];
        integral_sum += 0.5 * (f_vals[i] + f_vals[i + 1]) * du;
    }
    const double kappa = (3.0 * b_scaled * b_scaled) / (2.0 * (2.0 * tau + 1.0));
    const double a0 = -2.0 * rq.nu_prime;
    const double a2 = rq.nu_second - 2.0 * kappa * rq.nu_prime;
    // Last panel: exact quadratic endpoint model on [u_N, u_{N-1}]
    const double uL = u_grid[N];
    const double uR = u_grid[N-1];
    integral_sum += a0 * (uR - uL) + (a2 / 3.0) * (uR * uR * uR - uL * uL * uL);
    return integral_sum;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SECTION 4: Compilation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// [[Rcpp::export]]
double ou_fht_cdf(double t, double lambda, double theta, double sigma, double z0, double b_level, int num_steps) {
    ScaledParams p = prepare_ou_params(t, lambda, theta, sigma, z0, b_level);
    if (p.is_trivial) return p.cdf;

    const double v = 1.0 - std::exp(-p.t_scaled);
    const int N = (num_steps < 2) ? 2 : (num_steps % 2 ? num_steps + 1 : num_steps);

    // Block solver requires a uniform grid for its internal weights
    const double h_v = v / N;
    std::vector<double> v_grid(N + 1);
    for(int j = 0; j <= N; ++j) v_grid[j] = j * h_v;

	auto nu_vals = solve_nu_b_block_with_abel(v, p.z_scaled, p.b_scaled, N, v_grid);
    
    double integral_sum = integrate_cdf_backward_trapezoid(v, v_grid, nu_vals, p.z_scaled, p.b_scaled);
    
    // Formula for G(t,z) has a coefficient of 2
    return std::max(0.0, std::min(1.0, 2*integral_sum));
}

// [[Rcpp::export]]
double ou_fht_pdf_forward(double t, double lambda, double theta, double sigma, double z0, double b_level, int num_steps) {
    ScaledParams p = prepare_ou_params(t, lambda, theta, sigma, z0, b_level);
    if (p.is_trivial) return p.pdf;
	const double et  = std::exp(p.t_scaled);
	const double e2t = std::exp(2.0 * p.t_scaled);
	const double e2tm1 = e2t - 1.0;
	const double tau = 0.5 * (e2t - 1.0);
    const int N = (num_steps < 2) ? 2 : num_steps;

	std::vector<double> theta_grid_uniform, tau_grid_non_uniform, tau_grid_block, tau_grid_uniform, tau_grid_right_clustered, theta_grid_block;
	double theta_max;
	make_theta_and_tau_grids(p.t_scaled, N, tau, theta_grid_uniform, tau_grid_non_uniform, tau_grid_uniform, tau_grid_right_clustered,theta_max);
    tau_grid_non_uniform.back() = tau; // force last entry to be exactly tau_max
    std::vector<double> nu_vals(N + 1, 0.0);

	auto nu_f_vals = solve_nu_f_block_with_abel(theta_max, p.z_scaled, p.b_scaled, N, true, theta_grid_block);

	const double numA  = et * p.b_scaled - p.z_scaled;
    const double termA = - (numA * safe_exp(-((numA*numA)/e2tm1) + 2.0*p.t_scaled)
                     / std::sqrt(M_PI * e2tm1*e2tm1*e2tm1));

	const double termB = ( et*p.b_scaled + e2t/(std::sqrt(M_PI*e2tm1)) ) * nu_f_vals.back();

	const double prefactor = (e2t / std::sqrt(M_PI));
    
    const double integral_sum =  integrate_pdf_forward_theta_u(theta_max, theta_grid_block, nu_f_vals, p.b_scaled);
    
    const double termC = prefactor * integral_sum;
    
    const double g_scaled = termA - termB + termC;
			
    return lambda * g_scaled;
}

// [[Rcpp::export]]
Rcpp::List ou_fht_pdf_forward_debug(double t, double lambda, double theta, double sigma, double z0, double b_level, int num_steps) {
    ScaledParams p = prepare_ou_params(t, lambda, theta, sigma, z0, b_level);
    if (p.is_trivial){
		return Rcpp::List::create(
            Rcpp::Named("pdf") = p.pdf,
            Rcpp::Named("termA") = 0.0,
            Rcpp::Named("termB") = 0.0,
            Rcpp::Named("termC") = 0.0,
			Rcpp::Named("C_panel") = 0.0,
			Rcpp::Named("C_u") = 0.0,
            Rcpp::Named("nu_f_at_tau") = 0.0,
            Rcpp::Named("integral_val") = 0.0
        );
	}
	const double et  = std::exp(p.t_scaled);
	const double e2t = std::exp(2.0 * p.t_scaled);
	const double e2tm1 = e2t - 1.0;
	if (e2tm1 < 1e-12) {
        return Rcpp::List::create(
            Rcpp::Named("pdf") = 0.0,
            Rcpp::Named("termA") = 0.0,
            Rcpp::Named("termB") = 0.0,
            Rcpp::Named("termC") = 0.0,
			Rcpp::Named("C_panel") = 0.0,
			Rcpp::Named("C_u") = 0.0,
            Rcpp::Named("nu_f_at_tau") = 0.0,
            Rcpp::Named("integral_val") = 0.0
        );
    }

    const double tau = 0.5 * (e2t - 1.0);
    const int N = (num_steps < 2) ? 2 : num_steps;

	std::vector<double> theta_grid_uniform, tau_grid_non_uniform, tau_grid_block, tau_grid_uniform, tau_grid_right_clustered, theta_grid_block;
	double theta_max;
	make_theta_and_tau_grids(p.t_scaled, N, tau, theta_grid_uniform, tau_grid_non_uniform, tau_grid_uniform, tau_grid_right_clustered,theta_max);
    tau_grid_non_uniform.back() = tau; // force last entry to be exactly tau_max

	auto nu_f_vals_tau = solve_nu_f_block_with_abel(tau, p.z_scaled, p.b_scaled, N, false, tau_grid_block);
    auto nu_f_vals_theta_block = solve_nu_f_block_with_abel(theta_max, p.z_scaled, p.b_scaled, N, true, theta_grid_block);
    auto nu_f_vals_theta_grid = solve_nu_f_grid_with_abel(theta_grid_uniform, p.z_scaled, p.b_scaled, true);
	auto nu_f_vals_tau_grid = solve_nu_f_grid_with_abel(tau_grid_non_uniform,p.z_scaled, p.b_scaled, false);
	auto nu_f_vals_tau_grid_uniform = solve_nu_f_grid_with_abel(tau_grid_uniform,p.z_scaled, p.b_scaled, false);
    auto nu_f_vals_tau_grid_right = solve_nu_f_grid_with_abel(tau_grid_right_clustered,p.z_scaled, p.b_scaled, false);
    
    std::vector<double>tau_from_theta(theta_grid_block.size()), t_vals(theta_grid_block.size());
	for (size_t i=0; i<theta_grid_block.size(); ++i) {
          tau_from_theta[i] = tau_of_theta(theta_grid_block[i]);
          t_vals[i] = t_from_tau_theta(theta_grid_block[i], true);
    }
        
	const double numA  = et * p.b_scaled - p.z_scaled;
    const double termA = - (numA * safe_exp(-((numA*numA)/e2tm1) + 2.0*p.t_scaled)
                     / std::sqrt(M_PI * e2tm1*e2tm1*e2tm1));

	const double termB = ( et*p.b_scaled + e2t/(std::sqrt(M_PI*e2tm1)) );

    const double A_trap =  integrate_pdf_forward(tau, tau_grid_block, nu_f_vals_tau, p.b_scaled, "Richardson");
	const double B_trap =  integrate_pdf_forward(tau, tau_grid_non_uniform, nu_f_vals_tau_grid, p.b_scaled, "Richardson");
	const double C_trap =  integrate_pdf_forward(tau, tau_grid_right_clustered, nu_f_vals_tau_grid_right, p.b_scaled, "Richardson");
	const double A_u =  integrate_pdf_forward_u_robust(tau, tau_grid_block, nu_f_vals_tau, p.b_scaled);
	const double B_u =  integrate_pdf_forward_u_robust(tau, tau_grid_non_uniform, nu_f_vals_tau_grid, p.b_scaled);
	const double C_u =  integrate_pdf_forward_u_robust(tau, tau_grid_right_clustered, nu_f_vals_tau_grid_right, p.b_scaled);
    const double Theta_grid =  integrate_pdf_forward_theta_u(theta_max, theta_grid_uniform, nu_f_vals_theta_grid, p.b_scaled);
    const double Theta_block =  integrate_pdf_forward_theta_u(theta_max, theta_grid_block, nu_f_vals_theta_block, p.b_scaled);
	const double prefactor = (e2t / std::sqrt(8.0 * M_PI));
	const double prefactor_theta = (e2t / std::sqrt(M_PI));

	return Rcpp::List::create(
		Rcpp::Named("et") = et,
        Rcpp::Named("e2t") = e2t,
		Rcpp::Named("numA") = numA,
        Rcpp::Named("e2tm1") = e2tm1,
        Rcpp::Named("termA") = termA,
        Rcpp::Named("termB") = termB,
		Rcpp::Named("A_trap") = A_trap,
		Rcpp::Named("B_trap") = B_trap,
		Rcpp::Named("C_trap") = C_trap,
        Rcpp::Named("Theta_grid") = Theta_grid,
        Rcpp::Named("Theta_block") = Theta_block,
		Rcpp::Named("A_u") = A_u,
		Rcpp::Named("B_u") = B_u,
		Rcpp::Named("C_u") = C_u,
		Rcpp::Named("prefactor") = prefactor,
        Rcpp::Named("prefactor_theta") = prefactor_theta,
		Rcpp::Named("tau") = tau,
        Rcpp::Named("t_scaled") = p.t_scaled,
        Rcpp::Named("z_scaled") = p.z_scaled,
        Rcpp::Named("b_scaled") = p.b_scaled,
		Rcpp::Named("lambda") = lambda,
        Rcpp::Named("theta") = theta_max,
        Rcpp::Named("t_vals") = t_vals,
		Rcpp::Named("nu_f_vals_theta_block") = nu_f_vals_theta_block,
        Rcpp::Named("nu_f_vals_theta_grid") = nu_f_vals_theta_grid,
		Rcpp::Named("nu_f_vals_tau_grid") = nu_f_vals_tau_grid,
        Rcpp::Named("nu_f_vals_tau_block") = nu_f_vals_tau,
        Rcpp::Named("tau_grid_right_clustered") = tau_grid_right_clustered,
        Rcpp::Named("tau_grid_uniform") = tau_grid_uniform,
        Rcpp::Named("nu_f_vals_tau_grid_uniform") = nu_f_vals_tau_grid_uniform,
		Rcpp::Named("nu_f_vals_tau_grid_right") = nu_f_vals_tau_grid_right,
		Rcpp::Named("tau_vals_grid") = tau_grid_non_uniform,
        Rcpp::Named("theta_grid_block") = theta_grid_block,
        Rcpp::Named("theta_grid_uniform") = theta_grid_uniform
    );
}


// [[Rcpp::export]]
NumericVector ou_fht_pdf_forward_vec(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b_level, int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_pdf_forward(ti, lambda, theta, sigma, z0, b_level, num_steps) : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector ou_fht_cdf_vec(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b_level, int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_cdf(ti, lambda, theta, sigma, z0, b_level, num_steps) : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector simulate_ou_hit_times_std(int n,
                                        double lambda, double theta, double sigma,
                                        double z0, double b_level,
                                        double dt = 1e-3, double t_max = 10.0) {
  // 1. PARAMETER VALIDATION (same as original)
  if (lambda <= 0.0 || sigma <= 0.0 || dt <= 0.0 || t_max <= 0.0) {
    stop("simulate_ou_hit_times_std: invalid parameters.");
  }
  ScaledParams p = prepare_ou_params(t_max, lambda, theta, sigma, z0, b_level);
  const double lambda_scaled = 1, theta_scaled=0, sigma_scaled=1;
  // 3. SCALED TIME PARAMETERS
  // The simulation now runs in scaled time s = lambda * t
  const double ds = lambda_scaled * dt;        // Scaled time step
  const double s_max = lambda_scaled * t_max;  // Scaled max time
  const int steps = std::max(1, (int)std::ceil(s_max / ds)); // Same as std::ceil(t_max / dt)

  // 4. SIMULATION OF STANDARDIZED PROCESS
  // The standardized OU process is dZ = -Z ds + sqrt(2) dW'
  // Or, in the form your original code uses: dZ = -lambda_std * (Z - theta_std) ds + sigma_std dW
  // with lambda_std = 1.0, theta_std = 0.0, sigma_std = sqrt(2.0)
  // Let's use the exact transition from your original code with these new parameters.
  // The standard OU process used in papers is often dZ = -Z ds + dW, which means sigma_std=1.
  // We'll stick to that convention for simplicity, as it matches the scaling properties.
  // SDE: dZ = -Z ds + dW_s. This has lambda_std=1, theta_std=0, sigma_std=1.
    // The per-step variance is (sigma_std^2 / (2*lambda_std)) * (1 - e^(-2*lambda_std*ds))
  // Here: ( (sqrt(2))^2 / (2*1) ) * (1 - e^(-2*ds)) = (2/2) * (1 - e*e) = 1.0 * (1 - e*e)
  // Exact per-step coefficients for the *standardized* process
  const double e = std::exp(-lambda_scaled * ds); // e = exp(-1.0 * lambda * dt)
  const double var_step = (1.0 - e * e) / (2.0 * lambda_scaled);
  const double sd_step = std::sqrt(var_step);

  NumericVector tau(n);

  for (int j = 0; j < n; ++j) {
    double z = p.z_scaled; // Start with the scaled initial value
    double s = 0.0;      // Time variable is now the scaled time 's'
    bool hit = false;

    // Edge case: already at or beyond barrier at t=0
    if ((p.z_scaled - p.b_scaled) <= FPM_EPSILON) {
      tau[j] = 0.0;
      continue;
    }

    for (int k = 0; k < steps; ++k) {
      const double z_prev = z;
      
      // Exact transition for standardized process (theta_std = 0.0)
      z = z_prev * e + sd_step * R::rnorm(0.0, 1.0);
      const double s_next = s + ds;

      // Crossing detection for hitting from above.
      // The product check is robust and works for both directions.
      if ((z_prev - p.b_scaled) * (z - p.b_scaled) <= 0.0 && z_prev > p.b_scaled) {
        // Linear interpolation within the step
        const double dz = z - z_prev;
        double alpha = 0.5; // Default to midpoint if dz is zero
        if (std::abs(dz) > 0.0) {
          alpha = (p.b_scaled - z_prev) / dz;
          alpha = std::min(1.0, std::max(0.0, alpha));
        }
        
        // 5. UN-SCALE THE FINAL HIT TIME
        // The hit happened at scaled time s_hit = s + alpha * ds
        // We must convert it back to the original time domain.
        tau[j] = (s + alpha * ds) / lambda;
        hit = true;
        break;
      }

      s = s_next;
    }

    if (!hit) tau[j] = NA_REAL;
  }

  return tau;
}

inline double clamp_prob(double p) {
    // Tight but safe; keeps logit finite without biasing interior.
    const double eps = 1e-15;
    if (p <= eps) return eps;
    if (p >= 1.0 - eps) return 1.0 - eps;
    return p;
}

inline double logit_from_p(double p) {
    // log(p) - log1p(-p) is stable near 1
    p = clamp_prob(p);
    return std::log(p) - std::log1p(-p);
}

// Choose an adaptive step in scaled time ttbar = lambda * t.
// Shrinks with (1 - v) where v = 1 - exp(-ttbar), i.e., as you approach the endpoint.
inline double choose_tbar_step(double tbar) {
    const double v = 1.0 - std::exp(-tbar);
    const double s = std::max(1e-8, 1.0 - v);     // distance to endpoint in v-space
    // Bounds are conservative; you can tune if you like.
    const double h_min = 1e-5;
    const double h_max = 2.5e-2;
    double h = 0.35 * s;                          // scale with endpoint distance
    if (h < h_min) h = h_min;
    if (h > h_max) h = h_max;
    return h;
}

// Central 5-point derivative: f'(x0) ~ (-f_{+2}+8 f_{+1} - 8 f_{-1} + f_{-2})/(12 h)
template <typename T>
inline long double deriv5_central(long double fm2, long double fm1,
                                  long double fp1, long double fp2,
                                  long double h) {
    return (-fp2 + 8.0L*fp1 - 8.0L*fm1 + fm2) / (12.0L * h);
}

// Forward 5-point derivative for boundary: f'(x0) ~ (-25 f0 + 48 f1 - 36 f2 + 16 f3 - 3 f4)/(12 h)
template <typename T>
inline long double deriv5_forward(long double f0, long double f1,
                                  long double f2, long double f3,
                                  long double f4, long double h) {
    return (-25.0L*f0 + 48.0L*f1 - 36.0L*f2 + 16.0L*f3 - 3.0L*f4) / (12.0L * h);
}

// Thin wrapper to evaluate your existing CDF at physical time
inline double G_at(double t, double lambda, double theta, double sigma,
                   double z0, double b_level, int num_steps) {
    // Calls your already-correct backward CDF:
    return ou_fht_cdf(t, lambda, theta, sigma, z0, b_level, num_steps);
}

// [[Rcpp::export]]
double ou_fht_pdf_from_cdf_regularized(double t,
                                       double lambda, double theta, double sigma,
                                       double z0, double b_level,
                                       int num_steps=200)
{
    // Trivial/degenerate very small times: just fall back to central with tiny step
    if (t <= 0.0) return 0.0;

    const double tbar = lambda * t;
    const double hbar = choose_tbar_step(tbar);
    const double h    = hbar / lambda;  // physical-time step matching a step in tbar

    // Evaluate G at the needed offsets; prefer central stencil if we can go both sides
    const bool central_ok = (t > 2.0*h);

    long double Hprime_tbar = 0.0L;  // derivative of logit(G) w.r.t. tbar

    if (central_ok) {
        // central 5-point around t
        const double Gm2 = G_at(t - 2*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double Gm1 = G_at(t - 1*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double Gp1 = G_at(t + 1*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double Gp2 = G_at(t + 2*h, lambda, theta, sigma, z0, b_level, num_steps);

        const long double Hm2 = (long double)logit_from_p(Gm2);
        const long double Hm1 = (long double)logit_from_p(Gm1);
        const long double Hp1 = (long double)logit_from_p(Gp1);
        const long double Hp2 = (long double)logit_from_p(Gp2);

        Hprime_tbar = deriv5_central<long double>(Hm2, Hm1, Hp1, Hp2, (long double)hbar);
    } else {
        // too close to 0: forward 5-point at t
        const double G0 = G_at(t + 0*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double G1 = G_at(t + 1*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double G2 = G_at(t + 2*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double G3 = G_at(t + 3*h, lambda, theta, sigma, z0, b_level, num_steps);
        const double G4 = G_at(t + 4*h, lambda, theta, sigma, z0, b_level, num_steps);

        const long double H0 = (long double)logit_from_p(G0);
        const long double H1 = (long double)logit_from_p(G1);
        const long double H2 = (long double)logit_from_p(G2);
        const long double H3 = (long double)logit_from_p(G3);
        const long double H4 = (long double)logit_from_p(G4);

        Hprime_tbar = deriv5_forward<long double>(H0, H1, H2, H3, H4, (long double)hbar);
    }

    // Chain back: g(t) = lambda * G*(1-G) * d/dttbar logit(G)
    const double G0_center = G_at(t, lambda, theta, sigma, z0, b_level, num_steps);
    const double Gc = clamp_prob(G0_center);
    long double g = (long double)lambda * (long double)Gc * (1.0L - (long double)Gc) * Hprime_tbar;

    // Non-negativity guard (rare tiny negatives from numerical noise)
    if (g < 0.0L && g > -1e-12L) g = 0.0L;

    return (double)g;
}

// [[Rcpp::export]]
NumericVector ou_fht_pdf_from_cdf_regularized_vec(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b_level, int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_pdf_from_cdf_regularized(ti, lambda, theta, sigma, z0, b_level, num_steps) : NA_REAL;
  }
  return out;
}

#endif // OU_HITTING_TIME_H


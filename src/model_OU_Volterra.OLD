#ifndef OU_HITTING_TIME_H
#define OU_HITTING_TIME_H

#define _USE_MATH_DEFINES
#include "utility_functions.h"
#include "utils_reducible_diffusion.h"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <functional>
#include <limits>
#include <cstddef>
#include <Rcpp.h>
#include <quadmath.h>
using namespace Rcpp;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Section 1: Core Constants and Helper Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Centralized function for parameter validation, scaling, and handling edge cases.
inline RD_Params prepare_ou_params(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double pow) {
  RD_Params p;
    p.omega = 1.0; // default: hitting from above
    p.b0 = b0;
    p.c = std::sqrt(lambda) / sigma;
    p.t_scaled = lambda * t;
    p.z_scaled = p.c * (z0 - theta);
    p.b_scaled = p.c * (b0 - theta);
    p.binf = binf;
    p.fixed_b = true;
    p.mu = 0.0;
    p.sigma = sigma;
    p.z0 = z0;
    if (std::abs(p.binf-p.b0) > FPM_EPSILON) {
        p.fixed_b = false;
    }
    p.lambda = lambda;
    p.theta=theta;
    p.tau = tau;
    p.pow = pow;
    if (p.z_scaled < p.b_scaled) { // Hitting from below -> flip problem
        p.omega = -1.0;
        p.z_scaled = -p.z_scaled;
        p.b_scaled = -p.b_scaled;
        p.z0 = -p.z0;
    }
    return p;
}

double v_to_t(double v) {
    if (v >= 1.0) return std::numeric_limits<double>::infinity();
    return -std::log(1.0 - v);
}

double theta_to_t(double theta) {
    if (theta <= -1.0) return std::numeric_limits<double>::infinity();
    return std::log1p(theta);
}

inline double beta_from_v_raw(double v, const RD_Params& pars) {
    const double scale = std::max(0.0, 1.0 - v);
    if (scale <= 0.0) {
        return 0.0;
    }
    if (pars.fixed_b) {
        return scale * pars.b_scaled;
    }
    const double t = v_to_t(v) / pars.lambda;
    const double t_max = pars.t_scaled / pars.lambda;
    const double bt = exp_decay_scalar(t_max - t, pars.b0, pars.binf, pars.tau, pars.pow);
    const double bt_scaled = pars.omega * (pars.c * (bt - pars.theta));
    return scale * bt_scaled;
}

inline double beta_from_theta_raw(double theta, const RD_Params& pars) {
    const double scale = std::max(0.0, 1.0 + theta);
    if (scale <= 0.0) {
        return 0.0;
    }
    if (pars.fixed_b) {
        return scale * pars.b_scaled;
    }
    const double t = theta_to_t(theta) / pars.lambda;
    const double bt = exp_decay_scalar(t, pars.b0, pars.binf, pars.tau, pars.pow);
    const double bt_scaled = pars.omega * (pars.c * (bt - pars.theta));
    return scale * bt_scaled;
}



inline double beta_from_v(double v, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    if (pars.fixed_b) {
        const double scale = std::max(0.0, 1.0 - v);
        return scale * pars.b_scaled;
    }
    if (cache && !cache->empty()) {
        const double cached = cache->lookup(v);
        if (std::isfinite(cached)) {
            return cached;
        }
    }
    return beta_from_v_raw(v, pars);
}

inline double beta_from_theta(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    if (pars.fixed_b) {
        const double scale = std::max(0.0, 1.0 + theta);
        return scale * pars.b_scaled;
    }
    if (cache && !cache->empty()) {
        const double cached = cache->lookup(theta);
        if (std::isfinite(cached)) {
            return cached;
        }
    }
    return beta_from_theta_raw(theta, pars);
}

inline void effective_positions_tv(double v, double vp, const RD_Params& pars,
                                double& Delta, double& one_minus_vp,
                                const BoundaryDecayCache* cache = nullptr) {
  const double z_eff = pars.z_scaled * std::max(0.0, 1.0 - v);
  double beta_vp;
  if (pars.fixed_b) {
        beta_vp = pars.b_scaled * (1.0 - vp);
  } else {
        beta_vp = beta_from_v(vp, pars, cache);
  }
  Delta = z_eff - beta_vp;
  one_minus_vp = std::max(0.0, 1.0 - vp);
}

// KERNELS
// The non-singular part of the backwards Volterra kernel in v from page 11 Lipton & Kaushansky(2018) 
double kernel_backward(double v, double vp, const RD_Params& pars) {
    const double b = pars.b_scaled;
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

// Time-varying extension of above. replace b with (beta(v) - beta(v'))/dv
inline double kernel_backward_tv_core(double v, double vp, double beta_v, double beta_vp) {
    const double dv = v - vp;
    if (dv <= FPM_EPSILON) return 0.0;
    const double psi = beta_v - beta_vp;
    if (std::abs(psi) < 1e-15) return 0.0;
    const double s = 2.0 - v - vp;
    if (s <= FPM_EPSILON) return 0.0;

    const double m_eff = -psi / dv;
    const double xi = safe_exp(-(psi * psi) / (dv * s));

    const double term1 = (2.0 * m_eff / std::sqrt(M_PI));
    const double term2 = ((1.0 - vp) * xi) / (s * std::sqrt(s));
    return term1 * term2;
}

inline double kernel_backward_tv_cached(double v, double vp, const RD_Params& pars,
                                        const BoundaryDecayCache& cache) {
    if (pars.fixed_b) {
        return kernel_backward(v, vp, pars);
    }
    const double beta_v = beta_from_v(v, pars, &cache);
    const double beta_vp = beta_from_v(vp, pars, &cache);
    return kernel_backward_tv_core(v, vp, beta_v, beta_vp);
}

// The non-singular part of the backwards Volterra kernel in theta from page 11 Lipton & Kaushansky(2018)
double kernel_forward(double theta, double theta_p, const RD_Params& pars) {
    const double b = pars.b_scaled;
    if (std::abs(b) < FPM_EPSILON) return 0.0;

    const double s  = 1.0 + theta;
    const double sp = 1.0 + theta_p;
    const double sum = 2.0 + theta + theta_p;
	// Handle diagonal explicitly
    if (std::abs(s - sp) < 1e-15) {
        // K(tau,tau) = - b / (sqrt (2pi) sqrt s)
		if (s < FPM_EPSILON) return 0.0;
        return -(pars.beta_prime / std::sqrt(2.0 * M_PI)) / sqrt_pos(s); // in fixed case, beta' = b_scaled
    }
    if (theta <= theta_p + FPM_EPSILON) return 0.0;

    const double ratio = (theta - theta_p) / sum;
    const double expo  = safe_exp(-b*b*ratio);

    // Eq. 10 Lipton et al (2018), -sign when moving to RHS of equation
    return -(2.0*b/std::sqrt(M_PI)) * (expo * sp) / (sum * std::sqrt(sum));
}

// Time-varying extension of above. replace b with (beta(theta)-beta(theta_p))/dtheta
// Diagonal term uses the derivative of beta at theta. In the fixed case that term reduces to b_scaled, hence the change from above
double kernel_forward_tv_core(double theta, double theta_p, double beta_theta, double beta_theta_p,const RD_Params& pars) {
  const double dv = theta - theta_p;
  const double scale = 1.0 + theta;
    // Handle diagonal explicitly
    if (std::abs(dv) < 1e-15) {
        const double beta_prime = pars.beta_prime;
        if (scale < FPM_EPSILON) return 0.0;
        return -(beta_prime / std::sqrt(2.0 * M_PI)) / sqrt_pos(scale);
    }
    const double psi = beta_theta - beta_theta_p;
    const double s = 2.0 + theta + theta_p;
    if (s <= 1e-15) return 0.0;

    const double m_eff = psi / dv;
    const double xi = safe_exp(-(psi * psi) / (dv * s));

    const double term1 = (2.0 * m_eff / std::sqrt(M_PI));
    const double term2 = ((1.0 + theta_p) * xi) / (s * std::sqrt(s));

    return -term1 * term2;
}

inline double kernel_forward_tv_cached(double theta, double theta_p, const RD_Params& pars,
                                        const BoundaryDecayCache& cache) {
    if (pars.fixed_b) {
        return kernel_forward(theta, theta_p, pars);
    }
    const double beta_v = beta_from_theta(theta, pars, &cache);
    const double beta_vp = beta_from_theta(theta_p, pars, &cache);
    return kernel_forward_tv_core(theta, theta_p, beta_v, beta_vp, pars);
}

// Dummy function for backward kernel
double g_term_backward(double v, const RD_Params& pars) {
    return 1.0;
}

// Forward kernel g term from page 11 Lipton & Kaushansky(2018)
double g_term_forward(double theta, const RD_Params& pars) {
    const double z = pars.z_scaled, b = pars.b_scaled;
    const double s = 1.0 + theta;
	const double term1 = s*s-1.0;
    if (theta <= FPM_EPSILON) return 0.0;
    const double num = s*b - z;
	const double expo = safe_exp(-(num*num)/term1);
	const double denom = sqrt_pos(M_PI * term1);
	if (denom <= FPM_EPSILON) return 0.0;
    return -expo / denom;
}

// Time-varying version of the above, replace b with beta at theta (beta includes the (1+theta) factor already)
double g_term_forward_tv(double theta, const RD_Params& pars,
                                        const BoundaryDecayCache& cache) {
    if (pars.fixed_b) {
        return g_term_forward(theta, pars);
    }
    const double z = pars.z_scaled;
    const double s = 1.0 + theta;
	const double term1 = s*s-1.0;
    if (theta <= FPM_EPSILON) return 0.0;
    const double beta = beta_from_theta(theta, pars, &cache); 
    const double num = beta - z;
	const double expo = safe_exp(-(num*num)/term1);
	const double denom = sqrt_pos(M_PI * term1);
	if (denom <= FPM_EPSILON) return 0.0;
    return -expo / denom;
}

// Abel Approximations
// Analytical solution for nu_b and nu_f for small t (Abel approximation from Eq. 24 and 26) Lipton & Kaushansky(2018)
// When b is not fixed, we use the scaled b value at time t (shouldn't make a huge difference for small t but it's more mathematically correct)
double abel_approx_nu_b(double v, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    double b_eff = pars.b_scaled;
    if (!pars.fixed_b) {
        const double beta_v = beta_from_v(v, pars, cache);
        const double scale = std::max(FPM_EPSILON, 1.0 - v);
        b_eff = beta_v / scale;
    }
    if (std::abs(b_eff) < FPM_EPSILON) {
        return 1.0;
    }
    return 2.0 * std::exp((b_eff * b_eff * v) / 2.0) * normal_cdf(b_eff * std::sqrt(v));
}

double abel_approx_nu_f(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    if (theta < FPM_EPSILON) {
        return 0.0; // nu_f(0) = g(0) = 0
    }
    double b_eff = pars.b_scaled;
    if (!pars.fixed_b) {
        const double beta_theta = beta_from_theta(theta, pars, cache);
        const double scale = std::max(FPM_EPSILON, 1.0 + theta);
        b_eff = beta_theta / scale;
    }
    const double z = pars.z_scaled;
    const double num = b_eff * theta + z - b_eff;

    const double exp_arg1 = 0.5 * b_eff * b_eff * theta + b_eff * (z - b_eff);
    const double cdf_arg = -num / std::sqrt(theta);
    const double term1 = b_eff * safe_exp(exp_arg1) * normal_cdf(cdf_arg);

    const double exp_arg2 = -((b_eff - z) * (b_eff - z)) / (2 * theta);
    const double term2 = -safe_exp(exp_arg2) / std::sqrt(2.0 * M_PI * theta);

    return term1 + term2;
}

// Smooth Integrands
// Regularized integrands: factor out 1/sqrt(v - vp)^(1/2) so we can
// integrate with Stieltjes weights on the same uniform v'-grid.
// singularity reduces because the exponential term goes to zero at v'=v
inline double G_integrand_smooth(double v, double vp, const RD_Params& pars, double nu,
                                 const BoundaryDecayCache* cache = nullptr) {
    const double dv = v - vp;
    const double s = 2.0 - v - vp;
    double Delta, one_minus_vp;
    const double s_eff = (s > 1e-14) ? s : 1e-14;
    const double dv_eff = (dv > 1e-14) ? dv : 1e-14;
    effective_positions_tv(v, vp, pars, Delta, one_minus_vp, cache);
    const double expo = safe_exp(-Delta * Delta / (dv_eff * s_eff));

    const double denom = std::sqrt(M_PI) * (s_eff * std::sqrt(s_eff)); 
    const double num = Delta * expo * one_minus_vp * nu;
    if (denom < FPM_EPSILON) return 0.0;
    return num / denom;
}

// The smooth integrand of the final PDF integral after a u-transform.
// f(u) = 2 * H_smooth(θ, θ-u²) * (v(θ-u²) - v(θ)) * (1+θ-u²) / u²
// This function computes the value of f(u) given theta and theta'.
inline double regularized_pdf_integrand_theta(
    double theta, double theta_p,
    double nu_f_theta, double nu_f_p,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr)
{
    const double dv = theta - theta_p;
    if (dv <= 0) return 0.0;

    const double beta_theta = beta_from_theta(theta, pars, cache);
    const double beta_theta_p = beta_from_theta(theta_p, pars, cache);
    const double psi = beta_theta - beta_theta_p;

    // Kernel's smooth part, H_smooth(theta, theta_p)
    const double s = 2.0 + theta + theta_p;
    const double ratio = dv / s;
    const double psi_sq_over_dv = safe_div(psi * psi, dv * dv);
    const double H_smooth = (1.0 - 2.0 * psi_sq_over_dv * ratio) * safe_exp(-psi_sq_over_dv * ratio);

    const double jacobian_part = 1.0 + theta_p;

    const double denom = dv * s * sqrt_pos(s);
    if (denom < FPM_EPSILON) return 0.0;
    // NB 2* is from the jacobian of the u-transform
    return 2.0 * (H_smooth * nu_f_p - nu_f_theta) * jacobian_part / denom;
}

// Integration Functions
// High-order integrator for the final CDF calculation. Uses product intergration with
// trapezoidal weights adapted to the sqrt singularity at the endpoint v'=v.
inline double integrate_cdf_backward_trapezoid(
    double v, const std::vector<double>& v_grid, const std::vector<double>& nu_vals,
    const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    
    double integral_sum = 0.0;
    const int N = static_cast<int>(v_grid.size()) - 1;

    for (int i = 0; i < N; ++i) { // Loop over N intervals
        const double a = v_grid[i], b = v_grid[i+1];
        const double delta_v = b - a;
        if (delta_v < FPM_EPSILON) continue;

        const double H_a = G_integrand_smooth(v, a, pars, nu_vals[i], cache);
		double H_b = 0.0; 
        const double sqrt_v_minus_a = std::sqrt(std::max(0.0, v - a));
 
        // Final interval b==v and Hb=0
        double W_a, W_b;
        if (i == N - 1) {  // Check if we are in the final interval where b ~ v
			// The smooth part H_b at v'=v is zero due to the exponential term.
			 weights_panel_12(v, a, b, W_a, W_b);
        } else {
           // Use the standard product trapezoidal formula for all other intervals.
            H_b = G_integrand_smooth(v, b, pars, nu_vals[i+1], cache);
			weights_panel_12(v, a, b, W_a, W_b);
        }

        integral_sum += W_a * H_a + W_b * H_b;
    }
    
    return integral_sum;
}


// Function to integrate the PDF's singular term using a u-transform and
// high-order panel integration (Simpson's rule). This is designed to work
// directly on the solver's native uniform theta-grid.
double integrate_pdf_forward_theta_u(
    double theta_max,
    const std::vector<double>& theta_grid,
    const std::vector<double>& nu_f_vals,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr)
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
            pars, cache
        );
    }

    // Evaluate f(u_N) = f(0) using the derivative limit
    const double mid = 0.5 * (theta_grid[N-1] + theta_max);
	const RightQuad rq = right_end_quadratic(theta_grid, nu_f_vals, mid);
    const double nu_prime_theta = rq.nu_prime;
    const double s = 1.0 + theta_max;
    const double f_right_limit = -nu_prime_theta / std::sqrt(2.0 * s);

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
        const double nuU = quad_interp(theta_grid, nu_f_vals, i, thetaU);

        // f at left endpoint
        const double fL = regularized_pdf_integrand_theta(theta_max, thetaL,
                                                      nu_f_at_max, nuL, pars, cache);

        // f at right endpoint: exact limit on the final panel, standard elsewhere
        const double fR = (i == N - 1) ? f_right_limit :  regularized_pdf_integrand_theta(theta_max, thetaR,
                                          nu_f_at_max, nuR, pars, cache);

        // f at midpoint with quadratic ν and midpoint θ′ in the smooth factors
        const double fM = regularized_pdf_integrand_theta(theta_max, thetaU,
                                                      nu_f_at_max, nuU, pars, cache);

        // Simpson on this u-panel
        integral_sum += (h / 6.0) * (fL + 4.0 * fM + fR);
    }

    return integral_sum;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SECTION 4: Compilation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// [[Rcpp::export]]
double ou_fht_cdf(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double p, int num_steps) {
    if (t < 0.0 || lambda <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }
    RD_Params pars = prepare_ou_params(t, lambda, theta, sigma, z0, b0, binf, tau, p);
    if (pars.z_scaled <= pars.b_scaled + FPM_EPSILON) { // Already hit
         return 1.0;
    }

    if (std::abs(pars.b_scaled) < FPM_EPSILON) { // Analytical case for b_scaled = 0
        double sinh_t = std::sinh(pars.t_scaled);
        if (sinh_t < FPM_EPSILON) {
            return 0.0;
        } else {
            double common_arg = -pars.z_scaled * std::exp(-0.5 * pars.t_scaled) / std::sqrt(sinh_t);
            return 2.0 * normal_cdf(common_arg);
        }
    }

    const double v = 1.0 - std::exp(-pars.t_scaled);
    const int N = (num_steps < 2) ? 2 : (num_steps % 2 ? num_steps + 1 : num_steps);

    const double h_v = v / N;
    std::vector<double> v_grid(N + 1);
    for (int j = 0; j <= N; ++j) v_grid[j] = j * h_v;

    BoundaryDecayCache beta_cache;
    const BoundaryDecayCache* cache_ptr = nullptr;
    KernelFn kernel_fn;

    if (pars.fixed_b) {
        kernel_fn = KernelFn(kernel_backward);
    } else {
        beta_cache = make_boundary_decay_cache(v, N, pars, beta_from_v_raw);
        kernel_fn = KernelFn([&beta_cache](double vv, double vvp, const RD_Params& p) {
            return kernel_backward_tv_cached(vv, vvp, p, beta_cache);
        });
        cache_ptr = &beta_cache;
    }

        const BoundaryDecayCache* abel_cache = cache_ptr;
    AbelFn abel_fn = AbelFn([abel_cache](double vv, const RD_Params& p) {
        return abel_approx_nu_b(vv, p, abel_cache);
    });

    auto nu_vals = solve_nu_block_with_abel(v, pars, N, v_grid, kernel_fn, g_term_backward, abel_fn);
    
    const double integral_sum = integrate_cdf_backward_trapezoid(v, v_grid, nu_vals, pars, cache_ptr);
    
    // Formula for G(t,z) has a coefficient of 2
    return std::max(0.0, std::min(1.0, 2*integral_sum));
}

// [[Rcpp::export]]
double ou_fht_pdf_forward(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double p, int num_steps) {
    RD_Params pars = prepare_ou_params(t, lambda, theta, sigma, z0, b0, binf, tau, p);
    if (pars.z_scaled <= pars.b_scaled + FPM_EPSILON) { // Already hit
         return 1.0;
    }

    if (t < FPM_EPSILON || lambda <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }

    if (std::abs(pars.b_scaled) < FPM_EPSILON) { // Analytical case for b_scaled = 0
        double sinh_t = std::sinh(pars.t_scaled);
        if (sinh_t < FPM_EPSILON) {
            return 0.0;
        } else {
            double common_arg = -pars.z_scaled * std::exp(-0.5 * pars.t_scaled) / std::sqrt(sinh_t);
            const double exp_term = std::exp(-0.5 * common_arg * common_arg + 0.5 * pars.t_scaled);
            const double g_scaled = (pars.z_scaled * exp_term) / (std::sqrt(2.0 * M_PI) * std::pow(sinh_t, 1.5));
            return lambda * g_scaled;
        }
    }
	const double et  = std::exp(pars.t_scaled);
	const double e2t = std::exp(2.0 * pars.t_scaled);
	const double e2tm1 = e2t - 1.0;
    const int N = (num_steps < 2) ? 2 : num_steps;
    std::vector<double> theta_grid_block;
	double theta_max = std::exp(pars.t_scaled) - 1.0;
    
    const double h_t = theta_max / N;
    std::vector<double> t_grid(N + 1);
    for (int j = 0; j <= N; ++j) t_grid[j] = j * h_t;

    BoundaryDecayCache beta_cache;
    const BoundaryDecayCache* cache_ptr = nullptr;
    KernelFn kernel_fn;

    if (pars.fixed_b) {
        kernel_fn = KernelFn(kernel_forward);
    } else {
        beta_cache = make_boundary_decay_cache(theta_max, N, pars, beta_from_theta_raw);
        kernel_fn = KernelFn([&beta_cache](double tt, double ttp, const RD_Params& p) {
            return kernel_forward_tv_cached(tt, ttp, p, beta_cache);
        });
        cache_ptr = &beta_cache;
    }

    ForcingFn g_fn;
    if(pars.fixed_b) {
        g_fn = ForcingFn(g_term_forward);
    } else {
        g_fn =  ForcingFn([cache_ptr](double theta_val, const RD_Params &p) {
                return cache_ptr ? g_term_forward_tv(theta_val, p, *cache_ptr)
                                 : g_term_forward(theta_val, p);
                });
    }

	const BoundaryDecayCache* abel_cache_theta = cache_ptr;
	AbelFn abel_fn_theta = AbelFn([abel_cache_theta](double theta_val, const RD_Params& p) {
	    return abel_approx_nu_f(theta_val, p, abel_cache_theta);
	});
    const double scale = 1+theta_max;
    std::vector<double> beta_grid(N + 1);
    for (int j = 0; j <= N; ++j) beta_grid[j] = beta_from_theta(t_grid[j], pars, cache_ptr);
    const double beta = beta_grid.back();
    pars.beta_prime = beta; // constant-barrier default
    if (!pars.fixed_b) {
      const double mid = 0.5 * (t_grid[N] + t_grid.back()); // theta_grid is N+1 in size 
      const RightQuad rq = right_end_quadratic(t_grid, beta_grid,mid); 
      pars.beta_prime = rq.nu_prime;
    }
    auto nu_f_vals = solve_nu_block_with_abel(theta_max, pars, N, theta_grid_block, kernel_fn, g_fn, abel_fn_theta);

	const double numA  = beta - pars.z_scaled;
    const double termA = - (numA * safe_exp(-((numA*numA)/e2tm1) + 2.0*pars.t_scaled)
                     / std::sqrt(M_PI * e2tm1*e2tm1*e2tm1));

	const double termB = ( pars.beta_prime*scale + e2t/(std::sqrt(M_PI*e2tm1)) ) * nu_f_vals.back();

	const double prefactor = (e2t / std::sqrt(M_PI));
    
    const double integral_sum =  integrate_pdf_forward_theta_u(theta_max, theta_grid_block, nu_f_vals, pars, cache_ptr);
    
    const double termC = prefactor * integral_sum;
    
    const double g_scaled = termA - termB + termC;
			
    return lambda * g_scaled;
}

// [[Rcpp::export]]
NumericVector ou_fht_pdf_forward_vec(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b0, double binf, double tau, double pow, int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_pdf_forward(ti, lambda, theta, sigma, z0, b0, binf, tau, pow, num_steps) : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector ou_fht_cdf_vec(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b0, double binf, double tau, double pow, int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_cdf(ti, lambda, theta, sigma, z0, b0, binf, tau, pow, num_steps) : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector simulate_ou_hit_times_std(int n,
                                        double lambda, double theta, double sigma,
                                        double z0, double b0, double binf,
                                        double tau = 1.0, double p = 1.0,
                                        double dt = 1e-3, double t_max = 10.0) {
  // 1. PARAMETER VALIDATION (same as original)
if (lambda <= 0.0 || sigma <= 0.0 || dt <= 0.0 || t_max <= 0.0) {
    stop("simulate_ou_hit_times_std: invalid parameters.");
}
// Determine if we use fixed boundary or time-varying boundary
const bool use_fixed = (tau <= 0.0) || (p <= 0.0) || (std::abs(binf - b0) <= 0.0);
// Scaling constants
const double c = std::sqrt(lambda) / sigma;
double z0_scaled = c * (z0 - theta);

// Time grid (original time for boundary eval); step count
const int steps = std::max(1, (int)std::ceil(t_max / dt));
const double ds = lambda * dt;          // scaled time step
const double s_max = lambda * t_max;    // scaled max time
// Precompute standardized boundary trajectory B_k (k=0..steps)
std::vector<double> B_scaled(steps + 1);
if (use_fixed) {
    const double B0 = c * (b0 - theta);
    std::fill(B_scaled.begin(), B_scaled.end(), B0);
} else {
    // Build t-grid and decay, then scale
    NumericVector t_grid(steps + 1);
    for (int k = 0; k <= steps; ++k) t_grid[k] = k * dt;
    NumericVector b_t = exp_decay(t_grid, std::abs(b0), std::abs(binf), tau, p);
    for (int k = 0; k <= steps; ++k) B_scaled[k] = c * (b_t[k] - theta);
}
if (z0_scaled < B_scaled[0]) {
    z0_scaled = -z0_scaled;
    for (double &Bk : B_scaled) Bk = -Bk;
}

  const double lambda_scaled = 1, theta_scaled=0, sigma_scaled=1;


  const double e = std::exp(-lambda_scaled * ds); // e = exp(-1.0 * lambda * dt)
  const double var_step = (1.0 - e * e) / (2.0 * lambda_scaled);
  const double sd_step = std::sqrt(var_step);

  NumericVector X(n);

  for (int j = 0; j < n; ++j) {
    double Z = z0_scaled; // Start with the scaled initial value
    double s = 0.0;      // Time variable is now the scaled time 's'
    bool hit = false;

    // Edge case: already at or beyond barrier at t=0
    if (Z <= B_scaled[0] + FPM_EPSILON) {
      X[j] = 0.0;
      continue;
    }

    for (int k = 0; k < steps; ++k) {
      const double Z_prev = Z;
      
      // Exact transition for standardized process (theta_std = 0.0)
      Z = Z_prev * e + sd_step * R::rnorm(0.0, 1.0);
      const double s_next = s + ds;

      
      // Moving-boundary crossing detection via sign change of D = Z - B
      const double D_prev = Z_prev - B_scaled[k];
      const double D_next = Z - B_scaled[k + 1];
      
      if (D_prev > 0.0 && D_next <= 0.0) {
        // Linear interpolation in D across the step
        const double denom = (D_prev - D_next);
        double alpha = (denom > 0.0) ? (D_prev / denom) : 0.5;
        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;

        // Convert scaled hit time back to original time
        X[j] = (s + alpha * ds) / lambda;
        hit = true;
        break;
      }

      s = s_next;
    }


    if (!hit) X[j] = NA_REAL;
  }

  return X;
}

#endif // OU_HITTING_TIME_H

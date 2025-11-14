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
#include <memory>
#include <functional>
#include <limits>
#include <cstddef>
#include <utility>
#include <Rcpp.h>
#include <quadmath.h>
using namespace Rcpp;

NumericVector ou_fht_cdf_vec_fixed_zero_branch(NumericVector t,
                                               const RD_Params &pars);





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Core Constants and Helper Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline RD_Params prepare_ou_params(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double pow,
                                 BoundaryDecayFn boundary_fn = BoundaryDecayFn(),
                                 std::vector<double> boundary_params = {}) {
  RD_Params p{};
    p.sp_var = false;
    p.b0 = b0;
    p.c = std::sqrt(lambda) / sigma;
    p.t_scaled = lambda * t;

    const double upper_scaled = p.c * (z0 - theta);
    const double b_scaled = p.c * (b0 - theta);
    const double lower_scaled = p.c * (0.0 - theta);
    const double span = std::abs(upper_scaled - lower_scaled);
    p.omega = (upper_scaled >= b_scaled) ? 1.0 : -1.0;

    p.zU_scaled = upper_scaled;
    p.zL_scaled = p.zU_scaled;
    if (span > FPM_EPSILON) {
        p.sp_var = true;
        const double lo = std::min(lower_scaled, upper_scaled);
        const double hi = std::max(lower_scaled, upper_scaled);
        p.zL_scaled = lo;
        p.zU_scaled = hi;
    }
    p.b_scaled = b_scaled;
    p.binf = binf;
    p.fixed_b = true;
    p.sigma = sigma;
    if (std::abs(p.binf - p.b0) > FPM_EPSILON) {
        p.fixed_b = false;
    }
    p.lambda = lambda;
    p.theta = theta;
    p.tau = tau;
    p.pow = pow;
    p.boundary_params = std::move(boundary_params);
    if (!boundary_fn) {
        boundary_fn = p.fixed_b ? BoundaryDecayFn(fixed_boundary_decay) : BoundaryDecayFn(default_boundary_decay);
    }
    p.boundary_fn = std::move(boundary_fn);
    return p;
}

/**
 * @brief Calculates the core numerator term for the G-integral.
 * * This computes f(y) in the point-start case, or E[f(y)] in the averaged-start case,
 * where f(y) = (y - beta) * exp(-(y - beta)^2 / variance).
 * * @param f_term        (OUT) The resulting term, f(y) or E[f(y)].
 * @param one_minus_vp  (OUT) Helper term 1 - v'.
 */
inline void calculate_f_term(double v, double vp, double variance, const RD_Params& pars,
                             double& f_term, double& one_minus_vp,
                             const BoundaryDecayCache* cache = nullptr) {

    const double alpha = std::max(0.0, 1.0 - v);
    double beta_vp;
    if (pars.fixed_b) {
        beta_vp = pars.b_scaled * std::max(0.0, 1.0 - vp);
    } else {
        beta_vp = beta_from_v(vp, pars, cache);
    }
    one_minus_vp = std::max(0.0, 1.0 - vp);

    // --- Point-start case ---
    // (If sp_var is false OR variance is zero OR the z-span is zero)
    if (!pars.sp_var || variance <= FPM_EPSILON) {
        const double z_start = pars.zU_scaled;
        const double y_start = alpha * z_start;
        const double delta = y_start - beta_vp;
        
        // Calculate f(y) = delta * exp(-delta^2 / variance)
        f_term = (variance <= FPM_EPSILON) ? 0.0 : (delta * safe_exp(-(delta * delta) / variance));
    } 
    // --- Averaged-start case ---
    else {
        double z_lo = pars.zL_scaled; // The bug fix
        double z_hi = pars.zU_scaled;
        const double span = z_hi - z_lo;

        if (span <= FPM_EPSILON) {
            // This is still a point-start, just from zL/zH
            const double y_start = alpha * z_hi; 
            const double delta = y_start - beta_vp;
            f_term = delta * safe_exp(-(delta * delta) / variance); // var is > FPM_EPSILON
        } else {
            // Call our new, clean function to get E[f(y)]
            const double y_lo = alpha * z_lo;
            const double y_hi = alpha * z_hi;
            f_term = average_of_exp_deriv(beta_vp, variance, y_lo, y_hi);
        }
    }
}

// KERNELS
// Time-varying backward volterra kernel
inline double kernel_backward_tv_core(double v, double vp, const RD_Params& pars, const BoundaryDecayCache* cache) {
    const double dv = v - vp;
    if (dv <= FPM_EPSILON) {
      return 0.0;
	}
	const double omega = pars.omega;
    const double beta_v = beta_from_v(v, pars, cache);
    const double beta_vp = beta_from_v(vp, pars, cache);
    const double psi = beta_v - beta_vp;
    if (std::abs(psi) < 1e-15) return 0.0;
    const double s = 2.0 - v - vp;
    if (s <= FPM_EPSILON) return 0.0;
    
	// effective slope
    const double b_eff = psi / dv;
	// effective time for Gstar
	// This ensures exp(-psi^2/(dv*s)) matches exp(-delta^2/(2t)) in G*
	const double t_eff = 0.5 * dv * s; 
	// change of variables transform and jacobian
	const double jacobian  = (1.0 - vp) / (s * std::sqrt(s)); 
	// -omega because we move the kernel term across to the right hand side when solving (Eq 13 Lipton & Kaushansky 2020)
	return -omega * (2.0 * b_eff * jacobian * Gstar(t_eff, psi));
}

// Time-varying forward volterra kernel
// Diagonal term uses the derivative of beta at theta. In the fixed case that term reduces to b_scaled, hence the change from above
// Uses Gstar for Gaussian structure and separates Jacobian/scaling factors for clarity.
double kernel_forward_tv_core(double theta, double theta_p, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr, const util::AkimaSpline* spline = nullptr) {
    const double dv = theta - theta_p;
    const double scale = 1.0 + theta;
    const double omega = pars.omega;
    // Handle diagonal explicitly
    if (std::abs(dv) < 1e-15) {
        // Use the local derivative beta'(theta) for the diagonal limit.
        // For fixed barriers beta'(theta) = b_scaled (constant); for
        // time-varying barriers use a cache-aware robust local finite-difference
        // that does not use the cache (for accuracy).
        double beta_prime_local = beta_prime_at_theta(theta, pars, theta, cache);
        if (scale < FPM_EPSILON) return 0.0;
        return -omega*(beta_prime_local / std::sqrt(2.0 * M_PI)) / sqrt_pos(scale);
    }
	  // Compute beta -- prefer cache (exact) otherwise use highly accurate Akima Spline
	  const double beta_theta = beta_from_theta(theta, pars, cache, spline);
    const double beta_theta_p = beta_from_theta(theta_p, pars, cache, spline);
    const double psi = beta_theta - beta_theta_p; 
    const double s = 2.0 + theta + theta_p;
    if (s <= FPM_EPSILON) return 0.0;
	
	// effective slope
    const double b_eff = psi / dv;
	// effective time for Gstar
	// This ensures exp(-psi^2/(dv*s)) matches exp(-delta^2/(2t)) in G*
	const double t_eff = 0.5 * dv * s; 
	// change of variables transform and jacobian
	const double jacobian  = ((1.0 + theta_p)) / (s * std::sqrt(s)); 
	// -omega because we move the kernel term across to the right hand side when solving (Eq 13 Lipton & Kaushansky 2020)
	return -omega * (2.0 * b_eff * jacobian * Gstar(t_eff, psi));
}

// Dummy function for backward kernel
double g_term_backward(double v, const RD_Params& pars) {
    return 1.0;
}

// Refactored to be entirely in terms of standard heat kernels for clarity and simplicity 
inline double g_term_forward(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr, const util::AkimaSpline* spline = nullptr) {
  const double scale = 1.0 + theta;
  const double tau = 0.5*(scale*scale-1.0); // this is tau, which is the solver units
	const double beta = beta_from_theta(theta, pars, cache, spline); // beta already has "scale" factored in -- in fixed_b this ends up just being scale*b_scaled
  if (theta <= FPM_EPSILON || tau <= FPM_EPSILON) return 0.0;
  double z_lo = pars.zL_scaled, z_hi = pars.zU_scaled;
	const double span = z_hi - z_lo;
  if (span <= FPM_EPSILON) {
		const double delta = pars.zU_scaled - beta;
    return -Gstar(tau,delta);
  }
	
	// If uniform z, integrate across the range
	return -Gstar_Integral(tau, beta, z_lo, z_hi)/span;
}

// NEW function for the Gompertz model (or any model with p(y) ~ e^y)
inline double g_term_forward_exp(double theta, const RD_Params& pars, const BoundaryDecayCache* cache, const util::AkimaSpline* spline = nullptr) {
    const double scale = 1.0 + theta;
    const double tau = 0.5 * (scale * scale - 1.0);
    const double beta = beta_from_theta(theta, pars, cache, spline);
    if (theta <= FPM_EPSILON || tau <= FPM_EPSILON) return 0.0;

    // Handle point-start case (same as before)
    if (!pars.sp_var) {
        const double delta = pars.zU_scaled - beta;
        return -Gstar(tau, delta);
    }

    double y_lo = pars.zL_scaled;
    double y_hi = pars.zU_scaled;
    if (y_hi < y_lo) std::swap(y_lo, y_hi);

    const double norm_const = std::exp(y_hi) - std::exp(y_lo);
    if (norm_const <= FPM_EPSILON) return 0.0;

    // Completing the square introduces these terms
    const double exp_factor = std::exp(beta + tau / 2.0);
    const double new_mean = beta + tau;

    const double integral_G = Gstar_Integral(tau, new_mean, y_lo, y_hi);

    // The final averaged value, including the minus sign from the original formula
    return - (1.0 / norm_const) * exp_factor * integral_G;
}

// Abel Approximations
// Analytical solution for nu_b and nu_f for small t (Abel approximation from Eq. 24 and 26) Lipton & Kaushansky(2018)
// When b is not fixed, we use the scaled b value at time t (shouldn't make a huge difference for small t but it's more mathematically correct)
double abel_approx_nu_b(double v, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    double b_eff = pars.b_scaled;
    const double omega = pars.omega;
    const double b_signed = omega * b_eff;
    const double z_signed = omega * pars.zU_scaled;
    if (!pars.fixed_b) {
        const double beta_v = omega*beta_from_v(v, pars, cache);
        const double scale = std::max(FPM_EPSILON, 1.0 - v);
        b_eff = beta_v / scale;
    }
    if (std::abs(b_eff) < FPM_EPSILON) {
        return 1.0;
    }
    return 2.0 * std::exp((b_signed * b_signed * v) / 2.0) * gaussian_cdf(b_signed * std::sqrt(v));
}

// Updated to do uniform start better (needs checking)
double abel_approx_nu_f(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr, const util::AkimaSpline* spline = nullptr) {
    
	if (theta < FPM_EPSILON) {
        return 0.0; // nu_f(0) = g(0) = 0
  }

  const double beta = beta_from_theta(theta, pars, cache, spline);
  const double scale = 1.0 + theta;
  const double b = beta / scale;
  const double omega = pars.omega;
	
  const double b_signed = omega * b;
	const double z_lo = pars.zL_scaled;
  const double z_hi = pars.zU_scaled;
  const double span = z_hi - z_lo;
	
	// Generalise for both point and uniform versions
	auto term1_function = [=] (double z) {
		const double z_signed = omega * z;
		const double num = b_signed * theta + z_signed - b_signed;
		const double exp_arg1 = 0.5 * b_signed * b_signed * theta + b_signed * (z_signed - b_signed);
		const double cdf_arg = -num / std::sqrt(theta);
		return b_signed * safe_exp(exp_arg1) * gaussian_cdf(cdf_arg);
	};
	
	if (span<FPM_EPSILON) {
    const double z_signed = omega * z_hi;
		const double term1 = term1_function(z_hi);
		const double delta = b_signed - omega*z_signed;
		const double term2 = -Gstar(theta,delta); 
		return term1 + term2;
	}
	
	// Time-varying extension
	// Term 2 uses exact integral
	const double term2 = -Gstar_Integral(theta,b,z_lo,z_hi)/span;
	
	// --- Term 1: Approximate integral using 3-point Gauss–Legendre quadrature
    // GL weights and nodes for 3-point quadrature on [-1,1] using static const Rcpp::Function gauss_quad = statmodNS["gauss.quad"];
	static Rcpp::List gl_3 = gauss_quad(3, "legendre");
	const Rcpp::NumericVector& x = gl_3["nodes"];
	const Rcpp::NumericVector& w = gl_3["weights"];

  double term1_integral = 0.0;
  for (int i = 0; i < 3; ++i) {
    double z = ((z_hi - z_lo) * x[i] + (z_hi + z_lo)) / 2.0; // map [-1,1] → [z_lo,z_hi]
    term1_integral += w[i] * term1_function(z);
  }
  term1_integral *= (z_hi - z_lo) / 2.0; // scale by interval length
  const double term1 = term1_integral / span;

	return term1 + term2;  
}

double abel_approx_nu_f_exp(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr, const util::AkimaSpline* spline = nullptr) {

    if (theta < FPM_EPSILON) {
        return 0.0;
    }

    const double scale = 1.0 + theta;
    const double tau = 0.5 * (scale * scale - 1.0);
    const double beta = beta_from_theta(theta, pars, cache, spline);

    const double omega = pars.omega;

    double y_lo = pars.zL_scaled; // These are already log-transformed
    double y_hi = pars.zU_scaled;
    if (y_hi < y_lo) std::swap(y_lo, y_hi);

    // --- Same term1_function as before (now takes y as input) ---
    auto term1_function = [&](double y) {
        const double y_signed = omega * y;
        const double b_signed = omega * (beta / scale);
        const double num = b_signed * theta + y_signed - b_signed;
        const double exp_arg1 = 0.5 * b_signed * b_signed * theta + b_signed * (y_signed - b_signed);
        const double cdf_arg = -num / std::sqrt(theta);
        return b_signed * safe_exp(exp_arg1) * gaussian_cdf(cdf_arg);
    };

    // --- Handle point-start case (same logic as before) ---
    if (!pars.sp_var) {
        const double term1 = term1_function(y_hi);
        const double delta = (beta / scale) - y_hi;
        const double term2 = -Gstar(tau, delta); // Note: tau not theta here for Gstar
        return term1 + term2;
    }

    // --- Distributed start-point case ---

    // 1. Calculate Term 2 using our new analytic solution
    const double term2 = g_term_forward_exp(theta, pars, cache, spline);

    // 2. Calculate Term 1 using quadrature with the new PDF
    const double norm_const = std::exp(y_hi) - std::exp(y_lo);
    if (norm_const <= FPM_EPSILON) return term2; // term1 will be zero

    auto start_point_pdf = [&](double y) {
        return std::exp(y) / norm_const;
    };
    
    // Use the same 3-point Gauss-Legendre quadrature
    static Rcpp::List gl_3 = Rcpp::Function("statmod::gauss.quad")(3, "legendre");
	  const Rcpp::NumericVector& x = gl_3["nodes"];
	  const Rcpp::NumericVector& w = gl_3["weights"];

    double term1_integral = 0.0;
    for (int i = 0; i < 3; ++i) {
        double y = ((y_hi - y_lo) * x[i] + (y_hi + y_lo)) / 2.0; // map [-1,1] -> [y_lo, y_hi]
        // The new integrand is the original function multiplied by the new PDF
        term1_integral += w[i] * term1_function(y) * start_point_pdf(y);
    }
    term1_integral *= (y_hi - y_lo) / 2.0; // scale by interval length

    const double term1 = term1_integral; // This is the final averaged term1

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

    const double s_eff = (s > 1e-14) ? s : 1e-14;
    const double dv_eff = (dv > 1e-14) ? dv : 1e-14;
    const double variance = dv_eff * s_eff;

    double f_term; // This is f(y) or E[f(y)]
    double one_minus_vp;

    calculate_f_term(v, vp, variance, pars, f_term, one_minus_vp, cache);
    
    const double denom = std::sqrt(M_PI) * (s_eff * std::sqrt(s_eff)); 
    const double num = f_term * one_minus_vp * nu;           

    if (denom < FPM_EPSILON)
        return 0.0;
    return num / denom;
}

// The smooth integrand of the final PDF integral after a u-transform.
// This function computes the value of f(u) given theta and theta'.
inline double regularized_pdf_integrand(double theta, double theta_p,
                                        double nu_f_theta, double nu_f_p,
                                        const RD_Params &pars,
                                        const BoundaryDecayCache* cache = nullptr,
                                        const util::AkimaSpline* beta_spline = nullptr)
  {

    const double dv = theta - theta_p;
    const double scale = 1.0 + theta_p;
    const double s = 2.0 + theta + theta_p;
    if (dv <= 0) return 0.0;
    if (s <= 1e-15) return 0.0;

    const double beta_theta =   beta_from_theta(theta, pars, cache, beta_spline);
    const double beta_theta_p = beta_from_theta(theta_p, pars, cache, beta_spline);
    const double psi = beta_theta - beta_theta_p; // this is equivalent to b*dv in the fixed case
    const double t_eff = 0.5 * dv * s;

    // Get the normalized heat kernel value
    const double Gstar_val = Gstar(t_eff, psi); // accounts for one power of dv in denominator, as well as b_eff
    // This is the (1 - delta^2/t_eff) part characteristic of the 2nd derivative
    const double H2_part = (1.0 - (psi * psi) / t_eff);
    // Reconstruct the unnormalized gaussian 'xi' from Gstar
    // Gstar_val = xi / sqrt(2*pi*t_eff)  =>  xi = Gstar_val * sqrt(2*pi*t_eff)
    const double xi = Gstar_val * std::sqrt(2.0 * M_PI * t_eff);
    const double H_smooth = H2_part * xi;
    // Get the coordinate change jacobian
    // jacobian = (1.0 + theta_p) / (s * std::sqrt(s))
    const double jacobian = (1.0 + theta_p) / (s * std::sqrt(s));

    // NB 2* is from the jacobian of the u-transform
    return 2.0 * (H_smooth * nu_f_p - nu_f_theta) * (jacobian / dv);
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

/**
 * Computes the special PDF integral up to an ARBITRARY endpoint theta_k.
 * This function is self-contained and relies *only* on the spline.
 *
 * @param theta_k The arbitrary endpoint (your t_obs, transformed).
 * @param nu_spline The pre-built Akima spline for the kernel.
 * @param M The number of panels to use for this integration.
 * This can be small (e.g., 50-100) and fixed.
 * @return The value of the integral from 0 to theta_k.
 */
double integrate_pdf_forward(
    double theta_k,
    const util::AkimaSpline& nu_spline,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr,
    int M = 100,
    double nu_prime_at_k_override = std::numeric_limits<double>::quiet_NaN(),
    const util::AkimaSpline* beta_spline = nullptr)
{
    if (theta_k <= 0.0 || M <= 0) return 0.0;

    // Construct a uniform grid size M up to theta
    std::vector<double> temp_theta_grid(M + 1);
    std::vector<double> temp_u_grid(M + 1);
    const double h_theta = theta_k / M;
    OU_DEBUG_LOG(4, "integrate_pdf_forward: theta_k=" << theta_k << ", M=" << M);

    for (int i = 0; i <= M; ++i) {
        temp_theta_grid[i] = i * h_theta;
        temp_u_grid[i] = sqrt_pos(theta_k - temp_theta_grid[i]);
    }

    std::vector<double> theta_mid(M);
    for (int i = 0; i < M; ++i) {
        const double uL = temp_u_grid[i];
        const double uR = temp_u_grid[i + 1];
        const double uM = 0.5 * (uL + uR);
        theta_mid[i] = theta_k - uM * uM;
    }

    std::vector<double> nu_theta_grid(M + 1);
    for (int i = 0; i <= M; ++i) {
        nu_theta_grid[i] = nu_spline.interpolate(temp_theta_grid[i]);
    }
    std::vector<double> nu_mid(M);
    for (int i = 0; i < M; ++i) {
        nu_mid[i] = nu_spline.interpolate(theta_mid[i]);
    }

    const double nu_at_k = nu_theta_grid.back();
    const double nu_prime_at_k = std::isfinite(nu_prime_at_k_override)
                                   ? nu_prime_at_k_override
                                   : nu_spline.derivative(theta_k);
    const double s = 1.0 + theta_k;
    const double f_right_limit = -nu_prime_at_k / std::sqrt(2.0 * s);
    OU_DEBUG_LOG(5, "integrate_pdf_forward: nu_at_k=" << nu_at_k
                    << ", nu_prime_at_k=" << nu_prime_at_k
                    << ", s=" << s
                    << ", f_right_limit=" << f_right_limit);

    double integral_sum = 0.0;

    const double beta_theta_k = beta_from_theta(theta_k, pars, cache, beta_spline);

    for (int i = M - 1; i >= 0; --i) {
        const double uL = temp_u_grid[i];
        const double uR = temp_u_grid[i + 1];
        const double h_u = uL - uR;
        if (h_u == 0.0) continue;

        const double thetaL = temp_theta_grid[i];
        const double thetaR = temp_theta_grid[i + 1];
        const double theta_mid_val = theta_mid[i];

        const double nuL = nu_theta_grid[i];
        const double nuR = nu_theta_grid[i + 1];
        const double nuM = nu_mid[i];

        // f at left endpoint
        const double fL = regularized_pdf_integrand(theta_k, thetaL,
                                                      nu_at_k, nuL, pars, cache, beta_spline);

        // f at right endpoint: exact limit on the final panel, standard elsewhere
        const double fR = (i == M - 1) ? f_right_limit :  regularized_pdf_integrand(theta_k, thetaR,
                                          nu_at_k, nuR, pars, cache, beta_spline);

        // f at midpoint
        const double fM = regularized_pdf_integrand(theta_k, theta_mid_val,
                                                      nu_at_k, nuM, pars, cache, beta_spline);

        const double incr = (h_u / 6.0) * (fL + 4.0 * fM + fR);
        integral_sum += incr;
        OU_DEBUG_LOG(6, "integrate_pdf_forward panel i=" << i
                        << ": h_u=" << h_u
                        << ", fL=" << fL << ", fM=" << fM << ", fR=" << fR
                        << ", incr=" << incr);
    }
    OU_DEBUG_LOG(4, "integrate_pdf_forward: integral_sum=" << integral_sum);
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
  const double v_max = 1.0 - std::exp(-pars.t_scaled);
  pars.v_max = v_max;
  const double omega = pars.omega;
  const double b_at_t = (beta_from_v(v_max, pars) / (1 - v_max)) / pars.c + theta; // unscaled b at time t
    if (std::abs(b_at_t) < FPM_EPSILON) { // Already hit
         return 1.0;
    }

    /*
    if (pars.fixed_b && std::abs(pars.b_scaled) < FPM_EPSILON) { // Analytical
    case for b_scaled = 0 NumericVector t_vec = NumericVector::create(t);
        NumericVector cdf_vec = ou_fht_cdf_vec_fixed_zero_branch(t_vec, pars);
        const double cdf_val = cdf_vec[0];
        if (!R_finite(cdf_val)) return NA_REAL;
        return std::max(0.0, std::min(1.0, cdf_val));
    }
    */
    
    const int N = (num_steps < 2) ? 2 : num_steps;

    const double h_v = v_max / N;
    std::vector<double> v_grid(N + 1);
    for (int j = 0; j <= N; ++j) v_grid[j] = j * h_v;

    BoundaryDecayCache beta_cache;
    std::unique_ptr<util::AkimaSpline> beta_spline;
    const BoundaryDecayCache* cache_ptr = nullptr;
    KernelFn kernel_fn;

    if (pars.fixed_b) {
         kernel_fn = KernelFn([&beta_cache](double vv, double vvp, const RD_Params& p) {
            return kernel_backward_tv_core(vv, vvp, p, &beta_cache);
        });
    } else {
        std::vector<double> beta_vals(v_grid.size());
        beta_cache = make_boundary_decay_cache(v_grid, pars);
        cache_ptr = &beta_cache;
        for (size_t j = 0; j < v_grid.size(); ++j) {
            beta_vals[j] = beta_from_theta(v_grid[j], pars, cache_ptr);
        }
        beta_spline = std::make_unique<util::AkimaSpline>(v_grid, beta_vals);
        kernel_fn = KernelFn([&beta_cache](double vv, double vvp, const RD_Params& p) {
            return kernel_backward_tv_core(vv, vvp, p, &beta_cache);
        });
    }

        const BoundaryDecayCache* abel_cache = cache_ptr;
    AbelFn abel_fn = AbelFn([abel_cache](double vv, const RD_Params& p) {
        return abel_approx_nu_b(vv, p, abel_cache);
    });

    auto nu_vals = solve_nu_block_with_abel(v_max, pars, N, v_grid, kernel_fn, g_term_backward, abel_fn);

    const double integral_sum = integrate_cdf_backward_trapezoid(v_max, v_grid, nu_vals, pars, cache_ptr);
    const double cdf = omega * 2.0 * integral_sum;
    return std::max(0.0, std::min(cdf,1.0));
}

NumericVector ou_fht_pdf_vec_fixed_zero_branch(NumericVector t,
                                               const RD_Params &pars) {
  const int n = t.size();
  NumericVector out(n);
  const double lambda = pars.lambda;
  const double omega = pars.omega;
  const bool has_uniform_start =
        pars.sp_var && (std::abs(pars.zU_scaled - pars.zL_scaled) > FPM_EPSILON);
    for (int i = 0; i < n; ++i) {
      const double ti = t[i];
      if (!R_finite(ti) || ti <= 0.0) {
        out[i] = R_finite(ti) ? 0.0 : NA_REAL;
        continue;
      }

      const double t_scaled = lambda * ti;
      const double sinh_t = std::sinh(t_scaled);
      if (sinh_t < FPM_EPSILON) {
        out[i] = 0.0;
        continue;
      }

      const double prefactor =
          lambda * std::exp(0.5 * t_scaled) /
          (std::sqrt(2.0 * M_PI) * std::pow(sinh_t, 1.5));
      const double decay_coeff = 0.5 * std::exp(-t_scaled) / sinh_t;

      if (!has_uniform_start) {
        const double z_signed = omega * pars.zU_scaled;
        out[i] = prefactor * z_signed * std::exp(-decay_coeff * z_signed * z_signed);
        out[i] = std::max(0.0, out[i]);
        continue;
      }

      double z_lo = pars.zL_scaled;
      double z_hi = pars.zU_scaled;
      if (z_hi < z_lo) std::swap(z_lo, z_hi);
      const double span = z_hi - z_lo;

      if (span <= FPM_EPSILON) {
        const double z = 0.5 * (z_lo + z_hi);
        const double z_signed = omega * z;
        out[i] = prefactor * z_signed * std::exp(-decay_coeff * z_signed * z_signed);
        out[i] = std::max(0.0,out[i]);
        continue;
      }

      double averaged_z;
      if (decay_coeff <= FPM_EPSILON) {
        averaged_z = 0.5 * (z_lo + z_hi);
      } else {
        const double exp_lo = std::exp(-decay_coeff * z_lo * z_lo);
        const double exp_hi = std::exp(-decay_coeff * z_hi * z_hi);
        averaged_z = (0.5 / (decay_coeff * span)) * (exp_lo - exp_hi);
      }
      const double averaged_z_signed = omega * averaged_z;
      out[i] = prefactor * averaged_z_signed;
      out[i] = std::max(0.0,out[i]);
    }
    return out;
}

double ou_fht_pdf_forward(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double p, double steps_fineness, double min_steps) {
    RD_Params pars = prepare_ou_params(t, lambda, theta, sigma, z0, b0, binf, tau, p);
    const double omega = pars.omega;
    if (t < FPM_EPSILON || lambda <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }
    /*
    if (std::abs(pars.b_scaled) < FPM_EPSILON && pars.fixed_b) { // Analytical case for b_scaled = 0
        NumericVector t_vec(1);
        t_vec[0] = t;
        NumericVector pdf_vec = ou_fht_pdf_vec_fixed_zero_branch(t_vec, pars);
        return pdf_vec[0];
    }
    */
	  const double et  = std::exp(pars.t_scaled);
	  const double e2t = std::exp(2.0 * pars.t_scaled);
    const double e2tm1 = e2t - 1.0;
    double theta_max = et - 1.0;
    const double tau_max = 0.5 * (e2tm1);
    int N = calculate_num_steps(theta_max, steps_fineness, min_steps);
    std::vector<double> theta_grid_block;
    
    pars.theta_max = theta_max;
    OU_DEBUG_LOG(1, "ou_fht_pdf_forward: t=" << t << ", theta_max=" << theta_max
                    << ", tau_max=" << tau_max << ", N=" << N
                    << ", fixed_b=" << (pars.fixed_b?1:0));

    if (omega*(pars.zU_scaled - pars.b_scaled) < FPM_EPSILON) { 
         return 0.0;
    }
    double h_t = theta_max / N;
    std::vector<double> theta_uniform(N + 1);
    for (int j = 0; j <= N; ++j) theta_uniform[j] = j * h_t;

    // Build a beta spline on the solver uniform theta grid
    // Construct a cache that will hold beta and beta_prime values to save re-calculating
    BoundaryDecayCache beta_cache;
    const BoundaryDecayCache* cache_ptr = nullptr;
    KernelFn kernel_fn;
    ForcingFn g_fn;
    std::unique_ptr<util::AkimaSpline> beta_spline;
    if (pars.fixed_b) {
         kernel_fn = KernelFn([](double tt, double ttp, const RD_Params& p) {
          return kernel_forward_tv_core(tt, ttp, p);
        });
        g_fn = ForcingFn([](double theta_val, const RD_Params& p) {
          return g_term_forward(theta_val, p);
        });
    } else {
      std::vector<double> beta_vals(theta_uniform.size());
      beta_cache = make_boundary_decay_cache(theta_uniform, pars);
      cache_ptr = &beta_cache;
        for (size_t j = 0; j < theta_uniform.size(); ++j) {
            beta_vals[j] = beta_from_theta(theta_uniform[j], pars, cache_ptr);
        }
        beta_spline = std::make_unique<util::AkimaSpline>(theta_uniform, beta_vals);
        const double inv_h = (theta_max > 0.0) ? (static_cast<double>(N) / theta_max) : 0.0;
        kernel_fn = KernelFn([cache_ptr,&beta_spline](double tt, double ttp, const RD_Params& p) {
          return kernel_forward_tv_core(tt, ttp, p, cache_ptr, beta_spline.get());
        });
         g_fn = ForcingFn([cache_ptr, &beta_spline](double theta_val, const RD_Params &p) {
              return g_term_forward(theta_val, p, cache_ptr,beta_spline.get());
            });
      }
    

    AbelFn abel_fn_theta = AbelFn([cache_ptr,&beta_spline](double theta_val, const RD_Params& p) {
        return abel_approx_nu_f(theta_val, p, p.fixed_b ? nullptr : cache_ptr,p.fixed_b ? nullptr : beta_spline.get());
    });
    const double scale = 1+theta_max;
    const double beta_t = pars.fixed_b ? (1.0 + theta_max) * pars.b_scaled
                                       : beta_spline->interpolate(theta_max);
    const double beta_prime = beta_prime_at_theta(theta_max, pars, theta_max);

    auto nu_f_vals = solve_nu_block_with_abel(theta_max, pars, N, theta_grid_block, kernel_fn, g_fn, abel_fn_theta);
    const double nu_t = nu_f_vals.back();

    double nu_prime_at_max = 0.0;
    if (theta_grid_block.size() >= 3) {
        const double mid_nu = 0.5 * (theta_grid_block[N - 1] + theta_max);
        RightQuad rq_nu = right_end_quadratic(theta_grid_block, nu_f_vals, mid_nu);
        nu_prime_at_max = rq_nu.nu_prime;
    }

    util::AkimaSpline nu_spline(theta_grid_block, nu_f_vals);
    // Equations here are re-arranged to match general form in Lipton &
    // Kaushansky (2020) Equation 10 They are mathematically equivalent to the
    // original form in section 4.4, but easier to relate to the general case
    // Term 1: The integral from Eq. (10)
    int M_eval = std::max(80, std::min(240, 2 * N));
    const double integral_sum = integrate_pdf_forward(theta_max, nu_spline, pars, cache_ptr, M_eval, nu_prime_at_max, pars.fixed_b ? nullptr : beta_spline.get());
    const double termA = (1/std::sqrt(M_PI))*integral_sum;
    const double termB = -nu_t / std::sqrt(2.0 * M_PI * tau_max);
    const double beta_prime_term = omega*(-(beta_prime/scale) * nu_t);
    const double image_term = omega* (0.5* averaged_image_term(tau_max, beta_t, pars));
    const double g_scaled = termA + termB + beta_prime_term + image_term;
    const double pdf = lambda * g_scaled * e2t;
    OU_DEBUG_LOG(1, "ou_fht_pdf_forward terms: nu_t=" << nu_t
                    << ", nu_prime_at_max=" << nu_prime_at_max
                    << ", integral_sum=" << integral_sum
                    << ", termA=" << termA << ", termB=" << termB
                    << ", beta_prime_term=" << beta_prime_term
                    << ", image_term=" << image_term
                    << ", scale=" << scale << ", tau=" << tau_max
                    << ", beta_t=" << beta_t << ", b_scaled=" << pars.b_scaled
                    << ", zL_scaled=" << pars.zL_scaled << ", zU_scaled=" << pars.zU_scaled
                    << ", pdf=" << pdf);
    return std::max(0.0, pdf);
}

// [[Rcpp::export]]
NumericVector ou_fht_pdf_vec(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b0, double binf, double tau, double pow, double steps_fineness=0.005, double min_steps=300) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_pdf_forward(ti, lambda, theta, sigma, z0, b0, binf, tau, pow, steps_fineness, min_steps) : NA_REAL;
  }
  return out;
}

// Forward declaration for chunked variant used by the wrapper below
NumericVector ou_fht_pdf_vec_grid_chunked(NumericVector t,
                                          double lambda, double theta, double sigma,
                                          double z0, double b0, double binf,
                                          double tau, double pow, double steps_fineness, double min_steps,double min_t_fineness,
                                          double chunk_ratio, double chunk_base_panels,
                                          double chunk_max);

// Forward declaration for chunked CDF variant (used by the unchunked wrapper)
NumericVector ou_fht_cdf_vec_grid_chunked(NumericVector t,
                                          double lambda, double theta, double sigma,
                                          double z0, double b0, double binf,
                                          double tau, double pow, double steps_fineness, double min_steps,double min_t_fineness,
                                          double chunk_ratio, double chunk_base_panels,
                                          double chunk_max,
                                          double rt_resolution);

// [[Rcpp::export]]
NumericVector ou_fht_cdf_vec_grid(NumericVector t,
                                  double lambda, double theta, double sigma,
                                  double z0, double b0, double binf,
                                  double tau, double pow, double steps_fineness, double min_steps, double min_t_fineness,
                                  double rt_resolution = 0.02)
{
  const int n = t.size();
  // Determine horizon and earliest positive request time
  double t_max = 0.0;
  double t_min_pos = std::numeric_limits<double>::infinity();
  bool has_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i])) {
      if (t[i] > t_max) { t_max = t[i]; has_finite = true; }
      if (t[i] > 0.0 && t[i] < t_min_pos) t_min_pos = t[i];
    }
  }
  if (!has_finite || t_max <= 0.0 || lambda <= 0.0 || sigma <= 0.0) {
    NumericVector out(n, NA_REAL);
    for (int i = 0; i < n; ++i) out[i] = (R_finite(t[i]) && t[i] >= 0.0) ? 0.0 : NA_REAL;
    return out;
  }

  // Compute base panels for zero-compression handoff, mirroring the PDF wrapper
  int N = calculate_num_steps(std::exp(lambda*t_max)-1.0, steps_fineness, (int)min_steps);
  int P = 0;
  if (R_finite(t_min_pos) && t_min_pos > 0.0) {
    P = calculate_num_steps(std::exp(lambda*t_min_pos)-1.0, steps_fineness, (int)min_steps);
    if (P % 2 != 0) ++P;
    if (P > N - 2) P = std::max(2, N - 2);
    if ((N - P) % 2 != 0) { if (P > 2) --P; else ++N; }
  }
  const double ratio = 1.0;
  const double base_panels = static_cast<double>(std::max(2, N - P));
  const double max_chunks = 2.0;

  return ou_fht_cdf_vec_grid_chunked(t, lambda, theta, sigma,
                                     z0, b0, binf,
                                     tau, pow,
                                     steps_fineness, min_steps, min_t_fineness,
                                     ratio, base_panels, max_chunks,
                                     rt_resolution);
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


NumericVector ou_fht_cdf_vec_fixed_zero_branch(NumericVector t,
                                               const RD_Params &pars) {
  const int n = t.size();
  NumericVector out(n);
  const double lambda = pars.lambda;
  const double omega = pars.omega;
  const bool has_uniform_start =
      pars.sp_var && (std::abs(pars.zU_scaled - pars.zL_scaled) > FPM_EPSILON);

  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti)) {
      out[i] = NA_REAL;
      continue;
    }
    if (ti <= 0.0) {
      out[i] = 0.0;
      continue;
    }

    const double t_scaled = lambda * ti;
    const double sinh_t = std::sinh(t_scaled);
    if (sinh_t < FPM_EPSILON) {
      out[i] = 0.0;
      continue;
    }

    const double sqrt_sinh = std::sqrt(sinh_t);
    const double atten = std::exp(-0.5 * t_scaled) / sqrt_sinh;

    if (!has_uniform_start) {
      const double z_signed = omega * pars.zU_scaled;
      const double arg = -z_signed * atten;
      const double cdf = 2.0 * gaussian_cdf(arg);
      out[i] = std::min(1.0, std::max(0.0, cdf));
      continue;
    }

    double z_lo = pars.zL_scaled;
    double z_hi = pars.zU_scaled;
    if (z_hi < z_lo) std::swap(z_lo, z_hi);
    const double span = z_hi - z_lo;

    if (span <= FPM_EPSILON) {
      const double z_mid = 0.5 * (z_lo + z_hi);
      const double z_signed = omega * z_mid;
      const double arg = -z_signed * atten;
      const double cdf = 2.0 * gaussian_cdf(arg);
      out[i] = std::min(1.0, std::max(0.0, cdf));
      continue;
    }

    const double a = omega * atten;
    if (std::abs(a) <= FPM_EPSILON) {
      out[i] = 1.0;
      continue;
    }

    auto antideriv = [a](double z_val) {
      const double y = -a * z_val;
      const double Phi = gaussian_cdf(y);
      const double phi = gaussian_pdf(y);
      return z_val * Phi - (1.0 / a) * phi;
    };

    double cdf = (2.0 / span) * (antideriv(z_hi) - antideriv(z_lo));
    if (!R_finite(cdf)) {
      cdf = (cdf > 0.0) ? 1.0 : 0.0;
    }
    out[i] = std::min(1.0, std::max(0.0, cdf));
  }
  return out;
}



// [[Rcpp::export]]
NumericVector ou_fht_pdf_vec_closed_form(NumericVector t,
                                         double lambda, double theta, double sigma,
                                         double z0, double b0, double binf,
                                         double tau = 1.0, double pow = 1.0) {
  const int n = t.size();
  NumericVector out(n, NA_REAL);
  // Find max finite nonnegative time to set scaling
  double t_max = 0.0;
  bool any_pos = false;
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (R_finite(ti)) {
      if (ti > t_max) t_max = ti;
      if (ti >= 0.0) any_pos = true;
    }
  }
  if (!any_pos || !(lambda > 0.0) || !(sigma > 0.0)) {
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) out[i] = 0.0;
      else out[i] = NA_REAL;
    }
    return out;
  }

  RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, pow);
  // Require fixed, zero scaled barrier to apply this closed form
  if (!pars.fixed_b || std::abs(pars.b_scaled) > FPM_EPSILON) {
    stop("ou_fht_pdf_vec_closed_form: requires fixed barrier at theta (scaled zero barrier).");
  }
  return ou_fht_pdf_vec_fixed_zero_branch(t, pars);
}

// [[Rcpp::export]]
NumericVector ou_fht_cdf_vec_closed_form(NumericVector t,
                                         double lambda, double theta, double sigma,
                                         double z0, double b0, double binf,
                                         double tau = 1.0, double pow = 1.0) {
  const int n = t.size();
  NumericVector out(n, NA_REAL);
  double t_max = 0.0;
  bool any_pos = false;
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (R_finite(ti)) {
      if (ti > t_max) t_max = ti;
      if (ti >= 0.0) any_pos = true;
    }
  }
  if (!any_pos || !(lambda > 0.0) || !(sigma > 0.0)) {
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) out[i] = 0.0;
    }
    return out;
  }

  RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, pow);
  if (!pars.fixed_b || std::abs(pars.b_scaled) > FPM_EPSILON) {
    stop("ou_fht_cdf_vec_closed_form: requires fixed barrier at theta (scaled zero barrier).");
  }

  NumericVector res = ou_fht_cdf_vec_fixed_zero_branch(t, pars);
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i]) && t[i] < 0.0) res[i] = 0.0;
  }
  return res;
}

// [[Rcpp::export]]
NumericVector ou_fht_pdf_vec_grid(NumericVector t,
                                  double lambda, double theta, double sigma,
                                  double z0, double b0, double binf,
                                  double tau, double pow, double steps_fineness, double min_steps, double min_t_fineness)
{
  const int n = t.size();
  // Find horizon and earliest positive time
  double t_max = 0.0;
  double t_min_pos = std::numeric_limits<double>::infinity();
  bool has_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i])) {
      if (t[i] > t_max) { t_max = t[i]; has_finite = true; }
      if (t[i] > 0.0 && t[i] < t_min_pos) t_min_pos = t[i];
    }
  }
  if (!has_finite || t_max <= 0.0 || lambda <= 0.0 || sigma <= 0.0) {
    NumericVector out(n, NA_REAL);
    for (int i = 0; i < n; ++i) out[i] = (R_finite(t[i]) && t[i] >= 0.0) ? 0.0 : NA_REAL;
    return out;
  }
  double theta_max = std::exp(t_max*lambda) - 1.0;
  // Compute total panels and prelude panels using same step rule
  int N = calculate_num_steps(theta_max, steps_fineness, min_steps);
  int P = 0;
  if (R_finite(t_min_pos) && t_min_pos > 0.0) {
    P = calculate_num_steps(std::exp(t_min_pos*lambda)-1, min_t_fineness, min_steps);
    if (P % 2 != 0) ++P;
    if (P > N - 2) P = std::max(2, N - 2);
    if ((N - P) % 2 != 0) { if (P > 2) --P; else ++N; }
  }
  const double ratio = 1.0;
  const double base_panels = static_cast<double>(std::max(2, N - P));
  const double max_chunks = 2.0;
  return ou_fht_pdf_vec_grid_chunked(t, lambda, theta, sigma, z0, b0, binf,
                                     tau, pow, steps_fineness, min_steps,min_t_fineness,
                                     ratio, base_panels, max_chunks);
}


// Variant of the grid-based solver using the chunked Volterra kernel
// acceleration. Keeps the legacy function untouched for comparative testing.
// [[Rcpp::export]]
NumericVector ou_fht_pdf_vec_grid_chunked(NumericVector t, double lambda,
                                          double theta, double sigma, double z0,
                                          double b0, double binf, double tau,
                                          double pow, double steps_fineness,
                                          double min_steps,
                                          double min_t_fineness,
                                          double chunk_ratio, double chunk_base_panels,
                                          double chunk_max)
{
  const int n = t.size();
  NumericVector out(n, NA_REAL);

  double t_max = 0.0;
  bool has_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i]) && t[i] > t_max) {
      t_max = t[i];
      has_finite = true;
    }
  }
  if (!has_finite || t_max <= 0.0 || lambda <= 0.0 || sigma <= 0.0) {
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) {
        out[i] = 0.0;
      } else {
        out[i] = NA_REAL;
      }
    }
    return out;
  }
  // Identify earliest positive time to improve early-time resolution
  double t_min_pos = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i]) && t[i] > 0.0) {
      if (t[i] < t_min_pos) t_min_pos = t[i];
    }
  }
  RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, pow);
  
  // Construct a cache that will hold beta and beta_prime values to save re-calculating
  BoundaryDecayCache beta_cache;
  const BoundaryDecayCache* cache_ptr = nullptr;
  
  // Early exit when analytic solution is valid
  /*
  if (std::abs(pars.b_scaled) < FPM_EPSILON && pars.fixed_b) {
    NumericVector analytical = ou_fht_pdf_vec_fixed_zero_branch(t, pars);
    for (int i = 0; i < n; ++i) {
      out[i] = analytical[i];
    }
    return out;
  }
  */
  double theta_max = std::exp(pars.t_scaled) - 1.0;

  int N = calculate_num_steps(theta_max, steps_fineness, min_steps);
  Rcout<<"N: "<<N<<std::endl;
  // Chunking layout and early-time prelude options
  ChunkingOptions chunk_opts;
  chunk_opts.ratio = chunk_ratio;
  chunk_opts.base_panels = chunk_base_panels;
  chunk_opts.max_chunks = chunk_max;
  if (R_finite(t_min_pos) && t_min_pos > 0.0) {
    chunk_opts.theta_min = std::exp(lambda * t_min_pos) - 1.0;
    int P = calculate_num_steps(chunk_opts.theta_min, min_t_fineness, (int)min_steps);
    if (P % 2 != 0) ++P;
    chunk_opts.pre_min_panels = std::max(2, P);
  }
  OU_DEBUG_LOG(1, "ou_fht_pdf_vec_grid_chunked: t_max=" << t_max
                  << ", theta_max=" << theta_max
                  << ", N=" << N
                  << ", fixed_b=" << (pars.fixed_b?1:0)
                  << ", chunk_ratio=" << chunk_ratio
                  << ", chunk_base_panels=" << chunk_base_panels
                  << ", chunk_max=" << chunk_max
                  << ", t_min_pos=" << (R_finite(t_min_pos)?t_min_pos:NA_REAL));

  // Prebuild the chunked theta grid to align beta-spline; no horizon clamp
  std::vector<double> theta_grid_pre = build_chunked_theta_grid(theta_max, N, chunk_opts);
  
  KernelFn kernel_fn;
  ForcingFn g_fn;
  std::unique_ptr<util::AkimaSpline> beta_spline;

  if (pars.fixed_b) {
    kernel_fn = KernelFn([](double tt, double ttp, const RD_Params& p) {
      return kernel_forward_tv_core(tt, ttp, p);
    });
    g_fn = ForcingFn([](double theta_val, const RD_Params& p) {
      return g_term_forward(theta_val, p);
    });
  } else {
    // Build beta spline directly on the chunked theta grid
	// Evaluate beta once at each grid point -- slow but O(N) and saves to cache
    beta_cache = make_boundary_decay_cache(theta_grid_pre, pars);
    cache_ptr = &beta_cache;
	// Iterate over the same grid to get a vector of pre-cached beta_vals to create spline
	std::vector<double> beta_vals(theta_grid_pre.size());
    for (size_t j = 0; j < theta_grid_pre.size(); ++j) {
      beta_vals[j] = beta_from_theta(theta_grid_pre[j], pars, cache_ptr);
    }
    beta_spline = std::make_unique<util::AkimaSpline>(theta_grid_pre, beta_vals);
    // Kernel and forcing using spline; diagonal uses cached beta'
    kernel_fn = KernelFn([cache_ptr,&beta_spline](double tt, double ttp, const RD_Params& p) {
      return kernel_forward_tv_core(tt, ttp, p, cache_ptr, beta_spline.get());
    });
    g_fn = ForcingFn([cache_ptr,&beta_spline](double theta_val, const RD_Params &p) {
      return g_term_forward(theta_val, p, cache_ptr,beta_spline.get());
    });
  }

  AbelFn abel_fn_theta = AbelFn([cache_ptr,&beta_spline](double theta_val, const RD_Params& p) {
        return abel_approx_nu_f(theta_val, p, p.fixed_b ? nullptr : cache_ptr,p.fixed_b ? nullptr : beta_spline.get());
  });

  std::vector<double> theta_grid_block;
  // Beta spline already prepared; still compute a beta_prime diagnostic at the end
  const double beta_prime = beta_prime_at_theta(theta_max, pars, theta_max); // TODO replace this spline derivative with the precomputed accurate one from beta_prime_diag_ptr
  auto nu_f_vals = solve_nu_block_with_abel_chunked(theta_max, pars, N, theta_grid_block,
                                                    kernel_fn, g_fn, abel_fn_theta,
                                                    chunk_opts);
  util::AkimaSpline nu_spline(theta_grid_block, nu_f_vals);

  const int total_nodes = static_cast<int>(theta_grid_block.size());

  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti) || ti <= 0.0) {
      out[i] = R_finite(ti) ? 0.0 : NA_REAL;
      continue;
    }

    RD_Params pars_i = prepare_ou_params(ti, lambda, theta, sigma, z0, b0, binf, tau, pow);
    const double theta_eval = std::exp(pars_i.t_scaled) - 1.0;
    const double et  = std::exp(pars_i.t_scaled);
    const double e2t = et * et;
    const double tau_eval   = 0.5 * (e2t - 1.0);

    const double nu_t = nu_spline.interpolate(theta_eval);
    const double scale = 1.0 + theta_eval;

    const double beta_prime_i = beta_prime_at_theta(theta_eval, pars_i, theta_max); // TODO this should be pre-computed why are we recalculating it?
    const double beta_t = pars_i.fixed_b ? (1.0 + theta_eval) * pars_i.b_scaled
                                          : beta_spline->interpolate(theta_eval);

    double nu_prime_local = 0.0;
    if (total_nodes >= 3) {
      auto upper_idx = std::upper_bound(theta_grid_block.begin(), theta_grid_block.end(), theta_eval);
      int jpanel = static_cast<int>(upper_idx - theta_grid_block.begin()) - 1;
      if (jpanel < 1) jpanel = 1;
      if (jpanel > total_nodes - 2) jpanel = total_nodes - 2;
      int c_index = jpanel + 1;
      nu_prime_local = local_quad_derivative(theta_grid_block, nu_f_vals, c_index, theta_eval);
    }

    // Choose derivative override. Use right-end quadratic at the endpoint,
    // otherwise the local quadratic derivative around theta_eval.
    double nu_prime_override = nu_prime_local;
    const double last_theta_eval = theta_grid_block.back();
    if (std::abs(theta_eval - last_theta_eval) <= std::max(1e-12, 1e-12 * std::abs(last_theta_eval)) && total_nodes >= 3) {
      const double mid_theta = 0.5 * (theta_grid_block[total_nodes - 2] + last_theta_eval);
      RightQuad rq_end = right_end_quadratic(theta_grid_block, nu_f_vals, mid_theta);
      nu_prime_override = rq_end.nu_prime;
      OU_DEBUG_LOG(3, "endpoint derivative override: local=" << nu_prime_local
                      << ", right_end_quad=" << nu_prime_override);
    }

    // Adaptive M: proportional to evaluation horizon at ti
    int N_eval = calculate_num_steps(ti, steps_fineness, min_steps);
    int M_eval = std::max(80, std::min(240, 2 * N_eval));
    const double integral_sum = integrate_pdf_forward(
        theta_eval, nu_spline, pars_i, nullptr, M_eval, nu_prime_override,
        pars_i.fixed_b ? nullptr : beta_spline.get());
    const double termA = (1.0 / std::sqrt(M_PI)) * integral_sum;
    const double termB = -nu_t / std::sqrt(2.0 * M_PI * tau_eval);

    const double omega_i = pars_i.omega;
    const double beta_prime_term = omega_i * (-(beta_prime_i / scale) * nu_t);

    const double image_term = omega_i* (0.5* averaged_image_term(tau_eval, beta_t, pars_i));
    const double g_scaled = termA + termB + beta_prime_term + image_term;
    const double pdf = lambda * g_scaled * e2t;
    out[i] = std::max(0.0, pdf);
    const bool near_end = (std::abs(theta_eval - theta_grid_block.back()) <= std::max(1e-12, 1e-12 * std::abs(theta_grid_block.back())));
    if (OU_DEBUG_LEVEL >= 4 || near_end) {
      OU_DEBUG_LOG(1, "ou_fht_pdf_vec_grid_chunked eval: t=" << ti
                      << ", theta_eval=" << theta_eval
                      << ", M_eval=" << M_eval
                      << ", nu_t=" << nu_t
                      << ", nu_prime_used=" << nu_prime_override
                      << ", integral_sum=" << integral_sum
                      << ", termA=" << termA << ", termB=" << termB
                      << ", beta_prime_i=" << beta_prime_i
                      << ", image_term=" << image_term
                      << ", scale=" << scale << ", tau=" << tau_eval
                      << ", beta_t=" << beta_t << ", b_scaled=" << pars_i.b_scaled
                      << ", zL_scaled=" << pars_i.zL_scaled << ", zU_scaled=" << pars_i.zU_scaled
                      << ", pdf=" << out[i]);
    }
  }

  return out;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SECTION 5: CDF via chunked PDF integration
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline static double linear_interp(const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   const double xv) {
  const size_t n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return y[0];
  if (!(xv > x.front())) return y.front();
  if (!(xv < x.back())) return y.back();
  auto it = std::lower_bound(x.begin(), x.end(), xv);
  size_t idx = static_cast<size_t>(std::distance(x.begin(), it));
  if (idx == 0) return y[0];
  const double x1 = x[idx - 1];
  const double x2 = x[idx];
  const double y1 = y[idx - 1];
  const double y2 = y[idx];
  const double dx = x2 - x1;
  if (dx <= FPM_EPSILON) return y1;
  const double w = (xv - x1) / dx;
  return y1 + w * (y2 - y1);
}

// [[Rcpp::export]]
NumericVector ou_fht_cdf_vec_grid_chunked(NumericVector t, double lambda,
                                          double theta, double sigma, double z0,
                                          double b0, double binf, double tau,
                                          double pow,
                                          double steps_fineness = 0.005,
                                          double min_steps = 300,
                                          double min_t_fineness = 0.001,
                                          double chunk_ratio = 1.5,
                                          double chunk_base_panels = 500,
                                          double chunk_max = 12,
                                          double rt_resolution = 0.02) {
  const int n = t.size();
  NumericVector out(n, NA_REAL);

  // Parameter sanity and time scanning
  double t_max = 0.0;
  bool has_pos = false;
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (R_finite(ti)) {
      if (ti < 0.0) {
        out[i] = 0.0; // CDF at negative times is 0
      } else {
        if (ti > t_max) t_max = ti;
        has_pos = true;
      }
    }
  }

  if (!has_pos || t_max <= 0.0 || !(lambda > 0.0) || !(sigma > 0.0)) {
    // For non-negative finite times, CDF is 0 under invalid settings
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) out[i] = 0.0;
    }
    return out;
  }

    RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, pow);
    /*
    if (pars.fixed_b && std::abs(pars.b_scaled) <= FPM_EPSILON) {
      NumericVector closed_form = ou_fht_cdf_vec_fixed_zero_branch(t, pars);
      for (int i = 0; i < n; ++i) {
        out[i] = closed_form[i];
        if (R_finite(t[i]) && t[i] < 0.0) out[i] = 0.0;
      }
      return out;
    }
    */

  // Build evaluation grid
  std::vector<double> t_grid;
  const int steps = std::max(0, static_cast<int>(std::ceil(t_max / rt_resolution)));
  for (int k = 0; k <= steps; ++k) t_grid.push_back(k * rt_resolution);
  if (t_grid.empty() || (t_grid.back() + 1e-15) < t_max) t_grid.push_back(t_max);
  if (t_grid.front() > 0.0) t_grid.insert(t_grid.begin(), 0.0);

  // Evaluate PDF on the grid via the chunked forward solver
  NumericVector t_grid_nv(t_grid.begin(), t_grid.end());
  NumericVector pdf_grid_nv = ou_fht_pdf_vec_grid_chunked(t_grid_nv,
                                                          lambda, theta, sigma,
                                                          z0, b0, binf,
                                                          tau, pow,
                                                          steps_fineness, min_steps,min_t_fineness,
                                                          chunk_ratio, chunk_base_panels, chunk_max);
  OU_DEBUG_LOG(1, "ou_fht_cdf_vec_grid_chunked: t_max=" << t_max
                  << ", grid_points=" << t_grid.size());
  const size_t m = static_cast<size_t>(pdf_grid_nv.size());
  std::vector<double> pdf_grid(m, 0.0);
  for (size_t i = 0; i < m; ++i) {
    const double val = pdf_grid_nv[static_cast<R_xlen_t>(i)];
    pdf_grid[i] = (R_finite(val) && val > 0.0) ? val : 0.0;
  }

  // Trapezoidal cumulative integration to get CDF on grid
  std::vector<double> cdf_grid(m, 0.0);
  for (size_t i = 1; i < m; ++i) {
    const double dt = t_grid[i] - t_grid[i - 1];
    if (dt > 0.0) {
      const double avg = 0.5 * (pdf_grid[i] + pdf_grid[i - 1]);
      cdf_grid[i] = cdf_grid[i - 1] + dt * avg;
    } else {
      cdf_grid[i] = cdf_grid[i - 1];
    }
    if (cdf_grid[i] < cdf_grid[i - 1]) cdf_grid[i] = cdf_grid[i - 1];
    if (cdf_grid[i] > 1.0) cdf_grid[i] = 1.0;
  }

  // Prefer Akima spline interpolation for smooth CDF mapping (exact on-grid)
  util::AkimaSpline cdf_spline(t_grid, cdf_grid);
  
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti)) { out[i] = NA_REAL; continue; }
    if (ti <= 0.0) { out[i] = 0.0; continue; }
    out[i] = cdf_spline.interpolate(ti);
    OU_DEBUG_LOG(6, "ou_fht_cdf_vec_grid_chunked eval: t=" << ti << ", cdf=" << out[i]);
  }

  return out;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OU -> BM (time-changed) helper utilities for BB-based simulation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inline double ou_u_time(double t, double lambda, double sigma) {
  // U(t) = (sigma^2/(2*lambda)) * (exp(2 lambda t) - 1)
  return (sigma * sigma) * (std::exp(2.0 * lambda * t) - 1.0) / (2.0 * lambda);
}

inline double ou_t_from_u(double u, double lambda, double sigma) {
  // t(u) = (1/(2 lambda)) * log(1 + (2 lambda / sigma^2) u)
  const double arg = 1.0 + (2.0 * lambda / (sigma * sigma)) * std::max(u, 0.0);
  return (0.5 / lambda) * std::log(std::max(arg, 1e-300));
}

// Recursively refine hit time within [t0,t1] for an OU bridge by operating
// in the Brownian-bridge (u) domain with exact midpoint sampling.
static double locate_cross_time_ou_bridge(double t0, double t1,
                                          double Y0, double Y1,
                                          double L0, double L1,
                                          const RD_Params& pars,
                                          double lambda, double theta, double sigma,
                                          int depth = 0) {
  const int max_depth = 12;
  const double tol_t = 1e-9;
  const double u0 = ou_u_time(t0, lambda, sigma);
  const double u1 = ou_u_time(t1, lambda, sigma);
  const double d0 = Y0 - L0;
  const double d1 = Y1 - L1;

  // If endpoints straddle in transformed domain, use linear-in-u root
  if (d0 * d1 <= 0.0) {
    const double denom = (d0 - d1);
    double alpha_u = (std::abs(denom) > FPM_EPSILON) ? (d0 / denom) : 0.5;
    if (alpha_u < 0.0) alpha_u = 0.0;
    if (alpha_u > 1.0) alpha_u = 1.0;
    const double u_hit = u0 + alpha_u * (u1 - u0);
    const double t_hit = ou_t_from_u(u_hit, lambda, sigma);
    return std::min(std::max(t_hit, t0), t1);
  }

  if (depth >= max_depth || (t1 - t0) <= tol_t * (1.0 + t0)) {
    return 0.5 * (t0 + t1);
  }

  const double tm = 0.5 * (t0 + t1);
  const double um = ou_u_time(tm, lambda, sigma);
  // Midpoint of Brownian bridge in u-domain: Normal(mean=(Y0+Y1)/2, var=(u1-u0)/4)
  const double mean_m = 0.5 * (Y0 + Y1);
  const double sd_m = std::sqrt(std::max(0.25 * (u1 - u0), 0.0));
  const double Ym = mean_m + sd_m * R::rnorm(0.0, 1.0);
  const double Bm = evaluate_boundary_decay(tm, pars);
  const double em = std::exp(lambda * tm);
  const double Lm = em * (Bm - theta);

  const double dL0 = Y0 - L0;
  const double dLm = Ym - Lm;
  const double dL1 = Y1 - L1;

  if (dL0 * dLm <= 0.0) {
    return locate_cross_time_ou_bridge(t0, tm, Y0, Ym, L0, Lm, pars, lambda, theta, sigma, depth + 1);
  }
  if (dLm * dL1 <= 0.0) {
    return locate_cross_time_ou_bridge(tm, t1, Ym, Y1, Lm, L1, pars, lambda, theta, sigma, depth + 1);
  }

  // No immediate straddle: compute per-half crossing probabilities in u-domain
  const double p_left = std::exp(-2.0 * std::abs(dL0) * std::abs(dLm) / std::max(um - u0, FPM_EPSILON));
  const double p_right = std::exp(-2.0 * std::abs(dLm) * std::abs(dL1) / std::max(u1 - um, FPM_EPSILON));
  const double sum_p = p_left + p_right;
  if (sum_p <= 0.0) {
    // Degenerate; choose the half with larger relative approach to the barrier
    if (std::abs(dL0) + std::abs(dLm) < std::abs(dLm) + std::abs(dL1)) {
      return locate_cross_time_ou_bridge(t0, tm, Y0, Ym, L0, Lm, pars, lambda, theta, sigma, depth + 1);
    } else {
      return locate_cross_time_ou_bridge(tm, t1, Ym, Y1, Lm, L1, pars, lambda, theta, sigma, depth + 1);
    }
  }
  const double u_rand = R::runif(0.0, 1.0);
  if (u_rand < (p_left / sum_p)) {
    return locate_cross_time_ou_bridge(t0, tm, Y0, Ym, L0, Lm, pars, lambda, theta, sigma, depth + 1);
  } else {
    return locate_cross_time_ou_bridge(tm, t1, Ym, Y1, Lm, L1, pars, lambda, theta, sigma, depth + 1);
  }
}

// [[Rcpp::export]]
NumericVector simulate_ou_hit_times(int n,
                                    double lambda, double theta, double sigma,
                                    double z0, double b0, double binf,
                                    double tau = 1.0, double p = 1.0,
                                    double dt = 1e-3, double t_max = 10.0) {
  if (n <= 0) {
    stop("simulate_ou_hit_times: n must be positive.");
  }
  if (!R_finite(lambda) || !R_finite(theta) || !R_finite(sigma) ||
      !R_finite(z0) || !R_finite(b0) || !R_finite(binf) ||
      !R_finite(tau) || !R_finite(p)) {
    stop("simulate_ou_hit_times: parameters must be finite.");
  }
  if (lambda <= 0.0 || sigma <= 0.0 || dt <= 0.0 || t_max <= 0.0) {
    stop("simulate_ou_hit_times: invalid parameters.");
  }
  if (z0 < 0.0) {
    stop("simulate_ou_hit_times: z0 must be >= 0 to define Uniform(0, z0) start.");
  }

  RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, p);
  const bool use_fixed = pars.fixed_b;

  const int steps = std::max(1, static_cast<int>(std::ceil(t_max / dt)));
  std::vector<double> B(steps + 1);
  if (use_fixed) {
    std::fill(B.begin(), B.end(), b0);
  } else {
    for (int k = 0; k <= steps; ++k) {
      const double t_k = k * dt;
      B[k] = evaluate_boundary_decay(t_k, pars);
    }
  }

  const double lambda_dt = lambda * dt;
  const double exp_neg = std::exp(-lambda_dt);
  const double var_factor = (sigma * sigma) * (-std::expm1(-2.0 * lambda_dt)) / (2.0 * lambda);
  const double sd_step = std::sqrt(std::max(var_factor, 0.0));

  NumericVector T(n, NA_REAL);
  for (int j = 0; j < n; ++j) {
    const double X0 = (z0 > 0.0) ? R::runif(0.0, z0) : 0.0;
    const double B0 = B[0];

    if (std::abs(X0 - B0) <= FPM_EPSILON) {
      T[j] = 0.0;
      continue;
    }

    const double dir = (X0 < B0) ? 1.0 : -1.0; // 1 => approach from below, -1 => from above

    double X = X0;
    double t = 0.0;
    bool hit = false;

    for (int k = 0; k < steps; ++k) {
      const double X_prev = X;
      const double B_prev = B[k];

      X = theta + (X_prev - theta) * exp_neg + sd_step * R::rnorm(0.0, 1.0);
      const double B_next = B[k + 1];

      const double D_prev = X_prev - B_prev;
      const double D_next = X - B_next;

      if ((dir > 0.0 && D_prev < 0.0 && D_next >= 0.0) ||
          (dir < 0.0 && D_prev > 0.0 && D_next <= 0.0)) {
        const double num = B_prev - X_prev;
        const double den = (X - X_prev) - (B_next - B_prev);
        double alpha = 0.5;
        if (std::abs(den) > FPM_EPSILON) {
          alpha = num / den;
        }
        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;

        T[j] = t + alpha * dt;
        hit = true;
        break;
      }

      t += dt;
    }

    if (!hit) {
      T[j] = NA_REAL;
    }
  }

  return T;
}

// [[Rcpp::export]]
NumericVector simulate_ou_hit_times_bb(int n,
                                       double lambda, double theta, double sigma,
                                       double z0, double b0, double binf,
                                       double tau = 1.0, double p = 1.0,
                                       double dt = 1e-3, double t_max = 10.0,
                                       double p_tol = 1e-9,
                                       double eps_curv = 0.1,
                                       double adapt_factor = 32.0,
                                       bool adaptive = true) {
  if (n <= 0) {
    stop("simulate_ou_hit_times_bb: n must be positive.");
  }
  if (!R_finite(lambda) || !R_finite(theta) || !R_finite(sigma) ||
      !R_finite(z0) || !R_finite(b0) || !R_finite(binf) ||
      !R_finite(tau) || !R_finite(p)) {
    stop("simulate_ou_hit_times_bb: parameters must be finite.");
  }
  if (lambda <= 0.0 || sigma <= 0.0 || dt <= 0.0 || t_max <= 0.0) {
    stop("simulate_ou_hit_times_bb: invalid parameters.");
  }
  if (z0 < 0.0) {
    stop("simulate_ou_hit_times_bb: z0 must be >= 0 to define Uniform(0, z0) start.");
  }

  RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, p);
  const bool use_fixed = pars.fixed_b;

  // Adaptive stepping controls (user-tunable)
  if (!(p_tol > 0.0 && p_tol < 1.0) || !R_finite(p_tol)) p_tol = 1e-8;
  if (!(eps_curv >= 0.0) || !R_finite(eps_curv)) eps_curv = 0.1;
  if (!(adapt_factor >= 1.0) || !R_finite(adapt_factor)) adapt_factor = 64.0;
  // Auto-tighten curvature tolerance for time-varying boundaries
  if (!use_fixed) {
    eps_curv *= 0.5; // halve tolerance when B(t) is not fixed
  }
  const double dt_min_floor = 1e-8;   // absolute minimum step

  NumericVector T(n, NA_REAL);
  for (int j = 0; j < n; ++j) {
    const double X0 = (z0 > 0.0) ? R::runif(0.0, z0) : 0.0;
    double X = X0;
    double t = 0.0;
    bool hit = false;

    const double B0 = use_fixed ? b0 : evaluate_boundary_decay(0.0, pars);
    if (std::abs(X - B0) <= FPM_EPSILON) {
      T[j] = 0.0;
      continue;
    }

    if (adaptive) {
      // Main adaptive stepping loop
      while (t < t_max) {
        const double Bt = use_fixed ? b0 : evaluate_boundary_decay(t, pars);
        const double e_t = std::exp(lambda * t);
        const double Y0 = e_t * (X - theta);
        const double L0 = e_t * (Bt - theta);
        const double d = std::abs(L0 - Y0);

        // Candidate coarse step bound
        const double dt_max_coarse = std::min(t_max - t, adapt_factor * dt);
        if (dt_max_coarse <= dt_min_floor) break; // out of time

        // Δu implied by max coarse step
        const double U_t = ou_u_time(t, lambda, sigma);
        const double U_cap = ou_u_time(t + dt_max_coarse, lambda, sigma) - U_t;

        // Δu target from p_tol (if d==0, back off to tiny Δu)
        double Du_prob = (d > 0.0) ? (2.0 * d * d / std::max(std::log(1.0 / p_tol), 1e-12)) : (1e-12);
        if (!R_finite(Du_prob) || Du_prob <= 0.0) Du_prob = 1e-12;
        double Du = std::min(Du_prob, std::max(U_cap, 0.0));

        // Map Δu -> Δt and clamp within [dt_min, dt_max_coarse]
        const double dt_min_adapt = std::max(dt_min_floor, dt / 32.0);
        double t_next = ou_t_from_u(U_t + Du, lambda, sigma);
        if (!R_finite(t_next)) t_next = t + dt_max_coarse;
        double dt_step = std::min(std::max(t_next - t, dt_min_adapt), dt_max_coarse);
        t_next = t + dt_step;

        // Curvature check in L(t) = e^{λ t} (B(t) - θ)
        for (int shrink = 0; shrink < 10; ++shrink) {
          const double tm = 0.5 * (t + t_next);
          const double Lt = L0;
          const double e_t1 = std::exp(lambda * t_next);
          const double B1 = use_fixed ? b0 : evaluate_boundary_decay(t_next, pars);
          const double L1 = e_t1 * (B1 - theta);
          const double e_tm = std::exp(lambda * tm);
          const double Bm = use_fixed ? b0 : evaluate_boundary_decay(tm, pars);
          const double Lm = e_tm * (Bm - theta);
          const double Llin = 0.5 * (Lt + L1);
          const double curv = std::abs(Lm - Llin);
          const double Du_step = std::max(ou_u_time(t_next, lambda, sigma) - U_t, 0.0);
          const double scale = std::max(d, std::sqrt(Du_step));
          if (curv <= eps_curv * scale) {
            // Curvature ok; proceed
            break;
          }
          // Too curved: shrink step
          dt_step *= 0.5;
          if (dt_step <= dt_min_adapt) { dt_step = dt_min_adapt; break; }
          t_next = t + dt_step;
        }

        // Exact OU update to t_next
        const double lambda_dt = lambda * dt_step;
        const double phi = std::exp(-lambda_dt);
        const double var_step = (sigma * sigma) * (-std::expm1(-2.0 * lambda_dt)) / (2.0 * lambda);
        const double sd_step = std::sqrt(std::max(var_step, 0.0));
        const double X_prev = X;
        X = theta + (X_prev - theta) * phi + sd_step * R::rnorm(0.0, 1.0);
        const double B_next = use_fixed ? b0 : evaluate_boundary_decay(t_next, pars);

        const double D_prev = X_prev - Bt;
        const double D_next = X - B_next;

        // Endpoint straddle: linear interpolation for D across the step
        if ((D_prev <= 0.0 && D_next >= 0.0) || (D_prev >= 0.0 && D_next <= 0.0)) {
          const double denom = (D_prev - D_next);
          double alpha = (std::abs(denom) > FPM_EPSILON) ? (D_prev / denom) : 0.5;
          if (alpha < 0.0) alpha = 0.0;
          if (alpha > 1.0) alpha = 1.0;
          T[j] = t + alpha * dt_step;
          hit = true;
          break;
        }

        // Same-side endpoints: OU->BM Brownian-bridge crossing probability
        if (D_prev * D_next > 0.0) {
          const double e_tn = std::exp(lambda * t_next);
          const double Y1 = e_tn * (X - theta);
          const double L1 = e_tn * (B_next - theta);
          const double d0 = L0 - Y0;
          const double d1 = L1 - Y1;
          const double Du_step = std::max(ou_u_time(t_next, lambda, sigma) - U_t, FPM_EPSILON);
          const double p_cross = std::exp(-2.0 * std::abs(d0) * std::abs(d1) / Du_step);
          if (R::runif(0.0, 1.0) < p_cross) {
            const double t_hit = locate_cross_time_ou_bridge(
                t, t_next,
                Y0, Y1,
                L0, L1,
                pars, lambda, theta, sigma, 0);
            T[j] = t_hit;
            hit = true;
            break;
          }
        }

        // Advance
        t = t_next;
      }
    } else {
      // Fixed-step loop with dt, using OU->BM crossing probability and recursive refinement
      const int steps = std::max(1, static_cast<int>(std::ceil(t_max / dt)));
      std::vector<double> B(steps + 1);
      if (use_fixed) {
        std::fill(B.begin(), B.end(), b0);
      } else {
        for (int k = 0; k <= steps; ++k) {
          const double t_k = k * dt;
          B[k] = evaluate_boundary_decay(t_k, pars);
        }
      }

      const double lambda_dt = lambda * dt;
      const double phi = std::exp(-lambda_dt);
      const double var_step = (sigma * sigma) * (-std::expm1(-2.0 * lambda_dt)) / (2.0 * lambda);
      const double sd_step = std::sqrt(std::max(var_step, 0.0));
      std::vector<double> e1(steps + 1);
      std::vector<double> U(steps + 1);
      for (int k = 0; k <= steps; ++k) {
        const double tk = k * dt;
        const double el = std::exp(lambda * tk);
        e1[k] = el;
        U[k] = (sigma * sigma) * (el * el - 1.0) / (2.0 * lambda);
      }

      const double B0 = B[0];
      if (std::abs(X - B0) <= FPM_EPSILON) {
        T[j] = 0.0;
        continue;
      }

      for (int k = 0; k < steps; ++k) {
        const double X_prev = X;
        const double B_prev = B[k];
        X = theta + (X_prev - theta) * phi + sd_step * R::rnorm(0.0, 1.0);
        const double B_next = B[k + 1];

        const double D_prev = X_prev - B_prev;
        const double D_next = X - B_next;

        // Straddle across endpoints: linear interpolation in D
        if ((D_prev <= 0.0 && D_next >= 0.0) || (D_prev >= 0.0 && D_next <= 0.0)) {
          const double denom = (D_prev - D_next);
          double alpha = (std::abs(denom) > FPM_EPSILON) ? (D_prev / denom) : 0.5;
          if (alpha < 0.0) alpha = 0.0;
          if (alpha > 1.0) alpha = 1.0;
          T[j] = t + alpha * dt;
          hit = true;
          break;
        }

        // Same-side: BB crossing in u-domain + within-step refinement
        if (D_prev * D_next > 0.0) {
          const double Y_prev = e1[k] * (X_prev - theta);
          const double Y_next = e1[k + 1] * (X - theta);
          const double L_prev = e1[k] * (B_prev - theta);
          const double L_next = e1[k + 1] * (B_next - theta);
          const double d0 = L_prev - Y_prev;
          const double d1 = L_next - Y_next;
          const double Du = std::max(U[k + 1] - U[k], FPM_EPSILON);
          const double p_cross = std::exp(-2.0 * std::abs(d0) * std::abs(d1) / Du);
          if (R::runif(0.0, 1.0) < p_cross) {
            const double t0 = t;
            const double t1 = t + dt;
            const double t_hit = locate_cross_time_ou_bridge(
                t0, t1,
                Y_prev, Y_next,
                L_prev, L_next,
                pars, lambda, theta, sigma, 0);
            T[j] = t_hit;
            hit = true;
            break;
          }
        }

        t += dt;
      }
    }

    if (!hit) {
      T[j] = NA_REAL;
    }
  }

  return T;
}

inline double gompertz_log_boundary_decay(double t, const RD_Params& pars) {
  const double cap = exp_decay_scalar(t, pars.b0, pars.binf, pars.tau, pars.pow);
  if (!(cap > 0.0)) {
    stop("simulate_gompertz_hit_times_bb: carrying capacity must stay positive.");
  }
  return std::log(cap);
}

inline double gompertz_log_fixed_boundary_decay(double, const RD_Params& pars) {
  if (!(pars.b0 > 0.0)) {
    stop("simulate_gompertz_hit_times_bb: carrying capacity must stay positive.");
  }
  return std::log(pars.b0);
}

inline RD_Params prepare_gompertz_log_boundary_params(double k0,
                                                      double kinf,
                                                      double tau,
                                                      double pow) {
  RD_Params pars{};
  pars.b0 = k0;
  pars.binf = kinf;
  pars.tau = tau;
  pars.pow = pow;
  pars.fixed_b = (std::abs(kinf - k0) <= FPM_EPSILON) || tau <= 0.0 || pow <= 0.0;
  if (pars.fixed_b) {
    pars.boundary_fn = BoundaryDecayFn(gompertz_log_fixed_boundary_decay);
  } else {
    pars.boundary_fn = BoundaryDecayFn(gompertz_log_boundary_decay);
  }
  return pars;
}

// [[Rcpp::export]]
NumericVector simulate_gompertz_hit_times_bb(int n,
                                             double alpha,
                                             double beta,
                                             double z0,
                                             double k0,
                                             double kinf,
                                             double tau = 1.0,
                                             double pow = 1.0,
                                             double dt = 1e-3,
                                             double t_max = 10.0,
                                             double start_floor = 1e-3,
                                             double p_tol = 1e-9,
                                             double eps_curv = 0.1,
                                             double adapt_factor = 32.0,
                                             bool adaptive = true) {
  if (n <= 0) {
    stop("simulate_gompertz_hit_times_bb: n must be positive.");
  }
  if (!R_finite(alpha) || !R_finite(beta) || !R_finite(z0) ||
      !R_finite(k0) || !R_finite(kinf) || !R_finite(tau) ||
      !R_finite(pow) || !R_finite(dt) || !R_finite(t_max) ||
      !R_finite(start_floor)) {
    stop("simulate_gompertz_hit_times_bb: parameters must be finite.");
  }
  if (alpha <= 0.0 || beta <= 0.0 || dt <= 0.0 || t_max <= 0.0) {
    stop("simulate_gompertz_hit_times_bb: invalid parameters.");
  }
  if (start_floor <= 0.0) {
    stop("simulate_gompertz_hit_times_bb: start_floor must be positive and finite.");
  }
  if (z0 <= start_floor) {
    stop("simulate_gompertz_hit_times_bb: z0 must exceed start_floor.");
  }
  if (k0 <= 0.0 || kinf <= 0.0) {
    stop("simulate_gompertz_hit_times_bb: carrying capacity must stay positive.");
  }

  RD_Params log_boundary_pars = prepare_gompertz_log_boundary_params(k0, kinf, tau, pow);
  const bool use_fixed = log_boundary_pars.fixed_b;

  if (!(p_tol > 0.0 && p_tol < 1.0) || !R_finite(p_tol)) p_tol = 1e-8;
  if (!(eps_curv >= 0.0) || !R_finite(eps_curv)) eps_curv = 0.1;
  if (!(adapt_factor >= 1.0) || !R_finite(adapt_factor)) adapt_factor = 64.0;
  if (!use_fixed) {
    eps_curv *= 0.5;
  }
  const double dt_min_floor = 1e-8;
  const double beta2 = beta * beta;
  const double lambda = alpha;
  const double sigma = beta;
  const double log_k_inf = std::log(kinf);
  const double theta = log_k_inf - (beta2 / (2.0 * alpha));
  const double log_tol = std::max(std::log(1.0 / p_tol), 1e-12);

  NumericVector T(n, NA_REAL);
  for (int j = 0; j < n; ++j) {
    double X0 = (z0 > start_floor) ? R::runif(start_floor, z0) : start_floor;
    if (!R_finite(X0) || X0 <= 0.0) {
      X0 = start_floor;
    }
    double Y = std::log(X0);
    double t = 0.0;
    bool hit = false;

    const double logK0 = evaluate_boundary_decay(0.0, log_boundary_pars);
    if (std::abs(Y - logK0) <= FPM_EPSILON) {
      T[j] = 0.0;
      continue;
    }

    if (adaptive) {
      while (t < t_max) {
        const double logKt = use_fixed ? std::log(k0) : evaluate_boundary_decay(t, log_boundary_pars);
        const double e_t = std::exp(lambda * t);
        const double Y0 = e_t * (Y - theta);
        const double L0 = e_t * (logKt - theta);
        const double d = std::abs(L0 - Y0);

        const double dt_max_coarse = std::min(t_max - t, adapt_factor * dt);
        if (dt_max_coarse <= dt_min_floor) break;

        const double U_t = ou_u_time(t, lambda, sigma);
        const double U_cap = ou_u_time(t + dt_max_coarse, lambda, sigma) - U_t;
        double Du_prob = (d > 0.0) ? (2.0 * d * d / std::max(std::log(1.0 / p_tol), 1e-12)) : 1e-12;
        if (!R_finite(Du_prob) || Du_prob <= 0.0) Du_prob = 1e-12;
        double Du = std::min(Du_prob, std::max(U_cap, 0.0));

        const double dt_min_adapt = std::max(dt_min_floor, dt / 32.0);
        double t_next = ou_t_from_u(U_t + Du, lambda, sigma);
        if (!R_finite(t_next)) t_next = t + dt_max_coarse;
        double dt_step = std::min(std::max(t_next - t, dt_min_adapt), dt_max_coarse);
        t_next = t + dt_step;

        double logK_next = use_fixed ? std::log(k0) : evaluate_boundary_decay(t_next, log_boundary_pars);

        for (int shrink = 0; shrink < 10; ++shrink) {
          const double tm = 0.5 * (t + t_next);
          const double Lt = L0;
          const double e_t1 = std::exp(lambda * t_next);
          const double L1 = e_t1 * (logK_next - theta);
          const double e_tm = std::exp(lambda * tm);
          const double logKm = use_fixed ? std::log(k0) : evaluate_boundary_decay(tm, log_boundary_pars);
          const double Lm = e_tm * (logKm - theta);
          const double Llin = 0.5 * (Lt + L1);
          const double curv = std::abs(Lm - Llin);
          const double Du_step = std::max(ou_u_time(t_next, lambda, sigma) - U_t, 0.0);
          const double scale = std::max(d, std::sqrt(std::max(Du_step, 0.0)));
          if (curv <= eps_curv * scale) {
            break;
          }
          dt_step *= 0.5;
          if (dt_step <= dt_min_adapt) {
            dt_step = dt_min_adapt;
            t_next = t + dt_step;
            logK_next = use_fixed ? std::log(k0) : evaluate_boundary_decay(t_next, log_boundary_pars);
            break;
          }
          t_next = t + dt_step;
          logK_next = use_fixed ? std::log(k0) : evaluate_boundary_decay(t_next, log_boundary_pars);
        }

        const double lambda_dt = lambda * dt_step;
        const double phi = std::exp(-lambda_dt);
        const double var_step = (sigma * sigma) * (-std::expm1(-2.0 * lambda_dt)) / (2.0 * lambda);
        const double sd_step = std::sqrt(std::max(var_step, 0.0));
        const double Y_prev = Y;
        Y = theta + (Y_prev - theta) * phi + sd_step * R::rnorm(0.0, 1.0);

        const double D_prev = Y_prev - logKt;
        const double D_next = Y - logK_next;

        if ((D_prev <= 0.0 && D_next >= 0.0) || (D_prev >= 0.0 && D_next <= 0.0)) {
          const double denom = (D_prev - D_next);
          double alpha_dt = (std::abs(denom) > FPM_EPSILON) ? (D_prev / denom) : 0.5;
          if (alpha_dt < 0.0) alpha_dt = 0.0;
          if (alpha_dt > 1.0) alpha_dt = 1.0;
          T[j] = t + alpha_dt * dt_step;
          hit = true;
          break;
        }

        if (D_prev * D_next > 0.0) {
          const double e_tn = std::exp(lambda * t_next);
          const double Y1 = e_tn * (Y - theta);
          const double L1 = e_tn * (logK_next - theta);
          const double d0 = L0 - Y0;
          const double d1 = L1 - Y1;
          const double Du_step = std::max(ou_u_time(t_next, lambda, sigma) - U_t, FPM_EPSILON);
          const double p_cross = std::exp(-2.0 * std::abs(d0) * std::abs(d1) / Du_step);
          if (R::runif(0.0, 1.0) < p_cross) {
            const double t_hit = locate_cross_time_ou_bridge(
                t, t_next,
                Y0, Y1,
                L0, L1,
                log_boundary_pars, lambda, theta, sigma, 0);
            T[j] = t_hit;
            hit = true;
            break;
          }
        }

        t = t_next;
      }
    } else {
      const int steps = std::max(1, static_cast<int>(std::ceil(t_max / dt)));
      std::vector<double> logK(steps + 1);
      if (use_fixed) {
        std::fill(logK.begin(), logK.end(), std::log(k0));
      } else {
        for (int k = 0; k <= steps; ++k) {
          const double t_k = k * dt;
          logK[k] = evaluate_boundary_decay(t_k, log_boundary_pars);
        }
      }

      const double lambda_dt = lambda * dt;
      const double phi = std::exp(-lambda_dt);
      const double var_step = (sigma * sigma) * (-std::expm1(-2.0 * lambda_dt)) / (2.0 * lambda);
      const double sd_step = std::sqrt(std::max(var_step, 0.0));
      std::vector<double> e1(steps + 1);
      std::vector<double> U(steps + 1);
      for (int k = 0; k <= steps; ++k) {
        const double tk = k * dt;
        const double el = std::exp(lambda * tk);
        e1[k] = el;
        U[k] = (sigma * sigma) * (el * el - 1.0) / (2.0 * lambda);
      }

      const double logB0 = logK[0];
      if (std::abs(Y - logB0) <= FPM_EPSILON) {
        T[j] = 0.0;
        continue;
      }

      double t_cur = 0.0;
      for (int k = 0; k < steps; ++k) {
        const double Y_prev = Y;
        const double logB_prev = logK[k];
        Y = theta + (Y_prev - theta) * phi + sd_step * R::rnorm(0.0, 1.0);
        const double logB_next = logK[k + 1];

        const double D_prev = Y_prev - logB_prev;
        const double D_next = Y - logB_next;

        if ((D_prev <= 0.0 && D_next >= 0.0) || (D_prev >= 0.0 && D_next <= 0.0)) {
          const double denom = (D_prev - D_next);
          double alpha_dt = (std::abs(denom) > FPM_EPSILON) ? (D_prev / denom) : 0.5;
          if (alpha_dt < 0.0) alpha_dt = 0.0;
          if (alpha_dt > 1.0) alpha_dt = 1.0;
          T[j] = t_cur + alpha_dt * dt;
          hit = true;
          break;
        }

        if (D_prev * D_next > 0.0) {
          const double Y_prev_u = e1[k] * (Y_prev - theta);
          const double Y_next_u = e1[k + 1] * (Y - theta);
          const double L_prev = e1[k] * (logB_prev - theta);
          const double L_next = e1[k + 1] * (logB_next - theta);
          const double d0 = L_prev - Y_prev_u;
          const double d1 = L_next - Y_next_u;
          const double Du = std::max(U[k + 1] - U[k], FPM_EPSILON);
          const double p_cross = std::exp(-2.0 * std::abs(d0) * std::abs(d1) / Du);
          if (R::runif(0.0, 1.0) < p_cross) {
            const double t0 = t_cur;
            const double t1 = t_cur + dt;
            const double t_hit = locate_cross_time_ou_bridge(
                t0, t1,
                Y_prev_u, Y_next_u,
                L_prev, L_next,
                log_boundary_pars, lambda, theta, sigma, 0);
            T[j] = t_hit;
            hit = true;
            break;
          }
        }

        t_cur += dt;
      }
    }

    if (!hit) {
      T[j] = NA_REAL;
    }
  }

  return T;
}

#endif // OU_HITTING_TIME_H

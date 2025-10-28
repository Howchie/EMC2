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
#include <utility>
#include <Rcpp.h>
#include <quadmath.h>
using namespace Rcpp;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Section 1: Core Constants and Helper Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inline double averaged_shifted_gaussian(double beta, double variance, double z_lo, double z_hi) {
  if (variance <= FPM_EPSILON) {
    return 0.0;
  }
  if (z_hi < z_lo) {
    std::swap(z_lo, z_hi);
  }
  const double span = z_hi - z_lo;
  if (span <= FPM_EPSILON) {
    const double delta = z_lo - beta;
    return delta * safe_exp(-(delta * delta) / variance);
  }
  const double half_var = 0.5 * variance;
  const double scale = std::sqrt(2.0 * M_PI * half_var);
  const double exp_hi = Gstar(half_var, z_hi - beta) * scale;
  const double exp_lo = Gstar(half_var, z_lo - beta) * scale;
  return (variance / (2.0 * span)) * (exp_hi - exp_lo);
}

inline RD_Params prepare_ou_params(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double pow,
                                 BoundaryDecayFn boundary_fn = BoundaryDecayFn(),
                                 std::vector<double> boundary_params = {}) {
  RD_Params p{};
    p.sp_var = false;
    p.b0 = b0;
    p.c = std::sqrt(lambda) / sigma;
    p.t_scaled = lambda * t;

    const double z_scaled_raw = p.c * (z0 - theta);
    const double b_scaled_raw = p.c * (b0 - theta);
    const double lower_scaled_raw = p.c * (0.0 - theta);

    p.omega = (z_scaled_raw >= b_scaled_raw) ? 1.0 : -1.0;

    p.zU_scaled = z_scaled_raw; // these will both be zero for point start
    p.zL_scaled = p.zU_scaled;
    if (p.sp_var) {
        const double lo = std::min(lower_scaled_raw, z_scaled_raw);
        const double hi = std::max(lower_scaled_raw, z_scaled_raw);
        p.zL_scaled = lo;
        p.zU_scaled = hi;
    }
    p.b_scaled = b_scaled_raw;
    p.binf = binf;
    p.fixed_b = true;
    p.sigma = sigma;
    p.z0 = z0;
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
    const double bt = -pars.omega*evaluate_boundary_decay(t_max - t, pars);
    const double bt_scaled =  (pars.c * (bt - pars.theta));
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
    const double bt = -pars.omega * evaluate_boundary_decay(t, pars);
    const double bt_scaled = pars.c * (bt - pars.theta);
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

inline void effective_positions_tv(double v, double vp, double variance, const RD_Params& pars,
                                double& Delta, double& one_minus_vp,
                                const BoundaryDecayCache* cache = nullptr) {
  const double alpha = std::max(0.0, 1.0 - v);
  const double omega = pars.omega;
  double beta_vp;
  if (pars.fixed_b) {
        beta_vp = pars.b_scaled * std::max(0.0, 1.0 - vp);
  } else {
        beta_vp = beta_from_v(vp, pars, cache);
  }
  const double z_start = pars.zU_scaled;
  Delta = alpha * z_start - beta_vp;
  if (pars.sp_var) {
      double z_lo = pars.zU_scaled, z_hi = pars.zU_scaled;
      if ((pars.zU_scaled - pars.zL_scaled > FPM_EPSILON) && variance > FPM_EPSILON) {
          const double y_lo = alpha * z_lo;
          const double y_hi = alpha * z_hi;
          const double avg = averaged_shifted_gaussian(beta_vp, variance, y_lo, y_hi);
          Delta = -avg;
      }
  }
  one_minus_vp = std::max(0.0, 1.0 - vp);
}

// KERNELS
// The non-singular part of the backwards Volterra kernel in v from page 11 Lipton & Kaushansky(2018) 
double kernel_backward(double v, double vp, const RD_Params& pars) {
    const double b = pars.b_scaled;
    if (std::abs(b) < FPM_EPSILON) {
        return 0.0;
    }
    const double omega = pars.omega;
    const double beta_v = beta_from_v(v, pars);
    const double beta_vp = beta_from_v(vp, pars);
	const double dv = v - vp;
    const double s = 2.0 - v - vp;
	if (dv < 0.0 || s <= FPM_EPSILON) return 0.0;
    const double expo = safe_exp(-b*b*dv/s); // clamp extreme values
    const double term1 = (-2.0 * b / std::sqrt(M_PI)); // minus sign because b here is implicilty negative. In fixed b case, beta_v<beta_vp, and (beta_v - beta_vp)/dv = -b. In Lipton & Kaushansky (2020) they removed omega and assumed a hitting from above problem which canceled some signs, this is the more general form
    const double term2 = (((1.0 - vp) * expo) / (s * std::sqrt(s)));
    const double kernel = omega * (term1 * term2);
    return -kernel; // - because we move the kernel term across to the right hand side when solving (Eq 13 Lipton & Kaushansky 2020)
}

// Time-varying extension of above. replace b with (beta(v) - beta(v'))/dv
inline double kernel_backward_tv_core(double v, double vp, const RD_Params& pars, const BoundaryDecayCache& cache) {
    const double dv = v - vp;
    if (dv <= FPM_EPSILON)
      return 0.0;
    const double beta_v = beta_from_v(v, pars, &cache);
    const double beta_vp = beta_from_v(vp, pars, &cache);
    const double psi = beta_v - beta_vp;
    if (std::abs(psi) < 1e-15) return 0.0;
    const double s = 2.0 - v - vp;
    if (s <= FPM_EPSILON) return 0.0;
    const double omega = pars.omega;
    const double b_eff = psi / dv;
    const double xi = safe_exp(-(psi * psi) / (dv * s));

    const double term1 = (-2.0 * b_eff / std::sqrt(M_PI));
    const double term2 = ((1.0 - vp) * xi) / (s * std::sqrt(s));
    const double kernel = omega * (term1 * term2);
    return -kernel;
}

// The non-singular part of the backwards Volterra kernel in theta from page 11 Lipton & Kaushansky(2018)
double kernel_forward(double theta, double theta_p, const RD_Params& pars) {
    const double b = pars.b_scaled;
    if (std::abs(b) < FPM_EPSILON) return 0.0;
    const double omega = pars.omega;
    const double dv = theta - theta_p;
    const double s = 2.0 + theta + theta_p;
	// Handle diagonal explicitly
    if (std::abs(dv) < 1e-15) {
        // K(tau,tau) = - b / (sqrt (2pi) sqrt s)
        if (s < FPM_EPSILON)
          return 0.0;
        const double kernel_limit = -omega*(pars.beta_prime / std::sqrt(2.0 * M_PI) / sqrt_pos((1+theta)));
        return kernel_limit;// in fixed case, beta' = -b_scaled
    }
    if (theta <= theta_p + FPM_EPSILON) return 0.0;

    const double ratio = dv / s;
    const double xi  = safe_exp(-b*b*ratio);
    const double term1 = 2.0*b*xi  / std::sqrt(M_PI);
    const double term2 = (1.0 + theta_p)/ (s * std::sqrt(s)); // change of variables transform and jacobian
    // Eq. 10 Lipton et al (2018), -sign when moving to RHS of equation
    const double kernel = omega * (term1 * term2);
    return -kernel; // -omega because we move the kernel term across to the right hand side when solving (Eq 13 Lipton & Kaushansky 2020)
}

// Time-varying extension of above. replace b with (beta(theta)-beta(theta_p))/dtheta
// Diagonal term uses the derivative of beta at theta. In the fixed case that term reduces to b_scaled, hence the change from above
double kernel_forward_tv_core(double theta, double theta_p, const RD_Params& pars, const BoundaryDecayCache& cache) {
    const double dv = theta - theta_p;
    const double scale = 1.0 + theta;
    const double omega = pars.omega;
    const double beta_theta = beta_from_theta(theta, pars, &cache);
    const double beta_theta_p = beta_from_theta(theta_p, pars, &cache);
    // Handle diagonal explicitly
    if (std::abs(dv) < 1e-15) {
        // Use the local derivative beta'(theta) for the diagonal limit.
        // For fixed barriers beta'(theta) = b_scaled (constant); for
        // time-varying barriers we use the precomputed derivative from
        // the boundary cache at this theta.
        double beta_prime_local = pars.b_scaled; // fixed-barrier fallback
        if (!pars.fixed_b) {
            const double bp_cached = cache.lookup_prime(theta);
            if (std::isfinite(bp_cached)) {
                beta_prime_local = bp_cached;
            } else {
                // As a conservative fallback, approximate via a small
                // centered difference on the cached beta values if available.
                // Since the cache is built on a uniform grid [0, theta_max]
                // with spacing h, we can sample neighboring nodes.
                const double h = cache.h;
                const double tL = std::max(0.0, theta - h);
                const double tR = theta + h;
                const double bL = beta_from_theta(tL, pars, &cache);
                const double bR = beta_from_theta(tR, pars, &cache);
                const double denom = (tR - tL);
                beta_prime_local = (denom > FPM_EPSILON) ? ((bR - bL) / denom) : 0.0;
            }
        }
        if (scale < FPM_EPSILON) return 0.0;
        return -omega*(beta_prime_local / std::sqrt(2.0 * M_PI)) / sqrt_pos(scale);
    }
    const double psi = beta_theta - beta_theta_p; // this is equivalent to b*dv in the fixed case
    const double s = 2.0 + theta + theta_p;
    if (s <= 1e-15) return 0.0;

    const double b_eff = psi / dv;
    const double xi = safe_exp(-(psi * psi) / (dv * s)); 

    const double term1 = (2.0 * b_eff  * xi/ std::sqrt(M_PI));
    const double term2 = ((1.0 + theta_p)) / (s * std::sqrt(s)); // change of variables transform and jacobian
    const double kernel = omega * (term1 * term2);
    return -kernel; // -omega because we move the kernel term across to the right hand side when solving (Eq 13 Lipton & Kaushansky 2020)
}

// Dummy function for backward kernel
double g_term_backward(double v, const RD_Params& pars) {
    return 1.0;
}

// Forward kernel g term from page 11 Lipton & Kaushansky(2018)
double g_term_forward(double theta, const RD_Params& pars) {
    const double z = pars.zU_scaled;
    const double b = pars.b_scaled;
    const double scale = 1.0 + theta;
	const double tau_2 = scale*scale-1.0; // this is 2*tau, which is the solver units
    if (theta <= FPM_EPSILON) return 0.0;
    const double num = z - scale*b;
	const double expo = safe_exp(-(num*num)/tau_2);
	const double denom = sqrt_pos(M_PI * tau_2);
	if (denom <= FPM_EPSILON) return 0.0;
    return -expo / denom;
}

double g_term_forward_uniform(double theta, const RD_Params& pars) {
    const double scale = 1.0 + theta;
    const double tau_2 = scale*scale-1.0; // this is 2*tau, which is the solver units
    if (theta <= FPM_EPSILON || tau_2 <= FPM_EPSILON) return 0.0;
    double z_lo = pars.zU_scaled, z_hi = pars.zU_scaled;
    if (!(pars.zU_scaled - pars.zL_scaled > FPM_EPSILON)) {
        return g_term_forward(theta, pars);
    }
    const double span = z_hi - z_lo;
    if (span <= FPM_EPSILON) {
        return g_term_forward(theta, pars);
    }
    const double a = scale * pars.b_scaled;
    const double root = std::sqrt(tau_2);
    const double erf_hi = std::erf((z_hi - a) / root);
    const double erf_lo = std::erf((z_lo - a) / root);
    return -(erf_hi - erf_lo) / (2.0 * span);
}

// Time-varying version of the above, replace b with beta at theta (beta includes the (1+theta) factor already)
double g_term_forward_tv(double theta, const RD_Params& pars,
                                        const BoundaryDecayCache& cache) {
    //if (pars.fixed_b) {
    //    return g_term_forward(theta, pars);
    //}
    //TODO implement the averaging across 0-z0 for sp_var case
    const double z = pars.zU_scaled;
    const double scale = 1.0 + theta;
	const double tau_2 = scale*scale-1.0; // this is 2*tau, which is the solver units
    if (theta <= FPM_EPSILON) return 0.0;
    const double beta = beta_from_theta(theta, pars, &cache); // beta already has "scale" factored in
    const double num = beta - z;
	const double expo = safe_exp(-(num*num)/tau_2);
	const double denom = sqrt_pos(M_PI * tau_2);
	if (denom <= FPM_EPSILON) return 0.0;
    return -expo / denom;
}

inline double g_term_forward_tv_uniform(double theta, const RD_Params& pars, const BoundaryDecayCache& cache) {
    const double scale = 1.0 + theta;
    const double tau_2 = scale*scale-1.0; // this is 2*tau, which is the solver units
    if (theta <= FPM_EPSILON || tau_2 <= FPM_EPSILON) return 0.0;
    double z_lo = pars.zL_scaled, z_hi = pars.zU_scaled;
    if (!(pars.zU_scaled - pars.zL_scaled > FPM_EPSILON)) {
        return g_term_forward_tv(theta, pars, cache);
    }
    const double span = z_hi - z_lo;
    if (span <= FPM_EPSILON) {
        return g_term_forward_tv(theta, pars, cache);
    }
    const double beta = beta_from_theta(theta, pars, &cache); // beta already has "scale" factored in
    const double root = std::sqrt(tau_2);
    const double erf_hi = std::erf((z_hi - beta) / root);
    const double erf_lo = std::erf((z_lo - beta) / root);
    return -(erf_hi - erf_lo) / (2.0 * span);
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
    return 2.0 * std::exp((b_signed * b_signed * v) / 2.0) * normal_cdf(b_signed * std::sqrt(v));
}

double abel_approx_nu_f(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    if (theta < FPM_EPSILON) {
        return 0.0; // nu_f(0) = g(0) = 0
    }
    if (pars.sp_var) {
        if (cache) {
            return g_term_forward_tv_uniform(theta, pars, *cache);
        }
        return g_term_forward_uniform(theta, pars);
    }
    double b_eff = pars.b_scaled;
    if (!pars.fixed_b) {
        const double beta_theta = beta_from_theta(theta, pars, cache);
        const double scale = std::max(FPM_EPSILON, 1.0 + theta);
        b_eff = beta_theta / scale;
    }
    const double omega = pars.omega;
    const double b_signed = omega * b_eff;
    const double z_signed = omega * pars.zU_scaled;
    // TODO: implement the uniform start point case here
    const double num = b_signed * theta + z_signed - b_signed;

    const double exp_arg1 = 0.5 * b_signed * b_signed * theta + b_signed * (z_signed - b_signed);
    const double cdf_arg = -num / std::sqrt(theta);
    const double term1 = b_signed * safe_exp(exp_arg1) * normal_cdf(cdf_arg);

    const double exp_arg2 = -((b_signed - z_signed) * (b_signed - z_signed)) / (2 * theta);
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
    effective_positions_tv(v, vp, dv_eff * s_eff, pars, Delta, one_minus_vp,
                           cache);
    const double expo = safe_exp(-Delta * Delta / (dv_eff * s_eff));

    const double denom = std::sqrt(M_PI) * (s_eff * std::sqrt(s_eff)); 
    const double num = Delta * expo * one_minus_vp * nu;
    if (denom < FPM_EPSILON)
      return 0.0;
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
    const double scale = 1.0 + theta;
    const double beta_theta = beta_from_theta(theta, pars, cache);
    const double beta_theta_p = beta_from_theta(theta_p, pars, cache);

    const double psi = beta_theta - beta_theta_p; // this is equivalent to b*dv in the fixed case
    const double psi_sq = psi * psi;
    const double s = 2.0 + theta + theta_p;
    if (s <= 1e-15) return 0.0;

    const double b_eff = psi / dv;
    const double ratio = (psi * psi) / (dv * s);
    const double xi = safe_exp(-ratio);

    const double term1 = (2.0 * b_eff * xi / std::sqrt(M_PI));

    if (dv <= 0) return 0.0;

    // Kernel's smooth part, H_smooth(theta, theta_p)
    //const double ratio = dv / s;
    const double H_smooth = (1.0 - 2.0 * ratio) * xi;
    const double jacobian_part = 1.0 + theta_p;

    const double denom = dv * s * sqrt_pos(s); // dv here = u^2
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

    // Evaluate f(u_N) = f(0) using the derivative limit
    const double mid = 0.5 * (theta_grid[N-1] + theta_max);
	  const RightQuad rq = right_end_quadratic(theta_grid, nu_f_vals, mid);
    const double nu_prime_theta = rq.nu_prime;
    const double s = 1.0 + theta_max;
    const double f_right_limit = -nu_prime_theta / std::sqrt(2.0 * s);
    Rcout << "f_right_limit: " << f_right_limit << std::endl;
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
        double u2 = (u_grid[i] * u_grid[i]);
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
double integrate_pdf_forward_spline(
    double theta_k,
    const util::AkimaSpline& nu_spline,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr,
    int M = 100) // 100 panels is probably overkill, but FAST.
{
    if (theta_k <= 0.0) return 0.0;

    // 1. DYNAMICALLY create a small, temporary grid for THIS integration
    std::vector<double> temp_theta_grid(M + 1);
    std::vector<double> temp_u_grid(M + 1);
    const double h_theta = theta_k / M; // Step size for *this* integration

    for (int i = 0; i <= M; ++i) {
        temp_theta_grid[i] = i * h_theta;
        // The u-transform is now RELATIVE to theta_k
        temp_u_grid[i] = std::sqrt(std::max(0.0, theta_k - temp_theta_grid[i]));
    }

    // Get the ANALYTIC ENDPOINT at theta_k using the spline
    const double nu_at_k = nu_spline.interpolate(theta_k);
    const double nu_prime_at_k = nu_spline.derivative(theta_k);
    const double s = 1.0 + theta_k;
    const double f_right_limit = -nu_prime_at_k / std::sqrt(2.0 * s);
    Rcout << "f_right_limit: " << f_right_limit << std::endl;
    // Run Simpson's loop
    double integral_sum = 0.0;
    
    // Panels are [i, i+1] in temp_theta (equiv. [u_{i+1}, u_i] in u)
    for (int i = M - 1; i >= 0; --i) {
        const double uL = temp_u_grid[i];     // u(theta_i)
        const double uR = temp_u_grid[i + 1]; // u(theta_{i+1})
        const double h_u = uL - uR;           // panel width in u

        if (h_u == 0.0) continue; // Should only happen on first panel if i=M

        const double thetaL = temp_theta_grid[i];
        const double thetaR = temp_theta_grid[i + 1];

        // Midpoint in u-space
        const double uM = 0.5 * (uL + uR);
        const double thetaU = theta_k - uM * uM; // Midpoint in theta-space

        // Get nu values from the spline
        const double nuL = nu_spline.interpolate(thetaL);
        const double nuR = nu_spline.interpolate(thetaR);
        const double nuU = nu_spline.interpolate(thetaU);

        // Get f at left endpoint
        const double fL = regularized_pdf_integrand_theta(theta_k, thetaL,
                                                nu_at_k, nuL, pars, cache);

        // Get f at right endpoint
        // (Use the exact limit *only* on the final panel, i == M-1)
        const double fR = (i == M - 1) ? f_right_limit : 
                        regularized_pdf_integrand_theta(theta_k, thetaR,
                                                nu_at_k, nuR, pars, cache);

        // Get f at midpoint
        const double fM = regularized_pdf_integrand_theta(theta_k, thetaU,
                                                nu_at_k, nuU, pars, cache);

        // Simpson on this u-panel
        integral_sum += (h_u / 6.0) * (fL + 4.0 * fM + fR);
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
    const double v_max = 1.0 - std::exp(-pars.t_scaled);
    pars.v_max = v_max;
    const double omega = pars.omega;
    const double b_at_t = (beta_from_v(v_max, pars) / (1 - v_max)) / pars.c + theta; // unscaled b at time t
    if (std::abs(b_at_t) < FPM_EPSILON) { // Already hit
         return 1.0;
    } 

    // TODO implement the sp_var case here
    if (std::abs(pars.b_scaled) < FPM_EPSILON) { // Analytical case for b_scaled = 0
        double sinh_t = std::sinh(pars.t_scaled);
        if (sinh_t < FPM_EPSILON) {
            return 0.0;
        } else {
            double common_arg = -pars.zU_scaled * std::exp(-0.5 * pars.t_scaled) / std::sqrt(sinh_t);
            return 2.0 * normal_cdf(common_arg);
        }
    }
  
    const int N = (num_steps < 2) ? 2 : num_steps;

    const double h_v = v_max / N;
    std::vector<double> v_grid(N + 1);
    for (int j = 0; j <= N; ++j) v_grid[j] = j * h_v;

    BoundaryDecayCache beta_cache;
    const BoundaryDecayCache* cache_ptr = nullptr;
    KernelFn kernel_fn;

    if (pars.fixed_b) {
        kernel_fn = KernelFn(kernel_backward);
    } else {
        beta_cache = make_boundary_decay_cache(v_max, N, pars, beta_from_v_raw);
        kernel_fn = KernelFn([&beta_cache](double vv, double vvp, const RD_Params& p) {
            return kernel_backward_tv_core(vv, vvp, p, beta_cache);
        });
        cache_ptr = &beta_cache;
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

double ou_fht_pdf_forward(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double p, int num_steps) {
    RD_Params pars = prepare_ou_params(t, lambda, theta, sigma, z0, b0, binf, tau, p);
    const double omega = pars.omega;
    if (t < FPM_EPSILON || lambda <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }

    if (std::abs(pars.b_scaled) < FPM_EPSILON) { // Analytical case for b_scaled = 0
        double sinh_t = std::sinh(pars.t_scaled);
        if (sinh_t < FPM_EPSILON) {
            return 0.0;
        } else {
            // TOD0: implement the sp_var case here
            double common_arg = -pars.zU_scaled * std::exp(-0.5 * pars.t_scaled) / std::sqrt(sinh_t);
            const double exp_term = std::exp(-0.5 * common_arg * common_arg + 0.5 * pars.t_scaled);
            const double g_scaled = (pars.zU_scaled * exp_term) / (std::sqrt(2.0 * M_PI) * std::pow(sinh_t, 1.5));
            return lambda * g_scaled;
        }
    }
	  const double et  = std::exp(pars.t_scaled);
	  const double e2t = std::exp(2.0 * pars.t_scaled);
    const double e2tm1 = e2t - 1.0;
    const double theta_max = et - 1.0;
    const double tau_max = 0.5 * (e2tm1);
    const int N = (num_steps < 2) ? 2 : num_steps;
    std::vector<double> theta_grid_block;
	
    pars.theta_max = theta_max;

    if (omega*(pars.zU_scaled - pars.b_scaled) < FPM_EPSILON) { 
         return 0.0;
    }
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
            return kernel_forward_tv_core(tt, ttp, p, beta_cache);
        });
        cache_ptr = &beta_cache;
    }

    ForcingFn g_fn;
    if (pars.fixed_b) {
        g_fn = pars.sp_var ? ForcingFn(g_term_forward_uniform)
                           : ForcingFn(g_term_forward);
    } else {
        g_fn = ForcingFn([cache_ptr](double theta_val, const RD_Params &p) {
                if (cache_ptr) {
                    return p.sp_var ? g_term_forward_tv_uniform(theta_val, p, *cache_ptr)
                                     : g_term_forward_tv(theta_val, p, *cache_ptr);
                }
                return p.sp_var ? g_term_forward_uniform(theta_val, p)
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
    const double beta_t = beta_grid.back();
    pars.beta_prime = beta_t; // constant-barrier default
    double mid = 0.5 * (t_grid[N] + t_grid.back()); // theta_grid is N+1 in size 
    RightQuad rq = right_end_quadratic(t_grid, beta_grid,mid); 
    pars.beta_prime = rq.nu_prime;

    auto nu_f_vals = solve_nu_block_with_abel(theta_max, pars, N, theta_grid_block, kernel_fn, g_fn, abel_fn_theta);
    const double nu_t = nu_f_vals.back();
    // Endpoint derivative of nu at theta_max via local right-end quadratic fit
    double nu_prime_at_max = 0.0;
    if (theta_grid_block.size() >= 3) {
        const double mid_nu = 0.5 * (theta_grid_block[N - 1] + theta_max);
        RightQuad rq_nu = right_end_quadratic(theta_grid_block, nu_f_vals, mid_nu);
        nu_prime_at_max = rq_nu.nu_prime;
    }
    // Equations here are re-arranged to match general form in Lipton &
    // Kaushansky (2020) Equation 10 They are mathematically equivalent to the
    // original form in section 4.4, but easier to relate to the general case
    // Term 1: The integral from Eq. (10)
    const double integal_sum = integrate_pdf_forward_theta_u(theta_max, theta_grid_block, nu_f_vals, pars, cache_ptr);
    const double termA = (1/std::sqrt(M_PI))*integal_sum;
    // Term 2: -ν(t) / sqrt(2πt)
    const double termB = -nu_t / std::sqrt(2.0 * M_PI * tau_max);

    // The b'(t)ν(t) part
    const double beta_prime_term = omega*(-(pars.beta_prime/scale) * nu_t);
    const double image_term = omega* (0.5* averaged_image_term(tau_max, beta_t, pars));
    const double g_scaled = termA + termB + beta_prime_term + image_term;
    const double pdf = lambda * g_scaled * e2t;
    Rcout<<"nu_prime=" << nu_prime_at_max << std::endl;
    Rcout<<"PDF components: Integral=" << termA <<", TermB=" << termB <<", Nu=" << nu_t << ", Beta_prime_term=" << beta_prime_term << ", Image_term=" << image_term << ", PDF=" << pdf << std::endl;
    return std::max(0.0, pdf);
}

double ou_fht_pdf_spline(double t, double lambda, double theta, double sigma, double z0, double b0, double binf, double tau, double p, int num_steps) {
    RD_Params pars = prepare_ou_params(t, lambda, theta, sigma, z0, b0, binf, tau, p);
    const double omega = pars.omega;
    if (t < FPM_EPSILON || lambda <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }

    if (std::abs(pars.b_scaled) < FPM_EPSILON) { // Analytical case for b_scaled = 0
        double sinh_t = std::sinh(pars.t_scaled);
        if (sinh_t < FPM_EPSILON) {
            return 0.0;
        } else {
            // TOD0: implement the sp_var case here
            double common_arg = -pars.zU_scaled * std::exp(-0.5 * pars.t_scaled) / std::sqrt(sinh_t);
            const double exp_term = std::exp(-0.5 * common_arg * common_arg + 0.5 * pars.t_scaled);
            const double g_scaled = (pars.zU_scaled * exp_term) / (std::sqrt(2.0 * M_PI) * std::pow(sinh_t, 1.5));
            return lambda * g_scaled;
        }
    }
	  const double et  = std::exp(pars.t_scaled);
	  const double e2t = std::exp(2.0 * pars.t_scaled);
    const double e2tm1 = e2t - 1.0;
    const double theta_max = et - 1.0;
    const double tau_max = 0.5 * (e2tm1);
    const int N = (num_steps < 2) ? 2 : num_steps;
    std::vector<double> theta_grid_block;
	
    pars.theta_max = theta_max;

    if (omega*(pars.zU_scaled - pars.b_scaled) < FPM_EPSILON) { 
         return 0.0;
    }
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
            return kernel_forward_tv_core(tt, ttp, p, beta_cache);
        });
        cache_ptr = &beta_cache;
    }

    ForcingFn g_fn;
    if (pars.fixed_b) {
        g_fn = pars.sp_var ? ForcingFn(g_term_forward_uniform)
                           : ForcingFn(g_term_forward);
    } else {
        g_fn = ForcingFn([cache_ptr](double theta_val, const RD_Params &p) {
                if (cache_ptr) {
                    return p.sp_var ? g_term_forward_tv_uniform(theta_val, p, *cache_ptr)
                                     : g_term_forward_tv(theta_val, p, *cache_ptr);
                }
                return p.sp_var ? g_term_forward_uniform(theta_val, p)
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
    const double beta_t = beta_grid.back();
    pars.beta_prime = beta_t; // constant-barrier default
    util::AkimaSpline beta_spline(t_grid, beta_grid);
    pars.beta_prime = beta_spline.derivative(t_grid.back());

    auto nu_f_vals = solve_nu_block_with_abel(theta_max, pars, N, theta_grid_block, kernel_fn, g_fn, abel_fn_theta);
    util::AkimaSpline nu_spline(theta_grid_block, nu_f_vals);
    
    const double nu_t = nu_spline.interpolate(theta_max);
    const double nu_prime_at_max = nu_spline.derivative(theta_max);
    // Equations here are re-arranged to match general form in Lipton &
    // Kaushansky (2020) Equation 10 They are mathematically equivalent to the
    // original form in section 4.4, but easier to relate to the general case
    // Term 1: The integral from Eq. (10)
    const double integral_sum = integrate_pdf_forward_spline(theta_max, nu_spline, pars, cache_ptr, num_steps);
    const double termA = (1/std::sqrt(M_PI))*integral_sum;
    // Term 2: -ν(t) / sqrt(2πt)
    const double termB = -nu_t / std::sqrt(2.0 * M_PI * tau_max);

    // The b'(t)ν(t) part
    const double beta_prime_term = omega*(-(pars.beta_prime/scale) * nu_t);
    const double image_term = omega* (0.5* averaged_image_term(tau_max, beta_t, pars));
    const double g_scaled = termA + termB + beta_prime_term + image_term;
    const double pdf = lambda * g_scaled * e2t;
    Rcout<<"nu_prime=" << nu_prime_at_max << std::endl;
    Rcout<<"PDF components: Integral=" << termA <<", TermB=" << termB <<", Nu=" << nu_t << ", Beta_prime_term=" << beta_prime_term << ", Image_term=" << image_term << ", PDF=" << pdf <<", beta= "<<beta_t<<", beta_p= "<<pars.beta_prime<< std::endl;
    return std::max(0.0, pdf);
}

// [[Rcpp::export]]
NumericVector ou_fht_pdf_vec(NumericVector t,
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
NumericVector ou_fht_pdf_vec_spline(NumericVector t,
                             double lambda, double theta, double sigma,
                             double z0, double b0, double binf, double tau, double pow, int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? ou_fht_pdf_spline(ti, lambda, theta, sigma, z0, b0, binf, tau, pow, num_steps) : NA_REAL;
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

NumericVector ou_fht_pdf_vec_fixed_zero_branch(NumericVector t,
                                               const RD_Params &pars) {
  const int n = t.size() + 1;
  NumericVector out(n);
  const double lambda = pars.lambda;
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
        const double z = pars.zU_scaled;
        out[i] = prefactor * z * std::exp(-decay_coeff * z * z);
        continue;
      }

      double z_lo = pars.zL_scaled;
      double z_hi = pars.zU_scaled;
      if (z_hi < z_lo) std::swap(z_lo, z_hi);
      const double span = z_hi - z_lo;

      if (span <= FPM_EPSILON) {
        const double z = 0.5 * (z_lo + z_hi);
        out[i] = prefactor * z * std::exp(-decay_coeff * z * z);
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
      out[i] = prefactor * averaged_z;
    }
    return out;
}

// [[Rcpp::export]]
NumericVector ou_fht_pdf_vec_grid(NumericVector t,
                                  double lambda, double theta, double sigma,
                                  double z0, double b0, double binf,
                                  double tau, double pow)
{
  const int n = t.size();
  NumericVector out(n, NA_REAL);

  // Find maximum finite t (evaluation horizon). If none, return NA.
  double t_max = 0.0;
  bool has_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i]) && t[i] > t_max) { t_max = t[i]; has_finite = true; }
  }
  if (!has_finite || t_max <= 0.0 || lambda <= 0.0 || sigma <= 0.0) {
    // Degenerate or invalid inputs
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) out[i] = 0.0; else out[i] = NA_REAL;
    }
    return out;
  }

  // Prepare parameters at t_max (used only to configure scaling and boundary behaviour)
  RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, pow);
  const double omega = pars.omega;
  // Choose solver steps in the native theta variable (the block solver is uniform in theta)
  const double et_max  = std::exp(pars.t_scaled);
  const double theta_max = et_max - 1.0;
  int N = calculate_num_steps(t_max, 0.005, 200);
  Rcout << "num_steps: " << N << std::endl;

  // Special analytical case for b_scaled ~ 0 (fixed barrier at zero in scaled space)
  if (std::abs(pars.b_scaled) < FPM_EPSILON && pars.fixed_b) {
    return ou_fht_pdf_vec_fixed_zero_branch(t, pars);
  }

  // Build uniform theta-grid up to theta_max
  std::vector<double> theta_grid_block;
  const double h = theta_max / N;
  // prefill time grids
  std::vector<double> theta_grid(N+1);
  for (int j = 0; j <= N; ++j) {
    theta_grid[j] = j * h;
  }

  // Prepare kernel, forcing, Abel seed, and optional boundary cache
  BoundaryDecayCache beta_cache;
  const BoundaryDecayCache* cache_ptr = nullptr;
  KernelFn kernel_fn;
  if (pars.fixed_b) {
    kernel_fn = KernelFn(kernel_forward);
  } else {
    beta_cache = make_boundary_decay_cache(theta_max, N, pars, beta_from_theta_raw);
    kernel_fn = KernelFn([&beta_cache](double tt, double ttp, const RD_Params& p) {
      return kernel_forward_tv_core(tt, ttp, p, beta_cache);
    });
    cache_ptr = &beta_cache;
  }

  ForcingFn g_fn;
  if (pars.fixed_b) {
    g_fn = pars.sp_var ? ForcingFn(g_term_forward_uniform)
                       : ForcingFn(g_term_forward);
  } else {
    g_fn = ForcingFn([cache_ptr](double theta_val, const RD_Params &p) {
      if (cache_ptr) {
        return p.sp_var ? g_term_forward_tv_uniform(theta_val, p, *cache_ptr)
                        : g_term_forward_tv(theta_val, p, *cache_ptr);
      }
      return p.sp_var ? g_term_forward_uniform(theta_val, p)
                      : g_term_forward(theta_val, p);
    });
  }

  const BoundaryDecayCache* abel_cache_theta = cache_ptr;
  AbelFn abel_fn_theta = AbelFn([abel_cache_theta](double theta_val, const RD_Params& p) {
    return abel_approx_nu_f(theta_val, p, abel_cache_theta);
  });

  // Precompute beta-grid across the mesh for later local derivatives
  std::vector<double> beta_grid(N+1), betap_grid(N+1);
  for (int j = 0; j <= N; ++j) {
    beta_grid[j]  = beta_from_theta(theta_grid[j], pars, cache_ptr);
  }
  util::AkimaSpline beta_spline(theta_grid, beta_grid);
  pars.beta_prime = beta_spline.derivative(theta_grid.back());
  auto nu_f_vals = solve_nu_block_with_abel(theta_max, pars, N, theta_grid_block, kernel_fn, g_fn, abel_fn_theta);
  util::AkimaSpline nu_spline(theta_grid_block, nu_f_vals);
  // theta_grid is now filled with the uniform grid used by the solver, with a spline fit to it

  for (int i=0;i<n;++i) {
    const double ti = t[i];
    if (!R_finite(ti) || ti <= 0.0) {
      out[i] = R_finite(ti) ? 0.0 : NA_REAL;
      continue;
    }
    RD_Params pars_i = prepare_ou_params(ti, lambda, theta, sigma, z0, b0, binf, tau, pow);
    const double et  = std::exp(pars_i.t_scaled);
    const double e2t = et * et;
    const double theta_eval = et - 1.0;
    const double tau_eval   = 0.5 * (e2t - 1.0);
    
    
    const double nu_t = nu_spline.interpolate(theta_eval);
    const double scale = 1+theta_eval;
    pars_i.beta_prime = beta_spline.derivative(theta_eval);
    const double beta_t = beta_spline.interpolate(theta_eval);
    // Derivative diagnostics at theta_eval: spline vs local quadratic fit
    double nu_prime_spline = nu_spline.derivative(theta_eval);
    double nu_prime_local = 0.0;
    {
      const int nTheta = static_cast<int>(theta_grid_block.size());
      if (nTheta >= 3) {
        int jpanel = static_cast<int>(theta_eval / h);
        if (jpanel < 1) jpanel = 1;
        if (jpanel > nTheta - 2) jpanel = nTheta - 2;
        const int c_index = jpanel + 1; // triple (jpanel-1, jpanel, jpanel+1)
        nu_prime_local = local_quad_derivative(theta_grid_block, nu_f_vals, c_index);
      }
    }
    // Equations here are re-arranged to match general form in Lipton &
    // Kaushansky (2020) Equation 10 They are mathematically equivalent to the
    // original form in section 4.4, but easier to relate to the general case
    // Term 1: The integral from Eq. (10)
    // Match single-t accuracy: scale endpoint u-panels with theta_eval/h,
    // but keep a reasonable floor and cap for performance.
    int M_int = static_cast<int>(theta_eval / h);
    if (M_int < 200) M_int = 200;
    if (M_int > N)   M_int = N;
    const double integral_sum = integrate_pdf_forward_spline(theta_eval, nu_spline, pars_i, cache_ptr, M_int);
    const double termA = (1/std::sqrt(M_PI))*integral_sum;
    // Term 2: -ν(t) / sqrt(2πt)
    const double termB = -nu_t / std::sqrt(2.0 * M_PI * tau_eval);

    // The b'(t)ν(t) part
    const double omega_i = pars_i.omega;
    const double beta_prime_term = omega_i*(-(pars_i.beta_prime/scale) * nu_t);
    const double image_term = omega_i* (0.5* averaged_image_term(tau_eval, beta_t, pars_i));
    const double g_scaled = termA + termB + beta_prime_term + image_term;
    const double pdf = lambda * g_scaled * e2t;
    Rcout<<"nu_prime_spline="<< nu_prime_spline << ", nu_prime_local=" << nu_prime_local << std::endl;
    Rcout<<"PDF components: Integral=" << termA <<", TermB=" << termB <<", Nu=" << nu_t << ", Beta_prime_term=" << beta_prime_term << ", Image_term=" << image_term << ", PDF=" << pdf <<", Theta= "<<theta_eval<<", beta= "<<beta_t<<", beta_p= "<<pars_i.beta_prime << std::endl;
    out[i] = pdf;
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
RD_Params pars = prepare_ou_params(t_max, lambda, theta, sigma, z0, b0, binf, tau, p);
const bool use_fixed = pars.fixed_b;
const double c = pars.c;

// Time grid (original time for boundary eval); step count
const int steps = std::max(1, (int)std::ceil(t_max / dt));
const double ds = lambda * dt;          // scaled time step
// Precompute standardized boundary trajectory B_k (k=0..steps)
std::vector<double> B_scaled(steps + 1);
if (use_fixed) {
    std::fill(B_scaled.begin(), B_scaled.end(), pars.b_scaled);
} else {
    for (int k = 0; k <= steps; ++k) {
        const double t_k = k * dt;
        const double b_t = evaluate_boundary_decay(t_k, pars);
        B_scaled[k] = c * (b_t - pars.theta);
    }
}

  const double e = std::exp(-ds);
  const double var_step = 0.5 * (1.0 - std::exp(-2.0 * ds));
  const double sd_step = std::sqrt(var_step);

  NumericVector X(n, NA_REAL);

  for (int j = 0; j < n; ++j) {
    // Sample start point: X0 ~ Uniform(0, z0); point start at 0 if z0 == 0
    const double X0 = (z0 > 0.0) ? R::runif(0.0, z0) : 0.0;
    double Z = c * (X0 - pars.theta);
    const double B0 = B_scaled[0];
    double s = 0.0;      // Time variable is now the scaled time 's'
    bool hit = false;

    // Edge case: already at the barrier at t=0
    if (std::abs(Z - B0) <= FPM_EPSILON) {
      X[j] = 0.0;
      continue;
    }

    const double dir = (Z < B0) ? 1.0 : -1.0;

    for (int k = 0; k < steps; ++k) {
      const double Z_prev = Z;
      
      // Exact transition for standardized process (theta_std = 0.0)
      Z = Z_prev * e + sd_step * R::rnorm(0.0, 1.0);
      const double s_next = s + ds;

      
      // Moving-boundary crossing detection via sign change of D = Z - B
      const double D_prev = Z_prev - B_scaled[k];
      const double D_next = Z - B_scaled[k + 1];
      
      if ((dir > 0.0 && D_prev < 0.0 && D_next >= 0.0) ||
          (dir < 0.0 && D_prev > 0.0 && D_next <= 0.0)) {
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

#endif // OU_HITTING_TIME_H






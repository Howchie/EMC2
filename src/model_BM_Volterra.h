#ifndef BM_HITTING_TIME_H
#define BM_HITTING_TIME_H

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

inline RD_Params prepare_bm_params(double t,
                                   double mu,
                                   double sigma,
                                   double z0,
                                   double b0,
                                   double binf,
                                   double tau,
                                   double pow,
                                   BoundaryDecayFn boundary_fn = BoundaryDecayFn(),
                                   std::vector<double> boundary_params = {}) {
  if (sigma <= 0.0) {
    stop("prepare_bm_params: sigma must be positive.");
  }
  RD_Params pars{};
  pars.b0 = b0;
  pars.binf = binf;
  pars.mu = mu;
  pars.sigma = sigma;
  pars.z0 = z0;
  pars.tau = tau;
  pars.pow = pow;
  pars.t_scaled = t;
  pars.zL_scaled = 0.0/sigma;
  pars.zU_scaled = z0/sigma;
  pars.omega = 1.0;
  if (b0>=FPM_EPSILON) { // Hitting from above
    pars.omega = -1.0;
  }
  pars.fixed_b = (std::abs(binf - b0) <= FPM_EPSILON) || tau <= 0.0 || pow <= 0.0;
  if (pars.fixed_b) {
    pars.tau = 0.0;
    pars.pow = 0.0;
  }
  pars.sp_var = (z0 > 0);
  pars.boundary_params = std::move(boundary_params);
  if (!boundary_fn) {
    boundary_fn = pars.fixed_b ? BoundaryDecayFn(fixed_boundary_decay) : BoundaryDecayFn(default_boundary_decay);
  }
  pars.boundary_fn = std::move(boundary_fn);

  return pars;
}

inline double physical_boundary(double t, const RD_Params& pars) {
  return evaluate_boundary_decay(t, pars);
}

inline double beta_from_time_raw(double t, const RD_Params& pars) {
  const double b_t = -pars.omega*physical_boundary(t, pars);
  const double shifted = b_t - pars.mu * t;
  return shifted/pars.sigma;
}

inline double beta_from_time(double t,
                             const RD_Params& pars,
                             const BoundaryDecayCache* cache = nullptr) {
  if (cache && !cache->empty()) {
    const double cached = cache->lookup(t);
    if (std::isfinite(cached)) {
      return cached;
    }
  }
  return beta_from_time_raw(t, pars);
}

inline double boundary_from_time(double t,
                                 const RD_Params& pars,
                                 const BoundaryDecayCache* cache = nullptr) {
  if (pars.fixed_b) {
    return physical_boundary(t, pars);
  }
  if (cache && !cache->empty()) {
    const double beta_t = cache->lookup(t);
    if (std::isfinite(beta_t)) {
      return (beta_t * pars.sigma) + pars.mu * t;
    }
  }
  return physical_boundary(t, pars);
}

inline double beta_prime_from_time(double t,
                                       const RD_Params& pars,
                                       const BoundaryDecayCache* cache = nullptr) {
  if (cache && cache->has_derivatives()) {
    const double beta_prime = cache->lookup_prime(t);
    if (std::isfinite(beta_prime)) {
      return beta_prime;
    }
  }
  const double step = (cache && cache->h > 0.0) ? cache->h : std::max(1e-6, 0.5 * std::max(1e-6, std::abs(t)));
  const double forward = t + step;
  const double backward = (t > step) ? (t - step) : 0.0;
  double beta_prime = std::numeric_limits<double>::quiet_NaN();
  if (backward <= FPM_EPSILON) {
    const double beta_t = beta_from_time_raw(t, pars);
    const double beta_forward = beta_from_time_raw(forward, pars);
    const double dt = forward - t;
    beta_prime = (dt > 0.0) ? (beta_forward - beta_t) / dt : 0.0;
  } else {
    const double beta_forward = beta_from_time_raw(forward, pars);
    const double beta_backward = beta_from_time_raw(backward, pars);
    const double dt = forward - backward;
    beta_prime = (dt > 0.0) ? (beta_forward - beta_backward) / dt : 0.0;
  }
  if (!std::isfinite(beta_prime)) {
    beta_prime = beta_prime_from_time(std::max(t, 1e-6), pars);
  }
  return beta_prime;
}

inline double boundary_prime_from_time(double t, const RD_Params& pars,
                                       const BoundaryDecayCache* cache = nullptr) {
  if (pars.fixed_b) {
	return 0.0;
  }
  const double beta_prime = beta_prime_from_time(t,pars,cache);
  return (beta_prime * pars.sigma) + pars.mu;
}


inline double kernel_bm(double t,
                        double s,
                        const RD_Params& pars,
                        const BoundaryDecayCache* cache = nullptr) {
  const double dt = t - s;
  // omega = sgn(z - b(0)) -> sgn(0 - beta(0))
  const double omega = pars.omega;
    if (dt <= FPM_EPSILON) {
        const double beta_prime_t = beta_prime_from_time(t, pars);
        return -omega * (beta_prime_t / std::sqrt(2.0 * M_PI));
    }

    const double beta_t = beta_from_time(t, pars, cache);
    const double beta_s = beta_from_time(s, pars, cache);
    
    // Psi(t, t') = b(t) - b(t') in the paper's notation
    const double psi = beta_t - beta_s;
    // Xi(t, t') = exp(-Psi^2 / (2 * (t - t')))
    const double xi_on_sqrt_dt = Gstar(dt, psi)/std::sqrt(dt); // G_star returns xi/sqrt(2π(t-s)); dividing by sqrt(t-s) gives xi/sqrt(2π)*(1/(t-s))
  

    // The full kernel from Eq. (11), noting that the sqrt(t-s) is handled by
    // the solver
    return -omega * psi * xi_on_sqrt_dt;
}

inline double kernel_bm_backward(double t,      // Current time in the solver grid (from 0 to T)
                                 double s,      // Past time in the solver grid (from 0 to t)
                                 const RD_Params& pars,
                                 const BoundaryDecayCache* cache = nullptr) {
    
    // T is the fixed final time for the CDF calculation
    const double t_max = pars.t_scaled;

    // The backward method evaluates the boundary at time T-t, T-s, etc.
    const double time_T_minus_t = std::max(0.0, t_max - t);
    const double time_T_minus_s = std::max(0.0, t_max - s);
    const double dt = t - s; // same as dt
    if (dt <= FPM_EPSILON) {
        return 0.0;
    }
    
    const double beta_T_minus_t = beta_from_time_raw(time_T_minus_t, pars);
    const double beta_T_minus_s = beta_from_time_raw(time_T_minus_s, pars);

    // Psi(T-t, T-s) = b(T-t) - b(T-s) in the paper's notation
    const double psi = beta_T_minus_t - beta_T_minus_s;

    // Xi(T-t, T-s) = exp(-Psi^2 / (2 * (t - s)))
    const double xi_on_sqrt_dt = Gstar(dt, psi) / std::sqrt(dt); // same reasoning: extra 1/sqrt(t-s) turns G_star into the (t-s)^{-1} kernel factor

    const double omega = pars.omega; 

    return -omega * psi * xi_on_sqrt_dt;
}


inline double forcing_bm_point(double t,
                               const RD_Params& pars,
                               const BoundaryDecayCache* cache = nullptr) {
    if (t <= FPM_EPSILON) {
        // Avoid division by zero at t=0; forcing is undefined but limit is 0.
        return 0.0; 
    }
    const double beta_t = beta_from_time(t, pars, cache);
    
    // This is the right-hand side of Eq. (11), with z=0 and b(t)=beta(t)
    return -Gstar(t,pars.zL_scaled-beta_t);
}

inline double forcing_bm_uniform(double t,
                                 const RD_Params& pars,
                                 const BoundaryDecayCache* cache = nullptr) {
  if (t <= FPM_EPSILON) {
    return 0.0;                       // limit is zero as t→0
  }
  const double z_hi = std::max(pars.z0, 0.0);   // assumed uniform on [0, z0]
  if (z_hi <= FPM_EPSILON) {
    return forcing_bm_point(t, pars, cache);    // degenerates to point start
  }

  const double beta_t = beta_from_time(t, pars, cache);
  const double inv_sigma = 1.0 / pars.sigma;
  const double lower = -beta_t;                          // z = 0
  const double upper = z_hi * inv_sigma - beta_t;         // z = z0
  const double scale = std::sqrt(t);

  // G* averaged over z ∈ [0, z0] → σ/z0 × [Φ(upper/√t) − Φ(lower/√t)]
  const double cdf_upper = normal_cdf(upper / scale);
  const double cdf_lower = normal_cdf(lower / scale);
  const double averaged = (pars.sigma / z_hi) * (cdf_upper - cdf_lower);

  return -averaged;
}

// Unified forcing selector
inline double forcing_bm(double t,
                               const RD_Params& pars,
                               const BoundaryDecayCache* cache = nullptr) {
  if (!pars.sp_var) {
    // Point start at z0 = S.zL (optionally: assert S.zL ≈ pars.z0)
    return forcing_bm_point(t, pars, cache);
  } else {
    return forcing_bm_uniform(t, pars, cache);
  }
}

inline double abel_approx_bm(double t,
                             const RD_Params& pars,
                             const BoundaryDecayCache* cache = nullptr) {
  // Use the forcing itself as the small-t seed; the Volterra integral is higher order in t.
  return forcing_bm(std::max(t, 1e-12), pars, cache);
}

inline double abel_approx_bm_backward(double t,
                                      const RD_Params& pars,
                                      const BoundaryDecayCache* cache = nullptr) {
  (void)cache;
  if (t <= FPM_EPSILON) {
    return 1.0;
  }
  const double drift = pars.mu / pars.sigma;
  const double sqrt_t = std::sqrt(std::max(t, 0.0));
  const double arg = drift * sqrt_t;
  const double prefactor = safe_exp(0.5 * drift * drift * t);
  const double value = 2.0 * prefactor * normal_cdf(arg);
  return std::max(0.0, value);
}

// A regularized integrand for the BM
inline double regularized_pdf_integrand_bm(
    double t_max, double t_prime,
    double nu_at_max, double nu_at_prime,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr)
{
    const double dt = t_max - t_prime;
    if (dt <= FPM_EPSILON) return 0.0;

    const double beta_max = beta_from_time(t_max, pars, cache);
    const double beta_prime = beta_from_time(t_prime, pars, cache);
    const double psi = beta_max - beta_prime;

    // Kernel's smooth part, H_smooth(t, t')
    const double psi_sq = psi * psi;
    const double r = psi_sq / dt;
    const double e = safe_exp(-0.5 * r);
    const double H_smooth = std::fma(-r, e, e);   // e*(1 - r) with one rounding

    const double denom = std::sqrt(2.0 * M_PI) * dt; // 1/sqrt(dt) (i.e., u) cancels with the -2u in the jacobian
    if (denom < FPM_EPSILON) return 0.0;

    return 2.0 * (H_smooth * nu_at_prime - nu_at_max) / denom;
}

// 2. The high-order integrator for the BM PDF, mimicking your OU version
double integrate_pdf_forward_bm_u(
    double t_max,
    const std::vector<double>& t_grid,
    const std::vector<double>& nu_vals,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr)
{
    const int N = t_grid.size() - 1;
    if (N < 1) return 0.0;

    const double nu_at_max = nu_vals.back();
	
	std::vector<double> u_grid(N + 1);

    // The u-transform: u = sqrt(t_max - t_prime) => dt' = -2u du
    // Integral from t'=0..t_max becomes integral u=sqrt(t_max)..0
	
	for (int i = 0; i <= N; ++i) {
        u_grid[i] = sqrt_pos(t_max - t_grid[i]);
    }
	
	// Evaluate f(u_N) = f(0) using the derivative limit
    const double mid = 0.5 * (t_grid[N-1] + t_max);
	const RightQuad rq = right_end_quadratic(t_grid, nu_vals, mid);
    const double nu_prime = rq.nu_prime;
    const double beta_prime_t = pars.beta_prime;
    const double beta_prime_sq = beta_prime_t * beta_prime_t;
    const double f_right_limit = -2.0 * (nu_prime + 1.5 * beta_prime_sq * nu_at_max) / std::sqrt(2.0 * M_PI);
    
    double integral_sum = 0.0;
    
    // We integrate over u from u_max=sqrt(t_max) down to u_min=0
    for (int i = N - 1; i >= 0; --i) {
        // Endpoints of the interval in the original t-grid
        const double t_L = t_grid[i];
        const double t_R = t_grid[i+1];

        // Corresponding endpoints in the u-grid
        const double u_L = u_grid[i];
        const double u_R = u_grid[i+1];
        const double panel_h_u = u_L - u_R;

        // Midpoint in u-space
        const double u_M = 0.5 * (u_L + u_R);
        // Corresponding midpoint in t-space
        const double t_M = t_max - u_M * u_M;

        // Get nu values at the three points for Simpson's rule
        const double nu_L = nu_vals[i];
        const double nu_R = nu_vals[i+1];
        // Interpolate to find nu at the midpoint
        const double nu_M = quad_interp(t_grid, nu_vals, i, t_M);

        // Evaluate the regularized integrand f(t') at the three points (regularised integrand handles the Jacobian)
        const double f_L = regularized_pdf_integrand_bm(t_max, t_L, nu_at_max, nu_L, pars, cache);
        const double f_R = regularized_pdf_integrand_bm(t_max, t_R, nu_at_max, nu_R, pars, cache);
        const double f_M = regularized_pdf_integrand_bm(t_max, t_M, nu_at_max, nu_M, pars, cache);
		
		
		if (i == N - 1) {
			const double a = f_right_limit; // constant term                       // previous u
			const double f_u2 = regularized_pdf_integrand_bm(
				t_max, t_grid[N - 1], nu_at_max, nu_vals[N - 1], pars, cache);
			const double b = (f_L - a) / (u_L * u_L);                  // quadratic term
			integral_sum += a * u_L + (b / 3.0) * u_L * u_L * u_L;        // analytic tail
			continue; // skip Simpson for last panel
		}


        // Simpson's rule for this panel in u-space, including the Jacobian term (2u)
        integral_sum += (panel_h_u / 6.0) * ( f_L + 4.0 * f_M + f_R );
    }

    return integral_sum;
}

// 2b. Spline-based variant: identical numerics on-grid, but integrates using Akima spline
//     and the same regularized integrand as integrate_pdf_forward_bm_u. This mirrors
//     OU's integrate_pdf_forward_spline and allows per-time adaptive integration.
inline double integrate_pdf_forward_bm_spline(
    double t_k,
    const util::AkimaSpline& nu_spline,
    const RD_Params& pars,
    const BoundaryDecayCache* cache = nullptr,
    int M = 100,
    double nu_prime_at_k_override = std::numeric_limits<double>::quiet_NaN())
{
    if (t_k <= FPM_EPSILON) return 0.0;

    // Build a temporary uniform partition in t' and corresponding u-grid
    std::vector<double> t_grid(M + 1);
    std::vector<double> u_grid(M + 1);
    const double h_t = t_k / M;
    for (int i = 0; i <= M; ++i) {
        const double tp = i * h_t;
        t_grid[i] = tp;
        u_grid[i] = sqrt_pos(t_k - tp);
    }

    // Endpoint values via spline
    const double nu_at_k = nu_spline.interpolate(t_k);
    const double nu_prime_at_k = std::isfinite(nu_prime_at_k_override)
                                     ? nu_prime_at_k_override
                                     : nu_spline.derivative(t_k);

    // Right-end limit of the integrand in u-space for BM (matches integrate_pdf_forward_bm_u)
    const double beta_prime_t = pars.beta_prime;
    const double beta_prime_sq = beta_prime_t * beta_prime_t;
    const double f_right_limit = -2.0 * (nu_prime_at_k + 1.5 * beta_prime_sq * nu_at_k) / std::sqrt(2.0 * M_PI);

    double integral_sum = 0.0;
    for (int i = M - 1; i >= 0; --i) {
        const double tL = t_grid[i];
        const double tR = t_grid[i + 1];
        const double uL = u_grid[i];
        const double uR = u_grid[i + 1];
        const double h_u = uL - uR; // panel size in u
        if (h_u == 0.0) continue;

        const double uM = 0.5 * (uL + uR);
        const double tM = t_k - uM * uM;

        // nu values from spline
        const double nuL = nu_spline.interpolate(tL);
        const double nuR = nu_spline.interpolate(tR);
        const double nuM = nu_spline.interpolate(tM);

        // Regularized integrand at the three points
        const double fL = regularized_pdf_integrand_bm(t_k, tL, nu_at_k, nuL, pars, cache);
        const double fR = (i == M - 1)
                            ? f_right_limit
                            : regularized_pdf_integrand_bm(t_k, tR, nu_at_k, nuR, pars, cache);
        const double fM = regularized_pdf_integrand_bm(t_k, tM, nu_at_k, nuM, pars, cache);

        // Simpson on this u-panel
        integral_sum += (h_u / 6.0) * (fL + 4.0 * fM + fR);
    }

    return integral_sum;
}

double calculate_cdf_from_nu_backward(const std::vector<double>& t_grid,
                                      const std::vector<double>& nu_vec,
                                      const RD_Params& pars) {

  const double t = pars.t_scaled;
  const double omega = pars.omega;
  double integral_sum = 0.0;
  const int N = static_cast<int>(t_grid.size()) - 1;
    // Product integration over the solution
    for (size_t i = 0; i < N; ++i) {
        const double s1 = t_grid[i];
        const double s2 = t_grid[i+1];
        const double ds = s2 - s1;
        if (ds < FPM_EPSILON) continue;
        auto integrand = [&](double s, double nu_s) {
            const double T_minus_s = t - s;
            if (T_minus_s <= FPM_EPSILON) return 0.0;

            const double beta_T_minus_s = beta_from_time_raw(T_minus_s, pars);

            // The backward CDF integrand (Eq. 12) is exactly twice the forward
            // PDF's image term (Eq. 10). We can reuse the same function.
            const double image_term = averaged_image_term(T_minus_s, beta_T_minus_s, pars);

            return (image_term*std::sqrt(T_minus_s) * nu_s);
        };

        const double H_a = integrand(s1, nu_vec[i]);
		    double H_b = 0.0; 
 
        // Final interval b==v and Hb=0
        double W_a, W_b;
        if (i == N - 1) {  // Check if we are in the final interval where b ~ t
			    // The smooth part H_b at t'=t is zero due to the exponential term.
			    weights_panel_12(t, s1, s2, W_a, W_b);
        } else {
           // Use the standard product trapezoidal formula for all other intervals.
            H_b = integrand(s2, nu_vec[i+1]);
			      weights_panel_12(t, s1, s2, W_a, W_b);
        }

        integral_sum += W_a * H_a + W_b * H_b;
    }
    
    return std::max(0.0, std::min(1.0, omega * integral_sum));
}

double bm_fht_pdf(double t,
                         double mu,
                         double sigma,
                         double z0,
                         double b0,
                         double binf,
                         double tau,
                         double pow,
                         int num_steps,
                         BoundaryDecayFn boundary_fn = BoundaryDecayFn(),
                         std::vector<double> boundary_params = {}) {
  if (t <= 0.0) {
    return 0.0;
  }
  const int N = (num_steps < 2) ? 2 : num_steps;
  std::vector<double> t_grid(N + 1);
  const double h_t = t / N;
  for (int j = 0; j <= N; ++j) t_grid[j] = j * h_t;

  RD_Params pars = prepare_bm_params(t, mu, sigma, z0, b0, binf, tau, pow,
                                             std::move(boundary_fn),
                                             std::move(boundary_params));
  BoundaryDecayCache beta_cache;
  const bool use_cache = !pars.fixed_b || std::abs(pars.mu) > FPM_EPSILON;
  const BoundaryDecayCache* cache_ptr = nullptr;
  if (use_cache) {
    beta_cache = make_boundary_decay_cache(t, N, pars, beta_from_time_raw);
    cache_ptr = &beta_cache;
  }

  KernelFn kernel_fn = KernelFn([cache_ptr](double tt, double ss, const RD_Params& p) {
    return kernel_bm(tt, ss, p, cache_ptr);
  });

  ForcingFn forcing_fn = ForcingFn([cache_ptr](double tt, const RD_Params& p) {
    return forcing_bm(tt, p, cache_ptr);
  });

  AbelFn abel_fn = AbelFn([cache_ptr](double tt, const RD_Params &p) {
    return abel_approx_bm(tt, p, cache_ptr);
  });
  
  std::vector<double> beta_grid(N + 1);
  for (int j = 0; j <= N; ++j) beta_grid[j] = beta_from_time(t_grid[j], pars, cache_ptr);
  const double beta_t = beta_grid.back();
  pars.beta_prime = -pars.mu/pars.sigma; // constant-barrier default
  //if (!pars.fixed_b) {
      const double mid = 0.5 * (t_grid[N-1] + t_grid.back()); // theta_grid is N+1 in size 
      const RightQuad rq = right_end_quadratic(t_grid, beta_grid,mid); 
      pars.beta_prime = rq.nu_prime;
  //}
  auto nu_vec = solve_nu_block_with_abel(t, pars, N, t_grid, kernel_fn,
                                          forcing_fn, abel_fn);

  // Spline for nu for exact on-grid evaluation and stable derivatives
  util::AkimaSpline nu_spline(t_grid, nu_vec);

	  const double nu_t = nu_spline.interpolate(t);
    const double omega = pars.omega;

    // Term 1: The integral from Eq. (10) using spline-based integrator
    const double termA = 0.5 * integrate_pdf_forward_bm_spline(t, nu_spline, pars, cache_ptr, N);
    
    // Term 2: -ν(t) / sqrt(2πt)
    const double termB = -nu_t / std::sqrt(2.0 * M_PI * t);

    // Term 3: The standard inverse-Gaussian-like part
    // Image term from Eq. (10): ((z - b(t)) exp(-Psi_forward^2 /(2 t)))/(2
    // sqrt{2 pi t^3}).
    // When z0 > 0 we need to average this over the uniform start distribution by analytically integrating
    // the image term over z ∈ [0, z0].
    const double image_term = 0.5* averaged_image_term(t, beta_t, pars);
    const double termC = omega * (-pars.beta_prime * nu_t + image_term);


    const double density = termA + termB + termC;

  return density;
}

// Forward declarations for chunked BM grid solver
NumericVector bm_fht_pdf_vec_grid_chunked(NumericVector t,
                                          double mu, double sigma, double z0,
                                          double b0, double binf, double tau, double pow,
                                          double steps_fineness, double min_steps,
                                          double chunk_ratio, double chunk_base_panels,
                                          double chunk_max);

// [[Rcpp::export]]
NumericVector bm_fht_pdf_vec_grid(NumericVector t,
                                  double mu, double sigma, double z0,
                                  double b0, double binf, double tau, double pow,
                                  double steps_fineness, double min_steps)
{
  const int n = t.size();
  // Determine horizon and earliest positive time
  double t_max = 0.0;
  double t_min_pos = std::numeric_limits<double>::infinity();
  bool has_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i])) {
      if (t[i] > t_max) { t_max = t[i]; has_finite = true; }
      if (t[i] > 0.0 && t[i] < t_min_pos) t_min_pos = t[i];
    }
  }
  if (!has_finite || t_max <= 0.0 || sigma <= 0.0) {
    NumericVector out(n, NA_REAL);
    for (int i = 0; i < n; ++i) out[i] = (R_finite(t[i]) && t[i] >= 0.0) ? 0.0 : NA_REAL;
    return out;
  }

  int N = calculate_num_steps(t_max, steps_fineness, (int)min_steps);
  int P = 0;
  if (R_finite(t_min_pos) && t_min_pos > 0.0) {
    P = calculate_num_steps(t_min_pos, steps_fineness, (int)min_steps);
    if (P % 2 != 0) ++P;
    if (P > N - 2) P = std::max(2, N - 2);
    if ((N - P) % 2 != 0) { if (P > 2) --P; else ++N; }
  }
  const double ratio = 1.0; // zero compression path
  const double base_panels = static_cast<double>(std::max(2, N - P));
  const double max_chunks = 2.0;
  return bm_fht_pdf_vec_grid_chunked(t, mu, sigma, z0, b0, binf, tau, pow,
                                     steps_fineness, min_steps,
                                     ratio, base_panels, max_chunks);
}

// Variant of the grid-based solver using the chunked Volterra kernel
// acceleration, mirroring OU. Keeps legacy functions available.
// [[Rcpp::export]]
NumericVector bm_fht_pdf_vec_grid_chunked(NumericVector t,
                                          double mu, double sigma, double z0,
                                          double b0, double binf, double tau, double pow,
                                          double steps_fineness, double min_steps,
                                          double chunk_ratio, double chunk_base_panels,
                                          double chunk_max)
{
  const int n = t.size();
  NumericVector out(n, NA_REAL);

  double t_max = 0.0;
  bool has_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(t[i]) && t[i] > t_max) { t_max = t[i]; has_finite = true; }
  }
  if (!has_finite || t_max <= 0.0 || sigma <= 0.0) {
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) out[i] = 0.0; else out[i] = NA_REAL;
    }
    return out;
  }

  // Earliest positive time for prelude panels
  double t_min_pos = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (R_finite(ti) && ti > 0.0 && ti < t_min_pos) t_min_pos = ti;
  }

  RD_Params pars = prepare_bm_params(t_max, mu, sigma, z0, b0, binf, tau, pow);
  const bool use_cache = !pars.fixed_b || std::abs(mu) > FPM_EPSILON;
  BoundaryDecayCache beta_cache;
  BoundaryDecayCache* cache_ptr = nullptr;
  int N = calculate_num_steps(t_max, steps_fineness, (int)min_steps);
  if (use_cache) {
    beta_cache = make_boundary_decay_cache(t_max, N, pars, beta_from_time_raw);
    cache_ptr = &beta_cache;
  }

  KernelFn kernel_fn = KernelFn([cache_ptr](double tt, double ss, const RD_Params& p) {
    return kernel_bm(tt, ss, p, cache_ptr);
  });

  ForcingFn forcing_fn = ForcingFn([cache_ptr](double tt, const RD_Params& p) {
    return forcing_bm(tt, p, cache_ptr);
  });

  AbelFn abel_fn = AbelFn([cache_ptr](double tt, const RD_Params &p) {
    return abel_approx_bm(tt, p, cache_ptr);
  });

  // Chunking configuration (reuse theta helpers for generic t-grid)
  ChunkingOptions chunk_opts;
  chunk_opts.ratio = chunk_ratio;
  chunk_opts.base_panels = static_cast<int>(std::round(chunk_base_panels));
  chunk_opts.max_chunks = chunk_max;
  if (R_finite(t_min_pos) && t_min_pos > 0.0) {
    chunk_opts.theta_min = t_min_pos; // same semantics: position in the domain
    int P = calculate_num_steps(t_min_pos, steps_fineness, (int)min_steps);
    if (P % 2 != 0) ++P;
    chunk_opts.pre_min_panels = std::max(2, P);
  }

  std::vector<double> t_grid_block;
  auto t_grid_block_upfront = build_chunked_theta_grid(t_max, N, chunk_opts);
  std::vector<double> beta_grid(t_grid_block_upfront.size());
  for (size_t j = 0; j < beta_grid.size(); ++j) {
    beta_grid[j] = beta_from_time(t_grid_block_upfront[j], pars, cache_ptr);
  }
  util::AkimaSpline beta_spline(t_grid_block_upfront, beta_grid);
  pars.beta_prime = beta_spline.derivative(t_max);

  auto nu_vals = solve_nu_block_with_abel_chunked(t_max, pars, N, t_grid_block,
                                                  kernel_fn, forcing_fn, abel_fn,
                                                  chunk_opts);
  util::AkimaSpline nu_spline(t_grid_block, nu_vals);

  const int total_nodes = static_cast<int>(t_grid_block.size());

  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti)) { out[i] = NA_REAL; continue; }
    if (ti <= 0.0) { out[i] = 0.0; continue; }

    RD_Params pars_i = prepare_bm_params(ti, mu, sigma, z0, b0, binf, tau, pow);
    const double nu_ti = nu_spline.interpolate(ti);
    pars_i.beta_prime = beta_spline.derivative(ti);
    const double beta_ti = beta_spline.interpolate(ti);

    // For derivative override near the endpoint; compute a local quadratic derivative
    double deriv_override = std::numeric_limits<double>::quiet_NaN();
    if (total_nodes >= 3) {
      const double last_t = t_grid_block.back();
      if (std::abs(ti - last_t) <= std::max(1e-12, 1e-12 * std::abs(last_t))) {
        // Use right-end quadratic on the solved grid for stability
        const RightQuad rq = right_end_quadratic(t_grid_block, nu_vals, 0.5 * (t_grid_block[total_nodes-2] + last_t));
        deriv_override = rq.nu_prime;
      }
    }

    const double integral_sum = integrate_pdf_forward_bm_spline(ti, nu_spline, pars_i, cache_ptr, N, deriv_override);
    const double termA = 0.5 * integral_sum;
    const double termB = -nu_ti / std::sqrt(2.0 * M_PI * ti);
    const double omega_i = pars_i.omega;
    const double termC = omega_i * (-pars_i.beta_prime * nu_ti + 0.5 * averaged_image_term(ti, beta_ti, pars_i));
    const double density = termA + termB + termC;
    out[i] = std::max(0.0, density);
  }

  return out;
}

inline double bm_fht_cdf(double t,
                         double mu,
                         double sigma,
                         double z0,
                         double b0,
                         double binf,
                         double tau,
                         double pow,
                         int num_steps,
                         BoundaryDecayFn boundary_fn = BoundaryDecayFn(),
                         std::vector<double> boundary_params = {}) {
  if (t <= 0.0) {
    return 0.0;
  }
    const int N = (num_steps < 2) ? 2 : num_steps;
  std::vector<double> t_grid(N + 1);
  const double h_t = t / N;
  for (int j = 0; j <= N; ++j) t_grid[j] = j * h_t;

  RD_Params pars = prepare_bm_params(t, mu, sigma, z0, b0, binf, tau, pow,
                                             std::move(boundary_fn),
                                             std::move(boundary_params));
  BoundaryDecayCache beta_cache;
  const bool use_cache = !pars.fixed_b || std::abs(pars.mu) > FPM_EPSILON;
  const BoundaryDecayCache* cache_ptr = nullptr;
  if (use_cache) {
    beta_cache = make_boundary_decay_cache(t, N, pars, beta_from_time_raw);
    cache_ptr = &beta_cache;
  }

  KernelFn kernel_fn = KernelFn([cache_ptr](double tt, double ss, const RD_Params& p) {
    return kernel_bm_backward(tt, ss, p, cache_ptr);
  });

  ForcingFn forcing_fn = ForcingFn([](double, const RD_Params&) {
    return 1.0;
  });

  AbelFn abel_fn = AbelFn([cache_ptr](double tt, const RD_Params &p) {
    return abel_approx_bm_backward(tt, p, cache_ptr);
  });

  auto nu_vals = solve_nu_block_with_abel(t, pars, N, t_grid, kernel_fn,
                                          forcing_fn, abel_fn);
  const double cdf_val = calculate_cdf_from_nu_backward(t_grid, nu_vals, pars);
  return std::max(0.0, std::min(1.0, cdf_val));
}

// [[Rcpp::export]]
NumericVector bm_fht_cdf_vec_grid_chunked(NumericVector t,
                                          double mu, double sigma, double z0,
                                          double b0, double binf, double tau, double pow,
                                          double steps_fineness = 0.005,
                                          double min_steps = 300,
                                          double chunk_ratio = 1.5,
                                          double chunk_base_panels = 500,
                                          double chunk_max = 12,
                                          double rt_resolution = NA_REAL)
{
  const int n = t.size();
  NumericVector out(n, NA_REAL);

  // Determine horizon for building a uniform grid for trapezoidal accumulation
  double t_max = 0.0;
  for (int i = 0; i < n; ++i) if (R_finite(t[i]) && t[i] > t_max) t_max = t[i];
  if (!(sigma > 0.0) || t_max <= 0.0) {
    for (int i = 0; i < n; ++i) {
      if (R_finite(t[i]) && t[i] >= 0.0) out[i] = 0.0; else out[i] = NA_REAL;
    }
    return out;
  }

  // Uniform evaluation grid in t for CDF integration
  int N = calculate_num_steps(t_max, steps_fineness, (int)min_steps);
  std::vector<double> t_grid(N + 1);
  const double h = t_max / N;
  for (int j = 0; j <= N; ++j) t_grid[j] = j * h;

  // Evaluate PDF on this grid using the chunked forward solver (amortizes kernel solve)
  NumericVector t_grid_nv(t_grid.begin(), t_grid.end());
  NumericVector pdf_on_grid = bm_fht_pdf_vec_grid_chunked(t_grid_nv, mu, sigma, z0, b0, binf, tau, pow,
                                                          steps_fineness, min_steps,
                                                          chunk_ratio, chunk_base_panels, chunk_max);
  const size_t m = static_cast<size_t>(pdf_on_grid.size());
  std::vector<double> pdf_grid(m, 0.0);
  for (size_t i = 0; i < m; ++i) {
    const double val = pdf_on_grid[static_cast<R_xlen_t>(i)];
    pdf_grid[i] = (R_finite(val) && val > 0.0) ? val : 0.0;
  }

  // Trapezoidal accumulation to get CDF on the grid
  std::vector<double> cdf_grid(m, 0.0);
  for (size_t i = 1; i < m; ++i) {
    const double dt_ = t_grid[i] - t_grid[i - 1];
    if (dt_ > 0.0) {
      const double avg = 0.5 * (pdf_grid[i] + pdf_grid[i - 1]);
      cdf_grid[i] = cdf_grid[i - 1] + dt_ * avg;
    } else {
      cdf_grid[i] = cdf_grid[i - 1];
    }
    if (cdf_grid[i] < cdf_grid[i - 1]) cdf_grid[i] = cdf_grid[i - 1];
    if (cdf_grid[i] > 1.0) cdf_grid[i] = 1.0;
  }

  // Build Akima spline on either uniform grid or requested times
  std::vector<double> spline_x;
  std::vector<double> spline_y;
  if (R_finite(rt_resolution) && rt_resolution > 0.0) {
    spline_x = t_grid;
    spline_y = cdf_grid;
  } else {
    spline_x.reserve(static_cast<size_t>(n) + 1);
    for (int i = 0; i < n; ++i) {
      const double ti = t[i];
      if (R_finite(ti) && ti >= 0.0) spline_x.push_back(ti);
    }
    std::sort(spline_x.begin(), spline_x.end());
    spline_x.erase(std::unique(spline_x.begin(), spline_x.end(), [](double a, double b){ return std::abs(a-b) <= 1e-15; }), spline_x.end());
    spline_y.resize(spline_x.size());
    for (size_t i = 0; i < spline_x.size(); ++i) {
      spline_y[i] = linear_interp(t_grid, cdf_grid, spline_x[i]);
    }
  }

  bool use_akima = (spline_x.size() >= 5);
  std::unique_ptr<util::AkimaSpline> cdf_spline;
  if (use_akima) {
    try { cdf_spline = std::unique_ptr<util::AkimaSpline>(new util::AkimaSpline(spline_x, spline_y)); }
    catch (...) { use_akima = false; }
  }
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (!R_finite(ti)) { out[i] = NA_REAL; continue; }
    if (ti <= 0.0) { out[i] = 0.0; continue; }
    out[i] = use_akima ? cdf_spline->interpolate(ti) : linear_interp(t_grid, cdf_grid, ti);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector bm_fht_cdf_vec_grid(NumericVector t,
                                  double mu, double sigma, double z0,
                                  double b0, double binf, double tau, double pow,
                                  double steps_fineness, double min_steps,
                                  double rt_resolution = NA_REAL)
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
  if (!has_finite || t_max <= 0.0 || !(sigma > 0.0)) {
    NumericVector out(n, NA_REAL);
    for (int i = 0; i < n; ++i) out[i] = (R_finite(t[i]) && t[i] >= 0.0) ? 0.0 : NA_REAL;
    return out;
  }

  int N = calculate_num_steps(t_max, steps_fineness, (int)min_steps);
  int P = 0;
  if (R_finite(t_min_pos) && t_min_pos > 0.0) {
    P = calculate_num_steps(t_min_pos, steps_fineness, (int)min_steps);
    if (P % 2 != 0) ++P;
    if (P > N - 2) P = std::max(2, N - 2);
    if ((N - P) % 2 != 0) { if (P > 2) --P; else ++N; }
  }
  const double ratio = 1.0;
  const double base_panels = static_cast<double>(std::max(2, N - P));
  const double max_chunks = 2.0;

  return bm_fht_cdf_vec_grid_chunked(t, mu, sigma, z0, b0, binf, tau, pow,
                                     steps_fineness, min_steps,
                                     ratio, base_panels, max_chunks,
                                     rt_resolution);
}

// [[Rcpp::export]]
NumericVector bm_fht_pdf_vec(NumericVector t,
                             double mu,
                             double sigma,
                             double z0,
                             double b0,
                             double binf,
                             double tau,
                             double pow,
                             int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? bm_fht_pdf(ti, mu, sigma, z0, b0, binf, tau, pow, num_steps) : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector bm_fht_cdf_vec(NumericVector t,
                             double mu,
                             double sigma,
                             double z0,
                             double b0,
                             double binf,
                             double tau,
                             double pow,
                             int num_steps = 200) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    out[i] = R_finite(ti) ? bm_fht_cdf(ti, mu, sigma, z0, b0, binf, tau, pow, num_steps) : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector simulate_bm_hit_times(
    int n,
    double mu,    // drift
    double sigma, // diffusion SD (not variance)
    double z0, // upper bound of Uniform(0, z0) for start; if 0 => point start at 0
    double b0,   // boundary at t = 0
    double binf, // asymptote for decaying boundary; if == b0 use fixed boundary
    double tau = 1.0, // decay scale
    double pow = 1.0,   // decay shape
    double dt = 1e-3, double t_max = 10.0) {
  BoundaryDecayFn boundary_fn = BoundaryDecayFn();
  std::vector<double> boundary_params = {};
  if (n <= 0) stop("n must be positive.");
  if (!R_finite(mu) || !R_finite(sigma) || !R_finite(z0) || !R_finite(b0) || !R_finite(binf))
    stop("Non-finite parameter.");
  if (sigma <= 0.0) stop("sigma must be > 0 and is the diffusion SD.");
  if (dt <= 0.0)    stop("dt must be > 0.");
  if (t_max <= 0.0) stop("t_max must be > 0.");
  if (z0 < 0.0)     stop("z0 must be >= 0.");
// Determine if we use fixed boundary or time-varying boundary

  const int steps = std::max(1, (int)std::ceil(t_max / dt));
  const double sdt = std::sqrt(dt);
  const bool use_fixed = std::abs(b0 - binf) < 1e-12;

  // Precompute boundary trajectory B_k at t_k = k*dt
  // Precompute standardized boundary trajectory B_k (k=0..steps)
  RD_Params pars = prepare_bm_params(t_max, mu, sigma, z0, b0, binf, tau, pow);
  std::vector<double> B(steps + 1);
  if (use_fixed) {
      std::fill(B.begin(), B.end(), b0);
  } else {
      for (int k = 0; k <= steps; ++k) {
          const double t_k = k * dt;
          const double b_t = evaluate_boundary_decay(t_k, pars);
          B[k] = b_t;
      }
  }

  // We hit from below: require start < boundary at t=0 (or we hit at 0)
  Rcpp::NumericVector T(n, NA_REAL);

  for (int j = 0; j < n; ++j) {
    // Sample start point: X0 ~ Uniform(0, z0); point start at 0 if z0 == 0
    const double X0 = (z0 > 0.0) ? R::runif(0.0, z0) : 0.0;

    // If already at/above boundary at t=0, first-passage time is 0
    if (X0 >= B[0]) {
      T[j] = 0.0;
      continue;
    }

    double X = X0;
    double t = 0.0;
    bool hit = false;

    // Euler–Maruyama (exact for BM with drift on uniform grid increments)
    for (int k = 0; k < steps; ++k) {
      const double X_prev = X;
      const double B_prev = B[k];

      // advance one step
      X = X_prev + mu * dt + sigma * sdt * R::rnorm(0.0, 1.0);
      const double B_next = B[k + 1];

      // Detect crossing from below: D = X - B, look for D_prev < 0 and D_next >= 0
      const double D_prev = X_prev - B_prev;
      const double D_next = X - B_next;

      if (D_prev < 0.0 && D_next >= 0.0) {
        // Linear interpolation within step:
        // Solve for alpha in X_prev + alpha*(X - X_prev) = B_prev + alpha*(B_next - B_prev)
        const double num = (B_prev - X_prev);
        const double den = ( (X - X_prev) - (B_next - B_prev) );
        double alpha = 0.5; // fallback
        if (std::abs(den) > 0.0) {
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
      T[j] = NA_REAL; // no hit by t_max
    }
  }

  return T;
}

// [[Rcpp::export]]
Rcpp::NumericVector simulate_bm_hit_times_bb(
    int n,
    double mu,
    double sigma,
    double z0,
    double b0,
    double binf,
    double tau = 1.0,
    double pow = 1.0,
    double dt = 1e-3,
    double t_max = 10.0)
{
  if (n <= 0) stop("simulate_bm_hit_times_bb: n must be positive.");
  if (!R_finite(mu) || !R_finite(sigma) || !R_finite(z0) || !R_finite(b0) || !R_finite(binf) ||
      !R_finite(tau) || !R_finite(pow) || !R_finite(dt) || !R_finite(t_max)) {
    stop("simulate_bm_hit_times_bb: parameters must be finite.");
  }
  if (sigma <= 0.0 || dt <= 0.0 || t_max <= 0.0) {
    stop("simulate_bm_hit_times_bb: invalid parameters.");
  }
  if (z0 < 0.0) stop("simulate_bm_hit_times_bb: z0 must be >= 0.");

  RD_Params pars = prepare_bm_params(t_max, mu, sigma, z0, b0, binf, tau, pow);
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

  const double sdt = std::sqrt(dt);
  const double var_step = (sigma * sigma) * dt;

  NumericVector T(n, NA_REAL);
  for (int j = 0; j < n; ++j) {
    const double X0 = (z0 > 0.0) ? R::runif(0.0, z0) : 0.0;
    double X = X0;
    double t = 0.0;
    bool hit = false;

    if (std::abs(X - B[0]) <= FPM_EPSILON) { T[j] = 0.0; continue; }

    for (int k = 0; k < steps; ++k) {
      const double X_prev = X;
      const double B_prev = B[k];
      X = X_prev + mu * dt + sigma * sdt * R::rnorm(0.0, 1.0);
      const double B_next = B[k + 1];

      const double D_prev = X_prev - B_prev;
      const double D_next = X - B_next;

      // Straddle detection: linear interpolation
      if ((D_prev <= 0.0 && D_next >= 0.0) || (D_prev >= 0.0 && D_next <= 0.0)) {
        const double denom = (D_prev - D_next);
        double alpha = (std::abs(denom) > FPM_EPSILON) ? (D_prev / denom) : 0.5;
        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;
        T[j] = t + alpha * dt;
        hit = true;
        break;
      }

      // Brownian bridge correction when both endpoints on same side
      if (D_prev * D_next > 0.0) {
        const double a = std::abs(D_prev);
        const double b = std::abs(D_next);
        const double s2 = std::max(var_step, FPM_EPSILON);
        const double p_cross = std::exp(-2.0 * (a * b) / s2);
        if (R::runif(0.0, 1.0) < p_cross) {
          double alpha = R::rbeta(0.5, 0.5);
          if (!R_finite(alpha)) alpha = 0.5;
          if (alpha < 0.0) alpha = 0.0;
          if (alpha > 1.0) alpha = 1.0;
          T[j] = t + alpha * dt;
          hit = true;
          break;
        }
      }

      t += dt;
    }

    if (!hit) T[j] = NA_REAL;
  }

  return T;
}

#endif // BM_HITTING_TIME_H


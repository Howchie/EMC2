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
									  
	  const double nu_t = nu_vec.back();
    const double omega = pars.omega;
    const RightQuad rq_nu = right_end_quadratic(t_grid, nu_vec,mid);

    // Term 1: The integral from Eq. (10)
    const double termA = 0.5*integrate_pdf_forward_bm_u(t, t_grid, nu_vec, pars, cache_ptr);
    
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

#endif // BM_HITTING_TIME_H


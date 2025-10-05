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
#include <Rcpp.h>
#include <quadmath.h>

using namespace Rcpp;

inline double physical_boundary(double t, const RD_Params& pars) {
  if (pars.fixed_b) {
    return pars.b0;
  }
  return exp_decay_scalar(t, pars.b0, pars.binf, pars.tau, pars.pow);
}

inline double beta_from_time_raw(double t, const RD_Params& pars) {
  const double b_t = physical_boundary(t, pars);
  const double shifted = b_t - pars.mu * t - pars.z0;
  return pars.c * shifted;
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

inline double physical_boundary_derivative(double t, const RD_Params& pars) {
  if (pars.fixed_b) {
    return 0.0;
  }
  if (pars.tau <= 0.0 || pars.pow <= 0.0) {
    return 0.0;
  }
  const double x0 = std::abs(pars.b0);
  const double xinf = std::abs(pars.binf);
  const double amp = x0 - xinf;
  if (std::abs(amp) < FPM_EPSILON) {
    return 0.0;
  }
  if (t <= FPM_EPSILON) {
    return 0.0;
  }
  const double log_ratio = std::log(t) - std::log(pars.tau);
  const double s = safe_exp(pars.pow * log_ratio);
  const double exp_term = safe_exp(-s);
  const double ds_dt = (pars.pow / t) * s;
  const double db_dt = amp * exp_term * (-ds_dt);
  return db_dt;
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
      return (beta_t / pars.c) + pars.mu * t + pars.z0;
    }
  }
  return physical_boundary(t, pars);
}

inline double beta_prime_from_time(double t, const RD_Params& pars) {
  const double b_prime = physical_boundary_derivative(t, pars);
  return pars.c * (b_prime - pars.mu);
}

inline double boundary_prime_from_time(double t,
                                       const RD_Params& pars,
                                       const BoundaryDecayCache* cache = nullptr) {
  if (pars.fixed_b) {
    return 0.0;
  }
  if (cache && cache->has_derivatives()) {
    const double beta_prime = cache->lookup_prime(t);
    if (std::isfinite(beta_prime)) {
      return (beta_prime / pars.c) + pars.mu;
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
  return (beta_prime / pars.c) + pars.mu;
}

// Heat kernel G*(t, dbeta) = exp(-dbeta^2/(2t)) / sqrt(2 pi t)
inline double Gstar(double dt, double dbeta) {
  if (dt <= FPM_EPSILON) return 0.0;
  const double denom = std::sqrt(2.0 * M_PI * dt);
  return safe_exp(-(dbeta * dbeta) / (2.0 * dt)) / denom;
}

inline double kernel_bm(double t,
                        double s,
                        const RD_Params& pars,
                        const BoundaryDecayCache* cache = nullptr) {
    const double dt = t - s;
    if (dt <= FPM_EPSILON) {
        return 0.0;
    }

    const double beta_t = beta_from_time(t, pars, cache);
    const double beta_s = beta_from_time(s, pars, cache);
    
    // Psi(t, t') = b(t) - b(t') in the paper's notation
    const double psi = beta_t - beta_s;

    // Xi(t, t') = exp(-Psi^2 / (2 * (t - t')))
    const double xi = safe_exp(-(psi * psi) / (2.0 * dt));
    
    // omega = sgn(z - b(0)) -> sgn(0 - beta(0))
    // Since we assume z0 < b0, beta(0) is positive, so omega is -1.
    const double omega = -1.0; 

    // The full kernel from Eq. (11), noting that the sqrt(t-s) is handled by the solver
    return omega * psi * xi / std::sqrt(2.0 * M_PI * dt);
}

inline double kernel_bm_backward(double t,      // Current time in the solver grid (from 0 to T)
                                 double s,      // Past time in the solver grid (from 0 to t)
                                 const RD_Params& pars,
                                 const BoundaryDecayCache* cache = nullptr) {
    const double dt = t - s;
    if (dt <= FPM_EPSILON) {
        return 0.0;
    }
    
    // T is the fixed final time for the CDF calculation
    const double T = pars.t_scaled; 

    // The backward method evaluates the boundary at time T-t, T-s, etc.
    // The cache is built for forward-time, so we'll call the raw function here.
    const double beta_T_minus_t = beta_from_time_raw(T - t, pars);
    const double beta_T_minus_s = beta_from_time_raw(T - s, pars);
    
    // Psi(T-t, T-s) = b(T-t) - b(T-s) in the paper's notation
    const double psi = beta_T_minus_t - beta_T_minus_s;

    // Xi(T-t, T-s) = exp(-Psi^2 / (2 * (t - s)))
    // Note the denominator is (t-s), the difference in the integration variable
    const double xi = safe_exp(-(psi * psi) / (2.0 * dt));
    
    const double omega = -1.0; 

    return omega * psi * xi / std::sqrt(2.0 * M_PI * dt);
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
    return -Gstar(t,beta_t);
}

// Unified forcing selector
inline double forcing_bm(double t,
                               const RD_Params& pars,
                               const BoundaryDecayCache* cache = nullptr) {
  //if (!pars.sp_var) {
    // Point start at z0 = S.zL (optionally: assert S.zL ≈ pars.z0)
    return forcing_bm_point(t, pars, cache);
  //} else {
  //  return forcing_bm_uniform(t, pars, cache);
  //}
}

inline double abel_approx_bm(double t,
                             const RD_Params& pars,
                             const BoundaryDecayCache* cache = nullptr) {
  // Use the forcing itself as the small-t seed; the Volterra integral is higher order in t.
  return forcing_bm(std::max(t, 1e-12), pars, cache);
}

// 1. A regularized integrand for the BM, analogous to your theta version
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
    const double H_smooth = (1.0 - psi_sq / dt) * safe_exp(-psi_sq / (2.0 * dt));

    const double denom = std::sqrt(2.0 * M_PI) * dt * std::sqrt(dt);
    if (denom < FPM_EPSILON) return 0.0;

    return (H_smooth * nu_at_prime - nu_at_max) / denom;
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

    // The u-transform: u = sqrt(t_max - t_prime) => dt' = -2u du
    // Integral from t'=0..t_max becomes integral u=sqrt(t_max)..0
    
    double integral_sum = 0.0;
    
    // We integrate over u from u_max=sqrt(t_max) down to u_min=0
    for (int i = N - 1; i >= 0; --i) {
        // Endpoints of the interval in the original t-grid
        const double t_L = t_grid[i];
        const double t_R = t_grid[i+1];

        // Corresponding endpoints in the u-grid
        const double u_L = std::sqrt(t_max - t_L);
        const double u_R = std::sqrt(t_max - t_R);
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

        // Evaluate the regularized integrand f(t') at the three points
        const double f_L = regularized_pdf_integrand_bm(t_max, t_L, nu_at_max, nu_L, pars, cache);
        const double f_R = regularized_pdf_integrand_bm(t_max, t_R, nu_at_max, nu_R, pars, cache);
        const double f_M = regularized_pdf_integrand_bm(t_max, t_M, nu_at_max, nu_M, pars, cache);
        
        // Simpson's rule for this panel in u-space, including the Jacobian term (2u)
        integral_sum += (panel_h_u / 6.0) * ( (2.0 * u_L * f_L) + (8.0 * u_M * f_M) + (2.0 * u_R * f_R) );
    }

    return integral_sum;
}


// This is the function you need to call *after* solving for nu
double calculate_g_from_nu(double t, const std::vector<double>& t_grid, const std::vector<double>& nu_vec, const RD_Params& pars, const BoundaryDecayCache* cache) {
    if (t <= FPM_EPSILON) return 0.0;

    const double nu_t = nu_vec.back();
    const double beta_t = beta_from_time(t, pars, cache);
    const double b_prime_t = beta_prime_from_time(t, pars);
    const double omega = -1.0;
    const double z = pars.z_scaled;

    // Term 1: The new integral from Eq. (10)
    const double term1 = 0.5*integrate_pdf_forward_bm_u(t, t_grid, nu_vec, pars, cache);
    
    // Term 2: -ν(t) / sqrt(2πt)
    const double term2 = -nu_t / std::sqrt(2.0 * M_PI * t);

    // Term 3: The standard inverse-Gaussian-like part
    const double pdf_like_term = ((z - beta_t) * safe_exp(-((z-beta_t) * (z-beta_t)) / (2.0 * t))) / (2.0 * std::sqrt(M_PI * 2.0 * t * t * t));
    
    // Term 4: The b'(t)ν(t) part
    const double b_prime_term = -b_prime_t * nu_t;

    const double term3_and_4 = omega * (b_prime_term + pdf_like_term);

    return term1 + term2 + term3_and_4;
}

double calculate_cdf_from_nu_backward(const std::vector<double>& t_grid,
                                      const std::vector<double>& nu_vec,
                                      const RD_Params& pars) {
    
    const double T = pars.t_scaled;
    const double omega = -1.0;
    double integral = 0.0;

    // Simple trapezoidal integration over the solution
    for (size_t i = 0; i < t_grid.size() - 1; ++i) {
        const double s1 = t_grid[i];
        const double s2 = t_grid[i+1];
        const double ds = s2 - s1;

        auto integrand = [&](double s, double nu_s) {
            const double T_minus_s = T - s;
            if (T_minus_s <= FPM_EPSILON) return 0.0;

            const double beta_T_minus_s = beta_from_time_raw(T_minus_s, pars);

            // This is the integrand from Eq. (12) in transformed space (z=0)
            const double term1 = (0.0 - beta_T_minus_s); // (z - b(T-s))
            const double exp_term = safe_exp(-(term1 * term1) / (2.0 * T_minus_s));
            const double denom = std::sqrt(2.0 * M_PI * T_minus_s * T_minus_s * T_minus_s);

            return (term1 * exp_term / denom) * nu_s;
        };

        integral += (integrand(s1, nu_vec[i]) + integrand(s2, nu_vec[i+1])) * 0.5 * ds;
    }
    
    return std::max(0.0, std::min(1.0, omega * integral));
}

inline RD_Params prepare_bm_params(double t,
                                   double mu,
                                   double sigma,
                                   double z0,
                                   double b0,
                                   double binf,
                                   double tau,
                                   double pow) {
  if (sigma <= 0.0) {
    stop("prepare_bm_params: sigma must be positive.");
  }
  RD_Params pars{};
  pars.b0 = b0;
  pars.binf = binf;
  pars.mu = mu;
  pars.sigma = sigma;
  pars.z0 = z0;
  pars.c = 1.0 / sigma;
  pars.lambda = 1.0;
  pars.theta = 0.0;
  pars.tau = tau;
  pars.pow = pow;
  pars.omega = 1.0;
  pars.t_scaled = t;
  pars.z_scaled = 0.0;
  pars.b_scaled = (b0 - z0) * pars.c;
  pars.fixed_b = (std::abs(binf - b0) <= FPM_EPSILON) || tau <= 0.0 || pow <= 0.0;
  if (pars.fixed_b) {
    pars.tau = 0.0;
    pars.pow = 0.0;
  }
  pars.sp_var = false;
  if (z0 > 0) {
    pars.sp_var = true;
  }

  const double beta0 = (b0 - z0) * pars.c;
  if (beta0 <= FPM_EPSILON) {
    stop("prepare_bm_params: initial position must be below the boundary.");
  }

  return pars;
}

double bm_fht_pdf(double t,
                         double mu,
                         double sigma,
                         double z0,
                         double b0,
                         double binf,
                         double tau,
                         double pow,
                         int num_steps) {
  if (t <= 0.0) {
    return 0.0;
  }
  const int N = (num_steps < 2) ? 2 : num_steps;
  std::vector<double> t_grid(N + 1);
  const double h_t = t / N;
  for (int j = 0; j <= N; ++j) t_grid[j] = j * h_t;

  RD_Params pars = prepare_bm_params(t, mu, sigma, z0, b0, binf, tau, pow);
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

  auto nu = solve_nu_block_with_abel(t, pars, N, t_grid, kernel_fn,
                                          forcing_fn, abel_fn);
  double density = calculate_g_from_nu(t, t_grid, nu, pars, cache_ptr);
  return std::max(0.0, density);
}

inline double bm_fht_cdf(double t,
                         double mu,
                         double sigma,
                         double z0,
                         double b0,
                         double binf,
                         double tau,
                         double pow,
                         int num_steps) {
  if (t <= 0.0) {
    return 0.0;
  }
    const int N = (num_steps < 2) ? 2 : num_steps;
  std::vector<double> t_grid(N + 1);
  const double h_t = t / N;
  for (int j = 0; j <= N; ++j) t_grid[j] = j * h_t;

  RD_Params pars = prepare_bm_params(t, mu, sigma, z0, b0, binf, tau, pow);
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

  auto density = solve_nu_block_with_abel(t, pars, N, t_grid, kernel_fn,
                                          forcing_fn, abel_fn);
  const double integral = integrate_density(t_grid, density);
  return std::max(0.0, std::min(1.0, integral));
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

#endif // BM_HITTING_TIME_H


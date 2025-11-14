#ifndef utils_reducible_diffusion_h
#define utils_reducible_diffusion_h

#define _USE_MATH_DEFINES
#include "utility_functions.h"
#include <cmath>
#include <vector>
#include <map>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <functional>
#include <limits>
#include <cstddef>
#include <utility>
#include <iomanip>
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
// Threshold for switching to Abel approximation, based on the scaled time.
// The paper notes solutions are "visually indistinguishable" up to t=0.02.
const double SMALL_T_SCALED_THRESHOLD = 0.02;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lightweight debug utilities (toggle from R via ou_debug_set)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inline bool OU_DEBUG_ENABLED = false;   // header inline => one definition across TU
inline int  OU_DEBUG_LEVEL   = 1;

// Guarded logging macro
#define OU_DEBUG_LOG(level, expr) \
  do { \
    if (OU_DEBUG_ENABLED && OU_DEBUG_LEVEL >= (level)) { \
      Rcpp::Rcout << std::setprecision(12) << expr << std::endl; \
    } \
  } while(0)

// [[Rcpp::export]]
void ou_debug_set(bool enabled = true, int level = 1) {
  OU_DEBUG_ENABLED = enabled;
  OU_DEBUG_LEVEL   = level;
  Rcpp::Rcout << "[OU DEBUG] enabled=" << (OU_DEBUG_ENABLED?"true":"false")
              << ", level=" << OU_DEBUG_LEVEL << std::endl;
}
// Prepare and scale parameters for the OU hitting time problem.
struct RD_Params {
  double b0, t_scaled, zL_scaled, zU_scaled, b_scaled, binf, c, lambda, theta, tau, pow, omega, beta_prime;
  double mu, sigma, v_max, theta_max, time_scale, z0;
  bool fixed_b, sp_var;
  std::string model;
  BoundaryDecayFn boundary_fn;
  std::vector<double> boundary_params;
};

struct GompertzParams {
    double alpha;    // Rate of reversion (like lambda in OU)
    double beta;     // Volatility term (like sigma in OU)
    double K;        // Carrying capacity / reversion point
    double x0;       // Starting point of the process
    double b_gomp;   // Boundary for the Gompertz process
};

RD_Params transform_gompertz_to_ou(const GompertzParams& g_params, double t) {
    RD_Params ou_params{};

    // 1. Time is unchanged by this transformation.

    // 2. Mean-reversion strength (lambda) is directly equivalent to alpha.
    // The drift term -alpha*Y in the transformed SDE matches -lambda*Y.
    ou_params.lambda = g_params.alpha;

    // 3. Volatility (sigma) is directly equivalent to beta.
    // The noise term beta*dW in the transformed SDE matches sigma*dW.
    ou_params.sigma = g_params.beta;

    // 4. The reversion point (theta) is the most complex. It absorbs
    // the constant drift terms from Ito's lemma and the original SDE.
    // The constant drift is: alpha*ln(K) - beta^2/2, which must equal lambda*theta.
    // theta = (alpha*ln(K) - beta^2/2) / lambda
    // Since lambda = alpha, this simplifies.
    ou_params.theta = std::log(g_params.K) - (g_params.beta * g_params.beta) / (2.0 * g_params.alpha);

    // 5. The starting point (z0) and boundary (b0) are log-transformed.
    ou_params.z0 = std::log(g_params.x0);
    ou_params.b0 = std::log(g_params.b_gomp);

    return ou_params;
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
        // TODO implement safe fallbacks for <5 points instead of throwing error
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


// Re-worked to use maps to allow arbitrary additions to the cache
struct BoundaryDecayCache {
    mutable std::map<double, double> beta_map;
    mutable std::map<double, double> beta_prime_map;

    bool empty() const {
        return beta_map.empty();
    }

    bool has_derivatives() const {
        return !beta_prime_map.empty();
    }

	double lookup(double theta) const {
        auto it = beta_map.find(theta);
        return (it != beta_map.end()) ? it->second : std::numeric_limits<double>::quiet_NaN();
    }
 
	double lookup_prime(double theta) const {
        auto it = beta_prime_map.find(theta);
        return (it != beta_prime_map.end()) ? it->second : std::numeric_limits<double>::quiet_NaN();
    }
	
	void put_beta(double theta, double value) const {
        beta_map.emplace(theta, value);
    }

    void put_beta_prime(double theta, double value) const {
        beta_prime_map.emplace(theta, value);
    }
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

double theta_to_t(double theta) {
    if (theta <= -1.0) return std::numeric_limits<double>::infinity();
    return std::log1p(theta);
}

double v_to_t(double v) {
    if (v >= 1.0) return std::numeric_limits<double>::infinity();
    return -std::log(1.0 - v);
}

inline double beta_from_v(double v, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr) {
    const double scale = std::max(0.0, 1.0 - v);
    if (scale <= 0.0) {
        return 0.0;
    }
    if (pars.fixed_b) {
        return scale * pars.b_scaled;
    }
    if (cache && !cache->empty()) {
        const double cached = cache->lookup(v);
        if (std::isfinite(cached)) {
            return cached;
        }
    }
    const double t = v_to_t(v) / pars.lambda;
    const double t_max = pars.t_scaled / pars.lambda;
    const double bt = -pars.omega*evaluate_boundary_decay(t_max - t, pars);
    const double bt_scaled =  (pars.c * (bt - pars.theta));
    const double beta = scale * bt_scaled;
	if (cache) cache->put_beta(v, beta);
    return beta;
}

inline double beta_from_theta(double theta, const RD_Params& pars, const BoundaryDecayCache* cache = nullptr,const util::AkimaSpline* spline = nullptr) {
    const double scale = std::max(0.0, 1.0 + theta);
    if (theta <= 0.0) {
        return 0.0;
    }
	if (pars.fixed_b) {
        return scale * pars.b_scaled;
    }
    if (cache) {
        const double cached = cache->lookup(theta);
        if (std::isfinite(cached)) {
            return cached;
        }
    }
	
	// Try spline if provided
    if (spline) {
        const double interp = spline->interpolate(theta);
        if (std::isfinite(interp)) {
            return interp;
        }
	}

	// If not cached and no spline is built, compute the value
    const double t = theta_to_t(theta) / pars.lambda;
    const double bt = -pars.omega * evaluate_boundary_decay(t, pars); // always use a positive boundary to compute decay for convenience. Omega is negative when hitting from above, and we always initialise at z=0 therefore we can be confident in the sign of the raw boundary aligning here
    const double bt_scaled = pars.c * (bt - pars.theta);
	const double beta = scale * bt_scaled;
	if (cache) cache->put_beta(theta, beta);
    return beta;
}

// Utility: build a central-difference stencil around theta while respecting bounds.
inline bool make_central_difference_times(double theta,
                                          double lower_bound,
                                          double upper_bound,
                                          double& tL,
                                          double& tR) {
    const double rel_scale = 1.0 + std::abs(theta);
    double step = 1e-4 * rel_scale;
    if (step > 1e-3) step = 1e-3;
    if (step < 1e-7) step = 1e-7;

    tL = theta - step;
    tR = theta + step;
    if (tL < lower_bound) {
        tL = lower_bound;
        tR = std::min(upper_bound, tL + 2.0 * step);
    }
    if (tR > upper_bound) {
        tR = upper_bound;
        tL = std::max(lower_bound, tR - 2.0 * step);
    }
    return tR > tL;
}

// Utility: compute beta'(theta) at a single point using a robust central difference.
inline double beta_prime_at_theta(double theta, const RD_Params& pars, double theta_max, const BoundaryDecayCache* cache = nullptr) {
    if (pars.fixed_b) {
        return pars.b_scaled;
    }
	if (cache) {
        const double cached = cache->lookup_prime(theta);
        if (std::isfinite(cached)) {
            return cached;
        }
    }
	// Compute derivative via finite difference
    double tL = 0.0;
    double tR = 0.0;
    if (!make_central_difference_times(theta, 0.0, theta_max, tL, tR)) {
        return 0.0;
    }
    const double bL = beta_from_theta(tL, pars, cache);
    const double bR = beta_from_theta(tR, pars, cache);
    const double denom = tR - tL;
	const double beta_prime = (std::abs(denom) > FPM_EPSILON) ? ((bR - bL) / denom) : 0.0;
    if (cache) cache->put_beta_prime(theta, beta_prime);
    return beta_prime;
}

inline BoundaryDecayCache make_boundary_decay_cache(const std::vector<double>& theta_vals,
                                                    const RD_Params& pars) {
    BoundaryDecayCache cache;
    if (pars.fixed_b || theta_vals.empty()) {
        return cache;
    }
    const double theta_max = *std::max_element(theta_vals.begin(), theta_vals.end());
    for (double theta : theta_vals) {
        (void)beta_from_theta(theta, pars, &cache);
        (void)beta_prime_at_theta(theta, pars, theta_max, &cache);
    }
    return cache;
}

double local_quad_derivative(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  int j,
                                  double eval_point) {
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
    return f01 + f012 * ((eval_point - x1) + (eval_point - x0));
}

namespace util {

/**
 * @brief Computes the derivative at a point using a Savitzky-Golay filter.
 *
 * This method fits a polynomial of a given order to a window of data points
 * surrounding the target point using least-squares. The derivative of that
 * fitted polynomial at the target point is returned. This is robust to noise.
 * This implementation is for non-uniformly spaced data.
 *
 * @param x The vector of independent variable values (must be sorted).
 * @param y The vector of dependent variable values.
 * @param index The index of the point at which to compute the derivative.
 * @param window_size The number of points in the fitting window (must be odd, e.g., 5, 7).
 * @param poly_order The order of the polynomial to fit (e.g., 2 for quadratic).
 * @return The estimated derivative at x[index].
 */
double savitzky_golay_derivative(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 int index,
                                 double eval_point,
                                 int window_size,
                                 int poly_order)
{
    // Validation and basic setup
    const int n = static_cast<int>(x.size());
    if (n == 0 || y.size() != x.size()) return std::numeric_limits<double>::quiet_NaN();

    // Require odd window; clamp to available data if needed
    if (window_size < 3) window_size = 3;
    if ((window_size & 1) == 0) ++window_size; // make odd
    if (window_size > n) window_size = (n % 2 == 0 ? n - 1 : n);

    if (poly_order != 2 || poly_order >= window_size) {
        // Only quadratic supported for now
        return std::numeric_limits<double>::quiet_NaN();
    }

    const int half_window = window_size / 2;
    // Prefer the window centered around the actual evaluation location.
    int k = 0;
    {
        // k = index of the last x <= eval_point (clamped)
        auto it = std::upper_bound(x.begin(), x.end(), eval_point);
        k = static_cast<int>(std::distance(x.begin(), (it == x.begin() ? it : std::prev(it))));
        if (k < 0) k = 0;
        if (k > n - 1) k = n - 1;
    }

    // Near-boundary: fall back to a local 3-point derivative centered near eval_point.
    if (k < half_window || k > n - 1 - half_window) {
        if (n >= 3) {
            return local_quad_derivative(x, y, k, eval_point);
        }
        if (n > 1) {
            // fall back to one-sided difference around the nearest available segment
            if (k == 0) return (y[1] - y[0]) / (x[1] - x[0]);
            return (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
        }
        return 0.0;
    }

    // Choose a contiguous window that slides to the boundary when needed.
    // This ensures a proper right/left-sided fit near the ends.
    int start = k - half_window;
    int end = k + half_window;
    if (start < 0) { end -= start; start = 0; }
    if (end > n - 1) { start -= (end - (n - 1)); end = n - 1; }
    if (start < 0) start = 0; // final clamp

    // Center the polynomial at the actual evaluation point to avoid bias
    const double x_center = eval_point;

    // Build normal equations for quadratic LS fit: y ~= c0 + c1*(x-xc) + c2*(x-xc)^2
    double S0=0.0, S1=0.0, S2=0.0, S3=0.0, S4=0.0;
    double T0=0.0, T1=0.0, T2=0.0;
    for (int i = start; i <= end; ++i) {
        const double dx = x[i] - x_center;
        const double yi = y[i];
        const double dx2 = dx * dx;
        const double dx3 = dx2 * dx;
        const double dx4 = dx2 * dx2;

        S0 += 1.0;
        S1 += dx;
        S2 += dx2;
        S3 += dx3;
        S4 += dx4;

        T0 += yi;
        T1 += yi * dx;
        T2 += yi * dx2;
    }

    // Solve 3x3 via Cramer's rule. If ill-conditioned, fall back to local quadratic.
    const double det = S0*(S2*S4 - S3*S3) - S1*(S1*S4 - S2*S3) + S2*(S1*S3 - S2*S2);
    if (!std::isfinite(det) || std::abs(det) < 1e-15) {
        // Use a simple 3-point quadratic derivative around the nearest triple
        return local_quad_derivative(x, y, k, eval_point);
    }
    const double inv_det = 1.0 / det;

    // derivative at eval_point is c1 when centered at x_center = eval_point
    const double det1 = S0*(T1*S4 - S3*T2) - T0*(S1*S4 - S2*S3) + S2*(S1*T2 - S2*T1);
    const double c1 = det1 * inv_det;

    return c1;
}

}; // namespace util

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

// No merging in dense history mode

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
    return panel;
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
// No refined split needed in dense history mode
int seed_nu_on_grid(
    const std::vector<double>& grid,
    const RD_Params& pars,
    std::vector<double>& F,
    const AbelFn& abel,
    double t_cut = SMALL_T_SCALED_THRESHOLD)
{
    const int N = static_cast<int>(grid.size());
    int J = 0;
    OU_DEBUG_LOG(2, "seed_nu_on_grid: N=" << N << ", t_cut=" << t_cut);
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
    OU_DEBUG_LOG(2, "seed_nu_on_grid: seeded up to index=" << (J-1)
                    << ", first_unseeded_index=" << J
                    << ", last_seed_val_at=" << (J>0? grid[J-1] : grid[0]));
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
    OU_DEBUG_LOG(2, "seed_nu_on_nonuniform_grid: N=" << N << ", t_cut=" << t_cut);
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
        if (first_unseeded > 0) {
            --first_unseeded;
        } else {
            first_unseeded = 0;
        }
    }
    if (first_unseeded >= N) {
        first_unseeded = N - 1;
        if (first_unseeded % 2 != 0 && first_unseeded > 0) {
            --first_unseeded;
        }
    }
    info.first_unseeded = std::clamp(first_unseeded, 0, N - 1);
    OU_DEBUG_LOG(2, "seed_nu_on_nonuniform_grid: last_seeded_index=" << info.last_seeded
                    << ", first_unseeded_index=" << info.first_unseeded
                    << ", last_seed_val_at=" << (info.last_seeded>=0? grid[info.last_seeded] : 0.0));
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
    OU_DEBUG_LOG(1, "solve_nu_block_by_block_impl: t_max=" << t_max
                    << ", N=" << N << ", h=" << h
                    << ", k0_seeded=" << k0_seeded);

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
        OU_DEBUG_LOG(5, "block m=" << m
                        << ", t0=" << t0 << ", t1=" << t1 << ", t2=" << t2
                        << ", S1=" << S1 << ", S2=" << S2
                        << ", F1_pred=" << F1_pred << ", F2_pred=" << F2_pred
                        << ", F1_corr=" << F[j1] << ", F2_corr=" << F[j2]);
    }
    return F;
}

// Dense non-uniform history sum: sum quadratic product-integration over all
// native panels strictly before tn, up to (but excluding) panel index cut_index.
// This restores the exact O(N^2) history on a non-uniform grid.
inline double stieltjes_history_dense_nonuniform(
    int cut_index_exclusive,      // panels strictly before this node
    int max_known_index,          // last index i such that F[i] is known
    double tn,
    const std::vector<double>& grid,
    const std::vector<double>& F,
    const RD_Params& pars,
    const KernelFn& K)
{
    const int total_nodes = static_cast<int>(grid.size());
    if (total_nodes < 2) return 0.0;

    int last_panel = cut_index_exclusive - 1;
    if (last_panel > total_nodes - 2) {
        last_panel = total_nodes - 2;
    }
    if (last_panel < 0) {
        return 0.0;
    }

    double S = 0.0;
    for (int i = 0; i <= last_panel; ++i) {
        const double a = grid[i];
        const double b = grid[i + 1];
        if (b >= tn - FPM_EPSILON) break;

        const double Fa = F[i];
        const double Fb = F[i + 1];

        bool have_mid = false;
        double Fm = 0.5 * (Fa + Fb);

        // Interior panel midpoint uses (i-1, i, i+1) if all known
        if (i > 0 && (i + 1) <= max_known_index) {
            Fm = (-0.125) * F[i - 1] + 0.75 * F[i] + 0.375 * F[i + 1];
            have_mid = true;
        }
        // First panel midpoint uses (0,1,2) if all known
        else if (i == 0 && (i + 2) <= max_known_index) {
            Fm = 0.375 * F[i] + 0.75 * F[i + 1] - 0.125 * F[i + 2];
            have_mid = true;
        }

        const Panel P{a, 0.5 * (a + b), b, b - a};
        const StieltjesMoments M = stieltjes_moments(tn, P.a, P.b);

        if (have_mid) {
            const LagrangeWeights W = stieltjes_lagrange_weights(P, M);
            const double Qa = K(tn, a, pars) * Fa;
            const double Qm = K(tn, P.mid, pars) * Fm;
            const double Qb = K(tn, b, pars) * Fb;
            S += W.w0 * Qa + W.w1 * Qm + W.w2 * Qb;
        } else {
            // Fallback to endpoint-only rule if midpoint would peek forward
            const double Wa = (P.b * M.J0 - M.J1) / P.d;
            const double Wb = (M.J1 - P.a * M.J0) / P.d;
            S += Wa * (K(tn, a, pars) * Fa) + Wb * (K(tn, b, pars) * Fb);
        }
    }
    OU_DEBUG_LOG(6, "stieltjes_history_dense_nonuniform: tn=" << tn
                    << ", last_panel=" << last_panel << ", S=" << S);
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
    OU_DEBUG_LOG(1, "solve_nu_block_with_abel: t_max=" << t_max
                    << ", num_steps=" << num_steps
                    << ", seeded_up_to_index=" << k0);
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
    double ratio = 2.0;       // h_{m+1} = ratio * h_m
    int    base_panels = 128; // first chunk panel count
    int    max_chunks = 12;   // max number of chunks
    // Optional: ensure at least this many panels before theta_min
    double theta_min = 0.0;
    int    pre_min_panels = 0;
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

inline std::vector<ChunkSpec> build_chunked_theta_layout_with_prelude(
    double theta_max,
    int uniform_num_steps,
    const ChunkingOptions& opts)
{
    // Fallback if no prelude requested
    if (!(opts.theta_min > 0.0) || opts.pre_min_panels <= 0) {
        return build_chunked_theta_layout(theta_max, uniform_num_steps, opts);
    }

    std::vector<ChunkSpec> chunks;
    if (!(theta_max > 0.0) || uniform_num_steps <= 0) {
        return chunks;
    }

    int N = uniform_num_steps;
    if (N % 2 != 0) ++N;
    int P = opts.pre_min_panels;
    if (P % 2 != 0) ++P;
    // Leave at least 2 panels for the remainder
    if (P > N - 2) P = std::max(2, N - 2);

    const double theta_min = std::min(opts.theta_min, theta_max);

    // Prelude chunk [0, theta_min]
    {
        ChunkSpec pre;
        pre.start = 0.0;
        pre.end = theta_min;
        pre.num_steps = P;
        pre.h = (pre.end - pre.start) / static_cast<double>(P);
        pre.start_index = 0;
        chunks.push_back(pre);
    }

    // Remainder [theta_min, theta_max] with N_tail panels, geometric growth
    int N_tail = std::max(2, N - P);
    if (N_tail % 2 != 0) ++N_tail;

    const double ratio = std::max(1.0, opts.ratio);
    int panels_per_chunk = std::max(4, opts.base_panels);
    if (panels_per_chunk % 2 != 0) ++panels_per_chunk;
    int max_chunks = std::max(1, opts.max_chunks);
    const double base_h = (theta_max - theta_min) / static_cast<double>(N_tail);
    if (ratio > 1.0 && theta_min > 0.0) {
    const double pre_h = theta_min / static_cast<double>(P);
        if (pre_h > 0.0 && base_h > pre_h * (1.0 + 1e-12)) {
            const double growth_needed = base_h / pre_h;
            const double required = std::log(growth_needed) / std::log(ratio);
            if (std::isfinite(required)) {
                const int needed_chunks = static_cast<int>(std::ceil(required)) + 1;
                max_chunks = std::max(max_chunks, std::min(needed_chunks, 64));
                OU_DEBUG_LOG(1, "needed_chunks=" << needed_chunks
                                << ", max_chunks=" << max_chunks << std::endl);
            }
        }
    }
    const double pre_h = (theta_min > 0.0 && P > 0) ? theta_min / static_cast<double>(P) : base_h;
    double start = theta_min;
    int start_index = P;
    double current_h = base_h;
    if (ratio > 1.0 && max_chunks > 1) {
        const double scale = std::pow(ratio, max_chunks - 1);
        if (std::isfinite(scale) && scale > 0.0) {
            const double candidate_h = base_h / scale;
            if (candidate_h > 0.0) {
                current_h = std::max(candidate_h, pre_h);
            }
        }
    }
    for (int chunk_id = 0; chunk_id < max_chunks && start < theta_max - 1e-12; ++chunk_id) {
        int panels = panels_per_chunk;
        const double span = panels * current_h;
        double end = start + span;

        const bool is_last = (chunk_id == max_chunks - 1) || (end >= theta_max - 1e-12);
        if (is_last) {
            const double remaining = theta_max - start;
            if (remaining <= 0.0) break;
            panels = std::max(2, static_cast<int>(std::ceil(remaining / std::max(current_h, FPM_EPSILON))));
            if (panels % 2 != 0) ++panels;
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
        current_h = std::min(current_h * ratio, base_h);
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

    const auto chunks = (opts.theta_min > 0.0 && opts.pre_min_panels > 0)
        ? build_chunked_theta_layout_with_prelude(theta_max, uniform_num_steps, opts)
        : build_chunked_theta_layout(theta_max, uniform_num_steps, opts);
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

// Chunked theta-grid Volterra solver using dense (non-uniform) history sums.
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
{
    if (t_max <= 0.0) {
        grid.assign(1, 0.0);
        return std::vector<double>{G(0.0, pars)};
    }
    if (uniform_num_steps < 2) uniform_num_steps = 2;
    if (uniform_num_steps % 2 != 0) ++uniform_num_steps;

    auto chunks = (opts.theta_min > 0.0 && opts.pre_min_panels > 0)
        ? build_chunked_theta_layout_with_prelude(t_max, uniform_num_steps, opts)
        : build_chunked_theta_layout(t_max, uniform_num_steps, opts);
    if (chunks.empty()) {
        ChunkSpec chunk; chunk.start=0.0; chunk.end=t_max; chunk.h=t_max/ uniform_num_steps; chunk.num_steps=uniform_num_steps; chunk.start_index=0; chunks.push_back(chunk);
    }

    grid = build_chunked_theta_grid(t_max, uniform_num_steps, opts);
    const int total_nodes = static_cast<int>(grid.size());
    OU_DEBUG_LOG(1, "solve_nu_chunked_with_abel: t_max=" << t_max
                    << ", uniform_num_steps=" << uniform_num_steps
                    << ", total_nodes=" << total_nodes
                    << ", opts(ratio=" << opts.ratio
                    << ", base_panels=" << opts.base_panels
                    << ", max_chunks=" << opts.max_chunks
                    << ", theta_min=" << opts.theta_min
                    << ", pre_min_panels=" << opts.pre_min_panels << ")");
    for (size_t ci=0; ci<chunks.size(); ++ci) {
      const auto& c = chunks[ci];
      OU_DEBUG_LOG(3, "  chunk[" << ci << "]: start=" << c.start << ", end=" << c.end
                      << ", h=" << c.h << ", panels=" << c.num_steps
                      << ", start_index=" << c.start_index);
    }

    std::vector<double> F(total_nodes, 0.0);
    F[0] = G(grid[0], pars);

    SeedInfo seed = seed_nu_on_nonuniform_grid(grid, pars, F, abel, t_cut);
    int first_unseeded = std::min(seed.first_unseeded, total_nodes - 1);
    int last_seeded = std::min(std::max(seed.last_seeded, 0), total_nodes - 1);
    if (first_unseeded >= total_nodes - 1) return F;

    for (const auto& chunk : chunks) {
        const int chunk_start = chunk.start_index;
        const int chunk_end = chunk.end_index();
        if (chunk.num_steps < 2) continue;

        int first_block = std::max(chunk_start / 2, first_unseeded / 2);
        const int block_end = chunk_end / 2;

        const double h = chunk.h;
        double a1, b1, g1, a2, b2, g2;
        block_coeffs(h, 0.5 * h, a1, b1, g1);
        block_coeffs(2.0 * h, h, a2, b2, g2);

        for (int m = first_block; m < block_end; ++m) {
            const int j0 = 2 * m;
            const int j1 = j0 + 1;
            const int j2 = j0 + 2;
            if (j2 >= total_nodes) break;
            if (j0 < chunk_start || j2 > chunk_end) continue;
            if (j1 <= last_seeded && j2 <= last_seeded) continue;

            const double t0 = grid[j0];
            const double t1 = grid[j1];
            const double t2 = grid[j2];
            const double tmid = t0 + 0.5 * h;

            const double S1 = stieltjes_history_dense_nonuniform(j0, j0, t1, grid, F, pars, K);
            const double S2 = stieltjes_history_dense_nonuniform(j0, j0, t2, grid, F, pars, K);

            const double K1_0 = K(t1, t0, pars);
            const double K1_mid = K(t1, tmid, pars);
            const double K1_1 = K(t1, t1, pars);
            const double K2_0 = K(t2, t0, pars);
            const double K2_1 = K(t2, t1, pars);
            const double K2_2 = K(t2, t2, pars);
            const double G1 = G(t1, pars);
            const double G2 = G(t2, pars);

            double F1_pred = 0.0;
            if (j0 >= 1) {
                double Fmid_pred = (-0.125) * F[j0 - 1] + 0.75 * F[j0];
                if (j1 <= last_seeded) {
                    Fmid_pred += 0.375 * F[j1];
                }
                const double den1_pred = 1.0 - g1 * K1_1 - (3.0 / 8.0) * b1 * K1_mid;
                const double rhs1_pred = G1 + S1 + a1 * K1_0 * F[j0] + b1 * K1_mid * Fmid_pred;
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
            const double C = G1 + S1 + a1 * K1_0 * F[j0] + (3.0 / 8.0) * b1 * K1_mid * F[j0];
            const double D = -b2 * K2_1;
            const double E = 1.0 - g2 * K2_2;
            const double R = G2 + S2 + a2 * K2_0 * F[j0];

            const double det = A * E - B * D;
            const double inv_det = (std::abs(det) < FPM_EPSILON) ? (1.0 / FPM_EPSILON) : (1.0 / det);
            const double F1_corr = (E * C - B * R) * inv_det;
            const double F2_corr = (-D * C + A * R) * inv_det;

            F[j1] = F1_corr;
            F[j2] = F2_corr;

            last_seeded = std::max(last_seeded, j2);
            OU_DEBUG_LOG(5, "chunk step: j0=" << j0 << ", j1=" << j1 << ", j2=" << j2
                            << ", t0=" << t0 << ", t1=" << t1 << ", t2=" << t2
                            << ", S1=" << S1 << ", S2=" << S2
                            << ", F1=" << F[j1] << ", F2=" << F[j2]);
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

// Computes the uniform-in-X-space averaged image weight for the Gompertz case,
// which is equivalent to an exp(y)-weighted average in Y-space (OU/Wiener space).
// This is the direct replacement for your 'averaged_image_term' function.
inline double averaged_image_term_exp(double t, double beta_t, const RD_Params& pars) {
    if (t <= FPM_EPSILON) return 0.0;

    // Point-start case is identical to the original function
    if (!pars.sp_var) {
        const double num = pars.zU_scaled - beta_t; // zU_scaled is log(x0)
        return (num / t) * Gstar(t, num);
    }

    double y_lo = pars.zL_scaled;
    double y_hi = pars.zU_scaled;
    if (y_hi < y_lo) std::swap(y_lo, y_hi);
    
    // Normalization constant for the p(y) = e^y PDF
    const double norm_const = std::exp(y_hi) - std::exp(y_lo);
    if (norm_const <= FPM_EPSILON) {
        return 0.0;
    }

    // This term comes from completing the square in the exponent
    const double exp_factor = std::exp(beta_t + t / 2.0);

    // The new mean of the Gaussian after completing the square
    const double new_mean = beta_t + t;

    // We evaluate the Gaussian PDF and CDF with this new mean
    const double G_hi = Gstar(t, y_hi - new_mean);
    const double G_lo = Gstar(t, y_lo - new_mean);

    const double I_hi = Gstar_Integral(t, new_mean, -std::numeric_limits<double>::infinity(), y_hi);
    const double I_lo = Gstar_Integral(t, new_mean, -std::numeric_limits<double>::infinity(), y_lo);

    // The integral of (y - b) * exp(...) splits into two parts:
    // (new_mean - b) * integral[G*] - t * [G]_lo^hi
    const double term1 = (new_mean - beta_t) * (I_hi - I_lo);
    const double term2 = -t * (G_hi - G_lo);

    return (1.0 / norm_const) * exp_factor * (term1 + term2);
}

// Computes the average of f(y) = (y - beta) * exp(-(y - beta)^2 / v)
// over the interval [y_lo, y_hi]
inline double average_of_exp_deriv(double beta, double variance, double y_lo, double y_hi) {
    // This check handles all divide-by-zero or bad exp() calls
    if (variance <= FPM_EPSILON) {
        return 0.0;
    }

    // The core exponential term: exp(-(y - beta)^2 / v)
    auto exp_term = [beta, variance](double y) {
        const double delta = y - beta;
        return std::exp(-(delta * delta) / variance);
    };

    // The point-start function: f(y) = (y - beta) * exp_term(...)
    auto point_start_f = [beta, &exp_term](double y) {
        return (y - beta) * exp_term(y);
    };

    // The antiderivative: F(y) = -0.5 * v * exp_term(...)
    auto antideriv_F = [variance, &exp_term](double y) {
        return -0.5 * variance * exp_term(y);
    };


    // --- Main logic ---
    if (y_hi < y_lo) {
        std::swap(y_lo, y_hi);
    }
    const double span = y_hi - y_lo;

    // Point-start case: return f(y_lo)
    if (span <= FPM_EPSILON) {
        return point_start_f(y_lo);
    }

    // Averaged case: return (F(y_hi) - F(y_lo)) / span
    const double F_hi = antideriv_F(y_hi);
    const double F_lo = antideriv_F(y_lo);
    
    return (F_hi - F_lo) / span;
}


/**
 * @brief Determines the number of integration steps for a given t_max.
 *
 * This function calculates the number of steps for the Volterra solver.
 * It uses a logarithmic scaling law to balance accuracy and computational cost.
 * The number of steps grows with `t_max`, but slower than linearly,
 * which prevents excessively large step counts for large time horizons.
 * This is based on the observation that for the OU process, a finer grid
 * is more critical at earlier times.
 *
 * The scaling is designed to be equivalent to a linear scaling (`t_max / fineness`)
 * at `t_max = 1.0`, while growing logarithmically thereafter.
 *
 * It also enforces a minimum number of steps and ensures the final count is
 * even, as required by the underlying numerical solvers.
 *
 * @param t_max The maximum time (or theta) for the integration.
 * @param fineness A parameter controlling the grid density. A smaller value
 *                 leads to more steps. Conceptually, `1/fineness` is the
 *                 number of steps at `t_max = 1.0`.
 * @param min_steps The minimum number of steps to use.
 * @return The calculated (even) number of steps.
 */
// [[Rcpp::export]]
int calculate_num_steps(double t_max, double fineness = 0.01, int min_steps = 100) {
    // Logarithmic scaling: equivalent to linear at t_max=1, but grows slower.
    // The term `(std::exp(1.0) - 1.0)` is used to scale the argument of log1p
    // such that for t_max = 1, log1p(t_max * (e-1)) = log(e) = 1.
    const double e_minus_1 = std::exp(1.0) - 1.0;
    double ideal_steps_double = std::ceil(
        (1.0 / fineness) * std::log1p(t_max * e_minus_1)
    );

    // Enforce the minimum number of steps
    int num_steps = static_cast<int>(
        std::max(static_cast<double>(min_steps), ideal_steps_double)
    );

    // Ensure the number of steps is even
    if (num_steps % 2 != 0) {
        num_steps++;
    }

    num_steps = (num_steps < 2) ? 2 : num_steps;

    return num_steps;
}




#endif // utils_reducible_diffusion_h


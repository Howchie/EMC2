#ifndef gaussian_h
#define gaussian_h

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>
#include "utility_functions.h"
#include "wald_functions.h"  // pnorm_std / fast_norm_phi (USE_FAST_PNORM fast path)
using namespace Rcpp;

#ifndef FPM_EPSILON
#define FPM_EPSILON 1e-12
#endif

constexpr double sqrt2 = 1.4142135623730950488;
constexpr double sqrt_twoPI = 2.5066282746310005024;
constexpr double inv_sqrt2_pi = 0.39894228040143267794;
constexpr double log_twoPI = 1.83787706640934548356;
constexpr double inv3sqrt2pi = 0.13298076013381089265;
constexpr double fourPI = 12.566370614359172953;
constexpr double minus_inv_twoPI = -0.15915494309189533577;

// exp(x) as straight-line arithmetic, so that `#pragma omp simd` loops over
// quadrature nodes vectorise.  The package builds with -fno-finite-math-only,
// which keeps __FAST_MATH__ undefined and therefore keeps GCC from calling
// glibc's vector exp (libmvec): std::exp in a node loop stays one scalar call
// per node.  Cody-Waite reduction by ln 2 (split constant, fma, so the
// reduction is exact to well below an ulp for |k| < 2^11), a degree-12 Taylor
// polynomial on |r| <= ln2/2 (truncation < 2e-16), and the scale 2^k applied
// as two factors 2^k1 * 2^k2 so that results in the subnormal range round
// correctly instead of wrapping the exponent field.  Measured over 2e7 random
// arguments in [-745, 709]: max relative error 3.2e-16 (1.4 ulp) in the normal
// range and no mismatch against std::exp in the subnormal range; x < -1000
// returns 0 and x >= 709.79 overflows to +Inf, as std::exp does.
// It must be inlined into the node loop to vectorise: GCC declined to inline
// it into the (large) BVN kernels on its own, which left the `omp simd` loops
// as scalar calls slower than glibc's exp.
#if defined(__GNUC__)
#define EMC2_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define EMC2_ALWAYS_INLINE inline
#endif
EMC2_ALWAYS_INLINE double emc2_exp_simd(double x) {
  x = x < -1000.0 ? -1000.0 : (x > 710.0 ? 710.0 : x);
  const double k = std::rint(x * 1.4426950408889634074);
  double r = std::fma(-k, 6.93147180369123816490e-01, x);
  r = std::fma(-k, 1.90821492927058770002e-10, r);
  double p = 1.0 / 479001600.0;
  p = std::fma(p, r, 1.0 / 39916800.0);
  p = std::fma(p, r, 1.0 / 3628800.0);
  p = std::fma(p, r, 1.0 / 362880.0);
  p = std::fma(p, r, 1.0 / 40320.0);
  p = std::fma(p, r, 1.0 / 5040.0);
  p = std::fma(p, r, 1.0 / 720.0);
  p = std::fma(p, r, 1.0 / 120.0);
  p = std::fma(p, r, 1.0 / 24.0);
  p = std::fma(p, r, 1.0 / 6.0);
  p = std::fma(p, r, 0.5);
  p = std::fma(p, r, 1.0);
  p = std::fma(p, r, 1.0);
  // Biased exponents k_i + 1023 lie in [301, 1535] for the clamped range, so
  // adding 2^52 leaves each as the integer in the low mantissa bits.
  const double k1 = std::floor(0.5 * k);
  const double k2 = k - k1;
  const double m1 = (k1 + 1023.0) + 4503599627370496.0;
  const double m2 = (k2 + 1023.0) + 4503599627370496.0;
  uint64_t b1, b2;
  std::memcpy(&b1, &m1, sizeof(b1));
  std::memcpy(&b2, &m2, sizeof(b2));
  b1 = (b1 & 0xFFFFFFFFFFFFFull) << 52;
  b2 = (b2 & 0xFFFFFFFFFFFFFull) << 52;
  double s1, s2;
  std::memcpy(&s1, &b1, sizeof(s1));
  std::memcpy(&s2, &b2, sizeof(s2));
  return p * s1 * s2;
}

// The BVN node loops call EMC2_BVN_EXP.  Build with -DEMC2_BVN_SCALAR_EXP to
// fall back to std::exp there (benchmarking; the loops are then scalar).
#ifdef EMC2_BVN_SCALAR_EXP
#define EMC2_BVN_EXP(x) std::exp(x)
#else
#define EMC2_BVN_EXP(x) emc2_exp_simd(x)
#endif

// Useful Distribution Functions
inline double gaussian_pdf(double x, double mean = 0.0, double var = 1.0, bool log_p = false) {
  if (!(var > 0.0) || !emc2_isfinite(var)) return log_p ? R_NegInf : 0.0;;
  
  const double z = x - mean;

  if (log_p) {
    return -0.5 * (log_twoPI + std::log(var) + z * z / var);
  }

  return inv_sqrt2_pi / std::sqrt(var) * std::exp(-0.5 * z * z / var);
}

inline double gaussian_cdf(double x, double mean = 0.0, double var = 1.0, bool log_p=false) {
  if (!(var > 0.0) || !emc2_isfinite(var)) {
    if (log_p) {
      return (x < mean) ? R_NegInf : 0.0;
    }
    return (x < mean) ? 0.0 : 1.0;
  }
  return pnorm_std((x - mean) / std::sqrt(var), true, log_p);
}

// Heat kernel G*(t, delta) = N(delta | 0, t)
inline double Gstar(double var, double delta, bool log_p = false) {
  if (var <= 0.0) return log_p ? R_NegInf : 0.0;
  return gaussian_pdf(delta, 0.0, var, log_p);
}

// Erlang-n kill survival: log S_K^(n)(t) = log(exp(-lambda*t) * sum_{m=0}^{n-1} (lambda*t)^m / m!)
// n=3 is mixture: omega*E(1, lam) + (1-omega)*E(2, 2*lam)
inline double erlang_log_surv(double t, double lambda, int n, double omega = 1.0) {
  if (lambda <= 0.0 || t <= 0.0) return 0.0;
  if (t == R_PosInf) return R_NegInf;
  if (n == 3) {
    const double s1 = erlang_log_surv(t, lambda, 1);
    const double s2 = erlang_log_surv(t, 2.0 * lambda, 2);
    const double w = std::fmax(0.0, std::fmin(1.0, omega));
    if (w >= 1.0) return s1;
    if (w <= 0.0) return s2;
    return log_sum_exp(std::log(w) + s1, std::log1p(-w) + s2);
  }
  if (n <= 1) return -lambda * t;
  if (n == 2) return -lambda * t + std::log1p(lambda * t);

  // General n: use a loop for the sum to avoid factorials
  double sum_term = 1.0;
  double current_term = 1.0;
  double lt = lambda * t;
  for (int m = 1; m < n; ++m) {
    current_term *= (lt / m);
    sum_term += current_term;
  }
  return -lambda * t + std::log(sum_term);
}

// Erlang-n kill density: f_K^(n)(t)
inline double erlang_log_pdf(double t, double lambda, int n, double omega = 1.0) {
  if (lambda <= 0.0 || t < 0.0) return R_NegInf;
  if (t == R_PosInf) return R_NegInf;
  if (n == 3) {
    const double f1 = erlang_log_pdf(t, lambda, 1);
    const double f2 = erlang_log_pdf(t, 2.0 * lambda, 2);
    const double w = std::fmax(0.0, std::fmin(1.0, omega));
    if (w >= 1.0) return f1;
    if (w <= 0.0) return f2;
    return log_sum_exp(std::log(w) + f1, std::log1p(-w) + f2);
  }
  if (n <= 1) return std::log(lambda) - lambda * t;
  if (n == 2) return 2.0 * std::log(lambda) + std::log(t) - lambda * t;

  // General n: log(lambda^n * t^(n-1) / (n-1)! * exp(-lambda*t))
  return n * std::log(lambda) + (n - 1.0) * std::log(t) -
         std::lgamma(static_cast<double>(n)) - lambda * t;
}
// CDF of heat kernel N(mean, t) at x
inline double Gstar_CDF(double var, double mean, double x, bool log_p = false) {
  if (!(var > 0.0) || !emc2_isfinite(var)) {
    if (log_p) {
      return (x < mean) ? R_NegInf : 0.0;
    }
    return (x < mean) ? 0.0 : 1.0;
  }

  return gaussian_cdf(x, mean, var, log_p);
}

inline double log_pnorm_diff(double z_lo, double z_hi) {
  // log(Phi(z_hi) - Phi(z_lo)), assuming z_hi >= z_lo
  if (z_hi <= z_lo) return R_NegInf;

  // Both on left side: lower-tail difference is stable.
  if (z_hi <= 0.0) {
    const double log_hi = pnorm_std(z_hi, true, true);
    const double log_lo = pnorm_std(z_lo, true, true);
    return log_diff_exp(log_hi, log_lo);
  }

  // Both on right side: upper-tail difference is stable.
  if (z_lo >= 0.0) {
    const double log_lo_upper = pnorm_std(z_lo, false, true);
    const double log_hi_upper = pnorm_std(z_hi, false, true);
    return log_diff_exp(log_lo_upper, log_hi_upper);
  }

  // Direct probability subtraction is stable when the interval crosses zero.
  const double p_hi = pnorm_std(z_hi, true, false);
  const double p_lo = pnorm_std(z_lo, true, false);
  const double diff = p_hi - p_lo;

  return diff > 0.0 ? std::log(diff) : R_NegInf;
}

// Definite integral ∫_{x_lo}^{x_hi} N(mean, t) dx
inline double Gstar_Integral(double var, double mean,
                             double x_lo, double x_hi,
                             bool log_p = false) {
  if (var <= 0.0 ) {
    const double p = (mean > x_lo && mean <= x_hi) ? 1.0 : 0.0;
    return log_p ? std::log(p) : p;
  }

  if (x_hi <= x_lo) return log_p ? R_NegInf : 0.0;

  if (!log_p) {
    return Gstar_CDF(var, mean, x_hi) - Gstar_CDF(var, mean, x_lo);
  }

  const double sd = std::sqrt(var);
  const double z_lo = (x_lo - mean) / sd;
  const double z_hi = (x_hi - mean) / sd;

  return log_pnorm_diff(z_lo, z_hi);
}

/**
 * @brief Computes ∫[x_lo, x_hi] Φ(ax+c) dx
 * Where Φ is the standard Normal CDF.
 */
inline double integrate_gaussian_cdf(double a, double c, double x_lo, double x_hi) {
  if (std::abs(a) < FPM_EPSILON) {
    // Fallback: integral of a constant
    return Gstar_CDF(1.0, 0.0, c) * (x_hi - x_lo);
  }
  
  // Indefinite integral is (1/a) * [u*Φ(u) + φ(u)], where u = ax+c
  auto eval_integral = [&](double x){
    const double u = a * x + c;
    // Gstar_CDF is Φ(u), Gstar is φ(u) (with var=1)
    return (u * Gstar_CDF(1.0, 0.0, u) + Gstar(1.0, u)) / a;
  };
  
  return eval_integral(x_hi) - eval_integral(x_lo);
}

/**
 * @brief Computes ∫[x_lo, x_hi] exp(kx) * Φ(ax+c) dx
 */
inline double integrate_exp_times_normal_cdf(double k, double a, double c, double x_lo, double x_hi) {
  if (std::abs(k) < FPM_EPSILON) { // Fallback to non-exp integral
    return integrate_gaussian_cdf(a, c, x_lo, x_hi);
  }
  if (std::abs(a) < FPM_EPSILON) { // Fallback: integral of exp(kx) * C
    const double C = Gstar_CDF(1.0, 0.0, c);
    return C * (std::exp(k * x_hi) - std::exp(k * x_lo)) / k;
  }
  
  // Standard identity for this integral
  const double exp_factor = std::exp(k * k / (2.0 * a * a) - k * c / a);
  
  auto eval_integral = [&](double x){
    const double u = a * x + c;
    const double u_shifted = a * x + c - k / a;
    
    const double term1 = (std::exp(k * x) / k) * Gstar_CDF(1.0, 0.0, u);
    const double term2 = (1.0 / k) * exp_factor * Gstar_CDF(1.0, 0.0, u_shifted);
    
    return term1 - term2;
  };
  
  return eval_integral(x_hi) - eval_integral(x_lo);
}

inline double log_h_gaussian_cdf_integral(double u) {
  // h(u) = u * Phi(u) + phi(u)
  // Used by ∫ Phi(ax+c) dx.
  const double log_phi = Gstar(1.0, u, true);
  const double log_Phi = pnorm_std(u, true, true);

  if (u == 0.0) return log_phi;

  if (u > 0.0) {
    return log_sum_exp(std::log(u) + log_Phi, log_phi);
  }

  // u < 0: h(u) = phi(u) - |u| Phi(u)
  const double log_sub = std::log(-u) + log_Phi;
  const double ratio_log = log_sub - log_phi;

  if (ratio_log >= 0.0) {
    // Extreme left-tail fallback:
    // h(u) ≈ phi(u) / u^2
    return log_phi - 2.0 * std::log(-u);
  }

  return log_phi + log1m_exp(ratio_log);
}

inline double log_integrate_gaussian_cdf(double a, double c,
                                         double x_lo, double x_hi) {
  if (x_hi <= x_lo) return R_NegInf;

  if (std::abs(a) < FPM_EPSILON) {
    const double log_C = pnorm_std(c, true, true);
    return log_C + std::log(x_hi - x_lo);
  }

  auto eval_integral_log = [&](double x) {
    const double u = a * x + c;
    return make_signed_log(
      log_h_gaussian_cdf_integral(u) - std::log(std::abs(a)),
      a > 0.0 ? 1 : -1
    );
  };

  signed_log hi = eval_integral_log(x_hi);
  signed_log lo = eval_integral_log(x_lo);
  signed_log out = signed_log_sub(hi, lo);

  if (out.sign <= 0) return R_NegInf;
  return out.log_abs;
}

/**
 * @brief Computes log ∫[x_lo, x_hi] exp(kx) * Φ(ax+c) dx
 */
inline double log_integrate_exp_times_normal_cdf(double k, double a, double c,
                                                 double x_lo, double x_hi) {
  if (x_hi <= x_lo) return R_NegInf;

  if (std::abs(k) < FPM_EPSILON) {
    return log_integrate_gaussian_cdf(a, c, x_lo, x_hi);
  }

  if (std::abs(a) < FPM_EPSILON) {
    const double log_C = pnorm_std(c, true, true);

    if (k > 0.0) {
      return log_C +
        log_diff_exp(k * x_hi, k * x_lo) -
        std::log(k);
    } else {
      return log_C +
        log_diff_exp(k * x_lo, k * x_hi) -
        std::log(-k);
    }
  }

  const double log_abs_k = std::log(std::abs(k));
  const int sign_k = k > 0.0 ? 1 : -1;

  const double log_exp_factor =
    k * k / (2.0 * a * a) - k * c / a;

  auto eval_integral = [&](double x) {
    const double u = a * x + c;
    const double u_shifted = a * x + c - k / a;

    // term1 = exp(k*x) / k * Phi(u)
    signed_log term1 = make_signed_log(
      k * x - log_abs_k + pnorm_std(u, true, true),
      sign_k
    );

    // term2 = exp_factor / k * Phi(u_shifted)
    // eval = term1 - term2
    signed_log term2 = make_signed_log(
      log_exp_factor - log_abs_k + pnorm_std(u_shifted, true, true),
      sign_k
    );

    return signed_log_sub(term1, term2);
  };

  signed_log hi = eval_integral(x_hi);
  signed_log lo = eval_integral(x_lo);
  signed_log out = signed_log_sub(hi, lo);

  if (out.sign <= 0 || out.log_abs == R_NegInf) {
    // Rare numerical fallback. This keeps the sampler from detonating.
    const double val = integrate_exp_times_normal_cdf(k, a, c, x_lo, x_hi);
    return val > 0.0 ? std::log(val) : R_NegInf;
  }

  return out.log_abs;
}

// Bivariate Normal CDF Functions
struct LegendreHalfRule {
  std::vector<double> x;
  std::vector<double> w;
  std::vector<double> xp;
  std::vector<double> xn;
  std::vector<double> w_div_4pi;
  int pair_count = 0;
  bool has_center = false;
};

inline LegendreHalfRule make_legendre_half_rule(int order) {
  const Rcpp::List q = gauss_quad(order, "legendre");
  const Rcpp::NumericVector nodes = q["nodes"];
  const Rcpp::NumericVector weights = q["weights"];

  LegendreHalfRule rule;
  rule.pair_count = order / 2;
  rule.has_center = (order % 2) == 1;
  rule.x.reserve(rule.pair_count + (rule.has_center ? 1 : 0));
  rule.w.reserve(rule.pair_count + (rule.has_center ? 1 : 0));
  rule.xp.reserve(rule.pair_count + (rule.has_center ? 1 : 0));
  rule.xn.reserve(rule.pair_count + (rule.has_center ? 1 : 0));
  rule.w_div_4pi.reserve(rule.pair_count + (rule.has_center ? 1 : 0));

  // Keep only nonnegative abscissas, ordered from largest to smallest;
  // for odd rules append center node (x=0) last.
  constexpr double tol = 1e-14;
  for (int i = nodes.size() - 1; i >= 0; --i) {
    if (nodes[i] > tol) {
      const double xi = nodes[i];
      const double wi = weights[i];
      rule.x.push_back(xi);
      rule.w.push_back(wi);
      rule.xp.push_back(0.5 * xi + 0.5);
      rule.xn.push_back(-0.5 * xi + 0.5);
      rule.w_div_4pi.push_back(wi / (4.0 * M_PI));
    }
  }

  if (rule.has_center) {
    for (int i = 0; i < nodes.size(); ++i) {
      if (std::abs(nodes[i]) <= tol) {
        const double xi = 0.0;
        const double wi = weights[i];
        rule.x.push_back(xi);
        rule.w.push_back(wi);
        rule.xp.push_back(0.5);
        rule.xn.push_back(0.5);
        rule.w_div_4pi.push_back(wi / (4.0 * M_PI));
        break;
      }
    }
  }

  return rule;
}

inline const LegendreHalfRule& get_legendre_half_rule(int order) {
  switch (order) {
  case 5: {
    static const LegendreHalfRule rule = make_legendre_half_rule(5);
    return rule;
  }
  case 6: {
    static const LegendreHalfRule rule = make_legendre_half_rule(6);
    return rule;
  }
  case 8: {
    static const LegendreHalfRule rule = make_legendre_half_rule(8);
    return rule;
  }
  case 12: {
    static const LegendreHalfRule rule = make_legendre_half_rule(12);
    return rule;
  }
  case 16: {
    static const LegendreHalfRule rule = make_legendre_half_rule(16);
    return rule;
  }
  case 20: {
    static const LegendreHalfRule rule = make_legendre_half_rule(20);
    return rule;
  }
  case 24: {
    static const LegendreHalfRule rule = make_legendre_half_rule(24);
    return rule;
  }
  default:
    Rcpp::stop("Unsupported Legendre rule order: %d", order);
  }
  static const LegendreHalfRule unreachable{};
  return unreachable;
}

// The accurate bivariate-normal route is called repeatedly with the same rho
// while one correlated BAwL parameter cell is evaluated.  TVPACK's formula
// spends a noticeable fraction of that time rebuilding the rho-only part of
// its Legendre rule (asin, sin, and 1 - sin^2).  Keep the last rule in the
// calling thread, just as the RDMSWTN BVN seed keeps its per-law table.  The
// cache is deliberately one-entry: a likelihood call normally has one rho,
// while a small cache would add linear scans to a much hotter path.
struct BvnTvpackCache {
  double rho = R_NaN;
  int order = 0;
  int n = 0;
  int n_vec = 0;  // n rounded up to a multiple of 4; lanes past n are inert
  double scale = 0.0;
  std::array<double, 24> sn{};
  std::array<double, 24> den{};
  std::array<double, 24> weight{};

  void prepare(double r, int requested_order) {
    if (rho == r && order == requested_order) return;
    const LegendreHalfRule& rule = get_legendre_half_rule(requested_order);
    const double asr = std::asin(r);
    const double half = 0.5 * asr;
    n = static_cast<int>(rule.x.size()) * 2;
    n_vec = (n + 3) & ~3;
    for (int ix = n; ix < n_vec; ++ix) {
      sn[static_cast<size_t>(ix)] = 0.0;
      den[static_cast<size_t>(ix)] = 1.0;
      weight[static_cast<size_t>(ix)] = 0.0;
    }
    scale = asr / fourPI;
    for (int ix = 0; ix < static_cast<int>(rule.x.size()); ++ix) {
      const double xp = half * (1.0 + rule.x[ix]);
      const double xn = half * (1.0 - rule.x[ix]);
      const double sp = std::sin(xp);
      const double snv = std::sin(xn);
      sn[static_cast<size_t>(2 * ix)] = sp;
      sn[static_cast<size_t>(2 * ix + 1)] = snv;
      den[static_cast<size_t>(2 * ix)] = std::fma(-sp, sp, 1.0);
      den[static_cast<size_t>(2 * ix + 1)] = std::fma(-snv, snv, 1.0);
      weight[static_cast<size_t>(2 * ix)] = rule.w[ix];
      weight[static_cast<size_t>(2 * ix + 1)] = rule.w[ix];
    }
    rho = r;
    order = requested_order;
  }
};

inline BvnTvpackCache& bvn_tvpack_cache() {
  static thread_local BvnTvpackCache cache;
  return cache;
}

// For |rho| >= .925 the high-correlation TVPACK branch has a second
// rho-only setup: the transformed Legendre abscissas and their square-root
// denominators.  Cache those too.  The x-dependent exponentials remain in the
// hot loop, so this does not alter the numerical formula or its tail behaviour.
struct BvnTvpackHighCache {
  double rho = R_NaN;
  double half_sqrt = 0.0;
  int n = 0;
  std::array<double, 24> xs{};
  std::array<double, 24> rs{};
  std::array<double, 24> weight{};

  void prepare(double r) {
    if (rho == r) return;
    const LegendreHalfRule& rule = get_legendre_half_rule(20);
    half_sqrt = 0.5 * std::sqrt(1.0 - r * r);
    n = static_cast<int>(rule.x.size()) * 2;
    for (int ix = 0; ix < static_cast<int>(rule.x.size()); ++ix) {
      const double tp = half_sqrt * (1.0 + rule.x[ix]);
      const double tn = half_sqrt * (1.0 - rule.x[ix]);
      const double xp = tp * tp;
      const double xn = tn * tn;
      xs[static_cast<size_t>(2 * ix)] = xp;
      xs[static_cast<size_t>(2 * ix + 1)] = xn;
      rs[static_cast<size_t>(2 * ix)] = std::sqrt(1.0 - xp);
      rs[static_cast<size_t>(2 * ix + 1)] = std::sqrt(1.0 - xn);
      weight[static_cast<size_t>(2 * ix)] = half_sqrt * rule.w[ix];
      weight[static_cast<size_t>(2 * ix + 1)] = half_sqrt * rule.w[ix];
    }
    rho = r;
  }
};

inline BvnTvpackHighCache& bvn_tvpack_high_cache() {
  static thread_local BvnTvpackHighCache cache;
  return cache;
}

// Drezner's fast approximation has the same opportunity on its hot path.
// The five-point rule is fixed, but its rho-scaled abscissas, square-root
// denominators, and weights are not.  Trends on v make the CDF call count
// large while rho is usually a cell constant, so retain this small table per
// thread and per rho.
struct BvnDreznerCache {
  double rho = R_NaN;
  bool high = false;
  int n = 0;
  double rho_var = 0.0;
  double r3 = 0.0;
  std::array<double, 8> node{};
  std::array<double, 8> variance{};
  std::array<double, 8> inv_sqrt{};
  std::array<double, 8> denom{};
  std::array<double, 8> weight{};

  void prepare(double r) {
    if (rho == r) return;
    const LegendreHalfRule& rule = get_legendre_half_rule(5);
    high = std::fabs(r) >= 0.7;
    rho_var = 1.0 - r * r;
    n = rule.pair_count * 2 + (rule.has_center ? 1 : 0);
    if (high) {
      r3 = std::sqrt(rho_var);
      for (int ix = 0; ix < rule.pair_count; ++ix) {
        const double rp = r3 * rule.xp[ix];
        const double rn = r3 * rule.xn[ix];
        node[static_cast<size_t>(2 * ix)] = rp * rp;
        node[static_cast<size_t>(2 * ix + 1)] = rn * rn;
        variance[static_cast<size_t>(2 * ix)] =
          1.0 - node[static_cast<size_t>(2 * ix)];
        variance[static_cast<size_t>(2 * ix + 1)] =
          1.0 - node[static_cast<size_t>(2 * ix + 1)];
        denom[static_cast<size_t>(2 * ix)] =
          std::sqrt(variance[static_cast<size_t>(2 * ix)]);
        denom[static_cast<size_t>(2 * ix + 1)] =
          std::sqrt(variance[static_cast<size_t>(2 * ix + 1)]);
        inv_sqrt[static_cast<size_t>(2 * ix)] =
          1.0 / denom[static_cast<size_t>(2 * ix)];
        inv_sqrt[static_cast<size_t>(2 * ix + 1)] =
          1.0 / denom[static_cast<size_t>(2 * ix + 1)];
        weight[static_cast<size_t>(2 * ix)] = rule.w_div_4pi[ix];
        weight[static_cast<size_t>(2 * ix + 1)] = rule.w_div_4pi[ix];
      }
      if (rule.has_center) {
        const int ix = rule.pair_count;
        const double rp = r3 * rule.xp[ix];
        const size_t k = static_cast<size_t>(2 * ix);
        node[k] = rp * rp;
        variance[k] = 1.0 - node[k];
        denom[k] = std::sqrt(variance[k]);
        inv_sqrt[k] = 1.0 / denom[k];
        weight[k] = rule.w_div_4pi[ix];
      }
    } else {
      for (int ix = 0; ix < rule.pair_count; ++ix) {
        const double rp = r * rule.xp[ix];
        const double rn = r * rule.xn[ix];
        node[static_cast<size_t>(2 * ix)] = rp;
        node[static_cast<size_t>(2 * ix + 1)] = rn;
        variance[static_cast<size_t>(2 * ix)] = 1.0 - rp * rp;
        variance[static_cast<size_t>(2 * ix + 1)] = 1.0 - rn * rn;
        denom[static_cast<size_t>(2 * ix)] =
          std::sqrt(variance[static_cast<size_t>(2 * ix)]);
        denom[static_cast<size_t>(2 * ix + 1)] =
          std::sqrt(variance[static_cast<size_t>(2 * ix + 1)]);
        inv_sqrt[static_cast<size_t>(2 * ix)] =
          1.0 / denom[static_cast<size_t>(2 * ix)];
        inv_sqrt[static_cast<size_t>(2 * ix + 1)] =
          1.0 / denom[static_cast<size_t>(2 * ix + 1)];
        weight[static_cast<size_t>(2 * ix)] = rule.w_div_4pi[ix];
        weight[static_cast<size_t>(2 * ix + 1)] = rule.w_div_4pi[ix];
      }
      if (rule.has_center) {
        const int ix = rule.pair_count;
        const double rp = r * rule.xp[ix];
        const size_t k = static_cast<size_t>(2 * ix);
        node[k] = rp;
        variance[k] = 1.0 - rp * rp;
        denom[k] = std::sqrt(variance[k]);
        inv_sqrt[k] = 1.0 / denom[k];
        weight[k] = rule.w_div_4pi[ix];
      }
    }
    // Inert lanes n..7 for the fixed-width vector pass in norm_ucdf_2d_fast:
    // finite in both branches' formulas, zero weight, never summed.
    for (int ix = n; ix < 8; ++ix) {
      node[static_cast<size_t>(ix)] = 0.5;
      variance[static_cast<size_t>(ix)] = 0.75;
      denom[static_cast<size_t>(ix)] = std::sqrt(0.75);
      inv_sqrt[static_cast<size_t>(ix)] = 1.0 / std::sqrt(0.75);
      weight[static_cast<size_t>(ix)] = 0.0;
    }
    rho = r;
  }
};

inline BvnDreznerCache& bvn_drezner_cache() {
  static thread_local BvnDreznerCache cache;
  return cache;
}

// Marginal normal CDFs that norm_ucdf_2d_fast would otherwise evaluate itself,
// supplied by a caller that visits the same boundary at several corners (the
// correlated-BAwL boundary grid).  In the upper-orthant arguments (x1, x2):
//   g1 = gaussian_cdf(-x1), g2 = gaussian_cdf(-x2), h2 = gaussian_cdf(x2),
// the last used only on the high-|rho|, rho < 0 branch.  A flag that is false
// means "not supplied": the value is computed exactly as before, so a supplied
// value must be the bit-identical pnorm_std result of the same argument.
struct BvnFastMarginals {
  double g1 = 0.0, g2 = 0.0, h2 = 0.0;
  bool has_g1 = false, has_g2 = false, has_h2 = false;
};

/* Bivariate normal CDF approximation using 5-point Gauss-Legendre quadrature
 * with p-split refinement and tail parameter cutoffs.  `cache` must already be
 * prepared for rho. */
inline double norm_ucdf_2d_fast_core(double x1, double x2, double rho,
                                     const BvnDreznerCache& cache,
                                     const BvnFastMarginals& m)
{
  auto g1 = [&]() { return m.has_g1 ? m.g1 : gaussian_cdf(-x1); };
  double x12 = 0.5 * (x1*x1 + x2*x2);

  double out = 0;
  double x3;
  // Node exponentials are computed for all eight lanes in one vector pass
  // (lanes past cache.n are inert) and summed in the original node order.
  double e_a[8], e_b[8];
  if (cache.high) {
    const double r2 = cache.rho_var;
    const double r3 = cache.r3;
    const double x2_in = x2;
    if (rho < 0) {
      x2 = -x2;
    }
    x3 = x1*x2;
    double x7 = std::exp(-0.5 * x3);
    if (r2) {
      double x6 = std::fabs(x1 - x2);
      double x5 = 0.5 * x6*x6;
      x6 /= r3;
      double aa = 0.5 - 0.125*x3;
      double ab = 3. - 2.*aa*x5;
      out = inv3sqrt2pi * (
        x6 * ab * gaussian_cdf(-x6) -
          std::exp(-x5/r2) * std::fma(aa, r2, ab) * inv_sqrt2_pi
      );
#pragma omp simd
      for (int ix = 0; ix < 8; ++ix) {
        const double xs = cache.node[static_cast<size_t>(ix)];
        const double rs = cache.denom[static_cast<size_t>(ix)];
        e_a[ix] = EMC2_BVN_EXP(-x5/xs);
        e_b[ix] = EMC2_BVN_EXP(-x3/(1. + rs));
      }
      for (int ix = 0; ix < cache.n; ++ix) {
        const double xs = cache.node[static_cast<size_t>(ix)];
        const double inv_s = cache.inv_sqrt[static_cast<size_t>(ix)];
        out -= cache.weight[static_cast<size_t>(ix)] *
          e_a[ix] * (e_b[ix] * inv_s / x7 - 1. - aa * xs);
      }
    }
    if (rho > 0) {
      // gaussian_cdf(-max(x1, x2)) is whichever marginal has the larger x.
      const double gmax = (x1 >= x2)
        ? g1() : (m.has_g2 ? m.g2 : gaussian_cdf(-x2));
      out = std::fma(out, r3*x7, gmax);
    }
    else {
      // x2 was negated above, so gaussian_cdf(-x2) here is gaussian_cdf(x2_in).
      const double h2 = m.has_h2 ? m.h2 : gaussian_cdf(x2_in);
      out = std::fmax(0., g1() - h2) - out*r3*x7;
    }
    return out;
  }
  else {
    x3 = x1*x2;
#pragma omp simd
    for (int ix = 0; ix < 8; ++ix) {
      const double r1 = cache.node[static_cast<size_t>(ix)];
      const double rr2 = cache.variance[static_cast<size_t>(ix)];
      e_a[ix] = EMC2_BVN_EXP((r1*x3 - x12) / rr2);
    }
    for (int ix = 0; ix < cache.n; ++ix) {
      out += cache.weight[static_cast<size_t>(ix)] * e_a[ix] *
        cache.inv_sqrt[static_cast<size_t>(ix)];
    }
    const double g2 = m.has_g2 ? m.g2 : gaussian_cdf(-x2);
    return std::fma(out, rho, g1() * g2);
  }
}

inline double norm_ucdf_2d_fast(double x1, double x2, double rho)
{
  BvnDreznerCache& cache = bvn_drezner_cache();
  cache.prepare(rho);
  return norm_ucdf_2d_fast_core(x1, x2, rho, cache, BvnFastMarginals());
}

inline double norm_cdf_2d_fast(double x1, double x2, double rho)
{
  return norm_ucdf_2d_fast(-x1, -x2, rho);
}

/* Fast bivariate-normal approximation. */
constexpr const static double c1 = -1.0950081470333;
constexpr const static double c2 = -0.75651138383854;
inline double norm_cdf_2d_vfast(double x1, double x2, double rho)
{
  if (std::fabs(rho) <= std::numeric_limits<double>::epsilon()) {
    return gaussian_cdf(x1) * gaussian_cdf(x2);
  }
  
  double denom = std::sqrt(1 - rho * rho);
  double a = -rho / denom;
  double b = x1 / denom;
  double aq_plus_b = a*x2 + b;
  
  if (a > 0) {
    if (aq_plus_b >= 0) {
      double aa = a * a;
      double a_sq_c1 = aa*c1;
      double a_sq_c2 = aa*c2;
      double sqrt2b = sqrt2*b;
      double sqrt2x2 = sqrt2*x2;
      double sqrt_recpr_a_sq_c2 = std::sqrt(1. - a_sq_c2);
      double twicea_sqrt_recpr_a_sq_c2 = 2.*a*sqrt_recpr_a_sq_c2;
      double temp = 1. / (4. * sqrt_recpr_a_sq_c2);
      double t1 = a_sq_c1*c1 + 2.*b*b*c2;
      double t2 = 2.*sqrt2b*c1;
      double t3 = 4. - 4.*a_sq_c2;
      
      return
      0.5 * (std::erf(x2 / sqrt2) + std::erf(b / (sqrt2*a))) +
        temp
        * std::exp((t1 - t2) / t3)
        * (1. - std::erf((sqrt2b - a_sq_c1) / twicea_sqrt_recpr_a_sq_c2)) -
        temp
        * std::exp((t1 + t2) / t3)
        * (
            std::erf((sqrt2x2 - sqrt2x2*a_sq_c2 - sqrt2b*a*c2 - a*c1) / (2.*sqrt_recpr_a_sq_c2)) +
            std::erf((a_sq_c1 + sqrt2b) / twicea_sqrt_recpr_a_sq_c2)
        );
      
    }
    else {
      double sqrt2b = sqrt2*b;
      double sqrt2x2 = sqrt2*x2;
      double a_sq_c2 = a*a*c2;
      double recpr_a_sq_c2 = 1. - a_sq_c2;
      double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
      double a_c1 = a*c1;
      
      return
      (1. / (4. * sqrt_recpr_a_sq_c2)) *
        std::exp((a_c1*a_c1 - 2.*sqrt2b*c1 + 2*b*b*c2) / (4.*recpr_a_sq_c2)) *
        (1. + std::erf((sqrt2x2 - sqrt2x2*a_sq_c2 - sqrt2b*a*c2 + a_c1) / (2.*sqrt_recpr_a_sq_c2)));
    }
  }
  else {
    if (aq_plus_b >= 0) {
      double sqrt2b = sqrt2*b;
      double a_sq_c2 = a*a*c2;
      double recpr_a_sq_c2 = 1. - a_sq_c2;
      double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
      double a_c1 = a*c1;
      double sqrt2_x2 = sqrt2*x2;
      
      return
      0.5 + 0.5 * std::erf(x2 / sqrt2) -
        (1. / (4. * sqrt_recpr_a_sq_c2)) *
        std::exp((a_c1*a_c1 + 2.*sqrt2b*c1 + 2.*b*b*c2) / (4.*recpr_a_sq_c2)) *
        (1. + std::erf((sqrt2_x2 - sqrt2_x2*a_sq_c2 - sqrt2b*a*c2 - a_c1) / (2.*sqrt_recpr_a_sq_c2)));
    }
    else {
      double sqrt2a = sqrt2*a;
      double sqrt2b = sqrt2*b;
      double a_sq_c2 = a*a*c2;
      double recpr_a_sq_c2 = 1. - a_sq_c2;
      double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
      double a_c1 = a*c1;
      double temp = 1. / (4. * sqrt_recpr_a_sq_c2);
      double t1 = a_c1*a_c1 + 2.*b*b*c2;
      double t2 = 2.*sqrt2b*c1;
      double t3 = 4.*recpr_a_sq_c2;
      double sqrt2_x2 = sqrt2*x2;
      
      return
      0.5 - 0.5 * std::erf(b / sqrt2a) -
        temp
        * std::exp((t1 + t2) / t3)
        * (1. - std::erf((sqrt2b + a*a_c1) / (2.*a*sqrt_recpr_a_sq_c2))) +
        temp
        * std::exp((t1 - t2) / t3)
        * (
            std::erf((sqrt2_x2 - sqrt2_x2*a_sq_c2 - sqrt2b*a*c2 + a_c1) / (2.*sqrt_recpr_a_sq_c2)) +
            std::erf((sqrt2b - a*a_c1) / (2.*a*sqrt_recpr_a_sq_c2))
        );
    }
  }
} 

inline double norm_ucdf_2d(double x1, double x2, double rho)
{
  const int low_rho_order = 6;
  const int mid_rho_order = 12;
  const int high_rho_order = 20;
  const LegendreHalfRule& low_rule = get_legendre_half_rule(low_rho_order);
  const LegendreHalfRule& mid_rule = get_legendre_half_rule(mid_rho_order);
  const LegendreHalfRule& high_rule = get_legendre_half_rule(high_rho_order);

  double abs_rho = std::fabs(rho);
  
  double out = 0;
  double hk;
  if (abs_rho < 0.925) {
    if (abs_rho > std::numeric_limits<double>::epsilon()) {
      hk = x1 * x2;
      double hs = 0.5 * (x1*x1 + x2*x2);
      const LegendreHalfRule& active_rule =
        (abs_rho < 0.3) ? low_rule : ((abs_rho < 0.5) ? mid_rule : high_rule);
      BvnTvpackCache& cache = bvn_tvpack_cache();
      cache.prepare(rho, active_rule.x.size() == low_rule.x.size() ? 6 :
                    (active_rule.x.size() == mid_rule.x.size() ? 12 : 20));
      // Exponentials in one vector pass; the weighted sum keeps its
      // sequential order.
      double node_exp[24];
#pragma omp simd
      for (int ix = 0; ix < cache.n_vec; ++ix) {
        node_exp[ix] = EMC2_BVN_EXP(
          std::fma(cache.sn[static_cast<size_t>(ix)], hk, -hs) /
            cache.den[static_cast<size_t>(ix)]);
      }
      for (int ix = 0; ix < cache.n; ++ix) {
        out += cache.weight[static_cast<size_t>(ix)] * node_exp[ix];
      }
      out *= cache.scale;
    }
    out = std::fma(gaussian_cdf(-x1), gaussian_cdf(-x2), out);
  }
  else {
    
    if (rho < 0) {
      x2 = -x2;
    }
    if (abs_rho < 1) {
      hk = x1 * x2;
      double as = std::fma(-rho, rho, 1.);
      double a = std::sqrt(as);
      double b;
      double bs = (x1 - x2) * (x1 - x2);
      double c = std::fma(-0.125, hk, 0.5);
      double d = std::fma(-0.0625, hk, 0.75);
      double asr = -0.5 * (hk + bs / as);
      double rfdbs = std::fma(-d, bs, 5.)*(1./15.);
      if (asr > -100.) {
        out = a * std::exp(asr) * (1. - c * (bs - as) * rfdbs + 0.2*c*d*as*as);
      }
      if (hk > -100.) {
        b = std::sqrt(bs);
        out -= std::exp(-0.5 * hk) * sqrt_twoPI * gaussian_cdf(-b / a) * b * (1. - c*bs*rfdbs);
      }
      BvnTvpackHighCache& cache = bvn_tvpack_high_cache();
      cache.prepare(rho);
      // Both exponentials for every node in one vector pass (n = 20); the
      // asr > -100 screen and the sequential sum are applied afterwards.  A
      // screened-out node never enters the sum, so its second exponential
      // (whose argument is bounded by |hk|/2 there) is simply discarded.
      double node_asr[24], node_e1[24], node_e2[24];
#pragma omp simd
      for (int ix = 0; ix < cache.n; ++ix) {
        const double xs = cache.xs[static_cast<size_t>(ix)];
        const double rs = cache.rs[static_cast<size_t>(ix)];
        node_asr[ix] = -0.5 * (hk + bs / xs);
        node_e1[ix] = EMC2_BVN_EXP(node_asr[ix]);
        node_e2[ix] = EMC2_BVN_EXP(-hk*xs/(2.*(1.+rs)*(1.+rs)));
      }
      for (int ix = 0; ix < cache.n; ++ix) {
        const double xs = cache.xs[static_cast<size_t>(ix)];
        const double rs = cache.rs[static_cast<size_t>(ix)];
        if (node_asr[ix] > -100.) {
          out += cache.weight[static_cast<size_t>(ix)] * node_e1[ix] *
            (node_e2[ix]/rs - (1. + c*xs*std::fma(d, xs, 1.)));
        }
      }
      out *= minus_inv_twoPI;
    }
    if (rho > 0) {
      out += gaussian_cdf(-std::fmax(x1, x2));
    }
    else {
      out = -out;
      if (x2 > x1) {
        if (x1 < 0) {
          out += gaussian_cdf(x2) - gaussian_cdf(x1);
        }
        else {
          out += gaussian_cdf(-x1) - gaussian_cdf(-x2);
        }
      }
    }
  }
  return out;
}

inline double norm_cdf_2d(double x1, double x2, double rho)
{
  return norm_ucdf_2d(-x1, -x2, rho);
}

// Hybrid dispatch is based on probability magnitude rather than |x| or |y|.
// Since Phi2(x,y,rho) <= min(Phi(x), Phi(y)), the lower-tail screen routes
// sufficiently small probabilities to the stable tvpack path.
// Small negative approximation outputs are clamped to zero before log().
constexpr double EMC2_BVN_DREZNER_MIN_P = 1e-3;
constexpr double EMC2_BVN_MIN_Z = -3.090232306167813;  // qnorm(1e-3)

inline double norm_cdf_2d_hybrid(double x, double y, double rho) {
  double out;
  if (std::fabs(rho) > 0.9999 ||
      (x < EMC2_BVN_MIN_Z) || (y < EMC2_BVN_MIN_Z)) {
    out = norm_cdf_2d(x, y, rho);
  } else {
    out = norm_cdf_2d_fast(x, y, rho);
    if (!(out >= EMC2_BVN_DREZNER_MIN_P)) out = norm_cdf_2d(x, y, rho);
  }
  return out > 0.0 ? out : 0.0;
}

// The R-facing pbvn_* wrappers need external linkage so RcppExports.cpp can link
// against them, so they are defined out-of-line in gaussian.cpp (same pattern as
// wald_functions.cpp / cdf_fncs.cpp) rather than inline here.
double pbvn_tsay(double h, double k, double rho);
double pbvn_tvpack(double h, double k, double rho);
double pbvn_drezner(double h, double k, double rho);
double pbvn_hybrid(double h, double k, double rho);

#endif

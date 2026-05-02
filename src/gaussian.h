#ifndef gaussian_h
#define gaussian_h

#include <RcppArmadillo.h>
#include <algorithm>
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

// Useful Distribution Functions
inline double gaussian_pdf(double x, double mean = 0.0, double var = 1.0, bool log_p = false) {
  if (!(var > 0.0) || !std::isfinite(var)) return log_p ? R_NegInf : 0.0;;
  
  const double z = x - mean;

  if (log_p) {
    return -0.5 * (log_twoPI + std::log(var) + z * z / var);
  }

  return inv_sqrt2_pi / std::sqrt(var) * std::exp(-0.5 * z * z / var);
}

inline double gaussian_cdf(double x, double mean = 0.0, double var = 1.0, bool log_p=false) {
  if (!(var > 0.0) || !std::isfinite(var)) {
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
inline double erlang_log_surv(double t, double lambda, int n) {
  if (lambda <= 0.0 || t <= 0.0) return 0.0;
  if (t == R_PosInf) return R_NegInf;
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
inline double erlang_log_pdf(double t, double lambda, int n) {
  if (lambda <= 0.0 || t < 0.0) return R_NegInf;
  if (t == R_PosInf) return R_NegInf;
  if (n <= 1) return std::log(lambda) - lambda * t;
  if (n == 2) return 2.0 * std::log(lambda) + std::log(t) - lambda * t;
  
  // General n: log(lambda^n * t^(n-1) / (n-1)! * exp(-lambda*t))
  return n * std::log(lambda) + (n - 1.0) * std::log(t) - 
         std::lgamma(static_cast<double>(n)) - lambda * t;
}

// CDF of heat kernel N(mean, t) at x
inline double Gstar_CDF(double var, double mean, double x, bool log_p = false) {
  if (!(var > 0.0) || !std::isfinite(var)) {
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

  // Interval crosses zero, direct probability difference is usually fine.
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
// Functions distributed in this file taken from https://github.com/david-cortes/approxcdf
// Copyright 2022 David Cortes under BSD-3 license
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

/* Bivariate normal CDF.
 Algorithm:  Drezner (1978) 5-point GL with the p-split
 refined by Drezner & Wesolowsky (1990);
 parameter cut-off as in West (2004).
 Expected accuracy ~ 1e-6.
 */
double norm_ucdf_2d_fast(double x1, double x2, double rho)
{
  const int drezner_order = 5;
  const LegendreHalfRule& drezner_rule = get_legendre_half_rule(drezner_order);

  double x12 = 0.5 * (x1*x1 + x2*x2);
  
  double out = 0;
  double r1, x3;
  if (std::fabs(rho) >= 0.7) {
    double r2 = 1. - rho*rho;
    double r3 = std::sqrt(r2);
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
      double rr;
      double nr1, nrr, nr2;
      for (int ix = 0; ix < drezner_rule.pair_count; ix++) {
        r1 = r3 * drezner_rule.xp[ix];
        rr = r1*r1;
        r2 = std::sqrt(1. - rr);
        
        nr1 = r3 * drezner_rule.xn[ix];
        nrr = nr1*nr1;
        nr2 = std::sqrt(1. - nrr);
        
        out -= drezner_rule.w_div_4pi[ix] * (
          std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr) +
            std::exp(-x5/nrr) * (std::exp(-x3/(1. + nr2))/nr2/x7 - 1.- aa*nrr)
        );
      }
      if (drezner_rule.has_center) {
        const int center_ix = drezner_rule.pair_count;
        r1 = r3 * drezner_rule.xp[center_ix];
        rr = r1*r1;
        r2 = std::sqrt(1. - rr);
        out -= drezner_rule.w_div_4pi[center_ix] *
          std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr);
      }
    }
    if (rho > 0) {
      out = std::fma(out, r3*x7, gaussian_cdf(-std::fmax(x1, x2)));
    }
    else {
      out = std::fmax(0., gaussian_cdf(-x1) - gaussian_cdf(-x2)) - out*r3*x7;
    }
    return out;
  }
  else {
    x3 = x1*x2;
    double rr2;
    double nr1, nrr2;
    for (int ix = 0; ix < drezner_rule.pair_count; ix++) {
      r1 = rho * drezner_rule.xp[ix];
      rr2 = 1. - r1*r1;
      
      nr1 = rho * drezner_rule.xn[ix];
      nrr2 = 1. - nr1*nr1;
      
      out += drezner_rule.w_div_4pi[ix] * (
        std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2) +
          std::exp((nr1*x3 - x12) / nrr2) / std::sqrt(nrr2)
      );
    }
    if (drezner_rule.has_center) {
      const int center_ix = drezner_rule.pair_count;
      r1 = rho * drezner_rule.xp[center_ix];
      rr2 = 1. - r1*r1;
      out += drezner_rule.w_div_4pi[center_ix] *
        std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2);
    }
    return std::fma(out, rho, gaussian_cdf(-x1) * gaussian_cdf(-x2));
  }
}

double norm_cdf_2d_fast(double x1, double x2, double rho)
{
  return norm_ucdf_2d_fast(-x1, -x2, rho);
}

/* Tsay, Wen-Jen, and Peng-Hsuan Ke.
 "A simple approximation for the bivariate normal integral."
 Communications in Statistics-Simulation and Computation (2021): 1-14. */
constexpr const static double c1 = -1.0950081470333;
constexpr const static double c2 = -0.75651138383854;
double norm_cdf_2d_vfast(double x1, double x2, double rho)
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

[[gnu::flatten]]
double norm_ucdf_2d(double x1, double x2, double rho)
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
      double asr = std::asin(rho);
      double asr_half = 0.5 * asr;
      double sn1, sn2;
      
      const LegendreHalfRule& active_rule =
        (abs_rho < 0.3) ? low_rule : ((abs_rho < 0.5) ? mid_rule : high_rule);

      for (int ix = 0; ix < static_cast<int>(active_rule.x.size()); ix++) {
        sn1 = std::sin(asr_half * (1. + active_rule.x[ix]));
        sn2 = std::sin(asr_half * (1. - active_rule.x[ix]));
        out += active_rule.w[ix] * (
          std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
        );
      }
      out *= asr / fourPI;
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
      a *= 0.5;
      double xs;
      double rs;
      double temp;
      
      for (int ix = 0; ix < static_cast<int>(high_rule.x.size()); ix++) {
        temp = a * (1. + high_rule.x[ix]);
        xs = temp * temp;
        rs = std::sqrt(1. - xs);
        asr = -0.5 * (hk + bs / xs);
        if (asr > -100.) {
          out += a * high_rule.w[ix] * std::exp(asr) *
            (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
        }
        
        temp = a * (1. - high_rule.x[ix]);
        xs = temp * temp;
        rs = std::sqrt(1. - xs);
        asr = -0.5 * (hk + bs / xs);
        if (asr > -100.) {
          out += a * high_rule.w[ix] * std::exp(asr) *
            (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
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

double norm_cdf_2d(double x1, double x2, double rho)
{
  return norm_ucdf_2d(-x1, -x2, rho);
}

// [[Rcpp::export]]
double pbvn_tsay(double h, double k, double rho) {
  return norm_cdf_2d_vfast(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_tvpack(double h, double k, double rho) {
  return norm_cdf_2d(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_drezner(double h, double k, double rho) {
  return norm_cdf_2d_fast(h, k, rho);
}


#endif

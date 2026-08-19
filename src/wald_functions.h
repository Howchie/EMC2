#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include "utility_functions.h"
#include "composite_functions.h"
#include "gl_quad.h"

using namespace Rcpp;

constexpr double L_PI = 1.1447298858494001741434; 
constexpr double LOG_2PI = 1.83787706640934548356; 
constexpr double FAST_NORM_RT2PI = 2.506628274631000502415765284811;
constexpr double FAST_NORM_SPLIT = 7.07106781186547;
constexpr double FAST_NORM_N0 = 220.206867912376;
constexpr double FAST_NORM_N1 = 221.213596169931;
constexpr double FAST_NORM_N2 = 112.079291497871;
constexpr double FAST_NORM_N3 = 33.912866078383;
constexpr double FAST_NORM_N4 = 6.37396220353165;
constexpr double FAST_NORM_N5 = 0.700383064443688;
constexpr double FAST_NORM_N6 = 3.52624965998911e-02;
constexpr double FAST_NORM_M0 = 440.413735824752;
constexpr double FAST_NORM_M1 = 793.826512519948;
constexpr double FAST_NORM_M2 = 637.333633378831;
constexpr double FAST_NORM_M3 = 296.564248779674;
constexpr double FAST_NORM_M4 = 86.7807322029461;
constexpr double FAST_NORM_M5 = 16.064177579207;
constexpr double FAST_NORM_M6 = 1.75566716318264;
constexpr double FAST_NORM_M7 = 8.83883476483184e-02;
constexpr double LOG_SQRT_2PI = 0.91893853320467274178;

inline double fast_norm_phi(double x) {
  const double z = std::fabs(x);
  double c = 0.0;

  if (z <= 37.0) {
    const double e = std::exp(-z * z / 2.0);
    if (z < FAST_NORM_SPLIT) {
      const double n = (((((FAST_NORM_N6 * z + FAST_NORM_N5) * z + FAST_NORM_N4) * z + FAST_NORM_N3) * z + FAST_NORM_N2) * z + FAST_NORM_N1) * z + FAST_NORM_N0;
      const double d = ((((((FAST_NORM_M7 * z + FAST_NORM_M6) * z + FAST_NORM_M5) * z + FAST_NORM_M4) * z + FAST_NORM_M3) * z + FAST_NORM_M2) * z + FAST_NORM_M1) * z + FAST_NORM_M0;
      c = e * n / d;
    } else {
      const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
      c = e / (FAST_NORM_RT2PI * f);
    }
  }

  return x <= 0.0 ? c : 1.0 - c;
}

// log(P(Z > z)) for z >= 0, computed directly without materialising the probability.
inline double fast_log_upper_tail(double z) {
  if (z > 37.0) return R_NegInf;
  const double z2 = z * z;
  if (z >= FAST_NORM_SPLIT) {
    const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
    return -0.5 * z2 - LOG_SQRT_2PI - std::log(f);
  }
  const double n = (((((FAST_NORM_N6 * z + FAST_NORM_N5) * z + FAST_NORM_N4) * z + FAST_NORM_N3) * z + FAST_NORM_N2) * z + FAST_NORM_N1) * z + FAST_NORM_N0;
  const double d = ((((((FAST_NORM_M7 * z + FAST_NORM_M6) * z + FAST_NORM_M5) * z + FAST_NORM_M4) * z + FAST_NORM_M3) * z + FAST_NORM_M2) * z + FAST_NORM_M1) * z + FAST_NORM_M0;
  return -0.5 * z2 + std::log(n) - std::log(d);
}

inline double pnorm_std(double x, bool lower = true, bool log_p = false) {
#ifdef USE_FAST_PNORM
  if (log_p) {
    const bool want_tail = (x >= 0.0) == (!lower);
    if (want_tail) {
      double p = fast_norm_phi(x);
      if (!lower) p = 1.0 - p;
      if (p > 1e-10) return std::log(p);
      return lower ? fast_log_upper_tail(-x) : fast_log_upper_tail(x);
    }
    // Bulk side: fast_norm_phi(-|x|) gives the small lower-tail probability (≤ 0.5),
    // so log1p is stable — no cancellation.
    return std::log1p(-fast_norm_phi(-std::fabs(x)));
  }
  double p = fast_norm_phi(x);
  if (!lower) p = 1.0 - p;
  return p;
#else
  return R::pnorm(x, 0.0, 1.0, lower, log_p);
#endif
}

// --------------------------------------------------------------------------
// Shared robust normal log tails and natural/log adapters.
//
// pnorm_std's fast log-tail approximation deliberately returns -Inf beyond
// z = 37 (mirroring natural-scale underflow).  Log-space race kernels need
// those tails to remain representable, so extend them with the same
// continued-fraction asymptotic the fast path uses below 37:
//   log Q(z) = -z^2/2 - log(sqrt(2*pi)) - log(f(z)).
// This stays in cheap arithmetic instead of falling back to R::pnorm, which
// profiling showed dominating truncated-LBA likelihoods.  All log-space model
// fallbacks (LBA/BAwL, LNR, RDMSWTN, ...) share this one tail implementation.
// --------------------------------------------------------------------------
inline double pnorm_log_direct(double x, bool lower = true) {
  double out = pnorm_std(x, lower, true);
  if (out == R_NegInf && emc2_isfinite(x)) {
    const double z = lower ? -x : x;  // tail argument: result is log Q(z)
    if (z > 0.0) {
      const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
      return -0.5 * z * z - LOG_SQRT_2PI - std::log(f);
    }
    out = R::pnorm(x, 0.0, 1.0, lower, true);
  }
  return out;
}

// log Q(z) = log Phi(-z), always evaluated on the small-tail side.
// fast_norm_phi computes the tail probability there directly with full
// relative accuracy, whereas the upper-tail call materialises 1 - Phi and
// loses ~7 digits to the subtraction from 1 before the CF takes over.
inline double log_normal_upper_tail(double z) {
  return pnorm_log_direct(-z, true);
}

// Return log(Phi(hi) - Phi(lo)) without materialising either probability.
// Using the upper tails when both arguments are positive avoids subtracting
// two values that are both numerically equal to one.
inline double log_normal_interval(double lo, double hi) {
  if (!(hi > lo)) return R_NegInf;
  if (lo >= 0.0) {
    return log_diff_exp(log_normal_upper_tail(lo),
                        log_normal_upper_tail(hi));
  }
  return log_diff_exp(pnorm_log_direct(hi, true),
                      pnorm_log_direct(lo, true));
}

inline double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

// Natural-space twin of log_normal_interval(): Phi(hi) - Phi(lo), taken in
// whichever tail avoids subtracting two values that both round to one.  The
// larger of the two probabilities is returned in `scale` so the caller can
// check how much of it the difference retained -- a natural difference is only
// as good as that retained fraction, and there is no way to recover it after
// the subtraction has happened.
inline double normal_interval_nat(double lo, double hi, double &scale) {
  scale = 0.0;
  if (!(hi > lo) || ISNAN(lo) || ISNAN(hi)) return 0.0;
  double a, b;
  if (lo >= 0.0) {
    a = pnorm_std(lo, false, false);   // Q(lo) >= Q(hi)
    b = pnorm_std(hi, false, false);
  } else {
    a = pnorm_std(hi, true, false);    // Phi(hi) >= Phi(lo)
    b = pnorm_std(lo, true, false);
  }
  scale = a;
  return a - b;
}

// Natural-space twin of log_normal_q_antiderivative_abs(): M(z) = phi(z) - z Q(z).
inline double normal_q_antiderivative_abs_nat(double z) {
  if (z == R_PosInf || ISNAN(z)) return 0.0;
  if (z == R_NegInf) return R_PosInf;
  const double phi = dnormP(z, 0.0, 1.0, false);
  if (z > 0.0) {
    const double q = pnorm_std(z, false, false);
    return phi - z * q;
  } else if (z < 0.0) {
    const double q = pnorm_std(z, false, false);
    return phi + (-z) * q;
  }
  return phi;
}

// Natural-space twin of log_normal_q_interval(): integral_lo^hi Q(z) dz = M(lo) - M(hi).
inline double normal_q_interval_nat(double lo, double hi, double &scale) {
  scale = 0.0;
  if (!(hi > lo) || ISNAN(lo) || ISNAN(hi)) return 0.0;
  const double m_lo = normal_q_antiderivative_abs_nat(lo);
  const double m_hi = normal_q_antiderivative_abs_nat(hi);
  scale = m_lo;
  return m_lo - m_hi;
}

// Natural-space twin of log_normal_phi_antiderivative(): h(z) = z Phi(z) + phi(z).
inline double normal_phi_antiderivative_nat(double z) {
  if (z == R_PosInf) return R_PosInf;
  if (z == R_NegInf || ISNAN(z)) return 0.0;
  const double phi = dnormP(z, 0.0, 1.0, false);
  const double phi_cdf = pnorm_std(z, true, false);
  return z * phi_cdf + phi;
}

// Natural-space twin of log_normal_phi_integral(): integral_lo^hi Phi(z) dz = h(hi) - h(lo).
inline double normal_phi_integral_nat(double lo, double hi, double &scale) {
  scale = 0.0;
  if (!(hi > lo) || ISNAN(lo) || ISNAN(hi)) return 0.0;

  if (hi <= 0.0) {
    const double h_hi = normal_phi_antiderivative_nat(hi);
    const double h_lo = normal_phi_antiderivative_nat(lo);
    scale = h_hi;
    return h_hi - h_lo;
  }

  if (lo >= 0.0) {
    const double width = hi - lo;
    double q_scale = 0.0;
    const double q_int = normal_q_interval_nat(lo, hi, q_scale);
    scale = width;
    return width - q_int;
  }

  const double h_0 = normal_phi_antiderivative_nat(0.0);
  const double h_lo = normal_phi_antiderivative_nat(lo);
  const double left = h_0 - h_lo;

  const double right_width = hi;
  double right_q_scale = 0.0;
  const double right_q_int = normal_q_interval_nat(0.0, hi, right_q_scale);
  const double right = right_width - right_q_int;

  scale = h_0 + right_width;
  return left + right;
}

// Natural-space twin of log_lognormal_stoploss(): C(v) = E[(V - v)_+].
inline double lognormal_stoploss_nat(double v, double mu, double sigma, double &scale) {
  scale = 0.0;
  if (!(sigma > 0.0) || ISNAN(v)) return 0.0;
  if (v == R_PosInf || !emc2_isfinite(v)) return 0.0;
  const double M = std::exp(mu + 0.5 * sigma * sigma);
  if (!(v > 0.0)) {
    scale = M;
    return M - v;
  }
  const double x = (std::log(v) - mu) / sigma;
  const double Q1 = pnorm_std(x - sigma, false, false);
  const double Q2 = pnorm_std(x, false, false);
  const double term1 = M * Q1;
  const double term2 = v * Q2;
  scale = term1;
  return term1 - term2;
}

// Natural-space twin for stoploss differences: C(lo) - C(hi).
inline double lognormal_stoploss_interval_nat(double lo, double hi, double mu, double sigma, double &scale) {
  scale = 0.0;
  if (!(hi > lo) || ISNAN(lo) || ISNAN(hi)) return 0.0;
  double scale_lo = 0.0, scale_hi = 0.0;
  const double c_lo = lognormal_stoploss_nat(lo, mu, sigma, scale_lo);
  const double c_hi = lognormal_stoploss_nat(hi, mu, sigma, scale_hi);
  scale = c_lo;
  return c_lo - c_hi;
}


// Let Q(z) = 1 - Phi(z). The positive quantity
//   M(z) = phi(z) - z Q(z)
// is decreasing and satisfies
//   integral_lo^hi Q(z) dz = M(lo) - M(hi).
// This is the positive-integrand form of the LBA CDF and avoids the severe
// cancellation in the usual "1 + correction" expression at early times.
inline double log_normal_q_antiderivative_abs(double z) {
  const double log_phi = dnormP(z, 0.0, 1.0, true);
  const double log_q = log_normal_upper_tail(z);

  if (z > 0.0) {
    // M(z) = phi(z) * [1 - z Q(z) / phi(z)].
    const double log_ratio = std::log(z) + log_q - log_phi;
    if (log_ratio >= 0.0) {
      // Asymptotic M(z) ~ phi(z) / z^2.  This branch is reached only when
      // the two log terms are indistinguishable at machine precision.
      return log_phi - 2.0 * std::log(z);
    }
    return log_phi + log1m_exp(log_ratio);
  }
  if (z < 0.0) {
    // M(z) = phi(z) + (-z) Q(z).
    return log_sum_exp(log_phi, std::log(-z) + log_q);
  }
  return log_phi;
}

inline double log_normal_q_interval(double lo, double hi) {
  if (!(hi > lo)) return R_NegInf;
  return log_diff_exp(log_normal_q_antiderivative_abs(lo),
                      log_normal_q_antiderivative_abs(hi));
}

// log h(z), where h(z) = z Phi(z) + phi(z), the antiderivative of Phi.
// The negative tail is evaluated as a log difference; this is the same
// cancellation-sensitive part as the LBA CDF and must not be formed on the
// natural scale.
inline double log_normal_phi_antiderivative(double z) {
  if (z == R_PosInf) return R_PosInf;
  if (z == R_NegInf) return R_NegInf;

  const double log_phi = dnormP(z, 0.0, 1.0, true);
  const double log_phi_cdf = pnorm_log_direct(z, true);
  if (z > 0.0) {
    return log_sum_exp(std::log(z) + log_phi_cdf, log_phi);
  }
  if (z < 0.0) {
    const double log_ratio = std::log(-z) + log_phi_cdf - log_phi;
    if (log_ratio >= 0.0) {
      // Asymptotic h(z) ~ phi(z) / z^2.  This branch is reached only when
      // the two log terms are indistinguishable at machine precision.
      return log_phi - 2.0 * std::log(-z);
    }
    return log_phi + log1m_exp(log_ratio);
  }
  return log_phi;
}

// Return log( integral_lo^hi Phi(z) dz ).  The positive-CDF representation
// uses the lower-tail antiderivative on the left and the complementary-tail
// antiderivative on the right, avoiding subtraction of nearly equal values.
inline double log_normal_phi_integral(double lo, double hi) {
  if (!(hi > lo)) return R_NegInf;

  if (hi <= 0.0) {
    return log_diff_exp(log_normal_phi_antiderivative(hi),
                        log_normal_phi_antiderivative(lo));
  }

  auto log_positive_interval = [](double xlo, double xhi) {
    if (!(xhi > xlo)) return R_NegInf;
    const double log_width = std::log(xhi - xlo);
    const double log_q = log_normal_q_interval(xlo, xhi);
    return log_diff_exp(log_width, log_q);
  };

  if (lo >= 0.0) {
    return log_positive_interval(lo, hi);
  }

  // Split at zero so that neither endpoint calculation subtracts terms from
  // opposite tails.
  const double left = log_diff_exp(log_normal_phi_antiderivative(0.0),
                                   log_normal_phi_antiderivative(lo));
  const double right = log_positive_interval(0.0, hi);
  return log_sum_exp(left, right);
}

// Output adapter shared by all natural-fast / log-fallback wrappers: the
// calculation scale is chosen by the caller, this only converts at the
// boundary according to the requested output scale.
inline double return_from_log(double log_value, bool log_out) {
  return log_out ? log_value : std::exp(log_value);
}

// A signed log-scale product.  The coefficient is kept on its natural scale
// as an input, but the product itself is never formed in natural space.
inline signed_log signed_log_product(double coefficient, double log_factor) {
  if (coefficient == 0.0 || log_factor == R_NegInf) {
    return {R_NegInf, 0};
  }
  return make_signed_log(std::log(std::fabs(coefficient)) + log_factor,
                         coefficient > 0.0 ? 1 : -1);
}

inline double log_phi_std(double z) {
  return -0.5 * z * z - LOG_SQRT_2PI;
}

// Mills ratio R(z) = Q(z) / phi(z) for z >= 0.  Central region evaluates the
// tail probability directly (full relative accuracy from fast_norm_phi); the
// far region uses the same continued fraction as the fast log tail, so the
// two agree at the split.
inline double mills_ratio_std(double z) {
  if (z >= FAST_NORM_SPLIT) {
    const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
    return 1.0 / f;
  }
  return fast_norm_phi(-z) * std::exp(0.5 * z * z + LOG_SQRT_2PI);
}

// log C(v), where C(v) = E[(V - v)_+] for log V ~ N(mu, sigma^2): the
// undiscounted Black call / stop-loss price.  With M = e^{mu + sigma^2/2} and
// x = (log v - mu)/sigma,
//   C(v) = M Q(x - sigma) - v Q(x).
// The two terms agree to a relative sigma/x in the upper tail, so the raw
// difference is unusable there.  M phi(x - sigma) = v phi(x) regroups it as
//   C(v) = v phi(x) [R(x - sigma) - R(x)],     R = Mills ratio,
// where the subtraction is between two order-1/x quantities computed to full
// relative accuracy.  Past a retained fraction of 1e-12 the Mills difference is
// itself replaced by its factored asymptotic expansion (in which the
// subtraction has been performed symbolically), so this function keeps its
// digits and stays finite arbitrarily deep into the tail.  Both the phi and Q
// terms are analytic, never routed through a natural dnorm/pnorm, so nothing
// here underflows to NaN.
//
// C is the antiderivative of the lognormal survivor: C'(v) = -P(V >= v), which
// is why BAwD's start-point and frozen-mass integrals both reduce to it.
inline double log_mills_gap(double x, double sigma) {
  const double y = x - sigma;
  const double r_y = mills_ratio_std(y);   // valid for y of either sign
  const double r_x = mills_ratio_std(x);
  const double d = r_y - r_x;              // > 0: R is strictly decreasing
  if (d > 1e-12 * r_y) {
    return std::log(d);
  }
  // Retained fraction below 1e-12 means sigma/x < 1e-12, i.e. x is enormous.
  // R(z) = 1/z - 1/z^3 + 3/z^5 - ... differenced symbolically:
  //   R(y) - R(x) = sigma/(x y) - sigma (x^2 + x y + y^2)/(x^3 y^3) + ...
  // whose relative error is O(1/x^2) and therefore negligible wherever this
  // branch can be reached.
  if (!(x > 0.0) || !(y > 0.0)) return R_NegInf;
  const double t1 = sigma / (x * y);
  const double t2 = sigma * (x * x + x * y + y * y) /
                    (x * x * x * y * y * y);
  const double series = t1 - t2;
  if (!(series > 0.0)) return R_NegInf;
  return std::log(series);
}

inline double log_lognormal_stoploss(double v, double mu, double sigma) {
  if (!(sigma > 0.0)) return R_NegInf;
  const double log_M = mu + 0.5 * sigma * sigma;
  if (!(v > 0.0)) return log_M;  // C(0) = E[V] = M
  const double x = (std::log(v) - mu) / sigma;
  if (x <= 0.0) {
    // Lower half: Q(x - sigma) and Q(x) are both order one and not close, so
    // the plain log difference retains its digits.
    return log_diff_exp(log_M + pnorm_log_direct(x - sigma, false),
                        std::log(v) + pnorm_log_direct(x, false));
  }
  return std::log(v) + log_phi_std(x) + log_mills_gap(x, sigma);
}
// log Lambda_m(v) = log integral_v^infinity w^(-(m+1)) P(V >= w) dw
// for m > 0, m Lambda_m(v) =
// v^(-m) phi(x) [R(x) - R(x + m sigma)], with x = (log v - mu)/sigma.
// The Mills-ratio difference is evaluated directly until its retained
// fraction is too small, then differenced symbolically to avoid tail loss.
inline double log_lognormal_power_stoploss(double v, double mu, double sigma,
                                           double m) {
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  const double x = (std::log(v) - mu) / sigma;
  if (std::fabs(m) <= 1e-14) {
    return std::log(sigma) + log_normal_q_antiderivative_abs(x);
  }
  const double a = m * sigma;
  if (x <= 0.0) {
    const double first = -m * std::log(v) + pnorm_log_direct(x, false);
    const double second = -m * mu + 0.5 * a * a +
      pnorm_log_direct(x + a, false);
    const double gap = (m > 0.0)
      ? log_diff_exp(first, second)
      : log_diff_exp(second, first);
    return gap - std::log(std::fabs(m));
  }
  const double y = x + a;
  const double r_x = mills_ratio_std(x);
  const double r_y = mills_ratio_std(y);
  const double d = (m > 0.0) ? r_x - r_y : r_y - r_x;
  if (d > 1e-12 * std::fmax(r_x, r_y)) {
    return -m * std::log(v) + log_phi_std(x) + std::log(d) -
           std::log(std::fabs(m));
  }
  const double aa = std::fabs(a);
  const double t1 = aa / (x * y);
  const double t2 = aa * (x * x + x * y + y * y) /
                    (x * x * x * y * y * y);
  const double series = t1 - t2;
  if (!(series > 0.0)) return R_NegInf;
  return -m * std::log(v) + log_phi_std(x) + std::log(series) -
         std::log(std::fabs(m));
}

inline double log_lognormal_logratio_stoploss(double v, double mu,
                                              double sigma, double log_ell) {
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  const double x = (std::log(v) - mu) / sigma;
  const double y = x - sigma;
  const double lambda = std::log(v) - log_ell;
  const double gap = log_mills_gap(x, sigma);
  if (gap == R_NegInf) return R_NegInf;
  const signed_log t1 = make_signed_log(
    std::log(std::fabs(lambda - 1.0)) + gap,
    lambda >= 1.0 ? 1 : -1);
  const signed_log t2 = make_signed_log(
    std::log(sigma) + log_normal_q_antiderivative_abs(y) - log_phi_std(y), 1);
  const auto br = signed_log_add(t1, t2);
  return (br.sign <= 0) ? R_NegInf : std::log(v) + log_phi_std(x) + br.log_abs;
}


// Shared acceptance constants for the guarded natural race kernels.
// A natural probability difference is rejected once it retains less than
// EMC2_NAT_REL_CANCEL of its largest term; a natural CDF consumed as
// log(1 - F) is rejected within EMC2_CDF_SAT_MARGIN of one; signed-log
// brackets are trusted while they keep log(1e-6) of their largest term.
constexpr double EMC2_NAT_REL_CANCEL = 1e-10;
constexpr double EMC2_CDF_SAT_MARGIN = 1e-8;
constexpr double EMC2_LOG_CANCEL_MIN = -13.815510557964274;  // log(1e-6)

// Fast log-normal CDF/PDF wrappers.
// For sdlog <= 0 or non-finite sdlog, fall back to R's implementation to preserve semantics.
inline double plnorm_std(double x, double meanlog, double sdlog,
                         bool lower_tail = true, bool log_p = false) {
  if (x <= 0.0) {
    if (log_p) return lower_tail ? R_NegInf : 0.0;
    return lower_tail ? 0.0 : 1.0;
  }
  if (!emc2_isfinite(sdlog) || !(sdlog > 0.0)) {
    return R::plnorm(x, meanlog, sdlog, lower_tail, log_p);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return pnorm_std(z, lower_tail, log_p);
}

inline double dlnorm_std(double x, double meanlog, double sdlog,
                         bool log_p = false) {
  if (x <= 0.0) return log_p ? R_NegInf : 0.0;
  if (!emc2_isfinite(sdlog) || !(sdlog > 0.0)) {
    return R::dlnorm(x, meanlog, sdlog, log_p);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  const double log_pdf = -std::log(x) - std::log(sdlog) - 0.5 * z * z - LOG_SQRT_2PI;
  return log_p ? log_pdf : std::exp(log_pdf);
}

inline double lnorm_log_surv_std(double x, double meanlog, double sdlog) {
  if (x <= 0.0) return 0.0;
  if (!emc2_isfinite(sdlog) || !(sdlog > 0.0)) {
    const double cdf = R::plnorm(x, meanlog, sdlog, true, false);
    if (cdf >= 1.0) return R_NegInf;
    if (cdf <= 0.0) return 0.0;
    return std::log1p(-cdf);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  // Extended tail: keep the log survivor finite past the fast path's z = 37
  // cutoff instead of collapsing far-tail races to -Inf.
  return pnorm_log_direct(z, false);
}

inline double pigt0(double t, double k, double l){
  if (l == 0.0) return 2.0 * pnorm_std(-k / std::sqrt(t), true, false);
  double mu = k / l;
  double lambda = k * k;

  double p1 = pnorm_std(std::sqrt(lambda/t) * (1. + t/mu), false, true);
  double p2 = pnorm_std(std::sqrt(lambda/t) * (1. - t/mu), false, false);

  return std::exp(2.0 * k * l + p1) + p2;
}

inline double digt0(double t, double k, double l){
  double lambda = k * k;
  double e;
  if (l == 0.) {
    e = -.5 * lambda / t;
  } else {
    double mu = k / l;
    e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
  }
  return std::exp(e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI));
}

inline double pigt_impl(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
  if (t <= 0.) return 0.;
  if (a < threshold) return pigt0(t, k, l);

  const double sqt = std::sqrt(t);

  if (std::abs(l) < threshold) {
    const double p1 = pnorm_std((k + a) / sqt, true, false);
    const double p2 = pnorm_std((k - a) / sqt, true, false);
    const double diff_sq = (k + a) * (k + a) - (k - a) * (k - a);
    const double t1 = sqt * std::exp(-0.5 * (k - a) * (k - a) / t) * (-std::expm1(-0.5 * diff_sq / t)) / FAST_NORM_RT2PI;
    return (t1 + (k + a) * (1.0 - p1) - (k - a) * (1.0 - p2)) / a;
  }

  const double t1a = std::exp(- .5 * (k - a - t * l) * (k - a - t * l) / t);
  const double t1b = std::exp(- .5 * (a + k - t * l) * (a + k - t * l) / t);
  const double t1 = sqt * (t1a - t1b) / FAST_NORM_RT2PI;

  const double t2a = std::exp(2. * l * (k - a) + pnorm_std(- (k - a + t * l) / sqt, true, true));
  const double t2b = std::exp(2. * l * (k + a) + pnorm_std(- (k + a + t * l) / sqt, true, true));
  const double t2 = a + (t2b - t2a) / (2. * l);

  const double t4a = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
  const double t4b = 2. * pnorm_std((k - a) / sqt - sqt * l, true, false) - 1.;
  const double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;

  return .5 * (t4 + t2 + t1) / a;
}

inline double digt_impl(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10){
  if (t <= 0.) return 0.;
  if (a < threshold) return digt0(t, k, l);

  if (std::abs(l) < threshold) {
    const double diff_sq = (k + a) * (k + a) - (k - a) * (k - a);
    const double log_term = -0.5 * (k - a) * (k - a) / t + std::log(-std::expm1(-0.5 * diff_sq / t));
    return std::exp(-0.5 * (LOG_2PI + std::log(t)) + log_term - std::log(2.0 * a));
  }

  const double sqt = std::sqrt(t);

  const double t1a = - (a - k + t * l) * (a - k + t * l) / (2. * t);
  const double t1b = - (a + k - t * l) * (a + k - t * l) / (2. * t);
  const double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

  const double t2a = 2. * pnorm_std((- k + a) / sqt + sqt * l, true, false) - 1.;
  const double t2b = 2. * pnorm_std((k + a) / sqt - sqt * l, true, false) - 1.;
  const double t2 = 0.5 * l * (t2a + t2b);

  return (t1 + t2) / (2.0 * a);
}

// Global symbols for exported functions (defined with [[Rcpp::export]] in wald_functions.cpp)
double pigt(double t, double k, double l, double a, double threshold);
double digt(double t, double k, double l, double a, double threshold);

// --------------------------------------------------------------------------
// Stable point-start Wald (inverse Gaussian) log primitives, sigma = 1.
// Distance d > 0, drift mu.  These are the authoritative log evaluators
// behind the guarded natural k = 0 formulas below; they are implemented from
// the Wald-form expressions (not the legacy digt/pigt formulas).
// --------------------------------------------------------------------------
inline double wald_pt_log_pdf(double t, double d, double mu) {
  if (!(t > 0.0) || t == R_PosInf || !(d > 0.0)) return R_NegInf;
  const double z = (d - mu * t) / std::sqrt(t);
  return std::log(d) - 1.5 * std::log(t) + log_phi_std(z);
}

// log F(t) = log[ Phi((mu t - d)/sqrt(t)) + e^{2 mu d} Phi((-mu t - d)/sqrt(t)) ].
// Both terms are positive, so the log CDF is a plain log_sum_exp and stays
// finite far into the early-time tail where the natural CDF underflows.
inline double wald_pt_log_cdf(double t, double d, double mu) {
  if (!(t > 0.0)) return R_NegInf;
  if (!(d > 0.0)) return 0.0;  // start at/above threshold: absorbed at t = 0
  if (t == R_PosInf) return (mu >= 0.0) ? 0.0 : 2.0 * mu * d;
  const double sqt = std::sqrt(t);
  const double log_p1 = pnorm_log_direct((mu * t - d) / sqt, true);
  const double log_p2 = 2.0 * mu * d + pnorm_log_direct((-mu * t - d) / sqt, true);
  const double out = log_sum_exp(log_p1, log_p2);
  return out > 0.0 ? 0.0 : out;
}

// log S(t) = log[ Phi(z1) - e^{2 mu d} Phi(z2) ], z1 = (d - mu t)/sqrt(t),
// z2 = -(d + mu t)/sqrt(t).  In the far upper tail (mu > 0, z1 < 0) the two
// terms share the same Gaussian exponent exactly, so the difference is formed
// through Mills ratios: S = phi(z1) * (R(-z1) - R(-z2)).
inline double wald_pt_log_surv(double t, double d, double mu) {
  if (!(d > 0.0)) return R_NegInf;
  if (!(t > 0.0)) return 0.0;
  if (t == R_PosInf) return (mu >= 0.0) ? R_NegInf : log1m_exp(2.0 * mu * d);

  const double log_cdf = wald_pt_log_cdf(t, d, mu);
  if (log_cdf < -M_LN2) return log1m_exp(log_cdf);  // F < 1/2: well conditioned

  const double sqt = std::sqrt(t);
  const double z1 = (d - mu * t) / sqt;
  const double z2 = -(d + mu * t) / sqt;

  if (z1 >= 0.0) {
    // Exact regrouping S = [Phi(z1) - Phi(z2)] - expm1(2 mu d) Phi(z2):
    // the interval term is evaluated as a stable log interval and the
    // correction is O(mu d), so the small-d corner keeps its digits.
    const signed_log interval = make_signed_log(log_normal_interval(z2, z1), 1);
    const double e = 2.0 * mu * d;
    signed_log corr{R_NegInf, 0};
    if (e != 0.0) {
      const double log_em1 = (e > 0.0) ? e + log1m_exp(-e) : log1m_exp(e);
      corr = make_signed_log(log_em1 + pnorm_log_direct(z2, true),
                             (e > 0.0) ? -1 : 1);
    }
    const signed_log total = signed_log_add(interval, corr);
    if (total.sign > 0 && !ISNAN(total.log_abs))
      return std::fmin(total.log_abs, 0.0);
    return R_NegInf;  // survivor below representable resolution here
  }

  // z1 < 0 (mu > 0 past d/mu): Mills-ratio difference, positive by
  // monotonicity of R.
  const double x1 = -z1, x2 = -z2;  // 0 < x1 < x2
  const double r1 = mills_ratio_std(x1);
  const double r2 = mills_ratio_std(x2);
  const double dr = r1 - r2;
  if (dr > 0.0 && dr > 1e-6 * r1)
    return std::fmin(log_phi_std(z1) + std::log(dr), 0.0);
  // Subtraction lost its digits: leading order R(x) ~ 1/x gives
  // R(x1) - R(x2) ~ (x2 - x1)/(x1 x2) with x2 - x1 = 2 d / sqrt(t).
  return std::fmin(log_phi_std(z1) + std::log(2.0 * d) - 0.5 * std::log(t) -
                   std::log(x1) - std::log(x2), 0.0);
}

// --------------------------------------------------------------------------
// Wald SPV CDF/PDF without kill. The 4-argument versions are the sigma = 1
// primitives; the 5-argument versions keep physical-scale parameters at the
// call site and apply Brownian scale invariance inside this helper.
// Parameterised in canonical Wald form:
//   start x ~ U[0, A], threshold b, so distance d = b - x ~ U[b-A, b].
// The CDF is the closed-form integral of the point-Wald CDF over d.
// --------------------------------------------------------------------------
inline double pwald_k0(double t, double b, double mu, double A) {
  if (t <= 0.0 || b <= 0.0) return 0.0;
  if (A <= 1e-12) return pigt0(t, b, mu);
  if (A <= 0.0) return 0.0;

  const double d_raw_lo = b - A;
  const double d_lo = std::fmax(0.0, d_raw_lo);
  const double d_hi = b;
  const double immediate = (d_raw_lo < 0.0) ? std::fmin(-d_raw_lo, A) : 0.0;
  if (d_hi <= d_lo) return std::fmax(0.0, std::fmin(1.0, immediate / A));

  const double sqt = std::sqrt(t);
  auto phi_std = [](double z) {
    return std::exp(-0.5 * z * z) / FAST_NORM_RT2PI;
  };

  if (std::abs(mu) <= 1e-12) {
    auto antideriv_mu0 = [&](double d) {
      const double z = -d / sqt;
      return 2.0 * (d * pnorm_std(z, true, false) - sqt * phi_std(z));
    };
    const double p = (immediate + antideriv_mu0(d_hi) - antideriv_mu0(d_lo)) / A;
    return std::fmax(0.0, std::fmin(1.0, p));
  }

  auto antideriv = [&](double d) {
    const double z1 = (mu * t - d) / sqt;
    const double z2 = -(d + mu * t) / sqt;
    const double term1 = (d - mu * t) * pnorm_std(z1, true, false) - sqt * phi_std(z1);
    const double term2 = (std::exp(2.0 * mu * d + pnorm_std(z2, true, true)) -
                          pnorm_std(z1, true, false)) / (2.0 * mu);
    return term1 + term2;
  };

  const double p = (immediate + antideriv(d_hi) - antideriv(d_lo)) / A;
  return std::fmax(0.0, std::fmin(1.0, p));
}

inline double pwald_k0(double t, double b, double mu, double A, double s) {
  if (!(s > 0.0)) return 0.0;
  const double inv_s = 1.0 / s;
  return pwald_k0(t, b * inv_s, mu * inv_s, A * inv_s);
}

// Guarded natural-space SPV Wald density.  Returns true when `pdf` holds a
// value the caller may use directly (including exact zeros for out-of-support
// inputs); false when the natural formula under/overflowed or cancelled, in
// which case dwald_k0_log() is authoritative.
inline bool dwald_k0_natural(double t, double b, double mu, double A, double& pdf) {
  pdf = 0.0;
  if (t <= 0.0 || b <= 0.0) return true;  // exact zero
  const double b_lo  = b - A;
  if (A <= 1e-12) {
    pdf = digt0(t, b, mu);
    if (!emc2_isfinite(pdf)) return false;
    return pdf > 0.0;
  }
  const double sqt   = std::sqrt(t);
  const double inv_A = 1.0 / A;

  if (std::abs(mu) <= 1e-12) {
    const double t1 = M_SQRT1_2 * (std::exp(-0.5 * b_lo * b_lo / t)
                                  -std::exp(-0.5 * b    * b    / t))
                    / (std::sqrt(M_PI) * sqt);
    pdf = t1 * inv_A;
    if (!emc2_isfinite(pdf)) return false;
    return pdf > 0.0 || b_lo < 0.0;  // b_lo < 0: keep legacy continuation
  }

  const double t1a = -(b_lo - t * mu) * (b_lo - t * mu) / (2.0 * t);
  const double t1b = -(b    - t * mu) * (b    - t * mu) / (2.0 * t);
  const double t1  = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

  const double t2a = 2.0 * pnorm_std(-b_lo / sqt + sqt * mu, true, false) - 1.0;
  const double t2b = 2.0 * pnorm_std( b    / sqt - sqt * mu, true, false) - 1.0;
  const double t2  = 0.5 * mu * (t2a + t2b);

  const double sum = t1 + t2;
  if (!emc2_isfinite(sum)) return false;
  pdf = sum * inv_A;
  if (b_lo < 0.0) return true;  // outside supported b >= A regime: legacy value
  if (!(sum > 0.0)) return false;                    // underflow or sign loss
  if (sum <= EMC2_NAT_REL_CANCEL * (std::fabs(t1) + std::fabs(t2)))
    return false;                                    // cancellation
  return true;
}

// Authoritative log SPV Wald density.  In z-space the start-point average is
//   f(t) = [ mu * (Phi(z_hi) - Phi(z_lo)) + (phi(z_lo) - phi(z_hi))/sqrt(t) ] / A,
// z(d) = (d - mu t)/sqrt(t) at d = b - A and d = b: the same bracket structure
// as the LBA density, handled with the same signed-log arithmetic.
inline double dwald_k0_log(double t, double b, double mu, double A) {
  if (t <= 0.0 || t == R_PosInf || b <= 0.0) return R_NegInf;
  if (A <= 1e-12) return wald_pt_log_pdf(t, b, mu);
  const double d_lo = b - A;
  if (d_lo < 0.0) return R_NegInf;  // callers keep the legacy natural value here

  const double sqt = std::sqrt(t);
  const double z_lo = (d_lo - mu * t) / sqt;
  const double z_hi = (b - mu * t) / sqt;
  const double log_dphi = log_normal_interval(z_lo, z_hi);
  const signed_log lead = signed_log_product(mu, log_dphi);
  signed_log term2 = signed_log_sub(make_signed_log(log_phi_std(z_lo), 1),
                                    make_signed_log(log_phi_std(z_hi), 1));
  if (term2.sign != 0)
    term2 = make_signed_log(term2.log_abs - 0.5 * std::log(t), term2.sign);
  const signed_log total = signed_log_add(lead, term2);
  if (total.sign > 0 && !ISNAN(total.log_abs)) {
    const double max_term = std::max(lead.log_abs, term2.log_abs);
    if (total.log_abs > max_term + EMC2_LOG_CANCEL_MIN)
      return total.log_abs - std::log(A);
  }
  // The terms cancelled past the accuracy of the tail logs (the true density
  // is positive): midpoint-in-d limit of the start-point average.
  return wald_pt_log_pdf(t, b - 0.5 * A, mu);
}

inline double dwald_k0(double t, double b, double mu, double A) {
  double pdf;
  if (dwald_k0_natural(t, b, mu, A, pdf)) return pdf;
  return std::exp(dwald_k0_log(t, b, mu, A));
}

inline double dwald_k0(double t, double b, double mu, double A, double s) {
  if (!(s > 0.0)) return 0.0;
  const double inv_s = 1.0 / s;
  return dwald_k0(t, b * inv_s, mu * inv_s, A * inv_s);
}

inline double dwald_k0_log(double t, double b, double mu, double A, double s) {
  if (!(s > 0.0)) return R_NegInf;
  const double inv_s = 1.0 / s;
  return dwald_k0_log(t, b * inv_s, mu * inv_s, A * inv_s);
}

// --------------------------------------------------------------------------
// Authoritative log CDF / log survivor for the k = 0 SPV Wald.  The A > 0
// start-point average is formed in log space over the point-start primitives
// with a fixed Gauss-Legendre rule; this path is only entered when the
// natural closed form has saturated or underflowed, so the quadrature cost is
// exceptional, not per-trial.
// --------------------------------------------------------------------------
template <typename PtLogFn>
inline double wald_k0_log_avg_over_d(double d_lo, double d_hi, PtLogFn&& pt_log_fn) {
  // The point log integrand can fall by hundreds of log units across the
  // start range in the regimes that reach this fallback, so a single global
  // rule under-resolves the dominant endpoint.  Panelled Gauss-Legendre
  // keeps the per-panel variation manageable; cost is exceptional-path only.
  constexpr int n_panels = 4;
  const GLRule& gl = gl_get_rule(20);
  const double panel_w = (d_hi - d_lo) / n_panels;
  const double half = 0.5 * panel_w;
  double log_acc = R_NegInf;
  for (int p = 0; p < n_panels; ++p) {
    const double center = d_lo + (p + 0.5) * panel_w;
    for (size_t j = 0; j < gl.x.size(); ++j) {
      if (!(gl.w[j] > 0.0)) continue;
      const double lf = pt_log_fn(center + half * gl.x[j]);
      if (lf == R_NegInf) continue;
      log_acc = log_sum_exp(log_acc, std::log(gl.w[j]) + lf);
    }
  }
  // integral over [d_lo, d_hi] = half * sum over panels/nodes of w_j f_j
  return log_acc + std::log(half);
}

// Antiderivative core for the reflected component of the point-Wald CDF:
//   T(d) = e^{2 mu d} Phi(z2(d)) - Phi(z1(d)),
//   z1 = (mu t - d)/sqrt(t), z2 = -(d + mu t)/sqrt(t),
// satisfying dT/dd = 2 mu e^{2 mu d} Phi(z2(d)) (the phi terms cancel through
// the exact identity e^{2 mu d} phi(z2) = phi(z1)).  That same identity makes
// the two terms of T share one Gaussian exponent, so when both z's are in the
// lower tail T is regrouped through Mills ratios:
//   T = phi(z1) * (R(-z2) - R(-z1)).
inline signed_log wald_k0_cdf_exp_antideriv(double t, double d, double mu) {
  const double sqt = std::sqrt(t);
  const double z1 = (mu * t - d) / sqt;
  const double z2 = -(d + mu * t) / sqt;
  if (z1 < 0.0 && z2 < 0.0) {
    const double x1 = -z1, x2 = -z2;
    const double r1 = mills_ratio_std(x1);
    const double r2 = mills_ratio_std(x2);
    const double adr = std::fabs(r2 - r1);
    if (!(adr > 1e-6 * std::fmax(r1, r2))) {
      // Subtraction lost its digits: leading order R(x) ~ 1/x gives
      // |R(x2) - R(x1)| ~ |x1 - x2|/(x1 x2) with |x1 - x2| = 2 |mu| sqrt(t).
      if (mu == 0.0) return make_signed_log(R_NegInf, 0);
      return make_signed_log(log_phi_std(z1) + std::log(2.0 * std::fabs(mu)) +
                                 0.5 * std::log(t) - std::log(x1) - std::log(x2),
                             (mu > 0.0) ? -1 : 1);
    }
    return make_signed_log(log_phi_std(z1) + std::log(adr), (r2 > r1) ? 1 : -1);
  }
  return signed_log_sub(
      make_signed_log(2.0 * mu * d + pnorm_log_direct(z2, true), 1),
      make_signed_log(pnorm_log_direct(z1, true), 1));
}

// Closed-form log CDF for the A > 0 start-point average, built from the two
// positive components of the point CDF:
//   A * F_A(t) = int Phi((mu t - d)/sqrt(t)) dd
//              + int e^{2 mu d} Phi(z2(d)) dd    over d in [d_lo, d_hi].
// The first integral maps onto the stable log-int-Phi machinery shared with
// the LBA; the second integrates to (T(d_hi) - T(d_lo)) / (2 mu).  Returns
// NaN when an outer signed-log difference has cancelled past its accuracy;
// the caller then uses the panelled quadrature instead.
inline double wald_k0_log_cdf_closed(double t, double mu, double A,
                                     double d_lo, double d_hi) {
  const double sqt = std::sqrt(t);
  const double log_I1 = 0.5 * std::log(t) +
      log_normal_phi_integral((mu * t - d_hi) / sqt, (mu * t - d_lo) / sqt);
  if (ISNAN(log_I1)) return NA_REAL;

  double log_I2;
  if (std::fabs(mu) <= 1e-12) {
    // mu = 0: the reflected component equals the direct one exactly.
    log_I2 = log_I1;
  } else {
    const signed_log T_hi = wald_k0_cdf_exp_antideriv(t, d_hi, mu);
    const signed_log T_lo = wald_k0_cdf_exp_antideriv(t, d_lo, mu);
    const signed_log diff = signed_log_sub(T_hi, T_lo);
    // I2 = (T(d_hi) - T(d_lo)) / (2 mu) is positive: T is monotone in d with
    // slope of mu's sign.  Reject sign failures and differences that kept
    // less than log(1e-6) of their largest term.
    const int want_sign = (mu > 0.0) ? 1 : -1;
    if (diff.sign != want_sign || ISNAN(diff.log_abs)) return NA_REAL;
    const double max_term = std::fmax(T_hi.log_abs, T_lo.log_abs);
    if (emc2_isfinite(max_term) &&
        diff.log_abs <= max_term + EMC2_LOG_CANCEL_MIN)
      return NA_REAL;
    log_I2 = diff.log_abs - std::log(2.0 * std::fabs(mu));
  }
  const double out = log_sum_exp(log_I1, log_I2) - std::log(A);
  return ISNAN(out) ? NA_REAL : out;
}

inline double wald_k0_log_cdf(double t, double b, double mu, double A) {
  if (!(t > 0.0) || !(b > 0.0)) return R_NegInf;
  if (A <= 1e-12) return wald_pt_log_cdf(t, b, mu);
  const double d_raw_lo = b - A;
  const double d_lo = std::fmax(0.0, d_raw_lo);
  const double d_hi = b;
  const double immediate = (d_raw_lo < 0.0) ? std::fmin(-d_raw_lo, A) : 0.0;
  if (d_hi <= d_lo) {
    const double p = std::fmax(0.0, std::fmin(1.0, immediate / A));
    return p > 0.0 ? std::log(p) : R_NegInf;
  }
  double out = NA_REAL;
  if (emc2_isfinite(t)) {
    out = wald_k0_log_cdf_closed(t, mu, A, d_lo, d_hi);
  }
  if (ISNAN(out)) {
    // t = Inf or an outer difference cancelled: panelled quadrature over the
    // stable point log CDF.
    out = wald_k0_log_avg_over_d(d_lo, d_hi, [&](double d) {
      return wald_pt_log_cdf(t, d, mu);
    }) - std::log(A);
  }
  if (immediate > 0.0)
    out = log_sum_exp(out, std::log(immediate) - std::log(A));
  return std::fmin(out, 0.0);
}

inline double wald_k0_log_surv(double t, double b, double mu, double A) {
  if (!(b > 0.0)) return R_NegInf;
  if (!(t > 0.0)) return 0.0;
  if (A <= 1e-12) return wald_pt_log_surv(t, b, mu);
  const double d_raw_lo = b - A;
  const double d_lo = std::fmax(0.0, d_raw_lo);
  const double d_hi = b;
  if (d_hi <= d_lo) return R_NegInf;  // all starts at/above threshold
  // Starts with d <= 0 are absorbed immediately and contribute no survivor
  // mass, so the average runs over the positive-d part only.
  const double out = wald_k0_log_avg_over_d(d_lo, d_hi, [&](double d) {
    return wald_pt_log_surv(t, d, mu);
  }) - std::log(A);
  return std::fmin(out, 0.0);
}

inline double wald_k0_log_cdf(double t, double b, double mu, double A, double s) {
  if (!(s > 0.0)) return R_NegInf;
  const double inv_s = 1.0 / s;
  return wald_k0_log_cdf(t, b * inv_s, mu * inv_s, A * inv_s);
}

inline double wald_k0_log_surv(double t, double b, double mu, double A, double s) {
  if (!(s > 0.0)) return R_NegInf;
  const double inv_s = 1.0 / s;
  return wald_k0_log_surv(t, b * inv_s, mu * inv_s, A * inv_s);
}

#endif

#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <cstring>
#include "utility_functions.h"
#include "composite_functions.h"
#include "gl_quad.h"

using namespace Rcpp;


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

// Standard-normal fast path.  R::dnorm is an out-of-line call into Rmath that
// re-validates mean/sd on every invocation (~10.5 ns against ~4.9 ns for the
// expression below); all but a handful of the call sites in this package pass
// (0, 1).  Rmath's own |x| < 5 branch IS this expression and its log branch is
// -(M_LN_SQRT_2PI + x^2/2), so the values are bit-identical where this fires;
// |x| >= 5, non-standard mean/sd and NaN all defer to Rmath, which is where its
// split-argument high-precision tail actually matters.
inline double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  if (mean == 0.0 && sd == 1.0) {
    const double z2 = x * x;
    if (log) return -0.5 * z2 - LOG_SQRT_2PI;
    if (z2 < 25.0) return 0.398942280401432677939946059934 * std::exp(-0.5 * z2);
  }
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

// Lognormal launch parameterization.  The SAMPLED pair for the plain
// lognormal launch is the natural-scale arithmetic mean m = E[V] and the
// coefficient of variation cv = SD(V)/E[V]; every downstream primitive works
// in (mu, sigma) on log V.  The map is a bijection,
//   sigma^2 = log(1 + cv^2),   mu = log m - sigma^2/2,
// with inverse m = exp(mu + sigma^2/2), cv = sqrt(expm1(sigma^2)).  log1p
// keeps sigma accurate as cv -> 0, where sigma -> cv.
//
// The SPLIT-lognormal launch is NOT reparameterized: it keeps (mu, sigma,
// delta) with mu the exact median, because inverting (m, cv, delta) has no
// closed form.  Callers therefore convert only when the launch code is the
// plain lognormal one -- that is what `meancv` selects.
struct LaunchLogNormal { double mu; double sigma; };

inline LaunchLogNormal launch_lognormal_pair(double p1, double p2,
                                             bool meancv) {
  if (!meancv) return {p1, p2};
  if (!(p1 > 0.0) || !(p2 > 0.0) || !R_finite(p1) || !R_finite(p2))
    return {R_NaN, R_NaN};
  const double s2 = std::log1p(p2 * p2);
  return {std::log(p1) - 0.5 * s2, std::sqrt(s2)};
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

// Lognormal put primitive: int_0^v Phi((log w - mu)/sigma) dw, i.e. the
// integral of the launch CDF rather than its survivor.  Keeping this primitive
// in log space avoids the cancellation in 1 - stoploss/mean that would defeat
// the survivor-tail path.
inline double log_lognormal_put(double v, double mu, double sigma) {
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  const double x = (std::log(v) - mu) / sigma;
  const double a = std::log(v) + pnorm_log_direct(x, true);
  const double b = mu + 0.5 * sigma * sigma +
    pnorm_log_direct(x - sigma, true);
  if (!(a > b)) return R_NegInf;
  const double d = b - a;
  return a + log1m_exp(d);
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

// --------------------------------------------------------------------------
// Cached log-gamma and fast regularized incomplete gamma.
//
// Both are Weibull-launch hot-path primitives.  Every Weibull form needs
// log Gamma at an order fixed by the launch shape (1 + 1/shape for the scale
// and the first moment, 1/shape for the stop-loss), so the call is a per-row
// constant that R::lgammafn otherwise recomputes at ~45 ns a row; memoise it
// on that key, as model_BAwD.cpp does for the saturation solve.  Constant
// initialisation keeps the thread_local free of a guard variable, and key 0
// is never queried because the orders are positive.
// --------------------------------------------------------------------------
inline double lgammafn_cached(double a) {
  if (!(a > 0.0) || !emc2_isfinite(a)) return R::lgammafn(a);
  static constexpr int SLOTS = 128;  // power of two: the hash masks
  struct Slot { double key; double val; };
  static thread_local Slot tab[SLOTS] = {};
  std::uint64_t h;
  std::memcpy(&h, &a, sizeof(h));
  // murmur3 finalizer: orders separated only in low mantissa bits must not
  // land in one slot, which the raw bit pattern does not guarantee.
  h ^= h >> 33; h *= 0xff51afd7ed558ccdull;
  h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ull;
  h ^= h >> 33;
  Slot& s = tab[h & (SLOTS - 1)];
  if (s.key == a) return s.val;
  s.key = a;
  s.val = R::lgammafn(a);
  return s.val;
}

// P(a,z) and Q(a,z), the regularized lower and upper incomplete gamma.
// R::pgamma's Didonato-Morris path costs ~190 ns and dominated the Weibull
// launch the way R::pnorm once dominated the truncated LBA; the classical
// series/continued-fraction pair below holds ~3e-11 relative accuracy over
// a in [A_MIN, A_MAX] and the whole representable z range at ~50 ns.  The
// band is where that pair is both accurate and cheap: below it 1 - P eats
// the upper tail, above it the expansions need hundreds of terms near
// z = a.  Outside, R::pgamma is called, carrying the uniform asymptotics
// that extreme orders need.
//
// log z is an argument because every call site already holds it: z is always
// (v / scale)^shape and is formed in log space, so passing it in removes the
// only std::log an evaluation would otherwise need.
namespace inc_gamma {

constexpr int MAX_ITER = 400;
constexpr double REL_EPS = 1e-16;
// The continued fraction converges on |del - 1|, so its tolerance must stay
// clear of one ulp of 1.0 or the loop runs to MAX_ITER instead of stopping.
constexpr double CF_EPS = 1e-15;
constexpr double TINY = 1e-300;
constexpr double A_MIN = 1e-4;
constexpr double A_MAX = 1e3;

inline bool in_band(double a, double z) {
  return a >= A_MIN && a <= A_MAX && z > 0.0 && emc2_isfinite(z);
}

// log of the exp(-z) z^a / Gamma(a) factor shared by both expansions.  It is
// maximised at z = a, where Stirling makes it negative, so exp() cannot
// overflow.
inline double log_prefactor(double a, double z, double log_z) {
  return a * log_z - z - lgammafn_cached(a);
}

// P(a,z) = exp(log_prefactor) * lower_series(a,z), for z < a + 1.
inline double lower_series(double a, double z) {
  double ap = a, del = 1.0 / a, sum = del;
  for (int n = 0; n < MAX_ITER; ++n) {
    ap += 1.0;
    del *= z / ap;
    sum += del;
    if (std::fabs(del) <= std::fabs(sum) * REL_EPS) break;
  }
  return sum;
}

// Q(a,z) = exp(log_prefactor) * upper_cf(a,z), for z >= a + 1.  Modified
// Lentz evaluation of the Legendre continued fraction.
inline double upper_cf(double a, double z) {
  double b = z + 1.0 - a, c = 1.0 / TINY, d = 1.0 / b, h = d;
  for (int i = 1; i < MAX_ITER; ++i) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;  if (std::fabs(d) < TINY) d = TINY;
    c = b + an / c;  if (std::fabs(c) < TINY) c = TINY;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) <= CF_EPS) break;
  }
  return h;
}

}  // namespace inc_gamma

// Regularized lower incomplete gamma on the natural scale.
inline double gamma_p_fast(double a, double z, double log_z) {
  if (!inc_gamma::in_band(a, z)) return R::pgamma(z, a, 1.0, true, false);
  const double pref = std::exp(inc_gamma::log_prefactor(a, z, log_z));
  if (z < a + 1.0)
    return std::fmin(pref * inc_gamma::lower_series(a, z), 1.0);
  return std::fmax(1.0 - pref * inc_gamma::upper_cf(a, z), 0.0);
}

// log P(a,z) and log Q(a,z); either pointer may be null.  Each tail is taken
// from the expansion that computes it directly, so the complement is only
// ever formed on the side that keeps its digits.
inline void gamma_log_pq_fast(double a, double z, double log_z,
                              double* log_p, double* log_q) {
  if (!inc_gamma::in_band(a, z)) {
    if (log_p) *log_p = R::pgamma(z, a, 1.0, true, true);
    if (log_q) *log_q = R::pgamma(z, a, 1.0, false, true);
    return;
  }
  const double lpref = inc_gamma::log_prefactor(a, z, log_z);
  if (z < a + 1.0) {
    const double lp =
      std::fmin(lpref + std::log(inc_gamma::lower_series(a, z)), 0.0);
    if (log_p) *log_p = lp;
    if (log_q) *log_q = log1m_exp(lp);
    return;
  }
  const double lq =
    std::fmin(lpref + std::log(inc_gamma::upper_cf(a, z)), 0.0);
  if (log_q) *log_q = lq;
  if (log_p) *log_p = log1m_exp(lq);
}

// --------------------------------------------------------------------------
// Weibull launch primitives.  V ~ Weibull(shape, mean), V >= 0.  The
// corresponding conventional scale is mean / Gamma(1 + 1 / shape).
// The incomplete-gamma forms are kept in log space because all ballistic
// kernels consume stop-loss, put, and partial-moment differences.
// --------------------------------------------------------------------------
inline bool weibull_valid(double shape, double mean) {
  return shape > 0.0 && mean > 0.0 && emc2_isfinite(shape) &&
    emc2_isfinite(mean);
}

inline double weibull_log_scale(double shape, double mean) {
  if (!weibull_valid(shape, mean)) return R_NaN;
  return std::log(mean) - lgammafn_cached(1.0 + 1.0 / shape);
}

inline double weibull_scale_from_mean(double shape, double mean) {
  const double ls = weibull_log_scale(shape, mean);
  if (ISNAN(ls)) return R_NaN;
  return ls == R_NegInf ? 0.0 : (ls == R_PosInf ? R_PosInf : std::exp(ls));
}

inline double rweibull_mean(double shape, double mean) {
  if (!weibull_valid(shape, mean)) return R_NaN;
  const double u = R::unif_rand();
  const double log_x = weibull_log_scale(shape, mean) +
    std::log(-std::log1p(-u)) / shape;
  return std::exp(log_x);
}

inline double weibull_log_z(double v, double shape, double mean) {
  if (!(v > 0.0)) return (v == 0.0) ? R_NegInf : R_NaN;
  const double log_scale = weibull_log_scale(shape, mean);
  if (ISNAN(log_scale)) return R_NaN;
  const double lz = shape * (std::log(v) - log_scale);
  return lz;
}

inline double log_weibull_survivor(double v, double shape, double mean) {
  if (!weibull_valid(shape, mean)) return R_NegInf;
  if (!(v > 0.0)) return 0.0;
  const double lz = weibull_log_z(v, shape, mean);
  return (lz == R_PosInf) ? R_NegInf : -std::exp(lz);
}

inline double log_weibull_cdf(double v, double shape, double mean) {
  if (!weibull_valid(shape, mean) || !(v > 0.0)) return R_NegInf;
  const double lz = weibull_log_z(v, shape, mean);
  if (lz == R_PosInf) return 0.0;
  // In the extreme lower tail z may underflow even though log F is still a
  // perfectly representable number.  log(1 - exp(-z)) ~ log z there.
  if (lz < -36.0) return lz;
  const double z = std::exp(lz);
  if (!(z > 0.0)) return R_NegInf;
  return std::log(-std::expm1(-z));
}

inline double log_weibull_density(double v, double shape, double mean) {
  if (!weibull_valid(shape, mean) || !(v > 0.0)) return R_NegInf;
  const double log_scale = weibull_log_scale(shape, mean);
  const double lz = weibull_log_z(v, shape, mean);
  if (lz == R_PosInf) return R_NegInf;
  return std::log(shape) - log_scale +
    (shape - 1.0) * (std::log(v) - log_scale) - std::exp(lz);
}

// log Gamma(s,z), upper incomplete gamma, for arbitrary real s and z > 0.
// The regularized tail covers s > 0; the recurrence extends it to the
// negative orders needed by the power-stoploss terms.
inline double log_weibull_upper_gamma(double s, double z, double log_z) {
  if (!(z > 0.0) || !emc2_isfinite(s)) return R_NegInf;
  if (s > 0.0) {
    double lq;
    gamma_log_pq_fast(s, z, log_z, nullptr, &lq);
    return (lq > R_NegInf) ? lgammafn_cached(s) + lq : R_NegInf;
  }
  if (std::fabs(s) < 1e-12) {
    // E1(z) = Gamma(0,z).  Power series for z < 1 and a continued fraction
    // for the upper tail avoid both cancellation and underflow.
    if (z < 1.0) {
      double term = 1.0, sum = 0.0;
      for (int n = 1; n < 200; ++n) {
        term *= -z / static_cast<double>(n);
        const double add = term / static_cast<double>(n);
        sum += add;
        if (std::fabs(add) <= 2e-16 * std::fmax(1.0, std::fabs(sum))) break;
      }
      const double e1 = -0.5772156649015328606 - log_z - sum;
      return (e1 > 0.0) ? std::log(e1) : R_NegInf;
    }
    double b = z + 1.0, c = 1e300, d = 1.0 / b, h = d;
    for (int i = 1; i < 200; ++i) {
      const double a = -static_cast<double>(i * i);
      b += 2.0;
      d = a * d + b; if (std::fabs(d) < 1e-300) d = 1e-300;
      c = b + a / c; if (std::fabs(c) < 1e-300) c = 1e-300;
      d = 1.0 / d;
      const double del = d * c;
      h *= del;
      if (std::fabs(del - 1.0) < 2e-15) break;
    }
    const double e1 = std::exp(-z) * h;
    return (e1 > 0.0) ? std::log(e1) : R_NegInf;
  }
  const int nshift = static_cast<int>(std::floor(-s)) + 1;
  const double q = s + nshift;
  double log_g = log_weibull_upper_gamma(q, z, log_z);
  for (int j = nshift - 1; j >= 0; --j) {
    const double r = s + j;
    if (std::fabs(r) < 1e-12) {
      log_g = log_weibull_upper_gamma(0.0, z, log_z);
      continue;
    }
    const double log_term = r * log_z - z;
    const signed_log num = signed_log_sub(make_signed_log(log_g, 1),
                                          make_signed_log(log_term, 1));
    if (num.sign == 0) return R_NegInf;
    log_g = num.log_abs - std::log(std::fabs(r));
  }
  return log_g;
}

inline double log_weibull_stoploss(double v, double shape, double mean) {
  if (!weibull_valid(shape, mean)) return R_NegInf;
  const double log_mean = std::log(mean);
  const double log_scale = weibull_log_scale(shape, mean);
  if (!(v > 0.0)) return log_mean;
  const double lz = weibull_log_z(v, shape, mean);
  if (lz < -36.0) {
    // C(v) = E[V] - v + E[(v - V)_+].  The put term is negligible in this
    // tail, but retaining -v avoids a visible bias for large shape values.
    const double log_v = std::log(v);
    if (log_mean > log_v)
      return log_mean + log1m_exp(log_v - log_mean);
    return log_mean;
  }
  const double z = std::exp(lz);
  if (!emc2_isfinite(z)) return R_NegInf;
  return log_scale - std::log(shape) +
    log_weibull_upper_gamma(1.0 / shape, z, lz);
}

inline double log_weibull_put(double v, double shape, double mean) {
  if (!weibull_valid(shape, mean) || !(v > 0.0)) return R_NegInf;
  const double log_mean = std::log(mean);
  const double log_scale = weibull_log_scale(shape, mean);
  // The put primitive is E[(v - V)_+], so it grows as v - E[V] in the
  // upper tail.  Handle the endpoint before forming log(v) - log_mean.
  if (v == R_PosInf) return R_PosInf;
  const double lz = weibull_log_z(v, shape, mean);
  if (lz == R_PosInf || lz >= 700.0) {
    const double log_v = std::log(v);
    return (log_v > log_mean)
      ? log_v + log1m_exp(log_mean - log_v)
      : log_v;
  }
  if (lz < -36.0)
    return log_scale + (1.0 + 1.0 / shape) * lz - std::log(shape + 1.0);
  const double z = std::exp(lz);
  // For small z, integrate the CDF series directly; x F(x) minus a lower
  // incomplete gamma loses all digits when shape is large.
  if (z < 0.5) {
    const double r = 1.0 / shape;
    double term = std::exp((r + 1.0) * std::log(z)) / (shape + 1.0);
    double sum = term;
    for (int n = 2; n < 200; ++n) {
      term *= -z * (shape * (n - 1.0) + 1.0) /
        (static_cast<double>(n) * (shape * n + 1.0));
      sum += term;
      if (std::fabs(term) <= 2e-16 * std::fmax(1.0, std::fabs(sum))) break;
    }
    if (sum > 0.0 && emc2_isfinite(sum))
      return log_scale + std::log(sum);
    double direct = 0.0, powz = std::exp((r + 1.0) * std::log(z));
    for (int n = 1; n < 200; ++n) {
      const double add = ((n & 1) ? 1.0 : -1.0) * powz /
        (std::tgamma(n + 1.0) * (shape * n + 1.0));
      direct += add;
      powz *= z;
      if (std::fabs(add) <= 2e-16 * std::fmax(1.0, std::fabs(direct))) break;
    }
    return direct > 0.0 ? log_scale + std::log(direct) : R_NegInf;
  }
  const double x = v;
  const double log_f = log_weibull_cdf(v, shape, mean);
  const double a_shape = 1.0 + 1.0 / shape;
  double log_pl;
  gamma_log_pq_fast(a_shape, z, lz, &log_pl, nullptr);
  const double log_m = log_scale + (lgammafn_cached(a_shape) + log_pl);
  const signed_log out = signed_log_sub(
    make_signed_log(std::log(x) + log_f, 1), make_signed_log(log_m, 1));
  return out.sign > 0 ? out.log_abs : R_NegInf;
}

inline double log_weibull_mass_interval(double lo, double hi,
                                        double shape, double mean) {
  if (!(hi > lo) || !(lo >= 0.0) || !weibull_valid(shape, mean)) return R_NegInf;
  const double lzlo = (lo > 0.0) ? weibull_log_z(lo, shape, mean) : R_NegInf;
  const double lzhi = weibull_log_z(hi, shape, mean);
  if (!(lzhi > lzlo)) return R_NegInf;
  // If both z endpoints are tiny, differencing exp(-z) would erase the mass;
  // the leading term is z_hi - z_lo and is evaluated directly in log space.
  if (lzhi < -20.0)
    return log_diff_exp(lzhi, lzlo);
  const double zlo = (lzlo > R_NegInf) ? std::exp(lzlo) : 0.0;
  const double zhi = std::exp(lzhi);
  if (!(zhi > zlo)) return R_NegInf;
  return -zlo + log1m_exp(-(zhi - zlo));
}

inline double log_weibull_gamma_interval(double s, double zlo, double zhi,
                                         double log_zlo, double log_zhi) {
  if (!(s > 0.0) || !(zhi > zlo)) return R_NegInf;
  if (zhi < 1e-8) {
    const double a = s * log_zhi;
    const double b = (zlo > 0.0) ? s * log_zlo : R_NegInf;
    if (a > b) return log_diff_exp(a, b) - std::log(s);
  }
  double lp_lo, lq_lo, lp_hi, lq_hi;
  gamma_log_pq_fast(s, zlo, log_zlo, &lp_lo, &lq_lo);
  gamma_log_pq_fast(s, zhi, log_zhi, &lp_hi, &lq_hi);
  const double lower = (lp_hi > lp_lo) ? log_diff_exp(lp_hi, lp_lo) : R_NegInf;
  const double upper = (lq_lo > lq_hi) ? log_diff_exp(lq_lo, lq_hi) : R_NegInf;
  return lgammafn_cached(s) + std::fmax(lower, upper);
}

inline double log_weibull_first_interval(double lo, double hi,
                                         double shape, double mean) {
  if (!(hi > lo) || !(lo >= 0.0) || !weibull_valid(shape, mean)) return R_NegInf;
  const double lzlo = (lo > 0.0) ? weibull_log_z(lo, shape, mean) : R_NegInf;
  const double lzhi = weibull_log_z(hi, shape, mean);
  const double zlo = (lzlo > R_NegInf) ? std::exp(lzlo) : 0.0;
  const double zhi = std::exp(lzhi);
  return weibull_log_scale(shape, mean) +
    log_weibull_gamma_interval(1.0 + 1.0 / shape, zlo, zhi, lzlo, lzhi);
}

inline double log_weibull_power_stoploss(double v, double shape, double mean,
                                         double m) {
  if (!weibull_valid(shape, mean) || !(v > 0.0)) return R_NegInf;
  const double log_scale = weibull_log_scale(shape, mean);
  const double lz = weibull_log_z(v, shape, mean);
  if (lz < -36.0) {
    const double s = -m / shape;
    if (s > 1e-12)
      return -m * log_scale - std::log(shape) + lgammafn_cached(s);
    if (std::fabs(s) <= 1e-12)
      return -m * log_scale - std::log(shape) + std::log(-lz);
    // Gamma(s,z) ~ -z^s / s for s < 0 as z -> 0.
    return -m * log_scale - std::log(shape) + s * lz - std::log(-s);
  }
  const double z = std::exp(lz);
  if (!(z > 0.0) || !emc2_isfinite(z)) return R_NegInf;
  return -m * log_scale - std::log(shape) +
    log_weibull_upper_gamma(-m / shape, z, lz);
}

// Integral of log(w / ell) times the Weibull survivor.  A short numerical
// derivative of the same incomplete-gamma primitive is stable over the range
// used by the rho = 1 BAwD member and avoids a second special-function stack.
inline double log_weibull_logratio_stoploss(double v, double shape,
                                            double mean, double log_ell) {
  if (!weibull_valid(shape, mean) || !(v > 0.0)) return R_NegInf;
  const double h = 1e-5;
  const double lm = log_weibull_power_stoploss(v, shape, mean, -1.0 - h);
  const double lp = log_weibull_power_stoploss(v, shape, mean, -1.0 + h);
  const double l0 = log_weibull_power_stoploss(v, shape, mean, -1.0);
  if (!(lm > R_NegInf) || !(lp > R_NegInf) || !(l0 > R_NegInf)) return R_NegInf;
  const double deriv = (lp - lm) / (2.0 * h);
  const double coeff = -deriv - log_ell;
  return (coeff > 0.0 && emc2_isfinite(coeff)) ? l0 + std::log(coeff) : R_NegInf;
}


inline double dlnorm_std(double x, double meanlog, double sdlog, bool log_p = false);
inline double lnorm_log_surv_std(double x, double meanlog, double sdlog);
// --------------------------------------------------------------------------
// Continuous split-lognormal launch primitives.  Let Y = log(V) have split
// point c and widths sL,sR, with a=sL/(sL+sR).  Each side is a half-normal
// piece with total mass 2a and 2(1-a).  The location parameter mu is the
// median; c is solved from the closed-form median equations below.
// --------------------------------------------------------------------------
struct split_lognormal_shape {
  double c = 0.0;
  double sL = 0.0;
  double sR = 0.0;
  double a = 0.5;
};

inline bool split_lognormal_shape_params(double mu, double sigma, double delta,
                                         split_lognormal_shape& h) {
  if (!(sigma > 0.0) || !emc2_isfinite(mu) || !emc2_isfinite(delta))
    return false;
  if (delta == 0.0) {
    h.c = mu; h.sL = sigma; h.sR = sigma; h.a = 0.5;
    return true;
  }
  h.sL = sigma * std::exp(0.5 * delta);
  h.sR = sigma * std::exp(-0.5 * delta);
  if (!(h.sL > 0.0) || !(h.sR > 0.0) ||
      !emc2_isfinite(h.sL) || !emc2_isfinite(h.sR)) return false;
  h.a = 1.0 / (1.0 + std::exp(-delta));
  if (delta > 0.0)
    h.c = mu - h.sL * R::qnorm(1.0 / (4.0 * h.a), 0.0, 1.0, true, false);
  else
    h.c = mu - h.sR * R::qnorm((3.0 * h.sR - h.sL) /
                               (4.0 * h.sR), 0.0, 1.0, true, false);
  return emc2_isfinite(h.c);
}

inline double split_lognormal_log_side_tail(double x, double r,
                                            const split_lognormal_shape& h,
                                            bool left) {
  const double s = left ? h.sL : h.sR;
  const double mass = left ? h.a : (1.0 - h.a);
  const double base = std::log(2.0 * mass) + r * h.c + 0.5 * r * r * s * s;
  if (left) {
    if (!(x < h.c)) return R_NegInf;
    const double lo = (x - h.c - r * s * s) / s;
    const double hi = -r * s;
    return base + log_normal_interval(lo, hi);
  }
  const double lo = (x > h.c) ? (x - h.c - r * s * s) / s : -r * s;
  return base + log_normal_upper_tail(lo);
}

inline double split_lognormal_log_side_lower(double x, double r,
                                             const split_lognormal_shape& h,
                                             bool left) {
  const double s = left ? h.sL : h.sR;
  const double mass = left ? h.a : (1.0 - h.a);
  const double base = std::log(2.0 * mass) + r * h.c + 0.5 * r * r * s * s;
  if (left) {
    if (!(x < h.c))
      return base + pnorm_log_direct(-r * s, true);
    return base + pnorm_log_direct((x - h.c - r * s * s) / s, true);
  }
  if (!(x > h.c)) return R_NegInf;
  return base + log_normal_interval(
    -r * s, (x - h.c - r * s * s) / s);
}

inline double split_lognormal_log_moment_tail(double v, const split_lognormal_shape& h, double r) {
  const double lpL = split_lognormal_log_side_tail(std::log(v), r, h, true);
  const double lpR = split_lognormal_log_side_tail(std::log(v), r, h, false);
  return log_sum_exp(lpL, lpR);
}

inline double split_lognormal_log_moment_tail(double v, double mu,
                                              double sigma, double delta,
                                              double r) {
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return split_lognormal_log_moment_tail(v, h, r);
}

inline double log_split_lognormal_survivor(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, h, 0.0);
}

inline double log_split_lognormal_survivor(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return lnorm_log_surv_std(v, mu, sigma);
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 0.0);
}

inline double log_split_lognormal_first_partial_moment(double v, const split_lognormal_shape& h) {
  return split_lognormal_log_moment_tail(v, h, 1.0);
}

inline double log_split_lognormal_first_partial_moment(double v, double mu,
                                                       double sigma,
                                                       double delta) {
  if (delta == 0.0) {
    if (!(v > 0.0)) return mu + 0.5 * sigma * sigma;
    const double x = (std::log(v) - mu) / sigma;
    return mu + 0.5 * sigma * sigma +
      pnorm_log_direct(x - sigma, false);
  }
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 1.0);
}

inline double split_lognormal_log_total_moment(const split_lognormal_shape& h, double r) {
  return split_lognormal_log_moment_tail(0.0, h, r);
}

inline double split_lognormal_log_total_moment(double mu, double sigma,
                                               double delta, double r) {
  return split_lognormal_log_moment_tail(0.0, mu, sigma, delta, r);
}

inline double split_lognormal_log_stoploss_side(double x, double v,
                                                const split_lognormal_shape& h,
                                                bool left) {
  const double lm = split_lognormal_log_side_tail(x, 1.0, h, left);
  const double lp = split_lognormal_log_side_tail(x, 0.0, h, left);
  if (!(lm > R_NegInf) || !(lp > R_NegInf)) return R_NegInf;
  const signed_log out = signed_log_sub(
    make_signed_log(lm, 1), make_signed_log(std::log(v) + lp, 1));
  return out.sign > 0 ? out.log_abs : R_NegInf;
}
inline double log_split_lognormal_stoploss(double v, const split_lognormal_shape& h) {
  const double lm = split_lognormal_log_total_moment(h, 1.0);
  if (!(v > 0.0)) {
    if (v == 0.0) return lm;
    return log_sum_exp(lm, std::log(-v));
  }
  const double x = std::log(v);
  double ll = R_NegInf;
  if (x > h.c) {
    const double z = (x - h.c) / h.sR;
    const double r1 = mills_ratio_std(z - h.sR);
    const double r0 = mills_ratio_std(z);
    const double gap = r1 - r0;
    if (gap > 1e-12 * r1)
      ll = std::log(2.0 * (1.0 - h.a)) + std::log(v) +
        log_phi_std(z) + std::log(gap);
  }
  if (!(ll > R_NegInf))
    ll = log_sum_exp(split_lognormal_log_stoploss_side(x, v, h, true),
                     split_lognormal_log_stoploss_side(x, v, h, false));
  return ll;
}

inline double log_split_lognormal_stoploss(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return log_lognormal_stoploss(v, mu, sigma);
  if (!(sigma > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_stoploss(v, h);
}

inline double log_split_lognormal_put(double v, const split_lognormal_shape& h) {
  const double x = std::log(v);
  auto side_put = [&](bool left) {
    const double lp = split_lognormal_log_side_lower(x, 0.0, h, left);
    const double lm = split_lognormal_log_side_lower(x, 1.0, h, left);
    if (!(lp > R_NegInf) || !(lm > R_NegInf)) return R_NegInf;
    const signed_log out = signed_log_sub(
      make_signed_log(std::log(v) + lp, 1), make_signed_log(lm, 1));
    return out.sign > 0 ? out.log_abs : R_NegInf;
  };
  return log_sum_exp(side_put(true), side_put(false));
}

inline double log_split_lognormal_put(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) return log_lognormal_put(v, mu, sigma);
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_put(v, h);
}

inline double log_split_lognormal_power_stoploss(double v, const split_lognormal_shape& h, double m) {
  if (std::fabs(m) <= 1e-14) {
    const double x = std::log(v);
    auto yside = [&](bool left) {
      const double s = left ? h.sL : h.sR;
      const double mass = left ? h.a : (1.0 - h.a);
      const double z = (x - h.c) / s;
      if (left) {
        const double p = pnorm_std(0.0, false, false) -
          pnorm_std(z, false, false);
        const double q = dnormP(z) - dnormP(0.0);
        const double val = (h.c - x) * p + s * q;
        return val > 0.0 ? std::log(2.0 * mass * val) : R_NegInf;
      }
      const double zz = (x > h.c) ? z : 0.0;
      const double p = pnorm_std(zz, false, false);
      const double q = dnormP(zz);
      const double val = (h.c - x) * p + s * q;
      return val > 0.0 ? std::log(2.0 * mass * val) : R_NegInf;
    };
    return log_sum_exp(yside(true), yside(false));
  }
  const double lm = split_lognormal_log_moment_tail(v, h, -m);
  const double ls = log_split_lognormal_survivor(v, h);
  const signed_log out = (m > 0.0)
    ? signed_log_sub(make_signed_log(std::log(v) * (-m) + ls, 1),
                     make_signed_log(lm, 1))
    : signed_log_sub(make_signed_log(lm, 1),
                     make_signed_log(std::log(v) * (-m) + ls, 1));
  return out.sign > 0 ? out.log_abs - std::log(std::fabs(m)) : R_NegInf;
}

inline double log_split_lognormal_power_stoploss(double v, double mu,
                                                 double sigma, double delta,
                                                 double m) {
  if (delta == 0.0) return log_lognormal_power_stoploss(v, mu, sigma, m);
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_power_stoploss(v, h, m);
}

inline double log_split_lognormal_density(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return R_NegInf;
  const double x = std::log(v);
  const double s = (x <= h.c) ? h.sL : h.sR;
  const double mass = (x <= h.c) ? h.a : (1.0 - h.a);
  return std::log(2.0 * mass) - std::log(s) - std::log(v) -
    0.5 * (x - h.c) * (x - h.c) / (s * s) - LOG_SQRT_2PI;
}

inline double log_split_lognormal_density(double v, double mu, double sigma,
                                          double delta) {
  if (delta == 0.0) return dlnorm_std(v, mu, sigma, true);
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_density(v, h);
}
inline double log_split_lognormal_cdf(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return R_NegInf;
  const double x = std::log(v);
  if (x <= h.c) {
    const double z = (x - h.c) / h.sL;
    return std::fmin(std::log(2.0 * h.a) + pnorm_log_direct(z, true), 0.0);
  }
  const double z = (x - h.c) / h.sR;
  // F(x) = 2(1-a) Phi(z) + (2a - 1) on the right side.  The constant
  // term is negative when a < 1/2, so it must be subtracted rather than
  // passed to log(); using Q(z) here also makes the result decrease with x.
  const double lp = std::log(2.0 * (1.0 - h.a)) + pnorm_log_direct(z, true);
  if (!(lp > R_NegInf)) return R_NegInf;
  if (h.a >= 0.5) {
    const double lconst = std::log(2.0 * h.a - 1.0);
    return std::fmin((lconst > R_NegInf) ? log_sum_exp(lp, lconst) : lp, 0.0);
  }
  const double lconst = std::log1p(-2.0 * h.a);
  const double out = log_diff_exp(lp, lconst);
  return (out > R_NegInf) ? std::fmin(out, 0.0) : R_NegInf;
}

inline double log_split_lognormal_cdf(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) {
    if (!(v > 0.0) || !(sigma > 0.0)) return R_NegInf;
    return pnorm_log_direct((std::log(v) - mu) / sigma, true);
  }
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_cdf(v, h);
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
                         bool log_p) {
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

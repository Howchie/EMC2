#ifndef exgaussian_functions_h
#define exgaussian_functions_h

#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <cstring>
#include "composite_functions.h"
#include "wald_functions.h"
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

const double SIG_TAU_EPS = 1e-12;

// H(x) = exp(x^2 / 2) Phi(-x), evaluated without forming either a tiny
// Gaussian tail or the compensating exponential. This is the scaled normal
// tail used by the fast ex-Gaussian density path. The rational coefficients
// are Cody's normal-tail approximation (the same coefficients used in the
// R normal CDF implementation).
inline double exg_normal_tail_factor(double x, double* derivative = nullptr) {
  if (emc2_isinf(x)) {
    if (derivative) *derivative = -0.0;
    return 0.0;
  }
  if (x <= 0.67448975) {
    const double u = x * x;
    const double numerator = (((0.065682337918207449113 * u +
      2.2352520354606839287) * u + 161.02823106855587881) * u +
      1067.6894854603709582) * u + 18154.981253343561249;
    const double denominator = (((u + 47.20258190468824187) * u +
      976.09855173777669322) * u + 10260.932208618978205) * u +
      45507.789335026729956;
    const double central = x * numerator / denominator;
    const double h = std::exp(0.5 * u) * (0.5 - central);
    if (derivative) *derivative = x * h - 0.39894228040143267794;
    return h;
  }
  if (x <= 5.6568542494923801952) {
    constexpr double c[] = {1.0765576773720192317e-8, 0.39894151208813466764,
      8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
      2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124,
      9842.7148383839780218};
    constexpr double d[] = {1.0, 22.266688044328115691, 235.38790178262499861,
      1519.377599407554805, 6485.558298266760755, 18615.571640885098091,
      34900.952721145977266, 38912.003286093271411, 19685.429676859990727};
    double numerator = c[0], denominator = d[0], dn = 0.0, dd = 0.0;
    for (int i = 1; i < 9; ++i) {
      if (derivative) { dn = dn * x + numerator; dd = dd * x + denominator; }
      numerator = numerator * x + c[i];
      denominator = denominator * x + d[i];
    }
    const double h = numerator / denominator;
    if (derivative) *derivative = (dn - h * dd) / denominator;
    return h;
  }
  const double inv = 1.0 / x;
  const double u = inv * inv;
  constexpr double p[] = {0.02307344176494017303, 0.21589853405795699,
    0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
    2.9112874951168792e-5};
  constexpr double q[] = {1.0, 1.28426009614491121, 0.468238212480865118,
    0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5};
  double numerator = p[0], denominator = q[0], dn = 0.0, dd = 0.0;
  for (int i = 1; i < 6; ++i) {
    if (derivative) { dn = dn * u + numerator; dd = dd * u + denominator; }
    numerator = numerator * u + p[i];
    denominator = denominator * u + q[i];
  }
  const double r = numerator / denominator;
  if (derivative) {
    *derivative = -u * (0.39894228040143267794 - 3.0 * u * r -
      2.0 * u * u * (dn - r * dd) / denominator);
  }
  return (0.39894228040143267794 - u * r) * inv;
}

// log H(x).  In the central branch H is exp(x^2/2) times a rational factor, so
// the log is taken of the factor alone and the exponential is never formed.
inline double exg_log_normal_tail_factor(double x) {
  if (!emc2_isinf(x) && x <= 0.67448975) {
    const double u = x * x;
    const double numerator = (((0.065682337918207449113 * u +
      2.2352520354606839287) * u + 161.02823106855587881) * u +
      1067.6894854603709582) * u + 18154.981253343561249;
    const double denominator = (((u + 47.20258190468824187) * u +
      976.09855173777669322) * u + 10260.932208618978205) * u +
      45507.789335026729956;
    return 0.5 * u + std::log(0.5 - x * numerator / denominator);
  }
  return std::log(exg_normal_tail_factor(x));
}

inline double exg_normal_tail_difference(double x, double width) {
  if (!(width > 0.0)) return 0.0;
  if (width < 0.1 * (1.0 + std::max(0.0, x))) {
    constexpr double nodes[] = {0.095012509837637440185, 0.28160355077925891323,
      0.45801677765722738634, 0.61787624440264374845, 0.75540440835500303390,
      0.86563120238783174388, 0.94457502307323257608, 0.98940093499164993260};
    constexpr double weights[] = {0.18945061045506849629, 0.18260341504492358887,
      0.16915651939500238634, 0.14959598881657673208, 0.12462897125553387205,
      0.09515851168249278481, 0.06225352393864789286, 0.02715245941175409485};
    double sum = 0.0;
    for (int i = 0; i < 8; ++i) {
      const double left = x + width * (0.5 - 0.5 * nodes[i]);
      const double right = x + width * (0.5 + 0.5 * nodes[i]);
      double derivative_left, derivative_right;
      exg_normal_tail_factor(left, &derivative_left);
      exg_normal_tail_factor(right, &derivative_right);
      sum += weights[i] * (-derivative_left - derivative_right);
    }
    return 0.5 * width * sum;
  }
  return exg_normal_tail_factor(x) - exg_normal_tail_factor(x + width);
}

// log[exp(a^2/2 - a*y) Phi(y-a)]. The scaled-tail branch avoids the
// cancellation between a large exponential and a tiny Phi in the far-left
// density tail. Keeping the opposite branch in log space also handles very
// small scales without multiplying an overflow by an underflow.
inline double exg_log_weighted_tail(double y, double a, double log_a) {
  const double w = a - y;
  if (w >= 0.0) {
    if (emc2_isinf(w)) {
      // a can overflow even though log(a) and the final density are finite.
      // At this scale H(a-y) is 1/(sqrt(2*pi)*(a-y)) to machine precision.
      const double log_w = log_a + std::log1p(-y * std::exp(-log_a));
      return -0.5 * y * y - LOG_SQRT_2PI - log_w;
    }
    return -0.5 * y * y + exg_log_normal_tail_factor(w);
  }
  return a * (0.5 * a - y) + pnorm_log_direct(-w, true);
}

inline double dexg_prepared(const double x, const double mu,
                            const double inverse_sigma, const double tau,
                            const double ratio, const double log_ratio,
                            const double log_tau,
                            const bool log_d) {
  if (emc2_isinf(x)) return log_d ? R_NegInf : 0.0;
  const double y = (x - mu) * inverse_sigma;
  if (emc2_isinf(y)) return log_d ? R_NegInf : 0.0;
  const double log_out = exg_log_weighted_tail(y, ratio, log_ratio) - log_tau;
  if (emc2_isnan(log_out)) return NA_REAL;
  if (log_d) return log_out;

  const double w = ratio - y;
  if (w >= 0.0 && !emc2_isinf(ratio) && std::fabs(y) <= 37.0) {
    return std::exp(-0.5 * y * y) * exg_normal_tail_factor(w) / tau;
  }
  return std::exp(log_out);
}

// probability density function of ex-Gaussian distribution
inline double dexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {

  // minimal parameter validation - assuming all inputs are finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // protect against numerical issues due to extremely small sigma or tau values
  const double tau_p = std::max(tau, SIG_TAU_EPS);
  const double sig_p = std::max(sigma, SIG_TAU_EPS);
  return dexg_prepared(x, mu, 1.0 / sig_p, tau_p, sig_p / tau_p,
                       std::log(sig_p) - std::log(tau_p),
                       std::log(tau_p), log_d);
}

// Direct upper-tail log probability of the ex-Gaussian:
//   S(q) = Q(u) + exp(k) * Phi(u - sigma/tau),
//   u = (q - mu)/sigma, k = (mu - q)/tau + sigma^2/(2 tau^2).
// Both terms are non-negative, so log S is a plain log_sum_exp of two log
// terms and never round-trips through the lower-tail CDF (no exp, no
// rounding to one, no log1m).
inline double pexg_log_upper(double q, double mu, double sigma, double tau) {
  const double tau_p = std::max(tau, SIG_TAU_EPS);
  const double sig_p = std::max(sigma, SIG_TAU_EPS);
  const double u = (q - mu) / sig_p;
  const double log_q_term = pnorm_log_direct(u, false);
  const double log_weighted_tail = exg_log_weighted_tail(
    u, sig_p / tau_p, std::log(sig_p) - std::log(tau_p));
  return std::fmin(0.0, log_sum_exp(log_q_term, log_weighted_tail));
}

inline double pexg(double q, double mu, double sigma, double tau,
                   bool lower_tail, bool log_p);

// Parameter-only values for a truncated ex-Gaussian. Raw likelihood batches
// can reuse these when every active row in a particle shares the same tuple.
// Keep both natural and log quantities: ptexg's ordinary path uses natural
// subtraction where it is well conditioned and its tail fallback uses logs.
struct ExgTruncationCache {
  double mu = NA_REAL;
  double sigma = NA_REAL;
  double tau = NA_REAL;
  double inverse_sigma = NA_REAL;
  double tau_prepared = NA_REAL;
  double sigma_tau_ratio = NA_REAL;
  double log_sigma_tau_ratio = NA_REAL;
  double log_tau = NA_REAL;
  double lower = NA_REAL;
  double upper = NA_REAL;
  double lower_cdf = NA_REAL;
  double upper_cdf = NA_REAL;
  double lower_lcdf = NA_REAL;
  double upper_lcdf = NA_REAL;
  double lower_lsurv = NA_REAL;
  double upper_lsurv = NA_REAL;
  double natural_normalizer = NA_REAL;
  double log_normalizer = NA_REAL;
  bool valid = false;

  bool matches(double mu_, double sigma_, double tau_, double lower_,
               double upper_) const {
    return valid && std::memcmp(&mu, &mu_, sizeof(double)) == 0 &&
      std::memcmp(&sigma, &sigma_, sizeof(double)) == 0 &&
      std::memcmp(&tau, &tau_, sizeof(double)) == 0 &&
      std::memcmp(&lower, &lower_, sizeof(double)) == 0 &&
      std::memcmp(&upper, &upper_, sizeof(double)) == 0;
  }
};

inline ExgTruncationCache make_exg_truncation_cache(
    double mu, double sigma, double tau, double lower, double upper) {
  ExgTruncationCache cache;
  cache.mu = mu;
  cache.sigma = sigma;
  cache.tau = tau;
  const double sigma_p = std::max(sigma, SIG_TAU_EPS);
  cache.inverse_sigma = 1.0 / sigma_p;
  cache.tau_prepared = std::max(tau, SIG_TAU_EPS);
  cache.sigma_tau_ratio = sigma_p / cache.tau_prepared;
  cache.log_sigma_tau_ratio = std::log(sigma_p) - std::log(cache.tau_prepared);
  cache.log_tau = std::log(cache.tau_prepared);
  cache.lower = lower;
  cache.upper = upper;
  cache.lower_cdf = lower == R_NegInf ? 0.0
      : pexg(lower, mu, sigma, tau, true, false);
  cache.upper_cdf = upper == R_PosInf ? 1.0
      : pexg(upper, mu, sigma, tau, true, false);
  cache.lower_lcdf = lower == R_NegInf ? R_NegInf
                                        : pexg(lower, mu, sigma, tau, true, true);
  cache.upper_lcdf = upper == R_PosInf ? 0.0
                                        : pexg(upper, mu, sigma, tau, true, true);
  cache.lower_lsurv = lower == R_NegInf ? 0.0
      : (lower == R_PosInf ? R_NegInf : pexg_log_upper(lower, mu, sigma, tau));
  cache.upper_lsurv = upper == R_PosInf ? R_NegInf
      : (upper == R_NegInf ? 0.0 : pexg_log_upper(upper, mu, sigma, tau));
  cache.natural_normalizer = cache.upper_cdf - cache.lower_cdf;
  if (cache.lower_lcdf == R_NegInf) {
    cache.log_normalizer = cache.upper_lcdf;
  } else if (cache.lower_lcdf > -M_LN2) {
    cache.log_normalizer = log_diff_exp(cache.lower_lsurv, cache.upper_lsurv);
  } else {
    cache.log_normalizer = log_diff_exp(cache.upper_lcdf, cache.lower_lcdf);
  }
  cache.valid = true;
  return cache;
}

// GSL scalar callbacks revisit the same distribution parameters at many
// integration points. Keep one exact, thread-local truncation context for
// these routes; the full parameter and bound key makes reuse safe across
// particles, models and forked workers without a separate invalidation hook.
inline const ExgTruncationCache* exg_truncation_cache_get(
    ExgTruncationCache& cache, double mu, double sigma, double tau,
    double lower, double upper) {
  if (!R_FINITE(mu) || !R_FINITE(sigma) || !R_FINITE(tau) ||
      sigma <= 0.0 || tau <= 0.0 || ISNAN(lower) || ISNAN(upper) ||
      !(lower < upper) ||
      (lower == R_NegInf && upper == R_PosInf)) return nullptr;
  if (!cache.matches(mu, sigma, tau, lower, upper))
    cache = make_exg_truncation_cache(mu, sigma, tau, lower, upper);
  return &cache;
}

struct ExgTruncationMemo {
  ExgTruncationCache cache;

  const ExgTruncationCache* get(double mu, double sigma, double tau,
                                double lower, double upper) {
    return exg_truncation_cache_get(cache, mu, sigma, tau, lower, upper);
  }
};

inline ExgTruncationMemo& exg_scalar_truncation_memo() {
  static thread_local ExgTruncationMemo memo;
  return memo;
}

// cumulative distribution function of ex-Gaussian distribution
inline double pexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {

  // minimal parameter validation - assuming mu, sigma, and tau are all finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // handle infinite q
  if (emc2_isinf(q)) {
    double out = (q < 0.) ? 0. : 1.;
    if (!lower_tail) out = 1. - out;
    return log_p ? std::log(out) : out;
  }

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  double log_cdf_lower;
  const double y = (q - mu) / sig_p;
  if (y <= 0.0) {
    // F(q) = exp(-y^2/2) [H(-y) - H(sigma/tau-y)]. The scaled-tail
    // difference avoids subtracting two tiny CDF terms after separately
    // evaluating the cancelling exponential and normal tail.
    const double difference = exg_normal_tail_difference(-y, sig_p / tau_p);
    log_cdf_lower = difference > 0.0
      ? -0.5 * y * y + std::log(difference) : R_NegInf;
  } else {
    // Above the Gaussian mean, form F as 1-S from the direct positive
    // survivor terms. This avoids cancellation between Phi(y) and the
    // exponentially weighted Phi(y-sigma/tau) term.
    log_cdf_lower = log1m_exp(pexg_log_upper(q, mu, sig_p, tau_p));
  }

  double out;
  if (lower_tail) {
    out = log_p ? log_cdf_lower : std::exp(log_cdf_lower);
  } else {
    if (log_cdf_lower == R_NegInf) {
      out = log_p ? 0. : 1.;
    } else {
      double cdf_lower = std::exp(log_cdf_lower);
      if (cdf_lower >= 1. - 1e-8) {
        // Near saturation 1 - cdf has lost its digits; switch to the direct
        // upper-tail expression instead of rounding the survivor to zero.
        const double log_surv = pexg_log_upper(q, mu, sig_p, tau_p);
        out = log_p ? log_surv : std::exp(log_surv);
      } else {
        out = log_p ? log1m(cdf_lower) : -std::expm1(log_cdf_lower);
      }
    }
  }

  return(out);
}

// probability density function of truncated ex-Gaussian distribution
inline double dtexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = R_NegInf,
    const double upper = R_PosInf,
    const bool log_d = false,
    const ExgTruncationCache* cache = nullptr
) {

  if (lower == R_NegInf && upper == R_PosInf) {
    return dexg(x, mu, sigma, tau, log_d);
  }
  if (sigma <= 0. || tau <= 0.) return NA_REAL;
  if (lower >= upper) return NA_REAL;
  if (x <= lower || x >= upper) return log_d ? R_NegInf : 0.;

  const bool use_cache = cache != nullptr && cache->matches(mu, sigma, tau, lower, upper);
  const double x_ld = use_cache
    ? dexg_prepared(x, mu, cache->inverse_sigma, cache->tau_prepared,
                    cache->sigma_tau_ratio, cache->log_sigma_tau_ratio,
                    cache->log_tau, true)
    : dexg(x, mu, sigma, tau, true);
  if (x_ld == R_NegInf) return log_d ? R_NegInf : 0.;

  double log_normaliser;
  if (use_cache) {
    log_normaliser = cache->log_normalizer;
  } else {
    double lower_lcdf, upper_lcdf;
    if (lower == R_NegInf) {
      lower_lcdf = R_NegInf;
    } else {
      lower_lcdf = pexg(lower, mu, sigma, tau, true, true);
    }
    if (upper == R_PosInf) {
      upper_lcdf = 0.;
    } else {
      upper_lcdf = pexg(upper, mu, sigma, tau, true, true);
    }
    if (lower_lcdf == R_NegInf) {
      log_normaliser = upper_lcdf;
    } else if (lower_lcdf > -M_LN2) {
      // Both bounds sit in the upper half: log-CDF subtraction loses digits.
      const double ls_lo = pexg_log_upper(lower, mu, sigma, tau);
      const double ls_up = (upper == R_PosInf) ? R_NegInf
                                               : pexg_log_upper(upper, mu, sigma, tau);
      log_normaliser = log_diff_exp(ls_lo, ls_up);
    } else {
      log_normaliser = log_diff_exp(upper_lcdf, lower_lcdf);
    }
  }
  if (!(log_normaliser > R_NegInf) || ISNAN(log_normaliser)) {
    return log_d ? R_NegInf : 0.;
  }

  double log_out = x_ld - log_normaliser;
  return log_d ? log_out : std::exp(log_out);
}

// cumulative distribution function of truncated ex-Gaussian distribution
inline double ptexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = R_NegInf,
    const double upper = R_PosInf,
    const bool lower_tail = true,
    const bool log_p = false,
    const ExgTruncationCache* cache = nullptr
) {

  if (lower == R_NegInf && upper == R_PosInf) {
    return pexg(q, mu, sigma, tau, lower_tail, log_p);
  }
  if (sigma <= 0. || tau <= 0.) return NA_REAL;
  if (lower >= upper) return NA_REAL;
  if (q <= lower) {
    double out = lower_tail ? 0. : 1.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }
  if (q >= upper) {
    double out = lower_tail ? 1. : 0.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }

  const bool use_cache = cache != nullptr && cache->matches(mu, sigma, tau, lower, upper);
  double q_cdf = pexg(q, mu, sigma, tau);

  double lower_cdf, upper_cdf;
  if (use_cache) {
    lower_cdf = cache->lower_cdf;
    upper_cdf = cache->upper_cdf;
  } else {
    if (lower == R_NegInf) {
      lower_cdf = 0.;
    } else {
      lower_cdf = pexg(lower, mu, sigma, tau);
    }
    if (upper == R_PosInf) {
      upper_cdf = 1.;
    } else {
      upper_cdf = pexg(upper, mu, sigma, tau);
    }
  }

  // Natural fast path: ordinary bounds where the subtractions keep their
  // digits.  Otherwise fall back to matched log differences below.
  const double normaliser = use_cache ? cache->natural_normalizer : upper_cdf - lower_cdf;
  if (normaliser > 1e-10 * std::max(upper_cdf, 1e-300)) {
    double out;
    if (lower_tail) {
      out = (q_cdf - lower_cdf) / normaliser;
    } else {
      out = (upper_cdf - q_cdf) / normaliser;
    }
    out = std::max(0., std::min(1., out));
    // A zero here can only be numerator cancellation (q is strictly inside
    // the bounds), which matters when the caller wants a log tail.
    if (out > 0. || !log_p) {
      return log_p ? std::log(out) : out;
    }
  }

  // Log-space fallback: form numerator and normaliser as log differences on
  // the better-conditioned side (lower-tail logs on the left, survivor logs
  // on the right of the median).
  auto lcdf = [&](double x) {
    if (use_cache && x == lower) return cache->lower_lcdf;
    if (use_cache && x == upper) return cache->upper_lcdf;
    return pexg(x, mu, sigma, tau, true, true);
  };
  auto lsurv = [&](double x) {
    if (use_cache && x == lower) return cache->lower_lsurv;
    if (use_cache && x == upper) return cache->upper_lsurv;
    return (x == R_PosInf) ? R_NegInf : pexg_log_upper(x, mu, sigma, tau);
  };
  const double lower_lcdf = (lower == R_NegInf) ? R_NegInf : lcdf(lower);

  double log_norm;
  if (use_cache) {
    log_norm = cache->log_normalizer;
  } else if (lower_lcdf == R_NegInf) {
    log_norm = (upper == R_PosInf) ? 0.0 : lcdf(upper);
  } else if (lower_lcdf > -M_LN2) {
    log_norm = log_diff_exp(lsurv(lower), lsurv(upper));
  } else {
    log_norm = log_diff_exp((upper == R_PosInf) ? 0.0 : lcdf(upper), lower_lcdf);
  }
  if (!(log_norm > R_NegInf) || ISNAN(log_norm)) return NA_REAL;

  double log_num;
  if (lower_tail) {
    // F(q) - F(lower) = S(lower) - S(q)
    if (lower_lcdf == R_NegInf) {
      log_num = lcdf(q);
    } else if (lower_lcdf > -M_LN2) {
      log_num = log_diff_exp(lsurv(lower), lsurv(q));
    } else {
      log_num = log_diff_exp(lcdf(q), lower_lcdf);
    }
  } else {
    // F(upper) - F(q) = S(q) - S(upper)
    const double q_lcdf = lcdf(q);
    if (q_lcdf > -M_LN2) {
      log_num = log_diff_exp(lsurv(q), lsurv(upper));
    } else {
      log_num = log_diff_exp((upper == R_PosInf) ? 0.0 : lcdf(upper), q_lcdf);
    }
  }
  if (ISNAN(log_num)) log_num = R_NegInf;
  double log_out = log_num - log_norm;
  if (ISNAN(log_out) || log_out == R_NegInf) return log_p ? R_NegInf : 0.;
  if (log_out > 0.) log_out = 0.;
  return log_p ? log_out : std::exp(log_out);
}

// Scalar versions for GSL integration
inline double dtexg_scalar(double x, const double* par, void* /*ctx*/) {
  // par: mu, sigma, tau, ..., lb=8
  const double lower = par[8];
  const ExgTruncationCache* cache =
    exg_scalar_truncation_memo().get(par[0], par[1], par[2], lower, R_PosInf);
  return dtexg(x, par[0], par[1], par[2], lower, R_PosInf, false, cache);
}

inline double ptexg_scalar(double x, const double* par, void* /*ctx*/) {
  // par: mu, sigma, tau, ..., lb=8
  const double lower = par[8];
  const ExgTruncationCache* cache =
    exg_scalar_truncation_memo().get(par[0], par[1], par[2], lower, R_PosInf);
  return ptexg(x, par[0], par[1], par[2], lower, R_PosInf, true, false, cache);
}

double dexg_scalar(double x, const double* par, void* /*ctx*/);

double pexg_scalar(double x, const double* par, void* /*ctx*/);

#endif

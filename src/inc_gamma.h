#ifndef EMC2_INC_GAMMA_H
#define EMC2_INC_GAMMA_H

#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include "utility_functions.h"
#include "composite_functions.h"

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

// Stirling series of log Gamma(y) - [(y - 1/2) log y - y + log(2 pi)/2]
// through 1/y^13; below rounding for y >= 12.
inline double lgamma_stirling_series(double y) {
  const double r = 1.0 / y, r2 = r * r;
  return r * (1.0 / 12 + r2 * (-1.0 / 360 + r2 * (1.0 / 1260
    + r2 * (-1.0 / 1680 + r2 * (1.0 / 1188 + r2 * (-691.0 / 360360
    + r2 * (1.0 / 156)))))));
}

// log Gamma(a) for a >= 0.5, for orders that change per row and so miss any
// cache: the recurrence shifts to y >= 12, then the Stirling series.  Within
// ~1e-14 absolute of R::lgammafn (the cancellation against log(prod) near
// a = 1), at a third of its cost.
inline double lgamma_pos(double a) {
  double y = a, prod = 1.0;
  while (y < 12.0) { prod *= y; y += 1.0; }
  return (y - 0.5) * std::log(y) - y + 0.91893853320467274178 +
    lgamma_stirling_series(y) - std::log(prod);
}

inline double lgamma_fast(double a) {
  return (a >= 0.5 && emc2_isfinite(a)) ? lgamma_pos(a) : R::lgammafn(a);
}

// log B(a, b) for a, b >= 0.5.  With the larger argument l >= 12,
// lgamma(l) - lgamma(s + l) is the Stirling difference with
// log(s + l) = log l + log1p(s / l), so nothing cancels however large l is.
// Within ~2e-14 relative of R::lbeta at under half its cost.
inline double lbeta_fast(double a, double b) {
  const double sm = std::fmin(a, b), lg = std::fmax(a, b);
  if (!(sm >= 0.5) || !emc2_isfinite(lg)) return R::lbeta(a, b);
  if (lg < 12.0) return lgamma_pos(sm) + lgamma_pos(lg) - lgamma_pos(sm + lg);
  const double sl = sm + lg;
  return lgamma_pos(sm) - sm * std::log(lg) - (sl - 0.5) * std::log1p(sm / lg) +
    sm + lgamma_stirling_series(lg) - lgamma_stirling_series(sl);
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
inline double log_prefactor(double a, double z, double log_z,
                            double lgamma_a) {
  return a * log_z - z - lgamma_a;
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
  const double pref =
    std::exp(inc_gamma::log_prefactor(a, z, log_z, lgammafn_cached(a)));
  if (z < a + 1.0)
    return std::fmin(pref * inc_gamma::lower_series(a, z), 1.0);
  return std::fmax(1.0 - pref * inc_gamma::upper_cf(a, z), 0.0);
}

// log P(a,z) and log Q(a,z); either pointer may be null.  Each tail is taken
// from the expansion that computes it directly, so the complement is only
// ever formed on the side that keeps its digits.  lgamma_a is log Gamma(a),
// passed by callers that already hold it.
inline void gamma_log_pq_fast(double a, double z, double log_z,
                              double lgamma_a, double* log_p, double* log_q) {
  if (!inc_gamma::in_band(a, z)) {
    if (log_p) *log_p = R::pgamma(z, a, 1.0, true, true);
    if (log_q) *log_q = R::pgamma(z, a, 1.0, false, true);
    return;
  }
  const double lpref = inc_gamma::log_prefactor(a, z, log_z, lgamma_a);
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

inline void gamma_log_pq_fast(double a, double z, double log_z,
                              double* log_p, double* log_q) {
  gamma_log_pq_fast(a, z, log_z,
                    inc_gamma::in_band(a, z) ? lgammafn_cached(a) : 0.0,
                    log_p, log_q);
}

#endif  // EMC2_INC_GAMMA_H

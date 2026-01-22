#ifndef PROB_UTILS_H
#define PROB_UTILS_H

#include <cmath>
#include <limits>

// Generic probability/CDF clamping helpers shared across models.
// These are intentionally lightweight and header-only.

// Clamp finite values to [0, 1]. Non-finite values are passed through unchanged.
double clamp_prob01(double p) {
  if (!std::isfinite(p)) return p;
  if (p < 0.0) return 0.0;
  if (p > 1.0) return 1.0;
  return p;
}

// Alias for code that thinks in terms of CDFs.
double clamp_cdf01(double cdf) {
  return clamp_prob01(cdf);
}

// Variant used when downstream logic prefers treating non-finite CDFs as "finished".
double clamp_cdf01_or_one(double cdf) {
  if (!std::isfinite(cdf)) return 1.0;
  if (cdf <= 0.0) return 0.0;
  if (cdf >= 1.0) return 1.0;
  return cdf;
}

// CDF clamping variant used in race likelihood code:
// - non-finite => NaN (so downstream checks/logs treat as invalid)
// - values extremely close to 1 => clamp to 1 - eps
double clamp_cdf01_race(double cdf, double one_minus_eps = 1e-15) {
  if (!std::isfinite(cdf)) return std::numeric_limits<double>::quiet_NaN();
  if (cdf <= 0.0) return 0.0;
  if (cdf >= 1.0) return 1.0;
  if (cdf > 1.0 - one_minus_eps) return 1.0 - one_minus_eps;
  return cdf;
}

// Stable log(1 - p) for p in [0,1], with common edge-case handling:
// - p <= 0  => log(1) = 0
// - p >= 1  => log(0) = -Inf
// - p ~ 1   => clamp away from 1 to avoid log(0) in finite arithmetic
// - non-finite => -Inf (treat as impossible)
double safe_log1m_prob(double p, double one_minus_eps = 1e-15) {
  if (!std::isfinite(p)) return -std::numeric_limits<double>::infinity();
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return -std::numeric_limits<double>::infinity();
  if (p > 1.0 - one_minus_eps) p = 1.0 - one_minus_eps;
  return std::log1p(-p);
}

// Backward-compatible name used in race likelihood code.
double safe_log1m_race(double p) { return safe_log1m_prob(p); }

constexpr double kMinSurv = 1e-12;
double safe_surv_from_cdf(double cdf) {
    const double F = clamp_prob01(cdf);
    if (!R_FINITE(F)) return NA_REAL;
    const double s = std::fma(-F, 1.0, 1.0); // 1 - F
    return (s >= kMinSurv) ? s : kMinSurv;
}
double safe_surv_from_prod_cdf(double cdf1, double cdf2) {
    const double F1 = clamp_prob01(cdf1);
    const double F2 = clamp_prob01(cdf2);
    if (!R_FINITE(F1) || !R_FINITE(F2)) return NA_REAL;
    const double s = std::fma(-F1, F2, 1.0); // 1 - F1*F2 (single-rounding)
    return (s >= kMinSurv) ? s : kMinSurv;
}

#endif

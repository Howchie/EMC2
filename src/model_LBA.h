#ifndef lba_h
#define lba_h

// Shared LBA/BAwL model kernels and inline helpers.  Rcpp-exported wrappers
// are defined in model_LBA.cpp.
#include <Rcpp.h>
#include "utility_functions.h"
#include "wald_functions.h"  // pnorm_std() — fast normal CDF under USE_FAST_PNORM
#include "composite_functions.h"  // clamp_pos, safe_log
#include "gaussian.h"
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;

// dnormP, pnorm_log_direct, log_normal_upper_tail, log_normal_interval, the
// log ∫Phi machinery (log_normal_phi_integral and its antiderivatives), and
// the signed-log product helper are shared with the other race models and
// live in wald_functions.h next to pnorm_std().

constexpr double BAWL_K_EPS = 1e-10;
constexpr double BAWL_A_EPS = 1e-10;
constexpr double BAWL_NATURAL_Z_MAX = 7.5;
constexpr double BAWL_NATURAL_MIN_SPAN = 1e-6;
constexpr double BAWL_LOG_MIN_SPAN = 1e-8;
constexpr double BAWL_NATURAL_REL_TOL = 1e-10;
// Trust the signed-log PDF numerator only while it retains at least this
// (log-scale) fraction of its largest term; past that, cancellation has
// consumed the tail-log accuracy and the structural fallback is safer.
constexpr double BAWL_LOG_BRACKET_MIN = -13.815510557964274;  // log(1e-6)
// Normalizer floors match each model's legacy pmax(pnorm(v/sv), floor).
constexpr double LBA_DENOM_FLOOR = 1e-10;
// BAwL's 1e-300 value is a historical model contract, not a requirement for
// the log-space arithmetic: the log path can represent much smaller normal
// probabilities.  Keep it here because removing it changes BAwL values in
// the z < about -37 regime; the exact k = 0 LBA entry points intentionally
// continue to use LBA_DENOM_FLOOR.
constexpr double BAWL_DENOM_FLOOR = 1e-300;

// Launch-strength distribution selector (ContextForRaceModels::bawl_launch and
// the `launch` argument of the exported d/p functions -- keep them in sync).
// Deliberately the same values as BAWD_LAUNCH_* in model_BAwD.h, which is
// included at the bottom of this file and so cannot supply them here.
constexpr int BAWL_LAUNCH_NORMAL = 0;
constexpr int BAWL_LAUNCH_LOGNORMAL = 1;
constexpr int BAWL_LAUNCH_SPLITLOGNORMAL = 2;
// Acceptance modes for the guarded natural-space evaluators.
//   STRICT: value must stand on its own (feeds std::log directly).
//   RAW:    underflow to 0 is fine (caller floors at min_ll), but near-1
//           CDFs are rejected so log1p(-cdf) survivors keep their tail.
//   CLAMP:  as RAW, and near-1/saturated CDFs are clamped to 1 -- for
//           consumers that clamp to [0, 1] anyway (truncation normalisers,
//           GSL integrands), where the survivor tail is immaterial.
constexpr int BA_ACCEPT_STRICT = 0;
constexpr int BA_ACCEPT_RAW = 1;
constexpr int BA_ACCEPT_CLAMP = 2;

inline double log_positive_normalizer(double v, double sv, bool posdrift,
                                     double denom_floor = LBA_DENOM_FLOOR) {
  if (!posdrift) return 0.0;

  // Match the selected model's pmax(pnorm(v / sv), floor) normalization,
  // but apply the floor directly to its logarithm.
  const double log_denom = pnorm_log_direct(v / sv, true);
  const double log_floor = std::log(denom_floor);
  return (log_denom < log_floor) ? log_floor : log_denom;
}

// Log CDF of an unnormalised positive normal launch.  The lower-tail term at
// V = 0 is essential: dividing Phi((w-v)/sv) by Phi(v/sv) alone is the CDF of
// an *untruncated* normal with a rescaled mass, not the CDF conditional on V>0.
// Keeping this subtraction in log space is what makes the survivor paths
// correct for ordinary positive-drift fits as well as for far tails.
inline double log_normal_cdf_positive_raw(double w, double v, double sv) {
  if (!(sv > 0.0) || !(w > 0.0)) return R_NegInf;
  if (!emc2_isfinite(w)) return pnorm_log_direct(v / sv, true);
  const double lhi = pnorm_log_direct((w - v) / sv, true);
  const double llo = pnorm_log_direct(-v / sv, true);
  return log_diff_exp(lhi, llo);
}

inline double log_normal_cdf_positive(double w, double v, double sv,
                                      bool posdrift, double denom_floor) {
  if (!posdrift) {
    if (!(sv > 0.0) || ISNAN(w)) return R_NegInf;
    return pnorm_log_direct((w - v) / sv, true);
  }
  const double raw = log_normal_cdf_positive_raw(w, v, sv);
  return (raw > R_NegInf)
    ? raw - log_positive_normalizer(v, sv, true, denom_floor) : R_NegInf;
}

// Log integral of the same positive-normal CDF over a standardized interval.
// The direct difference is well conditioned in the central region; when both
// terms coincide at floating-point precision, the midpoint CDF gives the
// correct narrow-interval limit and avoids manufacturing a zero survivor.
inline double log_normal_phi_integral_positive_raw(double lo, double hi,
                                                   double v, double sv,
                                                   bool posdrift) {
  if (!(hi > lo)) return R_NegInf;
  if (!posdrift) return log_normal_phi_integral(lo, hi);
  const double base = pnorm_log_direct(-v / sv, true);
  const double integral = log_normal_phi_integral(lo, hi);
  const double sub = std::log(hi - lo) + base;
  const double out = log_diff_exp(integral, sub);
  if (out > R_NegInf) return out;
  const double mid = 0.5 * (lo + hi);
  const double point = log_diff_exp(pnorm_log_direct(mid, true), base);
  return point > R_NegInf ? std::log(hi - lo) + point : R_NegInf;
}

// Natural formulas are much cheaper, but are only used while their normal
// endpoints are central and the relevant subtraction is well-conditioned.
// The log-space branches remain authoritative for tails and cancellation.
inline bool natural_normalizer(double v, double sv, bool posdrift,
                               double &denom,
                               double denom_floor = LBA_DENOM_FLOOR) {
  if (!(sv > 0.0)) return false;
  if (!posdrift) {
    denom = 1.0;
    return true;
  }
  const double z = v / sv;
  denom = pnorm_std(z, true, false);
  // Below the selected model's floor, preserve the normalizer in log space.
  return denom > denom_floor;
}

inline bool natural_normal_interval_safe(double lo, double hi) {
  return hi > lo && std::fabs(lo) <= BAWL_NATURAL_Z_MAX &&
    std::fabs(hi) <= BAWL_NATURAL_Z_MAX;
}

// The CDF is often consumed as log(1 - F), so values too close to either
// endpoint must stay on the stable log path even when the probability itself
// is representable.  In particular, 1 - mean{Phi(z)} loses the low tail.
inline bool natural_cdf_safe(double p) {
  constexpr double cdf_margin = 1e-8;
  return p > cdf_margin && p < 1.0 - cdf_margin;
}

// --------------------------------------------------------------------------
// Shared BAwL core (standard LBA is the exact k = 0 member).
//
// For t > 0, define
//   c = (v - k*b/(1-exp(-k*t))) / sv,
//   m = k*exp(-k*t) / ((1-exp(-k*t))*sv).
// The hit probability conditional on start a ~ Unif(0, A) is Phi(c + m*a),
// so with span = m*A:
//   F(t) = integral_c^{c+span} Phi(z) dz / (span * denom),
//   f(t) = jacobian * [(b*m + c)*(Phi(c+span)-Phi(c)) + phi(c+span)-phi(c)]
//          / (A * sv * m^2 * denom).
// k = 0 uses its exact limit c = (v - b/t)/sv, m = 1/(t*sv), jacobian =
// 1/t^2, which recovers the standard LBA formulas identically.
//
// Each quantity pairs a guarded natural-space evaluation (fast; valid while
// the normal endpoints are central and subtractions are well-conditioned)
// with an authoritative log-space evaluation for tails and cancellation.
// --------------------------------------------------------------------------

inline void bawl_leak_factors(double kt, double &E, double &G) {
  // expm1 is needed only in the cancellation-prone small-kt regime.  Away
  // from it, derive G from the single exponential evaluation.
  if (kt < 1e-4) {
    G = -std::expm1(-kt);
    E = 1.0 - G;
  } else {
    E = std::exp(-kt);
    G = 1.0 - E;
  }
}

inline void bawl_threshold_terms(double t, double b, double sv, double k,
                                 double &c, double &m,
                                 double *log_jacobian = nullptr,
                                 double *jacobian = nullptr) {
  if (k <= BAWL_K_EPS) {
    // Exact k -> 0 limit of the general expressions below.
    c = -(b / t) / sv;
    m = (1.0 / t) / sv;
    if (log_jacobian != nullptr) *log_jacobian = -2.0 * std::log(t);
    if (jacobian != nullptr) *jacobian = 1.0 / (t * t);
    return;
  }

  const double kt = k * t;
  double E, G;
  bawl_leak_factors(kt, E, G);
  const double G_safe = clamp_pos(G, 1e-300);
  const double C1 = (k * b) / G_safe;
  const double C2 = (k * E) / G_safe;
  c = -C1 / sv;
  m = C2 / sv;
  if (log_jacobian != nullptr) {
    *log_jacobian = 2.0 * std::log(k) - kt - 2.0 * std::log(G_safe);
  }
  if (jacobian != nullptr) {
    *jacobian = (k * k * E) / (G_safe * G_safe);
  }
}

// Guarded natural-space CDF.  Returns true when `cdf` holds a value the
// caller may use under the given acceptance mode (see BA_ACCEPT_*): STRICT
// requires central normal endpoints and a value away from both 0 and 1;
// RAW additionally admits saturated tails that round down to 0; CLAMP also
// admits (and clamps) values at or near 1.
inline bool ba_natural_cdf(double t, double A, double b, double v, double sv,
                           double k, bool posdrift, double denom_floor,
                           int accept_mode, double &cdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  const auto accept = [accept_mode](double &p) {
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(p);
    if (accept_mode == BA_ACCEPT_RAW) return p < 1.0 - 1e-8;
    if (p > 1.0) p = 1.0;
    return true;
  };
  if (!(t > 0.0) || !(sv > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) {
    cdf = 0.0;
    return lenient;  // strict callers take the exact -Inf log path
  }
  if (t == R_PosInf && k <= BAWL_K_EPS) {
    cdf = posdrift ? 1.0 : pnorm_std(v / sv, true, false);
    if (!R_FINITE(cdf)) return false;
    return accept(cdf);
  }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;
  double c, m;
  bawl_threshold_terms(t, b, sv, k, c, m);
  c += v / sv;
  if (!(m > 0.0) || !emc2_isfinite(m)) return false;  // t = Inf with k > 0
  const double span = m * A;

  if (A <= BAWL_A_EPS || span < BAWL_NATURAL_MIN_SPAN) {
    // Point-mass start (a = 0 limit) or collapsed start range (midpoint).
    const double z = (A <= BAWL_A_EPS) ? c : c + 0.5 * span;
    if (!emc2_isfinite(z)) return false;
    if (!lenient && std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(z, true, false) / denom;
  } else {
    const double z_hi = c + span;
    if (!emc2_isfinite(c) || !emc2_isfinite(z_hi)) return false;
    if (!lenient && !natural_normal_interval_safe(c, z_hi)) return false;
    // H(z) = z Phi(z) + phi(z), so H(hi) - H(lo) integrates Phi.
    const double H_hi = z_hi * pnorm_std(z_hi, true, false) + dnormP(z_hi);
    const double H_lo = c * pnorm_std(c, true, false) + dnormP(c);
    const double integral = H_hi - H_lo;
    if (!(integral > 0.0)) {
      // Saturated or fully cancelled: mass below natural resolution.
      if (!lenient) return false;
      // With a plainly saturated interval the CDF is 1, not 0.
      cdf = (integral == 0.0 && c > BAWL_NATURAL_Z_MAX) ? 1.0 : 0.0;
      if (cdf == 1.0 && accept_mode != BA_ACCEPT_CLAMP) return false;
      return true;
    }
    if (!lenient) {
      const double scale = std::fabs(H_hi) + std::fabs(H_lo);
      if (integral <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;
    }
    cdf = integral / (span * denom);
  }

  if (!R_FINITE(cdf)) return false;
  if (cdf <= 0.0) {
    if (!lenient) return false;
    cdf = 0.0;
    return true;
  }
  return accept(cdf);
}

// Guarded natural-space PDF; acceptance semantics as for ba_natural_cdf.
// Near-1 has no meaning for a density, so RAW and CLAMP coincide; genuine
// positive cancellation always defers to the log path.
inline bool ba_natural_pdf(double t, double A, double b, double v, double sv,
                           double k, bool posdrift, double denom_floor,
                           int accept_mode, double &pdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!(t > 0.0) || t == R_PosInf || !(sv > 0.0) || !(A >= 0.0) ||
      !(b >= A) || !(b > 0.0)) {
    pdf = 0.0;
    return lenient;
  }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;
  double c, m, jacobian;
  bawl_threshold_terms(t, b, sv, k, c, m, nullptr, &jacobian);
  c += v / sv;
  if (!(m > 0.0) || !emc2_isfinite(m)) return false;
  if (!(jacobian > 0.0) || !(jacobian < R_PosInf)) return false;
  const double span = m * A;

  if (A <= BAWL_A_EPS || span < BAWL_NATURAL_MIN_SPAN) {
    const double z = (A <= BAWL_A_EPS) ? c : c + 0.5 * span;
    const double eff_b = (A <= BAWL_A_EPS) ? b : b - 0.5 * A;
    if (!emc2_isfinite(z)) return false;
    if (!lenient && std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    pdf = jacobian * eff_b * dnormP(z) / (sv * denom);
  } else {
    const double z_hi = c + span;
    if (!emc2_isfinite(c) || !emc2_isfinite(z_hi)) return false;
    if (!lenient && !natural_normal_interval_safe(c, z_hi)) return false;
    const double term1 = (b * m + c) *
      (pnorm_std(z_hi, true, false) - pnorm_std(c, true, false));
    const double term2 = dnormP(z_hi) - dnormP(c);
    const double bracket = term1 + term2;
    if (!(bracket > 0.0)) {
      if (!lenient) return false;
      pdf = 0.0;
      return true;
    }
    const double scale = std::fabs(term1) + std::fabs(term2);
    if (bracket <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;
    const double denom_scale = A * sv * m * m * denom;
    if (!(denom_scale > 0.0) || !(denom_scale < R_PosInf)) return false;
    pdf = jacobian * bracket / denom_scale;
  }

  if (!R_FINITE(pdf)) return false;
  if (pdf <= 0.0) {
    if (!lenient) return false;
    pdf = 0.0;
  }
  return true;
}

// Authoritative log-space CDF.
inline double log_ba_cdf(double t, double A, double b, double v, double sv,
                         double k, bool posdrift, double denom_floor) {
  if (!(t > 0.0) || !(sv > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0))
    return R_NegInf;
  if (t == R_PosInf && k <= BAWL_K_EPS)
    return posdrift ? 0.0 : pnorm_log_direct(v / sv, true);

  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  double c, m;
  bawl_threshold_terms(t, b, sv, k, c, m);
  c += v / sv;
  if (!(c > R_NegInf)) return R_NegInf;  // catches -Inf and NaN

  // Point-mass start; t = Inf with k > 0 lands here (m = 0) and keeps the
  // defective upper tail Phi((v - k*b)/sv).
  if (A <= BAWL_A_EPS || !(m > 0.0))
    return pnorm_log_direct(c, true) - log_denom;

  const double span = m * A;
  if (span > BAWL_LOG_MIN_SPAN && emc2_isfinite(span)) {
    const double log_integral = log_normal_phi_integral(c, c + span);
    const double out = log_integral - std::log(span) - log_denom;
    // NaN: both endpoints so far in one tail that the antiderivative
    // difference lost all digits; the midpoint below is the stable limit.
    if (!ISNAN(out)) return out >= 0.0 ? 0.0 : out;
  }
  const double out = pnorm_log_direct(c + 0.5 * span, true) - log_denom;
  return out >= 0.0 ? 0.0 : out;
}

// Authoritative log-space PDF.
inline double log_ba_pdf(double t, double A, double b, double v, double sv,
                         double k, bool posdrift, double denom_floor) {
  if (!(t > 0.0) || t == R_PosInf || !(sv > 0.0) || !(A >= 0.0) ||
      !(b >= A) || !(b > 0.0))
    return R_NegInf;

  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  double c, m, log_jacobian;
  bawl_threshold_terms(t, b, sv, k, c, m, &log_jacobian);
  c += v / sv;
  if (!(c > R_NegInf) || !(log_jacobian > R_NegInf)) return R_NegInf;

  if (A <= BAWL_A_EPS || !(m > 0.0))
    return log_jacobian + std::log(b) + dnormP(c, 0.0, 1.0, true) -
      std::log(sv) - log_denom;

  const double span = m * A;
  if (span > BAWL_LOG_MIN_SPAN && emc2_isfinite(span)) {
    const double hi = c + span;
    const double log_dphi = log_normal_interval(c, hi);
    if (!ISNAN(log_dphi)) {
      // Density numerator: (b*m + c) * DeltaPhi + Delta(phi), signed.
      const signed_log lead = signed_log_product(b * m + c, log_dphi);
      const signed_log phi_diff = signed_log_sub(
        make_signed_log(dnormP(hi, 0.0, 1.0, true), 1),
        make_signed_log(dnormP(c, 0.0, 1.0, true), 1));
      const signed_log bracket = signed_log_add(lead, phi_diff);
      const double max_term = std::max(lead.log_abs, phi_diff.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWL_LOG_BRACKET_MIN) {
        return log_jacobian - std::log(A) - std::log(sv) - log_denom -
          2.0 * std::log(m) + bracket.log_abs;
      }
      // The two terms cancelled past the accuracy of the tail logs (a sign
      // flip is always numerical: the true numerator is positive).  It is
      // m^2 * integral_0^A (b - a) phi(c + m a) da, with the factor (b - a)
      // confined to [b - A, b], so the midpoint-in-a form
      // m * (b - A/2) * DeltaPhi is exact in DeltaPhi and bounds the
      // weighting error by roughly a factor of two in the worst case.
      return log_jacobian - std::log(A) - std::log(sv) - log_denom -
        std::log(m) + std::log(b - 0.5 * A) + log_dphi;
    }
  }
  // Collapsed start range (or unusable DeltaPhi): midpoint limit in z.
  return log_jacobian + std::log(b - 0.5 * A) +
    dnormP(c + 0.5 * span, 0.0, 1.0, true) - std::log(sv) - log_denom;
}

// --------------------------------------------------------------------------
// Lognormal launch strength (log V ~ N(mu, sigma^2)).
//
// The crossing condition is the same one the Gaussian branch above uses, read
// in launch-strength rather than z units: with E = exp(-k t), G = 1 - E, an
// accumulator started at a ~ Unif(0, A) has crossed b by t iff
//   V >= W(a) = k (b - a E) / G,
// which is AFFINE in a.  Substituting w = W(a) therefore turns both integrals
// over (V, a) into integrals over w alone, whatever the launch distribution:
//   F(t) = (s / A) [C(w_lo) - C(w_hi)],      C(x) = E[(V - x)_+]
//   f(t) = (1 / (A E)) integral_{w_lo}^{w_hi} (w - k b) g(w) dw
// with s = G / (k E) (limit t at k = 0), w_hi = W(0), w_lo = W(A), and the
// weight (w - k b) = k E (b - a) / G non-negative because b >= A.
//
// For the lognormal launch C is the stop-loss price log_lognormal_stoploss()
// and the density integral is a partial expectation minus a probability, so
// every piece is closed form -- the same three primitives BAwD's lognormal
// branch uses.  The Gaussian branch keeps its own z-space formulation because
// its (b m + c) grouping is what makes its cancellation guard tractable.
// --------------------------------------------------------------------------

struct BawlLaunchGeom {
  double w_hi = 0.0;    // W(0) = k b / G      (b / t at k = 0)
  double w_lo = 0.0;    // W(A)
  double width = 0.0;   // w_hi - w_lo = A / s, formed without cancellation
  double s = 0.0;       // G / (k E); its k = 0 limit is t
  double log_s = R_NegInf;
  double E = 1.0;       // exp(-k t)
  double log_E = 0.0;   // -k t, exact once E itself has underflowed
  bool frozen = false;  // t = Inf: the live interval has collapsed to a point
  bool ok = false;
};

// log_jacobian of the Gaussian branch, expressed in the same quantities:
//   2 log k - k t - 2 log G = -2 log s - log E.
inline double bawl_geom_log_jacobian(const BawlLaunchGeom &g) {
  return -2.0 * g.log_s - g.log_E;
}

inline BawlLaunchGeom bawl_launch_geom(double t, double A, double b, double k) {
  BawlLaunchGeom g;
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) return g;

  if (k <= BAWL_K_EPS) {
    if (t == R_PosInf) {
      // Proper k = 0 limit: every positive launch strength eventually crosses.
      g.frozen = true;
      g.s = R_PosInf;
      g.log_s = R_PosInf;
      g.ok = true;
      return g;
    }
    g.s = t;
    g.log_s = std::log(t);
    g.w_hi = b / t;
    g.width = A / t;
    g.w_lo = (b - A) / t;
    g.ok = true;
    return g;
  }

  const double kt = k * t;
  if (!(kt < R_PosInf)) {
    // t = Inf with leak: only launches above k b ever finish, and the start
    // point no longer matters.  This is the defective upper tail.
    g.frozen = true;
    g.w_hi = k * b;
    g.w_lo = g.w_hi;
    g.E = 0.0;
    g.log_E = R_NegInf;
    g.s = R_PosInf;
    g.log_s = R_PosInf;
    g.ok = true;
    return g;
  }

  double E, G;
  bawl_leak_factors(kt, E, G);
  const double G_safe = clamp_pos(G, 1e-300);
  g.E = E;
  g.log_E = -kt;
  g.w_hi = (k * b) / G_safe;
  // b - A E carries no cancellation (b >= A >= A E), so both endpoints and the
  // width are formed directly rather than by differencing w_hi.
  g.w_lo = (k * (b - A * E)) / G_safe;
  g.width = (k * A * E) / G_safe;
  g.log_s = std::log(G_safe) - std::log(k) + kt;
  g.s = (E > 0.0) ? G_safe / (k * E) : R_PosInf;
  g.ok = emc2_isfinite(g.w_hi) && emc2_isfinite(g.w_lo);
  return g;
}

// Direct launch-CDF integral for the BAwL survivor.  Unlike complementing the
// CDF, this remains accurate when all start points have essentially crossed.
inline double log_bawl_surv_normal(double t, double A, double b, double v,
                                   double sv, double k, bool posdrift,
                                   double denom_floor) {
  if (!(t > 0.0) || !(sv > 0.0)) return R_NegInf;
  const BawlLaunchGeom g = bawl_launch_geom(t, A, b, k);
  if (!g.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  if (g.frozen || A <= BAWL_A_EPS || !(g.width > 0.0)) {
    const double w = g.frozen ? k * b : g.w_hi;
    return std::fmin(log_normal_cdf_positive(w, v, sv, posdrift, denom_floor), 0.0);
  }
  const double li = log_normal_phi_integral_positive_raw(
    (g.w_lo - v) / sv, (g.w_hi - v) / sv, v, sv, posdrift);
  if (!(li > R_NegInf)) return R_NegInf;
  const double out = g.log_s - std::log(A) + std::log(sv) + li - log_denom;
  return std::fmin(out, 0.0);
}

inline double log_bawl_surv_logn(double t, double A, double b, double mu,
                                 double sigma, double k, double delta) {
  if (!(t > 0.0) || !(sigma > 0.0)) return R_NegInf;
  const BawlLaunchGeom g = bawl_launch_geom(t, A, b, k);
  if (!g.ok) return R_NegInf;
  if (g.frozen || A <= BAWL_A_EPS || !(g.width > 0.0)) {
    const double w = g.frozen ? k * b : g.w_hi;
    return delta == 0.0
      ? std::fmin(pnorm_log_direct((std::log(w) - mu) / sigma, true), 0.0)
      : std::fmin(log_split_lognormal_cdf(w, mu, sigma, delta), 0.0);
  }
  const double pa = delta == 0.0
    ? log_lognormal_put(g.w_hi, mu, sigma)
    : log_split_lognormal_put(g.w_hi, mu, sigma, delta);
  const double pb = delta == 0.0
    ? log_lognormal_put(g.w_lo, mu, sigma)
    : log_split_lognormal_put(g.w_lo, mu, sigma, delta);
  if (!(pa - pb > BAWL_LOG_MIN_SPAN)) {
    const double w_mid = 0.5 * (g.w_hi + g.w_lo);
    return delta == 0.0
      ? pnorm_log_direct((std::log(w_mid) - mu) / sigma, true)
      : log_split_lognormal_cdf(w_mid, mu, sigma, delta);
  }
  const double out = g.log_s - std::log(A) + log_diff_exp(pa, pb);
  return std::fmin(out, 0.0);
}

// log P(V >= w) for the lognormal launch.
inline double bawl_logn_log_surv(double w, double mu, double sigma,
                                 double delta = 0.0) {
  if (!(w > 0.0)) return 0.0;
  if (!emc2_isfinite(w)) return R_NegInf;
  return delta == 0.0
    ? pnorm_log_direct((mu - std::log(w)) / sigma, true)
    : log_split_lognormal_survivor(w, mu, sigma, delta);
}

inline double log_bawl_cdf_logn(double t, double A, double b, double mu,
                                double sigma, double k, double delta) {
  if (!(sigma > 0.0)) return R_NegInf;
  const BawlLaunchGeom g = bawl_launch_geom(t, A, b, k);
  if (!g.ok) return R_NegInf;
  if (g.frozen)
    return (k <= BAWL_K_EPS) ? 0.0 :
      std::fmin(bawl_logn_log_surv(g.w_hi, mu, sigma, delta), 0.0);

  const double w_mid = 0.5 * (g.w_lo + g.w_hi);
  if (A <= BAWL_A_EPS || !(g.width > 0.0))
    return std::fmin(bawl_logn_log_surv(g.w_hi, mu, sigma, delta), 0.0);

  // int_{w_lo}^{w_hi} P(V >= w) dw = C(w_lo) - C(w_hi).
  const double lc_lo = delta == 0.0
    ? log_lognormal_stoploss(g.w_lo, mu, sigma)
    : log_split_lognormal_stoploss(g.w_lo, mu, sigma, delta);
  const double lc_hi = delta == 0.0
    ? log_lognormal_stoploss(g.w_hi, mu, sigma)
    : log_split_lognormal_stoploss(g.w_hi, mu, sigma, delta);
  if (lc_lo - lc_hi > 1e-8) {
    const double ld = log_diff_exp(lc_lo, lc_hi);
    if (!ISNAN(ld) && ld > R_NegInf) {
      const double out = g.log_s - std::log(A) + ld;
      if (!ISNAN(out)) return std::fmin(out, 0.0);
    }
  }
  // Interval too narrow (or endpoint-differenced away): (s/A) times the width
  // is exactly one, so the midpoint survivor is the stable limit.
  return std::fmin(bawl_logn_log_surv(w_mid, mu, sigma, delta), 0.0);
}

inline double log_bawl_pdf_logn(double t, double A, double b, double mu,
                                double sigma, double k, double delta) {
  if (!(sigma > 0.0) || t == R_PosInf) return R_NegInf;
  const BawlLaunchGeom g = bawl_launch_geom(t, A, b, k);
  if (!g.ok || g.frozen) return R_NegInf;
  const double log_jac = bawl_geom_log_jacobian(g);

  if (A <= BAWL_A_EPS || !(g.width > 0.0)) {
    const double lg = delta == 0.0
      ? dlnorm_std(g.w_hi, mu, sigma, true)
      : log_split_lognormal_density(g.w_hi, mu, sigma, delta);
    if (!(lg > R_NegInf)) return R_NegInf;
    return log_jac + std::log(b) + lg;
  }

  if (delta == 0.0) {
    // Numerator = E[V 1{w_lo < V < w_hi}] - k b P(w_lo < V < w_hi).
    const double d1_lo = (mu + sigma * sigma - std::log(g.w_lo)) / sigma;
    const double d1_hi = (mu + sigma * sigma - std::log(g.w_hi)) / sigma;
    const double d2_lo = (mu - std::log(g.w_lo)) / sigma;
    const double d2_hi = (mu - std::log(g.w_hi)) / sigma;
    const double log_prob = log_normal_interval(d2_hi, d2_lo);

    if (!ISNAN(log_prob)) {
      const double log_pexp = log_normal_interval(d1_hi, d1_lo);
      if (!ISNAN(log_pexp)) {
        const signed_log t1 =
          make_signed_log(mu + 0.5 * sigma * sigma + log_pexp, 1);
        const signed_log t2 = signed_log_product(k * b, log_prob);
        const signed_log bracket = signed_log_sub(t1, t2);
        const double max_term = std::fmax(t1.log_abs, t2.log_abs);
        if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
            bracket.log_abs > max_term + BAWL_LOG_BRACKET_MIN) {
          return bracket.log_abs - std::log(A) - g.log_E;
        }
      }
      if (log_prob > R_NegInf)
        return std::log(b - 0.5 * A) + log_prob - std::log(A) - g.log_s -
          g.log_E;
    }
    const double lg = dlnorm_std(0.5 * (g.w_lo + g.w_hi), mu, sigma, true);
    if (!(lg > R_NegInf)) return R_NegInf;
    return log_jac + std::log(b - 0.5 * A) + lg;
  }

  const double lp_lo = log_split_lognormal_survivor(g.w_lo, mu, sigma, delta);
  const double lp_hi = log_split_lognormal_survivor(g.w_hi, mu, sigma, delta);
  const double lm_lo = log_split_lognormal_first_partial_moment(
    g.w_lo, mu, sigma, delta);
  const double lm_hi = log_split_lognormal_first_partial_moment(
    g.w_hi, mu, sigma, delta);
  const double log_prob = log_diff_exp(lp_lo, lp_hi);
  if (log_prob > R_NegInf) {
    const double log_pexp = log_diff_exp(lm_lo, lm_hi);
    if (log_pexp > R_NegInf) {
      const signed_log t1 = make_signed_log(log_pexp, 1);
      const signed_log t2 = signed_log_product(k * b, log_prob);
      const signed_log bracket = signed_log_sub(t1, t2);
      const double max_term = std::fmax(t1.log_abs, t2.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWL_LOG_BRACKET_MIN) {
        return bracket.log_abs - std::log(A) - g.log_E;
      }
    }
    return std::log(b - 0.5 * A) + log_prob - std::log(A) - g.log_s -
      g.log_E;
  }
  const double lg = log_split_lognormal_density(
    0.5 * (g.w_lo + g.w_hi), mu, sigma, delta);
  if (!(lg > R_NegInf)) return R_NegInf;
  return log_jac + std::log(b - 0.5 * A) + lg;
}
inline bool bawl_natural_cdf_logn(double t, double A, double b, double mu,
                                  double sigma, double k, int accept_mode,
                                  double &cdf, double delta = 0.0) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (delta != 0.0) {
    const double lp = log_bawl_cdf_logn(t, A, b, mu, sigma, k, delta);
    cdf = (lp > R_NegInf) ? std::exp(lp) : 0.0;
    if (!R_FINITE(cdf)) return false;
    if (cdf <= 0.0) return lenient;
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(cdf);
    if (accept_mode == BA_ACCEPT_RAW) return cdf < 1.0 - 1e-8;
    if (cdf > 1.0) cdf = 1.0;
    return true;
  }
  const auto accept = [accept_mode](double &p) {
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(p);
    if (accept_mode == BA_ACCEPT_RAW) return p < 1.0 - 1e-8;
    if (p > 1.0) p = 1.0;
    return true;
  };
  if (!(sigma > 0.0)) return false;
  const BawlLaunchGeom g = bawl_launch_geom(t, A, b, k);
  if (!g.ok) {
    cdf = 0.0;
    return lenient;
  }
  if (g.frozen) {
    if (k <= BAWL_K_EPS) {
      cdf = 1.0;
      return accept_mode == BA_ACCEPT_CLAMP;
    }
    const double z = (mu - std::log(g.w_hi)) / sigma;
    if (!lenient && std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(z, true, false);
    return accept(cdf);
  }

  if (A <= BAWL_A_EPS || g.width < BAWL_NATURAL_MIN_SPAN * g.w_hi) {
    const double w = (A <= BAWL_A_EPS) ? g.w_hi : 0.5 * (g.w_lo + g.w_hi);
    if (!(w > 0.0)) return false;
    const double z = (mu - std::log(w)) / sigma;
    if (!emc2_isfinite(z)) return false;
    if (!lenient && std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(z, true, false);
  } else {
    double scale = 0.0;
    const double diff =
      lognormal_stoploss_interval_nat(g.w_lo, g.w_hi, mu, sigma, scale);
    if (!(diff > 0.0)) {
      if (!lenient) return false;
      cdf = 0.0;
      return true;
    }
    if (!lenient && diff <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale))
      return false;
    cdf = (g.s / A) * diff;
  }

  if (!R_FINITE(cdf)) return false;
  if (cdf <= 0.0) {
    if (!lenient) return false;
    cdf = 0.0;
    return true;
  }
  return accept(cdf);
}

inline bool bawl_natural_pdf_logn(double t, double A, double b, double mu,
                                  double sigma, double k, int accept_mode,
                                  double &pdf, double delta = 0.0) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (delta != 0.0) {
    const double lp = log_bawl_pdf_logn(t, A, b, mu, sigma, k, delta);
    pdf = (lp > R_NegInf) ? std::exp(lp) : 0.0;
    if (!R_FINITE(pdf)) return false;
    if (pdf <= 0.0 && !lenient) return false;
    return true;
  }
  if (!(sigma > 0.0) || t == R_PosInf) return false;
  const BawlLaunchGeom g = bawl_launch_geom(t, A, b, k);
  if (!g.ok || g.frozen) {
    pdf = 0.0;
    return lenient;
  }

  if (A <= BAWL_A_EPS || g.width < BAWL_NATURAL_MIN_SPAN * g.w_hi) {
    const double w = (A <= BAWL_A_EPS) ? g.w_hi : 0.5 * (g.w_lo + g.w_hi);
    const double eff_b = (A <= BAWL_A_EPS) ? b : b - 0.5 * A;
    if (!(w > 0.0)) return false;
    const double z = (mu - std::log(w)) / sigma;
    if (!lenient && std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    const double log_pdf = bawl_geom_log_jacobian(g) + std::log(eff_b) +
      dlnorm_std(w, mu, sigma, true);
    if (!emc2_isfinite(log_pdf)) return false;
    pdf = std::exp(log_pdf);
  } else {
    if (!(g.E > 0.0)) return false;
    const double M = std::exp(mu + 0.5 * sigma * sigma);
    if (!R_FINITE(M)) return false;
    const double d1_lo = (mu + sigma * sigma - std::log(g.w_lo)) / sigma;
    const double d1_hi = (mu + sigma * sigma - std::log(g.w_hi)) / sigma;
    const double d2_lo = (mu - std::log(g.w_lo)) / sigma;
    const double d2_hi = (mu - std::log(g.w_hi)) / sigma;
    if (!lenient && (!natural_normal_interval_safe(d1_hi, d1_lo) ||
                     !natural_normal_interval_safe(d2_hi, d2_lo)))
      return false;
    const double P1 = pnorm_std(d1_lo, true, false) - pnorm_std(d1_hi, true, false);
    const double P0 = pnorm_std(d2_lo, true, false) - pnorm_std(d2_hi, true, false);
    const double term1 = M * P1;
    const double term2 = k * b * P0;
    const double num = term1 - term2;
    if (!(num > 0.0)) {
      if (!lenient) return false;
      pdf = 0.0;
      return true;
    }
    const double scale = std::fabs(term1) + std::fabs(term2);
    if (num <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;
    pdf = num / (A * g.E);
  }

  if (!R_FINITE(pdf)) return false;
  if (pdf <= 0.0) {
    if (!lenient) return false;
    pdf = 0.0;
  }
  return true;
}

// Launch-distribution dispatch.  `p1`/`p2` are (v, sv) for the normal launch
// and (mu, sigma) for the lognormal one; they occupy the same kernel columns,
// exactly as in BAwD.  posdrift and the normalizer floor are meaningless for
// the lognormal launch (V > 0 by construction, so nothing is truncated).
inline bool ba_natural_cdf_launch(double t, double A, double b, double p1,
                                  double p2, double k, bool posdrift,
                                  double denom_floor, int accept_mode,
                                  double &cdf, int launch,
                                  double delta = 0.0) {
  if (launch == BAWL_LAUNCH_LOGNORMAL ||
      launch == BAWL_LAUNCH_SPLITLOGNORMAL)
    return bawl_natural_cdf_logn(t, A, b, p1, p2, k, accept_mode, cdf, delta);
  return ba_natural_cdf(t, A, b, p1, p2, k, posdrift, denom_floor,
                        accept_mode, cdf);
}

inline bool ba_natural_pdf_launch(double t, double A, double b, double p1,
                                  double p2, double k, bool posdrift,
                                  double denom_floor, int accept_mode,
                                  double &pdf, int launch,
                                  double delta = 0.0) {
  if (launch == BAWL_LAUNCH_LOGNORMAL ||
      launch == BAWL_LAUNCH_SPLITLOGNORMAL)
    return bawl_natural_pdf_logn(t, A, b, p1, p2, k, accept_mode, pdf, delta);
  return ba_natural_pdf(t, A, b, p1, p2, k, posdrift, denom_floor,
                        accept_mode, pdf);
}

inline double log_ba_cdf_launch(double t, double A, double b, double p1,
                                double p2, double k, bool posdrift,
                                double denom_floor, int launch,
                                double delta = 0.0) {
  if (launch == BAWL_LAUNCH_LOGNORMAL ||
      launch == BAWL_LAUNCH_SPLITLOGNORMAL)
    return log_bawl_cdf_logn(t, A, b, p1, p2, k, delta);
  return log_ba_cdf(t, A, b, p1, p2, k, posdrift, denom_floor);
}

inline double log_ba_surv_launch(double t, double A, double b, double p1,
                                 double p2, double k, bool posdrift,
                                 double denom_floor, int launch,
                                 double delta = 0.0) {
  if (launch == BAWL_LAUNCH_LOGNORMAL ||
      launch == BAWL_LAUNCH_SPLITLOGNORMAL)
    return log_bawl_surv_logn(t, A, b, p1, p2, k, delta);
  return log_bawl_surv_normal(t, A, b, p1, p2, k, posdrift, denom_floor);
}

inline double log_ba_pdf_launch(double t, double A, double b, double p1,
                                double p2, double k, bool posdrift,
                                double denom_floor, int launch,
                                double delta = 0.0) {
  if (launch == BAWL_LAUNCH_LOGNORMAL ||
      launch == BAWL_LAUNCH_SPLITLOGNORMAL)
    return log_bawl_pdf_logn(t, A, b, p1, p2, k, delta);
  return log_ba_pdf(t, A, b, p1, p2, k, posdrift, denom_floor);
}
inline double bawl_cdf_norm(double t, double A, double b, double v,
                            double sv, double k, bool posdrift, bool log_out,
                            double denom_floor = BAWL_DENOM_FLOOR,
                            int launch = BAWL_LAUNCH_NORMAL,
                            double delta = 0.0) {
  // The k = 0 process is proper whenever the launch strength is positive by
  // construction (posdrift, or a lognormal launch).  Its CDF at +Inf is
  // exactly one, so do not send the upper truncation bound through the strict
  // natural guard and then back through log space.  This is deliberately in
  // the shared BAwL wrapper: k = 0 can arrive here through pleakyba_norm(),
  // not only through the LBA() adapter.  Killed-clock callers use their own
  // deadline and upper-tail handling.
  if (t == R_PosInf && std::fabs(k) <= BAWL_K_EPS &&
      (posdrift || launch == BAWL_LAUNCH_LOGNORMAL ||
       launch == BAWL_LAUNCH_SPLITLOGNORMAL) && sv > 0.0 &&
      A >= 0.0 && b >= A && b > 0.0)
    return log_out ? 0.0 : 1.0;

  double cdf;
  if (ba_natural_cdf_launch(t, A, b, v, sv, k, posdrift, denom_floor,
                            BA_ACCEPT_STRICT, cdf, launch, delta))
    return log_out ? std::log(cdf) : cdf;
  return return_from_log(
    log_ba_cdf_launch(t, A, b, v, sv, k, posdrift, denom_floor, launch, delta),
    log_out);
}

inline double bawl_pdf_norm(double t, double A, double b, double v,
                            double sv, double k, bool posdrift, bool log_out,
                            double denom_floor = BAWL_DENOM_FLOOR,
                            int launch = BAWL_LAUNCH_NORMAL,
                            double delta = 0.0) {
  double pdf;
  if (ba_natural_pdf_launch(t, A, b, v, sv, k, posdrift, denom_floor,
                            BA_ACCEPT_STRICT, pdf, launch, delta))
    return log_out ? std::log(pdf) : pdf;
  return return_from_log(
    log_ba_pdf_launch(t, A, b, v, sv, k, posdrift, denom_floor, launch, delta),
    log_out);
}

inline double lba_k0_cdf_norm(double t, double A, double b, double v,
                              double sv, bool posdrift, bool log_out,
                              double denom_floor = LBA_DENOM_FLOOR) {
  return bawl_cdf_norm(t, A, b, v, sv, 0.0, posdrift, log_out, denom_floor);
}

inline double lba_k0_pdf_norm(double t, double A, double b, double v,
                              double sv, bool posdrift, bool log_out,
                              double denom_floor = LBA_DENOM_FLOOR) {
  return bawl_pdf_norm(t, A, b, v, sv, 0.0, posdrift, log_out, denom_floor);
}

// Natural-scale scalar evaluators for consumers that clamp to [0, 1] and
// tolerate tail saturation: truncation normalisers and GSL integrands.
// These accept near-0/near-1 natural values (the common case at truncation
// bounds) instead of re-deriving them through the log machinery, which
// profiling showed dominating truncated-LBA likelihoods.
inline double bawl_cdf_scalar_natural(double t, double A, double b, double v,
                                      double sv, double k, bool posdrift,
                                      double denom_floor = BAWL_DENOM_FLOOR,
                                      int launch = BAWL_LAUNCH_NORMAL,
                                      double delta = 0.0) {
  double cdf;
  if (ba_natural_cdf_launch(t, A, b, v, sv, k, posdrift, denom_floor,
                            BA_ACCEPT_CLAMP, cdf, launch, delta))
    return cdf;
  return std::exp(
    log_ba_cdf_launch(t, A, b, v, sv, k, posdrift, denom_floor, launch, delta));
}

inline double bawl_pdf_scalar_natural(double t, double A, double b, double v,
                                      double sv, double k, bool posdrift,
                                      double denom_floor = BAWL_DENOM_FLOOR,
                                      int launch = BAWL_LAUNCH_NORMAL,
                                      double delta = 0.0) {
  double pdf;
  if (ba_natural_pdf_launch(t, A, b, v, sv, k, posdrift, denom_floor,
                            BA_ACCEPT_CLAMP, pdf, launch, delta))
    return pdf;
  return std::exp(
    log_ba_pdf_launch(t, A, b, v, sv, k, posdrift, denom_floor, launch, delta));
}
double pleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift, bool log_out, int launch,
                     double delta);
double dleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift, bool log_out, int launch,
                     double delta);

// --------------------------------------------------------------------------
// Ballistic Accumulator with Leak (BAwL)
//
// Accumulator: x(t) = a*exp(-k*t) + (D/k)*(1 - exp(-k*t))
//   where a ~ Unif(0,A), D ~ N(v, sv^2) [optionally truncated D>0].
// Hits threshold b when x(t) >= b.
// k -> 0 limit recovers standard LBA exactly.
// --------------------------------------------------------------------------



// Killed-leaky BA density (hit + guess mixture):
//   f_R(t_eam) * S_K(t) * S_G(t)  +  f_G(t) * S_K(t) * S_R(t_eam)
// hit term = racer density surviving the kill and guess clocks; guess term =
// guess density surviving the kill clock and the racer (S_R = 1 - F_R).
// t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang clocks use raw t.
inline double dkilledleakyba_norm(double t, double v, double b, double A,
                                  double sv, double t0 = 0.0,
                                  double k = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
                                  bool posdrift = true, bool log_out = false,
                                  int kill_shape = 1, bool guess = false, double erlang_omega = 1.0,
                                  int launch = BAWL_LAUNCH_NORMAL, double delta = 0.0) {
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  const double t_eam = t - t0;
  const bool use_guess = guess && (lambda_g > 0.0);
  const bool use_kill = (lambda_k > 0.0);

  // Erlang survivals always use raw t
  const double log_sK = use_kill  ? erlang_log_surv(t, lambda_k, kill_shape, erlang_omega) : 0.0;
  const double log_sG = use_guess ? erlang_log_surv(t, lambda_g, kill_shape, erlang_omega) : 0.0;

  // When EAM hasn't started: f_EAM = 0, S_EAM = 1.
  if (t_eam <= 0.0) {
    if (!use_guess) return log_out ? R_NegInf : 0.0;
    // Only guess path contributes: f_G(t) * S_K(t) * S_R(0) = f_G(t) * S_K(t) * 1
    const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape, erlang_omega);
    const double log_f_guess = log_fG + log_sK;
    return log_out ? log_f_guess : std::exp(log_f_guess);
  }

  if (!use_guess && !use_kill) {
    // Natural-scale scalar consumers clamp; skip the strict wrapper's
    // near-boundary rejections (hot in truncation normalisers).
    if (!log_out) return bawl_pdf_scalar_natural(t_eam, A, b, v, sv, k, posdrift,
                                                 BAWL_DENOM_FLOOR, launch, delta);
    return dleakyba_norm(t_eam, A, b, v, sv, k, posdrift, log_out, launch, delta);
  }

  const double log_fR = dleakyba_norm(t_eam, A, b, v, sv, k, posdrift, true, launch, delta);
  const double log_f_hit = log_fR + log_sK + log_sG;
  if (!use_guess) return log_out ? log_f_hit : std::exp(log_f_hit);

  const double log_cdf_r = pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, true, launch, delta);
  const double log_sr = (log_cdf_r >= 0.0) ? R_NegInf : log1m_exp(log_cdf_r);
  const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape, erlang_omega);
  const double log_f_guess = log_fG + log_sK + log_sr;
  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double integrate_bawl_pdf_raw(double t, double v, double b, double A,
                                     double sv, double t0, double k,
                                     double lambda_g, double lambda_k,
                                     bool posdrift, int kill_shape, bool guess, double erlang_omega = 1.0,
                                     int launch = BAWL_LAUNCH_NORMAL, double delta = 0.0) {
  const double lower = guess ? 0.0 : t0;
  if (t != R_PosInf && t <= lower) return 0.0;

  std::function<double(double)> fn = [&](double u) -> double {
    return dkilledleakyba_norm(u, v, b, A, sv, t0, k, lambda_g, lambda_k,
                               posdrift, false, kill_shape, guess, erlang_omega,
                               launch, delta);
  };
  gsl_function F;
  F.function = [](double u, void* p) -> double {
    return (*static_cast<std::function<double(double)>*>(p))(u);
  };
  F.params = &fn;

  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(ws);
  double result = 0.0, err = 0.0;
  gsl_error_handler_t* old = gsl_set_error_handler_off();
  int status;
  if (t == R_PosInf) {
    status = gsl_integration_qagiu(&F, lower, 1e-8, 1e-5, 200, w, &result, &err);
  } else {
    status = gsl_integration_qags(&F, lower, t, 1e-8, 1e-5, 200, w, &result, &err);
  }
  gsl_set_error_handler(old);
  if (status != GSL_SUCCESS || !R_FINITE(result)) return 0.0;
  return std::max(0.0, std::min(1.0, result));
}

// Killed-leaky BA sub-CDF:
// P(T_R <= t, T_R < T_K) with kill_shape=1 (exponential) or kill_shape=2 (Erlang-2).
// With guess=true: mixture CDF = 1 - S_R(t_eam)*S_K(t)*S_G(t).
// t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
inline double pkilledleakyba_norm(double t, double v, double b, double A,
                                  double sv, double t0 = 0.0,
                                  double k = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
                                  bool posdrift = true, bool log_out = false,
                                  int kill_shape = 1, bool guess = false, double erlang_omega = 1.0,
                                  int launch = BAWL_LAUNCH_NORMAL, double delta = 0.0) {
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  const double t_eam = t - t0;
  const bool use_guess = guess && (lambda_g > 0.0);
  const bool use_kill = (lambda_k > 0.0);

  // When EAM hasn't started: P(EAM fires by t_eam) = 0, S_R = 1.
  if (t_eam <= 0.0) {
    if (!use_guess) return log_out ? R_NegInf : 0.0;
    if (!use_kill) {
      // Guess-only response process before EAM onset.
      const double log_sG = erlang_log_surv(t, lambda_g, kill_shape, erlang_omega);
      const double out = -std::expm1(log_sG);
      return log_out ? safe_log(out) : out;
    }
    // Before EAM onset, the only observed response is a guess before the kill clock:
    // integral_0^t f_G(u) S_K(u) du.  Kill wins are omissions, not CDF mass.
    std::function<double(double)> fn = [&](double u) -> double {
      return std::exp(erlang_log_pdf(u, lambda_g, kill_shape, erlang_omega) +
                      erlang_log_surv(u, lambda_k, kill_shape, erlang_omega));
    };
    gsl_function F;
    F.function = [](double u, void* p) -> double {
      return (*static_cast<std::function<double(double)>*>(p))(u);
    };
    F.params = &fn;
    static thread_local GslWorkspacePtr ws_pre_eam(nullptr, &gsl_integration_workspace_free);
    gsl_integration_workspace* w = ensure_gsl_workspace(ws_pre_eam);
    double out = 0.0, err = 0.0;
    gsl_error_handler_t* old = gsl_set_error_handler_off();
    gsl_integration_qags(&F, 0.0, t, 1e-8, 1e-5, 200, w, &out, &err);
    gsl_set_error_handler(old);
    out = std::max(0.0, std::min(1.0, out));
    return log_out ? safe_log(out) : out;
  }

  if (!use_guess && !use_kill) {
    if (!log_out) return bawl_cdf_scalar_natural(t_eam, A, b, v, sv, k, posdrift,
                                                 BAWL_DENOM_FLOOR, launch, delta);
    return pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, log_out, launch, delta);
  }

  if (use_guess && use_kill) {
    const double out = integrate_bawl_pdf_raw(
      t, v, b, A, sv, t0, k, lambda_g, lambda_k, posdrift, kill_shape, true,
      erlang_omega, launch, delta
    );
    return log_out ? safe_log(out) : out;
  }

  const double lambda = use_guess ? lambda_g : lambda_k;
  if (use_guess) {
    // 1 - S_R(t_eam) * S_G(t); erlang uses raw t
    const double log_cdf_r = pleakyba_norm(t_eam, A, b, v, sv, k, posdrift,
                                           true, launch, delta);
    const double log_sr = (log_cdf_r >= 0.0) ? R_NegInf : log1m_exp(log_cdf_r);
    const double log_sg = erlang_log_surv(t, lambda, kill_shape, erlang_omega);
    const double log_val = log1m_exp(log_sr + log_sg);
    return log_out ? log_val : std::exp(log_val);
  }
  if (sv <= 0.0 || b < A || b <= 0.0) return log_out ? R_NegInf : 0.0;

  // Pure kill path: integrate f_EAM(u - t0) * S_K(u) over raw time.
  // No closed normal-form primitive remains once leak and clock survival are combined.
  const double out = integrate_bawl_pdf_raw(
    t, v, b, A, sv, t0, k, 0.0, lambda, posdrift, kill_shape, false, erlang_omega,
    launch, delta
  );
  return log_out ? safe_log(out) : out;
}








#endif

#ifndef lba_h
#define lba_h

// This header may be included by exactly ONE translation unit (particle_ll.cpp,
// directly and via utils.h) because it defines [[Rcpp::export]] functions
// (pleakyba_norm, dleakyba_norm, …) that RcppExports links to. Non-exported
// free helpers here are
// marked `inline`; the exported ones must NOT be inline.
#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "wald_functions.h"  // pnorm_std() — fast normal CDF under USE_FAST_PNORM
#include "composite_functions.h"  // clamp_pos, safe_log
#include "gaussian.h"
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;

// Route through pnorm_std so USE_FAST_PNORM applies to LBA as well as RDM/Wald.
// pnorm(q, mean, sd) = pnorm_std((q - mean) / sd) for sd > 0.
inline double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log = false){
  return pnorm_std((q - mean) / sd, lower, log);
}

inline double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

// pnorm_std's fast log-tail approximation deliberately returns -Inf beyond
// its useful range. LBA needs those tails to remain representable, so fall
// back to R's direct log-tail implementation rather than materialising a
// probability and taking its log.
inline double pnorm_log_direct(double x, bool lower = true) {
  double out = pnorm_std(x, lower, true);
  if (out == R_NegInf && emc2_isfinite(x)) {
    out = R::pnorm(x, 0.0, 1.0, lower, true);
  }
  return out;
}

// Return log(Phi(hi) - Phi(lo)) without materialising either probability.
// Using the upper tails when both arguments are positive avoids subtracting
// two values that are both numerically equal to one.
inline double log_normal_interval(double lo, double hi) {
  if (!(hi > lo)) return R_NegInf;
  if (lo >= 0.0) {
    return log_diff_exp(pnorm_log_direct(lo, false),
                        pnorm_log_direct(hi, false));
  }
  return log_diff_exp(pnorm_log_direct(hi, true),
                      pnorm_log_direct(lo, true));
}

// Let Q(z) = 1 - Phi(z). The positive quantity
//   M(z) = phi(z) - z Q(z)
// is decreasing and satisfies
//   integral_lo^hi Q(z) dz = M(lo) - M(hi).
// This is the positive-integrand form of the LBA CDF and avoids the severe
// cancellation in the usual "1 + correction" expression at early times.
inline double log_normal_q_antiderivative_abs(double z) {
  const double log_phi = dnormP(z, 0.0, 1.0, true);
  const double log_q = pnorm_log_direct(z, false);

  if (z > 0.0) {
    // M(z) = phi(z) * [1 - z Q(z) / phi(z)].
    const double log_ratio = std::log(z) + log_q - log_phi;
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

// A signed log-scale product.  The coefficient is kept on its natural scale
// as an input, but the product itself is never formed in natural space.
inline signed_log signed_log_product(double coefficient, double log_factor) {
  if (coefficient == 0.0 || log_factor == R_NegInf) {
    return {R_NegInf, 0};
  }
  return make_signed_log(std::log(std::fabs(coefficient)) + log_factor,
                         coefficient > 0.0 ? 1 : -1);
}

constexpr double BAWL_K_EPS = 1e-10;
constexpr double BAWL_NATURAL_Z_MAX = 7.5;
constexpr double BAWL_NATURAL_MIN_SPAN = 1e-6;
constexpr double BAWL_NATURAL_REL_TOL = 1e-10;
constexpr double LBA_DENOM_FLOOR = 1e-10;
constexpr double BAWL_DENOM_FLOOR = 1e-300;

inline double log_positive_normalizer(double v, double sv, bool posdrift,
                                     double denom_floor = LBA_DENOM_FLOOR) {
  if (!posdrift) return 0.0;

  // Match the selected model's pmax(pnorm(v / sv), floor) normalization,
  // but apply the floor directly to its logarithm.
  const double log_denom = pnorm_log_direct(v / sv, true);
  const double log_floor = std::log(denom_floor);
  return (log_denom < log_floor) ? log_floor : log_denom;
}

inline double return_from_log(double log_value, bool log_out) {
  return log_out ? log_value : std::exp(log_value);
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

inline bool natural_probability(double p) {
  return p > 0.0 && p < 1.0;
}

// The CDF is often consumed as log(1 - F), so values too close to either
// endpoint must stay on the stable log path even when the probability itself
// is representable.  In particular, 1 - mean{Phi(z)} loses the low tail.
inline bool natural_cdf_safe(double p) {
  constexpr double cdf_margin = 1e-8;
  return p > cdf_margin && p < 1.0 - cdf_margin;
}

// Standard LBA is the exact k=0 limit of BAwL.  These guarded natural-space
// helpers are used by both public model entry points; false means that the
// caller must use the stable log-space calculation.
inline bool lba_k0_natural_cdf(double t, double A, double b, double v,
                               double sv, bool posdrift, double &cdf,
                               double denom_floor = LBA_DENOM_FLOOR) {
  if (!(sv > 0.0)) return false;
  if (t == R_PosInf) {
    cdf = posdrift ? 1.0 : pnorm_std(v / sv, true, false);
    return R_FINITE(cdf);
  }
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) return false;

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;
  const double zs = t * sv;
  const double cmz = b - t * v;
  if (!(zs > 0.0)) return false;

  if (A <= 1e-10) {
    const double z = cmz / zs;
    if (std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(-z, true, false) / denom;
    return natural_cdf_safe(cdf);
  }

  const double z_hi = cmz / zs;
  const double z_lo = (cmz - A) / zs;
  if (!natural_normal_interval_safe(z_lo, z_hi)) return false;
  const double span = z_hi - z_lo;
  if (!(span > 0.0)) return false;
  if (span < BAWL_NATURAL_MIN_SPAN) {
    cdf = pnorm_std(-(z_lo + 0.5 * span), true, false) / denom;
    return natural_cdf_safe(cdf);
  }

  // H(z) = z Phi(z) + phi(z), so H(hi)-H(lo) is the integral of Phi.
  const double phi_hi = dnormP(z_hi, 0.0, 1.0, false);
  const double phi_lo = dnormP(z_lo, 0.0, 1.0, false);
  const double H_hi = z_hi * pnorm_std(z_hi, true, false) + phi_hi;
  const double H_lo = z_lo * pnorm_std(z_lo, true, false) + phi_lo;
  const double integral = H_hi - H_lo;
  const double scale = std::fabs(H_hi) + std::fabs(H_lo);
  if (!(integral > 0.0) ||
      integral <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;

  const double avg_phi = integral / span;
  if (!natural_cdf_safe(1.0 - avg_phi)) return false;
  cdf = (1.0 - avg_phi) / denom;
  return natural_cdf_safe(cdf);
}

inline bool lba_k0_natural_pdf(double t, double A, double b, double v,
                               double sv, bool posdrift, double &pdf,
                               double denom_floor = LBA_DENOM_FLOOR) {
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) return false;

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;
  const double zs = t * sv;
  const double cmz = b - t * v;
  if (!(zs > 0.0)) return false;

  if (A <= 1e-10) {
    const double z = cmz / zs;
    if (std::fabs(z) > BAWL_NATURAL_Z_MAX) return false;
    const double scale = b / (t * t * sv * denom);
    pdf = scale * dnormP(z, 0.0, 1.0, false);
    return pdf > 0.0 && pdf < R_PosInf;
  }

  const double z_hi = cmz / zs;
  const double z_lo = (cmz - A) / zs;
  if (!natural_normal_interval_safe(z_lo, z_hi)) return false;
  const double span = z_hi - z_lo;
  if (!(span > 0.0)) return false;
  if (span < BAWL_NATURAL_MIN_SPAN) {
    const double z_mid = z_lo + 0.5 * span;
    const double scale = (b - 0.5 * A) / (t * t * sv * denom);
    pdf = scale * dnormP(z_mid, 0.0, 1.0, false);
    return pdf > 0.0 && pdf < R_PosInf;
  }

  const double dphi = pnorm_std(z_hi, true, false) -
    pnorm_std(z_lo, true, false);
  const double dnorm = dnormP(z_lo, 0.0, 1.0, false) -
    dnormP(z_hi, 0.0, 1.0, false);
  const double term1 = v * dphi;
  const double term2 = sv * dnorm;
  const double bracket = term1 + term2;
  const double scale = std::fabs(term1) + std::fabs(term2);
  if (!(bracket > 0.0) ||
      bracket <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;

  pdf = bracket / (A * denom);
  return pdf > 0.0 && pdf < R_PosInf;
}

inline double log_lba_k0_cdf_norm(double t, double A, double b, double v,
                                  double sv, bool posdrift = true,
                                  double denom_floor = LBA_DENOM_FLOOR) {
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0) || !(sv > 0.0))
    return R_NegInf;
  if (t == R_PosInf)
    return posdrift ? 0.0 : pnorm_log_direct(v / sv, true);
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  if (A <= 1e-10)
    return pnorm_log_direct((v - b / t) / sv, true) - log_denom;

  const double zs = t * sv;
  if (!(zs > 0.0)) return R_NegInf;
  const double z_hi = (b - t * v) / zs;
  const double z_lo = (b - A - t * v) / zs;
  return std::log(zs) - std::log(A) +
    log_normal_q_interval(z_lo, z_hi) - log_denom;
}

inline double log_lba_k0_pdf_norm(double t, double A, double b, double v,
                                  double sv, bool posdrift = true,
                                  double denom_floor = LBA_DENOM_FLOOR) {
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0) || !(sv > 0.0))
    return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  const double zs = t * sv;
  if (!(zs > 0.0)) return R_NegInf;

  if (A <= 1e-10) {
    const double z = (b - t * v) / zs;
    return std::log(b) - 2.0 * std::log(t) - std::log(sv) +
      dnormP(z, 0.0, 1.0, true) - log_denom;
  }

  const double z_hi = (b - t * v) / zs;
  const double z_lo = (b - A - t * v) / zs;
  const double log_dphi = log_normal_interval(z_lo, z_hi);
  signed_log bracket = signed_log_product(v, log_dphi);
  signed_log phi_diff = signed_log_sub(
    make_signed_log(dnormP(z_lo, 0.0, 1.0, true), 1),
    make_signed_log(dnormP(z_hi, 0.0, 1.0, true), 1));
  // The density numerator is v * ΔPhi + sv * Δphi.
  if (phi_diff.sign != 0) phi_diff.log_abs += std::log(sv);
  bracket = signed_log_add(bracket, phi_diff);
  if (bracket.sign <= 0 || bracket.log_abs == R_NegInf) return R_NegInf;
  return bracket.log_abs - std::log(A) - log_denom;
}

inline double lba_k0_cdf_norm(double t, double A, double b, double v,
                              double sv, bool posdrift, bool log_out,
                              double denom_floor = LBA_DENOM_FLOOR) {
  double cdf;
  if (lba_k0_natural_cdf(t, A, b, v, sv, posdrift, cdf, denom_floor))
    return log_out ? std::log(cdf) : cdf;
  return return_from_log(log_lba_k0_cdf_norm(t, A, b, v, sv, posdrift, denom_floor), log_out);
}

inline double lba_k0_pdf_norm(double t, double A, double b, double v,
                              double sv, bool posdrift, bool log_out,
                              double denom_floor = LBA_DENOM_FLOOR) {
  double pdf;
  if (lba_k0_natural_pdf(t, A, b, v, sv, posdrift, pdf, denom_floor))
    return log_out ? std::log(pdf) : pdf;
  return return_from_log(log_lba_k0_pdf_norm(t, A, b, v, sv, posdrift, denom_floor), log_out);
}

// Raw likelihoods may floor finite natural underflow without reconstructing a
// log tail that will be discarded by min_ll.  They still route genuine
// overflow, PDF cancellation, and near-one CDFs through the stable helpers.
inline bool lba_k0_raw_natural_cdf(double t, double A, double b, double v,
                                   double sv, bool posdrift, double &cdf) {
  if (!(sv > 0.0)) {
    cdf = 0.0;
    return true;
  }
  if (t == R_PosInf) {
    if (posdrift) {
      cdf = 1.0;
      return true;
    }
    cdf = pnorm_std(v / sv, true, false);
    return R_FINITE(cdf) && cdf < 1.0 - 1e-8;
  }
  if (!(t > 0.0) || !(sv > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) {
    cdf = 0.0;
    return true;
  }
  double denom = 1.0;
  if (posdrift) {
    denom = pnorm_std(v / sv, true, false);
    if (denom < 1e-10) denom = 1e-10;
  }
  if (A <= 1e-10) {
    cdf = pnorm_std((v - b / t) / sv, true, false) / denom;
  } else {
    const double zs = t * sv;
    const double cmz = b - t * v;
    const double z_hi = cmz / zs;
    const double z_lo = (cmz - A) / zs;
    cdf = (1.0 + (zs * (dnormP(z_lo) - dnormP(z_hi)) +
      (cmz - A) * pnorm_std(z_lo) - cmz * pnorm_std(z_hi)) / A) / denom;
  }
  // Near-one CDFs are used as log survivors and need the stable log path.
  return R_FINITE(cdf) && cdf < 1.0 - 1e-8;
}

inline bool lba_k0_raw_natural_pdf(double t, double A, double b, double v,
                                   double sv, bool posdrift, double &pdf) {
  if (!(t > 0.0) || !(sv > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) {
    pdf = 0.0;
    return true;
  }
  double denom = 1.0;
  if (posdrift) {
    denom = pnorm_std(v / sv, true, false);
    if (denom < 1e-10) denom = 1e-10;
  }
  if (A <= 1e-10) {
    pdf = dnormP(b / t, v, sv, false) * b / (t * t * denom);
  } else {
    const double zs = t * sv;
    const double cmz = b - t * v;
    const double z_hi = cmz / zs;
    const double z_lo = (cmz - A) / zs;
    const double dphi = pnorm_std(z_hi) - pnorm_std(z_lo);
    const double dnorm = dnormP(z_lo) - dnormP(z_hi);
    const double term1 = v * dphi;
    const double term2 = sv * dnorm;
    const double bracket = term1 + term2;
    const double scale = std::fabs(term1) + std::fabs(term2);
    if (bracket > 0.0 &&
        bracket <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) {
      return false;
    }
    pdf = bracket / (A * denom);
  }
  // A finite zero/negative result is already the raw likelihood floor; an
  // infinite result is the case that needs the stable log evaluator.
  return R_FINITE(pdf);
}

// --------------------------------------------------------------------------
// Shared log-space BAwL core
//
// For t > 0, define
//   c = (v - k*b/(1-exp(-k*t))) / sv,
//   m = k*exp(-k*t) / ((1-exp(-k*t))*sv).
// Then the hit probability conditional on start a is Phi(c + m*a).
// k=0 is evaluated by its exact limit, c=(v-b/t)/sv and m=1/(t*sv).
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
    c = 0.0;
    m = (1.0 / t) / sv;
    if (log_jacobian != nullptr) *log_jacobian = -2.0 * std::log(t);
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

// Guarded natural-space BAwL core for k>0.  The exact k=0 limit is delegated
// to the LBA helpers above so the shared evaluator does not duplicate model
// semantics while avoiding generic BAwL work in the nested case.
inline bool bawl_natural_cdf(double t, double A, double b, double v,
                             double sv, double k, bool posdrift,
                             double &cdf) {
  if (!(k > BAWL_K_EPS))
    return lba_k0_natural_cdf(t, A, b, v, sv, posdrift, cdf,
                              BAWL_DENOM_FLOOR);
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) return false;

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, BAWL_DENOM_FLOOR)) return false;
  double c, m;
  bawl_threshold_terms(t, b, sv, k, c, m);
  if (!(m > 0.0)) return false;
  c += v / sv;

  if (A <= 1e-10) {
    if (std::fabs(c) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(c, true, false) / denom;
    return natural_cdf_safe(cdf);
  }

  const double span = m * A;
  const double z_hi = c + span;
  if (!(span > 0.0) ||
      !natural_normal_interval_safe(c, z_hi)) return false;
  if (span < BAWL_NATURAL_MIN_SPAN) {
    cdf = pnorm_std(c + 0.5 * span, true, false) / denom;
    return natural_cdf_safe(cdf);
  }

  const double H_hi = z_hi * pnorm_std(z_hi, true, false) +
    dnormP(z_hi, 0.0, 1.0, false);
  const double H_lo = c * pnorm_std(c, true, false) +
    dnormP(c, 0.0, 1.0, false);
  const double integral = H_hi - H_lo;
  const double scale = std::fabs(H_hi) + std::fabs(H_lo);
  if (!(integral > 0.0) ||
      integral <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;

  cdf = integral / (A * m * denom);
  return natural_cdf_safe(cdf);
}

inline bool bawl_natural_pdf(double t, double A, double b, double v,
                             double sv, double k, bool posdrift,
                             double &pdf) {
  if (!(k > BAWL_K_EPS))
    return lba_k0_natural_pdf(t, A, b, v, sv, posdrift, pdf,
                              BAWL_DENOM_FLOOR);
  if (!(t > 0.0) || !(A >= 0.0) || !(b >= A) || !(b > 0.0)) return false;

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, BAWL_DENOM_FLOOR)) return false;
  double c, m, jacobian;
  bawl_threshold_terms(t, b, sv, k, c, m, nullptr, &jacobian);
  if (!(m > 0.0)) return false;
  c += v / sv;

  if (!(jacobian > 0.0) || !(jacobian < R_PosInf)) return false;
  if (A <= 1e-10) {
    if (std::fabs(c) > BAWL_NATURAL_Z_MAX) return false;
    pdf = jacobian * b * dnormP(c, 0.0, 1.0, false) /
      (sv * denom);
    return pdf > 0.0 && pdf < R_PosInf;
  }

  const double span = m * A;
  const double z_hi = c + span;
  if (!(span > 0.0) ||
      !natural_normal_interval_safe(c, z_hi)) return false;
  if (span < BAWL_NATURAL_MIN_SPAN) {
    const double z_mid = c + 0.5 * span;
    pdf = jacobian * (b - 0.5 * A) * dnormP(z_mid, 0.0, 1.0, false) /
      (sv * denom);
    return pdf > 0.0 && pdf < R_PosInf;
  }

  const double dphi = pnorm_std(z_hi, true, false) -
    pnorm_std(c, true, false);
  const double dnorm = dnormP(z_hi, 0.0, 1.0, false) -
    dnormP(c, 0.0, 1.0, false);
  const double term1 = (b * m + c) * dphi;
  const double term2 = dnorm;
  const double bracket = term1 + term2;
  const double scale = std::fabs(term1) + std::fabs(term2);
  if (!(bracket > 0.0) ||
      bracket <= BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) return false;

  const double denom_scale = A * sv * m * m * denom;
  if (!(denom_scale > 0.0) || !(denom_scale < R_PosInf)) return false;
  pdf = jacobian * bracket / denom_scale;
  return pdf > 0.0 && pdf < R_PosInf;
}

inline double log_bawl_cdf_norm(double t, double A, double b, double v,
                                double sv, double k, bool posdrift = true) {
  if (k <= BAWL_K_EPS)
    return log_lba_k0_cdf_norm(t, A, b, v, sv, posdrift,
                               BAWL_DENOM_FLOOR);
  if (t <= 0.0 || !(sv > 0.0) || !(b >= A) || !(b > 0.0)) return R_NegInf;

  const double log_denom = log_positive_normalizer(v, sv, posdrift,
                                                    BAWL_DENOM_FLOOR);
  double c, m;
  bawl_threshold_terms(t, b, sv, k, c, m);

  c = v / sv + c;

  if (A <= 1e-10 || !(m > 0.0)) {
    return pnorm_log_direct(c, true) - log_denom;
  }

  const double span = m * A;
  if (!(span > 1e-8) || !emc2_isfinite(span)) {
    // Midpoint evaluation is the continuous limit as the start range
    // collapses, and avoids manufacturing a difference of equal antiderives.
    const double zmid = c + 0.5 * span;
    return pnorm_log_direct(zmid, true) - log_denom;
  }

  const double log_integral = log_normal_phi_integral(c, c + span);
  const double out = log_integral - std::log(A) - std::log(m) - log_denom;
  return out >= 0.0 ? 0.0 : out;
}

inline double log_bawl_pdf_norm(double t, double A, double b, double v,
                                double sv, double k, bool posdrift = true) {
  if (k <= BAWL_K_EPS)
    return log_lba_k0_pdf_norm(t, A, b, v, sv, posdrift,
                               BAWL_DENOM_FLOOR);
  if (t <= 0.0 || !(sv > 0.0) || !(b >= A) || !(b > 0.0)) return R_NegInf;

  const double log_denom = log_positive_normalizer(v, sv, posdrift,
                                                    BAWL_DENOM_FLOOR);
  double c, m, log_jacobian;
  bawl_threshold_terms(t, b, sv, k, c, m, &log_jacobian);
  c = v / sv + c;

  if (!(log_jacobian > R_NegInf)) return R_NegInf;
  if (A <= 1e-10 || !(m > 0.0)) {
    return log_jacobian + std::log(b) + dnormP(c, 0.0, 1.0, true) -
      std::log(sv) - log_denom;
  }

  const double span = m * A;
  if (!(span > 1e-8) || !emc2_isfinite(span)) {
    const double zmid = c + 0.5 * span;
    return log_jacobian + std::log(b - 0.5 * A) +
      dnormP(zmid, 0.0, 1.0, true) - std::log(sv) - log_denom;
  }

  const double hi = c + span;
  const double log_dphi = log_normal_interval(c, hi);

  // m^2 * integral_0^A (b-a) phi(c+m*a) da
  //   = (b*m+c) [Phi(hi)-Phi(c)] + phi(hi)-phi(c).
  signed_log bracket = signed_log_product(b * m + c, log_dphi);
  const signed_log log_phi_diff = signed_log_sub(
    make_signed_log(dnormP(hi, 0.0, 1.0, true), 1),
    make_signed_log(dnormP(c, 0.0, 1.0, true), 1));
  bracket = signed_log_add(bracket, log_phi_diff);
  if (bracket.sign <= 0 || bracket.log_abs == R_NegInf) return R_NegInf;

  return log_jacobian - std::log(A) - std::log(sv) - log_denom -
    2.0 * std::log(m) + bracket.log_abs;
}

inline double bawl_cdf_norm(double t, double A, double b, double v,
                            double sv, double k, bool posdrift,
                            bool log_out) {
  if (k <= BAWL_K_EPS)
    return lba_k0_cdf_norm(t, A, b, v, sv, posdrift, log_out,
                           BAWL_DENOM_FLOOR);
  double cdf;
  if (bawl_natural_cdf(t, A, b, v, sv, k, posdrift, cdf))
    return log_out ? std::log(cdf) : cdf;
  return return_from_log(log_bawl_cdf_norm(t, A, b, v, sv, k, posdrift), log_out);
}

inline double bawl_pdf_norm(double t, double A, double b, double v,
                            double sv, double k, bool posdrift,
                            bool log_out) {
  if (k <= BAWL_K_EPS)
    return lba_k0_pdf_norm(t, A, b, v, sv, posdrift, log_out,
                           BAWL_DENOM_FLOOR);
  double pdf;
  if (bawl_natural_pdf(t, A, b, v, sv, k, posdrift, pdf))
    return log_out ? std::log(pdf) : pdf;
  return return_from_log(log_bawl_pdf_norm(t, A, b, v, sv, k, posdrift), log_out);
}

// --------------------------------------------------------------------------
// Ballistic Accumulator with Leak (BAwL)
//
// Accumulator: x(t) = a*exp(-k*t) + (D/k)*(1 - exp(-k*t))
//   where a ~ Unif(0,A), D ~ N(v, sv^2) [optionally truncated D>0].
// Hits threshold b when x(t) >= b.
// k -> 0 limit recovers standard LBA exactly.
// --------------------------------------------------------------------------

// CDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double pleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false) {
  // At infinite time, k=0 has the usual LBA limit.  For k>0 the accumulator
  // approaches D/k, so only drifts above k*b can finish; retain that
  // defective upper tail instead of returning one.
  if (t == R_PosInf && k <= 1e-10) {
    const double log_cdf = posdrift ? 0.0 : pnorm_log_direct(v / sv, true);
    return return_from_log(log_cdf, log_out);
  }
  return bawl_cdf_norm(t, A, b, v, sv, k, posdrift, log_out);
}

// PDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double dleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false) {
  return bawl_pdf_norm(t, A, b, v, sv, k, posdrift, log_out);
}

// Killed-leaky BA density (hit + guess mixture):
//   f_R(t_eam) * S_K(t) * S_G(t)  +  f_G(t) * S_K(t) * S_R(t_eam)
// hit term = racer density surviving the kill and guess clocks; guess term =
// guess density surviving the kill clock and the racer (S_R = 1 - F_R).
// t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang clocks use raw t.
inline double dkilledleakyba_norm(double t, double v, double b, double A,
                                  double sv, double t0 = 0.0,
                                  double k = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
                                  bool posdrift = true, bool log_out = false,
                                  int kill_shape = 1, bool guess = false, double erlang_omega = 1.0) {
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

  if (!use_guess && !use_kill) return dleakyba_norm(t_eam, A, b, v, sv, k, posdrift, log_out);

  const double log_fR = dleakyba_norm(t_eam, A, b, v, sv, k, posdrift, true);
  const double log_f_hit = log_fR + log_sK + log_sG;
  if (!use_guess) return log_out ? log_f_hit : std::exp(log_f_hit);

  const double log_cdf_r = pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, true);
  const double log_sr = (log_cdf_r >= 0.0) ? R_NegInf : log1m_exp(log_cdf_r);
  const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape, erlang_omega);
  const double log_f_guess = log_fG + log_sK + log_sr;
  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double integrate_bawl_pdf_raw(double t, double v, double b, double A,
                                     double sv, double t0, double k,
                                     double lambda_g, double lambda_k,
                                     bool posdrift, int kill_shape, bool guess, double erlang_omega = 1.0) {
  const double lower = guess ? 0.0 : t0;
  if (t != R_PosInf && t <= lower) return 0.0;

  std::function<double(double)> fn = [&](double u) -> double {
    return dkilledleakyba_norm(u, v, b, A, sv, t0, k, lambda_g, lambda_k,
                               posdrift, false, kill_shape, guess, erlang_omega);
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
                                  int kill_shape = 1, bool guess = false, double erlang_omega = 1.0) {
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

  if (!use_guess && !use_kill) return pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, log_out);

  if (use_guess && use_kill) {
    const double out = integrate_bawl_pdf_raw(
      t, v, b, A, sv, t0, k, lambda_g, lambda_k, posdrift, kill_shape, true, erlang_omega
    );
    return log_out ? safe_log(out) : out;
  }

  const double lambda = use_guess ? lambda_g : lambda_k;
  if (use_guess) {
    // 1 - S_R(t_eam) * S_G(t); erlang uses raw t
    const double log_cdf_r = pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, true);
    const double log_sr = (log_cdf_r >= 0.0) ? R_NegInf : log1m_exp(log_cdf_r);
    const double log_sg = erlang_log_surv(t, lambda, kill_shape, erlang_omega);
    const double log_val = log1m_exp(log_sr + log_sg);
    return log_out ? log_val : std::exp(log_val);
  }
  if (sv <= 0.0 || b < A || b <= 0.0) return log_out ? R_NegInf : 0.0;

  // Pure kill path: integrate f_EAM(u - t0) * S_K(u) over raw time.
  // No closed normal-form primitive remains once leak and clock survival are combined.
  const double out = integrate_bawl_pdf_raw(
    t, v, b, A, sv, t0, k, 0.0, lambda, posdrift, kill_shape, false, erlang_omega
  );
  return log_out ? safe_log(out) : out;
}

// [[Rcpp::export]]
NumericVector dkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    pdf[i] = dkilledleakyba_norm(t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
                                 pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
                                 posdrift, log_out, kill_shape, guess, omega);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    cdf[i] = pkilledleakyba_norm(t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
                                 pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
                                 posdrift, log_out, kill_shape, guess, omega);
  }
  return cdf;
}

// Vectorised R-callable wrappers (recycle scalar parameters).
// [[Rcpp::export]]
NumericVector dleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = dleakyba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i), posdrift);
  return pdf;
}

// [[Rcpp::export]]
NumericVector pleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = pleakyba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i), posdrift);
  return cdf;
}

#endif

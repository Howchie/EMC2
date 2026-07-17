#ifndef EMC2_BAWL_GEOMETRY_H
#define EMC2_BAWL_GEOMETRY_H

// Authoritative fixed-time LBA/leak geometry for the correlated BAwL routes
// (bawl_corr_exact_kernel_plan.md, design rule A/B).
//
// BAwLTimeGeometry owns every quantity that depends on (t, t0, A, b, k) but
// not on the drift mean or SD: leak factors, drift bounds L/U, piecewise
// survivor coefficients alpha/beta, and cause-density coefficients
// gamma0/gamma1.  It is built on the same helpers as the scalar kernels
// (bawl_leak_factors(), the exact k -> 0 limit in bawl_threshold_terms());
// there is no second implementation of exp(-k*tau) or the small-k*tau branch.
//
// BAwLPreparedRow composes a time geometry with a row's conditional drift
// scale for the shared-factor GH fallback: only the conditional mean varies
// over nodes, so the standardised threshold offset, span, and Jacobian are
// prepared once per row and particle.
//
// Like model_LBA.h (which this header requires for the shared numerical
// helpers and guard constants), it may be included by exactly one
// translation unit: particle_ll.cpp.

#include <cmath>
#include <cstdint>
#include "model_LBA.h"

enum class BAwLTimeStatus : uint8_t {
  valid,        // 0 < tau < Inf, interval start (A > 0)
  point_start,  // 0 < tau < Inf, point start (A <= BAWL_A_EPS)
  not_started,  // tau <= 0: survivor one, cause density zero
  infinite,     // t = +Inf: limiting indicator survivor only
  invalid       // malformed b/A (b <= 0 or b < A) or non-finite inputs
};

struct BAwLTimeGeometry {
  double tau = 0.0;
  double b = 0.0;
  double A = 0.0;
  double k = 0.0;
  double E = 1.0;       // exp(-k*tau); 1 in the exact k -> 0 limit
  double G = 0.0;       // 1 - E
  double C1 = 0.0;      // k*b/G   (b/tau at k -> 0);  drift upper bound U
  double C2 = 0.0;      // k*E/G   (1/tau at k -> 0)
  double L = 0.0;       // survivor 1 for V <= L
  double U = 0.0;       // survivor 0 for V >= U; point-start pinned drift v*
  double alpha = 0.0;   // survivor = alpha + beta*V on (L, U)
  double beta = 0.0;
  double gamma0 = 0.0;  // cause weight = gamma0 + gamma1*V on (L, U)
  double gamma1 = 0.0;
  double dv_dt = 0.0;   // point start: |dv*/dt| = C1*C2
  BAwLTimeStatus status = BAwLTimeStatus::invalid;
};

// Build the fixed-time geometry for one racer.  `b` is the absolute
// threshold B + A.  For t = +Inf the survivor degenerates to the indicator
// 1(V <= k*b) (1(V <= 0) for the k = 0 limit), represented as L = U with
// status `infinite`.
inline BAwLTimeGeometry bawl_time_geometry(double t, double t0, double A,
                                           double b, double k) {
  BAwLTimeGeometry g;
  g.A = A;
  g.b = b;
  g.k = k;
  if (!(A >= 0.0) || !(b > 0.0) || !(b >= A) || ISNAN(t) || ISNAN(t0) ||
      ISNAN(k)) {
    return g;  // invalid
  }
  const double tau = t - t0;
  g.tau = tau;
  if (ISNAN(tau) || !(tau > 0.0)) {
    g.status = BAwLTimeStatus::not_started;
    return g;
  }
  if (tau == R_PosInf) {
    const double v_limit = (k > BAWL_K_EPS) ? k * b : 0.0;
    g.E = 0.0;
    g.G = 1.0;
    g.C1 = v_limit;
    g.C2 = 0.0;
    g.L = v_limit;
    g.U = v_limit;
    g.status = BAwLTimeStatus::infinite;
    return g;
  }

  if (k <= BAWL_K_EPS) {
    // Exact k -> 0 limit (standard LBA), matching bawl_threshold_terms().
    g.E = 1.0;
    g.G = 0.0;
    g.C1 = b / tau;
    g.C2 = 1.0 / tau;
  } else {
    double E, G;
    bawl_leak_factors(k * tau, E, G);
    const double G_safe = clamp_pos(G, 1e-300);
    g.E = E;
    g.G = G;
    g.C1 = (k * b) / G_safe;
    g.C2 = (k * E) / G_safe;
  }
  g.U = g.C1;
  g.L = g.C1 - g.C2 * A;

  if (A <= BAWL_A_EPS) {
    // Point start: drift pinned at v* = U with time Jacobian |dv*/dt|.
    g.L = g.U;
    g.dv_dt = g.C1 * g.C2;
    g.status = BAwLTimeStatus::point_start;
    return g;
  }

  const double inv_C2A = 1.0 / (g.C2 * A);
  g.alpha = g.C1 * inv_C2A;
  g.beta = -inv_C2A;
  // Cause weight h(V) = gamma0 + gamma1*V on (L, U); the k = 0 member is
  // gamma0 = 0, gamma1 = 1/A exactly (E = 1).
  g.gamma1 = 1.0 / (g.E * A);
  g.gamma0 = -(g.k > BAWL_K_EPS ? g.k * b : 0.0) / (g.E * A);
  g.status = BAwLTimeStatus::valid;
  return g;
}

// --------------------------------------------------------------------------
// Prepared per-row conditional view for the shared-factor GH fallback.
//
// Conditional on a latent node z the row's drift is N(v0 + slope*z, sv^2)
// with sv the residual (conditional) SD.  Everything except the drift mean
// is node-invariant.  The endpoint evaluators below mirror the acceptance
// semantics of the raw BAwL kernels (BA_ACCEPT_RAW natural attempt, then the
// authoritative log-space branch) but return unnormalised quantities: the
// unrestricted log PDF f0 and log CDF F0, the positivity constant
// q = Phi(v_q / sv), and the two survivor modes 1 - F0 and q - F0.
// --------------------------------------------------------------------------

struct BAwLPreparedRow {
  // Factor model (node-invariant, per particle):
  double v0 = 0.0;      // marginal drift mean
  double slope = 0.0;   // v_q = v0 + slope * z
  double sv = 1.0;      // residual (conditional) drift SD
  double inv_sv = 1.0;
  // Time-geometry composition:
  double c0 = 0.0;      // -C1 / sv; node threshold offset c = c0 + v_q/sv
  double m = 0.0;       // C2 / sv
  double span = 0.0;    // m * A
  double jacobian = 0.0;
  double log_jacobian = R_NegInf;
  double b = 0.0;
  double A = 0.0;
  BAwLTimeStatus status = BAwLTimeStatus::invalid;
};

inline BAwLPreparedRow bawl_prepare_row(const BAwLTimeGeometry& g, double v0,
                                        double slope, double sv_residual) {
  BAwLPreparedRow r;
  r.v0 = v0;
  r.slope = slope;
  r.b = g.b;
  r.A = g.A;
  r.status = g.status;
  if (!(sv_residual > 0.0) || !emc2_isfinite(v0)) {
    r.status = BAwLTimeStatus::invalid;
    return r;
  }
  r.sv = sv_residual;
  r.inv_sv = 1.0 / sv_residual;
  if (g.status != BAwLTimeStatus::valid &&
      g.status != BAwLTimeStatus::point_start &&
      g.status != BAwLTimeStatus::infinite) {
    return r;
  }
  r.c0 = -g.C1 * r.inv_sv;
  r.m = g.C2 * r.inv_sv;
  r.span = r.m * g.A;
  if (g.status != BAwLTimeStatus::infinite) {
    // Same Jacobian as bawl_threshold_terms(): k^2 E / G^2, limit 1/tau^2.
    if (g.k <= BAWL_K_EPS) {
      r.jacobian = 1.0 / (g.tau * g.tau);
      r.log_jacobian = -2.0 * std::log(g.tau);
    } else {
      const double G_safe = clamp_pos(g.G, 1e-300);
      r.jacobian = (g.k * g.k * g.E) / (G_safe * G_safe);
      r.log_jacobian =
        2.0 * std::log(g.k) - g.k * g.tau - 2.0 * std::log(G_safe);
    }
  }
  return r;
}

// log q(z) = log Phi(v_q / sv): the row's positivity constant under joint
// positive drifts, evaluated with the shared log-tail helper.
inline double bawl_prepared_log_q(const BAwLPreparedRow& r, double vq) {
  if (r.status == BAwLTimeStatus::invalid) return R_NegInf;
  return pnorm_log_direct(vq * r.inv_sv, true);
}

// Unrestricted log PDF f0(t | v_q).  Natural attempt mirrors
// ba_natural_pdf() with denom = 1 under BA_ACCEPT_RAW (a saturated tail that
// rounds to zero is a genuine -Inf, matching the raw kernels); otherwise the
// authoritative log branch mirrors log_ba_pdf() with log_denom = 0.
inline double bawl_prepared_log_pdf(const BAwLPreparedRow& r, double vq) {
  if (r.status != BAwLTimeStatus::valid &&
      r.status != BAwLTimeStatus::point_start)
    return R_NegInf;
  const double c = r.c0 + vq * r.inv_sv;
  if (!(c > R_NegInf)) return R_NegInf;

  const bool point = r.status == BAwLTimeStatus::point_start ||
    r.span < BAWL_NATURAL_MIN_SPAN;
  if (point) {
    const bool true_point = r.status == BAwLTimeStatus::point_start;
    const double z = true_point ? c : c + 0.5 * r.span;
    const double eff_b = true_point ? r.b : r.b - 0.5 * r.A;
    if (!emc2_isfinite(z)) return R_NegInf;
    const double pdf = r.jacobian * eff_b * dnormP(z) * r.inv_sv;
    if (R_FINITE(pdf) && pdf > 0.0) return std::log(pdf);
    // Tail underflow: the log form stays representable.
    return r.log_jacobian + std::log(eff_b) + dnormP(z, 0.0, 1.0, true) -
      std::log(r.sv);
  }

  const double hi = c + r.span;
  if (emc2_isfinite(hi)) {
    const double term1 = (r.b * r.m + c) *
      (pnorm_std(hi, true, false) - pnorm_std(c, true, false));
    const double term2 = dnormP(hi) - dnormP(c);
    const double bracket = term1 + term2;
    if (bracket > 0.0) {
      const double scale = std::fabs(term1) + std::fabs(term2);
      if (bracket > BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) {
        const double denom_scale = r.A * r.sv * r.m * r.m;
        const double pdf = r.jacobian * bracket / denom_scale;
        if (R_FINITE(pdf) && pdf > 0.0) return std::log(pdf);
      }
    }
  }

  // Authoritative log-space branch (log_ba_pdf with log_denom = 0).
  if (r.span > BAWL_LOG_MIN_SPAN && emc2_isfinite(r.span)) {
    const double hi2 = c + r.span;
    const double log_dphi = log_normal_interval(c, hi2);
    if (!ISNAN(log_dphi)) {
      const signed_log lead = signed_log_product(r.b * r.m + c, log_dphi);
      const signed_log phi_diff = signed_log_sub(
        make_signed_log(dnormP(hi2, 0.0, 1.0, true), 1),
        make_signed_log(dnormP(c, 0.0, 1.0, true), 1));
      const signed_log bracket = signed_log_add(lead, phi_diff);
      const double max_term = std::max(lead.log_abs, phi_diff.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWL_LOG_BRACKET_MIN) {
        return r.log_jacobian - std::log(r.A) - std::log(r.sv) -
          2.0 * std::log(r.m) + bracket.log_abs;
      }
      return r.log_jacobian - std::log(r.A) - std::log(r.sv) -
        std::log(r.m) + std::log(r.b - 0.5 * r.A) + log_dphi;
    }
  }
  return r.log_jacobian + std::log(r.b - 0.5 * r.A) +
    dnormP(c + 0.5 * r.span, 0.0, 1.0, true) - std::log(r.sv);
}

// Unrestricted log CDF F0(t | v_q); mirrors log_ba_cdf() with log_denom = 0,
// with the same guarded natural fast path.  The t = +Inf and point-mass
// branches share the m = 0 limit Phi(c).
inline double bawl_prepared_log_cdf(const BAwLPreparedRow& r, double vq) {
  if (r.status == BAwLTimeStatus::invalid ||
      r.status == BAwLTimeStatus::not_started)
    return R_NegInf;
  const double c = r.c0 + vq * r.inv_sv;
  if (!(c > R_NegInf)) return R_NegInf;

  if (r.status == BAwLTimeStatus::point_start ||
      r.status == BAwLTimeStatus::infinite || !(r.m > 0.0)) {
    const double out = pnorm_log_direct(c, true);
    return out >= 0.0 ? 0.0 : out;
  }

  if (r.span < BAWL_NATURAL_MIN_SPAN) {
    const double out = pnorm_log_direct(c + 0.5 * r.span, true);
    return out >= 0.0 ? 0.0 : out;
  }

  const double hi = c + r.span;
  // Natural attempt: H(z) = z Phi(z) + phi(z) endpoints, accepted while the
  // subtraction retains relative accuracy (same guard as ba_natural_cdf).
  if (emc2_isfinite(hi)) {
    const double H_hi = hi * pnorm_std(hi, true, false) + dnormP(hi);
    const double H_lo = c * pnorm_std(c, true, false) + dnormP(c);
    const double integral = H_hi - H_lo;
    if (integral > 0.0) {
      const double scale = std::fabs(H_hi) + std::fabs(H_lo);
      if (integral > BAWL_NATURAL_REL_TOL * std::max(1.0, scale)) {
        const double cdf = integral / r.span;
        if (R_FINITE(cdf) && cdf > 0.0)
          return cdf >= 1.0 ? 0.0 : std::log(cdf);
      }
    } else if (integral == 0.0 && c > BAWL_NATURAL_Z_MAX) {
      return 0.0;  // saturated interval: CDF is 1
    }
  }

  if (r.span > BAWL_LOG_MIN_SPAN && emc2_isfinite(r.span)) {
    const double log_integral = log_normal_phi_integral(c, c + r.span);
    const double out = log_integral - std::log(r.span);
    if (!ISNAN(out)) return out >= 0.0 ? 0.0 : out;
  }
  const double out = pnorm_log_direct(c + 0.5 * r.span, true);
  return out >= 0.0 ? 0.0 : out;
}

// Log survivor.  joint_positive = false: log(1 - F0) (unrestricted-drift
// mode).  joint_positive = true: the unnormalised weighted survivor
// log(q - F0) used by the fused no-clock fallback; the caller subtracts the
// aggregated positivity normaliser exactly once per trial.
inline double bawl_prepared_log_survivor(const BAwLPreparedRow& r, double vq,
                                         bool joint_positive) {
  if (r.status == BAwLTimeStatus::invalid) return R_NaN;
  if (r.status == BAwLTimeStatus::not_started)
    return joint_positive ? bawl_prepared_log_q(r, vq) : 0.0;

  const double log_cdf = bawl_prepared_log_cdf(r, vq);
  if (!joint_positive) {
    if (log_cdf == R_NegInf) return 0.0;
    if (log_cdf >= 0.0) return R_NegInf;
    return log1m_exp(log_cdf);
  }
  const double log_q = bawl_prepared_log_q(r, vq);
  if (log_cdf == R_NegInf) return log_q;
  if (log_cdf >= log_q) return R_NegInf;
  return log_diff_exp(log_q, log_cdf);
}

#endif

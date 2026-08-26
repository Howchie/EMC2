#include "model_BAwR.h"
#include <algorithm>
#include "race_contract.h"
#include "utility_functions.h"
#include "col_registry.h"

double bawr_sat_time(const BawrGeom& g, double b_minus_z) {
  if (!(b_minus_z > 0.0)) return 0.0;
  const double log_u = (std::log1p(g.pw) + std::log(b_minus_z) -
                        g.log_kappa - std::log(g.pw)) / (g.pw + 1.0);
  return std::exp(log_u);
}

double bawr_log_critical_launch(const BawrGeom& g, double s) {
  if (!(s > 0.0)) return R_NegInf;
  return g.log_kappa + g.pw * std::log(s);
}
double bawr_critical_launch(const BawrGeom& g, double s) {
  const double l = bawr_log_critical_launch(g, s);
  return (l > R_NegInf) ? std::exp(l) : 0.0;
}

BawrGeom bawr_geometry(double A, double b, double kappa, double pw) {
  BawrGeom g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(kappa) ||
      !emc2_isfinite(pw)) return g;
  if (!(b > 0.0) || !(A >= 0.0) || !(b >= A) || !(kappa >= 0.0)) return g;
  // p = 0 is not a member of the family: the drive would lose a constant, so
  // the model degenerates to an LBA with a shifted launch and T_max diverges.
  if (!(pw > 0.0)) return g;
  g.ok = true; g.b = b; g.A = A; g.kappa = kappa; g.pw = pw;
  g.kappa_zero = (kappa <= BAWR_K_EPS);
  if (g.kappa_zero) {
    // No decay: the trace is z + V u, nothing saturates, exactly an LBA.
    g.log_kappa = R_NegInf;
    g.T_max = R_PosInf;
    g.T_sat_A = R_PosInf;
    g.V_c0 = 0.0;
    return g;
  }
  g.log_kappa = std::log(kappa);
  g.T_max = bawr_sat_time(g, b);
  g.T_sat_A = (A > 0.0) ? bawr_sat_time(g, b - A) : g.T_max;
  g.V_c0 = bawr_critical_launch(g, g.T_max);
  return g;
}


BawrAtU bawr_at_u(const BawrGeom& g, double u) {
  BawrAtU s;
  if (!g.ok || !(u > 0.0) || ISNAN(u)) return s;
  s.ok = true;
  const bool inf_u = (u == R_PosInf);
  s.q = u;

  if (g.kappa_zero) {
    s.Z = g.A;
    s.c = 0.0;
    // At u = Inf the required launch is zero, so the LBA finishes with
    // probability P(V > 0) -- the survivor evaluated at w = 0.
    s.w_hi = inf_u ? 0.0 : g.b / u;
    s.w_lo = inf_u ? 0.0 : (g.b - s.Z) / u;
    s.s_lo = g.T_sat_A;
    s.s_hi = g.T_max;
    return s;
  }

  s.s_lo = g.T_sat_A;
  s.s_hi = inf_u ? g.T_max : std::fmin(u, g.T_max);

  if (inf_u || u >= g.T_max) {
    s.Z = 0.0;
    s.saturated = true;
    s.partial = (g.A > 0.0);
    s.w_hi = g.V_c0;
    s.w_lo = g.V_c0;
    s.c = g.V_c0;
    return s;
  }

  // z*(u) = b [1 - (u / T_max)^(p+1)], via expm1 so that it stays accurate
  // (and correctly signed) as u approaches the endpoint, where the whole
  // live contribution is carried by this vanishing width.
  const double zstar =
    -g.b * std::expm1((g.pw + 1.0) * (std::log(u) - std::log(g.T_max)));
  s.saturated = !(zstar > 0.0);
  s.Z = std::fmin(std::fmax(zstar, 0.0), g.A);
  s.partial = (s.Z < g.A);

  s.c = std::exp(g.log_kappa + g.pw * std::log(u));
  s.w_hi = g.b / u + s.c / (g.pw + 1.0);
  // At an unclipped upper limit (Z = z*(u)) the required launch is exactly the
  // critical launch: (b - z*)/u = kappa p u^p/(p+1), so w_lo = kappa u^p = c
  // identically.  Setting it rather than recomputing it keeps the density's
  // integrand (w - c) vanishing EXACTLY at w_lo, which is what makes the
  // endpoint shutdown clean; the clipped branch has no such identity.
  s.w_lo = (s.Z < g.A) ? s.c : (g.b - g.A) / u + s.c / (g.pw + 1.0);
  return s;
}

double bawr_log_frozen_quad(const BawrGeom& g, double s_lo, double s_hi,
                                   double p1, double p2, bool logn,
                                   bool posdrift, double delta) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (logn && delta != 0.0 && !split_lognormal_shape_params(p1, p2, delta, h))
    return R_NegInf;
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_vc = bawr_log_critical_launch(g, s);
    const double log_jac = g.log_kappa + std::log(g.pw) + g.pw * std::log(s);
    const double log_surv = logn
      ? (delta == 0.0
           ? pnorm_log_direct((log_vc - p1) / p2, false)
           : log_split_lognormal_survivor(std::exp(log_vc), h))
      : pnorm_log_direct((p1 - std::exp(log_vc)) / p2, true);
    return log_jac + log_surv;
  };
  (void)posdrift;
  double mid = s_lo;
  const double log_ref = logn ? p1 : ((p1 > 0.0) ? std::log(p1) : R_NegInf);
  if (log_ref > R_NegInf) mid = std::exp((log_ref - g.log_kappa) / g.pw);
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES);
}

double bawr_log_frozen_normal(const BawrGeom& g, double s_lo,
                                     double s_hi, double v, double sv) {
  return bawr_log_frozen_quad(g, s_lo, s_hi, v, sv, false, true, 0.0);
}

double bawr_log_frozen_logn(const BawrGeom& g, double s_lo, double s_hi,
                                   double mu, double sigma, double delta) {
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  const double w_a = bawr_critical_launch(g, s_lo);
  const double w_b = bawr_critical_launch(g, s_hi);
  // w_a == 0 (A == b, so the highest start point is already at threshold) has
  // a finite integral that the power primitive cannot represent; quadrature.
  if (!(w_b > w_a) || !emc2_isfinite(w_b) || !(w_a > 0.0))
    return bawr_log_frozen_quad(g, s_lo, s_hi, mu, sigma, true, true, delta);

  // kappa^(-1/p) [ Lambda_m(w_a) - Lambda_m(w_b) ], m = -(p + 1)/p.
  // The split primitive carries an exact delta == 0 branch delegating to the
  // plain power stop-loss, so the lognormal numbers are bit-for-bit unchanged.
  const double m = -(g.pw + 1.0) / g.pw;
  const double lam_a =
    log_split_lognormal_power_stoploss(w_a, mu, sigma, delta, m);
  const double lam_b =
    log_split_lognormal_power_stoploss(w_b, mu, sigma, delta, m);
  if (lam_a - lam_b > BAWR_MIN_LOG_GAP) {
    const double d = log_diff_exp(lam_a, lam_b);
    if (emc2_isfinite(d)) return d - g.log_kappa / g.pw;
  }
  return bawr_log_frozen_quad(g, s_lo, s_hi, mu, sigma, true, true, delta);
}

double log_bawr_cdf_normal(double u, const BawrGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);

  if (g.A <= BAWR_A_EPS) {
    const double out = pnorm_log_direct((v - s.w_hi) / sv, true) - log_denom;
    return std::fmin(out, 0.0);
  }

  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double cl = (v - s.w_hi) / sv;
    const double ch = (v - s.w_lo) / sv;
    if (!emc2_isfinite(cl) && !emc2_isfinite(ch)) return R_NegInf;
    const double span = ch - cl;
    if (span > BAWR_MIN_SPAN && emc2_isfinite(span)) {
      const double li = log_normal_phi_integral(cl, ch);
      log_live = ISNAN(li)
        ? std::log(s.Z) + pnorm_log_direct(0.5 * (cl + ch), true)
        : li + std::log(sv) + std::log(s.q);
    } else {
      log_live = std::log(s.Z) + pnorm_log_direct(0.5 * (cl + ch), true);
    }
  }

  const double log_frozen = s.partial
    ? bawr_log_frozen_normal(g, s.s_lo, s.s_hi, v, sv) : R_NegInf;

  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A) -
    log_denom;
  if (ISNAN(out)) return R_NegInf;
  return std::fmin(out, 0.0);
}

double log_bawr_cdf_logn(double u, const BawrGeom& g, double mu,
                                double sigma, double delta) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const auto log_gbar = [&](double w) -> double {
    if (!(w > 0.0)) return 0.0;
    if (!emc2_isfinite(w)) return R_NegInf;
    return delta == 0.0
      ? pnorm_log_direct((mu - std::log(w)) / sigma, true)
      : log_split_lognormal_survivor(w, mu, sigma, delta);
  };

  if (g.A <= BAWR_A_EPS) return std::fmin(log_gbar(s.w_hi), 0.0);

  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double lc_lo = log_split_lognormal_stoploss(s.w_lo, mu, sigma, delta);
    const double lc_hi = log_split_lognormal_stoploss(s.w_hi, mu, sigma, delta);
    if (lc_lo - lc_hi > BAWR_MIN_LOG_GAP) {
      const double ld = log_diff_exp(lc_lo, lc_hi);
      log_live = ISNAN(ld) ? R_NegInf : std::log(s.q) + ld;
    }
    if (!(log_live > R_NegInf)) {
      log_live = std::log(s.Z) + log_gbar(0.5 * (s.w_lo + s.w_hi));
    }
  }

  const double log_frozen = s.partial
    ? bawr_log_frozen_logn(g, s.s_lo, s.s_hi, mu, sigma, delta) : R_NegInf;

  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  if (ISNAN(out)) return R_NegInf;
  return std::fmin(out, 0.0);
}

double bawr_log_frozen_surv_quad(const BawrGeom& g, double s_lo,
                                        double s_hi, double p1, double p2,
                                        bool logn, bool posdrift, double delta) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_vc = bawr_log_critical_launch(g, s);
    const double log_jac = g.log_kappa + std::log(g.pw) + g.pw * std::log(s);
    const double log_cdf = logn
      ? (delta == 0.0
           ? pnorm_log_direct((log_vc - p1) / p2, true)
           : log_split_lognormal_cdf(std::exp(log_vc), p1, p2, delta))
      : (posdrift
           ? log_normal_cdf_positive_raw(std::exp(log_vc), p1, p2)
           : pnorm_log_direct((std::exp(log_vc) - p1) / p2, true));
    return log_jac + log_cdf;
  };
  const double log_ref = logn ? p1 : ((p1 > 0.0) ? std::log(p1) : R_NegInf);
  const double mid = (log_ref > R_NegInf)
    ? std::exp((log_ref - g.log_kappa) / g.pw) : s_lo;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES);
}

double bawr_log_frozen_surv_normal(const BawrGeom& g, double s_lo,
                                          double s_hi, double v, double sv,
                                          bool posdrift) {
  return bawr_log_frozen_surv_quad(g, s_lo, s_hi, v, sv, false, posdrift, 0.0);
}
double bawr_log_frozen_surv_logn(const BawrGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma,
                                        double delta) {
  return bawr_log_frozen_surv_quad(g, s_lo, s_hi, mu, sigma, true, false, delta);
}

double log_bawr_surv_normal(double u, const BawrGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  if (u == R_PosInf && g.kappa_zero)
    return std::fmin(log_normal_cdf_positive(s.w_hi, v, sv, posdrift,
                                             denom_floor), 0.0);
  if (g.A <= BAWR_A_EPS)
    return std::fmin(log_normal_cdf_positive(s.w_hi, v, sv, posdrift,
                                             denom_floor), 0.0);
  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double li = log_normal_phi_integral_positive_raw(
      (s.w_lo - v) / sv, (s.w_hi - v) / sv, v, sv, posdrift);
    if (li > R_NegInf) log_live = li + std::log(sv) + std::log(s.q);
  }
  const double log_frozen = s.partial
    ? bawr_log_frozen_surv_normal(g, s.s_lo, s.s_hi, v, sv, posdrift)
    : R_NegInf;
  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A) - log_denom;
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double log_bawr_surv_logn(double u, const BawrGeom& g, double mu,
                                  double sigma, double delta) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const auto log_g = [&](double w) {
    return (w > 0.0 && emc2_isfinite(w))
      ? log_split_lognormal_cdf(w, mu, sigma, delta) : R_NegInf;
  };
  if (u == R_PosInf && g.kappa_zero) return std::fmin(log_g(s.w_hi), 0.0);
  if (g.A <= BAWR_A_EPS) return std::fmin(log_g(s.w_hi), 0.0);
  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double pa = log_split_lognormal_put(s.w_hi, mu, sigma, delta);
    const double pb = log_split_lognormal_put(s.w_lo, mu, sigma, delta);
    if (pa - pb > BAWR_MIN_LOG_GAP)
      log_live = std::log(s.q) + log_diff_exp(pa, pb);
    if (!(log_live > R_NegInf))
      log_live = std::log(s.Z) + log_g(0.5 * (s.w_hi + s.w_lo));
  }
  const double log_frozen = s.partial
    ? bawr_log_frozen_surv_logn(g, s.s_lo, s.s_hi, mu, sigma, delta) : R_NegInf;
  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double log_bawr_pdf_normal(double u, const BawrGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);

  if (g.A <= BAWR_A_EPS) {
    const double wgt = s.w_hi - s.c;
    if (!(wgt > 0.0)) return R_NegInf;
    return dnormP((v - s.w_hi) / sv, 0.0, 1.0, true) + std::log(wgt) -
      std::log(sv) - std::log(s.q) - log_denom;
  }
  const double cl = (v - s.w_hi) / sv, ch = (v - s.w_lo) / sv;
  const double span = ch - cl;
  const double log_scale = -std::log(g.A) - log_denom;
  if (span > BAWR_MIN_SPAN && emc2_isfinite(span)) {
    const double log_dphi = log_normal_interval(cl, ch);
    if (!ISNAN(log_dphi)) {
      // int w g dw = v DPhi + sv (phi(ch) - phi(cl)); then subtract c DPhi.
      const signed_log dphi = signed_log_sub(
        make_signed_log(dnormP(ch, 0.0, 1.0, true), 1),
        make_signed_log(dnormP(cl, 0.0, 1.0, true), 1));
      // v and c are kept as separate signed-log terms rather than differenced
      // in natural space: they are the same order of magnitude near the
      // endpoint, where their difference carries the whole density.
      const signed_log v_term = signed_log_product(v, log_dphi);
      const signed_log sv_term = (dphi.sign != 0)
        ? make_signed_log(std::log(sv) + dphi.log_abs, dphi.sign)
        : make_signed_log(R_NegInf, 0);
      const signed_log lead = signed_log_add(v_term, sv_term);
      const signed_log c_term = signed_log_product(-s.c, log_dphi);
      const signed_log bracket = signed_log_add(lead, c_term);
      const double max_term = std::fmax(lead.log_abs, c_term.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWR_LOG_BRACKET_MIN) {
        return bracket.log_abs + log_scale;
      }
      const double w_mid = 0.5 * (s.w_hi + s.w_lo);
      const double wgt = w_mid - s.c;
      if (!(wgt > 0.0)) return R_NegInf;
      return std::log(wgt) + log_dphi + log_scale;
    }
  }
  const double w_mid = 0.5 * (s.w_hi + s.w_lo);
  const double wgt = w_mid - s.c;
  if (!(wgt > 0.0)) return R_NegInf;
  return std::log(s.Z) + std::log(wgt) - std::log(sv) - std::log(s.q) +
    dnormP(0.5 * (cl + ch), 0.0, 1.0, true) + log_scale;
}

double log_bawr_pdf_logn(double u, const BawrGeom& g, double mu,
                                double sigma, double delta) {
  // Split launch: the density bracket assembles from the split survivor and
  // first-partial-moment primitives; the plain-lognormal d1/d2 reduction
  // below has no split analogue.  delta == 0 never enters this block.
  if (delta != 0.0) {
    if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf)
      return R_NegInf;
    const BawrAtU s = bawr_at_u(g, u);
    if (!s.ok || s.saturated) return R_NegInf;
    if (g.A <= BAWR_A_EPS) {
      const double wgt = s.w_hi - s.c;
      return wgt > 0.0
        ? log_split_lognormal_density(s.w_hi, mu, sigma, delta) +
            std::log(wgt) - std::log(s.q)
        : R_NegInf;
    }
    const double lSlo = log_split_lognormal_survivor(s.w_lo, mu, sigma, delta);
    const double lShi = log_split_lognormal_survivor(s.w_hi, mu, sigma, delta);
    const double lMlo = log_split_lognormal_first_partial_moment(
      s.w_lo, mu, sigma, delta);
    const double lMhi = log_split_lognormal_first_partial_moment(
      s.w_hi, mu, sigma, delta);
    const signed_log t1 = signed_log_sub(
      make_signed_log(lMlo, 1), make_signed_log(lMhi, 1));
    const signed_log sdiff = signed_log_sub(
      make_signed_log(lSlo, 1), make_signed_log(lShi, 1));
    signed_log t2 = signed_log_product(s.c, sdiff.log_abs);
    t2.sign *= sdiff.sign;
    const signed_log bracket = signed_log_sub(t1, t2);
    if (bracket.sign > 0 && emc2_isfinite(bracket.log_abs))
      return bracket.log_abs - std::log(g.A);
    const double w_mid = 0.5 * (s.w_hi + s.w_lo);
    const double wgt = w_mid - s.c;
    return wgt > 0.0
      ? std::log(wgt) +
          log_split_lognormal_density(w_mid, mu, sigma, delta) -
          std::log(s.q) - std::log(g.A)
      : R_NegInf;
  }
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;

  if (g.A <= BAWR_A_EPS) {
    const double wgt = s.w_hi - s.c;
    if (!(wgt > 0.0)) return R_NegInf;
    return dlnorm_std(s.w_hi, mu, sigma, true) + std::log(wgt) -
      std::log(s.q);
  }

  const double d1_lo = (mu + sigma * sigma - std::log(s.w_lo)) / sigma;
  const double d1_hi = (mu + sigma * sigma - std::log(s.w_hi)) / sigma;
  const double d2_lo = (mu - std::log(s.w_lo)) / sigma;
  const double d2_hi = (mu - std::log(s.w_hi)) / sigma;
  const double log_prob = log_normal_interval(d2_hi, d2_lo);
  const double log_scale = -std::log(g.A);

  if (!ISNAN(log_prob)) {
    const double log_pexp = log_normal_interval(d1_hi, d1_lo);
    if (!ISNAN(log_pexp)) {
      const signed_log t1 = make_signed_log(
        mu + 0.5 * sigma * sigma + log_pexp, 1);
      const signed_log t2 = signed_log_product(s.c, log_prob);
      const signed_log bracket = signed_log_sub(t1, t2);
      const double max_term = std::fmax(t1.log_abs, t2.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWR_LOG_BRACKET_MIN) {
        return bracket.log_abs + log_scale;
      }
    }
  }
  const double w_mid = 0.5 * (s.w_hi + s.w_lo);
  const double wgt = w_mid - s.c;
  if (!(wgt > 0.0)) return R_NegInf;
  if (!ISNAN(log_prob) && log_prob > R_NegInf)
    return std::log(wgt) + log_prob + log_scale;
  return std::log(s.Z) + std::log(wgt) + dlnorm_std(w_mid, mu, sigma, true) -
    std::log(s.q) + log_scale;
}

bool bawr_natural_cdf_normal(double u, const BawrGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &cdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  const auto accept = [accept_mode](double &p) {
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(p);
    if (accept_mode == BA_ACCEPT_RAW) return p < 1.0 - 1e-8;
    if (p > 1.0) p = 1.0;
    return true;
  };

  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) { cdf = 0.0; return lenient; }
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.q > 0.0)) { cdf = 0.0; return lenient; }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;

  if (g.A <= BAWR_A_EPS) {
    const double z_norm = (v - s.w_hi) / sv;
    if (!emc2_isfinite(z_norm)) return false;
    if (!lenient && std::fabs(z_norm) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(z_norm, true, false) / denom;
  } else {
    double live_part = 0.0;
    if (s.Z > 0.0) {
      const double cl = (v - s.w_hi) / sv;
      const double ch = (v - s.w_lo) / sv;
      const double span = ch - cl;
      if (span > BAWR_MIN_SPAN && emc2_isfinite(span)) {
        if (!lenient && !natural_normal_interval_safe(cl, ch)) return false;
        double scale = 0.0;
        const double integral = normal_phi_integral_nat(cl, ch, scale);
        if (!(integral > 0.0)) {
          if (!lenient) return false;
          live_part = 0.0;
        } else {
          if (!lenient &&
              integral <= BAWL_NATURAL_REL_TOL * std::fmax(1.0, scale))
            return false;
          live_part = integral * sv * s.q;
        }
      } else {
        const double mid_z = 0.5 * (cl + ch);
        if (!lenient && std::fabs(mid_z) > BAWL_NATURAL_Z_MAX) return false;
        live_part = s.Z * pnorm_std(mid_z, true, false);
      }
    }
    double frozen_part = 0.0;
    if (s.partial) {
      const double log_frozen = bawr_log_frozen_normal(g, s.s_lo, s.s_hi, v, sv);
      if (log_frozen > R_NegInf) frozen_part = std::exp(log_frozen);
    }
    cdf = (live_part + frozen_part) / (g.A * denom);
  }

  if (!R_FINITE(cdf)) return false;
  if (cdf <= 0.0) {
    if (!lenient) return false;
    cdf = 0.0;
    return true;
  }
  return accept(cdf);
}

bool bawr_natural_cdf_logn(double u, const BawrGeom& g, double mu,
                                  double sigma, int accept_mode, double &cdf,
                                  double delta) {
  // No natural fast path for the split launch yet: fall back to the log
  // kernel, whose split branch is exact at every acceptance level.
  if (delta != 0.0) {
    const double lp = log_bawr_cdf_logn(u, g, mu, sigma, delta);
    cdf = lp > R_NegInf ? std::exp(lp) : 0.0;
    return emc2_isfinite(cdf);
  }
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  const auto accept = [accept_mode](double &p) {
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(p);
    if (accept_mode == BA_ACCEPT_RAW) return p < 1.0 - 1e-8;
    if (p > 1.0) p = 1.0;
    return true;
  };

  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) { cdf = 0.0; return lenient; }
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.w_hi > 0.0) || !(s.w_lo > 0.0)) { cdf = 0.0; return lenient; }

  if (g.A <= BAWR_A_EPS) {
    const double z_norm = (mu - std::log(s.w_hi)) / sigma;
    if (!emc2_isfinite(z_norm)) return false;
    if (!lenient && std::fabs(z_norm) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(z_norm, true, false);
  } else {
    double live_part = 0.0;
    if (s.Z > 0.0) {
      double scale = 0.0;
      const double diff = lognormal_stoploss_interval_nat(s.w_lo, s.w_hi, mu,
                                                          sigma, scale);
      if (diff > 0.0 && scale > 0.0 && diff / scale > BAWR_MIN_LOG_GAP) {
        if (!lenient && diff <= BAWL_NATURAL_REL_TOL * std::fmax(1.0, scale))
          return false;
        live_part = s.q * diff;
      } else {
        const double w_mid = 0.5 * (s.w_lo + s.w_hi);
        if (!(w_mid > 0.0)) return false;
        const double z_mid = (mu - std::log(w_mid)) / sigma;
        if (!lenient && std::fabs(z_mid) > BAWL_NATURAL_Z_MAX) return false;
        live_part = s.Z * pnorm_std(z_mid, true, false);
      }
    }
    double frozen_part = 0.0;
    if (s.partial) {
      const double log_frozen = bawr_log_frozen_logn(g, s.s_lo, s.s_hi, mu,
                                                     sigma);
      if (log_frozen > R_NegInf) frozen_part = std::exp(log_frozen);
    }
    cdf = (live_part + frozen_part) / g.A;
  }

  if (!R_FINITE(cdf)) return false;
  if (cdf <= 0.0) {
    if (!lenient) return false;
    cdf = 0.0;
    return true;
  }
  return accept(cdf);
}

bool bawr_natural_pdf_normal(double u, const BawrGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &pdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) {
    pdf = 0.0;
    return lenient;
  }
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.q > 0.0)) { pdf = 0.0; return lenient; }
  if (s.saturated) { pdf = 0.0; return true; }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;

  if (g.A <= BAWR_A_EPS) {
    const double wgt = s.w_hi - s.c;
    if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
    const double z_norm = (v - s.w_hi) / sv;
    if (!emc2_isfinite(z_norm)) return false;
    if (!lenient && std::fabs(z_norm) > BAWL_NATURAL_Z_MAX) return false;
    pdf = wgt * dnormP(z_norm) / (sv * s.q * denom);
  } else {
    const double cl = (v - s.w_hi) / sv;
    const double ch = (v - s.w_lo) / sv;
    const double span = ch - cl;
    const double w_mid = 0.5 * (s.w_hi + s.w_lo);
    if (span > BAWR_MIN_SPAN && emc2_isfinite(span)) {
      if (!lenient && !natural_normal_interval_safe(cl, ch)) return false;
      double scale = 0.0;
      const double dphi = normal_interval_nat(cl, ch, scale);
      if (dphi > 0.0 && scale > 0.0 && dphi / scale > 1e-6) {
        const double term1 = (v - s.c) * dphi;
        const double term2 = sv * (dnormP(ch) - dnormP(cl));
        const double bracket = term1 + term2;
        const double max_scale = std::fmax(std::fabs(term1), std::fabs(term2));
        if (bracket > 0.0 && bracket > 1e-6 * max_scale) {
          pdf = bracket / (g.A * denom);
        } else {
          const double wgt = w_mid - s.c;
          if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
          pdf = wgt * dphi / (g.A * denom);
        }
      } else {
        const double wgt = w_mid - s.c;
        if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
        const double mid_z = 0.5 * (cl + ch);
        if (!lenient && std::fabs(mid_z) > BAWL_NATURAL_Z_MAX) return false;
        pdf = s.Z * wgt * dnormP(mid_z) / (g.A * sv * s.q * denom);
      }
    } else {
      const double wgt = w_mid - s.c;
      if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
      const double mid_z = 0.5 * (cl + ch);
      if (!lenient && std::fabs(mid_z) > BAWL_NATURAL_Z_MAX) return false;
      pdf = s.Z * wgt * dnormP(mid_z) / (g.A * sv * s.q * denom);
    }
  }

  if (!R_FINITE(pdf)) return false;
  if (pdf <= 0.0) {
    if (!lenient) return false;
    pdf = 0.0;
  }
  return true;
}

bool bawr_natural_pdf_logn(double u, const BawrGeom& g, double mu,
                                  double sigma, int accept_mode, double &pdf,
                                  double delta) {
  // As the CDF: no natural fast path for the split launch yet.
  if (delta != 0.0) {
    const double lp = log_bawr_pdf_logn(u, g, mu, sigma, delta);
    pdf = lp > R_NegInf ? std::exp(lp) : 0.0;
    return emc2_isfinite(pdf);
  }
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) {
    pdf = 0.0;
    return lenient;
  }
  const BawrAtU s = bawr_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.w_hi > 0.0) || !(s.w_lo > 0.0)) { pdf = 0.0; return lenient; }
  if (s.saturated) { pdf = 0.0; return true; }

  if (g.A <= BAWR_A_EPS) {
    const double wgt = s.w_hi - s.c;
    if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
    pdf = dlnorm_std(s.w_hi, mu, sigma, false) * wgt / s.q;
  } else {
    const double d1_lo = (mu + sigma * sigma - std::log(s.w_lo)) / sigma;
    const double d1_hi = (mu + sigma * sigma - std::log(s.w_hi)) / sigma;
    const double d2_lo = (mu - std::log(s.w_lo)) / sigma;
    const double d2_hi = (mu - std::log(s.w_hi)) / sigma;

    double scale_prob = 0.0;
    const double prob = normal_interval_nat(d2_hi, d2_lo, scale_prob);
    double scale_pexp = 0.0;
    const double pexp = normal_interval_nat(d1_hi, d1_lo, scale_pexp);

    const double M = std::exp(mu + 0.5 * sigma * sigma);
    const double t1 = M * pexp;
    const double t2 = s.c * prob;
    const double bracket = t1 - t2;
    const double w_mid = 0.5 * (s.w_hi + s.w_lo);

    if (bracket > 0.0 && prob > 0.0 && pexp > 0.0 &&
        (prob / scale_prob > 1e-6) && (pexp / scale_pexp > 1e-6)) {
      const double max_scale = std::fmax(std::fabs(t1), std::fabs(t2));
      if (!lenient &&
          bracket <= BAWL_NATURAL_REL_TOL * std::fmax(1.0, max_scale))
        return false;
      pdf = bracket / g.A;
    } else {
      const double wgt = w_mid - s.c;
      if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
      if (prob > 0.0 && scale_prob > 0.0 && prob / scale_prob > 1e-6) {
        pdf = wgt * prob / g.A;
      } else {
        pdf = s.Z * wgt * dlnorm_std(w_mid, mu, sigma, false) / (g.A * s.q);
      }
    }
  }

  if (!R_FINITE(pdf)) return false;
  if (pdf <= 0.0) {
    if (!lenient) return false;
    pdf = 0.0;
  }
  return true;
}

bool ba_natural_cdf_bawr(double u, double A, double b, double p1,
                                double p2, double kappa, double pw, int launch,
                                bool posdrift, double denom_floor,
                                int accept_mode, double &cdf, double delta) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  if (launch == BAWR_LAUNCH_LOGNORMAL ||
      launch == BAWR_LAUNCH_SPLITLOGNORMAL)
    return bawr_natural_cdf_logn(u, g, p1, p2, accept_mode, cdf, delta);
  return bawr_natural_cdf_normal(u, g, p1, p2, posdrift, denom_floor,
                                 accept_mode, cdf);
}

bool ba_natural_pdf_bawr(double u, double A, double b, double p1,
                                double p2, double kappa, double pw, int launch,
                                bool posdrift, double denom_floor,
                                int accept_mode, double &pdf, double delta) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  if (launch == BAWR_LAUNCH_LOGNORMAL ||
      launch == BAWR_LAUNCH_SPLITLOGNORMAL)
    return bawr_natural_pdf_logn(u, g, p1, p2, accept_mode, pdf, delta);
  return bawr_natural_pdf_normal(u, g, p1, p2, posdrift, denom_floor,
                                 accept_mode, pdf);
}

double bawr_log_cdf(double u, double A, double b, double p1, double p2,
                           double kappa, double pw, int launch, bool posdrift,
                           double denom_floor, double delta) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  if (launch == BAWR_LAUNCH_LOGNORMAL ||
      launch == BAWR_LAUNCH_SPLITLOGNORMAL)
    return log_bawr_cdf_logn(u, g, p1, p2, delta);
  return log_bawr_cdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawr_log_surv(double u, double A, double b, double p1, double p2,
                            double kappa, double pw, int launch, bool posdrift,
                            double denom_floor, double delta) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  if (launch == BAWR_LAUNCH_LOGNORMAL ||
      launch == BAWR_LAUNCH_SPLITLOGNORMAL)
    return log_bawr_surv_logn(u, g, p1, p2, delta);
  return log_bawr_surv_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawr_log_pdf(double u, double A, double b, double p1, double p2,
                           double kappa, double pw, int launch, bool posdrift,
                           double denom_floor, double delta) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  if (launch == BAWR_LAUNCH_LOGNORMAL ||
      launch == BAWR_LAUNCH_SPLITLOGNORMAL)
    return log_bawr_pdf_logn(u, g, p1, p2, delta);
  return log_bawr_pdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawr_cdf_norm(double t, double A, double b, double p1, double p2,
                            double kappa, double pw, int launch, bool posdrift,
                            bool log_out,
                            double denom_floor, double delta) {
  double cdf;
  if (ba_natural_cdf_bawr(t, A, b, p1, p2, kappa, pw, launch, posdrift,
                          denom_floor, BA_ACCEPT_STRICT, cdf, delta))
    return log_out ? std::log(cdf) : cdf;
  return return_from_log(
    bawr_log_cdf(t, A, b, p1, p2, kappa, pw, launch, posdrift, denom_floor,
                 delta),
    log_out);
}

double bawr_pdf_norm(double t, double A, double b, double p1, double p2,
                            double kappa, double pw, int launch, bool posdrift,
                            bool log_out,
                            double denom_floor, double delta) {
  double pdf;
  if (ba_natural_pdf_bawr(t, A, b, p1, p2, kappa, pw, launch, posdrift,
                          denom_floor, BA_ACCEPT_STRICT, pdf, delta))
    return log_out ? std::log(pdf) : pdf;
  return return_from_log(
    bawr_log_pdf(t, A, b, p1, p2, kappa, pw, launch, posdrift, denom_floor,
                 delta),
    log_out);
}

double bawr_cdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double kappa, double pw,
                                      int launch, bool posdrift, double delta) {
  double cdf;
  if (ba_natural_cdf_bawr(t, A, b, p1, p2, kappa, pw, launch, posdrift,
                          BAWR_DENOM_FLOOR, BA_ACCEPT_CLAMP, cdf, delta))
    return cdf;
  const double lp = bawr_log_cdf(t, A, b, p1, p2, kappa, pw, launch, posdrift,
                                 BAWR_DENOM_FLOOR, delta);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

double bawr_pdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double kappa, double pw,
                                      int launch, bool posdrift, double delta) {
  double pdf;
  if (ba_natural_pdf_bawr(t, A, b, p1, p2, kappa, pw, launch, posdrift,
                          BAWR_DENOM_FLOOR, BA_ACCEPT_CLAMP, pdf, delta))
    return pdf;
  const double lp = bawr_log_pdf(t, A, b, p1, p2, kappa, pw, launch, posdrift,
                                 BAWR_DENOM_FLOOR, delta);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}
// --------------------------------------------------------------------------
// R-callable entry points.  As for BAwD and BAwF these bypass
// ContextForRaceModels, so the launch distribution MUST be passed explicitly;
// R/model_BAwR.R derives both this argument and the c_name suffix from one
// `drift_distribution` value so dfun/pfun cannot silently disagree with the
// sampled likelihood.
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dbawr(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector kappa,
                    NumericVector pw, int launch = 1, bool posdrift = true,
                    bool log_out = false, NumericVector delta = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawr_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(kappa, i), pick(pw, i), launch,
                           posdrift, log_out, BAWR_DENOM_FLOOR, pick(delta, i));
  return out;
}

// [[Rcpp::export]]
NumericVector pbawr(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector kappa,
                    NumericVector pw, int launch = 1, bool posdrift = true,
                    bool log_out = false, NumericVector delta = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawr_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(kappa, i), pick(pw, i), launch,
                           posdrift, log_out, BAWR_DENOM_FLOOR, pick(delta, i));
  return out;
}

// [[Rcpp::export]]
double dbawr_norm(double t, double A, double b, double p1, double p2,
                  double kappa, double pw, int launch = 1,
                  bool posdrift = true, bool log_out = false,
                  double delta = 0.0) {
  return bawr_pdf_norm(t, A, b, p1, p2, kappa, pw, launch, posdrift, log_out,
                       BAWR_DENOM_FLOOR, delta);
}

// [[Rcpp::export]]
double pbawr_norm(double t, double A, double b, double p1, double p2,
                  double kappa, double pw, int launch = 1,
                  bool posdrift = true, bool log_out = false,
                  double delta = 0.0) {
  return bawr_cdf_norm(t, A, b, p1, p2, kappa, pw, launch, posdrift, log_out,
                       BAWR_DENOM_FLOOR, delta);
}

// Right endpoint of the supported decision-time window,
// T_max = [(p+1) b / (kappa p)]^(1/(p+1)); Inf only at kappa = 0.  Unlike
// BAwF's it DOES depend on b, which is the substantive difference between the
// two models, so it is exported and reported by Ttransform: a test or a
// summary that recomputed it in R would not be testing the same code the
// likelihood uses.
// [[Rcpp::export]]
double bawr_tmax(double A, double b, double kappa, double pw) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  return g.ok ? g.T_max : NA_REAL;
}

// Vectorised bawr_tmax for the R Ttransform, which derives T_max per
// accumulator row; a per-row .Call would dominate mapped_pars()/make_data().
// [[Rcpp::export]]
NumericVector bawr_tmax_vec(NumericVector A, NumericVector b,
                            NumericVector kappa, NumericVector pw) {
  // Sized by the LONGEST argument, not by A: Ttransform passes full columns,
  // but a scalar A with vector b (a caution sweep) would otherwise silently
  // return a single value.
  const int n = std::max(std::max(A.size(), b.size()),
                         std::max(kappa.size(), pw.size()));
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawrGeom g = bawr_geometry(pick(A, i), pick(b, i), pick(kappa, i),
                                     pick(pw, i));
    out[i] = g.ok ? g.T_max : NA_REAL;
  }
  return out;
}


// The critical launch strength at the lowest start point,
// V_c(0) = kappa T_max^p.  Launches below it never reach the threshold, so
// this is the observable that sets the omission rate.
// [[Rcpp::export]]
NumericVector bawr_vcrit_vec(NumericVector A, NumericVector b,
                             NumericVector kappa, NumericVector pw) {
  // Sized by the longest argument; see bawr_tmax_vec().
  const int n = std::max(std::max(A.size(), b.size()),
                         std::max(kappa.size(), pw.size()));
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawrGeom g = bawr_geometry(pick(A, i), pick(b, i), pick(kappa, i),
                                     pick(pw, i));
    out[i] = g.ok ? g.V_c0 : NA_REAL;
  }
  return out;
}

// --------------------------------------------------------------------------
// Race adapters: scalar, raw, and truncation kernels routed through
// particle_ll.cpp's function pointers.  Moved from utils.h.
// --------------------------------------------------------------------------
namespace {

// BAwR shares BAwD's launch context field, for the same reason BAwF does: the
// models never coexist in one adapter and the "0 = normal, 1 = lognormal,
// 2 = continuous split-lognormal" contract should have exactly one
// definition.  BAwR needs no rho: its decay shape is the sampled exponent
// `p`, not a fixed kernel index.
int bawr_launch_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->bawd_launch : BAWR_LAUNCH_LOGNORMAL;
}


}  // namespace

double dbawr_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = bawr_launch_of(ctx);
  const bool split = launch == BAWR_LAUNCH_SPLITLOGNORMAL;
  const int iv = emc2col::select_index(split, emc2col::bawrsplit::mu, emc2col::bawr::v);
  const int isv = emc2col::select_index(split, emc2col::bawrsplit::sigma, emc2col::bawr::sv);
  const int iB = emc2col::select_index(split, emc2col::bawrsplit::B, emc2col::bawr::B);
  const int iA = emc2col::select_index(split, emc2col::bawrsplit::A, emc2col::bawr::A);
  const int it0 = emc2col::select_index(split, emc2col::bawrsplit::t0, emc2col::bawr::t0);
  const int ikap = emc2col::select_index(split, emc2col::bawrsplit::kappa, emc2col::bawr::kappa);
  const int ipw = emc2col::select_index(split, emc2col::bawrsplit::p, emc2col::bawr::p);
  if (R_IsNA(par[iv])) return 0.0;
  const double tt = t - par[it0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  return bawr_pdf_scalar_natural(
    tt, par[iA], par[iB] + par[iA],
    par[iv], par[isv], par[ikap], par[ipw],
    launch, ctx ? ctx->use_posdrift : true,
    split ? par[emc2col::bawrsplit::delta] : 0.0);
}

double pbawr_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = bawr_launch_of(ctx);
  const bool split = launch == BAWR_LAUNCH_SPLITLOGNORMAL;
  const int iv = emc2col::select_index(split, emc2col::bawrsplit::mu, emc2col::bawr::v);
  const int isv = emc2col::select_index(split, emc2col::bawrsplit::sigma, emc2col::bawr::sv);
  const int iB = emc2col::select_index(split, emc2col::bawrsplit::B, emc2col::bawr::B);
  const int iA = emc2col::select_index(split, emc2col::bawrsplit::A, emc2col::bawr::A);
  const int it0 = emc2col::select_index(split, emc2col::bawrsplit::t0, emc2col::bawr::t0);
  const int ikap = emc2col::select_index(split, emc2col::bawrsplit::kappa, emc2col::bawr::kappa);
  const int ipw = emc2col::select_index(split, emc2col::bawrsplit::p, emc2col::bawr::p);
  if (R_IsNA(par[iv])) return 0.0;
  const double tt = t - par[it0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf is handled inside the kernel and returns F_max, not 1.
  return bawr_cdf_scalar_natural(
    tt, par[iA], par[iB] + par[iA],
    par[iv], par[isv], par[ikap], par[ipw],
    launch, ctx ? ctx->use_posdrift : true,
    split ? par[emc2col::bawrsplit::delta] : 0.0);
}

void dbawr_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawr_launch_of(ctx);
  const bool split = launch == BAWR_LAUNCH_SPLITLOGNORMAL;
  const double* p1_ = cols[emc2col::select_index(split, emc2col::bawrsplit::mu, emc2col::bawr::v)];
  const double* p2_ =
    cols[emc2col::select_index(split, emc2col::bawrsplit::sigma, emc2col::bawr::sv)];
  const double* B_ = cols[emc2col::select_index(split, emc2col::bawrsplit::B, emc2col::bawr::B)];
  const double* A_ = cols[emc2col::select_index(split, emc2col::bawrsplit::A, emc2col::bawr::A)];
  const double* t0_ = cols[emc2col::select_index(split, emc2col::bawrsplit::t0, emc2col::bawr::t0)];
  const double* kap_ =
    cols[emc2col::select_index(split, emc2col::bawrsplit::kappa, emc2col::bawr::kappa)];
  const double* pw_ = cols[emc2col::select_index(split, emc2col::bawrsplit::p, emc2col::bawr::p)];
  const double* delta_ = split ? cols[emc2col::bawrsplit::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double log_pdf = bawr_log_pdf(
      tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i], kap_[i], pw_[i], launch, pd,
      BAWR_DENOM_FLOOR, delta_ ? delta_[i] : 0.0);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

void pbawr_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawr_launch_of(ctx);
  const bool split = launch == BAWR_LAUNCH_SPLITLOGNORMAL;
  const double* p1_ = cols[emc2col::select_index(split, emc2col::bawrsplit::mu, emc2col::bawr::v)];
  const double* p2_ =
    cols[emc2col::select_index(split, emc2col::bawrsplit::sigma, emc2col::bawr::sv)];
  const double* B_ = cols[emc2col::select_index(split, emc2col::bawrsplit::B, emc2col::bawr::B)];
  const double* A_ = cols[emc2col::select_index(split, emc2col::bawrsplit::A, emc2col::bawr::A)];
  const double* t0_ = cols[emc2col::select_index(split, emc2col::bawrsplit::t0, emc2col::bawr::t0)];
  const double* kap_ =
    cols[emc2col::select_index(split, emc2col::bawrsplit::kappa, emc2col::bawr::kappa)];
  const double* pw_ = cols[emc2col::select_index(split, emc2col::bawrsplit::p, emc2col::bawr::p)];
  const double* delta_ = split ? cols[emc2col::bawrsplit::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(p1_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    double cdf = 0.0;
    if (ba_natural_cdf_bawr(tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i],
                            kap_[i], pw_[i], launch, pd, BAWR_DENOM_FLOOR,
                            BA_ACCEPT_RAW, cdf,
                            delta_ ? delta_[i] : 0.0)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = bawr_log_surv(
        tt, A_[i], B_[i] + A_[i], p1_[i], p2_[i], kap_[i], pw_[i], launch,
        pd, BAWR_DENOM_FLOOR, delta_ ? delta_[i] : 0.0);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

void bawr_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int /*n_par*/,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const int launch = bawr_launch_of(ctx);
  const bool split = launch == BAWR_LAUNCH_SPLITLOGNORMAL;
  const double* p1_ = cols[emc2col::select_index(split, emc2col::bawrsplit::mu, emc2col::bawr::v)];
  const double* p2_ =
    cols[emc2col::select_index(split, emc2col::bawrsplit::sigma, emc2col::bawr::sv)];
  const double* B_ = cols[emc2col::select_index(split, emc2col::bawrsplit::B, emc2col::bawr::B)];
  const double* A_ = cols[emc2col::select_index(split, emc2col::bawrsplit::A, emc2col::bawr::A)];
  const double* t0_ = cols[emc2col::select_index(split, emc2col::bawrsplit::t0, emc2col::bawr::t0)];
  const double* kap_ =
    cols[emc2col::select_index(split, emc2col::bawrsplit::kappa, emc2col::bawr::kappa)];
  const double* pw_ = cols[emc2col::select_index(split, emc2col::bawrsplit::p, emc2col::bawr::p)];
  const double* delta_ = split ? cols[emc2col::bawrsplit::delta] : nullptr;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(p1_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      double cdf = 0.0;
      if (ba_natural_cdf_bawr(tt, A_[r], B_[r] + A_[r], p1_[r], p2_[r],
                              kap_[r], pw_[r], launch, pd, BAWR_DENOM_FLOOR,
                              BA_ACCEPT_RAW, cdf,
                              delta_ ? delta_[r] : 0.0)) {
        if (cdf > 0.0) logS += std::log1p(-cdf);
      } else {
        const double ls = bawr_log_surv(
          tt, A_[r], B_[r] + A_[r], p1_[r], p2_[r], kap_[r], pw_[r], launch,
          pd, BAWR_DENOM_FLOOR, delta_ ? delta_[r] : 0.0);
        if (!(ls > R_NegInf) || ISNAN(ls)) { bad = true; break; }
        logS += ls;
      }
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
  (void)n_rows_total;
}

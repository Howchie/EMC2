#include "model_BAwF.h"

double bawf_log_H(const BawfGeom& g, double x) {
  return g.rho_inf ? x : g.rho * bawd_pk_log_tau(g.rho, x);
}

double bawf_log_Hp(const BawfGeom& g, double x) {
  return g.rho_inf ? x : (g.rho - 1.0) * bawd_pk_log_tau(g.rho, x);
}

double bawf_log_Hpp(const BawfGeom& g, double x) {
  return g.rho_inf ? x
                   : std::log1p(-1.0 / g.rho) +
                     (g.rho - 2.0) * bawd_pk_log_tau(g.rho, x);
}

double bawf_zrel(const BawfGeom& g, double x) {
  return std::exp(bawf_log_Hp(g, x)) * (1.0 - x / g.y_0);
}

double bawf_critical_launch(const BawfGeom& g, double s) {
  return g.k * g.b * std::exp(bawf_log_Hp(g, s));
}

double bawf_s_of_logratio(const BawfGeom& g, double log_ratio) {
  if (!(log_ratio > 0.0)) return 0.0;
  return g.rho_inf ? log_ratio
                   : g.rho * std::expm1(log_ratio / (g.rho - 1.0));
}

double bawf_newton_x(const BawfGeom& g, double c) {
  if (!(c > 0.0)) return g.y_0;      // z = 0: tangency at the far endpoint
  if (!(c < 1.0)) return 0.0;        // z >= b: already at threshold
  double lo = 0.0, hi = g.y_0;
  const double h0 = g.rho_inf ? 1.0 : (g.rho - 1.0) / g.rho;  // H''(0)
  double x = std::sqrt(2.0 * (1.0 - c) / h0);
  if (!(x > lo) || !(x < hi)) x = 0.5 * (lo + hi);
  for (int it = 0; it < 100; ++it) {
    const double f = bawf_zrel(g, x) - c;
    if (f > 0.0) lo = x; else hi = x;
    const double deriv = -x * std::exp(bawf_log_Hpp(g, x));
    double xn = (deriv < 0.0 && emc2_isfinite(deriv))
      ? x - f / deriv : 0.5 * (lo + hi);
    if (!(xn > lo) || !(xn < hi) || !emc2_isfinite(xn)) xn = 0.5 * (lo + hi);
    const bool done = std::fabs(xn - x) <= 1e-14 * std::fmax(1.0, xn);
    x = xn;
    if (done) break;
  }
  return x;
}

BawfGeom bawf_geometry(double A, double b, double k, double rho) {
  BawfGeom g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(k)) return g;
  if (!(b > 0.0) || !(A >= 0.0) || !(b >= A) || !(k >= 0.0)) return g;
  if (ISNAN(rho)) return g;
  // rho == 0.0 is the compileAttributes-safe FFI sentinel for Inf, exactly as
  // in bawd_geometry().
  g.rho_inf = (rho == 0.0) || (rho > 0.0 && !R_FINITE(rho));
  if (!g.rho_inf) {
    // rho = 1 has no finite endpoint and is not a member of this family.
    if (!(rho > 1.0)) return g;
    g.rho = rho;
    g.y_0 = rho / (rho - 1.0);
  } else {
    g.rho = R_PosInf;
    g.y_0 = 1.0;
  }
  g.ok = true; g.b = b; g.A = A; g.k = k;
  g.k_zero = (k <= BAWF_K_EPS);
  if (g.k_zero) {
    // No fading: the trace is z + V u, nothing saturates, and the model is
    // exactly an LBA.  y_A is left at y_0 so the frozen interval is empty.
    g.y_A = g.y_0;
    g.T_max = R_PosInf;
    g.T_sat_A = R_PosInf;
    g.V_c0 = 0.0;
    return g;
  }
  g.y_A = (A > 0.0) ? bawf_newton_x(g, A / b) : g.y_0;
  g.T_max = g.y_0 / k;
  g.T_sat_A = g.y_A / k;
  g.V_c0 = bawf_critical_launch(g, g.y_0);
  return g;
}

BawfAtU bawf_at_u(const BawfGeom& g, double u) {
  BawfAtU s;
  if (!g.ok || !(u > 0.0) || ISNAN(u)) return s;
  s.ok = true;
  const bool inf_u = (u == R_PosInf);
  s.q = u;

  if (g.k_zero) {
    s.Z = g.A;
    s.c = 0.0;
    // At u = Inf the required launch is zero, so the LBA finishes with
    // probability P(V > 0) -- the survivor evaluated at w = 0.
    s.w_hi = inf_u ? 0.0 : g.b / u;
    s.w_lo = inf_u ? 0.0 : (g.b - s.Z) / u;
    s.s_lo = g.y_A;
    s.s_hi = g.y_0;
    return s;
  }

  const double x = g.k * u;
  if (inf_u || x >= g.y_0) {
    s.Z = 0.0;
    s.saturated = true;
  } else {
    const double zstar = g.b * bawf_zrel(g, x);
    s.saturated = !(zstar > 0.0);
    s.Z = std::fmin(std::fmax(zstar, 0.0), g.A);
  }
  s.partial = (s.Z < g.A);
  s.s_lo = g.y_A;
  s.s_hi = inf_u ? g.y_0 : std::fmin(x, g.y_0);

  if (s.saturated) {
    s.w_hi = g.V_c0;
    s.w_lo = g.V_c0;
    s.c = g.V_c0;
  } else {
    const double bH = g.b * std::exp(bawf_log_H(g, x));
    s.w_hi = bH / u;
    s.w_lo = (bH - s.Z) / u;
    s.c = bawf_critical_launch(g, x);
  }
  return s;
}

double bawf_log_frozen_quad(const BawfGeom& g, double s_lo, double s_hi,
                                   double p1, double p2, bool logn,
                                   bool posdrift) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  const double log_kb = std::log(g.k) + std::log(g.b);
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_jac = std::log(s) + bawf_log_Hpp(g, s);
    // P(V >= k b H'(s)); in both launches the standardized argument is affine
    // in log H'(s) or in H'(s), so no exp/log round trip is needed for the
    // lognormal case.
    const double log_surv = logn
      ? pnorm_log_direct((log_kb + bawf_log_Hp(g, s) - p1) / p2, false)
      : pnorm_log_direct((p1 - bawf_critical_launch(g, s)) / p2, true);
    return log_jac + log_surv;
  };
  (void)posdrift;  // the normalizer is applied by the caller, as in BAwD
  // Split where the launch density turns over, i.e. where V_c(s) passes the
  // lognormal median / the normal mean.
  const double mid = logn
    ? bawf_s_of_logratio(g, p1 - log_kb)
    : ((p1 > 0.0) ? bawf_s_of_logratio(g, std::log(p1) - log_kb) : s_lo);
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log(g.b);
}

double bawf_log_frozen_normal(const BawfGeom& g, double s_lo,
                                     double s_hi, double v, double sv) {
  return bawf_log_frozen_quad(g, s_lo, s_hi, v, sv, false, true);
}

double bawf_log_frozen_logn(const BawfGeom& g, double s_lo, double s_hi,
                                   double mu, double sigma) {
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  const double w_a = bawf_critical_launch(g, s_lo);
  const double w_b = bawf_critical_launch(g, s_hi);
  if (!(w_b > w_a) || !emc2_isfinite(w_b) || !(w_a > 0.0))
    return bawf_log_frozen_quad(g, s_lo, s_hi, mu, sigma, true, true);
  const double log_kb = std::log(g.k) + std::log(g.b);

  if (g.rho_inf) {
    // (1/k) int log(w / (k b)) Gbar(w) dw.
    const double ta = log_lognormal_logratio_stoploss(w_a, mu, sigma, log_kb);
    const double tb = log_lognormal_logratio_stoploss(w_b, mu, sigma, log_kb);
    if (ta - tb > BAWF_MIN_LOG_GAP) {
      const double td = log_diff_exp(ta, tb);
      if (emc2_isfinite(td)) return td - std::log(g.k);
    }
    return bawf_log_frozen_quad(g, s_lo, s_hi, mu, sigma, true, true);
  }

  // (rho/k) [ (k b)^{-p} int w^p Gbar dw - int Gbar dw ], p = 1/(rho - 1).
  // int w^p Gbar dw is Lambda_m with -(m + 1) = p, i.e. m = -rho/(rho - 1).
  const double p = 1.0 / (g.rho - 1.0);
  const double m = -g.rho / (g.rho - 1.0);
  const double lc_a = log_lognormal_stoploss(w_a, mu, sigma);
  const double lc_b = log_lognormal_stoploss(w_b, mu, sigma);
  if (lc_a - lc_b > BAWF_MIN_LOG_GAP) {
    const double L0 = log_diff_exp(lc_a, lc_b);
    const double lam_a = log_lognormal_power_stoploss(w_a, mu, sigma, m);
    const double lam_b = log_lognormal_power_stoploss(w_b, mu, sigma, m);
    const double L1 = -p * log_kb + log_diff_exp(lam_a, lam_b);
    if (emc2_isfinite(L1)) {
      const double d = L0 - L1;   // w >= k b throughout, so L1 >= L0
      if (d < 0.0 && -std::expm1(d) > 1e-6)
        return std::log(g.rho) - std::log(g.k) + L1 + log1m_exp(d);
      if (!(L0 > R_NegInf))
        return std::log(g.rho) - std::log(g.k) + L1;
    }
  }
  return bawf_log_frozen_quad(g, s_lo, s_hi, mu, sigma, true, true);
}

double log_bawf_cdf_normal(double u, const BawfGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);

  if (g.A <= BAWF_A_EPS) {
    const double out = pnorm_log_direct((v - s.w_hi) / sv, true) - log_denom;
    return std::fmin(out, 0.0);
  }

  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double cl = (v - s.w_hi) / sv;
    const double ch = (v - s.w_lo) / sv;
    if (!emc2_isfinite(cl) && !emc2_isfinite(ch)) return R_NegInf;
    const double span = ch - cl;
    if (span > BAWF_MIN_SPAN && emc2_isfinite(span)) {
      const double li = log_normal_phi_integral(cl, ch);
      log_live = ISNAN(li)
        ? std::log(s.Z) + pnorm_log_direct(0.5 * (cl + ch), true)
        : li + std::log(sv) + std::log(s.q);
    } else {
      log_live = std::log(s.Z) + pnorm_log_direct(0.5 * (cl + ch), true);
    }
  }

  const double log_frozen = s.partial
    ? bawf_log_frozen_normal(g, s.s_lo, s.s_hi, v, sv) : R_NegInf;

  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A) -
    log_denom;
  if (ISNAN(out)) return R_NegInf;
  return std::fmin(out, 0.0);
}

double log_bawf_cdf_logn(double u, const BawfGeom& g, double mu,
                                double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const auto log_gbar = [&](double w) -> double {
    if (!(w > 0.0)) return 0.0;
    if (!emc2_isfinite(w)) return R_NegInf;
    return pnorm_log_direct((mu - std::log(w)) / sigma, true);
  };

  if (g.A <= BAWF_A_EPS) return std::fmin(log_gbar(s.w_hi), 0.0);

  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double lc_lo = log_lognormal_stoploss(s.w_lo, mu, sigma);
    const double lc_hi = log_lognormal_stoploss(s.w_hi, mu, sigma);
    if (lc_lo - lc_hi > BAWF_MIN_LOG_GAP) {
      const double ld = log_diff_exp(lc_lo, lc_hi);
      log_live = ISNAN(ld) ? R_NegInf : std::log(s.q) + ld;
    }
    if (!(log_live > R_NegInf)) {
      log_live = std::log(s.Z) + log_gbar(0.5 * (s.w_lo + s.w_hi));
    }
  }

  const double log_frozen = s.partial
    ? bawf_log_frozen_logn(g, s.s_lo, s.s_hi, mu, sigma) : R_NegInf;

  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  if (ISNAN(out)) return R_NegInf;
  return std::fmin(out, 0.0);
}

double bawf_log_frozen_surv_quad(const BawfGeom& g, double s_lo,
                                        double s_hi, double p1, double p2,
                                        bool logn, bool posdrift) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  const double log_kb = std::log(g.k) + std::log(g.b);
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_jac = std::log(s) + bawf_log_Hpp(g, s) + std::log(g.b);
    const double log_cdf = logn
      ? pnorm_log_direct((log_kb + bawf_log_Hp(g, s) - p1) / p2, true)
      : (posdrift
           ? log_normal_cdf_positive_raw(bawf_critical_launch(g, s), p1, p2)
           : pnorm_log_direct((bawf_critical_launch(g, s) - p1) / p2, true));
    return log_jac + log_cdf;
  };
  const double mid = logn
    ? bawf_s_of_logratio(g, p1 - log_kb)
    : ((p1 > 0.0) ? bawf_s_of_logratio(g, std::log(p1) - log_kb) : s_lo);
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES);
}

double bawf_log_frozen_surv_normal(const BawfGeom& g, double s_lo,
                                          double s_hi, double v, double sv,
                                          bool posdrift) {
  return bawf_log_frozen_surv_quad(g, s_lo, s_hi, v, sv, false, posdrift);
}

double bawf_log_frozen_surv_logn(const BawfGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma) {
  return bawf_log_frozen_surv_quad(g, s_lo, s_hi, mu, sigma, true, false);
}

double log_bawf_surv_normal(double u, const BawfGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  if (u == R_PosInf && g.k_zero)
    return std::fmin(log_normal_cdf_positive(s.w_hi, v, sv, posdrift,
                                             denom_floor), 0.0);
  if (g.A <= BAWF_A_EPS)
    return std::fmin(log_normal_cdf_positive(s.w_hi, v, sv, posdrift,
                                             denom_floor), 0.0);
  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double lo = (s.w_lo - v) / sv;
    const double hi = (s.w_hi - v) / sv;
    const double li = log_normal_phi_integral_positive_raw(lo, hi, v, sv,
                                                            posdrift);
    if (li > R_NegInf) log_live = li + std::log(sv) + std::log(s.q);
  }
  const double log_frozen = s.partial
    ? bawf_log_frozen_surv_normal(g, s.s_lo, s.s_hi, v, sv, posdrift)
    : R_NegInf;
  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A) - log_denom;
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double log_bawf_surv_logn(double u, const BawfGeom& g, double mu,
                                  double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const auto log_g = [&](double w) {
    return (w > 0.0 && emc2_isfinite(w))
      ? pnorm_log_direct((std::log(w) - mu) / sigma, true) : R_NegInf;
  };
  if (u == R_PosInf && g.k_zero) return std::fmin(log_g(s.w_hi), 0.0);
  if (g.A <= BAWF_A_EPS) return std::fmin(log_g(s.w_hi), 0.0);
  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double pa = log_lognormal_put(s.w_hi, mu, sigma);
    const double pb = log_lognormal_put(s.w_lo, mu, sigma);
    if (pa - pb > BAWF_MIN_LOG_GAP)
      log_live = std::log(s.q) + log_diff_exp(pa, pb);
    if (!(log_live > R_NegInf))
      log_live = std::log(s.Z) + log_g(0.5 * (s.w_hi + s.w_lo));
  }
  const double log_frozen = s.partial
    ? bawf_log_frozen_surv_logn(g, s.s_lo, s.s_hi, mu, sigma) : R_NegInf;
  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double log_bawf_pdf_normal(double u, const BawfGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);

  if (g.A <= BAWF_A_EPS) {
    const double wgt = s.w_hi - s.c;
    if (!(wgt > 0.0)) return R_NegInf;
    return dnormP((v - s.w_hi) / sv, 0.0, 1.0, true) + std::log(wgt) -
      std::log(sv) - std::log(s.q) - log_denom;
  }
  const double cl = (v - s.w_hi) / sv, ch = (v - s.w_lo) / sv;
  const double span = ch - cl;
  const double log_scale = -std::log(g.A) - log_denom;
  if (span > BAWF_MIN_SPAN && emc2_isfinite(span)) {
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
          bracket.log_abs > max_term + BAWF_LOG_BRACKET_MIN) {
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

double log_bawf_pdf_logn(double u, const BawfGeom& g, double mu,
                                double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;

  if (g.A <= BAWF_A_EPS) {
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
          bracket.log_abs > max_term + BAWF_LOG_BRACKET_MIN) {
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

bool bawf_natural_cdf_normal(double u, const BawfGeom& g, double v,
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
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.q > 0.0)) { cdf = 0.0; return lenient; }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;

  if (g.A <= BAWF_A_EPS) {
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
      if (span > BAWF_MIN_SPAN && emc2_isfinite(span)) {
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
      const double log_frozen = bawf_log_frozen_normal(g, s.s_lo, s.s_hi, v, sv);
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

bool bawf_natural_cdf_logn(double u, const BawfGeom& g, double mu,
                                  double sigma, int accept_mode, double &cdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  const auto accept = [accept_mode](double &p) {
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(p);
    if (accept_mode == BA_ACCEPT_RAW) return p < 1.0 - 1e-8;
    if (p > 1.0) p = 1.0;
    return true;
  };

  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) { cdf = 0.0; return lenient; }
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.w_hi > 0.0) || !(s.w_lo > 0.0)) { cdf = 0.0; return lenient; }

  if (g.A <= BAWF_A_EPS) {
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
      if (diff > 0.0 && scale > 0.0 && diff / scale > BAWF_MIN_LOG_GAP) {
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
      const double log_frozen = bawf_log_frozen_logn(g, s.s_lo, s.s_hi, mu,
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

bool bawf_natural_pdf_normal(double u, const BawfGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &pdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) {
    pdf = 0.0;
    return lenient;
  }
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.q > 0.0)) { pdf = 0.0; return lenient; }
  if (s.saturated) { pdf = 0.0; return true; }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;

  if (g.A <= BAWF_A_EPS) {
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
    if (span > BAWF_MIN_SPAN && emc2_isfinite(span)) {
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

bool bawf_natural_pdf_logn(double u, const BawfGeom& g, double mu,
                                  double sigma, int accept_mode, double &pdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) {
    pdf = 0.0;
    return lenient;
  }
  const BawfAtU s = bawf_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) ||
      !(s.w_hi > 0.0) || !(s.w_lo > 0.0)) { pdf = 0.0; return lenient; }
  if (s.saturated) { pdf = 0.0; return true; }

  if (g.A <= BAWF_A_EPS) {
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

bool ba_natural_cdf_bawf(double u, double A, double b, double p1,
                                double p2, double k, int launch, bool posdrift,
                                double rho, double denom_floor,
                                int accept_mode, double &cdf) {
  const BawfGeom g = bawf_geometry(A, b, k, rho);
  if (launch == BAWF_LAUNCH_LOGNORMAL)
    return bawf_natural_cdf_logn(u, g, p1, p2, accept_mode, cdf);
  return bawf_natural_cdf_normal(u, g, p1, p2, posdrift, denom_floor,
                                 accept_mode, cdf);
}

bool ba_natural_pdf_bawf(double u, double A, double b, double p1,
                                double p2, double k, int launch, bool posdrift,
                                double rho, double denom_floor,
                                int accept_mode, double &pdf) {
  const BawfGeom g = bawf_geometry(A, b, k, rho);
  if (launch == BAWF_LAUNCH_LOGNORMAL)
    return bawf_natural_pdf_logn(u, g, p1, p2, accept_mode, pdf);
  return bawf_natural_pdf_normal(u, g, p1, p2, posdrift, denom_floor,
                                 accept_mode, pdf);
}

double bawf_log_cdf(double u, double A, double b, double p1, double p2,
                           double k, int launch, bool posdrift, double rho,
                           double denom_floor) {
  const BawfGeom g = bawf_geometry(A, b, k, rho);
  if (launch == BAWF_LAUNCH_LOGNORMAL) return log_bawf_cdf_logn(u, g, p1, p2);
  return log_bawf_cdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawf_log_surv(double u, double A, double b, double p1, double p2,
                            double k, int launch, bool posdrift, double rho,
                            double denom_floor) {
  const BawfGeom g = bawf_geometry(A, b, k, rho);
  if (launch == BAWF_LAUNCH_LOGNORMAL) return log_bawf_surv_logn(u, g, p1, p2);
  return log_bawf_surv_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawf_log_pdf(double u, double A, double b, double p1, double p2,
                           double k, int launch, bool posdrift, double rho,
                           double denom_floor) {
  const BawfGeom g = bawf_geometry(A, b, k, rho);
  if (launch == BAWF_LAUNCH_LOGNORMAL) return log_bawf_pdf_logn(u, g, p1, p2);
  return log_bawf_pdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawf_cdf_norm(double t, double A, double b, double p1, double p2,
                            double k, int launch, bool posdrift, bool log_out,
                            double rho,
                            double denom_floor) {
  double cdf;
  if (ba_natural_cdf_bawf(t, A, b, p1, p2, k, launch, posdrift, rho,
                          denom_floor, BA_ACCEPT_STRICT, cdf))
    return log_out ? std::log(cdf) : cdf;
  return return_from_log(
    bawf_log_cdf(t, A, b, p1, p2, k, launch, posdrift, rho, denom_floor),
    log_out);
}

double bawf_pdf_norm(double t, double A, double b, double p1, double p2,
                            double k, int launch, bool posdrift, bool log_out,
                            double rho,
                            double denom_floor) {
  double pdf;
  if (ba_natural_pdf_bawf(t, A, b, p1, p2, k, launch, posdrift, rho,
                          denom_floor, BA_ACCEPT_STRICT, pdf))
    return log_out ? std::log(pdf) : pdf;
  return return_from_log(
    bawf_log_pdf(t, A, b, p1, p2, k, launch, posdrift, rho, denom_floor),
    log_out);
}

double bawf_cdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, int launch,
                                      bool posdrift, double rho) {
  double cdf;
  if (ba_natural_cdf_bawf(t, A, b, p1, p2, k, launch, posdrift, rho,
                          BAWF_DENOM_FLOOR, BA_ACCEPT_CLAMP, cdf))
    return cdf;
  const double lp = bawf_log_cdf(t, A, b, p1, p2, k, launch, posdrift, rho);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

double bawf_pdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, int launch,
                                      bool posdrift, double rho) {
  double pdf;
  if (ba_natural_pdf_bawf(t, A, b, p1, p2, k, launch, posdrift, rho,
                          BAWF_DENOM_FLOOR, BA_ACCEPT_CLAMP, pdf))
    return pdf;
  const double lp = bawf_log_pdf(t, A, b, p1, p2, k, launch, posdrift, rho);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

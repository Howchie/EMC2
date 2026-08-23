#include "model_BAwD.h"
#include <cmath>

// --------------------------------------------------------------------------
// BAwD analytic core (kept out of the header so downstream family headers
// retain only the shared geometry, constants, and quadrature templates).
// --------------------------------------------------------------------------

double bawd_em1my(double y) {
  if (std::fabs(y) < 1e-4) {
    // (y^2/2) (1 + y/3 + y^2/12 + y^3/60)
    return 0.5 * y * y * (1.0 + (y / 3.0) * (1.0 + (y / 4.0) * (1.0 + y / 5.0)));
  }
  return std::expm1(y) - y;
}

double bawd_newton_y(double c) {
  if (!(c > 0.0)) return 0.0;
  if (!emc2_isfinite(c)) return R_PosInf;
  // y ~ sqrt(2c) for small c; y ~ log(c + 1 + log(c + 1)) for large c.
  double y = (c < 1.0) ? std::sqrt(2.0 * c)
                       : std::log(c + 1.0 + std::log(c + 1.0));
  for (int it = 0; it < 60; ++it) {
    const double h = bawd_em1my(y) - c;
    const double hp = std::expm1(y);          // h'(y) = e^y - 1 > 0
    if (!(hp > 0.0) || !emc2_isfinite(hp)) break;
    double yn = y - h / hp;
    if (!(yn > 0.0)) yn = 0.5 * y;            // convexity keeps us positive
    const bool done = std::fabs(yn - y) <= 1e-14 * std::fmax(1.0, yn);
    y = yn;
    if (done) break;
  }
  return y;
}

double bawd_log_psi_prime(const BawdGeom& g, double x) {
  if (!g.rho_inf)
    return bawd_pk_log_psi_prime(g.rho, g.gamma,
                                 bawd_pk_log_tau(g.rho, x));
  return std::log(g.omg) - g.gamma * x + std::log(std::expm1(x));
}

double bawd_log_h(const BawdGeom& g, double x) {
  return g.rho_inf ? -x : -g.rho * bawd_pk_log_tau(g.rho, x);
}

double bawd_log_wrel(const BawdGeom& g, double x) {
  return g.rho_inf ? g.omg * x
                   : bawd_pk_log_wrel(g.rho, g.gamma,
                                       bawd_pk_log_tau(g.rho, x));
}

double bawd_kq(const BawdGeom& g, double x) {
  return g.rho_inf
    ? -std::expm1(-x)
    : bawd_pk_kq(g.rho, bawd_pk_log_tau(g.rho, x));
}

double bawd_kr(const BawdGeom& g, double x) {
  if (g.rho_inf) {
    if (g.gamma_one) return -std::expm1(-x);
    if (g.gamma_zero) return x;
    return -std::expm1(-g.gamma * x) / g.gamma;
  }
  return bawd_pk_kr(g.rho, g.gamma, bawd_pk_log_tau(g.rho, x));
}

double bawd_psi(double s, const BawdGeom& g) {
  if (!g.rho_inf)
    return bawd_pk_psi(g.rho, g.gamma, bawd_pk_log_tau(g.rho, s));
  if (g.gamma_zero) return bawd_em1my(s);
  if (s < 1e-4) {
    const double x = s;
    const double gg = g.gamma;
    return 0.5 * g.omg * x * x *
      (1.0 + x * (1.0 - 2.0 * gg) / 3.0 +
       x * x * (3.0 * gg * gg - 3.0 * gg + 1.0) / 12.0 +
       x * x * x * (-4.0 * gg * gg * gg + 6.0 * gg * gg - 4.0 * gg + 1.0) / 60.0);
  }
  return std::expm1(g.omg * s) + (g.omg / g.gamma) * std::expm1(-g.gamma * s);
}

double bawd_newton_s(double c, const BawdGeom& g) {
  if (!g.rho_inf) return bawd_pk_newton_x(g.rho, g.gamma, c);
  if (g.gamma_zero) return bawd_newton_y(c);
  if (!(c > 0.0)) return 0.0;
  if (!emc2_isfinite(c)) return R_PosInf;
  double s = (c < 1.0) ? std::sqrt(2.0 * c / g.omg)
    : std::log(c + g.omg / g.gamma + 1.0) / g.omg;
  for (int it = 0; it < 60; ++it) {
    const double h = bawd_psi(s, g) - c;
    const double hp = g.omg * std::exp(-g.gamma * s) * std::expm1(s);
    if (!(hp > 0.0) || !emc2_isfinite(hp)) break;
    double sn = s - h / hp;
    if (!(sn > 0.0)) sn = 0.5 * s;
    const bool done = std::fabs(sn - s) <= 1e-14 * std::fmax(1.0, sn);
    s = sn;
    if (done) break;
  }
  return s;
}

bool bawd_shape_flags(BawdGeom& g, double gamma, double rho) {
  if (!emc2_isfinite(gamma) || gamma < 0.0 || gamma > 1.0) return false;
  if (ISNAN(rho)) return false;
  g.rho_inf = (rho == 0.0) || (rho > 0.0 && !R_FINITE(rho));
  if (!g.rho_inf) {
    if (!(rho >= 1.0)) return false;
    g.rho = rho;
    g.rho_one = std::fabs(rho - 1.0) <= 1e-12;
    g.m_shape = rho * (1.0 - gamma);
    g.kq_inf = g.rho_one ? R_PosInf : rho / (rho - 1.0);
  } else {
    g.rho = R_PosInf;
    g.m_shape = R_PosInf;
    g.kq_inf = 1.0;
  }
  g.gamma = gamma;
  g.gamma_zero = gamma <= BAWD_GAMMA_EPS;
  g.gamma_one = gamma >= 1.0 - BAWD_GAMMA_EPS;
  g.omg = 1.0 - gamma;
  if (g.gamma_one) {
    g.frozen_alpha = R_PosInf;
    g.frozen_m = R_PosInf;
  } else if (g.rho_inf) {
    g.frozen_alpha = 1.0 / g.omg;
    g.frozen_m = g.gamma / g.omg;
  } else {
    g.frozen_alpha = (rho - 1.0) / (rho * g.omg);
    g.frozen_m = g.frozen_alpha - 1.0;
  }
  return true;
}

double bawd_ell_from_tmax(double b, double k, double Tmax,
                                 double gamma, double rho ) {
  if (ISNAN(b) || ISNAN(k) || ISNAN(Tmax)) return R_NaN;
  if (!(b > 0.0) || !(k > BAWD_K_EPS) || !(Tmax > 0.0) ||
      !emc2_isfinite(Tmax)) return 0.0;
  BawdGeom g;
  if (!bawd_shape_flags(g, gamma, rho) || g.gamma_one) return 0.0;
  const double theta = bawd_psi(k * Tmax, g);
  if (!(theta > 0.0) || !emc2_isfinite(theta)) return 0.0;
  return k * b / theta;
}

bool bawd_uses_tmax(double gamma) {
  return gamma < 1.0 - BAWD_GAMMA_EPS;
}

BawdGeom bawd_geometry(double A, double b, double k, double ell,
                              double gamma, double rho ) {
  BawdGeom g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(k) ||
      !emc2_isfinite(ell))
    return g;
  if (!(b > 0.0) || !(A >= 0.0) || !(b >= A) || !(k >= 0.0) ||
      !(ell >= 0.0)) return g;
  if (!bawd_shape_flags(g, gamma, rho)) return g;
  g.ok = true; g.b = b; g.A = A; g.k = k; g.ell = ell;
  g.k_zero = (k <= BAWD_K_EPS);
  g.ell_zero = (ell <= BAWD_ELL_EPS);
  if (g.k_zero || g.ell_zero || g.gamma_one) return g;
  g.y_0 = bawd_newton_s(k * b / ell, g);
  g.y_A = (b > A) ? bawd_newton_s(k * (b - A) / ell, g) : 0.0;
  g.T_max = g.y_0 / k;
  g.T_sat_A = g.y_A / k;
  g.V_c0 = !g.rho_inf
    ? bawd_critical_launch(g, g.y_0)
    : (g.gamma_zero
       ? k * b + ell * (1.0 + g.y_0)
       : ell * (1.0 + bawd_psi(g.y_0, g) -
                (g.omg / g.gamma) * std::expm1(-g.gamma * g.y_0)));
  return g;
}

double bawd_critical_launch(const BawdGeom& g, double s) {
  if (!g.rho_inf)
    return g.ell * std::exp(bawd_log_wrel(g, s));
  if (g.gamma_zero)
    return g.ell * (bawd_em1my(s) + 1.0 + s);
  return g.ell * (1.0 + bawd_psi(s, g) -
                  (g.omg / g.gamma) * std::expm1(-g.gamma * s));
}

double bawd_c_gamma(const BawdGeom& g, double u, double q) {
  if (g.gamma_one) return q;
  if (g.gamma_zero) return u;
  return -std::expm1(-g.gamma * g.k * u) / (g.gamma * g.k);
}

BawdAtU bawd_at_u(const BawdGeom& g, double u) {
  BawdAtU s;
  if (!g.ok || !(u > 0.0) || ISNAN(u)) return s;
  s.ok = true;
  const bool inf_u = (u == R_PosInf);

  if (!g.rho_inf) {
    if (g.k_zero) {
      s.q = u;
      s.E = 1.0;
      s.log_E = 0.0;
    } else if (inf_u) {
      s.q = g.kq_inf / g.k;
      s.E = 0.0;
      s.log_E = R_NegInf;
    } else {
      const double L = bawd_pk_log_tau(g.rho, g.k * u);
      s.log_E = -g.rho * L;
      s.E = std::exp(s.log_E);
      s.q = bawd_pk_kq(g.rho, L) / g.k;
    }
  } else if (g.k_zero) {
    s.q = u;
    s.E = 1.0;
    s.log_E = 0.0;
  } else if (inf_u) {
    s.q = 1.0 / g.k;
    s.E = 0.0;
    s.log_E = R_NegInf;
  } else {
    double E, G;
    bawl_leak_factors(g.k * u, E, G);
    s.E = E;
    s.log_E = -g.k * u;
    s.q = G / g.k;
  }
  if (g.k_zero) {
    s.E_g = 1.0; s.log_E_g = 0.0;
    s.E_rel = 1.0; s.log_E_rel = 0.0;
  } else if (g.gamma_zero) {
    s.E_g = 1.0; s.log_E_g = 0.0;
    s.E_rel = s.E; s.log_E_rel = s.log_E;
  } else if (inf_u) {
    s.E_g = 0.0; s.log_E_g = R_NegInf;
    s.E_rel = g.gamma_one ? 1.0 : 0.0;
    s.log_E_rel = g.gamma_one ? 0.0 : R_NegInf;
  } else if (!g.rho_inf) {
    s.log_E_g = g.gamma * s.log_E;
    s.E_g = std::exp(s.log_E_g);
    s.log_E_rel = g.omg * s.log_E;
    s.E_rel = std::exp(s.log_E_rel);
  } else {
    s.log_E_g = -g.gamma * g.k * u;
    s.E_g = std::exp(s.log_E_g);
    s.log_E_rel = -g.omg * g.k * u;
    s.E_rel = std::exp(s.log_E_rel);
  }

  if (g.k_zero || g.ell_zero || g.gamma_one) {
    s.Z = g.A;
    s.partial = false;
    s.saturated = false;
  } else {
    const double ku = g.k * u;
    if (inf_u || ku >= g.y_0) {
      s.Z = 0.0;
      s.saturated = true;
    } else {
      const double zstar = g.b - (g.ell / g.k) * bawd_psi(ku, g);
      s.saturated = !(zstar > 0.0);
      s.Z = std::fmin(std::fmax(zstar, 0.0), g.A);
    }
    s.partial = (s.Z < g.A);
    s.s_lo = g.y_A;
    s.s_hi = inf_u ? g.y_0 : std::fmin(ku, g.y_0);
  }
  if (!g.rho_inf && g.gamma_one) {
    s.w_hi = g.b / s.q + g.ell;
    s.w_lo = (g.b - s.Z) / s.q + g.ell;
  } else if (!g.rho_inf && !g.k_zero && !s.saturated) {
    const double ellu = g.ell_zero ? 0.0 : g.ell * bawd_kr(g, g.k * u) / g.k;
    s.w_hi = (g.b + ellu) / s.q;
    s.w_lo = (g.b - s.Z + ellu) / s.q;
  } else if (s.saturated) {
    s.w_hi = g.V_c0;
    s.w_lo = g.V_c0;
  } else if (g.k_zero) {
    // (b - z)/u + ell, exact at u = Inf (where it is just ell).
    s.w_hi = g.b / u + g.ell;
    s.w_lo = (g.b - s.Z) / u + g.ell;
  } else {
    const double ellu = g.ell_zero ? 0.0 : g.ell * bawd_c_gamma(g, u, s.q);
    s.w_hi = (g.b + ellu) / s.q;
    s.w_lo = (g.b - s.Z + ellu) / s.q;
  }
  return s;
}

double bawd_log_frozen_normal(const BawdGeom& g, double s_lo, double s_hi,
                                     double v, double sv) {
  if (!g.rho_inf) {
    if (!(s_hi > s_lo) || !(sv > 0.0)) return R_NegInf;
    const auto lf = [&](double x) -> double {
      return bawd_log_psi_prime(g, x) +
        pnorm_log_direct((v - bawd_critical_launch(g, x)) / sv, true);
    };
    const double mid = (v > 0.0)
      ? bawd_pk_peak_x(g.rho, g.gamma, std::log(v / g.ell)) : s_lo;
    return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
      std::log(g.ell) - std::log(g.k);
  }
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    return std::log(e1) - g.gamma * s +
      pnorm_log_direct((v - bawd_critical_launch(g, s)) / sv, true);
  };
  // ell e^((1-gamma)s) is where the normal factor passes its midpoint.
  const double mid = (v > 0.0) ? std::log(v / g.ell) / g.omg : s_lo;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log1p(-g.gamma) + std::log(g.ell) - std::log(g.k);
}

double bawd_log_frozen_logn_quad(const BawdGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma) {
  if (!g.rho_inf) {
    if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
    const double log_ell = std::log(g.ell);
    const auto lf = [&](double x) -> double {
      return bawd_log_psi_prime(g, x) +
        pnorm_log_direct((log_ell + bawd_log_wrel(g, x) - mu) / sigma,
                         false);
    };
    const double mid = bawd_pk_peak_x(g.rho, g.gamma, mu - log_ell);
    return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
      std::log(g.ell) - std::log(g.k);
  }
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  const double x_ell = (std::log(g.ell) - mu) / sigma;
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    // P(V >= ell e^((1-gamma)s)) = Q(x_ell + (1-gamma)s/sigma).
    return std::log(e1) - g.gamma * s +
      pnorm_log_direct(x_ell + g.omg * s / sigma, false);
  };
  const double mid = -sigma * x_ell / g.omg;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log1p(-g.gamma) + std::log(g.ell) - std::log(g.k);
}

double bawd_log_frozen_logn(const BawdGeom& g, double s_lo, double s_hi,
                                   double mu, double sigma) {
  if (!g.rho_inf) {
    if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
    const double w_a = bawd_critical_launch(g, s_lo);
    const double w_b = bawd_critical_launch(g, s_hi);
    if (!(w_b > w_a) || !emc2_isfinite(w_b))
      return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);
    const double log_ell = std::log(g.ell);
    if (g.rho_one) {
      const double ta = log_lognormal_logratio_stoploss(
        w_a, mu, sigma, log_ell);
      const double tb = log_lognormal_logratio_stoploss(
        w_b, mu, sigma, log_ell);
      const double gap = ta - tb;
      if (gap > BAWD_MIN_LOG_GAP) {
        const double td = log_diff_exp(ta, tb);
        if (emc2_isfinite(td))
          return -std::log(g.k) - std::log1p(-g.gamma) + td;
      }
      return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);
    }
    const double lc_a = log_lognormal_stoploss(w_a, mu, sigma);
    const double lc_b = log_lognormal_stoploss(w_b, mu, sigma);
    if (lc_a - lc_b > BAWD_MIN_LOG_GAP) {
      const double L0 = log_diff_exp(lc_a, lc_b);
      const double lam_a = log_lognormal_power_stoploss(
        w_a, mu, sigma, g.frozen_m);
      const double lam_b = log_lognormal_power_stoploss(
        w_b, mu, sigma, g.frozen_m);
      const double L1 = g.frozen_alpha * log_ell +
        log_diff_exp(lam_a, lam_b);
      if (emc2_isfinite(L0)) {
        const double d = L1 - L0;
        if (d < 0.0 && -std::expm1(d) > 1e-6)
          return std::log(g.rho) - std::log(g.k) -
            std::log(g.rho - 1.0) + L0 + log1m_exp(d);
        if (!(L1 > R_NegInf))
          return std::log(g.rho) - std::log(g.k) -
            std::log(g.rho - 1.0) + L0;
      }
    }
    return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);
  }
  const double w_a = bawd_critical_launch(g, s_lo);
  const double w_b = bawd_critical_launch(g, s_hi);
  if (!(w_b > w_a) || !emc2_isfinite(w_b))
    return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);

  const double lc_a = log_lognormal_stoploss(w_a, mu, sigma);
  const double lc_b = log_lognormal_stoploss(w_b, mu, sigma);
  if (lc_a - lc_b > BAWD_MIN_LOG_GAP) {
    if (g.gamma_zero) {
      const double L0 = log_diff_exp(lc_a, lc_b);
      const double x_a = (std::log(w_a) - mu) / sigma;
      const double x_b = (std::log(w_b) - mu) / sigma;
      const double L1 = std::log(g.ell) + std::log(sigma) +
        log_normal_q_interval(x_a, x_b);
      if (!ISNAN(L0) && emc2_isfinite(L0)) {
        const double d = L1 - L0;
        // The two frozen primitives have lost relative precision.
        if (d < 0.0 && -std::expm1(d) > 1e-6)
          return L0 + log1m_exp(d) - std::log(g.k);
        if (!(L1 > R_NegInf))  // ell I_{-1} underflowed: I_0 stands alone
          return L0 - std::log(g.k);
      }
    } else {
      const double L0 = log_diff_exp(lc_a, lc_b);
      const double lam_a = log_lognormal_power_stoploss(
        w_a, mu, sigma, g.frozen_m);
      const double lam_b = log_lognormal_power_stoploss(
        w_b, mu, sigma, g.frozen_m);
      const double L1 = g.frozen_alpha * std::log(g.ell) +
        log_diff_exp(lam_a, lam_b);
      if (emc2_isfinite(L0)) {
        const double d = L1 - L0;
        if (d < 0.0 && -std::expm1(d) > 1e-6)
          return L0 + log1m_exp(d) - std::log(g.k);
        if (!(L1 > R_NegInf))
          return L0 - std::log(g.k);
      }
    }
  }
  return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);
}

double log_bawd_cdf_normal(double u, const BawdGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);

  if (g.A <= BAWD_A_EPS) {
    // Point start: one required launch strength, frozen at V_c0 past T_max.
    const double out = pnorm_log_direct((v - s.w_hi) / sv, true) - log_denom;
    return std::fmin(out, 0.0);
  }

  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double cl = (v - s.w_hi) / sv;   // c
    const double ch = (v - s.w_lo) / sv;   // c + m Z
    if (!emc2_isfinite(cl) && !emc2_isfinite(ch)) return R_NegInf;
    const double span = ch - cl;           // = m Z
    if (span > BAWD_MIN_SPAN && emc2_isfinite(span)) {
      const double li = log_normal_phi_integral(cl, ch);
      // 1/m = sv q, so (1/m) int Phi dz is the live part.
      log_live = ISNAN(li)
        ? std::log(s.Z) + pnorm_log_direct(0.5 * (cl + ch), true)
        : li + std::log(sv) + std::log(s.q);
    } else {
      log_live = std::log(s.Z) + pnorm_log_direct(0.5 * (cl + ch), true);
    }
  }

  const double log_frozen = s.partial
    ? bawd_log_frozen_normal(g, s.s_lo, s.s_hi, v, sv) : R_NegInf;

  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A) -
    log_denom;
  if (ISNAN(out)) return R_NegInf;
  return std::fmin(out, 0.0);
}

double log_bawd_cdf_logn(double u, const BawdGeom& g, double mu,
                                double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const auto log_gbar = [&](double w) -> double {
    if (!(w > 0.0)) return 0.0;
    if (!emc2_isfinite(w)) return R_NegInf;
    return pnorm_log_direct((mu - std::log(w)) / sigma, true);
  };

  if (g.A <= BAWD_A_EPS) return std::fmin(log_gbar(s.w_hi), 0.0);

  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    // int_0^Z Gbar(V*(u,z)) dz = q [C(w_lo) - C(w_hi)], C the stop-loss price.
    const double lc_lo = log_lognormal_stoploss(s.w_lo, mu, sigma);
    const double lc_hi = log_lognormal_stoploss(s.w_hi, mu, sigma);
    if (lc_lo - lc_hi > BAWD_MIN_LOG_GAP) {
      const double ld = log_diff_exp(lc_lo, lc_hi);
      log_live = ISNAN(ld) ? R_NegInf : std::log(s.q) + ld;
    }
    if (!(log_live > R_NegInf)) {
      // Interval too narrow (or cancelled) for an endpoint difference; the
      // midpoint limit is then accurate to the square of its width.
      log_live = std::log(s.Z) + log_gbar(0.5 * (s.w_lo + s.w_hi));
    }
  }

  const double log_frozen = s.partial
    ? bawd_log_frozen_logn(g, s.s_lo, s.s_hi, mu, sigma) : R_NegInf;

  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  if (ISNAN(out)) return R_NegInf;
  return std::fmin(out, 0.0);
}

double bawd_log_frozen_surv_normal(const BawdGeom& g, double s_lo,
                                          double s_hi, double v, double sv,
                                          bool posdrift, double denom_floor) {
  if (!(s_hi > s_lo) || !(sv > 0.0)) return R_NegInf;
  (void)denom_floor;
  const auto log_launch_cdf = [&](double w) -> double {
    return posdrift
      ? log_normal_cdf_positive_raw(w, v, sv)
      : pnorm_log_direct((w - v) / sv, true);
  };
  if (!g.rho_inf) {
    const auto lf = [&](double x) -> double {
      return bawd_log_psi_prime(g, x) +
        log_launch_cdf(bawd_critical_launch(g, x));
    };
    const double mid = (v > 0.0 && g.ell > 0.0)
      ? bawd_pk_peak_x(g.rho, g.gamma, std::log(v / g.ell)) : s_lo;
    return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
      std::log(g.ell) - std::log(g.k);
  }
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    return std::log(e1) - g.gamma * s +
      log_launch_cdf(bawd_critical_launch(g, s));
  };
  const double mid = (v > 0.0) ? std::log(v / g.ell) / g.omg : s_lo;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log1p(-g.gamma) + std::log(g.ell) - std::log(g.k);
}

double bawd_log_frozen_surv_logn(const BawdGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma) {
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  if (!g.rho_inf) {
    const double log_ell = std::log(g.ell);
    const auto lf = [&](double x) -> double {
      return bawd_log_psi_prime(g, x) +
        pnorm_log_direct((log_ell + bawd_log_wrel(g, x) - mu) / sigma, true);
    };
    const double mid = bawd_pk_peak_x(g.rho, g.gamma, mu - log_ell);
    return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
      std::log(g.ell) - std::log(g.k);
  }
  const double x_ell = (std::log(g.ell) - mu) / sigma;
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    return std::log(e1) - g.gamma * s +
      pnorm_log_direct(x_ell + g.omg * s / sigma, true);
  };
  const double mid = -sigma * x_ell / g.omg;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log1p(-g.gamma) + std::log(g.ell) - std::log(g.k);
}

double log_bawd_surv_normal(double u, const BawdGeom& g, double v,
                                   double sv, bool posdrift,
                                   double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  if (u == R_PosInf && g.k_zero)
    return std::fmin(log_normal_cdf_positive(s.w_hi, v, sv, posdrift,
                                             denom_floor), 0.0);
  if (g.A <= BAWD_A_EPS)
    return std::fmin(log_normal_cdf_positive(s.w_hi, v, sv, posdrift,
                                             denom_floor), 0.0);
  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double li = log_normal_phi_integral_positive_raw(
      (s.w_lo - v) / sv, (s.w_hi - v) / sv, v, sv, posdrift);
    if (li > R_NegInf) log_live = li + std::log(sv) + std::log(s.q);
  }
  const double log_frozen = (s.partial && !g.gamma_one)
    ? bawd_log_frozen_surv_normal(g, s.s_lo, s.s_hi, v, sv, posdrift,
                                  denom_floor) : R_NegInf;
  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  return ISNAN(out) ? R_NegInf : std::fmin(out - log_denom, 0.0);
}

double log_bawd_surv_logn(double u, const BawdGeom& g, double mu,
                                 double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok) return R_NegInf;
  const auto log_g = [&](double w) {
    return (w > 0.0 && emc2_isfinite(w))
      ? pnorm_log_direct((std::log(w) - mu) / sigma, true) : R_NegInf;
  };
  if (u == R_PosInf && g.k_zero) return std::fmin(log_g(s.w_hi), 0.0);
  if (g.A <= BAWD_A_EPS) return std::fmin(log_g(s.w_hi), 0.0);
  double log_live = R_NegInf;
  if (s.Z > 0.0) {
    const double pa = log_lognormal_put(s.w_hi, mu, sigma);
    const double pb = log_lognormal_put(s.w_lo, mu, sigma);
    if (pa - pb > BAWD_MIN_LOG_GAP)
      log_live = std::log(s.q) + log_diff_exp(pa, pb);
    if (!(log_live > R_NegInf))
      log_live = std::log(s.Z) + log_g(0.5 * (s.w_hi + s.w_lo));
  }
  const double log_frozen = (s.partial && !g.gamma_one)
    ? bawd_log_frozen_surv_logn(g, s.s_lo, s.s_hi, mu, sigma) : R_NegInf;
  const double out = log_sum_exp(log_live, log_frozen) - std::log(g.A);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double log_bawd_pdf_normal(double u, const BawdGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);
  if (g.A <= BAWD_A_EPS) {
    const double log_wgt = g.gamma_zero
      ? (g.ell_zero ? std::log(s.w_hi) + s.log_E
                     : ((s.w_hi * s.E - g.ell) > 0.0 ? std::log(s.w_hi * s.E - g.ell)
                                                     : R_NegInf))
      : (g.ell_zero ? std::log(s.w_hi) + s.log_E
                     : ((s.w_hi * s.E_rel - g.ell) > 0.0
                        ? s.log_E_g + std::log(s.w_hi * s.E_rel - g.ell) : R_NegInf));
    if (!(log_wgt > R_NegInf)) return R_NegInf;
    return dnormP((v - s.w_hi) / sv, 0.0, 1.0, true) + log_wgt -
      std::log(sv) - std::log(s.q) - log_denom;
  }
  const double cl = (v - s.w_hi) / sv, ch = (v - s.w_lo) / sv;
  const double span = ch - cl;
  const double log_scale = -std::log(g.A) - log_denom;
  if (span > BAWD_MIN_SPAN && emc2_isfinite(span)) {
    const double log_dphi = log_normal_interval(cl, ch);
    if (!ISNAN(log_dphi)) {
      const signed_log dphi = signed_log_sub(
        make_signed_log(dnormP(ch, 0.0, 1.0, true), 1),
        make_signed_log(dnormP(cl, 0.0, 1.0, true), 1));
      const signed_log v_term = signed_log_product(v, log_dphi);
      const signed_log sv_term = (dphi.sign != 0)
        ? make_signed_log(std::log(sv) + dphi.log_abs, dphi.sign)
        : make_signed_log(R_NegInf, 0);
      signed_log lead = signed_log_add(v_term, sv_term);
      if (lead.sign != 0) lead.log_abs += s.log_E_rel;
      const signed_log clear_term = signed_log_product(-g.ell, log_dphi);
      const signed_log bracket = signed_log_add(lead, clear_term);
      const double max_term = std::fmax(lead.log_abs, clear_term.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWD_LOG_BRACKET_MIN) {
        return bracket.log_abs + s.log_E_g + log_scale;
      }
      const double w_mid = 0.5 * (s.w_hi + s.w_lo);
      const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
      if (!(wgt > 0.0)) return R_NegInf;
      return std::log(wgt) + log_dphi + log_scale;
    }
  }
  const double w_mid = 0.5 * (s.w_hi + s.w_lo);
  const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
  if (!(wgt > 0.0)) return R_NegInf;
  return std::log(s.Z) + std::log(wgt) - std::log(sv) - std::log(s.q) +
    dnormP(0.5 * (cl + ch), 0.0, 1.0, true) + log_scale;
}

double log_bawd_pdf_logn(double u, const BawdGeom& g, double mu,
                                double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;

  if (g.A <= BAWD_A_EPS) {
    const double log_wgt = g.gamma_zero
      ? (g.ell_zero ? std::log(s.w_hi) + s.log_E
                     : ((s.w_hi * s.E - g.ell) > 0.0 ? std::log(s.w_hi * s.E - g.ell)
                                                     : R_NegInf))
      : (g.ell_zero ? std::log(s.w_hi) + s.log_E
                     : ((s.w_hi * s.E_rel - g.ell) > 0.0
                        ? s.log_E_g + std::log(s.w_hi * s.E_rel - g.ell) : R_NegInf));
    if (!(log_wgt > R_NegInf)) return R_NegInf;
    return dlnorm_std(s.w_hi, mu, sigma, true) + log_wgt - std::log(s.q);
  }

  // f A = E * int w g(w) dw - ell * int g(w) dw over [w_lo, w_hi], i.e. the
  // lognormal partial expectation minus its probability -- both closed form.
  const double d1_lo = (mu + sigma * sigma - std::log(s.w_lo)) / sigma;
  const double d1_hi = (mu + sigma * sigma - std::log(s.w_hi)) / sigma;
  const double d2_lo = (mu - std::log(s.w_lo)) / sigma;
  const double d2_hi = (mu - std::log(s.w_hi)) / sigma;
  const double log_prob = log_normal_interval(d2_hi, d2_lo);
  const double log_scale = -std::log(g.A);

  if (!ISNAN(log_prob) && s.log_E_rel > R_NegInf) {
    const double log_pexp = log_normal_interval(d1_hi, d1_lo);
    if (!ISNAN(log_pexp)) {
      const signed_log t1 = make_signed_log(
        s.log_E_rel + mu + 0.5 * sigma * sigma + log_pexp, 1);
      const signed_log t2 = signed_log_product(g.ell, log_prob);
      const signed_log bracket = signed_log_sub(t1, t2);
      const double max_term = std::fmax(t1.log_abs, t2.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWD_LOG_BRACKET_MIN) {
        return bracket.log_abs + s.log_E_g + log_scale;
      }
    }
  }
  const double w_mid = 0.5 * (s.w_hi + s.w_lo);
  const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
  if (!(wgt > 0.0)) return R_NegInf;
  if (!ISNAN(log_prob) && log_prob > R_NegInf)
    return std::log(wgt) + log_prob + log_scale;
  return std::log(s.Z) + std::log(wgt) + dlnorm_std(w_mid, mu, sigma, true) -
    std::log(s.q) + log_scale;
}

bool bawd_natural_cdf_normal(double u, const BawdGeom& g, double v,
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

  if (!g.ok || !(sv > 0.0) || !(u > 0.0)) {
    cdf = 0.0;
    return lenient;
  }
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) || !(s.q > 0.0)) {
    cdf = 0.0;
    return lenient;
  }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;

  if (g.A <= BAWD_A_EPS) {
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
      if (span > BAWD_MIN_SPAN && emc2_isfinite(span)) {
        if (!lenient && !natural_normal_interval_safe(cl, ch)) return false;
        double scale = 0.0;
        const double integral = normal_phi_integral_nat(cl, ch, scale);
        if (!(integral > 0.0)) {
          if (!lenient) return false;
          live_part = 0.0;
        } else {
          if (!lenient && integral <= BAWL_NATURAL_REL_TOL * std::fmax(1.0, scale)) return false;
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
      const double log_frozen = bawd_log_frozen_normal(g, s.s_lo, s.s_hi, v, sv);
      if (log_frozen > R_NegInf) {
        frozen_part = std::exp(log_frozen);
      }
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

bool bawd_natural_cdf_logn(double u, const BawdGeom& g, double mu,
                                  double sigma, int accept_mode, double &cdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  const auto accept = [accept_mode](double &p) {
    if (accept_mode == BA_ACCEPT_STRICT) return natural_cdf_safe(p);
    if (accept_mode == BA_ACCEPT_RAW) return p < 1.0 - 1e-8;
    if (p > 1.0) p = 1.0;
    return true;
  };

  if (!g.ok || !(sigma > 0.0) || !(u > 0.0)) {
    cdf = 0.0;
    return lenient;
  }
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) || !(s.w_hi > 0.0) || !(s.w_lo > 0.0)) {
    cdf = 0.0;
    return lenient;
  }

  if (g.A <= BAWD_A_EPS) {
    const double z_norm = (mu - std::log(s.w_hi)) / sigma;
    if (!emc2_isfinite(z_norm)) return false;
    if (!lenient && std::fabs(z_norm) > BAWL_NATURAL_Z_MAX) return false;
    cdf = pnorm_std(z_norm, true, false);
  } else {
    double live_part = 0.0;
    if (s.Z > 0.0) {
      double scale = 0.0;
      const double diff = lognormal_stoploss_interval_nat(s.w_lo, s.w_hi, mu, sigma, scale);
      if (diff > 0.0 && scale > 0.0 && diff / scale > BAWD_MIN_LOG_GAP) {
        if (!lenient && diff <= BAWL_NATURAL_REL_TOL * std::fmax(1.0, scale)) return false;
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
      const double log_frozen = bawd_log_frozen_logn(g, s.s_lo, s.s_hi, mu, sigma);
      if (log_frozen > R_NegInf) {
        frozen_part = std::exp(log_frozen);
      }
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

bool bawd_natural_pdf_normal(double u, const BawdGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &pdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) {
    pdf = 0.0;
    return lenient;
  }
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) || !(s.q > 0.0)) {
    pdf = 0.0;
    return lenient;
  }
  if (s.saturated) {
    pdf = 0.0;
    return true;
  }

  double denom;
  if (!natural_normalizer(v, sv, posdrift, denom, denom_floor)) return false;

  if (g.A <= BAWD_A_EPS) {
    const double wgt = g.ell_zero ? s.E_g * (s.w_hi * s.E_rel)
                                 : s.E_g * (s.w_hi * s.E_rel - g.ell);
    if (!(wgt > 0.0)) {
      pdf = 0.0;
      return lenient;
    }
    const double z_norm = (v - s.w_hi) / sv;
    if (!emc2_isfinite(z_norm)) return false;
    if (!lenient && std::fabs(z_norm) > BAWL_NATURAL_Z_MAX) return false;
    pdf = wgt * dnormP(z_norm) / (sv * s.q * denom);
  } else {
    const double cl = (v - s.w_hi) / sv;
    const double ch = (v - s.w_lo) / sv;
    const double span = ch - cl;
    if (span > BAWD_MIN_SPAN && emc2_isfinite(span)) {
      if (!lenient && !natural_normal_interval_safe(cl, ch)) return false;
      double scale = 0.0;
      const double dphi = normal_interval_nat(cl, ch, scale);
      if (dphi > 0.0 && scale > 0.0 && dphi / scale > 1e-6) {
        const double dnorm_diff = dnormP(ch) - dnormP(cl);
        const double term1 = s.E_g * (v * s.E_rel - g.ell) * dphi;
        const double term2 = s.E_g * sv * s.E_rel * dnorm_diff;
        const double bracket = term1 + term2;
        const double max_scale = std::fmax(std::fabs(term1), std::fabs(term2));
        if (bracket > 0.0 && bracket > 1e-6 * max_scale) {
          pdf = bracket / (g.A * denom);
        } else {
          const double w_mid = 0.5 * (s.w_hi + s.w_lo);
          const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
          if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
          pdf = wgt * dphi / (g.A * denom);
        }
      } else {
        const double w_mid = 0.5 * (s.w_hi + s.w_lo);
        const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
        if (!(wgt > 0.0)) { pdf = 0.0; return lenient; }
        const double mid_z = 0.5 * (cl + ch);
        if (!lenient && std::fabs(mid_z) > BAWL_NATURAL_Z_MAX) return false;
        pdf = s.Z * wgt * dnormP(mid_z) / (g.A * sv * s.q * denom);
      }
    } else {
      const double w_mid = 0.5 * (s.w_hi + s.w_lo);
      const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
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

bool bawd_natural_pdf_logn(double u, const BawdGeom& g, double mu,
                                 double sigma, int accept_mode, double &pdf) {
  const bool lenient = accept_mode != BA_ACCEPT_STRICT;
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) {
    pdf = 0.0;
    return lenient;
  }
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || !emc2_isfinite(s.w_hi) || !emc2_isfinite(s.w_lo) || !(s.w_hi > 0.0) || !(s.w_lo > 0.0)) {
    pdf = 0.0;
    return lenient;
  }
  if (s.saturated) {
    pdf = 0.0;
    return true;
  }

  if (g.A <= BAWD_A_EPS) {
    const double wgt = g.ell_zero ? s.E_g * (s.w_hi * s.E_rel)
                                 : s.E_g * (s.w_hi * s.E_rel - g.ell);
    if (!(wgt > 0.0)) {
      pdf = 0.0;
      return lenient;
    }
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
    const double t1 = s.E_g * s.E_rel * M * pexp;
    const double t2 = s.E_g * g.ell * prob;
    const double bracket = t1 - t2;

    if (bracket > 0.0 && prob > 0.0 && pexp > 0.0 &&
        (prob / scale_prob > 1e-6) && (pexp / scale_pexp > 1e-6)) {
      const double max_scale = std::fmax(std::fabs(t1), std::fabs(t2));
      if (!lenient && bracket <= BAWL_NATURAL_REL_TOL * std::fmax(1.0, max_scale))
        return false;
      pdf = bracket / g.A;
    } else {
      const double w_mid = 0.5 * (s.w_hi + s.w_lo);
      const double wgt = s.E_g * (w_mid * s.E_rel - g.ell);
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

bool ba_natural_cdf_bawd(double u, double A, double b, double p1, double p2,
                                double k, double ell, int launch, bool posdrift,
                                double gamma, double rho, double denom_floor,
                                int accept_mode, double &cdf) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  if (launch == BAWD_LAUNCH_LOGNORMAL)
    return bawd_natural_cdf_logn(u, g, p1, p2, accept_mode, cdf);
  return bawd_natural_cdf_normal(u, g, p1, p2, posdrift, denom_floor, accept_mode, cdf);
}

bool ba_natural_pdf_bawd(double u, double A, double b, double p1, double p2,
                                double k, double ell, int launch, bool posdrift,
                                double gamma, double rho, double denom_floor,
                                int accept_mode, double &pdf) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  if (launch == BAWD_LAUNCH_LOGNORMAL)
    return bawd_natural_pdf_logn(u, g, p1, p2, accept_mode, pdf);
  return bawd_natural_pdf_normal(u, g, p1, p2, posdrift, denom_floor, accept_mode, pdf);
}

double bawd_log_cdf(double u, double A, double b, double p1, double p2,
                           double k, double ell, int launch, bool posdrift,
                           double gamma, double rho,
                           double denom_floor ) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  if (launch == BAWD_LAUNCH_LOGNORMAL)
    return log_bawd_cdf_logn(u, g, p1, p2);
  return log_bawd_cdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawd_log_surv(double u, double A, double b, double p1, double p2,
                            double k, double ell, int launch, bool posdrift,
                            double gamma, double rho,
                            double denom_floor ) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  if (launch == BAWD_LAUNCH_LOGNORMAL) return log_bawd_surv_logn(u, g, p1, p2);
  return log_bawd_surv_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawd_log_pdf(double u, double A, double b, double p1, double p2,
                           double k, double ell, int launch, bool posdrift,
                           double gamma, double rho,
                           double denom_floor ) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  if (launch == BAWD_LAUNCH_LOGNORMAL)
    return log_bawd_pdf_logn(u, g, p1, p2);
  return log_bawd_pdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

double bawd_cdf_norm(double t, double A, double b, double p1, double p2,
                            double k, double ell, int launch, bool posdrift,
                            bool log_out, double gamma, double rho,
                            double denom_floor ) {
  double cdf;
  if (ba_natural_cdf_bawd(t, A, b, p1, p2, k, ell, launch, posdrift, gamma,
                           rho, denom_floor, BA_ACCEPT_STRICT, cdf))
    return log_out ? std::log(cdf) : cdf;
  return return_from_log(
    bawd_log_cdf(t, A, b, p1, p2, k, ell, launch, posdrift, gamma, rho,
                 denom_floor),
    log_out);
}

double bawd_pdf_norm(double t, double A, double b, double p1, double p2,
                            double k, double ell, int launch, bool posdrift,
                            bool log_out, double gamma, double rho,
                            double denom_floor ) {
  double pdf;
  if (ba_natural_pdf_bawd(t, A, b, p1, p2, k, ell, launch, posdrift, gamma,
                           rho, denom_floor, BA_ACCEPT_STRICT, pdf))
    return log_out ? std::log(pdf) : pdf;
  return return_from_log(
    bawd_log_pdf(t, A, b, p1, p2, k, ell, launch, posdrift, gamma, rho,
                 denom_floor),
    log_out);
}

double bawd_cdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, double ell,
                                      int launch, bool posdrift, double gamma,
                                      double rho) {
  double cdf;
  if (ba_natural_cdf_bawd(t, A, b, p1, p2, k, ell, launch, posdrift, gamma,
                           rho, BAWD_DENOM_FLOOR, BA_ACCEPT_CLAMP, cdf))
    return cdf;
  const double lp = bawd_log_cdf(t, A, b, p1, p2, k, ell, launch, posdrift,
                                 gamma, rho);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

double bawd_pdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, double ell,
                                      int launch, bool posdrift, double gamma,
                                      double rho) {
  double pdf;
  if (ba_natural_pdf_bawd(t, A, b, p1, p2, k, ell, launch, posdrift, gamma,
                           rho, BAWD_DENOM_FLOOR, BA_ACCEPT_CLAMP, pdf))
    return pdf;
  const double lp = bawd_log_pdf(t, A, b, p1, p2, k, ell, launch, posdrift,
                                 gamma, rho);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

// --------------------------------------------------------------------------
// BAwD direct Rcpp wrappers.
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dbawd(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    NumericVector ell, int launch = 1, bool posdrift = true,
                    bool log_out = false, double gamma = 0.0,
                    double rho = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), pick(ell, i), launch,
                           posdrift, log_out, gamma, rho);
  return out;
}

// [[Rcpp::export]]
NumericVector pbawd(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    NumericVector ell, int launch = 1, bool posdrift = true,
                    bool log_out = false, double gamma = 0.0,
                    double rho = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), pick(ell, i), launch,
                           posdrift, log_out, gamma, rho);
  return out;
}

// [[Rcpp::export]]
double dbawd_norm(double t, double A, double b, double p1, double p2, double k,
                  double ell, int launch = 1, bool posdrift = true,
                  bool log_out = false, double gamma = 0.0,
                  double rho = 0.0) {
  return bawd_pdf_norm(t, A, b, p1, p2, k, ell, launch, posdrift, log_out,
                       gamma, rho);
}

// [[Rcpp::export]]
double pbawd_norm(double t, double A, double b, double p1, double p2, double k,
                  double ell, int launch = 1, bool posdrift = true,
                  bool log_out = false, double gamma = 0.0,
                  double rho = 0.0) {
  return bawd_cdf_norm(t, A, b, p1, p2, k, ell, launch, posdrift, log_out,
                       gamma, rho);
}

// [[Rcpp::export]]
double bawd_tmax(double A, double b, double k, double ell, double gamma = 0.0,
                 double rho = 0.0) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  return g.ok ? g.T_max : NA_REAL;
}

// [[Rcpp::export]]
NumericVector bawd_tmax_vec(NumericVector A, NumericVector b, NumericVector k,
                            NumericVector ell, NumericVector gamma = 0.0,
                            NumericVector rho = 0.0) {
  const int n = A.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawdGeom g = bawd_geometry(A[i], pick(b, i), pick(k, i),
                                     pick(ell, i), pick(gamma, i), pick(rho, i));
    out[i] = g.ok ? g.T_max : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector bawd_ell_vec(NumericVector Tmax, NumericVector b, NumericVector k,
                           NumericVector gamma = 0.0, NumericVector rho = 0.0) {
  const int n = Tmax.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_ell_from_tmax(pick(b, i), pick(k, i), Tmax[i],
                                pick(gamma, i), pick(rho, i));
  return out;
}

// [[Rcpp::export]]
double lognormal_stoploss_log(double v, double mu, double sigma) {
  return log_lognormal_stoploss(v, mu, sigma);
}

// [[Rcpp::export]]
double lognormal_power_stoploss_log(double v, double mu, double sigma, double m) {
  return log_lognormal_power_stoploss(v, mu, sigma, m);
}

// [[Rcpp::export]]
double lognormal_logratio_stoploss_log(double v, double mu, double sigma,
                                       double log_ell) {
  return log_lognormal_logratio_stoploss(v, mu, sigma, log_ell);
}

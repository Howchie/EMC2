#include "model_BAwD.h"
#include "model_BTAwL.h"
#include "race_contract.h"
#include "utility_functions.h"
#include "col_registry.h"

// ---------------------------------------------------------------------------
// BTAwL analytic core (kept out of the header so downstream model headers
// retain only the geometry, constants, and caller-visible declarations).
// ---------------------------------------------------------------------------

double btawl_h(double t, double k, double tau) {
  if (!(tau > 0.0) || t <= 0.0) return 0.0;
  if (t == R_PosInf) return k <= BTAWL_K_EPS ? tau : 0.0;
  const double x = t / tau;
  if (k <= BTAWL_K_EPS)
    return tau * (-std::expm1(-x) - x * std::exp(-x));
  const double d = k - 1.0 / tau;
  if (d == 0.0)
    return 0.5 * tau * x * x * std::exp(-x);
  const double y = d * t;
  // H = exp(-k t) / (tau d^2) * ((y - 1)e^y + 1). The expm1 form avoids
  // cancellation around the coincident-timescale case k = 1/tau.
  double numer;
  if (std::fabs(y) < 1e-4) {
    const double y2 = y * y;
    numer = y2 * (0.5 + y * (1.0 / 3.0 + y * (1.0 / 8.0 + y / 30.0)));
  } else {
    // Use the origin Taylor form for |y| <= 0.1 so the transition to the
    // coincident limit is smooth down to 1e-16 relative.
    if (std::fabs(y) < 0.1) {
      double term = 0.5 * y * y;
      numer = term;
      for (int m = 3; m <= 12; ++m) {
        term *= y / static_cast<double>(m);
        numer += (m - 1.0) * term;
      }
    } else {
      numer = (y - 1.0) * std::expm1(y) + y;
    }
  }
  const double h = std::exp(-k * t) * numer / (tau * d * d);
  return (h > 0.0 && emc2_isfinite(h)) ? h : 0.0;
}

double btawl_g(double t, double tau) {
  if (!(tau > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  return (t / tau) * std::exp(-t / tau);
}

double btawl_hp(double t, double k, double tau) {
  if (!(tau > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  return btawl_g(t, tau) - k * btawl_h(t, k, tau);
}

double btawl_tmax(double k, double tau) {
  if (!(k > BTAWL_K_EPS) || !(tau > 0.0)) return R_PosInf;

  // H has one maximum above tau.  Solve H'(t) = g(t) - k H(t) directly;
  // a fixed dimensionless bracket is only valid at k * tau = 1.
  double lo = tau;
  double hi = tau + 1.0 / k;
  if (!(hi > lo) || !emc2_isfinite(hi)) hi = 2.0 * lo;
  for (int it = 0; it < 80 && btawl_hp(hi, k, tau) > 0.0; ++it)
    hi += hi - lo;
  if (!(btawl_hp(lo, k, tau) >= 0.0) ||
      !(btawl_hp(hi, k, tau) <= 0.0) || !emc2_isfinite(hi))
    return R_PosInf;

  double x = std::fmin(std::fmax(tau + 1.0 / k, lo), hi);
  for (int it = 0; it < 80; ++it) {
    const double f = btawl_hp(x, k, tau);
    if (f > 0.0) lo = x; else hi = x;
    const double gg = btawl_g(x, tau);
    const double gp = gg * (1.0 / x - 1.0 / tau);
    const double d = gp - k * f;
    double xn = (emc2_isfinite(d) && d < 0.0) ? x - f / d
                                               : 0.5 * (lo + hi);
    if (!(xn > lo) || !(xn < hi) || !emc2_isfinite(xn))
      xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) <= 2e-15 * std::fmax(1.0, std::fabs(x)))
      return xn;
    x = xn;
  }
  return 0.5 * (lo + hi);
}

double btawl_tangent_time(double z, const BtawlGeom& g) {
  if (!(z > 0.0)) return g.t_max;
  if (g.k_zero) return R_PosInf;
  if (!(z < g.b)) return 0.0;

  // The tangency time is the zero of d V*(t,z) / dt.  Solve its numerator
  // directly; the previous Gamma reparameterization was not equivalent to
  // this condition and selected the wrong frozen/live seam.
  auto numerator = [&](double t) {
    const double h = btawl_h(t, g.k, g.tau);
    const double hp = btawl_hp(t, g.k, g.tau);
    const double E = std::exp(-g.k * t);
    return z * g.k * E * h - (g.b - z * E) * hp;
  };
  double lo = g.tau, hi = g.t_max;
  const double flo = numerator(lo), fhi = numerator(hi);
  if (!(flo < 0.0)) return lo;
  if (!(fhi >= 0.0)) return hi;
  for (int it = 0; it < 80; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (numerator(mid) > 0.0) hi = mid; else lo = mid;
  }
  return 0.5 * (lo + hi);
}

double btawl_tau_from_ttrans(double k, double Ttrans) {
  if (!(Ttrans > 0.0) || !emc2_isfinite(Ttrans)) return R_NaN;
  if (k <= BTAWL_K_EPS) return R_NaN;

  // T_max(k, tau) is strictly increasing in tau.  Bisect in log(tau) so
  // weak-leak endpoints do not spend half their iterations near zero.
  double hi = Ttrans;
  if (!(btawl_tmax(k, hi) >= Ttrans)) hi = std::nextafter(Ttrans, R_PosInf);
  double lo = hi;
  for (int it = 0; it < 128; ++it) {
    const double next = lo * 0.5;
    if (!(next > 0.0) || next == lo) break;
    lo = next;
    if (btawl_tmax(k, lo) < Ttrans) break;
  }
  double log_lo = std::log(lo), log_hi = std::log(hi);
  for (int it = 0; it < 80; ++it) {
    const double log_mid = 0.5 * (log_lo + log_hi);
    const double mid = std::exp(log_mid);
    const double tm = btawl_tmax(k, mid);
    if (!(tm > 0.0) || tm < Ttrans) log_lo = log_mid;
    else log_hi = log_mid;
  }
  return std::exp(0.5 * (log_lo + log_hi));
}

BtawlGeom btawl_geometry(double A, double b, double k, double tau) {
  BtawlGeom g;
  g.A = A; g.b = b; g.k = k; g.tau = tau;
  g.k_zero = (k <= BTAWL_K_EPS);
  g.ok = (tau > 0.0 && b > 0.0 && b > A && A >= 0.0 && k >= 0.0 &&
          emc2_isfinite(A) && emc2_isfinite(b) && emc2_isfinite(k) &&
          emc2_isfinite(tau));
  if (!g.ok) return g;
  g.t_max = g.k_zero ? R_PosInf : btawl_tmax(k, tau);
  g.h_max = g.k_zero ? tau : btawl_h(g.t_max, k, tau);
  g.s_lo = (g.k_zero || A <= BTAWL_A_EPS) ? g.t_max : btawl_tangent_time(A, g);
  return g;
}

BtawlGeom btawl_geometry_from_ttrans(double A, double b, double k,
                                     double Ttrans,
                                     double tau_hint) {
  BtawlGeom g;
  g.A = A; g.b = b; g.k = k;
  g.k_zero = (k <= BTAWL_K_EPS);
  g.t_max = Ttrans;
  const double tau = (!ISNAN(tau_hint) && tau_hint > 0.0)
    ? tau_hint : btawl_tau_from_ttrans(k, Ttrans);
  g.tau = tau;
  g.ok = (tau > 0.0 && b > 0.0 && b > A && A >= 0.0 && k >= 0.0 &&
          Ttrans > 0.0 && emc2_isfinite(A) && emc2_isfinite(b) &&
          emc2_isfinite(k) && emc2_isfinite(Ttrans));
  if (!g.ok) return g;
  g.h_max = g.k_zero ? tau : btawl_h(g.t_max, k, tau);
  g.s_lo = (g.k_zero || A <= BTAWL_A_EPS) ? g.t_max : btawl_tangent_time(A, g);
  return g;
}

BtawlGeom btawl_geometry_cached(ContextForRaceModels* ctx, bool endpoint,
                                double A, double b, double k, double clear,
                                double tau) {
  if (ctx != nullptr) {
    if (!ctx->btawl_cache)
      ctx->btawl_cache = std::make_shared<btawl::SolveCache>();
    btawl::SolveCache& cache = *ctx->btawl_cache;
    for (const btawl::SolveCacheEntry& entry : cache.geometry) {
      if (entry.endpoint == endpoint && entry.A == A && entry.b == b &&
          entry.k == k && entry.clear == clear &&
          (endpoint || entry.tau == tau))
        return entry.geom;
    }
    const BtawlGeom geom = endpoint
      ? btawl_geometry_from_ttrans(A, b, k, clear, tau)
      : btawl_geometry(A, b, k, tau);
    if (cache.geometry.size() >= btawl::SolveCache::max_entries)
      cache.geometry.pop_back();
    cache.geometry.insert(cache.geometry.begin(), {endpoint, A, b, k, clear, tau, geom});
    return geom;
  }
  return endpoint ? btawl_geometry_from_ttrans(A, b, k, clear, tau)
                  : btawl_geometry(A, b, k, tau);
}

void btawl_cache_new_particle(ContextForRaceModels* ctx) {
  if (ctx != nullptr && ctx->btawl_cache)
    ctx->btawl_cache->new_particle();
}

double btawl_log_surv_cached(ContextForRaceModels* ctx, double t,
                             const BtawlGeom& g, double clear,
                             double p1, double p2, int launch,
                             bool posdrift, double delta) {
  const auto evaluate = [&]() {
    if (launch == BTAWL_LAUNCH_WEIBULL)
      return btawl_log_eval(true, t, g, p1, p2, BTAWL_LAUNCH_WEIBULL, false, 0.0);
    return (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (launch == BTAWL_LAUNCH_SPLITLOGNORMAL
           ? btawl_log_eval(true, t, g, p1, p2, BTAWL_LAUNCH_SPLITLOGNORMAL, false, delta)
           : log_btawl_surv_logn(t, g, p1, p2))
      : log_btawl_surv_normal(t, g, p1, p2, posdrift);
  };
  if (ctx == nullptr || !emc2_isfinite(g.t_max) || !(t >= g.t_max))
    return evaluate();
  if (!ctx->btawl_cache)
    ctx->btawl_cache = std::make_shared<btawl::SolveCache>();
  btawl::SolveCache& cache = *ctx->btawl_cache;
  for (const btawl::PlateauCacheEntry& entry : cache.plateau_survivors) {
    if (entry.A == g.A && entry.b == g.b && entry.k == g.k &&
        entry.clear == clear && entry.tau == g.tau && entry.p1 == p1 &&
        entry.p2 == p2 && entry.launch == launch &&
        entry.posdrift == posdrift &&
        (launch != BTAWL_LAUNCH_SPLITLOGNORMAL || entry.delta == delta))
      return entry.value;
  }
  const double value = evaluate();
  if (cache.plateau_survivors.size() >= btawl::SolveCache::max_entries)
    cache.plateau_survivors.pop_back();
  cache.plateau_survivors.insert(cache.plateau_survivors.begin(),
                                 {g.A, g.b, g.k, clear, g.tau, p1, p2,
                                  value, launch, posdrift, delta});
  return value;
}

double btawl_vstar(double t, double z, const BtawlGeom& g) {
  if (t <= 0.0) return R_PosInf;
  const double h = btawl_h(t, g.k, g.tau);
  if (!(h > 0.0)) return R_PosInf;
  const double E = g.k_zero ? 1.0 : std::exp(-g.k * t);
  const double num = g.b - z * E;
  return (num > 0.0) ? num / h : 0.0;
}

double btawl_vstar_prime(double t, double z, const BtawlGeom& g) {
  if (t <= 0.0) return R_NegInf;
  const double h = btawl_h(t, g.k, g.tau);
  const double hp = btawl_hp(t, g.k, g.tau);
  if (!(h > 0.0)) return R_NegInf;
  const double E = g.k_zero ? 1.0 : std::exp(-g.k * t);
  const double num = z * g.k * E * h - (g.b - z * E) * hp;
  return num / (h * h);
}

double btawl_z_t(double t, const BtawlGeom& g) {
  if (!(t > 0.0) || g.k_zero) return g.A;
  if (t >= g.t_max) return 0.0;
  const double h = btawl_h(t, g.k, g.tau);
  const double hp = btawl_hp(t, g.k, g.tau);
  const double E = std::exp(-g.k * t);
  const double num = g.b * hp;
  const double den = E * (g.k * h + hp);
  if (!(den > 0.0)) return g.A;
  const double z = num / den;
  return (z < 0.0) ? 0.0 : ((z > g.A) ? g.A : z);
}

double btawl_tangent_z_prime(double t, const BtawlGeom& g) {
  if (!(t > 0.0) || g.k_zero) return R_NegInf;
  const double E = std::exp(-g.k * t);
  const double gg = btawl_g(t, g.tau);
  const double hp = btawl_hp(t, g.k, g.tau);
  const double gp = gg * (1.0 / t - 1.0 / g.tau);
  const double hpp = gp - g.k * hp;
  const double q = E * gg;
  const double qp = E * (gp - g.k * gg);
  const double num = g.b * (hpp * q - hp * qp);
  const double den = q * q;
  return (den > 0.0 && num != 0.0) ? num / den : R_NegInf;
}

double btawl_normal_denom(double v, double sv, bool posdrift) {
  if (!posdrift) return 1.0;
  double d = pnorm_std(v / sv, true, false);
  return std::fmax(d, BTAWL_DENOM_FLOOR);
}

double btawl_surv(double w, double p1, double p2, int launch,
                  bool posdrift, double delta) {
  if (!(w > 0.0)) {
    if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL || launch == BTAWL_LAUNCH_WEIBULL) return 1.0;
    return pnorm_std(p1 / p2, true, false) /
      btawl_normal_denom(p1, p2, posdrift);
  }
  if (!emc2_isfinite(w)) return 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
      const double ls = log_split_lognormal_survivor(w, p1, p2, delta);
      return (ls > R_NegInf) ? std::exp(ls) : 0.0;
    }
    return pnorm_std((p1 - std::log(w)) / p2, true, false);
  }
  if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double ls = log_weibull_survivor(w, p1, p2);
    return (ls > R_NegInf) ? std::exp(ls) : 0.0;
  }
  const double d = btawl_normal_denom(p1, p2, posdrift);
  return pnorm_std((p1 - w) / p2, true, false) / d;
}

double btawl_pdf_v(double w, double p1, double p2, int launch,
                   bool posdrift, double delta) {
  if (!(w > 0.0) || !emc2_isfinite(w) || !(p2 > 0.0)) return 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
      const double lp = log_split_lognormal_density(w, p1, p2, delta);
      return (lp > R_NegInf) ? std::exp(lp) : 0.0;
    }
    return dlnorm_std(w, p1, p2, false);
  }
  if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double lp = log_weibull_density(w, p1, p2);
    return (lp > R_NegInf) ? std::exp(lp) : 0.0;
  }
  return dnormP(w, p1, p2, false) / btawl_normal_denom(p1, p2, posdrift);
}

double btawl_J(double x) {
  return x * pnorm_std(x, true, false) + dnormP(x);
}

double btawl_live_cdf(double t, double zlo, double zhi,
                      const BtawlGeom& g, double p1, double p2,
                      int launch, bool posdrift, double delta) {
  if (!(zhi > zlo)) return 0.0;
  const double h = btawl_h(t, g.k, g.tau);
  if (!(h > 0.0)) return 0.0;
  const double a = g.b / h;
  const double q = (g.k_zero ? 1.0 : std::exp(-g.k * t)) / h;
  if (!(q > 0.0) || !emc2_isfinite(q))
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift, delta) * (zhi - zlo);
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    const double whi = a - q * zhi;
    const double wlo = a - q * zlo;
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
      const double clo = log_split_lognormal_stoploss(whi, p1, p2, delta);
      const double chi = log_split_lognormal_stoploss(wlo, p1, p2, delta);
      if (clo > chi) {
        const double ld = log_diff_exp(clo, chi);
        const double val = std::exp(ld - std::log(q));
        if (val >= 0.0 && emc2_isfinite(val)) return val;
      }
      return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                        launch, posdrift, delta) * (zhi - zlo);
    }
    double sc1 = 0.0, sc2 = 0.0;
    const double clo = lognormal_stoploss_nat(whi, p1, p2, sc1);
    const double chi = lognormal_stoploss_nat(wlo, p1, p2, sc2);
    const double val = (clo - chi) / q;
    if (val >= 0.0 && emc2_isfinite(val)) return val;
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift, delta) * (zhi - zlo);
  }
  if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double whi = a - q * zhi;
    const double wlo = a - q * zlo;
    const double c_hi = log_weibull_stoploss(whi, p1, p2);
    const double c_lo = log_weibull_stoploss(wlo, p1, p2);
    if (c_lo > c_hi) {
      const double ld = log_diff_exp(c_lo, c_hi);
      const double val = std::exp(ld - std::log(q));
      if (val >= 0.0 && emc2_isfinite(val)) return val;
    }
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift, delta) * (zhi - zlo);
  }
  const double c = (p1 - a) / p2;
  const double m = q / p2;
  const double den = btawl_normal_denom(p1, p2, posdrift);
  if (!(m > 1e-14))
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift, delta) * (zhi - zlo);
  const double val = (btawl_J(c + m * zhi) - btawl_J(c + m * zlo)) / m / den;
  return (val >= 0.0 && emc2_isfinite(val)) ? val : 0.0;
}

double btawl_live_pdf(double t, double zlo, double zhi,
                      const BtawlGeom& g, double p1, double p2,
                      int launch, bool posdrift, double delta) {
  if (!(zhi > zlo)) return 0.0;
  const double h = btawl_h(t, g.k, g.tau);
  const double hp = btawl_hp(t, g.k, g.tau);
  if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
  const double E = std::exp(-g.k * t);
  const double q = E / h;
  if (!(q > 0.0)) return 0.0;
  const double w_hi = btawl_vstar(t, zlo, g);
  // On the partial live/frozen seam zhi = z_t(t), the tangency identity is
  // exact: V*(t, z_t(t)) = k b / g(t). Use that common expression so the
  // live density and frozen quadrature share the same boundary at full
  // precision. Clamped zhi=A is an all-live endpoint, not this seam.
  const double w_lo = (zhi > 0.0 && zhi < g.A && t < g.t_max)
    ? g.k * g.b / btawl_g(t, g.tau)
    : btawl_vstar(t, zhi, g);
  const double d0 = g.b * hp / (h * h);
  const double d1 = E * (g.k * h + hp) / (h * h);
  const double c0 = d0 - d1 * (g.b / h) / q;
  const double c1 = d1 / q;
  double mass = 0.0, first = 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
      const double ls_lo = log_split_lognormal_survivor(w_lo, p1, p2, delta);
      const double ls_hi = log_split_lognormal_survivor(w_hi, p1, p2, delta);
      const double lm_lo = log_split_lognormal_first_partial_moment(w_lo, p1, p2, delta);
      const double lm_hi = log_split_lognormal_first_partial_moment(w_hi, p1, p2, delta);
      const double log_mass = (ls_lo > ls_hi) ? log_diff_exp(ls_lo, ls_hi) : R_NegInf;
      const double log_first = (lm_lo > lm_hi) ? log_diff_exp(lm_lo, lm_hi) : R_NegInf;
      mass = (log_mass > R_NegInf) ? std::exp(log_mass) : 0.0;
      first = (log_first > R_NegInf) ? std::exp(log_first) : 0.0;
    } else {
      const double x0 = (std::log(w_lo) - p1) / p2;
      const double x1 = (std::log(w_hi) - p1) / p2;
      mass = pnorm_std(x1, true, false) - pnorm_std(x0, true, false);
      const double M = std::exp(p1 + 0.5 * p2 * p2);
      first = M * (pnorm_std(x1 - p2, true, false) -
                   pnorm_std(x0 - p2, true, false));
    }
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double lm = log_weibull_mass_interval(w_lo, w_hi, p1, p2);
    const double l1 = log_weibull_first_interval(w_lo, w_hi, p1, p2);
    mass = lm > R_NegInf ? std::exp(lm) : 0.0;
    first = l1 > R_NegInf ? std::exp(l1) : 0.0;
  } else {
    const double x0 = (w_lo - p1) / p2;
    const double x1 = (w_hi - p1) / p2;
    const double den = btawl_normal_denom(p1, p2, posdrift);
    mass = (pnorm_std(x1, true, false) - pnorm_std(x0, true, false)) / den;
    first = (p1 * (pnorm_std(x1, true, false) - pnorm_std(x0, true, false)) +
             p2 * (dnormP(x0) - dnormP(x1))) / den;
  }
  const double out = (c0 * mass + c1 * first) / (g.A * q);
  return (out > 0.0 && emc2_isfinite(out)) ? out : 0.0;
}

double btawl_frozen_cdf(double t, double zlo, double zhi,
                        const BtawlGeom& g, double p1, double p2,
                        int launch, bool posdrift, double delta) {
  if (!(zhi > zlo) || g.k_zero) return 0.0;
  const double li = btawl_log_frozen(false, t, g, p1, p2, launch, posdrift, delta);
  if (!(li > R_NegInf)) return 0.0;
  const double ln = (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    ? log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : 0.0;
  const double out = li - ln;
  return (emc2_isfinite(out) && out > -745.0) ? std::exp(out) : 0.0;
}

double btawl_frozen_split(const BtawlGeom& g, double s_lo, double s_hi,
                          double target) {
  if (!(target > 0.0) || !(s_hi > s_lo) || !(g.k > 0.0) || !(g.b > 0.0))
    return 0.5 * (s_lo + s_hi);
  const double desc_lo = std::fmax(s_lo, g.tau);
  const double desc_hi = std::fmin(s_hi, g.t_max);
  if (!(desc_hi > desc_lo)) return 0.5 * (s_lo + s_hi);
  const double target_g = g.k * g.b / target;
  if (!(target_g > 0.0) || !emc2_isfinite(target_g)) return desc_hi;
  const double g_lo = btawl_g(desc_lo, g.tau);
  const double g_hi = btawl_g(desc_hi, g.tau);
  if (target_g >= g_lo) return desc_lo;
  if (target_g <= g_hi) return desc_hi;
  double lo = desc_lo, hi = desc_hi;
  const double tol = 2e-14 * std::fmax(1.0, std::fabs(desc_hi));
  double x = 0.5 * (lo + hi);
  for (int i = 0; i < 64; ++i) {
    const double val = btawl_g(x, g.tau) - target_g;
    if (val > 0.0) lo = x; else hi = x;
    if (hi - lo <= tol) break;
    const double gp = (1.0 / g.tau - x / (g.tau * g.tau)) * std::exp(-x / g.tau);
    double xn = (std::fabs(gp) > 1e-14) ? x - val / gp : 0.5 * (lo + hi);
    if (!(xn > lo) || !(xn < hi) || !emc2_isfinite(xn))
      xn = 0.5 * (lo + hi);
    x = xn;
  }
  return x;
}

double btawl_log_frozen(bool survivor, double t, const BtawlGeom& g,
                        double p1, double p2,
                        int launch, bool posdrift, double delta) {
  if (g.k_zero || !(p2 > 0.0)) return R_NegInf;
  const double s_lo = g.s_lo;
  const double s_hi = std::fmin(std::fmax(t, s_lo), g.t_max);
  if (!(s_hi > s_lo)) return R_NegInf;
  
  split_lognormal_shape h;
  if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
    if (!split_lognormal_shape_params(p1, p2, delta, h)) return R_NegInf;
  }

  const auto lf = [&](double s) -> double {
    const double zp = btawl_tangent_z_prime(s, g);
    const double gg = btawl_g(s, g.tau);
    const double vc = (gg > 0.0) ? g.k * g.b / gg : R_PosInf;
    if (!(zp < 0.0) || !(vc > 0.0) || !emc2_isfinite(vc)) return R_NegInf;
    const double lp = (launch == BTAWL_LAUNCH_WEIBULL)
      ? (survivor ? log_weibull_cdf(vc, p1, p2)
                  : log_weibull_survivor(vc, p1, p2))
      : (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (launch == BTAWL_LAUNCH_SPLITLOGNORMAL
           ? (survivor
                ? (delta == 0.0 ? pnorm_log_direct((std::log(vc) - p1) / p2, true) : log_split_lognormal_cdf(vc, h))
                : (delta == 0.0 ? lnorm_log_surv_std(vc, p1, p2) : log_split_lognormal_survivor(vc, h)))
           : pnorm_log_direct((survivor ? std::log(vc) - p1 : p1 - std::log(vc)) / p2, true))
      : (survivor
           ? (posdrift ? log_normal_cdf_positive_raw(vc, p1, p2)
                       : pnorm_log_direct((vc - p1) / p2, true))
           : pnorm_log_direct((p1 - vc) / p2, true));
    return std::log(-zp) + lp;
  };

  double target = p1;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    target = std::exp(p1);
  } else if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    target = (delta == 0.0) ? std::exp(p1) : std::exp(h.c);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    target = p2;
  }
  
  const double split = btawl_frozen_split(g, s_lo, s_hi, target);
  return bawd_log_gl_split(lf, s_lo, s_hi, split, BAWD_GL_NODES);
}

double btawl_log_live(bool survivor, double t, double zlo, double zhi,
                      const BtawlGeom& g, double p1, double p2,
                      int launch, bool posdrift, double delta) {
  if (!(zhi > zlo) || !(p2 > 0.0)) return R_NegInf;
  const double h = btawl_h(t, g.k, g.tau);
  const double E = g.k_zero ? 1.0 : std::exp(-g.k * t);
  const double q = E / h;
  if (!(h > 0.0) || !(q > 0.0) || !emc2_isfinite(q)) return R_NegInf;
  const double a = g.b / h;
  const double w_hi = a - q * zlo;
  const double w_lo = a - q * zhi;
  if (!(w_hi >= w_lo) || !(w_lo > 0.0)) return R_NegInf;
  double li = R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    const double pa = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (survivor ? log_split_lognormal_put(w_hi, p1, p2, delta)
                  : log_split_lognormal_stoploss(w_lo, p1, p2, delta))
      : (survivor ? log_lognormal_put(w_hi, p1, p2)
                  : log_lognormal_stoploss(w_lo, p1, p2));
    const double pb = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (survivor ? log_split_lognormal_put(w_lo, p1, p2, delta)
                  : log_split_lognormal_stoploss(w_hi, p1, p2, delta))
      : (survivor ? log_lognormal_put(w_lo, p1, p2)
                  : log_lognormal_stoploss(w_hi, p1, p2));
    if (pa - pb > 1e-10) li = log_diff_exp(pa, pb) - std::log(q);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double pa = survivor ? log_weibull_put(w_hi, p1, p2)
                               : log_weibull_stoploss(w_lo, p1, p2);
    const double pb = survivor ? log_weibull_put(w_lo, p1, p2)
                               : log_weibull_stoploss(w_hi, p1, p2);
    if (pa - pb > 1e-10) li = log_diff_exp(pa, pb) - std::log(q);
  } else {
    const double lo = survivor ? (w_lo - p1) / p2 : (p1 - w_hi) / p2;
    const double hi = survivor ? (w_hi - p1) / p2 : (p1 - w_lo) / p2;
    const double x = survivor
      ? log_normal_phi_integral_positive_raw(lo, hi, p1, p2, posdrift)
      : log_normal_phi_integral(lo, hi);
    if (x > R_NegInf) li = x + std::log(p2) - std::log(q);
  }
  if (!(li > R_NegInf)) {
    const double wm = 0.5 * (w_hi + w_lo);
    const double lm = (launch == BTAWL_LAUNCH_WEIBULL)
      ? (survivor ? log_weibull_cdf(wm, p1, p2) : log_weibull_survivor(wm, p1, p2))
      : (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (launch == BTAWL_LAUNCH_SPLITLOGNORMAL
           ? (survivor ? log_split_lognormal_cdf(wm, p1, p2, delta)
                       : log_split_lognormal_survivor(wm, p1, p2, delta))
           : pnorm_log_direct((survivor ? std::log(wm) - p1 : p1 - std::log(wm)) / p2, true))
      : (survivor ? log_normal_cdf_positive_raw(wm, p1, p2)
                  : pnorm_log_direct((p1 - wm) / p2, true));
    li = std::log(zhi - zlo) + lm;
  }
  return li;
}

double btawl_log_eval(bool survivor, double t, const BtawlGeom& g,
                      double p1, double p2, int launch, bool posdrift,
                      double delta) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return R_NegInf;
  // Once the transient state has passed its unique maximum, the CDF is the
  // finishing mass.  Compute it as the complement of the stable survivor so
  // the two evaluators cannot disagree at the endpoint.
  if (!survivor && !g.k_zero && t >= g.t_max) {
    const double ls = btawl_log_eval(true, g.t_max, g, p1, p2, launch,
                                     posdrift, delta);
    if (ls == R_NegInf) return 0.0;
    if (!(ls < 0.0) || !emc2_isfinite(ls)) return R_NegInf;
    const double out = -std::expm1(ls);
    return (out > 0.0 && emc2_isfinite(out)) ? std::log(out) : R_NegInf;
  }
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    if (!(w > 0.0)) return R_NegInf;
    const double l = (launch == BTAWL_LAUNCH_WEIBULL)
      ? (survivor ? log_weibull_survivor(w, p1, p2)
                  : log_weibull_cdf(w, p1, p2))
      : (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (launch == BTAWL_LAUNCH_SPLITLOGNORMAL
           ? (survivor ? log_split_lognormal_cdf(w, p1, p2, delta)
                       : log_split_lognormal_survivor(w, p1, p2, delta))
           : pnorm_log_direct((survivor ? std::log(w) - p1 : p1 - std::log(w)) / p2, true))
      : (survivor ? log_normal_cdf_positive_raw(w, p1, p2)
                  : pnorm_log_direct((p1 - w) / p2, true));
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double zcut = g.k_zero ? g.A : std::fmin(std::fmax(btawl_z_t(u, g), 0.0), g.A);
  const double ll = btawl_log_live(survivor, u, 0.0, zcut, g, p1, p2,
                                   launch, posdrift, delta);
  const double lf = (!g.k_zero && zcut < g.A)
    ? btawl_log_frozen(survivor, u, g, p1, p2, launch, posdrift, delta)
    : R_NegInf;
  double out = log_sum_exp(ll, lf) - std::log(g.A);
  if (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    out -= log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double log_btawl_cdf_normal(double t, const BtawlGeom& g, double v,
                            double sv, bool posdrift) {
  return btawl_log_eval(false, t, g, v, sv, BTAWL_LAUNCH_NORMAL, posdrift);
}

double log_btawl_cdf_logn(double t, const BtawlGeom& g, double mu,
                          double sigma, double delta) {
  return btawl_log_eval(false, t, g, mu, sigma,
                        delta == 0.0 ? BTAWL_LAUNCH_LOGNORMAL : BTAWL_LAUNCH_SPLITLOGNORMAL,
                        false, delta);
}

double log_btawl_cdf_weib(double t, const BtawlGeom& g, double shape,
                          double scale) {
  return btawl_log_eval(false, t, g, shape, scale,
                        BTAWL_LAUNCH_WEIBULL, false, 0.0);
}

double log_btawl_surv_normal(double t, const BtawlGeom& g, double v,
                             double sv, bool posdrift) {
  return btawl_log_eval(true, t, g, v, sv, BTAWL_LAUNCH_NORMAL, posdrift);
}

double log_btawl_surv_logn(double t, const BtawlGeom& g, double mu,
                           double sigma, double delta) {
  return btawl_log_eval(true, t, g, mu, sigma,
                        delta == 0.0 ? BTAWL_LAUNCH_LOGNORMAL : BTAWL_LAUNCH_SPLITLOGNORMAL,
                        false, delta);
}

double log_btawl_surv_weib(double t, const BtawlGeom& g, double shape,
                           double scale) {
  return btawl_log_eval(true, t, g, shape, scale,
                        BTAWL_LAUNCH_WEIBULL, false, 0.0);
}

double btawl_cdf(double t, double A, double b, double p1, double p2,
                 double k, double tau, int launch, bool posdrift,
                 double delta) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return 0.0;
  if (!g.k_zero && t >= g.t_max) {
    const double lp = btawl_log_eval(false, t, g, p1, p2, launch,
                                     posdrift, delta);
    return (lp > R_NegInf && emc2_isfinite(lp)) ? std::exp(lp) : 0.0;
  }
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (u == R_PosInf) u = g.k_zero ? R_PosInf : g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    return std::fmin(std::fmax(btawl_surv(w, p1, p2, launch, posdrift, delta), 0.0), 1.0);
  }
  // Starts below zcut are still on the rising/live limb. Higher starts have
  // already reached their tangency wall and contribute the frozen mass.
  double zcut = g.k_zero ? g.A : btawl_z_t(u, g);
  zcut = std::fmin(std::fmax(zcut, 0.0), g.A);
  const double live = (u == R_PosInf && !g.k_zero) ? 0.0 :
    btawl_live_cdf(u, 0.0, zcut, g, p1, p2, launch, posdrift, delta);
  const double frozen = (!g.k_zero && zcut < g.A)
    ? btawl_frozen_cdf(u, zcut, g.A, g, p1, p2, launch, posdrift, delta) : 0.0;
  const double out = (live + frozen) / g.A;
  return std::fmin(std::fmax(out, 0.0), 1.0);
}

double btawl_pdf(double t, double A, double b, double p1, double p2,
                 double k, double tau, int launch, bool posdrift,
                 double delta) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  const double lp = btawl_log_pdf_from_geom(t, g, p1, p2, launch, posdrift, delta);
  return (lp > R_NegInf && lp < 709.0) ? std::exp(lp) : 0.0;
}

double btawl_cdf_from_geom(double t, const BtawlGeom& g, double p1,
                           double p2, int launch, bool posdrift,
                           double delta) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return 0.0;
  if (!g.k_zero && t >= g.t_max) {
    const double lp = btawl_log_eval(false, t, g, p1, p2, launch,
                                     posdrift, delta);
    return (lp > R_NegInf && emc2_isfinite(lp)) ? std::exp(lp) : 0.0;
  }
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (u == R_PosInf) u = g.k_zero ? R_PosInf : g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    return std::fmin(std::fmax(btawl_surv(w, p1, p2, launch, posdrift, delta), 0.0), 1.0);
  }
  double zcut = g.k_zero ? g.A : btawl_z_t(u, g);
  zcut = std::fmin(std::fmax(zcut, 0.0), g.A);
  const double live = (u == R_PosInf && !g.k_zero) ? 0.0 :
    btawl_live_cdf(u, 0.0, zcut, g, p1, p2, launch, posdrift, delta);
  const double frozen = (!g.k_zero && zcut < g.A)
    ? btawl_frozen_cdf(u, zcut, g.A, g, p1, p2, launch, posdrift, delta) : 0.0;
  const double out = (live + frozen) / g.A;
  return std::fmin(std::fmax(out, 0.0), 1.0);
}

bool btawl_natural_cdf_from_geom(double t, const BtawlGeom& g,
                                 double p1, double p2, int launch,
                                 bool posdrift, double& cdf,
                                 double delta) {
  if (launch == BTAWL_LAUNCH_WEIBULL ||
      (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0))
    return false;
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return false;
  if (!g.k_zero && t >= g.t_max) {
    const double lp = btawl_log_eval(false, t, g, p1, p2, launch,
                                     posdrift, delta);
    if (!(lp > R_NegInf) || !emc2_isfinite(lp)) return false;
    cdf = std::exp(lp);
    return emc2_isfinite(cdf);
  }
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (u == R_PosInf) u = g.k_zero ? R_PosInf : g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    cdf = std::fmin(std::fmax(btawl_surv(w, p1, p2, launch, posdrift, delta), 0.0), 1.0);
    return true;
  }
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(u, g), 0.0), g.A);
  const double live = (u == R_PosInf && !g.k_zero) ? 0.0 :
    btawl_live_cdf(u, 0.0, zcut, g, p1, p2, launch, posdrift, delta);
  const double frozen = (!g.k_zero && zcut < g.A)
    ? btawl_frozen_cdf(u, zcut, g.A, g, p1, p2, launch, posdrift, delta) : 0.0;
  const double out = (live + frozen) / g.A;
  cdf = std::fmin(std::fmax(out, 0.0), 1.0);
  return true;
}

double btawl_pdf_from_geom(double t, const BtawlGeom& g, double p1,
                           double p2, int launch, bool posdrift,
                           double delta) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  if (!g.k_zero && !(t < g.t_max)) return 0.0;
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double hp = btawl_hp(t, g.k, g.tau);
    if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
    const double w = g.b / h;
    return btawl_pdf_v(w, p1, p2, launch, posdrift, delta) * g.b * hp / (h * h);
  }
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  return btawl_live_pdf(t, 0.0, zcut, g, p1, p2, launch, posdrift, delta);
}

double btawl_log_launch_pdf(double w, double p1, double p2, int launch,
                            bool posdrift, double delta) {
  if (!(w > 0.0) || !emc2_isfinite(w) || !(p2 > 0.0)) return R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      return log_split_lognormal_density(w, p1, p2, delta);
    return dlnorm_std(w, p1, p2, true);
  }
  if (launch == BTAWL_LAUNCH_WEIBULL)
    return log_weibull_density(w, p1, p2);
  const double ld = dnormP((w - p1) / p2, 0.0, 1.0, true) - std::log(p2);
  return posdrift ? ld - log_positive_normalizer(p1, p2, true,
                                                 BTAWL_DENOM_FLOOR) : ld;
}

bool btawl_natural_pdf_accepted(double t, const BtawlGeom& g,
                                double p1, double p2, int launch,
                                bool posdrift, double p_nat,
                                double delta) {
  (void)posdrift;
  if (launch == BTAWL_LAUNCH_WEIBULL ||
      (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0))
    return false;
  if (!(p_nat > 0.0) || !emc2_isfinite(p_nat)) return false;
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double w = (h > 0.0) ? g.b / h : R_PosInf;
    const double x = (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (std::log(w) - p1) / p2 : (w - p1) / p2;
    return emc2_isfinite(x) && std::fabs(x) <= BAWL_NATURAL_Z_MAX;
  }
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  if (!(zcut > 0.0)) return false;
  const double h = btawl_h(t, g.k, g.tau);
  const double q = (g.k_zero ? 1.0 : std::exp(-g.k * t)) / h;
  if (!(h > 0.0) || !(q > 0.0) || !emc2_isfinite(q)) return false;
  const double a = g.b / h;
  const double w_hi = a;
  const double w_lo = (zcut < g.A && t < g.t_max)
    ? g.k * g.b / btawl_g(t, g.tau) : a - q * zcut;
  if (!(w_hi > w_lo) || !(w_lo > 0.0)) return false;
  const double x0 = (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
    ? (std::log(w_lo) - p1) / p2 : (w_lo - p1) / p2;
  const double x1 = (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
    ? (std::log(w_hi) - p1) / p2 : (w_hi - p1) / p2;
  if (!emc2_isfinite(x0) || !emc2_isfinite(x1) ||
      std::fabs(x0) > BAWL_NATURAL_Z_MAX ||
      std::fabs(x1) > BAWL_NATURAL_Z_MAX)
    return false;
  const double shift = (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) ? p2 : 0.0;
  if (std::fabs(x0 - shift) > BAWL_NATURAL_Z_MAX ||
      std::fabs(x1 - shift) > BAWL_NATURAL_Z_MAX)
    return false;
  const double cdf_hi = pnorm_std(x1, true, false);
  const double cdf_lo = pnorm_std(x0, true, false);
  const double mass = cdf_hi - cdf_lo;
  const bool mass_ok = mass > BAWL_NATURAL_REL_TOL *
    std::fmax(1.0, std::fabs(cdf_hi) + std::fabs(cdf_lo));
  if (!mass_ok) return false;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    const double M = std::exp(p1 + 0.5 * p2 * p2);
    const double m_hi = pnorm_std(x1 - p2, true, false);
    const double m_lo = pnorm_std(x0 - p2, true, false);
    const double first = M * (m_hi - m_lo);
    return first > BAWL_NATURAL_REL_TOL *
      std::fmax(1.0, std::fabs(M * m_hi) + std::fabs(M * m_lo));
  }
  const double lead = p1 * mass;
  const double corr = p2 * (dnormP(x0) - dnormP(x1));
  const double first = lead + corr;
  return first > BAWL_NATURAL_REL_TOL *
    std::fmax(1.0, std::fabs(lead) + std::fabs(corr));
}

namespace {

double btawl_log_live_pdf_stable(double t, const BtawlGeom& g,
                                 double p1, double p2, int launch,
                                 bool posdrift, double delta) {
  if (g.A <= BTAWL_A_EPS || !(p2 > 0.0)) return R_NegInf;
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  if (!(zcut > 0.0)) return R_NegInf;
  const double h = btawl_h(t, g.k, g.tau);
  const double E = g.k_zero ? 1.0 : std::exp(-g.k * t);
  const double q = E / h;
  if (!(h > 0.0) || !(q > 0.0) || !emc2_isfinite(q)) return R_NegInf;
  const double whi = g.b / h;
  const double wlo = (zcut < g.A && t < g.t_max)
    ? g.k * g.b / btawl_g(t, g.tau) : whi - q * zcut;
  if (!(whi > wlo) || !(wlo > 0.0) || !emc2_isfinite(wlo)) return R_NegInf;

  const double hp = btawl_hp(t, g.k, g.tau);
  const double d0 = g.b * hp / (h * h);
  const double d1 = E * (g.k * h + hp) / (h * h);
  const double c0 = d0 - d1 * (g.b / h) / q;
  const double c1 = d1 / q;
  if (!(c1 > 0.0) || !emc2_isfinite(c0) || !emc2_isfinite(c1))
    return R_NegInf;

  signed_log mass{R_NegInf, 0};
  signed_log first{R_NegInf, 0};

  if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
    const double ls_lo = log_split_lognormal_survivor(wlo, p1, p2, delta);
    const double ls_hi = log_split_lognormal_survivor(whi, p1, p2, delta);
    const double lm_lo = log_split_lognormal_first_partial_moment(wlo, p1, p2, delta);
    const double lm_hi = log_split_lognormal_first_partial_moment(whi, p1, p2, delta);
    if (!(ls_lo > R_NegInf) || !(lm_lo > R_NegInf)) return R_NegInf;
    mass = signed_log_sub(make_signed_log(ls_lo, 1), make_signed_log(ls_hi, 1));
    first = signed_log_sub(make_signed_log(lm_lo, 1), make_signed_log(lm_hi, 1));
    if (mass.sign <= 0 || !(mass.log_abs > R_NegInf)) return R_NegInf;
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double lm = log_weibull_mass_interval(wlo, whi, p1, p2);
    const double l1 = log_weibull_first_interval(wlo, whi, p1, p2);
    if (!(lm > R_NegInf) || !(l1 > R_NegInf)) return R_NegInf;
    mass = make_signed_log(lm, 1);
    first = make_signed_log(l1, 1);
  } else {
    const double x0 = (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (std::log(wlo) - p1) / p2 : (wlo - p1) / p2;
    const double x1 = (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (std::log(whi) - p1) / p2 : (whi - p1) / p2;
    if (!(x1 > x0) || !emc2_isfinite(x0) || !emc2_isfinite(x1))
      return R_NegInf;

    const double lm0 = log_normal_interval(x0, x1);
    if (!(lm0 > R_NegInf)) return R_NegInf;
    mass = make_signed_log(lm0, 1);
    if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
      const double lm1 = p1 + 0.5 * p2 * p2 +
        log_normal_interval(x0 - p2, x1 - p2);
      first = make_signed_log(lm1, 1);
    } else {
      const signed_log lead = signed_log_product(p1, lm0);
      const signed_log phi_diff = signed_log_sub(
        make_signed_log(dnormP(x0, 0.0, 1.0, true), 1),
        make_signed_log(dnormP(x1, 0.0, 1.0, true), 1));
      signed_log tail = signed_log_product(p2, phi_diff.log_abs);
      tail.sign = phi_diff.sign == 0 ? 0 : tail.sign * phi_diff.sign;
      first = signed_log_add(lead, tail);
    }
  }
  const signed_log affine = signed_log_add(
    signed_log_product(c0, mass.log_abs),
    signed_log_product(c1 * static_cast<double>(first.sign), first.log_abs));
  if (affine.sign <= 0 || !(affine.log_abs > R_NegInf)) return R_NegInf;
  const double log_den = (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    ? log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : 0.0;
  const double out = affine.log_abs - std::log(g.A) - std::log(q) - log_den;
  return emc2_isfinite(out) ? out : R_NegInf;
}

} // namespace

double btawl_log_pdf_from_geom(double t, const BtawlGeom& g, double p1,
                               double p2, int launch, bool posdrift,
                               double delta) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return R_NegInf;
  if (!g.k_zero && !(t < g.t_max)) return R_NegInf;
  const double p_nat = btawl_pdf_from_geom(t, g, p1, p2, launch, posdrift, delta);
  if (p_nat > 1e-280 &&
      btawl_natural_pdf_accepted(t, g, p1, p2, launch, posdrift, p_nat, delta))
    return std::log(p_nat);
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double hp = btawl_hp(t, g.k, g.tau);
    if (!(h > 0.0) || !(hp > 0.0)) return R_NegInf;
    const double w = g.b / h;
    return btawl_log_launch_pdf(w, p1, p2, launch, posdrift, delta) +
      std::log(g.b) + std::log(hp) - 2.0 * std::log(h);
  }
  const double stable = btawl_log_live_pdf_stable(t, g, p1, p2,
                                                  launch, posdrift, delta);
  if (stable > R_NegInf) return stable;
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  if (!(zcut > 0.0)) return R_NegInf;
  const auto lf = [&](double z) -> double {
    const double w = btawl_vstar(t, z, g);
    const double vp = btawl_vstar_prime(t, z, g);
    return (vp < 0.0) ? btawl_log_launch_pdf(w, p1, p2, launch, posdrift, delta) +
      std::log(-vp) : R_NegInf;
  };
  const double li = bawd_log_gl_split(lf, 0.0, zcut, 0.5 * zcut,
                                      BAWD_GL_NODES);
  return (li > R_NegInf) ? li - std::log(g.A) : R_NegInf;
}

double btawl_log_cdf(double t, double A, double b, double p1, double p2,
                     double k, double tau, int launch, bool posdrift,
                     double delta) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  if (launch == BTAWL_LAUNCH_WEIBULL)
    return log_btawl_cdf_weib(t, g, p1, p2);
  return (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
    ? log_btawl_cdf_logn(t, g, p1, p2, delta)
    : log_btawl_cdf_normal(t, g, p1, p2, posdrift);
}

double btawl_log_surv(double t, double A, double b, double p1, double p2,
                      double k, double tau, int launch, bool posdrift,
                      double delta) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  if (launch == BTAWL_LAUNCH_WEIBULL)
    return log_btawl_surv_weib(t, g, p1, p2);
  return (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
    ? log_btawl_surv_logn(t, g, p1, p2, delta)
    : log_btawl_surv_normal(t, g, p1, p2, posdrift);
}

double btawl_log_pdf(double t, double A, double b, double p1, double p2,
                     double k, double tau, int launch, bool posdrift,
                     double delta) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  return btawl_log_pdf_from_geom(t, g, p1, p2, launch, posdrift, delta);
}

double btawl_cdf_log(double t, double A, double b, double p1, double p2,
                     double k, double tau, int launch, bool posdrift,
                     double delta) {
  return btawl_log_cdf(t, A, b, p1, p2, k, tau, launch, posdrift, delta);
}

double btawl_pdf_log(double t, double A, double b, double p1, double p2,
                     double k, double tau, int launch, bool posdrift,
                     double delta) {
  return btawl_log_pdf(t, A, b, p1, p2, k, tau, launch, posdrift, delta);
}

double btawl_cdf_chart(double t, double A, double b, double p1,
                       double p2, double k, double clear, int launch,
                       bool posdrift, bool endpoint_chart,
                       double tau_hint, double delta) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  return btawl_cdf_from_geom(t, g, p1, p2, launch, posdrift, delta);
}

double btawl_pdf_chart(double t, double A, double b, double p1,
                       double p2, double k, double clear, int launch,
                       bool posdrift, bool endpoint_chart,
                       double tau_hint, double delta) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  const double lp = btawl_log_pdf_from_geom(t, g, p1, p2, launch, posdrift, delta);
  return (lp > R_NegInf && lp < 709.0) ? std::exp(lp) : 0.0;
}

double btawl_log_surv_chart(double t, double A, double b, double p1,
                            double p2, double k, double clear,
                            int launch, bool posdrift,
                            bool endpoint_chart,
                            double tau_hint, double delta) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  if (launch == BTAWL_LAUNCH_WEIBULL)
    return log_btawl_surv_weib(t, g, p1, p2);
  return (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
    ? log_btawl_surv_logn(t, g, p1, p2, delta)
    : log_btawl_surv_normal(t, g, p1, p2, posdrift);
}

double btawl_log_pdf_chart(double t, double A, double b, double p1,
                           double p2, double k, double clear,
                           int launch, bool posdrift,
                           bool endpoint_chart,
                           double tau_hint, double delta) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  return btawl_log_pdf_from_geom(t, g, p1, p2, launch, posdrift, delta);
}

double btawl_hs(double t, double k, double tau) {
  if (!(tau > 0.0) || !(t > 0.0)) return 0.0;
  if (k <= BTAWL_K_EPS) {
    return t + tau * std::expm1(-t / tau);
  }
  const double kt = k * tau;
  if (std::fabs(kt - 1.0) < 1e-5) {
    const double d = k - 1.0 / tau;
    const double y = d * t;
    const double Ek = std::exp(-k * t);
    return (1.0 - Ek) / k - t * Ek * (1.0 - 0.5 * y + y * y / 6.0);
  }
  const double inv_k = 1.0 / k;
  const double inv_diff = 1.0 / (kt - 1.0);
  return inv_k * (1.0 - std::exp(-k * t)) -
    (tau * inv_diff) * (std::exp(-t / tau) - std::exp(-k * t));
}

double btawl_hs_p(double t, double k, double tau) {
  if (!(tau > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  if (k <= BTAWL_K_EPS) return -std::expm1(-t / tau);

  // h_s' = (exp(-t/tau) - exp(-k t)) / (k tau - 1).  Evaluating this as
  // (1 - exp(-t/tau)) - k h_s loses all significant digits once the
  // sustained input has reached its plateau.  Use expm1 around the
  // coincident-timescale pole and the direct difference elsewhere.
  const double kt = k * tau;
  const double y = (kt - 1.0) * (t / tau);
  const double E = std::exp(-k * t);
  if (std::fabs(y) < 1e-4) {
    const double series = 1.0 + y * (0.5 + y * (1.0 / 6.0 + y / 24.0));
    return (E > 0.0 && emc2_isfinite(E * series)) ? E * series : 0.0;
  }
  const double D = std::exp(-t / tau);
  const double out = (D - E) / (kt - 1.0);
  return (out > 0.0 && emc2_isfinite(out)) ? out : 0.0;
}

namespace {

// The sustained start-point interval has width q*A.  Once that width is
// below machine resolution, evaluating endpoint differences is both slower
// and less accurate than taking the analytic point-limit at w = b/h.
constexpr double BTAWL_SUSTAINED_REL_WIDTH = 1e-10;
constexpr double BTAWL_SUSTAINED_WEIB_REL_WIDTH = 1e-8;

inline bool btawl_sustained_interval_collapsed(double a, double q, double A) {
  if (!(q > 0.0) || !emc2_isfinite(q) || !(A > 0.0)) return true;
  const double width = q * A;
  return !emc2_isfinite(width) || width <=
    BTAWL_SUSTAINED_REL_WIDTH * std::fmax(1.0, std::fabs(a));
}

inline bool btawl_sustained_weib_interval_collapsed(double a, double q,
                                                    double A) {
  if (!(q > 0.0) || !emc2_isfinite(q) || !(A > 0.0)) return true;
  const double width = q * A;
  return !emc2_isfinite(width) || width <=
    BTAWL_SUSTAINED_WEIB_REL_WIDTH * std::fmax(1.0, std::fabs(a));
}

inline double btawl_sustained_point_log_pdf(double a, double d0,
                                            double p1, double p2,
                                            int launch, bool posdrift,
                                            double delta) {
  if (!(a > 0.0) || !(d0 > 0.0) || !emc2_isfinite(a) ||
      !emc2_isfinite(d0) || !(p2 > 0.0)) return R_NegInf;
  double lp = R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    lp = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? log_split_lognormal_density(a, p1, p2, delta)
      : dlnorm_std(a, p1, p2, true);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    lp = log_weibull_density(a, p1, p2);
  } else {
    lp = dnormP(a, p1, p2, true) -
      ((posdrift) ? log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : 0.0);
  }
  return (lp > R_NegInf && emc2_isfinite(lp)) ? lp + std::log(d0) : R_NegInf;
}

inline double btawl_sustained_point_cdf(double a, double p1, double p2,
                                        int launch, bool posdrift,
                                        double delta) {
  const double out = btawl_surv(a, p1, p2, launch, posdrift, delta);
  return (out >= 0.0 && emc2_isfinite(out)) ? out : 0.0;
}

} // namespace

double btawl_sustained_cdf(double t, double A, double b, double p1, double p2,
                           double k, double tau_s, int launch, bool posdrift,
                           double delta) {
  if (!(p2 > 0.0)) return 0.0;
  // The Weibull launch uses the log survivor complement to avoid cancellation
  // near its upper plateau.  Keep legacy normal/lognormal arithmetic intact.
  if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double lsurv = btawl_sustained_log_surv(t, A, b, p1, p2, k, tau_s,
                                                  launch, posdrift, delta);
    if (lsurv < 0.0 && emc2_isfinite(lsurv))
      return std::fmin(std::fmax(-std::expm1(lsurv), 0.0), 1.0);
  }
  if (t == R_PosInf) {
    const double w = (k <= BTAWL_K_EPS) ? 0.0 : (k * b);
    return btawl_surv(w, p1, p2, launch, posdrift, delta);
  }
  const double h = btawl_hs(t, k, tau_s);
  if (!(h > 0.0)) return 0.0;
  const double E = std::exp(-k * t);
  const double q = E / h;
  const double a = b / h;
  if (!(a > 0.0) || !emc2_isfinite(a)) return 0.0;
  if (A <= BTAWL_A_EPS || btawl_sustained_interval_collapsed(a, q, A))
    return btawl_sustained_point_cdf(a, p1, p2, launch, posdrift, delta);

  const double whi = a;
  const double wlo = a - q * A;
  if (!(wlo > 0.0) || !(whi > wlo))
    return btawl_sustained_point_cdf(a, p1, p2, launch, posdrift, delta);

  double log_num = R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    const double clo = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? log_split_lognormal_stoploss(wlo, p1, p2, delta)
      : log_lognormal_stoploss(wlo, p1, p2);
    const double chi = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? log_split_lognormal_stoploss(whi, p1, p2, delta)
      : log_lognormal_stoploss(whi, p1, p2);
    if (clo > chi) log_num = log_diff_exp(clo, chi);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double clo = log_weibull_stoploss(wlo, p1, p2);
    const double chi = log_weibull_stoploss(whi, p1, p2);
    if (clo > chi) log_num = log_diff_exp(clo, chi);
  } else {
    const double c = (p1 - a) / p2;
    const double m = q / p2;
    if (m > 1e-14 && emc2_isfinite(c) && emc2_isfinite(m)) {
      const double val = btawl_J(c + m * A) - btawl_J(c);
      if (val > 0.0 && emc2_isfinite(val))
        return std::fmin(std::fmax(val / (m * A) /
                                    btawl_normal_denom(p1, p2, posdrift), 0.0), 1.0);
    }
  }

  if (log_num > R_NegInf) {
    const double out = std::exp(log_num - std::log(q) - std::log(A));
    if (out >= 0.0 && emc2_isfinite(out))
      return std::fmin(out, 1.0);
  }
  return btawl_sustained_point_cdf(a, p1, p2, launch, posdrift, delta);
}

double btawl_sustained_log_pdf(double t, double A, double b, double p1,
                               double p2, double k, double tau_s,
                               int launch, bool posdrift, double delta) {
  if (!(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return R_NegInf;
  const double h = btawl_hs(t, k, tau_s);
  const double hp = btawl_hs_p(t, k, tau_s);
  if (!(h > 0.0) || !(hp > 0.0)) return R_NegInf;

  const double a = b / h;
  const double d0 = b * hp / (h * h);
  if (!(a > 0.0) || !(d0 > 0.0) || !emc2_isfinite(a) || !emc2_isfinite(d0))
    return R_NegInf;
  if (A <= BTAWL_A_EPS)
    return btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);

  const double E = std::exp(-k * t);
  const double q = E / h;
  if (btawl_sustained_interval_collapsed(a, q, A))
    return btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);

  const double w_hi = a;
  const double w_lo = a - q * A;
  if (!(w_lo > 0.0) || !(w_hi > w_lo))
    return btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);

  const double d1 = E * (k * h + hp) / (h * h);
  const double c0 = d0 - d1 * a / q;
  const double c1 = d1 / q;
  if (!emc2_isfinite(c0) || !emc2_isfinite(c1))
    return btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);

  double log_mass = R_NegInf;
  double log_first = R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    const double x0 = (std::log(w_lo) - p1) / p2;
    const double x1 = (std::log(w_hi) - p1) / p2;
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL && delta != 0.0) {
      const double f_lo = log_split_lognormal_cdf(w_lo, p1, p2, delta);
      const double f_hi = log_split_lognormal_cdf(w_hi, p1, p2, delta);
      const double m_lo = log_split_lognormal_first_partial_moment(w_lo, p1, p2, delta);
      const double m_hi = log_split_lognormal_first_partial_moment(w_hi, p1, p2, delta);
      if (f_hi > f_lo) log_mass = log_diff_exp(f_hi, f_lo);
      if (m_lo > m_hi) log_first = log_diff_exp(m_lo, m_hi);
    } else {
      log_mass = log_normal_interval(x0, x1);
      log_first = p1 + 0.5 * p2 * p2 +
        log_normal_interval(x0 - p2, x1 - p2);
    }
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    log_mass = log_weibull_mass_interval(w_lo, w_hi, p1, p2);
    log_first = log_weibull_first_interval(w_lo, w_hi, p1, p2);
  } else {
    const double x0 = (w_lo - p1) / p2;
    const double x1 = (w_hi - p1) / p2;
    const double interval = log_normal_interval(x0, x1);
    if (interval > R_NegInf) {
      const signed_log lead = signed_log_product(p1, interval);
      const signed_log phi_diff = signed_log_sub(
        make_signed_log(log_phi_std(x0), 1),
        make_signed_log(log_phi_std(x1), 1));
      const signed_log tail = signed_log_product(p2, phi_diff.log_abs);
      const signed_log first = signed_log_add(
        lead, make_signed_log(tail.log_abs, tail.sign * phi_diff.sign));
      if (first.sign > 0) log_first = first.log_abs;
      log_mass = interval;
      if (posdrift) {
        const double log_den = log_positive_normalizer(p1, p2, true,
                                                       BTAWL_DENOM_FLOOR);
        log_mass -= log_den;
        log_first -= log_den;
      }
    }
  }

  if (!(log_mass > R_NegInf) || !(log_first > R_NegInf))
    return btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);

  const signed_log affine = signed_log_add(
    signed_log_product(c0, log_mass), signed_log_product(c1, log_first));
  if (affine.sign <= 0 || !(affine.log_abs > R_NegInf))
    return btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);
  const double out = affine.log_abs - std::log(A) - std::log(q);
  return emc2_isfinite(out) ? out :
    btawl_sustained_point_log_pdf(a, d0, p1, p2, launch, posdrift, delta);
}

double btawl_sustained_pdf(double t, double A, double b, double p1, double p2,
                           double k, double tau_s, int launch, bool posdrift,
                           double delta) {
  const double lp = btawl_sustained_log_pdf(t, A, b, p1, p2, k, tau_s,
                                            launch, posdrift, delta);
  return (lp > R_NegInf && lp < 709.0) ? std::exp(lp) : 0.0;
}

double btawl_local_race_cdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift, double delta) {
  if (pi <= 1e-14) return btawl_cdf(t, A, b, p1, p2, k, tau_t, launch, posdrift, delta);
  if (pi >= 1.0 - 1e-14) return btawl_sustained_cdf(t, A, b, p1, p2, k, tau_s, launch, posdrift, delta);

  double p1_S = p1, p2_S = p2;
  double p1_T = p1, p2_T = p2;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    p1_S += std::log(pi);
    p1_T += std::log(1.0 - pi);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    p1_S = p1_T = p1;
    p2_S = p2 * pi;
    p2_T = p2 * (1.0 - pi);
  } else {
    p1_S *= pi;         p2_S *= pi;
    p1_T *= (1.0 - pi); p2_T *= (1.0 - pi);
  }

  const double F_T = btawl_cdf(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift, delta);
  const double F_S = btawl_sustained_cdf(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift, delta);
  return 1.0 - (1.0 - F_T) * (1.0 - F_S);
}

double btawl_local_race_pdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift, double delta) {
  if (pi <= 1e-14) return btawl_pdf(t, A, b, p1, p2, k, tau_t, launch, posdrift, delta);
  if (pi >= 1.0 - 1e-14)
    return btawl_sustained_pdf(t, A, b, p1, p2, k, tau_s, launch, posdrift, delta);

  double p1_S = p1, p2_S = p2;
  double p1_T = p1, p2_T = p2;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    p1_S += std::log(pi);
    p1_T += std::log(1.0 - pi);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    p1_S = p1_T = p1;
    p2_S = p2 * pi;
    p2_T = p2 * (1.0 - pi);
  } else {
    p1_S *= pi;         p2_S *= pi;
    p1_T *= (1.0 - pi); p2_T *= (1.0 - pi);
  }

  const double F_T = btawl_cdf(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift, delta);
  const double F_S = btawl_sustained_cdf(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift, delta);
  const double f_T = btawl_pdf(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift, delta);
  const double f_S = btawl_sustained_pdf(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift, delta);

  return f_T * (1.0 - F_S) + f_S * (1.0 - F_T);
}

double btawl_local_race_cdf_log(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, double tau_t, double pi,
                                int launch, bool posdrift, double delta) {
  const double p = btawl_local_race_cdf(t, A, b, p1, p2, k, tau_s, tau_t, pi, launch, posdrift, delta);
  return p > 0.0 ? std::log(p) : R_NegInf;
}

double btawl_local_race_pdf_log(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, double tau_t, double pi,
                                int launch, bool posdrift, double delta) {
  if (pi >= 1.0 - 1e-14)
    return btawl_sustained_log_pdf(t, A, b, p1, p2, k, tau_s,
                                   launch, posdrift, delta);
  const double p = btawl_local_race_pdf(t, A, b, p1, p2, k, tau_s, tau_t, pi, launch, posdrift, delta);
  return p > 0.0 ? std::log(p) : R_NegInf;
}

double btawl_sustained_log_surv(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, int launch, bool posdrift,
                                double delta) {
  if (!emc2_isfinite(b) || !(p2 > 0.0) || !(t > 0.0)) return R_NegInf;
  if (t == R_PosInf) {
    const double w = (k <= BTAWL_K_EPS) ? 0.0 : (k * b);
    if (!(w > 0.0)) {
      if (launch == BTAWL_LAUNCH_WEIBULL ||
          launch == BTAWL_LAUNCH_LOGNORMAL ||
          launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
        return 0.0;
      return posdrift ? 0.0 : pnorm_log_direct(p1 / p2, true);
    }
    if (launch == BTAWL_LAUNCH_WEIBULL)
      return log_weibull_cdf(w, p1, p2);
    if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
      if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
        return log_split_lognormal_cdf(w, p1, p2, delta);
      return pnorm_log_direct((std::log(w) - p1) / p2, true);
    }
    const double l = launch == BTAWL_LAUNCH_WEIBULL
      ? log_weibull_cdf(w, p1, p2)
      : log_normal_cdf_positive_raw(w, p1, p2);
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double h = btawl_hs(t, k, tau_s);
  if (!(h > 0.0)) return 0.0;
  const double E = std::exp(-k * t);
  const double q = E / h;
  if (!(q > 0.0) || !emc2_isfinite(q)) {
    const double w = b / h;
    if (launch == BTAWL_LAUNCH_WEIBULL)
      return log_weibull_cdf(w, p1, p2);
    if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
      if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
        return log_split_lognormal_cdf(w, p1, p2, delta);
      return pnorm_log_direct((std::log(w) - p1) / p2, true);
    }
    const double l = launch == BTAWL_LAUNCH_WEIBULL
      ? log_weibull_cdf(w, p1, p2)
      : log_normal_cdf_positive_raw(w, p1, p2);
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double a = b / h;
  if (A <= BTAWL_A_EPS ||
      (launch == BTAWL_LAUNCH_WEIBULL &&
       btawl_sustained_weib_interval_collapsed(a, q, A))) {
    const double w = b / h;
    if (!(w > 0.0)) return R_NegInf;
    const double l = (launch == BTAWL_LAUNCH_WEIBULL)
      ? log_weibull_cdf(w, p1, p2)
      : (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (launch == BTAWL_LAUNCH_SPLITLOGNORMAL
           ? log_split_lognormal_cdf(w, p1, p2, delta)
           : pnorm_log_direct((std::log(w) - p1) / p2, true))
      : log_normal_cdf_positive_raw(w, p1, p2);
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double whi = a;
  const double wlo = a - q * A;
  if (!(wlo > 0.0) || !(whi >= wlo)) return R_NegInf;

  double li = R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    const double pa = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? log_split_lognormal_put(whi, p1, p2, delta)
      : log_lognormal_put(whi, p1, p2);
    const double pb = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? log_split_lognormal_put(wlo, p1, p2, delta)
      : log_lognormal_put(wlo, p1, p2);
    if (pa - pb > 1e-10) li = log_diff_exp(pa, pb) - std::log(q);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    const double pa = log_weibull_put(whi, p1, p2);
    const double pb = log_weibull_put(wlo, p1, p2);
    if (pa - pb > 1e-10) li = log_diff_exp(pa, pb) - std::log(q);
  } else {
    const double lo = (wlo - p1) / p2;
    const double hi = (whi - p1) / p2;
    const double x = log_normal_phi_integral_positive_raw(lo, hi, p1, p2, posdrift);
    if (x > R_NegInf) li = x + std::log(p2) - std::log(q);
  }
  if (!(li > R_NegInf)) {
    const double wm = 0.5 * (whi + wlo);
    const double lm = (launch == BTAWL_LAUNCH_WEIBULL)
      ? log_weibull_cdf(wm, p1, p2)
      : (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? (launch == BTAWL_LAUNCH_SPLITLOGNORMAL
           ? log_split_lognormal_cdf(wm, p1, p2, delta)
           : pnorm_log_direct((std::log(wm) - p1) / p2, true))
      : log_normal_cdf_positive_raw(wm, p1, p2);
    li = std::log(A) + lm;
  }
  double out = li - std::log(A);
  if (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    out -= log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

double btawl_local_race_log_surv(double t, double A, double b, double p1,
                                 double p2, double k, double tau_s,
                                 double tau_t, double pi, int launch,
                                 bool posdrift, double delta) {
  if (pi <= 1e-14)
    return btawl_log_surv(t, A, b, p1, p2, k, tau_t, launch, posdrift, delta);
  if (pi >= 1.0 - 1e-14)
    return btawl_sustained_log_surv(t, A, b, p1, p2, k, tau_s, launch, posdrift, delta);

  double p1_S = p1, p2_S = p2;
  double p1_T = p1, p2_T = p2;
  if (launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL) {
    p1_S += std::log(pi);
    p1_T += std::log(1.0 - pi);
  } else if (launch == BTAWL_LAUNCH_WEIBULL) {
    p1_S = p1_T = p1;
    p2_S = p2 * pi;
    p2_T = p2 * (1.0 - pi);
  } else {
    p1_S *= pi;         p2_S *= pi;
    p1_T *= (1.0 - pi); p2_T *= (1.0 - pi);
  }

  const double ls_T = btawl_log_surv(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift, delta);
  const double ls_S = btawl_sustained_log_surv(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift, delta);

  if (!(ls_T > R_NegInf) || !(ls_S > R_NegInf)) return R_NegInf;
  return ls_T + ls_S;
}

// ============================================================
// BTAwL transient-member adapters. The full local race has a separate
// nine-column adapter below; both use p1=0 (v | mu), p2=1 (sv | sigma).
// ============================================================
namespace {

inline int btawl_launch_of(const ContextForRaceModels* ctx) {
  return ctx ? ctx->btawl_launch : BTAWL_LAUNCH_NORMAL;
}

inline bool btawl_uses_ttrans(const ContextForRaceModels* ctx) {
  return ctx && ctx->btawl_ttrans_chart;
}

inline double btawl_tau_of(const ContextForRaceModels* ctx, double clear,
                           double k) {
  if (!btawl_uses_ttrans(ctx)) return clear;
  if (ctx != nullptr) {
    for (int j = 0; j < ContextForRaceModels::btawl_tau_cache_size; ++j) {
      if (ctx->btawl_tau_cache_k[j] == k &&
          ctx->btawl_tau_cache_clear[j] == clear)
        return ctx->btawl_tau_cache_value[j];
    }
  }
  const double tau = btawl_tau_from_ttrans(k, clear);
  if (ctx != nullptr) {
    const int j = ctx->btawl_tau_cache_next;
    ctx->btawl_tau_cache_k[j] = k;
    ctx->btawl_tau_cache_clear[j] = clear;
    ctx->btawl_tau_cache_value[j] = tau;
    ctx->btawl_tau_cache_next =
      (j + 1) % ContextForRaceModels::btawl_tau_cache_size;
  }
  return tau;
}
} // namespace

double dbtawl_transient_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = split ? int(emc2col::btawlsplit_transient::mu) : int(emc2col::btawl_transient::v);
  const int isv = split ? int(emc2col::btawlsplit_transient::sigma) : int(emc2col::btawl_transient::sv);
  const int iB = split ? int(emc2col::btawlsplit_transient::B) : int(emc2col::btawl_transient::B);
  const int iA = split ? int(emc2col::btawlsplit_transient::A) : int(emc2col::btawl_transient::A);
  const int it0 = split ? int(emc2col::btawlsplit_transient::t0) : int(emc2col::btawl_transient::t0);
  const int ik = split ? int(emc2col::btawlsplit_transient::k) : int(emc2col::btawl_transient::k);
  const int iclear = split ? int(emc2col::btawlsplit_transient::clear) : int(emc2col::btawl_transient::clear);
  const double delta = split ? par[emc2col::btawlsplit_transient::delta] : 0.0;
  if (R_IsNA(par[iv])) return 0.0;
  const double tt = t - par[it0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[iB] + par[iA];
  if (btawl_uses_ttrans(ctx))
    return btawl_pdf_chart(tt, par[iA], b,
                           par[iv], par[isv],
                           par[ik], par[iclear],
                           launch, ctx ? ctx->use_posdrift : true, true,
                           btawl_tau_of(ctx, par[iclear], par[ik]), delta);
  return btawl_pdf(tt, par[iA], b,
                   par[iv], par[isv],
                   par[ik],
                   btawl_tau_of(ctx, par[iclear], par[ik]),
                   launch, ctx ? ctx->use_posdrift : true, delta);
}

double pbtawl_transient_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = split ? int(emc2col::btawlsplit_transient::mu) : int(emc2col::btawl_transient::v);
  const int isv = split ? int(emc2col::btawlsplit_transient::sigma) : int(emc2col::btawl_transient::sv);
  const int iB = split ? int(emc2col::btawlsplit_transient::B) : int(emc2col::btawl_transient::B);
  const int iA = split ? int(emc2col::btawlsplit_transient::A) : int(emc2col::btawl_transient::A);
  const int it0 = split ? int(emc2col::btawlsplit_transient::t0) : int(emc2col::btawl_transient::t0);
  const int ik = split ? int(emc2col::btawlsplit_transient::k) : int(emc2col::btawl_transient::k);
  const int iclear = split ? int(emc2col::btawlsplit_transient::clear) : int(emc2col::btawl_transient::clear);
  const double delta = split ? par[emc2col::btawlsplit_transient::delta] : 0.0;
  if (R_IsNA(par[iv])) return 0.0;
  const double tt = t - par[it0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[iB] + par[iA];
  if (btawl_uses_ttrans(ctx))
    return btawl_cdf_chart(tt, par[iA], b,
                           par[iv], par[isv],
                           par[ik], par[iclear],
                           launch, ctx ? ctx->use_posdrift : true, true,
                           btawl_tau_of(ctx, par[iclear], par[ik]), delta);
  return btawl_cdf(tt, par[iA], b,
                   par[iv], par[isv],
                   par[ik],
                   btawl_tau_of(ctx, par[iclear], par[ik]),
                   launch, ctx ? ctx->use_posdrift : true, delta);
}

void dbtawl_transient_raw(const double* rt, const double* const* cols, int n_rows,
                          const int* mask, const int* isok, double* out,
                          double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_transient::mu) : int(emc2col::btawl_transient::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_transient::sigma) : int(emc2col::btawl_transient::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_transient::B) : int(emc2col::btawl_transient::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_transient::A) : int(emc2col::btawl_transient::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_transient::t0) : int(emc2col::btawl_transient::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_transient::k) : int(emc2col::btawl_transient::k)];
  const double* clear = cols[split ? int(emc2col::btawlsplit_transient::clear) : int(emc2col::btawl_transient::clear)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_transient::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double tt = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(tt > 0.0)) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double tau = btawl_tau_of(ctx, clear[i], k[i]);
    const BtawlGeom g = btawl_geometry_cached(ctx, btawl_uses_ttrans(ctx),
                                              A[i], B[i] + A[i], k[i],
                                              clear[i], tau);
    const double lp = btawl_log_pdf_from_geom(tt, g, p1[i], p2[i], launch, pd,
                                              delta_ ? delta_[i] : 0.0);
    out[i] = (lp > R_NegInf && emc2_isfinite(lp))
      ? raw_log_value(lp, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

void pbtawl_transient_raw(const double* rt, const double* const* cols, int n_rows,
                          const int* mask, const int* isok, double* out,
                          double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_transient::mu) : int(emc2col::btawl_transient::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_transient::sigma) : int(emc2col::btawl_transient::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_transient::B) : int(emc2col::btawl_transient::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_transient::A) : int(emc2col::btawl_transient::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_transient::t0) : int(emc2col::btawl_transient::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_transient::k) : int(emc2col::btawl_transient::k)];
  const double* clear = cols[split ? int(emc2col::btawlsplit_transient::clear) : int(emc2col::btawl_transient::clear)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_transient::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(tt > 0.0)) { out[i] = 0.0; continue; }
    const double tau = btawl_tau_of(ctx, clear[i], k[i]);
    const BtawlGeom g = btawl_geometry_cached(ctx, btawl_uses_ttrans(ctx),
                                              A[i], B[i] + A[i], k[i],
                                              clear[i], tau);
    const double d_val = delta_ ? delta_[i] : 0.0;
    if (emc2_isfinite(g.t_max) && tt >= g.t_max) {
      const double ls = btawl_log_surv_cached(ctx, tt, g, clear[i], p1[i],
                                              p2[i], launch, pd, d_val);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
      continue;
    }
    double cdf = 0.0;
    if (btawl_natural_cdf_from_geom(tt, g, p1[i], p2[i], launch, pd, cdf, d_val)) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = btawl_log_surv_cached(ctx, tt, g, clear[i], p1[i],
                                              p2[i], launch, pd, d_val);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

void btawl_transient_logS_at_t(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_transient::mu) : int(emc2col::btawl_transient::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_transient::sigma) : int(emc2col::btawl_transient::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_transient::B) : int(emc2col::btawl_transient::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_transient::A) : int(emc2col::btawl_transient::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_transient::t0) : int(emc2col::btawl_transient::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_transient::k) : int(emc2col::btawl_transient::k)];
  const double* clear = cols[split ? int(emc2col::btawlsplit_transient::clear) : int(emc2col::btawl_transient::clear)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_transient::delta] : nullptr;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double ls = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = j * n_lR + kk;
      if (r >= n_rows_total || !isok_all[r] || R_IsNA(p1[r])) { bad = true; break; }
      const double tt = t - t0[r];
      if (!(tt > 0.0)) continue;
      const double tau = btawl_tau_of(ctx, clear[r], k[r]);
      const BtawlGeom g = btawl_geometry_cached(ctx, btawl_uses_ttrans(ctx),
                                                A[r], B[r] + A[r], k[r],
                                                clear[r], tau);
      const double lsr = btawl_log_surv_cached(ctx, tt, g, clear[r], p1[r],
                                               p2[r], launch, pd,
                                               delta_ ? delta_[r] : 0.0);
      if (!(lsr > R_NegInf) || ISNAN(lsr)) { bad = true; break; }
      ls += lsr;
    }
    logS_out[j] = bad ? R_NegInf : ls;
  }
}

// BTAwL sustained/transient local race.
double dbtawl_local_race_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = split ? int(emc2col::btawlsplit_local_race::mu) : int(emc2col::btawl_local_race::v);
  const int isv = split ? int(emc2col::btawlsplit_local_race::sigma) : int(emc2col::btawl_local_race::sv);
  const int iB = split ? int(emc2col::btawlsplit_local_race::B) : int(emc2col::btawl_local_race::B);
  const int iA = split ? int(emc2col::btawlsplit_local_race::A) : int(emc2col::btawl_local_race::A);
  const int it0 = split ? int(emc2col::btawlsplit_local_race::t0) : int(emc2col::btawl_local_race::t0);
  const int ik = split ? int(emc2col::btawlsplit_local_race::k) : int(emc2col::btawl_local_race::k);
  const int its = split ? int(emc2col::btawlsplit_local_race::tau_s) : int(emc2col::btawl_local_race::tau_s);
  const int iclear = split ? int(emc2col::btawlsplit_local_race::clear) : int(emc2col::btawl_local_race::clear);
  const int ipi = split ? int(emc2col::btawlsplit_local_race::pi) : int(emc2col::btawl_local_race::pi);
  const double delta = split ? par[emc2col::btawlsplit_local_race::delta] : 0.0;
  if (R_IsNA(par[iv])) return 0.0;
  const double tt = t - par[it0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[iB] + par[iA];
  const double tau_t = btawl_tau_of(ctx, par[iclear], par[ik]);
  if (par[ipi] <= 1e-14) {
    const BtawlGeom g = btawl_geometry_cached(
      ctx, btawl_uses_ttrans(ctx), par[iA], b,
      par[ik], par[iclear], tau_t);
    const double lp = btawl_log_pdf_from_geom(tt, g, par[iv], par[isv],
                                              launch, ctx ? ctx->use_posdrift : true, delta);
    return (lp > R_NegInf && lp < 709.0) ? std::exp(lp) : 0.0;
  }
  if (par[ipi] >= 1.0 - 1e-14) {
    const double lp = btawl_sustained_log_pdf(tt, par[iA], b,
                                              par[iv], par[isv], par[ik],
                                              par[its], launch,
                                              ctx ? ctx->use_posdrift : true,
                                              delta);
    return (lp > R_NegInf && lp < 709.0) ? std::exp(lp) : 0.0;
  }
  return btawl_local_race_pdf(tt, par[iA], b,
                              par[iv], par[isv],
                              par[ik], par[its],
                              tau_t, par[ipi],
                              launch, ctx ? ctx->use_posdrift : true, delta);
}

double pbtawl_local_race_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = split ? int(emc2col::btawlsplit_local_race::mu) : int(emc2col::btawl_local_race::v);
  const int isv = split ? int(emc2col::btawlsplit_local_race::sigma) : int(emc2col::btawl_local_race::sv);
  const int iB = split ? int(emc2col::btawlsplit_local_race::B) : int(emc2col::btawl_local_race::B);
  const int iA = split ? int(emc2col::btawlsplit_local_race::A) : int(emc2col::btawl_local_race::A);
  const int it0 = split ? int(emc2col::btawlsplit_local_race::t0) : int(emc2col::btawl_local_race::t0);
  const int ik = split ? int(emc2col::btawlsplit_local_race::k) : int(emc2col::btawl_local_race::k);
  const int its = split ? int(emc2col::btawlsplit_local_race::tau_s) : int(emc2col::btawl_local_race::tau_s);
  const int iclear = split ? int(emc2col::btawlsplit_local_race::clear) : int(emc2col::btawl_local_race::clear);
  const int ipi = split ? int(emc2col::btawlsplit_local_race::pi) : int(emc2col::btawl_local_race::pi);
  const double delta = split ? par[emc2col::btawlsplit_local_race::delta] : 0.0;
  if (R_IsNA(par[iv])) return 0.0;
  const double tt = t - par[it0];
  if (!(t > 0.0) || !(tt > 0.0)) return 0.0;
  const double b = par[iB] + par[iA];
  const double tau_t = btawl_tau_of(ctx, par[iclear], par[ik]);
  if (par[ipi] <= 1e-14) {
    const BtawlGeom g = btawl_geometry_cached(
      ctx, btawl_uses_ttrans(ctx), par[iA], b,
      par[ik], par[iclear], tau_t);
    return btawl_cdf_from_geom(tt, g, par[iv], par[isv],
                               launch, ctx ? ctx->use_posdrift : true, delta);
  }
  return btawl_local_race_cdf(tt, par[iA], b,
                              par[iv], par[isv],
                              par[ik], par[its],
                              tau_t, par[ipi],
                              launch, ctx ? ctx->use_posdrift : true, delta);
}

void dbtawl_local_race_raw(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok, double* out,
                           double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx); const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_local_race::mu) : int(emc2col::btawl_local_race::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_local_race::sigma) : int(emc2col::btawl_local_race::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_local_race::B) : int(emc2col::btawl_local_race::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_local_race::A) : int(emc2col::btawl_local_race::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_local_race::t0) : int(emc2col::btawl_local_race::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_local_race::k) : int(emc2col::btawl_local_race::k)];
  const double* ts = cols[split ? int(emc2col::btawlsplit_local_race::tau_s) : int(emc2col::btawl_local_race::tau_s)];
  const double* tt_clear = cols[split ? int(emc2col::btawlsplit_local_race::clear) : int(emc2col::btawl_local_race::clear)];
  const double* pi = cols[split ? int(emc2col::btawlsplit_local_race::pi) : int(emc2col::btawl_local_race::pi)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_local_race::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tau_t = btawl_tau_of(ctx, tt_clear[i], k[i]);
    const double d_val = delta_ ? delta_[i] : 0.0;
    double lp = R_NegInf;
    if (pi[i] <= 1e-14) {
      const BtawlGeom g = btawl_geometry_cached(ctx, btawl_uses_ttrans(ctx),
                                                A[i], B[i] + A[i], k[i],
                                                tt_clear[i], tau_t);
      lp = btawl_log_pdf_from_geom(u, g, p1[i], p2[i], launch, pd, d_val);
    } else {
      lp = btawl_local_race_pdf_log(u, A[i], B[i] + A[i], p1[i], p2[i], k[i],
                                    ts[i], tau_t, pi[i], launch, pd, d_val);
    }
    out[i] = (lp > R_NegInf && emc2_isfinite(lp))
      ? raw_log_value(lp, min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

void pbtawl_local_race_raw(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok, double* out,
                           double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const int launch = btawl_launch_of(ctx); const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_local_race::mu) : int(emc2col::btawl_local_race::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_local_race::sigma) : int(emc2col::btawl_local_race::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_local_race::B) : int(emc2col::btawl_local_race::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_local_race::A) : int(emc2col::btawl_local_race::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_local_race::t0) : int(emc2col::btawl_local_race::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_local_race::k) : int(emc2col::btawl_local_race::k)];
  const double* ts = cols[split ? int(emc2col::btawlsplit_local_race::tau_s) : int(emc2col::btawl_local_race::tau_s)];
  const double* tt_clear = cols[split ? int(emc2col::btawlsplit_local_race::clear) : int(emc2col::btawl_local_race::clear)];
  const double* pi = cols[split ? int(emc2col::btawlsplit_local_race::pi) : int(emc2col::btawl_local_race::pi)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_local_race::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = 0.0; continue; }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) { out[i] = 0.0; continue; }
    const double tau_t = btawl_tau_of(ctx, tt_clear[i], k[i]);
    const double d_val = delta_ ? delta_[i] : 0.0;
    double cdf = 0.0;
    bool cdf_ready = false;
    if (pi[i] <= 1e-14) {
      const BtawlGeom g = btawl_geometry_cached(ctx, btawl_uses_ttrans(ctx),
                                                A[i], B[i] + A[i], k[i],
                                                tt_clear[i], tau_t);
      if (emc2_isfinite(g.t_max) && u >= g.t_max) {
        const double ls = btawl_log_surv_cached(ctx, u, g, tt_clear[i],
                                                p1[i], p2[i], launch, pd, d_val);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
        continue;
      }
      cdf_ready = btawl_natural_cdf_from_geom(u, g, p1[i], p2[i], launch, pd, cdf, d_val);
      if (!cdf_ready) {
        const double ls = btawl_log_surv_cached(ctx, u, g, tt_clear[i],
                                                p1[i], p2[i], launch, pd, d_val);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
        continue;
      }
    } else {
      cdf = btawl_local_race_cdf(u, A[i], B[i] + A[i], p1[i], p2[i], k[i],
                                 ts[i], tau_t, pi[i], launch, pd, d_val);
      cdf_ready = true;
    }
    if (cdf_ready && emc2_isfinite(cdf) && cdf >= 0.0 && cdf < 1.0 - 1e-8) {
      out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = btawl_local_race_log_surv(u, A[i], B[i] + A[i], p1[i], p2[i], k[i],
                                                  ts[i], tau_t, pi[i], launch, pd, d_val);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

void btawl_local_race_logS_at_t(double t, const double* const* cols,
                                int n_rows_total, int n_lR, int /*n_par*/,
                                const int* trunc_mask, int n_unique_trials,
                                const int* isok_all, void* ctx_, double* out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx); const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_local_race::mu) : int(emc2col::btawl_local_race::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_local_race::sigma) : int(emc2col::btawl_local_race::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_local_race::B) : int(emc2col::btawl_local_race::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_local_race::A) : int(emc2col::btawl_local_race::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_local_race::t0) : int(emc2col::btawl_local_race::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_local_race::k) : int(emc2col::btawl_local_race::k)];
  const double* ts = cols[split ? int(emc2col::btawlsplit_local_race::tau_s) : int(emc2col::btawl_local_race::tau_s)];
  const double* tt_clear = cols[split ? int(emc2col::btawlsplit_local_race::clear) : int(emc2col::btawl_local_race::clear)];
  const double* pi = cols[split ? int(emc2col::btawlsplit_local_race::pi) : int(emc2col::btawl_local_race::pi)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_local_race::delta] : nullptr;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double ls = 0.0; bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = j * n_lR + kk;
      if (r >= n_rows_total || !isok_all[r] || R_IsNA(p1[r])) { bad = true; break; }
      const double u = t - t0[r]; if (!(u > 0.0)) continue;
      const double tau_t = btawl_tau_of(ctx, tt_clear[r], k[r]);
      const double d_val = delta_ ? delta_[r] : 0.0;
      double lsr = R_NegInf;
      if (pi[r] <= 1e-14) {
        const BtawlGeom g = btawl_geometry_cached(ctx, btawl_uses_ttrans(ctx),
                                                  A[r], B[r] + A[r], k[r],
                                                  tt_clear[r], tau_t);
        lsr = btawl_log_surv_cached(ctx, u, g, tt_clear[r], p1[r], p2[r],
                                    launch, pd, d_val);
      } else {
        lsr = btawl_local_race_log_surv(u, A[r], B[r] + A[r], p1[r], p2[r], k[r],
                                        ts[r], tau_t, pi[r], launch, pd, d_val);
      }
      if (!(lsr > R_NegInf) || ISNAN(lsr)) { bad = true; break; }
      ls += lsr;
    }
    out[j] = bad ? R_NegInf : ls;
  }
}

double dbtawl_sustained_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = split ? int(emc2col::btawlsplit_sustained::mu) : int(emc2col::btawl_sustained::v);
  const int isv = split ? int(emc2col::btawlsplit_sustained::sigma) : int(emc2col::btawl_sustained::sv);
  const int iB = split ? int(emc2col::btawlsplit_sustained::B) : int(emc2col::btawl_sustained::B);
  const int iA = split ? int(emc2col::btawlsplit_sustained::A) : int(emc2col::btawl_sustained::A);
  const int it0 = split ? int(emc2col::btawlsplit_sustained::t0) : int(emc2col::btawl_sustained::t0);
  const int ik = split ? int(emc2col::btawlsplit_sustained::k) : int(emc2col::btawl_sustained::k);
  const int its = split ? int(emc2col::btawlsplit_sustained::tau_s) : int(emc2col::btawl_sustained::tau_s);
  const double delta = split ? par[emc2col::btawlsplit_sustained::delta] : 0.0;
  const double tt = t - par[it0];
  if (R_IsNA(par[iv]) || !(t > 0.0) || !(tt > 0.0))
    return 0.0;
  return btawl_sustained_pdf(tt, par[iA],
                             par[iB] + par[iA],
                             par[iv],
                             par[isv],
                             par[ik],
                             par[its], launch,
                             ctx ? ctx->use_posdrift : true, delta);
}

double pbtawl_sustained_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = split ? int(emc2col::btawlsplit_sustained::mu) : int(emc2col::btawl_sustained::v);
  const int isv = split ? int(emc2col::btawlsplit_sustained::sigma) : int(emc2col::btawl_sustained::sv);
  const int iB = split ? int(emc2col::btawlsplit_sustained::B) : int(emc2col::btawl_sustained::B);
  const int iA = split ? int(emc2col::btawlsplit_sustained::A) : int(emc2col::btawl_sustained::A);
  const int it0 = split ? int(emc2col::btawlsplit_sustained::t0) : int(emc2col::btawl_sustained::t0);
  const int ik = split ? int(emc2col::btawlsplit_sustained::k) : int(emc2col::btawl_sustained::k);
  const int its = split ? int(emc2col::btawlsplit_sustained::tau_s) : int(emc2col::btawl_sustained::tau_s);
  const double delta = split ? par[emc2col::btawlsplit_sustained::delta] : 0.0;
  const double tt = t - par[it0];
  if (R_IsNA(par[iv]) || !(t > 0.0) || !(tt > 0.0))
    return 0.0;
  return btawl_sustained_cdf(tt, par[iA],
                             par[iB] + par[iA],
                             par[iv],
                             par[isv],
                             par[ik],
                             par[its], launch,
                             ctx ? ctx->use_posdrift : true, delta);
}

void dbtawl_sustained_raw(const double* rt, const double* const* cols,
                          int n_rows, const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_sustained::mu) : int(emc2col::btawl_sustained::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_sustained::sigma) : int(emc2col::btawl_sustained::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_sustained::B) : int(emc2col::btawl_sustained::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_sustained::A) : int(emc2col::btawl_sustained::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_sustained::t0) : int(emc2col::btawl_sustained::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_sustained::k) : int(emc2col::btawl_sustained::k)];
  const double* ts = cols[split ? int(emc2col::btawlsplit_sustained::tau_s) : int(emc2col::btawl_sustained::tau_s)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_sustained::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    const double p = btawl_sustained_pdf(u, A[i], B[i] + A[i], p1[i], p2[i],
                                         k[i], ts[i], launch, pd,
                                         delta_ ? delta_[i] : 0.0);
    out[i] = (p > 0.0 && emc2_isfinite(p))
      ? raw_log_value(std::log(p), min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

void pbtawl_sustained_raw(const double* rt, const double* const* cols,
                          int n_rows, const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_sustained::mu) : int(emc2col::btawl_sustained::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_sustained::sigma) : int(emc2col::btawl_sustained::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_sustained::B) : int(emc2col::btawl_sustained::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_sustained::A) : int(emc2col::btawl_sustained::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_sustained::t0) : int(emc2col::btawl_sustained::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_sustained::k) : int(emc2col::btawl_sustained::k)];
  const double* ts = cols[split ? int(emc2col::btawlsplit_sustained::tau_s) : int(emc2col::btawl_sustained::tau_s)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_sustained::delta] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || R_IsNA(p1[i])) { out[i] = 0.0; continue; }
    const double u = rt[i] - t0[i];
    if (!(rt[i] > 0.0) || !(u > 0.0)) { out[i] = 0.0; continue; }
    const double d_val = delta_ ? delta_[i] : 0.0;
    const double cdf = btawl_sustained_cdf(u, A[i], B[i] + A[i], p1[i], p2[i],
                                           k[i], ts[i], launch, pd, d_val);
    if (emc2_isfinite(cdf) && cdf >= 0.0 && cdf < 1.0 - 1e-8) {
      out[i] = cdf > 0.0 ? std::log1p(-cdf) : 0.0;
    } else {
      const double ls = btawl_sustained_log_surv(u, A[i], B[i] + A[i], p1[i], p2[i],
                                                 k[i], ts[i], launch, pd, d_val);
      out[i] = (ls > R_NegInf && emc2_isfinite(ls))
        ? ls : raw_log_zero(min_ll, floor_raw);
    }
  }
}

void btawl_sustained_logS_at_t(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int /*n_par*/,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int launch = btawl_launch_of(ctx);
  const bool pd = ctx ? ctx->use_posdrift : true;
  const bool split = launch == BTAWL_LAUNCH_SPLITLOGNORMAL;
  const double* p1 = cols[split ? int(emc2col::btawlsplit_sustained::mu) : int(emc2col::btawl_sustained::v)];
  const double* p2 = cols[split ? int(emc2col::btawlsplit_sustained::sigma) : int(emc2col::btawl_sustained::sv)];
  const double* B = cols[split ? int(emc2col::btawlsplit_sustained::B) : int(emc2col::btawl_sustained::B)];
  const double* A = cols[split ? int(emc2col::btawlsplit_sustained::A) : int(emc2col::btawl_sustained::A)];
  const double* t0 = cols[split ? int(emc2col::btawlsplit_sustained::t0) : int(emc2col::btawl_sustained::t0)];
  const double* k = cols[split ? int(emc2col::btawlsplit_sustained::k) : int(emc2col::btawl_sustained::k)];
  const double* ts = cols[split ? int(emc2col::btawlsplit_sustained::tau_s) : int(emc2col::btawl_sustained::tau_s)];
  const double* delta_ = split ? cols[emc2col::btawlsplit_sustained::delta] : nullptr;
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double ls = 0.0; bool bad = false;
    for (int kk = 0; kk < n_lR; ++kk) {
      const int r = j * n_lR + kk;
      if (r >= n_rows_total || !isok_all[r] || R_IsNA(p1[r])) { bad = true; break; }
      const double u = t - t0[r]; if (!(u > 0.0)) continue;
      const double lsr = btawl_sustained_log_surv(u, A[r], B[r] + A[r], p1[r], p2[r],
                                                  k[r], ts[r], launch, pd,
                                                  delta_ ? delta_[r] : 0.0);
      if (!(lsr > R_NegInf) || ISNAN(lsr)) { bad = true; break; }
      ls += lsr;
    }
    out[j] = bad ? R_NegInf : ls;
  }
}

// [[Rcpp::export]]
NumericVector dbtawl_transient(NumericVector t, NumericVector A, NumericVector b,
                              NumericVector p1, NumericVector p2, NumericVector k,
                              NumericVector tau, int launch = 1,
                              bool posdrift = true, bool log_out = false,
                              NumericVector delta = 0.0) {
  const int n = t.size(); NumericVector out(n);
  ContextForRaceModels cache_ctx;
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double Ai = pick(A, i), bi = pick(b, i), ki = pick(k, i), ti = pick(tau, i);
    const BtawlGeom g = btawl_geometry_cached(&cache_ctx, false, Ai, bi, ki, ti, ti);
    const double lp = btawl_log_pdf_from_geom(t[i], g, pick(p1,i), pick(p2,i),
                                              launch, posdrift, pick(delta, i));
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pbtawl_transient(NumericVector t, NumericVector A, NumericVector b,
                              NumericVector p1, NumericVector p2, NumericVector k,
                              NumericVector tau, int launch = 1,
                              bool posdrift = true, bool log_out = false,
                              NumericVector delta = 0.0) {
  const int n = t.size(); NumericVector out(n);
  ContextForRaceModels cache_ctx;
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double Ai = pick(A, i), bi = pick(b, i), ki = pick(k, i), ti = pick(tau, i);
    const BtawlGeom g = btawl_geometry_cached(&cache_ctx, false, Ai, bi, ki, ti, ti);
    const double lp = (launch == BTAWL_LAUNCH_WEIBULL)
      ? log_btawl_cdf_weib(t[i], g, pick(p1,i), pick(p2,i))
      : ((launch == BTAWL_LAUNCH_LOGNORMAL || launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      ? log_btawl_cdf_logn(t[i], g, pick(p1,i), pick(p2,i), pick(delta, i))
      : log_btawl_cdf_normal(t[i], g, pick(p1,i), pick(p2,i), posdrift));
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_transient_log_surv_vec(NumericVector t, NumericVector A,
                                          NumericVector b, NumericVector p1,
                                          NumericVector p2, NumericVector k,
                                          NumericVector tau, int launch = 1,
                                          bool posdrift = true,
                                          NumericVector delta = 0.0) {
  const int n = std::max({t.size(), A.size(), b.size(), p1.size(), p2.size(),
                          k.size(), tau.size(), delta.size()});
  NumericVector out(n);
  ContextForRaceModels cache_ctx;
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double Ai = pick(A, i), bi = pick(b, i), ki = pick(k, i), ti = pick(tau, i);
    const BtawlGeom g = btawl_geometry_cached(&cache_ctx, false, Ai, bi, ki, ti, ti);
    out[i] = btawl_log_surv_cached(&cache_ctx, pick(t, i), g, ti,
                                   pick(p1, i), pick(p2, i), launch, posdrift,
                                   pick(delta, i));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_local_race_log_surv_vec(NumericVector t, NumericVector A,
                                           NumericVector b, NumericVector p1,
                                           NumericVector p2, NumericVector k,
                                           NumericVector tau_s, NumericVector tau_t,
                                           NumericVector pi, int launch = 1,
                                           bool posdrift = true,
                                           NumericVector delta = 0.0) {
  const int n = std::max({t.size(), A.size(), b.size(), p1.size(), p2.size(),
                          k.size(), tau_s.size(), tau_t.size(), pi.size(),
                          delta.size()});
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i)
    out[i] = btawl_local_race_log_surv(pick(t, i), pick(A, i), pick(b, i), pick(p1, i),
                                      pick(p2, i), pick(k, i), pick(tau_s, i),
                                      pick(tau_t, i), pick(pi, i), launch, posdrift,
                                      pick(delta, i));
  return out;
}

// [[Rcpp::export]]
NumericVector dbtawl_local_race(NumericVector t, NumericVector A, NumericVector b,
                               NumericVector p1, NumericVector p2, NumericVector k,
                               NumericVector tau_s, NumericVector tau_t,
                               NumericVector pi, int launch = 1,
                               bool posdrift = true, bool log_out = false,
                               NumericVector delta = 0.0) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_local_race_pdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i), pick(p2,i),
                                               pick(k,i), pick(tau_s,i), pick(tau_t,i), pick(pi,i),
                                               launch, posdrift, pick(delta, i));
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pbtawl_local_race(NumericVector t, NumericVector A, NumericVector b,
                               NumericVector p1, NumericVector p2, NumericVector k,
                               NumericVector tau_s, NumericVector tau_t,
                               NumericVector pi, int launch = 1,
                               bool posdrift = true, bool log_out = false,
                               NumericVector delta = 0.0) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_local_race_cdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i), pick(p2,i),
                                               pick(k,i), pick(tau_s,i), pick(tau_t,i), pick(pi,i),
                                               launch, posdrift, pick(delta, i));
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_tmax_vec(NumericVector k, NumericVector tau) {
  const int n = std::max(k.size(), tau.size()); NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = btawl_tmax(k.size() == 1 ? k[0] : k[i],
                                                 tau.size() == 1 ? tau[0] : tau[i]);
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_tau_vec(NumericVector k, NumericVector Ttrans) {
  const int n = std::max(k.size(), Ttrans.size()); NumericVector out(n);
  for (int i = 0; i < n; ++i)
    out[i] = btawl_tau_from_ttrans(k.size() == 1 ? k[0] : k[i],
                                   Ttrans.size() == 1 ? Ttrans[0] : Ttrans[i]);
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_vcrit_vec(NumericVector k, NumericVector tau,
                              NumericVector b) {
  const int n = std::max({k.size(), tau.size(), b.size()}); NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ki = k.size() == 1 ? k[0] : k[i];
    const double ti = tau.size() == 1 ? tau[0] : tau[i];
    const double bi = b.size() == 1 ? b[0] : b[i];
    const double tm = btawl_tmax(ki, ti);
    const double h = btawl_h(tm, ki, ti);
    out[i] = (h > 0.0 && emc2_isfinite(h)) ? bi / h : R_PosInf;
  }
  return out;
}

#ifndef EMC2_MODEL_BTAWL_H
#define EMC2_MODEL_BTAWL_H

// BTAwL: ballistic transient accumulator with leak.
//
//   dX/du = V (u/tau) exp(-u/tau) - k X,  X(0) = z, z ~ U(0,A)
//
// The transient response is affine in the launch point and the required
// launch V*(u,z) is therefore affine in z on the live part of the CDF.  The
// live start-point average is evaluated in closed form (normal and lognormal
// launches).  The already-saturated start region is a one-dimensional
// quadrature in tangency time; there is no grid or PDE solve.
//
// This header is included by model_LBA.h, after the shared normal/lognormal
// primitives and the GSL helpers have been declared.

#include <cmath>
#include <functional>
#include <algorithm>
#include <limits>

// Keep BTAwL's launch selector and numerical thresholds identical to BAwL.
constexpr int BTAWL_LAUNCH_NORMAL = BAWL_LAUNCH_NORMAL;
constexpr int BTAWL_LAUNCH_LOGNORMAL = BAWL_LAUNCH_LOGNORMAL;
constexpr double BTAWL_K_EPS = BAWL_K_EPS;
constexpr double BTAWL_A_EPS = BAWL_A_EPS;
constexpr double BTAWL_DENOM_FLOOR = BAWL_DENOM_FLOOR;

struct BtawlGeom {
  bool ok = false;
  bool k_zero = false;
  double A = 0.0, b = 0.0, k = 0.0, tau = 0.0;
  double t_max = R_PosInf;
  double h_max = 0.0;
  // Tangency time for the upper start point z = A.  Frozen quadratures use
  // this as their lower endpoint; keeping it with the geometry avoids solving
  // the same monotone equation for every evaluation.
  double s_lo = R_PosInf;
};

inline double btawl_h(double t, double k, double tau) {
  if (!(tau > 0.0) || t <= 0.0) return 0.0;
  if (t == R_PosInf) return k <= BTAWL_K_EPS ? tau : 0.0;
  const double x = t / tau;
  if (k <= BTAWL_K_EPS)
    return tau * (-std::expm1(-x) - x * std::exp(-x));
  const double d = k - 1.0 / tau;
  const double y = d * t;
  // H = exp(-k t) / (tau d^2) * ((y - 1)e^y + 1).  The expm1 form avoids
  // cancellation around the coincident-timescale case k = 1/tau.
  double numer;
  if (std::fabs(y) < 1e-4) {
    const double y2 = y * y;
    numer = y2 * (0.5 + y * (1.0 / 3.0 + y * (1.0 / 8.0 + y / 30.0)));
  } else {
    // Use the original two-decay expression away from the coincident
    // timescale.  Rewriting it with expm1(y) is attractive near y = 0 but
    // can overflow for k > 1/tau even though H itself is tiny and finite.
    const double e_tau = std::exp(-t / tau);
    const double e_k = std::exp(-k * t);
    return (e_tau * (d * t - 1.0) + e_k) / (tau * d * d);
  }
  if (!(numer > 0.0)) {
    // Exact limit, also used as a guard if fast-math rounds a tiny positive
    // numerator down to zero.
    return 0.5 * t * t / tau * std::exp(-t / tau);
  }
  return std::exp(-k * t) * numer / (tau * d * d);
}

inline double btawl_g(double t, double tau) {
  if (!(tau > 0.0) || !(t >= 0.0) || t == R_PosInf) return 0.0;
  return (t / tau) * std::exp(-t / tau);
}

inline double btawl_hp(double t, double k, double tau) {
  const double h = btawl_h(t, k, tau);
  return btawl_g(t, tau) - k * h;
}

inline double btawl_tmax(double k, double tau) {
  if (!(k > BTAWL_K_EPS) || !(tau > 0.0)) return R_PosInf;
  // phi(u) = H'(u) has its unique root above tau.  The analytic seed is
  // already close in both the stiff and diffuse regimes; retain a bracket so
  // the Newton step can never jump across the root.
  double lo = tau;
  double hi = tau + 1.0 / k;
  if (!(hi > lo) || !emc2_isfinite(hi)) hi = 2.0 * lo;
  for (int i = 0; i < 80 && btawl_hp(hi, k, tau) > 0.0; ++i) hi = hi + (hi - lo);
  if (!(btawl_hp(lo, k, tau) >= 0.0) || !(btawl_hp(hi, k, tau) <= 0.0) ||
      !emc2_isfinite(hi)) return R_PosInf;
  double x = std::fmin(std::fmax(tau * (1.0 + 1.0 / (k * tau)), lo), hi);
  for (int i = 0; i < 80; ++i) {
    const double f = btawl_hp(x, k, tau);
    if (f > 0.0) lo = x; else hi = x;
    const double gg = btawl_g(x, tau);
    const double gp = (x > 0.0) ? gg * (1.0 / x - 1.0 / tau) : 0.0;
    const double d = gp - k * f;
    double xn = (emc2_isfinite(d) && d < 0.0) ? x - f / d : 0.5 * (lo + hi);
    if (!(xn > lo) || !(xn < hi) || !emc2_isfinite(xn)) xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) <= 2e-15 * std::fmax(1.0, std::fabs(x))) return xn;
    x = xn;
  }
  return 0.5 * (lo + hi);
}

inline double btawl_tangent_time(double z, const BtawlGeom& g);

// Invert the endpoint chart Ttrans = T_max(k, tau).  The endpoint is strictly
// increasing in tau; the bracket (0, Ttrans) is guaranteed by the transient
// geometry, so a safeguarded bisection is sufficient and has no branch
// ambiguity in the stiff k*tau regime.
inline double btawl_tau_from_ttrans(double k, double Ttrans) {
  if (ISNAN(k) || ISNAN(Ttrans)) return R_NaN;
  if (!(k > BTAWL_K_EPS) || !(Ttrans > 0.0) || !emc2_isfinite(Ttrans))
    return (Ttrans == R_PosInf && k <= BTAWL_K_EPS) ? R_PosInf : R_NaN;
  // The inverse is monotone in tau.  Bisect in log(tau), rather than from
  // denorm_min, so the bracket does not spend half its iterations walking
  // through decades that cannot contain the root.
  double hi = Ttrans;
  if (!(btawl_tmax(k, hi) >= Ttrans)) hi = std::nextafter(Ttrans, R_PosInf);
  double lo = hi;
  for (int i = 0; i < 128; ++i) {
    const double next = lo * 0.5;
    if (!(next > 0.0) || next == lo) break;
    lo = next;
    if (btawl_tmax(k, lo) < Ttrans) break;
  }
  double log_lo = std::log(lo), log_hi = std::log(hi);
  for (int i = 0; i < 80; ++i) {
    const double log_mid = 0.5 * (log_lo + log_hi);
    const double mid = std::exp(log_mid);
    const double tm = btawl_tmax(k, mid);
    if (!(tm > 0.0) || tm < Ttrans) log_lo = log_mid; else log_hi = log_mid;
  }
  return std::exp(0.5 * (log_lo + log_hi));
}

inline BtawlGeom btawl_geometry(double A, double b, double k, double tau) {
  BtawlGeom g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(k) ||
      !emc2_isfinite(tau) || !(A >= 0.0) || !(b > 0.0) || !(b >= A) ||
      !(k >= 0.0) || !(tau > 0.0)) return g;
  g.ok = true; g.A = A; g.b = b; g.k = k; g.tau = tau;
  g.k_zero = k <= BTAWL_K_EPS;
  g.t_max = btawl_tmax(k, tau);
  g.h_max = g.k_zero ? tau : btawl_h(g.t_max, k, tau);
  g.s_lo = g.k_zero ? R_PosInf : btawl_tangent_time(A, g);
  return g;
}

// Endpoint-chart geometry: the R-side transform has already solved
// Ttrans = T_max(k, tau).  Reuse that solved value instead of running the
// forward endpoint root finder on every density/CDF evaluation.
inline BtawlGeom btawl_geometry_from_ttrans(double A, double b, double k,
                                            double Ttrans,
                                            double tau_hint = R_NaN) {
  BtawlGeom g;
  const double tau = (emc2_isfinite(tau_hint) && tau_hint > 0.0)
    ? tau_hint : btawl_tau_from_ttrans(k, Ttrans);
  if (!emc2_isfinite(tau) || !(tau > 0.0)) return g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(k) ||
      !(A >= 0.0) || !(b > 0.0) || !(b >= A) || !(k > BTAWL_K_EPS)) return g;
  g.ok = true; g.A = A; g.b = b; g.k = k; g.tau = tau;
  g.k_zero = false;
  g.t_max = Ttrans;
  g.h_max = btawl_h(Ttrans, k, tau);
  g.s_lo = btawl_tangent_time(A, g);
  return g;
}

inline double btawl_vstar(double t, double z, const BtawlGeom& g) {
  const double h = btawl_h(t, g.k, g.tau);
  if (!(h > 0.0)) return R_PosInf;
  // In the k = 0 asymptote, exp(-k * Inf) is the indeterminate 0*Inf
  // expression numerically, but the start point does not decay and E = 1.
  const double E = g.k_zero ? 1.0 : std::exp(-g.k * t);
  return (g.b - z * E) / h;
}

// Derivative of V*(t,z).  It is negative on the live/rising limb.
inline double btawl_vstar_prime(double t, double z, const BtawlGeom& g) {
  const double h = btawl_h(t, g.k, g.tau);
  if (!(h > 0.0)) return 0.0;
  const double E = std::exp(-g.k * t);
  const double hp = btawl_g(t, g.tau) - g.k * h;
  const double num = g.b - z * E;
  return (z * g.k * E * h - num * hp) / (h * h);
}

// Start point at which the tangency time equals t.  Gamma starts at one,
// rises to its maximum at u = tau, and then decreases to zero at t_max.  The
// clamp below is therefore load-bearing: Gamma is not globally monotone, and
// a Newton solve seeded from that assumption can select the wrong branch.
inline double btawl_z_t(double t, const BtawlGeom& g) {
  if (g.k_zero || !(t > 0.0) || !(t < g.t_max)) return 0.0;
  const double h = btawl_h(t, g.k, g.tau);
  const double hp = btawl_hp(t, g.k, g.tau);
  // At tangency, dV*/dt = 0 gives z E g(t) = b H'(t), where
  // E = exp(-k t) and g is the transient input (not H itself).
  const double den = std::exp(-g.k * t) * btawl_g(t, g.tau);
  if (!(hp > 0.0) || !(den > 0.0)) return 0.0;
  return std::fmin(std::fmax(g.b * hp / den, 0.0), g.b);
}

inline double btawl_tangent_time(double z, const BtawlGeom& g) {
  if (g.k_zero) return R_PosInf;
  if (!(z > 0.0)) return g.t_max;
  if (z >= g.b) return 0.0;
  double lo = 0.0, hi = g.t_max;
  for (int i = 0; i < 90; ++i) {
    const double mid = 0.5 * (lo + hi);
    const double h = btawl_h(mid, g.k, g.tau);
    const double hp = btawl_hp(mid, g.k, g.tau);
    const double d = z * std::exp(-g.k * mid) * btawl_g(mid, g.tau) -
      g.b * hp;
    if (d < 0.0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

inline double btawl_tangent_z_prime(double t, const BtawlGeom& g) {
  if (!(t > 0.0) || !(t < g.t_max)) return R_NegInf;
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

inline double btawl_normal_denom(double v, double sv, bool posdrift) {
  if (!posdrift) return 1.0;
  double d = pnorm_std(v / sv, true, false);
  return std::fmax(d, BTAWL_DENOM_FLOOR);
}

inline double btawl_surv(double w, double p1, double p2, int launch,
                         bool posdrift) {
  if (!(w > 0.0)) {
    if (launch == BTAWL_LAUNCH_LOGNORMAL) return 1.0;
    return pnorm_std(p1 / p2, true, false) /
      btawl_normal_denom(p1, p2, posdrift);
  }
  if (!emc2_isfinite(w)) return 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL)
    return pnorm_std((p1 - std::log(w)) / p2, true, false);
  const double d = btawl_normal_denom(p1, p2, posdrift);
  return pnorm_std((p1 - w) / p2, true, false) / d;
}

inline double btawl_pdf_v(double w, double p1, double p2, int launch,
                          bool posdrift) {
  if (!(w > 0.0) || !emc2_isfinite(w) || !(p2 > 0.0)) return 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL)
    return dlnorm_std(w, p1, p2, false);
  return dnormP(w, p1, p2, false) / btawl_normal_denom(p1, p2, posdrift);
}

inline double btawl_J(double x) {
  return x * pnorm_std(x, true, false) + dnormP(x);
}

inline double btawl_live_cdf(double t, double zlo, double zhi,
                             const BtawlGeom& g, double p1, double p2,
                             int launch, bool posdrift) {
  if (!(zhi > zlo)) return 0.0;
  const double h = btawl_h(t, g.k, g.tau);
  if (!(h > 0.0)) return 0.0;
  const double a = g.b / h;
  const double q = (g.k_zero ? 1.0 : std::exp(-g.k * t)) / h;
  if (!(q > 0.0) || !emc2_isfinite(q))
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift) * (zhi - zlo);
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    double sc1 = 0.0, sc2 = 0.0;
    const double whi = a - q * zhi;
    const double wlo = a - q * zlo;
    const double clo = lognormal_stoploss_nat(whi, p1, p2, sc1);
    const double chi = lognormal_stoploss_nat(wlo, p1, p2, sc2);
    const double val = (clo - chi) / q;
    if (val >= 0.0 && emc2_isfinite(val)) return val;
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift) * (zhi - zlo);
  }
  const double c = (p1 - a) / p2;
  const double m = q / p2;
  const double den = btawl_normal_denom(p1, p2, posdrift);
  if (!(m > 1e-14))
    return btawl_surv(btawl_vstar(t, 0.5 * (zlo + zhi), g), p1, p2,
                      launch, posdrift) * (zhi - zlo);
  const double val = (btawl_J(c + m * zhi) - btawl_J(c + m * zlo)) / m / den;
  return (val >= 0.0 && emc2_isfinite(val)) ? val : 0.0;
}

inline double btawl_live_pdf(double t, double zlo, double zhi,
                             const BtawlGeom& g, double p1, double p2,
                             int launch, bool posdrift) {
  if (!(zhi > zlo)) return 0.0;
  const double h = btawl_h(t, g.k, g.tau);
  const double hp = btawl_hp(t, g.k, g.tau);
  if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
  const double E = std::exp(-g.k * t);
  const double q = E / h;
  if (!(q > 0.0)) return 0.0;
  const double w_hi = btawl_vstar(t, zlo, g);
  // On the partial live/frozen seam zhi = z_t(t), the tangency identity is
  // exact: V*(t, z_t(t)) = k b / g(t).  Use that common expression so the
  // live density and frozen quadrature share the same boundary at full
  // precision.  Clamped zhi=A is an all-live endpoint, not this seam.
  const double w_lo = (zhi > 0.0 && zhi < g.A && t < g.t_max)
    ? g.k * g.b / btawl_g(t, g.tau)
    : btawl_vstar(t, zhi, g);
  const double d0 = g.b * hp / (h * h);
  const double d1 = E * (g.k * h + hp) / (h * h);
  const double c0 = d0 - d1 * (g.b / h) / q;
  const double c1 = d1 / q;
  double mass = 0.0, first = 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    const double x0 = (std::log(w_lo) - p1) / p2;
    const double x1 = (std::log(w_hi) - p1) / p2;
    mass = pnorm_std(x1, true, false) - pnorm_std(x0, true, false);
    const double M = std::exp(p1 + 0.5 * p2 * p2);
    first = M * (pnorm_std(x1 - p2, true, false) -
                 pnorm_std(x0 - p2, true, false));
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

inline double btawl_log_frozen(bool survivor, double t, const BtawlGeom& g,
                               double p1, double p2,
                               int launch, bool posdrift);

inline double btawl_frozen_cdf(double t, double zlo, double zhi,
                               const BtawlGeom& g, double p1, double p2,
                               int launch, bool posdrift) {
  if (!(zhi > zlo) || g.k_zero) return 0.0;
  // Integrate in tangency time rather than repeatedly solving for the
  // tangency point at every quadrature node.  `btawl_log_frozen` uses the
  // closed Jacobian dz/ds and the same log-Gauss--Kronrod splitter as the
  // other ballistic families, so the conversion back to natural space is
  // just one exponentiation.  The normal-launch truncation normalizer is
  // applied here because the log helper intentionally leaves it to the
  // caller, just as `btawl_log_eval` does.
  const double li = btawl_log_frozen(false, t, g, p1, p2, launch, posdrift);
  if (!(li > R_NegInf)) return 0.0;
  const double ln = (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    ? log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : 0.0;
  const double out = li - ln;
  return (emc2_isfinite(out) && out > -745.0) ? std::exp(out) : 0.0;
}

// Launch-density turnover on the descending branch of g(s) = (s/tau)e^{-s/tau}.
// The frozen quadrature is concentrated near this point in the stiff endpoint
// window, so an arithmetic midpoint is a particularly poor split.
inline double btawl_frozen_split(const BtawlGeom& g, double s_lo, double s_hi,
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
  for (int i = 0; i < 64; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (btawl_g(mid, g.tau) > target_g) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

inline double btawl_log_frozen(bool survivor, double t, const BtawlGeom& g,
                               double p1, double p2,
                               int launch, bool posdrift) {
  if (g.k_zero || !(p2 > 0.0)) return R_NegInf;
  const double s_lo = g.s_lo;
  const double s_hi = std::fmin(std::fmax(t, s_lo), g.t_max);
  if (!(s_hi > s_lo)) return R_NegInf;
  const auto lf = [&](double s) -> double {
    const double zp = btawl_tangent_z_prime(s, g);
    const double gg = btawl_g(s, g.tau);
    const double vc = (gg > 0.0) ? g.k * g.b / gg : R_PosInf;
    if (!(zp < 0.0) || !(vc > 0.0) || !emc2_isfinite(vc)) return R_NegInf;
    const double lp = launch == BTAWL_LAUNCH_LOGNORMAL
      ? pnorm_log_direct((survivor ? std::log(vc) - p1 : p1 - std::log(vc)) / p2, true)
      : (survivor
           ? (posdrift ? log_normal_cdf_positive_raw(vc, p1, p2)
                       : pnorm_log_direct((vc - p1) / p2, true))
                  : pnorm_log_direct((p1 - vc) / p2, true));
    return std::log(-zp) + lp;
  };
  const double target = launch == BTAWL_LAUNCH_LOGNORMAL
    ? std::exp(p1) : p1;
  const double split = btawl_frozen_split(g, s_lo, s_hi, target);
  return bawd_log_gl_split(lf, s_lo, s_hi, split, BAWD_GL_NODES);
}

inline double btawl_log_live(bool survivor, double t, double zlo, double zhi,
                             const BtawlGeom& g, double p1, double p2,
                             int launch, bool posdrift) {
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
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    const double pa = survivor ? log_lognormal_put(w_hi, p1, p2)
                               : log_lognormal_stoploss(w_lo, p1, p2);
    const double pb = survivor ? log_lognormal_put(w_lo, p1, p2)
                               : log_lognormal_stoploss(w_hi, p1, p2);
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
    const double lm = launch == BTAWL_LAUNCH_LOGNORMAL
      ? pnorm_log_direct((survivor ? std::log(wm) - p1 : p1 - std::log(wm)) / p2, true)
      : (survivor ? log_normal_cdf_positive_raw(wm, p1, p2)
                  : pnorm_log_direct((p1 - wm) / p2, true));
    li = std::log(zhi - zlo) + lm;
  }
  return li;
}

inline double btawl_log_eval(bool survivor, double t, const BtawlGeom& g,
                             double p1, double p2, int launch, bool posdrift) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return R_NegInf;
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    if (!(w > 0.0)) return R_NegInf;
    const double l = launch == BTAWL_LAUNCH_LOGNORMAL
      ? pnorm_log_direct((survivor ? std::log(w) - p1 : p1 - std::log(w)) / p2, true)
      : (survivor ? log_normal_cdf_positive_raw(w, p1, p2)
                  : pnorm_log_direct((p1 - w) / p2, true));
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double zcut = g.k_zero ? g.A : std::fmin(std::fmax(btawl_z_t(u, g), 0.0), g.A);
  const double ll = btawl_log_live(survivor, u, 0.0, zcut, g, p1, p2,
                                   launch, posdrift);
  const double lf = (!g.k_zero && zcut < g.A)
    ? btawl_log_frozen(survivor, u, g, p1, p2, launch, posdrift)
    : R_NegInf;
  double out = log_sum_exp(ll, lf) - std::log(g.A);
  if (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    out -= log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

inline double log_btawl_cdf_normal(double t, const BtawlGeom& g, double v,
                                   double sv, bool posdrift) {
  return btawl_log_eval(false, t, g, v, sv, BTAWL_LAUNCH_NORMAL, posdrift);
}
inline double log_btawl_cdf_logn(double t, const BtawlGeom& g, double mu,
                                 double sigma) {
  return btawl_log_eval(false, t, g, mu, sigma, BTAWL_LAUNCH_LOGNORMAL, false);
}
inline double log_btawl_surv_normal(double t, const BtawlGeom& g, double v,
                                    double sv, bool posdrift) {
  return btawl_log_eval(true, t, g, v, sv, BTAWL_LAUNCH_NORMAL, posdrift);
}
inline double log_btawl_surv_logn(double t, const BtawlGeom& g, double mu,
                                  double sigma) {
  return btawl_log_eval(true, t, g, mu, sigma, BTAWL_LAUNCH_LOGNORMAL, false);
}

inline double btawl_cdf(double t, double A, double b, double p1, double p2,
                        double k, double tau, int launch, bool posdrift) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return 0.0;
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (u == R_PosInf) u = g.k_zero ? R_PosInf : g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    return std::fmin(std::fmax(btawl_surv(w, p1, p2, launch, posdrift), 0.0), 1.0);
  }
  // Starts below zcut are still on the rising/live limb.  Higher starts have
  // already reached their tangency wall and contribute the frozen mass.
  double zcut = g.k_zero ? g.A : btawl_z_t(u, g);
  zcut = std::fmin(std::fmax(zcut, 0.0), g.A);
  const double live = (u == R_PosInf && !g.k_zero) ? 0.0 :
    btawl_live_cdf(u, 0.0, zcut, g, p1, p2, launch, posdrift);
  const double frozen = (!g.k_zero && zcut < g.A)
    ? btawl_frozen_cdf(u, zcut, g.A, g, p1, p2, launch, posdrift) : 0.0;
  const double out = (live + frozen) / g.A;
  return std::fmin(std::fmax(out, 0.0), 1.0);
}

inline double btawl_pdf(double t, double A, double b, double p1, double p2,
                        double k, double tau, int launch, bool posdrift) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  if (!g.k_zero && !(t < g.t_max)) return 0.0;
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double hp = btawl_hp(t, g.k, g.tau);
    if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
    const double w = g.b / h;
    return btawl_pdf_v(w, p1, p2, launch, posdrift) * g.b * hp / (h * h);
  }
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  return btawl_live_pdf(t, 0.0, zcut, g, p1, p2, launch, posdrift);
}

inline double btawl_cdf_from_geom(double t, const BtawlGeom& g, double p1,
                                  double p2, int launch, bool posdrift) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return 0.0;
  double u = t;
  if (!g.k_zero && u >= g.t_max) u = g.t_max;
  if (u == R_PosInf) u = g.k_zero ? R_PosInf : g.t_max;
  if (g.A <= BTAWL_A_EPS) {
    const double w = (u == R_PosInf) ? g.b / g.h_max : btawl_vstar(u, 0.0, g);
    return std::fmin(std::fmax(btawl_surv(w, p1, p2, launch, posdrift), 0.0), 1.0);
  }
  double zcut = g.k_zero ? g.A : btawl_z_t(u, g);
  zcut = std::fmin(std::fmax(zcut, 0.0), g.A);
  const double live = (u == R_PosInf && !g.k_zero) ? 0.0 :
    btawl_live_cdf(u, 0.0, zcut, g, p1, p2, launch, posdrift);
  const double frozen = (!g.k_zero && zcut < g.A)
    ? btawl_frozen_cdf(u, zcut, g.A, g, p1, p2, launch, posdrift) : 0.0;
  return std::fmin(std::fmax((live + frozen) / g.A, 0.0), 1.0);
}

inline bool btawl_natural_cdf_from_geom(double t, const BtawlGeom& g,
                                        double p1, double p2, int launch,
                                        bool posdrift, double& cdf) {
  cdf = btawl_cdf_from_geom(t, g, p1, p2, launch, posdrift);
  // BA_ACCEPT_RAW deliberately rejects near-one CDFs: log1p(-cdf) has lost
  // the survivor tail there.  A rounded zero is harmless and is accepted.
  return emc2_isfinite(cdf) && cdf >= 0.0 && cdf < 1.0 - 1e-8;
}

inline double btawl_pdf_from_geom(double t, const BtawlGeom& g, double p1,
                                  double p2, int launch, bool posdrift) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  if (!g.k_zero && !(t < g.t_max)) return 0.0;
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double hp = btawl_hp(t, g.k, g.tau);
    if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
    const double w = g.b / h;
    return btawl_pdf_v(w, p1, p2, launch, posdrift) * g.b * hp / (h * h);
  }
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  return btawl_live_pdf(t, 0.0, zcut, g, p1, p2, launch, posdrift);
}

inline double btawl_log_launch_pdf(double w, double p1, double p2, int launch,
                                   bool posdrift) {
  if (!(w > 0.0) || !emc2_isfinite(w) || !(p2 > 0.0)) return R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL)
    return dlnorm_std(w, p1, p2, true);
  const double ld = dnormP((w - p1) / p2, 0.0, 1.0, true) - std::log(p2);
  return posdrift ? ld - log_positive_normalizer(p1, p2, true,
                                                  BTAWL_DENOM_FLOOR) : ld;
}

inline bool btawl_natural_pdf_accepted(double t, const BtawlGeom& g,
                                        double p1, double p2, int launch,
                                        bool posdrift, double p_nat) {
  (void)posdrift;
  if (!(p_nat > 0.0) || !emc2_isfinite(p_nat)) return false;
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double w = (h > 0.0) ? g.b / h : R_PosInf;
    const double x = launch == BTAWL_LAUNCH_LOGNORMAL
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
  const double x0 = launch == BTAWL_LAUNCH_LOGNORMAL
    ? (std::log(w_lo) - p1) / p2 : (w_lo - p1) / p2;
  const double x1 = launch == BTAWL_LAUNCH_LOGNORMAL
    ? (std::log(w_hi) - p1) / p2 : (w_hi - p1) / p2;
  if (!emc2_isfinite(x0) || !emc2_isfinite(x1) ||
      std::fabs(x0) > BAWL_NATURAL_Z_MAX ||
      std::fabs(x1) > BAWL_NATURAL_Z_MAX)
    return false;
  // The first-moment primitive uses the shifted interval as well.  Requiring
  // both differences to remain resolved is the BA_ACCEPT_RAW cancellation
  // contract; a merely positive natural density is not enough.
  const double shift = (launch == BTAWL_LAUNCH_LOGNORMAL) ? p2 : 0.0;
  if (std::fabs(x0 - shift) > BAWL_NATURAL_Z_MAX ||
      std::fabs(x1 - shift) > BAWL_NATURAL_Z_MAX)
    return false;
  const double mass = pnorm_std(x1, true, false) - pnorm_std(x0, true, false);
  return mass > BAWL_NATURAL_REL_TOL *
    std::fmax(1.0, std::fabs(pnorm_std(x1, true, false)) +
                     std::fabs(pnorm_std(x0, true, false)));
}

inline double btawl_log_pdf_from_geom(double t, const BtawlGeom& g, double p1,
                                      double p2, int launch, bool posdrift) {
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return R_NegInf;
  if (!g.k_zero && !(t < g.t_max)) return R_NegInf;
  const double p_nat = btawl_pdf_from_geom(t, g, p1, p2, launch, posdrift);
  if (p_nat > 1e-280 &&
      btawl_natural_pdf_accepted(t, g, p1, p2, launch, posdrift, p_nat))
    return std::log(p_nat);
  if (g.A <= BTAWL_A_EPS) {
    const double h = btawl_h(t, g.k, g.tau);
    const double hp = btawl_hp(t, g.k, g.tau);
    if (!(h > 0.0) || !(hp > 0.0)) return R_NegInf;
    const double w = g.b / h;
    return btawl_log_launch_pdf(w, p1, p2, launch, posdrift) +
      std::log(g.b) + std::log(hp) - 2.0 * std::log(h);
  }
  const double zcut = g.k_zero ? g.A :
    std::fmin(std::fmax(btawl_z_t(t, g), 0.0), g.A);
  if (!(zcut > 0.0)) return R_NegInf;
  const auto lf = [&](double z) -> double {
    const double w = btawl_vstar(t, z, g);
    const double vp = btawl_vstar_prime(t, z, g);
    return (vp < 0.0) ? btawl_log_launch_pdf(w, p1, p2, launch, posdrift) +
      std::log(-vp) : R_NegInf;
  };
  const double li = bawd_log_gl_split(lf, 0.0, zcut, 0.5 * zcut,
                                      BAWD_GL_NODES);
  return (li > R_NegInf) ? li - std::log(g.A) : R_NegInf;
}

inline double btawl_log_cdf(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  return launch == BTAWL_LAUNCH_LOGNORMAL
    ? log_btawl_cdf_logn(t, g, p1, p2)
    : log_btawl_cdf_normal(t, g, p1, p2, posdrift);
}
inline double btawl_log_surv(double t, double A, double b, double p1, double p2,
                             double k, double tau, int launch, bool posdrift) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  return launch == BTAWL_LAUNCH_LOGNORMAL
    ? log_btawl_surv_logn(t, g, p1, p2)
    : log_btawl_surv_normal(t, g, p1, p2, posdrift);
}
inline double btawl_log_pdf(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift) {
  const BtawlGeom g = btawl_geometry(A, b, k, tau);
  return btawl_log_pdf_from_geom(t, g, p1, p2, launch, posdrift);
}
inline double btawl_cdf_log(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift) {
  return btawl_log_cdf(t, A, b, p1, p2, k, tau, launch, posdrift);
}
inline double btawl_pdf_log(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift) {
  return btawl_log_pdf(t, A, b, p1, p2, k, tau, launch, posdrift);
}

inline double btawl_cdf_chart(double t, double A, double b, double p1,
                              double p2, double k, double clear, int launch,
                              bool posdrift, bool endpoint_chart,
                              double tau_hint = R_NaN) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  return btawl_cdf_from_geom(t, g, p1, p2, launch, posdrift);
}
inline double btawl_pdf_chart(double t, double A, double b, double p1,
                              double p2, double k, double clear, int launch,
                              bool posdrift, bool endpoint_chart,
                              double tau_hint = R_NaN) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  return btawl_pdf_from_geom(t, g, p1, p2, launch, posdrift);
}
inline double btawl_log_surv_chart(double t, double A, double b, double p1,
                                   double p2, double k, double clear,
                                   int launch, bool posdrift,
                                   bool endpoint_chart,
                                   double tau_hint = R_NaN) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  return launch == BTAWL_LAUNCH_LOGNORMAL
    ? log_btawl_surv_logn(t, g, p1, p2)
    : log_btawl_surv_normal(t, g, p1, p2, posdrift);
}
inline double btawl_log_pdf_chart(double t, double A, double b, double p1,
                                  double p2, double k, double clear,
                                  int launch, bool posdrift,
                                  bool endpoint_chart,
                                  double tau_hint = R_NaN) {
  const BtawlGeom g = endpoint_chart
    ? btawl_geometry_from_ttrans(A, b, k, clear, tau_hint)
    : btawl_geometry(A, b, k, clear);
  return btawl_log_pdf_from_geom(t, g, p1, p2, launch, posdrift);
}

// ---------------------------------------------------------------------------
// Shared-strength sustained + transient extension.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Shared-strength sustained + transient extension.
// ---------------------------------------------------------------------------

inline double btawl_hs(double t, double k, double tau) {
  if (!(tau > 0.0) || t <= 0.0) return 0.0;
  if (t == R_PosInf) return k <= BTAWL_K_EPS ? R_PosInf : 1.0 / k;
  if (k <= BTAWL_K_EPS)
    return t - tau * (-std::expm1(-t / tau));
  const double d = k - 1.0 / tau;
  const double y = d * t;
  double ratio;
  if (std::fabs(y) < 1e-4) {
    ratio = t * (1.0 + y * (0.5 + y * (1.0 / 6.0 + y / 24.0)));
  } else {
    return (1.0 - std::exp(-k * t)) / k - (std::exp(-t / tau) - std::exp(-k * t)) / d;
  }
  return (1.0 - std::exp(-k * t)) / k - std::exp(-k * t) * ratio;
}

inline double btawl_hs_p(double t, double k, double tau) {
  return (1.0 - std::exp(-t / tau)) - k * btawl_hs(t, k, tau);
}

inline double btawl_sustained_cdf(double t, double A, double b, double p1, double p2,
                                  double k, double tau_s, int launch, bool posdrift) {
  if (!(p2 > 0.0)) return 0.0;
  if (t == R_PosInf) {
    const double w = (k <= BTAWL_K_EPS) ? 0.0 : (k * b);
    return btawl_surv(w, p1, p2, launch, posdrift);
  }
  const double h = btawl_hs(t, k, tau_s);
  if (!(h > 0.0)) return 0.0;
  const double E = std::exp(-k * t);
  const double q = E / h;
  if (!(q > 0.0) || !emc2_isfinite(q)) return 0.0;
  
  if (A <= BTAWL_A_EPS) return btawl_surv(b / h, p1, p2, launch, posdrift);
  
  const double a = b / h;
  const double whi = a;
  const double wlo = a - q * A;
  if (!(wlo > 0.0) || !(whi >= wlo)) return 0.0;
  
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    double s_lo = 0.0, s_hi = 0.0;
    const double val = (lognormal_stoploss_nat(wlo, p1, p2, s_lo) -
                        lognormal_stoploss_nat(whi, p1, p2, s_hi));
    return (val >= 0.0 && emc2_isfinite(val)) ? (val / q) / A : 0.0;
  }
  
  const double c = (p1 - a) / p2;
  const double m = q / p2;
  if (!(m > 1e-14)) return btawl_surv(0.5 * (wlo + whi), p1, p2, launch, posdrift);
  const double den = btawl_normal_denom(p1, p2, posdrift);
  const double val = (btawl_J(c + m * A) - btawl_J(c)) / m / den;
  return (val >= 0.0 && emc2_isfinite(val)) ? val / A : 0.0;
}

inline double btawl_sustained_pdf(double t, double A, double b, double p1, double p2,
                                  double k, double tau_s, int launch, bool posdrift) {
  if (!(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  const double h = btawl_hs(t, k, tau_s);
  const double hp = btawl_hs_p(t, k, tau_s);
  if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
  
  if (A <= BTAWL_A_EPS) {
     const double v = b / h;
     return btawl_pdf_v(v, p1, p2, launch, posdrift) * (b * hp / (h * h));
  }
  
  const double E = std::exp(-k * t);
  const double q = E / h;
  const double w_hi = b / h;
  const double w_lo = (b - A * E) / h;
  
  const double d0 = b * hp / (h * h);
  const double d1 = E * (k * h + hp) / (h * h);
  const double c0 = d0 - d1 * (b / h) / q;
  const double c1 = d1 / q;
  double mass = 0.0, first = 0.0;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    const double x0 = (std::log(w_lo) - p1) / p2;
    const double x1 = (std::log(w_hi) - p1) / p2;
    mass = pnorm_std(x1, true, false) - pnorm_std(x0, true, false);
    const double M = std::exp(p1 + 0.5 * p2 * p2);
    first = M * (pnorm_std(x1 - p2, true, false) - pnorm_std(x0 - p2, true, false));
  } else {
    const double x0 = (w_lo - p1) / p2;
    const double x1 = (w_hi - p1) / p2;
    const double den = btawl_normal_denom(p1, p2, posdrift);
    mass = (pnorm_std(x1, true, false) - pnorm_std(x0, true, false)) / den;
    first = (p1 * (pnorm_std(x1, true, false) - pnorm_std(x0, true, false)) +
             p2 * (dnormP(x0) - dnormP(x1))) / den;
  }
  const double out = (c0 * mass + c1 * first) / (A * q);
  return (out > 0.0 && emc2_isfinite(out)) ? out : 0.0;
}

inline double btawl_mix_cdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift) {
  if (pi <= 1e-14) return btawl_cdf(t, A, b, p1, p2, k, tau_t, launch, posdrift);
  if (pi >= 1.0 - 1e-14) return btawl_sustained_cdf(t, A, b, p1, p2, k, tau_s, launch, posdrift);
  
  double p1_S = p1, p2_S = p2;
  double p1_T = p1, p2_T = p2;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    p1_S += std::log(pi);
    p1_T += std::log(1.0 - pi);
  } else {
    p1_S *= pi;      p2_S *= pi;
    p1_T *= (1.0 - pi); p2_T *= (1.0 - pi);
  }
  
  const double F_T = btawl_cdf(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift);
  const double F_S = btawl_sustained_cdf(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift);
  return 1.0 - (1.0 - F_T) * (1.0 - F_S);
}

inline double btawl_mix_pdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift) {
  if (pi <= 1e-14) return btawl_pdf(t, A, b, p1, p2, k, tau_t, launch, posdrift);
  if (pi >= 1.0 - 1e-14) return btawl_sustained_pdf(t, A, b, p1, p2, k, tau_s, launch, posdrift);
  
  double p1_S = p1, p2_S = p2;
  double p1_T = p1, p2_T = p2;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    p1_S += std::log(pi);
    p1_T += std::log(1.0 - pi);
  } else {
    p1_S *= pi;      p2_S *= pi;
    p1_T *= (1.0 - pi); p2_T *= (1.0 - pi);
  }
  
  const double F_T = btawl_cdf(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift);
  const double F_S = btawl_sustained_cdf(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift);
  const double f_T = btawl_pdf(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift);
  const double f_S = btawl_sustained_pdf(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift);
  
  return f_T * (1.0 - F_S) + f_S * (1.0 - F_T);
}

inline double btawl_mix_cdf_log(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, double tau_t, double pi,
                                int launch, bool posdrift) {
  const double p = btawl_mix_cdf(t, A, b, p1, p2, k, tau_s, tau_t, pi, launch, posdrift);
  return p > 0.0 ? std::log(p) : R_NegInf;
}

inline double btawl_mix_pdf_log(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, double tau_t, double pi,
                                int launch, bool posdrift) {
  const double p = btawl_mix_pdf(t, A, b, p1, p2, k, tau_s, tau_t, pi, launch, posdrift);
  return p > 0.0 ? std::log(p) : R_NegInf;
}

inline double btawl_sustained_log_surv(double t, double A, double b, double p1, double p2,
                                       double k, double tau_s, int launch, bool posdrift) {
  if (!emc2_isfinite(b) || !(p2 > 0.0) || !(t > 0.0)) return R_NegInf;
  if (t == R_PosInf) {
    const double w = (k <= BTAWL_K_EPS) ? 0.0 : (k * b);
    if (!(w > 0.0)) return R_NegInf;
    if (launch == BTAWL_LAUNCH_LOGNORMAL)
      return pnorm_log_direct((std::log(w) - p1) / p2, true);
    const double l = log_normal_cdf_positive_raw(w, p1, p2);
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double h = btawl_hs(t, k, tau_s);
  if (!(h > 0.0)) return 0.0;
  const double E = std::exp(-k * t);
  const double q = E / h;
  if (!(q > 0.0) || !emc2_isfinite(q)) {
    const double w = b / h;
    if (launch == BTAWL_LAUNCH_LOGNORMAL)
      return pnorm_log_direct((std::log(w) - p1) / p2, true);
    const double l = log_normal_cdf_positive_raw(w, p1, p2);
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  if (A <= BTAWL_A_EPS) {
    const double w = b / h;
    if (!(w > 0.0)) return R_NegInf;
    const double l = launch == BTAWL_LAUNCH_LOGNORMAL
      ? pnorm_log_direct((std::log(w) - p1) / p2, true)
      : log_normal_cdf_positive_raw(w, p1, p2);
    return (launch == BTAWL_LAUNCH_NORMAL && posdrift)
      ? l - log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR) : l;
  }
  const double a = b / h;
  const double whi = a;
  const double wlo = a - q * A;
  if (!(wlo > 0.0) || !(whi >= wlo)) return R_NegInf;
  
  double li = R_NegInf;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    const double pa = log_lognormal_put(whi, p1, p2);
    const double pb = log_lognormal_put(wlo, p1, p2);
    if (pa - pb > 1e-10) li = log_diff_exp(pa, pb) - std::log(q);
  } else {
    const double lo = (wlo - p1) / p2;
    const double hi = (whi - p1) / p2;
    const double x = log_normal_phi_integral_positive_raw(lo, hi, p1, p2, posdrift);
    if (x > R_NegInf) li = x + std::log(p2) - std::log(q);
  }
  if (!(li > R_NegInf)) {
    const double wm = 0.5 * (whi + wlo);
    const double lm = launch == BTAWL_LAUNCH_LOGNORMAL
      ? pnorm_log_direct((std::log(wm) - p1) / p2, true)
      : log_normal_cdf_positive_raw(wm, p1, p2);
    li = std::log(A) + lm;
  }
  double out = li - std::log(A);
  if (launch == BTAWL_LAUNCH_NORMAL && posdrift)
    out -= log_positive_normalizer(p1, p2, true, BTAWL_DENOM_FLOOR);
  return ISNAN(out) ? R_NegInf : std::fmin(out, 0.0);
}

inline double btawl_mix_log_surv(double t, double A, double b, double p1,
                                 double p2, double k, double tau_s,
                                 double tau_t, double pi, int launch,
                                 bool posdrift) {
  if (pi <= 1e-14)
    return btawl_log_surv(t, A, b, p1, p2, k, tau_t, launch, posdrift);
  if (pi >= 1.0 - 1e-14)
    return btawl_sustained_log_surv(t, A, b, p1, p2, k, tau_s, launch, posdrift);
  
  double p1_S = p1, p2_S = p2;
  double p1_T = p1, p2_T = p2;
  if (launch == BTAWL_LAUNCH_LOGNORMAL) {
    p1_S += std::log(pi);
    p1_T += std::log(1.0 - pi);
  } else {
    p1_S *= pi;      p2_S *= pi;
    p1_T *= (1.0 - pi); p2_T *= (1.0 - pi);
  }
  
  const double ls_T = btawl_log_surv(t, A, b, p1_T, p2_T, k, tau_t, launch, posdrift);
  const double ls_S = btawl_sustained_log_surv(t, A, b, p1_S, p2_S, k, tau_s, launch, posdrift);
  
  if (!(ls_T > R_NegInf) || !(ls_S > R_NegInf)) return R_NegInf;
  return ls_T + ls_S;
}
// [[Rcpp::export]]
NumericVector dbtawl(NumericVector t, NumericVector A, NumericVector b,
                     NumericVector p1, NumericVector p2, NumericVector k,
                     NumericVector tau, int launch = 1,
                     bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_pdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i),
                                    pick(p2,i), pick(k,i), pick(tau,i), launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pbtawl(NumericVector t, NumericVector A, NumericVector b,
                     NumericVector p1, NumericVector p2, NumericVector k,
                     NumericVector tau, int launch = 1,
                     bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_cdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i),
                                    pick(p2,i), pick(k,i), pick(tau,i), launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_log_surv_vec(NumericVector t, NumericVector A,
                                 NumericVector b, NumericVector p1,
                                 NumericVector p2, NumericVector k,
                                 NumericVector tau, int launch = 1,
                                 bool posdrift = true) {
  const int n = std::max({t.size(), A.size(), b.size(), p1.size(), p2.size(),
                          k.size(), tau.size()});
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i)
    out[i] = btawl_log_surv(pick(t, i), pick(A, i), pick(b, i), pick(p1, i),
                            pick(p2, i), pick(k, i), pick(tau, i), launch,
                            posdrift);
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_mix_log_surv_vec(NumericVector t, NumericVector A,
                                     NumericVector b, NumericVector p1,
                                     NumericVector p2, NumericVector k,
                                     NumericVector tau_s, NumericVector tau_t,
                                     NumericVector pi, int launch = 1,
                                     bool posdrift = true) {
  const int n = std::max({t.size(), A.size(), b.size(), p1.size(), p2.size(),
                          k.size(), tau_s.size(), tau_t.size(), pi.size()});
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i)
    out[i] = btawl_mix_log_surv(pick(t, i), pick(A, i), pick(b, i), pick(p1, i),
                                pick(p2, i), pick(k, i), pick(tau_s, i),
                                pick(tau_t, i), pick(pi, i), launch, posdrift);
  return out;
}

// [[Rcpp::export]]
NumericVector dbtawlmix(NumericVector t, NumericVector A, NumericVector b,
                        NumericVector p1, NumericVector p2, NumericVector k,
                        NumericVector tau_s, NumericVector tau_t,
                        NumericVector pi, int launch = 1,
                        bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_mix_pdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i), pick(p2,i),
                                        pick(k,i), pick(tau_s,i), pick(tau_t,i), pick(pi,i),
                                        launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pbtawlmix(NumericVector t, NumericVector A, NumericVector b,
                        NumericVector p1, NumericVector p2, NumericVector k,
                        NumericVector tau_s, NumericVector tau_t,
                        NumericVector pi, int launch = 1,
                        bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_mix_cdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i), pick(p2,i),
                                        pick(k,i), pick(tau_s,i), pick(tau_t,i), pick(pi,i),
                                        launch, posdrift);
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

#endif

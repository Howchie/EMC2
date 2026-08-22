#ifndef EMC2_MODEL_BTAWL_H
#define EMC2_MODEL_BTAWL_H

// BTAwL: ballistic transient accumulator with leak.
//
//   dX/du = V (u/tau) exp(-u/tau) - k X,  X(0) = z, z ~ U(0,A)
//
// The transient response is affine in the launch point and the required
// launch V*(u,z) is therefore affine in z on the live part of the CDF.  The
// live start-point average is evaluated in closed form (normal and lognormal
// launches).  The only non-elementary piece is the already-saturated start
// region; it is a one-dimensional adaptive integral over z, with no grid or
// PDE solve.
//
// This header is included by model_LBA.h, after the shared normal/lognormal
// primitives and the GSL helpers have been declared.

#include <cmath>
#include <functional>
#include <algorithm>

constexpr int BTAWL_LAUNCH_NORMAL = 0;
constexpr int BTAWL_LAUNCH_LOGNORMAL = 1;
constexpr double BTAWL_K_EPS = 1e-10;
constexpr double BTAWL_A_EPS = 1e-10;
constexpr double BTAWL_DENOM_FLOOR = BAWL_DENOM_FLOOR;

struct BtawlGeom {
  bool ok = false;
  bool k_zero = false;
  double A = 0.0, b = 0.0, k = 0.0, tau = 0.0;
  double t_max = R_PosInf;
  double h_max = 0.0;
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
  double hi = std::fmax(tau, 1.0 / k);
  if (!(hi > 0.0) || !emc2_isfinite(hi)) return R_PosInf;
  hi *= 2.0;
  for (int i = 0; i < 80 && btawl_hp(hi, k, tau) > 0.0; ++i)
    hi *= 2.0;
  if (!(btawl_hp(hi, k, tau) <= 0.0) || !emc2_isfinite(hi)) return R_PosInf;
  double lo = 0.0;
  for (int i = 0; i < 100; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (btawl_hp(mid, k, tau) > 0.0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
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

// Start point at which the tangency time equals t.  The value is monotone
// decreasing from b at onset to zero at t_max.
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
  const double w_lo = btawl_vstar(t, zhi, g);
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

struct BtawlFrozenCtx {
  const BtawlGeom* g;
  double p1, p2;
  int launch;
  bool posdrift;
};

inline double btawl_frozen_integrand(double z, void* ptr) {
  const BtawlFrozenCtx& c = *static_cast<BtawlFrozenCtx*>(ptr);
  const double tz = btawl_tangent_time(z, *c.g);
  const double v = (tz > 0.0 && tz < R_PosInf)
    ? c.g->k * c.g->b / btawl_g(tz, c.g->tau) : R_PosInf;
  return btawl_surv(v, c.p1, c.p2, c.launch, c.posdrift);
}

inline double btawl_frozen_cdf(double zlo, double zhi, const BtawlGeom& g, double p1,
                               double p2, int launch, bool posdrift) {
  if (!(zhi > zlo) || g.k_zero) return 0.0;
  BtawlFrozenCtx ctx{&g, p1, p2, launch, posdrift};
  gsl_function F;
  F.function = &btawl_frozen_integrand;
  F.params = &ctx;
  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(ws);
  double out = 0.0, err = 0.0;
  gsl_error_handler_t* old = gsl_set_error_handler_off();
  const int status = gsl_integration_qags(&F, zlo, zhi, 1e-9, 1e-7, 256,
                                          w, &out, &err);
  gsl_set_error_handler(old);
  if (status != GSL_SUCCESS || !emc2_isfinite(out)) return 0.0;
  return std::fmax(0.0, out);
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
    ? btawl_frozen_cdf(zcut, g.A, g, p1, p2, launch, posdrift) : 0.0;
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

inline double btawl_cdf_log(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift) {
  const double p = btawl_cdf(t, A, b, p1, p2, k, tau, launch, posdrift);
  return p > 0.0 ? std::log(p) : R_NegInf;
}
inline double btawl_pdf_log(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift) {
  const double p = btawl_pdf(t, A, b, p1, p2, k, tau, launch, posdrift);
  return p > 0.0 ? std::log(p) : R_NegInf;
}

// ---------------------------------------------------------------------------
// Shared-strength sustained + transient extension.
// ---------------------------------------------------------------------------

struct BtawlMixGeom {
  bool ok = false;
  double A = 0.0, b = 0.0, k = 0.0, tau_s = 0.0, tau_t = 0.0, pi = 0.0;
  double t_peak = R_PosInf;
  double h_peak = 0.0;
  double h_inf = R_PosInf;
  bool finite_peak = false;
};

inline double btawl_hs(double t, double k, double tau) {
  if (!(tau > 0.0) || t <= 0.0) return 0.0;
  if (t == R_PosInf) return k <= BTAWL_K_EPS ? R_PosInf : 1.0 / k;
  if (k <= BTAWL_K_EPS)
    return t - tau * (-std::expm1(-t / tau));
  const double d = k - 1.0 / tau;
  const double y = d * t;
  double ratio;
  if (std::fabs(y) < 1e-4) {
    // (exp(y) - 1) / d, written as t * expm1(y) / y.
    ratio = t * (1.0 + y * (0.5 + y * (1.0 / 6.0 + y / 24.0)));
  } else {
    // The direct two-exponential form avoids overflow when d*t is large and
    // positive, while the series branch above handles cancellation at d = 0.
    return (1.0 - std::exp(-k * t)) / k -
      (std::exp(-t / tau) - std::exp(-k * t)) / d;
  }
  return (1.0 - std::exp(-k * t)) / k - std::exp(-k * t) * ratio;
}

inline double btawl_hmix(double t, double k, double tau_s, double tau_t,
                         double pi) {
  return pi * btawl_hs(t, k, tau_s) + (1.0 - pi) * btawl_h(t, k, tau_t);
}

inline double btawl_hmix_p(double t, double k, double tau_s, double tau_t,
                           double pi) {
  return pi * (1.0 - std::exp(-t / tau_s)) +
    (1.0 - pi) * btawl_g(t, tau_t) - k * btawl_hmix(t, k, tau_s, tau_t, pi);
}

inline double btawl_mix_peak(double k, double tau_s, double tau_t, double pi) {
  if (!(k > BTAWL_K_EPS) || !(tau_s > 0.0) || !(tau_t > 0.0) || !(pi >= 0.0))
    return R_PosInf;
  double hi = 2.0 * std::fmax(std::fmax(tau_s, tau_t), 1.0 / k);
  if (!(btawl_hmix_p(hi, k, tau_s, tau_t, pi) < 0.0)) return R_PosInf;
  double lo = 0.0;
  for (int i = 0; i < 100; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (btawl_hmix_p(mid, k, tau_s, tau_t, pi) > 0.0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

inline BtawlMixGeom btawl_mix_geometry(double A, double b, double k,
                                       double tau_s, double tau_t, double pi) {
  BtawlMixGeom g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(k) ||
      !emc2_isfinite(tau_s) || !emc2_isfinite(tau_t) || !emc2_isfinite(pi) ||
      !(A >= 0.0) || !(b > 0.0) || !(b >= A) || !(k >= 0.0) ||
      !(tau_s > 0.0) || !(tau_t > 0.0) || !(pi >= 0.0) || !(pi <= 1.0)) return g;
  g.ok = true; g.A = A; g.b = b; g.k = k; g.tau_s = tau_s; g.tau_t = tau_t; g.pi = pi;
  g.h_inf = k <= BTAWL_K_EPS ? (pi > 0.0 ? R_PosInf : tau_t) : pi / k;
  g.t_peak = btawl_mix_peak(k, tau_s, tau_t, pi);
  g.finite_peak = emc2_isfinite(g.t_peak);
  g.h_peak = g.finite_peak ? btawl_hmix(g.t_peak, k, tau_s, tau_t, pi) :
    (g.h_inf == R_PosInf ? 0.0 : g.h_inf);
  return g;
}

inline double btawl_mix_vstar(double t, double z, const BtawlMixGeom& g) {
  const double h = btawl_hmix(t, g.k, g.tau_s, g.tau_t, g.pi);
  if (!(h > 0.0)) return R_PosInf;
  return (g.b - z * std::exp(-g.k * t)) / h;
}

inline double btawl_mix_d(double t, double z, const BtawlMixGeom& g) {
  const double h = btawl_hmix(t, g.k, g.tau_s, g.tau_t, g.pi);
  const double hp = btawl_hmix_p(t, g.k, g.tau_s, g.tau_t, g.pi);
  const double input = g.pi * (1.0 - std::exp(-t / g.tau_s)) +
    (1.0 - g.pi) * btawl_g(t, g.tau_t);
  return z * std::exp(-g.k * t) * input - g.b * hp;
}

// Minimum of V*(u,z) over the available history.  The derivative has at most
// two sign changes for the positive sustained/transient mixture.  We bracket
// all negative-to-positive changes on a logarithmic time scale, so an
// overshooting transient and the later sustained asymptote are both retained.
inline double btawl_mix_min_v(double t, double z, const BtawlMixGeom& g) {
  const double scale = std::fmax(std::fmax(g.tau_s, g.tau_t),
                                 g.k > BTAWL_K_EPS ? 1.0 / g.k : g.tau_t);
  double hi = t;
  if (!(hi > 0.0) || hi == R_PosInf)
    hi = g.finite_peak ? std::fmax(g.t_peak * 16.0, scale * 16.0) : scale * 1e4;
  hi = std::fmax(hi, scale * 1e-8);
  const double lo = hi * 1e-10;
  double best = btawl_mix_vstar(hi, z, g);
  const double h_inf = g.h_inf;
  if (t == R_PosInf && h_inf > 0.0 && emc2_isfinite(h_inf))
    best = std::fmin(best, g.b / h_inf);
  else if (t == R_PosInf && h_inf == R_PosInf)
    best = 0.0;
  double prev_t = lo;
  double prev_d = btawl_mix_d(prev_t, z, g);
  const int n = 96;
  for (int i = 1; i <= n; ++i) {
    const double frac = static_cast<double>(i) / n;
    const double cur_t = lo * std::pow(hi / lo, frac);
    const double cur_d = btawl_mix_d(cur_t, z, g);
    if (prev_d < 0.0 && cur_d >= 0.0) {
      double a = prev_t, b = cur_t;
      for (int it = 0; it < 70; ++it) {
        const double m = 0.5 * (a + b);
        if (btawl_mix_d(m, z, g) < 0.0) a = m; else b = m;
      }
      best = std::fmin(best, btawl_mix_vstar(0.5 * (a + b), z, g));
    }
    prev_t = cur_t; prev_d = cur_d;
    best = std::fmin(best, btawl_mix_vstar(cur_t, z, g));
  }
  if (t < R_PosInf) best = std::fmin(best, btawl_mix_vstar(t, z, g));
  return best;
}

struct BtawlMixIntCtx { const BtawlMixGeom* g; double t, p1, p2; int launch; bool posdrift; };
inline double btawl_mix_cdf_integrand(double z, void* ptr) {
  const BtawlMixIntCtx& c = *static_cast<BtawlMixIntCtx*>(ptr);
  return btawl_surv(btawl_mix_min_v(c.t, z, *c.g), c.p1, c.p2,
                    c.launch, c.posdrift);
}

inline double btawl_mix_pdf_integrand(double z, void* ptr) {
  const BtawlMixIntCtx& c = *static_cast<BtawlMixIntCtx*>(ptr);
  const BtawlMixGeom& g = *c.g;
  if (!(c.t > 0.0) || c.t == R_PosInf) return 0.0;
  const double d = btawl_mix_d(c.t, z, g);
  if (!(d < 0.0)) return 0.0;
  // A descending limb contributes density only once it has become the
  // running minimum of V*(u,z).  This matters when a strong transient
  // overshoots the later sustained asymptote: the interval between the two
  // minima is frozen, not live.
  const double t_prev = c.t * (1.0 - 1e-8);
  if (t_prev > 0.0 && btawl_mix_vstar(c.t, z, g) >
      btawl_mix_min_v(t_prev, z, g) * (1.0 + 1e-7)) return 0.0;
  const double h = btawl_hmix(c.t, g.k, g.tau_s, g.tau_t, g.pi);
  const double hp = btawl_hmix_p(c.t, g.k, g.tau_s, g.tau_t, g.pi);
  const double E = std::exp(-g.k * c.t);
  const double vp = (z * g.k * E * h - (g.b - z * E) * hp) / (h * h);
  return btawl_pdf_v(btawl_mix_vstar(c.t, z, g), c.p1, c.p2,
                     c.launch, c.posdrift) * (-vp);
}

inline double btawl_mix_integrate(double (*fn)(double, void*), double A,
                                  BtawlMixIntCtx& ctx) {
  if (!(A > 0.0)) return 0.0;
  gsl_function F; F.function = fn; F.params = &ctx;
  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(ws);
  double out = 0.0, err = 0.0;
  gsl_error_handler_t* old = gsl_set_error_handler_off();
  const int status = gsl_integration_qags(&F, 0.0, A, 1e-8, 1e-6, 256, w, &out, &err);
  gsl_set_error_handler(old);
  return (status == GSL_SUCCESS && emc2_isfinite(out)) ? std::fmax(0.0, out) : 0.0;
}

inline double btawl_mix_cdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift) {
  // Exact nesting: route pi = 0 through the original transient-only kernel,
  // preserving every last bit of its closed-form evaluation.
  if (pi <= 1e-14)
    return btawl_cdf(t, A, b, p1, p2, k, tau_t, launch, posdrift);
  const BtawlMixGeom g = btawl_mix_geometry(A, b, k, tau_s, tau_t, pi);
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0)) return 0.0;
  if (A <= BTAWL_A_EPS) {
    // With both channels present, H_mix can overshoot its sustained
    // asymptote.  First passage is governed by the minimum required launch
    // over the whole history, not by b/H_mix(t) on a later declining limb.
    const double vreq = btawl_mix_min_v(t, 0.0, g);
    return btawl_surv(vreq, p1, p2, launch, posdrift);
  }
  BtawlMixIntCtx ctx{&g, t, p1, p2, launch, posdrift};
  return std::fmin(1.0, btawl_mix_integrate(&btawl_mix_cdf_integrand, A, ctx) / A);
}

inline double btawl_mix_pdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift) {
  if (pi <= 1e-14)
    return btawl_pdf(t, A, b, p1, p2, k, tau_t, launch, posdrift);
  const BtawlMixGeom g = btawl_mix_geometry(A, b, k, tau_s, tau_t, pi);
  if (!g.ok || !(p2 > 0.0) || !(t > 0.0) || t == R_PosInf) return 0.0;
  if (A <= BTAWL_A_EPS) {
    const double h = btawl_hmix(t, k, tau_s, tau_t, pi);
    const double hp = btawl_hmix_p(t, k, tau_s, tau_t, pi);
    if (!(h > 0.0) || !(hp > 0.0)) return 0.0;
    // A mixture can in principle have a transient peak followed by a dip and
    // a later sustained rise.  Only a rising limb that is still the running
    // minimum of V*(u, 0) contributes first-passage density.
    const double t_prev = t * (1.0 - 1e-8);
    if (t_prev > 0.0 && btawl_mix_vstar(t, 0.0, g) >
        btawl_mix_min_v(t_prev, 0.0, g) * (1.0 + 1e-7)) return 0.0;
    return btawl_pdf_v(b / h, p1, p2, launch, posdrift) * b * hp / (h * h);
  }
  BtawlMixIntCtx ctx{&g, t, p1, p2, launch, posdrift};
  return btawl_mix_integrate(&btawl_mix_pdf_integrand, A, ctx) / A;
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

#endif

#ifndef EMC2_MODEL_BAWD_H
#define EMC2_MODEL_BAWD_H

// ---------------------------------------------------------------------------
// BAwD: the ballistic accumulator with drive decay
//       ("thresholded transient-drive ballistic accumulator")
//
//   U(t) = V e^{-k t}                      transient drive, launch strength V
//   X(t) = z + (V/k)(1 - e^{-k t}) - ell t start point z ~ U(0, A), static
//   threshold b = B + A
//
// The state does NOT leak (that is BAwL); the DRIVE decays while a constant
// clearance ell opposes it.  The trajectory therefore rises to a peak at
// log(V/ell)/k and falls, so the required launch strength for a crossing at u,
//
//   V*(u, z) = (b - z + ell u) / q(u),     q(u) = (1 - e^{-k u}) / k,
//
// is U-shaped rather than monotone.  Consequences, all of which the code below
// has to represent exactly:
//
//   * RTs have a HARD right endpoint T_max, the time at which V* is minimised
//     for z = 0.  Beyond it the CDF is flat and the density is exactly zero.
//   * Start points saturate one at a time: z = A first (at T_sat_A), z = 0 last
//     (at T_max).  In between, the CDF is a live closed-form part over
//     z in [0, Z(u)] plus a FROZEN part over z in (Z(u), A] whose hit
//     probability no longer changes.
//   * F(Inf) = F_max < 1.  The never-finish mass 1 - F_max is intrinsic.
//
// The flat-CDF contract is load-bearing for censoring and truncation: for every
// t >= T_max the CDF must return the same F_max bit-for-bit (evaluate the frozen
// branch, never extrapolate the live formula), and t = Inf must return F_max
// rather than 1.  An upper censoring bound beyond T_max is then well posed --
// the response mass above it is F_max - F(UC) = 0 and the omission mass is
// 1 - F_max.
//
// Two launch-strength distributions share all of the geometry above and differ
// only in the two integrals over V (ctx->bawd_launch selects between them):
//
//   BAWD_LAUNCH_NORMAL     V ~ N(v, sv^2), optionally truncated to V > 0.
//                          Nests BAwL exactly at ell = 0, A = 0, and LBA with
//                          drift v - ell at k = 0 (untruncated).  Its frozen
//                          integral has no elementary form -> Gauss-Legendre.
//   BAWD_LAUNCH_LOGNORMAL  log V ~ N(mu, sigma^2).  V > 0 by construction, so
//                          no truncation normalizer exists at all, and every
//                          regime is closed form in Phi/phi.  This is the
//                          default because its ordinary path needs no
//                          quadrature.
//
// Unlike BAwL this file has a SINGLE authoritative evaluation path (log space)
// rather than a guarded natural formula plus a log fallback.  The three
// regimes x two launch distributions already give six code paths; duplicating
// them in natural space is where the subtle bugs would live, and the log
// primitives it is built on (log_normal_phi_integral, log_normal_interval,
// log_normal_q_interval, log_lognormal_stoploss) are the same ones BAwL's log
// branch uses.  Natural-scale consumers exponentiate at the boundary.
//
// Like model_LBA.h (which it requires for the shared guard constants and the
// normalizer helpers) this header may be included by exactly one translation
// unit: particle_ll.cpp, directly and via utils.h.
// ---------------------------------------------------------------------------

#include <cmath>
#include "model_LBA.h"
#include "gl_quad.h"

constexpr double BAWD_K_EPS = 1e-10;
constexpr double BAWD_A_EPS = 1e-10;
constexpr double BAWD_ELL_EPS = 1e-12;
// Below this standardised start-point span the endpoint difference carries no
// information and the midpoint-in-z limit is used instead.
constexpr double BAWD_MIN_SPAN = 1e-8;
// Relative width below which a positive endpoint difference is replaced by its
// midpoint limit (which is then accurate to the square of the width).
constexpr double BAWD_MIN_LOG_GAP = 1e-7;
// Trust a signed-log numerator while it retains this (log-scale) fraction of
// its largest term; shared with BAwL.
constexpr double BAWD_LOG_BRACKET_MIN = BAWL_LOG_BRACKET_MIN;
// BAwD reuses BAwL's normalizer floor so that the ell = 0, A = 0 member
// reproduces pleakyba() exactly rather than to within a floor difference.
constexpr double BAWD_DENOM_FLOOR = BAWL_DENOM_FLOOR;

constexpr int BAWD_GL_NODES = 24;      // per panel, frozen integral
constexpr int BAWD_GL_MAX_NODES = 64;

// Launch-strength distribution selector (ContextForRaceModels::bawd_launch and
// the `launch` argument of the exported d/p functions -- keep them in sync).
constexpr int BAWD_LAUNCH_NORMAL = 0;
constexpr int BAWD_LAUNCH_LOGNORMAL = 1;

// --------------------------------------------------------------------------
// Geometry (launch-distribution free)
// --------------------------------------------------------------------------

// e^y - 1 - y, accurate for small y where expm1(y) - y cancels.
inline double bawd_em1my(double y) {
  if (std::fabs(y) < 1e-4) {
    // (y^2/2) (1 + y/3 + y^2/12 + y^3/60)
    return 0.5 * y * y * (1.0 + (y / 3.0) * (1.0 + (y / 4.0) * (1.0 + y / 5.0)));
  }
  return std::expm1(y) - y;
}

// Solve e^y - 1 - y = c for y >= 0.  This is -W_{-1}(-e^{-(1+c)}) computed by
// Newton on a convex increasing function, which needs no branch selection: the
// Lambert-W closed form would require choosing W_{-1} over W_0 at every call
// site and is less robust near c = 0.
inline double bawd_newton_y(double c) {
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

// Row geometry: everything that depends on (A, b, k, ell) but not on u or on
// the launch distribution.  Two Newton solves per row, never per quadrature
// node, so the truncation and censoring paths that evaluate the CDF at several
// times for one row reuse them.
struct BawdGeom {
  bool ok = false;
  bool k_zero = false;      // exact LBA-with-drift-(V - ell) limit
  bool ell_zero = false;    // no clearance: no saturation, T_max = Inf
  double b = 0.0, A = 0.0, k = 0.0, ell = 0.0;
  double y_A = 0.0;         // saturation exponent for z = A (earliest)
  double y_0 = 0.0;         // saturation exponent for z = 0 (latest)
  double T_sat_A = R_PosInf;
  double T_max = R_PosInf;
  double V_c0 = R_PosInf;   // critical launch strength at z = 0
};

inline BawdGeom bawd_geometry(double A, double b, double k, double ell) {
  BawdGeom g;
  if (!emc2_isfinite(A) || !emc2_isfinite(b) || !emc2_isfinite(k) ||
      !emc2_isfinite(ell))
    return g;
  if (!(b > 0.0) || !(A >= 0.0) || !(b >= A) || !(k >= 0.0) || !(ell >= 0.0))
    return g;
  g.ok = true;
  g.b = b; g.A = A; g.k = k; g.ell = ell;
  g.k_zero = (k <= BAWD_K_EPS);
  g.ell_zero = (ell <= BAWD_ELL_EPS);
  if (g.k_zero || g.ell_zero) return g;  // nothing ever saturates

  g.y_0 = bawd_newton_y(k * b / ell);
  g.y_A = (b > A) ? bawd_newton_y(k * (b - A) / ell) : 0.0;
  g.T_max = g.y_0 / k;
  g.T_sat_A = g.y_A / k;
  // V_c(z) = ell e^y with e^y - 1 - y = k(b - z)/ell, so
  // ell e^y = k(b - z) + ell(1 + y) exactly -- no exp() to overflow, and the
  // form makes the BAwL-like V > k*b condition visible at ell -> 0.
  g.V_c0 = k * b + ell * (1.0 + g.y_0);
  return g;
}

inline double bawd_critical_launch(const BawdGeom& g, double y) {
  // ell e^y with e^y = (e^y - 1 - y) + 1 + y.
  return g.ell * (bawd_em1my(y) + 1.0 + y);
}

// Per-time quantities.  `Z` splits the start-point range into the live part
// [0, Z] and the already-saturated (frozen) part (Z, A].
struct BawdAtU {
  bool ok = false;
  bool saturated = false;   // u >= T_max: Z = 0, only frozen mass remains
  bool partial = false;     // Z < A: a frozen contribution exists
  double q = 0.0;           // (1 - e^{-k u}) / k, = u at k = 0
  double E = 1.0;           // e^{-k u}
  // log E kept alongside E: e^{-k u} underflows to zero for k u > 745 while the
  // density it multiplies is still representable in log space, so every term
  // that carries a bare factor of E is formed as log(coef) + log_E instead.
  double log_E = 0.0;
  double Z = 0.0;
  double w_hi = 0.0;        // V*(u, 0), required launch at the lowest start
  double w_lo = 0.0;        // V*(u, Z)
  double s_lo = 0.0;        // frozen-integral exponent limits (y coordinates)
  double s_hi = 0.0;
};

inline BawdAtU bawd_at_u(const BawdGeom& g, double u) {
  BawdAtU s;
  if (!g.ok || !(u > 0.0) || ISNAN(u)) return s;
  s.ok = true;
  const bool inf_u = (u == R_PosInf);

  if (g.k_zero) {
    s.q = u;                 // Inf at u = Inf; only ever used via w below
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

  if (g.k_zero || g.ell_zero) {
    s.Z = g.A;
    s.partial = false;
    s.saturated = false;
  } else {
    const double ku = g.k * u;
    if (inf_u || ku >= g.y_0) {
      s.Z = 0.0;
      s.saturated = true;
    } else {
      // Saturation is z*(u) <= 0, NOT Z <= 0: with a point start (A = 0) the
      // clamped Z is identically zero while the accumulator is still live.
      const double zstar = g.b - (g.ell / g.k) * bawd_em1my(ku);
      s.saturated = !(zstar > 0.0);
      s.Z = std::fmin(std::fmax(zstar, 0.0), g.A);
    }
    s.partial = (s.Z < g.A);
    s.s_lo = g.y_A;
    s.s_hi = inf_u ? g.y_0 : std::fmin(ku, g.y_0);
  }

  if (s.saturated) {
    s.w_hi = g.V_c0;
    s.w_lo = g.V_c0;
  } else if (g.k_zero) {
    // (b - z)/u + ell, exact at u = Inf (where it is just ell).
    s.w_hi = g.b / u + g.ell;
    s.w_lo = (g.b - s.Z) / u + g.ell;
  } else {
    const double ellu = g.ell_zero ? 0.0 : g.ell * u;
    s.w_hi = (g.b + ellu) / s.q;
    s.w_lo = (g.b - s.Z + ellu) / s.q;
  }
  return s;
}

// --------------------------------------------------------------------------
// Positive-integrand log quadrature.  One routine serves the truncated-normal
// frozen integral (its ordinary path) and the lognormal frozen fallback, so
// there is a single place where node counts and panel splitting live.
// --------------------------------------------------------------------------

template <typename LogFun>
inline double bawd_log_gl(const LogFun& log_f, double a, double b, int n) {
  if (!(b > a) || !emc2_isfinite(a) || !emc2_isfinite(b)) return R_NegInf;
  if (n > BAWD_GL_MAX_NODES) n = BAWD_GL_MAX_NODES;
  const GLRule& r = gl_get_rule(n);
  const double c1 = 0.5 * (b - a), c2 = 0.5 * (b + a);
  double terms[BAWD_GL_MAX_NODES];
  double best = R_NegInf;
  for (int i = 0; i < n; ++i) {
    const double lf = log_f(c1 * r.x[i] + c2);
    terms[i] = (r.w[i] > 0.0 && lf > R_NegInf && !ISNAN(lf))
      ? (std::log(r.w[i]) + lf) : R_NegInf;
    if (terms[i] > best) best = terms[i];
  }
  if (!(best > R_NegInf)) return R_NegInf;
  double acc = 0.0;
  for (int i = 0; i < n; ++i)
    if (terms[i] > R_NegInf) acc += std::exp(terms[i] - best);
  return std::log(c1) + best + std::log(acc);
}

// Two panels split at `mid` when it lies strictly inside; the split point is
// where the integrand's normal factor turns over, which is where a single
// fixed rule is least accurate.
template <typename LogFun>
inline double bawd_log_gl_split(const LogFun& log_f, double a, double b,
                                double mid, int n) {
  if (mid > a && mid < b) {
    return log_sum_exp(bawd_log_gl(log_f, a, mid, n),
                       bawd_log_gl(log_f, mid, b, n));
  }
  return bawd_log_gl(log_f, a, b, n);
}

// --------------------------------------------------------------------------
// Frozen contribution: the mass of start points that have already saturated.
//
//   Psi = (ell/k) int_{s_lo}^{s_hi} (e^s - 1) P(V >= ell e^s) ds
//
// Written in the exponent s = k * T_sat(z) rather than in z, so no root
// finding appears inside the integral: the substitution
// z = b - (ell/k)(e^s - 1 - s) is explicit.
// --------------------------------------------------------------------------

inline double bawd_log_frozen_normal(const BawdGeom& g, double s_lo, double s_hi,
                                     double v, double sv) {
  if (!(s_hi > s_lo) || !(sv > 0.0)) return R_NegInf;
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    return std::log(e1) +
      pnorm_log_direct((v - bawd_critical_launch(g, s)) / sv, true);
  };
  // ell e^s = v is where the normal factor passes its midpoint.
  const double mid = (v > 0.0) ? std::log(v / g.ell) : s_lo;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log(g.ell) - std::log(g.k);
}

inline double bawd_log_frozen_logn_quad(const BawdGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma) {
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  const double x_ell = (std::log(g.ell) - mu) / sigma;
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    // P(V >= ell e^s) = Q(x_ell + s/sigma).
    return std::log(e1) + pnorm_log_direct(x_ell + s / sigma, false);
  };
  const double mid = -sigma * x_ell;   // ell e^s = median(V)
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log(g.ell) - std::log(g.k);
}

// Closed form for the lognormal launch.  Substituting w = ell e^s turns Psi
// into (1/k) [I_0 - ell I_{-1}] over the critical-launch interval, where
//   I_0    = int Gbar(w) dw       = C(w_a) - C(w_b)   (the stop-loss price)
//   I_{-1} = int Gbar(w)/w dw     = sigma [D(x_a) - D(x_b)]
// and D(x) = phi(x) - x Q(x) is log_normal_q_antiderivative_abs().  Both
// pieces therefore land on primitives that already exist.
//
// I_0 - ell I_{-1} = int Gbar(w)(1 - ell/w) dw >= 0 because every critical
// launch satisfies w >= ell, but it is ill conditioned exactly when the
// interval concentrates at w = ell, and no rearrangement of two separately
// evaluated terms recovers relative precision there.  That case falls back to
// the combined positive integrand above, which cancels nowhere.
inline double bawd_log_frozen_logn(const BawdGeom& g, double s_lo, double s_hi,
                                   double mu, double sigma) {
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  const double w_a = bawd_critical_launch(g, s_lo);
  const double w_b = bawd_critical_launch(g, s_hi);
  if (!(w_b > w_a) || !emc2_isfinite(w_b))
    return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);

  const double lc_a = log_lognormal_stoploss(w_a, mu, sigma);
  const double lc_b = log_lognormal_stoploss(w_b, mu, sigma);
  if (lc_a - lc_b > BAWD_MIN_LOG_GAP) {
    const double L0 = log_diff_exp(lc_a, lc_b);
    const double x_a = (std::log(w_a) - mu) / sigma;
    const double x_b = (std::log(w_b) - mu) / sigma;
    const double L1 = std::log(g.ell) + std::log(sigma) +
      log_normal_q_interval(x_a, x_b);
    if (!ISNAN(L0) && emc2_isfinite(L0)) {
      const double d = L1 - L0;
      // Retained fraction of I_0 - ell I_{-1}; below 1e-6 the closed pair has
      // lost its digits (plan section 2c, layer 3).
      if (d < 0.0 && -std::expm1(d) > 1e-6)
        return L0 + log1m_exp(d) - std::log(g.k);
      if (!(L1 > R_NegInf))  // ell I_{-1} underflowed: I_0 stands alone
        return L0 - std::log(g.k);
    }
  }
  return bawd_log_frozen_logn_quad(g, s_lo, s_hi, mu, sigma);
}

// --------------------------------------------------------------------------
// log CDF
// --------------------------------------------------------------------------

inline double log_bawd_cdf_normal(double u, const BawdGeom& g, double v,
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

inline double log_bawd_cdf_logn(double u, const BawdGeom& g, double mu,
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

// --------------------------------------------------------------------------
// log PDF
//
// The density numerator is signed: unlike BAwL's (b m + c) the lead coefficient
// (v E - ell) genuinely goes negative near T_max, so the cancellation guard and
// its midpoint-in-z fallback are required, not optional.  Both fallbacks use
//   f A = int_0^Z phi/g(V*) * (-dV*/du) dz,   -dV*/du = (V* E - ell)/q,
// whose weight is non-negative on [0, Z] and vanishes exactly at Z when the
// row is partially saturated.
// --------------------------------------------------------------------------

inline double log_bawd_pdf_normal(double u, const BawdGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor) {
  if (!g.ok || !(sv > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;   // exactly zero past T_max
  const double log_denom = log_positive_normalizer(v, sv, posdrift, denom_floor);

  if (g.A <= BAWD_A_EPS) {
    // The weight is V*(u,0) E - ell.  With no clearance E factors out cleanly,
    // which keeps the deep tail representable once E itself has underflowed.
    const double log_wgt = g.ell_zero
      ? std::log(s.w_hi) + s.log_E
      : ((s.w_hi * s.E - g.ell) > 0.0 ? std::log(s.w_hi * s.E - g.ell)
                                      : R_NegInf);
    if (!(log_wgt > R_NegInf)) return R_NegInf;
    return dnormP((v - s.w_hi) / sv, 0.0, 1.0, true) + log_wgt -
      std::log(sv) - std::log(s.q) - log_denom;
  }

  const double cl = (v - s.w_hi) / sv;
  const double ch = (v - s.w_lo) / sv;
  const double span = ch - cl;
  const double log_scale = -std::log(g.A) - log_denom;

  if (span > BAWD_MIN_SPAN && emc2_isfinite(span)) {
    const double log_dphi = log_normal_interval(cl, ch);
    if (!ISNAN(log_dphi)) {
      // Numerator (v E - ell) DeltaPhi + sv E Delta(phi), regrouped as
      //   E (v DeltaPhi + sv Delta(phi)) - ell DeltaPhi
      // so the single factor of E is applied in log space and does not
      // underflow the whole term once k u exceeds ~745.
      const signed_log dphi = signed_log_sub(
        make_signed_log(dnormP(ch, 0.0, 1.0, true), 1),
        make_signed_log(dnormP(cl, 0.0, 1.0, true), 1));
      const signed_log v_term = signed_log_product(v, log_dphi);
      const signed_log sv_term = (dphi.sign != 0)
        ? make_signed_log(std::log(sv) + dphi.log_abs, dphi.sign)
        : make_signed_log(R_NegInf, 0);
      signed_log lead = signed_log_add(v_term, sv_term);
      if (lead.sign != 0) lead.log_abs += s.log_E;
      const signed_log clear_term = signed_log_product(-g.ell, log_dphi);
      const signed_log bracket = signed_log_add(lead, clear_term);
      const double max_term = std::fmax(lead.log_abs, clear_term.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWD_LOG_BRACKET_MIN) {
        return bracket.log_abs + log_scale;
      }
      // Cancelled past the accuracy of the tail logs (the true numerator is
      // positive).  Midpoint in z: the weight (V* E - ell) is monotone and
      // bounded on [0, Z], and DeltaPhi is exact.
      const double w_mid = 0.5 * (s.w_hi + s.w_lo);
      const double wgt = w_mid * s.E - g.ell;
      if (!(wgt > 0.0)) return R_NegInf;
      return std::log(wgt) + log_dphi + log_scale;
    }
  }
  // Collapsed start range: point limit in z, weighted by the span.
  const double w_mid = 0.5 * (s.w_hi + s.w_lo);
  const double wgt = w_mid * s.E - g.ell;
  if (!(wgt > 0.0)) return R_NegInf;
  return std::log(s.Z) + std::log(wgt) - std::log(sv) - std::log(s.q) +
    dnormP(0.5 * (cl + ch), 0.0, 1.0, true) + log_scale;
}

inline double log_bawd_pdf_logn(double u, const BawdGeom& g, double mu,
                                double sigma) {
  if (!g.ok || !(sigma > 0.0) || !(u > 0.0) || u == R_PosInf) return R_NegInf;
  const BawdAtU s = bawd_at_u(g, u);
  if (!s.ok || s.saturated) return R_NegInf;

  if (g.A <= BAWD_A_EPS) {
    const double log_wgt = g.ell_zero
      ? std::log(s.w_hi) + s.log_E
      : ((s.w_hi * s.E - g.ell) > 0.0 ? std::log(s.w_hi * s.E - g.ell)
                                      : R_NegInf);
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

  if (!ISNAN(log_prob) && s.log_E > R_NegInf) {
    const double log_pexp = log_normal_interval(d1_hi, d1_lo);
    if (!ISNAN(log_pexp)) {
      const signed_log t1 = make_signed_log(
        s.log_E + mu + 0.5 * sigma * sigma + log_pexp, 1);
      const signed_log t2 = signed_log_product(g.ell, log_prob);
      const signed_log bracket = signed_log_sub(t1, t2);
      const double max_term = std::fmax(t1.log_abs, t2.log_abs);
      if (bracket.sign > 0 && !ISNAN(bracket.log_abs) &&
          bracket.log_abs > max_term + BAWD_LOG_BRACKET_MIN) {
        return bracket.log_abs + log_scale;
      }
    }
  }
  // Midpoint in w: (w_mid E - ell) * P(w_lo <= V <= w_hi).
  const double w_mid = 0.5 * (s.w_hi + s.w_lo);
  const double wgt = w_mid * s.E - g.ell;
  if (!(wgt > 0.0)) return R_NegInf;
  if (!ISNAN(log_prob) && log_prob > R_NegInf)
    return std::log(wgt) + log_prob + log_scale;
  return std::log(s.Z) + std::log(wgt) + dlnorm_std(w_mid, mu, sigma, true) -
    std::log(s.q) + log_scale;
}

// --------------------------------------------------------------------------
// Dispatch and output wrappers.  `p1`/`p2` are (v, sv) for the normal launch
// and (mu, sigma) for the lognormal one; they occupy the same kernel columns.
// --------------------------------------------------------------------------

inline double bawd_log_cdf(double u, double A, double b, double p1, double p2,
                           double k, double ell, int launch, bool posdrift,
                           double denom_floor = BAWD_DENOM_FLOOR) {
  const BawdGeom g = bawd_geometry(A, b, k, ell);
  if (launch == BAWD_LAUNCH_LOGNORMAL)
    return log_bawd_cdf_logn(u, g, p1, p2);
  return log_bawd_cdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

inline double bawd_log_pdf(double u, double A, double b, double p1, double p2,
                           double k, double ell, int launch, bool posdrift,
                           double denom_floor = BAWD_DENOM_FLOOR) {
  const BawdGeom g = bawd_geometry(A, b, k, ell);
  if (launch == BAWD_LAUNCH_LOGNORMAL)
    return log_bawd_pdf_logn(u, g, p1, p2);
  return log_bawd_pdf_normal(u, g, p1, p2, posdrift, denom_floor);
}

inline double bawd_cdf_norm(double t, double A, double b, double p1, double p2,
                            double k, double ell, int launch, bool posdrift,
                            bool log_out,
                            double denom_floor = BAWD_DENOM_FLOOR) {
  return return_from_log(
    bawd_log_cdf(t, A, b, p1, p2, k, ell, launch, posdrift, denom_floor),
    log_out);
}

inline double bawd_pdf_norm(double t, double A, double b, double p1, double p2,
                            double k, double ell, int launch, bool posdrift,
                            bool log_out,
                            double denom_floor = BAWD_DENOM_FLOOR) {
  return return_from_log(
    bawd_log_pdf(t, A, b, p1, p2, k, ell, launch, posdrift, denom_floor),
    log_out);
}

// Natural-scale scalar evaluators for consumers that clamp to [0, 1] and
// tolerate tail saturation: truncation normalisers and GSL integrands.
inline double bawd_cdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, double ell,
                                      int launch, bool posdrift) {
  const double lp = bawd_log_cdf(t, A, b, p1, p2, k, ell, launch, posdrift);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

inline double bawd_pdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, double ell,
                                      int launch, bool posdrift) {
  const double lp = bawd_log_pdf(t, A, b, p1, p2, k, ell, launch, posdrift);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

// --------------------------------------------------------------------------
// R-callable entry points.  These bypass ContextForRaceModels entirely, so the
// launch distribution MUST be passed explicitly; R/model_BAwD.R derives both
// this argument and the c_name suffix from one `drift_distribution` value so
// dfun/pfun cannot silently disagree with the sampled likelihood.
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dbawd(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    NumericVector ell, int launch = 1, bool posdrift = true,
                    bool log_out = false) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), pick(ell, i), launch,
                           posdrift, log_out);
  return out;
}

// [[Rcpp::export]]
NumericVector pbawd(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    NumericVector ell, int launch = 1, bool posdrift = true,
                    bool log_out = false) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), pick(ell, i), launch,
                           posdrift, log_out);
  return out;
}

// [[Rcpp::export]]
double dbawd_norm(double t, double A, double b, double p1, double p2, double k,
                  double ell, int launch = 1, bool posdrift = true,
                  bool log_out = false) {
  return bawd_pdf_norm(t, A, b, p1, p2, k, ell, launch, posdrift, log_out);
}

// [[Rcpp::export]]
double pbawd_norm(double t, double A, double b, double p1, double p2, double k,
                  double ell, int launch = 1, bool posdrift = true,
                  bool log_out = false) {
  return bawd_cdf_norm(t, A, b, p1, p2, k, ell, launch, posdrift, log_out);
}

// Right endpoint of the supported decision-time window (Inf when ell = 0 or
// k = 0, where nothing saturates).  Exported because the flat-CDF contract and
// the support constraint are both stated in terms of it, and a test that
// recomputed it in R would not be testing the same quantity.
// [[Rcpp::export]]
double bawd_tmax(double A, double b, double k, double ell) {
  const BawdGeom g = bawd_geometry(A, b, k, ell);
  return g.ok ? g.T_max : NA_REAL;
}

// log E[(V - v)_+] for log V ~ N(mu, sigma^2); exposed so the cancellation
// layer can be checked against a high-precision reference from R.
// [[Rcpp::export]]
double lognormal_stoploss_log(double v, double mu, double sigma) {
  return log_lognormal_stoploss(v, mu, sigma);
}

#endif  // EMC2_MODEL_BAWD_H

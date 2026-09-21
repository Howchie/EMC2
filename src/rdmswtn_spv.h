#ifndef RDMSWTN_SPV_H
#define RDMSWTN_SPV_H

#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include "utility_functions.h"
#include "composite_functions.h"
#include "wald_functions.h"
#include "gaussian.h"
#include "gl_quad.h"

// ==========================================================================
// Wald first passage averaged over a uniform distance D ~ U[b - A, b] and a
// Gaussian drift V ~ N(mu, sv^2), the drift law either truncated to V > 0
// (renormalised by Z = Phi(mu/sv)) or full-Gaussian (defective, never
// renormalised).  The no-clock decision CDF and survivor at decision time t
// come from the distance recurrence of Math/RDMSWTN_Derivation.tex.  At a
// fixed distance r, with c = s sqrt(t),
//
//   T(r) = E[1{V>0} Phi((Vt - r)/c)],       U(r) = E[1{V>0} Phi((r - Vt)/c)],
//   H(r) = E[1{V>0} exp(2rV/s^2) Phi(-(r + Vt)/c)],
//
// so F(r) = (T + H)/Z and S(r) = (U - H)/Z (full-Gaussian: drop the
// indicator, Z = 1).  T, U and H close under d/dr together with four
// Gaussian functions of r, every coefficient affine in r, so on a step
// r = r0 + h x their Taylor coefficients obey a two-term recurrence and the
// step average is the coefficient sum weighted by 1/(n + 1).
//
// Numerics.
//  * Every step is carried on its own log scale, so a CDF far below the
//    smallest double at tiny t, or a survivor far below it at long t, keeps
//    its digits.
//  * F integrates T + H downwards from the top of the range.  S integrates U
//    upwards from the bottom and H downwards from the top, then takes the
//    difference of the two integrals: every carried quantity grows in its
//    direction of travel and H's homogeneous mode decays downwards, and S
//    keeps its relative accuracy when F is close to one.
//  * Step sizes adapt to the local rates and curvatures of the Gaussian
//    factors and to the stiffness of H; each step's series runs until its
//    remainder, including what the auxiliary factors can still feed into
//    T, U or H, is below 1e-17 of the step average.
//  * F is decreasing and S increasing in r, so when the integrand is steep
//    the end of the range carrying less than 1e-16 of the mass is skipped
//    and that bound is verified afterwards from the endpoint value.
//  * Initial values: ordinary normal CDFs for full-Gaussian drift.  Under
//    positive truncation, bivariate normal probabilities where the tvpack
//    series keeps relative accuracy (m >= 0, correlation below 0.925, the
//    univariate tail argument above -4, no cancellation), otherwise the same
//    quantity written as a posterior-normal expectation of a Mills ratio: in
//    closed form when the drift truncation is negligible, as that closed
//    form less a one-sided integral when the truncated mass is small, and
//    otherwise integrated around the mode of its log-concave integrand.
// ==========================================================================

namespace rdmswtn_spv {

#ifdef RDMSWTN_SPV_STATS
struct Stats { long seeds_bvn = 0, seeds_post = 0, steps = 0, terms = 0, mills = 0; };
inline Stats& stats() { static Stats st; return st; }
#define RDMSWTN_SPV_COUNT(field, n) (rdmswtn_spv::stats().field += (n))
#else
#define RDMSWTN_SPV_COUNT(field, n) ((void)0)
#endif

// Per-step exponent budgets.
constexpr double kDecayBudget = 7.0;    // |rate h| for a factor shrinking along the step
constexpr double kGrowBudget  = 30.0;   // |rate h| for a factor growing along the step
constexpr double kStiffBudget = 30.0;   // (a + beta r)|h| for H's decaying homogeneous mode
constexpr double kCurvBudget  = 6.0;    // curvature * h^2 of a Gaussian factor
constexpr double kNegligible  = -45.0;  // log weight below which a factor is dropped
constexpr int    kMaxTerms    = 160;
constexpr double kTermTol     = 1e-17;
constexpr int    kMaxSteps    = 200000;
// A tail of the distance range is skipped when a Gaussian proxy puts it
// kSkipDrop below the dominant end; the skip is kept only if monotonicity
// bounds the discarded integral by kSkipTol of the retained one.
constexpr double kSkipDrop    = 55.0;
constexpr double kSkipTol     = 1e-16;
// Initial-value routing.
constexpr double kBvnRhoMax   = 0.925;  // tvpack's low-correlation branch
constexpr double kBvnMinArg   = -4.0;   // tvpack keeps ~3e-13 relative accuracy above this tail argument
constexpr double kBvnMaxRatio = 0.95;
constexpr int    kPostNodes   = 24;     // Gauss-Legendre nodes per posterior panel
constexpr double kPostPanel   = 12.0;   // width of the first posterior panel
constexpr double kPostDrop    = 40.0;   // log-integrand drop that ends a posterior window

// Normal tail functions at full double precision, through the Mills ratio
// R(x) = Phi(-x)/phi(x).  On [-3, 26) log R comes from a Chebyshev table
// (degree 12 on half-unit intervals, about 4e-16 absolute) built once from
// extended-precision erfc; beyond 26 the Laplace continued fraction is used,
// and below -3 the complement of the rational fast_norm_phi, whose relative
// error there is irrelevant next to 1.
inline double mills_cf_inverse(double x) {  // 1/R(x) for x >= 26
  double v = x;
  for (int k = 14; k >= 1; --k) v = x + k / v;
  return v;
}

struct MillsTable {
  static constexpr double lo = -3.0, hi = 26.0;
  static constexpr int per = 2, deg = 12, n = static_cast<int>((hi - lo) * per);
  double c[n][deg];
  MillsTable() {
    const long double hw = 0.5L / per, pi = 3.141592653589793238462643383279502884L;
    for (int k = 0; k < n; ++k) {
      const long double ctr = lo + (k + 0.5L) / per;
      long double f[deg];
      for (int j = 0; j < deg; ++j) {
        const long double x = ctr + hw * std::cos(pi * (j + 0.5L) / deg);
        f[j] = 0.5L * std::log(pi / 2) + 0.5L * x * x + std::log(std::erfc(x / std::sqrt(2.0L)));
      }
      for (int m = 0; m < deg; ++m) {
        long double acc = 0;
        for (int j = 0; j < deg; ++j) acc += f[j] * std::cos(pi * m * (j + 0.5L) / deg);
        c[k][m] = static_cast<double>(acc * (m == 0 ? 1.0L : 2.0L) / deg);
      }
    }
  }
  double operator()(double x) const {   // lo <= x < hi
    const double y = (x - lo) * per;
    const int k = std::min(static_cast<int>(y), n - 1);
    const double u = 2.0 * (y - k) - 1.0, u2 = 2.0 * u;
    const double* a = c[k];
    double b1 = 0.0, b2 = 0.0;
    for (int m = deg - 1; m >= 1; --m) {
      const double b0 = std::fma(u2, b1, a[m] - b2);
      b2 = b1;
      b1 = b0;
    }
    return std::fma(u, b1, a[0] - b2);
  }
};
inline const MillsTable& mills_table() { static const MillsTable t; return t; }

// log R(x) for either sign of x.
inline double log_mills(double x) {
  RDMSWTN_SPV_COUNT(mills, 1);
  if (x >= MillsTable::hi) return -std::log(mills_cf_inverse(x));
  if (x >= MillsTable::lo) return mills_table()(x);
  return std::log1p(-fast_norm_phi(x)) + 0.5 * x * x + LOG_SQRT_2PI;
}

inline double log_Phi(double x) {
  if (x <= 0.0) return log_phi_std(x) + log_mills(-x);
  return std::log1p(-std::exp(log_phi_std(x) + log_mills(x)));
}

// Normal hazard phi(x)/Phi(-x) = 1/R(x).
inline double normal_hazard(double x) { return std::exp(-log_mills(x)); }

struct Law {
  double t, mu, s, sv;
  bool pos;
  double s2, sv2, Q, sqrtQ, S, S2, c, c2, mut, rho, u, k, a, beta;
  double log_delta;  // log(2 sv phi(mu/sv) / s^2); -Inf for full-Gaussian drift
  double logZ;       // log Phi(mu/sv); 0 for full-Gaussian drift
  double log_c, log_S;
  // tvpack's low-correlation series for Phi2(., .; rho) at this law's rho:
  // Phi2(x, y) = Phi(x) Phi(y) + sum_j bw_j exp((bs_j x y - (x^2 + y^2)/2) bi_j).
  // Built only when positive-truncation initial values can use it.
  int bn;
  double bs[20], bi[20], bw[20];
};

inline Law make_law(double t, double mu, double s, double sv, bool pos) {
  Law L;
  L.t = t; L.mu = mu; L.s = s; L.sv = sv; L.pos = pos;
  L.s2 = s * s;
  L.sv2 = sv * sv;
  L.Q = L.s2 + L.sv2 * t;
  L.sqrtQ = std::sqrt(L.Q);
  L.S2 = t * L.Q;
  L.S = std::sqrt(L.S2);
  L.c2 = L.s2 * t;
  L.c = std::sqrt(L.c2);
  L.mut = mu * t;
  L.rho = sv * std::sqrt(t) / L.sqrtQ;
  L.u = sv / (s * L.sqrtQ);
  L.k = 1.0 + 2.0 * L.sv2 * t / L.s2;
  L.a = 2.0 * mu / L.s2;
  L.beta = 4.0 * L.sv2 / (L.s2 * L.s2);
  L.log_delta = pos ? std::log(2.0 * sv / L.s2) + log_phi_std(mu / sv) : R_NegInf;
  L.logZ = pos ? log_Phi(mu / sv) : 0.0;
  L.log_c = std::log(L.c);
  L.log_S = std::log(L.S);
  L.bn = 0;
  if (pos && mu >= 0.0 && L.rho < kBvnRhoMax &&
      L.rho > std::numeric_limits<double>::epsilon()) {
    const int order = (L.rho < 0.3) ? 6 : ((L.rho < 0.5) ? 12 : 20);
    const LegendreHalfRule& rule = get_legendre_half_rule(order);
    const double asr = std::asin(L.rho);
    const double half = 0.5 * asr;
    const double scale = asr / (4.0 * M_PI);
    // Every angle below is at most asin(kBvnRhoMax)/2 < 0.6, where these
    // Taylor series are exact to double precision.
    auto sin_cos = [](double a, double& sa, double& ca) {
      const double z = a * a;
      sa = a * (1.0 + z * (-1.0 / 6 + z * (1.0 / 120 + z * (-1.0 / 5040 + z * (1.0 / 362880 +
               z * (-1.0 / 39916800 + z * (1.0 / 6227020800.0 + z * (-1.0 / 1307674368000.0))))))));
      ca = 1.0 + z * (-0.5 + z * (1.0 / 24 + z * (-1.0 / 720 + z * (1.0 / 40320 + z * (-1.0 / 3628800 +
               z * (1.0 / 479001600 + z * (-1.0 / 87178291200.0 + z * (1.0 / 20922789888000.0))))))));
    };
    double sh, ch;
    sin_cos(half, sh, ch);
    // sin(half (1 -+ x)) = sin(half) cos(half x) -+ cos(half) sin(half x).
    const int m = static_cast<int>(rule.x.size());
    const double* xs = rule.x.data();
    const double* ws = rule.w.data();
    for (int ix = 0; ix < m; ++ix) {
      double sa, ca;
      sin_cos(half * xs[ix], sa, ca);
      L.bs[2 * ix] = sh * ca - ch * sa;
      L.bs[2 * ix + 1] = sh * ca + ch * sa;
      L.bw[2 * ix] = L.bw[2 * ix + 1] = scale * ws[ix];
    }
    L.bn = 2 * m;
    for (int j = 0; j < L.bn; ++j) L.bi[j] = 1.0 / std::fma(-L.bs[j], L.bs[j], 1.0);
  }
  return L;
}

// Phi2(x, y; rho) at the law's correlation (rho < kBvnRhoMax).
inline double bvn(const Law& L, double x, double y) {
  const double hk = x * y, hs = 0.5 * (x * x + y * y);
  double acc = 0.0;
  for (int j = 0; j < L.bn; ++j)
    acc += L.bw[j] * std::exp(std::fma(L.bs[j], hk, -hs) * L.bi[j]);
  return std::fma(pnorm_std(x), pnorm_std(y), acc);
}

// log of  int_{zlo}^inf phi(z) R(x0 + sl z) dz  for |sl| < 1.
// (log R)'' lies in (0, 1), so the log integrand f has f'' in
// (-1, -(1 - sl^2)): it is strictly concave.
//
// Without truncation the integral is closed form: R(x) = int_0^inf
// exp(-x u - u^2/2) du, so E[R(x0 + sl Z)] = R(x0/k)/k with k = sqrt(1 - sl^2).
// When f is still rising at zlo, concavity bounds the mass below zlo by
// exp(f(zlo))/f'(zlo).  If that is below exp(-kPostDrop) of the closed form,
// the closed form is returned; if it is below half, the integral is the
// closed form minus the mass below zlo, which is the same kind of integral
// mirrored (z -> -z) and has its maximum at its boundary.
//
// Otherwise the mode is found by safeguarded Newton and the integration
// window runs to where f has fallen by kPostDrop, or to zlo.  Each side of
// the mode is integrated by Gauss-Legendre with nodes in proportion to its
// length (every feature is at least of order one wide because |f''| <= 1); a
// side longer than kPostPanel (possible only as sl^2 -> 1, when it can
// stretch over hundreds of units with its curvature concentrated next to the
// mode) is cut into panels that double in width away from the mode.
inline double log_posterior_mills(double x0, double sl, double zlo,
                                  int nodes = kPostNodes,
                                  double first_panel = kPostPanel,
                                  bool allow_closed_form = true) {
  constexpr double kLaplaceReach = 9.0;
  auto f = [&](double z) { return log_phi_std(z) + log_mills(x0 + sl * z); };
  auto df = [&](double z, double& d2) {
    const double x = x0 + sl * z;
    const double lam = normal_hazard(x);
    d2 = std::fmin(-1.0 + sl * sl * (1.0 - lam * (lam - x)), -1e-12);
    return -z + sl * (x - lam);
  };

  const double om = std::fma(-sl, sl, 1.0);
  if (allow_closed_form && om > 0.0) {
    const double k = std::sqrt(om);
    const double full = log_mills(x0 / k) - std::log(k);
    if (!(zlo > R_NegInf)) return full;
    const double x = x0 + sl * zlo;
    const double lm = log_mills(x);
    const double slope = -zlo + sl * (x - std::exp(-lm));
    if (slope > 0.0) {
      const double log_bound = log_phi_std(zlo) + lm - std::log(slope);
      if (log_bound <= full - kPostDrop) return full;
      if (log_bound < full - M_LN2) {
        // The mass below zlo, int_{-inf}^{zlo} = I(x0, -sl, -zlo), starts at
        // its own maximum: one side of quadrature, and less than half of full.
        const double log_below = log_posterior_mills(x0, -sl, -zlo, nodes, first_panel, false);
        return full + std::log1p(-std::exp(log_below - full));
      }
    }
  }

  // Mode.  f' is decreasing, so a bracket [lo, hi] with f'(lo) > 0 > f'(hi)
  // is kept and any Newton step leaving it is replaced by bisection or, while
  // the bracket is open, by an expanding step.  The mode only places the
  // quadrature, so it is not refined past 1e-7.
  double d2 = -1.0;
  double lo = zlo, hi = R_PosInf;
  double z = std::max(zlo, sl * (x0 - normal_hazard(x0)));
  bool at_boundary = false;
  double g = 0.0;
  for (int it = 0; it < 200; ++it) {
    g = df(z, d2);
    if (g <= 0.0 && z <= zlo) { at_boundary = true; break; }
    if (g > 0.0) lo = z; else hi = z;
    double zn = z - g / d2;
    if (!(zn > lo && zn < hi)) {
      if (R_FINITE(hi) && R_FINITE(lo)) zn = 0.5 * (lo + hi);
      else if (!R_FINITE(hi)) zn = z + std::max(1.0, 2.0 * std::fabs(zn - z));
      else zn = std::max(zlo, z - std::max(1.0, 2.0 * std::fabs(zn - z)));
    }
    if (std::fabs(zn - z) <= 1e-7 * (1.0 + std::fabs(z))) break;
    z = zn;
  }
  const double z_mode = at_boundary ? zlo : z;
  const double g_mode = at_boundary ? g : 0.0;
  const double fmax = f(z_mode);
  const double sd = 1.0 / std::sqrt(-d2);
  const double target = fmax - kPostDrop;

  // Window end on one side: Newton on f = target.  From either side of the
  // root the concave f makes Newton land beyond it and then converge back.
  auto window_end = [&](double dir, double ze) {
    for (int it = 0; it < 60; ++it) {
      double dd;
      const double diff = f(ze) - target;
      if (std::fabs(diff) < 1e-3) break;
      const double slope = df(ze, dd);
      double zn = (slope * dir < 0.0) ? ze - diff / slope : ze + dir * (ze - z_mode);
      if ((zn - z_mode) * dir <= 0.0) zn = z_mode + 0.5 * (ze - z_mode);
      if (dir < 0.0 && zn <= zlo) {
        if (f(zlo) >= target) return zlo;
        zn = 0.5 * (ze + zlo);
      }
      ze = zn;
    }
    return ze;
  };
  double reach = kLaplaceReach * sd;
  if (g_mode < 0.0) reach = std::min(reach, kPostDrop / -g_mode);
  double zR = z_mode + reach;
  if (f(zR) > target + 2.0) zR = window_end(1.0, zR);
  double zL = z_mode;
  if (z_mode > zlo) {
    zL = std::max(zlo, z_mode - kLaplaceReach * sd);
    if (zL > zlo && f(zL) > target + 2.0) zL = window_end(-1.0, zL);
  }

  auto panel = [&](double a, double b, int n) {
    if (!(b > a)) return 0.0;
    const GLRule& gl = gl_get_rule(n);
    const double mid = 0.5 * (a + b), half = 0.5 * (b - a);
    double acc = 0.0;
    for (int j = 0; j < n; ++j)
      acc += gl.w[j] * std::exp(f(mid + half * gl.x[j]) - fmax);
    return half * acc;
  };
  auto side = [&](double to) {
    const double dir = (to >= z_mode) ? 1.0 : -1.0;
    const double len = dir * (to - z_mode);
    if (!(len > 0.0)) return 0.0;
    if (len <= first_panel) {
      const int n = std::min(nodes, 8 + 3 * static_cast<int>(std::ceil(len)));
      return panel(std::min(z_mode, to), std::max(z_mode, to), n);
    }
    double total = 0.0, from = z_mode, w = first_panel;
    while (dir * (to - from) > 0.0) {
      const double next = (dir * (to - from) > w) ? from + dir * w : to;
      total += panel(std::min(from, next), std::max(from, next), nodes);
      from = next;
      w *= 2.0;
    }
    return total;
  };
  return fmax + std::log(side(zL) + side(zR));
}

// ---- functions of the distance r -----------------------------------------

inline double w_of(const Law& L, double r) {
  return (L.sv2 * r + L.mu * L.s2) / (L.s * L.sv * L.sqrtQ);
}
inline double log_G(const Law& L, double r) {
  return log_phi_std((r - L.mut) / L.S) - L.log_S;
}
// Posterior drift mean given the fixed-time position r.
inline double post_mean(const Law& L, double r) {
  return (L.sv2 * r + L.mu * L.s2) / L.Q;
}
inline bool bvn_route(const Law& L) {
  return L.pos && L.mu >= 0.0 && L.rho < kBvnRhoMax;
}

// ---- initial values (log scale) --------------------------------------------
// Posterior forms: with V* ~ N(post_mean, (s sv)^2/Q) and z its standardised
// value (V* > 0  <=>  z > -w),
//   T = c G E[R(xT - rho z)],  U = c G E[R(-xT + rho z)],  H = c G E[R(xH + rho z)],
// xT = (r - post_mean t)/c, xH = (r + post_mean t)/c.

inline double seed_log_H(const Law& L, double r) {
  const double ell = 2.0 * r * L.mu / L.s2 + 2.0 * r * r * L.sv2 / (L.s2 * L.s2);
  const double mup = L.mu + 2.0 * r * L.sv2 / L.s2;
  const double h2 = (-mup * L.t - r) / L.S;
  if (!L.pos) return ell + log_Phi(h2);
  if (bvn_route(L) && h2 >= kBvnMinArg) {
    const double P = fast_norm_phi(h2);
    {
      const double q = bvn(L, -mup / L.sv, h2) / P;
      if (q >= 0.0 && q < kBvnMaxRatio) {
        RDMSWTN_SPV_COUNT(seeds_bvn, 1);
        return ell + std::log(P) + std::log1p(-q);
      }
    }
  }
  const double xH = (r + post_mean(L, r) * L.t) / L.c;
  RDMSWTN_SPV_COUNT(seeds_post, 1);
  return L.log_c + log_G(L, r) + log_posterior_mills(xH, L.rho, -w_of(L, r));
}

// Posterior-form T (want_T) or U.  The representation whose Mills argument is
// positive over the posterior bulk is evaluated; its complement in Z is taken
// by subtraction only when the complement is the larger part.
inline double post_log_TU(const Law& L, double r, bool want_T) {
  const double pre = L.log_c + log_G(L, r);
  const double zlo = -w_of(L, r);
  const double xT = (r - post_mean(L, r) * L.t) / L.c;
  const bool native_T = (xT >= 0.0);
  auto form = [&](bool T) {
    return pre + (T ? log_posterior_mills(xT, -L.rho, zlo)
                    : log_posterior_mills(-xT, L.rho, zlo));
  };
  if (native_T == want_T) return form(want_T);
  const double other = form(!want_T);
  if (other < L.logZ - M_LN2) return log_diff_exp(L.logZ, other);
  return form(want_T);
}

inline double seed_log_T(const Law& L, double r) {
  const double h1 = (L.mut - r) / L.S;
  if (!L.pos) return log_Phi(h1);
  if (bvn_route(L) && h1 >= kBvnMinArg) {
    RDMSWTN_SPV_COUNT(seeds_bvn, 1);
    return std::log(bvn(L, L.mu / L.sv, h1));
  }
  RDMSWTN_SPV_COUNT(seeds_post, 1);
  return post_log_TU(L, r, true);
}

inline double seed_log_U(const Law& L, double r) {
  const double h1 = (L.mut - r) / L.S;
  if (!L.pos) return log_Phi(-h1);
  if (bvn_route(L) && -h1 >= kBvnMinArg) {
    const double P = fast_norm_phi(-h1);
    {
      const double q = bvn(L, -L.mu / L.sv, -h1) / P;
      if (q >= 0.0 && q < kBvnMaxRatio) {
        RDMSWTN_SPV_COUNT(seeds_bvn, 1);
        return std::log(P) + std::log1p(-q);
      }
    }
  }
  RDMSWTN_SPV_COUNT(seeds_post, 1);
  return post_log_TU(L, r, false);
}

// ---- the recurrence ---------------------------------------------------------

// Gaussian factors at an expansion point, on the log scale, with their
// log-derivatives in r.  J = G Phi(w), K = G phi(w), M = delta Phi(-r/c),
// N = delta phi(r/c); full-Gaussian drift has J = G and no K, M, N.
struct Aux {
  double lJ, lK, lM, lN, w0;
  double rG, rJ, rK, rM, rN, rH;
};

inline Aux aux_at(const Law& L, double r0) {
  Aux x;
  const double lG = log_G(L, r0);
  x.rG = -(r0 - L.mut) / L.S2;
  x.rH = L.a + L.beta * r0;
  if (L.pos) {
    x.w0 = w_of(L, r0);
    x.lJ = lG + log_Phi(x.w0);
    x.lK = lG + log_phi_std(x.w0);
    x.rJ = x.rG + L.u * normal_hazard(-x.w0);
    x.rK = x.rG - L.u * x.w0;
    const double zc = r0 / L.c;
    x.lM = L.log_delta + log_Phi(-zc);
    x.lN = L.log_delta + log_phi_std(zc);
    x.rM = -normal_hazard(zc) / L.c;
    x.rN = -r0 / L.c2;
  } else {
    x.w0 = 0.0;
    x.lJ = lG;
    x.lK = x.lM = x.lN = R_NegInf;
    x.rJ = x.rG;
    x.rK = x.rM = x.rN = 0.0;
  }
  return x;
}

// Step-size limit contributed by one Gaussian factor of log weight lw (its
// reach into the integrated quantities, relative to the carried scale), log
// rate `rate`, curvature `curv`, travelling in direction dir.  A factor that
// is negligible and shrinking is dropped; a negligible growing one either
// stays negligible over the step or is resolved as an active factor.
struct FactorPlan { double hmax; double h_inactive; };

inline FactorPlan plan_factor(double lw, double rate, double curv, double dir) {
  FactorPlan p{R_PosInf, R_PosInf};
  if (lw == R_NegInf) return p;
  const double along = dir * rate;
  const double slack = std::max(0.0, -lw);
  double h_active = R_PosInf;
  if (along > 0.0) h_active = kGrowBudget / along;
  else if (along < 0.0) h_active = (kDecayBudget + slack) / (-along);
  if (curv > 0.0) h_active = std::min(h_active, std::sqrt((kCurvBudget + 2.0 * slack) / curv));
  if (lw < kNegligible) {
    if (along <= 0.0) return p;               // negligible and not growing: dropped
    p.h_inactive = (kNegligible - lw) / along;
    p.hmax = std::max(p.h_inactive, h_active);
    return p;
  }
  p.h_inactive = 0.0;
  p.hmax = h_active;
  return p;
}

enum class Carry { TH, H, U };

// Integrate the carried quantities over the range, one adaptive step at a
// time.  TH and H travel downwards from r_from (seeds at r_from); U travels
// upwards.  Returns log of the integral of T + H, H or U.  Full-Gaussian drift
// re-seeds from closed forms at every step; positive truncation carries the
// values (and re-seeds H where its homogeneous mode would grow downwards,
// which needs m < 0).
template <Carry C>
inline double sweep(const Law& L, double r_from, double r_to,
                    double log_seed1, double log_seed2) {
  constexpr bool hasT = (C == Carry::TH);
  constexpr bool hasH = (C != Carry::U);
  constexpr bool hasU = (C == Carry::U);
  const double dir = hasU ? 1.0 : -1.0;
  const double one_k = hasH ? 1.0 + L.k : 1.0;

  double kappa, X1, X2 = 0.0;   // carried values are exp(kappa) * (X1, X2)
  auto set_scale = [&](double l1, double l2) {
    if (hasT) {
      kappa = log_sum_exp(l1, l2);
      X1 = std::exp(l1 - kappa);
      X2 = std::exp(l2 - kappa);
    } else {
      kappa = l1;
      X1 = 1.0;
    }
  };
  set_scale(log_seed1, log_seed2);

  double r0 = r_from;
  double logI = R_NegInf;
  for (int step = 0; step < kMaxSteps; ++step) {
    const double rem = dir < 0.0 ? r0 - r_to : r_to - r0;
    if (!(rem > 1e-15 * (1.0 + std::fabs(r_to)))) break;

    if (step > 0) {
      if (!L.pos) {
        if (hasT) set_scale(seed_log_T(L, r0), seed_log_H(L, r0));
        else if (hasH) set_scale(seed_log_H(L, r0), 0.0);
        else set_scale(seed_log_U(L, r0), 0.0);
      } else if (hasH && L.a + L.beta * r0 < 0.0) {
        if (hasT) set_scale(kappa + std::log(X1), seed_log_H(L, r0));
        else set_scale(seed_log_H(L, r0), 0.0);
      }
    }
    if (!R_FINITE(kappa)) return kappa == R_NegInf ? R_NegInf : R_NaN;

    const Aux ax = aux_at(L, r0);
    const double lwJ = ax.lJ - kappa + std::log(one_k * L.S);
    const double lwK = ax.lK - kappa + std::log(one_k * L.u * L.S2);
    const double lwM = hasH ? ax.lM - kappa + L.log_S : R_NegInf;
    const double lwN = hasH ? ax.lN - kappa + 2.0 * L.log_S - L.log_c : R_NegInf;

    const FactorPlan pJ = plan_factor(lwJ, ax.rJ, 1.0 / L.S2, dir);
    const FactorPlan pK = plan_factor(lwK, ax.rK, 1.0 / L.S2 + L.u * L.u, dir);
    const FactorPlan pM = plan_factor(lwM, ax.rM, 1.0 / L.c2, dir);
    const FactorPlan pN = plan_factor(lwN, ax.rN, 1.0 / L.c2, dir);
    double habs = std::min({rem, pJ.hmax, pK.hmax, pM.hmax, pN.hmax});
    if (hasH) {
      const double along = dir * ax.rH;
      if (along < 0.0) habs = std::min(habs, kStiffBudget / (-along));
      else if (along > 0.0) habs = std::min(habs, kGrowBudget / along);
      if (L.beta > 0.0) habs = std::min(habs, std::sqrt(2.0 * kStiffBudget / L.beta));
    }
    if (!(habs > 0.0)) return R_NaN;

    for (int attempt = 0; attempt < 40; ++attempt) {
      const double h = dir * habs;
      const double J0 = (habs > pJ.h_inactive) ? std::exp(ax.lJ - kappa) : 0.0;
      const double K0 = (habs > pK.h_inactive) ? std::exp(ax.lK - kappa) : 0.0;
      const double M0 = (hasH && habs > pM.h_inactive) ? std::exp(ax.lM - kappa) : 0.0;
      const double N0 = (hasH && habs > pN.h_inactive) ? std::exp(ax.lN - kappa) : 0.0;

      const double a0 = h * ax.rH, a1 = L.beta * h * h;
      const double l0 = h * ax.rG, l1 = -h * h / L.S2;
      const double hu = h * L.u;
      const double m0 = l0 - hu * ax.w0, m1 = l1 - hu * hu;
      const double e0 = h * ax.rN, e1 = -h * h / L.c2;
      const double hc = h / L.c, kh = L.k * h;
      const double wJ = habs * one_k, wK = habs * habs * one_k * L.u;
      const double wN = habs * habs / L.c;
      const double peak = std::max({hasH ? std::fabs(a0) : 0.0, std::fabs(l0),
                                    std::fabs(m0), std::fabs(e0)});
      const int nmin = static_cast<int>(std::ceil(peak)) + 3;

      double Yn = X1, Hn = hasT ? X2 : X1;   // Y: T or U (unused for C == H)
      double Jn = J0, Kn = K0, Mn = M0, Nn = N0;
      double Hp = 0.0, Jp = 0.0, Kp = 0.0, Np = 0.0;
      double sum, Yend = Yn, Hend = Hn;
      if (hasT) sum = Yn + Hn;
      else if (hasH) sum = Hn;
      else sum = Yn;
      int quiet = 0;
      bool converged = false;
      for (int n = 0; n < kMaxTerms; ++n) {
        const double inv = 1.0 / (n + 1);
        const double nY = (hasU ? h : -h) * Jn * inv;
        const double nH = hasH ? (a0 * Hn + a1 * Hp - kh * Jn + h * Mn) * inv : 0.0;
        const double nJ = (l0 * Jn + l1 * Jp + hu * Kn) * inv;
        const double nK = (m0 * Kn + m1 * Kp) * inv;
        const double nM = hasH ? -hc * Nn * inv : 0.0;
        const double nN = hasH ? (e0 * Nn + e1 * Np) * inv : 0.0;
        Hp = Hn; Jp = Jn; Kp = Kn; Np = Nn;
        Yn = nY; Hn = nH; Jn = nJ; Kn = nK; Mn = nM; Nn = nN;
        const double inv2 = 1.0 / (n + 2);
        double lead;
        if (hasT) { sum += (Yn + Hn) * inv2; Yend += Yn; Hend += Hn; lead = std::fabs(Yn) + std::fabs(Hn); }
        else if (hasH) { sum += Hn * inv2; Hend += Hn; lead = std::fabs(Hn); }
        else { sum += Yn * inv2; Yend += Yn; lead = std::fabs(Yn); }
        const double mag = lead + wJ * std::fabs(Jn) + wK * std::fabs(Kn) +
                           habs * std::fabs(Mn) + wN * std::fabs(Nn);
        if (n + 1 >= nmin && mag <= kTermTol * (1.0 + std::fabs(sum))) {
          if (++quiet >= 2) { converged = true; RDMSWTN_SPV_COUNT(terms, n + 1); break; }
        } else {
          quiet = 0;
        }
      }
      if (!converged) { habs *= 0.5; continue; }
      if (!(sum > 0.0)) return R_NaN;

      logI = log_sum_exp(logI, kappa + std::log(habs * sum));
      if (hasT) {
        const double tot = Yend + Hend;
        if (!(tot > 0.0)) return R_NaN;
        kappa += std::log(tot);
        X1 = Yend / tot;
        X2 = Hend / tot;
      } else {
        const double v = hasH ? Hend : Yend;
        if (!(v > 0.0)) return R_NaN;
        kappa += std::log(v);
        X1 = 1.0;
      }
      RDMSWTN_SPV_COUNT(steps, 1);
      r0 += h;
      break;
    }
  }
  return logI;
}

// Skip bounds.  A Gaussian proxy for the fixed-time position locates where the
// integrand has fallen kSkipDrop below its dominant end.
inline double cdf_top(const Law& L, double r_lo, double r_hi) {
  auto P = [&](double r) { return log_Phi((L.mut - r) / L.S); };
  const double p_lo = P(r_lo);
  if (!(p_lo - P(r_hi) > kSkipDrop)) return r_hi;
  double lo = r_lo, hi = r_hi;
  for (int it = 0; it < 80; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (p_lo - P(mid) > kSkipDrop) hi = mid; else lo = mid;
  }
  return hi;
}

inline double surv_bottom(const Law& L, double r_lo, double r_hi) {
  auto P = [&](double r) { return log_Phi((r - L.mut) / L.S); };
  const double p_hi = P(r_hi);
  if (!(p_hi - P(r_lo) > kSkipDrop)) return r_lo;
  double lo = r_lo, hi = r_hi;
  for (int it = 0; it < 80; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (p_hi - P(mid) > kSkipDrop) lo = mid; else hi = mid;
  }
  return lo;
}

// log of the integral of T + H over [r_lo, r_hi].
inline double integral_F(const Law& L, double r_lo, double r_hi) {
  double r_top = cdf_top(L, r_lo, r_hi);
  for (;;) {
    const double lT = seed_log_T(L, r_top), lH = seed_log_H(L, r_top);
    const double logI = sweep<Carry::TH>(L, r_top, r_lo, lT, lH);
    if (r_top >= r_hi || ISNAN(logI)) return logI;
    if (log_sum_exp(lT, lH) + std::log(r_hi - r_top) <= std::log(kSkipTol) + logI)
      return logI;
    r_top = r_hi;
  }
}

// log of the integral of U - H over [r_lo, r_hi]; NaN when the difference is
// below the resolution of the two integrals.
inline double integral_S(const Law& L, double r_lo, double r_hi) {
  double r_bot = surv_bottom(L, r_lo, r_hi);
  const double lH_top = seed_log_H(L, r_hi);
  for (;;) {
    const double lU = seed_log_U(L, r_bot);
    const double logIU = sweep<Carry::U>(L, r_bot, r_hi, lU, 0.0);
    const double logIH = sweep<Carry::H>(L, r_hi, r_bot, lH_top, 0.0);
    if (ISNAN(logIU) || ISNAN(logIH) || !(logIU > logIH)) return R_NaN;
    const double logS = log_diff_exp(logIU, logIH);
    if (r_bot <= r_lo) return logS;
    if (lU + std::log(r_bot - r_lo) <= std::log(kSkipTol) + logS) return logS;
    r_bot = r_lo;
  }
}

// log of the integral of d J(d) over [r_lo, r_hi]; the decision density
// averaged over the distance range is that integral / (A t Z).  d J(d) is
// log-concave for d > 0, so the range is split at its mode and each part is
// integrated towards the mode, along which the integrand grows.  J and K are
// re-seeded in closed form at the start of every step.
inline double integral_dJ(const Law& L, double r_lo, double r_hi) {
  auto slope = [&](double r, double& d_slope) {
    const double w = w_of(L, r);
    double rJ = -(r - L.mut) / L.S2, drJ = -1.0 / L.S2;
    if (L.pos) {
      const double lam = normal_hazard(-w);   // phi(w) / Phi(w)
      rJ += L.u * lam;
      drJ -= L.u * L.u * lam * (lam + w);
    }
    d_slope = drJ - 1.0 / (r * r);
    return rJ + 1.0 / r;
  };
  double r_mode;
  double ds;
  if (slope(r_hi, ds) >= 0.0) {
    r_mode = r_hi;
  } else if (slope(r_lo, ds) <= 0.0) {
    r_mode = r_lo;
  } else {
    double lo = r_lo, hi = r_hi, r = 0.5 * (r_lo + r_hi);
    for (int it = 0; it < 100; ++it) {
      const double g = slope(r, ds);
      if (g > 0.0) lo = r; else hi = r;
      double rn = r - g / ds;
      if (!(rn > lo && rn < hi)) rn = 0.5 * (lo + hi);
      if (std::fabs(rn - r) <= 1e-10 * (1.0 + std::fabs(r))) { r = rn; break; }
      r = rn;
    }
    r_mode = r;
  }

  auto part = [&](double r_from, double r_to) {
    const double dir = (r_to >= r_from) ? 1.0 : -1.0;
    double r0 = r_from, logI = R_NegInf;
    for (int step = 0; step < kMaxSteps; ++step) {
      const double rem = dir * (r_to - r0);
      if (!(rem > 1e-15 * (1.0 + std::fabs(r_to)))) break;
      const Aux ax = aux_at(L, r0);
      if (ax.lJ == R_NegInf) return R_NegInf;
      const double lwK = L.pos ? ax.lK - ax.lJ + std::log(L.u * L.S2) : R_NegInf;
      const FactorPlan pJ = plan_factor(0.0, ax.rJ, 1.0 / L.S2, dir);
      const FactorPlan pK = plan_factor(lwK, ax.rK, 1.0 / L.S2 + L.u * L.u, dir);
      double habs = std::min({rem, pJ.hmax, pK.hmax});
      if (!(habs > 0.0)) return R_NaN;
      bool done = false;
      for (int attempt = 0; attempt < 40 && !done; ++attempt) {
        const double h = dir * habs;
        const double K0 = (habs > pK.h_inactive) ? std::exp(ax.lK - ax.lJ) : 0.0;
        const double l0 = h * ax.rG, l1 = -h * h / L.S2;
        const double hu = h * L.u;
        const double m0 = l0 - hu * ax.w0, m1 = l1 - hu * hu;
        const double wK = habs * habs * L.u * (std::fabs(r0) + habs);
        const int nmin = static_cast<int>(std::ceil(std::max(std::fabs(l0), std::fabs(m0)))) + 3;
        double Jn = 1.0, Kn = K0, Jp = 0.0, Kp = 0.0;
        double sum = r0 + 0.5 * h;
        int quiet = 0;
        for (int n = 0; n < kMaxTerms; ++n) {
          const double inv = 1.0 / (n + 1);
          const double nJ = (l0 * Jn + l1 * Jp + hu * Kn) * inv;
          const double nK = (m0 * Kn + m1 * Kp) * inv;
          Jp = Jn; Kp = Kn; Jn = nJ; Kn = nK;
          sum += Jn * (r0 / (n + 2) + h / (n + 3));
          const double mag = std::fabs(Jn) * (std::fabs(r0) + habs) + wK * std::fabs(Kn);
          if (n + 1 >= nmin && mag <= kTermTol * std::fabs(sum)) {
            if (++quiet >= 2) { done = true; break; }
          } else {
            quiet = 0;
          }
        }
        if (!done) { habs *= 0.5; continue; }
        if (!(sum > 0.0)) return R_NaN;
        logI = log_sum_exp(logI, ax.lJ + std::log(habs * sum));
        r0 += h;
      }
      if (!done) return R_NaN;
    }
    return logI;
  };
  return log_sum_exp(part(r_lo, r_mode), part(r_hi, r_mode));
}

// Full-Gaussian eventual hit and miss masses at distance r and averaged over
// [r_lo, r_hi].  With P = Phi(-mu/sv) and Hinf(r) = E[exp(2rV/s^2); V < 0] =
// exp(ell) Phi(-mu_p/sv):
//   Hinf' = (a + beta r) Hinf - delta,
//   D = P - Hinf,  D' = (a + beta r)(D - P) + delta,
// delta = 2 sv phi(mu/sv)/s^2.  Both are re-seeded at every step.
inline void inf_masses(double mu, double sv, double s, double r_lo, double r_hi,
                       double& log_hit, double& log_miss) {
  const double s2 = s * s, sv2 = sv * sv;
  const double a = 2.0 * mu / s2, beta = 4.0 * sv2 / (s2 * s2);
  const double z = mu / sv;
  const double logP = log_Phi(-z);
  const double delta = 2.0 * sv / s2 * std::exp(log_phi_std(z));
  // delta - a P = (2 sv/s^2)(phi(z) - z Phi(-z)), formed without cancellation.
  const double delta_minus_aP = (z >= 0.0)
    ? 2.0 * sv / s2 * std::exp(log_phi_std(z)) * mills_ratio_decrement(z)
    : delta - a * std::exp(logP);
  auto log_hinf = [&](double r) {
    const double ell = 2.0 * r * mu / s2 + 2.0 * r * r * sv2 / (s2 * s2);
    return ell + log_Phi(-(mu + 2.0 * r * sv2 / s2) / sv);
  };
  auto log_d = [&](double r) {
    const double lh = log_hinf(r);
    return (lh < logP) ? logP + log1m_exp(lh - logP) : R_NegInf;
  };

  const double span = r_hi - r_lo;
  if (!(span > 1e-12 * std::max(1.0, r_hi))) {
    log_hit = log_sum_exp(log_Phi(z), log_hinf(r_hi));
    log_miss = log_d(r_hi);
    return;
  }

  double logIH = R_NegInf, logID = R_NegInf;
  double r0 = r_hi;
  for (int step = 0; step < kMaxSteps; ++step) {
    const double rem = r0 - r_lo;
    if (!(rem > 1e-15 * (1.0 + std::fabs(r_lo)))) break;
    const double rate = a + beta * r0;
    double habs = rem;
    if (rate != 0.0) habs = std::min(habs, kDecayBudget / std::fabs(rate));
    habs = std::min(habs, std::sqrt(kCurvBudget / beta));
    const double h = -habs;
    const double a0 = h * rate, a1 = beta * h * h;
    const double peak = std::fabs(a0);
    const int nmin = static_cast<int>(std::ceil(peak)) + 3;

    // Hinf series, scaled by its value at r0.
    const double kH = log_hinf(r0);
    {
      double Hn = 1.0, Hp = 0.0, sum = 1.0;
      int quiet = 0;
      for (int n = 0; n < kMaxTerms; ++n) {
        double nH = a0 * Hn + a1 * Hp;
        if (n == 0) nH -= h * delta * std::exp(-kH);
        nH /= (n + 1);
        Hp = Hn; Hn = nH;
        sum += Hn / (n + 2);
        if (n + 1 >= nmin && std::fabs(Hn) <= kTermTol * (1.0 + std::fabs(sum))) {
          if (++quiet >= 2) break;
        } else quiet = 0;
      }
      if (sum > 0.0) logIH = log_sum_exp(logIH, kH + std::log(habs * sum));
    }
    // D series, scaled by its value at r0.
    const double kD = log_d(r0);
    if (kD > R_NegInf) {
      const double f0 = h * (delta_minus_aP - beta * r0 * std::exp(logP)) * std::exp(-kD);
      const double f1 = -beta * h * h * std::exp(logP - kD);
      double Dn = 1.0, Dp = 0.0, sum = 1.0;
      int quiet = 0;
      for (int n = 0; n < kMaxTerms; ++n) {
        double nD = a0 * Dn + a1 * Dp;
        if (n == 0) nD += f0;
        if (n == 1) nD += f1;
        nD /= (n + 1);
        Dp = Dn; Dn = nD;
        sum += Dn / (n + 2);
        if (n + 1 >= std::max(nmin, 3) &&
            std::fabs(Dn) <= kTermTol * (1.0 + std::fabs(sum))) {
          if (++quiet >= 2) break;
        } else quiet = 0;
      }
      if (sum > 0.0) logID = log_sum_exp(logID, kD + std::log(habs * sum));
    }
    r0 += h;
  }
  const double log_span = std::log(span);
  log_hit = log_sum_exp(log_Phi(z), logIH - log_span);
  log_miss = logID - log_span;
}

// ---- public entry points ----------------------------------------------------

// log F(t) of the no-clock decision time: t is decision time (RT - t0), b the
// upper distance, A the uniform range (A = 0 is a fixed distance), s > 0 the
// diffusion coefficient, sv > 0 the drift SD.
inline double log_cdf(double t, double mu, double b, double A, double s,
                      double sv, bool pos) {
  if (ISNAN(t) || ISNAN(mu) || ISNAN(b) || ISNAN(A) || ISNAN(s) || ISNAN(sv))
    return R_NaN;
  if (!(s > 0.0) || !(sv > 0.0) || A < 0.0) return R_NaN;
  if (!(t > 0.0)) return R_NegInf;
  if (!(b > 0.0)) return 0.0;
  const double r_hi = b, r_lo = std::max(0.0, b - A);
  if (!R_FINITE(t)) {
    if (pos) return 0.0;
    double lh, lm;
    inf_masses(mu, sv, s, r_lo, r_hi, lh, lm);
    return std::fmin(lh, 0.0);
  }
  const Law L = make_law(t, mu, s, sv, pos);
  const double span = r_hi - r_lo;
  double out;
  if (!(span > 1e-12 * std::max(1.0, r_hi))) {
    out = log_sum_exp(seed_log_T(L, r_hi), seed_log_H(L, r_hi)) - L.logZ;
  } else {
    out = integral_F(L, r_lo, r_hi) - std::log(span) - L.logZ;
  }
  return ISNAN(out) ? out : std::fmin(out, 0.0);
}

// log f(t) of the no-clock decision time for a uniform distance range
// (A > 0): the distance average of the fixed-distance density (d/t) J(d)/Z,
// i.e. (1/(A t Z)) E[Y 1{r_lo <= Y <= r_hi, V > 0}] for the fixed-time
// position Y ~ N(mu t, S^2).  That first moment is a Gaussian-rectangle
// expression; it is used unless its terms cancel below kDensityKeep of their
// scale (bivariate-normal probabilities carry absolute, not relative, error),
// in which case the distance recurrence gives the integral directly.
constexpr double kDensityKeep = 1e-4;

inline double log_density(double t, double mu, double b, double A, double s,
                          double sv, bool pos) {
  if (ISNAN(t) || ISNAN(mu) || ISNAN(b) || ISNAN(A) || ISNAN(s) || ISNAN(sv))
    return R_NaN;
  if (!(s > 0.0) || !(sv > 0.0) || !(A > 0.0)) return R_NaN;
  if (!(t > 0.0) || !R_FINITE(t) || !(b > 0.0)) return R_NegInf;
  const double r_hi = b, r_lo = std::max(0.0, b - A);
  if (!(r_hi > r_lo)) return R_NegInf;
  const Law L = make_law(t, mu, s, sv, pos);

  const double a_lo = (r_lo - L.mut) / L.S, a_hi = (r_hi - L.mut) / L.S;
  const double phi_lo = std::exp(log_phi_std(a_lo));
  const double phi_hi = std::exp(log_phi_std(a_hi));
  // Phi(a_hi) - Phi(a_lo) from whichever tail keeps it relative.
  const double dPhi = (a_lo >= 0.0)
    ? std::exp(log_Phi(-a_lo)) * -std::expm1(log_Phi(-a_hi) - log_Phi(-a_lo))
    : std::exp(log_Phi(a_hi)) * -std::expm1(log_Phi(a_lo) - log_Phi(a_hi));
  double num, scale;
  if (!pos) {
    num = L.mut * dPhi + L.S * (phi_lo - phi_hi);
    scale = std::fabs(L.mut) * dPhi + L.S * (phi_lo + phi_hi);
  } else {
    // (Y, V) standardised have correlation rho; V > 0 is Z_V > gamma.
    const double gamma = -mu / sv;
    const double sq = std::sqrt(std::fmax(1e-300, 1.0 - L.rho * L.rho));
    const double phi_g = std::exp(log_phi_std(gamma));
    auto Phi2 = [&](double x) {
      return bvn_route(L) ? bvn(L, x, gamma) : norm_cdf_2d(x, gamma, L.rho);
    };
    auto Gm = [&](double y, double phi_y) {
      return -phi_y * pnorm_std((L.rho * y - gamma) / sq) +
             L.rho * phi_g * pnorm_std((y - L.rho * gamma) / sq);
    };
    const double p_rect = dPhi - (Phi2(a_hi) - Phi2(a_lo));
    num = L.mut * p_rect + L.S * (Gm(a_hi, phi_hi) - Gm(a_lo, phi_lo));
    scale = std::fabs(L.mut) + L.S;
  }
  if (num > kDensityKeep * scale && R_FINITE(num))
    return std::log(num) - std::log(A) - std::log(t) - L.logZ;
  return integral_dJ(L, r_lo, r_hi) - std::log(A) - std::log(t) - L.logZ;
}

// log S(t) = log(1 - F(t)) with relative accuracy when F is close to one.
inline double log_surv(double t, double mu, double b, double A, double s,
                       double sv, bool pos) {
  if (ISNAN(t) || ISNAN(mu) || ISNAN(b) || ISNAN(A) || ISNAN(s) || ISNAN(sv))
    return R_NaN;
  if (!(s > 0.0) || !(sv > 0.0) || A < 0.0) return R_NaN;
  if (!(t > 0.0)) return 0.0;
  if (!(b > 0.0)) return R_NegInf;
  const double r_hi = b, r_lo = std::max(0.0, b - A);
  if (!R_FINITE(t)) {
    if (pos) return R_NegInf;
    double lh, lm;
    inf_masses(mu, sv, s, r_lo, r_hi, lh, lm);
    return std::fmin(lm, 0.0);
  }
  const Law L = make_law(t, mu, s, sv, pos);
  const double span = r_hi - r_lo;
  const bool point = !(span > 1e-12 * std::max(1.0, r_hi));
  auto from_cdf = [](double lF) {
    if (ISNAN(lF)) return lF;
    return (lF >= 0.0) ? R_NegInf : log1m_exp(lF);
  };

  // Full-Gaussian F at the top distance (exact for that law, a proxy under
  // truncation) decides whether the survivor is the small side.
  const double h1 = (L.mut - r_hi) / L.S;
  const double ell = 2.0 * r_hi * mu / L.s2 + 2.0 * r_hi * r_hi * L.sv2 / (L.s2 * L.s2);
  const double h2 = (-(mu + 2.0 * r_hi * L.sv2 / L.s2) * t - r_hi) / L.S;
  const double proxy = log_sum_exp(log_Phi(h1), ell + log_Phi(h2));

  if (point) {
    // S = (U - H)/Z when the survivor is the small side (U > Z/2 whenever
    // S > 1/2, so a misjudged proxy costs at most a factor two), otherwise
    // 1 - (T + H)/Z.
    const double lH = seed_log_H(L, r_hi);
    if (proxy > -M_LN2) {
      const double lU = seed_log_U(L, r_hi);
      if (lU > lH) return std::fmin(log_diff_exp(lU, lH) - L.logZ, 0.0);
      return from_cdf(log_sum_exp(seed_log_T(L, r_hi), lH) - L.logZ);
    }
    const double lF = log_sum_exp(seed_log_T(L, r_hi), lH) - L.logZ;
    if (!(lF > -M_LN2)) return from_cdf(lF);
    const double lU = seed_log_U(L, r_hi);
    if (lU > lH) return std::fmin(log_diff_exp(lU, lH) - L.logZ, 0.0);
    return from_cdf(lF);
  }

  const double log_span = std::log(span);
  double lF = R_NaN;
  if (!(proxy > -M_LN2)) {
    lF = std::fmin(integral_F(L, r_lo, r_hi) - log_span - L.logZ, 0.0);
    if (!(lF > std::log(0.9))) return from_cdf(lF);
  }
  const double lIS = integral_S(L, r_lo, r_hi);
  if (!ISNAN(lIS)) return std::fmin(lIS - log_span - L.logZ, 0.0);
  if (ISNAN(lF)) lF = std::fmin(integral_F(L, r_lo, r_hi) - log_span - L.logZ, 0.0);
  return from_cdf(lF);
}

// Full-Gaussian eventual hit mass and never-finish mass (t = Inf).
inline double log_hit_mass_fg(double mu, double sv, double s, double b, double A) {
  double lh, lm;
  inf_masses(mu, sv, s, std::max(0.0, b - A), b, lh, lm);
  return std::fmin(lh, 0.0);
}
inline double log_miss_mass_fg(double mu, double sv, double s, double b, double A) {
  double lh, lm;
  inf_masses(mu, sv, s, std::max(0.0, b - A), b, lh, lm);
  return std::fmin(lm, 0.0);
}

}  // namespace rdmswtn_spv

#endif

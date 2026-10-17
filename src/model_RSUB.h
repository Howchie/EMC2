#ifndef EMC2_MODEL_RSUB_H
#define EMC2_MODEL_RSUB_H

// Racing subordinators (Math/RDMSWTN_Derivation.tex).
//
// Accumulator evidence G(u) is a nondecreasing pure-jump Levy process with
// E[G(u)] = v u and Var[G(u)] = c v u, c = s^2, and the accumulator finishes
// at the first u with G(u) >= d.  Because G never decreases,
// F_d(u) = P(G(u) >= d), and with the distance uniform on [B, B + A] the
// averaged CDF, survivor and density are stop-loss differences of G(u) at
// the two ends (m = v u):
//
//   F = [H(B) - H(B+A)]/A,   H(d) = E(G - d)+,
//   S = [L(B+A) - L(B)]/A,   L(d) = E(d - G)+ = H(d) + d - m,
//   f = v [J(B) - J(B+A)]/A, J(d) = dH/dm,  Jc = 1 - J.
//
// Two evidence laws share this structure (namespaces ig and cp):
//
// ig, "dense": the inverse-Gaussian subordinator, G(u) the first time a
// Brownian motion with unit diffusion and drift 1/s reaches (v/s) u.  With
// z = s sqrt(d), alpha = (m - d)/z and beta = (m + d)/z,
//
//   S_d = Phi(-alpha) + phi(alpha) R(beta),
//   f_d = (2v/z) phi(alpha) [D(beta) + (m/z) R(beta)],
//   H = z phi(alpha) [D(-alpha) - D(beta)],  L = z phi(alpha) [D(alpha) - D(beta)],
//   Jc = phi(alpha) [R(alpha) - R(beta) + (beta - alpha) D(beta)],
//
// with R the standard normal Mills ratio, D = -R' = 1 - xR and E = -D' =
// R - xD.  Every bracket is an integral of D or E plus positive terms.  The
// finish time is (s/v) max_{0<=w<=d} (W_w + w/s).
//
// cp, "packets": compound Poisson evidence, events at rate 2v/c carrying
// exponential packets of mean c/2.  With a = 2m/c, b = 2d/c and
// Delta = N - K for independent N ~ Pois(a), K ~ Pois(b),
//
//   F_d = P(Delta >= 1),  f_d = (2v/c) P(Delta = 0),
//   H = (c/2) E[Delta+],  L = (c/2) E[Delta-],  J = P(Delta >= 0),
//
// with the Skellam law P(Delta = k) = e^{-(sqrt a - sqrt b)^2} (a/b)^{k/2}
// e^{-z} I_|k|(z), z = 2 sqrt(ab).  The finish time is Gamma(1 + K, 2v/c).
//
// The difference of two ends is taken in whichever complementary form (J or
// Jc, H or L) cancels least.  Where the ends still agree too closely for their
// difference to keep its digits, it is replaced by Gauss-Legendre quadrature
// of the point-start quantity over the start-point interval.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "utility_functions.h"
#include "composite_functions.h"
#include "wald_functions.h"
#include "gl_quad.h"
#include "model_OPR.h"

namespace rsub {

// Two ends whose ratio exceeds this are differenced by quadrature instead;
// below it a subtraction loses at most two digits.
constexpr double kCancel = 0.99;
// The same for Mills-ratio differences inside an end, whose operands already
// carry the cancellation of D = 1 - xR and E = R - xD.
constexpr double kMillsCancel = 0.9;
constexpr int kGLNodes = 8;

// v = 0 is valid: the accumulator never finishes.  cv = 0 is a fixed rate.
inline bool pars_ok(double v, double B, double A, double s, double cv = 0.0) {
  return R_FINITE(v) && v >= 0.0 && R_FINITE(B) && B >= 0.0 &&
         R_FINITE(A) && A >= 0.0 && B + A > 0.0 && R_FINITE(s) && s > 0.0 &&
         R_FINITE(cv) && cv >= 0.0;
}

inline double log_diff(double lx, double ly) {  // log(e^lx - e^ly), lx > ly
  return lx + std::log1p(-std::exp(ly - lx));
}

// min(x, 0) for a log probability, keeping NaN (std::fmin would drop it).
inline double clamp_log_prob(double x) { return x > 0.0 ? 0.0 : x; }

// log(e^lx - e^ly) for ly below lx, which keeps -Inf operands.
inline double log_less(double lx, double ly) {
  if (ly == R_NegInf || lx == R_NegInf) return lx;
  return lx + std::log1p(-std::exp(ly - lx));
}

inline double log_add(double lx, double ly) {
  if (lx < ly) std::swap(lx, ly);
  if (ly == R_NegInf) return lx;
  return lx + std::log1p(std::exp(ly - lx));
}

struct EndJ { double lJ, lJc; };
struct EndHL { double lH, lL; };

// The start-range differences in linear arithmetic, under the rules of the
// log-scale path below; NaN where that path has to settle them.  From the
// stop-loss ends, log S (want_S) or log F:
inline double lin_cdf_surv(double H1, double L1, double H2, double L2, double A, bool want_S) {
  const double S = (L2 - L1) / A, F = (H1 - H2) / A;
  if (L1 < kCancel * L2 && S < 0.5) return want_S ? std::log(S) : std::log1p(-S);
  if (H2 < kCancel * H1 && F < 0.5) return want_S ? std::log1p(-F) : std::log(F);
  return R_NaN;
}

// From the J ends (NaN for a form not formed), log of J(B) - J(B+A) =
// Jc(B+A) - Jc(B), in the form whose larger term is smaller when both are
// at hand and in the preferred one otherwise; a form missing at either end
// leaves the other untried.
inline double lin_J_diff(double J1, double Jc1, double J2, double Jc2, bool prefer_Jc) {
  const bool both = !ISNAN(J1) && !ISNAN(J2) && !ISNAN(Jc1) && !ISNAN(Jc2);
  const bool use_Jc = both ? Jc2 < J1 : prefer_Jc;
  if (use_Jc ? Jc1 < kCancel * Jc2 : J2 < kCancel * J1)
    return std::log(use_Jc ? Jc2 - Jc1 : J1 - J2);
  if (both && (use_Jc ? J2 < kCancel * J1 : Jc1 < kCancel * Jc2))
    return std::log(use_Jc ? J1 - J2 : Jc2 - Jc1);
  return R_NaN;
}

// ===========================================================================
// Inverse-Gaussian subordinator
// ===========================================================================
namespace ig {

// The Mills ratio R(x) = Phi(-x)/phi(x) and its successive negative
// derivatives D = 1 - xR and E = R - xD.  Below x = 1.5 they are formed
// directly, which costs a few ulps in E.  Above it the subtractions are
// exact in the tails of the Laplace continued fraction R = 1/(x + c1),
// c_k = k/(x + c_{k+1}), as D = c1 R and E = c1 c2 R; more generally
// R^(n) = (-1)^n R c_1 ... c_n.  On [1.5, 32) those derivatives are tabulated
// at a grid of nodes and R, D and E are Taylor series from the nearest node.
struct Mills { double R, D, E; };

// R(z) = sqrt(2 pi) n(z)/d(z) for 0 <= z < FAST_NORM_SPLIT, the rational of
// fast_norm_phi with its exponential cancelled.
inline double mills_rational(double z) {
  const double n = (((((FAST_NORM_N6 * z + FAST_NORM_N5) * z + FAST_NORM_N4) * z + FAST_NORM_N3) * z + FAST_NORM_N2) * z + FAST_NORM_N1) * z + FAST_NORM_N0;
  const double d = ((((((FAST_NORM_M7 * z + FAST_NORM_M6) * z + FAST_NORM_M5) * z + FAST_NORM_M4) * z + FAST_NORM_M3) * z + FAST_NORM_M2) * z + FAST_NORM_M1) * z + FAST_NORM_M0;
  return FAST_NORM_RT2PI * n / d;
}

// R(x) for x < 4; below zero R(x) = 1/phi(x) - R(-x), which loses at most a
// bit because R(-x) < R(0) < 1/phi(x) / 2.
inline double mills_small(double x) {
  if (x >= 0.0) return mills_rational(x);
  const double z = -x;
  const double Rz = z < FAST_NORM_SPLIT ? mills_rational(z) : mills_ratio_std(z);
  return FAST_NORM_RT2PI * std::exp(0.5 * z * z) - Rz;
}

// Taylor table: nodes x_j on [1.5, 8) at spacing 1/8 and on [8, 32) at
// spacing 1/2, and for each node the coefficients (-1)^n P_{n+i}/n!,
// i = 0, 1, 2 (R, D, E), with P_n = R c_1 ... c_n from a fraction deep enough
// to be exact at the lowest node.  Within half a spacing of a node the series
// terms fall by a factor of at least ten per order (c_n is at most about
// n/x), so the alternating sums lose nothing.
constexpr double kTabLo = 1.5, kTabMid = 8.0, kTabHi = 32.0;
constexpr double kTabH0 = 0.125, kTabH1 = 0.5;
constexpr int kTabN0 = 52, kTabN1 = 48;  // nodes below / above kTabMid
constexpr int kTabDeg = 15;

struct MillsTable { double c[kTabN0 + kTabN1 + 1][3][kTabDeg + 1]; };

inline double mills_tab_node(int j) {
  return j <= kTabN0 ? kTabLo + j * kTabH0 : kTabMid + (j - kTabN0) * kTabH1;
}

inline const MillsTable& mills_table() {
  static const MillsTable* tab = [] {
    auto* t = new MillsTable;
    constexpr int K = 4000;
    for (int j = 0; j <= kTabN0 + kTabN1; ++j) {
      const double x = mills_tab_node(j);
      double ck[kTabDeg + 3 + 1];
      double c = 0.0;
      for (int k = K; k >= 1; --k) {
        c = k / (x + c);
        if (k <= kTabDeg + 3) ck[k] = c;
      }
      double P[kTabDeg + 3];
      P[0] = 1.0 / (x + ck[1]);
      for (int n = 1; n < kTabDeg + 3; ++n) P[n] = P[n - 1] * ck[n];
      double fact = 1.0;
      for (int n = 0; n <= kTabDeg; ++n) {
        if (n > 0) fact *= n;
        const double sg = (n & 1) ? -1.0 : 1.0;
        for (int i = 0; i < 3; ++i) t->c[j][i][n] = sg * P[n + i] / fact;
      }
    }
    return t;
  }();
  return *tab;
}

// For kTabLo <= x < kTabHi - kTabH1/2.
inline Mills mills_tab(double x) {
  const int j = x < kTabMid - 0.5 * kTabH0
    ? static_cast<int>((x - kTabLo) / kTabH0 + 0.5)
    : kTabN0 + static_cast<int>((x - kTabMid) / kTabH1 + 0.5);
  const double dx = x - mills_tab_node(j);
  const double (*c)[kTabDeg + 1] = mills_table().c[j];
  double R = c[0][kTabDeg], D = c[1][kTabDeg], E = c[2][kTabDeg];
  for (int n = kTabDeg - 1; n >= 0; --n) {
    R = R * dx + c[0][n];
    D = D * dx + c[1][n];
    E = E * dx + c[2][n];
  }
  return {R, D, E};
}

// Above the table the fraction is entered at depth K with the tail seeded at
// the fixed point of c = (K + 1)/(x + c); each depth keeps R, D and E within
// an ulp or two over its range of x.
inline Mills mills_seq(double x) {
  if (x < kTabLo) {
    const double R = mills_small(x);
    const double D = 1.0 - x * R;
    return {R, D, R - x * D};
  }
  if (x < kTabHi - 0.5 * kTabH1) return mills_tab(x);
  if (!R_FINITE(x)) return {0.0, 0.0, 0.0};
  const int K = x < 40.0 ? 9 : 7;
  const double tail = 2.0 * (K + 1) / (x + std::sqrt(x * x + 4.0 * (K + 1)));
  double c2;
  if (x < 20.0) {
    // c2 = 2/(x + 3/(x + ... K/(x + tail))) from its convergents by Wallis'
    // forward recurrences, which only multiply and add; with every element
    // positive they are stable, and over this range they stay below ~1e40.
    double A0 = 1.0, A1 = 0.0, B0 = 0.0, B1 = 1.0;  // index j - 2, j - 1
    for (int j = 1; j < K; ++j) {
      const double b = (j == K - 1) ? x + tail : x, a = j + 1.0;
      const double A = b * A1 + a * A0, B = b * B1 + a * B0;
      A0 = A1; A1 = A; B0 = B1; B1 = B;
    }
    c2 = A1 / B1;
  } else {
    double c = tail;
    for (int k = K; k >= 2; --k) c = k / (x + c);
    c2 = c;
  }
  const double c1 = 1.0 / (x + c2);
  const double R = 1.0 / (x + c1);
  return {R, c1 * R, c1 * c2 * R};
}

// int_a^b D (order 1, = R(a) - R(b)) or int_a^b E (order 2, = D(a) - D(b)),
// a < b, from the Mills values already formed at both ends; by 8-point
// Gauss-Legendre where the ends agree too closely to subtract.  The width
// w = b - a is passed in its own closed form, since b - a formed from the
// rounded ends loses digits exactly when they are close.
inline double mills_int(int order, double a, const Mills& ma, double b, const Mills& mb,
                        double w) {
  if (!(w > 0.0)) return 0.0;
  const double xa = order == 1 ? ma.R : ma.D;
  const double xb = order == 1 ? mb.R : mb.D;
  if (xb < kMillsCancel * xa) return xa - xb;
  const GLRule& g = gl_get_rule(kGLNodes);
  const double h = 0.5 * w, c = a + h;
  double sum = 0.0;
  for (int i = 0; i < kGLNodes; ++i) {
    const Mills mi = mills_seq(c + h * g.x[i]);
    sum += g.w[i] * (order == 1 ? mi.D : mi.E);
  }
  return h * sum;
}

inline double log_mills_int(int order, double a, double b, double w) {
  if (!(w > 0.0)) return R_NegInf;
  return std::log(mills_int(order, a, mills_seq(a), b, mills_seq(b), w));
}

// Coordinates of the distance-d first passage at mean evidence m = v u.
// wd = beta - alpha = 2d/z and wm = beta + alpha = 2m/z.
struct Pt {
  double z, alpha, beta, lphi, wd, wm;
};

inline Pt point(double m, double d, double s) {
  Pt p;
  p.z = s * std::sqrt(d);
  p.alpha = (m - d) / p.z;
  p.beta = (m + d) / p.z;
  p.wd = 2.0 * d / p.z;
  p.wm = 2.0 * m / p.z;
  p.lphi = -0.5 * p.alpha * p.alpha - LOG_SQRT_2PI;
  return p;
}

struct Kernel {
  static double log_pdf_point(double u, double v, double d, double s) {
    const double m = v * u;
    const Pt p = point(m, d, s);
    const Mills mb = mills_seq(p.beta);
    return std::log(2.0 * v / p.z) + p.lphi + std::log(mb.D + (m / p.z) * mb.R);
  }

  static double log_surv_point(double u, double v, double d, double s) {
    const Pt p = point(v * u, d, s);
    const double Rb = mills_seq(p.beta).R;
    if (p.alpha > 0.0) return p.lphi + std::log(mills_seq(p.alpha).R + Rb);
    return std::log(pnorm_std(-p.alpha) + std::exp(p.lphi) * Rb);
  }

  // F = phi(alpha) int_{-alpha}^{beta} D.  At alpha >= 0 the integral is only
  // taken where S >= 1/2; the interval then starts at a small negative value,
  // where D = 1 - xR is still formed without loss.
  static double log_cdf_point(double u, double v, double d, double s) {
    const Pt p = point(v * u, d, s);
    if (p.alpha >= 0.0) {
      const double ls = log_surv_point(u, v, d, s);
      if (ls < -M_LN2) return std::log(-std::expm1(ls));
    }
    return p.lphi + log_mills_int(1, -p.alpha, p.beta, p.wm);
  }

  // The J form that log_pdf differences first: Jc when the start range is
  // mostly passed (alpha summed over its ends above zero, since J ~
  // Phi(alpha)), J otherwise.
  static bool prefer_Jc(double m, double B, double A, double s) {
    const double d2 = B + A;
    return (m - B) / (s * std::sqrt(B)) + (m - d2) / (s * std::sqrt(d2)) > 0.0;
  }

  // Each end in the complementary forms asked for (want: 1 J, 2 Jc), NaN for
  // a form not formed.  Both are positive sums,
  //   J  = phi(alpha) [int_{-alpha}^{beta} D + 2E(beta) + (alpha + beta) D(beta)],
  //   Jc = phi(alpha) [int_alpha^beta D + (beta - alpha) D(beta)],
  // and at alpha >= 0, where Jc <= 1/2 because J >= Phi(alpha), J is taken as
  // 1 - Jc.  Beyond alpha = -kAlphaMax, Jc is 1 - J.
  static constexpr double kAlphaMax = 30.0;

  static EndJ end_J(double m, double d, double s, int want = 3) {
    if (d <= 0.0) return {0.0, R_NegInf};  // J(0) = 1
    const Pt p = point(m, d, s);
    const Mills mb = mills_seq(p.beta);
    EndJ e = {R_NaN, R_NaN};
    if (p.alpha >= 0.0) {
      e.lJc = p.lphi + std::log(mills_int(1, p.alpha, mills_seq(p.alpha), p.beta, mb, p.wd) +
                                p.wd * mb.D);
      e.lJ = std::log1p(-std::exp(e.lJc));
      return e;
    }
    if ((want & 1) || p.alpha <= -kAlphaMax) {
      e.lJ = p.lphi + std::log(mills_int(1, -p.alpha, mills_seq(-p.alpha), p.beta, mb, p.wm) +
                               2.0 * mb.E + p.wm * mb.D);
    }
    if (want & 2) {
      e.lJc = p.alpha > -kAlphaMax
        ? p.lphi + std::log(mills_int(1, p.alpha, mills_seq(p.alpha), p.beta, mb, p.wd) +
                            p.wd * mb.D)
        : std::log(-std::expm1(clamp_log_prob(e.lJ)));
    }
    return e;
  }

  // H = z phi(alpha) int_{-alpha}^{beta} E and L = z phi(alpha) int_alpha^beta E
  // differ by H - L = m - d, so the smaller is formed as an integral and the
  // larger as a positive sum.
  static EndHL end_HL(double m, double d, double s) {
    if (d <= 0.0) return {std::log(m), R_NegInf};  // H(0) = m, L(0) = 0
    const Pt p = point(m, d, s);
    const Mills mb = mills_seq(p.beta);
    const double lzphi = std::log(p.z) + p.lphi;
    EndHL e;
    if (p.alpha >= 0.0) {
      e.lL = lzphi + std::log(mills_int(2, p.alpha, mills_seq(p.alpha), p.beta, mb, p.wd));
      e.lH = log_add(std::log(m - d), e.lL);
    } else {
      e.lH = lzphi + std::log(mills_int(2, -p.alpha, mills_seq(-p.alpha), p.beta, mb, p.wm));
      e.lL = log_add(std::log(d - m), e.lH);
    }
    return e;
  }

  static void ends_J(double m, double d1, double d2, double s, int want, EndJ* e1, EndJ* e2) {
    *e1 = end_J(m, d1, s, want);
    *e2 = end_J(m, d2, s, want);
  }
  static void ends_HL(double m, double d1, double d2, double s, EndHL* e1, EndHL* e2) {
    *e1 = end_HL(m, d1, s);
    *e2 = end_HL(m, d2, s);
  }

  // The same ends in linear arithmetic, which saves the logarithms of every
  // term; only for |alpha| < kAlphaMax, where phi(alpha) > 1e-196.
  static bool lin_HL(double m, double d, double s, double* H, double* L) {
    if (d <= 0.0) { *H = m; *L = 0.0; return true; }
    const Pt p = point(m, d, s);
    if (!(std::fabs(p.alpha) < kAlphaMax)) return false;
    const Mills mb = mills_seq(p.beta);
    const double zphi = p.z * std::exp(p.lphi);
    if (p.alpha >= 0.0) {
      *L = zphi * mills_int(2, p.alpha, mills_seq(p.alpha), p.beta, mb, p.wd);
      *H = *L + (m - d);
    } else {
      *H = zphi * mills_int(2, -p.alpha, mills_seq(-p.alpha), p.beta, mb, p.wm);
      *L = *H + (d - m);
    }
    return true;
  }

  static bool lin_J(double m, double d, double s, int want, double* J, double* Jc) {
    if (d <= 0.0) { *J = 1.0; *Jc = 0.0; return true; }
    const Pt p = point(m, d, s);
    if (!(std::fabs(p.alpha) < kAlphaMax)) return false;
    const Mills mb = mills_seq(p.beta);
    const double phi = std::exp(p.lphi);
    *J = *Jc = R_NaN;
    if (p.alpha >= 0.0 || (want & 2)) {
      *Jc = phi * (mills_int(1, p.alpha, mills_seq(p.alpha), p.beta, mb, p.wd) +
                   p.wd * mb.D);
    }
    if (p.alpha >= 0.0) {
      *J = 1.0 - *Jc;
    } else if (want & 1) {
      *J = phi * (mills_int(1, -p.alpha, mills_seq(-p.alpha), p.beta, mb, p.wm) +
                  2.0 * mb.E + p.wm * mb.D);
    }
    return true;
  }

  // log M(d), M = E[min(G, d)] = m - H = d - L, as the positive sum
  //   M = m phi(alpha) int_alpha^beta D + d F_d
  // (the partial mean E[G; G < d] = m [Phi(-alpha) - phi(alpha) R(beta)] plus
  // d P(G >= d)), so it keeps its digits where H is close to m.
  static double log_end_M(double m, double d, double s) {
    if (d <= 0.0) return R_NegInf;
    const Pt p = point(m, d, s);
    const Mills mb = mills_seq(p.beta);
    double lI;  // log [Phi(-alpha) - phi(alpha) R(beta)]
    if (p.alpha > -kAlphaMax) {
      lI = p.lphi + std::log(mills_int(1, p.alpha, mills_seq(p.alpha), p.beta, mb, p.wd));
    } else {
      lI = std::log(pnorm_std(-p.alpha) - std::exp(p.lphi) * mb.R);
    }
    const double lFd = p.alpha >= 0.0
      ? std::log1p(-std::exp(p.lphi) * (mills_seq(p.alpha).R + mb.R))
      : p.lphi + std::log(mills_int(1, -p.alpha, mills_seq(-p.alpha), p.beta, mb, p.wm));
    return log_add(std::log(m) + lI, std::log(d) + lFd);
  }

  static double fast_cdf_surv(double m, double B, double A, double s, bool want_S) {
    double H1, L1, H2, L2;
    if (!lin_HL(m, B, s, &H1, &L1) || !lin_HL(m, B + A, s, &H2, &L2)) return R_NaN;
    return lin_cdf_surv(H1, L1, H2, L2, A, want_S);
  }

  static double fast_pdf(double m, double v, double B, double A, double s) {
    const int want = prefer_Jc(m, B, A, s) ? 2 : 1;
    double J1, Jc1, J2, Jc2;
    if (!lin_J(m, B, s, want, &J1, &Jc1) || !lin_J(m, B + A, s, want, &J2, &Jc2))
      return R_NaN;
    double n = lin_J_diff(J1, Jc1, J2, Jc2, want == 2);
    if (ISNAN(n)) {  // form the other one where it was skipped
      lin_J(m, B, s, 3, &J1, &Jc1);
      lin_J(m, B + A, s, 3, &J2, &Jc2);
      n = lin_J_diff(J1, Jc1, J2, Jc2, want == 2);
    }
    return std::log(v / A) + n;
  }

  // Distance d = B + A U, then (s/v) times the maximum over [0, d] of
  // W + w/s, drawn from its endpoint and the Brownian-bridge maximum.
  static double draw(double v, double B, double A, double s) {
    const double d = B + A * R::unif_rand();
    const double x = d / s + std::sqrt(d) * R::norm_rand();
    const double e = -2.0 * d * std::log(R::unif_rand());
    const double r = std::sqrt(x * x + e);
    // (x + r)/2, as e/(2(r - x)) where x < 0 would cancel against r.
    const double mx = x >= 0.0 ? 0.5 * (x + r) : 0.5 * e / (r - x);
    return s * mx / v;
  }
};

}  // namespace ig

// ===========================================================================
// Compound Poisson with exponential packets
// ===========================================================================
namespace cp {

constexpr double kTermEps = 1e-18;
constexpr double kSideFull = -1.0005003335835335e-03;  // log(1 - 1e-3)
constexpr double kSideFullLin = 1.0 - 1e-3;
constexpr double kRescale = 1e250;
constexpr double kLogRescale = 575.6462732485114;  // log(1e250)

// sqrt(a) - sqrt(b) without the cancellation of the root difference.
inline double root_gap(double a, double b) {
  const double r = std::sqrt(a) + std::sqrt(b);
  return r > 0.0 ? (a - b) / r : 0.0;
}

// Upper tails T(k) = P(X >= k), k = 1..n, of X ~ Pois(lam) into out[0..n-1],
// in positive arithmetic: as 1 - CDF(k - 1) when the bulk lies beyond n (the
// CDF there is at most about one half), otherwise as tail sums of the pmf
// formed relative to Pois(1; lam), whose terms beyond n are summed to
// rounding first.
inline void pois_tails(double lam, int n, double* out) {
  if (!(lam > 0.0)) { for (int k = 0; k < n; ++k) out[k] = 0.0; return; }
  if (lam >= n) {
    double p = std::exp(-lam), cdf = 0.0;
    for (int k = 1; k <= n; ++k) { cdf += p; out[k - 1] = 1.0 - cdf; p *= lam / k; }
    return;
  }
  double w = 1.0;
  for (int k = 1; k <= n; ++k) { out[k - 1] = w; w *= lam / (k + 1); }
  double rest = 0.0;  // sum_{j > n} w_j, terms falling by lam/(j + 1) < 1
  for (int j = n + 1; w > 1e-18 * (rest + out[0]); ++j) { rest += w; w *= lam / (j + 1); }
  const double p1 = lam * std::exp(-lam);
  for (int k = n; k >= 1; --k) { rest += out[k - 1]; out[k - 1] = p1 * rest; }
}

// log(e^-z I_0(z)): the power series sum_k (z^2/4)^k / k!^2 below z = 25,
// the Hankel expansion above, where its smallest term is ~e^{-2z}.
inline double log_i0_scaled(double z) {
  if (z >= 25.0) return opr::log_bessel_i_scaled_hankel(0.0, z);
  const double q = 0.25 * z * z;
  double term = 1.0, sum = 1.0;
  for (int k = 1; k < 200; ++k) {
    term *= q / (static_cast<double>(k) * k);
    sum += term;
    if (term < 1e-17 * sum) break;
  }
  return std::log(sum) - z;
}

// One side of the Skellam law: q_k = P(Delta = sgn k) for k >= 0, with
// q_{k+1}/q_k = sqrt(t) rho_k, rho_k = I_{k+1}(z)/I_k(z), and t = a/b for
// the upper side, b/a for the lower.  Returns sum_{k>=1} q_k and
// sum_{k>=1} k q_k, both relative to q_0, and nrm with
// e^-z I_0(z) = 1/(1 + nrm); NaN where the sums fail.
struct Side { double s1, sk, nrm; };

// The ratios come from the stable backward recurrence
// rho_k = 1/(2(k+1)/z + rho_{k+1}), seeded at Amos' approximation far enough
// above the last index used, kstop, that the seed's error has decayed by
// prod rho_k^2 ~ exp(-(K^2 - kstop^2)/z).  kstop comes from the Gaussian
// tail of the Skellam law and kmode is the index of the side's largest term.
// This form divides at every step and is kept for z too small for the
// unnormalised recurrence below to stay in range.
inline Side side_sums_ratio(double z, double t, double kmode, double kstop) {
  static thread_local std::vector<double> rho;
  const double hst = std::sqrt(t);
  int K = static_cast<int>(std::sqrt(kstop * kstop + 40.0 * z) + 10.0);
  for (int attempt = 0; attempt < 8; ++attempt, K *= 2) {
    if (static_cast<int>(rho.size()) < K + 1) rho.resize(K + 1);
    const double kk = K + 1.0;
    double r = z / (kk + std::sqrt(kk * kk + z * z));
    for (int k = K; k >= 0; --k) {
      r = 1.0 / (2.0 * (k + 1.0) / z + r);
      rho[k] = r;
    }
    double L = 0.0, term = 1.0, s0 = 1.0, s1 = 0.0, sk = 0.0, prev = 1.0;
    int k = 0;
    bool done = false;
    for (; k < K; ++k) {
      term *= hst * rho[k];
      const double kn = k + 1.0;
      s0 += term; s1 += term; sk += kn * term;
      if (kn > kmode && term <= prev && kn * term < kTermEps * sk) { done = true; break; }
      prev = term;
      if (s0 > kRescale) {
        term /= kRescale; s0 /= kRescale; s1 /= kRescale; sk /= kRescale;
        prev /= kRescale; L += kLogRescale;
      }
    }
    if (done || term == 0.0) {
      const double sc = std::exp(L);
      return {sc * s1, sc * sk, std::expm1(-log_i0_scaled(z))};
    }
  }
  return {R_NaN, R_NaN, R_NaN};
}

constexpr double kMillerMinZ = 1e-30;
constexpr double kMillerBig = 0x1p664;     // ~1e200
constexpr double kMillerSmall = 0x1p-664;

// The same recurrence run on unnormalised y_k proportional to I_k(z),
// y_{k-1} = (2k/z) y_k + y_{k+1}, from the same seed and depth (raised so
// that it also covers the Bessel terms themselves, which fall like
// e^{-k^2/2z}).  The side sums are accumulated by Horner's rule in the same
// downward pass, so no step divides, and I_0 is normalised by
// I_0 + 2 sum_{k>=1} I_k = e^z.  The depth is doubled while the top terms
// are not negligible.
inline Side side_sums(double z, double t, double kmode, double kstop) {
  if (!(z >= kMillerMinZ)) return side_sums_ratio(z, t, kmode, kstop);
  const double h = std::sqrt(t), tz = 2.0 / z;
  int K = static_cast<int>(std::sqrt(std::fmax(kstop * kstop, 80.0 * z) + 40.0 * z) + 10.0);
  if (static_cast<double>(K) <= kmode) K = static_cast<int>(kmode) + 10;
  for (int attempt = 0; attempt < 8; ++attempt, K *= 2) {
    const double kk = K + 1.0;
    double yb = z / (kk + std::sqrt(kk * kk + z * z)), ya = 1.0;  // y_{k+2}, y_{k+1}
    double P = 0.0, W = 0.0, N = 0.0, yK = 0.0;
    int nsc = 0;
    for (int k = K; k >= 1; --k) {
      const double yk = (k + 1.0) * tz * ya + yb;
      if (k == K) yK = yk;
      P = P * h + yk;       // sum_{j>=k} y_j h^{j-k}
      W = W * h + k * yk;   // sum_{j>=k} j y_j h^{j-k}
      N += yk;
      yb = ya; ya = yk;
      if (W > kMillerBig) {  // W >= y_k, P and N / K
        ya *= kMillerSmall; yb *= kMillerSmall; P *= kMillerSmall;
        W *= kMillerSmall; N *= kMillerSmall; ++nsc;
      }
    }
    const double y0 = tz * ya + yb;
    const double s1 = h * P / y0, sk = h * W / y0;
    // y_K in the final scale: the top terms of both sums must be negligible.
    const double yKf = nsc > 0 ? std::ldexp(yK, -664 * nsc) : yK;
    const double top = K * yKf, lim = kTermEps * h * W;
    if (yKf < 1e-17 * (y0 + 2.0 * N) &&
        (top < lim || top * std::pow(h, K) < lim)) {  // h <= 1 on the far side
      return {s1, sk, 2.0 * N / y0};
    }
  }
  return {R_NaN, R_NaN, R_NaN};
}

// The same side by direct Poisson sums where both pmfs stay in range: with
// outer mean p, inner mean q and C_q, g_q(k) = E(k - N_q)+ = sum_{j<k} C_q(j)
// the inner CDF and stop-loss,
//   t0 = sum_k Pois(k; p) C_q(k),  t1 = sum_k Pois(k; p) C_q(k - 1),
//   tE = sum_k Pois(k; p) g_q(k),
// which are P(Delta <= 0), P(Delta <= -1) and E[Delta-] for Delta = N_q - N_p.
// Every recurrence is forward and positive.  Pois(k; p), C_q and g_q are
// log-concave in k, so their step ratios only fall, and the tail beyond k is
// at most the current summand times r/(1 - r), with r the product of the
// next step ratios; that is tested every fourth step.
// Used while the outer (smaller) mean is at most kDirectMaxOuter, where it
// takes fewer steps than the Bessel recurrence, and e^-q stays far above the
// underflow threshold; a sum that lands near it hands over to the recurrence.
// Only the sums named in the mask W (kNT0, kNT1, kNTE) are formed and tested,
// the others are NaN.
constexpr double kDirectMaxOuter = 60.0;
constexpr double kDirectMaxInner = 650.0;
constexpr double kDirectTiny = 1e-290;
constexpr int kDirectMaxK = 512;
struct Near { double t0, t1, tE; };
constexpr int kNT0 = 1, kNT1 = 2, kNTE = 4, kNAll = 7;

inline bool direct_ok(double a, double b) {
  return a < b ? a <= kDirectMaxOuter && b <= kDirectMaxInner
               : b <= kDirectMaxOuter && a <= kDirectMaxInner;
}

template <int W>
inline bool near_ok(const Near& n) {  // false for NaN
  return (!(W & kNT0) || n.t0 > kDirectTiny) && (!(W & kNT1) || n.t1 > kDirectTiny) &&
         (!(W & kNTE) || n.tE > kDirectTiny);
}

inline const double* inv_int_table() {
  static const std::vector<double> inv = [] {
    std::vector<double> v(kDirectMaxK + 1, 0.0);
    for (int k = 1; k <= kDirectMaxK; ++k) v[k] = 1.0 / k;
    return v;
  }();
  return inv.data();
}

// Whether the wanted sums have converged at index k (see above).  The step
// ratio of the t0 and t1 summands is at most r = p/(k+1) C_q(k+1)/C_q(k), and
// that of the tE summands p/(k+1) (1 + C/g); each tail test
// term r/(1 - r) < eps sum is taken multiplied through by 1 - r.
template <int W>
inline bool direct_done(int k, double p, double q, const double* inv, double pp,
                        double pq, double C, double g, double t0, double t1, double tE) {
  const double pk = p * inv[k + 1];
  if (W & (kNT0 | kNT1)) {
    const double r = pk * (1.0 + q * inv[k + 1] * pq / C);
    const double t = (W & kNT1) ? t1 : t0;
    if (!(r < 1.0) || !(pp * C * r < 0.5 * kTermEps * t * (1.0 - r))) return false;
  }
  if (W & kNTE) {
    const double r = pk * (1.0 + C / g);
    if (!(r < 1.0) || !(pp * g * r < 0.5 * kTermEps * tE * (1.0 - r))) return false;
  }
  return true;
}

inline Near near_out(int W, double t0, double t1, double tE) {
  return {(W & kNT0) ? t0 : R_NaN, (W & kNT1) ? t1 : R_NaN, (W & kNTE) ? tE : R_NaN};
}

template <int W>
inline Near near_direct(double p, double q) {
  const double* inv = inv_int_table();
  double pp = std::exp(-p), pq = std::exp(-q);
  double C = pq, g = 0.0;
  double t0 = pp * C, t1 = 0.0, tE = 0.0;
  for (int k = 1; k < kDirectMaxK; ++k) {
    pp *= p * inv[k];
    pq *= q * inv[k];
    if (W & kNTE) g += C;
    if (W & kNT1) t1 += pp * C;
    C += pq;
    if (W & kNT0) t0 += pp * C;
    if (W & kNTE) tE += pp * g;
    if ((k & 3) == 0 && direct_done<W>(k, p, q, inv, pp, pq, C, g, t0, t1, tE)) {
      return near_out(W, t0, t1, tE);
    }
  }
  return {R_NaN, R_NaN, R_NaN};
}

// Two sides summed in one pass, so that their independent recurrences
// overlap; the pass runs until both have converged.
template <int W>
inline void near_direct2(double p1, double q1, double p2, double q2, Near* o1, Near* o2) {
  const double* inv = inv_int_table();
  double pp1 = std::exp(-p1), pq1 = std::exp(-q1), pp2 = std::exp(-p2), pq2 = std::exp(-q2);
  double C1 = pq1, g1 = 0.0, C2 = pq2, g2 = 0.0;
  double t01 = pp1 * C1, t11 = 0.0, tE1 = 0.0, t02 = pp2 * C2, t12 = 0.0, tE2 = 0.0;
  for (int k = 1; k < kDirectMaxK; ++k) {
    pp1 *= p1 * inv[k]; pq1 *= q1 * inv[k];
    pp2 *= p2 * inv[k]; pq2 *= q2 * inv[k];
    if (W & kNTE) { g1 += C1; g2 += C2; }
    if (W & kNT1) { t11 += pp1 * C1; t12 += pp2 * C2; }
    C1 += pq1; C2 += pq2;
    if (W & kNT0) { t01 += pp1 * C1; t02 += pp2 * C2; }
    if (W & kNTE) { tE1 += pp1 * g1; tE2 += pp2 * g2; }
    if ((k & 3) == 0 && direct_done<W>(k, p1, q1, inv, pp1, pq1, C1, g1, t01, t11, tE1) &&
        direct_done<W>(k, p2, q2, inv, pp2, pq2, C2, g2, t02, t12, tE2)) {
      *o1 = near_out(W, t01, t11, tE1);
      *o2 = near_out(W, t02, t12, tE2);
      return;
    }
  }
  *o1 = *o2 = {R_NaN, R_NaN, R_NaN};
}

// The sums skel needs for need: t0 always (it decides whether the far side is
// summed too), t1 for the probabilities and tE for the stop-losses.
constexpr int kSkelJ = 1, kSkelE = 2, kSkelSF = 4;

inline Near near_direct_need(double p, double q, int need) {
  if (!(need & kSkelE)) return near_direct<kNT0 | kNT1>(p, q);
  if (!(need & (kSkelJ | kSkelSF))) return near_direct<kNT0 | kNTE>(p, q);
  return near_direct<kNAll>(p, q);
}

inline void near_direct2_need(double p1, double q1, double p2, double q2, int need,
                              Near* o1, Near* o2) {
  if (!(need & kSkelE)) near_direct2<kNT0 | kNT1>(p1, q1, p2, q2, o1, o2);
  else if (!(need & (kSkelJ | kSkelSF))) near_direct2<kNT0 | kNTE>(p1, q1, p2, q2, o1, o2);
  else near_direct2<kNAll>(p1, q1, p2, q2, o1, o2);
}

inline bool near_ok_need(const Near& n, int need) {
  if (!(need & kSkelE)) return near_ok<kNT0 | kNT1>(n);
  if (!(need & (kSkelJ | kSkelSF))) return near_ok<kNT0 | kNTE>(n);
  return near_ok<kNAll>(n);
}

// The near side for large z by an expansion in 1/z.  With big >= small the
// near quantities P(N_big <= K_small), P(N_big <= K_small - 1) and
// E(K_small - N_big)+ are tails beyond big of the mixed gamma laws
// Gamma(1 + K) and Gamma(K), whose densities are e^{-x-small} I_0(w) and
// e^{-x-small} sqrt(small/x) I_1(w), w = 2 sqrt(small x).  With
// x = (sqrt(big) + t)^2, delta = sqrt(big) - sqrt(small), tau = t/sqrt(big)
// and I~_nu(w) = e^{-w} I_nu(w),
//
//   P0 = 2 sqrt(big)       int_0^inf e^{-(delta+t)^2} (1 + tau) I~_0(z(1 + tau)) dt,
//   P1 = 2 sqrt(small)     int_0^inf e^{-(delta+t)^2} I~_1(z(1 + tau)) dt,
//   E  = 2 sqrt(small) big int_0^inf e^{-(delta+t)^2} (2 tau + tau^2) I~_1(z(1 + tau)) dt.
//
// Hankel's expansion I~_nu(w) = (2 pi w)^{-1/2} sum_k C_k(nu) w^{-k} and the
// binomial series of (1 + tau)^{1/2-k} or (1 + tau)^{-1/2-k} reduce each to
// the moments int_0^inf e^{-(delta+t)^2} tau^j dt = e^{-delta^2} R(x) nu_j/sqrt(2),
// x = sqrt(2) delta, with nu_j = (U_j/U_0) (2 big)^{-j/2} and
// U_j(x) = int_0^inf s^j e^{-xs-s^2/2} ds: U_0, U_1, U_2 are the Mills R, D, E
// and U_{j+1} = j U_{j-1} - x U_j.  That recurrence is unstable upward for
// large x, but the error it carries into nu_j shrinks like (delta/sqrt(big))^j
// relative to nu_0, so the sums keep their digits.  The Hankel remainder is
// of order e^{-2z} and the binomial series' of order e^{-big}.  nu_j is
// truncated by its bounds Gamma((j+1)/2)/(sqrt(pi) big^{j/2}) (its value at
// x = 0) and j! (1 + 1/x^2) (2 big)^{-j/2}/x^j.  The products
// C_k(nu) binom(p_k, j) are tabulated, so each k is a dot product over j.
// Used where the sums above would need the Bessel recurrence; there
// big >= kDirectMaxOuter.
constexpr double kAsymMinZ = 40.0;
constexpr int kAsymK = 24, kAsymJ = 48;
constexpr double kAsymTol = 1e-17;

struct AsymTable { double c0[kAsymK][kAsymJ + 1], c1[kAsymK][kAsymJ + 1]; };

inline const AsymTable& asym_table() {
  static const AsymTable* tab = [] {
    auto* t = new AsymTable;
    double C0 = 1.0, C1 = 1.0;  // C_k(0), C_k(1)
    for (int k = 0; k < kAsymK; ++k) {
      const double p0 = 0.5 - k, p1 = -0.5 - k;
      double b0 = 1.0, b1 = 1.0;
      for (int j = 0; j <= kAsymJ; ++j) {
        if (j > 0) { b0 *= (p0 - j + 1) / j; b1 *= (p1 - j + 1) / j; }
        t->c0[k][j] = C0 * b0;
        t->c1[k][j] = C1 * b1;
      }
      const double q = (2.0 * k + 1) * (2.0 * k + 1);
      C0 *= q / (8.0 * (k + 1));
      C1 *= (q - 4.0) / (8.0 * (k + 1));
    }
    return t;
  }();
  return *tab;
}

struct NearLog { double l0, l1, lE; };

// big >= small, z = 2 sqrt(big small) >= kAsymMinZ, big >= kDirectMaxOuter.
// l1 and lE are formed only when asked for (NaN otherwise).
inline NearLog near_asym(double big, double small, bool want1, bool wantE) {
  const double z = 2.0 * std::sqrt(big * small);
  const double sb = std::sqrt(big), ss = std::sqrt(small);
  const double delta = (big - small) / (sb + ss);
  const double x = M_SQRT2 * delta, r = 1.0 / std::sqrt(2.0 * big), ib = 1.0 / big;
  const ig::Mills ml = ig::mills_seq(x);
  double nu[kAsymJ + 3];
  nu[0] = 1.0; nu[1] = r * ml.D / ml.R; nu[2] = r * r * ml.E / ml.R;
  double G0 = 1.0, G1 = 1.0 / std::sqrt(M_PI * big), L = 1.0 + 1.0 / (x * x);
  int J = 0;
  for (;; ++J) {
    if (J >= 3) nu[J] = r * ((J - 1) * r * nu[J - 2] - x * nu[J - 1]);
    if (G0 < kAsymTol || L < kAsymTol || J >= kAsymJ) break;
    const double G2 = G0 * 0.5 * (J + 1) * ib;
    G0 = G1; G1 = G2;
    L *= (J + 1) * r / x;
  }
  if (J < 2) J = 2;
  double eps[kAsymJ + 1];
  if (wantE) {
    for (int j = J + 1; j <= J + 2; ++j) nu[j] = r * ((j - 1) * r * nu[j - 2] - x * nu[j - 1]);
    for (int j = 0; j <= J; ++j) eps[j] = 2.0 * nu[j + 1] + nu[j + 2];
  }
  const AsymTable& tab = asym_table();
  const double iz = 1.0 / z;
  double s0 = 0.0, s1 = 0.0, sE = 0.0, zk = 1.0;
  for (int k = 0; k < kAsymK; ++k, zk *= iz) {
    const double* c0 = tab.c0[k];
    const double* c1 = tab.c1[k];
    double a0 = 0.0, a1 = 0.0, aE = 0.0;
    for (int j = 0; j <= J; ++j) a0 += c0[j] * nu[j];
    if (want1) for (int j = 0; j <= J; ++j) a1 += c1[j] * nu[j];
    if (wantE) for (int j = 0; j <= J; ++j) aE += c1[j] * eps[j];
    s0 += zk * a0; s1 += zk * a1; sE += zk * aE;
    if (zk * std::fabs(c0[0]) < kAsymTol && zk * std::fabs(c1[0]) < kAsymTol) break;
  }
  const double lC = -0.5 * std::log(4.0 * M_PI * z) - delta * delta + std::log(ml.R);
  return {lC + std::log(2.0 * sb * s0),
          want1 ? lC + std::log(2.0 * ss * s1) : R_NaN,
          wantE ? lC + std::log(2.0 * ss * big * sE) : R_NaN};
}

// The Skellam quantities of one (a, b), each on its small side; only the
// pairs asked for in need are formed, NaN otherwise.  lo_direct, when given,
// is the near side already summed directly for the same need.
struct Skel { double lS, lF, lJ, lJc, lEpos, lEneg; };

inline Skel skel(double a, double b, int need, const Near* lo_direct = nullptr) {
  const double z = 2.0 * std::sqrt(a * b);
  const double sd = root_gap(a, b);
  const double gap = std::fabs(a - b);
  // Terms beyond k with (k -/+ gap)^2 - gap^2 > 82 (a + b) are below e^-41
  // of the largest.
  const double tail2 = 82.0 * (a + b);
  // The side away from the mean has decreasing terms; the other is only
  // summed when the first carries nearly all the mass, so that its complement
  // would lose more than three digits.
  const bool up_mean = a >= b;
  double near0, near1 = R_NaN, far0 = R_NaN, far1 = R_NaN, nearE = R_NaN, farE = R_NaN;
  const bool want1 = ((need & kSkelJ) && up_mean) || ((need & kSkelSF) && !up_mean);
  const bool want0 = ((need & kSkelJ) && !up_mean) || ((need & kSkelSF) && up_mean);
  auto lg0 = [](double x) { return x > 0.0 ? std::log(x) : R_NegInf; };
  Near lo_d = {R_NaN, R_NaN, R_NaN};
  NearLog lo_a;
  if (direct_ok(a, b)) {
    lo_d = lo_direct ? *lo_direct : near_direct_need(up_mean ? b : a, up_mean ? a : b, need);
  }
  if (near_ok_need(lo_d, need)) {
    const Near& lo = lo_d;
    near0 = std::log(lo.t0);
    const Near hi = near0 > kSideFull ? near_direct<kNAll>(up_mean ? a : b, up_mean ? b : a)
                                      : Near{R_NaN, R_NaN, R_NaN};
    const bool use_hi = !ISNAN(hi.t0);
    if (want1) {
      near1 = lg0(lo.t1);
      far0 = use_hi ? std::log(hi.t0) : log1m_exp(clamp_log_prob(near1));  // excludes 0
    }
    if (want0) far1 = use_hi ? lg0(hi.t1) : log1m_exp(clamp_log_prob(near0));
    if (need & kSkelE) {
      nearE = lg0(lo.tE);
      farE = use_hi ? lg0(hi.tE) : log_add(std::log(gap), nearE);
    }
  } else if (z >= kAsymMinZ && std::fmax(a, b) >= kDirectMaxOuter &&
             (lo_a = near_asym(up_mean ? a : b, up_mean ? b : a, want1, need & kSkelE),
              lo_a.l0 <= kSideFull)) {
    // The far side as complements, as when the near side is summed directly.
    near0 = lo_a.l0;
    if (want1) {
      near1 = lo_a.l1;
      far0 = log1m_exp(clamp_log_prob(near1));
    }
    if (want0) far1 = log1m_exp(clamp_log_prob(near0));
    if (need & kSkelE) {
      nearE = lo_a.lE;
      farE = log_add(std::log(gap), nearE);
    }
  } else {
    const Side lo = side_sums(z, up_mean ? b / a : a / b, 0.0, std::sqrt(gap * gap + tail2) - gap);
    const double lp0 = -sd * sd - std::log1p(lo.nrm);
    // lo is P(Delta <= 0), P(Delta <= -1), E[Delta-] when a >= b, and
    // P(Delta >= 0), P(Delta >= 1), E[Delta+] otherwise; far is the other side.
    near0 = lp0 + std::log1p(lo.s1);
    const bool use_hi = near0 > kSideFull;
    const Side hi = use_hi ? side_sums(z, up_mean ? a / b : b / a, gap, gap + std::sqrt(tail2))
                           : Side{R_NaN, R_NaN, R_NaN};
    auto lg = [&](double x) { return x > 0.0 ? lp0 + std::log(x) : R_NegInf; };
    if (want1) {
      near1 = lg(lo.s1);
      far0 = use_hi ? lp0 + std::log1p(hi.s1) : log1m_exp(clamp_log_prob(near1));  // excludes 0
    }
    if (want0) far1 = use_hi ? lg(hi.s1) : log1m_exp(clamp_log_prob(near0));
    if (need & kSkelE) {
      nearE = lg(lo.sk);
      farE = use_hi ? lg(hi.sk) : log_add(std::log(gap), nearE);
    }
  }
  Skel o;
  if (up_mean) {
    o.lS = near0; o.lJc = near1; o.lEneg = nearE;
    o.lJ = far0; o.lF = far1; o.lEpos = farE;
  } else {
    o.lJ = near0; o.lF = near1; o.lEpos = nearE;
    o.lS = far0; o.lJc = far1; o.lEneg = farE;
  }
  return o;
}

// The Skellam quantities at two distances sharing one a, with the direct
// near sides summed together.
inline void skel2(double a, double b1, double b2, int need, Skel* o1, Skel* o2) {
  if (direct_ok(a, b1) && direct_ok(a, b2)) {
    Near n1, n2;
    near_direct2_need(a >= b1 ? b1 : a, a >= b1 ? a : b1, a >= b2 ? b2 : a, a >= b2 ? a : b2,
                      need, &n1, &n2);
    *o1 = skel(a, b1, need, &n1);
    *o2 = skel(a, b2, need, &n2);
  } else {
    *o1 = skel(a, b1, need);
    *o2 = skel(a, b2, need);
  }
}

struct Kernel {
  static double log_pdf_point(double u, double v, double d, double s) {
    const double c = s * s, a = 2.0 * v * u / c, b = 2.0 * d / c;
    const double sd = root_gap(a, b);
    return std::log(2.0 * v / c) - sd * sd + log_i0_scaled(2.0 * std::sqrt(a * b));
  }
  static double log_surv_point(double u, double v, double d, double s) {
    const double c = s * s;
    return skel(2.0 * v * u / c, 2.0 * d / c, kSkelSF).lS;
  }
  static double log_cdf_point(double u, double v, double d, double s) {
    const double c = s * s;
    return skel(2.0 * v * u / c, 2.0 * d / c, kSkelSF).lF;
  }
  static bool prefer_Jc(double, double, double, double) { return false; }

  // log M(d), M = E[min(G, d)] = m - H = (c/2) E[min(N, K)]
  //   = (c/2) sum_{k>=1} P(N >= k) P(K >= k),
  // a positive series cut where the smaller mean's tail is negligible; NaN
  // where both means are large (the quadrature then settles the ends).
  static double log_end_M(double m, double d, double s) {
    if (d <= 0.0) return R_NegInf;
    const double c = s * s, a = 2.0 * m / c, b = 2.0 * d / c;
    const double lo = std::fmin(a, b);
    if (!(lo <= 60.0)) return R_NaN;
    const int n = static_cast<int>(std::ceil(lo + 12.0 * std::sqrt(lo) + 30.0));
    double ta[200], tb[200];
    pois_tails(a, n, ta);
    pois_tails(b, n, tb);
    double sum = 0.0;
    for (int k = n - 1; k >= 0; --k) sum += ta[k] * tb[k];
    return sum > 0.0 ? std::log(0.5 * c * sum) : R_NegInf;
  }
  // Linear ends from the direct sums, when both ends have them and the start
  // range does not begin at zero.  With up = (a >= b) the near side is
  // P(Delta <= 0), P(Delta <= -1), E[Delta-]; otherwise P(Delta >= 0),
  // P(Delta >= 1), E[Delta+].  The far probabilities are complements, kept
  // only while the near side leaves them three digits (as in skel); the far
  // stop-loss is the near one plus |a - b|.
  template <int W>
  static bool lin_near(double m, double B, double A, double s, double* a,
                       double* b1, double* b2, Near* n1, Near* n2) {
    if (!(B > 0.0)) return false;
    const double c = s * s;
    *a = 2.0 * m / c; *b1 = 2.0 * B / c; *b2 = 2.0 * (B + A) / c;
    if (!direct_ok(*a, *b1) || !direct_ok(*a, *b2)) return false;
    near_direct2<W>(*a >= *b1 ? *b1 : *a, *a >= *b1 ? *a : *b1,
                    *a >= *b2 ? *b2 : *a, *a >= *b2 ? *a : *b2, n1, n2);
    return near_ok<W>(*n1) && near_ok<W>(*n2);
  }

  static double fast_cdf_surv(double m, double B, double A, double s, bool want_S) {
    double a, b1, b2;
    Near n1, n2;
    if (!lin_near<kNTE>(m, B, A, s, &a, &b1, &b2, &n1, &n2)) return R_NaN;
    const double hc = 0.5 * s * s;
    auto HL = [&](double b, const Near& n, double* H, double* L) {
      if (a >= b) { *L = hc * n.tE; *H = *L + hc * (a - b); }
      else        { *H = hc * n.tE; *L = *H + hc * (b - a); }
    };
    double H1, L1, H2, L2;
    HL(b1, n1, &H1, &L1);
    HL(b2, n2, &H2, &L2);
    return lin_cdf_surv(H1, L1, H2, L2, A, want_S);
  }

  static double fast_pdf(double m, double v, double B, double A, double s) {
    double a, b1, b2;
    Near n1, n2;
    if (!lin_near<kNT0 | kNT1>(m, B, A, s, &a, &b1, &b2, &n1, &n2)) return R_NaN;
    if (n1.t0 > kSideFullLin || n2.t0 > kSideFullLin) return R_NaN;
    auto J = [&](double b, const Near& n, double* Jv, double* Jc) {
      if (a >= b) { *Jc = n.t1; *Jv = 1.0 - n.t1; }
      else        { *Jv = n.t0; *Jc = 1.0 - n.t0; }
    };
    double J1, Jc1, J2, Jc2;
    J(b1, n1, &J1, &Jc1);
    J(b2, n2, &J2, &Jc2);
    return std::log(v / A) + lin_J_diff(J1, Jc1, J2, Jc2, false);
  }
  // Both forms come from one side sum; want is ignored.
  static EndJ end_J(double m, double d, double s, int /*want*/ = 3) {
    if (d <= 0.0) return {0.0, R_NegInf};
    const double c = s * s;
    const Skel k = skel(2.0 * m / c, 2.0 * d / c, kSkelJ);
    return {k.lJ, k.lJc};
  }
  static EndHL end_HL(double m, double d, double s) {
    if (d <= 0.0) return {std::log(m), R_NegInf};
    const double c = s * s, lth = std::log(0.5 * c);
    const Skel k = skel(2.0 * m / c, 2.0 * d / c, kSkelE);
    return {lth + k.lEpos, lth + k.lEneg};
  }
  // Both ends of the start range, summed together.
  static void ends_J(double m, double d1, double d2, double s, int want, EndJ* e1, EndJ* e2) {
    if (d1 <= 0.0) { *e1 = end_J(m, d1, s, want); *e2 = end_J(m, d2, s, want); return; }
    const double c = s * s;
    Skel k1, k2;
    skel2(2.0 * m / c, 2.0 * d1 / c, 2.0 * d2 / c, kSkelJ, &k1, &k2);
    *e1 = {k1.lJ, k1.lJc};
    *e2 = {k2.lJ, k2.lJc};
  }
  static void ends_HL(double m, double d1, double d2, double s, EndHL* e1, EndHL* e2) {
    if (d1 <= 0.0) { *e1 = end_HL(m, d1, s); *e2 = end_HL(m, d2, s); return; }
    const double c = s * s, lth = std::log(0.5 * c);
    Skel k1, k2;
    skel2(2.0 * m / c, 2.0 * d1 / c, 2.0 * d2 / c, kSkelE, &k1, &k2);
    *e1 = {lth + k1.lEpos, lth + k1.lEneg};
    *e2 = {lth + k2.lEpos, lth + k2.lEneg};
  }
  // The packets needed to cross d are 1 + Pois(2d/c); the finish time is the
  // arrival of that many events at rate 2v/c.
  static double draw(double v, double B, double A, double s) {
    const double c = s * s;
    const double d = B + A * R::unif_rand();
    const double n = 1.0 + R::rpois(2.0 * d / c);
    return R::rgamma(n, c / (2.0 * v));
  }
};

}  // namespace cp

// ===========================================================================
// Gamma rate mixtures (Math/gamma_mixtures.tex, secs. 3 and 4)
//
// The trial rate V ~ Gamma(beta, rate beta/v), beta = 1/cv^2, is drawn once
// per trial.  The rate is a pure time scale, so in mean-evidence units
// m = V u the mixing law is Gamma(beta, rate lam = beta/mbar), mbar = v u.
// The averaged density is v E[V' ...] with V' size-biased, which is the same
// law with shape beta + 1; the J ends, which give the density, are therefore
// mixed at shape beta + 1 and the H, L ends and point survivors at beta.
// ===========================================================================

// ---------------------------------------------------------------------------
// Packets.  Given V the decision is the arrival of event 1 + K, K ~ Pois(b),
// b = 2d/c, at rate 2V/c; over V the event count by time u is negative
// binomial with P(N >= a) = I_q(a, p), q = abar/(p + abar), abar = 2 mbar/c
// (the rate, and so q, is the same at shapes beta and beta + 1).  With
// h_a = q^a (1-q)^p / B(a, p) and g_a = h_a / a = I_q(a, p) - I_q(a + 1, p),
// and C, T the Poisson(b) CDF and upper tail,
//
//   F = sum_{a>=1} g_a C(a-1),              S = (1-q)^beta + sum_a g_a T(a),
//   H = (c/2) sum_a g_a E(a - K)+,          L = (c/2) [(1-q)^beta b + sum_a g_a E(K - a)+],
//   J = e^-b + sum_a g'_a P(1 <= K <= a),   Jc = (1-q)^{beta+1}(1 - e^-b) + sum_a g'_a T(a+1),
//   f = (v/mbar) sum_{a>=1} a P(K = a-1) g_a,
//
// (g' at shape beta + 1).  Every weight is log-concave in a, rising
// (C, E(a-K)+, P(1 <= K <= a), a P(K = a-1) up to its mode) or falling
// (T, E(K-a)+), and is built by additions in its rising direction: the
// rising sums run up from a = 1, the falling ones down from a cut a_hi above
// the peak.  For p >= 1 the summands are log-concave too, so each sweep stops
// when its current ratio bounds a geometric remainder below rounding; for
// p < 1 g is log-convex and a downward sweep runs to a = 1.
// ---------------------------------------------------------------------------
namespace cpg {

constexpr double kEps = 1e-17;
constexpr int kMaxTerms = 1000000;

struct QG { double q, log_q, log_omq; };

// q and 1 - q for mean count abar at shape p.
inline QG qg(double abar, double p) {
  const double t = p + abar;
  return {abar / t, std::log(abar) - std::log(t), std::log(p) - std::log(t)};
}

// Rising weights: init() returns log W_1 and next(a) the ratio W_{a+1}/W_a.
// Internal state is kept relative to the first weight and renormalised.
struct WCprev {          // W_a = C(a - 1)
  double b, pi = 1.0, c = 1.0;
  explicit WCprev(double b_) : b(b_) {}
  double init() const { return -b; }
  double next(double a) {
    pi *= b / a;
    const double cn = c + pi, r = cn / c;
    c = cn;
    if (c > cp::kRescale) { c /= cp::kRescale; pi /= cp::kRescale; }
    return r;
  }
};
struct WStop {           // W_a = E(a - K)+ = sum_{n<a} C(n)
  double b, pi = 1.0, c = 1.0, l = 1.0;
  explicit WStop(double b_) : b(b_) {}
  double init() const { return -b; }
  double next(double a) {
    pi *= b / a;
    c += pi;
    const double ln = l + c, r = ln / l;
    l = ln;
    if (l > cp::kRescale) { l /= cp::kRescale; c /= cp::kRescale; pi /= cp::kRescale; }
    return r;
  }
};
struct WCstar {          // W_a = P(1 <= K <= a)
  double b, pi = 1.0, c = 1.0;
  explicit WCstar(double b_) : b(b_) {}
  double init() const { return std::log(b) - b; }
  double next(double a) {
    pi *= b / (a + 1.0);
    const double cn = c + pi, r = cn / c;
    c = cn;
    if (c > cp::kRescale) { c /= cp::kRescale; pi /= cp::kRescale; }
    return r;
  }
};
struct WDens {           // W_a = a P(K = a - 1)
  double b;
  explicit WDens(double b_) : b(b_) {}
  double init() const { return -b; }
  double next(double a) { return (a + 1.0) / a * b / a; }
};

// log sum_{a>=1} g_a W_a at shape p, from g_1 = p q (1-q)^p.  The ratios of
// a log-concave W only fall, and those of g fall to q (p >= 1) or rise to it
// (p < 1), so max(g ratio, q) times the current W ratio bounds every later
// summand ratio.
template <class W>
inline double up_sum(const QG& z, double p, W& w) {
  double L = std::log(p) + z.log_q + p * z.log_omq + w.init();
  if (!(L > R_NegInf)) return R_NegInf;
  double t = 1.0, sum = 1.0;
  for (double a = 1.0;; a += 1.0) {
    if (a > kMaxTerms) return R_NaN;
    const double gr = z.q * (a + p) / (a + 1.0);
    const double r = gr * w.next(a);
    const double rb = (p >= 1.0) ? r : r / gr * std::max(gr, z.q);
    t *= r;
    sum += t;
    if (rb < 1.0 && t * rb < kEps * sum * (1.0 - rb)) break;
    if (sum > cp::kRescale) { t /= cp::kRescale; sum /= cp::kRescale; L += cp::kLogRescale; }
  }
  return L + std::log(sum);
}

// Falling weights from a_hi > b: start() sets the state at a_hi from the
// forward series of the Poisson tail, t_k = b^k a!/(a + k)!, and returns
// log W_{a_hi}; prev(a) returns W_{a-1}/W_a.  Every falling weight has
// W_{a+1}/W_a <= b/(a + 1) above b, which bounds the sum above a_hi.
inline void tail_series(double a, double b, double* s0, double* s1, double* sj) {
  double t = 1.0, k = 0.0;
  *s0 = 1.0; *s1 = 0.0; *sj = 0.0;
  for (int i = 0; i < 100000; ++i) {
    k += 1.0;
    t *= b / (a + k);
    *s1 += t;
    *sj += k * t;
    if (t < 1e-18 * *s1) break;
  }
  *s0 += *s1;
}
struct WT {              // W_a = T(a)
  double b, p = 1.0, T = 1.0;
  explicit WT(double b_) : b(b_) {}
  double start(double a) {
    double s0, s1, sj;
    tail_series(a, b, &s0, &s1, &sj);
    p = 1.0; T = s0;
    return R::dpois(a, b, 1) + std::log(s0);
  }
  double prev(double a) {
    p *= a / b;
    const double Tn = T + p, r = Tn / T;
    T = Tn;
    if (T > cp::kRescale) { T /= cp::kRescale; p /= cp::kRescale; }
    return r;
  }
};
struct WT1 {             // W_a = T(a + 1)
  double b, p = 1.0, T = 1.0;
  explicit WT1(double b_) : b(b_) {}
  double start(double a) {
    double s0, s1, sj;
    tail_series(a, b, &s0, &s1, &sj);
    p = 1.0; T = s1;
    return R::dpois(a, b, 1) + std::log(s1);
  }
  double prev(double a) {
    const double Tn = T + p, r = Tn / T;
    T = Tn;
    p *= a / b;
    if (T > cp::kRescale) { T /= cp::kRescale; p /= cp::kRescale; }
    return r;
  }
};
struct WU {              // W_a = E(K - a)+ = sum_{n>a} T(n)
  double b, p = 1.0, T = 1.0, U = 1.0;
  explicit WU(double b_) : b(b_) {}
  double start(double a) {
    double s0, s1, sj;
    tail_series(a, b, &s0, &s1, &sj);
    p = 1.0; T = s0; U = sj;
    return R::dpois(a, b, 1) + std::log(sj);
  }
  double prev(double a) {
    const double Un = U + T, r = Un / U;
    U = Un;
    p *= a / b;
    T += p;
    if (U > cp::kRescale) { U /= cp::kRescale; T /= cp::kRescale; p /= cp::kRescale; }
    return r;
  }
};

// log sum_{a>=1} g_a W_a for a falling W, swept down from a_hi.  a_hi lies
// where the summand ratio has fallen below one half (g's ratio q(a+p)/(a+1)
// times b/(a+1)) and some way beyond; the remainder above it is checked
// against the finished sum, and the sweep is redone from twice as far up in
// the rare case that it is not negligible.
template <class W>
inline double down_sum(const QG& z, double p, double b, const W& w0) {
  const double qb = z.q * b;
  const double a1 = std::max(qb + std::sqrt(qb * qb + 2.0 * qb * (p + 1.0)), b);
  double a_hi = std::ceil(a1 + 10.0 * std::sqrt(a1) + 30.0);
  for (int attempt = 0; attempt < 8; ++attempt, a_hi *= 2.0) {
    if (a_hi > kMaxTerms) return R_NaN;
    W w = w0;
    const double L0 = a_hi * z.log_q + p * z.log_omq - lbeta_fast(a_hi, p) -
                      std::log(a_hi) + w.start(a_hi);
    if (!(L0 > R_NegInf)) return R_NegInf;
    double L = L0, t = 1.0, sum = 1.0;
    for (double a = a_hi; a > 1.0; a -= 1.0) {
      const double r = a / (z.q * (a - 1.0 + p)) * w.prev(a);
      t *= r;
      sum += t;
      if (p >= 1.0 && r < 1.0 && t * r < kEps * sum * (1.0 - r)) break;
      if (sum > cp::kRescale) { t /= cp::kRescale; sum /= cp::kRescale; L += cp::kLogRescale; }
    }
    // The summands above a_hi fall at least by ru per step.
    const double ru = std::max(z.q * (a_hi + p) / (a_hi + 1.0), p >= 1.0 ? 0.0 : z.q) *
                      b / (a_hi + 1.0);
    const double lsum = L + std::log(sum);
    if (ru < 0.5 && L0 + std::log(ru / (1.0 - ru)) < lsum + std::log(kEps)) return lsum;
  }
  return R_NaN;
}

struct Kernel {
  double beta = 1.0;

  // Point survivor and CDF: F is the rising sum, S the falling one plus
  // (1-q)^beta; whichever the mean count makes small is summed and the other
  // is its complement once it is below one half.
  static double log_F(const QG& z, double beta, double b) {
    WCprev w(b);
    return up_sum(z, beta, w);
  }
  static double log_S(const QG& z, double beta, double b) {
    return log_add(beta * z.log_omq, down_sum(z, beta, b, WT(b)));
  }
  double log_surv_point(double u, double v, double d, double s) const {
    const double c = s * s, abar = 2.0 * v * u / c, b = 2.0 * d / c;
    const QG z = qg(abar, beta);
    if (abar < b + 1.0) {
      const double lF = log_F(z, beta, b);
      if (lF < -M_LN2) return std::log1p(-std::exp(lF));
    }
    return log_S(z, beta, b);
  }
  double log_cdf_point(double u, double v, double d, double s) const {
    const double c = s * s, abar = 2.0 * v * u / c, b = 2.0 * d / c;
    const QG z = qg(abar, beta);
    if (abar >= b + 1.0) {
      const double lS = log_S(z, beta, b);
      if (lS < -M_LN2) return std::log1p(-std::exp(lS));
    }
    return log_F(z, beta, b);
  }
  double log_pdf_point(double u, double v, double d, double s) const {
    const double c = s * s, m = v * u, b = 2.0 * d / c;
    const QG z = qg(2.0 * m / c, beta);
    WDens w(b);
    return std::log(v / m) + up_sum(z, beta, w);
  }
  bool prefer_Jc(double, double, double, double) const { return false; }
  double log_end_M(double, double, double) const { return R_NaN; }
  double fast_cdf_surv(double, double, double, double, bool) const { return R_NaN; }
  double fast_pdf(double, double, double, double, double) const { return R_NaN; }

  // J and Jc at shape beta + 1: J is small while the size-biased mean count
  // abar (beta + 1)/beta is below b.
  EndJ end_J(double m, double d, double s, int /*want*/ = 3) const {
    if (d <= 0.0) return {0.0, R_NegInf};
    const double c = s * s, abar = 2.0 * m / c, b = 2.0 * d / c, p = beta + 1.0;
    const QG z = qg(abar * p / beta, p);
    auto lJ = [&] { WCstar w(b); return log_add(-b, up_sum(z, p, w)); };
    auto lJc = [&] {
      return log_add(p * z.log_omq + std::log(-std::expm1(-b)), down_sum(z, p, b, WT1(b)));
    };
    EndJ e;
    if (abar * p / beta < b) {
      e.lJ = lJ();
      if (!(e.lJ >= -M_LN2)) { e.lJc = std::log1p(-std::exp(e.lJ)); return e; }
      e.lJc = lJc();
      return e;
    }
    e.lJc = lJc();
    e.lJ = (e.lJc < -M_LN2) ? std::log1p(-std::exp(e.lJc)) : lJ();
    return e;
  }

  // H - L = mbar - d, so the smaller is summed and the larger added to it.
  EndHL end_HL(double m, double d, double s) const {
    if (d <= 0.0) return {std::log(m), R_NegInf};
    const double c = s * s, abar = 2.0 * m / c, b = 2.0 * d / c;
    const double lth = std::log(0.5 * c);
    const QG z = qg(abar, beta);
    EndHL e;
    if (m < d) {
      WStop w(b);
      e.lH = lth + up_sum(z, beta, w);
      e.lL = log_add(e.lH, std::log(d - m));
    } else {
      e.lL = lth + log_add(beta * z.log_omq + std::log(b), down_sum(z, beta, b, WU(b)));
      e.lH = (m > d) ? log_add(e.lL, std::log(m - d)) : e.lL;
    }
    return e;
  }
  void ends_J(double m, double d1, double d2, double s, int want, EndJ* e1, EndJ* e2) const {
    *e1 = end_J(m, d1, s, want);
    *e2 = end_J(m, d2, s, want);
  }
  void ends_HL(double m, double d1, double d2, double s, EndHL* e1, EndHL* e2) const {
    *e1 = end_HL(m, d1, s);
    *e2 = end_HL(m, d2, s);
  }
};

// The trial rate, then the fixed-rate packet arrivals at that rate.
inline double draw(double v, double B, double A, double s, double cv) {
  const double V = R::rgamma(1.0 / (cv * cv), v * cv * cv);
  return cp::Kernel::draw(V, B, A, s);
}

}  // namespace cpg

// ---------------------------------------------------------------------------
// Inverse-Gaussian evidence.  There is no closed form for real beta, so each
// mixed quantity is a Gauss rule over m against Gamma(p, rate lam) of the
// stable conditional forms of namespace ig.  Every conditional form is
// phi(alpha) times a slowly varying bracket of Mills functions on one side of
// m = d (alpha = (m - d)/z, z = s sqrt(d)) and a constant or linear limit on
// the other, so each quantity is split at d into Gaussian-tailed integrals
// and closed-form Gamma pieces:
//
//   S = P_p(lam d) - int_0^d g F + int_d^inf g S,    F <= 1/2 below d,
//   F = Q_p(lam d) - int_d^inf g S + int_0^d g F,    S <= 1/2 above d,
//   H = int_0^d g H + E(m - d)+ + int_d^inf g L,     L = H + d - mbar,
//   J = Q_p(lam d) + int_0^d g J - int_d^inf g Jc,   Jc = 1 - J,
//   f = v int_0^inf g' j,                            j = f_d / V,
//
// so that no subtraction loses more than a bit or two.  Below d, J reaches
// 1/2 only when z is comparable to d; there J is no longer small and the
// bracket varies slowly, and Jc is integrated directly instead.
//
// The Gaussian-tailed integrands are m^{p-1} e^{-lam m - (m-d)^2/(2 z^2)}
// times their bracket, so the nodes are placed on that envelope: its maximum
// on the side, where it is completed-square Gaussian times a power, and the
// two points kDrop below it.  Where the lower point reaches the origin the
// power is non-analytic there and goes into a Gauss-Jacobi weight.  This
// keeps relative accuracy where the Gamma law is narrow against z (large
// beta) and in the far tails, where a window fixed by z alone misses the
// product of a Gamma tail and a Gaussian tail.
// ---------------------------------------------------------------------------
namespace igg {

constexpr int kNodes = 32;
constexpr double kDrop = 38.0;

// Gauss-Jacobi rule on [-1, 1] for the weight (1 + u)^c, c > -1, by
// Golub-Welsch: the eigenvalues of the Jacobi matrix of the Jacobi
// polynomials P^(0, c), and the squared first components of its
// eigenvectors, found by implicit QL that carries only those components.
struct JRule { double c = R_NaN; double x[kNodes], w[kNodes]; };

inline bool jacobi_rule(double c, JRule& r) {
  constexpr int n = kNodes;
  double d[n], e[n], z0[n];
  for (int k = 0; k < n; ++k) {
    const double s = 2.0 * k + c;
    d[k] = (k == 0) ? c / (c + 2.0) : c * c / (s * (s + 2.0));
    if (k > 0) e[k - 1] = 2.0 * k * (k + c) / (s * std::sqrt((s + 1.0) * (s - 1.0)));
    z0[k] = (k == 0) ? 1.0 : 0.0;
  }
  e[n - 1] = 0.0;
  for (int l = 0; l < n; ++l) {
    int iter = 0, m;
    do {
      for (m = l; m < n - 1; ++m) {
        const double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
        if (std::fabs(e[m]) <= 1e-16 * dd) break;
      }
      if (m != l) {
        if (iter++ == 60) return false;
        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double rr = std::hypot(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? rr : -rr));
        double s = 1.0, cc = 1.0, p = 0.0;
        int i;
        for (i = m - 1; i >= l; --i) {
          double f = s * e[i];
          const double b = cc * e[i];
          rr = std::hypot(f, g);
          e[i + 1] = rr;
          if (rr == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
          s = f / rr; cc = g / rr;
          g = d[i + 1] - p;
          rr = (d[i] - g) * s + 2.0 * cc * b;
          p = s * rr;
          d[i + 1] = g + p;
          g = cc * rr - b;
          f = z0[i + 1];
          z0[i + 1] = s * z0[i] + cc * f;
          z0[i] = cc * z0[i] - s * f;
        }
        if (rr == 0.0 && i >= l) continue;
        d[l] -= p; e[l] = g; e[m] = 0.0;
      }
    } while (m != l);
  }
  const double mu0 = std::exp((c + 1.0) * M_LN2) / (c + 1.0);
  for (int k = 0; k < n; ++k) { r.x[k] = d[k]; r.w[k] = mu0 * z0[k] * z0[k]; }
  r.c = c;
  return true;
}

// The rules recur for the shapes of the rows of a condition; a small
// round-robin cache keyed on c serves them.
inline const JRule* jrule(double c) {
  constexpr int kSlots = 8;
  static thread_local JRule cache[kSlots];
  static thread_local int next = 0;
  for (int k = 0; k < kSlots; ++k) if (cache[k].c == c) return &cache[k];
  JRule& e = cache[next];
  next = (next + 1) % kSlots;
  if (!jacobi_rule(c, e)) { e.c = R_NaN; return nullptr; }
  return &e;
}

// Gamma(p, rate lam) on m, with lnorm = p log lam - log Gamma(p).
struct Mix { double p, lam, lnorm; };
inline Mix mix(double p, double lam) {
  return {p, lam, p * std::log(lam) - lgamma_fast(p)};
}

// log int_a^b Gamma(m; p, lam) phi(m) dm, phi = phi_std((m - d)/z) br(m) when
// gauss, br(m) otherwise.  Nodes sit on the envelope
// e(m) = pe log m - lam m - [gauss] (m - d)^2/(2 z^2), pe = max(p - 1, 0),
// between its maximum on [a, b] and the points kDrop below it; for p < 1 the
// power is left out of the envelope, which only falls faster with it.
template <class Fn>
inline double log_integral(const Mix& g, double a, double b, bool gauss, double d,
                           double z, Fn br) {
  const double pe = std::max(g.p - 1.0, 0.0), lam = g.lam;
  const double iz2 = gauss ? 1.0 / (z * z) : 0.0;
  auto env = [&](double m) {
    return (pe > 0.0 ? pe * std::log(m) : 0.0) - lam * m - 0.5 * (m - d) * (m - d) * iz2;
  };
  auto denv = [&](double m) { return pe / m - lam - (m - d) * iz2; };
  // Unconstrained maximum: the positive root of m^2 - mu m - pe z^2 = 0,
  // mu = d - lam z^2 (gauss), or pe/lam.
  double m0;
  if (gauss) {
    const double mu = d - lam * z * z;
    const double disc = std::sqrt(mu * mu + 4.0 * pe * z * z);
    m0 = (mu >= 0.0) ? 0.5 * (mu + disc) : 2.0 * pe * z * z / (disc - mu);
  } else {
    m0 = pe / lam;
  }
  const double ms = std::min(std::max(m0, a), b);
  const double es = env(ms);
  const double curv = (ms > 0.0 ? pe / (ms * ms) : 0.0) + iz2;
  const double sig = curv > 0.0 ? 1.0 / std::sqrt(curv) : 1.0 / lam;
  const double target = es - kDrop;
  const double step0 = (curv > 0.0) ? sig * std::sqrt(2.0 * kDrop) : kDrop / lam;
  // Drop points by Newton from the far side, where it converges monotonically
  // on the concave envelope.
  double hi;
  if (R_FINITE(b) && env(b) >= target) {
    hi = b;
  } else {
    double x = ms + step0;
    for (int i = 0; i < 60 && env(x) > target; ++i) x = ms + 2.0 * (x - ms);
    if (R_FINITE(b)) x = std::min(x, b);
    for (int i = 0; i < 8; ++i) {
      const double xn = x - (env(x) - target) / denv(x);
      if (!(xn > ms && xn < x)) break;
      if (x - xn < 1e-4 * sig) { x = xn; break; }
      x = xn;
    }
    hi = x;
  }
  double lo = a;
  if (ms > a && env(a) < target) {
    double x = ms - step0;
    for (int i = 0; i < 60 && x > a && env(x) > target; ++i) x = ms - 2.0 * (ms - x);
    // log m falls faster than the quadratic the step assumed: approach a
    // geometrically until the start lies below the target.
    for (int i = 0; i < 60 && !(x > a && env(x) <= target); ++i) {
      const double t = std::ldexp(1.0, -(i + 1));
      x = a + t * (ms - a);
    }
    if (x > a && env(x) <= target) {
      for (int i = 0; i < 8; ++i) {
        const double xn = x - (env(x) - target) / denv(x);
        if (!(xn > x && xn < ms)) break;
        if (xn - x < 1e-4 * sig) { x = xn; break; }
        x = xn;
      }
      lo = x;
    }
  }
  if (!(hi > lo)) return R_NegInf;

  double ex[kNodes], bv[kNodes];
  double emax = R_NegInf;
  const double lgau = gauss ? -LOG_SQRT_2PI : 0.0;
  const JRule* jr = nullptr;
  // Jacobi carries m^{p-1} exactly; the rest of the integrand must vary by
  // little enough over [0, hi] for the polynomial degree of the rule.
  auto rest = [&](double m) { return -lam * m - 0.5 * (m - d) * (m - d) * iz2; };
  const double mr = std::min(std::max(gauss ? d - lam * z * z : 0.0, 0.0), hi);
  const bool smooth_rest = rest(mr) - std::min(rest(0.0), rest(hi)) <= 50.0;
  if (a == 0.0 && g.p != 1.0 && (g.p < 1.0 || (lo < 0.1 * hi && smooth_rest)))
    jr = jrule(g.p - 1.0);
  if (jr != nullptr) {
    // int_0^hi m^{p-1} f(m) dm = (hi/2)^p int (1 + u)^{p-1} f(hi (1 + u)/2) du
    const double lsc = g.p * std::log(0.5 * hi);
    for (int i = 0; i < kNodes; ++i) {
      const double m = 0.5 * hi * (1.0 + jr->x[i]);
      const double al = gauss ? (m - d) / z : 0.0;
      bv[i] = (m > 0.0 && jr->w[i] > 0.0) ? br(m) : 0.0;
      ex[i] = g.lnorm + lsc + std::log(jr->w[i]) - lam * m - 0.5 * al * al + lgau;
      if (bv[i] > 0.0 && ex[i] > emax) emax = ex[i];
    }
  } else {
    static thread_local const GLRule& gl = gl_get_rule(kNodes);
    const double h = 0.5 * (hi - lo), c = 0.5 * (hi + lo);
    const double lh = std::log(h);
    for (int i = 0; i < kNodes; ++i) {
      const double m = c + h * gl.x[i];
      const double al = gauss ? (m - d) / z : 0.0;
      bv[i] = (m > 0.0) ? br(m) : 0.0;
      ex[i] = g.lnorm + lh + std::log(gl.w[i]) + (g.p - 1.0) * std::log(m) - lam * m -
              0.5 * al * al + lgau;
      if (bv[i] > 0.0 && ex[i] > emax) emax = ex[i];
    }
  }
  if (!(emax > R_NegInf)) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < kNodes; ++i)
    if (bv[i] > 0.0) sum += std::exp(ex[i] - emax) * bv[i];
  return sum > 0.0 ? emax + std::log(sum) : R_NegInf;
}

// Gamma stop-losses of X ~ Gamma(p, 1), each in the regime where it is the
// small side.  For y <= p, E(y - X)+ = sum_{j>=1} j y^{p+j} e^{-y}/Gamma(p+j+1).
inline double log_stoploss_below(double p, double y) {
  if (!(y > 0.0)) return R_NegInf;
  const double l1 = (p + 1.0) * std::log(y) - y - lgamma_fast(p + 2.0);
  double t = 1.0, sum = 1.0;
  for (double j = 1.0; j < 1e6; j += 1.0) {
    const double r = (j + 1.0) / j * y / (p + j + 1.0);
    t *= r;
    sum += t;
    if (r < 1.0 && t * r < 1e-17 * sum * (1.0 - r)) break;
  }
  return l1 + std::log(sum);
}

// For y >= p, E(X - y)+ = y^p e^{-y}/Gamma(p) (1 - T)/(y + 1 - p - T), with
// T = 1(1-p)/(y+3-p - 2(2-p)/(y+5-p - ...)) the tail of Legendre's fraction
// for Gamma(p, y).  1 - T > 2/3 (T <= 0 for p >= 1), so nothing cancels.
inline double log_stoploss_above(double p, double y) {
  constexpr double tiny = 1e-300;
  double f = tiny, C = f, D = 0.0;
  for (int k = 1; k < 2000; ++k) {
    const double ak = (k == 1) ? (1.0 - p) : -(k * (k - p));
    const double bk = y + 2.0 * k + 1.0 - p;
    D = bk + ak * D;
    if (std::fabs(D) < tiny) D = tiny;
    C = bk + ak / C;
    if (std::fabs(C) < tiny) C = tiny;
    D = 1.0 / D;
    const double del = C * D;
    f *= del;
    if (std::fabs(del - 1.0) < 1e-15) break;
  }
  const double T = (std::fabs(f) < 1e-290) ? 0.0 : f;
  return p * std::log(y) - y - lgamma_fast(p) + std::log(1.0 - T) -
         std::log(y + 1.0 - p - T);
}

struct Kernel {
  double beta = 1.0;

  // F and S below and above d at shape beta; lP, lQ are log P_beta(lam d)
  // and log Q_beta(lam d).
  void point_FS(double m, double d, double s, double* lF, double* lS) const {
    const double z = s * std::sqrt(d), lam = beta / m, y = lam * d;
    const Mix g = mix(beta, lam);
    const double lIF = log_integral(g, 0.0, d, true, d, z, [&](double mm) {
      const double al = (mm - d) / z, bb = (mm + d) / z;
      return ig::mills_int(1, -al, ig::mills_seq(-al), bb, ig::mills_seq(bb), 2.0 * mm / z);
    });
    const double lIS = log_integral(g, d, R_PosInf, true, d, z, [&](double mm) {
      const double al = (mm - d) / z, bb = (mm + d) / z;
      return ig::mills_seq(al).R + ig::mills_seq(bb).R;
    });
    const double lP = R::pgamma(y, beta, 1.0, 1, 1), lQ = R::pgamma(y, beta, 1.0, 0, 1);
    *lS = log_add(log_less(lP, lIF), lIS);
    *lF = log_add(log_less(lQ, lIS), lIF);
  }
  double log_surv_point(double u, double v, double d, double s) const {
    double lF, lS;
    point_FS(v * u, d, s, &lF, &lS);
    return lS;
  }
  double log_cdf_point(double u, double v, double d, double s) const {
    double lF, lS;
    point_FS(v * u, d, s, &lF, &lS);
    return lF;
  }
  // f = v E'[(2/z) phi(alpha) (D(beta) + (m/z) R(beta))] at shape beta + 1.
  double log_pdf_point(double u, double v, double d, double s) const {
    const double m = v * u, z = s * std::sqrt(d);
    const Mix g = mix(beta + 1.0, beta / m);
    return std::log(v) + log_integral(g, 0.0, R_PosInf, true, d, z, [&](double mm) {
      const ig::Mills mb = ig::mills_seq((mm + d) / z);
      return (2.0 / z) * (mb.D + (mm / z) * mb.R);
    });
  }
  bool prefer_Jc(double, double, double, double) const { return false; }
  double log_end_M(double, double, double) const { return R_NaN; }
  double fast_cdf_surv(double, double, double, double, bool) const { return R_NaN; }
  double fast_pdf(double, double, double, double, double) const { return R_NaN; }

  EndJ end_J(double m, double d, double s, int /*want*/ = 3) const {
    if (d <= 0.0) return {0.0, R_NegInf};
    const double p = beta + 1.0, z = s * std::sqrt(d), lam = beta / m, y = lam * d;
    const double wd = 2.0 * d / z;
    const Mix g = mix(p, lam);
    auto Jc_br = [&](double mm) {
      const double al = (mm - d) / z, bb = (mm + d) / z;
      const ig::Mills mb = ig::mills_seq(bb);
      return ig::mills_int(1, al, ig::mills_seq(al), bb, mb, wd) + wd * mb.D;
    };
    const double lIJcp = log_integral(g, d, R_PosInf, true, d, z, Jc_br);
    const double lP = R::pgamma(y, p, 1.0, 1, 1), lQ = R::pgamma(y, p, 1.0, 0, 1);
    const double lQm = log_less(lQ, lIJcp);  // int_d^inf g J
    // J at m = d bounds J below d.
    const ig::Mills m0 = ig::mills_seq(0.0), mw = ig::mills_seq(wd);
    const double Jd = M_1_SQRT_2PI * (ig::mills_int(1, 0.0, m0, wd, mw, wd) + 2.0 * mw.E +
                                      wd * mw.D);
    EndJ e;
    if (Jd <= 0.6) {
      const double lIJm = log_integral(g, 0.0, d, true, d, z, [&](double mm) {
        const double al = (mm - d) / z, bb = (mm + d) / z, wm = 2.0 * mm / z;
        const ig::Mills mb = ig::mills_seq(bb);
        return ig::mills_int(1, -al, ig::mills_seq(-al), bb, mb, wm) + 2.0 * mb.E + wm * mb.D;
      });
      e.lJ = log_add(lQm, lIJm);
      e.lJc = log_add(log_less(lP, lIJm), lIJcp);
    } else {
      const double lIJcm = log_integral(g, 0.0, d, false, d, z, [&](double mm) {
        const double al = (mm - d) / z;
        return M_1_SQRT_2PI * std::exp(-0.5 * al * al) * Jc_br(mm);
      });
      e.lJc = log_add(lIJcm, lIJcp);
      e.lJ = (e.lJc < -M_LN2) ? std::log1p(-std::exp(e.lJc))
                              : log_add(lQm, log_less(lP, lIJcm));
    }
    return e;
  }

  // H - L = mbar - d; the smaller is the two Gaussian-tailed integrals plus
  // the Gamma stop-loss on its own small side.
  EndHL end_HL(double m, double d, double s) const {
    if (d <= 0.0) return {std::log(m), R_NegInf};
    const double z = s * std::sqrt(d), lam = beta / m, y = lam * d;
    const Mix g = mix(beta, lam);
    const double lIH = log_integral(g, 0.0, d, true, d, z, [&](double mm) {
      const double al = (mm - d) / z, bb = (mm + d) / z;
      return z * ig::mills_int(2, -al, ig::mills_seq(-al), bb, ig::mills_seq(bb), 2.0 * mm / z);
    });
    const double lIL = log_integral(g, d, R_PosInf, true, d, z, [&](double mm) {
      const double al = (mm - d) / z, bb = (mm + d) / z;
      return z * ig::mills_int(2, al, ig::mills_seq(al), bb, ig::mills_seq(bb), 2.0 * d / z);
    });
    const double lI = log_add(lIH, lIL);
    EndHL e;
    if (m < d) {
      e.lH = log_add(lI, log_stoploss_above(beta, y) - std::log(lam));
      e.lL = log_add(e.lH, std::log(d - m));
    } else {
      e.lL = log_add(lI, log_stoploss_below(beta, y) - std::log(lam));
      e.lH = (m > d) ? log_add(e.lL, std::log(m - d)) : e.lL;
    }
    return e;
  }
  void ends_J(double m, double d1, double d2, double s, int want, EndJ* e1, EndJ* e2) const {
    *e1 = end_J(m, d1, s, want);
    *e2 = end_J(m, d2, s, want);
  }
  void ends_HL(double m, double d1, double d2, double s, EndHL* e1, EndHL* e2) const {
    *e1 = end_HL(m, d1, s);
    *e2 = end_HL(m, d2, s);
  }
};

// The trial rate, then the fixed-rate running maximum at that rate.
inline double draw(double v, double B, double A, double s, double cv) {
  const double V = R::rgamma(1.0 / (cv * cv), v * cv * cv);
  return ig::Kernel::draw(V, B, A, s);
}

}  // namespace igg

// ===========================================================================
// Uniform start point, shared by both evidence laws
// ===========================================================================

// log of (1/A) int_B^{B+A} exp(fn(d)) dd by Gauss-Legendre in w = sqrt(d),
// which absorbs the 1/sqrt(d) of the IG point density at B = 0.  The weights
// are normalised by the rule's own integral of 2w, so that a constant
// integrand averages to itself whatever the rounding of the square roots.
template <class Fn>
inline double log_gl_average(double B, double A, Fn fn) {
  const GLRule& g = gl_get_rule(kGLNodes);
  const double w0 = std::sqrt(B), w1 = std::sqrt(B + A);
  const double h = 0.5 * (w1 - w0), c = 0.5 * (w1 + w0);
  double lmax = R_NegInf;
  double lv[kGLNodes];
  for (int i = 0; i < kGLNodes; ++i) {
    const double w = c + h * g.x[i];
    lv[i] = fn(w * w) + std::log(w * g.w[i] / (2.0 * c));  // weights sum to one
    if (lv[i] > lmax) lmax = lv[i];
  }
  if (lmax == R_NegInf) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < kGLNodes; ++i) sum += std::exp(lv[i] - lmax);
  return lmax + std::log(sum);
}

// The same average, accepted only when the rule on the whole interval agrees
// with the rule on its two halves (split in w); NaN otherwise.
template <class Fn>
inline double log_gl_checked(double B, double A, Fn fn) {
  const double whole = log_gl_average(B, A, fn);
  const double wm = 0.5 * (std::sqrt(B) + std::sqrt(B + A));
  const double dm = wm * wm, A1 = dm - B, A2 = (B + A) - dm;
  if (!(A1 > 0.0) || !(A2 > 0.0)) return R_NaN;
  const double halves = log_add(std::log(A1 / A) + log_gl_average(B, A1, fn),
                                std::log(A2 / A) + log_gl_average(dm, A2, fn));
  return std::fabs(whole - halves) < 1e-10 * std::fmax(1.0, std::fabs(halves))
    ? halves : R_NaN;
}

// log(e^lx - e^ly) for two ends, or NaN when they agree too closely for the
// difference to keep its digits.
inline double log_end_diff(double lx, double ly) {
  if (ly == R_NegInf) return lx;
  if (!(lx > ly) || ly - lx > std::log(kCancel)) return R_NaN;
  return log_diff(lx, ly);
}

// Differences that amplify the rounding of their operands by more than this
// (log scale) are replaced by quadrature where it is available.
constexpr double kLogAmpMax = 4.1588830833596715;  // log(64)

// log F from A F = M(B+A) - M(B), M = E[min(G, d)]; NaN where the kernel has
// no M form.  *lamp receives the log of the factor M(B+A) / (A F) by which
// the difference amplifies the rounding of its ends.
template <class K>
inline double end_M_diff(const K& k, double m, double B, double A, double s, double lA,
                         double* lamp) {
  *lamp = R_PosInf;
  const double l2 = k.log_end_M(m, B + A, s);
  if (ISNAN(l2)) return R_NaN;
  const double l1 = k.log_end_M(m, B, s);
  if (!(l2 > l1)) return R_NaN;
  const double n = log_diff(l2, l1);
  *lamp = l2 - n;
  return clamp_log_prob(n - lA);
}

// log F and log S.  A S = L(B+A) - L(B) and A F = H(B) - H(B+A); each side
// keeps relative precision below one half and the other is its complement.
// A null lF or lS skips the quadratures that only that output needs.
template <class K>
inline void log_cdf_surv(const K& k, double u, double v, double B, double A, double s,
                         double* lF, double* lS) {
  double lF_, lS_;
  const bool want_F = lF != nullptr, want_S = lS != nullptr;
  if (!want_F) lF = &lF_;
  if (!want_S) lS = &lS_;
  A = (B + A) - B;  // the width of the interval the ends actually span
  if (A <= 0.0) {
    if (!want_S) {
      *lF = k.log_cdf_point(u, v, B, s);
      return;
    }
    *lS = k.log_surv_point(u, v, B, s);
    if (want_F)
      *lF = (*lS < -M_LN2) ? std::log(-std::expm1(*lS)) : k.log_cdf_point(u, v, B, s);
    return;
  }
  const double m = v * u, lA = std::log(A);
  if (want_F != want_S) {
    const double q = k.fast_cdf_surv(m, B, A, s, want_S);
    if (!ISNAN(q)) { *(want_S ? lS : lF) = q; return; }
  }
  EndHL e1, e2;
  k.ends_HL(m, B, B + A, s, &e1, &e2);
  double nS = log_end_diff(e2.lL, e1.lL), nF = log_end_diff(e1.lH, e2.lH);
  if (!ISNAN(nS)) nS -= lA;
  if (!ISNAN(nF)) nF -= lA;
  if (!ISNAN(nS) && nS < -M_LN2) { *lS = nS; *lF = std::log(-std::expm1(nS)); return; }
  if (!ISNAN(nF) && nF < -M_LN2) { *lF = nF; *lS = std::log(-std::expm1(nF)); return; }
  if (!ISNAN(nS)) {
    *lS = clamp_log_prob(nS);
    *lF = ISNAN(nF) ? std::log(-std::expm1(*lS)) : nF;
    // Where the H ends cancel, F is the complement of S or A F = M(B+A) - M(B)
    // with M = E[min(G, d)] = m - H, whose ends need not agree when the H
    // ends do (M(0) = 0), whichever amplifies rounding less.  Where both lose
    // more than kLogAmpMax, the point CDF, which then varies little over the
    // interval, is averaged instead.  S >= 1/2 is retaken as the complement
    // of F, which is the more accurate of the two and keeps F + S = 1
    // whichever output was asked for.
    if (ISNAN(nF)) {
      // log of 1 / (F (1 - L(B)/L(B+A))), the complement's amplification
      const double lcomp = -(*lF + std::log1p(-std::exp(e1.lL - e2.lL)));
      double lamp;
      double lq = end_M_diff(k, m, B, A, s, lA, &lamp);
      if (!(lamp < lcomp)) lq = R_NaN;
      if (std::fmin(lamp, lcomp) > kLogAmpMax) {
        const double lg = log_gl_checked(B, A, [&](double d) {
          return k.log_cdf_point(u, v, d, s); });
        if (!ISNAN(lg)) lq = lg;
      }
      if (!ISNAN(lq)) { *lF = lq; *lS = std::log1p(-std::exp(lq)); }
    }
    return;
  }
  if (!want_S && !ISNAN(nF)) { *lF = clamp_log_prob(nF); return; }
  // The survivor ends agree too closely: average the point survivor, whose
  // integrand then varies little over the interval.
  *lS = clamp_log_prob(log_gl_average(B, A, [&](double d) {
    return k.log_surv_point(u, v, d, s); }));
  if (!ISNAN(nF)) *lF = clamp_log_prob(nF);
  else if (*lS < -M_LN2) *lF = std::log(-std::expm1(*lS));
  else if (want_F) {
    double lamp;
    *lF = end_M_diff(k, m, B, A, s, lA, &lamp);
    if (!(lamp <= kLogAmpMax))
      *lF = log_gl_average(B, A, [&](double d) { return k.log_cdf_point(u, v, d, s); });
  }
}

template <class K>
inline double log_surv(const K& k, double u, double v, double B, double A, double s) {
  if (!(u > 0.0) || v == 0.0) return 0.0;
  if (!R_FINITE(u)) return R_NegInf;
  double lS;
  log_cdf_surv(k, u, v, B, A, s, nullptr, &lS);
  return clamp_log_prob(lS);
}

template <class K>
inline double log_cdf(const K& k, double u, double v, double B, double A, double s) {
  if (!(u > 0.0) || v == 0.0) return R_NegInf;
  if (!R_FINITE(u)) return 0.0;
  double lF;
  log_cdf_surv(k, u, v, B, A, s, &lF, nullptr);
  return clamp_log_prob(lF);
}

template <class K>
inline double log_surv(double u, double v, double B, double A, double s) {
  return log_surv(K{}, u, v, B, A, s);
}
template <class K>
inline double log_cdf(double u, double v, double B, double A, double s) {
  return log_cdf(K{}, u, v, B, A, s);
}

// f = v [J(B) - J(B+A)] / A = v [Jc(B+A) - Jc(B)] / A.  Where an end carries
// both forms the one whose larger term is smaller is differenced; otherwise
// the kernel's preferred form is formed first and the other only when the
// first agrees too closely at the two ends.
template <class K>
inline double log_pdf(const K& k, double u, double v, double B, double A, double s) {
  if (!(u > 0.0) || !R_FINITE(u) || v == 0.0) return R_NegInf;
  A = (B + A) - B;
  if (A <= 0.0) return k.log_pdf_point(u, v, B, s);
  const double m = v * u;
  const double fast = k.fast_pdf(m, v, B, A, s);
  if (!ISNAN(fast)) return fast;
  const int want = k.prefer_Jc(m, B, A, s) ? 2 : 1;
  EndJ e1, e2;
  k.ends_J(m, B, B + A, s, want, &e1, &e2);
  const bool use_Jc = (ISNAN(e1.lJ) || ISNAN(e2.lJ) || ISNAN(e1.lJc) || ISNAN(e2.lJc))
    ? want == 2 : e2.lJc < e1.lJ;
  double nJ = use_Jc ? log_end_diff(e2.lJc, e1.lJc) : log_end_diff(e1.lJ, e2.lJ);
  if (ISNAN(nJ)) {
    if (ISNAN(use_Jc ? e1.lJ : e1.lJc) || ISNAN(use_Jc ? e2.lJ : e2.lJc)) {
      k.ends_J(m, B, B + A, s, use_Jc ? 1 : 2, &e1, &e2);
    }
    nJ = use_Jc ? log_end_diff(e1.lJ, e2.lJ) : log_end_diff(e2.lJc, e1.lJc);
  }
  if (!ISNAN(nJ)) return std::log(v) + nJ - std::log(A);
  return log_gl_average(B, A, [&](double d) { return k.log_pdf_point(u, v, d, s); });
}

template <class K>
inline double log_pdf(double u, double v, double B, double A, double s) {
  return log_pdf(K{}, u, v, B, A, s);
}

template <class K>
inline double draw_decision(double v, double B, double A, double s) {
  if (!pars_ok(v, B, A, s) || v == 0.0) return R_PosInf;
  return K::draw(v, B, A, s);
}

}  // namespace rsub

double drsub_scalar(double t, const double* par, void* ctx_);
double prsub_scalar(double t, const double* par, void* ctx_);
void drsub_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void prsub_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void rsub_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);

#endif

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

// Two ends whose ratio exceeds this are differenced by quadrature instead.
constexpr double kCancel = 0.75;
constexpr int kGLNodes = 8;
// log of the smallest probability taken as the complement of its partner.
constexpr double kComplementFloor = -9.210340371976184;  // log(1e-4)

inline bool pars_ok(double v, double B, double A, double s) {
  return R_FINITE(v) && v > 0.0 && R_FINITE(B) && B >= 0.0 &&
         R_FINITE(A) && A >= 0.0 && B + A > 0.0 && R_FINITE(s) && s > 0.0;
}

inline double log_diff(double lx, double ly) {  // log(e^lx - e^ly), lx > ly
  return lx + std::log1p(-std::exp(ly - lx));
}

inline double log_add(double lx, double ly) {
  if (lx < ly) std::swap(lx, ly);
  if (ly == R_NegInf) return lx;
  return lx + std::log1p(std::exp(ly - lx));
}

struct EndJ { double lJ, lJc; };
struct EndHL { double lH, lL; };

// ===========================================================================
// Inverse-Gaussian subordinator
// ===========================================================================
namespace ig {

// The Mills ratio R(x) = Phi(-x)/phi(x) and its successive negative
// derivatives D = 1 - xR and E = R - xD.  Below x = 4 they are formed
// directly, which costs at most a few hundred ulps in E; above it the
// subtractions are exact in the tails of the Laplace continued fraction
// R = 1/(x + c1), c_k = k/(x + c_{k+1}), as D = c1 R and E = c1 c2 R.
struct Mills { double R, D, E; };

inline Mills mills_seq(double x) {
  if (x < 4.0) {
    const double R = mills_ratio_std(x);
    const double D = 1.0 - x * R;
    return {R, D, R - x * D};
  }
  if (!R_FINITE(x)) return {0.0, 0.0, 0.0};
  const int K = x < 6.0 ? 40 : (x < 10.0 ? 24 : (x < 20.0 ? 16 : 10));
  double c = 0.0, c2 = 0.0;
  for (int k = K; k >= 1; --k) {
    if (k == 1) c2 = c;
    c = k / (x + c);
  }
  const double R = 1.0 / (x + c);
  return {R, c * R, c * c2 * R};
}

// log of int_a^b D (order 1, = R(a) - R(b)) or int_a^b E (order 2,
// = D(a) - D(b)), a < b.
inline double log_mills_int(int order, double a, double b) {
  if (!(b > a)) return R_NegInf;
  const Mills ma = mills_seq(a), mb = mills_seq(b);
  const double xa = order == 1 ? ma.R : ma.D;
  const double xb = order == 1 ? mb.R : mb.D;
  if (xb < kCancel * xa) return std::log(xa - xb);
  const GLRule& g = gl_get_rule(kGLNodes);
  const double h = 0.5 * (b - a), c = 0.5 * (b + a);
  double sum = 0.0;
  for (int i = 0; i < kGLNodes; ++i) {
    const Mills mi = mills_seq(c + h * g.x[i]);
    sum += g.w[i] * (order == 1 ? mi.D : mi.E);
  }
  return std::log(h * sum);
}

// Coordinates of the distance-d first passage at mean evidence m = v u.
struct Pt {
  double z, alpha, beta, lphi;
};

inline Pt point(double m, double d, double s) {
  Pt p;
  p.z = s * std::sqrt(d);
  p.alpha = (m - d) / p.z;
  p.beta = (m + d) / p.z;
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
    return p.lphi + log_mills_int(1, -p.alpha, p.beta);
  }

  // Each end is carried in both complementary forms, each written as a
  // positive sum wherever its terms stay in range (|alpha| < kAlphaMax) and
  // as the complement of the other form beyond.
  static constexpr double kAlphaMax = 30.0;

  static EndJ end_J(double m, double d, double s) {
    if (d <= 0.0) return {0.0, R_NegInf};  // J(0) = 1
    const Pt p = point(m, d, s);
    const Mills mb = mills_seq(p.beta);
    EndJ e;
    // J = phi(alpha) [int_{-alpha}^{beta} D + 2E(beta) + (alpha + beta) D(beta)]
    //   = Phi(alpha) + phi(alpha) [E(beta) + alpha D(beta)]   (alpha >= 0).
    if (p.alpha < 0.0) {
      e.lJ = p.lphi + log_add(log_mills_int(1, -p.alpha, p.beta),
                              std::log(2.0 * mb.E + (p.alpha + p.beta) * mb.D));
    } else {
      e.lJ = log_add(pnorm_std(p.alpha, true, true),
                     p.lphi + std::log(mb.E + p.alpha * mb.D));
    }
    // Jc = phi(alpha) [int_alpha^beta D + (beta - alpha) D(beta)].
    if (p.alpha > -kAlphaMax) {
      e.lJc = p.lphi + log_add(log_mills_int(1, p.alpha, p.beta),
                               std::log((p.beta - p.alpha) * mb.D));
    } else {
      e.lJc = std::log(-std::expm1(std::fmin(e.lJ, 0.0)));
    }
    return e;
  }

  static EndHL end_HL(double m, double d, double s) {
    if (d <= 0.0) return {std::log(m), R_NegInf};  // H(0) = m, L(0) = 0
    const Pt p = point(m, d, s);
    const double lz = std::log(p.z);
    EndHL e;
    // H = z phi(alpha) int_{-alpha}^{beta} E,  L = z phi(alpha) int_alpha^beta E.
    e.lH = p.alpha < kAlphaMax ? lz + p.lphi + log_mills_int(2, -p.alpha, p.beta)
                               : R_NaN;
    e.lL = p.alpha > -kAlphaMax ? lz + p.lphi + log_mills_int(2, p.alpha, p.beta)
                                : R_NaN;
    if (ISNAN(e.lH)) e.lH = log_add(std::log(m - d), e.lL);
    if (ISNAN(e.lL)) e.lL = log_add(std::log(d - m), e.lH);
    return e;
  }

  // Distance d = B + A U, then (s/v) times the maximum over [0, d] of
  // W + w/s, drawn from its endpoint and the Brownian-bridge maximum.
  static double draw(double v, double B, double A, double s) {
    const double d = B + A * R::unif_rand();
    const double x = d / s + std::sqrt(d) * R::norm_rand();
    const double mx = 0.5 * (x + std::sqrt(x * x - 2.0 * d * std::log(R::unif_rand())));
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
constexpr double kRescale = 1e250;
constexpr double kLogRescale = 575.6462732485114;  // log(1e250)

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
// the upper side, b/a for the lower.  Returns the logs of sum_{k>=0} q_k,
// sum_{k>=1} q_k and sum_{k>=1} k q_k.
struct Side { double s0, s1, sk; };

// The ratios come from the stable backward recurrence
// rho_k = 1/(2(k+1)/z + rho_{k+1}), seeded at Amos' approximation far enough
// above the last index used, kstop, that the seed's error has decayed by
// prod rho_k^2 ~ exp(-(K^2 - kstop^2)/z).  kstop comes from the Gaussian
// tail of the Skellam law and kmode is the index of the side's largest term.
inline Side side_sums(double lp0, double z, double t, double kmode, double kstop) {
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
    double L = lp0, term = 1.0, s0 = 1.0, s1 = 0.0, sk = 0.0, prev = 1.0;
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
      return {L + std::log(s0), s1 > 0.0 ? L + std::log(s1) : R_NegInf,
              sk > 0.0 ? L + std::log(sk) : R_NegInf};
    }
  }
  return {R_NaN, R_NaN, R_NaN};
}

// The Skellam quantities of one (a, b), each on its small side.
struct Skel { double lS, lF, lJ, lJc, lEpos, lEneg; };

inline Skel skel(double a, double b) {
  const double z = 2.0 * std::sqrt(a * b);
  const double sd = std::sqrt(a) - std::sqrt(b);
  const double lp0 = -sd * sd + log_i0_scaled(z);
  const double gap = std::fabs(a - b);
  // Terms beyond k with (k -/+ gap)^2 - gap^2 > 82 (a + b) are below e^-41
  // of the largest.
  const double tail2 = 82.0 * (a + b);
  Skel o;
  // The side away from the mean has decreasing terms; the other is only
  // summed when the first carries nearly all the mass, so that its complement
  // would lose more than three digits.
  const bool up_mean = a >= b;
  const Side lo = side_sums(lp0, z, up_mean ? b / a : a / b, 0.0,
                            std::sqrt(gap * gap + tail2) - gap);
  Side hi = {R_NaN, R_NaN, R_NaN};
  if (lo.s0 > kSideFull)
    hi = side_sums(lp0, z, up_mean ? a / b : b / a, gap, gap + std::sqrt(tail2));
  // lo is P(Delta <= 0), P(Delta <= -1), E[Delta-] when a >= b, and
  // P(Delta >= 0), P(Delta >= 1), E[Delta+] otherwise.
  double near0, near1, nearE, far0, far1, farE;
  near0 = lo.s0; near1 = lo.s1; nearE = lo.sk;
  if (ISNAN(hi.s0)) {
    far0 = log1m_exp(std::fmin(near1, 0.0));  // complement excludes 0
    far1 = log1m_exp(std::fmin(near0, 0.0));
    farE = log_add(std::log(gap), nearE);
  } else {
    far0 = hi.s0; far1 = hi.s1; farE = hi.sk;
  }
  if (up_mean) {
    o.lS = near0; o.lJc = near1; o.lEneg = nearE;
    o.lJ = far0; o.lF = far1; o.lEpos = farE;
  } else {
    o.lJ = near0; o.lF = near1; o.lEpos = nearE;
    o.lS = far0; o.lJc = far1; o.lEneg = farE;
  }
  return o;
}

struct Kernel {
  static double log_pdf_point(double u, double v, double d, double s) {
    const double c = s * s, a = 2.0 * v * u / c, b = 2.0 * d / c;
    const double sd = std::sqrt(a) - std::sqrt(b);
    return std::log(2.0 * v / c) - sd * sd + log_i0_scaled(2.0 * std::sqrt(a * b));
  }
  static double log_surv_point(double u, double v, double d, double s) {
    const double c = s * s;
    return skel(2.0 * v * u / c, 2.0 * d / c).lS;
  }
  static double log_cdf_point(double u, double v, double d, double s) {
    const double c = s * s;
    return skel(2.0 * v * u / c, 2.0 * d / c).lF;
  }
  static EndJ end_J(double m, double d, double s) {
    if (d <= 0.0) return {0.0, R_NegInf};
    const double c = s * s;
    const Skel k = skel(2.0 * m / c, 2.0 * d / c);
    return {k.lJ, k.lJc};
  }
  static EndHL end_HL(double m, double d, double s) {
    if (d <= 0.0) return {std::log(m), R_NegInf};
    const double c = s * s, lth = std::log(0.5 * c);
    const Skel k = skel(2.0 * m / c, 2.0 * d / c);
    return {lth + k.lEpos, lth + k.lEneg};
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

// log F and log S.  A S = L(B+A) - L(B) and A F = H(B) - H(B+A); each side
// keeps relative precision below one half and the other is its complement.
template <class K>
inline void log_cdf_surv(double u, double v, double B, double A, double s,
                         double* lF, double* lS) {
  A = (B + A) - B;  // the width of the interval the ends actually span
  if (A <= 0.0) {
    *lS = K::log_surv_point(u, v, B, s);
    *lF = (*lS < -M_LN2) ? std::log(-std::expm1(*lS)) : K::log_cdf_point(u, v, B, s);
    return;
  }
  const double m = v * u, lA = std::log(A);
  const EndHL e1 = K::end_HL(m, B, s), e2 = K::end_HL(m, B + A, s);
  double nS = log_end_diff(e2.lL, e1.lL), nF = log_end_diff(e1.lH, e2.lH);
  if (!ISNAN(nS)) nS -= lA;
  if (!ISNAN(nF)) nF -= lA;
  if (!ISNAN(nS) && nS < -M_LN2) { *lS = nS; *lF = std::log(-std::expm1(nS)); return; }
  if (!ISNAN(nF) && nF < -M_LN2) { *lF = nF; *lS = std::log(-std::expm1(nF)); return; }
  if (!ISNAN(nS)) {
    *lS = std::fmin(nS, 0.0);
    *lF = ISNAN(nF) ? std::log(-std::expm1(*lS)) : nF;
    // The H ends agree closely, so the point CDF varies little over the
    // interval; average it where the complement has lost too many digits.
    if (ISNAN(nF) && *lF < kComplementFloor) {
      const double lq = log_gl_checked(B, A, [&](double d) {
        return K::log_cdf_point(u, v, d, s); });
      if (!ISNAN(lq)) *lF = lq;
    }
    return;
  }
  // The survivor ends agree too closely: average the point survivor, whose
  // integrand then varies little over the interval.
  *lS = std::fmin(log_gl_average(B, A, [&](double d) {
    return K::log_surv_point(u, v, d, s); }), 0.0);
  if (!ISNAN(nF)) *lF = std::fmin(nF, 0.0);
  else if (*lS < -M_LN2) *lF = std::log(-std::expm1(*lS));
  else *lF = log_gl_average(B, A, [&](double d) { return K::log_cdf_point(u, v, d, s); });
}

template <class K>
inline double log_surv(double u, double v, double B, double A, double s) {
  if (!(u > 0.0)) return 0.0;
  if (!R_FINITE(u)) return R_NegInf;
  double lF, lS;
  log_cdf_surv<K>(u, v, B, A, s, &lF, &lS);
  return std::fmin(lS, 0.0);
}

template <class K>
inline double log_cdf(double u, double v, double B, double A, double s) {
  if (!(u > 0.0)) return R_NegInf;
  if (!R_FINITE(u)) return 0.0;
  double lF, lS;
  log_cdf_surv<K>(u, v, B, A, s, &lF, &lS);
  return std::fmin(lF, 0.0);
}

// f = v [J(B) - J(B+A)] / A = v [Jc(B+A) - Jc(B)] / A.
template <class K>
inline double log_pdf(double u, double v, double B, double A, double s) {
  if (!(u > 0.0) || !R_FINITE(u)) return R_NegInf;
  A = (B + A) - B;
  if (A <= 0.0) return K::log_pdf_point(u, v, B, s);
  const double m = v * u;
  const EndJ e1 = K::end_J(m, B, s), e2 = K::end_J(m, B + A, s);
  const double nJ = (e2.lJc < e1.lJ) ? log_end_diff(e2.lJc, e1.lJc)
                                     : log_end_diff(e1.lJ, e2.lJ);
  if (!ISNAN(nJ)) return std::log(v) + nJ - std::log(A);
  return log_gl_average(B, A, [&](double d) { return K::log_pdf_point(u, v, d, s); });
}

template <class K>
inline double draw_decision(double v, double B, double A, double s) {
  if (!pars_ok(v, B, A, s)) return R_PosInf;
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

#ifndef EMC2_MODEL_OPR_H
#define EMC2_MODEL_OPR_H

// Opponent Poisson race (Math/CounterAndQuorumModels.tex, sec:opponent).
//
// Accumulator evidence X(u) = N+(u) - N-(u), with independent Poisson streams
// of rates a = lambda+ > 0 and b = lambda- >= 0, finishing at the continuous
// criterion K > 0.  The first-passage subdensity is
//
//   f(u) = (K/u) (a/b)^(K/2) exp(-(a+b)u) I_K(2u sqrt(ab))
//        = (K/u) sum_j Pois(j; bu) g(K+j; au),   g(s; x) = x^s e^-x / Gamma(s+1),
//
// and the survivor, from F = H(2au; 2K, 2bu) + (a/b)^K H(2bu; 2K+2, 2au) with
// H the noncentral chi-squared CDF written as a Poisson mixture of regularized
// incomplete gammas, is
//
//   S(u) = A - B,  A = sum_j Pois(j; bu) Q(K+j, au),
//                  B = (a/b)^K sum_j Pois(j; au) P(K+1+j, bu).
//
// Both sums are formed by one-directional recurrences that only add positive
// terms: Q(s+1) = Q(s) + g(s) upward in j for A, P(s) = P(s+1) + g(s)
// downward in j for B.  The summands are log-concave in j, so each sweep
// starts beyond its peak on the side the recurrence leaves and stops once the
// terms past the peak are below rounding.  Everything is carried on a running
// log scale, so neither sum under- or overflows in the far tail.
//
// Below its median the CDF is the positive sum
// F = sum_j Pois(j; bu) P(K+j, au) + B, formed by the same downward
// recurrence as B, so log F keeps its relative digits at the leading edge
// where 1 - S would round to zero.
//
// Where z = 2u sqrt(ab) is large enough for the Hankel expansion of I_K,
// that expansion is integrated term by term from u to Inf instead,
//
//   S(u) = K (a/b)^(K/2) e^{-cu} (2 pi z)^{-1/2} sum_k (-1)^k a_k(K) z^-k e^{cu} E_{3/2+k}(cu),
//
// c = (sqrt a - sqrt b)^2 and E_p the generalised exponential integral: a
// fixed handful of terms however large the rates, and no A - B cancellation
// near a = b.  Where b u is too small to move S or F by a rounding unit the
// accumulator is the Gamma(K, a) first passage.
//
// The survivor and CDF sums are used only for a >= b.  For a < b the completion
// probability is h = (a/b)^K and f(u; a, b) = h f(u; b, a), so
// S(u; a, b) = (1 - h) + h S(u; b, a) keeps the defective atom exact and
// the subtraction A - B always in the proper regime.
//
// b = 0 is the Gamma(K, a) first passage.
//
// The u-free quantities of an accumulator live in State, which the race
// adapters memoise per row (StateMemo) so that the rows of a trial and of a
// condition share them.

#include <Rcpp.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <utility>
#include "utility_functions.h"
#include "composite_functions.h"
#include "inc_gamma.h"

namespace opr {

constexpr double kTermEps = 1e-18;       // relative cutoff for a series term
constexpr double kRescale = 1e200;       // running-scale renormalisation point
constexpr double kLogRescale = 460.51701859880916;  // log(1e200)
constexpr double kWidth = 9.0;           // sweep margin, in Poisson SDs
constexpr int kMaxTerms = 200000;

// log j! for integer j.
inline double log_fact(double j) {
  static const std::array<double, 256> table = [] {
    std::array<double, 256> t{};
    for (int i = 0; i < 256; ++i) t[i] = R::lgammafn(i + 1.0);
    return t;
  }();
  return (j < 256.0) ? table[static_cast<int>(j)] : lgamma_fast(j + 1.0);
}

// log Pois(j; lam) and log g(s; x) = log(x^s e^-x / Gamma(s + 1)), the
// latter given lg_s1 = log Gamma(s + 1).  Large arguments go through R's
// saddle-point forms, which avoid the cancellation between -lam, j log lam
// and log j! that costs ~|lam| eps absolute.
constexpr double kDirect = 25.0;
inline double log_pois(double j, double lam, double log_lam) {
  if (j < kDirect && lam < kDirect) return -lam + j * log_lam - log_fact(j);
  return R::dpois(j, lam, 1);
}
inline double log_g(double s, double x, double log_x, double lg_s1) {
  if (s < kDirect && x < kDirect) return s * log_x - x - lg_s1;
  return R::dgamma(x, s + 1.0, 1.0, 1);
}

// Valid rates and criterion.  K > 0, a > 0, b >= 0, all finite.
inline bool pars_ok(double K, double a, double b) {
  return K > 0.0 && a > 0.0 && b >= 0.0 && R_FINITE(K) && R_FINITE(a) &&
         R_FINITE(b);
}

// log Gamma at the shapes K + n.  They repeat across the rows of a condition
// and change per row under a trend on K, so a direct-mapped cache is filled
// by the Stirling form: a hit costs a hash, a miss little more than the
// uncached evaluation.
inline double lgam_shape(double s) {
  if (!(s >= 0.5) || !R_FINITE(s)) return R::lgammafn(s);
  static constexpr int SLOTS = 256;  // power of two: the hash masks
  struct Slot { double key; double val; };
  static thread_local Slot tab[SLOTS] = {};
  std::uint64_t h;
  std::memcpy(&h, &s, sizeof(h));
  h ^= h >> 33; h *= 0xff51afd7ed558ccdull;
  h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ull;
  h ^= h >> 33;
  Slot& e = tab[h & (SLOTS - 1)];
  if (e.key == s) return e.val;
  e.key = s;
  e.val = lgamma_pos(s);
  return e.val;
}

// The u-free part of one accumulator.
struct State {
  bool ok = false;
  double K = 0.0, a = 0.0, b = 0.0;
  double log_K = 0.0, log_a = 0.0, log_b = R_NegInf;
  double lgK = 0.0;      // log Gamma(K)
  double log_h = 0.0;    // log P(U < Inf)
};

inline State make_state(double K, double a, double b) {
  State s;
  s.K = K; s.a = a; s.b = b;
  s.ok = pars_ok(K, a, b);
  if (!s.ok) return s;
  s.log_K = std::log(K);
  s.log_a = std::log(a);
  s.log_b = (b > 0.0) ? std::log(b) : R_NegInf;
  s.lgK = lgam_shape(K);
  s.log_h = (a >= b) ? 0.0 : K * (s.log_a - s.log_b);
  return s;
}

// The same accumulator with its rates exchanged.
inline State swapped(const State& s) {
  State w = s;
  std::swap(w.a, w.b);
  std::swap(w.log_a, w.log_b);
  w.log_h = (w.a >= w.b) ? 0.0 : w.K * (w.log_a - w.log_b);
  return w;
}

// log Gamma(K + n) for integer n >= 0, from lgK at the shapes that recur.
inline double lgam_K(const State& s, double n) {
  if (n == 0.0) return s.lgK;
  if (n == 1.0) return s.lgK + s.log_K;
  return lgam_shape(s.K + n);
}

// log sum_j Pois(j; lam) g(K+j; x) by a peak-centred sweep; the terms are
// T_j proportional to (lam x)^j / (j! Gamma(K+j+1)).
inline double log_density_series(const State& s, double lam, double x,
                                 double log_lam, double log_x) {
  const double K = s.K;
  const double lx = lam * x;
  const double jp = std::floor(0.5 * (std::sqrt(K * K + 4.0 * lx) - K));
  const double log_tp = log_pois(jp, lam, log_lam) +
                        log_g(K + jp, x, log_x, lgam_K(s, jp + 1.0));
  double sum = 1.0;
  double term = 1.0;
  int n = 0;
  bool done = false;
  for (double j = jp; n < kMaxTerms; j += 1.0, ++n) {
    term *= lx / ((j + 1.0) * (K + j + 1.0));
    sum += term;
    if (term < kTermEps * sum) { done = true; break; }
  }
  if (!done) return R_NaN;
  term = 1.0;
  n = 0;
  for (double j = jp; j > 0.0; j -= 1.0, ++n) {
    if (n == kMaxTerms) return R_NaN;
    term *= j * (K + j) / lx;
    sum += term;
    if (term < kTermEps * sum) break;
  }
  return log_tp + std::log(sum);
}

// log(e^-z I_nu(z)) by the Hankel expansion, for z large against nu^2.
inline double log_bessel_i_scaled_hankel(double nu, double z) {
  const double mu = 4.0 * nu * nu;
  double term = 1.0, sum = 1.0, prev = 1.0;
  for (int k = 1; k < 64; ++k) {
    const double odd = 2.0 * k - 1.0;
    term *= -(mu - odd * odd) / (8.0 * k * z);
    const double at = std::fabs(term);
    if (at >= prev) break;
    sum += term;
    if (at < 1e-17 * std::fabs(sum)) break;
    prev = at;
  }
  return std::log(sum) - 0.5 * std::log(2.0 * M_PI * z);
}

// Bessel argument above which the Hankel expansion replaces the series.
inline double hankel_threshold(double K) { return std::fmax(150.0, 4.0 * K * K); }

// (sqrt a - sqrt b)^2, without the cancellation of the root difference.
inline double rate_gap(double a, double b) {
  const double r = (a - b) / (std::sqrt(a) + std::sqrt(b));
  return r * r;
}

// e^x E_p(x) for p = 3/2 + k, k = 0..n-1, x >= 0.  Each recurrence is run in
// its stable direction: e^x E_{p+1} = (1 - x e^x E_p)/p upward while p > x,
// and its inverse downward while p < x, from a seed at p ~ x by the Legendre
// continued fraction (or, for x < 1, from E_{3/2} in closed form).
inline void scaled_expint_half(double x, int n, double* out) {
  if (x == 0.0) {
    for (int k = 0; k < n; ++k) out[k] = 1.0 / (0.5 + k);
    return;
  }
  int k0 = 0;
  if (x < 1.0) {
    // e^x E_{3/2}(x) = 2 - 2 sqrt(pi x) e^x erfc(sqrt x).
    const double erfcx = 2.0 * std::exp(x) * R::pnorm(-std::sqrt(2.0 * x), 0.0, 1.0, 1, 0);
    out[0] = 2.0 - 2.0 * std::sqrt(M_PI * x) * erfcx;
  } else {
    k0 = static_cast<int>(std::floor(x - 1.5));
    if (k0 < 0) k0 = 0;
    if (k0 > n - 1) k0 = n - 1;
    const double p = 1.5 + k0;
    double b = x + p, c = 1e300, d = 1.0 / b, h = d;
    for (int i = 1; i < 500; ++i) {
      const double an = -i * (p - 1.0 + i);
      b += 2.0;
      d = 1.0 / (an * d + b);
      c = b + an / c;
      const double del = c * d;
      h *= del;
      if (std::fabs(del - 1.0) < 3e-16) break;  // clear of one ulp of 1.0
    }
    out[k0] = h;
  }
  for (int k = k0 - 1; k >= 0; --k) out[k] = (1.0 - (1.5 + k) * out[k + 1]) / x;
  for (int k = k0 + 1; k < n; ++k) out[k] = (1.0 - x * out[k - 1]) / (0.5 + k);
}

// log S(u) of a proper accumulator (a >= b > 0) from the Hankel expansion of
// I_K integrated over (u, Inf), for z = 2u sqrt(ab) above hankel_threshold(K),
// where it holds at every later time too.  The series is cut where the Hankel
// terms themselves stop falling or drop below rounding; the exponential
// integrals only shrink with k, so the integrated terms are smaller still.
inline double log_surv_hankel(double u, const State& s) {
  const double K = s.K;
  const double z = 2.0 * u * std::sqrt(s.a * s.b);
  const double cu = rate_gap(s.a, s.b) * u;
  constexpr int kMax = 64;
  double t[kMax];
  const double mu = 4.0 * K * K;
  t[0] = 1.0;
  int n = 1;
  double tsum = 1.0;
  for (; n < kMax; ++n) {
    const double odd = 2.0 * n - 1.0;
    const double next = -t[n - 1] * (mu - odd * odd) / (8.0 * n * z);
    if (std::fabs(next) >= std::fabs(t[n - 1])) break;
    t[n] = next;
    tsum += next;
    if (std::fabs(next) < 1e-17 * std::fabs(tsum)) { ++n; break; }
  }
  double e[kMax];
  scaled_expint_half(cu, n, e);
  double sum = 0.0;
  for (int k = n - 1; k >= 0; --k) sum += t[k] * e[k];
  if (!(sum > 0.0)) return R_NaN;
  return s.log_K + 0.5 * K * (s.log_a - s.log_b) - cu -
         0.5 * std::log(2.0 * M_PI * z) + std::log(sum);
}

// True where b u is too small to change S or F of the Gamma(K, a) first
// passage by a rounding unit: a negative event before u has probability
// ~bu and moves the tail by a relative ~bu(1 + au).
inline bool gamma_limit(double u, const State& s) {
  return s.b * u * (1.0 + s.a * u) < 1e-18;
}

// log f(u) for the decision time u > 0.
inline double log_pdf(double u, const State& s) {
  if (!(u > 0.0)) return R_NegInf;
  if (u == R_PosInf) return R_NegInf;
  const double x = s.a * u;
  const double log_u = std::log(u);
  if (s.b == 0.0) {
    return s.log_a + log_g(s.K - 1.0, x, s.log_a + log_u, s.lgK);
  }
  const double lam = s.b * u;
  const double z = 2.0 * u * std::sqrt(s.a * s.b);
  if (z > hankel_threshold(s.K)) {
    return s.log_K - log_u + 0.5 * s.K * (s.log_a - s.log_b) -
           rate_gap(s.a, s.b) * u + log_bessel_i_scaled_hankel(s.K, z);
  }
  return s.log_K - log_u +
         log_density_series(s, lam, x, s.log_b + log_u, s.log_a + log_u);
}

// log of one summand of A, Pois(j; lam) Q(K+j, x).
inline double log_term_A(double j, const State& s, double lam, double x,
                         double log_lam, double log_x) {
  double lq;
  gamma_log_pq_fast(s.K + j, x, log_x, lgam_K(s, j), nullptr, &lq);
  return log_pois(j, lam, log_lam) + lq;
}

// Start of the upward sweep for A.  The summand peaks between the Poisson
// mode lam and jq, the peak of Pois(j; lam) g(K-1+j; x) that it follows
// while Q(K+j, x) is small.  When that range is wide the peak is located by
// an integer ternary search (the summand is log-concave in j), and the start
// is stepped down until its term is negligible against the peak.
inline double start_A(const State& s, double lam, double x, double log_lam,
                      double log_x) {
  const double K0 = s.K - 1.0;
  const double jq = 0.5 * (std::sqrt(K0 * K0 + 4.0 * lam * x) - K0);
  double lo = std::floor(lam), hi = std::ceil(std::fmax(lam, jq)) + 1.0;
  if (hi - lo < 64.0)
    return std::fmax(0.0, std::floor(lam - kWidth * std::sqrt(lam) - 5.0));
  while (hi - lo > 2.0) {
    const double m1 = std::floor(lo + (hi - lo) / 3.0);
    const double m2 = std::ceil(hi - (hi - lo) / 3.0);
    if (log_term_A(m1, s, lam, x, log_lam, log_x) <
        log_term_A(m2, s, lam, x, log_lam, log_x)) lo = m1;
    else hi = m2;
  }
  const double jp = lo;
  const double lp = log_term_A(jp, s, lam, x, log_lam, log_x);
  const double step = kWidth * std::sqrt(std::fmax(jp, 1.0)) + 5.0;
  double j0 = std::fmax(0.0, std::floor(jp - step));
  while (j0 > 0.0 && log_term_A(j0, s, lam, x, log_lam, log_x) > lp - 45.0)
    j0 = std::fmax(0.0, std::floor(j0 - step));
  return j0;
}

// A sum on a log scale: A = exp(L) * sum.
struct ScaledSum { double L; double sum; };

// A = sum_j Pois(j; lam) Q(K+j, x), upward from below its peak.
inline ScaledSum sum_A(const State& s, double lam, double x, double log_lam,
                       double log_x) {
  const double K = s.K;
  const double j0 = start_A(s, lam, x, log_lam, log_x);
  const double s0 = K + j0;
  const double lg0 = lgam_K(s, j0);
  double log_q0;
  gamma_log_pq_fast(s0, x, log_x, lg0, nullptr, &log_q0);
  const double log_s0 = (j0 == 0.0) ? s.log_K : std::log(s0);
  const double log_g0 = log_g(s0, x, log_x, lg0 + log_s0);
  const double top = std::fmax(log_q0, log_g0);
  double L = log_pois(j0, lam, log_lam) + top;
  if (!R_FINITE(L)) return {R_NegInf, 0.0};
  double tA = 1.0, tD = 1.0;
  if (log_q0 < top) tA = std::exp(log_q0 - top); else tD = std::exp(log_g0 - top);
  double sum = tA;
  double prev = tA;
  const double lx = lam * x;
  int n = 0;
  for (double j = j0;; j += 1.0, ++n) {
    if (n == kMaxTerms) return {R_NaN, R_NaN};  // never reached the bulk
    const double inv = 1.0 / (j + 1.0);
    tA = lam * inv * (tA + tD);
    tD *= lx * inv / (K + j + 1.0);
    sum += tA;
    if (j + 1.0 > lam && tA <= prev && tA < kTermEps * sum) break;
    prev = tA;
    if (sum > kRescale) {
      tA /= kRescale; tD /= kRescale; sum /= kRescale; prev /= kRescale;
      L += kLogRescale;
    }
  }
  return {L, sum};
}

// sum_j Pois(j; m_p) P(K+n0+j, z), downward from above its peak.  The peak
// lies below both the Poisson mode m_p and the peak of
// Pois(j; m_p) g(K+n0-1+j; z).
inline ScaledSum sum_P(const State& s, double n0, double m_p, double log_mp,
                       double z, double log_z) {
  const double o = s.K + n0;
  const double jg = 0.5 * (std::sqrt(o * o + 4.0 * m_p * z) - o);
  const double m = std::fmin(m_p, jg);
  const double j1 = std::ceil(m + kWidth * std::sqrt(m) + 5.0);
  // P(o+j1, z) and g(o-1+j1; z) share log Gamma(o+j1).
  const double lg1 = lgam_K(s, j1 + n0);
  double log_p1;
  gamma_log_pq_fast(o + j1, z, log_z, lg1, &log_p1, nullptr);
  const double log_e1 = log_g(o - 1.0 + j1, z, log_z, lg1);
  const double top = std::fmax(log_p1, log_e1);
  double L = log_pois(j1, m_p, log_mp) + top;
  if (!R_FINITE(L)) return {R_NegInf, 0.0};
  double tB = 1.0, tE = 1.0;
  if (log_p1 < top) tB = std::exp(log_p1 - top); else tE = std::exp(log_e1 - top);
  double sum = tB;
  double prev = tB;
  const double inv_m = 1.0 / m_p;
  const double inv_mz = 1.0 / (m_p * z);
  int n = 0;
  for (double j = j1; j > 0.0; j -= 1.0, ++n) {
    if (n == kMaxTerms) return {R_NaN, R_NaN};
    tB = j * inv_m * (tB + tE);
    tE *= j * (o - 1.0 + j) * inv_mz;
    sum += tB;
    if (j - 1.0 < m && tB <= prev && tB < kTermEps * sum) break;
    prev = tB;
    if (sum > kRescale) {
      tB /= kRescale; tE /= kRescale; sum /= kRescale; prev /= kRescale;
      L += kLogRescale;
    }
  }
  return {L, sum};
}

// The Gamma(K, a) first passage: log Q(K, au) (upper) or log P(K, au).
inline double log_gamma_tail(double u, const State& s, bool upper) {
  double lp, lq;
  gamma_log_pq_fast(s.K, s.a * u, s.log_a + std::log(u), s.lgK,
                    upper ? nullptr : &lp, upper ? &lq : nullptr);
  return upper ? lq : lp;
}

// log of (a/b)^K B, B = sum_j Pois(j; au) P(K+1+j, bu).
inline double log_B_term(const State& s, double x, double log_x, double lam,
                         double log_lam) {
  const ScaledSum B = sum_P(s, 1.0, x, log_x, lam, log_lam);
  if (ISNAN(B.L)) return R_NaN;
  return (B.L > R_NegInf)
    ? s.K * (s.log_a - s.log_b) + B.L + std::log(B.sum) : R_NegInf;
}

// log S(u) = log(A - B) from the sums, given log((a/b)^K B).
inline double log_surv_sums(const State& s, double x, double log_x, double lam,
                            double log_lam, double lB) {
  const ScaledSum A = sum_A(s, lam, x, log_lam, log_x);
  if (ISNAN(A.L) || ISNAN(lB)) return R_NaN;
  if (!(A.L > R_NegInf) || !(A.sum > 0.0)) return R_NegInf;
  // S = exp(A.L) (A.sum - exp(lB - A.L)).
  const double d = A.sum - ((lB > R_NegInf) ? std::exp(lB - A.L) : 0.0);
  if (!(d > 0.0)) return R_NegInf;
  return A.L + std::log(d);
}

// log S(u) of a proper accumulator (a >= b > 0), u finite and positive.
inline double log_surv_proper(double u, const State& s) {
  if (gamma_limit(u, s)) return log_gamma_tail(u, s, true);
  if (2.0 * u * std::sqrt(s.a * s.b) > hankel_threshold(s.K))
    return log_surv_hankel(u, s);
  const double x = s.a * u, lam = s.b * u;
  const double log_u = std::log(u);
  const double log_x = s.log_a + log_u, log_lam = s.log_b + log_u;
  return log_surv_sums(s, x, log_x, lam, log_lam,
                       log_B_term(s, x, log_x, lam, log_lam));
}

// log F(u) of a proper accumulator (a >= b > 0), u finite and positive.
// Below F = 1/2 the positive sum F = sum_j Pois(j; bu) P(K+j, au) + (a/b)^K B
// keeps relative accuracy however small F is; above it 1 - S is the more
// accurate.  In the Hankel range F is never small (the walk has spread
// sqrt(z) >= 2K about a nonnegative drift), so it is 1 - S there.
inline double log_cdf_proper(double u, const State& s) {
  if (gamma_limit(u, s)) return log_gamma_tail(u, s, false);
  if (2.0 * u * std::sqrt(s.a * s.b) > hankel_threshold(s.K))
    return log1m_exp(log_surv_hankel(u, s));
  const double x = s.a * u, lam = s.b * u;
  const double log_u = std::log(u);
  const double log_x = s.log_a + log_u, log_lam = s.log_b + log_u;
  const ScaledSum C = sum_P(s, 0.0, lam, log_lam, x, log_x);
  const double lB = log_B_term(s, x, log_x, lam, log_lam);
  if (ISNAN(C.L) || ISNAN(lB)) return R_NaN;
  const double lC = (C.L > R_NegInf) ? C.L + std::log(C.sum) : R_NegInf;
  const double lF = log_sum_exp(lC, lB);
  if (lF < -M_LN2) return lF;
  return log1m_exp(log_surv_sums(s, x, log_x, lam, log_lam, lB));
}

// log S(u) = log P(U > u), including the never-finish atom.
inline double log_surv(double u, const State& s) {
  if (!(u > 0.0)) return 0.0;
  if (s.b == 0.0) {
    if (u == R_PosInf) return R_NegInf;
    double lq;
    const double x = s.a * u;
    gamma_log_pq_fast(s.K, x, s.log_a + std::log(u), s.lgK, nullptr, &lq);
    return lq;
  }
  if (s.a >= s.b) {
    if (u == R_PosInf) return R_NegInf;
    return log_surv_proper(u, s);
  }
  const double l_atom = log1m_exp(s.log_h);
  if (u == R_PosInf) return l_atom;
  return log_sum_exp(l_atom, s.log_h + log_surv_proper(u, swapped(s)));
}

// log F(u) = log P(U <= u).
inline double log_cdf(double u, const State& s) {
  if (!(u > 0.0)) return R_NegInf;
  if (u == R_PosInf) return s.log_h;
  if (s.b == 0.0) {
    double lp;
    gamma_log_pq_fast(s.K, s.a * u, s.log_a + std::log(u), s.lgK, &lp, nullptr);
    return lp;
  }
  if (s.a < s.b) return s.log_h + log_cdf_proper(u, swapped(s));
  return log_cdf_proper(u, s);
}

// Scalar forms for callers that hold no State.
inline double log_pdf(double u, double K, double a, double b) {
  return log_pdf(u, make_state(K, a, b));
}
inline double log_surv(double u, double K, double a, double b) {
  return log_surv(u, make_state(K, a, b));
}
inline double log_cdf(double u, double K, double a, double b) {
  return log_cdf(u, make_state(K, a, b));
}

// u with log S(u) = lv for a proper accumulator (a >= b > 0), by safeguarded
// Newton in y = log u on g(y) = log S(e^y) - lv, which decreases from 0 to
// -Inf with g'(y) = -u f(u) / S(u).
inline double inv_surv_proper(double lv, const State& s) {
  const double K = s.K, a = s.a, b = s.b;
  double y;
  if (a > b) {
    // Lognormal match to the mean K / (a - b) and variance K (a + b) / (a - b)^3.
    const double d = a - b;
    const double s2 = std::log1p((a + b) / (K * d));
    y = std::log(K / d) - 0.5 * s2 +
        std::sqrt(s2) * R::qnorm(lv, 0.0, 1.0, 0, 1);
  } else {
    y = std::log(K / a);
  }
  auto g = [&](double yy) { return log_surv_proper(std::exp(yy), s) - lv; };
  double lo = y, hi = y;
  double glo = g(lo), ghi = glo;
  for (int i = 0; glo <= 0.0 && i < 200; ++i) { lo -= 1.0; glo = g(lo); }
  for (int i = 0; ghi >= 0.0 && i < 200; ++i) { hi += 1.0; ghi = g(hi); }
  if (!(glo > 0.0) || !(ghi < 0.0)) return R_NaN;
  y = (glo < -ghi) ? lo : hi;
  for (int it = 0; it < 100; ++it) {
    const double u = std::exp(y);
    const double ls = log_surv_proper(u, s);
    const double gy = ls - lv;
    if (gy > 0.0) lo = y; else hi = y;
    const double slope = -std::exp(y + log_pdf(u, s) - ls);
    double y_new = y - gy / slope;
    if (!(y_new > lo && y_new < hi) || !R_FINITE(y_new)) y_new = 0.5 * (lo + hi);
    if (std::fabs(y_new - y) < 1e-12 || hi - lo < 1e-12) return std::exp(y_new);
    y = y_new;
  }
  return std::exp(y);
}

// One decision time; +Inf is a never-finishing accumulator.  Uses R's RNG.
inline double draw_decision(double K, double a, double b) {
  const State s = make_state(K, a, b);
  if (!s.ok) return R_PosInf;
  if (b == 0.0) return R::rgamma(K, 1.0 / a);
  if (a < b) {
    // With probability h the accumulator finishes, with the law of the
    // proper accumulator whose rates are exchanged.
    if (!(std::log(R::unif_rand()) < s.log_h)) return R_PosInf;
    return inv_surv_proper(std::log(R::unif_rand()), swapped(s));
  }
  return inv_surv_proper(std::log(R::unif_rand()), s);
}

// Rates (a, b) from the sampled chart.  Chart 0: (v, rho) with
// a = (1 + rho) v, b = rho v.  Chart 1: (lambda_p, lambda_m) directly.
// Chart 2: (v, lambda0) with a = v + lambda0, b = lambda0.
inline void chart_rates(int chart, double c1, double c2, double* a, double* b) {
  if (chart == 0) { *a = (1.0 + c2) * c1; *b = c2 * c1; }
  else if (chart == 2) { *a = c1 + c2; *b = c2; }
  else { *a = c1; *b = c2; }
}

// ===========================================================================
// Gamma rate mixture (Math/gamma_mixtures.tex, sec. 2), ratio chart only.
//
// The net rate V ~ Gamma(beta, rate beta/v) is drawn once per trial, with
// beta = 1/cv^2, and lambda+ = (1 + rho) V, lambda- = rho V.  Each event of
// the superposed stream is positive with probability p+ = (1 + rho)/(1 + 2 rho)
// whatever V is, so the walk is a homogeneous counter of rate (1 + 2 rho) V
// with the random criterion K + 2J, J the negative events before absorption:
//
//   w_j = P(J = j) = K Gamma(K + 2j) / (j! Gamma(K + j + 1)) p+^(K+j) p-^j,
//   w_{j+1}/w_j = (K+2j)(K+2j+1) / ((j+1)(K+j+1)) p+ p-  ->  4 p+ p-.
//
// Integrating the Gamma(K + 2j) waiting time over V gives beta-prime laws
// with shape a = K + 2j, tail beta and scale r = beta/((1 + 2 rho) v), so with
// q = x/(r + x) and h_a = q^a (1-q)^beta / B(a, beta),
//
//   f = (1/x) sum_j w_j h_{K+2j},
//   F = sum_j w_j I_q(K+2j, beta),   S = sum_j w_j I_{1-q}(beta, K+2j).
//
// I_q(a, beta) - I_q(a+2, beta) = g_a = (h_a/a)(1 + q(a + beta)/(a + 1)), so
// F = sum_i g_{K+2i} W_i with W_i = w_0 + ... + w_i the cumulative weight:
// a forward sum of positive terms that needs no incomplete beta at all.  S is
// the forward sum of w_j S_j with S_{j+1} = S_j + g_{K+2j} from
// S_0 = I_{1-q}(beta, K), one library call.  Every series stops when a
// geometric bound on its remainder is below rounding: h_{a+2}/h_a =
// q^2 (a+beta)(a+beta+1)/(a(a+1)) falls with a, and the weight ratio stays
// below max(its current value, 4 p+ p-) (it falls to the limit while above it
// and rises to it from below).
//
// The F series converges at about q^2 per term and the S series at about the
// weight ratio, so F is summed while the beta-prime median of the first
// criterion is ahead, S otherwise, and each output is the complement of the
// other only where that leaves it at least one half.
// ===========================================================================
namespace gm {

constexpr double kEps = 1e-17;

struct GState {
  bool ok = false;
  double K = 0.0, beta = 0.0;
  double r = 0.0;        // beta-prime scale beta/((1 + 2 rho) v)
  double ppm = 0.0;      // p+ p-
  double w_lim = 0.0;    // 4 p+ p-
  double log_w0 = 0.0;   // K log p+
  double lbK = 0.0;      // log B(K, beta)
  double v = 0.0, rho = 0.0;
};

// w_0 = p+^K must stay representable: K log(1/p+) <= K log 2.
inline GState make_gstate(double K, double v, double rho, double cv) {
  GState g;
  if (!(K > 0.0 && v > 0.0 && rho >= 0.0 && cv > 0.0) || !R_FINITE(K) ||
      !R_FINITE(v) || !R_FINITE(rho) || !R_FINITE(cv)) return g;
  const double beta = 1.0 / (cv * cv);
  const double tot = 1.0 + 2.0 * rho;
  const double pm = rho / tot;
  g.K = K; g.beta = beta; g.v = v; g.rho = rho;
  g.r = beta / (tot * v);
  g.ppm = (1.0 - pm) * pm;
  g.w_lim = 4.0 * g.ppm;
  g.log_w0 = K * std::log1p(-pm);
  g.lbK = lbeta_fast(K, beta);
  g.ok = g.log_w0 > -700.0 && R_FINITE(g.lbK) && g.r > 0.0 && R_FINITE(g.r);
  return g;
}

// q, 1 - q and their logs at decision time x, each without cancellation.
struct QPt { double q, omq, log_q, log_omq; };
inline QPt qpoint(double x, double r) {
  const double rx = r + x, lrx = std::log(rx);
  return {x / rx, r / rx, std::log(x) - lrx, std::log(r) - lrx};
}

inline double log_h(double a, const QPt& p, double beta, double lb) {
  return a * p.log_q + beta * p.log_omq - lb;
}

// log f(x).
inline double log_pdf(double x, const GState& g) {
  if (!(x > 0.0) || !R_FINITE(x)) return R_NegInf;
  const QPt p = qpoint(x, g.r);
  const double K = g.K, b = g.beta, q2 = p.q * p.q;
  double L = g.log_w0 + log_h(K, p, b, g.lbK);
  double t = 1.0, sum = 1.0, a = K;
  for (int j = 0;; ++j, a += 2.0) {
    if (j == kMaxTerms) return R_NaN;
    const double wr = a * (a + 1.0) / ((j + 1.0) * (K + j + 1.0)) * g.ppm;
    const double hr = q2 * (a + b) * (a + b + 1.0) / (a * (a + 1.0));
    const double rb = std::max(wr, g.w_lim) * hr;
    if (rb < 1.0 && t * rb < kEps * sum * (1.0 - rb)) break;
    t *= wr * hr;
    sum += t;
    if (sum > kRescale) { t /= kRescale; sum /= kRescale; L += kLogRescale; }
  }
  return L + std::log(sum) - std::log(x);
}

// log F = log sum_i g_{K+2i} W_i, with h carried relative to exp(L).  The
// remainder after term i is at most sum_{l>i} g_l <= g_i gam/(1 - gam), gam
// bounding every later ratio g_{l+2}/g_l (W <= 1).
inline double log_F_series(const QPt& p, const GState& g) {
  const double K = g.K, b = g.beta, q = p.q, q2 = q * q;
  double L = log_h(K, p, b, g.lbK);
  if (!(L > R_NegInf)) return R_NegInf;
  double h = 1.0, w = std::exp(g.log_w0), W = w, sum = 0.0, a = K;
  for (int i = 0;; ++i, a += 2.0) {
    if (i == kMaxTerms) return R_NaN;
    const double br = 1.0 + q * (a + b) / (a + 1.0);
    const double gi = h / a * br;
    sum += gi * W;
    const double hr = q2 * (a + b) * (a + b + 1.0) / (a * (a + 1.0));
    const double gam = (b >= 1.0) ? hr : hr * (1.0 + q) / br;
    if (gam < 1.0 && gi * gam < kEps * sum * (1.0 - gam)) break;
    h *= hr;
    w *= a * (a + 1.0) / ((i + 1.0) * (K + i + 1.0)) * g.ppm;
    W += w;
    if (sum > kRescale || h > kRescale) {
      h /= kRescale; sum /= kRescale; L += kLogRescale;
    }
  }
  return L + std::log(sum);
}

// log S = log sum_j w_j S_j, with S_j and h carried relative to exp(L).
// The remainder after term j is at most Smax sum_{l>j} w_l, Smax bounding
// the later S_l (by one, or by S_{j+1} plus the geometric g tail).
inline double log_S_series(const QPt& p, const GState& g) {
  const double K = g.K, b = g.beta, q = p.q, q2 = q * q;
  const double lS0 = R::pbeta(p.omq, b, K, 1, 1);
  const double lh0 = log_h(K, p, b, g.lbK);
  if (ISNAN(lS0)) return R_NaN;
  double L = std::max(lS0, lh0);
  if (!(L > R_NegInf)) return R_NegInf;
  double S = std::exp(lS0 - L), h = std::exp(lh0 - L);
  double one = std::exp(-L);  // S <= 1 on this scale (may be Inf)
  double w = std::exp(g.log_w0), sum = 0.0, a = K;
  for (int j = 0;; ++j, a += 2.0) {
    if (j == kMaxTerms) return R_NaN;
    sum += w * S;
    const double br = 1.0 + q * (a + b) / (a + 1.0);
    const double gj = h / a * br;
    const double hr = q2 * (a + b) * (a + b + 1.0) / (a * (a + 1.0));
    const double wr = a * (a + 1.0) / ((j + 1.0) * (K + j + 1.0)) * g.ppm;
    S += gj;
    const double wb = std::max(wr, g.w_lim);
    if (wb < 1.0) {
      const double gam = (b >= 1.0) ? hr : hr * (1.0 + q) / br;
      double smax = one;
      if (gam < 1.0) smax = std::min(smax, S + gj * gam / (1.0 - gam));
      if (w * wb * smax < kEps * sum * (1.0 - wb)) break;
    }
    h *= hr;
    w *= wr;
    if (sum > kRescale || S > kRescale) {
      S /= kRescale; h /= kRescale; sum /= kRescale; one /= kRescale;
      L += kLogRescale;
    }
  }
  return L + std::log(sum);
}

// Whether the beta-prime median of the first criterion K is still ahead,
// where F <= I_q(K, beta) is the small side.
inline bool early(const QPt& p, const GState& g) {
  return p.q * (g.K + g.beta) <= g.K;
}

inline double log_surv(double x, const GState& g) {
  if (!(x > 0.0)) return 0.0;
  if (!R_FINITE(x)) return R_NegInf;
  const QPt p = qpoint(x, g.r);
  if (early(p, g)) {
    const double lF = log_F_series(p, g);
    if (lF < -M_LN2) return std::log1p(-std::exp(lF));
  }
  return log_S_series(p, g);
}

inline double log_cdf(double x, const GState& g) {
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return 0.0;
  const QPt p = qpoint(x, g.r);
  if (!early(p, g)) {
    const double lS = log_S_series(p, g);
    if (lS < -M_LN2) return std::log1p(-std::exp(lS));
  }
  return log_F_series(p, g);
}

// The trial rate, then the fixed-rate decision time at that rate.
inline double draw_decision(double K, double v, double rho, double cv) {
  const double V = R::rgamma(1.0 / (cv * cv), v * cv * cv);
  return opr::draw_decision(K, (1.0 + rho) * V, rho * V);
}

}  // namespace gm

// One accumulator row: the fixed-rate State, or the Gamma mixture when
// cv > 0 (ratio chart only).
struct Acc {
  bool mix = false;
  State s;
  gm::GState g;
  bool ok() const { return mix ? g.ok : s.ok; }
};

inline double log_pdf(double u, const Acc& c) {
  return c.mix ? gm::log_pdf(u, c.g) : log_pdf(u, c.s);
}
inline double log_surv(double u, const Acc& c) {
  return c.mix ? gm::log_surv(u, c.g) : log_surv(u, c.s);
}
inline double log_cdf(double u, const Acc& c) {
  return c.mix ? gm::log_cdf(u, c.g) : log_cdf(u, c.s);
}

// Two-slot memo of row states keyed on the sampled (K, c1, c2, cv), so that
// the accumulators of a trial, which alternate, and the rows of a condition
// reuse theirs.  NaN keys never match, so the first call always derives.
struct StateMemo {
  struct Slot {
    double K = R_NaN, c1 = R_NaN, c2 = R_NaN, cv = R_NaN;
    Acc s;
  };
  Slot slot[2];
  int cur = 0;
  int chart = 0;

  explicit StateMemo(int chart_) : chart(chart_) {}

  const Acc& get(double K, double c1, double c2, double cv = 0.0) {
    Slot& s0 = slot[cur];
    if (K == s0.K && c1 == s0.c1 && c2 == s0.c2 && cv == s0.cv) return s0.s;
    cur ^= 1;
    Slot& s1 = slot[cur];
    if (K == s1.K && c1 == s1.c1 && c2 == s1.c2 && cv == s1.cv) return s1.s;
    s1.K = K; s1.c1 = c1; s1.c2 = c2; s1.cv = cv;
    s1.s.mix = (cv != 0.0);
    if (s1.s.mix) {
      // Only the ratio chart separates the criterion from the rate.
      s1.s.g = (chart == 0) ? gm::make_gstate(K, c1, c2, cv) : gm::GState();
    } else {
      double a, b;
      chart_rates(chart, c1, c2, &a, &b);
      s1.s.s = make_state(K, a, b);
    }
    return s1.s;
  }
};

}  // namespace opr

// Race adapters (model_OPR.cpp).  Column layout: K = 0, then the two rate
// coordinates, t0 = 3, and cv = 4 under the ratio chart.
double dopr_scalar(double t, const double* par, void* ctx_);
double popr_scalar(double t, const double* par, void* ctx_);
void dopr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_);
void popr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_);
void opr_logS_at_t(double t, const double* const* cols,
                   int n_rows_total, int n_lR, int n_par,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* ctx_, double* logS_out);

#endif

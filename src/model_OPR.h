#ifndef EMC2_MODEL_OPR_H
#define EMC2_MODEL_OPR_H

// Opponent Poisson race (Math/Opponent_Poisson_Race.md).
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
// The survivor sums are used only for a >= b.  For a < b the completion
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
constexpr double kRescale = 1e250;       // running-scale renormalisation point
constexpr double kLogRescale = 575.6462732485114;  // log(1e250)
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
  for (double j = jp; n < kMaxTerms; j += 1.0, ++n) {
    term *= lx / ((j + 1.0) * (K + j + 1.0));
    sum += term;
    if (term < kTermEps * sum) break;
  }
  term = 1.0;
  n = 0;
  for (double j = jp; j > 0.0 && n < kMaxTerms; j -= 1.0, ++n) {
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
    const double d = std::sqrt(s.a) - std::sqrt(s.b);
    return s.log_K - log_u + 0.5 * s.K * (s.log_a - s.log_b) -
           d * d * u + log_bessel_i_scaled_hankel(s.K, z);
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
  for (double j = j0; n < kMaxTerms; j += 1.0, ++n) {
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

// sum_j Pois(j; x) P(K+1+j, lam), downward from above its peak.  The peak
// lies below both the Poisson mode x and the peak of
// Pois(j; x) g(K+1+j; lam).
inline ScaledSum sum_B(const State& s, double lam, double x, double log_lam,
                       double log_x) {
  const double K = s.K;
  const double K1 = K + 1.0;
  const double jg = 0.5 * (std::sqrt(K1 * K1 + 4.0 * x * lam) - K1);
  const double m = std::fmin(x, jg);
  const double j1 = std::ceil(m + kWidth * std::sqrt(m) + 5.0);
  // P(K+1+j1, lam) and g(K+j1; lam) share log Gamma(K+1+j1).
  const double lg1 = lgam_K(s, j1 + 1.0);
  double log_p1;
  gamma_log_pq_fast(K1 + j1, lam, log_lam, lg1, &log_p1, nullptr);
  const double log_e1 = log_g(K + j1, lam, log_lam, lg1);
  const double top = std::fmax(log_p1, log_e1);
  double L = log_pois(j1, x, log_x) + top;
  if (!R_FINITE(L)) return {R_NegInf, 0.0};
  double tB = 1.0, tE = 1.0;
  if (log_p1 < top) tB = std::exp(log_p1 - top); else tE = std::exp(log_e1 - top);
  double sum = tB;
  double prev = tB;
  const double inv_x = 1.0 / x;
  const double inv_xl = 1.0 / (x * lam);
  for (double j = j1; j > 0.0; j -= 1.0) {
    tB = j * inv_x * (tB + tE);
    tE *= j * (K + j) * inv_xl;
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

// log S(u) of a proper accumulator (a >= b > 0), u finite and positive.
inline double log_surv_proper(double u, const State& s) {
  const double x = s.a * u, lam = s.b * u;
  const double log_u = std::log(u);
  const double log_x = s.log_a + log_u, log_lam = s.log_b + log_u;
  const ScaledSum A = sum_A(s, lam, x, log_lam, log_x);
  if (!(A.L > R_NegInf) || !(A.sum > 0.0)) return R_NegInf;
  const ScaledSum B = sum_B(s, lam, x, log_lam, log_x);
  // S = A - B = exp(A.L) (A.sum - B.sum exp(K log(a/b) + B.L - A.L)).
  const double rB = (B.L > R_NegInf)
    ? B.sum * std::exp(s.K * (s.log_a - s.log_b) + B.L - A.L) : 0.0;
  const double d = A.sum - rB;
  if (!(d > 0.0)) return R_NegInf;
  return A.L + std::log(d);
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
  if (s.a < s.b) return s.log_h + log1m_exp(log_surv(u, swapped(s)));
  return log1m_exp(log_surv(u, s));
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

// Two-slot memo of row states keyed on the sampled (K, c1, c2), so that the
// accumulators of a trial, which alternate, and the rows of a condition reuse
// theirs.  NaN keys never match, so the first call always derives.
struct StateMemo {
  struct Slot {
    double K = R_NaN, c1 = R_NaN, c2 = R_NaN;
    State s;
  };
  Slot slot[2];
  int cur = 0;
  int chart = 0;

  explicit StateMemo(int chart_) : chart(chart_) {}

  const State& get(double K, double c1, double c2) {
    Slot& s0 = slot[cur];
    if (K == s0.K && c1 == s0.c1 && c2 == s0.c2) return s0.s;
    cur ^= 1;
    Slot& s1 = slot[cur];
    if (K == s1.K && c1 == s1.c1 && c2 == s1.c2) return s1.s;
    double a, b;
    chart_rates(chart, c1, c2, &a, &b);
    s1.K = K; s1.c1 = c1; s1.c2 = c2;
    s1.s = make_state(K, a, b);
    return s1.s;
  }
};

}  // namespace opr

// Race adapters (model_OPR.cpp).  Column layout: K = 0, then the two rate
// coordinates, t0 = 3.
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

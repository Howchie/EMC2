#ifndef EMC2_MODEL_FRQ_H
#define EMC2_MODEL_FRQ_H

// ============================================================================
// FRQ -- the Finite Reservoir Quorum process.
//
// GENERATIVE STORY.  One accumulator draws on N potential evidence units and
// responds once K distinct units have registered.  Each unit is available
// independently with probability p, and an available unit registers at a delay
// with mean rate lambda; clock = 0 is exponential registration and clock > 0
// integrates an independent Gamma-distributed unit rate.  An unavailable unit
// never registers at all.  The
// decision time is the K-th order statistic of the effective arrival times,
// which is +Inf whenever fewer than K units happen to be available.  That is
// the model's intrinsic omission mechanism, and it is NOT an additive
// contaminant: p also reshapes the finite-RT distribution.
//
// CLOSED FORM.  By time x exactly q(x) = p G(x) of the reservoir has
// registered, so the registered count is Binomial(N, q(x)) and a response has
// occurred iff that count is at least K.  Hence with A = K and R = N - K + 1
// (named `alpha` and `beta` here, because R is a reserved data column in EMC2
// and because U ~ Beta(alpha, beta) is the exact latent representation):
//
//   F(x) = I_{q(x)}(alpha, beta)
//   f(x) = p g(x) q(x)^{alpha-1} (1 - q(x))^{beta-1} / B(alpha, beta)
//   h    = F(Inf) = I_p(alpha, beta)          eventual completion probability
//
// For exponential registration G(x) = 1 - e^{-lambda x}, so its density is
// g(x) = G'(x) = lambda e^{-lambda x}.  With unit-rate CV c = clock > 0, the
// integrated registration CDF and density are
//   G_c(x) = 1 - (1 + c^2 lambda x)^(-1/c^2)
//   g_c(x) = lambda (1 + c^2 lambda x)^(-1/c^2 - 1).
// Nothing here needs a first-passage solve or a quadrature: the entire
// likelihood is four Rmath special functions.  The continuous relaxation is
// just alpha, beta > 0 rather than positive integers; every formula above is
// already valid there, because nothing above uses integrality -- only that
// I_q(alpha, beta) and the Beta density exist, which they do for any
// alpha, beta > 0. The fitted model keeps K = alpha >= 1 (R/model_FRQ.R),
// reached through an offset link K = 1 + exp(x): below one the density is
// singular at t0 and the likelihood is unbounded as t0 -> min(rt).  beta needs
// no such guard: 1 - p is solved from its own quantile, so beta < 1 (a rate
// SD above the mean rate) keeps every digit of the defect.
//
// COUNTER READING.  With phi(x) = q(x) / (1 - q(x)), the CDF is P(G_K <= W phi(x)) for independent
// G_K ~ Gamma(alpha) and W ~ Gamma(beta): a Poisson counter run on the
// operational clock phi, whose trial intensity W is Gamma distributed, and
// whose threshold is K = alpha events.  clock fixes the shape of the clock
// (1: constant rate, < 1 accelerating, > 1 decelerating) and phi'(0) = p lambda
// for every clock.
//
// PARAMETERISATION.  The fitted coordinates are (K, v, sv, h, t0, plo,
// clock), read in that counter language:
//   K   = alpha                 event threshold
//   v   = beta p lambda         mean initial input rate, E[W] phi'(0)
//   sv  = v / sqrt(beta)        between-trial SD of the counter's input rate
//   h   = F(Inf)                eventual completion probability
// The kernel derives beta = (v/sv)^2, p = I^{-1}_h(alpha, beta) and
// lambda = v / (beta p).  h = 1 gives p = 1 and a proper accumulator.
//
// sv = 0 is the exact beta -> Inf limit at fixed (K, v, h), a counter with no
// between-trial rate variability.  With h = 1 it is F(x) = P_K(v x), the
// regularised lower incomplete gamma (the classic Poisson counter; the clock
// drops out).  With h < 1 the total evidence on a trial is finite: E =
// P_K^{-1}(h) expected events, the clock runs at lambda = v / E, and
// F(x) = P_K(E G(x)).
//
// LOW STATE (plo).  On each trial each accumulator independently commits at
// the low threshold K_lo with probability plo, and at K otherwise.  Both
// thresholds run on the same counter (same W, clock and availability), so
//
//   F(x) = (1 - plo) I_{q(x)}(K, beta) + plo I_{q(x)}(K_lo, beta),
//
// and likewise for the density, the survivor and the cv = 0 limit with P_K.
// K_lo is fixed at FRQ_K_LO = 1, where the low component is elementary:
// I_q(1, beta) = 1 - (1 - q)^beta and P_1(X) = 1 - e^{-X}.  Its density at
// onset is plo v, so the low state sets the leading edge of the law.  The
// draw is per accumulator, so the race still factors over accumulators.
//
// The completion probability is a mixture too,
//
//   h = (1 - plo) I_p(K, beta) + plo I_p(K_lo, beta),
//
// so for plo > 0 and h < 1 the availability p (or E at cv = 0) is a
// one-dimensional monotone root, bracketed by the two pure quantiles and
// solved by safeguarded Newton (frq_mix_solve).  At h = 1, p = 1 exactly and
// nothing is solved.  plo = 0 is the base model, bit for bit.
//
// NUMERICS.  Everything is evaluated on the log scale.  The one real hazard is
// 1 - q(x) when q is close to one: forming it by subtraction loses every digit.
// It is available in closed form without cancellation,
//
//   1 - q(x) = (1 - p) + p e^{-lambda x},
//
// and that expression is what feeds both the density's (beta-1) log term and
// the survivor.  The survivor itself uses the reflection
// 1 - I_q(alpha, beta) = I_{1-q}(beta, alpha) rather than a subtraction, so a
// deep upper tail stays accurate all the way out to its defective limit
// S(Inf) = 1 - h.
//
// This header is pure inline Rmath and is safe to include from any number of
// translation units; the R-callable entry points live in model_FRQ.cpp, which
// is why they are not here.
// ============================================================================

#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <utility>
#include "inc_gamma.h"

using namespace Rcpp;

// Shapes below this are treated as degenerate.  The fitted model keeps
// alpha, beta >= 1 (R/model_FRQ.R bounds), so this floor is only a guard for
// the exported dfrq/pfrq/frq_rate entry points, which take shapes directly.
constexpr double FRQ_SHAPE_MIN = 1e-8;
// Below this rate CV (beta > 1e8) the beta branch differs from its
// beta -> Inf limit by O(cv^2), about 5e-9 in the CDF, while the incomplete
// beta and its quantile start losing digits to the huge shape; the exact
// Poisson-counter branch takes over.
constexpr double FRQ_CV_POISSON = 1e-4;

// Threshold of the low state.  Fixed rather than a fitted column; every
// evaluator takes it from FrqLow, so freeing it later only changes where the
// value comes from.
constexpr double FRQ_K_LO = 1.0;


// log(1 - e^y) for y <= 0.
inline double frq_log1mexp(double y) {
  return (y > -M_LN2) ? std::log(-std::expm1(y)) : std::log1p(-std::exp(y));
}

inline double frq_logaddexp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  const double m = std::max(a, b);
  return m + std::log1p(std::exp(-std::abs(a - b)));
}

// log Q_{a,b}(u), the log beta quantile.  For a < 1 the quantile underflows
// long before u does, so the small-x side uses the leading term of
// I_x(a, b) = x^a / (a B(a, b)) (1 + a(1-b) x/(a+1) + ...), exact to double
// precision once that series correction is below rounding.
inline double frq_log_qbeta(double u, double a, double b, double lbeta_ab) {
  if (!(u > 0.0)) return R_NegInf;
  if (!(u < 1.0)) return 0.0;
  const double ly0 = (std::log(u) + std::log(a) + lbeta_ab) / a;
  if (ly0 < std::log(1e-17 * a * (a + 1.0) / (std::abs(1.0 - b) + 1.0)))
    return ly0;
  const double x = R::qbeta(u, a, b, 1, 0);
  return (x > 0.0) ? std::log(x) : ly0;
}

// ---------------------------------------------------------------------------
// The low state.  An inactive FrqLow (plo = 0) must cost nothing: every call
// site branches on `active` rather than mixing with a zero weight.
// ---------------------------------------------------------------------------

inline double frq_log1pexp(double x) {
  return (x > 0.0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

struct FrqLow {
  bool ok = true;           // false only for a rejected plo
  bool active = false;      // plo > 0
  double plo = 0.0;
  double log_plo = R_NegInf;
  double log_1m_plo = 0.0;
  double K_lo = FRQ_K_LO;
  double K_lo_minus_1 = FRQ_K_LO - 1.0;
  double lbeta_lo = 0.0;    // lbeta(K_lo, beta), beta branch
  double lgamma_lo = 0.0;   // lgamma(K_lo), cv = 0 branch
};

inline FrqLow frq_low_make(double plo, double K_lo = FRQ_K_LO) {
  FrqLow L;
  if (ISNAN(plo) || !(plo >= 0.0) || !(plo <= 1.0) ||
      !(K_lo >= FRQ_SHAPE_MIN) || !R_FINITE(K_lo)) { L.ok = false; return L; }
  L.K_lo = K_lo;
  L.K_lo_minus_1 = K_lo - 1.0;
  if (plo == 0.0) return L;                   // base model; `active` stays false
  L.active = true;
  L.plo = plo;
  L.log_plo = std::log(plo);
  L.log_1m_plo = std::log1p(-plo);            // -Inf at plo = 1
  L.lgamma_lo = lgamma_fast(K_lo);
  return L;
}

// log[(1 - plo) e^main + plo e^low].
inline double frq_low_mix(double main, double low, const FrqLow& L) {
  return frq_logaddexp(L.log_1m_plo + main, L.log_plo + low);
}

// log P(K_lo, X) or log Q(K_lo, X); closed form at K_lo = 1.
inline double frq_low_gamma_tail(double X, double log_X, const FrqLow& L,
                                 bool lower) {
  if (L.K_lo == 1.0) return lower ? frq_log1mexp(-X) : -X;
  double out;
  if (lower) gamma_log_pq_fast(L.K_lo, X, log_X, L.lgamma_lo, &out, nullptr);
  else gamma_log_pq_fast(L.K_lo, X, log_X, L.lgamma_lo, nullptr, &out);
  return out;
}

// ---------------------------------------------------------------------------
// Geometric criterion excess (omega).  The criterion is
//   J = Y_{K-1} + (1 + omega) Y_1,   Y_{K-1} ~ Gamma(K - 1), Y_1 ~ Exp(1),
// so the last criterion unit is stretched by 1 + omega; the low state's
// criterion is (1 + omega) Y_1 (K_lo = 1).  With r = omega / (1 + omega),
// s = 1 - r and the counter's cumulative intensity M,
//   P(J <= M)   = P_K(M) - r^{1-K} e^{-s M} P_K(r M),
//   density J   = s r^{1-K} e^{-s M} P_{K-1}(r M),
// and averaged over W ~ Gamma(beta) at q = phi / (1 + phi),
//   F = I_q(K, beta) - T,   S = I_{1-q}(beta, K) + T,
//   T = r^{1-K} [(1 - q) / (1 - r q)]^beta I_{rq}(K, beta),
//   dF/dq = s r^{1-K} beta (1 - q)^{beta-1} (1 - r q)^{-beta-1}
//           I_{rq}(K - 1, beta + 1).
// The survivor and density are sums of positive terms; the CDF loses about
// log10(1 + omega) digits to its difference.  omega = 0 is the base model and
// never reaches these paths.
// ---------------------------------------------------------------------------

// Below this K - 1 the (K - 1)-shape factor is one (P_0 = I(0, .) = 1).
constexpr double FRQ_OMEGA_UNIT_K = 1e-10;

struct FrqOmega {
  bool ok = true;           // false only for a rejected omega
  bool active = false;      // omega > 0
  double omega = 0.0;
  double log_r = R_NegInf;  // log(omega / (1 + omega))
  double log_s = 0.0;       // log(1 / (1 + omega))
  double s = 1.0;
  double log1p_omega = 0.0;
  bool main_unit = true;    // K - 1 below FRQ_OMEGA_UNIT_K
  double lb_main1 = 0.0;    // lbeta(K - 1, beta + 1), or lgamma(K - 1) at cv = 0
  bool lo_unit = true;      // K_lo - 1 below FRQ_OMEGA_UNIT_K
  double lb_lo1 = 0.0;
};

inline FrqOmega frq_omega_make(double omega) {
  FrqOmega w;
  if (ISNAN(omega) || !(omega >= 0.0) || !R_FINITE(omega)) { w.ok = false; return w; }
  if (omega == 0.0) return w;
  w.active = true;
  w.omega = omega;
  w.log1p_omega = std::log1p(omega);
  w.log_s = -w.log1p_omega;
  w.log_r = std::log(omega) - w.log1p_omega;
  w.s = 1.0 / (1.0 + omega);
  return w;
}

// log F, log S and the log density with respect to the clock coordinate (q in
// the beta branch, the cumulative intensity X at cv = 0).
struct FrqOmComp {
  double lF = R_NaN;
  double lS = R_NaN;
  double lf = R_NaN;
};
enum : int { FRQ_OM_F = 1, FRQ_OM_S = 2, FRQ_OM_FD = 4 };

inline FrqOmComp frq_omega_mix_beta(double log_q, double log_omq, double K,
                                    double beta, double lb, const FrqLow& L,
                                    const FrqOmega& w, int want);
inline FrqOmComp frq_omega_mix_pois(double X, double log_X, double K,
                                    double lgK, const FrqLow& L,
                                    const FrqOmega& w, int want);

// ---------------------------------------------------------------------------
// Compiled state, split by what each half actually depends on.
//
// The evidence-level availability p is the expensive half: it comes from the
// fitted completion probability h by an incomplete-beta INVERSION, an iterative
// solve costing roughly fourteen density evaluations, plus a root solve when
// plo > 0.  It is a function of (K, cv, h, plo) ONLY, with cv = sv / v the rate CV; clock does not enter
// it.  The rate half is one log and a reciprocal.
//
// Keeping them on separate keys matters whenever the rate varies faster than
// the shapes.  Because cv = sv / v, a per-trial v moves the key only when
// sv > 0 and h < 1: at h = 1 there is no quantile, and at sv = 0 the key is
// cv = 0 whatever v is.  FrqAvail is what the shared cache
// stores; FrqPars is what the evaluators read, assembled by pasting a cached
// FrqAvail together with the row's rate.
//
// In both halves `ok` is false for any combination that does not define a
// proper accumulator; every evaluator then returns log 0 for the density and
// lets the caller floor the trial, which is what the rest of the race
// machinery expects from an invalid row.
// ---------------------------------------------------------------------------

// The (K, cv, h, plo) half -- everything the beta quantile touches, with
// cv = sv / v.  `poisson` marks the cv -> 0 limit, where beta is infinite and the kernel is
// the incomplete gamma P_K(E G(x)) (h < 1) or P_K(v x) (h = 1).  E is the
// expected total number of events available on a trial: beta p / (1 - p) for
// cv > 0, the root of the mixed completion probability in the limit, and Inf
// whenever h = 1.
struct FrqAvail {
  double alpha = 0.0;
  double beta = 0.0;
  bool poisson = false;
  double E = 0.0;
  double log_E = 0.0;
  double lgamma_K = 0.0;
  double p = 0.0;
  double log_p = 0.0;
  double om_p = 0.0;
  double log_om_p = R_NegInf;
  double h = 0.0;
  double lbeta_val = 0.0;
  double log_p_minus_lbeta = 0.0;
  double alpha_minus_1 = 0.0;
  double beta_minus_1 = 0.0;
  FrqLow lo;
  FrqOmega om;
  bool ok = false;
};

struct FrqPars {
  double alpha = 0.0;
  double beta = 0.0;
  bool poisson = false;
  double E = 0.0;
  double log_E = 0.0;
  double lgamma_K = 0.0;
  double v = 0.0;
  double log_v = 0.0;
  double p = 0.0;
  double log_p = 0.0;
  double om_p = 0.0;
  double log_om_p = R_NegInf;
  double h = 0.0;
  double lambda = 0.0;
  double log_lambda = 0.0;
  double clock = 0.0;
  double cv2 = 0.0;
  double cv2_lambda = 0.0;
  double inv_cv2 = 0.0;
  double lbeta_val = 0.0;
  double log_p_minus_lbeta = 0.0;
  double alpha_minus_1 = 0.0;
  double beta_minus_1 = 0.0;
  FrqLow lo;
  FrqOmega om;
  bool ok = false;
};

inline double frq_log_ibeta(double x, double omx, double log_x, double log_omx,
                            double a, double b, double lbeta_val);

// The beta quantile of h as the pair (log p, log(1 - p)).  1 - p comes from
// the complementary quantile (1 - B ~ Beta(b, a)), since 1 - qbeta(h, a, b)
// loses every digit once p is near one, which is the normal case for b < 1.
// p itself is re-solved only when it is the small side.  `closed` takes the
// elementary inverse at a = 1, which the low state uses; the main threshold
// keeps the quantile path.
inline bool frq_beta_quantile_logs(double h, double a, double b, double lb,
                                   double* log_p, double* log_om_p,
                                   bool closed = false) {
  if (closed && a == 1.0) {             // I_p(1, b) = 1 - (1 - p)^b
    *log_om_p = std::log1p(-h) / b;
    *log_p = frq_log1mexp(*log_om_p);
  } else {
    *log_om_p = frq_log_qbeta(1.0 - h, b, a, lb);
    if (ISNAN(*log_om_p) || !(*log_om_p < 0.0)) return false;
    *log_p = (*log_om_p > -M_LN2) ? frq_log_qbeta(h, a, b, lb)
                                  : frq_log1mexp(*log_om_p);
  }
  return !ISNAN(*log_p) && *log_p > R_NegInf && *log_p <= 0.0 &&
         *log_om_p < 0.0;
}

// The mixed completion probability m(t) and its slope at a solver coordinate t:
// log m, log(1 - m) and log dm/dt.  Only the side the residual uses is formed.
struct FrqMixEval {
  double lm;
  double l1m;
  double lD;
};

// Root of m(t) = h for an increasing m, bracketed by [t_a, t_b] and started at
// t_start.  The residual is log m - log h on the lower half and
// log(1 - h) - log(1 - m) on the upper half, so a completion probability within
// rounding of either end keeps its digits.  Newton steps that leave the
// bracket, or fail to shrink the residual, are replaced by bisection.
template <class Eval>
inline bool frq_mix_solve(double h, double t_a, double t_b, double t_start,
                          Eval ev, double* t_out) {
  const bool lower = (h <= 0.5);
  const double lh = std::log(h), l1h = std::log1p(-h);
  double a = std::min(t_a, t_b), b = std::max(t_a, t_b);
  double t = (t_start > a && t_start < b) ? t_start : 0.5 * (a + b);
  if (!(b > a)) { *t_out = a; return R_FINITE(a); }
  for (int it = 0; it < 100; ++it) {
    const FrqMixEval e = ev(t, lower);
    const double r = lower ? e.lm - lh : l1h - e.l1m;
    if (ISNAN(r)) return false;
    if (std::abs(r) <= 1e-14) break;
    if (r < 0.0) a = t; else b = t;
    const double dr = std::exp(e.lD - (lower ? e.lm : e.l1m));
    double tn = t - r / dr;
    if (!(tn > a && tn < b)) tn = 0.5 * (a + b);
    if (std::abs(tn - t) <= 1e-15 * (1.0 + std::abs(t))) { t = tn; break; }
    t = tn;
  }
  *t_out = t;
  return R_FINITE(t);
}

inline void frq_derive_avail_into(FrqAvail& s, double K, double cv, double h,
                                  double plo = 0.0, double omega = 0.0) {
  s = FrqAvail();
  if (ISNAN(K) || ISNAN(cv) || ISNAN(h)) return;
  if (!(K >= FRQ_SHAPE_MIN) || !R_FINITE(K)) return;
  if (!(cv >= 0.0) || !R_FINITE(cv)) return;
  if (!(h > 0.0) || !(h <= 1.0) || !R_FINITE(h)) return;
  s.lo = frq_low_make(plo);
  if (!s.lo.ok) return;
  s.om = frq_omega_make(omega);
  if (!s.om.ok) return;
  // The excess criterion splits off one unit, Y_{K-1}, so it needs K >= 1.
  if (s.om.active && !(K >= 1.0)) return;
  const FrqLow& L = s.lo;
  FrqOmega& W = s.om;
  W.main_unit = !(K - 1.0 >= FRQ_OMEGA_UNIT_K);
  W.lo_unit = !(L.K_lo - 1.0 >= FRQ_OMEGA_UNIT_K);
  s.alpha = K;
  s.alpha_minus_1 = K - 1.0;
  s.h = h;

  if (cv < FRQ_CV_POISSON) {
    // beta -> Inf at fixed (K, v, h).  The availability p -> 0 while the
    // expected total evidence E = beta p / (1 - p) tends to the root of
    // h = (1 - plo) P_K(E) + plo P_{K_lo}(E).
    s.poisson = true;
    s.beta = R_PosInf;
    s.lgamma_K = lgamma_fast(K);
    if (W.active) {
      if (!W.main_unit) W.lb_main1 = lgamma_fast(K - 1.0);
      if (!W.lo_unit) W.lb_lo1 = lgamma_fast(L.K_lo - 1.0);
    }
    if (h == 1.0) {
      s.E = R_PosInf; s.log_E = R_PosInf;
      s.p = 1.0; s.log_p = 0.0; s.om_p = 0.0; s.log_om_p = R_NegInf;
    } else if (W.active) {
      // The excess criterion lies between Y_K and (1 + omega) Y_K, so the root
      // lies in [E_0, (1 + omega) E_0] for each component's omega = 0 root E_0.
      double ta = R_PosInf, tb = R_NegInf;
      if (L.plo < 1.0) {
        const double E0 = R::qgamma(h, K, 1.0, 1, 0);
        if (!(E0 > 0.0) || !R_FINITE(E0)) return;
        ta = std::min(ta, std::log(E0)); tb = std::max(tb, std::log(E0));
      }
      if (L.active) {
        const double E0 = (L.K_lo == 1.0) ? -std::log1p(-h)
                                          : R::qgamma(h, L.K_lo, 1.0, 1, 0);
        if (!(E0 > 0.0) || !R_FINITE(E0)) return;
        ta = std::min(ta, std::log(E0)); tb = std::max(tb, std::log(E0));
      }
      tb += W.log1p_omega;
      const double lgK = s.lgamma_K;
      auto ev = [&](double t, bool lower) {
        const FrqOmComp c = frq_omega_mix_pois(
          std::exp(t), t, K, lgK, L, W,
          (lower ? FRQ_OM_F : FRQ_OM_S) | FRQ_OM_FD);
        FrqMixEval e;
        e.lm = c.lF; e.l1m = c.lS; e.lD = c.lf + t;
        return e;
      };
      double t;
      if (!frq_mix_solve(h, ta, tb, 0.5 * (ta + tb), ev, &t)) return;
      s.E = std::exp(t); s.log_E = t;
      if (!(s.E > 0.0) || !R_FINITE(s.E)) return;
      s.p = 0.0; s.log_p = R_NegInf; s.om_p = 1.0; s.log_om_p = 0.0;
    } else {
      double E = (L.plo == 1.0) ? R_NaN : R::qgamma(h, K, 1.0, 1, 0);
      if (L.active) {
        const double E_lo = (L.K_lo == 1.0) ? -std::log1p(-h)
                                            : R::qgamma(h, L.K_lo, 1.0, 1, 0);
        if (L.plo == 1.0) {
          E = E_lo;
        } else {
          if (!(E > 0.0) || !R_FINITE(E) || !(E_lo > 0.0)) return;
          const double lgK = s.lgamma_K;
          auto ev = [&](double t, bool lower) {
            const double X = std::exp(t);
            FrqMixEval e;
            double lmain, llow;
            if (lower) {
              gamma_log_pq_fast(K, X, t, lgK, &lmain, nullptr);
              llow = frq_low_gamma_tail(X, t, L, true);
              e.lm = frq_low_mix(lmain, llow, L);
            } else {
              gamma_log_pq_fast(K, X, t, lgK, nullptr, &lmain);
              llow = frq_low_gamma_tail(X, t, L, false);
              e.l1m = frq_low_mix(lmain, llow, L);
            }
            e.lD = frq_low_mix(K * t - lgK, L.K_lo * t - L.lgamma_lo, L) - X;
            return e;
          };
          const double ta = std::log(E), tb = std::log(E_lo);
          double t;
          if (!frq_mix_solve(h, ta, tb, ta + L.plo * (tb - ta), ev, &t)) return;
          E = std::exp(t);
        }
      }
      if (!(E > 0.0) || !R_FINITE(E)) return;
      s.E = E; s.log_E = std::log(E);
      s.p = 0.0; s.log_p = R_NegInf; s.om_p = 1.0; s.log_om_p = 0.0;
    }
    s.ok = true;
    return;
  }

  const double beta = 1.0 / (cv * cv);
  if (!(beta >= FRQ_SHAPE_MIN) || !R_FINITE(beta)) return;
  if (L.active) s.lo.lbeta_lo = lbeta_fast(L.K_lo, beta);
  // h = 1 maps exactly to p = 1 for any shape and plo; skipping the quantile
  // there matters because the proper boundary is the common case.  Both sides
  // are held as logs, so a 1 - p that underflows still enters log(1 - q(x))
  // exactly.
  const double lb = lbeta_fast(K, beta);
  if (W.active) {
    if (!W.main_unit) W.lb_main1 = lbeta_fast(K - 1.0, beta + 1.0);
    if (!W.lo_unit) W.lb_lo1 = lbeta_fast(L.K_lo - 1.0, beta + 1.0);
  }
  double log_p = 0.0, log_om_p = R_NegInf;
  if (h < 1.0 && W.active) {
    // Solved in t = logit p.  The excess criterion lies between Y_K and
    // (1 + omega) Y_K, which brackets each component's root between its
    // omega = 0 logit t_0 and t_0 + log(1 + omega).
    double ta = R_PosInf, tb = R_NegInf;
    double lp0, lomp0;
    if (L.plo < 1.0) {
      if (!frq_beta_quantile_logs(h, K, beta, lb, &lp0, &lomp0)) return;
      ta = std::min(ta, lp0 - lomp0); tb = std::max(tb, lp0 - lomp0);
    }
    if (L.active) {
      if (!frq_beta_quantile_logs(h, L.K_lo, beta, s.lo.lbeta_lo, &lp0, &lomp0,
                                  true)) return;
      ta = std::min(ta, lp0 - lomp0); tb = std::max(tb, lp0 - lomp0);
    }
    tb += W.log1p_omega;
    auto ev = [&](double t, bool lower) {
      const double lp = -frq_log1pexp(-t), lomp = -frq_log1pexp(t);
      const FrqOmComp c = frq_omega_mix_beta(
        lp, lomp, K, beta, lb, s.lo, W,
        (lower ? FRQ_OM_F : FRQ_OM_S) | FRQ_OM_FD);
      FrqMixEval e;
      e.lm = c.lF; e.l1m = c.lS; e.lD = c.lf + lp + lomp;
      return e;
    };
    double t;
    if (!frq_mix_solve(h, ta, tb, 0.5 * (ta + tb), ev, &t)) return;
    log_p = -frq_log1pexp(-t);
    log_om_p = -frq_log1pexp(t);
    if (ISNAN(log_p) || !(log_p > R_NegInf) || log_p > 0.0) return;
  } else if (h < 1.0) {
    if (L.plo == 1.0) {
      if (!frq_beta_quantile_logs(h, L.K_lo, beta, s.lo.lbeta_lo, &log_p,
                                  &log_om_p, true)) return;
    } else {
      if (!frq_beta_quantile_logs(h, K, beta, lb, &log_p, &log_om_p)) return;
      if (L.active) {
        // Solved in t = logit p, held as the pair (log p, log(1 - p)).
        double lp_lo, lomp_lo;
        if (!frq_beta_quantile_logs(h, L.K_lo, beta, s.lo.lbeta_lo, &lp_lo,
                                    &lomp_lo, true)) return;
        const double lb_lo = s.lo.lbeta_lo;
        auto ev = [&](double t, bool lower) {
          const double lp = -frq_log1pexp(-t), lomp = -frq_log1pexp(t);
          const double p = std::exp(lp), omp = std::exp(lomp);
          FrqMixEval e;
          if (lower) {
            e.lm = frq_low_mix(frq_log_ibeta(p, omp, lp, lomp, K, beta, lb),
                               frq_log_ibeta(p, omp, lp, lomp, L.K_lo, beta,
                                             lb_lo), L);
          } else {
            e.l1m = frq_low_mix(frq_log_ibeta(omp, p, lomp, lp, beta, K, lb),
                                frq_log_ibeta(omp, p, lomp, lp, beta, L.K_lo,
                                              lb_lo), L);
          }
          e.lD = frq_low_mix(K * lp - lb, L.K_lo * lp - lb_lo, L) + beta * lomp;
          return e;
        };
        const double ta = log_p - log_om_p, tb = lp_lo - lomp_lo;
        double t;
        if (!frq_mix_solve(h, ta, tb, ta + L.plo * (tb - ta), ev, &t)) return;
        log_p = -frq_log1pexp(-t);
        log_om_p = -frq_log1pexp(t);
      }
    }
    if (ISNAN(log_p) || !(log_p > R_NegInf) || log_p > 0.0) return;
  }

  s.beta = beta;
  s.p = std::exp(log_p);
  s.log_p = log_p;
  s.om_p = std::exp(log_om_p);
  s.log_om_p = log_om_p;
  s.log_E = (h == 1.0) ? R_PosInf : std::log(beta) + log_p - log_om_p;
  s.E = std::exp(s.log_E);
  s.lbeta_val = lb;
  s.log_p_minus_lbeta = s.log_p - s.lbeta_val;
  s.beta_minus_1 = beta - 1.0;
  s.ok = true;
}

inline FrqAvail frq_derive_avail(double K, double cv, double h,
                                 double plo = 0.0, double omega = 0.0) {
  FrqAvail s;
  frq_derive_avail_into(s, K, cv, h, plo, omega);
  return s;
}

// Paste the availability half into the full state.  Runs only when one of
// (K, cv, h, plo) actually changes.
inline void frq_copy_avail(FrqPars& s, const FrqAvail& av) {
  s.alpha = av.alpha;
  s.beta = av.beta;
  s.poisson = av.poisson;
  s.E = av.E;
  s.log_E = av.log_E;
  s.lgamma_K = av.lgamma_K;
  s.p = av.p;
  s.log_p = av.log_p;
  s.om_p = av.om_p;
  s.log_om_p = av.log_om_p;
  s.h = av.h;
  s.lbeta_val = av.lbeta_val;
  s.log_p_minus_lbeta = av.log_p_minus_lbeta;
  s.alpha_minus_1 = av.alpha_minus_1;
  s.beta_minus_1 = av.beta_minus_1;
  s.lo = av.lo;
  s.om = av.om;
}

// Fill the rate half from the mean initial input rate v.  The registration
// rate follows from the availability half: lambda = v / (beta p), or v / E in
// the cv = 0 limit with h < 1.  With cv = 0 and h = 1 the clock drops out of
// the likelihood and lambda is unused.  On rejection only `ok` is meaningful.
inline void frq_set_rate(FrqPars& s, bool avail_ok, double v, double clock) {
  s.ok = false;
  if (!avail_ok) return;
  if (ISNAN(v) || ISNAN(clock)) return;
  if (!(v > 0.0) || !R_FINITE(v)) return;
  if (!(clock >= 0.0) || !R_FINITE(clock)) return;
  const double cv2 = clock * clock;
  if (!R_FINITE(cv2)) return;
  const bool clock_free = s.poisson && s.h == 1.0;
  const double lambda = clock_free ? 0.0
    : (s.poisson ? v / s.E : v / (s.beta * s.p));
  if (!clock_free && (!(lambda > 0.0) || !R_FINITE(lambda))) return;
  s.v = v;
  s.log_v = std::log(v);
  s.lambda = lambda;
  s.log_lambda = clock_free ? R_NegInf : std::log(lambda);
  s.clock = clock;
  s.cv2 = cv2;
  s.cv2_lambda = cv2 * lambda;
  s.inv_cv2 = (cv2 > 0.0) ? (1.0 / cv2) : 0.0;
  s.ok = true;
}

// The fitted rate variability is the SD sv of the trial input rate; the
// beta-quantile half works in its CV, sv / v, which is 1/sqrt(beta).
inline double frq_rate_cv(double sv, double v) {
  return (sv == 0.0) ? 0.0 : sv / v;
}

// One-shot derive, for callers that hold no memo.
inline FrqPars frq_derive(double K, double v, double sv, double h,
                          double plo = 0.0, double clock = 1.0,
                          double omega = 0.0) {
  FrqPars s;
  const FrqAvail av = frq_derive_avail(K, frq_rate_cv(sv, v), h, plo, omega);
  frq_copy_avail(s, av);
  frq_set_rate(s, av.ok, v, clock);
  return s;
}

// Shared exact-key cache for the availability half.  The likelihood clears
// this once per particle, so entries are shared across density, survivor,
// scalar and truncation callbacks without accumulating proposals from earlier
// particles.  The key is exactly the beta quantile's arguments, so a rate that
// varies per trial does not multiply the number of inversions.
struct FrqCacheKey {
  std::uint64_t bits[5];

  bool operator==(const FrqCacheKey& other) const {
    for (int i = 0; i < 5; ++i) {
      if (bits[i] != other.bits[i]) return false;
    }
    return true;
  }
};

struct FrqCacheKeyHash {
  std::size_t operator()(const FrqCacheKey& key) const {
    std::size_t h = 0;
    for (int i = 0; i < 5; ++i) {
      h ^= std::hash<std::uint64_t>{}(key.bits[i]) +
           static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
           (h << 6) + (h >> 2);
    }
    return h;
  }
};

inline std::uint64_t frq_cache_bits(double x) {
  // Match the scalar memo's == semantics: +0 and -0 are the same key.
  if (x == 0.0) return 0;
  std::uint64_t bits = 0;
  std::memcpy(&bits, &x, sizeof(bits));
  return bits;
}

struct FrqCache {
  std::unordered_map<FrqCacheKey, FrqAvail, FrqCacheKeyHash> entries;

  FrqCache() { entries.reserve(256); }

  void clear() { entries.clear(); }

  // std::unordered_map is node-based, so a reference handed out here survives
  // later insertions; only clear() invalidates it, and that runs between
  // particles, outside any evaluator call.
  const FrqAvail& get(double K, double cv, double h, double plo,
                      double omega = 0.0) {
    const FrqCacheKey key{{frq_cache_bits(K), frq_cache_bits(cv),
                           frq_cache_bits(h), frq_cache_bits(plo),
                           frq_cache_bits(omega)}};
    const auto found = entries.find(key);
    if (found != entries.end()) return found->second;
    const auto inserted = entries.emplace(
      key, frq_derive_avail(K, cv, h, plo, omega));
    return inserted.first->second;
  }
};

// q(x) = p G_c(x) together with its cancellation-free complement
// 1 - q(x) = (1 - p) + p S_c(x).  x == Inf gives the saturated pair
// (p, 1 - p); the CDF and survivor return h and 1 - h there directly, because
// 1 - p can underflow when beta is small.
struct FrqState {
  double q;
  double omq;
  double log_q;
  double log_omq;
  double log_g;
};

inline double frq_registration_inv(double ratio, const FrqPars& s) {
  if (s.cv2 == 0.0)
    return -std::log1p(-ratio) / s.lambda;
  return std::expm1(-s.cv2 * std::log1p(-ratio)) / s.cv2_lambda;
}

// The same inverse given the log unit survivor log(1 - G(x)), for draws whose
// G(x) is too close to one to be held directly.
inline double frq_registration_inv_logsurv(double log_sr, const FrqPars& s) {
  if (s.cv2 == 0.0)
    return -log_sr / s.lambda;
  return std::expm1(-s.cv2 * log_sr) / s.cv2_lambda;
}

inline FrqState frq_state(double x, const FrqPars& s) {
  FrqState z;
  if (!R_FINITE(x)) {           // x == +Inf: fully saturated reservoir
    z.q = s.p;
    z.omq = s.om_p;
    z.log_q = s.log_p;
    z.log_omq = s.log_om_p;
    z.log_g = R_NegInf;
    return z;
  }
  double log_sr;
  if (s.cv2 == 0.0) {
    log_sr = -s.lambda * x;
    z.log_g = s.log_lambda + log_sr;
  } else {
    const double lx = std::log1p(s.cv2_lambda * x);
    log_sr = -lx * s.inv_cv2;
    // Split the exponent to avoid multiplying a potentially huge 1/cv^2 by a
    // small log1p term; lx/cv^2 has a finite exponential-limit value.
    z.log_g = s.log_lambda - lx * s.inv_cv2 - lx;
  }
  // 1 - q = (1 - p) + p S(x) is kept in logs: for beta < 1 it underflows
  // while (1 - q)^beta is still of order one.
  z.log_q = s.log_p + frq_log1mexp(log_sr);
  z.log_omq = frq_logaddexp(s.log_om_p, s.log_p + log_sr);
  z.q = std::exp(z.log_q);
  z.omq = std::exp(z.log_omq);
  return z;
}

// ---------------------------------------------------------------------------
// Regularized incomplete beta log-evaluator: log I_x(a, b), with omx = 1 - x.
// Evaluates Gauss's continued fraction in the convergent regime
// x < (a + 1)/(a + b + 2) and uses the symmetry I_x(a, b) = 1 - I_{1-x}(b, a)
// elsewhere. Special-cases a == 1 and b == 1 algebraically.
//
// The fraction f = 1 + d_1/(1 + d_2/(1 + ...)) is evaluated by the forward
// recurrence on its convergents, A_n = A_{n-1} + d_n A_{n-2} (B likewise),
// which needs no division per step; the pair is rescaled when B leaves
// [1e-100, 1e100], and convergence is |A_n B_{n-1} - A_{n-1} B_n| small
// against A_n B_{n-1}.
// ---------------------------------------------------------------------------
inline double frq_ibeta_cf_log(double x, double a, double b, double log_x,
                               double log_omx, double lbeta_val) {
  const double EPS = 1e-15;
  const int MAXIT = 200;
  const double apb = a + b;

  double Am2 = 1.0, Bm2 = 0.0, Am1 = 1.0, Bm1 = 1.0;
  const double d1 = -apb * x / (a + 1.0);
  double A = Am1 + d1 * Am2, B = Bm1 + d1 * Bm2;
  Am2 = Am1; Bm2 = Bm1; Am1 = A; Bm1 = B;
  for (int m = 1; m <= MAXIT; ++m) {
    const double m2 = 2.0 * m;
    const double d_even = (m * (b - m) * x) / ((a + m2 - 1.0) * (a + m2));
    A = Am1 + d_even * Am2; B = Bm1 + d_even * Bm2;
    Am2 = Am1; Bm2 = Bm1; Am1 = A; Bm1 = B;
    const double d_odd = -((a + m) * (apb + m) * x) / ((a + m2) * (a + m2 + 1.0));
    A = Am1 + d_odd * Am2; B = Bm1 + d_odd * Bm2;
    const double diff = A * Bm1 - Am1 * B;
    Am2 = Am1; Bm2 = Bm1; Am1 = A; Bm1 = B;
    if (std::abs(diff) <= EPS * std::abs(A * Bm2)) break;
    const double scale = std::abs(B);
    if (scale > 1e100 || scale < 1e-100) {
      const double r = 1.0 / scale;
      Am1 *= r; Bm1 *= r; Am2 *= r; Bm2 *= r;
    }
  }

  const double log_pref = a * log_x + b * log_omx - std::log(a) - lbeta_val;
  return log_pref - std::log(A / B);
}

// The prefactor takes log x and log(1 - x) directly: either side may underflow
// while the other rounds to one.
inline double frq_log_ibeta(double x, double omx, double log_x, double log_omx,
                            double a, double b, double lbeta_val) {
  if (log_x == R_NegInf) return R_NegInf;
  if (log_omx == R_NegInf) return 0.0;

  // Exact algebraic special cases
  if (b == 1.0) return a * log_x;
  if (a == 1.0) return frq_log1mexp(b * log_omx);   // I_x(1, b) = 1 - omx^b

  const double thresh = (a + 1.0) / (a + b + 2.0);
  if (x < thresh) {
    return frq_ibeta_cf_log(x, a, b, log_x, log_omx, lbeta_val);
  } else {
    // Symmetry: I_x(a, b) = 1 - I_{omx}(b, a)
    const double log_w = frq_ibeta_cf_log(omx, b, a, log_omx, log_x, lbeta_val);
    return (log_w < 0.0) ? frq_log1mexp(log_w) : R_NegInf;
  }
}

// ---------------------------------------------------------------------------
// Excess-criterion components (see FrqOmega).  One criterion shape Kc per
// call; the mixers combine the main and low thresholds with plo.
// ---------------------------------------------------------------------------

inline FrqOmComp frq_omega_beta_comp(double log_q, double log_omq, double Kc,
                                     double beta, double lb_c, bool unit_c,
                                     double lb_c1, const FrqOmega& w,
                                     int want) {
  FrqOmComp out;
  const double q = std::exp(log_q), omq = std::exp(log_omq);
  const double log_rq = w.log_r + log_q, rq = std::exp(log_rq);
  // 1 - r q = s + r (1 - q), with no cancellation.
  const double log_omrq = frq_logaddexp(w.log_s, w.log_r + log_omq);
  const double omrq = std::exp(log_omrq);
  const double lead = (1.0 - Kc) * w.log_r;
  if (want & (FRQ_OM_F | FRQ_OM_S)) {
    const double T = lead + beta * (log_omq - log_omrq) +
      frq_log_ibeta(rq, omrq, log_rq, log_omrq, Kc, beta, lb_c);
    out.lS = frq_logaddexp(frq_log_ibeta(omq, q, log_omq, log_q, beta, Kc, lb_c),
                           T);
    if (want & FRQ_OM_F) {
      const double lI = frq_log_ibeta(q, omq, log_q, log_omq, Kc, beta, lb_c);
      out.lF = (T < lI) ? lI + frq_log1mexp(T - lI) : frq_log1mexp(out.lS);
    }
  }
  if (want & FRQ_OM_FD) {
    double lf = w.log_s + lead + std::log(beta) - (beta + 1.0) * log_omrq;
    if (beta != 1.0) lf += (beta - 1.0) * log_omq;
    if (!unit_c)
      lf += frq_log_ibeta(rq, omrq, log_rq, log_omrq, Kc - 1.0, beta + 1.0,
                          lb_c1);
    out.lf = lf;
  }
  return out;
}

// log P(a, z) and log Q(a, z); closed form at a = 1.
inline void frq_omega_pgamma(double a, double z, double log_z, double lga,
                             double* lp, double* lq) {
  if (a == 1.0) {
    if (lp) *lp = frq_log1mexp(-z);
    if (lq) *lq = -z;
    return;
  }
  gamma_log_pq_fast(a, z, log_z, lga, lp, lq);
}

inline FrqOmComp frq_omega_pois_comp(double X, double log_X, double Kc,
                                     double lg_c, bool unit_c, double lg_c1,
                                     const FrqOmega& w, int want) {
  FrqOmComp out;
  const double lead = (1.0 - Kc) * w.log_r;
  if (!(X > 0.0)) {
    out.lF = R_NegInf;
    out.lS = 0.0;
    out.lf = unit_c ? w.log_s + lead : R_NegInf;
    return out;
  }
  const double log_rX = w.log_r + log_X, rX = std::exp(log_rX);
  if (want & (FRQ_OM_F | FRQ_OM_S)) {
    double lPr, lP, lQ;
    frq_omega_pgamma(Kc, rX, log_rX, lg_c, &lPr, nullptr);
    frq_omega_pgamma(Kc, X, log_X, lg_c, &lP, &lQ);
    const double T = lead - w.s * X + lPr;
    out.lS = frq_logaddexp(lQ, T);
    if (want & FRQ_OM_F)
      out.lF = (T < lP) ? lP + frq_log1mexp(T - lP) : frq_log1mexp(out.lS);
  }
  if (want & FRQ_OM_FD) {
    double lf = w.log_s + lead - w.s * X;
    if (!unit_c) {
      double lP1;
      frq_omega_pgamma(Kc - 1.0, rX, log_rX, lg_c1, &lP1, nullptr);
      lf += lP1;
    }
    out.lf = lf;
  }
  return out;
}

inline FrqOmComp frq_omega_mix(const FrqOmComp& m, const FrqOmComp& l,
                               const FrqLow& L, int want) {
  FrqOmComp out;
  if (want & FRQ_OM_F) out.lF = frq_low_mix(m.lF, l.lF, L);
  if (want & FRQ_OM_S) out.lS = frq_low_mix(m.lS, l.lS, L);
  if (want & FRQ_OM_FD) out.lf = frq_low_mix(m.lf, l.lf, L);
  return out;
}

inline FrqOmComp frq_omega_mix_beta(double log_q, double log_omq, double K,
                                    double beta, double lb, const FrqLow& L,
                                    const FrqOmega& w, int want) {
  const FrqOmComp m = frq_omega_beta_comp(log_q, log_omq, K, beta, lb,
                                          w.main_unit, w.lb_main1, w, want);
  if (!L.active) return m;
  const FrqOmComp l = frq_omega_beta_comp(log_q, log_omq, L.K_lo, beta,
                                          L.lbeta_lo, w.lo_unit, w.lb_lo1, w,
                                          want);
  return frq_omega_mix(m, l, L, want);
}

inline FrqOmComp frq_omega_mix_pois(double X, double log_X, double K,
                                    double lgK, const FrqLow& L,
                                    const FrqOmega& w, int want) {
  const FrqOmComp m = frq_omega_pois_comp(X, log_X, K, lgK, w.main_unit,
                                          w.lb_main1, w, want);
  if (!L.active) return m;
  const FrqOmComp l = frq_omega_pois_comp(X, log_X, L.K_lo, L.lgamma_lo,
                                          w.lo_unit, w.lb_lo1, w, want);
  return frq_omega_mix(m, l, L, want);
}

// ---------------------------------------------------------------------------
// Log density / CDF / survivor of the DECISION time x = rt - t0.  Callers are
// responsible for subtracting t0; x <= 0 is handled here so that the race
// kernels can pass a negative x without branching.
// ---------------------------------------------------------------------------

// The cv = 0 limit evaluates through its gamma argument X(x): v x when h = 1,
// E G(x) when h < 1 (so X(Inf) = E).
struct FrqPoissonArg {
  double X;
  double log_dX;
};

inline FrqPoissonArg frq_poisson_arg(double x, const FrqPars& s) {
  FrqPoissonArg a;
  if (s.h == 1.0) {
    a.X = s.v * x;
    a.log_dX = s.log_v;
    return a;
  }
  if (!R_FINITE(x)) {
    a.X = s.E;
    a.log_dX = R_NegInf;
    return a;
  }
  double G, log_g;
  if (s.cv2 == 0.0) {
    G = -std::expm1(-s.lambda * x);
    log_g = s.log_lambda - s.lambda * x;
  } else {
    const double lx = std::log1p(s.cv2_lambda * x);
    G = -std::expm1(-lx * s.inv_cv2);
    log_g = s.log_lambda - lx * s.inv_cv2 - lx;
  }
  a.X = s.E * G;
  a.log_dX = s.log_E + log_g;
  return a;
}

inline double frq_poisson_log_pdf(double x, const FrqPars& s) {
  const FrqPoissonArg a = frq_poisson_arg(x, s);
  if (!(a.X > 0.0) || !R_FINITE(a.X)) return R_NegInf;
  const double log_X = std::log(a.X);
  double lp;
  if (s.lo.active) {
    const double lmain = s.alpha_minus_1 * log_X - s.lgamma_K;
    const double llow = (s.lo.K_lo_minus_1 != 0.0)
      ? s.lo.K_lo_minus_1 * log_X - s.lo.lgamma_lo : -s.lo.lgamma_lo;
    lp = frq_low_mix(lmain, llow, s.lo) - a.X + a.log_dX;
  } else {
    lp = s.alpha_minus_1 * log_X - a.X - s.lgamma_K + a.log_dX;
  }
  return ISNAN(lp) ? R_NegInf : lp;
}

// Only the requested tail is evaluated.
inline double frq_poisson_log_tail(double x, const FrqPars& s, bool lower) {
  const FrqPoissonArg a = frq_poisson_arg(x, s);
  const double log_X = std::log(a.X);
  double out;
  if (lower) gamma_log_pq_fast(s.alpha, a.X, log_X, s.lgamma_K, &out, nullptr);
  else gamma_log_pq_fast(s.alpha, a.X, log_X, s.lgamma_K, nullptr, &out);
  if (s.lo.active)
    out = frq_low_mix(out, frq_low_gamma_tail(a.X, log_X, s.lo, lower), s.lo);
  return out;
}

inline double frq_log_pdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return R_NegInf;   // no atom at infinity in the density
  if (s.om.active) {
    double lp;
    if (s.poisson) {
      const FrqPoissonArg a = frq_poisson_arg(x, s);
      if (!(a.X > 0.0) || !R_FINITE(a.X)) return R_NegInf;
      lp = frq_omega_mix_pois(a.X, std::log(a.X), s.alpha, s.lgamma_K, s.lo,
                              s.om, FRQ_OM_FD).lf + a.log_dX;
    } else {
      const FrqState z = frq_state(x, s);
      lp = frq_omega_mix_beta(z.log_q, z.log_omq, s.alpha, s.beta, s.lbeta_val,
                              s.lo, s.om, FRQ_OM_FD).lf + s.log_p + z.log_g;
    }
    return ISNAN(lp) ? R_NegInf : lp;
  }
  if (s.poisson) return frq_poisson_log_pdf(x, s);
  const FrqState z = frq_state(x, s);
  double lp;
  if (s.lo.active) {
    // Both components share p g(x) (1 - q)^{beta - 1}; only q^{k - 1} / B(k, beta)
    // differs.
    double lmain = -s.lbeta_val, llow = -s.lo.lbeta_lo;
    if (s.alpha_minus_1 != 0.0) lmain += s.alpha_minus_1 * z.log_q;
    if (s.lo.K_lo_minus_1 != 0.0) llow += s.lo.K_lo_minus_1 * z.log_q;
    lp = s.log_p + z.log_g + frq_low_mix(lmain, llow, s.lo);
  } else {
    lp = s.log_p_minus_lbeta + z.log_g;
    if (s.alpha_minus_1 != 0.0) lp += s.alpha_minus_1 * z.log_q;
  }
  if (s.beta_minus_1 != 0.0) lp += s.beta_minus_1 * z.log_omq;
  return ISNAN(lp) ? R_NegInf : lp;
}

inline double frq_log_cdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return std::log(s.h);   // F(Inf) = h by construction
  if (s.om.active) {
    double out;
    if (s.poisson) {
      const FrqPoissonArg a = frq_poisson_arg(x, s);
      out = frq_omega_mix_pois(a.X, std::log(a.X), s.alpha, s.lgamma_K, s.lo,
                               s.om, FRQ_OM_F).lF;
    } else {
      const FrqState z = frq_state(x, s);
      out = frq_omega_mix_beta(z.log_q, z.log_omq, s.alpha, s.beta,
                               s.lbeta_val, s.lo, s.om, FRQ_OM_F).lF;
    }
    return ISNAN(out) ? R_NegInf : out;
  }
  if (s.poisson) {
    const double out = frq_poisson_log_tail(x, s, true);
    return ISNAN(out) ? R_NegInf : out;
  }
  const FrqState z = frq_state(x, s);
  double out = frq_log_ibeta(z.q, z.omq, z.log_q, z.log_omq, s.alpha, s.beta,
                             s.lbeta_val);
  if (s.lo.active)
    out = frq_low_mix(out, frq_log_ibeta(z.q, z.omq, z.log_q, z.log_omq,
                                         s.lo.K_lo, s.beta, s.lo.lbeta_lo),
                      s.lo);
  return ISNAN(out) ? R_NegInf : out;
}

// log S(x) = log(1 - F(x)), via I_{1-q}(beta, alpha) so that the defective
// limit log(1 - h) is reached without ever subtracting two near-equal numbers.
// An unstarted accumulator (x <= 0) has survivor one, i.e. log S = 0.  The low
// component reflects the same way, I_{1-q}(beta, K_lo) = (1 - q)^beta at
// K_lo = 1.
inline double frq_log_surv_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NaN;             // caller distinguishes bad from certain
  if (!(x > 0.0)) return 0.0;
  if (!R_FINITE(x)) return std::log1p(-s.h);
  if (s.om.active) {
    double out;
    if (s.poisson) {
      const FrqPoissonArg a = frq_poisson_arg(x, s);
      out = frq_omega_mix_pois(a.X, std::log(a.X), s.alpha, s.lgamma_K, s.lo,
                               s.om, FRQ_OM_S).lS;
    } else {
      const FrqState z = frq_state(x, s);
      out = frq_omega_mix_beta(z.log_q, z.log_omq, s.alpha, s.beta,
                               s.lbeta_val, s.lo, s.om, FRQ_OM_S).lS;
    }
    return ISNAN(out) ? R_NaN : out;
  }
  if (s.poisson) {
    const double out = frq_poisson_log_tail(x, s, false);
    return ISNAN(out) ? R_NaN : out;
  }
  const FrqState z = frq_state(x, s);
  double out = frq_log_ibeta(z.omq, z.q, z.log_omq, z.log_q, s.beta, s.alpha,
                             s.lbeta_val);
  if (s.lo.active)
    out = frq_low_mix(out, frq_log_ibeta(z.omq, z.q, z.log_omq, z.log_q,
                                         s.beta, s.lo.K_lo, s.lo.lbeta_lo),
                      s.lo);
  return ISNAN(out) ? R_NaN : out;
}

// ---------------------------------------------------------------------------
// Natural-scale wrappers.  Used by the truncation normalisers and by the
// scalar pdf1/cdf1 pointers, both of which want a probability in [0, 1] and
// tolerate saturation rather than wanting an -Inf.
// ---------------------------------------------------------------------------

inline double frq_pdf_natural_dt(double x, const FrqPars& s) {
  const double lp = frq_log_pdf_dt(x, s);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

inline double frq_cdf_natural_dt(double x, const FrqPars& s) {
  const double lp = frq_log_cdf_dt(x, s);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

// ---------------------------------------------------------------------------
// Accelerating clock (the accel chart).  The FRQ state runs on the operational
// time A(u) = (e^{a u} - 1) / a (A(u) = u at a = 0), so its input intensity
// is multiplied by e^{a u}:
//
//   F(u) = F_s(A(u)),   f(u) = f_s(A(u)) e^{a u},   S(u) = S_s(A(u)).
//
// A(Inf) = Inf, so h, p and the defective mass are those of the state itself.
// a = 0 is a branch, not a limit, so it reproduces the state bit for bit.
// ---------------------------------------------------------------------------

inline bool frq_accel_ok(double a) { return R_FINITE(a) && a >= 0.0; }

inline double frq_accel_warp(double u, double a) {
  if (a == 0.0 || !R_FINITE(u)) return u;
  return std::expm1(a * u) / a;
}

// Physical time of an operational time x; +Inf stays an omission.
inline double frq_accel_unwarp(double x, double a) {
  if (a == 0.0 || !R_FINITE(x)) return x;
  return std::log1p(a * x) / a;
}

inline double frq_accel_log_pdf_dt(double u, const FrqPars& s, double a) {
  if (!frq_accel_ok(a)) return R_NegInf;
  if (a == 0.0) return frq_log_pdf_dt(u, s);
  if (!(u > 0.0) || !R_FINITE(u)) return R_NegInf;
  const double lp = frq_log_pdf_dt(frq_accel_warp(u, a), s);
  return (lp > R_NegInf) ? lp + a * u : R_NegInf;
}

inline double frq_accel_log_cdf_dt(double u, const FrqPars& s, double a) {
  if (!frq_accel_ok(a)) return R_NegInf;
  return frq_log_cdf_dt(frq_accel_warp(u, a), s);
}

inline double frq_accel_log_surv_dt(double u, const FrqPars& s, double a) {
  if (!frq_accel_ok(a)) return R_NaN;
  return frq_log_surv_dt(frq_accel_warp(u, a), s);
}

// Front-end memo for every call site.  It keeps the two halves on separate
// keys, so a row that moves only the rate -- a trend or a covariate on v,
// against constant shapes -- redoes one log rather than one beta quantile.
// The likelihood adapters hand it the per-particle FrqCache; the exported
// entry points run without one and memoise locally.
// At h = 1 there is no beta quantile (p = 1), so a valid proper state depends
// on K only through alpha and lgamma(K) (sv = 0) or lbeta(K, beta) (sv > 0),
// and a K that changes per row (a trend) updates those instead of
// re-deriving.  The rate half is K-free.
// The excess criterion's K - 1 terms are not refreshed here, so an active
// omega always re-derives.
inline bool frq_proper_k_update(const FrqAvail& a, double K) {
  return a.ok && a.h == 1.0 && !a.om.active && K >= FRQ_SHAPE_MIN &&
         R_FINITE(K);
}

template <class S>
inline void frq_set_proper_k(S& s, double K, double lg) {
  s.alpha = K;
  s.alpha_minus_1 = K - 1.0;
  if (s.poisson) {
    s.lgamma_K = lg;
  } else {
    s.lbeta_val = lg;
    s.log_p_minus_lbeta = s.log_p - lg;
  }
}

// lgamma(K) or lbeta(K, beta), whichever the proper state carries.
inline double frq_proper_k_term(const FrqAvail& a, double K) {
  return a.poisson ? lgamma_fast(K) : lbeta_fast(K, a.beta);
}

struct FrqMemo {
  FrqCache* shared = nullptr;
  // Two availability states, keyed on (K, sv/v, h, plo), so rows that
  // alternate between two keys (the accumulators of one trial) find theirs
  // without a re-derive.  The memo owns its copies rather than pointing into
  // the shared cache, so its lifetime is independent of when that is cleared.
  struct Slot {
    double K = R_NaN, cv = R_NaN, h = R_NaN, plo = R_NaN, omega = R_NaN;
    FrqAvail avail;
    bool matches(double k, double c, double hh, double d, double o) const {
      return k == K && c == cv && hh == h && d == plo && o == omega;
    }
  };
  Slot slot[2];
  int cur = 0;
  double v = R_NaN, clock = R_NaN;
  FrqPars value;

  explicit FrqMemo(FrqCache* shared_cache = nullptr) : shared(shared_cache) {}

  const FrqPars& get(double k, double vv, double sv, double hh,
                     double d = 0.0, double cu = 1.0, double o = 0.0) {
    return get_cv(k, vv, frq_rate_cv(sv, vv), hh, d, cu, o);
  }

  // The same state from the rate CV c = sv / v directly (the cv chart).
  const FrqPars& get_cv(double k, double vv, double c, double hh,
                        double d = 0.0, double cu = 1.0, double o = 0.0) {
    // NaN never compares equal, so a NaN parameter always re-derives, and the
    // NaN-initialised keys make the first call always miss.
    Slot& s0 = slot[cur];
    if (!s0.matches(k, c, hh, d, o)) {
      Slot& s1 = slot[cur ^ 1];
      if (c == s0.cv && hh == s0.h && d == s0.plo && o == s0.omega &&
          frq_proper_k_update(s0.avail, k)) {
        s0.K = k;
        const double lg = frq_proper_k_term(s0.avail, k);
        frq_set_proper_k(s0.avail, k, lg);
        frq_set_proper_k(value, k, lg);
        if (vv == v && cu == clock) return value;
      } else if (s1.matches(k, c, hh, d, o)) {
        cur ^= 1;
        frq_copy_avail(value, s1.avail);
      } else {
        cur ^= 1;
        if (c == s1.cv && hh == s1.h && d == s1.plo && o == s1.omega &&
            frq_proper_k_update(s1.avail, k)) {
          s1.K = k;
          frq_set_proper_k(s1.avail, k, frq_proper_k_term(s1.avail, k));
        } else {
          s1.K = k; s1.cv = c; s1.h = hh; s1.plo = d; s1.omega = o;
          // h = 1 needs no beta quantile, so the derive is cheaper than a
          // shared-cache lookup and insert.
          if (shared != nullptr && hh != 1.0)
            s1.avail = shared->get(k, c, hh, d, o);
          else frq_derive_avail_into(s1.avail, k, c, hh, d, o);
        }
        frq_copy_avail(value, s1.avail);
      }
    } else if (vv == v && cu == clock) {
      return value;
    }
    v = vv; clock = cu;
    frq_set_rate(value, slot[cur].avail.ok, vv, cu);
    return value;
  }
};

// Race adapters used by the particle likelihood dispatch.  Definitions live
// in model_FRQ.cpp so this header remains safe for model_rng.cpp.  The
// `_cv` set reads the rate CV (column "cv") where the default reads sv; the
// `_accel` set reads the rate CV and the clock acceleration (column "accel")
// in the clock's position.  The layouts are otherwise identical.
double dfrq_scalar(double t, const double* par, void* /*ctx_*/);
double pfrq_scalar(double t, const double* par, void* /*ctx_*/);
void dfrq_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_);
void pfrq_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_);
void frq_logS_at_t(double t, const double* const* cols,
                   int /*n_rows_total*/, int n_lR, int /*n_par*/,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* /*ctx_*/, double* logS_out);
double dfrq_cv_scalar(double t, const double* par, void* /*ctx_*/);
double pfrq_cv_scalar(double t, const double* par, void* /*ctx_*/);
void dfrq_cv_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_);
void pfrq_cv_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_);
void frq_cv_logS_at_t(double t, const double* const* cols,
                      int /*n_rows_total*/, int n_lR, int /*n_par*/,
                      const int* trunc_mask, int n_unique_trials,
                      const int* isok_all, void* /*ctx_*/, double* logS_out);
double dfrq_accel_scalar(double t, const double* par, void* /*ctx_*/);
double pfrq_accel_scalar(double t, const double* par, void* /*ctx_*/);
void dfrq_accel_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_);
void pfrq_accel_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_);
void frq_accel_logS_at_t(double t, const double* const* cols,
                         int /*n_rows_total*/, int n_lR, int /*n_par*/,
                         const int* trunc_mask, int n_unique_trials,
                         const int* isok_all, void* /*ctx_*/, double* logS_out);

#endif  // EMC2_MODEL_FRQ_H

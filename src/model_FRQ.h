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
// PARAMETERISATION.  The fitted coordinates are (K, v, sv, h, t0, delta,
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
// THRESHOLD VARIABILITY (delta).  The quorum a trial actually demands is not
// fixed.  Let U be the latent quorum percentile, so that the baseline model is
// P(U <= u) = I_u(alpha, beta), and shift its LOG-ODDS by a trialwise state
// eps ~ Uniform(-delta, delta):
//
//   P(U <= u | eps) = logit^{-1}[ logit I_u(alpha, beta) + eps ].
//
// Positive eps makes small effective criteria more likely (a liberal trial),
// negative eps makes large ones more likely (a conservative trial), so delta is
// symmetric between-trial variability in caution on an odds scale, with
// SD(eps) = delta/sqrt(3).  Marginalising eps is a closed-form integral,
// because the antiderivative of logit^{-1} is the softplus:
//
//   H(z) = [log(1 + e^{x+delta}) - log(1 + e^{x-delta})] / (2 delta),  x = logit z
//
// and the model becomes F(x) = H(I_{q(x)}(alpha, beta)).  Three properties make
// this cost almost nothing here:
//
//   * H is a SCALAR POST-MAP on the existing CDF.  It does not touch the
//     incomplete beta, so the a == 1 / b == 1 fast branches below and the
//     hoisted lbeta all stay live.
//   * H is symmetric, H(1 - z) = 1 - H(z), so the survivor needs no new code
//     path: it is the same H evaluated at logit(1 - z) = -x.  This composes
//     with the reflection identity already used for I.
//   * H is analytically invertible,
//
//       H^{-1}(y) = sinh(delta y) / [sinh(delta y) + sinh(delta (1 - y))],
//
//     so deriving p from h remains stable without a second time-scale anchor.
//     Solve H(z) = y for e^x directly: (1+e^{x+d})/(1+e^{x-d}) = e^{2dy} gives
//     e^x = (e^{2dy} - 1)/(e^d - e^{2dy-d}) = sinh(dy)/sinh(d(1-y)).
//
// The derivative needed by the density is likewise elementary,
//
//   H'(z) = (sinh(delta)/delta) / (1 + 2 z (1-z) (cosh(delta) - 1)),
//
// which is > 1 at both endpoints and < 1 at z = 1/2: mass is pulled off middle
// quorum percentiles toward both extremes.  delta = 0 recovers the base model
// exactly (H'(z) = 1) and is the default; the delta > 0 density does have to
// evaluate the incomplete beta, which the base density does not.
//
// A caveat that is NOT visible from the formulas: eps is marginalised inside
// each accumulator's own CDF, so this is per-accumulator criterion noise (like
// LBA's sv), not a single caution state shared by a trial's whole race.  A
// shared eps would not factor through f_i prod_j S_j and would need quadrature.
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

// Beyond this point sinh(delta) and the log-odds inversion lose finite-double
// guarantees.  The fitted model is much more conservative (delta <= 6), but
// the exported entry points also accept direct values and must reject rather
// than manufacture NaNs for extreme inputs.
constexpr double FRQ_DELTA_MAX = 700.0;

// xlogy(a, y) = a log y, with the a == 0 branch returning 0 rather than the
// NaN of 0 * (-Inf).  alpha == 1 or beta == 1 hits this on every call.
inline double frq_xlogy(double a, double y) {
  if (a == 0.0) return 0.0;
  if (y <= 0.0) return (a > 0.0) ? R_NegInf : R_PosInf;
  return a * std::log(y);
}

// ---------------------------------------------------------------------------
// The threshold-variability generator H_delta and its inverse/derivative.  All
// three collapse to the identity at delta = 0, which is the default, so an
// inactive FrqH must cost nothing: every call site branches on `active` rather
// than evaluating an identity map.
// ---------------------------------------------------------------------------

// log sinh(u) for u > 0, without overflowing at large u.  Below 1/2 the direct
// evaluation is already exact (sinh is faithful down to the denormals); above
// it, factor out e^u so the exponential never exceeds one.
inline double frq_logsinh(double u) {
  if (!(u > 0.0)) return R_NegInf;
  if (u < 0.5) return std::log(std::sinh(u));
  return u + std::log1p(-std::exp(-2.0 * u)) - 0.6931471805599453;
}

struct FrqH {
  bool ok = true;         // false only for a rejected delta
  bool active = false;    // delta > 0
  double delta = 0.0;
  double two_delta = 0.0;

  double two_sinh = 0.0;    // 2 sinh(delta)
  double exp_mdelta = 1.0;  // e^{-delta}
  double log_sinhc = 0.0;   // log(sinh(delta)/delta) = log H'(0)
  double c2 = 0.0;          // 2(cosh(delta) - 1), formed as 4 sinh^2(delta/2)
};
inline FrqH frq_h_make(double delta) {
  FrqH H;
  if (ISNAN(delta) || delta < 0.0 || !R_FINITE(delta) ||
      delta > FRQ_DELTA_MAX) { H.ok = false; return H; }
  if (delta == 0.0) return H;                 // identity; `active` stays false
  H.active = true;
  H.delta = delta;
  H.two_delta = 2.0 * delta;
  H.two_sinh = 2.0 * std::sinh(delta);
  H.exp_mdelta = std::exp(-delta);
  // 2(cosh d - 1) suffers total cancellation for small d; the half-angle
  // identity 2(cosh d - 1) = 4 sinh^2(d/2) does not.
  const double sh = std::sinh(0.5 * delta);
  H.c2 = 4.0 * sh * sh;
  // sinh(d)/d - 1 = d^2/6 + d^4/120 + ...  Taking the log of the ratio directly
  // would round the whole small-delta signal away.
  H.log_sinhc = (delta < 0.1)
    ? std::log1p(delta * delta * ((1.0 / 6.0) + delta * delta * (1.0 / 120.0)))
    : std::log(std::sinh(delta) / delta);
  return H;
}

// H_delta(sigma(x)) on the natural scale, for x <= 0 only (so the result is in
// (0, 1/2]).  The defining difference of softplus terms cancels catastrophically
// in the lower tail, so it is rearranged as an exact log1p:
//
//   (1 + e^{x+d}) / (1 + e^{x-d}) = 1 + 2 sinh(d) / (e^{-x} + e^{-d}),
//
// which also removes the overflow, since e^{-x} >= 1 replaces e^{x+d}.
inline double frq_h_lower_val(double x, const FrqH& H) {
  if (x < -700.0) return 0.0;                 // e^{-x} would overflow; H is 0
  const double r = H.two_sinh / (std::exp(-x) + H.exp_mdelta);
  return std::log1p(r) / H.two_delta;
}

// log H_delta(sigma(x)) for any x.  Above zero it uses the exact symmetry
// H(1 - z) = 1 - H(z) so that only the well-conditioned lower branch is ever
// evaluated -- the same trick as the incomplete beta's reflection, and the
// reason the survivor below needs no separate generator.
inline double frq_h_log(double x, const FrqH& H) {
  if (x <= 0.0) {
    // H ~ e^x sinh(d)/d as x -> -Inf, which stays representable long after the
    // value itself underflows.
    if (x < -700.0) return x + H.log_sinhc;
    return std::log(frq_h_lower_val(x, H));
  }
  const double lo = frq_h_lower_val(-x, H);   // = 1 - H(sigma(x)), at most 1/2
  return std::log1p(-lo);
}

// log H'(z) from log z and log(1 - z).  z(1-z) <= 1/4 always, so the exp is
// safe, and it underflows to zero in exactly the tails where H' -> sinh(d)/d.
inline double frq_h_logderiv(double lz, double lomz, const FrqH& H) {
  return H.log_sinhc - std::log1p(std::exp(lz + lomz) * H.c2);
}

// H^{-1}(y) = sinh(d y) / [sinh(d y) + sinh(d (1-y))], written as a logistic of
// the log ratio so that large delta cannot overflow the individual sinh terms
// and a near-one y keeps its complement.
inline double frq_h_inv(double y, const FrqH& H) {
  if (!H.active) return y;
  if (!(y > 0.0)) return 0.0;
  if (!(y < 1.0)) return 1.0;
  const double lr = frq_logsinh(H.delta * (1.0 - y)) - frq_logsinh(H.delta * y);
  if (lr > 0.0) { const double e = std::exp(-lr); return e / (1.0 + e); }
  return 1.0 / (1.0 + std::exp(lr));
}

// ---------------------------------------------------------------------------
// Compiled state, split by what each half actually depends on.
//
// The evidence-level availability p is the expensive half: it comes from the
// fitted completion probability h by an incomplete-beta INVERSION, an iterative
// solve costing roughly fourteen density evaluations.  It is a function of
// (K, cv, h, delta) ONLY, with cv = sv / v the rate CV; clock does not enter
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

// The (K, cv, h, delta) half -- everything the beta quantile touches, with
// cv = sv / v.  `poisson` marks the cv -> 0 limit, where beta is infinite and the kernel is
// the incomplete gamma P_K(E G(x)) (h < 1) or P_K(v x) (h = 1).  E is the
// expected total number of events available on a trial: beta p / (1 - p) for
// cv > 0, P_K^{-1}(h) in the limit, and Inf whenever h = 1.
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
  double h = 0.0;
  double lbeta_val = 0.0;
  double log_p_minus_lbeta = 0.0;
  double alpha_minus_1 = 0.0;
  double beta_minus_1 = 0.0;
  FrqH hh;
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
  FrqH hh;
  bool ok = false;
};

inline FrqAvail frq_derive_avail(double K, double cv, double h,
                                 double delta = 0.0) {
  FrqAvail s;
  if (ISNAN(K) || ISNAN(cv) || ISNAN(h)) return s;
  if (!(K >= FRQ_SHAPE_MIN) || !R_FINITE(K)) return s;
  if (!(cv >= 0.0) || !R_FINITE(cv)) return s;
  if (!(h > 0.0) || !(h <= 1.0) || !R_FINITE(h)) return s;
  s.hh = frq_h_make(delta);
  if (!s.hh.ok) return s;
  s.alpha = K;
  s.alpha_minus_1 = K - 1.0;
  s.h = h;

  if (cv < FRQ_CV_POISSON) {
    // beta -> Inf at fixed (K, v, h).  The availability p -> 0 while the
    // expected total evidence E = beta p / (1 - p) tends to the gamma quantile
    // of the (threshold-variability adjusted) completion probability.
    s.poisson = true;
    s.beta = R_PosInf;
    s.lgamma_K = R::lgammafn(K);
    if (h == 1.0) {
      s.E = R_PosInf; s.log_E = R_PosInf;
      s.p = 1.0; s.log_p = 0.0; s.om_p = 0.0;
    } else {
      const double E = R::qgamma(frq_h_inv(h, s.hh), K, 1.0, 1, 0);
      if (!(E > 0.0) || !R_FINITE(E)) return s;
      s.E = E; s.log_E = std::log(E);
      s.p = 0.0; s.log_p = R_NegInf; s.om_p = 1.0;
    }
    s.ok = true;
    return s;
  }

  const double beta = 1.0 / (cv * cv);
  if (!(beta >= FRQ_SHAPE_MIN) || !R_FINITE(beta)) return s;
  // h = 1 maps exactly to p = 1 for any shape or threshold-variability map;
  // skipping qbeta there matters because the proper boundary is the common
  // case.  Below one, 1 - p comes from the complementary quantile
  // (1 - B ~ Beta(beta, K)), since 1 - qbeta(z, K, beta) loses every digit once
  // p is near one, which is the normal case for beta < 1.  H^{-1} is symmetric,
  // so its complement is H^{-1}(1 - h).  p itself is re-solved only when it is
  // the small side.  An om_p that underflows to zero is still a valid state:
  // it only enters 1 - q(x) next to p S(x), and the tails at x = Inf return
  // h and 1 - h directly.
  double p = 1.0, om_p = 0.0;
  if (h < 1.0) {
    om_p = R::qbeta(frq_h_inv(1.0 - h, s.hh), beta, K, 1, 0);
    if (!(om_p >= 0.0) || !(om_p < 1.0)) return s;
    p = (om_p > 0.5) ? R::qbeta(frq_h_inv(h, s.hh), K, beta, 1, 0)
                     : 1.0 - om_p;
    if (!(p > 0.0) || !(p <= 1.0)) return s;
  }

  s.beta = beta;
  s.p = p;
  s.log_p = std::log(p);
  s.om_p = om_p;
  s.E = (h == 1.0) ? R_PosInf : beta * p / om_p;
  s.log_E = std::log(s.E);
  s.lbeta_val = R::lbeta(K, beta);
  s.log_p_minus_lbeta = s.log_p - s.lbeta_val;
  s.beta_minus_1 = beta - 1.0;
  s.ok = true;
  return s;
}

// Paste the availability half into the full state.  Runs only when one of
// (K, cv, h, delta) actually changes.
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
  s.h = av.h;
  s.lbeta_val = av.lbeta_val;
  s.log_p_minus_lbeta = av.log_p_minus_lbeta;
  s.alpha_minus_1 = av.alpha_minus_1;
  s.beta_minus_1 = av.beta_minus_1;
  s.hh = av.hh;
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
                          double delta = 0.0, double clock = 1.0) {
  FrqPars s;
  const FrqAvail av = frq_derive_avail(K, frq_rate_cv(sv, v), h, delta);
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
  std::uint64_t bits[4];

  bool operator==(const FrqCacheKey& other) const {
    for (int i = 0; i < 4; ++i) {
      if (bits[i] != other.bits[i]) return false;
    }
    return true;
  }
};

struct FrqCacheKeyHash {
  std::size_t operator()(const FrqCacheKey& key) const {
    std::size_t h = 0;
    for (int i = 0; i < 4; ++i) {
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
  const FrqAvail& get(double K, double cv, double h, double delta) {
    const FrqCacheKey key{{frq_cache_bits(K), frq_cache_bits(cv),
                           frq_cache_bits(h), frq_cache_bits(delta)}};
    const auto found = entries.find(key);
    if (found != entries.end()) return found->second;
    const auto inserted = entries.emplace(
      key, frq_derive_avail(K, cv, h, delta));
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
  double log_g;
};

inline double frq_registration_inv(double ratio, const FrqPars& s) {
  if (s.cv2 == 0.0)
    return -std::log1p(-ratio) / s.lambda;
  return std::expm1(-s.cv2 * std::log1p(-ratio)) / s.cv2_lambda;
}

// The same inverse given the unit survivor sr = 1 - G(x), for draws whose
// G(x) is too close to one to be held directly.
inline double frq_registration_inv_surv(double sr, const FrqPars& s) {
  if (s.cv2 == 0.0)
    return -std::log(sr) / s.lambda;
  return std::expm1(-s.cv2 * std::log(sr)) / s.cv2_lambda;
}

inline FrqState frq_state(double x, const FrqPars& s) {
  FrqState z;
  if (!R_FINITE(x)) {           // x == +Inf: fully saturated reservoir
    z.q = s.p;
    z.omq = s.om_p;
    z.log_g = R_NegInf;
    return z;
  }
  if (s.cv2 == 0.0) {
    const double e = std::exp(-s.lambda * x);
    z.q = s.p * (-std::expm1(-s.lambda * x));
    z.omq = s.om_p + s.p * e;
    z.log_g = s.log_lambda - s.lambda * x;
  } else {
    const double lx = std::log1p(s.cv2_lambda * x);
    const double log_sr = -lx * s.inv_cv2;
    const double sr = std::exp(log_sr);
    z.q = s.p * (-std::expm1(log_sr));
    z.omq = s.om_p + s.p * sr;
    // Split the exponent to avoid multiplying a potentially huge 1/cv^2 by a
    // small log1p term; lx/cv^2 has a finite exponential-limit value.
    z.log_g = s.log_lambda - lx * s.inv_cv2 - lx;
  }
  return z;
}

// ---------------------------------------------------------------------------
// Regularized incomplete beta log-evaluator: log I_x(a, b), with omx = 1 - x.
// Uses Modified Lentz's method on Gauss's continued fraction in the convergent
// regime x < (a + 1)/(a + b + 2) and the symmetry transformation
// I_x(a, b) = 1 - I_{1-x}(b, a) elsewhere. Special-cases a == 1 and b == 1
// algebraically to eliminate all loops.
// ---------------------------------------------------------------------------
inline double frq_ibeta_cf_log(double x, double a, double b, double omx, double lbeta_val) {
  const double FPMIN = 1e-300;
  const double EPS = 1e-15;
  const int MAXIT = 200;

  double f = 1.0;
  double c = 1.0;
  double d = 0.0;

  const double a1 = - (a + b) * x / (a + 1.0);
  d = 1.0 + a1 * d;
  if (std::abs(d) < FPMIN) d = FPMIN;
  c = 1.0 + a1 / c;
  if (std::abs(c) < FPMIN) c = FPMIN;
  d = 1.0 / d;
  double delta = c * d;
  f *= delta;

  for (int m = 1; m <= MAXIT; ++m) {
    const int m2 = 2 * m;
    // Step 2m (even)
    const double a_even = (m * (b - m) * x) / ((a + m2 - 1.0) * (a + m2));
    d = 1.0 + a_even * d;
    if (std::abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + a_even / c;
    if (std::abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    delta = c * d;
    f *= delta;

    // Step 2m + 1 (odd)
    const double a_odd = - ((a + m) * (a + b + m) * x) / ((a + m2) * (a + m2 + 1.0));
    d = 1.0 + a_odd * d;
    if (std::abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + a_odd / c;
    if (std::abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    delta = c * d;
    f *= delta;

    if (std::abs(delta - 1.0) < EPS) break;
  }

  const double log_pref = a * std::log(x) + b * std::log(omx) - std::log(a) - lbeta_val;
  return log_pref - std::log(f);
}

inline double frq_log_ibeta(double x, double omx, double a, double b, double lbeta_val) {
  if (x <= 0.0) return R_NegInf;
  if (x >= 1.0 || omx <= 0.0) return 0.0;

  // Exact algebraic special cases
  if (b == 1.0) {
    return a * std::log(x);
  }
  if (a == 1.0) {
    // I_x(1, b) = 1 - (1-x)^b = 1 - omx^b = -expm1(b * log(omx))
    const double z = -std::expm1(b * std::log(omx));
    return (z > 0.0) ? std::log(z) : R_NegInf;
  }

  const double thresh = (a + 1.0) / (a + b + 2.0);
  if (x < thresh) {
    return frq_ibeta_cf_log(x, a, b, omx, lbeta_val);
  } else {
    // Symmetry: I_x(a, b) = 1 - I_{omx}(b, a)
    const double log_w = frq_ibeta_cf_log(omx, b, a, x, lbeta_val);
    const double w = std::exp(log_w);
    return (w < 1.0) ? std::log1p(-w) : R_NegInf;
  }
}

// ---------------------------------------------------------------------------
// Log density / CDF / survivor of the DECISION time x = rt - t0.  Callers are
// responsible for subtracting t0; x <= 0 is handled here so that the race
// kernels can pass a negative x without branching.
// ---------------------------------------------------------------------------

// log I_{q}(alpha, beta) and log I_{1-q}(beta, alpha) -- i.e. log z and
// log(1 - z) for the untransformed CDF at x.  The threshold-variability paths
// need BOTH (their logit is the difference, their product the H' argument), so
// the pair is formed once and shared rather than recomputed per evaluator.
struct FrqLogBeta {
  double lz;
  double lomz;
};

inline FrqLogBeta frq_log_beta_pair(const FrqState& z, const FrqPars& s) {
  FrqLogBeta o;
  o.lz   = frq_log_ibeta(z.q, z.omq, s.alpha, s.beta, s.lbeta_val);
  o.lomz = frq_log_ibeta(z.omq, z.q, s.beta, s.alpha, s.lbeta_val);
  return o;
}

// The cv = 0 limit evaluates through its gamma argument X(x): v x when h = 1,
// E G(x) when h < 1 (so X(Inf) = E and P_K(E) = h before the H map).
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

inline FrqLogBeta frq_poisson_pair(const FrqPoissonArg& a, const FrqPars& s) {
  FrqLogBeta o;
  o.lz   = R::pgamma(a.X, s.alpha, 1.0, 1, 1);
  o.lomz = R::pgamma(a.X, s.alpha, 1.0, 0, 1);
  return o;
}

inline double frq_poisson_log_pdf(double x, const FrqPars& s) {
  const FrqPoissonArg a = frq_poisson_arg(x, s);
  double lp = R::dgamma(a.X, s.alpha, 1.0, 1) + a.log_dX;
  if (s.hh.active) {
    const FrqLogBeta b = frq_poisson_pair(a, s);
    lp += frq_h_logderiv(b.lz, b.lomz, s.hh);
  }
  return ISNAN(lp) ? R_NegInf : lp;
}

inline double frq_poisson_log_tail(double x, const FrqPars& s, bool lower) {
  const FrqLogBeta b = frq_poisson_pair(frq_poisson_arg(x, s), s);
  if (s.hh.active)
    return lower ? frq_h_log(b.lz - b.lomz, s.hh)
                 : frq_h_log(b.lomz - b.lz, s.hh);
  return lower ? b.lz : b.lomz;
}

inline double frq_log_pdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return R_NegInf;   // no atom at infinity in the density
  if (s.poisson) return frq_poisson_log_pdf(x, s);
  const FrqState z = frq_state(x, s);
  double lp = s.log_p_minus_lbeta + z.log_g;
  lp += frq_xlogy(s.alpha_minus_1, z.q);
  lp += frq_xlogy(s.beta_minus_1, z.omq);
  // Chain rule through the threshold-variability map.  Unlike the CDF, the base
  // density does not need the incomplete beta at all, so this is the one place
  // where delta > 0 costs real work: H' is a function of z = I_q(alpha, beta).
  if (s.hh.active) {
    const FrqLogBeta b = frq_log_beta_pair(z, s);
    lp += frq_h_logderiv(b.lz, b.lomz, s.hh);
  }
  return ISNAN(lp) ? R_NegInf : lp;
}

inline double frq_log_cdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return std::log(s.h);   // F(Inf) = h by construction
  if (s.poisson) {
    const double out = frq_poisson_log_tail(x, s, true);
    return ISNAN(out) ? R_NegInf : out;
  }
  const FrqState z = frq_state(x, s);
  double out;
  if (s.hh.active) {
    const FrqLogBeta b = frq_log_beta_pair(z, s);
    out = frq_h_log(b.lz - b.lomz, s.hh);          // logit z = log z - log(1-z)
  } else {
    out = frq_log_ibeta(z.q, z.omq, s.alpha, s.beta, s.lbeta_val);
  }
  return ISNAN(out) ? R_NegInf : out;
}

// log S(x) = log(1 - F(x)), via I_{1-q}(beta, alpha) so that the defective
// limit log(1 - h) is reached without ever subtracting two near-equal numbers.
// An unstarted accumulator (x <= 0) has survivor one, i.e. log S = 0.
//
// With delta > 0 the same guarantee survives for free: H is symmetric, so
// 1 - H(z) = H(1 - z) and the survivor is the SAME generator evaluated at
// logit(1 - z) = -logit z.  No subtraction is introduced anywhere.
inline double frq_log_surv_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NaN;             // caller distinguishes bad from certain
  if (!(x > 0.0)) return 0.0;
  if (!R_FINITE(x)) return std::log1p(-s.h);
  if (s.poisson) {
    const double out = frq_poisson_log_tail(x, s, false);
    return ISNAN(out) ? R_NaN : out;
  }
  const FrqState z = frq_state(x, s);
  double out;
  if (s.hh.active) {
    const FrqLogBeta b = frq_log_beta_pair(z, s);
    out = frq_h_log(b.lomz - b.lz, s.hh);
  } else {
    out = frq_log_ibeta(z.omq, z.q, s.beta, s.alpha, s.lbeta_val);
  }
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

// Front-end memo for every call site.  It keeps the two halves on separate
// keys, so a row that moves only the rate -- a trend or a covariate on v,
// against constant shapes -- redoes one log rather than one beta quantile.
// The likelihood adapters hand it the per-particle FrqCache; the exported
// entry points run without one and memoise locally.
struct FrqMemo {
  FrqCache* shared = nullptr;
  double K = R_NaN, cv = R_NaN, h = R_NaN, delta = R_NaN;
  double v = R_NaN, clock = R_NaN;
  // The memo owns its copy rather than holding a pointer into the shared
  // cache: the copy runs only when the shape key changes -- the same branch
  // that would otherwise pay a hash lookup or a beta quantile -- and it keeps
  // the memo's lifetime independent of when the cache is cleared.
  FrqAvail avail;
  FrqPars value;

  explicit FrqMemo(FrqCache* shared_cache = nullptr) : shared(shared_cache) {}

  const FrqPars& get(double k, double vv, double sv, double hh,
                     double d = 0.0, double cu = 1.0) {
    // NaN never compares equal, so a NaN parameter always re-derives, and the
    // NaN-initialised keys make the first call always miss.
    const double c = frq_rate_cv(sv, vv);
    const bool same_avail = k == K && c == cv && hh == h && d == delta;
    if (same_avail && vv == v && cu == clock) return value;
    if (!same_avail) {
      K = k; cv = c; h = hh; delta = d;
      avail = (shared != nullptr) ? shared->get(k, c, hh, d)
                                  : frq_derive_avail(k, c, hh, d);
      frq_copy_avail(value, avail);
    }
    v = vv; clock = cu;
    frq_set_rate(value, avail.ok, vv, cu);
    return value;
  }
};

// Race adapters used by the particle likelihood dispatch.  Definitions live
// in model_FRQ.cpp so this header remains safe for model_rng.cpp.
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

#endif  // EMC2_MODEL_FRQ_H

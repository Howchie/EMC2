#ifndef EMC2_MODEL_FRQ_H
#define EMC2_MODEL_FRQ_H

// ============================================================================
// FRQ -- the Finite Reservoir Quorum process.
//
// GENERATIVE STORY.  One accumulator draws on N potential evidence units and
// responds once K distinct units have registered.  Each unit is available
// independently with probability p, and an available unit registers at an
// Exponential(lambda) latency; an unavailable unit never registers at all.  The
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
// g(x) = G'(x) = lambda e^{-lambda x}.
// Nothing here needs a first-passage solve or a quadrature: the entire
// likelihood is four Rmath special functions.  The continuous relaxation is
// just alpha, beta > 0 rather than positive integers; every formula above is
// already valid there, because nothing above uses integrality -- only that
// I_q(alpha, beta) and the Beta density exist, which they do for any
// alpha, beta > 0.  The FITTED model bounds both
// shapes below at 1 (R/model_FRQ.R): the sub-one corner is mathematically
// fine but the (h, tau) coordinates are not representable there in double
// precision, so it is not offered.
//
// PARAMETERISATION.  Raw (p, lambda) are poor estimation coordinates, so the
// exposed parameters are (alpha, beta, h, tau, t0) where h is the eventual
// completion probability and tau is the CONDITIONAL MEDIAN decision time,
// P(T <= tau | T < Inf) = 1/2.  frq_derive() inverts that exactly:
//
//   p      = I^{-1}_{H^{-1}(h)}(alpha, beta)
//   u_r    = I^{-1}_{H^{-1}(r h)}(alpha, beta)         with r = FRQ_QUANTILE
//   lambda = -log1p(-u_r / p) / tau
//
// (H is the identity at the default delta = 0, so those reduce to qbeta(h, .)
// and qbeta(r h, .).)
//
// r is hard-coded at 1/2.  Exposing it would add a parameter that is not
// identified separately from tau: (r, tau) enter the likelihood only through
// the single point they pin on F, so any change in r can be absorbed exactly
// by a change in tau.
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
//     so (h, tau) keep their exact meanings for STILL exactly two qbeta calls.
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

using namespace Rcpp;

// Conditional-quantile level defining tau.  1/2 = the conditional median.
// Changing this changes the meaning of tau in every fitted model, so it is a
// constant rather than an argument.
constexpr double FRQ_QUANTILE = 0.5;

// Shapes below this are treated as degenerate.  The fitted model keeps
// alpha, beta >= 1 (R/model_FRQ.R bounds), so this floor is only a guard for
// the exported dfrq/pfrq/frq_rate entry points, which take shapes directly.
constexpr double FRQ_SHAPE_MIN = 1e-8;

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
// (alpha, beta, h, tau, delta) -> (p, lambda).  `ok` is false for any parameter
// combination that does not define a proper accumulator; every evaluator then
// returns log 0 for the density and lets the caller floor the trial, which is
// what the rest of the race machinery expects from an invalid row.
// ---------------------------------------------------------------------------
struct FrqPars {
  double alpha = 0.0;
  double beta = 0.0;
  double p = 0.0;
  double lambda = 0.0;
  double lbeta_val = 0.0;
  FrqH hh;
  bool ok = false;
};

inline FrqPars frq_derive(double alpha, double beta, double h, double tau,
                          double delta = 0.0) {
  FrqPars s;
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(h) || ISNAN(tau)) return s;
  if (!(alpha >= FRQ_SHAPE_MIN) || !R_FINITE(alpha)) return s;
  if (!(beta >= FRQ_SHAPE_MIN) || !R_FINITE(beta)) return s;
  if (!(h > 0.0) || !(h <= 1.0)) return s;
  if (!(tau > 0.0) || !R_FINITE(tau)) return s;
  s.hh = frq_h_make(delta);
  if (!s.hh.ok) return s;

  s.alpha = alpha;
  s.beta = beta;
  // p is the saturation level of q(x): the registered fraction the reservoir
  // can ever reach.  h = H(I_p(alpha, beta)) inverts to it exactly, in two
  // steps: undo the threshold-variability map, then undo the incomplete beta.
  s.p = R::qbeta(frq_h_inv(h, s.hh), alpha, beta, /*lower_tail=*/1, /*log_p=*/0);
  // If h is defective, p must remain strictly below one.  qbeta can round a
  // representable but extreme h to one; accepting that value would silently
  // turn the requested defective distribution into a proper one.
  if (ISNAN(s.p) || !(s.p > 0.0) || !(s.p <= 1.0) ||
      (h < 1.0 && !(s.p < 1.0))) return s;

  // The r-quantile of the CONDITIONAL distribution sits at F = r h, so on the
  // Beta scale it is the H^{-1}(r h)-quantile.  It is strictly below p whenever
  // r < 1 (H is strictly increasing), which is what makes the log1p below
  // finite.
  const double u_r = R::qbeta(frq_h_inv(FRQ_QUANTILE * h, s.hh), alpha, beta, 1, 0);
  if (ISNAN(u_r) || !(u_r >= 0.0)) return s;
  const double gap = s.p - u_r;
  // A vanishing gap means the two quantiles have collapsed into each other in
  // double precision -- lambda is then unresolvable rather than merely large,
  // so reject instead of returning a garbage rate.
  if (!(gap > 0.0)) return s;

  s.lambda = -std::log1p(-u_r / s.p) / tau;
  if (ISNAN(s.lambda) || !(s.lambda > 0.0) || !R_FINITE(s.lambda)) return s;
  s.lbeta_val = R::lbeta(alpha, beta);
  s.ok = true;
  return s;
}

// q(x) = p (1 - e^{-lambda x}) together with its cancellation-free complement
// 1 - q(x) = (1 - p) + p e^{-lambda x}.  x == Inf gives the saturated pair
// (p, 1 - p), which is exactly what makes F(Inf) = h fall out of the ordinary
// CDF code path rather than needing a special case.
struct FrqState {
  double q;
  double omq;
};

inline FrqState frq_state(double x, const FrqPars& s) {
  FrqState z;
  if (!R_FINITE(x)) {           // x == +Inf: fully saturated reservoir
    z.q = s.p;
    z.omq = 1.0 - s.p;
    return z;
  }
  const double e = std::exp(-s.lambda * x);
  z.q = s.p * (-std::expm1(-s.lambda * x));
  z.omq = (1.0 - s.p) + s.p * e;
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

inline double frq_log_pdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return R_NegInf;   // no atom at infinity in the density
  const FrqState z = frq_state(x, s);
  double lp = std::log(s.p) + std::log(s.lambda) - s.lambda * x
              - s.lbeta_val;
  lp += frq_xlogy(s.alpha - 1.0, z.q);
  lp += frq_xlogy(s.beta - 1.0, z.omq);
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

// A single-entry memo for frq_derive().  The two qbeta inversions are by far
// the most expensive part of an FRQ likelihood evaluation, and the batch
// kernels walk a compressed parameter matrix in which consecutive rows very
// often share (alpha, beta, h, tau) -- shapes are usually constant across
// accumulators, and h/tau only vary across design cells.  Keyed on exact bit
// equality, so it can never change an answer.
struct FrqMemo {
  double alpha = R_NaN, beta = R_NaN, h = R_NaN, tau = R_NaN, delta = R_NaN;
  FrqPars value;

  const FrqPars& get(double a, double b, double hh, double t, double d = 0.0) {
    if (a == alpha && b == beta && hh == h && t == tau && d == delta) return value;
    alpha = a; beta = b; h = hh; tau = t; delta = d;
    value = frq_derive(a, b, hh, t, d);
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

#ifndef EMC2_MODEL_FRQ_H
#define EMC2_MODEL_FRQ_H

// ============================================================================
// FRQ -- the Finite Reservoir Quorum process (Math/FRQ.tex).
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
// already valid there (Math/FRQ.tex Sec. 12).  The FITTED model bounds both
// shapes below at 1 (R/model_FRQ.R): the sub-one corner is mathematically
// fine but the (h, tau) coordinates are not representable there in double
// precision, so it is not offered.
//
// PARAMETERISATION.  Raw (p, lambda) are poor estimation coordinates, so the
// exposed parameters are (alpha, beta, h, tau, t0) where h is the eventual
// completion probability and tau is the CONDITIONAL MEDIAN decision time,
// P(T <= tau | T < Inf) = 1/2.  frq_derive() inverts that exactly:
//
//   p      = I^{-1}_h(alpha, beta)                     = qbeta(h, alpha, beta)
//   u_r    = I^{-1}_{r h}(alpha, beta)                 with r = FRQ_QUANTILE
//   lambda = -log1p(-u_r / p) / tau
//
// r is hard-coded at 1/2 per Sec. 8; exposing it would add a parameter that is
// not identified separately from tau.
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
// constant rather than an argument; see Math/FRQ.tex Sec. 8.
constexpr double FRQ_QUANTILE = 0.5;

// Shapes below this are treated as degenerate.  The fitted model keeps
// alpha, beta >= 1 (R/model_FRQ.R bounds), so this floor is only a guard for
// the exported dfrq/pfrq/frq_rate entry points, which take shapes directly.
constexpr double FRQ_SHAPE_MIN = 1e-8;

// xlogy(a, y) = a log y, with the a == 0 branch returning 0 rather than the
// NaN of 0 * (-Inf).  alpha == 1 or beta == 1 hits this on every call.
inline double frq_xlogy(double a, double y) {
  if (a == 0.0) return 0.0;
  if (y <= 0.0) return (a > 0.0) ? R_NegInf : R_PosInf;
  return a * std::log(y);
}

// ---------------------------------------------------------------------------
// (alpha, beta, h, tau) -> (p, lambda).  `ok` is false for any parameter
// combination that does not define a proper accumulator; every evaluator then
// returns log 0 for the density and lets the caller floor the trial, which is
// what the rest of the race machinery expects from an invalid row.
// ---------------------------------------------------------------------------
struct FrqPars {
  double alpha = 0.0;
  double beta = 0.0;
  double p = 0.0;
  double lambda = 0.0;
  bool ok = false;
};

inline FrqPars frq_derive(double alpha, double beta, double h, double tau) {
  FrqPars s;
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(h) || ISNAN(tau)) return s;
  if (!(alpha >= FRQ_SHAPE_MIN) || !R_FINITE(alpha)) return s;
  if (!(beta >= FRQ_SHAPE_MIN) || !R_FINITE(beta)) return s;
  if (!(h > 0.0) || !(h <= 1.0)) return s;
  if (!(tau > 0.0) || !R_FINITE(tau)) return s;

  s.alpha = alpha;
  s.beta = beta;
  // p is the saturation level of q(x): the registered fraction the reservoir
  // can ever reach.  h = I_p(alpha, beta) inverts to it exactly.
  s.p = R::qbeta(h, alpha, beta, /*lower_tail=*/1, /*log_p=*/0);
  if (ISNAN(s.p) || !(s.p > 0.0) || !(s.p <= 1.0)) return s;

  // The r-quantile of the CONDITIONAL distribution sits at F = r h, so on the
  // Beta scale it is the (r h)-quantile.  It is strictly below p whenever
  // r < 1, which is what makes the log1p below finite.
  const double u_r = R::qbeta(FRQ_QUANTILE * h, alpha, beta, 1, 0);
  if (ISNAN(u_r) || !(u_r >= 0.0)) return s;
  const double gap = s.p - u_r;
  // A vanishing gap means the two quantiles have collapsed into each other in
  // double precision -- lambda is then unresolvable rather than merely large,
  // so reject instead of returning a garbage rate.
  if (!(gap > 0.0)) return s;

  s.lambda = -std::log1p(-u_r / s.p) / tau;
  if (ISNAN(s.lambda) || !(s.lambda > 0.0) || !R_FINITE(s.lambda)) return s;
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
// Log density / CDF / survivor of the DECISION time x = rt - t0.  Callers are
// responsible for subtracting t0; x <= 0 is handled here so that the race
// kernels can pass a negative x without branching.
// ---------------------------------------------------------------------------

inline double frq_log_pdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  if (!R_FINITE(x)) return R_NegInf;   // no atom at infinity in the density
  const FrqState z = frq_state(x, s);
  double lp = std::log(s.p) + std::log(s.lambda) - s.lambda * x
              - R::lbeta(s.alpha, s.beta);
  lp += frq_xlogy(s.alpha - 1.0, z.q);
  lp += frq_xlogy(s.beta - 1.0, z.omq);
  return ISNAN(lp) ? R_NegInf : lp;
}

inline double frq_log_cdf_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NegInf;
  if (!(x > 0.0)) return R_NegInf;
  const FrqState z = frq_state(x, s);
  const double out = R::pbeta(z.q, s.alpha, s.beta, /*lower_tail=*/1, /*log_p=*/1);
  return ISNAN(out) ? R_NegInf : out;
}

// log S(x) = log(1 - F(x)), via I_{1-q}(beta, alpha) so that the defective
// limit log(1 - h) is reached without ever subtracting two near-equal numbers.
// An unstarted accumulator (x <= 0) has survivor one, i.e. log S = 0.
inline double frq_log_surv_dt(double x, const FrqPars& s) {
  if (!s.ok) return R_NaN;             // caller distinguishes bad from certain
  if (!(x > 0.0)) return 0.0;
  const FrqState z = frq_state(x, s);
  const double out = R::pbeta(z.omq, s.beta, s.alpha, /*lower_tail=*/1, /*log_p=*/1);
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
  double alpha = R_NaN, beta = R_NaN, h = R_NaN, tau = R_NaN;
  FrqPars value;

  const FrqPars& get(double a, double b, double hh, double t) {
    if (a == alpha && b == beta && hh == h && t == tau) return value;
    alpha = a; beta = b; h = hh; tau = t;
    value = frq_derive(a, b, hh, t);
    return value;
  }
};

#endif  // EMC2_MODEL_FRQ_H

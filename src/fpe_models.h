#ifndef fpe_models_h
#define fpe_models_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The four leaky-accumulation models, expressed for the Fokker-Planck solver in
// fpe_solver.h.  Each model supplies five things:
//
//   x_lo(t)            lower edge of the computational domain
//   x_hi(t)            upper edge (the absorbing barrier)
//   B()                diffusion coefficient (constant in x -- see below)
//   drift(x, t)        A(x,t)
//   length(t)          L(t) = x_hi(t) - x_lo(t)
//   length_prime(t)    L'(t)
//   static_op()        true when the operator is time-invariant (fixed bounds)
//
// x_lo is a function of t rather than a constant because the bounded model can
// collapse BOTH barriers (see FPE_ModelBoundedOU).  A moving lower edge adds a
// frame term -x_lo'(t)/L to the transformed drift; each model folds that into
// its own atil_affine, so there is no separate x_lo_prime accessor to keep in
// step with it.  For the one-boundary models x_lo is constant and there is no
// such term at all.
//
// Only TWO drift laws are needed, because the geometric models reduce by a state
// log-transform Y = log X (Ito):
//
//   BM        A = mu                       B = sigma
//   OU        A = -lambda*(x - theta)      B = sigma
//   GBM       Y = log X  =>  BM with A = mu - sigma^2/2,  boundary log a(t)
//   Gompertz  Y = log X  =>  OU with lambda = alpha, sigma = beta,
//                            theta = log(K) - beta^2/(2 alpha),  boundary log a(t)
//
// (The Gompertz map is the same one as transform_gompertz_to_ou() in
// utils_reducible_diffusion.h:119; restated here so this translation unit does
// not have to include the Volterra headers, which own non-inline free functions.)
//
// Parameterisation is deliberately identical to the Volterra entry points, so
// parameters transfer 1:1:
//   * boundary  b(t) = binf + (|b0| - |binf|) * exp(-(t/tau)^pow)
//   * z0        UPPER limit of a Uniform(0, z0) start point; z0 == 0 is a point
//               start at 0.  (Geometric models offset by start_floor, matching
//               simulate_gbm_hit_times_bb / simulate_gompertz_hit_times_bb.)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <cmath>
#include <algorithm>
#include "fpe_solver.h"

namespace fpe {

#ifndef FPE_EPS
#define FPE_EPS 1e-12
#endif

// ---------------------------------------------------------------------------
// Collapsing boundary -- a self-contained module.
//
// This is deliberately the ONLY place that knows the functional form of b(t).
// It answers three questions and nothing else: what is the boundary at time t,
// what is its derivative, and is it constant (which enables the one-time
// factorisation in fpe_solve).  It knows no parameter NAMES: the model
// definition owns the mapping from named parameters to (kind, b0, binf, p1, p2)
// and does that mapping in its Ttransform.  Keeping it ignorant is what lets
// the same module serve the solver and the simulator, and what lets a new
// collapse form be added without touching either.
//
// FPE_BND_WEIBULL is the form the package already uses, matching
// exp_decay_scalar() (utils_reducible_diffusion.h:56) bit for bit:
//   b(t) = binf + (|b0| - |binf|) * exp(-(t/p1)^p2)
// The other kinds are declared but not yet implemented; adding one is a case in
// the two switches below and a p_types entry in the model definition.
// ---------------------------------------------------------------------------
enum FPE_BoundaryKind {
  FPE_BND_FIXED = 0,
  FPE_BND_WEIBULL = 1,                 // b = binf + (b0-binf)*exp(-(t/p1)^p2)
  FPE_BND_EXPONENTIAL = 2,             // b = binf + (b0-binf)*exp(-t/p1)
  FPE_BND_LINEAR_ADDITIVE = 3,         // b = binf + (b0-binf)*max(0, 1 - t/p1)
  FPE_BND_LINEAR_MULTIPLICATIVE = 4    // b = binf + (b0-binf)/(1 + t/p1)
};

struct FPE_Boundary {
  int kind = FPE_BND_FIXED;
  double b0 = 1.0, binf = 1.0, tau = 0.0, pw = 0.0;
  bool fixed = true;
  bool log_state = false;      // solve in Y = log X, so the barrier is log a(t)

  // Generic entry point: (kind, b0, binf, p1, p2).  A form that degenerates to a
  // constant reports fixed = true, so callers never have to special-case it.
  void set_kind(int kind_, double b0_, double binf_, double p1, double p2,
                bool log_state_) {
    kind = kind_;
    b0 = std::abs(b0_);
    binf = std::abs(binf_);
    tau = p1;
    pw = p2;
    log_state = log_state_;
    switch (kind) {
      case FPE_BND_WEIBULL:
        fixed = (std::abs(binf - b0) <= FPE_EPS) || tau <= 0.0 || pw <= 0.0;
        break;
      case FPE_BND_EXPONENTIAL:
      case FPE_BND_LINEAR_ADDITIVE:
      case FPE_BND_LINEAR_MULTIPLICATIVE:
        // Same degeneracy test minus the shape exponent, which these forms do
        // not have.  Reporting fixed = true here is not a shortcut: it restores
        // the one-time factorisation for a collapse that does not collapse, so
        // a sampler that wanders to binf == b0 does not silently pay 3x.
        fixed = (std::abs(binf - b0) <= FPE_EPS) || tau <= 0.0;
        break;
      case FPE_BND_FIXED:
      default:
        fixed = true;
        binf = b0;
        break;
    }
  }

  // Back-compatible Weibull entry point used by the validation exports.
  void set(double b0_, double binf_, double tau_, double pow_, bool log_state_) {
    set_kind(FPE_BND_WEIBULL, b0_, binf_, tau_, pow_, log_state_);
  }

  // physical boundary b(t) (before any log transform).  Called once per TIME
  // STEP, not per cell, so the switch costs nothing.
  double b(double t) const {
    if (fixed) return b0;
    switch (kind) {
      case FPE_BND_WEIBULL: {
        const double amp = b0 - binf;
        if (!(t > 0.0)) return b0;
        const double s = std::exp(pw * (std::log(t) - std::log(tau)));
        return binf + amp * std::exp(-s);
      }
      case FPE_BND_EXPONENTIAL: {
        if (!(t > 0.0)) return b0;
        return binf + (b0 - binf) * std::exp(-t / tau);
      }
      case FPE_BND_LINEAR_ADDITIVE: {
        if (!(t > 0.0)) return b0;
        // Held at binf past t = tau rather than continuing down through it.  A
        // boundary that keeps falling would cross the start-point range and
        // then the lower domain face, which is not a model, it is a bug.
        const double f = 1.0 - t / tau;
        return binf + (b0 - binf) * ((f > 0.0) ? f : 0.0);
      }
      case FPE_BND_LINEAR_MULTIPLICATIVE: {
        if (!(t > 0.0)) return b0;
        return binf + (b0 - binf) / (1.0 + t / tau);
      }
      default:
        return b0;
    }
  }

  // db/dt.  Weibull: -amp * (p/t) * s * exp(-s),  s = (t/tau)^p
  double b_prime(double t) const {
    if (fixed) return 0.0;
    switch (kind) {
      case FPE_BND_WEIBULL: {
        const double amp = b0 - binf;
        const double tt = std::max(t, 1e-9);
        const double s = std::exp(pw * (std::log(tt) - std::log(tau)));
        if (!std::isfinite(s)) return 0.0;
        const double d = -amp * (pw / tt) * s * std::exp(-s);
        return std::isfinite(d) ? d : 0.0;
      }
      case FPE_BND_EXPONENTIAL: {
        if (!(t > 0.0)) return -(b0 - binf) / tau;
        const double d = -(b0 - binf) / tau * std::exp(-t / tau);
        return std::isfinite(d) ? d : 0.0;
      }
      case FPE_BND_LINEAR_ADDITIVE:
        // Kinked at t = tau.  CN is second order on each smooth piece and the
        // kink costs at most one step's worth of order there, which is well
        // inside the discretisation error the graded schedule already carries.
        return (t < tau) ? -(b0 - binf) / tau : 0.0;
      case FPE_BND_LINEAR_MULTIPLICATIVE: {
        if (!(t > 0.0)) return -(b0 - binf) / tau;
        const double denom = 1.0 + t / tau;
        const double d = -(b0 - binf) / (tau * denom * denom);
        return std::isfinite(d) ? d : 0.0;
      }
      default:
        return 0.0;
    }
  }

  double a(double t) const {
    const double v = b(t);
    return log_state ? std::log(std::max(v, 1e-300)) : v;
  }
  double a_prime(double t) const {
    if (fixed) return 0.0;
    const double v = b(t);
    return log_state ? b_prime(t) / std::max(v, 1e-300) : b_prime(t);
  }
};

// ---------------------------------------------------------------------------
// Constant-drift model: BM directly, GBM after Y = log X.
// ---------------------------------------------------------------------------
struct FPE_ModelBM {
  // One absorbing boundary (xi = 1); xi = 0 is the no-flux far field.
  static constexpr bool lower_absorbing = false;
  static constexpr bool symmetric_mesh  = false;

  double A = 0.0;          // mu, or mu - sigma^2/2 in log state
  double sigma = 1.0;
  double xlo = -1.0;
  FPE_Boundary bnd;

  double x_lo(double /*t*/) const { return xlo; }
  double x_hi(double t) const { return bnd.a(t); }
  double B() const { return sigma; }
  double drift(double /*x*/, double /*t*/) const { return A; }
  double length(double t) const { return bnd.a(t) - xlo; }
  double length_prime(double t) const { return bnd.a_prime(t); }
  bool static_op() const { return bnd.fixed; }

  // Atil(xi,t) = ( A - xi*L'(t) ) / L(t)
  void atil_affine(double /*t*/, double L, double Lp, double& a0, double& a1) const {
    a0 = A / L;
    a1 = -Lp / L;
  }
};

// ---------------------------------------------------------------------------
// Mean-reverting model: OU directly, Gompertz after Y = log X.
// ---------------------------------------------------------------------------
// Parameterised as dX = (v - lambda X) dt + sigma dW, i.e. INPUT and LEAK, not
// (lambda, theta).  These are the same process wherever both are defined --
// v = lambda*theta -- but v is regular at lambda = 0 and theta is not.  At
// lambda = 0 the drift is simply the constant v and the coefficients below
// become FPE_ModelBM's exactly, so the Brownian race is a genuine member of this
// model rather than a dispatch branch to some other kernel.  That is what keeps
// the k = 0 comparison against the analytic Wald an honest cross-validation.
struct FPE_ModelOU {
  static constexpr bool lower_absorbing = false;
  static constexpr bool symmetric_mesh  = false;

  double v = 0.0;              // input rate; theta = v / lambda where defined
  double lambda = 1.0;         // leak
  double sigma = 1.0;
  double xlo = -1.0;
  FPE_Boundary bnd;

  // Convenience for callers that think in (lambda, theta) -- the validation
  // entry points and the Gompertz reduction.
  void set_lambda_theta(double lambda_, double theta_) {
    lambda = lambda_;
    v = lambda_ * theta_;
  }
  double theta() const { return (lambda > FPE_EPS) ? v / lambda : 0.0; }

  double x_lo(double /*t*/) const { return xlo; }
  double x_hi(double t) const { return bnd.a(t); }
  double B() const { return sigma; }
  double drift(double x, double /*t*/) const { return v - lambda * x; }
  double length(double t) const { return bnd.a(t) - xlo; }
  double length_prime(double t) const { return bnd.a_prime(t); }
  bool static_op() const { return bnd.fixed; }

  // Atil(xi,t) = ( v - lambda*(x_lo + xi*L) - xi*L'(t) ) / L(t)
  void atil_affine(double /*t*/, double L, double Lp, double& a0, double& a1) const {
    a0 = (v - lambda * xlo) / L;
    a1 = (-lambda * L - Lp) / L;
  }
};

// ---------------------------------------------------------------------------
// Bounded OU -- the Smith & Ratcliff (2004) two-choice OU diffusion.
//
//   dX = ( v + beta*(anchor - X) ) dt + sigma dW,   X(0) = z,
//   absorbing at BOTH 0 and a; the barrier hit determines the response.
//
// This is the DDM with leak.  `beta` is the decay rate; at beta = 0 the drift
// is the constant v and the coefficients below become FPE_ModelBM's exactly, so
// the model reduces to the Wiener diffusion and can be checked against the
// package's own Navarro-Fuss DDM as an oracle -- while still travelling this
// solver's code path rather than that one's.
//
// `anchor` is the point the decay pulls toward.  Smith & Ratcliff take it to be
// the STARTING POINT z (their Appendix p. 42, and footnote 2 argues for it over
// the alternative), which is the default here.  Note that zero is NOT an
// available alternative in this parameterisation -- it is one of the response
// boundaries -- so the only other sensible anchor is the midpoint a/2.  Keeping
// it a separate field rather than hard-coding z costs one addition in
// atil_affine and leaves that choice open; it also matters for cost, because an
// anchor that does not depend on z leaves the operator independent of the start
// point and collapses across-trial start-point variability into a single solve.
//
// COLLAPSING BOUNDS.  Both barriers move, symmetrically, toward the MIDPOINT:
//
//   upper(t) = mid + h(t),   lower(t) = mid - h(t),   mid = a/2,
//
// with the half-separation h(t) supplied by the same FPE_Boundary module the
// race models use, so every collapse form (Weibull, exponential, the two linear
// ones) is available here for free and behaves identically.  h(0) = a/2, so the
// separation runs from a down to the asymptote.
//
// Collapsing toward the midpoint rather than toward a fixed level is the only
// choice that leaves the model well defined: an asymmetric collapse would have
// to specify which barrier moves, and a collapse toward one of the barriers is
// not a collapse at all.  It also keeps the geometry symmetric under the
// response relabelling (v, Z) -> (-v, 1-Z), which is the property the
// zero-drift/unbiased test in test-fpe-bounded-ou.R checks.
//
// The half-separation is floored at FPE_BOU_MIN_SEP * a rather than being
// allowed to reach zero.  Nothing in the SPACE discretisation cares -- the mesh
// lives on the normalised xi in [0,1] and follows the barriers down -- but
// D(t) = (B/L)^2/2 grows like 1/L^2, so the time step stops resolving the flux
// long before the width does.  By the time the barriers are 2% of a apart the
// remaining survivor crosses in (0.02a/s)^2 of a second, so the floor changes
// where that last sliver of mass is absorbed by well under a millisecond.
// ---------------------------------------------------------------------------
constexpr double FPE_BOU_MIN_SEP = 0.02;

struct FPE_ModelBoundedOU {
  static constexpr bool lower_absorbing = true;
  static constexpr bool symmetric_mesh  = true;

  double v = 0.0;              // drift rate (xi in Smith & Ratcliff)
  double beta = 0.0;           // decay / leak; 0 => Wiener diffusion
  double anchor = 0.0;         // point the decay pulls toward (default: z)
  double sigma = 1.0;          // within-trial sd (s)
  double mid = 0.5;            // (fixed) midpoint the barriers collapse toward
  FPE_Boundary bnd;            // HALF-separation h(t); barriers are mid +/- h

  // The one place the (a, a_inf) user parameterisation meets the (mid, h)
  // internal one.  a is the separation at t = 0 and a_inf its asymptote, both in
  // the same units, so a design reads as "a collapses from a to a_inf".
  void set_separation(double a, int kind, double a_inf, double tau, double pw) {
    mid = 0.5 * a;
    if (kind == FPE_BND_FIXED) {
      // 0.5*a is exact in binary and 2*(0.5*a) == a, so a fixed bound reproduces
      // the pre-collapse geometry (x_lo = 0, x_hi = a, L = a) bit for bit.
      bnd.set_kind(FPE_BND_FIXED, mid, mid, 0.0, 0.0, false);
    } else {
      const double lo = FPE_BOU_MIN_SEP * a;
      bnd.set_kind(kind, mid, 0.5 * std::max(std::abs(a_inf), lo), tau, pw,
                   false);
    }
  }

  double half(double t) const { return bnd.a(t); }
  double x_lo(double t) const { return mid - bnd.a(t); }
  double x_hi(double t) const { return mid + bnd.a(t); }
  double B() const { return sigma; }
  double drift(double x, double /*t*/) const { return v + beta * (anchor - x); }
  double length(double t) const { return 2.0 * bnd.a(t); }
  double length_prime(double t) const { return 2.0 * bnd.a_prime(t); }
  bool static_op() const { return bnd.fixed; }

  // A(x) = (v + beta*anchor) - beta*x, so with x = x_lo(t) + xi*L(t),
  //   Atil(xi,t) = ( v + beta*(anchor - x_lo) - x_lo' - xi*(beta*L + L') ) / L
  // The -x_lo' term is the frame velocity of the lower edge and is what a fixed
  // lower bound does not have.  Here x_lo' = -h' = -L'/2 exactly, so it enters as
  // +L'/2: as the barriers close in, the moving frame adds an outward drift on
  // the lower side, which is the correct sign -- mass that is standing still in
  // the lab frame is moving toward xi = 0 relative to a rising floor.
  void atil_affine(double /*t*/, double L, double Lp, double& a0, double& a1) const {
    a0 = (v + beta * (anchor - (mid - 0.5 * L)) + 0.5 * Lp) / L;
    a1 = (-beta * L - Lp) / L;
  }
};

// ---------------------------------------------------------------------------
// Domain lower edge.
//
// Sized as (lowest reachable MEAN) - FPE_NSD * (diffusive sd), NOT as a blanket
// multiple of |A|*t_max + sigma*sqrt(t_max).  The distinction matters a great
// deal: the whole domain is spanned by a uniform mesh, so every unit of slack at
// the bottom is resolution stolen from the absorbing barrier, where all of the
// accuracy actually lives.  With the naive sizing, BM(mu=1, sigma=1, t_max=2)
// gave L = 28.3 and a cdf error of 2e-3 at M = 256.
//
// Drift is only allowed to *lower* the floor (min(0, A*t_max)); an upward drift
// does not license a deeper domain.  FPE_NSD = 6 truncates ~1e-9 of the mass.
// ---------------------------------------------------------------------------
constexpr double FPE_NSD = 6.0;

// z_min is the lowest reachable start point IN SOLVER COORDINATES (already
// log-transformed for the geometric models).
inline double fpe_x_lo_bm(double z_min, double A, double sigma, double t_max) {
  return z_min + std::min(0.0, A * t_max) - FPE_NSD * sigma * std::sqrt(t_max);
}

// (v, lambda) form, so that lambda = 0 is a regular point: the OU mean travels
// monotonically from z_min toward theta = v/lambda, so the lowest mean ever
// attained is min(z_min, theta) -- and as lambda -> 0 that limit is the BM
// expression min(z_min, z_min + v*t_max), which is what the second branch is.
inline double fpe_x_lo_ou(double z_min, double v, double lambda,
                          double sigma, double t_max) {
  const double var = (lambda > FPE_EPS)
    ? (1.0 - std::exp(-2.0 * lambda * t_max)) / (2.0 * lambda)
    : t_max;
  const double lo_mean = (lambda > FPE_EPS)
    ? std::min(z_min, v / lambda)
    : z_min + std::min(0.0, v * t_max);
  return lo_mean - FPE_NSD * sigma * std::sqrt(var);
}

// Number of mesh cells the seeding Gaussian's sd should span.  Smaller => the
// solver takes over sooner => less reliance on the frozen-coefficient formula.
constexpr double FPE_SEED_CELLS = 4.0;

// Normal probability of [a, b], branched on sign so the two erfc values are
// never subtracted when both are near 1.
inline double fpe_norm_mass(double a, double b, double mean, double sd) {
  const double za = (a - mean) / (sd * M_SQRT2);
  const double zb = (b - mean) / (sd * M_SQRT2);
  if (za >= 0.0)      return 0.5 * (std::erfc(za) - std::erfc(zb));   // both above
  else if (zb <= 0.0) return 0.5 * (std::erfc(-zb) - std::erfc(-za)); // both below
  return 1.0 - 0.5 * std::erfc(zb) - 0.5 * std::erfc(-za);            // straddling
}

// ---------------------------------------------------------------------------
// Initial condition, uniform start on [z_lo, z_hi] (z_lo == z_hi => point).
//
// Uniform: exact cell-overlap averaging, so the discontinuous edges of the
// uniform density do not cost an order of accuracy.
//
// Point:   do NOT seed a delta.  Start at a short t_seed from the analytic
//          method-of-images absorbed Gaussian in a frozen-coefficient local
//          frame, which already carries the mass absorbed before t_seed:
//
//   p(x,t) = phi(x; z + A t, B^2 t)
//            - exp(2 A (a - z)/B^2) * phi(x; 2a - z + A t, B^2 t)
//
//          NB the image mean is 2a - z + A t (reflect the START point about the
//          barrier, then drift), not 2a - (z + A t).
//
// Returns t0; fills q (length M) with the normalised sub-density q = L * p.
// ---------------------------------------------------------------------------
// `t_seed_force` > 0 overrides the seed time chosen below.  Lane-batched callers
// need it: a batch marches on ONE clock, so every lane has to start from the
// same t, and the analytic seed is equally valid at any sufficiently small time.
// The batch picks the smallest t_seed any of its lanes would have chosen, so no
// lane is seeded later -- and therefore less accurately -- than it would have
// been on its own.
template <class Model>
inline double fpe_seed(const Model& m, double z_lo, double z_hi,
                       const FPE_Mesh& g, double t_max, std::vector<double>& q,
                       double* absorbed_lower = nullptr,
                       double t_seed_force = -1.0) {
  const int M = g.M;
  q.assign(M, 0.0);
  if (absorbed_lower != nullptr) *absorbed_lower = 0.0;

  if (z_hi - z_lo > FPE_EPS) {
    // ---- uniform start ----
    // The start point is Uniform in the PHYSICAL state, so for the geometric
    // models (which solve in Y = log X) it is NOT uniform in the solver
    // coordinate -- its density there is e^y/(Zhi - Zlo).  Getting this wrong is
    // invisible over a narrow range and catastrophic over a wide one: it cost
    // 1.6e-2 on GBM (start 1 -> 1.5) and 2.3e-1 on Gompertz (start 1e-3 -> 0.5).
    // Taking the cell mass in the physical variable handles both exactly.
    const double L = m.length(0.0);
    const bool ls = m.bnd.log_state;
    const double Zlo = ls ? std::exp(z_lo) : z_lo;
    const double Zhi = ls ? std::exp(z_hi) : z_hi;
    const double span = Zhi - Zlo;
    const double xl0 = m.x_lo(0.0);
    for (int i = 0; i < M; ++i) {
      const double xa = xl0 + g.xf[i] * L;
      const double xb = xl0 + g.xf[i + 1] * L;
      const double Xa = ls ? std::exp(xa) : xa;
      const double Xb = ls ? std::exp(xb) : xb;
      const double ov = std::min(Xb, Zhi) - std::max(Xa, Zlo);
      if (ov > 0.0) q[i] = (ov / span) * g.rdx[i];  // cell mass / dx, q = L*p
    }
    return 0.0;
  }

  // ---- point start: short analytic warm-up ----
  const double z = z_hi;
  const double L0 = m.length(0.0);
  const double Bc = m.B();
  const double Ac = m.drift(z, 0.0);

  // Resolve the Gaussian by FPE_SEED_CELLS cells, and keep it clear of the
  // barrier.  t_seed wants to be as SMALL as possible -- everything before it is
  // answered by the frozen-coefficient approximation rather than by the solver,
  // and that is where the error lives (rel. error at t = 0.05 was 17% when
  // t_seed = 0.022).  Exact cell-average seeding below is what buys the licence
  // to shrink it: with midpoint sampling a 2-cell-wide Gaussian is badly
  // mis-integrated, with cell masses it is not.
  // The mesh is graded, so "FPE_SEED_CELLS cells wide" has to be measured with
  // the cell width where the start point actually sits, not with a nominal 1/M.
  //
  // It is capped at the width the UNIFORM mesh would have had, because grading
  // must never enlarge t_seed: everything before t_seed is answered by the
  // frozen-coefficient formula rather than by the solver, and a start point that
  // lands in the coarse far field (Gompertz, whose start range is 500x wide in
  // the physical state) would otherwise pay for the grading instead of gaining
  // from it.  Under-resolving the Gaussian is the cheaper error of the two --
  // the seed integrates exact cell masses, not midpoint samples.
  const double xi_z = (L0 > 0.0) ? (z - m.x_lo(0.0)) / L0 : 0.0;
  double s_target = FPE_SEED_CELLS * std::min(g.local_dx(xi_z), g.h) * L0;
  const double gap = m.x_hi(0.0) - z;
  if (gap > 0.0) s_target = std::min(s_target, 0.25 * gap);
  // With a second absorbing barrier the seed has to stay clear of BOTH, or the
  // frozen-coefficient Gaussian leaks mass through the lower one and the image
  // series below is truncated where it is not yet small.
  if constexpr (Model::lower_absorbing) {
    const double gap_lo = z - m.x_lo(0.0);
    if (gap_lo > 0.0) s_target = std::min(s_target, 0.25 * gap_lo);
  }
  double t_seed = (s_target / Bc) * (s_target / Bc);
  t_seed = std::min(t_seed, 0.25 * t_max);
  // The gaps above are measured at t = 0, but the seed is evaluated at t_seed,
  // so a boundary that collapses fast could invalidate them before the warm-up
  // is over.  Keeping t_seed well inside the collapse time scale makes the two
  // measurements the same to first order; t_seed is ~1 ms and tau is tens of ms
  // at the very least, so in practice this never binds.
  if (!m.bnd.fixed && m.bnd.tau > 0.0)
    t_seed = std::min(t_seed, 0.1 * m.bnd.tau);
  if (!(t_seed > 0.0)) t_seed = 1e-6;
  if (t_seed_force > 0.0) t_seed = t_seed_force;

  const double L = m.length(t_seed);
  const double a = m.x_hi(t_seed);
  const double sd = Bc * std::sqrt(t_seed);
  const double mean = z + Ac * t_seed;
  const double img_mean = 2.0 * a - z + Ac * t_seed;
  double img_w = std::exp(2.0 * Ac * (a - z) / (Bc * Bc));
  if (!std::isfinite(img_w)) img_w = 0.0;

  // Second image, reflected about the LOWER barrier.  Same construction as the
  // upper one -- reflect the start point about the barrier, then drift -- with
  // the sign of the distance reversed in the weight.  The exact two-barrier
  // solution is an infinite image series; both barriers are at least
  // 4*sd away by the t_seed cap above, so every image past this pair sits
  // >= 8 sd from the domain and contributes < 1e-14.
  const double xlo_s = m.x_lo(t_seed);
  double lo_mean = 0.0, lo_w = 0.0;
  if constexpr (Model::lower_absorbing) {
    lo_mean = 2.0 * xlo_s - z + Ac * t_seed;
    lo_w = std::exp(2.0 * Ac * (xlo_s - z) / (Bc * Bc));
    if (!std::isfinite(lo_w)) lo_w = 0.0;
  }

  // q_i = (cell mass) / dx_i, since q = L*p and the physical cell width is dx*L.
  for (int i = 0; i < M; ++i) {
    const double xa = xlo_s + g.xf[i] * L;
    const double xb = xlo_s + g.xf[i + 1] * L;
    double mass = fpe_norm_mass(xa, xb, mean, sd)
                  - img_w * fpe_norm_mass(xa, xb, img_mean, sd);
    if constexpr (Model::lower_absorbing) {
      mass -= lo_w * fpe_norm_mass(xa, xb, lo_mean, sd);
    }
    if (!(mass > 0.0)) mass = 0.0;
    q[i] = mass * g.rdx[i];
  }

  // The mass absorbed at the LOWER barrier before t_seed.  fpe_solve knows the
  // total absorbed exactly (1 - sum dx*q) but not how it splits, so hand it the
  // lower share and let the upper take the remainder -- the same division of
  // labour the time march uses.  This is not negligible: the t_seed cap leaves
  // the barrier 4 sd away, so the single-barrier hitting probability is ~3e-5,
  // which would otherwise be silently attributed to the wrong response.
  if constexpr (Model::lower_absorbing) {
    if (absorbed_lower != nullptr) {
      const double d = z - xlo_s;
      if (d > 0.0 && sd > 0.0) {
        const double u = (-Ac * t_seed - d) / (sd * M_SQRT2);
        const double w = (-d + Ac * t_seed) / (sd * M_SQRT2);
        double wt = std::exp(-2.0 * Ac * d / (Bc * Bc));
        if (!std::isfinite(wt)) wt = 0.0;
        double p = 0.5 * std::erfc(-u) + wt * 0.5 * std::erfc(-w);
        if (!(p > 0.0)) p = 0.0;
        *absorbed_lower = std::min(1.0, p);
      } else {
        *absorbed_lower = 0.0;
      }
    }
  }
  return t_seed;
}

// ---------------------------------------------------------------------------
// Convenience: seed, march, and report.
// ---------------------------------------------------------------------------
// Default mesh grading: ratio of far-field to barrier cell width.  See the
// sweep in the implementation log; 1.0 recovers the uniform mesh.
constexpr double FPE_GRADE = 8.0;

// ...but NOT for the two-boundary model, which wants a uniform mesh.  Grading
// pays off when most of the domain is empty far field; the bounded model has no
// far field at all -- the domain is exactly [0,a] and the density is O(1) across
// it -- so cells moved toward the barriers are taken from where the solution
// actually varies.  Measured (beta = 0 against the Navarro-Fuss DDM, nx = 384,
// nt = 3000, worst |dlog f| over the central mass of the response):
//
//   grade        1        4        8
//   sigma=1.00  1.22e-4  1.24e-4  1.91e-4
//   sigma=0.50  5.52e-3  8.18e-3  1.08e-2
//   sigma=0.25  1.69e-1  2.36e-1  3.08e-1
//
// Uniform wins at every setting and the gap widens as the domain gets wide in
// units of sigma.  The symmetric map is still worth having for a caller that
// asks for grading, but it is not the default here.
constexpr double FPE_GRADE_BOUNDED = 1.0;

// Default time grading.  1.0 (uniform) is kept for the validation entry points
// in fpe_diffusion.cpp, whose value lies in being directly comparable to the
// Volterra solver and the closed forms; a likelihood caller, which sets t_max
// from the largest RT in a data set and therefore has a large dt to spend,
// should pass FPE_TGRADE.  See FPE_TimeSchedule for the measurements.
constexpr double FPE_TGRADE = 32.0;

template <class Model>
inline FPE_Result fpe_run(const Model& m, double z_lo, double z_hi,
                          double t_max, int M, int nt,
                          double grade = FPE_GRADE, double tgrade = 1.0) {
  FPE_Mesh g;
  g.build(M, grade, Model::symmetric_mesh);
  std::vector<double> q0;
  double absorbed_lower = 0.0;
  const double t0 = fpe_seed(m, z_lo, z_hi, g, t_max, q0, &absorbed_lower);
  return fpe_solve(m, q0, t0, t_max, g, std::max(1, nt), tgrade, absorbed_lower);
}

} // namespace fpe

#endif // fpe_models_h

#ifndef fpe_models_h
#define fpe_models_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The four leaky-accumulation models, expressed for the Fokker-Planck solver in
// fpe_solver.h.  Each model supplies five things:
//
//   x_lo()             fixed lower edge of the computational domain
//   B()                diffusion coefficient (constant in x -- see below)
//   drift(x, t)        A(x,t)
//   length(t)          L(t) = a(t) - x_lo
//   length_prime(t)    L'(t) = a'(t)
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
// Collapsing boundary.  Same Weibull-style decay as exp_decay_scalar()
// (utils_reducible_diffusion.h:56), with the derivative in closed form.
// ---------------------------------------------------------------------------
struct FPE_Boundary {
  double b0 = 1.0, binf = 1.0, tau = 0.0, pw = 0.0;
  bool fixed = true;
  bool log_state = false;      // solve in Y = log X, so the barrier is log a(t)

  void set(double b0_, double binf_, double tau_, double pow_, bool log_state_) {
    b0 = std::abs(b0_);
    binf = std::abs(binf_);
    tau = tau_;
    pw = pow_;
    fixed = (std::abs(binf - b0) <= FPE_EPS) || tau <= 0.0 || pw <= 0.0;
    log_state = log_state_;
  }

  // physical boundary b(t) (before any log transform)
  double b(double t) const {
    if (fixed) return b0;
    const double amp = b0 - binf;
    if (!(t > 0.0)) return b0;
    const double s = std::exp(pw * (std::log(t) - std::log(tau)));
    return binf + amp * std::exp(-s);
  }

  // db/dt = -amp * (p/t) * s * exp(-s),  s = (t/tau)^p
  double b_prime(double t) const {
    if (fixed) return 0.0;
    const double amp = b0 - binf;
    const double tt = std::max(t, 1e-9);
    const double s = std::exp(pw * (std::log(tt) - std::log(tau)));
    if (!std::isfinite(s)) return 0.0;
    const double d = -amp * (pw / tt) * s * std::exp(-s);
    return std::isfinite(d) ? d : 0.0;
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
  double A = 0.0;          // mu, or mu - sigma^2/2 in log state
  double sigma = 1.0;
  double xlo = -1.0;
  FPE_Boundary bnd;

  double x_lo() const { return xlo; }
  double B() const { return sigma; }
  double drift(double /*x*/, double /*t*/) const { return A; }
  double length(double t) const { return bnd.a(t) - xlo; }
  double length_prime(double t) const { return bnd.a_prime(t); }

  // Atil(xi,t) = ( A - xi*L'(t) ) / L(t)
  void atil_affine(double /*t*/, double L, double Lp, double& a0, double& a1) const {
    a0 = A / L;
    a1 = -Lp / L;
  }
};

// ---------------------------------------------------------------------------
// Mean-reverting model: OU directly, Gompertz after Y = log X.
// ---------------------------------------------------------------------------
struct FPE_ModelOU {
  double lambda = 1.0, theta = 0.0, sigma = 1.0;
  double xlo = -1.0;
  FPE_Boundary bnd;

  double x_lo() const { return xlo; }
  double B() const { return sigma; }
  double drift(double x, double /*t*/) const { return -lambda * (x - theta); }
  double length(double t) const { return bnd.a(t) - xlo; }
  double length_prime(double t) const { return bnd.a_prime(t); }

  // Atil(xi,t) = ( -lambda*(x_lo + xi*L - theta) - xi*L'(t) ) / L(t)
  void atil_affine(double /*t*/, double L, double Lp, double& a0, double& a1) const {
    a0 = -lambda * (xlo - theta) / L;
    a1 = (-lambda * L - Lp) / L;
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

inline double fpe_x_lo_ou(double z_min, double lambda, double theta,
                          double sigma, double t_max) {
  const double var = (lambda > FPE_EPS)
    ? (1.0 - std::exp(-2.0 * lambda * t_max)) / (2.0 * lambda)
    : t_max;
  // the OU mean travels monotonically from z_min toward theta, so the lowest
  // mean ever attained is min(z_min, theta)
  return std::min(z_min, theta) - FPE_NSD * sigma * std::sqrt(var);
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
template <class Model>
inline double fpe_seed(const Model& m, double z_lo, double z_hi, int M,
                       double t_max, std::vector<double>& q) {
  const double h = 1.0 / M;
  q.assign(M, 0.0);

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
    for (int i = 0; i < M; ++i) {
      const double xa = m.x_lo() + (i * h) * L;
      const double xb = m.x_lo() + ((i + 1) * h) * L;
      const double Xa = ls ? std::exp(xa) : xa;
      const double Xb = ls ? std::exp(xb) : xb;
      const double ov = std::min(Xb, Zhi) - std::max(Xa, Zlo);
      if (ov > 0.0) q[i] = (ov / span) / h;   // cell mass / h, since q = L*p
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
  double s_target = FPE_SEED_CELLS * h * L0;
  const double gap = m.bnd.a(0.0) - z;
  if (gap > 0.0) s_target = std::min(s_target, 0.25 * gap);
  double t_seed = (s_target / Bc) * (s_target / Bc);
  t_seed = std::min(t_seed, 0.25 * t_max);
  if (!(t_seed > 0.0)) t_seed = 1e-6;

  const double L = m.length(t_seed);
  const double a = m.bnd.a(t_seed);
  const double sd = Bc * std::sqrt(t_seed);
  const double mean = z + Ac * t_seed;
  const double img_mean = 2.0 * a - z + Ac * t_seed;
  double img_w = std::exp(2.0 * Ac * (a - z) / (Bc * Bc));
  if (!std::isfinite(img_w)) img_w = 0.0;

  // q_i = (cell mass) / h, since q = L*p and the cell width is h*L.
  for (int i = 0; i < M; ++i) {
    const double xa = m.x_lo() + (i * h) * L;
    const double xb = m.x_lo() + ((i + 1) * h) * L;
    double mass = fpe_norm_mass(xa, xb, mean, sd)
                  - img_w * fpe_norm_mass(xa, xb, img_mean, sd);
    if (!(mass > 0.0)) mass = 0.0;
    q[i] = mass / h;
  }
  return t_seed;
}

// ---------------------------------------------------------------------------
// Convenience: seed, march, and report.
// ---------------------------------------------------------------------------
template <class Model>
inline FPE_Result fpe_run(const Model& m, double z_lo, double z_hi,
                          double t_max, int M, int nt) {
  std::vector<double> q0;
  const double t0 = fpe_seed(m, z_lo, z_hi, M, t_max, q0);
  return fpe_solve(m, q0, t0, t_max, M, std::max(1, nt));
}

} // namespace fpe

#endif // fpe_models_h

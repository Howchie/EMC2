#ifndef fpe_race_h
#define fpe_race_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Making a grid-valued PDE solver behave like a pointwise race kernel.
//
// EMC2's race contract is pointwise: every kernel is handed one row per
// (trial x accumulator) and asked for f(rt | par) or S(rt | par) at that row's
// own rt.  A Fokker-Planck march is the opposite shape -- it amortises O(M*nt)
// work over a WHOLE time grid and answers every rt at once, so evaluating it
// once per row would cost several hundred times more than necessary.
//
// The bridge is a solve cache keyed on the parameter tuple.  The number of
// solves then equals the number of DISTINCT parameter rows, not the number of
// trials: a forstmann-style design (v~lM, s~lM, B~E+lR, A~1, t0~1) has 12 of
// them covering 524 compressed rows, so a subject costs ~12 solves rather than
// ~524.  A trend model makes every trial distinct and degenerates to one solve
// per row -- correct, just slow; see the cost message in make_emc().
//
// The cache is also what makes CENSORING affordable.  The common unknown-winner
// path (log_surv_cm, particle_ll.cpp) needs scalar CDF values at the censoring
// and truncation bounds, on parameter rows the batch kernels have already
// solved.  Those calls become interpolations instead of solves.
//
// Everything here is model-agnostic except rou_solve(); adding the BM, GBM or
// Gompertz race means another solve function and another Key layout.
//
// NOTE ON t0: it is deliberately NOT part of the cache key.  t0 shifts the time
// axis (tt = rt - t0) and does not change the solve, so two design cells that
// differ only in t0 share one march.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include "fpe_models.h"

namespace fperace {

// Log of the smallest density/survivor worth representing.  Grids are stored in
// LOG space and interpolated there, because both quantities are close to
// exponential in t over a grid interval -- linear interpolation of log S is
// dramatically better than of S once S is small, and that is precisely the
// regime a race likelihood spends its time in.
constexpr double LOG_FLOOR = -700.0;

inline double safe_log(double x) {
  return (x > 0.0 && std::isfinite(x)) ? std::max(std::log(x), LOG_FLOOR)
                                       : LOG_FLOOR;
}

// ---------------------------------------------------------------------------
// Resolution.  dt_target rather than a fixed nt: t_max is set by the data (and
// is extended when a scalar/GSL call needs a longer horizon), so pinning nt
// would make dt -- and therefore the accuracy -- a function of the horizon, and
// a re-solve to a longer horizon would silently change the likelihood of rows
// that had already been evaluated.  Pinning dt keeps re-solves consistent.
// ---------------------------------------------------------------------------
// Defaults set by the step-10 sweep (race_ou_integration_plan.md section 8):
// nx is the binding constraint and dt is nearly free once the time grid is
// graded, so the budget goes to nx.  At (384, 2e-3) the likelihood carries
// 1.1e-3 nats/trial of discretisation bias against a converged reference, and a
// fit is within 2e-4 nats/trial of the reference fit.
struct FPE_Grid {
  int nx = 384;
  double dt_target = 2e-3;
  int nt_min = 256;
  int nt_max = 4096;
  double grade = fpe::FPE_GRADE;
  double tgrade = fpe::FPE_TGRADE;

  int nt_for(double t_max) const {
    const double want = std::ceil(t_max / std::max(dt_target, 1e-6));
    const int n = (want < static_cast<double>(nt_max)) ? static_cast<int>(want)
                                                       : nt_max;
    return std::max(n, nt_min);
  }
};

// ---------------------------------------------------------------------------
// Cache key.  Bit-exact comparison is the right test: these doubles come from
// the same ParamTable columns for every row in a design cell, so equal cells
// compare equal and unequal cells must not share a solve.
// ---------------------------------------------------------------------------
struct Key {
  double v = 0.0, k = 0.0, b = 0.0, A = 0.0;

  // Boundary.  For FPE_BND_FIXED the shape fields are forced to zero by
  // rou_key(), so every fixed-bound row compares equal on them and a model that
  // never collapses behaves exactly as it did before collapse existed.
  int bkind = fpe::FPE_BND_FIXED;
  double binf = 0.0, tau = 0.0, pw = 0.0;

  bool operator==(const Key& o) const {
    return v == o.v && k == o.k && b == o.b && A == o.A &&
           bkind == o.bkind && binf == o.binf && tau == o.tau && pw == o.pw;
  }
  bool finite() const {
    return std::isfinite(v) && std::isfinite(k) && std::isfinite(b) &&
           std::isfinite(A) && std::isfinite(binf) && std::isfinite(tau) &&
           std::isfinite(pw);
  }
};

// Collapse spec on the MODEL's scale, before s-scaling.  Binf is the asymptotic
// boundary measured from zero (unlike B, which is measured from the top of the
// start-point range).  The solver's FPE_Boundary knows no parameter names; this
// is the one place that maps them.
struct BndSpec {
  int kind = fpe::FPE_BND_FIXED;
  double Binf = 0.0, tau = 0.0, pw = 0.0;
};

// Build the key for one accumulator's parameters.  s is scaled out exactly as
// the RDM does (v/s, B/s, A/s, with k unchanged -- k has units of 1/time and is
// invariant under the state rescaling Y = X/s), which is what makes the k = 0
// member reproduce the analytic Wald parameterisation term for term.  Returns
// false for a row this model cannot represent, which the callers turn into a
// zero density / full survivor.
//
// It lives here rather than with the kernels so that the C++ likelihood path
// and the R-facing dROU/pROU cannot drift apart in how they form a key.
inline bool rou_key(double v, double k, double B, double A, double s,
                    const BndSpec& bs, Key& out) {
  if (!(s > 0.0) || !std::isfinite(s)) return false;
  const double inv_s = 1.0 / s;
  out.v = v * inv_s;
  out.k = (k > 0.0 && std::isfinite(k)) ? k : 0.0;
  const double AA = (A > 0.0 && std::isfinite(A)) ? A * inv_s : 0.0;
  out.A = AA;
  out.b = B * inv_s + AA;                  // b = B + A, both already scaled

  // The collapse asymptote is Binf itself, measured from ZERO -- not Binf + A.
  // The bound may therefore descend into the start-point range [0, A], which is
  // the point: a bound that meets the start point is a forced response, and
  // that is a state the model should be able to express.  Binf is bounded away
  // from 0 only so the solver's domain [x_lo, b(t)] cannot collapse.
  out.bkind = bs.kind;
  if (bs.kind == fpe::FPE_BND_FIXED) {
    out.binf = 0.0; out.tau = 0.0; out.pw = 0.0;   // canonical, so keys compare
  } else {
    if (!(bs.Binf >= 0.0) || !std::isfinite(bs.Binf) ||
        !(bs.tau > 0.0) || !std::isfinite(bs.tau)) return false;
    out.binf = bs.Binf * inv_s;
    out.tau = bs.tau;
    // Only the Weibull form has a shape exponent; zero it otherwise so that two
    // rows differing in an unused column still share one solve.
    out.pw = (bs.kind == fpe::FPE_BND_WEIBULL) ? bs.pw : 0.0;
    if (bs.kind == fpe::FPE_BND_WEIBULL &&
        (!(bs.pw > 0.0) || !std::isfinite(bs.pw))) return false;
  }
  return out.finite() && out.b > 0.0 && out.v > 0.0;
}

// Fixed-boundary overload -- the common case, and every pre-collapse call site.
inline bool rou_key(double v, double k, double B, double A, double s, Key& out) {
  static const BndSpec fixed;
  return rou_key(v, k, B, A, s, fixed, out);
}

struct Entry {
  Key key;
  double t_max = 0.0;
  std::vector<double> t, log_pdf, log_S;
};

struct SolveCache {
  FPE_Grid grid;
  // Boundary form for this model instance, set once when the adapter is built.
  // It lives here rather than on ContextForRaceModels because it is a property
  // of the solve, and because every path that can reach a solve already holds
  // the cache.
  int bnd_kind = fpe::FPE_BND_FIXED;
  std::vector<Entry> e;
  std::vector<int> row_group;   // scratch: row -> index into e, or -1

  // Cleared once per particle.  Keys are exact, so a stale entry could never be
  // returned for the wrong parameters; the clear exists to bound memory, since
  // a run visits thousands of particles.
  void new_particle() { e.clear(); }
};

// ---------------------------------------------------------------------------
// One racing-OU march.
//
//   dX = (v - k X) dt + s dW,   X_0 ~ U(0, A),   absorbed at b = B + A.
//
// s is scaled out exactly as the RDM does (Y = X/s carries v/s, b/s, A/s and
// the same k), so the k = 0 member reproduces the analytic Wald parameterisation
// term for term.
// ---------------------------------------------------------------------------
inline void rou_solve(const Key& p, double t_max, const FPE_Grid& gr, Entry& out) {
  fpe::FPE_ModelOU m;
  m.v = p.v;
  m.lambda = p.k;
  m.sigma = 1.0;                                  // s already divided out
  m.bnd.set_kind(p.bkind, p.b,
                 (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                 p.tau, p.pw, false);
  m.xlo = fpe::fpe_x_lo_ou(0.0, p.v, p.k, 1.0, t_max);

  const fpe::FPE_Result r =
      fpe::fpe_run(m, 0.0, p.A, t_max, gr.nx, gr.nt_for(t_max), gr.grade,
                   gr.tgrade);

  const size_t n = r.t.size();
  out.key = p;
  out.t_max = t_max;
  out.t = r.t;
  out.log_pdf.resize(n);
  out.log_S.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.log_pdf[i] = safe_log(r.pdf[i]);
    out.log_S[i] = safe_log(r.surv[i]);
  }
}

// Find the entry for `p` covering at least t_need, solving or extending as
// required.  Returns an INDEX, not a pointer: solving can reallocate `e`.
inline int cache_get(SolveCache& C, const Key& p, double t_need) {
  if (!(t_need > 0.0)) t_need = 1e-3;
  for (size_t i = 0; i < C.e.size(); ++i) {
    if (!(C.e[i].key == p)) continue;
    if (C.e[i].t_max >= t_need) return static_cast<int>(i);
    // Extend in place.  dt is pinned by FPE_Grid, so the values on the shared
    // part of the horizon move only at the discretisation-error level.
    rou_solve(p, t_need, C.grid, C.e[i]);
    return static_cast<int>(i);
  }
  C.e.push_back(Entry());
  rou_solve(p, t_need, C.grid, C.e.back());
  return static_cast<int>(C.e.size()) - 1;
}

// ---------------------------------------------------------------------------
// Interpolation, in log space.
//
// Below the grid the density is zero (the march is seeded at t_seed > 0, so the
// grid does not reach back to 0) but the survivor is not -- cache_get's grid
// carries the mass absorbed before t_seed.  This is the same asymmetry
// fpe::grid_lookup handles with zero_below_grid.
// ---------------------------------------------------------------------------
inline double interp_log(const std::vector<double>& tg,
                         const std::vector<double>& vg, double t,
                         bool floor_below) {
  const size_t n = tg.size();
  if (n == 0) return LOG_FLOOR;
  if (t <= tg.front()) return floor_below ? LOG_FLOOR : vg.front();
  if (t >= tg.back()) return vg.back();
  const size_t j = static_cast<size_t>(
      std::lower_bound(tg.begin(), tg.end(), t) - tg.begin());
  const double lo = vg[j - 1], hi = vg[j];
  // A floored neighbour carries no information; do not average garbage into a
  // live value, take the live one.
  if (lo <= LOG_FLOOR) return hi;
  if (hi <= LOG_FLOOR) return lo;
  const double w = (t - tg[j - 1]) / (tg[j] - tg[j - 1]);
  return lo + w * (hi - lo);
}

inline double entry_log_pdf(const Entry& en, double tt) {
  return interp_log(en.t, en.log_pdf, tt, true);
}
inline double entry_log_S(const Entry& en, double tt) {
  return interp_log(en.t, en.log_S, tt, false);
}

// ---------------------------------------------------------------------------
// Reference simulator: dX = (v - k X) dt + s dW, X_0 ~ U(0,A), absorbed at b.
//
// Written here rather than reusing simulate_ou_hit_times_bb() because that one
// is parameterised by (lambda, theta) and rejects lambda <= 0, so it cannot
// reach the k = 0 member at all -- which is the one case with an independent
// analytic answer to check against.
//
// Exact OU transition per step (no Euler discretisation of the drift), with the
// k -> 0 limits taken through expm1 so k = 0 is a regular value:
//   phi   = exp(-k dt)
//   mean  = X phi + v * (1 - phi)/k        -> X + v dt
//   var   = s^2 (1 - phi^2)/(2k)           -> s^2 dt
// plus the Brownian-bridge correction for a crossing that happens strictly
// inside a step, which is what stops the hit times being biased late:
//   P(cross | X_t, X_{t+dt}) = exp(-2 (b - X_t)(b - X_{t+dt}) / (s^2 dt))
// ---------------------------------------------------------------------------
// A COLLAPSING bound needs nothing new here beyond evaluating b at both ends of
// the step.  Over one step the barrier is treated as linear, and the crossing
// probability for a Brownian bridge against a linear barrier is the same
// expression with the two endpoint gaps -- so the correction stays exact to the
// order the rest of the step already is, rather than degrading to "compare
// against b at one end", which biases hit times late exactly when the bound is
// moving fastest.
template <class Unif01, class Norm01>
inline double rou_hit_time_bnd(double v, double k, double s, double A,
                               const fpe::FPE_Boundary& bnd, double dt,
                               double t_max, Unif01 runif, Norm01 rnorm) {
  double X = (A > 0.0) ? A * runif() : 0.0;
  double b = bnd.b(0.0);
  if (X >= b) return 0.0;

  const double phi = std::exp(-k * dt);
  const double m1 = -std::expm1(-k * dt);            // 1 - phi
  const double m2 = -std::expm1(-2.0 * k * dt);      // 1 - phi^2
  const double drift_gain = (k > 1e-10) ? (m1 / k) : dt;
  const double var = (k > 1e-10) ? (s * s * m2 / (2.0 * k)) : (s * s * dt);
  const double sd = std::sqrt(std::max(var, 0.0));
  const double inv_2var_bb = 2.0 / (s * s * dt);     // local BB uses s^2 dt

  double t = 0.0;
  while (t < t_max) {
    const double X1 = X * phi + v * drift_gain + sd * rnorm();
    const double b1 = bnd.fixed ? b : bnd.b(t + dt);
    t += dt;
    if (X1 >= b1) {
      // Where the path and the barrier, both linear over the step, meet.
      const double d0 = b - X, d1 = b1 - X1;
      const double frac = d0 / std::max(d0 - d1, 1e-300);
      return t - dt + std::min(std::max(frac, 0.0), 1.0) * dt;
    }
    // Did the path cross and come back inside the step?
    const double pc = std::exp(-(b - X) * (b1 - X1) * inv_2var_bb);
    if (runif() < pc) return t - dt + runif() * dt;
    X = X1;
    b = b1;
  }
  return std::numeric_limits<double>::infinity();  // never finished in horizon
}

// Fixed-boundary overload, for the k = 0 / k > 0 cross-validation.
template <class Unif01, class Norm01>
inline double rou_hit_time(double v, double k, double s, double A, double b,
                           double dt, double t_max, Unif01 runif, Norm01 rnorm) {
  fpe::FPE_Boundary bnd;
  bnd.set_kind(fpe::FPE_BND_FIXED, b, b, 0.0, 0.0, false);
  return rou_hit_time_bnd(v, k, s, A, bnd, dt, t_max, runif, rnorm);
}

} // namespace fperace

#endif // fpe_race_h

#ifndef fpe_solver_h
#define fpe_solver_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Model-agnostic Fokker-Planck first-passage-time solver.
//
// Solves the forward Kolmogorov equation for the sub-density of paths that have
// not yet been absorbed at a (possibly moving) upper boundary a(t), and reads the
// first-passage density off as the probability flux through that boundary.
//
//   dp/dt = -d/dx [ A(x,t) p ] + (1/2) d2/dx2 [ B^2 p ],   p(a(t),t) = 0
//
// Cost is O(M * N_t) -- LINEAR in max_rt -- as against the O(N_t^2) history sum
// of the Volterra formulation in utils_reducible_diffusion.h.
//
// Design notes
// ------------
// * The model is a TEMPLATE PARAMETER, never a std::function.  The Volterra path
//   pays ~330 ns per kernel evaluation almost entirely because KernelFn is an
//   un-inlinable indirect call inside its innermost loop; do not repeat that.
// * Cell-centred FINITE VOLUME, so the FPT density *is* the numerical face flux
//   at xi = 1.  No separate one-sided derivative of the solution is needed, and
//   the cdf comes free from the mass that has left the domain.
// * Interior faces use the Scharfetter-Gummel exponentially-fitted flux, which
//   degrades gracefully to upwinding at high Peclet number, so small sigma (or a
//   fast-collapsing bound) does not produce oscillations.
// * Crank-Nicolson with a Rannacher (backward-Euler) start, since CN rings on the
//   near-singular initial data of a point start.
//
// No Rcpp, no model knowledge here.  See fpe_models.h for the four models.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace fpe {

// Bernoulli function Bern(z) = z / (exp(z) - 1), with Bern(0) = 1.
inline double bern(double z) {
  if (std::abs(z) < 1e-8) return 1.0 - 0.5 * z;   // series, avoids 0/0
  if (z > 700.0)  return 0.0;                     // exp overflows; limit is 0
  if (z < -700.0) return -z;                      // expm1(z) -> -1
  return z / std::expm1(z);
}

// Same, but reusing an exp(z) that the caller already has.  Used by the fast
// path in build_op(), which generates exp(P_j) along a face row by recurrence
// instead of calling expm1() once per face.  The |z| > 1e-3 guard keeps the
// cancellation in (e - 1) below ~1e-13 relative; smaller z takes the series.
inline double bern_from_exp(double z, double e) {
  if (std::abs(z) < 1e-3) return 1.0 - 0.5 * z + z * z / 12.0;
  if (!(e > 0.0)) return -z;                      // underflowed: exp(z) -> 0
  if (!std::isfinite(e)) return 0.0;              // overflowed
  return z / (e - 1.0);
}

// Bern(-z) = Bern(z) + z.
//   z/(e^z - 1) + z = z*e^z/(e^z - 1) = -(-z)/(e^{-z} - 1)
// so the upwind and downwind face weights cost one transcendental between them,
// not two.
inline double bern_neg(double bz, double z) { return bz + z; }

// Tridiagonal system (Thomas algorithm), split into factor and solve.
//   sub[i] * x[i-1] + dia[i] * x[i] + sup[i] * x[i+1] = rhs[i]
// sub[0] and sup[n-1] are ignored.
//
// The split exists because the left-hand side is CONSTANT for the whole march
// whenever the boundary is fixed (see FPE_Op / static_op below), so the O(M)
// elimination can be done once instead of once per step.  It also replaces the
// per-step divisions with multiplies by a stored reciprocal.
struct FPE_Tri {
  std::vector<double> sub, dia, sup, cprime, dinv;

  void resize(int n) {
    if (static_cast<int>(dia.size()) != n) {
      sub.resize(n); dia.resize(n); sup.resize(n);
      cprime.resize(n); dinv.resize(n);
    }
  }

  void factor() {
    const int n = static_cast<int>(dia.size());
    if (n == 0) return;
    double denom = dia[0];
    dinv[0]   = 1.0 / denom;
    cprime[0] = sup[0] * dinv[0];
    for (int i = 1; i < n; ++i) {
      denom = dia[i] - sub[i] * cprime[i - 1];
      if (denom == 0.0) denom = std::numeric_limits<double>::min();
      dinv[i]   = 1.0 / denom;
      cprime[i] = (i + 1 < n) ? sup[i] * dinv[i] : 0.0;
    }
  }

  void solve(const std::vector<double>& rhs, std::vector<double>& out) const {
    const int n = static_cast<int>(dia.size());
    if (n == 0) return;
    out.resize(n);
    out[0] = rhs[0] * dinv[0];
    for (int i = 1; i < n; ++i) {
      out[i] = (rhs[i] - sub[i] * out[i - 1]) * dinv[i];
    }
    for (int i = n - 2; i >= 0; --i) {
      out[i] -= cprime[i] * out[i + 1];
    }
  }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Graded mesh on the normalised domain xi in [0,1].
//
// Cells are uniform in a stretched coordinate eta and clustered toward xi = 1:
//
//   xi(eta) = 1 - sinh(c (1 - eta)) / sinh(c)
//
// `grade` is the ratio of the far-field cell width to the barrier cell width,
// and c = acosh(grade).  grade = 1 gives back the uniform mesh exactly.
//
// This matters because essentially all of the accuracy lives in the last few
// percent of the domain: the absorbing barrier is at xi = 1 and (once the domain
// is sized by fpe_x_lo_*) the start point sits just below it, while the lower
// 80-90% of the domain holds the far tail and almost no mass.  A uniform mesh
// spends its cells in proportion to length rather than to where the solution
// actually varies.
//
// Finite volume takes a non-uniform mesh natively -- the balance
// dq_i/dt = -(F_{i+1} - F_i)/dx_i is exact for any cell widths -- but every
// stencil coefficient becomes cell-dependent, so all the geometry that used to
// be the single scalar h is precomputed here, once per solve, and the time march
// only ever touches D(t) and the two affine drift coefficients.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct FPE_Mesh {
  int M = 0;
  bool uniform = true;
  std::vector<double> xf;    // M+1 face positions, xf[0] = 0, xf[M] = 1
  std::vector<double> xc;    // M cell centres
  std::vector<double> dx;    // M cell widths
  std::vector<double> rdx;   // 1 / dx
  std::vector<double> rdc;   // 1 / (xc[j] - xc[j-1]), interior faces j = 1..M-1
  std::vector<double> pu;    // xc[j] - xc[j-1]        )  P_j = (a0*pu + a1*pv)/D
  std::vector<double> pv;    // xf[j] * (that spacing) )
  double cA = 0.0, cB = 0.0; // absorbing face: F_M = D*(cA*q[M-1] - cB*q[M-2])
  double h = 0.0;            // uniform mesh only

  void build(int M_, double grade) {
    M = std::max(M_, 2);
    h = 1.0 / M;
    uniform = !(grade > 1.0 + 1e-9);

    xf.resize(M + 1); xc.resize(M); dx.resize(M); rdx.resize(M);
    rdc.assign(M + 1, 0.0); pu.assign(M + 1, 0.0); pv.assign(M + 1, 0.0);

    if (uniform) {
      for (int j = 0; j <= M; ++j) xf[j] = j * h;
    } else {
      const double c  = std::acosh(grade);
      const double sh = std::sinh(c);
      for (int j = 0; j <= M; ++j) {
        xf[j] = 1.0 - std::sinh(c * (1.0 - j * h)) / sh;
      }
    }
    xf[0] = 0.0; xf[M] = 1.0;

    for (int i = 0; i < M; ++i) {
      dx[i]  = xf[i + 1] - xf[i];
      rdx[i] = 1.0 / dx[i];
      xc[i]  = 0.5 * (xf[i] + xf[i + 1]);
    }
    for (int j = 1; j < M; ++j) {
      const double d = xc[j] - xc[j - 1];
      rdc[j] = 1.0 / d;
      pu[j]  = d;
      pv[j]  = xf[j] * d;
    }

    // Absorbing face at xi = 1.  q(1) = 0, so fit u(s) = c1 s + c2 s^2 through
    // the two nearest centres, s = 1 - xi, and take F_M = D * u'(0) = D * c1:
    //   c1 = (q1 s2^2 - q2 s1^2) / (s1 s2 (s2 - s1))
    // which for a uniform mesh (s1 = h/2, s2 = 3h/2) collapses to the familiar
    // (9 q[M-1] - q[M-2]) / (3h).
    const double s1 = 1.0 - xc[M - 1];
    const double s2 = 1.0 - xc[M - 2];
    cA = s2 / (s1 * (s2 - s1));
    cB = s1 / (s2 * (s2 - s1));
  }

  // Width of the cell containing xi -- used to size the analytic seed.
  double local_dx(double xi) const {
    const int j = static_cast<int>(
      std::upper_bound(xf.begin(), xf.end(), xi) - xf.begin()) - 1;
    return dx[std::min(std::max(j, 0), M - 1)];
  }
};

// The spatial operator L in dq/dt = L q, as a tridiagonal matrix, plus the
// diffusivity D(t) needed to turn the solution into a boundary flux.
struct FPE_Op {
  std::vector<double> sub, diag, sup;
  double D = 0.0;

  // Every entry is overwritten by build_op(), so only grow -- never re-zero.
  // build_op() runs once per time step, so a gratuitous 3*M zero-fill there is
  // a measurable fraction of the whole solve.
  void resize(int M) {
    if (static_cast<int>(diag.size()) != M) {
      sub.resize(M); diag.resize(M); sup.resize(M);
    }
  }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Build L at time t.
//
// The moving boundary is normalised away with L(t) = a(t) - x_lo,
// xi = (x - x_lo)/L(t) in [0,1], q(xi,t) = L(t) p(x,t) (so mass is preserved):
//
//   Atil(xi,t) = ( A(x_lo + xi*L, t) - xi * L'(t) ) / L(t)
//   Btil(t)    = B / L(t)
//   D(t)       = 0.5 * Btil^2
//
// D is INDEPENDENT OF xi for all four supported models (B is constant in x after
// the log-state transform for GBM/Gompertz), so the flux is exactly
// F = Atil*q - D*dq/dxi with no d(D)/dxi correction term.  That is the reason
// this scheme stays a plain three-point stencil.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// On a UNIFORM mesh every supported model has Atil affine in xi, so the face
// Peclet numbers
//   P_j = Atil(j*h, t) * h / D = P1 + dP*(j-1)
// form an arithmetic sequence and exp(P_j) can be generated by one multiply per
// face.  That removes the two expm1() calls per face, which measured as the
// single dominant cost of the whole solve (13.3 ms -> see the implementation log
// in fpe_fht_solver_plan.md).  The recurrence drifts by only ~M*eps; it is
// abandoned for the safe path if the row spans an implausible Peclet range.
//
// A graded mesh breaks the arithmetic progression (the face spacing varies), so
// it pays one exp() per face.  That is the price of the grading and it is a good
// trade -- see the log.
template <class Model>
inline void build_op(const Model& m, double t, const FPE_Mesh& g, FPE_Op& op) {
  const int M      = g.M;
  const double Lb  = m.length(t);       // a(t) - x_lo
  const double Lp  = m.length_prime(t);
  const double Bt  = m.B() / Lb;
  const double D   = 0.5 * Bt * Bt;
  const double invD = 1.0 / D;
  op.D = D;
  op.resize(M);

  double a0, a1;                        // Atil(xi,t) = a0 + a1*xi
  m.atil_affine(t, Lb, Lp, a0, a1);

  const double P1 = (a0 + a1 * g.h) * g.h * invD;   // uniform mesh, face j = 1
  const double dP = a1 * g.h * g.h * invD;
  const bool uni = g.uniform && std::isfinite(P1) && std::isfinite(dP)
                   && std::abs(dP) * M < 20.0;
  double Pu = P1;
  double e  = uni ? std::exp(P1) : 0.0;
  const double er = uni ? std::exp(dP) : 0.0;

  // Face j sits at xi = xf[j], j = 1..M-1, between cells j-1 and j:
  //   F_j = (D/dc_j) * ( bern(-P_j) * q[j-1] - bern(P_j) * q[j] ),  P_j = Atil*dc_j/D
  // wm = bern(-P_j) multiplies q[j-1]; wp = bern(P_j) multiplies q[j].
  // Row i needs face i (computed on the previous iteration) and face i+1, so a
  // single forward sweep carries the previous face rather than caching arrays.
  // The same face weights enter the two adjacent rows divided by DIFFERENT cell
  // widths, so the carried quantity is D/dc_j and 1/dx_i is applied per row.
  double wm_prev = 0.0, wp_prev = 0.0, gf_prev = 0.0;

  // dq_i/dt = -( F_{i+1} - F_i ) / dx_i
  for (int i = 0; i < M; ++i) {
    double a_sub = 0.0, a_dia = 0.0, a_sup = 0.0;
    const double rdx = g.rdx[i];

    if (i > 0) {                      // inflow face F_i is an interior face
      a_sub += gf_prev * wm_prev * rdx;
      a_dia -= gf_prev * wp_prev * rdx;
    }
    // else: F_0 = 0, the no-flux far-field face.  This is what makes
    //       CDF = 1 - sum(dx*q) exact.

    if (i < M - 1) {                  // outflow face F_{i+1} is an interior face
      const int j = i + 1;
      double P, wp;
      if (uni) {
        P  = Pu;
        wp = bern_from_exp(P, e);
        Pu += dP;
        e  *= er;
      } else {
        P  = (a0 * g.pu[j] + a1 * g.pv[j]) * invD;
        wp = bern(P);
      }
      const double wm = bern_neg(wp, P);
      const double gf = D * g.rdc[j];
      a_dia -= gf * wm * rdx;
      a_sup += gf * wp * rdx;
      wm_prev = wm;
      wp_prev = wp;
      gf_prev = gf;
    } else {
      // Absorbing face at xi = 1.  q(1) = 0, so advection carries nothing and
      // the flux is purely diffusive; cA/cB are the one-sided weights built by
      // FPE_Mesh::build().
      a_dia -= D * g.cA * rdx;
      a_sub += D * g.cB * rdx;
    }

    op.sub[i]  = a_sub;
    op.diag[i] = a_dia;
    op.sup[i]  = a_sup;
  }
}

// Absorbing-face flux = the first-passage density g(t).
inline double flux_out(const FPE_Op& op, const std::vector<double>& q,
                       const FPE_Mesh& g) {
  const int M = g.M;
  if (M < 2) return 0.0;
  return op.D * (g.cA * q[M - 1] - g.cB * q[M - 2]);
}

// y = (I + c*L) q
inline void apply_shifted(const FPE_Op& op, double c,
                          const std::vector<double>& q, std::vector<double>& y) {
  const int M = static_cast<int>(q.size());
  y.resize(M);
  for (int i = 0; i < M; ++i) {
    double v = q[i] + c * op.diag[i] * q[i];
    if (i > 0)     v += c * op.sub[i] * q[i - 1];
    if (i < M - 1) v += c * op.sup[i] * q[i + 1];
    y[i] = v;
  }
}

struct FPE_Result {
  std::vector<double> t, pdf, cdf;
  // Survivor, reported SEPARATELY from cdf rather than left to the caller as
  // 1 - cdf.  S = sum(dx*q) is a sum of positive quantities, so it keeps full
  // RELATIVE accuracy far into the tail (measured |dlog S| <= 1.4e-2 down to
  // S = 7.5e-14), whereas 1 - cdf loses every significant digit once cdf -> 1.
  // A race likelihood multiplies loser survivors and lives or dies on that.
  std::vector<double> surv;
  double flux_mass_mismatch = 0.0;   // instability detector, NOT an error bound
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Time-step schedule.
//
// A uniform dt = T/nt spends its steps in proportion to elapsed time, but the
// solution's time variation is concentrated at the START: the first-passage
// density climbs from ~0 to its mode within the first few tens of ms and then
// decays smoothly over the remaining second or more.  With T set by the largest
// RT in a data set, a uniform grid leaves the rising flank badly under-resolved,
// and the error there is a pure function of dt (measured, OU lambda=4 theta=2
// b0=1, M=256; relative pdf error at t = 0.05 s):
//
//   dt      9.8ms   4.9ms   2.4ms   1.2ms   0.61ms
//   err      83%     26%      8%     2.6%    1.3%
//
// -- identically at T = 0.25 and T = 2.5, i.e. it tracks dt alone and not the
// number of steps.  That is the same argument that grades the SPATIAL mesh
// toward the barrier (see FPE_Mesh), applied to the time axis.
//
// Rather than a smoothly stretched grid, dt is held constant within blocks and
// DOUBLES between them.  Crank-Nicolson is second order for any step sequence,
// but a fixed boundary makes the left-hand side depend on dt alone (see
// static_op / set_lhs below), so a piecewise-constant dt keeps the one-time
// O(M) factorisation -- it is redone once per block, i.e. log2(tgrade) times for
// the whole march, instead of once per step.
//
// `tgrade` is the ratio of the last block's step to the first's; tgrade <= 1
// recovers the uniform grid exactly.  `nt` remains the total step budget.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct FPE_TimeSchedule {
  std::vector<double> dt;      // step size, per block
  std::vector<int> steps;      // step count, per block
};

inline FPE_TimeSchedule fpe_time_schedule(double T, int nt, double tgrade) {
  FPE_TimeSchedule s;
  const int n_steps = std::max(nt, 1);
  if (!(tgrade > 1.0 + 1e-9) || n_steps < 4) {
    s.dt.push_back(T / n_steps);
    s.steps.push_back(n_steps);
    return s;
  }
  int n_blk = 1 + static_cast<int>(std::floor(std::log(tgrade) / std::log(2.0) + 1e-9));
  n_blk = std::max(1, std::min(n_blk, n_steps / 2));
  const int per = std::max(2, n_steps / n_blk);
  // sum_b per * dt0 * 2^b = per * dt0 * (2^n_blk - 1) = T
  const double dt0 = T / (per * (std::ldexp(1.0, n_blk) - 1.0));
  s.dt.reserve(n_blk);
  s.steps.reserve(n_blk);
  for (int b = 0; b < n_blk; ++b) {
    s.dt.push_back(std::ldexp(dt0, b));
    s.steps.push_back(per);
  }
  return s;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// March q from t0 to t_max.
//
//   q0        cell-averaged sub-density at t0 on the normalised grid, already
//             carrying whatever mass was absorbed before t0 (see the seeding in
//             fpe_models.h) -- i.e. sum(dx*q0) = 1 - CDF(t0).
//   nt        total step budget.  The first two steps are replaced by four
//             backward-Euler half-steps (Rannacher), which is what keeps CN from
//             ringing on a point start.
//   tgrade    ratio of the last time step to the first (see FPE_TimeSchedule);
//             tgrade <= 1 gives the uniform grid dt = (t_max-t0)/nt.
//
// Two independent routes to the cdf are computed:
//   mass flavour:  CDF = 1 - sum(dx*q)         (exact, given F_0 = 0)
//   flux flavour:  CDF = CDF(t0) + trapz(g)
// Their maximum absolute difference is returned as flux_mass_mismatch.
//
// READ THIS BEFORE TRUSTING flux_mass_mismatch.  It is an INSTABILITY detector,
// not a convergence diagnostic, and it must never be used as a runtime
// accept/reject gate on accuracy.  Both routes are functionals of the same q, so
// any discretisation error common to the two cancels out of their difference.
// Measured (OU lambda=4 theta=2 sigma=1 b0=1, M=256, nt=512, uniform grid):
// mismatch reports 4.3e-8 while the true max CDF error against a converged
// reference is 1.4e-3 -- five orders of magnitude apart.  What the mismatch does
// catch is a march that has gone unstable or lost mass through a mis-built
// operator, which is worth having; it says nothing about resolution.
// Convergence must be established by refining (M, nt) and comparing solves.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class Model>
inline FPE_Result fpe_solve(const Model& m, const std::vector<double>& q0,
                            double t0, double t_max, const FPE_Mesh& g, int nt,
                            double tgrade = 1.0) {
  FPE_Result res;
  const int M = g.M;
  const FPE_TimeSchedule sched = fpe_time_schedule(t_max - t0, nt, tgrade);

  std::vector<double> q = q0, rhs;
  FPE_Tri tri;
  tri.resize(M);
  FPE_Op op_a, op_b;
  FPE_Op* op_old = &op_a;               // swapped, never copied
  FPE_Op* op_new = &op_b;

  // With a fixed boundary L(t), L'(t), D and the affine drift coefficients are
  // all constant, so the spatial operator never changes -- and since a Rannacher
  // half-step (c = step = dt/2) and a Crank-Nicolson full step (c = step/2 =
  // dt/2) use the SAME c, the left-hand side matrix is constant for as long as
  // dt is.  Factorise it once per BLOCK of the time schedule (once for the whole
  // march on a uniform grid; log2(tgrade) times on a graded one).
  const bool stat = m.static_op();

  auto set_lhs = [&](double c, const FPE_Op& o) {
    for (int i = 0; i < M; ++i) {
      tri.sub[i] = -c * o.sub[i];
      tri.dia[i] = 1.0 - c * o.diag[i];
      tri.sup[i] = -c * o.sup[i];
    }
    tri.factor();
  };

  int n_steps = 0;
  for (size_t b = 0; b < sched.steps.size(); ++b) n_steps += sched.steps[b];
  const int n_rann = std::min(2, n_steps);   // full steps replaced by BE halves
  const int n_lev  = 1 + 2 * n_rann + (n_steps - n_rann);
  res.t.reserve(n_lev);
  res.pdf.reserve(n_lev);
  res.cdf.reserve(n_lev);
  res.surv.reserve(n_lev);

  double t = t0;
  build_op(m, t, g, *op_old);

  double mass = 0.0;
  for (int i = 0; i < M; ++i) mass += g.dx[i] * q[i];

  // Clamp the reported cdf to [0,1] and keep it non-decreasing.  Once the
  // sub-density has decayed to the far tail, sum(dx*q) is a sum of quantities at
  // the 1e-17 level and 1 - mass is not reliably inside [0,1] nor monotone; a
  // consumer forming log(1 - cdf) or diff(cdf) must not have to defend against
  // that.  The RAW value is retained for the mismatch diagnostic so that
  // clamping cannot hide an unstable march.
  const double cdf0 = std::min(1.0, std::max(0.0, 1.0 - mass));
  double cdf_prev = cdf0;
  double cdf_flux = cdf0;
  double surv_prev = std::min(1.0, std::max(0.0, mass));
  // The absorbing-face flux is a one-sided difference of two nearly equal cell
  // values, so in the far tail it goes slightly NEGATIVE (measured -3e-16 where
  // the true density is 1.5e-17).  log() of that is NaN, which would silently
  // poison a likelihood; floor it at the source.
  double g_prev = flux_out(*op_old, q, g);
  if (!(g_prev > 0.0)) g_prev = 0.0;

  res.t.push_back(t);
  res.pdf.push_back(g_prev);
  res.cdf.push_back(cdf0);
  res.surv.push_back(surv_prev);
  res.flux_mass_mismatch = 0.0;

  // step_kind: true = backward Euler, false = Crank-Nicolson
  auto do_step = [&](double step, bool backward_euler) {
    const double t_new = t + step;
    if (!stat) build_op(m, t_new, g, *op_new);

    if (backward_euler) {
      rhs = q;
    } else {
      apply_shifted(*op_old, 0.5 * step, q, rhs);
    }

    // (I - c*L^{n+1}) q^{n+1} = rhs
    if (!stat) set_lhs(backward_euler ? step : 0.5 * step, *op_new);
    tri.solve(rhs, q);

    t = t_new;
    if (!stat) std::swap(op_old, op_new);

    double s = 0.0;
    for (int i = 0; i < M; ++i) s += g.dx[i] * q[i];
    const double cdf_raw = 1.0 - s;
    double cdf_mass = std::min(1.0, std::max(0.0, cdf_raw));
    if (cdf_mass < cdf_prev) cdf_mass = cdf_prev;
    cdf_prev = cdf_mass;

    // Survivor from the mass route directly, clamped and forced non-increasing.
    // NOT 1 - cdf_mass: that would throw away exactly the relative accuracy this
    // vector exists to preserve.
    double s_keep = std::min(1.0, std::max(0.0, s));
    if (s_keep > surv_prev) s_keep = surv_prev;
    surv_prev = s_keep;

    double gt = flux_out(*op_old, q, g);
    if (!(gt > 0.0)) gt = 0.0;
    cdf_flux += 0.5 * step * (g_prev + gt);
    g_prev = gt;

    res.t.push_back(t);
    res.pdf.push_back(gt);
    res.cdf.push_back(cdf_mass);
    res.surv.push_back(s_keep);
    res.flux_mass_mismatch =
      std::max(res.flux_mass_mismatch, std::abs(cdf_raw - cdf_flux));
  };

  int done = 0;
  for (size_t b = 0; b < sched.dt.size(); ++b) {
    const double dtb = sched.dt[b];
    // One factorisation per block: a Rannacher half-step (c = step = dt/2) and a
    // CN full step (c = step/2 = dt/2) share the same c, so the whole block runs
    // off a single elimination.
    if (stat) set_lhs(0.5 * dtb, *op_old);
    for (int k = 0; k < sched.steps[b]; ++k, ++done) {
      if (done < n_rann) {
        do_step(0.5 * dtb, true);
        do_step(0.5 * dtb, true);
      } else {
        do_step(dtb, false);
      }
    }
  }

  return res;
}

// Linear interpolation of a solution grid onto requested times.  Mirrors the
// role of rd_lookup_grid_value() in utils_reducible_diffusion.h:2146, restated
// here so this TU does not have to include the Volterra headers.
// `zero_below_grid` matters for the PDF.  A point start seeds the march at a
// short t_seed > 0 (fpe_seed), so the grid does not reach back to t = 0, and the
// default clamp would hand back the flux AT t_seed for every earlier time --
// i.e. a positive constant where the true density is ~0.  The CDF wants the
// clamp (cdf.front() is the mass genuinely absorbed before t_seed); the PDF does
// not.  t_seed is ~1 ms at M = 256, so this is a small correction, but it is
// the difference between a density that rises from zero and one that does not.
inline double grid_lookup(const std::vector<double>& tg,
                          const std::vector<double>& vg, double t,
                          bool zero_below_grid = false) {
  const int n = static_cast<int>(tg.size());
  if (n == 0) return 0.0;
  if (t <= tg.front()) {
    return (zero_below_grid && t < tg.front()) ? 0.0 : vg.front();
  }
  if (t >= tg.back())  return vg.back();
  const int j = static_cast<int>(
    std::lower_bound(tg.begin(), tg.end(), t) - tg.begin());
  const double w = (t - tg[j - 1]) / (tg[j] - tg[j - 1]);
  return vg[j - 1] + w * (vg[j] - vg[j - 1]);
}

} // namespace fpe

#endif // fpe_solver_h

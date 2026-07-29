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

// Solve a tridiagonal system in place-free form (Thomas algorithm).
//   sub[i] * x[i-1] + diag[i] * x[i] + sup[i] * x[i+1] = rhs[i]
// sub[0] and sup[n-1] are ignored.  `cprime` is scratch of length n.
inline void thomas_solve(const std::vector<double>& sub,
                         const std::vector<double>& diag,
                         const std::vector<double>& sup,
                         const std::vector<double>& rhs,
                         std::vector<double>& out,
                         std::vector<double>& cprime) {
  const int n = static_cast<int>(diag.size());
  if (n == 0) return;
  out.resize(n);
  cprime.resize(n);

  double denom = diag[0];
  cprime[0] = sup[0] / denom;
  out[0]    = rhs[0] / denom;
  for (int i = 1; i < n; ++i) {
    denom = diag[i] - sub[i] * cprime[i - 1];
    if (denom == 0.0) denom = std::numeric_limits<double>::min();
    cprime[i] = (i + 1 < n) ? sup[i] / denom : 0.0;
    out[i]    = (rhs[i] - sub[i] * out[i - 1]) / denom;
  }
  for (int i = n - 2; i >= 0; --i) {
    out[i] -= cprime[i] * out[i + 1];
  }
}

// The spatial operator L in dq/dt = L q, as a tridiagonal matrix, plus the
// diffusivity D(t) needed to turn the solution into a boundary flux.
struct FPE_Op {
  std::vector<double> sub, diag, sup;
  double D = 0.0;

  void resize(int M) { sub.assign(M, 0.0); diag.assign(M, 0.0); sup.assign(M, 0.0); }
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
template <class Model>
inline void build_op(const Model& m, double t, int M, double h, FPE_Op& op) {
  const double Lb  = m.length(t);       // a(t) - x_lo
  const double Lp  = m.length_prime(t);
  const double Bt  = m.B() / Lb;
  const double D   = 0.5 * Bt * Bt;
  const double r   = D / (h * h);
  op.D = D;
  op.resize(M);

  // Face Peclet-weighted coefficients.  Face j sits at xi = j*h, j = 1..M-1.
  //   F_j = (D/h) * ( bern(-P_j) * q[j-1] - bern(P_j) * q[j] ),  P_j = Atil_j*h/D
  // Cache them so each face is evaluated once, not twice.
  //   wm[j] = bern(-P_j) multiplies q[j-1];  wp[j] = bern(P_j) multiplies q[j].
  // (allocated on the caller's stack via static thread_local scratch)
  static thread_local std::vector<double> wm, wp;
  wm.assign(M + 1, 0.0);
  wp.assign(M + 1, 0.0);
  for (int j = 1; j <= M - 1; ++j) {
    const double xi = j * h;
    const double x  = m.x_lo() + xi * Lb;
    const double Atil = (m.drift(x, t) - xi * Lp) / Lb;
    const double P = Atil * h / D;
    wm[j] = bern(-P);
    wp[j] = bern(P);
  }

  // dq_i/dt = -( F_{i+1} - F_i ) / h
  for (int i = 0; i < M; ++i) {
    double a_sub = 0.0, a_dia = 0.0, a_sup = 0.0;

    if (i > 0) {                      // inflow face F_i is an interior face
      a_sub += r * wm[i];
      a_dia -= r * wp[i];
      // NB sign: +F_i/h contributes +(D/h^2)(wm[i] q[i-1] - wp[i] q[i])
    }
    // else: F_0 = 0, the no-flux far-field face.  This is what makes
    //       CDF = 1 - h*sum(q) exact.

    if (i < M - 1) {                  // outflow face F_{i+1} is an interior face
      a_dia -= r * wm[i + 1];
      a_sup += r * wp[i + 1];
      // -F_{i+1}/h = -(D/h^2)(wm[i+1] q[i] - wp[i+1] q[i+1])
    } else {
      // Absorbing face at xi = 1.  q(1) = 0, so advection carries nothing and
      // the flux is purely diffusive:  F_M = D * dq/ds |_{s=0}, s = 1 - xi.
      // 2nd-order one-sided through (s=0, 0), (s=h/2, q[M-1]), (s=3h/2, q[M-2]):
      //   (-8*u0 + 9*u1 - u2) / (3h)  =>  (9 q[M-1] - q[M-2]) / (3h)
      a_dia -= r * 3.0;
      a_sub += r / 3.0;
    }

    op.sub[i]  = a_sub;
    op.diag[i] = a_dia;
    op.sup[i]  = a_sup;
  }
}

// Absorbing-face flux = the first-passage density g(t).
inline double flux_out(const FPE_Op& op, const std::vector<double>& q, double h) {
  const int M = static_cast<int>(q.size());
  if (M < 2) return 0.0;
  return op.D * (9.0 * q[M - 1] - q[M - 2]) / (3.0 * h);
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
  double flux_mass_mismatch = 0.0;   // per-solve error estimate, see below
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// March q from t0 to t_max.
//
//   q0        cell-averaged sub-density at t0 on the normalised grid, already
//             carrying whatever mass was absorbed before t0 (see the seeding in
//             fpe_models.h) -- i.e. h*sum(q0) = 1 - CDF(t0).
//   nt        number of Crank-Nicolson steps of size dt = (t_max-t0)/nt.  The
//             first two are replaced by four backward-Euler half-steps
//             (Rannacher), which is what keeps CN from ringing on a point start.
//
// Two independent routes to the cdf are computed:
//   mass flavour:  CDF = 1 - h*sum(q)          (exact, given F_0 = 0)
//   flux flavour:  CDF = CDF(t0) + trapz(g)
// Their maximum absolute difference is returned as flux_mass_mismatch.  It is a
// genuine convergence diagnostic: the mass identity ALONE is a tautology of the
// discretisation and checks nothing, but the two agree only when the scheme has
// resolved the solution.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class Model>
inline FPE_Result fpe_solve(const Model& m, const std::vector<double>& q0,
                            double t0, double t_max, int M, int nt) {
  FPE_Result res;
  const int n_steps = std::max(nt, 1);
  const double h  = 1.0 / M;
  const double dt = (t_max - t0) / n_steps;

  std::vector<double> q = q0, rhs, cprime;
  std::vector<double> sub(M), dia(M), sup(M);
  FPE_Op op_old, op_new;

  const int n_rann = std::min(2, n_steps);   // full steps replaced by BE halves
  const int n_lev  = 1 + 2 * n_rann + (n_steps - n_rann);
  res.t.reserve(n_lev);
  res.pdf.reserve(n_lev);
  res.cdf.reserve(n_lev);

  double t = t0;
  build_op(m, t, M, h, op_old);

  double mass = 0.0;
  for (int i = 0; i < M; ++i) mass += q[i];
  mass *= h;

  const double cdf0 = 1.0 - mass;
  double cdf_flux = cdf0;
  double g_prev = flux_out(op_old, q, h);

  res.t.push_back(t);
  res.pdf.push_back(g_prev);
  res.cdf.push_back(cdf0);
  res.flux_mass_mismatch = 0.0;

  // step_kind: true = backward Euler, false = Crank-Nicolson
  auto do_step = [&](double step, bool backward_euler) {
    const double t_new = t + step;
    build_op(m, t_new, M, h, op_new);

    if (backward_euler) {
      rhs = q;
    } else {
      apply_shifted(op_old, 0.5 * step, q, rhs);
    }
    const double c = backward_euler ? step : 0.5 * step;

    // (I - c*L^{n+1}) q^{n+1} = rhs
    for (int i = 0; i < M; ++i) {
      sub[i] = -c * op_new.sub[i];
      dia[i] = 1.0 - c * op_new.diag[i];
      sup[i] = -c * op_new.sup[i];
    }
    thomas_solve(sub, dia, sup, rhs, q, cprime);

    t = t_new;
    op_old = op_new;

    double s = 0.0;
    for (int i = 0; i < M; ++i) s += q[i];
    const double cdf_mass = 1.0 - h * s;

    const double g = flux_out(op_old, q, h);
    cdf_flux += 0.5 * step * (g_prev + g);
    g_prev = g;

    res.t.push_back(t);
    res.pdf.push_back(g);
    res.cdf.push_back(cdf_mass);
    res.flux_mass_mismatch =
      std::max(res.flux_mass_mismatch, std::abs(cdf_mass - cdf_flux));
  };

  for (int k = 0; k < n_rann; ++k) {         // Rannacher start
    do_step(0.5 * dt, true);
    do_step(0.5 * dt, true);
  }
  for (int k = n_rann; k < n_steps; ++k) {   // Crank-Nicolson
    do_step(dt, false);
  }

  return res;
}

// Linear interpolation of a solution grid onto requested times.  Mirrors the
// role of rd_lookup_grid_value() in utils_reducible_diffusion.h:2146, restated
// here so this TU does not have to include the Volterra headers.
inline double grid_lookup(const std::vector<double>& tg,
                          const std::vector<double>& vg, double t) {
  const int n = static_cast<int>(tg.size());
  if (n == 0) return 0.0;
  if (t <= tg.front()) return vg.front();
  if (t >= tg.back())  return vg.back();
  const int j = static_cast<int>(
    std::lower_bound(tg.begin(), tg.end(), t) - tg.begin());
  const double w = (t - tg[j - 1]) / (tg[j] - tg[j - 1]);
  return vg[j - 1] + w * (vg[j] - vg[j - 1]);
}

} // namespace fpe

#endif // fpe_solver_h

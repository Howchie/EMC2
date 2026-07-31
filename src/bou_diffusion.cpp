// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R entry points for the bounded OU (Smith & Ratcliff 2004) -- the DDM with
// leak.  Density, cdf and simulator.
//
// Shares fpe_bou.h with the likelihood kernels: same key construction, same
// cache, same quadrature over across-trial variability, so the R-side density
// and the sampled likelihood cannot drift apart.  This is the arrangement
// rou_diffusion.cpp uses for the racing OU, and for the same reason.
//
// Response coding follows the DDM throughout: 1 = lower, 2 = upper, matching
// R's first and second factor levels, and the density/cdf are DEFECTIVE (they
// sum to 1 across the two responses, not each alone).
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// RcppArmadillo, not Rcpp: fpe_bou.h pulls in gh_quad.h, which includes
// RcppArmadillo.h, and that header refuses to be reached after plain Rcpp.h.
#include <RcppArmadillo.h>
#include "fpe_race.h"
#include "fpe_bou.h"

using namespace Rcpp;

namespace {

// Map the R-facing (Z, SZ) proportions onto absolute (z, sz).
//
// SZ arrives RAW -- Ttransform is an R-side step and does not run on this path,
// so the widening is done here, exactly as d_DDM_Wien_raw does it from its own
// raw SZ column (model_DDM.h:122-123):
//
//   sz_relative = 2 * SZ * min(Z, 1 - Z),   sz = a * sz_relative
//
// which is what the expression below computes, since z = Z*a and a - z = (1-Z)*a.
// The start point is then uniform on z +/- sz/2, centred, as in the DDM.
//
// NB the R-level dDDM()/pDDM() wrappers are NOT the same contract: they hand
// their SZ column straight to WienR as `sw`, because by the time they are called
// Ttransform has already widened it.  Comparing this function against dDDM()
// therefore requires passing dDDM the widened value, or the two disagree by
// ~0.6% at Z = 0.45 for reasons that have nothing to do with the solver.
inline void bou_zsz(double a, double Z, double SZ, double& z, double& sz) {
  z = Z * a;
  sz = 2.0 * SZ * std::min(z, a - z);
}

fpebou::SolveCache& scratch_cache() {
  static fpebou::SolveCache C;
  return C;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Defective density and cdf for a vector of (rt, R) under per-row parameters.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List bou_pdf_cdf_vec(NumericVector rt, IntegerVector R,
                           NumericVector v, NumericVector a, NumericVector Z,
                           NumericVector sv, NumericVector SZ,
                           NumericVector t0, NumericVector st0,
                           NumericVector s, NumericVector beta,
                           int bkind = 0,
                           NumericVector aInf = NumericVector::create(),
                           NumericVector tau = NumericVector::create(),
                           NumericVector pw = NumericVector::create(),
                           int nx = 512, double dt_target = 5e-4,
                           double grade = 1.0, double tgrade = 32.0,
                           int n_sv = 7, int n_sz = 7, int n_st0 = 7,
                           bool anchor_at_z = true, double anchor_fix = NA_REAL,
                           bool want_cdf = true) {
  const int n = rt.size();
  NumericVector pdf(n), cdf(n);

  fpebou::SolveCache& C = scratch_cache();
  C.grid.nx = nx;
  C.grid.dt_target = dt_target;
  C.grid.grade = grade;
  C.grid.tgrade = tgrade;
  C.grid.n_sv = n_sv;
  C.grid.n_sz = n_sz;
  C.grid.n_st0 = n_st0;
  C.new_particle();

  auto at = [&](const NumericVector& x, int i) {
    return x.size() == 1 ? x[0] : x[i];
  };
  // Collapse columns are absent under a fixed bound and are never read there.
  auto atc = [&](const NumericVector& x, int i) {
    return x.size() == 0 ? 0.0 : (x.size() == 1 ? x[0] : x[i]);
  };

  // Pass 1: the horizon.  Every solve is taken to the largest decision time any
  // row will ask about, so a sorted rt vector does not re-solve each quadrature
  // node as it walks (see SolveCache::t_horizon).  st0 shifts the query time
  // down by at most half its width, so the largest decision time uses the
  // smallest t0 in the range.
  double t_h = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!R_finite(rt[i])) continue;
    const double td = rt[i] - (at(t0, i) - 0.5 * std::max(0.0, at(st0, i)));
    if (td > t_h) t_h = td;
  }
  C.t_horizon = std::max(t_h * 1.02, 1e-3);

  // Pass 2.
  for (int i = 0; i < n; ++i) {
    const double rti = rt[i];
    if (!R_finite(rti)) continue;
    const double ai = at(a, i);
    double zi, szi;
    bou_zsz(ai, at(Z, i), at(SZ, i), zi, szi);

    const fpebou::BouMix m = fpebou::bou_mix_t0(
      C, rti, at(v, i), at(sv, i), ai, zi, szi, at(s, i), at(beta, i),
      anchor_at_z, R_finite(anchor_fix) ? anchor_fix : 0.5 * ai,
      at(t0, i), at(st0, i),
      bkind, atc(aInf, i), atc(tau, i), atc(pw, i), want_cdf);

    const bool up = (R[i] == 2);
    pdf[i] = up ? m.d_up : m.d_lo;
    if (want_cdf) cdf[i] = up ? m.F_up : m.F_lo;
  }
  return Rcpp::List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                            _["n_solves"] = static_cast<double>(C.n_entries));
}

// ---------------------------------------------------------------------------
// Simulator.
//
// Exact OU transition per step -- no Euler discretisation of the drift -- with
// the k -> 0 limits taken through expm1 so beta = 0 is a regular value rather
// than a special case, and a Brownian-bridge correction at BOTH barriers for
// paths that cross and return within a step.  Written in the same shape as
// fperace::rou_hit_time_bnd (fpe_race.h:1238), extended to two barriers.
//
// Without the bridge correction a simulator systematically over-estimates RTs
// and under-estimates error rates, which would make it useless as a check on
// the solver -- the one job it has here.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::DataFrame rbou_cpp(int n, NumericVector v, NumericVector a,
                         NumericVector Z, NumericVector sv, NumericVector SZ,
                         NumericVector t0, NumericVector st0,
                         NumericVector s, NumericVector beta,
                         int bkind = 0,
                         NumericVector aInf = NumericVector::create(),
                         NumericVector tau = NumericVector::create(),
                         NumericVector pw = NumericVector::create(),
                         double dt = 1e-4, double t_max = 30.0,
                         bool anchor_at_z = true, double anchor_fix = NA_REAL) {
  NumericVector rt(n);
  IntegerVector resp(n);

  auto at = [&](const NumericVector& x, int i) {
    return x.size() == 1 ? x[0] : x[i];
  };
  auto atc = [&](const NumericVector& x, int i) {
    return x.size() == 0 ? 0.0 : (x.size() == 1 ? x[0] : x[i]);
  };

  const double phi_eps = 1e-10;
  Rcpp::RNGScope scope;

  for (int i = 0; i < n; ++i) {
    const double ai = at(a, i), si = at(s, i), bi = at(beta, i);
    double zi, szi;
    bou_zsz(ai, at(Z, i), at(SZ, i), zi, szi);

    // Across-trial draws.  These are the same three sources the solver
    // integrates over, drawn here rather than integrated.
    const double vi = at(v, i) + at(sv, i) * ((at(sv, i) > 0.0) ? R::norm_rand() : 0.0);
    const double z0 = (szi > 0.0) ? (zi - 0.5 * szi + szi * unif_rand()) : zi;
    // Uniform(t0, t0 + st0), the LOWER-EDGE convention -- the same one the
    // solver's t0 quadrature uses (fpe_bou.h:bou_mix_t0) and the one the
    // package's DDM uses.  A centred draw here would put the simulator half an
    // st0 ahead of the density it is meant to be checking.
    const double t0i = at(t0, i) +
      ((at(st0, i) > 0.0) ? (at(st0, i) * unif_rand()) : 0.0);

    if (!(ai > 0.0) || !(si > 0.0) || !(z0 > 0.0) || !(z0 < ai)) {
      rt[i] = NA_REAL; resp[i] = NA_INTEGER; continue;
    }

    const double anc = anchor_at_z ? z0
                                   : (R_finite(anchor_fix) ? anchor_fix : 0.5 * ai);
    const double v_eff = vi + bi * anc;    // drift = v_eff - beta*X

    const double phi = std::exp(-bi * dt);
    const double m1 = -std::expm1(-bi * dt);
    const double m2 = -std::expm1(-2.0 * bi * dt);
    const double drift_gain = (bi > phi_eps) ? (m1 / bi) : dt;
    const double var = (bi > phi_eps) ? (si * si * m2 / (2.0 * bi))
                                      : (si * si * dt);
    const double sd = std::sqrt(std::max(var, 0.0));
    const double inv_2var_bb = 2.0 / (si * si * dt);

    // Barrier geometry, taken from the SOLVER'S own model object rather than
    // re-derived here, so the midpoint convention and the FPE_BOU_MIN_SEP floor
    // cannot drift apart between the simulator and the density it checks.
    fpe::FPE_ModelBoundedOU geo;
    geo.set_separation(ai, bkind, atc(aInf, i), atc(tau, i), atc(pw, i));

    double X = z0, t = 0.0;
    double hit = NA_REAL;
    int which = NA_INTEGER;
    double hi0 = geo.x_hi(0.0), lo0 = geo.x_lo(0.0);

    while (t < t_max) {
      const double X1 = X * phi + v_eff * drift_gain + sd * R::norm_rand();
      const double hi1 = geo.x_hi(t + dt), lo1 = geo.x_lo(t + dt);
      t += dt;

      // Direct crossings first.  If the step ends outside both barriers -- only
      // possible for a very coarse dt -- take whichever it passed first by
      // linear interpolation.  The barrier is interpolated linearly too, so the
      // crossing solves X + f*(X1-X) = b0 + f*(b1-b0); with fixed bounds the
      // b1-b0 terms vanish and this is the expression it replaces.
      const bool up_hit = (X1 >= hi1);
      const bool lo_hit = (X1 <= lo1);
      if (up_hit || lo_hit) {
        double f_up = 2.0, f_lo = 2.0;
        if (up_hit) f_up = (hi0 - X) / std::max((X1 - X) - (hi1 - hi0), 1e-300);
        if (lo_hit) f_lo = (X - lo0) / std::max((X - X1) + (lo1 - lo0), 1e-300);
        const double f = std::min(f_up, f_lo);
        hit = t - dt + std::min(std::max(f, 0.0), 1.0) * dt;
        which = (f_up <= f_lo) ? 2 : 1;
        break;
      }

      // Crossed and came back within the step.  Both barriers are tested; when
      // both fire (rare, and only when dt is coarse relative to a/s) the earlier
      // sampled time wins, which is the same tie-break the direct case uses.
      // The bridge probability against a LINEAR boundary is the same expression
      // with the barrier evaluated at each end of the step, which is exactly why
      // the collapse costs nothing here.
      const double pu = std::exp(-(hi0 - X) * (hi1 - X1) * inv_2var_bb);
      const double pl = std::exp(-(X - lo0) * (X1 - lo1) * inv_2var_bb);
      const bool cu = (unif_rand() < pu);
      const bool cl = (unif_rand() < pl);
      if (cu || cl) {
        const double tu = cu ? unif_rand() : 2.0;
        const double tl = cl ? unif_rand() : 2.0;
        const double f = std::min(tu, tl);
        hit = t - dt + f * dt;
        which = (tu <= tl) ? 2 : 1;
        break;
      }
      X = X1;
      hi0 = hi1; lo0 = lo1;
    }

    if (R_finite(hit)) { rt[i] = hit + t0i; resp[i] = which; }
    else               { rt[i] = R_PosInf; resp[i] = NA_INTEGER; }
  }

  return Rcpp::DataFrame::create(_["rt"] = rt, _["R"] = resp);
}

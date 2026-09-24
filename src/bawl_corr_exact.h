#ifndef EMC2_BAWL_CORR_EXACT_H
#define EMC2_BAWL_CORR_EXACT_H

// Exact two-racer component for correlated LBA/BAwL trials.
//
// The fixed-time survivor is affine in the drift on each of at most two
// intervals, and a cause density is affine on its winner interval.  The
// bivariate normal rectangle moments below therefore contain every term needed
// by the pair component.  This header deliberately has no Rcpp objects: the
// likelihood route can keep the ParamTable column pointers in its hot path.

#include "bawl_corr_counters.h"
#include "bawl_geometry.h"
#include "gl_quad.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

enum class BAwLCorrMomentStatus : uint8_t {
  ok = 0,
  zero_mass,
  unstable,
  invalid
};

struct BvnRectMoments {
  double p = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m12 = 0.0;
  // Absolute error bounds, split by how each part enters a weighted term:
  // err_p for p (every moment carries mu * p), err_c1/err_c2 for the
  // derivative parts m1 - mu1 p and m2 - mu2 p, and err_k for the derivative
  // part of m12 beyond p s12.  bawl_corr_add_rect_terms() combines them with
  // the term's coefficients; see there.
  double err_p = 0.0;
  double err_c1 = 0.0;
  double err_c2 = 0.0;
  double err_k = 0.0;
  // Moment coefficients of p: m1 = mu1 p + ..., m2 = mu2 p + ...,
  // m12 = (mu1 mu2 + s12) p + ...
  double mu1 = 0.0;
  double mu2 = 0.0;
  double s12 = 0.0;
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
};

struct BvnBoundaryGrid {
  // Four boundaries are enough for the pair survivor and pair cause grids;
  // unused entries are ignored.  The arrays also make the common 3x3 case
  // allocation-free.
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  double y[4] = {0.0, 0.0, 0.0, 0.0};
  double cdf[4][4] = {{0.0}};
  double cdf_err[4][4] = {{0.0}};
  double dx[4][4] = {{0.0}};
  double dy[4][4] = {{0.0}};
  double dxx[4][4] = {{0.0}};
  double dyy[4][4] = {{0.0}};
  double dxy[4][4] = {{0.0}};
  uint8_t nx = 0;
  uint8_t ny = 0;
  double rho = 0.0;
  double sd1 = 1.0;
  double sd2 = 1.0;
  double mu1 = 0.0;
  double mu2 = 0.0;
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
};

// Beyond this many standard deviations a normal tail holds Q(9) < 1.2e-19.
// The survivor-level shortcut below drops such a tail only where the result
// it is dropped from is bounded well away from zero, so the error is relative
// as well as absolute.  (Zeroing far-tail corner CDFs is NOT safe: when the
// whole survivor lives deep in a correlated tail, those corners are the
// answer -- see the "deep correlated tails" regression test.)  A lower
// truncation bound shortly after t0 is the common case: tau is small, so each
// racer's survivor boundaries sit many SDs above its drift mean and S(LT) is
// one to double precision -- yet it was costing as much as S(UT), through a
// full grid whose far-tail densities also fall off dnormP's fast branch.
constexpr double kBawlCorrFarZ = 9.0;
constexpr double kBawlCorrFarTail = 1.2e-19;  // >= Q(kBawlCorrFarZ)
// Relative error bound of the univariate conditional survivor (a Gaussian
// partial expectation whose two parts cancel by at most z^2 < 1500 before
// phi(z) underflows).
constexpr double kBawlCorrUnivRelErr = 1e-13;

// Absolute error bounds for one corner of norm_cdf_2d_hybrid, measured
// against mvtnorm's TVPACK (abseps 1e-15) on random corners with x, y in
// [-40, 10] and |rho| < .9999.  tvpack never exceeded 3.3e-16.  Drezner's
// five-point rule depends on |rho|: max 2.6e-13 below .3, 3.2e-10 below .5,
// 1.1e-7 below .7, 3.8e-7 below .9, 1.9e-8 below .99 and 6.3e-11 above; the
// bounds below keep a factor ~3 over those maxima.  Both are ABSOLUTE.
// Deep in a negatively correlated tail tvpack cancels Phi(x)Phi(y) against
// an integral of the same size and returns noise at up to ~1e-2 of that
// product, e.g. e^-348 where the true corner is e^-1889.  The exact pair
// route therefore carries an absolute error bound (BAwLCorrPairResult::noise)
// alongside every value, and a value that does not clear it is escalated
// rather than being trusted to relative precision.
inline double bawl_corr_drezner_abs_err(double r) {
  const double a = std::fabs(r);
  if (a < 0.3) return 1e-12;
  if (a < 0.5) return 1e-9;
  if (a < 0.7) return 3e-7;
  if (a < 0.9) return 1e-6;
  if (a < 0.99) return 6e-8;
  return 2e-10;
}
constexpr double kBvnTvpackAbsErr = 1e-15;
// Relative error of the pnorm/dnorm-built corner derivatives, and of the
// double sums that combine corners.
constexpr double kBvnRelErr = 1e-15;
// The likelihood trusts an exact value only when it exceeds its error bound
// by this factor, i.e. when its relative error is at most 1e-3.
constexpr double kBawlCorrResolveRatio = 1e3;

// norm_cdf_2d_hybrid, bit for bit, plus the error bound of the branch that
// produced the value.  `precise` skips the Drezner branch: the same exact
// algebra on tvpack corners is the cheap second tier for a value that the
// hybrid's Drezner corner bound leaves unresolved.
//
// `marg`, when given, holds the boundary grid's per-axis marginals for this
// corner (see BvnFastMarginals; in Drezner's upper-orthant arguments -x, -y
// they are Phi(x), Phi(y) and Phi(-y)), and `cache` the grid's Drezner table,
// so neither is re-derived per corner.
inline double bawl_corr_bvn_cdf(double x, double y, double r, double& abs_err,
                                bool precise = false,
                                const BvnFastMarginals* marg = nullptr,
                                BvnDreznerCache* cache = nullptr) {
  if (!precise &&
      !(std::fabs(r) > 0.9999 || x < EMC2_BVN_MIN_Z || y < EMC2_BVN_MIN_Z)) {
    double fast;
    if (marg != nullptr && cache != nullptr) {
      cache->prepare(r);
      fast = norm_ucdf_2d_fast_core(-x, -y, r, *cache, *marg);
    } else {
      fast = norm_cdf_2d_fast(x, y, r);
    }
    if (fast >= EMC2_BVN_DREZNER_MIN_P) {
      abs_err = bawl_corr_drezner_abs_err(r);
      return fast;
    }
  }
  const double out = norm_cdf_2d(x, y, r);
  abs_err = kBvnTvpackAbsErr;
  return out > 0.0 ? out : 0.0;
}

inline void bawl_corr_kahan_add(double value, double& sum, double& correction) {
  const double y = value - correction;
  const double t = sum + y;
  correction = (t - sum) - y;
  sum = t;
}

struct BvnCornerValues {
  double cdf = 0.0;
  double cdf_err = 0.0;  // absolute error bound on cdf
  double dx = 0.0;
  double dy = 0.0;
  double dxx = 0.0;
  double dyy = 0.0;
  double dxy = 0.0;
  bool valid = true;
};

// dnormP(x) and pnorm_std(x) together.  In dnormP's fast branch (|x| < 5)
// both evaluate the identical exp(-x^2/2) -- the same double argument, since
// -0.5*(x*x) and -(|x|*|x|)/2 are exact rescalings of one product -- so it is
// computed once.  Same arithmetic as the two separate calls; everywhere else
// it simply makes them.
inline void bawl_corr_phi_Phi(double x, double& phi, double& Phi) {
#ifdef USE_FAST_PNORM
  const double z = std::fabs(x);
  if (z < 5.0) {
    const double e = std::exp(-z * z / 2.0);
    phi = 0.398942280401432677939946059934 * e;
    const double n = (((((FAST_NORM_N6 * z + FAST_NORM_N5) * z + FAST_NORM_N4) * z +
                        FAST_NORM_N3) * z + FAST_NORM_N2) * z + FAST_NORM_N1) * z +
                     FAST_NORM_N0;
    const double d = ((((((FAST_NORM_M7 * z + FAST_NORM_M6) * z + FAST_NORM_M5) * z +
                         FAST_NORM_M4) * z + FAST_NORM_M3) * z + FAST_NORM_M2) * z +
                      FAST_NORM_M1) * z + FAST_NORM_M0;
    const double c = e * n / d;   // z < 5 < FAST_NORM_SPLIT: rational branch
    Phi = x <= 0.0 ? c : 1.0 - c;
    return;
  }
#endif
  phi = dnormP(x);
  Phi = pnorm_std(x, true, false);
}

// fx_pre/fy_pre, when finite, are dnormP(x)/dnormP(y) supplied by a caller
// that visits the same boundary at several corners (the boundary grid below),
// so each boundary density is evaluated once rather than once per corner.
//
// px_pre/py_pre, when finite, are pnorm_std(x)/pnorm_std(y) and qy_pre is
// pnorm_std(-y), from the same per-axis pass; `dcache` is the grid's Drezner
// table (nullptr: look it up per call, as a standalone corner does).
inline BvnCornerValues bawl_corr_bvn_corner(
    double x, double y, double rho,
    double r_override = R_NaN, double sr_override = R_NaN,
    double fx_pre = R_NaN, double fy_pre = R_NaN, bool precise = false,
    double px_pre = R_NaN, double py_pre = R_NaN, double qy_pre = R_NaN,
    BvnDreznerCache* dcache = nullptr) {
  BvnCornerValues out;
  if (x == R_NegInf || y == R_NegInf) return out;
  if (x == R_PosInf && y == R_PosInf) {
    out.cdf = 1.0;
    return out;
  }
  if (x == R_PosInf) {
    out.cdf = R_FINITE(py_pre) ? py_pre : pnorm_std(y, true, false);
    out.cdf_err = kBvnRelErr * out.cdf;
    out.dy = R_FINITE(fy_pre) ? fy_pre : dnormP(y);
    out.dyy = -y * out.dy;
    return out;
  }
  if (y == R_PosInf) {
    out.cdf = R_FINITE(px_pre) ? px_pre : pnorm_std(x, true, false);
    out.cdf_err = kBvnRelErr * out.cdf;
    out.dx = R_FINITE(fx_pre) ? fx_pre : dnormP(x);
    out.dxx = -x * out.dx;
    return out;
  }
  if (!R_FINITE(x) || !R_FINITE(y) || !R_FINITE(rho)) {
    out.valid = false;
    return out;
  }

  const double r = R_FINITE(r_override)
    ? r_override : std::fmax(-1.0 + 1e-12, std::fmin(1.0 - 1e-12, rho));
  const double sr = R_FINITE(sr_override)
    ? sr_override : std::sqrt(std::fmax(1.0 - r * r, 1e-24));
  const double fx = R_FINITE(fx_pre) ? fx_pre : dnormP(x);
  const double fy = R_FINITE(fy_pre) ? fy_pre : dnormP(y);
  if (std::fabs(r) <= 1e-15) {
    const double px = R_FINITE(px_pre) ? px_pre : pnorm_std(x, true, false);
    const double py = R_FINITE(py_pre) ? py_pre : pnorm_std(y, true, false);
    out.cdf = px * py;
    out.cdf_err = kBvnRelErr * out.cdf;
    out.dx = fx * py;
    out.dy = fy * px;
    out.dxx = -x * fx * py;
    out.dyy = -y * fy * px;
    out.dxy = fx * fy;
    return out;
  }

  BvnFastMarginals marg;
  marg.g1 = px_pre; marg.has_g1 = R_FINITE(px_pre);
  marg.g2 = py_pre; marg.has_g2 = R_FINITE(py_pre);
  marg.h2 = qy_pre; marg.has_h2 = R_FINITE(qy_pre);
  out.cdf = bawl_corr_bvn_cdf(x, y, r, out.cdf_err, precise, &marg,
                              dcache != nullptr ? dcache : &bvn_drezner_cache());
  const double qx = (y - r * x) / sr;
  const double qy = (x - r * y) / sr;
  // The same conditional CDFs feed both the first and second boundary
  // derivatives.  Calling pnorm_std once per argument matters here: a
  // finite two-racer cause rectangle visits four finite corners, and the
  // lower-tail hybrid already makes the BVN CDF itself relatively expensive.
  double fxq, fyq, pxq, pyq;
  bawl_corr_phi_Phi(qx, fxq, pxq);
  bawl_corr_phi_Phi(qy, fyq, pyq);
  out.dx = fx * pxq;
  out.dy = fy * pyq;
  // d/dy of d/dx Phi2 is phi(x) phi((y-rho*x)/s) / s.  Using the
  // symmetric conditional argument here is tempting, but it is paired with
  // the wrong marginal density away from the diagonal and corrupts M12 on
  // narrow asymmetric rectangles.
  out.dxy = fx * fxq / sr;
  out.dxx = fx * (-x * pxq - r * fxq / sr);
  out.dyy = fy * (-y * pyq - r * fyq / sr);
  if (!R_FINITE(out.cdf) || !R_FINITE(out.dx) || !R_FINITE(out.dy) ||
      !R_FINITE(out.dxx) || !R_FINITE(out.dyy) || !R_FINITE(out.dxy)) {
    out.valid = false;
  }
  return out;
}

inline bool bawl_corr_add_boundary(double* values, uint8_t& n, double value) {
  for (uint8_t i = 0; i < n; ++i) {
    if (values[i] == value) return true;
  }
  if (n >= 4) return false;
  values[n++] = value;
  return true;
}

inline BvnBoundaryGrid bawl_corr_make_boundary_grid(
    double mu1, double sd1, double mu2, double sd2, double rho,
    const double* x_bounds, int n_x_bounds,
    const double* y_bounds, int n_y_bounds, bool precise = false) {
  BvnBoundaryGrid g;
  g.mu1 = mu1; g.mu2 = mu2; g.sd1 = sd1; g.sd2 = sd2; g.rho = rho;
  g.status = BAwLCorrMomentStatus::ok;
  for (int i = 0; i < n_x_bounds; ++i) {
    const double value = x_bounds[i];
    if (!bawl_corr_add_boundary(g.x, g.nx, value)) {
      g.status = BAwLCorrMomentStatus::invalid;
      return g;
    }
  }
  for (int i = 0; i < n_y_bounds; ++i) {
    const double value = y_bounds[i];
    if (!bawl_corr_add_boundary(g.y, g.ny, value)) {
      g.status = BAwLCorrMomentStatus::invalid;
      return g;
    }
  }
  double r_grid = R_NaN;
  double sr_grid = R_NaN;
  if (R_FINITE(rho)) {
    r_grid = std::fmax(-1.0 + 1e-12, std::fmin(1.0 - 1e-12, rho));
    sr_grid = std::sqrt(std::fmax(1.0 - r_grid * r_grid, 1e-24));
  }
  // Standardised boundaries, their marginal densities and marginal CDFs are
  // per axis, not per corner: a 2x2 finite grid otherwise evaluates each of
  // them twice.  Phi(-y) is needed only on Drezner's high negative-rho branch.
  // A corner whose x or y is -Inf returns before using any of them, so those
  // axes are skipped.
  const bool need_qy = R_FINITE(r_grid) && r_grid <= -0.7;
  double xs[4], ys[4], fxs[4], fys[4], pxs[4], pys[4], qys[4];
  for (uint8_t i = 0; i < g.nx; ++i) {
    const bool inf = g.x[i] == R_NegInf || g.x[i] == R_PosInf;
    xs[i] = inf ? g.x[i] : (g.x[i] - mu1) / sd1;
    const bool fin = R_FINITE(xs[i]);
    fxs[i] = fin ? dnormP(xs[i]) : R_NaN;
    pxs[i] = fin ? pnorm_std(xs[i], true, false) : R_NaN;
  }
  for (uint8_t j = 0; j < g.ny; ++j) {
    const bool inf = g.y[j] == R_NegInf || g.y[j] == R_PosInf;
    ys[j] = inf ? g.y[j] : (g.y[j] - mu2) / sd2;
    const bool fin = R_FINITE(ys[j]);
    fys[j] = fin ? dnormP(ys[j]) : R_NaN;
    pys[j] = fin ? pnorm_std(ys[j], true, false) : R_NaN;
    qys[j] = (fin && need_qy) ? pnorm_std(-ys[j], true, false) : R_NaN;
  }
  BvnDreznerCache& dcache = bvn_drezner_cache();
  for (uint8_t i = 0; i < g.nx; ++i) {
    for (uint8_t j = 0; j < g.ny; ++j) {
      const BvnCornerValues corner = bawl_corr_bvn_corner(
          xs[i], ys[j], rho, r_grid, sr_grid, fxs[i], fys[j], precise,
          pxs[i], pys[j], qys[j], &dcache);
      if (!corner.valid) g.status = BAwLCorrMomentStatus::unstable;
      g.cdf[i][j] = corner.cdf;
      g.cdf_err[i][j] = corner.cdf_err;
      g.dx[i][j] = corner.dx;
      g.dy[i][j] = corner.dy;
      g.dxx[i][j] = corner.dxx;
      g.dyy[i][j] = corner.dyy;
      g.dxy[i][j] = corner.dxy;
    }
  }
  if (bawl_corr_counters_active()) {
    bawl_corr_counters().bvn_corner_evaluations +=
      static_cast<long long>(g.nx) * static_cast<long long>(g.ny);
  }
  return g;
}

inline int bawl_corr_boundary_index(const double* x, uint8_t n, double value) {
  for (uint8_t i = 0; i < n; ++i) if (x[i] == value) return static_cast<int>(i);
  return -1;
}

inline BvnRectMoments bawl_corr_rect_from_grid(const BvnBoundaryGrid& g,
                                               double lo1, double hi1,
                                               double lo2, double hi2) {
  BvnRectMoments out;
  if (!(hi1 > lo1) || !(hi2 > lo2)) {
    out.status = BAwLCorrMomentStatus::zero_mass;
    return out;
  }
  const int il = bawl_corr_boundary_index(g.x, g.nx, lo1);
  const int ih = bawl_corr_boundary_index(g.x, g.nx, hi1);
  const int jl = bawl_corr_boundary_index(g.y, g.ny, lo2);
  const int jh = bawl_corr_boundary_index(g.y, g.ny, hi2);
  if (il < 0 || ih < 0 || jl < 0 || jh < 0 || ih <= il || jh <= jl) {
    out.status = BAwLCorrMomentStatus::invalid;
    return out;
  }

  auto signed_corners = [&](auto member) {
    double sum = 0.0, correction = 0.0;
    bawl_corr_kahan_add(member(ih, jh), sum, correction);
    bawl_corr_kahan_add(-member(il, jh), sum, correction);
    bawl_corr_kahan_add(-member(ih, jl), sum, correction);
    bawl_corr_kahan_add(member(il, jl), sum, correction);
    return sum;
  };
  const double p = signed_corners([&](int i, int j){ return g.cdf[i][j]; });
  const double px = signed_corners([&](int i, int j){ return g.dx[i][j]; });
  const double py = signed_corners([&](int i, int j){ return g.dy[i][j]; });
  const double pxx = signed_corners([&](int i, int j){ return g.dxx[i][j]; });
  const double pyy = signed_corners([&](int i, int j){ return g.dyy[i][j]; });
  const double pxy = signed_corners([&](int i, int j){ return g.dxy[i][j]; });
  // Error bounds.  The corner CDFs carry their algorithm's absolute bound;
  // the derivative corners are relative-accurate, so their corner sums are
  // bounded by kBvnRelErr times the sum of magnitudes.  In standardised
  // boundary units m1 - mu1 p = -sd1 (px + rho py),
  // m2 - mu2 p = -sd2 (rho px + py), and the derivative part of m12 is
  // sd1 sd2 (rho (pxx + pyy) + (1 + rho^2) pxy).
  {
    auto abs_corners = [&](auto member) {
      return std::fabs(member(ih, jh)) + std::fabs(member(il, jh)) +
        std::fabs(member(ih, jl)) + std::fabs(member(il, jl));
    };
    out.mu1 = g.mu1;
    out.mu2 = g.mu2;
    out.s12 = g.rho * g.sd1 * g.sd2;
    const double ar = std::fabs(g.rho);
    const double sx = abs_corners([&](int i, int j){ return g.dx[i][j]; });
    const double sy = abs_corners([&](int i, int j){ return g.dy[i][j]; });
    out.err_p = abs_corners([&](int i, int j){ return g.cdf_err[i][j]; }) +
      kBvnRelErr * abs_corners([&](int i, int j){ return g.cdf[i][j]; });
    out.err_c1 = kBvnRelErr * g.sd1 * (sx + ar * sy);
    out.err_c2 = kBvnRelErr * g.sd2 * (ar * sx + sy);
    out.err_k = kBvnRelErr * g.sd1 * g.sd2 *
      (ar * (abs_corners([&](int i, int j){ return g.dxx[i][j]; }) +
             abs_corners([&](int i, int j){ return g.dyy[i][j]; })) +
       (1.0 + g.rho * g.rho) *
         abs_corners([&](int i, int j){ return g.dxy[i][j]; }));
  }
  if (!R_FINITE(p) || !R_FINITE(px) || !R_FINITE(py) || !R_FINITE(pxx) ||
      !R_FINITE(pyy) || !R_FINITE(pxy)) {
    out.status = BAwLCorrMomentStatus::unstable;
    return out;
  }
  if (p <= 0.0) {
    // A zero corner subtraction with nonzero boundary derivatives is the
    // characteristic tiny-rectangle cancellation that the numeric-pair
    // route is meant to catch.  The rectangle's mass is bounded by its
    // boundary line densities times the finite standardized span over which
    // the corner CDFs still resolve, so when every boundary derivative is
    // this small the bound sits far below both double resolution of the
    // surrounding sums and the global 1e-10 likelihood floor: that is a
    // genuine tail rectangle, not a cancellation, and must not push the
    // whole trial onto the numeric route.
    const double boundary_scale = std::fabs(px) + std::fabs(py) +
      std::fabs(pxx) + std::fabs(pyy) + std::fabs(pxy);
    out.status = (p > -1e-12 && boundary_scale <= 1e-14)
      ? BAwLCorrMomentStatus::zero_mass : BAwLCorrMomentStatus::unstable;
    return out;
  }
  if (p > 1.0 + 1e-9) {
    out.status = BAwLCorrMomentStatus::unstable;
    return out;
  }
  out.p = std::fmin(1.0, std::fmax(0.0, p));

  // g is the derivative with respect to the marginal means.  The rectangle
  // derivatives above are with respect to standardised upper boundaries.
  const double g1 = -px / g.sd1;
  const double g2 = -py / g.sd2;
  const double s11 = g.sd1 * g.sd1;
  const double s22 = g.sd2 * g.sd2;
  const double s12 = g.rho * g.sd1 * g.sd2;
  out.m1 = g.mu1 * out.p + s11 * g1 + s12 * g2;
  out.m2 = g.mu2 * out.p + s12 * g1 + s22 * g2;

  const double h11 = pxx / (g.sd1 * g.sd1);
  const double h22 = pyy / (g.sd2 * g.sd2);
  const double h12 = pxy / (g.sd1 * g.sd2);
  const double sh11 = s11 * h11 + s12 * h12;
  const double sh12 = s11 * h12 + s12 * h22;
  const double sh21 = s12 * h11 + s22 * h12;
  const double sh22 = s12 * h12 + s22 * h22;
  const double k12 = out.p * s12 + sh11 * s12 + sh12 * s22;
  const double c1 = out.m1 - g.mu1 * out.p;
  const double c2 = out.m2 - g.mu2 * out.p;
  out.m12 = g.mu1 * g.mu2 * out.p + g.mu1 * c2 + g.mu2 * c1 + k12;
  if (!R_FINITE(out.m1) || !R_FINITE(out.m2) || !R_FINITE(out.m12)) {
    out.status = BAwLCorrMomentStatus::unstable;
    return out;
  }
  (void)sh21; (void)sh22; // retained above to document Sigma H Sigma.
  out.status = (g.status == BAwLCorrMomentStatus::unstable)
    ? BAwLCorrMomentStatus::unstable : BAwLCorrMomentStatus::ok;
  return out;
}

inline BvnRectMoments bawl_corr_bvn_rect_moments(
    double mu1, double sd1, double mu2, double sd2, double rho,
    double lo1, double hi1, double lo2, double hi2) {
  BvnRectMoments out;
  if (!(sd1 > 0.0) || !(sd2 > 0.0) || !R_FINITE(mu1) || !R_FINITE(mu2) ||
      !R_FINITE(sd1) || !R_FINITE(sd2) || !(hi1 > lo1) || !(hi2 > lo2) ||
      !R_FINITE(rho) || std::fabs(rho) > 1.0) {
    out.status = BAwLCorrMomentStatus::invalid;
    return out;
  }
  if (std::fabs(rho) > 1.0 - 1e-10) {
    out.status = BAwLCorrMomentStatus::unstable;
    return out;
  }
  const double x_bounds[2] = {lo1, hi1};
  const double y_bounds[2] = {lo2, hi2};
  BvnBoundaryGrid grid = bawl_corr_make_boundary_grid(
      mu1, sd1, mu2, sd2, rho, x_bounds, 2, y_bounds, 2);
  return bawl_corr_rect_from_grid(grid, lo1, hi1, lo2, hi2);
}

struct BAwLCorrRegion {
  double lo = 0.0;
  double hi = 0.0;
  double c0 = 0.0;
  double c1 = 0.0;
};

struct BAwLCorrRegionList {
  std::array<BAwLCorrRegion, 2> value{};
  uint8_t n = 0;

  bool empty() const { return n == 0; }
  void push_back(const BAwLCorrRegion& region) {
    if (n < value.size()) value[static_cast<size_t>(n++)] = region;
  }
  const BAwLCorrRegion* begin() const { return value.data(); }
  const BAwLCorrRegion* end() const { return value.data() + n; }
};

inline bool bawl_corr_domain_lower(bool positive, double& d) {
  d = positive ? 0.0 : R_NegInf;
  return true;
}

inline BAwLCorrRegionList bawl_corr_survivor_regions(
    const BAwLTimeGeometry& g, bool positive) {
  BAwLCorrRegionList out;
  double d = 0.0;
  bawl_corr_domain_lower(positive, d);
  if (g.status == BAwLTimeStatus::invalid) return out;
  if (g.status == BAwLTimeStatus::not_started) {
    out.push_back({d, R_PosInf, 1.0, 0.0});
  } else if (g.status == BAwLTimeStatus::infinite) {
    if (g.L > d) out.push_back({d, g.L, 1.0, 0.0});
  } else if (g.status == BAwLTimeStatus::point_start) {
    if (g.U > d) out.push_back({d, g.U, 1.0, 0.0});
  } else {
    if (g.L > d) out.push_back({d, g.L, 1.0, 0.0});
    const double lo = std::fmax(d, g.L);
    if (g.U > lo) out.push_back({lo, g.U, g.alpha, g.beta});
  }
  return out;
}

inline bool bawl_corr_cause_interval(const BAwLTimeGeometry& g, bool positive,
                                     BAwLCorrRegion& out) {
  if (g.status != BAwLTimeStatus::valid) return false;
  const double d = positive ? 0.0 : R_NegInf;
  const double lo = std::fmax(d, g.L);
  if (!(g.U > lo)) return false;
  out = {lo, g.U, g.gamma0, g.gamma1};
  return true;
}

inline BvnRectMoments bawl_corr_zero_rect() {
  BvnRectMoments out;
  out.status = BAwLCorrMomentStatus::zero_mass;
  return out;
}

inline BvnRectMoments bawl_corr_rect_grid_for_regions(
    const BvnBoundaryGrid& grid, const BAwLCorrRegion& r1,
    const BAwLCorrRegion& r2) {
  return bawl_corr_rect_from_grid(grid, r1.lo, r1.hi, r2.lo, r2.hi);
}

struct BAwLCorrPairResult {
  double value = 0.0;
  // Absolute error bound on value.  The numeric route is trusted (zero);
  // the exact route accumulates its corner error bounds here.
  double noise = 0.0;
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
};

// Adds the weighted rectangle term (a0 + a1 V1)(b0 + b1 V2) and its error
// bound.  Expanding the moments, the term is
//   Kp p + a1 (b0 + b1 mu2) c1 + b1 (a0 + a1 mu1) c2 + a1 b1 k
// with Kp = a0 b0 + a1 b0 mu1 + a0 b1 mu2 + a1 b1 (mu1 mu2 + s12), so the
// p error enters once through Kp.  Kp is the weight at the drift means and
// can be far smaller than the separate coefficients (on a thin region the
// weights cancel), which is why the parts are bounded separately rather
// than through |a0 b0| + |a1 b0 mu1| + ...  A zero-mass rectangle adds no
// value but keeps its bound: its true mass is only known to within it.
inline bool bawl_corr_add_rect_terms(const BvnRectMoments& m,
                                     double a0, double a1,
                                     double b0, double b1,
                                     double& sum, double& correction,
                                     double& noise) {
  if (m.status == BAwLCorrMomentStatus::zero_mass ||
      m.status == BAwLCorrMomentStatus::ok) {
    const double kp = a0 * b0 + a1 * b0 * m.mu1 + a0 * b1 * m.mu2 +
      a1 * b1 * (m.mu1 * m.mu2 + m.s12);
    const double bound = std::fabs(kp) * m.err_p +
      std::fabs(a1 * (b0 + b1 * m.mu2)) * m.err_c1 +
      std::fabs(b1 * (a0 + a1 * m.mu1)) * m.err_c2 +
      std::fabs(a1 * b1) * m.err_k;
    if (R_FINITE(bound)) noise += bound;
    else noise = R_PosInf;
  }
  if (m.status == BAwLCorrMomentStatus::zero_mass) return true;
  if (m.status != BAwLCorrMomentStatus::ok) return false;
  const double t0 = a0 * b0 * m.p, t1 = a1 * b0 * m.m1,
               t2 = a0 * b1 * m.m2, t3 = a1 * b1 * m.m12;
  const double term = t0 + t1 + t2 + t3;
  if (!R_FINITE(term)) return false;
  noise += kBvnRelErr * (std::fabs(t0) + std::fabs(t1) + std::fabs(t2) +
                         std::fabs(t3));
  bawl_corr_kahan_add(term, sum, correction);
  return true;
}

// Before either racer has started (t <= t0 for both, e.g. a lower truncation
// bound below t0) the pair survivor is one for every drift pair: under
// posdrift the numerator and the orthant normalizer are the same probability.
// Returning it directly skips a grid whose corners are all infinite, and
// under posdrift also avoids dividing two independently rounded evaluations
// of that probability.
inline bool bawl_corr_pair_not_started(const BAwLTimeGeometry& g1,
                                       const BAwLTimeGeometry& g2,
                                       double normalizer,
                                       BAwLCorrPairResult& out) {
  if (g1.status != BAwLTimeStatus::not_started ||
      g2.status != BAwLTimeStatus::not_started) return false;
  if (!(normalizer > 0.0) || !R_FINITE(normalizer)) {
    out.status = BAwLCorrMomentStatus::invalid;
    return true;
  }
  out.value = 1.0;
  out.status = BAwLCorrMomentStatus::ok;
  return true;
}

inline double bawl_corr_conditional_survivor(
    const BAwLTimeGeometry& loser, bool positive, double cond_mu,
    double cond_sd, double loser_given_winner);

// Whether a racer's fixed-time survivor is one for (all but a sliver of) its
// drift distribution: exactly, before t0; or up to Q(kBawlCorrFarZ) when the
// geometry's lower survivor edge (S(v) = 1 below it) sits more than
// kBawlCorrFarZ SDs above mu.
enum class BAwLCorrSurely : uint8_t { no, far, exact };

inline BAwLCorrSurely bawl_corr_racer_surely_surviving(
    const BAwLTimeGeometry& g, double mu, double sd) {
  double edge;
  switch (g.status) {
  case BAwLTimeStatus::not_started: return BAwLCorrSurely::exact;
  case BAwLTimeStatus::valid:       edge = g.L; break;
  case BAwLTimeStatus::point_start: edge = g.U; break;
  case BAwLTimeStatus::infinite:    edge = g.L; break;
  default: return BAwLCorrSurely::no;
  }
  return (sd > 0.0 && R_FINITE(mu) && R_FINITE(edge) &&
          (edge - mu) / sd > kBawlCorrFarZ)
    ? BAwLCorrSurely::far : BAwLCorrSurely::no;
}

// Survivor-level far-tail reduction.  When racer k surely survives,
//   E[S1 S2] = E[S_other] - E[S_other (1 - S_k)],
// and the dropped term is at most P(V_k above its edge) = Q(kBawlCorrFarZ).
// That bound is absolute, so each case is accepted only where the result is
// bounded away from zero and the error is therefore also relative:
//   * both surviving: the pair survivor is >= 1 - 2Q, so it is one (under
//     posdrift the ratio is >= 1 - 2Q/q, so the orthant q must not be tiny);
//   * one surviving, without joint positivity: the other racer's marginal
//     survivor, kept only when it is >= 1e-3 unless the surviving racer has
//     not started (then S_k == 1 and the reduction is exact).  A tiny marginal
//     survivor is exactly where strong negative correlation can put the
//     other racer's surviving drifts in this racer's dropped tail.
// Under posdrift the remaining racer is still conditioned on the orthant, so
// the one-surviving case keeps the grid.
inline bool bawl_corr_pair_survival_far(
    const BAwLTimeGeometry& g1, const BAwLTimeGeometry& g2,
    double mu1, double sd1, double mu2, double sd2, bool positive,
    double normalizer, BAwLCorrPairResult& out) {
  const BAwLCorrSurely s1 = bawl_corr_racer_surely_surviving(g1, mu1, sd1);
  const BAwLCorrSurely s2 = bawl_corr_racer_surely_surviving(g2, mu2, sd2);
  if (s1 == BAwLCorrSurely::no && s2 == BAwLCorrSurely::no) return false;
  if (!(normalizer > 0.0) || !R_FINITE(normalizer)) return false;
  if (s1 != BAwLCorrSurely::no && s2 != BAwLCorrSurely::no) {
    if (positive && normalizer < 1e-3) return false;
    out.value = 1.0;
    out.noise = 2.0 * kBawlCorrFarTail / normalizer;
    out.status = BAwLCorrMomentStatus::ok;
    return true;
  }
  if (positive) return false;
  const bool first = s1 != BAwLCorrSurely::no;
  const double value = first
    ? bawl_corr_conditional_survivor(g2, false, mu2, sd2, 0.0)
    : bawl_corr_conditional_survivor(g1, false, mu1, sd1, 0.0);
  if (!R_FINITE(value)) return false;
  const BAwLCorrSurely kept = first ? s1 : s2;
  if (kept == BAwLCorrSurely::far && !(value >= 1e-3)) return false;
  out.value = std::fmin(1.0, value);
  out.noise = (kept == BAwLCorrSurely::far ? kBawlCorrFarTail : 0.0) +
    kBawlCorrUnivRelErr * out.value;
  out.status = (out.value > 0.0) ? BAwLCorrMomentStatus::ok
                                 : BAwLCorrMomentStatus::zero_mass;
  return true;
}

inline BAwLCorrPairResult bawl_corr_pair_survival_exact(
    const BAwLTimeGeometry& g1, const BAwLTimeGeometry& g2,
    double mu1, double sd1, double mu2, double sd2, double rho,
    bool positive, double normalizer, bool precise = false) {
  BAwLCorrPairResult out;
  if (bawl_corr_pair_not_started(g1, g2, normalizer, out)) return out;
  if (bawl_corr_pair_survival_far(g1, g2, mu1, sd1, mu2, sd2, positive,
                                  normalizer, out)) return out;
  const auto r1 = bawl_corr_survivor_regions(g1, positive);
  const auto r2 = bawl_corr_survivor_regions(g2, positive);
  if (r1.empty() || r2.empty()) {
    out.status = BAwLCorrMomentStatus::zero_mass;
    return out;
  }
  double x[4], y[4];
  int nx = 0, ny = 0;
  for (const auto& r : r1) { x[nx++] = r.lo; x[nx++] = r.hi; }
  for (const auto& r : r2) { y[ny++] = r.lo; y[ny++] = r.hi; }
  BvnBoundaryGrid grid = bawl_corr_make_boundary_grid(
      mu1, sd1, mu2, sd2, rho, x, nx, y, ny, precise);
  if (grid.status == BAwLCorrMomentStatus::invalid) {
    out.status = grid.status; return out;
  }
  double sum = 0.0, correction = 0.0, noise = 0.0;
  for (const auto& a : r1) for (const auto& b : r2) {
    const BvnRectMoments m = bawl_corr_rect_grid_for_regions(grid, a, b);
    if (!bawl_corr_add_rect_terms(m, a.c0, a.c1, b.c0, b.c1,
                                  sum, correction, noise)) {
      out.status = (m.status == BAwLCorrMomentStatus::unstable)
        ? BAwLCorrMomentStatus::unstable : BAwLCorrMomentStatus::invalid;
      return out;
    }
  }
  if (!(normalizer > 0.0) || !R_FINITE(normalizer) || !R_FINITE(sum)) {
    out.status = BAwLCorrMomentStatus::invalid; return out;
  }
  if (sum < -1e-9 || sum > normalizer + 1e-8) {
    out.status = BAwLCorrMomentStatus::unstable; return out;
  }
  out.value = std::fmax(0.0, sum) / normalizer;
  out.noise = noise / normalizer;
  out.status = (out.value > 0.0) ? BAwLCorrMomentStatus::ok
                                 : BAwLCorrMomentStatus::zero_mass;
  return out;
}

inline double bawl_corr_univ_interval_prob(double mu, double sd,
                                           double lo, double hi) {
  if (!(hi > lo)) return 0.0;
  const double zl = (lo == R_NegInf) ? R_NegInf : (lo - mu) / sd;
  const double zh = (hi == R_PosInf) ? R_PosInf : (hi - mu) / sd;
  if (zh <= zl) return 0.0;
  if (zl == R_NegInf) return pnorm_std(zh, true, false);
  if (zh == R_PosInf) return pnorm_std(zl, false, false);
  const double lp = log_normal_interval(zl, zh);
  return R_FINITE(lp) ? std::exp(lp) : 0.0;
}

inline void bawl_corr_univ_interval_moments(double mu, double sd,
                                            double lo, double hi,
                                            double& p, double& m1) {
  p = bawl_corr_univ_interval_prob(mu, sd, lo, hi);
  if (!(p > 0.0)) { m1 = 0.0; return; }
  const double zl = (lo == R_NegInf) ? R_NegInf : (lo - mu) / sd;
  const double zh = (hi == R_PosInf) ? R_PosInf : (hi - mu) / sd;
  const double fl = (zl == R_NegInf || zl == R_PosInf) ? 0.0 : dnormP(zl);
  const double fh = (zh == R_NegInf || zh == R_PosInf) ? 0.0 : dnormP(zh);
  m1 = mu * p + sd * (fl - fh);
  if (!R_FINITE(m1)) m1 = mu * p;
}

inline void bawl_corr_univ_interval_moments3(double mu, double sd,
                                             double lo, double hi,
                                             double& p, double& m1,
                                             double& m2, double& m3) {
  bawl_corr_univ_interval_moments(mu, sd, lo, hi, p, m1);
  if (!(p > 0.0)) { m2 = 0.0; m3 = 0.0; return; }
  const double zl = (lo == R_NegInf) ? R_NegInf : (lo - mu) / sd;
  const double zh = (hi == R_PosInf) ? R_PosInf : (hi - mu) / sd;
  const double fl = (zl == R_NegInf || zl == R_PosInf) ? 0.0 : dnormP(zl);
  const double fh = (zh == R_NegInf || zh == R_PosInf) ? 0.0 : dnormP(zh);
  const double ez2 = p +
    ((zl == R_NegInf || zl == R_PosInf) ? 0.0 : zl * fl) -
    ((zh == R_NegInf || zh == R_PosInf) ? 0.0 : zh * fh);
  const double ez3 =
    ((zl == R_NegInf || zl == R_PosInf) ? 0.0 : (zl * zl + 2.0) * fl) -
    ((zh == R_NegInf || zh == R_PosInf) ? 0.0 : (zh * zh + 2.0) * fh);
  m2 = mu * mu * p + 2.0 * mu * sd * (fl - fh) + sd * sd * ez2;
  m3 = mu * mu * mu * p + 3.0 * mu * mu * sd * (fl - fh) +
    3.0 * mu * sd * sd * ez2 + sd * sd * sd * ez3;
  if (!R_FINITE(m2)) m2 = mu * m1;
  if (!R_FINITE(m3)) m3 = mu * m2;
}

inline double bawl_corr_conditional_survivor(
    const BAwLTimeGeometry& loser, bool positive, double cond_mu,
    double cond_sd, double loser_given_winner) {
  const auto regions = bawl_corr_survivor_regions(loser, positive);
  double sum = 0.0, correction = 0.0;
  for (const auto& r : regions) {
    double p = 0.0, m = 0.0;
    bawl_corr_univ_interval_moments(cond_mu, cond_sd, r.lo, r.hi, p, m);
    const double term = r.c0 * p + r.c1 * m;
    if (term < -1e-10 || !R_FINITE(term)) return R_NaN;
    bawl_corr_kahan_add(term, sum, correction);
  }
  (void)loser_given_winner;
  return std::fmax(0.0, sum);
}

inline double bawl_corr_marginal_density(double v, double mu, double sd) {
  const double z = (v - mu) / sd;
  return dnormP(z) / sd;
}

inline BAwLCorrPairResult bawl_corr_pair_cause_exact(
    const BAwLTimeGeometry& gw, const BAwLTimeGeometry& gl,
    double muw, double sdw, double mul, double sdl, double rho,
    bool positive, double normalizer, bool precise = false) {
  BAwLCorrPairResult out;
  double d = 0.0;
  bawl_corr_domain_lower(positive, d);
  if (gw.status == BAwLTimeStatus::point_start) {
    const double vstar = gw.U;
    if (!(vstar > d) || !(gw.dv_dt > 0.0)) {
      out.status = BAwLCorrMomentStatus::zero_mass; return out;
    }
    const double cond_mu = mul + rho * (sdl / sdw) * (vstar - muw);
    const double cond_sd = sdl * std::sqrt(std::fmax(1.0 - rho * rho, 1e-24));
    const double e = bawl_corr_conditional_survivor(gl, positive, cond_mu,
                                                     cond_sd, vstar);
    const double f = bawl_corr_marginal_density(vstar, muw, sdw);
    const double value = f * gw.dv_dt * e / normalizer;
    if (!R_FINITE(value) || value < -1e-12) {
      out.status = BAwLCorrMomentStatus::unstable; return out;
    }
    out.value = std::fmax(0.0, value);
    out.noise = kBawlCorrUnivRelErr * out.value;
    out.status = out.value > 0.0 ? BAwLCorrMomentStatus::ok
                                 : BAwLCorrMomentStatus::zero_mass;
    return out;
  }
  BAwLCorrRegion wr;
  if (!bawl_corr_cause_interval(gw, positive, wr)) {
    out.status = BAwLCorrMomentStatus::zero_mass; return out;
  }
  const auto lr = bawl_corr_survivor_regions(gl, positive);
  if (lr.empty()) { out.status = BAwLCorrMomentStatus::zero_mass; return out; }
  const double x[2] = {wr.lo, wr.hi};
  double y[4];
  int ny = 0;
  for (const auto& r : lr) { y[ny++] = r.lo; y[ny++] = r.hi; }
  BvnBoundaryGrid grid = bawl_corr_make_boundary_grid(
      muw, sdw, mul, sdl, rho, x, 2, y, ny, precise);
  if (grid.status == BAwLCorrMomentStatus::invalid) {
    out.status = grid.status; return out;
  }
  double sum = 0.0, correction = 0.0, noise = 0.0;
  for (const auto& r : lr) {
    const BvnRectMoments m = bawl_corr_rect_from_grid(
        grid, wr.lo, wr.hi, r.lo, r.hi);
    if (!bawl_corr_add_rect_terms(m, wr.c0, wr.c1, r.c0, r.c1,
                                  sum, correction, noise)) {
      out.status = BAwLCorrMomentStatus::unstable; return out;
    }
  }
  if (!(normalizer > 0.0) || !R_FINITE(sum)) {
    out.status = BAwLCorrMomentStatus::invalid; return out;
  }
  if (sum < -1e-9) { out.status = BAwLCorrMomentStatus::unstable; return out; }
  out.value = std::fmax(0.0, sum) / normalizer;
  out.noise = noise / normalizer;
  out.status = out.value > 0.0 ? BAwLCorrMomentStatus::ok
                               : BAwLCorrMomentStatus::zero_mass;
  return out;
}

// The deterministic one-dimensional route is deliberately independent of the
// rectangle derivative algebra.  It is the route for pairs the exact grid
// cannot resolve (see BAwLCorrPairResult::noise) and the scalar oracle in
// tests, so it must stay accurate to relative precision wherever the mass
// is.  A window fixed in marginal SD units cannot promise that: strong
// correlation puts a whole survivor 17 SD out, and the cause density at an
// RT deep in the winner's tail lives 26+ SD out.  The mass is located
// instead.  Every pair integrand is log-concave in the first racer's
// standardised drift x -- a normal density, times a weight that is affine
// and positive on its region, times the other racer's conditional survivor,
// a Gaussian smoothing of a log-concave clamp (Prekopa) -- so it is unimodal
// with monotone flanks and {g > -Inf} is an interval:
//   1. scan a grid (plus caller breakpoints) and golden-section the mode;
//   2. bisect each flank for where g has fallen kBawlCorrNumericDrop nats;
//   3. integrate exp(g - g_mode) between them with adaptive Gauss-Kronrod
//      (G7/K15) panels that meet at the mode and at the breakpoints.
// Beyond |x| = kBawlCorrNumericZMax the normal factor alone is below e^-800.
constexpr double kBawlCorrNumericZMax = 40.0;
constexpr double kBawlCorrNumericDrop = 36.0;
constexpr int kBawlCorrNumericScan = 17;
constexpr int kBawlCorrNumericMaxBreaks = 4;
constexpr int kBawlCorrNumericMaxPanels = 64;
constexpr double kBawlCorrNumericRelTol = 1e-9;

struct BAwLCorrGKPanel {
  double a = 0.0, b = 0.0, value = 0.0, err = 0.0;
};

template <typename F>
inline BAwLCorrGKPanel bawl_corr_gk15(F f, double a, double b) {
  static constexpr double xgk[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.0};
  static constexpr double wgk[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
  static constexpr double wg[4] = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
  const double c = 0.5 * (a + b);
  const double h = 0.5 * (b - a);
  const double fc = f(c);
  double rk = wgk[7] * fc;
  double rg = wg[3] * fc;
  for (int j = 0; j < 7; ++j) {
    const double dx = h * xgk[j];
    const double pair = f(c - dx) + f(c + dx);
    rk += wgk[j] * pair;
    if (j % 2 == 1) rg += wg[j / 2] * pair;
  }
  BAwLCorrGKPanel out;
  out.a = a; out.b = b;
  out.value = rk * h;
  out.err = std::fabs(rk - rg) * h;
  return out;
}

// Returns the integral over standardised x of exp(log_integrand(v, x)),
// v = mu + sd x, for v in (lo, hi).  `breaks` are standardised x values
// where the integrand may change sharply (they need not lie inside).
template <typename LogIntegrand>
inline double bawl_corr_numeric_integral(double mu, double sd, double lo,
                                         double hi, LogIntegrand log_integrand,
                                         const double* breaks = nullptr,
                                         int n_breaks = 0) {
  double zl = (lo == R_NegInf) ? -kBawlCorrNumericZMax : (lo - mu) / sd;
  double zh = (hi == R_PosInf) ? kBawlCorrNumericZMax : (hi - mu) / sd;
  zl = std::fmax(zl, -kBawlCorrNumericZMax);
  zh = std::fmin(zh, kBawlCorrNumericZMax);
  if (!(zh > zl)) return 0.0;
  long long n_eval = 0;
  auto g = [&](double x) {
    ++n_eval;
    const double value = log_integrand(mu + sd * x, x);
    return ISNAN(value) ? R_NegInf : value;
  };
  struct EvalCount {
    long long& n;
    ~EvalCount() {
      if (bawl_corr_counters_active()) {
        ++bawl_corr_counters().numeric_pair_integrals;
        bawl_corr_counters().numeric_pair_integrand_evaluations += n;
      }
    }
  } eval_count{n_eval};

  // 1. Scan, then golden-section the mode inside the best point's
  // neighbours.  The best point stays in the bracket, so an interior -Inf
  // probe can never lose the mode.
  constexpr int kMaxScan = kBawlCorrNumericScan + kBawlCorrNumericMaxBreaks;
  std::array<double, kMaxScan> xs{}, gs{};
  int n = 0;
  for (int i = 0; i < kBawlCorrNumericScan; ++i)
    xs[static_cast<size_t>(n++)] =
      zl + (zh - zl) * i / (kBawlCorrNumericScan - 1.0);
  xs[static_cast<size_t>(n - 1)] = zh;
  for (int i = 0; i < n_breaks && i < kBawlCorrNumericMaxBreaks; ++i)
    if (breaks[i] > zl && breaks[i] < zh) xs[static_cast<size_t>(n++)] = breaks[i];
  std::sort(xs.begin(), xs.begin() + n);
  int best = -1;
  for (int i = 0; i < n; ++i) {
    gs[static_cast<size_t>(i)] = g(xs[static_cast<size_t>(i)]);
    if (gs[static_cast<size_t>(i)] > R_NegInf &&
        (best < 0 || gs[static_cast<size_t>(i)] > gs[static_cast<size_t>(best)]))
      best = i;
  }
  if (best < 0) return 0.0;
  const int ia = std::max(best - 1, 0), ib = std::min(best + 1, n - 1);
  double a = xs[static_cast<size_t>(ia)], ga = gs[static_cast<size_t>(ia)];
  double b = xs[static_cast<size_t>(ib)], gb = gs[static_cast<size_t>(ib)];
  double xm = xs[static_cast<size_t>(best)];
  double gm = gs[static_cast<size_t>(best)];
  // The mode only places a panel cut and sets the flank level, so it is
  // refined until the bracket is flat to 0.01 nats (a concave g cannot then
  // rise much above gm inside it) or negligibly narrow.
  const double x_tol = 1e-6 * (zh - zl);
  for (int it = 0; it < 60 && b - a > x_tol; ++it) {
    if (gm - ga < 0.01 && gm - gb < 0.01) break;
    const double left = xm - a, right = b - xm;
    const double u = (left > right) ? xm - 0.3819660112501051 * left
                                    : xm + 0.3819660112501051 * right;
    const double gu = g(u);
    if (gu > gm) {
      if (u < xm) { b = xm; gb = gm; } else { a = xm; ga = gm; }
      xm = u; gm = gu;
    } else {
      if (u < xm) { a = u; ga = gu; } else { b = u; gb = gu; }
    }
  }
  if (!R_FINITE(gm)) return gm == R_PosInf ? R_PosInf : 0.0;

  // 2. Flanks.  g is nondecreasing left of the mode and nonincreasing right
  // of it, so the scan brackets each level crossing; bisection keeps the
  // outer (below-level) end, erring toward a wider window.
  const double level = gm - kBawlCorrNumericDrop;
  auto flank = [&](bool left_side) {
    double inner = xm, outer = left_side ? zl : zh;
    // xs[0] and xs[n - 1] are exactly zl and zh.
    if (!(gs[static_cast<size_t>(left_side ? 0 : n - 1)] < level)) return outer;
    for (int i = 0; i < n; ++i) {
      const double x = xs[static_cast<size_t>(i)];
      if (left_side ? (x >= xm) : (x <= xm)) continue;
      const bool below = gs[static_cast<size_t>(i)] < level;
      if (left_side) {
        if (below) outer = std::fmax(outer, x); else inner = std::fmin(inner, x);
      } else {
        if (below) outer = std::fmin(outer, x); else inner = std::fmax(inner, x);
      }
    }
    if (left_side ? !(outer < inner) : !(outer > inner)) return outer;
    // The outer end is kept, so precision here only trims the window.
    for (int it = 0; it < 8; ++it) {
      const double mid = 0.5 * (inner + outer);
      if (g(mid) < level) outer = mid; else inner = mid;
    }
    return outer;
  };
  const double wa = flank(true);
  const double wb = flank(false);

  // 3. Adaptive G7/K15 on exp(g - gm), panels split at the mode and at the
  // breakpoints inside the window.
  auto f = [&](double x) { return std::exp(g(x) - gm); };
  std::array<double, 2 + kBawlCorrNumericMaxBreaks> cuts{};
  int n_cuts = 0;
  cuts[static_cast<size_t>(n_cuts++)] = wa;
  if (xm > wa && xm < wb) cuts[static_cast<size_t>(n_cuts++)] = xm;
  for (int i = 0; i < n_breaks && i < kBawlCorrNumericMaxBreaks; ++i)
    if (breaks[i] > wa && breaks[i] < wb) cuts[static_cast<size_t>(n_cuts++)] = breaks[i];
  std::sort(cuts.begin(), cuts.begin() + n_cuts);
  std::array<BAwLCorrGKPanel, kBawlCorrNumericMaxPanels> panels{};
  int n_panels = 0;
  for (int i = 0; i < n_cuts; ++i) {
    const double pa = cuts[static_cast<size_t>(i)];
    const double pb = (i + 1 < n_cuts) ? cuts[static_cast<size_t>(i + 1)] : wb;
    if (pb > pa) panels[static_cast<size_t>(n_panels++)] = bawl_corr_gk15(f, pa, pb);
  }
  if (n_panels == 0) return 0.0;
  for (;;) {
    double total = 0.0, err = 0.0;
    int worst = 0;
    for (int i = 0; i < n_panels; ++i) {
      total += panels[static_cast<size_t>(i)].value;
      err += panels[static_cast<size_t>(i)].err;
      if (panels[static_cast<size_t>(i)].err > panels[static_cast<size_t>(worst)].err)
        worst = i;
    }
    if (!(err > kBawlCorrNumericRelTol * total) ||
        n_panels >= kBawlCorrNumericMaxPanels) {
      return (total > 0.0 && R_FINITE(total)) ? std::exp(gm) * total : 0.0;
    }
    const BAwLCorrGKPanel w = panels[static_cast<size_t>(worst)];
    const double mid = 0.5 * (w.a + w.b);
    if (!(mid > w.a && mid < w.b)) {
      panels[static_cast<size_t>(worst)].err = 0.0;
      continue;
    }
    panels[static_cast<size_t>(worst)] = bawl_corr_gk15(f, w.a, mid);
    panels[static_cast<size_t>(n_panels++)] = bawl_corr_gk15(f, mid, w.b);
  }
}

// Where the other racer's conditional mean crosses one of its survivor
// region boundaries, in the first racer's standardised units.  The
// conditional survivor bends there over sqrt(1 - rho^2)/|rho| SD, which a
// panel only resolves as a breakpoint once correlation is extreme.
inline int bawl_corr_numeric_breaks(const BAwLTimeGeometry& other,
                                    bool positive, double mu_other,
                                    double sd_other, double rho,
                                    double* breaks) {
  if (!(std::fabs(rho) > 0.99)) return 0;
  int n = 0;
  const auto regions = bawl_corr_survivor_regions(other, positive);
  for (const auto& r : regions) {
    for (double edge : {r.lo, r.hi}) {
      if (!R_FINITE(edge) || n >= kBawlCorrNumericMaxBreaks) continue;
      bool seen = false;
      const double x = (edge - mu_other) / (rho * sd_other);
      for (int i = 0; i < n; ++i) seen = seen || breaks[i] == x;
      if (!seen && R_FINITE(x)) breaks[n++] = x;
    }
  }
  return n;
}

inline double bawl_corr_log_conditional_survivor(
    const BAwLTimeGeometry& loser, bool positive, double cond_mu,
    double cond_sd) {
  const double value = bawl_corr_conditional_survivor(
      loser, positive, cond_mu, cond_sd, 0.0);
  return value > 0.0 && R_FINITE(value) ? std::log(value) : R_NegInf;
}

inline BAwLCorrPairResult bawl_corr_pair_survival_numeric(
    const BAwLTimeGeometry& gw, const BAwLTimeGeometry& gl,
    double muw, double sdw, double mul, double sdl, double rho,
    bool positive, double normalizer) {
  BAwLCorrPairResult out;
  if (bawl_corr_pair_not_started(gw, gl, normalizer, out)) return out;
  if (bawl_corr_pair_survival_far(gw, gl, muw, sdw, mul, sdl, positive,
                                  normalizer, out)) return out;
  const auto wr = bawl_corr_survivor_regions(gw, positive);
  if (wr.empty()) { out.status = BAwLCorrMomentStatus::zero_mass; return out; }
  double breaks[kBawlCorrNumericMaxBreaks];
  const int n_breaks = bawl_corr_numeric_breaks(gl, positive, mul, sdl, rho,
                                                breaks);
  double total = 0.0;
  for (const auto& r : wr) {
    const double piece = bawl_corr_numeric_integral(
        muw, sdw, r.lo, r.hi,
        [&](double v, double x) {
          const double cm = mul + rho * (sdl / sdw) * (v - muw);
          const double cs = sdl * std::sqrt(std::fmax(1.0 - rho * rho, 1e-24));
          const double e = bawl_corr_log_conditional_survivor(gl, positive, cm, cs);
          const double sw = r.c0 + r.c1 * v;
          if (!(sw > 0.0)) return R_NegInf;
          return log_phi_std(x) + std::log(sw) + e;
        }, breaks, n_breaks);
    total += piece;
  }
  if (!(normalizer > 0.0) || !R_FINITE(total)) {
    out.status = BAwLCorrMomentStatus::unstable; return out;
  }
  const double value = total / normalizer;
  if (!R_FINITE(value) || value < 0.0 || value > 1.0 + 1e-8) {
    out.status = BAwLCorrMomentStatus::unstable; return out;
  }
  out.value = std::fmin(1.0, std::fmax(0.0, value));
  out.status = (out.value > 0.0)
    ? BAwLCorrMomentStatus::ok : BAwLCorrMomentStatus::zero_mass;
  return out;
}

inline BAwLCorrPairResult bawl_corr_pair_cause_numeric(
    const BAwLTimeGeometry& gw, const BAwLTimeGeometry& gl,
    double muw, double sdw, double mul, double sdl, double rho,
    bool positive, double normalizer) {
  BAwLCorrPairResult out;
  double d = positive ? 0.0 : R_NegInf;
  if (gw.status == BAwLTimeStatus::point_start) {
    const double v = gw.U;
    if (!(v > d)) { out.status = BAwLCorrMomentStatus::zero_mass; return out; }
    const double cm = mul + rho * (sdl / sdw) * (v - muw);
    const double cs = sdl * std::sqrt(std::fmax(1.0 - rho * rho, 1e-24));
    const double e = bawl_corr_conditional_survivor(gl, positive, cm, cs, v);
    out.value = bawl_corr_marginal_density(v, muw, sdw) * gw.dv_dt * e / normalizer;
    out.status = out.value > 0.0 ? BAwLCorrMomentStatus::ok
                                 : BAwLCorrMomentStatus::zero_mass;
    return out;
  }
  BAwLCorrRegion wr;
  if (!bawl_corr_cause_interval(gw, positive, wr)) {
    out.status = BAwLCorrMomentStatus::zero_mass; return out;
  }
  double breaks[kBawlCorrNumericMaxBreaks];
  const int n_breaks = bawl_corr_numeric_breaks(gl, positive, mul, sdl, rho,
                                                breaks);
  const double value = bawl_corr_numeric_integral(
      muw, sdw, wr.lo, wr.hi,
      [&](double v, double x) {
        const double cm = mul + rho * (sdl / sdw) * (v - muw);
        const double cs = sdl * std::sqrt(std::fmax(1.0 - rho * rho, 1e-24));
        const double e = bawl_corr_log_conditional_survivor(gl, positive, cm, cs);
        const double h = gw.gamma0 + gw.gamma1 * v;
        return (h > 0.0 && e > R_NegInf)
          ? log_phi_std(x) + std::log(h) + e : R_NegInf;
      }, breaks, n_breaks);
  out.value = value / normalizer;
  out.status = out.value > 0.0 ? BAwLCorrMomentStatus::ok
                               : BAwLCorrMomentStatus::zero_mass;
  return out;
}

inline double bawl_corr_pair_positive_normalizer(double mu1, double sd1,
                                                  double mu2, double sd2,
                                                  double rho, bool positive) {
  if (!positive) return 1.0;
  if (std::fabs(rho) >= 1.0 - 1e-10) {
    // The exact rectangle derivative layer intentionally declines this
    // nearly singular regime.  Use the continuous degenerate-Gaussian limit
    // for its normalizer so the pair can still take the numeric route rather
    // than falling through immediately as an invalid trial.
    const double a = -mu1 / sd1;
    if (rho >= 0.0) {
      const double lower = std::fmax(a, -mu2 / sd2);
      return pnorm_std(lower, false, false);
    }
    const double upper = mu2 / sd2;
    if (!(upper > a)) return 0.0;
    const double lp = log_normal_interval(a, upper);
    return R_FINITE(lp) ? std::exp(lp) : 0.0;
  }
  // The pair normalizer is just the positive orthant probability.  Do not
  // build the full rectangle-moment grid here: its first moment derivatives
  // are needed by the exact survivor/cause integrals, but not by Z.
  const double p = norm_cdf_2d_hybrid(mu1 / sd1, mu2 / sd2, rho);
  return (R_FINITE(p) && p > 0.0) ? p : R_NaN;
}

#endif

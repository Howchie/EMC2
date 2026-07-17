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
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
};

struct BvnBoundaryGrid {
  // Four boundaries are enough for the pair survivor and pair cause grids;
  // unused entries are ignored.  The arrays also make the common 3x3 case
  // allocation-free.
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  double y[4] = {0.0, 0.0, 0.0, 0.0};
  double cdf[4][4] = {{0.0}};
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

inline void bawl_corr_kahan_add(double value, double& sum, double& correction) {
  const double y = value - correction;
  const double t = sum + y;
  correction = (t - sum) - y;
  sum = t;
}

struct BvnCornerValues {
  double cdf = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dxx = 0.0;
  double dyy = 0.0;
  double dxy = 0.0;
  bool valid = true;
};

inline BvnCornerValues bawl_corr_bvn_corner(double x, double y, double rho) {
  BvnCornerValues out;
  if (x == R_NegInf || y == R_NegInf) return out;
  if (x == R_PosInf && y == R_PosInf) {
    out.cdf = 1.0;
    return out;
  }
  if (x == R_PosInf) {
    out.cdf = pnorm_std(y, true, false);
    out.dy = dnormP(y);
    out.dyy = -y * out.dy;
    return out;
  }
  if (y == R_PosInf) {
    out.cdf = pnorm_std(x, true, false);
    out.dx = dnormP(x);
    out.dxx = -x * out.dx;
    return out;
  }
  if (!R_FINITE(x) || !R_FINITE(y) || !R_FINITE(rho)) {
    out.valid = false;
    return out;
  }

  const double r = std::fmax(-1.0 + 1e-12, std::fmin(1.0 - 1e-12, rho));
  const double one_minus_r2 = std::fmax(1.0 - r * r, 1e-24);
  const double sr = std::sqrt(one_minus_r2);
  if (std::fabs(r) <= 1e-15) {
    const double px = pnorm_std(x, true, false);
    const double py = pnorm_std(y, true, false);
    const double fx = dnormP(x);
    const double fy = dnormP(y);
    out.cdf = px * py;
    out.dx = fx * py;
    out.dy = fy * px;
    out.dxx = -x * fx * py;
    out.dyy = -y * fy * px;
    out.dxy = fx * fy;
    return out;
  }

  out.cdf = norm_cdf_2d(x, y, r);
  const double qx = (y - r * x) / sr;
  const double qy = (x - r * y) / sr;
  const double fx = dnormP(x);
  const double fy = dnormP(y);
  const double fxq = dnormP(qx);
  const double fyq = dnormP(qy);
  out.dx = fx * pnorm_std(qx, true, false);
  out.dy = fy * pnorm_std(qy, true, false);
  // d/dy of d/dx Phi2 is phi(x) phi((y-rho*x)/s) / s.  Using the
  // symmetric conditional argument here is tempting, but it is paired with
  // the wrong marginal density away from the diagonal and corrupts M12 on
  // narrow asymmetric rectangles.
  out.dxy = fx * fxq / sr;
  out.dxx = fx * (-x * pnorm_std(qx, true, false) - r * fxq / sr);
  out.dyy = fy * (-y * pnorm_std(qy, true, false) - r * fyq / sr);
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
    const std::vector<double>& x_bounds,
    const std::vector<double>& y_bounds) {
  BvnBoundaryGrid g;
  g.mu1 = mu1; g.mu2 = mu2; g.sd1 = sd1; g.sd2 = sd2; g.rho = rho;
  g.status = BAwLCorrMomentStatus::ok;
  for (double value : x_bounds) {
    if (!bawl_corr_add_boundary(g.x, g.nx, value)) {
      g.status = BAwLCorrMomentStatus::invalid;
      return g;
    }
  }
  for (double value : y_bounds) {
    if (!bawl_corr_add_boundary(g.y, g.ny, value)) {
      g.status = BAwLCorrMomentStatus::invalid;
      return g;
    }
  }
  for (uint8_t i = 0; i < g.nx; ++i) {
    const double xs = (g.x[i] == R_NegInf || g.x[i] == R_PosInf)
      ? g.x[i] : (g.x[i] - mu1) / sd1;
    for (uint8_t j = 0; j < g.ny; ++j) {
      const double ys = (g.y[j] == R_NegInf || g.y[j] == R_PosInf)
        ? g.y[j] : (g.y[j] - mu2) / sd2;
      const BvnCornerValues corner = bawl_corr_bvn_corner(xs, ys, rho);
      if (!corner.valid) g.status = BAwLCorrMomentStatus::unstable;
      g.cdf[i][j] = corner.cdf;
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
  BvnBoundaryGrid grid = bawl_corr_make_boundary_grid(
      mu1, sd1, mu2, sd2, rho, {lo1, hi1}, {lo2, hi2});
  return bawl_corr_rect_from_grid(grid, lo1, hi1, lo2, hi2);
}

struct BAwLCorrRegion {
  double lo = 0.0;
  double hi = 0.0;
  double c0 = 0.0;
  double c1 = 0.0;
};

inline bool bawl_corr_domain_lower(bool positive, double& d) {
  d = positive ? 0.0 : R_NegInf;
  return true;
}

inline std::vector<BAwLCorrRegion> bawl_corr_survivor_regions(
    const BAwLTimeGeometry& g, bool positive) {
  std::vector<BAwLCorrRegion> out;
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
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
};

inline bool bawl_corr_add_rect_terms(const BvnRectMoments& m,
                                     double a0, double a1,
                                     double b0, double b1,
                                     double& sum, double& correction) {
  if (m.status == BAwLCorrMomentStatus::zero_mass) return true;
  if (m.status != BAwLCorrMomentStatus::ok) return false;
  const double term = a0 * b0 * m.p + a1 * b0 * m.m1 +
                      a0 * b1 * m.m2 + a1 * b1 * m.m12;
  if (!R_FINITE(term)) return false;
  bawl_corr_kahan_add(term, sum, correction);
  return true;
}

inline BAwLCorrPairResult bawl_corr_pair_survival_exact(
    const BAwLTimeGeometry& g1, const BAwLTimeGeometry& g2,
    double mu1, double sd1, double mu2, double sd2, double rho,
    bool positive, double normalizer) {
  BAwLCorrPairResult out;
  const auto r1 = bawl_corr_survivor_regions(g1, positive);
  const auto r2 = bawl_corr_survivor_regions(g2, positive);
  if (r1.empty() || r2.empty()) {
    out.status = BAwLCorrMomentStatus::zero_mass;
    return out;
  }
  std::vector<double> x, y;
  for (const auto& r : r1) { x.push_back(r.lo); x.push_back(r.hi); }
  for (const auto& r : r2) { y.push_back(r.lo); y.push_back(r.hi); }
  BvnBoundaryGrid grid = bawl_corr_make_boundary_grid(
      mu1, sd1, mu2, sd2, rho, x, y);
  if (grid.status == BAwLCorrMomentStatus::invalid) {
    out.status = grid.status; return out;
  }
  double sum = 0.0, correction = 0.0;
  for (const auto& a : r1) for (const auto& b : r2) {
    const BvnRectMoments m = bawl_corr_rect_grid_for_regions(grid, a, b);
    if (!bawl_corr_add_rect_terms(m, a.c0, a.c1, b.c0, b.c1,
                                  sum, correction)) {
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
    bool positive, double normalizer) {
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
  std::vector<double> x = {wr.lo, wr.hi};
  std::vector<double> y;
  for (const auto& r : lr) { y.push_back(r.lo); y.push_back(r.hi); }
  BvnBoundaryGrid grid = bawl_corr_make_boundary_grid(
      muw, sdw, mul, sdl, rho, x, y);
  if (grid.status == BAwLCorrMomentStatus::invalid) {
    out.status = grid.status; return out;
  }
  double sum = 0.0, correction = 0.0;
  for (const auto& r : lr) {
    const BvnRectMoments m = bawl_corr_rect_from_grid(
        grid, wr.lo, wr.hi, r.lo, r.hi);
    if (!bawl_corr_add_rect_terms(m, wr.c0, wr.c1, r.c0, r.c1,
                                  sum, correction)) {
      out.status = BAwLCorrMomentStatus::unstable; return out;
    }
  }
  if (!(normalizer > 0.0) || !R_FINITE(sum)) {
    out.status = BAwLCorrMomentStatus::invalid; return out;
  }
  if (sum < -1e-9) { out.status = BAwLCorrMomentStatus::unstable; return out; }
  out.value = std::fmax(0.0, sum) / normalizer;
  out.status = out.value > 0.0 ? BAwLCorrMomentStatus::ok
                               : BAwLCorrMomentStatus::zero_mass;
  return out;
}

// The deterministic one-dimensional route is deliberately independent of the
// rectangle derivative algebra.  It is used both as the unstable-pair route
// and as the scalar oracle in tests.  Integrating on [-12, 12] in marginal
// standard-deviation units leaves less than 1e-32 of Gaussian mass outside the
// window; all production pair bounds are finite except for the survivor's
// lower/upper normal tail, which is represented by this window.
template <typename LogIntegrand>
inline double bawl_corr_numeric_integral(double mu, double sd, double lo,
                                         double hi, LogIntegrand log_integrand) {
  double zl = (lo == R_NegInf) ? -12.0 : (lo - mu) / sd;
  double zh = (hi == R_PosInf) ? 12.0 : (hi - mu) / sd;
  zl = std::fmax(zl, -12.0);
  zh = std::fmin(zh, 12.0);
  if (!(zh > zl)) return 0.0;
  const GLRule& rule = gl_get_rule(64);
  double max_log = R_NegInf;
  std::array<double, 64> logs{};
  int n = 0;
  for (int i = 0; i < 64; ++i) {
    const double x = 0.5 * (zh - zl) * rule.x[i] + 0.5 * (zh + zl);
    const double lw = std::log(0.5 * (zh - zl) * rule.w[i]);
    const double value = log_integrand(mu + sd * x, x);
    logs[static_cast<size_t>(n++)] = lw + value;
    max_log = std::max(max_log, logs[static_cast<size_t>(n - 1)]);
  }
  if (!R_FINITE(max_log)) return 0.0;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += std::exp(logs[static_cast<size_t>(i)] - max_log);
  return std::exp(max_log) * sum;
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
  const auto wr = bawl_corr_survivor_regions(gw, positive);
  if (wr.empty()) { out.status = BAwLCorrMomentStatus::zero_mass; return out; }
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
        });
    total += piece;
  }
  if (!(normalizer > 0.0) || !R_FINITE(total)) {
    out.status = BAwLCorrMomentStatus::unstable; return out;
  }
  out.value = total / normalizer;
  out.status = (out.value > 0.0 && R_FINITE(out.value))
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
  const double value = bawl_corr_numeric_integral(
      muw, sdw, wr.lo, wr.hi,
      [&](double v, double x) {
        const double cm = mul + rho * (sdl / sdw) * (v - muw);
        const double cs = sdl * std::sqrt(std::fmax(1.0 - rho * rho, 1e-24));
        const double e = bawl_corr_log_conditional_survivor(gl, positive, cm, cs);
        const double h = gw.gamma0 + gw.gamma1 * v;
        return (h > 0.0 && e > R_NegInf)
          ? log_phi_std(x) + std::log(h) + e : R_NegInf;
      });
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
  const double p = norm_cdf_2d(mu1 / sd1, mu2 / sd2, rho);
  return (R_FINITE(p) && p > 0.0) ? p : R_NaN;
}

#endif

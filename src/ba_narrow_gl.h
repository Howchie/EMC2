#ifndef EMC2_BA_NARROW_GL_H
#define EMC2_BA_NARROW_GL_H

// ---------------------------------------------------------------------------
// Narrow-interval fallback shared by the ballistic launch kernels.
//
// Every ballistic model reduces its live start-point average to an integral
// over launch strengths w in [w_lo, w_hi], with w affine in the start point.
// The closed forms difference stop-losses, puts or partial moments at the two
// ends of that interval, and lose digits in proportion to how little the
// launch density varies across it; a density bracket of the form
// E[V; I] - c P(I) additionally cancels by the ratio of w to its small
// crossing weight (about exp(k t) late in a leaky race).  Where log g varies
// by less than BA_NARROW_VAR across the interval, the launch CDF, survivor and
// density are all smooth on it, so the average is evaluated directly with a
// Gauss-Legendre rule, whose error there is far below that cancellation.
// ---------------------------------------------------------------------------

#include <cmath>
#include <Rcpp.h>
#include "gl_quad.h"
#include "utility_functions.h"

constexpr double BA_NARROW_VAR = 0.5;
constexpr int BA_NARROW_GL_N = 6;

// Variation of log g across [w_lo, w_hi], estimated from its values at the
// ends and midpoint.  Twice the larger half-interval change: exact for a
// linear log density, and conservative for the unimodal launch densities.
template <class LG>
inline double ba_log_density_variation(double w_lo, double w_hi, LG lg) {
  if (!(w_hi > w_lo)) return 0.0;
  const double a = lg(w_lo), m = lg(0.5 * (w_lo + w_hi)), c = lg(w_hi);
  if (!emc2_isfinite(a) || !emc2_isfinite(m) || !emc2_isfinite(c))
    return R_PosInf;
  return 2.0 * std::fmax(std::fabs(m - a), std::fabs(c - m));
}

template <class LG>
inline bool ba_narrow_interval(double w_lo, double w_hi, LG lg) {
  return ba_log_density_variation(w_lo, w_hi, lg) < BA_NARROW_VAR;
}

// log of the average over u in [0, 1] of exp(lf(w(u), u)), with
// w(u) = w_hi - u (w_hi - w_lo); u is the start point as a fraction of its
// range, so a model's start-point weight can be formed from it exactly.
template <class LF>
inline double ba_log_mean_gl(double w_lo, double w_hi, LF lf,
                             int n = BA_NARROW_GL_N) {
  const GLRule &rule = gl_get_rule(n);
  double l[32];
  double M = R_NegInf;
  for (int i = 0; i < n; ++i) {
    const double u = 0.5 * (1.0 + rule.x[i]);
    l[i] = std::log(0.5 * rule.w[i]) + lf(w_hi - u * (w_hi - w_lo), u);
    if (l[i] > M) M = l[i];
  }
  if (!(M > R_NegInf) || !emc2_isfinite(M)) return M;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += std::exp(l[i] - M);
  return M + std::log(sum);
}

#endif  // EMC2_BA_NARROW_GL_H

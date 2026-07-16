#ifndef EMC2_ADAPTIVE_GK_H
#define EMC2_ADAPTIVE_GK_H

// Shared-node adaptive Gauss--Kronrod integration for positive vector-valued
// integrands.  The rule below is the 21-point GSL rule (qk21.c).  Unlike
// gsl_integration_qag, this wrapper refines one common interval partition for
// all targets.  That matters for BAwLcorr: every target is a trial, and the
// conditional likelihood can be evaluated for all trials in one pass at a
// latent-factor node.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

struct AdaptiveLogVectorResult {
  std::vector<double> log_integral;
  std::vector<double> log_error;
  std::size_t intervals = 0;
  std::size_t eval_points = 0;
  bool converged = false;
};

inline double adaptive_log_add_exp(double a, double b) {
  if (a == -std::numeric_limits<double>::infinity()) return b;
  if (b == -std::numeric_limits<double>::infinity()) return a;
  if (!std::isfinite(a) || !std::isfinite(b)) return std::max(a, b);
  const double hi = std::max(a, b);
  return hi + std::log1p(std::exp(std::min(a, b) - hi));
}

inline double adaptive_log_sum_exp(const std::vector<double>& x) {
  double out = -std::numeric_limits<double>::infinity();
  for (double v : x) out = adaptive_log_add_exp(out, v);
  return out;
}

namespace emc2_adaptive_gk_detail {

// GSL integration/qk21.c, evaluated with 80-decimal arithmetic by Fullerton.
static constexpr double xgk[11] = {
  0.995657163025808080735527280689003,
  0.973906528517171720077964012084452,
  0.930157491355708226001207180059508,
  0.865063366688984510732096688423493,
  0.780817726586416897063717578345042,
  0.679409568299024406234327365114874,
  0.562757134668604683339000099272694,
  0.433395394129247190799265943165784,
  0.294392862701460198131126603103866,
  0.148874338981631210884826001129720,
  0.0
};

static constexpr double wg[5] = {
  0.066671344308688137593568809893332,
  0.149451349150580593145776339657697,
  0.219086362515982043995534934228163,
  0.269266719309996355091226921569469,
  0.295524224714752870173892994651338
};

static constexpr double wgk[11] = {
  0.011694638867371874278064396062192,
  0.032558162307964727478818972459390,
  0.054755896574351996031381300244580,
  0.075039674810919952767043140916190,
  0.093125454583697605535065465083366,
  0.109387158802297641899210590325805,
  0.123491976262065851077958109831074,
  0.134709217311473325928054001771707,
  0.142775938577060080797094273138717,
  0.147739104901338491374841515972068,
  0.149445554002916905664936468389821
};

struct Interval {
  double a;
  double b;
  std::vector<double> log_result;
  std::vector<double> log_error;
  double priority = -std::numeric_limits<double>::infinity();
};

template <typename Evaluator>
inline Interval evaluate_interval(double a, double b, std::size_t n_targets,
                                  Evaluator& evaluator) {
  Interval out{a, b, std::vector<double>(n_targets,
                                          -std::numeric_limits<double>::infinity()),
                std::vector<double>(n_targets,
                                    -std::numeric_limits<double>::infinity())};
  if (!(b > a)) return out;

  const double center = 0.5 * (a + b);
  const double half = 0.5 * (b - a);
  std::vector<double> values(21 * n_targets,
                             -std::numeric_limits<double>::infinity());
  std::vector<double> at_x(n_targets);
  std::vector<double> plus(n_targets);
  auto evaluate_at = [&](std::size_t node, double x) {
    evaluator(x, at_x);
    for (std::size_t j = 0; j < n_targets; ++j) {
      values[node * n_targets + j] = std::isfinite(at_x[j]) ? at_x[j]
                                                            : -std::numeric_limits<double>::infinity();
    }
  };

  // Node 20 is the center.  Odd Kronrod nodes also carry Gauss weights;
  // even nodes extend the Gauss rule.
  evaluate_at(20, center);
  for (int i = 0; i < 5; ++i) {
    const int odd = 2 * i + 1;
    const int even = 2 * i;
    const double d_odd = half * xgk[odd];
    const double d_even = half * xgk[even];
    evaluate_at(static_cast<std::size_t>(odd), center - d_odd);
    // Store the positive-side value in the paired slot after the negative
    // side.  The exact ordering is private to this rule evaluation.
    evaluator(center + d_odd, plus);
    for (std::size_t j = 0; j < n_targets; ++j)
      values[(10 + i) * n_targets + j] = std::isfinite(plus[j]) ? plus[j]
                                                                  : -std::numeric_limits<double>::infinity();
    evaluate_at(static_cast<std::size_t>(even), center - d_even);
    evaluator(center + d_even, plus);
    for (std::size_t j = 0; j < n_targets; ++j)
      values[(15 + i) * n_targets + j] = std::isfinite(plus[j]) ? plus[j]
                                                                  : -std::numeric_limits<double>::infinity();
  }

  // The first five slots contain the negative odd/even evaluations, slots
  // 10--14 and 15--19 the corresponding positive evaluations.  Recombine
  // them using a local scale so narrow likelihoods do not underflow.
  for (std::size_t j = 0; j < n_targets; ++j) {
    double scale = -std::numeric_limits<double>::infinity();
    for (std::size_t node = 0; node < 21; ++node)
      scale = std::max(scale, values[node * n_targets + j]);
    if (!std::isfinite(scale)) continue;

    double kronrod = wgk[10] * std::exp(values[20 * n_targets + j] - scale);
    double gauss = 0.0;
    for (int i = 0; i < 5; ++i) {
      const int odd = 2 * i + 1;
      const int even = 2 * i;
      const double neg_odd = std::exp(values[odd * n_targets + j] - scale);
      const double pos_odd = std::exp(values[(10 + i) * n_targets + j] - scale);
      const double neg_even = std::exp(values[even * n_targets + j] - scale);
      const double pos_even = std::exp(values[(15 + i) * n_targets + j] - scale);
      kronrod += wgk[odd] * (neg_odd + pos_odd);
      kronrod += wgk[even] * (neg_even + pos_even);
      gauss += wg[i] * (neg_odd + pos_odd);
    }
    const double result = half * kronrod;
    const double error = half * std::fabs(kronrod - gauss);
    if (result > 0.0) out.log_result[j] = scale + std::log(result);
    if (error > 0.0) out.log_error[j] = scale + std::log(error);
  }

  for (std::size_t j = 0; j < n_targets; ++j) {
    if (std::isfinite(out.log_error[j]) && std::isfinite(out.log_result[j])) {
      out.priority = std::max(out.priority, out.log_error[j] - out.log_result[j]);
    }
  }
  return out;
}

} // namespace emc2_adaptive_gk_detail

// Integrate n_targets nonnegative functions supplied in log-space over
// x in [0,1]. The caller maps x to its preferred latent variable.
//
// `initial_breaks` (interior points of (0,1), sorted) seeds the partition.
// A narrow peak whose neighbourhood evaluates to zero at every node of a
// coarse rule is invisible to the error estimator, so the caller should
// supply breaks that keep node spacing below the narrowest feature it needs
// resolved; adaptivity only sharpens what at least one node has seen. `limit`
// is an interval count, not an evaluation count: a split adds two 21-point
// evaluations (42 shared callback points) for every target in the vector.
template <typename Evaluator>
inline AdaptiveLogVectorResult adaptive_gk21_log_vector(
    std::size_t n_targets, Evaluator evaluator,
    double abs_tol = 1e-12, double rel_tol = 1e-8, std::size_t limit = 128,
    const std::vector<double>& initial_breaks = {}) {
  using emc2_adaptive_gk_detail::Interval;
  if (!(abs_tol > 0.0) || !(rel_tol > 0.0))
    throw std::invalid_argument("adaptive GK tolerances must be positive");

  AdaptiveLogVectorResult out;
  out.log_integral.assign(n_targets, -std::numeric_limits<double>::infinity());
  out.log_error.assign(n_targets, -std::numeric_limits<double>::infinity());
  if (n_targets == 0) {
    out.converged = true;
    return out;
  }

  const double eps = 1e-14;
  std::vector<double> edges;
  edges.reserve(initial_breaks.size() + 2);
  edges.push_back(eps);
  for (double b : initial_breaks) {
    if (b > edges.back() && b < 1.0 - eps) edges.push_back(b);
  }
  edges.push_back(1.0 - eps);
  std::vector<Interval> intervals;
  intervals.reserve(edges.size() - 1);
  for (std::size_t e = 0; e + 1 < edges.size(); ++e) {
    intervals.emplace_back(emc2_adaptive_gk_detail::evaluate_interval(
        edges[e], edges[e + 1], n_targets, evaluator));
  }
  out.eval_points = 21 * intervals.size();

  auto recompute_totals = [&]() {
    for (std::size_t j = 0; j < n_targets; ++j) {
      out.log_integral[j] = -std::numeric_limits<double>::infinity();
      out.log_error[j] = -std::numeric_limits<double>::infinity();
      for (const Interval& interval : intervals) {
        out.log_integral[j] = adaptive_log_add_exp(out.log_integral[j],
                                                   interval.log_result[j]);
        out.log_error[j] = adaptive_log_add_exp(out.log_error[j],
                                                interval.log_error[j]);
      }
    }
  };

  auto converged = [&]() {
    for (std::size_t j = 0; j < n_targets; ++j) {
      if (!std::isfinite(out.log_integral[j])) continue;
      const double log_abs = std::log(abs_tol);
      const double log_rel = std::log(rel_tol) + out.log_integral[j];
      const double log_tol = adaptive_log_add_exp(log_abs, log_rel);
      if (out.log_error[j] > log_tol) return false;
    }
    return true;
  };

  recompute_totals();
  while (!converged() && intervals.size() < limit) {
    std::size_t split = 0;
    double priority = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < intervals.size(); ++i) {
      if (intervals[i].priority > priority) {
        priority = intervals[i].priority;
        split = i;
      }
    }
    if (!std::isfinite(priority)) break;
    const double a = intervals[split].a;
    const double b = intervals[split].b;
    const double mid = 0.5 * (a + b);
    if (!(mid > a && mid < b)) break;
    Interval left = emc2_adaptive_gk_detail::evaluate_interval(
        a, mid, n_targets, evaluator);
    Interval right = emc2_adaptive_gk_detail::evaluate_interval(
        mid, b, n_targets, evaluator);
    intervals[split] = std::move(left);
    intervals.emplace_back(std::move(right));
    out.eval_points += 42;
    recompute_totals();
  }
  out.intervals = intervals.size();
  out.converged = converged();
  return out;
}

#endif

#ifndef EMC2_QUAD_TEMPLATES_H
#define EMC2_QUAD_TEMPLATES_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <utility>

#include "gl_quad.h"
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "utility_functions.h"
#include "composite_functions.h"

// Log-space counterpart of integrate_positive_drift_quad below: accumulates
// log(w_j) + log_kernel(drift_j) with log_sum_exp, so node under/overflow in
// the natural accumulation cannot zero out a representable mixture.  Only
// entered after the natural quadrature has been rejected; log_out remains a
// pure output-scale choice at the call sites.
template <typename LogKernelFn>
inline double integrate_positive_drift_quad_log(double mu_drift, double sv,
                                                LogKernelFn&& log_kernel_fn,
                                                int n_gauss_nodes = 20) {
  if (!(sv > 1e-10) || !emc2_isfinite(sv)) {
    return log_kernel_fn(mu_drift);
  }

  const double lower_p_raw = pnorm_std(-mu_drift / sv, true, false);
  const double upper_p = std::nextafter(1.0, 0.0);
  const double lower_p = std::fmax(0.0, std::fmin(lower_p_raw, upper_p));
  const double width = upper_p - lower_p;
  if (!(width > 0.0)) return R_NegInf;

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  double log_acc = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    if (!(gl.w[j] > 0.0)) continue;
    const double p = lower_p + 0.5 * width * (gl.x[j] + 1.0);
    const double drift = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    const double lk = log_kernel_fn(drift);
    if (lk == R_NegInf || ISNAN(lk)) continue;
    log_acc = log_sum_exp(log_acc, std::log(gl.w[j]) + lk);
  }
  return log_acc - M_LN2;
}

template <typename KernelFn>
inline double integrate_positive_drift_quad(double mu_drift, double sv,
                                            KernelFn&& kernel_fn,
                                            int n_gauss_nodes = 20) {
  if (!(sv > 1e-10) || !emc2_isfinite(sv)) {
    return kernel_fn(mu_drift);
  }

  const double lower_p_raw = pnorm_std(-mu_drift / sv, true, false);
  const double upper_p = std::nextafter(1.0, 0.0);
  const double lower_p = std::fmax(0.0, std::fmin(lower_p_raw, upper_p));
  const double width = upper_p - lower_p;
  if (!(width > 0.0)) return 0.0;

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;

  double acc = 0.0;
  for (int j = 0; j < n_nodes; ++j) {
    const double p = lower_p + 0.5 * width * (nodes[j] + 1.0);
    const double drift = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    acc += weights[j] * kernel_fn(drift);
  }

  const double result = 0.5 * acc;
  return emc2_isfinite(result) ? result : 0.0;
}

template <typename DensityFn>
inline double integrate_density_gl20_finite(double t_upper, DensityFn&& density_fn) {
  if (!(t_upper > 0.0)) return 0.0;
  const GLRule& gl = gl_get_rule(20);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;
  double acc = 0.0;
  for (int j = 0; j < static_cast<int>(nodes.size()); ++j) {
    const double u = 0.5 * t_upper * (nodes[j] + 1.0);
    acc += weights[j] * density_fn(u);
  }
  return 0.5 * t_upper * acc;
}

template <typename DensityFn>
inline double integrate_density_gl20_infinite(double rate_scale, DensityFn&& density_fn) {
  const double scale = std::fmax(rate_scale, 1e-8);
  const GLRule& gl = gl_get_rule(20);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;
  double acc = 0.0;
  for (int j = 0; j < static_cast<int>(nodes.size()); ++j) {
    const double q = 0.5 * (nodes[j] + 1.0);
    const double qq = std::fmin(1.0 - 1e-12, std::fmax(1e-15, q));
    const double t = -std::log1p(-qq) / scale;
    const double jac = 1.0 / (scale * (1.0 - qq));
    acc += weights[j] * density_fn(t) * jac;
  }
  return 0.5 * acc;
}

template <typename DensityFn>
inline double integrate_density_adaptive_finite(double t_upper, DensityFn&& density_fn) {
  if (!(t_upper > 0.0)) return 0.0;

  using Fn = std::decay_t<DensityFn>;
  Fn fn = std::forward<DensityFn>(density_fn);
  struct Adapter {
    const Fn* fn;
  } adapter{&fn};

  gsl_function F;
  F.function = +[](double x, void* p) -> double {
    const auto* adapter_ptr = static_cast<const Adapter*>(p);
    return (*(adapter_ptr->fn))(x);
  };
  F.params = &adapter;

  GslIntegrationControls ctl = default_gsl_controls();
  ctl.try_qng_first_finite = true;
  ctl.qag_key = GSL_INTEG_GAUSS21;
  ctl.rel_tol = 1e-6;

  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws, ctl.retry_limit);

  double result = 0.0;
  double err = 0.0;
  int status = GSL_EFAILED;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

  size_t neval = 0;
  status = gsl_integration_qng(&F, 0.0, t_upper, ctl.abs_tol, ctl.rel_tol,
                               &result, &err, &neval);
  if (status != GSL_SUCCESS || !emc2_isfinite(result)) {
    status = gsl_integration_qag(&F, 0.0, t_upper,
                                 ctl.abs_tol, ctl.rel_tol,
                                 ctl.limit, ctl.qag_key,
                                 workspace, &result, &err);
  }
  if (status != GSL_SUCCESS || !emc2_isfinite(result)) {
    status = gsl_integration_qags(&F, 0.0, t_upper,
                                  ctl.retry_abs_tol, ctl.retry_rel_tol,
                                  ctl.retry_limit, workspace, &result, &err);
  }

  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS || !emc2_isfinite(result) || result < 0.0) {
    return integrate_density_gl20_finite(t_upper, fn);
  }
  return result;
}

template <typename DensityFn>
inline double integrate_density_adaptive_infinite(double rate_scale, DensityFn&& density_fn) {
  using Fn = std::decay_t<DensityFn>;
  Fn fn = std::forward<DensityFn>(density_fn);
  struct Adapter {
    const Fn* fn;
  } adapter{&fn};

  gsl_function F;
  F.function = +[](double x, void* p) -> double {
    const auto* adapter_ptr = static_cast<const Adapter*>(p);
    return (*(adapter_ptr->fn))(x);
  };
  F.params = &adapter;

  GslIntegrationControls ctl = default_gsl_controls();
  ctl.rel_tol = 1e-6;
  if (rate_scale > 1e-8 && emc2_isfinite(rate_scale)) {
    ctl.abs_tol = std::min(ctl.abs_tol, 1e-10 / rate_scale);
  }

  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws, ctl.retry_limit);

  double result = 0.0;
  double err = 0.0;
  int status = GSL_EFAILED;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

  status = gsl_integration_qagiu(&F, 0.0, ctl.abs_tol, ctl.rel_tol,
                                 ctl.limit, workspace, &result, &err);
  if (status != GSL_SUCCESS || !emc2_isfinite(result)) {
    status = gsl_integration_qagiu(&F, 0.0, ctl.retry_abs_tol, ctl.retry_rel_tol,
                                   ctl.retry_limit, workspace, &result, &err);
  }

  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS || !emc2_isfinite(result) || result < 0.0) {
    return integrate_density_gl20_infinite(rate_scale, fn);
  }
  return result;
}

constexpr int BAWD_GL_NODES = 24;      // per panel, frozen integral
constexpr int BAWD_GL_MAX_NODES = 64;

template <typename LogFun>
inline double bawd_log_gl(const LogFun& log_f, double a, double b, int n) {
  if (!(b > a) || !emc2_isfinite(a) || !emc2_isfinite(b)) return R_NegInf;
  if (n > BAWD_GL_MAX_NODES) n = BAWD_GL_MAX_NODES;
  const GLRule& r = gl_get_rule(n);
  const double c1 = 0.5 * (b - a), c2 = 0.5 * (b + a);
  double terms[BAWD_GL_MAX_NODES];
  double best = R_NegInf;
  for (int i = 0; i < n; ++i) {
    const double lf = log_f(c1 * r.x[i] + c2);
    terms[i] = (r.w[i] > 0.0 && lf > R_NegInf && !ISNAN(lf))
      ? (std::log(r.w[i]) + lf) : R_NegInf;
    if (terms[i] > best) best = terms[i];
  }
  if (!(best > R_NegInf)) return R_NegInf;
  double acc = 0.0;
  for (int i = 0; i < n; ++i)
    if (terms[i] > R_NegInf) acc += std::exp(terms[i] - best);
  return std::log(c1) + best + std::log(acc);
}

// Two panels split at `mid` when it lies strictly inside; the split point is
// where the integrand's normal factor turns over, which is where a single
// fixed rule is least accurate.
template <typename LogFun>
inline double bawd_log_gl_split(const LogFun& log_f, double a, double b,
                                double mid, int n) {
  if (mid > a && mid < b) {
    return log_sum_exp(bawd_log_gl(log_f, a, mid, n),
                       bawd_log_gl(log_f, mid, b, n));
  }
  return bawd_log_gl(log_f, a, b, n);
}

#endif  // EMC2_QUAD_TEMPLATES_H

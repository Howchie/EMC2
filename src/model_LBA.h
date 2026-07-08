#ifndef lba_h
#define lba_h

// This header may be included by exactly ONE translation unit (particle_ll.cpp,
// directly and via utils.h) because it defines [[Rcpp::export]] functions
// (dlba, plba, …) that RcppExports links to. Non-exported free helpers here are
// marked `inline`; the exported ones must NOT be inline.
#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "wald_functions.h"  // pnorm_std() — fast normal CDF under USE_FAST_PNORM
#include "composite_functions.h"  // clamp_pos, safe_log
#include "gaussian.h"
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;

// Route through pnorm_std so USE_FAST_PNORM applies to LBA as well as RDM/Wald.
// pnorm(q, mean, sd) = pnorm_std((q - mean) / sd) for sd > 0.
inline double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log = false){
  return pnorm_std((q - mean) / sd, lower, log);
}

inline double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

inline double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true, bool log_out = false){
  if (t == R_PosInf) {
    const double cdf = posdrift ? 1.0 : pnormP(0.0, v, sv, false, false);
    return log_out ? std::log(cdf) : cdf;
  }

  double denom = 1.;
  if (posdrift) {
    denom = pnormP(v / sv, 0., 1., true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double cdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;
    double xx = cmz - A;
    double cz = cmz / zs;
    double cz_max = xx / zs;
    cdf = (1. + (zs * (dnormP(cz_max, 0., 1., false) - dnormP(cz, 0., 1., false))
                   + xx * pnormP(cz_max, 0., 1., true, false) - cmz * pnormP(cz, 0., 1., true, false))/A) / denom;
  } else {
    cdf = pnormP(b / t, v, sv, false, false) / denom;
  }

  if (cdf < 0.) {
    return log_out ? R_NegInf : 0.0;
  } else if (cdf > 1.){
    return log_out ? 0.0 : 1.0;
  }
  return log_out ? std::log(cdf) : cdf;
}

inline double dlba_norm(double t, double A,double b, double v, double sv,
                 bool posdrift = true, bool log_out = false){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP(v / sv, 0., 1., true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double pdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;;
    double cz = cmz / zs;
    double cz_max = (cmz - A) / zs;
    pdf = (v * (pnormP(cz, 0., 1., true, false) - pnormP(cz_max, 0., 1., true, false)) +
      sv * (dnormP(cz_max, 0., 1., false) - dnormP(cz, 0., 1., false))) / (A * denom);
  } else {
    pdf = dnormP(b / t, v, sv, false) * b / (t * t * denom);
  }

  if (pdf < 0.) {
    return log_out ? R_NegInf : 0.0;
  }
  return log_out ? std::log(pdf) : pdf;
}

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true, bool log_out = false)

{
  int n = t.size();
  NumericVector pdf(n);

  for (int i = 0; i < n; i++){
    pdf[i] = dlba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift, log_out);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true, bool log_out = false)

{
  int n = t.size();
  NumericVector cdf(n);

  for (int i = 0; i < n; i++){
    cdf[i] = plba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift, log_out);
  }
  return cdf;
}

// --------------------------------------------------------------------------
// Ballistic Accumulator with Leak (BAwL)
//
// Accumulator: x(t) = a*exp(-k*t) + (D/k)*(1 - exp(-k*t))
//   where a ~ Unif(0,A), D ~ N(v, sv^2) [optionally truncated D>0].
// Hits threshold b when x(t) >= b.
// k -> 0 limit recovers standard LBA exactly.
// --------------------------------------------------------------------------

// Stable exp(-kt) and 1-exp(-kt).
inline void leak_terms(double k, double t, double &E, double &G) {
  E = std::exp(-k * t);
  G = -std::expm1(-k * t);  // = 1 - exp(-kt), numerically stable for small kt
}

// H(z) = z * Phi(z) + phi(z)
inline double Hfun(double z) {
  return z * pnorm_std(z) + dnormP(z);
}

// CDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double pleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false) {
  const double eps        = 1e-10;
  const double denom_floor = 1e-300;

  if (t <= 0.0 || !(sv > 0.0) || !(b >= A) || !(b > 0.0)) {
    return log_out ? R_NegInf : 0.0;
  }

  double denom = 1.0;
  if (posdrift) {
    denom = pnorm_std(v / sv);
    denom = clamp_pos(denom, denom_floor);
  }

  // k -> 0: reduces to standard LBA
  if (std::fabs(k) < eps) {
    const double c = (v - (b / t)) / sv;
    const double m = (1.0 / t) / sv;
    if (A < eps) {
      double F = pnorm_std(c) / denom;
      return log_out ? safe_log(F) : F;
    }
    if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
      double F = pnorm_std(c + 0.5 * m * A) / denom;
      return log_out ? safe_log(F) : F;
    }
    double zA = c + m * A;
    double F = (Hfun(zA) - Hfun(c)) / (A * m);
    F /= denom;
    if (!(F >= 0.0) || std::isnan(F)) F = 0.0;
    return log_out ? safe_log(F) : F;
  }

  // General k > 0
  double E, G;
  leak_terms(k, t, E, G);
  G = clamp_pos(G, 1e-300);
  const double C1 = (k * b) / G;
  const double C2 = (k * E) / G;
  const double c  = (v - C1) / sv;
  const double m  = C2 / sv;

  if (A < eps) {
    double F = pnorm_std(c) / denom;
    return log_out ? safe_log(F) : F;
  }
  if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
    double F = pnorm_std(c + 0.5 * m * A) / denom;
    return log_out ? safe_log(F) : F;
  }
  double zA = c + m * A;
  double F = (Hfun(zA) - Hfun(c)) / (A * m);
  F /= denom;
  if (!(F >= 0.0) || std::isnan(F)) F = 0.0;
  return log_out ? safe_log(F) : F;
}

// PDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double dleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false) {
  const double eps        = 1e-10;
  const double denom_floor = 1e-300;

  if (t <= 0.0 || !(sv > 0.0) || !(b >= A) || !(b > 0.0)) {
    return log_out ? R_NegInf : 0.0;
  }

  double denom = 1.0;
  if (posdrift) {
    denom = pnorm_std(v / sv);
    denom = clamp_pos(denom, denom_floor);
  }

  // k -> 0: reduces to standard LBA
  if (std::fabs(k) < eps) {
    const double c = (v - (b / t)) / sv;
    const double m = (1.0 / t) / sv;
    const double B = 1.0 / (t * t);

    if (A < eps) {
      double f = B * b * (dnormP(c) / sv) / denom;
      if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
      return log_out ? safe_log(f) : f;
    }
    if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
      double zmid = c + 0.5 * m * A;
      double f = B * (b - 0.5 * A) * (dnormP(zmid) / sv) / denom;
      if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
      return log_out ? safe_log(f) : f;
    }
    double zA = c + m * A;
    double dphi = pnorm_std(zA) - pnorm_std(c);
    double ddnorm = dnormP(zA) - dnormP(c);
    double bracket = (b * m + c) * dphi + ddnorm;
    double f = (B / (A * sv * m * m)) * bracket / denom;
    if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
    return log_out ? safe_log(f) : f;
  }

  // General k > 0
  double E, G;
  leak_terms(k, t, E, G);
  G = clamp_pos(G, 1e-300);
  const double C1 = (k * b) / G;
  const double C2 = (k * E) / G;
  const double c  = (v - C1) / sv;
  const double m  = C2 / sv;
  const double B  = (k * k * E) / (G * G);  // Jacobian base

  if (A < eps) {
    double f = B * b * (dnormP(c) / sv) / denom;
    if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
    return log_out ? safe_log(f) : f;
  }
  if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
    double zmid = c + 0.5 * m * A;
    double f = B * (b - 0.5 * A) * (dnormP(zmid) / sv) / denom;
    if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
    return log_out ? safe_log(f) : f;
  }
  double zA = c + m * A;
  double dphi  = pnorm_std(zA) - pnorm_std(c);
  double ddnorm = dnormP(zA) - dnormP(c);
  double bracket = (b * m + c) * dphi + ddnorm;
  double f = (B / (A * sv * m * m)) * bracket / denom;
  if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
  return log_out ? safe_log(f) : f;
}

// Killed-leaky BA density (hit + guess mixture):
//   f_R(t_eam) * S_K(t) * S_G(t)  +  f_G(t) * S_K(t) * S_R(t_eam)
// hit term = racer density surviving the kill and guess clocks; guess term =
// guess density surviving the kill clock and the racer (S_R = 1 - F_R).
// t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang clocks use raw t.
inline double dkilledleakyba_norm(double t, double v, double b, double A,
                                  double sv, double t0 = 0.0,
                                  double k = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
                                  bool posdrift = true, bool log_out = false,
                                  int kill_shape = 1, bool guess = false, double erlang_omega = 1.0) {
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  const double t_eam = t - t0;
  const bool use_guess = guess && (lambda_g > 0.0);
  const bool use_kill = (lambda_k > 0.0);

  // Erlang survivals always use raw t
  const double log_sK = use_kill  ? erlang_log_surv(t, lambda_k, kill_shape, erlang_omega) : 0.0;
  const double log_sG = use_guess ? erlang_log_surv(t, lambda_g, kill_shape, erlang_omega) : 0.0;

  // When EAM hasn't started: f_EAM = 0, S_EAM = 1.
  if (t_eam <= 0.0) {
    if (!use_guess) return log_out ? R_NegInf : 0.0;
    // Only guess path contributes: f_G(t) * S_K(t) * S_R(0) = f_G(t) * S_K(t) * 1
    const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape, erlang_omega);
    const double log_f_guess = log_fG + log_sK;
    return log_out ? log_f_guess : std::exp(log_f_guess);
  }

  if (!use_guess && !use_kill) return dleakyba_norm(t_eam, A, b, v, sv, k, posdrift, log_out);

  const double log_fR = dleakyba_norm(t_eam, A, b, v, sv, k, posdrift, true);
  const double log_f_hit = log_fR + log_sK + log_sG;
  if (!use_guess) return log_out ? log_f_hit : std::exp(log_f_hit);

  const double cdf_r = pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, false);
  const double log_sr = std::log1p(-std::max(0.0, std::min(1.0, cdf_r)));
  const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape, erlang_omega);
  const double log_f_guess = log_fG + log_sK + log_sr;
  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double integrate_bawl_pdf_raw(double t, double v, double b, double A,
                                     double sv, double t0, double k,
                                     double lambda_g, double lambda_k,
                                     bool posdrift, int kill_shape, bool guess, double erlang_omega = 1.0) {
  const double lower = guess ? 0.0 : t0;
  if (t != R_PosInf && t <= lower) return 0.0;

  std::function<double(double)> fn = [&](double u) -> double {
    return dkilledleakyba_norm(u, v, b, A, sv, t0, k, lambda_g, lambda_k,
                               posdrift, false, kill_shape, guess, erlang_omega);
  };
  gsl_function F;
  F.function = [](double u, void* p) -> double {
    return (*static_cast<std::function<double(double)>*>(p))(u);
  };
  F.params = &fn;

  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(ws);
  double result = 0.0, err = 0.0;
  gsl_error_handler_t* old = gsl_set_error_handler_off();
  int status;
  if (t == R_PosInf) {
    status = gsl_integration_qagiu(&F, lower, 1e-8, 1e-5, 200, w, &result, &err);
  } else {
    status = gsl_integration_qags(&F, lower, t, 1e-8, 1e-5, 200, w, &result, &err);
  }
  gsl_set_error_handler(old);
  if (status != GSL_SUCCESS || !R_FINITE(result)) return 0.0;
  return std::max(0.0, std::min(1.0, result));
}

// Killed-leaky BA sub-CDF:
// P(T_R <= t, T_R < T_K) with kill_shape=1 (exponential) or kill_shape=2 (Erlang-2).
// With guess=true: mixture CDF = 1 - S_R(t_eam)*S_K(t)*S_G(t).
// t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
inline double pkilledleakyba_norm(double t, double v, double b, double A,
                                  double sv, double t0 = 0.0,
                                  double k = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
                                  bool posdrift = true, bool log_out = false,
                                  int kill_shape = 1, bool guess = false, double erlang_omega = 1.0) {
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  const double t_eam = t - t0;
  const bool use_guess = guess && (lambda_g > 0.0);
  const bool use_kill = (lambda_k > 0.0);

  // When EAM hasn't started: P(EAM fires by t_eam) = 0, S_R = 1.
  if (t_eam <= 0.0) {
    if (!use_guess) return log_out ? R_NegInf : 0.0;
    if (!use_kill) {
      // Guess-only response process before EAM onset.
      const double log_sG = erlang_log_surv(t, lambda_g, kill_shape, erlang_omega);
      const double out = -std::expm1(log_sG);
      return log_out ? safe_log(out) : out;
    }
    // Before EAM onset, the only observed response is a guess before the kill clock:
    // integral_0^t f_G(u) S_K(u) du.  Kill wins are omissions, not CDF mass.
    std::function<double(double)> fn = [&](double u) -> double {
      return std::exp(erlang_log_pdf(u, lambda_g, kill_shape, erlang_omega) +
                      erlang_log_surv(u, lambda_k, kill_shape, erlang_omega));
    };
    gsl_function F;
    F.function = [](double u, void* p) -> double {
      return (*static_cast<std::function<double(double)>*>(p))(u);
    };
    F.params = &fn;
    static thread_local GslWorkspacePtr ws_pre_eam(nullptr, &gsl_integration_workspace_free);
    gsl_integration_workspace* w = ensure_gsl_workspace(ws_pre_eam);
    double out = 0.0, err = 0.0;
    gsl_error_handler_t* old = gsl_set_error_handler_off();
    gsl_integration_qags(&F, 0.0, t, 1e-8, 1e-5, 200, w, &out, &err);
    gsl_set_error_handler(old);
    out = std::max(0.0, std::min(1.0, out));
    return log_out ? safe_log(out) : out;
  }

  if (!use_guess && !use_kill) return pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, log_out);

  if (use_guess && use_kill) {
    const double out = integrate_bawl_pdf_raw(
      t, v, b, A, sv, t0, k, lambda_g, lambda_k, posdrift, kill_shape, true, erlang_omega
    );
    return log_out ? safe_log(out) : out;
  }

  const double lambda = use_guess ? lambda_g : lambda_k;
  if (use_guess) {
    // 1 - S_R(t_eam) * S_G(t); erlang uses raw t
    const double cdf_r = pleakyba_norm(t_eam, A, b, v, sv, k, posdrift, false);
    const double log_sr = std::log1p(-std::max(0.0, std::min(1.0, cdf_r)));
    const double log_sg = erlang_log_surv(t, lambda, kill_shape, erlang_omega);
    const double log_val = std::log1p(-std::exp(log_sr + log_sg));
    return log_out ? log_val : std::exp(log_val);
  }
  if (sv <= 0.0 || b < A || b <= 0.0) return log_out ? R_NegInf : 0.0;

  // Pure kill path: integrate f_EAM(u - t0) * S_K(u) over raw time.
  // No closed normal-form primitive remains once leak and clock survival are combined.
  const double out = integrate_bawl_pdf_raw(
    t, v, b, A, sv, t0, k, 0.0, lambda, posdrift, kill_shape, false, erlang_omega
  );
  return log_out ? safe_log(out) : out;
}

// [[Rcpp::export]]
NumericVector dkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    pdf[i] = dkilledleakyba_norm(t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
                                 pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
                                 posdrift, log_out, kill_shape, guess, omega);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    cdf[i] = pkilledleakyba_norm(t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
                                 pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
                                 posdrift, log_out, kill_shape, guess, omega);
  }
  return cdf;
}

// Vectorised R-callable wrappers (recycle scalar parameters).
// [[Rcpp::export]]
NumericVector dleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = dleakyba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i), posdrift);
  return pdf;
}

// [[Rcpp::export]]
NumericVector pleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = pleakyba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i), posdrift);
  return cdf;
}

#endif

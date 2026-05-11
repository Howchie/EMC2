#ifndef lba_h
#define lba_h

#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "wald_functions.h"  // pnorm_std() — fast normal CDF under USE_FAST_PNORM
#include "composite_functions.h"  // clamp_pos, safe_log
#include "gaussian.h"

using namespace Rcpp;

// Route through pnorm_std so USE_FAST_PNORM applies to LBA as well as RDM/Wald.
// pnorm(q, mean, sd) = pnorm_std((q - mean) / sd) for sd > 0.
double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log = false){
  return pnorm_std((q - mean) / sd, lower, log);
}

double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true, bool log_out = false){
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

double dlba_norm(double t, double A,double b, double v, double sv,
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

// Killed-leaky BA density: f_hit(t_eam) * S_K(t) + f_K(t) * S_R(t_eam)
// t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
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
  const Rcpp::List& gl = get_gl20();
  const Rcpp::NumericVector nodes = gl["nodes"];
  const Rcpp::NumericVector weights = gl["weights"];
  double acc = 0.0;

  if (t == R_PosInf) {
    // Clocked defective mass over [0, Inf).  Scale the transform by the
    // active Erlang rates because the clock survival terms set the tail scale.
    const double rate = std::max(1e-8, (guess ? lambda_g : 0.0) + lambda_k);
    for (int j = 0; j < nodes.size(); ++j) {
      const double q = std::min(1.0 - 1e-12, std::max(1e-15, 0.5 * (nodes[j] + 1.0)));
      const double u = -std::log1p(-q) / rate;
      const double jac = 1.0 / (rate * (1.0 - q));
      acc += weights[j] * dkilledleakyba_norm(
        u, v, b, A, sv, t0, k, lambda_g, lambda_k, posdrift, false, kill_shape, guess, erlang_omega
      ) * jac;
    }
    double out = 0.5 * acc;
    return std::max(0.0, std::min(1.0, out));
  }

  for (int j = 0; j < nodes.size(); ++j) {
    const double u = 0.5 * t * (nodes[j] + 1.0);
    acc += weights[j] * dkilledleakyba_norm(
      u, v, b, A, sv, t0, k, lambda_g, lambda_k, posdrift, false, kill_shape, guess, erlang_omega
    );
  }
  double out = 0.5 * t * acc;
  return std::max(0.0, std::min(1.0, out));
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
    const Rcpp::List& gl = get_gl20();
    const Rcpp::NumericVector nodes = gl["nodes"];
    const Rcpp::NumericVector weights = gl["weights"];
    double acc = 0.0;
    for (int j = 0; j < nodes.size(); ++j) {
      const double u = 0.5 * t * (nodes[j] + 1.0);
      const double log_fG = erlang_log_pdf(u, lambda_g, kill_shape, erlang_omega);
      const double log_sK = erlang_log_surv(u, lambda_k, kill_shape, erlang_omega);
      acc += weights[j] * std::exp(log_fG + log_sK);
    }
    double out = 0.5 * t * acc;
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
    const double log_sk = erlang_log_surv(t, lambda, kill_shape, erlang_omega);
    const double log_val = std::log1p(-std::exp(log_sr + log_sk));
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

// --------------------------------------------------------------------------
// Killed LBA: standard LBA hits truncated by an exponential expiry clock
// T_kill ~ Exp(k).
// --------------------------------------------------------------------------

// PDF of killed LBA: f_k(t) = f_lba(t) * exp(-k*t)
// [[Rcpp::export]]
double dkilledlba_norm(double t, double A, double b, double v, double sv, double k,
                       bool posdrift = true, bool log_out = false) {
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  if (k <= 0.0) return dlba_norm(t, A, b, v, sv, posdrift, log_out);

  double log_pdf = dlba_norm(t, A, b, v, sv, posdrift, true) - k * t;
  return log_out ? log_pdf : std::exp(log_pdf);
}

// Sub-CDF of killed LBA: P(T_lba < t AND T_lba < T_kill)
// Uses 20-point Gauss-Legendre quadrature over the drift distribution.
// [[Rcpp::export]]
double pkilledlba_norm(double t, double A, double b, double v, double sv, double k,
                       bool posdrift = true, bool log_out = false) {
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  if (k <= 0.0) return plba_norm(t, A, b, v, sv, posdrift, log_out);

  const double d0 = b - A;
  const double eps = 1e-10;

  auto get_subcdf_fixed_v = [&](double drift) {
    if (drift <= 0.0) return 0.0;
    double t_start = d0 / drift;
    if (t <= t_start) return 0.0;
    double t_end = b / drift;
    double t_upper = std::min(t, t_end);
    if (t_upper <= t_start) return 0.0;
    // Integral_{t_start}^{t_upper} (drift/A) * exp(-k*tau) dtau
    return (drift / (A * k)) * (std::exp(-k * t_start) - std::exp(-k * t_upper));
  };

  double cdf = 0.0;
  if (sv < eps) {
    cdf = get_subcdf_fixed_v(v);
  } else {
    // Integrate over drift distribution
    const Rcpp::List& gl = get_gl20();
    const Rcpp::NumericVector nodes = gl["nodes"];
    const Rcpp::NumericVector weights = gl["weights"];
    const int n_nodes = nodes.size();

    double alpha = pnorm_std(-v / sv); // P(V < 0)

    for (int j = 0; j < n_nodes; j++) {
      double u = 0.5 * (nodes[j] + 1.0);
      double p = alpha + (1.0 - alpha) * u;
      // Clamp p for stability
      p = std::max(1e-15, std::min(1.0 - 1e-15, p));
      double drift_j = v + sv * R::qnorm(p, 0.0, 1.0, true, false);
      cdf += weights[j] * get_subcdf_fixed_v(drift_j);
    }
    cdf *= 0.5;
    if (!posdrift) {
      cdf *= (1.0 - alpha);
    }
    // Note: if posdrift=true, we should NOT divide by (1-alpha) here because
    // the integral itself is over the truncated density.
    // Wait, if posdrift=true, g(r) = phi(r)/ (1-alpha).
    // The integral over p from 0 to 1 gives the integral over the truncated density.
    // Yes, 0.5 * sum(w * f(r)) integrates over the truncated space.
  }

  if (cdf < 0.0) cdf = 0.0;
  if (cdf > 1.0) cdf = 1.0;
  return log_out ? std::log(cdf) : cdf;
}

// [[Rcpp::export]]
NumericVector dkilledlba(NumericVector t,
                         NumericVector A, NumericVector b,
                         NumericVector v, NumericVector sv, NumericVector k,
                         bool posdrift = true, bool log_out = false) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = dkilledlba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i), posdrift, log_out);
  return pdf;
}

// [[Rcpp::export]]
NumericVector pkilledlba(NumericVector t,
                         NumericVector A, NumericVector b,
                         NumericVector v, NumericVector sv, NumericVector k,
                         bool posdrift = true, bool log_out = false) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = pkilledlba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i), posdrift, log_out);
  return cdf;
}

#endif

#ifndef lba_h
#define lba_h

#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "wald_functions.h"  // pnorm_std() — fast normal CDF under USE_FAST_PNORM
#include "composite_functions.h"  // clamp_pos, safe_log

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
                 bool posdrift = true){
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
    return 0.;
  } else if (cdf > 1.){
    return 1.;
  }
  return cdf;
}

double dlba_norm(double t, double A,double b, double v, double sv,
                 bool posdrift = true){
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
    return 0.;
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector pdf(n);

  for (int i = 0; i < n; i++){
    pdf[i] = dlba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector cdf(n);

  for (int i = 0; i < n; i++){
    cdf[i] = plba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift);
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


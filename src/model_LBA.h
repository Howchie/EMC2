#ifndef lba_h
#define lba_h

#include <RcppArmadillo.h>
#include "utility_functions.h"

using namespace Rcpp;

double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log = false){
  return R::pnorm(q, mean, sd, lower, log);
}

double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

// ZH new math for the Ballistic Accumulator with Leak
inline void leak_terms(double k, double t, double &E, double &G) {
  // E = exp(-k t)
  E = std::exp(-k * t);
  // G = 1 - E, but stable for small kt:
  // expm1(x) returns exp(x)-1, so -expm1(-kt) = 1-exp(-kt)
  G = -std::expm1(-k * t);
}
// H(z) = z*Pnorm(z) + dnorm(z)
double Hfun(double z) { return z * pnormP(z) + dnormP(z); }


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

NumericVector dlba_c1(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector plba_c1(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool posdrift = true){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), posdrift);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool posdrift = true){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), posdrift);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
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
// Closed-form leaky 1D ballistic accumulator
// x(t) = a*exp(-k t) + (D/k)*(1-exp(-k t)) , hit when x(t) >= b
// a ~ Unif(0,A), D ~ N(v, sv^2); optionally truncated to D>0.
//
// Defective: F(infty) = P(D >= k b) (or conditional on D>0 if truncated).
// --------------------------------------------------------------------------

// [[Rcpp::export]]
double pleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false) {

  const double eps  = 1e-10;
  const double denom_floor = 1e-300;


  // Truncation normalizer for D>0
  double denom = 1.0;
  if (posdrift) {
    denom = pnormP(v / sv);
    denom = clamp_pos(denom, denom_floor);
  }

  // k = 0 branch (no leak): v_req(a,t) = (b-a)/t
  // this should reduce to the standard lba! 
  if (std::fabs(k) < eps) {
    // z(a,t) = (v - (b-a)/t)/sv = c + m a
    const double c = (v - (b / t)) / sv;
    const double m = (1.0 / t) / sv;

    if (A < eps) {
      double F = pnormP(c) / denom;
      return log_out ? safe_log(F) : F;
    }

    // If slope is tiny, avoid division; use midpoint approx (accurate when m*A small)
    if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
      double F = pnormP(c + 0.5 * m * A) / denom;
      return log_out ? safe_log(F) : F;
    }

    double z0 = c;
    double zA = c + m * A;
    double F = (Hfun(zA) - Hfun(z0)) / (A * m);
    F /= denom;
    if (!(F >= 0.0) || std::isnan(F)) F = 0.0;
    return log_out ? safe_log(F) : F;
  }

  // General k > 0
  double E, G;
  leak_terms(k, t, E, G);
  // Guard extremely small G (should only happen if kt tiny; handled by k~0 already)
  G = clamp_pos(G, 1e-300);

  const double C1 = (k * b) / G;
  const double C2 = (k * E) / G;

  const double c = (v - C1) / sv;
  const double m = C2 / sv;

  if (A < eps) {
    double F = pnormP(c) / denom;
    return log_out ? safe_log(F) : F;
  }

  // Small-slope / large-time limit: m -> 0, d_req becomes independent of a
  if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
    // midpoint is a slightly better approximation than pnormP(c) alone when mA is small but not microscopic
    double F = pnormP(c + 0.5 * m * A) / denom;
    return log_out ? safe_log(F) : F;
  }

  double z0 = c;
  double zA = c + m * A;

  double F = (Hfun(zA) - Hfun(z0)) / (A * m);
  F /= denom;

  if (!(F >= 0.0) || std::isnan(F)) F = 0.0;
  return log_out ? safe_log(F) : F;
}

// [[Rcpp::export]]
double dleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false) {

  const double eps  = 1e-10;
  const double denom_floor = 1e-300;

  // Truncation normalizer for D>0
  double denom = 1.0;
  if (posdrift) {
    denom = pnormP(v / sv);
    denom = clamp_pos(denom, denom_floor);
  }

  // k -> 0 branch (no leak):
  // d_req(a,t) = (b-a)/t
  // |d d_req/dt| = (b-a)/t^2
  // should reduce to the standard LBA
  if (std::fabs(k) < eps) {
    const double c = (v - (b / t)) / sv;
    const double m = (1.0 / t) / sv;
    const double B = 1.0 / (t * t); // replaces k^2 E / (1-E)^2

    if (A < eps) {
      // a=0 -> (b-a)=b, z=c, f_D = (1/sv)dnormP(c)
      double f = B * b * (dnormP(c) / sv);
      f /= denom;
      if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
      return log_out ? safe_log(f) : f;
    }

    // small-slope approx
    if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
      // z approx at midpoint; average(b-a)=b-A/2
      double zmid = c + 0.5 * m * A;
      double f = (B * (b - 0.5 * A) * (dnormP(zmid) / sv));
      f /= denom;
      if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
      return log_out ? safe_log(f) : f;
    }

    // exact closed form:
    // f(t) = B/(A*sv*m^2) * [ (b m + c)(pnormP(zA)-pnormP(z0)) + (dnormP(zA)-dnormP(z0)) ]
    double z0 = c;
    double zA = c + m * A;

    double dpnormP = pnormP(zA) - pnormP(z0);
    double ddnormP = dnormP(zA) - dnormP(z0);

    double bracket = (b * m + c) * dpnormP + ddnormP;
    double f = (B / (A * sv * m * m)) * bracket;

    f /= denom;
    if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
    return log_out ? safe_log(f) : f;
  }

  // General k > 0
  double E, G;
  leak_terms(k, t, E, G);
  G = clamp_pos(G, 1e-300);

  const double C1 = (k * b) / G;
  const double C2 = (k * E) / G;

  const double c = (v - C1) / sv;
  const double m = C2 / sv;

  // Jacobian base B(t) = k^2 E / (1-E)^2
  const double B = (k * k * E) / (G * G);

  if (A < eps) {
    // a=0 -> (b-a)=b
    double f = B * b * (dnormP(c) / sv);
    f /= denom;
    if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
    return log_out ? safe_log(f) : f;
  }

  // small-slope / large-time approx (m -> 0)
  if (std::fabs(m) < eps || std::fabs(m * A) < 1e-8) {
    double zmid = c + 0.5 * m * A;
    double f = (B * (b - 0.5 * A) * (dnormP(zmid) / sv));
    f /= denom;
    if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
    return log_out ? safe_log(f) : f;
  }

  // exact closed form:
  double z0 = c;
  double zA = c + m * A;

  double dpnormP = pnormP(zA) - pnormP(z0);
  double ddnormP = dnormP(zA) - dnormP(z0);

  double bracket = (b * m + c) * dpnormP + ddnormP;
  double f = (B / (A * sv * m * m)) * bracket;

  f /= denom;
  if (!(f >= 0.0) || std::isnan(f)) f = 0.0;
  return log_out ? safe_log(f) : f;
}


// [[Rcpp::export]]
NumericVector dleakyba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv, NumericVector k,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector pdf(n);

  for (int i = 0; i < n; i++){
    pdf[i] = dleakyba_norm(t[i], A[i], b[i], v[i], sv[i], k[i], posdrift);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pleakyba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv, NumericVector k,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector cdf(n);

  for (int i = 0; i < n; i++){
    cdf[i] = pleakyba_norm(t[i], A[i], b[i], v[i], sv[i], k[i], posdrift);
  }
  return cdf;
}

NumericVector dleakyba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool posdrift = true){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4, k = 5
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = dleakyba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), pars(i, 5), posdrift);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericVector pleakyba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool posdrift = true){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4, k = 5
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = pleakyba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), pars(i, 5), posdrift);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}


#endif



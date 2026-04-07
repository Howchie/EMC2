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
  // v = 0, sv = 1, B = 2, A = 3, t0 = 4
  const int n_rows = rts.size();
  const int n_out = sum(idx);
  NumericVector out(n_out);

  const double* rt = rts.begin();
  const double* v  = &pars(0, 0);
  const double* sv = &pars(0, 1);
  const double* B  = &pars(0, 2);
  const double* A  = &pars(0, 3);
  const double* t0 = &pars(0, 4);
  int* idx_ptr = LOGICAL(idx);
  int* ok_ptr = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < n_rows; ++i) {
    if (!idx_ptr[i]) continue;
    if (R_IsNA(v[i])) {
      out[k++] = 0.0;
      continue;
    }
    const double t_eff = rt[i] - t0[i];
    if (t_eff > 0.0 && ok_ptr[i]) {
      out[k++] = dlba_norm(t_eff, A[i], B[i] + A[i], v[i], sv[i], posdrift);
    } else {
      out[k++] = min_ll;
    }
  }
  return out;
}

NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool posdrift = true){
  // v = 0, sv = 1, B = 2, A = 3, t0 = 4
  const int n_rows = rts.size();
  const int n_out = sum(idx);
  NumericVector out(n_out);

  const double* rt = rts.begin();
  const double* v  = &pars(0, 0);
  const double* sv = &pars(0, 1);
  const double* B  = &pars(0, 2);
  const double* A  = &pars(0, 3);
  const double* t0 = &pars(0, 4);
  int* idx_ptr = LOGICAL(idx);
  int* ok_ptr = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < n_rows; ++i) {
    if (!idx_ptr[i]) continue;
    if (R_IsNA(v[i])) {
      out[k++] = 0.0;
      continue;
    }
    const double t_eff = rt[i] - t0[i];
    if (t_eff > 0.0 && ok_ptr[i]) {
      out[k++] = plba_norm(t_eff, A[i], B[i] + A[i], v[i], sv[i], posdrift);
    } else {
      out[k++] = min_ll;
    }
  }
  return out;
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

#endif


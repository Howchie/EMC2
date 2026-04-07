#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <RcppArmadillo.h>
#include "wald_functions.h"

using namespace Rcpp;

NumericVector drdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  const int n_rows = rts.size();
  NumericVector out(sum(idx));
  const double* rt = rts.begin();
  const double* v  = &pars(0, 0);
  const double* B  = &pars(0, 1);
  const double* A  = &pars(0, 2);
  const double* t0 = &pars(0, 3);
  const double* s  = &pars(0, 4);
  int* idx_ptr = LOGICAL(idx);
  int* ok_ptr = LOGICAL(is_ok);
  int k = 0;
  for(int i = 0; i < n_rows; i++){
    if(idx_ptr[i]){
      if(R_IsNA(v[i])){
        out[k] = 0.0;
      } else if((rt[i] - t0[i] > 0.0) && ok_ptr[i]){
        const double inv_s = 1.0 / s[i];
        out[k] = digt(rt[i] - t0[i], B[i] * inv_s + 0.5 * A[i] * inv_s, v[i] * inv_s, 0.5 * A[i] * inv_s);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector prdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  const int n_rows = rts.size();
  NumericVector out(sum(idx));
  const double* rt = rts.begin();
  const double* v  = &pars(0, 0);
  const double* B  = &pars(0, 1);
  const double* A  = &pars(0, 2);
  const double* t0 = &pars(0, 3);
  const double* s  = &pars(0, 4);
  int* idx_ptr = LOGICAL(idx);
  int* ok_ptr = LOGICAL(is_ok);
  int k = 0;
  for(int i = 0; i < n_rows; i++){
    if(idx_ptr[i]){
      if(R_IsNA(v[i])){
        out[k] = 0.0;
      } else if((rt[i] - t0[i] > 0.0) && ok_ptr[i]){
        const double inv_s = 1.0 / s[i];
        out[k] = pigt(rt[i] - t0[i], B[i] * inv_s + 0.5 * A[i] * inv_s, v[i] * inv_s, 0.5 * A[i] * inv_s);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}


// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return cdf;
}

#endif

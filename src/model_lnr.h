#ifndef lnr_h
#define lnr_h

#include <RcppArmadillo.h>
#include "utility_functions.h"

using namespace Rcpp;

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  // 0 = m, 1 = s, 2 = t0
  const int n_rows = rts.size();
  NumericVector out(sum(idx));
  const double* rt = rts.begin();
  const double* m  = &pars(0, 0);
  const double* s  = &pars(0, 1);
  const double* t0 = &pars(0, 2);
  int* idx_ptr = LOGICAL(idx);
  int* ok_ptr = LOGICAL(is_ok);
  int k = 0;
  for(int i = 0; i < n_rows; i++){
    if(idx_ptr[i]){
      if(R_IsNA(m[i])){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rt[i] - t0[i] > 0.0) && ok_ptr[i]){
        out[k] = R::plnorm(rt[i] - t0[i], m[i], s[i], TRUE, FALSE);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  const int n_rows = rts.size();
  NumericVector out(sum(idx));
  const double* rt = rts.begin();
  const double* m  = &pars(0, 0);
  const double* s  = &pars(0, 1);
  const double* t0 = &pars(0, 2);
  int* idx_ptr = LOGICAL(idx);
  int* ok_ptr = LOGICAL(is_ok);
  int k = 0;
  for(int i = 0; i < n_rows; i++){
    if(idx_ptr[i]){
      if(R_IsNA(m[i])){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rt[i] - t0[i] > 0.0) && ok_ptr[i]){
        out[k] = R::dlnorm(rt[i] - t0[i], m[i], s[i], FALSE);
      } else{
        out[k] = min_ll;
      }
      k++;
    }

  }

  return(out);
}


#endif

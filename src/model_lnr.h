#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "utility_functions.h"

using namespace Rcpp;

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i]){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = log_out ? R::plnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), true, true) : R::plnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), true, false);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }

  return(out);
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool log_out=false){
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i]){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = log_out ? R::dlnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), true) : R::dlnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), false);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }

  }

  return(out);
}


#endif

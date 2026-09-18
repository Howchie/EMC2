// [[Rcpp::depends(EMC2)]]
#include <Rcpp.h>
#include "EMC2/userfun.hpp"

Rcpp::NumericVector serialized_custom_trend(
    Rcpp::NumericMatrix trend_pars,
    Rcpp::NumericMatrix input) {
  const int n = input.nrow();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = trend_pars(i, 0) * input(i, 0);
  }
  return out;
}

// [[Rcpp::export]]
SEXP EMC2_make_serialized_custom_trend_ptr();
EMC2_MAKE_PTR(serialized_custom_trend);

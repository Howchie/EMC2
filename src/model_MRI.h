#ifndef mri_H
#define mri_H

#include <Rcpp.h>
#include <cmath>

// Definitions live in src/model_MRI.cpp. Default min_ll arguments are declared
// here only; the definitions in the cpp inherit them.

double c_log_likelihood_MRI(Rcpp::NumericMatrix pars, Rcpp::NumericVector y, Rcpp::LogicalVector is_ok,
                            int n, int m,
                            double min_ll = std::log(1e-10));

double c_log_likelihood_MRI_white(Rcpp::NumericMatrix pars, Rcpp::NumericVector y, Rcpp::LogicalVector is_ok,
                                  int n, int m,
                                  double min_ll = std::log(1e-10));

// Returns the first column of `data` whose name is not one of the excluded
// bookkeeping columns ("subjects", "run", "time", "trials"); an empty vector
// if none qualifies.
Rcpp::NumericVector extract_y(Rcpp::DataFrame data);

#endif

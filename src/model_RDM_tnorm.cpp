#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// forward declarations of density functions implemented in model_RDM.h
double digt_tnorm(double t, double k = 1., double mu = 1., double sd2 = 1.,
                  double a = .1, double threshold = 1e-10);
double pigt_tnorm(double t, double k = 1., double mu = 1., double sd2 = 1.,
                  double a = .1, int steps = 100);

// Implementation of the truncated normal drift versions

NumericVector dWald_tnorm(NumericVector t, NumericVector mu,
                          NumericVector sd2, NumericVector B,
                          NumericVector A, NumericVector t0);
NumericVector pWald_tnorm(NumericVector t, NumericVector mu,
                          NumericVector sd2, NumericVector B,
                          NumericVector A, NumericVector t0);

// include function definitions from header (not inline)

// [[Rcpp::export]]
NumericVector dWald_tnorm(NumericVector t, NumericVector mu,
                          NumericVector sd2, NumericVector B,
                          NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt_tnorm(t[i], B[i] + .5 * A[i], mu[i], sd2[i], .5 * A[i]);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pWald_tnorm(NumericVector t, NumericVector mu,
                          NumericVector sd2, NumericVector B,
                          NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt_tnorm(t[i], B[i] + .5 * A[i], mu[i], sd2[i], .5 * A[i]);
    }
  }
  return cdf;
}

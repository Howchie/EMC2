#ifndef lba_h
#define lba_h

#include <Rcpp.h>
#include "utility_functions.h"

using namespace Rcpp;

double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log_out = false){
  return R::pnorm(q, mean, sd, lower, log_out);
}

double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log_out = false){
  return R::dnorm(x, mean, sd, log_out);
}

// [[Rcpp::export]]
double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true, bool log_out=false){
  double denom = 1.0;
  if (posdrift) {
    denom = pnormP(v / sv, 0.0, 1.0, true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double cdf_val;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;
    double xx = cmz - A;
    double cz = cmz / zs;
    double cz_max = xx / zs;
    cdf_val = (1. + (zs * (dnormP(cz_max, 0.0, 1.0, false) - dnormP(cz, 0.0, 1.0, false))
                   + xx * pnormP(cz_max, 0.0, 1.0, true, false) - cmz * pnormP(cz, 0.0, 1.0, true, false))/A) / denom;
  } else {
    cdf_val = pnormP(b / t, v, sv, false, false) / denom;
  }

	if (cdf_val < 0.0 || std::isnan(cdf_val)) {
		return log_out ? R_NegInf : 0;
	}
	if (cdf_val > 1.0) {
		return log_out ? 0 : 1;
	}
	return log_out ? std::log(cdf_val) : cdf_val;
}

// [[Rcpp::export]]
double dlba_norm(double t, double A,double b, double v, double sv,
                 bool posdrift = true, bool log_out=false){
  double denom = 1.0;
  if (posdrift) {
    denom = pnormP(v / sv, 0.0, 1.0, true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double pdf_val;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;;
    double cz = cmz / zs;
    double cz_max = (cmz - A) / zs;
    pdf_val = (v * (pnormP(cz, 0.0, 1.0, true, false) - pnormP(cz_max, 0.0, 1.0, true, false)) +
      sv * (dnormP(cz_max, 0.0, 1.0, false) - dnormP(cz, 0.0, 1.0, false))) / (A * denom);
  } else {
    pdf_val = dnormP(b / t, v, sv, false) * b / (t * t * denom);
  }

  if (pdf_val < 0.0) {
    return log_out ? R_NegInf : 0;
  }
  return log_out ? std::log(pdf_val) : pdf_val;
}

// [[Rcpp::export]]
NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool use_posdrift = true, bool log_out=false){ // Added use_posdrift
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i]){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i])){
        // Pass use_posdrift to dlba_norm
        out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), use_posdrift,log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
      }
      k++;
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool use_posdrift = true, bool log_out=false){ // Added use_posdrift
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i]){
      if(NumericVector::is_na(pars(i,0))){ // for RACE
        out[k] = log_out ? R_NegInf : 0; 
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i])){
        // Pass use_posdrift to plba_norm
        out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), use_posdrift,log_out);
      } else{
        out[k] = log_out ? min_ll : std::exp(min_ll);
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
#endif
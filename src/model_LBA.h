#ifndef lba_h
#define lba_h

#include <Rcpp.h>
#include "utility_functions.h" // Assuming this provides ifelse or similar if needed

using namespace Rcpp;

// Scalar versions (can be kept for reference or specific uses if any)
// Or remove if dlba_norm_vec and plba_norm_vec fully supersede their usage
double pnormP_scalar(double q, double mean = 0.0, double sd = 1.0,
                     bool lower = true, bool log = false){
  return R::pnorm(q, mean, sd, lower, log);
}

double dnormP_scalar(double x, double mean = 0.0, double sd = 1.0,
                     bool log = false){
  return R::dnorm(x, mean, sd, log);
}


NumericVector plba_norm(NumericVector t, NumericVector A, NumericVector b, // b is B+A from original
                            NumericVector v, NumericVector sv, bool posdrift = true) {
  int n = t.size();
  NumericVector denom(n, 1.0);

  if (posdrift) {
    denom = pnorm_vec(v / sv);
	for (int i = 0; i < n; ++i) {
      if (denom[i] < 1e-10) denom[i] = 1e-10;
    }
  }

  NumericVector pdf(n);
  // Check A for branching math
  LogicalVector tiny_A = (A <= 1e-10);
  NumericVector cdf(n);


  // Case 1: A > 1e-10
  // Calculations for elements where cond_A_gt_eps is true
  NumericVector zs = t * sv;

  NumericVector cmz = b - t * v; // b is already B+A effectively, so b_val in original code is b here.
                                     // The original 'b' parameter for plba_norm is B+A.
                                     // The 'xx' term used 'cmz - A'. Here, if b = B+A, then b = b - A.
                                     // So, cmz = (B+A) - t*v. xx = (B+A) - t*v - A = B - t*v.
                                     // Let's use 'b' as the parameter name for (B+A) as in original plba_norm call structure.
                                     // And 'A_param' for A.
                                     // So, cmz = b - t*v. xx = b - A_param - t*v.
  NumericVector xx = cmz - A; // A is A_param
  NumericVector cz = cmz / zs;
  NumericVector cz_max = xx / zs;

  NumericVector dnorm_cz_max = dnorm_vec(cz_max);
  NumericVector dnorm_cz = dnorm_vec(cz);
  NumericVector pnorm_cz_max = pnorm_vec(cz_max);
  NumericVector pnorm_cz = pnorm_vec(cz);

  NumericVector cdf_A = (1.0 + (zs * (dnorm_cz_max - dnorm_cz) +
                           xx * pnorm_cz_max - cmz * pnorm_cz) / A) / denom;

  // Case 2: A <= 1e-10
  // Ensure t is not zero
  NumericVector b_div_t = b / t; // Here b is B+A. If A is small, b approx B.
  NumericVector cdf_noA = pnorm_vec(b_div_t, v, sv, false, false) / denom; // pnorm_vec is vectorized, handles v, sv as vectors

  // Combine using ifelse logic
  for(int i=0; i<n; ++i) {
    if(tiny_A[i]) {
      cdf[i] = cdf_noA[i];
    } else {
      cdf[i] = cdf_A[i];
    }
    // Clamp values
    if (cdf[i] < 0.0) cdf[i] = 0.0;
    else if (cdf[i] > 1.0) cdf[i] = 1.0;
  }
  return cdf;
}


NumericVector dlba_norm(NumericVector t, NumericVector A, NumericVector b, // b is B+A from original
                            NumericVector v, NumericVector sv, bool posdrift = true) {
  int n = t.size();
  NumericVector denom(n, 1.0);

  if (posdrift) {
    denom = pnorm_vec(v / sv);
	for (int i = 0; i < n; ++i) {
      if (denom[i] < 1e-10) denom[i] = 1e-10;
    }
  }

  NumericVector pdf(n);
  LogicalVector tiny_A = (A <= 1e-10);

  // Case 1: A > 1e-10
  NumericVector zs = t * sv;
  NumericVector cmz = b - t * v; // b is B+A
  NumericVector cz = cmz / zs;
  NumericVector cz_max = (cmz - A) / zs;

  NumericVector pnorm_cz = pnorm_vec(cz);
  NumericVector pnorm_cz_max = pnorm_vec(cz_max);
  NumericVector dnorm_cz = dnorm_vec(cz);
  NumericVector dnorm_cz_max = dnorm_vec(cz_max);

  NumericVector pdf_A = (v * (pnorm_cz - pnorm_cz_max) +
                           sv * (dnorm_cz_max - dnorm_cz)) / (A * denom);

  // Case 2: A <= 1e-10
  NumericVector t_sq = t * t;
  NumericVector b_div_t = b / t; // b is B+A. If A is small, b approx B.
  NumericVector pdf_noA = dnorm_vec(b_div_t, v, sv, false) * b / (t_sq * denom);

  // Combine
  for(int i=0; i<n; ++i) {
    if(tiny_A[i]) {
      pdf[i] = pdf_noA[i];
    } else {
      pdf[i] = pdf_A[i];
    }
    if (pdf[i] < 0.0) pdf[i] = 0.0;
  }
  return pdf;
}


// Keep scalar versions if they are used elsewhere or for clarity/testing
/* double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP_scalar(v / sv, 0., 1., true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double cdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v; // b here is B+A
    double xx = cmz - A;    // xx is (B+A)-A = B
    double cz = cmz / zs;
    double cz_max = xx / zs;
    cdf = (1. + (zs * (dnormP_scalar(cz_max, 0., 1., false) - dnormP_scalar(cz, 0., 1., false))
                   + xx * pnormP_scalar(cz_max, 0., 1., true, false) - cmz * pnormP_scalar(cz, 0., 1., true, false))/A) / denom;
  } else {
    cdf = pnormP_scalar(b / t, v, sv, false, false) / denom; // b here is B+A
  }

  if (cdf < 0.) {
    return 0.;
  } else if (cdf > 1.){
    return 1.;
  }
  return cdf;
} */

/* double dlba_norm(double t, double A,double b, double v, double sv, // b here is B+A
                 bool posdrift = true){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP_scalar(v / sv, 0., 1., true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double pdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v; // b is B+A
    double cz = cmz / zs;
    double cz_max = (cmz - A) / zs; // (B+A-A)/zs = B/zs
    pdf = (v * (pnormP_scalar(cz, 0., 1., true, false) - pnormP_scalar(cz_max, 0., 1., true, false)) +
      sv * (dnormP_scalar(cz_max, 0., 1., false) - dnormP_scalar(cz, 0., 1., false))) / (A * denom);
  } else {
    pdf = dnormP_scalar(b / t, v, sv, false) * b / (t * t * denom); // b is B+A
  }

  if (pdf < 0.) {
    return 0.;
  }
  return pdf;
} */


NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool use_posdrift = true){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n_rts = rts.length();
  NumericVector out(n_rts);
  std::fill(out.begin(), out.end(), NA_REAL); // pre-fill NA_REAL for probs

  NumericVector arg_t = rts-pars(_, 4);
  NumericVector arg_A = pars(_, 3);
  NumericVector arg_b = pars(_, 2)+pars(_, 3); // This will be B + A
  NumericVector arg_v = pars(_, 0);
  NumericVector arg_sv = pars(_, 1);
  Rcout<<arg_t<<arg_v;
  NumericVector computed_values = dlba_norm(arg_t, arg_A, arg_b, arg_v, arg_sv, use_posdrift);
  for(int i = 0; i < n_rts; i++){
    if(idx[i]){
      if(NumericVector::is_na(pars(i,0))){
        out[i] = 0;
      } else if((arg_t[i]> 0) && (is_ok[i])){
        out[i] = computed_values[i];
      } else{
        out[i] = min_ll;
      }
    } 
  }
  return out;
}

NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok, bool use_posdrift = true){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n_rts = rts.length();
  NumericVector out(n_rts);
  std::fill(out.begin(), out.end(), NA_REAL); // pre-fill NA_REAL for probs

  NumericVector arg_t = rts-pars(_, 4);
  NumericVector arg_A = pars(_, 3);
  NumericVector arg_b = pars(_, 2)+pars(_, 3); // This will be B + A
  NumericVector arg_v = pars(_, 0);
  NumericVector arg_sv = pars(_, 1);
  Rcout<<arg_t<<arg_v;
  NumericVector computed_values = plba_norm(arg_t, arg_A, arg_b, arg_v, arg_sv, use_posdrift);
  for(int i = 0; i < n_rts; i++){
    if(idx[i]){
      if(NumericVector::is_na(pars(i,0))){
        out[i] = 0;
      } else if((arg_t[i]> 0) && (is_ok[i])){
        out[i] = computed_values[i];
      } else{
        out[i] = min_ll;
      }
    } 
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv, // Here 'b' is B+A for consistency with dlba_norm
                   bool posdrift = true)
{

  if (t.length() == 0) return NumericVector(0);
  return dlba_norm(t, A, b, v, sv, posdrift);
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv, // Here 'b' is B+A
                   bool posdrift = true)
{

  if (t.length() == 0) return NumericVector(0);
  return plba_norm(t, A, b, v, sv, posdrift);
}

#endif


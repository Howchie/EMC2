#include "model_LBA.h"
#include "model_BAwD.h"
#include "model_BAwF.h"
#include "model_BAwR.h"
#include "model_BTAwL.h"
// CDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double pleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false,
                     int launch = 0, double delta = 0.0) {
  // At infinite time k = 0 has the usual LBA limit; for k > 0 only drifts
  // above k*b can finish, and the m = 0 point limit inside the evaluators
  // retains that defective upper tail instead of returning one.
  return bawl_cdf_norm(t, A, b, v, sv, k, posdrift, log_out,
                       BAWL_DENOM_FLOOR, launch, delta);
}
// PDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double dleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false,
                     int launch = 0, double delta = 0.0) {
  return bawl_pdf_norm(t, A, b, v, sv, k, posdrift, log_out,
                       BAWL_DENOM_FLOOR, launch, delta);
}
// [[Rcpp::export]]
NumericVector dkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0,
                             int launch = 0, NumericVector delta = 0.0) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    pdf[i] = dkilledleakyba_norm(
      t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
      pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
      posdrift, log_out, kill_shape, guess, omega, launch, pick(delta, i));
  }
  return pdf;
}
// [[Rcpp::export]]
NumericVector pkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0,
                             int launch = 0, NumericVector delta = 0.0) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    cdf[i] = pkilledleakyba_norm(
      t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
      pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
      posdrift, log_out, kill_shape, guess, omega, launch, pick(delta, i));
  }
  return cdf;
}
// Vectorised R-callable wrappers (recycle scalar parameters).
// [[Rcpp::export]]
NumericVector dleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true, int launch = 0,
                       NumericVector delta = 0.0) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = dleakyba_norm(
      t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i),
      posdrift, false, launch, pick(delta, i));
  return pdf;
}
// [[Rcpp::export]]
NumericVector pleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true, int launch = 0,
                       NumericVector delta = 0.0) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = pleakyba_norm(
      t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i),
      posdrift, false, launch, pick(delta, i));
  return cdf;
}
// Standard LBA (exact k = 0 member with the legacy LBA normalizer floor),
// restored so the R-side dfun/pfun agree exactly with the C++ likelihood
// kernels, which also use LBA_DENOM_FLOOR for this model.
// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b,
                   NumericVector v, NumericVector sv,
                   bool posdrift = true, bool log_out = false) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = lba_k0_pdf_norm(t[i], pick(A,i), pick(b,i), pick(v,i),
                             pick(sv,i), posdrift, log_out);
  return pdf;
}
// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b,
                   NumericVector v, NumericVector sv,
                   bool posdrift = true, bool log_out = false) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = lba_k0_cdf_norm(t[i], pick(A,i), pick(b,i), pick(v,i),
                             pick(sv,i), posdrift, log_out);
  return cdf;
}

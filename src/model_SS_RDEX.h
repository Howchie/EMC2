#ifndef ss_rdex_h
#define ss_rdex_h

#include <cmath>
#include <vector>
#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "wald_functions.h"
#include "exgaussian_functions.h"
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------
// HYBRID WALD / EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around Wald log PDF for race function
NumericVector rdex_go_lpdf(
    // single RT, broadcast to repped vector (one per accumulator)
    NumericVector rt,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // accumulator index: for which accumulator(s) should log density be computed?
    LogicalVector idx,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  
  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;
  
  NumericVector out(n_acc_selected);
  int k = 0;
  
  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;
    
    double dt_i = rt[i] - pars(i, 3);
    double log_d = R_NegInf;
    if (dt_i > 0.) {
      log_d = std::log(
        digt(
          dt_i,
          (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
          pars(i, 0) / pars(i, 4),
          .5 * (pars(i, 2) / pars(i, 4))
        )
      );
    }
    
    out[k] = std::isfinite(log_d) ? log_d : R_NegInf;
    
    k++;
  }
  
  return(out);
}

// wrapper around Wald log complementary CDF for race function
NumericVector rdex_go_lccdf(
    // single RT, broadcast to repped vector (one per accumulator)
    NumericVector rt,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // accumulator index: for which accumulator(s) should survivor probability be computed?
    LogicalVector idx,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  
  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;
  
  NumericVector out(n_acc_selected);
  int k = 0;
  
  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;
    
    double dt_i = rt[i] - pars(i, 3);
    double log_s = 0.; // log(1)
    if (dt_i > 0.) {
      log_s = log1m(
        pigt(
          dt_i,
          (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
          pars(i, 0) / pars(i, 4),
          .5 * (pars(i, 2) / pars(i, 4))
        )
      );
    }
    
    out[k] = std::isfinite(log_s) ? log_s : R_NegInf;
    
    k++;
  }
  
  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_rdex_go_lpdf(
    // single RT
    double RT,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  const int n_acc = pars.nrow();
  double logpdf_winner = 0.0;
  double logsurv_losers = 0.0;
  bool winner_found = false;
  
  for (int i = 0; i < n_acc; ++i) {
    double s = pars(i, 4);
    double alpha = (pars(i, 1) / s) + .5 * (pars(i, 2) / s);
    double nu    =  pars(i, 0) / s;
    double gamma = .5 * (pars(i, 2) / s);
    double t0    =  pars(i, 3);
    double dt_i = RT - t0;

    if (winner[i]) {
      double ld = R_NegInf;
      if (dt_i > 0.) {
        ld = std::log(digt(dt_i, alpha, nu, gamma));
      }
      logpdf_winner += std::isfinite(ld) ? ld : min_ll;
      winner_found = true;
    } else {
      double ls = 0.; // log(1)
      if (dt_i > 0.) {
        ls = log1m(pigt(dt_i, alpha, nu, gamma));
      }
      logsurv_losers += std::isfinite(ls) ? ls : min_ll;
    }
  }
  if (!winner_found) return min_ll;
  double out = logpdf_winner + logsurv_losers;
  return std::isfinite(out) ? out : min_ll;
}

// go vs stop race log likelihood function for the case of stop trials with
// failed inhibition (i.e., stop process losing), not accounting for trigger
// failure and go failure.
double ss_rdex_stop_fail_lpdf(
    // single RT
    double RT,
    // single stop signal delay
    double SSD,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  // obtain the go race process log likelihood
  double go_lprob = ss_rdex_go_lpdf(RT, pars, winner, min_ll);
  // obtain the survivor log probability of the stop process
  // NB SSD subtracted from observed RT to get stop finish time
  // stop pars: muS=5, sigS=6, tauS=7, lbS=10
  double stop_survivor_lprob = ptexg(
    RT - SSD, pars(0, 5), pars(0, 6), pars(0, 7), pars(0, 10), R_PosInf, false, true
  );
  if (!std::isfinite(stop_survivor_lprob)) {
    stop_survivor_lprob = min_ll;
  }
  // final output of race model is summed log likelihood
  return go_lprob + stop_survivor_lprob;
}

// function to compute the stop success integral, not accounting for trigger
// failure and go failure
static inline double ss_rdex_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll,
    double upper = R_PosInf,
    int max_subdiv = 30,
    double abs_tol = 1e-5,
    double rel_tol = 1e-4,
    double k_sigma = 6.0,
    double k_tau = 12.0
) {
  gsl_function F;
  struct Wrapper {
    NumericMatrix pars;
    double SSD;
  } w_params = {pars, SSD};
  
  auto integrand = [](double x, void* p) -> double {
    Wrapper* wp = static_cast<Wrapper*>(p);
    const int n_acc = wp->pars.nrow();
    // Stop process (muS=5, sigS=6, tauS=7, lbS=10)
    double muS = wp->pars(0, 5);
    double sigS = wp->pars(0, 6);
    double tauS = wp->pars(0, 7);
    double lbS  = wp->pars(0, 10);
    double fS = dtexg(x, muS, sigS, tauS, lbS, R_PosInf, false);
    if (fS <= 0.0) return 0.0;
    
    double S_go_all = 1.0;
    for (int i = 0; i < n_acc; ++i) {
      double s = wp->pars(i, 4);
      double alpha = (wp->pars(i, 1) / s) + .5 * (wp->pars(i, 2) / s);
      double nu    =  wp->pars(i, 0) / s;
      double gamma = .5 * (wp->pars(i, 2) / s);
      double t0    =  wp->pars(i, 3);
      double dt_i = (x + wp->SSD) - t0;
      double Si = 1.0;
      if (dt_i > 0.) {
        Si = 1.0 - pigt(dt_i, alpha, nu, gamma);
      }
      S_go_all *= Si;
      if (S_go_all <= 0.0) return 0.0;
    }
    return fS * S_go_all;
  };

  F.function = integrand;
  F.params = &w_params;

  double muS = pars(0, 5);
  double sigS = pars(0, 6);
  double tauS = pars(0, 7);
  double lbS  = pars(0, 10);
  double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  // Use emc2_isfinite / emc2_isinf (not std:: versions) — -ffast-math breaks them
  double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lbS)) ub = lbS + 1e-12;

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, max_subdiv);
  double res, err;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

  int status;
  if (emc2_isinf(ub)) {
    status = gsl_integration_qagiu(&F, lbS, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  } else {
    status = gsl_integration_qags(&F, lbS, ub, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  }

  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// ----------------------------------------------------------------------------
// OLD STUFF BELOW, KEPT FOR TESTING R CODE
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector pEXG_RDEX(NumericVector q,
                        double mu = 5., double sigma = 1., double tau = 1.,
                        bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  if (tau <= 0 || sigma <= 0) {
    NumericVector cdf(n, NA_REAL);
    return cdf;
  }
  
  NumericVector cdf(n);
  if (sigma < 1e-4){
    for (int i = 0; i < n; i++){
      cdf[i] = R::pexp(q[i] - mu, tau, lower_tail, log_p);
    }
    return cdf;
  }
  
  for (int i = 0; i < n; i++){
    if (!traits::is_infinite<REALSXP>(q[i])){
      if (tau > .05 * sigma){
        double z_i = q[i] - mu - (sigma * sigma) / tau;
        double mu_term = mu + (sigma * sigma / tau);
        cdf[i] = R::pnorm((q[i] - mu) / sigma, 0., 1., true, false) - std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) + (mu_term * mu_term - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
      } else {
        cdf[i] = R::pnorm(q[i], mu, sigma, true, false);
      }
    } else {
      if (q[i] < 0) {
        cdf[i] = 0.;
      } else {
        cdf[i] = 1.;
      }
    }
  }
  if (!lower_tail){
    for(int i = 0; i < n; i++){
      cdf[i] = 1. - cdf[i];
    }
  }
  if (log_p){
    for(int i = 0; i < n; i++){
      cdf[i] = std::log(cdf[i]);
    }
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector dEXG_RDEX(NumericVector x,
                        double mu = 5., double sigma = 1., double tau = 1.,
                        bool log_d = false) {
  int n = x.size();
  if (tau <= 0 || sigma <= 0) {
    NumericVector pdf(n, NA_REAL);
    return pdf;
  }
  
  NumericVector pdf(n);
  if (sigma < 1e-4){
    for (int i = 0; i < n; i++){
      pdf[i] = R::dexp(x[i] - mu, tau, log_d);
    }
    return pdf;
  }
  
  for (int i = 0; i < n; i++){
    if (tau > .05 * sigma){
      double z_i = x[i] - mu - (sigma * sigma) / tau;
      pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + std::log(R::pnorm(z_i / sigma, 0., 1., true, false));
    } else {
      pdf[i] = R::dnorm(x[i], mu, sigma, true);
    }
  }
  if (!log_d){
    for(int i = 0; i < n; i++){
      pdf[i] = std::exp(pdf[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector dWald_RDEX_old(NumericVector t, double v,
                             double B, double A, double t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0;
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B + .5 * A, v, .5 * A);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    pdf[i] = 0.;
    if (t[i] > 0.) {
      pdf[i] = digt(t[i], (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald_RDEX_old(NumericVector t, double v,
                             double B, double A, double t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0;
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B + .5 * A, v, .5 * A);
    }
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector pWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    cdf[i] = 0.;
    if (t[i] > 0.) {
      cdf[i] = pigt(t[i], (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
    }
  }
  return cdf;
}


// [[Rcpp::export]]
NumericVector pTEXG_RDEX(
    NumericVector q, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool lower_tail = true, bool log_p = false
) {
  int n = q.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector cdf(n, NA_REAL);
    return cdf;
  }
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    cdf[i] = ptexg(q[i], mu, sigma, tau, lb, R_PosInf, lower_tail, log_p);
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector dTEXG_RDEX(
    NumericVector x, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool log_d = false
) {
  int n = x.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector pdf(n, NA_REAL);
    return pdf;
  }
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    pdf[i] = dtexg(x[i], mu, sigma, tau, lb, R_PosInf, log_d);
  }
  return pdf;
}



// [[Rcpp::export]]
NumericVector dRDEXrace_old(NumericMatrix dt,
                            double mu, double sigma, double tau,
                            NumericVector v, NumericVector B, NumericVector A,
                            NumericVector t0, bool exgWinner = true){
  int n = v.size();
  NumericVector out(dt.ncol());
  if (exgWinner){
    out = dEXG_RDEX(dt(0, _), mu, sigma, tau, false);
    out = out * (1. - pWald_RDEX_old(dt(1, _), v[0], B[0], A[0], t0[0]));
  } else {
    out = dWald_RDEX_old(dt(0, _), v[0], B[0], A[0], t0[0]);
    out = out * (1. - pEXG_RDEX(dt(1, _), mu, sigma, tau));
  }
  for (int i = 1; i < n; i++){
    out = out * (1. - pWald_RDEX_old(dt(i + 1, _), v[i], B[i], A[i], t0[i]));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dRDEXrace(
    NumericMatrix dt,
    double mu, double sigma, double tau, double lb,
    NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector s,
    bool exgWinner = true
) {
  int n = v.size();
  NumericVector out(dt.ncol());
  if (exgWinner) {
    out = dTEXG_RDEX(dt(0, _), mu, sigma, tau, lb);
    out = out * (1. - pWald_RDEX(dt(1, _), v[0], B[0], A[0], t0[0], s[0]));
  } else {
    out = dWald_RDEX(dt(0, _), v[0], B[0], A[0], t0[0], s[0]);
    out = out * (1. - pTEXG_RDEX(dt(1, _), mu, sigma, tau, lb));
  }
  for (int i = 1; i < n; i++){
    out = out * (1. - pWald_RDEX(dt(i + 1, _), v[i], B[i], A[i], t0[i], s[i]));
  }
  return out;
}


// [[Rcpp::export]]
NumericVector stopfn_rdex_old(NumericVector t, int n_acc,
                              double mu, double sigma, double tau,
                              NumericVector v, NumericVector B, NumericVector A,
                              NumericVector t0, double SSD){
  NumericVector tmp( (n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace_old(dt, mu, sigma, tau, v, B, A, t0);
}

// [[Rcpp::export]]
NumericVector stopfn_rdex(
    NumericVector t, int n_acc,
    double mu, double sigma, double tau, double lb,
    NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector s,
    double SSD
) {
  NumericVector tmp((n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace(dt, mu, sigma, tau, lb, v, B, A, t0, s);
}

#endif
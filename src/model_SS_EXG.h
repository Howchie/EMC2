#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <vector>
#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "wald_functions.h"
#include "exgaussian_functions.h"
#include "gsl_utils.h"
#include "gl_quad.h"
#include "ss_exg_analytic.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;

// ----------------------------------------------------------------------------
// TRUNCATED EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around truncated ex-Gaussian log PDF for race function
NumericVector texg_go_lpdf(
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
    
    // input args: x, mu, sigma, tau, exg_lb, upper = Inf, log_d = TRUE
    double log_d = dtexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, true
    );
    out[k] = emc2_isfinite(log_d) ? log_d : R_NegInf;
    
    k++;
  }
  
  return(out);
}

// wrapper around truncated ex-Gaussian log complementary CDF for race function
NumericVector texg_go_lccdf(
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
    
    // input args: q, mu, sigma, tau, exg_lb, upper = Inf, lower_tail = FALSE, log_p = TRUE
    double log_s = ptexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true
    );
    out[k] = emc2_isfinite(log_s) ? log_s : R_NegInf;
    
    k++;
  }
  
  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_texg_go_lpdf(
    // single RT
    double RT,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  // Manual loop for racing likelihood (much faster than 1D integrator)
  const int n_acc = pars.nrow();
  double logpdf_winner = 0.0;
  double logsurv_losers = 0.0;
  bool winner_found = false;
  
  for (int i = 0; i < n_acc; ++i) {
    if (winner[i]) {
      // muG=0, sigG=1, tauG=2, lbG=8
      double ld = dtexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, true);
      logpdf_winner += emc2_isfinite(ld) ? ld : min_ll;
      winner_found = true;
    } else {
      double ls = ptexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true);
      logsurv_losers += emc2_isfinite(ls) ? ls : min_ll;
    }
  }
  if (!winner_found) return min_ll;
  double out = logpdf_winner + logsurv_losers;
  return emc2_isfinite(out) ? out : min_ll;
}

// go vs stop race log likelihood function for the case of stop trials with
// failed inhibition (i.e., stop process losing), not accounting for trigger
// failure and go failure.
double ss_texg_stop_fail_lpdf(
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
  double go_lprob = ss_texg_go_lpdf(RT, pars, winner, min_ll);
  // obtain the survivor log probability of the stop process
  // NB SSD subtracted from observed RT to get stop finish time
  // pars(0, 3..5, 9) are stop parameters
  double stop_survivor_lprob = ptexg(
    RT - SSD, pars(0, 3), pars(0, 4), pars(0, 5), pars(0, 9), R_PosInf, false, true
  );
  if (!emc2_isfinite(stop_survivor_lprob)) {
    stop_survivor_lprob = min_ll;
  }
  // final output of race model is summed log likelihood
  return go_lprob + stop_survivor_lprob;
}

// function to compute the stop success integral, not accounting for trigger
// failure and go failure
static inline double ss_texg_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll,
    double upper = R_PosInf,
    int max_subdiv = 100,
    double abs_tol = 1e-8,
    double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA,
    double k_tau = SS_WINDOW_K_TAU
) {
  gsl_function F;
  // Let's pass the NumericMatrix itself via a struct.
  struct Wrapper {
    NumericMatrix pars;
    double SSD;
  } w_params = {pars, SSD};
  
  auto integrand = [](double x, void* p) -> double {
    Wrapper* wp = static_cast<Wrapper*>(p);
    const int n_acc = wp->pars.nrow();
    // Stop process (muS=3, sigS=4, tauS=5, lbS=9)
    double muS = wp->pars(0, 3);
    double sigS = wp->pars(0, 4);
    double tauS = wp->pars(0, 5);
    double lbS  = wp->pars(0, 9);
    double fS = dtexg(x, muS, sigS, tauS, lbS, R_PosInf, false);
    if (fS <= 0.0) return 0.0;
    
    double S_go_all = 1.0;
    for (int i = 0; i < n_acc; ++i) {
      double muG = wp->pars(i, 0);
      double sigG = wp->pars(i, 1);
      double tauG = wp->pars(i, 2);
      double lbG = wp->pars(i, 8);
      double Si = ptexg(x + wp->SSD, muG, sigG, tauG, lbG, R_PosInf, false, false);
      S_go_all *= Si;
      if (S_go_all <= 0.0) return 0.0;
    }
    return fS * S_go_all;
  };

  F.function = integrand;
  F.params = &w_params;

  double muS = pars(0, 3);
  double sigS = pars(0, 4);
  double tauS = pars(0, 5);
  double lbS  = pars(0, 9);
  double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  // Use emc2_isfinite / emc2_isinf (not std:: versions) — -ffast-math breaks them
  upper = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(upper > lbS)) return min_ll;

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, max_subdiv);
  double res, err;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

  int status;
  if (emc2_isinf(upper)) {
    status = gsl_integration_qagiu(&F, lbS, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  } else {
    status = gsl_integration_qags(&F, lbS, upper, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  }

  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}


// ----------------------------------------------------------------------------
// REGULAR EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around ex-Gaussian log PDF for race function
NumericVector exg_go_lpdf(
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
    
    // input args: x, mu, sigma, tau, log_d = TRUE
    double log_d = dexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), true);
    out[k] = emc2_isfinite(log_d) ? log_d : R_NegInf;
    
    k++;
  }
  
  return(out);
}

// wrapper around ex-Gaussian log complementary CDF for race function
NumericVector exg_go_lccdf(
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
    
    // input args: q, mu, sigma, tau, lower_tail = FALSE, log_p = TRUE
    double log_s = pexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), false, true);
    out[k] = emc2_isfinite(log_s) ? log_s : R_NegInf;
    
    k++;
  }
  
  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_exg_go_lpdf(
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
    if (winner[i]) {
      // muG=0, sigG=1, tauG=2
      double ld = dexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), true);
      logpdf_winner += emc2_isfinite(ld) ? ld : min_ll;
      winner_found = true;
    } else {
      double ls = pexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), false, true);
      logsurv_losers += emc2_isfinite(ls) ? ls : min_ll;
    }
  }
  if (!winner_found) return min_ll;
  double out = logpdf_winner + logsurv_losers;
  return emc2_isfinite(out) ? out : min_ll;
}

// go vs stop race log likelihood function for the case of stop trials with
// failed inhibition (i.e., stop process losing), not accounting for trigger
// failure and go failure.
double ss_exg_stop_fail_lpdf(
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
  double go_lprob = ss_exg_go_lpdf(RT, pars, winner, min_ll);
  // obtain the survivor log probability of the stop process
  // NB SSD subtracted from observed RT to get stop finish time
  // pars(0, 3..5) are stop parameters
  double stop_survivor_lprob = pexg(
    RT - SSD, pars(0, 3), pars(0, 4), pars(0, 5), false, true
  );
  if (!emc2_isfinite(stop_survivor_lprob)) {
    stop_survivor_lprob = min_ll;
  }
  // final output of race model is summed log likelihood
  return go_lprob + stop_survivor_lprob;
}

// function to compute the stop success integral, not accounting for trigger
// failure and go failure
static inline double ss_exg_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll,
    double upper = R_PosInf,
    int max_subdiv = 100,
    double abs_tol = 1e-8,
    double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA,
    double k_tau = SS_WINDOW_K_TAU
) {
  gsl_function F;
  struct Wrapper {
    NumericMatrix pars;
    double SSD;
  } w_params = {pars, SSD};
  
  auto integrand = [](double x, void* p) -> double {
    Wrapper* wp = static_cast<Wrapper*>(p);
    const int n_acc = wp->pars.nrow();
    // Stop process (muS=3, sigS=4, tauS=5)
    double muS = wp->pars(0, 3);
    double sigS = wp->pars(0, 4);
    double tauS = wp->pars(0, 5);
    double fS = dexg(x, muS, sigS, tauS, false);
    if (fS <= 0.0) return 0.0;
    
    double S_go_all = 1.0;
    for (int i = 0; i < n_acc; ++i) {
      double muG = wp->pars(i, 0);
      double sigG = wp->pars(i, 1);
      double tauG = wp->pars(i, 2);
      double Si = pexg(x + wp->SSD, muG, sigG, tauG, false, false);
      S_go_all *= Si;
      if (S_go_all <= 0.0) return 0.0;
    }
    return fS * S_go_all;
  };

  F.function = integrand;
  F.params = &w_params;

  double muS = pars(0, 3);
  double sigS = pars(0, 4);
  double tauS = pars(0, 5);
  double halfw = k_sigma * sigS + k_tau * tauS;
  double a = muS - halfw;
  double b = muS + halfw;

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, max_subdiv);
  double res, err;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  
  int status = gsl_integration_qags(&F, a, b, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  
  gsl_set_error_handler(old_handler);
  
  if (status != GSL_SUCCESS || !emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// ----------------------------------------------------------------------------
// OLD STUFF BELOW, KEPT FOR TESTING R CODE
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dexg_c(
    const NumericVector x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = dexg(x[i], mu, sigma, tau, log_d);
  }
  return(out);
}

NumericVector pexg_c(
    const NumericVector q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {
  int n = q.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = pexg(q[i], mu, sigma, tau, lower_tail, log_p);
  }
  return(out);
}

NumericVector protect_finite(const NumericVector& x, double min_ll) {
  NumericVector out = clone(x);
  for (int i = 0; i < out.size(); i++) {
    if (!R_FINITE(out[i])) {
      out[i] = min_ll;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dEXGrace(
    NumericMatrix dt, NumericVector mu, NumericVector sigma, NumericVector tau,
    double min_ll
){
  int n = mu.size();
  NumericVector log_out(dt.ncol());
  log_out = protect_finite(
    dexg_c(dt(0, _), mu[0], sigma[0], tau[0], true),
    min_ll
  );
  for (int i = 1; i < n; i++) {
    log_out += protect_finite(
      pexg_c(dt(i, _), mu[i], sigma[i], tau[i], false, true),
      min_ll
    );
  }
  return exp(log_out);
}

// [[Rcpp::export]]
NumericVector stopfn_exg(
    NumericVector t, NumericVector mu, NumericVector sigma, NumericVector tau,
    double SSD, double min_ll
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace(dt, mu, sigma, tau, min_ll);
}

// old

// [[Rcpp::export]]
NumericVector pEXG_old(NumericVector q,
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
      double z_i = q[i] - mu - (sigma * sigma) / tau;
      double mu_term = mu + (sigma * sigma / tau);
      cdf[i] = pnorm_std((q[i] - mu) / sigma) - std::exp(pnorm_std(z_i / sigma, true, true) + (mu_term * mu_term - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
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
NumericVector dEXG_old(NumericVector x,
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
    //    if (tau > .05 * sigma){
    double z_i = x[i] - mu - (sigma * sigma) / tau;
    pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + pnorm_std(z_i / sigma, true, true);
    //    } else {
    //      pdf[i] = R::dnorm(x[i], mu, sigma, true);
    //    }
  }
  if (!log_d){
    for(int i = 0; i < n; i++){
      pdf[i] = std::exp(pdf[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector dEXGrace_old(NumericMatrix dt, NumericVector mu, NumericVector sigma, NumericVector tau){
  int n = mu.size();
  NumericVector out(dt.ncol());
  out = dEXG_old(dt(0, _), mu[0], sigma[0], tau[0], false);
  for (int i = 1; i < n; i++){
    out = out * pEXG_old(dt(i, _), mu[i], sigma[i], tau[i], false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_exg_old(NumericVector t,
                             NumericVector mu, NumericVector sigma, NumericVector tau,
                             double SSD){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace_old(dt, mu, sigma, tau);
}


// [[Rcpp::export]]
NumericVector pTEXG_vec(
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
NumericVector dTEXG_vec(
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
NumericVector dTEXGrace(
    NumericMatrix dt,
    NumericVector mu, NumericVector sigma, NumericVector tau, NumericVector lb 
){
  int n = mu.size();
  int n_points = dt.ncol();
  NumericVector out(n_points);
  out = dTEXG_vec(dt(0, _), mu[0], sigma[0], tau[0], lb[0]);
  for (int i = 1; i < n; i++){
    out = out * pTEXG_vec(dt(i, _), mu[i], sigma[i], tau[i], lb[i], false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_texg(
    NumericVector t,
    NumericVector mu, NumericVector sigma, NumericVector tau, NumericVector lb,
    double SSD
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dTEXGrace(dt, mu, sigma, tau, lb);
}

// ----------------------------------------------------------------------------
// GAUSS-LEGENDRE VARIANTS OF THE STOP-SUCCESS INTEGRALS  (added 2026-06)
// Reuse the same integrand internals as the adaptive (qags/qagiu) versions
// above; integrate over the SAME finite window with a fixed n-node GL rule.
// See gl_quad.h. These are a fast alternative to the GSL adaptive routines.
// ----------------------------------------------------------------------------

// --- truncated ex-Gaussian stop-success, GL version (LIVE PATH twin) ---
static inline double ss_texg_stop_success_lpdf_gl(
    double SSD,
    NumericMatrix pars,
    double min_ll,
    double upper = R_PosInf,
    int n_nodes = 64,
    double k_sigma = SS_WINDOW_K_SIGMA,
    double k_tau = SS_WINDOW_K_TAU
) {
  struct Wrapper { NumericMatrix pars; double SSD; } w_params = {pars, SSD};
  auto integrand = [](double x, void* p) -> double {
    Wrapper* wp = static_cast<Wrapper*>(p);
    const int n_acc = wp->pars.nrow();
    double muS = wp->pars(0,3), sigS = wp->pars(0,4), tauS = wp->pars(0,5);
    double lbS = wp->pars(0,9);
    double fS = dtexg(x, muS, sigS, tauS, lbS, R_PosInf, false);
    if (fS <= 0.0) return 0.0;
    double S_go_all = 1.0;
    for (int i = 0; i < n_acc; ++i) {
      double Si = ptexg(x + wp->SSD, wp->pars(i,0), wp->pars(i,1),
                        wp->pars(i,2), wp->pars(i,8), R_PosInf, false, false);
      S_go_all *= Si;
      if (S_go_all <= 0.0) return 0.0;
    }
    return fS * S_go_all;
  };
  const double muS = pars(0,3), sigS = pars(0,4), tauS = pars(0,5);
  const double lbS = pars(0,9);
  const double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  // lbS = -Inf (untruncated stop) needs a finite GL window: mirror the upper
  // heuristic on the lower side
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lo)) return min_ll;
  double res = gl_integrate(integrand, &w_params, lo, ub, n_nodes);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// --- "auto"/"analytic" dispatch (PLAN, Andrew 2026-06-11) ---------------------
// n_go == 1 with an infinite upper limit: exact closed form (full-line or
// truncated; see ss_exg_analytic.h). Guard trips, kinked domains, finite upper
// limits, and n_go >= 2: GL with an auto node bump for tight stop densities.
static inline double ss_texg_stop_success_lpdf_autodisp(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int n_nodes = 64,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU
) {
  if (pars.nrow() == 1 && !emc2_isfinite(upper)) {
    bool ok = false;
    double p = ss_texg_stop_success_analytic1(
      SSD, pars(0,0), pars(0,1), pars(0,2), pars(0,8),
      pars(0,3), pars(0,4), pars(0,5), pars(0,9), ok);
    if (ok) return std::log(p);
  }
  const double muS = pars(0,3), sigS = pars(0,4), tauS = pars(0,5);
  const double lbS = pars(0,9);
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  const double ub = emc2_isfinite(upper) ? upper
                                         : muS + k_sigma * sigS + k_tau * tauS;
  const int n_eff = gl_auto_nodes(n_nodes, lo, ub, sigS);
  return ss_texg_stop_success_lpdf_gl(SSD, pars, min_ll, upper, n_eff,
                                      k_sigma, k_tau);
}

// --- LIVE entry point: matches the ss_stop_success_fn signature so
//     particle_ll.cpp can point the function pointer at it. Reads the
//     process-global stop_method_config(), which calc_ll_manager() sets from
//     SSEXG(stop_method =, stop_n_nodes =) once per likelihood call.
//     "integrate" is the original qags route, numerically untouched. ---
static inline double ss_texg_stop_success_lpdf_live(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int max_subdiv = 100,
    double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU
) {
  const StopMethodConfig& c = stop_method_config();
  switch (c.method) {
  case STOP_METHOD_INTEGRATE:
    return ss_texg_stop_success_lpdf(SSD, pars, min_ll, upper, max_subdiv,
                                     abs_tol, rel_tol, k_sigma, k_tau);
  case STOP_METHOD_GL:
    return ss_texg_stop_success_lpdf_gl(SSD, pars, min_ll, upper, c.n_nodes,
                                        k_sigma, k_tau);
  default:   // STOP_METHOD_AUTO / STOP_METHOD_ANALYTIC
    return ss_texg_stop_success_lpdf_autodisp(SSD, pars, min_ll, upper,
                                              c.n_nodes, k_sigma, k_tau);
  }
}

// --- stop_method config setter/getter (called from calc_ll_manager and tests).
//     Process-global; see gl_quad.h for the threading rationale. ---
// [[Rcpp::export]]
void emc2_set_stop_method(std::string method = "auto", int n_nodes = 64) {
  StopMethodConfig& c = stop_method_config();
  if (method == "auto")           c.method = STOP_METHOD_AUTO;
  else if (method == "integrate") c.method = STOP_METHOD_INTEGRATE;
  else if (method == "gl")        c.method = STOP_METHOD_GL;
  else if (method == "analytic")  c.method = STOP_METHOD_ANALYTIC;
  else Rcpp::stop("unknown stop_method '" + method + "'");
  if (n_nodes < 2) Rcpp::stop("stop_n_nodes must be >= 2");
  c.n_nodes = n_nodes;
}
// [[Rcpp::export]]
Rcpp::List emc2_get_stop_method() {
  const StopMethodConfig& c = stop_method_config();
  static const char* method_names[4] = {"auto", "integrate", "gl", "analytic"};
  return Rcpp::List::create(
    Rcpp::_["method"]  = std::string(method_names[c.method]),
    Rcpp::_["n_nodes"] = c.n_nodes);
}

// --- which branch would "auto" take for this trial? (test/benchmark aid) ---
// [[Rcpp::export]]
std::string ss_texg_stop_success_auto_branch(
    double SSD, NumericMatrix pars, double upper = -1.0
) {
  if (upper <= 0.0) upper = R_PosInf;   // sentinel: <=0 means auto/Inf
  if (pars.nrow() != 1) return "gl_ngo2plus";
  if (emc2_isfinite(upper)) return "gl_finite_upper";
  bool ok = false;
  ss_texg_stop_success_analytic1(
    SSD, pars(0,0), pars(0,1), pars(0,2), pars(0,8),
    pars(0,3), pars(0,4), pars(0,5), pars(0,9), ok);
  if (!ok) return "gl_guard";
  const bool full_line = (pars(0,9) == R_NegInf) && (pars(0,8) == R_NegInf);
  return full_line ? "analytic_fullline" : "analytic_trunc";
}

// --- Rcpp test wrappers: return the INTEGRAL VALUE (not log) for one trial.
//     method = "integrate" -> pure qags (config-independent),
//              "gl"        -> explicit fixed Gauss-Legendre with n_nodes,
//              "analytic"/"auto" -> the auto dispatch rule with n_nodes as the
//                                   GL fallback base,
//              "live"      -> the live dispatcher (honours the config set by
//                             emc2_set_stop_method(), i.e. the real plumbing).
//     Used by testthat / benchmarking from R. Not on the fitting hot path. ---
// NB: k_sigma/k_tau defaults are literal copies of SS_WINDOW_K_SIGMA/K_TAU
// (gl_quad.h) — compileAttributes copies defaults verbatim into the generated
// R wrapper, so named constants cannot be used in this exported signature.
// max_subdiv/abs_tol/rel_tol defaults match the LIVE call sites in
// particle_ll.cpp (100, 1e-8, 1e-6) so wrapper-default results equal live.
// [[Rcpp::export]]
double ss_texg_stop_success_value(
    double SSD, NumericMatrix pars, std::string method = "integrate",
    double upper = -1.0, int n_nodes = 64,
    double k_sigma = 8.0, double k_tau = 16.0,
    int max_subdiv = 100, double abs_tol = 1e-8, double rel_tol = 1e-6
) {
  const double min_ll = -1e10;
  if (upper <= 0.0) upper = R_PosInf;   // sentinel: <=0 means auto/Inf
  double lp;
  if (method == "gl") {
    lp = ss_texg_stop_success_lpdf_gl(SSD, pars, min_ll, upper, n_nodes,
                                      k_sigma, k_tau);
  } else if (method == "auto" || method == "analytic") {
    lp = ss_texg_stop_success_lpdf_autodisp(SSD, pars, min_ll, upper, n_nodes,
                                            k_sigma, k_tau);
  } else if (method == "live") {
    lp = ss_texg_stop_success_lpdf_live(SSD, pars, min_ll, upper, max_subdiv,
                                        abs_tol, rel_tol, k_sigma, k_tau);
  } else {
    lp = ss_texg_stop_success_lpdf(SSD, pars, min_ll, upper, max_subdiv,
                                   abs_tol, rel_tol, k_sigma, k_tau);
  }
  return (lp <= min_ll) ? 0.0 : std::exp(lp);
}

#endif

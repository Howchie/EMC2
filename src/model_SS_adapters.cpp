#include "model_SS_adapters.h"
#include "utility_functions.h"
#include "race_contract.h"
#include "utils.h"
#include "col_registry.h"
#include "ss_raw.h"
#include "ss_exg_analytic.h"
#include "model_SS_EXG.h"
#include "model_SS_RDEX.h"
#include "gsl_utils.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace Rcpp;
// SS helper pointer types
using ss_go_pdf_fn = NumericVector (*)(NumericVector, NumericMatrix, LogicalVector, double);
using ss_stop_surv_fn = double (*)(double, NumericMatrix);
using ss_stop_success_fn = double (*)(double, NumericMatrix, double, double, int, double, double, double, double);

// Model-specific stop survivor wrappers (columns per src/col_registry.h)
static inline double stop_logsurv_texg_fn(double q, NumericMatrix P) {
  namespace tc = emc2col::ss_texg;
  return ptexg(q, P(0, tc::muS), P(0, tc::sigmaS), P(0, tc::tauS),
               P(0, tc::exgS_lb), R_PosInf, false, true);
}
static inline double stop_logsurv_rdex_fn(double q, NumericMatrix P) {
  namespace rc = emc2col::ss_rdex;
  return ptexg(q, P(0, rc::muS), P(0, rc::sigmaS), P(0, rc::tauS),
               P(0, rc::exgS_lb), R_PosInf, false, true);
}

static inline double sum_log_terms(const NumericVector& x) {
  double out = 0.0;
  for (int i = 0; i < x.size(); ++i) out += x[i];
  return out;
}

static inline double ss_log_go_density_for_winner(
    double rt,
    const NumericMatrix& P,
    const LogicalVector& is_go,
    int winner_idx,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    double min_ll
) {
  const int n_acc = P.nrow();
  NumericVector rt_go(n_acc, rt);
  LogicalVector go_win_mask(n_acc, false);
  LogicalVector go_loss_mask(n_acc, false);
  int go_loss_count = 0;
  for (int i = 0; i < n_acc; ++i) {
    const bool go_i = (is_go[i] == TRUE);
    if (!go_i) continue;
    if (i == winner_idx) go_win_mask[i] = true;
    else {
      go_loss_mask[i] = true;
      ++go_loss_count;
    }
  }

  NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
  double out = (lw.size() > 0) ? sum_log_terms(lw) : R_NegInf;
  if (!R_FINITE(out)) return R_NegInf;

  if (go_loss_count == 0) return out;

  NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
  return out + sum_log_terms(ls);
}

static inline double ss_log_st_density_for_winner(
    double rt,
    double SSD,
    const NumericMatrix& P,
    const LogicalVector& is_st,
    int winner_idx,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    double min_ll
) {
  const int n_acc = P.nrow();
  double rt_st_val = rt - SSD;
  if (rt_st_val < 0.0) rt_st_val = 0.0;
  NumericVector rt_st(n_acc, rt_st_val);
  LogicalVector st_win_mask(n_acc, false);
  LogicalVector st_loss_mask(n_acc, false);
  int st_loss_count = 0;
  for (int i = 0; i < n_acc; ++i) {
    const bool st_i = (is_st[i] == TRUE);
    if (!st_i) continue;
    if (i == winner_idx) st_win_mask[i] = true;
    else {
      st_loss_mask[i] = true;
      ++st_loss_count;
    }
  }

  NumericVector lw = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
  double out = (lw.size() > 0) ? sum_log_terms(lw) : R_NegInf;
  if (!R_FINITE(out)) return R_NegInf;

  if (st_loss_count == 0) return out;

  NumericVector ls = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
  return out + sum_log_terms(ls);
}

static inline double ss_trial_log_response_density(
    double rt,
    const NumericMatrix& P,
    const NumericVector& lR_trial,
    const LogicalVector& is_go,
    const LogicalVector& is_st,
    double SSD,
    bool stop_signal_presented,
    double tf,
    double gf,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    double min_ll,
    int response_code = NA_INTEGER
) {
  const int n_acc = P.nrow();
  const bool response_known = response_code != NA_INTEGER;
  int n_accG = 0;
  int n_accST = 0;
  for (int i = 0; i < n_acc; ++i) {
    if (is_go[i] == TRUE) ++n_accG;
    if (is_st[i] == TRUE) ++n_accST;
  }

  double log_go_any = R_NegInf;
  for (int i = 0; i < n_acc; ++i) {
    if (is_go[i] != TRUE) continue;
    if (response_known && lR_trial[i] != response_code) continue;

    const double go_lprob = ss_log_go_density_for_winner(
      rt, P, is_go, i, go_lpdf_ptr, go_lccdf_ptr, min_ll
    );
    if (!R_FINITE(go_lprob)) continue;

    double log_term = log1m(gf) + go_lprob;
    if (stop_signal_presented && rt > SSD) {
      double rt_eff = rt - SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = stop_logsurv_ptr(rt_eff, P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;

      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        NumericVector rt_st(n_acc, rt_eff);
        LogicalVector st_loss_mask(n_acc, false);
        int st_loss_count = 0;
        for (int j = 0; j < n_acc; ++j) {
          st_loss_mask[j] = (is_st[j] == TRUE);
          if (st_loss_mask[j] == TRUE) ++st_loss_count;
        }
        if (st_loss_count > 0) {
          NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
          st_loss_sum = sum_log_terms(ls_st);
        }
      }

      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      log_term = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
    }

    log_go_any = log_sum_exp(log_go_any, log_term);
  }

  if (!stop_signal_presented || n_accST == 0) return log_go_any;

  double log_st_any = R_NegInf;
  NumericMatrix P_go = submat_rcpp(P, is_go);
  for (int i = 0; i < n_acc; ++i) {
    if (is_st[i] != TRUE) continue;
    if (response_known && lR_trial[i] != response_code) continue;
    if (rt <= SSD) continue;

    const double st_base = ss_log_st_density_for_winner(
      rt, SSD, P, is_st, i, go_lpdf_ptr, go_lccdf_ptr, min_ll
    );
    if (!R_FINITE(st_base)) continue;

    double go_loss_sum = 0.0;
    if (n_accG > 0) {
      NumericVector rt_go(n_acc, rt);
      NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
      go_loss_sum = sum_log_terms(ls_go);
    }

    double rt_eff = rt - SSD;
    if (rt_eff < 0.0) rt_eff = 0.0;
    double log_pstop = stop_success_ptr(SSD, P_go, min_ll, rt_eff,
                                        100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
    if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

    const double term_gf = std::log(gf) + st_base;
    const double term_stop_win = log1m(gf) + log_pstop + st_base;
    const double term_stop_lose =
      log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
    const double log_term =
      log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));

    log_st_any = log_sum_exp(log_st_any, log_term);
  }

  return log_sum_exp(log_go_any, log_st_any);
}

struct SsLcIntegrandData {
  NumericMatrix P;
  NumericVector lR_trial;
  LogicalVector is_go;
  LogicalVector is_st;
  double SSD;
  bool stop_signal_presented;
  double tf;
  double gf;
  ss_go_pdf_fn go_lpdf_ptr;
  ss_go_pdf_fn go_lccdf_ptr;
  ss_stop_surv_fn stop_logsurv_ptr;
  ss_stop_success_fn stop_success_ptr;
  double min_ll;
  int response_code;
};

static double ss_lc_integrand(double x, void* params) {
  SsLcIntegrandData* d = static_cast<SsLcIntegrandData*>(params);
  const double logf = ss_trial_log_response_density(
    x, d->P, d->lR_trial, d->is_go, d->is_st, d->SSD, d->stop_signal_presented,
    d->tf, d->gf, d->go_lpdf_ptr, d->go_lccdf_ptr, d->stop_logsurv_ptr,
    d->stop_success_ptr, d->min_ll, d->response_code
  );
  if (!R_FINITE(logf)) return 0.0;
  return std::exp(logf);
}

static inline double ss_integrate_lc_response_mass(
    double lower,
    double upper,
    const NumericMatrix& P,
    const NumericVector& lR_trial,
    const LogicalVector& is_go,
    const LogicalVector& is_st,
    double SSD,
    bool stop_signal_presented,
    double tf,
    double gf,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    double min_ll,
    int response_code = NA_INTEGER
) {
  if (!(upper > lower)) return min_ll;

  SsLcIntegrandData data{
    P, lR_trial, is_go, is_st, SSD, stop_signal_presented, tf, gf,
    go_lpdf_ptr, go_lccdf_ptr, stop_logsurv_ptr, stop_success_ptr, min_ll,
    response_code
  };

  gsl_function F;
  F.function = &ss_lc_integrand;
  F.params = &data;

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, 200);
  double res = 0.0;
  double err = 0.0;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags(&F, lower, upper, 1e-8, 1e-6, 200, workspace, &res, &err);
  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !R_FINITE(res) || res <= 0.0) return min_ll;
  return std::log(res);
}
double c_log_likelihood_ss(
    NumericMatrix pars,
    DataFrame data,
    const int n_trials,
    IntegerVector expand,
    double min_ll,
    LogicalVector is_ok,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    int idx_tf,
    int idx_gf
) {

  
  // initialise local variables
  const int n_out = expand.length();
  if (is_true(all(!is_ok))) {
    return static_cast<double>(n_out) * min_ll;
  }
  NumericVector lls(n_trials);
  // extract data
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector LT = get_col_with_default(data, "LT", 0.0);
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? as<IntegerVector>(data["lI"]) : IntegerVector(lR.size(), 2);
  
  // dimensional expectations: pars has one row per accumulator per trial
  
  // compute log likelihoods (generalized, matching R's log_likelihood_race_ss)
  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();
  
  NumericVector tt(n_acc);
  auto log_surv_mask = [&](double t, const NumericMatrix& Pcur,
                           const LogicalVector& mask) -> double {
                             tt.fill(t);
                             NumericVector ls = go_lccdf_ptr(tt, Pcur, mask, min_ll);
                             double out = 0.0;
                             for (int i = 0; i < ls.size(); ++i) {
                               out += ls[i];
                             }
                             return out;
                           };
  // n_trials equals data rows grouped by accumulators
  for (int trial = 0; trial < n_trials; trial++) {
    if (is_ok[trial * n_acc] != 1) { lls[trial] = min_ll; continue; }
    
    int start_row = trial * n_acc;
    int end_row   = (trial + 1) * n_acc - 1;
    // basic bounds are guaranteed by correct n_trials passed into this function
    NumericMatrix P = pars(Range(start_row, end_row), _);
    NumericVector lR_trial = lR[Range(start_row, end_row)];
    IntegerVector lI_trial = lI[Range(start_row, end_row)];
    LogicalVector is_go(n_acc, true), is_st(n_acc, false);
    // determine go/ST accumulators if present
    if (has_lI) {
      int go_code = max(lI_trial);
      for (int i = 0; i < n_acc; i++) {
        is_go[i] = (lI_trial[i] == go_code);
        is_st[i] = !is_go[i];
      }
    } else {
      for (int i = 0; i < n_acc; i++) is_st[i] = false;
    }
    int n_accG = sum(is_go);
    int n_accST = sum(is_st);
    
    double tf = P(0, idx_tf);
    double gf = P(0, idx_gf);
    
    double rt = RT[start_row];
    bool response_observed = R[start_row] != NA_INTEGER;
    bool stop_signal_presented = emc2_isfinite(SSD[start_row]);
    // Added UC handling
    double uc = UC[start_row];
    // Identify whether observed response is GO or ST (when response observed)
    bool response_is_go = false;
    if (response_observed) {
      int r_obs = R[start_row];
      for (int i = 0; i < n_acc; i++) {
        if (lR[start_row + i] == r_obs) {
          response_is_go = is_go[i];
          break;
        }
      }
    }
    
    // Build rt vectors for go and st contexts
    NumericVector rt_go(n_acc, rt);
    double rt_st_val = rt - SSD[start_row];
    if (rt_st_val < 0.0) rt_st_val = 0.0;
    NumericVector rt_st(n_acc, rt_st_val);
    
    // GO masks for current trial
    LogicalVector win_mask = winner[Range(start_row, end_row)];
    LogicalVector go_win_mask(n_acc); // winner and go
    LogicalVector go_loss_mask(n_acc);
    for (int i = 0; i < n_acc; i++) {
      go_win_mask[i] = (win_mask[i] && is_go[i]);
      go_loss_mask[i] = (!win_mask[i] && is_go[i]);
    }

    bool lower_censored = (rt == R_NegInf);
    if (lower_censored) {
      const double upper = LC[start_row];
      lls[trial] = ss_integrate_lc_response_mass(
        LT[start_row], upper, P, lR_trial, is_go, is_st, SSD[start_row], stop_signal_presented,
        tf, gf, go_lpdf_ptr, go_lccdf_ptr, stop_logsurv_ptr, stop_success_ptr,
        min_ll, response_observed ? R[start_row] : NA_INTEGER
      );
      continue;
    }
    
    if (!response_observed) {
      bool has_deadline = R_FINITE(uc) && !Rcpp::NumericVector::is_na(uc);
      if (!stop_signal_presented) {
        // No response
        if (!has_deadline) {
          // intrinsic NR
          lls[trial] = std::log(gf);
        } else {
          // deadline-censored GO: gf + (1-gf)*S_go(UC)
          // GO trial, no response by UC: gf + (1-gf)*S_go(UC)
          double logS_go = (n_accG > 0) ? log_surv_mask(uc, P, is_go) : 0.0; // product of all go accumulators not having finished
          lls[trial] = log_sum_exp(std::log(gf), log1m(gf) + logS_go);
        }
        continue;
      }
      // stop trial, no response observed
      if (!has_deadline) {
        if (n_accST == 0) {
          // Stop trial, no ST accumulators: gf + (1-gf)*(1-tf)*pStop
          NumericMatrix P_go = submat_rcpp(P, is_go);
          double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, R_PosInf,
                                              100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
          if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
          double comp1 = std::log(gf);
          double comp2 = log1m(gf) + log1m(tf) + log_pstop;
          lls[trial] = log_sum_exp(comp1, comp2);
        } else {
          // With an ST accumulator present, a stop-trial non-response (with no deadline)
          // can only happen via trigger-failure AND go-failure: tf * gf.
          lls[trial] = std::log(gf) + std::log(tf);
        }
        continue;
      }
      
      // Deadline-censored stop trial: "stop win by UC OR Go/St unfinished by UC"
      // GO survivor by UC
      double logS_go = (n_accG > 0) ? log_surv_mask(uc, P, is_go) : 0.0;
      
      // stop survivor by UC (duration scale: UC-SSD, clamp at 0)
      double uc_eff = uc-SSD[start_row];
      if (!R_FINITE(uc_eff) || uc_eff <= 0.0) uc_eff = 0.0;
      bool stop_can_act = R_FINITE(uc_eff) && (uc_eff > 0.0);
      
      double log_pstop = R_NegInf;
      double logS_stop = 0.0; // log(1)
      
      // pStop(UC): stop finishes before GO and before UC-SSD
      // integral from SSD to UC (integrand function adds SSD to go accumulators so we use uc here)
      NumericMatrix P_go = submat_rcpp(P, is_go);
      if (stop_can_act) {
        logS_stop = stop_logsurv_ptr(uc_eff, P);
        log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, uc_eff, 100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      }
      
      // Triggered, go-not-failed core no-response by UC:
      // p = pStop(UC) + S_go(UC)*S_stop(UC-SSD)
      double log_core_trig = log_sum_exp(log_pstop, logS_go + logS_stop);
      
      if (n_accST == 0) {
        // No ST accumulators:
        // p = gf + (1-gf) * [ tf*S_go + (1-tf)*core_trig ]
        double log_no_nogf = log_mix(tf, logS_go, log_core_trig);
        lls[trial] = log_sum_exp(std::log(gf), log1m(gf) + log_no_nogf);
        continue;
      }
      
      // Stop-triggered case - 
      double logS_st = log_surv_mask(uc_eff, P, is_st);
      double log_trig = logS_st + log_sum_exp(std::log(gf), log1m(gf) + log_core_trig);
      
      double log_tfbranch = log_sum_exp(std::log(gf), log1m(gf) + logS_go);
      
      lls[trial] = log_mix(tf, log_tfbranch, log_trig);
      continue;
    }
    
    // Response observed
    if (!stop_signal_presented) {
      // GO trial with response: (1-gf) * GO race ll
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : R_NegInf;
      if (!R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) {
          go_lprob += ls[i];
        }
      }
      lls[trial] = log1m(gf) + go_lprob;
      continue;
    }
    
    // Stop trial with response
    if (response_is_go) {
      // GO wins on stop trial: (1-gf) * [ tf * go + (1-tf) * (go + stop_surv + st_loss) ]
      // go race ll at observed rt
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : R_NegInf;
      if (!R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) {
          go_lprob += ls[i];
        }
      }
      // stop survivor at observed rt (rt - SSD)
      double rt_eff = rt - SSD[start_row];
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = stop_logsurv_ptr(rt_eff, P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;
      // ST losers survivors (if any)
      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (is_st[i] && !win_mask[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) {
          st_loss_sum += ls_st[i];
        }
      }
      double comp_tf = go_lprob; // only go
      double comp_notf = go_lprob + log_stop_surv + st_loss_sum; // fair race (stop loses)
      lls[trial] = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
      continue;
    } else {
      // ST wins on stop trial
      // ST winner log pdf at rt - SSD
      LogicalVector st_win_mask(n_acc);
      for (int i = 0; i < n_acc; i++) st_win_mask[i] = (win_mask[i] && is_st[i]);
      NumericVector lw_st = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
      double st_winner_logpdf = (lw_st.size() > 0) ? sum(lw_st) : R_NegInf;
      if (!R_FINITE(st_winner_logpdf)) st_winner_logpdf = R_NegInf;
      // ST losers survivors
      double st_loss_sum = 0.0;
      if (n_accST > 1) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (!win_mask[i] && is_st[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) {
          st_loss_sum += ls_st[i];
        }
      }
      // GO losers survivors
      double go_loss_sum = 0.0;
      if (n_accG > 0) {
        NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
        for (int i = 0; i < ls_go.size(); i++) {
          go_loss_sum += ls_go[i];
        }
      }
      // Stop success probability up to observed rt (only go racers influence integral)
      NumericMatrix P_go = submat_rcpp(P, is_go);
      // stop survivor at observed rt (rt - SSD)
      double rt_eff = rt - SSD[start_row];
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, rt_eff,
                                          100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
      if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      
      double st_base = st_winner_logpdf + st_loss_sum;
      double term_gf = std::log(gf) + st_base; // go failure -> only ST race
      double term_stop_win = log1m(gf) + log_pstop + st_base; // stop beats go -> only ST race
      double term_stop_lose = log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum; // all race, no stop win
      lls[trial] = log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
      continue;
    }
  }
  lls[is_na(lls)] = min_ll;
  lls[is_infinite(lls)] = min_ll;
  lls[lls < min_ll] = min_ll;
  double sum_ll = 0.0;
  for (int i = 0; i < n_out; ++i) {
    sum_ll += lls[expand[i] - 1]; // expand created in 1-based R
  }
  return sum_ll;
}

// ============================================================================
// Raw-buffer stop-signal path.
//
// c_log_likelihood_ss above re-derived the whole trial structure from Rcpp
// columns for every trial of every particle (Range() submatrix copies,
// LogicalVector masks, NumericVector rt broadcasts). Everything data-fixed is
// now computed ONCE per calc_ll_oo call into SsSharedState, and the
// per-particle evaluator reads ParamTable base columns directly through the
// canonical kernel blocks defined in ss_raw.h (see SsRawModel). Numerically
// identical to c_log_likelihood_ss, which remains as the fallback (and
// readable reference) when the ParamTable column mapping is unavailable.
// ============================================================================
SsSharedState build_ss_shared_state(const Rcpp::DataFrame& data,
                                           int n_trials_ss, int n_acc) {
  SsSharedState ss;
  ss.n_trials = n_trials_ss;
  ss.n_acc = n_acc;
  if (n_trials_ss <= 0 || n_acc <= 0) return ss;
  if (!data.containsElementNamed("rt") || !data.containsElementNamed("R") ||
      !data.containsElementNamed("SSD") || !data.containsElementNamed("lR") ||
      !data.containsElementNamed("winner")) {
    return ss;
  }

  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  IntegerVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector LT = get_col_with_default(data, "LT", 0.0);
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  const bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? Rcpp::as<IntegerVector>(data["lI"])
                            : IntegerVector(lR.size(), 2);
  if (RT.size() != n_trials_ss * n_acc) return ss;

  const int n_rows = n_trials_ss * n_acc;
  ss.kind.resize(n_trials_ss);
  ss.rt.resize(n_trials_ss); ss.SSD.resize(n_trials_ss);
  ss.UC.resize(n_trials_ss); ss.LC.resize(n_trials_ss); ss.LT.resize(n_trials_ss);
  ss.resp_code.resize(n_trials_ss);
  ss.n_go.resize(n_trials_ss); ss.n_st.resize(n_trials_ss);
  ss.first_go_row.resize(n_trials_ss);
  ss.is_go_row.resize(n_rows); ss.is_st_row.resize(n_rows); ss.win_row.resize(n_rows);
  ss.lR_row.resize(n_rows);

  for (int trial = 0; trial < n_trials_ss; ++trial) {
    const int start = trial * n_acc;
    const double rt = RT[start];
    const double ssd = SSD[start];
    const double uc = UC[start];
    const int r_obs = R[start];
    const bool response_observed = (r_obs != NA_INTEGER);
    const bool stop_presented = emc2_isfinite(ssd);
    const bool has_deadline = R_FINITE(uc) && !Rcpp::NumericVector::is_na(uc);

    ss.rt[trial] = rt; ss.SSD[trial] = ssd; ss.UC[trial] = uc;
    ss.LC[trial] = LC[start]; ss.LT[trial] = LT[start];
    ss.resp_code[trial] = r_obs;

    // go/ST partition exactly as before: within-trial max lI code marks go
    int go_code = lI[start];
    for (int i = 1; i < n_acc; ++i) {
      if (lI[start + i] > go_code) go_code = lI[start + i];
    }
    int ng = 0, nst = 0, first_go = -1;
    for (int i = 0; i < n_acc; ++i) {
      const bool go_i = has_lI ? (lI[start + i] == go_code) : true;
      ss.is_go_row[start + i] = go_i ? 1 : 0;
      ss.is_st_row[start + i] = go_i ? 0 : 1;
      ss.win_row[start + i] = (winner[start + i] == TRUE) ? 1 : 0;
      ss.lR_row[start + i] = lR[start + i];
      if (go_i) { ++ng; if (first_go < 0) first_go = i; }
      else ++nst;
    }
    ss.n_go[trial] = ng; ss.n_st[trial] = nst;
    ss.first_go_row[trial] = (first_go < 0) ? 0 : first_go;

    bool response_is_go = false;
    if (response_observed) {
      for (int i = 0; i < n_acc; ++i) {
        if (lR[start + i] == r_obs) { response_is_go = ss.is_go_row[start + i] != 0; break; }
      }
    }

    uint8_t k;
    if (rt == R_NegInf) {
      k = SS_TRIAL_LC_CENSORED;
    } else if (!response_observed) {
      if (!stop_presented) k = has_deadline ? SS_TRIAL_GO_NR_DEADLINE : SS_TRIAL_GO_NR_FREE;
      else                 k = has_deadline ? SS_TRIAL_STOP_NR_DEADLINE : SS_TRIAL_STOP_NR_FREE;
    } else if (!stop_presented) {
      k = SS_TRIAL_GO_RESP;
    } else {
      k = response_is_go ? SS_TRIAL_STOP_GO_RESP : SS_TRIAL_STOP_ST_RESP;
    }
    ss.kind[trial] = k;
  }

  ss.valid = true;
  return ss;
}

// Per-trial-per-particle evaluation context: canonical parameter blocks plus
// the data-fixed masks for this trial. All raw pointers; nothing allocates.
struct SsTrialEvalCtx {
  const SsRawModel* M;
  const double* acc;      // n_acc * SS_ACC_STRIDE canonical blocks (trial order)
  const double* acc_go;   // n_go * SS_ACC_STRIDE, go rows compacted in order
  const uint8_t* is_go;   // per trial row
  const uint8_t* is_st;
  const int* lR_codes;
  int n_acc, n_go, n_st;
  double stop_r0[4];      // stop block from trial row 0   (stop survivor)
  double stop_g0[4];      // stop block from first go row  (stop-success integral)
  double tf, gf;
  double SSD;
  bool stop_presented;
  int response_code;      // NA_INTEGER when unconstrained (LC integration)
  double min_ll;
};

static inline SsStopCtx ss_make_stop_ctx(const SsTrialEvalCtx& c) {
  SsStopCtx sc;
  sc.SSD = c.SSD;
  sc.acc_go = c.acc_go;
  sc.n_go = c.n_go;
  sc.muS = c.stop_g0[0]; sc.sigS = c.stop_g0[1];
  sc.tauS = c.stop_g0[2]; sc.lbS = c.stop_g0[3];
  sc.acc_surv = c.M->acc_surv;
  return sc;
}

// log survivor of the stop process at duration t (raw twin of
// stop_logsurv_texg_fn / stop_logsurv_rdex_fn: both models race a TEXG stop).
static inline double ss_stop_logsurv_raw(double t, const double* stop4) {
  return ptexg(t, stop4[0], stop4[1], stop4[2], stop4[3], R_PosInf, false, true);
}

// Sum of protected log survivors over a masked subset of accumulators, at
// time t. Matches the texg_go_lccdf/rdex_go_lccdf + sum pattern: each element
// is floored to -Inf when non-finite.
static inline double ss_masked_lsurv_sum_raw(const SsTrialEvalCtx& c, double t,
                                             const uint8_t* mask) {
  double out = 0.0;
  for (int i = 0; i < c.n_acc; ++i) {
    if (!mask[i]) continue;
    const double v = c.M->acc_lsurv(t, c.acc + SS_ACC_STRIDE * i);
    out += emc2_isfinite(v) ? v : R_NegInf;
  }
  return out;
}

// Full trial response density at rt, marginalizing over candidate winners.
// Raw twin of ss_trial_log_response_density (used by the LC-censoring path).
static double ss_trial_log_response_density_raw(double rt, const SsTrialEvalCtx& c) {
  const SsRawModel& M = *c.M;
  const bool response_known = (c.response_code != NA_INTEGER);

  double log_go_any = R_NegInf;
  for (int i = 0; i < c.n_acc; ++i) {
    if (!c.is_go[i]) continue;
    if (response_known && c.lR_codes[i] != c.response_code) continue;

    // winner-i go race density (twin of ss_log_go_density_for_winner)
    const double lw = M.acc_lpdf(rt, c.acc + SS_ACC_STRIDE * i);
    double go_lprob = emc2_isfinite(lw) ? lw : R_NegInf;
    if (!R_FINITE(go_lprob)) continue;
    for (int j = 0; j < c.n_acc; ++j) {
      if (!c.is_go[j] || j == i) continue;
      const double ls = M.acc_lsurv(rt, c.acc + SS_ACC_STRIDE * j);
      go_lprob += emc2_isfinite(ls) ? ls : R_NegInf;
    }
    if (!R_FINITE(go_lprob)) continue;

    double log_term = log1m(c.gf) + go_lprob;
    if (c.stop_presented && rt > c.SSD) {
      double rt_eff = rt - c.SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = ss_stop_logsurv_raw(rt_eff, c.stop_r0);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = c.min_ll;
      double st_loss_sum = 0.0;
      if (c.n_st > 0) {
        st_loss_sum = ss_masked_lsurv_sum_raw(c, rt_eff, c.is_st);
      }
      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      log_term = log1m(c.gf) + log_mix(c.tf, comp_tf, comp_notf);
    }
    log_go_any = log_sum_exp(log_go_any, log_term);
  }

  if (!c.stop_presented || c.n_st == 0) return log_go_any;

  double log_st_any = R_NegInf;
  for (int i = 0; i < c.n_acc; ++i) {
    if (!c.is_st[i]) continue;
    if (response_known && c.lR_codes[i] != c.response_code) continue;
    if (rt <= c.SSD) continue;

    double rt_st = rt - c.SSD;
    if (rt_st < 0.0) rt_st = 0.0;
    // winner-i ST race density (twin of ss_log_st_density_for_winner)
    const double lw = M.acc_lpdf(rt_st, c.acc + SS_ACC_STRIDE * i);
    double st_base = emc2_isfinite(lw) ? lw : R_NegInf;
    if (!R_FINITE(st_base)) continue;
    for (int j = 0; j < c.n_acc; ++j) {
      if (!c.is_st[j] || j == i) continue;
      const double ls = M.acc_lsurv(rt_st, c.acc + SS_ACC_STRIDE * j);
      st_base += emc2_isfinite(ls) ? ls : R_NegInf;
    }
    if (!R_FINITE(st_base)) continue;

    double go_loss_sum = 0.0;
    if (c.n_go > 0) {
      go_loss_sum = ss_masked_lsurv_sum_raw(c, rt, c.is_go);
    }
    double rt_eff = rt - c.SSD;
    if (rt_eff < 0.0) rt_eff = 0.0;
    const SsStopCtx sc = ss_make_stop_ctx(c);
    double log_pstop = ss_stop_success_raw_live(sc, c.min_ll, rt_eff,
                                                M.clamp_empty_window,
                                                M.try_analytic1);
    if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

    const double term_gf = std::log(c.gf) + st_base;
    const double term_stop_win = log1m(c.gf) + log_pstop + st_base;
    const double term_stop_lose =
      log1m(c.gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
    const double log_term =
      log1m(c.tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));

    log_st_any = log_sum_exp(log_st_any, log_term);
  }

  return log_sum_exp(log_go_any, log_st_any);
}

static double ss_lc_integrand_raw(double x, void* params) {
  const SsTrialEvalCtx* c = static_cast<const SsTrialEvalCtx*>(params);
  const double logf = ss_trial_log_response_density_raw(x, *c);
  if (!R_FINITE(logf)) return 0.0;
  return std::exp(logf);
}

// Raw twin of ss_integrate_lc_response_mass (same GSL setup/tolerances).
static inline double ss_integrate_lc_response_mass_raw(double lower, double upper,
                                                       const SsTrialEvalCtx& c) {
  if (!(upper > lower)) return c.min_ll;

  gsl_function F;
  F.function = &ss_lc_integrand_raw;
  F.params = const_cast<SsTrialEvalCtx*>(&c);

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, 200);
  double res = 0.0;
  double err = 0.0;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags(&F, lower, upper, 1e-8, 1e-6, 200, workspace, &res, &err);
  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !R_FINITE(res) || res <= 0.0) return c.min_ll;
  return std::log(res);
}
// Raw-path per-particle stop-signal log likelihood. `cols` are pointers to the
// ParamTable base columns in p_types order (the same layout the materialized
// matrix had). Branch bodies mirror c_log_likelihood_ss trial by trial.
double c_log_likelihood_ss_pt(
    const double* const* cols,
    const SsSharedState& ss,
    const SsRawModel& M,
    const LogicalVector& is_ok,
    const IntegerVector& expand,
    double min_ll,
    SsRawWorkspace& ws
) {
  const int n_out = expand.length();
  if (is_true(all(!is_ok))) {
    return static_cast<double>(n_out) * min_ll;
  }
  const int n_acc = ss.n_acc;
  const int n_trials = ss.n_trials;
  ws.acc.resize(static_cast<size_t>(SS_ACC_STRIDE) * n_acc);
  ws.acc_go.resize(static_cast<size_t>(SS_ACC_STRIDE) * n_acc);
  ws.lls.assign(static_cast<size_t>(n_trials), min_ll);

  for (int trial = 0; trial < n_trials; ++trial) {
    const int start = trial * n_acc;
    if (is_ok[start] != TRUE) { ws.lls[trial] = min_ll; continue; }

    // canonical parameter blocks for this trial/particle
    double* acc = ws.acc.data();
    for (int i = 0; i < n_acc; ++i) {
      M.fill_acc(cols, start + i, acc + SS_ACC_STRIDE * i);
    }
    double* acc_go = ws.acc_go.data();
    int k = 0;
    for (int i = 0; i < n_acc; ++i) {
      if (!ss.is_go_row[start + i]) continue;
      std::memcpy(acc_go + SS_ACC_STRIDE * k, acc + SS_ACC_STRIDE * i,
                  SS_ACC_STRIDE * sizeof(double));
      ++k;
    }

    SsTrialEvalCtx c;
    c.M = &M;
    c.acc = acc;
    c.acc_go = acc_go;
    c.is_go = ss.is_go_row.data() + start;
    c.is_st = ss.is_st_row.data() + start;
    c.lR_codes = ss.lR_row.data() + start;
    c.n_acc = n_acc;
    c.n_go = ss.n_go[trial];
    c.n_st = ss.n_st[trial];
    M.fill_stop(cols, start, c.stop_r0);
    M.fill_stop(cols, start + ss.first_go_row[trial], c.stop_g0);
    c.tf = cols[M.idx_tf][start];
    c.gf = cols[M.idx_gf][start];
    c.SSD = ss.SSD[trial];
    c.stop_presented = emc2_isfinite(c.SSD);
    c.response_code = ss.resp_code[trial];
    c.min_ll = min_ll;

    const uint8_t* win = ss.win_row.data() + start;
    const double rt = ss.rt[trial];
    double ll = min_ll;

    switch (ss.kind[trial]) {

    case SS_TRIAL_LC_CENSORED: {
      ll = ss_integrate_lc_response_mass_raw(ss.LT[trial], ss.LC[trial], c);
      break;
    }

    case SS_TRIAL_GO_NR_FREE: {
      ll = std::log(c.gf);
      break;
    }

    case SS_TRIAL_GO_NR_DEADLINE: {
      const double logS_go = (c.n_go > 0)
        ? ss_masked_lsurv_sum_raw(c, ss.UC[trial], c.is_go) : 0.0;
      ll = log_sum_exp(std::log(c.gf), log1m(c.gf) + logS_go);
      break;
    }

    case SS_TRIAL_STOP_NR_FREE: {
      if (c.n_st == 0) {
        const SsStopCtx sc = ss_make_stop_ctx(c);
        double log_pstop = ss_stop_success_raw_live(sc, min_ll, R_PosInf,
                                                    M.clamp_empty_window,
                                                    M.try_analytic1);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
        ll = log_sum_exp(std::log(c.gf),
                         log1m(c.gf) + log1m(c.tf) + log_pstop);
      } else {
        // NR with an ST accumulator and no deadline: tf AND gf must both fail
        ll = std::log(c.gf) + std::log(c.tf);
      }
      break;
    }

    case SS_TRIAL_STOP_NR_DEADLINE: {
      const double uc = ss.UC[trial];
      const double logS_go = (c.n_go > 0)
        ? ss_masked_lsurv_sum_raw(c, uc, c.is_go) : 0.0;

      double uc_eff = uc - c.SSD;
      if (!R_FINITE(uc_eff) || uc_eff <= 0.0) uc_eff = 0.0;
      const bool stop_can_act = R_FINITE(uc_eff) && (uc_eff > 0.0);

      double log_pstop = R_NegInf;
      double logS_stop = 0.0;   // log(1)
      if (stop_can_act) {
        logS_stop = ss_stop_logsurv_raw(uc_eff, c.stop_r0);
        const SsStopCtx sc = ss_make_stop_ctx(c);
        log_pstop = ss_stop_success_raw_live(sc, min_ll, uc_eff,
                                             M.clamp_empty_window,
                                             M.try_analytic1);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      }

      const double log_core_trig = log_sum_exp(log_pstop, logS_go + logS_stop);

      if (c.n_st == 0) {
        const double log_no_nogf = log_mix(c.tf, logS_go, log_core_trig);
        ll = log_sum_exp(std::log(c.gf), log1m(c.gf) + log_no_nogf);
      } else {
        const double logS_st = ss_masked_lsurv_sum_raw(c, uc_eff, c.is_st);
        const double log_trig =
          logS_st + log_sum_exp(std::log(c.gf), log1m(c.gf) + log_core_trig);
        const double log_tfbranch =
          log_sum_exp(std::log(c.gf), log1m(c.gf) + logS_go);
        ll = log_mix(c.tf, log_tfbranch, log_trig);
      }
      break;
    }

    case SS_TRIAL_GO_RESP:
    case SS_TRIAL_STOP_GO_RESP: {
      // go race density with the observed (data) winner
      double go_lprob = 0.0;
      bool any_win = false;
      for (int i = 0; i < n_acc; ++i) {
        if (!(win[i] && c.is_go[i])) continue;
        const double lw = M.acc_lpdf(rt, acc + SS_ACC_STRIDE * i);
        go_lprob += emc2_isfinite(lw) ? lw : R_NegInf;
        any_win = true;
      }
      if (!any_win || !R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (c.n_go > 1) {
        for (int i = 0; i < n_acc; ++i) {
          if (!(!win[i] && c.is_go[i])) continue;
          const double ls = M.acc_lsurv(rt, acc + SS_ACC_STRIDE * i);
          go_lprob += emc2_isfinite(ls) ? ls : R_NegInf;
        }
      }

      if (ss.kind[trial] == SS_TRIAL_GO_RESP) {
        ll = log1m(c.gf) + go_lprob;
        break;
      }

      // stop trial, go response: mixture over trigger failure
      double rt_eff = rt - c.SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = ss_stop_logsurv_raw(rt_eff, c.stop_r0);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;
      double st_loss_sum = 0.0;
      if (c.n_st > 0) {
        for (int i = 0; i < n_acc; ++i) {
          if (!(c.is_st[i] && !win[i])) continue;
          const double ls = M.acc_lsurv(rt_eff, acc + SS_ACC_STRIDE * i);
          st_loss_sum += emc2_isfinite(ls) ? ls : R_NegInf;
        }
      }
      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      ll = log1m(c.gf) + log_mix(c.tf, comp_tf, comp_notf);
      break;
    }

    case SS_TRIAL_STOP_ST_RESP: {
      double rt_st = rt - c.SSD;
      if (rt_st < 0.0) rt_st = 0.0;
      // ST winner log pdf at rt - SSD
      double st_winner_logpdf = 0.0;
      bool any_win = false;
      for (int i = 0; i < n_acc; ++i) {
        if (!(win[i] && c.is_st[i])) continue;
        const double lw = M.acc_lpdf(rt_st, acc + SS_ACC_STRIDE * i);
        st_winner_logpdf += emc2_isfinite(lw) ? lw : R_NegInf;
        any_win = true;
      }
      if (!any_win || !R_FINITE(st_winner_logpdf)) st_winner_logpdf = R_NegInf;
      // ST losers survivors
      double st_loss_sum = 0.0;
      if (c.n_st > 1) {
        for (int i = 0; i < n_acc; ++i) {
          if (!(!win[i] && c.is_st[i])) continue;
          const double ls = M.acc_lsurv(rt_st, acc + SS_ACC_STRIDE * i);
          st_loss_sum += emc2_isfinite(ls) ? ls : R_NegInf;
        }
      }
      // GO losers survivors (all go accumulators, at raw rt)
      double go_loss_sum = 0.0;
      if (c.n_go > 0) {
        go_loss_sum = ss_masked_lsurv_sum_raw(c, rt, c.is_go);
      }
      double rt_eff = rt - c.SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      const SsStopCtx sc = ss_make_stop_ctx(c);
      double log_pstop = ss_stop_success_raw_live(sc, min_ll, rt_eff,
                                                  M.clamp_empty_window,
                                                  M.try_analytic1);
      if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

      const double st_base = st_winner_logpdf + st_loss_sum;
      const double term_gf = std::log(c.gf) + st_base;
      const double term_stop_win = log1m(c.gf) + log_pstop + st_base;
      const double term_stop_lose =
        log1m(c.gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
      ll = log1m(c.tf) +
        log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
      break;
    }
    }

    // same clamp as the vectorized epilogue of c_log_likelihood_ss
    if (!R_FINITE(ll) || ll < min_ll) ll = min_ll;
    ws.lls[trial] = ll;
  }

  double sum_ll = 0.0;
  for (int i = 0; i < n_out; ++i) {
    sum_ll += ws.lls[expand[i] - 1];   // expand created in 1-based R
  }
  return sum_ll;
}
SSModelAdapter resolve_ss_adapter(const std::string& type_std) {
  SSModelAdapter a;
  const bool is_exg = type_std.find("EXG") != std::string::npos;
  a.go_lpdf_ptr      = is_exg ? texg_go_lpdf          : rdex_go_lpdf;
  a.go_lccdf_ptr     = is_exg ? texg_go_lccdf          : rdex_go_lccdf;
  a.stop_logsurv_ptr = is_exg ? stop_logsurv_texg_fn   : stop_logsurv_rdex_fn;
  a.stop_success_ptr = is_exg ? ss_texg_stop_success_lpdf_live : ss_rdex_stop_success_lpdf_live;
  a.idx_tf           = is_exg ? emc2col::ss_texg::tf : emc2col::ss_rdex::tf;
  a.idx_gf           = is_exg ? emc2col::ss_texg::gf : emc2col::ss_rdex::gf;
  return a;
}

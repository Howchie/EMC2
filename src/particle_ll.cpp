#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "model_SS_EXG.h"
#include "model_SS_RDEX.h"
#include "composite_functions.h"
#include "trend.h"
#include "utils.h"
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include <cmath>
#include <string>
#include <memory>
#include <array>

using namespace Rcpp;
//todo check whether all the includes are still required

NumericMatrix get_pars_c_wrapper_oo_core(NumericMatrix particle_matrix,
                                         DataFrame data,
                                         NumericVector constants,
                                         List designs,
                                         List bounds,
                                         List transforms,
                                         List pretransforms,
                                         Rcpp::Nullable<Rcpp::List> trend,
                                         bool return_kernel_matrix,
                                         bool return_all_pars,
                                         IntegerVector kernel_output_codes);

// [[Rcpp::export]]
Rcpp::NumericMatrix do_transform(Rcpp::NumericMatrix pars, Rcpp::List transform) {
  // Build the specs for these parameters
  std::vector<TransformSpec> specs = make_transform_specs(pars, transform);
  // Apply transformation in place and return
  return c_do_transform(pars, specs);
}


NumericMatrix c_map_p(NumericVector p_vector,
                      CharacterVector p_types,
                      List designs,
                      int n_trials,
                      DataFrame data,
                      List trend,
                      const std::vector<TransformSpec>& full_specs,
                      bool drop_trend_pars=true) {
  
  // Extract information about trends
  const bool has_trend = (trend.length() > 0);
  bool premap = false;
  bool pretransform = false;
  CharacterVector trend_names;
  if (has_trend) {
    trend_names = trend.names();
    for (int i = 0; i < trend.size(); ++i) {
      List cur = trend[i];
      std::string ph = Rcpp::as<std::string>(cur["phase"]);
      if (ph == "premap") premap = true;
      if (ph == "pretransform") pretransform = true;
    }
  }
  
  const int n_params = p_types.size();
  NumericMatrix pars(n_trials, n_params);
  colnames(pars) = p_types;

  // Build a fast lookup for available particle coefficients.
  // Some design columns (e.g., contrast coefficients like v_lMd) may be
  // absent from p_vector in certain call paths and should contribute zero.
  std::unordered_map<std::string, double> p_lookup;
  CharacterVector pv_names = p_vector.names();
  if (pv_names.size() > 0) {
    p_lookup.reserve(pv_names.size());
    for (int i = 0; i < pv_names.size(); ++i) {
      p_lookup[Rcpp::as<std::string>(pv_names[i])] = p_vector[i];
    }
  }
  
  // Prepare trend parameter columns when needed
  NumericMatrix trend_pars;
  LogicalVector trend_index(n_params, FALSE);
  CharacterVector trend_pnames;
  if (has_trend && (premap || pretransform)) {
    // Fill in trend columns first so that they can be used in premapped trend
    // This function also applies transformations to the trend parameters
    // to ensure real-lines support.
    // The pre-transform trends are also included here, they are used after map_p
    // but they need to be transformed already (before the other trend parameters)
    // are transformed.
    trend_pars = build_trend_columns_from_design(p_vector, p_types, designs, n_trials, trend, full_specs);
    trend_pnames = colnames(trend_pars);
    trend_index = contains_multiple(p_types, trend_pnames);
  }
  
  // Map non-trend parameters from designs, applying premap trends if requested
  for (int i = 0; i < n_params; i++) {
    if (trend_index[i] == TRUE) continue; // skip trend parameters here
    NumericMatrix cur_design = designs[i];
    CharacterVector cur_names = colnames(cur_design);
    for (int j = 0; j < cur_design.ncol(); j++) {
      std::string cur_name = Rcpp::as<std::string>(cur_names[j]);
      auto it_coef = p_lookup.find(cur_name);
      if (it_coef == p_lookup.end()) continue;
      NumericVector p_mult_design(n_trials, it_coef->second);
      if (has_trend && premap) {
        p_mult_design = apply_premap_trends(data, trend, trend_names, cur_name, p_mult_design, trend_pars, p_vector);
      }
      p_mult_design = p_mult_design * cur_design(_, j);
      LogicalVector bad = is_na(p_mult_design) | is_nan(p_mult_design);
      p_mult_design[bad] = 0;
      pars(_, i) = pars(_, i) + p_mult_design;
    }
  }
  
  // If using pretransform trends, copy the pre-transformed trend cols into pars by name
  if (has_trend && pretransform) {
    // Only fill columns for pretransform entries
    CharacterVector tf_names = collect_trend_param_names_phase(trend, "pretransform");
    NumericMatrix trend_pars_tf = (tf_names.size() > 0) ? submat_rcpp_col_by_names(trend_pars, tf_names) : NumericMatrix(n_trials, 0);
    fill_trend_columns_for_pretransform(pars, p_types, trend_pars_tf);
  }
  
  // If premap, trend parameter columns are not part of the final matrix
  if (has_trend && premap && drop_trend_pars) {
    CharacterVector names_premap = collect_trend_param_names_phase(trend, "premap");
    if (names_premap.size() > 0) {
      pars = submat_rcpp_col(pars, !contains_multiple(p_types, names_premap));
    }
  }
  return pars;
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, const std::vector<PreTransformSpec>& p_specs,
                              CharacterVector p_types, List designs, int n_trials, DataFrame data, List trend,
                              const std::vector<TransformSpec>& full_specs,
                              bool drop_trend_pars=true){
  const bool has_trend = (trend.length() > 0);
  bool pretransform = false;
  bool posttransform = false;
  if (has_trend) {
    for (int i = 0; i < trend.size(); ++i) {
      List cur = trend[i];
      std::string ph = Rcpp::as<std::string>(cur["phase"]);
      if (ph == "pretransform") pretransform = true;
      if (ph == "posttransform") posttransform = true;
    }
  }
  NumericVector p_vector_updtd(clone(p_vector));
  p_vector_updtd = c_do_pre_transform(p_vector_updtd, p_specs);
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, full_specs, drop_trend_pars);
  // // Check if pretransform trend applies
  if(pretransform){
    pars = prep_trend_phase(data, trend, pars, "pretransform");
  }
  std::vector<TransformSpec> t_specs = make_transform_specs_from_full(pars, p_types, full_specs);
  pars = c_do_transform(pars, t_specs);
  // Check if posttransform trend applies
  if(posttransform){
    // Build trend parameter columns once (transformed) and pass override
    NumericMatrix trend_pars_all = build_trend_columns_from_design(p_vector_updtd, p_types, designs, n_trials, trend, full_specs);
    CharacterVector names_post = collect_trend_param_names_phase(trend, "posttransform");
    NumericMatrix trend_pars_post = (names_post.size() > 0) ? submat_rcpp_col_by_names(trend_pars_all, names_post) : NumericMatrix(n_trials, 0);
    pars = prep_trend_phase_with_pars(data, trend, pars, "posttransform", trend_pars_post);
  }
  // ok is calculated afterwards and Ttransform applied in the function
  return(pars);
}


// SS helper pointer types
using ss_go_pdf_fn = NumericVector (*)(NumericVector, NumericMatrix, LogicalVector, double);
using ss_stop_surv_fn = double (*)(double, NumericMatrix);
using ss_stop_success_fn = double (*)(double, NumericMatrix, double, double, int, double, double, double, double);

// Model-specific stop survivor wrappers (read fixed columns)
static inline double stop_logsurv_texg_fn(double q, NumericMatrix P) {
  // EXG stop: muS=3, sigmaS=4, tauS=5, exgS_lb=9
  return ptexg(q, P(0, 3), P(0, 4), P(0, 5), P(0, 9), R_PosInf, false, true);
}
static inline double stop_logsurv_rdex_fn(double q, NumericMatrix P) {
  // RDEX stop: muS=5, sigmaS=6, tauS=7, exgS_lb=10
  return ptexg(q, P(0, 5), P(0, 6), P(0, 7), P(0, 10), R_PosInf, false, true);
}

struct RaceModelAdapter {
  RacePdfFun model_dfun_ptr = nullptr;
  RaceCdfFun model_pfun_ptr = nullptr;
  RacePdf1Fun pdf1_ptr = nullptr;
  RaceCdf1Fun cdf1_ptr = nullptr;
  ContextForRaceModels ctx;
};

static inline RaceModelAdapter resolve_race_model_adapter(const std::string& type_std,
                                                          const std::string& caller) {
  RaceModelAdapter out;
  out.ctx.min_lik_for_pdf = std::exp(std::log(1e-10));
  out.ctx.use_posdrift = true;
  out.ctx.gng = false;

  if (type_std.find("LBA") != std::string::npos) {
    out.model_dfun_ptr = &lba_dfun_adapter;
    out.model_pfun_ptr = &lba_pfun_adapter;
    out.pdf1_ptr = &dlba_scalar;
    out.cdf1_ptr = &plba_scalar;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("RDM") != std::string::npos) {
    out.model_dfun_ptr = &rdm_dfun_adapter;
    out.model_pfun_ptr = &rdm_pfun_adapter;
    out.pdf1_ptr = &drdm_scalar;
    out.cdf1_ptr = &prdm_scalar;
  } else if (type_std.find("LNR") != std::string::npos) {
    out.model_dfun_ptr = &lnr_dfun_adapter;
    out.model_pfun_ptr = &lnr_pfun_adapter;
    out.pdf1_ptr = &dlnr_scalar;
    out.cdf1_ptr = &plnr_scalar;
  } else {
    Rcpp::stop("Unsupported race model type string in %s: %s", caller.c_str(), type_std.c_str());
  }

  if (type_std.find("GNG") != std::string::npos) {
    out.ctx.gng = true;
  }
  return out;
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
    NumericVector lls_expanded(n_out, min_ll);
    return(sum(lls_expanded));
  }
  NumericVector lls(n_trials);
  NumericVector lls_expanded(n_out);
  // extract data
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? as<IntegerVector>(data["lI"]) : IntegerVector(lR.size(), 2);
  
  // dimensional expectations: pars has one row per accumulator per trial
  
  // compute log likelihoods (generalized, matching R's log_likelihood_race_ss)
  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();
  
  auto log_surv_mask = [&](double t, const NumericMatrix& Pcur,
                           const LogicalVector& mask) -> double {
                             NumericVector tt(n_acc);
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
    bool stop_signal_presented = std::isfinite(SSD[start_row]);
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
                                              100, 1e-8, 1e-6, 8.0, 16.0);
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
        log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, uc_eff, 100, 1e-8, 1e-6, 8.0, 16.0);
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
                                          100, 1e-8, 1e-6, 8.0, 16.0);
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
  lls_expanded = c_expand(lls, expand); // decompress
  return(sum(lls_expanded));
}

double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials, IntegerVector expand,
                            double min_ll, LogicalVector is_ok, bool gng){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector LT = get_col_with_default(data, "LT", 0.0);
  NumericVector UT = get_col_with_default(data, "UT", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  NumericVector lls(n_trials);
  NumericVector lls_exp(n_out);
  LogicalVector finite = !is_na(rts) & is_finite(rts);
  LogicalVector ok_finite = is_ok & finite;
  LogicalVector ok_nonfinite = is_ok & !finite;

  auto log_interval_mass = [&](double logF_low, double logF_high) -> double {
    if (!R_FINITE(logF_high)) return R_NegInf;
    if (logF_low == R_NegInf) return logF_high;
    return log_diff_exp(logF_high, logF_low);
  };
  
  if (is_true(any(ok_finite))) {
    lls = d_DDM_Wien(rts, R, pars, ok_finite); // Inf rts will get passed back with R_NegInf ll here
  } else {
    lls.fill(R_NegInf);
  }

  // Apply truncation normaliser to finite trials for standard DDM only.
  // (For DDMGNG we preserve the existing semantics and do not apply truncation.)
  if (!gng && is_true(any(ok_finite))) {
    IntegerVector R1_for_cdf(n_trials, 1);
    IntegerVector R2_for_cdf(n_trials, 2);
    NumericVector logF_LT_1 = p_DDM_Wien(LT, R1_for_cdf, pars, ok_finite);
    NumericVector logF_LT_2 = p_DDM_Wien(LT, R2_for_cdf, pars, ok_finite);
    NumericVector logF_UT_1 = p_DDM_Wien(UT, R1_for_cdf, pars, ok_finite);
    NumericVector logF_UT_2 = p_DDM_Wien(UT, R2_for_cdf, pars, ok_finite);
    for (int i = 0; i < n_trials; ++i) {
      if (!ok_finite[i]) continue;
      if (LT[i] == 0.0 && !R_FINITE(UT[i])) continue; // no truncation
      const double logF_LT_all = log_sum_exp(logF_LT_1[i], logF_LT_2[i]);
      const double logF_UT_all = R_FINITE(UT[i]) ? log_sum_exp(logF_UT_1[i], logF_UT_2[i]) : 0.0; // log(1)
      const double logZ = log_interval_mass(logF_LT_all, logF_UT_all);
      if (R_FINITE(logZ)) lls[i] -= logZ;
      else lls[i] = min_ll;
    }
  }

  if (is_true(any(ok_nonfinite))) {
    if (gng) {
      // Preserve current right-censoring semantics for go/no-go, and add LC-aware
      // handling for rt == -Inf / NA without introducing truncation correction.
      int go_idx = -1; // 1-based index into the R factor levels
      SEXP lev_sexp = R.attr("levels");
      if (lev_sexp != R_NilValue) {
        CharacterVector levs(lev_sexp);
        int nogo_idx = -1;
        for (int j = 0; j < levs.size(); ++j) {
          if (Rcpp::as<std::string>(levs[j]) == "nogo") { nogo_idx = j + 1; break; }
        }
        if (nogo_idx != -1) {
          for (int j = 0; j < levs.size(); ++j) {
            const int idx = j + 1;
            if (idx != nogo_idx) { go_idx = idx; break; }
          }
        }
      }
      if (go_idx == -1) {
        go_idx = 1;
        for (int i = 0; i < n_trials; ++i) {
          if (finite[i] && R[i] != NA_INTEGER) {
            go_idx = R[i];
            break;
          }
        }
      }

      IntegerVector R_go_for_cdf(n_trials, go_idx);
      NumericVector logcdf_U_go = p_DDM_Wien(UC, R_go_for_cdf, pars, ok_nonfinite);
      NumericVector logcdf_L_go = p_DDM_Wien(LC, R_go_for_cdf, pars, ok_nonfinite);

      for (int i = 0; i < n_trials; ++i) {
        if (!ok_nonfinite[i]) continue;
        if (rts[i] == R_PosInf) {
          // Keep existing right-censoring behavior: no go-boundary hit before UC.
          lls[i] = log1m_exp(logcdf_U_go[i]);
        } else if (rts[i] == R_NegInf) {
          // LC support: go-boundary hit by LC.
          lls[i] = logcdf_L_go[i];
        } else {
          // NA non-response/censored coding can represent left or right censoring.
          lls[i] = log_sum_exp(logcdf_L_go[i], log1m_exp(logcdf_U_go[i]));
        }
      }
    } else {
      // Standard DDM: full LC/UC/LT/UT interval handling for non-finite RTs.
      IntegerVector R1_for_cdf(n_trials, 1);
      IntegerVector R2_for_cdf(n_trials, 2);
      NumericVector logF_LT_1 = p_DDM_Wien(LT, R1_for_cdf, pars, ok_nonfinite);
      NumericVector logF_LT_2 = p_DDM_Wien(LT, R2_for_cdf, pars, ok_nonfinite);
      NumericVector logF_LC_1 = p_DDM_Wien(LC, R1_for_cdf, pars, ok_nonfinite);
      NumericVector logF_LC_2 = p_DDM_Wien(LC, R2_for_cdf, pars, ok_nonfinite);
      NumericVector logF_UC_1 = p_DDM_Wien(UC, R1_for_cdf, pars, ok_nonfinite);
      NumericVector logF_UC_2 = p_DDM_Wien(UC, R2_for_cdf, pars, ok_nonfinite);
      NumericVector logF_UT_1 = p_DDM_Wien(UT, R1_for_cdf, pars, ok_nonfinite);
      NumericVector logF_UT_2 = p_DDM_Wien(UT, R2_for_cdf, pars, ok_nonfinite);

      for (int i = 0; i < n_trials; ++i) {
        if (!ok_nonfinite[i]) continue;

        const bool r_known = (R[i] != NA_INTEGER);
        const int r_idx = r_known ? R[i] : 0;

        const double logF_LT_all = log_sum_exp(logF_LT_1[i], logF_LT_2[i]);
        const double logF_LC_all = log_sum_exp(logF_LC_1[i], logF_LC_2[i]);
        const double logF_UC_all = log_sum_exp(logF_UC_1[i], logF_UC_2[i]);
        const double logF_UT_all = R_FINITE(UT[i]) ? log_sum_exp(logF_UT_1[i], logF_UT_2[i]) : 0.0; // log(1)

        if (rts[i] == R_NegInf) {
          // Left-censored interval (LT, LC]
          if (!r_known) {
            lls[i] = log_interval_mass(logF_LT_all, logF_LC_all);
          } else if (r_idx == 1) {
            lls[i] = log_interval_mass(logF_LT_1[i], logF_LC_1[i]);
          } else {
            lls[i] = log_interval_mass(logF_LT_2[i], logF_LC_2[i]);
          }
        } else if (rts[i] == R_PosInf) {
          // Right-censored interval (UC, UT]
          if (!r_known) {
            lls[i] = log_interval_mass(logF_UC_all, logF_UT_all);
          } else if (r_idx == 1) {
            lls[i] = log_interval_mass(logF_UC_1[i], logF_UT_1[i]);
          } else {
            lls[i] = log_interval_mass(logF_UC_2[i], logF_UT_2[i]);
          }
        } else {
          // rt == NA: union of left and right censored intervals
          if (!r_known) {
            const double ll_left = log_interval_mass(logF_LT_all, logF_LC_all);
            const double ll_right = log_interval_mass(logF_UC_all, logF_UT_all);
            lls[i] = log_sum_exp(ll_left, ll_right);
          } else if (r_idx == 1) {
            const double ll_left = log_interval_mass(logF_LT_1[i], logF_LC_1[i]);
            const double ll_right = log_interval_mass(logF_UC_1[i], logF_UT_1[i]);
            lls[i] = log_sum_exp(ll_left, ll_right);
          } else {
            const double ll_left = log_interval_mass(logF_LT_2[i], logF_LC_2[i]);
            const double ll_right = log_interval_mass(logF_UC_2[i], logF_UT_2[i]);
            lls[i] = log_sum_exp(ll_left, ll_right);
          }
        }
      }
    }
  }
  lls_exp = c_expand(lls, expand); // decompress
  // lls_exp = lls;
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return(sum(lls_exp));
}

// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, List bounds, List transforms, List pretransforms,
                      CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue){
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  NumericMatrix pars(n_trials, p_types.length());
  p_vector.names() = p_names;
  LogicalVector is_ok(n_trials);
  
  // Once (outside the main loop over particles):
  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<PreTransformSpec> p_specs;
  std::vector<BoundSpec> bound_specs;
  std::vector<TransformSpec> full_t_specs; // precomputed transform specs for p_types
  Rcpp::List trend_list;
  if (!trend.isNull()) {
    trend_list = Rcpp::List(trend);
  } else {
    trend_list = Rcpp::List::create();
  }
  std::string type_std = Rcpp::as<std::string>(Rcpp::wrap(type));
  if(type_std.find("DDM") != std::string::npos){
    bool gng = (type_std.find("GNG") != std::string::npos);
    IntegerVector expand = data.attr("expand");
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      p_vector.names() = p_names;
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
      // Precompute specs
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok, gng);
    }
  } else if(type == "MRI" || type == "MRI_AR1"){
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      p_vector.names() = p_names;
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
      // Precompute specs
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      if(type == "MRI"){
        lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      } else{
        lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
    }
  } else if(type_std.find("SS") != std::string::npos){
      IntegerVector expand = data.attr("expand");
      NumericVector lR = data["lR"];
      int n_lR = unique(lR).length();
      int n_trials_ss = (n_lR > 0) ? (n_trials / n_lR) : n_trials;
      // Pick function pointers and indices based on type
      bool is_exg = type_std.find("EXG") != std::string::npos;
      ss_go_pdf_fn go_lpdf_ptr = is_exg ? texg_go_lpdf : rdex_go_lpdf;
      ss_go_pdf_fn go_lccdf_ptr = is_exg ? texg_go_lccdf : rdex_go_lccdf;
      ss_stop_surv_fn stop_logsurv_ptr = is_exg ? stop_logsurv_texg_fn : stop_logsurv_rdex_fn;
      ss_stop_success_fn stop_success_ptr = is_exg ? ss_texg_stop_success_lpdf : ss_rdex_stop_success_lpdf;
      int idx_tf = is_exg ? 6 : 8;
      int idx_gf = is_exg ? 7 : 9;
      for (int i = 0; i < n_particles; ++i) {
        p_vector = p_matrix(i, _);
        p_vector.names() = p_names;
        if(i == 0){
          p_specs = make_pretransform_specs(p_vector, pretransforms);
          // Precompute transform specs for all p_types using a one-time dummy
          NumericMatrix dummy(1, p_types.size());
          colnames(dummy) = p_types;
          full_t_specs = make_transform_specs(dummy, transforms);
        }
        pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
        if (i == 0) {                            // first particle only, just to get colnames
          bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
        }
        is_ok = c_do_bound(pars, bound_specs);
        is_ok = lr_all(is_ok, n_lR); // reduce to per-trial ok
        lls[i] = c_log_likelihood_ss(pars, data, n_trials_ss, expand, min_ll, is_ok,
                                     go_lpdf_ptr, go_lccdf_ptr,
                                     stop_logsurv_ptr, stop_success_ptr,
                                     idx_tf, idx_gf);
      }
    } else {
    // ZH I have adjusted a lot here. I made new race model pointers that help branch things like intrinsic omissions LBA, and the new GNG models
    // The "1" versions handle scalar likelihoods that do not call Rcpp functions,which speeds up the GSL integrations used in censoring, truncation, and GNG
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Standardize incoming model type string from R
    // For race models: e.g. "LBA", "LBAIO", "LBAGNG", "RDM", "RDMGNG", "LNR", "LNRGNG"

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);
    
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    
    bool all_finite_trials = false;
    if (data.hasAttribute("emc2_all_finite_trials")) {
      Rcpp::LogicalVector v = data.attr("emc2_all_finite_trials");
      if (v.size() == 1 && v[0] != NA_LOGICAL) {
        all_finite_trials = v[0];
      }
    } else if (n_trials > 0 && n_lR > 0 && (n_trials % n_lR) == 0) {
      // Fallback for direct calc_ll callers that bypass R-side cache setup.
      all_finite_trials = true;
      NumericVector rts = data["rt"];
      IntegerVector R_idx = data["R"];
      for (int start_row_idx = 0; start_row_idx < n_trials; start_row_idx += n_lR) {
        const double rt_j = rts[start_row_idx];
        const int R_j_idx = R_idx[start_row_idx];
        if (!(R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx != NA_INTEGER)) {
          all_finite_trials = false;
          break;
        }
      }
    }
      
    for (int i = 0; i < n_particles; ++i) {
        p_vector = p_matrix(i, Rcpp::_);
        p_vector.names() = p_names;
        if (i == 0) {
          p_specs = make_pretransform_specs(p_vector, pretransforms);
          // Precompute transform specs for all p_types using a one-time dummy
          NumericMatrix dummy(1, p_types.size());
          colnames(dummy) = p_types;
          full_t_specs = make_transform_specs(dummy, transforms);
        }
        pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
        if (i == 0) {                            // first particle only, just to get colnames
          bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
        }
        is_ok = c_do_bound(pars, bound_specs);
        is_ok = lr_all(is_ok, n_lR);
        lls[i] = c_log_likelihood_race(pars, data,
                                       adapter.model_dfun_ptr, adapter.model_pfun_ptr,
                                       adapter.pdf1_ptr, adapter.cdf1_ptr,
                                       n_trials,
                                       winner, expand, min_ll, is_ok, n_lR,
                                       &adapter.ctx,
                                       all_finite_trials);
      }
  }
  return(lls);
}

// [[Rcpp::export]]
NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue) {
  const int n_particles = particle_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  LogicalVector is_ok(n_trials);
  NumericMatrix pars;

  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs;
  CharacterVector p_names = colnames(particle_matrix);
  std::string type_std = Rcpp::as<std::string>(Rcpp::wrap(type));

  NumericMatrix one_particle(1, particle_matrix.ncol());
  colnames(one_particle) = p_names;
  IntegerVector kernel_output_codes = IntegerVector::create(1);
  auto pars_for_particle = [&](int i) -> NumericMatrix {
    for (int j = 0; j < particle_matrix.ncol(); ++j) {
      one_particle(0, j) = particle_matrix(i, j);
    }
    return get_pars_c_wrapper_oo_core(one_particle, data, constants, designs, bounds, transforms,
                                      pretransforms, trend, false, false, kernel_output_codes);
  };

  if (type_std.find("DDM") != std::string::npos) {
    bool gng = (type_std.find("GNG") != std::string::npos);
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok, gng);
    }
  } else if (type == "MRI" || type == "MRI_AR1") {
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      if (type == "MRI") {
        lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      } else {
        lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
    }
  } else if (type_std.find("SS") != std::string::npos) {
    IntegerVector expand = data.attr("expand");
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    int n_trials_ss = (n_lR > 0) ? (n_trials / n_lR) : n_trials;
    bool is_exg = type_std.find("EXG") != std::string::npos;
    ss_go_pdf_fn go_lpdf_ptr = is_exg ? texg_go_lpdf : rdex_go_lpdf;
    ss_go_pdf_fn go_lccdf_ptr = is_exg ? texg_go_lccdf : rdex_go_lccdf;
    ss_stop_surv_fn stop_logsurv_ptr = is_exg ? stop_logsurv_texg_fn : stop_logsurv_rdex_fn;
    ss_stop_success_fn stop_success_ptr = is_exg ? ss_texg_stop_success_lpdf : ss_rdex_stop_success_lpdf;
    int idx_tf = is_exg ? 6 : 8;
    int idx_gf = is_exg ? 7 : 9;

    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      lls[i] = c_log_likelihood_ss(pars, data, n_trials_ss, expand, min_ll, is_ok,
                                   go_lpdf_ptr, go_lccdf_ptr,
                                   stop_logsurv_ptr, stop_success_ptr,
                                   idx_tf, idx_gf);
    }
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll_oo");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);

    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    bool all_finite_trials = false;
    if (data.hasAttribute("emc2_all_finite_trials")) {
      Rcpp::LogicalVector v = data.attr("emc2_all_finite_trials");
      if (v.size() == 1 && v[0] != NA_LOGICAL) {
        all_finite_trials = v[0];
      }
    } else if (n_trials > 0 && n_lR > 0 && (n_trials % n_lR) == 0) {
      // Fallback for direct calc_ll_oo callers that bypass R-side cache setup.
      all_finite_trials = true;
      NumericVector rts = data["rt"];
      IntegerVector R_idx = data["R"];
      for (int start_row_idx = 0; start_row_idx < n_trials; start_row_idx += n_lR) {
        const double rt_j = rts[start_row_idx];
        const int R_j_idx = R_idx[start_row_idx];
        if (!(R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx != NA_INTEGER)) {
          all_finite_trials = false;
          break;
        }
      }
    }

    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      lls[i] = c_log_likelihood_race(pars, data,
                                     adapter.model_dfun_ptr, adapter.model_pfun_ptr,
                                     adapter.pdf1_ptr, adapter.cdf1_ptr,
                                     n_trials,
                                     winner, expand, min_ll, is_ok, n_lR,
                                     &adapter.ctx,
                                     all_finite_trials);
    }
  }
  return lls;
}



// [[Rcpp::export]]
NumericMatrix get_pars_c_wrapper(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
                                 List designs, List bounds, List transforms, List pretransforms,
                                 CharacterVector p_types,
                                 Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                 bool return_kernel_matrix = false,
                                 bool drop_trend_pars = true){
  if (return_kernel_matrix && trend.isNull()) {
    stop("return_kernel_matrix=TRUE requires a non-NULL 'trend' specification");
  }
  if (return_kernel_matrix) {
    return get_pars_c_wrapper_oo_core(p_matrix, data, constants, designs, bounds, transforms,
                                      pretransforms, trend, true, false, IntegerVector::create(1));
  }

  const int n_trials = data.nrow();
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  p_vector.names() = p_names;
  NumericMatrix pars;

  std::vector<PreTransformSpec> p_specs;
  std::vector<TransformSpec> full_t_specs;

  p_vector = p_matrix(0, _);
  p_specs = make_pretransform_specs(p_vector, pretransforms);
  NumericMatrix dummy(1, p_types.size());
  colnames(dummy) = p_types;
  full_t_specs = make_transform_specs(dummy, transforms);

  Rcpp::List trend_list;
  if (!trend.isNull()) {
    trend_list = Rcpp::List(trend);
  } else {
    trend_list = Rcpp::List::create();
  }

  pars = get_pars_matrix(p_vector, constants,
                         p_specs, p_types,
                         designs, n_trials, data,
                         trend_list, full_t_specs, drop_trend_pars);
  return pars;
}

// [[Rcpp::export]]
NumericMatrix get_pars_c_wrapper_oo(NumericMatrix particle_matrix,
                                    DataFrame data,
                                    NumericVector constants,
                                    List designs,
                                    List bounds,
                                    List transforms,
                                    List pretransforms,
                                    Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                    bool return_kernel_matrix = false,
                                    bool return_all_pars = false,
                                    IntegerVector kernel_output_codes = 1) {
  return get_pars_c_wrapper_oo_core(particle_matrix, data, constants, designs, bounds, transforms,
                                    pretransforms, trend, return_kernel_matrix,
                                    return_all_pars, kernel_output_codes);
}

// gsl adapter for integrals - uses scalar, Rcpp-independent functions for speed
double gsl_f_race_scalar(double t, void* p) {
  auto* P = static_cast<gsl_race_params_scalar*>(p);
  if (t <= 0.0) return 0.0;
  const int w = P->winner_idx0;
  if (w < 0 || w >= P->n_lR) return 0.0;
  if (!P->isok[w]) return 0.0;
  
  const double* par_w = P->pars + static_cast<size_t>(w) * P->n_par;
  double out = P->pdf1(t, par_w, P->ctx);
  if (!(out > 0.0) || !std::isfinite(out)) return 0.0;
  
  for (int j = 0; j < P->n_lR; ++j) {
    if (j == w) continue;
    if (!P->isok[j]) return 0.0;
    const double* par_j = P->pars + static_cast<size_t>(j) * P->n_par;
    const double cdf = P->cdf1(t, par_j, P->ctx);
    if (!(cdf > 0.0)) continue;
    if (cdf >= 1.0 || !std::isfinite(cdf)) return 0.0;
    out *= (1.0 - cdf);
    if (!(out > 0.0) || !std::isfinite(out)) return 0.0;
  }
  return out;
}

// Log probability of intrinsic omission
inline double log_pIO_rowmajor(const double* pars_rowmajor,
                                                 const int* isok_int,
                                                 int n_lR,
                                                 int n_par) {
  double log_p = 0.0;
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    const double v = par_k[0];
    const double sv = par_k[1];
    if (!std::isfinite(v) || !std::isfinite(sv) || sv <= 0.0) return R_NegInf;
    const double ll = R::pnorm(0.0, v, sv, 1, 1);
    if (!std::isfinite(ll)) return R_NegInf;
    log_p += ll;
  }
  return log_p;
}

// Log survivor and cdf of the race at time t:
//   log S(t) = sum_k log(1 - F_k(t))
//
// "rowmajor" here means per-trial parameters for the n_lR accumulators are
// packed as a contiguous buffer:
//   pars_rowmajor[k * n_par + c]
//
// This avoids repeatedly indexing into an Rcpp::NumericMatrix inside tight
// loops and matches the representation needed for the scalar/GSL integrands.

inline double log_survivor_rowmajor(double t,
                                    const double* pars_rowmajor,
                                    const int* isok_int,
                                    int n_lR,
                                    int n_par,
                                    RaceCdf1Fun cdf1,
                                    void* ctx) {
  if (t == R_PosInf) {
    auto* race_ctx = static_cast<ContextForRaceModels*>(ctx);
    if (race_ctx && !race_ctx->use_posdrift) {
      return log_pIO_rowmajor(pars_rowmajor, isok_int, n_lR, n_par);
    }
    return R_NegInf;
  }
  double logS = 0.0;
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = safe_log1m_race(cdf);
    if (!std::isfinite(ll)) return R_NegInf;
    logS += ll;
  }
  return logS;
}

inline double log_cdf_rowmajor(double t,
                               const double* pars_rowmajor,
                               const int* isok_int,
                               int n_lR,
                               int n_par,
                               RaceCdf1Fun cdf1,
                               void* ctx) {
  if (t == R_PosInf) {
      auto* race_ctx = static_cast<ContextForRaceModels*>(ctx);
      if (race_ctx && !race_ctx->use_posdrift) {
          double logC = 0.0;
          for (int k = 0; k < n_lR; ++k) {
            if (!isok_int[k]) return R_NegInf;
            const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
            const double v = par_k[0];
            const double sv = par_k[1];
            if (!std::isfinite(v) || !std::isfinite(sv) || sv <= 0.0) return R_NegInf;
            const double lp_neg = R::pnorm(0.0, v, sv, 1, 1);
            const double ll = log1m_exp(lp_neg);
            if (!std::isfinite(ll)) return R_NegInf;
            logC += ll;
          }
          return logC;
     }
    return 0.0;
  }
  double logC = 0.0;
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = std::log(cdf);
    if (!std::isfinite(ll)) return R_NegInf;
    logC += ll;
  }
  return logC;
}

inline double log_min_density_rowmajor(double t,
                                       const double* pars_rowmajor,
                                       const int* isok_int,
                                       int n_lR,
                                       int n_par,
                                       RacePdf1Fun pdf1,
                                       RaceCdf1Fun cdf1,
                                       void* ctx,
                                       std::vector<double>& logS_k) {
  // Log density of the minimum (no known winner) at time t:
  //   f_min(t) = sum_k f_k(t) * prod_{j != k} (1 - F_j(t))
  //
  // Work in log space:
  //   log f_min(t) = logsumexp_k [ log f_k(t) + sum_{j != k} log(1 - F_j(t)) ].
  //
  // logS_k is a per-accumulator scratch buffer (reused across calls) holding
  // log(1 - F_k(t)) for the current t, to avoid reallocations in the "other
  // trial" loop.
  if (!(t > 0.0) || !std::isfinite(t)) return R_NegInf;
  double logS_all = 0.0;
  if (logS_k.size() < static_cast<size_t>(n_lR)) {
    logS_k.assign(static_cast<size_t>(n_lR), R_NegInf);
  } else {
    std::fill(logS_k.begin(), logS_k.begin() + n_lR, R_NegInf);
  }
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = safe_log1m_race(cdf);
    if (!std::isfinite(ll)) return R_NegInf;
    logS_k[static_cast<size_t>(k)] = ll;
    logS_all += ll;
  }
  double out = R_NegInf;
  for (int k = 0; k < n_lR; ++k) {
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    const double pdf = pdf1(t, par_k, ctx);
    if (!(pdf > 0.0) || !std::isfinite(pdf)) continue;
    const double term = std::log(pdf) + (logS_all - logS_k[static_cast<size_t>(k)]);
    out = log_sum_exp(out, term);
  }
  return out;
}


double integrate_for_kth_winner_rowmajor_cpp(
    int k_winner_idx, // 1-based
    const double* pars_rowmajor,
    const int* isok_int,
    double low,
    double upp,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    int n_lR_j,
    int n_par,
    const GslIntegrationControls& gsl_ctl,
    void* model_specific_context,
    gsl_integration_workspace* w) {
  
  // Integrate the k-th winner density over an interval [low, upp] for a single
  // unique trial, using rowmajor buffers (raw pointers). This is the fast path
  // used by truncation/censoring normalisers and by go/no-go branches; it avoids
  // Rcpp object traffic inside the GSL callback.
  if (low >= upp && !(low == 0 && upp == R_PosInf)) return R_NegInf;
  if (k_winner_idx < 1 || k_winner_idx > n_lR_j) return R_NegInf;
  if (w == nullptr) Rcpp::stop("integrate_for_kth_winner_rowmajor_cpp: GSL workspace is null.");
  
  gsl_function F;
  gsl_race_params_scalar params_struct;
  params_struct.pars = pars_rowmajor;
  params_struct.n_lR = n_lR_j;
  params_struct.n_par = n_par;
  params_struct.winner_idx0 = k_winner_idx - 1;
  params_struct.isok = isok_int;
  params_struct.pdf1 = pdf1;
  params_struct.cdf1 = cdf1;
  params_struct.ctx = model_specific_context;
  
  F.function = &gsl_f_race_scalar;
  F.params = &params_struct;
  
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status;
  double result = 0.0;
  double error = 0.0;
  auto run_integral = [&](double abs_tol, double rel_tol, size_t limit) -> int {
    if (upp == R_PosInf) {
      if (low < 0) low = 0; // QAGIU requires a >= 0
      if (low >= R_PosInf) {
        result = 0.0;
        error = 0.0;
        return GSL_SUCCESS;
      }
      return gsl_integration_qagiu(&F, low, abs_tol, rel_tol, limit, w, &result, &error);
    }
    return gsl_integration_qags(&F, low, upp, abs_tol, rel_tol, limit, w, &result, &error);
  };
  status = run_integral(gsl_ctl.abs_tol, gsl_ctl.rel_tol, gsl_ctl.limit);
  // Retry once with tighter controls if the fast pass fails.
  if (status != GSL_SUCCESS) {
    status = run_integral(gsl_ctl.retry_abs_tol, gsl_ctl.retry_rel_tol, gsl_ctl.retry_limit);
  }
  
  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS) return R_NegInf;
  if (!(result > 0.0) || !R_finite(result)) return R_NegInf;
  return std::log(result);
}

double get_trunc_normaliser_rowmajor_cpp(const double* pars_rowmajor,
                                         const int* isok_int,
                                         RacePdf1Fun pdf1,
                                         RaceCdf1Fun cdf1,
                                         double LT,
                                         double UT,
                                         int n_lR,
                                         int n_par,
                                         const GslIntegrationControls& gsl_ctl,
                                         void* model_specific_context,
                                         GslWorkspacePtr& workspace) {
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
  const double logS_LT = log_survivor_rowmajor(LT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
  if (!R_FINITE(logS_LT)) return R_NegInf;
  if (UT == R_PosInf) return logS_LT;
  
  const double logS_UT = log_survivor_rowmajor(UT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
  double logP = log_diff_exp(logS_LT, logS_UT);
  if (R_FINITE(logP) && logP > log_prob_eps) return logP;
  
  // Only falls back to GSL integration if the analytic trick fails (e.g. due to catastrophic cancellation)
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
  double log_total = R_NegInf;
  for (int k_win = 1; k_win <= n_lR; ++k_win) {
    const double log_k = integrate_for_kth_winner_rowmajor_cpp(k_win,
                                                               pars_rowmajor,
                                                               isok_int,
                                                               LT,
                                                               UT,
                                                               pdf1,
                                                               cdf1,
                                                               n_lR,
                                                               n_par,
                                                               gsl_ctl,
                                                               model_specific_context,
                                                               w);
    log_total = log_sum_exp(log_total, log_k);
  }
  if (R_FINITE(log_total) && log_total > log_prob_eps) return log_total;
  return R_NegInf;
}

// Main C++ function for censored/truncated race likelihood calculation
// This function is now the unified entry point for all race models (LBA, RDM, LNR),
// whether they are standard or explicitly handling censoring/truncation.
// It uses batching for finite RTs and iterative processing for others (censored/NA RTs).
double c_log_likelihood_race(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, covering all dadm rows for that particle
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions, structured for all accumulators
    RacePdfFun model_dfun,                  // Pointer to the model's PDF adapter function
    RaceCdfFun model_pfun,                  // Pointer to the model's CDF adapter function
    RacePdf1Fun pdf1,                       // Scalar PDF (for GSL)
    RaceCdf1Fun cdf1,                       // Scalar CDF (for GSL)
    const int n_trials,
    LogicalVector winner,
    Rcpp::IntegerVector expand,  // Vector for expanding unique LLs to full trial count
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector isok,   // Parameter validity for each row in 'pars' matrix
    int n_lR,                              // Number of accumulators in the race (must be > 0 if data exists)
    void* model_context_for_funcs,          // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
    bool all_finite_trials             // Data-only hint: all trials finite/in-bounds/known response
) {
  
  // Only allocate a GSL workspace if/when numerical integration is needed.
  GslWorkspacePtr workspace(nullptr, &gsl_integration_workspace_free);
  
  // Fetch censoring and truncation values from dadm columns. These are passed across all rows for ease of access, but should be identical at least at the subject level as they don't correspond to a data entry. Attributes probably a better fit here but clunky.
  Rcpp::NumericVector LT = get_col_with_default(dadm, "LT", 0.0);
  Rcpp::NumericVector UT = get_col_with_default(dadm, "UT", R_PosInf);
  Rcpp::NumericVector LC = get_col_with_default(dadm, "LC", 0.0);
  Rcpp::NumericVector UC = get_col_with_default(dadm, "UC", R_PosInf);
  // Fast default controls for MCMC throughput, with an automatic stricter retry.
  // Retry settings match the old behaviour (rel_tol=1e-7, limit=1000).
  const GslIntegrationControls gsl_ctl{
    1e-10, // abs_tol
    1e-5,  // rel_tol
    200,   // max subintervals
    0.0,   // retry_abs_tol
    1e-7,  // retry_rel_tol
    1000   // retry max subintervals
  };
  int n_lR_j = n_lR;
  Rcpp::NumericVector rts_dadm = dadm["rt"];
  Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
  NumericVector lds(n_trials, min_ll); // initialise at min_ll
  // If a RACE column exists, set parameters of accumulators not present on a
  // given trial to NA so the density functions return zero for them. This
  // mirrors logic from the old c_log_likelihood_race implementation.
  bool has_RACE_col = dadm.containsElementNamed("RACE");
  Rcpp::IntegerVector RACE;
  Rcpp::LogicalVector RACE_mask;
  if (has_RACE_col) {
    if (dadm.hasAttribute("RACE_nacc_by_row") && dadm.hasAttribute("RACE_mask")) {
      RACE = dadm.attr("RACE_nacc_by_row");
      RACE_mask = dadm.attr("RACE_mask");
      for (int row = 0; row < pars.nrow(); ++row) {
        if (!RACE_mask[row]) std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
      }
    } else {
      has_RACE_col = false;
    }
  }
  if (n_trials == 0) return 0.0; // No data, no likelihood
  
  if (n_lR <= 0) Rcpp::stop("c_log_likelihood_race: n_lR must be positive and correctly determined before this call.");
  if (n_trials % n_lR != 0) Rcpp::stop("c_log_likelihood_race: dadm nrows not a multiple of n_lR.");
  
  // Here we check for a pC parameter corresponding to probability of contaminant OMISSION.
  // Index the column (if it exists)
  bool use_pC = false;
  int pc_col = -1;
  Rcpp::List dimnames = pars.attr("dimnames");
  Rcpp::CharacterVector colnames = as<Rcpp::CharacterVector>(dimnames[1]);
  for (int j = 0; j < colnames.size(); ++j) {
    if (as<std::string>(colnames[j]) == "pContaminant") {
      pc_col = j;
      use_pC = true;
      break;
    }
  }
  // Also check that pC is not all zeroes (i.e. a model that has pC optional but is not using it here)
  if (use_pC) {
    bool all_zero = true;
    for (int i = 0; i < pars.nrow(); ++i) {
      if (pars(i, pc_col) != 0.0) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      use_pC = false;
    }
  }
  
  int n_unique_trials = n_trials / n_lR;
  Rcpp::NumericVector ll_unique(n_unique_trials, min_ll);
  Rcpp::NumericVector pC_values(n_unique_trials);
  if (use_pC) {
    for (int j = 0; j < n_unique_trials; ++j) {
      pC_values[j] = pars(j * n_lR, pc_col);
    }
  }
  
  // Parameter matrix and validity vector checks
  if (pars.nrow() != n_trials) {
    Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_trials: " << n_trials << std::endl;
    Rcpp::stop("c_log_likelihood_race: pars matrix dimensions do not match total dadm rows.");
  }
  if (isok.size() != pars.nrow()) {
    Rcpp::stop("c_log_likelihood_race: isok size does not match "
                 "pars matrix rows.");
  }
  if (winner.size() != pars.nrow()) {
    Rcpp::stop("c_log_likelihood_race: isok size does not match pars matrix rows.");
  }
  const int n_par = pars.ncol();
  std::vector<double> pars_rowmajor_buffer(static_cast<size_t>(n_lR) * n_par);
  std::vector<int> isok_int_buffer(static_cast<size_t>(n_lR), 0);
  std::vector<double> logS_k_buffer(static_cast<size_t>(n_lR), R_NegInf);
  auto fill_trial_buffers = [&](int start_row_idx, int n_lR_j) {
    for (int k = 0; k < n_lR_j; ++k) {
      isok_int_buffer[k] = isok[start_row_idx + k] ? 1 : 0;
      for (int c = 0; c < n_par; ++c) {
        pars_rowmajor_buffer[static_cast<size_t>(k) * n_par + c] = pars(start_row_idx + k, c);
      }
    }
  };
  
  // Fast path hint: computed once per calc_ll call (data-only), to avoid re-scanning per particle.
  const bool use_full_finite_batch = all_finite_trials;
  Rcpp::LogicalVector finite_rt_mask;
  std::vector<int> finite_rt_unique_trial_indices;
  std::vector<int> other_unique_trial_indices;
  if (!use_full_finite_batch) {
    const bool has_partition_attrs =
      dadm.hasAttribute("finite_rt_mask") &&
      dadm.hasAttribute("finite_rt_unique_trial_indices") &&
      dadm.hasAttribute("other_unique_trial_indices");

    if (has_partition_attrs) {
      finite_rt_mask = dadm.attr("finite_rt_mask");
      Rcpp::IntegerVector finite_attr = dadm.attr("finite_rt_unique_trial_indices");
      Rcpp::IntegerVector other_attr = dadm.attr("other_unique_trial_indices");
      finite_rt_unique_trial_indices.assign(finite_attr.begin(), finite_attr.end());
      other_unique_trial_indices.assign(other_attr.begin(), other_attr.end());
    } else {
      // Fallback path for direct calc_ll/calc_ll_oo callers that bypass
      // .cache_ll_data_attrs(): derive finite/other unique-trial partitions.
      finite_rt_mask = Rcpp::LogicalVector(n_trials, false);
      for (int unique_trial_idx = 0; unique_trial_idx < n_unique_trials; ++unique_trial_idx) {
        const int start_row_idx = unique_trial_idx * n_lR;
        const int n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
        const double rt_j = rts_dadm[start_row_idx];
        const int R_j_idx = R_idxs_dadm[start_row_idx];
        if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx != NA_INTEGER) {
          finite_rt_unique_trial_indices.push_back(unique_trial_idx);
          for (int k = 0; k < n_lR_j; ++k) finite_rt_mask[start_row_idx + k] = true;
        } else {
          other_unique_trial_indices.push_back(unique_trial_idx);
        }
      }
    }
  }
  double log_Z_this = 0;  // Default inv_Z if no truncation. Should never be used but here as a precaution.
  // Scratch interval bounds for per-trial branches below. Keeping these as
  // mutable locals makes the censoring/truncation logic read like "set bounds,
  // then integrate / take survivor", rather than threading LTj/UTj/LCj/UCj
  // through every helper call.
  double lower_for_trial = 0;
  double upper_for_trial = R_PosInf;
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  ContextForRaceModels dense_ctx = *ctx;
  dense_ctx.min_lik_for_pdf = 0.0;
  void* dense_ctx_ptr = static_cast<void*>(&dense_ctx);
  bool gng=ctx->gng;
  bool posdrift=ctx->use_posdrift;
  
  const bool has_finite_batch = use_full_finite_batch || (finite_rt_unique_trial_indices.size() > 0);
  if (has_finite_batch) {
    Rcpp::LogicalVector idx_win;
    Rcpp::LogicalVector idx_loss;
    bool any_win = false;
    bool any_loss = false;
    if (use_full_finite_batch && !has_RACE_col) {
      idx_win = winner;
      any_win = true;
      if (n_lR > 1) {
        idx_loss = !winner;
        any_loss = true;
      }
    } else {
      idx_win = Rcpp::LogicalVector(n_trials, false);
      idx_loss = Rcpp::LogicalVector(n_trials, false);
      for (int i = 0; i < n_trials; ++i) {
        if (has_RACE_col && !RACE_mask[i]) continue;
        if (!use_full_finite_batch && !finite_rt_mask[i]) continue;
        if (winner[i]) { idx_win[i] = true; any_win = true; }
        else { idx_loss[i] = true; any_loss = true; }
      }
    }
    
    if (any_win) {
      Rcpp::NumericVector win_pdf = model_dfun(rts_dadm, pars, idx_win, isok, dense_ctx_ptr);
      int out_i = 0;
      for (int i = 0; i < n_trials; ++i) {
        if (!idx_win[i]) continue;
        const double pdf = win_pdf[out_i++];
        lds[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : R_NegInf;
      }
    }
    
    if (n_lR > 1 && any_loss) {
      Rcpp::NumericVector loss_cdf = model_pfun(rts_dadm, pars, idx_loss, isok, dense_ctx_ptr);
      int out_i = 0;
      for (int i = 0; i < n_trials; ++i) {
        if (!idx_loss[i]) continue;
        double cdf = loss_cdf[out_i++];
        cdf = clamp_cdf01_race(cdf);
        lds[i] = safe_log1m_race(cdf);
      }
    }

    // Apply truncation correction and calculate log-likelihood for each trial in the batch
    const int n_finite_unique = use_full_finite_batch
    ? n_unique_trials
    : static_cast<int>(finite_rt_unique_trial_indices.size());
    for (int i = 0; i < n_finite_unique; ++i) {
      // When all unique trials are finite, iterate in order 0..n_unique_trials-1.
      // Otherwise, iterate over the precomputed subset of finite-RT unique trials.
      int unique_trial_idx = use_full_finite_batch
      ? i
      : finite_rt_unique_trial_indices[static_cast<size_t>(i)];
      int start_row_idx = unique_trial_idx * n_lR;
      n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
      log_Z_this = 0.0;
      const double LTj = LT[start_row_idx];
      const double UTj = UT[start_row_idx];
      if (LTj != 0.0 || UTj != R_PosInf) { // truncation active
        fill_trial_buffers(start_row_idx, n_lR_j);
        log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer.data(),
                                                       isok_int_buffer.data(),
                                                       pdf1,
                                                       cdf1,
                                                       LTj,
                                                       UTj,
                                                       n_lR_j,
                                                       n_par,
                                                       gsl_ctl,
                                                       model_context_for_funcs,
                                                       workspace);
      }
      double current_trial_ll_sum = 0.0;
      for(int k=0; k < n_lR_j; ++k) {
        int dadm_row_idx = start_row_idx + k; // Correct index into lds
        current_trial_ll_sum += lds[dadm_row_idx];
      }
      
      if ((LT[start_row_idx] != 0 || UT[start_row_idx] != R_PosInf) && (NumericVector::is_na(log_Z_this) || !R_FINITE(log_Z_this))) {
        // If truncation active but inv_Z is bad, probability is effectively zero. Could mean the probability of observing an untruncated RT is zero, which would be bad.
        ll_unique[unique_trial_idx] = min_ll;
        continue;
      }
      else if (LT[start_row_idx] != 0 || UT[start_row_idx] != R_PosInf) { // Truncation active and inv_Z is good
        current_trial_ll_sum -=  log_Z_this;
      }
      
      ll_unique[unique_trial_idx] = current_trial_ll_sum;
      ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]);
    }
    
  }
  // --- Process other trials (Infinite RTs, NA RTs, or finite RTs outside truncation) ---
  // These trials require individual processing, often involving numerical integration for censored intervals.
  double current_ll_val;
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
  for (size_t i = 0; i < other_unique_trial_indices.size(); ++i) {
    int unique_trial_idx = other_unique_trial_indices[i];
    int start_row_idx = unique_trial_idx*n_lR;
    log_Z_this = 0.0;
    n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
    const double rt_j = rts_dadm[start_row_idx];
    const int R_j_idx = R_idxs_dadm[start_row_idx];
    const double LTj = LT[start_row_idx];
    const double UTj = UT[start_row_idx];
    const double LCj = LC[start_row_idx];
    const double UCj = UC[start_row_idx];
    const bool has_trunc = (LTj != 0.0 || UTj != R_PosInf);
    
    const bool needs_model = (rt_j == R_NegInf) || (rt_j == R_PosInf) || Rcpp::NumericVector::is_na(rt_j) ||
      (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER);
    if (!needs_model) {
      ll_unique[unique_trial_idx] = min_ll;
      continue;
    }
    
    fill_trial_buffers(start_row_idx, n_lR_j);
    
    auto integrate_interval = [&](int k_winner_1based, double low, double upp) -> double {
      gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
      return integrate_for_kth_winner_rowmajor_cpp(k_winner_1based,
                                                   pars_rowmajor_buffer.data(),
                                                   isok_int_buffer.data(),
                                                   low,
                                                   upp,
                                                   pdf1,
                                                   cdf1,
                                                   n_lR_j,
                                                   n_par,
                                                   gsl_ctl,
                                                   model_context_for_funcs,
                                                   w);
    };
    
    auto log_survivor = [&](double t) -> double {
      return log_survivor_rowmajor(t,
                                   pars_rowmajor_buffer.data(),
                                   isok_int_buffer.data(),
                                   n_lR_j,
                                   n_par,
                                   cdf1,
                                   model_context_for_funcs);
    };

    current_ll_val = R_NegInf;
    if (rt_j == R_NegInf) { 
      lower_for_trial = LTj;
      upper_for_trial = LCj;
      if (R_j_idx == NA_INTEGER) {
        if (gng) {
          // For go/no-go, rt == -Inf means an observed response below LC.
          // The no-go accumulator cannot be the winner for this event, so we
          // skip that winner index when summing winner-specific densities.
          int k_nogo = -1; // 1-based
          int n_true = 0;
          for (int k = 0; k < n_lR_j; ++k) {
            if (winner[start_row_idx + k]) { n_true++; k_nogo = k + 1; }
          }
          if (n_true != 1) k_nogo = -1; // fallback: keep normal behavior
          for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
            if (k_win == k_nogo) continue;
            current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
          }
        } else {
          // Left-censoring with unknown winner: probability the *minimum* finishes in (LT, LC].
          // This is S(LT) - S(LC), where S(t)=P(min > t)=prod_k (1 - F_k(t)).
          // Note: this works for both posdrift and intrinsic-omission (posdrift=FALSE) cases;
          // any mass at +Inf cancels in the survivor difference.
          const double logP = log_diff_exp(log_survivor(lower_for_trial), log_survivor(upper_for_trial));
          if (R_FINITE(logP) && logP > log_prob_eps) {
            current_ll_val = logP;
          } else {
            for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
              current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
            }
          }
        }
      } else {
        current_ll_val = integrate_interval(R_j_idx, lower_for_trial, upper_for_trial);
      }
    } else if (rt_j == R_PosInf) {
      if (gng) {
        lower_for_trial = LTj;
        upper_for_trial = UCj;
        int k_nogo = -1; // 1-based
        int n_true = 0;
        for (int k = 0; k < n_lR_j; ++k) {
          if (winner[start_row_idx + k]) { n_true++; k_nogo = k + 1; }
        }
        // Robust fallback: for a single-accumulator go/no-go trial with no
        // observed response, winner may be all FALSE in dadm. In that case the
        // only valid no-go accumulator is index 1.
        if (n_true == 0 && n_lR_j == 1 && R_j_idx == NA_INTEGER) {
          k_nogo = 1;
          n_true = 1;
        }
        if (n_true != 1) Rcpp::stop("No winner identified in go/no-go withheld response");
        const double logA = integrate_interval(k_nogo, lower_for_trial, upper_for_trial); // probability of no-go winning
        const double logB = log_survivor(upper_for_trial); // For LBAIO this includes "all accumulators never finished" 
        current_ll_val = log_sum_exp(logA, logB);
      } else {
        lower_for_trial = UCj;
        upper_for_trial = UTj;
        if (n_lR_j == 1) {
          const double cdf = cdf1(lower_for_trial, pars_rowmajor_buffer.data(), model_context_for_funcs);
          current_ll_val = safe_log1m_race(clamp_cdf01_race(cdf));
        } else if (R_j_idx == NA_INTEGER) {
          const double logP = posdrift
          ? log_diff_exp(log_survivor(lower_for_trial), log_survivor(upper_for_trial))
          : log_survivor(lower_for_trial); // For LBAIO this includes "all accumulators never finished" 
          if (R_FINITE(logP) && logP > log_prob_eps) {
            current_ll_val = logP;
          } else { // numerical integration fallback integrates across each winner case
            for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
              current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
            }
            // For the fallback we need to add in pIO
            if (!posdrift) { // should only ever be false for models with intrinsic omissions
              const double log_p_I = log_pIO_rowmajor(pars_rowmajor_buffer.data(),
                                                                        isok_int_buffer.data(),
                                                                        n_lR_j,
                                                                        n_par);
              current_ll_val = log_sum_exp(current_ll_val, log_p_I);
            }
          }
        } else {
          current_ll_val = integrate_interval(R_j_idx, lower_for_trial, upper_for_trial);
        }
      }
    } else if (Rcpp::NumericVector::is_na(rt_j)) {
      const double lower1 = LTj;
      const double upper1 = LCj;
      const double lower2 = UCj;
      const double upper2 = UTj;
      if (R_j_idx != NA_INTEGER) {
        current_ll_val = log_sum_exp(integrate_interval(R_j_idx, lower1, upper1),
                                     integrate_interval(R_j_idx, lower2, upper2));
      } else {
        const double logP1 = log_diff_exp(log_survivor(lower1), log_survivor(upper1));
        const double logP2 = posdrift ? log_diff_exp(log_survivor(lower2), log_survivor(upper2)) : log_survivor(lower2);
        const double logPsum = log_sum_exp(logP1, logP2);
        if (R_FINITE(logPsum) && logPsum > log_prob_eps) {
          current_ll_val = logPsum;
        } else {
          for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
            const double ll_L_k = integrate_interval(k_win, lower1, upper1);
            const double ll_U_k = integrate_interval(k_win, lower2, upper2);
            current_ll_val = log_sum_exp(current_ll_val, log_sum_exp(ll_L_k, ll_U_k));
          }
          if (!posdrift) {
            const double log_p_I = log_pIO_rowmajor(pars_rowmajor_buffer.data(), isok_int_buffer.data(), n_lR_j, n_par);
            current_ll_val = log_sum_exp(current_ll_val, log_p_I);
          }
        }
      }
    } else if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER) {
      current_ll_val = log_min_density_rowmajor(rt_j,
                                                pars_rowmajor_buffer.data(),
                                                isok_int_buffer.data(),
                                                n_lR_j,
                                                n_par,
                                                pdf1,
                                                cdf1,
                                                model_context_for_funcs,
                                                logS_k_buffer);
    }
    
    if (current_ll_val > min_ll && has_trunc) {
      log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer.data(),
                                                     isok_int_buffer.data(),
                                                     pdf1,
                                                     cdf1,
                                                     LTj,
                                                     UTj,
                                                     n_lR_j,
                                                     n_par,
                                                     gsl_ctl,
                                                     model_context_for_funcs,
                                                     workspace);
      if (!R_FINITE(log_Z_this)) current_ll_val = min_ll;
      else current_ll_val -= log_Z_this;
    }
    
    ll_unique[unique_trial_idx] = std::max(min_ll, current_ll_val);
  }
  
  
  // --- Summation of log-likelihoods for all unique trials ---
  double total_ll = 0;
  if (expand.length() > 0) { // If an expansion vector is provided (e.g. from non-compressed dadm)
    if (use_pC) { // multiply all likelihoods by the appropriate pC adjustment if it's present
      for (int j = 0; j < n_unique_trials; ++j) {
        double pC = pC_values[j];
        double log1m_pC = std::log1p(-pC); // log(1-pC) todo: check if one of the safe helpers is appropriate here
        int start_row_idx = j * n_lR;
        double rt_j = rts_dadm[start_row_idx];
        
        if (rt_j == R_PosInf) {
          double term1 = std::log(pC);
          double term2 = log1m_pC + ll_unique[j]; //pIO + pD
          ll_unique[j] = log_sum_exp(term1, term2); // pC * pIO * pD
        } else {
          ll_unique[j] = log1m_pC + ll_unique[j];    // pRT * 1-pC
        }
        
      }
    }
    for (int i = 0; i < expand.length(); ++i) {
      total_ll += ll_unique[expand[i] - 1];
    }
  } else { // Default: sum ll_unique directly (each unique trial counted once)
    for (int j = 0; j < n_unique_trials; ++j) {
      if (use_pC) {
        double pC = pC_values[j];
        double log1m_pC = std::log1p(-pC);
        int start_row_idx = j * n_lR;
        double rt_j = rts_dadm[start_row_idx];
        
        if (R_FINITE(rt_j)) {
          ll_unique[j] = log1m_pC + ll_unique[j];
        } else {
          double term1 = std::log(pC);
          double term2 = log1m_pC + ll_unique[j];
          ll_unique[j] = log_sum_exp(term1, term2);
        }
      }
      total_ll += ll_unique[j];
    }
  }
  
  return total_ll;
}

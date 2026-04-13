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
#include "ParamTable.h"
#include "TrendEngine.h"
#include "transform_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include <cmath>
#include <string>
#include <memory>
#include <array>

using namespace Rcpp;

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

static NumericMatrix get_pars_matrix_oo_fast(ParamTable& param_table,
                                             const Rcpp::List& designs,
                                             TrendRuntime* trend_runtime,
                                             const std::vector<TransformSpec>& full_specs,
                                             const Rcpp::CharacterVector& keep_names) {
  if (trend_runtime) trend_runtime->reset_all_kernels();

  const int n_designs = designs.size();
  LogicalVector map_next(n_designs, false);
  std::unordered_set<std::string> transform_next;
  std::unordered_set<std::string> empty_set;

  const auto& premap_set = trend_runtime ? trend_runtime->premap_trend_params() : empty_set;
  const auto& pretransform_set = trend_runtime ? trend_runtime->pretransform_trend_params() : empty_set;

  if (trend_runtime && trend_runtime->has_premap()) {
    map_next = trend_runtime->premap_design_mask(designs);
    param_table.map_from_designs(designs, map_next);

    const auto& specs_premap = trend_runtime->premap_specs();
    if (!specs_premap.empty()) {
      c_do_transform_pt(param_table, specs_premap);
    }

    for (std::size_t i = 0; i < trend_runtime->premap_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->premap_ops[i], param_table);
    }
  }

  for (int i = 0; i < n_designs; ++i) {
    bool is_premap = (trend_runtime && trend_runtime->has_premap()) ? map_next[i] : false;
    map_next[i] = !is_premap;
  }
  param_table.map_from_designs(designs, map_next);

  if (trend_runtime && trend_runtime->has_pretransform()) {
    const auto& specs_pretransform = trend_runtime->pretransform_specs();
    if (!specs_pretransform.empty()) {
      c_do_transform_pt(param_table, specs_pretransform);
    }

    for (std::size_t i = 0; i < trend_runtime->pretransform_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->pretransform_ops[i], param_table);
    }
  }

  transform_next = param_names_excluding(param_table, { &premap_set, &pretransform_set });
  auto postmap_specs = filter_specs_by_param_set(param_table, full_specs, transform_next);
  c_do_transform_pt(param_table, postmap_specs);

  if (trend_runtime && trend_runtime->has_posttransform()) {
    for (std::size_t i = 0; i < trend_runtime->posttransform_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->posttransform_ops[i], param_table);
    }
  }

  return param_table.materialize_by_param_names(keep_names);
}

static bool ddm_data_all_finite_untruncated(const Rcpp::DataFrame& data,
                                            const int n_trials) {
  Rcpp::NumericVector rts = data["rt"];
  Rcpp::IntegerVector R = data["R"];

  for (int i = 0; i < n_trials; ++i) {
    if (!R_FINITE(rts[i]) || Rcpp::NumericVector::is_na(rts[i])) return false;
    if (R[i] == NA_INTEGER) return false;
  }

  if (data.containsElementNamed("LT")) {
    Rcpp::NumericVector LT = data["LT"];
    for (int i = 0; i < n_trials; ++i) {
      if (LT[i] != 0.0) return false;
    }
  }

  if (data.containsElementNamed("UT")) {
    Rcpp::NumericVector UT = data["UT"];
    for (int i = 0; i < n_trials; ++i) {
      if (R_FINITE(UT[i])) return false;
    }
  }

  return true;
}

static bool race_data_all_finite_untruncated(const Rcpp::DataFrame& data,
                                             const int n_trials,
                                             const int n_lR) {
  if (n_trials <= 0 || n_lR <= 0 || (n_trials % n_lR) != 0) return false;

  Rcpp::NumericVector rts = data["rt"];
  Rcpp::IntegerVector R = data["R"];

  for (int start = 0; start < n_trials; start += n_lR) {
    if (!R_FINITE(rts[start]) || Rcpp::NumericVector::is_na(rts[start]) || rts[start] <= 0.0) {
      return false;
    }
    if (R[start] == NA_INTEGER) return false;
  }

  if (data.containsElementNamed("LT")) {
    Rcpp::NumericVector LT = data["LT"];
    for (int start = 0; start < n_trials; start += n_lR) {
      if (LT[start] != 0.0) return false;
    }
  }

  if (data.containsElementNamed("UT")) {
    Rcpp::NumericVector UT = data["UT"];
    for (int start = 0; start < n_trials; start += n_lR) {
      if (R_FINITE(UT[start])) return false;
    }
  }

  return true;
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
  RacePdf1Fun pdf1_ptr = nullptr;
  RaceCdf1Fun cdf1_ptr = nullptr;
  RaceRawFun model_dfun_raw = nullptr;   // fast-path: write log-density to pre-allocated buffer
  RaceRawFun model_pfun_raw = nullptr;   // fast-path: write log-survivor to pre-allocated buffer
  RaceLogSAtTFun logS_at_t_ptr = nullptr; // batch: log-survivor at scalar t for truncation norms
  ContextForRaceModels ctx;
};

static inline RaceModelAdapter resolve_race_model_adapter(const std::string& type_std,
                                                          const std::string& caller) {
  RaceModelAdapter out;
  out.ctx.min_lik_for_pdf = std::exp(std::log(1e-10));
  out.ctx.use_posdrift = true;
  out.ctx.gng = false;

  if (type_std.find("LBA") != std::string::npos) {
    out.pdf1_ptr = &dlba_scalar;
    out.cdf1_ptr = &plba_scalar;
    out.model_dfun_raw = &dlba_raw;
    out.model_pfun_raw = &plba_raw;
    out.logS_at_t_ptr = &lba_logS_at_t;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("RDM") != std::string::npos) {
    out.pdf1_ptr = &drdm_scalar;
    out.cdf1_ptr = &prdm_scalar;
    out.model_dfun_raw = &drdm_raw;
    out.model_pfun_raw = &prdm_raw;
    out.logS_at_t_ptr = &rdm_logS_at_t;
  } else if (type_std.find("LNR") != std::string::npos) {
    out.pdf1_ptr = &dlnr_scalar;
    out.cdf1_ptr = &plnr_scalar;
    out.model_dfun_raw = &dlnr_raw;
    out.model_pfun_raw = &plnr_raw;
    out.logS_at_t_ptr = &lnr_logS_at_t;
  } else {
    Rcpp::stop("Unsupported race model type string in %s: %s", caller.c_str(), type_std.c_str());
  }

  if (type_std.find("GNG") != std::string::npos) {
    out.ctx.gng = true;
  }
  return out;
}

// Pre-computed per-data state for c_log_likelihood_race.
// Built once outside the particle loop; reused across particles to eliminate
// per-particle R-heap allocations and repeated attribute/column reads.
struct RaceSharedState {
  bool valid = false;
  // Pre-read censoring/truncation bounds (Rcpp vectors keep memory alive)
  Rcpp::NumericVector LT_vec, UT_vec, LC_vec, UC_vec;
  // Pre-computed finite/other trial partition
  Rcpp::LogicalVector finite_mask;       // length n_trials; shared ref from dadm attr
  std::vector<int>    finite_unique_idx; // indices of finite unique trials
  std::vector<int>    other_unique_idx;  // indices of other (non-finite) unique trials
  // Pre-allocated mutable scratch buffers (length n_trials)
  std::vector<double> lds_buf;   // log-density; NOT re-initialised between particles
  std::vector<int>    idx_win;   // data-fixed winner mask  (finite rows only)
  std::vector<int>    idx_loss;  // data-fixed loser mask   (finite rows only)
  std::vector<int>    isok_buf;  // per-particle validity   (re-filled each call)
  bool any_win  = false;
  bool any_loss = false;
  // pContaminant column: -2=not yet searched, -1=absent, >=0=column index
  int  pc_col   = -2;
};

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
                            double min_ll, LogicalVector is_ok, bool gng,
                            bool all_finite_untruncated = false,
                            Rcpp::NumericVector* trial_ll_out = nullptr){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  const double* rt_ptr = rts.begin();
  const int* R_ptr = R.begin();
  const double* pars_cm = pars.begin();

  if (all_finite_untruncated) {
    std::vector<double> lls(static_cast<size_t>(n_trials), min_ll);
    std::vector<int> mask(static_cast<size_t>(n_trials), 1);
    std::vector<int> ok_int(static_cast<size_t>(n_trials), 0);
    for (int i = 0; i < n_trials; ++i) ok_int[static_cast<size_t>(i)] = is_ok[i] ? 1 : 0;
    d_DDM_Wien_raw(rt_ptr, R_ptr, pars_cm, n_trials, pars.ncol(),
                   mask.data(), ok_int.data(), lls.data(), min_ll);
    const double* lls_ptr = lls.data();
    const int* expand_ptr = expand.begin();
    double sum_ll = 0.0;

    if (trial_ll_out != nullptr) {
      for (int i = 0; i < n_out; ++i) {
        const int idx = expand_ptr[i] - 1; // expand is 1-based
        double v = lls_ptr[idx];
        if (!std::isfinite(v) || v < min_ll) {
          v = min_ll;
        }
        (*trial_ll_out)[i] = v;
        sum_ll += v;
      }
    } else {
      #pragma omp simd reduction(+:sum_ll)
      for (int i = 0; i < n_out; ++i) {
        const int idx = expand_ptr[i] - 1; // expand is 1-based
        double v = lls_ptr[idx];
        if (!std::isfinite(v) || v < min_ll) {
          v = min_ll;
        }
        sum_ll += v;
      }
    }

    return sum_ll;
  }

  NumericVector LT = get_col_with_default(data, "LT", 0.0);
  NumericVector UT = get_col_with_default(data, "UT", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  std::vector<double> lls(static_cast<size_t>(n_trials), R_NegInf);
  std::vector<int> finite_mask(static_cast<size_t>(n_trials), 0);
  std::vector<int> ok_int(static_cast<size_t>(n_trials), 0);
  std::vector<int> nonfinite_idx;
  bool any_ok_finite = false;
  bool any_ok_nonfinite = false;
  for (int i = 0; i < n_trials; ++i) {
    ok_int[static_cast<size_t>(i)] = is_ok[i] ? 1 : 0;
    const bool finite_i = R_FINITE(rts[i]);
    if (finite_i) {
      finite_mask[static_cast<size_t>(i)] = 1;
      if (is_ok[i]) any_ok_finite = true;
    } else {
      nonfinite_idx.push_back(i);
      if (is_ok[i]) any_ok_nonfinite = true;
    }
  }

  auto log_interval_mass = [&](double logF_low, double logF_high) -> double {
    if (!R_FINITE(logF_high)) return R_NegInf;
    if (logF_low == R_NegInf) return logF_high;
    return log_diff_exp(logF_high, logF_low);
  };
  
  if (any_ok_finite) {
    d_DDM_Wien_raw(rt_ptr, R_ptr, pars_cm, n_trials, pars.ncol(),
                   finite_mask.data(), ok_int.data(), lls.data(), min_ll);
  }

  // Apply truncation normaliser to finite trials for standard DDM only.
  // (For DDMGNG we preserve the existing semantics and do not apply truncation.)
  if (!gng && any_ok_finite) {
    // Early-out: skip all 4 expensive p_DDM_Wien (CDF) calls if no trial has truncation.
    // p_DDM_Wien with sv/SZ nonzero uses pdiff() with Neval=6000 numerical integration —
    // calling it unconditionally even when LT==0 and UT==Inf for all trials was the main
    // cause of ~100x slowdown vs the upstream SIMD branch for DDM/WDM.
    bool needs_truncation = false;
    for (int i = 0; i < n_trials; ++i) {
      if (is_ok[i] && finite_mask[static_cast<size_t>(i)] &&
          (LT[i] != 0.0 || R_FINITE(UT[i]))) {
        needs_truncation = true;
        break;
      }
    }
    if (needs_truncation) {
      // Extra scan: detect whether all truncated trials have LT==0 (very common —
      // avoids 2 of the 4 p_DDM_Wien CDF calls since pwiener/pdiff at t=0 always
      // returns R_NegInf when t0 > 0, so logF_LT_all is unconditionally R_NegInf).
      // Likewise detect whether any truncated trial has finite UT; if none do, the
      // UT CDFs are superfluous (log F_any(Inf) = 0 is used analytically).
      bool all_LT_zero   = true;
      bool any_finite_UT = false;
      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i] || !finite_mask[static_cast<size_t>(i)]) continue;
        if (LT[i] == 0.0 && !R_FINITE(UT[i])) continue; // trial not truncated
        if (LT[i] != 0.0) all_LT_zero   = false;
        if (R_FINITE(UT[i])) any_finite_UT = true;
      }

      LogicalVector ok_finite(n_trials);
      for (int i = 0; i < n_trials; ++i) {
        ok_finite[i] = is_ok[i] && finite_mask[static_cast<size_t>(i)];
      }
      IntegerVector R1_for_cdf(n_trials, 1);
      IntegerVector R2_for_cdf(n_trials, 2);

      // Compute only the CDF vectors that are actually required.
      NumericVector logF_LT_1, logF_LT_2, logF_UT_1, logF_UT_2;
      if (!all_LT_zero) {
        logF_LT_1 = p_DDM_Wien(LT, R1_for_cdf, pars, ok_finite);
        logF_LT_2 = p_DDM_Wien(LT, R2_for_cdf, pars, ok_finite);
      }
      if (any_finite_UT) {
        logF_UT_1 = p_DDM_Wien(UT, R1_for_cdf, pars, ok_finite);
        logF_UT_2 = p_DDM_Wien(UT, R2_for_cdf, pars, ok_finite);
      }

      for (int i = 0; i < n_trials; ++i) {
        if (!ok_finite[i]) continue;
        if (LT[i] == 0.0 && !R_FINITE(UT[i])) continue; // no truncation for this trial
        // logF_LT_all: R_NegInf when LT==0 (CDF(0-t0)<0→R_NegInf), otherwise from vector.
        const double logF_LT_all = all_LT_zero ? R_NegInf
            : log_sum_exp(logF_LT_1[i], logF_LT_2[i]);
        // logF_UT_all: 0 (= log 1) when UT==Inf, otherwise from vector.
        const double logF_UT_all = R_FINITE(UT[i])
            ? log_sum_exp(logF_UT_1[i], logF_UT_2[i]) : 0.0;
        const double logZ = log_interval_mass(logF_LT_all, logF_UT_all);
        if (R_FINITE(logZ)) lls[i] -= logZ;
        else lls[i] = min_ll;
      }
    } // end if (needs_truncation)
  }

  if (any_ok_nonfinite) {
    const int n_nonfinite = static_cast<int>(nonfinite_idx.size());
    NumericMatrix pars_nf(n_nonfinite, pars.ncol());
    LogicalVector ok_nf(n_nonfinite);
    for (int row = 0; row < n_nonfinite; ++row) {
      const int src = nonfinite_idx[static_cast<size_t>(row)];
      ok_nf[row] = is_ok[src];
      for (int col = 0; col < pars.ncol(); ++col) {
        pars_nf(row, col) = pars(src, col);
      }
    }
    Rcpp::List pars_dn = pars.attr("dimnames");
    if (!Rf_isNull(pars_dn)) {
      Rcpp::List nf_dn = Rcpp::List::create(R_NilValue, pars_dn[1]);
      pars_nf.attr("dimnames") = nf_dn;
    }
    auto compact_bound = [&](NumericVector source) {
      NumericVector out(n_nonfinite);
      for (int row = 0; row < n_nonfinite; ++row) {
        out[row] = source[nonfinite_idx[static_cast<size_t>(row)]];
      }
      return out;
    };
    NumericVector LT_nf = compact_bound(LT);
    NumericVector UT_nf = compact_bound(UT);
    NumericVector LC_nf = compact_bound(LC);
    NumericVector UC_nf = compact_bound(UC);

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
          if (finite_mask[static_cast<size_t>(i)] && R[i] != NA_INTEGER) {
            go_idx = R[i];
            break;
          }
        }
      }

      IntegerVector R_go_for_cdf(n_nonfinite, go_idx);
      NumericVector logcdf_U_go = p_DDM_Wien(UC_nf, R_go_for_cdf, pars_nf, ok_nf);
      NumericVector logcdf_L_go = p_DDM_Wien(LC_nf, R_go_for_cdf, pars_nf, ok_nf);

      for (int row = 0; row < n_nonfinite; ++row) {
        const int i = nonfinite_idx[static_cast<size_t>(row)];
        if (!is_ok[i]) continue;
        if (rts[i] == R_PosInf) {
          // Keep existing right-censoring behavior: no go-boundary hit before UC.
          lls[static_cast<size_t>(i)] = log1m_exp(logcdf_U_go[row]);
        } else if (rts[i] == R_NegInf) {
          // LC support: go-boundary hit by LC.
          lls[static_cast<size_t>(i)] = logcdf_L_go[row];
        } else {
          // NA non-response/censored coding can represent left or right censoring.
          lls[static_cast<size_t>(i)] = log_sum_exp(logcdf_L_go[row], log1m_exp(logcdf_U_go[row]));
        }
      }
    } else {
      // Standard DDM: classify non-finite trials by RT type and R value to avoid
      // computing CDFs that are never used.
      //
      //   has_left_cens  → need (LT, LC] CDFs
      //   has_right_cens → need (UC, UT] CDFs
      //   has_na_rt      → need both sides
      //   all_r_known + uniform R → compute only one boundary's CDFs
      bool has_left_cens  = false;  // rt == R_NegInf
      bool has_right_cens = false;  // rt == R_PosInf
      bool has_na_rt      = false;  // rt == NA
      bool all_r_known    = true;
      bool has_r1_only    = true;   // every r_known trial has R == 1
      bool has_r2_only    = true;   // every r_known trial has R == 2
      bool any_finite_UT  = false;  // at least one non-finite trial has finite UT

      for (int row = 0; row < n_nonfinite; ++row) {
        if (!ok_nf[row]) continue;
        const int src = nonfinite_idx[static_cast<size_t>(row)];
        const double rt = rts[src];
        if      (rt == R_NegInf) has_left_cens  = true;
        else if (rt == R_PosInf) has_right_cens = true;
        else                      has_na_rt      = true;
        if (R[src] == NA_INTEGER) {
          all_r_known = false; has_r1_only = false; has_r2_only = false;
        } else {
          if (R[src] != 1) has_r1_only = false;
          if (R[src] != 2) has_r2_only = false;
        }
        if (R_FINITE(UT_nf[row])) any_finite_UT = true;
      }

      // Which CDFs are required?
      const bool need_left  = has_left_cens  || has_na_rt;
      const bool need_right = has_right_cens || has_na_rt;
      // Boundary-specific: skip a boundary when all r_known trials share the other.
      // need_r1=false only when every trial has R==2; need_r2=false only when R==1.
      const bool need_r1 = !all_r_known || !has_r2_only;
      const bool need_r2 = !all_r_known || !has_r1_only;
      // UT CDFs: only needed for the right side AND when some trial has finite UT.
      // When UT == Inf the log-CDF is log(1) = 0 and needs no function call.
      const bool need_UT = need_right && any_finite_UT;

      IntegerVector R1_for_cdf(n_nonfinite, 1);
      IntegerVector R2_for_cdf(n_nonfinite, 2);

      NumericVector logF_LT_1, logF_LT_2, logF_LC_1, logF_LC_2;
      NumericVector logF_UC_1, logF_UC_2, logF_UT_1, logF_UT_2;

      if (need_left) {
        if (need_r1) logF_LT_1 = p_DDM_Wien(LT_nf, R1_for_cdf, pars_nf, ok_nf);
        if (need_r2) logF_LT_2 = p_DDM_Wien(LT_nf, R2_for_cdf, pars_nf, ok_nf);
        if (need_r1) logF_LC_1 = p_DDM_Wien(LC_nf, R1_for_cdf, pars_nf, ok_nf);
        if (need_r2) logF_LC_2 = p_DDM_Wien(LC_nf, R2_for_cdf, pars_nf, ok_nf);
      }
      if (need_right) {
        if (need_r1) logF_UC_1 = p_DDM_Wien(UC_nf, R1_for_cdf, pars_nf, ok_nf);
        if (need_r2) logF_UC_2 = p_DDM_Wien(UC_nf, R2_for_cdf, pars_nf, ok_nf);
        if (need_UT) {
          if (need_r1) logF_UT_1 = p_DDM_Wien(UT_nf, R1_for_cdf, pars_nf, ok_nf);
          if (need_r2) logF_UT_2 = p_DDM_Wien(UT_nf, R2_for_cdf, pars_nf, ok_nf);
        }
      }

      // Safe indexed access into a possibly-empty vector.
      // Returns R_NegInf when the vector was not computed (boundary not needed).
      auto lF = [](const NumericVector& v, int row) -> double {
        return v.size() > 0 ? v[row] : R_NegInf;
      };
      // UT log-CDF: log(1) = 0 when UT == Inf, otherwise from pre-computed vector.
      auto lF_UT = [&](const NumericVector& v, int row) -> double {
        return R_FINITE(UT_nf[row]) ? (v.size() > 0 ? v[row] : R_NegInf) : 0.0;
      };

      for (int row = 0; row < n_nonfinite; ++row) {
        const int i = nonfinite_idx[static_cast<size_t>(row)];
        if (!is_ok[i]) continue;

        const bool r_known = (R[i] != NA_INTEGER);
        const int  r_idx   = r_known ? R[i] : 0;

        if (rts[i] == R_NegInf) {
          // Left-censored: probability in (LT, LC]
          if (!r_known) {
            const double logF_LT_all = log_sum_exp(lF(logF_LT_1, row), lF(logF_LT_2, row));
            const double logF_LC_all = log_sum_exp(lF(logF_LC_1, row), lF(logF_LC_2, row));
            lls[static_cast<size_t>(i)] = log_interval_mass(logF_LT_all, logF_LC_all);
          } else if (r_idx == 1) {
            lls[static_cast<size_t>(i)] = log_interval_mass(lF(logF_LT_1, row), lF(logF_LC_1, row));
          } else {
            lls[static_cast<size_t>(i)] = log_interval_mass(lF(logF_LT_2, row), lF(logF_LC_2, row));
          }
        } else if (rts[i] == R_PosInf) {
          // Right-censored: probability in (UC, UT]
          if (!r_known) {
            const double logF_UC_all = log_sum_exp(lF(logF_UC_1, row), lF(logF_UC_2, row));
            const double logF_UT_all = R_FINITE(UT[i])
                ? log_sum_exp(lF_UT(logF_UT_1, row), lF_UT(logF_UT_2, row)) : 0.0;
            lls[static_cast<size_t>(i)] = log_interval_mass(logF_UC_all, logF_UT_all);
          } else if (r_idx == 1) {
            lls[static_cast<size_t>(i)] = log_interval_mass(lF(logF_UC_1, row), lF_UT(logF_UT_1, row));
          } else {
            lls[static_cast<size_t>(i)] = log_interval_mass(lF(logF_UC_2, row), lF_UT(logF_UT_2, row));
          }
        } else {
          // rt == NA: union of left (LT, LC] and right (UC, UT] intervals
          if (!r_known) {
            const double logF_LT_all = log_sum_exp(lF(logF_LT_1, row), lF(logF_LT_2, row));
            const double logF_LC_all = log_sum_exp(lF(logF_LC_1, row), lF(logF_LC_2, row));
            const double logF_UC_all = log_sum_exp(lF(logF_UC_1, row), lF(logF_UC_2, row));
            const double logF_UT_all = R_FINITE(UT[i])
                ? log_sum_exp(lF_UT(logF_UT_1, row), lF_UT(logF_UT_2, row)) : 0.0;
            lls[static_cast<size_t>(i)] = log_sum_exp(
                log_interval_mass(logF_LT_all, logF_LC_all),
                log_interval_mass(logF_UC_all, logF_UT_all));
          } else if (r_idx == 1) {
            lls[static_cast<size_t>(i)] = log_sum_exp(
                log_interval_mass(lF(logF_LT_1, row), lF(logF_LC_1, row)),
                log_interval_mass(lF(logF_UC_1, row), lF_UT(logF_UT_1, row)));
          } else {
            lls[static_cast<size_t>(i)] = log_sum_exp(
                log_interval_mass(lF(logF_LT_2, row), lF(logF_LC_2, row)),
                log_interval_mass(lF(logF_UC_2, row), lF_UT(logF_UT_2, row)));
          }
        }
      }
    }
  }
  const double* lls_ptr = lls.data();
  const int* expand_ptr = expand.begin();
  double sum_ll = 0.0;

  if (trial_ll_out != nullptr) {
    for (int i = 0; i < n_out; ++i) {
      const int idx = expand_ptr[i] - 1; // expand is 1-based
      double v = lls_ptr[idx];
      if (!std::isfinite(v) || v < min_ll) {
        v = min_ll;
      }
      (*trial_ll_out)[i] = v;
      sum_ll += v;
    }
  } else {
    #pragma omp simd reduction(+:sum_ll)
    for (int i = 0; i < n_out; ++i) {
      const int idx = expand_ptr[i] - 1; // expand is 1-based
      double v = lls_ptr[idx];
      if (!std::isfinite(v) || v < min_ll) {
        v = min_ll;
      }
      sum_ll += v;
    }
  }

  return sum_ll;
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
    const bool all_finite_untruncated = ddm_data_all_finite_untruncated(data, n_trials);
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
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok,
                                    gng, all_finite_untruncated);
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
                                       adapter.pdf1_ptr, adapter.cdf1_ptr,
                                       n_trials,
                                       winner, expand, min_ll, is_ok, n_lR,
                                       &adapter.ctx,
                                       all_finite_trials,
                                       adapter.model_dfun_raw,
                                       adapter.model_pfun_raw,
                                       adapter.logS_at_t_ptr);
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

  const bool is_ddm_type = type_std.find("DDM") != std::string::npos;
  const bool is_mri_type = (type == "MRI" || type == "MRI_AR1");
  const bool is_ss_type = type_std.find("SS") != std::string::npos;
  // For now, SS models (SSEXG, SSRDEX) use Ttransform in R which isn't easily
  // handled by the cached-spec fast path.
  const bool use_pt_mapping = is_ddm_type || (!is_mri_type && !is_ss_type);

  NumericMatrix one_particle(1, particle_matrix.ncol());
  colnames(one_particle) = p_names;
  IntegerVector kernel_output_codes = IntegerVector::create(1);
  auto pars_for_particle_generic = [&](int i) -> NumericMatrix {
    for (int j = 0; j < particle_matrix.ncol(); ++j) {
      one_particle(0, j) = particle_matrix(i, j);
    }
    return get_pars_c_wrapper_oo_core(one_particle, data, constants, designs, bounds, transforms,
                                      pretransforms, trend, false, false, kernel_output_codes);
  };

  NumericMatrix particle_matrix_pt;
  ParamTable param_table_template;
  std::vector<TransformSpec> transform_specs_pt;
  std::unique_ptr<TrendPlan> trend_plan;
  std::unique_ptr<TrendRuntime> trend_runtime;
  Rcpp::CharacterVector keep_names;
  std::vector<int> pm_col_to_base_idx;

  if (use_pt_mapping) {
    std::vector<TransformSpec> pre_specs = make_transform_specs(particle_matrix, pretransforms);
    particle_matrix_pt = c_do_transform(particle_matrix, pre_specs);

    const bool has_constants = !(constants.size() == 1 &&
                                 Rcpp::NumericVector::is_na(constants[0]));
    if (has_constants) {
      particle_matrix_pt = add_constants_columns(particle_matrix_pt, constants);
    }

    NumericVector p_vector = particle_matrix_pt(0, _);
    p_vector.attr("names") = colnames(particle_matrix_pt);
    param_table_template = ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials);
    transform_specs_pt = make_transform_specs_for_paramtable(param_table_template, transforms);

    Rcpp::CharacterVector pm_names = colnames(particle_matrix_pt);
    pm_col_to_base_idx.assign(pm_names.size(), -1);
    for (int j = 0; j < pm_names.size(); ++j) {
      std::string nm = Rcpp::as<std::string>(pm_names[j]);
      auto it = param_table_template.name_to_base_idx.find(nm);
      if (it != param_table_template.name_to_base_idx.end()) {
        pm_col_to_base_idx[j] = it->second;
      }
    }

    if (!trend.isNull()) {
      trend_plan.reset(new TrendPlan(trend, data));
      trend_runtime.reset(new TrendRuntime(*trend_plan));
      trend_runtime->bind_all_ops_to_paramtable(param_table_template);
      trend_runtime->init_cached_specs(param_table_template, transform_specs_pt);

      Rcpp::CharacterVector dnames = designs.names();
      const auto& trend_params = trend_runtime->all_trend_params();
      keep_names = names_excluding(dnames, { &trend_params });
    } else {
      keep_names = designs.names();
    }
  }

  TrendRuntime* trend_runtime_ptr = trend_runtime ? trend_runtime.get() : nullptr;

  if (is_ddm_type) {
    bool gng = (type_std.find("GNG") != std::string::npos);
    IntegerVector expand = data.attr("expand");
    const bool all_finite_untruncated = ddm_data_all_finite_untruncated(data, n_trials);
    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) {
        param_table_template.fill_from_particle_row(particle_matrix_pt, i, pm_col_to_base_idx);
      }
      pars = get_pars_matrix_oo_fast(param_table_template,
                                     designs,
                                     trend_runtime_ptr,
                                     transform_specs_pt,
                                     keep_names);
      if (i == 0) {
        bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      }
      is_ok = c_do_bound_pt(param_table_template, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok,
                                    gng, all_finite_untruncated);
    }
  } else if (is_mri_type) {
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle_generic(i);
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
  } else if (is_ss_type) {
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
      pars = pars_for_particle_generic(i);
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

    bool has_RACE_col_fp = data.containsElementNamed("RACE");
    bool has_RACE_attrs_fp = false;
    Rcpp::IntegerVector RACE_fp;
    Rcpp::LogicalVector RACE_mask_fp;
    if (has_RACE_col_fp &&
        data.hasAttribute("RACE_nacc_by_row") &&
        data.hasAttribute("RACE_mask")) {
      RACE_fp = data.attr("RACE_nacc_by_row");
      RACE_mask_fp = data.attr("RACE_mask");
      has_RACE_attrs_fp = (RACE_fp.size() == n_trials && RACE_mask_fp.size() == n_trials);
    }

    // --- Fast path: all-finite RTs, no truncation ---
    // Pre-allocates buffers once outside the particle loop; each particle writes
    // log-density/log-survivor directly into raw double arrays, bypassing all
    // Rcpp vector allocations inside c_log_likelihood_race.
    // Variable-accumulator RACE designs are supported when cache attrs are
    // present: inactive accumulator rows are masked out and not included in the
    // unique-trial sum.
    // pContaminant is handled inline: log(1-pC) shift for finite RTs (no-op when pC=0).
    const bool use_raw_fast_path =
      race_data_all_finite_untruncated(data, n_trials, n_lR) &&
      (!has_RACE_col_fp || has_RACE_attrs_fp) &&
      adapter.model_pfun_raw != nullptr &&
      adapter.model_dfun_raw != nullptr;

    const int n_unique_fp = n_trials / n_lR;
    std::vector<double> lds_buf;
    std::vector<double> ll_uniq_buf;
    std::vector<int> winner_int_buf;
    std::vector<int> loser_int_buf;
    std::vector<int> isok_int_fp;
    NumericVector rts_fp_hold;
    const double* rt_ptr = nullptr;
    const int* expand_ptr = expand.begin();
    const int n_exp = expand.length();
    int fast_pc_col = -2;  // -2: not yet searched; -1: not present in pars

    if (use_raw_fast_path) {
      lds_buf.resize(n_trials);
      ll_uniq_buf.resize(n_unique_fp);
      winner_int_buf.resize(n_trials);
      loser_int_buf.resize(n_trials);
      isok_int_fp.resize(n_trials);
      // winner/loser masks are data-fixed; fill once
      for (int j = 0; j < n_trials; ++j) {
        const bool active = !has_RACE_col_fp || RACE_mask_fp[j];
        winner_int_buf[j] = (active && winner[j]) ? 1 : 0;
        loser_int_buf[j]  = (active && n_lR > 1 && !winner[j]) ? 1 : 0;
      }
      rts_fp_hold = data["rt"];
      rt_ptr = rts_fp_hold.begin();
    }

    // --- Build shared state for the mixed (non-raw-fast) path ---
    // Pre-compute all data-fixed structures once so c_log_likelihood_race can skip
    // per-particle R-heap allocations, column reads, and attribute lookups.
    RaceSharedState race_shared;
    if (!use_raw_fast_path && adapter.model_dfun_raw != nullptr && !all_finite_trials) {
      const bool has_part =
        data.hasAttribute("finite_rt_mask") &&
        data.hasAttribute("finite_rt_unique_trial_indices") &&
        data.hasAttribute("other_unique_trial_indices");
      if (has_part) {
        race_shared.finite_mask = data.attr("finite_rt_mask");
        Rcpp::IntegerVector fa = data.attr("finite_rt_unique_trial_indices");
        Rcpp::IntegerVector oa = data.attr("other_unique_trial_indices");
        race_shared.finite_unique_idx.assign(fa.begin(), fa.end());
        race_shared.other_unique_idx.assign(oa.begin(), oa.end());

        // Pre-read censoring/truncation bounds once
        race_shared.LT_vec = get_col_with_default(data, "LT", 0.0);
        race_shared.UT_vec = get_col_with_default(data, "UT", R_PosInf);
        race_shared.LC_vec = get_col_with_default(data, "LC", 0.0);
        race_shared.UC_vec = get_col_with_default(data, "UC", R_PosInf);

        // Pre-allocate mutable scratch buffers
        race_shared.lds_buf.resize(static_cast<size_t>(n_trials)); // no init needed
        race_shared.idx_win.assign(static_cast<size_t>(n_trials), 0);
        race_shared.idx_loss.assign(static_cast<size_t>(n_trials), 0);
        race_shared.isok_buf.resize(static_cast<size_t>(n_trials));

        // Fill data-fixed winner/loser masks (finite trials only)
        const Rcpp::LogicalVector& fmask = race_shared.finite_mask;
        for (int j = 0; j < n_trials; ++j) {
          if (!fmask[j]) continue;
          if (has_RACE_col_fp && !RACE_mask_fp[j]) continue;
          if (winner[j]) {
            race_shared.idx_win[static_cast<size_t>(j)] = 1;
            race_shared.any_win = true;
          } else if (n_lR > 1) {
            race_shared.idx_loss[static_cast<size_t>(j)] = 1;
            race_shared.any_loss = true;
          }
        }
        race_shared.valid = true;
      }
    }

    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) {
        param_table_template.fill_from_particle_row(particle_matrix_pt, i, pm_col_to_base_idx);
      }
      pars = get_pars_matrix_oo_fast(param_table_template,
                                     designs,
                                     trend_runtime_ptr,
                                     transform_specs_pt,
                                     keep_names);
      if (i == 0) {
        bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      }
      is_ok = c_do_bound_pt(param_table_template, bound_specs);
      is_ok = lr_all(is_ok, n_lR);

      if (use_raw_fast_path) {
        // Locate pContaminant column once from first particle's pars dimnames
        if (fast_pc_col == -2) {
          fast_pc_col = -1;
          Rcpp::List dn = pars.attr("dimnames");
          if (!Rf_isNull(dn)) {
            Rcpp::CharacterVector cn = Rcpp::as<Rcpp::CharacterVector>(dn[1]);
            for (int j = 0; j < cn.size(); ++j) {
              if (Rcpp::as<std::string>(cn[j]) == "pContaminant") { fast_pc_col = j; break; }
            }
          }
        }

        // Fill per-particle isok buffer
        for (int j = 0; j < n_trials; ++j) isok_int_fp[j] = is_ok[j] ? 1 : 0;

        const double* pars_cm = pars.begin();  // column-major: col j at pars_cm + j*n_trials

        // Log-density for winner rows
        adapter.model_dfun_raw(rt_ptr, pars_cm, n_trials,
                               winner_int_buf.data(), isok_int_fp.data(),
                               lds_buf.data(), min_ll, &adapter.ctx);

        // Log-survivor for loser rows (writes into same buffer, different slots)
        if (n_lR > 1) {
          adapter.model_pfun_raw(rt_ptr, pars_cm, n_trials,
                                 loser_int_buf.data(), isok_int_fp.data(),
                                 lds_buf.data(), min_ll, &adapter.ctx);
        }

        // Sum n_lR rows per unique trial; apply pC correction (no-op when pC=0).
        // For all-finite RTs: ll_trial = log(1-pC) + sum_k(lds[k])
        for (int j = 0; j < n_unique_fp; ++j) {
          double s = 0.0;
          const int base = j * n_lR;
          const int n_lR_j = has_RACE_col_fp ? RACE_fp[base] : n_lR;
          for (int k = 0; k < n_lR_j; ++k) s += lds_buf[base + k];
          s = s < min_ll ? min_ll : s;
          if (fast_pc_col >= 0) {
            const double pC = pars_cm[fast_pc_col * n_trials + base];
            if (pC != 0.0) s += std::log1p(-pC);
          }
          ll_uniq_buf[j] = s;
        }

        // Expand unique-trial LLs and accumulate (1-based expand indices)
        double total_ll = 0.0;
        #pragma omp simd reduction(+:total_ll)
        for (int ei = 0; ei < n_exp; ++ei) {
          total_ll += ll_uniq_buf[expand_ptr[ei] - 1];
        }
        lls[i] = total_ll;
      } else {
        lls[i] = c_log_likelihood_race(pars, data,
                                       adapter.pdf1_ptr, adapter.cdf1_ptr,
                                       n_trials,
                                       winner, expand, min_ll, is_ok, n_lR,
                                       &adapter.ctx,
                                       all_finite_trials,
                                       adapter.model_dfun_raw,
                                       adapter.model_pfun_raw,
                                       adapter.logS_at_t_ptr,
                                       race_shared.valid ? &race_shared : nullptr);
      }
    }
  }
  return lls;
}

// [[Rcpp::export]]
NumericMatrix calc_ll_oo_pw(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                            List designs, String type, List bounds, List transforms, List pretransforms,
                            CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue) {
  const int n_particles = particle_matrix.nrow();
  const int n_trials = data.nrow();
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
    const int n_out = expand.length();
    const bool all_finite_untruncated = ddm_data_all_finite_untruncated(data, n_trials);
    NumericMatrix result(n_particles, n_out);
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      is_ok = c_do_bound(pars, bound_specs);
      NumericVector row_vec(n_out);
      c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok,
                           gng, all_finite_untruncated, &row_vec);
      result(i, _) = row_vec;
    }
    return result;
  } else if (type == "MRI" || type == "MRI_AR1" ||
             type_std.find("SS") != std::string::npos ||
             type_std.find("SOFTMAX") != std::string::npos) {
    Rcpp::stop("calc_ll_oo_pw: not implemented for model type '%s'", type_std.c_str());
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll_oo_pw");
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

    const int n_out_race = (expand.length() > 0) ? expand.length() : (n_trials / n_lR);
    NumericMatrix result(n_particles, n_out_race);
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      NumericVector row_vec(n_out_race);
      c_log_likelihood_race(pars, data,
                            adapter.pdf1_ptr, adapter.cdf1_ptr,
                            n_trials,
                            winner, expand, min_ll, is_ok, n_lR,
                            &adapter.ctx,
                            all_finite_trials,
                            adapter.model_dfun_raw,
                            adapter.model_pfun_raw,
                            adapter.logS_at_t_ptr,
                            nullptr,
                            &row_vec);
      result(i, _) = row_vec;
    }
    return result;
  }

  return NumericMatrix(0, 0);
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
  // When LT == 0, every proper distribution has CDF(0) = 0 → S(0) = 1 → log = 0.
  // Skip the n_lR scalar cdf1 calls in this common case.
  double logS_LT;
  if (LT == 0.0) {
    logS_LT = 0.0;
  } else {
    logS_LT = log_survivor_rowmajor(LT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
    if (!R_FINITE(logS_LT)) return R_NegInf;
  }
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
    RacePdf1Fun pdf1,                       // Scalar PDF (for GSL)
    RaceCdf1Fun cdf1,                       // Scalar CDF (for GSL)
    const int n_trials,
    LogicalVector winner,
    Rcpp::IntegerVector expand,  // Vector for expanding unique LLs to full trial count
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector isok,   // Parameter validity for each row in 'pars' matrix
    int n_lR,                              // Number of accumulators in the race (must be > 0 if data exists)
    void* model_context_for_funcs,          // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
    bool all_finite_trials,             // Data-only hint: all trials finite/in-bounds/known response
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,            // batch log-survivor at scalar t (for truncation norms)
    RaceSharedState* shared,             // optional pre-computed per-data state (nullptr = compute per-call)
    NumericVector* trial_ll_out
) {

  // Only allocate a GSL workspace if/when numerical integration is needed.
  GslWorkspacePtr workspace(nullptr, &gsl_integration_workspace_free);

  const bool use_shared = (shared != nullptr && shared->valid);

  // Fetch censoring and truncation columns only when we actually need them.
  // When all_finite_trials=true the caller guarantees LT=0, UT=Inf, no censoring —
  // deferring these reads avoids O(n_trials) allocation+fill per particle call.
  // When use_shared=true the pre-read vectors are reused directly (no allocation).
  const bool may_need_ct = !all_finite_trials;
  Rcpp::NumericVector LT, UT, LC, UC;
  if (may_need_ct) {
    if (use_shared) {
      LT = shared->LT_vec;   // shared reference, no copy
      UT = shared->UT_vec;
      LC = shared->LC_vec;
      UC = shared->UC_vec;
    } else {
      LT = get_col_with_default(dadm, "LT", 0.0);
      UT = get_col_with_default(dadm, "UT", R_PosInf);
      LC = get_col_with_default(dadm, "LC", 0.0);
      UC = get_col_with_default(dadm, "UC", R_PosInf);
    }
  }
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
  // Use pre-allocated std::vector when shared state is available; avoids per-particle
  // R-heap allocation.  The buffer does NOT need re-initialisation between particles:
  // dfun_raw/pfun_raw overwrite every finite-trial slot before it is read, and
  // other-trial slots are computed independently in the "other trials" loop below.
  std::vector<double> lds_local;
  double* lds_ptr;
  if (use_shared && static_cast<int>(shared->lds_buf.size()) == n_trials) {
    lds_ptr = shared->lds_buf.data();
  } else {
    lds_local.assign(static_cast<size_t>(n_trials), min_ll);
    lds_ptr = lds_local.data();
  }
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
  // The column index is data-structure-fixed so we cache it in shared state after the
  // first search.  The "all zeroes" check still runs per particle (values change).
  bool use_pC = false;
  int pc_col = -1;
  if (use_shared && shared->pc_col != -2) {
    // Fast path: use cached result from a previous particle call
    pc_col = shared->pc_col;
    use_pC = (pc_col >= 0);
  } else {
    Rcpp::List dimnames = pars.attr("dimnames");
    Rcpp::CharacterVector colnames = as<Rcpp::CharacterVector>(dimnames[1]);
    for (int j = 0; j < colnames.size(); ++j) {
      if (as<std::string>(colnames[j]) == "pContaminant") {
        pc_col = j;
        use_pC = true;
        break;
      }
    }
    if (use_shared) shared->pc_col = pc_col; // cache for subsequent particles
  }
  // Per-particle check: pC column present but all values zero → treat as absent
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
  // Use std::vector to avoid per-particle R-heap allocation overhead.
  std::vector<double> ll_unique(static_cast<size_t>(n_unique_trials), min_ll);
  std::vector<double> pC_values(static_cast<size_t>(n_unique_trials), 0.0);
  // Raw pointer into pars matrix: col-major layout, element (row,col) = pars_cm_ptr[col*n_trials+row]
  const double* pars_cm_ptr = pars.begin();
  if (use_pC) {
    for (int j = 0; j < n_unique_trials; ++j) {
      pC_values[static_cast<size_t>(j)] = pars_cm_ptr[static_cast<size_t>(pc_col) * n_trials + j * n_lR];
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

  // fill_trial_buffers: copies one trial's params into row-major scratch (for GSL/rowmajor helpers).
  // Uses raw column-major pointer to avoid Rcpp subscript overhead.
  auto fill_trial_buffers = [&](int start_row_idx, int n_lR_j) {
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start_row_idx + k;
      isok_int_buffer[static_cast<size_t>(k)] = isok[row] ? 1 : 0;
      for (int c = 0; c < n_par; ++c)
        pars_rowmajor_buffer[static_cast<size_t>(k) * n_par + c] =
            pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
    }
  };

  // Fast path hint: computed once per calc_ll call (data-only), to avoid re-scanning per particle.
  const bool use_full_finite_batch = all_finite_trials;
  Rcpp::LogicalVector finite_rt_mask;
  std::vector<int> finite_rt_unique_trial_indices;
  std::vector<int> other_unique_trial_indices;
  if (!use_full_finite_batch) {
    if (use_shared && static_cast<int>(shared->finite_mask.size()) == n_trials) {
      // Fast path: re-use pre-read partition from shared state (no attr lookup, no copy)
      finite_rt_mask = shared->finite_mask;
      finite_rt_unique_trial_indices = shared->finite_unique_idx;
      other_unique_trial_indices     = shared->other_unique_idx;
    } else {
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
  
  // log_surv_cm: log S(t) = sum_k log(1-F_k(t)) read directly from column-major pars.
  // Used by the "other trials" path for analytical survivor calls — avoids the
  // row-major copy that fill_trial_buffers performs.
  auto log_surv_cm = [&](double t, int start_row_idx, int n_lR_j) -> double {
    double logS = 0.0;
    double par_buf[16];
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start_row_idx + k;
      if (!isok[row]) return R_NegInf;
      for (int c = 0; c < n_par; ++c)
        par_buf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
      const double Fk = cdf1(t, par_buf, model_context_for_funcs);
      if (Fk >= 1.0) return R_NegInf;
      logS += std::log1p(-Fk);
    }
    return logS;
  };

  const bool has_finite_batch = use_full_finite_batch || (finite_rt_unique_trial_indices.size() > 0);
  if (has_finite_batch) {
    bool any_win = false;
    bool any_loss = false;
    // Local fallback storage (used when shared state is not available)
    std::vector<int> idx_win_int_local;
    std::vector<int> idx_loss_int_local;
    std::vector<int> isok_int_all_local;
    // Pointers resolved below (point into shared or local storage)
    const int* idx_win_ptr  = nullptr;
    const int* idx_loss_ptr = nullptr;
    int*       isok_ptr     = nullptr;

    const bool use_shared_bufs = use_shared &&
      static_cast<int>(shared->idx_win.size()) == n_trials;

    if (use_shared_bufs) {
      // Per-particle: only fill isok (winner/loser masks are data-fixed)
      for (int i = 0; i < n_trials; ++i)
        shared->isok_buf[static_cast<size_t>(i)] = isok[i] ? 1 : 0;
      idx_win_ptr  = shared->idx_win.data();
      idx_loss_ptr = shared->idx_loss.data();
      isok_ptr     = shared->isok_buf.data();
      any_win      = shared->any_win;
      any_loss     = shared->any_loss;
    } else {
      idx_win_int_local.assign(static_cast<size_t>(n_trials), 0);
      idx_loss_int_local.assign(static_cast<size_t>(n_trials), 0);
      isok_int_all_local.assign(static_cast<size_t>(n_trials), 0);
      for (int i = 0; i < n_trials; ++i) {
        isok_int_all_local[static_cast<size_t>(i)] = isok[i] ? 1 : 0;
        if (has_RACE_col && !RACE_mask[i]) continue;
        if (!use_full_finite_batch && !finite_rt_mask[i]) continue;
        if (winner[i]) {
          idx_win_int_local[static_cast<size_t>(i)] = 1;
          any_win = true;
        } else if (n_lR > 1) {
          idx_loss_int_local[static_cast<size_t>(i)] = 1;
          any_loss = true;
        }
      }
      idx_win_ptr  = idx_win_int_local.data();
      idx_loss_ptr = idx_loss_int_local.data();
      isok_ptr     = isok_int_all_local.data();
    }

    const double* rt_ptr = rts_dadm.begin();
    const double* pars_cm = pars.begin();
    if (any_win) {
      model_dfun_raw(rt_ptr, pars_cm, n_trials,
                     idx_win_ptr, isok_ptr,
                     lds_ptr, min_ll, dense_ctx_ptr);
    }
    if (n_lR > 1 && any_loss) {
      model_pfun_raw(rt_ptr, pars_cm, n_trials,
                     idx_loss_ptr, isok_ptr,
                     lds_ptr, min_ll, dense_ctx_ptr);
    }

    // --- Pre-compute truncation normalisers in batch when possible ---
    // When may_need_ct AND logS_at_t is available, we can replace per-trial
    // fill_trial_buffers + get_trunc_normaliser_rowmajor_cpp with a single
    // column-major pass: eliminates n_lR*n_par row copies per truncated trial.
    //
    // Fast-batch conditions: LT == 0 everywhere (logS_LT = 0 analytically),
    // UT is uniform across all truncated trials, and logS_at_t_ptr is provided.
    // Trials that fail the fast path fall back to the per-trial scalar route.
    bool batch_trunc_done = false;
    std::vector<double> logZ_batch;  // length n_unique_trials, filled when batch_trunc_done
    const int n_finite_unique = use_full_finite_batch
        ? n_unique_trials
        : static_cast<int>(finite_rt_unique_trial_indices.size());

    if (may_need_ct && logS_at_t != nullptr) {
      // Pass 1: scan for truncated trials; check uniformity of LT/UT.
      double uniform_UT = R_PosInf;
      bool all_LT_zero = true;
      bool uniform_UT_ok = true;
      bool any_trunc = false;
      std::vector<int> trunc_mask(static_cast<size_t>(n_unique_trials), 0);

      for (int i = 0; i < n_finite_unique; ++i) {
        const int j = use_full_finite_batch
            ? i : finite_rt_unique_trial_indices[static_cast<size_t>(i)];
        const int start = j * n_lR;
        const double LTj = LT[start];
        const double UTj = UT[start];
        if (LTj != 0.0 || UTj != R_PosInf) {
          trunc_mask[static_cast<size_t>(j)] = 1;
          any_trunc = true;
          if (LTj != 0.0) all_LT_zero = false;
          if (uniform_UT == R_PosInf) uniform_UT = UTj;
          else if (UTj != uniform_UT) uniform_UT_ok = false;
        }
      }

      if (any_trunc && all_LT_zero && uniform_UT_ok && uniform_UT != R_PosInf) {
        // Pass 2: batch-compute log S(UT) for all truncated trials at once.
        // logS_LT = 0 (LT==0), so logZ = log(1 - exp(logS_UT)) = log_diff_exp(0, logS_UT).
        // Reuse shared isok buffer when available (already filled above); else build locally.
        std::vector<int> isok_all_int_local;
        const int* isok_all_ptr_trunc;
        if (isok_ptr != nullptr) {
          isok_all_ptr_trunc = isok_ptr;  // already filled per-particle above
        } else {
          isok_all_int_local.resize(static_cast<size_t>(n_trials));
          for (int r = 0; r < n_trials; ++r)
            isok_all_int_local[static_cast<size_t>(r)] = isok[r] ? 1 : 0;
          isok_all_ptr_trunc = isok_all_int_local.data();
        }

        std::vector<double> logS_UT_vec(static_cast<size_t>(n_unique_trials), 0.0);
        logS_at_t(uniform_UT,
                  pars.begin(), n_trials, n_lR, n_par,
                  trunc_mask.data(), n_unique_trials,
                  isok_all_ptr_trunc,
                  model_context_for_funcs,
                  logS_UT_vec.data());

        logZ_batch.assign(static_cast<size_t>(n_unique_trials), 0.0);
        for (int j = 0; j < n_unique_trials; ++j) {
          if (!trunc_mask[static_cast<size_t>(j)]) continue;
          const double logS_UT_j = logS_UT_vec[static_cast<size_t>(j)];
          if (!R_FINITE(logS_UT_j)) { logZ_batch[static_cast<size_t>(j)] = R_NegInf; continue; }
          // logZ = log(S(LT=0) - S(UT)) = log(1 - S(UT)) = log_diff_exp(0, logS_UT)
          const double logP = log_diff_exp(0.0, logS_UT_j);
          // Use the same epsilon threshold as the scalar fallback path.
          static const double kLogProbEps = std::log(std::numeric_limits<double>::epsilon());
          if (R_FINITE(logP) && logP > kLogProbEps) {
            logZ_batch[static_cast<size_t>(j)] = logP;
          } else {
            // Catastrophic cancellation: fall back to per-trial GSL path below.
            // Signal with NaN so the per-trial check recomputes it.
            logZ_batch[static_cast<size_t>(j)] = std::numeric_limits<double>::quiet_NaN();
          }
        }
        batch_trunc_done = true;
      }
    }

    // Apply truncation correction and calculate log-likelihood for each trial in the batch
    for (int i = 0; i < n_finite_unique; ++i) {
      // When all unique trials are finite, iterate in order 0..n_unique_trials-1.
      // Otherwise, iterate over the precomputed subset of finite-RT unique trials.
      const int unique_trial_idx = use_full_finite_batch
          ? i : finite_rt_unique_trial_indices[static_cast<size_t>(i)];
      const int start_row_idx = unique_trial_idx * n_lR;
      n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
      log_Z_this = 0.0;

      if (may_need_ct) {
        const double LTj = LT[start_row_idx];
        const double UTj = UT[start_row_idx];
        if (LTj != 0.0 || UTj != R_PosInf) { // truncation active
          bool need_scalar = true;
          if (batch_trunc_done) {
            const double bz = logZ_batch[static_cast<size_t>(unique_trial_idx)];
            if (!std::isnan(bz)) {       // NaN flags a fallback case
              log_Z_this = bz;
              need_scalar = false;
            }
          }
          if (need_scalar) {
            fill_trial_buffers(start_row_idx, n_lR_j);
            log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer.data(),
                                                           isok_int_buffer.data(),
                                                           pdf1, cdf1,
                                                           LTj, UTj,
                                                           n_lR_j, n_par,
                                                           gsl_ctl,
                                                           model_context_for_funcs,
                                                           workspace);
          }
          double current_trial_ll_sum = 0.0;
          for (int k = 0; k < n_lR_j; ++k)
            current_trial_ll_sum += lds_ptr[start_row_idx + k];
          if (NumericVector::is_na(log_Z_this) || !R_FINITE(log_Z_this)) {
            ll_unique[unique_trial_idx] = min_ll;
          } else {
            ll_unique[unique_trial_idx] = std::max(min_ll, current_trial_ll_sum - log_Z_this);
          }
          continue;
        }
      }

      // No truncation (or all_finite_trials guarantees none): plain sum.
      double current_trial_ll_sum = 0.0;
      for (int k = 0; k < n_lR_j; ++k)
        current_trial_ll_sum += lds_ptr[start_row_idx + k];
      ll_unique[unique_trial_idx] = std::max(min_ll, current_trial_ll_sum);
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
    
    // Lazy row-major buffer fill — only when GSL integration or row-major
    // helper functions are actually needed.  Pure survivor / CDF calls use
    // log_surv_cm which reads column-major directly.
    bool buffers_filled = false;
    auto ensure_buffers = [&]() {
      if (!buffers_filled) {
        fill_trial_buffers(start_row_idx, n_lR_j);
        buffers_filled = true;
      }
    };

    auto integrate_interval = [&](int k_winner_1based, double low, double upp) -> double {
      ensure_buffers();
      gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
      return integrate_for_kth_winner_rowmajor_cpp(k_winner_1based,
                                                   pars_rowmajor_buffer.data(),
                                                   isok_int_buffer.data(),
                                                   low, upp, pdf1, cdf1,
                                                   n_lR_j, n_par, gsl_ctl,
                                                   model_context_for_funcs, w);
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
          // Left-censoring with unknown winner: S(LT) - S(LC), S(t)=prod_k(1-F_k(t)).
          const double logP = log_diff_exp(log_surv_cm(lower_for_trial, start_row_idx, n_lR_j),
                                           log_surv_cm(upper_for_trial, start_row_idx, n_lR_j));
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
        const double logA = integrate_interval(k_nogo, lower_for_trial, upper_for_trial);
        const double logB = log_surv_cm(upper_for_trial, start_row_idx, n_lR_j); // incl. LBAIO mass at +Inf
        current_ll_val = log_sum_exp(logA, logB);
      } else {
        lower_for_trial = UCj;
        upper_for_trial = UTj;
        if (R_j_idx == NA_INTEGER) {
          const double logP = posdrift
            ? log_diff_exp(log_surv_cm(lower_for_trial, start_row_idx, n_lR_j),
                           log_surv_cm(upper_for_trial, start_row_idx, n_lR_j))
            : log_surv_cm(lower_for_trial, start_row_idx, n_lR_j); // LBAIO: incl. intrinsic-omission mass
          if (R_FINITE(logP) && logP > log_prob_eps) {
            current_ll_val = logP;
          } else { // numerical integration fallback
            for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
              current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
            }
            if (!posdrift) { // intrinsic-omission models only
              ensure_buffers();
              const double log_p_I = log_pIO_rowmajor(pars_rowmajor_buffer.data(),
                                                      isok_int_buffer.data(),
                                                      n_lR_j, n_par);
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
        const double logP1 = log_diff_exp(log_surv_cm(lower1, start_row_idx, n_lR_j),
                                          log_surv_cm(upper1, start_row_idx, n_lR_j));
        const double logP2 = posdrift
          ? log_diff_exp(log_surv_cm(lower2, start_row_idx, n_lR_j),
                         log_surv_cm(upper2, start_row_idx, n_lR_j))
          : log_surv_cm(lower2, start_row_idx, n_lR_j);
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
            ensure_buffers();
            const double log_p_I = log_pIO_rowmajor(pars_rowmajor_buffer.data(), isok_int_buffer.data(), n_lR_j, n_par);
            current_ll_val = log_sum_exp(current_ll_val, log_p_I);
          }
        }
      }
    } else if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER) {
      ensure_buffers();
      current_ll_val = log_min_density_rowmajor(rt_j,
                                                pars_rowmajor_buffer.data(),
                                                isok_int_buffer.data(),
                                                n_lR_j, n_par, pdf1, cdf1,
                                                model_context_for_funcs,
                                                logS_k_buffer);
    }

    if (current_ll_val > min_ll && has_trunc) {
      ensure_buffers();
      log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer.data(),
                                                     isok_int_buffer.data(),
                                                     pdf1, cdf1,
                                                     LTj, UTj, n_lR_j, n_par,
                                                     gsl_ctl,
                                                     model_context_for_funcs,
                                                     workspace);
      if (!R_FINITE(log_Z_this)) current_ll_val = min_ll;
      else current_ll_val -= log_Z_this;
    }
    
    ll_unique[unique_trial_idx] = std::max(min_ll, current_ll_val);
  }
  
  
  // --- Summation of log-likelihoods for all unique trials ---
  // pC modification (if any) is kept in a separate pass so the final summation
  // loop is a pure reduction — allowing #pragma omp simd to vectorize it.
  double total_ll = 0;
  if (expand.length() > 0) { // non-compressed dadm: sum via expand index vector
    if (use_pC) {
      for (int j = 0; j < n_unique_trials; ++j) {
        double pC = pC_values[j];
        double log1m_pC = std::log1p(-pC);
        int start_row_idx = j * n_lR;
        double rt_j = rts_dadm[start_row_idx];
        if (rt_j == R_PosInf) {
          ll_unique[j] = log_sum_exp(std::log(pC), log1m_pC + ll_unique[j]);
        } else {
          ll_unique[j] = log1m_pC + ll_unique[j];
        }
      }
    }
    const double* ll_ptr = ll_unique.data();
    const int* ex_ptr = expand.begin();
    const int n_exp = expand.length();
    if (trial_ll_out != nullptr) {
      for (int i = 0; i < n_exp; ++i) {
        const double v = ll_ptr[ex_ptr[i] - 1];
        (*trial_ll_out)[i] = v;
        total_ll += v;
      }
    } else {
      // Gather + reduce: indirect indexing prevents full SIMD gather, but the
      // reduction itself still benefits from the pragma.
      #pragma omp simd reduction(+:total_ll)
      for (int i = 0; i < n_exp; ++i) {
        total_ll += ll_ptr[ex_ptr[i] - 1];
      }
    }
  } else { // compressed dadm: each unique trial counted once
    if (use_pC) {
      for (int j = 0; j < n_unique_trials; ++j) {
        double pC = pC_values[j];
        double log1m_pC = std::log1p(-pC);
        int start_row_idx = j * n_lR;
        double rt_j = rts_dadm[start_row_idx];
        if (R_FINITE(rt_j)) {
          ll_unique[j] = log1m_pC + ll_unique[j];
        } else {
          ll_unique[j] = log_sum_exp(std::log(pC), log1m_pC + ll_unique[j]);
        }
      }
    }
    const double* ll_ptr = ll_unique.data();
    if (trial_ll_out != nullptr) {
      for (int j = 0; j < n_unique_trials; ++j) {
        const double v = ll_ptr[j];
        (*trial_ll_out)[j] = v;
        total_ll += v;
      }
    } else {
      // Pure sequential reduction — fully vectorizable.
      #pragma omp simd reduction(+:total_ll)
      for (int j = 0; j < n_unique_trials; ++j) {
        total_ll += ll_ptr[j];
      }
    }
  }
  
  return total_ll;
}

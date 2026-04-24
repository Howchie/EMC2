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
#include "wald_functions.h"
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
#include <cstring>
#include <mutex>

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

NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend);

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

static void update_pt_only(ParamTable& param_table,
                           const Rcpp::List& designs,
                           TrendRuntime* trend_runtime,
                           const std::vector<TransformSpec>& full_specs) {
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

  if (type_std.find("RDMSWTN") != std::string::npos) {
    // Must be checked before "RDM" since "RDMSWTN" contains "RDM"
    out.pdf1_ptr       = &drdmswtn_scalar;
    out.cdf1_ptr       = &prdmswtn_scalar;
    out.model_dfun_raw = &drdmswtn_raw;
    out.model_pfun_raw = &prdmswtn_raw;
    out.logS_at_t_ptr  = &rdmswtn_logS_at_t;
    out.ctx.t0_index   = 3;
  } else if (type_std.find("GBM") != std::string::npos) {
    // Must be checked before "RDM" since "RDMGBM" contains "RDM"
    out.pdf1_ptr       = &drdmgbm_scalar;
    out.cdf1_ptr       = &prdmgbm_scalar;
    out.model_dfun_raw = &drdmgbm_raw;
    out.model_pfun_raw = &prdmgbm_raw;
    out.logS_at_t_ptr  = &rdmgbm_logS_at_t;
    out.ctx.t0_index   = 3;
  } else if (type_std.find("BAwL") != std::string::npos) {
    out.pdf1_ptr       = &dbawl_scalar;
    out.cdf1_ptr       = &pbawl_scalar;
    out.model_dfun_raw = &dbawl_raw;
    out.model_pfun_raw = &pbawl_raw;
    out.logS_at_t_ptr  = &bawl_logS_at_t;
    out.ctx.t0_index   = 4;
    // Leaky ballistic accumulators can have defective upper tails (never-finish
    // mass) even when posdrift=TRUE.
    out.ctx.defective_upper_tail = true;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("LBA") != std::string::npos) {
    out.pdf1_ptr = &dlba_scalar;
    out.cdf1_ptr = &plba_scalar;
    out.model_dfun_raw = &dlba_raw;
    out.model_pfun_raw = &plba_raw;
    out.logS_at_t_ptr = &lba_logS_at_t;
    out.ctx.t0_index = 4;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
      out.ctx.defective_upper_tail = true;
    }
  } else if (type_std.find("RDM") != std::string::npos) {
    out.pdf1_ptr = &drdm_scalar;
    out.cdf1_ptr = &prdm_scalar;
    out.model_dfun_raw = &drdm_raw;
    out.model_pfun_raw = &prdm_raw;
    out.logS_at_t_ptr = &rdm_logS_at_t;
    out.ctx.t0_index = 3;
  } else if (type_std.find("LNR") != std::string::npos) {
    out.pdf1_ptr = &dlnr_scalar;
    out.cdf1_ptr = &plnr_scalar;
    out.model_dfun_raw = &dlnr_raw;
    out.model_pfun_raw = &plnr_raw;
    out.logS_at_t_ptr = &lnr_logS_at_t;
    out.ctx.t0_index = 2;
  } else {
    Rcpp::stop("Unsupported race model type string in %s: %s", caller.c_str(), type_std.c_str());
  }

  if (type_std.find("GNG") != std::string::npos) {
    out.ctx.gng = true;
  }
  return out;
}

// Pre-computed per-data state for likelihood functions (Race, DDM, etc.).
// Built once outside the particle loop; reused across particles to eliminate
// per-particle R-heap allocations and repeated attribute/column reads.
struct ModelSharedState {
  bool valid = false;
  // Pre-read censoring/truncation bounds (Rcpp vectors keep memory alive)
  Rcpp::NumericVector LT_vec, UT_vec, LC_vec, UC_vec;
  // Pre-computed finite/other trial partition
  Rcpp::LogicalVector finite_mask;       // length n_trials; shared ref from dadm attr
  std::vector<int>    finite_mask_int;   // 0/1 representation (for DDM/SIMD)
  std::vector<int>    finite_unique_idx; // indices of finite unique trials
  std::vector<int>    other_unique_idx;  // indices of other (non-finite) unique trials
  // Pre-allocated mutable scratch buffers (length n_trials)
  std::vector<double> res_buf;   // log-density or result; NOT re-initialised between particles
  std::vector<int>    idx_win;   // data-fixed winner mask  (finite rows only)
  std::vector<int>    idx_loss;  // data-fixed loser mask   (finite rows only)
  std::vector<int>    ok_int_buf;  // per-particle validity   (re-filled each call)
  bool any_win  = false;
  bool any_loss = false;
  // pContaminant column: -2=not yet searched, -1=absent, >=0=column index
  int  pc_col   = -2;

  // DDM-specific data
  std::vector<double> logF_LT_1, logF_LT_2, logF_UT_1, logF_UT_2;
  // Pre-allocated scratch buffers for nonfinite/trunc path; avoids per-particle R heap.
  std::vector<double> lF_LC_1_buf, lF_LC_2_buf, lF_UC_1_buf, lF_UC_2_buf;
  std::vector<int>    R1_int_buf, R2_int_buf;  // constant all-1 / all-2 response vectors
  std::vector<int>    all_ones_int_buf;        // constant all-1 mask
  bool any_ok_finite = false;
  bool any_ok_nonfinite = false;
  SEXP shared_R_levels = R_NilValue;     // Cached response levels for Go/No-go logic
};

struct RaceSharedState : ModelSharedState {};
using DDMSharedState = ModelSharedState;

struct LogicalRulesSharedState {
  bool valid = false;
  int n_trials = 0;
  int n_acc = 0;
  int n_unique_trials = 0;
  std::vector<int> idxA;
  std::vector<int> idxB;
  std::vector<int> idxnA;
  std::vector<int> idxnB;
  // Rule code: 1=OR, 2=AND, 3=XOR, 4=ID
  std::vector<int> rule_code;
  // Response code: 1=yes, 2=no, 3=NN, 4=AN, 5=NB, 6=AB
  std::vector<int> resp_code;
  std::vector<double> rt_unique;
  // Per-row RTs (length n_trials) for raw dfun/pfun kernels.
  std::vector<double> rt_by_row;
};

struct RedundantTargetSharedState {
  bool valid = false;
  int n_trials = 0;
  int n_acc = 0;
  int n_unique_trials = 0;
  bool has_nogo = false;
  std::vector<int> idxA;
  std::vector<int> idxB;
  std::vector<int> idxNogo;
  std::vector<int> cond_code;  // 1=A, 2=B, 3=AB
  // Response code per unique trial: 0=missing/omission, 1=go response, 2=nogo.
  std::vector<int> resp_code;
  std::vector<double> rt_unique;
  // Per-unique-trial bounds for go/no-go omission handling.
  std::vector<double> LT_unique;
  std::vector<double> UC_unique;
  // Per-row RTs (length n_trials) for raw dfun/pfun kernels.
  std::vector<double> rt_by_row;
};

static LogicalRulesSharedState build_logicalrules_shared_state(const Rcpp::DataFrame& dadm,
                                                               int n_trials,
                                                               int n_acc);
static RedundantTargetSharedState build_redundant_target_shared_state(const Rcpp::DataFrame& dadm,
                                                                      int n_trials,
                                                                      int n_acc);

static double local_race_helper(double t,
                                const double* par_target,
                                const double* par_nontarget,
                                int n_par,
                                ContextForRaceModels* ctx,
                                RacePdf1Fun pdf1,
                                RaceCdf1Fun cdf1,
                                const GslIntegrationControls& gsl_ctl,
                                gsl_integration_workspace* w,
                                double* pars_2buf,
                                int* isok_2buf);

static double c_log_likelihood_logicalrules(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& ok_params,
    int n_acc,
    ContextForRaceModels* model_ctx,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    const LogicalRulesSharedState& shared,
    Rcpp::NumericVector* trial_ll_out = nullptr);

static double c_log_likelihood_redundant_target_race(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& ok_params,
    int n_acc,
    ContextForRaceModels* model_ctx,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    const RedundantTargetSharedState& shared,
    Rcpp::NumericVector* trial_ll_out = nullptr);

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
    // Use R_FINITE (not std::isfinite) — -ffast-math breaks std::isfinite for Inf values
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
  double sum_ll = 0.0;
  for (int i = 0; i < n_out; ++i) {
    sum_ll += lls[expand[i] - 1]; // expand created in 1-based R
  }
  return sum_ll;
}

// Raw-buffer variant for DDM to skip materialization and allocations.
// Handles truncation and censoring with high numerical stability.
double c_log_likelihood_DDM_pt(const double* pars_cm,
                               const double* rt_ptr,
                               const int* R_ptr,
                               const int n_trials,
                               const int* expand_ptr,
                               const int n_out,
                               double min_ll,
                               const int* is_ok,
                               bool gng,
                               bool all_finite_untruncated,
                               const std::vector<int>& p_idx,
                               ModelSharedState* shared,
                               Rcpp::NumericVector* trial_ll_out = nullptr) {

  // 1. Fast Path: All RTs finite, no truncation, no censoring
  if (all_finite_untruncated) {
    if (shared->all_ones_int_buf.size() != static_cast<size_t>(n_trials)) {
      shared->all_ones_int_buf.assign(n_trials, 1);
    }
    d_DDM_Wien_raw(rt_ptr, R_ptr, pars_cm, n_trials, (int)p_idx.size(),
                   shared->all_ones_int_buf.data(), is_ok, shared->res_buf.data(), min_ll, p_idx);
    
    const double* lls_ptr = shared->res_buf.data();
    double total_ll = 0.0;
    if (trial_ll_out != nullptr) {
      for (int i = 0; i < n_out; ++i) {
        const int idx = expand_ptr[i] - 1;
        double v = lls_ptr[idx];
        if (!R_FINITE(v) || v < min_ll) v = min_ll;
        (*trial_ll_out)[i] = v;
        total_ll += v;
      }
    } else {
      #pragma omp simd reduction(+:total_ll)
      for (int i = 0; i < n_out; ++i) {
        const int idx = expand_ptr[i] - 1;
        double v = lls_ptr[idx];
        if (!R_FINITE(v) || v < min_ll) v = min_ll;
        total_ll += v;
      }
    }
    return total_ll;
  }

  // 2. Comprehensive Path: Handles truncation, censoring, and non-finite RTs
  const double* LT = shared->LT_vec.begin();
  const double* UT = shared->UT_vec.begin();
  
  std::fill(shared->res_buf.begin(), shared->res_buf.end(), min_ll);
  shared->any_ok_finite = false;
  shared->any_ok_nonfinite = false;
  
  for (int i = 0; i < n_trials; ++i) {
    if (R_FINITE(rt_ptr[i])) {
      shared->finite_mask_int[static_cast<size_t>(i)] = 1;
      if (is_ok[i]) shared->any_ok_finite = true;
    } else {
      shared->finite_mask_int[static_cast<size_t>(i)] = 0;
      if (is_ok[i]) shared->any_ok_nonfinite = true;
    }
  }

  // Helper for stable log-difference of CDFs
  auto clamp_log_cdf = [&](double logF) -> double {
    if (!R_FINITE(logF)) return R_NegInf;
    if (logF > 0.0) return 0.0;
    return logF;
  };
  auto log_sum_cdf = [&](double logF_a, double logF_b) -> double {
    return clamp_log_cdf(log_sum_exp(clamp_log_cdf(logF_a), clamp_log_cdf(logF_b)));
  };
  auto log_interval_mass = [&](double logF_low, double logF_high) -> double {
    const double lo = clamp_log_cdf(logF_low);
    const double hi = clamp_log_cdf(logF_high);
    if (hi == R_NegInf) return R_NegInf;
    if (lo == R_NegInf) return hi;
    // Allow slight numerical noise where lo > hi due to integration error
    if (hi <= lo + 1e-14) return R_NegInf; 
    return log_diff_exp(hi, lo);
  };

  // Pre-calculate truncation normalizers (logZ) for all trials
  std::vector<double> logZ(static_cast<size_t>(n_trials), 0.0);
  if (!gng) {
    // These buffers are also consumed by the censored (non-finite RT) branch
    // below. Even when truncation is inactive, keep them sized with stable
    // defaults corresponding to LT=0 and UT=Inf:
    //   logF(LT=0) = -Inf, logF(UT=Inf) = 0
    if (shared->logF_LT_1.size() != static_cast<size_t>(n_trials)) shared->logF_LT_1.resize(n_trials);
    if (shared->logF_LT_2.size() != static_cast<size_t>(n_trials)) shared->logF_LT_2.resize(n_trials);
    if (shared->logF_UT_1.size() != static_cast<size_t>(n_trials)) shared->logF_UT_1.resize(n_trials);
    if (shared->logF_UT_2.size() != static_cast<size_t>(n_trials)) shared->logF_UT_2.resize(n_trials);
    std::fill(shared->logF_LT_1.begin(), shared->logF_LT_1.end(), R_NegInf);
    std::fill(shared->logF_LT_2.begin(), shared->logF_LT_2.end(), R_NegInf);
    std::fill(shared->logF_UT_1.begin(), shared->logF_UT_1.end(), 0.0);
    std::fill(shared->logF_UT_2.begin(), shared->logF_UT_2.end(), 0.0);

    bool any_LT = false, any_UT_finite = false;
    for (int i = 0; i < n_trials; ++i) {
      if (is_ok[i]) {
        if (LT[i] != 0.0) any_LT = true;
        if (R_FINITE(UT[i])) any_UT_finite = true;
      }
      if (any_LT && any_UT_finite) break;
    }
    const bool any_trunc = any_LT || any_UT_finite;
    if (any_trunc) {
      // Use pre-allocated buffers from shared state (avoids per-particle R heap)
      const int* R1_ptr = shared->R1_int_buf.empty()
          ? nullptr : shared->R1_int_buf.data();
      const int* R2_ptr = shared->R2_int_buf.empty()
          ? nullptr : shared->R2_int_buf.data();
      const int* ones_ptr = shared->all_ones_int_buf.empty()
          ? nullptr : shared->all_ones_int_buf.data();
      // Fallback local buffers for the c_log_likelihood_DDM (non-pt) compatibility path
      Rcpp::IntegerVector R1_local, R2_local;
      std::vector<int> all_ones_local;
      if (!R1_ptr || !R2_ptr || !ones_ptr) {
        R1_local = Rcpp::IntegerVector(n_trials, 1);
        R2_local = Rcpp::IntegerVector(n_trials, 2);
        all_ones_local.assign(n_trials, 1);
        R1_ptr  = R1_local.begin();
        R2_ptr  = R2_local.begin();
        ones_ptr = all_ones_local.data();
      }
      // Only call p_DDM at LT when some trial has LT > 0 — avoids p_DDM(0,...) calls.
      // Only call p_DDM at UT when some trial has finite UT — avoids p_DDM(Inf,...) calls
      // which trigger expensive/degenerate Wiener evaluations at t=Inf.
      if (any_LT) {
        p_DDM_Wien_raw(shared->LT_vec.begin(), R1_ptr, pars_cm, n_trials, (int)p_idx.size(),
                       ones_ptr, is_ok, shared->logF_LT_1.data(), R_NegInf, p_idx);
        p_DDM_Wien_raw(shared->LT_vec.begin(), R2_ptr, pars_cm, n_trials, (int)p_idx.size(),
                       ones_ptr, is_ok, shared->logF_LT_2.data(), R_NegInf, p_idx);
      }
      if (any_UT_finite) {
        p_DDM_Wien_raw(shared->UT_vec.begin(), R1_ptr, pars_cm, n_trials, (int)p_idx.size(),
                       ones_ptr, is_ok, shared->logF_UT_1.data(), R_NegInf, p_idx);
        p_DDM_Wien_raw(shared->UT_vec.begin(), R2_ptr, pars_cm, n_trials, (int)p_idx.size(),
                       ones_ptr, is_ok, shared->logF_UT_2.data(), R_NegInf, p_idx);
      }

      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i]) continue;
        double hi_val = R_FINITE(UT[i]) ? log_sum_cdf(shared->logF_UT_1[i], shared->logF_UT_2[i]) : 0.0;
        double lo_val = (LT[i] == 0.0) ? R_NegInf : log_sum_cdf(shared->logF_LT_1[i], shared->logF_LT_2[i]);
        logZ[i] = log_interval_mass(lo_val, hi_val);
      }
    }
  }

  // Calculate finite RT densities
  if (shared->any_ok_finite) {
  d_DDM_Wien_raw(rt_ptr, R_ptr, pars_cm, n_trials, (int)p_idx.size(),
                 shared->finite_mask_int.data(), is_ok, shared->res_buf.data(), min_ll, p_idx);
  if (!gng) {

      for (int i = 0; i < n_trials; ++i) {
        if (is_ok[i] && shared->finite_mask_int[i]) {
          if (R_FINITE(logZ[i])) {
            if (shared->res_buf[i] > min_ll) shared->res_buf[i] -= logZ[i];
          }
          else shared->res_buf[i] = min_ll;
        }
      }
    }
  }

  // Calculate non-finite RT masses (censored)
    if (shared->any_ok_nonfinite) {
    if (gng) {
      // Go/No-go DDM (untruncated)
      if (shared->all_ones_int_buf.size() != static_cast<size_t>(n_trials)) {
        shared->all_ones_int_buf.assign(n_trials, 1);
      }
      if (shared->R1_int_buf.size() != static_cast<size_t>(n_trials)) {
        shared->R1_int_buf.assign(n_trials, 1);
      }
      int* R_go_ptr = shared->R1_int_buf.data();
      SEXP lev_sexp = shared->shared_R_levels;
      if (lev_sexp != R_NilValue) {
        Rcpp::CharacterVector levs(lev_sexp);
        int nogo_idx = -1;
        for (int j = 0; j < levs.size(); ++j) {
          if (Rcpp::as<std::string>(levs[j]) == "nogo") { nogo_idx = j + 1; break; }
        }
        if (nogo_idx != -1) {
          int go_idx = 1;
          for (int j = 0; j < levs.size(); ++j) {
            if (j + 1 != nogo_idx) { go_idx = j + 1; break; }
          }
          std::fill(shared->R1_int_buf.begin(), shared->R1_int_buf.end(), go_idx);
          R_go_ptr = shared->R1_int_buf.data();
        }
      }
      Rcpp::NumericVector logcdf_U(n_trials), logcdf_L(n_trials);
      p_DDM_Wien_raw(shared->UC_vec.begin(), R_go_ptr, pars_cm, n_trials, (int)p_idx.size(),
                     shared->all_ones_int_buf.data(), is_ok, logcdf_U.begin(), R_NegInf, p_idx);
      p_DDM_Wien_raw(shared->LC_vec.begin(), R_go_ptr, pars_cm, n_trials, (int)p_idx.size(),
                     shared->all_ones_int_buf.data(), is_ok, logcdf_L.begin(), R_NegInf, p_idx);
      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i] || shared->finite_mask_int[i]) continue;
        if (rt_ptr[i] == R_PosInf) shared->res_buf[i] = log1m_exp(logcdf_U[i]);
        else if (rt_ptr[i] == R_NegInf) shared->res_buf[i] = logcdf_L[i];
        else shared->res_buf[i] = log_sum_exp(logcdf_L[i], log1m_exp(logcdf_U[i]));
      }
    } else {
      // Standard DDM censored trials — use pre-allocated shared buffers to avoid
      // per-particle R heap allocation (critical for LT path where n_nonfinite > 0).
      const bool have_shared_bufs = (shared != nullptr &&
                                     !shared->R1_int_buf.empty() &&
                                     (int)shared->R1_int_buf.size() >= n_trials);
      // Fallback allocations (only used when called outside calc_ll_oo, e.g. standalone)
      Rcpp::IntegerVector R1_fb, R2_fb;
      Rcpp::NumericVector lF_LC_1_fb, lF_LC_2_fb, lF_UC_1_fb, lF_UC_2_fb;
      std::vector<int> all_ones_fb;
      int*    R1_ptr;    int*    R2_ptr;    int*    ones_ptr;
      double* lF_LC_1;  double* lF_LC_2;
      double* lF_UC_1;  double* lF_UC_2;
      if (have_shared_bufs) {
        R1_ptr   = shared->R1_int_buf.data();
        R2_ptr   = shared->R2_int_buf.data();
        ones_ptr = shared->all_ones_int_buf.data();
        lF_LC_1  = shared->lF_LC_1_buf.data();
        lF_LC_2  = shared->lF_LC_2_buf.data();
        lF_UC_1  = shared->lF_UC_1_buf.data();
        lF_UC_2  = shared->lF_UC_2_buf.data();
      } else {
        R1_fb = Rcpp::IntegerVector(n_trials, 1);
        R2_fb = Rcpp::IntegerVector(n_trials, 2);
        lF_LC_1_fb = Rcpp::NumericVector(n_trials);
        lF_LC_2_fb = Rcpp::NumericVector(n_trials);
        lF_UC_1_fb = Rcpp::NumericVector(n_trials);
        lF_UC_2_fb = Rcpp::NumericVector(n_trials);
        all_ones_fb.assign(n_trials, 1);
        R1_ptr   = R1_fb.begin();
        R2_ptr   = R2_fb.begin();
        ones_ptr = all_ones_fb.data();
        lF_LC_1  = lF_LC_1_fb.begin();
        lF_LC_2  = lF_LC_2_fb.begin();
        lF_UC_1  = lF_UC_1_fb.begin();
        lF_UC_2  = lF_UC_2_fb.begin();
      }
      p_DDM_Wien_raw(shared->LC_vec.begin(), R1_ptr, pars_cm, n_trials, (int)p_idx.size(),
                     ones_ptr, is_ok, lF_LC_1, R_NegInf, p_idx);
      p_DDM_Wien_raw(shared->LC_vec.begin(), R2_ptr, pars_cm, n_trials, (int)p_idx.size(),
                     ones_ptr, is_ok, lF_LC_2, R_NegInf, p_idx);
      p_DDM_Wien_raw(shared->UC_vec.begin(), R1_ptr, pars_cm, n_trials, (int)p_idx.size(),
                     ones_ptr, is_ok, lF_UC_1, R_NegInf, p_idx);
      p_DDM_Wien_raw(shared->UC_vec.begin(), R2_ptr, pars_cm, n_trials, (int)p_idx.size(),
                     ones_ptr, is_ok, lF_UC_2, R_NegInf, p_idx);

      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i] || shared->finite_mask_int[i]) continue;
        const bool r_known = (R_ptr[i] != NA_INTEGER);
        const int r_idx = r_known ? R_ptr[i] : 0;
        double mass = R_NegInf;
        if (rt_ptr[i] == R_NegInf) {
          if (!r_known) mass = log_sum_cdf(log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]), log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]));
          else if (r_idx == 1) mass = log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]);
          else mass = log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]);
        } else if (rt_ptr[i] == R_PosInf) {
          // P(UC < RT < UT).
          // When UT is finite, logF_UT_r = log F_r(UT) (computed above), and the
          // interval formula is correct: (F_1(UT)-F_1(UC)) + (F_2(UT)-F_2(UC)).
          // When UT=Inf, logF_UT_r was initialised to 0 (representing F_r(Inf)=1),
          // but the DDM CDF is DEFECTIVE: F_1(Inf)+F_2(Inf)=1, not F_r(Inf)=1
          // individually.  Using logF_UT_r=0 gives (1-F_1)+(1-F_2)=2-F_total >= 1,
          // which log_sum_cdf clamps to 0 — every censored trial contributes mass=1
          // regardless of parameters.  The correct mass is 1-F_total(UC).
          if (!r_known) {
            if (R_FINITE(UT[i])) {
              mass = log_sum_cdf(log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i]),
                                 log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]));
            } else {
              // P(RT > UC) = 1 - F_1(UC) - F_2(UC)
              mass = log_diff_exp(0.0, log_sum_exp(lF_UC_1[i], lF_UC_2[i]));
            }
          } else if (r_idx == 1) {
            mass = log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i]);
          } else {
            mass = log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]);
          }
        } else { // NA
          double m1 = (!r_known) ? log_sum_cdf(log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]), log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]))
                                 : (r_idx == 1) ? log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]) : log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]);
          // Same defective-CDF fix as the R_PosInf branch above.
          double m2;
          if (!r_known) {
            m2 = R_FINITE(UT[i])
              ? log_sum_cdf(log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i]),
                            log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]))
              : log_diff_exp(0.0, log_sum_exp(lF_UC_1[i], lF_UC_2[i]));
          } else {
            m2 = (r_idx == 1) ? log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i])
                               : log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]);
          }
          mass = log_sum_exp(m1, m2);
        }
        if (R_FINITE(logZ[i])) shared->res_buf[i] = mass - logZ[i];
        else shared->res_buf[i] = min_ll;
      }
    }
  }

  // 3. Accumulate results
  const double* res_ptr = shared->res_buf.data();
  double total_ll = 0.0;
  if (trial_ll_out != nullptr) {
    for (int i = 0; i < n_out; ++i) {
      double v = res_ptr[expand_ptr[i] - 1];
      if (!R_FINITE(v) || v < min_ll) v = min_ll;
      (*trial_ll_out)[i] = v;
      total_ll += v;
    }
  } else {
    #pragma omp simd reduction(+:total_ll)
    for (int i = 0; i < n_out; ++i) {
      double v = res_ptr[expand_ptr[i] - 1];
      if (!R_FINITE(v) || v < min_ll) v = min_ll;
      total_ll += v;
    }
  }
  return total_ll;
}

// Compatibility wrapper for old calc_ll path
double c_log_likelihood_DDM(Rcpp::NumericMatrix pars, Rcpp::DataFrame data,
                            const int n_trials, Rcpp::IntegerVector expand,
                            double min_ll, Rcpp::LogicalVector is_ok, bool gng,
                            bool all_finite_untruncated = false,
                            Rcpp::NumericVector* trial_ll_out = nullptr) {
  ModelSharedState shared;
  shared.LT_vec = get_col_with_default(data, "LT", 0.0);
  shared.UT_vec = get_col_with_default(data, "UT", R_PosInf);
  shared.LC_vec = get_col_with_default(data, "LC", 0.0);
  shared.UC_vec = get_col_with_default(data, "UC", R_PosInf);
  shared.finite_mask_int.assign(static_cast<size_t>(n_trials), 0);
  shared.res_buf.assign(static_cast<size_t>(n_trials), min_ll);
  shared.ok_int_buf.assign(static_cast<size_t>(n_trials), 0);
  shared.all_ones_int_buf.assign(static_cast<size_t>(n_trials), 1);
  shared.R1_int_buf.assign(static_cast<size_t>(n_trials), 1);
  shared.R2_int_buf.assign(static_cast<size_t>(n_trials), 2);
  Rcpp::IntegerVector R_col = data["R"];
  shared.shared_R_levels = R_col.attr("levels");
  shared.valid = true;

  Rcpp::NumericVector rts = data["rt"];
  Rcpp::IntegerVector R = data["R"];
  const int* expand_ptr = expand.begin();
  const int n_out = expand.length();
  
  std::vector<int> ok_int(n_trials);
  for(int i=0; i<n_trials; ++i) ok_int[i] = is_ok[i] ? 1 : 0;

  // Use a fixed identity mapping for the old path
  std::vector<int> p_idx = {0, 1, 2, 3, 4, 5, 6, 7}; 

  return c_log_likelihood_DDM_pt(pars.begin(), rts.begin(), R.begin(),
                                 n_trials, expand_ptr, n_out, min_ll,
                                 ok_int.data(), gng, all_finite_untruncated,
                                 p_idx, &shared, trial_ll_out);
}

// ---------------------------------------------------------------------------
// Shared helpers used by calc_ll, calc_ll_oo, and calc_ll_oo_pw
// ---------------------------------------------------------------------------

// Read the emc2_all_finite_trials data attribute; fall back to a direct scan
// when the attribute is absent (e.g. direct C++ callers that bypass R-side cache).
static inline bool read_all_finite_trials_attr(const Rcpp::DataFrame& data,
                                               int n_trials, int n_lR) {
  if (data.hasAttribute("emc2_all_finite_trials")) {
    Rcpp::LogicalVector v = data.attr("emc2_all_finite_trials");
    if (v.size() == 1 && v[0] != NA_LOGICAL) return (bool)v[0];
  }
  if (n_trials > 0 && n_lR > 0 && (n_trials % n_lR) == 0) {
    return race_data_all_finite_untruncated(data, n_trials, n_lR);
  }
  return false;
}

// Resolve function pointers and fixed column indices for SS (EXG or RDEX).
struct SSModelAdapter {
  ss_go_pdf_fn    go_lpdf_ptr;
  ss_go_pdf_fn    go_lccdf_ptr;
  ss_stop_surv_fn stop_logsurv_ptr;
  ss_stop_success_fn stop_success_ptr;
  int idx_tf;
  int idx_gf;
};

static inline SSModelAdapter resolve_ss_adapter(const std::string& type_std) {
  SSModelAdapter a;
  const bool is_exg = type_std.find("EXG") != std::string::npos;
  a.go_lpdf_ptr      = is_exg ? texg_go_lpdf          : rdex_go_lpdf;
  a.go_lccdf_ptr     = is_exg ? texg_go_lccdf          : rdex_go_lccdf;
  a.stop_logsurv_ptr = is_exg ? stop_logsurv_texg_fn   : stop_logsurv_rdex_fn;
  a.stop_success_ptr = is_exg ? ss_texg_stop_success_lpdf : ss_rdex_stop_success_lpdf;
  a.idx_tf           = is_exg ? 6 : 8;
  a.idx_gf           = is_exg ? 7 : 9;
  return a;
}

static LogicalRulesSharedState build_logicalrules_shared_state(const Rcpp::DataFrame& dadm,
                                                               int n_trials,
                                                               int n_acc) {
  LogicalRulesSharedState out;
  out.n_trials = n_trials;
  out.n_acc = n_acc;

  if (n_trials <= 0 || n_acc <= 0 || (n_trials % n_acc) != 0) {
    Rcpp::stop("LogicalRules likelihood expects n_trials > 0 and n_trials %% n_acc == 0.");
  }
  if (!dadm.containsElementNamed("rt") ||
      !dadm.containsElementNamed("lR") ||
      !dadm.containsElementNamed("R") ||
      !dadm.containsElementNamed("LogicalRule")) {
    Rcpp::stop("LogicalRules likelihood requires dadm columns: rt, lR, R, LogicalRule.");
  }

  Rcpp::NumericVector rts = dadm["rt"];
  SEXP role_sexp = dadm["lR"];
  SEXP resp_sexp = dadm["R"];
  SEXP rule_sexp = dadm["LogicalRule"];
  if (rts.size() != n_trials ||
      Rf_length(role_sexp) != n_trials ||
      Rf_length(resp_sexp) != n_trials ||
      Rf_length(rule_sexp) != n_trials) {
    Rcpp::stop("LogicalRules likelihood: dadm column lengths must match n_trials.");
  }

  const bool role_is_factor = Rf_inherits(role_sexp, "factor");
  const bool resp_is_factor = Rf_inherits(resp_sexp, "factor");
  const bool rule_is_factor = Rf_inherits(rule_sexp, "factor");

  Rcpp::IntegerVector role_code;
  Rcpp::IntegerVector resp_code;
  Rcpp::IntegerVector rule_code;
  Rcpp::CharacterVector role_chr;
  Rcpp::CharacterVector resp_chr;
  Rcpp::CharacterVector rule_chr;

  auto level_code = [](const Rcpp::CharacterVector& levels, const std::string& target) -> int {
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == target) return i + 1;
    }
    return -1;
  };

  int codeA = -1, codeB = -1, codenA = -1, codenB = -1;
  int codeYes = -1, codeNo = -1;
  int codeNN = -1, codeAN = -1, codeNB = -1, codeAB = -1;
  int codeAND = -1, codeOR = -1, codeXOR = -1, codeID = -1;

  if (role_is_factor) {
    role_code = Rcpp::IntegerVector(role_sexp);
    Rcpp::CharacterVector lev = role_code.attr("levels");
    codeA = level_code(lev, "A");
    codeB = level_code(lev, "B");
    codenA = level_code(lev, "n_A");
    codenB = level_code(lev, "n_B");
    if (codeA < 0 || codeB < 0 || codenA < 0 || codenB < 0) {
      Rcpp::stop("LogicalRules likelihood: lR levels must include A, B, n_A, n_B.");
    }
  } else {
    role_chr = Rcpp::CharacterVector(role_sexp);
  }

  if (resp_is_factor) {
    resp_code = Rcpp::IntegerVector(resp_sexp);
    Rcpp::CharacterVector lev = resp_code.attr("levels");
    codeYes = level_code(lev, "yes");
    codeNo = level_code(lev, "no");
    codeNN = level_code(lev, "NN");
    codeAN = level_code(lev, "AN");
    codeNB = level_code(lev, "NB");
    codeAB = level_code(lev, "AB");
    if (codeYes < 0 && codeNo < 0 && codeNN < 0 && codeAN < 0 && codeNB < 0 && codeAB < 0) {
      Rcpp::stop("LogicalRules likelihood: response levels must include binary (yes/no) and/or ID (NN/AN/NB/AB) codes.");
    }
  } else {
    resp_chr = Rcpp::CharacterVector(resp_sexp);
  }

  if (rule_is_factor) {
    rule_code = Rcpp::IntegerVector(rule_sexp);
    Rcpp::CharacterVector lev = rule_code.attr("levels");
    codeAND = level_code(lev, "AND");
    codeOR = level_code(lev, "OR");
    codeXOR = level_code(lev, "XOR");
    codeID = level_code(lev, "ID");
    if (codeAND < 0 && codeOR < 0 && codeXOR < 0 && codeID < 0) {
      Rcpp::stop("LogicalRules likelihood: LogicalRule levels must include at least one of AND, OR, XOR, ID.");
    }
  } else {
    rule_chr = Rcpp::CharacterVector(rule_sexp);
  }

  out.n_unique_trials = n_trials / n_acc;
  out.idxA.assign(out.n_unique_trials, -1);
  out.idxB.assign(out.n_unique_trials, -1);
  out.idxnA.assign(out.n_unique_trials, -1);
  out.idxnB.assign(out.n_unique_trials, -1);
  out.rule_code.assign(out.n_unique_trials, 0);
  out.resp_code.assign(out.n_unique_trials, 0);
  out.rt_unique.assign(out.n_unique_trials, NA_REAL);
  out.rt_by_row.assign(rts.begin(), rts.end());

  for (int j = 0; j < out.n_unique_trials; ++j) {
    const int start = j * n_acc;
    out.rt_unique[static_cast<size_t>(j)] = rts[start];

    int rcode = 0;
    if (rule_is_factor) {
      const int rc = rule_code[start];
      if (rc == codeOR) rcode = 1;
      else if (rc == codeAND) rcode = 2;
      else if (rc == codeXOR) rcode = 3;
      else if (rc == codeID) rcode = 4;
    } else {
      const std::string rv = Rcpp::as<std::string>(rule_chr[start]);
      if (rv == "OR") rcode = 1;
      else if (rv == "AND") rcode = 2;
      else if (rv == "XOR") rcode = 3;
      else if (rv == "ID") rcode = 4;
    }
    if (rcode == 0) {
      Rcpp::stop("LogicalRules likelihood: each unique trial must have rule AND, OR, XOR, or ID.");
    }
    out.rule_code[static_cast<size_t>(j)] = rcode;

    int resp = 0;
    if (resp_is_factor) {
      const int rc = resp_code[start];
      if (rc == codeYes) resp = 1;
      else if (rc == codeNo) resp = 2;
      else if (rc == codeNN) resp = 3;
      else if (rc == codeAN) resp = 4;
      else if (rc == codeNB) resp = 5;
      else if (rc == codeAB) resp = 6;
    } else {
      const std::string rv = Rcpp::as<std::string>(resp_chr[start]);
      if (rv == "yes") resp = 1;
      else if (rv == "no") resp = 2;
      else if (rv == "NN") resp = 3;
      else if (rv == "AN") resp = 4;
      else if (rv == "NB") resp = 5;
      else if (rv == "AB") resp = 6;
    }
    if (resp == 0) {
      Rcpp::stop("LogicalRules likelihood: response must be yes/no/NN/AN/NB/AB on each unique trial.");
    }
    out.resp_code[static_cast<size_t>(j)] = resp;

    for (int k = 0; k < n_acc; ++k) {
      const int idx = start + k;
      if (role_is_factor) {
        const int rc = role_code[idx];
        if (rc == codeA) out.idxA[static_cast<size_t>(j)] = idx;
        else if (rc == codeB) out.idxB[static_cast<size_t>(j)] = idx;
        else if (rc == codenA) out.idxnA[static_cast<size_t>(j)] = idx;
        else if (rc == codenB) out.idxnB[static_cast<size_t>(j)] = idx;
      } else {
        const std::string rv = Rcpp::as<std::string>(role_chr[idx]);
        if (rv == "A") out.idxA[static_cast<size_t>(j)] = idx;
        else if (rv == "B") out.idxB[static_cast<size_t>(j)] = idx;
        else if (rv == "n_A") out.idxnA[static_cast<size_t>(j)] = idx;
        else if (rv == "n_B") out.idxnB[static_cast<size_t>(j)] = idx;
      }
    }

    if (out.idxA[static_cast<size_t>(j)] < 0 ||
        out.idxB[static_cast<size_t>(j)] < 0 ||
        out.idxnA[static_cast<size_t>(j)] < 0 ||
        out.idxnB[static_cast<size_t>(j)] < 0) {
      Rcpp::stop("LogicalRules likelihood: each unique trial requires A, B, n_A, n_B accumulators.");
    }
  }

  out.valid = true;
  return out;
}

static RedundantTargetSharedState build_redundant_target_shared_state(const Rcpp::DataFrame& dadm,
                                                                      int n_trials,
                                                                      int n_acc) {
  RedundantTargetSharedState out;
  out.n_trials = n_trials;
  out.n_acc = n_acc;
  if (n_trials <= 0 || n_acc <= 0 || (n_trials % n_acc) != 0) {
    Rcpp::stop("RedundantTarget likelihood expects n_trials > 0 and n_trials %% n_acc == 0.");
  }
  if (!dadm.containsElementNamed("rt") || !dadm.containsElementNamed("lR")) {
    Rcpp::stop("RedundantTarget likelihood requires dadm columns: rt and lR.");
  }

  std::string stim_col;
  if (dadm.containsElementNamed("S")) stim_col = "S";
  else if (dadm.containsElementNamed("stimulus")) stim_col = "stimulus";
  else if (dadm.containsElementNamed("condition")) stim_col = "condition";
  else {
    Rcpp::stop("RedundantTarget likelihood requires trial-wise stimulus column `S` (or `stimulus`/`condition`).");
  }

  Rcpp::NumericVector rts = dadm["rt"];
  Rcpp::NumericVector LT = get_col_with_default(dadm, "LT", 0.0);
  Rcpp::NumericVector UC = get_col_with_default(dadm, "UC", R_PosInf);
  SEXP role_sexp = dadm["lR"];
  SEXP stim_sexp = dadm[stim_col];
  const bool has_resp_col = dadm.containsElementNamed("R");
  SEXP resp_sexp = has_resp_col ? static_cast<SEXP>(dadm["R"]) : R_NilValue;
  if (rts.size() != n_trials || LT.size() != n_trials || UC.size() != n_trials ||
      Rf_length(role_sexp) != n_trials || Rf_length(stim_sexp) != n_trials ||
      (has_resp_col && Rf_length(resp_sexp) != n_trials)) {
    Rcpp::stop("RedundantTarget likelihood: dadm column lengths must match n_trials.");
  }

  const bool role_is_factor = Rf_inherits(role_sexp, "factor");
  const bool stim_is_factor = Rf_inherits(stim_sexp, "factor");
  const bool resp_is_factor = has_resp_col && Rf_inherits(resp_sexp, "factor");
  Rcpp::IntegerVector role_code;
  Rcpp::CharacterVector role_chr;
  Rcpp::IntegerVector stim_code;
  Rcpp::CharacterVector stim_chr;
  Rcpp::IntegerVector resp_code;
  Rcpp::CharacterVector resp_chr;

  auto level_code = [](const Rcpp::CharacterVector& levels, const std::string& target) -> int {
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == target) return i + 1;
    }
    return -1;
  };

  int codeA = -1, codeB = -1, codeNogo = -1;
  int stimA = -1, stimB = -1, stimAB = -1;
  if (role_is_factor) {
    role_code = Rcpp::IntegerVector(role_sexp);
    Rcpp::CharacterVector lev = role_code.attr("levels");
    codeA = level_code(lev, "A");
    codeB = level_code(lev, "B");
    if (codeA < 0 || codeB < 0) {
      Rcpp::stop("RedundantTarget likelihood: lR levels must include A and B.");
    }
    codeNogo = level_code(lev, "nogo");
  } else {
    role_chr = Rcpp::CharacterVector(role_sexp);
  }

  if (stim_is_factor) {
    stim_code = Rcpp::IntegerVector(stim_sexp);
    Rcpp::CharacterVector lev = stim_code.attr("levels");
    stimA = level_code(lev, "A");
    stimB = level_code(lev, "B");
    stimAB = level_code(lev, "AB");
    if (stimAB < 0) {
      const int sBA = level_code(lev, "BA");
      if (sBA > 0) stimAB = sBA;
    }
    if (stimA < 0 || stimB < 0 || stimAB < 0) {
      Rcpp::stop("RedundantTarget likelihood: stimulus levels must include A, B, AB.");
    }
  } else {
    stim_chr = Rcpp::CharacterVector(stim_sexp);
  }

  int respNogo = -1;
  if (has_resp_col) {
    if (resp_is_factor) {
      resp_code = Rcpp::IntegerVector(resp_sexp);
      Rcpp::CharacterVector lev = resp_code.attr("levels");
      respNogo = level_code(lev, "nogo");
    } else {
      resp_chr = Rcpp::CharacterVector(resp_sexp);
    }
  }

  out.n_unique_trials = n_trials / n_acc;
  out.idxA.assign(out.n_unique_trials, -1);
  out.idxB.assign(out.n_unique_trials, -1);
  out.idxNogo.assign(out.n_unique_trials, -1);
  out.cond_code.assign(out.n_unique_trials, 0);
  out.resp_code.assign(out.n_unique_trials, 0);
  out.rt_unique.assign(out.n_unique_trials, NA_REAL);
  out.LT_unique.assign(out.n_unique_trials, 0.0);
  out.UC_unique.assign(out.n_unique_trials, R_PosInf);
  out.rt_by_row.assign(rts.begin(), rts.end());
  out.has_nogo = false;

  for (int j = 0; j < out.n_unique_trials; ++j) {
    const int start = j * n_acc;
    out.rt_unique[static_cast<size_t>(j)] = rts[start];
    out.LT_unique[static_cast<size_t>(j)] = LT[start];
    out.UC_unique[static_cast<size_t>(j)] = UC[start];

    int cond = 0;
    if (stim_is_factor) {
      const int sc = stim_code[start];
      if (sc == stimA) cond = 1;
      else if (sc == stimB) cond = 2;
      else if (sc == stimAB) cond = 3;
    } else {
      std::string s = Rcpp::as<std::string>(stim_chr[start]);
      if (s == "A") cond = 1;
      else if (s == "B") cond = 2;
      else if (s == "AB" || s == "BA" || s == "A+B" || s == "B+A") cond = 3;
    }
    if (cond == 0) {
      Rcpp::stop("RedundantTarget likelihood: each unique trial must have stimulus A, B, or AB.");
    }
    out.cond_code[static_cast<size_t>(j)] = cond;

    int resp = 0;
    if (has_resp_col) {
      if (resp_is_factor) {
        const int rc = resp_code[start];
        if (rc == NA_INTEGER) resp = 0;
        else if (respNogo > 0 && rc == respNogo) resp = 2;
        else resp = 1;
      } else {
        if (resp_chr[start] == NA_STRING) resp = 0;
        else {
          const std::string rv = Rcpp::as<std::string>(resp_chr[start]);
          if (rv == "nogo") resp = 2;
          else resp = 1;
        }
      }
    } else {
      resp = (R_FINITE(out.rt_unique[static_cast<size_t>(j)]) && out.rt_unique[static_cast<size_t>(j)] > 0.0) ? 1 : 0;
    }
    out.resp_code[static_cast<size_t>(j)] = resp;

    for (int k = 0; k < n_acc; ++k) {
      const int idx = start + k;
      if (role_is_factor) {
        const int rc = role_code[idx];
        if (rc == codeA) out.idxA[static_cast<size_t>(j)] = idx;
        else if (rc == codeB) out.idxB[static_cast<size_t>(j)] = idx;
        else if (codeNogo > 0 && rc == codeNogo) out.idxNogo[static_cast<size_t>(j)] = idx;
      } else {
        const std::string r = Rcpp::as<std::string>(role_chr[idx]);
        if (r == "A") out.idxA[static_cast<size_t>(j)] = idx;
        else if (r == "B") out.idxB[static_cast<size_t>(j)] = idx;
        else if (r == "nogo") out.idxNogo[static_cast<size_t>(j)] = idx;
      }
    }
    if (out.idxA[static_cast<size_t>(j)] < 0 || out.idxB[static_cast<size_t>(j)] < 0) {
      Rcpp::stop("RedundantTarget likelihood: each unique trial requires A and B accumulators.");
    }
    if (out.idxNogo[static_cast<size_t>(j)] >= 0) out.has_nogo = true;
  }

  out.valid = true;
  return out;
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
    // Canonical DDM path lives in calc_ll_oo; keep legacy entry point as a wrapper.
    return calc_ll_oo(p_matrix, data, constants, designs, type, bounds, transforms,
                      pretransforms, p_types, min_ll, trend);
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
      SSModelAdapter ssa = resolve_ss_adapter(type_std);
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
                                     ssa.go_lpdf_ptr, ssa.go_lccdf_ptr,
                                     ssa.stop_logsurv_ptr, ssa.stop_success_ptr,
                                     ssa.idx_tf, ssa.idx_gf);
      }
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    const bool is_logicalrules = (type_std.find("LogicalRules") != std::string::npos);
    const bool is_redundant_target = (type_std.find("RedundantTarget") != std::string::npos);

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);

    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    LogicalRulesSharedState logicalrules_shared;
    RedundantTargetSharedState redundant_target_shared;
    if (is_logicalrules) {
      logicalrules_shared = build_logicalrules_shared_state(data, n_trials, n_lR);
    } else if (is_redundant_target) {
      redundant_target_shared = build_redundant_target_shared_state(data, n_trials, n_lR);
    }

    const bool all_finite_trials = read_all_finite_trials_attr(data, n_trials, n_lR);
      
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
        if (is_logicalrules) {
          lls[i] = c_log_likelihood_logicalrules(pars, expand, min_ll, is_ok, n_lR,
                                                 &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                                 adapter.model_dfun_raw, adapter.model_pfun_raw,
                                                 logicalrules_shared, nullptr);
        } else if (is_redundant_target) {
          lls[i] = c_log_likelihood_redundant_target_race(pars, expand, min_ll, is_ok, n_lR,
                                                          &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                                          adapter.model_dfun_raw, adapter.model_pfun_raw,
                                                          redundant_target_shared, nullptr);
        } else {
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
  const bool use_pt_mapping = !is_mri_type;

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
    }
    keep_names = p_types;
  }

  TrendRuntime* trend_runtime_ptr = trend_runtime ? trend_runtime.get() : nullptr;

  ModelSharedState ddm_shared;
  std::vector<int> ddm_p_idx;
  bool ddm_raw_ready = true;
  if (is_ddm_type) {
      ddm_shared.LT_vec = get_col_with_default(data, "LT", 0.0);
      ddm_shared.UT_vec = get_col_with_default(data, "UT", R_PosInf);
      ddm_shared.LC_vec = get_col_with_default(data, "LC", 0.0);
      ddm_shared.UC_vec = get_col_with_default(data, "UC", R_PosInf);
      ddm_shared.finite_mask_int.resize(n_trials);
      ddm_shared.res_buf.resize(n_trials);
      ddm_shared.ok_int_buf.resize(n_trials);
      // Pre-allocate DDM scratch buffers (avoids per-particle R heap in trunc/cens paths)
      ddm_shared.lF_LC_1_buf.resize(n_trials);
      ddm_shared.lF_LC_2_buf.resize(n_trials);
      ddm_shared.lF_UC_1_buf.resize(n_trials);
      ddm_shared.lF_UC_2_buf.resize(n_trials);
      ddm_shared.R1_int_buf.assign(n_trials, 1);
      ddm_shared.R2_int_buf.assign(n_trials, 2);
      ddm_shared.all_ones_int_buf.assign(n_trials, 1);
      Rcpp::IntegerVector R_col = data["R"];
      ddm_shared.shared_R_levels = R_col.attr("levels");
      ddm_shared.valid = true;
      
      const std::vector<std::string> ddm_names = {"v", "a", "sv", "t0", "st0", "s", "Z", "SZ"};
      for(const auto& nm : ddm_names) {
        auto it = param_table_template.name_to_base_idx.find(nm);
        int idx = (it != param_table_template.name_to_base_idx.end()) ? it->second : -1;
        ddm_p_idx.push_back(idx);
        if (idx < 0) ddm_raw_ready = false;
      }
  }

  if (is_ddm_type) {
    bool gng = (type_std.find("GNG") != std::string::npos);
    IntegerVector expand = data.attr("expand");
    NumericVector rts = data["rt"];
    IntegerVector R = data["R"];
    const double* rt_ptr = rts.begin();
    const int* R_ptr = R.begin();
    const int* expand_ptr = expand.begin();
    const int n_out = expand.length();
    const bool all_finite_untruncated = ddm_data_all_finite_untruncated(data, n_trials);
    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) {
        param_table_template.fill_from_particle_row(particle_matrix_pt, i, pm_col_to_base_idx);
      }
      update_pt_only(param_table_template, designs, trend_runtime_ptr, transform_specs_pt);
      if (i == 0) {
        bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      }
      is_ok = c_do_bound_pt(param_table_template, bound_specs);
      if (ddm_raw_ready) {
        for(int j = 0; j < n_trials; ++j) ddm_shared.ok_int_buf[j] = is_ok[j] ? 1 : 0;
        lls[i] = c_log_likelihood_DDM_pt(param_table_template.base.begin(),
                                        rt_ptr, R_ptr, n_trials, expand_ptr, n_out,
                                        min_ll, ddm_shared.ok_int_buf.data(), gng,
                                        all_finite_untruncated, ddm_p_idx, &ddm_shared);
      } else {
        pars = param_table_template.materialize_by_param_names(keep_names);
        lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok,
                                      gng, all_finite_untruncated);
      }
    }
  }
 else if (is_mri_type) {
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
    SSModelAdapter ssa = resolve_ss_adapter(type_std);

    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) {
        param_table_template.fill_from_particle_row(particle_matrix_pt, i, pm_col_to_base_idx);
      }
      update_pt_only(param_table_template, designs, trend_runtime_ptr, transform_specs_pt);
      if (i == 0) {
        bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      }
      pars = param_table_template.materialize_by_param_names(keep_names);
      is_ok = c_do_bound_pt(param_table_template, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      lls[i] = c_log_likelihood_ss(pars, data, n_trials_ss, expand, min_ll, is_ok,
                                   ssa.go_lpdf_ptr, ssa.go_lccdf_ptr,
                                   ssa.stop_logsurv_ptr, ssa.stop_success_ptr,
                                   ssa.idx_tf, ssa.idx_gf);
    }
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    const bool is_logicalrules = (type_std.find("LogicalRules") != std::string::npos);
    const bool is_redundant_target = (type_std.find("RedundantTarget") != std::string::npos);

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll_oo");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);

    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();

    if (is_logicalrules) {
      LogicalRulesSharedState logicalrules_shared =
        build_logicalrules_shared_state(data, n_trials, n_lR);

      for (int i = 0; i < n_particles; ++i) {
        if (i > 0) {
          param_table_template.fill_from_particle_row(particle_matrix_pt, i, pm_col_to_base_idx);
        }
        update_pt_only(param_table_template, designs, trend_runtime_ptr, transform_specs_pt);
        if (i == 0) {
          bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
        }
        is_ok = c_do_bound_pt(param_table_template, bound_specs);
        is_ok = lr_all(is_ok, n_lR);
        pars = param_table_template.materialize_by_param_names(keep_names);
        lls[i] = c_log_likelihood_logicalrules(pars, expand, min_ll, is_ok, n_lR,
                                               &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                               adapter.model_dfun_raw, adapter.model_pfun_raw,
                                               logicalrules_shared, nullptr);
      }
      return lls;
    }
    if (is_redundant_target) {
      RedundantTargetSharedState redundant_target_shared =
        build_redundant_target_shared_state(data, n_trials, n_lR);

      for (int i = 0; i < n_particles; ++i) {
        if (i > 0) {
          param_table_template.fill_from_particle_row(particle_matrix_pt, i, pm_col_to_base_idx);
        }
        update_pt_only(param_table_template, designs, trend_runtime_ptr, transform_specs_pt);
        if (i == 0) {
          bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
        }
        is_ok = c_do_bound_pt(param_table_template, bound_specs);
        is_ok = lr_all(is_ok, n_lR);
        pars = param_table_template.materialize_by_param_names(keep_names);
        lls[i] = c_log_likelihood_redundant_target_race(pars, expand, min_ll, is_ok, n_lR,
                                                       &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                                       adapter.model_dfun_raw, adapter.model_pfun_raw,
                                                       redundant_target_shared, nullptr);
      }
      return lls;
    }
    const bool all_finite_trials = read_all_finite_trials_attr(data, n_trials, n_lR);

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
      all_finite_trials &&  // already scanned above; avoids a redundant O(n) pass
      (!has_RACE_col_fp || has_RACE_attrs_fp) &&
      adapter.model_pfun_raw != nullptr &&
      adapter.model_dfun_raw != nullptr;

    const int n_unique_fp = n_trials / n_lR;
    std::vector<double> res_buf;
    std::vector<double> ll_uniq_buf;
    std::vector<int> winner_int_buf;
    std::vector<int> loser_int_buf;
    std::vector<int> isok_int_fp;
    NumericVector rts_fp_hold;
    const double* rt_ptr = nullptr;
    const int* expand_ptr = expand.begin();
    const int n_exp = expand.length();

    // Pre-allocated staging buffer: copies base columns into p_types order so
    // raw model functions (dlba_raw etc.) see the column layout they expect.
    // Populated once per particle from param_table_template.base.
    std::vector<int>    race_base_col_order;  // base_idx for each keep_names entry
    std::vector<double> race_staging_buf;     // n_trials * keep_names.size()
    int fast_pc_staging_col = -1;             // keep_names position of pContaminant

    if (use_raw_fast_path) {
      res_buf.resize(n_trials);
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

      // Build column-order mapping: keep_names[j] -> base column index
      const int n_kn = keep_names.size();
      race_base_col_order.resize(n_kn, -1);
      for (int j = 0; j < n_kn; ++j) {
        std::string nm = Rcpp::as<std::string>(keep_names[j]);
        auto it = param_table_template.name_to_base_idx.find(nm);
        if (it != param_table_template.name_to_base_idx.end()) {
          race_base_col_order[j] = it->second;
        }
        if (nm == "pContaminant") fast_pc_staging_col = j;
      }
      race_staging_buf.resize(static_cast<size_t>(n_trials) * n_kn);
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
        race_shared.res_buf.resize(static_cast<size_t>(n_trials)); // no init needed
        race_shared.idx_win.assign(static_cast<size_t>(n_trials), 0);
        race_shared.idx_loss.assign(static_cast<size_t>(n_trials), 0);
        race_shared.ok_int_buf.resize(static_cast<size_t>(n_trials));

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
      
      update_pt_only(param_table_template, designs, trend_runtime_ptr, transform_specs_pt);
      if (i == 0) {
        bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      }

      is_ok = c_do_bound_pt(param_table_template, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      if (use_raw_fast_path) {
        // Fill per-particle isok buffer
        for (int j = 0; j < n_trials; ++j) isok_int_fp[j] = is_ok[j] ? 1 : 0;

        // Copy base columns into staging buffer in p_types order so raw model
        // functions (dlba_raw, drdm_raw, etc.) see the expected column layout.
        const int n_kn = (int)race_base_col_order.size();
        for (int col = 0; col < n_kn; ++col) {
          const int bidx = race_base_col_order[col];
          if (bidx < 0) continue;
          std::memcpy(race_staging_buf.data() + col * n_trials,
                      &param_table_template.base(0, bidx),
                      static_cast<size_t>(n_trials) * sizeof(double));
        }
        const double* pars_cm = race_staging_buf.data();
        adapter.ctx.mode_hint = 0;
        // Set once-per-particle mode hints so raw kernels can skip per-row
        // variability checks in common zero-variability cases.
        if (type_std.find("RDMSWTN") != std::string::npos) {
          const double* sv_col = pars_cm + 5 * n_trials;
          bool sv_zero = true;
          for (int j = 0; j < n_trials; ++j) {
            if (!isok_int_fp[j]) continue;
            const double sv = sv_col[j];
            if (std::isfinite(sv) && std::fabs(sv) > 1e-10) { sv_zero = false; break; }
          }
          adapter.ctx.mode_hint = sv_zero ? 1 : 2;
        }

        // Log-density for winner rows
        adapter.model_dfun_raw(rt_ptr, pars_cm, n_trials,
                               winner_int_buf.data(), isok_int_fp.data(),
                               res_buf.data(), min_ll, &adapter.ctx);

        // Log-survivor for loser rows (writes into same buffer, different slots)
        if (n_lR > 1) {
          adapter.model_pfun_raw(rt_ptr, pars_cm, n_trials,
                                 loser_int_buf.data(), isok_int_fp.data(),
                                 res_buf.data(), min_ll, &adapter.ctx);
        }

        // Sum n_lR rows per unique trial; apply pC correction (no-op when pC=0).
        for (int j = 0; j < n_unique_fp; ++j) {
          double s = 0.0;
          const int base = j * n_lR;
          const int n_lR_j = has_RACE_col_fp ? RACE_fp[base] : n_lR;
          for (int k = 0; k < n_lR_j; ++k) s += res_buf[base + k];
          s = s < min_ll ? min_ll : s;
          if (fast_pc_staging_col >= 0) {
            const double pC = pars_cm[fast_pc_staging_col * n_trials + base];
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
        pars = param_table_template.materialize_by_param_names(keep_names);
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
    const bool is_logicalrules = (type_std.find("LogicalRules") != std::string::npos);
    const bool is_redundant_target = (type_std.find("RedundantTarget") != std::string::npos);

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll_oo_pw");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);

    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    const bool all_finite_trials = read_all_finite_trials_attr(data, n_trials, n_lR);

    const int n_out_race = (expand.length() > 0) ? expand.length() : (n_trials / n_lR);
    NumericMatrix result(n_particles, n_out_race);
    LogicalRulesSharedState logicalrules_shared;
    RedundantTargetSharedState redundant_target_shared;
    if (is_logicalrules) {
      logicalrules_shared = build_logicalrules_shared_state(data, n_trials, n_lR);
    } else if (is_redundant_target) {
      redundant_target_shared = build_redundant_target_shared_state(data, n_trials, n_lR);
    }
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle(i);
      if (i == 0) bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      NumericVector row_vec(n_out_race);
      if (is_logicalrules) {
        c_log_likelihood_logicalrules(pars, expand, min_ll, is_ok, n_lR,
                                      &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                      adapter.model_dfun_raw, adapter.model_pfun_raw,
                                      logicalrules_shared, &row_vec);
      } else if (is_redundant_target) {
        c_log_likelihood_redundant_target_race(pars, expand, min_ll, is_ok, n_lR,
                                               &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                               adapter.model_dfun_raw, adapter.model_pfun_raw,
                                               redundant_target_shared, &row_vec);
      } else {
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
      }
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
  if (!(out > 0.0) || !emc2_isfinite(out)) return 0.0;
  
  for (int j = 0; j < P->n_lR; ++j) {
    if (j == w) continue;
    if (!P->isok[j]) return 0.0;
    const double* par_j = P->pars + static_cast<size_t>(j) * P->n_par;
    const double cdf = P->cdf1(t, par_j, P->ctx);
    if (!(cdf > 0.0)) continue;
    if (cdf >= 1.0 || !emc2_isfinite(cdf)) return 0.0;
    out *= (1.0 - cdf);
    if (!(out > 0.0) || !emc2_isfinite(out)) return 0.0;
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
    if (!emc2_isfinite(v) || !emc2_isfinite(sv) || sv <= 0.0) return R_NegInf;
    const double ll = pnorm_std(-v / sv, true, true); 
    if (!emc2_isfinite(ll)) return R_NegInf;
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
    if (race_ctx && race_ctx->defective_upper_tail) {
      double logS = 0.0;
      for (int k = 0; k < n_lR; ++k) {
        if (!isok_int[k]) return R_NegInf;
        const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
        double cdf_inf = cdf1(R_PosInf, par_k, ctx);
        cdf_inf = clamp_cdf01_race(cdf_inf);
        const double ll = safe_log1m_race(cdf_inf);
        if (!std::isfinite(ll)) return R_NegInf;
        logS += ll;
      }
      return logS;
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
    if (!emc2_isfinite(ll)) return R_NegInf;
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
      if (race_ctx && race_ctx->defective_upper_tail) {
          double logC = 0.0;
          for (int k = 0; k < n_lR; ++k) {
            if (!isok_int[k]) return R_NegInf;
            const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
            double cdf_inf = cdf1(R_PosInf, par_k, ctx);
            cdf_inf = clamp_cdf01_race(cdf_inf);
            const double ll = std::log(cdf_inf);
            if (!emc2_isfinite(ll)) return R_NegInf;
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
    if (!emc2_isfinite(ll)) return R_NegInf;
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
  if (!(t > 0.0) || !emc2_isfinite(t)) return R_NegInf;
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
    if (!emc2_isfinite(ll)) return R_NegInf;
    logS_k[static_cast<size_t>(k)] = ll;
    logS_all += ll;
  }
  double out = R_NegInf;
  for (int k = 0; k < n_lR; ++k) {
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    const double pdf = pdf1(t, par_k, ctx);
    if (!(pdf > 0.0) || !emc2_isfinite(pdf)) continue;
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
    if (gsl_ctl.try_qng_first_finite) {
      size_t neval = 0;
      const int qng_status = gsl_integration_qng(&F, low, upp, abs_tol, rel_tol,
                                                 &result, &error, &neval);
      if (qng_status == GSL_SUCCESS && R_FINITE(result)) return GSL_SUCCESS;
      return gsl_integration_qag(&F, low, upp, abs_tol, rel_tol, limit, gsl_ctl.qag_key,
                                 w, &result, &error);
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
  if (!(result > 0.0) || !R_FINITE(result)) return R_NegInf;
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

static inline void copy_par_row_colmajor(const double* pars_cm_ptr,
                                         int n_rows,
                                         int n_par,
                                         int row_idx,
                                         double* out_row) {
  for (int c = 0; c < n_par; ++c) {
    out_row[c] = pars_cm_ptr[static_cast<size_t>(c) * n_rows + row_idx];
  }
}

static inline bool row_equal_colmajor(const double* pars_cm_ptr,
                                      int n_rows,
                                      int n_par,
                                      int row_a,
                                      int row_b) {
  for (int c = 0; c < n_par; ++c) {
    if (pars_cm_ptr[static_cast<size_t>(c) * n_rows + row_a] !=
        pars_cm_ptr[static_cast<size_t>(c) * n_rows + row_b]) {
      return false;
    }
  }
  return true;
}

static double local_race_helper(double t,
                                const double* par_target,
                                const double* par_nontarget,
                                int n_par,
                                ContextForRaceModels* ctx,
                                RacePdf1Fun pdf1,
                                RaceCdf1Fun cdf1,
                                const GslIntegrationControls& gsl_ctl,
                                gsl_integration_workspace* w,
                                double* pars_2buf,
                                int* isok_2buf) {
  if (!(t > 0.0) || !R_FINITE(t)) return 0.0;
  if (ctx == nullptr || pdf1 == nullptr || cdf1 == nullptr || w == nullptr) return NA_REAL;

  const int t0_idx = ctx->t0_index;
  double t0_tgt = (t0_idx >= 0 && t0_idx < n_par) ? par_target[t0_idx] : 0.0;
  double t0_nt  = (t0_idx >= 0 && t0_idx < n_par) ? par_nontarget[t0_idx] : 0.0;
  if (t0_tgt < 0.0) t0_tgt = 0.0;
  if (t0_nt < 0.0) t0_nt = 0.0;

  if (t <= t0_tgt) return 0.0;

  double term1 = 0.0;
  const double split = std::min(t, t0_nt);
  if (split > t0_tgt) {
    double c_hi = cdf1(split, par_target, ctx);
    double c_lo = cdf1(t0_tgt, par_target, ctx);
    c_hi = clamp_cdf01_race(c_hi);
    c_lo = clamp_cdf01_race(c_lo);
    if (emc2_isnan(c_hi) || emc2_isnan(c_lo)) return NA_REAL;
    term1 = c_hi - c_lo;
    if (term1 < 0.0) term1 = 0.0;
    if (term1 > 1.0) term1 = 1.0;
  }

  const double low_limit = std::max(t0_tgt, t0_nt);
  if (t <= low_limit) return term1;

  std::memcpy(pars_2buf, par_target, static_cast<size_t>(n_par) * sizeof(double));
  std::memcpy(pars_2buf + n_par, par_nontarget, static_cast<size_t>(n_par) * sizeof(double));
  isok_2buf[0] = 1;
  isok_2buf[1] = 1;

  const double log_g = integrate_for_kth_winner_rowmajor_cpp(
    1,
    pars_2buf,
    isok_2buf,
    low_limit,
    t,
    pdf1,
    cdf1,
    2,
    n_par,
    gsl_ctl,
    ctx,
    w);
  if (!R_FINITE(log_g)) return NA_REAL;
  double g = term1 + std::exp(log_g);
  if (g < 0.0) g = 0.0;
  if (g > 1.0) g = 1.0;
  return g;
}

static double c_log_likelihood_logicalrules(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& ok_params,
    int n_acc,
    ContextForRaceModels* model_ctx,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    const LogicalRulesSharedState& shared,
    Rcpp::NumericVector* trial_ll_out) {
  if (!shared.valid || n_acc <= 0 || pdf1 == nullptr || cdf1 == nullptr || model_ctx == nullptr) {
    Rcpp::stop("c_log_likelihood_logicalrules: invalid logical-rules configuration.");
  }
  const int n_trials = shared.n_trials;
  if (pars.nrow() != n_trials || ok_params.size() != n_trials) {
    Rcpp::stop("c_log_likelihood_logicalrules: pars/ok_params sizes must match shared data.");
  }
  if (n_trials == 0) return 0.0;

  static const double kMinSurv = 1e-300;
  const int n_par = pars.ncol();
  const int n_unique_trials = shared.n_unique_trials;
  const double* pars_cm_ptr = pars.begin();
  std::vector<double> ll_unique(static_cast<size_t>(n_unique_trials), min_ll);
  const bool use_raw_local =
    (model_dfun_raw != nullptr &&
     model_pfun_raw != nullptr &&
     static_cast<int>(shared.rt_by_row.size()) == n_trials);

  GslIntegrationControls gsl_ctl = default_gsl_controls();
  gsl_ctl.try_qng_first_finite = true;
  gsl_ctl.qag_key = GSL_INTEG_GAUSS21;
  gsl_ctl.rel_tol = 1e-4;  // tighter than default unnecessary for log-likelihood accuracy

  // Change 3: skip the raw batch (and all integration) if every parameter row is invalid.
  bool any_ok = false;
  for (int i = 0; i < n_trials && !any_ok; ++i) any_ok = (bool)ok_params[i];
  if (!any_ok) {
    const int n_out = expand.length() > 0 ? (int)expand.length() : n_unique_trials;
    if (trial_ll_out != nullptr) for (int i = 0; i < n_out; ++i) (*trial_ll_out)[i] = min_ll;
    return static_cast<double>(n_out) * min_ll;
  }

  std::vector<double> logf_all;
  std::vector<double> logS_all;
  std::vector<int> all_mask;
  std::vector<int> isok_int_all;
  if (use_raw_local) {
    logf_all.resize(static_cast<size_t>(n_trials), min_ll);
    logS_all.resize(static_cast<size_t>(n_trials), min_ll);
    all_mask.assign(static_cast<size_t>(n_trials), 1);
    isok_int_all.resize(static_cast<size_t>(n_trials), 0);
    for (int i = 0; i < n_trials; ++i) isok_int_all[static_cast<size_t>(i)] = ok_params[i] ? 1 : 0;
    model_dfun_raw(shared.rt_by_row.data(), pars_cm_ptr, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logf_all.data(), min_ll, model_ctx);
    model_pfun_raw(shared.rt_by_row.data(), pars_cm_ptr, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logS_all.data(), min_ll, model_ctx);
  }

  // ---- Gauss-Legendre batch integration pre-pass ----
  // Replaces per-trial local_race_helper (serial GSL) with N_GL vectorised batch sweeps.
  // GA_no_gl[j] = P(nontarget-A wins channel-A sub-race before RT_j)  [linear probability]
  // GB_no_gl[j] = same for channel B, or = GA_no_gl[j] when channels are equal.
  // OR/AND rules use GA_no_gl directly; XOR/ID rules derive GA_yes = 1 - GA_no_gl - S_dec_A.
  const bool use_gl_pass = use_raw_local;
  std::vector<double> GA_no_gl, GB_no_gl;
  std::vector<bool> ch_eq_vec;

  if (use_gl_pass) {
    static constexpr int N_GL = 31;
    static std::once_flag s_gl_flag;
    static std::vector<double> s_gl_x, s_gl_w;
    std::call_once(s_gl_flag, []() {
      // Newton iteration on the Legendre polynomial P_N to find GL nodes/weights.
      // No LAPACK required; converges to machine precision in ~4 iterations per root.
      s_gl_x.resize(N_GL); s_gl_w.resize(N_GL);
      const int M = (N_GL + 1) / 2;
      for (int i = 1; i <= M; ++i) {
        double x = std::cos(M_PI * (i - 0.25) / (N_GL + 0.5));
        double p1, p2, p3, dp, dx;
        do {
          p1 = 1.0; p2 = 0.0;
          for (int j = 1; j <= N_GL; ++j) {
            p3 = p2; p2 = p1;
            p1 = ((2*j - 1) * x * p2 - (j - 1) * p3) / j;
          }
          dp = N_GL * (x * p1 - p2) / (x * x - 1.0);
          dx = p1 / dp;
          x -= dx;
        } while (std::abs(dx) > 1e-15);
        s_gl_x[i - 1]       = -x;
        s_gl_x[N_GL - i]    =  x;
        const double w = 2.0 / ((1.0 - x * x) * dp * dp);
        s_gl_w[i - 1] = s_gl_w[N_GL - i] = w;
      }
    });
    const double* gl_x = s_gl_x.data();
    const double* gl_w = s_gl_w.data();

    // Per-trial GL parameters: half-length h_j, midpoint mid_j, channels-equal flag.
    const int t0_col = model_ctx->t0_index;
    std::vector<double> gl_h(n_unique_trials), gl_mid(n_unique_trials);
    ch_eq_vec.resize(n_unique_trials);
    bool any_unequal = false;
    for (int j = 0; j < n_unique_trials; ++j) {
      const int iA  = shared.idxA[j],  inA = shared.idxnA[j];
      const int iB  = shared.idxB[j],  inB = shared.idxnB[j];
      const double t0 = (t0_col >= 0 && t0_col < n_par)
                        ? pars_cm_ptr[static_cast<size_t>(t0_col) * n_trials + iA] : 0.0;
      const double RT = shared.rt_unique[j];
      const double lo = std::max(0.0, t0);
      gl_h[j]   = (RT > lo) ? (RT - lo) * 0.5 : 0.0;
      gl_mid[j] = (RT + lo) * 0.5;
      ch_eq_vec[j] = row_equal_colmajor(pars_cm_ptr, n_trials, n_par, iA, iB)
                  && row_equal_colmajor(pars_cm_ptr, n_trials, n_par, inA, inB);
      if (!ch_eq_vec[j]) any_unequal = true;
    }

    // Compact column-major parameter matrices (n_unique_trials rows) for each accumulator role.
    // Avoids iterating over the full n_trials parameter matrix with sparse masks.
    const size_t cpt = static_cast<size_t>(n_unique_trials * n_par);
    std::vector<double> pars_nA(cpt), pars_A(cpt), pars_nB, pars_B;
    std::vector<int> isok_nA(n_unique_trials), isok_A(n_unique_trials), isok_nB, isok_B;
    std::vector<int> mask_ne;
    if (any_unequal) {
      pars_nB.resize(cpt); pars_B.resize(cpt);
      isok_nB.resize(n_unique_trials); isok_B.resize(n_unique_trials);
      mask_ne.assign(n_unique_trials, 0);
    }
    for (int p = 0; p < n_par; ++p) {
      const double* src = pars_cm_ptr + static_cast<size_t>(p) * n_trials;
      double* d_nA = pars_nA.data() + static_cast<size_t>(p) * n_unique_trials;
      double* d_A  = pars_A.data()  + static_cast<size_t>(p) * n_unique_trials;
      for (int j = 0; j < n_unique_trials; ++j) {
        d_nA[j] = src[shared.idxnA[j]];
        d_A[j]  = src[shared.idxA[j]];
      }
      if (any_unequal) {
        double* d_nB = pars_nB.data() + static_cast<size_t>(p) * n_unique_trials;
        double* d_B  = pars_B.data()  + static_cast<size_t>(p) * n_unique_trials;
        for (int j = 0; j < n_unique_trials; ++j) {
          d_nB[j] = src[shared.idxnB[j]];
          d_B[j]  = src[shared.idxB[j]];
        }
      }
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      isok_nA[j] = isok_int_all[shared.idxnA[j]];
      isok_A[j]  = isok_int_all[shared.idxA[j]];
      if (any_unequal) {
        isok_nB[j] = isok_int_all[shared.idxnB[j]];
        isok_B[j]  = isok_int_all[shared.idxB[j]];
        mask_ne[j] = ch_eq_vec[j] ? 0 : 1;
      }
    }
    std::vector<int> all1(n_unique_trials, 1);

    // GL time buffer and batch output buffers (compact, n_unique_trials rows).
    std::vector<double> gl_rt(n_unique_trials);
    std::vector<double> lf_nA(n_unique_trials, min_ll), lS_A(n_unique_trials, min_ll);
    // B buffers pre-initialised to min_ll; mask=0 rows are never written → stay at min_ll.
    std::vector<double> lf_nB(any_unequal ? n_unique_trials : 1, min_ll);
    std::vector<double> lS_B(any_unequal ? n_unique_trials : 1, min_ll);

    GA_no_gl.assign(n_unique_trials, 0.0);
    GB_no_gl.assign(n_unique_trials, 0.0);

    for (int k = 0; k < N_GL; ++k) {
      const double xi = gl_x[k], wt = gl_w[k];

      for (int j = 0; j < n_unique_trials; ++j)
        gl_rt[j] = gl_mid[j] + gl_h[j] * xi;

      // Channel A: f_nA (winner density) × S_A (loser survivor)
      model_dfun_raw(gl_rt.data(), pars_nA.data(), n_unique_trials,
                     all1.data(), isok_nA.data(), lf_nA.data(), min_ll, model_ctx);
      model_pfun_raw(gl_rt.data(), pars_A.data(),  n_unique_trials,
                     all1.data(), isok_A.data(),  lS_A.data(),  min_ll, model_ctx);
      #pragma omp simd
      for (int j = 0; j < n_unique_trials; ++j)
        GA_no_gl[j] += wt * gl_h[j] * std::exp(lf_nA[j] + lS_A[j]);

      if (any_unequal) {
        // Channel B: only for trials where B params differ from A (mask_ne).
        // ch_eq rows stay at min_ll in lf_nB/lS_B → exp contribution ≈ 0.
        model_dfun_raw(gl_rt.data(), pars_nB.data(), n_unique_trials,
                       mask_ne.data(), isok_nB.data(), lf_nB.data(), min_ll, model_ctx);
        model_pfun_raw(gl_rt.data(), pars_B.data(),  n_unique_trials,
                       mask_ne.data(), isok_B.data(),  lS_B.data(),  min_ll, model_ctx);
        #pragma omp simd
        for (int j = 0; j < n_unique_trials; ++j)
          GB_no_gl[j] += wt * gl_h[j] * std::exp(lf_nB[j] + lS_B[j]);
      }
    }
    // Channels-equal trials: B sub-race is identical to A — no separate integration needed.
    for (int j = 0; j < n_unique_trials; ++j)
      if (ch_eq_vec[j]) GB_no_gl[j] = GA_no_gl[j];
  }

  static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace_tls, gsl_ctl.retry_limit);
  std::vector<double> parA(static_cast<size_t>(n_par));
  std::vector<double> parB(static_cast<size_t>(n_par));
  std::vector<double> parnA(static_cast<size_t>(n_par));
  std::vector<double> parnB(static_cast<size_t>(n_par));
  std::vector<double> pars_2buf(static_cast<size_t>(2 * n_par));
  int isok_2buf[2] = {1, 1};

  for (int j = 0; j < n_unique_trials; ++j) {
    const int idxA = shared.idxA[static_cast<size_t>(j)];
    const int idxB = shared.idxB[static_cast<size_t>(j)];
    const int idxnA = shared.idxnA[static_cast<size_t>(j)];
    const int idxnB = shared.idxnB[static_cast<size_t>(j)];

    if (idxA < 0 || idxB < 0 || idxnA < 0 || idxnB < 0 ||
        !ok_params[idxA] || !ok_params[idxB] || !ok_params[idxnA] || !ok_params[idxnB]) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    const double t = shared.rt_unique[static_cast<size_t>(j)];
    if (!(t > 0.0) || !R_FINITE(t)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxA, parA.data());
    copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxB, parB.data());
    copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxnA, parnA.data());
    copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxnB, parnB.data());

    auto eval_pdf_cdf = [&](int idx, const double* p_row, double& f_out, double& F_out) -> bool {
      if (use_raw_local) {
        const double lf = logf_all[static_cast<size_t>(idx)];
        const double ls = logS_all[static_cast<size_t>(idx)];
        if (!R_FINITE(lf) || !R_FINITE(ls)) return false;
        const double f = (lf <= min_ll) ? 0.0 : std::exp(lf);
        double S = (ls <= min_ll) ? 0.0 : std::exp(ls);
        if (!(S >= 0.0) || !R_FINITE(S)) return false;
        if (S > 1.0) S = 1.0;
        const double F = clamp_cdf01_race(1.0 - S);
        if (!R_FINITE(f) || !R_FINITE(F) || emc2_isnan(F)) return false;
        f_out = f;
        F_out = F;
        return true;
      }
      const double f = pdf1(t, p_row, model_ctx);
      if (!R_FINITE(f) || f < 0.0) return false;
      double F = cdf1(t, p_row, model_ctx);
      F = clamp_cdf01_race(F);
      if (!R_FINITE(F) || emc2_isnan(F)) return false;
      f_out = f;
      F_out = F;
      return true;
    };

    double fA = 0.0, fB = 0.0, fnA = 0.0, fnB = 0.0;
    double FA = 0.0, FB = 0.0, FnA = 0.0, FnB = 0.0;
    if (!eval_pdf_cdf(idxA, parA.data(), fA, FA) ||
        !eval_pdf_cdf(idxB, parB.data(), fB, FB) ||
        !eval_pdf_cdf(idxnA, parnA.data(), fnA, FnA) ||
        !eval_pdf_cdf(idxnB, parnB.data(), fnB, FnB)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    const double one_m_FA = std::max(kMinSurv, 1.0 - FA);
    const double one_m_FB = std::max(kMinSurv, 1.0 - FB);
    const double one_m_FnA = std::max(kMinSurv, 1.0 - FnA);
    const double one_m_FnB = std::max(kMinSurv, 1.0 - FnB);

    const double gA_yes = fA * one_m_FnA;
    const double gB_yes = fB * one_m_FnB;
    const double gA_no = fnA * one_m_FA;
    const double gB_no = fnB * one_m_FB;

    const double S_dec_A = one_m_FA * one_m_FnA;
    const double S_dec_B = one_m_FB * one_m_FnB;

    const int rule_code = shared.rule_code[static_cast<size_t>(j)];
    const int resp_code = shared.resp_code[static_cast<size_t>(j)];
    const bool is_or_rule = (rule_code == 1);
    const bool is_and_rule = (rule_code == 2);
    const bool is_xor_rule = (rule_code == 3);
    const bool is_id_rule = (rule_code == 4);
    const bool channels_equal = use_gl_pass ? ch_eq_vec[static_cast<size_t>(j)]
        : (row_equal_colmajor(pars_cm_ptr, n_trials, n_par, idxA, idxB) &&
           row_equal_colmajor(pars_cm_ptr, n_trials, n_par, idxnA, idxnB));

    double p_j = 0.0;
    if (is_or_rule || is_and_rule) {
      double GA_no, GB_no;
      if (use_gl_pass) {
        GA_no = std::max(0.0, std::min(GA_no_gl[static_cast<size_t>(j)], 1.0));
        GB_no = channels_equal ? GA_no
                               : std::max(0.0, std::min(GB_no_gl[static_cast<size_t>(j)], 1.0));
      } else {
        GA_no = local_race_helper(t, parnA.data(), parA.data(), n_par, model_ctx,
                                  pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        GB_no = GA_no;
        if (!channels_equal)
          GB_no = local_race_helper(t, parnB.data(), parB.data(), n_par, model_ctx,
                                    pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        if (!R_FINITE(GA_no) || !R_FINITE(GB_no)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
      const double s_GA_yes = std::min(1.0, std::max(kMinSurv, GA_no + S_dec_A));
      const double s_GB_yes = std::min(1.0, std::max(kMinSurv, GB_no + S_dec_B));
      if (is_or_rule) {
        if (resp_code == 1) p_j = gA_yes * s_GB_yes + gB_yes * s_GA_yes; // yes
        else if (resp_code == 2) p_j = gA_no * GB_no + gB_no * GA_no;    // no
      } else { // AND
        const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
        const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
        const double s_GA_no = std::min(1.0, std::max(kMinSurv, GA_yes + S_dec_A));
        const double s_GB_no = std::min(1.0, std::max(kMinSurv, GB_yes + S_dec_B));
        if (resp_code == 1) p_j = gA_yes * GB_yes + gB_yes * GA_yes;      // yes
        else if (resp_code == 2) p_j = gA_no * s_GB_no + gB_no * s_GA_no; // no
      }
    } else if (is_xor_rule || is_id_rule) {
      double GA_yes, GB_yes;
      if (use_gl_pass) {
        // Derive GA_yes from GA_no_gl: GA_yes + GA_no + S_dec_A = 1.
        const double GA_no = std::max(0.0, std::min(GA_no_gl[static_cast<size_t>(j)], 1.0));
        const double GB_no = channels_equal ? GA_no
                                            : std::max(0.0, std::min(GB_no_gl[static_cast<size_t>(j)], 1.0));
        GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
        GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
      } else {
        GA_yes = local_race_helper(t, parA.data(), parnA.data(), n_par, model_ctx,
                                   pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        GB_yes = GA_yes;
        if (!channels_equal)
          GB_yes = local_race_helper(t, parB.data(), parnB.data(), n_par, model_ctx,
                                     pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        if (!R_FINITE(GA_yes) || !R_FINITE(GB_yes)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
      const double GA_no = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_yes));
      const double GB_no = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_yes));
      const double dAB = gA_yes * GB_yes + gB_yes * GA_yes;
      const double dAN = gA_yes * GB_no + gB_no * GA_yes;
      const double dNB = gA_no * GB_yes + gB_yes * GA_no;
      const double dNN = gA_no * GB_no + gB_no * GA_no;

      if (is_id_rule) {
        if (resp_code == 6) p_j = dAB;
        else if (resp_code == 4) p_j = dAN;
        else if (resp_code == 5) p_j = dNB;
        else if (resp_code == 3) p_j = dNN;
      } else { // XOR
        if (resp_code == 1) p_j = dAN + dNB; // yes
        else if (resp_code == 2) p_j = dAB + dNN; // no
      }
    } else {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    if (!(p_j > 0.0) || !R_FINITE(p_j)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
    } else {
      const double ll = std::log(p_j);
      ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
    }
  }

  double sum_ll = 0.0;
  if (expand.length() > 0) {
    const int n_out = expand.length();
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_logicalrules: trial_ll_out size mismatch (expand path).");
    }
    for (int i = 0; i < n_out; ++i) {
      double val = ll_unique[static_cast<size_t>(expand[i] - 1)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = val;
      sum_ll += val;
    }
  } else {
    const int n_out = n_unique_trials;
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_logicalrules: trial_ll_out size mismatch (compressed path).");
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      double val = ll_unique[static_cast<size_t>(j)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = val;
      sum_ll += val;
    }
  }
  return sum_ll;
}

// Port of redundant-target race likelihood adapted to the refactored scalar
// adapters and column-major parameter storage used in calc_ll_oo.
// This function is intentionally kept standalone for now (wrapper wiring can
// be added via a dedicated c_name tag when model wrappers are ready).
static double c_log_likelihood_redundant_target_race(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& ok_params,
    int n_acc,
    ContextForRaceModels* model_ctx,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    const RedundantTargetSharedState& shared,
    Rcpp::NumericVector* trial_ll_out) {
  if (!shared.valid || n_acc <= 0 || pdf1 == nullptr || cdf1 == nullptr || model_ctx == nullptr) {
    Rcpp::stop("c_log_likelihood_redundant_target_race: invalid redundant-target configuration.");
  }
  const int n_trials = shared.n_trials;
  if (pars.nrow() != n_trials || ok_params.size() != n_trials) {
    Rcpp::stop("c_log_likelihood_redundant_target_race: pars/ok_params sizes must match shared data.");
  }
  if (n_trials == 0) return 0.0;

  static const double kMinSurv = 1e-12;
  const int n_par = pars.ncol();
  const int n_unique_trials = shared.n_unique_trials;
  const double* pars_cm_ptr = pars.begin();
  const bool use_raw_local =
    (model_dfun_raw != nullptr &&
     model_pfun_raw != nullptr &&
     static_cast<int>(shared.rt_by_row.size()) == n_trials);
  std::vector<double> ll_unique(static_cast<size_t>(n_unique_trials), min_ll);

  GslIntegrationControls gsl_ctl = default_gsl_controls();
  gsl_ctl.try_qng_first_finite = true;
  gsl_ctl.qag_key = GSL_INTEG_GAUSS21;
  gsl_ctl.rel_tol = 1e-4;  // tighter than default unnecessary for log-likelihood accuracy

  // Change 3: skip the raw batch (and all integration) if every parameter row is invalid.
  bool any_ok = false;
  for (int i = 0; i < n_trials && !any_ok; ++i) any_ok = (bool)ok_params[i];
  if (!any_ok) {
    const int n_out = expand.length() > 0 ? (int)expand.length() : n_unique_trials;
    if (trial_ll_out != nullptr) for (int i = 0; i < n_out; ++i) (*trial_ll_out)[i] = min_ll;
    return static_cast<double>(n_out) * min_ll;
  }

  std::vector<double> logf_all;
  std::vector<double> logS_all;
  std::vector<int> all_mask;
  std::vector<int> isok_int_all;
  if (use_raw_local) {
    logf_all.resize(static_cast<size_t>(n_trials), min_ll);
    logS_all.resize(static_cast<size_t>(n_trials), min_ll);
    all_mask.assign(static_cast<size_t>(n_trials), 1);
    isok_int_all.resize(static_cast<size_t>(n_trials), 0);
    for (int i = 0; i < n_trials; ++i) isok_int_all[static_cast<size_t>(i)] = ok_params[i] ? 1 : 0;
    model_dfun_raw(shared.rt_by_row.data(), pars_cm_ptr, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logf_all.data(), min_ll, model_ctx);
    model_pfun_raw(shared.rt_by_row.data(), pars_cm_ptr, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logS_all.data(), min_ll, model_ctx);
  }

  // eval_pdf_cdf and safe_surv_at only capture read-only shared state; thread-safe.
  auto eval_pdf_cdf = [&](double t, int idx, const double* p_row, double& f_out, double& F_out) -> bool {
    if (use_raw_local && R_FINITE(t) && t > 0.0) {
      const double lf = logf_all[static_cast<size_t>(idx)];
      const double ls = logS_all[static_cast<size_t>(idx)];
      if (!R_FINITE(lf) || !R_FINITE(ls)) return false;
      const double f = (lf <= min_ll) ? 0.0 : std::exp(lf);
      double S = (ls <= min_ll) ? 0.0 : std::exp(ls);
      if (!(S >= 0.0) || !R_FINITE(S)) return false;
      if (S > 1.0) S = 1.0;
      const double F = clamp_cdf01_race(1.0 - S);
      if (!R_FINITE(f) || !R_FINITE(F) || emc2_isnan(F)) return false;
      f_out = f;
      F_out = F;
      return true;
    }
    const double f = pdf1(t, p_row, model_ctx);
    if (!R_FINITE(f) || f < 0.0) return false;
    double F = cdf1(t, p_row, model_ctx);
    F = clamp_cdf01_race(F);
    if (!R_FINITE(F) || emc2_isnan(F)) return false;
    f_out = f;
    F_out = F;
    return true;
  };

  auto safe_surv_at = [&](const double* p_row, double t_eval) -> double {
    double F = cdf1(t_eval, p_row, model_ctx);
    F = clamp_cdf01_race(F);
    if (!R_FINITE(F) || emc2_isnan(F)) return 0.0;
    return std::max(kMinSurv, 1.0 - F);
  };

  static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace_tls, gsl_ctl.retry_limit);
  std::vector<double> parA(static_cast<size_t>(n_par));
  std::vector<double> parB(static_cast<size_t>(n_par));
  std::vector<double> parN(static_cast<size_t>(n_par));
  std::vector<double> pars_rowmajor_buf(static_cast<size_t>(3 * n_par));
  int isok_buf[3] = {1, 1, 1};
  for (int j = 0; j < n_unique_trials; ++j) {
    const int idxA = shared.idxA[static_cast<size_t>(j)];
    const int idxB = shared.idxB[static_cast<size_t>(j)];
    const int idxN = shared.idxNogo[static_cast<size_t>(j)];
    const bool use_nogo = (idxN >= 0);

    if (idxA < 0 || idxB < 0 || !ok_params[idxA] || !ok_params[idxB] ||
        (use_nogo && !ok_params[idxN])) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxA, parA.data());
    copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxB, parB.data());
    if (use_nogo) copy_par_row_colmajor(pars_cm_ptr, n_trials, n_par, idxN, parN.data());

    const double t = shared.rt_unique[static_cast<size_t>(j)];
    const int cond = shared.cond_code[static_cast<size_t>(j)];
    const int resp = shared.resp_code[static_cast<size_t>(j)];
    const double LTj = shared.LT_unique[static_cast<size_t>(j)];
    const double UCj = shared.UC_unique[static_cast<size_t>(j)];
    const bool finite_rt = (R_FINITE(t) && t > 0.0);

    double p = 0.0;
    if (cond < 1 || cond > 3) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    if (finite_rt) {
      if (resp != 1) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
        continue;
      }

      double fA = 0.0, FA = 0.0, fB = 0.0, FB = 0.0, fN = 0.0, FN = 0.0;
      if (!eval_pdf_cdf(t, idxA, parA.data(), fA, FA) || !eval_pdf_cdf(t, idxB, parB.data(), fB, FB)) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
        continue;
      }
      if (use_nogo && !eval_pdf_cdf(t, idxN, parN.data(), fN, FN)) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
        continue;
      }

      const double S_A = std::max(kMinSurv, 1.0 - FA);
      const double S_B = std::max(kMinSurv, 1.0 - FB);
      const double S_N = use_nogo ? std::max(kMinSurv, 1.0 - FN) : 1.0;
      if (cond == 1) {
        p = fA * S_N;
      } else if (cond == 2) {
        p = fB * S_N;
      } else { // AB
        p = (fA * S_B + fB * S_A) * S_N;
      }
    } else {
      // Omission / no-response path: allow missing response code and explicit "nogo".
      if (resp != 0 && resp != 2) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
        continue;
      }

      const double S_A_inf = safe_surv_at(parA.data(), R_PosInf);
      const double S_B_inf = safe_surv_at(parB.data(), R_PosInf);
      if (!use_nogo) {
        if (cond == 1) p = S_A_inf;
        else if (cond == 2) p = S_B_inf;
        else p = S_A_inf * S_B_inf;
      } else {
        if (!R_FINITE(UCj)) {
          Rcpp::stop("RedundantTarget likelihood: finite UC is required for no-go/omission trials with a nogo accumulator.");
        }
        const double low = std::max(0.0, LTj);
        const double upp = UCj;
        const double S_A_uc = safe_surv_at(parA.data(), upp);
        const double S_B_uc = safe_surv_at(parB.data(), upp);
        const double S_N_uc = safe_surv_at(parN.data(), upp);
        double p_none = 0.0;
        if (cond == 1) p_none = S_A_uc * S_N_uc;
        else if (cond == 2) p_none = S_B_uc * S_N_uc;
        else p_none = S_A_uc * S_B_uc * S_N_uc;

        // Integrate probability that nogo wins before UC (and after LT if truncated).
        std::memcpy(pars_rowmajor_buf.data(), parN.data(), static_cast<size_t>(n_par) * sizeof(double));
        int n_active = 1;
        if (cond == 1) {
          std::memcpy(pars_rowmajor_buf.data() + n_par, parA.data(), static_cast<size_t>(n_par) * sizeof(double));
          n_active = 2;
        } else if (cond == 2) {
          std::memcpy(pars_rowmajor_buf.data() + n_par, parB.data(), static_cast<size_t>(n_par) * sizeof(double));
          n_active = 2;
        } else {
          std::memcpy(pars_rowmajor_buf.data() + n_par, parA.data(), static_cast<size_t>(n_par) * sizeof(double));
          std::memcpy(pars_rowmajor_buf.data() + 2 * n_par, parB.data(), static_cast<size_t>(n_par) * sizeof(double));
          n_active = 3;
        }
        for (int k = 0; k < n_active; ++k) isok_buf[k] = 1;
        double log_nogo_win = R_NegInf;
        if (upp > low) {
          log_nogo_win = integrate_for_kth_winner_rowmajor_cpp(
            1, pars_rowmajor_buf.data(), isok_buf, low, upp, pdf1, cdf1,
            n_active, n_par, gsl_ctl, model_ctx, w);
        }
        const double p_nogo_win = R_FINITE(log_nogo_win) ? std::exp(log_nogo_win) : 0.0;
        p = p_nogo_win + p_none;
      }
    }

    if (!(p > 0.0) || !R_FINITE(p)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
    } else {
      const double ll = std::log(p);
      ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
    }
  }

  double sum_ll = 0.0;
  if (expand.length() > 0) {
    const int n_out = expand.length();
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_redundant_target_race: trial_ll_out size mismatch (expand path).");
    }
    for (int i = 0; i < n_out; ++i) {
      double val = ll_unique[static_cast<size_t>(expand[i] - 1)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = val;
      sum_ll += val;
    }
  } else {
    const int n_out = n_unique_trials;
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_redundant_target_race: trial_ll_out size mismatch (compressed path).");
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      double val = ll_unique[static_cast<size_t>(j)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = val;
      sum_ll += val;
    }
  }
  return sum_ll;
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
    void* model_context_for_funcs,          // Context for model functions (e.g., contains posdrift for LBA)
    bool all_finite_trials,             // Data-only hint: all trials finite/in-bounds/known response
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,            // batch log-survivor at scalar t (for truncation norms)
    RaceSharedState* shared,             // optional pre-computed per-data state (nullptr = compute per-call)
    NumericVector* trial_ll_out
) {

  // Reuse one workspace per thread across particles/model fits.
  static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
  GslWorkspacePtr& workspace = workspace_tls;

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
  GslIntegrationControls gsl_ctl = default_gsl_controls();
  gsl_ctl.try_qng_first_finite = true;   // try fixed-point QNG before adaptive QAG on finite intervals
  gsl_ctl.qag_key = GSL_INTEG_GAUSS21;  // fallback rule when QNG fails
  gsl_ctl.rel_tol = 1e-4;               // sufficient accuracy for log-likelihood; reduces QAG subdivisions
  int n_lR_j = n_lR;
  Rcpp::NumericVector rts_dadm = dadm["rt"];
  Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
  // Use pre-allocated std::vector when shared state is available; avoids per-particle
  // R-heap allocation.  The buffer does NOT need re-initialisation between particles:
  // dfun_raw/pfun_raw overwrite every finite-trial slot before it is read, and
  // other-trial slots are computed independently in the "other trials" loop below.
  std::vector<double> lds_local;
  double* lds_ptr;
  if (use_shared && static_cast<int>(shared->res_buf.size()) == n_trials) {
    lds_ptr = shared->res_buf.data();
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
  const bool defective_upper_tail = ctx->defective_upper_tail;
  
  // log_surv_cm: log S(t) = sum_k log(1-F_k(t)) read directly from column-major pars.
  // Used by the "other trials" path for analytical survivor calls — avoids the
  // row-major copy that fill_trial_buffers performs.
  auto log_surv_cm = [&](double t, int start_row_idx, int n_lR_j) -> double {
    if (t == R_PosInf) {
      auto* race_ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
      // Proper race models have S(Inf)=0 => log S(Inf) = -Inf.
      // Defective-tail models (e.g., LBAIO / BAwL) retain mass at +Inf.
      if (race_ctx && race_ctx->defective_upper_tail) {
        double log_p = 0.0;
        double par_buf_inf[16];
        for (int k = 0; k < n_lR_j; ++k) {
          const int row = start_row_idx + k;
          if (!isok[row]) return R_NegInf;
          for (int c = 0; c < n_par; ++c)
            par_buf_inf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
          double Fk_inf = cdf1(R_PosInf, par_buf_inf, model_context_for_funcs);
          Fk_inf = clamp_cdf01_race(Fk_inf);
          const double ll = safe_log1m_race(Fk_inf);
          if (!emc2_isfinite(ll)) return R_NegInf;
          log_p += ll;
        }
        return log_p;
      }
      return R_NegInf;
    }

    double logS = 0.0;
    double par_buf[16];
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start_row_idx + k;
      if (!isok[row]) return R_NegInf;
      for (int c = 0; c < n_par; ++c)
        par_buf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
      double Fk = cdf1(t, par_buf, model_context_for_funcs);
      Fk = clamp_cdf01_race(Fk);
      const double ll = safe_log1m_race(Fk);
      if (!R_FINITE(ll)) return R_NegInf;
      logS += ll;
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
        shared->ok_int_buf[static_cast<size_t>(i)] = isok[i] ? 1 : 0;
      idx_win_ptr  = shared->idx_win.data();
      idx_loss_ptr = shared->idx_loss.data();
      isok_ptr     = shared->ok_int_buf.data();
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
      // Pass 1: scan for truncated trials; check uniformity of LT and UT separately.
      // We can batch the normaliser whenever all truncated finite-RT trials share the
      // same LT value AND the same UT value (each may be 0/Inf trivially).
      double uniform_LT = -1.0;   // -1 = not yet observed (LT >= 0 always)
      double uniform_UT = -1.0;   // -1 = not yet observed
      bool uniform_LT_ok = true;
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
          if (uniform_LT < 0.0) uniform_LT = LTj;
          else if (LTj != uniform_LT) uniform_LT_ok = false;
          if (uniform_UT < 0.0) uniform_UT = UTj;
          else if (UTj != uniform_UT) uniform_UT_ok = false;
        }
      }
      // Resolve sentinels: if no truncated trial was seen, treat as trivial.
      if (uniform_LT < 0.0) uniform_LT = 0.0;
      if (uniform_UT < 0.0) uniform_UT = R_PosInf;

      if (any_trunc && uniform_LT_ok && uniform_UT_ok && !has_RACE_col) {
        // Pass 2: batch-compute logZ = log_diff_exp(logS(LT), logS(UT)) for all
        // truncated trials simultaneously.
        // RACE models are excluded: lba/rdm/lnr_logS_at_t always loops over the
        // global n_lR, but RACE trials may have fewer active accumulators (inactive
        // rows are NA-masked).  The per-trial scalar path handles n_lR_j correctly.
        //
        // Trivial endpoints:
        //   LT == 0   → logS(LT) = 0  (S(0)=1 for any proper distribution)
        //   UT == Inf → logS(UT) = -Inf (S(∞)=0 for any proper distribution)
        // For non-trivial endpoints we call logS_at_t in batch.
        //
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

        // logS_LT: 0 when LT==0 (trivial), otherwise computed in batch.
        std::vector<double> logS_LT_vec(static_cast<size_t>(n_unique_trials), 0.0);
        if (uniform_LT != 0.0) {
          logS_at_t(uniform_LT,
                    pars.begin(), n_trials, n_lR, n_par,
                    trunc_mask.data(), n_unique_trials,
                    isok_all_ptr_trunc,
                    model_context_for_funcs,
                    logS_LT_vec.data());
        }

        // logS_UT: -Inf when UT==Inf (trivial for proper distributions), else computed.
        std::vector<double> logS_UT_vec(static_cast<size_t>(n_unique_trials), R_NegInf);
        if (uniform_UT != R_PosInf) {
          logS_at_t(uniform_UT,
                    pars.begin(), n_trials, n_lR, n_par,
                    trunc_mask.data(), n_unique_trials,
                    isok_all_ptr_trunc,
                    model_context_for_funcs,
                    logS_UT_vec.data());
        }

        // logZ = log(S(LT) - S(UT)) = log_diff_exp(logS_LT, logS_UT)
        // log_diff_exp handles R_NegInf inputs correctly:
        //   log_diff_exp(0, -Inf)  = 0   (LT=0,  UT=Inf → P=1)
        //   log_diff_exp(x, -Inf)  = x   (UT=Inf  → P=S(LT))
        //   log_diff_exp(0, y)     = log(1-exp(y)) (LT=0 → original UT-only case)
        // NaN from logS_at_t (bad params) propagates through log_diff_exp as NA_REAL,
        // which is caught by emc2_isnan and sent to the per-trial fallback.
        logZ_batch.assign(static_cast<size_t>(n_unique_trials), 0.0);
        static const double kLogProbEps = std::log(std::numeric_limits<double>::epsilon());
        for (int j = 0; j < n_unique_trials; ++j) {
          if (!trunc_mask[static_cast<size_t>(j)]) continue;
          const double logS_LT_j = logS_LT_vec[static_cast<size_t>(j)];
          const double logS_UT_j = logS_UT_vec[static_cast<size_t>(j)];
          const double logP = log_diff_exp(logS_LT_j, logS_UT_j);
          if (R_FINITE(logP) && logP > kLogProbEps) {
            logZ_batch[static_cast<size_t>(j)] = logP;
          } else {
            // NaN (bad params), catastrophic cancellation, or P≈0:
            // fall back to the per-trial scalar/GSL path.
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
            if (!emc2_isnan(bz)) {       // NaN flags a fallback case
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
          bool hit_min_ll = false;
          for (int k = 0; k < n_lR_j; ++k) {
            double v = lds_ptr[start_row_idx + k];
            if (v <= min_ll) hit_min_ll = true;
            current_trial_ll_sum += v;
          }
          if (NumericVector::is_na(log_Z_this) || !R_FINITE(log_Z_this)) {
            ll_unique[unique_trial_idx] = min_ll;
          } else if (hit_min_ll) {
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
          // Fast exit: when no upper censoring (UCj=Inf) and model has no
          // defective upper tail, P(RT=Inf)=0. Avoid log_diff_exp(-Inf,-Inf)
          // and degenerate integration bounds [Inf,Inf].
          if (!defective_upper_tail && !R_FINITE(lower_for_trial)) {
            // current_ll_val stays at min_ll
          } else {
            const double logP = (!defective_upper_tail)
              ? log_diff_exp(log_surv_cm(lower_for_trial, start_row_idx, n_lR_j),
                             log_surv_cm(upper_for_trial, start_row_idx, n_lR_j))
              : log_surv_cm(lower_for_trial, start_row_idx, n_lR_j); // include intrinsic never-finish mass
            if (R_FINITE(logP) && logP > log_prob_eps) {
              current_ll_val = logP;
            } else { // numerical integration fallback
              for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
                current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
              }
              if (defective_upper_tail) { // include defective upper-tail mass
                const double log_p_I = log_surv_cm(R_PosInf, start_row_idx, n_lR_j);
                current_ll_val = log_sum_exp(current_ll_val, log_p_I);
              }
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
        const double logP2 = (!defective_upper_tail)
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
          if (defective_upper_tail) {
            const double log_p_I = log_surv_cm(R_PosInf, start_row_idx, n_lR_j);
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
        double log1m_pC = log1m(pC); // todo double check right helper
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

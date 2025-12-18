#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "model_SS.h"
#include "trend.h"
#include "utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include <cmath>
#include <string>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <iomanip>
#include <array>
#include <limits>
using namespace Rcpp;

// Forward declaration: stop-signal likelihood with deadline censoring (rt==Inf uses UC).
double c_log_likelihood_ss_cens_trunc(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::DataFrame& dadm,
    int n_trials,
    const Rcpp::LogicalVector& winner,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& isok,
    int n_lR,
    const std::string& ss_type);

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
                      const std::vector<TransformSpec>& full_specs) {

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
      String cur_name(cur_names[j]);
      NumericVector p_mult_design(n_trials, p_vector[cur_name]);
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
  if (has_trend && premap) {
    CharacterVector names_premap = collect_trend_param_names_phase(trend, "premap");
    if (names_premap.size() > 0) {
      pars = submat_rcpp_col(pars, !contains_multiple(p_types, names_premap));
    }
  }
  return pars;
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, const std::vector<PreTransformSpec>& p_specs,
                              CharacterVector p_types, List designs, int n_trials, DataFrame data, List trend,
                              const std::vector<TransformSpec>& full_specs){
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
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, full_specs);
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

double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials, IntegerVector expand,
                            double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector lls(n_trials);
  NumericVector lls_exp(n_out);
  lls = d_DDM_Wien(rts, R, pars, is_ok);
  lls_exp = c_expand(lls, expand); // decompress
  // lls_exp = lls;
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return(sum(lls_exp));
}

// ZH - this has been replaced by the more robust cens/log function
double c_log_likelihood_race(NumericMatrix pars, DataFrame data,
                             NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                             NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                             const int n_trials, LogicalVector winner, IntegerVector expand,
                             double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  CharacterVector R = data["R"];
  NumericVector lR = data["lR"];
  NumericVector lds_exp(n_out);
  const int n_lR = unique(lR).length();
  if(sum(contains(data.names(), "RACE")) == 1){
    NumericVector NACC = data["RACE"];
    CharacterVector vals_NACC = NACC.attr("levels");
    for(int x = 0; x < pars.nrow(); x++){
      // subtract 1 because R is 1 coded
      if(lR[x] > atoi(vals_NACC[NACC[x]-1])){
        pars(x,0) = NA_REAL;
      }
    }
  }
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll), is_ok)); //first for compressed
  lds[winner] = win;
  if(n_lR > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok)); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds)] = min_ll;

  if(n_lR > 1){
    // LogicalVector winner_exp = c_bool_expand(winner, expand);
    NumericVector ll_out = lds[winner];
    NumericVector lds_los = lds[!winner];
    if(n_lR == 2){
      ll_out = ll_out + lds_los;
    } else{
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_lR -1), (z+1) * (n_lR -1) -1)]);
      }
    }

    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    ll_out = c_expand(ll_out, expand); // decompress
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    lds_exp = c_expand(lds, expand); // decompress
    return(sum(lds_exp));
  }
}

// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
            List designs, String type, List bounds, List transforms, List pretransforms,
            CharacterVector p_types, double min_ll, List trend){
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

  if(type == "DDM"){
    IntegerVector expand = data.attr("expand");
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      // Precompute specs
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type == "MRI" || type == "MRI_AR1"){
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
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
  } else{
    // ZH I have adjusted a lot here. I made new race model pointers that help branch things like intrinsic omissions LBA, and the new GNG models
    // The "1" versions handle scalar likelihoods that do not call Rcpp functions,which speeds up the GSL integrations used in censoring, truncation, and GNG 
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Standardize incoming model type string from R
    // For race models: e.g. "LBA", "LBAIO", "LBAGNG", "RDM", "RDMGNG", "LNR", "LNRGNG"
    // For stop-signal models (C++): e.g. "SSexG", "SShybrid"
    std::string type_std = Rcpp::as<std::string>(Rcpp::wrap(type));
    const bool is_ss = (type_std.find("SS") != std::string::npos);

 	RacePdfFun model_dfun_ptr = nullptr;
 	RaceCdfFun model_pfun_ptr = nullptr;
 	RacePdf1Fun pdf1_ptr = nullptr;
 	RaceCdf1Fun cdf1_ptr = nullptr;
 	ContextForRaceModels current_model_ctx;
 	current_model_ctx.min_lik_for_pdf = std::exp(min_ll);
 	current_model_ctx.use_posdrift = true; // Default: LBA uses posdrift, RDM/LNR ignore this context field.
 	current_model_ctx.gng = false; // Default: Do not use go/no-go likelihoods by default
    if (!is_ss) {
      // Determine adapter functions and specific LBA posdrift setting based on type_std
      if (type_std.find("LBA") != std::string::npos) {
        model_dfun_ptr = &lba_dfun_adapter;
        model_pfun_ptr = &lba_pfun_adapter;
        pdf1_ptr = &dlba_scalar;
        cdf1_ptr = &plba_scalar;
        // Check for the 'IO' (Implicit Omissions / no posdrift) flag in the original type_std
        if (type_std.find("IO") != std::string::npos) {
          current_model_ctx.use_posdrift = false;
        }
      } else if (type_std == "RDM") {
        model_dfun_ptr = &rdm_dfun_adapter;
        model_pfun_ptr = &rdm_pfun_adapter;
        pdf1_ptr = &drdm_scalar;
        cdf1_ptr = &prdm_scalar;
      } else if (type_std.find("LNR") != std::string::npos) {
        model_dfun_ptr = &lnr_dfun_adapter;
        model_pfun_ptr = &lnr_pfun_adapter;
        pdf1_ptr = &dlnr_scalar;
        cdf1_ptr = &plnr_scalar;
      } else {
        Rcpp::stop("Unsupported race model type string in calc_ll: " + type_std);
      }
      if (type_std.find("GNG") != std::string::npos) {
        current_model_ctx.gng = true;
      }
    }
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    for (int i = 0; i < n_particles; ++i) {
        p_vector = p_matrix(i, Rcpp::_);
        if (i == 0) {
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      if (is_ss) {
        lls[i] = c_log_likelihood_ss_cens_trunc(pars, data,
                                               n_trials, winner, expand,
                                               min_ll, is_ok, n_lR, type_std);
      } else {
        lls[i] = c_log_likelihood_race_cens_trunc(pars, data,
                                                 model_dfun_ptr, model_pfun_ptr,
                                                 pdf1_ptr, cdf1_ptr,
                                                 n_trials,
                                                 winner, expand, min_ll, is_ok, n_lR,
                                                 &current_model_ctx);
      }
    }
  }
  return(lls);
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


bool row_is_finite(const Rcpp::NumericMatrix& mat, int row) {
    for (int j = 0; j < mat.ncol(); ++j) {
        if (!R_finite(mat(row, j))) return false;
    }
    return true;
}

// Numerical integration helper using GSL
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
	  LogicalVector isok,
    double low,
    double upp,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    int n_lR_j,
    int n_par,
    double epsilon,
    void* model_specific_context,
    gsl_integration_workspace* w) {

    if (low >= upp && !(low == 0 && upp == R_PosInf)) return R_NegInf;
    if (k_winner_idx < 1 || k_winner_idx > (n_lR_j)) {
      Rcpp::Rcerr << "Warning: Invalid k_winner_idx in integrate_for_kth_winner_cpp: " << k_winner_idx << std::endl;
      return R_NegInf;
    }
	
    if (w == nullptr) {
      Rcpp::stop("integrate_for_kth_winner_cpp: GSL workspace is null.");
    }

    std::vector<double> pars_rowmajor(static_cast<size_t>(n_lR_j) * n_par);
    for (int r = 0; r < n_lR_j; ++r) {
      for (int c = 0; c < n_par; ++c) {
        pars_rowmajor[static_cast<size_t>(r) * n_par + c] = p_all_acc(r, c);
      }
    }

    std::vector<int> isok_int(n_lR_j, 0);
    for (int r = 0; r < n_lR_j; ++r) isok_int[r] = isok[r] ? 1 : 0;

    gsl_function F;
    gsl_race_params_scalar params_struct;
    params_struct.pars = pars_rowmajor.data();
    params_struct.n_lR = n_lR_j;
    params_struct.n_par = n_par;
    params_struct.winner_idx0 = k_winner_idx - 1;
    params_struct.isok = isok_int.data();
    params_struct.pdf1 = pdf1;
    params_struct.cdf1 = cdf1;
    params_struct.ctx = model_specific_context;

    F.function = &gsl_f_race_scalar;
    F.params = &params_struct;

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status;
    double result, error;
    if (upp == R_PosInf) {
        if (low < 0) low = 0; // QAGIU requires a >= 0
        if (low >= R_PosInf) {
             result = 0; status = GSL_SUCCESS;
        } else {
            status = gsl_integration_qagiu(&F, low, 0, epsilon, 1000, w, &result, &error);
        }
    } else {
        status = gsl_integration_qags(&F, low, upp, 0, epsilon, 1000, w, &result, &error);
    }

    gsl_set_error_handler(old_handler);

    if (status != GSL_SUCCESS) {
        // On integration failure, treat as zero probability on this interval
        return R_NegInf;
    }

    if (!(result > 0) || !R_finite(result)) {
      return R_NegInf;
    }
    return std::log(result);
}

// New common normaliser
double get_trunc_normaliser_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
	  LogicalVector isok_trial, 
    double LT, double UT,
    int n_lR,
    int n_par,
    double epsilon,
    void* model_specific_context,
    gsl_integration_workspace* w){
    double total_log_mass = R_NegInf; 
    for(int k=1;k<=n_lR;++k) {
        double log_mass_k  = integrate_for_kth_winner_cpp(k,
                                                          p_all_acc,
                                                          isok_trial,
                                                          LT,
                                                          UT,
                                                          pdf1,
                                                          cdf1,
                                                          n_lR,
                                                          n_par,
                                                          epsilon,
                                                          model_specific_context,
                                                          w);
		    // Accumulate the log masses using the log-sum-exp trick
        total_log_mass = log_sum_exp(total_log_mass, log_mass_k);
	}
    return total_log_mass;
}

inline void append_par_block(std::ostringstream& oss,
                            int start_row,
                            const Rcpp::NumericMatrix& pars,
                            int n_lR,
                            int n_par) {
  for (int r = 0; r < n_lR; ++r) {
    for (int c = 0; c < n_par; ++c) {
      const double v = pars(start_row + r, c);
      if (Rcpp::NumericVector::is_na(v)) oss << "NA";
      else oss << v;
      oss << ',';
    }
    oss << '|';
  }
}

inline std::string make_logZ_key(int start_row,
                                 double low,
                                 double upp,
                                 const Rcpp::NumericMatrix& pars,
                                 int n_lR,
                                 int n_par) {
  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  oss << std::setprecision(17);
  oss << "Z|" << low << '|' << upp << '|';
  append_par_block(oss, start_row, pars, n_lR, n_par);
  return oss.str();
}

inline std::string make_winner_integral_key(int start_row,
                                            double low,
                                            double upp,
                                            int k_winner_idx,
                                            const Rcpp::NumericMatrix& pars,
                                            int n_lR,
                                            int n_par) {
  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  oss << std::setprecision(17);
  oss << "I|" << k_winner_idx << '|' << low << '|' << upp << '|';
  append_par_block(oss, start_row, pars, n_lR, n_par);
  return oss.str();
}

// Main C++ function for censored/truncated race likelihood calculation
// This function is now the unified entry point for all race models (LBA, RDM, LNR),
// whether they are standard or explicitly handling censoring/truncation.
// It uses batching for finite RTs and iterative processing for others (censored/NA RTs).
double c_log_likelihood_race_cens_trunc(
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
    void* model_context_for_funcs           // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
) {
  // Fetch censoring and truncation values from dadm columns. These are passed across all rows for ease of access, but should be identical at least at the subject level as they don't correspond to a data entry. Attributes probably a better fit here but clunky.
	auto workspace = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>(
	    gsl_integration_workspace_alloc(1000),
	    &gsl_integration_workspace_free);
	if (!workspace) Rcpp::stop("Failed to allocate GSL integration workspace.");
	Rcpp::NumericVector LT = dadm["LT"];
	Rcpp::NumericVector UT = dadm["UT"]; 
	Rcpp::NumericVector LC = dadm["LC"];    
	Rcpp::NumericVector UC = dadm["UC"];
	std::unordered_map<std::string,double> log_integral_cache; // truncation + go/no-go caching
	double integration_epsilon = 1e-7; // Tolerance for GSL integration
	int n_lR_j = n_lR;
	const int n_out = expand.length();
	Rcpp::NumericVector rts_dadm = dadm["rt"];
	Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
	Rcpp::IntegerVector lR_dadm = dadm["lR"];
	NumericVector lds(n_trials,min_ll); // initialise at min_ll
	NumericVector lds_exp(n_out);
	// If a RACE column exists, set parameters of accumulators not present on a
  // given trial to NA so the density functions return zero for them. This
  // mirrors logic from the old c_log_likelihood_race implementation.
	IntegerVector RACE(n_trials, n_lR);
	Rcpp::LogicalVector RACE_mask(n_trials, true); // Mask for all dadm rows
	if (dadm.containsElementNamed("RACE")) {
		// factor codes (1-based) for each row
		Rcpp::IntegerVector race_idx = dadm["RACE"];
		// character levels ("2", "3", ...)
		Rcpp::CharacterVector race_levels = race_idx.attr("levels");
		for (int row = 0; row < pars.nrow(); ++row) {
			// how many accumulators for this trial
			int n_lR_this_trial = std::stoi(
			  Rcpp::as<std::string>(race_levels[ race_idx[row] - 1 ])
			);
			RACE[row]=n_lR_this_trial;
			// lR_dadm is the (1-based) index of *this* accumulator on the trial
			if (lR_dadm[row] > n_lR_this_trial) {
			  // accumulator not present - blank its parameter row
			  std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
			  RACE_mask[row]=false;
			}
		}
	}
	if (n_trials == 0) return 0.0; // No data, no likelihood

	if (n_lR <= 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: n_lR must be positive and correctly determined before this call.");
	if (n_trials % n_lR != 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: dadm nrows not a multiple of n_lR.");

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
	Rcpp::NumericVector ll_unique(n_unique_trials);
	//ll_unique.fill(min_ll); // Initialize all to min_ll
	Rcpp::NumericVector pC_values(n_unique_trials);
	if (use_pC) {
		for (int j = 0; j < n_unique_trials; ++j) {
			pC_values[j] = pars(j * n_lR, pc_col);
		}
	}

	// Parameter matrix and validity vector checks
	if (pars.nrow() != n_trials) {
		Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_trials: " << n_trials << std::endl;
		Rcpp::stop("c_log_likelihood_race_cens_trunc: pars matrix dimensions do not match total dadm rows.");
	}
	if (isok.size() != pars.nrow()) {
		Rcpp::stop("c_log_likelihood_race_cens_trunc: isok size does not match "
				"pars matrix rows.");
	}
	if (winner.size() != pars.nrow()) {
		Rcpp::stop("c_log_likelihood_race_cens_trunc: isok size does not match pars matrix rows.");
	}

	Rcpp::LogicalVector finite_rt_mask(n_trials, false); // Mask for all dadm rows
	std::vector<int> finite_dadm_rows_indices; // Stores actual dadm/pars row indices for finite RTs
	finite_dadm_rows_indices.reserve(n_trials); // Reserve space
	double log_Z_this = 0;  // Default inv_Z if no truncation. Should never be used but here as a precaution.
	double lower_for_trial = 0; 
	double upper_for_trial = R_PosInf;
	ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
	// If posdrift=true we'll add the probability of an intrinsic omission in relevant sections
	bool posdrift=ctx->use_posdrift;
	bool gng=ctx->gng;
	// Batch process finite rt trials
	std::vector<int> finite_rt_unique_trial_indices; // Stores *unique trial* indices for finite RTs
	std::vector<int> other_unique_trial_indices;   // Stores *unique trial* indices for other cases
  finite_rt_unique_trial_indices.reserve(n_unique_trials);
  other_unique_trial_indices.reserve(n_unique_trials);
  // Categorize unique trials based on RT properties and parameter validity

  for (int j = 0; j < n_unique_trials; ++j) {
      int start_row_idx = j * n_lR; // Starting row in dadm/pars for this unique trial
      double rt_j = rts_dadm[start_row_idx]; // RT for this unique trial
      int R_j_idx = R_idxs_dadm[start_row_idx]; // Winner index (1-based from R factor)
		  int n_lR_j = RACE[start_row_idx];
      // Criteria for batch processing: finite, positive RT, within truncation bounds [LT, UT], and known winner
      if (R_FINITE(rt_j) && rt_j > 0 && rt_j >= LT[start_row_idx] && rt_j <= UT[start_row_idx] && R_j_idx != NA_INTEGER) {
			  finite_rt_unique_trial_indices.push_back(j); // Store unique trial index
        for(int k_acc = 0; k_acc < n_lR_j; ++k_acc) {
            int dadm_row_idx = start_row_idx + k_acc;
				    finite_rt_mask[dadm_row_idx] = true;
            finite_dadm_rows_indices.push_back(dadm_row_idx);
        }
        }  else {
            other_unique_trial_indices.push_back(j); // All other cases (NA, Inf, -Inf, outside bounds, or unknown winner with finite RT)
        }
    }
	
	  Rcpp::NumericMatrix finite_rt_pars;
    if (!finite_dadm_rows_indices.empty()) {
        finite_rt_pars = Rcpp::NumericMatrix(finite_dadm_rows_indices.size(), pars.ncol());
        Rcpp::CharacterVector par_colnames = Rcpp::colnames(pars);
        Rcpp::colnames(finite_rt_pars) = par_colnames;
        for(size_t i = 0; i < finite_dadm_rows_indices.size(); ++i) {
            finite_rt_pars.row(i) = pars.row(finite_dadm_rows_indices[i]);
        }
    }
    // Important: rts_dadm, isok, winner used for indexing lds calculation
    // should be sliced by finite_rt_mask, NOT finite_dadm_rows_indices directly
    // if the model functions expect a dense rt vector and pars matrix.
    // The current model functions (e.g. dlba_c) take full rts/pars and an 'idx' mask.
    // So, we pass the full rts_dadm, full pars, and finite_rt_mask as 'idx'.
    // However, the 'pars' argument to model_dfun/pfun should be the one corresponding to the RTs selected by finite_rt_mask.
    // This means finite_rt_pars is the correct parameter matrix to pass if the RTs are rts_dadm[finite_rt_mask].
	  Rcpp::NumericVector rts_for_finite_batch = rts_dadm[finite_rt_mask];
    Rcpp::LogicalVector winner_for_finite_batch = winner[finite_rt_mask];
    Rcpp::LogicalVector isok_for_finite_batch = isok[finite_rt_mask];
		  if (rts_for_finite_batch.size() > 0) { // Only calculate if there are any finite RTs
			  Rcpp::NumericVector win = log(model_dfun(rts_for_finite_batch, finite_rt_pars, winner_for_finite_batch, isok_for_finite_batch, model_context_for_funcs));
		    // Place results into the correct positions in 'lds'
	        int current_lds_idx = 0;
	        for(int i = 0; i < n_trials; ++i) {
	            if(finite_rt_mask[i] && winner[i]) {
	                lds[i] =win[current_lds_idx++];
	            }
	        }
			    // Prepare survivor functions for losers, one matrix per accumulator
			    if (n_lR > 1) {
				    Rcpp::NumericVector loss_cdf = model_pfun(rts_for_finite_batch,
				                                         finite_rt_pars,
				                                         !winner_for_finite_batch,
				                                         isok_for_finite_batch,
				                                         model_context_for_funcs);

            // Compute log(1 - CDF) directly to avoid log/exp round-trips.
            Rcpp::NumericVector loss(loss_cdf.size());
            for (int j = 0; j < loss_cdf.size(); ++j) {
              double cdf = loss_cdf[j];
              if (!R_finite(cdf)) { loss[j] = R_NegInf; continue; }
              if (cdf <= 0.0) { loss[j] = 0.0; continue; }           // log(1)
              if (cdf >= 1.0) { loss[j] = R_NegInf; continue; }       // log(0)
              if (cdf > 1.0 - 1e-15) cdf = 1.0 - 1e-15;
              loss[j] = std::log1p(-cdf);
            }
            
            current_lds_idx = 0; // Reset for iterating through loss_contributions
            for(int i = 0; i < n_trials; ++i) {
              if(finite_rt_mask[i] && !winner[i]) {
                double tmp_ll = loss[current_lds_idx++];
                if (tmp_ll < min_ll) tmp_ll = min_ll; // floor ll
                lds[i] = tmp_ll;
              }
            }
          }
        // Apply truncation correction and calculate log-likelihood for each trial in the batch
          for (size_t i = 0; i < finite_rt_unique_trial_indices.size(); ++i) {
              int unique_trial_idx = finite_rt_unique_trial_indices[i];
              int start_row_idx = unique_trial_idx*(n_lR);
              n_lR_j = RACE[start_row_idx];
              Rcpp::NumericMatrix p_all_acc_for_trial = pars(Rcpp::Range(start_row_idx, start_row_idx + n_lR_j - 1), Rcpp::_);
              Rcpp::LogicalVector isok_for_trial = isok[Rcpp::Range(start_row_idx, start_row_idx + n_lR_j - 1)];
              lower_for_trial = LT[start_row_idx];
              upper_for_trial = UT[start_row_idx];
              if (LT[start_row_idx]!=0 || UT[start_row_idx] != R_PosInf) { // if truncation
                std::string key = make_logZ_key(start_row_idx, lower_for_trial, upper_for_trial, pars, n_lR_j, pars.ncol());
                auto hit = log_integral_cache.find(key);
                if (hit != log_integral_cache.end()) { // check if we've integrated this exact parameter set already and re-use result if we have
                  log_Z_this = hit->second;
                } else { // otherwise compute the integral and cache it
                    log_Z_this  = get_trunc_normaliser_cpp(p_all_acc_for_trial,
                                                      pdf1,
                                                      cdf1,
                                                      isok_for_trial,
                                                      lower_for_trial,
                                                      upper_for_trial,
                                                      n_lR_j,
                                                      pars.ncol(),
                                                      integration_epsilon,
                                                      model_context_for_funcs,
                                                      workspace.get());
                    log_integral_cache.emplace(std::move(key), log_Z_this);
                  }
                }
                double current_trial_ll_sum = 0.0;
                for(int k=0; k < n_lR_j; ++k) {
                  int dadm_row_idx = start_row_idx + k; // Correct index into lds
                  current_trial_ll_sum += lds[dadm_row_idx];
                }

                if ((LT[start_row_idx] != 0 || UT[start_row_idx] != R_PosInf) && (NumericVector::is_na(log_Z_this) || !R_FINITE(log_Z_this))) {
                              // If truncation active but inv_Z is bad, probability is effectively zero. Could mean the probability of observing an untruncated RT is zero, which would be bad.
                              ll_unique[unique_trial_idx] = min_ll;
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
	// TODO update Case 3 to handle go/no-go withheld response
	double current_ll_val;
    for (size_t i = 0; i < other_unique_trial_indices.size(); ++i) {
		int unique_trial_idx = other_unique_trial_indices[i];
		int start_row_idx = unique_trial_idx*n_lR;
		n_lR_j = RACE[start_row_idx];
		Rcpp::NumericMatrix p_all_acc_for_trial = pars(Rcpp::Range(start_row_idx, start_row_idx + n_lR_j - 1), Rcpp::_);
        Rcpp::LogicalVector isok_for_trial = isok[Rcpp::Range(start_row_idx, start_row_idx + n_lR_j - 1)];
		Rcpp::LogicalVector winner_for_trial = winner[Rcpp::Range(start_row_idx, start_row_idx + n_lR_j - 1)];
        double rt_j = rts_dadm[start_row_idx];
        int R_j_idx = R_idxs_dadm[start_row_idx];
        current_ll_val = R_NegInf;
		bool slow_censor_closed_form = false; // true when we compute P(Tmin > UC) via survivor product (already includes intrinsic omissions when posdrift=false)

		auto integrate_cached = [&](int k_winner_1based, double low, double upp) -> double {
		  std::string key = make_winner_integral_key(start_row_idx,
		                                             low,
		                                             upp,
		                                             k_winner_1based,
		                                             pars,
		                                             n_lR_j,
		                                             pars.ncol());
		  auto hit = log_integral_cache.find(key);
		  if (hit != log_integral_cache.end()) return hit->second;

		  double log_val = integrate_for_kth_winner_cpp(k_winner_1based,
		                                               p_all_acc_for_trial,
		                                               isok_for_trial,
		                                               low,
		                                               upp,
		                                               pdf1,
		                                               cdf1,
		                                               n_lR_j,
		                                               pars.ncol(),
		                                               integration_epsilon,
		                                               model_context_for_funcs,
		                                               workspace.get());
		  log_integral_cache.emplace(std::move(key), log_val);
		  return log_val;
		};
        // Case 2: Fast censoring (RT=-Inf). Probability is integral from LT to LC.
        if (rt_j == R_NegInf) {
			lower_for_trial = LT[start_row_idx];
			upper_for_trial = LC[start_row_idx];
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                current_ll_val = integrate_cached(R_j_idx, lower_for_trial, upper_for_trial);
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
                    double log_l_k = integrate_cached(k_win, lower_for_trial, upper_for_trial);
					current_ll_val =  log_sum_exp(current_ll_val, log_l_k);
					}
            }
        // Case 3: Slow censoring (RT = Inf). Probability is integral from UC to UT.
		// ZH added Go/No-Go Race logic here.
        } else if (rt_j == R_PosInf) {
			if (gng) {
				const double D = UC[start_row_idx];
				lower_for_trial = LT[start_row_idx]; // maybe needs to be 0?
				upper_for_trial = D; // the censor value is the upper limit of a winning no-go, or the lower limit of a losing go/no-go
				// Identify which accumulator is the NOGO channel from winner_for_trial tag
				int k_nogo = -1; // 1-based index for integrate_for_kth_winner_cpp
				int n_true = 0;
				for (int k = 0; k < n_lR_j; ++k) {
				  if (winner_for_trial[k]) { n_true++; k_nogo = k + 1; }
				}
				if (n_true != 1) {
				  Rcpp::stop("No winner identified in go/no-go withheld response");
				} 
				// Term A: No-Go accumulator wins before the deadline D and before the Go accumulator
				const double logA = integrate_cached(k_nogo, lower_for_trial, D);
				
				// Term B: Neither accumulator finishes before the deadline = product of survivors			
				double logB = 0.0;
				std::vector<double> par_row(static_cast<size_t>(pars.ncol()));
				for (int k = 0; k < n_lR_j; ++k) {
					if (!isok_for_trial[k]) { logB = R_NegInf; break; }
					for (int c = 0; c < pars.ncol(); ++c) par_row[static_cast<size_t>(c)] = p_all_acc_for_trial(k, c);
					double cdf = cdf1(D, par_row.data(), model_context_for_funcs);
					if (!std::isfinite(cdf) || cdf >= 1.0) { logB = R_NegInf; break; }
					if (cdf <= 0.0) continue;
					if (cdf > 1.0 - 1e-15) cdf = 1.0 - 1e-15;
					logB += std::log1p(-cdf);
				}
				current_ll_val = log_sum_exp(logA, logB);
			
				} else {
					lower_for_trial = UC[start_row_idx];
					upper_for_trial = UT[start_row_idx];

					// Fast path: when the winner is unknown and there is no upper truncation (UT = +Inf),
					// the observation {rt == Inf} corresponds to the event {T_min > UC}. This can be
					// computed exactly as a product of survivor functions, avoiding numerical integration.
					// For LBAIO (posdrift=false) this already includes intrinsic-omission mass.
					if (R_j_idx == NA_INTEGER && upper_for_trial == R_PosInf) {
						double logS = 0.0;
						std::vector<double> par_row(static_cast<size_t>(pars.ncol()));
						for (int k = 0; k < n_lR_j; ++k) {
							if (!isok_for_trial[k]) { logS = R_NegInf; break; }
							for (int c = 0; c < pars.ncol(); ++c) par_row[static_cast<size_t>(c)] = p_all_acc_for_trial(k, c);
							double cdf = cdf1(lower_for_trial, par_row.data(), model_context_for_funcs);
							if (!std::isfinite(cdf)) { logS = R_NegInf; break; }
							if (cdf <= 0.0) continue;                   // log(1)
							if (cdf >= 1.0) { logS = R_NegInf; break; } // log(0)
							if (cdf > 1.0 - 1e-15) cdf = 1.0 - 1e-15;
							logS += std::log1p(-cdf);
						}
						current_ll_val = logS;
						slow_censor_closed_form = true;
					} else if (n_lR_j == 1) {
						// Simplified case for single accumulator
						if (!isok_for_trial[0]) {
							current_ll_val = R_NegInf;
						} else {
							std::vector<double> par_row(static_cast<size_t>(pars.ncol()));
							for (int c = 0; c < pars.ncol(); ++c) par_row[static_cast<size_t>(c)] = p_all_acc_for_trial(0, c);
							double cdf = cdf1(lower_for_trial, par_row.data(), model_context_for_funcs);
							if (!std::isfinite(cdf) || cdf >= 1.0) {
								current_ll_val = R_NegInf;
							} else if (cdf <= 0.0) {
								current_ll_val = 0.0; // log(1)
							} else {
								if (cdf > 1.0 - 1e-15) cdf = 1.0 - 1e-15;
								current_ll_val = std::log1p(-cdf);
							}
						}
						if (R_j_idx == NA_INTEGER && upper_for_trial == R_PosInf) {
							slow_censor_closed_form = true;
						}
					} else {
						if (R_j_idx != NA_INTEGER) { // Response (winner) is known
							current_ll_val = integrate_cached(R_j_idx, lower_for_trial, upper_for_trial);
						} else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
                  double log_l_k = integrate_cached(k_win, lower_for_trial, upper_for_trial);
                  current_ll_val =  log_sum_exp(current_ll_val, log_l_k);
                }
              }
				    }
			  }
        // Case 4: Missing RT (NA). Probability is sum of integral from LT to LC and UC to UT (i.e., outside the observation window but within truncation).
      } else if (NumericVector::is_na(rt_j)) {
        double lower_for_trial1 = LT[start_row_idx];
        double upper_for_trial1 = LC[start_row_idx];			
        double lower_for_trial2 = UC[start_row_idx];
        double upper_for_trial2 = UT[start_row_idx];
        if (R_j_idx != NA_INTEGER) { // Response (winner) is known
            double ll_L = integrate_cached(R_j_idx, lower_for_trial1, upper_for_trial1);
            double ll_U = integrate_cached(R_j_idx, lower_for_trial2, upper_for_trial2);
            current_ll_val = log_sum_exp(ll_L, ll_U);
        } else { // Response (winner) is unknown; sum probabilities over all possible winners
            for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
                double ll_L_k = integrate_cached(k_win, lower_for_trial1, upper_for_trial1);
                double ll_U_k = integrate_cached(k_win, lower_for_trial2, upper_for_trial2);
                double ll_k_sum = log_sum_exp(ll_L_k, ll_U_k);
                current_ll_val = log_sum_exp(current_ll_val,ll_k_sum); 
            }
        }
      }
		lower_for_trial = LT[start_row_idx];
		upper_for_trial = UT[start_row_idx];
    // Else: rt_j is 0 or negative finite (but not -Inf), current_prob_val remains 0.
		// Apply trial-level truncation-correction i.e. probability of observing ANY trial in the truncation window
		if (current_ll_val > min_ll && (LT[start_row_idx]!=0 || UT[start_row_idx] != R_PosInf) ) {
			std::string key = make_logZ_key(start_row_idx, lower_for_trial, upper_for_trial, pars, n_lR_j, pars.ncol());
			auto hit = log_integral_cache.find(key);
			if (hit != log_integral_cache.end()) { // fast path check if we've integrated this exact parameter set already and re-use result if we have
				log_Z_this = hit->second;
			} else { // otherwise compute the integral and cache it
				log_Z_this  = get_trunc_normaliser_cpp(p_all_acc_for_trial,
				                                       pdf1,
				                                       cdf1,
				                                       isok_for_trial,
				                                       lower_for_trial,
				                                       upper_for_trial,
				                                       n_lR_j,
				                                       pars.ncol(),
				                                       integration_epsilon,
				                                       model_context_for_funcs,
				                                       workspace.get());
				log_integral_cache.emplace(std::move(key), log_Z_this);
			}
			
			if (!NumericVector::is_na(log_Z_this) && R_FINITE(log_Z_this)){current_ll_val -= log_Z_this;}
		}
		
		current_ll_val = std::max(min_ll, current_ll_val); // Ensure probability is not negative
        ll_unique[unique_trial_idx] = current_ll_val;
        ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]); // Ensure not less than min_ll
    }
    // --- Summation of log-likelihoods for all unique trials ---
    double total_ll = 0;
    if (expand.length() > 0) { // If an expansion vector is provided (e.g. from non-compressed dadm)
		if (use_pC) { // multiply all likelihoods by the appropriate pC adjustment if it's present
            for (int j = 0; j < n_unique_trials; ++j) {
                double pC = pC_values[j];
                double log1m_pC = std::log1p(-pC); // log(1-pC)
                int start_row_idx = j * n_lR;
                double rt_j = rts_dadm[start_row_idx];

                if (R_FINITE(rt_j)) {
                    ll_unique[j] = log1m_pC + ll_unique[j]; // pRT * 1-pC
                } else {
                    double term1 = std::log(pC);
                    double term2 = log1m_pC + ll_unique[j];
                    ll_unique[j] = log_sum_exp(term1, term2); // pC * pIO * pCens
                }
            }
        }
        total_ll = sum(c_expand(ll_unique,expand));
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

namespace {

enum class SsGoDist { EXG, RDM };

inline int find_col_idx(const Rcpp::CharacterVector& cols, const std::string& name) {
  for (int i = 0; i < cols.size(); ++i) {
    if (Rcpp::as<std::string>(cols[i]) == name) return i;
  }
  return -1;
}

inline double safe_log1m(double p) {
  if (!std::isfinite(p)) return R_NegInf;
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return R_NegInf;
  if (p > 1.0 - 1e-15) p = 1.0 - 1e-15;
  return std::log1p(-p);
}

inline double clamp_cdf01(double cdf) {
  if (!std::isfinite(cdf)) return 1.0; // treat as finished (survivor=0)
  if (cdf <= 0.0) return 0.0;
  if (cdf >= 1.0) return 1.0;
  return cdf;
}

struct gsl_ss_stopwin_params {
  const double* go_pars; // row-major: n_go * n_go_par
  int n_go;
  int n_go_par;
  RaceCdf1Fun go_cdf1;
  void* go_ctx;
  double SSD;
  double muS;
  double sigmaS;
  double tauS;
};

double gsl_f_ss_stopwin(double s, void* p) {
  auto* P = static_cast<gsl_ss_stopwin_params*>(p);
  const double f_stop = ss_exgauss_pdf(s, P->muS, P->sigmaS, P->tauS);
  if (!(f_stop > 0.0) || !std::isfinite(f_stop)) return 0.0;

  const double t_abs = s + P->SSD;
  double logS = 0.0;
  for (int k = 0; k < P->n_go; ++k) {
    const double* par_k = P->go_pars + static_cast<size_t>(k) * P->n_go_par;
    double cdf = P->go_cdf1(t_abs, par_k, P->go_ctx);
    cdf = clamp_cdf01(cdf);
    const double ll = safe_log1m(cdf);
    if (!std::isfinite(ll)) return 0.0;
    logS += ll;
  }
  const double surv_prod = std::exp(logS);
  if (!(surv_prod > 0.0) || !std::isfinite(surv_prod)) return 0.0;
  return f_stop * surv_prod;
}

struct gsl_ss_stoptrunc_params {
  const double* go_pars; // row-major: n_go * n_go_par
  int n_go;
  int n_go_par;
  RacePdf1Fun go_pdf1;
  RaceCdf1Fun go_cdf1;
  void* go_ctx;
  double tf;
  double SSD;
  double muS;
  double sigmaS;
  double tauS;
};

double gsl_f_ss_stoptrunc(double t, void* p) {
  auto* P = static_cast<gsl_ss_stoptrunc_params*>(p);
  if (!(t > 0.0) || !std::isfinite(t)) return 0.0;

  double S_all = 1.0;
  double sum_pdf_over_S = 0.0;
  for (int k = 0; k < P->n_go; ++k) {
    const double* par_k = P->go_pars + static_cast<size_t>(k) * P->n_go_par;
    double cdf = P->go_cdf1(t, par_k, P->go_ctx);
    cdf = clamp_cdf01(cdf);
    if (cdf >= 1.0) return 0.0;
    const double S_k = 1.0 - cdf;
    if (!(S_k > 0.0) || !std::isfinite(S_k)) return 0.0;
    S_all *= S_k;
    if (!(S_all > 0.0) || !std::isfinite(S_all)) return 0.0;

    double pdf = P->go_pdf1(t, par_k, P->go_ctx);
    if (!(pdf > 0.0) || !std::isfinite(pdf)) continue;
    sum_pdf_over_S += pdf / S_k;
  }
  const double f_min_go = S_all * sum_pdf_over_S;
  if (!(f_min_go > 0.0) || !std::isfinite(f_min_go)) return 0.0;

  double S_stop = ss_stop_surv_abs(t, P->SSD, P->muS, P->sigmaS, P->tauS);
  if (!std::isfinite(S_stop)) S_stop = 0.0;
  S_stop = ss_clamp_prob01(S_stop);

  const double mix = P->tf + (1.0 - P->tf) * S_stop;
  if (!(mix > 0.0) || !std::isfinite(mix)) return 0.0;
  return f_min_go * mix;
}

double integrate_ss_stopwin_prob(double upper_duration,
                                 const double* go_pars,
                                 int n_go,
                                 int n_go_par,
                                 RaceCdf1Fun go_cdf1,
                                 void* go_ctx,
                                 double SSD,
                                 double muS,
                                 double sigmaS,
                                 double tauS,
                                 double epsilon,
                                 gsl_integration_workspace* w) {
  if (w == nullptr) Rcpp::stop("integrate_ss_stopwin_prob: GSL workspace is null.");
  if (n_go <= 0) return 0.0;
  if (!R_finite(upper_duration)) upper_duration = R_PosInf;
  if (upper_duration <= 0.0) return 0.0;

  gsl_function F;
  gsl_ss_stopwin_params P;
  P.go_pars = go_pars;
  P.n_go = n_go;
  P.n_go_par = n_go_par;
  P.go_cdf1 = go_cdf1;
  P.go_ctx = go_ctx;
  P.SSD = SSD;
  P.muS = muS;
  P.sigmaS = sigmaS;
  P.tauS = tauS;
  F.function = &gsl_f_ss_stopwin;
  F.params = &P;

  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status;
  double result = 0.0;
  double error = 0.0;

  if (upper_duration == R_PosInf) {
    status = gsl_integration_qagiu(&F, 0.0, 0, epsilon, 1000, w, &result, &error);
  } else {
    status = gsl_integration_qags(&F, 0.0, upper_duration, 0, epsilon, 1000, w, &result, &error);
  }

  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS) return 0.0;
  if (!(result > 0.0) || !std::isfinite(result)) return 0.0;
  if (result > 1.0) result = 1.0;
  return result;
}

double integrate_ss_stoptrunc_Z(double LT,
                                double UT,
                                const double* go_pars,
                                int n_go,
                                int n_go_par,
                                RacePdf1Fun go_pdf1,
                                RaceCdf1Fun go_cdf1,
                                void* go_ctx,
                                double tf,
                                double SSD,
                                double muS,
                                double sigmaS,
                                double tauS,
                                double epsilon,
                                gsl_integration_workspace* w) {
  if (w == nullptr) Rcpp::stop("integrate_ss_stoptrunc_Z: GSL workspace is null.");
  if (n_go <= 0) return 0.0;
  if (LT >= UT && !(LT == 0.0 && UT == R_PosInf)) return 0.0;

  gsl_function F;
  gsl_ss_stoptrunc_params P;
  P.go_pars = go_pars;
  P.n_go = n_go;
  P.n_go_par = n_go_par;
  P.go_pdf1 = go_pdf1;
  P.go_cdf1 = go_cdf1;
  P.go_ctx = go_ctx;
  P.tf = tf;
  P.SSD = SSD;
  P.muS = muS;
  P.sigmaS = sigmaS;
  P.tauS = tauS;
  F.function = &gsl_f_ss_stoptrunc;
  F.params = &P;

  double low = LT;
  double upp = UT;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status;
  double result = 0.0;
  double error = 0.0;
  if (upp == R_PosInf) {
    if (low < 0.0) low = 0.0;
    if (low >= R_PosInf) {
      result = 0.0;
      status = GSL_SUCCESS;
    } else {
      status = gsl_integration_qagiu(&F, low, 0, epsilon, 1000, w, &result, &error);
    }
  } else {
    status = gsl_integration_qags(&F, low, upp, 0, epsilon, 1000, w, &result, &error);
  }
  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS) return 0.0;
  if (!(result > 0.0) || !std::isfinite(result)) return 0.0;
  return result;
}

} // namespace

double c_log_likelihood_ss_cens_trunc(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::DataFrame& dadm,
    int n_trials,
    const Rcpp::LogicalVector& winner,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& isok,
    int n_lR,
    const std::string& ss_type) {

  if (n_trials == 0) return 0.0;
  if (n_lR <= 0) Rcpp::stop("c_log_likelihood_ss_cens_trunc: n_lR must be positive.");
  if (n_trials % n_lR != 0) Rcpp::stop("c_log_likelihood_ss_cens_trunc: dadm nrows not a multiple of n_lR.");

  SsGoDist go_dist;
  if (ss_type == "SSexG") {
    go_dist = SsGoDist::EXG;
  } else if (ss_type == "SShybrid") {
    go_dist = SsGoDist::RDM;
  } else {
    Rcpp::stop("Unsupported stop-signal model type string in calc_ll: " + ss_type);
  }

  const Rcpp::CharacterVector pcols = Rcpp::colnames(pars);
  const int col_gf = find_col_idx(pcols, "gf");
  const int col_tf = find_col_idx(pcols, "tf");
  const int col_muS = find_col_idx(pcols, "muS");
  const int col_sigmaS = find_col_idx(pcols, "sigmaS");
  const int col_tauS = find_col_idx(pcols, "tauS");
  if (col_gf < 0 || col_tf < 0 || col_muS < 0 || col_sigmaS < 0 || col_tauS < 0) {
    Rcpp::stop("c_log_likelihood_ss_cens_trunc: required stop-signal parameter columns not found in pars.");
  }

  std::array<int, 5> go_cols{};
  int n_go_par = 0;
  RacePdf1Fun go_pdf1 = nullptr;
  RaceCdf1Fun go_cdf1 = nullptr;
  void* go_ctx = nullptr;

  if (go_dist == SsGoDist::EXG) {
    const int col_mu = find_col_idx(pcols, "mu");
    const int col_sigma = find_col_idx(pcols, "sigma");
    const int col_tau = find_col_idx(pcols, "tau");
    if (col_mu < 0 || col_sigma < 0 || col_tau < 0) {
      Rcpp::stop("c_log_likelihood_ss_cens_trunc: SSexG requires go columns mu/sigma/tau.");
    }
    go_cols = {col_mu, col_sigma, col_tau, -1, -1};
    n_go_par = 3;
    go_pdf1 = &ss_exg_go_pdf1;
    go_cdf1 = &ss_exg_go_cdf1;
  } else {
    const int col_v = find_col_idx(pcols, "v");
    const int col_B = find_col_idx(pcols, "B");
    const int col_A = find_col_idx(pcols, "A");
    const int col_t0 = find_col_idx(pcols, "t0");
    const int col_s = find_col_idx(pcols, "s");
    if (col_v < 0 || col_B < 0 || col_A < 0 || col_t0 < 0 || col_s < 0) {
      Rcpp::stop("c_log_likelihood_ss_cens_trunc: SShybrid requires go columns v/B/A/t0/s.");
    }
    go_cols = {col_v, col_B, col_A, col_t0, col_s};
    n_go_par = 5;
    go_pdf1 = &drdm_scalar;
    go_cdf1 = &prdm_scalar;
  }

  if (!dadm.containsElementNamed("SSD")) {
    Rcpp::stop("c_log_likelihood_ss_cens_trunc: dadm must contain SSD column.");
  }
  const Rcpp::NumericVector RT = dadm["rt"];
  const Rcpp::NumericVector SSDv = dadm["SSD"];
  const bool has_LT = dadm.containsElementNamed("LT");
  const bool has_UT = dadm.containsElementNamed("UT");
  const Rcpp::NumericVector LT = has_LT ? Rcpp::NumericVector(dadm["LT"]) : Rcpp::NumericVector();
  const Rcpp::NumericVector UT = has_UT ? Rcpp::NumericVector(dadm["UT"]) : Rcpp::NumericVector();
  const bool has_UC = dadm.containsElementNamed("UC");
  const Rcpp::NumericVector UC = has_UC ? Rcpp::NumericVector(dadm["UC"]) : Rcpp::NumericVector();

  // Determine which response accumulators are GO (lI==2) vs stop-triggered (ST; lI==1).
  std::vector<int> go_offsets;
  std::vector<int> st_offsets;
  if (dadm.containsElementNamed("lI")) {
    const Rcpp::IntegerVector lI = dadm["lI"];
    for (int k = 0; k < n_lR; ++k) {
      const int code = lI[k];
      if (code == 2) go_offsets.push_back(k);
      else if (code == 1) st_offsets.push_back(k);
      else {
        Rcpp::stop("c_log_likelihood_ss_cens_trunc: unexpected lI code (expected 1=ST, 2=GO).");
      }
    }
  } else {
    // No lI provided: treat all accumulators as GO and none as ST.
    go_offsets.reserve(static_cast<size_t>(n_lR));
    for (int k = 0; k < n_lR; ++k) go_offsets.push_back(k);
  }
  const int n_go = static_cast<int>(go_offsets.size());
  const int n_st = static_cast<int>(st_offsets.size());

  auto workspace = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>(
      gsl_integration_workspace_alloc(1000),
      &gsl_integration_workspace_free);
  if (!workspace) Rcpp::stop("c_log_likelihood_ss_cens_trunc: Failed to allocate GSL integration workspace.");

  const double epsilon = 1e-7;
  std::unordered_map<std::string, double> cache_prob;
  std::unordered_map<std::string, double> cache_logZ;

  const int n_unique = n_trials / n_lR;
  Rcpp::NumericVector ll_unique(n_unique, min_ll);

  auto build_resp_pars = [&](int start_row, std::vector<double>& resp_pars_out) {
    resp_pars_out.resize(static_cast<size_t>(n_lR) * static_cast<size_t>(n_go_par));
    for (int k = 0; k < n_lR; ++k) {
      const int row = start_row + k;
      for (int j = 0; j < n_go_par; ++j) {
        resp_pars_out[static_cast<size_t>(k) * n_go_par + j] =
            pars(row, go_cols[static_cast<size_t>(j)]);
      }
    }
  };

  auto build_goonly_pars = [&](const std::vector<double>& resp_pars_in,
                              std::vector<double>& goonly_pars_out) {
    goonly_pars_out.resize(static_cast<size_t>(n_go) * static_cast<size_t>(n_go_par));
    for (int i = 0; i < n_go; ++i) {
      const int k0 = go_offsets[static_cast<size_t>(i)];
      for (int j = 0; j < n_go_par; ++j) {
        goonly_pars_out[static_cast<size_t>(i) * n_go_par + j] =
            resp_pars_in[static_cast<size_t>(k0) * n_go_par + j];
      }
    }
  };

  auto surv_prod = [&](double t, const std::vector<double>& pars_pack, int n_acc) -> double {
    double logS = 0.0;
    for (int k = 0; k < n_acc; ++k) {
      const double* par_k = pars_pack.data() + static_cast<size_t>(k) * n_go_par;
      double cdf = go_cdf1(t, par_k, go_ctx);
      cdf = clamp_cdf01(cdf);
      const double ll = safe_log1m(cdf);
      if (!std::isfinite(ll)) return 0.0;
      logS += ll;
    }
    const double out = std::exp(logS);
    if (!std::isfinite(out)) return 0.0;
    return out;
  };

  auto surv_prod_subset = [&](double t,
                             const std::vector<double>& resp_pars,
                             const std::vector<int>& offsets) -> double {
    double logS = 0.0;
    for (int off : offsets) {
      const double* par_k = resp_pars.data() + static_cast<size_t>(off) * n_go_par;
      double cdf = go_cdf1(t, par_k, go_ctx);
      cdf = clamp_cdf01(cdf);
      const double ll = safe_log1m(cdf);
      if (!std::isfinite(ll)) return 0.0;
      logS += ll;
    }
    const double out = std::exp(logS);
    if (!std::isfinite(out)) return 0.0;
    return out;
  };

  auto log_joint_for_winner_subset = [&](double t,
                                        const std::vector<double>& resp_pars,
                                        const std::vector<int>& offsets,
                                        int winner_off) -> double {
    const double* par_w = resp_pars.data() + static_cast<size_t>(winner_off) * n_go_par;
    double pdf_w = go_pdf1(t, par_w, go_ctx);
    if (!(pdf_w > 0.0) || !std::isfinite(pdf_w)) return R_NegInf;
    double logv = std::log(pdf_w);
    for (int off : offsets) {
      if (off == winner_off) continue;
      const double* par_k = resp_pars.data() + static_cast<size_t>(off) * n_go_par;
      double cdf = go_cdf1(t, par_k, go_ctx);
      cdf = clamp_cdf01(cdf);
      const double ll = safe_log1m(cdf);
      if (!std::isfinite(ll)) return R_NegInf;
      logv += ll;
    }
    return logv;
  };

  auto log_race_subset = [&](double t,
                             const std::vector<double>& resp_pars,
                             const std::vector<int>& offsets,
                             int observed_winner_idx0) -> double {
    if (offsets.empty()) return R_NegInf;
    if (observed_winner_idx0 >= 0) {
      // Only valid if the observed winner is in this subset.
      bool in_subset = false;
      for (int off : offsets) if (off == observed_winner_idx0) { in_subset = true; break; }
      if (!in_subset) return R_NegInf;
      return log_joint_for_winner_subset(t, resp_pars, offsets, observed_winner_idx0);
    }
    // Unknown winner: sum over subset.
    double out = R_NegInf;
    for (int off : offsets) {
      out = log_sum_exp(out, log_joint_for_winner_subset(t, resp_pars, offsets, off));
    }
    return out;
  };

  // Cache keys using only the go-only parameter pack plus stop parameters.
  auto make_stopwin_key = [&](double upper_dur,
                             double SSD,
                             double muS,
                             double sigmaS,
                             double tauS,
                             const std::vector<double>& goonly_pars) -> std::string {
    std::ostringstream oss;
    oss.setf(std::ios::scientific);
    oss << std::setprecision(17);
    oss << "SSP|" << upper_dur << '|' << SSD << '|' << muS << '|' << sigmaS << '|' << tauS << '|';
    for (double v : goonly_pars) {
      if (Rcpp::NumericVector::is_na(v)) oss << "NA";
      else oss << v;
      oss << ',';
    }
    return oss.str();
  };

  auto make_stoptrunc_key = [&](double LTj,
                               double UTj,
                               double tf,
                               double SSD,
                               double muS,
                               double sigmaS,
                               double tauS,
                               const std::vector<double>& goonly_pars) -> std::string {
    std::ostringstream oss;
    oss.setf(std::ios::scientific);
    oss << std::setprecision(17);
    oss << "SSZ|" << LTj << '|' << UTj << '|' << tf << '|' << SSD << '|' << muS << '|' << sigmaS << '|' << tauS << '|';
    for (double v : goonly_pars) {
      if (Rcpp::NumericVector::is_na(v)) oss << "NA";
      else oss << v;
      oss << ',';
    }
    return oss.str();
  };

  for (int j = 0; j < n_unique; ++j) {
    const int start_row = j * n_lR;
    if (!isok[start_row]) {
      ll_unique[j] = min_ll;
      continue;
    }

    const double rt_j = RT[start_row];
    const double SSD = SSDv[start_row];
    const bool is_stop = R_finite(SSD);
    const double gf = pars(start_row, col_gf);
    const double tf = pars(start_row, col_tf);
    const double muS = pars(start_row, col_muS);
    const double sigmaS = pars(start_row, col_sigmaS);
    const double tauS = pars(start_row, col_tauS);

    std::vector<double> resp_pars;
    build_resp_pars(start_row, resp_pars);
    std::vector<double> goonly_pars;
    build_goonly_pars(resp_pars, goonly_pars);

    if (rt_j == R_PosInf) {
      if (!has_UC) Rcpp::stop("c_log_likelihood_ss_cens_trunc: rt contains Inf but dadm has no UC column.");
      const double UCj = UC[start_row];
      const double S_go_UC = (n_go > 0) ? surv_prod(UCj, goonly_pars, n_go) : 1.0;

      double p_nr = 0.0;
      if (!is_stop) {
        p_nr = gf + (1.0 - gf) * S_go_UC;
      } else {
        const double S_stop_UC = ss_clamp_prob01(ss_stop_surv_abs(UCj, SSD, muS, sigmaS, tauS));
        const double upper_dur = UCj - SSD;

        double pStop_UC = 0.0;
        if (n_go > 0 && (1.0 - gf) > 0.0) {
          const std::string key = make_stopwin_key(upper_dur, SSD, muS, sigmaS, tauS, goonly_pars);
          auto hit = cache_prob.find(key);
          if (hit != cache_prob.end()) {
            pStop_UC = hit->second;
          } else {
            pStop_UC = integrate_ss_stopwin_prob(upper_dur,
                                                 goonly_pars.data(), n_go, n_go_par,
                                                 go_cdf1, go_ctx,
                                                 SSD, muS, sigmaS, tauS,
                                                 epsilon, workspace.get());
            cache_prob.emplace(key, pStop_UC);
          }
        } else if (n_go == 0) {
          pStop_UC = 1.0;
        }
        pStop_UC = ss_clamp_prob01(pStop_UC);

        // ST survivor by UC if stop triggers (deadline precedes trigger => ST cannot respond yet).
        double S_st_UC = 1.0;
        if (n_st > 0) {
          double t_eff = UCj - SSD;
          if (!R_finite(t_eff) || t_eff <= 0.0) t_eff = 0.0;
          S_st_UC = surv_prod_subset(t_eff, resp_pars, st_offsets);
        }

        const double p_nr_if_tf = gf + (1.0 - gf) * S_go_UC;
        const double p_nr_if_trig = S_st_UC * (gf + (1.0 - gf) * (pStop_UC + S_go_UC * S_stop_UC));
        p_nr = tf * p_nr_if_tf + (1.0 - tf) * p_nr_if_trig;
      }

      ll_unique[j] = (p_nr > 0.0 && std::isfinite(p_nr)) ? std::max(min_ll, std::log(p_nr)) : min_ll;
      continue;
    }

    if (Rcpp::NumericVector::is_na(rt_j)) {
      double p_nr = 0.0;
      if (!is_stop) {
        p_nr = gf;
      } else {
        if (n_st > 0) {
          // With stop-triggered accumulators, a true intrinsic non-response requires both go failure and trigger failure.
          p_nr = gf * tf;
        } else {
          const std::string key = make_stopwin_key(R_PosInf, SSD, muS, sigmaS, tauS, goonly_pars);
          auto hit = cache_prob.find(key);
          double pStop = 0.0;
          if (hit != cache_prob.end()) {
            pStop = hit->second;
          } else {
            pStop = integrate_ss_stopwin_prob(R_PosInf,
                                              goonly_pars.data(), n_go, n_go_par,
                                              go_cdf1, go_ctx,
                                              SSD, muS, sigmaS, tauS,
                                              epsilon, workspace.get());
            cache_prob.emplace(key, pStop);
          }
          pStop = ss_clamp_prob01(pStop);
          p_nr = gf + (1.0 - gf) * (1.0 - tf) * pStop;
        }
      }
      ll_unique[j] = (p_nr > 0.0 && std::isfinite(p_nr)) ? std::max(min_ll, std::log(p_nr)) : min_ll;
      continue;
    }

    if (!R_FINITE(rt_j) || !(rt_j > 0.0)) {
      ll_unique[j] = min_ll;
      continue;
    }

    const double LTj = has_LT ? LT[start_row] : 0.0;
    const double UTj = has_UT ? UT[start_row] : R_PosInf;
    if (rt_j < LTj || rt_j > UTj) {
      ll_unique[j] = min_ll;
      continue;
    }

    int observed_winner_idx0 = -1;
    int n_true = 0;
    for (int k = 0; k < n_lR; ++k) {
      if (winner[start_row + k] == TRUE) {
        n_true++;
        observed_winner_idx0 = k;
      }
    }
    if (n_true != 1) observed_winner_idx0 = -1; // unknown

    double log_like = R_NegInf;
    if (!is_stop) {
      // GO trial: only GO accumulators race.
      const double log_race_go = log_race_subset(rt_j, resp_pars, go_offsets, observed_winner_idx0);
      if (std::isfinite(log_race_go) && (1.0 - gf) > 0.0) {
        log_like = std::log1p(-gf) + log_race_go;
      }
    } else {
      // Stop trial: GO responses are a mixture over trigger failure; ST responses require triggering.
      bool winner_is_go = false;
      bool winner_is_st = false;
      if (observed_winner_idx0 >= 0) {
        for (int off : go_offsets) if (off == observed_winner_idx0) { winner_is_go = true; break; }
        for (int off : st_offsets) if (off == observed_winner_idx0) { winner_is_st = true; break; }
      }

      // If we can't classify the winner, fall back to summing GO+ST response masses.
      const bool unknown_winner = (observed_winner_idx0 < 0);

      // GO response mass (only GO racers, at absolute rt).
      double log_go_mass = R_NegInf;
      if (n_go > 0 && (unknown_winner || winner_is_go)) {
        const double log_race_go = log_race_subset(rt_j, resp_pars, go_offsets, unknown_winner ? -1 : observed_winner_idx0);
        if (std::isfinite(log_race_go) && (1.0 - gf) > 0.0) {
          const double log_base = std::log1p(-gf) + log_race_go;
          // tf branch: stop doesn't trigger => GO only
          if (tf > 0.0) log_go_mass = log_sum_exp(log_go_mass, std::log(tf) + log_base);
          // triggered branch: stop triggers => require stop survivor and ST survivors
          if ((1.0 - tf) > 0.0) {
            const double S_stop_rt = ss_clamp_prob01(ss_stop_surv_abs(rt_j, SSD, muS, sigmaS, tauS));
            double logS_stop = (S_stop_rt > 0.0) ? std::log(S_stop_rt) : R_NegInf;
            double logS_st = 0.0;
            if (n_st > 0) {
              const double t_eff = rt_j - SSD;
              const double S_st = surv_prod_subset(t_eff, resp_pars, st_offsets);
              logS_st = (S_st > 0.0) ? std::log(S_st) : R_NegInf;
            }
            if (std::isfinite(logS_stop) && std::isfinite(logS_st)) {
              log_go_mass = log_sum_exp(log_go_mass, std::log1p(-tf) + log_base + logS_stop + logS_st);
            }
          }
        }
      }

      // ST response mass (only if stop triggers).
      double log_st_mass = R_NegInf;
      if (n_st > 0 && (unknown_winner || winner_is_st) && (1.0 - tf) > 0.0) {
        const double t_eff = rt_j - SSD;
        const double log_race_st = log_race_subset(t_eff, resp_pars, st_offsets, unknown_winner ? -1 : observed_winner_idx0);
        if (std::isfinite(log_race_st)) {
          // pStop up to rt (duration upper = rt-SSD), depends only on GO racers.
          double pStop_rt = 0.0;
          if (n_go == 0) {
            pStop_rt = 1.0;
          } else if ((1.0 - gf) > 0.0) {
            const double upper_dur = rt_j - SSD;
            const std::string key = make_stopwin_key(upper_dur, SSD, muS, sigmaS, tauS, goonly_pars);
            auto hit = cache_prob.find(key);
            if (hit != cache_prob.end()) {
              pStop_rt = hit->second;
            } else {
              pStop_rt = integrate_ss_stopwin_prob(upper_dur,
                                                   goonly_pars.data(), n_go, n_go_par,
                                                   go_cdf1, go_ctx,
                                                   SSD, muS, sigmaS, tauS,
                                                   epsilon, workspace.get());
              cache_prob.emplace(key, pStop_rt);
            }
          }
          pStop_rt = ss_clamp_prob01(pStop_rt);
          const double log_pStop = (pStop_rt > 0.0) ? std::log(pStop_rt) : R_NegInf;
          const double log1m_pStop = safe_log1m(pStop_rt);

          // GO survivor by rt (only matters when stop doesn't beat go).
          double logS_go_rt = 0.0;
          if (n_go > 0) {
            const double S_go_rt = surv_prod(rt_j, goonly_pars, n_go);
            logS_go_rt = (S_go_rt > 0.0) ? std::log(S_go_rt) : R_NegInf;
          }

          // gf branch: only ST races
          if (gf > 0.0) {
            log_st_mass = log_sum_exp(log_st_mass, std::log(gf) + log_race_st);
          }
          // (1-gf) branch: mixture over pStop vs go surviving
          if ((1.0 - gf) > 0.0) {
            double log_mix = log_pStop;
            if (std::isfinite(log1m_pStop) && std::isfinite(logS_go_rt)) {
              log_mix = log_sum_exp(log_mix, log1m_pStop + logS_go_rt);
            }
            if (std::isfinite(log_mix)) {
              log_st_mass = log_sum_exp(log_st_mass, std::log1p(-gf) + log_race_st + log_mix);
            }
          }
          if (std::isfinite(log_st_mass)) {
            log_st_mass += std::log1p(-tf);
          }
        }
      }

      log_like = log_sum_exp(log_go_mass, log_st_mass);
    }

    if (LTj != 0.0 || UTj != R_PosInf) {
      double logZ = R_NegInf;
      if (!is_stop) {
        const double S_lt = (n_go > 0) ? surv_prod(LTj, goonly_pars, n_go) : 1.0;
        const double S_ut = (UTj == R_PosInf) ? 0.0 : ((n_go > 0) ? surv_prod(UTj, goonly_pars, n_go) : 0.0);
        const double diff = std::max(S_lt - S_ut, std::numeric_limits<double>::min());
        const double Z = (1.0 - gf) * diff;
        logZ = (Z > 0.0) ? std::log(Z) : R_NegInf;
      } else {
        if (n_st > 0) {
          Rcpp::stop("c_log_likelihood_ss_cens_trunc: truncation normalisation is not yet supported when stop-triggered accumulators are present.");
        }
        const std::string key = make_stoptrunc_key(LTj, UTj, tf, SSD, muS, sigmaS, tauS, goonly_pars);
        auto hit = cache_logZ.find(key);
        if (hit != cache_logZ.end()) {
          logZ = hit->second;
        } else {
          const double Zcore = integrate_ss_stoptrunc_Z(LTj, UTj,
                                                       goonly_pars.data(), n_go, n_go_par,
                                                       go_pdf1, go_cdf1, go_ctx,
                                                       tf, SSD, muS, sigmaS, tauS,
                                                       epsilon, workspace.get());
          const double Z = (1.0 - gf) * std::max(Zcore, std::numeric_limits<double>::min());
          logZ = (Z > 0.0) ? std::log(Z) : R_NegInf;
          cache_logZ.emplace(key, logZ);
        }
      }

      if (!std::isfinite(logZ)) {
        ll_unique[j] = min_ll;
        continue;
      }
      log_like -= logZ;
    }

    ll_unique[j] = std::isfinite(log_like) ? std::max(min_ll, log_like) : min_ll;
  }

  if (expand.length() > 0) {
    return sum(c_expand(ll_unique, expand));
  }
  double total = 0.0;
  for (int i = 0; i < ll_unique.size(); ++i) total += ll_unique[i];
  return total;
}

#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "trend.h"
#include "model_race_cens_trunc.h"

// ---- START: Context Structs and Static Adapters for Model Functions ----

// Context struct for LBA & RDM (can be shared if context is similar)
struct ContextForRaceModels {
    double min_lik_for_pdf;
    bool use_posdrift; // Added for LBA models, true by default or for non-LBA models
};
// For LNR, context might be simpler or could reuse above if only min_lik_for_pdf is needed.

// Static adapter for LBA dfun
static Rcpp::NumericVector lba_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            bool log_p,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    Rcpp::LogicalVector idx(pars.nrow(), true); // LBA default
    // Pass use_posdrift from context to dlba_c
    return dlba_c(rt, pars, idx, ctx->min_lik_for_pdf, is_ok, ctx->use_posdrift);
}

// Static adapter for LBA pfun
static Rcpp::NumericVector lba_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            bool log_p,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    Rcpp::LogicalVector idx(pars.nrow(), true); // LBA default
    // Pass use_posdrift from context to plba_c
    return plba_c(rt, pars, idx, ctx->min_lik_for_pdf, is_ok, ctx->use_posdrift);
}

// Static adapter for RDM dfun
static Rcpp::NumericVector rdm_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            bool log_p,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    Rcpp::LogicalVector idx(pars.nrow(), true); // RDM default
    return drdm_c(rt, pars, idx, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for RDM pfun
static Rcpp::NumericVector rdm_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            bool log_p,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    Rcpp::LogicalVector idx(pars.nrow(), true); // RDM default
    return prdm_c(rt, pars, idx, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LNR dfun
static Rcpp::NumericVector lnr_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            bool log_p,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context); // Reusing same context struct
    Rcpp::LogicalVector posdrift_ignored; // LNR's dlnr_c doesn't use posdrift
    return dlnr_c(rt, pars, posdrift_ignored, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LNR pfun
static Rcpp::NumericVector lnr_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            bool log_p,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context); // Reusing same context struct
    Rcpp::LogicalVector posdrift_ignored;
    return plnr_c(rt, pars, posdrift_ignored, ctx->min_lik_for_pdf, is_ok);
}

// ---- END: Context Structs and Static Adapters ----

using namespace Rcpp;

LogicalVector c_do_bound(NumericMatrix pars,
                              const std::vector<BoundSpec>& specs)
{
  int nrows = pars.nrow();
  LogicalVector result(nrows, true);

  // For each parameter that has bounds
  for (size_t j = 0; j < specs.size(); j++) {
    const BoundSpec& bs = specs[j];
    int col_idx   = bs.col_idx;
    double min_v  = bs.min_val;
    double max_v  = bs.max_val;
    bool has_exc  = bs.has_exception;
    double exc_val= bs.exception_val;

    // Check each row
    for (int i = 0; i < nrows; i++) {
      double val = pars(i, col_idx);
      bool ok = (val > min_v && val < max_v);
      if (!ok && has_exc) {
        // If out of range, see if exception matches
        ok = (val == exc_val);
      }
      // Merge with existing result (like result = result & ok_col)
      if (result[i] && !ok) {
        result[i] = false;
      }
    }
  }
  return result;
}

NumericVector c_do_pre_transform(NumericVector p_vector,
                                 const std::vector<PreTransformSpec>& specs)
{
  for (size_t i = 0; i < specs.size(); i++) {
    const PreTransformSpec& s = specs[i];
    double val = p_vector[s.index];

    switch (s.code) {
    case PTF_EXP: {
      // lower + exp(real)
      p_vector[s.index] = s.lower + std::exp(val);
      break;
    }
    case PTF_PNORM: {
      double range = s.upper - s.lower;
      // lower + range * Φ(real)
      p_vector[s.index] = s.lower +
        range * R::pnorm(val, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
      break;
    }
    default:
      // no transform
      break;
    }
  }
  return p_vector;
}

NumericMatrix c_do_transform(NumericMatrix pars,
                             const std::vector<TransformSpec>& specs)
{
  int nrow = pars.nrow();

  for (size_t j = 0; j < specs.size(); j++) {
    const TransformSpec& sp = specs[j];
    int          col_idx = sp.col_idx;
    TransformCode c      = sp.code;
    double        lw     = sp.lower;
    double        up     = sp.upper;

    switch (c) {
    case EXP: {
      for (int i = 0; i < nrow; i++) {
      // lower + exp(real)
      pars(i, col_idx) = lw + std::exp(pars(i, col_idx));
    }
      break;
    }
    case PNORM: {
      double range = up - lw;
      for (int i = 0; i < nrow; i++) {
        // lower + range * Φ(real)
        pars(i, col_idx) = lw +
          range * R::pnorm(pars(i, col_idx), 0.0, 1.0,
                           /*lower_tail=*/1, /*log_p=*/0);
      }
      break;
    }
    case IDENTITY:
    default:
      // do nothing
      break;
    }
  }
  return pars;
}


NumericMatrix c_map_p(NumericVector p_vector,
                      CharacterVector p_types,
                      List designs,
                      int n_trials,
                      DataFrame data,
                      List trend,
                      List transforms) {

  // Extract information about trends
  bool has_trend = (trend.length() > 0); // or another condition
  bool premap = false;
  bool pretransform = false;
  CharacterVector trend_names;
  // If trend has these flags
  if (has_trend) {
    premap = trend.attr("premap");
    pretransform = trend.attr("pretransform");
    trend_names = trend.names();
  }
  NumericVector p_mult_design;
  int n_params = p_types.size();
  NumericMatrix pars(n_trials, n_params);
  colnames(pars) = p_types;
  NumericMatrix trend_pars;
  // Identify trend parameters if any
  CharacterVector trend_pnames;
  LogicalVector trend_index(n_params, FALSE);
  if (has_trend && (premap || pretransform)) {
    // first loop over trends to get all trend pnames
    // But only for trends that are premap or pretransform
    for(unsigned int q = 0; q < trend.length(); q++){
      List cur_trend = trend[q];
      trend_pnames = c_add_charvectors(trend_pnames, as<CharacterVector>(cur_trend["trend_pnames"]));
      // Takes care of shared parameters
      trend_pnames = unique(trend_pnames);
    }
    // index which p_types are trends
    LogicalVector trend_index = contains_multiple(p_types,trend_pnames);
    for(unsigned int j = 0; j < trend_index.length(); j ++){
      // If we are a trend parameter:
      if(trend_index[j] == TRUE){
        NumericMatrix cur_design_trend = designs[j];
        CharacterVector cur_names_trend = colnames(cur_design_trend);
        // Take the current design and loop over columns
        // Multiply by design matrix
        for(int k = 0; k < cur_design_trend.ncol(); k ++){
          String cur_name_trend(cur_names_trend[k]);
          p_mult_design =  p_vector[cur_name_trend] * cur_design_trend(_, k);
          p_mult_design[is_nan(p_mult_design)] = 0;
          pars(_, j) = pars(_, j) + p_mult_design;
        }
      }
    }
    trend_pars = submat_rcpp_col_by_names(pars, trend_pnames);
    std::vector<TransformSpec> t_specs = make_transform_specs(trend_pars, transforms);
    trend_pars = c_do_transform(trend_pars, t_specs);
    trend_index = contains_multiple(p_types, trend_pnames);
  }
  for(int i = 0, t = 0; i < n_params; i++){
    if(trend_index[i] == FALSE){
      NumericMatrix cur_design = designs[i];
      CharacterVector cur_names = colnames(cur_design);

      for(int j = 0; j < cur_design.ncol(); j ++){
        String cur_name(cur_names[j]);
        NumericVector p_mult_design(n_trials, p_vector[cur_name]);
        // at this point we're multiplying by specific parameters (e.g. v_lMd)
        // So first apply trend to this parameter, then multiply by design matrix;
        if(has_trend && premap){
          // Check if trend is on current parameter
          LogicalVector cur_has_trend = contains(trend_names, cur_name);
          // This is a bit tricky and arguable.
          // Here we first fill a p_mult_design vector, then apply a trend then multiply with design matrix
          // Arguably you could also multiply parameter with design matrix and then apply trend
          // But that results in weird effects that if a parameter is set at 0, it could no longer be 0 post-trend
          for(unsigned int w = 0; w < cur_has_trend.length(); w ++){
            if(cur_has_trend[w] == TRUE){ // if so apply trend
              List cur_trend = trend[cur_name];
              CharacterVector cur_trend_pnames = cur_trend["trend_pnames"];
              p_mult_design = run_trend_rcpp(data, cur_trend, p_mult_design,
                                             submat_rcpp_col_by_names(trend_pars,cur_trend_pnames));

            }
          }
        }
        p_mult_design = p_mult_design * cur_design(_, j);
        p_mult_design[is_nan(p_mult_design)] = 0;
        pars(_, i) = pars(_, i) + p_mult_design;
      };
    } else if(pretransform){
      // These trends aren't applied here, but rather after mapping,
      // But they are transformed here already, so input them here.
      pars(_, i) = trend_pars(_, t);
      t++;
    }
  };
  if(has_trend && premap){
    pars = submat_rcpp_col(pars, !contains_multiple(p_types, trend_pnames));
  }
  return(pars);
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, List transforms, const std::vector<PreTransformSpec>& p_specs,
                              CharacterVector p_types, List designs, int n_trials, DataFrame data, List trend){
  bool has_trend = (trend.length() > 0);
  bool pretransform = false;
  bool posttransform = false;
  // If trend has these flags
  if (has_trend) {
    pretransform = trend.attr("premap"); // Note: R code used trend.attr("premap") for pretransform flag. Keeping consistent.
    posttransform = trend.attr("posttransform");
  }
  NumericVector p_vector_updtd(clone(p_vector));
  CharacterVector par_names = p_vector_updtd.names();
  p_vector_updtd = c_do_pre_transform(p_vector_updtd, p_specs);
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, transforms);
  // // Check if pretransform trend applies
  if(pretransform){ // automatically only applies if trend
    pars = prep_trend(data, trend, pars);
  }
  std::vector<TransformSpec> t_specs = make_transform_specs(pars, transforms);
  pars = c_do_transform(pars, t_specs);
  // Check if posttransform trend applies
  if(posttransform){ // automatically only applies if trend
    pars = prep_trend(data, trend, pars);
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

// This is the old standard race likelihood - will be commented out after refactor
double c_log_likelihood_race(NumericMatrix pars, DataFrame data,
                             NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector, bool), // Added bool for posdrift
                             NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector, bool), // Added bool for posdrift
                             const int n_trials, LogicalVector winner, IntegerVector expand,
                             double min_ll, LogicalVector is_ok, bool use_posdrift_arg = true){ // Added use_posdrift_arg
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  CharacterVector R = data["R"]; // Not directly used for R factor levels, but for n_acc via levels(dadm$R) in R
  NumericVector lR = data["lR"]; // Factor indicating accumulator index
  NumericVector lds_exp(n_out);
  const int n_acc = unique(lR).length();

  if(sum(contains(data.names(), "RACE")) == 1){ // Handling for subset of accumulators
    NumericVector NACC = data["RACE"];
    CharacterVector vals_NACC = NACC.attr("levels");
    for(int x = 0; x < pars.nrow(); x++){
      if(lR[x] > atoi(vals_NACC[NACC[x]-1])){ // -1 because NACC values are 1-based R factor levels
        pars(x,0) = NA_REAL; // Mark parameters for non-existent accumulators as NA
      }
    }
  }
  // For LBA, dfun/pfun are dlba_c/plba_c which now take use_posdrift_arg
  // For RDM/LNR, their dfun/pfun (e.g. drdm_c) don't take use_posdrift_arg.
  // This path needs to be smarter or dfun/pfun need a uniform signature (e.g. via adapters or optional args)
  // For now, assuming this function is primarily for LBA if posdrift is relevant.
  // This whole function will be commented out.
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll), is_ok, use_posdrift_arg));
  lds[winner] = win;
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok, use_posdrift_arg));
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll; // Avoid -Inf from log(0) if pfun result was min_ll for CDF
    lds[!winner] = loss;
  }
  lds[is_na(lds)] = min_ll; // Catch any other NAs

  if(n_acc > 1){
    NumericVector ll_out = lds[winner]; // Likelihoods for winning accumulators
    NumericVector lds_los = lds[!winner]; // Likelihoods for losing accumulators
    if(n_acc == 2){
      ll_out = ll_out + lds_los; // Simple sum if only one loser
    } else{
      // Summing log-likelihoods of (n_acc-1) losers for each winning accumulator's trial
      for(int z = 0; z < ll_out.length(); z++){ // ll_out.length() is n_unique_trials if winner is per unique trial
        ll_out[z] = ll_out[z] + sum(lds_los[Rcpp::Range( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }

    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    ll_out = c_expand(ll_out, expand); // Decompress to full trial count
    return(sum(ll_out));
  } else{ // Single accumulator case
    lds_exp = c_expand(lds, expand); // Decompress
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    return(sum(lds_exp));
  }
}

// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
            List designs, String type_rcpp, List bounds, List transforms, List pretransforms,
            CharacterVector p_types, double min_ll, List trend, SEXP model_r_object = R_NilValue){
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow(); // Total rows in dadm (n_unique_trials * n_accumulators)
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  p_vector.names() = p_names;
  NumericMatrix pars(n_trials, p_types.length()); // Will be resized per particle if n_trials changes (it shouldn't here)
  LogicalVector is_ok(n_trials);                 // Same as above

  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<PreTransformSpec> p_specs;
  std::vector<BoundSpec> bound_specs;
  std::string type_std(type_rcpp);

  if(type_std == "DDM"){
    IntegerVector expand = data.attr("expand");
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, Rcpp::_);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type_std == "MRI" || type_std == "MRI_AR1"){
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, Rcpp::_);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      if(type_std == "MRI"){
        lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      } else{
        lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
    }
  } else { // This block now handles ALL race models (LBA, RDM, LNR) using the _CENS_TRUNC machinery
    IntegerVector expand = data.attr("expand");
    std::string base_model_type = type_std;
    bool is_cens_trunc_explicitly = false;
    if (type_std.length() > 11 && type_std.substr(type_std.length() - 11) == "_CENS_TRUNC") {
        base_model_type = type_std.substr(0, type_std.find("_CENS_TRUNC"));
        is_cens_trunc_explicitly = true; // Flag if it was originally a _CENS_TRUNC type
    }

    // Check if the base_model_type is one of the known race models
    if (base_model_type != "LBA" && base_model_type != "RDM" && base_model_type != "LNR") {
        Rcpp::stop("Unsupported model type in race model block: " + type_std);
    }

    RacePdfFun model_dfun_ptr = nullptr;
    RaceCdfFun model_pfun_ptr = nullptr;
    ContextForRaceModels current_model_ctx;
    current_model_ctx.min_lik_for_pdf = exp(min_ll);
    current_model_ctx.use_posdrift = true; // Default to true, override for LBA if specified in model object

    if (base_model_type == "LBA") {
        model_dfun_ptr = &lba_dfun_adapter;
        model_pfun_ptr = &lba_pfun_adapter;
        if (model_r_object != R_NilValue) {
            Rcpp::List model_props(model_r_object);
            if (model_props.containsElementNamed("posdrift")) {
                current_model_ctx.use_posdrift = Rcpp::as<bool>(model_props["posdrift"]);
            } else {
                // Rcpp::Rcout << "Warning: LBA model from R object does not contain 'posdrift' element. Defaulting to true." << std::endl;
            }
        } else if (is_cens_trunc_explicitly || type_std == "LBA") { // Only warn if we expected to find it
             // Rcpp::Rcout << "Warning: R model object not passed to calc_ll for LBA model type '" << type_std << "'. Defaulting posdrift to true." << std::endl;
        }
    } else if (base_model_type == "RDM") {
        model_dfun_ptr = &rdm_dfun_adapter;
        model_pfun_ptr = &rdm_pfun_adapter;
    } else if (base_model_type == "LNR") {
        model_dfun_ptr = &lnr_dfun_adapter;
        model_pfun_ptr = &lnr_pfun_adapter;
    }
    // No else needed here due to the check above

    int n_acc = 0;
    if (data.containsElementNamed("lR")) {
        Rcpp::NumericVector lR_col = data["lR"]; // Use NumericVector for lR as in c_log_likelihood_race
        if (lR_col.length() > 0) {
            Rcpp::NumericVector lR_unique_vals = Rcpp::unique(lR_col);
            n_acc = lR_unique_vals.size();
        }
    }
    if (n_acc == 0 && data.nrows() > 0) { // Fallback or error if n_acc couldn't be determined
         Rcpp::Rcerr << "Warning: calc_ll could not determine n_acc from data[\"lR\"] for model type " << type_std
                  << ". This may lead to errors or incorrect likelihoods." << std::endl;
         // Setting n_acc to 1 if completely undetermined and data exists, as a last resort, though likely wrong.
         // A better approach might be to stop if n_acc is vital and cannot be found.
         // For now, this matches the old fallback logic's spirit if it were to use data.nrows() / unique_conditions.
         // If n_trials is total rows in dadm, and dadm has n_unique_cond * n_acc rows,
         // then n_acc = n_trials / n_unique_cond. This is hard to get here without n_unique_cond.
         // The c_log_likelihood_race_cens_trunc will re-calculate n_acc based on data["lR"] or stop.
         // So, we can pass a placeholder or let it be determined inside.
         // Let's ensure n_acc is determined before the loop or c_log_likelihood_race_cens_trunc will handle it.
         // The c_log_likelihood_race_cens_trunc has its own n_acc determination.
    }
    // n_acc determination is primarily handled inside c_log_likelihood_race_cens_trunc now.
    // We pass a placeholder n_acc (e.g. 0 or a guess) if needed, or rely on its internal logic.
    // The n_acc argument to c_log_likelihood_race_cens_trunc is critical.
    // Let's re-evaluate n_acc for the call:
    if (data.containsElementNamed("lR")) {
        Rcpp::NumericVector lR_col_for_nacc = data["lR"];
        if (lR_col_for_nacc.length() > 0) {
            n_acc = Rcpp::unique(lR_col_for_nacc).size();
        }
    }
    if (n_acc == 0 && data.nrows() > 0) { // If still zero
        Rcpp::stop("calc_ll: Failed to determine n_acc for race model from 'lR' column.");
    }


    for (int i = 0; i < n_particles; ++i) {
        p_vector = p_matrix(i, Rcpp::_);
        if (i == 0) {
            p_specs = make_pretransform_specs(p_vector, pretransforms);
        }
        pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);

        is_ok.resize(pars.nrow());
        if (i == 0 || bound_specs.empty()) {
            bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
        }
        is_ok = c_do_bound(pars, bound_specs);

        lls[i] = c_log_likelihood_race_cens_trunc(pars, data,
                                                  model_dfun_ptr, model_pfun_ptr,
                                                  min_ll, is_ok, n_acc, expand, // Pass determined n_acc
                                                  &current_model_ctx);
    }

    /* // Old standard race model block - now commented out
    } else { // This was the original "else" for standard LBA/RDM/LNR
        LogicalVector winner = data["winner"];
        // The function pointers for dfun/pfun in the old c_log_likelihood_race did not take a context or a separate posdrift.
        // Their C++ implementations (dlba_c etc.) would need to be adapted or wrapped if we were to keep this path.
        // For LBA, posdrift would need to be passed.
        bool model_use_posdrift_std = true;
        if (type_std == "LBA" && model_r_object != R_NilValue) {
            Rcpp::List model_props(model_r_object);
            if (model_props.containsElementNamed("posdrift")) {
                model_use_posdrift_std = Rcpp::as<bool>(model_props["posdrift"]);
            }
        }

        // Original function pointers for c_log_likelihood_race
        // These would need to be updated to match the new signatures if this block were active
        // e.g., dlba_c now needs a posdrift argument.
        // NumericVector (*dfun_std)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector, bool);
        // NumericVector (*pfun_std)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector, bool);

        // This block is complex to adapt due to differing function signatures for LBA vs RDM/LNR regarding posdrift.
        // The unified path using adapters and context is cleaner.
        // ... (original logic for selecting dfun/pfun and calling c_log_likelihood_race) ...
        // For example, for LBA:
        // if(type_std == "LBA"){
        //   lls[i] = c_log_likelihood_race(pars, data, dlba_c, plba_c, n_trials, winner, expand, min_ll, is_ok, model_use_posdrift_std);
        // } else if ...
        // This demonstrates why unifying is simpler.
    }
    */
  }
  return(lls);
}


// Batched version of f_race_integrand_cpp for finite RT trials
// rts_batch: vector of RTs for the batch
// pars_all_trials_ordered: A List of NumericMatrix, where each element is an ordered parameter matrix
//                          for a trial (winner first), matching the order in rts_batch.
// model_dfun, model_pfun: function pointers to the model's density and CDF
// n_acc: number of accumulators
// model_specific_context: context pointer for model functions
Rcpp::NumericVector f_race_integrand_batch_cpp(
    const Rcpp::NumericVector& rts_batch,
    const Rcpp::List& pars_all_trials_ordered, // List of NumericMatrix, each n_acc x n_params
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    void* model_specific_context) {

    int n_batch_trials = rts_batch.size();
    if (n_batch_trials == 0) {
        return Rcpp::NumericVector(0);
    }
    if (pars_all_trials_ordered.size() != n_batch_trials) {
        Rcpp::stop("f_race_integrand_batch_cpp: rts_batch size and pars_all_trials_ordered size mismatch.");
    }

    Rcpp::NumericVector results(n_batch_trials);

    for (int i = 0; i < n_batch_trials; ++i) {
        double t = rts_batch[i];
        Rcpp::NumericMatrix p_trial_this_winner_first = Rcpp::as<Rcpp::NumericMatrix>(pars_all_trials_ordered[i]);

        if (p_trial_this_winner_first.nrow() != n_acc) {
             Rcpp::stop("f_race_integrand_batch_cpp: Parameter matrix for a trial has incorrect number of rows for n_acc.");
        }

        // Calculate PDF for the winner
        Rcpp::NumericMatrix p_winner = p_trial_this_winner_first.row(0);
        Rcpp::NumericVector pdf_winner_vec = model_dfun(Rcpp::NumericVector::create(t), p_winner, Rcpp::LogicalVector::create(true), false, model_specific_context);
        double pdf_winner = pdf_winner_vec[0];

        if (R_IsNA(pdf_winner) || !R_finite(pdf_winner) || pdf_winner < 0) {
            results[i] = 0.0;
            continue;
        }

        double survivor_losers = 1.0;
        if (n_acc > 1) {
            for (int k = 1; k < n_acc; ++k) { // Losers are in rows 1 to n_acc-1
                Rcpp::NumericMatrix p_loser_k = p_trial_this_winner_first.row(k);
                Rcpp::NumericVector cdf_loser_k_vec = model_pfun(Rcpp::NumericVector::create(t), p_loser_k, Rcpp::LogicalVector::create(true), false, model_specific_context);
                double cdf_loser_k = cdf_loser_k_vec[0];

                if (R_IsNA(cdf_loser_k) || !R_finite(cdf_loser_k)) {
                    survivor_losers = 0.0; // If any CDF is problematic, product becomes zero
                    break;
                }
                double s_loser_k = 1.0 - cdf_loser_k;
                s_loser_k = (s_loser_k < 0.0) ? 0.0 : s_loser_k; // Clamp at 0
                s_loser_k = (s_loser_k > 1.0) ? 1.0 : s_loser_k; // Clamp at 1 (should not be needed for 1-CDF if CDF is valid)
                survivor_losers *= s_loser_k;
                if (survivor_losers == 0.0) break; // Optimization
            }
        }

        double final_val = pdf_winner * survivor_losers;
        if (R_IsNA(final_val) || !R_finite(final_val) || final_val < 0) {
            results[i] = 0.0;
        } else {
            results[i] = final_val;
        }
    }
    return results;
}

// ---- START: New code for censored/truncated race models ----
// #include "model_race_cens_trunc.h" // Include the new header // MOVED TO TOP
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling

// Struct to pass parameters to GSL integrand
struct gsl_race_params {
    const Rcpp::NumericMatrix* p_trial_this_winner_first;
    RacePdfFun model_dfun;
    RaceCdfFun model_pfun;
    int n_acc;
    void* model_specific_context; // To pass context to model_dfun/model_pfun
};

// GSL-compatible adapter for f_race_integrand_cpp
double gsl_f_race_adapter(double t, void *p) {
    gsl_race_params* params = static_cast<gsl_race_params*>(p);

    if (t < 0) return 0.0; // RTs are non-negative

    Rcpp::NumericVector t_vec(1); // model_dfun/pfun expect NumericVector for rt
    t_vec[0] = t;

    const Rcpp::NumericMatrix& p_mat = *(params->p_trial_this_winner_first);

    Rcpp::NumericMatrix p_winner_mat(1, p_mat.ncol());
    p_winner_mat.row(0) = p_mat.row(0);

    Rcpp::LogicalVector ok_win(1, true); // Assuming parameters passed here are valid for this specific call context
    Rcpp::NumericVector pdf_winner_vec = params->model_dfun(t_vec, p_winner_mat, ok_win, false, params->model_specific_context);
    double pdf_winner = pdf_winner_vec[0];

    if (!R_finite(pdf_winner) || pdf_winner < 0) pdf_winner = 0.0;

    if (params->n_acc > 1) {
        double survivor_losers = 1.0;
        for (int i = 1; i < params->n_acc; ++i) {
            Rcpp::NumericMatrix p_loser_i_mat(1, p_mat.ncol());
            p_loser_i_mat.row(0) = p_mat.row(i);
            Rcpp::LogicalVector ok_los(1, true);
            Rcpp::NumericVector cdf_loser_i_vec = params->model_pfun(t_vec, p_loser_i_mat, ok_los, false, params->model_specific_context);
            double cdf_loser_i = cdf_loser_i_vec[0];

            double s_loser_i = 1.0 - cdf_loser_i;
            if (!R_finite(s_loser_i) || s_loser_i < 0) s_loser_i = 0.0;
            if (s_loser_i > 1.0) s_loser_i = 1.0; // Should not happen with valid CDF
            survivor_losers *= s_loser_i;
            if (survivor_losers == 0) break; // Optimization
        }
        pdf_winner *= survivor_losers;
    }
    return (R_finite(pdf_winner) && pdf_winner > 0) ? pdf_winner : 0.0;
}


// Helper to order parameters (winner first)
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    int k_idx, // 1-based index
    int n_acc) {

    if (k_idx < 1 || k_idx > n_acc) {
        Rcpp::stop("order_pars_for_winner_cpp: k_idx out of bounds.");
    }

    Rcpp::NumericMatrix ordered_pars(n_acc, p_all_acc.ncol());
    ordered_pars.row(0) = p_all_acc.row(k_idx - 1);

    int current_row = 1;
    for (int i = 0; i < n_acc; ++i) {
        if (i == (k_idx - 1)) continue;
        ordered_pars.row(current_row++) = p_all_acc.row(i);
    }
    return ordered_pars;
}


// Numerical integration helper using GSL
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    double low,
    double upp,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    double epsilon = 1e-7,
    void* model_specific_context = nullptr) {

    if (low >= upp && !(low == 0 && upp == R_PosInf)) return 0.0;
    if (k_winner_idx < 1 || k_winner_idx > n_acc) {
      // Rcpp::Rcerr << "Warning: Invalid k_winner_idx in integrate_for_kth_winner_cpp: " << k_winner_idx << std::endl;
      return 0.0;
    }
    Rcpp::NumericMatrix pars_ordered = order_pars_for_winner_cpp(p_all_acc, k_winner_idx, n_acc);

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function F;
    gsl_race_params params_struct;
    params_struct.p_trial_this_winner_first = &pars_ordered;
    params_struct.model_dfun = model_dfun;
    params_struct.model_pfun = model_pfun;
    params_struct.n_acc = n_acc;
    params_struct.model_specific_context = model_specific_context;

    F.function = &gsl_f_race_adapter;
    F.params = &params_struct;

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status;
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
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
        // Rcpp::Rcerr << "GSL integration failed: " << gsl_strerror(status) << std::endl;
        return 0.0;
    }
    return (result > 0 && R_finite(result)) ? result : 0.0;
}

// Truncation correction factor helper - uses GSL via integrate_for_kth_winner_cpp
double get_trunc_corr_factor_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    double LT, double UT,
    int n_acc,
    double epsilon = 1e-7,
    void* model_specific_context = nullptr) {

    if (!(LT > 0 || UT < R_PosInf)) return 1.0; // No truncation if LT=0 and UT=Inf
    if (k_winner_idx < 1 || k_winner_idx > n_acc) return NA_REAL; // Invalid winner

    double prob_untruncated = integrate_for_kth_winner_cpp(k_winner_idx, p_all_acc, 0, R_PosInf, model_dfun, model_pfun, n_acc, epsilon, model_specific_context);
    double prob_truncated_interval = integrate_for_kth_winner_cpp(k_winner_idx, p_all_acc, LT, UT, model_dfun, model_pfun, n_acc, epsilon, model_specific_context);

    if (prob_truncated_interval > 1e-12) {
        if (!R_finite(prob_untruncated) || prob_untruncated < 0) return NA_REAL;
        return prob_untruncated / prob_truncated_interval;
    } else {
        if (prob_untruncated > 1e-12) return NA_REAL; // Density exists but not in [LT,UT] -> infinite correction
        return 1.0; // Both are zero or negligible
    }
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
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector& ok_params,   // Parameter validity for each row in 'pars' matrix
    int n_acc,                              // Number of accumulators in the race (must be > 0 if data exists)
    const Rcpp::IntegerVector& expand_vec,  // Vector for expanding unique LLs to full trial count
    void* model_context_for_funcs           // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
) {
    // Fetch censoring and truncation values from dadm attributes, with defaults
    double LT = 0.0, UT = R_PosInf, LC = 0.0, UC = R_PosInf;
    if (dadm.hasAttribute("LT") && !Rf_isNull(dadm.attr("LT"))) LT = Rcpp::as<double>(dadm.attr("LT"));
    if (dadm.hasAttribute("UT") && !Rf_isNull(dadm.attr("UT"))) UT = Rcpp::as<double>(dadm.attr("UT"));
    if (dadm.hasAttribute("LC") && !Rf_isNull(dadm.attr("LC"))) LC = Rcpp::as<double>(dadm.attr("LC"));
    if (dadm.hasAttribute("UC") && !Rf_isNull(dadm.attr("UC"))) UC = Rcpp::as<double>(dadm.attr("UC"));

    double integration_epsilon = 1e-7; // Tolerance for GSL integration

    Rcpp::NumericVector rts_dadm = dadm["rt"];
    Rcpp::IntegerVector R_idxs_dadm = dadm["R"];

    int n_total_dadm_rows = dadm.nrows();
    if (n_total_dadm_rows == 0) return 0.0; // No data, no likelihood

    if (n_acc <= 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: n_acc must be positive and correctly determined before this call.");
    if (n_total_dadm_rows % n_acc != 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: dadm nrows not a multiple of n_acc.");

    int n_unique_trials = n_total_dadm_rows / n_acc;
    Rcpp::NumericVector ll_unique(n_unique_trials);
    ll_unique.fill(min_ll); // Initialize all to min_ll

    // Parameter matrix and validity vector checks
    if (pars.nrow() != n_total_dadm_rows) {
      Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_total_dadm_rows: " << n_total_dadm_rows << std::endl;
      Rcpp::stop("c_log_likelihood_race_cens_trunc: pars matrix dimensions do not match total dadm rows.");
    }
    if (ok_params.size() != pars.nrow()) {
       Rcpp::stop("c_log_likelihood_race_cens_trunc: ok_params size does not match pars matrix rows.");
    }

    std::vector<int> finite_rt_unique_trial_indices; // Stores original indices of unique trials for batching
    std::vector<int> other_unique_trial_indices;   // Stores original indices for iterative processing

    // Categorize unique trials based on RT properties and parameter validity
    for (int j = 0; j < n_unique_trials; ++j) {
        int first_row_in_dadm_for_trial = j * n_acc; // Starting row in dadm/pars for this unique trial

        // Check parameter validity for the entire block of accumulators for this unique trial
        bool params_ok_for_this_unique_trial = true;
        for (int k_acc = 0; k_acc < n_acc; ++k_acc) {
            if (!ok_params[first_row_in_dadm_for_trial + k_acc]) {
                params_ok_for_this_unique_trial = false;
                break;
            }
        }
        if (!params_ok_for_this_unique_trial) {
            // ll_unique[j] is already min_ll by initialization, so just continue
            continue;
        }

        double rt_j = rts_dadm[first_row_in_dadm_for_trial]; // RT for this unique trial
        int R_j_idx = R_idxs_dadm[first_row_in_dadm_for_trial]; // Winner index (1-based from R factor)

        // Criteria for batch processing: finite, positive RT, within truncation bounds [LT, UT], and known winner
        if (R_FINITE(rt_j) && rt_j > 0 && rt_j >= LT && rt_j <= UT && R_j_idx != NA_INTEGER) {
            finite_rt_unique_trial_indices.push_back(j);
        } else {
            other_unique_trial_indices.push_back(j); // All other cases (NA, Inf, -Inf, outside bounds, or unknown winner with finite RT)
        }
    }

    // --- Part 1: Batch process finite RT trials (Observed RTs within truncation bounds) ---
    // This section handles trials where RT is finite, positive, within truncation bounds (LT, UT),
    // and the response (winner) is known. These are suitable for batch calculation of densities.
    if (!finite_rt_unique_trial_indices.empty()) {
        Rcpp::NumericVector batch_rts(finite_rt_unique_trial_indices.size());
        Rcpp::List batch_pars_ordered(finite_rt_unique_trial_indices.size());
        std::vector<int> batch_R_idxs(finite_rt_unique_trial_indices.size());
        std::vector<Rcpp::NumericMatrix> batch_pars_unord(finite_rt_unique_trial_indices.size());

        // Collect data for batch processing
        for (size_t i = 0; i < finite_rt_unique_trial_indices.size(); ++i) {
            int unique_trial_idx = finite_rt_unique_trial_indices[i];
            int first_row_in_dadm = unique_trial_idx * n_acc;

            batch_rts[i] = rts_dadm[first_row_in_dadm];
            batch_R_idxs[i] = R_idxs_dadm[first_row_in_dadm];

            Rcpp::NumericMatrix p_all_acc_for_trial(n_acc, pars.ncol());
            for(int k_acc=0; k_acc < n_acc; ++k_acc) {
                p_all_acc_for_trial.row(k_acc) = pars.row(first_row_in_dadm + k_acc);
            }
            batch_pars_unord[i] = p_all_acc_for_trial;
            batch_pars_ordered[i] = order_pars_for_winner_cpp(p_all_acc_for_trial, batch_R_idxs[i], n_acc);
        }

        // Calculate probability densities for the entire batch of finite RT trials
        Rcpp::NumericVector batch_prob_densities = f_race_integrand_batch_cpp(
            batch_rts, batch_pars_ordered, model_dfun, model_pfun, n_acc, model_context_for_funcs
        );

        // Apply truncation correction and calculate log-likelihood for each trial in the batch
        for (size_t i = 0; i < finite_rt_unique_trial_indices.size(); ++i) {
            int unique_trial_idx = finite_rt_unique_trial_indices[i];
            double prob_density = batch_prob_densities[i];

            if (prob_density > std::numeric_limits<double>::epsilon()) {
                double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(
                    batch_R_idxs[i], batch_pars_unord[i], model_dfun, model_pfun,
                    LT, UT, n_acc, integration_epsilon, model_context_for_funcs
                );

                if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf <= 0) {
                    ll_unique[unique_trial_idx] = min_ll;
                } else {
                    ll_unique[unique_trial_idx] = std::log(prob_density) + std::log(trunc_cf);
                }
            } else {
                ll_unique[unique_trial_idx] = min_ll;
            }
             ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]);
        }
    }

    // --- Part 2: Process other trials (Infinite RTs, NA RTs, or finite RTs outside truncation) ---
    // These trials require individual processing, often involving numerical integration for censored intervals.
    for (int unique_trial_idx : other_unique_trial_indices) {
        // Parameters for these trials are already deemed valid if they are in other_unique_trial_indices
        // (because bad param trials were skipped in the categorization loop, ll_unique[j] remains min_ll for them).

        int first_row_in_dadm = unique_trial_idx * n_acc;
        Rcpp::NumericMatrix pars_condition_j_all_acc(n_acc, pars.ncol());
        for(int k_acc=0; k_acc < n_acc; ++k_acc) {
            pars_condition_j_all_acc.row(k_acc) = pars.row(first_row_in_dadm + k_acc);
        }

        double rt_j = rts_dadm[first_row_in_dadm];
        int R_j_idx = R_idxs_dadm[first_row_in_dadm];
        double current_prob_val = 0.0;

        // Case 1: Finite RT but strictly outside truncation bounds (LT, UT).
        // These trials have zero probability density according to the truncated distribution.
        if (R_FINITE(rt_j) && rt_j > 0 && (rt_j < LT || rt_j > UT)) {
             current_prob_val = 0;
        // Case 2: Fast censoring (RT = -Inf). Probability is integral from LT to LC.
        } else if (rt_j == R_NegInf) {
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                current_prob_val = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                if (current_prob_val > 0) { // Apply truncation correction if probability is non-zero
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon, model_context_for_funcs);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf <= 0) current_prob_val = 0; else current_prob_val *= trunc_cf;
                }
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                    if (p_k > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon, model_context_for_funcs);
                        if (!ISNAN(trunc_cf_k) && R_FINITE(trunc_cf_k) && trunc_cf_k >=0) current_prob_val += (p_k * trunc_cf_k);
                    }
                }
            }
        // Case 3: Slow censoring (RT = Inf). Probability is integral from UC to UT.
        } else if (rt_j == R_PosInf) {
             if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                current_prob_val = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                 if (current_prob_val > 0) { // Apply truncation correction
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon, model_context_for_funcs);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf <= 0) current_prob_val = 0; else current_prob_val *= trunc_cf;
                }
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                     if (p_k > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon, model_context_for_funcs);
                        if (!ISNAN(trunc_cf_k) && R_FINITE(trunc_cf_k) && trunc_cf_k >=0) current_prob_val += (p_k * trunc_cf_k);
                    }
                }
            }
        // Case 4: Missing RT (NA). Probability is sum of integral from LT to LC and UC to UT (i.e., outside the observation window but within truncation).
        } else if (ISNAN(rt_j)) {
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                double p_L = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                double p_U = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                current_prob_val = p_L + p_U;
                 if (current_prob_val > 0) { // Apply truncation correction
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon, model_context_for_funcs);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf <= 0) current_prob_val = 0; else current_prob_val *= trunc_cf;
                }
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_L_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                    double p_U_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                    double p_k_sum = p_L_k + p_U_k;
                    if (p_k_sum > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon, model_context_for_funcs);
                        if (!ISNAN(trunc_cf_k) && R_FINITE(trunc_cf_k) && trunc_cf_k >=0) current_prob_val += (p_k_sum * trunc_cf_k);
                    }
                }
            }
        }
        // Else: rt_j is 0 or negative finite (but not -Inf), current_prob_val remains 0.

        current_prob_val = std::max(0.0, current_prob_val); // Ensure probability is not negative
        ll_unique[unique_trial_idx] = (current_prob_val > std::numeric_limits<double>::epsilon()) ? std::log(current_prob_val) : min_ll;
        ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]); // Ensure not less than min_ll
    }

    // --- Summation of log-likelihoods for all unique trials ---
    double total_ll = 0;
    if (expand_vec.length() > 0) { // If an expansion vector is provided (e.g. from non-compressed dadm)
        for (int i = 0; i < expand_vec.length(); ++i) {
            int unique_idx = expand_vec[i] - 1; // expand_vec is 1-based from R
            if (unique_idx >= 0 && unique_idx < n_unique_trials) {
                total_ll += ll_unique[unique_idx];
            } else { // Should not happen with valid expand_vec
                total_ll += min_ll;
            }
        }
    } else if (dadm.containsElementNamed("N")) { // Alternative: if dadm has trial counts 'N' per unique condition
        Rcpp::NumericVector trial_Ns_dadm = dadm["N"];
        if (trial_Ns_dadm.length() == n_unique_trials) {
            for (int j = 0; j < n_unique_trials; ++j) {
                total_ll += ll_unique[j] * trial_Ns_dadm[j];
            }
        } else { // Fallback if N column size doesn't match (should not happen with correct dadm)
             for (int j = 0; j < n_unique_trials; ++j) {
                total_ll += ll_unique[j];
            }
        }
    } else { // Default: sum ll_unique directly (each unique trial counted once)
        for (int j = 0; j < n_unique_trials; ++j) {
            total_ll += ll_unique[j];
        }
    }
    return total_ll;
}


// ---- Rcpp Exported Wrapper for Testing ----
// This wrapper is for testing the c_log_likelihood_race_cens_trunc from R.
// It will use existing C++ pdf/cdf functions for a known model (e.g., LBA)
// The true integration will be calc_ll calling c_log_likelihood_race_cens_trunc directly
// with proper C++ function pointers and context derived from the R model object.
// [[Rcpp::export]]
Rcpp::NumericVector test_c_loglik_cens_trunc_wrapper_R(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame dadm,
    std::string model_type_str,    // e.g., "LBA", "RDM", "LNR" to select hardcoded C++ funcs
    double min_ll,
    Rcpp::LogicalVector ok_params,
    int n_acc,
    Rcpp::List R_model_obj_list // Pass the R model list for posdrift extraction
) {
    RacePdfFun model_dfun_ptr = nullptr;
    RaceCdfFun model_pfun_ptr = nullptr;
    void* current_model_context_ptr = nullptr;
    ContextForRaceModels test_model_ctx;

    test_model_ctx.min_lik_for_pdf = exp(min_ll);
    test_model_ctx.use_posdrift = true; // Default

    if (model_type_str == "LBA") { // Renamed from LBA_test for clarity
        model_dfun_ptr = &lba_dfun_adapter;
        model_pfun_ptr = &lba_pfun_adapter;
        if (R_model_obj_list.containsElementNamed("posdrift")) {
            test_model_ctx.use_posdrift = Rcpp::as<bool>(R_model_obj_list["posdrift"]);
        }
    } else if (model_type_str == "RDM") { // Renamed from RDM_test
        model_dfun_ptr = &rdm_dfun_adapter;
        model_pfun_ptr = &rdm_pfun_adapter;
        // RDM doesn't use posdrift, context.use_posdrift is ignored
    } else if (model_type_str == "LNR") { // Renamed from LNR_test
        model_dfun_ptr = &lnr_dfun_adapter;
        model_pfun_ptr = &lnr_pfun_adapter;
        // LNR doesn't use posdrift
    }
    else {
        Rcpp::warning("Model type '%s' not recognized for C++ test wrapper. Likelihood calculation will likely fail.", model_type_str.c_str());
        // Return sum of min_ll for all trials indicated by expand_vec or dadm structure
        Rcpp::IntegerVector expand_vec_dummy;
        if (dadm.hasAttribute("expand")) expand_vec_dummy = dadm.attr("expand");
        else expand_vec_dummy = Rcpp::seq(1, dadm.nrows() / (n_acc > 0 ? n_acc : 1) ); // Simplistic fallback for expand

        double total_min_ll = 0;
        if (expand_vec_dummy.length() > 0) {
             for (int i = 0; i < expand_vec_dummy.length(); ++i) total_min_ll += min_ll;
        } else if (dadm.nrows() > 0 && n_acc > 0) {
             total_min_ll = min_ll * (dadm.nrows()/n_acc) ;
        } else if (dadm.nrows() == 0) {
            total_min_ll = 0; // No data, no likelihood
        } else { // dadm.nrows > 0 but n_acc = 0 (error case)
            total_min_ll = R_NaN; // Indicate error
        }
        return Rcpp::NumericVector::create(total_min_ll);
    }
    current_model_context_ptr = &test_model_ctx;


    Rcpp::IntegerVector expand_vec_test;
    if (dadm.hasAttribute("expand")) {
      expand_vec_test = dadm.attr("expand");
    } else { // Create a simple 1-to-N unique trials expand_vec if not present
      if (n_acc > 0 && dadm.nrows() > 0) {
        int n_unique_for_test = dadm.nrows() / n_acc;
        expand_vec_test = Rcpp::seq(1, n_unique_for_test);
      }
      // if n_acc is 0 or no rows, expand_vec_test remains empty, handled by c_log_likelihood_race_cens_trunc
    }


    double result_ll = c_log_likelihood_race_cens_trunc(
        pars, dadm,
        model_dfun_ptr, model_pfun_ptr,
        min_ll, ok_params, n_acc,
        expand_vec_test,
        current_model_context_ptr
    );
    return Rcpp::NumericVector::create(result_ll);
}

// ---- END: New code for censored/truncated race models ----
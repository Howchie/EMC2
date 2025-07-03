#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "trend.h"
#include "utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include <string>

// For LNR, context might be simpler or could reuse above if only min_lik_for_pdf is needed.


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

// [[Rcpp::export]]
NumericMatrix get_pars_matrix_rcpp(NumericVector p_vector, NumericVector constants,
                                   List transforms, List pretransforms,
                                   CharacterVector p_types, List designs,
                                   int n_trials, DataFrame data, List trend) {
  std::vector<PreTransformSpec> p_specs = make_pretransform_specs(p_vector, pretransforms);
  return get_pars_matrix(p_vector, constants, transforms, p_specs,
                         p_types, designs, n_trials, data, trend);
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
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return(sum(lls_exp));
}

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
  const int n_acc = unique(lR).length();
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
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok)); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds)] = min_ll;

  if(n_acc > 1){
    // LogicalVector winner_exp = c_bool_expand(winner, expand);
    NumericVector ll_out = lds[winner];
    NumericVector lds_los = lds[!winner];
    if(n_acc == 2){
      ll_out = ll_out + lds_los;
    } else{
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
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
  p_vector.names() = p_names;
  NumericMatrix pars(n_trials, p_types.length());
  LogicalVector is_ok(n_trials);

  // Once (outside the main loop over particles):
  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<PreTransformSpec> p_specs;
  std::vector<BoundSpec> bound_specs;

  if(type == "DDM"){
    IntegerVector expand = data.attr("expand");
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
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
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
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
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    if(type == "LBA"){
      dfun = dlba_c;
      pfun = plba_c;
    } else if(type == "RDM"){
      dfun = drdm_c;
      pfun = prdm_c;
    } else if(type == "RDMSWTN"){
      dfun = drdmswtn_c;
      pfun = prdmswtn_c;
    } else{
      dfun = dlnr_c;
      pfun = plnr_c;
    }
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    for (int i = 0; i < n_particles; ++i) {
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
    }
  }
  return(lls);
}

/*
// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
            List designs, String type_rcpp, List bounds, List transforms, List pretransforms,
            CharacterVector p_types, double min_ll, List trend){ // SEXP model_r_object removed
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow(); // Total rows in dadm (n_unique_trials * n_accumulators)
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  
  CharacterVector p_names = colnames(p_matrix);
  NumericMatrix pars(n_trials, p_types.length());
  p_vector.names() = p_names;
  LogicalVector is_ok(n_trials);

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
  } else { // This block now handles ALL race models (LBA, RDM, LNR, and their variants)
    IntegerVector expand = data.attr("expand");
	LogicalVector winner = data["winner"];
    RacePdfFun model_dfun_ptr = nullptr;
    RaceCdfFun model_pfun_ptr = nullptr;
    ContextForRaceModels current_model_ctx;
    current_model_ctx.min_lik_for_pdf = exp(min_ll);
    current_model_ctx.use_posdrift = true; // Default: LBA uses posdrift, RDM/LNR ignore this context field.

    // Determine adapter functions and specific LBA posdrift setting based on type_std
    // The type_std string is expected to be e.g. "LBA", "LBA_IO", "RDM", "LNR".

    if (type_std.find("LBA") != std::string::npos) {
        model_dfun_ptr = &lba_dfun_adapter;
        model_pfun_ptr = &lba_pfun_adapter;
        // Check for the 'IO' (Implicit Omissions / no posdrift) flag in the original type_std
        if (type_std.find("IO") != std::string::npos) {
            current_model_ctx.use_posdrift = false;
        }
    } else if (type_std=="RDM") {
        model_dfun_ptr = &rdm_dfun_adapter;
        model_pfun_ptr = &rdm_pfun_adapter;
    } else if (type_std=="RDM_SWTN") {
        model_dfun_ptr = &rdmswtn_dfun_adapter;
        model_pfun_ptr = &rdmswtn_pfun_adapter;
    } else if (type_std.find("LNR") != std::string::npos) {
        model_dfun_ptr = &lnr_dfun_adapter;
        model_pfun_ptr = &lnr_pfun_adapter;
    } else {
        Rcpp::stop("Unsupported race model type string in calc_ll: " + type_std);
    }

    int n_acc = 0;
    if (data.containsElementNamed("lR")) {
        Rcpp::NumericVector lR_col_for_nacc = data["lR"];
        if (lR_col_for_nacc.length() > 0) {
            n_acc = Rcpp::unique(lR_col_for_nacc).size();
        }
    }
    if (n_acc == 0 && data.nrows() > 0) { // If still zero and there's data
        Rcpp::stop("calc_ll: Failed to determine n_acc for race model from 'lR' column for type: " + type_std);
    }
    if (n_acc == 0 && data.nrows() == 0) { // No data, n_acc can be 0. Loop over particles won't run if n_trials is 0.
        // This is fine, likelihood will be sum over 0 particles if n_trials is 0.
    }

    for (int i = 0; i < n_particles; ++i) {
        p_vector = p_matrix(i, Rcpp::_);
        if (i == 0) {
            p_specs = make_pretransform_specs(p_vector, pretransforms);
        }
        pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
        if (i == 0 || bound_specs.empty()) {
            bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
        }
        is_ok = c_do_bound(pars, bound_specs);
        is_ok = lr_all(is_ok, n_acc);
        lls[i] = c_log_likelihood_race_cens_trunc(pars, data,
                                                  model_dfun_ptr, model_pfun_ptr,n_trials,
                                                  winner, expand, min_ll, is_ok, n_acc, 
                                                  &current_model_ctx);
    }
  }
  return(lls);
}
*/


// GSL-compatible adapter for f_race_integrand_cpp
double gsl_f_race_adapter(double t, void *p) {
    gsl_race_params* params = static_cast<gsl_race_params*>(p);
	// ① get the params block
    // ② recover the model-specific context
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(params->model_specific_context);
	const int min_ll = ctx->min_lik_for_pdf;
    if (t < 0) return 0.0; // RTs are non-negative

    Rcpp::NumericVector t_vec(1); // model_dfun/pfun expect NumericVector for rt
    t_vec[0] = t;

    const Rcpp::NumericMatrix& p_mat = *(params->p_trial);
	const Rcpp::LogicalVector& winner = params->winner;
	const Rcpp::LogicalVector& isok = params->isok;
	NumericVector pp(p_mat.size());
	NumericVector pp_out(p_mat.size());
    NumericVector win = params->model_dfun(t_vec, p_mat, isok, winner, ctx);
	pp[winner]=win;
	double prod = 0.0;
    if (params->n_acc > 1) {
        Rcpp::NumericVector loss = 1-params->model_pfun(t_vec, p_mat, isok, !winner, ctx);
        loss[is_na(loss)] = std::exp(min_ll);
		pp[!winner] = loss;
    }
	for (int i = 0; i < p_mat.size(); ++i) {
        double term;
        term = pp[i];
        if (Rcpp::NumericVector::is_na(term) || term <= 0.0) {
            term = std::exp(min_ll);                  // floor for log-lik
		}
        prod *= term;
    }

    return prod;   
	
}


// Helper to order parameters (winner first)
/*
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    int k_idx, // 1-based index
    int n_acc) {
    if (k_idx < 1 || k_idx > n_acc) {
        Rcpp::stop("order_pars_for_winner_cpp: k_idx out of bounds.");
    }

    Rcpp::NumericMatrix ordered_pars(n_acc, p_all_acc.ncol());
	std::fill(ordered_pars.begin(), ordered_pars.end(), NA_REAL);
    ordered_pars.row(0) = p_all_acc.row(k_idx - 1);

    int current_row = 1;
    for (int i = 0; i < n_acc; ++i) {
        if (i == (k_idx - 1)) continue;
        ordered_pars.row(current_row++) = p_all_acc.row(i);
    }
    return ordered_pars;
} */


// Numerical integration helper using GSL
// TODO - reinstate ordering - this one shouldn't use the winner vec
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
	LogicalVector isok,
    double low,
    double upp,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    double epsilon = 1e-7,
    void* model_specific_context = nullptr) {

	LogicalVector winner(isok.size(), false);           // initialise all FALSE
    winner[k_winner_idx-1] = true; 
    if (low >= upp && !(low == 0 && upp == R_PosInf)) return 0.0;
    if (k_winner_idx < 1 || k_winner_idx > n_acc) {
      // Rcpp::Rcerr << "Warning: Invalid k_winner_idx in integrate_for_kth_winner_cpp: " << k_winner_idx << std::endl;
      return 0.0;
    }

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function F;
    gsl_race_params params_struct;
    params_struct.p_trial = &p_all_acc;
	params_struct.winner = winner;
	params_struct.isok = isok;
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

/* Deprecated
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
} */

// New common normaliser (ZH: I think the normalise by kth-winner is wrong?)
double get_trunc_normaliser_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
	LogicalVector isok_trial, 
    double LT, double UT,
    int n_acc,
    double epsilon,
    void* model_specific_context){
    double mass = 0.0;
    for(int k=1;k<=n_acc;++k)
        mass += integrate_for_kth_winner_cpp(k,p_all_acc, isok_trial, LT, UT, model_dfun, model_pfun, n_acc, epsilon, model_specific_context);
    return (mass>0.0) ? 1.0/mass : NA_REAL;
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
	const int n_trials,
	LogicalVector winner,
	Rcpp::IntegerVector expand,  // Vector for expanding unique LLs to full trial count
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector isok,   // Parameter validity for each row in 'pars' matrix
    int n_acc,                              // Number of accumulators in the race (must be > 0 if data exists)
    void* model_context_for_funcs           // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
) {
    // Fetch censoring and truncation values from dadm attributes, with defaults
    double LT = 0.0, UT = R_PosInf, LC = 0.0, UC = R_PosInf;
    if (dadm.hasAttribute("LT") && !Rf_isNull(dadm.attr("LT"))) LT = Rcpp::as<double>(dadm.attr("LT"));
    if (dadm.hasAttribute("UT") && !Rf_isNull(dadm.attr("UT"))) UT = Rcpp::as<double>(dadm.attr("UT"));
    if (dadm.hasAttribute("LC") && !Rf_isNull(dadm.attr("LC"))) LC = Rcpp::as<double>(dadm.attr("LC"));
    if (dadm.hasAttribute("UC") && !Rf_isNull(dadm.attr("UC"))) UC = Rcpp::as<double>(dadm.attr("UC"));

    double integration_epsilon = 1e-8; // Tolerance for GSL integration
	
	const int n_out = expand.length();
    Rcpp::NumericVector rts_dadm = dadm["rt"];
    Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
    Rcpp::IntegerVector lR_dadm = dadm["lR"];
	NumericVector lds(n_trials,min_ll); // initialise at min_ll
	NumericVector lds_exp(n_out);
	// If a RACE column exists, set parameters of accumulators not present on a
    // given trial to NA so the density functions return zero for them. This
    // mirrors logic from the old c_log_likelihood_race implementation.
	if (dadm.containsElementNamed("RACE")) {
		// factor codes (1-based) for each row
		Rcpp::IntegerVector race_idx = dadm["RACE"];
		// character levels (“2”, “3”, …)
		Rcpp::CharacterVector race_levels = race_idx.attr("levels");

		for (int row = 0; row < pars.nrow(); ++row) {

			// how many accumulators for this trial
			int n_acc_this_trial = std::stoi(
				Rcpp::as<std::string>(race_levels[ race_idx[row] - 1 ])
			);

			// lR_dadm is the (1-based) index of *this* accumulator on the trial
			if (lR_dadm[row] > n_acc_this_trial) {
				// accumulator not present → blank its parameter row
				std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
			}
		}
	}
    if (n_trials == 0) return 0.0; // No data, no likelihood

    if (n_acc <= 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: n_acc must be positive and correctly determined before this call.");
    if (n_trials % n_acc != 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: dadm nrows not a multiple of n_acc.");

    int n_unique_trials = n_trials / n_acc;
    Rcpp::NumericVector ll_unique(n_unique_trials);
    //ll_unique.fill(min_ll); // Initialize all to min_ll

    // Parameter matrix and validity vector checks
    if (pars.nrow() != n_trials) {
      Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_trials: " << n_trials << std::endl;
      Rcpp::stop("c_log_likelihood_race_cens_trunc: pars matrix dimensions do not match total dadm rows.");
    }
    if (isok.size() != pars.nrow()) {
       Rcpp::stop("c_log_likelihood_race_cens_trunc: isok size does not match pars matrix rows.");
    }

    Rcpp::LogicalVector finite_rt_mask(n_trials, false); // Mask for all dadm rows
	std::vector<int> finite_dadm_rows_indices; // Stores actual dadm/pars row indices for finite RTs
    finite_dadm_rows_indices.reserve(n_trials); // Reserve space
	// Common truncation factor, compute once -- seems WRONG
	double inv_Z = 1;  // Default inv_Z if no truncation
	std::vector<int> finite_rt_unique_trial_indices; // Stores *unique trial* indices for finite RTs
	std::vector<int> other_unique_trial_indices;   // Stores *unique trial* indices for other cases
    finite_rt_unique_trial_indices.reserve(n_unique_trials);
    other_unique_trial_indices.reserve(n_unique_trials);
    // Categorize unique trials based on RT properties and parameter validity

    for (int j = 0; j < n_unique_trials; ++j) {
        int first_row_in_dadm_for_trial = j * n_acc; // Starting row in dadm/pars for this unique trial
        double rt_j = rts_dadm[first_row_in_dadm_for_trial]; // RT for this unique trial
        int R_j_idx = R_idxs_dadm[first_row_in_dadm_for_trial]; // Winner index (1-based from R factor)
        // Criteria for batch processing: finite, positive RT, within truncation bounds [LT, UT], and known winner
        if (R_FINITE(rt_j) && rt_j > 0 && rt_j >= LT && rt_j <= UT && R_j_idx != NA_INTEGER) {
			finite_rt_unique_trial_indices.push_back(j); // Store unique trial index
            for(int k_acc = 0; k_acc < n_acc; ++k_acc) {
                int dadm_row_idx = first_row_in_dadm_for_trial + k_acc;
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
		Rcpp::NumericVector win = log(model_dfun(rts_for_finite_batch, finite_rt_pars, isok_for_finite_batch, winner_for_finite_batch, model_context_for_funcs));
	// Place results into the correct positions in 'lds'
        int current_lds_idx = 0;
        for(int i = 0; i < n_trials; ++i) {
            if(finite_rt_mask[i] && winner[i]) {
                lds[i] =win[current_lds_idx++];
            }
        }
		lds[winner & finite_rt_mask]=win;
		// Prepare survivor functions for losers, one matrix per accumulator
		if (n_acc > 1) {
			Rcpp::NumericVector loss = log(1-model_pfun(rts_for_finite_batch, finite_rt_pars, isok_for_finite_batch, !winner_for_finite_batch, model_context_for_funcs));
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
			int start_row_idx = unique_trial_idx*(n_acc);
            Rcpp::NumericMatrix p_all_acc_for_trial = pars(Rcpp::Range(start_row_idx, start_row_idx + n_acc - 1), Rcpp::_);
            Rcpp::LogicalVector isok_for_trial = isok[Rcpp::Range(start_row_idx, start_row_idx + n_acc - 1)];
			if (LT!=0 || UT != R_PosInf) {
						inv_Z = get_trunc_normaliser_cpp(p_all_acc_for_trial,model_dfun, model_pfun,isok_for_trial,LT, UT, n_acc,integration_epsilon,model_context_for_funcs);
			}
			double current_trial_ll_sum = 0.0;
			bool any_valid_component = false;
			for(int k=0; k < n_acc; ++k) {
				int dadm_row_idx = start_row_idx + k; // Correct index into lds
				if (lds[dadm_row_idx] > min_ll - 1.0) { // Check if component is not essentially min_ll
					any_valid_component = true;
				}
				current_trial_ll_sum += lds[dadm_row_idx];
			}

            if (!any_valid_component) { // All components were min_ll or less
                 ll_unique[unique_trial_idx] = min_ll;
			} 
			else if ((LT != 0 || UT != R_PosInf) && (NumericVector::is_na(inv_Z) || !R_FINITE(inv_Z) || inv_Z <= 0)) {
                    // If truncation active but inv_Z is bad, probability is effectively zero
                    ll_unique[unique_trial_idx] = min_ll;
			} 
			else if (LT != 0 || UT != R_PosInf) { // Truncation active and inv_Z is good
                        current_trial_ll_sum += std::log(inv_Z);
			}
			// If no truncation, inv_Z is not added.
			ll_unique[unique_trial_idx] = current_trial_ll_sum;
			ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]);
		}
	}

    // --- Process other trials (Infinite RTs, NA RTs, or finite RTs outside truncation) ---
    // These trials require individual processing, often involving numerical integration for censored intervals.
    for (size_t i = 0; i < other_unique_trial_indices.size(); ++i) {
		int unique_trial_idx = other_unique_trial_indices[i];
		int start_row_idx = unique_trial_idx*n_acc;
		Rcpp::NumericMatrix p_all_acc_for_trial = pars(Rcpp::Range(start_row_idx, start_row_idx + n_acc - 1), Rcpp::_);
        Rcpp::LogicalVector isok_for_trial = isok[Rcpp::Range(start_row_idx, start_row_idx + n_acc - 1)];
		Rcpp::LogicalVector winner_for_trial = winner[Rcpp::Range(start_row_idx, start_row_idx + n_acc - 1)];

        double rt_j = rts_dadm[start_row_idx];
        int R_j_idx = R_idxs_dadm[start_row_idx];
        double current_prob_val = 0.0;

        // Case 1: Finite RT but strictly outside truncation bounds (LT, UT).
        // These trials have zero probability density according to the truncated distribution. Shouldn't actually be possible but checking just in case.
        if (R_FINITE(rt_j) && rt_j > 0 && (rt_j < LT || rt_j > UT)) {
             current_prob_val = 0;
        // Case 2: Fast censoring (RT = -Inf). Probability is integral from LT to LC.
        } else if (rt_j == R_NegInf) {
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                current_prob_val = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
					current_prob_val *= p_k; 
					}
            }
        // Case 3: Slow censoring (RT = Inf). Probability is integral from UC to UT.
        } else if (rt_j == R_PosInf) {
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                current_prob_val = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
					current_prob_val *= p_k; 
                }
			}
        // Case 4: Missing RT (NA). Probability is sum of integral from LT to LC and UC to UT (i.e., outside the observation window but within truncation).
        } else if (NumericVector::is_na(rt_j)) {
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                double p_L = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                double p_U = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                current_prob_val = p_L + p_U;
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_L_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                    double p_U_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon, model_context_for_funcs);
                    double p_k_sum = p_L_k + p_U_k;
                    current_prob_val *= p_k_sum; 
                }
            }
        }
        // Else: rt_j is 0 or negative finite (but not -Inf), current_prob_val remains 0.
		// Apply trial-level truncation-correction i.e. probability of observing ANY trial in the truncation window
		if (current_prob_val > 0 && (LT!=0 || UT != R_PosInf) ) {
			inv_Z = get_trunc_normaliser_cpp(p_all_acc_for_trial,model_dfun, model_pfun,isok_for_trial,LT, UT, n_acc,integration_epsilon,model_context_for_funcs);
		}
        if (!NumericVector::is_na(inv_Z) && R_FINITE(inv_Z) && inv_Z >=0){current_prob_val += (current_prob_val * inv_Z);}
        
		current_prob_val = std::max(0.0, current_prob_val); // Ensure probability is not negative
        ll_unique[unique_trial_idx] = (current_prob_val > std::numeric_limits<double>::epsilon()) ? std::log(current_prob_val) : min_ll;
        ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]); // Ensure not less than min_ll
    }

    // --- Summation of log-likelihoods for all unique trials ---
    double total_ll = 0;
    if (expand.length() > 0) { // If an expansion vector is provided (e.g. from non-compressed dadm)
        total_ll = sum(c_expand(ll_unique,expand));
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

// ---- END: New code for censored/truncated race models ----
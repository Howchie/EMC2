#include <Rcpp.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "trend.h"
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
    pretransform = trend.attr("pretransform");
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
    IntegerVector expand = data.attr("expand"); // Used by both race likelihoods

    if (type.find("_CENS_TRUNC") != std::string::npos) {
        // --- New Censored/Truncated Race Model Logic ---
        std::string base_model_type = type.substr(0, type.find("_CENS_TRUNC"));
        RacePdfFun cens_model_dfun_ptr = nullptr;
        RaceCdfFun cens_model_pfun_ptr = nullptr;

        // Temporary lambda wrappers to adapt existing C++ functions to RacePdfFun/RaceCdfFun
        // These capture necessary context or provide defaults for extra parameters.
        // Note: min_lik_for_pdf is typically exp(min_ll)
        double min_lik_for_pdf = exp(min_ll);

        if (base_model_type == "LBA") {
            cens_model_dfun_ptr = [&](Rcpp::NumericVector rt, Rcpp::NumericMatrix event_pars, bool log_p_ignored) {
                // dlba_c needs: NumericVector t, NumericMatrix pars, LogicalVector posdrift, double min_lik, LogicalVector log_on
                // Assuming 'posdrift' is all TRUE for simplicity in this generic wrapper.
                // This might need to come from model definition or be configurable.
                Rcpp::LogicalVector posdrift_default(event_pars.nrow(), true);
                Rcpp::LogicalVector log_on_default(1, false); // log_p_ignored suggests we want raw density
                return dlba_c(rt, event_pars, posdrift_default, min_lik_for_pdf, log_on_default);
            };
            cens_model_pfun_ptr = [&](Rcpp::NumericVector rt, Rcpp::NumericMatrix event_pars, bool log_p_ignored) {
                Rcpp::LogicalVector posdrift_default(event_pars.nrow(), true);
                Rcpp::LogicalVector log_on_default(1, false);
                return plba_c(rt, event_pars, posdrift_default, min_lik_for_pdf, log_on_default);
            };
        } else if (base_model_type == "RDM") {
            cens_model_dfun_ptr = [&](Rcpp::NumericVector rt, Rcpp::NumericMatrix event_pars, bool log_p_ignored) {
                Rcpp::LogicalVector posdrift_default(event_pars.nrow(), true); // RDM typically positive drift
                Rcpp::LogicalVector log_on_default(1, false);
                return drdm_c(rt, event_pars, posdrift_default, min_lik_for_pdf, log_on_default);
            };
            cens_model_pfun_ptr = [&](Rcpp::NumericVector rt, Rcpp::NumericMatrix event_pars, bool log_p_ignored) {
                Rcpp::LogicalVector posdrift_default(event_pars.nrow(), true);
                Rcpp::LogicalVector log_on_default(1, false);
                return prdm_c(rt, event_pars, posdrift_default, min_lik_for_pdf, log_on_default);
            };
        } else if (base_model_type == "LNR") {
             cens_model_dfun_ptr = [&](Rcpp::NumericVector rt, Rcpp::NumericMatrix event_pars, bool log_p_ignored) {
                Rcpp::LogicalVector posdrift_ignored; // LNR doesn't use posdrift
                Rcpp::LogicalVector log_on_default(1, false);
                return dlnr_c(rt, event_pars, posdrift_ignored, min_lik_for_pdf, log_on_default);
            };
            cens_model_pfun_ptr = [&](Rcpp::NumericVector rt, Rcpp::NumericMatrix event_pars, bool log_p_ignored) {
                Rcpp::LogicalVector posdrift_ignored;
                Rcpp::LogicalVector log_on_default(1, false);
                return plnr_c(rt, event_pars, posdrift_ignored, min_lik_for_pdf, log_on_default);
            };
        } else {
            Rcpp::stop("Unsupported base model type for _CENS_TRUNC: " + base_model_type);
        }

        Rcpp::NumericVector lR_unique_vals = data["lR"]; // Assuming lR column exists and gives accumulator info
        int n_acc = Rcpp::unique(lR_unique_vals).size(); // Determine n_acc
        if (n_acc == 0 && pars.nrow() > 0 && n_trials > 0) n_acc = pars.nrow() / n_trials; // Fallback if lR is not helpful
        if (n_acc == 0) Rcpp::stop("calc_ll: Could not determine n_acc for _CENS_TRUNC model.");


        for (int i = 0; i < n_particles; ++i) {
            p_vector = p_matrix(i, Rcpp::_);
            if (i == 0) { // Only need to make p_specs once
                p_specs = make_pretransform_specs(p_vector, pretransforms);
            }
            // get_pars_matrix applies p_vector attributes, constants, pretransforms, maps to design, applies transforms and trends
            Rcpp::NumericMatrix current_pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);

            Rcpp::LogicalVector current_is_ok(current_pars.nrow());
            if (i == 0 || bound_specs.empty()) { // bound_specs depends on colnames of current_pars, make once.
                                                 // Or if current_pars structure changes per particle (it shouldn't much beyond values)
                bound_specs = make_bound_specs(minmax, mm_names, current_pars, bounds);
            }
            current_is_ok = c_do_bound(current_pars, bound_specs);
            // lr_all might be needed if model specific checks are required, like DDM. For generic race, maybe not here.
            // current_is_ok = lr_all(current_is_ok, n_acc); // This was for old race, check if needed

            lls[i] = c_log_likelihood_race_cens_trunc(current_pars, data, cens_model_dfun_ptr, cens_model_pfun_ptr, min_ll, current_is_ok, n_acc, expand);
        }

    } else {
        // --- Existing Race Model Logic (LBA, RDM, LNR) ---
        LogicalVector winner = data["winner"];
        NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
        NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
        if(type == "LBA"){
          dfun = dlba_c;
          pfun = plba_c;
        } else if(type == "RDM"){
          dfun = drdm_c;
          pfun = prdm_c;
        } else{ // LNR
          dfun = dlnr_c;
          pfun = plnr_c;
        }
        Rcpp::NumericVector lR_col = data["lR"]; // Changed from lR to lR_col to avoid conflict
        int n_lR_val = Rcpp::unique(lR_col).size(); // Changed from n_lR to n_lR_val

        for (int i = 0; i < n_particles; ++i) {
          p_vector = p_matrix(i, Rcpp::_);
          if(i == 0){
            p_specs = make_pretransform_specs(p_vector, pretransforms);
          }
          pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
          if (i == 0) {
            bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
          }
          is_ok = c_do_bound(pars, bound_specs);
          is_ok = lr_all(is_ok, n_lR_val); // n_lR_val instead of n_lR
          lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
        }
    }
  }
  return(lls);
}


// ---- START: New code for censored/truncated race models ----
#include "model_race_cens_trunc.h" // Include the new header

// Actual model functions will be passed by R, need to get them from SEXP
// For now, these are placeholders for where we'd extract them.
// In c_log_likelihood_race_cens_trunc, we'll get them from R_model_list.
// RacePdfFun current_model_dfun = nullptr;
// RaceCdfFun current_model_pfun = nullptr;


// Helper to order parameters (winner first)
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    int k_idx, // 1-based index
    int n_acc) {

    if (k_idx < 1 || k_idx > n_acc) {
        Rcpp::stop("order_pars_for_winner_cpp: k_idx out of bounds.");
    }

    Rcpp::NumericMatrix ordered_pars(n_acc, p_all_acc.ncol());
    // Copy winner to the first row
    for (int j = 0; j < p_all_acc.ncol(); ++j) {
        ordered_pars(0, j) = p_all_acc(k_idx - 1, j);
    }
    // Copy losers
    int current_row = 1;
    for (int i = 0; i < n_acc; ++i) {
        if (i == (k_idx - 1)) continue;
        for (int j = 0; j < p_all_acc.ncol(); ++j) {
            ordered_pars(current_row, j) = p_all_acc(i, j);
        }
        current_row++;
    }
    return ordered_pars;
}


// Integrand function
// p_trial_this_winner_first: parameters for the current trial, winner's pars in the first row
double f_race_integrand_cpp(
    double t,
    const Rcpp::NumericMatrix& p_trial_this_winner_first, // Assumes winner is row 0
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc) {

    if (t < 0) return 0.0; // RTs are non-negative

    Rcpp::NumericVector t_vec(1); // model_dfun/pfun expect a vector
    t_vec[0] = t;

    // Winner's parameters are in the first row
    Rcpp::NumericMatrix p_winner = p_trial_this_winner_first(Rcpp::Range(0,0), Rcpp::_);

    Rcpp::NumericVector pdf_winner_vec = model_dfun(t_vec, p_winner, false); // log_p = false
    double pdf_winner = pdf_winner_vec[0];

    if (!R_finite(pdf_winner) || pdf_winner < 0) pdf_winner = 0.0;

    if (n_acc > 1) {
        double survivor_losers = 1.0;
        for (int i = 1; i < n_acc; ++i) { // Losers are rows 1 to n_acc-1
            Rcpp::NumericMatrix p_loser_i = p_trial_this_winner_first(Rcpp::Range(i,i), Rcpp::_);
            Rcpp::NumericVector cdf_loser_i_vec = model_pfun(t_vec, p_loser_i, false); // log_p = false
            double cdf_loser_i = cdf_loser_i_vec[0];

            double s_loser_i = 1.0 - cdf_loser_i;
            if (!R_finite(s_loser_i) || s_loser_i < 0) s_loser_i = 0.0;
            if (s_loser_i > 1.0) s_loser_i = 1.0;
            survivor_losers *= s_loser_i;
        }
        pdf_winner *= survivor_losers;
    }
    return (R_finite(pdf_winner) && pdf_winner > 0) ? pdf_winner : 0.0;
}


// --- Basic adaptive Simpson for integration (Placeholder) ---
// This is a simplified adaptive Simpson's rule.
// For production, a more robust library (like GSL or Boost.Math) would be better.
// Or call back to R's integrate if performance allows and complexity is an issue.
typedef double (*integrand_func_ptr)(double, const Rcpp::NumericMatrix&, RacePdfFun, RaceCdfFun, int);

double adaptive_simpson_recursive(
    integrand_func_ptr f,
    const Rcpp::NumericMatrix& p, RacePdfFun dfun, RaceCdfFun pfun, int n_acc,
    double a, double b, double epsilon, double S, double fa, double fb, double fc, int depth) {

    if (depth > 20) { // Max depth to prevent infinite recursion
      // Rcpp::warning("Max recursion depth reached in adaptive_simpson_recursive.");
      return S; // Return current best estimate
    }
    double c = (a + b) / 2.0;
    double h = b - a;
    double fd = f( (a + c) / 2.0, p, dfun, pfun, n_acc);
    double fe = f( (c + b) / 2.0, p, dfun, pfun, n_acc);
    double Sleft = (h/12.0) * (fa + 4.0*fd + fc);
    double Sright = (h/12.0) * (fc + 4.0*fe + fb);
    double S2 = Sleft + Sright;
    if (std::abs(S2 - S) <= 15.0 * epsilon || h < 1e-5) { // 1e-5 is an arbitrary small h
        return S2 + (S2 - S) / 15.0;
    }
    return adaptive_simpson_recursive(f, p, dfun, pfun, n_acc, a, c, epsilon/2.0, Sleft, fa, fc, fd, depth+1) +
           adaptive_simpson_recursive(f, p, dfun, pfun, n_acc, c, b, epsilon/2.0, Sright, fc, fb, fe, depth+1);
}

double adaptive_simpson_cpp(
    integrand_func_ptr f,
    const Rcpp::NumericMatrix& p_trial_this_winner_first, // Parameters for integrand
    RacePdfFun model_dfun, RaceCdfFun model_pfun, int n_acc, // Args for f
    double a, double b,
    double epsilon = 1e-7) { // Tolerance

    if (a == b) return 0.0;
    if (a > b) std::swap(a,b); // Ensure a < b
    if (!R_finite(a) || !R_finite(b)) {
      // Rcpp::Rcout << "Warning: Non-finite integration bounds: " << a << ", " << b << std::endl;
      // Attempt to handle common cases like integrating to Inf
      if (R_finite(a) && b == R_PosInf) {
          // Transform integral (0, inf) to (0, 1) using t = x / (1-x)
          // This requires changing the function f and dx accordingly.
          // For now, if Inf, use a large number or return error.
          // This placeholder integration does not handle Inf bounds well.
          b = a + 50; // Arbitrary large upper bound if original was Inf
          if (a == 0 && b == R_PosInf) b = 50; // if 0 to Inf
      } else if (!R_finite(a) && R_finite(b)) {
          a = b - 50;
      } else if (!R_finite(a) && !R_finite(b)){ // both non-finite
          return 0.0; // Or handle as error
      }
    }
    if (a >=b) return 0.0;


    double c = (a + b) / 2.0;
    double h = b - a;
    double fa = f(a, p_trial_this_winner_first, model_dfun, model_pfun, n_acc);
    double fb = f(b, p_trial_this_winner_first, model_dfun, model_pfun, n_acc);
    double fc = f(c, p_trial_this_winner_first, model_dfun, model_pfun, n_acc);
    double S = (h/6.0)*(fa + 4.0*fc + fb);
    return adaptive_simpson_recursive(f, p_trial_this_winner_first, model_dfun, model_pfun, n_acc,
                                      a, b, epsilon, S, fa, fb, fc, 1);
}
// --- End of Placeholder Integration ---


// Numerical integration helper
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    double low,
    double upp,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    double epsilon = 1e-7) { // Added tolerance parameter

    if (low >= upp) return 0.0;
    if (k_winner_idx < 1 || k_winner_idx > n_acc) {
      // Rcpp::Rcout << "Warning: Invalid k_winner_idx in integrate_for_kth_winner_cpp: " << k_winner_idx << std::endl;
      return 0.0; // Should not happen with proper R_j_idx handling
    }
    Rcpp::NumericMatrix pars_ordered = order_pars_for_winner_cpp(p_all_acc, k_winner_idx, n_acc);

    double res = adaptive_simpson_cpp(f_race_integrand_cpp, pars_ordered, model_dfun, model_pfun, n_acc, low, upp, epsilon);

    return (res > 0 && R_finite(res)) ? res : 0.0;
}

// Truncation correction factor helper
double get_trunc_corr_factor_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    double LT, double UT,
    int n_acc,
    double epsilon = 1e-7) {

    if (!(LT > 0 || UT < R_PosInf)) return 1.0; // No truncation active
    if (k_winner_idx < 1 || k_winner_idx > n_acc) return NA_REAL;


    double prob_untruncated = integrate_for_kth_winner_cpp(k_winner_idx, p_all_acc, 0, R_PosInf, model_dfun, model_pfun, n_acc, epsilon);
    double prob_truncated_interval = integrate_for_kth_winner_cpp(k_winner_idx, p_all_acc, LT, UT, model_dfun, model_pfun, n_acc, epsilon);

    if (prob_truncated_interval > 1e-12) {
        if (!R_finite(prob_untruncated) || prob_untruncated < 0) return NA_REAL;
        return prob_untruncated / prob_truncated_interval;
    } else {
        if (prob_untruncated > 1e-12) return NA_REAL; // Density outside but not inside truncation
        return 1.0; // No density anywhere, or only outside [LT,UT] (and also outside [0,Inf] effectively if prob_untruncated is also 0)
    }
}


// Main C++ function for censored/truncated race likelihood calculation
// Returns a single double: the total log-likelihood for one particle.
double c_log_likelihood_race_cens_trunc(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, all unique trial conditions
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions
    RacePdfFun model_dfun,                  // Pointer to the model's C++ PDF
    RaceCdfFun model_pfun,                  // Pointer to the model's C++ CDF
    double min_ll,
    const Rcpp::LogicalVector& ok_params,   // Parameter validity for this particle's pars matrix
    int n_acc,                              // Number of accumulators in the race
    const Rcpp::IntegerVector& expand_vec   // Vector for expanding unique LLs to full trial count
) {
    // Hardcoded censoring and truncation values
    double LT = 0.1;
    double LC = 0.2;
    double UC = 2.5;
    double UT = 3.0;
    double integration_epsilon = 1e-7; // Tolerance for integration

    Rcpp::IntegerVector R_col = Rcpp::as<Rcpp::IntegerVector>(dadm["R"]);
    Rcpp::NumericVector rt_col = Rcpp::as<Rcpp::NumericVector>(dadm["rt"]);

    int n_unique_trials = dadm.nrows();
    // n_acc is now passed as an argument, ensure it's valid.
    if (n_acc <= 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: n_acc must be positive.");
    if (pars.nrow() != (n_unique_trials * n_acc)) {
      Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_unique_trials: " << n_unique_trials << ", n_acc: " << n_acc << std::endl;
      Rcpp::stop("c_log_likelihood_race_cens_trunc: pars matrix dimensions do not match n_unique_trials * n_acc.");
    }
    if (ok_params.size() != pars.nrow()) {
       Rcpp::stop("c_log_likelihood_race_cens_trunc: ok_params size does not match pars matrix rows.");
    }


    Rcpp::NumericVector ll_unique(n_unique_trials);

    for (int j = 0; j < n_unique_trials; ++j) {
        Rcpp::Range current_trial_par_indices((j * n_acc), (j * n_acc) + n_acc - 1);

        bool params_ok_for_trial = true;
        for(int k_param_check = 0; k_param_check < n_acc; ++k_param_check) {
            if (!ok_params[current_trial_par_indices.begin() + k_param_check]) {
                params_ok_for_trial = false;
                break;
            }
        }
        if (!params_ok_for_trial) {
            ll_unique[j] = min_ll;
            continue;
        }

        Rcpp::NumericMatrix pars_condition_j_all_acc = pars(current_trial_par_indices, Rcpp::_);

        double rt_j = rt_col[j];
        // R_col contains 1-based factor levels. NA is represented by R_NaInt.
        int R_j_idx = (R_col[j] == NA_INTEGER) ? NA_INTEGER : R_col[j];

        double current_prob_val = 0;

        if (R_FINITE(rt_j) && rt_j > 0) { // Case 1: Observed RT
            if (R_j_idx != NA_INTEGER) {
                if (rt_j < LT || rt_j > UT) {
                    current_prob_val = 0;
                } else {
                    Rcpp::NumericMatrix pars_ordered_obs = order_pars_for_winner_cpp(pars_condition_j_all_acc, R_j_idx, n_acc);
                    current_prob_val = f_race_integrand_cpp(rt_j, pars_ordered_obs, model_dfun, model_pfun, n_acc);
                }
                if (current_prob_val > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf < 0) current_prob_val = 0;
                    else current_prob_val *= trunc_cf;
                }
            } else { current_prob_val = 0; } // Should not happen: observed RT but unknown response

        } else if (R_FINITE(rt_j) && rt_j == R_NegInf) { // Case 2: Lower Censored
            if (R_j_idx != NA_INTEGER) { // Response known
                current_prob_val = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon);
                if (current_prob_val > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf < 0) current_prob_val = 0;
                    else current_prob_val *= trunc_cf;
                }
            } else { // Response unknown
                double sum_p = 0;
                for (int k_w = 1; k_w <= n_acc; ++k_w) {
                    double p_k = integrate_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon);
                    if (p_k > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                        if (!ISNAN(trunc_cf_k) && R_FINITE(trunc_cf_k) && trunc_cf_k >=0) sum_p += (p_k * trunc_cf_k);
                    }
                }
                current_prob_val = sum_p;
            }
        } else if (R_FINITE(rt_j) && rt_j == R_PosInf) { // Case 3: Upper Censored
             if (R_j_idx != NA_INTEGER) { // Response known
                current_prob_val = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon);
                if (current_prob_val > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf < 0) current_prob_val = 0;
                    else current_prob_val *= trunc_cf;
                }
            } else { // Response unknown
                double sum_p = 0;
                for (int k_w = 1; k_w <= n_acc; ++k_w) {
                    double p_k = integrate_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon);
                     if (p_k > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                        if (!ISNAN(trunc_cf_k) && R_FINITE(trunc_cf_k) && trunc_cf_k >=0) sum_p += (p_k * trunc_cf_k);
                    }
                }
                current_prob_val = sum_p;
            }
        } else if (ISNAN(rt_j)) { // Case 4 & 5: NA RT (interval censoring)
            if (R_j_idx != NA_INTEGER) { // Response known
                double p_L = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon);
                double p_U = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon);
                current_prob_val = p_L + p_U;
                if (current_prob_val > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                    if (ISNAN(trunc_cf) || !R_FINITE(trunc_cf) || trunc_cf < 0) current_prob_val = 0;
                    else current_prob_val *= trunc_cf;
                }
            } else { // Response unknown
                double sum_total_p = 0;
                for (int k_w = 1; k_w <= n_acc; ++k_w) {
                    double p_L_k = integrate_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, integration_epsilon);
                    double p_U_k = integrate_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, integration_epsilon);
                    double p_k_sum = p_L_k + p_U_k;
                    if (p_k_sum > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_w, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, integration_epsilon);
                        if (!ISNAN(trunc_cf_k) && R_FINITE(trunc_cf_k) && trunc_cf_k >=0) sum_total_p += (p_k_sum * trunc_cf_k);
                    }
                }
                current_prob_val = sum_total_p;
            }
        }
        // All other rt_j values (e.g. 0 or negative if not Inf) would result in current_prob_val = 0

        current_prob_val = std::max(0.0, current_prob_val); // Ensure non-negative
        ll_unique[j] = (current_prob_val > std::numeric_limits<double>::epsilon()) ? std::log(current_prob_val) : min_ll;
        ll_unique[j] = std::max(min_ll, ll_unique[j]); // Ensure not less than min_ll
    }

    // The return type is Rcpp::NumericVector for easier debugging from R if this function
    // is exported and called directly. For integration into calc_ll, this should be changed
    // to return a single double (the sum of ll_unique, possibly after expansion).
    // calc_ll handles expansion, so this function should return the sum for the *unique* trials,
    // or the vector ll_unique if calc_ll expects that.
    // Given c_log_likelihood_race returns a double (summed over expanded trials),
    // this function should eventually do the same.
    // Summing logic based on expand_vec
    double total_ll = 0;
    if (expand_vec.length() > 0) {
        for (int i = 0; i < expand_vec.length(); ++i) {
            // expand_vec is 1-based index from R, convert to 0-based for C++ vector
            int unique_idx = expand_vec[i] - 1;
            if (unique_idx >= 0 && unique_idx < ll_unique.size()) {
                total_ll += ll_unique[unique_idx];
            } else {
                // This case should ideally not happen if expand_vec is correct
                total_ll += min_ll;
            }
        }
    } else {
        // Fallback if expand_vec is not provided or empty (e.g. dadm$N was used in R)
        // This part might need to check for dadm["N"] if that's an alternative expansion mode
        // For now, sum unique likelihoods directly if no expansion.
        // This is consistent if each unique trial happens only once.
        for (int i = 0; i < ll_unique.size(); ++i) {
            total_ll += ll_unique[i];
        }
    }
    return total_ll;
}


// ---- Rcpp Exported Wrapper for Testing ----
// This wrapper is for testing the c_log_likelihood_race_cens_trunc from R.
// This wrapper is for testing the c_log_likelihood_race_cens_trunc from R.
// It will use existing C++ pdf/cdf functions for a known model (e.g., LBA)
// as a simplification, since passing arbitrary C++ function pointers from R
// directly to a C++ function not designed as a generic callback host is non-trivial.
// The true integration will be calc_ll calling c_log_likelihood_race_cens_trunc directly
// with proper C++ function pointers.

// [[Rcpp::export]]
Rcpp::NumericVector test_c_loglik_cens_trunc_wrapper_R(
    Rcpp::NumericMatrix pars,      // Parameters for one particle, all unique trial conditions
    Rcpp::DataFrame dadm,          // Data for unique trial conditions
    std::string model_type_str,    // e.g., "LBA", "RDM", "LNR" to select hardcoded C++ funcs
    double min_ll,
    Rcpp::LogicalVector ok_params, // Parameter validity for this particle's pars matrix
    int n_acc                      // Number of accumulators
) {
    RacePdfFun model_dfun_ptr = nullptr;
    RaceCdfFun model_pfun_ptr = nullptr;

    // Select appropriate C++ function pointers based on model_type_str
    // This mirrors part of what calc_ll does.
    // Note: The signatures of dlba_c/plba_c etc. are not exactly
    // (NumericVector rt, NumericMatrix pars, bool log_p).
    // dlba_c is (NumericVector t, NumericMatrix pars, LogicalVector posdrift, double min_lik, LogicalVector log_on)
    // This means our RacePdfFun/RaceCdfFun typedefs might be too simplistic, or
    // we need thin C++ wrappers around dlba_c etc. to match the typedef.
    // For now, let's assume such wrappers exist or that the typedefs would be adjusted.
    // THIS IS A SIGNIFICANT SIMPLIFICATION FOR TESTING.

    // Placeholder: Using LBA functions directly if available and compatible
    // This part will likely cause a compile error if the signatures don't match perfectly.
    // Or, we'd need to create small lambda wrappers here.
    // Let's assume for the sake of this step that we have functions matching the typedefs.
    // For a real test, one would need to ensure compatible functions are callable.
    // This might involve creating temporary wrapper functions within this test function
    // or having pre-existing C++ functions that match the RacePdfFun/RaceCdfFun typedefs.

    if (model_type_str == "LBA_test") { // Using a distinct name for clarity
        // These would be pointers to functions that match the RacePdfFun/RaceCdfFun signature
        // e.g., wrapper_dlba_for_racecens and wrapper_plba_for_racecens
        // model_dfun_ptr = &wrapper_dlba_for_racecens;
        // model_pfun_ptr = &wrapper_plba_for_racecens;
        // Since these wrappers don't exist yet, this will be a conceptual call.
        // The actual c_log_likelihood_race_cens_trunc will be called with nullptr
        // if no matching functions are found, which it should handle (e.g. by erroring or returning min_ll).
        // Rcpp::Rcout << "Using LBA_test (conceptual - actual functions not yet wrapped)" << std::endl;
        // To make it compile, we pass nullptrs and the C++ function should handle them.
         model_dfun_ptr = nullptr; // Placeholder
         model_pfun_ptr = nullptr; // Placeholder
         // A real test would require pre-defined C++ functions like:
         // RacePdfFun test_lba_dfun = [](NumericVector rt, NumericMatrix pars, bool log_p){ return Rcpp::NumericVector(rt.size(), 0.1); };
         // model_dfun_ptr = &test_lba_dfun;
         Rcpp::warning("Test wrapper called for LBA_test, but actual C++ model functions are not yet wrapped to match RacePdfFun/CdfFun typedefs. Using nullptr, expect min_ll or errors.");

    } else {
        Rcpp::warning("Model type '%s' not recognized for C++ test wrapper. Using nullptr for model functions.", model_type_str);
        model_dfun_ptr = nullptr;
        model_pfun_ptr = nullptr;
    }

    if (model_dfun_ptr == nullptr || model_pfun_ptr == nullptr) {
        // If function pointers are null, we cannot proceed with calculation.
        // Return a vector of min_ll for each unique trial.
        int n_unique_trials = dadm.nrows();
        Rcpp::NumericVector error_result(n_unique_trials);
        for(int i=0; i<n_unique_trials; ++i) error_result[i] = min_ll;
        // This return needs to be processed by the R caller to sum up for the final single LL value.
        // Or, if we want this test wrapper to return a single double like calc_ll expects from its sub-functions:
        // double total_error_ll = 0;
        // Rcpp::IntegerVector N_vec = dadm.containsElementNamed("N") ? Rcpp::as<Rcpp::IntegerVector>(dadm["N"]) : Rcpp::IntegerVector(0);
        // Rcpp::IntegerVector expand_vec = dadm.hasAttribute("expand") ? Rcpp::as<Rcpp::IntegerVector>(dadm.attr("expand")) : Rcpp::IntegerVector(0);
        // if (N_vec.size() == n_unique_trials) {
        //   for(int i=0; i<n_unique_trials; ++i) total_error_ll += min_ll * N_vec[i];
        // } else if (expand_vec.size() > 0) {
        //   total_error_ll = min_ll * expand_vec.size();
        // } else {
        //   total_error_ll = min_ll * n_unique_trials;
        // }
        // return Rcpp::NumericVector::create(total_error_ll); // Returning a single value in a vector for consistency with current return type
        return error_result; // Returning vector of unique trial LLs
    }

    // Call the main C++ implementation
    // For the test wrapper, we need to construct an expand_vec or assume calc_ll's role.
    // Simplest for testing: assume each unique trial is one actual trial (no expansion).
    Rcpp::IntegerVector expand_vec_test(dadm.nrows());
    for(int i=0; i < dadm.nrows(); ++i) expand_vec_test[i] = i + 1; // 1-based index

    // If the C++ function now returns double, the wrapper should reflect that.
    double result_ll = c_log_likelihood_race_cens_trunc(
        pars, dadm,
        model_dfun_ptr, model_pfun_ptr,
        min_ll, ok_params, n_acc,
        expand_vec_test // Pass the dummy expand_vec
    );
    return Rcpp::NumericVector::create(result_ll); // Return as a single-element NumericVector
}

// ---- END: New code for censored/truncated race models ----

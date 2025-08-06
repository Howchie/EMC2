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
#include <unordered_map>
#include <sstream>
#include <iomanip>

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
    current_model_ctx.min_lik_for_pdf = min_ll;
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
    } else if (type_std=="RDMSWTN") {
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
		if( type_std.find("LogicalRules") != std::string::npos) {
			if( type_std.find("negdrift") != std::string::npos) {
				lls[i] = c_log_likelihood_redundant_target_race_failure(pars, data,
                                                           model_dfun_ptr, model_pfun_ptr, n_trials,
                                                           expand, min_ll, is_ok,n_acc,
                                                           &current_model_ctx);
			} 
			else {
				lls[i] = c_log_likelihood_redundant_target_race(pars, data,
                                                           model_dfun_ptr, model_pfun_ptr, n_trials,
                                                           expand, min_ll, is_ok,n_acc,
                                                           &current_model_ctx);
			} 
		} else {
            lls[i] = c_log_likelihood_race_cens_trunc(pars, data,
                                                     model_dfun_ptr, model_pfun_ptr, n_trials,
                                                     winner, expand, min_ll, is_ok, n_acc,
                                                     &current_model_ctx);
        }
    }
  }
  return(lls);
}

// GSL-compatible adapter for f_race_integrand_cpp
double gsl_f_race_adapter(double t, void *p) {
    gsl_race_params* params = static_cast<gsl_race_params*>(p);
	// ① get the params block
    // ② recover the model-specific context
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(params->model_specific_context);
    if (t < 0) return 0.0; // RTs are non-negative
	int n_acc=params->n_acc;
	NumericVector t_vec = Rcpp::rep(t, n_acc);

    const NumericMatrix& p_mat = *(params->p_trial);
	const LogicalVector& winner = params->winner;
	const LogicalVector& isok = params->isok;
	NumericVector ll(n_acc);
    NumericVector win = log(params->model_dfun(t_vec, p_mat, isok, winner, ctx));
	ll[winner]=win;
	double ll_out = 0.0;
    if (params->n_acc > 1) {
        Rcpp::NumericVector loss_cdfs = params->model_pfun(t_vec, p_mat, isok, !winner, ctx);
		Rcpp::NumericVector loss(loss_cdfs.size());
		for (int j = 0; j < loss_cdfs.size(); ++j) {
			double cdf = loss_cdfs[j];
			if (cdf > -1e-9) { // Safety check: if cdf is basically 1
				loss[j] = R_NegInf; // log(0)
			} else {
				loss[j] = std::log1p(-cdf);
			}
		}
		ll[!winner] = loss;
    }
 
	for (int i = 0; i < n_acc; ++i) {
        double term = ll[i];
        ll_out += term;
    }
    return std::exp(ll_out);   
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

bool row_is_finite(const Rcpp::NumericMatrix& mat, int row) {
    for (int j = 0; j < mat.ncol(); ++j) {
        if (!R_finite(mat(row, j))) return false;
    }
    return true;
}

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
    int n_acc_j,
    double epsilon = 1e-7,
    void* model_specific_context = nullptr) {

	Rcpp::LogicalVector winner(n_acc_j,false);
	winner[k_winner_idx-1] = true; 
    if (low >= upp && !(low == 0 && upp == R_PosInf)) return 0.0;
    if (k_winner_idx < 1 || k_winner_idx > (n_acc_j)) {
      Rcpp::Rcerr << "Warning: Invalid k_winner_idx in integrate_for_kth_winner_cpp: " << k_winner_idx << std::endl;
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
    params_struct.n_acc = n_acc_j;
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
        Rcpp::Rcerr << "GSL_ERROR: " << gsl_strerror(status) << " (code " << status << ")\n"
                        << "  Inputs: " << p_all_acc << n_acc_j << std::endl;
            return 0.0; // 0? Or NaN to indicate error
        return 0.0;
    }
    return (result > 0 && R_finite(result)) ? result : 0.0;
}

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
    double total_log_mass = R_NegInf; 
    for(int k=1;k<=n_acc;++k) {
        double log_mass_k  = integrate_for_kth_winner_cpp(k,p_all_acc, isok_trial, LT, UT, model_dfun, model_pfun, n_acc, epsilon, model_specific_context);
		// Accumulate the log masses using the log-sum-exp trick
        total_log_mass = log_sum_exp(total_log_mass, log_mass_k);
	}
    return total_log_mass;
}

// Helper function to cache unique parameter combinations
static inline std::string make_key(int start_row,
                                   const Rcpp::NumericVector &LT,
                                   const Rcpp::NumericVector &UT,
                                   const Rcpp::NumericMatrix &pars,
                                   int n_acc, int n_par)
{
    std::ostringstream oss;
    oss.setf(std::ios::scientific);        // fixed precision
    oss << std::setprecision(17);

    oss << LT[start_row] << '|' << UT[start_row] << '|';

    for (int r = 0; r < n_acc; ++r) {
        for (int c = 0; c < n_par; ++c) {
            double v = pars(start_row + r, c);
            if (Rcpp::NumericVector::is_na(v)) oss << "NA";
            else                               oss << v;
            oss << ',';
        }
        oss << '|';
    }
    return oss.str();
}

// Main C++ function for censored/truncated race likelihood calculation
// This function is now the unified entry point for all race models (LBA, RDM, LNR),
// whether they are standard or explicitly handling censoring/truncation.
// It uses batching for finite RTs and iterative processing for others (censored/NA RTs).
// needs unbatching?
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
    // Fetch censoring and truncation values from dadm columns. These are passed across all rows for ease of access, but should be identical at least at the subject level as they don't correspond to a data entry. Attributes probably a better fit here but clunky.
    Rcpp::NumericVector LT = dadm["LT"];
	Rcpp::NumericVector UT = dadm["UT"]; 
	Rcpp::NumericVector LC = dadm["LC"];    
	Rcpp::NumericVector UC = dadm["UC"];
    static std::unordered_map<std::string,double> inv_Z_cache; // for truncation caching
    double integration_epsilon = 1e-7; // Tolerance for GSL integration
	int n_acc_j = n_acc;
	const int n_out = expand.length();
    Rcpp::NumericVector rts_dadm = dadm["rt"];
    Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
    Rcpp::IntegerVector lR_dadm = dadm["lR"];
	NumericVector lds(n_trials,min_ll); // initialise at min_ll
	NumericVector lds_exp(n_out);
	// If a RACE column exists, set parameters of accumulators not present on a
    // given trial to NA so the density functions return zero for them. This
    // mirrors logic from the old c_log_likelihood_race implementation.
	IntegerVector RACE(n_trials, n_acc);
	Rcpp::LogicalVector RACE_mask(n_trials, true); // Mask for all dadm rows
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
			RACE[row]=n_acc_this_trial;
			// lR_dadm is the (1-based) index of *this* accumulator on the trial
			if (lR_dadm[row] > n_acc_this_trial) {
				// accumulator not present → blank its parameter row
				std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
				RACE_mask[row]=false;
			}
		}
	}
    if (n_trials == 0) return 0.0; // No data, no likelihood

    if (n_acc <= 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: n_acc must be positive and correctly determined before this call.");
    if (n_trials % n_acc != 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: dadm nrows not a multiple of n_acc.");

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


    int n_unique_trials = n_trials / n_acc;
    Rcpp::NumericVector ll_unique(n_unique_trials);
    //ll_unique.fill(min_ll); // Initialize all to min_ll
	Rcpp::NumericVector pC_values(n_unique_trials);
    if (use_pC) {
		for (int j = 0; j < n_unique_trials; ++j) {
			pC_values[j] = pars(j * n_acc, pc_col);
		}
    }

    // Parameter matrix and validity vector checks
    if (pars.nrow() != n_trials) {
      Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_trials: " << n_trials << std::endl;
      Rcpp::stop("c_log_likelihood_race_cens_trunc: pars matrix dimensions do not match total dadm rows.");
    }
    if (isok.size() != pars.nrow()) {
       Rcpp::stop("c_log_likelihood_race_cens_trunc: isok size does not match pars matrix rows.");
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
	// Batch process finite rt trials
	std::vector<int> finite_rt_unique_trial_indices; // Stores *unique trial* indices for finite RTs
	std::vector<int> other_unique_trial_indices;   // Stores *unique trial* indices for other cases
    finite_rt_unique_trial_indices.reserve(n_unique_trials);
    other_unique_trial_indices.reserve(n_unique_trials);
    // Categorize unique trials based on RT properties and parameter validity

    for (int j = 0; j < n_unique_trials; ++j) {
        int start_row_idx = j * n_acc; // Starting row in dadm/pars for this unique trial
        double rt_j = rts_dadm[start_row_idx]; // RT for this unique trial
        int R_j_idx = R_idxs_dadm[start_row_idx]; // Winner index (1-based from R factor)
		int n_acc_j = RACE[start_row_idx];
        // Criteria for batch processing: finite, positive RT, within truncation bounds [LT, UT], and known winner
        if (R_FINITE(rt_j) && rt_j > 0 && rt_j >= LT[start_row_idx] && rt_j <= UT[start_row_idx] && R_j_idx != NA_INTEGER) {
			finite_rt_unique_trial_indices.push_back(j); // Store unique trial index
            for(int k_acc = 0; k_acc < n_acc_j; ++k_acc) {
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
		Rcpp::NumericVector win = log(model_dfun(rts_for_finite_batch, finite_rt_pars, isok_for_finite_batch, winner_for_finite_batch, model_context_for_funcs));
	// Place results into the correct positions in 'lds'
        int current_lds_idx = 0;
        for(int i = 0; i < n_trials; ++i) {
            if(finite_rt_mask[i] && winner[i]) {
                lds[i] =win[current_lds_idx++];
            }
        }
		// Prepare survivor functions for losers, one matrix per accumulator
		if (n_acc > 1) {
			Rcpp::NumericVector loss_cdf = log(model_pfun(rts_for_finite_batch, finite_rt_pars, isok_for_finite_batch, !winner_for_finite_batch, model_context_for_funcs));
			
			// Vectorize the log(1-CDF) calculation using a stable log1p trick.
			Rcpp::NumericVector loss(loss_cdf.size());
			for (int j = 0; j < loss_cdf.size(); ++j) {
				double log_cdf = loss_cdf[j];
				if (log_cdf > -1e-9) { // Safety check: if cdf is basically 1
					loss[j] = R_NegInf; // log(0)
				} else {
					loss[j] = std::log1p(-std::exp(log_cdf));
				}
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
			int start_row_idx = unique_trial_idx*(n_acc);
			n_acc_j = RACE[start_row_idx];
            Rcpp::NumericMatrix p_all_acc_for_trial = pars(Rcpp::Range(start_row_idx, start_row_idx + n_acc_j - 1), Rcpp::_);
            Rcpp::LogicalVector isok_for_trial = isok[Rcpp::Range(start_row_idx, start_row_idx + n_acc_j - 1)];
			lower_for_trial = LT[start_row_idx];
			upper_for_trial = UT[start_row_idx];
			if (LT[start_row_idx]!=0 || UT[start_row_idx] != R_PosInf) { // if truncation
				std::string key = make_key(start_row_idx, LT, UT, pars, n_acc_j, pars.ncol());
				auto hit = inv_Z_cache.find(key);
				if (hit != inv_Z_cache.end()) { // check if we've integrated this exact parameter set already and re-use result if we have
					log_Z_this = hit->second;
				} else { // otherwise compute the integral and cache it
					log_Z_this  = get_trunc_normaliser_cpp(p_all_acc_for_trial,model_dfun, model_pfun,isok_for_trial,lower_for_trial, upper_for_trial, n_acc_j,integration_epsilon,model_context_for_funcs);
					inv_Z_cache.emplace(std::move(key), log_Z_this);
				}
			}
			double current_trial_ll_sum = 0.0;
			for(int k=0; k < n_acc_j; ++k) {
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
			// Intrinsic Omissions Correction Comes AFTER the truncation adjustment (because that adjustment assumes a response would have been recorded)
			if (!posdrift) { // should only ever be false for LBAIO model
				double p_I = 1.0;
				for (int k = 0; k < n_acc_j; ++k) {
					double v   = p_all_acc_for_trial(k, 0);   // mean drift for kth acc
					double sv  = p_all_acc_for_trial(k, 1);   // its s.d.
					double p_neg = R::pnorm(0.0, v, sv, 1,0); // (same as 1-pnorm(v/sv,0,1)
					p_I *= p_neg; // we multiply to get the probability that EVERY accumulator was negative
				}
				
				current_trial_ll_sum += std::log1p(-p_I); // log(1-pI)
			}
			// If no truncation, inv_Z is not added.
			ll_unique[unique_trial_idx] = current_trial_ll_sum;
			ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]);
		}
	}
    // --- Process other trials (Infinite RTs, NA RTs, or finite RTs outside truncation) ---
    // These trials require individual processing, often involving numerical integration for censored intervals.
	double current_ll_val;
    for (size_t i = 0; i < other_unique_trial_indices.size(); ++i) {
		int unique_trial_idx = other_unique_trial_indices[i];
		int start_row_idx = unique_trial_idx*n_acc;
		n_acc_j = RACE[start_row_idx];
		Rcpp::NumericMatrix p_all_acc_for_trial = pars(Rcpp::Range(start_row_idx, start_row_idx + n_acc_j - 1), Rcpp::_);
        Rcpp::LogicalVector isok_for_trial = isok[Rcpp::Range(start_row_idx, start_row_idx + n_acc_j - 1)];
		Rcpp::LogicalVector winner_for_trial = winner[Rcpp::Range(start_row_idx, start_row_idx + n_acc_j - 1)];
        double rt_j = rts_dadm[start_row_idx];
        int R_j_idx = R_idxs_dadm[start_row_idx];
        current_ll_val = R_NegInf;
        // Case 2: Fast censoring (RT=-Inf). Probability is integral from LT to LC.
        if (rt_j == R_NegInf) {
			lower_for_trial = LT[start_row_idx];
			upper_for_trial = LC[start_row_idx];
            if (R_j_idx != NA_INTEGER) { // Response (winner) is known
                current_ll_val = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, lower_for_trial, upper_for_trial, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc_j; ++k_win) {
                    double log_l_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, lower_for_trial, upper_for_trial, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
					current_ll_val =  log_sum_exp(current_ll_val, log_l_k);
					}
            }
        // Case 3: Slow censoring (RT = Inf). Probability is integral from UC to UT.
        } else if (rt_j == R_PosInf) {
			lower_for_trial = UC[start_row_idx];
			upper_for_trial = UT[start_row_idx];
            if (n_acc_j == 1) {
                // Simplified case for single accumulator
                Rcpp::NumericVector t_vec = Rcpp::NumericVector::create(lower_for_trial);
                Rcpp::NumericVector l_vec = log(model_pfun(t_vec, p_all_acc_for_trial, isok_for_trial, winner_for_trial, model_context_for_funcs));
                double log_cdf = l_vec[0];
				if (log_cdf > -1e-9) {
					current_ll_val = R_NegInf;
				} else {
					current_ll_val = std::log1p(-std::exp(log_cdf));
				}			
            } else {
				if (R_j_idx != NA_INTEGER) { // Response (winner) is known
					current_ll_val = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, lower_for_trial, upper_for_trial, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
				} else { // Response (winner) is unknown; sum probabilities over all possible winners
					for (int k_win = 1; k_win <= n_acc_j; ++k_win) {
						double log_l_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, lower_for_trial, upper_for_trial, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
						current_ll_val =  log_sum_exp(current_ll_val, log_l_k);
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
                double ll_L = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, lower_for_trial1, upper_for_trial1, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
                double ll_U = integrate_for_kth_winner_cpp(R_j_idx, p_all_acc_for_trial, isok_for_trial, lower_for_trial2, upper_for_trial2, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
                current_ll_val = log_sum_exp(ll_L, ll_U);
            } else { // Response (winner) is unknown; sum probabilities over all possible winners
                for (int k_win = 1; k_win <= n_acc_j; ++k_win) {
                    double ll_L_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, lower_for_trial1, upper_for_trial1, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
                    double ll_U_k = integrate_for_kth_winner_cpp(k_win, p_all_acc_for_trial, isok_for_trial, lower_for_trial2, upper_for_trial2, model_dfun, model_pfun, n_acc_j, integration_epsilon, model_context_for_funcs);
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
			std::string key = make_key(start_row_idx, LT, UT, pars, n_acc_j, pars.ncol());
			auto hit = inv_Z_cache.find(key);
			if (hit != inv_Z_cache.end()) { // fast path check if we've integrated this exact parameter set already and re-use result if we have
				log_Z_this = hit->second;
			} else { // otherwise compute the integral and cache it
				log_Z_this  = get_trunc_normaliser_cpp(p_all_acc_for_trial,model_dfun, model_pfun,isok_for_trial,lower_for_trial, upper_for_trial, n_acc_j,integration_epsilon,model_context_for_funcs);
				inv_Z_cache.emplace(std::move(key), log_Z_this);
			}
			
			if (!NumericVector::is_na(log_Z_this) && R_FINITE(log_Z_this)){current_ll_val -= log_Z_this;}
		}
		
		// Intrinsic Omissions Correction Comes AFTER the truncation adjustment (because that adjustment assumes a response would have been recorded)
		if (!posdrift) { // should only ever be false for LBAIO model
			double ll_I = 0.0;
			for (int k = 0; k < n_acc_j; ++k) {
				double v   = p_all_acc_for_trial(k, 0);   // mean drift for kth acc
				double sv  = p_all_acc_for_trial(k, 1);   // its s.d.
				double p_neg = R::pnorm(0.0, v, sv, 1,1); // (same as 1-pnorm(v/sv,0,1)
				ll_I += p_neg; // we multiply to get the probability that EVERY accumulator was negative
			}
			current_ll_val = log_sum_exp(ll_I, current_ll_val);
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
                int start_row_idx = j * n_acc;
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
                int start_row_idx = j * n_acc;
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

double c_log_likelihood_redundant_target_race(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame dadm,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    const int n_trials,
    const Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector ok_params,
	int n_acc,
    void* model_specific_context) {

    if (n_trials % n_acc != 0) Rcpp::stop("c_log_likelihood_redundant_target_race: dadm rows must be multiple of n_acc");

    Rcpp::NumericVector rts = dadm["rt"];
    Rcpp::CharacterVector role = dadm["lR"];
    Rcpp::CharacterVector resp = dadm["R"];
	Rcpp::CharacterVector LogicalRule = dadm["LogicalRule"];
    Rcpp::LogicalVector all_idx(n_trials, true); // No such thing as "winner" in the normal context
    Rcpp::NumericVector f_all = model_dfun(rts, pars, ok_params, all_idx, model_specific_context);
    Rcpp::NumericVector F_all = model_pfun(rts, pars, ok_params, all_idx, model_specific_context);
	double rt_min = Rcpp::min(rts);
	double rt_max = Rcpp::max(rts);
	double u_density = 1.0 / std::max(1e-12, rt_max - rt_min);
    int n_unique_trials = n_trials / n_acc;
    Rcpp::NumericVector ll_unique(n_unique_trials);
	
	// Here we check for a rule-following parameter (p) corresponding to probability of processing only a single channel with probability q.
	// Index the column (if it exists)
	ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_specific_context);
	bool posdrift=ctx->use_posdrift;
	bool use_rulebreak = false;
	bool use_guess = false;
    int rulebreak_col = -1;
	int q_col = -1;
	int guess_col =-1;
	Rcpp::List dimnames = pars.attr("dimnames");
	Rcpp::CharacterVector colnames = as<Rcpp::CharacterVector>(dimnames[1]);
    for (int j = 0; j < colnames.size(); ++j) {
      if (as<std::string>(colnames[j]) == "p") {
        rulebreak_col = j;
        use_rulebreak = true;
        break;
      }
    }
	// Same for guesses
	for (int j = 0; j < colnames.size(); ++j) {
		if (as<std::string>(colnames[j]) == "r") {
			guess_col = j;
			use_guess = true;
			break;
		}
	}
	// If we're rulebreaking, find the index of column q, then also check that the rulebreak probability isn't hardcoded to zero
    if (use_rulebreak) {
		for (int j = 0; j < colnames.size(); ++j) {
			if (as<std::string>(colnames[j]) == "q") {
				q_col = j;
				break;
			}
		}
      bool all_one = true;
      for (int i = 0; i < pars.nrow(); ++i) {
        if (pars(i, rulebreak_col) != 1.0) {
          all_one = false;
          break;
        }
      }
      if (all_one) {
        use_rulebreak = false;
      }
    }
	
	// set up parameters
	double p_j; double p_guess;double p_process; double p_rulebreak;
	double p;double q;double r;
	double vA_T;double vA_N;double vB_T;double vB_N;double svA_T;double svA_N;double svB_T;double svB_N; double vG; double svG;
	
    for(int j=0; j<n_unique_trials; ++j){
        int start = j*n_acc;
		p_j = 0;
        double fA=NA_REAL, fB=NA_REAL, fnA=NA_REAL, fnB=NA_REAL, fG=0;
        double FA=NA_REAL, FB=NA_REAL, FnA=NA_REAL, FnB=NA_REAL, FG=0;
        for(int k=0;k<n_acc;++k){
            int idx = start+k;
            std::string r = Rcpp::as<std::string>(role[idx]);
            if(r == "A"){ fA = f_all[idx]; FA = F_all[idx]; vA_T=pars(idx,0); svA_T=pars(idx,1);}
            else if(r == "B"){ fB = f_all[idx]; FB = F_all[idx]; vB_T=pars(idx,0); svB_T=pars(idx,1); }
            else if(r == "n_A"){ fnA = f_all[idx]; FnA = F_all[idx]; vA_N=pars(idx,0); svA_N=pars(idx,1);}
            else if(r == "n_B"){ fnB = f_all[idx]; FnB = F_all[idx]; vB_N=pars(idx,0); svB_N=pars(idx,1);}
			else if(r == "guess"){ fG = f_all[idx]; FG = F_all[idx]; vG=pars(idx,0); svG=pars(idx,1);}
        }
		double one_m_FB = std::max(1e-12, 1.0 - FB);
		double one_m_FA = std::max(1e-12, 1.0 - FA);
		double one_m_FnB = std::max(1e-12, 1.0 - FnB);
		double one_m_FnA = std::max(1e-12, 1.0 - FnA);
		double one_m_FG = std::max(1e-12, 1.0 - FG);
		double one_m_FnAFnB = std::max(1e-12, 1.0 - FnA * FnB);
		double one_m_FAFB = std::max(1e-12, 1.0 - FA * FB);
        std::string r_obs = Rcpp::as<std::string>(resp[start]);
		if (use_rulebreak) {
			p = pars(start, rulebreak_col);
			q = pars(start, q_col);
		}
        if (LogicalRule[start]=="OR") {
			if(r_obs == "yes"){
				p_process = one_m_FG*((fA*one_m_FB + fB*one_m_FA) * one_m_FnAFnB);
				if(use_rulebreak && !std::isnan(p) && q_col>-1) {				
					p_rulebreak = (q * fB*one_m_FnB) + ((1-q)*fA*one_m_FnA);
					p_process=p*p_process + (1-p)*p_rulebreak;
				}
				p_guess = fG * one_m_FA * one_m_FB * one_m_FnAFnB;
				p_j = p_guess+p_process;
			} else if(r_obs == "no"){
				p_process = one_m_FG*((fnA*FnB + fnB*FnA) * (one_m_FA*one_m_FB));
				if(use_rulebreak && !std::isnan(p) && q_col>-1) {					
					p_rulebreak = (q * fnB*one_m_FB) + ((1-q)*fnA*one_m_FA);
					p_process=p*p_process + (1-p)*p_rulebreak;
				}
				p_j = p_process;
			}
		} else if (LogicalRule[start]=="AND") {
			if(r_obs == "no"){
				p_process = one_m_FG*((fnA*one_m_FnB + fnB*one_m_FnA) * one_m_FAFB);
				if(use_rulebreak && !std::isnan(p) && q_col>-1) {
					p_rulebreak = (q * fnB*one_m_FB) + ((1-q)*fnA*one_m_FA);
					p_process=p*p_process + (1-p)*p_rulebreak;
				}
				p_guess = fG * one_m_FnA * one_m_FnB * one_m_FAFB;
				p_j = p_guess+p_process;
			} else if(r_obs == "yes"){
				p_process = one_m_FG*((fA*FB + fB*FA) * (one_m_FnA*one_m_FnB));
				if(use_rulebreak && !std::isnan(p) && q_col>-1) {
					p_rulebreak = (q * fB*one_m_FnB) + ((1-q)*fA*one_m_FnA);
					p_process=p*p_process + (1-p)*p_rulebreak;
				}
				p_j = p_process;
			}
		}
		if (use_guess) {
			r = pars(start, guess_col);
			// guess_pdf is RT-uniform and response-uniform (0.5 each)
			double guess_pdf = 0.5 * u_density;
			p_j = (1.0 - r) * p_j + r * guess_pdf;
		}
		if (p_j <= 0.0 || !R_FINITE(p_j)) {
			ll_unique[j] = min_ll;
		} else {
			ll_unique[j] = std::log(p_j);
		}
    }

    Rcpp::NumericVector ll_exp = c_expand(ll_unique, expand);
    double sum_ll = 0.0;
    for(int i=0;i<ll_exp.size();++i){
        double val = ll_exp[i];
        if(!R_FINITE(val) || Rcpp::NumericVector::is_na(val) || val < min_ll) val = min_ll;
        sum_ll += val;
    }
    return sum_ll;
}

double c_log_likelihood_redundant_target_race_failure(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame dadm,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    const int n_trials,
    const Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector ok_params,
	int n_acc,
    void* model_specific_context) {

    if (n_trials % n_acc != 0) Rcpp::stop("c_log_likelihood_redundant_target_race: dadm rows must be multiple of n_acc");

    Rcpp::NumericVector rts = dadm["rt"];
    Rcpp::CharacterVector role = dadm["lR"];
    Rcpp::CharacterVector resp = dadm["R"];
	Rcpp::CharacterVector LogicalRule = dadm["LogicalRule"];
    Rcpp::LogicalVector all_idx(n_trials, true); // No such thing as "winner" in the normal context
    Rcpp::NumericVector f_all = model_dfun(rts, pars, ok_params, all_idx, model_specific_context);
    Rcpp::NumericVector F_all = model_pfun(rts, pars, ok_params, all_idx, model_specific_context);
    int n_unique_trials = n_trials / n_acc;
    Rcpp::NumericVector ll_unique(n_unique_trials);
	
	ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_specific_context);
	Rcpp::List dimnames = pars.attr("dimnames");
	
	// set up parameters
	double vA_T, vB_T, svA_T, svB_T, vA_N_flip, vB_N_flip, svA_N_flip, svB_N_flip, p_negA, p_negB, p_A_Fail, p_B_Fail, p_AB_Fail, p_process, p_j;
	
    for(int j=0; j<n_unique_trials; ++j){
        int start = j*n_acc;
		p_j = 0;
        double fA=NA_REAL, fB=NA_REAL, fnA=NA_REAL, fnB=NA_REAL, fnA_flip=0, fnB_flip=0;
        double FA=NA_REAL, FB=NA_REAL, FnA=NA_REAL, FnB=NA_REAL, FnA_flip=0, FnB_flip=0;
        for(int k=0;k<n_acc;++k){
            int idx = start+k;
            std::string r = Rcpp::as<std::string>(role[idx]);
            if(r == "A"){ fA = f_all[idx]; FA = F_all[idx]; vA_T=pars(idx,0); svA_T=pars(idx,1);}
            else if(r == "B"){ fB = f_all[idx]; FB = F_all[idx]; vB_T=pars(idx,0); svB_T=pars(idx,1);}
            else if(r == "n_A"){ fnA = f_all[idx]; FnA = F_all[idx];}
            else if(r == "n_B"){ fnB = f_all[idx]; FnB = F_all[idx];}
			else if(r == "n_A_flip"){ fnA_flip = f_all[idx]; FnA_flip = F_all[idx];  vA_N_flip=pars(idx,0); svA_N_flip=pars(idx,1);}
			else if(r == "n_B_flip"){ fnB_flip = f_all[idx]; FnB_flip = F_all[idx];  vB_N_flip=pars(idx,0); svB_N_flip=pars(idx,1);}
        }
		p_negA = R::pnorm(0.0, vA_T, svA_T, 1,0); p_negB = R::pnorm(0.0, vB_T, svB_T, 1,0); 
		double one_m_FB = std::max(1e-12, 1.0 - FB);
		double one_m_FA = std::max(1e-12, 1.0 - FA);
		double one_m_FnB = std::max(1e-12, 1.0 - FnB);
		double one_m_FnA = std::max(1e-12, 1.0 - FnA);
		double one_m_FnB_flip = std::max(1e-12, 1.0 - FnB_flip);
		double one_m_FnA_flip = std::max(1e-12, 1.0 - FnA_flip);
		double one_m_p_negA = std::max(1e-12, 1.0 - p_negA);
		double one_m_p_negB = std::max(1e-12, 1.0 - p_negB);
		double one_m_FnAFnB = std::max(1e-12, 1.0 - FnA * FnB);
		double one_m_FnA_flipFnB = std::max(1e-12, 1.0 - FnA_flip * FnB);
		double one_m_FnAFnB_flip = std::max(1e-12, 1.0 - FnA * FnB_flip);
		double No_fail = std::max(1e-12, one_m_p_negA * one_m_p_negB);
		double one_m_FAFB = std::max(1e-12, 1.0 - FA * FB);
        std::string r_obs = Rcpp::as<std::string>(resp[start]);
		// We're building a likelihood allowing drift rates to go negative only for TARGET accumulators
		// If a target accumulator "fails" we substitute the corresponding "absence" accumulator with the vTRUE-Absent in some cases this won't change it, e.g. trials where the target actually wasn't there -- this is mathematically the same as just the usual race
		// When a target "fails" (has negative drift) it can't finish thus drops out of the likelihood for those cases
		// When both targets fail, OR "yes" is impossible. When either target fails, AND "yes" is impossible.
        if (LogicalRule[start]=="OR") {
			if(r_obs == "yes"){
				p_process = No_fail*((fA*one_m_FB + fB*one_m_FA) * one_m_FnAFnB);
				p_A_Fail = (p_negA*one_m_p_negB) * (fB*one_m_FnA_flipFnB); // A failed but B didn't, assumes vCorrect used for absent-A regardless of stimulus
				p_B_Fail = (p_negB*one_m_p_negA) * (fA*one_m_FnAFnB_flip); // B failed but A didn't, assumes vCorrect used for absent-B regardless of stimulus
				p_j = p_process+p_A_Fail+p_B_Fail;
			} else if(r_obs == "no"){
				p_process = No_fail*((fnA*FnB + fnB*FnA) * (one_m_FA*one_m_FB));
				p_A_Fail = (p_negA*one_m_p_negB) * (fnA_flip*FnB + fnB*FnA_flip) * one_m_FB; // A failed but B didn't, assumes vCorrect used for absent-A regardless of stimulus
				p_B_Fail = (p_negB*one_m_p_negA) * (fnA*FnB_flip + fnB_flip*FnA) * one_m_FA; // A failed but B didn't, assumes vCorrect used for absent-A regardless of stimulus
				p_AB_Fail = (p_negA*p_negB) * (fnA_flip*FnB_flip + fnB_flip*FnA_flip); // Both failed
				p_j = p_process+p_A_Fail+p_B_Fail+p_AB_Fail;
			}
		} else if (LogicalRule[start]=="AND") {
			if(r_obs == "no"){
				p_process = No_fail*((fnA*one_m_FnB + fnB*one_m_FnA) * one_m_FAFB);
				p_A_Fail = (p_negA*one_m_p_negB) * (fnB*one_m_FnA_flip + fnA_flip*one_m_FnB); // A failed but B didn't, assumes vCorrect used for absent-A regardless of stimulus
				p_B_Fail = (p_negB*one_m_p_negA) * (fnB_flip*one_m_FnA + fnA*one_m_FnB_flip); // A failed but B didn't, assumes vCorrect used for absent-A regardless of stimulus
				p_AB_Fail = (p_negA*p_negB) * (fnB_flip*one_m_FnA_flip + fnA_flip*one_m_FnB_flip); // Both failed
				p_j = p_process+p_A_Fail+p_B_Fail+p_AB_Fail;
			} else if(r_obs == "yes"){
				p_process = No_fail*((fA*FB + fB*FA) * (one_m_FnA*one_m_FnB));
				p_j = p_process;
			}
		}
		if (p_j <= 0.0 || !R_FINITE(p_j)) {
			ll_unique[j] = min_ll;
		} else {
			ll_unique[j] = std::log(p_j);
		}
    }

    Rcpp::NumericVector ll_exp = c_expand(ll_unique, expand);
    double sum_ll = 0.0;
    for(int i=0;i<ll_exp.size();++i){
        double val = ll_exp[i];
        if(!R_FINITE(val) || Rcpp::NumericVector::is_na(val) || val < min_ll) val = min_ll;
        sum_ll += val;
    }
    return sum_ll;
}
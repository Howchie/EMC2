#include "trend.h"

// Now accepts the full parameter matrix `pars_full` so we can use par_input columns as inputs too.
// Passes all inputs (covariates + par_input) to kernel in one call; kernel sums across columns.
// [[Rcpp::export]]
NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param,
                             NumericMatrix trend_pars, NumericMatrix pars_full,
                             bool return_kernel = false) {
  String kernel = as<String>(trend["kernel"]);
  String base = as<String>(trend["base"]);
  // Extract optional custom pointer attribute if present
  SEXP custom_ptr = R_NilValue;
  if (kernel == "custom") {
    custom_ptr = trend.attr("custom_ptr");
  }
  CharacterVector covnames;
  if (trend.containsElementNamed("covariate") && !Rf_isNull(trend["covariate"])) {
    covnames = trend["covariate"];
  } else {
    covnames = CharacterVector(0);
  }
  CharacterVector par_input;
  if (trend.containsElementNamed("par_input") && !Rf_isNull(trend["par_input"])) {
    par_input = trend["par_input"];
  } else {
    par_input = CharacterVector(0);
  }

  bool ffill_na;
  if (trend.containsElementNamed("ffill_na") && !Rf_isNull(trend["ffill_na"])) {
    ffill_na = as<bool>(trend["ffill_na"]);
  } else {
    ffill_na = false;
  }

  // Initialize output vector with zeros
  int n_trials = param.length();
  NumericVector out(n_trials, 0.0);

  int n_base_pars = 0;
  if (base == "lin" || base == "exp_lin" || base == "centered") {
    n_base_pars = 1;
  }
  // Check for covariate maps
  CharacterVector map_names;
  List map_list;
  int n_maps = 0;
  bool has_maps = trend.containsElementNamed("map") && !Rf_isNull(trend["map"]);
  if (has_maps) {
    map_list = trend["map"];
    map_names = map_list.names();
    n_maps = map_names.size();
    n_base_pars = n_maps * n_base_pars; // find n maps
  }

  // extract kernel parameters
  int n_trend_pars = trend_pars.ncol();
  NumericMatrix kernel_pars(n_trials, n_trend_pars - n_base_pars);
  for (int j = n_base_pars; j < n_trend_pars; j++) {
    kernel_pars(_, j - n_base_pars) = trend_pars(_, j);
  }

  // Build input matrix from covariates and par_input columns
  int n_cov = covnames.size();
  // Keep only par_input columns that actually exist in pars_full (if provided)
  CharacterVector pars_full_names = colnames(pars_full);
  std::vector<std::string> par_in_keep;
  for (int i = 0; i < par_input.size(); i++) {
    std::string nm = Rcpp::as<std::string>(par_input[i]);
    bool found = false;
    for (int j = 0; j < pars_full_names.size(); j++) {
      if (nm == Rcpp::as<std::string>(pars_full_names[j])) {
        found = true;
        break;
      }
    }
    if (found) par_in_keep.push_back(nm);
  }
  int n_par_in = (int)par_in_keep.size();

  int n_inputs = n_cov + n_par_in;
  // If no inputs provided, keep out as zeros (base will handle accordingly)
  if (n_inputs > 0) {
    NumericMatrix input_all(n_trials, n_inputs);
    // Fill covariate columns first
    for (int i = 0; i < n_cov; i++) {
      String cur_cov = covnames[i];
      NumericVector covariate = as<NumericVector>(data[cur_cov]);
      input_all(_, i) = covariate;
    }
    // Then par_input columns
    if (n_par_in > 0) {
      CharacterVector par_in_keep_cv(n_par_in);
      for (int i = 0; i < n_par_in; i++) par_in_keep_cv[i] = par_in_keep[i];
      NumericMatrix pin = submat_rcpp_col_by_names(pars_full, par_in_keep_cv);
      for (int j = 0; j < n_par_in; j++) {
        input_all(_, n_cov + j) = pin(_, j);
      }
    }

    // Build optional first-level mask and call kernel (handles compression/expansion internally)
    LogicalVector first_level;
    if (trend.containsElementNamed("at") && !Rf_isNull(trend["at"])) {
      String at_name = trend["at"];
      SEXP at_col = data[at_name];
      if (!Rf_inherits(at_col, "factor")) stop("'at' column must be a factor");
      IntegerVector f = as<IntegerVector>(at_col);
      first_level = (f == 1);
    } else {
      first_level = LogicalVector(0);
    }

    // Run kernel
    NumericMatrix kernel_out = run_kernel_rcpp(kernel_pars, kernel, input_all, custom_ptr, first_level, ffill_na);
    if (return_kernel) return kernel_out;

    int n_rows = kernel_out.nrow();
    int n_cols = kernel_out.ncol();

    // Outer loop over maps
    List covariate_maps = data.attr("covariate_maps");
    int n_loops = has_maps ? n_maps : 1;
    for (int map_n = 0; map_n < n_loops; ++map_n) {
      NumericMatrix covariate_map;
      if (has_maps) {
        std::string map_name = as<std::string>(map_names[map_n]);
        if (Rf_isNull(covariate_maps[map_name])) {
          stop("No map named '%s' found in covariate_maps", map_name.c_str());
        }
        covariate_map = as<NumericMatrix>(covariate_maps[map_name]);
      }

      // Loop over columns
      for (int c = 0; c < n_cols; ++c) {
        // Loop over rows
        for (int r = 0; r < n_rows; ++r) {
          double contrib = kernel_out(r, c);

          // Multiply by covariate map if it exists
          if (has_maps) contrib *= covariate_map(r, c);

          // Apply base parameter
          if (base == "lin" || base == "exp_lin") {
            contrib *= trend_pars(r, map_n);
          } else if (base == "centered") {
            contrib = (contrib - 0.5) * trend_pars(r, map_n);
          } // "add" leaves contrib unchanged

          // Accumulate into output
          out(r) += contrib;
        }
      }
    }
  }

  // Add parameter to obtain final summed input
  if (base == "lin") {
    out += param;
  } else if (base == "exp_lin") {
    out += exp(param);
  } else if (base == "centered") {
    out += param;
  } else if (base == "add") {
    out += param;
  }
  return out;
}

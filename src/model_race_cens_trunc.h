#ifndef MODEL_RACE_CENS_TRUNC_H
#define MODEL_RACE_CENS_TRUNC_H

#include <RcppArmadillo.h>
#include "utility_functions.h" // For any existing utilities we might use

// Function pointer types for model's dfun and pfun
// Added LogicalVector is_ok for bound checking and void* context for extra args
typedef Rcpp::NumericVector (*RacePdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
                                          Rcpp::LogicalVector is_ok,
                                          bool log_p,
                                          void* model_specific_context);
typedef Rcpp::NumericVector (*RaceCdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
                                          Rcpp::LogicalVector is_ok,
                                          bool log_p,
                                          void* model_specific_context);

// Forward declaration for the main exported function
// Returns a single double: the total log-likelihood for one particle.
double c_log_likelihood_race_cens_trunc(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, all unique trial conditions
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions
    RacePdfFun model_dfun,                  // Pointer to the model's C++ PDF
    RaceCdfFun model_pfun,                  // Pointer to the model's C++ CDF
    double min_ll,
    const Rcpp::LogicalVector& ok_params,   // Parameter validity for this particle's pars matrix
    int n_acc,                              // Number of accumulators in the race
    const Rcpp::IntegerVector& expand_vec,  // Vector for expanding unique LLs to full trial count
    void* model_context_for_funcs           // Context for model_dfun/model_pfun
) {
    // Extract needed columns from dadm
    Rcpp::NumericVector rts = dadm["rt"];
    Rcpp::IntegerVector R_idxs = dadm["R"]; // Assuming 'R' column holds 1-based winner index or NA
    Rcpp::NumericVector trial_Ns = dadm["N"]; // Number of times each unique trial condition occurs

    // Truncation and censoring limits (assuming these are attributes of dadm or passed differently)
    // For now, access them as attributes if they exist, otherwise use defaults.
    double LT = dadm.hasAttribute("LT") ? Rcpp::as<double>(dadm.attr("LT")) : 0.0;
    double UT = dadm.hasAttribute("UT") ? Rcpp::as<double>(dadm.attr("UT")) : R_PosInf;
    double LC = dadm.hasAttribute("LC") ? Rcpp::as<double>(dadm.attr("LC")) : 0.0;
    double UC = dadm.hasAttribute("UC") ? Rcpp::as<double>(dadm.attr("UC")) : R_PosInf;

    int n_unique_trials = dadm.nrows() / n_acc; // Assuming dadm is structured with n_acc rows per unique trial condition
    Rcpp::NumericVector ll_unique(n_unique_trials, min_ll); // Initialize with min_ll

    std::vector<int> finite_rt_trial_indices;
    std::vector<int> other_trial_indices; // For NA, Inf, -Inf RTs

    for (int j = 0; j < n_unique_trials; ++j) {
        int first_row_for_trial = j * n_acc;
        double rt_j = rts[first_row_for_trial]; // Assuming RT is same for all acc rows of a unique trial

        if (R_IsNA(rt_j) || !R_finite(rt_j)) {
            other_trial_indices.push_back(j);
        } else {
            finite_rt_trial_indices.push_back(j);
        }
    }

    double total_ll = 0.0;

    // --- Part 1: Batch process finite RT trials (Observed RTs) ---
    if (!finite_rt_trial_indices.empty()) {
        // TODO: Implement vectorized calculation for finite RT trials
        // This will involve:
        // 1. Creating sub-matrices/vectors of pars and dadm for these trials.
        // 2. Calling a vectorized version of f_race_integrand_cpp (or similar logic).
        // 3. Applying truncation correction.
        // For now, placeholder:
        for (int trial_idx : finite_rt_trial_indices) {
            // This is a simplified placeholder, actual logic will be more complex
            // and use the existing helper functions after adaptation.
            int first_row_for_trial = trial_idx * n_acc;
            if (!ok_params[first_row_for_trial]) { // Check ok_params for the first row of this trial block
                ll_unique[trial_idx] = min_ll;
                continue;
            }

            Rcpp::NumericMatrix pars_condition_j_all_acc(n_acc, pars.ncol());
            for(int i = 0; i < n_acc; ++i) {
                pars_condition_j_all_acc.row(i) = pars.row(first_row_for_trial + i);
            }

            double rt_j = rts[first_row_for_trial];
            int R_j_idx = R_idxs[first_row_for_trial]; // 1-based winner index

            if (R_IsNA(R_j_idx)) { // Should not happen for finite RT if data is clean
                 ll_unique[trial_idx] = min_ll;
                 continue;
            }

            double prob_density = 0.0;
            if (rt_j >= LT && rt_j <= UT) { // Check if RT is within truncation bounds
                Rcpp::NumericMatrix pars_ordered_obs = order_pars_for_winner_cpp(pars_condition_j_all_acc, R_j_idx, n_acc);
                prob_density = f_race_integrand_cpp(rt_j, pars_ordered_obs, model_dfun, model_pfun, n_acc, model_context_for_funcs);

                if (prob_density > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, 1e-6, model_context_for_funcs);
                    if (R_IsNA(trunc_cf) || !R_finite(trunc_cf) || trunc_cf < 0) {
                        prob_density = 0;
                    } else {
                        prob_density *= trunc_cf;
                    }
                }
            }

            if (prob_density > std::numeric_limits<double>::epsilon()) {
                ll_unique[trial_idx] = std::log(prob_density);
            } else {
                ll_unique[trial_idx] = min_ll;
            }
        }
    }

    // --- Part 2: Process other trials (Infinite RTs, NA RTs - Censored/Missing) ---
    // This part will largely follow the structure of the original R loop, using helper functions.
    for (int trial_idx : other_trial_indices) {
        int first_row_for_trial = trial_idx * n_acc;
        if (!ok_params[first_row_for_trial]) { // Check ok_params for the first row
            ll_unique[trial_idx] = min_ll;
            continue;
        }

        Rcpp::NumericMatrix pars_condition_j_all_acc(n_acc, pars.ncol());
        for(int i = 0; i < n_acc; ++i) {
            pars_condition_j_all_acc.row(i) = pars.row(first_row_for_trial + i);
        }

        double rt_j = rts[first_row_for_trial]; // This will be Inf, -Inf, or NA
        int R_j_idx = R_idxs[first_row_for_trial]; // Winner index (1-based) or NA

        double prob_mass = 0.0;
        double epsilon_integration = 1e-6; // Tolerance for integration

        if (rt_j == R_NegInf) { // Fast censoring (rt = -Inf)
            if (!R_IsNA(R_j_idx)) {
                prob_mass = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                if (prob_mass > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, epsilon_integration, model_context_for_funcs);
                    if (R_IsNA(trunc_cf) || !R_finite(trunc_cf) || trunc_cf < 0) prob_mass = 0; else prob_mass *= trunc_cf;
                }
            } else { // Winner unknown
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                    if (p_k > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, epsilon_integration, model_context_for_funcs);
                        if (!R_IsNA(trunc_cf_k) && R_finite(trunc_cf_k) && trunc_cf_k >=0) prob_mass += (p_k * trunc_cf_k);
                    }
                }
            }
        } else if (rt_j == R_PosInf) { // Slow censoring (rt = Inf)
             if (!R_IsNA(R_j_idx)) {
                prob_mass = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                 if (prob_mass > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, epsilon_integration, model_context_for_funcs);
                    if (R_IsNA(trunc_cf) || !R_finite(trunc_cf) || trunc_cf < 0) prob_mass = 0; else prob_mass *= trunc_cf;
                }
            } else { // Winner unknown
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double p_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                     if (p_k > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, epsilon_integration, model_context_for_funcs);
                        if (!R_IsNA(trunc_cf_k) && R_finite(trunc_cf_k) && trunc_cf_k >=0) prob_mass += (p_k * trunc_cf_k);
                    }
                }
            }
        } else if (R_IsNA(rt_j)) { // Missing RT (NA) - implies censoring at both ends
            if (!R_IsNA(R_j_idx)) {
                double prob_L = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                double prob_U = integrate_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                prob_mass = prob_L + prob_U;
                 if (prob_mass > 0) {
                    double trunc_cf = get_trunc_corr_factor_for_kth_winner_cpp(R_j_idx, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, epsilon_integration, model_context_for_funcs);
                    if (R_IsNA(trunc_cf) || !R_finite(trunc_cf) || trunc_cf < 0) prob_mass = 0; else prob_mass *= trunc_cf;
                }
            } else { // Winner unknown
                for (int k_win = 1; k_win <= n_acc; ++k_win) {
                    double prob_L_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, LT, LC, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                    double prob_U_k = integrate_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, UC, UT, model_dfun, model_pfun, n_acc, epsilon_integration, model_context_for_funcs);
                    double p_k_sum = prob_L_k + prob_U_k;
                    if (p_k_sum > 0) {
                        double trunc_cf_k = get_trunc_corr_factor_for_kth_winner_cpp(k_win, pars_condition_j_all_acc, model_dfun, model_pfun, LT, UT, n_acc, epsilon_integration, model_context_for_funcs);
                        if (!R_IsNA(trunc_cf_k) && R_finite(trunc_cf_k) && trunc_cf_k >=0) prob_mass += (p_k_sum * trunc_cf_k);
                    }
                }
            }
        }

        if (prob_mass > std::numeric_limits<double>::epsilon()) {
            ll_unique[trial_idx] = std::log(prob_mass);
        } else {
            ll_unique[trial_idx] = min_ll;
        }
    }

    // --- Combine LLs from all unique trials, considering expansion/weighting ---
    if (expand_vec.length() > 0) { // If expand_vec is provided and not empty
        for (int i = 0; i < expand_vec.length(); ++i) {
            // expand_vec is 1-based from R, adjust to 0-based for C++ vector access
            int unique_idx = expand_vec[i] - 1;
            if (unique_idx >= 0 && unique_idx < n_unique_trials) {
                 // Ensure ll_unique[unique_idx] is not min_ll before adding, to prevent sum(lots of min_ll)
                 // if (ll_unique[unique_idx] > min_ll - std::numeric_limits<double>::epsilon()) {
                    total_ll += ll_unique[unique_idx];
                 // } else {
                 //    total_ll += min_ll; // Add min_ll if the unique trial's ll was min_ll
                 // }
            } else {
                // Handle error or warning: expand_vec contains invalid index
                // For now, just skip or add min_ll
                total_ll += min_ll;
            }
        }
    } else if (trial_Ns.length() == n_unique_trials) { // If N (weights) are provided
        for (int j = 0; j < n_unique_trials; ++j) {
            // if (ll_unique[j] > min_ll - std::numeric_limits<double>::epsilon()) {
                total_ll += ll_unique[j] * trial_Ns[j];
            // } else {
            //    total_ll += min_ll * trial_Ns[j];
            // }
        }
    } else { // Default: sum ll_unique directly (each unique trial counted once)
        for (int j = 0; j < n_unique_trials; ++j) {
            // if (ll_unique[j] > min_ll - std::numeric_limits<double>::epsilon()) {
                total_ll += ll_unique[j];
            // } else {
            //    total_ll += min_ll;
            // }
        }
    }

    // Ensure total_ll is not worse than n_total_trials * min_ll if all fail
    // This depends on how expand_vec or N are used.
    // If expand_vec is used, n_total_trials = expand_vec.length().
    // If N is used, n_total_trials = sum(N).
    // If neither, n_total_trials = n_unique_trials.
    // For now, this simple check might be okay, but a more robust lower bound might be needed.
    // double min_total_ll_possible = (expand_vec.length() > 0 ? expand_vec.length() : (trial_Ns.length() == n_unique_trials ? Rcpp::sum(trial_Ns) : n_unique_trials)) * min_ll;
    // if (total_ll < min_total_ll_possible && min_total_ll_possible < 0) total_ll = min_total_ll_possible;


    return total_ll;
}


// --- Internal Helper Function Declarations (not directly exported to R via Rcpp::export) ---

// Integrand function equivalent to R's f_race_integrand
// p_trial_this_winner_first: (n_accumulators x n_model_params)
// model_dfun, model_pfun: function pointers to the model's density and CDF
double f_race_integrand_cpp(
    double t,
    const Rcpp::NumericMatrix& p_trial_this_winner_first, // Parameters for the winner (1st row) and losers
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    void* model_specific_context // Added context
);

// Numerical integration helper
// k_winner_idx: 1-based index of the accumulator treated as winner
// p_all_acc: parameters for all accumulators for the current trial condition
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    double low,
    double upp,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    double epsilon,
    void* model_specific_context
);

// Truncation correction factor helper
// Returns the factor itself (not log)
double get_trunc_corr_factor_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    double LT, double UT,
    int n_acc,
    double epsilon,
    void* model_specific_context
);

// Helper to order parameters (winner first)
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc, // n_acc x n_params
    int k_idx, // 1-based index of the winner
    int n_acc
);


#endif // MODEL_RACE_CENS_TRUNC_H
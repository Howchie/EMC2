// Top-level hUVSD log-likelihood, moved verbatim out of particle_ll.cpp.
//
// The scalar kernel log_likelihood_huvsd_single stays inline in model_SDT.h;
// only the per-data-frame driver moves here. Signal/yes level resolution,
// parameter column lookup, is_ok censoring and the min_ll floor are unchanged.

#include "model_SDT.h"

#include <string>

double c_log_likelihood_huvsd(Rcpp::NumericMatrix pars, Rcpp::DataFrame data,

                              const int n_trials, Rcpp::IntegerVector expand,
                              double min_ll, Rcpp::LogicalVector is_ok,
                              Rcpp::NumericVector* trial_ll_out) {
    Rcpp::IntegerVector S = data["S"];
    Rcpp::IntegerVector R = data["R"];
    Rcpp::CharacterVector S_levels = S.attr("levels");
    Rcpp::CharacterVector R_levels = R.attr("levels");
    
    int signal_level = 2; 
    for(int i=0; i<S_levels.size(); ++i) {
        std::string lev = Rcpp::as<std::string>(S_levels[i]);
        if(lev == "Signal" || lev == "S" || lev == "signal") { signal_level = i+1; break; }
    }
    int yes_level = 2;
    for(int i=0; i<R_levels.size(); ++i) {
        std::string lev = Rcpp::as<std::string>(R_levels[i]);
        if(lev == "Yes" || lev == "Y" || lev == "yes" || lev == "Hit" || lev == "hit") { yes_level = i+1; break; }
    }

    Rcpp::CharacterVector p_types = colnames(pars);
    int d_idx = -1, c_idx = -1, sd_idx = -1;
    for(int j=0; j<p_types.size(); ++j) {
        if(p_types[j] == "d") d_idx = j;
        else if(p_types[j] == "c") c_idx = j;
        else if(p_types[j] == "sd") sd_idx = j;
    }
    if(d_idx == -1 || c_idx == -1 || sd_idx == -1) Rcpp::stop("hUVSD model requires parameters d, c, and sd");

    const int n_out = (expand.length() > 0) ? expand.length() : n_trials;
    double total_ll = 0.0;
    for(int j=0; j<n_out; ++j) {
        int row = (expand.length() > 0) ? (expand[j] - 1) : j;
        double ll = min_ll;
        if(is_ok[row]) {
            ll = log_likelihood_huvsd_single(pars(row, d_idx), pars(row, c_idx), pars(row, sd_idx),
                                             S[row] == signal_level, R[row] == yes_level, min_ll);
        }
        if(trial_ll_out) (*trial_ll_out)[j] = ll;
        total_ll += ll;
    }
    return total_ll;
}

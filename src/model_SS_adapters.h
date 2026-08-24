#ifndef EMC2_MODEL_SS_ADAPTERS_H
#define EMC2_MODEL_SS_ADAPTERS_H

#include <RcppArmadillo.h>
#include <cstdint>
#include <string>
#include <vector>

enum SsTrialKind : uint8_t {
  SS_TRIAL_LC_CENSORED = 0,  // rt == -Inf: integrate response mass over [LT, LC)
  SS_TRIAL_GO_NR_FREE,       // go trial, no response, no deadline
  SS_TRIAL_GO_NR_DEADLINE,   // go trial, no response by UC
  SS_TRIAL_STOP_NR_FREE,     // stop trial, no response, no deadline
  SS_TRIAL_STOP_NR_DEADLINE, // stop trial, no response by UC
  SS_TRIAL_GO_RESP,          // response on a go trial
  SS_TRIAL_STOP_GO_RESP,     // go response on a stop trial
  SS_TRIAL_STOP_ST_RESP      // ST (or unmatched) response on a stop trial
};

struct SsSharedState {
  bool valid = false;
  int n_trials = 0;                       // unique trials
  int n_acc = 0;                          // accumulators per trial
  std::vector<uint8_t> kind;              // SsTrialKind per trial
  std::vector<double> rt, SSD, UC, LC, LT;
  std::vector<int> resp_code;             // R code, NA_INTEGER if no response
  std::vector<int> n_go, n_st, first_go_row;   // first_go_row: offset in trial
  std::vector<uint8_t> is_go_row, is_st_row, win_row;  // n_trials * n_acc
  std::vector<int> lR_row;                // lR codes, n_trials * n_acc
};

SsSharedState build_ss_shared_state(const Rcpp::DataFrame& data,
                                    int n_trials_ss, int n_acc);

struct SsRawWorkspace {
  std::vector<double> acc;     // n_acc * SS_ACC_STRIDE
  std::vector<double> acc_go;  // n_go  * SS_ACC_STRIDE
  std::vector<double> lls;     // per unique trial
};

double c_log_likelihood_ss(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame data,
    const int n_trials,
    Rcpp::IntegerVector expand,
    double min_ll,
    Rcpp::LogicalVector is_ok,
    Rcpp::NumericVector (*go_lpdf_ptr)(Rcpp::NumericVector, Rcpp::NumericMatrix,
                                       Rcpp::LogicalVector, double),
    Rcpp::NumericVector (*go_lccdf_ptr)(Rcpp::NumericVector, Rcpp::NumericMatrix,
                                        Rcpp::LogicalVector, double),
    double (*stop_logsurv_ptr)(double, Rcpp::NumericMatrix),
    double (*stop_success_ptr)(double, Rcpp::NumericMatrix, double, double, int,
                               double, double, double, double),
    int idx_tf,
    int idx_gf);

struct SsRawModel;
double c_log_likelihood_ss_pt(
    const double* const* cols,
    const SsSharedState& ss,
    const SsRawModel& M,
    const Rcpp::LogicalVector& is_ok,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    SsRawWorkspace& ws);

struct SSModelAdapter {
  Rcpp::NumericVector (*go_lpdf_ptr)(Rcpp::NumericVector, Rcpp::NumericMatrix,
                                     Rcpp::LogicalVector, double);
  Rcpp::NumericVector (*go_lccdf_ptr)(Rcpp::NumericVector, Rcpp::NumericMatrix,
                                      Rcpp::LogicalVector, double);
  double (*stop_logsurv_ptr)(double, Rcpp::NumericMatrix);
  double (*stop_success_ptr)(double, Rcpp::NumericMatrix, double, double, int,
                             double, double, double, double);
  int idx_tf;
  int idx_gf;
};

SSModelAdapter resolve_ss_adapter(const std::string& type_std);

#endif

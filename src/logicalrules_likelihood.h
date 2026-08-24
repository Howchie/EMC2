#ifndef EMC2_LOGICALRULES_LIKELIHOOD_H
#define EMC2_LOGICALRULES_LIKELIHOOD_H

#include <Rcpp.h>
#include <cstddef>
#include <vector>

#include "race_contract.h"

struct GuessKernel;

struct LogicalRulesCellKey {
  int idxA, idxB, idxnA, idxnB, idxNogo;
  int rule_code, cond_code;
  double LT, UT;

  bool operator==(const LogicalRulesCellKey& o) const {
    return idxA == o.idxA && idxB == o.idxB && idxnA == o.idxnA && idxnB == o.idxnB &&
           idxNogo == o.idxNogo && rule_code == o.rule_code && cond_code == o.cond_code &&
           LT == o.LT && UT == o.UT;
  }
};

struct LogicalRulesSharedState {
  bool valid = false;
  int n_trials = 0;
  int n_acc = 0;
  int n_unique_trials = 0;
  bool has_nogo = false;
  std::vector<int> idxA;
  std::vector<int> idxB;
  std::vector<int> idxnA;
  std::vector<int> idxnB;
  std::vector<int> idxNogo;
  // Rule code: 1=OR, 2=AND, 3=XOR, 4=ID, 5=OR_DETECTION_ANALYTIC, 6=OR_DETECTION_GNG
  std::vector<int> rule_code;
  // Stimulus code for detection rules: 0=NN, 1=AN/A, 2=NB/B, 3=AB
  std::vector<int> cond_code;
  // Response code: 0=missing, 1=yes, 2=no/nogo, 3=NN, 4=AN, 5=NB, 6=AB
  std::vector<int> resp_code;
  std::vector<double> rt_unique;
  std::vector<double> LT_unique;
  std::vector<double> UT_unique;
  std::vector<double> LC_unique;
  std::vector<double> UC_unique;
  // Per-row RTs (length n_trials) for raw dfun/pfun kernels.
  std::vector<double> rt_by_row;

  // Pre-computed design cells and trial-to-cell mapping
  int n_cells = 0;
  std::vector<int> cell_id;                     // length n_unique_trials
  std::vector<LogicalRulesCellKey> cells;       // length n_cells
  std::vector<unsigned char> ch_eq_cell;        // per unique trial, channel equality

  // Per-particle cache for truncation normalisers and capacity denominators.
  mutable std::vector<double> cell_log_z;       // length n_cells
  mutable std::vector<double> cell_log_den_cap; // length n_cells

  void clear_particle_cache() const {
    cell_log_z.assign(static_cast<size_t>(n_cells), NA_REAL);
    cell_log_den_cap.assign(static_cast<size_t>(n_cells), NA_REAL);
  }
};

LogicalRulesSharedState build_logicalrules_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_acc,
    bool capacity = false);

double c_log_likelihood_logicalrules(
    const double* const* pars_cols,
    int n_par,
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
    Rcpp::NumericVector* trial_ll_out = nullptr,
    int kappa_col = -1,
    int tau_col = -1,
    int pc_col = -1,
    int pg_col = -1,
    const GuessKernel* guess = nullptr);

Rcpp::List lr_capacity_counter_values();
void lr_capacity_counters_reset();

#endif // EMC2_LOGICALRULES_LIKELIHOOD_H

#ifndef EMC2_CORRELATED_LIKELIHOOD_H
#define EMC2_CORRELATED_LIKELIHOOD_H

#include <Rcpp.h>
#include <cstdint>
#include <functional>
#include <vector>

#include "race_contract.h"
#include "contaminant_mixture.h"
#include "drift_factor.h"

struct ParamTable;
struct RaceSharedState;

enum class BAwLCorrRoute : uint8_t {
  ordinary = 0,       // 0/1 nonzero loadings: independent race machinery
  exact_pair,         // two loaded rows, no clocks, exact rectangle kernel
  numeric_pair,       // exact route reported unstable: direct 1-D integration
  gh_no_clock,        // shared-factor GH with the fused no-clock evaluator
  gh_generic_clock,   // shared-factor GH through the generic race evaluator
  invalid             // malformed rows: trial value floors to min_ll
};

struct BAwLCorrTrialLayout {
  BAwLCorrRoute route = BAwLCorrRoute::gh_no_clock;
  int n_active = 0;             // isok && RACE-active rows
  int n_loaded = 0;             // active rows with |effective rho| > threshold
  int loaded_row[2] = {-1, -1}; // row indices, meaningful while n_loaded <= 2
  int winner_row = -1;
  bool winner_loaded = false;
  bool any_bad_row = false;     // active row failed isok/rho validity
};

struct CorrDriftSharedState {
  bool valid = false;
  int n_trials = 0;
  int n_lR = 0;
  int n_unique = 0;
  int n_out = 0;
  int n_par = 0;

  // R holders keeping backing memory alive
  Rcpp::NumericVector rt, LT, UT, LC, UC;
  Rcpp::IntegerVector race_nacc;
  Rcpp::LogicalVector race_mask;
  bool has_RACE = false;
  bool has_truncation = false;

  // Data-fixed arrays
  std::vector<unsigned char> role_correct;    // per row, lM role
  std::vector<int> winner_row;                // per unique trial (-1 if none)
  std::vector<unsigned char> has_trunc_trial; // per unique trial
  std::vector<int> j_to_i;                    // unique trial -> output row

  // Direct ParamTable column pointers in keep_names (p_types) order; the
  // base storage is refilled in place per particle, so these stay valid.
  std::vector<const double*> cols;
  int pc_col = -1;    // pContaminant column, -1 when absent
  int pg_col = -1;    // pGuess column, -1 when absent
  GuessKernel guess;  // uniform guess kernel; inactive when pGuess is unused
  int rho_col = -1;

  // Where the model keeps the drift mean and SD, and whether it has an
  // affine-in-drift survivor.  BAwL does, so it gets the exact rectangle and
  // fused fast routes; a Wald-kernel model does not, and takes the generic
  // shared-factor route for every loaded trial.  See drift_factor.h.
  DriftFactorColumns factor_cols;
  bool generic_only = false;

  // Reused per-particle scratch
  std::vector<double> effective_rho;
  std::vector<BAwLCorrTrialLayout> layout;

  // Per-particle memoization of parameter-cell quantities.  Unique trials
  // multiply design cells by distinct RTs, so the positive-drift orthant
  // normalizer q_AB and the truncation normalizer log Z — which depend only
  // on the cell parameters (and the data-constant LT/UT window) — are
  // otherwise recomputed once per unique trial instead of once per cell.
  // Exact double keys, linear scan, capped size; cleared for each particle
  // because the ParamTable base storage is refilled in place.  Mutable so the
  // const-view trial evaluators can memoize.
  struct QabCacheEntry { double key[5]; double value; };
  struct ZCacheEntry { double key[15]; double log_z; };
  mutable std::vector<QabCacheEntry> qab_cache;
  mutable std::vector<ZCacheEntry> z_cache;
};

CorrDriftSharedState build_corr_drift_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names);

double c_log_likelihood_corr_drift(
    CorrDriftSharedState& cshared,
    Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    const int n_trials,
    Rcpp::LogicalVector winner,
    Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector isok,
    int n_lR,
    void* model_context_for_funcs,
    bool all_finite_trials,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,
    RaceSharedState* shared,
    Rcpp::NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize);

struct RDMSWTNCorrTrialLayout {
  int pair_row[2] = {-1, -1};
  int n_nonzero = 0;
  int winner_row = -1;
  bool ordinary = true;
  bool params_ok = true;
};

struct RDMSWTNCorrSharedState {
  bool valid = false;
  int n_trials = 0;
  int n_lR = 0;
  int n_unique = 0;
  int n_out = 0;
  int n_par = 0;
  int rho_col = -1;
  int sv_col = -1;    // sv column, needed by the posdrift = FALSE gate
  int pc_col = -1;
  int pg_col = -1;    // pGuess column, -1 when absent
  GuessKernel guess;  // uniform guess kernel; inactive when pGuess is unused
  bool has_RACE = false;
  Rcpp::NumericVector rt, LT, UT, LC, UC;
  Rcpp::IntegerVector lR_codes;
  Rcpp::IntegerVector race_nacc;
  Rcpp::LogicalVector race_mask;
  std::vector<int> winner_row;
  std::vector<const double*> cols;
  std::vector<RDMSWTNCorrTrialLayout> layout;
};

RDMSWTNCorrSharedState build_rdmswtn_corr_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names);

double c_log_likelihood_rdmswtn_correlated(
    RDMSWTNCorrSharedState& cshared, Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, const int n_trials,
    Rcpp::LogicalVector winner, Rcpp::IntegerVector expand, double min_ll,
    const Rcpp::LogicalVector isok, int n_lR,
    void* model_context_for_funcs, bool all_finite_trials,
    RaceRawFun model_dfun_raw, RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t, RaceSharedState* shared,
    Rcpp::NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize);

#endif // EMC2_CORRELATED_LIKELIHOOD_H

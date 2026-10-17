#ifndef EMC2_FRQ_SHARED_LIKELIHOOD_H
#define EMC2_FRQ_SHARED_LIKELIHOOD_H

#include <Rcpp.h>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "race_contract.h"
#include "contaminant_mixture.h"

struct ParamTable;

// FRQ(rate_sharing = "shared"): one Gamma rate multiplier per trial, shared by
// every accumulator (src/frq_shared.h).  The race does not factorise, so it
// has its own trial evaluator rather than the raw-kernel product.
struct FrqSharedState {
  bool valid = false;
  bool accel_chart = false;
  int n_trials = 0;
  int n_lR = 0;
  int n_unique = 0;
  int n_out = 0;
  int n_par = 0;
  int pc_col = -1;
  int pg_col = -1;
  int omega_col = -1;
  GuessKernel guess;
  bool has_RACE = false;
  double rt_resolution = 0.0;
  Rcpp::NumericVector rt, LT, UT, LC, UC;
  Rcpp::IntegerVector race_nacc;
  Rcpp::LogicalVector race_mask;
  std::vector<int> winner_row;
  std::vector<const double*> cols;

  // Truncation normaliser per parameter tuple, cleared per particle: unique
  // trials repeat design cells, and the normaliser depends on the cell only.
  struct ZKeyHash {
    std::size_t operator()(const std::vector<std::uint64_t>& k) const {
      std::size_t h = 0;
      for (std::uint64_t b : k)
        h ^= std::hash<std::uint64_t>{}(b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      return h;
    }
  };
  std::unordered_map<std::vector<std::uint64_t>, double, ZKeyHash> z_cache;
};

FrqSharedState build_frq_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names);

double c_log_likelihood_frq_shared(
    FrqSharedState& s, const Rcpp::DataFrame& dadm, double min_ll,
    const Rcpp::LogicalVector& isok, const Rcpp::IntegerVector& expand,
    ContextForRaceModels* ctx, Rcpp::NumericVector* trial_ll_out);

#endif  // EMC2_FRQ_SHARED_LIKELIHOOD_H

#ifndef EMC2_RACE_INTEGRANDS_H
#define EMC2_RACE_INTEGRANDS_H

#include "race_contract.h"
#include "gsl_utils.h"

#include <cstddef>
#include <unordered_map>
#include <vector>

// A generic race endpoint cache.  It groups complete trial parameter blocks so
// the model-specific logS_at_t callback only evaluates one representative per
// exact parameter key.  The grouping is built once per particle and reused for
// both truncation endpoints.
struct RaceEndpointGroupCache {
  std::vector<int> group_id;
  std::vector<int> representative;
  std::unordered_map<std::size_t, std::vector<int>> hash_groups;
  std::vector<double> compact_cols;
  std::vector<int> compact_isok;
  std::vector<int> compact_mask;
  std::vector<const double*> compact_col_ptrs;
  int n_included = 0;
  bool prepared = false;

  void new_particle() {
    group_id.clear();
    representative.clear();
    hash_groups.clear();
    compact_cols.clear();
    compact_isok.clear();
    compact_mask.clear();
    compact_col_ptrs.clear();
    n_included = 0;
    prepared = false;
  }
};

void race_endpoint_prepare_groups(
    RaceEndpointGroupCache& cache, const double* const* cols,
    int n_unique_trials, int n_lR, int n_par,
    const int* include_mask, const int* isok, int skip_a, int skip_b);

double log_survivor_rowmajor(double t,
                             const double* pars_rowmajor,
                             const int* isok_int,
                             int n_lR,
                             int n_par,
                             RaceCdf1Fun cdf1,
                             void* ctx);

double log_cdf_rowmajor(double t,
                        const double* pars_rowmajor,
                        const int* isok_int,
                        int n_lR,
                        int n_par,
                        RaceCdf1Fun cdf1,
                        void* ctx);

double log_min_density_rowmajor(double t,
                                const double* pars_rowmajor,
                                const int* isok_int,
                                int n_lR,
                                int n_par,
                                RacePdf1Fun pdf1,
                                RaceCdf1Fun cdf1,
                                void* ctx,
                                double* logS_k);

double integrate_for_kth_winner_rowmajor_cpp(
    int k_winner_idx,
    const double* pars_rowmajor,
    const int* isok_int,
    double low,
    double upp,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    int n_lR_j,
    int n_par,
    const GslIntegrationControls& gsl_ctl,
    void* model_specific_context,
    gsl_integration_workspace* w);

double get_trunc_normaliser_rowmajor_cpp(const double* pars_rowmajor,
                                         const int* isok_int,
                                         RacePdf1Fun pdf1,
                                         RaceCdf1Fun cdf1,
                                         double LT,
                                         double UT,
                                         int n_lR,
                                         int n_par,
                                         const GslIntegrationControls& gsl_ctl,
                                         void* model_specific_context,
                                         GslWorkspacePtr& workspace);

#endif // EMC2_RACE_INTEGRANDS_H

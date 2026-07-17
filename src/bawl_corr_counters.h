#ifndef EMC2_BAWL_CORR_COUNTERS_H
#define EMC2_BAWL_CORR_COUNTERS_H

// Test/benchmark observability for the correlated BAwL likelihood routes
// (bawl_corr_exact_kernel_plan.md).  Counters are enabled only when the
// environment variable EMC2_BAWLCORR_COUNTERS is set to a non-empty value
// other than "0"; the flag is read once per likelihood invocation, and the
// counters are deliberately plain (non-atomic) since the particle likelihood
// runs the correlated route single-threaded.

#include <cstdlib>
#include <cstring>

struct BAwLCorrCounters {
  // Per-trial route selection (unique trials, per particle evaluation).
  long long ordinary_zero_rho_trials = 0;
  long long ordinary_single_loaded_trials = 0;
  long long exact_pair_trials = 0;
  long long exact_pair_pair_winner_trials = 0;
  long long exact_pair_independent_winner_trials = 0;
  long long exact_pair_point_start_trials = 0;
  long long numeric_pair_trials = 0;
  long long gh_no_clock_trials = 0;
  long long gh_generic_clock_trials = 0;

  // Loaded (nonzero-loading) dimension of the positivity denominator.
  long long loaded_dimension_0 = 0;
  long long loaded_dimension_1 = 0;
  long long loaded_dimension_2 = 0;
  long long loaded_dimension_3plus = 0;

  // Work volume.
  long long prepared_rows = 0;
  long long fused_node_evaluations = 0;
  long long bvn_corner_evaluations = 0;

  // Adaptive-GH behaviour.
  long long analytic_center_eligible_trials = 0;
  long long analytic_center_success_trials = 0;
  long long survivor_scan_trials = 0;
  long long scan_refinement_trials = 0;

  // Legacy-route work volume (baseline observability): per-trial node
  // evaluations in the scan/fine numerator passes and the quadrature
  // denominator sweep, plus how many trials required that sweep at all.
  long long scan_node_evaluations = 0;
  long long fine_node_evaluations = 0;
  long long den_quadrature_trials = 0;
  long long den_quadrature_node_evaluations = 0;

  void reset() { *this = BAwLCorrCounters(); }
};

inline BAwLCorrCounters& bawl_corr_counters() {
  static BAwLCorrCounters counters;
  return counters;
}

inline bool bawl_corr_counters_enabled() {
  const char* raw = std::getenv("EMC2_BAWLCORR_COUNTERS");
  return raw != nullptr && *raw != '\0' && std::strcmp(raw, "0") != 0;
}

// Cached copy of the enabled flag, set once per likelihood invocation, so
// header-level work counters (e.g. BVN corner evaluations) do not re-read
// the environment on hot paths.
inline bool& bawl_corr_counters_active() {
  static bool active = false;
  return active;
}

#endif

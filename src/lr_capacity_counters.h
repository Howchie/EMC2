#ifndef EMC2_LR_CAPACITY_COUNTERS_H
#define EMC2_LR_CAPACITY_COUNTERS_H

// Test/benchmark observability for the LogicalRules correlated-capacity
// likelihood routes. Counters are
// enabled only when the environment variable EMC2_LRCAP_COUNTERS is set to a
// non-empty value other than "0"; the flag is read once per likelihood
// invocation, and the counters are deliberately plain (non-atomic) since the
// particle likelihood runs the LogicalRules route single-threaded.

#include <cstdlib>
#include <cstring>

struct LrCapacityCounters {
  // Per-trial route selection (unique trials, per particle evaluation).
  // Ordinary trials are single-target trials, trials without capacity
  // parameters, and AB trials sitting exactly at kappa = 1, tau = 0.
  long long ordinary_trials = 0;
  long long capacity_detection_trials = 0;
  long long capacity_choice_trials = 0;
  long long invalid_trials = 0;

  // Capacity trials evaluated with the exact single-node tau = 0 shortcut
  // (kappa != 1 mean effect only; the integrand is factor-invariant).
  long long tau_zero_trials = 0;

  // Work volume.
  long long factor_node_evaluations = 0;   // conditional integrand evaluations
  long long channel_gl_evaluations = 0;    // prepared target endpoint evals in GL sums

  void reset() { *this = LrCapacityCounters(); }
};

inline LrCapacityCounters& lr_capacity_counters() {
  static LrCapacityCounters counters;
  return counters;
}

inline bool lr_capacity_counters_enabled() {
  const char* raw = std::getenv("EMC2_LRCAP_COUNTERS");
  return raw != nullptr && *raw != '\0' && std::strcmp(raw, "0") != 0;
}

#endif

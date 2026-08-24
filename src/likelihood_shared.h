#ifndef EMC2_LIKELIHOOD_SHARED_H
#define EMC2_LIKELIHOOD_SHARED_H

// ---------------------------------------------------------------------------
// Shared per-likelihood-call state contract, extracted verbatim from
// particle_ll.cpp.
//
// This header is the single definition owner for the per-data state objects
// the race and DDM likelihood kernels operate on:
//
//   DDMEndpointCacheKey / DDMEndpointCacheKeyHash / DDMEndpointCache
//       Wiener endpoint-CDF memoisation keyed on the full normalised
//       parameter row of a single trial.
//   ModelSharedState
//       Pre-computed per-data state reused across particles: censoring /
//       truncation bounds, finite/non-finite partition, scratch buffers,
//       contaminant/guess columns.
//   RaceSharedState / DDMSharedState
//       Family aliases over ModelSharedState; both families deliberately
//       share one layout so particle_ll.cpp and likelihood_ddm.cpp see
//       identical types.
//
// Field order, default initializers, types, comments and methods below are
// part of the contract and must stay identical for every translation unit
// that includes this header.  Includes are limited to the contract headers
// that layout needs (RcppArmadillo, race_integrands.h, contaminant_mixture.h,
// standard containers/types).
// ---------------------------------------------------------------------------

#include <RcppArmadillo.h>

// Included first because composite_functions.h (via contaminant_mixture.h)
// consumes its emc2_isfinite fast-math predicates; this keeps the header
// self-contained so any including translation unit sees the same layout.
#include "utility_functions.h"

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <vector>

#include "race_integrands.h"     // RaceEndpointGroupCache (ModelSharedState member)
#include "contaminant_mixture.h" // GuessKernel (ModelSharedState member)


struct DDMEndpointCacheKey {
  double rt = 0.0;
  double v = 0.0;
  double a = 0.0;
  double sv = 0.0;
  double t0 = 0.0;
  double st0 = 0.0;
  double s = 0.0;
  double Z = 0.0;
  double sz = 0.0;
  int response = 0;

  bool operator==(const DDMEndpointCacheKey& other) const {
    return rt == other.rt && v == other.v && a == other.a &&
           sv == other.sv && t0 == other.t0 && st0 == other.st0 &&
           s == other.s && Z == other.Z && sz == other.sz &&
           response == other.response;
  }
};

struct DDMEndpointCacheKeyHash {
  static inline void mix(std::size_t& h, std::size_t x) {
    h ^= x + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
         (h << 6) + (h >> 2);
  }

  std::size_t operator()(const DDMEndpointCacheKey& key) const noexcept {
    std::size_t h = 0;
    const std::hash<double> hd;
    mix(h, hd(key.rt));
    mix(h, hd(key.v));
    mix(h, hd(key.a));
    mix(h, hd(key.sv));
    mix(h, hd(key.t0));
    mix(h, hd(key.st0));
    mix(h, hd(key.s));
    mix(h, hd(key.Z));
    mix(h, hd(key.sz));
    mix(h, std::hash<int>{}(key.response));
    return h;
  }
};

struct DDMEndpointCache {
  std::unordered_map<DDMEndpointCacheKey, double, DDMEndpointCacheKeyHash> values;

  // ParamTable columns are refilled in place for every particle.  Drop the
  // entries, but retain the map's bucket allocation for the next particle.
  void new_particle() { values.clear(); }
};


// Pre-computed per-data state for likelihood functions (Race, DDM, etc.).
// Built once outside the particle loop; reused across particles to eliminate
// per-particle R-heap allocations and repeated attribute/column reads.
struct ModelSharedState {
  bool valid = false;
  // Pre-read censoring/truncation bounds (Rcpp vectors keep memory alive)
  Rcpp::NumericVector LT_vec, UT_vec, LC_vec, UC_vec;
  // Pre-computed finite/other trial partition
  Rcpp::LogicalVector finite_mask;       // length n_trials; shared ref from dadm attr
  std::vector<int>    finite_mask_int;   // 0/1 representation (for DDM/SIMD)
  std::vector<int>    finite_unique_idx; // indices of finite unique trials
  std::vector<int>    other_unique_idx;  // indices of other (non-finite) unique trials
  std::vector<int>    active_nogo_trial_mask; // per-unique-trial 0/1 nogo-active dispatch flag
  // Pre-allocated mutable scratch buffers (length n_trials)
  std::vector<double> res_buf;   // log-density or result; NOT re-initialised between particles
  std::vector<int>    idx_win;   // data-fixed winner mask  (finite rows only)
  std::vector<int>    idx_loss;  // data-fixed loser mask   (finite rows only)
  std::vector<int>    ok_int_buf;  // per-particle validity   (re-filled each call)
  bool any_win  = false;
  bool any_loss = false;
  // pContaminant column: -2=not yet searched, -1=absent, >=0=column index
  int  pc_col   = -2;
  // pGuess column and the uniform guess kernel (window resolved once in R by
  // resolve_guess_window(); see src/contaminant_mixture.h).
  int  pg_col   = -2;
  GuessKernel guess;
  // Direct column pointers for the DDM path, which does not carry a keep_names
  // index for the trailing nuisance columns.  nullptr when the model omits them.
  const double* pc_ptr = nullptr;
  const double* pg_ptr = nullptr;
  int time_code = -1;
  int nogo_code = -1;
  std::vector<int> idx_time_only;   // time-accumulator mask (finite rows only)
  std::vector<int> n_resp;          // per-unique-trial guessable accumulator count (excl. time, nogo)
  std::vector<double> alt_res_buf;  // scratch for timed-race f_T and S_W

  // DDM-specific data
  std::vector<double> logF_LT_1, logF_LT_2, logF_UT_1, logF_UT_2;
  DDMEndpointCache ddm_endpoint_cache;
  RaceEndpointGroupCache race_endpoint_cache;
  // Pre-allocated scratch buffers for nonfinite/trunc path; avoids per-particle R heap.
  std::vector<double> lF_LC_1_buf, lF_LC_2_buf, lF_UC_1_buf, lF_UC_2_buf;
  std::vector<int>    R1_int_buf, R2_int_buf;  // constant all-1 / all-2 response vectors
  std::vector<int>    all_ones_int_buf;        // constant all-1 mask
  bool any_ok_finite = false;
  bool any_ok_nonfinite = false;
  SEXP shared_R_levels = R_NilValue;     // Cached response levels for Go/No-go logic
};

struct RaceSharedState : ModelSharedState {};
using DDMSharedState = ModelSharedState;


#endif // EMC2_LIKELIHOOD_SHARED_H

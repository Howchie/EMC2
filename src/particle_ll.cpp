#include "utility_functions.h"
#include "race_dispatch.h"
#include "model_lnr.h"
#include "model_EXG.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_SDT.h"
#include "model_MRI.h"
#include "model_SS_adapters.h"
#include "ss_raw.h"
#include "composite_functions.h"
#include "trend.h"
#include "utils.h"
#include "model_PCOUNTER.h"
#include "wald_functions.h"
#include "gsl_utils.h"
#include "ParamTable.h"
#include "col_registry.h"
#include "TrendEngine.h"
#include "transform_utils.h"
#include "gh_quad.h"
#include "bawl_corr_counters.h"
#include "lr_capacity_counters.h"
#include "logicalrules_likelihood.h"
#include "correlated_likelihood.h"
#include "marginal_likelihood.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include "bawl_geometry.h"
#include "bawl_corr_exact.h"
#include "drift_factor.h"
#include "contaminant_mixture.h"
#include <cmath>
#include <string>
#include <memory>
#include <array>
#include <algorithm>   // std::max, for the scratch-buffer widths
#include <vector>
#include <cstring>
#include <cstdlib>
#include <mutex>
#include <functional>
#include <unordered_map>
#include <cstdint>

using namespace Rcpp;

// Count accumulators in [start, start+n) that are neither the time accumulator
// nor the nogo accumulator.  Used to determine the number of guessable responses
// for the timed-race and Erlang-guess likelihood paths.
static inline int count_resp_accumulators(const int* lR_codes, int start, int n,
                                          int time_code, int nogo_code) {
  int count = 0;
  for (int k = 0; k < n; ++k) {
    const int code = lR_codes[start + k];
    if (code != time_code && code != nogo_code) count++;
  }
  return count;
}

NumericMatrix get_pars_c_wrapper_oo_core(NumericMatrix particle_matrix,
                                         DataFrame data,
                                         NumericVector constants,
                                         List designs,
                                         List bounds,
                                         List transforms,
                                         List pretransforms,
                                         Rcpp::Nullable<Rcpp::List> trend,
                                         bool return_kernel_matrix,
                                         bool return_all_pars,
                                         IntegerVector kernel_output_codes);

NumericVector get_pars_c_batch_wrapper_oo_core(NumericMatrix particle_matrix,
                                               DataFrame data,
                                               NumericVector constants,
                                               List designs,
                                               List bounds,
                                               List transforms,
                                               List pretransforms,
                                               Rcpp::Nullable<Rcpp::List> trend,
                                               bool return_kernel_matrix,
                                               bool return_all_pars,
                                               IntegerVector kernel_output_codes);

NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
                         Rcpp::Nullable<Rcpp::List> marginalise);


static void update_pt_only(ParamTable& param_table,
                           const Rcpp::List& designs,
                           TrendRuntime* trend_runtime,
                           const std::vector<TransformSpec>& full_specs,
                           const Rcpp::LogicalVector* invariant_design_mask = nullptr,
                           const std::unordered_set<std::string>* invariant_param_names = nullptr,
                           const Rcpp::LogicalVector* planned_map_next = nullptr,
                           const std::vector<TransformSpec>* planned_postmap = nullptr) {
  // Planned lane: with no trend runtime the premap / pretransform /
  // posttransform blocks below are all inert, and the design mask and postmap
  // spec list depend only on names (designs, split set, invariant set) — never
  // on parameter values.  The caller resolves them once per likelihood call,
  // which removes the per-particle set-of-strings rebuild and the
  // STRSXP -> std::string conversion filter_specs_by_param_set does per spec.
  if (!trend_runtime && planned_map_next && planned_postmap) {
    param_table.map_from_designs(designs, *planned_map_next);
    c_do_transform_pt(param_table, *planned_postmap);
    return;
  }

  if (trend_runtime) trend_runtime->reset_all_kernels();

  const int n_designs = designs.size();
  LogicalVector map_next(n_designs, false);
  std::unordered_set<std::string> transform_next;
  std::unordered_set<std::string> empty_set;

  const auto& premap_set = trend_runtime ? trend_runtime->premap_trend_params() : empty_set;
  const auto& pretransform_set = trend_runtime ? trend_runtime->pretransform_trend_params() : empty_set;

  if (trend_runtime && trend_runtime->has_premap()) {
    map_next = trend_runtime->premap_design_mask(designs);
    param_table.map_from_designs(designs, map_next);

    const auto& specs_premap = trend_runtime->premap_specs();
    if (!specs_premap.empty()) {
      c_do_transform_pt(param_table, specs_premap);
    }

    for (std::size_t i = 0; i < trend_runtime->premap_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->premap_ops[i], param_table);
    }
  }

  for (int i = 0; i < n_designs; ++i) {
    bool is_premap = (trend_runtime && trend_runtime->has_premap()) ? map_next[i] : false;
    const bool is_invariant = (invariant_design_mask && i < invariant_design_mask->size())
      ? static_cast<bool>((*invariant_design_mask)[i]) : false;
    map_next[i] = (!is_premap) && (!is_invariant);
  }
  param_table.map_from_designs(designs, map_next);

  if (trend_runtime && trend_runtime->has_pretransform()) {
    const auto& specs_pretransform = trend_runtime->pretransform_specs();
    if (!specs_pretransform.empty()) {
      c_do_transform_pt(param_table, specs_pretransform);
    }

    for (std::size_t i = 0; i < trend_runtime->pretransform_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->pretransform_ops[i], param_table);
    }
  }

  const auto split_set = param_table.split_transform_params();
  transform_next = param_names_excluding(param_table, { &premap_set, &pretransform_set });
  for (const auto& nm : split_set) transform_next.erase(nm);
  if (invariant_param_names && !invariant_param_names->empty()) {
    for (const auto& nm : *invariant_param_names) transform_next.erase(nm);
  }
  auto postmap_specs = filter_specs_by_param_set(param_table, full_specs, transform_next);
  c_do_transform_pt(param_table, postmap_specs);

  if (trend_runtime && trend_runtime->has_posttransform()) {
    for (std::size_t i = 0; i < trend_runtime->posttransform_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->posttransform_ops[i], param_table);
    }
  }
}

static bool ddm_data_all_finite_untruncated(const Rcpp::DataFrame& data,
                                            const int n_trials) {
  Rcpp::NumericVector rts = data["rt"];
  Rcpp::IntegerVector R = data["R"];

  for (int i = 0; i < n_trials; ++i) {
    if (!R_FINITE(rts[i]) || Rcpp::NumericVector::is_na(rts[i])) return false;
    if (R[i] == NA_INTEGER) return false;
  }

  if (data.containsElementNamed("LT")) {
    Rcpp::NumericVector LT = data["LT"];
    for (int i = 0; i < n_trials; ++i) {
      if (LT[i] != 0.0) return false;
    }
  }

  if (data.containsElementNamed("UT")) {
    Rcpp::NumericVector UT = data["UT"];
    for (int i = 0; i < n_trials; ++i) {
      if (R_FINITE(UT[i])) return false;
    }
  }

  return true;
}

static bool race_data_all_finite_untruncated(const Rcpp::DataFrame& data,
                                             const int n_trials,
                                             const int n_lR) {
  if (n_trials <= 0 || n_lR <= 0 || (n_trials % n_lR) != 0) return false;

  Rcpp::NumericVector rts = data["rt"];
  Rcpp::IntegerVector R = data["R"];

  for (int start = 0; start < n_trials; start += n_lR) {
    if (!R_FINITE(rts[start]) || Rcpp::NumericVector::is_na(rts[start]) || rts[start] <= 0.0) {
      return false;
    }
    if (R[start] == NA_INTEGER) return false;
  }

  if (data.containsElementNamed("LT")) {
    Rcpp::NumericVector LT = data["LT"];
    for (int start = 0; start < n_trials; start += n_lR) {
      if (LT[start] != 0.0) return false;
    }
  }

  if (data.containsElementNamed("UT")) {
    Rcpp::NumericVector UT = data["UT"];
    for (int start = 0; start < n_trials; start += n_lR) {
      if (R_FINITE(UT[start])) return false;
    }
  }

  return true;
}


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

static inline void race_endpoint_hash_mix(std::size_t& h, std::size_t x) {
  h ^= x + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
       (h << 6) + (h >> 2);
}

static inline std::size_t race_endpoint_double_hash(double x) {
  std::uint64_t bits = 0;
  if (x != 0.0) std::memcpy(&bits, &x, sizeof(bits));
  return std::hash<std::uint64_t>{}(bits);
}

static inline bool race_endpoint_block_equal(
    const double* const* cols, int n_lR, int n_par,
    const int* isok, int row_a, int row_b, int skip_a, int skip_b) {
  for (int k = 0; k < n_lR; ++k) {
    const int a = row_a + k;
    const int b = row_b + k;
    if ((isok[a] != 0) != (isok[b] != 0)) return false;
    if (!isok[a]) continue;
    for (int c = 0; c < n_par; ++c) {
      if (c == skip_a || c == skip_b) continue;
      if (cols[c][a] != cols[c][b]) return false;
    }
  }
  return true;
}

static inline void race_endpoint_prepare_groups(
    RaceEndpointGroupCache& cache, const double* const* cols,
    int n_unique_trials, int n_lR, int n_par,
    const int* include_mask, const int* isok, int skip_a, int skip_b) {
  if (cache.prepared) return;
  cache.group_id.assign(static_cast<size_t>(n_unique_trials), -1);

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!include_mask[j]) continue;
    ++cache.n_included;
    const int start = j * n_lR;
    std::size_t h = 0;
    for (int k = 0; k < n_lR; ++k) {
      const int row = start + k;
      race_endpoint_hash_mix(h, std::hash<int>{}(isok[row] ? 1 : 0));
      if (!isok[row]) continue;
      for (int c = 0; c < n_par; ++c) {
        if (c == skip_a || c == skip_b) continue;
        race_endpoint_hash_mix(h, race_endpoint_double_hash(cols[c][row]));
      }
    }

    int group = -1;
    auto& candidates = cache.hash_groups[h];
    for (int candidate : candidates) {
      if (race_endpoint_block_equal(cols, n_lR, n_par, isok, start,
                                    cache.representative[static_cast<size_t>(candidate)] * n_lR,
                                    skip_a, skip_b)) {
        group = candidate;
        break;
      }
    }
    if (group < 0) {
      group = static_cast<int>(cache.representative.size());
      cache.representative.push_back(j);
      candidates.push_back(group);
    }
    cache.group_id[static_cast<size_t>(j)] = group;
  }

  const int n_groups = static_cast<int>(cache.representative.size());
  if (n_groups == cache.n_included) {
    cache.prepared = true;
    return;
  }
  cache.compact_cols.assign(static_cast<size_t>(n_groups) * n_lR * n_par, 0.0);
  cache.compact_isok.assign(static_cast<size_t>(n_groups) * n_lR, 0);
  for (int g = 0; g < n_groups; ++g) {
    const int source = cache.representative[static_cast<size_t>(g)] * n_lR;
    const int target = g * n_lR;
    for (int k = 0; k < n_lR; ++k) {
      cache.compact_isok[static_cast<size_t>(target + k)] = isok[source + k];
      for (int c = 0; c < n_par; ++c) {
        cache.compact_cols[static_cast<size_t>(c) * n_groups * n_lR + target + k] =
          cols[c][source + k];
      }
    }
  }
  cache.compact_col_ptrs.resize(static_cast<size_t>(n_par));
  for (int c = 0; c < n_par; ++c)
    cache.compact_col_ptrs[static_cast<size_t>(c)] =
      cache.compact_cols.data() + static_cast<size_t>(c) * n_groups * n_lR;
  cache.compact_mask.assign(static_cast<size_t>(n_groups), 1);
  cache.prepared = true;
}

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

double c_log_likelihood_race(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    const int n_trials,
    LogicalVector winner,
    Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector isok,
    int n_lR,
    void* model_context_for_funcs,
    bool all_finite_trials,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,
    RaceSharedState* shared = nullptr,
    NumericVector* trial_ll_out = nullptr,
    bool apply_truncation_correction = true);



// The Wiener DDM's entries in a DDMAdapter.  Thin shims that drop the model
// context the analytic model has no use for; they inline away, so the default
// adapter runs exactly the code this kernel used to call directly.
inline void ddm_wien_d_raw(const double* rts, const int* Rs,
                           const double* const* cols, int n_rows,
                           const int* mask, const int* is_ok,
                           double* out, double floor_,
                           ContextForDDMModels* /*ctx*/) {
  d_DDM_Wien_raw(rts, Rs, cols, n_rows, mask, is_ok, out, floor_);
}
inline void ddm_wien_p_raw(const double* rts, const int* Rs,
                           const double* const* cols, int n_rows,
                           const int* mask, const int* is_ok,
                           double* out, double floor_,
                           ContextForDDMModels* /*ctx*/) {
  p_DDM_Wien_raw(rts, Rs, cols, n_rows, mask, is_ok, out, floor_);
}

// The adapter every existing DDM caller gets when it passes none.
inline const DDMAdapter& ddm_wien_adapter() {
  static const DDMAdapter a = [] {
    DDMAdapter x;
    x.d_raw = &ddm_wien_d_raw;
    x.p_raw = &ddm_wien_p_raw;
    x.col_spec = emc2col::ddm::spec();
    x.endpoint_cdf_cache = true;
    return x;
  }();
  return a;
}

static inline bool ddm_wiener_endpoint_key(const double* rts,
                                           const int* Rs,
                                           const double* const* cols,
                                           int i,
                                           DDMEndpointCacheKey& key) {
  const double v = cols[emc2col::ddm::v][i];
  const double a = cols[emc2col::ddm::a][i];
  const double sv = cols[emc2col::ddm::sv][i];
  const double t0 = cols[emc2col::ddm::t0][i];
  const double st0 = cols[emc2col::ddm::st0][i];
  const double s = cols[emc2col::ddm::s][i];
  const double Z = cols[emc2col::ddm::Z][i];
  const double SZ = cols[emc2col::ddm::SZ][i];
  if (!R_FINITE(rts[i]) || !R_FINITE(v) || !R_FINITE(a) ||
      !R_FINITE(sv) || !R_FINITE(t0) || !R_FINITE(st0) ||
      !R_FINITE(s) || !R_FINITE(Z) || !R_FINITE(SZ) || s == 0.0) {
    return false;
  }
  const double sz = (Z < (1.0 - Z)) ? 2.0 * SZ * Z
                                    : 2.0 * SZ * (1.0 - Z);
  if (!R_FINITE(sz)) return false;
  key.rt = rts[i];
  key.v = v / s;
  key.a = a / s;
  key.sv = sv / s;
  key.t0 = t0;       // DDM integration depends on the effective time.
  key.st0 = st0;
  key.s = s;         // Retain the raw scale in the cache key.
  key.Z = Z;
  key.sz = sz;
  key.response = Rs[i];
  return R_FINITE(key.v) && R_FINITE(key.a) && R_FINITE(key.sv);
}

static inline void ddm_wien_p_raw_cached(
    const double* rts, const int* Rs, const double* const* cols, int n_rows,
    const int* mask, const int* is_ok, double* out, double min_ll,
    const DDMAdapter& ker, ContextForDDMModels* kctx,
    DDMEndpointCache& cache) {
  const double* sv = cols[emc2col::ddm::sv];
  const double* st0 = cols[emc2col::ddm::st0];
  const double* SZ = cols[emc2col::ddm::SZ];
  bool any_numeric = false;
  for (int i = 0; i < n_rows; ++i) {
    if (mask[i] && is_ok[i] &&
        (sv[i] != 0.0 || SZ[i] != 0.0 || st0[i] != 0.0)) {
      any_numeric = true;
      break;
    }
  }
  if (!any_numeric) {
    ker.p_raw(rts, Rs, cols, n_rows, mask, is_ok, out, min_ll, kctx);
    return;
  }
  std::vector<int> direct_mask(static_cast<size_t>(n_rows), 0);
  bool have_direct = false;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!is_ok[i] || (sv[i] == 0.0 && SZ[i] == 0.0 && st0[i] == 0.0)) {
      direct_mask[static_cast<size_t>(i)] = 1;
      have_direct = true;
      continue;
    }

    DDMEndpointCacheKey key;
    if (!ddm_wiener_endpoint_key(rts, Rs, cols, i, key)) {
      direct_mask[static_cast<size_t>(i)] = 1;
      have_direct = true;
      continue;
    }

    const auto found = cache.values.find(key);
    if (found != cache.values.end()) {
      out[i] = found->second;
      continue;
    }

    std::array<const double*, emc2col::ddm::N_REQ> row_cols{};
    for (int c = 0; c < emc2col::ddm::N_REQ; ++c)
      row_cols[static_cast<size_t>(c)] = cols[c] + i;
    const int row_mask = 1;
    double value = min_ll;
    ker.p_raw(rts + i, Rs + i, row_cols.data(), 1, &row_mask,
              is_ok + i, &value, min_ll, kctx);
    cache.values.emplace(key, value);
    out[i] = value;
  }

  if (have_direct) {
    ker.p_raw(rts, Rs, cols, n_rows, direct_mask.data(), is_ok, out,
              min_ll, kctx);
  }
}

// Collapsing-bound variants of the bounded OU, selected by
// BOU(boundary_collapse=).  As for ROU, the c_name suffix carries the FORM only;
// the shape parameters (aInf/tau/pw) are ordinary optional columns the design
// system estimates like any other.  The "B" prefix on each suffix is what keeps
// these from colliding with the ROU_B* names, which are matched the same way.
inline int bou_bnd_kind_from_type(const std::string& type_std) {
  if (type_std.find("BOU_BWEIB") != std::string::npos)
    return fpe::FPE_BND_WEIBULL;
  if (type_std.find("BOU_BEXP") != std::string::npos)
    return fpe::FPE_BND_EXPONENTIAL;
  if (type_std.find("BOU_BLIN_MULT") != std::string::npos)
    return fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
  if (type_std.find("BOU_BLIN_ADD") != std::string::npos)
    return fpe::FPE_BND_LINEAR_ADDITIVE;
  return fpe::FPE_BND_FIXED;
}

// Raw-buffer variant for DDM to skip materialization and allocations.
// Handles truncation and censoring with high numerical stability.
//
// `ker` supplies the only two model-specific operations (see DDMAdapter in
// utils.h); it defaults to the Wiener pair, so every existing call site is
// unchanged in behaviour.
double c_log_likelihood_DDM_pt(const double* const* cols,
                               const double* rt_ptr,
                               const int* R_ptr,
                               const int n_trials,
                               const int* expand_ptr,
                               const int n_out,
                               double min_ll,
                               const int* is_ok,
                               bool gng,
                               bool all_finite_untruncated,
                               ModelSharedState* shared,
                               Rcpp::NumericVector* trial_ll_out = nullptr,
                               const DDMAdapter* ker_in = nullptr) {

  // `cols` and n_trials are the same for every call below, so bind them once and
  // let the call sites read as they did before.
  const DDMAdapter& ker = (ker_in != nullptr) ? *ker_in : ddm_wien_adapter();
  ContextForDDMModels* kctx = const_cast<ContextForDDMModels*>(&ker.ctx);
  if (shared != nullptr) shared->ddm_endpoint_cache.new_particle();
  auto d_raw = [&](const double* rts, const int* Rs, const int* mask,
                   const int* ok_, double* out, double floor_) {
    ker.d_raw(rts, Rs, cols, n_trials, mask, ok_, out, floor_, kctx);
  };
  auto p_raw = [&](const double* rts, const int* Rs, const int* mask,
                   const int* ok_, double* out, double floor_) {
    if (ker.endpoint_cdf_cache && shared != nullptr) {
      ddm_wien_p_raw_cached(rts, Rs, cols, n_trials, mask, ok_, out,
                            floor_, ker, kctx, shared->ddm_endpoint_cache);
    } else {
      ker.p_raw(rts, Rs, cols, n_trials, mask, ok_, out, floor_, kctx);
    }
  };

  // Contaminant mixture (pContaminant omission + pGuess uniform outlier).
  // Both are no-ops at 0, which is their default, so a DDM design that does not
  // name them is bit-for-bit unchanged.  See src/contaminant_mixture.h.
  const double* pc_ptr = shared->pc_ptr;
  const double* pg_ptr = shared->pg_ptr;
  const GuessKernel& gk = shared->guess;
  const bool use_mix = (pc_ptr != nullptr) || (pg_ptr != nullptr && gk.active());
  auto apply_mix = [&](double ll, int i) -> double {
    const double pC = (pc_ptr != nullptr) ? pc_ptr[i] : 0.0;
    const double pG = (pg_ptr != nullptr) ? pg_ptr[i] : 0.0;
    if (pC == 0.0 && pG == 0.0) return ll;
    return mix_contaminants_rt(ll, pC, pG, gk, rt_ptr[i], R_ptr[i] != NA_INTEGER);
  };

  // 1. Fast Path: All RTs finite, no truncation, no censoring
  if (all_finite_untruncated) {
    if (shared->all_ones_int_buf.size() != static_cast<size_t>(n_trials)) {
      shared->all_ones_int_buf.assign(n_trials, 1);
    }
    d_raw(rt_ptr, R_ptr,
          shared->all_ones_int_buf.data(), is_ok, shared->res_buf.data(), min_ll);
    if (use_mix) {
      for (int i = 0; i < n_trials; ++i)
        if (is_ok[i]) shared->res_buf[i] = apply_mix(shared->res_buf[i], i);
    }

    const double* lls_ptr = shared->res_buf.data();
    double total_ll = 0.0;
    if (expand_ptr == nullptr) {
      for (int i = 0; i < n_trials; ++i) {
        double v = lls_ptr[i];
        if (!R_FINITE(v) || v < min_ll) v = min_ll;
        if (trial_ll_out != nullptr) (*trial_ll_out)[i] = v;
        total_ll += v;
      }
    } else if (trial_ll_out != nullptr) {
      for (int i = 0; i < n_out; ++i) {
        const int idx = expand_ptr[i] - 1;
        double v = lls_ptr[idx];
        if (!R_FINITE(v) || v < min_ll) v = min_ll;
        (*trial_ll_out)[i] = v;
        total_ll += v;
      }
    } else {
      #pragma omp simd reduction(+:total_ll)
      for (int i = 0; i < n_out; ++i) {
        const int idx = expand_ptr[i] - 1;
        double v = lls_ptr[idx];
        if (!R_FINITE(v) || v < min_ll) v = min_ll;
        total_ll += v;
      }
    }
    return total_ll;
  }

  // 2. Comprehensive Path: Handles truncation, censoring, and non-finite RTs
  const double* LT = shared->LT_vec.begin();
  const double* UT = shared->UT_vec.begin();
  
  std::fill(shared->res_buf.begin(), shared->res_buf.end(), min_ll);
  shared->any_ok_finite = false;
  shared->any_ok_nonfinite = false;
  
  for (int i = 0; i < n_trials; ++i) {
    if (R_FINITE(rt_ptr[i])) {
      shared->finite_mask_int[static_cast<size_t>(i)] = 1;
      if (is_ok[i]) shared->any_ok_finite = true;
    } else {
      shared->finite_mask_int[static_cast<size_t>(i)] = 0;
      if (is_ok[i]) shared->any_ok_nonfinite = true;
    }
  }

  // Helper for stable log-difference of CDFs
  auto clamp_log_cdf = [&](double logF) -> double {
    if (!R_FINITE(logF)) return R_NegInf;
    if (logF > 0.0) return 0.0;
    return logF;
  };
  auto log_sum_cdf = [&](double logF_a, double logF_b) -> double {
    return clamp_log_cdf(log_sum_exp(clamp_log_cdf(logF_a), clamp_log_cdf(logF_b)));
  };
  auto log_interval_mass = [&](double logF_low, double logF_high) -> double {
    const double lo = clamp_log_cdf(logF_low);
    const double hi = clamp_log_cdf(logF_high);
    if (hi == R_NegInf) return R_NegInf;
    if (lo == R_NegInf) return hi;
    // Allow slight numerical noise where lo > hi due to integration error
    if (hi <= lo + 1e-14) return R_NegInf; 
    return log_diff_exp(hi, lo);
  };

  // Pre-calculate truncation normalizers (logZ) for all trials
  std::vector<double> logZ(static_cast<size_t>(n_trials), 0.0);
  if (!gng) {
    // These buffers are also consumed by the censored (non-finite RT) branch
    // below. Even when truncation is inactive, keep them sized with stable
    // defaults corresponding to LT=0 and UT=Inf:
    //   logF(LT=0) = -Inf, logF(UT=Inf) = 0
    if (shared->logF_LT_1.size() != static_cast<size_t>(n_trials)) shared->logF_LT_1.resize(n_trials);
    if (shared->logF_LT_2.size() != static_cast<size_t>(n_trials)) shared->logF_LT_2.resize(n_trials);
    if (shared->logF_UT_1.size() != static_cast<size_t>(n_trials)) shared->logF_UT_1.resize(n_trials);
    if (shared->logF_UT_2.size() != static_cast<size_t>(n_trials)) shared->logF_UT_2.resize(n_trials);
    std::fill(shared->logF_LT_1.begin(), shared->logF_LT_1.end(), R_NegInf);
    std::fill(shared->logF_LT_2.begin(), shared->logF_LT_2.end(), R_NegInf);
    std::fill(shared->logF_UT_1.begin(), shared->logF_UT_1.end(), 0.0);
    std::fill(shared->logF_UT_2.begin(), shared->logF_UT_2.end(), 0.0);

    bool any_LT = false, any_UT_finite = false;
    for (int i = 0; i < n_trials; ++i) {
      if (is_ok[i]) {
        if (LT[i] != 0.0) any_LT = true;
        if (R_FINITE(UT[i])) any_UT_finite = true;
      }
      if (any_LT && any_UT_finite) break;
    }
    const bool any_trunc = any_LT || any_UT_finite;
    if (any_trunc) {
      // Use pre-allocated buffers from shared state (avoids per-particle R heap)
      const int* R1_ptr = shared->R1_int_buf.empty()
          ? nullptr : shared->R1_int_buf.data();
      const int* R2_ptr = shared->R2_int_buf.empty()
          ? nullptr : shared->R2_int_buf.data();
      const int* ones_ptr = shared->all_ones_int_buf.empty()
          ? nullptr : shared->all_ones_int_buf.data();
      // Fallback local buffers for the c_log_likelihood_DDM (non-pt) compatibility path
      Rcpp::IntegerVector R1_local, R2_local;
      std::vector<int> all_ones_local;
      if (!R1_ptr || !R2_ptr || !ones_ptr) {
        R1_local = Rcpp::IntegerVector(n_trials, 1);
        R2_local = Rcpp::IntegerVector(n_trials, 2);
        all_ones_local.assign(n_trials, 1);
        R1_ptr  = R1_local.begin();
        R2_ptr  = R2_local.begin();
        ones_ptr = all_ones_local.data();
      }
      // Only call p_DDM at LT when some trial has LT > 0 — avoids p_DDM(0,...) calls.
      // Only call p_DDM at UT when some trial has finite UT — avoids p_DDM(Inf,...) calls
      // which trigger expensive/degenerate Wiener evaluations at t=Inf.
      if (any_LT) {
        std::vector<int> LT_mask(n_trials, 0);
        for (int i = 0; i < n_trials; ++i) {
          if (is_ok[i] && LT[i] > 0.0) LT_mask[i] = 1;
        }
        p_raw(shared->LT_vec.begin(), R1_ptr,
              LT_mask.data(), is_ok, shared->logF_LT_1.data(), R_NegInf);
        p_raw(shared->LT_vec.begin(), R2_ptr,
              LT_mask.data(), is_ok, shared->logF_LT_2.data(), R_NegInf);
      }
      if (any_UT_finite) {
        std::vector<int> UT_mask(n_trials, 0);
        for (int i = 0; i < n_trials; ++i) {
          if (is_ok[i] && R_FINITE(UT[i])) UT_mask[i] = 1;
        }
        p_raw(shared->UT_vec.begin(), R1_ptr,
              UT_mask.data(), is_ok, shared->logF_UT_1.data(), R_NegInf);
        p_raw(shared->UT_vec.begin(), R2_ptr,
              UT_mask.data(), is_ok, shared->logF_UT_2.data(), R_NegInf);
      }

      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i]) continue;
        double hi_val = R_FINITE(UT[i]) ? log_sum_cdf(shared->logF_UT_1[i], shared->logF_UT_2[i]) : 0.0;
        double lo_val = (LT[i] == 0.0) ? R_NegInf : log_sum_cdf(shared->logF_LT_1[i], shared->logF_LT_2[i]);
        logZ[i] = log_interval_mass(lo_val, hi_val);
      }
    }
  }

  // Calculate finite RT densities
  if (shared->any_ok_finite) {
  d_raw(rt_ptr, R_ptr,
        shared->finite_mask_int.data(), is_ok, shared->res_buf.data(), min_ll);
  if (!gng) {

      for (int i = 0; i < n_trials; ++i) {
        if (is_ok[i] && shared->finite_mask_int[i]) {
          if (R_FINITE(logZ[i])) {
            if (shared->res_buf[i] > min_ll) shared->res_buf[i] -= logZ[i];
          }
          else shared->res_buf[i] = min_ll;
        }
      }
    }
  }

  // Calculate non-finite RT masses (censored)
    if (shared->any_ok_nonfinite) {
    if (gng) {
      // Go/No-go DDM (untruncated)
      if (shared->all_ones_int_buf.size() != static_cast<size_t>(n_trials)) {
        shared->all_ones_int_buf.assign(n_trials, 1);
      }
      if (shared->R1_int_buf.size() != static_cast<size_t>(n_trials)) {
        shared->R1_int_buf.assign(n_trials, 1);
      }
      int* R_go_ptr = shared->R1_int_buf.data();
      SEXP lev_sexp = shared->shared_R_levels;
      if (lev_sexp != R_NilValue) {
        Rcpp::CharacterVector levs(lev_sexp);
        int nogo_idx = -1;
        for (int j = 0; j < levs.size(); ++j) {
          if (Rcpp::as<std::string>(levs[j]) == "nogo") { nogo_idx = j + 1; break; }
        }
        if (nogo_idx != -1) {
          int go_idx = 1;
          for (int j = 0; j < levs.size(); ++j) {
            if (j + 1 != nogo_idx) { go_idx = j + 1; break; }
          }
          std::fill(shared->R1_int_buf.begin(), shared->R1_int_buf.end(), go_idx);
          R_go_ptr = shared->R1_int_buf.data();
        }
      }
      Rcpp::NumericVector logcdf_U(n_trials), logcdf_L(n_trials);
      p_raw(shared->UC_vec.begin(), R_go_ptr,
            shared->all_ones_int_buf.data(), is_ok, logcdf_U.begin(), R_NegInf);
      p_raw(shared->LC_vec.begin(), R_go_ptr,
            shared->all_ones_int_buf.data(), is_ok, logcdf_L.begin(), R_NegInf);
      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i] || shared->finite_mask_int[i]) continue;
        if (rt_ptr[i] == R_PosInf) shared->res_buf[i] = log1m_exp(logcdf_U[i]);
        else if (rt_ptr[i] == R_NegInf) shared->res_buf[i] = logcdf_L[i];
        else shared->res_buf[i] = log_sum_exp(logcdf_L[i], log1m_exp(logcdf_U[i]));
      }
    } else {
      // Standard DDM censored trials — use pre-allocated shared buffers to avoid
      // per-particle R heap allocation (critical for LT path where n_nonfinite > 0).
      const bool have_shared_bufs = (shared != nullptr &&
                                     !shared->R1_int_buf.empty() &&
                                     (int)shared->R1_int_buf.size() >= n_trials);
      // Fallback allocations (only used when called outside calc_ll_oo, e.g. standalone)
      Rcpp::IntegerVector R1_fb, R2_fb;
      Rcpp::NumericVector lF_LC_1_fb, lF_LC_2_fb, lF_UC_1_fb, lF_UC_2_fb;
      std::vector<int> all_ones_fb;
      int*    R1_ptr;    int*    R2_ptr;    int*    ones_ptr;
      double* lF_LC_1;  double* lF_LC_2;
      double* lF_UC_1;  double* lF_UC_2;
      if (have_shared_bufs) {
        R1_ptr   = shared->R1_int_buf.data();
        R2_ptr   = shared->R2_int_buf.data();
        ones_ptr = shared->all_ones_int_buf.data();
        lF_LC_1  = shared->lF_LC_1_buf.data();
        lF_LC_2  = shared->lF_LC_2_buf.data();
        lF_UC_1  = shared->lF_UC_1_buf.data();
        lF_UC_2  = shared->lF_UC_2_buf.data();
      } else {
        R1_fb = Rcpp::IntegerVector(n_trials, 1);
        R2_fb = Rcpp::IntegerVector(n_trials, 2);
        lF_LC_1_fb = Rcpp::NumericVector(n_trials);
        lF_LC_2_fb = Rcpp::NumericVector(n_trials);
        lF_UC_1_fb = Rcpp::NumericVector(n_trials);
        lF_UC_2_fb = Rcpp::NumericVector(n_trials);
        all_ones_fb.assign(n_trials, 1);
        R1_ptr   = R1_fb.begin();
        R2_ptr   = R2_fb.begin();
        ones_ptr = all_ones_fb.data();
        lF_LC_1  = lF_LC_1_fb.begin();
        lF_LC_2  = lF_LC_2_fb.begin();
        lF_UC_1  = lF_UC_1_fb.begin();
        lF_UC_2  = lF_UC_2_fb.begin();
      }
      
      std::vector<int> nonfinite_mask(n_trials, 0);
      for (int i = 0; i < n_trials; ++i) {
        if (!shared->finite_mask_int[i] && is_ok[i]) nonfinite_mask[i] = 1;
      }
      
      p_raw(shared->LC_vec.begin(), R1_ptr,
            nonfinite_mask.data(), is_ok, lF_LC_1, R_NegInf);
      p_raw(shared->LC_vec.begin(), R2_ptr,
            nonfinite_mask.data(), is_ok, lF_LC_2, R_NegInf);
      p_raw(shared->UC_vec.begin(), R1_ptr,
            nonfinite_mask.data(), is_ok, lF_UC_1, R_NegInf);
      p_raw(shared->UC_vec.begin(), R2_ptr,
            nonfinite_mask.data(), is_ok, lF_UC_2, R_NegInf);

      for (int i = 0; i < n_trials; ++i) {
        if (!is_ok[i] || shared->finite_mask_int[i]) continue;
        const bool r_known = (R_ptr[i] != NA_INTEGER);
        const int r_idx = r_known ? R_ptr[i] : 0;
        double mass = R_NegInf;
        if (rt_ptr[i] == R_NegInf) {
          if (!r_known) mass = log_sum_cdf(log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]), log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]));
          else if (r_idx == 1) mass = log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]);
          else mass = log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]);
        } else if (rt_ptr[i] == R_PosInf) {
          // P(UC < RT < UT).
          // When UT is finite, logF_UT_r = log F_r(UT) (computed above), and the
          // interval formula is correct: (F_1(UT)-F_1(UC)) + (F_2(UT)-F_2(UC)).
          // When UT=Inf, logF_UT_r was initialised to 0 (representing F_r(Inf)=1),
          // but the DDM CDF is DEFECTIVE: F_1(Inf)+F_2(Inf)=1, not F_r(Inf)=1
          // individually.  Using logF_UT_r=0 gives (1-F_1)+(1-F_2)=2-F_total >= 1,
          // which log_sum_cdf clamps to 0 — every censored trial contributes mass=1
          // regardless of parameters.  The correct mass is 1-F_total(UC).
          if (!r_known) {
            if (R_FINITE(UT[i])) {
              mass = log_sum_cdf(log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i]),
                                 log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]));
            } else {
              // P(RT > UC) = 1 - F_1(UC) - F_2(UC)
              mass = log_diff_exp(0.0, log_sum_exp(lF_UC_1[i], lF_UC_2[i]));
            }
          } else if (r_idx == 1) {
            mass = log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i]);
          } else {
            mass = log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]);
          }
        } else { // NA
          double m1 = (!r_known) ? log_sum_cdf(log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]), log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]))
                                 : (r_idx == 1) ? log_interval_mass(shared->logF_LT_1[i], lF_LC_1[i]) : log_interval_mass(shared->logF_LT_2[i], lF_LC_2[i]);
          // Same defective-CDF fix as the R_PosInf branch above.
          double m2;
          if (!r_known) {
            m2 = R_FINITE(UT[i])
              ? log_sum_cdf(log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i]),
                            log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]))
              : log_diff_exp(0.0, log_sum_exp(lF_UC_1[i], lF_UC_2[i]));
          } else {
            m2 = (r_idx == 1) ? log_interval_mass(lF_UC_1[i], shared->logF_UT_1[i])
                               : log_interval_mass(lF_UC_2[i], shared->logF_UT_2[i]);
          }
          mass = log_sum_exp(m1, m2);
        }
        if (R_FINITE(logZ[i])) shared->res_buf[i] = mass - logZ[i];
        else shared->res_buf[i] = min_ll;
      }
    }
  }

  // Contaminant mixture, applied after truncation renormalisation and after the
  // censored-interval branches -- exactly where the race kernels apply it, so
  // pGuess is the guess proportion among *retained* trials.  A guess can never
  // be censored or truncated away (the window is [max(LT,LC), min(UC,UT)]), so
  // the censored branches need no guess term of their own.
  if (use_mix) {
    for (int i = 0; i < n_trials; ++i)
      if (is_ok[i]) shared->res_buf[i] = apply_mix(shared->res_buf[i], i);
  }

  // 3. Accumulate results
  const double* res_ptr = shared->res_buf.data();
  double total_ll = 0.0;
  if (expand_ptr == nullptr) {
    for (int i = 0; i < n_trials; ++i) {
      double v = res_ptr[i];
      if (!R_FINITE(v) || v < min_ll) v = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = v;
      total_ll += v;
    }
  } else if (trial_ll_out != nullptr) {
    for (int i = 0; i < n_out; ++i) {
      double v = res_ptr[expand_ptr[i] - 1];
      if (!R_FINITE(v) || v < min_ll) v = min_ll;
      (*trial_ll_out)[i] = v;
      total_ll += v;
    }
  } else {
    #pragma omp simd reduction(+:total_ll)
    for (int i = 0; i < n_out; ++i) {
      double v = res_ptr[expand_ptr[i] - 1];
      if (!R_FINITE(v) || v < min_ll) v = min_ll;
      total_ll += v;
    }
  }
  return total_ll;
}

// Compatibility wrapper for old calc_ll path
double c_log_likelihood_DDM(Rcpp::NumericMatrix pars, Rcpp::DataFrame data,
                            const int n_trials, Rcpp::IntegerVector expand,
                            double min_ll, Rcpp::LogicalVector is_ok, bool gng,
                            bool all_finite_untruncated = false,
                            Rcpp::NumericVector* trial_ll_out = nullptr) {
  ModelSharedState shared;
  shared.LT_vec = get_col_with_default(data, "LT", 0.0);
  shared.UT_vec = get_col_with_default(data, "UT", R_PosInf);
  shared.LC_vec = get_col_with_default(data, "LC", 0.0);
  shared.UC_vec = get_col_with_default(data, "UC", R_PosInf);
  shared.finite_mask_int.assign(static_cast<size_t>(n_trials), 0);
  shared.res_buf.assign(static_cast<size_t>(n_trials), min_ll);
  shared.ok_int_buf.assign(static_cast<size_t>(n_trials), 0);
  shared.all_ones_int_buf.assign(static_cast<size_t>(n_trials), 1);
  shared.R1_int_buf.assign(static_cast<size_t>(n_trials), 1);
  shared.R2_int_buf.assign(static_cast<size_t>(n_trials), 2);
  shared.lF_LC_1_buf.assign(static_cast<size_t>(n_trials), R_NegInf);
  shared.lF_LC_2_buf.assign(static_cast<size_t>(n_trials), R_NegInf);
  shared.lF_UC_1_buf.assign(static_cast<size_t>(n_trials), R_NegInf);
  shared.lF_UC_2_buf.assign(static_cast<size_t>(n_trials), R_NegInf);
  Rcpp::IntegerVector R_col = data["R"];
  shared.shared_R_levels = R_col.attr("levels");
  shared.valid = true;

  Rcpp::NumericVector rts = data["rt"];
  Rcpp::IntegerVector R = data["R"];
  const int* expand_ptr = (expand.length() > 0) ? expand.begin() : nullptr;
  const int n_out_val = (expand.length() > 0) ? expand.length() : n_trials;

  std::vector<int> ok_int(n_trials);
  for(int i=0; i<n_trials; ++i) ok_int[i] = is_ok[i] ? 1 : 0;

  // The materialized fallback may retain a non-canonical storage order. Resolve
  // the canonical DDM pointer array from column names once, just as the raw
  // ParamTable path does.
  const emc2col::ColSpec ddm_spec = emc2col::ddm::spec();
  Rcpp::CharacterVector par_names = Rcpp::colnames(pars);
  if (static_cast<int>(par_names.size()) != pars.ncol()) {
    Rcpp::stop("c_log_likelihood_DDM: materialized parameter matrix has no column names.");
  }
  const double* ddm_cols[emc2col::ddm::N_REQ];
  for (int k = 0; k < ddm_spec.n_required; ++k) {
    int col = -1;
    for (int j = 0; j < pars.ncol(); ++j) {
      if (Rcpp::as<std::string>(par_names[j]) == ddm_spec.names[k]) {
        col = j;
        break;
      }
    }
    if (col < 0) {
      Rcpp::stop("c_log_likelihood_DDM: missing required parameter column '%s'.",
                 ddm_spec.names[k]);
    }
    ddm_cols[k] = pars.begin() + static_cast<size_t>(col) * n_trials;
  }
  // Trailing nuisance columns, resolved by name from the same matrix.
  for (int j = 0; j < pars.ncol(); ++j) {
    const std::string nm = Rcpp::as<std::string>(par_names[j]);
    if (nm == "pContaminant") shared.pc_ptr = pars.begin() + static_cast<size_t>(j) * n_trials;
    else if (nm == "pGuess")  shared.pg_ptr = pars.begin() + static_cast<size_t>(j) * n_trials;
  }
  shared.guess = resolve_guess_kernel(data);
  shared.guess.pg_col = (shared.pg_ptr != nullptr) ? 0 : -1;

  return c_log_likelihood_DDM_pt(ddm_cols, rts.begin(), R.begin(),
                                 n_trials, expand_ptr, n_out_val, min_ll,
                                 ok_int.data(), gng, all_finite_untruncated,
                                 &shared, trial_ll_out);}

// ---------------------------------------------------------------------------
// Shared helpers used by calc_ll, calc_ll_oo, and calc_ll_oo_pw
// ---------------------------------------------------------------------------

// Read the emc2_all_finite_trials data attribute; fall back to a direct scan
// when the attribute is absent (e.g. direct C++ callers that bypass R-side cache).
static inline bool read_all_finite_trials_attr(const Rcpp::DataFrame& data,
                                               int n_trials, int n_lR) {
  if (data.hasAttribute("emc2_all_finite_trials")) {
    Rcpp::LogicalVector v = data.attr("emc2_all_finite_trials");
    if (v.size() == 1 && v[0] != NA_LOGICAL) return (bool)v[0];
  }
  if (n_trials > 0 && n_lR > 0 && (n_trials % n_lR) == 0) {
    return race_data_all_finite_untruncated(data, n_trials, n_lR);
  }
  return false;
}


double c_log_likelihood_huvsd(Rcpp::NumericMatrix pars, Rcpp::DataFrame data,

                              const int n_trials, Rcpp::IntegerVector expand,
                              double min_ll, Rcpp::LogicalVector is_ok,
                              Rcpp::NumericVector* trial_ll_out = nullptr) {
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



// ---------------------------------------------------------------------------
// Per-call ParamTable machinery shared by calc_ll_oo and calc_ll_oo_pw:
// pretransformed particle matrix, template table, transform specs, optional
// trend runtime, invariant-parameter bookkeeping, and bound specs.
// prepare(i) is the per-particle prologue (refill from particle row i,
// design mapping + transforms, then bounds). i == 0 runs on the template as
// built (must fully compute every parameter, including transforms of
// constants); i > 0 may skip invariant designs/transforms.
struct PtMapper {
  NumericMatrix particle_matrix_pt;
  ParamTable table;
  std::vector<TransformSpec> transform_specs;
  std::unique_ptr<TrendPlan> trend_plan;
  std::unique_ptr<TrendRuntime> trend_runtime;
  Rcpp::CharacterVector keep_names;
  std::vector<int> pm_col_to_base_idx;
  Rcpp::LogicalVector invariant_design_mask;
  std::unordered_set<std::string> invariant_param_names;
  bool use_invariants = false;
  std::vector<int> invariant_base_idx_vec;
  std::vector<BoundSpec> bound_specs;
  Rcpp::List designs;
  Rcpp::List bounds;
  Rcpp::NumericMatrix minmax;
  Rcpp::CharacterVector mm_names;

  // Reusable materialization: one NumericMatrix per likelihood call, refilled
  // per particle (replaces a per-particle materialize_by_param_names R-heap
  // allocation on the LogicalRules / mixed-race / pw hot paths).
  Rcpp::NumericMatrix mat_buf;
  std::vector<int> mat_base_idx;
  bool mat_ready = false;

  // Once-per-likelihood-call plans for the three per-particle steps whose cost
  // otherwise scales with parameter-table *width* even when the extra columns
  // are held constant (mapping/transform selection, base refill, bounds).
  bool plan_ready = false;
  Rcpp::LogicalVector plan_map_next;              // designs to remap for i > 0
  std::vector<TransformSpec> plan_postmap;        // transforms to apply for i > 0
  std::vector<int> plan_zero_base_idx;            // base cols to clear per particle
  std::vector<std::pair<int,int>> plan_fill_pairs;// (particle col, base col)
  std::vector<BoundSpec> bound_specs_variant;     // bounds that can change per particle
  Rcpp::LogicalVector bound_seed;                 // invariant bounds' verdict
  bool bound_plan_ready = false;

  // Build the i > 0 plans.  Only valid without a trend runtime, which is also
  // the only case in which use_invariants is ever set.
  void build_plan() {
    const int n_designs = designs.size();
    plan_map_next = Rcpp::LogicalVector(n_designs, true);
    for (int i = 0; i < n_designs; ++i) {
      if (use_invariants && i < invariant_design_mask.size() &&
          static_cast<bool>(invariant_design_mask[i])) {
        plan_map_next[i] = false;
      }
    }

    std::unordered_set<std::string> transform_next = param_names_excluding(table, {});
    for (const auto& nm : table.split_transform_params()) transform_next.erase(nm);
    if (use_invariants) {
      for (const auto& nm : invariant_param_names) transform_next.erase(nm);
    }
    plan_postmap = filter_specs_by_param_set(table, transform_specs, transform_next);

    // Refill plan: invariant base columns are neither zeroed nor refilled, so
    // their natural-scale values survive from the template particle.
    const int n_base = table.base.ncol();
    std::vector<char> is_inv_base(static_cast<size_t>(n_base), 0);
    if (use_invariants) {
      for (int bidx : invariant_base_idx_vec) {
        if (bidx >= 0 && bidx < n_base) is_inv_base[static_cast<size_t>(bidx)] = 1;
      }
    }
    plan_zero_base_idx.clear();
    for (int j = 0; j < n_base; ++j) {
      if (!is_inv_base[static_cast<size_t>(j)]) plan_zero_base_idx.push_back(j);
    }
    plan_fill_pairs.clear();
    for (std::size_t j = 0; j < pm_col_to_base_idx.size(); ++j) {
      const int bidx = pm_col_to_base_idx[j];
      if (bidx < 0 || is_inv_base[static_cast<size_t>(bidx)]) continue;
      plan_fill_pairs.emplace_back(static_cast<int>(j), bidx);
    }

    plan_ready = true;
  }

  // Split the bound specs into the ones that can change from particle to
  // particle and the ones that cannot, and cache the latter's verdict.
  void build_bound_plan() {
    bound_specs_variant.clear();
    std::vector<BoundSpec> invariant_specs;
    for (const auto& bs : bound_specs) {
      bool inv = false;
      if (use_invariants && bs.col_idx >= 0 && bs.col_idx < table.base_names.size()) {
        inv = invariant_param_names.count(
                Rcpp::as<std::string>(table.base_names[bs.col_idx])) > 0;
      }
      if (inv) invariant_specs.push_back(bs);
      else bound_specs_variant.push_back(bs);
    }
    bound_seed = invariant_specs.empty()
      ? Rcpp::LogicalVector(table.n_trials, true)
      : c_do_bound_pt(table, invariant_specs);
    bound_plan_ready = true;
  }

  Rcpp::LogicalVector prepare(int i) {
    const bool planned = !trend_runtime && plan_ready;
    if (i > 0) {
      if (planned) {
        table.fill_from_particle_row_planned(particle_matrix_pt, i,
                                             plan_zero_base_idx, plan_fill_pairs);
      } else {
        table.fill_from_particle_row(particle_matrix_pt, i,
                                     pm_col_to_base_idx,
                                     invariant_base_idx_vec);
      }
    }
    const bool skip_inv = (i > 0) && use_invariants;
    update_pt_only(table, designs, trend_runtime ? trend_runtime.get() : nullptr,
                   transform_specs,
                   skip_inv ? &invariant_design_mask : nullptr,
                   skip_inv ? &invariant_param_names : nullptr,
                   (planned && i > 0) ? &plan_map_next : nullptr,
                   (planned && i > 0) ? &plan_postmap : nullptr);
    if (i == 0) {
      bound_specs = make_bound_specs_pt(minmax, mm_names, table, bounds);
      Rcpp::LogicalVector ok = c_do_bound_pt(table, bound_specs);
      if (!trend_runtime) {
        if (!plan_ready) build_plan();
        if (!bound_plan_ready) build_bound_plan();
      }
      return ok;
    }
    if (bound_plan_ready) {
      return c_do_bound_pt_from(table, bound_specs_variant, bound_seed);
    }
    return c_do_bound_pt(table, bound_specs);
  }

  Rcpp::NumericMatrix materialize_reusable() {
    if (!mat_ready) {
      const int k = keep_names.size();
      mat_base_idx.resize(k);
      for (int j = 0; j < k; ++j) {
        // base_index_for throws for unknown names, exactly like
        // materialize_by_param_names did on this path.
        mat_base_idx[j] = table.base_index_for(Rcpp::as<std::string>(keep_names[j]));
      }
      mat_buf = Rcpp::NumericMatrix(table.n_trials, k);
      Rcpp::colnames(mat_buf) = keep_names;
      mat_ready = true;
    }
    table.materialize_into(mat_buf, mat_base_idx);
    return mat_buf;
  }
};

static PtMapper make_pt_mapper(NumericMatrix particle_matrix, DataFrame data,
                               NumericVector constants, List designs, List bounds,
                               List transforms, List pretransforms,
                               Rcpp::Nullable<Rcpp::List> trend,
                               CharacterVector p_types,
                               int n_trials, int n_particles) {
  PtMapper m;
  m.designs = designs;
  m.bounds = bounds;
  m.minmax = Rcpp::as<Rcpp::NumericMatrix>(bounds["minmax"]);
  m.mm_names = colnames(m.minmax);
  m.keep_names = p_types;

  std::vector<TransformSpec> pre_specs = make_transform_specs(particle_matrix, pretransforms);
  m.particle_matrix_pt = c_do_transform(particle_matrix, pre_specs);

  const bool has_constants = !(constants.size() == 1 &&
                               Rcpp::NumericVector::is_na(constants[0]));
  if (has_constants) {
    m.particle_matrix_pt = add_constants_columns(m.particle_matrix_pt, constants);
  }

  NumericVector p_vector = m.particle_matrix_pt(0, Rcpp::_);
  p_vector.attr("names") = colnames(m.particle_matrix_pt);
  m.table = ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials, transforms);
  m.transform_specs = make_transform_specs_for_paramtable(m.table, transforms);

  Rcpp::CharacterVector pm_names = colnames(m.particle_matrix_pt);
  m.pm_col_to_base_idx.assign(pm_names.size(), -1);
  for (int j = 0; j < pm_names.size(); ++j) {
    std::string nm = Rcpp::as<std::string>(pm_names[j]);
    auto it = m.table.name_to_base_idx.find(nm);
    if (it != m.table.name_to_base_idx.end()) {
      m.pm_col_to_base_idx[j] = it->second;
    }
  }

  if (!trend.isNull()) {
    m.trend_plan.reset(new TrendPlan(trend, data));
    m.trend_runtime.reset(new TrendRuntime(*m.trend_plan));
    m.trend_runtime->bind_all_ops_to_paramtable(m.table);
    m.trend_runtime->init_cached_specs(m.table, m.transform_specs);
  }

  // Invariant-parameter optimization: designs whose coefficients are not
  // sampled map to the same natural-scale column for every particle, so
  // particles after the first can skip re-mapping/re-transforming them.
  if (n_particles > 1 && !m.trend_runtime) {
    Rcpp::CharacterVector p_names = colnames(particle_matrix);
    std::unordered_set<std::string> sampled_coef_names;
    sampled_coef_names.reserve(p_names.size());
    for (int j = 0; j < p_names.size(); ++j) {
      sampled_coef_names.insert(Rcpp::as<std::string>(p_names[j]));
    }

    m.invariant_design_mask = Rcpp::LogicalVector(designs.size(), false);
    Rcpp::CharacterVector design_names = designs.names();
    for (int i = 0; i < designs.size(); ++i) {
      Rcpp::NumericMatrix d = designs[i];
      Rcpp::CharacterVector dcols = colnames(d);
      bool invariant = true;
      for (int c = 0; c < dcols.size(); ++c) {
        if (sampled_coef_names.find(Rcpp::as<std::string>(dcols[c])) != sampled_coef_names.end()) {
          invariant = false;
          break;
        }
      }
      m.invariant_design_mask[i] = invariant;
      if (invariant) {
        std::string pnm;
        if (design_names.size() == designs.size()) pnm = Rcpp::as<std::string>(design_names[i]);
        else if (i < m.keep_names.size()) pnm = Rcpp::as<std::string>(m.keep_names[i]);
        if (!pnm.empty()) m.invariant_param_names.insert(pnm);
      }
    }
    m.use_invariants = !m.invariant_param_names.empty();
  }

  // Precompute base column indices for invariant parameters so fill_from_particle_row
  // can preserve their natural-scale values across reset_base_to_zero.
  if (m.use_invariants) {
    for (const auto& nm : m.invariant_param_names) {
      auto it = m.table.name_to_base_idx.find(nm);
      if (it != m.table.name_to_base_idx.end())
        m.invariant_base_idx_vec.push_back(it->second);
    }
  }
  return m;
}

// DDM per-data shared state + canonical base-column pointers for the raw
// c_log_likelihood_DDM_pt kernel. Returns false (raw path unusable) when any
// canonical DDM parameter is missing from the table.
static bool init_ddm_shared_state(DataFrame data, int n_trials,
                                  const ParamTable& table,
                                  ModelSharedState& shared,
                                  std::vector<const double*>& cols,
                                  const emc2col::ColSpec* spec_in = nullptr,
                                  const Rcpp::CharacterVector* keep_names_in = nullptr) {
  shared.LT_vec = get_col_with_default(data, "LT", 0.0);
  shared.UT_vec = get_col_with_default(data, "UT", R_PosInf);
  shared.LC_vec = get_col_with_default(data, "LC", 0.0);
  shared.UC_vec = get_col_with_default(data, "UC", R_PosInf);
  shared.finite_mask_int.resize(n_trials);
  shared.res_buf.resize(n_trials);
  shared.ok_int_buf.resize(n_trials);
  // Pre-allocate DDM scratch buffers (avoids per-particle R heap in trunc/cens paths)
  shared.lF_LC_1_buf.resize(n_trials);
  shared.lF_LC_2_buf.resize(n_trials);
  shared.lF_UC_1_buf.resize(n_trials);
  shared.lF_UC_2_buf.resize(n_trials);
  shared.R1_int_buf.assign(n_trials, 1);
  shared.R2_int_buf.assign(n_trials, 2);
  shared.all_ones_int_buf.assign(n_trials, 1);
  Rcpp::IntegerVector R_col = data["R"];
  shared.shared_R_levels = R_col.attr("levels");
  shared.valid = true;

  // Contaminant nuisance columns.  Both are trailing p_types, so they sit past
  // the spec's canonical prefix and are resolved by name, never positionally.
  // The ParamTable base storage is fixed for this likelihood call (only its
  // values are refilled per particle), so these addresses stay valid.
  shared.pc_ptr = nullptr;
  shared.pg_ptr = nullptr;
  {
    auto itc = table.name_to_base_idx.find("pContaminant");
    if (itc != table.name_to_base_idx.end())
      shared.pc_ptr = table.base.begin() + static_cast<size_t>(itc->second) * n_trials;
    auto itg = table.name_to_base_idx.find("pGuess");
    if (itg != table.name_to_base_idx.end())
      shared.pg_ptr = table.base.begin() + static_cast<size_t>(itg->second) * n_trials;
  }
  shared.guess = resolve_guess_kernel(data);
  shared.guess.pg_col = (shared.pg_ptr != nullptr) ? 0 : -1;

  bool raw_ready = true;
  // The column set is the model's, not the family's: a two-boundary model with
  // extra parameters (bounded OU's leak) declares its own spec.
  const emc2col::ColSpec ddm_spec =
    (spec_in != nullptr) ? *spec_in : emc2col::ddm::spec();

  // Resolve the OPTIONAL trailing columns too, not just the required prefix.
  // A model variant may declare columns past N_REQ -- the bounded OU's
  // aInf/tau/pw for a collapsing bound -- and its kernel indexes them
  // positionally, so a `cols` sized to n_required alone is an out-of-bounds
  // read, not a missing feature.  This mirrors what the race path already does
  // (see the keep_names walk in calc_ll_oo); the required prefix has already
  // been checked by validate_col_prefix, so walking keep_names positionally
  // agrees with the spec on the first n_required entries by construction.
  if (keep_names_in != nullptr) {
    const int n_kn = keep_names_in->size();
    cols.assign(std::max(n_kn, ddm_spec.n_required), nullptr);
    for (int j = 0; j < n_kn; ++j) {
      std::string nm = Rcpp::as<std::string>((*keep_names_in)[j]);
      auto it = table.name_to_base_idx.find(nm);
      if (it != table.name_to_base_idx.end()) {
        cols[j] = table.base.begin() + static_cast<size_t>(it->second) * n_trials;
      } else if (j < ddm_spec.n_required) {
        raw_ready = false;
      }
    }
    for (int j = 0; j < ddm_spec.n_required; ++j)
      if (cols[j] == nullptr) raw_ready = false;
    return raw_ready;
  }

  for (int j = 0; j < ddm_spec.n_required; ++j) {
    auto it = table.name_to_base_idx.find(ddm_spec.names[j]);
    int idx = (it != table.name_to_base_idx.end()) ? it->second : -1;
    if (idx < 0) {
      raw_ready = false;
    } else {
      // ParamTable base storage is fixed for this likelihood call; only its
      // values are refilled per particle, so these addresses are reusable.
      cols.push_back(table.base.begin() + static_cast<size_t>(idx) * n_trials);
    }
  }
  return raw_ready;
}

// Data-fixed shared state for the mixed (non-raw-fast) race path: partition
// attrs, censoring/truncation bounds, winner/loser masks, scratch buffers.
// Leaves shared.valid=false when the partition attributes are absent.
static RaceSharedState build_race_shared_state(DataFrame data, int n_trials, int n_lR,
                                               const LogicalVector& winner,
                                               const Rcpp::IntegerVector& lR_code_vec,
                                               int time_code, int nogo_code,
                                               bool has_RACE_col,
                                               const Rcpp::IntegerVector& RACE_nacc,
                                               const Rcpp::LogicalVector& RACE_mask) {
  RaceSharedState race_shared;
  const int n_unique = n_trials / n_lR;
  const bool has_part =
    data.hasAttribute("finite_rt_mask") &&
    data.hasAttribute("finite_rt_unique_trial_indices") &&
    data.hasAttribute("other_unique_trial_indices") &&
    data.hasAttribute("active_nogo_trial_mask");
  if (!has_part) return race_shared;

  race_shared.finite_mask = data.attr("finite_rt_mask");
  Rcpp::IntegerVector fa = data.attr("finite_rt_unique_trial_indices");
  Rcpp::IntegerVector oa = data.attr("other_unique_trial_indices");
  Rcpp::LogicalVector ng = data.attr("active_nogo_trial_mask");
  race_shared.finite_unique_idx.assign(fa.begin(), fa.end());
  race_shared.other_unique_idx.assign(oa.begin(), oa.end());
  race_shared.active_nogo_trial_mask.assign(ng.begin(), ng.end());

  // Pre-read censoring/truncation bounds once
  race_shared.LT_vec = get_col_with_default(data, "LT", 0.0);
  race_shared.UT_vec = get_col_with_default(data, "UT", R_PosInf);
  race_shared.LC_vec = get_col_with_default(data, "LC", 0.0);
  race_shared.UC_vec = get_col_with_default(data, "UC", R_PosInf);

  // Pre-allocate mutable scratch buffers
  race_shared.res_buf.resize(static_cast<size_t>(n_trials)); // no init needed
  race_shared.idx_win.assign(static_cast<size_t>(n_trials), 0);
  race_shared.idx_loss.assign(static_cast<size_t>(n_trials), 0);
  race_shared.ok_int_buf.resize(static_cast<size_t>(n_trials));

  // Fill data-fixed winner/loser masks (finite trials only)
  const Rcpp::LogicalVector& fmask = race_shared.finite_mask;
  race_shared.time_code = time_code;
  race_shared.nogo_code = nogo_code;
  if (time_code != -1) {
    race_shared.idx_time_only.assign(static_cast<size_t>(n_trials), 0);
    race_shared.alt_res_buf.resize(static_cast<size_t>(n_trials));
  }
  const bool needs_n_resp_shared = (time_code != -1);
  if (needs_n_resp_shared) {
    race_shared.n_resp.assign(static_cast<size_t>(n_unique), 0);
  }
  for (int j = 0; j < n_trials; ++j) {
    if (!fmask[j]) continue;
    if (has_RACE_col && !RACE_mask[j]) continue;
    if (winner[j]) {
      race_shared.idx_win[static_cast<size_t>(j)] = 1;
      race_shared.any_win = true;
    } else if (n_lR > 1) {
      race_shared.idx_loss[static_cast<size_t>(j)] = 1;
      race_shared.any_loss = true;
    }
    if (time_code != -1 && lR_code_vec[j] == time_code) {
      race_shared.idx_time_only[static_cast<size_t>(j)] = 1;
    }
  }
  if (needs_n_resp_shared) {
    for (int j = 0; j < n_unique; ++j) {
      const int start = j * n_lR;
      const int n_lR_curr = has_RACE_col ? RACE_nacc[start] : n_lR;
      race_shared.n_resp[static_cast<size_t>(j)] = count_resp_accumulators(
          lR_code_vec.begin(), start, n_lR_curr, time_code, nogo_code);
    }
  }
  race_shared.valid = true;
  return race_shared;
}

// [[Rcpp::export]]
NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                         Rcpp::Nullable<Rcpp::List> marginalise = R_NilValue) {
  // -------------------------------------------------------------------------
  // Optional coherent marginalization of a shared subject-level parameter
  // (Stage 1: the non-decision time t0 in race/GNG models). When `marginalise`
  // is supplied we integrate ONE shared t0 out of the COMPLETE subject
  // likelihood via Gauss-Legendre quadrature on the log-t0 (sampled) axis:
  //     Lbar(theta) = int p(t0|eta) L(t0, theta) dt0,   log t0 ~ N(mu, sigma)
  // The SAME t0 is inserted into every trial at each node (coherent, not the
  // trialwise prod-of-integrals). Each node overwrites the sampled t0 column
  // and re-calls THIS function with marginalise off, so the existing kernel is
  // reused verbatim for every race model (no per-model code). The integration
  // interval is clipped to where the likelihood is live: below the t0 lower
  // bound c_do_bound floors the whole subject to min_ll; above min(rt) a
  // response is infeasible. When `marginalise` is null this is a single branch
  // test and the byte-identical existing path runs (zero cost to other models).
  if (marginalise.isNotNull()) {
    return calc_ll_oo_marginal(particle_matrix, data, constants, designs, type, bounds,
                               transforms, pretransforms, p_types, min_ll, trend,
                               Rcpp::List(marginalise));
  }
  const int n_particles = particle_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  LogicalVector is_ok(n_trials);
  NumericMatrix pars;

  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs;
  CharacterVector p_names = colnames(particle_matrix);
  std::string type_std = Rcpp::as<std::string>(Rcpp::wrap(type));

  // BOU is the DDM with leak: a different kernel behind the SAME two-boundary
  // likelihood, so it takes the DDM path and swaps only the two primitives (see
  // DDMAdapter).  `type` here is the model's c_name, so it arrives as "BOU".
  const bool is_bou_type = type_std.find("BOU") != std::string::npos;
  const bool is_ddm_type = is_bou_type ||
                           type_std.find("DDM") != std::string::npos;
  const bool is_mri_type = (type == "MRI" || type == "MRI_AR1");
  const bool is_ss_type = is_stop_signal_type(type_std);
  const bool use_pt_mapping = !is_mri_type;

  NumericMatrix one_particle(1, particle_matrix.ncol());
  colnames(one_particle) = p_names;
  IntegerVector kernel_output_codes = IntegerVector::create(1);
  auto pars_for_particle_generic = [&](int i) -> NumericMatrix {
    for (int j = 0; j < particle_matrix.ncol(); ++j) {
      one_particle(0, j) = particle_matrix(i, j);
    }
    return get_pars_c_wrapper_oo_core(one_particle, data, constants, designs, bounds, transforms,
                                      pretransforms, trend, false, false, kernel_output_codes);
  };

  PtMapper pt;
  if (use_pt_mapping) {
    pt = make_pt_mapper(particle_matrix, data, constants, designs, bounds,
                        transforms, pretransforms, trend, p_types,
                        n_trials, n_particles);
  }
  ParamTable& param_table_template = pt.table;
  Rcpp::CharacterVector& keep_names = pt.keep_names;
  auto prepare_particle = [&](int i) -> Rcpp::LogicalVector { return pt.prepare(i); };

  ModelSharedState ddm_shared;
  bool ddm_raw_ready = true;
  std::vector<const double*> ddm_cols;
  // The adapter carries the model's two primitives, its column spec and its
  // solve cache.  Built once per likelihood call; the Wiener default leaves the
  // DDM on exactly the code it had before.
  DDMAdapter ddm_adapter = ddm_wien_adapter();
  if (is_bou_type) {
    ddm_adapter.d_raw = &bou::d_BOU_raw;
    ddm_adapter.p_raw = &bou::p_BOU_raw;
    ddm_adapter.col_spec = emc2col::bou::spec();
    ddm_adapter.endpoint_cdf_cache = false;
    ddm_adapter.ctx.bou_cache = std::make_shared<fpebou::SolveCache>();
    bou::bou_configure(*ddm_adapter.ctx.bou_cache);
    ddm_adapter.ctx.bnd_kind = bou_bnd_kind_from_type(type_std);
  }
  if (is_ddm_type) {
    emc2col::validate_col_prefix(keep_names, ddm_adapter.col_spec);
    ddm_raw_ready = init_ddm_shared_state(data, n_trials, param_table_template,
                                          ddm_shared, ddm_cols,
                                          &ddm_adapter.col_spec, &keep_names);
  }

  if (is_ddm_type) {
    bool gng = (type_std.find("GNG") != std::string::npos);
    IntegerVector expand = data.attr("expand");
    NumericVector rts = data["rt"];
    IntegerVector R = data["R"];
    const double* rt_ptr = rts.begin();
    const int* R_ptr = R.begin();
    const int* expand_ptr = expand.begin();
    const int n_out = expand.length();
    const bool all_finite_untruncated = ddm_data_all_finite_untruncated(data, n_trials);
    for (int i = 0; i < n_particles; ++i) {
      is_ok = prepare_particle(i);
      // Every particle has different parameters, so solves from the previous one
      // are dead weight -- and the horizon must be re-derived rather than
      // inherited.  Keys are exact, so this is about cost and staleness of the
      // horizon, not about ever returning a wrong entry.
      if (ddm_adapter.ctx.bou_cache) {
        ddm_adapter.ctx.bou_cache->new_particle();
        ddm_adapter.ctx.bou_cache->t_horizon = 0.0;
      }
      if (ddm_raw_ready) {
        for(int j = 0; j < n_trials; ++j) ddm_shared.ok_int_buf[j] = is_ok[j] ? 1 : 0;
        lls[i] = c_log_likelihood_DDM_pt(ddm_cols.data(),
                                        rt_ptr, R_ptr, n_trials, expand_ptr, n_out,
                                        min_ll, ddm_shared.ok_int_buf.data(), gng,
                                        all_finite_untruncated, &ddm_shared,
                                        nullptr, &ddm_adapter);
      } else {
        pars = param_table_template.materialize_by_param_names(keep_names);
        lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok,
                                      gng, all_finite_untruncated);
      }
    }
  } else if (type_std == "hUVSD") {
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      is_ok = prepare_particle(i);
      pars = param_table_template.materialize_by_param_names(keep_names);
      lls[i] = c_log_likelihood_huvsd(pars, data, n_trials, expand, min_ll, is_ok);
    }
  }
 else if (is_mri_type) {
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for (int i = 0; i < n_particles; ++i) {
      pars = pars_for_particle_generic(i);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      if (type == "MRI") {
        lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      } else {
        lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
    }
  } else if (is_ss_type) {
    IntegerVector expand = data.attr("expand");
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    int n_trials_ss = (n_lR > 0) ? (n_trials / n_lR) : n_trials;
    SSModelAdapter ssa = resolve_ss_adapter(type_std);
    emc2col::validate_col_prefix(keep_names,
                                 type_std.find("EXG") != std::string::npos
                                     ? emc2col::ss_texg::spec()
                                     : emc2col::ss_rdex::spec());

    // Raw path: data-fixed trial structure once, then read ParamTable columns
    // directly per particle. Falls back to the materialized NumericMatrix path
    // if the ParamTable is missing a p_types column (should not happen for
    // SSEXG/SSRDEX, but the fallback keeps the reference implementation live).
    const bool ss_is_exg = type_std.find("EXG") != std::string::npos;
    const SsRawModel& ss_raw_model = ss_is_exg ? ss_texg_raw_model()
                                               : ss_rdex_raw_model();
    SsSharedState ss_shared = build_ss_shared_state(data, n_trials_ss, n_lR);
    std::vector<int> ss_col_base_idx(keep_names.size(), -1);
    // Test/benchmark hook: force the materialized fallback so the two paths
    // can be compared from R (see test-ss-raw-path.R).
    bool ss_raw_ready = ss_shared.valid &&
                        (std::getenv("EMC2_SS_FORCE_MATERIALIZE") == nullptr);
    for (int j = 0; j < keep_names.size(); ++j) {
      auto it = param_table_template.name_to_base_idx.find(
          Rcpp::as<std::string>(keep_names[j]));
      if (it == param_table_template.name_to_base_idx.end()) {
        ss_raw_ready = false;
        break;
      }
      ss_col_base_idx[j] = it->second;
    }
    // Base column addresses are particle-invariant (prepare_particle refills
    // values in place), so the pointer array is built once.
    std::vector<const double*> ss_cols(ss_col_base_idx.size(), nullptr);
    if (ss_raw_ready) {
      for (size_t j = 0; j < ss_col_base_idx.size(); ++j) {
        ss_cols[j] = &param_table_template.base(0, ss_col_base_idx[j]);
      }
    }
    SsRawWorkspace ss_ws;

    for (int i = 0; i < n_particles; ++i) {
      is_ok = prepare_particle(i);
      is_ok = lr_all(is_ok, n_lR);
      if (ss_raw_ready) {
        lls[i] = c_log_likelihood_ss_pt(ss_cols.data(), ss_shared, ss_raw_model,
                                        is_ok, expand, min_ll, ss_ws);
      } else {
        pars = param_table_template.materialize_by_param_names(keep_names);
        lls[i] = c_log_likelihood_ss(pars, data, n_trials_ss, expand, min_ll, is_ok,
                                     ssa.go_lpdf_ptr, ssa.go_lccdf_ptr,
                                     ssa.stop_logsurv_ptr, ssa.stop_success_ptr,
                                     ssa.idx_tf, ssa.idx_gf);
      }
    }
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    const bool is_logicalrules = (type_std.find("LogicalRules") != std::string::npos);

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll_oo");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);
    if (adapter.col_spec.names != nullptr) {
      emc2col::validate_col_prefix(keep_names, adapter.col_spec);
    }
    configure_corr_drift_context(adapter, keep_names, "calc_ll_oo");
    configure_rdmswtn_corr_context(adapter, keep_names, "calc_ll_oo");
    if (is_logicalrules && adapter.ctx.rdmswtn_correlated) {
      Rcpp::stop("calc_ll_oo: correlated RDMSWTN logical-rule races are not supported.");
    }

    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    const Rcpp::IntegerVector lR_code_vec(static_cast<SEXP>(lR));
    const Rcpp::CharacterVector lR_levels = lR_code_vec.attr("levels");
    int time_code = -1;
    int nogo_code = -1;
    for (int j = 0; j < lR_levels.size(); ++j) {
      std::string lev = Rcpp::as<std::string>(lR_levels[j]);
      if (lev == "time") time_code = j + 1;
      else if (lev == "nogo") nogo_code = j + 1;
    }
    adapter.ctx.time_code = time_code;
    adapter.ctx.nogo_code = nogo_code;

    if (is_logicalrules) {
      int kappa_col = -1;
      int tau_col = -1;
      int pc_col = -1;
      int pg_col = -1;
      for (int j = 0; j < keep_names.size(); ++j) {
        const std::string nm = Rcpp::as<std::string>(keep_names[j]);
        if (nm == "kappa") kappa_col = j;
        else if (nm == "tau") tau_col = j;
        else if (nm == "pContaminant") pc_col = j;
        else if (nm == "pGuess") pg_col = j;
      }
      const GuessKernel lr_guess = resolve_guess_kernel(data);
      const bool capacity = kappa_col >= 0 || tau_col >= 0;
      LogicalRulesSharedState logicalrules_shared =
        build_logicalrules_shared_state(data, n_trials, n_lR, capacity);
      if (capacity && lr_capacity_counters_enabled())
        lr_capacity_counters().reset();

      // Column pointers into ParamTable base in p_types order; base column
      // addresses are particle-invariant, so this replaces the per-particle
      // materialization copy. Padded with nullptr slots for optional columns.
      const int n_par_lr = keep_names.size();
      std::vector<const double*> lr_cols(std::max(n_par_lr, 16), nullptr);
      for (int j = 0; j < n_par_lr; ++j) {
        lr_cols[j] = &param_table_template.base(
            0, param_table_template.base_index_for(Rcpp::as<std::string>(keep_names[j])));
      }
      for (int i = 0; i < n_particles; ++i) {
        is_ok = prepare_particle(i);
        // LogicalRules evaluates validity at the stimulus/rule level.  Do
        // not collapse the row-wise flags across the fixed-role scaffold:
        // an inactive detector row is deliberately allowed to be irrelevant
        // for an A-only or B-only detection trial.
        lls[i] = c_log_likelihood_logicalrules(lr_cols.data(), n_par_lr, expand, min_ll, is_ok, n_lR,
                                               &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                               adapter.model_dfun_raw, adapter.model_pfun_raw,
                                               logicalrules_shared, nullptr,
                                               kappa_col, tau_col, pc_col,
                                               pg_col, &lr_guess);
      }
      return lls;
    }
    const bool all_finite_trials = read_all_finite_trials_attr(data, n_trials, n_lR);

    bool has_RACE_col_fp = data.containsElementNamed("RACE");
    bool has_RACE_attrs_fp = false;
    Rcpp::IntegerVector RACE_fp;
    Rcpp::LogicalVector RACE_mask_fp;
    if (has_RACE_col_fp &&
        data.hasAttribute("RACE_nacc_by_row") &&
        data.hasAttribute("RACE_mask")) {
      RACE_fp = data.attr("RACE_nacc_by_row");
      RACE_mask_fp = data.attr("RACE_mask");
      has_RACE_attrs_fp = (RACE_fp.size() == n_trials && RACE_mask_fp.size() == n_trials);
    }

    // --- Fast path: all-finite RTs, no truncation ---
    // Pre-allocates buffers once outside the particle loop; each particle writes
    // log-density/log-survivor directly into raw double arrays, bypassing all
    // Rcpp vector allocations inside c_log_likelihood_race.
    // Variable-accumulator RACE designs are supported when cache attrs are
    // present: inactive accumulator rows are masked out and not included in the
    // unique-trial sum.
    // pContaminant is handled inline: log(1-pC) shift for finite RTs (no-op when pC=0).
    const bool use_raw_fast_path =
      all_finite_trials &&  // already scanned above; avoids a redundant O(n) pass
      !adapter.ctx.corr_drift_active && !adapter.ctx.rdmswtn_correlated &&
      (!has_RACE_col_fp || has_RACE_attrs_fp) &&
      adapter.model_pfun_raw != nullptr &&
      adapter.model_dfun_raw != nullptr;

    const int n_unique_fp = n_trials / n_lR;
    std::vector<double> res_buf;
    std::vector<double> ll_uniq_buf;
    std::vector<int> winner_int_buf;
    std::vector<int> loser_int_buf;
    std::vector<int> isok_int_fp;
    NumericVector rts_fp_hold;
    const double* rt_ptr = nullptr;
    const int* expand_ptr = expand.begin();
    const int n_exp = expand.length();

    // Column pointers into ParamTable base, in p_types order (the layout the
    // raw kernels expect per src/col_registry.h). Base column addresses are
    // particle-invariant — prepare_particle refills values in place — so the
    // pointer array is built once, replacing the old per-particle staging copy.
    std::vector<const double*> race_cols;     // keep_names.size() pointers
    int fast_pc_col = -1;                     // keep_names position of pContaminant
    int fast_pg_col = -1;                     // keep_names position of pGuess
    GuessKernel fast_guess;                   // uniform guess kernel (data-fixed)

    std::vector<int> time_win_int_buf;
    std::vector<double> alt_res_buf_fp;
    std::vector<int> n_resp_fp;

    if (use_raw_fast_path) {
      res_buf.resize(n_trials);
      ll_uniq_buf.resize(n_unique_fp);
      winner_int_buf.resize(n_trials);
      loser_int_buf.resize(n_trials);
      isok_int_fp.resize(n_trials);
      // Allocate auxiliary buffers only where actually needed.
      // alt_res_buf_fp: needed for time-race (f_T, S_W) and global_kill (log S_K at winner).
      // n_resp_fp:      needed for time-race only.
      // time_win_int_buf: needed for time-race only.
      if (time_code != -1) {
        time_win_int_buf.assign(n_trials, 0);
        alt_res_buf_fp.resize(n_trials);
        n_resp_fp.assign(n_unique_fp, 0);
      } else if (adapter.ctx.is_global_kill) {
        alt_res_buf_fp.resize(n_trials);
      }
      // winner/loser masks are data-fixed; fill once
      for (int j = 0; j < n_trials; ++j) {
        const bool active = !has_RACE_col_fp || RACE_mask_fp[j];
        winner_int_buf[j] = (active && winner[j]) ? 1 : 0;
        loser_int_buf[j]  = (active && n_lR > 1 && !winner[j]) ? 1 : 0;
        if (time_code != -1 && active && lR_code_vec[j] == time_code) {
          time_win_int_buf[j] = 1;
        }
      }
      if (!n_resp_fp.empty()) {
        for (int j = 0; j < n_unique_fp; ++j) {
          const int start = j * n_lR;
          const int n_lR_curr = has_RACE_col_fp ? RACE_fp[start] : n_lR;
          n_resp_fp[j] = count_resp_accumulators(lR_code_vec.begin(), start, n_lR_curr,
                                                 time_code, nogo_code);
        }
      }
      rts_fp_hold = data["rt"];
      rt_ptr = rts_fp_hold.begin();

      // Resolve keep_names[j] -> base column pointer once; kernels index the
      // leading columns positionally, so a missing required column is fatal.
      // Pad with nullptr slots: kernels may fetch (not dereference) pointers
      // for optional trailing columns that this model variant lacks.
      const int n_kn = keep_names.size();
      race_cols.assign(std::max(n_kn, 16), nullptr);
      for (int j = 0; j < n_kn; ++j) {
        std::string nm = Rcpp::as<std::string>(keep_names[j]);
        auto it = param_table_template.name_to_base_idx.find(nm);
        if (it != param_table_template.name_to_base_idx.end()) {
          race_cols[j] = &param_table_template.base(0, it->second);
        } else if (j < adapter.col_spec.n_required) {
          Rcpp::stop("calc_ll_oo: required parameter column '%s' missing from ParamTable.",
                     nm.c_str());
        }
        if (nm == "pContaminant") fast_pc_col = j;
        else if (nm == "pGuess") fast_pg_col = j;
      }
      fast_guess = resolve_guess_kernel(data);
      fast_guess.pg_col = fast_pg_col;
    }

    // --- Build shared state for the mixed (non-raw-fast) path ---
    // Pre-compute all data-fixed structures once so c_log_likelihood_race can skip
    // per-particle R-heap allocations, column reads, and attribute lookups.
    RaceSharedState race_shared;
    if (!use_raw_fast_path && adapter.model_dfun_raw != nullptr && !all_finite_trials) {
      race_shared = build_race_shared_state(data, n_trials, n_lR, winner, lR_code_vec,
                                            time_code, nogo_code,
                                            has_RACE_col_fp, RACE_fp, RACE_mask_fp);
    }

    // Correlated BAwL: data-fixed layout, role mapping, and direct
    // ParamTable column pointers resolved once outside the particle loop.
    // The generic-clock fallback materialises lazily via this callback.
    CorrDriftSharedState corr_drift_shared;
    if (adapter.ctx.corr_drift_active) {
      corr_drift_shared = build_corr_drift_shared_state(
          data, n_trials, n_lR, winner, expand, adapter.ctx,
          param_table_template, keep_names);
    }
    RDMSWTNCorrSharedState rdmswtn_corr_shared;
    if (adapter.ctx.rdmswtn_correlated) {
      rdmswtn_corr_shared = build_rdmswtn_corr_shared_state(
          data, n_trials, n_lR, winner, expand, adapter.ctx,
          param_table_template, keep_names);
    }
    const std::function<Rcpp::NumericMatrix()> corr_materialize =
        [&pt]() { return pt.materialize_reusable(); };

    for (int i = 0; i < n_particles; ++i) {
      is_ok = prepare_particle(i);
      is_ok = lr_all(is_ok, n_lR);
      // PDE-backed models cache one solve per distinct parameter tuple; the
      // tuples change with the particle, so drop them here.  Keys are exact, so
      // this is about bounding memory, not about correctness.
      if (adapter.ctx.fpe_cache) adapter.ctx.fpe_cache->new_particle();
      if (adapter.ctx.rlf_cache) adapter.ctx.rlf_cache->new_particle();
      if (use_raw_fast_path) {
        // Fill per-particle isok buffer
        for (int j = 0; j < n_trials; ++j) isok_int_fp[j] = is_ok[j] ? 1 : 0;

        const double* const* pars_cols = race_cols.data();
        adapter.ctx.mode_hint = 0;
        adapter.ctx.kill_active = adapter.ctx.is_local_guess ||
                                  adapter.ctx.is_global_kill ||
                                  adapter.ctx.is_local_kill ||
                                  adapter.ctx.is_local_kill_guess;
        // Set once-per-particle mode hints so raw kernels can skip per-row
        // variability checks in common zero-variability cases.
        if (type_std.find("RDMSWTN") != std::string::npos) {
          const double* sv_col = pars_cols[emc2col::rdmswtn::sv];
          bool sv_zero = true;
          bool lambda_active = false;
          for (int j = 0; j < n_trials; ++j) {
            if (!isok_int_fp[j]) continue;
            const double sv = sv_col[j];
            if (emc2_isfinite(sv) && std::fabs(sv) > 1e-10) { sv_zero = false; break; }
          }
          // For any Erlang variant, disable kill bookkeeping when all lambdas are zero
          // so kernels fall through to standard Wald without computing mixtures.
          const bool any_erlang = adapter.ctx.is_global_kill || adapter.ctx.is_local_kill ||
                                  adapter.ctx.is_local_guess || adapter.ctx.is_local_kill_guess;
          if (any_erlang) {
            const double* lambda_g_col = pars_cols[emc2col::rdmswtn::mG];
            const double* lambda_k_col = pars_cols[emc2col::rdmswtn::mK];
            for (int j = 0; j < n_trials; ++j) {
              if (!isok_int_fp[j]) continue;
              const double lg = erlang_lambda_from_mean(lambda_g_col[j], adapter.ctx.kill_shape);
              const double lk = erlang_lambda_from_mean(lambda_k_col[j], adapter.ctx.kill_shape);
              if ((emc2_isfinite(lg) && lg > 1e-12) || (emc2_isfinite(lk) && lk > 1e-12)) {
                lambda_active = true;
                break;
              }
            }
            adapter.ctx.kill_active = lambda_active;
          } else {
            // No erlang type active: disable kill so scalar kernels don't read
            // past the end of the parameter block (col 7 doesn't exist for "none").
            adapter.ctx.kill_active = false;
          }
          adapter.ctx.mode_hint = sv_zero ? 1 : 2;
        }

        const bool floor_raw_log_lik_prev = adapter.ctx.floor_raw_log_lik;
        if (time_code != -1) adapter.ctx.floor_raw_log_lik = false;

        // Log-density for winner rows
        adapter.model_dfun_raw(rt_ptr, pars_cols, n_trials,
                               winner_int_buf.data(), isok_int_fp.data(),
                               res_buf.data(), min_ll, &adapter.ctx);

        // Log-survivor for loser rows (writes into same buffer, different slots)
        if (n_lR > 1) {
          adapter.model_pfun_raw(rt_ptr, pars_cols, n_trials,
                                 loser_int_buf.data(), isok_int_fp.data(),
                                 res_buf.data(), min_ll, &adapter.ctx);
        }

        if (time_code != -1) {
          // Get f_T
          adapter.model_dfun_raw(rt_ptr, pars_cols, n_trials,
                                 time_win_int_buf.data(), isok_int_fp.data(),
                                 alt_res_buf_fp.data(), min_ll, &adapter.ctx);
          // Get S_W
          adapter.model_pfun_raw(rt_ptr, pars_cols, n_trials,
                                 winner_int_buf.data(), isok_int_fp.data(),
                                 alt_res_buf_fp.data(), min_ll, &adapter.ctx);
          adapter.ctx.floor_raw_log_lik = floor_raw_log_lik_prev;
        } else if (adapter.ctx.is_global_kill && adapter.ctx.kill_active) {
          // Pre-fill log S_K at the winner row for each trial.
          const double* mean_k_ptr = pars_cols[adapter.ctx.mean_k_index];
          for (int j = 0; j < n_trials; ++j) {
            if (!isok_int_fp[j] || !winner_int_buf[j]) continue;
            const double tt = rt_ptr[j];
            alt_res_buf_fp[j] = (tt > 0.0)
                ? erlang_log_surv(tt, erlang_lambda_from_mean(mean_k_ptr[j], adapter.ctx.kill_shape), adapter.ctx.kill_shape)
                : min_ll;
          }
        }
        if (time_code == -1) adapter.ctx.floor_raw_log_lik = floor_raw_log_lik_prev;
        // Sum n_lR rows per unique trial; apply pC correction (no-op when pC=0).
        for (int j = 0; j < n_unique_fp; ++j) {
          double s = 0.0;
          const int base = j * n_lR;
          const int n_lR_j = has_RACE_col_fp ? RACE_fp[base] : n_lR;
          for (int k = 0; k < n_lR_j; ++k) s += res_buf[base + k];

          if (adapter.ctx.is_global_kill && adapter.ctx.kill_active) {
            int idx_W = -1;
            for (int k = 0; k < n_lR_j; ++k) { if (winner_int_buf[base + k]) { idx_W = k; break; } }
            if (idx_W != -1) {
              s += alt_res_buf_fp[base + idx_W]; // log S_K at winner RT
            }
          } else if (time_code != -1 && n_resp_fp[j] > 0) {
            int idx_W = -1, idx_T = -1;
            for (int k = 0; k < n_lR_j; ++k) {
              if (winner_int_buf[base + k]) idx_W = k;
              if (time_win_int_buf[base + k]) idx_T = k;
            }
            if (idx_W != -1 && idx_T != -1) {
              double s_T = 0.0;
              for (int k = 0; k < n_lR_j; ++k) {
                if (k == idx_W) s_T += alt_res_buf_fp[base + idx_W];
                else if (k == idx_T) s_T += alt_res_buf_fp[base + idx_T];
                else s_T += res_buf[base + k];
              }
              s = log_sum_exp(s, s_T - std::log(n_resp_fp[j]));
            }
          }
          s = s < min_ll ? min_ll : s;
          // Contaminant mixture.  This path is reached only when every trial has
          // a finite rt and a known R (see .cache_ll_data_attrs), so no trial is
          // an omission and every trial is guess-eligible; mix_contaminants()
          // still owns the arithmetic so the sites cannot drift apart.
          if (fast_pc_col >= 0 || (fast_pg_col >= 0 && fast_guess.active())) {
            const double pC = (fast_pc_col >= 0) ? pars_cols[fast_pc_col][base] : 0.0;
            // A guess must not lift an invalid-parameter trial off min_ll; the
            // pC branch stays unconditional, exactly as it was before pGuess.
            const double pG = (fast_pg_col >= 0 && fast_guess.active() && isok_int_fp[base])
                                ? pars_cols[fast_pg_col][base] : 0.0;
            if (pC != 0.0 || pG != 0.0)
              s = mix_contaminants(s, pC, pG, fast_guess.log_g, false);
          }
          ll_uniq_buf[j] = s;
        }

        // Expand unique-trial LLs and accumulate (1-based expand indices)
        double total_ll = 0.0;
        #pragma omp simd reduction(+:total_ll)
        for (int ei = 0; ei < n_exp; ++ei) {
          total_ll += ll_uniq_buf[expand_ptr[ei] - 1];
        }
        lls[i] = total_ll;
      } else if (adapter.ctx.corr_drift_active) {
        lls[i] = c_log_likelihood_corr_drift(
            corr_drift_shared, data, adapter.pdf1_ptr, adapter.cdf1_ptr,
            n_trials, winner, expand, min_ll, is_ok, n_lR, &adapter.ctx,
            all_finite_trials, adapter.model_dfun_raw, adapter.model_pfun_raw,
            adapter.logS_at_t_ptr,
            race_shared.valid ? &race_shared : nullptr, nullptr,
            corr_materialize);
      } else if (adapter.ctx.rdmswtn_correlated) {
        lls[i] = c_log_likelihood_rdmswtn_correlated(
            rdmswtn_corr_shared, data, adapter.pdf1_ptr, adapter.cdf1_ptr,
            n_trials, winner, expand, min_ll, is_ok, n_lR, &adapter.ctx,
            all_finite_trials, adapter.model_dfun_raw, adapter.model_pfun_raw,
            adapter.logS_at_t_ptr,
            race_shared.valid ? &race_shared : nullptr, nullptr,
            corr_materialize);
      } else {
        pars = pt.materialize_reusable();
        lls[i] = c_log_likelihood_race(pars, data,
                                       adapter.pdf1_ptr, adapter.cdf1_ptr,
                                       n_trials,
                                       winner, expand, min_ll, is_ok, n_lR,
                                       &adapter.ctx,
                                       all_finite_trials,
                                       adapter.model_dfun_raw,
                                       adapter.model_pfun_raw,
                                       adapter.logS_at_t_ptr,
                                       race_shared.valid ? &race_shared : nullptr);
      }
    }
  }
  return lls;
}

// [[Rcpp::export]]
NumericMatrix calc_ll_oo_pw(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                            List designs, String type, List bounds, List transforms, List pretransforms,
                            CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue) {
  const int n_particles = particle_matrix.nrow();
  const int n_trials = data.nrow();
  LogicalVector is_ok(n_trials);
  NumericMatrix pars;
  std::string type_std = Rcpp::as<std::string>(Rcpp::wrap(type));

  if (type == "MRI" || type == "MRI_AR1" ||
      is_stop_signal_type(type_std) ||
      type_std.find("SOFTMAX") != std::string::npos) {
    Rcpp::stop("calc_ll_oo_pw: not implemented for model type '%s'", type_std.c_str());
  }

  // Same per-call ParamTable machinery as calc_ll_oo: one template build, then
  // a cheap refill + remap per particle (instead of the full mapping wrapper).
  PtMapper pt = make_pt_mapper(particle_matrix, data, constants, designs, bounds,
                               transforms, pretransforms, trend, p_types,
                               n_trials, n_particles);
  ParamTable& param_table_template = pt.table;
  Rcpp::CharacterVector& keep_names = pt.keep_names;

  const bool is_bou_pw = type_std.find("BOU") != std::string::npos;
  if (is_bou_pw || type_std.find("DDM") != std::string::npos) {
    bool gng = (type_std.find("GNG") != std::string::npos);
    // Same adapter arrangement as calc_ll_oo; see there.
    DDMAdapter ddm_adapter = ddm_wien_adapter();
    if (is_bou_pw) {
      ddm_adapter.d_raw = &bou::d_BOU_raw;
      ddm_adapter.p_raw = &bou::p_BOU_raw;
      ddm_adapter.col_spec = emc2col::bou::spec();
      ddm_adapter.endpoint_cdf_cache = false;
      ddm_adapter.ctx.bou_cache = std::make_shared<fpebou::SolveCache>();
      bou::bou_configure(*ddm_adapter.ctx.bou_cache);
      ddm_adapter.ctx.bnd_kind = bou_bnd_kind_from_type(type_std);
    }
    IntegerVector expand = data.attr("expand");
    const int n_out = (expand.length() > 0) ? expand.length() : n_trials;
    const bool all_finite_untruncated = ddm_data_all_finite_untruncated(data, n_trials);
    ModelSharedState ddm_shared;
    emc2col::validate_col_prefix(keep_names, ddm_adapter.col_spec);
    std::vector<const double*> ddm_cols;
    const bool ddm_raw_ready = init_ddm_shared_state(data, n_trials, param_table_template,
                                                     ddm_shared, ddm_cols,
                                                     &ddm_adapter.col_spec, &keep_names);
    NumericVector rts = data["rt"];
    IntegerVector R = data["R"];
    const int* expand_ptr = (expand.length() > 0) ? expand.begin() : nullptr;
    NumericMatrix result(n_particles, n_out);
    for (int i = 0; i < n_particles; ++i) {
      is_ok = pt.prepare(i);
      if (ddm_adapter.ctx.bou_cache) {
        ddm_adapter.ctx.bou_cache->new_particle();
        ddm_adapter.ctx.bou_cache->t_horizon = 0.0;
      }
      NumericVector row_vec(n_out);
      if (ddm_raw_ready) {
        for (int j = 0; j < n_trials; ++j) ddm_shared.ok_int_buf[j] = is_ok[j] ? 1 : 0;
        c_log_likelihood_DDM_pt(ddm_cols.data(), rts.begin(), R.begin(),
                                n_trials, expand_ptr, n_out, min_ll,
                                ddm_shared.ok_int_buf.data(), gng, all_finite_untruncated,
                                &ddm_shared, &row_vec, &ddm_adapter);
      } else {
        pars = param_table_template.materialize_by_param_names(keep_names);
        c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok,
                             gng, all_finite_untruncated, &row_vec);
      }
      result(i, _) = row_vec;
    }
    return result;
  } else if (type_std == "hUVSD") {
    IntegerVector expand = data.attr("expand");
    const int n_out = (expand.length() > 0) ? expand.length() : n_trials;
    NumericMatrix result(n_particles, n_out);
    for (int i = 0; i < n_particles; ++i) {
      is_ok = pt.prepare(i);
      pars = param_table_template.materialize_by_param_names(keep_names);
      NumericVector row_vec(n_out);
      c_log_likelihood_huvsd(pars, data, n_trials, expand, min_ll, is_ok, &row_vec);
      result(i, _) = row_vec;
    }
    return result;
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    const bool is_logicalrules = (type_std.find("LogicalRules") != std::string::npos);

    RaceModelAdapter adapter = resolve_race_model_adapter(type_std, "calc_ll_oo_pw");
    adapter.ctx.min_lik_for_pdf = std::exp(min_ll);
    if (adapter.col_spec.names != nullptr) {
      emc2col::validate_col_prefix(keep_names, adapter.col_spec);
    }
    configure_corr_drift_context(adapter, keep_names, "calc_ll_oo_pw");
    configure_rdmswtn_corr_context(adapter, keep_names, "calc_ll_oo_pw");
    if (is_logicalrules && adapter.ctx.rdmswtn_correlated) {
      Rcpp::stop("calc_ll_oo_pw: correlated RDMSWTN logical-rule races are not supported.");
    }

    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    const Rcpp::IntegerVector lR_code_vec(static_cast<SEXP>(lR));
    const Rcpp::CharacterVector lR_levels = lR_code_vec.attr("levels");
    int time_code = -1;
    int nogo_code = -1;
    for (int j = 0; j < lR_levels.size(); ++j) {
      std::string lev = Rcpp::as<std::string>(lR_levels[j]);
      if (lev == "time") time_code = j + 1;
      else if (lev == "nogo") nogo_code = j + 1;
    }
    adapter.ctx.time_code = time_code;
    adapter.ctx.nogo_code = nogo_code;

    const bool all_finite_trials = read_all_finite_trials_attr(data, n_trials, n_lR);
    const int n_out_race = (expand.length() > 0) ? expand.length() : (n_trials / n_lR);
    NumericMatrix result(n_particles, n_out_race);

    if (is_logicalrules) {
      int kappa_col = -1;
      int tau_col = -1;
      int pc_col = -1;
      int pg_col = -1;
      for (int j = 0; j < keep_names.size(); ++j) {
        const std::string nm = Rcpp::as<std::string>(keep_names[j]);
        if (nm == "kappa") kappa_col = j;
        else if (nm == "tau") tau_col = j;
        else if (nm == "pContaminant") pc_col = j;
        else if (nm == "pGuess") pg_col = j;
      }
      const GuessKernel lr_guess = resolve_guess_kernel(data);
      const bool capacity = kappa_col >= 0 || tau_col >= 0;
      LogicalRulesSharedState logicalrules_shared =
        build_logicalrules_shared_state(data, n_trials, n_lR, capacity);
      if (capacity && lr_capacity_counters_enabled())
        lr_capacity_counters().reset();
      const int n_par_lr = keep_names.size();
      std::vector<const double*> lr_cols(std::max(n_par_lr, 16), nullptr);
      for (int j = 0; j < n_par_lr; ++j) {
        lr_cols[j] = &param_table_template.base(
            0, param_table_template.base_index_for(Rcpp::as<std::string>(keep_names[j])));
      }
      for (int i = 0; i < n_particles; ++i) {
        is_ok = pt.prepare(i);
        // Keep parameter validity row-wise; the LogicalRules evaluator
        // selects only the active roles for each stimulus condition.
        NumericVector row_vec(n_out_race);
        c_log_likelihood_logicalrules(lr_cols.data(), n_par_lr, expand, min_ll, is_ok, n_lR,
                                      &adapter.ctx, adapter.pdf1_ptr, adapter.cdf1_ptr,
                                      adapter.model_dfun_raw, adapter.model_pfun_raw,
                                      logicalrules_shared, &row_vec,
                                      kappa_col, tau_col, pc_col,
                                      pg_col, &lr_guess);
        result(i, _) = row_vec;
      }
      return result;
    }

    bool has_RACE_col = data.containsElementNamed("RACE");
    Rcpp::IntegerVector RACE_nacc;
    Rcpp::LogicalVector RACE_mask;
    bool has_RACE_attrs = false;
    if (has_RACE_col &&
        data.hasAttribute("RACE_nacc_by_row") &&
        data.hasAttribute("RACE_mask")) {
      RACE_nacc = data.attr("RACE_nacc_by_row");
      RACE_mask = data.attr("RACE_mask");
      has_RACE_attrs = (RACE_nacc.size() == n_trials && RACE_mask.size() == n_trials);
    }

    RaceSharedState race_shared;
    if (adapter.model_dfun_raw != nullptr && !all_finite_trials && (!has_RACE_col || has_RACE_attrs)) {
      race_shared = build_race_shared_state(data, n_trials, n_lR, winner, lR_code_vec,
                                            time_code, nogo_code,
                                            has_RACE_col, RACE_nacc, RACE_mask);
    }

    CorrDriftSharedState corr_drift_shared;
    if (adapter.ctx.corr_drift_active) {
      corr_drift_shared = build_corr_drift_shared_state(
          data, n_trials, n_lR, winner, expand, adapter.ctx,
          param_table_template, keep_names);
    }
    RDMSWTNCorrSharedState rdmswtn_corr_shared;
    if (adapter.ctx.rdmswtn_correlated) {
      rdmswtn_corr_shared = build_rdmswtn_corr_shared_state(
          data, n_trials, n_lR, winner, expand, adapter.ctx,
          param_table_template, keep_names);
    }
    const std::function<Rcpp::NumericMatrix()> corr_materialize =
        [&pt]() { return pt.materialize_reusable(); };

    for (int i = 0; i < n_particles; ++i) {
      is_ok = pt.prepare(i);
      is_ok = lr_all(is_ok, n_lR);
      if (adapter.ctx.fpe_cache) adapter.ctx.fpe_cache->new_particle();
      if (adapter.ctx.rlf_cache) adapter.ctx.rlf_cache->new_particle();
      NumericVector row_vec(n_out_race);
      if (adapter.ctx.corr_drift_active) {
        c_log_likelihood_corr_drift(
            corr_drift_shared, data, adapter.pdf1_ptr, adapter.cdf1_ptr,
            n_trials, winner, expand, min_ll, is_ok, n_lR, &adapter.ctx,
            all_finite_trials, adapter.model_dfun_raw, adapter.model_pfun_raw,
            adapter.logS_at_t_ptr,
            race_shared.valid ? &race_shared : nullptr, &row_vec,
            corr_materialize);
      } else if (adapter.ctx.rdmswtn_correlated) {
        c_log_likelihood_rdmswtn_correlated(
            rdmswtn_corr_shared, data, adapter.pdf1_ptr, adapter.cdf1_ptr,
            n_trials, winner, expand, min_ll, is_ok, n_lR, &adapter.ctx,
            all_finite_trials, adapter.model_dfun_raw, adapter.model_pfun_raw,
            adapter.logS_at_t_ptr,
            race_shared.valid ? &race_shared : nullptr, &row_vec,
            corr_materialize);
      } else {
        pars = pt.materialize_reusable();
        c_log_likelihood_race(pars, data,
                              adapter.pdf1_ptr, adapter.cdf1_ptr,
                              n_trials,
                              winner, expand, min_ll, is_ok, n_lR,
                              &adapter.ctx,
                              all_finite_trials,
                              adapter.model_dfun_raw,
                              adapter.model_pfun_raw,
                              adapter.logS_at_t_ptr,
                              race_shared.valid ? &race_shared : nullptr,
                              &row_vec);
      }
      result(i, _) = row_vec;
    }
    return result;
  }
}



// [[Rcpp::export]]
NumericMatrix get_pars_c_wrapper_oo(NumericMatrix particle_matrix,
                                    DataFrame data,
                                    NumericVector constants,
                                    List designs,
                                    List bounds,
                                    List transforms,
                                    List pretransforms,
                                    Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                    bool return_kernel_matrix = false,
                                    bool return_all_pars = false,
                                    IntegerVector kernel_output_codes = 1) {
  return get_pars_c_wrapper_oo_core(particle_matrix, data, constants, designs, bounds, transforms,
                                    pretransforms, trend, return_kernel_matrix,
                                    return_all_pars, kernel_output_codes);
}

// Map a batch of particles while keeping the design and transform objects on
// the C++ side.  The scalar wrapper above is intentionally retained for the
// public low-level API; mapped posterior summaries use this batch form to
// avoid one R <-> C++ transition per draw.
//
// [[Rcpp::export]]
NumericVector get_pars_c_batch_wrapper_oo(NumericMatrix particle_matrix,
                                          DataFrame data,
                                          NumericVector constants,
                                          List designs,
                                          List bounds,
                                          List transforms,
                                          List pretransforms,
                                          Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                          bool return_kernel_matrix = false,
                                          bool return_all_pars = false,
                                          IntegerVector kernel_output_codes = 1) {
  return get_pars_c_batch_wrapper_oo_core(
    particle_matrix, data, constants, designs, bounds, transforms, pretransforms,
    trend, return_kernel_matrix, return_all_pars, kernel_output_codes);
}

// gsl adapter for integrals - uses scalar, Rcpp-independent functions for speed
double gsl_f_race_scalar(double t, void* p) {
  auto* P = static_cast<gsl_race_params_scalar*>(p);
  if (t <= 0.0) return 0.0;
  const int w = P->winner_idx0;
  if (w < 0 || w >= P->n_lR) return 0.0;
  if (!P->isok[w]) return 0.0;

  const double* par_w = P->pars + static_cast<size_t>(w) * P->n_par;
  double out = P->pdf1(t, par_w, P->ctx);
  if (!(out > 0.0) || !emc2_isfinite(out)) return 0.0;

  for (int j = 0; j < P->n_lR; ++j) {
    if (j == w) continue;
    if (!P->isok[j]) return 0.0;
    const double* par_j = P->pars + static_cast<size_t>(j) * P->n_par;
    const double cdf_raw = P->cdf1(t, par_j, P->ctx);
    if (!R_FINITE(cdf_raw)) return 0.0;
    if (!(cdf_raw > 0.0)) continue;  // S_j(t)=1, skip multiply
    const double cdf = cdf_raw >= 1.0 ? 1.0 - 1e-15 : cdf_raw;
    out *= (1.0 - cdf);
    if (!(out > 0.0) || !emc2_isfinite(out)) return 0.0;
  }
  return out;
}

// Complete race integrand (winner pdf x loser survivors) as a log sum, for
// the shifted-integrand fallback.  Uses the same natural scalar callbacks as
// gsl_f_race_scalar, but the product is never formed on the natural scale.
double log_race_integrand_scalar(double t, const gsl_race_params_scalar* P) {
  if (t <= 0.0) return R_NegInf;
  const int w = P->winner_idx0;
  if (w < 0 || w >= P->n_lR) return R_NegInf;
  if (!P->isok[w]) return R_NegInf;

  const double* par_w = P->pars + static_cast<size_t>(w) * P->n_par;
  const double pdf = P->pdf1(t, par_w, P->ctx);
  if (!(pdf > 0.0) || !emc2_isfinite(pdf)) return R_NegInf;
  double log_out = std::log(pdf);

  for (int j = 0; j < P->n_lR; ++j) {
    if (j == w) continue;
    if (!P->isok[j]) return R_NegInf;
    const double* par_j = P->pars + static_cast<size_t>(j) * P->n_par;
    const double cdf_raw = P->cdf1(t, par_j, P->ctx);
    if (!R_FINITE(cdf_raw)) return R_NegInf;
    if (!(cdf_raw > 0.0)) continue;  // S_j(t)=1
    if (cdf_raw >= 1.0) return R_NegInf;
    log_out += std::log1p(-cdf_raw);
  }
  return log_out;
}

// Shifted natural integrand for GSL: exp(log integrand - log_scale).  Entered
// only after the natural integrand has failed its guard (the natural product
// under/overflowed while the log integrand is representable).
double gsl_f_race_scalar_logshift(double t, void* p) {
  auto* P = static_cast<gsl_race_params_scalar*>(p);
  const double log_out = log_race_integrand_scalar(t, P);
  if (log_out == R_NegInf || ISNAN(log_out)) return 0.0;
  // The shift keeps the integrand O(1); cap the exponent defensively so a
  // poorly-placed pilot cannot hand GSL an Inf.
  return std::exp(std::fmin(log_out - P->log_scale, 700.0));
}

// Log survivor and cdf of the race at time t:
//   log S(t) = sum_k log(1 - F_k(t))
//
// "rowmajor" here means per-trial parameters for the n_lR accumulators are
// packed as a contiguous buffer:
//   pars_rowmajor[k * n_par + c]
//
// This avoids repeatedly indexing into an Rcpp::NumericMatrix inside tight
// loops and matches the representation needed for the scalar/GSL integrands.

inline double log_survivor_rowmajor(double t,
                                    const double* pars_rowmajor,
                                    const int* isok_int,
                                    int n_lR,
                                    int n_par,
                                    RaceCdf1Fun cdf1,
                                    void* ctx) {
  if (t == R_PosInf) {
    auto* race_ctx = static_cast<ContextForRaceModels*>(ctx);
    if (race_ctx && race_ctx->defective_upper_tail) {
      double logS = 0.0;
      for (int k = 0; k < n_lR; ++k) {
        if (!isok_int[k]) return R_NegInf;
        const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
        double cdf_inf = cdf1(R_PosInf, par_k, ctx);
        cdf_inf = clamp_cdf01_race(cdf_inf);
        const double ll = safe_log1m_race(cdf_inf);
        if (!emc2_isfinite(ll)) return R_NegInf;
        logS += ll;
      }
      return logS;
    }
    return R_NegInf;
  }
  double logS = 0.0;
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = safe_log1m_race(cdf);
    if (!emc2_isfinite(ll)) return R_NegInf;
    logS += ll;
  }
  return logS;
}

inline double log_cdf_rowmajor(double t,
                               const double* pars_rowmajor,
                               const int* isok_int,
                               int n_lR,
                               int n_par,
                               RaceCdf1Fun cdf1,
                               void* ctx) {
  if (t == R_PosInf) {
      auto* race_ctx = static_cast<ContextForRaceModels*>(ctx);
      if (race_ctx && race_ctx->defective_upper_tail) {
          double logC = 0.0;
          for (int k = 0; k < n_lR; ++k) {
            if (!isok_int[k]) return R_NegInf;
            const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
            double cdf_inf = cdf1(R_PosInf, par_k, ctx);
            cdf_inf = clamp_cdf01_race(cdf_inf);
            const double ll = std::log(cdf_inf);
            if (!emc2_isfinite(ll)) return R_NegInf;
            logC += ll;
          }
          return logC;
     }
    return 0.0;
  }
  double logC = 0.0;
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = std::log(cdf);
    if (!emc2_isfinite(ll)) return R_NegInf;
    logC += ll;
  }
  return logC;
}

inline double log_min_density_rowmajor(double t,
                                       const double* pars_rowmajor,
                                       const int* isok_int,
                                       int n_lR,
                                       int n_par,
                                       RacePdf1Fun pdf1,
                                       RaceCdf1Fun cdf1,
                                       void* ctx,
                                       double* logS_k) {
  // Log density of the minimum (no known winner) at time t:
  //   f_min(t) = sum_k f_k(t) * prod_{j != k} (1 - F_j(t))
  //
  // Work in log space:
  //   log f_min(t) = logsumexp_k [ log f_k(t) + sum_{j != k} log(1 - F_j(t)) ].
  //
  // logS_k is a per-accumulator scratch buffer (reused across calls) holding
  // log(1 - F_k(t)) for the current t, to avoid reallocations in the "other
  // trial" loop.
  if (!(t > 0.0) || !emc2_isfinite(t)) return R_NegInf;
  double logS_all = 0.0;
  std::fill(logS_k, logS_k + n_lR, R_NegInf);
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = safe_log1m_race(cdf);
    if (!emc2_isfinite(ll)) return R_NegInf;
    logS_k[static_cast<size_t>(k)] = ll;
    logS_all += ll;
  }
  double out = R_NegInf;
  for (int k = 0; k < n_lR; ++k) {
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    const double pdf = pdf1(t, par_k, ctx);
    if (!(pdf > 0.0) || !emc2_isfinite(pdf)) continue;
    const double term = std::log(pdf) + (logS_all - logS_k[static_cast<size_t>(k)]);
    out = log_sum_exp(out, term);
  }
  return out;
}


double integrate_for_kth_winner_rowmajor_cpp(
    int k_winner_idx, // 1-based
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
    gsl_integration_workspace* w) {
  
  // Integrate the k-th winner density over an interval [low, upp] for a single
  // unique trial, using rowmajor buffers (raw pointers). This is the fast path
  // used by truncation/censoring normalisers and by go/no-go branches; it avoids
  // Rcpp object traffic inside the GSL callback.
  if (low >= upp && !(low == 0 && upp == R_PosInf)) return R_NegInf;
  if (k_winner_idx < 1 || k_winner_idx > n_lR_j) return R_NegInf;
  if (w == nullptr) Rcpp::stop("integrate_for_kth_winner_rowmajor_cpp: GSL workspace is null.");

  // Start at the integrand's support edge, not at the caller's nominal bound.
  // The winner's density is identically zero below its own non-decision time,
  // so an interval like the go/no-go withheld branch's [LT, UC] = [0, UC] hands
  // GSL a long dead stretch.  A Gauss-Kronrod panel covering only that stretch
  // returns zero WITH a near-zero error estimate, so the adaptive routine can
  // declare convergence before it ever resolves the peak further right --
  // silently returning a value ~1% low at isolated parameter values (which the
  // trial counts of a compressed dadm then amplify into tens of nats).  Clipping
  // to t0 removes the dead zone entirely, and is strictly cheaper.  t0_index
  // is -1 only for models whose t0 is not a pure additive shift.
  {
    const ContextForRaceModels* rctx =
      static_cast<const ContextForRaceModels*>(model_specific_context);
    if (rctx != nullptr && rctx->t0_index >= 0 && rctx->t0_index < n_par) {
      const double t0_w =
        pars_rowmajor[static_cast<size_t>(k_winner_idx - 1) * n_par + rctx->t0_index];
      if (R_FINITE(t0_w) && t0_w > low) {
        low = t0_w;
        if (low >= upp) return R_NegInf;   // support lies outside the interval
      }
    }
  }
  
  gsl_function F;
  gsl_race_params_scalar params_struct;
  params_struct.pars = pars_rowmajor;
  params_struct.n_lR = n_lR_j;
  params_struct.n_par = n_par;
  params_struct.winner_idx0 = k_winner_idx - 1;
  params_struct.isok = isok_int;
  params_struct.pdf1 = pdf1;
  params_struct.cdf1 = cdf1;
  params_struct.ctx = model_specific_context;
  
  F.function = &gsl_f_race_scalar;
  F.params = &params_struct;
  
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status;
  double result = 0.0;
  double error = 0.0;
  auto run_integral = [&](double abs_tol, double rel_tol, size_t limit) -> int {
    if (upp == R_PosInf) {
      if (low < 0) low = 0; // QAGIU requires a >= 0
      if (low >= R_PosInf) {
        result = 0.0;
        error = 0.0;
        return GSL_SUCCESS;
      }
      return gsl_integration_qagiu(&F, low, abs_tol, rel_tol, limit, w, &result, &error);
    }
    if (gsl_ctl.try_qng_first_finite) {
      size_t neval = 0;
      const int qng_status = gsl_integration_qng(&F, low, upp, abs_tol, rel_tol,
                                                 &result, &error, &neval);
      if (qng_status == GSL_SUCCESS && R_FINITE(result)) return GSL_SUCCESS;
      return gsl_integration_qag(&F, low, upp, abs_tol, rel_tol, limit, gsl_ctl.qag_key,
                                 w, &result, &error);
    }
    return gsl_integration_qags(&F, low, upp, abs_tol, rel_tol, limit, w, &result, &error);
  };
  status = run_integral(gsl_ctl.abs_tol, gsl_ctl.rel_tol, gsl_ctl.limit);
  // Retry once with tighter controls if the fast pass fails.
  if (status != GSL_SUCCESS) {
    status = run_integral(gsl_ctl.retry_abs_tol, gsl_ctl.retry_rel_tol, gsl_ctl.retry_limit);
  }

  if (status == GSL_SUCCESS && result > 0.0 && R_FINITE(result)) {
    gsl_set_error_handler(old_handler);
    return std::log(result);
  }

  // Natural-integrand failure: the winner-pdf x loser-survivor product
  // under/overflowed pointwise even though the log integrand may still be
  // representable.  Pilot the log integrand to pick a shift, then integrate
  // the shifted natural integrand g(t) = exp(log g(t) - log_scale) so the
  // GSL API still consumes natural values.
  double max_log = R_NegInf;
  {
    const double lo = std::fmax(low, 0.0);
    if (upp == R_PosInf) {
      static const double offsets[] = {0.05, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0};
      const double base = std::fmax(lo, 0.0);
      for (double off : offsets) {
        const double lg = log_race_integrand_scalar(base + off, &params_struct);
        if (lg > max_log) max_log = lg;
      }
    } else if (upp > lo) {
      static const double fracs[] = {0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98};
      for (double f : fracs) {
        const double lg = log_race_integrand_scalar(lo + f * (upp - lo), &params_struct);
        if (lg > max_log) max_log = lg;
      }
    }
  }
  if (!R_FINITE(max_log)) {
    gsl_set_error_handler(old_handler);
    return R_NegInf;
  }

  params_struct.log_scale = max_log;
  F.function = &gsl_f_race_scalar_logshift;
  status = run_integral(gsl_ctl.abs_tol, gsl_ctl.rel_tol, gsl_ctl.limit);
  if (status != GSL_SUCCESS) {
    status = run_integral(gsl_ctl.retry_abs_tol, gsl_ctl.retry_rel_tol, gsl_ctl.retry_limit);
  }
  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS) return R_NegInf;
  if (!(result > 0.0) || !R_FINITE(result)) return R_NegInf;
  return max_log + std::log(result);
}

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
                                         GslWorkspacePtr& workspace) {
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());

  // Global kill: the shared clock is applied as a per-trial factor at the
  // likelihood-assembly sites, NOT inside cdf1 (apply_lk_to_racers is false),
  // so the per-racer survivor product below would normalise by
  //   P(no racer crossed by t)
  // instead of the quantity the numerator is conditioned on,
  //   P(no response by t) = 1 - integral_0^t f_race(u) S_K(u) du.
  // Those differ by the paths where the kill fired first, which are also
  // non-responses; multiplying in S_K(t) does NOT fix it.
  //
  // With one accumulator the correct survivor is the killed sub-CDF, obtained
  // by switching the kill back on for the racer -- the same trick the
  // single-accumulator omission branch uses.  With several accumulators a
  // single shared clock does not factorise into per-racer survivors, so no
  // product form exists to correct and the race integral would be required.
  auto* race_ctx = static_cast<ContextForRaceModels*>(model_specific_context);
  const bool global_kill_trunc =
      race_ctx && race_ctx->has_global_kill() && ((LT != 0.0) || R_FINITE(UT));
  if (global_kill_trunc && n_lR != 1) {
    Rcpp::stop("erlang_type = \"global_kill\" does not support truncation "
               "(LT/UT) with more than one accumulator: a shared kill clock "
               "does not factorise into per-accumulator survivors, so the "
               "truncation normaliser would be wrong. Use "
               "erlang_type = \"local_kill\" (equivalent when the "
               "accumulators are independent) or remove the truncation "
               "bounds.");
  }
  // Restores apply_lk_to_racers even if a kernel throws.
  struct KillOnScope {
    ContextForRaceModels* ctx = nullptr;
    bool saved = false;
    ~KillOnScope() { if (ctx) ctx->apply_lk_to_racers = saved; }
  } kill_scope;
  if (global_kill_trunc) {
    kill_scope.ctx = race_ctx;
    kill_scope.saved = race_ctx->apply_lk_to_racers;
    race_ctx->apply_lk_to_racers = true;
  }

  // When LT == 0, every proper distribution has CDF(0) = 0 → S(0) = 1 → log = 0.
  // Skip the n_lR scalar cdf1 calls in this common case.
  double logS_LT;
  if (LT == 0.0) {
    logS_LT = 0.0;
  } else {
    logS_LT = log_survivor_rowmajor(LT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
    if (!R_FINITE(logS_LT)) return R_NegInf;
  }
  if (UT == R_PosInf) {
    // UT == +Inf means no upper truncation. A defective model's intrinsic
    // omission mass sits at +Inf, i.e. >= LT, so it stays inside the retained
    // window [LT, Inf): the truncation normaliser only ever excludes finite
    // density in [0, LT). This matches the batch path (which leaves
    // logS_UT == R_NegInf for UT==Inf, defective or not).
    return logS_LT;
  }

  const double logS_UT = log_survivor_rowmajor(UT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
  // The retained sample space under [LT, UT]. make_missing() cuts only FINITE
  // RTs outside the window (R/make_data.R), so an intrinsic never-finish
  // outcome survives upper truncation: it is a T = +Inf atom, not a slow
  // response that got cut off. For a defective model the retained mass is
  //   {LT <= T <= UT} u {T = +Inf}
  // and the atom P(T = Inf) = S(Inf) = prod_j (1 - h_j) must be added back,
  // otherwise the finite density is renormalised to one on its own while the
  // omission score keeps its (undivided) atom -- the conditional distribution
  // then integrates to more than one, badly so when F(UT) < S(Inf). Proper
  // models have S(Inf) = 0 and log_defective_atom stays -Inf, so this is a
  // no-op for them.
  const double log_defective_atom =
      (race_ctx && race_ctx->defective_upper_tail)
        ? log_survivor_rowmajor(R_PosInf, pars_rowmajor, isok_int, n_lR, n_par,
                                cdf1, model_specific_context)
        : R_NegInf;
  double logP = log_sum_exp(log_diff_exp(logS_LT, logS_UT), log_defective_atom);
  if (R_FINITE(logP) && logP > log_prob_eps) return logP;

  // Only falls back to GSL integration if the analytic trick fails (e.g. due to catastrophic cancellation)
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
  double log_total = R_NegInf;
  for (int k_win = 1; k_win <= n_lR; ++k_win) {
    const double log_k = integrate_for_kth_winner_rowmajor_cpp(k_win,
                                                               pars_rowmajor,
                                                               isok_int,
                                                               LT,
                                                               UT,
                                                               pdf1,
                                                               cdf1,
                                                               n_lR,
                                                               n_par,
                                                               gsl_ctl,
                                                               model_specific_context,
                                                               w);
    log_total = log_sum_exp(log_total, log_k);
  }
  // The integrals cover the finite window only; add the retained +Inf atom for
  // the same reason as the analytic branch above.
  log_total = log_sum_exp(log_total, log_defective_atom);
  if (R_FINITE(log_total) && log_total > log_prob_eps) return log_total;
  return R_NegInf;
}


// Main C++ function for censored/truncated race likelihood calculation
// This function is now the unified entry point for all race models (LBA, RDM, LNR),
// whether they are standard or explicitly handling censoring/truncation.
// It uses batching for finite RTs and iterative processing for others (censored/NA RTs).
double c_log_likelihood_race(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, covering all dadm rows for that particle
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions, structured for all accumulators
    RacePdf1Fun pdf1,                       // Scalar PDF (for GSL)
    RaceCdf1Fun cdf1,                       // Scalar CDF (for GSL)
    const int n_trials,
    LogicalVector winner,
    Rcpp::IntegerVector expand,  // Vector for expanding unique LLs to full trial count
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector isok,   // Parameter validity for each row in 'pars' matrix
    int n_lR,                              // Number of accumulators in the race (must be > 0 if data exists)
    void* model_context_for_funcs,          // Context for model functions (e.g., contains posdrift for LBA)
    bool all_finite_trials,             // Data-only hint: all trials finite/in-bounds/known response
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,            // batch log-survivor at scalar t (for truncation norms)
    RaceSharedState* shared,             // optional pre-computed per-data state (nullptr = compute per-call)
    NumericVector* trial_ll_out,
    bool apply_truncation_correction
) {

  // Reuse one workspace per thread across particles/model fits.
  static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
  GslWorkspacePtr& workspace = workspace_tls;

  const bool use_shared = (shared != nullptr && shared->valid);

  // Fetch censoring and truncation columns only when we actually need them.
  // When all_finite_trials=true the caller guarantees LT=0, UT=Inf, no censoring —
  // deferring these reads avoids O(n_trials) allocation+fill per particle call.
  // When use_shared=true the pre-read vectors are reused directly (no allocation).
  const bool may_need_ct = !all_finite_trials;
  Rcpp::NumericVector LT, UT, LC, UC;
  if (may_need_ct) {
    if (use_shared) {
      LT = shared->LT_vec;   // shared reference, no copy
      UT = shared->UT_vec;
      LC = shared->LC_vec;
      UC = shared->UC_vec;
    } else {
      LT = get_col_with_default(dadm, "LT", 0.0);
      UT = get_col_with_default(dadm, "UT", R_PosInf);
      LC = get_col_with_default(dadm, "LC", 0.0);
      UC = get_col_with_default(dadm, "UC", R_PosInf);
    }
  }
  // Default controls with a stricter retry for difficult integrals.
  GslIntegrationControls gsl_ctl = default_gsl_controls();
  gsl_ctl.try_qng_first_finite = true;   // try fixed-point QNG before adaptive QAG on finite intervals
  gsl_ctl.qag_key = GSL_INTEG_GAUSS21;  // fallback rule when QNG fails
  // 1e-4 was enough for the integral itself but not for the log-likelihood SURFACE:
  // the adaptive error estimate is optimistic on a peaked race integrand, so isolated
  // parameter values came out ~1e-3 off and, once multiplied by a compressed dadm's
  // trial counts, left visible roughness for any t0-scanning method to trip over.
  gsl_ctl.rel_tol = 1e-5;
  int n_lR_j = n_lR;
  Rcpp::NumericVector rts_dadm = dadm["rt"];
  Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
  // Use pre-allocated std::vector when shared state is available; avoids per-particle
  // R-heap allocation.  The buffer does NOT need re-initialisation between particles:
  // dfun_raw/pfun_raw overwrite every finite-trial slot before it is read, and
  // other-trial slots are computed independently in the "other trials" loop below.
  std::vector<double> lds_local;
  double* lds_ptr;
  if (use_shared && static_cast<int>(shared->res_buf.size()) == n_trials) {
    lds_ptr = shared->res_buf.data();
  } else {
    lds_local.assign(static_cast<size_t>(n_trials), min_ll);
    lds_ptr = lds_local.data();
  }
  // If a RACE column exists, set parameters of accumulators not present on a
  // given trial to NA so the density functions return zero for them. This
  // mirrors logic from the old c_log_likelihood_race implementation.
  bool has_RACE_col = dadm.containsElementNamed("RACE");
  Rcpp::IntegerVector RACE;
  Rcpp::LogicalVector RACE_mask;
  if (has_RACE_col) {
    if (dadm.hasAttribute("RACE_nacc_by_row") && dadm.hasAttribute("RACE_mask")) {
      RACE = dadm.attr("RACE_nacc_by_row");
      RACE_mask = dadm.attr("RACE_mask");
      for (int row = 0; row < pars.nrow(); ++row) {
        if (!RACE_mask[row]) std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
      }
    } else {
      has_RACE_col = false;
    }
  }
  if (n_trials == 0) return 0.0; // No data, no likelihood
  
  if (n_lR <= 0) Rcpp::stop("c_log_likelihood_race: n_lR must be positive and correctly determined before this call.");
  if (n_trials % n_lR != 0) Rcpp::stop("c_log_likelihood_race: dadm nrows not a multiple of n_lR.");
  
  // Here we check for a pC parameter corresponding to probability of contaminant OMISSION.
  // The column index is data-structure-fixed so we cache it in shared state after the
  // first search.  The "all zeroes" check still runs per particle (values change).
  bool use_pC = false;
  int pc_col = -1;
  if (use_shared && shared->pc_col != -2) {
    // Fast path: use cached result from a previous particle call
    pc_col = shared->pc_col;
    use_pC = (pc_col >= 0);
  } else {
    Rcpp::List dimnames = pars.attr("dimnames");
    Rcpp::CharacterVector colnames = as<Rcpp::CharacterVector>(dimnames[1]);
    for (int j = 0; j < colnames.size(); ++j) {
      if (as<std::string>(colnames[j]) == "pContaminant") {
        pc_col = j;
        use_pC = true;
        break;
      }
    }
    if (use_shared) shared->pc_col = pc_col; // cache for subsequent particles
  }
  // Per-particle check: pC column present but all values zero → treat as absent
  if (use_pC) {
    bool all_zero = true;
    for (int i = 0; i < pars.nrow(); ++i) {
      if (pars(i, pc_col) != 0.0) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      use_pC = false;
    }
  }

  // pGuess: the same column-cache and all-zeroes shortcut, run in parallel with
  // pContaminant so the two can never disagree about which trials they touch.
  // The uniform guess kernel itself is data-fixed (resolve_guess_window() in R).
  GuessKernel guess;
  bool use_pG = false;
  int pg_col = -1;
  if (use_shared && shared->pg_col != -2) {
    pg_col = shared->pg_col;
    guess = shared->guess;
  } else {
    Rcpp::List dimnames = pars.attr("dimnames");
    Rcpp::CharacterVector colnames = as<Rcpp::CharacterVector>(dimnames[1]);
    for (int j = 0; j < colnames.size(); ++j) {
      if (as<std::string>(colnames[j]) == "pGuess") { pg_col = j; break; }
    }
    guess = resolve_guess_kernel(dadm);
    guess.pg_col = pg_col;
    if (use_shared) { shared->pg_col = pg_col; shared->guess = guess; }
  }
  use_pG = guess.active();
  if (use_pG) {
    bool all_zero = true;
    for (int i = 0; i < pars.nrow(); ++i) {
      if (pars(i, pg_col) != 0.0) { all_zero = false; break; }
    }
    if (all_zero) use_pG = false;
  }

  int n_unique_trials = n_trials / n_lR;
  // Use std::vector to avoid per-particle R-heap allocation overhead.
  std::vector<double> ll_unique(static_cast<size_t>(n_unique_trials), min_ll);
  // Raw pointer into pars matrix: col-major layout, element (row,col) = pars_cm_ptr[col*n_trials+row]
  const double* pars_cm_ptr = pars.begin();

  // Parameter matrix and validity vector checks
  if (pars.nrow() != n_trials) {
    Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_trials: " << n_trials << std::endl;
    Rcpp::stop("c_log_likelihood_race: pars matrix dimensions do not match total dadm rows.");
  }
  if (isok.size() != pars.nrow()) {
    Rcpp::stop("c_log_likelihood_race: isok size does not match "
                 "pars matrix rows.");
  }
  if (winner.size() != pars.nrow()) {
    Rcpp::stop("c_log_likelihood_race: isok size does not match pars matrix rows.");
  }
  const int n_par = pars.ncol();
  // Scratch buffers.
  //
  // These were fixed stack arrays, which forced a size limit and made the
  // frame ~280 kB.  `static thread_local` vectors -- the same idiom this
  // function already uses for the GSL workspace -- keep their capacity for the
  // life of the process, so the allocation is paid once rather than per call,
  // the frame stays small, and there is no limit to hit.  That matters: n_par
  // is the number of model parameter *types* (v, B, A, t0, ... for a race),
  // which a trend model extends one name at a time via update_model_trend().
  //
  // `.assign()` rather than `.resize()`: resize leaves previously-used slots
  // holding the last call's values, and every buffer here is written only up
  // to the current trial's width.  Zeroing keeps a read past that width
  // deterministic instead of quietly plausible.
  static thread_local std::vector<const double*> cols_view_tls;
  static thread_local std::vector<double> pars_rowmajor_tls;
  static thread_local std::vector<int> isok_int_tls;
  static thread_local std::vector<double> logS_k_tls;
  cols_view_tls.assign(static_cast<size_t>(n_par), nullptr);
  pars_rowmajor_tls.assign(static_cast<size_t>(n_lR) * n_par, 0.0);
  isok_int_tls.assign(static_cast<size_t>(n_lR), 0);
  logS_k_tls.assign(static_cast<size_t>(n_lR), R_NegInf);
  const double** cols_view = cols_view_tls.data();
  double* pars_rowmajor_buffer = pars_rowmajor_tls.data();
  int* isok_int_buffer = isok_int_tls.data();
  double* logS_k_buffer = logS_k_tls.data();

  // Column-pointer view of the materialized pars matrix for the raw kernels
  // (this path keeps the matrix: the RACE NA-fill above writes into it).
  for (int c = 0; c < n_par; ++c) {
    cols_view[c] = pars_cm_ptr + static_cast<size_t>(c) * n_trials;
  }

  // One accumulator row of parameters, staged for the scalar cdf1/pdf1
  // kernels.  Six branches below need one; they are in mutually exclusive
  // paths, but two of them hand the buffer to a lambda that outlives the
  // filling loop, so each site gets its own slot rather than sharing one.
  // Padded to at least 64 because a kernel may *fetch* (not dereference) an
  // optional trailing column this model variant does not carry -- the fixed
  // arrays these replace gave that slack implicitly.
  const size_t par_row_width = static_cast<size_t>(std::max(n_par, 64));

  // fill_trial_buffers: copies one trial's params into row-major scratch (for GSL/rowmajor helpers).
  // Uses raw column-major pointer to avoid Rcpp subscript overhead.
  auto fill_trial_buffers = [&](int start_row_idx, int n_lR_j) {
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start_row_idx + k;
      isok_int_buffer[static_cast<size_t>(k)] = isok[row] ? 1 : 0;
      for (int c = 0; c < n_par; ++c)
        pars_rowmajor_buffer[static_cast<size_t>(k) * n_par + c] =
            pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
    }
  };

  // Fast path hint: computed once per calc_ll call (data-only), to avoid re-scanning per particle.
  const bool use_full_finite_batch = all_finite_trials;
  Rcpp::LogicalVector finite_rt_mask;
  std::vector<int> finite_rt_unique_trial_indices;
  std::vector<int> other_unique_trial_indices;
  std::vector<int> active_nogo_trial_mask;
  if (!use_full_finite_batch) {
    if (use_shared && static_cast<int>(shared->finite_mask.size()) == n_trials) {
      // Fast path: re-use pre-read partition from shared state (no attr lookup, no copy)
      finite_rt_mask = shared->finite_mask;
      finite_rt_unique_trial_indices = shared->finite_unique_idx;
      other_unique_trial_indices     = shared->other_unique_idx;
      active_nogo_trial_mask         = shared->active_nogo_trial_mask;
    } else {
      const bool has_partition_attrs =
        dadm.hasAttribute("finite_rt_mask") &&
        dadm.hasAttribute("finite_rt_unique_trial_indices") &&
        dadm.hasAttribute("other_unique_trial_indices") &&
        dadm.hasAttribute("active_nogo_trial_mask");

      if (has_partition_attrs) {
        finite_rt_mask = dadm.attr("finite_rt_mask");
        Rcpp::IntegerVector finite_attr = dadm.attr("finite_rt_unique_trial_indices");
        Rcpp::IntegerVector other_attr = dadm.attr("other_unique_trial_indices");
        Rcpp::LogicalVector nogo_attr = dadm.attr("active_nogo_trial_mask");
        finite_rt_unique_trial_indices.assign(finite_attr.begin(), finite_attr.end());
        other_unique_trial_indices.assign(other_attr.begin(), other_attr.end());
        active_nogo_trial_mask.assign(nogo_attr.begin(), nogo_attr.end());
      } else {
        // Fallback path for direct calc_ll/calc_ll_oo callers that bypass
        // .cache_ll_data_attrs(): derive finite/other unique-trial partitions.
        finite_rt_mask = Rcpp::LogicalVector(n_trials, false);
        active_nogo_trial_mask.assign(static_cast<size_t>(n_unique_trials), 0);
        const SEXP lR_sexp = dadm["lR"];
        const Rcpp::IntegerVector lR_code(lR_sexp);
        const Rcpp::CharacterVector lR_levels = lR_code.attr("levels");
        int nogo_code = -1;
        for (int i = 0; i < lR_levels.size(); ++i) {
          if (Rcpp::as<std::string>(lR_levels[i]) == "nogo") {
            nogo_code = i + 1;
            break;
          }
        }
        for (int unique_trial_idx = 0; unique_trial_idx < n_unique_trials; ++unique_trial_idx) {
          const int start_row_idx = unique_trial_idx * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
          const double rt_j = rts_dadm[start_row_idx];
          const int R_j_idx = R_idxs_dadm[start_row_idx];
          if (nogo_code > 0) {
            bool has_active_nogo = false;
            for (int k = 0; k < n_lR_j; ++k) {
              const int row = start_row_idx + k;
              if (has_RACE_col && static_cast<int>(RACE_mask.size()) == n_trials && !RACE_mask[row]) continue;
              if (lR_code[row] == nogo_code) { has_active_nogo = true; break; }
            }
            active_nogo_trial_mask[static_cast<size_t>(unique_trial_idx)] = has_active_nogo ? 1 : 0;
          }
          if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx != NA_INTEGER) {
            finite_rt_unique_trial_indices.push_back(unique_trial_idx);
            for (int k = 0; k < n_lR_j; ++k) finite_rt_mask[start_row_idx + k] = true;
          } else {
            other_unique_trial_indices.push_back(unique_trial_idx);
          }
        }
      }
    }
  }
  double log_Z_this = 0;  // Default inv_Z if no truncation. Should never be used but here as a precaution.
  // Scratch interval bounds for per-trial branches below. Keeping these as
  // mutable locals makes the censoring/truncation logic read like "set bounds,
  // then integrate / take survivor", rather than threading LTj/UTj/LCj/UCj
  // through every helper call.
  double lower_for_trial = 0;
  double upper_for_trial = R_PosInf;
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  int time_code = ctx->time_code;
  int nogo_code = ctx->nogo_code;
  // Use pre-resolved codes from shared state when available.
  if (use_shared) {
    if (time_code == -2) time_code = shared->time_code;
    if (nogo_code == -2) nogo_code = shared->nogo_code;
  }
  // Fall back to scanning the factor levels when still unresolved.
  if (time_code == -2 || nogo_code == -2) {
    const SEXP lR_sexp = dadm["lR"];
    if (Rf_inherits(lR_sexp, "factor")) {
      Rcpp::IntegerVector lR_code_tmp(lR_sexp);
      Rcpp::CharacterVector lR_levels_tmp = lR_code_tmp.attr("levels");
      if (time_code == -2) time_code = -1;
      if (nogo_code == -2) nogo_code = -1;
      for (int j = 0; j < lR_levels_tmp.size(); ++j) {
        std::string lev = Rcpp::as<std::string>(lR_levels_tmp[j]);
        if (lev == "time") time_code = j + 1;
        else if (lev == "nogo") nogo_code = j + 1;
      }
    } else {
      if (time_code == -2) time_code = -1;
      if (nogo_code == -2) nogo_code = -1;
    }
    if (ctx->time_code == -2) ctx->time_code = time_code;
    if (ctx->nogo_code == -2) ctx->nogo_code = nogo_code;
  }
  const bool has_time = (time_code > 0);

  // For any Erlang variant, disable kill bookkeeping when lambda is identically
  // zero for this particle so kernels use the standard (non-mixture) forms.
  const bool any_erlang_ctx = ctx && (ctx->mean_k_index >= 0 || ctx->mean_g_index >= 0) &&
    (ctx->is_global_kill || ctx->is_local_kill || ctx->is_local_guess || ctx->is_local_kill_guess);
  if (any_erlang_ctx) {
    bool lambda_active = false;
    const double* mean_k_ptr = (ctx->mean_k_index >= 0)
      ? (pars_cm_ptr + static_cast<size_t>(ctx->mean_k_index) * n_trials) : nullptr;
    const double* mean_g_ptr = (ctx->mean_g_index >= 0)
      ? (pars_cm_ptr + static_cast<size_t>(ctx->mean_g_index) * n_trials) : nullptr;
    for (int row = 0; row < n_trials; ++row) {
      if (!isok[row]) continue;
      const double lk = mean_k_ptr ? erlang_lambda_from_mean(mean_k_ptr[row], ctx->kill_shape) : 0.0;
      const double lg = mean_g_ptr ? erlang_lambda_from_mean(mean_g_ptr[row], ctx->kill_shape) : 0.0;
      if ((emc2_isfinite(lk) && lk > 1e-12) || (emc2_isfinite(lg) && lg > 1e-12)) {
        lambda_active = true; break;
      }
    }
    ctx->kill_active = lambda_active;
  } else if (ctx) {
    // If erlang column indices are registered but no erlang type flag is active
    // (e.g. RDMSWTN with erlang_type="none"), disable kill to prevent scalar
    // kernels from reading past the parameter block at par[mean_k_index].
    const bool has_erlang_indices = (ctx->mean_k_index >= 0 || ctx->mean_g_index >= 0);
    const bool has_erlang_type = ctx->is_global_kill || ctx->is_local_kill ||
                                 ctx->is_local_guess  || ctx->is_local_kill_guess;
    ctx->kill_active = !has_erlang_indices || has_erlang_type;
  }

  ContextForRaceModels dense_ctx = *ctx;
  dense_ctx.min_lik_for_pdf = 0.0;
  void* dense_ctx_ptr = static_cast<void*>(&dense_ctx);
  const bool defective_upper_tail = ctx->defective_upper_tail;
  const bool global_omission_active = (ctx && ctx->is_global_kill && ctx->kill_active);

  // Cache 20-point Gauss-Legendre nodes/weights for the global-kill t=Inf branch.
  // This avoids per-trial Rcpp List/Vector allocation and repeated node transforms.
  std::array<double, 20> gl_u{};
  std::array<double, 20> gl_w_half{};
  std::array<double, 20> gl_jac_erlang2{};
  std::array<double, 20> gl_log1m_u{};
  bool gl20_ready = false;
  if (global_omission_active) {
    const int n_nodes = 20;
    const GLRule& gl = gl_get_rule(n_nodes);
    const std::vector<double>& gl_nodes = gl.x;
    const std::vector<double>& gl_weights = gl.w;
    for (int j = 0; j < n_nodes; ++j) {
      const double u = 0.5 * (gl_nodes[j] + 1.0);
      gl_u[static_cast<size_t>(j)] = u;
      gl_w_half[static_cast<size_t>(j)] = 0.5 * gl_weights[j];
      gl_log1m_u[static_cast<size_t>(j)] = std::log1p(-u);
      gl_jac_erlang2[static_cast<size_t>(j)] = -gl_log1m_u[static_cast<size_t>(j)];
    }
    gl20_ready = true;
  }

  // Cache log_surv_cm(+Inf) per unique trial across all models/branches.
  std::vector<double> global_log_surv_inf_by_trial;
  std::vector<uint8_t> global_log_surv_inf_ready;
  global_log_surv_inf_by_trial.assign(static_cast<size_t>(n_unique_trials), R_NaN);
  global_log_surv_inf_ready.assign(static_cast<size_t>(n_unique_trials), 0);
  
  // log_surv_cm: log S(t) = sum_k log(1-F_k(t)) read directly from column-major pars.
  // Used by the "other trials" path for analytical survivor calls — avoids the
  // row-major copy that fill_trial_buffers performs.
  auto log_surv_cm = [&](double t, int start_row_idx, int n_lR_j) -> double {
    auto* race_ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
    const int unique_trial_idx = start_row_idx / n_lR;
    if (t == R_PosInf) {
      if (global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)]) {
        return global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)];
      }
      // Proper race models have S(Inf)=0 => log S(Inf) = -Inf.
      // Defective-tail models (e.g., LBAIO / BAwL) retain mass at +Inf.
      if (race_ctx && race_ctx->is_global_kill && race_ctx->kill_active) {
        // Global omission probability: Integral f_K(u) * S_race(u) du on [0, Inf).
        // Use Gauss-Legendre quadrature (20 nodes) with mapping u -> t = -1/lam * log(1-u).
        const double* mean_k_ptr = pars_cm_ptr + race_ctx->mean_k_index * n_trials;
        const double lam = erlang_lambda_from_mean(mean_k_ptr[start_row_idx],
                                                   race_ctx->kill_shape);
        if (lam <= 0.0) {
          global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = R_NegInf;
          global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
          return R_NegInf;
        }
        // Single-accumulator global-kill case is analytic:
        // P(omission) = 1 - P(hit by Inf) where P(hit by Inf) is model CDF at +Inf.
        if (n_lR_j == 1 && isok[start_row_idx]) {
          static thread_local std::vector<double> par_buf_inf_tls1;
          par_buf_inf_tls1.assign(par_row_width, 0.0);
          double* par_buf_inf = par_buf_inf_tls1.data();
          for (int c = 0; c < n_par; ++c)
            par_buf_inf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + start_row_idx];
          ContextForRaceModels* race_ctx_local = static_cast<ContextForRaceModels*>(model_context_for_funcs);
          const bool old_apply = race_ctx_local->apply_lk_to_racers;
          race_ctx_local->apply_lk_to_racers = true;
          double F_inf = cdf1(R_PosInf, par_buf_inf, race_ctx_local);
          race_ctx_local->apply_lk_to_racers = old_apply;
          F_inf = clamp_cdf01_race(F_inf);
          const double ans = safe_log1m_race(F_inf);
          global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = ans;
          global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
          return ans;
        }
        // Build row-major trial buffer once; avoids re-packing parameter rows per node.
        fill_trial_buffers(start_row_idx, n_lR_j);
        if (!gl20_ready) {
          global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = R_NegInf;
          global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
          return R_NegInf;
        }
        const int n_nodes = 20;
        const bool is_erlang2 = (race_ctx->kill_shape >= 2);
        const double inv_lam = 1.0 / lam;
        int active_count = 0;
        std::array<const double*, 64> par_rows{};
        for (int acc = 0; acc < n_lR_j; ++acc) {
          if (!isok_int_buffer[static_cast<size_t>(acc)]) continue;
          par_rows[static_cast<size_t>(active_count++)] =
            pars_rowmajor_buffer + static_cast<size_t>(acc) * n_par;
        }
        if (active_count == 0) {
          global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = R_NegInf;
          global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
          return R_NegInf;
        }
        double prob = 0.0;
        for (int j = 0; j < n_nodes; ++j) {
          const double tt = -gl_log1m_u[static_cast<size_t>(j)] * inv_lam;
          // Calculate S_race(tt)
          double logS_race_tt = 0.0;
          for (int acc = 0; acc < active_count; ++acc) {
            const double* par_row = par_rows[static_cast<size_t>(acc)];
            double Fk = cdf1(tt, par_row, model_context_for_funcs);
            Fk = clamp_cdf01_race(Fk);
            const double ll = safe_log1m_race(Fk);
            if (!R_FINITE(ll)) { logS_race_tt = R_NegInf; break; }
            logS_race_tt += ll;
          }
          if (!R_FINITE(logS_race_tt)) continue;
          const double sk_jac = is_erlang2
            ? gl_jac_erlang2[static_cast<size_t>(j)]
            : 1.0;
          prob += gl_w_half[static_cast<size_t>(j)] * std::exp(logS_race_tt) * sk_jac;
        }
        const double ans = (prob > 0.0) ? std::log(prob) : R_NegInf;
        global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = ans;
        global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
        return ans;
      }
      if (race_ctx && race_ctx->defective_upper_tail) {
        if (n_lR_j == 1 && isok[start_row_idx]) {
          static thread_local std::vector<double> par_buf_inf_tls2;
          par_buf_inf_tls2.assign(par_row_width, 0.0);
          double* par_buf_inf = par_buf_inf_tls2.data();
          for (int c = 0; c < n_par; ++c)
            par_buf_inf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + start_row_idx];
          double F_inf = cdf1(R_PosInf, par_buf_inf, model_context_for_funcs);
          F_inf = clamp_cdf01_race(F_inf);
          const double ans = safe_log1m_race(F_inf);
          global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = ans;
          global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
          return ans;
        }
        double log_p = 0.0;
        static thread_local std::vector<double> par_buf_inf_tls3;
        par_buf_inf_tls3.assign(par_row_width, 0.0);
        double* par_buf_inf = par_buf_inf_tls3.data();
        for (int k = 0; k < n_lR_j; ++k) {
          const int row = start_row_idx + k;
          if (!isok[row]) return R_NegInf;
          for (int c = 0; c < n_par; ++c)
            par_buf_inf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
          double Fk_inf = cdf1(R_PosInf, par_buf_inf, model_context_for_funcs);
          Fk_inf = clamp_cdf01_race(Fk_inf);
          const double ll = safe_log1m_race(Fk_inf);
          if (!emc2_isfinite(ll)) return R_NegInf;
          log_p += ll;
        }
        global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = log_p;
        global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
        return log_p;
      }
      global_log_surv_inf_by_trial[static_cast<size_t>(unique_trial_idx)] = R_NegInf;
      global_log_surv_inf_ready[static_cast<size_t>(unique_trial_idx)] = 1;
      return R_NegInf;
    }

    if (n_lR_j == 1) {
      if (!isok[start_row_idx]) return R_NegInf;
      static thread_local std::vector<double> par_buf_tls4;
      par_buf_tls4.assign(par_row_width, 0.0);
      double* par_buf = par_buf_tls4.data();
      for (int c = 0; c < n_par; ++c)
        par_buf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + start_row_idx];
      double Fk = cdf1(t, par_buf, model_context_for_funcs);
      Fk = clamp_cdf01_race(Fk);
      double logS = safe_log1m_race(Fk);
      if (!R_FINITE(logS)) return R_NegInf;
      if (race_ctx && race_ctx->is_global_kill && race_ctx->kill_active) {
        const double* mean_k_ptr = pars_cm_ptr + race_ctx->mean_k_index * n_trials;
        if (t > 0.0) {
          logS += erlang_log_surv(
            t,
            erlang_lambda_from_mean(mean_k_ptr[start_row_idx], race_ctx->kill_shape),
            race_ctx->kill_shape
          );
        }
      }
      return logS;
    }

    double logS = 0.0;
    static thread_local std::vector<double> par_buf_tls5;
    par_buf_tls5.assign(par_row_width, 0.0);
    double* par_buf = par_buf_tls5.data();
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start_row_idx + k;
      if (!isok[row]) return R_NegInf;
      for (int c = 0; c < n_par; ++c)
        par_buf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + row];
      double Fk = cdf1(t, par_buf, model_context_for_funcs);
      Fk = clamp_cdf01_race(Fk);
      const double ll = safe_log1m_race(Fk);
      if (!R_FINITE(ll)) return R_NegInf;
      logS += ll;
    }
    if (race_ctx && race_ctx->is_global_kill && race_ctx->kill_active) {
      const double* mean_k_ptr = pars_cm_ptr + race_ctx->mean_k_index * n_trials;
      if (t > 0.0) {
        logS += erlang_log_surv(
          t,
          erlang_lambda_from_mean(mean_k_ptr[start_row_idx], race_ctx->kill_shape),
          race_ctx->kill_shape
        );
      }
    }
    return logS;
  };

  const bool has_finite_batch = use_full_finite_batch || (finite_rt_unique_trial_indices.size() > 0);
  std::vector<int> winner_row_by_trial;
  std::vector<double> global_log_sk_by_trial;
  if (ctx && ctx->is_global_kill && ctx->kill_active) {
    winner_row_by_trial.assign(static_cast<size_t>(n_unique_trials), -1);
    global_log_sk_by_trial.assign(static_cast<size_t>(n_unique_trials), 0.0);
    const double* mean_k_ptr = pars_cm_ptr + static_cast<size_t>(ctx->mean_k_index) * n_trials;
    for (int j = 0; j < n_unique_trials; ++j) {
      const int start = j * n_lR;
      const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
      int idx_w = -1;
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (winner[row]) { idx_w = row; break; }
      }
      winner_row_by_trial[static_cast<size_t>(j)] = idx_w;
      if (idx_w >= 0) {
        const double tt = rts_dadm[idx_w];
        global_log_sk_by_trial[static_cast<size_t>(j)] =
          (tt > 0.0)
            ? erlang_log_surv(
                tt,
                erlang_lambda_from_mean(mean_k_ptr[idx_w], ctx->kill_shape),
                ctx->kill_shape
              )
            : min_ll;
      }
    }
  }

  if (has_finite_batch) {
    bool any_win = false;
    bool any_loss = false;
    // Local fallback storage (used when shared state is not available)
    std::vector<int> idx_win_int_local;
    std::vector<int> idx_loss_int_local;
    std::vector<int> isok_int_all_local;
    // Pointers resolved below (point into shared or local storage)
    const int* idx_win_ptr  = nullptr;
    const int* idx_loss_ptr = nullptr;
    int*       isok_ptr     = nullptr;

    const bool use_shared_bufs = use_shared &&
      static_cast<int>(shared->idx_win.size()) == n_trials;

    if (use_shared_bufs) {
      // Per-particle: only fill isok (winner/loser masks are data-fixed)
      for (int i = 0; i < n_trials; ++i)
        shared->ok_int_buf[static_cast<size_t>(i)] = isok[i] ? 1 : 0;
      idx_win_ptr  = shared->idx_win.data();
      idx_loss_ptr = shared->idx_loss.data();
      isok_ptr     = shared->ok_int_buf.data();
      any_win      = shared->any_win;
      any_loss     = shared->any_loss;
    } else {
      idx_win_int_local.assign(static_cast<size_t>(n_trials), 0);
      idx_loss_int_local.assign(static_cast<size_t>(n_trials), 0);
      isok_int_all_local.assign(static_cast<size_t>(n_trials), 0);
      for (int i = 0; i < n_trials; ++i) {
        isok_int_all_local[static_cast<size_t>(i)] = isok[i] ? 1 : 0;
        if (has_RACE_col && !RACE_mask[i]) continue;
        if (!use_full_finite_batch && !finite_rt_mask[i]) continue;
        if (winner[i]) {
          idx_win_int_local[static_cast<size_t>(i)] = 1;
          any_win = true;
        } else if (n_lR > 1) {
          idx_loss_int_local[static_cast<size_t>(i)] = 1;
          any_loss = true;
        }
      }
      idx_win_ptr  = idx_win_int_local.data();
      idx_loss_ptr = idx_loss_int_local.data();
      isok_ptr     = isok_int_all_local.data();
    }

    const double* rt_ptr = rts_dadm.begin();
    const double* const* pars_cm = cols_view;
    const bool dense_floor_raw_log_lik_prev = dense_ctx.floor_raw_log_lik;
    if (has_time) dense_ctx.floor_raw_log_lik = false;
    if (any_win) {
      model_dfun_raw(rt_ptr, pars_cm, n_trials,
                     idx_win_ptr, isok_ptr,
                     lds_ptr, min_ll, dense_ctx_ptr);
    }
    if (n_lR > 1 && any_loss) {
      model_pfun_raw(rt_ptr, pars_cm, n_trials,
                     idx_loss_ptr, isok_ptr,
                     lds_ptr, min_ll, dense_ctx_ptr);
    }

    std::vector<int> idx_time_only_local;
    const int* idx_time_only_ptr = nullptr;
    std::vector<double> alt_lds_local;
    double* alt_lds_ptr = nullptr;
    std::vector<int> n_resp_local;
    const int* n_resp_ptr = nullptr;

    const bool needs_n_resp = has_time;
    if (needs_n_resp) {
      if (use_shared && !shared->n_resp.empty()) {
        n_resp_ptr = shared->n_resp.data();
      } else {
        n_resp_local.assign(static_cast<size_t>(n_unique_trials), 0);
        const SEXP lR_sexp = dadm["lR"];
        Rcpp::IntegerVector lR_code_vec_int(lR_sexp);
        for (int j = 0; j < n_unique_trials; ++j) {
          const int start = j * n_lR;
          const int n_lR_curr = has_RACE_col ? RACE[start] : n_lR;
          n_resp_local[static_cast<size_t>(j)] = count_resp_accumulators(
              lR_code_vec_int.begin(), start, n_lR_curr, time_code, nogo_code);
        }
        n_resp_ptr = n_resp_local.data();
      }
    }

    if (has_time) {
      if (use_shared && !shared->idx_time_only.empty()) {
        idx_time_only_ptr = shared->idx_time_only.data();
        alt_lds_ptr = shared->alt_res_buf.data();
        // n_resp_ptr already set above
      } else {
        idx_time_only_local.assign(n_trials, 0);
        alt_lds_local.assign(n_trials, min_ll);
        const SEXP lR_sexp = dadm["lR"];
        Rcpp::IntegerVector lR_code_vec_int(lR_sexp);
        for (int j = 0; j < n_trials; ++j) {
          const bool active = !has_RACE_col || RACE_mask[j];
          if (active && lR_code_vec_int[j] == time_code) idx_time_only_local[j] = 1;
        }
        // n_resp already filled in the needs_n_resp block above
        idx_time_only_ptr = idx_time_only_local.data();
        alt_lds_ptr = alt_lds_local.data();
      }
      model_dfun_raw(rt_ptr, pars_cm, n_trials,
                     idx_time_only_ptr, isok_ptr,
                     alt_lds_ptr, min_ll, dense_ctx_ptr);
      model_pfun_raw(rt_ptr, pars_cm, n_trials,
                     idx_win_ptr, isok_ptr,
                     alt_lds_ptr, min_ll, dense_ctx_ptr);
    }
    dense_ctx.floor_raw_log_lik = dense_floor_raw_log_lik_prev;

    // --- Pre-compute truncation normalisers in batch when possible ---
    // When may_need_ct AND logS_at_t is available, we can replace per-trial
    // fill_trial_buffers + get_trunc_normaliser_rowmajor_cpp with a single
    // column-major pass: eliminates n_lR*n_par row copies per truncated trial.
    //
    // Fast-batch conditions: LT == 0 everywhere (logS_LT = 0 analytically),
    // UT is uniform across all truncated trials, and logS_at_t_ptr is provided.
    // Trials that fail the fast path fall back to the per-trial scalar route.
    bool batch_trunc_done = false;
    std::vector<double> logZ_batch;  // length n_unique_trials, filled when batch_trunc_done
    const int n_finite_unique = use_full_finite_batch
        ? n_unique_trials
        : static_cast<int>(finite_rt_unique_trial_indices.size());

    RaceEndpointGroupCache local_endpoint_cache;
    RaceEndpointGroupCache* endpoint_cache =
      use_shared ? &shared->race_endpoint_cache : &local_endpoint_cache;
    endpoint_cache->new_particle();
    auto batch_logS_at_t = [&](double t, const std::vector<int>& include_mask,
                               std::vector<double>& out) {
      // RACE-column trials have varying active accumulator counts and remain on
      // the existing scalar-safe path below.  Fixed-width race trials can be
      // compacted by complete parameter-block key before calling the model
      // callback, so repeated design cells pay one endpoint calculation.
      if (has_RACE_col || logS_at_t == nullptr) return;
      const int skip_pc = use_shared ? shared->pc_col : -1;
      const int skip_pg = use_shared ? shared->pg_col : -1;
      race_endpoint_prepare_groups(
          *endpoint_cache, cols_view, n_unique_trials, n_lR, n_par,
          include_mask.data(), isok_ptr, skip_pc, skip_pg);
      const int n_groups = static_cast<int>(endpoint_cache->representative.size());
      out.assign(static_cast<size_t>(n_unique_trials), R_NegInf);
      if (n_groups == 0) return;
      if (n_groups == endpoint_cache->n_included) {
        logS_at_t(t, cols_view, n_trials, n_lR, n_par,
                  include_mask.data(), n_unique_trials, isok_ptr,
                  model_context_for_funcs, out.data());
        return;
      }
      std::vector<double> compact_out(static_cast<size_t>(n_groups), R_NegInf);
      logS_at_t(t, endpoint_cache->compact_col_ptrs.data(),
                n_groups * n_lR, n_lR, n_par,
                endpoint_cache->compact_mask.data(), n_groups,
                endpoint_cache->compact_isok.data(),
                model_context_for_funcs, compact_out.data());
      for (int j = 0; j < n_unique_trials; ++j) {
        const int group = endpoint_cache->group_id[static_cast<size_t>(j)];
        if (group >= 0) out[static_cast<size_t>(j)] = compact_out[static_cast<size_t>(group)];
      }
    };

    if (apply_truncation_correction && may_need_ct && logS_at_t != nullptr && !has_RACE_col) {
      // Pass 1: scan for truncated trials; check uniformity of LT and UT separately.
      // We can batch the normaliser whenever all truncated finite-RT trials share the
      // same LT value AND the same UT value (each may be 0/Inf trivially).
      double uniform_LT = -1.0;   // -1 = not yet observed (LT >= 0 always)
      double uniform_UT = -1.0;   // -1 = not yet observed
      bool uniform_LT_ok = true;
      bool uniform_UT_ok = true;
      bool any_trunc = false;
      std::vector<int> trunc_mask(static_cast<size_t>(n_unique_trials), 0);

      for (int i = 0; i < n_finite_unique; ++i) {
        const int j = use_full_finite_batch
            ? i : finite_rt_unique_trial_indices[static_cast<size_t>(i)];
        const int start = j * n_lR;
        const double LTj = LT[start];
        const double UTj = UT[start];
        if (LTj != 0.0 || UTj != R_PosInf) {
          trunc_mask[static_cast<size_t>(j)] = 1;
          any_trunc = true;
          if (uniform_LT < 0.0) uniform_LT = LTj;
          else if (LTj != uniform_LT) uniform_LT_ok = false;
          if (uniform_UT < 0.0) uniform_UT = UTj;
          else if (UTj != uniform_UT) uniform_UT_ok = false;
        }
      }
      // Resolve sentinels: if no truncated trial was seen, treat as trivial.
      if (uniform_LT < 0.0) uniform_LT = 0.0;
      if (uniform_UT < 0.0) uniform_UT = R_PosInf;

      // Global kill is excluded: logS_at_t zeroes the kill rate for the racers
      // (utils.h), so the batched survivor is the unkilled one.  The per-trial
      // scalar path routes through get_trunc_normaliser_rowmajor_cpp, which
      // switches the kill back on for the single-accumulator case and rejects
      // the multi-accumulator one.
      // A defective model under FINITE upper truncation needs the retained
      // T = +Inf atom added to the normaliser (see
      // get_trunc_normaliser_rowmajor_cpp).  Rather than ask every model's
      // logS_at_t adapter to be correct at t = +Inf, these trials are sent to
      // the per-trial scalar route, whose log_survivor_rowmajor() has an
      // explicit defective +Inf branch.  UT == Inf is unaffected and stays
      // batched: there the atom is already inside S(LT).
      const bool defective_finite_UT =
          defective_upper_tail && uniform_UT != R_PosInf;
      if (any_trunc && uniform_LT_ok && uniform_UT_ok && !has_RACE_col &&
          !global_omission_active && !defective_finite_UT) {
        // Pass 2: batch-compute logZ = log_diff_exp(logS(LT), logS(UT)) for all
        // truncated trials simultaneously.
        // RACE models are excluded: lba/rdm/lnr_logS_at_t always loops over the
        // global n_lR, but RACE trials may have fewer active accumulators (inactive
        // rows are NA-masked).  The per-trial scalar path handles n_lR_j correctly.
        //
        // Trivial endpoints:
        //   LT == 0   → logS(LT) = 0  (S(0)=1 for any proper distribution)
        //   UT == Inf → logS(UT) = -Inf (S(∞)=0 for any proper distribution)
        // For non-trivial endpoints we call logS_at_t in batch.
        //
        // logS_LT: 0 when LT==0 (trivial), otherwise computed in batch.
        std::vector<double> logS_LT_vec(static_cast<size_t>(n_unique_trials), 0.0);
        if (uniform_LT != 0.0) {
          batch_logS_at_t(uniform_LT, trunc_mask, logS_LT_vec);
        }

        // logS_UT: -Inf when UT==Inf (trivial for proper distributions), else computed.
        std::vector<double> logS_UT_vec(static_cast<size_t>(n_unique_trials), R_NegInf);
        if (uniform_UT != R_PosInf) {
          batch_logS_at_t(uniform_UT, trunc_mask, logS_UT_vec);
        }

        // logZ = log(S(LT) - S(UT)) = log_diff_exp(logS_LT, logS_UT)
        // log_diff_exp handles R_NegInf inputs correctly:
        //   log_diff_exp(0, -Inf)  = 0   (LT=0,  UT=Inf → P=1)
        //   log_diff_exp(x, -Inf)  = x   (UT=Inf  → P=S(LT))
        //   log_diff_exp(0, y)     = log(1-exp(y)) (LT=0 → original UT-only case)
        // NaN from logS_at_t (bad params) propagates through log_diff_exp as NA_REAL,
        // which is caught by emc2_isnan and sent to the per-trial fallback.
        logZ_batch.assign(static_cast<size_t>(n_unique_trials), 0.0);
        static const double kLogProbEps = std::log(std::numeric_limits<double>::epsilon());
        for (int j = 0; j < n_unique_trials; ++j) {
          if (!trunc_mask[static_cast<size_t>(j)]) continue;
          const double logS_LT_j = logS_LT_vec[static_cast<size_t>(j)];
          const double logS_UT_j = logS_UT_vec[static_cast<size_t>(j)];
          const double logP = log_diff_exp(logS_LT_j, logS_UT_j);
          if (R_FINITE(logP) && logP > kLogProbEps) {
            logZ_batch[static_cast<size_t>(j)] = logP;
          } else {
            // NaN (bad params), catastrophic cancellation, or P≈0:
            // fall back to the per-trial scalar/GSL path.
            logZ_batch[static_cast<size_t>(j)] = std::numeric_limits<double>::quiet_NaN();
          }
        }
        batch_trunc_done = true;
      }
    }

    // Apply truncation correction and calculate log-likelihood for each trial in the batch
    for (int i = 0; i < n_finite_unique; ++i) {
      // When all unique trials are finite, iterate in order 0..n_unique_trials-1.
      // Otherwise, iterate over the precomputed subset of finite-RT unique trials.
      const int unique_trial_idx = use_full_finite_batch
          ? i : finite_rt_unique_trial_indices[static_cast<size_t>(i)];
      const int start_row_idx = unique_trial_idx * n_lR;
      n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
      const double rt_j = rts_dadm[start_row_idx];
      log_Z_this = 0.0;

      if (apply_truncation_correction && may_need_ct) {
        const double LTj = LT[start_row_idx];
        const double UTj = UT[start_row_idx];
        if (LTj != 0.0 || UTj != R_PosInf) { // truncation active
          bool need_scalar = true;
          if (batch_trunc_done) {
            const double bz = logZ_batch[static_cast<size_t>(unique_trial_idx)];
            if (!emc2_isnan(bz)) {       // NaN flags a fallback case
              log_Z_this = bz;
              need_scalar = false;
            }
          }
          if (need_scalar) {
            fill_trial_buffers(start_row_idx, n_lR_j);
            log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer,
                                                           isok_int_buffer,
                                                           pdf1, cdf1,
                                                           LTj, UTj,
                                                           n_lR_j, n_par,
                                                           gsl_ctl,
                                                           model_context_for_funcs,
                                                           workspace);
          }
          double current_trial_ll_sum = 0.0;
          bool hit_min_ll = false;
          for (int k = 0; k < n_lR_j; ++k) {
            double v = lds_ptr[start_row_idx + k];
            if (!has_time && v <= min_ll) hit_min_ll = true;
            current_trial_ll_sum += v;
          }
          if (has_time && n_resp_ptr[unique_trial_idx] > 0) {
            int idx_W = -1, idx_T = -1;
            for (int k = 0; k < n_lR_j; ++k) {
              if (idx_win_ptr[start_row_idx + k]) idx_W = k;
              if (idx_time_only_ptr[start_row_idx + k]) idx_T = k;
            }
            if (idx_W != -1 && idx_T != -1) {
              double s_T = 0.0;
              for (int k = 0; k < n_lR_j; ++k) {
                if (k == idx_W) s_T += alt_lds_ptr[start_row_idx + idx_W];
                else if (k == idx_T) s_T += alt_lds_ptr[start_row_idx + idx_T];
                else s_T += lds_ptr[start_row_idx + k];
              }
              current_trial_ll_sum = log_sum_exp(current_trial_ll_sum, s_T - std::log(n_resp_ptr[unique_trial_idx]));
            }
          } else if (ctx && ctx->is_global_kill && ctx->kill_active) {
            current_trial_ll_sum += global_log_sk_by_trial[static_cast<size_t>(unique_trial_idx)];
          }
          if (NumericVector::is_na(log_Z_this) || !R_FINITE(log_Z_this)) {
            ll_unique[unique_trial_idx] = min_ll;
          } else if (hit_min_ll) {
            ll_unique[unique_trial_idx] = min_ll;
          } else {
            ll_unique[unique_trial_idx] = std::max(min_ll, current_trial_ll_sum - log_Z_this);
          }
          continue;
        }
      }

      // No truncation (or all_finite_trials guarantees none): plain sum.
      double current_trial_ll_sum = 0.0;
      for (int k = 0; k < n_lR_j; ++k)
        current_trial_ll_sum += lds_ptr[start_row_idx + k];
      if (has_time && n_resp_ptr[unique_trial_idx] > 0) {
        int idx_W = -1, idx_T = -1;
        for (int k = 0; k < n_lR_j; ++k) {
          if (idx_win_ptr[start_row_idx + k]) idx_W = k;
          if (idx_time_only_ptr[start_row_idx + k]) idx_T = k;
        }
        if (idx_W != -1 && idx_T != -1) {
          double s_T = 0.0;
          for (int k = 0; k < n_lR_j; ++k) {
            if (k == idx_W) s_T += alt_lds_ptr[start_row_idx + idx_W];
            else if (k == idx_T) s_T += alt_lds_ptr[start_row_idx + idx_T];
            else s_T += lds_ptr[start_row_idx + k];
          }
          current_trial_ll_sum = log_sum_exp(current_trial_ll_sum, s_T - std::log(n_resp_ptr[unique_trial_idx]));
        }
      } else if (ctx && ctx->is_global_kill && ctx->kill_active) {
        current_trial_ll_sum += global_log_sk_by_trial[static_cast<size_t>(unique_trial_idx)];
      }
      ll_unique[unique_trial_idx] = std::max(min_ll, current_trial_ll_sum);
    }

  }
  // --- Process other trials (Infinite RTs, NA RTs, or finite RTs outside truncation) ---
  // These trials require individual processing, often involving numerical integration for censored intervals.
  double current_ll_val;
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
  for (size_t i = 0; i < other_unique_trial_indices.size(); ++i) {
    int unique_trial_idx = other_unique_trial_indices[i];
    int start_row_idx = unique_trial_idx*n_lR;
    log_Z_this = 0.0;
    n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
    const double rt_j = rts_dadm[start_row_idx];
    const int R_j_idx = R_idxs_dadm[start_row_idx];
    const double LTj = LT[start_row_idx];
    const double UTj = UT[start_row_idx];
    const double LCj = LC[start_row_idx];
    const double UCj = UC[start_row_idx];
    const bool has_trunc = apply_truncation_correction &&
      (LTj != 0.0 || UTj != R_PosInf);
    const bool trial_has_active_nogo =
      (static_cast<size_t>(unique_trial_idx) < active_nogo_trial_mask.size()) &&
      (active_nogo_trial_mask[static_cast<size_t>(unique_trial_idx)] != 0);
    
    const bool needs_model = (rt_j == R_NegInf) || (rt_j == R_PosInf) || Rcpp::NumericVector::is_na(rt_j) ||
      (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER);
    if (!needs_model) {
      ll_unique[unique_trial_idx] = min_ll;
      continue;
    }
    
    // Lazy row-major buffer fill — only when GSL integration or row-major
    // helper functions are actually needed.  Pure survivor / CDF calls use
    // log_surv_cm which reads column-major directly.
    bool buffers_filled = false;
    auto ensure_buffers = [&]() {
      if (!buffers_filled) {
        fill_trial_buffers(start_row_idx, n_lR_j);
        buffers_filled = true;
      }
    };

    auto integrate_interval = [&](int k_winner_1based, double low, double upp) -> double {
      ensure_buffers();
      gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
      return integrate_for_kth_winner_rowmajor_cpp(k_winner_1based,
                                                   pars_rowmajor_buffer,
                                                   isok_int_buffer,
                                                   low, upp, pdf1, cdf1,
                                                   n_lR_j, n_par, gsl_ctl,
                                                   model_context_for_funcs, w);
    };

    int idx_T_1based = -1;
    int n_resp_j = 0;
    if (has_time) {
      const SEXP lR_sexp = dadm["lR"];
      Rcpp::IntegerVector lR_code_vec_int(lR_sexp);
      for (int k = 0; k < n_lR_j; ++k) {
        int code = lR_code_vec_int[start_row_idx + k];
        if (code == time_code) idx_T_1based = k + 1;
        else if (code != nogo_code) n_resp_j++;
      }
    }

    auto integrate_timed = [&](int k_winner_1based, double low, double upp) -> double {
      double ll = integrate_interval(k_winner_1based, low, upp);
      if (idx_T_1based != -1 && k_winner_1based != idx_T_1based && n_resp_j > 0) {
        double ll_T = integrate_interval(idx_T_1based, low, upp);
        ll = log_sum_exp(ll, ll_T - std::log(n_resp_j));
      }
      return ll;
    };

    current_ll_val = R_NegInf;
    if (rt_j == R_NegInf) {
      lower_for_trial = LTj;
      upper_for_trial = LCj;
      if (R_j_idx == NA_INTEGER) {
        if (trial_has_active_nogo) {
          // For go/no-go, rt == -Inf means an observed response below LC.
          // The no-go accumulator cannot be the winner for this event, so we
          // skip that winner index when summing winner-specific densities.
          int k_nogo = -1; // 1-based
          int n_true = 0;
          for (int k = 0; k < n_lR_j; ++k) {
            if (winner[start_row_idx + k]) { n_true++; k_nogo = k + 1; }
          }
          if (n_true != 1) k_nogo = -1; // fallback: keep normal behavior
          for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
            if (k_win == k_nogo) continue;
            current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
          }
        } else {
          // Left-censoring with unknown winner: S(LT) - S(LC), S(t)=prod_k(1-F_k(t)).
          const double logP = log_diff_exp(log_surv_cm(lower_for_trial, start_row_idx, n_lR_j),
                                           log_surv_cm(upper_for_trial, start_row_idx, n_lR_j));
          if (R_FINITE(logP) && logP > log_prob_eps) {
            current_ll_val = logP;
          } else {
            for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
              current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
            }
          }
        }
      } else {
        current_ll_val = integrate_timed(R_j_idx, lower_for_trial, upper_for_trial);
      }
    } else if (rt_j == R_PosInf) {
      if (trial_has_active_nogo) {
        lower_for_trial = LTj;
        upper_for_trial = UCj;
        int k_nogo = -1; // 1-based
        int n_true = 0;
        for (int k = 0; k < n_lR_j; ++k) {
          if (winner[start_row_idx + k]) { n_true++; k_nogo = k + 1; }
        }
        // Robust fallback: for a single-accumulator go/no-go trial with no
        // observed response, winner may be all FALSE in dadm. In that case the
        // only valid no-go accumulator is index 1.
        if (n_true == 0 && n_lR_j == 1 && R_j_idx == NA_INTEGER) {
          k_nogo = 1;
          n_true = 1;
        }
        if (n_true != 1) Rcpp::stop("No winner identified in go/no-go withheld response");
        const double logA = integrate_interval(k_nogo, lower_for_trial, upper_for_trial);
        const double logB = log_surv_cm(upper_for_trial, start_row_idx, n_lR_j); // incl. LBAIO mass at +Inf
        current_ll_val = log_sum_exp(logA, logB);
      } else {
        lower_for_trial = UCj;
        upper_for_trial = UTj;
        if (R_j_idx == NA_INTEGER) {
          if (global_omission_active && n_lR_j == 1 && !trial_has_active_nogo) {
            // Single-accumulator analytic branch.
            static thread_local std::vector<double> par_buf_tls6;
            par_buf_tls6.assign(par_row_width, 0.0);
            double* par_buf = par_buf_tls6.data();
            for (int c = 0; c < n_par; ++c)
              par_buf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + start_row_idx];
            const auto logS_single = [&](double t) -> double {
              ContextForRaceModels* ctx_loc = static_cast<ContextForRaceModels*>(model_context_for_funcs);
              const bool old_apply = ctx_loc->apply_lk_to_racers;
              ctx_loc->apply_lk_to_racers = true;
              double F = cdf1(t, par_buf, ctx_loc);
              ctx_loc->apply_lk_to_racers = old_apply;
              F = clamp_cdf01_race(F);
              return safe_log1m_race(F);
            };
            const double logP = (!defective_upper_tail)
              ? log_diff_exp(logS_single(lower_for_trial), logS_single(upper_for_trial))
              : logS_single(lower_for_trial);
            if (R_FINITE(logP) && logP > log_prob_eps) {
              current_ll_val = logP;
              goto apply_trial_trunc;
            }
          }
          // Fast exit: when no upper censoring (UCj=Inf) and model has no
          // defective upper tail, P(RT=Inf)=0. Avoid log_diff_exp(-Inf,-Inf)
          // and degenerate integration bounds [Inf,Inf].
          if (!defective_upper_tail && !R_FINITE(lower_for_trial)) {
            // current_ll_val stays at min_ll
          } else {
            const double logP = (!defective_upper_tail)
              ? log_diff_exp(log_surv_cm(lower_for_trial, start_row_idx, n_lR_j),
                             log_surv_cm(upper_for_trial, start_row_idx, n_lR_j))
              : log_surv_cm(lower_for_trial, start_row_idx, n_lR_j); // include intrinsic never-finish mass
            if (R_FINITE(logP) && logP > log_prob_eps) {
              current_ll_val = logP;
            } else { // numerical integration fallback
              for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
                current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
              }
              if (defective_upper_tail) { // include defective upper-tail mass
                const double log_p_I = log_surv_cm(R_PosInf, start_row_idx, n_lR_j);
                current_ll_val = log_sum_exp(current_ll_val, log_p_I);
              }
            }
          }
        } else {
          current_ll_val = integrate_timed(R_j_idx, lower_for_trial, upper_for_trial);
          if (defective_upper_tail && n_lR_j == 1) {
            current_ll_val = log_sum_exp(current_ll_val,
                                         log_surv_cm(R_PosInf, start_row_idx, n_lR_j));
          }
        }
      }
    } else if (Rcpp::NumericVector::is_na(rt_j)) {
      const double lower1 = LTj;
      const double upper1 = LCj;
      const double lower2 = UCj;
      const double upper2 = UTj;
      if (R_j_idx != NA_INTEGER) {
        current_ll_val = log_sum_exp(integrate_timed(R_j_idx, lower1, upper1),
                                     integrate_timed(R_j_idx, lower2, upper2));
      } else {
        if (global_omission_active && n_lR_j == 1 && !trial_has_active_nogo) {
          // Single-accumulator analytic branch for missing RT interval union.
          static thread_local std::vector<double> par_buf_tls7;
          par_buf_tls7.assign(par_row_width, 0.0);
          double* par_buf = par_buf_tls7.data();
          for (int c = 0; c < n_par; ++c)
            par_buf[c] = pars_cm_ptr[static_cast<size_t>(c) * n_trials + start_row_idx];
          const auto logS_single = [&](double t) -> double {
            ContextForRaceModels* ctx_loc = static_cast<ContextForRaceModels*>(model_context_for_funcs);
            const bool old_apply = ctx_loc->apply_lk_to_racers;
            ctx_loc->apply_lk_to_racers = true;
            double F = cdf1(t, par_buf, ctx_loc);
            ctx_loc->apply_lk_to_racers = old_apply;
            F = clamp_cdf01_race(F);
            return safe_log1m_race(F);
          };
          const double logP1 = log_diff_exp(logS_single(lower1), logS_single(upper1));
          const double logP2 = (!defective_upper_tail)
            ? log_diff_exp(logS_single(lower2), logS_single(upper2))
            : logS_single(lower2);
          const double logPsum = log_sum_exp(logP1, logP2);
          if (R_FINITE(logPsum) && logPsum > log_prob_eps) {
            current_ll_val = logPsum;
            goto apply_trial_trunc;
          }
        }
        const double logP1 = log_diff_exp(log_surv_cm(lower1, start_row_idx, n_lR_j),
                                          log_surv_cm(upper1, start_row_idx, n_lR_j));
        const double logP2 = (!defective_upper_tail)
          ? log_diff_exp(log_surv_cm(lower2, start_row_idx, n_lR_j),
                         log_surv_cm(upper2, start_row_idx, n_lR_j))
          : log_surv_cm(lower2, start_row_idx, n_lR_j);
        const double logPsum = log_sum_exp(logP1, logP2);
        if (R_FINITE(logPsum) && logPsum > log_prob_eps) {
          current_ll_val = logPsum;
        } else {
          for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
            const double ll_L_k = integrate_interval(k_win, lower1, upper1);
            const double ll_U_k = integrate_interval(k_win, lower2, upper2);
            current_ll_val = log_sum_exp(current_ll_val, log_sum_exp(ll_L_k, ll_U_k));
          }
          if (defective_upper_tail) {
            const double log_p_I = log_surv_cm(R_PosInf, start_row_idx, n_lR_j);
            current_ll_val = log_sum_exp(current_ll_val, log_p_I);
          }
        }
      }
    } else if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER) {
      ensure_buffers();
      current_ll_val = log_min_density_rowmajor(rt_j,
                                                pars_rowmajor_buffer,
                                                isok_int_buffer,
                                                n_lR_j, n_par, pdf1, cdf1,
                                                model_context_for_funcs,
                                                logS_k_buffer);
    }

apply_trial_trunc:
    if (current_ll_val > min_ll && has_trunc) {
      ensure_buffers();
      log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer,
                                                     isok_int_buffer,
                                                     pdf1, cdf1,
                                                     LTj, UTj, n_lR_j, n_par,
                                                     gsl_ctl,
                                                     model_context_for_funcs,
                                                     workspace);
      if (!R_FINITE(log_Z_this)) current_ll_val = min_ll;
      else current_ll_val -= log_Z_this;
    }
    
    ll_unique[unique_trial_idx] = std::max(min_ll, current_ll_val);
  }
  
  
  // --- Summation of log-likelihoods for all unique trials ---
  // pC modification (if any) is kept in a separate pass so the final summation
  // loop is a pure reduction — allowing #pragma omp simd to vectorize it.
  //
  // A contaminant is an intrinsic omission that only manifests as a
  // never-responded (rt == +Inf) trial, so ONLY those trials pick up the
  // log(pC) mass; every trial is down-weighted by log(1 - pC). This must match
  // the R reference (likelihood.R:372-374), which adds the omission mass solely
  // when rt is +Inf — a left-censored (-Inf) or missing (NA) rt does NOT. One
  // lambda applied in both the expand and compressed branches so they cannot
  // drift apart again.
  //
  // pGuess rides in the same lambda: a uniform density over the guess window,
  // for finite-rt trials only.  The one branch that reaches a finite rt with an
  // UNKNOWN R (log_min_density_rowmajor, above) takes the window density
  // without the n_resp division, which mix_contaminants_rt() handles.
  auto apply_pC = [&](int j) {
    const double rt_j = rts_dadm[j * n_lR];
    const double pC = use_pC ? pars_cm_ptr[static_cast<size_t>(pc_col) * n_trials + j * n_lR] : 0.0;
    const double pG = use_pG ? pars_cm_ptr[static_cast<size_t>(pg_col) * n_trials + j * n_lR] : 0.0;
    ll_unique[j] = mix_contaminants_rt(ll_unique[j], pC, pG,
                                       guess, rt_j,
                                       R_idxs_dadm[j * n_lR] != NA_INTEGER);
  };
  const bool use_mix = use_pC || use_pG;
  double total_ll = 0;
  if (expand.length() > 0) { // non-compressed dadm: sum via expand index vector
    if (use_mix) {
      for (int j = 0; j < n_unique_trials; ++j) apply_pC(j);
    }
    const double* ll_ptr = ll_unique.data();
    const int* ex_ptr = expand.begin();
    const int n_exp = expand.length();
    if (trial_ll_out != nullptr) {
      for (int i = 0; i < n_exp; ++i) {
        const double v = ll_ptr[ex_ptr[i] - 1];
        (*trial_ll_out)[i] = v;
        total_ll += v;
      }
    } else {
      // Gather + reduce: indirect indexing prevents full SIMD gather, but the
      // reduction itself still benefits from the pragma.
      #pragma omp simd reduction(+:total_ll)
      for (int i = 0; i < n_exp; ++i) {
        total_ll += ll_ptr[ex_ptr[i] - 1];
      }
    }
  } else { // compressed dadm: each unique trial counted once
    if (use_mix) {
      for (int j = 0; j < n_unique_trials; ++j) apply_pC(j);
    }
    const double* ll_ptr = ll_unique.data();
    if (trial_ll_out != nullptr) {
      for (int j = 0; j < n_unique_trials; ++j) {
        const double v = ll_ptr[j];
        (*trial_ll_out)[j] = v;
        total_ll += v;
      }
    } else {
      // Pure sequential reduction — fully vectorizable.
      #pragma omp simd reduction(+:total_ll)
      for (int j = 0; j < n_unique_trials; ++j) {
        total_ll += ll_ptr[j];
      }
    }
  }
  
  return total_ll;
}


// Test/benchmark observability accessors for the correlated BAwL route
// counters (see src/bawl_corr_counters.h).  Counting is active only while
// EMC2_BAWLCORR_COUNTERS is set; the accessors themselves always work.
// [[Rcpp::export]]
Rcpp::List bawl_corr_counter_values() {
  const BAwLCorrCounters& c = bawl_corr_counters();
  return Rcpp::List::create(
      Rcpp::Named("ordinary_zero_rho_trials") = (double)c.ordinary_zero_rho_trials,
      Rcpp::Named("ordinary_single_loaded_trials") = (double)c.ordinary_single_loaded_trials,
      Rcpp::Named("exact_pair_trials") = (double)c.exact_pair_trials,
      Rcpp::Named("exact_pair_pair_winner_trials") = (double)c.exact_pair_pair_winner_trials,
      Rcpp::Named("exact_pair_independent_winner_trials") = (double)c.exact_pair_independent_winner_trials,
      Rcpp::Named("exact_pair_point_start_trials") = (double)c.exact_pair_point_start_trials,
      Rcpp::Named("numeric_pair_trials") = (double)c.numeric_pair_trials,
      Rcpp::Named("gh_no_clock_trials") = (double)c.gh_no_clock_trials,
      Rcpp::Named("gh_generic_clock_trials") = (double)c.gh_generic_clock_trials,
      Rcpp::Named("unstable_pair_floored_trials") = (double)c.unstable_pair_floored_trials,
      Rcpp::Named("loaded_dimension_0") = (double)c.loaded_dimension_0,
      Rcpp::Named("loaded_dimension_1") = (double)c.loaded_dimension_1,
      Rcpp::Named("loaded_dimension_2") = (double)c.loaded_dimension_2,
      Rcpp::Named("loaded_dimension_3plus") = (double)c.loaded_dimension_3plus,
      Rcpp::Named("prepared_rows") = (double)c.prepared_rows,
      Rcpp::Named("fused_node_evaluations") = (double)c.fused_node_evaluations,
      Rcpp::Named("bvn_corner_evaluations") = (double)c.bvn_corner_evaluations,
      Rcpp::Named("analytic_center_eligible_trials") = (double)c.analytic_center_eligible_trials,
      Rcpp::Named("analytic_center_success_trials") = (double)c.analytic_center_success_trials,
      Rcpp::Named("survivor_scan_trials") = (double)c.survivor_scan_trials,
      Rcpp::Named("scan_refinement_trials") = (double)c.scan_refinement_trials,
      Rcpp::Named("scan_node_evaluations") = (double)c.scan_node_evaluations,
      Rcpp::Named("fine_node_evaluations") = (double)c.fine_node_evaluations,
      Rcpp::Named("den_quadrature_trials") = (double)c.den_quadrature_trials,
      Rcpp::Named("den_quadrature_node_evaluations") = (double)c.den_quadrature_node_evaluations);
}

// [[Rcpp::export]]
void bawl_corr_counters_reset() {
  bawl_corr_counters().reset();
}


// Numerical probes for the exact probability layer.  They are intentionally
// small and allocation-free on the C++ side; testthat uses them to compare the
// rectangle derivatives and pair component against an independent oracle.
// [[Rcpp::export]]
Rcpp::List bawl_corr_bvn_rect_probe(double mu1, double sd1, double mu2,
                                    double sd2, double rho, double lo1,
                                    double hi1, double lo2, double hi2) {
  const BvnRectMoments m = bawl_corr_bvn_rect_moments(
      mu1, sd1, mu2, sd2, rho, lo1, hi1, lo2, hi2);
  return Rcpp::List::create(
      Rcpp::Named("p") = m.p, Rcpp::Named("m1") = m.m1,
      Rcpp::Named("m2") = m.m2, Rcpp::Named("m12") = m.m12,
      Rcpp::Named("status") = static_cast<int>(m.status));
}

// [[Rcpp::export]]
Rcpp::List bawl_corr_pair_probe(double t, double t01, double A1, double B1,
                                double k1, double v1, double sv1,
                                double t02, double A2, double B2, double k2,
                                double v2, double sv2, double rho,
                                bool posdrift = false, bool numeric = false) {
  const BAwLTimeGeometry g1 = bawl_time_geometry(t, t01, A1, B1 + A1, k1);
  const BAwLTimeGeometry g2 = bawl_time_geometry(t, t02, A2, B2 + A2, k2);
  const double D = bawl_corr_pair_positive_normalizer(v1, sv1, v2, sv2,
                                                       rho, posdrift);
  const BAwLCorrPairResult s = numeric
    ? bawl_corr_pair_survival_numeric(g1, g2, v1, sv1, v2, sv2, rho,
                                      posdrift, D)
    : bawl_corr_pair_survival_exact(g1, g2, v1, sv1, v2, sv2, rho,
                                    posdrift, D);
  const BAwLCorrPairResult c1 = numeric
    ? bawl_corr_pair_cause_numeric(g1, g2, v1, sv1, v2, sv2, rho,
                                   posdrift, D)
    : bawl_corr_pair_cause_exact(g1, g2, v1, sv1, v2, sv2, rho,
                                 posdrift, D);
  const BAwLCorrPairResult c2 = numeric
    ? bawl_corr_pair_cause_numeric(g2, g1, v2, sv2, v1, sv1, rho,
                                   posdrift, D)
    : bawl_corr_pair_cause_exact(g2, g1, v2, sv2, v1, sv1, rho,
                                 posdrift, D);
  return Rcpp::List::create(
      Rcpp::Named("normalizer") = D,
      Rcpp::Named("survival") = s.value,
      Rcpp::Named("cause1") = c1.value,
      Rcpp::Named("cause2") = c2.value,
      Rcpp::Named("survival_status") = static_cast<int>(s.status),
      Rcpp::Named("cause1_status") = static_cast<int>(c1.status),
      Rcpp::Named("cause2_status") = static_cast<int>(c2.status));
}

// Test probes for the shared correlated-BAwL geometry (T1).  `B` is the
// relative threshold; the geometry consumes b = B + A like the raw kernels.
// [[Rcpp::export]]
Rcpp::List bawl_time_geometry_probe(double t, double t0, double A, double B,
                                    double k) {
  const BAwLTimeGeometry g = bawl_time_geometry(t, t0, A, B + A, k);
  return Rcpp::List::create(
      Rcpp::Named("status") = static_cast<int>(g.status),
      Rcpp::Named("tau") = g.tau,
      Rcpp::Named("b") = g.b,
      Rcpp::Named("E") = g.E,
      Rcpp::Named("G") = g.G,
      Rcpp::Named("C1") = g.C1,
      Rcpp::Named("C2") = g.C2,
      Rcpp::Named("L") = g.L,
      Rcpp::Named("U") = g.U,
      Rcpp::Named("alpha") = g.alpha,
      Rcpp::Named("beta") = g.beta,
      Rcpp::Named("gamma0") = g.gamma0,
      Rcpp::Named("gamma1") = g.gamma1,
      Rcpp::Named("dv_dt") = g.dv_dt);
}

// Prepared conditional endpoints at latent node z under the factor model
// v_q = v + sign(rho) sv sqrt(|rho|) z, sv_q = sv max(sqrt(1-|rho|), 1e-12).
// Returns unnormalised unrestricted values plus the positivity constant so
// R tests can reassemble both drift modes.
// [[Rcpp::export]]
Rcpp::NumericVector bawl_prepared_endpoints_probe(
    double t, double t0, double A, double B, double k, double v, double sv,
    double rho, double z, bool joint_positive) {
  const DriftFactorLoading l = drift_factor_loading(sv, rho);
  const double sv_res = l.sv_res;
  const BAwLTimeGeometry g = bawl_time_geometry(t, t0, A, B + A, k);
  const BAwLPreparedRow r = bawl_prepare_row(g, v, l.slope, sv_res);
  const double vq = r.v0 + r.slope * z;
  return Rcpp::NumericVector::create(
      Rcpp::Named("vq") = vq,
      Rcpp::Named("sv_res") = sv_res,
      Rcpp::Named("log_pdf") = bawl_prepared_log_pdf(r, vq),
      Rcpp::Named("log_cdf") = bawl_prepared_log_cdf(r, vq),
      Rcpp::Named("log_q") = bawl_prepared_log_q(r, vq),
      Rcpp::Named("log_surv") =
        bawl_prepared_log_survivor(r, vq, joint_positive));
}

// Test-only accessor for the shared Gauss-Legendre cache (gl_quad.h). Lets R
// unit tests pin gl_get_rule() nodes/weights against statmod::gauss.quad.
// [[Rcpp::export]]
Rcpp::List gl_rule_nodes_weights(int n) {
  const GLRule& r = gl_get_rule(n);
  return Rcpp::List::create(
      Rcpp::Named("nodes")   = Rcpp::NumericVector(r.x.begin(), r.x.end()),
      Rcpp::Named("weights") = Rcpp::NumericVector(r.w.begin(), r.w.end()));
}

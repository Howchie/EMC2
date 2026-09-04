#include "gh_quad.h"
#include "utility_functions.h"

// ---------------------------------------------------------------------------
// DDM/BOU likelihood kernels, extracted verbatim from particle_ll.cpp.
//
// Definitions for the declarations in likelihood_ddm.h.  The shared-state
// contract (ModelSharedState and friends) lives in likelihood_shared.h; the
// adapter/context contracts stay in race_contract.h.  Private helpers
// (ddm_wien_d_raw, ddm_wien_p_raw, ddm_wiener_endpoint_key,
// ddm_wien_p_raw_cached) keep internal linkage exactly as they had in
// particle_ll.cpp.
//
// Not moved here: race_data_all_finite_untruncated, read_all_finite_trials_attr
// and build_race_shared_state remain in particle_ll.cpp with the race
// orchestration they serve.
// ---------------------------------------------------------------------------

#include "likelihood_shared.h"
#include "likelihood_ddm.h"

// The legacy vectorised wrappers in model_DDM.h are owned by particle_ll.cpp;
// this TU needs only the _raw kernels, so it opts out of re-emitting them.
#define EMC2_MODEL_DDM_NO_LEGACY_WRAPPERS
#include "model_DDM.h"            // d_DDM_Wien_raw / p_DDM_Wien_raw
#include "col_registry.h"         // emc2col::ddm column indices / ColSpec
#include "ParamTable.h"           // ParamTable (init_ddm_shared_state)
#include "fpe_models.h"           // fpe::FPE_BND_* (bou_bnd_kind_from_type)
#include "contaminant_mixture.h"  // mix_contaminants_rt / resolve_guess_kernel
#include "composite_functions.h"  // log_sum_exp / log_diff_exp / log1m_exp
#include "utils.h"                // get_col_with_default

#include <Rcpp.h>

#include <algorithm>  // std::fill, std::max
#include <array>      // std::array (single-row kernel dispatch)
#include <string>
#include <vector>

using namespace Rcpp;

bool ddm_data_all_finite_untruncated(const Rcpp::DataFrame& data,
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

// The Wiener DDM's entries in a DDMAdapter.  Thin shims that drop the model
// context the analytic model has no use for; they inline away, so the default
// adapter runs exactly the code this kernel used to call directly.
static void ddm_wien_d_raw(const double* rts, const int* Rs,
                           const double* const* cols, int n_rows,
                           const int* mask, const int* is_ok,
                           double* out, double floor_,
                           ContextForDDMModels* /*ctx*/) {
  d_DDM_Wien_raw(rts, Rs, cols, n_rows, mask, is_ok, out, floor_);
}
static void ddm_wien_p_raw(const double* rts, const int* Rs,
                           const double* const* cols, int n_rows,
                           const int* mask, const int* is_ok,
                           double* out, double floor_,
                           ContextForDDMModels* /*ctx*/) {
  p_DDM_Wien_raw(rts, Rs, cols, n_rows, mask, is_ok, out, floor_);
}

// The adapter every existing DDM caller gets when it passes none.
const DDMAdapter& ddm_wien_adapter() {
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
int bou_bnd_kind_from_type(const std::string& type_std) {
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

bool bou_anchor_at_z_from_type(const std::string& type_std) {
  return type_std.rfind("BOU", 0) == 0 &&
         type_std.size() >= 6 &&
         type_std.compare(type_std.size() - 6, 6, "_START") == 0;
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
                               Rcpp::NumericVector* trial_ll_out,
                               const DDMAdapter* ker_in) {

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
  const EndpointQueryPlan* previous_endpoint_queries =
    kctx != nullptr ? kctx->endpoint_queries : nullptr;
  EndpointQueryPlan endpoint_plan;
  if (kctx != nullptr && shared != nullptr) {
    endpoint_plan.LT = shared->LT_vec.begin();
    endpoint_plan.UT = shared->UT_vec.begin();
    endpoint_plan.LC = shared->LC_vec.begin();
    endpoint_plan.UC = shared->UC_vec.begin();
    endpoint_plan.n_rows = n_trials;
    kctx->endpoint_queries = &endpoint_plan;
  }
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
  if (kctx != nullptr && shared != nullptr)
    kctx->endpoint_queries = previous_endpoint_queries;
  return total_ll;
}

// Compatibility wrapper for old calc_ll path
double c_log_likelihood_DDM(Rcpp::NumericMatrix pars, Rcpp::DataFrame data,
                            const int n_trials, Rcpp::IntegerVector expand,
                            double min_ll, Rcpp::LogicalVector is_ok, bool gng,
                            bool all_finite_untruncated,
                            Rcpp::NumericVector* trial_ll_out) {
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

// DDM per-data shared state + canonical base-column pointers for the raw
// c_log_likelihood_DDM_pt kernel. Returns false (raw path unusable) when any
// canonical DDM parameter is missing from the table.
bool init_ddm_shared_state(DataFrame data, int n_trials,
                                  const ParamTable& table,
                                  ModelSharedState& shared,
                                  std::vector<const double*>& cols,
                                  const emc2col::ColSpec* spec_in,
                                  const Rcpp::CharacterVector* keep_names_in) {
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

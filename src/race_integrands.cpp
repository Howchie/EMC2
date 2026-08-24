#include "utility_functions.h"
#include "race_contract.h"
#include "race_integrands.h"
#include "gsl_utils.h"
#include <gsl/gsl_errno.h>

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>

#include "composite_functions.h"

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

void race_endpoint_prepare_groups(
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

double log_survivor_rowmajor(double t,
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

double log_cdf_rowmajor(double t,
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

double log_min_density_rowmajor(double t,
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

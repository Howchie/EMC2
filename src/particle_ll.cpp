#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_SDT.h"
#include "model_MRI.h"
#include "model_SS_EXG.h"
#include "model_SS_RDEX.h"
#include "composite_functions.h"
#include "trend.h"
#include "utils.h"
#include "wald_functions.h"
#include "gsl_utils.h"
#include "ParamTable.h"
#include "col_registry.h"
#include "TrendEngine.h"
#include "transform_utils.h"
#include "gh_quad.h"
#include "bawl_corr_counters.h"
#include "lr_capacity_counters.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include "bawl_geometry.h"
#include "bawl_corr_exact.h"
#include "contaminant_mixture.h"
#include <cmath>
#include <string>
#include <memory>
#include <array>
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
                           const std::unordered_set<std::string>* invariant_param_names = nullptr) {
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


// SS helper pointer types
using ss_go_pdf_fn = NumericVector (*)(NumericVector, NumericMatrix, LogicalVector, double);
using ss_stop_surv_fn = double (*)(double, NumericMatrix);
using ss_stop_success_fn = double (*)(double, NumericMatrix, double, double, int, double, double, double, double);

// Model-specific stop survivor wrappers (columns per src/col_registry.h)
static inline double stop_logsurv_texg_fn(double q, NumericMatrix P) {
  namespace tc = emc2col::ss_texg;
  return ptexg(q, P(0, tc::muS), P(0, tc::sigmaS), P(0, tc::tauS),
               P(0, tc::exgS_lb), R_PosInf, false, true);
}
static inline double stop_logsurv_rdex_fn(double q, NumericMatrix P) {
  namespace rc = emc2col::ss_rdex;
  return ptexg(q, P(0, rc::muS), P(0, rc::sigmaS), P(0, rc::tauS),
               P(0, rc::exgS_lb), R_PosInf, false, true);
}

static inline double sum_log_terms(const NumericVector& x) {
  double out = 0.0;
  for (int i = 0; i < x.size(); ++i) out += x[i];
  return out;
}

static inline double ss_log_go_density_for_winner(
    double rt,
    const NumericMatrix& P,
    const LogicalVector& is_go,
    int winner_idx,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    double min_ll
) {
  const int n_acc = P.nrow();
  NumericVector rt_go(n_acc, rt);
  LogicalVector go_win_mask(n_acc, false);
  LogicalVector go_loss_mask(n_acc, false);
  int go_loss_count = 0;
  for (int i = 0; i < n_acc; ++i) {
    const bool go_i = (is_go[i] == TRUE);
    if (!go_i) continue;
    if (i == winner_idx) go_win_mask[i] = true;
    else {
      go_loss_mask[i] = true;
      ++go_loss_count;
    }
  }

  NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
  double out = (lw.size() > 0) ? sum_log_terms(lw) : R_NegInf;
  if (!R_FINITE(out)) return R_NegInf;

  if (go_loss_count == 0) return out;

  NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
  return out + sum_log_terms(ls);
}

static inline double ss_log_st_density_for_winner(
    double rt,
    double SSD,
    const NumericMatrix& P,
    const LogicalVector& is_st,
    int winner_idx,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    double min_ll
) {
  const int n_acc = P.nrow();
  double rt_st_val = rt - SSD;
  if (rt_st_val < 0.0) rt_st_val = 0.0;
  NumericVector rt_st(n_acc, rt_st_val);
  LogicalVector st_win_mask(n_acc, false);
  LogicalVector st_loss_mask(n_acc, false);
  int st_loss_count = 0;
  for (int i = 0; i < n_acc; ++i) {
    const bool st_i = (is_st[i] == TRUE);
    if (!st_i) continue;
    if (i == winner_idx) st_win_mask[i] = true;
    else {
      st_loss_mask[i] = true;
      ++st_loss_count;
    }
  }

  NumericVector lw = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
  double out = (lw.size() > 0) ? sum_log_terms(lw) : R_NegInf;
  if (!R_FINITE(out)) return R_NegInf;

  if (st_loss_count == 0) return out;

  NumericVector ls = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
  return out + sum_log_terms(ls);
}

static inline double ss_trial_log_response_density(
    double rt,
    const NumericMatrix& P,
    const NumericVector& lR_trial,
    const LogicalVector& is_go,
    const LogicalVector& is_st,
    double SSD,
    bool stop_signal_presented,
    double tf,
    double gf,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    double min_ll,
    int response_code = NA_INTEGER
) {
  const int n_acc = P.nrow();
  const bool response_known = response_code != NA_INTEGER;
  int n_accG = 0;
  int n_accST = 0;
  for (int i = 0; i < n_acc; ++i) {
    if (is_go[i] == TRUE) ++n_accG;
    if (is_st[i] == TRUE) ++n_accST;
  }

  double log_go_any = R_NegInf;
  for (int i = 0; i < n_acc; ++i) {
    if (is_go[i] != TRUE) continue;
    if (response_known && lR_trial[i] != response_code) continue;

    const double go_lprob = ss_log_go_density_for_winner(
      rt, P, is_go, i, go_lpdf_ptr, go_lccdf_ptr, min_ll
    );
    if (!R_FINITE(go_lprob)) continue;

    double log_term = log1m(gf) + go_lprob;
    if (stop_signal_presented && rt > SSD) {
      double rt_eff = rt - SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = stop_logsurv_ptr(rt_eff, P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;

      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        NumericVector rt_st(n_acc, rt_eff);
        LogicalVector st_loss_mask(n_acc, false);
        int st_loss_count = 0;
        for (int j = 0; j < n_acc; ++j) {
          st_loss_mask[j] = (is_st[j] == TRUE);
          if (st_loss_mask[j] == TRUE) ++st_loss_count;
        }
        if (st_loss_count > 0) {
          NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
          st_loss_sum = sum_log_terms(ls_st);
        }
      }

      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      log_term = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
    }

    log_go_any = log_sum_exp(log_go_any, log_term);
  }

  if (!stop_signal_presented || n_accST == 0) return log_go_any;

  double log_st_any = R_NegInf;
  NumericMatrix P_go = submat_rcpp(P, is_go);
  for (int i = 0; i < n_acc; ++i) {
    if (is_st[i] != TRUE) continue;
    if (response_known && lR_trial[i] != response_code) continue;
    if (rt <= SSD) continue;

    const double st_base = ss_log_st_density_for_winner(
      rt, SSD, P, is_st, i, go_lpdf_ptr, go_lccdf_ptr, min_ll
    );
    if (!R_FINITE(st_base)) continue;

    double go_loss_sum = 0.0;
    if (n_accG > 0) {
      NumericVector rt_go(n_acc, rt);
      NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
      go_loss_sum = sum_log_terms(ls_go);
    }

    double rt_eff = rt - SSD;
    if (rt_eff < 0.0) rt_eff = 0.0;
    double log_pstop = stop_success_ptr(SSD, P_go, min_ll, rt_eff,
                                        100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
    if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

    const double term_gf = std::log(gf) + st_base;
    const double term_stop_win = log1m(gf) + log_pstop + st_base;
    const double term_stop_lose =
      log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
    const double log_term =
      log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));

    log_st_any = log_sum_exp(log_st_any, log_term);
  }

  return log_sum_exp(log_go_any, log_st_any);
}

struct SsLcIntegrandData {
  NumericMatrix P;
  NumericVector lR_trial;
  LogicalVector is_go;
  LogicalVector is_st;
  double SSD;
  bool stop_signal_presented;
  double tf;
  double gf;
  ss_go_pdf_fn go_lpdf_ptr;
  ss_go_pdf_fn go_lccdf_ptr;
  ss_stop_surv_fn stop_logsurv_ptr;
  ss_stop_success_fn stop_success_ptr;
  double min_ll;
  int response_code;
};

static double ss_lc_integrand(double x, void* params) {
  SsLcIntegrandData* d = static_cast<SsLcIntegrandData*>(params);
  const double logf = ss_trial_log_response_density(
    x, d->P, d->lR_trial, d->is_go, d->is_st, d->SSD, d->stop_signal_presented,
    d->tf, d->gf, d->go_lpdf_ptr, d->go_lccdf_ptr, d->stop_logsurv_ptr,
    d->stop_success_ptr, d->min_ll, d->response_code
  );
  if (!R_FINITE(logf)) return 0.0;
  return std::exp(logf);
}

static inline double ss_integrate_lc_response_mass(
    double lower,
    double upper,
    const NumericMatrix& P,
    const NumericVector& lR_trial,
    const LogicalVector& is_go,
    const LogicalVector& is_st,
    double SSD,
    bool stop_signal_presented,
    double tf,
    double gf,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    double min_ll,
    int response_code = NA_INTEGER
) {
  if (!(upper > lower)) return min_ll;

  SsLcIntegrandData data{
    P, lR_trial, is_go, is_st, SSD, stop_signal_presented, tf, gf,
    go_lpdf_ptr, go_lccdf_ptr, stop_logsurv_ptr, stop_success_ptr, min_ll,
    response_code
  };

  gsl_function F;
  F.function = &ss_lc_integrand;
  F.params = &data;

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, 200);
  double res = 0.0;
  double err = 0.0;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags(&F, lower, upper, 1e-8, 1e-6, 200, workspace, &res, &err);
  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !R_FINITE(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

struct RaceModelAdapter {
  RacePdf1Fun pdf1_ptr = nullptr;
  RaceCdf1Fun cdf1_ptr = nullptr;
  RaceRawFun model_dfun_raw = nullptr;   // fast-path: write log-density to pre-allocated buffer
  RaceRawFun model_pfun_raw = nullptr;   // fast-path: write log-survivor to pre-allocated buffer
  RaceLogSAtTFun logS_at_t_ptr = nullptr; // batch: log-survivor at scalar t for truncation norms
  emc2col::ColSpec col_spec = {nullptr, 0, ""}; // kernel column contract (col_registry.h)
  ContextForRaceModels ctx;
};

static inline RaceModelAdapter resolve_race_model_adapter(const std::string& type_std,
                                                          const std::string& caller) {
  RaceModelAdapter out;
  out.ctx.min_lik_for_pdf = std::exp(std::log(1e-10));
  out.ctx.use_posdrift = true;
  out.ctx.gng = false;
  out.ctx.kill_shape = (type_std.find("_EMIX") != std::string::npos) ? 3 :
                       ((type_std.find("_E2") != std::string::npos) ? 2 : 1);

  // Erlang process flag resolution
  out.ctx.is_local_guess = (type_std.find("_LOCAL_GUESS") != std::string::npos);
  out.ctx.is_global_kill = (type_std.find("_GLOBAL_KILL") != std::string::npos);
  out.ctx.is_local_kill  = (type_std.find("_LOCAL_KILL")  != std::string::npos);
  out.ctx.is_local_kill_guess = (type_std.find("_LOCAL_KILL_GUESS") != std::string::npos);

  // Backward compatibility: bare _GLOBAL maps to global_kill
  if (!out.ctx.is_local_guess && !out.ctx.is_global_kill && !out.ctx.is_local_kill && !out.ctx.is_local_kill_guess) {
    if (type_std.find("_GLOBAL") != std::string::npos) out.ctx.is_global_kill = true;
  }

  // A bare BAwL/RDM-with-timers type has timer columns in its parameter
  // contract, but no active Erlang process.  Keep that distinction in the
  // scalar and raw paths; treating the default mG/mK values as active clocks
  // changes the likelihood and needlessly evaluates the clock mixture for
  // every proposal.  The suffixed variants enable the process below.
  out.ctx.kill_active = out.ctx.is_local_guess || out.ctx.is_global_kill ||
                        out.ctx.is_local_kill || out.ctx.is_local_kill_guess;

  out.ctx.apply_lk_to_racers = !out.ctx.is_global_kill;

  if (out.ctx.is_global_kill) out.ctx.defective_upper_tail = true;

  if (type_std.find("RLF") != std::string::npos) {
    out.pdf1_ptr       = &drlf_scalar;
    out.cdf1_ptr       = &prlf_scalar;
    out.model_dfun_raw = &drlf_raw;
    out.model_pfun_raw = &prlf_raw;
    out.logS_at_t_ptr  = &rlf_logS_at_t;
    out.col_spec       = emc2col::rlf::spec();
    out.ctx.t0_index   = emc2col::rlf::t0;
    out.ctx.rlf_cache = std::make_shared<rlf::SolveCache>();
    rlf_configure_grid(out.ctx.rlf_cache->grid);
  } else if (type_std.find("GOM") != std::string::npos ||
             type_std.find("GOMP") != std::string::npos) {
    // Gompertz growth is an OU after Y = log(X).  The adapter keeps the
    // physical alpha/beta/K columns and gomp_key() performs the reduction; its
    // log_state flag makes the FPE seed integrate A uniformly on X rather than
    // incorrectly treating the transformed start range as uniform in Y.
    out.pdf1_ptr       = &dgomp_scalar;
    out.cdf1_ptr       = &pgomp_scalar;
    out.model_dfun_raw = &dgomp_raw;
    out.model_pfun_raw = &pgomp_raw;
    out.logS_at_t_ptr  = &gomp_logS_at_t;
    out.col_spec       = emc2col::gompertz::spec();
    out.ctx.t0_index   = emc2col::gompertz::t0;
    out.ctx.defective_upper_tail = false;
    out.ctx.fpe_cache = std::make_shared<fperace::SolveCache>();
    rou_configure_cache(*out.ctx.fpe_cache);
    if (type_std.find("_BWEIB") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_WEIBULL;
    else if (type_std.find("_BEXP") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_EXPONENTIAL;
    else if (type_std.find("_BLIN_MULT") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
    else if (type_std.find("_BLIN_ADD") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_ADDITIVE;
  } else if (type_std.find("ROU") != std::string::npos) {
    // Ordered FIRST deliberately.  Dispatch here is by substring, so a key that
    // is a substring of a later one must be tested first; "ROU" collides with
    // nothing today, and testing it first is what keeps a future addition from
    // silently capturing it.
    out.pdf1_ptr       = &drou_scalar;
    out.cdf1_ptr       = &prou_scalar;
    out.model_dfun_raw = &drou_raw;
    out.model_pfun_raw = &prou_raw;
    out.logS_at_t_ptr  = &rou_logS_at_t;
    // With non-zero diffusion and a finite upper boundary, the OU hits the
    // boundary almost surely.  A subthreshold equilibrium creates long finite
    // survival and an approximately constant late hazard, not a point mass at
    // infinity.  Preserve a genuine global-kill defect if a suffixed variant
    // requests that clock.
    out.ctx.defective_upper_tail = out.ctx.is_global_kill;
    out.ctx.fpe_cache = std::make_shared<fperace::SolveCache>();
    rou_configure_cache(*out.ctx.fpe_cache);
    // Parameterisation, selected by ROU(parameterization=).  It changes only
    // which columns the kernels read and how they map onto (v, k, s); the
    // solve, the cache and the boundary handling are untouched.  Tested before
    // the boundary suffix because the two compose: ROUCURV_BEXP is both.
    if (type_std.find("ROUCURV") != std::string::npos) {
      out.ctx.fpe_cache->par_kind = fperace::ROU_PAR_CURVATURE;
      out.col_spec     = emc2col::rou_curv::spec();
      out.ctx.t0_index = emc2col::rou_curv::t0;
    } else if (type_std.find("ROUEQ") != std::string::npos) {
      out.ctx.fpe_cache->par_kind = fperace::ROU_PAR_EQUILIBRIUM;
      out.col_spec     = emc2col::rou_eq::spec();
      out.ctx.t0_index = emc2col::rou_eq::t0;
    } else {
      out.ctx.fpe_cache->par_kind = fperace::ROU_PAR_RATE;
      out.col_spec     = emc2col::rou::spec();
      out.ctx.t0_index = emc2col::rou::t0;
    }
    // Collapsing-bound variants, selected by ROU(boundary_collapse=).  The
    // suffix carries the FORM only; the shape parameters are ordinary optional
    // columns (Binf/tau/pw) that the design system estimates like any other.
    // Matched on the bare suffix rather than on "ROU_B*" so that the collapse
    // forms are available under every parameterisation; we are already inside
    // the ROU branch, so there is nothing else these can capture.
    if (type_std.find("_BWEIB") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_WEIBULL;
    else if (type_std.find("_BEXP") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_EXPONENTIAL;
    else if (type_std.find("_BLIN_MULT") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
    else if (type_std.find("_BLIN_ADD") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_ADDITIVE;
  } else if (type_std.find("RDMSWTN_TT") != std::string::npos) {
    // The time-changed model must precede generic RDMSWTN substring dispatch.
    out.pdf1_ptr       = &drdmswtn_tt_scalar;
    out.cdf1_ptr       = &prdmswtn_tt_scalar;
    out.model_dfun_raw = &drdmswtn_tt_raw;
    out.model_pfun_raw = &prdmswtn_tt_raw;
    out.logS_at_t_ptr  = &rdmswtn_tt_logS_at_t;
    out.col_spec       = emc2col::rdmswtn_tt::spec();
    out.ctx.t0_index   = emc2col::rdmswtn_tt::t0;
    out.ctx.defective_upper_tail = true;
    out.ctx.rdmswtn_correlated =
      (type_std.find("_CORR") != std::string::npos);
    if (type_std.find("_IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("RDMSWTN") != std::string::npos) {
    // Must be checked before "RDM" since "RDMSWTN" contains "RDM"
    out.pdf1_ptr       = &drdmswtn_scalar;
    out.cdf1_ptr       = &prdmswtn_scalar;
    out.model_dfun_raw = &drdmswtn_raw;
    out.model_pfun_raw = &prdmswtn_raw;
    out.logS_at_t_ptr  = &rdmswtn_logS_at_t;
    out.col_spec       = emc2col::rdmswtn::spec();
    out.ctx.t0_index   = emc2col::rdmswtn::t0;
    out.ctx.mean_g_index = emc2col::rdmswtn::mG;
    out.ctx.mean_k_index = emc2col::rdmswtn::mK;
    out.ctx.erlang_omega_index = (out.ctx.kill_shape == 3) ? emc2col::rdmswtn::omega : -1;
    out.ctx.defective_upper_tail = true;
    out.ctx.rdmswtn_correlated =
      (type_std.find("_CORR") != std::string::npos);
    if (type_std.find("_IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("GBM") != std::string::npos) {
    // Must be checked before "RDM" since "RDMGBM" contains "RDM"
    out.pdf1_ptr       = &drdmgbm_scalar;
    out.cdf1_ptr       = &prdmgbm_scalar;
    out.model_dfun_raw = &drdmgbm_raw;
    out.model_pfun_raw = &prdmgbm_raw;
    out.logS_at_t_ptr  = &rdmgbm_logS_at_t;
    out.col_spec       = emc2col::rdmgbm::spec();
    out.ctx.t0_index     = emc2col::rdmgbm::t0;
    out.ctx.mean_g_index = emc2col::rdmgbm::mG;
    out.ctx.mean_k_index = emc2col::rdmgbm::mK;
    out.ctx.erlang_omega_index = (out.ctx.kill_shape == 3) ? emc2col::rdmgbm::omega : -1;
    out.ctx.defective_upper_tail = true;
  } else if (type_std.find("BAwD") != std::string::npos) {
    // Dispatch is by substring, and "BAwD" is a substring of nothing here and
    // contains neither "BAwL" nor "LBA", so placement relative to those is
    // safe either way; it sits next to BAwL for readability.
    out.pdf1_ptr       = &dbawd_scalar;
    out.cdf1_ptr       = &pbawd_scalar;
    out.model_dfun_raw = &dbawd_raw;
    out.model_pfun_raw = &pbawd_raw;
    out.logS_at_t_ptr  = &bawd_logS_at_t;
    // The two launch distributions share column POSITIONS and differ only in
    // the names validate_col_prefix() enforces, so t0_index is common.
    const bool bawd_logn = (type_std.find("_LOGN") != std::string::npos);
    out.col_spec = bawd_logn ? emc2col::bawd_logn::spec() : emc2col::bawd::spec();
    out.ctx.t0_index = emc2col::bawd::t0;
    out.ctx.bawd_launch = bawd_logn ? BAWD_LAUNCH_LOGNORMAL : BAWD_LAUNCH_NORMAL;
    // Always defective: a decaying drive with a constant clearance has a hard
    // right endpoint and intrinsic never-finish mass regardless of posdrift.
    out.ctx.defective_upper_tail = true;
    // posdrift is meaningless for the lognormal launch (V > 0 by construction);
    // BAwD() refuses posdrift = FALSE there rather than silently ignoring it.
    if (!bawd_logn && type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("BAwL") != std::string::npos) {
    out.pdf1_ptr       = &dbawl_scalar;
    out.cdf1_ptr       = &pbawl_scalar;
    out.model_dfun_raw = &dbawl_raw;
    out.model_pfun_raw = &pbawl_raw;
    out.logS_at_t_ptr  = &bawl_logS_at_t;
    out.col_spec       = emc2col::bawl::spec();
    out.ctx.t0_index   = emc2col::bawl::t0;
    out.ctx.mean_g_index = emc2col::bawl::mG;
    out.ctx.mean_k_index = emc2col::bawl::mK;
    out.ctx.erlang_omega_index = (out.ctx.kill_shape == 3) ? emc2col::bawl::omega : -1;
    out.ctx.bawl_correlated = (type_std.find("_CORR") != std::string::npos);
    // Leaky ballistic accumulators can have defective upper tails (never-finish
    // mass) even when posdrift=TRUE.
    out.ctx.defective_upper_tail = true;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
    if (type_std.find("_E2") != std::string::npos) out.ctx.kill_shape = 2;
    if (type_std.find("_EMIX") != std::string::npos) out.ctx.kill_shape = 3;
  } else if (type_std.find("LBA") != std::string::npos) {
    // Standard LBA is the exact k=0, no-clock member of the shared BAwL
    // family.  Keep the five-column LBA contract and force the optional BAwL
    // parameters off in the adapter rather than adding hidden columns.
    out.pdf1_ptr = &dbawl_scalar;
    out.cdf1_ptr = &pbawl_scalar;
    out.model_dfun_raw = &dbawl_raw;
    out.model_pfun_raw = &pbawl_raw;
    out.logS_at_t_ptr = &bawl_logS_at_t;
    out.col_spec     = emc2col::lba::spec();
    out.ctx.t0_index = emc2col::lba::t0;
    out.ctx.bawl_k_fixed_zero = true;
    out.ctx.bawl_clocks_fixed_off = true;
    out.ctx.defective_upper_tail = false;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
      out.ctx.defective_upper_tail = true;
    }
  } else if (type_std.find("RDM") != std::string::npos) {
    out.pdf1_ptr = &drdm_scalar;
    out.cdf1_ptr = &prdm_scalar;
    out.model_dfun_raw = &drdm_raw;
    out.model_pfun_raw = &prdm_raw;
    out.logS_at_t_ptr = &rdm_logS_at_t;
    out.col_spec     = emc2col::rdm::spec();
    out.ctx.t0_index = emc2col::rdm::t0;
  } else if (type_std.find("REXG") != std::string::npos) {
    out.pdf1_ptr = &dexg_scalar;
    out.cdf1_ptr = &pexg_scalar;
    out.model_dfun_raw = &drexg_raw;
    out.model_pfun_raw = &prexg_raw;
    out.logS_at_t_ptr = &rexg_logS_at_t;
    out.col_spec     = emc2col::rexg::spec();
    out.ctx.t0_index = -1;
  } else if (type_std.find("LNR") != std::string::npos) {
    out.pdf1_ptr = &dlnr_scalar;
    out.cdf1_ptr = &plnr_scalar;
    out.model_dfun_raw = &dlnr_raw;
    out.model_pfun_raw = &plnr_raw;
    out.logS_at_t_ptr = &lnr_logS_at_t;
    out.col_spec     = emc2col::lnr::spec();
    out.ctx.t0_index = emc2col::lnr::t0;
  } else if (type_std.find("PCOUNTER") != std::string::npos) {
    out.pdf1_ptr = &dpcounter_scalar;
    out.cdf1_ptr = &ppcounter_scalar;
    out.model_dfun_raw = &dpcounter_raw;
    out.model_pfun_raw = &ppcounter_raw;
    out.logS_at_t_ptr = &pcounter_logS_at_t;
    out.col_spec     = emc2col::pcounter::spec();
    out.ctx.t0_index = -1;
  } else {
    Rcpp::stop("Unsupported race model type string in %s: %s", caller.c_str(), type_std.c_str());
  }

  if (type_std.find("GNG") != std::string::npos) {
    out.ctx.gng = true;
  }
  return out;
}

static inline void configure_bawl_corr_context(RaceModelAdapter& adapter,
                                               const Rcpp::CharacterVector& keep_names,
                                               const std::string& caller) {
  if (!adapter.ctx.bawl_correlated) return;
  adapter.ctx.bawl_rho_index = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == "rho") {
      adapter.ctx.bawl_rho_index = j;
      break;
    }
  }
  if (adapter.ctx.bawl_rho_index < 0) {
    Rcpp::stop("%s: correlated BAwL requires a parameter column named 'rho'.",
               caller.c_str());
  }
}

static inline void configure_rdmswtn_corr_context(
    RaceModelAdapter& adapter, const Rcpp::CharacterVector& keep_names,
    const std::string& caller) {
  if (!adapter.ctx.rdmswtn_correlated) return;
  adapter.ctx.rdmswtn_rho_index = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == "rho") {
      adapter.ctx.rdmswtn_rho_index = j;
      break;
    }
  }
  if (adapter.ctx.rdmswtn_rho_index < 0) {
    Rcpp::stop("%s: correlated RDMSWTN requires a parameter column named 'rho'.",
               caller.c_str());
  }
  if (adapter.ctx.kill_active) {
    Rcpp::stop("%s: correlated RDMSWTN does not support guess or kill clocks.",
               caller.c_str());
  }
}

static inline bool is_stop_signal_type(const std::string& type_std) {
  return type_std == "SSEXG" || type_std == "SSRDEX";
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

static LogicalRulesSharedState build_logicalrules_shared_state(const Rcpp::DataFrame& dadm,
                                                               int n_trials,
                                                               int n_acc,
                                                               bool capacity = false);

static double local_race_helper(double t,
                                const double* par_target,
                                const double* par_nontarget,
                                int n_par,
                                ContextForRaceModels* ctx,
                                RacePdf1Fun pdf1,
                                RaceCdf1Fun cdf1,
                                const GslIntegrationControls& gsl_ctl,
                                gsl_integration_workspace* w,
                                double* pars_2buf,
                                int* isok_2buf);

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

// --------------------------------------------------------------------------
// Correlated BAwL shared state and canonical trial layout
// (bawl_corr_exact_kernel_plan.md, design rule C).  Built once per likelihood
// call outside the particle loop; the per-particle classifier below is the
// sole authority for loaded-row discovery and positivity dimension.
// --------------------------------------------------------------------------

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

struct BAwLCorrSharedState {
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

static BAwLCorrSharedState build_bawl_corr_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names);

double c_log_likelihood_bawl_correlated(
    BAwLCorrSharedState& cshared,
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
    RaceSharedState* shared,
    NumericVector* trial_ll_out,
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

static RDMSWTNCorrSharedState build_rdmswtn_corr_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names);

double c_log_likelihood_rdmswtn_correlated(
    RDMSWTNCorrSharedState& cshared, Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, const int n_trials,
    LogicalVector winner, Rcpp::IntegerVector expand, double min_ll,
    const Rcpp::LogicalVector isok, int n_lR,
    void* model_context_for_funcs, bool all_finite_trials,
    RaceRawFun model_dfun_raw, RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t, RaceSharedState* shared,
    NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize);

static double c_log_likelihood_logicalrules(
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

static inline bool eval_pdf_cdf_race_scalar(
    bool use_raw_local,
    int idx,
    double t,
    const double* p_row,
    double min_ll,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const std::vector<double>& logf_all,
    const std::vector<double>& logS_all,
    double& f_out,
    double& F_out);

static inline double safe_surv_at_race_scalar(
    const double* p_row,
    double t_eval,
    double min_surv,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx);

static double logicalrules_detection_trial_ll(
    int rule_code,
    int cond,
    int resp,
    double t,
    double LTj,
    double UTj,
    double LCj,
    double UCj,
    int idxA,
    int idxB,
    int idxN,
    int n_par,
    double min_ll,
    double min_surv,
    double log_Z,
    bool use_raw_local,
    const std::vector<double>& logf_all,
    const std::vector<double>& logS_all,
    const double* parA,
    const double* parB,
    const double* parN,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& pars_rowmajor_buf);

double c_log_likelihood_ss(
    NumericMatrix pars,
    DataFrame data,
    const int n_trials,
    IntegerVector expand,
    double min_ll,
    LogicalVector is_ok,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    int idx_tf,
    int idx_gf
) {

  
  // initialise local variables
  const int n_out = expand.length();
  if (is_true(all(!is_ok))) {
    return static_cast<double>(n_out) * min_ll;
  }
  NumericVector lls(n_trials);
  // extract data
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector LT = get_col_with_default(data, "LT", 0.0);
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? as<IntegerVector>(data["lI"]) : IntegerVector(lR.size(), 2);
  
  // dimensional expectations: pars has one row per accumulator per trial
  
  // compute log likelihoods (generalized, matching R's log_likelihood_race_ss)
  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();
  
  NumericVector tt(n_acc);
  auto log_surv_mask = [&](double t, const NumericMatrix& Pcur,
                           const LogicalVector& mask) -> double {
                             tt.fill(t);
                             NumericVector ls = go_lccdf_ptr(tt, Pcur, mask, min_ll);
                             double out = 0.0;
                             for (int i = 0; i < ls.size(); ++i) {
                               out += ls[i];
                             }
                             return out;
                           };
  // n_trials equals data rows grouped by accumulators
  for (int trial = 0; trial < n_trials; trial++) {
    if (is_ok[trial * n_acc] != 1) { lls[trial] = min_ll; continue; }
    
    int start_row = trial * n_acc;
    int end_row   = (trial + 1) * n_acc - 1;
    // basic bounds are guaranteed by correct n_trials passed into this function
    NumericMatrix P = pars(Range(start_row, end_row), _);
    NumericVector lR_trial = lR[Range(start_row, end_row)];
    IntegerVector lI_trial = lI[Range(start_row, end_row)];
    LogicalVector is_go(n_acc, true), is_st(n_acc, false);
    // determine go/ST accumulators if present
    if (has_lI) {
      int go_code = max(lI_trial);
      for (int i = 0; i < n_acc; i++) {
        is_go[i] = (lI_trial[i] == go_code);
        is_st[i] = !is_go[i];
      }
    } else {
      for (int i = 0; i < n_acc; i++) is_st[i] = false;
    }
    int n_accG = sum(is_go);
    int n_accST = sum(is_st);
    
    double tf = P(0, idx_tf);
    double gf = P(0, idx_gf);
    
    double rt = RT[start_row];
    bool response_observed = R[start_row] != NA_INTEGER;
    bool stop_signal_presented = emc2_isfinite(SSD[start_row]);
    // Added UC handling
    double uc = UC[start_row];
    // Identify whether observed response is GO or ST (when response observed)
    bool response_is_go = false;
    if (response_observed) {
      int r_obs = R[start_row];
      for (int i = 0; i < n_acc; i++) {
        if (lR[start_row + i] == r_obs) {
          response_is_go = is_go[i];
          break;
        }
      }
    }
    
    // Build rt vectors for go and st contexts
    NumericVector rt_go(n_acc, rt);
    double rt_st_val = rt - SSD[start_row];
    if (rt_st_val < 0.0) rt_st_val = 0.0;
    NumericVector rt_st(n_acc, rt_st_val);
    
    // GO masks for current trial
    LogicalVector win_mask = winner[Range(start_row, end_row)];
    LogicalVector go_win_mask(n_acc); // winner and go
    LogicalVector go_loss_mask(n_acc);
    for (int i = 0; i < n_acc; i++) {
      go_win_mask[i] = (win_mask[i] && is_go[i]);
      go_loss_mask[i] = (!win_mask[i] && is_go[i]);
    }

    bool lower_censored = (rt == R_NegInf);
    if (lower_censored) {
      const double upper = LC[start_row];
      lls[trial] = ss_integrate_lc_response_mass(
        LT[start_row], upper, P, lR_trial, is_go, is_st, SSD[start_row], stop_signal_presented,
        tf, gf, go_lpdf_ptr, go_lccdf_ptr, stop_logsurv_ptr, stop_success_ptr,
        min_ll, response_observed ? R[start_row] : NA_INTEGER
      );
      continue;
    }
    
    if (!response_observed) {
      bool has_deadline = R_FINITE(uc) && !Rcpp::NumericVector::is_na(uc);
      if (!stop_signal_presented) {
        // No response
        if (!has_deadline) {
          // intrinsic NR
          lls[trial] = std::log(gf);
        } else {
          // deadline-censored GO: gf + (1-gf)*S_go(UC)
          // GO trial, no response by UC: gf + (1-gf)*S_go(UC)
          double logS_go = (n_accG > 0) ? log_surv_mask(uc, P, is_go) : 0.0; // product of all go accumulators not having finished
          lls[trial] = log_sum_exp(std::log(gf), log1m(gf) + logS_go);
        }
        continue;
      }
      // stop trial, no response observed
      if (!has_deadline) {
        if (n_accST == 0) {
          // Stop trial, no ST accumulators: gf + (1-gf)*(1-tf)*pStop
          NumericMatrix P_go = submat_rcpp(P, is_go);
          double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, R_PosInf,
                                              100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
          if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
          double comp1 = std::log(gf);
          double comp2 = log1m(gf) + log1m(tf) + log_pstop;
          lls[trial] = log_sum_exp(comp1, comp2);
        } else {
          // With an ST accumulator present, a stop-trial non-response (with no deadline)
          // can only happen via trigger-failure AND go-failure: tf * gf.
          lls[trial] = std::log(gf) + std::log(tf);
        }
        continue;
      }
      
      // Deadline-censored stop trial: "stop win by UC OR Go/St unfinished by UC"
      // GO survivor by UC
      double logS_go = (n_accG > 0) ? log_surv_mask(uc, P, is_go) : 0.0;
      
      // stop survivor by UC (duration scale: UC-SSD, clamp at 0)
      double uc_eff = uc-SSD[start_row];
      if (!R_FINITE(uc_eff) || uc_eff <= 0.0) uc_eff = 0.0;
      bool stop_can_act = R_FINITE(uc_eff) && (uc_eff > 0.0);
      
      double log_pstop = R_NegInf;
      double logS_stop = 0.0; // log(1)
      
      // pStop(UC): stop finishes before GO and before UC-SSD
      // integral from SSD to UC (integrand function adds SSD to go accumulators so we use uc here)
      NumericMatrix P_go = submat_rcpp(P, is_go);
      if (stop_can_act) {
        logS_stop = stop_logsurv_ptr(uc_eff, P);
        log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, uc_eff, 100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      }
      
      // Triggered, go-not-failed core no-response by UC:
      // p = pStop(UC) + S_go(UC)*S_stop(UC-SSD)
      double log_core_trig = log_sum_exp(log_pstop, logS_go + logS_stop);
      
      if (n_accST == 0) {
        // No ST accumulators:
        // p = gf + (1-gf) * [ tf*S_go + (1-tf)*core_trig ]
        double log_no_nogf = log_mix(tf, logS_go, log_core_trig);
        lls[trial] = log_sum_exp(std::log(gf), log1m(gf) + log_no_nogf);
        continue;
      }
      
      // Stop-triggered case - 
      double logS_st = log_surv_mask(uc_eff, P, is_st);
      double log_trig = logS_st + log_sum_exp(std::log(gf), log1m(gf) + log_core_trig);
      
      double log_tfbranch = log_sum_exp(std::log(gf), log1m(gf) + logS_go);
      
      lls[trial] = log_mix(tf, log_tfbranch, log_trig);
      continue;
    }
    
    // Response observed
    if (!stop_signal_presented) {
      // GO trial with response: (1-gf) * GO race ll
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : R_NegInf;
      if (!R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) {
          go_lprob += ls[i];
        }
      }
      lls[trial] = log1m(gf) + go_lprob;
      continue;
    }
    
    // Stop trial with response
    if (response_is_go) {
      // GO wins on stop trial: (1-gf) * [ tf * go + (1-tf) * (go + stop_surv + st_loss) ]
      // go race ll at observed rt
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : R_NegInf;
      if (!R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) {
          go_lprob += ls[i];
        }
      }
      // stop survivor at observed rt (rt - SSD)
      double rt_eff = rt - SSD[start_row];
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = stop_logsurv_ptr(rt_eff, P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;
      // ST losers survivors (if any)
      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (is_st[i] && !win_mask[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) {
          st_loss_sum += ls_st[i];
        }
      }
      double comp_tf = go_lprob; // only go
      double comp_notf = go_lprob + log_stop_surv + st_loss_sum; // fair race (stop loses)
      lls[trial] = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
      continue;
    } else {
      // ST wins on stop trial
      // ST winner log pdf at rt - SSD
      LogicalVector st_win_mask(n_acc);
      for (int i = 0; i < n_acc; i++) st_win_mask[i] = (win_mask[i] && is_st[i]);
      NumericVector lw_st = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
      double st_winner_logpdf = (lw_st.size() > 0) ? sum(lw_st) : R_NegInf;
      if (!R_FINITE(st_winner_logpdf)) st_winner_logpdf = R_NegInf;
      // ST losers survivors
      double st_loss_sum = 0.0;
      if (n_accST > 1) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (!win_mask[i] && is_st[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) {
          st_loss_sum += ls_st[i];
        }
      }
      // GO losers survivors
      double go_loss_sum = 0.0;
      if (n_accG > 0) {
        NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
        for (int i = 0; i < ls_go.size(); i++) {
          go_loss_sum += ls_go[i];
        }
      }
      // Stop success probability up to observed rt (only go racers influence integral)
      NumericMatrix P_go = submat_rcpp(P, is_go);
      // stop survivor at observed rt (rt - SSD)
      double rt_eff = rt - SSD[start_row];
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, rt_eff,
                                          100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
      if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      
      double st_base = st_winner_logpdf + st_loss_sum;
      double term_gf = std::log(gf) + st_base; // go failure -> only ST race
      double term_stop_win = log1m(gf) + log_pstop + st_base; // stop beats go -> only ST race
      double term_stop_lose = log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum; // all race, no stop win
      lls[trial] = log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
      continue;
    }
  }
  lls[is_na(lls)] = min_ll;
  lls[is_infinite(lls)] = min_ll;
  lls[lls < min_ll] = min_ll;
  double sum_ll = 0.0;
  for (int i = 0; i < n_out; ++i) {
    sum_ll += lls[expand[i] - 1]; // expand created in 1-based R
  }
  return sum_ll;
}

// ============================================================================
// Raw-buffer stop-signal path.
//
// c_log_likelihood_ss above re-derived the whole trial structure from Rcpp
// columns for every trial of every particle (Range() submatrix copies,
// LogicalVector masks, NumericVector rt broadcasts). Everything data-fixed is
// now computed ONCE per calc_ll_oo call into SsSharedState, and the
// per-particle evaluator reads ParamTable base columns directly through the
// canonical kernel blocks defined in ss_raw.h (see SsRawModel). Numerically
// identical to c_log_likelihood_ss, which remains as the fallback (and
// readable reference) when the ParamTable column mapping is unavailable.
// ============================================================================

enum SsTrialKind : uint8_t {
  SS_TRIAL_LC_CENSORED = 0,  // rt == -Inf: integrate response mass over [LT, LC)
  SS_TRIAL_GO_NR_FREE,       // go trial, no response, no deadline
  SS_TRIAL_GO_NR_DEADLINE,   // go trial, no response by UC
  SS_TRIAL_STOP_NR_FREE,     // stop trial, no response, no deadline
  SS_TRIAL_STOP_NR_DEADLINE, // stop trial, no response by UC
  SS_TRIAL_GO_RESP,          // response on a go trial
  SS_TRIAL_STOP_GO_RESP,     // go response on a stop trial
  SS_TRIAL_STOP_ST_RESP      // ST (or unmatched) response on a stop trial
};

struct SsSharedState {
  bool valid = false;
  int n_trials = 0;                       // unique trials
  int n_acc = 0;                          // accumulators per trial
  std::vector<uint8_t> kind;              // SsTrialKind per trial
  std::vector<double> rt, SSD, UC, LC, LT;
  std::vector<int> resp_code;             // R code, NA_INTEGER if no response
  std::vector<int> n_go, n_st, first_go_row;   // first_go_row: offset in trial
  std::vector<uint8_t> is_go_row, is_st_row, win_row;  // n_trials * n_acc
  std::vector<int> lR_row;                // lR codes, n_trials * n_acc
};

static SsSharedState build_ss_shared_state(const Rcpp::DataFrame& data,
                                           int n_trials_ss, int n_acc) {
  SsSharedState ss;
  ss.n_trials = n_trials_ss;
  ss.n_acc = n_acc;
  if (n_trials_ss <= 0 || n_acc <= 0) return ss;
  if (!data.containsElementNamed("rt") || !data.containsElementNamed("R") ||
      !data.containsElementNamed("SSD") || !data.containsElementNamed("lR") ||
      !data.containsElementNamed("winner")) {
    return ss;
  }

  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  IntegerVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector LT = get_col_with_default(data, "LT", 0.0);
  NumericVector UC = get_col_with_default(data, "UC", R_PosInf);
  NumericVector LC = get_col_with_default(data, "LC", 0.0);
  const bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? Rcpp::as<IntegerVector>(data["lI"])
                            : IntegerVector(lR.size(), 2);
  if (RT.size() != n_trials_ss * n_acc) return ss;

  const int n_rows = n_trials_ss * n_acc;
  ss.kind.resize(n_trials_ss);
  ss.rt.resize(n_trials_ss); ss.SSD.resize(n_trials_ss);
  ss.UC.resize(n_trials_ss); ss.LC.resize(n_trials_ss); ss.LT.resize(n_trials_ss);
  ss.resp_code.resize(n_trials_ss);
  ss.n_go.resize(n_trials_ss); ss.n_st.resize(n_trials_ss);
  ss.first_go_row.resize(n_trials_ss);
  ss.is_go_row.resize(n_rows); ss.is_st_row.resize(n_rows); ss.win_row.resize(n_rows);
  ss.lR_row.resize(n_rows);

  for (int trial = 0; trial < n_trials_ss; ++trial) {
    const int start = trial * n_acc;
    const double rt = RT[start];
    const double ssd = SSD[start];
    const double uc = UC[start];
    const int r_obs = R[start];
    const bool response_observed = (r_obs != NA_INTEGER);
    const bool stop_presented = emc2_isfinite(ssd);
    const bool has_deadline = R_FINITE(uc) && !Rcpp::NumericVector::is_na(uc);

    ss.rt[trial] = rt; ss.SSD[trial] = ssd; ss.UC[trial] = uc;
    ss.LC[trial] = LC[start]; ss.LT[trial] = LT[start];
    ss.resp_code[trial] = r_obs;

    // go/ST partition exactly as before: within-trial max lI code marks go
    int go_code = lI[start];
    for (int i = 1; i < n_acc; ++i) {
      if (lI[start + i] > go_code) go_code = lI[start + i];
    }
    int ng = 0, nst = 0, first_go = -1;
    for (int i = 0; i < n_acc; ++i) {
      const bool go_i = has_lI ? (lI[start + i] == go_code) : true;
      ss.is_go_row[start + i] = go_i ? 1 : 0;
      ss.is_st_row[start + i] = go_i ? 0 : 1;
      ss.win_row[start + i] = (winner[start + i] == TRUE) ? 1 : 0;
      ss.lR_row[start + i] = lR[start + i];
      if (go_i) { ++ng; if (first_go < 0) first_go = i; }
      else ++nst;
    }
    ss.n_go[trial] = ng; ss.n_st[trial] = nst;
    ss.first_go_row[trial] = (first_go < 0) ? 0 : first_go;

    bool response_is_go = false;
    if (response_observed) {
      for (int i = 0; i < n_acc; ++i) {
        if (lR[start + i] == r_obs) { response_is_go = ss.is_go_row[start + i] != 0; break; }
      }
    }

    uint8_t k;
    if (rt == R_NegInf) {
      k = SS_TRIAL_LC_CENSORED;
    } else if (!response_observed) {
      if (!stop_presented) k = has_deadline ? SS_TRIAL_GO_NR_DEADLINE : SS_TRIAL_GO_NR_FREE;
      else                 k = has_deadline ? SS_TRIAL_STOP_NR_DEADLINE : SS_TRIAL_STOP_NR_FREE;
    } else if (!stop_presented) {
      k = SS_TRIAL_GO_RESP;
    } else {
      k = response_is_go ? SS_TRIAL_STOP_GO_RESP : SS_TRIAL_STOP_ST_RESP;
    }
    ss.kind[trial] = k;
  }

  ss.valid = true;
  return ss;
}

// Per-trial-per-particle evaluation context: canonical parameter blocks plus
// the data-fixed masks for this trial. All raw pointers; nothing allocates.
struct SsTrialEvalCtx {
  const SsRawModel* M;
  const double* acc;      // n_acc * SS_ACC_STRIDE canonical blocks (trial order)
  const double* acc_go;   // n_go * SS_ACC_STRIDE, go rows compacted in order
  const uint8_t* is_go;   // per trial row
  const uint8_t* is_st;
  const int* lR_codes;
  int n_acc, n_go, n_st;
  double stop_r0[4];      // stop block from trial row 0   (stop survivor)
  double stop_g0[4];      // stop block from first go row  (stop-success integral)
  double tf, gf;
  double SSD;
  bool stop_presented;
  int response_code;      // NA_INTEGER when unconstrained (LC integration)
  double min_ll;
};

static inline SsStopCtx ss_make_stop_ctx(const SsTrialEvalCtx& c) {
  SsStopCtx sc;
  sc.SSD = c.SSD;
  sc.acc_go = c.acc_go;
  sc.n_go = c.n_go;
  sc.muS = c.stop_g0[0]; sc.sigS = c.stop_g0[1];
  sc.tauS = c.stop_g0[2]; sc.lbS = c.stop_g0[3];
  sc.acc_surv = c.M->acc_surv;
  return sc;
}

// log survivor of the stop process at duration t (raw twin of
// stop_logsurv_texg_fn / stop_logsurv_rdex_fn: both models race a TEXG stop).
static inline double ss_stop_logsurv_raw(double t, const double* stop4) {
  return ptexg(t, stop4[0], stop4[1], stop4[2], stop4[3], R_PosInf, false, true);
}

// Sum of protected log survivors over a masked subset of accumulators, at
// time t. Matches the texg_go_lccdf/rdex_go_lccdf + sum pattern: each element
// is floored to -Inf when non-finite.
static inline double ss_masked_lsurv_sum_raw(const SsTrialEvalCtx& c, double t,
                                             const uint8_t* mask) {
  double out = 0.0;
  for (int i = 0; i < c.n_acc; ++i) {
    if (!mask[i]) continue;
    const double v = c.M->acc_lsurv(t, c.acc + SS_ACC_STRIDE * i);
    out += emc2_isfinite(v) ? v : R_NegInf;
  }
  return out;
}

// Full trial response density at rt, marginalizing over candidate winners.
// Raw twin of ss_trial_log_response_density (used by the LC-censoring path).
static double ss_trial_log_response_density_raw(double rt, const SsTrialEvalCtx& c) {
  const SsRawModel& M = *c.M;
  const bool response_known = (c.response_code != NA_INTEGER);

  double log_go_any = R_NegInf;
  for (int i = 0; i < c.n_acc; ++i) {
    if (!c.is_go[i]) continue;
    if (response_known && c.lR_codes[i] != c.response_code) continue;

    // winner-i go race density (twin of ss_log_go_density_for_winner)
    const double lw = M.acc_lpdf(rt, c.acc + SS_ACC_STRIDE * i);
    double go_lprob = emc2_isfinite(lw) ? lw : R_NegInf;
    if (!R_FINITE(go_lprob)) continue;
    for (int j = 0; j < c.n_acc; ++j) {
      if (!c.is_go[j] || j == i) continue;
      const double ls = M.acc_lsurv(rt, c.acc + SS_ACC_STRIDE * j);
      go_lprob += emc2_isfinite(ls) ? ls : R_NegInf;
    }
    if (!R_FINITE(go_lprob)) continue;

    double log_term = log1m(c.gf) + go_lprob;
    if (c.stop_presented && rt > c.SSD) {
      double rt_eff = rt - c.SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = ss_stop_logsurv_raw(rt_eff, c.stop_r0);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = c.min_ll;
      double st_loss_sum = 0.0;
      if (c.n_st > 0) {
        st_loss_sum = ss_masked_lsurv_sum_raw(c, rt_eff, c.is_st);
      }
      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      log_term = log1m(c.gf) + log_mix(c.tf, comp_tf, comp_notf);
    }
    log_go_any = log_sum_exp(log_go_any, log_term);
  }

  if (!c.stop_presented || c.n_st == 0) return log_go_any;

  double log_st_any = R_NegInf;
  for (int i = 0; i < c.n_acc; ++i) {
    if (!c.is_st[i]) continue;
    if (response_known && c.lR_codes[i] != c.response_code) continue;
    if (rt <= c.SSD) continue;

    double rt_st = rt - c.SSD;
    if (rt_st < 0.0) rt_st = 0.0;
    // winner-i ST race density (twin of ss_log_st_density_for_winner)
    const double lw = M.acc_lpdf(rt_st, c.acc + SS_ACC_STRIDE * i);
    double st_base = emc2_isfinite(lw) ? lw : R_NegInf;
    if (!R_FINITE(st_base)) continue;
    for (int j = 0; j < c.n_acc; ++j) {
      if (!c.is_st[j] || j == i) continue;
      const double ls = M.acc_lsurv(rt_st, c.acc + SS_ACC_STRIDE * j);
      st_base += emc2_isfinite(ls) ? ls : R_NegInf;
    }
    if (!R_FINITE(st_base)) continue;

    double go_loss_sum = 0.0;
    if (c.n_go > 0) {
      go_loss_sum = ss_masked_lsurv_sum_raw(c, rt, c.is_go);
    }
    double rt_eff = rt - c.SSD;
    if (rt_eff < 0.0) rt_eff = 0.0;
    const SsStopCtx sc = ss_make_stop_ctx(c);
    double log_pstop = ss_stop_success_raw_live(sc, c.min_ll, rt_eff,
                                                M.clamp_empty_window,
                                                M.try_analytic1);
    if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

    const double term_gf = std::log(c.gf) + st_base;
    const double term_stop_win = log1m(c.gf) + log_pstop + st_base;
    const double term_stop_lose =
      log1m(c.gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
    const double log_term =
      log1m(c.tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));

    log_st_any = log_sum_exp(log_st_any, log_term);
  }

  return log_sum_exp(log_go_any, log_st_any);
}

static double ss_lc_integrand_raw(double x, void* params) {
  const SsTrialEvalCtx* c = static_cast<const SsTrialEvalCtx*>(params);
  const double logf = ss_trial_log_response_density_raw(x, *c);
  if (!R_FINITE(logf)) return 0.0;
  return std::exp(logf);
}

// Raw twin of ss_integrate_lc_response_mass (same GSL setup/tolerances).
static inline double ss_integrate_lc_response_mass_raw(double lower, double upper,
                                                       const SsTrialEvalCtx& c) {
  if (!(upper > lower)) return c.min_ll;

  gsl_function F;
  F.function = &ss_lc_integrand_raw;
  F.params = const_cast<SsTrialEvalCtx*>(&c);

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, 200);
  double res = 0.0;
  double err = 0.0;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qags(&F, lower, upper, 1e-8, 1e-6, 200, workspace, &res, &err);
  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !R_FINITE(res) || res <= 0.0) return c.min_ll;
  return std::log(res);
}

struct SsRawWorkspace {
  std::vector<double> acc;     // n_acc * SS_ACC_STRIDE
  std::vector<double> acc_go;  // n_go  * SS_ACC_STRIDE
  std::vector<double> lls;     // per unique trial
};

// Raw-path per-particle stop-signal log likelihood. `cols` are pointers to the
// ParamTable base columns in p_types order (the same layout the materialized
// matrix had). Branch bodies mirror c_log_likelihood_ss trial by trial.
static double c_log_likelihood_ss_pt(
    const double* const* cols,
    const SsSharedState& ss,
    const SsRawModel& M,
    const LogicalVector& is_ok,
    const IntegerVector& expand,
    double min_ll,
    SsRawWorkspace& ws
) {
  const int n_out = expand.length();
  if (is_true(all(!is_ok))) {
    return static_cast<double>(n_out) * min_ll;
  }
  const int n_acc = ss.n_acc;
  const int n_trials = ss.n_trials;
  ws.acc.resize(static_cast<size_t>(SS_ACC_STRIDE) * n_acc);
  ws.acc_go.resize(static_cast<size_t>(SS_ACC_STRIDE) * n_acc);
  ws.lls.assign(static_cast<size_t>(n_trials), min_ll);

  for (int trial = 0; trial < n_trials; ++trial) {
    const int start = trial * n_acc;
    if (is_ok[start] != TRUE) { ws.lls[trial] = min_ll; continue; }

    // canonical parameter blocks for this trial/particle
    double* acc = ws.acc.data();
    for (int i = 0; i < n_acc; ++i) {
      M.fill_acc(cols, start + i, acc + SS_ACC_STRIDE * i);
    }
    double* acc_go = ws.acc_go.data();
    int k = 0;
    for (int i = 0; i < n_acc; ++i) {
      if (!ss.is_go_row[start + i]) continue;
      std::memcpy(acc_go + SS_ACC_STRIDE * k, acc + SS_ACC_STRIDE * i,
                  SS_ACC_STRIDE * sizeof(double));
      ++k;
    }

    SsTrialEvalCtx c;
    c.M = &M;
    c.acc = acc;
    c.acc_go = acc_go;
    c.is_go = ss.is_go_row.data() + start;
    c.is_st = ss.is_st_row.data() + start;
    c.lR_codes = ss.lR_row.data() + start;
    c.n_acc = n_acc;
    c.n_go = ss.n_go[trial];
    c.n_st = ss.n_st[trial];
    M.fill_stop(cols, start, c.stop_r0);
    M.fill_stop(cols, start + ss.first_go_row[trial], c.stop_g0);
    c.tf = cols[M.idx_tf][start];
    c.gf = cols[M.idx_gf][start];
    c.SSD = ss.SSD[trial];
    c.stop_presented = emc2_isfinite(c.SSD);
    c.response_code = ss.resp_code[trial];
    c.min_ll = min_ll;

    const uint8_t* win = ss.win_row.data() + start;
    const double rt = ss.rt[trial];
    double ll = min_ll;

    switch (ss.kind[trial]) {

    case SS_TRIAL_LC_CENSORED: {
      ll = ss_integrate_lc_response_mass_raw(ss.LT[trial], ss.LC[trial], c);
      break;
    }

    case SS_TRIAL_GO_NR_FREE: {
      ll = std::log(c.gf);
      break;
    }

    case SS_TRIAL_GO_NR_DEADLINE: {
      const double logS_go = (c.n_go > 0)
        ? ss_masked_lsurv_sum_raw(c, ss.UC[trial], c.is_go) : 0.0;
      ll = log_sum_exp(std::log(c.gf), log1m(c.gf) + logS_go);
      break;
    }

    case SS_TRIAL_STOP_NR_FREE: {
      if (c.n_st == 0) {
        const SsStopCtx sc = ss_make_stop_ctx(c);
        double log_pstop = ss_stop_success_raw_live(sc, min_ll, R_PosInf,
                                                    M.clamp_empty_window,
                                                    M.try_analytic1);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
        ll = log_sum_exp(std::log(c.gf),
                         log1m(c.gf) + log1m(c.tf) + log_pstop);
      } else {
        // NR with an ST accumulator and no deadline: tf AND gf must both fail
        ll = std::log(c.gf) + std::log(c.tf);
      }
      break;
    }

    case SS_TRIAL_STOP_NR_DEADLINE: {
      const double uc = ss.UC[trial];
      const double logS_go = (c.n_go > 0)
        ? ss_masked_lsurv_sum_raw(c, uc, c.is_go) : 0.0;

      double uc_eff = uc - c.SSD;
      if (!R_FINITE(uc_eff) || uc_eff <= 0.0) uc_eff = 0.0;
      const bool stop_can_act = R_FINITE(uc_eff) && (uc_eff > 0.0);

      double log_pstop = R_NegInf;
      double logS_stop = 0.0;   // log(1)
      if (stop_can_act) {
        logS_stop = ss_stop_logsurv_raw(uc_eff, c.stop_r0);
        const SsStopCtx sc = ss_make_stop_ctx(c);
        log_pstop = ss_stop_success_raw_live(sc, min_ll, uc_eff,
                                             M.clamp_empty_window,
                                             M.try_analytic1);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      }

      const double log_core_trig = log_sum_exp(log_pstop, logS_go + logS_stop);

      if (c.n_st == 0) {
        const double log_no_nogf = log_mix(c.tf, logS_go, log_core_trig);
        ll = log_sum_exp(std::log(c.gf), log1m(c.gf) + log_no_nogf);
      } else {
        const double logS_st = ss_masked_lsurv_sum_raw(c, uc_eff, c.is_st);
        const double log_trig =
          logS_st + log_sum_exp(std::log(c.gf), log1m(c.gf) + log_core_trig);
        const double log_tfbranch =
          log_sum_exp(std::log(c.gf), log1m(c.gf) + logS_go);
        ll = log_mix(c.tf, log_tfbranch, log_trig);
      }
      break;
    }

    case SS_TRIAL_GO_RESP:
    case SS_TRIAL_STOP_GO_RESP: {
      // go race density with the observed (data) winner
      double go_lprob = 0.0;
      bool any_win = false;
      for (int i = 0; i < n_acc; ++i) {
        if (!(win[i] && c.is_go[i])) continue;
        const double lw = M.acc_lpdf(rt, acc + SS_ACC_STRIDE * i);
        go_lprob += emc2_isfinite(lw) ? lw : R_NegInf;
        any_win = true;
      }
      if (!any_win || !R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (c.n_go > 1) {
        for (int i = 0; i < n_acc; ++i) {
          if (!(!win[i] && c.is_go[i])) continue;
          const double ls = M.acc_lsurv(rt, acc + SS_ACC_STRIDE * i);
          go_lprob += emc2_isfinite(ls) ? ls : R_NegInf;
        }
      }

      if (ss.kind[trial] == SS_TRIAL_GO_RESP) {
        ll = log1m(c.gf) + go_lprob;
        break;
      }

      // stop trial, go response: mixture over trigger failure
      double rt_eff = rt - c.SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = ss_stop_logsurv_raw(rt_eff, c.stop_r0);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;
      double st_loss_sum = 0.0;
      if (c.n_st > 0) {
        for (int i = 0; i < n_acc; ++i) {
          if (!(c.is_st[i] && !win[i])) continue;
          const double ls = M.acc_lsurv(rt_eff, acc + SS_ACC_STRIDE * i);
          st_loss_sum += emc2_isfinite(ls) ? ls : R_NegInf;
        }
      }
      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      ll = log1m(c.gf) + log_mix(c.tf, comp_tf, comp_notf);
      break;
    }

    case SS_TRIAL_STOP_ST_RESP: {
      double rt_st = rt - c.SSD;
      if (rt_st < 0.0) rt_st = 0.0;
      // ST winner log pdf at rt - SSD
      double st_winner_logpdf = 0.0;
      bool any_win = false;
      for (int i = 0; i < n_acc; ++i) {
        if (!(win[i] && c.is_st[i])) continue;
        const double lw = M.acc_lpdf(rt_st, acc + SS_ACC_STRIDE * i);
        st_winner_logpdf += emc2_isfinite(lw) ? lw : R_NegInf;
        any_win = true;
      }
      if (!any_win || !R_FINITE(st_winner_logpdf)) st_winner_logpdf = R_NegInf;
      // ST losers survivors
      double st_loss_sum = 0.0;
      if (c.n_st > 1) {
        for (int i = 0; i < n_acc; ++i) {
          if (!(!win[i] && c.is_st[i])) continue;
          const double ls = M.acc_lsurv(rt_st, acc + SS_ACC_STRIDE * i);
          st_loss_sum += emc2_isfinite(ls) ? ls : R_NegInf;
        }
      }
      // GO losers survivors (all go accumulators, at raw rt)
      double go_loss_sum = 0.0;
      if (c.n_go > 0) {
        go_loss_sum = ss_masked_lsurv_sum_raw(c, rt, c.is_go);
      }
      double rt_eff = rt - c.SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      const SsStopCtx sc = ss_make_stop_ctx(c);
      double log_pstop = ss_stop_success_raw_live(sc, min_ll, rt_eff,
                                                  M.clamp_empty_window,
                                                  M.try_analytic1);
      if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

      const double st_base = st_winner_logpdf + st_loss_sum;
      const double term_gf = std::log(c.gf) + st_base;
      const double term_stop_win = log1m(c.gf) + log_pstop + st_base;
      const double term_stop_lose =
        log1m(c.gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
      ll = log1m(c.tf) +
        log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
      break;
    }
    }

    // same clamp as the vectorized epilogue of c_log_likelihood_ss
    if (!R_FINITE(ll) || ll < min_ll) ll = min_ll;
    ws.lls[trial] = ll;
  }

  double sum_ll = 0.0;
  for (int i = 0; i < n_out; ++i) {
    sum_ll += ws.lls[expand[i] - 1];   // expand created in 1-based R
  }
  return sum_ll;
}

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
  key.s = s;         // Retain the raw scale in the exact key as requested.
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

// Resolve function pointers and fixed column indices for SS (EXG or RDEX).
struct SSModelAdapter {
  ss_go_pdf_fn    go_lpdf_ptr;
  ss_go_pdf_fn    go_lccdf_ptr;
  ss_stop_surv_fn stop_logsurv_ptr;
  ss_stop_success_fn stop_success_ptr;
  int idx_tf;
  int idx_gf;
};

static inline SSModelAdapter resolve_ss_adapter(const std::string& type_std) {
  SSModelAdapter a;
  const bool is_exg = type_std.find("EXG") != std::string::npos;
  a.go_lpdf_ptr      = is_exg ? texg_go_lpdf          : rdex_go_lpdf;
  a.go_lccdf_ptr     = is_exg ? texg_go_lccdf          : rdex_go_lccdf;
  a.stop_logsurv_ptr = is_exg ? stop_logsurv_texg_fn   : stop_logsurv_rdex_fn;
  a.stop_success_ptr = is_exg ? ss_texg_stop_success_lpdf_live : ss_rdex_stop_success_lpdf_live;
  a.idx_tf           = is_exg ? emc2col::ss_texg::tf : emc2col::ss_rdex::tf;
  a.idx_gf           = is_exg ? emc2col::ss_texg::gf : emc2col::ss_rdex::gf;
  return a;
}

static LogicalRulesSharedState build_logicalrules_shared_state(const Rcpp::DataFrame& dadm,
                                                               int n_trials,
                                                               int n_acc,
                                                               bool capacity) {
  LogicalRulesSharedState out;
  out.n_trials = n_trials;
  out.n_acc = n_acc;

  if (n_trials <= 0 || n_acc <= 0 || (n_trials % n_acc) != 0) {
    Rcpp::stop("LogicalRules likelihood expects n_trials > 0 and n_trials %% n_acc == 0.");
  }
  if (!dadm.containsElementNamed("rt") ||
      !dadm.containsElementNamed("lR") ||
      !dadm.containsElementNamed("R") ||
      !dadm.containsElementNamed("LogicalRule")) {
    Rcpp::stop("LogicalRules likelihood requires dadm columns: rt, lR, R, LogicalRule.");
  }

  Rcpp::NumericVector rts = dadm["rt"];
  SEXP role_sexp = dadm["lR"];
  SEXP resp_sexp = dadm["R"];
  SEXP rule_sexp = dadm["LogicalRule"];
  Rcpp::NumericVector LT = get_col_with_default(dadm, "LT", 0.0);
  Rcpp::NumericVector UT = get_col_with_default(dadm, "UT", R_PosInf);
  Rcpp::NumericVector LC = get_col_with_default(dadm, "LC", 0.0);
  Rcpp::NumericVector UC = get_col_with_default(dadm, "UC", R_PosInf);
  const bool has_stim_col = dadm.containsElementNamed("S") ||
                            dadm.containsElementNamed("stimulus") ||
                            dadm.containsElementNamed("condition");
  SEXP stim_sexp = R_NilValue;
  bool stim_is_factor = false;
  Rcpp::IntegerVector stim_code;
  Rcpp::CharacterVector stim_chr;
  int stimNN = -1, stimA = -1, stimB = -1, stimAB = -1;
  if (has_stim_col) {
    if (dadm.containsElementNamed("S")) stim_sexp = dadm["S"];
    else if (dadm.containsElementNamed("stimulus")) stim_sexp = dadm["stimulus"];
    else stim_sexp = dadm["condition"];
    stim_is_factor = Rf_inherits(stim_sexp, "factor");
  }
  if (rts.size() != n_trials ||
      Rf_length(role_sexp) != n_trials ||
      Rf_length(resp_sexp) != n_trials ||
      Rf_length(rule_sexp) != n_trials) {
    Rcpp::stop("LogicalRules likelihood: dadm column lengths must match n_trials.");
  }

  const bool role_is_factor = Rf_inherits(role_sexp, "factor");
  const bool resp_is_factor = Rf_inherits(resp_sexp, "factor");
  const bool rule_is_factor = Rf_inherits(rule_sexp, "factor");

  Rcpp::IntegerVector role_code;
  Rcpp::IntegerVector resp_code;
  Rcpp::IntegerVector rule_code;
  Rcpp::CharacterVector role_chr;
  Rcpp::CharacterVector resp_chr;
  Rcpp::CharacterVector rule_chr;

  auto level_code = [](const Rcpp::CharacterVector& levels, const std::string& target) -> int {
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == target) return i + 1;
    }
    return -1;
  };

  int codeA = -1, codeB = -1, codenA = -1, codenB = -1, codeNogo = -1;
  int codeYes = -1, codeNo = -1;
  int codeNN = -1, codeAN = -1, codeNB = -1, codeAB = -1;
  int codeAND = -1, codeOR = -1, codeXOR = -1, codeID = -1, codeORDet = -1, codeORGng = -1;

  if (role_is_factor) {
    role_code = Rcpp::IntegerVector(role_sexp);
    Rcpp::CharacterVector lev = role_code.attr("levels");
    codeA = level_code(lev, "A");
    codeB = level_code(lev, "B");
    codenA = level_code(lev, "n_A");
    codenB = level_code(lev, "n_B");
    codeNogo = level_code(lev, "nogo");
    if (codeA < 0 || codeB < 0) {
      Rcpp::stop("LogicalRules likelihood: lR levels must include A and B.");
    }
  } else {
    role_chr = Rcpp::CharacterVector(role_sexp);
  }

  if (resp_is_factor) {
    resp_code = Rcpp::IntegerVector(resp_sexp);
    Rcpp::CharacterVector lev = resp_code.attr("levels");
    codeYes = level_code(lev, "yes");
    codeNo = level_code(lev, "no");
    codeNN = level_code(lev, "NN");
    codeAN = level_code(lev, "AN");
    codeNB = level_code(lev, "NB");
    codeAB = level_code(lev, "AB");
    if (codeYes < 0 && codeNo < 0 && codeNN < 0 && codeAN < 0 && codeNB < 0 && codeAB < 0) {
      Rcpp::stop("LogicalRules likelihood: response levels must include binary (yes/no) and/or ID (NN/AN/NB/AB) codes.");
    }
  } else {
    resp_chr = Rcpp::CharacterVector(resp_sexp);
  }

  if (rule_is_factor) {
    rule_code = Rcpp::IntegerVector(rule_sexp);
    Rcpp::CharacterVector lev = rule_code.attr("levels");
    codeAND = level_code(lev, "AND");
    codeOR = level_code(lev, "OR");
    codeXOR = level_code(lev, "XOR");
    codeID = level_code(lev, "ID");
    codeORDet = level_code(lev, "OR_DETECTION_ANALYTIC");
    codeORGng = level_code(lev, "OR_DETECTION_GNG");
    if (codeAND < 0 && codeOR < 0 && codeXOR < 0 && codeID < 0 &&
        codeORDet < 0 && codeORGng < 0) {
      Rcpp::stop("LogicalRules likelihood: LogicalRule levels must include AND, OR, XOR, ID, OR_DETECTION_ANALYTIC, or OR_DETECTION_GNG.");
    }
  } else {
    rule_chr = Rcpp::CharacterVector(rule_sexp);
  }

  if (has_stim_col && stim_is_factor) {
    stim_code = Rcpp::IntegerVector(stim_sexp);
    Rcpp::CharacterVector lev = stim_code.attr("levels");
    stimNN = level_code(lev, "NN");
    stimA = level_code(lev, "AN");
    if (stimA < 0) stimA = level_code(lev, "A");
    stimB = level_code(lev, "NB");
    if (stimB < 0) stimB = level_code(lev, "B");
    stimAB = level_code(lev, "AB");
  } else if (has_stim_col) {
    stim_chr = Rcpp::CharacterVector(stim_sexp);
  }

  out.n_unique_trials = n_trials / n_acc;
  out.idxA.assign(out.n_unique_trials, -1);
  out.idxB.assign(out.n_unique_trials, -1);
  out.idxnA.assign(out.n_unique_trials, -1);
  out.idxnB.assign(out.n_unique_trials, -1);
  out.idxNogo.assign(out.n_unique_trials, -1);
  out.rule_code.assign(out.n_unique_trials, 0);
  out.cond_code.assign(out.n_unique_trials, 0);
  out.resp_code.assign(out.n_unique_trials, 0);
  out.rt_unique.assign(out.n_unique_trials, NA_REAL);
  out.LT_unique.assign(out.n_unique_trials, 0.0);
  out.UT_unique.assign(out.n_unique_trials, R_PosInf);
  out.LC_unique.assign(out.n_unique_trials, 0.0);
  out.UC_unique.assign(out.n_unique_trials, R_PosInf);
  out.rt_by_row.assign(rts.begin(), rts.end());

  for (int j = 0; j < out.n_unique_trials; ++j) {
    const int start = j * n_acc;
    out.rt_unique[static_cast<size_t>(j)] = rts[start];
    out.LT_unique[static_cast<size_t>(j)] = LT[start];
    out.UT_unique[static_cast<size_t>(j)] = UT[start];
    out.LC_unique[static_cast<size_t>(j)] = LC[start];
    out.UC_unique[static_cast<size_t>(j)] = UC[start];

    int rcode = 0;
    if (rule_is_factor) {
      const int rc = rule_code[start];
      if (rc == codeOR) rcode = 1;
      else if (rc == codeAND) rcode = 2;
      else if (rc == codeXOR) rcode = 3;
      else if (rc == codeID) rcode = 4;
      else if (rc == codeORDet) rcode = 5;
      else if (rc == codeORGng) rcode = 6;
    } else {
      const std::string rv = Rcpp::as<std::string>(rule_chr[start]);
      if (rv == "OR") rcode = 1;
      else if (rv == "AND") rcode = 2;
      else if (rv == "XOR") rcode = 3;
      else if (rv == "ID") rcode = 4;
      else if (rv == "OR_DETECTION_ANALYTIC") rcode = 5;
      else if (rv == "OR_DETECTION_GNG") rcode = 6;
    }
    if (rcode == 0) {
      Rcpp::stop("LogicalRules likelihood: each unique trial must have a supported rule.");
    }
    out.rule_code[static_cast<size_t>(j)] = rcode;

    int resp = 0;
    if (resp_is_factor) {
      const int rc = resp_code[start];
      if (rc == NA_INTEGER) resp = 0;
      else if (rc == codeYes) resp = 1;
      else if (rc == codeNo) resp = 2;
      else if (rc == codeNN) resp = 3;
      else if (rc == codeAN) resp = 4;
      else if (rc == codeNB) resp = 5;
      else if (rc == codeAB) resp = 6;
    } else {
      if (resp_chr[start] == NA_STRING) resp = 0;
      else {
        const std::string rv = Rcpp::as<std::string>(resp_chr[start]);
        if (rv == "yes") resp = 1;
        else if (rv == "no" || rv == "nogo") resp = 2;
        else if (rv == "NN") resp = 3;
        else if (rv == "AN") resp = 4;
        else if (rv == "NB") resp = 5;
        else if (rv == "AB") resp = 6;
      }
    }
    if (resp == 0 && rcode <= 4) {
      Rcpp::stop("LogicalRules likelihood: response must be yes/no/NN/AN/NB/AB on each unique trial.");
    }
    out.resp_code[static_cast<size_t>(j)] = resp;

    // Stimulus condition codes are mandatory for the detection rules and for
    // the capacity model (which needs to identify the redundant-target AB
    // condition for every rule).  Ordinary choice designs may carry an
    // unrelated S factor; their unmapped levels are recorded as -1 and never
    // consulted.
    // OR_DETECTION_GNG (rcode 6) is the full four-horse OR task with a withheld
    // (rt = +Inf) outcome in place of the overt "no"; it does not consult the
    // stimulus condition (the design drives the drifts, exactly as for OR), so
    // only the analytic detector (rcode 5) and the capacity model require it.
    int cond = (rcode == 5 || capacity) ? -1 : 0;
    if (rcode == 5 || capacity) {
      if (!has_stim_col) {
        Rcpp::stop(capacity
          ? "LogicalRules capacity requires stimulus column `S` (or `stimulus`/`condition`)."
          : "LogicalRules detection rules require stimulus column `S` (or `stimulus`/`condition`).");
      }
    }
    if (has_stim_col) {
      if (stim_is_factor) {
        const int sc = stim_code[start];
        if (sc == NA_INTEGER || sc == stimNN) cond = 0;
        else if (sc == stimA) cond = 1;
        else if (sc == stimB) cond = 2;
        else if (sc == stimAB) cond = 3;
        else cond = -1;
      } else {
        if (stim_chr[start] == NA_STRING) cond = 0;
        else {
          const std::string sv = Rcpp::as<std::string>(stim_chr[start]);
          if (sv == "NN" || sv == "none") cond = 0;
          else if (sv == "AN" || sv == "A") cond = 1;
          else if (sv == "NB" || sv == "B") cond = 2;
          else if (sv == "AB" || sv == "BA" || sv == "A+B" || sv == "B+A") cond = 3;
          else cond = -1;
        }
      }
    }
    if ((rcode == 5 || capacity) && (cond < 0 || cond > 3)) {
      Rcpp::stop(capacity
        ? "LogicalRules capacity requires stimulus NN/AN/NB/AB (or A/B/AB) on every trial."
        : "LogicalRules detection rules require stimulus NN/AN/NB/AB (or A/B/AB).");
    }
    out.cond_code[static_cast<size_t>(j)] = cond;

    for (int k = 0; k < n_acc; ++k) {
      const int idx = start + k;
      if (role_is_factor) {
        const int rc = role_code[idx];
        if (rc == codeA) out.idxA[static_cast<size_t>(j)] = idx;
        else if (rc == codeB) out.idxB[static_cast<size_t>(j)] = idx;
        else if (rc == codenA) out.idxnA[static_cast<size_t>(j)] = idx;
        else if (rc == codenB) out.idxnB[static_cast<size_t>(j)] = idx;
        else if (rc == codeNogo) out.idxNogo[static_cast<size_t>(j)] = idx;
      } else {
        const std::string rv = Rcpp::as<std::string>(role_chr[idx]);
        if (rv == "A") out.idxA[static_cast<size_t>(j)] = idx;
        else if (rv == "B") out.idxB[static_cast<size_t>(j)] = idx;
        else if (rv == "n_A") out.idxnA[static_cast<size_t>(j)] = idx;
        else if (rv == "n_B") out.idxnB[static_cast<size_t>(j)] = idx;
        else if (rv == "nogo") out.idxNogo[static_cast<size_t>(j)] = idx;
      }
    }
    if (out.idxNogo[static_cast<size_t>(j)] >= 0) out.has_nogo = true;

    if (rcode <= 4 || rcode == 6) {
      // OR/AND/XOR/ID and OR_DETECTION_GNG all use the four-horse subrace layout
      // (two target/nontarget channels).  GNG differs only in that the "no"
      // outcome is a withheld response (rt = +Inf) rather than an overt "no".
      if (out.idxA[static_cast<size_t>(j)] < 0 ||
          out.idxB[static_cast<size_t>(j)] < 0 ||
          out.idxnA[static_cast<size_t>(j)] < 0 ||
          out.idxnB[static_cast<size_t>(j)] < 0) {
        Rcpp::stop(rcode == 6
          ? "LogicalRules likelihood: OR_DETECTION_GNG trials require A, n_A, B, n_B accumulators."
          : "LogicalRules likelihood: OR/AND/XOR/ID trials require A, B, n_A, n_B accumulators.");
      }
    } else if (rcode == 5) {
      if (out.idxA[static_cast<size_t>(j)] < 0 || out.idxB[static_cast<size_t>(j)] < 0) {
        Rcpp::stop("LogicalRules likelihood: OR_DETECTION_ANALYTIC trials require A and B accumulators.");
      }
    }
  }

  const int n_unique_trials = out.n_unique_trials;
  out.cell_id.resize(n_unique_trials);
  std::vector<Rcpp::IntegerVector> dm_expands;
  if (dadm.hasAttribute("designs")) {
    Rcpp::List designs = dadm.attr("designs");
    for (int d = 0; d < designs.size(); ++d) {
      if (TYPEOF(designs[d]) == EXTPTRSXP) continue;
      Rcpp::RObject obj = designs[d];
      if (obj.hasAttribute("expand")) {
        dm_expands.push_back(obj.attr("expand"));
      }
    }
  }
  std::vector<int> acc_des_id(static_cast<size_t>(n_trials * n_acc), 0);
  std::vector<std::vector<int>> unique_profiles;
  for (int r = 0; r < n_trials * n_acc; ++r) {
    std::vector<int> profile;
    if (!dm_expands.empty()) {
      for (const auto& exp : dm_expands) {
        profile.push_back((r >= 0 && r < exp.size()) ? exp[r] : r);
      }
    } else {
      profile.push_back(r);
    }
    int found_p = -1;
    for (size_t u = 0; u < unique_profiles.size(); ++u) {
      if (unique_profiles[u] == profile) {
        found_p = static_cast<int>(u);
        break;
      }
    }
    if (found_p < 0) {
      found_p = static_cast<int>(unique_profiles.size());
      unique_profiles.push_back(profile);
    }
    acc_des_id[static_cast<size_t>(r)] = found_p;
  }
  auto get_des_idx = [&acc_des_id](int idx) -> int {
    if (idx >= 0 && idx < static_cast<int>(acc_des_id.size())) {
      return acc_des_id[static_cast<size_t>(idx)];
    }
    return idx;
  };
  for (int j = 0; j < n_unique_trials; ++j) {
    LogicalRulesCellKey k{
      get_des_idx(out.idxA[j]), get_des_idx(out.idxB[j]),
      get_des_idx(out.idxnA[j]), get_des_idx(out.idxnB[j]),
      get_des_idx(out.idxNogo[j]),
      out.rule_code[j], out.cond_code[j],
      out.LT_unique[j], out.UT_unique[j]
    };
    int found = -1;
    for (size_t c = 0; c < out.cells.size(); ++c) {
      if (out.cells[c] == k) {
        found = static_cast<int>(c);
        break;
      }
    }
    if (found < 0) {
      found = static_cast<int>(out.cells.size());
      out.cells.push_back(k);
    }
    out.cell_id[j] = found;
  }
  out.n_cells = static_cast<int>(out.cells.size());
  out.cell_log_z.assign(static_cast<size_t>(out.n_cells), NA_REAL);
  out.cell_log_den_cap.assign(static_cast<size_t>(out.n_cells), NA_REAL);

  out.valid = true;
  return out;
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

  Rcpp::LogicalVector prepare(int i) {
    if (i > 0) {
      table.fill_from_particle_row(particle_matrix_pt, i,
                                   pm_col_to_base_idx,
                                   invariant_base_idx_vec);
    }
    const bool skip_inv = (i > 0) && use_invariants;
    update_pt_only(table, designs, trend_runtime ? trend_runtime.get() : nullptr,
                   transform_specs,
                   skip_inv ? &invariant_design_mask : nullptr,
                   skip_inv ? &invariant_param_names : nullptr);
    if (i == 0) {
      bound_specs = make_bound_specs_pt(minmax, mm_names, table, bounds);
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

// Per-node grid of a t0 marginalization: the quadrature node values (sampled
// log-t0 scale) and the unnormalized per-particle log-terms
//   ell_ik = log quadrature_weight_ik + log p(x_ik|eta) + log L_i(x_ik)
// that BOTH the marginal-ll reducer and the storage-time weight accessor share,
// so the node grid and the reconstruction weights are guaranteed identical.
// Nodes are PER PARTICLE (np x K): each proposal gets a rule centred on its own
// conditional t0 posterior, so one particle's rule can never be dragged off its
// peak by the rest of a scattered batch.
struct MarginalGrid {
  Rcpp::NumericMatrix nodes;      // np x K: per-particle log-t0 node values
  Rcpp::NumericMatrix log_terms;  // np x K: ell_ik (softmax over k => p(t0=x_ik|y,theta,eta))
  Rcpp::NumericVector mode;       // np: fitted conditional mode (warm-start carrier)
  Rcpp::NumericVector sd;         // np: fitted Laplace scale  (warm-start carrier)
  double warm_used = 0.0;         // fraction of rows retaining the warm hint
  int pred_used = 0;              // 1 when the subset-regression hint replaced it
  int repaired = 0;               // particles whose rule missed and was refitted
};

// A composite GL rule lets the t0 marginal retain support over its complete
// feasible interval while concentrating most nodes around the response-informed
// posterior bump.  Unlike simply truncating to a centred window, the two tail
// panels retain the integral's mass and therefore keep the approximation
// coherent even for a broad or weakly identified t0 posterior.
struct MarginalRule {
  std::vector<double> x;
  std::vector<double> w;
};

static void marginal_append_gl_panel(MarginalRule& out, double lo, double hi,
                                     int n) {
  if (n <= 0 || !(hi > lo)) return;
  const GLRule& rule = gl_get_rule(n);
  const double half = 0.5 * (hi - lo);
  const double mid = 0.5 * (hi + lo);
  for (int k = 0; k < n; ++k) {
    out.x.push_back(mid + half * rule.x[static_cast<size_t>(k)]);
    out.w.push_back(half * rule.w[static_cast<size_t>(k)]);
  }
}

// One batched kernel pass: every particle is evaluated at its OWN column of
// node values through a single calc_ll_oo call (rows laid out r = k*np + i)
// rather than one call per node.  The per-call setup inside the kernel (design
// expansion, trial bookkeeping) is O(n_trials) and used to be re-paid for every
// node of every particle step.
static NumericMatrix marginal_eval_nodes(
    const NumericMatrix& particle_matrix, const NumericMatrix& X, int t0col,
    DataFrame data, NumericVector constants, List designs, String type,
    List bounds, List transforms, List pretransforms, CharacterVector p_types,
    double min_ll, Rcpp::Nullable<Rcpp::List> trend) {
  const int np = particle_matrix.nrow();
  const int n_par = particle_matrix.ncol();
  const int m = X.ncol();
  NumericMatrix out(np, m);
  if (np == 0 || m == 0) return out;
  CharacterVector p_names = colnames(particle_matrix);
  // Cap the rows per call so a wide grid on an uncompressed data set cannot
  // blow the expanded parameter matrix up.
  const int max_rows = 4096;
  int cols_per_call = m;
  if (np * m > max_rows) cols_per_call = std::max(1, max_rows / np);
  NumericMatrix pm(np * std::min(cols_per_call, m), n_par);
  colnames(pm) = p_names;
  for (int k0 = 0; k0 < m; k0 += cols_per_call) {
    const int mk = std::min(cols_per_call, m - k0);
    if (pm.nrow() != np * mk) {
      pm = NumericMatrix(np * mk, n_par);
      colnames(pm) = p_names;
    }
    for (int k = 0; k < mk; ++k) {
      for (int i = 0; i < np; ++i) {
        const int r = k * np + i;
        for (int j = 0; j < n_par; ++j) pm(r, j) = particle_matrix(i, j);
        pm(r, t0col) = X(i, k0 + k);
      }
    }
    NumericVector ll = calc_ll_oo(pm, data, constants, designs, type, bounds,
                                  transforms, pretransforms, p_types, min_ll,
                                  trend, R_NilValue);
    for (int k = 0; k < mk; ++k)
      for (int i = 0; i < np; ++i) out(i, k0 + k) = ll[k * np + i];
  }
  return out;
}

// Mode and scale of the parabola through three points of the log integrand
// g(x) = log L(x) + log p(x|eta).  Returns false unless the fit is a proper
// interior maximum, in which case the caller keeps the full feasible interval
// rather than a Laplace window.
static bool marginal_parabola(double x1, double g1, double x2, double g2,
                              double x3, double g3, double* mode, double* sd) {
  if (!(R_FINITE(g1) && R_FINITE(g2) && R_FINITE(g3))) return false;
  if (!(x1 < x2 && x2 < x3)) return false;
  const double d12 = (g2 - g1) / (x2 - x1);
  const double d23 = (g3 - g2) / (x3 - x2);
  const double curv = 2.0 * (d23 - d12) / (x3 - x1);   // ~ g''(mode)
  if (!R_FINITE(curv) || curv >= 0.0) return false;
  const double num = (x2 - x1) * (x2 - x1) * (g2 - g3) -
                     (x2 - x3) * (x2 - x3) * (g2 - g1);
  const double den = (x2 - x1) * (g2 - g3) - (x2 - x3) * (g2 - g1);
  if (!R_FINITE(den) || den == 0.0) return false;
  const double xm = x2 - 0.5 * num / den;
  if (!R_FINITE(xm)) return false;
  const double s = 1.0 / std::sqrt(-curv);
  if (!R_FINITE(s) || !(s > 0.0)) return false;
  *mode = xm;
  *sd = s;
  return true;
}

// Ridge-stabilised least squares of y on [1 X] (X is n x p, column-major), by
// Cholesky on the normal equations.  p is a handful of sampled parameters and n
// a few dozen probe particles, so this is negligible next to one kernel pass.
// Returns false when the system is not positive definite (rank-deficient design,
// e.g. a component whose columns are all constant across the batch).
static bool marginal_lsfit(const std::vector<double>& X, const std::vector<double>& y,
                           int n, int p, std::vector<double>& beta) {
  const int q = p + 1;
  if (n < q + 2) return false;
  std::vector<double> A(static_cast<size_t>(q) * q, 0.0), b(static_cast<size_t>(q), 0.0);
  auto xcol = [&](int j, int i) { return j == 0 ? 1.0 : X[static_cast<size_t>(j - 1) * n + i]; };
  for (int j = 0; j < q; ++j) {
    for (int k = j; k < q; ++k) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) sum += xcol(j, i) * xcol(k, i);
      A[static_cast<size_t>(j) * q + k] = A[static_cast<size_t>(k) * q + j] = sum;
    }
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += xcol(j, i) * y[static_cast<size_t>(i)];
    b[static_cast<size_t>(j)] = sum;
  }
  double tr = 0.0;
  for (int j = 0; j < q; ++j) tr += A[static_cast<size_t>(j) * q + j];
  const double ridge = 1e-8 * (tr / q) + 1e-12;
  for (int j = 0; j < q; ++j) A[static_cast<size_t>(j) * q + j] += ridge;
  // Cholesky, in place
  for (int j = 0; j < q; ++j) {
    double d = A[static_cast<size_t>(j) * q + j];
    for (int k = 0; k < j; ++k) d -= A[static_cast<size_t>(j) * q + k] * A[static_cast<size_t>(j) * q + k];
    if (!(d > 0.0) || !R_FINITE(d)) return false;
    d = std::sqrt(d);
    A[static_cast<size_t>(j) * q + j] = d;
    for (int i = j + 1; i < q; ++i) {
      double v = A[static_cast<size_t>(i) * q + j];
      for (int k = 0; k < j; ++k) v -= A[static_cast<size_t>(i) * q + k] * A[static_cast<size_t>(j) * q + k];
      A[static_cast<size_t>(i) * q + j] = v / d;
    }
  }
  beta.assign(static_cast<size_t>(q), 0.0);
  for (int i = 0; i < q; ++i) {            // forward
    double v = b[static_cast<size_t>(i)];
    for (int k = 0; k < i; ++k) v -= A[static_cast<size_t>(i) * q + k] * beta[static_cast<size_t>(k)];
    beta[static_cast<size_t>(i)] = v / A[static_cast<size_t>(i) * q + i];
  }
  for (int i = q - 1; i >= 0; --i) {       // back
    double v = beta[static_cast<size_t>(i)];
    for (int k = i + 1; k < q; ++k) v -= A[static_cast<size_t>(k) * q + i] * beta[static_cast<size_t>(k)];
    beta[static_cast<size_t>(i)] = v / A[static_cast<size_t>(i) * q + i];
  }
  for (int j = 0; j < q; ++j) if (!R_FINITE(beta[static_cast<size_t>(j)])) return false;
  return true;
}

// Coherent marginalization core (option 3): integrate ONE shared t0 out of the
// COMPLETE subject likelihood by a per-particle Laplace-centred composite
// Gauss-Legendre rule on the log-t0 (sampled) axis.  Given a warm start
// (`warm_mode`/`warm_sd`: the values fitted for this subject's accepted
// particle last iteration) a single probe round re-fits each particle at the
// right scale; otherwise a coarse pilot scan
// brackets each particle's own conditional mode, two parabolic refinements
// resolve its curvature (the conditional posterior narrows like 1/sqrt(n_trials)
// and can be orders of magnitude thinner than any prior-scaled window), and the
// core panel is then flanked by one tail node each side so the whole feasible
// interval stays integrated.  Reuses the registered kernel verbatim -- nodes
// overwrite the sampled t0 column and re-enter calc_ll_oo with marginalise off
// -- so it is model-agnostic (every race model, no per-model code).  The
// interval is clipped to where the likelihood is live: below t0's natural-scale
// lower bound c_do_bound floors the whole subject to n_trials*min_ll; above
// min(rt) the fastest response is infeasible.  Correctness reference is
// WorkingTests/marginal_t0_lib.R::marginal_ll_t0.
static MarginalGrid calc_ll_oo_marginal_core(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
    Rcpp::List marginalise) {
  const int np = particle_matrix.nrow();
  const int n_par = particle_matrix.ncol();

  // --- marginalization spec: p(t0|eta) = lognormal, log t0 ~ N(mu, sigma) ---
  const std::string mparam = Rcpp::as<std::string>(marginalise["param"]);
  const double mu    = Rcpp::as<double>(marginalise["mu"]);
  const double sigma = Rcpp::as<double>(marginalise["sigma"]);
  int    n_nodes = marginalise.containsElementNamed("n_nodes")
                     ? Rcpp::as<int>(marginalise["n_nodes"]) : 12;
  int    n_scan = marginalise.containsElementNamed("n_scan")
                     ? Rcpp::as<int>(marginalise["n_scan"]) : 7;
  int    n_refine = marginalise.containsElementNamed("n_refine")
                     ? Rcpp::as<int>(marginalise["n_refine"]) : 2;
  const bool adaptive = marginalise.containsElementNamed("adaptive")
                     ? Rcpp::as<bool>(marginalise["adaptive"]) : true;
  const double span = marginalise.containsElementNamed("span")
                     ? Rcpp::as<double>(marginalise["span"]) : 6.0;
  // Half-width of the core panel in Laplace SDs: 5 covers a Gaussian core to
  // ~1e-6 relative, and the flanking tail nodes carry the rest.
  const double window_sd = marginalise.containsElementNamed("window_sd")
                     ? Rcpp::as<double>(marginalise["window_sd"]) : 5.0;
  const double eps  = marginalise.containsElementNamed("eps")
                     ? Rcpp::as<double>(marginalise["eps"]) : 1e-3;
  // Warm start: the mode/scale fitted for this subject's accepted particle on
  // the previous iteration.  The conditional t0 posterior moves very little
  // between MCMC iterations, so a valid hint replaces the whole pilot scan (and
  // one refinement) with a single probe round at the right scale.
  const double warm_mode = marginalise.containsElementNamed("warm_mode")
                     ? Rcpp::as<double>(marginalise["warm_mode"]) : NA_REAL;
  const double warm_sd = marginalise.containsElementNamed("warm_sd")
                     ? Rcpp::as<double>(marginalise["warm_sd"]) : NA_REAL;
  // Probe at a few hint-SDs so the bracket still spans the peak when the batch
  // has drifted; the fit is re-scaled from those points anyway.
  const double warm_inflate = marginalise.containsElementNamed("warm_inflate")
                     ? Rcpp::as<double>(marginalise["warm_inflate"]) : 4.0;
  int n_refine_warm = marginalise.containsElementNamed("n_refine_warm")
                     ? Rcpp::as<int>(marginalise["n_refine_warm"]) : 2;
  // Within-iteration hint: pilot a 1/pred_frac subset, regress its modes on the
  // other sampled parameters, and bracket the rest of the batch from that.
  const bool predict_mode = marginalise.containsElementNamed("predict_mode")
                     ? Rcpp::as<bool>(marginalise["predict_mode"]) : true;
  int pred_frac = marginalise.containsElementNamed("pred_frac")
                     ? Rcpp::as<int>(marginalise["pred_frac"]) : 8;
  int n_refine_pred = marginalise.containsElementNamed("n_refine_pred")
                     ? Rcpp::as<int>(marginalise["n_refine_pred"]) : 2;
  if (n_nodes < 2) n_nodes = 2;
  if (n_scan < 3) n_scan = 3;
  if (n_refine < 0) n_refine = 0;
  if (n_refine_warm < 1) n_refine_warm = 1;
  if (n_refine_pred < 1) n_refine_pred = 1;
  if (pred_frac < 2) pred_frac = 2;

  // --- locate the sampled t0 coordinate in the proposal matrix ---
  CharacterVector p_names = colnames(particle_matrix);
  int t0col = -1;
  for (int j = 0; j < p_names.size(); ++j)
    if (Rcpp::as<std::string>(p_names[j]) == mparam) { t0col = j; break; }
  if (t0col < 0)
    Rcpp::stop("calc_ll_oo marginalise: '%s' is not a sampled column.", mparam);

  // --- lower clip: t0's natural-scale bound (minmax), on the log axis. ---
  double log_lo = R_NegInf;
  {
    NumericMatrix minmax = bounds["minmax"];
    CharacterVector mm_names = colnames(minmax);
    for (int j = 0; j < mm_names.size(); ++j)
      if (Rcpp::as<std::string>(mm_names[j]) == mparam) {
        const double lo = minmax(0, j);
        if (R_FINITE(lo) && lo > 0.0) log_lo = std::log(lo);
        break;
      }
  }
  // --- upper clip: feasibility ceiling log(min finite rt - eps). ---
  double log_ms = R_PosInf;
  {
    NumericVector rt = data["rt"];
    double mn = R_PosInf;
    for (int k = 0; k < rt.size(); ++k)
      if (R_FINITE(rt[k]) && rt[k] < mn) mn = rt[k];
    if (R_FINITE(mn)) {
      double top = mn - eps;
      if (top < std::numeric_limits<double>::min()) top = std::numeric_limits<double>::min();
      log_ms = std::log(top);
    }
  }

  const double lo = std::max(mu - span * sigma, log_lo);
  const double hi = std::min(mu + span * sigma, log_ms);

  MarginalGrid g;
  g.nodes = NumericMatrix(np, n_nodes);
  g.log_terms = NumericMatrix(np, n_nodes);
  g.mode = NumericVector(np, NA_REAL);
  g.sd = NumericVector(np, NA_REAL);
  if (!(hi > lo)) {                                  // collapsed interval => no mass
    std::fill(g.nodes.begin(), g.nodes.end(), NA_REAL);
    std::fill(g.log_terms.begin(), g.log_terms.end(), R_NegInf);
    return g;
  }
  if (np == 0) return g;

  // --- per-particle Laplace fit of g(x) = log L(x) + log p(x|eta) -----------
  const double def_mode = 0.5 * (lo + hi), def_sd = 0.25 * (hi - lo);
  std::vector<double> mode(static_cast<size_t>(np), def_mode);
  std::vector<double> sdev(static_cast<size_t>(np), def_sd);
  std::vector<char> fitted(static_cast<size_t>(np), 0);
  const double max_sd = 0.5 * (hi - lo);

  // One refinement round: probe g at (mode - h, mode, mode + h) with h the
  // current scale estimate and re-fit the parabola THERE.  A fit at the wrong
  // spacing badly mis-states the width of a razor-thin posterior (it shrinks
  // like 1/sqrt(n_trials)), so each round contracts the probe onto the peak.
  // Returns how many particles the round did not resolve CLEANLY: either the
  // parabola was not a proper maximum, or its vertex sits outside the middle of
  // the probe bracket, meaning the incoming centre/scale was an extrapolation
  // rather than a local fit.  The warm path uses that count to decide whether
  // its hint is still usable; the cold path ignores it (its own next round is
  // the correction).
  auto refine_on = [&](const NumericMatrix& PM, std::vector<double>& md,
                       std::vector<double>& sv, std::vector<char>& ft) -> int {
    const int nn = PM.nrow();
    NumericMatrix Xr(nn, 3);
    for (int i = 0; i < nn; ++i) {
      const double h = std::max(std::min(sv[static_cast<size_t>(i)], max_sd), 1e-8);
      double c = std::min(hi, std::max(lo, md[static_cast<size_t>(i)]));
      double a = c - h, b = c + h;
      if (a < lo) { a = lo; c = std::min(hi, a + h); b = std::min(hi, c + h); }
      if (b > hi) { b = hi; c = std::max(lo, b - h); a = std::max(lo, c - h); }
      Xr(i, 0) = a; Xr(i, 1) = c; Xr(i, 2) = b;
    }
    NumericMatrix gr = marginal_eval_nodes(PM, Xr, t0col, data, constants,
                                           designs, type, bounds, transforms,
                                           pretransforms, p_types, min_ll, trend);
    int bad = 0;
    for (int i = 0; i < nn; ++i) {
      double m_i, s_i;
      const double g1 = gr(i, 0) + R::dnorm(Xr(i, 0), mu, sigma, 1);
      const double g2 = gr(i, 1) + R::dnorm(Xr(i, 1), mu, sigma, 1);
      const double g3 = gr(i, 2) + R::dnorm(Xr(i, 2), mu, sigma, 1);
      ft[static_cast<size_t>(i)] = 0;
      if (marginal_parabola(Xr(i, 0), g1, Xr(i, 1), g2, Xr(i, 2), g3, &m_i, &s_i)) {
        // Let the window contract freely but grow only gradually, so one bad
        // fit cannot throw the rule back out to the pilot scale.
        const double cap = std::min(max_sd, 4.0 * sv[static_cast<size_t>(i)]);
        const double mc = std::min(hi, std::max(lo, m_i));
        const double dev = mc - Xr(i, 1);
        // A mode pinned on a boundary can never sit at the probe centre (the
        // bracket is shifted inwards to stay feasible), so a flush bracket
        // pointing further out is consistent, not a failed fit.
        const bool ok_local = std::fabs(dev) <= 0.5 * (Xr(i, 2) - Xr(i, 1)) ||
                              (Xr(i, 0) <= lo && dev < 0.0) ||
                              (Xr(i, 2) >= hi && dev > 0.0);
        if (!ok_local) ++bad;
        else ft[static_cast<size_t>(i)] = 1;
        md[static_cast<size_t>(i)] = mc;
        sv[static_cast<size_t>(i)] = std::min(cap, s_i);
      } else {
        ++bad;
      }
    }
    return bad;
  };
  auto refine_round = [&]() -> int { return refine_on(particle_matrix, mode, sdev, fitted); };

  // Coarse uniform scan of the whole feasible interval: brackets each particle's
  // own mode.  The log integrand is smooth and unimodal in log t0, so the
  // bracket is reliable even when the exponentiated peak is far thinner than the
  // scan spacing -- but it costs n_scan full passes over the batch, which is why
  // it is worth running on a subset and predicting the rest.
  auto pilot_on = [&](const NumericMatrix& PM, std::vector<double>& md,
                      std::vector<double>& sv, std::vector<char>& ft) {
    const int nn = PM.nrow();
    NumericMatrix Xs(nn, n_scan);
    std::vector<double> xs(static_cast<size_t>(n_scan));
    for (int k = 0; k < n_scan; ++k) {
      xs[static_cast<size_t>(k)] = lo + (hi - lo) * k / (n_scan - 1);
      for (int i = 0; i < nn; ++i) Xs(i, k) = xs[static_cast<size_t>(k)];
    }
    NumericMatrix gs = marginal_eval_nodes(PM, Xs, t0col, data, constants,
                                           designs, type, bounds, transforms,
                                           pretransforms, p_types, min_ll, trend);
    for (int k = 0; k < n_scan; ++k) {
      const double lp = R::dnorm(xs[static_cast<size_t>(k)], mu, sigma, 1);
      for (int i = 0; i < nn; ++i) gs(i, k) += lp;
    }
    for (int i = 0; i < nn; ++i) {
      int best = 0;
      double bv = R_NegInf;
      for (int k = 0; k < n_scan; ++k) if (gs(i, k) > bv) { bv = gs(i, k); best = k; }
      const int j = std::min(std::max(best, 1), n_scan - 2);
      double m_i, s_i;
      if (marginal_parabola(xs[static_cast<size_t>(j - 1)], gs(i, j - 1),
                            xs[static_cast<size_t>(j)],     gs(i, j),
                            xs[static_cast<size_t>(j + 1)], gs(i, j + 1), &m_i, &s_i)) {
        ft[static_cast<size_t>(i)] = 1;
        md[static_cast<size_t>(i)] = std::min(hi, std::max(lo, m_i));
        sv[static_cast<size_t>(i)] = std::min(max_sd, s_i);
      }
    }
  };

  // Particles whose window came from a HINT rather than from the pilot scan:
  // only those can be wrong in a way the pilot would fix, so only those are
  // worth verifying against the finished grid.
  std::vector<char> hinted(static_cast<size_t>(np), 0);
  if (adaptive && n_nodes >= 6) {
    bool warm_ok = R_FINITE(warm_mode) && R_FINITE(warm_sd) && warm_sd > 0.0 &&
                   warm_mode >= lo && warm_mode <= hi;
    if (warm_ok) {
      const double h0 = std::min(max_sd, std::max(warm_sd * warm_inflate, 1e-8));
      for (int i = 0; i < np; ++i) {
        mode[static_cast<size_t>(i)] = warm_mode;
        sdev[static_cast<size_t>(i)] = h0;
        fitted[static_cast<size_t>(i)] = 1;
      }
      // One probe round is enough when every particle's mode sits inside the
      // hinted bracket. Otherwise give the batch a second round at the rescaled
      // spacing -- unless the hint was so far off that most particles missed,
      // in which case chasing it costs more than the pilot it replaces.
      int bad = refine_round();
      if (bad > 0 && bad * 4 <= np) {
        for (int r = 1; r < n_refine_warm && bad > 0; ++r) bad = refine_round();
      }
      g.warm_used = static_cast<double>(np - bad) / static_cast<double>(np);
      if (bad == 0) {
        std::fill(hinted.begin(), hinted.end(), 1);
      } else if (bad * 4 <= np) {
        std::vector<int> bad_rows;
        for (int i = 0; i < np; ++i) {
          if (fitted[static_cast<size_t>(i)]) hinted[static_cast<size_t>(i)] = 1;
          else bad_rows.push_back(i);
        }
        const int nb = static_cast<int>(bad_rows.size());
        NumericMatrix Pb(nb, n_par);
        colnames(Pb) = p_names;
        for (int j = 0; j < nb; ++j)
          for (int c = 0; c < n_par; ++c) Pb(j, c) = particle_matrix(bad_rows[static_cast<size_t>(j)], c);
        std::vector<double> mb(static_cast<size_t>(nb), def_mode);
        std::vector<double> sb(static_cast<size_t>(nb), def_sd);
        std::vector<char> fb(static_cast<size_t>(nb), 0);
        pilot_on(Pb, mb, sb, fb);
        refine_on(Pb, mb, sb, fb);
        for (int j = 0; j < nb; ++j) {
          const int i = bad_rows[static_cast<size_t>(j)];
          mode[static_cast<size_t>(i)] = mb[static_cast<size_t>(j)];
          sdev[static_cast<size_t>(i)] = sb[static_cast<size_t>(j)];
          fitted[static_cast<size_t>(i)] = fb[static_cast<size_t>(j)];
        }
      } else {
        warm_ok = false;
        g.warm_used = 0;
        std::fill(mode.begin(), mode.end(), def_mode);
        std::fill(sdev.begin(), sdev.end(), def_sd);
        std::fill(fitted.begin(), fitted.end(), 0);
      }
    }
    if (!warm_ok) {
      // Within-iteration hint: the conditional t0 mode is a smooth function of
      // the OTHER sampled parameters, so pilot-scan a small spread-out subset of
      // the batch, regress its modes on those parameters, and use the fit to
      // bracket every remaining particle.  Unlike the cross-iteration warm start
      // this hint is per particle and never stale -- it tracks the batch it was
      // built from -- and it replaces n_scan full passes with n_scan*k/np.
      bool pred_ok = false;
      if (predict_mode) {
        // Predictors: every sampled column that actually varies here (a constant
        // column, e.g. another component's parameters, carries no information
        // and would make the design rank-deficient).
        std::vector<int> cols;
        for (int j = 0; j < n_par; ++j) {
          if (j == t0col) continue;
          double mn = particle_matrix(0, j), mx = mn;
          for (int i = 1; i < np; ++i) {
            const double v = particle_matrix(i, j);
            if (v < mn) mn = v; else if (v > mx) mx = v;
          }
          if (mx - mn > 0.0) cols.push_back(j);
        }
        const int p_pred = static_cast<int>(cols.size());
        const int k_sub = std::max(p_pred + 6, np / pred_frac);
        // Only worth it when the subset is a small fraction of the batch.
        if (p_pred > 0 && np >= 2 * k_sub) {
          NumericMatrix Ps(k_sub, n_par);
          colnames(Ps) = p_names;
          std::vector<int> S(static_cast<size_t>(k_sub));
          for (int j = 0; j < k_sub; ++j) {
            // Evenly spaced: the batch is laid out in proposal-type blocks, so
            // a spread subset spans all of them.
            S[static_cast<size_t>(j)] = static_cast<int>(
              (static_cast<double>(j) * np) / k_sub);
            for (int c = 0; c < n_par; ++c)
              Ps(j, c) = particle_matrix(S[static_cast<size_t>(j)], c);
          }
          std::vector<double> ms(static_cast<size_t>(k_sub), def_mode);
          std::vector<double> ss(static_cast<size_t>(k_sub), def_sd);
          std::vector<char> fs(static_cast<size_t>(k_sub), 0);
          pilot_on(Ps, ms, ss, fs);
          refine_on(Ps, ms, ss, fs);
          // Regress the resolved subset modes on their parameters.
          std::vector<double> yv;
          std::vector<int> rows;
          for (int j = 0; j < k_sub; ++j)
            if (fs[static_cast<size_t>(j)]) { rows.push_back(j); yv.push_back(ms[static_cast<size_t>(j)]); }
          const int n_ok = static_cast<int>(rows.size());
          std::vector<double> Xv(static_cast<size_t>(n_ok) * p_pred);
          for (int c = 0; c < p_pred; ++c)
            for (int i = 0; i < n_ok; ++i)
              Xv[static_cast<size_t>(c) * n_ok + i] = Ps(rows[static_cast<size_t>(i)], cols[static_cast<size_t>(c)]);
          std::vector<double> beta;
          if (marginal_lsfit(Xv, yv, n_ok, p_pred, beta)) {
            // Residual scale sets the probe half-width: the bracket has to be
            // wide enough to contain the prediction error, and the refinement
            // rounds contract from there onto each particle's own peak.
            // Robust scale: the proposal batch is a MIXTURE, and its few very
            // wide draws are outliers to the local linear fit.  A plain RMS lets
            // them set the probe width for the entire batch (and then the whole
            // hint gets rejected as uninformative); a median absolute residual
            // sizes the probe for the bulk, and the verification pass below
            // catches the handful the regression could not place.
            std::vector<double> absres(static_cast<size_t>(n_ok));
            for (int i = 0; i < n_ok; ++i) {
              double f = beta[0];
              for (int c = 0; c < p_pred; ++c) f += beta[static_cast<size_t>(c + 1)] * Xv[static_cast<size_t>(c) * n_ok + i];
              absres[static_cast<size_t>(i)] = std::fabs(yv[static_cast<size_t>(i)] - f);
            }
            std::sort(absres.begin(), absres.end());
            const double rms = 1.4826 * absres[absres.size() / 2];
            std::vector<double> sorted;
            for (int j = 0; j < k_sub; ++j) if (fs[static_cast<size_t>(j)]) sorted.push_back(ss[static_cast<size_t>(j)]);
            std::sort(sorted.begin(), sorted.end());
            const double med_sd = sorted.empty() ? def_sd : sorted[sorted.size() / 2];
            const double h0 = std::min(max_sd, std::max(std::max(2.5 * rms, 4.0 * med_sd), 1e-8));
            // A prediction no sharper than the interval itself is worthless.
            if (R_FINITE(rms) && h0 < 0.5 * (hi - lo)) {
              for (int i = 0; i < np; ++i) {
                double f = beta[0];
                for (int c = 0; c < p_pred; ++c)
                  f += beta[static_cast<size_t>(c + 1)] * particle_matrix(i, cols[static_cast<size_t>(c)]);
                if (!R_FINITE(f)) { f = def_mode; }
                mode[static_cast<size_t>(i)] = std::min(hi, std::max(lo, f));
                sdev[static_cast<size_t>(i)] = h0;
                fitted[static_cast<size_t>(i)] = 1;
              }
              std::fill(hinted.begin(), hinted.end(), 1);
              // Keep the subset's own resolved fits: they are exact, not predicted.
              for (int j = 0; j < k_sub; ++j) if (fs[static_cast<size_t>(j)]) {
                hinted[static_cast<size_t>(S[static_cast<size_t>(j)])] = 0;
                const int i = S[static_cast<size_t>(j)];
                mode[static_cast<size_t>(i)] = ms[static_cast<size_t>(j)];
                sdev[static_cast<size_t>(i)] = ss[static_cast<size_t>(j)];
              }
              pred_ok = true;
              g.pred_used = 1;
            }
          }
        }
      }
      int bad = 0;
      if (pred_ok) {
        for (int r = 0; r < n_refine_pred; ++r)
          bad = refine_on(particle_matrix, mode, sdev, fitted);
        // The prediction is only a bracket; if it failed to bracket a
        // meaningful share of the batch, fall back to the full pilot.
        if (bad * 4 > np) {
          pred_ok = false;
          g.pred_used = 0;
          std::fill(hinted.begin(), hinted.end(), 0);
          std::fill(mode.begin(), mode.end(), def_mode);
          std::fill(sdev.begin(), sdev.end(), def_sd);
          std::fill(fitted.begin(), fitted.end(), 0);
        }
      }
      if (!pred_ok) {
        std::fill(hinted.begin(), hinted.end(), 0);
        pilot_on(particle_matrix, mode, sdev, fitted);
        for (int r = 0; r < n_refine; ++r) refine_round();
      }
    }
  }

  // --- per-particle composite rule: [lo, core) core [core, hi] --------------
  // Records whether the rule actually spends a node on each tail, which is what
  // makes the finished grid self-checking below.
  std::vector<char> has_left(static_cast<size_t>(np), 0), has_right(static_cast<size_t>(np), 0);
  std::vector<int> core_lo(static_cast<size_t>(np), 0), core_hi(static_cast<size_t>(np), n_nodes - 1);
  NumericMatrix logw(np, n_nodes);   // quadrature weights, so g = log_terms - logw
  auto build_rules = [&](const std::vector<int>& idxs) {
    for (size_t t = 0; t < idxs.size(); ++t) {
      const int i = idxs[t];
      double clo = lo, chi = hi;
      if (fitted[static_cast<size_t>(i)]) {
        const double half = window_sd * sdev[static_cast<size_t>(i)];
        clo = std::max(lo, mode[static_cast<size_t>(i)] - half);
        chi = std::min(hi, mode[static_cast<size_t>(i)] + half);
        if (!(chi > clo)) { clo = lo; chi = hi; }
        g.mode[i] = mode[static_cast<size_t>(i)];
        g.sd[i] = sdev[static_cast<size_t>(i)];
      } else {
        g.mode[i] = NA_REAL;
        g.sd[i] = NA_REAL;
      }
      const int n_left  = (clo > lo && n_nodes >= 4) ? 1 : 0;
      const int n_right = (chi < hi && n_nodes >= 4) ? 1 : 0;
      if (n_left == 0) clo = lo;
      if (n_right == 0) chi = hi;
      has_left[static_cast<size_t>(i)] = static_cast<char>(n_left);
      has_right[static_cast<size_t>(i)] = static_cast<char>(n_right);
      core_lo[static_cast<size_t>(i)] = n_left;
      core_hi[static_cast<size_t>(i)] = n_nodes - 1 - n_right;
      MarginalRule rule;
      marginal_append_gl_panel(rule, lo, clo, n_left);
      marginal_append_gl_panel(rule, clo, chi, n_nodes - n_left - n_right);
      marginal_append_gl_panel(rule, chi, hi, n_right);
      for (int k = 0; k < n_nodes; ++k) {
        const double xk = rule.x[static_cast<size_t>(k)];
        g.nodes(i, k) = xk;
        logw(i, k) = std::log(rule.w[static_cast<size_t>(k)]);
        g.log_terms(i, k) = logw(i, k) + R::dnorm(xk, mu, sigma, 1);
      }
    }
  };
  // Evaluate the kernel on a set of rows at their own nodes and add it in.
  auto eval_rows = [&](const std::vector<int>& idxs) {
    const int nb = static_cast<int>(idxs.size());
    if (nb == 0) return;
    NumericMatrix Pb(nb, n_par), Xb(nb, n_nodes);
    colnames(Pb) = p_names;
    for (int j = 0; j < nb; ++j) {
      for (int c = 0; c < n_par; ++c) Pb(j, c) = particle_matrix(idxs[static_cast<size_t>(j)], c);
      for (int k = 0; k < n_nodes; ++k) Xb(j, k) = g.nodes(idxs[static_cast<size_t>(j)], k);
    }
    NumericMatrix lb = marginal_eval_nodes(Pb, Xb, t0col, data, constants, designs, type,
                                           bounds, transforms, pretransforms, p_types,
                                           min_ll, trend);
    for (int j = 0; j < nb; ++j)
      for (int k = 0; k < n_nodes; ++k)
        g.log_terms(idxs[static_cast<size_t>(j)], k) += lb(j, k);
  };

  std::vector<int> all_rows(static_cast<size_t>(np));
  for (int i = 0; i < np; ++i) all_rows[static_cast<size_t>(i)] = i;
  build_rules(all_rows);

  // --- one batched pass over the whole np x K grid --------------------------
  NumericMatrix llf = marginal_eval_nodes(particle_matrix, g.nodes, t0col, data, constants,
                                          designs, type, bounds, transforms,
                                          pretransforms, p_types, min_ll, trend);
  for (int k = 0; k < n_nodes; ++k)
    for (int i = 0; i < np; ++i) g.log_terms(i, k) += llf(i, k);

  // --- verify the finished rule, and repair whatever it missed --------------
  // The mode search can be fooled: a hinted bracket that never contained the
  // peak, or a parabola fitted on a shoulder, yields a window that is either
  // off the peak or far too wide for it.  The finished grid detects both for
  // free.  Reading the integrand back off its own nodes (g = log_terms - log w),
  // a rule that resolved its peak has an interior maximum whose curvature
  // implies a width comparable to the node spacing.  A maximum at an end node,
  // or an implied width far BELOW the spacing (the tell-tale of a wide window
  // straddling a razor-thin peak, which silently loses most of the mass), means
  // the rule is not integrating what it should.  Those particles -- and only
  // those -- are refitted from the full pilot scan and re-integrated, so a hint
  // can cost time but never accuracy.
  if (adaptive && n_nodes >= 6) {
    std::vector<int> hurt;
    for (int i = 0; i < np; ++i) {
      // A window that already came from the pilot scan has nothing to fall back
      // to: refitting it would repeat the identical deterministic search.
      if (!hinted[static_cast<size_t>(i)]) continue;
      int best = 0;
      double bv = R_NegInf;
      for (int k = 0; k < n_nodes; ++k) if (g.log_terms(i, k) > bv) { bv = g.log_terms(i, k); best = k; }
      if (!R_FINITE(bv)) continue;
      const int clo_i = core_lo[static_cast<size_t>(i)], chi_i = core_hi[static_cast<size_t>(i)];
      // Peak at an end node, or outside the core panel: the window is off-peak.
      if (best <= clo_i || best >= chi_i) { hurt.push_back(i); continue; }
      const double x1 = g.nodes(i, best - 1), x2 = g.nodes(i, best), x3 = g.nodes(i, best + 1);
      double m_i, s_i;
      if (!marginal_parabola(x1, g.log_terms(i, best - 1) - logw(i, best - 1),
                             x2, g.log_terms(i, best) - logw(i, best),
                             x3, g.log_terms(i, best + 1) - logw(i, best + 1), &m_i, &s_i)) {
        hurt.push_back(i);
        continue;
      }
      // Gauss-Legendre nodes cluster towards the panel edges, so on a CORRECT
      // rule (core = mode +/- window_sd sigma) the spacing either side of the
      // peak is already a couple of sigma; only a width far below that means the
      // nodes are straddling a peak they never resolved.
      if (s_i < 0.2 * 0.5 * (x3 - x1)) hurt.push_back(i);   // peak under-resolved
    }
    g.repaired = static_cast<int>(hurt.size());
    if (!hurt.empty()) {
      const int nb = g.repaired;
      NumericMatrix Pb(nb, n_par);
      colnames(Pb) = p_names;
      for (int j = 0; j < nb; ++j)
        for (int c = 0; c < n_par; ++c) Pb(j, c) = particle_matrix(hurt[static_cast<size_t>(j)], c);
      std::vector<double> mb(static_cast<size_t>(nb), def_mode);
      std::vector<double> sb(static_cast<size_t>(nb), def_sd);
      std::vector<char> fb(static_cast<size_t>(nb), 0);
      pilot_on(Pb, mb, sb, fb);
      for (int r = 0; r < n_refine; ++r) refine_on(Pb, mb, sb, fb);
      for (int j = 0; j < nb; ++j) {
        const int i = hurt[static_cast<size_t>(j)];
        mode[static_cast<size_t>(i)] = mb[static_cast<size_t>(j)];
        sdev[static_cast<size_t>(i)] = sb[static_cast<size_t>(j)];
        fitted[static_cast<size_t>(i)] = fb[static_cast<size_t>(j)];
      }
      build_rules(hurt);
      eval_rows(hurt);
    }
  }
  return g;
}

// Marginal log-likelihood: log-sum-exp of the shared core's per-node terms.
// (Not Rcpp-exported: internal helper, reached only via calc_ll_oo.)
static NumericVector calc_ll_oo_marginal(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
    Rcpp::List marginalise) {
  MarginalGrid g = calc_ll_oo_marginal_core(particle_matrix, data, constants, designs, type,
      bounds, transforms, pretransforms, p_types, min_ll, trend, marginalise);
  const int np = g.log_terms.nrow(), nn = g.log_terms.ncol();
  NumericVector out(np);
  for (int i = 0; i < np; ++i) {
    double m = R_NegInf;
    for (int k = 0; k < nn; ++k) if (g.log_terms(i, k) > m) m = g.log_terms(i, k);
    if (!R_FINITE(m)) { out[i] = R_NegInf; continue; }
    double s = 0.0;
    for (int k = 0; k < nn; ++k) s += std::exp(g.log_terms(i, k) - m);
    out[i] = m + std::log(s);
  }
  return out;
}

// Storage-time accessor for the t0 reconstruction (reconstruct-at-storage design,
// plan "use the posterior as sampled"). Returns the quadrature node grid and the
// unnormalized per-particle log-terms; the sampler forms w = softmax_k(log_terms)
// and draws one node per subject to write into the stored alpha[t0], so the
// posterior carries a valid t0 draw and predict()/make_data() stay unchanged.
// Nodes are on the sampled (log-t0) scale, matching alpha.
// [[Rcpp::export]]
List calc_ll_oo_marginal_nodes(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::List marginalise,
    Rcpp::Nullable<Rcpp::List> trend = R_NilValue) {
  MarginalGrid g = calc_ll_oo_marginal_core(particle_matrix, data, constants, designs, type,
      bounds, transforms, pretransforms, p_types, min_ll, trend, marginalise);
  const int np = g.log_terms.nrow(), nn = g.log_terms.ncol();
  NumericVector ll(np);
  for (int i = 0; i < np; ++i) {
    double m = R_NegInf;
    for (int k = 0; k < nn; ++k) if (g.log_terms(i, k) > m) m = g.log_terms(i, k);
    if (!R_FINITE(m)) { ll[i] = R_NegInf; continue; }
    double s = 0.0;
    for (int k = 0; k < nn; ++k) s += std::exp(g.log_terms(i, k) - m);
    ll[i] = m + std::log(s);
  }
  return List::create(Rcpp::Named("nodes") = g.nodes,
                      Rcpp::Named("log_terms") = g.log_terms,
                      Rcpp::Named("ll") = ll,
                      Rcpp::Named("mode") = g.mode,
                      Rcpp::Named("sd") = g.sd,
                      Rcpp::Named("warm_used") = g.warm_used,
                      Rcpp::Named("pred_used") = g.pred_used,
                      Rcpp::Named("repaired") = g.repaired);
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
    configure_bawl_corr_context(adapter, keep_names, "calc_ll_oo");
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
      !adapter.ctx.bawl_correlated && !adapter.ctx.rdmswtn_correlated &&
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
    BAwLCorrSharedState bawl_corr_shared;
    if (adapter.ctx.bawl_correlated) {
      bawl_corr_shared = build_bawl_corr_shared_state(
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
      } else if (adapter.ctx.bawl_correlated) {
        lls[i] = c_log_likelihood_bawl_correlated(
            bawl_corr_shared, data, adapter.pdf1_ptr, adapter.cdf1_ptr,
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
    configure_bawl_corr_context(adapter, keep_names, "calc_ll_oo_pw");
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

    BAwLCorrSharedState bawl_corr_shared;
    if (adapter.ctx.bawl_correlated) {
      bawl_corr_shared = build_bawl_corr_shared_state(
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
      if (adapter.ctx.bawl_correlated) {
        c_log_likelihood_bawl_correlated(
            bawl_corr_shared, data, adapter.pdf1_ptr, adapter.cdf1_ptr,
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
  // to t0 removes the dead zone entirely, and is strictly cheaper.  t0_index is
  // -1 for models whose t0 is not a pure shift (e.g. REXG), which skips this.
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
    // logS_UT == R_NegInf for UT==Inf, defective or not). Upper-truncating a
    // defective distribution is conceptually ill-posed (an omission is a
    // never-finish outcome, not a slow response that got cut off) and is
    // intentionally not supported, so no S(Inf) is subtracted here.
    return logS_LT;
  }
  
  const double logS_UT = log_survivor_rowmajor(UT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
  double logP = log_diff_exp(logS_LT, logS_UT);
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
  if (R_FINITE(log_total) && log_total > log_prob_eps) return log_total;
  return R_NegInf;
}

static inline void copy_par_row_colmajor(const double* const* cols,
                                         int n_par,
                                         int row_idx,
                                         double* out_row) {
  for (int c = 0; c < n_par; ++c) {
    out_row[c] = cols[c][row_idx];
  }
}

static inline bool row_equal_colmajor(const double* const* cols,
                                      int n_par,
                                      int row_a,
                                      int row_b) {
  for (int c = 0; c < n_par; ++c) {
    if (cols[c][row_a] != cols[c][row_b]) {
      return false;
    }
  }
  return true;
}

struct LogicalRulesScratch {
  std::vector<double> ll_unique;
  std::vector<double> logf_all;
  std::vector<double> logS_all;
  std::vector<int> all_mask;
  std::vector<int> isok_int_all;
  std::vector<double> GA_no_gl;
  std::vector<double> GB_no_gl;
  std::vector<unsigned char> ch_eq_vec;
  std::vector<unsigned char> pair_ok_A;
  std::vector<unsigned char> pair_ok_B;
  std::vector<double> gl_h_A;
  std::vector<double> gl_h_B;
  std::vector<double> gl_lo_A;
  std::vector<double> gl_lo_B;
  std::vector<double> pars_nA;
  std::vector<double> pars_A;
  std::vector<double> pars_nB;
  std::vector<double> pars_B;
  std::vector<int> isok_nA;
  std::vector<int> isok_A;
  std::vector<int> isok_nB;
  std::vector<int> isok_B;
  // Auxiliary GL sweeps (G_no at censor/truncation time points), slot-major
  // flat layout [slot * n_unique_trials + j]; see LrAuxSlot.
  std::vector<double> aux_t;
  std::vector<double> aux_hA;
  std::vector<double> aux_hB;
  std::vector<double> aux_GA;
  std::vector<double> aux_GB;
  std::vector<int> aux_mask;
  std::vector<int> isok_nA_base;
  std::vector<int> isok_A_base;
  std::vector<int> isok_nB_base;
  std::vector<int> isok_B_base;
  std::vector<int> all1;
  std::vector<double> tmpA1;
  std::vector<double> tmpA2;
  std::vector<double> tmpB1;
  std::vector<double> tmpB2;
  std::vector<double> parA;
  std::vector<double> parB;
  std::vector<double> parnA;
  std::vector<double> parnB;
  std::vector<double> pars_2buf;
};

// Auxiliary time-point slots for the logical-rules GL pre-pass: G_no is
// batch-evaluated at these per-trial times so the trial loop's truncation
// normaliser and censored branches need no per-trial GSL integration on the
// raw-kernel path.
enum LrAuxSlot {
  LR_AUX_LT = 0,   // truncation lower bound (Z), OR/AND only
  LR_AUX_UT,       // truncation upper bound (Z + upper-censor subtraction)
  LR_AUX_CLO,      // lower-censor window start max(0, LT)
  LR_AUX_CHI,      // lower-censor window end LC
  LR_AUX_UC,       // upper-censor bound, OR/AND only
  LR_AUX_N
};

static inline void resize_assign_double(std::vector<double>& x, size_t n, double value) {
  x.assign(n, value);
}

static inline void resize_assign_int(std::vector<int>& x, size_t n, int value) {
  x.assign(n, value);
}

static double local_race_helper(double t,
                                const double* par_target,
                                const double* par_nontarget,
                                int n_par,
                                ContextForRaceModels* ctx,
                                RacePdf1Fun pdf1,
                                RaceCdf1Fun cdf1,
                                const GslIntegrationControls& gsl_ctl,
                                gsl_integration_workspace* w,
                                double* pars_2buf,
                                int* isok_2buf) {
  if (!(t > 0.0) || !R_FINITE(t)) return 0.0;
  if (ctx == nullptr || pdf1 == nullptr || cdf1 == nullptr) return NA_REAL;

  const int t0_idx = ctx->t0_index;
  double t0_tgt = (t0_idx >= 0 && t0_idx < n_par) ? par_target[t0_idx] : 0.0;
  double t0_nt  = (t0_idx >= 0 && t0_idx < n_par) ? par_nontarget[t0_idx] : 0.0;
  if (t0_tgt < 0.0) t0_tgt = 0.0;
  if (t0_nt < 0.0) t0_nt = 0.0;

  if (t <= t0_tgt) return 0.0;

  double term1 = 0.0;
  const double split = std::min(t, t0_nt);
  if (split > t0_tgt) {
    double c_hi = cdf1(split, par_target, ctx);
    double c_lo = cdf1(t0_tgt, par_target, ctx);
    c_hi = clamp_cdf01_race(c_hi);
    c_lo = clamp_cdf01_race(c_lo);
    if (emc2_isnan(c_hi) || emc2_isnan(c_lo)) return NA_REAL;
    term1 = c_hi - c_lo;
    if (term1 < 0.0) term1 = 0.0;
    if (term1 > 1.0) term1 = 1.0;
  }

  const double low_limit = std::max(t0_tgt, t0_nt);
  if (t <= low_limit) return term1;

  const GLRule& gl = gl_get_rule(31);
  const double h = 0.5 * (t - low_limit);
  double sum_g = 0.0;
  for (size_t i = 0; i < gl.x.size(); ++i) {
    const double s = low_limit + h * (1.0 + gl.x[i]);
    const double f = pdf1(s, par_target, ctx);
    if (!R_FINITE(f) || f <= 0.0) continue;
    double F = clamp_cdf01_race(cdf1(s, par_nontarget, ctx));
    if (emc2_isnan(F)) return NA_REAL;
    const double S = std::max(0.0, 1.0 - F);
    sum_g += gl.w[i] * h * f * S;
  }

  double g = term1 + sum_g;
  if (g < 0.0) g = 0.0;
  if (g > 1.0) g = 1.0;
  return g;
}

// P(target strictly wins its two-accumulator channel by time t), i.e.
// int_0^t f_target(u) S_nontarget(u) du.  Unlike local_race_helper this
// accepts t = +Inf, integrating the tail with Gauss-Legendre quadrature.
// Returns a probability in [0, 1], or NA_REAL on integration failure.
static double lr_channel_yes_prob(double t,
                                  const double* par_target,
                                  const double* par_nontarget,
                                  int n_par,
                                  ContextForRaceModels* ctx,
                                  RacePdf1Fun pdf1,
                                  RaceCdf1Fun cdf1,
                                  const GslIntegrationControls& gsl_ctl,
                                  gsl_integration_workspace* w,
                                  double* pars_2buf,
                                  int* isok_2buf) {
  if (!(t > 0.0)) return 0.0;
  if (R_FINITE(t))
    return local_race_helper(t, par_target, par_nontarget, n_par, ctx,
                             pdf1, cdf1, gsl_ctl, w, pars_2buf, isok_2buf);
  if (ctx == nullptr || pdf1 == nullptr || cdf1 == nullptr)
    return NA_REAL;

  const int t0_idx = ctx->t0_index;
  double t0_tgt = (t0_idx >= 0 && t0_idx < n_par) ? par_target[t0_idx] : 0.0;
  double t0_nt  = (t0_idx >= 0 && t0_idx < n_par) ? par_nontarget[t0_idx] : 0.0;
  if (t0_tgt < 0.0) t0_tgt = 0.0;
  if (t0_nt < 0.0) t0_nt = 0.0;

  // Early window [t0_tgt, t0_nt] where only the target can finish (nontarget
  // survivor is exactly one): a plain CDF gap of the target.
  double term1 = 0.0;
  if (t0_nt > t0_tgt) {
    double c_hi = clamp_cdf01_race(cdf1(t0_nt, par_target, ctx));
    double c_lo = clamp_cdf01_race(cdf1(t0_tgt, par_target, ctx));
    if (emc2_isnan(c_hi) || emc2_isnan(c_lo)) return NA_REAL;
    term1 = c_hi - c_lo;
    if (term1 < 0.0) term1 = 0.0;
    if (term1 > 1.0) term1 = 1.0;
  }

  const double low_limit = std::max(t0_tgt, t0_nt);
  const GLRule& gl = gl_get_rule(31);
  double tail = 0.0;
  for (size_t i = 0; i < gl.x.size(); ++i) {
    const double u = 0.5 * (1.0 + gl.x[i]);
    const double om_u = 1.0 - u;
    if (om_u <= 1e-12) continue;
    const double s = low_limit + u / om_u;
    const double jacobian = 0.5 * gl.w[i] / (om_u * om_u);
    const double f = pdf1(s, par_target, ctx);
    if (!R_FINITE(f) || f <= 0.0) continue;
    double F = clamp_cdf01_race(cdf1(s, par_nontarget, ctx));
    if (emc2_isnan(F)) return NA_REAL;
    const double S = std::max(0.0, 1.0 - F);
    tail += jacobian * f * S;
  }

  double g = term1 + tail;
  if (g < 0.0) g = 0.0;
  if (g > 1.0) g = 1.0;
  return g;
}

static inline bool eval_pdf_cdf_race_scalar(
    bool use_raw_local,
    int idx,
    double t,
    const double* p_row,
    double min_ll,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const std::vector<double>& logf_all,
    const std::vector<double>& logS_all,
    double& f_out,
    double& F_out) {
  if (use_raw_local) {
    const double lf = logf_all[static_cast<size_t>(idx)];
    const double ls = logS_all[static_cast<size_t>(idx)];
    if (!R_FINITE(lf) || !R_FINITE(ls)) return false;
    const double f = (lf <= min_ll) ? 0.0 : std::exp(lf);
    double S = (ls <= min_ll) ? 0.0 : std::exp(ls);
    if (!(S >= 0.0) || !R_FINITE(S)) return false;
    if (S > 1.0) S = 1.0;
    const double F = clamp_cdf01_race(1.0 - S);
    if (!R_FINITE(f) || !R_FINITE(F) || emc2_isnan(F)) return false;
    f_out = f;
    F_out = F;
    return true;
  }
  const double f = pdf1(t, p_row, model_ctx);
  if (!R_FINITE(f) || f < 0.0) return false;
  double F = cdf1(t, p_row, model_ctx);
  F = clamp_cdf01_race(F);
  if (!R_FINITE(F) || emc2_isnan(F)) return false;
  f_out = f;
  F_out = F;
  return true;
}

static inline double safe_surv_at_race_scalar(
    const double* p_row,
    double t_eval,
    double min_surv,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx) {
  double F = cdf1(t_eval, p_row, model_ctx);
  F = clamp_cdf01_race(F);
  if (!R_FINITE(F) || emc2_isnan(F)) return 0.0;
  return std::max(min_surv, 1.0 - F);
}

// ---------------------------------------------------------------------------
// Logical-rules channel-resolution probabilities at a fixed time t.
//
// A "channel" is the 2-accumulator sub-race between a target ("yes") detector
// and a nontarget ("no") detector. At time t the channel is in exactly one of
// three states, with G_yes + G_no + S_dec = 1:
//   G_yes(t) = P(target won the sub-race at some u <= t)     channel said yes
//   G_no(t)  = P(nontarget won at some u <= t)               channel said no
//   S_dec(t) = P(neither finished by t) = S_tgt(t)*S_ntgt(t) undecided
// ---------------------------------------------------------------------------
struct LrChannelState {
  double G_yes = 0.0;
  double G_no = 0.0;
  double S_dec = 1.0;   // t = 0 state: undecided with certainty
};

// P(no overt response by t) for the OR/AND/XOR/ID rules. An overt response is
// NOT the first accumulator event: OR responds at the first channel-yes or the
// second channel-no; AND at the second yes or the first no; XOR/ID only once
// BOTH channels have resolved. Derived from channel independence, and
// consistent with the finite-rt densities used in the main trial loop
// (e.g. d/dt[1 - Q_A*Q_B] = gA_yes*Q_B + gB_yes*Q_A with Q = G_no + S_dec).
static inline double lr_no_response_prob(int rule_code,
                                         const LrChannelState& A,
                                         const LrChannelState& B) {
  if (rule_code == 3 || rule_code == 4) {  // XOR/ID: waiting iff either channel undecided
    return 1.0 - (1.0 - A.S_dec) * (1.0 - B.S_dec);
  }
  if (rule_code == 1) {                    // OR: no yes yet, minus both-no (responded "no")
    const double Q_A = A.G_no + A.S_dec;
    const double Q_B = B.G_no + B.S_dec;
    return Q_A * Q_B - A.G_no * B.G_no;
  }
  // AND: no no yet, minus both-yes (responded "yes")
  const double Q_A = A.G_yes + A.S_dec;
  const double Q_B = B.G_yes + B.S_dec;
  return Q_A * Q_B - A.G_yes * B.G_yes;
}

// CDF of (response identity, RT <= t) for OR/AND/XOR/ID. resp_code follows
// LogicalRulesSharedState (0=unknown, 1=yes, 2=no, 3=NN, 4=AN, 5=NB, 6=AB);
// resp_code 0 marginalizes over identities. Returns NA_REAL for a response
// code that is impossible under the rule.
static inline double lr_response_cdf(int rule_code, int resp_code,
                                     const LrChannelState& A,
                                     const LrChannelState& B) {
  if (resp_code == 0) return 1.0 - lr_no_response_prob(rule_code, A, B);
  switch (rule_code) {
  case 1:  // OR: yes at first channel-yes; no once both channels said no
    if (resp_code == 1) return 1.0 - (A.G_no + A.S_dec) * (B.G_no + B.S_dec);
    if (resp_code == 2) return A.G_no * B.G_no;
    break;
  case 2:  // AND: yes once both channels said yes; no at first channel-no
    if (resp_code == 1) return A.G_yes * B.G_yes;
    if (resp_code == 2) return 1.0 - (A.G_yes + A.S_dec) * (B.G_yes + B.S_dec);
    break;
  case 3:  // XOR: yes = exactly one channel yes; both channels resolved
    if (resp_code == 1) return A.G_yes * B.G_no + A.G_no * B.G_yes;
    if (resp_code == 2) return A.G_yes * B.G_yes + A.G_no * B.G_no;
    break;
  case 4:  // ID: identity = the resolved (yes/no) pair
    if (resp_code == 6) return A.G_yes * B.G_yes;
    if (resp_code == 4) return A.G_yes * B.G_no;
    if (resp_code == 5) return A.G_no  * B.G_yes;
    if (resp_code == 3) return A.G_no  * B.G_no;
    break;
  }
  return NA_REAL;
}

// P(no overt response by t) for the legacy analytic-detection helper. The
// four-horse GNG route is handled by lr4_N_at below; it has no separate nogo
// accumulator. Used for upper-censor masses and the truncation normaliser.
// t must be finite; t <= 0 returns 1.
static double lr_detection_no_response_prob(
    int rule_code,
    int cond,
    double t,
    const double* parA,
    const double* parB,
    const double* parN,
    int n_par,
    double min_surv,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& pars_rowmajor_buf,
    bool& ok) {
  ok = true;
  if (!(t > 0.0)) return 1.0;
  if (!R_FINITE(t)) { ok = false; return 1.0; }
  if (cond == 0 && rule_code == 5) return 1.0;  // no detector: no response possible
  const double S_A = (cond == 1 || cond == 3)
    ? safe_surv_at_race_scalar(parA, t, min_surv, cdf1, model_ctx) : 1.0;
  const double S_B = (cond == 2 || cond == 3)
    ? safe_surv_at_race_scalar(parB, t, min_surv, cdf1, model_ctx) : 1.0;
  if (rule_code == 5) return S_A * S_B;

  // rule 6, cond 0: only the nogo accumulator races, so no overt response can
  // ever occur — S_N(t) + P(nogo won by t) = 1 identically.
  if (cond == 0) return 1.0;

  // rule 6: joint survivor of the active set, plus nogo-win mass on [0, t]
  const double S_N = safe_surv_at_race_scalar(parN, t, min_surv, cdf1, model_ctx);
  const double p_none = S_A * S_B * S_N;

  const int t0_idx = model_ctx->t0_index;
  const double t0_N = (t0_idx >= 0 && t0_idx < n_par) ? parN[t0_idx] : 0.0;
  const double lo = std::max(0.0, t0_N);
  const double h = (t > lo) ? 0.5 * (t - lo) : 0.0;
  double p_nogo_win = 0.0;
  if (h > 0.0) {
    const GLRule& gl = gl_get_rule(31);
    for (size_t i = 0; i < gl.x.size(); ++i) {
      const double s = lo + h * (1.0 + gl.x[i]);
      const double fN = pdf1(s, parN, model_ctx);
      if (!R_FINITE(fN) || fN <= 0.0) continue;
      double SA = 1.0, SB = 1.0;
      if (cond == 1 || cond == 3) SA = safe_surv_at_race_scalar(parA, s, min_surv, cdf1, model_ctx);
      if (cond == 2 || cond == 3) SB = safe_surv_at_race_scalar(parB, s, min_surv, cdf1, model_ctx);
      p_nogo_win += gl.w[i] * h * fN * SA * SB;
    }
  }

  return std::min(1.0, p_none + p_nogo_win);
}

static double logicalrules_detection_trial_ll(
    int rule_code,
    int cond,
    int resp,
    double t,
    double LTj,
    double UTj,
    double LCj,
    double UCj,
    int idxA,
    int idxB,
    int idxN,
    int n_par,
    double min_ll,
    double min_surv,
    double log_Z,
    bool use_raw_local,
    const std::vector<double>& logf_all,
    const std::vector<double>& logS_all,
    const double* parA,
    const double* parB,
    const double* parN,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& pars_rowmajor_buf) {
  const bool finite_rt = (R_FINITE(t) && t > 0.0);
  double p_j = 0.0;

  if (finite_rt) {
    if (resp != 1) return min_ll;
    // A detection stimulus activates only the detector(s) named by `cond`.
    // The fixed-role design still carries both A and B rows so that every
    // trial has the same layout, but an inactive row must not be evaluated:
    // besides being unnecessary, an invalid inactive parameter could
    // otherwise floor an otherwise valid A-only or B-only trial.
    const bool activeA = (cond == 1 || cond == 3);
    const bool activeB = (cond == 2 || cond == 3);
    if (!activeA && !activeB) return min_ll;
    double fA = 0.0, FA = 0.0, fB = 0.0, FB = 0.0;
    if (activeA && !eval_pdf_cdf_race_scalar(
          use_raw_local, idxA, t, parA, min_ll, pdf1, cdf1, model_ctx,
          logf_all, logS_all, fA, FA)) return min_ll;
    if (activeB && !eval_pdf_cdf_race_scalar(
          use_raw_local, idxB, t, parB, min_ll, pdf1, cdf1, model_ctx,
          logf_all, logS_all, fB, FB)) return min_ll;
    const double S_A = activeA ? std::max(min_surv, 1.0 - FA) : 1.0;
    const double S_B = activeB ? std::max(min_surv, 1.0 - FB) : 1.0;
    if (rule_code == 5) {
      if (cond == 1) p_j = fA;
      else if (cond == 2) p_j = fB;
      else p_j = fA * S_B + fB * S_A;
    } else {
      double fN = 0.0, FN = 0.0;
      if (!eval_pdf_cdf_race_scalar(
            use_raw_local, idxN, t, parN, min_ll, pdf1, cdf1, model_ctx,
            logf_all, logS_all, fN, FN)) return min_ll;
      const double S_N = std::max(min_surv, 1.0 - FN);
      if (cond == 1) p_j = fA * S_N;
      else if (cond == 2) p_j = fB * S_N;
      else p_j = (fA * S_B + fB * S_A) * S_N;
    }
  } else {
    // Upper censor (+Inf) is an omission: recorded resp may be missing (0) or
    // an explicit "no" (2). Lower censor (-Inf) is an overt detection response
    // that beat LC: resp may be missing (0) or "yes" (1) — the mass is the
    // same either way, as these rules have a single overt response type.
    const bool upper_censored = (t > 0.0);
    if (upper_censored ? (resp != 0 && resp != 2)
                       : (resp != 0 && resp != 1)) return min_ll;
    if (rule_code == 5) {
      // Upper censor (rt = +Inf): neither detector fired by UC. With a finite
      // upper truncation bound (and a proper upper tail) the response is known
      // to land in (UC, UT], so N(UT) is subtracted; the trial mass is then
      // normalised by Z = N(LT) - N(UT) in the caller.
      // Lower censor (rt = -Inf): at least one detector fired before LC.
      if (t > 0.0) {
        if (!R_FINITE(UCj))
          Rcpp::stop("LogicalRules OR_DETECTION_ANALYTIC requires finite UC for upper-censored trials.");
        bool okd = true;
        p_j = lr_detection_no_response_prob(rule_code, cond, UCj, parA, parB, parN,
                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                            gsl_ctl, w, pars_rowmajor_buf, okd);
        if (okd && R_FINITE(UTj) && !model_ctx->defective_upper_tail) {
          bool ok_ut = true;
          const double N_UT = lr_detection_no_response_prob(rule_code, cond, UTj, parA, parB, parN,
                                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                                            gsl_ctl, w, pars_rowmajor_buf, ok_ut);
          okd = okd && ok_ut;
          p_j = std::max(0.0, p_j - N_UT);
        }
        if (!okd) return min_ll;
      } else {
        const double lo = std::max(0.0, LTj);
        const double hi = std::max(lo, LCj);
        if (cond == 0) {
          p_j = 0.0;
        } else if (cond == 1) {
          const double S_A_lo = safe_surv_at_race_scalar(parA, lo, min_surv, cdf1, model_ctx);
          const double S_A_hi = safe_surv_at_race_scalar(parA, hi, min_surv, cdf1, model_ctx);
          p_j = std::max(0.0, S_A_lo - S_A_hi);
        } else if (cond == 2) {
          const double S_B_lo = safe_surv_at_race_scalar(parB, lo, min_surv, cdf1, model_ctx);
          const double S_B_hi = safe_surv_at_race_scalar(parB, hi, min_surv, cdf1, model_ctx);
          p_j = std::max(0.0, S_B_lo - S_B_hi);
        } else {
          const double S_A_lo = safe_surv_at_race_scalar(parA, lo, min_surv, cdf1, model_ctx);
          const double S_B_lo = safe_surv_at_race_scalar(parB, lo, min_surv, cdf1, model_ctx);
          const double S_A_hi = safe_surv_at_race_scalar(parA, hi, min_surv, cdf1, model_ctx);
          const double S_B_hi = safe_surv_at_race_scalar(parB, hi, min_surv, cdf1, model_ctx);
          p_j = std::max(0.0, S_A_lo * S_B_lo - S_A_hi * S_B_hi);
        }
      }
    } else {
      if (t > 0.0) {
        // Upper censor: no overt response by UC = nogo won by UC or nobody
        // fired by UC. The nogo-win mass now accumulates from 0 (not LT):
        // truncation conditioning is handled by the caller's normaliser Z,
        // and N(UT) is subtracted when the window's upper bound is finite.
        if (!R_FINITE(UCj))
          Rcpp::stop("LogicalRules OR_DETECTION_GNG requires finite UC for upper-censored trials.");
        bool okd = true;
        p_j = lr_detection_no_response_prob(rule_code, cond, UCj, parA, parB, parN,
                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                            gsl_ctl, w, pars_rowmajor_buf, okd);
        if (okd && R_FINITE(UTj) && !model_ctx->defective_upper_tail) {
          bool ok_ut = true;
          const double N_UT = lr_detection_no_response_prob(rule_code, cond, UTj, parA, parB, parN,
                                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                                            gsl_ctl, w, pars_rowmajor_buf, ok_ut);
          okd = okd && ok_ut;
          p_j = std::max(0.0, p_j - N_UT);
        }
        if (!okd) return min_ll;
      } else {
        // Lower censor: an overt (go) response occurred in [LT, LC]. Only the
        // go detectors active under this stimulus can produce it, and the
        // detector must win the race against the nogo accumulator (and any
        // other go detector) inside the window — a nogo finish is a withheld
        // response, never an observable one (mirrors the standard-race
        // go/no-go convention for rt == -Inf). Previously this used the joint
        // first-event mass, which wrongly counted nogo wins as responses
        // (e.g. cond == 0, where no overt response is possible at all).
        const double lo = std::max(0.0, LTj);
        const double hi = std::max(lo, LCj);
        if (cond == 0 || !(hi > lo)) {
          p_j = 0.0;
        } else {
          const GLRule& gl = gl_get_rule(31);
          const double h = 0.5 * (hi - lo);
          double p_go_win = 0.0;
          for (size_t i = 0; i < gl.x.size(); ++i) {
            const double s = lo + h * (1.0 + gl.x[i]);
            const double SN = safe_surv_at_race_scalar(parN, s, min_surv, cdf1, model_ctx);
            if (cond == 1 || cond == 3) {
              const double fA = pdf1(s, parA, model_ctx);
              if (R_FINITE(fA) && fA > 0.0) {
                const double SB = (cond == 3) ? safe_surv_at_race_scalar(parB, s, min_surv, cdf1, model_ctx) : 1.0;
                p_go_win += gl.w[i] * h * fA * SB * SN;
              }
            }
            if (cond == 2 || cond == 3) {
              const double fB = pdf1(s, parB, model_ctx);
              if (R_FINITE(fB) && fB > 0.0) {
                const double SA = (cond == 3) ? safe_surv_at_race_scalar(parA, s, min_surv, cdf1, model_ctx) : 1.0;
                p_go_win += gl.w[i] * h * fB * SA * SN;
              }
            }
          }
          p_j = p_go_win;
        }
      }
    }
  }

  if (!(p_j > 0.0) || !R_FINITE(p_j)) return min_ll;
  const double ll = std::log(p_j) - log_Z;
  return (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
}

// ===========================================================================
// LogicalRules correlated-target capacity routes.
//
// On a redundant-target (AB) trial one latent standard-normal factor z adds
// a shared shift to both target drift means: V_i | z ~
// N(v_i + kappa + tau*z, sv_i^2) for the
// A and B target rows; nontarget/nogo racers never load on the factor.  The
// logical rule is assembled conditional on z from normalized conditional
// channel quantities and the complete trial mass is integrated over z with a
// scan/recentred Gauss-Hermite pair (gh_quad.h).  Positive-drift models
// weight each node by q_A(z)*q_B(z) and divide once by the closed-form
// bivariate orthant probability q_AB (bawl_corr_pair_positive_normalizer);
// independent racers keep their ordinary posdrift-normalized kernels, whose
// constant q factors cancel between numerator and denominator.
// LogicalRulesLBA is plain LBA, i.e. the exact k = 0 member of the BAwL
// geometry, so all conditional target endpoints come from BAwLPreparedRow.
// ===========================================================================

// GL order for the conditional channel integrals; matches the batched GL
// pre-pass order used by the ordinary route.
static constexpr int LRCAP_N_GL = 31;

// Node-invariant conditional drift view of one loaded target row:
// V | z ~ N(v0 + slope*z, sv^2) with v0 = v + kappa and slope = tau.
struct LrCapTargetView {
  double v0 = 0.0, slope = 0.0, sv = 1.0;
  double t0 = 0.0, A = 0.0, b = 0.0;
  bool valid = false;
};

static inline LrCapTargetView lrcap_target_view(const double* const* cols,
                                                int row, double kappa,
                                                double tau) {
  LrCapTargetView tv;
  const double v = cols[emc2col::lba::v][row];
  const double sv = cols[emc2col::lba::sv][row];
  const double B = cols[emc2col::lba::B][row];
  const double A = cols[emc2col::lba::A][row];
  tv.t0 = cols[emc2col::lba::t0][row];
  tv.A = A;
  tv.b = B + A;
  tv.v0 = v + kappa;
  tv.slope = tau;
  tv.sv = sv;
  tv.valid = R_FINITE(v) && sv > 0.0 && A >= 0.0 && tv.b > 0.0 && tv.b >= A &&
    R_FINITE(tv.v0) && R_FINITE(tv.slope) &&
    R_FINITE(tv.t0) && tv.t0 >= 0.0;
  return tv;
}

static inline BAwLPreparedRow lrcap_prepare_at(const LrCapTargetView& tv,
                                               double t) {
  return bawl_prepare_row(bawl_time_geometry(t, tv.t0, tv.A, tv.b, 0.0),
                          tv.v0, tv.slope, tv.sv);
}

// log q(z) = log Phi(v_q / sv): the target's positivity constant, shared by
// every evaluation time of the row (it depends only on the drift model).
static inline double lrcap_log_q(const LrCapTargetView& tv, double z,
                                 bool posdrift) {
  if (!posdrift) return 0.0;
  return pnorm_log_direct((tv.v0 + tv.slope * z) / tv.sv, true);
}

// Normalized conditional survivor of a prepared target row at latent z.  The
// caller supplies log_q once per (row, z); a numerically hopeless node (or an
// invalid prepared row) reports ok = false and the node contributes zero
// integrand mass rather than failing the trial.
static inline double lrcap_cond_survivor(const BAwLPreparedRow& prep, double z,
                                         double log_q, bool posdrift,
                                         bool& ok) {
  if (prep.status == BAwLTimeStatus::invalid) { ok = false; return 0.0; }
  const double vq = prep.v0 + prep.slope * z;
  const double ls = bawl_prepared_log_survivor(prep, vq, posdrift);
  if (ISNAN(ls)) { ok = false; return 0.0; }
  if (ls == R_NegInf) return 0.0;
  double S = std::exp(ls - log_q);
  if (!R_FINITE(S)) { ok = false; return 0.0; }
  if (S < 0.0) S = 0.0;
  if (S > 1.0) S = 1.0;
  return S;
}

// Normalized conditional density of a prepared target row at latent z.
static inline double lrcap_cond_pdf(const BAwLPreparedRow& prep, double z,
                                    double log_q, bool& ok) {
  if (prep.status == BAwLTimeStatus::invalid) { ok = false; return 0.0; }
  const double vq = prep.v0 + prep.slope * z;
  const double lf = bawl_prepared_log_pdf(prep, vq);
  if (ISNAN(lf)) { ok = false; return 0.0; }
  if (lf == R_NegInf) return 0.0;
  const double f = std::exp(lf - log_q);
  if (!R_FINITE(f) || f < 0.0) { ok = false; return 0.0; }
  return f;
}

// One channel (loaded target vs independent nontarget) prepared at one
// evaluation time.  Everything except the conditional target mean is
// node-invariant: the independent racer's density/survivor, the prepared
// target geometry at the endpoint, and — when the channel-no CDF is needed —
// the GL abscissa weights wt_i*h*f_ind(s_i) with per-abscissa prepared target
// geometry (plan phase 4's two mandatory caches).
struct LrCapChannelPoint {
  bool built = false;
  bool at_zero = false;              // tt <= 0: channel still undecided
  BAwLPreparedRow prep_t;
  double f_ind_t = 0.0;              // independent density at tt
  double S_ind_t = 1.0;              // independent survivor at tt
  int n_gl = 0;
  std::array<double, LRCAP_N_GL> gl_w{};
  std::array<BAwLPreparedRow, LRCAP_N_GL> gl_prep{};
};

static void lrcap_build_channel_point(LrCapChannelPoint& cp,
                                      const LrCapTargetView& tv,
                                      const double* par_ind, double tt,
                                      bool need_G, bool need_f_ind,
                                      RacePdf1Fun pdf1, RaceCdf1Fun cdf1,
                                      ContextForRaceModels* ctx,
                                      double min_surv) {
  cp = LrCapChannelPoint();
  cp.built = true;
  if (!(tt > 0.0)) { cp.at_zero = true; return; }
  cp.prep_t = lrcap_prepare_at(tv, tt);
  cp.S_ind_t = safe_surv_at_race_scalar(par_ind, tt, min_surv, cdf1, ctx);
  if (need_f_ind) {
    const double f = pdf1(tt, par_ind, ctx);
    cp.f_ind_t = (R_FINITE(f) && f > 0.0) ? f : 0.0;
  }
  if (need_G) {
    const int t0_idx = ctx->t0_index;
    const double t0n = (t0_idx >= 0) ? par_ind[t0_idx] : 0.0;
    const double lo = std::max(0.0, t0n);
    const double h = (tt > lo) ? 0.5 * (tt - lo) : 0.0;
    if (h > 0.0) {
      const GLRule& gl = gl_get_rule(LRCAP_N_GL);
      cp.n_gl = LRCAP_N_GL;
      for (int i = 0; i < LRCAP_N_GL; ++i) {
        const double s = lo + h * (1.0 + gl.x[static_cast<size_t>(i)]);
        const double f = pdf1(s, par_ind, ctx);
        cp.gl_w[static_cast<size_t>(i)] =
          (R_FINITE(f) && f > 0.0) ? gl.w[static_cast<size_t>(i)] * h * f : 0.0;
        cp.gl_prep[static_cast<size_t>(i)] = lrcap_prepare_at(tv, s);
      }
    }
  }
}

// Complete conditional channel quantities at the point's time given latent z.
struct LrCapChanVals {
  double f_t = 0.0;   // conditional target density
  double S_t = 0.0;   // conditional target survivor
  double f_n = 0.0;   // independent density
  double S_n = 1.0;   // independent survivor
  double G_no = 0.0;  // P(nontarget won the channel by tt | z)
  bool ok = true;
};

static LrCapChanVals lrcap_chan_vals(const LrCapChannelPoint& cp, double z,
                                     double log_q, bool posdrift,
                                     bool need_pdf, bool count) {
  LrCapChanVals out;
  if (cp.at_zero) { out.S_t = 1.0; return out; }
  out.S_t = lrcap_cond_survivor(cp.prep_t, z, log_q, posdrift, out.ok);
  if (!out.ok) return out;
  out.f_n = cp.f_ind_t;
  out.S_n = cp.S_ind_t;
  if (need_pdf) {
    out.f_t = lrcap_cond_pdf(cp.prep_t, z, log_q, out.ok);
    if (!out.ok) return out;
  }
  double G_no = 0.0;
  for (int i = 0; i < cp.n_gl; ++i) {
    const double w = cp.gl_w[static_cast<size_t>(i)];
    if (w == 0.0) continue;
    G_no += w * lrcap_cond_survivor(cp.gl_prep[static_cast<size_t>(i)], z,
                                    log_q, posdrift, out.ok);
    if (!out.ok) return out;
  }
  if (count && cp.n_gl > 0)
    lr_capacity_counters().channel_gl_evaluations += cp.n_gl;
  out.G_no = std::min(1.0, std::max(0.0, G_no));
  return out;
}

static inline LrChannelState lrcap_state_from_vals(const LrCapChanVals& v) {
  LrChannelState st;
  st.S_dec = v.S_t * v.S_n;
  st.G_no = v.G_no;
  st.G_yes = std::min(1.0, std::max(0.0, 1.0 - st.S_dec - st.G_no));
  return st;
}

// One detection evaluation point: prepared A/B target rows plus (for GNG) the
// nogo survivor and the z-dependent GL integral caches.  `gl_w` carries
// wt_i*h*f_N(s_i) for the nogo-win mass, or wt_i*h*S_N(s_i) for the go-win
// window mass, depending on which evaluator consumes the point.
struct LrCapDetPoint {
  bool built = false;
  bool at_zero = false;
  BAwLPreparedRow prepA, prepB;
  double S_N = 1.0;
  int n_gl = 0;
  std::array<double, LRCAP_N_GL> gl_w{};
  std::array<BAwLPreparedRow, LRCAP_N_GL> gl_prepA{};
  std::array<BAwLPreparedRow, LRCAP_N_GL> gl_prepB{};
};

// N(tt) point for OR_DETECTION_ANALYTIC / OR_DETECTION_GNG: survivors at tt,
// and for GNG the nogo-win mass cache over [max(0, t0_N), tt].
static void lrcap_build_det_N_point(LrCapDetPoint& dp,
                                    const LrCapTargetView& tvA,
                                    const LrCapTargetView& tvB,
                                    const double* parN, bool has_nogo,
                                    double tt, bool need_nogo_win,
                                    RacePdf1Fun pdf1, RaceCdf1Fun cdf1,
                                    ContextForRaceModels* ctx,
                                    double min_surv) {
  dp = LrCapDetPoint();
  dp.built = true;
  if (!(tt > 0.0)) { dp.at_zero = true; return; }
  dp.prepA = lrcap_prepare_at(tvA, tt);
  dp.prepB = lrcap_prepare_at(tvB, tt);
  if (has_nogo) {
    dp.S_N = safe_surv_at_race_scalar(parN, tt, min_surv, cdf1, ctx);
    if (need_nogo_win) {
      const int t0_idx = ctx->t0_index;
      const double t0n = (t0_idx >= 0) ? parN[t0_idx] : 0.0;
      const double lo = std::max(0.0, t0n);
      const double h = (tt > lo) ? 0.5 * (tt - lo) : 0.0;
      if (h > 0.0) {
        const GLRule& gl = gl_get_rule(LRCAP_N_GL);
        dp.n_gl = LRCAP_N_GL;
        for (int i = 0; i < LRCAP_N_GL; ++i) {
          const double s = lo + h * (1.0 + gl.x[static_cast<size_t>(i)]);
          const double f = pdf1(s, parN, ctx);
          dp.gl_w[static_cast<size_t>(i)] =
            (R_FINITE(f) && f > 0.0) ? gl.w[static_cast<size_t>(i)] * h * f : 0.0;
          dp.gl_prepA[static_cast<size_t>(i)] = lrcap_prepare_at(tvA, s);
          dp.gl_prepB[static_cast<size_t>(i)] = lrcap_prepare_at(tvB, s);
        }
      }
    }
  }
}

// Go-win window point for the GNG lower-censor mass over [lo, hi]:
// integral of [f_A(s|z) S_B(s|z) + f_B(s|z) S_A(s|z)] * S_N(s) ds.
static void lrcap_build_det_gowin_point(LrCapDetPoint& dp,
                                        const LrCapTargetView& tvA,
                                        const LrCapTargetView& tvB,
                                        const double* parN, bool has_nogo,
                                        double lo, double hi,
                                        RaceCdf1Fun cdf1,
                                        ContextForRaceModels* ctx,
                                        double min_surv) {
  dp = LrCapDetPoint();
  dp.built = true;
  const double h = (hi > lo) ? 0.5 * (hi - lo) : 0.0;
  if (!(h > 0.0)) { dp.at_zero = true; return; }
  const GLRule& gl = gl_get_rule(LRCAP_N_GL);
  dp.n_gl = LRCAP_N_GL;
  for (int i = 0; i < LRCAP_N_GL; ++i) {
    const double s = lo + h * (1.0 + gl.x[static_cast<size_t>(i)]);
    const double SN = has_nogo
      ? safe_surv_at_race_scalar(parN, s, min_surv, cdf1, ctx) : 1.0;
    dp.gl_w[static_cast<size_t>(i)] = gl.w[static_cast<size_t>(i)] * h * SN;
    dp.gl_prepA[static_cast<size_t>(i)] = lrcap_prepare_at(tvA, s);
    dp.gl_prepB[static_cast<size_t>(i)] = lrcap_prepare_at(tvB, s);
  }
}

// Conditional detection no-response probability at a prepared N point.
static double lrcap_det_N_at(const LrCapDetPoint& dp, int rule_code, double z,
                             double lqA, double lqB, bool posdrift, bool& ok,
                             bool count) {
  if (dp.at_zero) return 1.0;
  const double SA = lrcap_cond_survivor(dp.prepA, z, lqA, posdrift, ok);
  if (!ok) return 1.0;
  const double SB = lrcap_cond_survivor(dp.prepB, z, lqB, posdrift, ok);
  if (!ok) return 1.0;
  if (rule_code == 5) return SA * SB;
  double p = SA * SB * dp.S_N;
  double nogo_win = 0.0;
  for (int i = 0; i < dp.n_gl; ++i) {
    const double w = dp.gl_w[static_cast<size_t>(i)];
    if (w == 0.0) continue;
    const double SAi = lrcap_cond_survivor(dp.gl_prepA[static_cast<size_t>(i)],
                                           z, lqA, posdrift, ok);
    if (!ok) return 1.0;
    const double SBi = lrcap_cond_survivor(dp.gl_prepB[static_cast<size_t>(i)],
                                           z, lqB, posdrift, ok);
    if (!ok) return 1.0;
    nogo_win += w * SAi * SBi;
  }
  if (count && dp.n_gl > 0)
    lr_capacity_counters().channel_gl_evaluations += 2LL * dp.n_gl;
  return std::min(1.0, p + nogo_win);
}

// Conditional go-win window mass at a prepared go-win point.
static double lrcap_det_gowin_at(const LrCapDetPoint& dp, double z,
                                 double lqA, double lqB, bool posdrift,
                                 bool& ok, bool count) {
  if (dp.at_zero) return 0.0;
  double mass = 0.0;
  for (int i = 0; i < dp.n_gl; ++i) {
    const double w = dp.gl_w[static_cast<size_t>(i)];
    if (w == 0.0) continue;
    const double fAi = lrcap_cond_pdf(dp.gl_prepA[static_cast<size_t>(i)], z,
                                      lqA, ok);
    if (!ok) return 0.0;
    const double SAi = lrcap_cond_survivor(dp.gl_prepA[static_cast<size_t>(i)],
                                           z, lqA, posdrift, ok);
    if (!ok) return 0.0;
    const double fBi = lrcap_cond_pdf(dp.gl_prepB[static_cast<size_t>(i)], z,
                                      lqB, ok);
    if (!ok) return 0.0;
    const double SBi = lrcap_cond_survivor(dp.gl_prepB[static_cast<size_t>(i)],
                                           z, lqB, posdrift, ok);
    if (!ok) return 0.0;
    mass += w * (fAi * SBi + fBi * SAi);
  }
  if (count && dp.n_gl > 0)
    lr_capacity_counters().channel_gl_evaluations += 2LL * dp.n_gl;
  return mass;
}

// Shared-factor integration: a fixed scan pass locates each integrand's mass
// (and is the fallback value), then a recentred pass refines it
// (gh_quad.h; scan-based centring per the capacity plan, phase 5).
template <typename LogMass>
static double lrcap_integrate_factor(LogMass&& log_mass, int n_scan,
                                     int n_fine, bool count,
                                     bool bounded_mass = false) {
  // Censoring masses and truncation-window masses are probabilities.  A
  // quadrature approximation with a positive log mass is therefore a
  // numerical failure, not a plausible large likelihood.  This is important
  // for the recentered pass: a narrow positive-drift integrand can be missed
  // by the scan and then badly over-weighted by the refined rule.
  static constexpr double kLogProbabilityTol = 1e-8;
  const double sqrt2 = std::sqrt(2.0);
  const GHRule& scan_rule = gh_rule(n_scan);
  std::array<double, 256> g{};
  std::array<double, 256> zv{};
  double scan_lse = R_NegInf;
  for (int q = 0; q < n_scan; ++q) {
    const double z = sqrt2 * scan_rule.x[static_cast<size_t>(q)];
    zv[static_cast<size_t>(q)] = z;
    const double w = gh_standard_normal_weight(scan_rule, q);
    const double v = (R_FINITE(w) && w > 0.0)
      ? std::log(w) + log_mass(z) : R_NegInf;
    g[static_cast<size_t>(q)] = v;
    scan_lse = log_sum_exp(scan_lse, v);
  }
  if (count) lr_capacity_counters().factor_node_evaluations += n_scan;
  if (!R_FINITE(scan_lse)) return R_NegInf;
  const AGHCenter c = agh_center_from_scan(g.data(), zv.data(), n_scan);
  const GHRule& fine_rule = gh_rule(n_fine);
  double fine_lse = R_NegInf;
  for (int q = 0; q < n_fine; ++q) {
    const double w = fine_rule.w[static_cast<size_t>(q)];
    if (!(R_FINITE(w) && w > 0.0)) continue;
    const double lgw = std::log(w) +
      fine_rule.x[static_cast<size_t>(q)] * fine_rule.x[static_cast<size_t>(q)];
    const double z = c.mu + c.sigma * sqrt2 * fine_rule.x[static_cast<size_t>(q)];
    fine_lse = log_sum_exp(fine_lse,
                           agh_log_weight(lgw, c.sigma, z) + log_mass(z));
  }
  if (count) lr_capacity_counters().factor_node_evaluations += n_fine;
  if (bounded_mass) {
    const bool scan_ok = R_FINITE(scan_lse) &&
      scan_lse <= kLogProbabilityTol;
    const bool fine_ok = R_FINITE(fine_lse) &&
      fine_lse <= kLogProbabilityTol;
    if (fine_ok) return fine_lse;
    // The fixed scan is a deliberately conservative fallback when the
    // recentered rule violates the probability bound.  If both rules fail,
    // report an invalid mass so the caller floors this particle/trial rather
    // than exposing an artificial positive likelihood to the sampler.
    return scan_ok ? scan_lse : R_NegInf;
  }
  return R_FINITE(fine_lse) ? fine_lse : scan_lse;
}

// Complete per-trial likelihood for one capacity (AB, kappa/tau active)
// trial.  Mirrors the ordinary logical-rules branches conditional on the
// factor: the rule is assembled inside the integral, censoring masses and the
// truncation normaliser are z-integrated, and the final contract is
// log num - log den with den the closed-form q_AB (untruncated posdrift) or
// the z-integrated truncation window.
static double lrcap_trial_ll(
    int rule_code, int cond, int resp,
    double t, double LTj, double UTj, double LCj, double UCj,
    double kappa, double tau,
    const double* const* pars_cols, int n_par,
    int idxA, int idxB,
    const double* parnA, const double* parnB, const double* parN,
    bool has_channels, bool has_nogo,
    double min_ll, double min_surv,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, ContextForRaceModels* ctx,
    int n_scan, int n_fine, bool count,
    int cell_idx, std::vector<double>* cell_log_den_cap,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& gng_scratch) {
  if (cond != 3) return min_ll;  // capacity trials are AB by construction
  const bool posdrift = ctx->use_posdrift;
  // GNG (rule 6) is the full four-horse OR task whose "no" outcome is a
  // withheld response (rt = +Inf); it shares OR's finite-RT go density but
  // routes its censor/withheld/truncation masses through the conditional
  // "no overt go by t" probability N_GNG(t | z) = Q_A(t|z) Q_B(t|z).
  const bool is_gng = (rule_code == 6);
  const LrCapTargetView tvA = lrcap_target_view(pars_cols, idxA, kappa, tau);
  const LrCapTargetView tvB = lrcap_target_view(pars_cols, idxB, kappa, tau);
  if (!tvA.valid || !tvB.valid) return min_ll;
  const bool finite_rt = R_FINITE(t) && t > 0.0;
  const bool upper_censored = !finite_rt && t == R_PosInf;
  const bool lower_censored = !finite_rt && t == R_NegInf;
  if (!finite_rt && !upper_censored && !lower_censored) return min_ll;

  // Closed-form joint positivity normalizer for the pair (untruncated
  // denominator); the quadrature identity int qA(z) qB(z) phi(z) dz is kept
  // only as a test cross-check.
  double log_qAB = 0.0;
  if (posdrift) {
    if (tau == 0.0) {
      log_qAB = pnorm_log_direct(tvA.v0 / tvA.sv, true) +
        pnorm_log_direct(tvB.v0 / tvB.sv, true);
    } else {
      const double sdA = std::sqrt(tvA.sv * tvA.sv + tvA.slope * tvA.slope);
      const double sdB = std::sqrt(tvB.sv * tvB.sv + tvB.slope * tvB.slope);
      const double rho = (tvA.slope * tvB.slope) / (sdA * sdB);
      const double qAB = bawl_corr_pair_positive_normalizer(
          tvA.v0, sdA, tvB.v0, sdB, rho, true);
      if (!R_FINITE(qAB) || !(qAB > 0.0)) return min_ll;
      log_qAB = std::log(std::fmin(qAB, 1.0));
    }
    if (!R_FINITE(log_qAB)) return min_ll;
  }

  const bool is_detection = (rule_code == 5);
  const bool need_G_for_N = (rule_code == 1 || rule_code == 2);
  const bool has_trunc = (LTj != 0.0 || R_FINITE(UTj));
  const bool subtract_UT = R_FINITE(UTj) && !ctx->defective_upper_tail;

  // The capacity denominator is independent of the observed RT. Cache it
  // by cell_idx, leaving the RT-specific event numerator on the per-trial path.
  double cached_log_den = R_NegInf;
  bool z_cache_hit = false;
  if (has_trunc && cell_log_den_cap != nullptr && cell_idx >= 0 &&
      cell_idx < static_cast<int>(cell_log_den_cap->size())) {
    if (!ISNAN((*cell_log_den_cap)[cell_idx])) {
      z_cache_hit = true;
      cached_log_den = (*cell_log_den_cap)[cell_idx];
      if (cached_log_den == R_NegInf) return min_ll;
    }
  }

  // Response gates mirror the ordinary branches; an incompatible recorded
  // response is zero mass for every factor value.
  if (is_detection) {
    if (finite_rt && resp != 1) return min_ll;
    if (upper_censored && resp != 0 && resp != 2) return min_ll;
    if (lower_censored && resp != 0 && resp != 1) return min_ll;
    if (upper_censored && !R_FINITE(UCj)) {
      Rcpp::stop("LogicalRules OR_DETECTION_ANALYTIC requires finite UC for upper-censored trials.");
    }
  } else if (is_gng) {
    // Overt (go) responses are "yes"; a withheld response is rt = +Inf with a
    // missing or "no" identity. UC may be +Inf (the withheld mass extends to
    // infinity when there is no upper deadline).
    if (finite_rt && resp != 1) return min_ll;
    if (upper_censored && resp != 0 && resp != 2) return min_ll;
    if (lower_censored && resp != 0 && resp != 1) return min_ll;
  } else {
    if (!has_channels) return min_ll;
    if (upper_censored && !R_FINITE(UCj)) return min_ll;
  }

  // ---- Node-invariant evaluation-point caches -----------------------------
  // Detection points.
  LrCapDetPoint det_t, det_uc, det_ut, det_lt, det_clo, det_chi, det_gowin;
  // Choice channel points.
  LrCapChannelPoint A_t, B_t, A_uc, B_uc, A_ut, B_ut, A_lt, B_lt,
    A_clo, B_clo, A_chi, B_chi;

  const double clo = std::max(0.0, LTj);
  const double chi = std::max(clo, LCj);

  if (is_detection) {
    if (finite_rt) {
      lrcap_build_det_N_point(det_t, tvA, tvB, parN, has_nogo, t, false,
                              pdf1, cdf1, ctx, min_surv);
    } else if (upper_censored) {
      lrcap_build_det_N_point(det_uc, tvA, tvB, parN, has_nogo, UCj,
                              rule_code == 6, pdf1, cdf1, ctx, min_surv);
    } else {
      if (rule_code == 5) {
        lrcap_build_det_N_point(det_clo, tvA, tvB, parN, false, clo, false,
                                pdf1, cdf1, ctx, min_surv);
        lrcap_build_det_N_point(det_chi, tvA, tvB, parN, false, chi, false,
                                pdf1, cdf1, ctx, min_surv);
      } else {
        lrcap_build_det_gowin_point(det_gowin, tvA, tvB, parN, has_nogo,
                                    clo, chi, cdf1, ctx, min_surv);
      }
    }
    const bool need_ut_for_event = upper_censored && subtract_UT;
    if (has_trunc && (!z_cache_hit || need_ut_for_event)) {
      if (!z_cache_hit && LTj > 0.0) {
        lrcap_build_det_N_point(det_lt, tvA, tvB, parN, has_nogo, LTj,
                                rule_code == 6, pdf1, cdf1, ctx, min_surv);
      }
      if (R_FINITE(UTj) && (!z_cache_hit || need_ut_for_event)) {
        lrcap_build_det_N_point(det_ut, tvA, tvB, parN, has_nogo, UTj,
                                rule_code == 6, pdf1, cdf1, ctx, min_surv);
      }
    }
  } else {
    // GNG builds only the finite-RT go-density channel points; its censor,
    // withheld and truncation masses are conditional "no overt go by t"
    // probabilities evaluated on the fly by gng_N_at (which handles UC = +Inf).
    if (finite_rt) {
      lrcap_build_channel_point(A_t, tvA, parnA, t, true, true,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_t, tvB, parnB, t, true, true,
                                pdf1, cdf1, ctx, min_surv);
    } else if (!is_gng && upper_censored) {
      lrcap_build_channel_point(A_uc, tvA, parnA, UCj, need_G_for_N, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_uc, tvB, parnB, UCj, need_G_for_N, false,
                                pdf1, cdf1, ctx, min_surv);
    } else if (!is_gng) {
      lrcap_build_channel_point(A_clo, tvA, parnA, clo, true, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_clo, tvB, parnB, clo, true, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(A_chi, tvA, parnA, chi, true, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_chi, tvB, parnB, chi, true, false,
                                pdf1, cdf1, ctx, min_surv);
    }
    const bool need_ut_for_event = upper_censored && subtract_UT;
    if (!is_gng && has_trunc && (!z_cache_hit || need_ut_for_event)) {
      if (!z_cache_hit && LTj > 0.0) {
        lrcap_build_channel_point(A_lt, tvA, parnA, LTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
        lrcap_build_channel_point(B_lt, tvB, parnB, LTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
      }
      if (R_FINITE(UTj) && (!z_cache_hit || need_ut_for_event)) {
        lrcap_build_channel_point(A_ut, tvA, parnA, UTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
        lrcap_build_channel_point(B_ut, tvB, parnB, UTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
      }
    }
  }

  // Conditional "no overt go response by tt" probability for GNG at latent z,
  // N_GNG(tt | z) = Q_A(tt|z) Q_B(tt|z), Q_X = 1 - GX_yes(tt|z).  The
  // conditional target drift mean is materialized per node (v -> v0 + slope*z)
  // and the channel win probability is integrated with the shared scalar GSL
  // helper, which also handles tt = +Inf (the withheld-to-infinity case).
  double bufA_row[64];
  double bufB_row[64];
  int gng_isok2[2] = {1, 1};
  // gng_N_at is used both for the censor/withheld event masses (!finite_rt) and
  // for the truncation-window denominator (any finite-RT GNG trial with LT/UT),
  // so initialize its per-node buffers whenever the trial is GNG.
  if (is_gng) {
    for (int p = 0; p < n_par; ++p) {
      bufA_row[p] = pars_cols[p][idxA];
      bufB_row[p] = pars_cols[p][idxB];
    }
    gng_scratch.resize(static_cast<size_t>(2 * n_par));
  }
  auto gng_N_at = [&](double tt, double z, bool& ok) -> double {
    ok = true;
    if (!(tt > 0.0)) return 1.0;
    bufA_row[emc2col::lba::v] = tvA.v0 + tvA.slope * z;
    bufB_row[emc2col::lba::v] = tvB.v0 + tvB.slope * z;
    const double gA = lr_channel_yes_prob(tt, bufA_row, parnA, n_par, ctx,
                                          pdf1, cdf1, gsl_ctl, w,
                                          gng_scratch.data(), gng_isok2);
    const double gB = lr_channel_yes_prob(tt, bufB_row, parnB, n_par, ctx,
                                          pdf1, cdf1, gsl_ctl, w,
                                          gng_scratch.data(), gng_isok2);
    if (!R_FINITE(gA) || !R_FINITE(gB)) { ok = false; return 1.0; }
    const double QA = std::min(1.0, std::max(0.0, 1.0 - gA));
    const double QB = std::min(1.0, std::max(0.0, 1.0 - gB));
    return QA * QB;
  };

  // Conditional no-response probability for the choice rules at a prepared
  // channel-point pair (XOR/ID need only the channel survivors).
  auto choice_N_at = [&](const LrCapChannelPoint& cpA,
                         const LrCapChannelPoint& cpB, double z,
                         double lqA, double lqB, bool& ok) -> double {
    if (cpA.at_zero && cpB.at_zero) return 1.0;
    const LrCapChanVals va = lrcap_chan_vals(cpA, z, lqA, posdrift, false, count);
    if (!va.ok) { ok = false; return 1.0; }
    const LrCapChanVals vb = lrcap_chan_vals(cpB, z, lqB, posdrift, false, count);
    if (!vb.ok) { ok = false; return 1.0; }
    if (rule_code == 3 || rule_code == 4) {
      const double SdA = va.S_t * va.S_n;
      const double SdB = vb.S_t * vb.S_n;
      return 1.0 - (1.0 - SdA) * (1.0 - SdB);
    }
    return lr_no_response_prob(rule_code, lrcap_state_from_vals(va),
                               lrcap_state_from_vals(vb));
  };

  // ---- Conditional event mass (numerator integrand, includes the node
  // positivity weight qA(z)*qB(z)) ------------------------------------------
  auto log_event_mass = [&](double z) -> double {
    const double lqA = lrcap_log_q(tvA, z, posdrift);
    const double lqB = lrcap_log_q(tvB, z, posdrift);
    if (!R_FINITE(lqA) || !R_FINITE(lqB)) return R_NegInf;
    bool ok = true;
    double p = 0.0;
    if (is_detection) {
      if (finite_rt) {
        const double fA = lrcap_cond_pdf(det_t.prepA, z, lqA, ok);
        if (!ok) return R_NegInf;
        const double SA = lrcap_cond_survivor(det_t.prepA, z, lqA, posdrift, ok);
        if (!ok) return R_NegInf;
        const double fB = lrcap_cond_pdf(det_t.prepB, z, lqB, ok);
        if (!ok) return R_NegInf;
        const double SB = lrcap_cond_survivor(det_t.prepB, z, lqB, posdrift, ok);
        if (!ok) return R_NegInf;
        p = fA * std::max(min_surv, SB) + fB * std::max(min_surv, SA);
        if (rule_code == 6) p *= det_t.S_N;
      } else if (upper_censored) {
        p = lrcap_det_N_at(det_uc, rule_code, z, lqA, lqB, posdrift, ok, count);
        if (ok && subtract_UT) {
          const double N_UT = lrcap_det_N_at(det_ut, rule_code, z, lqA, lqB,
                                             posdrift, ok, count);
          p = std::max(0.0, p - N_UT);
        }
        if (!ok) return R_NegInf;
      } else {
        if (rule_code == 5) {
          const double N_lo = lrcap_det_N_at(det_clo, 5, z, lqA, lqB,
                                             posdrift, ok, count);
          if (!ok) return R_NegInf;
          const double N_hi = lrcap_det_N_at(det_chi, 5, z, lqA, lqB,
                                             posdrift, ok, count);
          if (!ok) return R_NegInf;
          p = std::max(0.0, N_lo - N_hi);
        } else {
          p = lrcap_det_gowin_at(det_gowin, z, lqA, lqB, posdrift, ok, count);
          if (!ok) return R_NegInf;
        }
      }
    } else {
      if (finite_rt) {
        const LrCapChanVals va = lrcap_chan_vals(A_t, z, lqA, posdrift, true, count);
        if (!va.ok) return R_NegInf;
        const LrCapChanVals vb = lrcap_chan_vals(B_t, z, lqB, posdrift, true, count);
        if (!vb.ok) return R_NegInf;
        const double one_m_FA = std::max(min_surv, va.S_t);
        const double one_m_FB = std::max(min_surv, vb.S_t);
        const double one_m_FnA = std::max(min_surv, va.S_n);
        const double one_m_FnB = std::max(min_surv, vb.S_n);
        const double gA_yes = va.f_t * one_m_FnA;
        const double gB_yes = vb.f_t * one_m_FnB;
        const double gA_no = va.f_n * one_m_FA;
        const double gB_no = vb.f_n * one_m_FB;
        const double S_dec_A = one_m_FA * one_m_FnA;
        const double S_dec_B = one_m_FB * one_m_FnB;
        const double GA_no = va.G_no;
        const double GB_no = vb.G_no;
        if (rule_code == 1 || is_gng) {
          // GNG shares OR's go-response density (first channel-yes); its "no"
          // outcome is withheld (rt = +Inf), never a finite-RT event, so only
          // resp == 1 (gated above) contributes here.
          const double s_GA_yes = std::min(1.0, std::max(min_surv, GA_no + S_dec_A));
          const double s_GB_yes = std::min(1.0, std::max(min_surv, GB_no + S_dec_B));
          if (resp == 1) p = gA_yes * s_GB_yes + gB_yes * s_GA_yes;
          else if (resp == 2 && !is_gng) p = gA_no * GB_no + gB_no * GA_no;
        } else if (rule_code == 2) {
          const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
          const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
          const double s_GA_no = std::min(1.0, std::max(min_surv, GA_yes + S_dec_A));
          const double s_GB_no = std::min(1.0, std::max(min_surv, GB_yes + S_dec_B));
          if (resp == 1) p = gA_yes * GB_yes + gB_yes * GA_yes;
          else if (resp == 2) p = gA_no * s_GB_no + gB_no * s_GA_no;
        } else {
          const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
          const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
          const double dAB = gA_yes * GB_yes + gB_yes * GA_yes;
          const double dAN = gA_yes * GB_no + gB_no * GA_yes;
          const double dNB = gA_no * GB_yes + gB_yes * GA_no;
          const double dNN = gA_no * GB_no + gB_no * GA_no;
          if (rule_code == 4) {
            if (resp == 6) p = dAB;
            else if (resp == 4) p = dAN;
            else if (resp == 5) p = dNB;
            else if (resp == 3) p = dNN;
          } else {
            if (resp == 1) p = dAN + dNB;
            else if (resp == 2) p = dAB + dNN;
          }
        }
      } else if (is_gng) {
        // Withheld (rt = +Inf): no overt go by UC, N_GNG(UC|z) minus N_GNG(UT|z)
        // when the window has a finite upper bound; lower censor: go mass in
        // [clo, chi] = N_GNG(clo|z) - N_GNG(chi|z).
        if (upper_censored) {
          p = gng_N_at(UCj, z, ok);
          if (ok && subtract_UT) {
            const double N_UT = gng_N_at(UTj, z, ok);
            p = std::max(0.0, p - N_UT);
          }
          if (!ok) return R_NegInf;
        } else {
          const double N_lo = gng_N_at(clo, z, ok);
          if (!ok) return R_NegInf;
          const double N_hi = gng_N_at(chi, z, ok);
          if (!ok) return R_NegInf;
          p = std::max(0.0, N_lo - N_hi);
        }
      } else if (upper_censored) {
        p = choice_N_at(A_uc, B_uc, z, lqA, lqB, ok);
        if (ok && subtract_UT) {
          const double N_UT = choice_N_at(A_ut, B_ut, z, lqA, lqB, ok);
          p = std::max(0.0, p - N_UT);
        }
        if (!ok) return R_NegInf;
      } else {
        const LrCapChanVals va_lo = lrcap_chan_vals(A_clo, z, lqA, posdrift, false, count);
        const LrCapChanVals vb_lo = lrcap_chan_vals(B_clo, z, lqB, posdrift, false, count);
        const LrCapChanVals va_hi = lrcap_chan_vals(A_chi, z, lqA, posdrift, false, count);
        const LrCapChanVals vb_hi = lrcap_chan_vals(B_chi, z, lqB, posdrift, false, count);
        if (!va_lo.ok || !vb_lo.ok || !va_hi.ok || !vb_hi.ok) return R_NegInf;
        const double cdf_hi = lr_response_cdf(rule_code, resp,
                                              lrcap_state_from_vals(va_hi),
                                              lrcap_state_from_vals(vb_hi));
        const double cdf_lo = lr_response_cdf(rule_code, resp,
                                              lrcap_state_from_vals(va_lo),
                                              lrcap_state_from_vals(vb_lo));
        if (ISNAN(cdf_hi) || ISNAN(cdf_lo)) return R_NegInf;
        p = std::max(0.0, cdf_hi - cdf_lo);
      }
    }
    if (!(p > 0.0) || !R_FINITE(p)) return R_NegInf;
    return std::log(p) + lqA + lqB;
  };

  // ---- Conditional truncation-window mass (denominator integrand) ---------
  auto log_trunc_mass = [&](double z) -> double {
    const double lqA = lrcap_log_q(tvA, z, posdrift);
    const double lqB = lrcap_log_q(tvB, z, posdrift);
    if (!R_FINITE(lqA) || !R_FINITE(lqB)) return R_NegInf;
    bool ok = true;
    double N_LT = 1.0;
    double N_UT = 0.0;
    if (is_detection) {
      if (LTj > 0.0) {
        N_LT = lrcap_det_N_at(det_lt, rule_code, z, lqA, lqB, posdrift, ok, count);
        if (!ok) return R_NegInf;
      }
      if (R_FINITE(UTj)) {
        N_UT = lrcap_det_N_at(det_ut, rule_code, z, lqA, lqB, posdrift, ok, count);
        if (!ok) return R_NegInf;
      }
    } else if (is_gng) {
      if (LTj > 0.0) {
        N_LT = gng_N_at(LTj, z, ok);
        if (!ok) return R_NegInf;
      }
      if (R_FINITE(UTj)) {
        N_UT = gng_N_at(UTj, z, ok);
        if (!ok) return R_NegInf;
      }
    } else {
      if (LTj > 0.0) {
        N_LT = choice_N_at(A_lt, B_lt, z, lqA, lqB, ok);
        if (!ok) return R_NegInf;
      }
      if (R_FINITE(UTj)) {
        N_UT = choice_N_at(A_ut, B_ut, z, lqA, lqB, ok);
        if (!ok) return R_NegInf;
      }
    }
    const double Z = std::min(1.0, N_LT - N_UT);
    if (!(Z > 0.0) || !R_FINITE(Z)) return R_NegInf;
    return std::log(Z) + lqA + lqB;
  };

  double log_num = R_NegInf;
  double log_den = 0.0;
  if (tau == 0.0) {
    // The factor is degenerate: one node at z = 0 is the exact integral, and
    // the positivity weight cancels exactly against the rho = 0 orthant.
    if (count) ++lr_capacity_counters().tau_zero_trials;
    log_num = log_event_mass(0.0);
    const double log_w0 = lrcap_log_q(tvA, 0.0, posdrift) +
      lrcap_log_q(tvB, 0.0, posdrift);
    log_den = has_trunc
      ? (z_cache_hit ? cached_log_den : log_trunc_mass(0.0))
      : log_w0;
  } else {
    // All censored event masses are probabilities; finite-RT event densities
    // are not bounded by one and must retain the ordinary adaptive result.
    log_num = lrcap_integrate_factor(log_event_mass, n_scan, n_fine, count,
                                     !finite_rt);
    log_den = has_trunc
      ? (z_cache_hit
          ? cached_log_den
          : lrcap_integrate_factor(log_trunc_mass, n_scan, n_fine, count,
                                   true))
      : log_qAB;
  }
  if (has_trunc && !z_cache_hit && cell_log_den_cap != nullptr && cell_idx >= 0 &&
      cell_idx < static_cast<int>(cell_log_den_cap->size())) {
    const bool valid_z = R_FINITE(log_den);
    (*cell_log_den_cap)[cell_idx] = valid_z ? log_den : R_NegInf;
  }
  if (!R_FINITE(log_den)) return min_ll;
  const double ll = log_num - log_den;
  if (!finite_rt && (!R_FINITE(log_num) || log_num > 1e-8 ||
                     log_den > 1e-8 || ll > 1e-8)) return min_ll;
  return (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
}

static double c_log_likelihood_logicalrules(
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
    Rcpp::NumericVector* trial_ll_out,
    int kappa_col,
    int tau_col,
    int pc_col,
    int pg_col,
    const GuessKernel* guess) {
  if (!shared.valid || n_acc <= 0 || pdf1 == nullptr || cdf1 == nullptr || model_ctx == nullptr) {
    Rcpp::stop("c_log_likelihood_logicalrules: invalid logical-rules configuration.");
  }
  const int n_trials = shared.n_trials;
  if (ok_params.size() != n_trials) {
    Rcpp::stop("c_log_likelihood_logicalrules: ok_params size must match shared data.");
  }
  if (n_trials == 0) return 0.0;
  // ParamTable columns are refilled in place for each particle, so no cache
  // entry may survive this likelihood call.
  shared.clear_particle_cache();
  if (model_ctx->fpe_cache) model_ctx->fpe_cache->new_particle();
  if (model_ctx->rlf_cache) model_ctx->rlf_cache->new_particle();
  if (n_par > 64) {
    Rcpp::stop("c_log_likelihood_logicalrules: at most 64 parameter columns are supported.");
  }
  const bool capacity = kappa_col >= 0 || tau_col >= 0;
  if (capacity && (kappa_col < 0 || tau_col < 0 ||
                   kappa_col >= n_par || tau_col >= n_par)) {
    Rcpp::stop("c_log_likelihood_logicalrules: invalid LogicalRules capacity columns.");
  }
  const bool count_capacity = capacity && lr_capacity_counters_enabled();
  const int lrcap_scan_n = capacity
    ? emc2_quad_nodes("EMC2_LRCAP_SCAN_N", 12) : 0;
  const int lrcap_fine_n = capacity
    ? emc2_quad_nodes("EMC2_LRCAP_FINE_N", 12) : 0;

  static const double kMinSurv = 1e-300;
  const int n_unique_trials = shared.n_unique_trials;
  static thread_local LogicalRulesScratch scratch;
  resize_assign_double(scratch.ll_unique, static_cast<size_t>(n_unique_trials), min_ll);
  std::vector<double>& ll_unique = scratch.ll_unique;
  const bool use_raw_local =
    (model_dfun_raw != nullptr &&
     model_pfun_raw != nullptr &&
     static_cast<int>(shared.rt_by_row.size()) == n_trials);

  GslIntegrationControls gsl_ctl = default_gsl_controls();
  gsl_ctl.try_qng_first_finite = true;
  gsl_ctl.qag_key = GSL_INTEG_GAUSS21;
  gsl_ctl.rel_tol = 1e-5;  // 1e-4 let the adaptive rule's optimistic error estimate through

  // Change 3: skip the raw batch (and all integration) if every parameter row is invalid.
  bool any_ok = false;
  for (int i = 0; i < n_trials && !any_ok; ++i) any_ok = (bool)ok_params[i];
  if (!any_ok) {
    const int n_out = expand.length() > 0 ? (int)expand.length() : n_unique_trials;
    if (trial_ll_out != nullptr) for (int i = 0; i < n_out; ++i) (*trial_ll_out)[i] = min_ll;
    return static_cast<double>(n_out) * min_ll;
  }

  std::vector<double>& logf_all = scratch.logf_all;
  std::vector<double>& logS_all = scratch.logS_all;
  std::vector<int>& all_mask = scratch.all_mask;
  std::vector<int>& isok_int_all = scratch.isok_int_all;
  if (use_raw_local) {
    resize_assign_double(logf_all, static_cast<size_t>(n_trials), min_ll);
    resize_assign_double(logS_all, static_cast<size_t>(n_trials), min_ll);
    resize_assign_int(all_mask, static_cast<size_t>(n_trials), 1);
    isok_int_all.resize(static_cast<size_t>(n_trials));
    for (int i = 0; i < n_trials; ++i) isok_int_all[static_cast<size_t>(i)] = ok_params[i] ? 1 : 0;
    model_dfun_raw(shared.rt_by_row.data(), pars_cols, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logf_all.data(), min_ll, model_ctx);
    model_pfun_raw(shared.rt_by_row.data(), pars_cols, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logS_all.data(), min_ll, model_ctx);
  }

  // ---- Gauss-Legendre batch integration pre-pass ----
  // Replaces per-trial local_race_helper (serial GSL) with N_GL vectorised batch sweeps.
  // GA_no_gl[j] = P(nontarget-A wins channel-A sub-race before RT_j)  [linear probability]
  // GB_no_gl[j] = same for channel B, or = GA_no_gl[j] when channels are equal.
  // OR/AND rules use GA_no_gl directly; XOR/ID rules derive GA_yes = 1 - GA_no_gl - S_dec_A.
  const bool use_gl_pass = use_raw_local;
  std::vector<double>& GA_no_gl = scratch.GA_no_gl;
  std::vector<double>& GB_no_gl = scratch.GB_no_gl;
  std::vector<unsigned char>& ch_eq_vec = scratch.ch_eq_vec;

  if (use_gl_pass) {
    const GLRule& gl_rule31 = gl_get_rule(31);   // shared cache, gl_quad.h

    // Per-trial GL parameters for GA_no = P(n_A wins before RT) and
    // GB_no = P(n_B wins before RT).  The integration lower bound is the
    // density accumulator's t0, not its opponent's t0; otherwise unequal t0
    // pairs lose early mass where the opponent survivor is exactly one.
    const int t0_col = model_ctx->t0_index;
    std::vector<double>& gl_h_A = scratch.gl_h_A;
    std::vector<double>& gl_h_B = scratch.gl_h_B;
    std::vector<double>& gl_lo_A = scratch.gl_lo_A;
    std::vector<double>& gl_lo_B = scratch.gl_lo_B;
    std::vector<unsigned char>& pair_ok_A = scratch.pair_ok_A;
    std::vector<unsigned char>& pair_ok_B = scratch.pair_ok_B;
    gl_h_A.resize(static_cast<size_t>(n_unique_trials));
    gl_h_B.resize(static_cast<size_t>(n_unique_trials));
    gl_lo_A.resize(static_cast<size_t>(n_unique_trials));
    gl_lo_B.resize(static_cast<size_t>(n_unique_trials));
    ch_eq_vec.resize(static_cast<size_t>(n_unique_trials));
    pair_ok_A.resize(static_cast<size_t>(n_unique_trials));
    pair_ok_B.resize(static_cast<size_t>(n_unique_trials));
    bool any_unequal = false;
    for (int j = 0; j < n_unique_trials; ++j) {
      const int iA  = shared.idxA[j],  inA = shared.idxnA[j];
      const int iB  = shared.idxB[j],  inB = shared.idxnB[j];
      const double t0nA = (t0_col >= 0 && t0_col < n_par && inA >= 0)
                        ? pars_cols[t0_col][inA] : 0.0;
      const double t0nB = (t0_col >= 0 && t0_col < n_par && inB >= 0)
                        ? pars_cols[t0_col][inB] : 0.0;
      const double RT = shared.rt_unique[j];
      const double loA = std::max(0.0, t0nA);
      const double loB = std::max(0.0, t0nB);
      gl_lo_A[static_cast<size_t>(j)] = loA;
      gl_lo_B[static_cast<size_t>(j)] = loB;
      gl_h_A[static_cast<size_t>(j)] = (RT > loA) ? (RT - loA) * 0.5 : 0.0;
      gl_h_B[static_cast<size_t>(j)] = (RT > loB) ? (RT - loB) * 0.5 : 0.0;
      ch_eq_vec[j] = (inA >= 0 && inB >= 0) &&
                     row_equal_colmajor(pars_cols, n_par, iA, iB) &&
                     row_equal_colmajor(pars_cols, n_par, inA, inB);
      if (!ch_eq_vec[j]) any_unequal = true;
      pair_ok_A[static_cast<size_t>(j)] =
        (iA >= 0 && inA >= 0 && ok_params[iA] && ok_params[inA] &&
         RT > 0.0 && R_FINITE(RT) && gl_h_A[static_cast<size_t>(j)] > 0.0) ? 1 : 0;
      pair_ok_B[static_cast<size_t>(j)] =
        (iB >= 0 && inB >= 0 && ok_params[iB] && ok_params[inB] &&
         RT > 0.0 && R_FINITE(RT) && gl_h_B[static_cast<size_t>(j)] > 0.0) ? 1 : 0;
    }

    // Compact column-major parameter matrices (n_unique_trials rows) for each accumulator role.
    // Avoids iterating over the full n_trials parameter matrix with sparse masks.
    const size_t cpt = static_cast<size_t>(n_unique_trials * n_par);
    std::vector<double>& pars_nA = scratch.pars_nA;
    std::vector<double>& pars_A = scratch.pars_A;
    std::vector<double>& pars_nB = scratch.pars_nB;
    std::vector<double>& pars_B = scratch.pars_B;
    std::vector<int>& isok_nA = scratch.isok_nA;
    std::vector<int>& isok_A = scratch.isok_A;
    std::vector<int>& isok_nB = scratch.isok_nB;
    std::vector<int>& isok_B = scratch.isok_B;
    pars_nA.resize(cpt);
    pars_A.resize(cpt);
    isok_nA.resize(static_cast<size_t>(n_unique_trials));
    isok_A.resize(static_cast<size_t>(n_unique_trials));
    if (any_unequal) {
      pars_nB.resize(cpt);
      pars_B.resize(cpt);
      isok_nB.resize(static_cast<size_t>(n_unique_trials));
      isok_B.resize(static_cast<size_t>(n_unique_trials));
    }
    for (int p = 0; p < n_par; ++p) {
      const double* src = pars_cols[p];
      double* d_nA = pars_nA.data() + static_cast<size_t>(p) * n_unique_trials;
      double* d_A  = pars_A.data()  + static_cast<size_t>(p) * n_unique_trials;
      for (int j = 0; j < n_unique_trials; ++j) {
        const int inA_j = shared.idxnA[j];
        d_nA[j] = (inA_j >= 0) ? src[inA_j] : 0.0;
        d_A[j]  = src[shared.idxA[j]];
      }
      if (any_unequal) {
        double* d_nB = pars_nB.data() + static_cast<size_t>(p) * n_unique_trials;
        double* d_B  = pars_B.data()  + static_cast<size_t>(p) * n_unique_trials;
        for (int j = 0; j < n_unique_trials; ++j) {
          const int inB_j = shared.idxnB[j];
          d_nB[j] = (inB_j >= 0) ? src[inB_j] : 0.0;
          d_B[j]  = src[shared.idxB[j]];
        }
      }
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      isok_nA[j] = pair_ok_A[static_cast<size_t>(j)] ? isok_int_all[shared.idxnA[j]] : 0;
      isok_A[j]  = pair_ok_A[static_cast<size_t>(j)] ? isok_int_all[shared.idxA[j]] : 0;
      if (any_unequal) {
        isok_nB[j] = (!ch_eq_vec[j] && pair_ok_B[static_cast<size_t>(j)]) ? isok_int_all[shared.idxnB[j]] : 0;
        isok_B[j]  = (!ch_eq_vec[j] && pair_ok_B[static_cast<size_t>(j)]) ? isok_int_all[shared.idxB[j]] : 0;
      }
    }

    // Column-pointer views of the compact role matrices for the raw kernels.
    // Slots beyond n_par stay nullptr: kernels may fetch (not dereference)
    // optional trailing columns this model variant lacks.
    const double* colsp_nA[64] = {nullptr};
    const double* colsp_A[64]  = {nullptr};
    const double* colsp_nB[64] = {nullptr};
    const double* colsp_B[64]  = {nullptr};
    for (int p = 0; p < n_par; ++p) {
      colsp_nA[p] = pars_nA.data() + static_cast<size_t>(p) * n_unique_trials;
      colsp_A[p]  = pars_A.data()  + static_cast<size_t>(p) * n_unique_trials;
      if (any_unequal) {
        colsp_nB[p] = pars_nB.data() + static_cast<size_t>(p) * n_unique_trials;
        colsp_B[p]  = pars_B.data()  + static_cast<size_t>(p) * n_unique_trials;
      }
    }

    resize_assign_double(GA_no_gl, static_cast<size_t>(n_unique_trials), 0.0);
    resize_assign_double(GB_no_gl, static_cast<size_t>(n_unique_trials), 0.0);

    std::vector<int>& all1 = scratch.all1;
    resize_assign_int(all1, static_cast<size_t>(n_unique_trials), 1);
    std::vector<double>& gl_rt = scratch.tmpA1;
    std::vector<double>& lf_nA = scratch.tmpA2;
    std::vector<double>& lS_A = scratch.tmpB1;
    std::vector<double>& lf_nB = scratch.tmpB2;
    std::vector<double>& lS_B = scratch.parA;
    gl_rt.resize(static_cast<size_t>(n_unique_trials));
    resize_assign_double(lf_nA, static_cast<size_t>(n_unique_trials), min_ll);
    resize_assign_double(lS_A, static_cast<size_t>(n_unique_trials), min_ll);
    resize_assign_double(lf_nB, static_cast<size_t>(n_unique_trials), min_ll);
    resize_assign_double(lS_B, static_cast<size_t>(n_unique_trials), min_ll);

    for (size_t k = 0; k < gl_rule31.x.size(); ++k) {
      const double xi = gl_rule31.x[k], wt = gl_rule31.w[k];
      #pragma omp simd
      for (int j = 0; j < n_unique_trials; ++j) {
        gl_rt[static_cast<size_t>(j)] = gl_lo_A[static_cast<size_t>(j)]
                                      + gl_h_A[static_cast<size_t>(j)] * (1.0 + xi);
      }
      model_dfun_raw(gl_rt.data(), colsp_nA, n_unique_trials,
                     all1.data(), isok_nA.data(), lf_nA.data(), min_ll, model_ctx);
      model_pfun_raw(gl_rt.data(), colsp_A,  n_unique_trials,
                     all1.data(), isok_A.data(),  lS_A.data(),  min_ll, model_ctx);
      #pragma omp simd
      for (int j = 0; j < n_unique_trials; ++j)
        GA_no_gl[j] += wt * gl_h_A[j] * std::exp(lf_nA[j] + lS_A[j]);
      if (any_unequal) {
        #pragma omp simd
        for (int j = 0; j < n_unique_trials; ++j) {
          gl_rt[static_cast<size_t>(j)] = gl_lo_B[static_cast<size_t>(j)]
                                        + gl_h_B[static_cast<size_t>(j)] * (1.0 + xi);
        }
        model_dfun_raw(gl_rt.data(), colsp_nB, n_unique_trials,
                       all1.data(), isok_nB.data(), lf_nB.data(), min_ll, model_ctx);
        model_pfun_raw(gl_rt.data(), colsp_B,  n_unique_trials,
                       all1.data(), isok_B.data(),  lS_B.data(),  min_ll, model_ctx);
        #pragma omp simd
        for (int j = 0; j < n_unique_trials; ++j)
          GB_no_gl[j] += wt * gl_h_B[j] * std::exp(lf_nB[j] + lS_B[j]);
      }
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      if (ch_eq_vec[static_cast<size_t>(j)]) GB_no_gl[static_cast<size_t>(j)] = GA_no_gl[static_cast<size_t>(j)];
    }

    // ---- Auxiliary GL sweeps: G_no at censoring/truncation time points ----
    // Same machinery as the RT sweep above, evaluated at the per-trial times
    // the trial loop needs: the truncation normaliser endpoints (LT, UT), the
    // lower-censor window ends, and UC for upper-censored OR/AND trials.
    // Values left NaN are recomputed by the scalar GSL helper in the trial
    // loop, so a slot missed here costs speed, never correctness.
    std::vector<double>& aux_t  = scratch.aux_t;
    std::vector<double>& aux_hA = scratch.aux_hA;
    std::vector<double>& aux_hB = scratch.aux_hB;
    std::vector<double>& aux_GA = scratch.aux_GA;
    std::vector<double>& aux_GB = scratch.aux_GB;
    std::vector<int>& aux_mask  = scratch.aux_mask;
    const size_t aux_len = static_cast<size_t>(LR_AUX_N) * n_unique_trials;
    aux_t.assign(aux_len, NA_REAL);
    aux_hA.assign(aux_len, 0.0);
    aux_hB.assign(aux_len, 0.0);
    aux_GA.assign(aux_len, NA_REAL);
    aux_GB.assign(aux_len, NA_REAL);
    aux_mask.assign(aux_len, 0);
    bool slot_any[LR_AUX_N] = {false, false, false, false, false};

    std::vector<int>& isok_nA_b = scratch.isok_nA_base;
    std::vector<int>& isok_A_b  = scratch.isok_A_base;
    std::vector<int>& isok_nB_b = scratch.isok_nB_base;
    std::vector<int>& isok_B_b  = scratch.isok_B_base;
    resize_assign_int(isok_nA_b, static_cast<size_t>(n_unique_trials), 0);
    resize_assign_int(isok_A_b,  static_cast<size_t>(n_unique_trials), 0);
    resize_assign_int(isok_nB_b, static_cast<size_t>(n_unique_trials), 0);
    resize_assign_int(isok_B_b,  static_cast<size_t>(n_unique_trials), 0);

    bool any_aux = false;
    for (int j = 0; j < n_unique_trials; ++j) {
      const int rc = shared.rule_code[static_cast<size_t>(j)];
      if (rc == 5) continue;   // analytic detection keeps its own route; GNG (6) uses four horses
      const int iA = shared.idxA[static_cast<size_t>(j)], inA = shared.idxnA[static_cast<size_t>(j)];
      const int iB = shared.idxB[static_cast<size_t>(j)], inB = shared.idxnB[static_cast<size_t>(j)];
      const bool base_ok =
        iA >= 0 && inA >= 0 && iB >= 0 && inB >= 0 &&
        ok_params[iA] && ok_params[inA] && ok_params[iB] && ok_params[inB];
      if (!base_ok) continue;
      isok_nA_b[static_cast<size_t>(j)] = 1;
      isok_A_b[static_cast<size_t>(j)]  = 1;
      if (!ch_eq_vec[static_cast<size_t>(j)]) {
        isok_nB_b[static_cast<size_t>(j)] = 1;
        isok_B_b[static_cast<size_t>(j)]  = 1;
      }

      const bool needs_G = (rc == 1 || rc == 2 || rc == 6);   // XOR/ID N(t) is survivor-only
      const double rtj = shared.rt_unique[static_cast<size_t>(j)];
      const double LTj = shared.LT_unique[static_cast<size_t>(j)];
      const double UTj = shared.UT_unique[static_cast<size_t>(j)];
      const bool truncated = (LTj != 0.0 || R_FINITE(UTj));

      double slot_t[LR_AUX_N] = {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL};
      if (needs_G && truncated && LTj > 0.0) slot_t[LR_AUX_LT] = LTj;
      if (needs_G && R_FINITE(UTj))          slot_t[LR_AUX_UT] = UTj;
      if (rtj == R_NegInf) {   // lower censor: response CDFs need G for all rules <= 4
        const double lo = std::max(0.0, LTj);
        const double hi = std::max(lo, shared.LC_unique[static_cast<size_t>(j)]);
        if (lo > 0.0) slot_t[LR_AUX_CLO] = lo;
        if (hi > 0.0) slot_t[LR_AUX_CHI] = hi;
      } else if (rtj == R_PosInf && needs_G) {
        const double UCj = shared.UC_unique[static_cast<size_t>(j)];
        if (R_FINITE(UCj)) slot_t[LR_AUX_UC] = UCj;
      }
      for (int s = 0; s < LR_AUX_N; ++s) {
        if (ISNAN(slot_t[s])) continue;
        const size_t o = static_cast<size_t>(s) * n_unique_trials + j;
        aux_t[o] = slot_t[s];
        aux_hA[o] = (slot_t[s] > gl_lo_A[static_cast<size_t>(j)])
                  ? (slot_t[s] - gl_lo_A[static_cast<size_t>(j)]) * 0.5 : 0.0;
        aux_hB[o] = (slot_t[s] > gl_lo_B[static_cast<size_t>(j)])
                  ? (slot_t[s] - gl_lo_B[static_cast<size_t>(j)]) * 0.5 : 0.0;
        aux_mask[o] = 1;
        aux_GA[o] = 0.0;
        aux_GB[o] = 0.0;
        slot_any[s] = true;
        any_aux = true;
      }
    }

    if (any_aux) {
      for (int s = 0; s < LR_AUX_N; ++s) {
        if (!slot_any[s]) continue;
        const double* hA = aux_hA.data() + static_cast<size_t>(s) * n_unique_trials;
        const double* hB = aux_hB.data() + static_cast<size_t>(s) * n_unique_trials;
        const int* msk   = aux_mask.data() + static_cast<size_t>(s) * n_unique_trials;
        double* GA = aux_GA.data() + static_cast<size_t>(s) * n_unique_trials;
        double* GB = aux_GB.data() + static_cast<size_t>(s) * n_unique_trials;
        for (size_t k = 0; k < gl_rule31.x.size(); ++k) {
          const double xi = gl_rule31.x[k], wt = gl_rule31.w[k];
          #pragma omp simd
          for (int j = 0; j < n_unique_trials; ++j) {
            gl_rt[static_cast<size_t>(j)] = gl_lo_A[static_cast<size_t>(j)]
                                          + hA[j] * (1.0 + xi);
          }
          model_dfun_raw(gl_rt.data(), colsp_nA, n_unique_trials,
                         msk, isok_nA_b.data(), lf_nA.data(), min_ll, model_ctx);
          model_pfun_raw(gl_rt.data(), colsp_A,  n_unique_trials,
                         msk, isok_A_b.data(),  lS_A.data(),  min_ll, model_ctx);
          for (int j = 0; j < n_unique_trials; ++j) {
            if (msk[j]) GA[j] += wt * hA[j] * std::exp(lf_nA[j] + lS_A[j]);
          }
          if (any_unequal) {
            #pragma omp simd
            for (int j = 0; j < n_unique_trials; ++j) {
              gl_rt[static_cast<size_t>(j)] = gl_lo_B[static_cast<size_t>(j)]
                                            + hB[j] * (1.0 + xi);
            }
            model_dfun_raw(gl_rt.data(), colsp_nB, n_unique_trials,
                           msk, isok_nB_b.data(), lf_nB.data(), min_ll, model_ctx);
            model_pfun_raw(gl_rt.data(), colsp_B,  n_unique_trials,
                           msk, isok_B_b.data(),  lS_B.data(),  min_ll, model_ctx);
            for (int j = 0; j < n_unique_trials; ++j) {
              if (msk[j] && !ch_eq_vec[static_cast<size_t>(j)]) {
                GB[j] += wt * hB[j] * std::exp(lf_nB[j] + lS_B[j]);
              }
            }
          }
        }
        for (int j = 0; j < n_unique_trials; ++j) {
          if (msk[j] && ch_eq_vec[static_cast<size_t>(j)]) GB[j] = GA[j];
        }
      }
    }
  }

  static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace_tls, gsl_ctl.retry_limit);
  std::vector<double>& parA = scratch.parA;
  std::vector<double>& parB = scratch.parB;
  std::vector<double>& parnA = scratch.parnA;
  std::vector<double>& parnB = scratch.parnB;
  std::vector<double>& pars_2buf = scratch.pars_2buf;
  std::vector<double>& parN = scratch.tmpA1;
  std::vector<double>& pars_3buf = scratch.tmpA2;
  parA.resize(static_cast<size_t>(n_par));
  parB.resize(static_cast<size_t>(n_par));
  parnA.resize(static_cast<size_t>(n_par));
  parnB.resize(static_cast<size_t>(n_par));
  parN.resize(static_cast<size_t>(n_par));
  pars_2buf.resize(static_cast<size_t>(2 * n_par));
  int isok_2buf[2] = {1, 1};

  for (int j = 0; j < n_unique_trials; ++j) {
    const int idxA = shared.idxA[static_cast<size_t>(j)];
    const int idxB = shared.idxB[static_cast<size_t>(j)];
    const int idxnA = shared.idxnA[static_cast<size_t>(j)];
    const int idxnB = shared.idxnB[static_cast<size_t>(j)];

    const int rule_code = shared.rule_code[static_cast<size_t>(j)];
    // GNG (rule 6) is the full four-horse OR task whose "no" outcome is a
    // withheld response (rt = +Inf); it shares OR's subrace machinery here.
    const bool is_gng = (rule_code == 6);
    // Parameter-validity gate. OR/AND/XOR/ID and GNG need all four accumulators;
    // the analytic detector (rule 5) needs only valid A and B. Without this,
    // censored branches would evaluate survivors/densities with out-of-bounds
    // parameters.
    const bool pars_valid = (rule_code <= 4 || is_gng)
      ? (idxA >= 0 && idxB >= 0 && idxnA >= 0 && idxnB >= 0 &&
         ok_params[idxA] && ok_params[idxB] && ok_params[idxnA] && ok_params[idxnB])
      : (idxA >= 0 && idxB >= 0 && ok_params[idxA] && ok_params[idxB]);
    if (!pars_valid) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      if (count_capacity) ++lr_capacity_counters().invalid_trials;
      continue;
    }

    const double t = shared.rt_unique[static_cast<size_t>(j)];
    copy_par_row_colmajor(pars_cols, n_par, idxA, parA.data());
    copy_par_row_colmajor(pars_cols, n_par, idxB, parB.data());
    if (idxnA >= 0) copy_par_row_colmajor(pars_cols, n_par, idxnA, parnA.data());
    if (idxnB >= 0) copy_par_row_colmajor(pars_cols, n_par, idxnB, parnB.data());

    const int cell_idx = shared.cell_id[static_cast<size_t>(j)];

    // Capacity is a trial-level extension.  It is active only for redundant
    // target (AB) conditions and only when the particle's capacity effect is
    // non-degenerate.  Every other trial deliberately continues through the
    // legacy evaluator below, preserving its batched fast path and exact
    // baseline behaviour.
    if (capacity) {
      const double kappa = pars_cols[kappa_col][idxA];
      const double tau = pars_cols[tau_col][idxA];
      const bool shared_capacity =
        R_FINITE(kappa) && R_FINITE(tau) && tau >= 0.0 &&
        pars_cols[kappa_col][idxB] == kappa &&
        pars_cols[tau_col][idxB] == tau;
      if (!shared_capacity) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
        if (count_capacity) ++lr_capacity_counters().invalid_trials;
        continue;
      }
      const bool pair_active = shared.cond_code[static_cast<size_t>(j)] == 3 &&
        (kappa != 0.0 || tau != 0.0);
      if (pair_active) {
        if (count_capacity) {
          if (rule_code == 5) ++lr_capacity_counters().capacity_detection_trials;
          else ++lr_capacity_counters().capacity_choice_trials;
        }
        ll_unique[static_cast<size_t>(j)] = lrcap_trial_ll(
          rule_code,
          shared.cond_code[static_cast<size_t>(j)],
          shared.resp_code[static_cast<size_t>(j)],
          t,
          shared.LT_unique[static_cast<size_t>(j)],
          shared.UT_unique[static_cast<size_t>(j)],
          shared.LC_unique[static_cast<size_t>(j)],
          shared.UC_unique[static_cast<size_t>(j)],
          kappa, tau,
          pars_cols, n_par, idxA, idxB,
          parnA.data(), parnB.data(), parN.data(),
          rule_code <= 4 || is_gng, false,
          min_ll, kMinSurv,
          pdf1, cdf1, model_ctx,
          lrcap_scan_n, lrcap_fine_n,
          count_capacity, cell_idx, &shared.cell_log_den_cap,
          gsl_ctl, w, pars_3buf);
        continue;
      }
      if (count_capacity) ++lr_capacity_counters().ordinary_trials;
    }

    // Channel-parameter equality (A vs B), shared by the truncation
    // normaliser and the censoring branches below.
    const bool ch_eq_j = use_gl_pass ? ch_eq_vec[static_cast<size_t>(j)]
      : ((rule_code <= 4 || is_gng) &&
         row_equal_colmajor(pars_cols, n_par, idxA, idxB) &&
         row_equal_colmajor(pars_cols, n_par, idxnA, idxnB));

    // Sub-race win probability G_no at time tt for one channel. The GL
    // pre-pass batch-computes these at the aux slots (LT/UT/censor bounds);
    // slot values are used when present and matching tt, otherwise fall back
    // to the scalar GSL helper (non-raw path, or an uncovered combination).
    auto lr4_G_no_at = [&](double tt, int aux_slot, bool useB, bool& okf) -> double {
      okf = true;
      if (R_FINITE(tt) && use_gl_pass && aux_slot >= 0) {
        const size_t o = static_cast<size_t>(aux_slot) * n_unique_trials + j;
        const double g = (useB && !ch_eq_j) ? scratch.aux_GB[o] : scratch.aux_GA[o];
        if (!emc2_isnan(g) && !emc2_isnan(scratch.aux_t[o]) &&
            std::fabs(scratch.aux_t[o] - tt) < 1e-12) {
          return std::min(1.0, std::max(0.0, g));
        }
      }
      // tt = +Inf (GNG withheld / eventual N): P(nontarget eventually wins the
      // channel) via the QAGIU-capable helper; local_race_helper reports 0 for
      // an infinite horizon, so route infinite tt through lr_channel_yes_prob.
      const double g = R_FINITE(tt)
        ? local_race_helper(tt,
                            useB ? parnB.data() : parnA.data(),
                            useB ? parB.data()  : parA.data(),
                            n_par, model_ctx, pdf1, cdf1, gsl_ctl, w,
                            pars_2buf.data(), isok_2buf)
        : lr_channel_yes_prob(tt,
                            useB ? parnB.data() : parnA.data(),
                            useB ? parB.data()  : parA.data(),
                            n_par, model_ctx, pdf1, cdf1, gsl_ctl, w,
                            pars_2buf.data(), isok_2buf);
      okf = R_FINITE(g);
      return okf ? std::min(1.0, std::max(0.0, g)) : 0.0;
    };

    // Full channel state at tt (survivors are cheap scalar calls; G_no from
    // the batch/fallback above; G_yes by complement).
    auto lr4_state_at = [&](double tt, int aux_slot, bool useB, bool& okf) -> LrChannelState {
      LrChannelState s{0.0, 0.0, 1.0};
      okf = true;
      if (!(tt > 0.0)) return s;
      const double* pT = useB ? parB.data()  : parA.data();
      const double* pN = useB ? parnB.data() : parnA.data();
      const double S_dec = safe_surv_at_race_scalar(pT, tt, kMinSurv, cdf1, model_ctx)
                         * safe_surv_at_race_scalar(pN, tt, kMinSurv, cdf1, model_ctx);
      const double G_no = lr4_G_no_at(tt, aux_slot, useB, okf);
      if (!okf) return s;
      s.S_dec = S_dec;
      s.G_no  = G_no;
      s.G_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec - G_no));
      return s;
    };

    // N(t) = P(no overt rule response by t) for OR/AND/XOR/ID. XOR/ID need
    // only the channel survivors; OR/AND additionally need G_no.
    auto lr4_N_at = [&](double tt, int aux_slot, bool& okf) -> double {
      okf = true;
      if (!(tt > 0.0)) return 1.0;
      if (rule_code == 3 || rule_code == 4) {
        const double SdA =
          safe_surv_at_race_scalar(parA.data(),  tt, kMinSurv, cdf1, model_ctx) *
          safe_surv_at_race_scalar(parnA.data(), tt, kMinSurv, cdf1, model_ctx);
        const double SdB = ch_eq_j ? SdA
          : safe_surv_at_race_scalar(parB.data(),  tt, kMinSurv, cdf1, model_ctx) *
            safe_surv_at_race_scalar(parnB.data(), tt, kMinSurv, cdf1, model_ctx);
        return 1.0 - (1.0 - SdA) * (1.0 - SdB);
      }
      LrChannelState Ast = lr4_state_at(tt, aux_slot, false, okf);
      if (!okf) return 1.0;
      LrChannelState Bst = Ast;
      if (!ch_eq_j) {
        Bst = lr4_state_at(tt, aux_slot, true, okf);
        if (!okf) return 1.0;
      }
      if (is_gng) {
        // GNG: "no overt go response by t" = neither channel has said yes yet.
        // A completed channel-no is a withheld outcome here, not an overt
        // response, so (unlike OR) it stays inside N. tt may be +Inf (the
        // eventual withheld mass), which lr4_state_at handles via lr4_G_no_at.
        const double QA = std::min(1.0, std::max(0.0, Ast.G_no + Ast.S_dec));
        const double QB = std::min(1.0, std::max(0.0, Bst.G_no + Bst.S_dec));
        return QA * QB;
      }
      return lr_no_response_prob(rule_code, Ast, Bst);
    };

    // --- Truncation normaliser -------------------------------------------
    // Z = P(overt response RT in [LT, UT]) under the rule's outcome logic:
    // Z = N(LT) - N(UT), where the UT = Inf case is DEFINED as Z = N(LT)
    const double LTj_tr = shared.LT_unique[static_cast<size_t>(j)];
    const double UTj_tr = shared.UT_unique[static_cast<size_t>(j)];
    const bool has_trunc = (LTj_tr != 0.0 || R_FINITE(UTj_tr));
    double log_Z_j = 0.0;
    if (has_trunc) {
      if (!ISNAN(shared.cell_log_z[cell_idx])) {
        log_Z_j = shared.cell_log_z[cell_idx];
        if (log_Z_j == R_NegInf) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      } else {
        const int cond_j = shared.cond_code[static_cast<size_t>(j)];
        bool z_ok = true;
        double N_LT = 1.0;
        double N_UT = 0.0;
        if (rule_code == 5) {
          if (LTj_tr > 0.0) {
            N_LT = lr_detection_no_response_prob(rule_code, cond_j, LTj_tr,
                                                 parA.data(), parB.data(), parN.data(),
                                                 n_par, kMinSurv, pdf1, cdf1, model_ctx,
                                                 gsl_ctl, w, pars_3buf, z_ok);
          }
          if (z_ok && R_FINITE(UTj_tr)) {
            N_UT = lr_detection_no_response_prob(rule_code, cond_j, UTj_tr,
                                                  parA.data(), parB.data(), parN.data(),
                                                  n_par, kMinSurv, pdf1, cdf1, model_ctx,
                                                  gsl_ctl, w, pars_3buf, z_ok);
          }
        } else {
          if (LTj_tr > 0.0) N_LT = lr4_N_at(LTj_tr, LR_AUX_LT, z_ok);
          if (z_ok && R_FINITE(UTj_tr)) N_UT = lr4_N_at(UTj_tr, LR_AUX_UT, z_ok);
        }
        const double Z = std::min(1.0, N_LT - N_UT);
        const bool valid_z = z_ok && R_FINITE(Z) && (Z > 1e-12);
        log_Z_j = valid_z ? std::log(Z) : R_NegInf;
        shared.cell_log_z[cell_idx] = log_Z_j;
        if (!valid_z) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
    }

    if (rule_code == 5) {
      ll_unique[static_cast<size_t>(j)] = logicalrules_detection_trial_ll(
        rule_code,
        shared.cond_code[static_cast<size_t>(j)],
        shared.resp_code[static_cast<size_t>(j)],
        t,
        shared.LT_unique[static_cast<size_t>(j)],
        shared.UT_unique[static_cast<size_t>(j)],
        shared.LC_unique[static_cast<size_t>(j)],
        shared.UC_unique[static_cast<size_t>(j)],
        idxA,
        idxB,
        shared.idxNogo[static_cast<size_t>(j)],
        n_par,
        min_ll,
        kMinSurv,
        log_Z_j,
        use_raw_local,
        logf_all,
        logS_all,
        parA.data(),
        parB.data(),
        parN.data(),
        pdf1,
        cdf1,
        model_ctx,
        gsl_ctl,
        w,
        pars_3buf);
      continue;
    }

    if (!(t > 0.0) || !R_FINITE(t)) {
      double p_j = 0.0;
      const int rule_code_nc = shared.rule_code[static_cast<size_t>(j)];
      if (!R_FINITE(t) && t > 0.0) {
        // Upper censor (+Inf): probability no overt response fires by UC.
        // Semantics depend on rule because each rule has a different "first
        // response" trigger (see lr_no_response_prob):
        //   XOR/ID : response fires when both channels resolve  → omit iff either channel undecided
        //   OR     : yes = first yes; no = second no           → omit = P(no yes by UC) − P(both no by UC)
        //   AND    : yes = second yes; no = first no           → omit = P(no no by UC) − P(both yes by UC)
        // With a finite upper truncation bound (and a proper upper tail) the
        // response is known to land in (UC, UT], so N(UT) is subtracted; the
        // trial is then normalised by Z below.
        // GNG: rt = +Inf is the withheld response — no overt go by UC, where
        // N_GNG(t) = Q_A(t) Q_B(t) (a completed channel-no is withheld, not an
        // overt response). UC may be +Inf (the withheld-to-infinity mass), so
        // GNG does not require a finite UC (unlike OR/AND/XOR/ID, whose rt=+Inf
        // omission is only defined against a finite censoring bound).
        const double UCj = shared.UC_unique[static_cast<size_t>(j)];
        bool ok_uc = true;
        if (is_gng) {
          p_j = lr4_N_at(UCj, R_FINITE(UCj) ? LR_AUX_UC : -1, ok_uc);
        } else {
          if (!R_FINITE(UCj)) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
          p_j = lr4_N_at(UCj, LR_AUX_UC, ok_uc);
        }
        if (ok_uc && R_FINITE(UTj_tr) && !model_ctx->defective_upper_tail) {
          bool ok_ut = true;
          const double N_UT_uc = lr4_N_at(UTj_tr, LR_AUX_UT, ok_ut);
          ok_uc = ok_uc && ok_ut;
          p_j = std::max(0.0, p_j - N_UT_uc);
        }
        if (!ok_uc) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
      } else if (!R_FINITE(t) && t < 0.0) {
        // Lower censor (-Inf): an overt rule response occurred in [LT, LC],
        // with identity resp_code when recorded (0 = unknown). The response is
        // NOT the first accumulator event — each rule triggers on its own
        // channel-outcome pattern — so the mass is a difference of rule
        // response CDFs built from the channel states at the window ends
        // (previously this used S_joint(LT) - S_joint(LC), the first-event
        // mass, which overcounts whenever one channel resolves without an
        // overt response having occurred yet).
        const double lo = std::max(0.0, shared.LT_unique[static_cast<size_t>(j)]);
        const double hi = std::max(lo,  shared.LC_unique[static_cast<size_t>(j)]);
        const int resp_lc = shared.resp_code[static_cast<size_t>(j)];
        if (is_gng) {
          // Lower censor: an overt go response fired in [lo, hi]. The go-CDF is
          // 1 - N_GNG, so the window mass is N_GNG(lo) - N_GNG(hi).
          bool ok_lc = true;
          const double N_lo = lr4_N_at(lo, LR_AUX_CLO, ok_lc);
          const double N_hi = ok_lc ? lr4_N_at(hi, LR_AUX_CHI, ok_lc) : 1.0;
          if (!ok_lc) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
          p_j = std::max(0.0, N_lo - N_hi);
          if (!(p_j > 0.0) || !R_FINITE(p_j)) {
            ll_unique[static_cast<size_t>(j)] = min_ll;
          } else {
            const double ll = std::log(p_j) - log_Z_j;
            ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
          }
          continue;
        }
        bool states_ok = true;
        LrChannelState A_lo = lr4_state_at(lo, LR_AUX_CLO, false, states_ok);
        LrChannelState A_hi, B_lo, B_hi;
        if (states_ok) A_hi = lr4_state_at(hi, LR_AUX_CHI, false, states_ok);
        if (states_ok) {
          if (ch_eq_j) {
            B_lo = A_lo;
            B_hi = A_hi;
          } else {
            B_lo = lr4_state_at(lo, LR_AUX_CLO, true, states_ok);
            if (states_ok) B_hi = lr4_state_at(hi, LR_AUX_CHI, true, states_ok);
          }
        }
        if (!states_ok) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
        const double cdf_hi = lr_response_cdf(rule_code_nc, resp_lc, A_hi, B_hi);
        const double cdf_lo = lr_response_cdf(rule_code_nc, resp_lc, A_lo, B_lo);
        if (ISNAN(cdf_hi) || ISNAN(cdf_lo)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
        p_j = std::max(0.0, cdf_hi - cdf_lo);
      }
      if (!(p_j > 0.0) || !R_FINITE(p_j)) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
      } else {
        const double ll = std::log(p_j) - log_Z_j;
        ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
      }
      continue;
    }

    double fA = 0.0, fB = 0.0, fnA = 0.0, fnB = 0.0;
    double FA = 0.0, FB = 0.0, FnA = 0.0, FnB = 0.0;
    if (!eval_pdf_cdf_race_scalar(use_raw_local, idxA, t, parA.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fA, FA) ||
        !eval_pdf_cdf_race_scalar(use_raw_local, idxB, t, parB.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fB, FB) ||
        !eval_pdf_cdf_race_scalar(use_raw_local, idxnA, t, parnA.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fnA, FnA) ||
        !eval_pdf_cdf_race_scalar(use_raw_local, idxnB, t, parnB.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fnB, FnB)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    const double one_m_FA = std::max(kMinSurv, 1.0 - FA);
    const double one_m_FB = std::max(kMinSurv, 1.0 - FB);
    const double one_m_FnA = std::max(kMinSurv, 1.0 - FnA);
    const double one_m_FnB = std::max(kMinSurv, 1.0 - FnB);

    const double gA_yes = fA * one_m_FnA;
    const double gB_yes = fB * one_m_FnB;
    const double gA_no = fnA * one_m_FA;
    const double gB_no = fnB * one_m_FB;

    const double S_dec_A = one_m_FA * one_m_FnA;
    const double S_dec_B = one_m_FB * one_m_FnB;

    const int rule_code_lr = shared.rule_code[static_cast<size_t>(j)];
    const int resp_code = shared.resp_code[static_cast<size_t>(j)];
    // GNG shares OR's finite-RT go-response density; only its "no" outcome
    // differs (withheld, handled above), so a finite RT is always a "yes".
    const bool is_or_rule = (rule_code_lr == 1 || rule_code_lr == 6);
    const bool is_and_rule = (rule_code_lr == 2);
    const bool is_xor_rule = (rule_code_lr == 3);
    const bool is_id_rule = (rule_code_lr == 4);
    const bool channels_equal = use_gl_pass ? ch_eq_vec[static_cast<size_t>(j)]
        : (row_equal_colmajor(pars_cols, n_par, idxA, idxB) &&
           row_equal_colmajor(pars_cols, n_par, idxnA, idxnB));

    double p_j = 0.0;
    if (is_or_rule || is_and_rule) {
      double GA_no, GB_no;
      if (use_gl_pass) {
        GA_no = std::max(0.0, std::min(GA_no_gl[static_cast<size_t>(j)], 1.0));
        GB_no = channels_equal ? GA_no
                               : std::max(0.0, std::min(GB_no_gl[static_cast<size_t>(j)], 1.0));
      } else {
        GA_no = local_race_helper(t, parnA.data(), parA.data(), n_par, model_ctx,
                                  pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        GB_no = GA_no;
        if (!channels_equal)
          GB_no = local_race_helper(t, parnB.data(), parB.data(), n_par, model_ctx,
                                    pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        if (!R_FINITE(GA_no) || !R_FINITE(GB_no)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
      const double s_GA_yes = std::min(1.0, std::max(kMinSurv, GA_no + S_dec_A));
      const double s_GB_yes = std::min(1.0, std::max(kMinSurv, GB_no + S_dec_B));
      if (is_or_rule) {
        if (resp_code == 1) p_j = gA_yes * s_GB_yes + gB_yes * s_GA_yes; // yes
        else if (resp_code == 2 && rule_code_lr != 6) p_j = gA_no * GB_no + gB_no * GA_no; // no
      } else { // AND
        const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
        const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
        const double s_GA_no = std::min(1.0, std::max(kMinSurv, GA_yes + S_dec_A));
        const double s_GB_no = std::min(1.0, std::max(kMinSurv, GB_yes + S_dec_B));
        if (resp_code == 1) p_j = gA_yes * GB_yes + gB_yes * GA_yes;      // yes
        else if (resp_code == 2) p_j = gA_no * s_GB_no + gB_no * s_GA_no; // no
      }
    } else if (is_xor_rule || is_id_rule) {
      double GA_yes, GB_yes;
      if (use_gl_pass) {
        // Derive GA_yes from GA_no_gl: GA_yes + GA_no + S_dec_A = 1.
        const double GA_no = std::max(0.0, std::min(GA_no_gl[static_cast<size_t>(j)], 1.0));
        const double GB_no = channels_equal ? GA_no
                                            : std::max(0.0, std::min(GB_no_gl[static_cast<size_t>(j)], 1.0));
        GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
        GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
      } else {
        GA_yes = local_race_helper(t, parA.data(), parnA.data(), n_par, model_ctx,
                                   pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        GB_yes = GA_yes;
        if (!channels_equal)
          GB_yes = local_race_helper(t, parB.data(), parnB.data(), n_par, model_ctx,
                                     pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        if (!R_FINITE(GA_yes) || !R_FINITE(GB_yes)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
      const double GA_no = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_yes));
      const double GB_no = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_yes));
      const double dAB = gA_yes * GB_yes + gB_yes * GA_yes;
      const double dAN = gA_yes * GB_no + gB_no * GA_yes;
      const double dNB = gA_no * GB_yes + gB_yes * GA_no;
      const double dNN = gA_no * GB_no + gB_no * GA_no;

      if (is_id_rule) {
        if (resp_code == 6) p_j = dAB;
        else if (resp_code == 4) p_j = dAN;
        else if (resp_code == 5) p_j = dNB;
        else if (resp_code == 3) p_j = dNN;
      } else { // XOR
        if (resp_code == 1) p_j = dAN + dNB; // yes
        else if (resp_code == 2) p_j = dAB + dNN; // no
      }
    } else {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    if (!(p_j > 0.0) || !R_FINITE(p_j)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
    } else {
      const double ll = std::log(p_j) - log_Z_j;
      ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
    }
  }

  // Contaminant mixture: pC for intrinsic omissions (+Inf RT), pGuess for the
  // uniform outlier on observed RTs.  Both read the first row of each trial.
  // See src/contaminant_mixture.h -- the arithmetic lives there so this and the
  // seven other application sites cannot drift apart.
  bool use_pC = (pc_col >= 0);
  if (use_pC) {
    bool all_zero = true;
    for (int j = 0; j < n_unique_trials; ++j) {
      if (pars_cols[pc_col][j * n_acc] != 0.0) { all_zero = false; break; }
    }
    if (all_zero) use_pC = false;
  }

  const GuessKernel gk = (guess != nullptr) ? *guess : GuessKernel();
  bool use_pG = (pg_col >= 0) && gk.active();
  if (use_pG) {
    bool all_zero = true;
    for (int j = 0; j < n_unique_trials; ++j) {
      if (pars_cols[pg_col][j * n_acc] != 0.0) { all_zero = false; break; }
    }
    if (all_zero) use_pG = false;
  }

  if (use_pC || use_pG) {
    for (int j = 0; j < n_unique_trials; ++j) {
      const size_t sj = static_cast<size_t>(j);
      const double pC = use_pC ? pars_cols[pc_col][sj * n_acc] : 0.0;
      const double pG = use_pG ? pars_cols[pg_col][sj * n_acc] : 0.0;
      // resp_code 0 is a missing response; every other code names a response.
      ll_unique[sj] = mix_contaminants_rt(ll_unique[sj], pC, pG, gk,
                                          shared.rt_unique[sj],
                                          shared.resp_code[sj] != 0);
    }
  }

  double sum_ll = 0.0;
  if (expand.length() > 0) {
    const int n_out = expand.length();
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_logicalrules: trial_ll_out size mismatch (expand path).");
    }
    for (int i = 0; i < n_out; ++i) {
      double val = ll_unique[static_cast<size_t>(expand[i] - 1)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = val;
      sum_ll += val;
    }
  } else {
    const int n_out = n_unique_trials;
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_logicalrules: trial_ll_out size mismatch (compressed path).");
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      double val = ll_unique[static_cast<size_t>(j)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = val;
      sum_ll += val;
    }
  }
  return sum_ll;
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
  // Fast default controls for MCMC throughput, with an automatic stricter retry.
  // Retry settings match the old behaviour (rel_tol=1e-7, limit=1000).
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
  // Several branches below stage one accumulator row (or one pointer per
  // accumulator) in fixed stack arrays of 1024 entries; guard the sizes here so
  // a future wide model fails loudly instead of overflowing the stack.
  if (n_par > 1024 || n_lR > 1024) {
    Rcpp::stop("c_log_likelihood_race: at most 1024 parameter columns and 1024 accumulators are supported.");
  }
  if (n_par * n_lR > 32768) {
    Rcpp::stop("c_log_likelihood_race: n_par * n_lR exceeds 32768, which would overflow the safe stack size.");
  }
  // Column-pointer view of the materialized pars matrix for the raw kernels
  // (this path keeps the matrix: the RACE NA-fill above writes into it).
  const double* cols_view[1024] = {nullptr};
  for (int c = 0; c < n_par; ++c) {
    cols_view[c] = pars_cm_ptr + static_cast<size_t>(c) * n_trials;
  }
  double pars_rowmajor_buffer[32768];
  int isok_int_buffer[1024] = {0};
  double logS_k_buffer[1024];
  for (int i = 0; i < n_lR; ++i) {
    logS_k_buffer[i] = R_NegInf;
  }

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
          double par_buf_inf[64];
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
          double par_buf_inf[64];
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
        double par_buf_inf[64];
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
      double par_buf[64];
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
    double par_buf[64];
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

      if (any_trunc && uniform_LT_ok && uniform_UT_ok && !has_RACE_col) {
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
            double par_buf[64];
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
          double par_buf[64];
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

// Build the data-fixed correlated BAwL state once per likelihood call: lM
// role mapping, RACE masks, truncation windows, winner rows, unique-trial
// expansion, and direct ParamTable column pointers in keep_names order.
// The fast OO mapper intentionally stops at the design-matrix/natural-scale
// mapping and does not invoke the model's R Ttransform.  The role mapping
// resolved here mirrors BAwLcorr's role conversion for the likelihood path;
// R Ttransform remains the corresponding conversion for mapped parameter
// displays and simulation.
static BAwLCorrSharedState build_bawl_corr_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names) {
  BAwLCorrSharedState s;
  if (n_lR <= 0 || n_trials < 0 || n_trials % n_lR != 0) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: invalid race dimensions.");
  }
  s.n_trials = n_trials;
  s.n_lR = n_lR;
  s.n_unique = n_trials / n_lR;
  s.n_out = (expand.length() > 0) ? expand.length() : s.n_unique;
  s.n_par = keep_names.size();
  s.rho_col = ctx.bawl_rho_index;
  if (s.rho_col < 0 || s.rho_col >= s.n_par) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: rho column is outside pars.");
  }

  // lM role mapping.
  bool lM_is_logical = false;
  Rcpp::LogicalVector lM_logical;
  Rcpp::IntegerVector lM_factor;
  int lM_true_code = -1;
  if (!dadm.containsElementNamed("lM")) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: BAwLcorr requires lM from matchfun.");
  }
  SEXP lM_sexp = dadm["lM"];
  if (TYPEOF(lM_sexp) == LGLSXP && Rf_length(lM_sexp) == n_trials) {
    lM_logical = Rcpp::LogicalVector(lM_sexp);
    lM_is_logical = true;
  } else if (TYPEOF(lM_sexp) == INTSXP && Rf_length(lM_sexp) == n_trials) {
    lM_factor = Rcpp::IntegerVector(lM_sexp);
    Rcpp::CharacterVector levels = lM_factor.attr("levels");
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == "TRUE") {
        lM_true_code = i + 1;
        break;
      }
    }
    if (lM_true_code < 0 && levels.size() == 0) lM_true_code = 1;
  }
  if (!lM_is_logical && lM_true_code < 0) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: lM must be a logical or TRUE/FALSE factor.");
  }
  s.role_correct.assign(static_cast<size_t>(n_trials), 0);
  for (int r = 0; r < n_trials; ++r) {
    const bool missing = lM_is_logical
      ? (lM_logical[r] == NA_LOGICAL)
      : (lM_factor[r] == NA_INTEGER);
    if (missing) {
      Rcpp::stop("c_log_likelihood_bawl_correlated: lM contains missing role values.");
    }
    const bool correct = lM_is_logical
      ? static_cast<bool>(lM_logical[r])
      : (lM_factor[r] == lM_true_code);
    s.role_correct[static_cast<size_t>(r)] = correct ? 1 : 0;
  }

  // RACE masks (attribute-cached variable-accumulator designs only).
  s.has_RACE = dadm.containsElementNamed("RACE");
  if (s.has_RACE && dadm.hasAttribute("RACE_nacc_by_row") &&
      dadm.hasAttribute("RACE_mask")) {
    s.race_nacc = dadm.attr("RACE_nacc_by_row");
    s.race_mask = dadm.attr("RACE_mask");
    if (s.race_nacc.size() != n_trials || s.race_mask.size() != n_trials)
      s.has_RACE = false;
  } else {
    s.has_RACE = false;
  }

  // Truncation windows.  emc2_all_finite_trials guarantees LT = 0/UT = Inf,
  // so these are all zero whenever the caller reports all-finite data.
  s.rt = dadm["rt"];
  s.LT = get_col_with_default(dadm, "LT", 0.0);
  s.UT = get_col_with_default(dadm, "UT", R_PosInf);
  s.LC = get_col_with_default(dadm, "LC", R_PosInf);
  s.UC = get_col_with_default(dadm, "UC", 0.0);
  s.has_trunc_trial.assign(static_cast<size_t>(s.n_unique), 0);
  for (int j = 0; j < s.n_unique; ++j) {
    const int start = j * n_lR;
    if (s.LT[start] != 0.0 || s.UT[start] != R_PosInf) {
      s.has_trunc_trial[static_cast<size_t>(j)] = 1;
      s.has_truncation = true;
    }
  }

  // Winner rows and unique-trial expansion map.
  s.winner_row.assign(static_cast<size_t>(s.n_unique), -1);
  for (int r = 0; r < n_trials; ++r) {
    if (winner[r]) s.winner_row[static_cast<size_t>(r / n_lR)] = r;
  }
  s.j_to_i.assign(static_cast<size_t>(s.n_unique), -1);
  if (expand.length() > 0) {
    for (int i = s.n_out - 1; i >= 0; --i)
      s.j_to_i[static_cast<size_t>(expand[i] - 1)] = i;
  } else {
    for (int j = 0; j < s.n_unique; ++j) s.j_to_i[static_cast<size_t>(j)] = j;
  }

  // Direct column pointers; base storage is refilled in place per particle.
  // base_index_for throws for unknown names exactly like materialisation.
  s.cols.assign(static_cast<size_t>(std::max(s.n_par, 16)), nullptr);
  for (int j = 0; j < s.n_par; ++j) {
    const std::string nm = Rcpp::as<std::string>(keep_names[j]);
    s.cols[static_cast<size_t>(j)] = &table.base(0, table.base_index_for(nm));
    if (nm == "pContaminant") s.pc_col = j;
    else if (nm == "pGuess") s.pg_col = j;
  }
  s.guess = resolve_guess_kernel(dadm);
  s.guess.pg_col = s.pg_col;

  s.effective_rho.assign(static_cast<size_t>(n_trials), R_NegInf);
  s.layout.assign(static_cast<size_t>(s.n_unique), BAwLCorrTrialLayout());
  s.valid = true;
  return s;
}

// Canonical per-particle trial classification (plan design rule C).  Counts
// only isok && RACE-active rows; an exact-zero loading is an independent
// singleton, not a loaded dimension.  This classifier is the sole authority
// for loaded-row discovery and the positivity dimension.
static void bawl_corr_classify_particle(BAwLCorrSharedState& s,
                                        const Rcpp::LogicalVector& isok,
                                        bool no_clock_eligible) {
  constexpr double kLoadingEps = 1e-14;
  const int n_lR = s.n_lR;
  for (int j = 0; j < s.n_unique; ++j) {
    BAwLCorrTrialLayout& L = s.layout[static_cast<size_t>(j)];
    L = BAwLCorrTrialLayout();
    const int start = j * n_lR;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : n_lR;
    L.winner_row = s.winner_row[static_cast<size_t>(j)];
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (s.has_RACE && !s.race_mask[row]) continue;
      ++L.n_active;
      const double rho = s.effective_rho[static_cast<size_t>(row)];
      if (!isok[row] || !R_FINITE(rho) || std::fabs(rho) > 1.0) {
        L.any_bad_row = true;
        continue;
      }
      if (std::fabs(rho) > kLoadingEps) {
        if (L.n_loaded < 2) L.loaded_row[L.n_loaded] = row;
        ++L.n_loaded;
        if (row == L.winner_row) L.winner_loaded = true;
      }
    }
    if (L.any_bad_row) {
      L.route = BAwLCorrRoute::invalid;
    } else if (L.n_loaded <= 1) {
      L.route = BAwLCorrRoute::ordinary;
    } else {
      // A two-loaded no-clock trial is the exact component route.  Trials
      // with a larger factor dimension retain the fused GH fallback; active
      // clocks retain the generic race semantics even when only two rows are
      // loaded.
      L.route = no_clock_eligible
        ? ((L.n_loaded == 2) ? BAwLCorrRoute::exact_pair
                             : BAwLCorrRoute::gh_no_clock)
        : BAwLCorrRoute::gh_generic_clock;
    }
  }
}

struct BAwLCorrPairData {
  int row1 = -1;
  int row2 = -1;
  double mu1 = 0.0;
  double mu2 = 0.0;
  double sd1 = 1.0;
  double sd2 = 1.0;
  double rho = 0.0;
  double normalizer = 1.0;
  bool positive = false;
};

static inline double bawl_corr_row_k(const double* const* cols, int row,
                                     const ContextForRaceModels* ctx) {
  return ctx->bawl_k_fixed_zero ? 0.0 : cols[emc2col::bawl::k][row];
}

static inline BAwLTimeGeometry bawl_corr_row_geometry(
    double t, const double* const* cols, int row,
    const ContextForRaceModels* ctx) {
  const double A = cols[emc2col::bawl::A][row];
  const double B = cols[emc2col::bawl::B][row];
  return bawl_time_geometry(t, cols[emc2col::bawl::t0][row], A, B + A,
                            bawl_corr_row_k(cols, row, ctx));
}

static inline bool bawl_corr_make_pair_data(
    const BAwLCorrTrialLayout& layout, const BAwLCorrSharedState& s,
    const double* const* cols, const ContextForRaceModels* ctx,
    BAwLCorrPairData& out) {
  if (layout.n_loaded != 2 || layout.loaded_row[0] < 0 ||
      layout.loaded_row[1] < 0) return false;
  out.row1 = layout.loaded_row[0];
  out.row2 = layout.loaded_row[1];
  out.mu1 = cols[emc2col::bawl::v][out.row1];
  out.mu2 = cols[emc2col::bawl::v][out.row2];
  out.sd1 = cols[emc2col::bawl::sv][out.row1];
  out.sd2 = cols[emc2col::bawl::sv][out.row2];
  const double r1 = s.effective_rho[static_cast<size_t>(out.row1)];
  const double r2 = s.effective_rho[static_cast<size_t>(out.row2)];
  if (!(out.sd1 > 0.0) || !(out.sd2 > 0.0) || !R_FINITE(out.mu1) ||
      !R_FINITE(out.mu2) || !R_FINITE(r1) || !R_FINITE(r2) ||
      std::fabs(r1) > 1.0 || std::fabs(r2) > 1.0) return false;
  out.rho = (r1 * r2 < 0.0 ? -1.0 : 1.0) *
    std::sqrt(std::fabs(r1 * r2));
  out.positive = ctx->use_posdrift;
  const double key[5] = {out.mu1, out.sd1, out.mu2, out.sd2, out.rho};
  for (const auto& e : s.qab_cache) {
    if (e.key[0] == key[0] && e.key[1] == key[1] && e.key[2] == key[2] &&
        e.key[3] == key[3] && e.key[4] == key[4]) {
      out.normalizer = e.value;
      return R_FINITE(out.normalizer) && out.normalizer > 0.0;
    }
  }
  out.normalizer = bawl_corr_pair_positive_normalizer(
      out.mu1, out.sd1, out.mu2, out.sd2, out.rho, out.positive);
  if (!R_FINITE(out.normalizer) || !(out.normalizer > 0.0)) return false;
  if (s.qab_cache.size() < 64) {
    s.qab_cache.push_back({{key[0], key[1], key[2], key[3], key[4]},
                           out.normalizer});
  }
  return R_FINITE(out.normalizer) && out.normalizer > 0.0;
}

static inline double bawl_corr_single_log_survival(
    double t, int row, const double* const* cols,
    const ContextForRaceModels* ctx) {
  const double t0 = cols[emc2col::bawl::t0][row];
  const double tau = t - t0;
  if (ISNAN(t) || ISNAN(t0) || !(t > 0.0) || !(tau > 0.0)) return 0.0;
  const double A = cols[emc2col::bawl::A][row];
  const double B = cols[emc2col::bawl::B][row];
  const double v = cols[emc2col::bawl::v][row];
  const double sv = cols[emc2col::bawl::sv][row];
  const double k = bawl_corr_row_k(cols, row, ctx);
  const double log_cdf = log_ba_cdf(
      tau, A, B + A, v, sv, k, ctx->use_posdrift,
      ctx->bawl_k_fixed_zero ? LBA_DENOM_FLOOR : BAWL_DENOM_FLOOR);
  if (log_cdf == R_NegInf) return 0.0;
  if (log_cdf >= 0.0) return R_NegInf;
  return log1m_exp(log_cdf);
}

static inline double bawl_corr_single_log_density(
    double t, int row, const double* const* cols,
    const ContextForRaceModels* ctx) {
  if (!(t > 0.0) || !R_FINITE(t)) return R_NegInf;
  const double t0 = cols[emc2col::bawl::t0][row];
  const double tau = t - t0;
  if (!(tau > 0.0) || !R_FINITE(tau)) return R_NegInf;
  const double A = cols[emc2col::bawl::A][row];
  const double B = cols[emc2col::bawl::B][row];
  const double v = cols[emc2col::bawl::v][row];
  const double sv = cols[emc2col::bawl::sv][row];
  const double k = bawl_corr_row_k(cols, row, ctx);
  return log_ba_pdf(tau, A, B + A, v, sv, k, ctx->use_posdrift,
                    ctx->bawl_k_fixed_zero ? LBA_DENOM_FLOOR : BAWL_DENOM_FLOOR);
}

static inline double bawl_corr_log_pair_survival(
    const BAwLCorrPairData& pair, double t, const double* const* cols,
    const ContextForRaceModels* ctx, bool numeric,
    BAwLCorrMomentStatus* status_out = nullptr) {
  const BAwLTimeGeometry g1 = bawl_corr_row_geometry(t, cols, pair.row1, ctx);
  const BAwLTimeGeometry g2 = bawl_corr_row_geometry(t, cols, pair.row2, ctx);
  const bool use_numeric = numeric || std::fabs(pair.rho) >= 1.0 - 1e-10;
  const BAwLCorrPairResult result = use_numeric
    ? bawl_corr_pair_survival_numeric(g1, g2, pair.mu1, pair.sd1, pair.mu2,
                                      pair.sd2, pair.rho, pair.positive,
                                      pair.normalizer)
    : bawl_corr_pair_survival_exact(g1, g2, pair.mu1, pair.sd1, pair.mu2,
                                    pair.sd2, pair.rho, pair.positive,
                                    pair.normalizer);
  if (status_out != nullptr) *status_out = result.status;
  if (!(result.value > 0.0) || !R_FINITE(result.value)) return R_NegInf;
  return std::log(result.value);
}

static inline double bawl_corr_log_pair_cause(
    const BAwLCorrPairData& pair, int winner_row, double t,
    const double* const* cols, const ContextForRaceModels* ctx, bool numeric,
    BAwLCorrMomentStatus* status_out = nullptr) {
  const bool first = winner_row == pair.row1;
  const int rw = first ? pair.row1 : pair.row2;
  const int rl = first ? pair.row2 : pair.row1;
  const BAwLTimeGeometry gw = bawl_corr_row_geometry(t, cols, rw, ctx);
  const BAwLTimeGeometry gl = bawl_corr_row_geometry(t, cols, rl, ctx);
  const bool use_numeric = numeric || std::fabs(pair.rho) >= 1.0 - 1e-10;
  const BAwLCorrPairResult result = use_numeric
    ? bawl_corr_pair_cause_numeric(
        gw, gl, cols[emc2col::bawl::v][rw], cols[emc2col::bawl::sv][rw],
        cols[emc2col::bawl::v][rl], cols[emc2col::bawl::sv][rl], pair.rho,
        pair.positive, pair.normalizer)
    : bawl_corr_pair_cause_exact(
        gw, gl, cols[emc2col::bawl::v][rw], cols[emc2col::bawl::sv][rw],
        cols[emc2col::bawl::v][rl], cols[emc2col::bawl::sv][rl], pair.rho,
        pair.positive, pair.normalizer);
  if (status_out != nullptr) *status_out = result.status;
  if (!(result.value > 0.0) || !R_FINITE(result.value)) return R_NegInf;
  return std::log(result.value);
}

static inline bool bawl_corr_row_active(const BAwLCorrSharedState& s, int row) {
  return !s.has_RACE || s.race_mask[row];
}

static inline double bawl_corr_log_component_survival(
    const BAwLCorrSharedState& s, const BAwLCorrPairData& pair, int j,
    double t, const double* const* cols, const ContextForRaceModels* ctx,
    bool numeric, BAwLCorrMomentStatus* status_out = nullptr) {
  BAwLCorrMomentStatus pair_status = BAwLCorrMomentStatus::ok;
  double out = bawl_corr_log_pair_survival(pair, t, cols, ctx, numeric,
                                           &pair_status);
  if (status_out != nullptr) *status_out = pair_status;
  if (pair_status == BAwLCorrMomentStatus::unstable ||
      pair_status == BAwLCorrMomentStatus::invalid) return R_NaN;
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!bawl_corr_row_active(s, row) || row == pair.row1 || row == pair.row2) continue;
    const double ls = bawl_corr_single_log_survival(t, row, cols, ctx);
    if (!R_FINITE(ls)) return R_NegInf;
    out += ls;
  }
  return out;
}

static inline double bawl_corr_log_component_cause(
    const BAwLCorrSharedState& s, const BAwLCorrPairData& pair, int j,
    int winner_row, double t, const double* const* cols,
    const ContextForRaceModels* ctx, bool numeric,
    BAwLCorrMomentStatus* status_out = nullptr) {
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  if (winner_row == pair.row1 || winner_row == pair.row2) {
    BAwLCorrMomentStatus ps = BAwLCorrMomentStatus::ok;
    double out = bawl_corr_log_pair_cause(pair, winner_row, t, cols, ctx,
                                          numeric, &ps);
    if (status_out != nullptr) *status_out = ps;
    if (ps == BAwLCorrMomentStatus::unstable ||
        ps == BAwLCorrMomentStatus::invalid) return R_NaN;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (!bawl_corr_row_active(s, row) || row == pair.row1 || row == pair.row2) continue;
      const double ls = bawl_corr_single_log_survival(t, row, cols, ctx);
      if (!R_FINITE(ls)) return R_NegInf;
      out += ls;
    }
    return out;
  }
  if (!bawl_corr_row_active(s, winner_row)) return R_NegInf;
  double out = bawl_corr_single_log_density(t, winner_row, cols, ctx);
  if (!R_FINITE(out)) return out;
  BAwLCorrMomentStatus ps = BAwLCorrMomentStatus::ok;
  const double lp = bawl_corr_log_pair_survival(pair, t, cols, ctx, numeric, &ps);
  if (status_out != nullptr) *status_out = ps;
  if (ps == BAwLCorrMomentStatus::unstable || ps == BAwLCorrMomentStatus::invalid)
    return R_NaN;
  if (!R_FINITE(lp)) return R_NegInf;
  out += lp;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!bawl_corr_row_active(s, row) || row == winner_row ||
        row == pair.row1 || row == pair.row2) continue;
    const double ls = bawl_corr_single_log_survival(t, row, cols, ctx);
    if (!R_FINITE(ls)) return R_NegInf;
    out += ls;
  }
  return out;
}

template <typename LogIntegrand>
static double bawl_corr_integrate_log_interval(double lo, double hi,
                                               LogIntegrand log_integrand) {
  if (!(hi > lo) && !(lo == 0.0 && hi == R_PosInf)) return R_NegInf;
  const GLRule& rule = gl_get_rule(64);
  double max_log = R_NegInf;
  std::array<double, 64> values{};
  for (int i = 0; i < 64; ++i) {
    const double x = rule.x[i];
    double t = 0.0;
    double log_jac = 0.0;
    if (hi == R_PosInf) {
      const double u = 0.5 * (x + 1.0);
      const double one = 1.0 - u;
      t = lo + u / one;
      log_jac = -std::log(2.0) - 2.0 * std::log(one);
    } else {
      t = 0.5 * (hi - lo) * x + 0.5 * (hi + lo);
      log_jac = std::log(0.5 * (hi - lo));
    }
    values[static_cast<size_t>(i)] = log_jac + std::log(rule.w[i]) +
      log_integrand(t);
    max_log = std::max(max_log, values[static_cast<size_t>(i)]);
  }
  if (!R_FINITE(max_log)) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < 64; ++i)
    sum += std::exp(values[static_cast<size_t>(i)] - max_log);
  return max_log + std::log(sum);
}

struct BAwLCorrExactTrialResult {
  double log_likelihood = R_NegInf;
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
  bool used_numeric = false;
  bool point_start = false;
};

static BAwLCorrExactTrialResult bawl_corr_exact_trial_loglik(
    const BAwLCorrSharedState& s, const BAwLCorrTrialLayout& layout, int j,
    const double* const* cols, const ContextForRaceModels* ctx,
    const Rcpp::IntegerVector& response, const Rcpp::LogicalVector& winner,
    double min_ll, bool force_numeric = false) {
  BAwLCorrExactTrialResult out;
  BAwLCorrPairData pair;
  if (!bawl_corr_make_pair_data(layout, s, cols, ctx, pair)) return out;
  const int start = j * s.n_lR;
  const int response_code = response[start];
  const bool known = response_code != NA_INTEGER;
  const int known_row = s.winner_row[static_cast<size_t>(j)];
  const double rt = s.rt[start];
  const double LT = s.LT[start];
  const double UT = s.UT[start];
  const double LC = s.LC[start];
  const double UC = s.UC[start];
  bool numeric = force_numeric;
  BAwLCorrMomentStatus st = BAwLCorrMomentStatus::ok;
  auto absorb_status = [&](BAwLCorrMomentStatus value) {
    if (value == BAwLCorrMomentStatus::invalid) st = value;
    else if (value == BAwLCorrMomentStatus::unstable &&
             st != BAwLCorrMomentStatus::invalid) st = value;
  };
  auto component_survival_mode = [&](double t, bool numeric_survival,
                                     BAwLCorrMomentStatus* status) {
    BAwLCorrMomentStatus local = BAwLCorrMomentStatus::ok;
    const double value = bawl_corr_log_component_survival(
        s, pair, j, t, cols, ctx, numeric_survival, &local);
    if (status != nullptr) *status = local;
    absorb_status(local);
    return value;
  };
  auto component_survival = [&](double t, BAwLCorrMomentStatus* status) {
    return component_survival_mode(t, numeric, status);
  };
  auto component_cause = [&](int row, double t, BAwLCorrMomentStatus* status) {
    BAwLCorrMomentStatus local = BAwLCorrMomentStatus::ok;
    const double value = bawl_corr_log_component_cause(
        s, pair, j, row, t, cols, ctx, numeric, &local);
    if (status != nullptr) *status = local;
    absorb_status(local);
    return value;
  };

  auto finite_cause = [&](double t, int row, BAwLCorrMomentStatus* st) {
    if (row >= 0) return component_cause(row, t, st);
    double total = R_NegInf;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
    for (int k = 0; k < n_lR_j; ++k) {
      const int r = start + k;
      if (!bawl_corr_row_active(s, r)) continue;
      total = log_sum_exp(total, component_cause(r, t, st));
    }
    return total;
  };
  auto interval_cause = [&](double lo, double hi, int row,
                            BAwLCorrMomentStatus* st) {
    return bawl_corr_integrate_log_interval(lo, hi, [&](double t) {
      return finite_cause(t, row, st);
    });
  };
  auto survival_difference = [&](double lo, double hi,
                                 BAwLCorrMomentStatus* st) {
    const double a = (lo == 0.0) ? 0.0 : component_survival(lo, st);
    if (hi == R_PosInf) return a;
    const double b = component_survival(hi, st);
    if (a == R_NegInf) return R_NegInf;
    const double z = log_diff_exp(a, b);
    return z;
  };
  auto stable_survival_difference = [&](double lo, double hi,
                                        BAwLCorrMomentStatus* st) {
    // The exact rectangle formula obtains a survivor by subtracting four BVN
    // CDF corners.  In a rare positive orthant those corners can agree to all
    // available digits while leaving a small, positive residual that looks
    // valid.  Dividing by that residual as an LT/UT normaliser then rewards
    // the parameter point with an arbitrarily large finite likelihood.
    //
    // The deterministic pair route integrates one marginal drift against the
    // conditional survivor of the other.  It is independent of the corner
    // derivative algebra and is evaluated only for the data-window
    // normaliser, which is cached by parameter cell below.  Event densities
    // retain the closed-form route.  Because both quantities use the same
    // positive-orthant normaliser, q cancels exactly in their log ratio.
    const double a = (lo == 0.0)
      ? 0.0 : component_survival_mode(lo, true, st);
    if (hi == R_PosInf) return a;
    const double b = component_survival_mode(hi, true, st);
    if (a == R_NegInf) return R_NegInf;
    return log_diff_exp(a, b);
  };

  double log_value = R_NegInf;
  if (R_FINITE(rt) && rt > 0.0) {
    log_value = finite_cause(rt, known ? known_row : -1, &st);
  } else if (rt == R_NegInf) {
    if (known) log_value = interval_cause(LT, LC, known_row, &st);
    else log_value = survival_difference(LT, LC, &st);
  } else if (rt == R_PosInf) {
    if (known) log_value = interval_cause(UC, UT, known_row, &st);
    else log_value = survival_difference(UC, UT, &st);
  } else if (Rcpp::NumericVector::is_na(rt)) {
    if (known) {
      log_value = log_sum_exp(interval_cause(LT, LC, known_row, &st),
                              interval_cause(UC, UT, known_row, &st));
    } else {
      log_value = log_sum_exp(survival_difference(LT, LC, &st),
                              survival_difference(UC, UT, &st));
    }
  }
  if (st == BAwLCorrMomentStatus::unstable || st == BAwLCorrMomentStatus::invalid) {
    out.status = st;
    return out;
  }

  if (s.has_trunc_trial[static_cast<size_t>(j)]) {
    // Z depends only on the parameter cell and the data-constant truncation
    // window, never on rt, so unique trials sharing a cell share Z.  Only
    // pure pair trials are cacheable: a singleton survivor's parameters
    // would otherwise need to join the key.
    const bool z_cacheable = !numeric && layout.n_active == 2;
    double zkey[15] = {0.0};
    const double* z_hit = nullptr;
    if (z_cacheable) {
      const double* t0c = cols[emc2col::bawl::t0];
      const double* Ac = cols[emc2col::bawl::A];
      const double* Bc = cols[emc2col::bawl::B];
      zkey[0] = t0c[pair.row1]; zkey[1] = Ac[pair.row1]; zkey[2] = Bc[pair.row1];
      zkey[3] = bawl_corr_row_k(cols, pair.row1, ctx);
      zkey[4] = pair.mu1; zkey[5] = pair.sd1;
      zkey[6] = t0c[pair.row2]; zkey[7] = Ac[pair.row2]; zkey[8] = Bc[pair.row2];
      zkey[9] = bawl_corr_row_k(cols, pair.row2, ctx);
      zkey[10] = pair.mu2; zkey[11] = pair.sd2;
      zkey[12] = pair.rho; zkey[13] = LT; zkey[14] = UT;
      for (const auto& e : s.z_cache) {
        if (std::memcmp(e.key, zkey, sizeof(zkey)) == 0) {
          z_hit = &e.log_z;
          break;
        }
      }
    }
    if (z_hit != nullptr) {
      log_value -= *z_hit;
    } else {
      BAwLCorrMomentStatus zst = BAwLCorrMomentStatus::ok;
      const double log_z = stable_survival_difference(LT, UT, &zst);
      if (zst == BAwLCorrMomentStatus::unstable || zst == BAwLCorrMomentStatus::invalid ||
          !R_FINITE(log_z)) {
        out.status = zst == BAwLCorrMomentStatus::ok
          ? BAwLCorrMomentStatus::unstable : zst;
        return out;
      }
      if (z_cacheable && s.z_cache.size() < 64) {
        BAwLCorrSharedState::ZCacheEntry e;
        std::memcpy(e.key, zkey, sizeof(zkey));
        e.log_z = log_z;
        s.z_cache.push_back(e);
      }
      log_value -= log_z;
    }
  }
  // Contaminant mixture; see src/contaminant_mixture.h.
  if (s.pc_col >= 0 || (s.pg_col >= 0 && s.guess.active())) {
    const double pC = (s.pc_col >= 0) ? cols[static_cast<size_t>(s.pc_col)][start] : 0.0;
    const double pG = (s.pg_col >= 0 && s.guess.active())
                        ? cols[static_cast<size_t>(s.pg_col)][start] : 0.0;
    log_value = mix_contaminants_rt(log_value, pC, pG, s.guess, rt, known);
  }
  out.log_likelihood = (!R_FINITE(log_value) || log_value < min_ll)
    ? min_ll : log_value;
  out.status = BAwLCorrMomentStatus::ok;
  out.used_numeric = numeric;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  for (int k = 0; k < n_lR_j; ++k) {
    const int r = start + k;
    if (bawl_corr_row_active(s, r) &&
        cols[emc2col::bawl::A][r] <= BAWL_A_EPS) out.point_start = true;
  }
  (void)winner;
  return out;
}

// Correlated BAwL likelihood using a single shared Gaussian factor.  The
// factor is integrated out with two batched Gauss-Hermite passes: a fixed
// wide scan locates each trial's integrand mass, then a small rule is
// re-centered per trial on those moments (see gh_quad.h).  Broad central
// trials reuse the scan integral; narrow or tail-peaked trials get the
// recentered pass.  This replaces the old fixed 40/80/200-node schedule with
// a default 12-node scan plus a selective 12-node refinement.
// Conditional on a factor value the accumulators are independent, so two
// evaluators supply the node integrands: a raw batched evaluator for the
// common all-finite untruncated design (no allocation, kernels called
// column-wise, z-invariants hoisted), and the ordinary race likelihood for
// everything else: we alter only v and sv, turn off its truncation
// correction, and let the existing code handle winners, omissions, clocks,
// censoring, contaminants, RACE masks, and expansion.  Truncation is
// normalised as
//   log int p(data | z) phi(z) dz - log int Z(z) phi(z) dz,
// rather than by averaging already-normalised node likelihoods.
double c_log_likelihood_bawl_correlated(
    BAwLCorrSharedState& cshared,
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
    RaceSharedState* shared,
    NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize) {
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  if (ctx == nullptr || !ctx->bawl_correlated || ctx->bawl_rho_index < 0 ||
      !cshared.valid || cshared.n_trials != n_trials || cshared.n_lR != n_lR) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: invalid correlated BAwL context.");
  }

  const int n_unique = cshared.n_unique;
  const int n_out = cshared.n_out;
  if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: trial_ll_out size mismatch.");
  }
  const int n_par = cshared.n_par;
  if (n_par > 64 || n_lR > 64) {
    Rcpp::stop("c_log_likelihood_bawl_correlated: at most 64 parameter columns and 64 accumulators are supported.");
  }

  // Role-mapped effective loadings; roles were resolved once in the shared
  // state, so per particle this is one pass over the raw rho column.
  const double* rho_col_ptr = cshared.cols[static_cast<size_t>(cshared.rho_col)];
  std::vector<double>& effective_rho = cshared.effective_rho;
  for (int r = 0; r < n_trials; ++r) {
    const double rho_cell = rho_col_ptr[r];
    effective_rho[static_cast<size_t>(r)] =
      cshared.role_correct[static_cast<size_t>(r)] ? std::fabs(rho_cell)
                                                   : rho_cell;
  }
  // Cell-keyed memoization is only valid within one particle: the ParamTable
  // base storage the column pointers view is refilled in place per particle.
  cshared.qab_cache.clear();
  cshared.z_cache.clear();

  const bool count_routes = bawl_corr_counters_enabled();
  bawl_corr_counters_active() = count_routes;
  BAwLCorrCounters& route_counters = bawl_corr_counters();

  // If all active loadings are zero, use the exact existing path.  Apart from
  // being cheaper, this keeps rho=0 bit-for-bit equivalent to independent
  // BAwL.  This is the only route below that consumes a materialised matrix
  // besides the generic-clock node evaluator.
  bool any_loading = false;
  for (int j = 0; j < cshared.n_unique && !any_loading; ++j) {
    const int start = j * n_lR;
    const int n_lR_j = cshared.has_RACE ? cshared.race_nacc[start] : n_lR;
    for (int k = 0; k < n_lR_j; ++k) {
      const int r = start + k;
      if (cshared.has_RACE && !cshared.race_mask[r]) continue;
      if (!isok[r]) continue;
      const double rho = effective_rho[static_cast<size_t>(r)];
      if (!R_FINITE(rho) || std::fabs(rho) > 1e-14) {
        any_loading = true;
        break;
      }
    }
  }
  if (!any_loading) {
    if (count_routes)
      route_counters.ordinary_zero_rho_trials += n_unique;
    return c_log_likelihood_race(materialize(), dadm, pdf1, cdf1, n_trials,
                                  winner, expand, min_ll, isok, n_lR,
                                  model_context_for_funcs, all_finite_trials,
                                  model_dfun_raw, model_pfun_raw, logS_at_t,
                                  shared, trial_ll_out, true);
  }

  const bool has_RACE_col = cshared.has_RACE;
  const Rcpp::IntegerVector& RACE = cshared.race_nacc;
  const Rcpp::LogicalVector& RACE_mask = cshared.race_mask;
  const Rcpp::NumericVector& LT = cshared.LT;
  const Rcpp::NumericVector& UT = cshared.UT;
  const bool has_truncation = !all_finite_trials && cshared.has_truncation;
  const std::vector<unsigned char>& has_trunc_trial = cshared.has_trunc_trial;

  // Raw batched conditional evaluator eligibility for the common design:
  // all-finite RTs, no truncation, no global-kill/timer/nogo machinery
  // outside the kernels.  Everything else falls back to the generic race
  // evaluator below.
  const bool use_fast_node_eval =
      all_finite_trials && !has_truncation &&
      model_dfun_raw != nullptr && model_pfun_raw != nullptr &&
      !ctx->has_global_kill() && !ctx->kill_active && !ctx->gng &&
      ctx->time_code == -1 && ctx->nogo_code == -1 &&
      (!dadm.containsElementNamed("RACE") || has_RACE_col);

  // Canonical per-particle trial classification.  Numerator routing still
  // follows the legacy scan/refine machinery below; the layout is the single
  // authority for loaded-row counts (positivity dimension) and future route
  // dispatch.
  const bool no_clock_model = !ctx->gng && !ctx->kill_active &&
    !ctx->has_global_kill() && ctx->time_code == -1 && ctx->nogo_code == -1;
  // The exact pair component is data-fixed and pointer-based, so it remains
  // eligible even when the caller disables the all-finite GH hint (that hint
  // controls only the fused node evaluator).  This also lets nonfinite and
  // truncated no-clock trials use the same component path.
  bawl_corr_classify_particle(cshared, isok, no_clock_model);
  if (count_routes) {
    for (int j = 0; j < n_unique; ++j) {
      const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
      if (L.n_loaded == 0) ++route_counters.loaded_dimension_0;
      else if (L.n_loaded == 1) ++route_counters.loaded_dimension_1;
      else if (L.n_loaded == 2) ++route_counters.loaded_dimension_2;
      else ++route_counters.loaded_dimension_3plus;
      if (L.route == BAwLCorrRoute::ordinary) {
        if (L.n_loaded == 0) ++route_counters.ordinary_zero_rho_trials;
        else ++route_counters.ordinary_single_loaded_trials;
      }
    }
  }

  // Exact pair trials are evaluated before the shared-factor fallback.  The
  // arrays are also masks for the node evaluator, so mixed particles retain
  // per-trial routing rather than promoting an exact pair because another
  // trial needs GH.  `R` is data-fixed and is only used to distinguish a
  // known response from an unknown winner/censoring observation.
  std::vector<unsigned char> exact_done(static_cast<size_t>(n_unique), 0);
  std::vector<double> exact_ll(static_cast<size_t>(n_unique), min_ll);

  // A trial with at most one loaded racer has no cross-racer dependence left:
  // the shared factor can be integrated out analytically and the trial is
  // exactly the ordinary independent BAwL race.  Evaluate those trials once
  // through the existing race likelihood, then keep them out of every GH
  // node.  This is deliberately done at the unique-trial level so the
  // ordinary path retains its complete handling of omissions, censoring,
  // truncation, Erlang clocks, contaminants, RACE masks, and expansion.
  //
  // The ordinary call receives a mask containing only ordinary trials.  Rows
  // belonging to correlated trials are invalid for that call and their
  // returned values are ignored; this avoids constructing a second DataFrame
  // or changing the data-fixed expansion/index bookkeeping.
  bool any_ordinary_trial = false;
  Rcpp::LogicalVector ordinary_isok(n_trials, false);
  for (int j = 0; j < n_unique; ++j) {
    const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
    if (L.route != BAwLCorrRoute::ordinary) continue;
    any_ordinary_trial = true;
    const int start = j * n_lR;
    const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (!has_RACE_col || RACE_mask[row]) ordinary_isok[row] = isok[row];
    }
  }
  if (any_ordinary_trial) {
    Rcpp::NumericMatrix pars_ordinary = materialize();
    Rcpp::NumericVector ordinary_ll(n_unique, min_ll);
    Rcpp::IntegerVector ordinary_expand(0);
    // c_log_likelihood_race updates a few context flags while inspecting the
    // active rows.  Isolate those updates from the correlated remainder,
    // whose GH evaluator must retain the original particle context.
    ContextForRaceModels ordinary_ctx = *ctx;
    c_log_likelihood_race(
        pars_ordinary, dadm, pdf1, cdf1, n_trials, winner, ordinary_expand,
        min_ll, ordinary_isok, n_lR, &ordinary_ctx, all_finite_trials,
        model_dfun_raw, model_pfun_raw, logS_at_t, shared, &ordinary_ll, true);
    for (int j = 0; j < n_unique; ++j) {
      const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
      if (L.route != BAwLCorrRoute::ordinary) continue;
      exact_done[static_cast<size_t>(j)] = 1;
      exact_ll[static_cast<size_t>(j)] = ordinary_ll[j];
    }
  }

  const Rcpp::IntegerVector response = dadm["R"];
  for (int j = 0; j < n_unique; ++j) {
    BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
    if (L.route != BAwLCorrRoute::exact_pair) continue;
    BAwLCorrExactTrialResult exact = bawl_corr_exact_trial_loglik(
        cshared, L, j, cshared.cols.data(), ctx, response, winner, min_ll, false);
    if (exact.status == BAwLCorrMomentStatus::unstable) {
      L.route = BAwLCorrRoute::numeric_pair;
      exact = bawl_corr_exact_trial_loglik(
          cshared, L, j, cshared.cols.data(), ctx, response, winner, min_ll, true);
    }
    if (exact.status == BAwLCorrMomentStatus::ok) {
      exact_done[static_cast<size_t>(j)] = 1;
      exact_ll[static_cast<size_t>(j)] = exact.log_likelihood;
      if (count_routes) {
        if (L.route == BAwLCorrRoute::numeric_pair)
          ++route_counters.numeric_pair_trials;
        else {
          ++route_counters.exact_pair_trials;
          if (L.winner_loaded) ++route_counters.exact_pair_pair_winner_trials;
          else ++route_counters.exact_pair_independent_winner_trials;
          if (exact.point_start) ++route_counters.exact_pair_point_start_trials;
        }
      }
    } else {
      // Both the exact rectangle and its numeric fallback found this trial's
      // geometry too degenerate to trust (e.g. a near-zero residual SD driven
      // by a runaway sv/rho combination).  That verdict is a property of the
      // parameter point, not of the kernel: routing it into the generic GH
      // quadrature would silently hand a numerically pathological integrand
      // to code with no equivalent stability check, which can return a large
      // finite value instead of failing.  Floor it directly instead of
      // giving the sampler a route to exploit that hole.
      exact_done[static_cast<size_t>(j)] = 1;
      exact_ll[static_cast<size_t>(j)] = min_ll;
      if (count_routes) ++route_counters.unstable_pair_floored_trials;
    }
  }

  int n_fallback_trials = 0;
  for (int j = 0; j < n_unique; ++j)
    if (!exact_done[static_cast<size_t>(j)]) ++n_fallback_trials;
  if (n_fallback_trials == 0) {
    // The exact component has already applied truncation, contamination, and
    // the final floor.  Do not enter the GH machinery when every unique trial
    // is a pair: this is the common canonical path and should pay no fallback
    // node or prepared-geometry work.
    double total_ll = 0.0;
    if (expand.length() > 0) {
      for (int i = 0; i < n_out; ++i) {
        const double value = exact_ll[static_cast<size_t>(expand[i] - 1)];
        if (trial_ll_out != nullptr) (*trial_ll_out)[i] = value;
        total_ll += value;
      }
    } else {
      for (int j = 0; j < n_unique; ++j) {
        const double value = exact_ll[static_cast<size_t>(j)];
        if (trial_ll_out != nullptr) (*trial_ll_out)[j] = value;
        total_ll += value;
      }
    }
    return total_ll;
  }

  const bool joint_posdrift = ctx->use_posdrift;
  // Numerator/denominator integrals per unique trial; expanded rows share
  // their unique trial's value and are written out at the end.
  std::vector<double> log_num(static_cast<size_t>(n_unique), R_NegInf);
  // Joint positive-drift conditioning needs a denominator for every trial:
  // integral phi(z) prod_k q_k(z) Z+(z) dz.  Without joint conditioning the
  // untruncated denominator is exactly one and only truncated trials need it.
  std::vector<double> log_den(static_cast<size_t>(n_unique),
                              joint_posdrift ? R_NegInf : 0.0);
  if (has_truncation && !joint_posdrift) {
    for (int j = 0; j < n_unique; ++j) {
      if (has_trunc_trial[static_cast<size_t>(j)])
        log_den[static_cast<size_t>(j)] = R_NegInf;
    }
  }

  // Positive-drift denominators without truncation have a closed form: the
  // Phi argument of racer k is a_k + b_k z, so
  //   int phi(z) prod_k Phi(a_k + b_k z) dz
  //     = Phi_K(a_k / sqrt(1 + b_k^2); rho_kl = b_k b_l / sqrt((1+b_k^2)(1+b_l^2))),
  // the one-factor orthant probability.  K <= 2 covers the dominant designs
  // exactly (norm_cdf_2d); larger K and truncated trials keep the quadrature
  // denominator below.  At high |rho| this integrand is phi against a steep
  // sigmoid, which fixed and moment-matched Hermite rules resolve worst, so
  // the closed form is an accuracy fix, not just a shortcut.
  // The positivity dimension is the layout's loaded-row count: exact-zero
  // loadings contribute the same constant q to numerator and denominator and
  // are cancelled algebraically on both sides, so an independent PM racer
  // neither raises the closed-form dimension nor forces the quadrature
  // denominator sweep.
  std::vector<unsigned char> den_by_quadrature(static_cast<size_t>(n_unique), 0);
  {
    const double* v0c = cshared.cols[static_cast<size_t>(emc2col::bawl::v)];
    const double* sv0c = cshared.cols[static_cast<size_t>(emc2col::bawl::sv)];
    for (int j = 0; j < n_unique; ++j) {
      if (has_trunc_trial[static_cast<size_t>(j)]) {
        if (exact_done[static_cast<size_t>(j)]) continue;
        den_by_quadrature[static_cast<size_t>(j)] = 1;
        continue;
      }
      if (!joint_posdrift) continue;  // denominator is exactly one
      const BAwLCorrTrialLayout& Lj = cshared.layout[static_cast<size_t>(j)];
      if (exact_done[static_cast<size_t>(j)]) continue;
      if (Lj.any_bad_row) {
        log_den[static_cast<size_t>(j)] = R_NegInf;
        continue;
      }
      if (Lj.n_loaded > 2) {
        den_by_quadrature[static_cast<size_t>(j)] = 1;
        continue;
      }
      double h[2] = {0.0, 0.0};
      double b[2] = {0.0, 0.0};
      bool bad = false;
      for (int m = 0; m < Lj.n_loaded; ++m) {
        const int row = Lj.loaded_row[m];
        const double rho = effective_rho[static_cast<size_t>(row)];
        const double sv = sv0c[row];
        const double svq =
          sv * std::fmax(std::sqrt(std::fmax(0.0, 1.0 - std::fabs(rho))), 1e-12);
        const double mu = v0c[row];
        if (!(svq > 0.0) || !R_FINITE(mu)) { bad = true; break; }
        const double slope =
          ((rho < 0.0) ? -1.0 : 1.0) * sv * std::sqrt(std::fabs(rho));
        const double bk = slope / svq;
        const double hk = (mu / svq) / std::sqrt(1.0 + bk * bk);
        if (!R_FINITE(bk) || !R_FINITE(hk)) { bad = true; break; }
        h[m] = hk;
        b[m] = bk;
      }
      if (bad) {
        log_den[static_cast<size_t>(j)] = R_NegInf;
        continue;
      }
      if (Lj.n_loaded == 0) {
        log_den[static_cast<size_t>(j)] = 0.0;
      } else if (Lj.n_loaded == 1) {
        log_den[static_cast<size_t>(j)] = pnorm_log_direct(h[0], true);
      } else {
        const double rho12 =
          b[0] * b[1] / std::sqrt((1.0 + b[0] * b[0]) * (1.0 + b[1] * b[1]));
        const double p = norm_cdf_2d(h[0], h[1], rho12);
        log_den[static_cast<size_t>(j)] =
          (R_FINITE(p) && p > 0.0) ? std::log(std::fmin(p, 1.0)) : R_NegInf;
      }
    }
  }

  // Positivity reweighting is evaluated on the node-shifted Gaussian means
  // and residual SDs.  RACE-inactive rows do not participate in the product,
  // and zero-loading rows are cancelled: their constant q appears in both
  // the numerator and denominator integrands, so it is omitted from both.
  auto log_positive_trial = [&](const Rcpp::NumericMatrix& pars_q,
                                const Rcpp::LogicalVector& isok_q,
                                int j) -> double {
    const int start = j * n_lR;
    const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
    double out = 0.0;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (has_RACE_col && !RACE_mask[row]) continue;
      if (!isok_q[row]) return R_NegInf;
      const double sd = pars_q(row, emc2col::bawl::sv);
      const double mu = pars_q(row, emc2col::bawl::v);
      if (!(sd > 0.0) || !R_FINITE(mu)) return R_NegInf;
      const double rho = effective_rho[static_cast<size_t>(row)];
      if (std::fabs(rho) <= 1e-14) continue;  // cancelled constant
      const double log_q = pnorm_log_direct(mu / sd, true);
      if (!R_FINITE(log_q)) return R_NegInf;
      out += log_q;
    }
    return out;
  };

  // Both node evaluators report per-unique-trial integrands; expanded output
  // rows share their unique trial's integral, so expansion happens once at
  // the end rather than per node.  j_to_i lets the generic evaluator read one
  // representative expanded row per unique trial.
  const std::vector<int>& j_to_i = cshared.j_to_i;

  // The generic-clock node evaluator is the only remaining consumer of a
  // materialised parameter matrix; produce it once per particle and only
  // when that route is actually selected.
  Rcpp::NumericMatrix pars_generic;
  if (!use_fast_node_eval) pars_generic = materialize();
  // Per-trial log(1 - pG) for the generic node evaluator, and the matching
  // suppression of pGuess inside the nested race call -- same split as the fast
  // path above: the constant factor is integrated, the guess mass is not.
  std::vector<double> generic_log1m_pg(static_cast<size_t>(n_unique), 0.0);
  if (!use_fast_node_eval && cshared.pg_col >= 0 && cshared.guess.active()) {
    for (int j = 0; j < n_unique; ++j) {
      const double pG = pars_generic(j * n_lR, cshared.pg_col);
      if (pG != 0.0) generic_log1m_pg[static_cast<size_t>(j)] = log1m(pG);
    }
    for (int r = 0; r < n_trials; ++r) pars_generic(r, cshared.pg_col) = 0.0;
  }

  // Step 4: the no-clock route consumes prepared BAwL geometry directly;
  // prepared conditional survivors retain scalar natural-branch semantics,
  // including the saturated-CDF boundary and separate q accounting.
  const double* fast_rt = nullptr;
  const double* fast_v0 = nullptr;
  const double* fast_sv0 = nullptr;
  const double* fast_B0 = nullptr;
  const double* fast_A0 = nullptr;
  const double* fast_t00 = nullptr;
  const double* fast_k0 = nullptr;
  std::vector<double> fast_res, fast_log1m_pc;
  std::vector<double> fast_const_ll;
  std::vector<BAwLPreparedRow> fast_prepared;
  std::vector<int> fast_winner, fast_loser, fast_const_winner, fast_const_loser, fast_ok;
  std::vector<int> fine_winner, fine_loser;
  bool use_fine_masks = false;
  bool has_fast_const = false;
  if (use_fast_node_eval) {
    fast_rt = cshared.rt.begin();
    fast_v0 = cshared.cols[static_cast<size_t>(emc2col::bawl::v)];
    fast_sv0 = cshared.cols[static_cast<size_t>(emc2col::bawl::sv)];
    fast_B0 = cshared.cols[static_cast<size_t>(emc2col::bawl::B)];
    fast_A0 = cshared.cols[static_cast<size_t>(emc2col::bawl::A)];
    fast_t00 = cshared.cols[static_cast<size_t>(emc2col::bawl::t0)];
    fast_k0 = ctx->bawl_k_fixed_zero
      ? nullptr : cshared.cols[static_cast<size_t>(emc2col::bawl::k)];
    fast_res.assign(static_cast<size_t>(n_trials), R_NegInf);
    fast_prepared.assign(static_cast<size_t>(n_trials), BAwLPreparedRow());
    fast_winner.resize(static_cast<size_t>(n_trials));
    fast_loser.resize(static_cast<size_t>(n_trials));
    fast_const_winner.assign(static_cast<size_t>(n_trials), 0);
    fast_const_loser.assign(static_cast<size_t>(n_trials), 0);
    fast_ok.resize(static_cast<size_t>(n_trials));
    for (int r = 0; r < n_trials; ++r) {
      if (exact_done[static_cast<size_t>(r / n_lR)]) {
        fast_ok[static_cast<size_t>(r)] = 0;
        fast_winner[static_cast<size_t>(r)] = 0;
        fast_loser[static_cast<size_t>(r)] = 0;
        fast_const_winner[static_cast<size_t>(r)] = 0;
        fast_const_loser[static_cast<size_t>(r)] = 0;
        continue;
      }
      const double rho = effective_rho[static_cast<size_t>(r)];
      const bool ok_r = isok[r] && R_FINITE(rho) && std::fabs(rho) <= 1.0;
      fast_ok[static_cast<size_t>(r)] = ok_r ? 1 : 0;
      const double magnitude = std::fabs(rho);
      const double loading = std::sqrt(magnitude);
      const double residual = std::sqrt(std::fmax(0.0, 1.0 - magnitude));
      const double sign = (rho < 0.0) ? -1.0 : 1.0;
      const bool active = !has_RACE_col || RACE_mask[r];
      // A valid exact-zero loading is independent of the shared factor. Its
      // race contribution is evaluated once below and added to every node;
      // only nonzero-loading rows remain in the node masks.
      const bool independent = ok_r && magnitude <= 1e-14;
      if (active && independent) has_fast_const = true;
      fast_winner[static_cast<size_t>(r)] =
        (active && winner[r] && !independent) ? 1 : 0;
      fast_loser[static_cast<size_t>(r)] =
        (active && n_lR > 1 && !winner[r] && !independent) ? 1 : 0;
      fast_const_winner[static_cast<size_t>(r)] =
        (active && winner[r] && independent) ? 1 : 0;
      fast_const_loser[static_cast<size_t>(r)] =
        (active && n_lR > 1 && !winner[r] && independent) ? 1 : 0;

      if (active && ok_r && !independent) {
        // All fixed-time LBA/leak work is prepared once per particle.  The
        // node loop below only updates v_q and evaluates the normal endpoints
        // through BAwLPreparedRow; it never rebuilds tau, leak factors, or
        // Jacobians and never stages a conditional parameter table.
        const double kval = ctx->bawl_k_fixed_zero ? 0.0 : fast_k0[r];
        const BAwLTimeGeometry g = bawl_time_geometry(
            fast_rt[r], fast_t00[r], fast_A0[r], fast_B0[r] + fast_A0[r], kval);
        fast_prepared[static_cast<size_t>(r)] = bawl_prepare_row(
            g, fast_v0[r], sign * fast_sv0[r] * loading,
            fast_sv0[r] * std::fmax(residual, 1e-12));
        if (count_routes) ++route_counters.prepared_rows;
      }
    }

    fast_const_ll.assign(static_cast<size_t>(n_unique), 0.0);
    if (has_fast_const) {
      // Evaluate exact-zero-loading race terms once.  These are ordinary
      // singleton likelihood factors, so their own positive-drift
      // normalisers remain in the singleton kernel and are not part of the
      // shared-factor denominator.
      model_dfun_raw(fast_rt, cshared.cols.data(), n_trials,
                   fast_const_winner.data(), fast_ok.data(), fast_res.data(),
                   R_NegInf, ctx);
      for (int j = 0; j < n_unique; ++j) {
        const int start = j * n_lR;
        const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
        double value = 0.0;
        for (int k = 0; k < n_lR_j; ++k) {
          const int row = start + k;
          if (has_RACE_col && !RACE_mask[row]) continue;
          if (fast_const_winner[static_cast<size_t>(row)])
            value += fast_res[static_cast<size_t>(row)];
        }
        fast_const_ll[static_cast<size_t>(j)] = value;
      }
      if (n_lR > 1) {
        model_pfun_raw(fast_rt, cshared.cols.data(), n_trials,
                     fast_const_loser.data(), fast_ok.data(), fast_res.data(),
                     R_NegInf, ctx);
        for (int j = 0; j < n_unique; ++j) {
          const int start = j * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
          for (int k = 0; k < n_lR_j; ++k) {
            const int row = start + k;
            if (has_RACE_col && !RACE_mask[row]) continue;
            if (fast_const_loser[static_cast<size_t>(row)])
              fast_const_ll[static_cast<size_t>(j)] += fast_res[static_cast<size_t>(row)];
          }
        }
      }
    }

    // The generic evaluator inherits the contaminant shift from the race
    // path; with all-finite RTs that is a per-trial log(1 - pC) added to the
    // numerator integrand.
    //
    // pGuess is split in two here, because unlike pC's log(1 - pC) the guess is
    // NOT a constant factor and must not be integrated over z.  Only its
    // (1 - pG) down-weight -- which IS a constant factor, and therefore rides
    // through both the quadrature and the joint-positivity division unchanged --
    // goes into the integrand; the flat guess MASS is added once to the finished
    // trial likelihood in the accumulation loop below.  The split is exact.
    fast_log1m_pc.assign(static_cast<size_t>(n_unique), 0.0);
    if (cshared.pc_col >= 0) {
      const double* pc_ptr = cshared.cols[static_cast<size_t>(cshared.pc_col)];
      for (int j = 0; j < n_unique; ++j) {
        const double pC = pc_ptr[j * n_lR];
        if (pC != 0.0) fast_log1m_pc[static_cast<size_t>(j)] = std::log1p(-pC);
      }
    }
    if (cshared.pg_col >= 0 && cshared.guess.active()) {
      const double* pg_ptr = cshared.cols[static_cast<size_t>(cshared.pg_col)];
      for (int j = 0; j < n_unique; ++j) {
        const double pG = pg_ptr[j * n_lR];
        if (pG != 0.0) fast_log1m_pc[static_cast<size_t>(j)] += log1m(pG);
      }
    }
  }

  // Evaluators take one latent value per unique trial so the refined pass
  // can center each trial's rule on its own integrand; a shared node is just
  // a constant z_by_trial.
  auto evaluate_node_fast = [&](const double* z_by_trial, bool need_num,
                                std::vector<double>& node_num,
                                std::vector<double>& node_den) {
    node_num.assign(static_cast<size_t>(n_unique), R_NegInf);
    node_den.assign(static_cast<size_t>(n_unique),
                    joint_posdrift ? R_NegInf : 0.0);
    const std::vector<int>& winner_mask = use_fine_masks ? fine_winner : fast_winner;
    const std::vector<int>& loser_mask = use_fine_masks ? fine_loser : fast_loser;
    for (int j = 0; j < n_unique; ++j) {
      if (exact_done[static_cast<size_t>(j)]) continue;
      const int start = j * n_lR;
      const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
      double s = has_fast_const ? fast_const_ll[static_cast<size_t>(j)] : 0.0;
      double lp = 0.0;
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (has_RACE_col && !RACE_mask[row]) continue;
        if (fast_const_winner[static_cast<size_t>(row)] ||
            fast_const_loser[static_cast<size_t>(row)]) continue;
        const double rho = effective_rho[static_cast<size_t>(row)];
        const bool loaded = fast_ok[static_cast<size_t>(row)] &&
          std::fabs(rho) > 1e-14;
        if (!loaded) {
          if (need_num && (winner_mask[static_cast<size_t>(row)] ||
                           loser_mask[static_cast<size_t>(row)]))
            s = R_NegInf;
          if (joint_posdrift) lp = R_NegInf;
          continue;
        }
        const BAwLPreparedRow& prepared = fast_prepared[static_cast<size_t>(row)];
        const double vq = prepared.v0 + prepared.slope * z_by_trial[j];
        if (joint_posdrift) {
          const double log_q = bawl_prepared_log_q(prepared, vq);
          if (!R_FINITE(log_q)) lp = R_NegInf;
          else if (lp != R_NegInf) lp += log_q;
        }
        if (need_num) {
          double log_component = R_NegInf;
          if (winner_mask[static_cast<size_t>(row)]) {
            log_component = bawl_prepared_log_pdf(prepared, vq);
          } else if (loser_mask[static_cast<size_t>(row)]) {
            log_component = bawl_prepared_log_survivor(
                prepared, vq, joint_posdrift);
          }
          s += log_component;
        }
      }
      if (need_num)
        node_num[static_cast<size_t>(j)] =
          // Loaded prepared components are already unnormalised: a winner
          // contributes f0 and a loser contributes q-F0.  The q product is
          // reserved for the positivity denominator, so adding lp here
          // would double-count every loaded loser.
          s + fast_log1m_pc[static_cast<size_t>(j)];
      if (joint_posdrift) node_den[static_cast<size_t>(j)] = lp;
    }
  };

  // Generic conditional evaluator: one shared-race pass per node vector.
  // need_num/need_den let the refinement passes skip the parts they do not
  // consume (the denominator pass needs no race likelihood, only the
  // positivity product and truncation normalisers).
  auto evaluate_node = [&](const double* z_by_trial, bool need_num,
                           bool need_den, std::vector<double>& node_num,
                           std::vector<double>& node_den) {
    node_num.assign(static_cast<size_t>(n_unique), R_NegInf);
    node_den.assign(static_cast<size_t>(n_unique),
                    joint_posdrift ? R_NegInf : 0.0);
    Rcpp::NumericMatrix pars_q = Rcpp::clone(pars_generic);
    Rcpp::LogicalVector isok_q = Rcpp::clone(isok);

    for (int r = 0; r < n_trials; ++r) {
      if (exact_done[static_cast<size_t>(r / n_lR)]) isok_q[r] = false;
    }

    // Conditional factor model preserving each marginal sv:
    // v_q = v + sign(rho) sv sqrt(|rho|) z,
    // sv_q = sv sqrt(1-|rho|).
    for (int r = 0; r < n_trials; ++r) {
      const double rho = effective_rho[static_cast<size_t>(r)];
      if (!isok_q[r]) continue;
      if (!R_FINITE(rho) || std::fabs(rho) > 1.0) {
        isok_q[r] = false;
        continue;
      }
      const double sv = pars_q(r, emc2col::bawl::sv);
      const double magnitude = std::fabs(rho);
      const double loading = std::sqrt(magnitude);
      const double residual = std::sqrt(std::fmax(0.0, 1.0 - magnitude));
      const double sign = (rho < 0.0) ? -1.0 : 1.0;
      pars_q(r, emc2col::bawl::v) += sign * sv * loading * z_by_trial[r / n_lR];
      // Keep the conditional SD strictly positive at the exact boundary. The
      // public bound stays just inside +/-1, but this avoids undefined scalar
      // kernels for hand-built parameter matrices at rho == +/-1.
      pars_q(r, emc2col::bawl::sv) = sv * std::fmax(residual, 1e-12);
    }

    std::vector<double> log_pos(static_cast<size_t>(n_unique), 0.0);
    if (joint_posdrift) {
      for (int j = 0; j < n_unique; ++j) {
        if (!exact_done[static_cast<size_t>(j)])
          log_pos[static_cast<size_t>(j)] = log_positive_trial(pars_q, isok_q, j);
      }
    }

    ContextForRaceModels node_ctx = *ctx;
    if (need_num) {
      Rcpp::NumericVector node_ll(n_out);
      c_log_likelihood_race(pars_q, dadm, pdf1, cdf1, n_trials,
                            winner, expand, R_NegInf, isok_q, n_lR,
                            &node_ctx, all_finite_trials,
                            model_dfun_raw, model_pfun_raw, logS_at_t,
                            shared, &node_ll, false);

      // Expanded rows of one unique trial share the integrand; report once.
      for (int j = 0; j < n_unique; ++j) {
        if (exact_done[static_cast<size_t>(j)]) continue;
        const int i = j_to_i[static_cast<size_t>(j)];
        if (i < 0) continue;
        node_num[static_cast<size_t>(j)] =
          log_pos[static_cast<size_t>(j)] + node_ll[i] +
          generic_log1m_pg[static_cast<size_t>(j)];
      }
    }

    if (!need_den) return;
    if (has_truncation) {
      // The factor must remain inside the truncation normaliser:
      //   Z = int [S_race(LT | z) - S_race(UT | z)] phi(z) dz.
      // Do not replace this with an unconditional race normaliser.  However,
      // conditional on a GH node the accumulators are independent, so the
      // conditional race survivors can be evaluated in a batch.  The old
      // implementation repacked every trial into row-major storage and called
      // get_trunc_normaliser_rowmajor_cpp once per trial and node.
      GslIntegrationControls gsl_ctl = default_gsl_controls();
      gsl_ctl.try_qng_first_finite = true;
      gsl_ctl.qag_key = GSL_INTEG_GAUSS21;
      gsl_ctl.rel_tol = 1e-5;
      static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
      GslWorkspacePtr& workspace = workspace_tls;

      std::vector<int> trunc_mask(static_cast<size_t>(n_unique), 0);
      std::vector<int> isok_q_int(static_cast<size_t>(n_trials), 0);
      double uniform_lt = 0.0;
      double uniform_ut = R_PosInf;
      bool have_trunc = false;
      bool uniform_lt_ok = true;
      bool uniform_ut_ok = true;
      for (int j = 0; j < n_unique; ++j) {
        const int start = j * n_lR;
        const double lt = LT[start];
        const double ut = UT[start];
        if (!has_trunc_trial[static_cast<size_t>(j)]) continue;
        trunc_mask[static_cast<size_t>(j)] = 1;
        if (!have_trunc) {
          uniform_lt = lt;
          uniform_ut = ut;
          have_trunc = true;
        } else {
          if (lt != uniform_lt) uniform_lt_ok = false;
          if (ut != uniform_ut) uniform_ut_ok = false;
        }
      }

      if (have_trunc) {
        for (int r = 0; r < n_trials; ++r) {
          isok_q_int[static_cast<size_t>(r)] = isok_q[r] ? 1 : 0;
        }

        std::vector<double> log_s_lt(static_cast<size_t>(n_unique), R_NegInf);
        std::vector<double> log_s_ut(static_cast<size_t>(n_unique), R_NegInf);
        std::array<const double*, 64> node_cols{};
        const double* pars_q_ptr = pars_q.begin();
        for (int c = 0; c < n_par; ++c) {
          node_cols[static_cast<size_t>(c)] =
            pars_q_ptr + static_cast<size_t>(c) * n_trials;
        }

        // Evaluate one conditional race survivor without constructing a
        // row-major trial buffer.  This is the fallback for variable-accumulator
        // RACE layouts, which the existing column-major batch callback cannot
        // represent because each trial can have a different active prefix.
        auto log_surv_trial = [&](double t, int j) -> double {
          const int start = j * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
          double log_s = 0.0;
          std::array<double, 64> par_row{};
          for (int k = 0; k < n_lR_j; ++k) {
            const int row = start + k;
            if (!isok_q[row] || (has_RACE_col && !RACE_mask[row])) return R_NegInf;
            for (int c = 0; c < n_par; ++c) par_row[static_cast<size_t>(c)] = pars_q(row, c);
            double cdf = cdf1(t, par_row.data(), &node_ctx);
            cdf = clamp_cdf01_race(cdf);
            const double log_s_k = safe_log1m_race(cdf);
            if (!emc2_isfinite(log_s_k)) return R_NegInf;
            log_s += log_s_k;
          }
          return log_s;
        };

        // Batch the common fixed-race case through the model-specific callback.
        // This is still conditional on the current GH node, so the shared
        // factor is fully represented in the final log-sum-exp below.
        auto log_surv_batch = [&](double t, std::vector<double>& out) {
          if (!has_RACE_col && logS_at_t != nullptr) {
            logS_at_t(t, node_cols.data(), n_trials, n_lR, n_par,
                      trunc_mask.data(), n_unique, isok_q_int.data(),
                      &node_ctx, out.data());
            return;
          }
          for (int j = 0; j < n_unique; ++j) {
            if (trunc_mask[static_cast<size_t>(j)]) {
              out[static_cast<size_t>(j)] = log_surv_trial(t, j);
            }
          }
        };

        if (uniform_lt_ok) {
          if (uniform_lt == 0.0) {
            for (int j = 0; j < n_unique; ++j) {
              if (trunc_mask[static_cast<size_t>(j)]) log_s_lt[static_cast<size_t>(j)] = 0.0;
            }
          } else {
            log_surv_batch(uniform_lt, log_s_lt);
          }
        } else {
          for (int j = 0; j < n_unique; ++j) {
            if (trunc_mask[static_cast<size_t>(j)]) {
              log_s_lt[static_cast<size_t>(j)] = log_surv_trial(LT[j * n_lR], j);
            }
          }
        }
        if (uniform_ut_ok && uniform_ut != R_PosInf) {
          log_surv_batch(uniform_ut, log_s_ut);
        } else if (!uniform_ut_ok) {
          for (int j = 0; j < n_unique; ++j) {
            if (trunc_mask[static_cast<size_t>(j)] && UT[j * n_lR] != R_PosInf) {
              log_s_ut[static_cast<size_t>(j)] = log_surv_trial(UT[j * n_lR], j);
            }
          }
        }

        std::vector<double> rowmajor(static_cast<size_t>(n_lR) * n_par);
        std::vector<int> ok_int(static_cast<size_t>(n_lR));
        const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
        for (int j = 0; j < n_unique; ++j) {
          if (!trunc_mask[static_cast<size_t>(j)]) continue;
          const int start = j * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
          const double lt = LT[start];
          const double ut = UT[start];
          double log_z = R_NegInf;
          bool need_scalar_fallback = false;

          if (R_FINITE(log_s_lt[static_cast<size_t>(j)])) {
            if (ut == R_PosInf) {
              // For BAwL's defective tail, +Inf is part of the retained
              // [LT, Inf) window and must not be subtracted.
              log_z = log_s_lt[static_cast<size_t>(j)];
            } else {
              log_z = log_diff_exp(log_s_lt[static_cast<size_t>(j)],
                                   log_s_ut[static_cast<size_t>(j)]);
              need_scalar_fallback =
                !(R_FINITE(log_z) && log_z > log_prob_eps);
            }
          }

          if (need_scalar_fallback) {
            for (int k = 0; k < n_lR_j; ++k) {
              const int row = start + k;
              ok_int[static_cast<size_t>(k)] =
                (isok_q[row] && (!has_RACE_col || RACE_mask[row])) ? 1 : 0;
              for (int c = 0; c < n_par; ++c) {
                rowmajor[static_cast<size_t>(k) * n_par + c] = pars_q(row, c);
              }
            }
            log_z = get_trunc_normaliser_rowmajor_cpp(
                rowmajor.data(), ok_int.data(), pdf1, cdf1, lt, ut,
                n_lR_j, n_par, gsl_ctl, &node_ctx, workspace);
          }

          if (has_trunc_trial[static_cast<size_t>(j)]) {
            node_den[static_cast<size_t>(j)] =
              log_pos[static_cast<size_t>(j)] + log_z;
          }
        }
      }
    }

    if (joint_posdrift) {
      // Untruncated trials still require the positivity normaliser.  For a
      // truncated trial the block above supplied log_z; do not add a second
      // unit-mass term for it.
      for (int j = 0; j < n_unique; ++j) {
        if (has_trunc_trial[static_cast<size_t>(j)]) continue;
        node_den[static_cast<size_t>(j)] = log_pos[static_cast<size_t>(j)];
      }
    }
  };

  auto eval_selected = [&](const double* z_by_trial, bool need_num,
                           bool need_den, std::vector<double>& num,
                           std::vector<double>& den) {
    if (use_fast_node_eval) {
      evaluate_node_fast(z_by_trial, need_num, num, den);
      if (count_routes)
        route_counters.fused_node_evaluations += n_fallback_trials;
    } else {
      evaluate_node(z_by_trial, need_num, need_den, num, den);
    }
  };

  // Two batched passes of Gauss-Hermite nodes.
  //
  // Scan pass: a fixed shared rule whose node masses give each trial's
  // integrand mean and SD.  The scan is also used directly for broad,
  // central trials; only difficult trials pay for the recentered pass.
  //
  // Refined pass: a smaller rule re-centered per trial on those moments.
  // A weak-drift trial whose winner-implied peak sits in a tail or is narrow
  // gets its nodes placed inside that peak.  Numerator and denominator
  // integrands peak in different places (the denominator mass sits near the
  // prior mode), so each gets its own centering; the denominator refinement
  // needs no race pass.
  double max_abs_rho = 0.0;
  for (int r = 0; r < n_trials; ++r) {
    if (!isok[r]) continue;
    const double rho = effective_rho[static_cast<size_t>(r)];
    if (R_FINITE(rho) && std::fabs(rho) <= 1.0)
      max_abs_rho = std::fmax(max_abs_rho, std::fabs(rho));
  }
  // The difficult part of the integral is the narrow factor-conditioned
  // likelihood near |rho| -> 1. Moderate correlations are smooth enough for
  // a 10-node scan/refinement pair; retain the original 12-node tier for the
  // sharp cases. Environment overrides remain available for convergence
  // tests and explicit user tuning.
  const int scan_default = (max_abs_rho <= 0.6) ? 10 : 12;
  const int fine_default = (max_abs_rho <= 0.6) ? 10 : 12;
  const GHRule& scan_rule = gh_rule(
      emc2_quad_nodes("EMC2_BAWLCORR_SCAN_N", scan_default));
  const GHRule& fine_rule = gh_rule(
      emc2_quad_nodes("EMC2_BAWLCORR_FINE_N", fine_default));
  const int n_scan = static_cast<int>(scan_rule.x.size());
  const int n_fine = static_cast<int>(fine_rule.x.size());
  const double sqrt2 = std::sqrt(2.0);
  bool any_quad_den = false;
  for (int j = 0; j < n_unique; ++j) {
    if (den_by_quadrature[static_cast<size_t>(j)]) { any_quad_den = true; break; }
  }

  if (count_routes) {
    long long n_quad_den = 0;
    for (int j = 0; j < n_unique; ++j)
      if (den_by_quadrature[static_cast<size_t>(j)]) ++n_quad_den;
    route_counters.den_quadrature_trials += n_quad_den;
    if (use_fast_node_eval) route_counters.gh_no_clock_trials += n_fallback_trials;
    else route_counters.gh_generic_clock_trials += n_fallback_trials;
    route_counters.scan_node_evaluations +=
      static_cast<long long>(n_scan) * n_fallback_trials;
  }

  std::vector<double> node_num;
  std::vector<double> node_den;
  std::vector<double> z_by_trial(static_cast<size_t>(n_unique));
  std::vector<double> z_scan(static_cast<size_t>(n_scan));
  std::vector<unsigned char> refine_num(static_cast<size_t>(n_unique), 0);
  for (int j = 0; j < n_unique; ++j)
    refine_num[static_cast<size_t>(j)] = exact_done[static_cast<size_t>(j)] ? 0 : 1;
  // Column-major (node fastest) log node masses from the scan pass.
  std::vector<double> scan_num(static_cast<size_t>(n_scan) * n_unique, R_NegInf);
  std::vector<double> scan_den(any_quad_den
      ? static_cast<size_t>(n_scan) * n_unique : 0, R_NegInf);
  for (int q = 0; q < n_scan; ++q) {
    const double z = sqrt2 * scan_rule.x[static_cast<size_t>(q)];
    z_scan[static_cast<size_t>(q)] = z;
    std::fill(z_by_trial.begin(), z_by_trial.end(), z);
    eval_selected(z_by_trial.data(), true, any_quad_den, node_num, node_den);
    const double log_w = std::log(gh_standard_normal_weight(scan_rule, q));
    for (int j = 0; j < n_unique; ++j) {
      scan_num[static_cast<size_t>(j) * n_scan + q] =
        log_w + node_num[static_cast<size_t>(j)];
      if (any_quad_den)
        scan_den[static_cast<size_t>(j) * n_scan + q] =
          log_w + node_den[static_cast<size_t>(j)];
    }
  }

  std::vector<AGHCenter> center_num(static_cast<size_t>(n_unique));
  std::vector<AGHCenter> center_den(
      any_quad_den ? static_cast<size_t>(n_unique) : 0);
  for (int j = 0; j < n_unique; ++j) {
    center_num[static_cast<size_t>(j)] = agh_center_from_scan(
        &scan_num[static_cast<size_t>(j) * n_scan], z_scan.data(), n_scan);
    if (any_quad_den)
      center_den[static_cast<size_t>(j)] = agh_center_from_scan(
          &scan_den[static_cast<size_t>(j) * n_scan], z_scan.data(), n_scan);
  }

  // A loaded known winner supplies an analytic proposal for the latent
  // factor.  Its marginal winner drift has density
  // (gamma0 + gamma1 V) N(V | mu, sv), restricted to the cause interval.  The
  // shared incomplete-normal moments turn that into the mean/variance of z;
  // survivor-only and unknown-winner integrands deliberately retain the scan
  // proposal because they have no winner tilt.
  std::vector<unsigned char> analytic_center_trial(static_cast<size_t>(n_unique), 0);
  if (use_fast_node_eval) {
    for (int j = 0; j < n_unique; ++j) {
      if (exact_done[static_cast<size_t>(j)]) continue;
      const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
      if (!L.winner_loaded || L.winner_row < 0) continue;
      const int row = L.winner_row;
      const BAwLTimeGeometry g = bawl_corr_row_geometry(
          cshared.rt[row], cshared.cols.data(), row, ctx);
      BAwLCorrRegion cause_region;
      if (!bawl_corr_cause_interval(g, joint_posdrift, cause_region)) continue;
      double p = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
      bawl_corr_univ_interval_moments3(
          cshared.cols[emc2col::bawl::v][row],
          cshared.cols[emc2col::bawl::sv][row], cause_region.lo,
          cause_region.hi, p, m1, m2, m3);
      const double den = g.gamma0 * p + g.gamma1 * m1;
      const double first = g.gamma0 * m1 + g.gamma1 * m2;
      const double second = g.gamma0 * m2 + g.gamma1 * m3;
      if (!(den > 0.0) || !R_FINITE(den) || !R_FINITE(first) ||
          !R_FINITE(second)) continue;
      const double mu = cshared.cols[emc2col::bawl::v][row];
      const double sd = cshared.cols[emc2col::bawl::sv][row];
      const double rho = effective_rho[static_cast<size_t>(row)];
      const double ev = first / den;
      const double vv = std::fmax(second / den - ev * ev, 0.0);
      const double loading = std::sqrt(std::fabs(rho));
      const double ez = (rho < 0.0 ? -1.0 : 1.0) * loading * (ev - mu) / sd;
      const double vz = 1.0 - std::fabs(rho) +
        std::fabs(rho) * vv / (sd * sd);
      if (R_FINITE(ez) && R_FINITE(vz) && vz > 0.0) {
        center_num[static_cast<size_t>(j)].mu = ez;
        center_num[static_cast<size_t>(j)].sigma =
          std::fmin(std::fmax(1.3 * std::sqrt(vz), 0.3), 1.25);
        analytic_center_trial[static_cast<size_t>(j)] = 1;
        if (count_routes) {
          ++route_counters.analytic_center_eligible_trials;
          ++route_counters.analytic_center_success_trials;
        }
      }
    }
  }

  auto fine_log_gauss_weights = [](const GHRule& rule) {
    std::vector<double> lgw(rule.x.size());
    for (size_t q = 0; q < rule.x.size(); ++q) {
      lgw[q] = std::log(rule.w[q]) + rule.x[q] * rule.x[q];
    }
    return lgw;
  };
  const std::vector<double> fine_lgw = fine_log_gauss_weights(fine_rule);

  std::vector<double> lw_by_trial(static_cast<size_t>(n_unique));
  // The raw scan is already a valid quadrature estimate when its mass is
  // broad and well inside the scan rule.  Keep the refinement masks separate
  // so easy trials do not re-run the BAwL kernels in the second pass; narrow
  // or tail-peaked trials retain the recentered pass for accuracy.
  if (use_fast_node_eval && !any_quad_den) {
    for (int j = 0; j < n_unique; ++j) {
      if (exact_done[static_cast<size_t>(j)]) continue;
      if (analytic_center_trial[static_cast<size_t>(j)]) {
        refine_num[static_cast<size_t>(j)] = 1;
        continue;
      }
      const AGHCenter& c = center_num[static_cast<size_t>(j)];
      const double* g = &scan_num[static_cast<size_t>(j) * n_scan];
      double peak = R_NegInf;
      int peak_i = -1;
      for (int q = 0; q < n_scan; ++q) {
        if (g[q] > peak) { peak = g[q]; peak_i = q; }
      }
      const bool at_scan_edge = peak_i == 0 || peak_i == n_scan - 1;
      // These thresholds are deliberately conservative: the normal scan is
      // reused only when it has a comfortably wide, central mass.  The fine
      // pass remains the fallback for the high-rho/narrow cases that motivated
      // the adaptive rule in the first place.
      const bool well_resolved = R_FINITE(peak) && !at_scan_edge &&
        std::fabs(c.mu) < 0.75 && c.sigma > 0.45;
      refine_num[static_cast<size_t>(j)] = well_resolved ? 0 : 1;
      if (!well_resolved) continue;
      double scan_ll = R_NegInf;
      for (int q = 0; q < n_scan; ++q)
        scan_ll = log_sum_exp(scan_ll, g[q]);
      log_num[static_cast<size_t>(j)] = scan_ll;
    }
    fine_winner = fast_winner;
    fine_loser = fast_loser;
    for (int j = 0; j < n_unique; ++j) {
      if (refine_num[static_cast<size_t>(j)]) continue;
      const int start = j * n_lR;
      const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
      for (int k = 0; k < n_lR_j; ++k) {
        fine_winner[static_cast<size_t>(start + k)] = 0;
        fine_loser[static_cast<size_t>(start + k)] = 0;
      }
    }
    use_fine_masks = true;
  }

  if (count_routes) {
    long long n_refined = 0;
    for (int j = 0; j < n_unique; ++j)
      if (refine_num[static_cast<size_t>(j)]) ++n_refined;
    route_counters.scan_refinement_trials += n_refined;
    route_counters.fine_node_evaluations +=
      static_cast<long long>(n_fine) * n_refined;
  }

  // The per-node raw evaluator consults these masks only during the optional
  // refinement pass; the scan still evaluates every active racer.
  for (int q = 0; q < n_fine; ++q) {
    const double x = sqrt2 * fine_rule.x[static_cast<size_t>(q)];
    for (int j = 0; j < n_unique; ++j) {
      const AGHCenter& c = center_num[static_cast<size_t>(j)];
      const double z = c.mu + c.sigma * x;
      z_by_trial[static_cast<size_t>(j)] = z;
      lw_by_trial[static_cast<size_t>(j)] =
        agh_log_weight(fine_lgw[static_cast<size_t>(q)], c.sigma, z);
    }
    eval_selected(z_by_trial.data(), true, false, node_num, node_den);
    for (int j = 0; j < n_unique; ++j) {
      if (!refine_num[static_cast<size_t>(j)]) continue;
      log_num[static_cast<size_t>(j)] = log_sum_exp(
          log_num[static_cast<size_t>(j)],
          lw_by_trial[static_cast<size_t>(j)] + node_num[static_cast<size_t>(j)]);
    }
  }

  if (any_quad_den) {
    // Quadrature denominators (truncated trials, or more than two racers)
    // carry the steep positivity sigmoid and the truncation window, so they
    // get a denser re-centered rule.  These evaluations skip the race pass
    // entirely and are cheap relative to the numerator sweeps.
    const GHRule& den_rule = gh_rule(64);
    const int n_den = static_cast<int>(den_rule.x.size());
    const std::vector<double> den_lgw = fine_log_gauss_weights(den_rule);
    if (count_routes) {
      long long n_quad_den = 0;
      for (int j = 0; j < n_unique; ++j)
        if (den_by_quadrature[static_cast<size_t>(j)]) ++n_quad_den;
      route_counters.den_quadrature_node_evaluations +=
        static_cast<long long>(n_den) * n_quad_den;
    }
    for (int q = 0; q < n_den; ++q) {
      const double x = sqrt2 * den_rule.x[static_cast<size_t>(q)];
      for (int j = 0; j < n_unique; ++j) {
        const AGHCenter& c = center_den[static_cast<size_t>(j)];
        const double z = c.mu + c.sigma * x;
        z_by_trial[static_cast<size_t>(j)] = z;
        lw_by_trial[static_cast<size_t>(j)] =
          agh_log_weight(den_lgw[static_cast<size_t>(q)], c.sigma, z);
      }
      eval_selected(z_by_trial.data(), false, true, node_num, node_den);
      for (int j = 0; j < n_unique; ++j) {
        if (!den_by_quadrature[static_cast<size_t>(j)]) continue;
        log_den[static_cast<size_t>(j)] = log_sum_exp(
            log_den[static_cast<size_t>(j)],
            lw_by_trial[static_cast<size_t>(j)] + node_den[static_cast<size_t>(j)]);
      }
    }
  }

  // The flat guess mass for quadrature trials.  Exact trials already carry the
  // full mixture (bawl_corr_exact_trial_loglik); quadrature trials carry only
  // its (1 - pG) factor, so the mass is added here, where `value` is a finished
  // trial likelihood rather than an integrand.  The guess term picks up
  // log(1 - pC) because the integrand it is being summed with already has it:
  //   (1-pC)[(1-pG) L + pG g]  ==  [(1-pC)(1-pG) L] + pG (1-pC) g.
  // A non-finite RT gets log_g = -Inf, leaving just the (1 - pG) factor already
  // in the integrand -- which is the correct treatment for an omission.
  const bool guess_quad =
    cshared.pg_col >= 0 && cshared.guess.active();
  auto add_guess_mass = [&](double value, int j) -> double {
    if (!guess_quad || exact_done[static_cast<size_t>(j)]) return value;
    const int start = j * n_lR;
    const double pG = cshared.cols[static_cast<size_t>(cshared.pg_col)][start];
    if (pG <= 0.0) return value;
    const double rt_j = cshared.rt[start];
    if (!R_FINITE(rt_j) || rt_j <= 0.0) return value;
    const double pC = (cshared.pc_col >= 0)
      ? cshared.cols[static_cast<size_t>(cshared.pc_col)][start] : 0.0;
    return log_sum_exp(value,
                       std::log(pG) + log1m(pC) + cshared.guess.log_g);
  };

  double total_ll = 0.0;
  if (expand.length() > 0) {
    for (int i = 0; i < n_out; ++i) {
      const int j = expand[i] - 1;
      double value = exact_done[static_cast<size_t>(j)]
        ? exact_ll[static_cast<size_t>(j)]
        : add_guess_mass(log_num[static_cast<size_t>(j)] -
                         log_den[static_cast<size_t>(j)], j);
      if (!R_FINITE(value) || value < min_ll) value = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = value;
      total_ll += value;
    }
  } else {
    for (int j = 0; j < n_unique; ++j) {
      double value = exact_done[static_cast<size_t>(j)]
        ? exact_ll[static_cast<size_t>(j)]
        : add_guess_mass(log_num[static_cast<size_t>(j)] -
                         log_den[static_cast<size_t>(j)], j);
      if (!R_FINITE(value) || value < min_ll) value = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = value;
      total_ll += value;
    }
  }
  return total_ll;
}

static RDMSWTNCorrSharedState build_rdmswtn_corr_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names) {
  RDMSWTNCorrSharedState s;
  if (n_lR <= 0 || n_trials <= 0 || n_trials % n_lR != 0) {
    Rcpp::stop("RDMSWTNcorr: invalid race dimensions.");
  }
  s.n_trials = n_trials;
  s.n_lR = n_lR;
  s.n_unique = n_trials / n_lR;
  s.n_out = expand.length() > 0 ? expand.length() : s.n_unique;
  s.n_par = keep_names.size();
  s.rho_col = ctx.rdmswtn_rho_index;
  if (s.rho_col < 0 || s.rho_col >= s.n_par) {
    Rcpp::stop("RDMSWTNcorr: rho column is outside the parameter table.");
  }
  s.rt = dadm["rt"];
  s.LT = get_col_with_default(dadm, "LT", 0.0);
  s.UT = get_col_with_default(dadm, "UT", R_PosInf);
  s.LC = get_col_with_default(dadm, "LC", R_PosInf);
  s.UC = get_col_with_default(dadm, "UC", 0.0);
  s.lR_codes = Rcpp::IntegerVector(dadm["lR"]);

  s.has_RACE = dadm.containsElementNamed("RACE");
  if (s.has_RACE && dadm.hasAttribute("RACE_nacc_by_row") &&
      dadm.hasAttribute("RACE_mask")) {
    s.race_nacc = dadm.attr("RACE_nacc_by_row");
    s.race_mask = dadm.attr("RACE_mask");
    if (s.race_nacc.size() != n_trials || s.race_mask.size() != n_trials)
      s.has_RACE = false;
  } else {
    s.has_RACE = false;
  }

  s.winner_row.assign(static_cast<size_t>(s.n_unique), -1);
  for (int r = 0; r < n_trials; ++r) {
    if (winner[r]) s.winner_row[static_cast<size_t>(r / n_lR)] = r;
  }
  s.cols.assign(static_cast<size_t>(std::max(s.n_par, 16)), nullptr);
  for (int c = 0; c < s.n_par; ++c) {
    const std::string nm = Rcpp::as<std::string>(keep_names[c]);
    s.cols[static_cast<size_t>(c)] =
      &table.base(0, table.base_index_for(nm));
    if (nm == "pContaminant") s.pc_col = c;
    else if (nm == "pGuess") s.pg_col = c;
  }
  s.guess = resolve_guess_kernel(dadm);
  s.guess.pg_col = s.pg_col;
  s.layout.assign(static_cast<size_t>(s.n_unique),
                  RDMSWTNCorrTrialLayout());
  s.valid = true;
  return s;
}

static inline bool rdmswtn_corr_row_active(
    const RDMSWTNCorrSharedState& s, int row) {
  return !s.has_RACE || s.race_mask[row];
}

static inline double rdmswtn_corr_eval_selected(
    double t, int row, const double* const* cols, int n_par,
    RacePdf1Fun fn, ContextForRaceModels* ctx) {
  if (n_par > 32) return R_NaN;
  std::array<double, 32> par{};
  for (int c = 0; c < n_par; ++c) par[static_cast<size_t>(c)] = cols[c][row];
  return fn(t, par.data(), ctx);
}

static inline double rdmswtn_corr_logcdf(
    double t, int row, const double* const* cols, int n_par,
    RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  if (!(t > 0.0)) return R_NegInf;
  const double value = rdmswtn_corr_eval_selected(
      t, row, cols, n_par, cdf1, ctx);
  if (ISNAN(value)) return R_NaN;
  if (!(value > 0.0)) return R_NegInf;
  if (value >= 1.0) return 0.0;
  return std::log(value);
}

static inline double rdmswtn_corr_logdensity(
    double t, int row, const double* const* cols, int n_par,
    RacePdf1Fun pdf1, ContextForRaceModels* ctx) {
  if (!(t > 0.0) || !R_FINITE(t)) return R_NegInf;
  const double value = rdmswtn_corr_eval_selected(
      t, row, cols, n_par, pdf1, ctx);
  if (ISNAN(value)) return R_NaN;
  return value > 0.0 ? std::log(value) : R_NegInf;
}

static inline double rdmswtn_corr_z_from_logcdf(double log_f) {
  if (log_f == R_NegInf) return R_NegInf;
  if (log_f >= 0.0) return R_PosInf;
  return R::qnorm(log_f, 0.0, 1.0, 1, 1);
}

template <typename LogIntegrand>
static double rdmswtn_corr_integrate_log(double lo, double hi,
                                         LogIntegrand log_integrand) {
  if (!(hi > lo) || lo == R_PosInf) return R_NegInf;
  const GLRule& rule = gl_get_rule(64);
  std::array<double, 64> values{};
  double max_log = R_NegInf;
  for (int i = 0; i < 64; ++i) {
    const double x = rule.x[static_cast<size_t>(i)];
    double point, log_jac;
    if (hi == R_PosInf) {
      const double u = 0.5 * (x + 1.0);
      const double one = 1.0 - u;
      point = lo + u / one;
      log_jac = -M_LN2 - 2.0 * std::log(one);
    } else {
      point = 0.5 * (hi - lo) * x + 0.5 * (hi + lo);
      log_jac = std::log(0.5 * (hi - lo));
    }
    const double value = log_jac +
      std::log(rule.w[static_cast<size_t>(i)]) + log_integrand(point);
    values[static_cast<size_t>(i)] = value;
    max_log = std::max(max_log, value);
  }
  if (!R_FINITE(max_log)) return R_NegInf;
  double sum = 0.0;
  for (double value : values) sum += std::exp(value - max_log);
  return max_log + std::log(sum);
}

static double rdmswtn_corr_log_bvn_upper(double z1, double z2, double rho) {
  if (z1 == R_NegInf)
    return z2 == R_NegInf ? 0.0 : log_normal_upper_tail(z2);
  if (z2 == R_NegInf) return log_normal_upper_tail(z1);
  if (z1 == R_PosInf || z2 == R_PosInf) return R_NegInf;
  const double p = norm_cdf_2d(-z1, -z2, rho);
  if (R_FINITE(p) && p > 1e-280) return std::log(std::fmin(1.0, p));

  // Direct conditional-normal integration is a tail fallback for probabilities
  // smaller than the natural-scale BVN routine can represent.
  if (z2 > z1) std::swap(z1, z2);
  const double sd = std::sqrt(std::fmax(1e-16, 1.0 - rho * rho));
  return rdmswtn_corr_integrate_log(z1, R_PosInf, [&](double x) {
    return -0.5 * x * x - 0.91893853320467274178 +
      log_normal_upper_tail((z2 - rho * x) / sd);
  });
}

static inline double rdmswtn_corr_log_pair_survival(
    double t, int row1, int row2, double rho, const double* const* cols,
    int n_par, RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  const double lf1 = rdmswtn_corr_logcdf(t, row1, cols, n_par, cdf1, ctx);
  const double lf2 = rdmswtn_corr_logcdf(t, row2, cols, n_par, cdf1, ctx);
  if (ISNAN(lf1) || ISNAN(lf2)) return R_NaN;
  return rdmswtn_corr_log_bvn_upper(
      rdmswtn_corr_z_from_logcdf(lf1),
      rdmswtn_corr_z_from_logcdf(lf2), rho);
}

static inline double rdmswtn_corr_log_single_survival(
    double t, int row, const double* const* cols, int n_par,
    RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  const double lf = rdmswtn_corr_logcdf(t, row, cols, n_par, cdf1, ctx);
  if (ISNAN(lf)) return R_NaN;
  if (lf == R_NegInf) return 0.0;
  if (lf >= 0.0) return R_NegInf;
  return log1m_exp(lf);
}

static double rdmswtn_corr_log_component_survival(
    const RDMSWTNCorrSharedState& s, int j,
    const RDMSWTNCorrTrialLayout& layout, double t, double rho,
    const double* const* cols, RaceCdf1Fun cdf1,
    ContextForRaceModels* ctx) {
  double out = rdmswtn_corr_log_pair_survival(
      t, layout.pair_row[0], layout.pair_row[1], rho, cols,
      s.n_par, cdf1, ctx);
  if (ISNAN(out)) return out;
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!rdmswtn_corr_row_active(s, row) ||
        row == layout.pair_row[0] || row == layout.pair_row[1]) continue;
    const double ls = rdmswtn_corr_log_single_survival(
        t, row, cols, s.n_par, cdf1, ctx);
    if (ISNAN(ls)) return ls;
    out += ls;
  }
  return out;
}

static double rdmswtn_corr_log_component_cause(
    const RDMSWTNCorrSharedState& s, int j,
    const RDMSWTNCorrTrialLayout& layout, int winner_row, double t,
    double rho, const double* const* cols, RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  const int pair1 = layout.pair_row[0], pair2 = layout.pair_row[1];
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  if (!rdmswtn_corr_row_active(s, winner_row)) return R_NegInf;
  double out = rdmswtn_corr_logdensity(
      t, winner_row, cols, s.n_par, pdf1, ctx);
  if (!R_FINITE(out)) return out;

  if (winner_row == pair1 || winner_row == pair2) {
    const int loser = winner_row == pair1 ? pair2 : pair1;
    const double lfw = rdmswtn_corr_logcdf(
        t, winner_row, cols, s.n_par, cdf1, ctx);
    const double lfl = rdmswtn_corr_logcdf(
        t, loser, cols, s.n_par, cdf1, ctx);
    if (ISNAN(lfw) || ISNAN(lfl)) return R_NaN;
    const double zw = rdmswtn_corr_z_from_logcdf(lfw);
    const double zl = rdmswtn_corr_z_from_logcdf(lfl);
    const double sd = std::sqrt(std::fmax(1e-16, 1.0 - rho * rho));
    out += log_normal_upper_tail((zl - rho * zw) / sd);
  } else {
    out += rdmswtn_corr_log_pair_survival(
        t, pair1, pair2, rho, cols, s.n_par, cdf1, ctx);
  }

  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!rdmswtn_corr_row_active(s, row) || row == winner_row ||
        row == pair1 || row == pair2) continue;
    out += rdmswtn_corr_log_single_survival(
        t, row, cols, s.n_par, cdf1, ctx);
  }
  return out;
}

static double rdmswtn_corr_trial_loglik(
    const RDMSWTNCorrSharedState& s, int j,
    const RDMSWTNCorrTrialLayout& layout,
    ContextForRaceModels* ctx, const Rcpp::IntegerVector& response,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, double min_ll) {
  if (!layout.params_ok) return min_ll;
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  const double rho =
    s.cols[static_cast<size_t>(s.rho_col)][layout.pair_row[0]];
  const bool known = response[start] != NA_INTEGER;
  const int known_row = layout.winner_row;
  const double rt = s.rt[start];
  const double LT = s.LT[start], UT = s.UT[start];
  const double LC = s.LC[start], UC = s.UC[start];
  if (known && known_row < 0) return min_ll;

  int time_row = -1, nogo_row = -1, n_resp = 0;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!rdmswtn_corr_row_active(s, row)) continue;
    const int code = s.lR_codes[row];
    if (code == ctx->time_code) time_row = row;
    else if (code == ctx->nogo_code) nogo_row = row;
    else ++n_resp;
  }

  auto single_cause = [&](int row, double t) {
    return rdmswtn_corr_log_component_cause(
        s, j, layout, row, t, rho, s.cols.data(), pdf1, cdf1, ctx);
  };
  auto finite_cause = [&](int row, double t) {
    if (row < 0) {
      double total = R_NegInf;
      for (int k = 0; k < n_lR_j; ++k) {
        const int active_row = start + k;
        if (!rdmswtn_corr_row_active(s, active_row)) continue;
        const int code = s.lR_codes[active_row];
        if (code == ctx->time_code || code == ctx->nogo_code) continue;
        double cause = single_cause(active_row, t);
        if (time_row >= 0 && n_resp > 0) {
          cause = log_sum_exp(
              cause, single_cause(time_row, t) - std::log(n_resp));
        }
        total = log_sum_exp(total, cause);
      }
      return total;
    }
    double value = single_cause(row, t);
    if (time_row >= 0 && time_row != row && n_resp > 0 &&
        s.lR_codes[row] != ctx->nogo_code) {
      value = log_sum_exp(value,
                          single_cause(time_row, t) - std::log(n_resp));
    }
    return value;
  };
  auto interval_cause = [&](double lo, double hi, int row) {
    return rdmswtn_corr_integrate_log(lo, hi, [&](double t) {
      return finite_cause(row, t);
    });
  };
  auto survivor = [&](double t) {
    if (t == 0.0) return 0.0;
    return rdmswtn_corr_log_component_survival(
        s, j, layout, t, rho, s.cols.data(), cdf1, ctx);
  };
  auto survivor_difference = [&](double lo, double hi) {
    const double a = lo == 0.0 ? 0.0 : survivor(lo);
    if (hi == R_PosInf) return a;
    return log_diff_exp(a, survivor(hi));
  };

  double value = R_NegInf;
  if (R_FINITE(rt) && rt > 0.0) {
    value = finite_cause(known ? known_row : -1, rt);
  } else if (rt == R_NegInf) {
    if (!known && nogo_row >= 0) {
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (!rdmswtn_corr_row_active(s, row) || row == nogo_row) continue;
        value = log_sum_exp(value, interval_cause(LT, LC, row));
      }
    } else {
      value = known ? interval_cause(LT, LC, known_row)
                    : survivor_difference(LT, LC);
    }
  } else if (rt == R_PosInf) {
    if (nogo_row >= 0) {
      value = log_sum_exp(interval_cause(LT, UC, nogo_row), survivor(UC));
    } else {
      value = known ? interval_cause(UC, UT, known_row)
                    : survivor_difference(UC, UT);
    }
  } else if (Rcpp::NumericVector::is_na(rt)) {
    if (known) {
      value = log_sum_exp(interval_cause(LT, LC, known_row),
                          interval_cause(UC, UT, known_row));
    } else {
      value = log_sum_exp(survivor_difference(LT, LC),
                          survivor_difference(UC, UT));
    }
  }

  if (LT != 0.0 || UT != R_PosInf) {
    const double log_z = survivor_difference(LT, UT);
    if (!R_FINITE(log_z)) return min_ll;
    value -= log_z;
  }
  // Contaminant mixture; see src/contaminant_mixture.h.
  if (s.pc_col >= 0) {
    const double p_c = s.cols[static_cast<size_t>(s.pc_col)][start];
    if (!R_FINITE(p_c) || p_c < 0.0 || p_c >= 1.0) return min_ll;
  }
  if (s.pc_col >= 0 || (s.pg_col >= 0 && s.guess.active())) {
    const double pC = (s.pc_col >= 0) ? s.cols[static_cast<size_t>(s.pc_col)][start] : 0.0;
    const double pG = (s.pg_col >= 0 && s.guess.active())
                        ? s.cols[static_cast<size_t>(s.pg_col)][start] : 0.0;
    value = mix_contaminants_rt(value, pC, pG, s.guess, rt, known);
  }
  return (!R_FINITE(value) || value < min_ll) ? min_ll : value;
}

double c_log_likelihood_rdmswtn_correlated(
    RDMSWTNCorrSharedState& s, Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, const int n_trials,
    LogicalVector winner, Rcpp::IntegerVector expand, double min_ll,
    const Rcpp::LogicalVector isok, int n_lR,
    void* model_context_for_funcs, bool all_finite_trials,
    RaceRawFun model_dfun_raw, RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t, RaceSharedState* shared,
    NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize) {
  auto* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  if (ctx == nullptr || !ctx->rdmswtn_correlated ||
      ctx->rdmswtn_rho_index < 0 || !s.valid ||
      s.n_trials != n_trials || s.n_lR != n_lR) {
    Rcpp::stop("RDMSWTNcorr: invalid correlated likelihood context.");
  }
  const double* rho_col = s.cols[static_cast<size_t>(s.rho_col)];
  bool any_correlated = false;
  bool any_ordinary = false;
  Rcpp::LogicalVector ordinary_isok(n_trials, false);

  for (int j = 0; j < s.n_unique; ++j) {
    RDMSWTNCorrTrialLayout& layout = s.layout[static_cast<size_t>(j)];
    layout = RDMSWTNCorrTrialLayout();
    layout.winner_row = s.winner_row[static_cast<size_t>(j)];
    const int start = j * n_lR;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : n_lR;
    double first_rho = 0.0;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (!rdmswtn_corr_row_active(s, row)) continue;
      const double rho = rho_col[row];
      if (!R_FINITE(rho) || std::fabs(rho) > 1.0) {
        Rcpp::stop("RDMSWTNcorr requires finite natural-scale rho values in [-1, 1].");
      }
      if (!ctx->use_posdrift && std::fabs(rho) > 1e-12) {
        Rcpp::stop("RDMSWTNcorr with posdrift = FALSE is supported only when every active rho is zero.");
      }
      if (!isok[row]) layout.params_ok = false;
      if (std::fabs(rho) <= 1e-12) continue;
      if (layout.n_nonzero < 2)
        layout.pair_row[layout.n_nonzero] = row;
      ++layout.n_nonzero;
      if (layout.n_nonzero == 1) first_rho = rho;
      else if (std::fabs(rho - first_rho) > 1e-12) {
        Rcpp::stop("RDMSWTNcorr requires the two participating rows to have the same signed nonzero rho; all opted-out rows must have rho = 0.");
      }
    }
    if (layout.n_nonzero > 2) {
      Rcpp::stop("RDMSWTNcorr supports one pair per trial: at most two active rows may have nonzero rho, and all other rows must have rho = 0.");
    }
    layout.ordinary = layout.n_nonzero <= 1;
    if (layout.ordinary) {
      any_ordinary = true;
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (rdmswtn_corr_row_active(s, row)) ordinary_isok[row] = isok[row];
      }
    } else {
      any_correlated = true;
    }
  }

  // Exact nesting and the cheapest route for an all-zero (or otherwise
  // unpaired) particle.
  if (!any_correlated) {
    return c_log_likelihood_race(
        materialize(), dadm, pdf1, cdf1, n_trials, winner, expand, min_ll,
        isok, n_lR, model_context_for_funcs, all_finite_trials,
        model_dfun_raw, model_pfun_raw, logS_at_t, shared, trial_ll_out, true);
  }

  std::vector<double> unique_ll(static_cast<size_t>(s.n_unique), min_ll);
  if (any_ordinary) {
    Rcpp::NumericVector ordinary_ll(s.n_unique, min_ll);
    Rcpp::IntegerVector no_expand(0);
    ContextForRaceModels ordinary_ctx = *ctx;
    ordinary_ctx.rdmswtn_correlated = false;
    c_log_likelihood_race(
        materialize(), dadm, pdf1, cdf1, n_trials, winner, no_expand, min_ll,
        ordinary_isok, n_lR, &ordinary_ctx, all_finite_trials,
        model_dfun_raw, model_pfun_raw, logS_at_t, shared, &ordinary_ll, true);
    for (int j = 0; j < s.n_unique; ++j) {
      if (s.layout[static_cast<size_t>(j)].ordinary)
        unique_ll[static_cast<size_t>(j)] = ordinary_ll[j];
    }
  }

  const Rcpp::IntegerVector response = dadm["R"];
  for (int j = 0; j < s.n_unique; ++j) {
    const RDMSWTNCorrTrialLayout& layout =
      s.layout[static_cast<size_t>(j)];
    if (!layout.ordinary) {
      unique_ll[static_cast<size_t>(j)] = rdmswtn_corr_trial_loglik(
          s, j, layout, ctx, response, pdf1, cdf1, min_ll);
    }
  }

  double total = 0.0;
  if (expand.length() > 0) {
    for (int i = 0; i < s.n_out; ++i) {
      const double value = unique_ll[static_cast<size_t>(expand[i] - 1)];
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = value;
      total += value;
    }
  } else {
    for (int j = 0; j < s.n_unique; ++j) {
      const double value = unique_ll[static_cast<size_t>(j)];
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = value;
      total += value;
    }
  }
  return total;
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

// Test/benchmark observability accessors for the LogicalRules correlated
// capacity route.  Counting is active only while EMC2_LRCAP_COUNTERS is set;
// the accessors themselves always return the current counters.
// [[Rcpp::export]]
Rcpp::List lr_capacity_counter_values() {
  const LrCapacityCounters& c = lr_capacity_counters();
  return Rcpp::List::create(
      Rcpp::Named("ordinary_trials") = (double)c.ordinary_trials,
      Rcpp::Named("capacity_detection_trials") = (double)c.capacity_detection_trials,
      Rcpp::Named("capacity_choice_trials") = (double)c.capacity_choice_trials,
      Rcpp::Named("invalid_trials") = (double)c.invalid_trials,
      Rcpp::Named("tau_zero_trials") = (double)c.tau_zero_trials,
      Rcpp::Named("factor_node_evaluations") = (double)c.factor_node_evaluations,
      Rcpp::Named("channel_gl_evaluations") = (double)c.channel_gl_evaluations);
}

// [[Rcpp::export]]
void lr_capacity_counters_reset() {
  lr_capacity_counters().reset();
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
  const double magnitude = std::fabs(rho);
  const double slope =
    ((rho < 0.0) ? -1.0 : 1.0) * sv * std::sqrt(magnitude);
  const double sv_res =
    sv * std::fmax(std::sqrt(std::fmax(0.0, 1.0 - magnitude)), 1e-12);
  const BAwLTimeGeometry g = bawl_time_geometry(t, t0, A, B + A, k);
  const BAwLPreparedRow r = bawl_prepare_row(g, v, slope, sv_res);
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

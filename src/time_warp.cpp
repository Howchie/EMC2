#include "time_warp.h"

#include <algorithm>
#include <vector>

void configure_time_warp_context(RaceModelAdapter& adapter,
                                 const Rcpp::CharacterVector& keep_names,
                                 const std::string& caller) {
  int idx = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == "eta") {
      idx = j;
      break;
    }
  }
  if (idx < 0) return;                       // no warp column
  if (!adapter.ctx.tw.supported) {
    Rcpp::stop("%s: the design supplies an 'eta' column but this model does not "
               "implement the operational-time warp (Math/ballistic-time.md).",
               caller.c_str());
  }
  if (adapter.ctx.t0_index < 0) {
    Rcpp::stop("%s: the operational-time warp requires a purely additive t0 "
               "column.", caller.c_str());
  }

  adapter.ctx.tw.eta_index  = idx;
  adapter.ctx.tw.base_pdf1  = adapter.pdf1_ptr;
  adapter.ctx.tw.base_cdf1  = adapter.cdf1_ptr;
  adapter.ctx.tw.base_d_raw = adapter.model_dfun_raw;
  adapter.ctx.tw.base_p_raw = adapter.model_pfun_raw;
  adapter.ctx.tw.base_logS  = adapter.logS_at_t_ptr;

  if (adapter.pdf1_ptr       != nullptr) adapter.pdf1_ptr       = &tw_pdf1;
  if (adapter.cdf1_ptr       != nullptr) adapter.cdf1_ptr       = &tw_cdf1;
  if (adapter.model_dfun_raw != nullptr) adapter.model_dfun_raw = &tw_d_raw;
  if (adapter.model_pfun_raw != nullptr) adapter.model_pfun_raw = &tw_p_raw;
  if (adapter.logS_at_t_ptr  != nullptr) adapter.logS_at_t_ptr  = &tw_logS_at_t;
}

double tw_pdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double t0  = par[ctx->t0_index];
  const double eta = par[tw.eta_index];
  const double u   = t - t0;
  if (eta == 0.0 || !(u > 0.0) || !R_FINITE(u))
    return tw.base_pdf1(t, par, ctx_);        // includes t == +Inf
  const double p = tw.base_pdf1(t0 + emc2tw::fwd(u, eta), par, ctx_);
  if (!(p > 0.0) || !R_FINITE(p)) return 0.0;
  const double log_density = std::log(p) + emc2tw::log_jac(u, eta);
  return R_FINITE(log_density) ? std::exp(log_density) : 0.0;
}

double tw_cdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double t0  = par[ctx->t0_index];
  const double eta = par[tw.eta_index];
  const double u   = t - t0;
  if (eta == 0.0 || !(u > 0.0) || !R_FINITE(u))
    return tw.base_cdf1(t, par, ctx_);        // defect at +Inf reaches base
  return tw.base_cdf1(t0 + emc2tw::fwd(u, eta), par, ctx_);
}

void tw_d_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double* t0_  = cols[ctx->t0_index];
  const double* eta_ = cols[tw.eta_index];

  bool any = false;
  if (eta_ != nullptr && t0_ != nullptr) {
    for (int i = 0; i < n_rows && !any; ++i) {
      if (mask[i] && eta_[i] != 0.0) any = true;
    }
  }
  if (!any) {                                  // exact parent path
    tw.base_d_raw(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
    return;
  }

  static thread_local std::vector<double> rt_buf;
  static thread_local std::vector<double> lj_buf;
  rt_buf.assign(rt, rt + n_rows);
  lj_buf.assign(static_cast<size_t>(n_rows), 0.0);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double eta = eta_[i];
    if (eta == 0.0) continue;
    const double u = rt[i] - t0_[i];
    if (!(u > 0.0) || !R_FINITE(u)) continue;  // guards match base clock
    rt_buf[i] = t0_[i] + emc2tw::fwd(u, eta);
    lj_buf[i] = emc2tw::log_jac(u, eta);
  }

  const bool floor_raw = ctx->floor_raw_log_lik;
  ctx->floor_raw_log_lik = false;              // Jacobian precedes flooring
  tw.base_d_raw(rt_buf.data(), cols, n_rows, mask, isok, out, min_ll, ctx_);
  ctx->floor_raw_log_lik = floor_raw;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    out[i] = raw_log_value(out[i] + lj_buf[i], min_ll, floor_raw);
  }
}

void tw_p_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double* t0_  = cols[ctx->t0_index];
  const double* eta_ = cols[tw.eta_index];

  bool any = false;
  if (eta_ != nullptr && t0_ != nullptr) {
    for (int i = 0; i < n_rows && !any; ++i) {
      if (mask[i] && eta_[i] != 0.0) any = true;
    }
  }
  if (!any) {                                  // exact parent path
    tw.base_p_raw(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
    return;
  }

  static thread_local std::vector<double> rt_buf;
  rt_buf.assign(rt, rt + n_rows);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double eta = eta_[i];
    if (eta == 0.0) continue;
    const double u = rt[i] - t0_[i];
    if (!(u > 0.0) || !R_FINITE(u)) continue;
    rt_buf[i] = t0_[i] + emc2tw::fwd(u, eta);
  }
  tw.base_p_raw(rt_buf.data(), cols, n_rows, mask, isok, out, min_ll, ctx_);
}

void tw_logS_at_t(double t, const double* const* cols, int n_rows_total,
                  int n_lR, int n_par, const int* trunc_mask,
                  int n_unique_trials, const int* isok_all, void* ctx_,
                  double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double* t0_  = cols[ctx->t0_index];
  const double* eta_ = cols[tw.eta_index];

  bool any = false;
  if (eta_ != nullptr && t0_ != nullptr && R_FINITE(t)) {
    for (int r = 0; r < n_rows_total && !any; ++r) {
      if (eta_[r] != 0.0) any = true;
    }
  }
  if (!any) {
    tw.base_logS(t, cols, n_rows_total, n_lR, n_par, trunc_mask,
                 n_unique_trials, isok_all, ctx_, logS_out);
    return;
  }

  static thread_local std::vector<double> t0_buf;
  static thread_local std::vector<const double*> cols_buf;
  t0_buf.assign(t0_, t0_ + n_rows_total);
  for (int r = 0; r < n_rows_total; ++r) {
    const double eta = eta_[r];
    if (eta == 0.0) continue;
    const double u = t - t0_[r];
    if (!(u > 0.0)) continue;                  // preserve not-started branch
    t0_buf[r] = t - emc2tw::fwd(u, eta);       // t - t0_buf[r] == c_eta(u)
  }

  // Kernels may fetch optional trailing pointers for this model variant.  Do
  // not read beyond the compacted n_par source columns while padding them.
  const int n_slots = std::max(n_par, 16);
  cols_buf.assign(static_cast<size_t>(n_slots), nullptr);
  for (int j = 0; j < n_par; ++j) cols_buf[j] = cols[j];
  cols_buf[ctx->t0_index] = t0_buf.data();

  tw.base_logS(t, cols_buf.data(), n_rows_total, n_lR, n_par, trunc_mask,
               n_unique_trials, isok_all, ctx_, logS_out);
}

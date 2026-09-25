#include "time_warp.h"

#include <algorithm>
#include <vector>

using emc2tw::clock_fwd;
using emc2tw::clock_is_identity;
using emc2tw::clock_is_valid;
using emc2tw::clock_log_jac;

void configure_time_warp_context(RaceModelAdapter& adapter,
                                 const Rcpp::CharacterVector& keep_names,
                                 const std::string& caller) {
  TimeWarpPlan& tw = adapter.ctx.tw;
  int idx = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == tw.par_name) {
      idx = j;
      break;
    }
  }
  if (idx < 0) {                             // no clock column
    if (tw.required) {
      Rcpp::stop("%s: this model requires the clock parameter column '%s'.",
                 caller.c_str(), tw.par_name);
    }
    return;
  }
  if (!tw.supported) {
    Rcpp::stop("%s: the design supplies an '%s' column but this model does not "
               "implement the operational-time warp (Math/ballistic-time.md).",
               caller.c_str(), tw.par_name);
  }
  if (adapter.ctx.t0_index < 0) {
    Rcpp::stop("%s: the operational-time warp requires a purely additive t0 "
               "column.", caller.c_str());
  }

  tw.par_index  = idx;
  tw.base_pdf1  = adapter.pdf1_ptr;
  tw.base_cdf1  = adapter.cdf1_ptr;
  tw.base_d_raw = adapter.model_dfun_raw;
  tw.base_p_raw = adapter.model_pfun_raw;
  tw.base_logS  = adapter.logS_at_t_ptr;

  if (adapter.pdf1_ptr       != nullptr) adapter.pdf1_ptr       = &tw_pdf1;
  if (adapter.cdf1_ptr       != nullptr) adapter.cdf1_ptr       = &tw_cdf1;
  if (adapter.model_dfun_raw != nullptr) adapter.model_dfun_raw = &tw_d_raw;
  if (adapter.model_pfun_raw != nullptr) adapter.model_pfun_raw = &tw_p_raw;
  if (adapter.logS_at_t_ptr  != nullptr) adapter.logS_at_t_ptr  = &tw_logS_at_t;
}

// An invalid clock parameter (only possible for the RDMSWTN clocks; see
// clock_is_valid) makes the row impossible: zero density and zero survivor,
// so the trial takes the likelihood floor however the row is used.

double tw_pdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double t0 = par[ctx->t0_index];
  const double p  = par[tw.par_index];
  const double x  = t - t0;
  if (clock_is_identity(tw.clock, p) || !(x > 0.0) || !R_FINITE(x))
    return tw.base_pdf1(t, par, ctx_);        // includes t == +Inf
  if (!clock_is_valid(tw.clock, p)) return 0.0;
  const double d = tw.base_pdf1(t0 + clock_fwd(tw.clock, x, p), par, ctx_);
  if (!(d > 0.0) || !R_FINITE(d)) return 0.0;
  const double log_density = std::log(d) + clock_log_jac(tw.clock, x, p);
  return R_FINITE(log_density) ? std::exp(log_density) : 0.0;
}

double tw_cdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double t0 = par[ctx->t0_index];
  const double p  = par[tw.par_index];
  const double x  = t - t0;
  if (clock_is_identity(tw.clock, p) || !(x > 0.0))
    return tw.base_cdf1(t, par, ctx_);
  if (!clock_is_valid(tw.clock, p)) return 1.0;
  // At t = +Inf the unbounded clocks map to +Inf (the defect reaches the
  // base); the exhaustion clock maps to its frozen budget.
  return tw.base_cdf1(t0 + clock_fwd(tw.clock, x, p), par, ctx_);
}

void tw_d_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok, double* out,
              double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const int clock = tw.clock;
  const double* t0_ = cols[ctx->t0_index];
  const double* p_  = cols[tw.par_index];

  bool any = false;
  if (p_ != nullptr && t0_ != nullptr) {
    for (int i = 0; i < n_rows && !any; ++i) {
      if (mask[i] && !clock_is_identity(clock, p_[i])) any = true;
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
    const double p = p_[i];
    if (clock_is_identity(clock, p)) continue;
    const double x = rt[i] - t0_[i];
    if (!(x > 0.0) || !R_FINITE(x)) continue;  // guards match base clock
    if (!clock_is_valid(clock, p)) {
      lj_buf[i] = R_NegInf;
      continue;
    }
    rt_buf[i] = t0_[i] + clock_fwd(clock, x, p);
    lj_buf[i] = clock_log_jac(clock, x, p);
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
  const int clock = tw.clock;
  const double* t0_ = cols[ctx->t0_index];
  const double* p_  = cols[tw.par_index];

  bool any = false;
  if (p_ != nullptr && t0_ != nullptr) {
    for (int i = 0; i < n_rows && !any; ++i) {
      if (mask[i] && !clock_is_identity(clock, p_[i])) any = true;
    }
  }
  if (!any) {                                  // exact parent path
    tw.base_p_raw(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
    return;
  }

  static thread_local std::vector<double> rt_buf;
  rt_buf.assign(rt, rt + n_rows);
  bool any_invalid = false;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double p = p_[i];
    if (clock_is_identity(clock, p)) continue;
    const double x = rt[i] - t0_[i];
    // x = +Inf is kept: the unbounded clocks return +Inf, and the exhaustion
    // clock freezes at its budget.
    if (!(x > 0.0)) continue;
    if (!clock_is_valid(clock, p)) {
      any_invalid = true;
      continue;
    }
    rt_buf[i] = t0_[i] + clock_fwd(clock, x, p);
  }
  tw.base_p_raw(rt_buf.data(), cols, n_rows, mask, isok, out, min_ll, ctx_);
  if (!any_invalid) return;
  const bool floor_raw = raw_floor_log_lik(ctx_);
  for (int i = 0; i < n_rows; ++i) {
    if (mask[i] && !clock_is_valid(clock, p_[i]) && rt[i] - t0_[i] > 0.0)
      out[i] = raw_log_zero(min_ll, floor_raw);
  }
}

// log S at a common t on a clock with a finite budget, t = +Inf: the rows'
// operational times differ (tau_r/2), so the common-t kernel cannot take them.
// Evaluate each row's survivor through the raw kernel instead.
static void tw_logS_at_inf_rowwise(const double* const* cols,
                                   int n_rows_total, int n_lR,
                                   const int* trunc_mask, int n_unique_trials,
                                   const int* isok_all, void* ctx_,
                                   double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const double* t0_ = cols[ctx->t0_index];
  const double* p_  = cols[tw.par_index];

  static thread_local std::vector<double> rt_buf;
  static thread_local std::vector<double> out_buf;
  static thread_local std::vector<int> mask_buf;
  rt_buf.assign(static_cast<size_t>(n_rows_total), R_PosInf);
  out_buf.assign(static_cast<size_t>(n_rows_total), 0.0);
  mask_buf.assign(static_cast<size_t>(n_rows_total), 0);
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    for (int k = 0; k < n_lR; ++k) {
      const int r = j * n_lR + k;
      mask_buf[r] = 1;
      const double p = p_[r];
      if (clock_is_valid(tw.clock, p))
        rt_buf[r] = t0_[r] + clock_fwd(tw.clock, R_PosInf, p);
    }
  }
  const bool floor_raw = ctx->floor_raw_log_lik;
  ctx->floor_raw_log_lik = false;
  tw.base_p_raw(rt_buf.data(), cols, n_rows_total, mask_buf.data(), isok_all,
                out_buf.data(), R_NegInf, ctx_);
  ctx->floor_raw_log_lik = floor_raw;

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double logS = 0.0;
    for (int k = 0; k < n_lR; ++k) {
      const int r = j * n_lR + k;
      const double v = out_buf[r];
      if (!isok_all[r] || !clock_is_valid(tw.clock, p_[r]) || ISNAN(v) ||
          v == R_NegInf) {
        logS = R_NegInf;
        break;
      }
      logS += std::fmin(v, 0.0);
    }
    logS_out[j] = logS;
  }
}

void tw_logS_at_t(double t, const double* const* cols, int n_rows_total,
                  int n_lR, int n_par, const int* trunc_mask,
                  int n_unique_trials, const int* isok_all, void* ctx_,
                  double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const TimeWarpPlan& tw = ctx->tw;
  const int clock = tw.clock;
  const double* t0_ = cols[ctx->t0_index];
  const double* p_  = cols[tw.par_index];

  bool any = false;
  if (p_ != nullptr && t0_ != nullptr) {
    for (int r = 0; r < n_rows_total && !any; ++r) {
      if (!clock_is_identity(clock, p_[r])) any = true;
    }
  }
  if (any && !R_FINITE(t)) {
    // The unbounded clocks send t = +Inf to +Inf, which is the parent's own
    // survivor at infinity.  A finite budget needs the row-wise evaluation.
    if (clock != emc2tw::CLOCK_EXHAUSTION) any = false;
    else {
      tw_logS_at_inf_rowwise(cols, n_rows_total, n_lR, trunc_mask,
                             n_unique_trials, isok_all, ctx_, logS_out);
      return;
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
  bool any_invalid = false;
  for (int r = 0; r < n_rows_total; ++r) {
    const double p = p_[r];
    if (clock_is_identity(clock, p)) continue;
    const double x = t - t0_[r];
    if (!(x > 0.0)) continue;                  // preserve not-started branch
    if (!clock_is_valid(clock, p)) {
      any_invalid = true;
      continue;
    }
    t0_buf[r] = t - clock_fwd(clock, x, p);    // t - t0_buf[r] == q(x)
  }

  // Kernels may fetch optional trailing pointers for this model variant.  Do
  // not read beyond the compacted n_par source columns while padding them.
  const int n_slots = std::max(n_par, 16);
  cols_buf.assign(static_cast<size_t>(n_slots), nullptr);
  for (int j = 0; j < n_par; ++j) cols_buf[j] = cols[j];
  cols_buf[ctx->t0_index] = t0_buf.data();

  tw.base_logS(t, cols_buf.data(), n_rows_total, n_lR, n_par, trunc_mask,
               n_unique_trials, isok_all, ctx_, logS_out);
  if (!any_invalid) return;
  for (int r = 0; r < n_rows_total; ++r) {
    const int j = r / n_lR;
    if (j < n_unique_trials && trunc_mask[j] && t - t0_[r] > 0.0 &&
        !clock_is_valid(clock, p_[r]))
      logS_out[j] = R_NegInf;
  }
}

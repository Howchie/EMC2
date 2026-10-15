#include "state_mixture.h"

#include <algorithm>
#include <cmath>
#include <vector>

// Natural-scale weight of the alternative state.  Weights at or below 1e-12
// (including NA) leave the accumulator in its base state without paying for
// the second kernel.
static inline double mix_weight(double p) {
  return (p > 1e-12) ? std::fmin(p, 1.0) : 0.0;
}

static inline double mix_log_combine(double log_base, double log_alt, double p) {
  const double a = std::log1p(-p) + log_base;
  const double b = std::log(p) + log_alt;
  const double m = std::fmax(a, b);
  if (m == R_NegInf) return R_NegInf;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

static inline double mix_alt_value(double base, double twin, int sign) {
  return sign == 0 ? twin : base + sign * twin;
}

// Copy one scalar parameter row with every swapped column in its alternative
// state.
static const double* mix_alt_row(const double* par, const MixturePlan& mx) {
  static thread_local std::vector<double> buf;
  buf.assign(par, par + mx.n_par);
  for (const auto& s : mx.swaps)
    buf[s.col] = mix_alt_value(par[s.col], par[s.twin], s.sign);
  return buf.data();
}

// Column-pointer array with every swapped column in its alternative state.  A
// replacing twin is used in place; a gap twin is applied into a scratch
// column.  The array is padded the same way the race paths pad theirs, so
// kernels may fetch (not dereference) optional trailing slots.
static const double* const* mix_alt_cols(const double* const* cols, int n_rows,
                                         const MixturePlan& mx) {
  static thread_local std::vector<const double*> buf;
  static thread_local std::vector<std::vector<double>> scratch;
  const int n_slots = std::max(mx.n_par, 16);
  buf.assign(static_cast<size_t>(n_slots), nullptr);
  for (int j = 0; j < mx.n_par; ++j) buf[j] = cols[j];
  if (scratch.size() < mx.swaps.size()) scratch.resize(mx.swaps.size());
  for (size_t k = 0; k < mx.swaps.size(); ++k) {
    const auto& s = mx.swaps[k];
    if (s.sign == 0) { buf[s.col] = cols[s.twin]; continue; }
    std::vector<double>& col = scratch[k];
    col.resize(static_cast<size_t>(n_rows));
    const double* b = cols[s.col];
    const double* t = cols[s.twin];
    for (int r = 0; r < n_rows; ++r) col[r] = b[r] + s.sign * t[r];
    buf[s.col] = col.data();
  }
  return buf.data();
}

template <typename F>
static double mix_scalar(double t, const double* par, void* ctx_, F base) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const MixturePlan& mx = ctx->mix;
  const double p = mix_weight(par[mx.p_index]);
  if (p == 0.0) return base(t, par, ctx_);
  const int hint = ctx->mode_hint;
  ctx->mode_hint = 0;                          // alternative row may change branch
  const double alt = base(t, mix_alt_row(par, mx), ctx_);
  ctx->mode_hint = hint;
  if (p == 1.0) return alt;
  return (1.0 - p) * base(t, par, ctx_) + p * alt;
}

double mix_pdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return mix_scalar(t, par, ctx_, ctx->mix.base_pdf1);
}

double mix_cdf1(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return mix_scalar(t, par, ctx_, ctx->mix.base_cdf1);
}

// Shared body of the raw density and survivor wrappers: both outputs are
// logs of quantities that mix linearly in the natural scale.
static void mix_raw(RaceRawFun base, const double* rt, const double* const* cols,
                    int n_rows, const int* mask, const int* isok, double* out,
                    double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const MixturePlan& mx = ctx->mix;
  const double* p_ = cols[mx.p_index];

  bool any = false;
  for (int i = 0; i < n_rows && !any; ++i)
    if (mask[i] && mix_weight(p_[i]) > 0.0) any = true;
  if (!any) {                                  // exact parent path
    base(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
    return;
  }

  static thread_local std::vector<int> mask0, mask1;
  static thread_local std::vector<double> alt;
  mask0.assign(static_cast<size_t>(n_rows), 0);
  mask1.assign(static_cast<size_t>(n_rows), 0);
  alt.assign(static_cast<size_t>(n_rows), R_NegInf);
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double p = mix_weight(p_[i]);
    mask0[i] = p < 1.0;
    mask1[i] = p > 0.0;
  }

  const bool floor_raw = ctx->floor_raw_log_lik;
  const int hint = ctx->mode_hint;
  ctx->floor_raw_log_lik = false;              // mix before flooring
  base(rt, cols, n_rows, mask0.data(), isok, out, min_ll, ctx_);
  ctx->mode_hint = 0;
  base(rt, mix_alt_cols(cols, n_rows, mx), n_rows, mask1.data(), isok, alt.data(),
       min_ll, ctx_);
  ctx->mode_hint = hint;
  ctx->floor_raw_log_lik = floor_raw;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const double p = mix_weight(p_[i]);
    double v;
    if (p == 0.0) v = out[i];
    else if (p == 1.0) v = alt[i];
    else v = mix_log_combine(out[i], alt[i], p);
    out[i] = raw_log_value(v, min_ll, floor_raw);
  }
}

void mix_d_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok, double* out,
               double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  mix_raw(ctx->mix.base_d_raw, rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}

void mix_p_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok, double* out,
               double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  mix_raw(ctx->mix.base_p_raw, rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}

// The batch survivor is a per-trial product over accumulators, which does not
// factor through a mixture of whole rows; when any included row mixes, each
// accumulator's survivor is evaluated through the mixed raw survivor at t.
void mix_logS_at_t(double t, const double* const* cols, int n_rows_total,
                   int n_lR, int n_par, const int* trunc_mask,
                   int n_unique_trials, const int* isok_all, void* ctx_,
                   double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const MixturePlan& mx = ctx->mix;
  const double* p_ = cols[mx.p_index];

  bool any = false;
  for (int j = 0; j < n_unique_trials && !any; ++j) {
    if (!trunc_mask[j]) continue;
    for (int k = 0; k < n_lR; ++k)
      if (mix_weight(p_[j * n_lR + k]) > 0.0) { any = true; break; }
  }
  if (!any || mx.base_p_raw == nullptr) {
    mx.base_logS(t, cols, n_rows_total, n_lR, n_par, trunc_mask,
                 n_unique_trials, isok_all, ctx_, logS_out);
    return;
  }

  static thread_local std::vector<double> rt_buf;
  static thread_local std::vector<double> out_buf;
  static thread_local std::vector<int> mask_buf;
  rt_buf.assign(static_cast<size_t>(n_rows_total), t);
  out_buf.assign(static_cast<size_t>(n_rows_total), 0.0);
  mask_buf.assign(static_cast<size_t>(n_rows_total), 0);
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    for (int k = 0; k < n_lR; ++k) mask_buf[j * n_lR + k] = 1;
  }
  const bool floor_raw = ctx->floor_raw_log_lik;
  ctx->floor_raw_log_lik = false;
  mix_p_raw(rt_buf.data(), cols, n_rows_total, mask_buf.data(), isok_all,
            out_buf.data(), R_NegInf, ctx_);
  ctx->floor_raw_log_lik = floor_raw;

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double logS = 0.0;
    for (int k = 0; k < n_lR; ++k) {
      const int r = j * n_lR + k;
      const double v = out_buf[r];
      if (!isok_all[r] || ISNAN(v) || v == R_NegInf) {
        logS = R_NegInf;
        break;
      }
      logS += std::fmin(v, 0.0);
    }
    logS_out[j] = logS;
  }
}

void configure_state_mixture_context(RaceModelAdapter& adapter,
                                     const Rcpp::CharacterVector& keep_names,
                                     bool is_logicalrules,
                                     const std::string& caller) {
  MixturePlan& mx = adapter.ctx.mix;
  const int n = keep_names.size();
  std::vector<std::string> names(static_cast<size_t>(n));
  for (int j = 0; j < n; ++j) names[j] = Rcpp::as<std::string>(keep_names[j]);

  int p_index = -1;
  std::vector<MixturePlan::Swap> swaps;
  static const std::pair<const char*, int> suffixes[] = {
    {"_mix", 0}, {"_below", -1}, {"_above", 1}};
  for (int j = 0; j < n; ++j) {
    const std::string& nm = names[j];
    if (nm == "pmix") { p_index = j; continue; }
    int sign = 2;
    std::string suffix;
    for (const auto& sx : suffixes) {
      const std::string sf = sx.first;
      if (nm.size() > sf.size() &&
          nm.compare(nm.size() - sf.size(), sf.size(), sf) == 0) {
        sign = sx.second;
        suffix = sf;
        break;
      }
    }
    if (sign == 2) continue;
    const std::string target = nm.substr(0, nm.size() - suffix.size());
    const auto it = std::find(names.begin(), names.end(), target);
    if (it == names.end()) {
      Rcpp::stop("%s: mixture column '%s' has no base column '%s'.",
                 caller.c_str(), nm.c_str(), target.c_str());
    }
    if (target == "pmix" || target == "pContaminant" || target == "pGuess") {
      Rcpp::stop("%s: '%s' cannot have a mixture state.", caller.c_str(),
                 target.c_str());
    }
    const int col = static_cast<int>(it - names.begin());
    for (const auto& sw : swaps) {
      if (sw.col == col)
        Rcpp::stop("%s: '%s' has more than one mixture twin.", caller.c_str(),
                   target.c_str());
    }
    swaps.push_back({col, j, sign});
  }
  if (p_index < 0 && swaps.empty()) return;
  if (p_index < 0 || swaps.empty()) {
    Rcpp::stop("%s: a state mixture needs both a 'pmix' column and at least one "
               "'<par>_mix', '<par>_below' or '<par>_above' column.", caller.c_str());
  }
  if (is_logicalrules || adapter.ctx.corr_drift_active ||
      adapter.ctx.rdmswtn_correlated) {
    Rcpp::stop("%s: the state mixture needs independent accumulators; it is not "
               "available for LogicalRules or correlated races.", caller.c_str());
  }

  mx.p_index    = p_index;
  mx.n_par      = n;
  mx.swaps      = swaps;
  mx.base_pdf1  = adapter.pdf1_ptr;
  mx.base_cdf1  = adapter.cdf1_ptr;
  mx.base_d_raw = adapter.model_dfun_raw;
  mx.base_p_raw = adapter.model_pfun_raw;
  mx.base_logS  = adapter.logS_at_t_ptr;

  if (adapter.pdf1_ptr       != nullptr) adapter.pdf1_ptr       = &mix_pdf1;
  if (adapter.cdf1_ptr       != nullptr) adapter.cdf1_ptr       = &mix_cdf1;
  if (adapter.model_dfun_raw != nullptr) adapter.model_dfun_raw = &mix_d_raw;
  if (adapter.model_pfun_raw != nullptr) adapter.model_pfun_raw = &mix_p_raw;
  if (adapter.logS_at_t_ptr  != nullptr) adapter.logS_at_t_ptr  = &mix_logS_at_t;
}

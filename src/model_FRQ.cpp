// R-callable entry points for FRQ.
//
// These live in their own translation unit rather than in model_FRQ.h because
// the header is also included by model_rng.cpp for the simulator, and a header
// carrying [[Rcpp::export]] definitions can only be seen by one TU.  They
// bypass ContextForRaceModels entirely and so take the parameters directly.

#include "utility_functions.h"
#include "model_FRQ.h"
#include "race_contract.h"

// ---------------------------------------------------------------------------
// R-callable entry points.  All take the DECISION time (rt - t0); R/model_FRQ.R
// subtracts t0 before calling, exactly as it does for the other race models.
// Arguments recycle against `t`, as in dbawd/pbawd: each must be length 1 or
// the full length, and anything else is an error (see frq_check_len).
// ---------------------------------------------------------------------------

// Only scalars and full-length vectors recycle.  A partial-length argument is
// an error rather than R-style wrap-around: the loops index directly, so
// silently accepting one would read past the end of the vector.
static inline void frq_check_len(const NumericVector& x, int n, const char* nm) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("FRQ: `%s` must be length 1 or length %d, not %d.",
               nm, n, static_cast<int>(x.size()));
}

static inline double frq_pick(const NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

static inline FrqCache* frq_cache_from_context(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr && ctx->frq_cache
    ? ctx->frq_cache.get() : nullptr;
}

// [[Rcpp::export]]
NumericVector dfrq(NumericVector t, NumericVector K, NumericVector v,
                   NumericVector sv, NumericVector h,
                   NumericVector plo = NumericVector::create(0.0),
                   bool log_out = false,
                   NumericVector clock = NumericVector::create(1.0),
                   NumericVector accel = NumericVector::create(0.0),
                   NumericVector omega = NumericVector::create(0.0)) {
  const int n = t.size();
  frq_check_len(K, n, "K");         frq_check_len(v, n, "v");
  frq_check_len(sv, n, "sv");       frq_check_len(h, n, "h");
  frq_check_len(plo, n, "plo");
  frq_check_len(clock, n, "clock");
  frq_check_len(accel, n, "accel");
  frq_check_len(omega, n, "omega");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqCache cache;
  FrqMemo memo(&cache);
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(K, i), pick(v, i), pick(sv, i),
                                pick(h, i), pick(plo, i), pick(clock, i),
                                pick(omega, i));
    const double lp = frq_accel_log_pdf_dt(t[i], s, pick(accel, i));
    out[i] = log_out ? lp : ((lp > R_NegInf) ? std::exp(lp) : 0.0);
  }
  return out;
}

// `lower_tail = false` returns the SURVIVOR 1 - F, which tends to the
// defective mass 1 - h rather than to zero.  It is a separate branch rather
// than 1 - pfrq(...) because the whole point of the reflection identity is to
// avoid that subtraction.
// [[Rcpp::export]]
NumericVector pfrq(NumericVector t, NumericVector K, NumericVector v,
                   NumericVector sv, NumericVector h,
                   NumericVector plo = NumericVector::create(0.0),
                   bool lower_tail = true, bool log_out = false,
                   NumericVector clock = NumericVector::create(1.0),
                   NumericVector accel = NumericVector::create(0.0),
                   NumericVector omega = NumericVector::create(0.0)) {
  const int n = t.size();
  frq_check_len(K, n, "K");         frq_check_len(v, n, "v");
  frq_check_len(sv, n, "sv");       frq_check_len(h, n, "h");
  frq_check_len(plo, n, "plo");
  frq_check_len(clock, n, "clock");
  frq_check_len(accel, n, "accel");
  frq_check_len(omega, n, "omega");
  NumericVector out(n);
  auto pick = frq_pick;
  FrqCache cache;
  FrqMemo memo(&cache);
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(K, i), pick(v, i), pick(sv, i),
                                pick(h, i), pick(plo, i), pick(clock, i),
                                pick(omega, i));
    const double a = pick(accel, i);
    double lp;
    if (lower_tail) {
      lp = frq_accel_log_cdf_dt(t[i], s, a);
    } else {
      lp = frq_accel_log_surv_dt(t[i], s, a);
      if (ISNAN(lp)) lp = R_NaN;
    }
    if (ISNAN(lp)) { out[i] = R_NaN; continue; }
    out[i] = log_out ? lp : ((lp > R_NegInf) ? std::exp(lp) : 0.0);
  }
  return out;
}

// Map the fitted coordinates to the generative ones: the Beta shape
// beta = (v/sv)^2, the availability p, the registration rate lambda and the
// expected total evidence E (Inf when h = 1), and log(1 - p).  Used by
// Ttransform and the R reference simulator.
// [[Rcpp::export]]
NumericMatrix frq_rate(NumericVector K, NumericVector v, NumericVector sv,
                       NumericVector h,
                       NumericVector plo = NumericVector::create(0.0),
                       NumericVector clock = NumericVector::create(1.0),
                       NumericVector omega = NumericVector::create(0.0)) {
  // Length is the longest argument, not K's: a scalar K against a vector v has
  // to give one row per v.
  const int n = static_cast<int>(std::max({K.size(), v.size(), sv.size(),
                                           h.size(), plo.size(), clock.size(),
                                           omega.size()}));
  frq_check_len(K, n, "K");         frq_check_len(v, n, "v");
  frq_check_len(sv, n, "sv");       frq_check_len(h, n, "h");
  frq_check_len(plo, n, "plo");
  frq_check_len(clock, n, "clock");
  frq_check_len(omega, n, "omega");
  NumericMatrix out(n, 5);
  auto pick = frq_pick;
  FrqMemo memo;
  for (int i = 0; i < n; ++i) {
    const FrqPars& s = memo.get(pick(K, i), pick(v, i), pick(sv, i),
                                pick(h, i), pick(plo, i), pick(clock, i),
                                pick(omega, i));
    out(i, 0) = s.ok ? s.beta : NA_REAL;
    out(i, 1) = s.ok ? s.p : NA_REAL;
    out(i, 2) = s.ok ? ((s.poisson && s.h == 1.0) ? NA_REAL : s.lambda) : NA_REAL;
    out(i, 3) = s.ok ? s.E : NA_REAL;
    out(i, 4) = s.ok ? s.log_om_p : NA_REAL;
  }
  colnames(out) = CharacterVector::create("beta", "p", "lambda", "Etot",
                                          "log1mp");
  return out;
}

// ============================================================
// FRQ (Finite Reservoir Quorum) adapters
// Column layout: K=0, v=1, sv=2, h=3, t0=4, plo=5, clock=6; the cv chart
// (FRQ_CV) holds the rate CV sv / v in column 2, the gamma chart (FRQ_GAMMA)
// gamma = sv^2 / v in column 2, and the accel chart (FRQ_ACCEL) the rate CV in
// column 2 and the clock acceleration in column 6.  The power chart (FRQ_POWER)
// holds cv1 in column 2 and its exponent cvpow in column 7, with
// cv = cv1 v^(-cvpow).
// plo is the probability of the low-threshold state (zero by default) and
// clock the registration clock shape (one by default, a constant-rate counter).
//
// The upper tail may be defective: an accumulator terminates with probability
// h = I_p(alpha, beta), and the leftover mass 1 - h sits at t = +Inf.  At
// h = 1 the survivor reaches zero and the same paths become proper.
//
// Every entry point resolves its row through FrqMemo, which splits the state
// on (K, sv/v, h, plo) -- the beta-quantile half, cached per particle --
// and (v, clock), which is recomputed per row for the cost of one log.
// ============================================================

// Charts sharing the column layout: column 2 holds sv (FRQ_SV), the rate CV
// (FRQ_CV_CHART, FRQ_ACCEL_CHART), gamma = sv^2 / v (FRQ_GAMMA_CHART) or cv1
// (FRQ_POWER_CHART), and column 6 holds the clock shape, or the clock
// acceleration under the accel chart (whose clock shape is one).  The power
// chart's exponent cvpow (column 7) reaches frq_memo_row as `pw`; the other
// charts pass zero.
enum FrqChart : int { FRQ_SV = 0, FRQ_CV_CHART = 1, FRQ_ACCEL_CHART = 2,
                      FRQ_GAMMA_CHART = 3, FRQ_POWER_CHART = 4 };

template <int CH>
inline const FrqPars& frq_memo_row(FrqMemo& memo, double k, double v,
                                   double x, double h, double pl, double cu,
                                   double om, double pw = 0.0) {
  if (CH == FRQ_SV) return memo.get(k, v, x, h, pl, cu, om);
  if (CH == FRQ_GAMMA_CHART)
    return memo.get_cv(k, v, frq_gamma_cv(x, v), h, pl, cu, om);
  if (CH == FRQ_POWER_CHART)
    return memo.get_cv(k, v, frq_power_cv(x, v, pw), h, pl, cu, om);
  return memo.get_cv(k, v, x, h, pl, CH == FRQ_ACCEL_CHART ? 1.0 : cu, om);
}

// The excess-criterion column is resolved BY NAME (configure_frq_omega_context
// in race_dispatch.cpp), so a design without one reads omega = 0.
static inline int frq_omega_col(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  return ctx != nullptr ? ctx->frq_omega_index : -1;
}

template <int CH>
inline double frq_chart_log_pdf(double u, const FrqPars& s, double cu) {
  return CH == FRQ_ACCEL_CHART ? frq_accel_log_pdf_dt(u, s, cu)
                               : frq_log_pdf_dt(u, s);
}

template <int CH>
inline double frq_chart_log_surv(double u, const FrqPars& s, double cu) {
  return CH == FRQ_ACCEL_CHART ? frq_accel_log_surv_dt(u, s, cu)
                               : frq_log_surv_dt(u, s);
}

template <int CH>
inline double frq_chart_log_cdf(double u, const FrqPars& s, double cu) {
  return CH == FRQ_ACCEL_CHART ? frq_accel_log_cdf_dt(u, s, cu)
                               : frq_log_cdf_dt(u, s);
}

// The scalar entry points are called once per evaluation point, often several
// times on one row (bin edges and midpoint), so they share a per-thread memo.
// Its keys are the parameter values themselves, so a hit is valid across
// particles; only a miss consults the per-particle cache.
template <int CH>
static inline const FrqPars& frq_scalar_row(const double* par, void* ctx_) {
  static thread_local FrqMemo memo;
  memo.shared = frq_cache_from_context(ctx_);
  const int io = frq_omega_col(ctx_);
  return frq_memo_row<CH>(memo, par[emc2col::frq::K], par[emc2col::frq::v],
                          par[emc2col::frq::sv], par[emc2col::frq::h],
                          par[emc2col::frq::plo], par[emc2col::frq::clock],
                          io >= 0 ? par[io] : 0.0,
                          CH == FRQ_POWER_CHART ? par[emc2col::frq_power::cvpow] : 0.0);
}

template <int CH>
static double dfrq_scalar_impl(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq::K]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  const double lp = frq_chart_log_pdf<CH>(tt, frq_scalar_row<CH>(par, ctx_),
                                          par[emc2col::frq::clock]);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

template <int CH>
static double pfrq_scalar_impl(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::frq::K]) ||
      !R_FINITE(par[emc2col::frq::t0])) return 0.0;
  const double tt = t - par[emc2col::frq::t0];
  if (t <= 0.0 || tt <= 0.0) return 0.0;
  // tt == Inf deliberately reaches the kernel: the CDF there is h, not one
  // when h < 1, and is one at the proper h = 1 boundary.
  const double lp = frq_chart_log_cdf<CH>(tt, frq_scalar_row<CH>(par, ctx_),
                                          par[emc2col::frq::clock]);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

template <int CH>
static void dfrq_raw_impl(const double* rt, const double* const* cols, int n_rows,
                          const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* k_  = cols[emc2col::frq::K];
  const double* v_  = cols[emc2col::frq::v];
  const double* sv_ = cols[emc2col::frq::sv];
  const double* h_  = cols[emc2col::frq::h];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* pl_ = cols[emc2col::frq::plo];
  const double* cu_ = cols[emc2col::frq::clock];
  const double* pw_ = CH == FRQ_POWER_CHART ? cols[emc2col::frq_power::cvpow] : nullptr;
  const int io = frq_omega_col(ctx_);
  const double* om_ = io >= 0 ? cols[io] : nullptr;
  FrqMemo memo(frq_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(k_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    const FrqPars& s = frq_memo_row<CH>(memo, k_[i], v_[i], sv_[i], h_[i],
                                        pl_[i], cu_[i], om_ ? om_[i] : 0.0,
                                        pw_ ? pw_[i] : 0.0);
    const double log_pdf = frq_chart_log_pdf<CH>(tt, s, cu_[i]);
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw)
      : raw_log_zero(min_ll, floor_raw);
  }
}

template <int CH>
static void pfrq_raw_impl(const double* rt, const double* const* cols, int n_rows,
                          const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* k_  = cols[emc2col::frq::K];
  const double* v_  = cols[emc2col::frq::v];
  const double* sv_ = cols[emc2col::frq::sv];
  const double* h_  = cols[emc2col::frq::h];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* pl_ = cols[emc2col::frq::plo];
  const double* cu_ = cols[emc2col::frq::clock];
  const double* pw_ = CH == FRQ_POWER_CHART ? cols[emc2col::frq_power::cvpow] : nullptr;
  const int io = frq_omega_col(ctx_);
  const double* om_ = io >= 0 ? cols[io] : nullptr;
  FrqMemo memo(frq_cache_from_context(ctx_));
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    // A loser that cannot be evaluated contributes a survivor of one, matching
    // every other race adapter: the trial is failed by its winner's density.
    if (R_IsNA(k_[i]) || !isok[i] || !R_FINITE(t0_[i])) {
      out[i] = 0.0; continue;
    }
    const double tt = rt[i] - t0_[i];
    if (tt <= 0.0 || rt[i] <= 0.0) { out[i] = 0.0; continue; }
    const FrqPars& s = frq_memo_row<CH>(memo, k_[i], v_[i], sv_[i], h_[i],
                                        pl_[i], cu_[i], om_ ? om_[i] : 0.0,
                                        pw_ ? pw_[i] : 0.0);
    const double log_surv = frq_chart_log_surv<CH>(tt, s, cu_[i]);
    if (ISNAN(log_surv)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = (log_surv > R_NegInf) ? log_surv
                                   : raw_log_zero(min_ll, floor_raw);
  }
}

template <int CH>
static void frq_logS_at_t_impl(double t, const double* const* cols, int n_lR,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_,
                               double* logS_out) {
  const double* k_  = cols[emc2col::frq::K];
  const double* v_  = cols[emc2col::frq::v];
  const double* sv_ = cols[emc2col::frq::sv];
  const double* h_  = cols[emc2col::frq::h];
  const double* t0_ = cols[emc2col::frq::t0];
  const double* pl_ = cols[emc2col::frq::plo];
  const double* cu_ = cols[emc2col::frq::clock];
  const double* pw_ = CH == FRQ_POWER_CHART ? cols[emc2col::frq_power::cvpow] : nullptr;
  const int io = frq_omega_col(ctx_);
  const double* om_ = io >= 0 ? cols[io] : nullptr;
  FrqMemo memo(frq_cache_from_context(ctx_));
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(k_[r]) || !R_FINITE(t0_[r])) {
        bad = true; break;
      }
      const double tt = t - t0_[r];
      if (tt <= 0.0) continue;  // not started: survivor one
      const FrqPars& s = frq_memo_row<CH>(memo, k_[r], v_[r], sv_[r], h_[r],
                                          pl_[r], cu_[r], om_ ? om_[r] : 0.0,
                                        pw_ ? pw_[r] : 0.0);
      const double log_surv = !R_FINITE(t)
        ? (s.ok ? frq_chart_log_surv<CH>(R_PosInf, s, cu_[r]) : R_NaN)
        : frq_chart_log_surv<CH>(tt, s, cu_[r]);
      if (ISNAN(log_surv)) { bad = true; break; }
      logS += log_surv;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

double dfrq_scalar(double t, const double* par, void* ctx_) {
  return dfrq_scalar_impl<FRQ_SV>(t, par, ctx_);
}
double pfrq_scalar(double t, const double* par, void* ctx_) {
  return pfrq_scalar_impl<FRQ_SV>(t, par, ctx_);
}
void dfrq_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  dfrq_raw_impl<FRQ_SV>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void pfrq_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_) {
  pfrq_raw_impl<FRQ_SV>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void frq_logS_at_t(double t, const double* const* cols,
                   int /*n_rows_total*/, int n_lR, int /*n_par*/,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* ctx_, double* logS_out) {
  frq_logS_at_t_impl<FRQ_SV>(t, cols, n_lR, trunc_mask, n_unique_trials,
                            isok_all, ctx_, logS_out);
}

double dfrq_cv_scalar(double t, const double* par, void* ctx_) {
  return dfrq_scalar_impl<FRQ_CV_CHART>(t, par, ctx_);
}
double pfrq_cv_scalar(double t, const double* par, void* ctx_) {
  return pfrq_scalar_impl<FRQ_CV_CHART>(t, par, ctx_);
}
void dfrq_cv_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_) {
  dfrq_raw_impl<FRQ_CV_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void pfrq_cv_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_) {
  pfrq_raw_impl<FRQ_CV_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void frq_cv_logS_at_t(double t, const double* const* cols,
                      int /*n_rows_total*/, int n_lR, int /*n_par*/,
                      const int* trunc_mask, int n_unique_trials,
                      const int* isok_all, void* ctx_, double* logS_out) {
  frq_logS_at_t_impl<FRQ_CV_CHART>(t, cols, n_lR, trunc_mask, n_unique_trials,
                           isok_all, ctx_, logS_out);
}

double dfrq_gamma_scalar(double t, const double* par, void* ctx_) {
  return dfrq_scalar_impl<FRQ_GAMMA_CHART>(t, par, ctx_);
}
double pfrq_gamma_scalar(double t, const double* par, void* ctx_) {
  return pfrq_scalar_impl<FRQ_GAMMA_CHART>(t, par, ctx_);
}
void dfrq_gamma_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_) {
  dfrq_raw_impl<FRQ_GAMMA_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void pfrq_gamma_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_) {
  pfrq_raw_impl<FRQ_GAMMA_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void frq_gamma_logS_at_t(double t, const double* const* cols,
                         int /*n_rows_total*/, int n_lR, int /*n_par*/,
                         const int* trunc_mask, int n_unique_trials,
                         const int* isok_all, void* ctx_, double* logS_out) {
  frq_logS_at_t_impl<FRQ_GAMMA_CHART>(t, cols, n_lR, trunc_mask, n_unique_trials,
                                      isok_all, ctx_, logS_out);
}

double dfrq_power_scalar(double t, const double* par, void* ctx_) {
  return dfrq_scalar_impl<FRQ_POWER_CHART>(t, par, ctx_);
}
double pfrq_power_scalar(double t, const double* par, void* ctx_) {
  return pfrq_scalar_impl<FRQ_POWER_CHART>(t, par, ctx_);
}
void dfrq_power_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_) {
  dfrq_raw_impl<FRQ_POWER_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void pfrq_power_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_) {
  pfrq_raw_impl<FRQ_POWER_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void frq_power_logS_at_t(double t, const double* const* cols,
                         int /*n_rows_total*/, int n_lR, int /*n_par*/,
                         const int* trunc_mask, int n_unique_trials,
                         const int* isok_all, void* ctx_, double* logS_out) {
  frq_logS_at_t_impl<FRQ_POWER_CHART>(t, cols, n_lR, trunc_mask, n_unique_trials,
                                      isok_all, ctx_, logS_out);
}

double dfrq_accel_scalar(double t, const double* par, void* ctx_) {
  return dfrq_scalar_impl<FRQ_ACCEL_CHART>(t, par, ctx_);
}
double pfrq_accel_scalar(double t, const double* par, void* ctx_) {
  return pfrq_scalar_impl<FRQ_ACCEL_CHART>(t, par, ctx_);
}
void dfrq_accel_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_) {
  dfrq_raw_impl<FRQ_ACCEL_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void pfrq_accel_raw(const double* rt, const double* const* cols, int n_rows,
                    const int* mask, const int* isok,
                    double* out, double min_ll, void* ctx_) {
  pfrq_raw_impl<FRQ_ACCEL_CHART>(rt, cols, n_rows, mask, isok, out, min_ll, ctx_);
}
void frq_accel_logS_at_t(double t, const double* const* cols,
                         int /*n_rows_total*/, int n_lR, int /*n_par*/,
                         const int* trunc_mask, int n_unique_trials,
                         const int* isok_all, void* ctx_, double* logS_out) {
  frq_logS_at_t_impl<FRQ_ACCEL_CHART>(t, cols, n_lR, trunc_mask, n_unique_trials,
                                      isok_all, ctx_, logS_out);
}

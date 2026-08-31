#include <Rcpp.h>
#include "model_BAwL.h"
#include "model_LBA.h"
#include "race_contract.h"
#include "timer_helpers.h"
#include "utility_functions.h"
#include "col_registry.h"

// ============================================================
// BAwL (Ballistic Accumulator with Leak + killing/guessing) adapters
// Column layout: v=0, sv=1, B=2, A=3, t0=4, k=5, mG=6 (guess-clock mean), mK=7 (kill-clock mean)
// ============================================================
// BAwL(parameterization = "ratio") samples r = k / (mean launch strength) in
// k's column, so the leak has to be rebuilt per row.  The clamp at zero is
// deliberate: the normal launch's v is identity-transformed and unbounded, and
// a negative leak is not a model, so a non-positive mean drift collapses
// continuously to the LBA limit (k = 0).  This is about the mean-drift
// PARAMETER, not the per-trial draws -- posdrift is untouched here, so an IO
// model with a positive mean keeps its full leak.  Mirrored exactly by the
// Ttransform of BAwL() in R/model_LBA.R.
static inline double bawl_leak(const ContextForRaceModels* ctx,
                               double k_or_r, double mean_launch) {
  if (ctx == nullptr || !ctx->bawl_ratio_chart) return k_or_r;
  return (mean_launch > 0.0) ? k_or_r * mean_launch : 0.0;
}

double dbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const bool split = lau == BAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = emc2col::select_index(split, emc2col::bawlsplit::mu, emc2col::bawl::v);
  const int isv = emc2col::select_index(split, emc2col::bawlsplit::sigma, emc2col::bawl::sv);
  const int iB = emc2col::select_index(split, emc2col::bawlsplit::B, emc2col::bawl::B);
  const int iA = emc2col::select_index(split, emc2col::bawlsplit::A, emc2col::bawl::A);
  const int it0 = emc2col::select_index(split, emc2col::bawlsplit::t0, emc2col::bawl::t0);
  const int ik = emc2col::select_index(split, emc2col::bawlsplit::k, emc2col::bawl::k);
  const int imG = emc2col::select_index(split, emc2col::bawlsplit::mG, emc2col::bawl::mG);
  const int imK = emc2col::select_index(split, emc2col::bawlsplit::mK, emc2col::bawl::mK);
  const double delta = split ? par[emc2col::bawlsplit::delta] : 0.0;
  // The ratio chart divides by the launch MEAN: column 0 for the normal (v)
  // and lognormal (mean) launches, column 1 for the Weibull pair (shape, mean).
  const int imean = (lau == BAWL_LAUNCH_WEIBULL) ? isv : iv;
  if (R_IsNA(par[iv])) return 0.0;
  const double t0_val = par[it0];
  const double tt = t - t0_val;
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    if (tt <= 0.0 || t <= 0.0) return 0.0;
    return bawl_pdf_scalar_natural(
      tt, par[iA], par[iB] + par[iA], par[iv], par[isv], 0.0,
      ctx->use_posdrift, LBA_DENOM_FLOOR, lau, delta);
  }
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const int ks = ctx ? ctx->kill_shape : 1;
  const double k_val = (ctx && ctx->bawl_k_fixed_zero) ? 0.0
                                                       : bawl_leak(ctx, par[ik], par[imean]);
  const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[imG], ks) : 0.0);
  const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active && ctx->apply_lk_to_racers) ?
      erlang_lambda_from_mean(par[imK], ks) : 0.0);
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  return dkilledleakyba_norm(
    t, par[iv], par[iB] + par[iA], par[iA], par[isv], t0_val, k_val, lg, lk,
    ctx->use_posdrift, false, ks, local_guess, omega, lau, delta);
}

double pbawl_scalar(double t, const double* par, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const bool split = lau == BAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = emc2col::select_index(split, emc2col::bawlsplit::mu, emc2col::bawl::v);
  const int isv = emc2col::select_index(split, emc2col::bawlsplit::sigma, emc2col::bawl::sv);
  const int iB = emc2col::select_index(split, emc2col::bawlsplit::B, emc2col::bawl::B);
  const int iA = emc2col::select_index(split, emc2col::bawlsplit::A, emc2col::bawl::A);
  const int it0 = emc2col::select_index(split, emc2col::bawlsplit::t0, emc2col::bawl::t0);
  const int ik = emc2col::select_index(split, emc2col::bawlsplit::k, emc2col::bawl::k);
  const int imG = emc2col::select_index(split, emc2col::bawlsplit::mG, emc2col::bawl::mG);
  const int imK = emc2col::select_index(split, emc2col::bawlsplit::mK, emc2col::bawl::mK);
  const double delta = split ? par[emc2col::bawlsplit::delta] : 0.0;
  // The ratio chart divides by the launch MEAN: column 0 for the normal (v)
  // and lognormal (mean) launches, column 1 for the Weibull pair (shape, mean).
  const int imean = (lau == BAWL_LAUNCH_WEIBULL) ? isv : iv;
  if (R_IsNA(par[iv])) return 0.0;
  const double t0_val = par[it0];
  const double tt = t - t0_val;
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    if (tt <= 0.0 || t <= 0.0) return 0.0;
    if (tt == R_PosInf) {
      if (lau == BAWL_LAUNCH_NORMAL)
        return (ctx && ctx->use_posdrift) ? 1.0 :
          pnorm_std(par[iv] / par[isv], true, false);
      return 1.0;
    }
    return bawl_cdf_scalar_natural(
      tt, par[iA], par[iB] + par[iA], par[iv], par[isv], 0.0,
      ctx->use_posdrift, LBA_DENOM_FLOOR, lau, delta);
  }
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const int ks = ctx ? ctx->kill_shape : 1;
  const double k_val = (ctx && ctx->bawl_k_fixed_zero) ? 0.0
                                                       : bawl_leak(ctx, par[ik], par[imean]);
  const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active) ? erlang_lambda_from_mean(par[imG], ks) : 0.0);
  const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
    ((ctx && ctx->kill_active && ctx->apply_lk_to_racers) ?
      erlang_lambda_from_mean(par[imK], ks) : 0.0);
  const double omega = erlang_omega_for_shape(ks, par, ctx ? ctx->erlang_omega_index : -1);
  const bool erl = (lg > 1e-12 || lk > 1e-12);
  if (tt <= 0.0 && !erl) return 0.0;
  if (t <= 0.0) return 0.0;
  return pkilledleakyba_norm(
    t, par[iv], par[iB] + par[iA], par[iA], par[isv], t0_val, k_val, lg, lk,
    ctx->use_posdrift, false, ks, local_guess, omega, lau, delta);
}

void dbawl_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
  // Lognormal launch strengths occupy the v/sv slots as (mu, sigma); posdrift
  // and the normalizer floor are then inert inside the kernels.
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const bool split = lau == BAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = emc2col::select_index(split, emc2col::bawlsplit::mu, emc2col::bawl::v);
  const int isv = emc2col::select_index(split, emc2col::bawlsplit::sigma, emc2col::bawl::sv);
  const int iB = emc2col::select_index(split, emc2col::bawlsplit::B, emc2col::bawl::B);
  const int iA = emc2col::select_index(split, emc2col::bawlsplit::A, emc2col::bawl::A);
  const int it0 = emc2col::select_index(split, emc2col::bawlsplit::t0, emc2col::bawl::t0);
  const int ik = emc2col::select_index(split, emc2col::bawlsplit::k, emc2col::bawl::k);
  const int imG = emc2col::select_index(split, emc2col::bawlsplit::mG, emc2col::bawl::mG);
  const int imK = emc2col::select_index(split, emc2col::bawlsplit::mK, emc2col::bawl::mK);
  const int iomega = emc2col::select_index(split, emc2col::bawlsplit::omega, emc2col::bawl::omega);
  const double* v_  = cols[iv];
  const double* sv_ = cols[isv];
  const double* B_  = cols[iB];
  const double* A_  = cols[iA];
  const double* t0_ = cols[it0];
  const double* delta_ = split ? cols[emc2col::bawlsplit::delta] : nullptr;
  // The ratio chart divides by the launch MEAN: v_ for the normal and
  // lognormal launches, the second column for the Weibull pair (shape, mean).
  const double* mean_ = (lau == BAWL_LAUNCH_WEIBULL) ? sv_ : v_;
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    // Exact LBA member: retain the shared numerical kernel, but skip the
    // killed-clock wrapper and all optional-column bookkeeping per row.
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) {
        out[i] = raw_log_zero(min_ll, floor_raw);
        continue;
      }
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0 || rt[i] <= 0.0) {
        out[i] = raw_log_zero(min_ll, floor_raw);
        continue;
      }
      const double b_i = B_[i] + A_[i];
      double pdf = 0.0;
      if (ba_natural_pdf_launch(tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd,
                                LBA_DENOM_FLOOR, BA_ACCEPT_RAW, pdf, lau,
                                delta_ ? delta_[i] : 0.0)) {
        out[i] = (pdf > 0.0)
          ? raw_log_value(std::log(pdf), min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      } else {
        const double log_pdf = log_ba_pdf_launch(
          tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd, LBA_DENOM_FLOOR, lau,
          delta_ ? delta_[i] : 0.0);
        out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
          ? raw_log_value(log_pdf, min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      }
    }
    return;
  }
  const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* k_  = (ctx && ctx->bawl_k_fixed_zero) ? nullptr : cols[ik];
  const double* lg_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[imG];
  const double* lk_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[imK];
  const double* omega_ = (ks == 3) ? cols[iomega] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double kval = (ctx && ctx->bawl_k_fixed_zero) ? 0.0
                                                        : bawl_leak(ctx, k_[i], mean_[i]);
    const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[i], ks));
    const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[i], ks));
    const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(ks);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (rt[i] <= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    if (!erl) {
      // The no-clock BAwL member is evaluated directly in natural space in
      // the scalar adapter.  Keep that fast path in the raw callback too;
      // routing every conditional GH node through dkilledleakyba_norm(...,
      // log_out=true) needlessly enters the strict log wrapper and is costly
      // in the correlated model.  The log evaluator remains the tail fallback.
      double pdf = 0.0;
      if (ba_natural_pdf_launch(tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval,
                                pd, BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, pdf, lau,
                                delta_ ? delta_[i] : 0.0)) {
        out[i] = (pdf > 0.0)
          ? raw_log_value(std::log(pdf), min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      } else {
        const double log_pdf = log_ba_pdf_launch(
          tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval, pd,
          BAWL_DENOM_FLOOR, lau, delta_ ? delta_[i] : 0.0);
        out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
          ? raw_log_value(log_pdf, min_ll, floor_raw)
          : raw_log_zero(min_ll, floor_raw);
      }
      continue;
    }
    const double log_pdf = dkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, kval, lg, lk,
      pd, true, ks, local_guess, omega, lau, delta_ ? delta_[i] : 0.0
    );
    out[i] = (log_pdf > R_NegInf && emc2_isfinite(log_pdf))
      ? raw_log_value(log_pdf, min_ll, floor_raw) : raw_log_zero(min_ll, floor_raw);
  }
}

void pbawl_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const bool pd          = ctx->use_posdrift;
  // Lognormal launch strengths occupy the v/sv slots as (mu, sigma); posdrift
  // and the normalizer floor are then inert inside the kernels.
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const bool split = lau == BAWL_LAUNCH_SPLITLOGNORMAL;
  const int iv = emc2col::select_index(split, emc2col::bawlsplit::mu, emc2col::bawl::v);
  const int isv = emc2col::select_index(split, emc2col::bawlsplit::sigma, emc2col::bawl::sv);
  const int iB = emc2col::select_index(split, emc2col::bawlsplit::B, emc2col::bawl::B);
  const int iA = emc2col::select_index(split, emc2col::bawlsplit::A, emc2col::bawl::A);
  const int it0 = emc2col::select_index(split, emc2col::bawlsplit::t0, emc2col::bawl::t0);
  const int ik = emc2col::select_index(split, emc2col::bawlsplit::k, emc2col::bawl::k);
  const int imG = emc2col::select_index(split, emc2col::bawlsplit::mG, emc2col::bawl::mG);
  const int imK = emc2col::select_index(split, emc2col::bawlsplit::mK, emc2col::bawl::mK);
  const int iomega = emc2col::select_index(split, emc2col::bawlsplit::omega, emc2col::bawl::omega);
  const double* v_  = cols[iv];
  const double* sv_ = cols[isv];
  const double* B_  = cols[iB];
  const double* A_  = cols[iA];
  const double* t0_ = cols[it0];
  const double* delta_ = split ? cols[emc2col::bawlsplit::delta] : nullptr;
  // The ratio chart divides by the launch MEAN: v_ for the normal and
  // lognormal launches, the second column for the Weibull pair (shape, mean).
  const double* mean_ = (lau == BAWL_LAUNCH_WEIBULL) ? sv_ : v_;
  if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
    // Exact LBA member; pfun writes the log-survivor expected by the raw
    // likelihood path, while the launch dispatch supplies stable tails.
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (R_IsNA(v_[i]) || !isok[i]) {
        out[i] = 0.0;
        continue;
      }
      const double tt = rt[i] - t0_[i];
      if (tt <= 0.0 || rt[i] <= 0.0) {
        out[i] = 0.0;
        continue;
      }
      if (tt == R_PosInf) {
        const double ls = log_ba_surv_launch(
          tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], 0.0, pd,
          LBA_DENOM_FLOOR, lau, delta_ ? delta_[i] : 0.0);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
        continue;
      }
      const double b_i = B_[i] + A_[i];
      double cdf = 0.0;
      if (ba_natural_cdf_launch(
            tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd,
            LBA_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, lau,
            delta_ ? delta_[i] : 0.0)) {
        out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
      } else {
        const double ls = log_ba_surv_launch(
          tt, A_[i], b_i, v_[i], sv_[i], 0.0, pd, LBA_DENOM_FLOOR, lau,
          delta_ ? delta_[i] : 0.0);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
      }
    }
    return;
  }
  const bool local_guess = ctx->is_local_guess || ctx->is_local_kill_guess;
  const int ks = ctx ? ctx->kill_shape : 1;
  const double* k_  = (ctx && ctx->bawl_k_fixed_zero) ? nullptr : cols[ik];
  const double* lg_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[imG];
  const double* lk_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[imK];
  const double* omega_ = (ks == 3) ? cols[iomega] : nullptr;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double t0_i = t0_[i];
    const double tt = rt[i] - t0_i;
    const double kval = (ctx && ctx->bawl_k_fixed_zero) ? 0.0
                                                        : bawl_leak(ctx, k_[i], mean_[i]);
    const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[i], ks));
    const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
      ((!ctx->kill_active || !ctx->apply_lk_to_racers) ? 0.0 : erlang_lambda_from_mean(lk_[i], ks));
    const double omega = (ks == 3 && omega_ != nullptr) ? std::fmax(0.0, std::fmin(1.0, omega_[i])) :
                         erlang_omega_for_shape(ks);
    const bool erl = (lg > 1e-12 || lk > 1e-12);
    if (tt <= 0.0 && !erl) { out[i] = 0.0; continue; }
    if (rt[i] <= 0.0) { out[i] = 0.0; continue; }
    if (!erl) {
      double cdf = 0.0;
      if (ba_natural_cdf_launch(
            tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval, pd,
            BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, lau,
            delta_ ? delta_[i] : 0.0)) {
        out[i] = (cdf > 0.0) ? std::log1p(-cdf) : 0.0;
      } else {
        const double ls = log_ba_surv_launch(
          tt, A_[i], B_[i] + A_[i], v_[i], sv_[i], kval, pd,
          BAWL_DENOM_FLOOR, lau, delta_ ? delta_[i] : 0.0);
        out[i] = (ls > R_NegInf && emc2_isfinite(ls))
          ? ls : raw_log_zero(min_ll, floor_raw);
      }
      continue;
    }
    const double log_cdf = pkilledleakyba_norm(
      rt[i], v_[i], B_[i] + A_[i], A_[i], sv_[i], t0_i, kval, lg, lk,
      pd, true, ks, local_guess, omega, lau, delta_ ? delta_[i] : 0.0
    );
    if (!R_FINITE(log_cdf)) { out[i] = 0.0; continue; }
    if (log_cdf >= 0.0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    out[i] = log1m_exp(log_cdf);
  }
}

// BAwL: column layout v=0, sv=1, B=2, A=3, t0=4, k=5, mG=6 (guess-clock mean), mK=7 (kill-clock mean)
void bawl_logS_at_t(double t, const double* const* cols,
                    int /*n_rows_total*/, int n_lR, int /*n_par*/,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  const bool pd = ctx->use_posdrift;
  const int lau = ctx ? ctx->bawl_launch : BAWL_LAUNCH_NORMAL;
  const bool split = lau == BAWL_LAUNCH_SPLITLOGNORMAL;
  const int ks = ctx ? ctx->kill_shape : 1;
  const int iv = emc2col::select_index(split, emc2col::bawlsplit::mu, emc2col::bawl::v);
  const int isv = emc2col::select_index(split, emc2col::bawlsplit::sigma, emc2col::bawl::sv);
  const int iB = emc2col::select_index(split, emc2col::bawlsplit::B, emc2col::bawl::B);
  const int iA = emc2col::select_index(split, emc2col::bawlsplit::A, emc2col::bawl::A);
  const int it0 = emc2col::select_index(split, emc2col::bawlsplit::t0, emc2col::bawl::t0);
  const int ik = emc2col::select_index(split, emc2col::bawlsplit::k, emc2col::bawl::k);
  const int imG = emc2col::select_index(split, emc2col::bawlsplit::mG, emc2col::bawl::mG);
  const int imK = emc2col::select_index(split, emc2col::bawlsplit::mK, emc2col::bawl::mK);
  const int iomega = emc2col::select_index(split, emc2col::bawlsplit::omega, emc2col::bawl::omega);
  const double* v_ = cols[iv];
  const double* sv_ = cols[isv];
  const double* B_ = cols[iB];
  const double* A_ = cols[iA];
  const double* t0_ = cols[it0];
  const double* delta_ = split ? cols[emc2col::bawlsplit::delta] : nullptr;
  // The ratio chart divides by the launch MEAN: v_ for the normal and
  // lognormal launches, the second column for the Weibull pair (shape, mean).
  const double* mean_ = (lau == BAWL_LAUNCH_WEIBULL) ? sv_ : v_;
  const double* k_ = (ctx && ctx->bawl_k_fixed_zero) ? nullptr : cols[ik];
  const double* lg_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[imG];
  const double* lk_ = (ctx && ctx->bawl_clocks_fixed_off) ? nullptr : cols[imK];
  const double* omega_ = (ks == 3) ? cols[iomega] : nullptr;

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int kk = 0; kk < n_lR && !bad; ++kk) {
      const int r = start + kk;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double t0_r = t0_[r];
      const double tt = t - t0_r;
      const double delta_r = delta_ ? delta_[r] : 0.0;
      if (ctx && ctx->bawl_k_fixed_zero && ctx->bawl_clocks_fixed_off) {
        if (tt <= 0.0) continue;
        double cdf = 0.0;
        if (ba_natural_cdf_launch(
              tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], 0.0, pd,
              LBA_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, lau, delta_r)) {
          if (cdf > 0.0) logS += std::log1p(-cdf);
        } else {
          const double ls = log_ba_surv_launch(
            tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], 0.0, pd,
            LBA_DENOM_FLOOR, lau, delta_r);
          if (!(ls > R_NegInf)) { bad = true; break; }
          logS += ls;
        }
        continue;
      }
      const double kval = (ctx && ctx->bawl_k_fixed_zero) ? 0.0
                                                          : bawl_leak(ctx, k_[r], mean_[r]);
      const double lg = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
        ((!ctx->kill_active) ? 0.0 : erlang_lambda_from_mean(lg_[r], ks));
      const double lk = (ctx && ctx->bawl_clocks_fixed_off) ? 0.0 :
        ((!ctx->kill_active || !ctx->apply_lk_to_racers) ?
          0.0 : erlang_lambda_from_mean(lk_[r], ks));
      const double omega = (ks == 3 && omega_ != nullptr) ?
        std::fmax(0.0, std::fmin(1.0, omega_[r])) : erlang_omega_for_shape(ks);
      const bool local_guess = ctx && (ctx->is_local_guess || ctx->is_local_kill_guess);
      const bool erl = (lg > 1e-12 || lk > 1e-12);
      if (tt <= 0.0) {
        if (!erl) continue;
        const double log_cdf = pkilledleakyba_norm(
          t, v_[r], B_[r] + A_[r], A_[r], sv_[r], t0_r, kval, lg, lk,
          pd, true, ks, local_guess, omega, lau, delta_r);
        if (!R_FINITE(log_cdf)) continue;
        if (log_cdf >= 0.0) { bad = true; break; }
        logS += log1m_exp(log_cdf);
        continue;
      }
      if (!erl) {
        double cdf = 0.0;
        if (ba_natural_cdf_launch(
              tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], kval, pd,
              BAWL_DENOM_FLOOR, BA_ACCEPT_RAW, cdf, lau, delta_r)) {
          if (cdf > 0.0) logS += std::log1p(-cdf);
        } else {
          const double ls = log_ba_surv_launch(
            tt, A_[r], B_[r] + A_[r], v_[r], sv_[r], kval, pd,
            BAWL_DENOM_FLOOR, lau, delta_r);
          if (!(ls > R_NegInf) || ISNAN(ls)) { bad = true; break; }
          logS += ls;
        }
        continue;
      }
      const double log_cdf = pkilledleakyba_norm(
        t, v_[r], B_[r] + A_[r], A_[r], sv_[r], t0_r, kval, lg, lk,
        pd, true, ks, local_guess, omega, lau, delta_r);
      if (log_cdf >= 0.0) { bad = true; break; }
      if (R_FINITE(log_cdf)) logS += log1m_exp(log_cdf);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

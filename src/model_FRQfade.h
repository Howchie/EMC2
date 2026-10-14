#ifndef EMC2_MODEL_FRQ_FADE_H
#define EMC2_MODEL_FRQ_FADE_H

#include "model_FRQ.h"

// ============================================================================
// FRQfade -- the FRQ counter with a fading registration opportunity.
//
// The proper FRQ accumulator (h = 1) runs on the warped time
//
//   A(u) = (1 - e^{-kappa u}) / kappa          (A(u) = u at kappa = 0),
//
// i.e. its input intensity is multiplied by e^{-kappa u}.  With F_1, f_1 and
// S_1 the FRQ functions at h = 1 for the same (K, v, sv, plo, clock),
//
//   F(u) = F_1(A(u)),   f(u) = f_1(A(u)) e^{-kappa u},   S(u) = S_1(A(u)).
//
// Operational time is capped at 1/kappa, so the accumulator never finishes
// with probability S_1(1/kappa): h = F_1(1/kappa) is derived, not fitted.
// The base sits at h = 1, so it needs no beta quantile at all.
// ============================================================================

struct FrqFadePars {
  const FrqPars* base = nullptr;
  double kappa = 0.0;
  double horizon = R_PosInf;   // 1/kappa, the operational-time cap
  bool ok = false;
};

inline double frqfade_warp(double u, const FrqFadePars& s) {
  if (s.kappa == 0.0) return u;
  if (!R_FINITE(u)) return s.horizon;
  return -std::expm1(-s.kappa * u) / s.kappa;
}

// Front-end memo: the base state goes through FrqMemo (and the per-particle
// FrqCache when the likelihood supplies one), kappa only sets the warp.
struct FrqFadeMemo {
  FrqMemo base;
  FrqFadePars value;

  explicit FrqFadeMemo(FrqCache* shared_cache = nullptr) : base(shared_cache) {}

  const FrqFadePars& get(double K, double v, double sv, double kappa,
                         double plo = 0.0, double clock = 1.0) {
    value.base = &base.get(K, v, sv, 1.0, plo, clock);
    value.ok = false;
    if (!value.base->ok) return value;
    if (ISNAN(kappa) || !(kappa >= 0.0) || !R_FINITE(kappa)) return value;
    value.kappa = kappa;
    value.horizon = (kappa > 0.0) ? 1.0 / kappa : R_PosInf;
    value.ok = true;
    return value;
  }
};

inline double frqfade_log_pdf_dt(double u, const FrqFadePars& s) {
  if (!s.ok) return R_NegInf;
  if (!(u > 0.0) || !R_FINITE(u)) return R_NegInf;
  const double lp = frq_log_pdf_dt(frqfade_warp(u, s), *s.base) - s.kappa * u;
  return ISNAN(lp) ? R_NegInf : lp;
}

inline double frqfade_log_cdf_dt(double u, const FrqFadePars& s) {
  if (!s.ok) return R_NegInf;
  if (!(u > 0.0)) return R_NegInf;
  return frq_log_cdf_dt(frqfade_warp(u, s), *s.base);
}

inline double frqfade_log_surv_dt(double u, const FrqFadePars& s) {
  if (!s.ok) return R_NaN;             // caller distinguishes bad from certain
  if (!(u > 0.0)) return 0.0;
  return frq_log_surv_dt(frqfade_warp(u, s), *s.base);
}

inline double frqfade_pdf_natural_dt(double u, const FrqFadePars& s) {
  const double lp = frqfade_log_pdf_dt(u, s);
  return (lp > R_NegInf) ? std::exp(lp) : 0.0;
}

inline double frqfade_cdf_natural_dt(double u, const FrqFadePars& s) {
  const double lp = frqfade_log_cdf_dt(u, s);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return (out > 1.0) ? 1.0 : out;
}

// Physical decision time of a base (unwarped) decision time x: +Inf, an
// omission, when x lies beyond the operational-time cap.
inline double frqfade_unwarp(double x, const FrqFadePars& s) {
  if (s.kappa == 0.0) return x;
  if (!(x < s.horizon)) return R_PosInf;
  return -std::log1p(-s.kappa * x) / s.kappa;
}

// Expected total number of events a trial can deliver, E[W] phi(1/kappa):
// v / kappa in the sv = 0 limit, beta G / (1 - G) at the cap otherwise, and
// Inf at kappa = 0.
inline double frqfade_total_evidence(const FrqFadePars& s) {
  if (!s.ok) return R_NaN;
  if (s.kappa == 0.0) return R_PosInf;
  const FrqPars& b = *s.base;
  if (b.poisson) return b.v * s.horizon;
  const FrqState z = frq_state(s.horizon, b);
  return b.beta * std::exp(z.log_q - z.log_omq);
}

double dfrqfade_scalar(double t, const double* par, void* ctx_);
double pfrqfade_scalar(double t, const double* par, void* ctx_);
void dfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_);
void pfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_);
void frqfade_logS_at_t(double t, const double* const* cols,
                       int n_rows_total, int n_lR, int n_par,
                       const int* trunc_mask, int n_unique_trials,
                       const int* isok_all, void* ctx_, double* logS_out);

#endif

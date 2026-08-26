#ifndef EMC2_TIME_WARP_H
#define EMC2_TIME_WARP_H

#include <Rcpp.h>
#include <cmath>
#include <string>

#include "race_contract.h"
#include "race_dispatch.h"

namespace emc2tw {

// s = c_eta(u) = ((1+u)^omega - 1)/omega, omega = exp(eta).
// The log1p/expm1 forms avoid direct exponentiation's overflow and cancellation regions.
inline double fwd(double u, double eta) {
  if (eta == 0.0) return u;                  // exact parent, bitwise
  if (!(u > 0.0)) return u;                  // u <= 0 and NaN pass through
  if (!R_FINITE(u)) return u;                // c(+Inf) = +Inf for every omega > 0
  if (std::isnan(eta)) return R_NaN;         // malformed eta propagates
  const double omega = std::exp(eta);
  if (!R_FINITE(omega)) return R_PosInf;     // omega -> Inf
  if (omega < 1e-12) return std::log1p(u);   // omega -> 0
  const double x = omega * std::log1p(u);
  if (x > 709.0) return R_PosInf;            // (1+u)^omega overflows
  return std::expm1(x) / omega;
}

// log c'_eta(u) = (omega - 1) log(1 + u).
inline double log_jac(double u, double eta) {
  if (eta == 0.0) return 0.0;
  if (!(u > 0.0) || !R_FINITE(u)) return 0.0;
  if (std::isnan(eta)) return R_NaN;
  const double omega = std::exp(eta);
  if (!R_FINITE(omega)) return R_PosInf;
  return (omega - 1.0) * std::log1p(u);
}

// u = c_eta^{-1}(s) = (1 + omega s)^(1/omega) - 1.
inline double inv(double s, double eta) {
  if (eta == 0.0) return s;
  if (!(s > 0.0)) return s;
  if (!R_FINITE(s)) return s;
  if (std::isnan(eta)) return R_NaN;
  const double omega = std::exp(eta);
  if (!R_FINITE(omega)) return 0.0;          // omega -> Inf
  if (omega < 1e-12) return (s > 709.0) ? R_PosInf : std::expm1(s);
  const double os = omega * s;
  // When omega*s overflows, log1p(omega*s) == log(omega) + log(s)
  // to well past double precision.
  const double y = (R_FINITE(os) ? std::log1p(os)
                                 : (std::log(omega) + std::log(s))) / omega;
  if (y > 709.0) return R_PosInf;
  return std::expm1(y);
}

}  // namespace emc2tw

// Generic warp wrappers.  Signatures match RacePdf1Fun / RaceCdf1Fun /
// RaceRawFun / RaceLogSAtTFun so they can be dropped into RaceModelAdapter.
double tw_pdf1(double t, const double* par, void* ctx_);
double tw_cdf1(double t, const double* par, void* ctx_);
void   tw_d_raw(const double* rt, const double* const* cols, int n_rows,
                const int* mask, const int* isok, double* out,
                double min_ll, void* ctx_);
void   tw_p_raw(const double* rt, const double* const* cols, int n_rows,
                const int* mask, const int* isok, double* out,
                double min_ll, void* ctx_);
void   tw_logS_at_t(double t, const double* const* cols, int n_rows_total,
                    int n_lR, int n_par, const int* trunc_mask,
                    int n_unique_trials, const int* isok_all, void* ctx_,
                    double* logS_out);

// Resolve the eta column by name and install the wrappers.  No-op when the
// design has no eta; hard error when it has one but the model does not support
// the warp.  Mirrors configure_corr_drift_context().
void configure_time_warp_context(RaceModelAdapter& adapter,
                                 const Rcpp::CharacterVector& keep_names,
                                 const std::string& caller);

#endif  // EMC2_TIME_WARP_H

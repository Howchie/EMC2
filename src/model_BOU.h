#ifndef model_BOU_h
#define model_BOU_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Likelihood primitives for the bounded OU -- the DDM with leak.
//
// These are the two operations c_log_likelihood_DDM_pt needs from a two-boundary
// model (see DDMAdapter in utils.h): the LOG defective density and the LOG
// defective cdf of a response, vectorised over trials.  Everything else the
// kernel does -- truncation, censoring, go/no-go, the all-finite fast path -- is
// model independent and is reused verbatim.
//
// Both are shaped exactly like d_DDM_Wien_raw / p_DDM_Wien_raw (model_DDM.h:53,
// 133), including the conventions that are easy to get wrong:
//   * Rs: 1 = lower, 2 = upper
//   * rows with mask == 0 are skipped and `out` is left untouched
//   * SZ arrives RAW and is widened here by 2*SZ*min(Z,1-Z), because Ttransform
//     is an R-side step that does not run on this path
//   * st0 is Uniform(t0, t0+st0), the lower-edge form
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <RcppArmadillo.h>
#include <cmath>
#include "col_registry.h"
#include "fpe_bou.h"

namespace bou {

// Pull one row's parameters out of the column base pointers.
struct Row {
  double v, a, sv, t0, st0, s, Z, SZ, beta, z, sz;
  bool ok;
};

inline Row bou_row(const double* const* cols, int i) {
  Row r;
  r.v    = cols[emc2col::bou::v][i];
  r.a    = cols[emc2col::bou::a][i];
  r.sv   = cols[emc2col::bou::sv][i];
  r.t0   = cols[emc2col::bou::t0][i];
  r.st0  = cols[emc2col::bou::st0][i];
  r.s    = cols[emc2col::bou::s][i];
  r.Z    = cols[emc2col::bou::Z][i];
  r.SZ   = cols[emc2col::bou::SZ][i];
  r.beta = cols[emc2col::bou::beta][i];
  r.z = r.Z * r.a;
  // Same widening d_DDM_Wien_raw does from its own raw SZ column.
  r.sz = 2.0 * r.SZ * std::min(r.z, r.a - r.z);
  r.ok = (r.a > 0.0) && (r.s > 0.0) && (r.z > 0.0) && (r.z < r.a) &&
         R_FINITE(r.v) && R_FINITE(r.beta) && (r.beta >= 0.0) &&
         R_FINITE(r.t0);
  return r;
}

// The horizon has to be known before the first solve, or every quadrature node
// is re-solved as later rows arrive (see fpebou::SolveCache::t_horizon).  Grown
// monotonically: a later call with shorter times reuses the existing solves
// rather than shrinking the horizon and invalidating them.
inline void bou_set_horizon(fpebou::SolveCache& C, const double* rts,
                            const double* const* cols, int n_rows,
                            const int* mask, const int* is_ok) {
  double h = C.t_horizon;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i] || !is_ok[i]) continue;
    const double t = rts[i];
    if (!R_FINITE(t)) continue;
    const double st0 = cols[emc2col::bou::st0][i];
    const double td = t - (cols[emc2col::bou::t0][i] -
                           0.5 * std::max(0.0, st0));
    if (td > h) h = td;
  }
  C.t_horizon = std::max(h * 1.02, 1e-3);
}

inline void bou_eval(const double* rts, const int* Rs,
                     const double* const* cols, int n_rows,
                     const int* mask, const int* is_ok,
                     double* out, double floor_,
                     ContextForDDMModels* ctx, bool want_cdf) {
  if (ctx == nullptr || !ctx->bou_cache) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = floor_;
    return;
  }
  fpebou::SolveCache& C = *ctx->bou_cache;
  bou_set_horizon(C, rts, cols, n_rows, mask, is_ok);

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!is_ok[i]) { out[i] = floor_; continue; }
    const double t = rts[i];
    if (!R_FINITE(t)) {
      // A cdf at t = +Inf is the response's total probability; the kernel avoids
      // asking, but answer sensibly rather than solving to an infinite horizon.
      out[i] = floor_;
      continue;
    }
    const Row r = bou_row(cols, i);
    if (!r.ok) { out[i] = floor_; continue; }

    const fpebou::BouMix m = fpebou::bou_mix_t0(
      C, t, r.v, r.sv, r.a, r.z, r.sz, r.s, r.beta,
      /*anchor_at_z=*/true, 0.5 * r.a, r.t0, r.st0,
      ctx->bnd_kind, 0.0, 0.0, 0.0, want_cdf);

    const bool up = (Rs[i] == 2);
    const double val = want_cdf ? (up ? m.F_up : m.F_lo)
                                : (up ? m.d_up : m.d_lo);
    out[i] = (val > 0.0) ? std::log(val) : R_NegInf;
    if (!R_FINITE(out[i]) && !want_cdf) out[i] = R_NegInf;
  }
}

inline void d_BOU_raw(const double* rts, const int* Rs,
                      const double* const* cols, int n_rows,
                      const int* mask, const int* is_ok,
                      double* out, double floor_, ContextForDDMModels* ctx) {
  bou_eval(rts, Rs, cols, n_rows, mask, is_ok, out, floor_, ctx, false);
}

inline void p_BOU_raw(const double* rts, const int* Rs,
                      const double* const* cols, int n_rows,
                      const int* mask, const int* is_ok,
                      double* out, double floor_, ContextForDDMModels* ctx) {
  bou_eval(rts, Rs, cols, n_rows, mask, is_ok, out, floor_, ctx, true);
}

// Read the solver's grid settings from R options once per adapter build, never
// per row -- the same arrangement rou_configure_grid uses (model_ROU.h:35).
inline void bou_configure(fpebou::SolveCache& C) {
  Rcpp::Environment base("package:base");
  Rcpp::Function getOption = base["getOption"];
  auto opt_i = [&](const char* nm, int def) {
    SEXP v = getOption(nm, Rcpp::wrap(def));
    return Rf_isNull(v) ? def : Rcpp::as<int>(v);
  };
  auto opt_d = [&](const char* nm, double def) {
    SEXP v = getOption(nm, Rcpp::wrap(def));
    return Rf_isNull(v) ? def : Rcpp::as<double>(v);
  };
  C.grid.nx        = opt_i("emc2.bou_nx", 384);
  C.grid.dt_target = opt_d("emc2.bou_dt", 5e-4);
  C.grid.grade     = opt_d("emc2.bou_grade", fpe::FPE_GRADE_BOUNDED);
  C.grid.tgrade    = opt_d("emc2.bou_tgrade", fpe::FPE_TGRADE);
  C.grid.n_sv      = opt_i("emc2.bou_n_sv", 7);
  C.grid.n_sz      = opt_i("emc2.bou_n_sz", 7);
  C.grid.n_st0     = opt_i("emc2.bou_n_st0", 7);
}

} // namespace bou

#endif // model_BOU_h

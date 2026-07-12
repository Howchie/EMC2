#ifndef ss_raw_h
#define ss_raw_h

// ---------------------------------------------------------------------------
// Raw-buffer stop-signal kernels shared by SSEXG and SSRDEX.
//
// The stop-signal hot path used to pass Rcpp NumericMatrix/LogicalVector
// objects into every per-trial-per-particle call, and the stop-success
// integrands re-read matrix elements (and re-derived Wald parameters) at every
// quadrature node. This header defines the allocation-free replacements:
//
//  * a canonical per-accumulator parameter block of 4 doubles
//      TEXG go/ST accumulator:  {mu, sigma, tau, lb}
//      RDEX go/ST accumulator:  {alpha, nu, gamma, t0}   (s already divided out)
//    plus a canonical stop block {muS, sigS, tauS, lbS} (both models race a
//    truncated ex-Gaussian stop process);
//  * SsRawModel: per-model function pointers on those blocks (log-pdf,
//    log-survivor, natural-scale survivor, block fill from ParamTable columns);
//  * SsStopCtx + the stop-success integral strategies (integrate / GL / auto /
//    live), numerically identical twins of the NumericMatrix versions in
//    model_SS_EXG.h / model_SS_RDEX.h — those are now thin wrappers over these
//    (so the exported ss_*_stop_success_value() test entry points exercise this
//    code).
//
// Model-specific edge behaviour is preserved exactly:
//   TEXG integrate/GL return min_ll on an empty window; RDEX clamps the upper
//   bound to lo + 1e-12 (clamp_empty_window). Only TEXG has the n_go == 1,
//   upper == Inf closed form (try_analytic1).
// ---------------------------------------------------------------------------

#include <cmath>
#include <vector>
#include "utility_functions.h"
#include "composite_functions.h"
#include "exgaussian_functions.h"
#include "wald_functions.h"
#include "gsl_utils.h"
#include "gl_quad.h"
#include "ss_exg_analytic.h"
#include "col_registry.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// number of doubles per canonical accumulator block
constexpr int SS_ACC_STRIDE = 4;

// --- per-model raw kernels on canonical blocks ------------------------------

using ss_acc_fill_fn  = void (*)(const double* const* cols, int row, double* acc4);
using ss_stop_fill_fn = void (*)(const double* const* cols, int row, double* stop4);
using ss_acc_lpdf_fn  = double (*)(double t, const double* acc4);   // log density
using ss_acc_lsurv_fn = double (*)(double t, const double* acc4);   // log survivor
using ss_acc_surv_fn  = double (*)(double t, const double* acc4);   // survivor

struct SsRawModel {
  ss_acc_fill_fn  fill_acc;
  ss_stop_fill_fn fill_stop;
  ss_acc_lpdf_fn  acc_lpdf;
  ss_acc_lsurv_fn acc_lsurv;
  ss_acc_surv_fn  acc_surv;
  bool clamp_empty_window;  // RDEX: ub -> lo + 1e-12; TEXG: return min_ll
  bool try_analytic1;       // TEXG only: n_go==1, upper==Inf closed form
  int idx_tf;
  int idx_gf;
};

// --- stop-success integral ---------------------------------------------------

// Everything the integrand needs, extracted once per call: quadrature nodes do
// pure arithmetic on pre-fetched doubles (no matrix reads, no re-derivation).
struct SsStopCtx {
  double SSD;
  const double* acc_go;   // canonical go blocks, SS_ACC_STRIDE per accumulator
  int n_go;
  double muS, sigS, tauS, lbS;   // stop block (first GO row, as before)
  ss_acc_surv_fn acc_surv;
};

// Integrand: f_stop(x) * prod_i S_go_i(x + SSD). Matches the lambda bodies of
// the NumericMatrix versions term for term (same early exits, same order).
inline double ss_stop_success_integrand_raw(double x, void* p) {
  const SsStopCtx* c = static_cast<const SsStopCtx*>(p);
  const double fS = dtexg(x, c->muS, c->sigS, c->tauS, c->lbS, R_PosInf, false);
  if (fS <= 0.0) return 0.0;
  double S_go_all = 1.0;
  for (int i = 0; i < c->n_go; ++i) {
    S_go_all *= c->acc_surv(x + c->SSD, c->acc_go + SS_ACC_STRIDE * i);
    if (S_go_all <= 0.0) return 0.0;
  }
  return fS * S_go_all;
}

// Adaptive GSL route (stop_method = "integrate"). Twin of
// ss_texg_stop_success_lpdf / ss_rdex_stop_success_lpdf.
inline double ss_stop_success_raw_integrate(
    const SsStopCtx& c, double min_ll, double upper,
    bool clamp_empty_window,
    int max_subdiv = 100, double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU
) {
  const double ub_heur = c.muS + k_sigma * c.sigS + k_tau * c.tauS;
  double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > c.lbS)) {
    if (!clamp_empty_window) return min_ll;
    ub = c.lbS + 1e-12;
  }

  gsl_function F;
  F.function = &ss_stop_success_integrand_raw;
  F.params = const_cast<SsStopCtx*>(&c);

  static thread_local GslWorkspacePtr ws_ptr(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws_ptr, max_subdiv);
  double res, err;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status;
  if (emc2_isinf(ub)) {
    status = gsl_integration_qagiu(&F, c.lbS, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  } else {
    status = gsl_integration_qags(&F, c.lbS, ub, abs_tol, rel_tol, max_subdiv, workspace, &res, &err);
  }
  gsl_set_error_handler(old_handler);

  if (status != GSL_SUCCESS || !emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// Fixed Gauss-Legendre route (stop_method = "gl"). Twin of
// ss_texg_stop_success_lpdf_gl / ss_rdex_stop_success_lpdf_gl.
inline double ss_stop_success_raw_gl(
    const SsStopCtx& c, double min_ll, double upper, int n_nodes,
    bool clamp_empty_window,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU
) {
  const double ub_heur = c.muS + k_sigma * c.sigS + k_tau * c.tauS;
  // lbS = -Inf (untruncated stop) needs a finite GL window: mirror the upper
  // heuristic on the lower side
  const double lo = emc2_isfinite(c.lbS) ? c.lbS
                                         : c.muS - k_sigma * c.sigS - k_tau * c.tauS;
  double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lo)) {
    if (!clamp_empty_window) return min_ll;
    ub = lo + 1e-12;
  }
  double res = gl_integrate(&ss_stop_success_integrand_raw,
                            const_cast<SsStopCtx*>(&c), lo, ub, n_nodes);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// "auto" dispatch. Twin of ss_texg_stop_success_lpdf_autodisp /
// ss_rdex_stop_success_lpdf_autodisp: TEXG n_go == 1 with infinite upper limit
// takes the closed form; otherwise GL with the tight-stop-density node bump.
inline double ss_stop_success_raw_auto(
    const SsStopCtx& c, double min_ll, double upper, int n_nodes,
    bool clamp_empty_window, bool try_analytic1,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU
) {
  if (try_analytic1 && c.n_go == 1 && !emc2_isfinite(upper)) {
    bool ok = false;
    const double* g = c.acc_go;
    double p = ss_texg_stop_success_analytic1(
        c.SSD, g[0], g[1], g[2], g[3], c.muS, c.sigS, c.tauS, c.lbS, ok);
    if (ok) return std::log(p);
  }
  const double lo = emc2_isfinite(c.lbS) ? c.lbS
                                         : c.muS - k_sigma * c.sigS - k_tau * c.tauS;
  const double ub = emc2_isfinite(upper) ? upper
                                         : c.muS + k_sigma * c.sigS + k_tau * c.tauS;
  const int n_eff = gl_auto_nodes(n_nodes, lo, ub, c.sigS);
  return ss_stop_success_raw_gl(c, min_ll, upper, n_eff, clamp_empty_window,
                                k_sigma, k_tau);
}

// LIVE entry point: honours the process-global stop_method_config() exactly
// like ss_texg_stop_success_lpdf_live / ss_rdex_stop_success_lpdf_live.
inline double ss_stop_success_raw_live(
    const SsStopCtx& c, double min_ll, double upper,
    bool clamp_empty_window, bool try_analytic1
) {
  const StopMethodConfig& cfg = stop_method_config();
  switch (cfg.method) {
  case STOP_METHOD_INTEGRATE:
    return ss_stop_success_raw_integrate(c, min_ll, upper, clamp_empty_window);
  case STOP_METHOD_GL:
    return ss_stop_success_raw_gl(c, min_ll, upper, cfg.n_nodes,
                                  clamp_empty_window);
  default:   // STOP_METHOD_AUTO / STOP_METHOD_ANALYTIC
    return ss_stop_success_raw_auto(c, min_ll, upper, cfg.n_nodes,
                                    clamp_empty_window, try_analytic1);
  }
}

// --- TEXG raw kernels (SSEXG live path) --------------------------------------
// Column order per emc2col::ss_texg (src/col_registry.h).

inline void ss_texg_fill_acc_raw(const double* const* cols, int row, double* a) {
  namespace tc = emc2col::ss_texg;
  a[0] = cols[tc::mu][row];
  a[1] = cols[tc::sigma][row];
  a[2] = cols[tc::tau][row];
  a[3] = cols[tc::exg_lb][row];
}
inline void ss_texg_fill_stop_raw(const double* const* cols, int row, double* s) {
  namespace tc = emc2col::ss_texg;
  s[0] = cols[tc::muS][row];
  s[1] = cols[tc::sigmaS][row];
  s[2] = cols[tc::tauS][row];
  s[3] = cols[tc::exgS_lb][row];
}
inline double ss_texg_acc_lpdf_raw(double t, const double* a) {
  return dtexg(t, a[0], a[1], a[2], a[3], R_PosInf, true);
}
inline double ss_texg_acc_lsurv_raw(double t, const double* a) {
  return ptexg(t, a[0], a[1], a[2], a[3], R_PosInf, false, true);
}
inline double ss_texg_acc_surv_raw(double t, const double* a) {
  return ptexg(t, a[0], a[1], a[2], a[3], R_PosInf, false, false);
}

inline const SsRawModel& ss_texg_raw_model() {
  static const SsRawModel m = {
    &ss_texg_fill_acc_raw, &ss_texg_fill_stop_raw,
    &ss_texg_acc_lpdf_raw, &ss_texg_acc_lsurv_raw, &ss_texg_acc_surv_raw,
    /*clamp_empty_window=*/false, /*try_analytic1=*/true,
    /*idx_tf=*/emc2col::ss_texg::tf, /*idx_gf=*/emc2col::ss_texg::gf
  };
  return m;
}

// --- RDEX raw kernels (SSRDEX live path) -------------------------------------
// Column order per emc2col::ss_rdex (src/col_registry.h). The canonical block
// pre-divides by s exactly as the NumericMatrix code did per evaluation:
//   alpha = B/s + A/(2s), nu = v/s, gamma = A/(2s), t0.

inline void ss_rdex_fill_acc_raw(const double* const* cols, int row, double* a) {
  namespace rc = emc2col::ss_rdex;
  const double s = cols[rc::s][row];
  a[0] = (cols[rc::B][row] / s) + .5 * (cols[rc::A][row] / s);  // alpha
  a[1] = cols[rc::v][row] / s;                                  // nu
  a[2] = .5 * (cols[rc::A][row] / s);                           // gamma
  a[3] = cols[rc::t0][row];
}
inline void ss_rdex_fill_stop_raw(const double* const* cols, int row, double* s) {
  namespace rc = emc2col::ss_rdex;
  s[0] = cols[rc::muS][row];
  s[1] = cols[rc::sigmaS][row];
  s[2] = cols[rc::tauS][row];
  s[3] = cols[rc::exgS_lb][row];
}
inline double ss_rdex_acc_lpdf_raw(double t, const double* a) {
  const double dt = t - a[3];
  if (!(dt > 0.)) return R_NegInf;
  return std::log(digt_impl(dt, a[0], a[1], a[2]));
}
inline double ss_rdex_acc_lsurv_raw(double t, const double* a) {
  const double dt = t - a[3];
  if (!(dt > 0.)) return 0.0;  // log(1)
  return log1m(pigt_impl(dt, a[0], a[1], a[2]));
}
inline double ss_rdex_acc_surv_raw(double t, const double* a) {
  const double dt = t - a[3];
  if (!(dt > 0.)) return 1.0;
  return 1.0 - pigt_impl(dt, a[0], a[1], a[2]);
}

inline const SsRawModel& ss_rdex_raw_model() {
  static const SsRawModel m = {
    &ss_rdex_fill_acc_raw, &ss_rdex_fill_stop_raw,
    &ss_rdex_acc_lpdf_raw, &ss_rdex_acc_lsurv_raw, &ss_rdex_acc_surv_raw,
    /*clamp_empty_window=*/true, /*try_analytic1=*/false,
    /*idx_tf=*/emc2col::ss_rdex::tf, /*idx_gf=*/emc2col::ss_rdex::gf
  };
  return m;
}

#endif

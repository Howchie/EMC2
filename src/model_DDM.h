#ifndef DDM_wien
#define DDM_wien

#include <RcppArmadillo.h>
using namespace Rcpp;
#include "utility_functions.h"
#include "pdf_fncs.h"
#include "cdf_fncs.h"
#include "fncs_seven.h"
#include "tools.h"
#include "ddm_functions_inline.h"
#include <cstdlib> // getenv

NumericVector d_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, LogicalVector is_ok){
  int Epsflag = 1;
  double eps = 5e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0; //the type of integration method to choose.
  //0 = "v", 1 = "a", 2= "sv", 3 = "t0", 4 = "st0", 5 = "s", 6 = "Z", 7 = "SZ",
  int N = rts.length();
  NumericVector out(N);
  for(int i = 0; i < N; i++){
    if(is_ok[i] == FALSE){
      out[i] = R_NegInf;
    } else{
      // we divide v, a and sv by s to introduce the scaling parameter s
      double pm = (Rs[i]==1) ? -1 : 1;
      // if sz and st0 are zero we can use simple and fast dwiener function
      if(pars(i,7) == 0 && pars(i, 4) == 0){
        double new_rt = rts[i] - pars(i,3);
        if(new_rt > 0){
          out[i] = dwiener_inline(new_rt*pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 6), pars(i, 2)/pars(i,5), eps, K, Epsflag);
        } else{
          out[i] = 	R_NegInf;
        }
      } else{ // otherwise use ddiff function with integration
        double Rval;
        double Rerr;
        double sz = (pars(i,6) < (1 - pars(i,6))) ? 2*pars(i,7)*pars(i,6) : 2*pars(i,7)*(1-pars(i,6));
        ddiff(choice, rts[i], pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 3), pars(i, 6), sz, pars(i, 2)/pars(i,5), pars(i,4), eps, K, Epsflag, Neval, &Rval, &Rerr);
        out[i] = (Rval > 0.0 && R_FINITE(Rval)) ? log(Rval) : R_NegInf;
      }
      // PDF log-probabilities should always be finite; clamp any numerical
      // pathologies to log(0) so caller-level min_ll flooring can handle them.
      if (!R_FINITE(out[i])) out[i] = R_NegInf;
    }
  }
  return(out);
}

inline void d_DDM_Wien_raw(const double* rts, const int* Rs, const double* pars_cm,
                           int n_rows, int n_pars, const int* mask,
                           const int* is_ok, double* out, double min_ll,
                           const std::vector<int>& p_idx) {
  (void)n_pars;
  int Epsflag = 1;
  double eps = 5e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0;

  const double* v_   = pars_cm + p_idx[0] * n_rows;
  const double* a_   = pars_cm + p_idx[1] * n_rows;
  const double* sv_  = pars_cm + p_idx[2] * n_rows;
  const double* t0_  = pars_cm + p_idx[3] * n_rows;
  const double* st0_ = pars_cm + p_idx[4] * n_rows;
  const double* s_   = pars_cm + p_idx[5] * n_rows;
  const double* Z_   = pars_cm + p_idx[6] * n_rows;
  const double* SZ_  = pars_cm + p_idx[7] * n_rows;

  // Check if SZ and st0 are zero for ALL trials in this batch.
  // If so, we can run a much faster vectorized loop.
  bool all_wdm = true;
  for (int i = 0; i < n_rows; ++i) {
    if (mask[i] && is_ok[i]) {
      if (SZ_[i] != 0.0 || st0_[i] != 0.0) {
        all_wdm = false;
        break;
      }
    }
  }

  if (all_wdm) {
    #pragma omp simd
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (!is_ok[i]) {
        out[i] = min_ll;
        continue;
      }
      const double pm = (Rs[i] == 1) ? -1.0 : 1.0;
      const double new_rt = rts[i] - t0_[i];
      if (new_rt > 0.0) {
        out[i] = dwiener_inline(new_rt * pm, a_[i] / s_[i], v_[i] / s_[i],
                                Z_[i], sv_[i] / s_[i], eps, K, Epsflag);
      } else {
        out[i] = min_ll;
      }
    }
  } else {
    // Mixed or full-DDM path: logic remains same but we still use inlined funcs
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      if (!is_ok[i]) {
        out[i] = min_ll;
        continue;
      }

      const double pm = (Rs[i] == 1) ? -1.0 : 1.0;
      if (SZ_[i] == 0.0 && st0_[i] == 0.0) {
        const double new_rt = rts[i] - t0_[i];
        if (new_rt > 0.0) {
          out[i] = dwiener_inline(new_rt * pm, a_[i] / s_[i], v_[i] / s_[i],
                                  Z_[i], sv_[i] / s_[i], eps, K, Epsflag);
        } else {
          out[i] = min_ll;
        }
      } else {
        double Rval;
        double Rerr;
        const double sz = (Z_[i] < (1.0 - Z_[i])) ? 2.0 * SZ_[i] * Z_[i]
                                                  : 2.0 * SZ_[i] * (1.0 - Z_[i]);
        ddiff(choice, rts[i], pm, a_[i] / s_[i], v_[i] / s_[i], t0_[i],
              Z_[i], sz, sv_[i] / s_[i], st0_[i], eps, K, Epsflag, Neval,
              &Rval, &Rerr);
        out[i] = (Rval > 0.0 && R_FINITE(Rval)) ? std::log(Rval) : min_ll;
      }
    }
  }
}

inline void p_DDM_Wien_raw(const double* rts, const int* Rs, const double* pars_cm,
                           int n_rows, int n_pars, const int* mask,
                           const int* is_ok, double* out, double min_ll,
                           const std::vector<int>& p_idx) {
  (void)n_pars;
  int Epsflag = 1;
  double eps = 5e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0;

  const double* v_   = pars_cm + p_idx[0] * n_rows;
  const double* a_   = pars_cm + p_idx[1] * n_rows;
  const double* sv_  = pars_cm + p_idx[2] * n_rows;
  const double* t0_  = pars_cm + p_idx[3] * n_rows;
  const double* st0_ = pars_cm + p_idx[4] * n_rows;
  const double* s_   = pars_cm + p_idx[5] * n_rows;
  const double* Z_   = pars_cm + p_idx[6] * n_rows;
  const double* SZ_  = pars_cm + p_idx[7] * n_rows;

  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!is_ok[i]) {
      out[i] = min_ll;
      continue;
    }

    const double new_rt = rts[i] - t0_[i];
    if (new_rt <= 0.0) {
      out[i] = min_ll;
    } else if (sv_[i] == 0.0 && SZ_[i] == 0.0 && st0_[i] == 0.0) {
      double v = v_[i] / s_[i];
      double w = Z_[i];
      if (Rs[i] != 1) { // upper boundary
        v = -v;
        w = 1.0 - w;
      }
      out[i] = pwiener_inline(new_rt, a_[i] / s_[i], v, w, eps, K, Epsflag);
    } else {
      double Rval;
      double Rerr;
      const double pm = (Rs[i] == 1) ? -1.0 : 1.0;
      const double sz = (Z_[i] < (1.0 - Z_[i])) ? 2.0 * SZ_[i] * Z_[i]
                                                : 2.0 * SZ_[i] * (1.0 - Z_[i]);
      pdiff(choice, rts[i], pm, a_[i] / s_[i], v_[i] / s_[i], t0_[i],
            Z_[i], sz, sv_[i] / s_[i], st0_[i], eps, K, Epsflag, Neval,
            &Rval, &Rerr);
      out[i] = (Rval > 0.0 && R_FINITE(Rval)) ? std::log(Rval) : min_ll;
    }
    
    if (!R_FINITE(out[i])) {
      out[i] = min_ll;
    } else if (out[i] > 0.0) {
      out[i] = 0.0;
    }
  }
}

NumericVector p_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, LogicalVector is_ok){
  static int ddm_debug_prints_left = 20; // shared across calls; keep noise bounded
  const bool ddm_debug = (std::getenv("EMC2_DEBUG_DDM") != nullptr);
  int Epsflag = 1;
  double eps = 5e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0; //the type of integration method to choose.
  //0 = "v", 1 = "a", 2= "sv", 3 = "t0", 4 = "st0", 5 = "s", 6 = "Z", 7 = "SZ",
  int N = rts.length();
  NumericVector out(N);
  for(int i = 0; i < N; i++){
    if(is_ok[i] == FALSE){
      out[i] = R_NegInf;
    } else{
      // we divide v, a and sv by s to introduce the scaling parameter s
      double new_rt = rts[i] - pars(i,3);
      if(new_rt <= 0){
        out[i] = R_NegInf;
      } else if(pars(i,2) == 0 && pars(i,7) == 0 && pars(i, 4) == 0){
        // if sv, sz and st0 are zero we can use simple and fast pwiener function
        // NOTE: Unlike dwiener(), pwiener() does not accept a signed time argument
        // to indicate the boundary. To get the upper-bound CDF, we reflect the
        // diffusion (v -> -v, w -> 1-w) and evaluate the lower-bound CDF.
        double v = pars(i, 0)/pars(i,5);
        double w = pars(i, 6);
        if (Rs[i] != 1) { // upper boundary
          v = -v;
          w = 1.0 - w;
        }
        out[i] = pwiener_inline(new_rt, pars(i, 1)/pars(i,5), v, w, eps, K, Epsflag);
        if (ddm_debug && ddm_debug_prints_left-- > 0 && !R_FINITE(out[i])) {
          Rcpp::Rcout << "[EMC2_DEBUG_DDM] pwiener returned non-finite logcdf: "
                      << out[i]
                      << " new_rt=" << new_rt
                      << " a=" << (pars(i, 1)/pars(i,5))
                      << " v=" << v
                      << " w=" << w
                      << " Rs=" << Rs[i]
                      << "\n";
        }
      } else{ // otherwise use pdiff function with integration
        double Rval;
        double Rerr;
        double pm = (Rs[i]==1) ? -1 : 1;
        double sz = (pars(i,6) < (1 - pars(i,6))) ? 2*pars(i,7)*pars(i,6) : 2*pars(i,7)*(1-pars(i,6));
        pdiff(choice, rts[i], pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 3), pars(i, 6), sz, pars(i, 2)/pars(i,5), pars(i,4), eps, K, Epsflag, Neval, &Rval, &Rerr);
        out[i] = (Rval > 0.0 && R_FINITE(Rval)) ? log(Rval) : R_NegInf;
      }
      // CDF log-probabilities should always be finite and <= 0; guard against
      // occasional numerical spillover in backend routines.
      if (!R_FINITE(out[i])) {
        out[i] = R_NegInf;
      } else if (out[i] > 0.0) {
        out[i] = 0.0;
      }
    }
  }
  return(out);
}


#endif

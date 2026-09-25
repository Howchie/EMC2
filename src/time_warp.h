#ifndef EMC2_TIME_WARP_H
#define EMC2_TIME_WARP_H

#include <Rcpp.h>
#include <cmath>
#include <string>

#include "race_contract.h"
#include "race_dispatch.h"
#include "time_clock.h"

// The clock family itself (emc2tw::) lives in time_clock.h.

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

// Resolve the clock column (TimeWarpPlan::par_name) by name and install the
// wrappers.  An optional clock (the ballistic eta) is a no-op when the design
// lacks it and a hard error when the model does not support it; a required
// clock (RDMSWTN_TT tau, RDMSWTN_UT u) must be present.  Mirrors
// configure_corr_drift_context().
void configure_time_warp_context(RaceModelAdapter& adapter,
                                 const Rcpp::CharacterVector& keep_names,
                                 const std::string& caller);

#endif  // EMC2_TIME_WARP_H

#ifndef EMC2_MODEL_RDMSWTN_GAMMA_H
#define EMC2_MODEL_RDMSWTN_GAMMA_H

#include "race_contract.h"
// RDMSWTN with Gamma-distributed drift (RDMSWTN(drift_dist = "gamma")):
// race-model adapter entry points.  Kernel in rdmswtn_gamma.h.
double drdmswtn_gdrift_scalar(double t, const double* par, void* ctx_);
double prdmswtn_gdrift_scalar(double t, const double* par, void* ctx_);
void drdmswtn_gdrift_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_);
void prdmswtn_gdrift_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_);
void rdmswtn_gdrift_logS_at_t(double t, const double* const* cols,
                              int n_rows_total, int n_lR, int n_par,
                              const int* trunc_mask, int n_unique_trials,
                              const int* isok_all, void* ctx_, double* logS_out);

// Decision-time log density / CDF / survivor (t already net of t0).
double gdrift_log_density(double t, double v, double b, double s, double sv,
                          double lf_floor = R_NegInf);
double gdrift_log_cdf(double t, double v, double b, double s, double sv);
double gdrift_log_surv(double t, double v, double b, double s, double sv);

#endif

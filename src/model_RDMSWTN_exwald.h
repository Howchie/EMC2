#ifndef EMC2_MODEL_RDMSWTN_EXWALD_H
#define EMC2_MODEL_RDMSWTN_EXWALD_H

#include "race_contract.h"
// RDMSWTN with exponential hesitation (RDMSWTN(exwald = TRUE)): the Wald
// first passage plus an independent Exp(mean tau) delay, sv = 0.
double drdmswtn_exwald_scalar(double t, const double* par, void* ctx_);
double prdmswtn_exwald_scalar(double t, const double* par, void* ctx_);
void drdmswtn_exwald_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_);
void prdmswtn_exwald_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_);
void rdmswtn_exwald_logS_at_t(double t, const double* const* cols,
                              int n_rows_total, int n_lR, int n_par,
                              const int* trunc_mask, int n_unique_trials,
                              const int* isok_all, void* ctx_, double* logS_out);

// Decision-time log density / CDF / survivor (t already net of t0).  The
// distance to threshold is uniform on [B, B + A]; tau = 0 is the Wald.
double exwald_log_density(double t, double v, double B, double A, double s, double tau);
double exwald_log_cdf(double t, double v, double B, double A, double s, double tau);
double exwald_log_surv(double t, double v, double B, double A, double s, double tau);

#endif

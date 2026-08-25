#ifndef EMC2_MODEL_GBM_H
#define EMC2_MODEL_GBM_H

#include "race_contract.h"

// GBM race-model adapter entry points.  Definitions live in model_GBM.cpp so
// utils.h remains a declaration-only integration point for these adapters.
double drdmgbm_scalar(double t, const double* par, void* ctx_);
double prdmgbm_scalar(double t, const double* par, void* ctx_);
void drdmgbm_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_);
void prdmgbm_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok,
                 double* out, double min_ll, void* ctx_);
void rdmgbm_logS_at_t(double t, const double* const* cols,
                      int n_rows_total, int n_lR, int n_par,
                      const int* trunc_mask, int n_unique_trials,
                      const int* isok_all, void* ctx_, double* logS_out);

#endif

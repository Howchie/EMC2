#ifndef EMC2_MODEL_BAWL_H
#define EMC2_MODEL_BAWL_H

#include "race_contract.h"

// BAwL/LBA race-model adapter entry points.
// Definitions live in model_BAwL.cpp so utils.h remains a
// declaration-only integration point for these adapters.
double dbawl_scalar(double t, const double* par, void* ctx_);
double pbawl_scalar(double t, const double* par, void* ctx_);
void dbawl_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void pbawl_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void bawl_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);

#endif

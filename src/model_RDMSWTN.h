#ifndef EMC2_MODEL_RDMSWTN_H
#define EMC2_MODEL_RDMSWTN_H

#include "race_contract.h"
// RDMSWTN race-model adapter entry points.  Definitions live in
// model_RDMSWTN.cpp so utils.h remains a declaration-only integration point
// for these adapters.
double drdmswtn_scalar(double t, const double* par, void* ctx_);
double prdmswtn_scalar(double t, const double* par, void* ctx_);
void drdmswtn_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_);
void prdmswtn_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_);
void rdmswtn_logS_at_t(double t, const double* const* cols,
                       int n_rows_total, int n_lR, int n_par,
                       const int* trunc_mask, int n_unique_trials,
                       const int* isok_all, void* ctx_, double* logS_out);

// RDMSWTN_TT race-model adapter entry points.  Definitions live in
// model_RDMSWTN.cpp so utils.h remains a declaration-only integration point
// for these adapters.
double drdmswtn_tt_scalar(double t, const double* par, void* ctx_);
double prdmswtn_tt_scalar(double t, const double* par, void* ctx_);
void drdmswtn_tt_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_);
void prdmswtn_tt_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_);
void rdmswtn_tt_logS_at_t(double t, const double* const* cols,
                          int n_rows_total, int n_lR, int n_par,
                          const int* trunc_mask, int n_unique_trials,
                          const int* isok_all, void* ctx_, double* logS_out);

#endif

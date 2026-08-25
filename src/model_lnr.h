#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "utility_functions.h"
#include "lnr_functions.h"

using namespace Rcpp;

// LNR race-model kernel entry points.  Definitions live in model_lnr.cpp so
// that utils.h (the shared analytic-kernel header) and this header can be
// compiled independently; the function-pointer types come from race_contract.h.
double dlnr_scalar(double t, const double* par, void* ctx_);
double plnr_scalar(double t, const double* par, void* ctx_);
void dlnr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_);
void plnr_raw(const double* rt, const double* const* cols, int n_rows,
              const int* mask, const int* isok,
              double* out, double min_ll, void* ctx_);
void lnr_logS_at_t(double t, const double* const* cols,
                   int n_rows_total, int n_lR, int n_par,
                   const int* trunc_mask, int n_unique_trials,
                   const int* isok_all, void* ctx_, double* logS_out);
#endif

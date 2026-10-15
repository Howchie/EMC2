#ifndef EMC2_STATE_MIXTURE_H
#define EMC2_STATE_MIXTURE_H

#include <Rcpp.h>
#include <string>

#include "race_contract.h"
#include "race_dispatch.h"

// Generic two-state accumulator mixture wrappers (see MixturePlan in
// race_contract.h).  Signatures match the adapter entry points they replace.
double mix_pdf1(double t, const double* par, void* ctx_);
double mix_cdf1(double t, const double* par, void* ctx_);
void   mix_d_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok, double* out,
                 double min_ll, void* ctx_);
void   mix_p_raw(const double* rt, const double* const* cols, int n_rows,
                 const int* mask, const int* isok, double* out,
                 double min_ll, void* ctx_);
void   mix_logS_at_t(double t, const double* const* cols, int n_rows_total,
                     int n_lR, int n_par, const int* trunc_mask,
                     int n_unique_trials, const int* isok_all, void* ctx_,
                     double* logS_out);

// Resolve `pmix` and the `<par>_mix` columns by name and install the wrappers.
// A no-op when the design has neither; an error when it has only one of them
// or when the race path does not evaluate accumulators independently.  Must
// run after configure_time_warp_context() so the mixture is the outer layer.
void configure_state_mixture_context(RaceModelAdapter& adapter,
                                     const Rcpp::CharacterVector& keep_names,
                                     bool is_logicalrules,
                                     const std::string& caller);

#endif  // EMC2_STATE_MIXTURE_H

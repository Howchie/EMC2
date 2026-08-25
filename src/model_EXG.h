#ifndef EMC2_MODEL_EXG_H
#define EMC2_MODEL_EXG_H

void drexg_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void prexg_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void rexg_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);

#endif

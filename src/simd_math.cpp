// Built with -ffinite-math-only (see Makevars.in), which is what lets the
// compiler use the C library's vector math variants.  Nothing here tests for
// Inf or NaN, and callers only pass finite arguments in each function's domain.
#include <math.h>
#include "simd_math.h"

void emc2_simd_exp(double* x, int n) {
#pragma omp simd
  for (int i = 0; i < n; ++i) x[i] = exp(x[i]);
}

void emc2_simd_log(double* x, int n) {
#pragma omp simd
  for (int i = 0; i < n; ++i) x[i] = log(x[i]);
}

void emc2_simd_log1p(double* x, int n) {
#pragma omp simd
  for (int i = 0; i < n; ++i) x[i] = log1p(x[i]);
}

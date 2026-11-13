#ifndef EMC2_SIMD_MATH_H
#define EMC2_SIMD_MATH_H

// In-place elementwise exp / log / log1p over dense arrays.  The definitions
// live in simd_math.cpp, which is built with finite-only math so the compiler
// can call the vector variants of these functions where the C library has
// them (glibc's libmvec); elsewhere they are plain loops.  Callers pass only
// finite arguments inside each function's domain and keep every other case on
// their scalar route.
void emc2_simd_exp(double* x, int n);
void emc2_simd_log(double* x, int n);
void emc2_simd_log1p(double* x, int n);

#endif  // EMC2_SIMD_MATH_H

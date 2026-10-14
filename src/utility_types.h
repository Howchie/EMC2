#ifndef EMC2_UTILITY_TYPES_H
#define EMC2_UTILITY_TYPES_H

#include <cmath>
#include <limits>

// Shared struct/enum type definitions used by utility_functions.h and transform_utils.h.
// No function bodies here — only types.

// For do_bounds
// A parameter may declare several permitted values outside its [min,max]
// interval (e.g. BTAwL's `pi` is allowed to sit on either endpoint, 0 for a
// pure transient race and 1 for a pure sustained one), so the exception is a
// small set rather than a single value.
#define EMC2_BOUND_MAX_EXCEPTIONS 4
struct BoundSpec {
  int col_idx;
  double min_val;
  double max_val;
  int n_exception;
  double exception_val[EMC2_BOUND_MAX_EXCEPTIONS];
};

// For transforms
enum TransformCode {
  IDENTITY = 0,
  EXP      = 1,
  PNORM    = 2
};

struct TransformSpec {
  int col_idx;
  TransformCode code;
  double lower;
  double upper;
  // Natural values at a sampled -Inf / +Inf; NaN means lower / upper. They
  // differ for a side linked to a declared bound, whose plain endpoint (0, or
  // 1 for a pnorm upper side) stays reachable so that log(0) and qnorm(0)
  // still give exception values such as t0 = 0 or pGuess = 0.
  double lower_inf = std::numeric_limits<double>::quiet_NaN();
  double upper_inf = std::numeric_limits<double>::quiet_NaN();
};

// For pretransform
enum PreTFCode { PTF_EXP = 1, PTF_PNORM = 2, PTF_NONE = 0 };

struct PreTransformSpec {
  int index;
  PreTFCode code;
  double lower;
  double upper;
};

#endif

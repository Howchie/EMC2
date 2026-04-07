#ifndef EMC2_UTILITY_TYPES_H
#define EMC2_UTILITY_TYPES_H

// Shared struct/enum type definitions used by utility_functions.h and transform_utils.h.
// No function bodies here — only types.

// For do_bounds
struct BoundSpec {
  int col_idx;
  double min_val;
  double max_val;
  bool has_exception;
  double exception_val;
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

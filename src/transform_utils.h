#ifndef EMC2_TRANSFORM_UTILS_H
#define EMC2_TRANSFORM_UTILS_H

#include <RcppArmadillo.h>
#include <unordered_set>
#include "ParamTable.h"

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

std::vector<TransformSpec>
make_transform_specs_for_paramtable(const ParamTable& pt,
                                    const Rcpp::List& transform);

std::vector<TransformSpec>
make_transform_specs_for_paramtable_from_full(
  const ParamTable& pt,
  const Rcpp::CharacterVector& full_names,
  const std::vector<TransformSpec>& full_specs);

void c_do_transform_pt(ParamTable& pt,
                       const std::vector<TransformSpec>& specs);

std::vector<TransformSpec> filter_specs_by_param_set(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  const std::unordered_set<std::string>& allowed);

std::vector<TransformSpec> complement_specs_for_premap(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  const std::unordered_set<std::string>& premap_set);

#endif

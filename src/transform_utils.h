#ifndef EMC2_TRANSFORM_UTILS_H
#define EMC2_TRANSFORM_UTILS_H

#include <RcppArmadillo.h>
#include <unordered_set>
#include "ParamTable.h"
#include "utility_types.h"

// TransformCode, TransformSpec, PreTransformSpec, BoundSpec are defined in utility_types.h

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

Rcpp::LogicalVector c_do_bound_pt(const ParamTable& pt, const std::vector<BoundSpec>& specs);

std::vector<BoundSpec> make_bound_specs_pt(Rcpp::NumericMatrix minmax,
                                           Rcpp::CharacterVector minmax_colnames,
                                           const ParamTable& pt,
                                           Rcpp::List bound);

std::vector<TransformSpec> filter_specs_by_param_set(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  const std::unordered_set<std::string>& allowed);

std::vector<TransformSpec> complement_specs_for_premap(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  const std::unordered_set<std::string>& premap_set);

std::vector<TransformSpec> complement_specs_for_phases(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  std::initializer_list<const std::unordered_set<std::string>*> sets);

#endif

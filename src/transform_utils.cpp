#include "transform_utils.h"
#include "wald_functions.h"
#include <unordered_map>

using namespace Rcpp;

std::vector<TransformSpec>
make_transform_specs_for_paramtable(const ParamTable& pt,
                                    const List& transform)
{
  CharacterVector func = transform["func"];
  NumericVector lowervec = transform["lower"];
  NumericVector uppervec = transform["upper"];

  std::unordered_map<std::string, TransformCode> codeMap;
  CharacterVector fnames = func.names();
  for (int i = 0; i < func.size(); i++) {
    std::string name = as<std::string>(fnames[i]);
    std::string f = as<std::string>(func[i]);
    if (f == "exp") codeMap[name] = EXP;
    else if (f == "pnorm") codeMap[name] = PNORM;
    else codeMap[name] = IDENTITY;
  }

  std::unordered_map<std::string, std::pair<double,double>> boundMap;
  {
    CharacterVector ln = lowervec.names();
    for (int i = 0; i < lowervec.size(); i++) {
      boundMap[as<std::string>(ln[i])].first = lowervec[i];
    }
    CharacterVector un = uppervec.names();
    for (int i = 0; i < uppervec.size(); i++) {
      boundMap[as<std::string>(un[i])].second = uppervec[i];
    }
  }

  std::vector<TransformSpec> specs;
  specs.reserve(pt.base_names.size());

  for (int col_idx = 0; col_idx < pt.base_names.size(); ++col_idx) {
    std::string name = as<std::string>(pt.base_names[col_idx]);
    TransformSpec s;
    s.col_idx = col_idx;

    auto itc = codeMap.find(name);
    s.code = (itc != codeMap.end()) ? itc->second : IDENTITY;

    auto itb = boundMap.find(name);
    if (itb != boundMap.end()) {
      s.lower = itb->second.first;
      s.upper = itb->second.second;
    } else {
      s.lower = 0.0;
      s.upper = 1.0;
    }

    specs.push_back(s);
  }

  return specs;
}

std::vector<TransformSpec>
make_transform_specs_for_paramtable_from_full(
  const ParamTable& pt,
  const CharacterVector& full_names,
  const std::vector<TransformSpec>& full_specs)
{
  std::unordered_map<std::string, int> name_to_full_idx;
  for (int k = 0; k < full_names.size(); ++k) {
    name_to_full_idx[as<std::string>(full_names[k])] = k;
  }

  std::vector<TransformSpec> out;
  out.reserve(pt.base_names.size());

  for (int j = 0; j < pt.base_names.size(); ++j) {
    std::string nm = as<std::string>(pt.base_names[j]);
    auto it = name_to_full_idx.find(nm);
    if (it != name_to_full_idx.end()) {
      TransformSpec s = full_specs[it->second];
      s.col_idx = j;
      out.push_back(s);
    }
  }
  return out;
}

Rcpp::LogicalVector c_do_bound_pt(const ParamTable& pt,
                                  const std::vector<BoundSpec>& specs)
{
  using Rcpp::LogicalVector;
  using Rcpp::NumericMatrix;

  const NumericMatrix& base = pt.base;
  const int nrows = base.nrow();

  LogicalVector result(nrows, true);

  for (std::size_t j = 0; j < specs.size(); ++j) {
    const BoundSpec& bs = specs[j];
    const int col_idx   = bs.col_idx;
    const double min_v  = bs.min_val;
    const double max_v  = bs.max_val;
    const bool has_exc  = bs.has_exception;
    const double exc_val= bs.exception_val;

    for (int i = 0; i < nrows; ++i) {
      const double val = base(i, col_idx);
      bool ok = (val > min_v && val < max_v);
      if (!ok && has_exc) {
        ok = (val == exc_val);
      }
      if (result[i] && !ok) {
        result[i] = false;
      }
    }
  }

  return result;
}

void c_do_transform_pt(ParamTable& pt,
                       const std::vector<TransformSpec>& specs)
{
  const int nrow = pt.n_trials;

  for (size_t j = 0; j < specs.size(); ++j) {
    const TransformSpec& sp = specs[j];
    const int col_idx = sp.col_idx;
    const TransformCode c = sp.code;
    const double lw = sp.lower;
    const double up = sp.upper;

    double* col = &pt.base(0, col_idx);

    switch (c) {
    case EXP: {
      for (int i = 0; i < nrow; ++i) col[i] = lw + std::exp(col[i]);
      break;
    }
    case PNORM: {
      const double range = up - lw;
      for (int i = 0; i < nrow; ++i) {
        col[i] = lw + range * pnorm_std(col[i]);
      }
      break;
    }
    case IDENTITY:
    default:
      break;
    }
  }
}

std::vector<TransformSpec> filter_specs_by_param_set(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  const std::unordered_set<std::string>& allowed)
{
  std::vector<TransformSpec> out;
  out.reserve(full_specs.size());

  for (const auto& sp : full_specs) {
    int base_idx = sp.col_idx;
    std::string nm = Rcpp::as<std::string>(pt.base_names[base_idx]);
    if (allowed.find(nm) != allowed.end()) out.push_back(sp);
  }
  return out;
}

std::vector<TransformSpec> complement_specs_for_premap(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  const std::unordered_set<std::string>& premap_set)
{
  return complement_specs_for_phases(pt, full_specs, { &premap_set });
}

std::vector<TransformSpec> complement_specs_for_phases(
  const ParamTable& pt,
  const std::vector<TransformSpec>& full_specs,
  std::initializer_list<const std::unordered_set<std::string>*> sets)
{
  std::vector<TransformSpec> out;
  out.reserve(full_specs.size());

  for (const auto& sp : full_specs) {
    int base_idx = sp.col_idx;
    std::string nm = Rcpp::as<std::string>(pt.base_names[base_idx]);
    bool excluded = false;
    for (const auto* s_ptr : sets) {
      if (s_ptr->find(nm) != s_ptr->end()) {
        excluded = true;
        break;
      }
    }
    if (!excluded) out.push_back(sp);
  }
  return out;
}

// Same logic as make_bound_specs but indexed into ParamTable base columns
std::vector<BoundSpec> make_bound_specs_pt(Rcpp::NumericMatrix minmax,
                                           Rcpp::CharacterVector minmax_colnames,
                                           const ParamTable& pt,
                                           Rcpp::List bound)
{
  using namespace Rcpp;
  using std::string;

  bool has_exception =
    bound.containsElementNamed("exception") &&
    !Rf_isNull(bound["exception"]);

  std::unordered_map<string, double> exceptionMap;
  if (has_exception) {
    NumericVector except_vec = bound["exception"];
    CharacterVector except_names = except_vec.names();
    for (int i = 0; i < except_vec.size(); ++i) {
      exceptionMap[ as<string>(except_names[i]) ] = except_vec[i];
    }
  }

  const int ncols = minmax_colnames.size();
  std::vector<BoundSpec> specs(ncols);

  for (int j = 0; j < ncols; ++j) {
    string var_name = as<string>(minmax_colnames[j]);
    int base_idx = pt.base_index_for(var_name);

    BoundSpec s;
    s.col_idx = base_idx;
    s.min_val = minmax(0, j);
    s.max_val = minmax(1, j);

    auto it = exceptionMap.find(var_name);
    if (it != exceptionMap.end()) {
      s.has_exception = true;
      s.exception_val = it->second;
    } else {
      s.has_exception = false;
      s.exception_val = NA_REAL;
    }
    specs[j] = s;
  }
  return specs;
}

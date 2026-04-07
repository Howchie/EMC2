#include "transform_utils.h"
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
  std::vector<TransformSpec> out;
  out.reserve(pt.base_names.size());

  for (int j = 0; j < pt.base_names.size(); ++j) {
    std::string nm = as<std::string>(pt.base_names[j]);
    for (int k = 0; k < full_names.size(); ++k) {
      if (nm == as<std::string>(full_names[k])) {
        TransformSpec s = full_specs[k];
        s.col_idx = j;
        out.push_back(s);
        break;
      }
    }
  }
  return out;
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
        col[i] = lw + range * R::pnorm(col[i], 0.0, 1.0, 1, 0);
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
  std::vector<TransformSpec> out;
  out.reserve(full_specs.size());

  for (const auto& sp : full_specs) {
    int base_idx = sp.col_idx;
    std::string nm = Rcpp::as<std::string>(pt.base_names[base_idx]);
    if (premap_set.find(nm) == premap_set.end()) out.push_back(sp);
  }
  return out;
}

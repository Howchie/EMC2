#include <RcppArmadillo.h>
// #include <RcppArmadillo.h>
#include "ParamTable.h"
#include "TrendEngine.h"   // contains TrendEngine definition
#include "transform_utils.h"
#include "wald_functions.h"
#include <unordered_map>
using namespace Rcpp;

// Lightweight wrapper for R-facing API: holds both plan and runtime
struct TrendEngine {
  TrendPlan plan;
  TrendRuntime runtime;

  TrendEngine(Rcpp::Nullable<Rcpp::List> trend_,
              const Rcpp::DataFrame& data_)
    : plan(trend_, data_), runtime(plan) {}
};

static NumericVector apply_pretransforms_cpp(NumericVector p_vector, List pretransforms) {
  if (pretransforms.size() == 0) return p_vector;
  if (!pretransforms.containsElementNamed("func")) return p_vector;

  CharacterVector func = pretransforms["func"];
  NumericVector lowervec = pretransforms.containsElementNamed("lower") ? as<NumericVector>(pretransforms["lower"]) : NumericVector();
  NumericVector uppervec = pretransforms.containsElementNamed("upper") ? as<NumericVector>(pretransforms["upper"]) : NumericVector();

  std::unordered_map<std::string, std::string> f_map;
  CharacterVector f_names = func.names();
  for (int i = 0; i < func.size(); ++i) {
    f_map[as<std::string>(f_names[i])] = as<std::string>(func[i]);
  }

  std::unordered_map<std::string, double> lower_map, upper_map;
  CharacterVector ln = lowervec.names();
  for (int i = 0; i < lowervec.size(); ++i) lower_map[as<std::string>(ln[i])] = lowervec[i];
  CharacterVector un = uppervec.names();
  for (int i = 0; i < uppervec.size(); ++i) upper_map[as<std::string>(un[i])] = uppervec[i];

  CharacterVector p_names = p_vector.names();
  for (int i = 0; i < p_vector.size(); ++i) {
    std::string nm = as<std::string>(p_names[i]);
    auto fit = f_map.find(nm);
    if (fit == f_map.end()) continue;

    const std::string f = fit->second;
    const double lw = (lower_map.find(nm) != lower_map.end()) ? lower_map[nm] : 0.0;
    const double up = (upper_map.find(nm) != upper_map.end()) ? upper_map[nm] : 1.0;
    const double val = p_vector[i];

    if (f == "exp") {
      p_vector[i] = lw + std::exp(val);
    } else if (f == "pnorm") {
      const double range = up - lw;
      p_vector[i] = lw + range * pnorm_std(val);
    }
  }
  return p_vector;
}

static NumericVector add_constants_cpp(const NumericVector& p_vector, const NumericVector& constants) {
  if (constants.size() == 1 && NumericVector::is_na(constants[0])) return clone(p_vector);

  NumericVector out(p_vector.size() + constants.size());
  std::copy(p_vector.begin(), p_vector.end(), out.begin());
  std::copy(constants.begin(), constants.end(), out.begin() + p_vector.size());

  CharacterVector p_names = p_vector.names();
  CharacterVector c_names = constants.names();
  CharacterVector all_names(out.size());
  std::copy(p_names.begin(), p_names.end(), all_names.begin());
  std::copy(c_names.begin(), c_names.end(), all_names.begin() + p_vector.size());
  out.names() = all_names;
  return out;
}


static NumericMatrix get_pars_matrix_oo_core(ParamTable& param_table,
                                             const Rcpp::List& designs,
                                             TrendRuntime* trend_runtime,
                                             const std::vector<TransformSpec>& full_specs,
                                             const Rcpp::CharacterVector& keep_names,
                                             bool return_kernel_matrix,
                                             bool return_all_pars,
                                             const std::vector<int>& kernel_output_codes) {
  if (trend_runtime) trend_runtime->reset_all_kernels();

  const int n_designs = designs.size();
  if (!trend_runtime) {
    if (return_kernel_matrix) {
      Rcpp::stop("return_kernel_matrix=TRUE requires a non-NULL trend");
    }

    // Common fast path: no trend operations, so we can map all designs in one pass
    // and apply full transforms directly without any parameter-set filtering.
    LogicalVector include_all(n_designs, true);
    param_table.map_from_designs(designs, include_all);
    if (!full_specs.empty()) {
      const auto split_set = param_table.split_transform_params();
      if (split_set.empty()) {
        c_do_transform_pt(param_table, full_specs);
      } else {
        auto postmap_specs = complement_specs_for_phases(param_table, full_specs, { &split_set });
        c_do_transform_pt(param_table, postmap_specs);
      }
    }

    if (return_all_pars) return param_table.materialize();
    return param_table.materialize_by_param_names(keep_names);
  }
  LogicalVector map_next(n_designs, false);
  std::unordered_set<std::string> transform_next;
  std::unordered_set<std::string> empty_set;

  const auto& premap_set = (trend_runtime ? trend_runtime->premap_trend_params() : empty_set);
  const auto& pretransform_set = (trend_runtime ? trend_runtime->pretransform_trend_params() : empty_set);

  if (trend_runtime && trend_runtime->has_premap()) {
    map_next = trend_runtime->premap_design_mask(designs);
    param_table.map_from_designs(designs, map_next);

    transform_next = premap_set;
    if (!transform_next.empty()) {
      auto specs_premap = filter_specs_by_param_set(param_table, full_specs, transform_next);
      c_do_transform_pt(param_table, specs_premap);
    }

    for (std::size_t i = 0; i < trend_runtime->premap_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->premap_ops[i], param_table);
    }
  }

  for (int i = 0; i < n_designs; ++i) {
    bool is_premap = (trend_runtime && trend_runtime->has_premap()) ? map_next[i] : false;
    map_next[i] = !is_premap;
  }
  param_table.map_from_designs(designs, map_next);

  if (trend_runtime && trend_runtime->has_pretransform()) {
    transform_next = pretransform_set;
    if (!transform_next.empty()) {
      auto specs_pretransform = filter_specs_by_param_set(param_table, full_specs, transform_next);
      c_do_transform_pt(param_table, specs_pretransform);
    }
    for (std::size_t i = 0; i < trend_runtime->pretransform_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->pretransform_ops[i], param_table);
    }
  }

  const auto split_set = param_table.split_transform_params();
  transform_next = param_names_excluding(param_table, { &premap_set, &pretransform_set });
  for (const auto& nm : split_set) transform_next.erase(nm);
  auto postmap_specs = filter_specs_by_param_set(param_table, full_specs, transform_next);
  c_do_transform_pt(param_table, postmap_specs);

  if (trend_runtime && trend_runtime->has_posttransform()) {
    for (std::size_t i = 0; i < trend_runtime->posttransform_ops.size(); ++i) {
      trend_runtime->apply_base_for_op(trend_runtime->posttransform_ops[i], param_table);
    }
  }

  if (return_kernel_matrix) {
    if (!trend_runtime) {
      Rcpp::stop("return_kernel_matrix=TRUE requires a non-NULL trend");
    }
    std::vector<int> codes = kernel_output_codes;
    if (codes.empty()) codes.push_back(1);
    return trend_runtime->all_kernel_outputs(param_table, codes);
  }

  if (return_all_pars) return param_table.materialize();
  return param_table.materialize_by_param_names(keep_names);
}

NumericMatrix get_pars_c_wrapper_oo_core(NumericMatrix particle_matrix,
                                         DataFrame data,
                                         NumericVector constants,
                                         List designs,
                                         List bounds,
                                         List transforms,
                                         List pretransforms,
                                         Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                         bool return_kernel_matrix = false,
                                         bool return_all_pars = false,
                                         IntegerVector kernel_output_codes = 1) {
  (void)bounds;
  const int n_trials = data.nrow();
  if (particle_matrix.nrow() == 0) stop("particle_matrix must have at least one row");

  NumericVector p_vector = particle_matrix(0, _);
  p_vector.attr("names") = colnames(particle_matrix);
  p_vector = apply_pretransforms_cpp(p_vector, pretransforms);
  p_vector = add_constants_cpp(p_vector, constants);

  ParamTable param_table = ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials, transforms);
  std::vector<TransformSpec> full_specs = make_transform_specs_for_paramtable(param_table, transforms);

  std::unique_ptr<TrendPlan> trend_plan;
  std::unique_ptr<TrendRuntime> trend_runtime;
  Rcpp::CharacterVector keep_names;

  if (!trend.isNull()) {
    trend_plan.reset(new TrendPlan(trend, data));
    trend_runtime.reset(new TrendRuntime(*trend_plan));
    trend_runtime->bind_all_ops_to_paramtable(param_table);

    Rcpp::CharacterVector dnames = designs.names();
    const auto& trend_params = trend_runtime->all_trend_params();
    keep_names = names_excluding(dnames, { &trend_params });
  } else {
    keep_names = designs.names();
  }

  TrendRuntime* tr_ptr = trend_runtime ? trend_runtime.get() : nullptr;
  std::vector<int> kernel_codes;
  kernel_codes.reserve(kernel_output_codes.size());
  for (int i = 0; i < kernel_output_codes.size(); ++i) kernel_codes.push_back(kernel_output_codes[i]);

  return get_pars_matrix_oo_core(param_table, designs, tr_ptr, full_specs, keep_names,
                                 return_kernel_matrix, return_all_pars, kernel_codes);
}

// Factory: create a ParamTable from p_types and n_trials
// [[Rcpp::export]]
SEXP ParamTable_create_from_p_types(int n_trials,
                                    CharacterVector p_types) {
  ParamTable* pt = new ParamTable(ParamTable::from_p_types(n_trials, p_types));
  XPtr<ParamTable> ptr(pt, true); // true => R will delete pt when GC'ed
  return ptr;
}

// Factory: create from an existing matrix + names (useful for debugging)
// [[Rcpp::export]]
SEXP ParamTable_create_from_matrix(NumericMatrix base,
                                   CharacterVector names) {
  ParamTable* pt = new ParamTable(base, names);
  XPtr<ParamTable> ptr(pt, true);
  return ptr;
}

// Wrapper: call materialize()
// [[Rcpp::export]]
NumericMatrix ParamTable_materialize(SEXP pt_xptr) {
  XPtr<ParamTable> pt(pt_xptr);
  return pt->materialize();
}

// // Wrapper: drop columns by name
// // [[Rcpp::export]]
// void ParamTable_drop(SEXP pt_xptr, CharacterVector drop_names) {
//   XPtr<ParamTable> pt(pt_xptr);
//   pt->drop(drop_names);
// }

// Wrapper: set a column by name
// [[Rcpp::export]]
void ParamTable_set_column(SEXP pt_xptr,
                           std::string name,
                           NumericVector col) {
  XPtr<ParamTable> pt(pt_xptr);
  pt->set_column_by_name(name, col);
}

// Wrapper
// [[Rcpp::export]]
SEXP ParamTable_create_from_pvector_designs(NumericVector p_vector,
                                            List designs,
                                            int n_trials) {
  ParamTable* pt = new ParamTable(
    ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials)
  );
  XPtr<ParamTable> ptr(pt, true);
  return ptr;
}


// [[Rcpp::export]]
void ParamTable_map_designs(SEXP pt_xptr,
                            Rcpp::List designs,
                            Rcpp::LogicalVector include_param) {
  Rcpp::XPtr<ParamTable> pt(pt_xptr);
  pt->map_from_designs(designs, include_param);
}

// // [[Rcpp::export]]
// SEXP TrendEngine_create(List trend, DataFrame data) {
//   // allocate on heap and wrap in XPtr so R manages lifetime
//   TrendEngine* engine = new TrendEngine(trend, data);
//   Rcpp::XPtr<TrendEngine> xp(engine, true); // true => delete when GC'd
//   return xp;
// }
// [[Rcpp::export]]
SEXP TrendEngine_create(SEXP trend, Rcpp::DataFrame data) {
  Rcpp::Nullable<Rcpp::List> ntrend(trend);
  TrendEngine* engine = new TrendEngine(ntrend, data);
  Rcpp::XPtr<TrendEngine> xp(engine, true);
  return xp;
}

// [[Rcpp::export]]
LogicalVector TrendEngine_premap_mask(SEXP engine_ptr, List designs) {
  Rcpp::XPtr<TrendEngine> engine(engine_ptr);
  return engine->plan.premap_design_mask(designs);
}


// [[Rcpp::export]]
int ParamTable_get_n_trials(SEXP pt_xptr) {
  XPtr<ParamTable> pt(pt_xptr);
  return pt->n_trials;   // or a getter if you made n_trials private
}


// [[Rcpp::export]]
void ParamTable_do_transform_all(SEXP param_table_ptr, Rcpp::List transform) {
  // however you're storing/accessing ParamTable; here assume by reference.
  ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);

  auto specs = make_transform_specs_for_paramtable(pt, transform);
  c_do_transform_pt(pt, specs);
}


// [[Rcpp::export]]
void ParamTable_do_transform_premap(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
  ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  // 1) Build full transform specs for all active parameters
  auto full_specs = make_transform_specs_for_paramtable(pt, transform);

  // 2) Filter specs to only premap trend parameters
  auto premap_specs = filter_specs_by_param_set(pt, full_specs,
                                                te.plan.premap_trend_params);

  // 3) Apply transforms in-place to those columns
  c_do_transform_pt(pt, premap_specs);
}

// [[Rcpp::export]]
void ParamTable_do_transform_pretransform(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
  ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  // 1) Build full transform specs for all active parameters
  auto full_specs = make_transform_specs_for_paramtable(pt, transform);

  // 2) Filter specs to only premap trend parameters
  auto pretransform_specs = filter_specs_by_param_set(pt, full_specs,
                                                te.plan.pretransform_trend_params);

  // 3) Apply transforms in-place to those columns
  c_do_transform_pt(pt, pretransform_specs);
}

// [[Rcpp::export]]
void ParamTable_do_transform_posttransform(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
  ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  // 1) Build full transform specs for all active parameters
  auto full_specs = make_transform_specs_for_paramtable(pt, transform);

  // 2) Filter specs to only premap trend parameters
  auto posttransform_specs = filter_specs_by_param_set(pt, full_specs,
                                                      te.plan.posttransform_trend_params);

  // 3) Apply transforms in-place to those columns
  c_do_transform_pt(pt, posttransform_specs);
}


// [[Rcpp::export]]
void ParamTable_do_transform_postmap(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
  ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  // 1) Build full transform specs for all active parameters
  auto full_specs = make_transform_specs_for_paramtable(pt, transform);

  // 2) Filter specs to only premap trend parameters
  auto postmap_specs = complement_specs_for_premap(pt, full_specs,
                                                te.plan.premap_trend_params);


  // 3) Apply transforms in-place to those columns
  c_do_transform_pt(pt, postmap_specs);
}


// [[Rcpp::export]]
void ParamTable_bind_trendops(SEXP param_table_ptr, SEXP trend_engine_ptr) {
  ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  te.runtime.bind_all_ops_to_paramtable(pt);
}


// // [[Rcpp::export]]
// Rcpp::List TrendEngine_run_premap_kernels_debug(SEXP param_table_ptr,
//                                                 SEXP trend_engine_ptr) {
//   ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   const std::size_t n_ops = te.runtime.premap_ops.size();
//   Rcpp::List out(n_ops);
//   Rcpp::CharacterVector nms(n_ops);
//
//   for (std::size_t i = 0; i < n_ops; ++i) {
//     TrendOpRuntime& op = te.runtime.premap_ops[i];
//
//     // ensure kernel is run
//     te.runtime.run_kernels_for_op(op, pt);
//
//     const std::vector<double>& traj = op.kernel_ptr->get_output();
//     Rcpp::NumericVector rtraj(traj.begin(), traj.end());
//
//     out[i] = rtraj;
//     nms[i] = op.spec->target_param;
//   }
//
//   out.attr("names") = nms;
//   return out;
// }


// [[Rcpp::export]]
void TrendEngine_apply_premap_bases(SEXP param_table_ptr,
                                    SEXP trend_engine_ptr) {
  ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  auto& runtime = te.runtime;
  const std::size_t n_ops = runtime.premap_ops.size();

  for (std::size_t i = 0; i < n_ops; ++i) {
    TrendOpRuntime& op = runtime.premap_ops[i];
    runtime.apply_base_for_op(op, pt);
  }
}

// [[Rcpp::export]]
void TrendEngine_apply_pretransform_bases(SEXP param_table_ptr,
                                    SEXP trend_engine_ptr) {
  ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  auto& runtime = te.runtime;
  const std::size_t n_ops = runtime.pretransform_ops.size();

  Rcpp::Rcout << "has_pretransform = " << runtime.has_pretransform()
              << ", n_pretransform_ops = " << runtime.pretransform_ops.size() << "\n";

  for (std::size_t i = 0; i < n_ops; ++i) {
    TrendOpRuntime& op = runtime.pretransform_ops[i];
    const TrendOpSpec& spec = *op.spec;
    Rcpp::Rcout << "Pretransform op " << i
                << " target_param=" << spec.target_param
                << " base_type=" << spec.base_type
                << "\n";

    runtime.apply_base_for_op(op, pt);
  }
}

// [[Rcpp::export]]
void TrendEngine_apply_posttransform_bases(SEXP param_table_ptr,
                                          SEXP trend_engine_ptr) {
  ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
  TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);

  auto& runtime = te.runtime;
  const std::size_t n_ops = runtime.posttransform_ops.size();

  for (std::size_t i = 0; i < n_ops; ++i) {
    TrendOpRuntime& op = runtime.posttransform_ops[i];
    runtime.apply_base_for_op(op, pt);
  }
}






// // Helper: all parameter names in pt that are NOT in premap_trend_params
// std::unordered_set<std::string>
//   make_non_premap_param_set(const ParamTable& pt,
//                             const std::unordered_set<std::string>& premap_trend_params)
//   {
//     std::unordered_set<std::string> out;
//     const int p = pt.base_names.size();
//     out.reserve(p);
//
//     for (int j = 0; j < p; ++j) {
//       std::string nm = Rcpp::as<std::string>(pt.base_names[j]);
//       if (premap_trend_params.find(nm) == premap_trend_params.end()) {
//         out.insert(nm);
//       }
//     }
//     return out;
//   }

// // [[Rcpp::export]]
// Rcpp::NumericMatrix ParamTable_run_pipeline(Rcpp::NumericVector p_vector,
//                                             Rcpp::List designs,
//                                             int n_trials,
//                                             Rcpp::List trend,
//                                             Rcpp::DataFrame data,
//                                             Rcpp::List transform)
// {
//   using namespace Rcpp;
//
//   // 1) Create ParamTable from p_vector + designs
//   ParamTable pt = ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials);
//
//   // 2) Create TrendEngine
//   TrendEngine te(trend, data);
//   te.runtime.bind_all_ops_to_paramtable(pt);
//
//   // 3) Create TransformSpecs
//   auto full_specs = make_transform_specs_for_paramtable(pt, transform);
//
//   // 4) Mask for premap trend *design* entries
//   LogicalVector include_params;
//   if (te.runtime.has_premap()) {
//     include_params = te.plan.premap_design_mask(designs);
//   } else {
//     include_params = LogicalVector(designs.size(), false);
//   }
//
//   if (te.runtime.has_premap()) {
//     // 5) Map these premap-trend parameters
//     pt.map_from_designs(designs, include_params);
//
//     // 6) Filter specs to only premap trend parameters
//     auto premap_specs = filter_specs_by_param_set(pt, full_specs,
//                                                   te.plan.premap_trend_params);
//
//     // 7) Apply transforms in-place to those premap trend columns
//     c_do_transform_pt(pt, premap_specs);
//
//     // 9) Run premap kernels and apply bases (in-place on pt)
//     std::size_t n_ops = te.runtime.premap_ops.size();
//     for (std::size_t i = 0; i < n_ops; ++i) {
//       TrendOpRuntime& op = te.runtime.premap_ops[i];
//       te.runtime.apply_base_for_op(op, pt);  // will run kernel if needed and apply base in-place
//     }
//   }
//
//   // 10) Map designs for non-trend parameters (complement of include_params)
//   const int n_designs = designs.size();
//   LogicalVector include_nontrend(n_designs);
//   for (int i = 0; i < n_designs; ++i) {
//     include_nontrend[i] = !include_params[i];
//   }
//   pt.map_from_designs(designs, include_nontrend);
//
//   // 11) Apply pretransform trends
//   if (te.plan.has_pretransform()) {
//     std::size_t n_ops = te.runtime.pretransform_ops.size();
//     for (std::size_t i = 0; i < n_ops; ++i) {
//       TrendOpRuntime& op = te.runtime.pretransform_ops[i];
//       te.runtime.apply_base_for_op(op, pt);  // will run kernel if needed and apply base in-place
//     }
//   }
//
//   // 12) Apply transforms for non-trend parameters
//   auto nontrend_params = make_non_premap_param_set(pt, te.plan.premap_trend_params);
//   auto postmap_specs   = filter_specs_by_param_set(pt, full_specs, nontrend_params);
//   c_do_transform_pt(pt, postmap_specs);
//
//   // 13) Apply posttransform trends
//   if (te.plan.has_posttransform()) {
//     std::size_t n_ops = te.runtime.posttransform_ops.size();
//     for (std::size_t i = 0; i < n_ops; ++i) {
//       TrendOpRuntime& op = te.runtime.posttransform_ops[i];
//       te.runtime.apply_base_for_op(op, pt);  // will run kernel if needed and apply base in-place
//     }
//   }
//
//
//   // Only parameters from names(designs), but drop all trend_pnames
//   CharacterVector dnames = designs.names();
//   std::vector<std::string> keep;
//   keep.reserve(dnames.size());
//
//   for (int i = 0; i < dnames.size(); ++i) {
//     std::string nm = Rcpp::as<std::string>(dnames[i]);
//     if (te.plan.all_trend_params.find(nm) == te.plan.all_trend_params.end()) {
//       keep.push_back(nm);
//     }
//   }
//
//   CharacterVector keep_names(keep.size());
//   for (int i = 0; i < (int)keep.size(); ++i) {
//     keep_names[i] = keep[i];
//   }
//
//   return pt.materialize_by_param_names(keep_names);
//   // 14) Return final materialized matrix
//   // return pt.materialize();
// }

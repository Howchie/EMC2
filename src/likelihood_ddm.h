#ifndef EMC2_LIKELIHOOD_DDM_H
#define EMC2_LIKELIHOOD_DDM_H

// ---------------------------------------------------------------------------
// DDM/BOU likelihood kernels, extracted from particle_ll.cpp.
//
// Declares the entry points particle_ll.cpp's calc_ll / calc_ll_oo /
// calc_ll_oo_pw paths call; definitions (bodies moved verbatim) live in
// likelihood_ddm.cpp.  The shared per-data state these kernels operate on is
// owned by likelihood_shared.h.
//
// Contracts stay where they live:
//   * DDMAdapter / ContextForDDMModels come from race_contract.h and are not
//     duplicated here.
//   * ddm_wien_adapter() keeps its function-local static adapter, so the
//     function-pointer storage stays stable for every existing caller.
//   * bou_bnd_kind_from_type() keeps its suffix precedence (BWEIB > BEXP >
//     BLIN_MULT > BLIN_ADD > fixed).
//   * Default arguments live on these declarations only; the definitions in
//     likelihood_ddm.cpp do not repeat them.
// ---------------------------------------------------------------------------

#include "likelihood_shared.h"
#include "race_contract.h"       // DDMAdapter / ContextForDDMModels

#include <string>
#include <vector>

struct ParamTable;                      // ParamTable.h
namespace emc2col { struct ColSpec; }   // col_registry.h

// True when every trial is finite and untruncated/unconsored on the DDM path.
bool ddm_data_all_finite_untruncated(const Rcpp::DataFrame& data,
                                     const int n_trials);

// The adapter every existing DDM caller gets when it passes none: the Wiener
// raw-kernel pair plus the DDM column spec.
const DDMAdapter& ddm_wien_adapter();

// Collapsing-bound variant of the bounded OU, selected by
// BOU(boundary_collapse=).
int bou_bnd_kind_from_type(const std::string& type_std);

// True for the legacy BOU start-point anchor; the default BOU model uses the
// midpoint anchor so start-point variability can be seeded in one march.
bool bou_anchor_at_z_from_type(const std::string& type_std);

// Raw-buffer variant for DDM to skip materialization and allocations.
// `ker` supplies the two model-specific operations and defaults to the Wiener
// pair, so every existing call site is unchanged in behaviour.
double c_log_likelihood_DDM_pt(const double* const* cols,
                               const double* rt_ptr,
                               const int* R_ptr,
                               const int n_trials,
                               const int* expand_ptr,
                               const int n_out,
                               double min_ll,
                               const int* is_ok,
                               bool gng,
                               bool all_finite_untruncated,
                               ModelSharedState* shared,
                               Rcpp::NumericVector* trial_ll_out = nullptr,
                               const DDMAdapter* ker_in = nullptr);

// Compatibility wrapper for old calc_ll path (materialized parameter matrix).
double c_log_likelihood_DDM(Rcpp::NumericMatrix pars, Rcpp::DataFrame data,
                            const int n_trials, Rcpp::IntegerVector expand,
                            double min_ll, Rcpp::LogicalVector is_ok, bool gng,
                            bool all_finite_untruncated = false,
                            Rcpp::NumericVector* trial_ll_out = nullptr);

// DDM per-data shared state + canonical base-column pointers for the raw
// c_log_likelihood_DDM_pt kernel. Returns false (raw path unusable) when any
// canonical DDM parameter is missing from the table.
bool init_ddm_shared_state(Rcpp::DataFrame data, int n_trials,
                           const ParamTable& table,
                           ModelSharedState& shared,
                           std::vector<const double*>& cols,
                           const emc2col::ColSpec* spec_in = nullptr,
                           const Rcpp::CharacterVector* keep_names_in = nullptr);

#endif // EMC2_LIKELIHOOD_DDM_H

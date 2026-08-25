#ifndef EMC2_MARGINAL_LIKELIHOOD_H
#define EMC2_MARGINAL_LIKELIHOOD_H

#include <Rcpp.h>

// ---------------------------------------------------------------------------
// Coherent t0 marginalization, extracted verbatim from particle_ll.cpp.
//
// Integrates ONE shared subject-level parameter (the race/GNG non-decision
// time t0) out of the COMPLETE subject likelihood by a per-particle
// Laplace-centred composite Gauss-Legendre rule on the log-t0 (sampled)
// axis:
//     Lbar(theta) = int p(t0|eta) L(t0, theta) dt0,   log t0 ~ N(mu, sigma)
//
// The registered likelihood kernel is reused as a black box: node values
// overwrite the sampled t0 column and re-enter calc_ll_oo with marginalise
// disabled, so the machinery stays model-agnostic (no per-model code).
//
// The internal helpers (MarginalGrid/MarginalRule, the marginal_* probe and
// fit routines, and calc_ll_oo_marginal_core) live in marginal_likelihood.cpp
// with internal linkage; this header exposes only the two entry points below.
// ---------------------------------------------------------------------------

// Black-box evaluator the marginalization recurses into; defined in
// particle_ll.cpp.  Exact signature of its [[Rcpp::export]]ed definition.
Rcpp::NumericVector calc_ll_oo(Rcpp::NumericMatrix particle_matrix,
                               Rcpp::DataFrame data,
                               Rcpp::NumericVector constants,
                               Rcpp::List designs, Rcpp::String type,
                               Rcpp::List bounds, Rcpp::List transforms,
                               Rcpp::List pretransforms,
                               Rcpp::CharacterVector p_types, double min_ll,
                               Rcpp::Nullable<Rcpp::List> trend,
                               Rcpp::Nullable<Rcpp::List> marginalise);

// Marginal log-likelihood: log-sum-exp of the shared core's per-node terms.
// Not Rcpp-exported; external linkage so calc_ll_oo's delegation branch
// (particle_ll.cpp) can reach it when `marginalise` is supplied.
Rcpp::NumericVector calc_ll_oo_marginal(
    Rcpp::NumericMatrix particle_matrix, Rcpp::DataFrame data,
    Rcpp::NumericVector constants, Rcpp::List designs, Rcpp::String type,
    Rcpp::List bounds, Rcpp::List transforms, Rcpp::List pretransforms,
    Rcpp::CharacterVector p_types, double min_ll,
    Rcpp::Nullable<Rcpp::List> trend, Rcpp::List marginalise);

// Storage-time accessor for the t0 reconstruction; [[Rcpp::export]] remains
// on the definition in marginal_likelihood.cpp (with its trend = R_NilValue
// default), so codegen and registration are unchanged.
Rcpp::List calc_ll_oo_marginal_nodes(
    Rcpp::NumericMatrix particle_matrix, Rcpp::DataFrame data,
    Rcpp::NumericVector constants, Rcpp::List designs, Rcpp::String type,
    Rcpp::List bounds, Rcpp::List transforms, Rcpp::List pretransforms,
    Rcpp::CharacterVector p_types, double min_ll, Rcpp::List marginalise,
    Rcpp::Nullable<Rcpp::List> trend);

#endif

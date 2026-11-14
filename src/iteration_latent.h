// The latent-factor hierarchies' group state as the iteration's steps pass it
// to and from their group draws (src/gibbs_factor.cpp). `kind` is 0 for a
// factor hierarchy, 1 for the shrinkage-prior one, 2 for a structural
// equation model.
#ifndef EMC2_ITERATION_LATENT_H
#define EMC2_ITERATION_LATENT_H

#include <Rcpp.h>

namespace emc_iter {

// Thrown where a state is not in the form the native draws read.
struct LatentDecline {};

// The draw's starting state from a group state named as the draws return
// it (tmu, eta, lambda_untransf, ...), for n subjects, p parameters, k
// factors and n_cov covariates.
Rcpp::List latent_state(SEXP pars, int kind, int n, int p, int k, int n_cov);

// The draw's result in the form the hierarchy's own R draw returns it:
// `names` the parameters' names, `delta_form` the form the store gives the
// factors' shrinkage, `eta_rows` the subjects' names on the scores (NULL for
// none) and `alpha` the subjects' states.
Rcpp::List latent_pars(const Rcpp::List& state, int kind, SEXP names, SEXP delta_form,
                       SEXP eta_rows, SEXP alpha);

}  // namespace emc_iter

#endif

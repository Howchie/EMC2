// Access to the iteration loop's state (src/iteration_loop.cpp): R lists read
// by name and rebuilt as `x$name <- v` rebuilds them, the state environment's
// variables, and the clock the R steps read.
#ifndef EMC2_ITERATION_STATE_H
#define EMC2_ITERATION_STATE_H

#include <Rcpp.h>
#include <vector>

namespace emc_iter {

// x$name of a list; NULL when absent.
SEXP lget(SEXP x, const char* name);
// The list `x$name <- v` leaves (v not NULL): a shallow copy with the element
// replaced, or appended under that name.
SEXP lset(SEXP x, const char* name, SEXP v);
// A variable of the state environment; NULL when unset.
SEXP var(SEXP env, const char* name);
void set_var(SEXP env, const char* name, SEXP v);
// A scalar's value, `if_null` for NULL.
double num(SEXP x, double if_null);
// isTRUE().
bool flag(SEXP x);

// proc.time(): elapsed, and the process's own user plus system time.
struct Clock {
  double elapsed, cpu;
};
Clock clock_now();
// .emc_reserved_cost().
double reserved_cost(double elapsed, double cores);
// The store index of the block's iteration i, typed as start_iter + i.
SEXP iteration_index(SEXP cx, int i);

// The session's random stream (.Random.seed of the global environment).
SEXP global_seed();
void set_global_seed(SEXP seed);
// One uniform from the stream `seed`, as stats::runif(1) draws it there;
// `next` is the stream after the draw.
double stream_unif(SEXP seed, Rcpp::RObject& next);
// proposals[seq_len(p), , drop = FALSE] of a (p + 1)-row proposals matrix.
SEXP proposal_rows(SEXP proposals, int p);
// Whether a logical option is unset or TRUE.
bool option_on(const char* name);
// A function of the package's namespace, which the loop sets from the
// environment of a package function it is handed.
void set_package_namespace(SEXP ns);
Rcpp::Function package_function(const char* name);
// A variable of the package's namespace.
SEXP package_value(const char* name);

// A group move's context (.emc_shear_context(), .emc_collapsed_context()):
// the value last built for the slot while the sampler fields it was built
// from are unchanged, else NULL (the R function then builds it).
SEXP move_context(const char* slot, SEXP sampler);

// The native steps: TRUE when the step ran, FALSE when it left the state
// untouched for the R step.
bool step_group(SEXP st, SEXP cx, int i);
bool step_particle(SEXP st, SEXP cx, int i);
bool step_ensemble(SEXP st, SEXP cx, int i);
bool step_shear(SEXP st, SEXP cx, int i);
bool step_collapsed(SEXP st, SEXP cx, int i);
bool step_pooled(SEXP st, SEXP cx, int i);
bool step_hmc(SEXP st, SEXP cx, int i);

// .chol_factor() of a covariance matrix: its upper Cholesky factor, the
// inverse factor and the sum of the inverse's log diagonal. NULL when the
// matrix is not positive definite (R repairs it).
SEXP chol_factor(SEXP sigma);

// .emc_shear_draw_a_half(): the mixing scales given the covariance, drawn
// from the session's stream. False when the covariance's blocks do not
// factorise.
bool draw_a_half(SEXP V, SEXP spec, std::vector<double>& out);
// .emc_shear_log_prior_mu(): the group mean's log prior density up to its
// constant, with R's products. False when R would take its non-finite path.
bool log_prior_mu(SEXP mu, SEXP prior, double& out);

}  // namespace emc_iter

#endif

// The subjects' likelihoods of candidate rows for a group move's round, as
// .emc_wpool_subject_ll() (R/chain_pool.R) gives them: on the
// stored-argument path in this process where it serves every subject, else
// through that function and so the worker pool.
#ifndef EMC2_SUBJECT_LL_H
#define EMC2_SUBJECT_LL_H

#include <Rcpp.h>

struct SubjectLl {
  // The worker pool as the last round left it.
  Rcpp::RObject pool;
  // The context .emc_wpool_subject_ll() reads (the sampler), the subjects'
  // partition over the workers and the move's name.
  SEXP sampler, part;
  const char* source;
  // Whether the stored-argument path may serve a round in this process.
  bool in_process;
  // The function itself.
  SEXP dispatch;
};

// Whether .emc_wpool_subject_ll() would try the stored-argument path in this
// process for this pool and sampler.
bool emc_subject_ll_in_process(SEXP pool, SEXP sampler);

// rows[[s]]'s likelihoods under subject s's data, a column per candidate row.
// Without a source only the stored-argument path is tried, and FALSE says it
// did not serve every subject. The dispatch runs R code: the session's
// random stream is written back before it and read again after.
bool emc_subject_lls(const Rcpp::List& rows, const Rcpp::List& data, SEXP model, int n_cand,
                     Rcpp::NumericMatrix& ll, SubjectLl* src);

#endif

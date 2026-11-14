// The chain's side of an exchange with its worker pool (R/chain_pool.R): a
// round's requests out, the workers' replies back and checked, with the
// subjects no valid reply holds computed in this process. A pool that fails
// part-way is marked lost exactly as the R functions mark it.
#ifndef EMC2_POOL_EXCHANGE_H
#define EMC2_POOL_EXCHANGE_H

#include <Rcpp.h>

namespace emc_pool {

// .emc_lpt_partition(): the subjects over the workers, costliest first to
// the least loaded, with that dispatch order as the `order` attribute.
SEXP lpt_partition(const double* cost, int n, int n_workers);

// .emc_cores_now(): the cores the chain may use at this moment.
double cores_now(SEXP core_ctl, double base);

// .emc_wpool_subject_ll() over the shared subject queue: rows[[s]]'s
// likelihoods under subject s's data into the n x k matrix `ll`. FALSE,
// before anything is sent, when the pool does not serve the queue; the R
// function then handles the round.
bool subject_ll_queue(Rcpp::RObject& pool, SEXP rows, SEXP ctx, SEXP part, const char* source,
                      int k, double* ll);

// .emc_wpool_iter() over the shared subject queue.
struct ParticleRequest {
  SEXP ctx, part, alpha, population_mu, group_chol, pm_settings, prev_ll, seeds, extra,
      pooled, skip;
  int n_pars;
};
struct ParticleReply {
  Rcpp::RObject props, pm_settings, seeds, times, extra_ll, rejects, ensemble_rejects;
  double worker_cpu;
};
// FALSE, before anything is sent, when the pool does not serve the queue.
bool particle_queue(Rcpp::RObject& pool, const ParticleRequest& request, ParticleReply& reply);

// msgs[[w]] to worker w for each message present, and each worker's reply
// (NULL where it gave none).
Rcpp::List requests(Rcpp::RObject& pool, const Rcpp::List& msgs);

// The pool marked lost for the reason given.
void lost(Rcpp::RObject& pool, const char* why);

}  // namespace emc_pool

#endif

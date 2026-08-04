# ===========================================================================
# Flattened (chain x subject) particle scheduling on a persistent worker pool
# ===========================================================================
#
# The original design nests two levels of parallelism and puts a barrier at
# each:
#
#   auto_mclapply over chains            (cores_for_chains)
#     run_stage: for each iteration
#       mcmapply over subjects           (cores_per_chain)
#
# That has three costs.  (1) The inner mcmapply forks a fresh set of workers on
# *every* MCMC iteration.  (2) mc.preschedule splits subjects into contiguous
# equal-size chunks, so an uneven subject count or uneven per-subject cost
# leaves cores idle at the end of every iteration.  (3) The core budget of a
# chain is fixed when it is forked, so a chain that finishes its block early
# cannot hand its cores to a chain that is still running.
#
# Here the two levels are flattened into one: chains advance in lockstep, and
# the n_chains * n_subjects particle updates of a single iteration form one flat
# task list dispatched over a single, persistent, load-balanced pool.  Chains
# are independent given the group-level parameters, so this changes nothing
# statistically -- only who computes what, and when.
#
# Consequences worth knowing:
#   * Cores are never idle waiting on one straggler chain; a chain that is
#     cheap this iteration simply contributes fewer tasks to the shared queue.
#   * The Gibbs steps become serial across chains.  They are O(n_pars^3) matrix
#     draws, far cheaper than n_subjects likelihood evaluations, so this is a
#     good trade.
#   * Each (chain, subject) carries its own L'Ecuyer RNG stream, so results are
#     independent of how tasks happen to be scheduled and of how many cores are
#     used.  The nested-fork scheme did not have that property: the same fit on
#     1 and on 2 cores gave different draws.
#
# MEASURED RESULT: this loses to the nested path, and is OFF by default.
#
# On a 24-subject (150-450 trials each), 3-chain, 5-parameter RDM fit, 40 burn
# iterations, nested vs flat elapsed seconds:
#
#     cores    3      6     12        3      6     12
#     pf       30     30    30        100    100   100
#     nested   2.3    3.3   3.2       4.7    4.8   4.1
#     flat     4.9    4.6   4.5       7.9    6.0   5.3
#
# The reason is structural, and it is worth recording so this is not
# re-attempted blind.  With cores_for_chains = n_chains the nested path forks
# once per *block* and each chain then runs its whole 40-iteration block with
# no further communication at all.  Any flattened scheme must instead
# synchronise every iteration, because the Gibbs step needs all of a chain's
# subjects.  An iteration is only ~45-160 ms of per-chain work, which does not
# amortise a per-iteration fork or dispatch.
#
# Note the nested path does not scale past cores = n_chains either (3.3 s on 6
# cores, 3.2 s on 12, against 2.3 s on 3).  That is the same barrier seen from
# the other side: extra cores within a chain have nothing profitable to do at
# this granularity.  Idle cores here are not a scheduling failure that better
# dispatch can fix.
#
# Where this could still pay off, and why it is kept: models whose per-subject
# likelihood is expensive enough to dwarf the per-iteration barrier (the
# PDE-backed models, or very large trial counts), and any case where subject
# costs are wildly uneven.  Both were out of reach to test here.
#
# Turn it on with options(emc2.flat_parallel = TRUE).
# It requires RNGkind("L'Ecuyer-CMRG"), which it sets once if needed.

# Package-level home for worker state.  Fork workers inherit it for free, so the
# per-subject data is never serialised on Unix; PSOCK workers get one explicit
# copy per block.
.emc_runtime <- new.env(parent = emptyenv())

.emc_use_flat_parallel <- function() {
  isTRUE(getOption("emc2.flat_parallel", FALSE))
}

# --- pool lifecycle --------------------------------------------------------

# `static` is everything that is constant for a whole block.  It is installed
# into .emc_runtime before any worker starts, so forked workers inherit it
# rather than receiving it down a socket.
#
# Transport choice matters more than it looks.  A persistent socket cluster
# was the obvious design -- workers stay alive, nothing is re-forked -- but it
# loses badly here, because the barrier is per MCMC *iteration*: an iteration
# of this size is only ~130 ms of work, while every socket round trip costs
# ~30-40 ms (the delayed-ACK stall, measured both on clusterCall and on
# clusterApplyLB).  Three round trips per iteration is already more latency
# than there is work.  Measured on a 24-subject, 3-chain RDM fit at 3 cores:
# 14.5 s with a fork cluster and 3 chunks/worker, 5.5 s at 1 chunk/worker,
# against 2.4 s for the nested mclapply path.
#
# So on Unix the flattening is done over mclapply instead: no sockets, no
# persistent pool, and the per-iteration fork is a cost the nested path already
# pays (in fact it pays more of them -- one set of forks per chain).  What is
# kept is the part that actually matters: one flat, longest-first task list
# over all (chain, subject) pairs, so a chain with cheap subjects this
# iteration cannot strand cores.  Windows has no fork and falls back to the
# socket cluster.
.emc_pool_start <- function(n_workers, static) {
  .emc_runtime$static <- static
  .emc_runtime$iter <- NULL
  if (n_workers <= 1) return(NULL)

  if (Sys.info()[["sysname"]] != "Windows") {
    return(structure(list(n = as.integer(n_workers)), class = "emc_fork_pool"))
  }
  tryCatch({
    cl <- parallel::makePSOCKcluster(n_workers)
    parallel::clusterEvalQ(cl, suppressMessages(requireNamespace("EMC2")))
    parallel::clusterCall(cl, function(s) {
      assign("static", s, envir = get(".emc_runtime", envir = asNamespace("EMC2")))
      NULL
    }, static)
    cl
  }, error = function(e) NULL)
}

.emc_pool_stop <- function(pool) {
  if (!is.null(pool) && !inherits(pool, "emc_fork_pool")) {
    try(parallel::stopCluster(pool), silent = TRUE)
  }
  .emc_runtime$static <- NULL
  .emc_runtime$iter <- NULL
  invisible(NULL)
}

# State that changes every iteration (the group-level parameters each chain's
# subjects are conditioned on).  Under mclapply the forks simply inherit it, so
# this costs nothing; under the Windows socket cluster it rides along with the
# work chunks instead of being broadcast, because a clusterCall to every worker
# measured ~43 ms for a 7 kB payload.
.emc_pool_set_iter <- function(pool, iter_state) {
  .emc_runtime$iter <- iter_state
  invisible(NULL)
}

.emc_pool_apply <- function(pool, tasks, FUN, n_workers = 1L, iter_state = NULL) {
  if (is.null(pool) || length(tasks) == 0) {
    # Running in the master: each task installs its own RNG stream into
    # .Random.seed, which would otherwise leave the master's own stream (used
    # by the Gibbs steps and by the next block's stream derivation) sitting
    # wherever the last task left it.  Restoring it is what makes the serial
    # path agree with the pooled one.
    had_seed <- exists(".Random.seed", envir = globalenv())
    if (had_seed) saved <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) assign(".Random.seed", saved, envir = globalenv())
    }, add = TRUE)
    return(lapply(tasks, FUN))
  }
  if (inherits(pool, "emc_fork_pool")) {
    # Tasks arrive longest-first; mclapply's prescheduled split then deals them
    # round-robin, which on an LPT-ordered list is already well balanced.  The
    # per-task RNG stream is explicit, so mc.set.seed plays no part.
    return(parallel::mclapply(tasks, FUN, mc.cores = pool$n))
  }
  chunks <- lapply(.emc_chunk_tasks(tasks, n_workers),
                   function(ch) list(tasks = ch, iter = iter_state))
  unlist(parallel::clusterApplyLB(pool, chunks, .emc_run_chunk, FUN),
         recursive = FALSE)
}

.emc_run_chunk <- function(chunk, FUN) {
  .emc_runtime$iter <- chunk$iter
  lapply(chunk$tasks, FUN)
}

# Dispatching one (chain, subject) at a time costs a socket round trip per task,
# which is only worth paying when a task is much more expensive than the trip.
# Grouping the queue into a few chunks per worker amortises that while keeping
# most of the balancing: the task list arrives in longest-first order, so
# dealing it round-robin leaves every chunk with a comparable mix of big and
# small subjects, and there are still several chunks per worker for the
# load-balancer to even out.
.emc_chunk_tasks <- function(tasks, n_workers,
                             chunks_per_worker = getOption("emc2.chunks_per_worker", 3L)) {
  n <- length(tasks)
  if (n_workers <= 1L || n <= 1L) return(list(tasks))
  n_chunks <- min(n, as.integer(n_workers) * chunks_per_worker)
  unname(split(tasks, rep_len(seq_len(n_chunks), n)))
}

# --- worker entry points ---------------------------------------------------
#
# These are package-level functions on purpose: serialising a closure would drag
# its enclosing environment to every worker on every task, whereas a namespace
# function is sent by reference.

.emc_particle_task <- function(task) {
  st <- .emc_runtime$static
  it <- .emc_runtime$iter
  cc <- task$chain
  s <- task$subj
  assign(".Random.seed", task$seed, envir = globalenv())
  out <- safe_new_particle(
    s = s,
    data = st$data[[s]],
    pm_settings = task$pm_settings,
    eff_mu = st$eff_mu[[cc]][[s]], eff_var = st$eff_var[[cc]][[s]],
    chains_mu = st$chains_mu[[cc]][[s]], chains_var = st$chains_var[[cc]][[s]],
    prev_ll = task$prev_ll,
    parameters = it$pars[[cc]],
    model = st$model, stage = st$stage, type = st$type, tune = st$tune,
    marginalise = st$marginalise, r_cores = st$r_cores,
    chol_cache = st$chol_caches[[cc]][[s]],
    group_chol = it$group_chol[[cc]]
  )
  list(res = out, seed = get(".Random.seed", envir = globalenv()),
       chain = cc, subj = s)
}

.emc_start_task <- function(task) {
  st <- .emc_runtime$static
  it <- .emc_runtime$iter
  cc <- task$chain
  assign(".Random.seed", task$seed, envir = globalenv())
  out <- start_proposals(s = task$subj, parameters = it$startpoints[[cc]],
                         n_particles = st$particles,
                         pmwgs = st$samplers[[cc]], type = st$type,
                         r_cores = st$r_cores)
  list(res = out, seed = get(".Random.seed", envir = globalenv()),
       chain = cc, subj = task$subj)
}

# --- RNG streams -----------------------------------------------------------
#
# One independent L'Ecuyer stream per (chain, subject), advanced by the worker
# and stored back in the master.  This is what makes the sampler's output
# independent of the number of cores and of the scheduling order.

.emc_init_streams <- function(n_chains, n_subjects) {
  if (!identical(RNGkind()[1], "L'Ecuyer-CMRG")) {
    RNGkind("L'Ecuyer-CMRG")
  }
  if (!exists(".Random.seed", envir = globalenv())) stats::runif(1)
  seed <- get(".Random.seed", envir = globalenv())
  streams <- vector("list", n_chains)
  for (cc in seq_len(n_chains)) {
    streams[[cc]] <- vector("list", n_subjects)
    for (s in seq_len(n_subjects)) {
      seed <- parallel::nextRNGStream(seed)
      streams[[cc]][[s]] <- seed
    }
  }
  # Leave the master's own stream past every stream handed out.
  assign(".Random.seed", parallel::nextRNGStream(seed), envir = globalenv())
  streams
}

# --- task ordering ---------------------------------------------------------
#
# Longest-processing-time-first.  With load-balanced dispatch this is the
# classic 4/3-approximation to optimal makespan; dispatching in subject order
# instead lets a big subject start last and strand the whole iteration on it.
# Cost is proxied by the number of (compressed) dadm rows, which is what the
# likelihood loops over.
.emc_subject_order <- function(data) {
  cost <- vapply(data, function(d) {
    if (is.data.frame(d)) nrow(d) else sum(vapply(d, function(x)
      if (is.data.frame(x)) nrow(x) else 0L, numeric(1)))
  }, numeric(1))
  order(cost, decreasing = TRUE)
}

# --- flattened stage runner -------------------------------------------------

# Drop-in replacement for `auto_mclapply(emc_list, run_stages, ...)`: runs every
# chain of `emc_list` through one stage, scheduling all chains' subjects on a
# single shared pool.  Returns the list of updated samplers.
run_stages_flat <- function(emc_list, stage = "preburn", iter = 0,
                            verbose = TRUE, verboseProgress = TRUE,
                            particle_factor = 50, search_width = NULL,
                            n_workers = 1, r_cores = 1) {
  n_chains <- length(emc_list)
  particles <- round(particle_factor * sqrt(emc_list[[1]]$n_pars))

  # Generalisation of .particle_core_budget() to the flattened queue: with
  # fewer (chain, subject) tasks than workers the surplus workers would simply
  # sit idle, so hand them to the per-likelihood particle split instead.
  n_tasks <- n_chains * emc_list[[1]]$n_subjects
  if (n_workers > n_tasks && n_tasks >= 1L) {
    r_cores <- max(r_cores, n_workers %/% n_tasks)
    n_workers <- n_tasks
  }

  needs_init <- !vapply(emc_list, function(x) isTRUE(x$init), logical(1))
  if (any(needs_init)) {
    emc_list <- init_flat(emc_list, which(needs_init), particles = particles,
                          n_workers = n_workers, r_cores = r_cores)
  }
  if (iter == 0) return(emc_list)

  run_stage_flat(emc_list, stage = stage, iter = iter, particles = particles,
                 tune = list(search_width = search_width),
                 n_workers = n_workers, r_cores = r_cores,
                 verbose = verbose, verboseProgress = verboseProgress)
}

# Flattened equivalent of init(): start points for every (chain, subject) pair
# are drawn on one shared queue instead of one mclapply per chain.
init_flat <- function(emc_list, chain_idx, particles, n_workers, r_cores) {
  n_subjects <- emc_list[[1]]$n_subjects
  streams <- .emc_init_streams(length(emc_list), n_subjects)

  startpoints <- vector("list", length(emc_list))
  startpoints_comb <- vector("list", length(emc_list))
  for (cc in chain_idx) {
    pmwgs <- emc_list[[cc]]
    sp <- get_startpoints(pmwgs, start_mu = NULL, start_var = NULL, type = pmwgs$type)
    startpoints[[cc]] <- sp
    startpoints_comb[[cc]] <- sp
    if (any(pmwgs$nuisance)) {
      type_nuis <- pmwgs$sampler_nuis$type
      sp_nuis <- get_startpoints(pmwgs$sampler_nuis, start_mu = NULL,
                                 start_var = NULL, type = type_nuis)
      startpoints_comb[[cc]] <- merge_group_level(
        sp$tmu, sp_nuis$tmu, sp$tvar, sp_nuis$tvar, pmwgs$nuisance, sp$subj_mu)
      emc_list[[cc]]$sampler_nuis$samples <- fill_samples(
        samples = pmwgs$sampler_nuis$samples, group_level = sp_nuis, j = 1,
        proposals = NULL, n_pars = pmwgs$n_pars, type = type_nuis)
      emc_list[[cc]]$sampler_nuis$samples$idx <- 1
    }
  }

  static <- list(samplers = emc_list, particles = particles,
                 type = emc_list[[1]]$type, r_cores = r_cores)
  pool <- .emc_pool_start(n_workers, static)
  on.exit(.emc_pool_stop(pool), add = TRUE)
  .emc_pool_set_iter(pool, list(startpoints = startpoints_comb))

  ord <- .emc_subject_order(emc_list[[1]]$data)
  tasks <- list()
  for (s in ord) for (cc in chain_idx) {
    tasks[[length(tasks) + 1L]] <- list(chain = cc, subj = s, seed = streams[[cc]][[s]])
  }
  results <- .emc_pool_apply(pool, tasks, .emc_start_task, n_workers,
                             list(startpoints = startpoints_comb))

  n_pars <- emc_list[[1]]$n_pars
  props <- lapply(seq_along(emc_list), function(i) matrix(0, n_pars + 1L, n_subjects))
  for (r in results) {
    props[[r$chain]][, r$subj] <- c(r$res$proposal, r$res$ll)
  }
  for (cc in chain_idx) {
    emc_list[[cc]]$samples <- fill_samples(
      samples = emc_list[[cc]]$samples, group_level = startpoints[[cc]],
      proposals = props[[cc]], j = 1, n_pars = n_pars, type = emc_list[[cc]]$type)
    emc_list[[cc]]$init <- TRUE
  }
  emc_list
}

# Flattened equivalent of run_stage().  Chains advance in lockstep: the Gibbs
# steps run in the master, then every chain's subjects go onto one queue.
run_stage_flat <- function(samplers, stage, iter, particles, tune,
                           n_workers = 1, r_cores = 1,
                           verbose = TRUE, verboseProgress = TRUE) {
  n_chains <- length(samplers)
  n_pars <- samplers[[1]]$n_pars
  n_subjects <- samplers[[1]]$n_subjects
  nuisance <- samplers[[1]]$nuisance
  type <- samplers[[1]]$type
  marginalise <- samplers[[1]]$marginalise

  # Components come from the (shared) data, so tuning settings are common.
  tune$components <- attr(samplers[[1]]$data, "components")
  tune$shared_ll_idx <- attr(samplers[[1]]$data, "shared_ll_idx")
  tune <- check_tune_settings(tune, n_pars, stage, particles)

  pm_settings <- vector("list", n_chains)
  eff_mu <- eff_var <- chains_mu <- chains_var <- vector("list", n_chains)
  chol_caches <- vector("list", n_chains)
  start_iter <- integer(n_chains)

  marginal_idx_blk <- .marginal_par_idx(samplers[[1]]$par_names, marginalise)
  if (length(marginal_idx_blk) != length(tune$components)) {
    marginal_idx_blk <- rep(FALSE, length(tune$components))
  }
  idx_list_blk <- .component_idx_list(tune$components, marginal_idx_blk)

  for (cc in seq_len(n_chains)) {
    s <- samplers[[cc]]
    pms <- attr(s$samples, "pm_settings")
    if (is.null(pms)) {
      pms <- lapply(seq_len(n_subjects), function(x)
        vector("list", length(unique(tune$components))))
    }
    pm_settings[[cc]] <- lapply(pms, FUN = check_sampling_settings, stage = stage,
                                n_pars = n_pars, particles)
    start_iter[cc] <- s$samples$idx
    samplers[[cc]] <- extend_sampler(s, iter, stage)

    empty <- vector("list", n_subjects)
    eff_mu[[cc]]     <- if (is.null(s$eff_mu))     empty else s$eff_mu
    eff_var[[cc]]    <- if (is.null(s$eff_var))    empty else s$eff_var
    chains_mu[[cc]]  <- if (is.null(s$chains_mu))  empty else s$chains_mu
    chains_var[[cc]] <- if (is.null(s$chains_var)) empty else s$chains_var
    chol_caches[[cc]] <- lapply(seq_len(n_subjects), function(ss) {
      build_subject_chol_cache(chains_var[[cc]][[ss]], eff_var[[cc]][[ss]], idx_list_blk)
    })
    if (any(nuisance)) {
      samplers[[cc]]$sampler_nuis$samples$idx <- samplers[[cc]]$samples$idx
    }
  }

  streams <- .emc_init_streams(n_chains, n_subjects)
  static <- list(data = samplers[[1]]$data, model = samplers[[1]]$model,
                 stage = stage, type = type, tune = tune,
                 marginalise = marginalise, r_cores = r_cores,
                 eff_mu = eff_mu, eff_var = eff_var,
                 chains_mu = chains_mu, chains_var = chains_var,
                 chol_caches = chol_caches)
  pool <- .emc_pool_start(n_workers, static)
  on.exit(.emc_pool_stop(pool), add = TRUE)

  ord <- .emc_subject_order(samplers[[1]]$data)
  if (verboseProgress) pb <- accept_progress_bar(min = 0, max = iter)

  for (i in seq_len(iter)) {
    if (verboseProgress) {
      update_progress_bar(pb, i, extra = mean(accept_rate(samplers[[1]])))
    }

    pars_comb_list <- vector("list", n_chains)
    pars_list <- vector("list", n_chains)
    group_chols <- vector("list", n_chains)
    active <- rep(TRUE, n_chains)

    # Gibbs steps, one chain at a time in the master.  A numerical failure here
    # means no subject of that chain has been updated yet, so the safe move is
    # to repeat its previous iteration -- exactly as in run_stage().
    for (cc in seq_len(n_chains)) {
      j <- start_iter[cc] + i
      pars_attempt <- tryCatch(
        gibbs_step(samplers[[cc]],
                   samplers[[cc]]$samples$alpha[!nuisance, , j - 1], type),
        error = identity)
      if (inherits(pars_attempt, c("error", "try-error"))) {
        samplers[[cc]]$samples <- reject_sample_iteration(samplers[[cc]]$samples, j)
        if (any(nuisance)) {
          samplers[[cc]]$sampler_nuis$samples <-
            reject_sample_iteration(samplers[[cc]]$sampler_nuis$samples, j)
        }
        active[cc] <- FALSE
        next
      }
      pars <- pars_comb <- pars_attempt
      if (any(nuisance)) {
        pars_nuis_attempt <- tryCatch(
          gibbs_step(samplers[[cc]]$sampler_nuis,
                     samplers[[cc]]$samples$alpha[nuisance, , j - 1],
                     samplers[[cc]]$sampler_nuis$type),
          error = identity)
        if (inherits(pars_nuis_attempt, c("error", "try-error"))) {
          samplers[[cc]]$samples <- reject_sample_iteration(samplers[[cc]]$samples, j)
          samplers[[cc]]$sampler_nuis$samples <-
            reject_sample_iteration(samplers[[cc]]$sampler_nuis$samples, j)
          active[cc] <- FALSE
          next
        }
        pars_nuis <- pars_nuis_attempt
        pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu, pars$tvar,
                                       pars_nuis$tvar, nuisance, pars$subj_mu)
        pars_comb$alpha <- samplers[[cc]]$samples$alpha[, , j - 1]
        samplers[[cc]]$sampler_nuis$samples <- fill_samples(
          samples = samplers[[cc]]$sampler_nuis$samples, group_level = pars_nuis,
          j = j, proposals = NULL, n_pars = n_pars,
          type = samplers[[cc]]$sampler_nuis$type)
        samplers[[cc]]$sampler_nuis$samples$idx <- j
      }
      pars_list[[cc]] <- pars
      pars_comb_list[[cc]] <- pars_comb
      group_var_it <- tryCatch(
        .apply_marginal_group_var(get_group_level(pars_comb, 1L, type)$var,
                                  marginal_idx_blk, marginalise),
        error = function(e) NULL)
      group_chols[[cc]] <- build_group_chol_cache(group_var_it, idx_list_blk)
    }

    iter_state <- list(pars = pars_comb_list, group_chol = group_chols)
    .emc_pool_set_iter(pool, iter_state)

    # One flat queue over every (chain, subject) still in play, longest first.
    tasks <- vector("list", sum(active) * n_subjects)
    k <- 0L
    for (ss in ord) for (cc in seq_len(n_chains)) {
      if (!active[cc]) next
      k <- k + 1L
      tasks[[k]] <- list(chain = cc, subj = ss,
                         pm_settings = pm_settings[[cc]][[ss]],
                         prev_ll = samplers[[cc]]$samples$subj_ll[ss, start_iter[cc] + i - 1],
                         seed = streams[[cc]][[ss]])
    }
    length(tasks) <- k
    results <- .emc_pool_apply(pool, tasks, .emc_particle_task, n_workers, iter_state)

    props <- lapply(seq_len(n_chains), function(x) matrix(0, n_pars + 1L, n_subjects))
    for (r in results) {
      props[[r$chain]][, r$subj] <- c(r$res$proposal, r$res$ll)
      pm_settings[[r$chain]][[r$subj]] <- r$res$pm_settings
      streams[[r$chain]][[r$subj]] <- r$seed
    }

    for (cc in seq_len(n_chains)) {
      if (!active[cc]) next
      j <- start_iter[cc] + i
      samplers[[cc]]$samples <- fill_samples(
        samples = samplers[[cc]]$samples, group_level = pars_list[[cc]],
        proposals = props[[cc]], j = j, n_pars = n_pars, type = type)
    }
  }

  for (cc in seq_len(n_chains)) {
    attr(samplers[[cc]]$samples, "pm_settings") <- pm_settings[[cc]]
  }
  if (verboseProgress) close(pb)
  samplers
}

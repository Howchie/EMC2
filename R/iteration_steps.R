# One sampler iteration as a sequence of steps over a shared state, for the
# configurations .emc_steps_ok() admits. Each step reads and writes the state
# environment `st` (the sampler, the group draw, the subjects' states, their
# proposal settings and random streams, the block's records) given the
# block's constants `cx`; run_stage() runs the same iteration inline for
# every other configuration. The steps are the units the native loop
# (src/iteration_loop.cpp) runs or replaces.

# Standard hierarchies without nuisance or pooled parameters, outside the
# profiled and progress-bar runs, with no likelihood pool of a single subject.
.emc_steps_ok <- function(pmwgs, nuisance, sampler_profile, verbose_progress, wpool,
                          llpool) {
  isTRUE(getOption("emc2.iteration_steps", TRUE)) &&
    identical(pmwgs$type, "standard") && !any(nuisance) && is.null(pmwgs$pooled) &&
    !sampler_profile && !verbose_progress && is.null(llpool) && !is.null(wpool) &&
    pmwgs$n_subjects > 1L
}

# The group draw given the subjects' last states, the ensemble pool members
# drawn from the same stream, and the factors of the group covariance.
.emc_step_group <- function(st, cx, i) {
  j <- cx$start_iter + i
  if (cx$persistent_workers) st$pool_iteration <- st$pool_iteration + 1L
  if (cx$alloc_on) {
    st$el_iter <- proc.time()[["elapsed"]]
    st$cores_it <- if (is.null(cx$core_ctl)) cx$n_cores else
      max(1L, cx$core_ctl$total %/% cx$core_ctl$n_slots)
  }
  pmwgs <- st$pmwgs
  assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
  alpha_prev <- matrix(pmwgs$samples$alpha[, , j - 1], nrow = dim(pmwgs$samples$alpha)[1],
                       dimnames = dimnames(pmwgs$samples$alpha)[1:2])
  pars <- tryCatch(gibbs_step(pmwgs, alpha_prev, pmwgs$type), error = identity)
  if (inherits(pars, c("error", "try-error"))) {
    .emc_reject_record(pars, "gibbs")
    .emc_gibbs_abort(pars, stage = cx$stage, iteration = j)
  }
  ens_q <- pmwgs$ensemble_q
  ens_pools <- .emc_ensemble_draw_pools(pmwgs, cx$stage,
                                        if (is.null(ens_q)) cx$chains_mu else ens_q$mu,
                                        if (is.null(ens_q)) cx$chains_var else ens_q$var,
                                        factors = st$ens_factors)
  if (!is.null(ens_pools) && is.null(st$ens_factors))
    st$ens_factors <- list(L = ens_pools$L, logdet = ens_pools$logdet)
  pmwgs$rng$gibbs <- get(".Random.seed", envir = globalenv())
  group_var_it <- tryCatch(get_group_level(pars, 1L, pmwgs$type)$var,
                           error = function(e) NULL)
  if (is.null(st$group_chol_it) || !identical(st$group_chol_it$ref, group_var_it))
    st$group_chol_it <- build_group_chol_cache(group_var_it, cx$idx_list_blk,
                                               cx$marginal_idx_blk,
                                               include_full = cx$include_full)
  st$pmwgs <- pmwgs
  st$pars <- pars
  st$ens_pools <- ens_pools
  invisible(NULL)
}

# The subjects' particle steps, with the allocation's skips drawn from each
# subject's own stream, and the ensemble candidates' likelihoods.
.emc_step_particle <- function(st, cx, i) {
  j <- cx$start_iter + i
  pmwgs <- st$pmwgs
  wpool_streams <- st$streams
  hmc_plan <- st$hmc_plan
  part_skip <- NULL
  hmc_pick <- NULL
  q_part <- if (cx$alloc_on) pmwgs$alloc$particle$rate
  if (length(q_part) == pmwgs$n_subjects && any(q_part < 1)) {
    part_skip <- logical(pmwgs$n_subjects)
    r_now <- .emc_hmc_plan_rates(hmc_plan, pmwgs$n_subjects)
    r_alloc <- pmwgs$alloc$last$rate
    if (length(r_alloc) != pmwgs$n_subjects) r_alloc <- r_now
    .emc_with_preserved_rng(for (s in which(q_part < 1)) {
      if (r_alloc[s] > 0 && !(r_now[s] > 0)) next
      assign(".Random.seed", wpool_streams[[s]], envir = globalenv())
      u <- stats::runif(1)
      wpool_streams[[s]] <- get(".Random.seed", envir = globalenv())
      part_skip[s] <- u >= q_part[s]
      if (r_now[s] > 0 && r_now[s] < 1) {
        if (is.null(hmc_pick)) hmc_pick <- rep(NA, pmwgs$n_subjects)
        hmc_pick[s] <- u >= 1 - r_now[s]
      }
    })
  }
  part_started <- proc.time()[["elapsed"]]
  part_cpu_share <- 1
  wpool <- st$wpool
  wpool_part <- st$wpool_part
  wpool_cost <- st$wpool_cost
  ens_pools <- st$ens_pools
  wpool <- .emc_wpool_recycle(
    wpool, if (cx$persistent_workers) st$pool_iteration else i,
    cx$recycle_every, cx$wpool_ctx)
  chain_cores <- .emc_cores_now(cx$core_ctl, cx$n_cores)
  target_workers <- min(
    pmwgs$n_subjects,
    max(1L, chain_cores %/% max(1L, cx$pool_budget$likelihood))
  )
  n_before <- wpool$n
  wpool <- .emc_wpool_grow(wpool, target_workers, cx$wpool_ctx)
  if (cx$persistent_workers) wpool$iterations <- st$pool_iteration
  if (wpool$n != n_before || !is.null(part_skip)) {
    wpool_part <- .emc_lpt_partition(
      if (is.null(part_skip)) wpool_cost else wpool_cost * !part_skip, wpool$n)
  }
  it <- .emc_wpool_iter(wpool, cx$wpool_ctx, wpool_part, st$pars,
                        st$group_chol_it, st$pm_settings,
                        pmwgs$samples$subj_ll[, j - 1], wpool_streams,
                        extra = ens_pools$rows,
                        pooled = NULL, skip = part_skip)
  if (!is.null(ens_pools)) ens_pools$ll <- it$extra_ll
  proposals <- it$props
  wpool$alive <- it$alive
  if (!is.null(it$rejects))
    st$stage_rejects$particle <- st$stage_rejects$particle + it$rejects
  if (!is.null(it$ensemble_rejects))
    .emc_reject_add(it$ensemble_rejects, "ensemble")
  ran <- if (is.null(part_skip)) rep(TRUE, length(it$times)) else !part_skip
  if (any(it$times[ran] > 0)) {
    fresh <- pmax(it$times, min(it$times[ran & it$times > 0]))
    wpool_cost[ran] <- fresh[ran]
  }
  wpool_part <- .emc_lpt_partition(wpool_cost, wpool$n)
  if (is.numeric(it$worker_cpu) && it$worker_cpu > 0)
    part_cpu_share <- min(1, sum(it$times) / it$worker_cpu)
  if (cx$alloc_on) {
    pmwgs <- .emc_alloc_record_particle(
      pmwgs, matrix(pmwgs$samples$alpha[, , j - 1L], nrow = pmwgs$n_pars),
      proposals[seq_len(pmwgs$n_pars), , drop = FALSE],
      if (is.null(part_skip)) rep(TRUE, pmwgs$n_subjects) else !part_skip,
      part_cpu_share * .emc_reserved_cost(proc.time()[["elapsed"]] - part_started,
                                          st$cores_it),
      it$times)
  }
  st$pmwgs <- pmwgs
  st$proposals <- proposals
  st$pm_settings <- it$pm_settings
  st$streams <- it$seeds
  st$wpool <- wpool
  st$wpool_part <- wpool_part
  st$wpool_cost <- wpool_cost
  st$ens_pools <- ens_pools
  st$hmc_pick <- hmc_pick
  invisible(NULL)
}

# The ensemble group update, after the subject step.
.emc_step_ensemble <- function(st, cx, i) {
  ens <- .emc_ensemble_move(st$pmwgs, st$pars, st$proposals, cx$stage, st$wpool,
                            cx$wpool_ctx, st$wpool_part, st$ens_pools,
                            ll_cache_validated = TRUE)
  pmwgs <- ens$sampler
  if (cx$alloc_on) pmwgs <- .emc_alloc_record_group(pmwgs, ens$pars)
  st$pmwgs <- pmwgs
  st$pars <- ens$pars
  st$proposals <- ens$proposals
  st$wpool <- ens$wpool
  invisible(NULL)
}

# The shear move of the same group state.
.emc_step_shear <- function(st, cx, i) {
  st$shear_cpu <- st$shear_el <- st$shear_round_cpu <- 0
  st$shear_rounds <- 0L
  if (!.emc_shear_active(st$pmwgs)) return(invisible(NULL))
  shear_el0 <- proc.time()[["elapsed"]]
  sh <- .emc_shear_move(st$pmwgs, st$pars, st$proposals, cx$stage, st$wpool,
                        st$wpool_part)
  st$pmwgs <- sh$sampler
  st$pars <- sh$pars
  st$proposals <- sh$proposals
  st$wpool <- sh$wpool
  st$shear_rounds <- if (is.null(sh$rounds)) 0L else sh$rounds
  shear_round_el <- if (is.null(sh$round_el)) 0 else sh$round_el
  st$shear_el <- proc.time()[["elapsed"]] - shear_el0
  if (cx$alloc_on) {
    st$shear_cpu <- .emc_reserved_cost(st$shear_el, st$cores_it)
    st$shear_round_cpu <- .emc_reserved_cost(shear_round_el, st$cores_it)
  }
  invisible(NULL)
}

# The collapsed move of the whole group state.
.emc_step_collapsed <- function(st, cx, i) {
  st$coll_cpu <- st$coll_el <- st$coll_round_cpu <- 0
  st$coll_moves <- 0L
  if (!.emc_collapsed_active(st$pmwgs)) return(invisible(NULL))
  coll_el0 <- proc.time()[["elapsed"]]
  moves_before <- st$pmwgs$collapsed$plan$block$n
  cm <- .emc_collapsed_move(st$pmwgs, st$pars, st$proposals, cx$stage, st$wpool,
                            st$wpool_part)
  st$pmwgs <- cm$sampler
  st$pars <- cm$pars
  st$proposals <- cm$proposals
  st$wpool <- cm$wpool
  moves_after <- st$pmwgs$collapsed$plan$block$n
  if (length(moves_before) && length(moves_after))
    st$coll_moves <- as.integer(moves_after - moves_before)
  st$coll_el <- proc.time()[["elapsed"]] - coll_el0
  if (cx$alloc_on) {
    st$coll_cpu <- .emc_reserved_cost(st$coll_el, st$cores_it)
    st$coll_round_cpu <- .emc_reserved_cost(if (is.null(cm$round_el)) 0 else cm$round_el,
                                            st$cores_it)
  }
  invisible(NULL)
}

# Subject HMC, the last update of the iteration, from the final group state.
.emc_step_hmc <- function(st, cx, i) {
  for (kind in c("ensemble", "shear", "collapsed")) {
    delta <- .emc_reject_delta(st$rejects[[kind]], kind)
    st$rejects[[kind]] <- .emc_reject_counts(kind)
    st$stage_rejects[[kind]] <- st$stage_rejects[[kind]] + delta
  }
  pmwgs <- st$pmwgs
  proposals <- st$proposals
  hmc_tasks <- NULL
  hmc_el <- 0
  if (cx$alloc_on)
    pmwgs <- .emc_alloc_record_subjects(pmwgs,
                                        proposals[seq_len(pmwgs$n_pars), , drop = FALSE],
                                        st$pars$tmu)
  if (!is.null(st$hmc_plan)) {
    hmc_started <- .emc_cpu_self()
    hmc_el0 <- proc.time()[["elapsed"]]
    hr <- .emc_hmc_round(st$hmc_plan, proposals, st$pars, pmwgs$type, st$streams,
                         cx$data, pmwgs$model, wpool = st$wpool,
                         n_cores = cx$pool_budget$subject,
                         r_cores = cx$pool_budget$likelihood,
                         par_names = cx$par_names, pick = st$hmc_pick,
                         ll_cache_validated = TRUE)
    proposals <- hr$proposals
    st$streams <- hr$streams
    hmc_tasks <- hr$tasks
    st$hmc_plan <- hr$plan
    if (!is.null(hr$pool)) st$wpool$alive <- hr$pool$alive
    hmc_master <- .emc_cpu_self() - hmc_started
    hmc_el <- proc.time()[["elapsed"]] - hmc_el0
    attempted <- !vapply(hmc_tasks, is.null, logical(1))
    if (any(attempted)) {
      for (s in which(attempted)) {
        hmc_tasks[[s]]$work <- hmc_tasks[[s]]$cpu
        hmc_tasks[[s]]$cpu <- hmc_tasks[[s]]$cpu + hmc_master / sum(attempted)
      }
      # The round's wall time on the chain's cores, shared by CPU.
      if (cx$alloc_on) {
        cpu_s <- vapply(hmc_tasks[attempted], function(x) x$cpu, numeric(1))
        share <- if (sum(cpu_s) > 0) cpu_s / sum(cpu_s) else
          rep(1 / sum(attempted), sum(attempted))
        for (k in seq_along(share))
          hmc_tasks[[which(attempted)[k]]]$cost <-
            .emc_reserved_cost(hmc_el, st$cores_it) * share[k]
      }
    }
    pmwgs$subject_hmc <- .emc_hmc_record_tasks(pmwgs$subject_hmc, hmc_tasks)
  }
  st$pmwgs <- pmwgs
  st$proposals <- proposals
  st$hmc_tasks <- hmc_tasks
  st$hmc_el <- hmc_el
  invisible(NULL)
}

# The iteration's draw into the sample store, and its costs into the
# allocation's records.
.emc_step_store <- function(st, cx, i) {
  j <- cx$start_iter + i
  pmwgs <- st$pmwgs
  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = st$pars,
                                proposals = st$proposals, j = j, n_pars = pmwgs$n_pars,
                                type = pmwgs$type)
  if (cx$alloc_on) {
    rest_cpu <- .emc_reserved_cost(proc.time()[["elapsed"]] - st$el_iter - st$shear_el -
                                     st$coll_el - st$hmc_el, st$cores_it)
    pmwgs <- .emc_alloc_record_cost(pmwgs, rest_cpu, st$hmc_tasks, moves = list(
      shear = c(cpu = st$shear_cpu, own = st$shear_round_cpu, n = st$shear_rounds),
      collapsed = c(cpu = st$coll_cpu, own = st$coll_round_cpu, n = st$coll_moves),
      hmc = .emc_reserved_cost(st$hmc_el, st$cores_it)))
  }
  delta <- .emc_reject_delta(st$rejects$gibbs, "gibbs")
  st$rejects$gibbs <- .emc_reject_counts("gibbs")
  st$stage_rejects$gibbs <- st$stage_rejects$gibbs + delta
  st$pmwgs <- pmwgs
  invisible(NULL)
}

.EMC_ITERATION_STEPS <- c("group", "particle", "ensemble", "shear", "collapsed", "hmc",
                          "store")

# Iterations `from` to cx$iter of the block as the step sequence.
.emc_steps_run <- function(st, cx, from = 1L) {
  for (i in seq.int(from, cx$iter, length.out = max(0L, cx$iter - from + 1L))) {
    .emc_step_group(st, cx, i)
    .emc_step_particle(st, cx, i)
    .emc_step_ensemble(st, cx, i)
    .emc_step_shear(st, cx, i)
    .emc_step_collapsed(st, cx, i)
    .emc_step_hmc(st, cx, i)
    .emc_step_store(st, cx, i)
  }
  invisible(NULL)
}

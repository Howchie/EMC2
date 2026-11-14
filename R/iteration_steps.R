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

# What the native steps take ready-made for the block: the group draw's
# blocks (NULL with group designs, which the R update handles), the ensemble
# pools' inputs ("off" when no pools are drawn), the native particle kernel's
# tuning constants (NULL when the kernel does not cover the block), and
# whether the group covariance has one factor over all the parameters.
.emc_steps_constants <- function(pmwgs, stage, chains_mu, chains_var, idx_list_blk,
                                 marginal_idx_blk, include_full, ctx,
                                 persistent_workers) {
  p <- pmwgs$n_pars
  n <- pmwgs$n_subjects
  out <- list()
  if (is.null(pmwgs$gd) && !length(pmwgs$group_designs)) {
    is_blocked <- pmwgs$is_blocked[!pmwgs$nuisance]
    par_group <- pmwgs$par_group[!pmwgs$nuisance]
    if (is.null(is_blocked) || length(is_blocked) == 0L) is_blocked <- rep(TRUE, p)
    if (is.null(par_group)) par_group <- rep(1L, p)
    out$gibbs <- list(blocked = as.logical(is_blocked),
                      group = as.integer(factor(par_group, levels = unique(par_group))))
  }
  ens_q <- pmwgs$ensemble_q
  mu <- if (is.null(ens_q)) chains_mu else ens_q$mu
  var <- if (is.null(ens_q)) chains_var else ens_q$var
  cfg <- .emc_ensemble_options()
  out$pools <- if (stage == "preburn" || !.emc_ensemble_supported(pmwgs) ||
                   length(mu) != n || length(var) != n ||
                   any(vapply(mu, is.null, logical(1))) ||
                   any(vapply(var, is.null, logical(1)))) "off" else
    list(mu = lapply(mu, as.numeric), var = lapply(var, as.matrix), size = cfg$pool - 1L,
         rounds = cfg$rounds,
         tail_weight = cfg$tail_weight, tail_df = cfg$tail_df,
         tail_excess = cfg$tail_excess, par_names = pmwgs$par_names)
  # The particle steps run on the native kernel (.emc_native_particle_batch)
  # with one block of parameters, in this process.
  tune <- ctx$tune
  tune_in <- list(n0 = as.numeric(tune$n0),
                  step_target = if (is.null(tune$step_target)) step_target(1) else
                    tune$step_target,
                  mix_adapt = as.numeric(tune$mix_adapt),
                  target_ESS = as.numeric(tune$target_ESS),
                  ESS_scale = as.numeric(tune$ESS_scale),
                  max_particles = as.numeric(tune$max_particles))
  if (isTRUE(getOption("emc2.native_particle_block", TRUE)) &&
      is.null(ctx$marginalise) && length(unique(tune$components)) == 1L &&
      length(unique(tune$shared_ll_idx)) == 1L && !isTRUE(ctx$r_cores > 1L) &&
      !(!is.null(.emc_pool_state$ll_pool) && isTRUE(.emc_pool_state$ll_pool$alive)) &&
      !nzchar(Sys.getenv("EMC2_LL_TRACE")) && isTRUE(ctx$ll_cache_validated) &&
      !persistent_workers && all(lengths(tune_in) == 1L))
    out$particle <- list(tune_in = tune_in)
  out$shear_active <- .emc_shear_active(pmwgs)
  out$shear_step_range <- .EMC_SHEAR_LOG_STEP_RANGE
  out$chol_simple <- !isTRUE(include_full) && length(idx_list_blk) == 1L &&
    is.logical(idx_list_blk[[1L]]) && length(idx_list_blk[[1L]]) == p &&
    all(idx_list_blk[[1L]]) && !any(marginal_idx_blk)
  out
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
  st$pmwgs <- pmwgs
  .emc_steps_rejects(st)
  invisible(NULL)
}

# The iteration's rejected updates of each group move, into the stage's
# counts.
.emc_steps_rejects <- function(st) {
  for (kind in c("ensemble", "shear", "collapsed", "gibbs")) {
    now <- .emc_reject_counts(kind)
    if (identical(now, st$rejects[[kind]])) next
    st$stage_rejects[[kind]] <- st$stage_rejects[[kind]] + (now - st$rejects[[kind]])
    st$rejects[[kind]] <- now
  }
  invisible(NULL)
}

.EMC_ITERATION_STEPS <- c("group", "particle", "ensemble", "shear", "collapsed", "hmc",
                          "store")

# Iterations `from` to cx$iter of the block: the native loop
# (src/iteration_loop.cpp), which runs each step itself where it can and
# calls the step above where it cannot.
.emc_steps_run <- function(st, cx, from = 1L) {
  steps <- list(group = .emc_step_group, particle = .emc_step_particle,
                ensemble = .emc_step_ensemble, shear = .emc_step_shear,
                collapsed = .emc_step_collapsed, hmc = .emc_step_hmc,
                store = .emc_step_store)
  native <- getOption("emc2.native_steps", TRUE)
  native <- if (is.character(native)) .EMC_ITERATION_STEPS %in% native else
    rep(isTRUE(native), length(steps))
  emc_iteration_loop_cpp(st, cx, as.integer(from), steps, native, .emc_steps_rejects)
  invisible(NULL)
}

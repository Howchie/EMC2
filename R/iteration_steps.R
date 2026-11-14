# One sampler iteration as a sequence of steps over a shared state. Each step
# reads and writes the state environment `st` (the sampler, the group draw,
# the subjects' states, their proposal settings and random streams, the
# worker pools, the block's records) given the block's constants `cx`. The
# loop over a block's iterations (src/iteration_loop.cpp) runs each step in
# its native form where that form accepts the state, else through the
# functions here.

# The steps that have a native form.
.EMC_NATIVE_STEPS <- c(group = "group", particle = "particle", ensemble = "ensemble",
                       shear = "shear", collapsed = "collapsed", hmc = "hmc", store = "store")

# Whether the native forms of the steps cover the block: standard hierarchies
# and fits without a group level, without nuisance or pooled parameters,
# outside profiled runs. A step whose native form does not take the state
# (a likelihood pool of a single subject, no subject pool) runs in R.
.emc_steps_native_ok <- function(pmwgs, nuisance, sampler_profile, wpool, llpool) {
  (identical(pmwgs$type, "standard") || identical(pmwgs$type, "single")) &&
    !any(nuisance) && is.null(pmwgs$pooled) && !sampler_profile
}

# The state's own value of a block constant that a step may replace (the data
# and the pool context under pooled parameters), else the constant.
.emc_steps_now <- function(st, cx, name) {
  v <- st[[name]]
  if (is.null(v)) cx[[name]] else v
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
  # Without a group level the "group draw" is the prior itself.
  if (identical(pmwgs$type, "single"))
    out$single <- list(tmu = pmwgs$prior$theta_mu_mean, tvar = pmwgs$prior$theta_mu_var,
                       rows = pmwgs$par_names[!pmwgs$nuisance])
  else if (is.null(pmwgs$gd) && !length(pmwgs$group_designs)) {
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
  out$subject_cost <- .emc_subject_cost(ctx$data)
  out$hmc_warmup <- c(accept = .EMC_HMC_WARMUP_ACCEPT, moves = .EMC_HMC_WARMUP_MOVES,
                      eps_low = .EMC_HMC_EPS_RANGE[1], eps_high = .EMC_HMC_EPS_RANGE[2])
  out$shear_active <- .emc_shear_active(pmwgs)
  out$shear_step_range <- .EMC_SHEAR_LOG_STEP_RANGE
  out$chol_simple <- !isTRUE(include_full) && length(idx_list_blk) == 1L &&
    is.logical(idx_list_blk[[1L]]) && length(idx_list_blk[[1L]]) == p &&
    all(idx_list_blk[[1L]]) && !any(marginal_idx_blk)
  out
}

# The acceptance rate so far, on the block's progress bar.
.emc_step_progress <- function(st, cx, i) {
  update_progress_bar(cx$pb, i, extra = mean(accept_rate(st$pmwgs)))
  invisible(NULL)
}

# The group draw given the subjects' last states (with the nuisance
# parameters' own draw), the ensemble pool members drawn from the same
# stream, and the factors of the group covariance.
.emc_step_group <- function(st, cx, i) {
  j <- cx$start_iter + i
  prof <- cx$sampler_profile
  if (cx$persistent_workers) st$pool_iteration <- st$pool_iteration + 1L
  if (prof) st$prof <- list(started = proc.time()[["elapsed"]])
  if (cx$alloc_on) {
    # The allocation charges wall time on the chain's own share of the
    # cores: those it takes over when other chains finish their block would
    # otherwise idle, so they cost the fit nothing.
    st$el_iter <- proc.time()[["elapsed"]]
    st$cores_it <- if (is.null(cx$core_ctl)) cx$n_cores else
      max(1L, cx$core_ctl$total %/% cx$core_ctl$n_slots)
  }
  pmwgs <- st$pmwgs
  nuisance <- cx$nuisance
  assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
  gibbs_started <- if (prof) proc.time()[["elapsed"]]
  # Parameters x subjects, kept a matrix with one parameter or one subject.
  alpha_prev <- matrix(pmwgs$samples$alpha[, , j - 1], nrow = dim(pmwgs$samples$alpha)[1],
                       dimnames = dimnames(pmwgs$samples$alpha)[1:2])
  pars <- tryCatch(gibbs_step(pmwgs, if (cx$any_nuisance)
    alpha_prev[!nuisance, , drop = FALSE] else alpha_prev, pmwgs$type), error = identity)
  if (inherits(pars, c("error", "try-error"))) {
    .emc_reject_record(pars, "gibbs")
    .emc_gibbs_abort(pars, stage = cx$stage, iteration = j)
  }
  pars_comb <- pars
  if (cx$any_nuisance) {
    pars_nuis <- tryCatch(
      gibbs_step(pmwgs$sampler_nuis, alpha_prev[nuisance, , drop = FALSE],
                 pmwgs$sampler_nuis$type),
      error = identity)
    if (inherits(pars_nuis, c("error", "try-error"))) {
      .emc_reject_record(pars_nuis, "gibbs")
      .emc_gibbs_abort(pars_nuis, stage = cx$stage, iteration = j, nuisance = TRUE)
    }
    # pars_comb only feeds the particle step, so it carries the subject prior
    # covariance (SEM: subj_var, which excludes covariate spread).
    pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu,
                                   if (is.null(pars$subj_var)) pars$tvar else pars$subj_var,
                                   pars_nuis$tvar, nuisance, pars$subj_mu)
    pars_comb$alpha <- pmwgs$samples$alpha[, , j - 1]
    pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                               group_level = pars_nuis, j = j,
                                               proposals = NULL, n_pars = cx$n_pars,
                                               type = pmwgs$sampler_nuis$type)
    pmwgs$sampler_nuis$samples$idx <- j
  }
  # Ensemble pool members come from the group-step stream and depend only on
  # the frozen proposals, so their likelihoods ride along with the subject
  # step; the proposals' factors from the first draw serve the rest.
  ens_q <- pmwgs$ensemble_q
  ens_pools <- .emc_ensemble_draw_pools(pmwgs, cx$stage,
                                        if (is.null(ens_q)) cx$chains_mu else ens_q$mu,
                                        if (is.null(ens_q)) cx$chains_var else ens_q$var,
                                        factors = st$ens_factors)
  if (!is.null(ens_pools) && is.null(st$ens_factors))
    st$ens_factors <- list(L = ens_pools$L, logdet = ens_pools$logdet)
  pmwgs$rng$gibbs <- get(".Random.seed", envir = globalenv())
  if (prof) st$prof$gibbs <- proc.time()[["elapsed"]] - gibbs_started
  cache_started <- if (prof) proc.time()[["elapsed"]]
  group_var_it <- tryCatch(get_group_level(pars_comb, 1L, pmwgs$type)$var,
                           error = function(e) NULL)
  # The factors depend only on the group variance, which is fixed when there
  # is no group level.
  if (is.null(st$group_chol_it) || !identical(st$group_chol_it$ref, group_var_it))
    st$group_chol_it <- build_group_chol_cache(group_var_it, cx$idx_list_blk,
                                               cx$marginal_idx_blk,
                                               include_full = cx$include_full)
  if (prof) st$prof$group_cache <- proc.time()[["elapsed"]] - cache_started
  st$pmwgs <- pmwgs
  st$pars <- pars
  st$pars_comb <- if (cx$any_nuisance) pars_comb
  st$ens_pools <- ens_pools
  invisible(NULL)
}

# The subjects' particle steps, with the allocation's skips drawn from each
# subject's own stream, and the ensemble candidates' likelihoods.
.emc_step_particle <- function(st, cx, i) {
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
  .emc_step_particle_run(st, cx, i, part_skip, hmc_pick, wpool_streams)
}

# The particle step given the allocation's skips (NULL when no subject may
# skip), the joint selection with subject HMC and the subjects' streams after
# those draws: the steps on the worker pool, or in this process without one,
# and their records. `maintained`: the pool has had this iteration's recycling
# and growth.
.emc_step_particle_run <- function(st, cx, i, part_skip, hmc_pick, wpool_streams,
                                   maintained = FALSE) {
  j <- cx$start_iter + i
  prof <- cx$sampler_profile
  pmwgs <- st$pmwgs
  part_started <- proc.time()[["elapsed"]]
  part_cpu_share <- 1
  wpool <- st$wpool
  wpool_part <- st$wpool_part
  wpool_cost <- st$wpool_cost
  wpool_ctx <- .emc_steps_now(st, cx, "wpool_ctx")
  data <- .emc_steps_now(st, cx, "data")
  ens_pools <- st$ens_pools
  pm_settings <- st$pm_settings
  pars_comb <- if (is.null(st$pars_comb)) st$pars else st$pars_comb
  recycle_elapsed <- grow_elapsed <- 0
  pool_profile <- NULL
  it <- NULL
  if (!is.null(st$llpool)) {
    llpool <- .emc_wpool_recycle(
      st$llpool, if (cx$persistent_workers) st$pool_iteration else i,
      cx$recycle_every, wpool_ctx)
    llpool <- .emc_wpool_grow(llpool, .emc_cores_now(cx$core_ctl, cx$n_cores), wpool_ctx)
    if (cx$persistent_workers) llpool$iterations <- st$pool_iteration
    .emc_ll_pool_set(llpool, 1L, reset = FALSE)
    st$llpool <- llpool
  }
  if (!is.null(wpool)) {
    n_before <- wpool$n
    if (!maintained) {
      recycle_started <- if (prof) proc.time()[["elapsed"]]
      wpool <- .emc_wpool_recycle(
        wpool, if (cx$persistent_workers) st$pool_iteration else i,
        cx$recycle_every, wpool_ctx)
      if (prof) recycle_elapsed <- proc.time()[["elapsed"]] - recycle_started
      chain_cores <- .emc_cores_now(cx$core_ctl, cx$n_cores)
      target_workers <- min(
        pmwgs$n_subjects,
        max(1L, chain_cores %/% max(1L, cx$pool_budget$likelihood))
      )
      n_before <- wpool$n
      grow_started <- if (prof) proc.time()[["elapsed"]]
      wpool <- .emc_wpool_grow(wpool, target_workers, wpool_ctx)
      if (cx$persistent_workers) wpool$iterations <- st$pool_iteration
      if (prof) grow_elapsed <- proc.time()[["elapsed"]] - grow_started
    }
    if (wpool$n != n_before || !is.null(part_skip)) {
      wpool_part <- .emc_lpt_partition(
        if (is.null(part_skip)) wpool_cost else wpool_cost * !part_skip, wpool$n)
    }
    it <- .emc_wpool_iter(wpool, wpool_ctx, wpool_part, pars_comb,
                          st$group_chol_it, pm_settings,
                          pmwgs$samples$subj_ll[, j - 1], wpool_streams,
                          extra = ens_pools$rows,
                          pooled = pmwgs$pooled$value, skip = part_skip)
    if (!is.null(ens_pools)) ens_pools$ll <- it$extra_ll
    proposals <- it$props
    pm_settings <- it$pm_settings
    wpool_streams <- it$seeds
    wpool$alive <- it$alive
    pool_profile <- it$profile
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
    # The pass also computes the ensemble candidates: the particle steps'
    # share of it is their share of the workers' CPU.
    if (is.numeric(it$worker_cpu) && it$worker_cpu > 0)
      part_cpu_share <- min(1, sum(it$times) / it$worker_cpu)
  } else {
    core_budget <- .particle_core_budget(pmwgs$n_subjects, n_cores = cx$n_cores,
                                         r_cores = cx$r_cores,
                                         total_cores = cx$n_cores)
    st$core_budget <- core_budget
    tune <- cx$tune
    group_chol_it <- st$group_chol_it
    # Keep each result whole, including its subject tag and post-draw RNG state,
    # so a missing worker cannot shift or recycle another subject's result.
    raw_proposals <- auto_mclapply(seq_len(pmwgs$n_subjects), function(s) {
      assign(".Random.seed", wpool_streams[[s]], envir = globalenv())
      out <- if (isTRUE(part_skip[s])) {
        r <- reject_particle(pmwgs$samples$alpha[, s, j - 1L],
                             pmwgs$samples$subj_ll[s, j - 1L], pm_settings[[s]])
        r$subject <- s
        r
      } else safe_new_particle(
        s = s, data = data[[s]], pm_settings = pm_settings[[s]],
        eff_mu = cx$eff_mu[[s]], eff_var = cx$eff_var[[s]],
        chains_mu = cx$chains_mu[[s]], chains_var = cx$chains_var[[s]],
        prev_ll = pmwgs$samples$subj_ll[s, j - 1L], parameters = pars_comb,
        model = pmwgs$model, stage = cx$stage, type = pmwgs$type, tune = tune,
        marginalise = pmwgs$marginalise,
        r_cores = core_budget$likelihood, chol_cache = cx$chol_caches[[s]],
        group_chol = group_chol_it, ll_cache_validated = TRUE)
      out$rng_seed <- get(".Random.seed", envir = globalenv())
      if (!is.null(ens_pools))
        out$extra_ll <- .emc_ensemble_ll_candidate(ens_pools$rows[[s]], data[[s]],
                                                   pmwgs$model,
                                                   ll_cache_validated = TRUE)
      out
    }, mc.cores = core_budget$subject, mc.set.seed = FALSE)
    if (!is.null(ens_pools)) {
      ens_pools$ll <- vector("list", pmwgs$n_subjects)
      for (res in raw_proposals) if (is.list(res) && is.numeric(res$subject) &&
                                     length(res$subject) == 1L && !is.null(res$extra_ll))
        ens_pools$ll[[as.integer(res$subject)]] <- res$extra_ll
    }
    next_streams <- wpool_streams
    for (res in raw_proposals) {
      if (!is.list(res)) next
      s_res <- res$subject
      if (!is.numeric(s_res) || length(s_res) != 1L || is.na(s_res) ||
          s_res < 1L || s_res > pmwgs$n_subjects || is.null(res$rng_seed)) {
        next
      }
      next_streams[[as.integer(s_res)]] <- res$rng_seed
    }
    wpool_streams <- next_streams
    assembled <- .emc_assemble_proposals(
      raw_proposals, n_pars = pmwgs$n_pars, n_subjects = pmwgs$n_subjects,
      prev_alpha = matrix(pmwgs$samples$alpha[, , j - 1L], nrow = pmwgs$n_pars),
      prev_ll = pmwgs$samples$subj_ll[, j - 1L], prev_pm = pm_settings,
      stage = cx$stage, iteration = j)
    pm_settings <- assembled$pm_settings
    proposals <- assembled$proposals
  }
  if (cx$alloc_on) {
    pmwgs <- .emc_alloc_record_particle(
      pmwgs, matrix(pmwgs$samples$alpha[, , j - 1L], nrow = pmwgs$n_pars),
      proposals[seq_len(pmwgs$n_pars), , drop = FALSE],
      if (is.null(part_skip)) rep(TRUE, pmwgs$n_subjects) else !part_skip,
      part_cpu_share * .emc_reserved_cost(proc.time()[["elapsed"]] - part_started,
                                          st$cores_it),
      if (is.null(wpool)) NULL else it$times)
  }
  if (prof) {
    st$prof$recycle <- recycle_elapsed
    st$prof$grow <- grow_elapsed
    st$prof$pool <- pool_profile
  }
  st$pmwgs <- pmwgs
  st$proposals <- proposals
  st$pm_settings <- pm_settings
  st$streams <- wpool_streams
  st$wpool <- wpool
  st$wpool_part <- wpool_part
  st$wpool_cost <- wpool_cost
  st$ens_pools <- ens_pools
  st$hmc_pick <- hmc_pick
  invisible(NULL)
}

# The ensemble group update, after the subject step, so the group Gibbs
# rounds and the final reselection close the sweep.
.emc_step_ensemble <- function(st, cx, i) {
  started <- if (cx$sampler_profile) proc.time()[["elapsed"]]
  ens <- .emc_ensemble_move(st$pmwgs, st$pars, st$proposals, cx$stage, st$wpool,
                            .emc_steps_now(st, cx, "wpool_ctx"), st$wpool_part,
                            st$ens_pools, ll_cache_validated = TRUE)
  if (cx$sampler_profile) st$prof$ensemble <- proc.time()[["elapsed"]] - started
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
  prof <- cx$sampler_profile
  if (prof) st$prof$shear_raw <- st$prof$shear_build <- st$prof$shear_ll <- 0
  if (!.emc_shear_active(st$pmwgs)) return(invisible(NULL))
  shear_started <- if (prof) .emc_cpu_self()
  sll_shear <- if (prof) .emc_subject_ll_cpu()
  shear_el0 <- proc.time()[["elapsed"]]
  sh <- .emc_shear_move(st$pmwgs, st$pars, st$proposals, cx$stage, st$wpool,
                        st$wpool_part)
  st$pmwgs <- sh$sampler
  st$pars <- sh$pars
  st$proposals <- sh$proposals
  st$wpool <- sh$wpool
  st$shear_rounds <- if (is.null(sh$rounds)) 0L else sh$rounds
  shear_round_el <- if (is.null(sh$round_el)) 0 else sh$round_el
  if (prof) {
    st$prof$shear_build <- if (is.null(sh$build_el)) 0 else sh$build_el
    st$prof$shear_ll <- if (is.null(sh$ll_el)) 0 else sh$ll_el
    st$prof$shear_raw <- .emc_cpu_self() - shear_started + .emc_subject_ll_cpu() - sll_shear
  }
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

# Pooled coefficients given the final subject states (R/pooled.R); the data
# and the pool context carry their new values from here.
.emc_step_pooled <- function(st, cx, i) {
  pmwgs <- st$pmwgs
  wpool_ctx <- .emc_steps_now(st, cx, "wpool_ctx")
  pmv <- .emc_pooled_move(pmwgs, st$proposals, cx$stage, st$pars, st$group_chol_it,
                          st$pm_settings, st$streams, st$wpool, wpool_ctx,
                          st$wpool_part, n_cores = cx$n_cores,
                          ll_cache_validated = TRUE)
  pmwgs <- pmv$sampler
  st$proposals <- pmv$proposals
  st$pars <- pmv$pars
  if (!is.null(st$wpool) && !is.null(pmv$alive)) st$wpool$alive <- pmv$alive
  st$data <- pmwgs$data
  wpool_ctx$data <- pmwgs$data
  st$wpool_ctx <- wpool_ctx
  now <- .emc_reject_counts("pooled")
  delta <- now - st$rejects$pooled
  st$rejects$pooled <- now
  st$stage_rejects$pooled <- st$stage_rejects$pooled + delta
  if (cx$sampler_profile) st$prof$pooled_delta <- delta
  st$pmwgs <- pmwgs
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
                         .emc_steps_now(st, cx, "data"), pmwgs$model, wpool = st$wpool,
                         n_cores = cx$pool_budget$subject,
                         r_cores = cx$pool_budget$likelihood,
                         par_names = cx$par_names, pick = st$hmc_pick,
                         ll_cache_validated = TRUE)
    proposals <- hr$proposals
    st$streams <- hr$streams
    hmc_tasks <- hr$tasks
    st$hmc_plan <- hr$plan
    if (!is.null(st$wpool) && !is.null(hr$pool)) st$wpool$alive <- hr$pool$alive
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
  prof <- cx$sampler_profile
  fill_started <- if (prof) proc.time()[["elapsed"]]
  pmwgs <- st$pmwgs
  pars <- st$pars
  proposals <- st$proposals
  if (is.null(st$ens_pools) && pmwgs$type %in% c("factor", "infnt_factor", "SEM") &&
      !any(pmwgs$nuisance) && !any(as.logical(pmwgs$marginalised_idx), na.rm = TRUE)) {
    assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
    pars <- .emc_ensemble_refresh_scores(pars, pmwgs,
                                         proposals[seq_len(pmwgs$n_pars), , drop = FALSE])
    pmwgs$rng$gibbs <- get(".Random.seed", envir = globalenv())
    st$pars <- pars
  }
  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                proposals = proposals, j = j, n_pars = pmwgs$n_pars,
                                type = pmwgs$type)
  if (!is.null(pmwgs$pooled)) pmwgs$samples$pooled[, j] <- pmwgs$pooled$value
  if (cx$alloc_on) {
    rest_cpu <- .emc_reserved_cost(proc.time()[["elapsed"]] - st$el_iter - st$shear_el -
                                     st$coll_el - st$hmc_el, st$cores_it)
    pmwgs <- .emc_alloc_record_cost(pmwgs, rest_cpu, st$hmc_tasks, moves = list(
      shear = c(cpu = st$shear_cpu, own = st$shear_round_cpu, n = st$shear_rounds),
      collapsed = c(cpu = st$coll_cpu, own = st$coll_round_cpu, n = st$coll_moves),
      hmc = .emc_reserved_cost(st$hmc_el, st$cores_it)))
  }
  if (prof) st$prof$fill <- proc.time()[["elapsed"]] - fill_started
  st$pmwgs <- pmwgs
  .emc_steps_rejects(st)
  invisible(NULL)
}

# The iteration's rejected updates of each group move, into the stage's
# counts (and, for a profiled run, the iteration's own).
.emc_steps_rejects <- function(st) {
  for (kind in c("ensemble", "shear", "collapsed", "gibbs")) {
    now <- .emc_reject_counts(kind)
    if (!is.null(st$prof)) st$prof[[paste0(kind, "_delta")]] <- now - st$rejects[[kind]]
    if (identical(now, st$rejects[[kind]])) next
    st$stage_rejects[[kind]] <- st$stage_rejects[[kind]] + (now - st$rejects[[kind]])
    st$rejects[[kind]] <- now
  }
  invisible(NULL)
}

# The iteration's row of a profiled run.
.emc_step_profile <- function(st, cx, i) {
  j <- cx$start_iter + i
  pf <- st$prof
  pmwgs <- st$pmwgs
  wpool <- st$wpool
  pm_settings <- st$pm_settings
  hmc_tasks <- st$hmc_tasks
  iteration_total <- proc.time()[["elapsed"]] - pf$started
  pp <- if (is.null(pf$pool)) list() else pf$pool
  n_part <- vapply(pm_settings, function(x) {
    if (is.null(x$n_particles)) NA_real_ else as.numeric(x$n_particles)
  }, numeric(1))
  mem <- if (.emc_profile_expensive() &&
             (i == 1L || i == cx$iter || i %% cx$mem_every == 0L)) {
    .emc_profile_memory(if (is.null(wpool)) NULL else c(
      vapply(wpool$jobs, .emc_wpool_job_pid, integer(1)),
      .emc_wpool_job_pid(wpool$template$job)))
  } else list(bytes = NA_real_, method = NA_character_)
  route <- .emc_ll_route(wpool, st$llpool, cx$pool_budget, st$core_budget,
                         .emc_pool_state$ll_route)
  st$profile_rows[[i]] <- do.call(.emc_profile_row, c(list(
    iteration = j,
    stage = cx$stage,
    subjects = pmwgs$n_subjects,
    workers = if (is.null(wpool)) 0L else wpool$n,
    route = route,
    workers_active = pp$workers_active,
    particles_max = if (all(is.na(n_part))) NA_real_ else max(n_part, na.rm = TRUE),
    particles_sum = if (all(is.na(n_part))) NA_real_ else sum(n_part, na.rm = TRUE),
    total = iteration_total,
    gibbs = pf$gibbs,
    group_cache = pf$group_cache,
    recycle = pf$recycle,
    grow = pf$grow,
    startup = if (i == 1L) cx$startup_elapsed else NA_real_,
    particle = pp$elapsed,
    fill = pf$fill,
    shared_serialize = pp$shared_serialize,
    request_send = pp$send,
    response_wait = pp$receive,
    subject_cpu_max = pp$subject_cpu_max,
    subject_cpu_sum = pp$subject_cpu_sum,
    worker_cpu_max = pp$worker_cpu_max,
    worker_cpu_sum = pp$worker_cpu_sum,
    worker_elapsed_max = pp$worker_elapsed_max,
    worker_elapsed_sum = pp$worker_elapsed_sum,
    dispatch_first = pp$dispatch_first,
    dispatch_last = pp$dispatch_last,
    finish_first = pp$finish_first,
    finish_last = pp$finish_last,
    queue_delay_max = pp$queue_delay_max,
    queue_delay_mean = pp$queue_delay_mean,
    shared_bytes = pp$shared_bytes,
    private_bytes = pp$private_bytes,
    wire_bytes = pp$wire_bytes,
    private_memory = mem$bytes,
    private_memory_method = mem$method,
    particle_weight_ess = {
      we <- vapply(pm_settings, function(x) {
        if (is.null(x$weight_ess)) NA_real_ else as.numeric(x$weight_ess)
      }, numeric(1))
      if (all(is.na(we))) NA_real_ else mean(we, na.rm = TRUE)
    },
    fallback_workers = pp$fallback_workers,
    fallback_subjects = pp$fallback_subjects,
    degraded = pp$degraded,
    shear_rounds = st$shear_rounds, shear_cpu = pf$shear_raw,
    shear_elapsed = st$shear_el, shear_build = pf$shear_build, shear_ll = pf$shear_ll,
    ensemble_elapsed = pf$ensemble,
    collapsed_elapsed = st$coll_el, collapsed_moves = st$coll_moves,
    hmc_elapsed = st$hmc_el,
    hmc_task_cpu_max = max(c(0, vapply(hmc_tasks, function(x)
      if (is.null(x) || is.null(x$work)) 0 else x$work, numeric(1)))),
    hmc_attempts = sum(!vapply(hmc_tasks, is.null, logical(1))),
    hmc_accepted = sum(vapply(hmc_tasks, function(x) isTRUE(x$accepted), logical(1))),
    hmc_cpu = sum(vapply(hmc_tasks, function(x) if (is.null(x)) 0 else x$cpu, numeric(1)))),
    .emc_reject_row_fields(pp$rejects, "particle"),
    .emc_reject_row_fields(pf$gibbs_delta, "gibbs"),
    .emc_reject_row_fields(pf$ensemble_delta, "ensemble"),
    .emc_reject_row_fields(pf$shear_delta, "shear"),
    .emc_reject_row_fields(pf$collapsed_delta, "collapsed"),
    .emc_reject_row_fields(pf$pooled_delta, "pooled")))
  invisible(NULL)
}

# The block's iterations from `from` on: the steps in order, each in its
# native form where the block admits those and the form accepts the state.
.emc_steps_run <- function(st, cx, from = 1L) {
  steps <- c(if (cx$verbose_progress) list(progress = .emc_step_progress),
             list(group = .emc_step_group, particle = .emc_step_particle,
                  ensemble = .emc_step_ensemble, shear = .emc_step_shear,
                  collapsed = .emc_step_collapsed),
             if (cx$pooled) list(pooled = .emc_step_pooled),
             list(hmc = .emc_step_hmc, store = .emc_step_store),
             if (cx$sampler_profile) list(profile = .emc_step_profile))
  native <- getOption("emc2.native_steps", TRUE)
  native <- if (!cx$native_ok) rep(FALSE, length(steps)) else
    if (is.character(native)) names(steps) %in% native else
      rep(isTRUE(native), length(steps))
  if (from <= cx$iter)
    emc_iteration_loop_cpp(st, cx, as.integer(from), steps, native, .emc_steps_rejects)
  invisible(NULL)
}

get_stop_criteria <- function(stage, stop_criteria, type){
  if(is.null(stop_criteria)){
    if(stage == "preburn"){
      stop_criteria$iter <- 50
    }
    if(stage == "burn"){
      stop_criteria$mean_gd <- 1.1
      stop_criteria$omit_mpsrf <- TRUE
      if(type != "single"){
        stop_criteria$selection <- c("alpha", "mu")
      } else{
        stop_criteria$selection <- c("alpha")
      }

    }
    if(stage == "adapt"){
      stop_criteria$min_unique <- 150
    }
    if(stage == "sample"){
      stop_criteria$max_gd <- 1.1
      stop_criteria$omit_mpsrf <- TRUE
      if(type != "single"){
        stop_criteria$selection <- c("alpha", "mu")
      } else{
        stop_criteria$selection <- c("alpha")
      }
    }
  }
  if(!is.null(stop_criteria$max_gd) || !is.null(stop_criteria$mean_gd)){
    if(is.null(stop_criteria$selection)) stop_criteria$selection <- c('alpha', 'mu')
    .emc_check_gd_criteria(stop_criteria$selection, stop_criteria$mean_gd, stop_criteria$max_gd)
    if(is.null(stop_criteria$omit_mpsrf)) stop_criteria$omit_mpsrf <- TRUE
  }
  if(stage == "sample"){
    if(!is.null(stop_criteria$max_flat_loc)){
      if(is.null(stop_criteria$flat_selection)) stop_criteria$flat_selection <- c("alpha", "subj_ll", "theta_mu")
      if(is.null(stop_criteria$flat_p1)) stop_criteria$flat_p1 <- 1/3
      if(is.null(stop_criteria$flat_p2)) stop_criteria$flat_p2 <- 1/3
      if(stop_criteria$flat_p1 <= 0 || stop_criteria$flat_p1 > 1) stop("flat_p1 must be in (0, 1].")
      if(stop_criteria$flat_p2 <= 0 || stop_criteria$flat_p2 > 1) stop("flat_p2 must be in (0, 1].")
      if(!all(stop_criteria$flat_selection %in% c("alpha", "subj_ll", "theta_mu", "theta_var"))){
        stop("flat_selection must be one or more of: alpha, subj_ll, theta_mu, theta_var")
      }
    }
    if(!is.null(stop_criteria$max_sample_iter) && stop_criteria$max_sample_iter < 10){
      stop("max_sample_iter should be >= 10.")
    }
  }
  # One default, whether adapt criteria are omitted or given without min_unique.
  if(stage == "adapt" & is.null(stop_criteria$min_unique)) stop_criteria$min_unique <- 150
  if(stage != "adapt" & !is.null(stop_criteria$min_unique)) stop("min_unique only applicable for adapt stage, try min_es instead.")
  return(stop_criteria)
}

#' Fine-Tuned Model Estimation
#'
#' Although typically users will rely on ``fit``, this function can be used for more fine-tuned specification of estimation needs.
#' The function will throw an error if a stage is skipped,
#' the stages have to be run in order ("preburn", "burn", "adapt", "sample").
#' More details can be found in the ``fit`` help files (``?fit``).
#'
#' @param emc An emc object
#' @param stage A string. Indicates which stage is to be run, either `preburn`, `burn`, `adapt` or `sample`
#' @param search_width A double. Tunes target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability.
#' 1 is the default width, increases lead to broader search.
#' @param step_size An integer. After each step, the stopping requirements as
#' specified by `stop_criteria` are checked and proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the stop_criteria.
#' @param verboseProgress Logical. Whether to print a progress bar within each step or not. Will print one progress bar for each chain and only if cores_for_chains = 1.
#' @param fileName A string. If specified will autosave emc at this location on every iteration.
#' @param particle_factor An integer. `particle_factor` multiplied by the square root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain.
#' Parallelizes across participant calculations; for a single-participant fit,
#' the same budget is used across proposal likelihoods. Only available on
#' Linux or Mac OS.
#' For Windows, only parallelization across chains (``cores_for_chains``) is available.
#' @param cores_for_chains An integer. How many cores to use across chains.
#' Defaults to the number of chains. the total number of cores used is equal to ``cores_per_chain`` * ``cores_for_chains``.
#' @param max_tries An integer. How many times should it try to meet the finish
#' conditions as specified by stop_criteria? Defaults to 20. max_tries is ignored if the required number of iterations has not been reached yet.
#' @param n_blocks An integer. Number of blocks. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with a large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See ``?fit``.
#' @param thin A boolean. If `TRUE` will automatically thin the MCMC samples, closely matched to the ESS.
#' Can also be set to a double, in which case 1/thin of the chain will be removed (does not have to be an integer).
#' @param trim A boolean. If `TRUE` will automatically remove redundant samples (i.e. from preburn, burn, adapt).
#' @param r_cores An integer giving an explicit lower bound on the cores used
#' within one participant's proposal-likelihood calculation. This applies to
#' both R and registered C++ likelihoods. It defaults to 1; single-participant
#' fits automatically inherit a larger `cores_per_chain` budget. When nested
#' likelihood work would exceed the current global core budget, the effective
#' value is capped so the outer budget remains authoritative.
#' @param hmc_share A number in \[0, 1\]. The probability that an eligible subject's
#' update in an iteration is a Hamiltonian Monte Carlo transition on its conditional
#' posterior rather than its particle step. Subjects are eligible under a standard or
#' diagonal-gamma hierarchy with one likelihood per subject and a usable fit of it.
#' `0` runs the particle step alone.
#' @export
#' @return An emc object
#' @examples \donttest{
#' # First define a design
#' design_in <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1),
#'                            constants=c(s=log(1)))
#' # Then make the emc, we've omitted a prior here for brevity so default priors will be used.
#' emc <- make_emc(forstmann, design_in, compress = FALSE)
#'
#' # Now for example we can specify that we only want to run the "preburn" phase
#' # for MCMC 10 iterations
#' # emc <- run_emc(emc, stage = "preburn", stop_criteria = list(iter = 10), cores_for_chains = 1)
#'}

run_emc <- function(emc, stage, stop_criteria,
                    search_width = 1, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                    fileName = NULL,particle_factor=50, cores_per_chain = 1,
                    cores_for_chains = length(emc), max_tries = 20, n_blocks = 1,
                    thin = FALSE, trim = TRUE, r_cores=1, rhat_version = "old",
                    hmc_share = 0.5){
  .emc_check_hmc_share(hmc_share)
  emc <- restore_custom_kernel_pointers(emc, quiet = TRUE)
  emc <- restore_duplicates(emc, warn = TRUE)
  # Validate the immutable likelihood configuration once at the fit boundary;
  # per-subject lookup uses the saved fingerprint and small exact keys.
  has_lik_cache <- any(vapply(emc, function(chain)
    !is.null(chain$likelihood_cache_version) ||
      !is.null(chain$subject_fit_cache) ||
      !is.null(chain$laplace_check_cache), logical(1)))
  if (has_lik_cache)
    emc <- .emc_prepare_subject_fit_cache(emc, validate = TRUE)
  emc <- .emc_prepare_chain_rng(emc)
  .emc_ll_route_reset()
  if(Sys.info()[1] == "Windows" & cores_per_chain > 1) stop("only cores_for_chains can be set on Windows")
  if (verbose) message(paste0("Running ", stage, " stage"))
  total_iters_stage <- chain_n(emc)[,stage][1]
  if(stage != "preburn"){
    iter <- stop_criteria[["iter"]] + total_iters_stage
  } else{
    iter <- stop_criteria[["iter"]]
  }
  progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains, verbose, n_blocks = n_blocks, rhat_version = rhat_version)
  emc <- progress$emc
  progress <- progress[!names(progress) == 'emc']
  cur_thin <- ifelse(is.numeric(thin), thin, 1)
  # Stall detection, counted in tries.  It lives here rather than in
  # run_stage() because run_stage() only ever sees ONE try: a counter there
  # restarts every step_size iterations, so the limit could only be expressed
  # in iterations, which is too short for healthy burn-in on slow-moving
  # models.  Preburn is exempt: it proposes from the prior, not
  # an adapted covariance, so long rejection runs there are ordinary.
  stall_limit <- .emc_stall_try_limit()
  check_stall <- !identical(stage, "preburn") && is.finite(stall_limit)
  if (check_stall) {
    frozen_tries <- integer(length(emc))
    last_snap <- lapply(emc, function(chain) .emc_chain_snapshot(chain$samples))
  }
  # Chain runners fork their subject workers from clean templates shared
  # across blocks.
  pool_slots <- if (cores_for_chains >= length(emc)) length(emc)
                else if (cores_for_chains <= 1L) 1L else 0L
  pool_session <- NULL
  on.exit(.emc_wpool_session_stop(pool_session), add = TRUE)
  pool_session <- if (cores_per_chain > 1L && pool_slots > 0L) {
    .emc_wpool_session_start(
      list(data = emc[[1L]]$data,
           model = .emc_wpool_slim_model(emc[[1L]]$model)), pool_slots)
  } else NULL
  block_no <- 0L
  while(!progress$done){
    block_no <- block_no + 1L
    block_started <- proc.time()[["elapsed"]]
    emc <- reset_pm_settings(emc, stage)
    if(trim){
      emc <- fit_remove_samples(emc)
    }
    if (stage == "burn" && block_no == 1L &&
        chain_n(emc)[1L, "burn"] == 0L && get_last_stage(emc) == "preburn") {
      emc <- .emc_burn_mode_start(emc, cores_per_chain * cores_for_chains)
      # The mode starts are deliberately near one another. Give the group
      # hierarchy and subject proposals time to settle before using R-hat.
      # Without a group level the starts are approximate posterior draws, so
      # the ordinary R-hat rule applies from the first block.
      if (is.null(stop_criteria$iter) && emc[[1L]]$type != "single" &&
          any(vapply(emc, function(x) !is.null(x$burn_mode_index), logical(1)))) {
        stop_criteria$iter <- 300L
        iter <- total_iters_stage + stop_criteria$iter
      }
    } else if (stage %in% c("burn", "adapt") ||
               (stage == "sample" && (is.null(emc[[1L]]$sample_proposals_next) ||
                                      .emc_sample_rebuild_due(emc)))) {
      emc <- .emc_laplace_refresh(emc, cores_per_chain * cores_for_chains)
    }
    if(!is.null(progress$n_blocks)) n_blocks <- progress$n_blocks
    prop_started <- proc.time()[["elapsed"]]
    emc <- add_proposals(emc, stage, cores_per_chain*cores_for_chains, n_blocks,
                         block_iter = progress$step_size * max(1, cur_thin),
                         hmc_share = hmc_share)
    prop_elapsed <- proc.time()[["elapsed"]] - prop_started
    last_stage <- get_last_stage(emc)
    if(stage == "preburn"){
      sub_emc <- emc
    } else if(stage != last_stage){
      sub_emc <- subset(emc, filter = chain_n(emc)[1,last_stage]-1, stage = last_stage)
    } else{
      sub_emc <- subset(emc, filter = chain_n(emc)[1,stage] - 1, stage = stage)
    }
    t0 <- Sys.time()
    prepare_elapsed <- proc.time()[["elapsed"]] - block_started - prop_elapsed
    sample_started <- proc.time()[["elapsed"]]
    # Every stop rule waits for the iteration target, and production proposals
    # change only at their rebuild points, so a sample-stage block runs
    # straight to whichever comes first. Each block is a fork whose children
    # copy the session heap, so batching saves forks. Verbose runs need the
    # chain draws back at each step to report fresh diagnostics. Checkpoints
    # and thinning also keep per-step blocks. block_units counts the steps
    # the block stands for.
    block_units <- 1L
    if (stage == "sample" && !verbose && is.null(fileName) && isFALSE(thin) &&
        is.null(stop_criteria$max_sample_iter)) {
      remaining <- iter - progress$total_iters_stage
      to_rebuild <- emc[[1L]]$sample_proposals_next - emc[[1L]]$sample_transitions
      if (length(to_rebuild) != 1L || !is.finite(to_rebuild) || to_rebuild <= 0)
        to_rebuild <- Inf
      block_units <- max(1L, as.integer(floor(min(remaining, to_rebuild) /
                                                progress$step_size)))
      progress$step_size <- progress$step_size * block_units
    }
    stage_iter <- progress$step_size*max(1,cur_thin)
    core_ctl <- .emc_core_ctl(length(sub_emc), cores_per_chain,
                              cores_for_chains = cores_for_chains)
    pool_session <- .emc_wpool_session_check(pool_session)
    for (j in seq_along(sub_emc)) {
      attr(sub_emc[[j]], "emc_pool_slot") <-
        .emc_wpool_session_slot(pool_session, min(j, pool_slots))
    }
    sub_emc <- auto_mclapply(sub_emc, function(sampler, ...) {
                             outer_slot <- .emc_pool_state$session_slot
                             on.exit(.emc_pool_state$session_slot <- outer_slot,
                                     add = TRUE)
                             .emc_pool_state$session_slot <- attr(sampler, "emc_pool_slot")
                             attr(sampler, "emc_pool_slot") <- NULL
                             sampler <- .emc_set_sampler_rng(sampler)
                             run_stages(sampler, ...)
                           }, stage = stage, iter= stage_iter,
                             verbose=verbose,  verboseProgress = verboseProgress,
                             particle_factor=particle_factor,search_width=search_width,
                             n_cores=cores_per_chain, mc.cores = cores_for_chains,
                             r_cores = r_cores, core_ctl = core_ctl,
                             mc.preschedule = FALSE, mc.set.seed = FALSE)
    if (!is.null(core_ctl)) unlink(core_ctl$dir, recursive = TRUE)
    sample_elapsed <- proc.time()[["elapsed"]] - sample_started
    concat_started <- proc.time()[["elapsed"]]
    if(getOption("emc2.print_iteration_duration", FALSE)) { print(Sys.time()-t0) }
    class(sub_emc) <- "emc"
    .emc_iteration_record(sub_emc, block_no)
    if(cores_for_chains > 1) sub_emc <- pointer_reset_wrapper(sub_emc, emc)
    if(stage != 'preburn'){
      if(is.numeric(thin)){
        sub_emc <- subset(sub_emc, stage = c("preburn", "burn", "adapt", "sample"), thin = thin)
      } else if(thin){
        sub_emc <- auto_thin(sub_emc, stage = c("preburn", "burn", "adapt", "sample"))
        cur_thin <- progress$step_size/chain_n(sub_emc)[1,stage]
      }
      emc <- concat_emc(emc, sub_emc, progress$step_size, stage)
      if (stage == "sample")
        for (j in seq_along(emc))
          emc[[j]]$sample_transitions <-
            as.integer(emc[[j]]$sample_transitions + round(stage_iter))
    } else{
      emc <- sub_emc
    }
    if(stage == "sample" && !is.null(stop_criteria$max_sample_iter)){
      emc <- .trim_sample_stage(emc, stop_criteria$max_sample_iter)
    }
    concat_elapsed <- proc.time()[["elapsed"]] - concat_started
    check_started <- proc.time()[["elapsed"]]
    if (block_units > 1L) {
      progress$trys <- progress$trys + block_units - 1L
    }
    progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains,
                               verbose, progress, n_blocks, rhat_version = rhat_version)
    emc <- progress$emc
    progress <- progress[!names(progress) == 'emc']
    if (check_stall) {
      snap <- lapply(emc, function(chain) .emc_chain_snapshot(chain$samples))
      moved <- !mapply(identical, snap, last_snap)
      frozen_tries <- ifelse(moved, 0L, frozen_tries + block_units)
      last_snap <- snap
      if (any(frozen_tries >= stall_limit)) {
        .emc_stall_abort(stage = stage, iteration = chain_n(emc)[1, stage],
                         run_length = max(frozen_tries),
                         chains = which(frozen_tries >= stall_limit))
      }
    }
    check_elapsed <- proc.time()[["elapsed"]] - check_started
    save_started <- proc.time()[["elapsed"]]
    if(!is.null(fileName)){
      emc <- strip_duplicates(emc)
      fileName <- fix_fileName(fileName)
      class(emc) <- "emc"
      save(emc, file = fileName)
      emc <- restore_duplicates(emc)
    }
    save_elapsed <- proc.time()[["elapsed"]] - save_started
    .emc_block_record(
      block = block_no, stage = stage,
      block_iterations = as.integer(stage_iter),
      block_chains = length(emc),
      block_total = proc.time()[["elapsed"]] - block_started,
      block_prepare = prepare_elapsed,
      block_proposals = .emc_proposal_split("chain", prop_elapsed),
      block_eff_proposals = .emc_proposal_split("eff", prop_elapsed),
      block_sample = sample_elapsed,
      block_concat = concat_elapsed,
      block_check_progress = check_elapsed,
      block_save = save_elapsed)
    elapsed <- Sys.time() - t0
    if (verbose) {
      gd_message <- .emc_gd_message(progress$gd_stats)
      ess_message <- ""
      if (stage == "sample" && length(progress$curr_min_es) > 0L) {
        ess_message <- paste0(" | min ESS=", round(progress$curr_min_es))
      }
      flat_message <- ""
      if (length(progress$flat_max) > 0L) {
        flat_message <- paste0(" | flat=", round(progress$flat_max, 3))
      }
      current_iters <- progress$total_iters_stage
      target_iters <- iter
      rem <- estimate_remaining_total_time(
        stage = stage,
        tries_done = progress$trys,
        elapsed_dt = elapsed,
        max_tries = max_tries,
        current_iters = current_iters,
        target_iters = target_iters,
        step_size = progress$step_size
      )
      message(sprintf("[%s | try=%d | iters=%d%s%s%s] Duration: %s - ETA: %s-%s",
                      stage, progress$trys, progress$total_iters_stage,
                      gd_message, ess_message, flat_message,
                      format_duration(elapsed),
                      format_duration(rem$min_time),
                      format_duration(rem$max_time)))
    }
  }
  emc <- strip_duplicates(emc)
  class(emc) <- "emc"
  return(emc)
}

run_stages <- function(sampler, stage = "preburn", iter=0, verbose = TRUE, verboseProgress = TRUE,
                       particle_factor=50, search_width= NULL, n_cores=1, r_cores = 1,
                       core_ctl = NULL)
{
  sampler <- .emc_prepare_sampler_rng(sampler)
  if (!is.null(sampler$rng$gibbs)) {
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  }
  if (!isTRUE(sampler$init)) .emc_ll_route_reset()
  on.exit(.emc_core_release(core_ctl), add = TRUE)
  chain_cores <- .emc_cores_now(core_ctl, n_cores)
  particles <- round(particle_factor*sqrt(sampler$n_pars))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = chain_cores, r_cores = r_cores,
                    total_cores = chain_cores)
  }
  if (iter == 0) return(sampler)
  tune <- list(search_width = search_width)
  timer_dir <- getOption("emc2.chain_timer_dir", NULL)
  t_start <- if (is.null(timer_dir)) NULL else Sys.time()
  sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
                       n_cores = chain_cores, tune = tune, verbose = verbose,
                       verboseProgress = verboseProgress, r_cores = r_cores,
                       core_ctl = core_ctl)
  if (!is.null(t_start)) {
    saveRDS(as.numeric(difftime(Sys.time(), t_start, units = "secs")),
            file.path(timer_dir, paste0("chain_", Sys.getpid(), ".rds")))
  }
  return(sampler)
}

.emc_proposal_split <- function(which, total) {
  got <- .emc_profile_state$proposal_split
  if (is.null(got) || is.na(got[[which]])) return(NA_real_)
  got[[which]]
}

# Whether add_proposals() rebuilds the production proposals in this block.
.emc_sample_rebuild_due <- function(emc) {
  next_at <- emc[[1L]]$sample_proposals_next
  n_sample <- emc[[1L]]$sample_transitions
  !is.null(next_at) && !is.null(n_sample) && n_sample >= next_at
}

# Transition count of the next production rebuild after one at n_sample: the
# stage is refreshed once it holds 100, 400, 1600, ... transitions, so
# adaptation becomes increasingly rare.
.emc_next_sample_rebuild <- function(n_sample) {
  if (n_sample < 100) 100L else 4L * as.integer(n_sample)
}

add_proposals <- function(emc, stage, n_cores, n_blocks, block_iter = 100L,
                          hmc_share = 0.5){
  prof <- .emc_profile_enabled()
  split <- c(chain = 0, eff = 0)
  proposal_started <- if (prof) proc.time()[["elapsed"]] else NA_real_
  # Production proposal refreshes follow the transition counter, independent
  # of thinning and trimming. Mixture weights and scales remain fixed.
  sample_window <- NULL
  if (stage == "sample") {
    n_stored <- sum(emc[[1]]$samples$stage[seq_len(emc[[1]]$samples$idx)] == "sample")
    # A checkpoint without the counter resumes from the stored count.
    n_sample <- emc[[1]]$sample_transitions
    if (is.null(n_sample)) {
      n_sample <- as.integer(n_stored)
      for (j in seq_along(emc)) emc[[j]]$sample_transitions <- n_sample
    }
    next_at <- emc[[1]]$sample_proposals_next
    if (!is.null(next_at) && n_sample < next_at) {
      if (prof) .emc_profile_state$proposal_split <- split
      return(emc)
    }
    if (!is.null(next_at) && n_stored > 0)
      sample_window <- which(emc[[1]]$samples$stage == "sample" &
                               seq_along(emc[[1]]$samples$stage) <= emc[[1]]$samples$idx)
  }
  if(stage != "preburn"){
    emc <- create_chain_proposals(
      emc, samples_idx = sample_window, do_block = stage != "sample",
      stage = stage)
    emc <- .emc_check_likelihood_fits(emc, n_cores)
    # The shear plan uses only checked fits, so it is built after the checks.
    emc <- .emc_shear_set_plans(emc, sample_window, stage, n_cores)
    emc <- .emc_collapsed_set_plans(emc, sample_window, stage, n_cores)
    if(!is.null(n_blocks)){
      if(n_blocks > 1){
        # Each chain blocks on its own covariances.
        for(i in seq_along(emc)){
          attr(emc[[i]]$data, "components") <- sub_blocking(emc[i], n_blocks)
        }
      }
    }
  }
  eff_started <- if (prof) proc.time()[["elapsed"]] else NA_real_
  if(stage == "adapt"){
    emc <- create_eff_proposals(emc, n_cores, from_burn = TRUE)
  }
  if (stage == "burn") {
    for (j in seq_along(emc)) emc[[j]]$burn_eff_ready <- FALSE
    if (chain_n(emc)[1L, "burn"] >= 200L && !any(emc[[1L]]$nuisance)) {
      emc <- tryCatch(create_eff_proposals(emc, n_cores, from_burn = TRUE),
                      error = function(e) {
                        warning("Burn efficient proposal unavailable: ",
                                conditionMessage(e), call. = FALSE)
                        emc
                      })
      chain_ready <- function(chain) {
        isTRUE(chain$burn_eff_ready) &&
          length(chain$eff_mu) == chain$n_subjects &&
          length(chain$eff_var) == chain$n_subjects &&
          all(vapply(chain$eff_mu, function(x)
            length(x) == chain$n_pars && all(is.finite(x)), logical(1))) &&
          all(vapply(chain$eff_var, function(x)
            is.matrix(x) && nrow(x) == chain$n_pars &&
              ncol(x) == chain$n_pars &&
              .emc_empirical_covariance_ok(x), logical(1)))
      }
      # Each chain switches when its own proposals are ready.
      for (j in seq_along(emc)) emc[[j]]$burn_eff_ready <- chain_ready(emc[[j]])
    }
  }
  if(stage == "sample"){
    emc <- create_eff_proposals(emc, n_cores)
    next_at <- .emc_next_sample_rebuild(n_sample)
    for (j in seq_along(emc)) emc[[j]]$sample_proposals_next <- next_at
    block_iter <- next_at - n_sample
  }
  if (prof) split[["eff"]] <- proc.time()[["elapsed"]] - eff_started
  # Shear rounds, collapsed rate and subject HMC rates (R/move_allocation.R)
  # on the same rebuild clock, after the likelihood fits are checked and the
  # shear plans built; block_iter transitions run before the next.
  emc <- tryCatch(.emc_alloc_rebuild(emc, stage, n_cores, block_iter, hmc_share),
    error = function(e) {
      cls <- .emc_classify_failure(e)
      .emc_reject_record(e, "allocator")
      if (identical(.emc_failure_action(cls), "abort"))
        .emc_failure_abort(cls, "allocator", e)
      emc
    })
  if (prof) {
    split[["chain"]] <- proc.time()[["elapsed"]] - proposal_started - split[["eff"]]
    .emc_profile_state$proposal_split <- split
  }
  return(emc)
}

check_progress <- function (emc, stage, iter, stop_criteria,
                            max_tries, step_size, n_cores, verbose, progress = NULL,
                            n_blocks, rhat_version = "old", truncate_gd = TRUE)
{
  min_es <- stop_criteria$min_es
  if(is.null(min_es)) min_es <- 0
  selection <- stop_criteria$selection
  min_unique <- stop_criteria$min_unique
  total_iters_stage <- chain_n(emc)[, stage][1]
  if (is.null(progress)) {
    iters_total <- 0
    trys <- 0
  }
  else {
    iters_total <- progress$iters_total + step_size
    trys <- progress$trys + 1
  }
  gd <- check_gd(
    emc = emc,
    stage = stage,
    max_gd = stop_criteria[["max_gd"]],
    mean_gd = stop_criteria[["mean_gd"]],
    trys = trys,
    verbose = FALSE,
    iter = total_iters_stage,
    selection = selection,
    omit_mpsrf = stop_criteria[["omit_mpsrf"]],
    n_blocks = n_blocks,
    rhat_version = rhat_version,
    truncate = truncate_gd,
    n_cores = n_cores
  )
  iter_done <- ifelse(is.null(iter) || length(iter) == 0, TRUE, total_iters_stage >= iter)
  if (min_es == 0) {
    es_done <- TRUE
  } else if (total_iters_stage != 0) {
    class(emc) <- "emc"
    curr_min_es <- Inf
    for(select in selection){
      es <- ess_summary(emc, selection = select, stage = stage,
                         stat_only = TRUE, n_cores = n_cores)
      curr_min_es <- min(c(es, curr_min_es))
    }
    es_done <- ifelse(!emc[[1]]$init, FALSE, curr_min_es >
                        min_es)
  }
  else {
    es_done <- FALSE
  }
  trys_done <- ifelse(is.null(max_tries), FALSE, trys >= max_tries)
  if (stage == "adapt") {
    samples_merged <- merge_chains(emc)
    test_samples <- extract_samples(samples_merged, stage = "adapt",
                                    samples_merged$samples$idx, n_chains = length(emc))
    adapted <- test_adapted(emc[[1]], test_samples,
                            min_unique, n_cores, verbose)
    if (!isTRUE(adapted) && trys_done && iter_done) {
      # Out of tries. Sampling needs the efficient proposals, so only continue
      # if they can be built at all; the unique-value threshold is a quality
      # target, not a requirement.
      reason <- attr(adapted, "reason")
      if (identical(reason, "min_unique")) {
        buildable <- test_adapted(emc[[1]], test_samples, 0, n_cores, FALSE)
        if (isTRUE(buildable)) {
          warning("Adaptation reached max_tries before every subject had ",
                  min_unique, " unique values; efficient proposals are built ",
                  "from fewer draws. Consider more adapt iterations.", call. = FALSE)
          adapted <- TRUE
        } else reason <- attr(buildable, "reason")
      }
      if (!isTRUE(adapted)) {
        stop("Adaptation did not complete within max_tries: efficient proposals ",
             "could not be built (", reason, "). Increase the adapt iterations ",
             "or max_tries, or check the model for unidentified parameters.",
             call. = FALSE)
      }
    }

  }
  else {
    adapted <- TRUE
  }
  flat_done <- TRUE
  flat_out  <- NULL
  if(stage == "sample" && !is.null(stop_criteria$max_flat_loc)){
    flat_out <- .check_flatness_stop(
      emc = emc,
      max_flat_loc = stop_criteria$max_flat_loc,
      flat_selection = stop_criteria$flat_selection,
      flat_p1 = stop_criteria$flat_p1,
      flat_p2 = stop_criteria$flat_p2,
      stage = "sample"
    )
    flat_done <- flat_out$flat_done
  }

  done <- (es_done & iter_done & gd$gd_done & adapted & flat_done) | (trys_done & iter_done)
  if(es_done & gd$gd_done & adapted & flat_done & !iter_done){
    step_size <- min(step_size, abs(iter - total_iters_stage))[1]
  }
  if (trys_done & iter_done) {
    if(!(es_done & gd$gd_done & adapted & flat_done)){
      warning("Max tries reached. If this happens in burn-in while trying to get
            gelman diagnostics small enough, you might have a particularly hard model.
            Make sure your model is well specified. If so, you can run adapt and
            sample, if run for long enough, sample usually converges eventually.")
    }
  }

  log_entry <- list(
    stage     = stage,
    iter      = total_iters_stage,
    try       = trys,
    gd        = list(
      values    = gd$gd,
      stats     = gd$gd_stats,
      selection = selection,
      version   = rhat_version
    ),
    flat      = if (!is.null(flat_out)) list(
      max_flat         = flat_out$max_flat,
      max_by_selection = flat_out$max_by_selection,
      stats_by_selection = flat_out$stats_by_selection
    ) else NULL
  )
  out_emc <- gd$emc
  out_emc[[1]]$convergence_log <- c(out_emc[[1]]$convergence_log, list(log_entry))

  return(list(emc = out_emc, done = done, step_size = step_size,
              trys = trys, n_blocks = gd$n_blocks, gd = gd$gd,
              gd_stats = gd$gd_stats,
              total_iters_stage = total_iters_stage,
              curr_min_es = if (min_es > 0 && total_iters_stage != 0) curr_min_es else NULL,
              flat_max = if (!is.null(flat_out)) flat_out$max_flat else NULL))
}

# R-hat rules per selection. A scalar mean_gd or max_gd applies to every
# selection on its own; a named vector applies only to the selections it names,
# so e.g. mean_gd = c(alpha = 1.1), max_gd = c(mu = 1.1) mixes statistics.
.emc_gd_rules <- function(selection, mean_gd = NULL, max_gd = NULL) {
  pick <- function(x, s) {
    if (is.null(x)) return(NULL)
    if (is.null(names(x))) return(x[[1L]])
    if (s %in% names(x)) x[[s]] else NULL
  }
  rules <- lapply(selection, function(s) {
    out <- list()
    if (!is.null(m <- pick(mean_gd, s))) out$mean <- m
    if (!is.null(m <- pick(max_gd, s))) out$max <- m
    out
  })
  names(rules) <- selection
  rules[lengths(rules) > 0L]
}

.emc_check_gd_criteria <- function(selection, mean_gd, max_gd) {
  for (nm in c("mean_gd", "max_gd")) {
    x <- get(nm)
    if (is.null(x)) next
    if (!is.numeric(x) || anyNA(x)) stop(nm, " must be numeric.")
    if (is.null(names(x))) {
      if (length(x) != 1L) stop(nm, " must be a single number or a vector named by selection.")
    } else if (!all(names(x) %in% selection)) {
      stop(nm, " names ", paste(setdiff(names(x), selection), collapse = ", "),
           " are not in selection (", paste(selection, collapse = ", "), ").")
    }
  }
  invisible(TRUE)
}

# One row per rule: the statistic of that selection's R-hats against its threshold.
.emc_gd_stats <- function(gds, rules) {
  do.call(rbind, lapply(names(rules), function(s) {
    g <- gds[[s]]
    do.call(rbind, lapply(names(rules[[s]]), function(stat) {
      value <- if (!length(g) || !all(is.finite(g))) Inf else
        if (stat == "mean") mean(g) else max(g)
      data.frame(selection = s, statistic = stat, value = value,
                 threshold = rules[[s]][[stat]], stringsAsFactors = FALSE)
    }))
  }))
}

.emc_gd_message <- function(stats) {
  if (is.null(stats) || !nrow(stats)) return("")
  paste0(" | ", paste(sprintf("%s %s Rhat=%.3f", stats$selection, stats$statistic,
                               stats$value), collapse = ", "))
}

check_gd <- function(emc, stage, max_gd, mean_gd, omit_mpsrf, trys, verbose,
                     selection, iter, n_blocks = 1, rhat_version = "old",
                     truncate = TRUE, n_cores = 1L)
{
  get_gds <- function(emc, selection) {
    gds <- lapply(selection, function(select) {
      gd <- unlist(gd_summary.emc(emc, selection = select, stage = stage,
                                  omit_mpsrf = omit_mpsrf, stat = NULL,
                                  version = rhat_version, n_cores = n_cores))
      gd[is.na(gd)] <- Inf
      gd
    })
    names(gds) <- selection
    gds
  }
  if(is.null(max_gd) & is.null(mean_gd)) return(list(gd_done = TRUE, emc = emc, gd = NULL))
  if(!emc[[1]]$init | !stage %in% emc[[1]]$samples$stage)
    return(list(gd_done = FALSE, emc = emc, gd = NULL))
  if(is.null(omit_mpsrf)) omit_mpsrf <- TRUE
  rules <- .emc_gd_rules(selection, mean_gd, max_gd)
  if (!length(rules)) return(list(gd_done = TRUE, emc = emc, gd = NULL))
  # Selections without a rule are not computed.
  gds <- get_gds(emc, names(rules))
  stats <- .emc_gd_stats(gds, rules)
  # How far the worst rule is from passing; used to choose between the full
  # and the truncated stage.
  excess <- function(st) max(st$value / st$threshold)
  ok_gd <- all(stats$value < stats$threshold)
  if (!ok_gd && truncate) {
    n_remove <- round(chain_n(emc)[,stage][1]/3)
    samplers_short <- subset.emc(emc, filter=n_remove,stage=stage, keep_stages = TRUE)
    if (!is(samplers_short, "try-error")) {
      gds_short <- get_gds(samplers_short, names(rules))
      stats_short <- .emc_gd_stats(gds_short, rules)
      if (isTRUE(excess(stats_short) < excess(stats))) {
        gds <- gds_short
        stats <- stats_short
        emc <- samplers_short
      }
    }
    ok_gd <- all(stats$value < stats$threshold)
  }
  if(verbose) message(sub("^ \\| ", "", .emc_gd_message(stats)))
  if (!is.null(rules$alpha)) {
    emc <- set_tune_ess(emc, gds$alpha, rules$alpha$mean, rules$alpha$max)
  }
  class(emc) <- "emc"
  return(list(gd = unlist(gds, use.names = FALSE), gd_stats = stats,
              gd_done = ok_gd, emc = emc))
}

set_tune_ess <- function(emc, alpha_gd = NULL, mean_gd = NULL, max_gd = NULL){
  if(is.null(alpha_gd)){
    return(emc)
  }
  alpha_gd <- matrix(alpha_gd, emc[[1]]$n_subjects)
  if(!is.null(mean_gd)){
    mean_alpha_ok <- rowMeans(alpha_gd) < 1+(mean_gd-1)*.5
  } else{
    mean_alpha_ok <- rep(T, nrow(alpha_gd))
  }
  if(!is.null(max_gd)){
    max_alpha_ok <- alpha_gd[cbind(1:nrow(alpha_gd), max.col(alpha_gd, ties.method="first"))] < 1+(max_gd-1)*.5
  } else{
    max_alpha_ok <- rep(T, nrow(alpha_gd))
  }
  alpha_ok <- max_alpha_ok & mean_alpha_ok
  emc <- lapply(emc, function(x){
    pm_settings <- attr(x$samples, "pm_settings")
    pm_settings <- mapply(pm_settings, alpha_ok, FUN = function(y, z){
      for(i in 1:length(y)){
        y[[i]]$gd_good <- z
      }
      return(list(y))
    })
    attr(x$samples, "pm_settings") <- pm_settings
    return(x)
  })
  return(emc)
}


# Chains share nothing between them: each builds its proposals from its own
# draws, so chains are independent replicates, which R-hat assumes.
create_eff_proposals <- function(emc, n_cores, from_burn = FALSE){
  if (length(emc) > 1L) {
    for (i in seq_along(emc))
      emc[[i]] <- create_eff_proposals(list(emc[[i]]), n_cores, from_burn)[[1L]]
    return(emc)
  }
  samples_merged <- merge_chains(emc)
  if (from_burn && !any(emc[[1]]$nuisance)) {
    # Adapt tunes the production kernel, so it needs an efficient proposal
    # before any adapt draws exist: use the last 250 warm-up draws per chain.
    st <- samples_merged$samples$stage
    pool <- which(st %in% c("burn", "adapt") &
                    seq_along(st) <= samples_merged$samples$idx)
    keep <- tail(pool, 250L)
    test_samples <- filtered_samples(samples_merged, keep, type = samples_merged$type)
  } else
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 750, n_chains = length(emc))
  type <- emc[[1]]$type
  components <- attr(emc[[1]]$data, "components")
  for(i in 1:length(emc)){
    iteration = round(test_samples$iteration * i/length(emc))
    n_pars <- emc[[1]]$n_pars
    nuisance <- emc[[1]]$nuisance
    n_subjects <- emc[[1]]$n_subjects
    eff_mu <- matrix(0, nrow = n_pars, ncol = n_subjects)
    eff_var <- array(0, dim = c(n_pars, n_pars, n_subjects))
    # Regression of each conditional mean on the group mean it conditions on,
    # and the value used; new_particle() moves the mean with the group state.
    eff_B <- array(0, dim = c(n_pars, n_pars, n_subjects))
    eff_given <- matrix(NA_real_, nrow = n_pars, ncol = n_subjects)
    has_B <- TRUE
    place <- function(conds, rows) {
      for (sub in seq_len(n_subjects)) {
        cd <- conds[[sub]]
        eff_mu[rows, sub] <<- cd$eff_mu
        eff_var[rows, rows, sub] <<- cd$eff_var
        if (is.null(cd$B) || is.null(cd$given)) {
          has_B <<- FALSE
        } else {
          eff_B[rows, rows, sub] <<- cd$B
          eff_given[rows, sub] <<- cd$given
        }
      }
    }
    for(comp in unique(components)){
      idx <- comp == components
      nuis_idx <- nuisance[idx]
      if(any(nuis_idx)){
        type <- samples_merged$sampler_nuis$type
        conditionals <- .emc_rebuild_lapply(seq_len(n_subjects), function(sub)
          get_conditionals(sub, samples = test_samples,
                           n_pars = sum(idx[!nuisance]), iteration = iteration,
                           idx = idx[!nuisance], type = type),
          n_cores, cost_key = "eff_conditionals")
        conditionals_nuis <- .emc_rebuild_lapply(seq_len(n_subjects), function(sub)
          get_conditionals(sub, samples = test_samples$nuisance,
                           n_pars = sum(idx[nuisance]), iteration = iteration,
                           idx = idx[nuisance], type = type),
          n_cores, cost_key = "eff_conditionals_nuisance")
        place(conditionals, which(idx & !nuisance))
        place(conditionals_nuis, which(idx & nuisance))
      } else{
        # From warm-up draws the conditional can still be singular. Catch that
        # inside each forked job: a condition raised in a child reaches the
        # parent's handlers (and, under testthat's parallel reporter, its pipe).
        cond_fun <- if (from_burn) function(...) {
          tryCatch(get_conditionals(...), warning = function(w) NULL,
                   error = function(e) NULL)
        } else get_conditionals
        # The conditionals are deterministic; cheap ones run without forking.
        conditionals <- .emc_rebuild_lapply(seq_len(n_subjects), function(sub)
          cond_fun(sub, samples = test_samples,
                   n_pars = sum(idx[!nuisance]), iteration = iteration,
                   idx = idx[!nuisance], type = type),
          n_cores, cost_key = "eff_conditionals")
        if (from_burn && any(vapply(conditionals, is.null, logical(1)))) {
          # Not buildable yet: the adapt block uses the chain proposal in the
          # efficient slot, which keeps the production mixture's shape.
          for (j in seq_along(emc)) {
            emc[[j]]$eff_mu <- emc[[j]]$chains_mu
            emc[[j]]$eff_var <- emc[[j]]$chains_var
          }
          return(emc)
        }
        place(conditionals, which(idx & !nuisance))
      }

    }

    eff_mu <- lapply(seq_len(n_subjects), function(sub) {
      m <- eff_mu[, sub]
      if (has_B && all(is.finite(eff_B[, , sub]))) {
        attr(m, "B") <- matrix(eff_B[, , sub], n_pars, n_pars)
        attr(m, "given") <- eff_given[, sub]
      }
      m
    })
    eff_var <- lapply(apply(eff_var, 3, identity, simplify = F), function(v) {
      if (!.emc_empirical_covariance_ok(v)) .emc_regularize_covariance(v) else v
    })
    emc[[i]]$eff_mu <- eff_mu
    emc[[i]]$eff_var <- eff_var
    if (from_burn) emc[[i]]$burn_eff_ready <- TRUE
  }
  return(emc)
}


sub_blocking <- function(emc, n_blocks){
  covs <- lapply(emc, FUN = function(x){return(x$chains_var)})
  out <- array(0, dim = dim(covs[[1]][[1]]))
  n_summed <- 0L
  for(i in seq_along(covs)){
    cov_tmp <- covs[[i]]
    for(j in seq_along(cov_tmp)){
      out <- out + cov2cor(cov_tmp[[j]])
      n_summed <- n_summed + 1L
    }
  }
  if(n_summed == 0L) stop("sub_blocking: no chain covariance matrices to block on")
  out <- out/n_summed
  shared_ll_idx <- attr(emc[[1]]$data, "shared_ll_idx")
  min_comp <- 0
  components <- c()
  for(ll in unique(shared_ll_idx)){
    idx <- ll == shared_ll_idx
    n_members <- sum(idx)
    if(n_members <= 1L){
      sub_comps <- stats::setNames(rep(min_comp + 1L, n_members),
                                   dimnames(out)[[1]][idx])
    } else {
      distance <- as.dist(1 - abs(out[idx, idx, drop = FALSE]))
      clusts <- hclust(distance)
      sub_comps <- min_comp + cutree(clusts, k = min(n_blocks, n_members))
    }
    min_comp <- max(sub_comps)
    components <- c(components, sub_comps)
  }
  return(components)
}

# Whether S is an exactly symmetric, finite covariance whose eigenvalues and
# reciprocal condition number are safely positive.
.emc_empirical_covariance_ok <- function(S, subject = NULL, context = NULL) {
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) != ncol(S) || !nrow(S))
    return(FALSE)
  if (is.double(S)) {
    ok <- emc_covariance_ok_cpp(S)
    if (!is.null(ok)) return(ok)
  }
  if (!all(is.finite(S)) || !all(S == t(S))) return(FALSE)
  values <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) NULL)
  if (is.null(values) || !all(is.finite(values)) || min(values) < 1e-8)
    return(FALSE)
  rcond <- tryCatch(1 / kappa(S, exact = FALSE), error = function(e) NA_real_)
  is.finite(rcond) && rcond >= 1e-12
}

.emc_chain_proposal_window <- function(emc, samples_idx = NULL, stage) {
  if (is.null(samples_idx)) {
    idx <- emc[[1L]]$samples$idx
    samples_idx <- seq.int(round(idx - min(250, idx / 1.5)), idx)
  }
  if (stage == "burn") {
    starts <- vapply(emc, function(x) {
      if (is.null(x$burn_mode_index)) 1L else as.integer(x$burn_mode_index)
    }, integer(1))
    samples_idx <- samples_idx[samples_idx >= max(starts)]
  }
  samples_idx
}

create_chain_proposals <- function(emc, samples_idx = NULL, do_block = TRUE,
                                   stage){
  n_subjects <- emc[[1]]$n_subjects
  n_chains <- length(emc)
  n_pars <- emc[[1]]$n_pars
  stage <- match.arg(stage, c("burn", "adapt", "sample"))
  samples_idx <- .emc_chain_proposal_window(emc, samples_idx, stage)
  # Posterior moments of each subject over the chain's own window: the plain
  # mean and covariance. The covariance feeds an unscaled independence
  # proposal, whose efficiency depends on the correlations, so they are not
  # shrunk; a window that is not positive definite is regularised below.
  win_moments <- function(a) list(
    mu = lapply(seq_len(n_subjects), function(sub)
      rowMeans(matrix(a[, sub, ], nrow = n_pars))),
    cov = lapply(seq_len(n_subjects), function(sub)
      stats::cov(t(matrix(a[, sub, ], nrow = n_pars)))))
  for(j in 1:n_chains){
    wm <- win_moments(emc[[j]]$samples$alpha[, , samples_idx, drop = FALSE])
    win_mu <- wm$mu; win_cov <- wm$cov
    block_idx <- block_variance_idx(attr(emc[[j]]$data, "components"))
    old_chains_var <- emc[[j]]$chains_var
    old_chains_mu <- emc[[j]]$chains_mu
    old_ensemble_q <- emc[[j]]$ensemble_q
    burn_mode_var <- emc[[j]]$burn_mode_var
    burn_mode_mu <- emc[[j]]$burn_mode_mu
    # Subjects whose Laplace approximation is poor (a ridge or a support
    # boundary) take the covariance of their draws once there is a usable one.
    mode_stats <- emc[[j]]$burn_mode_stats
    laplace_poor <- .emc_laplace_poor(mode_stats)
    chains_var <- vector("list", n_subjects)
    chains_mu <- vector("list", n_subjects)
    window_informative <- window_ok <- logical(n_subjects)
    for(sub in 1:n_subjects){
      emp_covs <- win_cov[[sub]]
      chains_mu[[sub]] <- win_mu[[sub]]
      if(do_block) emp_covs[block_idx] <- 0
      informative <- is.matrix(emp_covs) &&
        nrow(emp_covs) == n_pars && ncol(emp_covs) == n_pars &&
        all(is.finite(emp_covs)) &&
        any(diag(emp_covs) > 0)
      raw_ok <- informative && .emc_empirical_covariance_ok(emp_covs, subject = sub)
      window_informative[sub] <- informative
      window_ok[sub] <- raw_ok
      old_ok <- length(old_chains_var) >= sub &&
        is.matrix(old_chains_var[[sub]]) &&
        nrow(old_chains_var[[sub]]) == n_pars &&
        ncol(old_chains_var[[sub]]) == n_pars &&
        .emc_empirical_covariance_ok(old_chains_var[[sub]], subject = sub)
      keep_old <- stage == "sample" && !raw_ok && old_ok
      if(informative && !raw_ok && !keep_old){
        emp_covs <- .emc_regularize_covariance(emp_covs)
      }
      if(informative && !keep_old &&
         (raw_ok || .emc_empirical_covariance_ok(emp_covs, subject = sub))){
        chains_var[[sub]] <- emp_covs
      } else if(old_ok){
        # An empty or unusable window should not replace a working proposal.
        chains_var[[sub]] <- old_chains_var[[sub]]
        if(length(old_chains_mu) >= sub &&
           length(old_chains_mu[[sub]]) == n_pars &&
           all(is.finite(old_chains_mu[[sub]]))) {
          chains_mu[[sub]] <- old_chains_mu[[sub]]
        }
      }
      if(length(chains_mu[[sub]]) != n_pars ||
         any(!is.finite(chains_mu[[sub]])))
        chains_mu[[sub]] <- win_mu[[sub]]
      if (length(burn_mode_var) >= sub &&
          is.matrix(burn_mode_var[[sub]]) &&
          !(isTRUE(laplace_poor[sub]) && raw_ok)) {
        mode_cov <- burn_mode_var[[sub]]
        if (do_block) mode_cov[block_idx] <- 0
        if (.emc_empirical_covariance_ok(mode_cov, subject = sub)) {
          chains_var[[sub]] <- mode_cov
          chains_mu[[sub]] <- if (stage != "burn" && length(burn_mode_mu) >= sub &&
                                  length(burn_mode_mu[[sub]]) == n_pars) {
            as.numeric(burn_mode_mu[[sub]])
          } else as.numeric(emc[[j]]$samples$alpha[, sub, emc[[j]]$samples$idx])
        }
      }
    }
    null_idx <- sapply(chains_var, is.null)
    if(all(null_idx)){
      mean_chains_var <- diag(n_pars) * .5
    } else{
      mean_chains_var <- Reduce(`+`, chains_var[!null_idx]) / sum(!null_idx)
    }
    if(!.emc_empirical_covariance_ok(mean_chains_var, context = "mean empirical proposal covariance")){
      mean_chains_var <- .emc_regularize_covariance(mean_chains_var)
    }
    if(!.emc_empirical_covariance_ok(mean_chains_var, context = "mean empirical proposal covariance"))
      mean_chains_var <- diag(n_pars) * .5
    if(any(null_idx)){
      for(q in 1:n_subjects){
        if(null_idx[q]) chains_var[[q]] <- mean_chains_var
      }
    }

    # A subject without a usable current or earlier covariance uses the pooled
    # one. Keep its adapted multiplier: resetting it to one at the sample
    # boundary can freeze a wide fit for the entire production stage.
    new_prop_var <- vapply(chains_var, function(v) mean(diag(v)), numeric(1))
    bad_var <- !is.finite(new_prop_var) | new_prop_var <= 0
    new_prop_var[bad_var] <- .5
    if(stage != "sample"){
      old_prop_var <- attr(emc[[j]], "prop_var")
    prop_var_ratio <- if(is.null(old_prop_var)) rep(1, n_subjects) else {
        if(length(old_prop_var) == 1L) old_prop_var <- rep(old_prop_var, n_subjects)
        if(length(old_prop_var) != n_subjects) rep(1, n_subjects) else {
          r <- old_prop_var / new_prop_var
          r[!is.finite(r) | r <= 0] <- 1
          r
        }
      }
      prop_var_ratio[null_idx | bad_var] <- 1
      emc[[j]] <- update_epsilon_scale(emc[[j]], prop_var_ratio)
    }
    attr(emc[[j]], "prop_var") <- new_prop_var
    emc[[j]]$chains_var <- chains_var
    emc[[j]]$chains_mu <- chains_mu
    emc[[j]]$burn_mode_var <- NULL
    # Pool proposal for the ensemble group update: the same window moments
    # (in burn, the Laplace proposal mode and a widened Laplace covariance).
    # The mean and covariance come from the same source, so a poor Laplace fit
    # is replaced by the window moments as a pair.
    ens_q <- lapply(seq_len(n_subjects), function(sub) {
      v <- win_cov[[sub]]
      ens_v <- v
      if (do_block) v[block_idx] <- 0
      informative <- window_informative[sub]
      raw_ok <- window_ok[sub]
      if (stage == "burn" && length(burn_mode_var) >= sub &&
          is.matrix(burn_mode_var[[sub]]) &&
          length(burn_mode_mu) >= sub &&
          length(burn_mode_mu[[sub]]) == n_pars &&
          !(isTRUE(laplace_poor[sub]) && raw_ok))
        return(list(mu = as.numeric(burn_mode_mu[[sub]]),
                    var = 3 * burn_mode_var[[sub]]))
      m <- win_mu[[sub]]
      old_mu_ok <- length(old_ensemble_q$mu) >= sub &&
        length(old_ensemble_q$mu[[sub]]) == n_pars &&
        all(is.finite(old_ensemble_q$mu[[sub]]))
      old_ok <- old_mu_ok && length(old_ensemble_q$var) >= sub &&
        is.matrix(old_ensemble_q$var[[sub]]) &&
        nrow(old_ensemble_q$var[[sub]]) == n_pars &&
        ncol(old_ensemble_q$var[[sub]]) == n_pars &&
        .emc_empirical_covariance_ok(old_ensemble_q$var[[sub]], subject = sub)
      old_q <- if (old_ok) list(mu = old_ensemble_q$mu[[sub]],
                                var = old_ensemble_q$var[[sub]])
      if (!all(is.finite(m))) {
        if (old_ok) return(old_q)
        return(list(mu = chains_mu[[sub]], var = chains_var[[sub]]))
      }
      ens_informative <- is.matrix(ens_v) && nrow(ens_v) == n_pars &&
        ncol(ens_v) == n_pars && all(is.finite(ens_v)) && any(diag(ens_v) > 0)
      ens_ok <- ens_informative && .emc_empirical_covariance_ok(ens_v, subject = sub)
      if (stage == "sample" && !ens_ok && old_ok) return(old_q)
      if (ens_informative && !ens_ok) ens_v <- .emc_regularize_covariance(ens_v)
      if (ens_informative && .emc_empirical_covariance_ok(ens_v, subject = sub))
        return(list(mu = m, var = ens_v))
      if (old_ok) return(old_q)
      list(mu = m, var = chains_var[[sub]])
    })
    emc[[j]]$ensemble_q <- list(mu = lapply(ens_q, `[[`, "mu"),
                                var = lapply(ens_q, `[[`, "var"))
    emc[[j]]$burn_mode_mu <- NULL
  }

  return(emc)
}

reset_pm_settings <- function(emc, stage){
  # Production does not tune particle settings. Keep the mixture learned in
  # adapt; check_mix() adds the new efficient-proposal component in run_stage().
  if(stage == "sample") return(emc)
  if(stage != get_last_stage(emc) || stage == "burn"){
    for(i in 1:length(emc)){
      pm_settings <- attr(emc[[i]]$samples, "pm_settings")
      attr(emc[[i]]$samples, "pm_settings") <- lapply(pm_settings, function(x){
        for(i in 1:length(x)){
          x[[i]]$proposal_counts <- rep(0, length(x[[i]]$proposal_counts))
          x[[i]]$acc_counts <- rep(0, length(x[[i]]$proposal_counts))
          x[[i]]$jump_counts <- NULL
          # Every burn block is a fresh adaptation window.  Leaving the
          # iteration counter above tune$n0 after the first failed block makes
          # the next block adapt from a single observation immediately.
          if(stage == "burn" || stage != get_last_stage(emc)){
            x[[i]]$iter <- 25
            x[[i]]$mix <- NULL
          }
          # Burn has met its stop rule, and sample freezes the particle
          # count, so adapt must tune it. gd_good uses a stricter R-hat than
          # the burn rule and is not rechecked in adapt, so it can otherwise
          # stay FALSE through adapt and sample.
          if(stage == "adapt") x[[i]]$gd_good <- TRUE
        }
        return(x)
      })
    }
  }
  return(emc)
}

update_epsilon_scale <- function(pmwgs, prop_var_ratio = NULL){
  scale_bounds <- c(.25, 4)
  epsilon_bounds <- c(.05, 5)
  pm_settings <- attr(pmwgs$samples, "pm_settings")
  n_subjects <- length(pm_settings)
  # prop_var_ratio is a per-subject variance ratio (length 1 is recycled),
  # while epsilon multiplies a Cholesky
  # factor in particle_draws().  Convert variance scaling to standard-
  # deviation scaling before applying it, and bound the one-block jump.
  scale <- rep(1, n_subjects)
  if(!is.null(prop_var_ratio)){
    r <- rep_len(as.numeric(prop_var_ratio), n_subjects)
    ok <- is.finite(r) & r > 0
    scale[ok] <- pmin(scale_bounds[2], pmax(scale_bounds[1], sqrt(r[ok])))
  }
  for(s in seq_len(n_subjects)){
    x <- pm_settings[[s]]
    for(i in seq_along(x)){
      epsilon <- x[[i]]$epsilon
      layout <- x[[i]]$step_layout
      if (length(layout) == length(epsilon)) {
        use_scale <- layout == "chains"
        epsilon[use_scale] <- epsilon[use_scale] * scale[s]
      }
      epsilon[!is.finite(epsilon) | epsilon <= 0] <- 1
      x[[i]]$epsilon <- pmin(epsilon_bounds[2],
                             pmax(epsilon_bounds[1], epsilon))
    }
    pm_settings[[s]] <- x
  }
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  return(pmwgs)
}

test_adapted <- function(sampler, test_samples, min_unique, n_cores_conditional = 1,
                         verbose = FALSE)
{

  first_par <- as.matrix(test_samples$alpha[1, , ])
  if(ncol(first_par) == 1) first_par <- t(first_par)
  n_unique_sub <- apply(first_par, 1, FUN = function(x) return(length(unique(x))))
  n_pars <- sampler$n_pars
  components <- attr(sampler$data, "components")
  nuisance <- sampler$nuisance
  if (length(n_unique_sub) != 0 & all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      for(comp in unique(components)){
        idx <- comp == components
        nuis_idx <- nuisance[idx]
        if(any(nuis_idx)){
          type <- sampler$sampler_nuis$type
          auto_mclapply(X = 1:sampler$n_subjects,
                        FUN = get_conditionals, samples = test_samples$nuisance,
                        n_pars = sum(idx[nuisance]), idx = idx[nuisance],
                        type = type,
                        mc.cores = n_cores_conditional)
        }
        auto_mclapply(X = 1:sampler$n_subjects,FUN = get_conditionals,samples = test_samples,
                      n_pars = sum(idx[!nuisance]), idx = idx[!nuisance], type = sampler$type,
                      mc.cores = n_cores_conditional)
      }
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      reason <- conditionMessage(attempt)
      if(verbose){
        message("Efficient proposals could not be built yet (", reason,
                "); continuing adaptation")
      }
      return(structure(FALSE, reason = reason))
    }
    else {
      if(verbose) message("Successfully adapted - stopping adaptation")
      return(TRUE)
    }
  } else{
    return(structure(FALSE, reason = "min_unique"))
  }
}

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#' Make an emc Object
#'
#' Creates an emc object by combining the data, prior,
#' and model specification into a `emc` object that is needed in `fit()`.
#'
#' @param data A data frame, or a list of data frames. Needs to have the variable `subjects` as participant identifier.
#' @param design A list with a pre-specified design, the output of `design()`.
#' @param model A model list. If none is supplied, the model specified in `design()` is used.
#' @param type A string indicating whether to run a `standard` group-level, `blocked`, `diagonal`, `factor`, or `single` (i.e., non-hierarchical) model.
#' @param n_chains An integer. Specifies the number of mcmc chains to be run (has to be more than 1 to compute `rhat`).
#' @param compress A Boolean, if `TRUE` (i.e., the default), the data is compressed to speed up likelihood calculations.
#' @param rt_resolution A double. Used for compression, response times will be binned based on this resolution.
#'   Each response time then stands for its bin: race models average the
#'   density over the part of the bin above `t0`, other models evaluate it at
#'   the bin midpoint. `NULL` uses the response times as recorded.
#' @param use_data Boolean. If `TRUE` (the default), parameters for factor
#'   combinations that are not present in `data` are omitted from the model.
#'   If `FALSE`, the complete factorial design is retained for backwards
#'   compatibility.
#' @param group_design A design for group-level mappings, made using `group_design()`.
#' @param par_groups A vector. Indicates which parameters are allowed to correlate. Could either be a list of character vectors of covariance blocks. Or
#' a numeric vector, e.g., `c(1,1,1,2,2)` means the covariances
#' of the first three and of the last two parameters are estimated as two separate blocks.
#' @param prior_list A named list containing the prior. Default prior created if `NULL`. For the default priors, see `?get_prior_{type}`.
#' @param ... Additional, optional arguments.
#' @details
#' Some likelihood backends attach internal, versioned cache attributes to the
#' processed per-subject data (dadm) to avoid repeating data-only scans during
#' sampling. These attributes are recomputed when missing or outdated and are not
#' part of the user-facing API.
#' @return An uninitialized emc object
#' @examples dat <- forstmann
#'
#' # function that takes the lR factor (named diff in the following function) and
#' # returns a logical defining the correct response for each stimulus. In this
#' # case the match is simply such that the S factor equals the latent response factor.
#' matchfun <- function(d)d$S==d$lR
#'
#' # design an "average and difference" contrast matrix
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
#'
#' # specify design
#' design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
#' formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#' contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#'
#' # specify priors
#' pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
#'            B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
#' psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
#'          B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
#' prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)
#'
#' # create emc object
#' LBABE <- make_emc(dat,design_LBABE,type="standard",  prior=prior_LBABE,
#'                   compress = FALSE)
#'
#' @export

make_emc <- function(data,design,model=NULL,
                    type="standard",
                    n_chains=3,compress=TRUE,rt_resolution=1/60,
                    use_data = TRUE,
                    prior_list = NULL, group_design = NULL,
                    par_groups=NULL, ...){
  type_was_missing <- missing(type)
  design <- .restore_custom_kernel_design(design, quiet = TRUE)
  n_factors <- NULL
  nuisance <- NULL
  nuisance_non_hyper <- NULL
  sem_settings <- NULL
  Lambda_mat <- NULL
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(!is.null(prior_list) & !is.null(prior_list$theta_mu_mean)){
    prior_list <- list(prior_list)
  }
  if(!is.null(prior_list) && type_was_missing){
    type <- attr(prior_list[[1]], "sampler_type")
    if (is.null(type)) type <- attr(prior_list[[1]], "type")
  }
  if(!is.null(group_design) && !type %in% c("standard", "diagonal", "blocked")){
    stop("group_design can only be used with standard, blocked or diagonal type")
  }
  if (length(use_data) != 1 || !is.logical(use_data) || is.na(use_data))
    stop("use_data must be a single TRUE/FALSE value")
  if(type != "single" && length(unique(data$subjects)) == 1){
    stop("can only use type = `single` when there's only one subject in the data")
  }

  if (!(type %in% c("standard","diagonal","blocked","factor","single", "infnt_factor", "SEM", "diagonal-gamma")))
    stop("type must be one of: standard,diagonal,blocked,factor,infnt_factor, single")

  if(!is.null(nuisance) & !is.null(nuisance_non_hyper)){
    stop("You can only specify nuisance OR nuisance_non_hyper")
  }
  if (is(data, "data.frame")) data <- list(data)
  data <- lapply(data,function(d){
    d$subjects <- factor(d$subjects)
    d <- d[order(d$subjects),]
    d <- add_trials(d)
    d
  })
  for (i in 1:length(data)) {
    if (inherits(design, "emc.design"))
      des <- design else des <- design[[i]]
    if (!any(names(data[[i]])=="LT") & !is.null(des$TC$LT))
      data[[i]]$LT <- des$TC$LT
    if (!any(names(data[[i]])=="LC") & !is.null(des$TC$LC))
      data[[i]]$LC <- des$TC$LC
    if (!any(names(data[[i]])=="UT") & !is.null(des$TC$UT))
      data[[i]]$UT <- des$TC$UT
    if (!any(names(data[[i]])=="UC") & !is.null(des$TC$UC))
      data[[i]]$UC <- des$TC$UC
  }
  if (!is.null(names(design)[1]) && names(design)[1]=="Flist"){
    design <- list(design)
  }
  checks <- sapply(design, function(x) is(x, "emc.design"))
  if (!all(checks)) stop("design must be a list of emc.design objects")
  if (length(design)!=length(data)){
    design <- rep(design,length(data))
  }
  for (i in seq_along(design)) {
    if (isTRUE(use_data)) {
      attr(design[[i]], "data") <- data[[i]]
    } else {
      attr(design[[i]], "data") <- NULL
    }
  }
  if (is.null(model)) model <- lapply(design,function(x){x$model})
  if (any(unlist(lapply(model,is.null))))
    stop("Must supply model if model is not in all design components")
  if (!is.null(names(model)[1]) && names(model)[1]=="type")
    model <- list(model)
  if (length(model)!=length(data))
    model <- rep(model,length(data))

  compress_passed <- compress
  compress <- rep(compress, length(model))
  has_delta_rule <- sapply(model, has_delta_rules)
  compress[has_delta_rule] <- FALSE
  if(compress_passed & any(has_delta_rule)) {
    if(length(model) == 1) message('Because the model contains a delta rule, data will not be compressed.')
    else message(paste0('Models ', which(has_delta_rule), ' contain a delta rule; the corresponding data will not be compressed.'))
  }

  no_bin <- !sapply(model, model_compress_ok)
  compress[no_bin] <- FALSE

  dadm_list <- vector(mode="list",length=length(data))
  dadm_names <- names(design)
  if (is.null(dadm_names)) dadm_names <- as.character(seq_along(dadm_list))
  names(dadm_list) <- dadm_names
  if (is.null(rt_resolution)) {
    rt_resolution <- rep(list(NULL), length(data))
  } else {
    rt_resolution <- as.list(rep(rt_resolution, length.out = length(data)))
  }
  if (any(no_bin)) {
    was_binned <- no_bin & !vapply(rt_resolution, is.null, logical(1))
    rt_resolution[no_bin] <- list(NULL)
    if (compress_passed || any(was_binned)) {
      what <- if (length(model) == 1) "The model's likelihood is" else
        paste0("The likelihood of model(s) ", paste(which(no_bin), collapse=", "), " is")
      message(what, " a cached grid solve, so binning response times saves no ",
              "computation while biasing the fit: using compress = FALSE and ",
              "rt_resolution = NULL.")
    }
  }
  for (i in 1:length(dadm_list)) {
    message("Processing data set ",i)
    if(is.null(attr(design[[i]], "custom_ll"))){
      dadm_list[[i]] <- design_model(data=data[[i]],design=design[[i]],
                                     compress=compress[[i]],model=model[[i]],rt_resolution=rt_resolution[[i]],
                                     drop_unobserved = isTRUE(use_data))
    } else{
      dadm_list[[i]] <- design_model_custom_ll(data = data[[i]],
                                               design = design[[i]],model=model[[i]])
    }
    dadm_list[[i]] <- .cache_ll_data_attrs(dadm_list[[i]], force_rebuild = TRUE)
  }
  class(design) <- "emc.design"
  prior_in <- merge_priors(prior_list)

  prior_in <- prior(design, type, update = prior_in, group_design = group_design, ...)
  attr(dadm_list[[1]], "prior") <- prior_in

  marginalise_flag <- attr(design[[1]], "marginalise")
  if (!is.null(marginalise_flag) && !type %in% c("single", "standard", "blocked", "diagonal")) {
    stop("marginalise is currently only supported for single, standard, blocked and diagonal types")
  }

  if (type %in% c("single", "infnt_factor", "diagonal-gamma")) {
    out <- pmwgs(dadm_list, type, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 n_factors = n_factors,
                 marginalise = marginalise_flag)
  } else if (type %in% c("standard", "blocked", "diagonal")) {
    if(type == "blocked"){
      if(is.null(par_groups)) stop("par_groups must be specified for blocked models")
    }
    if(type == "diagonal" && is.null(par_groups)){
      par_groups <- 1:length(sampled_pars(design))
    }
    if(type == "standard" && is.null(par_groups)){
      par_groups <- rep(1, length(sampled_pars(design)))
    }
    if(!is.null(par_groups)){
      if(is.character(par_groups)){
        par_groups <- list(par_groups)
      }
      if(is.list(par_groups)){
        par_names <- names(sampled_pars(design))
        new_par_groups <- rep(NA, length(par_names))
        for(i in 1:length(par_groups)){
          if(any(!par_groups[[i]] %in% par_names)) stop("Make sure you specified parameter names in par_groups correctly")
          new_par_groups[par_names %in% par_groups[[i]]] <- i
        }
        new_par_groups[is.na(new_par_groups)] <- (i+1):(i+sum(is.na(new_par_groups)))
        par_groups <- new_par_groups
      }
      if(length(par_groups) != length(sampled_pars(design))){
        stop("par_groups length does not match number of sampled parameters, make sure you specified par_groups correctly")
      }
      # Downstream code treats labels as 1..k; c(2, 2, 5, 5) must mean two blocks.
      par_groups <- match(par_groups, unique(par_groups))
    }
    if(type %in% c("diagonal", "blocked")) type <- "standard"
    out <- pmwgs(dadm_list, type, par_groups=par_groups,
                 group_design = group_design,
                 nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 marginalise = marginalise_flag)
  } else if (type == "factor") {
    out <- pmwgs(dadm_list, type, n_factors = n_factors,
                 nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 Lambda_mat = Lambda_mat)
  } else if (type == "SEM"){
    out <- pmwgs(dadm_list, type, sem_settings = sem_settings, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper)
  }
  out$model <- lapply(design, function(x) x$model)
  if(length(out$model) == 1) out$model <- out$model[[1]]
  out <- check_duplicate_designs(out)
  dadm_lists <- rep(list(out),n_chains)
  class(dadm_lists) <- "emc"
  return(dadm_lists)
}

fix_fileName <- function(x){
  ext <- substr(x, nchar(x)-5, nchar(x))
  if(ext != ".RData" & ext != ".Rdata"){
    return(paste0(x, ".RData"))
  } else{
    return(x)
  }
}

check_duplicate_designs <- function(out){
  if(is.data.frame(out$data[[1]])) return(out)
  for(i in 1:length(out$data)){
    designs <- lapply(out$data[[i]], function(y) attr(y, "designs"))
    duplicacy <- duplicated(designs)
    unq_idx <-   sapply(seq_along(designs), function(i) {
      for (j in seq_along(designs)) {
        if (identical(designs[[i]], designs[[j]])) return(j)
      }
    })
    for(j in 1:length(out$data[[i]])){
      if(is.null(designs[[j]])) next
      if(duplicacy[j]){
        attr(out$data[[i]][[j]], "designs") <- unq_idx[j]
      }
    }
  }
  return(out)
}

extractDadms <- function(dadms, names = NULL, par_names = NULL){
  if(is.null(names)) {
    names <- base::names(dadms)
    if (is.null(names)) names <- 1:length(dadms)
  }
  N_models <- length(dadms)
  sampled_names <- attr(dadms[[1]], "sampled_p_names")
  prior <- attr(dadms[[1]], "prior")
  subjects <- unique(factor(sapply(dadms, FUN = function(x) levels(x$subjects))))
  dadm_list <- dm_list(dadms[[1]])
  if (N_models == 1) {
    components <- setNames(rep.int(1L, length(sampled_names)), sampled_names)
  } else {
    total_dadm_list <- vector("list", length = N_models)
    k <- 1
    model_pars <- paste(names[1], sampled_names, sep = "|")
    dadm_list[as.character(which(!subjects %in% unique(dadms[[1]]$subjects)))] <- NA
    total_dadm_list[[1]] <- dadm_list
    components <- setNames(rep.int(1L, length(model_pars)), model_pars)
    for(dadm in dadms[-1]){
      k <- k + 1
      tmp_list <- vector("list", length = length(subjects))
      tmp_list[as.numeric(unique(dadm$subjects))] <- dm_list(dadm)
      total_dadm_list[[k]] <- tmp_list
      curr_pars <- attr(dadm, "sampled_p_names")
      model_pars <- paste(names[k], curr_pars, sep = "|")
      components <- c(components,
                      setNames(rep.int(as.integer(k), length(model_pars)), model_pars))
    }
    dadm_list <- do.call(mapply, c(list, total_dadm_list, SIMPLIFY = F))
  }
  if (!is.null(par_names)) {
    source_names <- base::names(components)
    if (is.null(source_names)) {
      stop("Parameter bookkeeping error: component metadata has no parameter names")
    }
    aligned <- components[match(par_names, source_names)]
    missing <- is.na(aligned)
    if (any(missing)) {
      model_idx <- if (N_models == 1L) {
        rep.int(1L, sum(missing))
      } else {
        match(sub("\\|.*$", "", par_names[missing]), names)
      }
      if (anyNA(model_idx)) {
        stop("Parameter bookkeeping error: cannot assign components to sampled parameters")
      }
      aligned[missing] <- model_idx
    }
    components <- setNames(as.integer(aligned), par_names)
  }
  attr(dadm_list, "components") <- components
  attr(dadm_list, "shared_ll_idx") <- components
  return(list(prior = prior,
              dadm_list = dadm_list, subjects = subjects))
}

.emc_core_ctl <- function(n_chains, cores_per_chain,
                          cores_for_chains = n_chains) {
  n_chains <- as.integer(n_chains)
  cores_per_chain <- as.integer(cores_per_chain)
  cores_for_chains <- max(1L, as.integer(cores_for_chains))
  if (n_chains <= 1 || cores_per_chain < 1 || cores_for_chains <= 1) return(NULL)
  if (Sys.info()[1] == "Windows") return(NULL)
  dir <- tempfile(paste0("emc_cores_", Sys.getpid(), "_"))
  if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) return(NULL)
  list(dir = dir, n_chains = n_chains,
       n_slots = cores_for_chains,
       total = cores_for_chains * cores_per_chain)
}

.emc_cores_now <- function(core_ctl, base) {
  if (is.null(core_ctl)) return(base)
  unfinished <- core_ctl$n_chains - length(list.files(core_ctl$dir))
  active <- min(core_ctl$n_slots, unfinished)
  if (active <= 0L) return(core_ctl$total)
  max(base, core_ctl$total %/% active)
}

.emc_core_release <- function(core_ctl) {
  if (!is.null(core_ctl)) {
    file.create(file.path(core_ctl$dir, paste0("done_", Sys.getpid())),
                showWarnings = FALSE)
  }
  invisible(NULL)
}

# A forked child must die with its caller: an interrupted mclapply SIGTERMs
# it, and its own workers follow through their death signal.  `fun` arms only
# in a fork, never in the process that made it, which mclapply() also runs
# `fun` in when it has one core.
.emc_armed <- function(fun) {
  owner <- Sys.getpid()
  function(...) {
    if (Sys.getpid() != owner &&
        identical(tryCatch(emc_arm_parent_death(), error = function(e) 0L), -1L))
      stop("the calling process exited before this child started")
    fun(...)
  }
}

# A forked child's first full garbage collection sweeps every page of its
# parent's small-object heap, dead objects included, and so takes a private
# copy of all of it.  Collecting in the parent first releases the pages that
# hold only garbage, about halving what each child copies once a long fit has
# churned its heap.  Throttled, since chains fork too.
.emc_fork_gc <- new.env(parent = emptyenv())
.emc_gc_before_fork <- function(n_forks) {
  if (n_forks < 3L) return(invisible(NULL))
  now <- proc.time()[["elapsed"]]
  last <- if (identical(.emc_fork_gc$pid, Sys.getpid())) .emc_fork_gc$at else -Inf
  if (now - last < 10) return(invisible(NULL))
  invisible(gc(full = TRUE))
  .emc_fork_gc$pid <- Sys.getpid()
  .emc_fork_gc$at <- proc.time()[["elapsed"]]
  invisible(NULL)
}

# Namespace functions are lazy-loaded, and a forked child that loads one keeps
# it only until it exits. Each child therefore reports which still-unloaded
# bindings it loaded, and the parent loads them, so that later forks inherit
# them instead of loading them again.
.emc_ns_unloaded <- function() {
  ns <- asNamespace("EMC2")
  nms <- ls(ns, all.names = TRUE, sorted = FALSE)
  nms[emc_ns_unforced(ns, nms)]
}

.emc_ns_reporting <- function(fun, unloaded) {
  force(fun)
  function(...) {
    value <- fun(...)
    ns <- asNamespace("EMC2")
    loaded <- unloaded[!emc_ns_unforced(ns, unloaded)]
    structure(list(value = value, loaded = loaded),
              class = "emc_ns_report")
  }
}

.emc_ns_load <- function(reports) {
  loaded <- unique(unlist(lapply(reports, function(r)
    if (inherits(r, "emc_ns_report")) r$loaded), use.names = FALSE))
  if (!length(loaded)) return(invisible(NULL))
  ns <- asNamespace("EMC2")
  loaded <- loaded[emc_ns_unforced(ns, loaded)]
  for (nm in loaded) invisible(get(nm, envir = ns, inherits = FALSE))
  invisible(NULL)
}

auto_mclapply <- function(X, FUN, mc.cores, ..., mc.preschedule = TRUE,
                          mc.set.seed = TRUE){
  if(mc.cores <= 1) return(lapply(X, FUN, ...))
  if(Sys.info()[1] == "Windows"){
    cluster <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    list_out <- parallel::parLapply(cl = cluster, X,FUN, ...)
  } else{
    .emc_gc_before_fork(mc.cores)
    unloaded <- .emc_ns_unloaded()
    run <- if (length(unloaded)) .emc_ns_reporting(FUN, unloaded) else FUN
    list_out <- parallel::mclapply(X, .emc_armed(run), mc.cores = mc.cores,
                                   mc.preschedule = mc.preschedule,
                                   mc.set.seed = mc.set.seed, ...)
    if (length(unloaded)) {
      .emc_ns_load(list_out)
      list_out <- lapply(list_out, function(r)
        if (inherits(r, "emc_ns_report")) r$value else r)
    }
    failed <- vapply(list_out, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      first <- list_out[[which(failed)[1]]]
      cond <- attr(first, "condition")
      if (is.null(cond)) stop(as.character(first), call. = FALSE) else stop(cond)
    }
  }
  return(list_out)
}

.EMC_SHARED_CHAIN_FIELDS <- c("data", "model", "prior", "gd", "XtX", "group_designs",
                             "par_names", "subjects", "n_pars", "n_subjects", "type",
                             "nuisance", "marginalise", "marginalised_idx",
                             "par_group", "is_blocked")

#' Share identical immutable configuration across chains
#' @param emc A list of EMC objects
#' @return The same list with identical immutable fields shared with the first chain
#' @noRd
strip_duplicates <- function(emc, incl_props = TRUE) {
  emc <- restore_duplicates(emc)
  # Only immutable configuration may be shared. All other fields (including
  # future dynamic fields) remain chain-local, even if currently identical.
  for (i in seq_along(emc)[-1]) {
    shared <- intersect(.EMC_SHARED_CHAIN_FIELDS, names(emc[[i]]))
    shared <- shared[vapply(shared, function(k) {
      identical(emc[[i]][[k]], emc[[1]][[k]])
    }, logical(1))]
    attr(emc[[i]], "emc_field_order") <- names(emc[[i]])
    emc[[i]][shared] <- NULL
    attr(emc[[i]], "emc_shared_fields") <- shared
  }
  # incl_props is retained for call compatibility. Proposal caches are dynamic
  # sampler state: deleting them makes autosave change the next transition.
  return(emc)
}

#' Restore shared immutable configuration without replacing chain-local state
#' @param emc A list of EMC objects with stripped duplicates
#' @param warn Whether to warn that an object without chain-local adaptive
#'   state cannot continue exactly (only relevant when sampling resumes)
#' @return The same list with shared fields restored and dynamic state preserved
#' @noRd
restore_duplicates <- function(emc, warn = FALSE) {
  if (length(emc) > 1) {
    for (i in 2:length(emc)) {
      shared <- attr(emc[[i]], "emc_shared_fields", exact = TRUE)
      if (!is.null(shared)) {
        for (k in shared) emc[[i]][k] <- emc[[1]][k]
        order <- attr(emc[[i]], "emc_field_order", exact = TRUE)
        if (!is.null(order)) {
          attrs <- attributes(emc[[i]])
          emc[[i]] <- emc[[i]][order]
          attrs$names <- order
          attributes(emc[[i]]) <- attrs
        }
        attr(emc[[i]], "emc_shared_fields") <- NULL
        attr(emc[[i]], "emc_field_order") <- NULL
      } else if (is.null(emc[[i]]$n_pars) && !is.null(emc[[i]]$samples)) {
        missing <- setdiff(.EMC_SHARED_CHAIN_FIELDS, names(emc[[i]]))
        for (k in missing) emc[[i]][k] <- emc[[1]][k]
        idx <- emc[[i]]$samples$idx
        if (is.null(emc[[i]]$init))
          emc[[i]]$init <- !is.null(idx) && idx >= 1L &&
            all(is.finite(emc[[i]]$samples$alpha[, , idx, drop = FALSE]))
        if (warn)
          warning("This object lacks chain-local adaptive state; exact continuation is unavailable.",
                  call. = FALSE)
      }
    }
  }
  return(emc)
}

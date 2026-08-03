resolve_marginalise_prior <- function(marginalise, prior) {
  if (is.null(marginalise)) return(NULL)
  param <- if (is.list(marginalise)) marginalise$param else marginalise
  if (!is.null(names(marginalise)) && "param" %in% names(marginalise)) {
    param <- marginalise[["param"]]
  } else if (is.character(marginalise) || is.atomic(marginalise)) {
    nm <- names(marginalise)
    if (!is.null(nm)) {
      idx_param <- which(!nm %in% c("n_nodes", "nodes"))
      if (length(idx_param) > 0) param <- marginalise[idx_param[1]] else param <- marginalise[1]
    } else {
      param <- marginalise[1]
    }
  }
  param <- as.character(param[1L])
  if (length(param) != 1L || is.na(param) || !nzchar(param)) {
    stop("marginalise must identify exactly one parameter")
  }
  mu <- prior$theta_mu_mean[[param]]
  if (is.null(mu) || !is.finite(mu)) {
    stop("Could not find a finite prior mean for marginalized parameter '", param, "'")
  }
  vv <- prior$theta_mu_var
  if (is.matrix(vv)) {
    vi <- if (!is.null(dimnames(vv)) && !is.null(rownames(vv))) {
      match(param, rownames(vv))
    } else match(param, names(prior$theta_mu_mean))
    if (is.na(vi)) stop("Could not find prior variance for marginalized parameter '", param, "'")
    variance <- vv[vi, vi]
  } else {
    variance <- if (!is.null(names(vv))) vv[[param]] else vv[[match(param, names(prior$theta_mu_mean))]]
  }
  sigma <- sqrt(variance)
  if (!is.finite(sigma) || sigma <= 0) {
    stop("Marginalized parameter '", param, "' needs a finite, positive prior SD")
  }
  n_nodes <- if (is.list(marginalise) && !is.null(marginalise$n_nodes)) {
    as.integer(marginalise$n_nodes)
  } else if (!is.null(names(marginalise)) && "n_nodes" %in% names(marginalise)) {
    as.integer(marginalise[["n_nodes"]])
  } else if (length(marginalise) > 1 && is.null(names(marginalise))) {
    if (suppressWarnings(!is.na(as.integer(marginalise[2])))) as.integer(marginalise[2]) else 12L
  } else 12L

  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes < 2L) {
    stop("marginalise$n_nodes must be an integer >= 2")
  }
  list(param = param, mu = unname(mu), sigma = unname(sigma), n_nodes = n_nodes)
}

# Build the t0 quadrature grid for one or more proposals: the node values
# (sampled t0 scale) and the np x K unnormalized log-terms. Shared by the
# particle step (which reduces it to the marginal ll AND reuses the accepted
# row for reconstruction) and the init/start path, so the grid is only ever
# computed once per likelihood evaluation.
compute_marginal_grid <- function(proposals, data, model, marginalise,
                                  r_cores = 1, warm = NULL) {
  model_spec <- if (is.function(model)) model() else model
  set_stop_method_from_model(model_spec)  # SS models: stop_method -> C++ config
  data <- .cache_ll_data_attrs(data)
  constants <- attr(data, "constants")
  if (is.null(constants)) constants <- NA
  designs <- .oo_expanded_designs(data)
  # Warm start: last iteration's accepted mode/scale for this subject. It only
  # seeds the probe, and the C++ side falls back to the full pilot scan if any
  # particle comes out unresolved, so a stale hint costs time, never accuracy.
  if (!is.null(warm) && all(is.finite(warm)) && warm[["sd"]] > 0) {
    marginalise$warm_mode <- unname(warm[["mode"]])
    marginalise$warm_sd <- unname(warm[["sd"]])
  }
  one <- function(props) {
    calc_ll_oo_marginal_nodes(
      props, data, constants = constants, designs = designs,
      type = model_spec$c_name, bounds = model_spec$bound,
      transforms = model_spec$transform, pretransforms = model_spec$pre_transform,
      p_types = names(model_spec$p_types), min_ll = log(1e-10),
      marginalise = marginalise, trend = model_spec$trend
    )
  }
  # Each particle now carries its own quadrature rule, so particles split across
  # cores exactly like the ordinary likelihood path in calc_ll_manager.
  if (r_cores <= 1 || nrow(proposals) <= r_cores) return(one(proposals))
  idx <- rep(1:r_cores, each = 1 + (nrow(proposals) %/% r_cores))[1:nrow(proposals)]
  parts <- auto_mclapply(1:r_cores, function(i) {
    one(proposals[idx == i, , drop = FALSE])
  }, mc.cores = r_cores)
  list(nodes = do.call(rbind, lapply(parts, `[[`, "nodes")),
       log_terms = do.call(rbind, lapply(parts, `[[`, "log_terms")),
       ll = unlist(lapply(parts, `[[`, "ll")),
       mode = unlist(lapply(parts, `[[`, "mode")),
       sd = unlist(lapply(parts, `[[`, "sd")),
       # fraction of splits that kept the hint (the backoff gates on a majority)
       warm_used = mean(unlist(lapply(parts, `[[`, "warm_used"))),
       pred_used = mean(unlist(lapply(parts, `[[`, "pred_used"))),
       repaired = sum(unlist(lapply(parts, `[[`, "repaired"))))
}

# Warm-start state carried on pm_settings between iterations: the Laplace fit
# of the accepted particle, or NULL when that particle had no usable fit.
marginal_warm_state <- function(grid, idx) {
  if (is.null(grid$mode) || is.null(grid$sd)) return(NULL)
  m <- grid$mode[idx]; s <- grid$sd[idx]
  if (!is.finite(m) || !is.finite(s) || s <= 0) return(NULL)
  c(mode = m, sd = s)
}

# A single hint can only serve the whole particle batch while the batch is
# tight: the conditional t0 mode moves with the other parameters, so a wide
# proposal cloud (preburn/burn) spreads the per-particle modes over hundreds of
# posterior SDs and the hint is rejected.  A rejected hint costs one wasted
# probe round, so back off geometrically after a failure and retry later -- the
# cloud contracts as the chain converges, and the hint then holds.
marginal_warm_backoff <- function(state, attempted, used) {
  pen <- if (is.null(state$penalty)) 0L else state$penalty
  wait <- if (is.null(state$backoff)) 0L else state$backoff
  if (attempted) {
    if (isTRUE(used)) { pen <- 0L; wait <- 0L }
    else { pen <- min(16L, max(1L, pen * 2L)); wait <- pen }
  } else wait <- max(0L, wait - 1L)
  list(penalty = pen, backoff = wait)
}

# Reduce a grid of log-terms (np x K) to the marginal log-likelihood per
# proposal by a numerically stable row-wise log-sum-exp.
marginal_ll_from_grid <- function(log_terms) {
  if (is.list(log_terms) && !is.null(log_terms$ll)) return(log_terms$ll)
  m <- apply(log_terms, 1L, max, na.rm = TRUE)
  m[!is.finite(m)] <- -Inf
  fin <- is.finite(m)
  res <- rep(-Inf, nrow(log_terms))
  if (any(fin)) {
    res[fin] <- m[fin] + log(rowSums(exp(log_terms[fin, , drop = FALSE] - m[fin]), na.rm = TRUE))
  }
  res
}

# Draw one t0 node ~ Categorical(softmax(log_terms_row)) for the stored alpha.
# `nodes` is that particle's OWN node row (rules are per particle).
# Falls back to the fixed prior when the row carries no finite mass.
draw_marginal_node <- function(nodes, log_terms_row, marginalise) {
  nodes <- as.numeric(nodes)
  if (length(nodes) == 0L || all(!is.finite(log_terms_row))) {
    return(stats::rnorm(1L, marginalise$mu, marginalise$sigma))
  }
  w <- exp(log_terms_row - max(log_terms_row[is.finite(log_terms_row)]))
  w[!is.finite(w)] <- 0
  if (!any(w > 0)) return(stats::rnorm(1L, marginalise$mu, marginalise$sigma))
  nodes[[sample.int(length(nodes), size = 1L, prob = w)]]
}

# Standalone reconstruction (init/no-reuse path): computes its own grid.
reconstruct_marginalised_particle <- function(particle, data, model, marginalise) {
  if (is.null(marginalise)) return(as.numeric(particle))
  p_names <- names(particle)
  if (is.null(p_names)) stop("A marginalized particle must have parameter names")
  particle <- as.numeric(particle)
  names(particle) <- p_names
  p_idx <- match(marginalise$param, p_names)
  if (is.na(p_idx)) stop("Marginalized parameter is missing from the particle")
  # A subject with no observations has no posterior node weights. Draw from the
  # fixed prior so the stored alpha still has the same shape as ordinary fits.
  if (is.null(data) || nrow(data) == 0L) {
    particle[p_idx] <- stats::rnorm(1L, marginalise$mu, marginalise$sigma)
    return(particle)
  }
  grid <- compute_marginal_grid(
    matrix(particle, nrow = 1L, dimnames = list(NULL, p_names)), data, model, marginalise)
  particle[p_idx] <- draw_marginal_node(grid$nodes[1L, ],
                                        as.numeric(grid$log_terms[1L, ]), marginalise)
  particle
}

pmwgs <- function(dadm, type, pars = NULL, prior = NULL,
                  nuisance = NULL, nuisance_non_hyper = NULL, marginalise = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  if(is.null(pars)) pars <- names(sampled_pars(attr(dadm[[1]], "prior")))
  if(is.null(prior)) prior <- attr(dadm[[1]], "prior")
  # The opt-in flag arrives as an explicit argument (from design(marginalise=)),
  # NOT tagged onto the dadm: downstream consumers (predict / make_data / IC)
  # read the reconstructed t0 like any other parameter and must never integrate.
  marginalise <- resolve_marginalise_prior(marginalise, prior)
  dadm <- extractDadms(dadm)

  dadm_list <-dadm$dadm_list
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  if(!is.null(nuisance) & !is.numeric(nuisance)) nuisance <- which(pars %in% nuisance)
  if(!is.null(nuisance_non_hyper) & !is.numeric(nuisance_non_hyper)) nuisance_non_hyper <- which(pars %in% nuisance_non_hyper)

  if (!is.null(marginalise) && any(pars %in% marginalise$param &
                                   is.element(seq_along(pars),
                                              unique(c(nuisance, nuisance_non_hyper))))) {
    stop("A marginalized parameter cannot be assigned to a nuisance sampler")
  }

  if(!is.null(nuisance_non_hyper)){
    is_nuisance <- is.element(seq_len(length(pars)), nuisance_non_hyper)
    nuis_type <- "single"
  } else if(!is.null(nuisance)) {
    is_nuisance <- is.element(seq_len(length(pars)), nuisance)
    nuis_type <- "diagonal"
  } else{
    is_nuisance <- rep(F, length(pars))
  }


  sampler_nuis <- NULL
  if(any(is_nuisance)){
    sampler_nuis <- list(
      samples = sample_store(dadm, pars, nuis_type, integrate = F,
                                                    is_nuisance = !is_nuisance, ...),
      n_subjects = length(subjects),
      n_pars = sum(is_nuisance),
      nuisance = rep(F, sum(is_nuisance)),
      type = nuis_type
    )
    if(nuis_type == "single") sampler_nuis$samples <- NULL
    sampler_nuis <- add_info(sampler_nuis, prior$prior_nuis, nuis_type, ...)
  }
  samples <- sample_store(dadm, pars, type, is_nuisance = is_nuisance, ...)
  sampler <- list(
    data = dadm_list,
    par_names = pars,
    subjects = subjects,
    n_pars = length(pars),
    nuisance = is_nuisance,
    n_subjects = length(subjects),
    samples = samples,
    sampler_nuis = sampler_nuis,
    type = type,
    marginalise = marginalise,
    marginalised_idx = pars %in% marginalise$param,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- add_info(sampler, prior, type, ...)
  return(sampler)
}

.particle_core_budget <- function(n_subjects, n_cores = 1L, r_cores = 1L) {
  n_subjects <- as.integer(n_subjects)
  n_cores <- max(1L, as.integer(n_cores))
  r_cores <- max(1L, as.integer(r_cores))
  if (n_subjects == 1L && n_cores > 1L) {
    return(list(subject = 1L, likelihood = max(r_cores, n_cores)))
  }
  list(subject = n_cores, likelihood = r_cores)
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000,
                 n_cores = 1, r_cores = 1) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  type <- pmwgs$type
  startpoints <-startpoints_comb <- get_startpoints(pmwgs, start_mu, start_var, type)
  if(any(pmwgs$nuisance)){
    type_nuis <- pmwgs$sampler_nuis$type
    startpoints_nuis <- get_startpoints(pmwgs$sampler_nuis, start_mu = NULL, start_var = NULL, type = type_nuis)
    startpoints_comb <- merge_group_level(startpoints$tmu, startpoints_nuis$tmu,
                                          startpoints$tvar, startpoints_nuis$tvar,
                                          pmwgs$nuisance, startpoints$subj_mu)
    pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                      group_level = startpoints_nuis,
                                                                      j = 1,
                                                                      proposals = NULL,
                                                                      n_pars = pmwgs$n_pars, type = type_nuis)
    pmwgs$sampler_nuis$samples$idx <- 1
  }
  # With one subject there is nothing to parallelise in the outer mclapply.
  # Spend the same per-chain core budget across that subject's proposal
  # likelihoods instead.  This matters especially for PDE-backed models, for
  # which a particle is an independent numerical solve.  The total process
  # count still respects n_cores; r_cores remains an explicit lower bound.
  core_budget <- .particle_core_budget(
    pmwgs$n_subjects, n_cores = n_cores, r_cores = r_cores
  )
  proposals <- parallel::mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,
                                  parameters = startpoints_comb, n_particles = particles,
                                  pmwgs = pmwgs, type = type,
                                  mc.cores = core_budget$subject,
                                  r_cores = core_budget$likelihood)
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))

  # Sample the mixture variables' initial values.

  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             j = 1, n_pars = pmwgs$n_pars, type = type)
  pmwgs$init <- TRUE
  return(pmwgs)
}

#' Initialize Chains
#'
#' Adds a set of start points to each chain. These start points are sampled from a user-defined multivariate
#' normal across subjects.
#'
#' @param emc An emc object made by `make_emc()`
#' @param start_mu A vector. Mean of multivariate normal used in proposal distribution
#' @param start_var A matrix. Variance covariance matrix of multivariate normal used in proposal distribution.
#' Smaller values will lead to less deviation around the mean.
#' @param cores_per_chain An integer. How many cores to use per chain.
#' Parallelizes across participant calculations; with one participant,
#' parallelizes its proposal likelihoods instead.
#' @param cores_for_chains An integer. How many cores to use to parallelize across chains. Default is the number of chains.
#' @param particles An integer. Number of starting values
#' @param ... optional additional arguments
#'
#' @return An emc object
#' @examples \donttest{
#' # Make a design and an emc object
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1),
#'                            constants=c(s=log(1)))
#'
#' DDMaE <- make_emc(forstmann, design_DDMaE, compress = FALSE)
#' # set up our mean starting points (same used across subjects).
#' mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'        t0=log(.2))
#' # Small variances to simulate start points from a tight range
#' var <- diag(0.05, length(mu))
#' # Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
#' # so 4*3 cores used.
#' DDMaE <- init_chains(DDMaE, start_mu = mu, start_var = var,
#'                      cores_per_chain = 1, cores_for_chains = 1, particles = 3)
#' # Afterwards we can just use fit
#' # DDMaE <- fit(DDMaE, cores_per_chain = 4)
#' }
#' @export
init_chains <- function(emc, start_mu = NULL, start_var = NULL, particles = 1000,
                        cores_per_chain=1,cores_for_chains = length(emc),
                        ...)
{
  dots <- add_defaults(list(...),r_cores=1)
  emc <- mclapply(emc,init,start_mu = start_mu, start_var = start_var,
           verbose = FALSE, particles = particles,r_cores=dots$r_cores,
           n_cores = cores_per_chain, mc.cores=cores_for_chains)
  class(emc) <- "emc"
  return(emc)
}

start_proposals <- function(s, parameters, n_particles, pmwgs, type, r_cores = 1){
  #Draw the first start point
  group_pars <- get_group_level(parameters, s, type)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha) # preserve par names
  data_s <- pmwgs$data[[which(pmwgs$subjects == s)]]
  marginalise <- pmwgs$marginalise
  if (!is.null(marginalise)) {
    # One grid: reduce to the marginal ll for start-point selection and reuse
    # the chosen particle's node weights to reconstruct its t0.
    grid <- compute_marginal_grid(proposals, data_s, pmwgs$model, marginalise,
                                  r_cores = r_cores)
    lw <- marginal_ll_from_grid(grid)
    weight <- exp(lw - max(lw))
    idx <- sample(x = n_particles, size = 1, prob = weight)
    proposal <- proposals[idx,]
    names(proposal) <- colnames(proposals)
    p_idx <- match(marginalise$param, names(proposal))
    proposal[p_idx] <- draw_marginal_node(grid$nodes[idx, ],
                                          as.numeric(grid$log_terms[idx, ]), marginalise)
    return(list(proposal = proposal, ll = lw[idx]))
  }
  lw <- calc_ll_manager(proposals, dadm = data_s, model = pmwgs$model, r_cores = r_cores)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}


check_tune_settings <- function(tune, n_pars, stage, particles){
  # Acceptance ratio tuning
  tune$alphaStar <- ifelse(stage == "sample", 2, 3)
  tune$p_accept <- set_p_accept(stage, tune$search_width)
  # Potential blocking settings
  if(is.null(tune$components)) tune$components <- rep(1, n_pars)
  if(is.null(tune$shared_ll_idx)) tune$shared_ll_idx <- tune$components
  # Tuning of number of particles, might be a bit arbitrary
  if(is.null(tune$target_ESS)) tune$target_ESS <- 2.5*sqrt(n_pars)
  if(is.null(tune$ESS_scale)) tune$ESS_scale <- .05
  if(is.null(tune$max_particles)) tune$max_particles <- particles*1.2
  # Mix tuning settings
  if(is.null(tune$mix_adapt)) tune$mix_adapt <- .05
  # After n0 all the tuning kicks in
  tune$n0 <- 25
  return(tune)
}

check_sampling_settings <- function(pm_settings, stage, n_pars, particles){
  for(i in 1:length(pm_settings)){
    # Mix settings
    pm_settings[[i]]$mix <- check_mix(pm_settings[[i]]$mix, stage)
    # For p_accept
    pm_settings[[i]]$epsilon <- check_epsilon(pm_settings[[i]]$epsilon, n_pars, pm_settings[[i]]$mix)
    # For mix and p_accept tuning
    pm_settings[[i]]$proposal_counts <- check_prop_performance(pm_settings[[i]]$proposal_counts, stage)
    pm_settings[[i]]$acc_counts <- check_prop_performance(pm_settings[[i]]$acc_counts, stage)
    # Setting particles
    if(is.null(pm_settings[[i]]$n_particles)) pm_settings[[i]]$n_particles <- particles
    if(is.null(pm_settings[[i]]$iter)) pm_settings[[i]]$iter <- 1
    if(is.null(pm_settings[[i]]$gd_good)) pm_settings[[i]]$gd_good <- FALSE
  }
  return(pm_settings)
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      n_cores = 1,
                      tune = NULL,
                      verbose = TRUE,
                      verboseProgress = TRUE,
                      r_cores = 1) {
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  tune$components <- attr(pmwgs$data, "components")
  tune$shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")

  pm_settings <- attr(pmwgs$samples, "pm_settings")
  # Intialize sampling tuning settings
  if(is.null(pm_settings)) pm_settings <- lapply(1:pmwgs$n_subjects, function(x) return(vector("list", length(unique(tune$components)))))
  tune <- check_tune_settings(tune, n_pars, stage, particles)
  pm_settings <- lapply(pm_settings, FUN = check_sampling_settings,  stage = stage, n_pars = n_pars, particles)

  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)

  # Add proposal distributions
  eff_mu <- pmwgs$eff_mu
  eff_var <- pmwgs$eff_var
  chains_var <- pmwgs$chains_var
  chains_mu <- pmwgs$chains_mu
  # Make sure that there's at least something to mapply over
  if(is.null(eff_mu)) eff_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(chains_mu)) chains_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(eff_var)) eff_var <- vector("list", pmwgs$n_subjects)
  if(is.null(chains_var)) chains_var <- vector("list", pmwgs$n_subjects)
  if (verboseProgress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  nuisance <- pmwgs$nuisance
  if(any(nuisance)){
    type <- pmwgs$sampler_nuis$type
    pmwgs$sampler_nuis$samples$idx <- pmwgs$samples$idx
  }
  block_idx <- block_variance_idx(tune$components)
  # Main iteration loop
  for (i in 1:iter) {
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    # Gibbs step. If a numerical failure occurs here, no subject update has
    # happened yet, so the safe rejection is to repeat the previous iteration.
    pars_attempt <- tryCatch(
      gibbs_step(pmwgs, pmwgs$samples$alpha[!nuisance,,j-1], pmwgs$type),
      error = identity
    )
    if (inherits(pars_attempt, c("error", "try-error"))) {
      pmwgs$samples <- reject_sample_iteration(pmwgs$samples, j)
      if(any(nuisance)){
        pmwgs$sampler_nuis$samples <- reject_sample_iteration(pmwgs$sampler_nuis$samples, j)
      }
      next
    }
    pars <- pars_comb <- pars_attempt
    if(any(nuisance)){
      pars_nuis_attempt <- tryCatch(
        gibbs_step(pmwgs$sampler_nuis, pmwgs$samples$alpha[nuisance,,j-1], pmwgs$sampler_nuis$type),
        error = identity
      )
      if (inherits(pars_nuis_attempt, c("error", "try-error"))) {
        pmwgs$samples <- reject_sample_iteration(pmwgs$samples, j)
        pmwgs$sampler_nuis$samples <- reject_sample_iteration(pmwgs$sampler_nuis$samples, j)
        next
      }
      pars_nuis <- pars_nuis_attempt
      pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu, pars$tvar, pars_nuis$tvar, nuisance, pars$subj_mu)
      pars_comb$alpha <- pmwgs$samples$alpha[,,j-1]
      pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                        group_level = pars_nuis,
                                                                        j = j,
                                                                        proposals = NULL,
                                                                        n_pars = n_pars, type = pmwgs$sampler_nuis$type)
      pmwgs$sampler_nuis$samples$idx <- j
    }
    # Particle step
    # A single-subject chain otherwise leaves cores_per_chain - 1 workers idle.
    # Route that existing budget into the independent proposal likelihoods.
    core_budget <- .particle_core_budget(
      pmwgs$n_subjects, n_cores = n_cores, r_cores = r_cores
    )
    proposals <- parallel::mcmapply(safe_new_particle, 1:pmwgs$n_subjects, data, pm_settings, eff_mu, eff_var,
                                    chains_mu, chains_var, pmwgs$samples$subj_ll[,j-1],
                                    MoreArgs = list(parameters = pars_comb,
                                                    model = pmwgs$model,
                                                    stage = stage,
                                                    type = pmwgs$type,
                                                    tune = tune,
                                                    marginalise = pmwgs$marginalise,
                                                    r_cores = core_budget$likelihood),
                                    mc.cores = core_budget$subject)
    pm_settings <- proposals[3,]
    proposals <- array(unlist(proposals[1:2,]), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, j = j, n_pars = pmwgs$n_pars, type = pmwgs$type)
  }
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  if (verboseProgress) close(pb)
  return(pmwgs)
}

reject_sample_iteration <- function(samples, j) {
  if (j <= 1) {
    samples$idx <- j
    return(samples)
  }
  for (nm in names(samples)) {
    obj <- samples[[nm]]
    d <- dim(obj)
    if (is.null(d)) next
    if (length(d) == 2 && d[2] >= j && d[1] != d[2]) {
      samples[[nm]][, j] <- samples[[nm]][, j - 1]
    } else if (length(d) == 3 && d[3] >= j) {
      samples[[nm]][, , j] <- samples[[nm]][, , j - 1]
    }
  }
  samples$idx <- j
  samples
}

safe_new_particle <- function (s, data, pm_settings, eff_mu = NULL,
                               eff_var = NULL, chains_mu = NULL,
                               chains_var = NULL, prev_ll,
                               parameters, model = NULL, stage,
                               type, tune, marginalise = NULL, r_cores = 1) {
  attempt <- tryCatch(
    new_particle(s, data, pm_settings, eff_mu, eff_var, chains_mu, chains_var,
                 prev_ll, parameters, model, stage, type, tune, marginalise, r_cores),
    error = identity
  )
  if (inherits(attempt, c("error", "try-error"))) {
    return(reject_particle(parameters$alpha[, s], prev_ll, pm_settings))
  }
  attempt
}

reject_particle <- function(subj_mu, prev_ll, pm_settings) {
  list(proposal = subj_mu, ll = prev_ll, pm_settings = pm_settings)
}


new_particle <- function (s, data, pm_settings, eff_mu = NULL,
                          eff_var = NULL, chains_mu = NULL,
                          chains_var = NULL, prev_ll,
                          parameters, model = NULL, stage,
                          type, tune, marginalise = NULL, r_cores = 1)
{
  group_pars <- get_group_level(parameters, s, type)
  unq_components <- unique(tune$components)
  proposal_out <- numeric(length(group_pars$mu))
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- parameters$alpha[,s]
  # Node grid of the accepted particle, captured during the likelihood step and
  # reused for reconstruction (avoids a second full marginal pass at storage).
  marg_nodes <- NULL
  marg_terms_row <- NULL
  marginal_idx <- rep(FALSE, length(subj_mu))
  if (!is.null(marginalise)) {
    if (is.null(names(subj_mu))) {
      names(subj_mu) <- rownames(parameters$alpha)[seq_along(subj_mu)]
    }
    marginal_idx <- names(subj_mu) %in% marginalise$param
    # Keep a placeholder column for design mapping, but hold the marginalized
    # coordinate out of every proposal draw and proposal-density evaluation.
    group_mu[marginal_idx] <- marginalise$mu
    group_var[marginal_idx, ] <- 0
    group_var[, marginal_idx] <- 0
    group_var[marginal_idx, marginal_idx] <- marginalise$sigma^2
  }
  out_lls <- numeric(length(unq_components))
  particle_multiplier <- 1
  # Set the proposals
  if(stage == "preburn"){
    Mus <- list(group_mu, subj_mu)
    Sigmas <- list(group_var, group_var)
    # For preburn use a lot of proposals, to increase initial search a bit
    particle_multiplier <- 2
  } else if(stage == "burn"){ # Burn
    Mus <- list(group_mu, subj_mu, subj_mu)
    Sigmas <- list(group_var, group_var, chains_var)
  } else if(stage == "adapt"){
    Mus <- list(group_mu, subj_mu, chains_mu)
    Sigmas <- list(group_var, chains_var, chains_var)
  } else{ # Sample
    Mus <- list(group_mu, subj_mu, chains_mu, eff_mu)
    Sigmas <- list(group_var, chains_var, chains_var, eff_var)
  }
  n_proposals <- length(Mus)
  for(i in unq_components){
    # Add 1 to epsilons such that prior/group-level proposals aren't scaled
    epsilons <- c(1, pm_settings[[i]]$epsilon)
    idx_full <- tune$components == i
    idx <- idx_full & !marginal_idx
    p_idx <- sum(idx)

    Rs <- vector("list", n_proposals)
    rootis <- vector("list", n_proposals)
    log_consts <- numeric(n_proposals)
    for (k in seq_len(n_proposals)) {
      if (p_idx > 0) {
        base_R <- tryCatch(chol(Sigmas[[k]][idx,idx,drop=FALSE]), error = function(e) NULL)
        if (is.null(base_R)) {
          base_R <- tryCatch(chol(Sigmas[[k]][idx,idx,drop=FALSE] + diag(1e-6, p_idx)), error = function(e) NULL)
        }
        if (is.null(base_R)) {
          base_R <- tryCatch(chol(Sigmas[[k]][idx,idx,drop=FALSE] + diag(1e-4, p_idx)), error = function(e) NULL)
        }
        if (is.null(base_R)) {
          cov_diag <- diag(Sigmas[[k]][idx,idx,drop=FALSE])
          cov_diag[is.na(cov_diag) | cov_diag <= 0] <- 1e-4
          base_R <- diag(sqrt(cov_diag), p_idx)
        }
        Rs[[k]] <- base_R * epsilons[k]
        rootis[[k]] <- backsolve(Rs[[k]], diag(p_idx))
        log_consts[[k]] <- sum(log(diag(rootis[[k]]))) - 0.5 * p_idx * log(2 * pi)
      }
    }

    # Draw new proposals for each component
    particle_numbers <- numbers_from_proportion(pm_settings[[i]]$mix, pm_settings[[i]]$n_particles*particle_multiplier)
    proposals <- vector("list", n_proposals +1)
    proposals[[1]] <- matrix(subj_mu[idx], nrow = 1L)
    for(j in 1:n_proposals){
      # Fill up the proposals
      if (p_idx > 0) {
        proposals[[j + 1]] <- particle_draws(
          particle_numbers[j], Mus[[j]][idx],
          Sigmas[[j]][idx,idx,drop=FALSE] * (epsilons[j]^2), R = Rs[[j]])
      } else {
        proposals[[j + 1]] <- matrix(numeric(0), nrow = particle_numbers[j], ncol = 0L)
      }
    }
    proposals <- do.call(rbind, proposals)

    # Non -used proposals (for prior calculations)
    # Rejoin new proposals with current MCMC values for other components
    if(any(!idx)){
      proposals_other <- do.call(rbind, rep(list(subj_mu[!idx]), nrow(proposals)))
      colnames(proposals_other) <- names(subj_mu)[!idx]
      colnames(proposals) <- names(subj_mu)[idx]
      proposals <- cbind(proposals, proposals_other)
      proposals <- proposals[,names(subj_mu),drop=FALSE]
    } else{
      colnames(proposals) <- names(subj_mu)
    }

    # Normally we assume that a component contains all the parameters to estimate the individual likelihood of a joint model
    # Sometimes we may also want to block within a model if it has very high dimensionality
    shared_idx <- tune$shared_ll_idx[idx_full][1]
    is_shared <- shared_idx == tune$shared_ll_idx

    # Calculate likelihoods
    if (!is.null(marginalise)) {
      # Compute the t0 quadrature grid ONCE: reduce it to the marginal ll for
      # the MH weights, and stash the grid so the accepted particle's node
      # weights reconstruct t0 without a second full marginal pass.
      warm_state <- pm_settings[[i]]$marg_warm
      warm_arg <- if (!is.null(warm_state) && warm_state$backoff <= 0L) {
        c(mode = warm_state$mode, sd = warm_state$sd)
      } else NULL
      # The within-iteration hint (pilot a subset, regress its modes on the other
      # parameters) is rejected outright when the proposal cloud is too wide for
      # the regression to bracket, and then its subset pilot is wasted work. Back
      # it off exactly like the warm hint: the cloud narrows as the chain
      # converges, and it starts paying off from then on.
      pred_state <- pm_settings[[i]]$marg_pred
      use_pred <- is.null(pred_state) || pred_state$backoff <= 0L
      marg_spec <- marginalise
      if (!use_pred) marg_spec$predict_mode <- FALSE
      marg_grid <- compute_marginal_grid(proposals[, is_shared, drop = FALSE],
                                         data, model, marg_spec,
                                         r_cores = r_cores, warm = warm_arg)
      lw <- marginal_ll_from_grid(marg_grid)
    } else if(tune$components[length(tune$components)] > 1){
      lw <- calc_ll_manager(proposals[,is_shared], dadm = data, model,
                            component = shared_idx, r_cores = r_cores)
    } else{
      lw <- calc_ll_manager(proposals[,is_shared], dadm = data, model,
                            r_cores = r_cores)
    }
    lw_total <- lw + prev_ll - lw[1] # make sure lls from other components are included
    # Prior density
    lp <- if (p_idx > 0) {
      fast_dmvnorm_rooti(x = proposals[,idx,drop=FALSE], mean = group_mu[idx],
                         rooti = rootis[[1]], log_const = log_consts[[1]])
    } else rep(0, nrow(proposals))
    if(length(unq_components) > 1){
      prior_density <- if (any(!marginal_idx)) {
        fast_dmvnorm(x = proposals[,!marginal_idx,drop=FALSE],
                     mean = group_mu[!marginal_idx],
                     sigma = group_var[!marginal_idx,!marginal_idx,drop=FALSE])
      } else rep(0, nrow(proposals))
    } else{
      prior_density <- lp
    }
    # Calculate mixture log-density using log-sum-exp for numerical stability
    log_mix_comps <- matrix(0, nrow = nrow(proposals), ncol = n_proposals)
    log_mix_comps[, 1] <- log(pm_settings[[i]]$mix[1]) + lp
    for (k in 2:n_proposals) {
      if (p_idx > 0) {
        log_mix_comps[, k] <- log(pm_settings[[i]]$mix[k]) + fast_dmvnorm_rooti(
          x = proposals[, idx, drop = FALSE], mean = Mus[[k]][idx],
          rooti = rootis[[k]], log_const = log_consts[[k]])
      } else {
        log_mix_comps[, k] <- log(pm_settings[[i]]$mix[k])
      }
    }
    max_log <- do.call(pmax, as.data.frame(log_mix_comps))
    lm <- max_log + log(rowSums(exp(log_mix_comps - max_log)))
    infnt_idx <- is.infinite(lm) | is.na(lm)
    if (any(infnt_idx)) {
      fin_vals <- lm[!infnt_idx]
      lm[infnt_idx] <- if (length(fin_vals) > 0) min(fin_vals) else -1e10
    }
    # Calculate weights and center
    l <- lw_total + prior_density - lm
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = sum(particle_numbers) + 1, size = 1, prob = weights)

    out_lls[i] <- lw[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,idx]
    if (!is.null(marginalise)) {
      marg_nodes <- as.numeric(marg_grid$nodes[idx_ll, ])
      marg_terms_row <- as.numeric(marg_grid$log_terms[idx_ll, ])
      warm_next <- marginal_warm_state(marg_grid, idx_ll)
      warm_next_state <- marginal_warm_backoff(
        warm_state, attempted = !is.null(warm_arg),
        used = isTRUE(marg_grid$warm_used >= 0.5))
      # The predictor only runs when the warm hint did not already carry the step.
      pred_next_state <- marginal_warm_backoff(
        pred_state, attempted = use_pred && !isTRUE(marg_grid$warm_used == 1),
        used = isTRUE(marg_grid$pred_used == 1))
    }
    pm_settings[[i]] <- update_pm_settings(pm_settings[[i]], idx_ll, weights, particle_numbers, tune, sum(idx))
    # Seed next iteration's quadrature from the particle that was accepted here.
    if (!is.null(marginalise)) {
      pm_settings[[i]]$marg_warm <- if (is.null(warm_next)) NULL else {
        c(as.list(warm_next), warm_next_state)
      }
      pm_settings[[i]]$marg_pred <- pred_next_state
    }
  }
  names(proposal_out) <- names(subj_mu)
  if (!is.null(marginalise)) {
    p_idx <- match(marginalise$param, names(proposal_out))
    if (!is.na(p_idx)) {
      proposal_out[p_idx] <- if (!is.null(marg_terms_row)) {
        draw_marginal_node(marg_nodes, marg_terms_row, marginalise)
      } else {
        stats::rnorm(1L, marginalise$mu, marginalise$sigma)
      }
    }
  }
  return(list(proposal = proposal_out, ll = sum(out_lls), pm_settings = pm_settings))
}


update_pm_settings <- function(pm_settings, chosen_idx, weights, particle_numbers,
                               tune, n_pars) {
  # 0) If we're past an initial burn-in, do the adaptation
  pm_settings$iter <- pm_settings$iter + 1
  if (pm_settings$iter > tune$n0) {
    # A) Update proposal_counts
    # -------------------------
    # Each proposal j was used "particle_numbers[j]" times
    pm_settings$proposal_counts <- pm_settings$proposal_counts + particle_numbers

    # B) Update acceptance counts: "percent better than old"
    # ------------------------------------------------------
    # weights[1] = old particle's weight
    # weights[2..(1+sum(particle_numbers))] = new draws' weights
    old_weight <- weights[1]

    # We'll parse out each proposal's chunk in weights[-1]
    # using the known counts in particle_numbers.
    # We're also tracking acceptance of group-level proposals, which is minorly wasteful
    offset <- 2  # start index in 'weights' for new proposals
    for (j in seq_along(particle_numbers)) {
      # The chunk of new weights for proposal j
      n_j <- particle_numbers[j]
      if (n_j > 0) {
        draws_j <- weights[offset:(offset + n_j - 1)]
        # Count how many draws_j exceed old_weight
        better_j <- sum(draws_j > old_weight)
        # Accumulate that in acceptance counts
        pm_settings$acc_counts[j] <- pm_settings$acc_counts[j] + better_j
        offset <- offset + n_j
      }
    }

    # C) Compute per-proposal acceptance rates
    # ----------------------------------------
    acc_rates <- ifelse(pm_settings$proposal_counts > 0, pm_settings$acc_counts / pm_settings$proposal_counts, 0)

    # .1 for preburn, .4 for burn and adapt and .6 for sample
    clamp_min <- ifelse(length(pm_settings$mix) == 2, .1, ifelse(length(pm_settings$mix) == 3, .4, .6))

    # D) Adapt epsilon via continuous approach
    # ----------------------------------------
    # pm_settings$epsilon is a vector, same length as pm_settings$mix
    # tune$p_accept is also a vector, e.g. c(0.2, 0.3, 0.6) for each proposal
    new_epsilon <- update_epsilon_continuous(
      epsilon   = pm_settings$epsilon,
      acceptance = acc_rates[-1],
      target     = tune$p_accept,
      iter       = pm_settings$iter,
      d          = n_pars,
      alphaStar  = tune$alphaStar,
      damp       = 100,          # Example
      clamp      = c(clamp_min, 5)    # Example range
    )
    pm_settings$epsilon <- new_epsilon

    # E) Adapt mixing weights based on acceptance vs. target
    # ------------------------------------------------------
    # If ratio_j > 1 (proposal j acceptance > target), its mix goes up;
    # if ratio_j < 1, mix goes down.

    if(length(pm_settings$mix) > 2){ # We're not in preburn
      eps_val <- 1e-12   # Avoid divide-by-zero

      # 1) Compute performance ~ (acceptance / old_mix), normalized
      performance <- (acc_rates + eps_val) / pm_settings$mix
      performance <- performance / sum(performance)

      # 2) Adjust the elements by expected acceptance (p_accept)
      # The first element is the group-level, so has no p_accept, just take the mean of the others
      # Bit hacky
      adj_factor <- c(mean(tune$p_accept), tune$p_accept)
      performance <- performance / adj_factor

      # 3) Normalize again
      performance <- performance / sum(performance)

      # 4) Blend with old mix: stable update
      new_mix <- (1 - tune$mix_adapt) * pm_settings$mix + tune$mix_adapt * performance

      # 5) Impose a floor, re-normalize
      new_mix <- pmax(new_mix, 0.02)
      new_mix <- new_mix / sum(new_mix)
      pm_settings$mix <- new_mix
    }


    # F) Adapt the number of particles (ESS logic)
    # -------------------------------------------------------
    # If length mix > 2, we're in sample stage
    # Only reduce number of particles when we're already converged
    if (length(pm_settings$mix) > 3 && pm_settings$gd_good) {
      ess <- sum(weights)^2 / sum(weights^2)
      desired_ess <- tune$target_ESS
      scale_factor <- (desired_ess / ess)^tune$ESS_scale
      new_num_particles <- round(pm_settings$n_particles * scale_factor)
      pm_settings$n_particles <- max(25, min(tune$max_particles, new_num_particles))
    }
  }

  return(pm_settings)
}



# Utility functions for sampling below ------------------------------------
update_epsilon_continuous <- function(
    epsilon,    # vector of current epsilons
    acceptance, # vector of acceptance rates, same length
    target,     # vector of target acceptance rates, same length
    iter,
    d,
    alphaStar,
    damp = 100,
    clamp = c(0.6, 4)
) {
  log_eps <- log(epsilon)
  # 2) define step size
  # We'll do one pass per element. If you want a single c_term, that's also fine.
  c_term <- (1 - 1/d)*sqrt(2*pi)*exp(alphaStar^2/2)/(2*alphaStar) + 1/(d*target*(1-target))
  step_size <- c_term / max(damp, iter)
  # 3) compute difference from target acceptance
  diff_accept <- acceptance - target
  # 4) update in log space (vectorized)
  log_eps_new <- log_eps + step_size * diff_accept
  # 5) exponentiate
  eps_new <- exp(log_eps_new)
  # 6) clamp
  eps_new <- pmin(clamp[2], pmax(eps_new, clamp[1]))

  return(eps_new)
}

update_epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c <- ((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta <- log(sqrt(epsilon2))
  Theta <- Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}


numbers_from_proportion <- function(mix_proportion, n_particles = 1000) {
  # Make sure each proposal has at least 1 particle
  return(pmax(1, rmultinom(1, n_particles, mix_proportion)))
}


particle_draws <- function(n, mu, covar, alpha = NULL, tau= NULL, R = NULL) {
  if (n <= 0) {
    return(matrix(numeric(0), nrow = 0, ncol = length(mu)))
  }
  if(is.null(alpha)){
    if (!is.null(R)) {
      p <- length(mu)
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      X <- Z %*% R
      return(sweep(X, 2, mu, "+"))
    } else {
      return(mvtnorm::rmvnorm(n, mu, covar))
    }
  }
}

extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. Big shout out to the rapply function
  sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
  if(any(sampler$nuisance)) sampler$sampler_nuis$samples <- rapply(sampler$sampler_nuis$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  return(sampler)
}

extend_obj <- function(obj, n_extend){
  old_dim <- dim(obj)
  n_dimensions <- length(old_dim)
  if(is.null(old_dim) | n_dimensions == 1) return(obj)
  if(n_dimensions == 2){
    if(nrow(obj) == ncol(obj)){
      if(nrow(obj) > 1){
        if(mean(abs(abs(rowSums(obj/max(obj))) - abs(colSums(obj/max(obj))))) < .01) return(obj)
      }
    }
  }
  new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
  extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
  extended[slice.index(extended,n_dimensions) <= old_dim[n_dimensions]] <- obj
  return(extended)
}

sample_store_base <- function(data, par_names, iters = 1, stage = "init", is_nuisance = rep(F, length(par_names)), ...) {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  samples <- list(
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

block_variance_idx <- function(components){
  vars_out <- matrix(0, length(components), length(components))
  for(i in unique(components)){
    idx <- i == components
    vars_out[idx,idx] <- NA
  }
  return(vars_out == 0)
}

fill_samples_base <- function(samples, group_level, proposals, j = 1, n_pars){
  # Fill samples both group level and random effects
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, j, n_pars)
  return(samples)
}



fill_samples_RE <- function(samples, proposals, j = 1, n_pars, ...){
  # Only for random effects, separated because group level sometimes differs.
  if(!is.null(proposals)){
    samples$alpha[, , j] <- proposals[1:n_pars,]
    samples$subj_ll[, j] <- proposals[n_pars + 1,]
    samples$idx <- j
  }
  return(samples)
}


set_p_accept <- function(stage, search_width){
  # Proposals distributions:
  # 1. Prior - unscaled: all stages
  # 2. Prev particle - scaled chain variance: all stages in preburn scaled by prior variance
  # 3. Chain mean - scaled chain variance: burn onwards
  # 4. Eff mean - scaled eff variance: sample onwards
  if(stage == "preburn") return(0.02 * (1/search_width))
  if(stage == "burn") return(c(0.02, 0.25)* (1/search_width))
  if(stage == "adapt") return(c(0.2, 0.25)* (1/search_width))
  if(stage == "sample") return(c(0.3, 0.3, 0.3)* (1/search_width))
}

get_default_mix <- function(stage){
  if (stage == "burn") {
    default_mix <- c(0.15, 0.35, 0.5)
  } else if(stage == "adapt"){
    default_mix <- c(0.1, 0.4, 0.4)
  } else if(stage == "sample"){
    default_mix <- c(0.05, 0.3, 0.3, 0.35)
  }  else{
    default_mix <- c(0.5, 0.5)
  }
  return(default_mix)
}


check_mix <- function(mix = NULL, stage) {
  default_mix <- get_default_mix(stage)
  if(is.null(mix)) mix <- default_mix
  if(length(mix) < length(default_mix)){
    if(length(default_mix) - length(mix) == 1){
      mix <- mix*(1-default_mix[length(default_mix)])
      mix <- c(mix, default_mix[length(default_mix)])
    } else{
      stop("mix settings are not compatible with stage")
    }
  }
  return(mix)
}

check_epsilon <- function(epsilon, n_pars, mix) {
  if (is.null(epsilon)) { # In preburn case, there's only one epsilon here
    if (n_pars > 15) {
      epsilon <- .5
    } else if (n_pars > 10) {
      epsilon <- .6
    } else {
      epsilon <- .7
    }
    # The first proposal is always unscaled
    # Every subsequent phase at most adds one epsilon
  } else if(length(epsilon) < (length(mix) -1)){
    epsilon <- c(epsilon, epsilon[length(epsilon)])
  }
  return(epsilon)
}

check_prop_performance <- function(prop_performance, stage){
  default_mix <- get_default_mix(stage)
  if(is.null(prop_performance) || stage == "adapt") prop_performance <- rep(0, length(default_mix))
  if(length(prop_performance) < length(default_mix)){
    if(length(default_mix) - length(prop_performance) == 1){
      n_total <- sum(prop_performance)
      prop_performance <- prop_performance*(1-default_mix[length(default_mix)])
      prop_performance <- c(prop_performance, default_mix[length(default_mix)]*n_total)
    } else{
      stop("prop_performance settings are not compatible with stage")
    }
  }
  return(round(prop_performance))
}

calc_ll_manager <- function(proposals, dadm, model, component = NULL, r_cores = 1,
                            marginalise = NULL){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component, r_cores = r_cores, marginalise = marginalise)
  } else{
    model <- model()
    dadm <- .cache_ll_data_attrs(dadm)
    # marginalise is threaded explicitly by the sampler; it is never inferred
    # from a dadm attribute, so predict/make_data/IC calls never integrate.
    if(is.null(model$c_name)){ # use the R implementation
      if (!is.null(marginalise)) {
        stop("marginalise requires a registered race-model likelihood")
      }
      lls <- unlist(
        auto_mclapply(1:nrow(proposals),
          function(i) calc_ll_R(proposals[i,], model=model, dadm = dadm),
         mc.cores=r_cores))
    } else {
      # SS models: push stop_method/stop_n_nodes into the process-global C++
      # config once per likelihood call (before any fork, so mclapply workers
      # inherit it; PSOCK workers run this themselves)
      set_stop_method_from_model(model)
      p_types <- names(model$p_types)
      designs <- .oo_expanded_designs(dadm)
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      if (nrow(proposals) <= r_cores) {
        lls <- calc_ll_oo(proposals, dadm, constants = constants, designs = designs,
                          type = model$c_name, bounds = model$bound,
                          transforms = model$transform, pretransforms = model$pre_transform,
                          p_types = p_types, min_ll = log(1e-10), trend = model$trend,
                          marginalise = marginalise)
      } else {
        idx <- rep(1:r_cores,each=1+(nrow(proposals) %/% r_cores))[1:nrow(proposals)]
        lls <- unlist(auto_mclapply(1:r_cores,function(i) {
          calc_ll_oo(proposals[idx==i,,drop=FALSE], dadm, constants = constants,
                     designs = designs, type = model$c_name, bounds = model$bound,
                     transforms = model$transform, pretransforms = model$pre_transform,
                     p_types = p_types, min_ll = log(1e-10),
                     trend = model$trend, marginalise = marginalise)
        },mc.cores=r_cores))
      }
    }
  }
  return(lls)
}

.waic_trial_row_groups <- function(dadm, model_type = NULL) {
  n_rows <- nrow(dadm)
  if (n_rows == 0) return(list(integer(0)))
  if (is.null(model_type)) model_type <- ""
  if (!("lR" %in% names(dadm)) || grepl("DDM", model_type)) {
    return(as.list(seq_len(n_rows)))
  }

  lR_codes <- as.integer(dadm$lR)
  first_code <- suppressWarnings(min(lR_codes, na.rm = TRUE))
  if (!is.finite(first_code)) return(as.list(seq_len(n_rows)))

  starts <- which(lR_codes == first_code)
  if (length(starts) == 0) return(as.list(seq_len(n_rows)))
  ends <- c(starts[-1] - 1L, n_rows)
  Map(seq.int, starts, ends)
}

.waic_subset_dadm <- function(dadm, rows) {
  out <- dadm[rows, , drop = FALSE]
  designs <- attr(dadm, "designs")
  if (!is.null(designs)) {
    out_designs <- lapply(designs, function(dm) {
      dm_out <- dm[rows, , drop = FALSE]
      dm_expand <- attr(dm, "expand")
      if (!is.null(dm_expand)) attr(dm_out, "expand") <- dm_expand[rows]
      dm_out
    })
    attr(out, "designs") <- out_designs
  }
  attr(out, "constants") <- attr(dadm, "constants")
  attr(out, "expand") <- 1L
  out
}

#' Calculate Pointwise Log-Likelihoods
#'
#' @param proposals A matrix of parameter proposals with dimensions `n_iter x n_pars`.
#' @param dadm A Data Augmented Design Matrix (DADM).
#' @param model A model function or object.
#' @param r_cores Number of cores to use (for non-C models).
#'
#' @return A matrix of pointwise log-likelihoods with dimensions `n_iter x n_trials`.
#' @export
calc_ll_pw <- function(proposals, dadm, model, r_cores = 1){
  model <- model()
  dadm <- .cache_ll_data_attrs(dadm)
  c_name <- model$c_name
  unsupported <- c("MRI", "MRI_AR1", "SSEXG", "SSRDEX", "SOFTMAX")
  if(!is.null(c_name) && c_name %in% unsupported)
    stop("WAIC is not supported for model type '", c_name, "'")

  if (!is.null(c_name) && exists("calc_ll_oo_pw", mode = "function", inherits = TRUE)) {
    p_types <- names(model$p_types)
    designs <- .oo_expanded_designs(dadm)
    constants <- attr(dadm, "constants")
    if (is.null(constants)) constants <- NA
    return(calc_ll_oo_pw(proposals, dadm, constants = constants, designs = designs,
                         type = c_name, model$bound, model$transform, model$pre_transform,
                         p_types = p_types, min_ll = log(1e-10), model$trend))
  }

  # Fallback for non-C models or environments without calc_ll_oo_pw.
  trial_groups <- .waic_trial_row_groups(dadm, c_name)
  if (length(trial_groups) == 0) {
    return(matrix(numeric(0), nrow = nrow(proposals), ncol = 0))
  }

  model_fun <- function() model
  ll_unique <- vapply(trial_groups, function(rows) {
    trial_dadm <- .waic_subset_dadm(dadm, rows)
    calc_ll_manager(proposals, trial_dadm, model = model_fun, r_cores = r_cores)
  }, numeric(nrow(proposals)))
  if (is.null(dim(ll_unique))) ll_unique <- matrix(ll_unique, ncol = 1)

  expand_map <- attr(dadm, "expand")
  if (is.null(expand_map)) return(ll_unique)
  expand_map <- as.integer(expand_map)
  if (anyNA(expand_map) || any(expand_map < 1L) || any(expand_map > ncol(ll_unique))) {
    stop("Invalid 'expand' mapping while computing WAIC pointwise log-likelihoods.")
  }
  ll_unique[, expand_map, drop = FALSE]
}

merge_group_level <- function(tmu, tmu_nuis, tvar, tvar_nuis, is_nuisance, subj_mu){
  n_pars <- length(is_nuisance)
  tmu_out <- numeric(n_pars)
  tmu_out[!is_nuisance] <- tmu
  tmu_out[is_nuisance] <- tmu_nuis
  tvar_out <- matrix(0, nrow = n_pars, ncol = n_pars)
  tvar_out[!is_nuisance, !is_nuisance] <- tvar

  subj_mu_out <- matrix(NA, ncol = ncol(subj_mu), nrow = length(tmu_out))
  subj_mu_out[is_nuisance,] <- do.call(cbind, rep(list(c(tmu_nuis)), ncol(subj_mu)))
  subj_mu_out[!is_nuisance,] <- subj_mu
  return(list(tmu = tmu_out, tvar = tvar_out, subj_mu = subj_mu_out))
}


#' Run a Group-level Model.
#'
#' Separate function for running only the group-level model. This can be useful in a
#' two-step analysis. Works similar in functionality to make_emc,
#' except also does the fitting and returns an emc object that works with
#' most posterior checking tests (but not the data generation/posterior predictives).
#'
#' @param prior an emc.prior object.
#' @param iter Number of MCMC samples to collect.
#' @inheritParams make_emc
#'
#' @returns an emc object with only group-level samples
#' @export
run_hyper <- function(type = "standard", data, prior = NULL, iter = 1000, n_chains =3, ...){
  args <- list(...)
  if(length(dim(data)) == 3){
    data_input <- data
    data <- as.data.frame(t(data_input[,,1]))
    data$subjects <- 1:nrow(data)
    iter <- dim(data_input)[3]
    is_mcmc <- T
    pars <- rownames(data_input)
  } else{
    data_input <- data[,colnames(data)!= "subjects"]
    is_mcmc <- F
    pars <- colnames(data_input)
  }
  emc <- list()
  for(j in 1:n_chains){
    samples <- sample_store(data = data ,par_names = pars, is_nuisance = rep(F, length(pars)), integrate = F, type = type, ...)
    subjects <- unique(data$subjects)
    sampler <- list(
      data = split(data, data$subjects),
      par_names = pars,
      subjects = subjects,
      n_pars = length(pars),
      nuisance = rep(F, length(pars)),
      n_subjects = length(subjects),
      samples = samples,
      init = TRUE
    )
    class(sampler) <- "pmwgs"
    sampler <- add_info(sampler, prior, type = type, ...)
    sampler$type <- type
    startpoints <- get_startpoints(sampler, start_mu = NULL, start_var = NULL, type = type)
    sampler$samples <- fill_samples(samples = sampler$samples, group_level = startpoints, proposals = NULL,
                                    j = 1, n_pars = sampler$n_pars, type = type)
    sampler$samples$idx <- 1
    sampler <- extend_sampler(sampler, iter-1, "sample")
    for(i in 2:iter){
      if(is_mcmc){
        group_pars <- gibbs_step(sampler, data_input[,,i], type = type)
      } else{
        group_pars <- gibbs_step(sampler, t(data_input), type = type)
      }
      sampler$samples$idx <- i
      sampler$samples <- fill_samples(samples = sampler$samples, group_level = group_pars, proposals = NULL,
                                                   j = i, n_pars = sampler$n_pars, type = type)
    }
    emc[[j]] <- sampler
  }
  emc[[1]]$type <- type
  class(emc) <- "emc"
  emc <- subset(emc, filter = 1)
  return(emc)
}

check_CR <- function(emc, p_vector, range = .2, N = 500){
  covs <- diag(length(p_vector)) * range
  props <- mvtnorm::rmvnorm(N, mean = p_vector, sigma = covs)
  model <- emc[[1]]$model
  if(is.null(model()$c_name)) stop("C not implemented yet for this model")
  dat <- emc[[1]]$data[[1]]
  modelRlist <- model()
  modelRlist$c_name <- NULL
  modelR <- function()return(modelRlist)
  t1 <- system.time(
    R <- calc_ll_manager(props, dat, modelR)
  )
  t2 <- system.time(
    C <- calc_ll_manager(props, dat, model)
  )
  print(paste0("C ", t1$elapsed/t2$elapsed, " times faster"))
  if(!identical(C, R)){
    warning("C and R results differ")
  }
  return(list(C = C, R = R))
}

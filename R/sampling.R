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

compute_marginal_grid <- function(proposals, data, model, marginalise,
                                  r_cores = 1, warm = NULL) {
  model_spec <- if (is.function(model)) model() else model
  set_stop_method_from_model(model_spec)
  data <- .cache_ll_data_attrs(data)
  constants <- attr(data, "constants")
  if (is.null(constants)) constants <- NA
  designs <- .oo_expanded_designs(data, expand = FALSE)
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
  if (r_cores <= 1 || nrow(proposals) <= r_cores) return(one(proposals))
  idx <- .split_work_indices(nrow(proposals), r_cores)
  parts <- auto_mclapply(1:r_cores, function(i) {
    one(proposals[idx == i, , drop = FALSE])
  }, mc.cores = r_cores)
  list(nodes = do.call(rbind, lapply(parts, `[[`, "nodes")),
       log_terms = do.call(rbind, lapply(parts, `[[`, "log_terms")),
       ll = unlist(lapply(parts, `[[`, "ll")),
       mode = unlist(lapply(parts, `[[`, "mode")),
       sd = unlist(lapply(parts, `[[`, "sd")),
       warm_used = mean(unlist(lapply(parts, `[[`, "warm_used"))),
       pred_used = mean(unlist(lapply(parts, `[[`, "pred_used"))),
       repaired = sum(unlist(lapply(parts, `[[`, "repaired"))))
}

.split_work_indices <- function(n, n_workers) {
  n <- as.integer(n)
  n_workers <- min(as.integer(n_workers), n)
  if (n <= 0L || n_workers <= 0L) return(integer(0))
  q <- n %/% n_workers
  rem <- n %% n_workers
  rep.int(seq_len(n_workers), q + (seq_len(n_workers) <= rem))
}

marginal_warm_state <- function(grid, idx) {
  if (is.null(grid$mode) || is.null(grid$sd)) return(NULL)
  m <- grid$mode[idx]; s <- grid$sd[idx]
  if (!is.finite(m) || !is.finite(s) || s <= 0) return(NULL)
  c(mode = m, sd = s)
}

marginal_warm_backoff <- function(state, attempted, used) {
  pen <- if (is.null(state$penalty)) 0L else state$penalty
  wait <- if (is.null(state$backoff)) 0L else state$backoff
  if (attempted) {
    if (isTRUE(used)) { pen <- 0L; wait <- 0L }
    else { pen <- min(16L, max(1L, pen * 2L)); wait <- pen }
  } else wait <- max(0L, wait - 1L)
  list(penalty = pen, backoff = wait)
}

marginal_ll_from_grid <- function(log_terms) {
  if (is.list(log_terms) && !is.null(log_terms$ll)) return(log_terms$ll)
  m <- rep(-Inf, nrow(log_terms))
  for (k in seq_len(ncol(log_terms))) {
    col <- log_terms[, k]
    col[is.na(col)] <- -Inf
    m <- pmax(m, col)
  }
  m[!is.finite(m)] <- -Inf
  fin <- is.finite(m)
  res <- rep(-Inf, nrow(log_terms))
  if (any(fin)) {
    res[fin] <- m[fin] + log(rowSums(exp(log_terms[fin, , drop = FALSE] - m[fin]), na.rm = TRUE))
  }
  res
}

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

pmwgs <- function(dadm, type, pars = NULL, prior = NULL,
                  nuisance = NULL, nuisance_non_hyper = NULL, marginalise = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  if(is.null(pars)) pars <- names(sampled_pars(attr(dadm[[1]], "prior")))
  if(is.null(prior)) prior <- attr(dadm[[1]], "prior")
  marginalise <- resolve_marginalise_prior(marginalise, prior)
  dadm <- extractDadms(dadm, par_names = pars)

  dadm_list <-dadm$dadm_list
  components <- .align_sampling_index(attr(dadm_list, "components"), pars,
                                       "components", remap = TRUE)
  shared_ll_idx <- .align_sampling_index(attr(dadm_list, "shared_ll_idx"), pars,
                                         "shared_ll_idx", remap = FALSE)
  attr(dadm_list, "components") <- components
  attr(dadm_list, "shared_ll_idx") <- shared_ll_idx
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

.emc_blas_threads <- function() {
  cached <- .emc_pool_state$blas_threads
  if (!is.null(cached)) return(cached)
  n <- 1L
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    got <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
    if (length(got) == 1L && !is.na(got) && got >= 1L) n <- as.integer(got)
  }
  .emc_pool_state$blas_threads <- n
  n
}

.emc_set_blas_threads <- function(n) {
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) return(NULL)
  n <- max(1L, as.integer(n))
  previous <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
  ok <- tryCatch({ RhpcBLASctl::blas_set_num_threads(n); TRUE },
                 error = function(e) FALSE)
  if (!ok) return(NULL)
  .emc_pool_state$blas_threads <- n
  previous
}

.emc_limit_blas_threads <- function(total_cores) {
  total_cores <- max(1L, as.integer(total_cores)[1L])
  requested <- max(1L, as.integer(.emc_blas_threads())[1L])
  target <- min(requested, total_cores)
  if (requested > target) {
    .emc_set_blas_threads(target)
    observed <- max(1L, as.integer(.emc_blas_threads())[1L])
    if (observed > target) {
      if (!isTRUE(.emc_pool_state$blas_warning)) {
        warning("BLAS thread width exceeds the chain core budget; using a serial allocation")
        .emc_pool_state$blas_warning <- TRUE
      }
      return(target)
    }
    return(observed)
  }
  requested
}

.particle_core_budget <- function(n_subjects, n_cores = 1L, r_cores = 1L,
                                 total_cores = NULL, blas_threads = NULL,
                                 allow_spare = FALSE) {
  n_subjects <- as.integer(n_subjects)
  n_cores <- max(1L, as.integer(n_cores))
  r_cores <- max(1L, as.integer(r_cores))
  if (is.null(total_cores)) {
    if (n_subjects == 1L && n_cores > 1L) {
      return(list(subject = 1L, likelihood = max(r_cores, n_cores)))
    }
    return(list(subject = n_cores, likelihood = r_cores))
  }
  total_cores <- max(1L, as.integer(total_cores))
  if (is.null(blas_threads)) {
    blas_threads <- .emc_limit_blas_threads(total_cores)
  } else {
    blas_threads <- min(total_cores, max(1L, as.integer(blas_threads)[1L]))
  }
  budget <- max(1L, total_cores %/% blas_threads)
  if (n_subjects == 1L && n_cores > 1L) {
    return(list(subject = 1L,
                likelihood = min(budget, max(r_cores, n_cores))))
  }
  likelihood <- min(r_cores, budget)
  subject <- max(1L, min(n_subjects, n_cores, budget %/% likelihood))
  if (isTRUE(allow_spare) || r_cores > 1L) {
    likelihood <- max(likelihood, budget %/% subject)
  }
  list(subject = subject, likelihood = likelihood)
}

.emc_ll_route <- function(wpool, llpool, pool_budget, core_budget = NULL,
                          ll_route = NULL) {
  if (!is.null(wpool) && isTRUE(wpool$alive)) {
    return("persistent_pool")
  }
  if (!is.null(llpool)) {
    if (!is.null(ll_route) && !is.null(ll_route$last_route)) {
      return(ll_route$last_route)
    }
    if (isTRUE(llpool$alive)) return("persistent_pool")
    return(if (pool_budget$likelihood > 1L) "nested_per_call" else "serial")
  }
  if (!is.null(wpool)) {
    return(if (pool_budget$likelihood > 1L) "nested_per_call" else "serial")
  }
  if (!is.null(core_budget) && core_budget$likelihood > 1L) {
    return("nested_per_call")
  }
  "serial"
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000,
                 n_cores = 1, r_cores = 1, total_cores = NULL) {
  blas_before <- .emc_blas_threads()
  on.exit(if (!is.null(blas_before)) .emc_set_blas_threads(blas_before), add = TRUE)
  pmwgs$samples <- emc_clone_sample_store(pmwgs$samples)
  if (!is.null(pmwgs$sampler_nuis$samples)) {
    pmwgs$sampler_nuis$samples <- emc_clone_sample_store(pmwgs$sampler_nuis$samples)
  }
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
  core_budget <- .particle_core_budget(
    pmwgs$n_subjects, n_cores = n_cores, r_cores = r_cores,
    total_cores = total_cores
  )
  start_streams <- .emc_subject_streams(pmwgs$n_subjects)
  proposals <- .emc_with_preserved_rng(
    parallel::mclapply(X=1:pmwgs$n_subjects,
                       FUN=function(s, ...) {
                         assign(".Random.seed", start_streams[[s]],
                                envir = globalenv())
                         start_proposals(s, ...)
                       },
                       parameters = startpoints_comb, n_particles = particles,
                       pmwgs = pmwgs, type = type,
                       mc.cores = core_budget$subject,
                       r_cores = core_budget$likelihood))
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))


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
           n_cores = cores_per_chain, total_cores = cores_per_chain,
           mc.cores=cores_for_chains)
  class(emc) <- "emc"
  return(emc)
}

start_proposals <- function(s, parameters, n_particles, pmwgs, type, r_cores = 1){
  group_pars <- get_group_level(parameters, s, type)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha)
  data_s <- pmwgs$data[[which(pmwgs$subjects == s)]]
  marginalise <- pmwgs$marginalise
  if (!is.null(marginalise)) {
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
  tune$alphaStar <- ifelse(stage == "sample", 2, 3)
  if (is.null(tune$search_width)) tune$search_width <- 1
  tune$p_accept <- set_p_accept(stage, tune$search_width)
  if(is.null(tune$components)) tune$components <- rep(1, n_pars)
  if(is.null(tune$shared_ll_idx)) tune$shared_ll_idx <- tune$components
  if(is.null(tune$target_ESS)) tune$target_ESS <- 2.5*sqrt(n_pars)
  if(is.null(tune$ESS_scale)) tune$ESS_scale <- .05
  if(is.null(tune$max_particles)) tune$max_particles <- particles*1.2
  if(is.null(tune$mix_adapt)) tune$mix_adapt <- .05
  tune$n0 <- 25
  return(tune)
}

check_sampling_settings <- function(pm_settings, stage, n_pars, particles){
  for(i in 1:length(pm_settings)){
    pm_settings[[i]]$mix <- check_mix(pm_settings[[i]]$mix, stage)
    pm_settings[[i]]$epsilon <- check_epsilon(pm_settings[[i]]$epsilon, n_pars, pm_settings[[i]]$mix)
    pm_settings[[i]]$proposal_counts <- check_prop_performance(pm_settings[[i]]$proposal_counts, stage)
    pm_settings[[i]]$acc_counts <- check_prop_performance(pm_settings[[i]]$acc_counts, stage)
    if(is.null(pm_settings[[i]]$n_particles)) pm_settings[[i]]$n_particles <- particles
    if(is.null(pm_settings[[i]]$iter)) pm_settings[[i]]$iter <- 1
    if(is.null(pm_settings[[i]]$gd_good)) pm_settings[[i]]$gd_good <- FALSE
  }
  return(pm_settings)
}

.align_sampling_index <- function(index, par_names, label, remap = FALSE) {
  n_pars <- length(par_names)
  if (is.null(index)) return(NULL)
  index_names <- names(index)
  if (!is.null(index_names) && length(index_names)) {
    if (anyDuplicated(index_names) || anyDuplicated(par_names)) {
      stop("Parameter bookkeeping error: duplicate names in ", label,
           " or pmwgs$par_names")
    }
    missing <- setdiff(par_names, index_names)
    extra <- setdiff(index_names, par_names)
    if (length(missing) || length(extra)) {
      stop("Parameter bookkeeping error: named ", label,
           " does not match pmwgs$par_names")
    }
    index <- index[match(par_names, index_names)]
  } else if (length(index) != n_pars) {
    stop("Parameter bookkeeping error: ", label, " has length ",
         length(index), " but pmwgs$par_names has length ", n_pars)
  }
  if (length(index) != n_pars) {
    stop("Parameter bookkeeping error: aligned ", label,
         " does not have covariance dimension ", n_pars)
  }
  labels <- as.numeric(index)
  if (anyNA(index) || any(!is.finite(labels)) ||
      any(labels != as.integer(labels)) || any(labels < 1L)) {
    stop("Parameter bookkeeping error: ", label,
         " must contain positive integer-like component labels")
  }
  out <- if (isTRUE(remap)) {
    match(as.integer(labels), unique(as.integer(labels)))
  } else {
    as.integer(labels)
  }
  names(out) <- par_names
  out
}



run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      n_cores = 1,
                      tune = NULL,
                      verbose = TRUE,
                      verboseProgress = TRUE,
                      r_cores = 1,
                      core_ctl = NULL) {
  blas_before <- .emc_blas_threads()
  on.exit(if (!is.null(blas_before)) .emc_set_blas_threads(blas_before), add = TRUE)
  n_pars <- pmwgs$n_pars
  par_names <- pmwgs$par_names
  if (is.null(par_names)) {
    alpha_dimnames <- dimnames(pmwgs$samples$alpha)
    par_names <- if (is.null(alpha_dimnames)) NULL else alpha_dimnames[[1L]]
  }
  if (is.null(par_names) || length(par_names) != n_pars) {
    stop("Parameter bookkeeping error: pmwgs$par_names does not match covariance dimension")
  }
  if (is.null(tune)) tune <- list()
  components <- attr(pmwgs$data, "components")
  shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")
  if (is.null(components)) components <- rep.int(1L, n_pars)
  components <- .align_sampling_index(components, par_names,
                                      "components", remap = TRUE)
  if (is.null(shared_ll_idx)) shared_ll_idx <- components
  shared_ll_idx <- .align_sampling_index(shared_ll_idx, par_names,
                                         "shared_ll_idx", remap = FALSE)
  tune$components <- components
  tune$shared_ll_idx <- shared_ll_idx

  pm_settings <- attr(pmwgs$samples, "pm_settings")
  if(is.null(pm_settings)) pm_settings <- lapply(1:pmwgs$n_subjects, function(x) return(vector("list", length(unique(tune$components)))))
  tune <- check_tune_settings(tune, n_pars, stage, particles)
  pm_settings <- lapply(pm_settings, FUN = check_sampling_settings,  stage = stage, n_pars = n_pars, particles)

  eff_mu <- pmwgs$eff_mu
  eff_var <- pmwgs$eff_var
  chains_var <- pmwgs$chains_var
  chains_mu <- pmwgs$chains_mu
  if(is.null(eff_mu)) eff_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(chains_mu)) chains_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(eff_var)) eff_var <- vector("list", pmwgs$n_subjects)
  if(is.null(chains_var)) chains_var <- vector("list", pmwgs$n_subjects)
  if (verboseProgress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }

  data <- pmwgs$data
  nuisance <- pmwgs$nuisance

  marginal_idx_blk <- .marginal_par_idx(par_names, marginalise = pmwgs$marginalise)
  if (length(marginal_idx_blk) != length(par_names)) {
    marginal_idx_blk <- rep(FALSE, length(par_names))
  }
  idx_list_blk <- .component_idx_list(tune$components, marginal_idx_blk)
  chol_caches <- lapply(seq_len(pmwgs$n_subjects), function(s) {
    build_subject_chol_cache(chains_var[[s]], eff_var[[s]], idx_list_blk)
  })

  pool_budget <- .particle_core_budget(pmwgs$n_subjects, n_cores = n_cores,
                                       r_cores = r_cores,
                                       total_cores = n_cores)
  wpool_ctx <- list(
    data = data, model = .emc_wpool_slim_model(pmwgs$model),
    stage = stage, type = pmwgs$type,
    tune = tune, marginalise = pmwgs$marginalise,
    r_cores = pool_budget$likelihood, n_pars = pmwgs$n_pars,
    eff_mu = eff_mu, eff_var = eff_var, chains_mu = chains_mu,
    chains_var = chains_var, chol_caches = chol_caches
  )
  sampler_profile <- .emc_profile_enabled()
  startup_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
  n_workers <- min(pool_budget$subject, pmwgs$n_subjects)
  wpool <- .emc_wpool_start(n_workers, wpool_ctx)
  if (!is.null(wpool)) {
    on.exit(.emc_wpool_stop(wpool), add = TRUE)
  } else if (n_workers <= 1L) {
    wpool <- list(n = 1L, dir = NULL, jobs = list(), wcs = list(),
                  rcs = list(), alive = FALSE)
  }
  if (!is.null(wpool)) {
    wpool_streams <- .emc_subject_streams(pmwgs$n_subjects)
    wpool_cost <- .emc_subject_cost(data)
    wpool_part <- .emc_lpt_partition(wpool_cost, wpool$n)
  }

  llpool <- NULL
  if (pmwgs$n_subjects == 1L && pool_budget$likelihood > 1L &&
      isTRUE(getOption("emc2.ll_pool", TRUE))) {
    llpool <- .emc_wpool_start(pool_budget$likelihood, wpool_ctx)
    if (!is.null(llpool)) {
      outer_ll <- list(pool = .emc_pool_state$ll_pool,
                       s = .emc_pool_state$ll_subject,
                       route = .emc_pool_state$ll_route)
      .emc_ll_pool_set(llpool, 1L, reset = is.null(outer_ll$route))
      on.exit({
        .emc_wpool_stop(llpool)
        .emc_pool_state$ll_pool <- outer_ll$pool
        .emc_pool_state$ll_subject <- outer_ll$s
        .emc_pool_state$ll_route <- outer_ll$route
      }, add = TRUE)
    }
  }

  startup_elapsed <- if (sampler_profile) {
    proc.time()[["elapsed"]] - startup_started
  } else NA_real_

  start_iter <- pmwgs$samples$idx
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  if(any(nuisance)) pmwgs$sampler_nuis$samples$idx <- pmwgs$samples$idx

  recycle_default <- .emc_wpool_recycle_default(wpool)
  recycle_every <- as.integer(getOption("emc2.worker_recycle", recycle_default))
  if (length(recycle_every) != 1L || is.na(recycle_every)) {
    recycle_every <- recycle_default
  }

  profile_rows <- if (sampler_profile) vector("list", iter) else NULL
  gibbs_rejects <- .emc_reject_counts("gibbs")
  stage_rejects <- list(
    particle = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                               .EMC_REJECT_CLASSES),
    gibbs = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                            .EMC_REJECT_CLASSES))
  mem_every <- max(1L, as.integer(
    getOption("emc2.sampler_profile_memory_every", 10L)))

  for (i in 1:iter) {
    iteration_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    gibbs_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    pars_attempt <- tryCatch(
      gibbs_step(pmwgs, pmwgs$samples$alpha[!nuisance,,j-1], pmwgs$type),
      error = identity
    )
    if (inherits(pars_attempt, c("error", "try-error"))) {
      cls <- .emc_reject_record(pars_attempt, "gibbs")
      if (identical(.emc_failure_action(cls), "abort")) {
        .emc_failure_abort(cls, "gibbs", pars_attempt)
      }
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
        cls <- .emc_reject_record(pars_nuis_attempt, "gibbs")
        if (identical(.emc_failure_action(cls), "abort")) {
          .emc_failure_abort(cls, "gibbs", pars_nuis_attempt)
        }
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
    gibbs_elapsed <- if (sampler_profile) {
      proc.time()[["elapsed"]] - gibbs_started
    } else NA_real_
    cache_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    group_var_it <- tryCatch(
      .apply_marginal_group_var(get_group_level(pars_comb, 1L, pmwgs$type)$var,
                                marginal_idx_blk, pmwgs$marginalise),
      error = function(e) NULL
    )
    group_chol_it <- build_group_chol_cache(group_var_it, idx_list_blk,
                                          marginal_idx_blk,
                                          include_full = length(unique(tune$components)) > 1L)
    cache_elapsed <- if (sampler_profile) {
      proc.time()[["elapsed"]] - cache_started
    } else NA_real_
    recycle_elapsed <- grow_elapsed <- if (sampler_profile) 0 else NA_real_
    pool_profile <- NULL
    if (!is.null(llpool)) {
      llpool <- .emc_wpool_recycle(llpool, i, recycle_every, wpool_ctx)
      llpool <- .emc_wpool_grow(llpool, .emc_cores_now(core_ctl, n_cores),
                                wpool_ctx)
      .emc_ll_pool_set(llpool, 1L, reset = FALSE)
    }
    if (!is.null(wpool)) {
      recycle_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
      wpool <- .emc_wpool_recycle(wpool, i, recycle_every, wpool_ctx)
      if (sampler_profile) {
        recycle_elapsed <- proc.time()[["elapsed"]] - recycle_started
      }
      chain_cores <- .emc_cores_now(core_ctl, n_cores)
      target_workers <- min(
        pmwgs$n_subjects,
        max(1L, chain_cores %/% max(1L, pool_budget$likelihood))
      )
      n_before <- wpool$n
      grow_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
      wpool <- .emc_wpool_grow(wpool, target_workers, wpool_ctx)
      if (sampler_profile) {
        grow_elapsed <- proc.time()[["elapsed"]] - grow_started
      }
      if (wpool$n != n_before) {
        wpool_part <- .emc_lpt_partition(wpool_cost, wpool$n)
      }
      it <- .emc_wpool_iter(wpool, wpool_ctx, wpool_part, pars_comb,
                            group_chol_it, pm_settings,
                            pmwgs$samples$subj_ll[, j - 1], wpool_streams)
      proposals <- it$props
      pm_settings <- it$pm_settings
      wpool_streams <- it$seeds
      wpool$alive <- it$alive
      pool_profile <- it$profile
      if (!is.null(it$rejects)) {
        stage_rejects$particle <- stage_rejects$particle + it$rejects
      }
      if (any(it$times > 0)) wpool_cost <- pmax(it$times, min(it$times[it$times > 0]))
      wpool_part <- .emc_lpt_partition(wpool_cost, wpool$n)
    } else {
    core_budget <- .particle_core_budget(pmwgs$n_subjects, n_cores = n_cores,
                                         r_cores = r_cores,
                                         total_cores = n_cores)
    proposals <- parallel::mcmapply(safe_new_particle, 1:pmwgs$n_subjects, data, pm_settings, eff_mu, eff_var,
                                    chains_mu, chains_var, pmwgs$samples$subj_ll[,j-1],
                                    MoreArgs = list(parameters = pars_comb,
                                                    model = pmwgs$model,
                                                    stage = stage,
                                                    type = pmwgs$type,
                                                    tune = tune,
                                                    marginalise = pmwgs$marginalise,
                                                    r_cores = core_budget$likelihood,
                                                    group_chol = group_chol_it),
                                    chol_cache = chol_caches,
                                    mc.cores = core_budget$subject)
    pm_settings <- proposals[3,]
    proposals <- array(unlist(proposals[1:2,]), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))
    }

    fill_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, j = j, n_pars = pmwgs$n_pars, type = pmwgs$type)
    fill_elapsed <- if (sampler_profile) {
      proc.time()[["elapsed"]] - fill_started
    } else NA_real_
    gibbs_delta <- .emc_reject_delta(gibbs_rejects, "gibbs")
    gibbs_rejects <- .emc_reject_counts("gibbs")
    stage_rejects$gibbs <- stage_rejects$gibbs + gibbs_delta
    if (sampler_profile) {
      iteration_total <- proc.time()[["elapsed"]] - iteration_started
      pp <- if (is.null(pool_profile)) list() else pool_profile
      n_part <- vapply(pm_settings, function(x) {
        if (is.null(x$n_particles)) NA_real_ else as.numeric(x$n_particles)
      }, numeric(1))
      mem <- if (.emc_profile_expensive() &&
                 (i == 1L || i == iter || i %% mem_every == 0L)) {
        .emc_profile_memory(if (is.null(wpool)) NULL else c(
          vapply(wpool$jobs, .emc_wpool_job_pid, integer(1)),
          .emc_wpool_job_pid(wpool$template$job)))
      } else list(bytes = NA_real_, method = NA_character_)
      route <- .emc_ll_route(wpool, llpool, pool_budget,
                             if (exists("core_budget", inherits = FALSE)) core_budget else NULL,
                             .emc_pool_state$ll_route)
      profile_rows[[i]] <- do.call(.emc_profile_row, c(list(
        iteration = j,
        stage = stage,
        subjects = pmwgs$n_subjects,
        workers = if (is.null(wpool)) 0L else wpool$n,
        route = route,
        workers_active = pp$workers_active,
        particles_max = if (all(is.na(n_part))) NA_real_ else max(n_part, na.rm = TRUE),
        particles_sum = if (all(is.na(n_part))) NA_real_ else sum(n_part, na.rm = TRUE),
        total = iteration_total,
        gibbs = gibbs_elapsed,
        group_cache = cache_elapsed,
        recycle = recycle_elapsed,
        grow = grow_elapsed,
        startup = if (i == 1L) startup_elapsed else NA_real_,
        particle = pp$elapsed,
        fill = fill_elapsed,
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
        degraded = pp$degraded),
        .emc_reject_row_fields(pp$rejects, "particle"),
        .emc_reject_row_fields(gibbs_delta, "gibbs")))
    }
  }
  .emc_failure_report_stage(stage_rejects, stage)
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  if (sampler_profile) {
    attr(pmwgs$samples, "sampler_profile") <- .emc_profile_bind(profile_rows)
  }
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
    if (length(d) == 2 && d[2] >= j) {
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
                               parameters = NULL, model = NULL, stage,
                               type, tune, marginalise = NULL, r_cores = 1,
                               chol_cache = NULL, group_chol = NULL,
                               current_alpha = NULL, population_mu = NULL,
                               population_var = NULL) {
  attempt <- tryCatch(
    new_particle(
      s = s, data = data, pm_settings = pm_settings,
      eff_mu = eff_mu, eff_var = eff_var,
      chains_mu = chains_mu, chains_var = chains_var,
      prev_ll = prev_ll, parameters = parameters, model = model,
      stage = stage, type = type, tune = tune, marginalise = marginalise,
      r_cores = r_cores, chol_cache = chol_cache, group_chol = group_chol,
      current_alpha = current_alpha, population_mu = population_mu,
      population_var = population_var
    ),
    error = identity
  )
  if (inherits(attempt, c("error", "try-error"))) {
    cls <- .emc_reject_record(attempt, "particle")
    if (identical(.emc_failure_action(cls), "abort")) {
      .emc_failure_abort(cls, "particle", attempt)
    }
    old <- if (is.null(current_alpha)) parameters$alpha[, s] else current_alpha
    return(reject_particle(old, prev_ll, pm_settings))
  }
  attempt
}

reject_particle <- function(subj_mu, prev_ll, pm_settings) {
  list(proposal = subj_mu, ll = prev_ll, pm_settings = pm_settings)
}



.chol_factor <- function(Sigma, p) {
  if (p <= 0) return(NULL)
  base_R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(base_R)) {
    base_R <- tryCatch(chol(Sigma + diag(1e-6, p)), error = function(e) NULL)
  }
  if (is.null(base_R)) {
    base_R <- tryCatch(chol(Sigma + diag(1e-4, p)), error = function(e) NULL)
  }
  if (is.null(base_R)) {
    cov_diag <- diag(Sigma)
    cov_diag[is.na(cov_diag) | cov_diag <= 0] <- 1e-4
    base_R <- diag(sqrt(cov_diag), p)
  }
  rooti <- backsolve(base_R, diag(p))
  list(R = base_R, rooti = rooti, sum_log_diag = sum(log(diag(rooti))), p = p)
}

.marginal_par_idx <- function(par_names, marginalise) {
  if (is.null(marginalise) || is.null(par_names)) {
    return(rep(FALSE, length(par_names)))
  }
  par_names %in% marginalise$param
}

.apply_marginal_group_var <- function(group_var, marginal_idx, marginalise) {
  if (!any(marginal_idx)) return(group_var)
  group_var[marginal_idx, ] <- 0
  group_var[, marginal_idx] <- 0
  group_var[marginal_idx, marginal_idx] <- marginalise$sigma^2
  group_var
}

.component_idx_list <- function(components, marginal_idx) {
  unq <- unique(components)
  out <- vector("list", max(unq))
  for (i in unq) out[[i]] <- (components == i) & !marginal_idx
  out
}

build_subject_chol_cache <- function(chains_var, eff_var, idx_list) {
  if (is.null(chains_var) && is.null(eff_var)) return(NULL)
  facs <- function(S) {
    if (is.null(S)) return(NULL)
    lapply(idx_list, function(idx) {
      if (is.null(idx)) return(NULL)
      .chol_factor(S[idx, idx, drop = FALSE], sum(idx))
    })
  }
  list(idx_list = idx_list,
       chains_ref = chains_var, chains = facs(chains_var),
       eff_ref = eff_var, eff = facs(eff_var))
}

build_group_chol_cache <- function(group_var, idx_list, marginal_idx = NULL,
                                   include_full = TRUE) {
  if (is.null(group_var)) return(NULL)
  full <- NULL
  if (isTRUE(include_full) && !is.null(marginal_idx) && any(!marginal_idx)) {
    keep <- !marginal_idx
    full <- fast_dmvnorm_factor(group_var[keep, keep, drop = FALSE])
    full$keep <- keep
  }
  list(idx_list = idx_list, ref = group_var, marginal_idx = marginal_idx,
       include_full = isTRUE(include_full),
       full = full,
       f = lapply(idx_list, function(idx) {
         if (is.null(idx)) return(NULL)
         .chol_factor(group_var[idx, idx, drop = FALSE], sum(idx))
       }))
}

new_particle <- function (s, data, pm_settings, eff_mu = NULL,
                          eff_var = NULL, chains_mu = NULL,
                          chains_var = NULL, prev_ll,
                          parameters, model = NULL, stage,
                          type, tune, marginalise = NULL, r_cores = 1,
                          chol_cache = NULL, group_chol = NULL,
                          current_alpha = NULL, population_mu = NULL,
                          population_var = NULL)
{
  direct_state <- !is.null(current_alpha) && !is.null(population_mu) &&
    !is.null(population_var)
  group_pars <- if (direct_state) {
    list(mu = population_mu, var = population_var)
  } else {
    get_group_level(parameters, s, type)
  }
  unq_components <- unique(tune$components)
  proposal_out <- numeric(length(group_pars$mu))
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- if (direct_state) current_alpha else parameters$alpha[,s]
  marg_nodes <- NULL
  marg_terms_row <- NULL
  marginal_idx <- rep(FALSE, length(subj_mu))
  if (!is.null(marginalise)) {
    if (is.null(names(subj_mu))) {
      alpha_names <- if (direct_state) names(current_alpha) else rownames(parameters$alpha)
      names(subj_mu) <- alpha_names[seq_along(subj_mu)]
    }
    marginal_idx <- .marginal_par_idx(names(subj_mu), marginalise)
    group_mu[marginal_idx] <- marginalise$mu
    group_var <- .apply_marginal_group_var(group_var, marginal_idx, marginalise)
  }
  out_lls <- numeric(length(unq_components))
  particle_multiplier <- 1
  if(stage == "preburn"){
    Mus <- list(group_mu, subj_mu)
    Sigmas <- list(group_var, group_var)
    sig_tags <- c("group", "group")
    particle_multiplier <- 2
  } else if(stage == "burn"){
    Mus <- list(group_mu, subj_mu, subj_mu)
    Sigmas <- list(group_var, group_var, chains_var)
    sig_tags <- c("group", "group", "chains")
  } else if(stage == "adapt"){
    Mus <- list(group_mu, subj_mu, chains_mu)
    Sigmas <- list(group_var, chains_var, chains_var)
    sig_tags <- c("group", "chains", "chains")
  } else{
    Mus <- list(group_mu, subj_mu, chains_mu, eff_mu)
    Sigmas <- list(group_var, chains_var, chains_var, eff_var)
    sig_tags <- c("group", "chains", "chains", "eff")
  }
  n_proposals <- length(Mus)
  cache_ok <- function(cache, ref_name, ref, idx) {
    !is.null(cache) && !is.null(cache[[ref_name]]) &&
      identical(cache[[ref_name]], ref) &&
      identical(cache$idx_list[[i]], idx)
  }
  for(i in unq_components){
    epsilons <- c(1, pm_settings[[i]]$epsilon)
    idx_full <- tune$components == i
    idx <- idx_full & !marginal_idx
    p_idx <- sum(idx)

    Rs <- vector("list", n_proposals)
    rootis <- vector("list", n_proposals)
    log_consts <- numeric(n_proposals)
    if (p_idx > 0) {
      bases <- list()
      for (tag in unique(sig_tags)) {
        k_first <- which(sig_tags == tag)[1]
        base <- switch(tag,
          group  = if (cache_ok(group_chol, "ref", group_var, idx)) group_chol$f[[i]],
          chains = if (cache_ok(chol_cache, "chains_ref", chains_var, idx)) chol_cache$chains[[i]],
          eff    = if (cache_ok(chol_cache, "eff_ref", eff_var, idx)) chol_cache$eff[[i]],
          NULL)
        if (is.null(base)) {
          base <- .chol_factor(Sigmas[[k_first]][idx, idx, drop = FALSE], p_idx)
        }
        bases[[tag]] <- base
      }
      const <- -0.5 * p_idx * log(2 * pi)
      for (k in seq_len(n_proposals)) {
        base <- bases[[sig_tags[k]]]
        Rs[[k]] <- base$R * epsilons[k]
        rootis[[k]] <- base$rooti / epsilons[k]
        log_consts[[k]] <- base$sum_log_diag - p_idx * log(epsilons[k]) + const
      }
    }

    particle_numbers <- numbers_from_proportion(pm_settings[[i]]$mix, pm_settings[[i]]$n_particles*particle_multiplier)
    proposals <- vector("list", n_proposals +1)
    proposals[[1]] <- matrix(subj_mu[idx], nrow = 1L)
    for(j in 1:n_proposals){
      if (p_idx > 0) {
        proposals[[j + 1]] <- particle_draws(
          particle_numbers[j], Mus[[j]][idx], covar = NULL, R = Rs[[j]])
      } else {
        proposals[[j + 1]] <- matrix(numeric(0), nrow = particle_numbers[j], ncol = 0L)
      }
    }
    proposals <- do.call(rbind, proposals)

    if(any(!idx)){
      proposals_full <- matrix(rep(subj_mu, each = nrow(proposals)),
                                nrow = nrow(proposals),
                                ncol = length(subj_mu))
      proposals_full[, idx] <- proposals
      colnames(proposals_full) <- names(subj_mu)
      proposals <- proposals_full
    } else{
      colnames(proposals) <- names(subj_mu)
    }

    shared_idx <- tune$shared_ll_idx[idx_full][1]
    is_shared <- shared_idx == tune$shared_ll_idx

    if (!is.null(marginalise)) {
      warm_state <- pm_settings[[i]]$marg_warm
      warm_arg <- if (!is.null(warm_state) && warm_state$backoff <= 0L) {
        c(mode = warm_state$mode, sd = warm_state$sd)
      } else NULL
      pred_state <- pm_settings[[i]]$marg_pred
      use_pred <- is.null(pred_state) || pred_state$backoff <= 0L
      marg_spec <- marginalise
      if (!use_pred) marg_spec$predict_mode <- FALSE
      marg_grid <- compute_marginal_grid(proposals[, is_shared, drop = FALSE],
                                         data, model, marg_spec,
                                         r_cores = r_cores, warm = warm_arg)
      lw <- marginal_ll_from_grid(marg_grid)
    } else if(tune$components[length(tune$components)] > 1){
      lw <- calc_ll_pooled(proposals[,is_shared], dadm = data, model,
                           component = shared_idx, r_cores = r_cores, s = s,
                           varying = idx[is_shared])
    } else{
      lw <- calc_ll_pooled(proposals[,is_shared], dadm = data, model,
                           r_cores = r_cores, s = s, varying = idx[is_shared])
    }
    lw_total <- lw + prev_ll - lw[1]
    lp <- if (p_idx > 0) {
      fast_dmvnorm_rooti(x = proposals[,idx,drop=FALSE], mean = group_mu[idx],
                         rooti = rootis[[1]], log_const = log_consts[[1]])
    } else rep(0, nrow(proposals))
    if(length(unq_components) > 1){
      prior_density <- if (any(!marginal_idx)) {
        gf <- group_chol$full
        if (!is.null(gf) && !is.null(group_chol$ref) &&
            identical(group_chol$ref, group_var) &&
            identical(gf$keep, !marginal_idx)) {
          if (isTRUE(gf$ok)) {
            fast_dmvnorm_rooti(x = proposals[, !marginal_idx, drop = FALSE],
                               mean = group_mu[!marginal_idx],
                               rooti = gf$rooti, log_const = gf$log_const)
          } else {
            rep(-Inf, nrow(proposals))
          }
        } else {
          fast_dmvnorm(x = proposals[,!marginal_idx,drop=FALSE],
                       mean = group_mu[!marginal_idx],
                       sigma = group_var[!marginal_idx,!marginal_idx,drop=FALSE])
        }
      } else rep(0, nrow(proposals))
    } else{
      prior_density <- lp
    }
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
    max_log <- log_mix_comps[, 1L]
    for (k in seq_len(n_proposals)[-1L]) {
      max_log <- pmax(max_log, log_mix_comps[, k])
    }
    lm <- max_log + log(rowSums(exp(log_mix_comps - max_log)))
    infnt_idx <- is.infinite(lm) | is.na(lm)
    if (any(infnt_idx)) {
      fin_vals <- lm[!infnt_idx]
      lm[infnt_idx] <- if (length(fin_vals) > 0) min(fin_vals) else -1e10
    }
    l <- lw_total + prior_density - lm
    weights <- exp(l - max(l))
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
      pred_next_state <- marginal_warm_backoff(
        pred_state, attempted = use_pred && !isTRUE(marg_grid$warm_used == 1),
        used = isTRUE(marg_grid$pred_used == 1))
    }
    pm_settings[[i]] <- update_pm_settings(pm_settings[[i]], idx_ll, weights, particle_numbers, tune, sum(idx))
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
  pm_settings$iter <- pm_settings$iter + 1
  if (pm_settings$iter > tune$n0) {
    pm_settings$proposal_counts <- pm_settings$proposal_counts + particle_numbers

    old_weight <- weights[1]

    offset <- 2
    for (j in seq_along(particle_numbers)) {
      n_j <- particle_numbers[j]
      if (n_j > 0) {
        draws_j <- weights[offset:(offset + n_j - 1)]
        better_j <- sum(draws_j > old_weight)
        pm_settings$acc_counts[j] <- pm_settings$acc_counts[j] + better_j
        offset <- offset + n_j
      }
    }

    acc_rates <- ifelse(pm_settings$proposal_counts > 0, pm_settings$acc_counts / pm_settings$proposal_counts, 0)

    clamp_min <- ifelse(length(pm_settings$mix) == 2, .1, ifelse(length(pm_settings$mix) == 3, .4, .6))

    new_epsilon <- update_epsilon_continuous(
      epsilon   = pm_settings$epsilon,
      acceptance = acc_rates[-1],
      target     = tune$p_accept,
      iter       = pm_settings$iter,
      d          = n_pars,
      alphaStar  = tune$alphaStar,
      damp       = 100,
      clamp      = c(clamp_min, 5)
    )
    pm_settings$epsilon <- new_epsilon


    if(length(pm_settings$mix) > 2){
      eps_val <- 1e-12

      performance <- (acc_rates + eps_val) / pm_settings$mix
      performance <- performance / sum(performance)

      adj_factor <- c(mean(tune$p_accept), tune$p_accept)
      performance <- performance / adj_factor

      performance <- performance / sum(performance)

      new_mix <- (1 - tune$mix_adapt) * pm_settings$mix + tune$mix_adapt * performance

      new_mix <- pmax(new_mix, 0.02)
      new_mix <- new_mix / sum(new_mix)
      pm_settings$mix <- new_mix
    }


    if (length(pm_settings$mix) > 3 && pm_settings$gd_good) {
      ess <- sum(weights)^2 / sum(weights^2)
      desired_ess <- tune$target_ESS
      scale_factor <- (desired_ess / ess)^tune$ESS_scale
      new_num_particles <- round(pm_settings$n_particles * scale_factor)
      pm_settings$n_particles <- max(25, min(tune$max_particles, new_num_particles))
    }
  }

  pm_settings$weight_ess <- if (length(weights)) {
    sum(weights)^2 / sum(weights^2)
  } else {
    NA_real_
  }
  return(pm_settings)
}



update_epsilon_continuous <- function(
    epsilon,
    acceptance,
    target,
    iter,
    d,
    alphaStar,
    damp = 100,
    clamp = c(0.6, 4)
) {
  log_eps <- log(epsilon)
  c_term <- (1 - 1/d)*sqrt(2*pi)*exp(alphaStar^2/2)/(2*alphaStar) + 1/(d*target*(1-target))
  step_size <- c_term / max(damp, iter)
  diff_accept <- acceptance - target
  log_eps_new <- log_eps + step_size * diff_accept
  eps_new <- exp(log_eps_new)
  eps_new <- pmin(clamp[2], pmax(eps_new, clamp[1]))

  return(eps_new)
}

numbers_from_proportion <- function(mix_proportion, n_particles = 1000) {
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
      return(X + rep(mu, each = n))
    } else {
      return(mvtnorm::rmvnorm(n, mu, covar))
    }
  }
}

extend_sampler <- function(sampler, n_samples, stage) {
  n_iter <- length(sampler$samples$stage)
  sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
  if(any(sampler$nuisance)) sampler$sampler_nuis$samples <- rapply(sampler$sampler_nuis$samples, f = function(x) extend_obj(x, n_samples, n_iter), how = "replace")
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples, n_iter), how = "replace")
  return(sampler)
}

is_iteration_array <- function(obj, n_iter){
  d <- dim(obj)
  if(is.null(d) || length(d) < 2L) return(FALSE)
  isTRUE(d[length(d)] == n_iter)
}

extend_obj <- function(obj, n_extend, n_iter = NULL){
  old_dim <- dim(obj)
  n_dimensions <- length(old_dim)
  if(is.null(old_dim) | n_dimensions == 1) return(obj)
  if(!is.null(n_iter) && !is_iteration_array(obj, n_iter)) return(obj)
  new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
  extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
  emc_copy_sample_prefix(extended, as.numeric(obj))
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
  emc_set_last_slice(samples$theta_mu, j, as.numeric(group_level$tmu))
  emc_set_last_slice(samples$theta_var, j, as.numeric(group_level$tvar))
  if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, j, n_pars)
  return(samples)
}



fill_samples_RE <- function(samples, proposals, j = 1, n_pars, ...){
  if(!is.null(proposals)){
    emc_set_particle_slice(samples$alpha, samples$subj_ll,
                           as.matrix(proposals), j, n_pars)
    samples$idx <- j
  }
  return(samples)
}


set_p_accept <- function(stage, search_width){
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
  if (is.null(epsilon)) {
    if (n_pars > 15) {
      epsilon <- .5
    } else if (n_pars > 10) {
      epsilon <- .6
    } else {
      epsilon <- .7
    }
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

calc_ll_pooled <- function(proposals, dadm, model, component = NULL, r_cores = 1,
                           s = NULL, varying = NULL){
  have_pool <- !is.null(s) && !is.null(.emc_pool_state$ll_subject) &&
    identical(s, .emc_pool_state$ll_subject) &&
    isTRUE(.emc_pool_state$ll_pool$alive)
  if (!have_pool) {
    route_started <- if (!is.null(.emc_pool_state$ll_route))
      proc.time()[["elapsed"]] else NA_real_
    out <- calc_ll_manager(proposals, dadm = dadm, model = model,
                           component = component, r_cores = r_cores,
                           varying = varying)
    if (is.finite(route_started)) {
      .emc_ll_route_record(FALSE,
                           proc.time()[["elapsed"]] - route_started,
                           if (is.matrix(proposals)) nrow(proposals) else 1L,
                           route = if (r_cores > 1L) "nested_per_call" else "serial")
    }
    return(out)
  }
  use_pool <- .emc_ll_route_use_pool()
  dynamic <- use_pool && .emc_ll_route_use_dynamic()
  started <- proc.time()[["elapsed"]]
  out <- if (use_pool) {
    .emc_wpool_ll(proposals, s, component, dynamic, varying)
  } else NULL
  pooled <- !is.null(out)
  if (!pooled) {
    out <- calc_ll_manager(proposals, dadm = dadm, model = model,
                           component = component, r_cores = 1L,
                           varying = varying)
  }
  el <- proc.time()[["elapsed"]] - started
  actual_route <- if (pooled) "persistent_pool" else "serial"
  .emc_ll_route_record(pooled, el,
                       if (is.matrix(proposals)) nrow(proposals) else 1L,
                       dynamic = pooled && isTRUE(.emc_pool_state$ll_dynamic),
                       route = actual_route)
  trace_file <- Sys.getenv("EMC2_LL_TRACE")
  if (nzchar(trace_file)) {
    st <- .emc_pool_state$ll_route
    cat(sprintf(paste("[ll] pid=%d n=%d mode=%s n_serial=%d n_pool=%d pooled=%s",
                      "split=%s n_static=%d n_dyn=%d dyn=%s %.3f s\n"),
                Sys.getpid(),
                if (is.matrix(proposals)) nrow(proposals) else 1L,
                st$mode, length(st$t_serial), length(st$t_pool), pooled,
                st$split, length(st$t_static), length(st$t_dynamic),
                pooled && isTRUE(.emc_pool_state$ll_dynamic), el),
        file = trace_file, append = TRUE)
  }
  out
}

calc_ll_manager <- function(proposals, dadm, model, component = NULL, r_cores = 1,
                            marginalise = NULL, varying = NULL){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component,
                                r_cores = r_cores, marginalise = marginalise,
                                varying = varying)
  } else{
    model <- model()
    dadm <- .cache_ll_data_attrs(dadm)
    if(is.null(model$c_name)){
      if (!is.null(marginalise)) {
        stop("marginalise requires a registered race-model likelihood")
      }
      lls <- unlist(
        auto_mclapply(1:nrow(proposals),
          function(i) calc_ll_R(proposals[i,], model=model, dadm = dadm),
         mc.cores=r_cores))
    } else {
      set_stop_method_from_model(model)
      p_types <- names(model$p_types)
      designs <- .oo_expanded_designs(dadm, expand = FALSE)
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      if (r_cores <= 1L || nrow(proposals) <= r_cores) {
        lls <- calc_ll_oo(proposals, dadm, constants = constants, designs = designs,
                          type = model$c_name, bounds = model$bound,
                          transforms = model$transform, pretransforms = model$pre_transform,
                          p_types = p_types, min_ll = log(1e-10), trend = model$trend,
                          marginalise = marginalise, varying = varying)
      } else {
        idx <- .split_work_indices(nrow(proposals), r_cores)
        lls <- unlist(auto_mclapply(1:r_cores,function(i) {
          calc_ll_oo(proposals[idx==i,,drop=FALSE], dadm, constants = constants,
                     designs = designs, type = model$c_name, bounds = model$bound,
                     transforms = model$transform, pretransforms = model$pre_transform,
                     p_types = p_types, min_ll = log(1e-10),
                     trend = model$trend, marginalise = marginalise,
                     varying = varying)
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
  subj_mu_out[is_nuisance,] <- matrix(c(tmu_nuis), nrow = sum(is_nuisance),
                                      ncol = ncol(subj_mu))
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

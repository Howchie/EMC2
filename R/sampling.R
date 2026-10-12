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
    if (suppressWarnings(!is.na(as.integer(marginalise[2])))) as.integer(marginalise[2]) else 5L
  } else 5L

  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes < 3L) {
    stop("marginalise$n_nodes must be an integer >= 3")
  }
  list(param = param, mu = unname(mu), sigma = unname(sigma), n_nodes = n_nodes)
}

compute_marginal_grid <- function(proposals, data, model, marginalise,
                                  r_cores = 1) {
  model_spec <- if (is.function(model)) model() else model
  set_stop_method_from_model(model_spec)
  data <- .cache_ll_data_attrs(data)
  constants <- attr(data, "constants")
  if (is.null(constants)) constants <- NA
  designs <- .oo_expanded_designs(data, expand = FALSE)
  one <- function(props, spec) {
    calc_ll_oo_marginal_nodes(
      props, data, constants = constants, designs = designs,
      type = model_spec$c_name, bounds = model_spec$bound,
      transforms = model_spec$transform, pretransforms = model_spec$pre_transform,
      p_types = names(model_spec$p_types), min_ll = log(1e-10),
      marginalise = spec, trend = model_spec$trend
    )
  }
  if (r_cores <= 1 || nrow(proposals) <= r_cores) return(one(proposals, marginalise))
  idx <- .split_work_indices(nrow(proposals), r_cores)
  parts <- auto_mclapply(1:r_cores, function(i) {
    spec <- marginalise
    for (f in c("mu", "centre", "probe_sd", "log_offset"))
      if (length(spec[[f]]) > 1L) spec[[f]] <- spec[[f]][idx == i]
    one(proposals[idx == i, , drop = FALSE], spec)
  }, mc.cores = r_cores)
  list(nodes = do.call(rbind, lapply(parts, `[[`, "nodes")),
       log_terms = do.call(rbind, lapply(parts, `[[`, "log_terms")),
       ll = unlist(lapply(parts, `[[`, "ll")),
       mode = unlist(lapply(parts, `[[`, "mode")),
       sd = unlist(lapply(parts, `[[`, "sd")),
       rounds = unlist(lapply(parts, `[[`, "rounds")),
       gl_used = sum(unlist(lapply(parts, `[[`, "gl_used"))),
       screened = sum(unlist(lapply(parts, `[[`, "screened"))))
}

.split_work_indices <- function(n, n_workers) {
  n <- as.integer(n)
  n_workers <- min(as.integer(n_workers), n)
  if (n <= 0L || n_workers <= 0L) return(integer(0))
  q <- n %/% n_workers
  rem <- n %% n_workers
  rep.int(seq_len(n_workers), q + (seq_len(n_workers) <= rem))
}

# Where a batch's peaks sat: a regression of each particle's located peak and
# log scale on the rest of its vector, which starts the next batch's rules.
# Coordinates the batch held fixed keep the subject's posterior slope
# (`predict`, a .marginal_conditional of the chains' moments) when there is
# one. NULL when the batch cannot identify the regression.
.marginal_state <- function(grid, proposals, cond, predict = NULL) {
  ok <- is.finite(grid$mode) & is.finite(grid$sd) & grid$sd > 0
  X <- proposals[ok, cond$o, drop = FALSE]
  if (!nrow(X)) return(NULL)
  x0 <- colMeans(X)
  D <- sweep(X, 2L, x0)
  vary <- colSums(D != 0) > 0
  if (nrow(X) <= sum(vary) + 1L) return(NULL)
  Z <- cbind(1, D[, vary, drop = FALSE])
  fm <- stats::lm.wfit(Z, grid$mode[ok], w = 1 / grid$sd[ok]^2)$coefficients
  fs <- stats::lm.fit(Z, log(grid$sd[ok]))$coefficients
  if (anyNA(fm) || anyNA(fs)) return(NULL)
  beta <- if (!is.null(predict)) predict$b else numeric(length(cond$o))
  lbeta <- numeric(length(cond$o))
  beta[vary] <- fm[-1L]
  lbeta[vary] <- fs[-1L]
  list(x0 = x0, m0 = fm[[1L]], beta = beta, l0 = fs[[1L]], lbeta = lbeta)
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

# The marginalised parameter's prior given the rest of a subject's vector:
# under the group normal, x_k | x_o ~ N(a + x_o b, sd^2).
.marginal_conditional <- function(group_mu, group_var, marginal_idx) {
  k <- which(marginal_idx)
  o <- which(!marginal_idx)
  group_mu <- as.numeric(group_mu)
  b <- numeric(length(o))
  v <- group_var[k, k]
  s_ok <- group_var[o, k]
  if (length(o) && any(s_ok != 0)) {
    b <- solve(group_var[o, o, drop = FALSE], s_ok)
    v <- v - sum(s_ok * b)
  }
  if (!is.finite(v) || v <= 0)
    stop(errorCondition(
      "group covariance is not positive definite; marginalised prior undefined",
      class = c("emc_group_prior_failure", "emc_failure"),
      emc_class = "numerical", emc_source = "particle"))
  list(k = k, o = o, a = group_mu[k] - sum(b * group_mu[o]), b = b, sd = sqrt(v))
}

.marginal_cond_mean <- function(cond, x) {
  as.numeric(cond$a + x[, cond$o, drop = FALSE] %*% cond$b)
}

# The quadrature spec for a batch of particles: each particle's own
# conditional prior mean and the common conditional SD. Each rule starts at
# the particle's predicted peak: from the subject's last batch (`state`, a
# .marginal_state) when there is one, else from its posterior moments
# (`predict`, a .marginal_conditional of the chains' mean and covariance).
.marginal_batch_spec <- function(marginalise, cond, proposals, state = NULL,
                                 predict = NULL) {
  marginalise$mu <- .marginal_cond_mean(cond, proposals)
  marginalise$sigma <- cond$sd
  if (!is.null(state)) {
    D <- sweep(proposals[, cond$o, drop = FALSE], 2L, state$x0)
    marginalise$centre <- as.numeric(state$m0 + D %*% state$beta)
    marginalise$probe_sd <- as.numeric(exp(state$l0 + D %*% state$lbeta))
  } else if (!is.null(predict)) {
    marginalise$centre <- .marginal_cond_mean(predict, proposals)
    marginalise$probe_sd <- predict$sd
  }
  marginalise
}

# Completes subject vector x with a draw of the marginalised parameter from
# its conditional posterior given the rest of x. A node drawn from the
# quadrature's discrete posterior starts one independence Metropolis step
# towards the continuous conditional, proposing from the grid's Laplace fit.
# Returns the completed vector and its log-likelihood.
.marginal_draw <- function(x, cond, nodes, log_terms, mode, sd, data, model) {
  k <- cond$k
  m <- .marginal_cond_mean(cond, matrix(x, nrow = 1L))
  fin <- is.finite(log_terms)
  cur <- if (length(nodes) && any(fin)) {
    w <- exp(log_terms - max(log_terms[fin]))
    w[!fin] <- 0
    nodes[[sample.int(length(nodes), size = 1L, prob = w)]]
  } else stats::rnorm(1L, m, cond$sd)
  laplace <- length(mode) == 1L && length(sd) == 1L &&
    is.finite(mode) && is.finite(sd) && sd > 0
  rows <- matrix(x, nrow = if (laplace) 2L else 1L, ncol = length(x),
                 byrow = TRUE, dimnames = list(NULL, names(x)))
  rows[1L, k] <- cur
  if (laplace) rows[2L, k] <- stats::rnorm(1L, mode, sd)
  ll <- as.numeric(calc_ll_manager(rows, dadm = data, model))
  j <- 1L
  if (laplace) {
    lt <- ll + stats::dnorm(rows[, k], m, cond$sd, log = TRUE) -
      stats::dnorm(rows[, k], mode, sd, log = TRUE)
    log_a <- lt[2L] - lt[1L]
    if (is.finite(log_a) && log(stats::runif(1L)) < log_a) j <- 2L
  }
  x[k] <- rows[j, k]
  list(x = x, ll = ll[j])
}

# The particle step targets an approximate joint density
#   pi~(x_o, x_k) = Q(x_o) N(x_o | group) g(x_k | x_o),
# where Q is the quadrature's marginal likelihood and g the reconstruction
# density: the grid's Laplace fit, or the conditional group prior where the
# grid has none. Each conditional-IS step moves x_o reversibly for pi~ and
# draws x_k from g, and a second accept/reject on
#   r = pi / pi~ = L(x_o, x_k) N(x_k | x_o, group) / (Q(x_o) g(x_k | x_o))
# makes the step reversible for the exact posterior (delayed acceptance).
.marginal_g <- function(m, cond_sd, mode, sd) {
  if (length(mode) == 1L && length(sd) == 1L && is.finite(mode) &&
      is.finite(sd) && sd > 0) list(mean = mode, sd = sd)
  else list(mean = m, sd = cond_sd)
}

.marginal_log_ratio <- function(ll, t, m, cond_sd, q, g) {
  ll + stats::dnorm(t, m, cond_sd, log = TRUE) - q -
    stats::dnorm(t, g$mean, g$sd, log = TRUE)
}

# One delayed-acceptance decision for the particle selected in a step.
# `proposals` holds the step's full vectors with the reference in row 1,
# `j` the selected row, `grid` and `q` the step's quadrature and marginal
# log-likelihoods, and (t_cur, ll_cur) the current marginalised value and
# complete log-likelihood. Returns whether to move, with the new value and
# complete log-likelihood.
.marginal_delayed_accept <- function(proposals, j, grid, q, cond, t_cur,
                                     ll_cur, data, model) {
  k <- cond$k
  m <- .marginal_cond_mean(cond, proposals[c(1L, j), , drop = FALSE])
  g_cur <- .marginal_g(m[1L], cond$sd, grid$mode[1L], grid$sd[1L])
  g_new <- .marginal_g(m[2L], cond$sd, grid$mode[j], grid$sd[j])
  t_new <- stats::rnorm(1L, g_new$mean, g_new$sd)
  row <- proposals[j, , drop = FALSE]
  row[, k] <- t_new
  ll_new <- as.numeric(calc_ll_manager(row, dadm = data, model))
  log_a <- .marginal_log_ratio(ll_new, t_new, m[2L], cond$sd, q[j], g_new) -
    .marginal_log_ratio(ll_cur, t_cur, m[1L], cond$sd, q[1L], g_cur)
  list(accept = !is.na(log_a) && log(stats::runif(1L)) < log_a,
       t = t_new, ll = ll_new)
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
  if (!is.null(attr(prior, "pooled")) &&
      (length(nuisance) || length(nuisance_non_hyper)))
    stop("pooled parameters cannot be combined with nuisance samplers")
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
  pooled_prior <- attr(prior, "pooled")
  if (!is.null(pooled_prior)) {
    sampler$pooled <- .emc_pooled_state(pooled_prior)
    .emc_pooled_unsupported(sampler)
    sampler$samples$pooled <- matrix(NA_real_, length(pooled_prior$mean),
                                     length(sampler$samples$stage),
                                     dimnames = list(names(pooled_prior$mean), NULL))
  }
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
  pmwgs <- .emc_prepare_sampler_rng(pmwgs)
  if (!is.null(pmwgs$sampler_nuis$samples)) {
    pmwgs$sampler_nuis$samples <- emc_clone_sample_store(pmwgs$sampler_nuis$samples)
  }
  type <- pmwgs$type
  pooled_start <- NULL
  if (!is.null(pmwgs$pooled) && !is.null(start_mu) && !is.null(names(start_mu))) {
    pn <- .emc_pooled_names(pmwgs)
    pooled_start <- start_mu[intersect(pn, names(start_mu))]
    keep <- !names(start_mu) %in% pn
    if (!is.null(start_var) && is.matrix(start_var) && nrow(start_var) == length(start_mu))
      start_var <- start_var[keep, keep, drop = FALSE]
    start_mu <- start_mu[keep]
  }
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
  start_streams <- pmwgs$rng$subjects
  raw_proposals <- .emc_with_preserved_rng(
    auto_mclapply(X=seq_len(pmwgs$n_subjects),
                  FUN=function(s) {
                    assign(".Random.seed", start_streams[[s]],
                           envir = globalenv())
                    out <- start_proposals(s, parameters = startpoints_comb,
                                           n_particles = particles,
                                           pmwgs = pmwgs, type = type,
                                           r_cores = core_budget$likelihood)
                    list(out = out,
                         seed = get(".Random.seed", envir = globalenv()))
                  },
                  mc.cores = core_budget$subject, mc.set.seed = FALSE))
  proposals <- lapply(raw_proposals, `[[`, "out")
  pmwgs$rng$subjects <- lapply(raw_proposals, `[[`, "seed")
  pooled_draws <- if (is.null(pmwgs$pooled)) NULL else
    vapply(proposals, function(x) as.numeric(x$pooled), numeric(length(pmwgs$pooled$value)))
  proposals <- lapply(proposals, `[`, c("proposal", "ll"))
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))
  if (!is.null(pmwgs$pooled)) {
    pin <- .emc_pooled_init(pmwgs, proposals, pooled_draws, pooled_start)
    pmwgs <- pin$sampler
    proposals <- pin$proposals
  }


  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             j = 1, n_pars = pmwgs$n_pars, type = type)
  if (!is.null(pmwgs$pooled)) pmwgs$samples$pooled[, 1] <- pmwgs$pooled$value
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
  emc <- restore_custom_kernel_pointers(emc, quiet = TRUE)
  emc <- .emc_prepare_chain_rng(emc)
  dots <- add_defaults(list(...),r_cores=1)
  emc <- auto_mclapply(emc, function(sampler, ...) {
           sampler <- .emc_set_sampler_rng(sampler)
           init(sampler, ...)
         }, start_mu = start_mu, start_var = start_var,
           verbose = FALSE, particles = particles, r_cores = dots$r_cores,
           n_cores = cores_per_chain, total_cores = cores_per_chain,
           mc.cores = cores_for_chains, mc.set.seed = FALSE)
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
    cond <- .marginal_conditional(group_pars$mu, group_pars$var,
                                  .marginal_par_idx(colnames(proposals), marginalise))
    # Selection is by likelihood alone, so particles far below the best can
    # keep their first rule.
    spec <- .marginal_batch_spec(marginalise, cond, proposals)
    spec$log_offset <- numeric(n_particles)
    grid <- compute_marginal_grid(proposals, data_s, pmwgs$model, spec,
                                  r_cores = r_cores)
    lw <- marginal_ll_from_grid(grid)
    weight <- exp(lw - max(lw))
    idx <- sample(x = n_particles, size = 1, prob = weight)
    proposal <- proposals[idx,]
    names(proposal) <- colnames(proposals)
    out <- .marginal_draw(proposal, cond, as.numeric(grid$nodes[idx, ]),
                          as.numeric(grid$log_terms[idx, ]), grid$mode[idx],
                          grid$sd[idx], data_s, pmwgs$model)
    return(list(proposal = out$x, ll = out$ll))
  }
  # Pooled coefficients are drawn per particle from their prior here; init()
  # pools the selected draws into one start value.
  pooled_draws <- NULL
  if (!is.null(pmwgs$pooled)) {
    st <- pmwgs$pooled
    pooled_draws <- matrix(stats::rnorm(n_particles * length(st$value),
                                        st$prior_mean, sqrt(st$prior_var)),
                           nrow = n_particles, byrow = TRUE,
                           dimnames = list(NULL, names(st$value)))
  }
  lw <- calc_ll_manager(cbind(proposals, pooled_draws), dadm = data_s,
                        model = pmwgs$model, r_cores = r_cores)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx],
              pooled = if (is.null(pooled_draws)) NULL else pooled_draws[idx, ]))
}


check_tune_settings <- function(tune, n_pars, stage, particles){
  if (is.null(tune$search_width)) tune$search_width <- 1
  tune$p_accept <- set_p_accept(stage, tune$search_width)
  # Preburn searches with larger steps; later stages target the random-walk
  # scale (update_local_steps).
  tune$step_target <- if (stage == "preburn") tune$p_accept[1L] else
    step_target(tune$search_width)
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
    # Squared jumps are fractional, so they are not rounded like the counts;
    # a new mixture layout starts them afresh.
    if (stage == "adapt" ||
        length(pm_settings[[i]]$jump_counts) != length(pm_settings[[i]]$mix))
      pm_settings[[i]]$jump_counts <- numeric(length(pm_settings[[i]]$mix))
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
  pmwgs <- .emc_prepare_sampler_rng(pmwgs)
  if (!is.null(pmwgs$rng$gibbs)) {
    assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
  }
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
  pmwgs <- .emc_pooled_sync(pmwgs)
  # Each block refreshes the pooled move's curvature once (R/pooled.R).
  if (!is.null(pmwgs$pooled))
    pmwgs$pooled$block <- (if (is.null(pmwgs$pooled$block)) 0L else pmwgs$pooled$block) + 1L
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
  n_comp <- length(unique(tune$components))
  if(is.null(pm_settings)) pm_settings <- lapply(1:pmwgs$n_subjects, function(x) return(vector("list", n_comp)))
  # n_blocks > 1 re-blocks the parameters after preburn (and may re-cluster
  # between blocks of iterations), so the per-block settings must follow. New
  # blocks inherit the subject's adapted settings with fresh acceptance counts.
  pm_settings <- lapply(pm_settings, function(x) {
    if (length(x) == n_comp) return(x)
    template <- x[[1L]]
    template$proposal_counts <- NULL
    template$acc_counts <- NULL
    template$jump_counts <- NULL
    rep(list(template), n_comp)
  })
  burn_eff_ready <- stage == "burn" && isTRUE(pmwgs$burn_eff_ready)
  proposal_stage <- if (burn_eff_ready) "adapt" else stage
  tune <- check_tune_settings(tune, n_pars, proposal_stage, particles)
  if (burn_eff_ready) {
    pm_settings <- lapply(pm_settings, function(x) {
      for (i in seq_along(x))
        if (is.null(x[[i]]$mix)) x[[i]]$mix <- c(0.1, 0.25, 0.5, 0.15)
      x
    })
  }
  pm_settings <- lapply(pm_settings, FUN = check_sampling_settings,
                        stage = proposal_stage, n_pars = n_pars, particles)

  eff_mu <- .emc_attach_lik(pmwgs$eff_mu, pmwgs$lik_laplace)
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
  wpool_streams <- pmwgs$rng$subjects

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
  # Optional moves (R/move_allocation.R): the subject HMC plan is fixed for
  # the whole block, and the allocation's measurements are recorded every
  # iteration while any optional move is on.
  hmc_plan <- if (stage == "preburn") NULL else .emc_hmc_worker_plan(pmwgs)
  alloc_on <- stage != "preburn" && .emc_alloc_wanted(pmwgs)
  wpool_ctx <- list(
    data = data, model = .emc_wpool_slim_model(pmwgs$model),
    stage = stage, type = pmwgs$type,
    tune = tune, marginalise = pmwgs$marginalise,
    r_cores = pool_budget$likelihood, n_pars = pmwgs$n_pars,
    eff_mu = eff_mu, eff_var = eff_var, chains_mu = chains_mu,
    chains_var = chains_var, chol_caches = chol_caches,
    hmc = hmc_plan, par_names = par_names
  )
  sampler_profile <- .emc_profile_enabled()
  startup_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
  n_workers <- min(pool_budget$subject, pmwgs$n_subjects)
  wpool_part <- NULL
  wpool <- .emc_wpool_start(n_workers, wpool_ctx, persistent = TRUE)
  if (!is.null(wpool)) {
    on.exit(.emc_wpool_stop(wpool), add = TRUE)
  } else if (n_workers <= 1L) {
    wpool <- list(n = 1L, dir = NULL, jobs = list(), wcs = list(),
                  rcs = list(), alive = FALSE)
  }
  if (!is.null(wpool)) {
    wpool_cost <- .emc_subject_cost(data)
    wpool_part <- .emc_lpt_partition(wpool_cost, wpool$n)
  }

  llpool <- NULL
  if (pmwgs$n_subjects == 1L && pool_budget$likelihood > 1L &&
      isTRUE(getOption("emc2.ll_pool", TRUE))) {
    llpool <- .emc_wpool_start(pool_budget$likelihood, wpool_ctx,
                               persistent = TRUE)
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
  ensemble_rejects <- .emc_reject_counts("ensemble")
  shear_rejects <- .emc_reject_counts("shear")
  collapsed_rejects <- .emc_reject_counts("collapsed")
  pooled_rejects <- .emc_reject_counts("pooled")
  stage_rejects <- list(
    particle = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                               .EMC_REJECT_CLASSES),
    gibbs = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                            .EMC_REJECT_CLASSES),
    ensemble = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                                 .EMC_REJECT_CLASSES),
    shear = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                            .EMC_REJECT_CLASSES),
    collapsed = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                                .EMC_REJECT_CLASSES),
    pooled = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
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
    if (alloc_on) {
      # The allocation charges wall time on the cores this chain holds now
      # (a chain's share grows when other chains finish their block).
      el_iter <- proc.time()[["elapsed"]]
      cores_it <- .emc_cores_now(core_ctl, n_cores)
    }

    assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
    gibbs_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    pars_attempt <- tryCatch(
      gibbs_step(pmwgs, pmwgs$samples$alpha[!nuisance,,j-1], pmwgs$type),
      error = identity
    )
    if (inherits(pars_attempt, c("error", "try-error"))) {
      .emc_reject_record(pars_attempt, "gibbs")
      .emc_gibbs_abort(pars_attempt, stage = stage, iteration = j)
    }
    pars <- pars_comb <- pars_attempt
    if(any(nuisance)){
      pars_nuis_attempt <- tryCatch(
        gibbs_step(pmwgs$sampler_nuis, pmwgs$samples$alpha[nuisance,,j-1], pmwgs$sampler_nuis$type),
        error = identity
      )
      if (inherits(pars_nuis_attempt, c("error", "try-error"))) {
        .emc_reject_record(pars_nuis_attempt, "gibbs")
        .emc_gibbs_abort(pars_nuis_attempt, stage = stage, iteration = j,
                         nuisance = TRUE)
      }
      pars_nuis <- pars_nuis_attempt
      # pars_comb only feeds the particle step, so it carries the subject prior
      # covariance (SEM: subj_var, which excludes covariate spread).
      pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu,
                                     if (is.null(pars$subj_var)) pars$tvar else pars$subj_var,
                                     pars_nuis$tvar, nuisance, pars$subj_mu)
      pars_comb$alpha <- pmwgs$samples$alpha[,,j-1]
      pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                        group_level = pars_nuis,
                                                                        j = j,
                                                                        proposals = NULL,
                                                                        n_pars = n_pars, type = pmwgs$sampler_nuis$type)
      pmwgs$sampler_nuis$samples$idx <- j
    }
    # Ensemble pool members come from the group-step stream and depend only on
    # the frozen proposals, so their likelihoods ride along with the subject
    # step below.
    ens_q <- pmwgs$ensemble_q
    ens_pools <- .emc_ensemble_draw_pools(pmwgs, stage,
                                          if (is.null(ens_q)) chains_mu else ens_q$mu,
                                          if (is.null(ens_q)) chains_var else ens_q$var,
                                          iteration = j)
    pmwgs$rng$gibbs <- get(".Random.seed", envir = globalenv())
    gibbs_elapsed <- if (sampler_profile) {
      proc.time()[["elapsed"]] - gibbs_started
    } else NA_real_
    cache_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    group_var_it <- tryCatch(get_group_level(pars_comb, 1L, pmwgs$type)$var,
                             error = function(e) NULL)
    group_chol_it <- build_group_chol_cache(group_var_it, idx_list_blk,
                                          marginal_idx_blk,
                                          include_full = length(unique(tune$components)) > 1L)
    cache_elapsed <- if (sampler_profile) {
      proc.time()[["elapsed"]] - cache_started
    } else NA_real_
    recycle_elapsed <- grow_elapsed <- if (sampler_profile) 0 else NA_real_
    pool_profile <- NULL
    # Subjects whose particle step the allocation (R/move_allocation.R) skips
    # this iteration, each drawn from the subject's own stream. A subject
    # that also takes HMC at rate r draws one uniform U for both kernels:
    # the particle step runs if U < q and HMC if U >= 1 - r, which keeps
    # both rates and runs at least one of them whenever q + r >= 1. A
    # subject whose allocated HMC rate has since dropped to zero (its warmup
    # ended) takes its particle step.
    part_skip <- NULL
    hmc_pick <- NULL
    q_part <- if (alloc_on) pmwgs$alloc$particle$rate
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
      if (wpool$n != n_before || !is.null(part_skip)) {
        wpool_part <- .emc_lpt_partition(
          if (is.null(part_skip)) wpool_cost else wpool_cost * !part_skip, wpool$n)
      }
      it <- .emc_wpool_iter(wpool, wpool_ctx, wpool_part, pars_comb,
                            group_chol_it, pm_settings,
                            pmwgs$samples$subj_ll[, j - 1], wpool_streams,
                            extra = ens_pools$rows,
                            pooled = pmwgs$pooled$value, skip = part_skip)
      if (!is.null(ens_pools)) ens_pools$ll <- it$extra_ll
      proposals <- it$props
      pm_settings <- it$pm_settings
      wpool_streams <- it$seeds
      wpool$alive <- it$alive
      pool_profile <- it$profile
      if (!is.null(it$rejects)) {
        stage_rejects$particle <- stage_rejects$particle + it$rejects
      }
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
    core_budget <- .particle_core_budget(pmwgs$n_subjects, n_cores = n_cores,
                                         r_cores = r_cores,
                                         total_cores = n_cores)
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
        eff_mu = eff_mu[[s]], eff_var = eff_var[[s]],
        chains_mu = chains_mu[[s]], chains_var = chains_var[[s]],
        prev_ll = pmwgs$samples$subj_ll[s, j - 1L], parameters = pars_comb,
        model = pmwgs$model, stage = stage, type = pmwgs$type, tune = tune,
        marginalise = pmwgs$marginalise,
        r_cores = core_budget$likelihood, chol_cache = chol_caches[[s]],
        group_chol = group_chol_it)
      out$rng_seed <- get(".Random.seed", envir = globalenv())
      if (!is.null(ens_pools))
        out$extra_ll <- .emc_ensemble_ll_candidate(ens_pools$rows[[s]], data[[s]],
                                                   pmwgs$model)
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
      stage = stage, iteration = j)
    pm_settings <- assembled$pm_settings
    proposals <- assembled$proposals
    }
    if (alloc_on) {
      pmwgs <- .emc_alloc_record_particle(
        pmwgs, matrix(pmwgs$samples$alpha[, , j - 1L], nrow = pmwgs$n_pars),
        proposals[seq_len(pmwgs$n_pars), , drop = FALSE],
        if (is.null(part_skip)) rep(TRUE, pmwgs$n_subjects) else !part_skip,
        part_cpu_share * .emc_reserved_cost(proc.time()[["elapsed"]] - part_started, cores_it),
        if (is.null(wpool)) NULL else it$times)
    }

    # Ensemble group update (R/ensemble_move.R): after the subject step, so the
    # group Gibbs rounds and the final reselection close the sweep.
    ens <- .emc_ensemble_move(pmwgs, pars, proposals, stage, wpool, wpool_ctx,
                              wpool_part, ens_pools)
    # Then the shear (R/shear_move.R) when enabled: the non-centred move of
    # the same group state.
    shear_cpu <- shear_raw <- shear_el <- shear_round_el <- shear_round_cpu <- 0
    shear_build_el <- shear_ll_el <- 0
    shear_rounds <- 0L
    if (alloc_on) ens$sampler <- .emc_alloc_record_group(ens$sampler, ens$pars)
    if (.emc_shear_active(ens$sampler)) {
      shear_started <- .emc_cpu_self()
      shear_el0 <- proc.time()[["elapsed"]]
      sll_shear <- .emc_subject_ll_cpu()
      # The rounds in K parts with a subject HMC round between parts, so
      # the subjects re-equilibrate before the next rounds (emc2.shear_interleave).
      K <- if (is.null(hmc_plan)) 1L else
        max(1L, as.integer(getOption("emc2.shear_interleave", 1L)))
      rounds_all <- ens$sampler$shear$plan$rounds
      after <- identical(getOption("emc2.shear_interleave_mode", "between"), "after")
      for (part in seq_len(K)) {
        if (K > 1L) ens$sampler$shear$plan$rounds <- if (after) rounds_all * (part == 1L) else rounds_all / K
        sh <- .emc_shear_move(ens$sampler, ens$pars, ens$proposals, stage,
                              ens$wpool, wpool_part)
        ens$sampler <- sh$sampler
        ens$pars <- sh$pars
        ens$proposals <- sh$proposals
        ens$wpool <- sh$wpool
        shear_rounds <- shear_rounds + (if (is.null(sh$rounds)) 0L else sh$rounds)
        shear_round_el <- shear_round_el + (if (is.null(sh$round_el)) 0 else sh$round_el)
        shear_build_el <- shear_build_el + (if (is.null(sh$build_el)) 0 else sh$build_el)
        shear_ll_el <- shear_ll_el + (if (is.null(sh$ll_el)) 0 else sh$ll_el)
        if (part < K) {
          hr <- .emc_hmc_round(hmc_plan, ens$proposals, ens$pars, pmwgs$type, wpool_streams,
                               data, pmwgs$model, wpool = ens$wpool,
                               n_cores = pool_budget$subject, r_cores = pool_budget$likelihood,
                               par_names = par_names)
          ens$proposals <- hr$proposals
          wpool_streams <- hr$streams
          hmc_plan <- hr$plan
          if (!is.null(ens$wpool) && !is.null(hr$pool)) ens$wpool$alive <- hr$pool$alive
        }
      }
      if (K > 1L) ens$sampler$shear$plan$rounds <- rounds_all
      shear_raw <- .emc_cpu_self() - shear_started + .emc_subject_ll_cpu() - sll_shear
      shear_el <- proc.time()[["elapsed"]] - shear_el0
      if (alloc_on) {
        shear_cpu <- .emc_reserved_cost(shear_el, cores_it)
        shear_round_cpu <- .emc_reserved_cost(shear_round_el, cores_it)
      }
    }
    # Then the collapsed move (R/collapsed_move.R) when enabled: the whole
    # group state with the subjects in their conditional frames.
    coll_cpu <- coll_el <- coll_round_cpu <- 0
    coll_moves <- 0L
    if (.emc_collapsed_active(ens$sampler)) {
      coll_el0 <- proc.time()[["elapsed"]]
      moves_before <- ens$sampler$collapsed$plan$block$n
      cm <- .emc_collapsed_move(ens$sampler, ens$pars, ens$proposals, stage,
                                ens$wpool, wpool_part)
      ens$sampler <- cm$sampler
      ens$pars <- cm$pars
      ens$proposals <- cm$proposals
      ens$wpool <- cm$wpool
      moves_after <- ens$sampler$collapsed$plan$block$n
      if (length(moves_before) && length(moves_after))
        coll_moves <- as.integer(moves_after - moves_before)
      coll_el <- proc.time()[["elapsed"]] - coll_el0
      if (alloc_on) {
        coll_cpu <- .emc_reserved_cost(coll_el, cores_it)
        coll_round_cpu <- .emc_reserved_cost(if (is.null(cm$round_el)) 0 else cm$round_el, cores_it)
      }
    }
    ensemble_delta <- .emc_reject_delta(ensemble_rejects, "ensemble")
    ensemble_rejects <- .emc_reject_counts("ensemble")
    stage_rejects$ensemble <- stage_rejects$ensemble + ensemble_delta
    shear_delta <- .emc_reject_delta(shear_rejects, "shear")
    shear_rejects <- .emc_reject_counts("shear")
    stage_rejects$shear <- stage_rejects$shear + shear_delta
    collapsed_delta <- .emc_reject_delta(collapsed_rejects, "collapsed")
    collapsed_rejects <- .emc_reject_counts("collapsed")
    stage_rejects$collapsed <- stage_rejects$collapsed + collapsed_delta
    pmwgs <- ens$sampler
    pars <- ens$pars
    proposals <- ens$proposals
    wpool <- ens$wpool

    # Pooled coefficients given the final subject states (R/pooled.R).
    pooled_delta <- NULL
    if (!is.null(pmwgs$pooled)) {
      pmv <- .emc_pooled_move(pmwgs, proposals, stage, pars, group_chol_it,
                              pm_settings, wpool_streams, wpool, wpool_ctx,
                              wpool_part, n_cores = n_cores)
      pmwgs <- pmv$sampler
      proposals <- pmv$proposals
      pars <- pmv$pars
      if (!is.null(wpool) && !is.null(pmv$alive)) wpool$alive <- pmv$alive
      data <- pmwgs$data
      wpool_ctx$data <- data
      pooled_delta <- .emc_reject_delta(pooled_rejects, "pooled")
      pooled_rejects <- .emc_reject_counts("pooled")
      stage_rejects$pooled <- stage_rejects$pooled + pooled_delta
    }

    # Subject HMC (R/subject_hmc.R): the last update of the iteration, from
    # the final group state.
    hmc_tasks <- NULL
    hmc_master <- hmc_el <- 0
    if (alloc_on)
      pmwgs <- .emc_alloc_record_subjects(pmwgs, proposals[seq_len(pmwgs$n_pars), , drop = FALSE],
                                          pars$tmu)
    if (!is.null(hmc_plan)) {
      hmc_started <- .emc_cpu_self()
      hmc_el0 <- proc.time()[["elapsed"]]
      hr <- .emc_hmc_round(hmc_plan, proposals, pars, pmwgs$type, wpool_streams,
                           data, pmwgs$model, wpool = wpool,
                           n_cores = pool_budget$subject, r_cores = pool_budget$likelihood,
                           par_names = par_names, pick = hmc_pick)
      proposals <- hr$proposals
      wpool_streams <- hr$streams
      hmc_tasks <- hr$tasks
      hmc_plan <- hr$plan
      if (!is.null(wpool) && !is.null(hr$pool)) wpool$alive <- hr$pool$alive
      hmc_master <- .emc_cpu_self() - hmc_started
      hmc_el <- proc.time()[["elapsed"]] - hmc_el0
      attempted <- !vapply(hmc_tasks, is.null, logical(1))
      if (any(attempted)) {
        for (s in which(attempted)) {
          hmc_tasks[[s]]$work <- hmc_tasks[[s]]$cpu
          hmc_tasks[[s]]$cpu <- hmc_tasks[[s]]$cpu + hmc_master / sum(attempted)
        }
        # The round's wall time on the chain's cores, shared by CPU.
        if (alloc_on) {
          cpu_s <- vapply(hmc_tasks[attempted], function(x) x$cpu, numeric(1))
          share <- if (sum(cpu_s) > 0) cpu_s / sum(cpu_s) else rep(1 / sum(attempted), sum(attempted))
          for (k in seq_along(share))
            hmc_tasks[[which(attempted)[k]]]$cost <- .emc_reserved_cost(hmc_el, cores_it) * share[k]
        }
      }
      pmwgs$subject_hmc <- .emc_hmc_record_tasks(pmwgs$subject_hmc, hmc_tasks)
    }

    fill_started <- if (sampler_profile) proc.time()[["elapsed"]] else NA_real_
    if (is.null(ens_pools) && pmwgs$type %in% c("factor", "infnt_factor", "SEM") &&
        !any(pmwgs$nuisance) && !any(as.logical(pmwgs$marginalised_idx), na.rm = TRUE)) {
      assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
      pars <- .emc_ensemble_refresh_scores(pars, pmwgs, proposals[seq_len(pmwgs$n_pars), , drop = FALSE])
      pmwgs$rng$gibbs <- get(".Random.seed", envir = globalenv())
    }
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, j = j, n_pars = pmwgs$n_pars, type = pmwgs$type)
    if (!is.null(pmwgs$pooled)) pmwgs$samples$pooled[, j] <- pmwgs$pooled$value
    if (alloc_on) {
      rest_cpu <- .emc_reserved_cost(proc.time()[["elapsed"]] - el_iter - shear_el - coll_el -
                                       hmc_el, cores_it)
      pmwgs <- .emc_alloc_record_cost(pmwgs, rest_cpu, hmc_tasks, moves = list(
        shear = c(cpu = shear_cpu, own = shear_round_cpu, n = shear_rounds),
        collapsed = c(cpu = coll_cpu, own = coll_round_cpu, n = coll_moves),
        hmc = .emc_reserved_cost(hmc_el, cores_it)))
    }
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
        degraded = pp$degraded,
        shear_rounds = shear_rounds, shear_cpu = shear_raw,
        shear_elapsed = shear_el, shear_build = shear_build_el, shear_ll = shear_ll_el,
        hmc_elapsed = hmc_el,
        hmc_task_cpu_max = max(c(0, vapply(hmc_tasks, function(x)
          if (is.null(x) || is.null(x$work)) 0 else x$work, numeric(1)))),
        hmc_attempts = sum(!vapply(hmc_tasks, is.null, logical(1))),
        hmc_accepted = sum(vapply(hmc_tasks, function(x) isTRUE(x$accepted), logical(1))),
        hmc_cpu = sum(vapply(hmc_tasks, function(x) if (is.null(x)) 0 else x$cpu, numeric(1)))),
        .emc_reject_row_fields(pp$rejects, "particle"),
        .emc_reject_row_fields(gibbs_delta, "gibbs"),
        .emc_reject_row_fields(ensemble_delta, "ensemble"),
        .emc_reject_row_fields(shear_delta, "shear"),
        .emc_reject_row_fields(collapsed_delta, "collapsed"),
        .emc_reject_row_fields(pooled_delta, "pooled")))
    }
  }
  .emc_failure_report_stage(stage_rejects, stage)
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  pmwgs$rng$subjects <- wpool_streams
  assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
  if (sampler_profile) {
    attr(pmwgs$samples, "sampler_profile") <- .emc_profile_bind(profile_rows)
  }
  if (verboseProgress) close(pb)
  return(pmwgs)
}

.emc_sample_iteration_equal <- function(samples, j) {
  if(is.null(samples$stage) || j <= 1L) return(FALSE)
  n_iter <- length(samples$stage)
  if(j > n_iter) return(FALSE)
  names_history <- names(samples)[vapply(samples, function(obj) {
    is_iteration_array(obj, n_iter)
  }, logical(1))]
  names_history <- setdiff(names_history, "stage")
  if(!length(names_history)) return(FALSE)
  same_slice <- function(obj) {
    d <- dim(obj)
    # seq_len, not TRUE: TRUE cannot index a dimension of extent zero.
    idx_prev <- lapply(d, seq_len)
    idx_now <- idx_prev
    idx_prev[[length(d)]] <- j - 1L
    idx_now[[length(d)]] <- j
    prev <- do.call(`[`, c(list(obj), idx_prev, list(drop = FALSE)))
    now <- do.call(`[`, c(list(obj), idx_now, list(drop = FALSE)))
    # Ignore per-iteration dimnames, if a caller supplied them, but require
    # exact numerical equality so slow movement is not mistaken for a stall.
    identical(unname(prev), unname(now))
  }
  all(vapply(samples[names_history], same_slice, logical(1)))
}

# The last iteration of every iteration-indexed field -- the same fields
# .emc_sample_iteration_equal() compares.  run_emc() takes one per try: two
# consecutive tries ending on bit-identical snapshots mean the chain did not
# move for the whole of the later try, because a continuous proposal does not
# return to the same bits by accident.
.emc_chain_snapshot <- function(samples) {
  n_iter <- length(samples$stage)
  if (!n_iter) return(NULL)
  nms <- names(samples)[vapply(samples, function(obj) {
    is_iteration_array(obj, n_iter)
  }, logical(1))]
  nms <- setdiff(nms, "stage")
  lapply(samples[nms], function(obj) {
    d <- dim(obj)
    idx <- lapply(d, seq_len)
    idx[[length(d)]] <- n_iter
    unname(do.call(`[`, c(list(obj), idx, list(drop = FALSE))))
  })
}

reject_sample_iteration <- function(samples, j) {
  if (j <= 1) {
    samples$idx <- j
    return(samples)
  }
  # `last_theta_var_inv` is the current Gibbs-state covariance inverse, not a
  # draw history.  It is a square 2-D matrix, so treating every 2-D object as
  # an iteration matrix silently overwrites its columns whenever a Gibbs
  # update is rejected.  The resulting repeated-column matrix is neither a
  # covariance nor a valid sample store and later makes `filter_obj()` try to
  # reshape (p x p) values to (p x n_iter).
  n_iter <- length(samples$stage)
  for (nm in names(samples)) {
    if (nm %in% .EMC_STATIC_SAMPLE_FIELDS) next
    obj <- samples[[nm]]
    d <- dim(obj)
    if (is.null(d)) next
    # Use the allocated history length, rather than just the current column
    # count.  This is the same rule used by extend_obj() and avoids mistaking
    # any other fixed matrix for a history when its dimensions happen to be
    # large enough to contain j.
    if (is_iteration_array(obj, n_iter) && d[length(d)] >= j) {
      if (length(d) == 2) {
        samples[[nm]][, j] <- samples[[nm]][, j - 1]
      } else if (length(d) == 3) {
        samples[[nm]][, , j] <- samples[[nm]][, , j - 1]
      }
    }
  }
  samples$idx <- j
  samples
}

# The subject's likelihood approximation (.emc_lik_canonical) rides on its
# efficient-proposal mean, so it reaches new_particle() on every path.
.emc_attach_lik <- function(eff_mu, lik) {
  if (is.null(eff_mu)) return(eff_mu)
  for (s in seq_along(eff_mu)) attr(eff_mu[[s]], "lik") <- NULL
  if (is.null(lik) || !isTRUE(getOption("emc2.laplace_conditional", TRUE))) return(eff_mu)
  for (s in seq_along(eff_mu)) {
    if (!is.null(eff_mu[[s]]) && length(lik) >= s && is.list(lik[[s]]) &&
        !is.null(lik[[s]]$check) && !isTRUE(.emc_laplace_poor(lik[[s]]$check)))
      attr(eff_mu[[s]], "lik") <- lik[[s]]
  }
  eff_mu
}

# Efficient proposal from the likelihood approximation and the current group
# prior N(mu, V): the Gaussian product, precision prec + V^-1 and mean
# (prec + V^-1)^-1 (lin + V^-1 mu). It is the conditional of alpha_s at this
# group state, and depends on the group state only, so the conditional-IS
# step stays invariant. NULL when there is no approximation or it is unusable.
.emc_laplace_conditional <- function(eff_mu, group_pars) {
  lik <- attr(eff_mu, "lik")
  if (is.null(lik)) return(NULL)
  mu <- as.numeric(group_pars$mu)
  V <- as.matrix(group_pars$var)
  p <- length(mu)
  if (!is.matrix(lik$prec) || nrow(lik$prec) != p || length(lik$lin) != p ||
      !identical(dim(V), c(p, p))) return(NULL)
  Vi <- tryCatch(chol2inv(chol(V)), error = function(e) NULL)
  if (is.null(Vi)) return(NULL)
  Q <- lik$prec + Vi
  R <- tryCatch(chol((Q + t(Q)) / 2), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  Sig <- chol2inv(R)
  m <- drop(Sig %*% (lik$lin + Vi %*% mu))
  if (!all(is.finite(m)) || !all(is.finite(Sig))) return(NULL)
  list(mu = m, var = (Sig + t(Sig)) / 2)
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
    attempt <- reject_particle(old, prev_ll, pm_settings)
  }
  # Tag the subject this result belongs to.  parallel::mclapply DROPS the
  # element of a worker that dies below R level, so position is not a reliable
  # identifier: without the tag the survivors shift into the wrong subjects'
  # slots.  Callers that index by name (the worker pool) ignore it.
  attempt$subject <- s
  attempt
}

# A batch likelihood that fails is retried one particle at a time, and a
# particle that still fails gets zero density (-Inf), as in the ensemble pool.
# Keeping the old state on every failure would stay put with a probability
# that depends on the state, which is not invariant; a zero-density candidate
# is. The reference particle is the current state, so if it fails the error is
# not a property of the candidates and is raised (the caller keeps the state).
.emc_particle_ll_rowwise <- function(err, proposals, data, model, component) {
  cls <- if (inherits(err, "emc_failure") && !is.null(err$emc_class))
    err$emc_class else .emc_classify_failure(err)
  if (identical(.emc_failure_action(cls), "abort")) stop(err)
  ll <- vapply(seq_len(nrow(proposals)), function(r) {
    v <- tryCatch(as.numeric(calc_ll_manager(
      proposals[r, , drop = FALSE], dadm = data, model = model,
      component = component, r_cores = 1L)), error = function(e) NA_real_)
    if (length(v) != 1L || is.na(v) || v == Inf) -Inf else v
  }, numeric(1))
  if (!is.finite(ll[1L])) stop(err)
  .emc_reject_record(err, "particle")
  ll
}

reject_particle <- function(subj_mu, prev_ll, pm_settings) {
  list(proposal = subj_mu, ll = prev_ll, pm_settings = pm_settings)
}

# Reassemble one iteration's particle results into the (n_pars + 1) x n_subjects
# proposal matrix, matching each result to its subject by the tag
# safe_new_particle() attaches.
#
# Matching by tag rather than by position is the whole point.  A worker that
# dies below R level is not an R error -- mclapply drops its element and warns
# "N function calls resulted in an error", nothing more.  The result therefore
# comes back SHORT.  The previous code did
#     array(unlist(proposals[1:2, ]), dim = c(n_pars + 1, n_subjects))
# which silently RECYCLED the survivors to fill the gap, so one dead worker
# handed another subject's parameters to the subject that lost, and the chain
# wedged with no error raised.  Here a lost subject repeats its own previous
# state -- exactly what safe_new_particle() does for an error it CAN catch --
# and no other subject is touched.
.emc_assemble_proposals <- function(raw, n_pars, n_subjects, prev_alpha,
                                    prev_ll, prev_pm, stage = NULL,
                                    iteration = NULL) {
  props <- matrix(NA_real_, nrow = n_pars + 1L, ncol = n_subjects)
  pm <- vector("list", n_subjects)
  seen <- logical(n_subjects)
  for (res in raw) {
    if (inherits(res, c("error", "try-error")) || !is.list(res)) next
    s <- res[["subject"]]
    if (!is.numeric(s) || length(s) != 1L || is.na(s) ||
        s < 1L || s > n_subjects || seen[[s]]) next
    prop <- res[["proposal"]]
    ll <- res[["ll"]]
    if (!is.numeric(prop) || length(prop) != n_pars ||
        !is.numeric(ll) || length(ll) != 1L) next
    props[, s] <- c(prop, ll)
    pm[[s]] <- res[["pm_settings"]]
    seen[[s]] <- TRUE
  }
  if (!all(seen)) {
    lost <- which(!seen)
    for (s in lost) {
      props[, s] <- c(prev_alpha[, s], prev_ll[[s]])
      pm[[s]] <- prev_pm[[s]]
    }
    # Under the strict policy this stops; otherwise it warns and the repeated
    # states above stand.
    .emc_worker_loss(lost, stage = stage, iteration = iteration)
  }
  list(proposals = props, pm_settings = pm)
}



# Cholesky factor for a proposal component. A matrix that does not factor is
# repaired (jitter, then its diagonal) and flagged: the repaired factor is a
# valid proposal as long as drawing and the density use it alike, but it is
# not the target's group prior, so the particle step never evaluates the prior
# with a repaired factor.
.chol_factor <- function(Sigma, p) {
  if (p <= 0) return(NULL)
  base_R <- tryCatch(chol(Sigma), error = function(e) NULL)
  repaired <- is.null(base_R)
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
  list(R = base_R, rooti = rooti, sum_log_diag = sum(log(diag(rooti))), p = p,
       repaired = repaired)
}

.marginal_par_idx <- function(par_names, marginalise) {
  if (is.null(marginalise) || is.null(par_names)) {
    return(rep(FALSE, length(par_names)))
  }
  par_names %in% marginalise$param
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

# Mean of the efficient proposal at the current group state. The proposal is
# the conditional of alpha_s given the group mean, fitted to earlier draws:
# its covariance is fixed and its mean is linear in the conditioning value, so
# attr(eff_mu, "B") and attr(eff_mu, "given") let it follow the group mean
# every iteration. It depends on the group state and not on alpha_s, so the
# conditional-IS step stays invariant. Without the attributes the stored mean
# is used as is.
.emc_eff_mean <- function(eff_mu, given_now) {
  B <- attr(eff_mu, "B")
  given <- attr(eff_mu, "given")
  if (is.null(eff_mu) || is.null(B) || is.null(given)) return(eff_mu)
  if (length(given_now) != length(given) || nrow(B) != length(eff_mu))
    return(eff_mu)
  delta <- given_now - given
  delta[!is.finite(delta)] <- 0
  out <- as.numeric(eff_mu) + drop(B %*% delta)
  if (any(!is.finite(out))) return(as.numeric(eff_mu))
  out
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
  # The efficient proposal conditions on the subject's prior mean (for SEM,
  # given its covariates: filtered_samples_SEM's conditional_mu), the same
  # value on the serial and the direct worker paths.
  lap <- .emc_laplace_conditional(eff_mu, group_pars)
  if (!is.null(lap)) {
    eff_mu <- lap$mu
    eff_var <- lap$var
  } else eff_mu <- .emc_eff_mean(eff_mu, as.numeric(group_pars$mu))
  unq_components <- unique(tune$components)
  proposal_out <- numeric(length(group_pars$mu))
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- if (direct_state) current_alpha else parameters$alpha[,s]
  marg_cond <- NULL
  marg_pred <- NULL
  marginal_idx <- rep(FALSE, length(subj_mu))
  if (!is.null(marginalise)) {
    if (is.null(names(subj_mu))) {
      alpha_names <- if (direct_state) names(current_alpha) else rownames(parameters$alpha)
      names(subj_mu) <- alpha_names[seq_along(subj_mu)]
    }
    marginal_idx <- .marginal_par_idx(names(subj_mu), marginalise)
    # The particles see the other parameters' marginal group prior and the
    # likelihood with the marginalised parameter integrated against its
    # conditional group prior; each step's selection is then accepted or
    # rejected jointly with a fresh value of the marginalised parameter
    # (.marginal_delayed_accept).
    marg_cond <- .marginal_conditional(group_mu, group_var, marginal_idx)
    t_cur <- subj_mu[[which(marginal_idx)]]
    ll_cur <- prev_ll
    # The subject's posterior regression of the marginalised parameter on the
    # others predicts where each particle's integrand peaks.
    marg_pred <- if (!is.null(chains_mu) && !is.null(chains_var))
      tryCatch(.marginal_conditional(chains_mu, chains_var, marginal_idx),
               error = function(e) NULL)
  }
  # Likelihood at the final state, per shared-likelihood component: blocks
  # that share one likelihood each end with its value at the state so far.
  ll_shared <- numeric(0)
  particle_multiplier <- 1
  if(stage == "preburn"){
    Mus <- list(group_mu, subj_mu)
    Sigmas <- list(group_var, group_var)
    sig_tags <- c("group", "group")
    state_centred <- c(FALSE, TRUE)
    particle_multiplier <- 2
  } else if(stage == "burn" && length(pm_settings[[1L]]$mix) < 4L){
    Mus <- list(group_mu, subj_mu, subj_mu)
    Sigmas <- list(group_var, group_var, chains_var)
    sig_tags <- c("group", "group", "chains")
    state_centred <- c(FALSE, TRUE, TRUE)
  } else if(stage == "adapt" &&
            (is.null(eff_mu) || length(pm_settings[[1]]$mix) < 4L)){
    Mus <- list(group_mu, subj_mu, chains_mu)
    Sigmas <- list(group_var, chains_var, chains_var)
    sig_tags <- c("group", "chains", "chains")
    state_centred <- c(FALSE, TRUE, FALSE)
  } else{
    Mus <- list(group_mu, subj_mu, chains_mu, eff_mu)
    Sigmas <- list(group_var, chains_var, chains_var, eff_var)
    sig_tags <- c("group", "chains", "chains", "eff")
    state_centred <- c(FALSE, TRUE, FALSE, FALSE)
  }
  n_proposals <- length(Mus)
  cache_ok <- function(cache, ref_name, ref, idx) {
    !is.null(cache) && !is.null(cache[[ref_name]]) &&
      identical(cache[[ref_name]], ref) &&
      identical(cache$idx_list[[i]], idx)
  }
  for(i in unq_components){
    epsilons <- c(1, pm_settings[[i]]$epsilon)
    # Independence proposals are not shrunk: a scale below one on a fitted
    # covariance removes their coverage in high dimension. Only the local
    # (state-centred) step size is tuned.
    epsilons <- rep_len(epsilons, n_proposals)
    epsilons[!state_centred] <- 1
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

    # A component centred on the current state makes the reference particle's
    # proposal density depend on itself, and conditional importance sampling is
    # then not invariant for the posterior. Such components are centred on an
    # auxiliary point z ~ N(current, S_k) instead, targeting the extended
    # density pi(x) N(z | x, S_k): given z the proposal no longer depends on the
    # reference particle, and the extra factor enters the weights.
    #
    # That factor would also penalise every independent-proposal particle far
    # from z, so the update runs as a composition of conditional-IS steps (each
    # invariant): the independent components (group, chains, eff) with plain
    # weights, then each state-centred component on its own around the winner.
    # Two state-centred components in one step would multiply their auxiliary
    # densities, so the tighter one would pin every particle to the current
    # state and each step size would set the other's acceptance.
    particle_numbers <- numbers_from_proportion(pm_settings[[i]]$mix, pm_settings[[i]]$n_particles*particle_multiplier)
    mix_i <- pm_settings[[i]]$mix
    split_steps <- p_idx > 0 && any(state_centred)
    steps <- if (split_steps) {
      c(if (any(!state_centred)) list(which(!state_centred)),
        as.list(which(state_centred)))
    } else list(seq_len(n_proposals))
    shared_idx <- tune$shared_ll_idx[idx_full][1]
    is_shared <- shared_idx == tune$shared_ll_idx
    ref <- subj_mu[idx]
    ref_ll <- if (!is.na(ll_shared[as.character(shared_idx)])) ll_shared[as.character(shared_idx)] else prev_ll
    better <- numeric(n_proposals)
    jump <- numeric(n_proposals)
    weight_ess <- 0
    # Squared jumps for the mixture weights are measured in one metric for
    # every component: the chains covariance where there is one, else the
    # group covariance.
    jump_rooti <- if (p_idx > 0) {
      if (!is.null(bases[["chains"]])) bases[["chains"]]$rooti else bases[[sig_tags[1L]]]$rooti
    } else NULL
    marg_seen <- list()
    for (ks in steps) {
      if (sum(particle_numbers[ks]) == 0) next
      mus_i <- lapply(Mus, function(m) m[idx])
      aux_k <- if (p_idx > 0) intersect(ks, which(state_centred)) else integer(0)
      for (k in aux_k) {
        mus_i[[k]] <- particle_draws(1L, ref, covar = NULL, R = Rs[[k]])[1L, ]
      }
      proposals <- vector("list", length(ks) + 1L)
      proposals[[1]] <- matrix(ref, nrow = 1L)
      for (j in seq_along(ks)) {
        k <- ks[j]
        if (p_idx > 0) {
          proposals[[j + 1]] <- particle_draws(
            particle_numbers[k], mus_i[[k]], covar = NULL, R = Rs[[k]])
        } else {
          proposals[[j + 1]] <- matrix(numeric(0), nrow = particle_numbers[k], ncol = 0L)
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

      if (is.null(marginalise)) {
        comp_arg <- if (length(unq_components) > 1L) shared_idx else NULL
        lw <- tryCatch(
          calc_ll_pooled(proposals[,is_shared], dadm = data, model,
                         component = comp_arg, r_cores = r_cores, s = s,
                         varying = idx[is_shared]),
          error = identity)
        if (inherits(lw, "error"))
          lw <- .emc_particle_ll_rowwise(lw, proposals[, is_shared, drop = FALSE],
                                         data, model, comp_arg)
      }
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
        # The group component's factor doubles as the target prior. If it had
        # to be repaired the group covariance is not positive definite and the
        # prior is undefined there: reject rather than evaluate a substitute.
        if (p_idx > 0 && isTRUE(bases[[sig_tags[1L]]]$repaired))
          stop(errorCondition(
            "group covariance is not positive definite; subject prior undefined",
            class = c("emc_group_prior_failure", "emc_failure"),
            emc_class = "numerical", emc_source = "particle"))
        prior_density <- lp
      }
      log_mix_comps <- matrix(0, nrow = nrow(proposals), ncol = length(ks))
      for (j in seq_along(ks)) {
        k <- ks[j]
        log_mix_comps[, j] <- log(mix_i[k]) + if (k == 1L) lp else if (p_idx > 0) {
          fast_dmvnorm_rooti(
            x = proposals[, idx, drop = FALSE], mean = mus_i[[k]],
            rooti = rootis[[k]], log_const = log_consts[[k]])
        } else 0
      }
      max_log <- apply(log_mix_comps, 1L, max)
      lm <- max_log + log(rowSums(exp(log_mix_comps - max_log)))
      infnt_idx <- is.infinite(lm) | is.na(lm)
      if (any(infnt_idx)) {
        fin_vals <- lm[!infnt_idx]
        lm[infnt_idx] <- if (length(fin_vals) > 0) min(fin_vals) else -1e10
      }
      aux_density <- numeric(nrow(proposals))
      for (k in aux_k) {
        aux_density <- aux_density + fast_dmvnorm_rooti(
          x = proposals[, idx, drop = FALSE], mean = mus_i[[k]],
          rooti = rootis[[k]], log_const = log_consts[[k]])
      }
      if (!is.null(marginalise)) {
        # The rest of each particle's log weight lets the quadrature leave
        # particles that cannot be selected on their first rule.
        marg_spec <- .marginal_batch_spec(marginalise, marg_cond, proposals,
                                          pm_settings[[i]]$marg_state, marg_pred)
        marg_spec$log_offset <- prior_density + aux_density - lm
        marg_grid <- compute_marginal_grid(proposals[, is_shared, drop = FALSE],
                                           data, model, marg_spec,
                                           r_cores = r_cores)
        lw <- marginal_ll_from_grid(marg_grid)
        marg_seen[[length(marg_seen) + 1L]] <- list(
          mode = marg_grid$mode, sd = marg_grid$sd, x = proposals)
      }
      lw_total <- lw + ref_ll - lw[1]
      l <- lw_total + prior_density + aux_density - lm
      weights <- exp(l - max(l))
      idx_ll <- sample(x = nrow(proposals), size = 1, prob = weights)
      if (!is.null(marginalise)) {
        # A rejected step keeps the current state, as if the reference
        # particle had been selected.
        proposals[, marginal_idx] <- t_cur
        da <- .marginal_delayed_accept(proposals, idx_ll, marg_grid, lw,
                                       marg_cond, t_cur, ll_cur, data, model)
        if (da$accept) {
          t_cur <- da$t
          ll_cur <- da$ll
        } else idx_ll <- 1L
      }
      # Diagnostics only: which component supplied the selected particle and
      # how far it moved the state.
      comp_of <- c(0L, rep(ks, particle_numbers[ks]))
      if (is.null(pm_settings[[i]]$sel[[stage]])) {
        pm_settings[[i]]$sel[[stage]] <- list(n = numeric(n_proposals),
                                              disp = numeric(n_proposals),
                                              steps = 0)
      } else if (length(pm_settings[[i]]$sel[[stage]]$n) < n_proposals) {
        extra <- n_proposals - length(pm_settings[[i]]$sel[[stage]]$n)
        pm_settings[[i]]$sel[[stage]]$n <- c(pm_settings[[i]]$sel[[stage]]$n,
                                               numeric(extra))
        pm_settings[[i]]$sel[[stage]]$disp <- c(pm_settings[[i]]$sel[[stage]]$disp,
                                                  numeric(extra))
      }
      pm_settings[[i]]$sel[[stage]]$steps <- pm_settings[[i]]$sel[[stage]]$steps + 1
      if (comp_of[idx_ll] > 0L) {
        kk <- comp_of[idx_ll]
        pm_settings[[i]]$sel[[stage]]$n[kk] <- pm_settings[[i]]$sel[[stage]]$n[kk] + 1
        pm_settings[[i]]$sel[[stage]]$disp[kk] <- pm_settings[[i]]$sel[[stage]]$disp[kk] +
          sqrt(sum((proposals[idx_ll, idx] - ref)^2))
        if (!is.null(jump_rooti))
          jump[kk] <- jump[kk] +
            sum((matrix(proposals[idx_ll, idx] - ref, nrow = 1L) %*% jump_rooti)^2)
      }

      offset <- 1L
      for (k in ks) {
        n_k <- particle_numbers[k]
        if (n_k > 0) {
          better[k] <- sum(weights[offset + seq_len(n_k)] > weights[1])
          offset <- offset + n_k
        }
      }
      # The particle count follows the weight ESS of the step holding the
      # group proposal (the independence step). A local step's weights say
      # nothing about how well the fitted proposals cover the posterior.
      if (1L %in% ks) weight_ess <- sum(weights)^2 / sum(weights^2)
      ref <- proposals[idx_ll, idx]
      ref_ll <- if (is.null(marginalise)) lw[idx_ll] else ll_cur
    }

    ll_shared[as.character(shared_idx)] <- ref_ll
    proposal_out[idx] <- ref
    # Gibbs over blocks: later blocks condition on this block's new values.
    subj_mu[idx] <- ref
    if (!is.null(marginalise)) {
      subj_mu[marginal_idx] <- t_cur
      # The next rules are placed by the peaks of every step's particles.
      marg_next <- .marginal_state(
        list(mode = unlist(lapply(marg_seen, `[[`, "mode")),
             sd = unlist(lapply(marg_seen, `[[`, "sd"))),
        do.call(rbind, lapply(marg_seen, `[[`, "x")), marg_cond, marg_pred)
    }
    # Proposal settings adapt only in the discarded stages. Retained samples
    # come from one fixed transition kernel.
    if (stage != "sample") {
      pm_settings[[i]] <- update_pm_settings(
        pm_settings[[i]], better, particle_numbers, weight_ess, tune, sum(idx),
        layout = ifelse(state_centred, sig_tags, "")[-1L], jump = jump)
    }
    # The state places the quadrature rules, so it is part of the likelihood
    # approximation: frozen in the sample stage like the other proposal
    # settings, so the approximate target does not follow the chain.
    if (!is.null(marginalise) && stage != "sample" && !is.null(marg_next)) {
      pm_settings[[i]]$marg_state <- marg_next
    }
  }
  names(proposal_out) <- names(subj_mu)
  ll_out <- sum(ll_shared)
  if (!is.null(marg_cond)) {
    # The stored log-likelihood is the complete vector's, like every other
    # stored alpha.
    proposal_out[marginal_idx] <- t_cur
    ll_out <- ll_cur
  }
  return(list(proposal = proposal_out, ll = ll_out, pm_settings = pm_settings))
}


# better[j]: component j's particles that outweighed the reference particle of
# the conditional-IS step they were drawn in; jump[j]: summed squared jump
# (in the chains metric) of the moves component j supplied; weight_ess: of the
# independence step.
# layout: per non-group component, the covariance a state-centred component
# steps on ("group", "chains") or "" for an independence component.
#
# Mixture weights. Every particle costs one likelihood evaluation, so a
# component is worth its expected squared jump per particle:
# jump_counts / proposal_counts. The weights move towards those values in
# proportion. This is one currency for every component, whether its moves come
# from an independence step or a local one, and it does not feed on itself: a
# larger weight gives a component more particles but lowers what each of them
# contributes, because they compete for the same selection.
update_pm_settings <- function(pm_settings, better, particle_numbers,
                               weight_ess, tune, n_pars, layout = NULL,
                               jump = NULL) {
  pm_settings$iter <- pm_settings$iter + 1
  if (pm_settings$iter > tune$n0) {
    pm_settings$proposal_counts <- pm_settings$proposal_counts + particle_numbers
    pm_settings$acc_counts <- pm_settings$acc_counts + better
    if (!is.null(jump)) {
      jc <- pm_settings$jump_counts
      if (length(jc) != length(jump)) jc <- numeric(length(jump))
      pm_settings$jump_counts <- jc + jump
    }

    pm_settings <- update_local_steps(pm_settings, better[-1L],
                                      particle_numbers[-1L], tune, n_pars,
                                      layout)

    jc <- pm_settings$jump_counts
    if (length(pm_settings$mix) > 2 && length(jc) == length(pm_settings$mix)) {
      per_particle <- ifelse(pm_settings$proposal_counts > 0,
                             jc / pm_settings$proposal_counts, 0)
      if (sum(per_particle) > 0) {
        performance <- per_particle / sum(per_particle)
        new_mix <- (1 - tune$mix_adapt) * pm_settings$mix + tune$mix_adapt * performance
        new_mix <- pmax(new_mix, 0.02)
        pm_settings$mix <- new_mix / sum(new_mix)
      }
    }

    if (length(pm_settings$mix) > 3 && pm_settings$gd_good && weight_ess > 0) {
      ess <- weight_ess
      desired_ess <- tune$target_ESS
      scale_factor <- (desired_ess / ess)^tune$ESS_scale
      new_num_particles <- round(pm_settings$n_particles * scale_factor)
      pm_settings$n_particles <- max(25, min(tune$max_particles, new_num_particles))
    }
  }

  pm_settings$weight_ess <- weight_ess
  return(pm_settings)
}

# Step size of each state-centred component. The controller steers the share
# of the component's particles that outweigh the reference particle in this
# iteration to tune$step_target. For a Gaussian conditional with a matched
# covariance that share is about Phi(-epsilon * sqrt(d / 2)), so the default
# target Phi(-2.38 / sqrt(2)) puts epsilon at the random-walk scale
# 2.38 / sqrt(d). The share is taken per iteration: a running mean since the
# block started lags the step it responds to and overshoots. A component whose
# covariance changed (a new layout) restarts at the random-walk scale.
update_local_steps <- function(pm_settings, better, particle_numbers, tune,
                               n_pars, layout = NULL) {
  n_comp <- length(particle_numbers)
  eps <- rep_len(pm_settings$epsilon, max(n_comp, length(pm_settings$epsilon)))
  rw_scale <- 2.38 / sqrt(n_pars)
  if (is.null(layout)) {
    local <- rep(TRUE, n_comp)
    fresh <- rep(FALSE, n_comp)
  } else {
    local <- nzchar(layout)
    old <- pm_settings$step_layout
    if (length(old) != n_comp) old <- rep(NA_character_, n_comp)
    fresh <- local & (is.na(old) | old != layout)
    pm_settings$step_layout <- layout
  }
  target <- tune$step_target
  if (is.null(target)) target <- step_target(1)
  upd <- which(local & !fresh & particle_numbers > 0)
  if (length(upd)) {
    eps[upd] <- update_epsilon_continuous(
      epsilon    = eps[upd],
      acceptance = better[upd] / particle_numbers[upd],
      target     = target,
      iter       = pm_settings$iter,
      d          = n_pars,
      alphaStar  = -stats::qnorm(target / 2),
      damp       = 100,
      # A step below a quarter of the random-walk scale moves the chain
      # at under 1/16 of the optimal squared jump distance.
      clamp      = c(rw_scale / 4, 5)
    )
  }
  eps[which(fresh)] <- rw_scale
  pm_settings$epsilon <- eps
  pm_settings
}

step_target <- function(search_width) stats::pnorm(-2.38 / sqrt(2)) / search_width

update_epsilon_continuous <- function(
    epsilon,
    acceptance,
    target,
    iter,
    d,
    alphaStar,
    damp = 100,
    clamp
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
  # Plain multinomial allocation: the particles are then iid draws from the
  # mixture whose density weights them. Forcing a minimum of one particle per
  # component (the old pmax(1, .)) over-samples low-weight components and
  # biases the stationary distribution.
  return(as.vector(rmultinom(1, n_particles, mix_proportion)))
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
  static <- lapply(.EMC_STATIC_SAMPLE_FIELDS, function(nm) sampler$samples[[nm]])
  names(static) <- .EMC_STATIC_SAMPLE_FIELDS
  if(any(sampler$nuisance) && !is.null(sampler$sampler_nuis$samples)) {
    nuis_static <- lapply(.EMC_STATIC_SAMPLE_FIELDS,
                          function(nm) sampler$sampler_nuis$samples[[nm]])
    names(nuis_static) <- .EMC_STATIC_SAMPLE_FIELDS
    sampler$sampler_nuis$samples <- rapply(
      sampler$sampler_nuis$samples,
      f = function(x) extend_obj(x, n_samples, n_iter), how = "replace")
    for (nm in .EMC_STATIC_SAMPLE_FIELDS) {
      sampler$sampler_nuis$samples[[nm]] <- nuis_static[[nm]]
    }
  }
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples, n_iter), how = "replace")
  for (nm in .EMC_STATIC_SAMPLE_FIELDS) sampler$samples[[nm]] <- static[[nm]]
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


# Preburn's local step target (its first entry): the share of local particles
# that should beat the reference while searching. Later stages target the
# random-walk scale (step_target); mixture weights have no target
# (update_pm_settings).
set_p_accept <- function(stage, search_width){
  if(stage == "preburn") return(0.02 * (1/search_width))
  if(stage == "burn") return(c(0.02, 0.25)* (1/search_width))
  if(stage == "adapt") return(c(0.2, 0.25, 0.25)* (1/search_width))
  if(stage == "sample") return(c(0.3, 0.3, 0.3)* (1/search_width))
}

get_default_mix <- function(stage){
  if (stage == "burn") {
    default_mix <- c(0.15, 0.35, 0.5)
  } else if(stage == "adapt"){
    default_mix <- c(0.05, 0.3, 0.3, 0.35)
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
    # A stage can add more than one component (adapt now tunes the full
    # production mixture), so pad to the full length.
    epsilon <- c(epsilon, rep(epsilon[length(epsilon)],
                              length(mix) - 1 - length(epsilon)))
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
  if (use_pool && is.data.frame(dadm)) {
    # Pool workers hold the data of their template: send the current pooled
    # values as columns.
    pc <- .emc_pooled_columns(proposals, dadm, varying)
    proposals <- pc$proposals
    varying <- pc$varying
  }
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
  model <- .restore_custom_model_container(model, quiet = TRUE)
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component,
                                r_cores = r_cores, marginalise = marginalise,
                                varying = varying)
  } else{
    model <- model()
    dadm <- .cache_ll_data_attrs(dadm)
    pc <- .emc_pooled_columns(proposals, dadm, varying)
    proposals <- pc$proposals
    varying <- pc$varying
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
  model <- .restore_custom_model_container(model, quiet = TRUE)
  model <- model()
  dadm <- .cache_ll_data_attrs(dadm)
  c_name <- model$c_name
  unsupported <- c("MRI", "MRI_AR1", "SSEXG", "SSRDEX", "SOFTMAX")
  if(!is.null(c_name) && c_name %in% unsupported)
    stop("WAIC is not supported for model type '", c_name, "'")

  proposals <- .emc_pooled_columns(proposals, dadm)$proposals
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

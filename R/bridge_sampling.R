make_info <- function(samples, type){
  info <- list(
    subjects = samples$subjects,
    par_names = samples$par_names,
    n_pars = samples$n_pars,
    n_subjects = samples$n_subjects,
    prior = samples$prior,
    type = type,
    model = samples$model
  )
  bridge_add_info(info, samples, type)
}

# Each subject's block of alpha columns, named for the likelihood.
bridge_subject_proposals <- function(proposals, info) {
  n_pars <- info$n_pars
  lapply(seq_len(info$n_subjects), function(i) {
    out <- proposals[, (i - 1L) * n_pars + seq_len(n_pars), drop = FALSE]
    colnames(out) <- info$par_names
    out
  })
}

# Data log-likelihood of each proposal row, summed over subjects.
bridge_ll <- function(proposals_list, data, info) {
  lws <- mapply(calc_ll_manager, proposals_list, data,
                MoreArgs = list(model = info$model), SIMPLIFY = FALSE)
  rowSums(do.call(cbind, lws))
}

# Unnormalised log posterior of each proposal row. `lw` supplies the data
# log-likelihood when it is already known.
h.unnormalized.posterior <- function(proposals, data, info, hyper_only, lw = NULL) {
  proposals_list <- bridge_subject_proposals(proposals, info)
  if (hyper_only) {
    lw <- 0
  } else if (is.null(lw)) {
    lw <- bridge_ll(proposals_list, data, info)
  }
  proposals_group <- proposals[, info$group_idx, drop = FALSE]
  lw + bridge_group_and_prior_and_jac(proposals_group, proposals_list, info, info$type)
}

bridge_proposals <- function(kind, draws, iter_idx, gen_samples, m, R, rows) {
  switch(kind,
    posterior = draws[iter_idx[rows], , drop = FALSE],
    reflected = sweep(-draws[iter_idx[rows], , drop = FALSE], 2L, 2 * m, "+"),
    negative = sweep(-(gen_samples[rows, , drop = FALSE] %*% R), 2L, m, "+"),
    positive = sweep(gen_samples[rows, , drop = FALSE] %*% R, 2L, m, "+")
  )
}

# Unnormalised log posterior of every row of each proposal kind. The kinds
# are cut into row chunks that share one queue over `cores` processes, so an
# uneven kind does not leave a process idle. `lw` holds the data
# log-likelihood of the posterior kind's rows when it is already known.
bridge_eval_kinds <- function(kinds, draws, iter_idx, gen_samples, data, m, R, info,
                              cores, hyper_only, block_size, lw = NULL) {
  n <- length(iter_idx)
  chunk <- block_size
  if (cores > 1L)
    chunk <- min(chunk, max(.BRIDGE_MIN_CHUNK, ceiling(n * length(kinds) / (4 * cores))))
  jobs <- expand.grid(first = seq.int(1L, n, by = chunk), kind = kinds,
                      stringsAsFactors = FALSE)
  values <- .emc_rebuild_lapply(seq_len(nrow(jobs)), function(j) {
    rows <- seq.int(jobs$first[j], min(n, jobs$first[j] + chunk - 1L))
    kind <- jobs$kind[j]
    proposals <- bridge_proposals(kind, draws, iter_idx, gen_samples, m, R, rows)
    h.unnormalized.posterior(proposals, data, info, hyper_only,
                             lw = if (kind == "posterior" && !is.null(lw)) lw[rows])
  }, cores)
  out <- lapply(kinds, function(kind) unlist(values[jobs$kind == kind], use.names = FALSE))
  names(out) <- kinds
  out
}

# Fewest proposal rows worth a likelihood call of their own.
.BRIDGE_MIN_CHUNK <- 32L

run.iterative.scheme <- function(l1, l2, log_r0, tol, criterion, silent, maxiter, neff) {
  lstar <- median(l1)
  n.1 <- length(l1)
  n.2 <- length(l2)
  if (is.null(neff)) neff <- n.1

  l1 <- l1 - lstar
  l2 <- l2 - lstar
  log_s1 <- log(neff / (neff + n.2))
  log_s2 <- log(n.2 / (neff + n.2))
  log_r <- log_r0
  previous_log_r <- log_r
  logml <- log_r + lstar
  i <- 1L
  criterion_val <- 1 + tol

  while (criterion_val > tol && i <= maxiter) {
    if (!silent) cat("Iteration: ", i, "\n", sep = "")
    previous_log_r <- log_r
    previous_logml <- logml
    log_r <- bridge_update(l1, l2, log_s1, log_s2 + log_r)
    logml <- log_r + lstar
    i <- i + 1L

    criterion_val <- switch(criterion,
      r = abs(expm1(previous_log_r - log_r)),
      logml = abs((logml - previous_logml) / logml)
    )
    if (is.na(criterion_val)) break
  }

  if (i > maxiter || is.na(criterion_val)) logml <- NA_real_
  list(logml = logml, restart = (previous_log_r + log_r) / 2)
}

bridge_draws <- function(emc, stage, info) {
  samples <- emc[[1L]]
  stage_idx <- lapply(emc, function(chain) which(chain$samples$stage %in% stage))
  n <- sum(lengths(stage_idx))
  p <- info$n_pars
  draws <- matrix(NA_real_, n, p * info$n_subjects + length(info$group_idx))
  offset <- 0L
  for (chain in seq_along(emc)) {
    idx <- stage_idx[[chain]]
    if (!length(idx)) next
    samples$samples <- emc[[chain]]$samples
    rows <- offset + seq_along(idx)
    for (subject in seq_len(info$n_subjects)) {
      cols <- (subject - 1L) * p + seq_len(p)
      alpha <- samples$samples$alpha[, subject, idx, drop = FALSE]
      draws[rows, cols] <- t(matrix(alpha, nrow = p))
    }
    if (length(info$group_idx)) {
      group <- bridge_add_group(matrix(numeric(), length(idx), 0L), samples, idx, samples$type)
      draws[rows, info$group_idx] <- group
    }
    offset <- offset + length(idx)
  }
  draws
}

# The data log-likelihood the sampler stored with each draw, in the row order
# of bridge_draws(). NULL when it is missing or is not the likelihood bridge
# sampling evaluates, which a spread of draws from every chain is recomputed
# to establish.
bridge_stored_ll <- function(emc, stage, draws, data, info) {
  stored <- lapply(emc, function(chain) {
    ll <- chain$samples$subj_ll
    if (is.null(ll)) return(NULL)
    colSums(ll[, chain$samples$stage %in% stage, drop = FALSE])
  })
  if (any(vapply(stored, is.null, logical(1)))) return(NULL)
  chain_rows <- lengths(stored)
  stored <- unlist(stored, use.names = FALSE)
  if (length(stored) != nrow(draws) || !all(is.finite(stored))) return(NULL)
  first <- cumsum(c(0L, chain_rows[-length(chain_rows)]))
  rows <- sort(unique(unlist(lapply(seq_along(chain_rows), function(chain) {
    n <- chain_rows[chain]
    first[chain] + round(seq(1, n, length.out = min(n, .BRIDGE_LL_CHECKS)))
  }))))
  check <- bridge_ll(bridge_subject_proposals(draws[rows, , drop = FALSE], info), data, info)
  agree <- abs(check - stored[rows]) <= 1e-8 * (1 + abs(check))
  if (!all(is.finite(check)) || !all(agree)) return(NULL)
  stored
}

# Draws per chain whose stored log-likelihood is recomputed.
.BRIDGE_LL_CHECKS <- 8L

# Mean and Cholesky factor of the covariance of the fitting draws.
bridge_fit <- function(draws, fit_idx) {
  fit <- draws[fit_idx, , drop = FALSE]
  m <- colMeans(fit)
  V <- crossprod(sweep(fit, 2L, m)) / (nrow(fit) - 1L)
  # Factorised as a correlation matrix, so that columns on different scales
  # cannot make a full-rank covariance look singular.
  s <- sqrt(diag(V))
  R <- if (nrow(fit) > ncol(fit))
    tryCatch(chol(V / tcrossprod(s)) * rep(s, each = length(s)), error = function(e) NULL)
  # Only a covariance that cannot be factorised is moved to a nearby
  # positive-definite matrix: the move floors eigenvalues relative to the
  # largest, which narrows the proposal along well-determined directions.
  if (is.null(R)) R <- chol(Matrix::nearPD(V, base.matrix = TRUE)$mat)
  list(mean = m, chol = R)
}

# The terms of one estimate that use no random numbers: the proposal density
# at the estimation draws and the unnormalised posterior at those draws and at
# their reflections through the proposal mean.
bridge_fixed_terms <- function(draws, data, info, split_idx, fit, cores, hyper_only,
                               lw = NULL) {
  m <- fit$mean
  R <- fit$chol
  iter_idx <- seq_len(nrow(draws))[-split_idx]
  n_iter <- length(iter_idx)
  p <- ncol(draws)
  block_size <- bridge_block_size(p)
  density_cols <- if (hyper_only) info$group_idx else seq_len(p)
  log_const <- -0.5 * length(density_cols) * log(2 * pi)
  q12 <- numeric(n_iter)
  for (first in seq.int(1L, n_iter, by = block_size)) {
    rows <- seq.int(first, min(n_iter, first + block_size - 1L))
    centered <- sweep(draws[iter_idx[rows], , drop = FALSE], 2L, m, "-")
    standard <- backsolve(R, t(centered), transpose = TRUE)
    q12[rows] <- log_const - 0.5 * colSums(standard[density_cols, , drop = FALSE]^2)
  }
  rm(centered, standard)
  q <- bridge_eval_kinds(c("posterior", "reflected"), draws, iter_idx, NULL, data, m, R,
                         info, cores, hyper_only, block_size, lw = lw[iter_idx])
  list(q11 = bridge_log_add(q$posterior, q$reflected), q12 = q12)
}

# Rows per block, capping the temporary proposal and whitening matrices.
bridge_block_size <- function(p) max(1L, floor(2^21 / p))

# One bridge estimate of the log marginal likelihood: proposal draws are
# generated here, and `fixed` carries bridge_fixed_terms() when the caller
# already has them.
bridge_sampling <- function(draws, data, info, n_eff, split_idx, fit,
                            cores_for_props, cores_per_prop, maxiter, r0, tol1, tol2, hyper_only,
                            fixed = NULL) {
  cores <- max(1L, as.integer(cores_for_props)) * max(1L, as.integer(cores_per_prop))
  if (is.null(fixed))
    fixed <- bridge_fixed_terms(draws, data, info, split_idx, fit, cores, hyper_only)
  m <- fit$mean
  R <- fit$chol
  iter_idx <- seq_len(nrow(draws))[-split_idx]
  n_iter <- length(iter_idx)
  p <- ncol(draws)
  gen_samples <- matrix(NA_real_, n_iter, p)

  block_size <- bridge_block_size(p)
  q22 <- numeric(n_iter)
  density_cols <- if (hyper_only) info$group_idx else seq_len(p)
  log_const <- -0.5 * length(density_cols) * log(2 * pi)
  for (first in seq.int(1L, n_iter, by = block_size)) {
    rows <- seq.int(first, min(n_iter, first + block_size - 1L))
    generated <- matrix(rnorm(as.double(length(rows)) * p), length(rows), byrow = TRUE)
    gen_samples[rows, ] <- generated
    q22[rows] <- log_const - 0.5 * rowSums(generated[, density_cols, drop = FALSE]^2)
  }

  rm(generated)
  q <- bridge_eval_kinds(c("negative", "positive"), draws, iter_idx, gen_samples, data, m, R,
                         info, cores, hyper_only, block_size)
  log_scale <- sum(log(diag(R))) - log(2)
  l1 <- log_scale + fixed$q11 - fixed$q12
  l2 <- log_scale + bridge_log_add(q$negative, q$positive) - q22
  result <- run.iterative.scheme(l1, l2, log(r0), tol1, "r", TRUE, maxiter, n_eff)
  if (is.na(result$logml)) {
    result <- run.iterative.scheme(l1, l2, result$restart, tol2, "logml", TRUE, maxiter, n_eff)
  }

  if (is.na(result$logml)) {
    stop("Bridge sampling did not converge, usually this means you need to fit longer to get more samples")
  }
  result$logml
}

#' Estimating Marginal Likelihoods Using WARP-III Bridge Sampling
#'
#' Uses bridge sampling that matches a proposal distribution to the first three moments
#' of the posterior distribution to get an accurate estimate of the marginal likelihood.
#' The marginal likelihood can be used for computing Bayes factors and posterior model probabilities.
#'
#'
#' If not enough posterior samples were collected using `fit()`,
#' bridge sampling can be unstable. It is recommended to run
#' `run_bridge_sampling()` several times with the ``repetitions`` argument
#'  and to examine how stable the results are.
#'
#' It can be difficult to converge bridge sampling for exceptionally large models,
#' because of a large number of subjects (> 100) and/or cognitive model parameters.
#'
#' For a practical introduction:
#'
#' Gronau, Q. F., Heathcote, A., & Matzke, D. (2020). Computing Bayes factors
#' for evidence-accumulation models using Warp-III bridge sampling.
#' *Behavior research methods*, 52(2), 918-937. doi.org/10.3758/s13428-019-01290-6
#'
#' For mathematical background:
#'
#' Meng, X.-L., & Wong, W. H. (1996). Simulating ratios of normalizing
#' constants via a simple identity: A theoretical exploration. *Statistica Sinica*,
#' 6, 831-860. http://www3.stat.sinica.edu.tw/statistica/j6n4/j6n43/j6n43.htm
#'
#' Meng, X.-L., & Schilling, S. (2002). Warp bridge sampling.
#' *Journal of Computational and Graphical Statistics*,
#' 11(3), 552-586. doi.org/10.1198/106186002457
#'
#' @param emc An emc object with a set of converged samples
#' @param stage A character indicating which stage to use, defaults to `sample`
#' @param filter An integer or vector. If integer, it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param repetitions An integer. How many times to repeat the bridge sampling scheme. Can help get an estimate of stability of the estimate.
#' @param cores_for_props Integer. Warp-III evaluates the posterior over 4 different proposal densities. Their draws are split into chunks of rows that are evaluated in parallel on ``cores_for_props`` * ``cores_per_prop`` cores.
#' @param cores_per_prop Integer. Multiplies ``cores_for_props``: the evaluations use ``cores_for_props`` * ``cores_per_prop`` cores in total. Must be 1 on Windows.
#' @param both_splits Boolean. Bridge sampling uses a proposal density and a target density. We can estimate the stability of our samples and therefore MLL estimate, by running 2 bridge sampling iterations
#' The first one uses the first half of the samples as the proposal and the second half as the target, the second run uses the opposite. If this is is set to ``FALSE``, it will only run bridge sampling once and
#' it will instead do an odd-even iterations split to get a more reasonable estimate for just one run.
#' @param ... Additional, optional more in-depth hyperparameters
#'
#' @return A vector of length repetitions which contains the marginal log likelihood estimates per repetition
#' @examples \donttest{
#' # After `fit` has converged on a specific model
#' # We can take those samples and calculate the marginal log-likelihood for them
#' MLL <- run_bridge_sampling(samples_LNR, cores_for_props = 1, both_splits = FALSE)
#' # This will run on 2*4 cores (since 4 is the default for ``cores_for_props``)
#' }
#' @export
#'
run_bridge_sampling <- function(emc, stage = "sample", filter = NULL, repetitions = 1, cores_for_props = 4,  cores_per_prop = 1, both_splits = TRUE, ...){
  emc <- restore_custom_kernel_pointers(emc, quiet = TRUE)
  if (!is.null(emc[[1]]$pooled))
    stop("Bridge sampling does not yet support pooled parameters")
  # Hyper parameters and dev options
  maxiter <- 5000
  r0 <- 1e-5
  tol1 <- 1e-10
  tol2 <- 1e-6
  hyper_only <- F
  # overwrite those that were supplied
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if (.Platform$OS.type == "windows" && cores_per_prop > 1)
    stop("only cores_for_props can be set on Windows")
  cores <- max(1L, as.integer(cores_for_props)) * max(1L, as.integer(cores_per_prop))
  # The stage and filter are applied here, once, for everything below.
  emc <- subset(emc, filter = filter, stage = stage)
  n_eff <- round(ess_summary(emc, selection = "alpha", stat = "median", stat_only = TRUE,
                             stage = stage, n_cores = cores) / 2)
  info <- make_info(emc[[1L]], emc[[1L]]$type)
  data <- emc[[1L]]$data
  draws <- bridge_draws(emc, stage, info)
  lw <- if (!hyper_only && repetitions) bridge_stored_ll(emc, stage, draws, data, info)
  rm(emc)
  n <- nrow(draws)
  if (!repetitions) return(numeric())
  if (n < 4L) stop("Bridge sampling requires at least four draws from the chosen stage")
  split1 <- if (both_splits) seq_len(round(n / 2)) else seq.int(1L, n, by = 2L)
  splits <- list(split1)
  fit_idx <- list(if (both_splits) split1 else head(split1, n - length(split1)))
  if (both_splits) splits[[2L]] <- fit_idx[[2L]] <- seq.int(length(split1) + 1L, n)
  # The fits and the terms that use no random numbers are shared by every
  # repetition.
  fits <- lapply(fit_idx, function(idx) bridge_fit(draws, idx))
  fixed <- lapply(seq_along(splits), function(s)
    bridge_fixed_terms(draws, data, info, splits[[s]], fits[[s]], cores, hyper_only, lw))
  mls <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    estimates <- vapply(seq_along(splits), function(s)
      bridge_sampling(draws, data, info, n_eff, splits[[s]], fits[[s]], cores_for_props,
                      cores_per_prop, maxiter, r0, tol1, tol2, hyper_only, fixed = fixed[[s]]),
      numeric(1))
    if (both_splits && abs(estimates[1L] - estimates[2L]) > 1) warning("First split and second split marginal likelihood estimates are off by 1 log point. \n This usually means that your MCMC chains aren't completely stable yet. \n Consider running the MCMC chain longer if you need more precise estimates (e.g. when comparing different priors)")
    mls[i] <- mean(estimates)
  }
  mls
}

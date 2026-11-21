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

bridge_eval_proposal <- function(kind, draws, iter_idx, gen_samples, data, m, R, info,
                                 cores_per_prop, hyper_only, block_size) {
  n <- length(iter_idx)
  out <- numeric(n)
  for (first in seq.int(1L, n, by = block_size)) {
    rows <- seq.int(first, min(n, first + block_size - 1L))
    proposals <- switch(kind,
      posterior = draws[iter_idx[rows], , drop = FALSE],
      reflected = sweep(-draws[iter_idx[rows], , drop = FALSE], 2L, 2 * m, "+"),
      negative = sweep(-(gen_samples[rows, , drop = FALSE] %*% R), 2L, m, "+"),
      positive = sweep(gen_samples[rows, , drop = FALSE] %*% R, 2L, m, "+")
    )
    out[rows] <- h.unnormalized.posterior(proposals, data, info, cores_per_prop, hyper_only)
  }
  out
}

eval.unnormalized.posterior <- function(draws, iter_idx, gen_samples, data, m, R, info,
                                        cores_for_props, cores_per_prop, hyper_only, block_size) {
  mls <- auto_mclapply(c("posterior", "reflected", "negative", "positive"), bridge_eval_proposal,
                       draws = draws, iter_idx = iter_idx, gen_samples = gen_samples, data = data,
                       m = m, R = R, info = info, cores_per_prop = cores_per_prop,
                       hyper_only = hyper_only, block_size = block_size, mc.cores = cores_for_props)
  list(q11 = bridge_log_add(mls[[1L]], mls[[2L]]),
       q21 = bridge_log_add(mls[[3L]], mls[[4L]]))
}

h.unnormalized.posterior <- function(proposals, data, info, n_cores, hyper_only) {
  n_pars <- info$n_pars
  n_subjects <- info$n_subjects
  proposals_list <- vector("list", length = n_subjects)
  for(i in seq_len(n_subjects)){
    tmp <- proposals[, ((i-1)*n_pars + 1):(i*n_pars), drop = FALSE]
    colnames(tmp) <- info$par_names
    proposals_list[[i]] <- tmp
  }
  if(hyper_only){
    lw <- 0
  } else{
    if(n_cores <= 1) {
      lws <- mapply(calc_ll_manager, proposals_list, data,
                    MoreArgs = list(model = info$model), SIMPLIFY = FALSE)
    } else {
      lws <- parallel::mcmapply(.emc_armed(calc_ll_manager), proposals_list, data,
                                MoreArgs = list(model = info$model),
                                mc.cores = n_cores, SIMPLIFY = FALSE)
    }
    lw <- rowSums(do.call(cbind, lws))
  }
  proposals_group <- proposals[, info$group_idx, drop = FALSE]
  gr_pr_jac <- bridge_group_and_prior_and_jac(proposals_group, proposals_list, info, info$type)
  return(lw + gr_pr_jac)
}

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

bridge_fit <- function(draws, fit_idx) {
  fit <- draws[fit_idx, , drop = FALSE]
  list(mean = colMeans(fit), chol = chol(Matrix::nearPD(var(fit), base.matrix = TRUE)$mat))
}

bridge_sampling <- function(draws, data, info, n_eff, split_idx, fit,
                            cores_for_props, cores_per_prop, maxiter, r0, tol1, tol2, hyper_only) {
  m <- fit$mean
  R <- fit$chol
  iter_idx <- seq_len(nrow(draws))[-split_idx]
  n_iter <- length(iter_idx)
  p <- ncol(draws)
  gen_samples <- matrix(NA_real_, n_iter, p)

  # Cap the size of temporary proposal and whitening matrices per worker.
  block_size <- max(1L, floor(2^21 / p))
  q12 <- q22 <- numeric(n_iter)
  density_cols <- if (hyper_only) info$group_idx else seq_len(p)
  log_const <- -0.5 * length(density_cols) * log(2 * pi)
  for (first in seq.int(1L, n_iter, by = block_size)) {
    rows <- seq.int(first, min(n_iter, first + block_size - 1L))
    generated <- matrix(rnorm(as.double(length(rows)) * p), length(rows), byrow = TRUE)
    gen_samples[rows, ] <- generated
    centered <- sweep(draws[iter_idx[rows], , drop = FALSE], 2L, m, "-")
    standard <- backsolve(R, t(centered), transpose = TRUE)
    q12[rows] <- log_const - 0.5 * colSums(standard[density_cols, , drop = FALSE]^2)
    q22[rows] <- log_const - 0.5 * rowSums(generated[, density_cols, drop = FALSE]^2)
  }

  rm(centered, standard, generated)
  q <- eval.unnormalized.posterior(draws, iter_idx, gen_samples, data, m, R, info,
                                   cores_for_props, cores_per_prop, hyper_only, block_size)
  log_scale <- sum(log(diag(R))) - log(2)
  l1 <- log_scale + q$q11 - q12
  l2 <- log_scale + q$q21 - q22
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
#' @param cores_for_props Integer. Warp-III evaluates the posterior over 4 different proposal densities. If you have the CPU, 4 cores will do this in parallel, 2 is also already helpful.
#' @param cores_per_prop Integer. Per density we can also parallelize across subjects. Eventual cores will be ``cores_for_props`` * ``cores_per_prop``. For efficiency users should prioritize cores_for_props being 4.
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
  emc <-subset(emc, filter = filter, stage = stage)
  n_eff <- round(ess_summary(emc, selection = "alpha", stat = "median", stat_only = TRUE, stage = stage, filter = filter)/2)
  if (.Platform$OS.type == "windows" && cores_per_prop > 1)
    stop("only cores_for_props can be set on Windows")
  info <- make_info(emc[[1L]], emc[[1L]]$type)
  data <- emc[[1L]]$data
  draws <- bridge_draws(emc, stage, info)
  rm(emc)
  n <- nrow(draws)
  if (!repetitions) return(numeric())
  if (n < 4L) stop("Bridge sampling requires at least four draws from the chosen stage")
  split1 <- if (both_splits) seq_len(round(n / 2)) else seq.int(1L, n, by = 2L)
  split2 <- if (both_splits) seq.int(length(split1) + 1L, n)
  fit_idx <- if (both_splits) split1 else head(split1, n - length(split1))
  fit1 <- bridge_fit(draws, fit_idx)
  mls <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    s1 <- bridge_sampling(draws, data, info, n_eff, split1, fit1, cores_for_props, cores_per_prop,
                          maxiter, r0, tol1, tol2, hyper_only)
    if (both_splits) {
      fit2 <- bridge_fit(draws, split2)
      s2 <- bridge_sampling(draws, data, info, n_eff, split2, fit2, cores_for_props, cores_per_prop,
                            maxiter, r0, tol1, tol2, hyper_only)
      rm(fit2)
      if (abs(s1 - s2) > 1) warning("First split and second split marginal likelihood estimates are off by 1 log point. \n This usually means that your MCMC chains aren't completely stable yet. \n Consider running the MCMC chain longer if you need more precise estimates (e.g. when comparing different priors)")
      s1 <- mean(c(s1, s2))
    }
    mls[i] <- s1
  }
  mls
}

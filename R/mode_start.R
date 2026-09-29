.emc_burn_subject_prior <- function(sampler, s) {
  store <- sampler$samples
  i <- store$idx
  p <- sampler$n_pars
  type <- sampler$type
  if (type == "single") {
    return(list(mu = sampler$prior$theta_mu_mean,
                var = sampler$prior$theta_mu_var))
  }
  if (type == "SEM") {
    last <- last_sample_SEM(store)
    pars <- list(tmu = last$mu, lambda = last$lambda, eta = last$eta,
                 B = last$B, K = last$K, G = last$G,
                 delta_inv = last$delta_inv, epsilon_inv = last$epsilon_inv)
    prior <- .emc_ensemble_subject_prior(pars, sampler,
                                         sampler$n_subjects, "marginal")
    return(list(mu = prior$mean[, s], var = prior$var))
  }
  mu <- store$theta_mu[, i]
  if (type == "standard") {
    group_designs <- sampler$gd
    if (is.null(group_designs)) {
      group_designs <- add_group_design(
        sampler$par_names[!sampler$nuisance], sampler$group_designs,
        sampler$n_subjects)
    }
    mu <- calculate_subject_means(group_designs, mu)[, s]
  }
  list(mu = as.numeric(mu),
       var = matrix(store$theta_var[, , i], p, p))
}

.emc_burn_ll_batch <- function(x, data, model, par_names, batch_size = 256L) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  out <- numeric(nrow(x))
  for (first in seq.int(1L, nrow(x), by = batch_size)) {
    last <- min(first + batch_size - 1L, nrow(x))
    proposals <- x[first:last, , drop = FALSE]
    colnames(proposals) <- par_names
    out[first:last] <- as.numeric(calc_ll_manager(
      proposals, dadm = data, model = model, r_cores = 1L))
  }
  out
}

.emc_burn_hessian <- function(x, prior_var, precision, ll_batch, step = NULL) {
  # Mixed differences with a step wider than the posterior along a sharply
  # identified direction come out indefinite, and shrinking them toward the
  # diagonal then collapses the covariance along every ridge. Try smaller
  # steps first; shrink only if no step gives a positive definite Hessian.
  if (is.null(step)) {
    # Gradient-difference Hessian (stats::optimHess) first: it stays positive
    # definite where the four-point mixed stencil does not.
    ll_hess <- tryCatch(stats::optimHess(x, function(z) {
      v <- ll_batch(z)
      if (length(v) != 1L || !is.finite(v)) stop("nonfinite likelihood")
      -v
    }), error = function(e) NULL)
    if (is.matrix(ll_hess) && all(is.finite(ll_hess))) {
      out <- .emc_hessian_covariance(precision + (ll_hess + t(ll_hess)) / 2,
                                     precision, prior_var, allow_shrink = FALSE)
      if (is.matrix(out)) return(out)
    }
    sd <- sqrt(diag(prior_var))
    for (m in c(0.003, 0.001, 3e-4, 1e-4)) {
      out <- .emc_burn_hessian_step(x, prior_var, precision, ll_batch,
                                    pmax(1e-6, m * sd), allow_shrink = FALSE)
      if (is.matrix(out)) return(out)
    }
    return(.emc_burn_hessian_step(x, prior_var, precision, ll_batch,
                                  pmax(1e-6, 0.003 * sd), allow_shrink = TRUE))
  }
  .emc_burn_hessian_step(x, prior_var, precision, ll_batch, step,
                         allow_shrink = TRUE)
}

.emc_burn_hessian_step <- function(x, prior_var, precision, ll_batch, step,
                                   allow_shrink = TRUE) {
  p <- length(x)
  # The prior can be much wider than a wide-design subject posterior. A step
  # of 1% of its SD can cross a sharp likelihood ridge and create spurious
  # negative cross-curvature even at a converged optimum.
  pairs <- which(upper.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
  n_eval <- 1L + 2L * p + 4L * nrow(pairs)
  proposals <- matrix(rep(x, each = n_eval), nrow = n_eval)
  row <- 1L
  for (k in seq_len(p)) {
    proposals[row + 1L, k] <- x[k] + step[k]
    proposals[row + 2L, k] <- x[k] - step[k]
    row <- row + 2L
  }
  if (nrow(pairs)) for (k in seq_len(nrow(pairs))) {
    a <- pairs[k, 1L]
    b <- pairs[k, 2L]
    proposals[row + 1L, c(a, b)] <- x[c(a, b)] + c(step[a], step[b])
    proposals[row + 2L, c(a, b)] <- x[c(a, b)] + c(step[a], -step[b])
    proposals[row + 3L, c(a, b)] <- x[c(a, b)] + c(-step[a], step[b])
    proposals[row + 4L, c(a, b)] <- x[c(a, b)] - c(step[a], step[b])
    row <- row + 4L
  }
  ll <- ll_batch(proposals)
  if (length(ll) != n_eval || any(!is.finite(ll))) return(NULL)
  hessian <- precision
  row <- 1L
  for (k in seq_len(p)) {
    hessian[k, k] <- hessian[k, k] -
      (ll[row + 1L] - 2 * ll[1L] + ll[row + 2L]) / step[k]^2
    row <- row + 2L
  }
  if (nrow(pairs)) for (k in seq_len(nrow(pairs))) {
    a <- pairs[k, 1L]
    b <- pairs[k, 2L]
    curvature <- -(ll[row + 1L] - ll[row + 2L] - ll[row + 3L] +
                     ll[row + 4L]) / (4 * step[a] * step[b])
    hessian[a, b] <- hessian[a, b] + curvature
    hessian[b, a] <- hessian[b, a] + curvature
    row <- row + 4L
  }
  .emc_hessian_covariance(hessian, precision, prior_var, allow_shrink)
}

.emc_hessian_covariance <- function(hessian, precision, prior_var,
                                    allow_shrink = TRUE) {
  p <- nrow(hessian)
  root <- t(chol(prior_var))
  # Mixed differences can still be unreliable near a likelihood boundary.
  # Shrink only those interactions toward the positive coordinate curvatures;
  # the Gaussian prior precision remains intact. This changes the proposal,
  # never the target density or its acceptance weights.
  likelihood_diag <- pmax(diag(hessian - precision), 0)
  reference <- precision + diag(likelihood_diag, p)
  eig <- NULL
  for (shrink in if (allow_shrink) c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1) else 0) {
    candidate <- (1 - shrink) * hessian + shrink * reference
    scaled <- crossprod(root, candidate %*% root)
    attempt <- tryCatch(eigen((scaled + t(scaled)) / 2, symmetric = TRUE),
                        error = function(e) NULL)
    if (!is.null(attempt) && all(is.finite(attempt$values)) &&
        min(attempt$values) > 1e-6) {
      eig <- attempt
      break
    }
  }
  if (is.null(eig)) return(NULL)
  inverse <- eig$vectors %*%
    ((1 / pmax(eig$values, 0.1)) * t(eig$vectors))
  covariance <- root %*% inverse %*% t(root)
  covariance <- (covariance + t(covariance)) / 2
  if (!.emc_empirical_covariance_ok(covariance))
    covariance <- .emc_regularize_covariance(covariance)
  if (!.emc_empirical_covariance_ok(covariance)) return(NULL)
  covariance
}

.emc_optimize_subject <- function(x0, prior, ll_batch, maxit = 2000L,
                                  proposal_prior_scale = 1) {
  p <- length(x0)
  root <- tryCatch(chol(prior$var), error = function(e) NULL)
  if (is.null(root)) return(list(reason = "prior covariance"))
  precision <- chol2inv(root)
  best <- new.env(parent = emptyenv())
  best$value <- Inf
  best$par <- x0
  objective <- function(x) {
    ll <- tryCatch(ll_batch(x), error = function(e) NA_real_)
    if (length(ll) != 1L || !is.finite(ll)) return(1e10)
    delta <- x - prior$mu
    value <- -ll + 0.5 * drop(crossprod(delta, precision %*% delta))
    if (is.finite(value) && value < best$value) {
      best$value <- value
      best$par <- x
    }
    value
  }
  f0 <- objective(x0)
  if (!is.finite(f0) || f0 >= 1e10) return(list(reason = "initial likelihood"))
  grad_step <- pmax(1e-5, 1e-3 * sqrt(diag(prior$var)))
  gradient <- function(x) {
    f <- objective(x)
    if (!is.finite(f) || f >= 1e10) stop("nonfinite central objective")
    steps <- grad_step
    ll_gradient <- rep(NA_real_, p)
    for (attempt in seq_len(5L)) {
      pending <- which(!is.finite(ll_gradient))
      if (!length(pending)) break
      proposals <- matrix(rep(x, each = 2L * length(pending)),
                          nrow = 2L * length(pending))
      for (i in seq_along(pending)) {
        k <- pending[i]
        proposals[2L * i - 1L, k] <- x[k] + steps[k]
        proposals[2L * i, k] <- x[k] - steps[k]
      }
      ll <- tryCatch(ll_batch(proposals), error = function(e)
        rep(NA_real_, nrow(proposals)))
      center_ll <- tryCatch(ll_batch(x), error = function(e) NA_real_)
      for (i in seq_along(pending)) {
        k <- pending[i]
        plus <- ll[2L * i - 1L]
        minus <- ll[2L * i]
        if (is.finite(plus) && is.finite(minus))
          ll_gradient[k] <- (plus - minus) / (2 * steps[k])
        else if (is.finite(plus))
          ll_gradient[k] <- (plus - center_ll) / steps[k]
        else if (is.finite(minus))
          ll_gradient[k] <- (center_ll - minus) / steps[k]
        else steps[k] <- steps[k] / 10
      }
    }
    if (any(!is.finite(ll_gradient))) stop("nonfinite likelihood gradient")
    -ll_gradient + drop(precision %*% (x - prior$mu))
  }
  opt <- tryCatch(stats::optim(x0, objective, gradient, method = "BFGS",
                               control = list(maxit = maxit, reltol = 1e-7)),
                  error = function(e) e)
  if (is.list(opt) && !inherits(opt, "error") &&
      all(is.finite(opt$par))) objective(opt$par)
  x <- best$par
  ll <- tryCatch(ll_batch(x), error = function(e) NA_real_)
  converged <- is.list(opt) && !inherits(opt, "error") &&
    identical(opt$convergence, 0L)
  grad <- tryCatch(gradient(x), error = function(e) rep(Inf, p))
  stationary <- all(is.finite(grad)) &&
    max(abs(grad * sqrt(diag(prior$var)))) < 0.1
  # The mode is found under the actual group prior. The proposal covariance
  # can use a wider one: directions the likelihood barely informs are then not
  # tied to the current group variance, which in a centred hierarchy would
  # otherwise let a shrinking group SD shrink the proposals with it.
  covariance <- if (p <= 40L && (converged || stationary)) {
    tryCatch(.emc_burn_hessian(x, prior$var * proposal_prior_scale,
                               precision / proposal_prior_scale, ll_batch),
             error = function(e) NULL)
  } else NULL
  list(alpha = x, ll = ll, covariance = covariance,
       gain = f0 - best$value, converged = converged,
       stationary = stationary,
       reason = if (!is.null(covariance)) "ok" else if (p > 40L) {
         "dimension limit"
       } else if (!converged && !stationary) "not stationary" else "Hessian")
}

.emc_burn_subject_mode <- function(sampler, s, proposal_prior_scale = 1) {
  store <- sampler$samples
  i <- store$idx
  p <- sampler$n_pars
  x0 <- as.numeric(store$alpha[, s, i])
  prior <- tryCatch(.emc_burn_subject_prior(sampler, s), error = function(e) NULL)
  if (is.null(prior) || length(prior$mu) != p ||
      any(!is.finite(prior$mu)) || any(!is.finite(prior$var))) return(NULL)
  par_names <- rownames(store$alpha)
  ll_batch <- function(x) .emc_burn_ll_batch(
    x, sampler$data[[s]], sampler$model, par_names)
  .emc_optimize_subject(x0, prior, ll_batch,
                        proposal_prior_scale = proposal_prior_scale)
}

# Optimisation probes the likelihood far from the posterior, where it can warn.
# In a forked child those conditions reach any handler inherited from the
# parent (testthat's parallel reporter writes them to its pipe from the child
# and corrupts the stream), and they carry no information here.
.emc_quiet_fork <- function(expr) suppressWarnings(suppressMessages(expr))

.emc_burn_mode_start <- function(emc, n_cores, verbose = FALSE) {
  if (!isTRUE(getOption("emc2.burn_mode_start", TRUE))) return(emc)
  if (any(emc[[1L]]$nuisance) || !is.null(emc[[1L]]$marginalise)) return(emc)
  jobs <- expand.grid(chain = seq_along(emc),
                      subject = seq_len(emc[[1L]]$n_subjects))
  results <- .emc_with_preserved_rng(auto_mclapply(
    seq_len(nrow(jobs)), function(k) {
      j <- jobs$chain[k]
      s <- jobs$subject[k]
      sampler <- emc[[j]]
      assign(".Random.seed", sampler$rng$subjects[[s]], envir = globalenv())
      result <- .emc_quiet_fork(.emc_burn_subject_mode(sampler, s))
      if (is.list(result) && is.matrix(result$covariance)) {
        # Spread the chains across the local posterior neighbourhood. Reject
        # proposals outside the likelihood support, keeping the mode as backup.
        root <- chol(3 * result$covariance)
        for (attempt in seq_len(8L)) {
          draw <- as.numeric(result$alpha + drop(stats::rnorm(length(result$alpha)) %*% root))
          ll <- tryCatch(.emc_burn_ll_batch(
            draw, sampler$data[[s]], sampler$model,
            rownames(sampler$samples$alpha)), error = function(e) NA_real_)
          if (length(ll) == 1L && is.finite(ll)) {
            result$start <- draw
            result$start_ll <- ll
            break
          }
        }
        if (is.null(result$start)) {
          result$start <- result$alpha
          result$start_ll <- result$ll
        }
      }
      list(result = result, seed = get(".Random.seed", envir = globalenv()))
    }, mc.cores = min(nrow(jobs), n_cores), mc.set.seed = FALSE))
  for (j in seq_along(emc)) {
    emc[[j]]$burn_mode_stats <- list(
      gain = rep(NA_real_, emc[[j]]$n_subjects),
      converged = rep(FALSE, emc[[j]]$n_subjects),
      hessian = rep(FALSE, emc[[j]]$n_subjects),
      reason = rep("optimizer error", emc[[j]]$n_subjects))
  }
  for (k in seq_len(nrow(jobs))) {
    job_result <- results[[k]]
    if (!is.list(job_result) || inherits(job_result, "try-error")) next
    j <- jobs$chain[k]
    s <- jobs$subject[k]
    emc[[j]]$rng$subjects[[s]] <- job_result$seed
    result <- job_result$result
    if (!is.list(result)) next
    if (!is.null(result$reason)) emc[[j]]$burn_mode_stats$reason[s] <- result$reason
    if (!is.null(result$gain)) emc[[j]]$burn_mode_stats$gain[s] <- result$gain
    if (!is.null(result$converged))
      emc[[j]]$burn_mode_stats$converged[s] <- result$converged
    if (is.null(result$covariance) || !is.matrix(result$covariance)) next
    i <- emc[[j]]$samples$idx
    emc[[j]]$samples$alpha[, s, i] <- result$start
    emc[[j]]$samples$subj_ll[s, i] <- result$start_ll
    emc[[j]]$burn_mode_stats$hessian[s] <- TRUE
    emc[[j]]$burn_mode_index <- i
    if (is.null(emc[[j]]$burn_mode_var))
      emc[[j]]$burn_mode_var <- vector("list", emc[[j]]$n_subjects)
    if (is.null(emc[[j]]$burn_mode_mu))
      emc[[j]]$burn_mode_mu <- vector("list", emc[[j]]$n_subjects)
    emc[[j]]$burn_mode_var[s] <- list(result$covariance)
    emc[[j]]$burn_mode_mu[s] <- list(result$alpha)
  }
  stats <- lapply(emc, `[[`, "burn_mode_stats")
  n_hessian <- sum(vapply(stats, function(x) sum(x$hessian), integer(1)))
  n_converged <- sum(vapply(stats, function(x) sum(x$converged), integer(1)))
  reasons <- table(unlist(lapply(stats, `[[`, "reason")))
  failures <- reasons[names(reasons) != "ok"]
  if (verbose || n_hessian < nrow(jobs)) {
    message("Burn mode starts: ", n_hessian, "/", nrow(jobs),
            " Hessians; ", n_converged, "/", nrow(jobs),
            " BFGS converged", if (n_hessian < nrow(jobs))
              paste0("; failures: ", paste(names(failures), failures,
                                           collapse = ", ")) else "")
  }
  emc
}

# Laplace proposal metric for the discarded stages: at each proposal refresh,
# find every chain's conditional mode for each subject under that chain's
# current group state and store the Laplace covariance (and mode). The chain
# state is not moved. create_chain_proposals() consumes these once.
.emc_laplace_refresh <- function(emc, n_cores, verbose = FALSE) {
  if (!isTRUE(getOption("emc2.laplace_refresh", TRUE))) return(emc)
  if (any(emc[[1L]]$nuisance) || !is.null(emc[[1L]]$marginalise)) return(emc)
  jobs <- expand.grid(chain = seq_along(emc),
                      subject = seq_len(emc[[1L]]$n_subjects))
  results <- .emc_with_preserved_rng(auto_mclapply(
    seq_len(nrow(jobs)), function(k) {
      tryCatch(.emc_quiet_fork(
        .emc_burn_subject_mode(emc[[jobs$chain[k]]], jobs$subject[k],
                               proposal_prior_scale = 4)),
        error = function(e) NULL)
    }, mc.cores = min(nrow(jobs), n_cores), mc.set.seed = FALSE))
  n_ok <- 0L
  for (k in seq_len(nrow(jobs))) {
    result <- results[[k]]
    if (!is.list(result) || !is.matrix(result$covariance)) next
    j <- jobs$chain[k]
    s <- jobs$subject[k]
    if (is.null(emc[[j]]$burn_mode_var))
      emc[[j]]$burn_mode_var <- vector("list", emc[[j]]$n_subjects)
    if (is.null(emc[[j]]$burn_mode_mu))
      emc[[j]]$burn_mode_mu <- vector("list", emc[[j]]$n_subjects)
    emc[[j]]$burn_mode_var[s] <- list(result$covariance)
    emc[[j]]$burn_mode_mu[s] <- list(result$alpha)
    n_ok <- n_ok + 1L
  }
  if (verbose || n_ok < nrow(jobs))
    message("Laplace refresh: ", n_ok, "/", nrow(jobs), " covariances")
  emc
}

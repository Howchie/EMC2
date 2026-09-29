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

.emc_burn_ll_batch <- function(x, data, model, par_names, batch_size = 256L,
                               component = NULL) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  out <- numeric(nrow(x))
  for (first in seq.int(1L, nrow(x), by = batch_size)) {
    last <- min(first + batch_size - 1L, nrow(x))
    proposals <- x[first:last, , drop = FALSE]
    colnames(proposals) <- par_names
    out[first:last] <- as.numeric(calc_ll_manager(
      proposals, dadm = data, model = model, component = component,
      r_cores = 1L))
  }
  out
}

# Likelihood Hessian (of -ll) assembled from separable blocks: a joint model's
# likelihood is a sum of per-model terms, each depending only on its own
# parameters, so the cross-model entries are exactly zero and each block needs
# only its own model's likelihood. blocks = NULL is the whole vector at once.
.emc_blockwise_ll_hessian <- function(x, ll_batch, blocks, fn, step = NULL) {
  if (is.null(blocks)) return(fn(x, ll_batch, step))
  out <- matrix(0, length(x), length(x))
  for (b in blocks) {
    h <- fn(x[b$idx], b$ll, step[b$idx])
    if (!is.matrix(h)) return(NULL)
    out[b$idx, b$idx] <- h
  }
  out
}

.emc_burn_hessian <- function(x, prior_var, precision, ll_batch, step = NULL,
                              blocks = NULL) {
  # Mixed differences with a step wider than the posterior along a sharply
  # identified direction come out indefinite, and shrinking them toward the
  # diagonal then collapses the covariance along every ridge. Try smaller
  # steps first; shrink only if no step gives a positive definite Hessian.
  if (is.null(step)) {
    # Gradient-difference Hessian first: it stays positive definite where the
    # four-point mixed stencil does not.
    ll_hess <- tryCatch(.emc_blockwise_ll_hessian(
      x, ll_batch, blocks,
      function(x, ll_batch, step) .emc_gradient_difference_hessian(x, ll_batch)),
      error = function(e) NULL)
    if (is.matrix(ll_hess) && all(is.finite(ll_hess))) {
      out <- .emc_hessian_covariance(precision + (ll_hess + t(ll_hess)) / 2,
                                     precision, prior_var, allow_shrink = FALSE)
      if (is.matrix(out)) return(out)
    }
    sd <- sqrt(diag(prior_var))
    for (m in c(0.003, 0.001, 3e-4, 1e-4)) {
      out <- .emc_burn_hessian_step(x, prior_var, precision, ll_batch,
                                    pmax(1e-6, m * sd), allow_shrink = FALSE,
                                    blocks = blocks)
      if (is.matrix(out)) return(out)
    }
    return(.emc_burn_hessian_step(x, prior_var, precision, ll_batch,
                                  pmax(1e-6, 0.003 * sd), allow_shrink = TRUE,
                                  blocks = blocks))
  }
  .emc_burn_hessian_step(x, prior_var, precision, ll_batch, step,
                         allow_shrink = TRUE, blocks = blocks)
}

# Hessian of -ll as the central difference of a central-difference gradient,
# the construction stats::optimHess() uses without an analytic gradient (both
# steps eps). optimHess() evaluates 4p^2 points one likelihood call at a time.
# Here they go through ll_batch, and each off-diagonal point is evaluated once:
# x + a e_i + b e_j (gradient j at outer step i) and x + b e_j + a e_i (gradient
# i at outer step j) coincide, leaving 2p^2 + 2p points. Returns NULL if any
# point is nonfinite.
.emc_gradient_difference_hessian <- function(x, ll_batch, eps = 1e-3,
                                             max_rows = 4096L) {
  p <- length(x)
  up <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)   # i <= j
  # Four points per pair, (a, b) = (+, +), (+, -), (-, +), (-, -), at
  # x + a e_i + b e_j; on the diagonal the two shifts add to coordinate i.
  n <- 4L * nrow(up)
  i <- rep(up[, 1L], each = 4L)
  j <- rep(up[, 2L], each = 4L)
  a <- rep(c(eps, eps, -eps, -eps), nrow(up))
  b <- rep(c(eps, -eps, eps, -eps), nrow(up))
  f <- numeric(n)
  for (first in seq.int(1L, n, by = max_rows)) {
    k <- first:min(first + max_rows - 1L, n)
    rows <- seq_along(k)
    proposals <- matrix(rep(x, each = length(k)), nrow = length(k))
    proposals[cbind(rows, i[k])] <- proposals[cbind(rows, i[k])] + a[k]
    proposals[cbind(rows, j[k])] <- proposals[cbind(rows, j[k])] + b[k]
    v <- -ll_batch(proposals)
    if (length(v) != length(k) || any(!is.finite(v))) return(NULL)
    f[k] <- v
  }
  f <- matrix(f, 4L)
  # fs[[ab]][i, j] = -ll(x + a e_i + b e_j); the lower triangle is the upper
  # one of the sign-swapped matrix.
  low <- lower.tri(diag(p))
  fs <- lapply(1:4, function(ab) {
    m <- matrix(NA_real_, p, p)
    m[up] <- f[ab, ]
    m
  })
  swap <- c(1L, 3L, 2L, 4L)
  for (ab in 1:4) fs[[ab]][low] <- t(fs[[swap[ab]]])[low]
  grad_plus <- (fs[[1L]] - fs[[2L]]) / (2 * eps)    # [i, j]: d_j at x + eps e_i
  grad_minus <- (fs[[3L]] - fs[[4L]]) / (2 * eps)
  hessian <- (grad_plus - grad_minus) / (2 * eps)
  (hessian + t(hessian)) / 2
}

.emc_burn_hessian_step <- function(x, prior_var, precision, ll_batch, step,
                                   allow_shrink = TRUE, blocks = NULL) {
  ll_hess <- .emc_blockwise_ll_hessian(x, ll_batch, blocks,
                                       .emc_stencil_ll_hessian, step)
  if (is.null(ll_hess)) return(NULL)
  .emc_hessian_covariance(precision + ll_hess, precision, prior_var,
                          allow_shrink)
}

# Four-point mixed-difference Hessian of -ll; NULL if any point is nonfinite.
.emc_stencil_ll_hessian <- function(x, ll_batch, step) {
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
  hessian <- matrix(0, p, p)
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
  hessian
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
                                  proposal_prior_scale = 1, blocks = NULL) {
  p <- length(x0)
  root <- tryCatch(chol(prior$var), error = function(e) NULL)
  if (is.null(root)) return(list(reason = "prior covariance"))
  precision <- chol2inv(root)
  best <- new.env(parent = emptyenv())
  best$value <- Inf
  best$par <- x0
  # optim() calls the gradient at the point it just evaluated, so keep the
  # last likelihood rather than recomputing it there.
  last <- new.env(parent = emptyenv())
  loglik <- function(x) {
    if (!identical(x, last$x)) {
      last$ll <- tryCatch(ll_batch(x), error = function(e) NA_real_)
      last$x <- x
    }
    last$ll
  }
  objective <- function(x) {
    ll <- loglik(x)
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
  # A separable likelihood's partial derivatives in block b need only block
  # b's likelihood, evaluated on block b's coordinates.
  groups <- if (is.null(blocks)) list(list(idx = seq_len(p), ll = NULL)) else blocks
  gradient <- function(x) {
    f <- objective(x)
    if (!is.finite(f) || f >= 1e10) stop("nonfinite central objective")
    ll_gradient <- rep(NA_real_, p)
    for (g in groups) {
      xg <- x[g$idx]
      ll_g <- if (is.null(g$ll)) ll_batch else g$ll
      # The centre is needed only for a one-sided difference.
      center_ll <- if (is.null(g$ll)) loglik(x) else NULL
      steps <- grad_step[g$idx]
      grad_g <- rep(NA_real_, length(xg))
      for (attempt in seq_len(8L)) {
        pending <- which(!is.finite(grad_g))
        if (!length(pending)) break
        proposals <- matrix(rep(xg, each = 2L * length(pending)),
                            nrow = 2L * length(pending))
        for (i in seq_along(pending)) {
          k <- pending[i]
          proposals[2L * i - 1L, k] <- xg[k] + steps[k]
          proposals[2L * i, k] <- xg[k] - steps[k]
        }
        ll <- tryCatch(ll_g(proposals), error = function(e)
          rep(NA_real_, nrow(proposals)))
        for (i in seq_along(pending)) {
          k <- pending[i]
          plus <- ll[2L * i - 1L]
          minus <- ll[2L * i]
          if (is.finite(plus) && is.finite(minus)) {
            grad_g[k] <- (plus - minus) / (2 * steps[k])
            next
          }
          if (attempt < 8L) {
            steps[k] <- steps[k] / 4
            next
          }
          if (is.null(center_ll))
            center_ll <- tryCatch(ll_g(xg), error = function(e) NA_real_)
          if (is.finite(center_ll))
            grad_g[k] <- if (is.finite(plus))
              (plus - center_ll) / steps[k] else if (is.finite(minus))
                (center_ll - minus) / steps[k] else NA_real_
        }
      }
      grad_step[g$idx] <<- pmin(grad_step[g$idx], steps)
      ll_gradient[g$idx] <- grad_g
    }
    if (any(!is.finite(ll_gradient))) stop("nonfinite likelihood gradient")
    -ll_gradient + drop(precision %*% (x - prior$mu))
  }
  # BFGS. nlminb, restarted BFGS and Pathfinder were compared on five fixtures
  # and none improved the starts it gives (benchmarks/mode-start/README.md).
  opt <- tryCatch(stats::optim(x0, objective, gradient, method = "BFGS",
                               control = list(maxit = maxit, reltol = 1e-7)),
                  error = function(e) e)
  if (is.list(opt) && !inherits(opt, "error") &&
      all(is.finite(opt$par))) objective(opt$par)
  x <- best$par
  ll <- loglik(x)
  converged <- is.list(opt) && !inherits(opt, "error") &&
    identical(opt$convergence, 0L)
  grad <- tryCatch(gradient(x), error = function(e) rep(Inf, p))
  stationary <- all(is.finite(grad)) &&
    max(abs(grad * sqrt(diag(prior$var)))) < 0.1
  # The mode is found under the actual group prior. The proposal covariance
  # can use a wider one: directions the likelihood barely informs are then not
  # tied to the current group variance, which in a centred hierarchy would
  # otherwise let a shrinking group SD shrink the proposals with it.
  covariance <- if (converged || stationary) {
    tryCatch(.emc_burn_hessian(x, prior$var * proposal_prior_scale,
                               precision / proposal_prior_scale, ll_batch,
                               blocks = blocks),
             error = function(e) NULL)
  } else NULL
  list(alpha = x, ll = ll, covariance = covariance,
       gain = f0 - best$value, converged = converged,
       stationary = stationary,
       reason = if (!is.null(covariance)) "ok" else
         if (!converged && !stationary) "not stationary" else "Hessian")
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
                        proposal_prior_scale = proposal_prior_scale,
                        blocks = .emc_joint_ll_blocks(sampler, s))
}

# One block per model of a joint fit, each with that model's likelihood on its
# own parameters (shared_ll_idx is the model index of every parameter). NULL
# when there is a single likelihood.
.emc_joint_ll_blocks <- function(sampler, s) {
  shared <- attr(sampler$data, "shared_ll_idx")
  if (is.null(shared) || length(unique(shared)) < 2L ||
      length(shared) != sampler$n_pars) return(NULL)
  par_names <- rownames(sampler$samples$alpha)
  lapply(unique(shared), function(k) {
    idx <- which(shared == k)
    list(idx = idx, ll = function(x) .emc_burn_ll_batch(
      x, sampler$data[[s]], sampler$model, par_names[idx], component = k))
  })
}

# A chain start drawn approximately from the subject's conditional posterior
# (given the current group state): importance-resample n_draws draws of the
# Laplace Gaussian N(mode, covariance), weighted posterior / Gaussian with
# Pareto-smoothed weights. The Pareto k of those weights says how well the
# Gaussian matches the posterior (above 0.7 it does not: typically a mode
# against the support boundary or a skewed ridge). NULL if fewer than two
# draws have a finite likelihood. k estimates the tail shape of the weights,
# a property of the Gaussian-posterior pair rather than of n_draws: at 64
# draws it missed 4 of 15 poor subjects on a p = 71 fixture, from 256 up to
# 4096 the counts were stable (2026-09-29). No n rescues k > 1 (weights of
# infinite variance), so n is fixed rather than scaled with p.
.emc_laplace_is_start <- function(mode, covariance, prior, ll_batch,
                                  n_draws = 256L) {
  p <- length(mode)
  z <- matrix(stats::rnorm(n_draws * p), p, n_draws)
  draws <- t(mode + t(chol(covariance)) %*% z)
  ll <- ll_batch(draws)
  # Normalising constants of both densities cancel in the weights.
  u <- forwardsolve(t(chol(prior$var)), t(draws) - prior$mu)
  log_ratio <- ll - 0.5 * colSums(u^2) + 0.5 * colSums(z^2)
  ok <- which(is.finite(log_ratio))
  if (length(ok) < 2L) return(NULL)
  if (diff(range(log_ratio[ok])) < 1e-3) {
    # Weights equal to 0.1%: the tail fit is meaningless (k is Inf for
    # constant ratios and arbitrary for rounding noise).
    w <- rep(1 / length(ok), length(ok))
    k <- 0
  } else {
    smoothed <- suppressWarnings(loo::psis(log_ratio[ok], r_eff = 1))
    w <- as.numeric(stats::weights(smoothed, log = FALSE, normalize = TRUE))
    k <- as.numeric(smoothed$diagnostics$pareto_k)
  }
  pick <- ok[sample.int(length(ok), 1L, prob = w)]
  # Draws outside the likelihood support have weight zero; k only sees the
  # others, so the share outside is reported alongside it.
  list(start = draws[pick, ], start_ll = ll[pick], pareto_k = k,
       outside = 1 - length(ok) / n_draws)
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
        ll_batch <- function(x) .emc_burn_ll_batch(
          x, sampler$data[[s]], sampler$model, rownames(sampler$samples$alpha))
        # Laplace importance resampling: always run for its diagnostic (how well
        # the Gaussian matches the conditional posterior); its draw is the
        # start only under emc2.mode_start_is, until burn-in comparisons show
        # it beats the spread start below.
        use_is <- isTRUE(getOption("emc2.mode_start_is", FALSE))
        seed <- get(".Random.seed", envir = globalenv())
        prior <- .emc_burn_subject_prior(sampler, s)
        is_start <- tryCatch(.emc_laplace_is_start(
          result$alpha, result$covariance, prior, ll_batch),
          error = function(e) NULL)
        if (!is.null(is_start)) {
          result$pareto_k <- is_start$pareto_k
          result$outside <- is_start$outside
        }
        # As a diagnostic only, it leaves the chain's random stream untouched.
        if (!use_is) assign(".Random.seed", seed, envir = globalenv())
        # A resampled draw is only as good as its weights: none above k = 0.7.
        start <- if (use_is && !is.null(is_start) && is_start$pareto_k <= 0.7)
          is_start
        if (is.null(start)) {
          # Spread the chains across the local posterior neighbourhood. Reject
          # proposals outside the likelihood support.
          root <- chol(3 * result$covariance)
          for (attempt in seq_len(8L)) {
            draw <- as.numeric(result$alpha +
                                 drop(stats::rnorm(length(result$alpha)) %*% root))
            ll <- tryCatch(ll_batch(draw), error = function(e) NA_real_)
            if (length(ll) == 1L && is.finite(ll)) {
              start <- list(start = draw, start_ll = ll)
              break
            }
          }
        }
        # The mode is the backup start.
        if (is.null(start)) start <- list(start = result$alpha, start_ll = result$ll)
        result$start <- start$start
        result$start_ll <- start$start_ll
      }
      list(result = result, seed = get(".Random.seed", envir = globalenv()))
    }, mc.cores = min(nrow(jobs), n_cores), mc.set.seed = FALSE))
  for (j in seq_along(emc)) {
    emc[[j]]$burn_mode_stats <- list(
      gain = rep(NA_real_, emc[[j]]$n_subjects),
      converged = rep(FALSE, emc[[j]]$n_subjects),
      hessian = rep(FALSE, emc[[j]]$n_subjects),
      reason = rep("optimizer error", emc[[j]]$n_subjects),
      pareto_k = rep(NA_real_, emc[[j]]$n_subjects),
      outside = rep(NA_real_, emc[[j]]$n_subjects))
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
    if (!is.null(result$pareto_k)) {
      emc[[j]]$burn_mode_stats$pareto_k[s] <- result$pareto_k
      emc[[j]]$burn_mode_stats$outside[s] <- result$outside
    }
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
  k <- unlist(lapply(stats, `[[`, "pareto_k"))
  outside <- unlist(lapply(stats, `[[`, "outside"))[is.finite(k)]
  k <- k[is.finite(k)]
  if (verbose || n_hessian < nrow(jobs)) {
    message("Burn mode starts: ", n_hessian, "/", nrow(jobs),
            " Hessians; ", n_converged, "/", nrow(jobs), " BFGS converged",
            if (n_hessian < nrow(jobs))
              paste0("; failures: ", paste(names(failures), failures,
                                           collapse = ", ")) else "",
            if (length(k)) paste0("; Laplace poor (Pareto k > 0.7 or over half",
                                  " its mass outside the support) for ",
                                  sum(k > 0.7 | outside > 0.5), "/",
                                  length(k)) else "")
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

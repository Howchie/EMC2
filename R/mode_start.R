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
    # The gradient-difference Hessian is the four-point stencil (diagonal at
    # twice the step) at a fixed step of 1e-3, the scale optimHess() uses;
    # what follows is a ladder of steps scaled to the prior SD.
    ll_hess <- tryCatch(.emc_blockwise_ll_hessian(
      x, ll_batch, blocks,
      function(x, ll_batch, step) .emc_gradient_difference_hessian(x, ll_batch)),
      error = function(e) NULL)
    if (is.matrix(ll_hess) && all(is.finite(ll_hess))) {
      ll_hess <- (ll_hess + t(ll_hess)) / 2
      out <- .emc_hessian_covariance(precision + ll_hess,
                                     precision, prior_var, allow_shrink = FALSE)
      if (is.matrix(out)) {
        attr(out, "ll_hessian") <- ll_hess
        return(out)
      }
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
# steps eps). Algebraically this is the four-point mixed stencil of
# .emc_stencil_ll_hessian() at step eps, with the diagonal at step 2 eps. optimHess() evaluates 4p^2 points one likelihood call at a time.
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
  out <- .emc_hessian_covariance(precision + ll_hess, precision, prior_var,
                                 allow_shrink)
  if (is.matrix(out)) attr(out, "ll_hessian") <- ll_hess
  out
}

# Four-point mixed-difference Hessian of -ll; NULL if any point is nonfinite.
.emc_stencil_ll_hessian <- function(x, ll_batch, step) {
  p <- length(x)
  # A step scaled to the prior can cross a sharp likelihood ridge and give
  # spurious negative cross-curvature at an optimum.
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
                                  proposal_prior_scale = 1, blocks = NULL,
                                  start_covariance = FALSE) {
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
  # The mode is found under the actual group prior. The proposal Gaussian uses
  # a wider one: directions the likelihood barely informs are then not tied to
  # the current group variance, which in a centred hierarchy would otherwise
  # let a shrinking group SD shrink the proposals with it. Its covariance is
  # (P / s + H)^-1 and its mean is that density's mode, one Newton step from
  # the mode x under P: the gradient of ll - (x - m)' (P / s) (x - m) / 2 at x
  # is (1 - 1 / s) P (x - m), exact when ll is quadratic. Mean and covariance
  # then describe the same density. The same likelihood Hessian gives the
  # conditional covariance under P itself (start_covariance), for chain starts
  # and the Laplace diagnostics.
  s_scale <- proposal_prior_scale
  covariance <- if (converged || stationary) {
    tryCatch(.emc_burn_hessian(x, prior$var * s_scale,
                               precision / s_scale, ll_batch,
                               blocks = blocks),
             error = function(e) NULL)
  } else NULL
  proposal_mode <- x
  start_cov <- NULL
  lik <- NULL
  if (is.matrix(covariance)) {
    ll_hess <- attr(covariance, "ll_hessian")
    attr(covariance, "ll_hessian") <- NULL
    if (is.matrix(ll_hess))
      lik <- tryCatch(.emc_lik_canonical(ll_hess, x, precision, prior$mu),
                      error = function(e) NULL)
    if (s_scale != 1) {
      shift <- drop(covariance %*% ((1 - 1 / s_scale) *
                                      drop(precision %*% (x - prior$mu))))
      candidate <- x + shift
      cand_ll <- if (all(is.finite(candidate)))
        tryCatch(ll_batch(candidate), error = function(e) NA_real_) else NA_real_
      if (length(cand_ll) == 1L && is.finite(cand_ll)) proposal_mode <- candidate
      if (isTRUE(start_covariance) && is.matrix(ll_hess))
        start_cov <- tryCatch(.emc_hessian_covariance(
          precision + ll_hess, precision, prior$var, allow_shrink = TRUE),
          error = function(e) NULL)
    } else if (isTRUE(start_covariance)) start_cov <- covariance
  }
  list(alpha = x, ll = ll, covariance = covariance,
       proposal_mode = proposal_mode,
       start_covariance = if (is.matrix(start_cov)) start_cov else NULL,
       lik = lik, gain = f0 - best$value, converged = converged,
       stationary = stationary,
       reason = if (!is.null(covariance)) "ok" else
         if (!converged && !stationary) "not stationary" else "Hessian")
}

# Gaussian approximation of the subject likelihood alone, in canonical form
# exp(lin' a - a' prec a / 2): prec is the likelihood Hessian at the
# conditional mode x (projected to positive semi-definite), and lin follows
# from the mode condition, where the likelihood gradient equals
# precision (x - mu). It does not involve the group prior, so new_particle()
# combines it with whatever group state is current.
.emc_lik_canonical <- function(ll_hess, x, precision, mu) {
  e <- eigen((ll_hess + t(ll_hess)) / 2, symmetric = TRUE)
  prec <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
  lin <- drop(prec %*% x + precision %*% (x - mu))
  if (!all(is.finite(prec)) || !all(is.finite(lin))) return(NULL)
  list(prec = prec, lin = lin)
}

# The covariance of the conditional posterior under the actual group prior
# when it was computed, else the proposal covariance.
.emc_start_cov <- function(r) {
  if (is.matrix(r$start_covariance)) r$start_covariance else r$covariance
}

.emc_burn_subject_mode <- function(sampler, s, proposal_prior_scale = 1,
                                   x0 = NULL, start_covariance = FALSE) {
  store <- sampler$samples
  i <- store$idx
  p <- sampler$n_pars
  if (is.null(x0)) x0 <- as.numeric(store$alpha[, s, i])
  prior <- tryCatch(.emc_burn_subject_prior(sampler, s), error = function(e) NULL)
  if (is.null(prior) || length(prior$mu) != p ||
      any(!is.finite(prior$mu)) || any(!is.finite(prior$var))) return(NULL)
  par_names <- rownames(store$alpha)
  ll_batch <- function(x) .emc_burn_ll_batch(
    x, sampler$data[[s]], sampler$model, par_names)
  .emc_optimize_subject(x0, prior, ll_batch,
                        proposal_prior_scale = proposal_prior_scale,
                        blocks = .emc_joint_ll_blocks(sampler, s),
                        start_covariance = start_covariance)
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
# against the support boundary or a skewed ridge). Always returns a status:
# "ok", or "support" when fewer than two draws have a finite weight (no start,
# k NA, outside the share of draws without one). The population k is a
# property of the Gaussian-posterior pair, but its estimate from n_draws draws
# is noisy, and the Gaussian assessed is the conditional Laplace fit, not the
# widened proposal the burn stage executes: k is a proxy for that proposal.
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
  if (length(ok) < 2L)
    return(list(start = NULL, start_ll = NULL, pareto_k = NA_real_,
                outside = 1 - length(ok) / n_draws, status = "support"))
  if (diff(range(log_ratio[ok])) < 1e-3) {
    # Equal weights: the tail fit is meaningless.
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
       outside = 1 - length(ok) / n_draws, status = "ok")
}

# The Laplace diagnostic of one fit attempt, total: a failed check is recorded
# as such rather than leaving an earlier result in place.
.emc_laplace_check <- function(is_result) {
  if (is.list(is_result) && !is.null(is_result$status))
    return(is_result[c("pareto_k", "outside", "status")])
  list(pareto_k = NA_real_, outside = NA_real_, status = "error")
}

# Subjects whose Laplace approximation should not be trusted: Pareto k above
# 0.7 (including Inf), over half the Gaussian mass outside the likelihood
# support, or no usable diagnostic at all (k NA).
.emc_laplace_poor <- function(stats) {
  k <- stats$pareto_k
  if (is.null(k)) return(logical(0))
  outside <- stats$outside
  if (length(outside) != length(k)) outside <- rep(NA_real_, length(k))
  is.na(k) | k > 0.7 | (!is.na(outside) & outside > 0.5)
}

# Optimisation probes the likelihood far from the posterior, where it can warn.
# In a forked child those conditions reach any handler inherited from the
# parent (testthat's parallel reporter writes them to its pipe from the child
# and corrupts the stream), and they carry no information here.
.emc_quiet_fork <- function(expr) suppressWarnings(suppressMessages(expr))

# Laplace log mass of a mode under a conditional posterior, up to a constant:
# log posterior at the mode plus half the log determinant of its covariance.
.emc_laplace_log_mass <- function(ll, x, prior, covariance) {
  u <- forwardsolve(t(chol(prior$var)), x - prior$mu)
  ll - 0.5 * sum(u^2) + sum(log(diag(chol(covariance))))
}

# A chain's mode moves to another chain's mode of the same subject when that
# has more Laplace mass, averaged over all chains' group states, by this many
# nats. Averaging matters: a chain stuck in a local mode has widened its own
# group variance around it, so under its own group state the trap still wins.
.emc_mode_share_gap <- 20

# Proposal covariances (and their means) use the group prior widened by this
# factor, at the burn mode start and at every refresh alike.
.emc_laplace_prior_scale <- 4

.emc_usable_mode <- function(r) is.list(r) && !inherits(r, "try-error") &&
  is.matrix(r$covariance)

# Laplace log mass of mode r averaged over a subject's group states (one prior
# per chain); -Inf for an unusable mode.
.emc_shared_mode_mass <- function(r, priors) {
  if (!.emc_usable_mode(r) || !is.finite(r$ll)) return(-Inf)
  m <- vapply(priors, function(prior) if (is.null(prior)) NA_real_ else
    tryCatch(.emc_laplace_log_mass(r$ll, r$alpha, prior, .emc_start_cov(r)),
             error = function(e) -Inf), numeric(1))
  if (all(is.na(m))) -Inf else mean(m, na.rm = TRUE)
}

# For each chain x subject job, the chain whose mode of that subject beats its
# own by more than .emc_mode_share_gap; NA to keep its own.
.emc_mode_share_sources <- function(jobs, results, priors) {
  out <- rep(NA_integer_, nrow(jobs))
  for (k in seq_len(nrow(jobs))) {
    if (is.null(priors[[k]])) next
    same <- which(jobs$subject == jobs$subject[k])
    subject_priors <- priors[same]
    own <- .emc_shared_mode_mass(results[[k]], subject_priors)
    others <- setdiff(same, k)
    cand <- vapply(others, function(kk)
      .emc_shared_mode_mass(results[[kk]], subject_priors), numeric(1))
    if (length(cand) && max(cand) > own + .emc_mode_share_gap)
      out[k] <- jobs$chain[others[which.max(cand)]]
  }
  out
}

.emc_burn_mode_start <- function(emc, n_cores) {
  if (!isTRUE(getOption("emc2.burn_mode_start", TRUE))) return(emc)
  if (any(emc[[1L]]$nuisance)) return(emc)
  # Subject modes are conditional on each chain's current pooled values.
  for (j in seq_along(emc)) emc[[j]] <- .emc_pooled_sync(emc[[j]])
  jobs <- expand.grid(chain = seq_along(emc),
                      subject = seq_len(emc[[1L]]$n_subjects))
  cores <- min(nrow(jobs), n_cores)
  # 1. Every chain's own conditional mode for every subject (no RNG used).
  results <- .emc_with_preserved_rng(auto_mclapply(
    seq_len(nrow(jobs)), function(k) {
      .emc_quiet_fork(.emc_burn_subject_mode(
        emc[[jobs$chain[k]]], jobs$subject[k],
        proposal_prior_scale = .emc_laplace_prior_scale,
        start_covariance = TRUE))
    }, mc.cores = cores, mc.set.seed = FALSE))
  # 2. Compare each subject's modes across chains.
  priors <- lapply(seq_len(nrow(jobs)), function(k) tryCatch(
    .emc_burn_subject_prior(emc[[jobs$chain[k]]], jobs$subject[k]),
    error = function(e) NULL))
  mass <- function(k, r) .emc_shared_mode_mass(
    r, priors[jobs$subject == jobs$subject[k]])
  switch_from <- if (.emc_independent_chains()) rep(NA_integer_, nrow(jobs)) else
    .emc_mode_share_sources(jobs, results, priors)
  # 3. Re-optimise from the better mode under the chain's own group state; keep
  # it only if it still wins by the same margin.
  todo <- which(!is.na(switch_from))
  if (length(todo)) {
    redo <- .emc_with_preserved_rng(auto_mclapply(todo, function(k) {
      src <- results[[which(jobs$chain == switch_from[k] &
                            jobs$subject == jobs$subject[k])]]
      .emc_quiet_fork(.emc_burn_subject_mode(
        emc[[jobs$chain[k]]], jobs$subject[k], x0 = src$alpha,
        proposal_prior_scale = .emc_laplace_prior_scale,
        start_covariance = TRUE))
    }, mc.cores = min(length(todo), n_cores), mc.set.seed = FALSE))
    for (m in seq_along(todo)) {
      k <- todo[m]
      new <- redo[[m]]
      if (!.emc_usable_mode(new)) next
      if (mass(k, new) > mass(k, results[[k]]) + .emc_mode_share_gap) {
        new$shared_from <- switch_from[k]
        results[[k]] <- new
      }
    }
  }
  # 4. Start draws, each on its chain's subject RNG stream.
  results <- .emc_with_preserved_rng(auto_mclapply(
    seq_len(nrow(jobs)), function(k) {
      j <- jobs$chain[k]
      s <- jobs$subject[k]
      sampler <- emc[[j]]
      result <- results[[k]]
      assign(".Random.seed", sampler$rng$subjects[[s]], envir = globalenv())
      if (.emc_usable_mode(result)) {
        ll_batch <- function(x) .emc_burn_ll_batch(
          x, sampler$data[[s]], sampler$model, rownames(sampler$samples$alpha))
        # Diagnostic probes consume the subject stream before the start draw.
        prior <- .emc_burn_subject_prior(sampler, s)
        is_start <- tryCatch(.emc_laplace_is_start(
          result$alpha, .emc_start_cov(result), prior, ll_batch),
          error = function(e) NULL)
        result[c("pareto_k", "outside", "laplace_status")] <-
          unname(.emc_laplace_check(is_start))
        # Spread the chains across the local posterior neighbourhood. Reject
        # proposals outside the likelihood support.
        start <- NULL
        root <- chol(3 * .emc_start_cov(result))
        for (attempt in seq_len(8L)) {
          draw <- as.numeric(result$alpha +
                               drop(stats::rnorm(length(result$alpha)) %*% root))
          ll <- tryCatch(ll_batch(draw), error = function(e) NA_real_)
          if (length(ll) == 1L && is.finite(ll)) {
            start <- list(start = draw, start_ll = ll)
            break
          }
        }
        # The mode is the backup start.
        if (is.null(start)) start <- list(start = result$alpha, start_ll = result$ll)
        result$start <- start$start
        result$start_ll <- start$start_ll
      }
      list(result = result, seed = get(".Random.seed", envir = globalenv()))
    }, mc.cores = cores, mc.set.seed = FALSE))
  for (j in seq_along(emc)) {
    emc[[j]]$burn_mode_stats <- list(
      gain = rep(NA_real_, emc[[j]]$n_subjects),
      converged = rep(FALSE, emc[[j]]$n_subjects),
      hessian = rep(FALSE, emc[[j]]$n_subjects),
      reason = rep("optimizer error", emc[[j]]$n_subjects),
      pareto_k = rep(NA_real_, emc[[j]]$n_subjects),
      outside = rep(NA_real_, emc[[j]]$n_subjects),
      laplace_status = rep("none", emc[[j]]$n_subjects),
      shared_from = rep(NA_integer_, emc[[j]]$n_subjects))
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
    if (!is.null(result$shared_from))
      emc[[j]]$burn_mode_stats$shared_from[s] <- result$shared_from
    if (!is.null(result$laplace_status)) {
      emc[[j]]$burn_mode_stats$pareto_k[s] <- result$pareto_k
      emc[[j]]$burn_mode_stats$outside[s] <- result$outside
      emc[[j]]$burn_mode_stats$laplace_status[s] <- result$laplace_status
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
    emc[[j]]$burn_mode_mu[s] <- list(result$proposal_mode)
    emc[[j]] <- .emc_store_lik(emc[[j]], s, result$lik)
  }
  emc
}

.emc_store_lik <- function(chain, s, lik) {
  if (is.null(lik)) return(chain)
  if (is.null(chain$lik_laplace))
    chain$lik_laplace <- vector("list", chain$n_subjects)
  chain$lik_laplace[s] <- list(lik)
  chain
}

# Check the Gaussian proposal from each likelihood fit and the current subject
# prior. The check belongs to that fit, including any change to its location.
.emc_check_likelihood_fits <- function(emc, n_cores) {
  has_fits <- vapply(emc, function(chain)
    any(vapply(chain$lik_laplace, is.list, logical(1))), logical(1))
  if (!any(has_fits)) return(emc)
  emc <- .emc_prepare_chain_rng(emc)
  jobs <- expand.grid(chain = seq_along(emc), subject = seq_len(emc[[1]]$n_subjects))
  present <- vapply(seq_len(nrow(jobs)), function(k) {
    lik <- emc[[jobs$chain[k]]]$lik_laplace
    length(lik) >= jobs$subject[k] && is.list(lik[[jobs$subject[k]]])
  }, logical(1))
  jobs <- jobs[present, , drop = FALSE]
  if (!nrow(jobs)) return(emc)
  results <- .emc_with_preserved_rng(auto_mclapply(seq_len(nrow(jobs)), function(k) {
    chain <- emc[[jobs$chain[k]]]
    s <- jobs$subject[k]
    fit <- chain$lik_laplace[[s]]
    .emc_quiet_fork(.emc_with_preserved_rng({
      check <- tryCatch({
        assign(".Random.seed", chain$rng$subjects[[s]], envir = globalenv())
        prior <- .emc_burn_subject_prior(chain, s)
        gaussian <- .emc_laplace_conditional(structure(numeric(chain$n_pars), lik = fit), prior)
        if (is.null(gaussian)) stop("unusable likelihood approximation")
        ll <- function(x) .emc_burn_ll_batch(x, chain$data[[s]], chain$model, chain$par_names)
        .emc_laplace_is_start(gaussian$mu, gaussian$var, prior, ll)
      }, error = function(e) NULL)
      list(check = .emc_laplace_check(check),
           seed = get(".Random.seed", envir = globalenv()))
    }))
  }, mc.cores = min(nrow(jobs), n_cores), mc.set.seed = FALSE))
  for (k in seq_len(nrow(jobs))) {
    result <- results[[k]]
    check <- if (is.list(result)) result$check else NULL
    if (!is.list(check)) check <- .emc_laplace_check(NULL)
    emc[[jobs$chain[k]]]$lik_laplace[[jobs$subject[k]]]$check <- check
    if (is.list(result) && !is.null(result$seed))
      emc[[jobs$chain[k]]]$rng$subjects[[jobs$subject[k]]] <- result$seed
  }
  emc
}

# Laplace proposal metric for the discarded stages: at each proposal refresh,
# find every chain's conditional mode for each subject under that chain's
# current group state and store the Laplace covariance (and mode). The chain
# state is not moved. create_chain_proposals() consumes these once.
.emc_laplace_refresh <- function(emc, n_cores) {
  if (!isTRUE(getOption("emc2.laplace_refresh", TRUE))) return(emc)
  if (any(emc[[1L]]$nuisance)) return(emc)
  for (j in seq_along(emc)) emc[[j]] <- .emc_pooled_sync(emc[[j]])
  emc <- .emc_prepare_chain_rng(emc)
  jobs <- expand.grid(chain = seq_along(emc),
                      subject = seq_len(emc[[1L]]$n_subjects))
  cost <- .emc_subject_cost(emc[[1L]]$data)
  # Each refresh checks the Laplace fit under the current group state.
  # Probe draws advance the subject stream before subsequent transitions.
  results <- .emc_with_preserved_rng(.emc_rebuild_lapply(
    seq_len(nrow(jobs)), function(k) {
      .emc_quiet_fork({
        sampler <- emc[[jobs$chain[k]]]
        s <- jobs$subject[k]
        assign(".Random.seed", sampler$rng$subjects[[s]], envir = globalenv())
        r <- tryCatch({
          r <- .emc_burn_subject_mode(sampler, s,
                                      proposal_prior_scale = .emc_laplace_prior_scale,
                                      start_covariance = TRUE)
          if (.emc_usable_mode(r)) {
            ll_batch <- function(x) .emc_burn_ll_batch(
              x, sampler$data[[s]], sampler$model, rownames(sampler$samples$alpha))
            is_check <- tryCatch(.emc_laplace_is_start(
              r$alpha, .emc_start_cov(r), .emc_burn_subject_prior(sampler, s),
              ll_batch), error = function(e) NULL)
            r[c("pareto_k", "outside", "laplace_status")] <-
              unname(.emc_laplace_check(is_check))
          }
          r
        }, error = function(e) NULL)
        list(result = r, seed = get(".Random.seed", envir = globalenv()))
      })
    }, n_cores = n_cores, order = order(cost[jobs$subject], decreasing = TRUE)))
  for (k in seq_len(nrow(jobs))) {
    job_result <- results[[k]]
    result <- if (is.list(job_result)) job_result$result else NULL
    j <- jobs$chain[k]
    s <- jobs$subject[k]
    if (is.list(job_result) && !is.null(job_result$seed))
      emc[[j]]$rng$subjects[[s]] <- job_result$seed
    n_sub <- emc[[j]]$n_subjects
    usable <- is.list(result) && is.matrix(result$covariance)
    if (is.null(emc[[j]]$burn_mode_var))
      emc[[j]]$burn_mode_var <- vector("list", n_sub)
    if (is.null(emc[[j]]$burn_mode_mu))
      emc[[j]]$burn_mode_mu <- vector("list", n_sub)
    ms <- emc[[j]]$burn_mode_stats
    if (length(ms$pareto_k) != n_sub) ms$pareto_k <- rep(NA_real_, n_sub)
    if (length(ms$outside) != n_sub) ms$outside <- rep(NA_real_, n_sub)
    if (length(ms$laplace_status) != n_sub)
      ms$laplace_status <- rep("none", n_sub)
    # A failed refresh drops the earlier curvature and its diagnostic: both
    # belonged to an older group state.
    emc[[j]]$burn_mode_var[s] <- list(if (usable) result$covariance)
    emc[[j]]$burn_mode_mu[s] <- list(if (usable) result$proposal_mode)
    # The likelihood approximation does not depend on the group state, so a
    # failed refresh keeps the earlier one.
    if (usable) emc[[j]] <- .emc_store_lik(emc[[j]], s, result$lik)
    chk <- if (usable && !is.null(result$laplace_status))
      result[c("pareto_k", "outside", "laplace_status")] else
        list(NA_real_, NA_real_, if (usable) "error" else "mode failed")
    ms$pareto_k[s] <- chk[[1]]
    ms$outside[s] <- chk[[2]]
    ms$laplace_status[s] <- chk[[3]]
    emc[[j]]$burn_mode_stats <- ms
  }
  emc
}

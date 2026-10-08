# Pooled parameters: design(pooled = ...).
#
# A pooled coefficient is shared by every subject and sampled once per
# iteration. The hierarchy never sees it: sampled_pars() leaves it out by
# default, so alpha, the group-level samplers, the particle step and the
# ensemble move work on the remaining parameters unchanged. The current value
# rides on every subject's data as attr(dadm, "pooled_values") (sampled scale)
# and .emc_pooled_columns() appends it to proposal matrices at the likelihood
# boundary, where the C++ code treats it like any other proposal column
# (pre-transform included).
#
# Sweep: subjects | pooled (particle step, ensemble), then the pooled move.
# The move is a fixed shear: a candidate phi' = phi + delta carries every
# subject and the group mean along the same slope c (alpha_s + c delta,
# mu + c delta), so every alpha_s - mu, and with it the subject prior, is
# unchanged. c and the step covariance come from a finite-difference
# curvature grid at the current state. They are refreshed once per block
# outside the sample stage, when the other proposals are rebuilt, and frozen
# in it. The map has unit Jacobian, so the target
# along it is prior(phi) prior(mu) prod_s L_s(alpha_s, phi). phi is drawn by
# conditional importance sampling with an auxiliary centre: z ~ N(phi, S),
# K candidates ~ N(z, S), select one of {phi, candidates} by that target.
# Given z the candidates do not depend on phi, so the step is invariant. The
# candidates are scored in one worker round trip, and the chosen one's
# per-subject likelihoods become the stored subject log-likelihoods.

.emc_pooled_options <- function() {
  k <- suppressWarnings(as.integer(getOption("emc2.pooled_candidates", 8L)))
  if (length(k) != 1L || is.na(k) || k < 1L) k <- 8L
  list(candidates = k, target = 0.5)
}

.emc_pooled_names <- function(sampler) {
  if (is.null(sampler$pooled)) character(0) else names(sampler$pooled$value)
}

# Proposal matrix plus the current pooled values for any pooled coefficient the
# caller did not supply. Supplied columns win, so post-hoc callers that pass
# per-draw pooled values are unaffected. `varying` (the C++ invariance mask)
# grows with the appended columns, which never vary within a call.
.emc_pooled_columns <- function(proposals, dadm, varying = NULL) {
  pooled <- attr(dadm, "pooled")
  if (is.null(pooled)) return(list(proposals = proposals, varying = varying))
  if (!is.matrix(proposals))
    proposals <- matrix(proposals, nrow = 1L,
                        dimnames = list(NULL, names(proposals)))
  missing <- setdiff(pooled, colnames(proposals))
  if (!length(missing)) return(list(proposals = proposals, varying = varying))
  vals <- attr(dadm, "pooled_values")
  if (is.null(vals) || !all(missing %in% names(vals)) || anyNA(vals[missing]))
    stop("No value for pooled parameter(s) ", paste(missing, collapse = ", "),
         ": supply them as columns of the parameter matrix")
  add <- matrix(rep(unname(vals[missing]), each = nrow(proposals)),
                nrow = nrow(proposals), dimnames = list(NULL, missing))
  if (!is.null(varying)) varying <- c(varying, rep(FALSE, length(missing)))
  list(proposals = cbind(proposals, add), varying = varying)
}

.emc_pooled_set_one <- function(dadm, values) {
  if (!is.null(values) && is.data.frame(dadm)) attr(dadm, "pooled_values") <- values
  dadm
}

.emc_pooled_set_data <- function(data, values) {
  if (is.null(values)) return(data)
  for (i in seq_along(data)) data[[i]] <- .emc_pooled_set_one(data[[i]], values)
  data
}

# Pooled values of the sampler's current iteration on every subject's data.
.emc_pooled_sync <- function(sampler) {
  if (is.null(sampler$pooled)) return(sampler)
  idx <- sampler$samples$idx
  if (!is.null(idx) && !is.null(sampler$samples$pooled) && idx >= 1L) {
    cur <- sampler$samples$pooled[, idx]
    if (!anyNA(cur)) sampler$pooled$value[] <- cur
  }
  sampler$data <- .emc_pooled_set_data(sampler$data, sampler$pooled$value)
  sampler
}

.emc_pooled_unsupported <- function(sampler) {
  if (is.null(sampler$pooled)) return(invisible(NULL))
  if (any(sampler$nuisance))
    stop("pooled parameters cannot be combined with nuisance samplers")
  if (!is.null(sampler$marginalise))
    stop("pooled parameters cannot be combined with marginalise")
  if (!is.null(sampler$group_designs))
    stop("pooled parameters cannot be combined with group_design")
  if (identical(sampler$type, "single"))
    stop("pooled parameters need a hierarchical type")
  invisible(NULL)
}

# Sampler state, added by pmwgs(). The proposal covariance starts at the prior
# variance until the first curvature refresh; the scale is tuned outside the
# sample stage only.
.emc_pooled_state <- function(prior) {
  nm <- names(prior$mean)
  list(value = stats::setNames(as.numeric(prior$mean), nm),
       prior_mean = prior$mean, prior_var = prior$var,
       cov = diag(as.numeric(prior$var), length(nm)),
       log_eps = log(0.1), n_adapt = 0L, shear = NULL,
       moves = c(preburn = 0, burn = 0, adapt = 0, sample = 0),
       tries = c(preburn = 0, burn = 0, adapt = 0, sample = 0))
}

.emc_pooled_log_prior <- function(x, st) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  colSums(stats::dnorm(t(x), st$prior_mean, sqrt(st$prior_var), log = TRUE))
}

# Log-likelihood of each row of a complete parameter matrix; the C++ mapping
# skips columns that are equal in every row. Every row is a candidate (the
# current state's likelihood is carried separately), so a candidate the model
# cannot evaluate gets zero density, as on the particle and ensemble paths. The
# failure policy still applies: a failed batch is classified, aborted under
# "strict" for non-numerical classes, and otherwise retried row by row so one
# bad candidate does not reject the rest. NA/NaN/+Inf values and a result of
# the wrong length are failures too, not silent -Inf. This may run in a forked
# worker, so the failures travel back as an attribute and are counted by
# .emc_pooled_record() in the master.
.emc_pooled_rows_ll <- function(rows, dadm, model) {
  K <- nrow(rows)
  failures <- list()
  fail <- function(cond) {
    cls <- if (inherits(cond, "emc_failure") && !is.null(cond$emc_class))
      cond$emc_class else .emc_classify_failure(cond)
    if (identical(.emc_failure_action(cls), "abort"))
      .emc_failure_abort(cls, "pooled", cond)
    failures[[length(failures) + 1L]] <<- cond
  }
  bad_value <- function(msg, cls)
    errorCondition(msg, class = c("emc_pooled_failure", "emc_failure"),
                   emc_class = cls, emc_source = "pooled")
  eval_rows <- function(r) {
    varying <- if (nrow(r) > 1L) apply(r, 2L, function(x) any(x != x[1L])) else NULL
    v <- as.numeric(calc_ll_manager(r, dadm = dadm, model = model,
                                    varying = varying))
    if (length(v) != nrow(r))
      stop(bad_value(sprintf("pooled likelihood returned %d values for %d rows",
                             length(v), nrow(r)), "programming"))
    v
  }
  ll <- tryCatch(eval_rows(rows), error = identity)
  if (inherits(ll, "error")) {
    fail(ll)
    ll <- if (K == 1L) -Inf else vapply(seq_len(K), function(r) {
      v <- tryCatch(eval_rows(rows[r, , drop = FALSE]), error = identity)
      if (inherits(v, "error")) { fail(v); -Inf } else v
    }, numeric(1))
  }
  odd <- is.na(ll) | ll == Inf
  if (any(odd)) {
    fail(bad_value(sprintf(
      "non-finite pooled likelihood (%d NA/NaN, %d +Inf) at %d candidate(s)",
      sum(is.na(ll)), sum(ll == Inf, na.rm = TRUE), sum(odd)), "numerical"))
    ll[odd] <- -Inf
  }
  if (length(failures)) attr(ll, "emc_failures") <- failures
  ll
}

# Count (and announce) the failures carried back by .emc_pooled_rows_ll().
.emc_pooled_record <- function(ll) {
  conds <- attr(ll, "emc_failures")
  for (cond in conds) .emc_reject_record(cond, "pooled")
  attr(ll, "emc_failures") <- NULL
  ll
}

# Start: each subject's start particles also carry pooled draws from the prior
# (as if the coefficients were per-subject). The pooled start is whichever
# subject's selected draw scores best when all subjects share it (start_mu
# overrides), and the subject log-likelihoods are taken at that value.
.emc_pooled_init <- function(pmwgs, proposals, pooled_draws, start_mu,
                             r_cores = 1) {
  st <- pmwgs$pooled
  nm <- names(st$value)
  n_pars <- pmwgs$n_pars
  cands <- if (!is.null(pooled_draws)) t(matrix(pooled_draws, nrow = length(nm))) else
    matrix(st$value, nrow = 1L)
  if (!is.null(start_mu) && !is.null(names(start_mu))) {
    given <- intersect(nm, names(start_mu))
    cands[, given] <- rep(start_mu[given], each = nrow(cands))
  }
  cands <- unique(cands)
  colnames(cands) <- nm
  ll <- vapply(seq_len(pmwgs$n_subjects), function(s) {
    a <- matrix(proposals[seq_len(n_pars), s], nrow(cands), n_pars, byrow = TRUE)
    rows <- cbind(a, cands)
    colnames(rows) <- c(pmwgs$par_names, nm)
    .emc_pooled_record(.emc_pooled_rows_ll(rows, pmwgs$data[[s]], pmwgs$model))
  }, numeric(nrow(cands)))
  ll <- matrix(ll, nrow = nrow(cands))
  best <- which.max(rowSums(ll))
  st$value[] <- cands[best, ]
  pmwgs$pooled <- st
  pmwgs$data <- .emc_pooled_set_data(pmwgs$data, st$value)
  proposals[n_pars + 1L, ] <- ll[best, ]
  list(sampler = pmwgs, proposals = proposals)
}

.emc_pooled_adapt <- function(st, stage, moved) {
  if (identical(stage, "sample")) return(st)
  st$n_adapt <- st$n_adapt + 1L
  gamma <- 1 / (1 + st$n_adapt / 10)^0.6
  st$log_eps <- st$log_eps + gamma * (as.numeric(moved) - .emc_pooled_options()$target)
  st
}

# Curvature grid: centre and axial rows estimate the diagonal, one forward
# corner per pair involving an alpha coordinate estimates its cross term, and
# central corners estimate the pooled-pooled cross terms. The alpha-alpha
# block is needed: the shear slope solves against it, and its diagonal alone
# can make the curvature along the move negative.
.emc_pooled_fd_rows <- function(x0, h, d, p) {
  D <- length(x0)
  rows <- list(x0)
  for (i in seq_len(D)) {
    plus <- minus <- x0
    plus[i] <- plus[i] + h[i]
    minus[i] <- minus[i] - h[i]
    rows <- c(rows, list(plus, minus))
  }
  for (i in seq_len(d)) for (j in setdiff(seq_len(D), seq_len(i))) {
    corner <- x0
    corner[i] <- corner[i] + h[i]
    corner[j] <- corner[j] + h[j]
    rows <- c(rows, list(corner))
  }
  if (p > 1L) {
    for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
      ii <- d + i
      jj <- d + j
      pp <- pm <- mp <- mm <- x0
      pp[ii] <- pp[ii] + h[ii]; pp[jj] <- pp[jj] + h[jj]
      pm[ii] <- pm[ii] + h[ii]; pm[jj] <- pm[jj] - h[jj]
      mp[ii] <- mp[ii] - h[ii]; mp[jj] <- mp[jj] + h[jj]
      mm[ii] <- mm[ii] - h[ii]; mm[jj] <- mm[jj] - h[jj]
      rows <- c(rows, list(pp, pm, mp, mm))
    }
  }
  out <- do.call(rbind, rows)
  colnames(out) <- names(x0)
  out
}

# Negative Hessian of one subject's log-likelihood from its grid values `fs`.
.emc_pooled_fd_hessian <- function(fs, hs, d, p) {
  D <- d + p
  Q <- matrix(0, D, D)
  ax <- function(i) 2L * i
  for (i in seq_len(D))
    Q[i, i] <- -(fs[ax(i)] - 2 * fs[1L] + fs[ax(i) + 1L]) / hs[i]^2
  k <- 2L * D + 2L
  for (i in seq_len(d)) for (j in setdiff(seq_len(D), seq_len(i))) {
    Q[i, j] <- Q[j, i] <- -(fs[k] - fs[ax(i)] - fs[ax(j)] + fs[1L]) / (hs[i] * hs[j])
    k <- k + 1L
  }
  if (p > 1L) {
    for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
      ii <- d + i
      jj <- d + j
      Q[ii, jj] <- Q[jj, ii] <-
        -(fs[k] - fs[k + 1L] - fs[k + 2L] + fs[k + 3L]) / (4 * hs[ii] * hs[jj])
      k <- k + 4L
    }
  }
  Q
}

# Shear slope c and step covariance from the summed curvature of the subject
# likelihoods (and the group-mean prior) at the current state: c minimises the
# curvature along the move, and the step covariance inverts the curvature left
# along it (the Schur complement of the pooled block). Away from the mode the
# summed curvature need not be positive definite; its eigenvalues are then
# taken in absolute value, which keeps every direction's scale.
.emc_pooled_curvature <- function(pmwgs, st, alpha, score, h = 1e-3) {
  nm <- names(st$value)
  p <- length(nm)
  d <- nrow(alpha)
  D <- d + p
  x0 <- lapply(seq_len(ncol(alpha)), function(s)
    stats::setNames(c(alpha[, s], st$value), c(pmwgs$par_names, nm)))
  steps <- lapply(x0, function(x) h * pmax(1, abs(x)))
  rows <- Map(.emc_pooled_fd_rows, x0, steps, MoreArgs = list(d = d, p = p))
  f <- score(rows)
  if (any(!is.finite(f))) return(st)
  A <- matrix(0, D, D)
  for (s in seq_len(ncol(f))) A <- A + .emc_pooled_fd_hessian(f[, s], steps[[s]], d, p)
  ia <- seq_len(d)
  ip <- d + seq_len(p)
  A[ia, ia] <- A[ia, ia] + .emc_pooled_mu_prec(pmwgs$prior, d)
  A[ip, ip] <- A[ip, ip] + diag(1 / as.numeric(st$prior_var), p)
  A <- (A + t(A)) / 2
  if (inherits(try(chol(A), silent = TRUE), "try-error")) {
    e <- eigen(A, symmetric = TRUE)
    lam <- pmax(abs(e$values), max(abs(e$values)) * .Machine$double.eps)
    A <- e$vectors %*% (lam * t(e$vectors))
  }
  shear <- if (d > 0L) -solve(A[ia, ia, drop = FALSE], A[ia, ip, drop = FALSE]) else
    matrix(0, 0L, p)
  curv <- A[ip, ip, drop = FALSE] + crossprod(A[ia, ip, drop = FALSE], shear)
  curv <- (curv + t(curv)) / 2
  st$shear <- shear
  st$cov <- solve(curv)
  if (st$n_adapt == 0L) st$log_eps <- 0
  st
}

.emc_pooled_mu_prec <- function(prior, d) {
  P <- prior$theta_mu_invar
  if (is.null(P)) {
    V <- prior$theta_mu_var
    P <- if (is.matrix(V)) tryCatch(solve(V), error = function(e) MASS::ginv(V)) else
      diag(1 / as.numeric(V), d)
  } else if (!is.matrix(P)) {
    P <- diag(as.numeric(P), d)
  }
  P[seq_len(d), seq_len(d), drop = FALSE]
}

.emc_pooled_log_prior_mu <- function(mu, prior) {
  d <- length(mu)
  r <- mu - as.numeric(prior$theta_mu_mean)[seq_len(d)]
  -0.5 * sum(r * (.emc_pooled_mu_prec(prior, d) %*% r))
}

# Log-likelihoods of every subject's candidate rows (a list of K-row matrices):
# a K x n_subjects matrix, scored on the worker pool when there is one (same
# transport as the particle step, rows_only mode).
.emc_pooled_score <- function(pmwgs, rows, alpha, wpool, wpool_ctx, wpool_part,
                              pars, group_chol, pm_settings, prev_ll, streams,
                              n_cores) {
  n <- pmwgs$n_subjects
  if (!is.null(wpool)) {
    pars$alpha <- alpha
    it <- .emc_wpool_iter(wpool, wpool_ctx, wpool_part, pars, group_chol,
                          pm_settings, prev_ll, streams, extra = rows,
                          pooled = pmwgs$pooled$value, rows_only = TRUE)
    ll <- it$extra_ll
    alive <- it$alive
  } else {
    ll <- auto_mclapply(seq_len(n), function(s) {
      .emc_pooled_rows_ll(rows[[s]], pmwgs$data[[s]], pmwgs$model)
    }, mc.cores = max(1L, min(n, n_cores)))
    alive <- NULL
  }
  K <- nrow(rows[[1L]])
  bad <- vapply(ll, function(x) !is.numeric(x) || length(x) != K, logical(1))
  # A failed worker reply is recomputed here, where a strict abort is raised.
  for (s in which(bad))
    ll[[s]] <- .emc_pooled_rows_ll(rows[[s]], pmwgs$data[[s]], pmwgs$model)
  ll <- lapply(ll, .emc_pooled_record)
  list(ll = matrix(unlist(ll), nrow = K), alive = alive)
}

.emc_pooled_move <- function(pmwgs, proposals, stage, pars, group_chol,
                             pm_settings, streams, wpool, wpool_ctx, wpool_part,
                             n_cores = 1L) {
  st <- pmwgs$pooled
  nm <- names(st$value)
  p <- length(nm)
  n <- pmwgs$n_subjects
  n_pars <- pmwgs$n_pars
  K <- .emc_pooled_options()$candidates
  alpha <- matrix(proposals[seq_len(n_pars), ], nrow = n_pars,
                  dimnames = list(pmwgs$par_names, NULL))
  prev_ll <- proposals[n_pars + 1L, ]
  score <- function(rows) {
    sc <- .emc_pooled_score(pmwgs, rows, alpha, wpool, wpool_ctx, wpool_part,
                            pars, group_chol, pm_settings, prev_ll, streams, n_cores)
    if (!is.null(wpool) && !is.null(sc$alive)) wpool$alive <<- sc$alive
    sc$ll
  }
  if (!identical(stage, "sample") && !isTRUE(st$curvature_block == st$block)) {
    st <- .emc_pooled_curvature(pmwgs, st, alpha, score)
    st$curvature_block <- st$block
  }
  if (!is.null(pmwgs$rng$gibbs))
    assign(".Random.seed", pmwgs$rng$gibbs, envir = globalenv())
  L <- t(chol(exp(2 * st$log_eps) * st$cov))
  z <- st$value + drop(L %*% stats::rnorm(p))
  cands <- t(z + L %*% matrix(stats::rnorm(p * K), p, K))
  colnames(cands) <- nm
  delta <- sweep(cands, 2L, st$value)
  shift <- if (is.null(st$shear)) matrix(0, K, n_pars) else delta %*% t(st$shear)
  rows <- lapply(seq_len(n), function(s) {
    out <- cbind(sweep(shift, 2L, alpha[, s], "+"), cands)
    colnames(out) <- c(pmwgs$par_names, nm)
    out
  })
  ll <- score(rows)
  lp <- .emc_pooled_log_prior(rbind(st$value, cands), st)
  if (!is.null(st$shear)) {
    mu0 <- as.numeric(pars$tmu)[seq_len(n_pars)]
    lp <- lp + apply(rbind(0, shift), 1L, function(m)
      .emc_pooled_log_prior_mu(mu0 + m, pmwgs$prior))
  }
  lw <- lp + c(sum(prev_ll), rowSums(ll))
  lw[!is.finite(lw)] <- -Inf
  pick <- if (all(lw == -Inf)) 1L else
    sample.int(K + 1L, 1L, prob = exp(lw - max(lw)))
  moved <- pick > 1L
  if (moved) {
    k <- pick - 1L
    st$value[] <- cands[k, ]
    proposals[seq_len(n_pars), ] <- alpha + shift[k, ]
    proposals[n_pars + 1L, ] <- ll[k, ]
    if (!is.null(st$shear)) {
      pars$tmu[seq_len(n_pars)] <- pars$tmu[seq_len(n_pars)] + shift[k, ]
      if (is.matrix(pars$subj_mu)) pars$subj_mu <- pars$subj_mu + shift[k, ]
      pars$alpha <- proposals[seq_len(n_pars), , drop = FALSE]
    }
  }
  st$tries[stage] <- st$tries[stage] + 1
  st$moves[stage] <- st$moves[stage] + moved
  st <- .emc_pooled_adapt(st, stage, moved)
  pmwgs$pooled <- st
  pmwgs$data <- .emc_pooled_set_data(pmwgs$data, st$value)
  pmwgs$rng$gibbs <- get(".Random.seed", envir = globalenv())
  list(sampler = pmwgs, proposals = proposals, pars = pars,
       alive = if (is.null(wpool)) NULL else wpool$alive)
}

# Output: pooled draws as extra rows of alpha (repeated over subjects) and of
# mu. Returns `x` unchanged for samplers without pooled parameters.
.emc_pooled_append <- function(x, chain, idx, selection) {
  pooled <- chain$samples$pooled
  if (is.null(pooled)) return(x)
  pv <- pooled[, idx, drop = FALSE]
  if (selection == "alpha") {
    d <- dim(x)
    out <- array(NA_real_, c(d[1] + nrow(pv), d[2], d[3]),
                 dimnames = list(c(dimnames(x)[[1]], rownames(pv)),
                                 dimnames(x)[[2]], dimnames(x)[[3]]))
    out[seq_len(d[1]), , ] <- x
    for (s in seq_len(d[2])) out[d[1] + seq_len(nrow(pv)), s, ] <- pv
    return(out)
  }
  rbind(x, pv)
}

# The hierarchical rows of an alpha array or mu matrix (pooled rows removed).
.emc_pooled_drop <- function(x, emc) {
  if (is.null(emc[[1]]$pooled)) return(x)
  keep <- dimnames(x)[[1]] %in% emc[[1]]$par_names
  if (length(dim(x)) == 3L) x[keep, , , drop = FALSE] else x[keep, , drop = FALSE]
}

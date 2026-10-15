# Collapsed group move: a move of the whole group state (mu and every free
# entry of V) with the subjects carried in their conditional frames.
#
# Target: pi(mu, V, a | y) propto p(mu) p(V) prod_s L_s(a_s) N(a_s | mu, V).
# Write each subject as a whitened residual in its conditional frame,
#   a_s = m_s + R_s^-1 r_s,
# where N(m_s, (R_s' R_s)^-1) is subject s's conditional posterior under its
# Gaussian likelihood approximation Lhat_s(a) = exp(-a' H_s a / 2 + lin_s' a)
# and prior N(mu, V) (the shear's frame, R/shear_move.R). In the chart
# (phi, r), with phi an unconstrained coordinate of (mu, V), the conditional
# of phi given r is
#   f(phi) propto p(mu) p(V) |dV/dphi| prod_s L_s(a_s) N(a_s | mu, V) / |R_s|
# and, with the Laplace-collapsed group posterior
#   phat(phi) propto p(mu) p(V) |dV/dphi| prod_s Int Lhat_s(a) N(a | mu, V) da,
# the ratio is f / phat = prod_s L_s(a_s) / Lhat_s(a_s) times a function of r
# alone. Under exact Gaussian likelihoods phi and r are independent and phi | r
# is phat: the funnel between the group variance and the subjects is gone,
# and so is the coupling between a log sd and its correlations that pins a
# one-coordinate shear orbit.
#
# The inner chain is compiled (src/collapsed_move.cpp) because it runs in the
# chain's own process while its likelihood workers wait.
#
# The approximations are not the particle step's Laplace fits: a Hessian at
# one conditional mode is accurate near that mode, but its errors add up over
# subjects into a collapsed posterior centred elsewhere, which the outer
# correction then rejects. The acceptance is governed by how much log w
# varies where the posterior lives, so each subject's log Lhat_s is the
# least-squares quadratic of its log likelihood over the window draws (every
# chain), projected to a concave quadratic. Without enough draws for the
# regression the Laplace fit is used.
#
# A move is a surrogate transition (Liu 2001, section 9.4.3): m Hamiltonian
# trajectories on phat, each accepted or rejected against phat, followed by
# one exact Metropolis-Hastings correction
#   min(1, w(phi') / w(phi)),  log w = sum_s [log L_s(a_s) - log Lhat_s(a_s)],
# evaluated at the transported subjects in one likelihood round. The inner
# chain is reversible with respect to phat, so its m-step kernel is too, and
# the correction makes the move reversible with respect to f. Hence every
# move leaves pi invariant for any approximation; a poor one lowers the
# acceptance only. The inner chain needs no likelihood calls.
#
# Chart: phi = (mu, the log-Cholesky factor of every inverse-Wishart block of
# V, the log sd of every gamma-prior parameter). It is unconstrained, so no
# proposal can leave the positive definite cone. The density of V is over its
# free entries (.emc_shear_log_prior_v, with a_half integrated out), and the
# Jacobian of a d x d block is prod_k L_kk^(d - k + 2) (V = L L' with log
# diagonal), of a gamma-prior variance exp(2 phi) it is exp(2 phi). a_half is
# redrawn from its conditional after the move, as in the shear.
#
# Tuning: the fits, the Hamiltonian mass is the window covariance of phi over every
# chain, shrunk towards its diagonal, and the step size is set on phat alone
# (no likelihood calls) by dual averaging from the chain's current state.
# Both change only at proposal rebuilds. Expected moves per iteration are set
# by the move allocation (R/move_allocation.R).
#
# Supported hierarchies are the shear's: every subject's prior mean is the
# group mean (standard without group designs, diagonal-gamma).

# Inner trajectories per move, trajectory length in whitened units, and the
# inner acceptance the step size targets.
.EMC_COLLAPSED_INNER <- 3L
.EMC_COLLAPSED_LENGTH <- 1.5
.EMC_COLLAPSED_TARGET <- 0.8
.EMC_COLLAPSED_PILOT <- 40L
.EMC_COLLAPSED_EPS_RANGE <- c(0.01, 1.5)

.emc_collapsed_active <- function(sampler) .emc_shear_supported(sampler)

# Positions of the chart's coordinates for covariance prior spec.
.emc_collapsed_chart <- function(spec, p) {
  pos <- p
  blocks <- lapply(spec$iw, function(b) {
    d <- length(b)
    lt <- which(lower.tri(diag(d), diag = TRUE))
    at <- pos + seq_along(lt)
    pos <<- pos + length(lt)
    # Position of L_kk in column-major lower-triangle order.
    dg <- vapply(seq_len(d), function(k) which(lt == (k - 1L) * d + k), integer(1))
    list(idx = b, d = d, lt = lt, at = at, diag = dg)
  })
  g <- spec$gamma_idx
  list(p = p, blocks = blocks, gamma = g, gamma_at = pos + seq_along(g),
       dim = pos + length(g))
}

.emc_collapsed_to_phi <- function(mu, V, chart) {
  phi <- numeric(chart$dim)
  phi[seq_len(chart$p)] <- mu
  for (b in chart$blocks) {
    L <- t(chol(V[b$idx, b$idx, drop = FALSE]))
    diag(L) <- log(diag(L))
    phi[b$at] <- L[b$lt]
  }
  if (length(chart$gamma))
    phi[chart$gamma_at] <- 0.5 * log(V[cbind(chart$gamma, chart$gamma)])
  phi
}

# mu, V, the Cholesky factor of each block and log |dV/dphi|.
.emc_collapsed_from_phi <- function(phi, chart) {
  p <- chart$p
  V <- matrix(0, p, p)
  Ls <- vector("list", length(chart$blocks))
  log_jac <- 0
  for (i in seq_along(chart$blocks)) {
    b <- chart$blocks[[i]]
    L <- matrix(0, b$d, b$d)
    L[b$lt] <- phi[b$at]
    ld <- diag(L)
    diag(L) <- exp(ld)
    V[b$idx, b$idx] <- tcrossprod(L)
    Ls[[i]] <- L
    log_jac <- log_jac + sum((b$d - seq_len(b$d) + 2) * ld)
  }
  if (length(chart$gamma)) {
    lg <- phi[chart$gamma_at]
    V[cbind(chart$gamma, chart$gamma)] <- exp(2 * lg)
    log_jac <- log_jac + sum(2 * lg)
  }
  list(mu = phi[seq_len(p)], V = V, L = Ls, log_jac = log_jac)
}

# The chart, approximations and priors as collapsed_surrogate_cpp() and
# collapsed_trajectories_cpp() read them (0-based positions, approximations
# of the subjects that have one).
.emc_collapsed_cpp_model <- function(chart, fits, prior, spec) {
  p <- chart$p
  has <- which(!vapply(fits, is.null, logical(1)))
  H <- array(0, c(p, p, length(has)))
  lin <- matrix(0, p, length(has))
  for (k in seq_along(has)) {
    H[, , k] <- fits[[has[k]]]$prec
    lin[, k] <- fits[[has[k]]]$lin
  }
  list(p = p, dim = chart$dim,
       block_idx = lapply(chart$blocks, function(b) b$idx - 1L),
       block_at = lapply(chart$blocks, function(b) b$at - 1L),
       gamma = as.integer(chart$gamma) - 1L, gamma_at = as.integer(chart$gamma_at) - 1L,
       H = H, lin = lin, m0 = as.numeric(prior$theta_mu_mean),
       P0 = as.matrix(prior$theta_mu_invar), uses_a_half = isTRUE(spec$uses_a_half),
       v = if (is.null(spec$v)) 2 else as.numeric(spec$v),
       A = if (is.null(spec$A)) rep(1, p) else rep_len(as.numeric(spec$A), p),
       shape = if (is.null(spec$shape)) rep(1, p) else rep_len(as.numeric(spec$shape), p),
       rate = if (is.null(spec$rate)) rep(1, p) else rep_len(as.numeric(spec$rate), p))
}

# log phat(phi) with its gradient (src/collapsed_move.cpp); lp is -Inf when
# phi gives no usable state.
.emc_collapsed_surrogate <- function(phi, ctx) {
  out <- collapsed_surrogate_cpp(as.numeric(phi), ctx$cpp)
  if (!is.null(out$grad)) out$grad <- as.numeric(out$grad)
  out
}

# A fitted quadratic can have flat curvature after its negative eigenvalues
# are projected away. A nonzero linear coefficient in that null space makes
# its Gaussian integral unbounded as the group variance grows. Remove those
# unsupported slopes in both empirical and mode-based fits.
.emc_collapsed_safe_fit <- function(fit) {
  if (is.null(fit)) return(NULL)
  H <- fit$prec
  lin <- fit$lin
  if (!is.matrix(H) || nrow(H) != ncol(H) || length(lin) != nrow(H) ||
      any(!is.finite(H)) || any(!is.finite(lin))) return(NULL)
  if (isTRUE(fit$.emc_collapsed_safe)) return(fit)
  H <- (H + t(H)) / 2
  e <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e) || any(!is.finite(e$values))) return(NULL)
  tol <- max(abs(e$values), .Machine$double.xmin) * 1e-10
  keep <- which(e$values > tol)
  U <- e$vectors[, keep, drop = FALSE]
  H <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
  if (ncol(U)) {
    lin <- drop(U %*% crossprod(U, lin))
    if (is.matrix(fit$cross) && nrow(fit$cross) == nrow(H) &&
        all(is.finite(fit$cross)))
      fit$cross <- U %*% crossprod(U, fit$cross)
  } else {
    lin[] <- 0
    if (is.matrix(fit$cross) && nrow(fit$cross) == nrow(H))
      fit$cross[,] <- 0
  }
  fit$prec <- (H + t(H)) / 2
  fit$lin <- lin
  fit$.emc_collapsed_safe <- TRUE
  fit
}

.emc_collapsed_context <- function(sampler) {
  spec <- .emc_shear_v_spec(sampler)
  fits <- sampler$collapsed$fits
  cacheable <- length(fits) == sampler$n_subjects
  if (!cacheable) fits <- .emc_shear_fits(sampler)
  fits <- lapply(fits, .emc_collapsed_safe_fit)
  base_fits <- fits
  fits <- .emc_collapsed_condition(fits, sampler$pooled$value)
  fits <- lapply(fits, .emc_collapsed_safe_fit)
  chart <- .emc_collapsed_chart(spec, sampler$n_pars)
  list(fits = fits, base_fits = base_fits, cacheable = cacheable,
       spec = spec, prior = sampler$prior, chart = chart,
       walls = .emc_shear_walls_of(sampler),
       cpp = .emc_collapsed_cpp_model(chart, fits, sampler$prior, spec))
}

# Least-squares quadratic of a subject's log likelihood y over its draws A
# (rows), in canonical form (prec, lin). Fitted in coordinates whitened by
# the draws' covariance, then projected to a concave quadratic with the
# linear part refitted. NULL without enough distinct draws.
.emc_collapsed_fit_subject <- function(A, y) {
  ok <- is.finite(y) & stats::complete.cases(A)
  A <- A[ok, , drop = FALSE]
  y <- y[ok]
  p <- ncol(A)
  iu <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  if (nrow(A) < 2L * (1L + p + nrow(iu))) return(NULL)
  abar <- colMeans(A)
  e <- eigen(stats::cov(A), symmetric = TRUE)
  if (any(!is.finite(e$values)) || any(e$values <= 1e-12 * max(e$values))) return(NULL)
  W <- t(e$vectors) / sqrt(e$values)
  Z <- sweep(A, 2L, abar) %*% t(W)
  X <- cbind(1, Z, Z[, iu[, 1L], drop = FALSE] * Z[, iu[, 2L], drop = FALSE])
  b <- tryCatch(qr.coef(qr(X), y), error = function(e) NULL)
  if (is.null(b) || anyNA(b)) return(NULL)
  # y ~ c + g'z - z'Bz / 2: a product coefficient is -B_ij off the diagonal
  # and -B_ii / 2 on it.
  cq <- b[-seq_len(p + 1L)]
  B <- matrix(0, p, p)
  B[iu] <- -cq
  B <- B + t(B)
  diag(B) <- -2 * cq[iu[, 1L] == iu[, 2L]]
  eb <- eigen(B, symmetric = TRUE)
  Bp <- eb$vectors %*% (pmax(eb$values, 0) * t(eb$vectors))
  g <- qr.coef(qr(cbind(1, Z)), y + 0.5 * rowSums((Z %*% Bp) * Z))[-1L]
  # PSD projection can create flat directions. Refitting an unconstrained
  # linear term after that projection gives those directions arbitrary slopes.
  tol <- max(abs(eb$values), .Machine$double.xmin) * 1e-10
  keep <- which(eb$values > tol)
  if (length(keep)) {
    U <- eb$vectors[, keep, drop = FALSE]
    g <- drop(U %*% crossprod(U, g))
  } else {
    g[] <- 0
  }
  H <- t(W) %*% Bp %*% W
  H <- (H + t(H)) / 2
  lin <- drop(H %*% abar + t(W) %*% g)
  if (!all(is.finite(H)) || !all(is.finite(lin))) return(NULL)
  list(prec = H, lin = lin)
}

# Every subject's fit from the window draws of every chain, else its Laplace
# fit. Each draw's likelihood was taken at that draw's pooled coefficients
# (design(pooled = )), so with pooled coefficients the quadratic is fitted in
# (alpha, pooled) jointly and keeps the cross block, from which
# .emc_collapsed_condition() gives the quadratic in alpha at the current
# pooled value: the one the move's likelihoods are evaluated at.
.emc_collapsed_fits <- function(emc, idx, n_cores = 1L) {
  ch <- emc[[1L]]
  laplace <- .emc_shear_fits(ch)
  p <- ch$n_pars
  pooled <- if (!is.null(ch$pooled) && is.matrix(ch$samples$pooled))
    do.call(rbind, lapply(emc, function(c) t(c$samples$pooled[, idx, drop = FALSE])))
  # Wide quadratic regressions justify forking; small regressions stay in
  # the caller to avoid paying process overhead for a few matrix operations.
  d <- p + if (is.null(pooled)) 0L else ncol(pooled)
  coefficients <- 1L + d + d * (d + 1L) / 2L
  draws <- length(idx) * length(emc)
  enough <- draws >= 2L * coefficients
  work <- draws * coefficients^2 * ch$n_subjects
  cores <- if (d >= 8L && enough && ch$n_subjects >= 8L && work >= 1e8)
    min(n_cores, ch$n_subjects %/% 2L) else 1L
  .emc_rebuild_lapply(seq_len(ch$n_subjects), function(s) {
    A <- do.call(rbind, lapply(emc, function(c) t(matrix(c$samples$alpha[, s, idx], p))))
    y <- unlist(lapply(emc, function(c) c$samples$subj_ll[s, idx]))
    fit <- tryCatch(.emc_collapsed_fit_subject(cbind(A, pooled), y), error = function(e) NULL)
    if (is.null(fit)) return(laplace[[s]])
    if (is.null(pooled)) return(fit)
    a <- seq_len(p)
    list(prec = fit$prec[a, a, drop = FALSE], lin = fit$lin[a],
         cross = fit$prec[a, -a, drop = FALSE])
  }, n_cores = cores)
}

# The fits as quadratics in alpha at pooled coefficients `value`: the joint
# quadratic -x'Hx/2 + lin'x in x = (alpha, pooled) has linear part
# lin_a - H_ap value in alpha.
.emc_collapsed_condition <- function(fits, value) {
  lapply(fits, function(f) {
    if (is.null(f$cross) || is.null(value)) return(f)
    list(prec = f$prec, lin = f$lin - drop(f$cross %*% as.numeric(value)),
         .emc_collapsed_safe = isTRUE(f$.emc_collapsed_safe))
  })
}

# log Lhat_s at each subject's state; 0 without a fit.
.emc_collapsed_log_lhat <- function(alpha, fits) {
  vapply(seq_along(fits), function(s) {
    fit <- fits[[s]]
    if (is.null(fit)) return(0)
    a <- alpha[, s]
    -0.5 * sum(a * (fit$prec %*% a)) + sum(fit$lin * a)
  }, numeric(1))
}

# Every chain's plan from the proposal window create_chain_proposals() uses:
# the mass from every chain's window draws, then each chain's step size on
# phat from its current state.
.emc_collapsed_set_plans <- function(emc, samples_idx, stage, n_cores = 1L) {
  if (!.emc_collapsed_active(emc[[1L]])) return(emc)
  idx <- .emc_chain_proposal_window(emc, samples_idx, stage)
  if (length(idx) < 2L) return(emc)
  fits <- .emc_collapsed_fits(emc, idx, n_cores)
  for (j in seq_along(emc)) emc[[j]]$collapsed$fits <- fits
  ctx <- .emc_collapsed_context(emc[[1L]])
  chart <- ctx$chart
  X <- do.call(rbind, lapply(emc, function(ch) t(vapply(idx, function(t)
    tryCatch(.emc_collapsed_to_phi(ch$samples$theta_mu[, t],
                                   as.matrix(ch$samples$theta_var[, , t]), chart),
             error = function(e) rep(NA_real_, chart$dim)), numeric(chart$dim)))))
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < 3L) return(emc)
  S <- stats::cov(X)
  lam <- min(1, chart$dim / nrow(X))
  M <- (1 - lam) * S + lam * diag(pmax(diag(S), 1e-8), chart$dim)
  C <- tryCatch(t(chol(M)), error = function(e) NULL)
  if (is.null(C)) return(emc)
  for (j in seq_along(emc)) emc[[j]] <- .emc_collapsed_set_plan(emc[[j]], C, stage)
  emc
}

.emc_collapsed_set_plan <- function(chain, C, stage) {
  ctx <- .emc_collapsed_context(chain)
  p <- chain$n_pars
  st <- chain$collapsed
  if (is.null(st)) st <- list()
  old <- st$plan
  idx <- chain$samples$idx
  phi <- tryCatch(.emc_collapsed_to_phi(chain$samples$theta_mu[, idx],
                                        as.matrix(chain$samples$theta_var[, , idx]),
                                        ctx$chart), error = function(e) NULL)
  if (is.null(phi)) return(chain)
  if (!is.finite(.emc_collapsed_surrogate(phi, ctx)$lp)) return(chain)
  if (!is.null(chain$rng$gibbs)) assign(".Random.seed", chain$rng$gibbs, envir = globalenv())
  eps0 <- if (is.list(old) && is.finite(old$eps)) old$eps else 0.5 * ctx$chart$dim^(-0.25)
  eps <- .emc_collapsed_tune_eps(phi, C, eps0, ctx)
  chain$rng$gibbs <- get(".Random.seed", envir = globalenv())
  zero <- .emc_collapsed_record(NULL, p)
  # Discount movement evidence across rebuilds so the exploration floor
  # re-probes moves at a diminishing rate.
  blk <- .emc_collapsed_record(if (is.list(old$block)) old$block, p)
  blk$iters <- .emc_alloc_block_iters(chain)
  same_stage <- identical(stage, old$stage)
  production_pool <- identical(stage, "sample") && same_stage
  prev <- if (production_pool || !identical(stage, "sample")) {
    .emc_collapsed_record(if (same_stage) old$pool else old$last, p)
  } else zero
  if (!same_stage && !identical(stage, "sample") && prev$n < 5L) prev <- zero
  if (production_pool || !identical(stage, "sample")) {
    d <- .EMC_ALLOC_EVIDENCE_DECAY
    pool <- list(E = d * prev$E + blk$E, n = d * prev$n + blk$n,
                 acc = d * prev$acc + blk$acc, moved = d * prev$moved + blk$moved,
                 N = d * prev$N + blk$N,
                 K = .emc_alloc_khist_sum(d * prev$K, blk$K),
                 iters = d * prev$iters + blk$iters)
    last <- pool
  } else {
    pool <- zero
    last <- blk
  }
  st$plan <- list(C = C, eps = eps, rounds = 0, dim = ctx$chart$dim, stage = stage,
                  last = last, pool = pool, block = zero)
  chain$collapsed <- st
  chain
}

# A block's move record: the acceptance-weighted squared jumps E of its n
# moves (acceptance sum acc, moved of them moving), and the crossproduct N of
# each iteration's net jump over its moves with K the iterations' move counts
# (one row), and the iterations the record spans (iters, zero where
# unknown). Fills in the fields a record made without them lacks.
.emc_collapsed_record <- function(rec, p) {
  if (is.null(rec)) rec <- list(E = matrix(0, 2L * p, 2L * p), n = 0L, acc = 0, moved = 0L)
  if (is.null(rec$N)) rec$N <- matrix(0, 2L * p, 2L * p)
  if (is.null(rec$K)) rec$K <- matrix(0L, 1L, 1L)
  if (is.null(rec$iters)) rec$iters <- 0
  rec
}

# Dual averaging (Hoffman & Gelman 2014) of the step size towards the target
# inner acceptance, on phat alone, from the chain's current state.
.emc_collapsed_tune_eps <- function(phi, C, eps0, ctx) {
  mu_da <- log(10 * eps0)
  h_bar <- 0
  log_eps <- log(eps0)
  log_bar <- log_eps
  for (k in seq_len(.EMC_COLLAPSED_PILOT)) {
    tr <- collapsed_trajectories_cpp(phi, C, exp(log_eps), 1L, .EMC_COLLAPSED_LENGTH, ctx$cpp)
    phi <- tr$phi
    h_bar <- (1 - 1 / (k + 10)) * h_bar + (.EMC_COLLAPSED_TARGET - tr$acc) / (k + 10)
    log_eps <- mu_da - sqrt(k) / 0.05 * h_bar
    log_eps <- min(log(.EMC_COLLAPSED_EPS_RANGE[2]), max(log(.EMC_COLLAPSED_EPS_RANGE[1]), log_eps))
    wk <- k^(-0.75)
    log_bar <- wk * log_eps + (1 - wk) * log_bar
  }
  min(.EMC_COLLAPSED_EPS_RANGE[2], max(.EMC_COLLAPSED_EPS_RANGE[1], exp(log_bar)))
}

# One move from x (group state, subjects and their likelihoods, worker pool).
.emc_collapsed_round <- function(x, plan, sampler, ctx) {
  p <- nrow(x$alpha)
  n <- ncol(x$alpha)
  x$acc <- 0
  x$jump <- numeric(2L * p)
  cur <- x$cur
  r <- matrix(vapply(seq_len(n), function(s) {
    drop(cur$fm$R[, , s] %*% (x$alpha[, s] - cur$fm$m[, s]))
  }, numeric(p)), p)
  phi0 <- .emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
  inner <- collapsed_trajectories_cpp(phi0, plan$C, plan$eps, .EMC_COLLAPSED_INNER,
                                      .EMC_COLLAPSED_LENGTH, ctx$cpp)
  if (!is.finite(inner$lp)) return(x)
  phi <- as.numeric(inner$phi)
  if (identical(phi, phi0)) return(x)
  st <- .emc_collapsed_from_phi(phi, ctx$chart)
  new <- .emc_shear_state(st$mu, st$V, ctx$fits, sampler$prior, ctx$spec)
  if (is.null(new) || !is.finite(new$lp)) return(x)
  walls <- ctx$walls
  # Subjects with a transport wall take the wall-aware map; the exact
  # correction then carries the map's Jacobian in place of the plain one and
  # the change of their whitened coordinate (see .emc_shear_wall_map).
  tr <- .emc_shear_transport(x$alpha, cur, new, walls)
  if (!all(tr$ok) || any(!is.finite(tr$alpha))) return(x)
  alpha_new <- tr$alpha
  jac <- if (is.null(walls)) numeric(n) else
    ifelse(walls$has, tr$jac - 0.5 * (tr$r1sq - colSums(r^2)), 0)
  pn <- sampler$par_names
  rows <- lapply(seq_len(n), function(s)
    matrix(alpha_new[, s], 1L, dimnames = list(NULL, pn)))
  evaluated <- .emc_wpool_subject_ll(x$wpool, rows, sampler, x$wpool_part,
                                    source = "collapsed")
  x$wpool <- evaluated$pool
  x$rounds <- x$rounds + 1L
  ll_new <- as.numeric(evaluated$ll)
  lhat_new <- .emc_collapsed_log_lhat(alpha_new, ctx$fits)
  lhat_old <- .emc_collapsed_log_lhat(x$alpha, ctx$fits)
  la <- sum(ll_new - lhat_new) - sum(x$ll - lhat_old) + sum(jac)
  acc <- if (is.finite(la) && all(is.finite(ll_new))) min(1, exp(la)) else 0
  x$acc <- acc
  x$jump <- c(new$mu - cur$mu, .emc_alloc_lsd(new$V) - .emc_alloc_lsd(cur$V))
  if (stats::runif(1L) < acc) {
    x$cur <- new
    x$alpha <- alpha_new
    x$ll <- ll_new
    x$moved <- TRUE
  }
  x
}

.emc_collapsed_move <- function(sampler, pars, proposals, stage, wpool = NULL,
                                wpool_part = NULL) {
  keep <- list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool,
               rounds = 0L)
  if (stage == "preburn" || !.emc_collapsed_active(sampler)) return(keep)
  st <- sampler$collapsed
  plan <- st$plan
  p <- sampler$n_pars
  if (!is.list(plan) || !is.matrix(plan$C) || !isTRUE(plan$rounds > 0)) return(keep)
  if (is.null(st$iterations)) st$iterations <- 0L
  if (is.null(st$rounds)) st$rounds <- 0L
  if (is.null(st$accepted)) st$accepted <- 0L
  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  # The allocated expected moves: the whole part, plus one more with the
  # fractional part as probability.
  whole <- floor(plan$rounds)
  n_rounds <- whole + (stats::runif(1L) < plan$rounds - whole)
  st$iterations <- st$iterations + 1L
  if (n_rounds == 0L) {
    sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
    sampler$collapsed <- st
    keep$sampler <- sampler
    return(keep)
  }
  alpha_start <- matrix(proposals[seq_len(p), ], p)
  ctx <- .emc_collapsed_context(sampler)
  # The base fits change at proposal rebuilds, but pooled conditioning changes
  # every iteration. Cache only the sanitized base fits and recondition them
  # against the current pooled values in each context.
  if (isTRUE(ctx$cacheable)) st$fits <- ctx$base_fits
  cur <- if (n_rounds > 0L && plan$dim == ctx$chart$dim)
    .emc_shear_state(as.numeric(pars$tmu), as.matrix(pars$tvar), ctx$fits,
                     sampler$prior, ctx$spec)
  x <- list(cur = cur, alpha = alpha_start,
            ll = proposals[p + 1L, ], rounds = 0L, wpool = wpool,
            wpool_part = wpool_part)
  any_moved <- FALSE
  round_el <- 0
  if (!is.null(cur) && is.finite(cur$lp)) {
    blk <- .emc_collapsed_record(plan$block, p)
    round_el0 <- proc.time()[["elapsed"]]
    for (k in seq_len(n_rounds)) {
      x$moved <- FALSE
      x <- .emc_collapsed_round(x, plan, sampler, ctx)
      blk$E <- blk$E + x$acc * tcrossprod(x$jump)
      blk$n <- blk$n + 1L
      blk$acc <- blk$acc + x$acc
      blk$moved <- blk$moved + as.integer(x$moved)
      any_moved <- any_moved || x$moved
    }
    round_el <- proc.time()[["elapsed"]] - round_el0
    # The net jump over this iteration's moves (their saturation,
    # R/move_allocation.R).
    d <- c(x$cur$mu - cur$mu, .emc_alloc_lsd(x$cur$V) - .emc_alloc_lsd(cur$V))
    if (n_rounds > 0L && all(is.finite(d))) {
      blk$N <- blk$N + tcrossprod(d)
      blk$K <- .emc_alloc_khist_add(blk$K, 1L, n_rounds)
    }
    plan$block <- blk
    st$plan <- plan
    st$rounds <- st$rounds + x$rounds
  }
  if (any_moved) {
    st$accepted <- st$accepted + 1L
    if (ctx$spec$uses_a_half)
      pars$a_half <- stats::setNames(.emc_shear_draw_a_half(x$cur$V, ctx$spec),
                                     names(pars$a_half))
    pars$tmu <- stats::setNames(x$cur$mu, names(pars$tmu))
    pars$tvar <- `dimnames<-`(x$cur$V, dimnames(pars$tvar))
    pars$tvinv <- `dimnames<-`(x$cur$Vi, dimnames(pars$tvar))
    # Every subject's prior mean is the group mean here (see
    # .emc_shear_supported); updates later in the iteration read it.
    if (is.matrix(pars$subj_mu)) pars$subj_mu[] <- x$cur$mu
    proposals[seq_len(p), ] <- x$alpha
    proposals[p + 1L, ] <- x$ll
    pars$alpha <- x$alpha
  }
  sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  sampler$collapsed <- st
  list(sampler = sampler, pars = pars, proposals = proposals, wpool = x$wpool,
       rounds = x$rounds, round_el = round_el)
}

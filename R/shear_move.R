# Shear group update: the ensemble update (R/ensemble_move.R) followed by a
# shear of the group state.
#
# Target: pi(mu, V, a | y) propto p(mu) p(V) prod_s L_s(a_s) N(a_s | mu, V).
# The ensemble update draws the group state given the subjects (the centred
# augmentation). In a direction the data inform weakly it can only move the
# group as far as the subjects already have, and vice versa: the group mean
# along a ridge, and the group sd of a parameter that each subject's data
# barely determine (the hierarchical funnel). The shear moves the group
# state while each subject keeps its standardised position in its own
# conditional (the non-centred augmentation). Alternating the two
# augmentations of one group state is the interweaving strategy of Yu & Meng
# (2011); each step leaves the posterior invariant, so the composition does.
#
# The shear map, for a step u = (u_mu, u_sd):
#   mu -> mu + u_mu,   V -> D V D,   D = diag(exp(u_sd)),
#   a_s -> m_s' + R_s'^-1 R_s (a_s - m_s),
# where N(m_s, (R_s' R_s)^-1) is subject s's conditional posterior under its
# Gaussian likelihood approximation (lik_laplace, precision H_s, see
# .emc_lik_canonical) and prior N(mu, V), and primes denote the moved group
# state. A well-informed subject barely moves; a prior-dominated one is
# shifted and rescaled with the group. Without a usable approximation
# H_s = 0, so a_s -> mu' + D (a_s - mu). The moved subjects depend only on
# the start and end group states, so the maps compose, T_v T_u = T_(u+v),
# and T_-u inverts T_u. The Jacobian of T_u is
#   J_u = exp(sum_i (g_i + 1) u_sd,i) prod_s |R_s| / |R_s'|
# (g_i is the size of the covariance block holding parameter i; the prior
# density of V is taken over the free entries of its blocks), and it composes
# the same way. The approximation shapes the map only; a wrong one makes the
# moves smaller, not invalid.
#
# A round moves along one coordinate: a mean direction d (u = t (d, 0)) or
# the log group sd of one parameter (u = t e_i). Restricted to the orbit
# t -> T_tv(x) the posterior has density f(t) = pi(T_tv x) J_tv(x), and the
# orbit through any point of it gives the same f up to a shift in t. A round
# is the multiple-proposal step of Tjelmeland (2004) and Calderhead (2014) on
# t: an auxiliary point t~ ~ N(0, c^2), N candidates t_j ~ N(t~, c^2), all
# evaluated in one likelihood round, and one of the N + 1 points (the
# current one is t_0 = 0) chosen with probability proportional to f(t_j).
# Given t~ the N + 1 points are exchangeable, so the choice is a Gibbs step
# on their index and leaves the posterior invariant. Per round, N = 8
# candidates halve the autocorrelation time of a single Metropolis step at
# about the same cost, because a likelihood call's fixed cost dwarfs its
# extra rows.
#
# Which coordinates. Under the approximations the centred Gibbs step moves
# the group mean as an autoregression with matrix
#   A = (n V^-1 + P0)^-1 V^-1 (sum_s Q_s^-1) V^-1,  Q_s = H_s + V^-1,
# whose eigenvectors are the mean directions (columns of dirs, each moving mu
# by one posterior sd of its direction under the approximations) and whose
# eigenvalues lambda_i predict their lag-1 autocorrelations; for a group sd,
# w_i = mean_s (Q_s^-1 V^-1)_ii is the prior share of each subject's
# conditional precision. How many rounds each coordinate gets per iteration
# is set by the move allocation (R/move_allocation.R) from what the rest of
# the iteration already mixes, at proposal rebuilds, and frozen in between.
# Mean steps keep V, and so every R_s: only the centres move.
#
# Under the Huang-Wand prior (standard type), V | a ~ IW(v + d - 1,
# 2 v diag(1 / a)) with a_k ~ IG(1/2, 1 / A_k^2). Holding a fixed would pin
# a scale step to the prior's conditional width in log sd whatever the data
# say. Given V the inverse-Wishart density factorises over the a_k, so a
# integrates out in closed form:
#   log p(V) = -(df + d + 1)/2 log|V| - (v + d)/2 sum_k log(v (V^-1)_kk + 1/A_k^2)
# (for an unblocked parameter, the d = 1 case on its variance). The rounds
# target the posterior with a integrated out, and a is then drawn from its
# conditional given the final V, the draw the Gibbs step makes: together a
# valid update of (V, a).
#
# Supported hierarchies: those in which every subject's prior mean is the
# group mean (standard without group designs, diagonal-gamma); the others
# keep the ensemble update alone.

# Candidates per round, and the share of rounds that should move (the
# optimum of the step scale c for N = 8 on a Gaussian orbit is flat around
# this value).
.EMC_SHEAR_CANDIDATES <- 8L
.EMC_SHEAR_MOVE_TARGET <- 0.68
# Initial step scale in posterior sds of the coordinate, and the bounds every
# adapted scale stays in.
.EMC_SHEAR_STEP0 <- 2.5
.EMC_SHEAR_LOG_STEP_RANGE <- log(c(0.05, 20))

.emc_shear_options <- function() {
  list(enabled = identical(getOption("emc2.group_move", "shear"), "shear"))
}

.emc_shear_supported <- function(sampler) {
  isTRUE(sampler$type %in% c("standard", "diagonal-gamma")) &&
    !any(sampler$nuisance) &&
    is.null(sampler$gd) &&
    length(sampler$group_designs) == 0L &&
    isTRUE(sampler$n_subjects >= 2L)
}

.emc_shear_active <- function(sampler) {
  .emc_shear_options()$enabled && .emc_shear_supported(sampler)
}

# Each subject's likelihood approximation if it is finite, else NULL. The
# particle step's Pareto check does not apply: it judges the Gaussian
# conditional at one group state, while the map is valid with any
# approximation and moves an informative subject far better with a rough one
# than with none (H_s = 0).
.emc_shear_fits <- function(sampler) {
  p <- sampler$n_pars
  lapply(seq_len(sampler$n_subjects), function(s) {
    fit <- if (length(sampler$lik_laplace) >= s) sampler$lik_laplace[[s]] else NULL
    ok <- is.list(fit) && is.matrix(fit$prec) && all(dim(fit$prec) == p) &&
      length(fit$lin) == p && all(is.finite(fit$prec)) && all(is.finite(fit$lin))
    if (ok) fit else NULL
  })
}

# Mean directions and scale weights at group covariance V (see the header).
# Columns of dirs move mu by one posterior sd of their direction under the
# approximations: mu = mu0 + t d with d' (n V^-1 + P0) d = 1 has collapsed
# precision 1 - lambda in t.
.emc_shear_plan <- function(fits, prior, V) {
  p <- nrow(V)
  n <- length(fits)
  Vi <- chol2inv(chol(V))
  S <- matrix(0, p, p)
  w <- numeric(p)
  for (fit in fits) {
    Qi <- chol2inv(chol(if (is.null(fit)) Vi else fit$prec + Vi))
    S <- S + Qi
    w <- w + rowSums(Qi * t(Vi))
  }
  Pc <- n * Vi + as.matrix(prior$theta_mu_invar)
  R <- chol((Pc + t(Pc)) / 2)
  Ri <- backsolve(R, diag(p))
  M <- t(Ri) %*% (Vi %*% S %*% Vi) %*% Ri
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  lam <- pmin(pmax(e$values, 0), 1 - 1e-6)
  dirs <- sweep(Ri %*% e$vectors, 2L, sqrt(1 - lam), "/")
  if (!all(is.finite(dirs))) return(NULL)
  list(dirs = dirs, lam = lam, w = pmin(1, pmax(0, w / n)))
}

# The move along coordinate k (1..p mean directions, p+1..2p log group sds)
# as a step direction in u = (u_mu, u_sd).
.emc_shear_coord_dir <- function(plan, k) {
  p <- nrow(plan$dirs)
  if (k <= p) c(plan$dirs[, k], numeric(p)) else replace(numeric(2L * p), k, 1)
}

# Every chain's plan from the proposal window create_chain_proposals() uses.
# Step scales and the last block's round records carry over to the new
# directions they match (|cosine| >= 0.9; a log group sd always matches
# itself); the round counts are set by the move allocation afterwards.
.emc_shear_set_plans <- function(emc, samples_idx, stage) {
  if (!.emc_shear_active(emc[[1L]])) return(emc)
  idx <- .emc_chain_proposal_window(emc, samples_idx, stage)
  if (!length(idx)) return(emc)
  for (j in seq_along(emc)) emc[[j]] <- .emc_shear_set_plan(emc[[j]], idx, stage)
  emc
}

.emc_shear_set_plan <- function(chain, idx, stage) {
  if (!.emc_shear_active(chain)) return(chain)
  n <- chain$n_subjects
  p <- chain$n_pars
  V <- apply(chain$samples$theta_var[, , idx, drop = FALSE], 1:2, mean)
  plan <- tryCatch(.emc_shear_plan(.emc_shear_fits(chain), chain$prior, V),
                   error = function(e) NULL)
  if (is.null(plan)) return(chain)
  st <- chain$shear
  if (is.null(st)) st <- list()
  old <- st$plan
  base <- c(rep(0, p), rep(-0.5 * log(2 * n), p))
  plan$log_step <- base + log(.EMC_SHEAR_STEP0)
  plan$rounds <- numeric(2L * p)
  plan$stage <- stage
  zero <- .emc_shear_record(NULL, p)
  recent <- plan$pool <- zero
  # Production blocks run settled kernels on a stationary chain, so their
  # rounds pool over the sample stage (from its first block on).
  pooled <- identical(stage, "sample") && identical(old$stage, "sample")
  if (is.list(old) && is.matrix(old$dirs) && all(dim(old$dirs) == p)) {
    scale <- c(rep(NA_real_, p), rep(1, p))
    match <- c(rep(NA_integer_, p), p + seq_len(p))
    for (k in seq_len(p)) {
      d <- plan$dirs[, k]
      cs <- abs(drop(crossprod(old$dirs, d))) /
        (sqrt(colSums(old$dirs^2)) * sqrt(sum(d^2)))
      m <- which.max(cs)
      if (cs[m] >= 0.9) {
        match[k] <- m
        scale[k] <- sqrt(sum(old$dirs[, m]^2) / sum(d^2))
      }
    }
    got <- which(!is.na(match))
    blk <- old$block
    plan$log_step[got] <- old$log_step[match[got]] + log(scale[got])
    carry <- function(rec) {
      out <- zero
      if (is.null(rec)) return(out)
      rec <- .emc_shear_record(rec, p)
      out$E[got] <- rec$E[match[got]] * scale[got]^2
      out$n[got] <- rec$n[match[got]]
      out$moved[got] <- rec$moved[match[got]]
      out$N[got] <- rec$N[match[got]] * scale[got]^2
      out$K <- matrix(0L, 2L * p, ncol(rec$K))
      out$K[got, ] <- rec$K[match[got], ]
      out
    }
    recent <- carry(blk)
    if (pooled) {
      prev <- carry(old$pool)
      plan$pool <- list(E = prev$E + recent$E, n = prev$n + recent$n,
                        moved = prev$moved + recent$moved, N = prev$N + recent$N,
                        K = .emc_alloc_khist_sum(prev$K, recent$K))
    }
  }
  # The allocation measures a coordinate's rounds from the pooled production
  # blocks, or else the last block alone. Step scales follow the last block.
  plan$last <- if (pooled) plan$pool else recent
  # Discarded stages adapt the scales round by round; production only here.
  if (identical(stage, "sample")) plan <- .emc_shear_tune_steps(plan, recent)
  plan$block <- zero
  st$plan <- plan
  chain$shear <- st
  chain
}

# A block's round record: per coordinate the selection-weighted squared
# moves E over n rounds (moved of them moving), and the net squared move N of
# its rounds per iteration with K the iterations' round counts. Fills in the
# net-move fields a record made without them lacks.
.emc_shear_record <- function(rec, p) {
  if (is.null(rec))
    rec <- list(E = numeric(2L * p), n = integer(2L * p), moved = integer(2L * p))
  if (is.null(rec$N)) rec$N <- numeric(2L * p)
  if (is.null(rec$K)) rec$K <- matrix(0L, 2L * p, 1L)
  rec
}

# A production rebuild moves each step scale towards the target share of
# moving rounds, from the last block's rounds `lst` (rounds within a block
# never change it).
.emc_shear_tune_steps <- function(plan, lst) {
  ok <- lst$n >= 10L
  if (any(ok)) {
    rate <- lst$moved[ok] / lst$n[ok]
    plan$log_step[ok] <- plan$log_step[ok] +
      pmax(-log(2), pmin(log(2), 2 * (rate - .EMC_SHEAR_MOVE_TARGET)))
  }
  plan$log_step <- pmin(.EMC_SHEAR_LOG_STEP_RANGE[2],
                        pmax(.EMC_SHEAR_LOG_STEP_RANGE[1], plan$log_step))
  plan
}

# Covariance prior as each variant's Gibbs step targets it: diagonal-gamma
# has a gamma prior on every precision; standard has an inverse-Wishart per
# blocked parameter group and a gamma prior on the precision of every
# unblocked parameter, given a_half, which is integrated out (see the header).
.emc_shear_v_spec <- function(sampler) {
  p <- sampler$n_pars
  prior <- sampler$prior
  if (identical(sampler$type, "diagonal-gamma")) {
    return(list(iw = list(), gamma_idx = seq_len(p), size = rep(1L, p),
                shape = rep_len(prior$shape, p), rate = rep_len(prior$rate, p),
                uses_a_half = FALSE))
  }
  is_blocked <- sampler$is_blocked[!sampler$nuisance]
  par_group <- sampler$par_group[!sampler$nuisance]
  if (is.null(is_blocked) || length(is_blocked) == 0L) is_blocked <- rep(TRUE, p)
  if (is.null(par_group)) par_group <- rep(1, p)
  blocked <- which(is_blocked)
  iw <- lapply(unique(par_group[blocked]), function(g) blocked[par_group[blocked] == g])
  size <- rep(1L, p)
  for (b in iw) size[b] <- length(b)
  list(iw = iw, gamma_idx = which(!is_blocked), size = size,
       v = prior$v, A = rep_len(prior$A, p), uses_a_half = TRUE)
}

# Log density of V over the free entries of its blocks.
.emc_shear_log_prior_v <- function(V, spec) {
  lp <- 0
  g <- spec$gamma_idx
  if (length(g)) {
    tau <- 1 / V[cbind(g, g)]
    # Density of the variance: the precision density times tau^2.
    lp <- lp + if (spec$uses_a_half) {
      sum((spec$v / 2 + 1) * log(tau) -
            (spec$v + 1) / 2 * log(spec$v * tau + 1 / spec$A[g]^2))
    } else {
      sum((spec$shape[g] + 1) * log(tau) - spec$rate[g] * tau)
    }
  }
  for (b in spec$iw) {
    d <- length(b)
    R <- chol(V[b, b, drop = FALSE])
    df <- spec$v + d - 1
    lp <- lp - (df + d + 1) * sum(log(diag(R))) -
      (spec$v + d) / 2 * sum(log(spec$v * diag(chol2inv(R)) + 1 / spec$A[b]^2))
  }
  lp
}

# a_half given V, as the standard Gibbs step draws it:
# IG((v + d_k) / 2, v (V^-1)_kk + 1 / A_k^2), with d_k the size of k's block.
.emc_shear_draw_a_half <- function(V, spec) {
  prec <- numeric(nrow(V))
  g <- spec$gamma_idx
  if (length(g)) prec[g] <- 1 / V[cbind(g, g)]
  for (b in spec$iw) prec[b] <- diag(chol2inv(chol(V[b, b, drop = FALSE])))
  1 / stats::rgamma(nrow(V), shape = (spec$v + spec$size) / 2,
                    rate = spec$v * prec + 1 / spec$A^2)
}

.emc_shear_log_prior_mu <- function(mu, prior) {
  d <- mu - prior$theta_mu_mean
  -0.5 * drop(crossprod(d, prior$theta_mu_invar %*% d))
}

# Everything the move needs about one group state (mu, V): its prior density
# and each subject's conditional frame (m_s, R_s), with R_s the upper Cholesky
# factor of the conditional precision H_s + V^-1.
.emc_shear_state <- function(mu, V, fits, prior, spec) {
  chol_V <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(chol_V)) return(NULL)
  Vi <- chol2inv(chol_V)
  Vi_mu <- drop(Vi %*% mu)
  frames <- lapply(fits, function(fit) {
    Q <- if (is.null(fit)) Vi else fit$prec + Vi
    R <- tryCatch(chol(Q), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    rhs <- if (is.null(fit)) Vi_mu else fit$lin + Vi_mu
    list(m = drop(backsolve(R, forwardsolve(t(R), rhs))), R = R,
         logdet = sum(log(diag(R))))
  })
  if (any(vapply(frames, is.null, logical(1)))) return(NULL)
  lp_v <- tryCatch(.emc_shear_log_prior_v(V, spec), error = function(e) NA_real_)
  list(mu = mu, V = V, Vi = Vi, chol_V = chol_V, frames = frames,
       logdet = sum(vapply(frames, `[[`, 0, "logdet")), lp_v = lp_v,
       lp = .emc_shear_log_prior_mu(mu, prior) + lp_v)
}

# The state after a mean step u_mu: V, and so every R_s, is unchanged and
# each centre moves by Q_s^-1 V^-1 u_mu.
.emc_shear_shift_state <- function(cur, u_mu, prior) {
  Vi_u <- drop(cur$Vi %*% u_mu)
  new <- cur
  new$frames <- lapply(cur$frames, function(f) {
    f$m <- f$m + drop(backsolve(f$R, forwardsolve(t(f$R), Vi_u)))
    f
  })
  new$mu <- cur$mu + u_mu
  new$lp <- .emc_shear_log_prior_mu(new$mu, prior) + cur$lp_v
  new
}

# sum_s log N(a_s | mu, V).
.emc_shear_log_hier <- function(alpha, state) {
  z <- backsolve(state$chol_V, alpha - state$mu, transpose = TRUE)
  -0.5 * sum(z^2) - ncol(alpha) * sum(log(diag(state$chol_V)))
}

# The map for a given u, from the group state cur and subjects alpha: the
# moved state, the moved subjects and log|det| of the map's Jacobian. NULL
# when the moved state is invalid.
.emc_shear_propose <- function(cur, u, alpha, fits, prior, spec) {
  p <- length(cur$mu)
  u_sd <- u[p + seq_len(p)]
  if (all(u_sd == 0)) {
    new <- .emc_shear_shift_state(cur, u[seq_len(p)], prior)
    if (!is.finite(new$lp)) return(NULL)
    shift <- vapply(seq_along(new$frames), function(s)
      new$frames[[s]]$m - cur$frames[[s]]$m, numeric(p))
    return(list(state = new, alpha = alpha + shift, log_jac = 0))
  }
  d <- exp(u_sd)
  new <- .emc_shear_state(cur$mu + u[seq_len(p)], cur$V * outer(d, d),
                          fits, prior, spec)
  if (is.null(new) || !is.finite(new$lp)) return(NULL)
  alpha_new <- matrix(vapply(seq_len(ncol(alpha)), function(s) {
    f0 <- cur$frames[[s]]; f1 <- new$frames[[s]]
    f1$m + drop(backsolve(f1$R, f0$R %*% (alpha[, s] - f0$m)))
  }, numeric(p)), p)
  list(state = new, alpha = alpha_new,
       log_jac = sum((spec$size + 1) * u_sd) + cur$logdet - new$logdet)
}

# What every round of one move shares: the usable fits and the covariance
# prior.
.emc_shear_context <- function(sampler) {
  list(fits = .emc_shear_fits(sampler), spec = .emc_shear_v_spec(sampler))
}

# A mean step's candidates lie on a line: with Delta_s = Q_s^-1 V^-1 d, the
# subjects move to a_s + t Delta_s and mu to mu + t d, so the hierarchical
# and mean-prior terms are quadratics in t. Returns, per candidate, the moved
# subjects and the change of those terms (V, and so the Jacobian, fixed).
.emc_shear_line <- function(cur, d, alpha, t, prior) {
  Vi_d <- drop(cur$Vi %*% d)
  Delta <- vapply(cur$frames, function(f)
    drop(backsolve(f$R, forwardsolve(t(f$R), Vi_d))), numeric(length(d)))
  Delta <- matrix(Delta, length(d))
  Z0 <- backsolve(cur$chol_V, alpha - cur$mu, transpose = TRUE)
  Zd <- backsolve(cur$chol_V, Delta - d, transpose = TRUE)
  P0 <- as.matrix(prior$theta_mu_invar)
  m0 <- cur$mu - prior$theta_mu_mean
  a1 <- sum(Z0 * Zd) + sum(m0 * (P0 %*% d))
  a2 <- sum(Zd^2) + sum(d * (P0 %*% d))
  list(Delta = Delta, delta = -(a1 * t + 0.5 * a2 * t^2))
}

# One round along step direction v with scale step (see the header); x
# carries the current state, its subject likelihoods and the worker pool.
# Adds E, the selection-weighted squared move in t, and whether it moved.
.emc_shear_round <- function(x, v, step, sampler, ctx) {
  t_aux <- step * stats::rnorm(1L)
  t <- t_aux + step * stats::rnorm(.EMC_SHEAR_CANDIDATES)
  p <- nrow(x$alpha)
  n <- ncol(x$alpha)
  mean_step <- all(v[p + seq_len(p)] == 0)
  if (mean_step) {
    d <- v[seq_len(p)]
    line <- .emc_shear_line(x$cur, d, x$alpha, t, sampler$prior)
    ok <- which(is.finite(line$delta))
    moved_alpha <- function(i) x$alpha + t[i] * line$Delta
    prior_delta <- line$delta
  } else {
    props <- lapply(t, function(tj)
      .emc_shear_propose(x$cur, tj * v, x$alpha, ctx$fits, sampler$prior, ctx$spec))
    ok <- which(!vapply(props, is.null, logical(1)))
    moved_alpha <- function(i) props[[i]]$alpha
    prior_delta <- vapply(seq_along(t), function(i) {
      pr <- props[[i]]
      if (is.null(pr)) return(-Inf)
      .emc_shear_log_hier(pr$alpha, pr$state) + pr$state$lp + pr$log_jac -
        x$cur_hier - x$cur$lp
    }, numeric(1))
  }
  lw <- rep(-Inf, length(t))
  ll_new <- NULL
  if (length(ok)) {
    pn <- sampler$par_names
    cand <- lapply(ok, moved_alpha)
    rows <- lapply(seq_len(n), function(s)
      matrix(vapply(cand, function(a) a[, s], numeric(p)), length(ok), byrow = TRUE,
             dimnames = list(NULL, pn)))
    evaluated <- .emc_wpool_subject_ll(x$wpool, rows, sampler, x$wpool_part,
                                      source = "shear")
    x$wpool <- evaluated$pool
    x$rounds <- x$rounds + 1L
    ll_new <- matrix(evaluated$ll, n, length(ok))
    v_ok <- colSums(ll_new) - sum(x$ll) + prior_delta[ok]
    lw[ok] <- ifelse(is.finite(v_ok), v_ok, -Inf)
  }
  w <- c(0, lw)
  w <- exp(w - max(w))
  w <- w / sum(w)
  j <- sample.int(length(w), 1L, prob = w) - 1L
  x$E <- sum(w[-1L] * t^2)
  x$moved <- j > 0L
  x$t <- if (x$moved) t[j] else 0
  if (x$moved) {
    if (mean_step) {
      x$cur <- .emc_shear_shift_state(x$cur, t[j] * d, sampler$prior)
      x$alpha <- moved_alpha(j)
    } else {
      x$cur <- props[[j]]$state
      x$alpha <- props[[j]]$alpha
    }
    x$cur_hier <- .emc_shear_log_hier(x$alpha, x$cur)
    x$ll <- ll_new[, match(j, ok)]
  }
  x
}

.emc_shear_move <- function(sampler, pars, proposals, stage, wpool = NULL,
                             wpool_part = NULL) {
  keep <- list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool)
  if (stage == "preburn" || !.emc_shear_active(sampler)) return(keep)
  st <- sampler$shear
  plan <- st$plan
  n <- sampler$n_subjects
  p <- sampler$n_pars
  if (!is.list(plan) || !is.matrix(plan$dirs) || any(dim(plan$dirs) != p) ||
      length(plan$rounds) != 2L * p || !any(plan$rounds > 0)) return(keep)
  if (is.null(st$iterations)) st$iterations <- st$rounds <- 0L
  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  # The allocated expected rounds per coordinate: the whole part, plus one
  # more with the fractional part as probability, in random order.
  whole <- floor(plan$rounds)
  n_rounds <- whole + (stats::runif(2L * p) < plan$rounds - whole)
  take <- rep.int(seq_len(2L * p), n_rounds)
  st$iterations <- st$iterations + 1L
  if (!length(take)) {
    sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
    sampler$shear <- st
    keep$sampler <- sampler
    return(keep)
  }
  ctx <- .emc_shear_context(sampler)
  alpha <- matrix(proposals[seq_len(p), ], p)
  cur <- .emc_shear_state(as.numeric(pars$tmu), as.matrix(pars$tvar), ctx$fits,
                          sampler$prior, ctx$spec)
  if (is.null(cur) || !is.finite(cur$lp)) {
    sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
    sampler$shear <- st
    keep$sampler <- sampler
    return(keep)
  }
  x <- list(cur = cur, cur_hier = .emc_shear_log_hier(alpha, cur),
            alpha = alpha, ll = proposals[p + 1L, ], rounds = 0L,
            wpool = wpool, wpool_part = wpool_part)
  adapt <- stage != "sample"
  any_moved <- FALSE
  blk <- .emc_shear_record(plan$block, p)
  net <- numeric(2L * p)
  for (k in take[sample.int(length(take))]) {
    x <- .emc_shear_round(x, .emc_shear_coord_dir(plan, k), exp(plan$log_step[k]),
                          sampler, ctx)
    net[k] <- net[k] + x$t
    blk$E[k] <- blk$E[k] + x$E
    blk$n[k] <- blk$n[k] + 1L
    blk$moved[k] <- blk$moved[k] + as.integer(x$moved)
    if (adapt) {
      plan$log_step[k] <- min(.EMC_SHEAR_LOG_STEP_RANGE[2], max(.EMC_SHEAR_LOG_STEP_RANGE[1],
        plan$log_step[k] + (x$moved - .EMC_SHEAR_MOVE_TARGET) / sqrt(blk$n[k])))
    }
    any_moved <- any_moved || x$moved
  }
  # Each coordinate's net move over its rounds in this iteration, by how many
  # it had (the rounds' saturation, R/move_allocation.R).
  ran <- which(n_rounds > 0L)
  blk$N[ran] <- blk$N[ran] + net[ran]^2
  for (kk in unique(n_rounds[ran]))
    blk$K <- .emc_alloc_khist_add(blk$K, ran[n_rounds[ran] == kk], kk)
  plan$block <- blk
  st$plan <- plan
  st$rounds <- st$rounds + x$rounds
  if (ctx$spec$uses_a_half)
    pars$a_half <- stats::setNames(.emc_shear_draw_a_half(x$cur$V, ctx$spec),
                                   names(pars$a_half))
  sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  if (any_moved) {
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
  sampler$shear <- st
  list(sampler = sampler, pars = pars, proposals = proposals, wpool = x$wpool,
       rounds = x$rounds)
}

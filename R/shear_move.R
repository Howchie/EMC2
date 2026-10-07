# Shear group update, the default group update: the ensemble update
# (R/ensemble_move.R) followed by a shear of the group state.
# options(emc2.group_move = "ensemble") keeps the ensemble update alone.
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
# H_s = 0, so a_s -> mu' + D (a_s - mu). For a given u the map is a
# bijection whose inverse is the same map with -u, and u and -u are equally
# likely, so the Metropolis-Hastings ratio is the posterior ratio times the
# Jacobian: sum_i (g_i + 1) u_sd,i for V (g_i is the size of the covariance
# block holding parameter i; the prior density of V is taken over the free
# entries of its blocks) and sum_s log|R_s| - log|R_s'| for the subjects. The
# approximation shapes the proposal only; a wrong one lowers acceptance but
# not validity.
#
# Where the steps go. Under the approximations, the centred Gibbs step moves
# the group mean as an autoregression with matrix
#   A = (n V^-1 + P0)^-1 V^-1 (sum_s Q_s^-1) V^-1,  Q_s = H_s + V^-1,
# whose eigenvalues lambda_i are the lag-1 autocorrelations of its
# directions; for a group sd, the prior share of each subject's conditional
# precision, w_i = mean_s (Q_s^-1 V^-1)_ii, plays that role. Each iteration
# takes a mean step (u_sd = 0) along direction i with probability lambda_i
# and a scale step (u_mu = 0) on parameter i with probability w_i, in random
# order, so a well-informed hierarchy gets few steps. The directions and
# probabilities are computed at proposal rebuilds from the group covariance
# over the chain's own window and frozen in between, so they never depend on
# the state being moved. Every direction and every scale has its own step
# size (directions in units of their posterior sd under the
# approximations), adapted by Robbins-Monro in the discarded stages only.
#
# When every subject has a usable approximation, each step is delayed
# acceptance (Christen & Fox 2005): a first stage on the approximations,
# which costs no likelihood evaluation, then the exact correction
# exp(sum_s r_s(a_s') - r_s(a_s)), r = ll - log L~; otherwise it is one
# exact stage.
#
# Under the Huang-Wand prior (standard type), V | a ~ IW(v + d - 1,
# 2 v diag(1 / a)) with a_k ~ IG(1/2, 1 / A_k^2). Holding a fixed would pin
# a scale step to the prior's conditional width in log sd whatever the data
# say. Given V the inverse-Wishart density factorises over the a_k, so a
# integrates out in closed form:
#   log p(V) = -(df + d + 1)/2 log|V| - (v + d)/2 sum_k log(v (V^-1)_kk + 1/A_k^2)
# (for an unblocked parameter, the d = 1 case on its variance). The steps
# target the posterior with a integrated out, and a is then drawn from its
# conditional given the final V, the draw the Gibbs step makes: together a
# valid update of (V, a).
#
# Supported hierarchies: those in which every subject's prior mean is the
# group mean (standard without group designs, diagonal-gamma); the others
# keep the ensemble update alone.

.emc_shear_options <- function() {
  list(enabled = identical(getOption("emc2.group_move", "shear"), "shear"),
       target = 0.44)
}

.emc_shear_supported <- function(sampler) {
  isTRUE(sampler$type %in% c("standard", "diagonal-gamma")) &&
    !any(sampler$nuisance) &&
    !any(as.logical(sampler$marginalised_idx), na.rm = TRUE) &&
    is.null(sampler$pooled) && is.null(sampler$gd) &&
    length(sampler$group_designs) == 0L &&
    isTRUE(sampler$n_subjects >= 2L)
}

.emc_shear_active <- function(sampler) {
  .emc_shear_options()$enabled && .emc_shear_supported(sampler)
}

# Each subject's likelihood approximation if it is usable (checked, and the
# check not poor), else NULL. The particle step uses the same ones.
.emc_shear_fits <- function(sampler) {
  p <- sampler$n_pars
  lapply(seq_len(sampler$n_subjects), function(s) {
    fit <- if (length(sampler$lik_laplace) >= s) sampler$lik_laplace[[s]] else NULL
    ok <- is.list(fit) && is.matrix(fit$prec) && all(dim(fit$prec) == p) &&
      length(fit$lin) == p && all(is.finite(fit$prec)) && all(is.finite(fit$lin)) &&
      is.list(fit$check) && !isTRUE(any(.emc_laplace_poor(fit$check)))
    if (ok) fit else NULL
  })
}

# Mean directions and step probabilities at group covariance V (see the
# header). Columns of dirs move mu by one posterior sd of their direction
# under the approximations: mu = mu0 + t d with d' (n V^-1 + P0) d = 1 has
# collapsed precision 1 - lambda in t.
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

# Plan from the mean group covariance over the chain's proposal window.
.emc_shear_set_plan <- function(chain, var_draws) {
  if (!.emc_shear_active(chain)) return(chain)
  if (length(dim(var_draws)) != 3L || dim(var_draws)[3] < 1L) return(chain)
  V <- apply(var_draws, 1:2, mean)
  plan <- tryCatch(.emc_shear_plan(.emc_shear_fits(chain), chain$prior, V),
                   error = function(e) NULL)
  if (is.null(plan)) return(chain)
  st <- chain$shear
  if (is.null(st)) st <- list()
  st$plan <- plan
  chain$shear <- st
  chain
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
       logdet = sum(vapply(frames, `[[`, 0, "logdet")),
       lp = .emc_shear_log_prior_mu(mu, prior) + lp_v)
}

# sum_s log N(a_s | mu, V).
.emc_shear_log_hier <- function(alpha, state) {
  z <- backsolve(state$chol_V, alpha - state$mu, transpose = TRUE)
  -0.5 * sum(z^2) - ncol(alpha) * sum(log(diag(state$chol_V)))
}

# The deterministic map for a given u, from the group state cur and subjects
# alpha: the moved state, the moved subjects and log|det| of the map's
# Jacobian. NULL when the moved state is invalid.
.emc_shear_propose <- function(cur, u, alpha, fits, prior, spec) {
  p <- length(cur$mu)
  u_sd <- u[p + seq_len(p)]
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

# Gaussian log-likelihood difference; subjects without a fit contribute none.
.emc_shear_lik_delta <- function(alpha, proposed, fits) {
  delta <- 0
  for (s in seq_len(ncol(alpha))) {
    fit <- fits[[s]]
    if (is.null(fit)) next
    a <- alpha[, s]
    b <- proposed[, s]
    delta <- delta + sum((b - a) * (fit$lin - drop(fit$prec %*% (a + b)) / 2))
  }
  if (is.finite(delta)) delta else NA_real_
}

# What every step of one move shares: the usable fits, the covariance prior,
# and whether to screen on the fits first (only when every subject has one).
.emc_shear_context <- function(sampler) {
  fits <- .emc_shear_fits(sampler)
  list(fits = fits, spec = .emc_shear_v_spec(sampler),
       delayed = !any(vapply(fits, is.null, logical(1))))
}

# One Metropolis-Hastings step of the map with step u; x carries the current
# state, its subject likelihoods and the worker pool.
.emc_shear_mh <- function(x, u, sampler, ctx) {
  x$accept_prob <- 0
  x$accepted <- FALSE
  prop <- .emc_shear_propose(x$cur, u, x$alpha, ctx$fits, sampler$prior, ctx$spec)
  if (is.null(prop)) return(x)
  pn <- sampler$par_names
  new_hier <- .emc_shear_log_hier(prop$alpha, prop$state)
  common <- new_hier - x$cur_hier + prop$state$lp - x$cur$lp + prop$log_jac
  approx_delta <- if (ctx$delayed) .emc_shear_lik_delta(x$alpha, prop$alpha, ctx$fits) else NA_real_
  delayed <- is.finite(common + approx_delta)
  probability <- function(log_ratio) if (is.na(log_ratio)) 0 else exp(min(0, log_ratio))
  screen <- if (delayed) probability(common + approx_delta) else 1
  if (delayed && stats::runif(1L) >= screen) return(x)
  rows <- lapply(seq_len(ncol(x$alpha)), function(s)
    matrix(prop$alpha[, s], 1L, dimnames = list(NULL, pn)))
  evaluated <- .emc_wpool_subject_ll(x$wpool, rows, sampler, x$wpool_part)
  ll_new <- evaluated$ll
  x$wpool <- evaluated$pool
  x$rounds <- x$rounds + 1L
  correction <- sum(ll_new - x$ll) + if (delayed) -approx_delta else common
  x$accept_prob <- screen * probability(correction)
  if (stats::runif(1L) < probability(correction)) {
    x$cur <- prop$state; x$cur_hier <- new_hier
    x$alpha <- prop$alpha; x$ll <- ll_new; x$accepted <- TRUE
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
  if (!is.list(plan) || !is.matrix(plan$dirs) || any(dim(plan$dirs) != p)) return(keep)
  cfg <- .emc_shear_options()
  # Mean steps start at one posterior sd of their direction; scale steps near
  # the posterior sd of a log group sd informed by n subjects, 1 / sqrt(2 n).
  if (length(st$dir_log_step) != p) {
    st$dir_log_step <- numeric(p)
    st$dir_n_adapt <- st$dir_proposed <- st$dir_accepted <- integer(p)
  }
  if (length(st$scale_log_step) != p) {
    st$scale_log_step <- rep(-0.5 * log(2 * n), p)
    st$scale_n_adapt <- st$scale_proposed <- st$scale_accepted <- integer(p)
  }
  if (is.null(st$iterations)) st$iterations <- st$rounds <- 0L
  ctx <- .emc_shear_context(sampler)
  alpha <- matrix(proposals[seq_len(p), ], p)
  cur <- .emc_shear_state(as.numeric(pars$tmu), as.matrix(pars$tvar), ctx$fits,
                          sampler$prior, ctx$spec)
  if (is.null(cur) || !is.finite(cur$lp)) return(keep)
  x <- list(cur = cur, cur_hier = .emc_shear_log_hier(alpha, cur),
            alpha = alpha, ll = proposals[p + 1L, ], rounds = 0L,
            wpool = wpool, wpool_part = wpool_part)

  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  adapt <- stage != "sample"
  zero <- numeric(p)
  any_accepted <- FALSE
  # Candidate steps: p mean directions, then p scale parameters, each taken
  # with its frozen probability, in random order. Each is a complete
  # Metropolis-Hastings step from the current state.
  take <- which(stats::runif(2L * p) < c(plan$lam, plan$w))
  for (k in take[sample.int(length(take))]) {
    if (k <= p) {
      u <- c(exp(st$dir_log_step[k]) * stats::rnorm(1L) * plan$dirs[, k], zero)
      x <- .emc_shear_mh(x, u, sampler, ctx)
      st$dir_proposed[k] <- st$dir_proposed[k] + 1L
      st$dir_accepted[k] <- st$dir_accepted[k] + as.integer(x$accepted)
      if (adapt) {
        st$dir_n_adapt[k] <- st$dir_n_adapt[k] + 1L
        st$dir_log_step[k] <- st$dir_log_step[k] +
          (x$accepted - cfg$target) / sqrt(st$dir_n_adapt[k])
      }
    } else {
      i <- k - p
      u_sd <- zero
      u_sd[i] <- exp(st$scale_log_step[i]) * stats::rnorm(1L)
      x <- .emc_shear_mh(x, c(zero, u_sd), sampler, ctx)
      st$scale_proposed[i] <- st$scale_proposed[i] + 1L
      st$scale_accepted[i] <- st$scale_accepted[i] + as.integer(x$accepted)
      if (adapt) {
        st$scale_n_adapt[i] <- st$scale_n_adapt[i] + 1L
        st$scale_log_step[i] <- st$scale_log_step[i] +
          (x$accepted - cfg$target) / sqrt(st$scale_n_adapt[i])
      }
    }
    any_accepted <- any_accepted || x$accepted
  }
  st$iterations <- st$iterations + 1L
  st$rounds <- st$rounds + x$rounds
  if (ctx$spec$uses_a_half)
    pars$a_half <- stats::setNames(.emc_shear_draw_a_half(x$cur$V, ctx$spec),
                                   names(pars$a_half))
  sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  if (any_accepted) {
    pars$tmu <- stats::setNames(x$cur$mu, names(pars$tmu))
    pars$tvar <- `dimnames<-`(x$cur$V, dimnames(pars$tvar))
    pars$tvinv <- `dimnames<-`(x$cur$Vi, dimnames(pars$tvar))
    proposals[seq_len(p), ] <- x$alpha
    proposals[p + 1L, ] <- x$ll
    pars$alpha <- x$alpha
  }
  sampler$shear <- st
  list(sampler = sampler, pars = pars, proposals = proposals, wpool = x$wpool)
}

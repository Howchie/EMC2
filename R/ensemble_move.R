# Ensemble group update (Neal 2011, arXiv:1101.0387; Shestopaloff & Neal 2013,
# arXiv:1305.0320, section 5), specialised to conditionally independent
# subjects.
#
# Target: pi(h, a | y) propto p(h) prod_s L_s(a_s) g(a_s | h), where h is the
# group state and g the subject prior. For every subject build a pool A_s of M
# states: the current a_s at a uniformly random position k_s, the other M - 1
# drawn independently from a frozen proposal q_s that does not depend on h. The
# augmented target
#   pi~(h, A, k) propto p(h) prod_s [ L_s(a_{s,k_s}) g(a_{s,k_s} | h) / M
#                                      prod_{j != k_s} q_s(a_sj) ]
# has pi(h, a | y) as its marginal. Conditional on the pools its full
# conditionals are
#   k_s | h, A  propto  L_s(a_sj) g(a_sj | h) / q_s(a_sj)          (reselection)
#   h   | k, A  propto  p(h) prod_s g(a_{s,k_s} | h)               (group Gibbs)
# and the second is exactly the conjugate group Gibbs step with the
# selected subject states. Alternating the two for several rounds is a Gibbs
# sampler on pi~(h, k | A): the group moves as if each subject could be any of
# its pool members, which is what lets it leave the centred Gibbs step's
# sigma / sqrt(n) neighbourhood. Only the pool costs likelihood evaluations
# (M - 1 per subject, one worker round trip); the rounds are hierarchy
# arithmetic. The final reselection under the final h sets the new subject
# states. Pools are discarded at the end of the sweep.
#
# No quantity is adapted inside the move. q_s is a two-component mixture built
# from the subject's posterior moments over the last warm-up window
# (sampler$ensemble_q, set with the chain proposals): a normal N(mu, V) and a
# defensive Student-t (Hesterberg 1995) with df 4 and SCALE matrix c V,
# c = 1 + 3 / sqrt(p). The t's covariance is 2 c V, not c V. The reselection
# is an independence-type step, so q_s lighter than the target in some
# direction makes rare draws in its tail win the reselection and then linger
# (Mengersen & Tweedie 1996); the t's polynomial tails guard against a
# covariance estimate that is too narrow. That is a robustness argument, not
# a containment proof. The radial spread of the t component does not
# concentrate: its squared Mahalanobis distance M2 (relative to V) satisfies
# M2 / p -> 4 c / chi^2_4 as p grows, with mean 2 c and infinite variance.
# Shrinking c with p is therefore a heuristic for how much of the pool the
# tail spends far from the bulk, chosen empirically; at p = 1, c = 4. q_s is
# rebuilt with the chain proposals: at every refresh in the discarded stages,
# and in sampling on the increasingly-rare schedule of add_proposals() (after
# 100, 400, 1600, ... production transitions).

# Pool size M, number of Gibbs rounds per sweep, and the defensive t
# component of q_s (mixture weight, degrees of freedom, scale excess).
.emc_ensemble_options <- function() {
  list(pool = 32L, rounds = 10L, tail_weight = 0.25, tail_df = 4, tail_excess = 3)
}

# Supported hierarchies: those whose group Gibbs step and subject prior are
# already used by the particle update. Nuisance coordinates have their own
# samplers and are not covered. A marginalised coordinate is stored as a
# draw from its conditional posterior, so the subject states the move sees
# are complete.
.emc_ensemble_supported <- function(sampler) {
  isTRUE(sampler$type %in% c("standard", "diagonal-gamma", "factor",
                             "infnt_factor", "SEM")) &&
    !any(sampler$nuisance) &&
    isTRUE(sampler$n_subjects >= 2L)
}

# One-iteration copy of the sample store holding the current group state, so
# gibbs_step() can read it through last_sample_<type>() and fill_samples() can
# write it without touching the chain's history (fill_samples writes in place).
.emc_ensemble_store <- function(samples) {
  n_iter <- length(samples$stage)
  idx <- samples$idx
  out <- samples
  hist <- setdiff(names(samples), .EMC_STATIC_SAMPLE_FIELDS)
  out[hist] <- lapply(samples[hist], function(x) {
    d <- dim(x)
    if (is.null(d) || d[length(d)] != n_iter) return(x)
    if (length(d) == 2L) x[, idx, drop = FALSE] else
      if (length(d) == 3L) x[, , idx, drop = FALSE] else x
  })
  out$idx <- 1L
  out$stage <- samples$stage[idx]
  out
}

# Pool proposal q_s from one Cholesky factor per subject:
# (1 - w) N(mu, V) + w t_df(mu, c V), c = 1 + tail_excess / sqrt(p), with
# Mahalanobis distances reused for both components.
.emc_ensemble_q <- function(mu, var, cfg) {
  L <- t(chol(var))
  p <- length(mu)
  list(mu = mu, L = L, logdet = 2 * sum(log(diag(L))), p = p,
       c = 1 + cfg$tail_excess / sqrt(p))
}

.emc_ensemble_log_q <- function(x, q, cfg) {
  z <- forwardsolve(q$L, t(x) - q$mu)
  .emc_ensemble_log_q_m2(colSums(z^2), q, cfg)
}

# The same density from squared Mahalanobis distances m2.
.emc_ensemble_log_q_m2 <- function(m2, q, cfg) {
  p <- q$p; nu <- cfg$tail_df; c <- q$c
  a <- log1p(-cfg$tail_weight) - 0.5 * (p * log(2 * pi) + q$logdet + m2)
  b <- log(cfg$tail_weight) + lgamma((nu + p) / 2) - lgamma(nu / 2) -
    0.5 * (p * log(nu * pi) + q$logdet + p * log(c)) -
    (nu + p) / 2 * log1p(m2 / (c * nu))
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

.emc_ensemble_draw_q <- function(n, q, cfg) {
  z <- matrix(stats::rnorm(n * q$p), q$p, n)
  tail <- stats::runif(n) < cfg$tail_weight
  scale <- rep(1, n)
  if (any(tail)) scale[tail] <- sqrt(q$c * cfg$tail_df /
                                       stats::rchisq(sum(tail), cfg$tail_df))
  x <- t(q$mu + (q$L %*% z) * rep(scale, each = q$p))
  attr(x, "m2") <- scale^2 * colSums(z^2)
  x
}

# Subject prior of every subject under one group state h, as a p x n matrix of
# means and one covariance shared by all subjects. For the latent-factor
# hierarchies h contains the factor scores eta, and the prior is conditional
# on them: a_s | h ~ N(mu + K x_s + Lambda eta_s, Sigma_err), the density the
# group Gibbs step conditions on. `weight = "marginal"` (SEM only, used by
# the burn mode start) integrates eta_s out instead:
#   a_s | h ~ N(mu + (K + Lambda B0^-1 G) x_s, Lambda B0^-1 Delta^-1 B0^-T Lambda' + Sigma_err).
# Factor and infinite-factor models have no covariates (K = G = 0, B = 0).
.emc_ensemble_subject_prior <- function(pars, sampler, n, weight = "conditional") {
  type <- sampler$type
  rep_mu <- function(m) matrix(as.numeric(m), length(m), n)
  if (type == "standard" && is.matrix(pars$subj_mu))
    return(list(mean = pars$subj_mu, var = pars$tvar))    # = get_group_level_standard
  if (type == "diagonal-gamma")                            # = get_group_level_base
    return(list(mean = rep_mu(pars$tmu), var = pars$tvar, dvar = diag(pars$tvar)))
  if (type %in% c("standard", "diagonal-gamma")) {
    gl <- lapply(seq_len(n), function(s) get_group_level(pars, s, type))
    return(list(mean = vapply(gl, function(g) as.numeric(g$mu), numeric(length(gl[[1L]]$mu))),
                var = gl[[1L]]$var))
  }
  mu <- as.numeric(pars$tmu)
  p <- length(mu)
  if (type == "SEM") {
    x <- sampler$sem_settings$covariates
    x <- if (is.null(x)) matrix(0, n, 0) else as.matrix(x)
    lam <- matrix(pars$lambda, p)
    nf <- ncol(lam)
    K <- matrix(pars$K, p, ncol(x))
    eps_var <- 1 / as.numeric(pars$epsilon_inv)
    if (weight == "conditional")
      return(list(mean = mu + K %*% t(x) + lam %*% t(matrix(pars$eta, n, nf)),
                  var = diag(eps_var, p), dvar = eps_var))
    B0_inv <- solve(diag(nf) - matrix(pars$B, nf, nf))
    G <- matrix(pars$G, nf, ncol(x))
    A <- lam %*% B0_inv
    return(list(mean = mu + (K + A %*% G) %*% t(x),
                var = A %*% solve(matrix(pars$delta_inv, nf, nf)) %*% t(A) + diag(eps_var, p)))
  }
  if (type == "factor") {
    lam <- matrix(pars$lambda_untransf, p)
    eps_var <- 1 / as.numeric(pars$sig_err_inv)
  } else {                                        # infnt_factor
    lam <- matrix(pars$lambda, p)
    eps_var <- 1 / as.numeric(pars$epsilon_inv)
  }
  list(mean = mu + lam %*% t(matrix(pars$eta, n, ncol(lam))), var = diag(eps_var, p),
       dvar = eps_var)
}

# Factor scores drawn from their full conditional given the subject states.
.emc_ensemble_refresh_scores <- function(pars, sampler, alpha) {
  type <- sampler$type
  if (!type %in% c("factor", "infnt_factor", "SEM")) return(pars)
  p <- sampler$n_pars
  n <- sampler$n_subjects
  loadings <- matrix(if (type == "factor") pars$lambda_untransf else pars$lambda, p)
  nf <- ncol(loadings)
  if (nf == 0L) return(pars)
  residual <- sweep(t(alpha), 2, as.numeric(pars$tmu))
  if (type == "SEM") {
    x <- sampler$sem_settings$covariates
    if (is.null(x)) x <- matrix(0, n, 0)
    x <- as.matrix(x)
    scores <- .sem_draw_eta(residual, x, matrix(pars$K, p, ncol(x)),
                           matrix(pars$G, nf, ncol(x)), matrix(pars$B, nf, nf),
                           matrix(pars$delta_inv, nf, nf), matrix(pars$lambda, p, nf),
                           diag(as.numeric(pars$epsilon_inv), p), nf, n)
  } else {
    error_precision <- if (type == "factor") pars$sig_err_inv else pars$epsilon_inv
    score_precision <- if (type == "factor") pars$psi_inv else rep(1, nf)
    weighted <- sweep(loadings, 1, error_precision, "*")
    covariance <- solve(diag(as.numeric(score_precision), nf) + crossprod(loadings, weighted))
    scores <- residual %*% weighted %*% covariance + .emc_rmvnorm(n, sigma = covariance)
  }
  pars$eta[,] <- scores
  pars
}

# Subject prior log densities of every pool member under one group state.
.emc_ensemble_log_g <- function(pool_stack, subj_of_row, pars, sampler, n) {
  sp <- .emc_ensemble_subject_prior(pars, sampler, n)
  if (!is.null(sp$dvar)) {                     # diagonal covariance: no factorisation
    z <- (t(pool_stack) - sp$mean[, subj_of_row, drop = FALSE]) / sqrt(sp$dvar)
    return(-0.5 * (length(sp$dvar) * log(2 * pi) + sum(log(sp$dvar)) + colSums(z^2)))
  }
  L <- t(chol(sp$var))
  z <- forwardsolve(L, t(pool_stack) - sp$mean[, subj_of_row, drop = FALSE])
  -0.5 * (nrow(L) * log(2 * pi) + 2 * sum(log(diag(L))) + colSums(z^2))
}

# One categorical draw per column of log-weights, all columns at once: the
# Gumbel-max trick, argmax_j (lw_j + G_j) with G_j iid standard Gumbel, is an
# exact draw with probabilities proportional to exp(lw_j).
.emc_ensemble_select_all <- function(lw) {
  g <- -log(-log(stats::runif(length(lw))))
  max.col(t(lw + g), ties.method = "first")
}

# Fresh pool members for every subject. They are drawn from the frozen q_s
# alone, never from the current state, so they can be drawn before the subject
# step and their likelihoods computed in the same worker round trip as the
# particles (run_stage passes them through .emc_wpool_iter()). Returns NULL
# when the update does not run (preburn, unsupported hierarchy, no proposals). The draws and their proposal
# densities come from one C++ call (src/ensemble_move.cpp). `factors` (the
# proposals' Cholesky factors and log determinants, from
# emc_ensemble_factors_cpp) can be passed in, since the proposals are frozen
# within a block; without them they are computed here.
.emc_ensemble_draw_pools <- function(sampler, stage, chains_mu, chains_var,
                                     factors = NULL) {
  cfg <- .emc_ensemble_options()
  if (stage == "preburn" || !.emc_ensemble_supported(sampler)) return(NULL)
  n <- sampler$n_subjects
  if (length(chains_mu) != n || length(chains_var) != n ||
      any(vapply(chains_mu, is.null, logical(1))) ||
      any(vapply(chains_var, is.null, logical(1))))
    return(NULL)
  if (is.null(factors)) factors <- emc_ensemble_factors_cpp(lapply(chains_var, as.matrix))
  d <- emc_ensemble_draw_pools_cpp(lapply(chains_mu, as.numeric), factors$L, factors$logdet,
                                   cfg$pool - 1L, cfg$tail_weight, cfg$tail_df,
                                   cfg$tail_excess, sampler$par_names)
  list(rows = d$rows, lq = d$lq, ll = NULL, L = factors$L, mu = d$mu,
       logdet = factors$logdet, c = d$c)
}

# Log-likelihoods of a matrix of candidate states for one subject. A failure
# is a zero-density candidate: it is counted under the failure policy and
# returned as -Inf (strict mode stops instead).
.emc_ensemble_ll_checked <- function(props, data, model,
                                     ll_cache_validated = FALSE) {
  ll <- tryCatch(as.numeric(calc_ll_manager(props, dadm = data, model = model,
                                             r_cores = 1L,
                                             ll_cache_validated = ll_cache_validated)),
    error = function(e) {
      cls <- .emc_classify_failure(e)
      stop(errorCondition(paste("Subject likelihood evaluation failed:",
        conditionMessage(e)), class = c("emc_ensemble_failure", "emc_failure"),
        emc_class = cls, emc_source = "ensemble", emc_cause = e))
    })
  if (length(ll) != nrow(props) || anyNA(ll) || any(ll == Inf)) {
    cls <- if (length(ll) != nrow(props)) "programming" else "numerical"
    stop(errorCondition(
      "Invalid subject likelihood: expected finite values or -Inf",
      class = c("emc_ensemble_failure", "emc_failure"),
      emc_class = cls, emc_source = "ensemble"))
  }
  ll
}

.emc_ensemble_ll_candidate <- function(props, data, model,
                                      ll_cache_validated = FALSE) {
  tryCatch(.emc_ensemble_ll_checked(props, data, model,
                                    ll_cache_validated = ll_cache_validated),
    error = function(e) {
      cls <- if (inherits(e, "emc_failure") && !is.null(e$emc_class))
        e$emc_class else .emc_classify_failure(e)
      .emc_reject_record(e, "ensemble")
      if (identical(.emc_failure_action(cls), "abort"))
        .emc_failure_abort(cls, "ensemble", e)
      rep(-Inf, nrow(props))
    })
}

# Pool log-likelihoods. The subject step normally returns them with its
# particles; any subject it did not cover is evaluated here.
.emc_ensemble_pool_ll <- function(sampler, pools,
                                 ll_cache_validated = FALSE) {
  n <- sampler$n_subjects
  M1 <- nrow(pools$rows[[1L]])
  ll <- pools$ll
  if (is.null(ll)) ll <- vector("list", n)
  for (s in seq_len(n)) if (is.null(ll[[s]]) || length(ll[[s]]) != M1) {
    ll[[s]] <- tryCatch(
      .emc_ensemble_ll_checked(pools$rows[[s]], sampler$data[[s]], sampler$model,
                               ll_cache_validated = ll_cache_validated),
      error = function(e) vapply(seq_len(M1), function(i)
        .emc_ensemble_ll_candidate(pools$rows[[s]][i, , drop = FALSE],
                                   sampler$data[[s]], sampler$model,
                                   ll_cache_validated = ll_cache_validated), numeric(1)))
  }
  ll_new <- matrix(unlist(ll), nrow = M1)
  if (anyNA(ll_new) || any(ll_new == Inf))
    stop("Invalid ensemble pool likelihood")
  ll_new
}

# The hierarchies whose Gibbs rounds run in C++: standard without group
# designs, and diagonal-gamma. The others take the R rounds below.
.emc_ensemble_rounds_in_cpp <- function(sampler) {
  if (sampler$type == "diagonal-gamma") return(TRUE)
  if (sampler$type != "standard") return(FALSE)
  is.null(sampler$gd) && (is.null(sampler$group_designs) || length(sampler$group_designs) == 0L)
}

# The Gibbs rounds on pi~(h, k | A) in C++ (src/ensemble_move.cpp): the
# group state from the package's conjugate update given the selected subject
# states, then every subject's pool position given the group state. Returns
# the final selection, the group state in the form gibbs_step() returns it,
# and the pools' weight ESS under the final state.
.emc_ensemble_rounds_cpp <- function(sampler, pars, pool_stack, base_w, k0, M, cfg) {
  n <- length(k0)
  p <- ncol(pool_stack)
  if (sampler$type == "standard") {
    is_blocked <- sampler$is_blocked[!sampler$nuisance]
    par_group <- sampler$par_group[!sampler$nuisance]
    if (is.null(is_blocked) || length(is_blocked) == 0L) is_blocked <- rep(TRUE, p)
    if (is.null(par_group)) par_group <- rep(1L, p)
    par_group <- as.integer(factor(par_group, levels = unique(par_group)))
    res <- emc_ensemble_rounds_cpp(pool_stack, base_w, k0, M, cfg$rounds, 0L,
                                   list(tvinv = pars$tvinv, a_half = as.numeric(pars$a_half)),
                                   sampler$prior, as.logical(is_blocked), par_group)
    pars <- list(tmu = res$tmu, tvar = res$tvar, tvinv = res$tvinv, a_half = res$a_half,
                 subj_mu = matrix(res$tmu, p, n), alpha = NULL)
  } else {
    res <- emc_ensemble_rounds_cpp(pool_stack, base_w, k0, M, cfg$rounds, 1L,
                                   list(tvinv = diag(pars$tvinv)), sampler$prior,
                                   logical(0), integer(0))
    tmu <- res$tmu
    names(tmu) <- sampler$par_names[!sampler$nuisance]
    pars <- list(tmu = tmu, tvar = diag(1 / res$tvinv, p), tvinv = diag(res$tvinv, p),
                 alpha = NULL)
  }
  list(selected = res$selected, pars = pars, wess = res$wess)
}

# The same rounds in R, through gibbs_step() and the subject prior of every
# pool member: the path for the latent-factor hierarchies and for standard
# hierarchies with group designs, and the reference the C++ rounds are
# tested against. The scratch store only needs the group state: gibbs_step()
# takes the subject states as an argument.
.emc_ensemble_rounds_r <- function(sampler, pars, pool_stack, base_w, k0, M, cfg) {
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  out <- if (!is.null(seed) && isTRUE(getOption("emc2.native_gibbs_factor", TRUE)))
    .emc_ensemble_rounds_latent(sampler, pars, pool_stack, base_w, k0, M, cfg)
  if (!is.null(out)) return(out)
  if (!is.null(seed)) assign(".Random.seed", seed, envir = globalenv())
  .emc_ensemble_rounds_r_form(sampler, pars, pool_stack, base_w, k0, M, cfg)
}

# The rounds of a factor or infinite-factor hierarchy in one native call
# (src/gibbs_factor.cpp), the same rounds on the same random stream; NULL
# for other hierarchies and where the native form declines.
.emc_ensemble_rounds_latent <- function(sampler, pars, pool_stack, base_w, k0, M, cfg) {
  type <- sampler$type
  n <- length(k0)
  if (!type %in% c("factor", "infnt_factor") || any(sampler$nuisance) ||
      n != sampler$n_subjects || ncol(pool_stack) != sampler$n_pars || !is.double(pool_stack))
    return(NULL)
  # The state as the rounds' first draw reads it back from the store.
  if (type == "factor") {
    inp <- .emc_gibbs_factor_inputs(sampler, list(
      mu = pars$tmu, eta = pars$eta, lambda = pars$lambda_untransf,
      psi_inv = pars$psi_inv, sig_err_inv = pars$sig_err_inv))
  } else {
    inp <- .emc_gibbs_infnt_factor_inputs(sampler, list(
      mu = pars$tmu, delta = pars$delta, lambda = pars$lambda, psi = pars$psi,
      epsilon_inv = pars$epsilon_inv))
    if (!is.null(inp)) inp$state$eta <- NULL
  }
  if (is.null(inp) || any(lengths(inp$state) == 0L)) return(NULL)
  res <- emc_ensemble_rounds_latent_cpp(pool_stack, as.numeric(base_w), as.integer(k0),
                                        as.integer(M), as.integer(cfg$rounds), inp$state,
                                        inp$free, inp$prior)
  if (is.null(res)) return(NULL)
  offset <- (seq_len(n) - 1L) * M
  alpha_k <- t(pool_stack[offset + res$conditioned, , drop = FALSE])
  pars <- if (type == "factor") .emc_gibbs_factor_pars(res$state, t(alpha_k)) else
    .emc_gibbs_infnt_factor_pars(res$state, alpha_k, t(alpha_k),
                                 last_sample_infnt_factor(
                                   .emc_ensemble_store(sampler$samples))$delta)
  list(selected = res$selected, pars = pars, wess = res$wess)
}

.emc_ensemble_rounds_r_form <- function(sampler, pars, pool_stack, base_w, k0, M, cfg) {
  n <- length(k0)
  n_pars <- ncol(pool_stack)
  subj_of_row <- rep(seq_len(n), each = M)
  offset <- (seq_len(n) - 1L) * M
  tmp <- sampler
  store <- .emc_ensemble_store(sampler$samples)
  selected <- k0
  wess <- numeric(n)
  for (r in seq_len(cfg$rounds)) {
    alpha_k <- t(pool_stack[offset + selected, , drop = FALSE])
    store <- fill_samples(samples = store, group_level = pars, proposals = NULL,
                          j = 1L, n_pars = n_pars, type = sampler$type)
    tmp$samples <- store
    pars <- gibbs_step(tmp, alpha_k, sampler$type)
    lw_all <- matrix(base_w + .emc_ensemble_log_g(pool_stack, subj_of_row, pars, sampler, n),
                     nrow = M)
    selected <- .emc_ensemble_select_all(lw_all)
    if (r == cfg$rounds) {
      w <- exp(sweep(lw_all, 2, apply(lw_all, 2, max)))
      wess <- colSums(w)^2 / colSums(w^2)
    }
  }
  list(selected = selected, pars = pars, wess = wess)
}

.emc_ensemble_move <- function(sampler, pars, proposals, stage, wpool, wpool_ctx,
                               wpool_part, pools,
                               ll_cache_validated = FALSE) {
  keep <- list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool)
  if (is.null(pools)) return(keep)
  out <- tryCatch(.emc_ensemble_move_impl(sampler, pars, proposals, stage,
                                          wpool, wpool_ctx, wpool_part, pools,
                                          ll_cache_validated = ll_cache_validated),
                  error = identity)
  if (!inherits(out, "error")) return(out)
  cls <- if (inherits(out, "emc_failure") && length(out$emc_class))
    out$emc_class else .emc_classify_failure(out)
  already_recorded <- inherits(out, "emc_failure") &&
    identical(out$emc_source, "ensemble")
  if (!already_recorded) .emc_reject_record(out, "ensemble")
  if (identical(.emc_failure_action(cls), "abort"))
    .emc_failure_abort(cls, "ensemble", out)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    keep$sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  keep
}

.emc_ensemble_move_impl <- function(sampler, pars, proposals, stage, wpool, wpool_ctx,
                                    wpool_part, pools,
                                    ll_cache_validated = FALSE) {
  keep <- list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool)
  if (is.null(pools)) return(keep)
  cfg <- .emc_ensemble_options()
  n <- sampler$n_subjects
  n_pars <- sampler$n_pars
  M <- cfg$pool
  pn <- sampler$par_names
  ll_new <- .emc_ensemble_pool_ll(sampler, pools,
                                  ll_cache_validated = ll_cache_validated)

  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  # Pools: the current state at a uniform position among the M - 1 fresh
  # draws. Stacked pools, subject s in rows (s - 1) * M + 1:M.
  k0 <- integer(n)
  if (n > 0L) k0[seq_len(n)] <- sample.int(M, n, replace = TRUE)
  st <- emc_ensemble_stack_cpp(pools$rows, pools$lq, ll_new,
                               proposals[seq_len(n_pars + 1L), , drop = FALSE], k0,
                               pools$L, pools$mu, pools$logdet, pools$c,
                               cfg$tail_weight, cfg$tail_df)
  pool_stack <- st$pool_stack
  colnames(pool_stack) <- pn
  offset <- (seq_len(n) - 1L) * M

  # Gibbs on pi~(h, k | A): h | k with the package's group step, then k | h.
  rounds <- if (.emc_ensemble_rounds_in_cpp(sampler)) .emc_ensemble_rounds_cpp else
    .emc_ensemble_rounds_r
  res <- rounds(sampler, pars, pool_stack, st$base_w, k0, M, cfg)
  selected <- res$selected
  pars <- res$pars
  wess <- res$wess

  moved <- selected != k0
  proposals[seq_len(n_pars), ] <- t(pool_stack[offset + selected, , drop = FALSE])
  proposals[n_pars + 1L, ] <- st$ll_stack[offset + selected]
  pars$alpha <- proposals[seq_len(n_pars), , drop = FALSE]

  st <- sampler$ensemble_stats
  if (is.null(st) || is.null(st$subj_moved))
    st <- list(sweeps = 0L, moved = 0, wess = 0, subj_moved = numeric(n), subj_wess = numeric(n))
  st$sweeps <- st$sweeps + 1L
  st$moved <- st$moved + mean(moved)
  st$wess <- st$wess + mean(wess) / M
  st$subj_moved <- st$subj_moved + moved
  st$subj_wess <- st$subj_wess + wess / M
  sampler$ensemble_stats <- st
  sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool)
}

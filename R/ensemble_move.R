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
# and the second is exactly the existing conjugate group Gibbs step with the
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
# (sampler$ensemble_q, set with the chain proposals): a normal and a defensive
# Student-t (Hesterberg 1995) with 4x the covariance. The reselection is an
# independence-type step, so q_s must not be lighter than the target in any
# direction (Mengersen & Tweedie 1996): a too-narrow q_s makes rare draws in
# its tail win the reselection and then linger. The window covariance of a
# slow-mixing subject can be 2x too narrow in SD, hence the 4x. It is rebuilt
# with the chain proposals: at every refresh in the discarded stages, and in
# sampling on the diminishing-adaptation schedule (100, 200, 400, ... draws).

.emc_ensemble_options <- function() {
  int_opt <- function(name, default) {
    x <- suppressWarnings(as.integer(getOption(name, default)))
    if (length(x) != 1L || is.na(x) || x < 1L) default else x
  }
  # Latent-factor hierarchies: weight pool members by the subject prior
  # conditional on the factor scores ("conditional") or with the scores
  # integrated out ("marginal"). Both are exact; "conditional" is the default
  # because it gave equal group-level efficiency and a better worst subject on
  # the factor and infinite-factor benchmarks (2026-09-27). Internal switch,
  # not a user setting.
  latent <- getOption("emc2.ensemble_latent_weight", "conditional")
  if (!identical(latent, "marginal")) latent <- "conditional"
  list(enabled = isTRUE(getOption("emc2.ensemble", TRUE)),
       pool = max(2L, int_opt("emc2.ensemble_pool", 32L)),
       rounds = int_opt("emc2.ensemble_rounds", 10L),
       cadence = int_opt("emc2.ensemble_cadence", 1L),
       latent_weight = latent,
       tail_weight = 0.25, tail_df = 4, tail_scale = 4)
}

# Supported hierarchies: those whose group Gibbs step and subject prior are
# already used by the particle update. Nuisance and marginalised coordinates
# have their own samplers and are not yet covered.
.emc_ensemble_supported <- function(sampler) {
  isTRUE(sampler$type %in% c("standard", "diagonal-gamma", "factor",
                             "infnt_factor", "SEM")) &&
    !any(sampler$nuisance) &&
    !any(as.logical(sampler$marginalised_idx), na.rm = TRUE) &&
    isTRUE(sampler$n_subjects >= 2L)
}

# One-iteration copy of the sample store holding the current group state, so
# gibbs_step() can read it through last_sample_<type>() and fill_samples() can
# write it without touching the chain's history (fill_samples writes in place).
.emc_ensemble_store <- function(samples) {
  n_iter <- length(samples$stage)
  idx <- samples$idx
  out <- lapply(samples, function(x) {
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
# (1 - w) N(mu, V) + w t_df(mu, c V), with Mahalanobis distances reused for
# both components.
.emc_ensemble_q <- function(mu, var, cfg) {
  L <- t(chol(var))
  list(mu = mu, L = L, logdet = 2 * sum(log(diag(L))), p = length(mu))
}

.emc_ensemble_log_q <- function(x, q, cfg) {
  z <- forwardsolve(q$L, t(x) - q$mu)
  .emc_ensemble_log_q_m2(colSums(z^2), q, cfg)
}

# The same density from squared Mahalanobis distances m2.
.emc_ensemble_log_q_m2 <- function(m2, q, cfg) {
  p <- q$p; nu <- cfg$tail_df; c <- cfg$tail_scale
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
  if (any(tail)) scale[tail] <- sqrt(cfg$tail_scale * cfg$tail_df /
                                       stats::rchisq(sum(tail), cfg$tail_df))
  x <- t(q$mu + (q$L %*% z) * rep(scale, each = q$p))
  attr(x, "m2") <- scale^2 * colSums(z^2)
  x
}

# Subject prior of every subject under one group state h, as a p x n matrix of
# means and one covariance shared by all subjects. For the latent-factor
# hierarchies h contains the factor scores eta, and there are two valid
# weights:
#   conditional: a_s | h ~ N(mu + K x_s + Lambda eta_s, Sigma_err), which is
#     the density the group Gibbs step conditions on;
#   marginal: eta_s integrated out, a_s | h ~ N(mu + (K + Lambda B0^-1 G) x_s,
#     Lambda B0^-1 Delta^-1 B0^-T Lambda' + Sigma_err). This leaves eta stale
#     after reselection, which is valid because every latent-factor Gibbs step
#     redraws eta from its full conditional before anything conditions on it
#     (partially collapsed Gibbs, van Dyk & Park 2008).
# Factor and infinite-factor models have no covariates (K = G = 0, B = 0).
.emc_ensemble_subject_prior <- function(pars, sampler, n, weight) {
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
  if (weight == "marginal") return(list(mean = rep_mu(mu), var = pars$tvar))
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

# Subject prior log densities of every pool member under one group state.
.emc_ensemble_log_g <- function(pool_stack, subj_of_row, pars, sampler, n, weight) {
  sp <- .emc_ensemble_subject_prior(pars, sampler, n, weight)
  if (!is.null(sp$dvar)) {                     # diagonal covariance: no factorisation
    z <- (t(pool_stack) - sp$mean[, subj_of_row, drop = FALSE]) / sqrt(sp$dvar)
    return(-0.5 * (length(sp$dvar) * log(2 * pi) + sum(log(sp$dvar)) + colSums(z^2)))
  }
  L <- t(chol(sp$var))
  z <- forwardsolve(L, t(pool_stack) - sp$mean[, subj_of_row, drop = FALSE])
  -0.5 * (nrow(L) * log(2 * pi) + 2 * sum(log(diag(L))) + colSums(z^2))
}

.emc_ensemble_select <- function(lw) {
  w <- exp(lw - max(lw))
  sample.int(length(w), 1L, prob = w)
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
# when the update does not run this iteration.
.emc_ensemble_draw_pools <- function(sampler, stage, chains_mu, chains_var,
                                     iteration = 0L) {
  cfg <- .emc_ensemble_options()
  if (!cfg$enabled || stage == "preburn" || !.emc_ensemble_supported(sampler) ||
      iteration %% cfg$cadence != 0L)
    return(NULL)
  n <- sampler$n_subjects
  if (length(chains_mu) != n || length(chains_var) != n ||
      any(vapply(chains_mu, is.null, logical(1))) ||
      any(vapply(chains_var, is.null, logical(1))))
    return(NULL)
  pn <- sampler$par_names
  q <- lapply(seq_len(n), function(s)
    .emc_ensemble_q(as.numeric(chains_mu[[s]]), as.matrix(chains_var[[s]]), cfg))
  lq <- vector("list", n)
  rows <- lapply(seq_len(n), function(s) {
    fresh <- .emc_ensemble_draw_q(cfg$pool - 1L, q[[s]], cfg)
    lq[[s]] <<- .emc_ensemble_log_q_m2(attr(fresh, "m2"), q[[s]], cfg)
    attr(fresh, "m2") <- NULL
    colnames(fresh) <- pn
    fresh
  })
  list(q = q, rows = rows, lq = lq, ll = NULL)
}

# Log-likelihoods of a matrix of candidate states for one subject. A failure
# is a zero-density candidate: it is counted under the failure policy and
# returned as -Inf (strict mode stops instead).
.emc_ensemble_ll_checked <- function(props, data, model, marginalise = NULL) {
  ll <- tryCatch(as.numeric(calc_ll_manager(props, dadm = data, model = model,
                                             r_cores = 1L, marginalise = marginalise)),
    error = function(e) {
      cls <- .emc_classify_failure(e)
      stop(errorCondition(paste("Ensemble pool likelihood evaluation failed:",
        conditionMessage(e)), class = c("emc_ensemble_failure", "emc_failure"),
        emc_class = cls, emc_source = "ensemble", emc_cause = e))
    })
  if (length(ll) != nrow(props) || anyNA(ll) || any(ll == Inf)) {
    cls <- if (length(ll) != nrow(props)) "programming" else "numerical"
    stop(errorCondition(
      "Invalid ensemble pool likelihood: expected finite values or -Inf",
      class = c("emc_ensemble_failure", "emc_failure"),
      emc_class = cls, emc_source = "ensemble"))
  }
  ll
}

.emc_ensemble_ll_candidate <- function(props, data, model, marginalise = NULL) {
  tryCatch(.emc_ensemble_ll_checked(props, data, model, marginalise),
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
.emc_ensemble_pool_ll <- function(sampler, pools) {
  n <- sampler$n_subjects
  M1 <- nrow(pools$rows[[1L]])
  ll <- pools$ll
  if (is.null(ll)) ll <- vector("list", n)
  for (s in seq_len(n)) if (is.null(ll[[s]]) || length(ll[[s]]) != M1)
    ll[[s]] <- .emc_ensemble_ll_candidate(pools$rows[[s]], sampler$data[[s]],
                                          sampler$model, sampler$marginalise)
  ll_new <- matrix(unlist(ll), nrow = M1)
  if (anyNA(ll_new) || any(ll_new == Inf))
    stop("Invalid ensemble pool likelihood")
  ll_new
}

.emc_ensemble_move <- function(sampler, pars, proposals, stage, wpool, wpool_ctx,
                               wpool_part, pools) {
  keep <- list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool)
  if (is.null(pools)) return(keep)
  cfg <- .emc_ensemble_options()
  n <- sampler$n_subjects
  n_pars <- sampler$n_pars
  M <- cfg$pool
  pn <- sampler$par_names
  ll_new <- .emc_ensemble_pool_ll(sampler, pools)

  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  # Pools: current state at a uniform position among the M - 1 fresh draws.
  # Stacked pools, subject s in rows (s - 1) * M + 1:M, with the current state
  # at a uniform position k0.
  pool_stack <- matrix(NA_real_, n * M, n_pars, dimnames = list(NULL, pn))
  ll_stack <- lq_stack <- numeric(n * M)
  k0 <- integer(n)
  for (s in seq_len(n)) {
    k0[s] <- sample.int(M, 1L)
    rows <- (s - 1L) * M + seq_len(M)
    cur <- rows[k0[s]]
    pool_stack[cur, ] <- proposals[seq_len(n_pars), s]
    pool_stack[rows[-k0[s]], ] <- pools$rows[[s]]
    ll_stack[cur] <- proposals[n_pars + 1L, s]
    ll_stack[rows[-k0[s]]] <- ll_new[, s]
    lq_stack[cur] <- .emc_ensemble_log_q(pool_stack[cur, , drop = FALSE], pools$q[[s]], cfg)
    lq_stack[rows[-k0[s]]] <- pools$lq[[s]]
  }
  base_w <- ll_stack - lq_stack                  # h-independent part of log w
  subj_of_row <- rep(seq_len(n), each = M)
  offset <- (seq_len(n) - 1L) * M

  # Gibbs on pi~(h, k | A): h | k with the package's group step, then k | h.
  # The scratch store only needs the group state: gibbs_step() takes the
  # subject states as an argument.
  tmp <- sampler
  store <- .emc_ensemble_store(sampler$samples)
  selected <- k0
  sel_alpha <- function() t(pool_stack[offset + selected, , drop = FALSE])
  wess <- numeric(n)
  for (r in seq_len(cfg$rounds)) {
    alpha_k <- sel_alpha()
    store <- fill_samples(samples = store, group_level = pars, proposals = NULL,
                          j = 1L, n_pars = n_pars, type = sampler$type)
    tmp$samples <- store
    pars <- gibbs_step(tmp, alpha_k, sampler$type)
    lw_all <- matrix(base_w + .emc_ensemble_log_g(pool_stack, subj_of_row, pars, sampler,
                                                  n, cfg$latent_weight), nrow = M)
    selected <- .emc_ensemble_select_all(lw_all)
    if (r == cfg$rounds) {
      w <- exp(sweep(lw_all, 2, apply(lw_all, 2, max)))
      wess <- colSums(w)^2 / colSums(w^2)
    }
  }

  moved <- selected != k0
  proposals[seq_len(n_pars), ] <- sel_alpha()
  proposals[n_pars + 1L, ] <- ll_stack[offset + selected]
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

#' Ensemble-move diagnostics
#'
#' Mean share of subjects whose state changed per ensemble sweep, and the mean
#' normalised weight ESS of the pools under the final group state.
#' @param emc An emc object.
#' @return A data frame with one row per chain.
#' @noRd
ensemble_diagnostics <- function(emc) {
  emc <- restore_duplicates(emc)
  do.call(rbind, lapply(seq_along(emc), function(i) {
    st <- emc[[i]]$ensemble_stats
    if (is.null(st) || !st$sweeps) return(data.frame(chain = i, sweeps = 0L,
                                                     moved = NA_real_, weight_ess = NA_real_))
    data.frame(chain = i, sweeps = st$sweeps, moved = st$moved / st$sweeps,
               weight_ess = st$wess / st$sweeps)
  }))
}

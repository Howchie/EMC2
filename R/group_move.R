# Chain-local adaptation for the exact interweaving transport. The supported
# kernels are the legacy move, regularized adaptive Metropolis (AM), and
# robust adaptive Metropolis (RAM).
#
# ram: Robust Adaptive Metropolis (Vihola 2012, Stat. Comput. 22:997-1008).
# AM adapts to the posterior covariance of theta, but the move is a Metropolis
# step on the conditional target of theta given the whitened subject
# residuals; RAM adapts the proposal factor S from the move's own acceptance
# probabilities, so it learns the geometry of the distribution the move
# actually targets. Update: S S' <- S (I + eta_t (alpha - a*) u u'/|u|^2) S'.
# To learn faster per move, each adapting move draws `ram_batch` directions
# u_k from the current state (the first is the Metropolis candidate; the rest
# are evaluated in the same likelihood batch and only feed adaptation) and
# applies them as consecutive RAM steps. RAM's fixed point is unchanged. Adaptation is
# confined to the discarded adapt stage and freezes before production.
.emc_group_move_config <- function() {
  cfg <- .emc_group_move_options()
  if (cfg$method == "legacy") return(cfg)
  raw_k <- getOption("emc2.group_move_proposals", 1L)
  if (length(raw_k) != 1L || !is.numeric(raw_k) || is.na(raw_k) || raw_k != 1)
    stop("AM group moves require emc2.group_move_proposals = 1")
  positive <- function(name, default, integer = FALSE) {
    x <- getOption(paste0("emc2.group_move_", name), default)
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0 ||
        (integer && x != floor(x))) stop("Invalid group move option: ", name)
    x
  }
  cfg$target <- positive("target", .234)
  if (cfg$target >= 1) stop("Group move target must be below one")
  cfg$warmup <- getOption("emc2.group_move_warmup", NULL)
  if (!is.null(cfg$warmup)) cfg$warmup <- positive("warmup", NULL, TRUE)
  cfg$settle <- getOption("emc2.group_move_settle", NULL)
  if (!is.null(cfg$settle)) cfg$settle <- positive("settle", NULL, TRUE)
  cfg$max_dim <- positive("max_dim", 256L, TRUE)
  if (cfg$method == "ram") cfg$ram_batch <- positive("ram_batch", 32L, TRUE)
  cfg
}

.emc_group_move_signature <- function(sampler, spec) {
  # The design attribute is helper metadata. It contains model closures and
  # process-local pointers reconstructed across stages and serialization.
  prior <- sampler$prior
  attr(prior, "design") <- NULL
  list(spec = spec, prior = prior, gd = sampler$gd,
       marginal = sampler$marginalised_idx, nuisance = sampler$nuisance)
}

.emc_group_move_initialize <- function(sampler, spec, pars, cfg) {
  d <- spec$n_move
  if (d > cfg$max_dim)
    stop("Dense group move exceeds emc2.group_move_max_dim; use legacy or disable explicitly")
  baseline <- .emc_group_move_gibbs_cov(
    spec, sampler$prior, as.matrix(pars$tvar), sampler$gd,
    .emc_group_move_marginal(sampler, spec))
  if (any(!is.finite(baseline))) stop("Nonfinite group move initialization metric")
  shape_factor <- t(chol(baseline)) * (2.38 / sqrt(d))
  ram <- identical(cfg$method, "ram")
  # Adaptation is paid for in discarded iterations. Batched RAM averages
  # ram_batch increments per move, so it needs far fewer moves than AM.
  warmup <- if (!is.null(cfg$warmup)) cfg$warmup else if (ram) {
    max(200L, 10L * d)
  } else max(1000L, 50L * d)
  list(schema_version = 1L,
       kernel_version = if (ram) "interweave-ram-2" else "interweave-am-2",
       signature = .emc_group_move_signature(sampler, spec), config = cfg,
       status = "adapting", baseline = baseline,
       shape_factor = shape_factor, factor = shape_factor,
       warmup = as.integer(warmup),
       settling_remaining = if (!is.null(cfg$settle)) cfg$settle else
         if (ram) max(50L, 2L * d) else max(100L, 5L * d),
       n_adapt = 0L, am_n = 0L, am_mean = numeric(d),
       am_M2 = matrix(0, d, d), am_log_scale = 0,
       attempts_warmup = 0L, accepts_warmup = 0L,
       sum_accept_prob_warmup = 0,
       attempts_sample = 0L, accepts_sample = 0L,
       sum_accept_prob_sample = 0, freeze_iteration = NULL)
}

.emc_group_move_coordinates <- function(cur, spec) {
  c(cur$tmu[spec$mu_idx], .5 * log(diag(cur$tvar)[spec$scale_main]))
}

.emc_group_move_adapt <- function(gm, a, coordinates) {
  gm$n_adapt <- gm$n_adapt + 1L
  gm$am_n <- gm$am_n + 1L
  diff <- coordinates - gm$am_mean
  gm$am_mean <- gm$am_mean + diff / gm$am_n
  gm$am_M2 <- gm$am_M2 + tcrossprod(diff, coordinates - gm$am_mean)

  # Shape updates are deliberately less frequent; scalar feedback is applied
  # to the proposal factor every opportunity so it is never stale.
  if (gm$am_n > 1L && gm$am_n %% 25L == 0L) {
    d <- length(coordinates)
    weight <- gm$am_n / (gm$am_n + 5*d)
    empirical <- gm$am_M2 / (gm$am_n - 1L)
    covariance <- weight * empirical + (1 - weight) * gm$baseline
    gm$shape_factor <- t(chol((covariance + t(covariance)) / 2)) *
      (2.38 / sqrt(d))
  }
  gain <- min(.5, (gm$n_adapt + 10)^(-2/3))
  gm$am_log_scale <- gm$am_log_scale + gain * (a - gm$config$target)
  gm$factor <- gm$shape_factor * exp(gm$am_log_scale)
  if (any(!is.finite(gm$factor))) stop("Nonfinite adapted group move factor")
  if (gm$n_adapt >= gm$warmup) {
    gm$status <- "frozen_settling"
    gm$freeze_iteration <- gm$attempts_warmup
  }
  gm
}

# Batched RAM update. U holds the K standard-normal directions (columns) and
# a their acceptance probabilities from the current state. Each direction is
# one ordinary RAM step (Vihola 2012, rank-one factor update).
.emc_group_move_adapt_ram <- function(gm, U, a) {
  gm$n_adapt <- gm$n_adapt + 1L
  if (is.null(gm$ram_steps)) gm$ram_steps <- 0L
  d <- nrow(gm$factor)
  # The gain decays per move, shared by that move's K steps: the gain
  # conditions (sum eta = Inf, sum eta^2 < Inf over steps) still hold.
  eta <- min(.5, d * (gm$n_adapt + 10)^(-2/3))
  for (k in which(is.finite(a))) {
    gm$ram_steps <- gm$ram_steps + 1L
    v <- U[, k] / sqrt(sum(U[, k]^2))
    b <- eta * (a[k] - gm$config$target)
    x <- as.numeric(gm$factor %*% v)
    S <- tcrossprod(gm$factor) + b * tcrossprod(x)   # S (I + b v v') S'
    gm$factor <- t(chol((S + t(S)) / 2))
  }
  if (any(!is.finite(gm$factor))) stop("Nonfinite adapted group move factor")
  if (gm$n_adapt >= gm$warmup) {
    gm$status <- "frozen_settling"
    gm$freeze_iteration <- gm$attempts_warmup
  }
  gm
}

.emc_group_move_ready <- function(sampler) {
  gm <- sampler[["group_move"]]
  is.null(gm) || gm$status %in% c("ready", "sampling", "unsupported", "disabled")
}

# Worker evaluations are strict so failures return to the master for a retry.
.emc_group_move_ll_checked <- function(props, data, model, marginalise = NULL) {
  ll <- tryCatch(as.numeric(calc_ll_manager(props, dadm = data, model = model,
                                  r_cores = 1L, marginalise = marginalise)),
    error = function(e) {
      cls <- .emc_classify_failure(e)
      stop(errorCondition(paste("Group move likelihood evaluation failed:",
        conditionMessage(e)), class = c("emc_group_move_failure", "emc_failure"),
        emc_class = cls, emc_source = "group_move", emc_cause = e))
    })
  if (length(ll) != nrow(props) || anyNA(ll) || any(ll == Inf)) {
    cls <- if (length(ll) != nrow(props)) "programming" else "numerical"
    cond <- errorCondition(
      "Invalid group move likelihood: expected finite values or -Inf",
      class = c("emc_group_move_failure", "emc_failure"),
      emc_class = cls, emc_source = "group_move")
    stop(cond)
  }
  ll
}

# A candidate error is a rejection under report/silent policies. Strict mode
# preserves the cause and stops on nonnumerical failures.
.emc_group_move_ll_candidate <- function(props, data, model,
                                          marginalise = NULL) {
  tryCatch(.emc_group_move_ll_checked(props, data, model, marginalise),
    error = function(e) {
      cls <- if (inherits(e, "emc_failure") && !is.null(e$emc_class))
        e$emc_class else .emc_classify_failure(e)
      .emc_reject_record(e, "group_move")
      if (identical(.emc_failure_action(cls), "abort"))
        .emc_failure_abort(cls, "group_move", e)
      rep(-Inf, nrow(props))
    })
}

.emc_group_move <- function(sampler, pars, proposals, stage, wpool, wpool_ctx,
                            wpool_part, move_cores = 1L) {
  if (is.null(sampler$group_move_config))
    sampler$group_move_config <- .emc_group_move_config()
  cfg <- sampler$group_move_config
  if (is.null(cfg$method) || !cfg$method %in% c("legacy", "am", "ram"))
    stop("Stored group-move method is no longer supported; start a fresh fit with 'legacy', 'am' or 'ram'")
  if (cfg$method == "legacy") return(.emc_group_move_legacy(
    sampler, pars, proposals, stage, wpool, wpool_ctx, wpool_part, move_cores))
  keep <- function() list(sampler = sampler, pars = pars,
                          proposals = proposals, wpool = wpool)
  if (!cfg$enabled) {
    sampler$group_move <- list(status = "disabled")
    return(keep())
  }
  spec <- .emc_group_move_spec(sampler, cfg$scale)
  if (is.null(spec)) {
    sampler$group_move <- list(status = "unsupported", reason = paste(
      "No verified group transport for", sampler$type,
      "with these active coordinates"))
    return(keep())
  }
  if (!is.matrix(proposals) || nrow(proposals) != sampler$n_pars + 1L ||
      ncol(proposals) != sampler$n_subjects)
    stop("Invalid group move proposal dimensions")

  gm <- sampler[["group_move"]]
  if (is.null(gm)) {
    if (stage == "sample")
      stop("Adaptive group moves require discarded warmup before sampling")
    if (any(sampler$samples$stage == "sample"))
      stop("Starting an adaptive group move from a production segment requires a fresh fit")
    gm <- .emc_group_move_initialize(sampler, spec, pars, cfg)
  }
  signature <- .emc_group_move_signature(sampler, spec)
  stored_signature <- gm$signature
  # Normalize signatures from checkpoints written before design metadata was
  # excluded. Keep the mathematical prior parameters under exact comparison.
  if (!is.null(stored_signature$prior))
    attr(stored_signature$prior, "design") <- NULL
  if (!identical(gm$schema_version, 1L) || !identical(gm$config, cfg) ||
      !identical(stored_signature, signature))
    stop("Group move schema changed or state is missing; start a new discarded warmup epoch")
  gm$signature <- signature
  if (stage == "sample" && !gm$status %in% c("ready", "sampling"))
    stop("Group move warmup/settling incomplete; continue discarded adapt stage")
  if (stage != "sample" && gm$status == "sampling")
    stop("Cannot retune an existing production segment; start a fresh fit")
  if (stage == "sample") gm$status <- "sampling"

  n <- sampler$n_subjects
  cur <- list(tmu = as.numeric(pars$tmu), tvar = as.matrix(pars$tvar),
              tvinv = if (is.null(pars$tvinv)) solve(pars$tvar) else as.matrix(pars$tvinv),
              subj_mu = if (is.null(pars$subj_mu))
                matrix(pars$tmu, spec$p, n) else as.matrix(pars$subj_mu),
              alpha = proposals[spec$main_full, , drop = FALSE])
  P <- .emc_group_move_mu_invar(spec, sampler$prior)
  lp <- function(x) .emc_group_move_log_prior(
    x$tmu, x$tvar, spec, sampler$prior, P, pars$a_half)
  lp0 <- lp(cur)
  ll0 <- proposals[sampler$n_pars + 1L, ]
  if (!is.finite(lp0) || any(!is.finite(ll0)))
    stop("Nonfinite current group move target")
  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  shadow <- identical(cfg$method, "ram") && stage == "adapt" &&
    gm$status == "adapting"
  U <- matrix(stats::rnorm(spec$n_move * if (shadow) cfg$ram_batch else 1L),
              spec$n_move)
  u <- U[, 1L]
  delta <- as.numeric(gm$factor %*% u)
  cand <- .emc_group_move_apply(delta, spec, cur)
  is_valid <- function(x) all(is.finite(unlist(x))) && all(diag(x$tvar) > 0)
  valid <- is_valid(cand)
  a <- 0
  accepted <- FALSE

  # RAM: extra directions from the current state, evaluated in the
  # candidate's likelihood batch; they only feed adaptation.
  if (shadow && ncol(U) > 1L) {
    D <- gm$factor %*% U[, -1L, drop = FALSE]
    probe_states <- lapply(seq_len(ncol(D)), function(j)
      .emc_group_move_apply(D[, j], spec, cur))
    probe_ok <- vapply(probe_states, is_valid, logical(1))
  } else {
    probe_states <- list()
    probe_ok <- logical(0)
  }
  rows <- c(if (valid) list(cand), probe_states[probe_ok])

  if (length(rows)) {
    props <- lapply(seq_len(n), function(s) {
      m <- matrix(proposals[seq_len(sampler$n_pars), s], length(rows),
                  sampler$n_pars, byrow = TRUE,
                  dimnames = list(NULL, sampler$par_names))
      for (r in seq_along(rows)) m[r, spec$main_full] <- rows[[r]]$alpha[, s]
      m
    })
    ctx <- wpool_ctx
    ctx$group_move_checked <- TRUE
    pooled <- .emc_wpool_group_move_ll(wpool, wpool_part, props, ctx)
    if (is.null(pooled)) {
      ll_all <- vapply(seq_len(n), function(s)
        .emc_group_move_ll_candidate(props[[s]], sampler$data[[s]],
          sampler$model, sampler$marginalise), numeric(length(rows)))
      ll_all <- matrix(ll_all, nrow = length(rows))
    } else {
      ll_all <- matrix(pooled$ll, nrow = length(rows))
      wpool <- pooled$pool
    }
    if (anyNA(ll_all) || any(ll_all == Inf))
      stop("Invalid group move worker likelihood")
  }

  jac_of <- function(dl) if (spec$n_s)
    sum(dl[spec$n_t + seq_len(spec$n_s)] * spec$jac_weight) else 0
  shadow_a <- numeric(0)
  if (length(probe_states)) {
    shadow_a <- rep(0, length(probe_states))
    sums <- rowSums(ll_all)[seq_len(sum(probe_ok)) + as.integer(valid)]
    lr <- sums - sum(ll0) + vapply(probe_states[probe_ok], lp, numeric(1)) - lp0 +
      apply(D[, probe_ok, drop = FALSE], 2L, jac_of)
    shadow_a[probe_ok] <- exp(pmin(0, ifelse(is.na(lr), -Inf, lr)))
  }

  if (valid) {
    lp1 <- lp(cand)
    if (is.na(lp1) || lp1 == Inf) stop("Invalid proposed group prior density")
    ll1 <- ll_all[1L, ]
    jac <- if (spec$n_s)
      sum(delta[spec$n_t + seq_len(spec$n_s)] * spec$jac_weight) else 0
    log_ratio <- sum(ll1) - sum(ll0) + lp1 - lp0 + jac
    if (is.na(log_ratio) || log_ratio == Inf)
      stop("Invalid group move target ratio")
    a <- exp(min(0, log_ratio))
  }
  accepted <- stats::runif(1L) < a
  if (accepted) {
    pars$tmu <- cand$tmu
    pars$tvar <- cand$tvar
    pars$tvinv <- cand$tvinv
    if (!is.null(pars$subj_mu)) pars$subj_mu <- cand$subj_mu
    if (!is.null(pars$alpha)) pars$alpha <- cand$alpha
    proposals[spec$main_full, ] <- cand$alpha
    proposals[sampler$n_pars + 1L, ] <- ll1
  }

  suffix <- if (stage == "sample") "sample" else "warmup"
  attempts <- paste0("attempts_", suffix)
  accepts <- paste0("accepts_", suffix)
  sum_prob <- paste0("sum_accept_prob_", suffix)
  gm[[attempts]] <- gm[[attempts]] + 1L
  gm[[accepts]] <- gm[[accepts]] + as.integer(accepted)
  gm[[sum_prob]] <- gm[[sum_prob]] + a

  # Preburn and burn transitions use the initial Gibbs metric but never enter
  # AM's empirical moments or Robbins-Monro scale. Only the adapt stage learns.
  if (shadow) {
    gm <- .emc_group_move_adapt_ram(gm, U, c(a, shadow_a))
  } else if (stage == "adapt" && gm$status == "adapting") {
    gm <- .emc_group_move_adapt(
      gm, a, .emc_group_move_coordinates(if (accepted) cand else cur, spec))
  } else if (stage == "adapt" && gm$status == "frozen_settling") {
    gm$settling_remaining <- gm$settling_remaining - 1L
    if (gm$settling_remaining <= 0L) gm$status <- "ready"
  }
  sampler$group_move <- gm
  sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  keep()
}

#' Inspect group-move adaptation and frozen sampling state
#'
#' Reports algorithm state, not posterior convergence. Use posterior diagnostics
#' for means, scales, correlations and subject parameters separately.
#' @param emc An EMC fit or list of sampler chains.
#' @return A data frame with one row per chain.
#' @export
group_move_diagnostics <- function(emc) {
  emc <- restore_duplicates(emc)
  rows <- lapply(seq_along(emc), function(i) {
    s <- emc[[i]]
    cfg <- s$group_move_config
    gm <- s[["group_move"]]
    legacy <- s$group_move_state
    value <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
    ratio <- function(a, n) if (is.null(n) || n == 0) NA_real_ else a/n
    data.frame(chain = i,
      method = if (is.null(cfg)) "uninitialized" else cfg$method,
      status = if (!is.null(gm$status)) gm$status else
        if (is.null(legacy)) "uninitialized" else "legacy",
      transition = "Metropolis",
      dimension = if (!is.null(gm$factor)) nrow(gm$factor) else value(legacy$dim),
      n_adapt = if (is.null(gm)) value(legacy$adapt_iter) else value(gm$n_adapt),
      settling_remaining = value(gm$settling_remaining),
      warmup_attempts = value(gm$attempts_warmup),
      warmup_acceptance = ratio(gm$accepts_warmup, gm$attempts_warmup),
      sample_attempts = value(gm$attempts_sample),
      sample_acceptance = ratio(gm$accepts_sample, gm$attempts_sample),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

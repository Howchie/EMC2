# Non-Gaussian and dimension stress tests for legacy and regular AM.
suppressPackageStartupMessages(pkgload::load_all(quiet = TRUE))
if (!requireNamespace("posterior", quietly = TRUE))
  stop("Install posterior for benchmark diagnostics")

args <- commandArgs(TRUE)
smoke <- length(args) > 0L && identical(args[[1L]], "smoke")
dimensions <- length(args) > 0L && identical(args[[1L]], "dimensions")
short_run <- smoke
warmup <- if (smoke) 220L else 6000L
preburn <- if (short_run) 25L else 500L
settle <- if (short_run) 20L else 500L
production <- if (short_run) 100L else 4000L
n_chains <- if (short_run) 2L else 4L
n_subjects <- 4L
methods <- if (dimensions) {
  c("disabled", "legacy", "am", "oracle")
} else c("disabled", "legacy", "am")
families <- if (smoke) "student_t" else if (dimensions) {
  c("gaussian32", "gaussian64")
} else c("student_t", "funnel", "two_mode")
seeds <- if (smoke) 17001L else if (dimensions) {
  c(47001L, 57001L, 67001L)
} else c(17001L, 27001L, 37001L)
result_set <- if (smoke) "stress-smoke" else if (dimensions) "dimensions" else "stress"
outdir <- file.path("benchmarks/group-move/results", result_set)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fingerprint <- paste(tools::md5sum(c("R/group_move.R", "R/sampling.R",
  "benchmarks/group-move/stress.R")), warmup, settle, production,
  n_chains, n_subjects, paste(methods, collapse = ","), paste(families, collapse = ","),
  sep = "|")

make_target <- function(family, seed) {
  set.seed(seed)
  if (startsWith(family, "gaussian")) {
    d <- as.integer(sub("^gaussian", "", family))
    Q <- qr.Q(qr(matrix(stats::rnorm(d * d), d)))
    eigenvalues <- exp(seq(log(0.01), log(100), length.out = d))
    covariance <- Q %*% diag(eigenvalues, d) %*% t(Q)
    precision <- Q %*% diag(1 / eigenvalues, d) %*% t(Q)
    log_normalizer <- -0.5 * (d * log(2 * pi) + sum(log(eigenvalues)))
    log_density <- function(x) log_normalizer -
      0.5 * drop(crossprod(x, precision %*% x))
    draw <- function() as.numeric(t(chol(covariance)) %*% stats::rnorm(d))
    transform <- function(x) as.numeric(t(Q) %*% x / sqrt(eigenvalues))
    list(d = d, log_density = log_density, draw = draw,
         covariance = covariance,
         transform = transform,
         description = paste0(d, "D rotated Gaussian; condition number 10,000"))
  } else if (family == "student_t") {
    d <- 8L
    Q <- qr.Q(qr(matrix(stats::rnorm(d * d), d)))
    scale <- Q %*% diag(exp(seq(log(0.1), log(10), length.out = d)), d) %*% t(Q)
    nu <- 5
    log_density <- function(x) {
      q <- drop(crossprod(x, solve(scale, x)))
      lgamma((nu + d) / 2) - lgamma(nu / 2) -
        d / 2 * log(nu * pi) - 0.5 * as.numeric(determinant(scale, logarithm = TRUE)$modulus) -
        (nu + d) / 2 * log1p(q / nu)
    }
    draw <- function() {
      as.numeric(t(chol(scale)) %*% stats::rnorm(d) /
                   sqrt(stats::rchisq(1, nu) / nu))
    }
    covariance <- nu / (nu - 2) * scale
    transform <- function(x) t(Q) %*% x / sqrt(diag(t(Q) %*% covariance %*% Q))
    list(d = d, log_density = log_density, draw = draw,
         transform = transform, description = "8D rotated Student-t, df=5")
  } else if (family == "funnel") {
    d <- 8L
    log_density <- function(x) {
      stats::dnorm(x[1], 0, 3, log = TRUE) -
        (d - 1) / 2 * (log(2 * pi) + x[1]) -
        0.5 * sum(x[-1]^2) * exp(-x[1])
    }
    draw <- function() {
      v <- stats::rnorm(1, 0, 3)
      c(v, stats::rnorm(d - 1) * exp(v / 2))
    }
    transform <- function(x) c(x[1] / 3, x[-1] / exp(2.25))
    list(d = d, log_density = log_density, draw = draw,
         transform = transform, description = "8D Neal funnel")
  } else {
    d <- 2L
    angle <- 0.43
    Q <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L)
    direction <- Q[, 1L]
    orthogonal <- Q[, 2L]
    means <- list(-6 * direction, 6 * direction)
    covariance <- Q %*% diag(c(0.5, 1), 2L) %*% t(Q)
    log_component <- function(x, mu) {
      z <- x - mu
      -log(2 * pi) - 0.5 * as.numeric(determinant(covariance, logarithm = TRUE)$modulus) -
        0.5 * drop(crossprod(z, solve(covariance, z)))
    }
    log_density <- function(x) {
      a <- log_component(x, means[[1L]])
      b <- log_component(x, means[[2L]])
      m <- max(a, b)
      m + log(exp(a - m) + exp(b - m)) - log(2)
    }
    draw_component <- function(sign) {
      as.numeric(sign * 6 * direction +
                   t(chol(covariance)) %*% stats::rnorm(d))
    }
    draw <- function() draw_component(if (stats::runif(1) < 0.5) -1 else 1)
    transform <- function(x) c(drop(crossprod(direction, x)) / sqrt(36.5),
                                drop(crossprod(orthogonal, x)))
    list(d = d, log_density = log_density, draw = draw,
         draw_component = draw_component,
         transform = transform, description = "2D rotated equal-weight Gaussian mixture; modes separated by 12")
  }
}

target_ll <- function(proposals, dadm, model, ...) {
  d <- dadm$d
  n <- dadm$n_subjects
  values <- apply(proposals[, seq_len(d), drop = FALSE], 1L, function(x) {
    # Cancel the Gaussian group-mean prior so the translated state has exactly
    # the stated target density. Every subject contributes one nth of this term.
    dadm$log_density(as.numeric(x)) + 0.5 * dadm$prior_precision * sum(x^2)
  })
  values / n
}

run_case <- function(family, seed, method) {
  target <- make_target(family, seed)
  d <- target$d
  check_point <- seq(-0.7, 0.9, length.out = d)
  check_data <- list(d = d, n_subjects = n_subjects,
                     prior_precision = 0.01,
                     log_density = target$log_density)
  check_props <- matrix(check_point, nrow = 1L,
                        dimnames = list(NULL, paste0("x", seq_len(d))))
  encoded_target <- n_subjects * target_ll(check_props, check_data, NULL) -
    0.5 * 0.01 * sum(check_point^2)
  if (!isTRUE(all.equal(encoded_target, target$log_density(check_point),
                        tolerance = 1e-12)))
    stop("Synthetic likelihood does not encode the declared target")
  set.seed(seed + 1000L)
  proposal_count <- if (method == "legacy") 4L else 1L
  configured_method <- if (method %in% c("disabled", "oracle")) "am" else method
  options(emc2.group_move_method = configured_method,
          emc2.group_move = method != "disabled",
          emc2.group_move_proposals = proposal_count,
          emc2.group_move_scale = FALSE,
          emc2.group_move_warmup = warmup,
          emc2.group_move_settle = settle)
  history <- array(NA_real_, c(d, n_chains, preburn + warmup + settle + production))
  legacy_tv_history <- if (method == "legacy")
    array(rep(diag(d), 250L), c(d, d, 250L)) else NULL
  states <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    set.seed(seed + ch)
    initial <- if (family == "two_mode")
      target$draw_component(if (ch %% 2L) -1 else 1) else target$draw()
    sampler <- list(type = "standard", n_pars = d, n_subjects = n_subjects,
      par_names = paste0("x", seq_len(d)), nuisance = rep(FALSE, d),
      marginalised_idx = rep(FALSE, d), is_blocked = rep(TRUE, d),
      par_group = rep(1L, d),
      prior = list(theta_mu_mean = rep(0, d), theta_mu_invar = diag(0.01, d)),
      data = rep(list(list(d = d, n_subjects = n_subjects,
                           prior_precision = 0.01,
                           log_density = target$log_density)), n_subjects))
    if (method == "oracle") {
      if (is.null(target$covariance)) stop("Oracle reference requires a Gaussian target")
      sampler$group_move_config <- EMC2:::.emc_group_move_config()
      spec <- EMC2:::.emc_group_move_spec(sampler, scale = FALSE)
      gm <- EMC2:::.emc_group_move_initialize(
        sampler, spec, list(tvar = diag(d)), sampler$group_move_config)
      oracle_factor <- t(chol(target$covariance)) * (2.38 / sqrt(d))
      gm$baseline <- target$covariance
      gm$shape_factor <- oracle_factor
      gm$factor <- oracle_factor
      gm$status <- "ready"
      sampler$group_move <- gm
    }
    states[[ch]] <- list(sampler = sampler,
      alpha = matrix(initial, d, n_subjects), rng = .Random.seed)
  }
  started <- proc.time()[["elapsed"]]
  total <- preburn + warmup + settle + production
  for (iteration in seq_len(total)) {
    for (ch in seq_len(n_chains)) {
      state <- states[[ch]]
      assign(".Random.seed", state$rng, globalenv())
      sampler <- state$sampler
      alpha <- state$alpha
      pars <- list(tmu = as.numeric(alpha[, 1L]), tvar = diag(d),
                   tvinv = diag(d))
      current_loglik <- (target$log_density(pars$tmu) +
                           0.5 * 0.01 * sum(pars$tmu^2)) / n_subjects
      likelihood_row <- rep(current_loglik, n_subjects)
      sampler$rng$gibbs <- .Random.seed
      stage <- if (method == "oracle") "sample" else
        if (iteration <= preburn) "preburn" else
          if (iteration <= preburn + warmup + settle) "adapt" else "sample"
      step <- EMC2:::.emc_group_move(sampler, pars,
        rbind(alpha, likelihood_row), stage, NULL, NULL, NULL)
      sampler <- step$sampler
      alpha <- step$proposals[seq_len(d), , drop = FALSE]
      states[[ch]] <- list(sampler = sampler, alpha = alpha,
                           rng = .Random.seed)
      history[, ch, iteration] <- step$pars$tmu
    }
    if (method == "legacy" && iteration <= preburn + warmup && iteration %% 100L == 0L) {
      ix <- max(1L, iteration - 249L):iteration
      chains <- lapply(seq_len(n_chains), function(ch) {
        s <- states[[ch]]$sampler
        s$samples <- list(theta_mu = matrix(history[, ch, ix], d),
          theta_var = legacy_tv_history)
        s
      })
      spec <- EMC2:::.emc_group_move_spec(chains[[1L]], FALSE)
      V <- EMC2:::.emc_group_move_covariance(chains, spec, seq_along(ix))
      if (!is.null(V)) for (ch in seq_len(n_chains))
        states[[ch]]$sampler <- EMC2:::.emc_group_move_install(
          states[[ch]]$sampler, V, spec, proposal_count)
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  transformed <- array(NA_real_, c(production, n_chains, d))
  for (i in seq_len(production)) for (ch in seq_len(n_chains))
    transformed[i, ch, ] <- target$transform(history[, ch, preburn + warmup + settle + i])
  metric_rows <- lapply(seq_len(d), function(j) {
    z <- transformed[, , j]
    data.frame(family = family, description = target$description,
      seed = seed, method = method, coordinate = j,
      rhat = posterior::rhat(z), bulk_ess = posterior::ess_bulk(z),
      tail_ess = posterior::ess_tail(z), mean_standardized = mean(z),
      variance_standardized = var(as.numeric(z)),
      elapsed_seconds = elapsed,
      likelihood_evaluations = if (method == "disabled") 0L else
        total * n_chains * n_subjects * proposal_count,
      ess_per_100k_likelihood_evaluations = if (method == "disabled") NA_real_ else
        posterior::ess_bulk(z) / (total * n_chains * n_subjects * proposal_count) * 1e5,
      stringsAsFactors = FALSE)
  })
  mode_occupancy <- if (family == "two_mode") {
    direction <- c(cos(0.43), sin(0.43))
    vapply(seq_len(n_chains), function(ch) {
      mean(vapply(seq_len(production), function(i) {
        sum(history[, ch, preburn + warmup + settle + i] * direction) > 0
      }, logical(1)))
    }, numeric(1))
  } else rep(NA_real_, n_chains)
  diagnostics <- do.call(rbind, lapply(seq_along(states), function(ch) {
    s <- states[[ch]]
    gm <- s$sampler$group_move
    legacy <- s$sampler$group_move_state
    adaptive_warmup <- if (length(gm$attempts_warmup) == 1L &&
        gm$attempts_warmup > 0) gm$accepts_warmup / gm$attempts_warmup else NULL
    adaptive_sample <- if (length(gm$attempts_sample) == 1L &&
        gm$attempts_sample > 0) gm$accepts_sample / gm$attempts_sample else NULL
    data.frame(chain = ch, method = method,
      warmup_acceptance = if (!is.null(adaptive_warmup)) adaptive_warmup else
        if (!is.null(legacy)) legacy$moves / max(legacy$attempts, 1) else NA_real_,
      sample_acceptance = if (!is.null(adaptive_sample)) adaptive_sample else NA_real_,
      adaptation_opportunities = if (length(gm$n_adapt) == 1L) gm$n_adapt else
        if (!is.null(legacy)) legacy$adapt_iter else 0,
      stringsAsFactors = FALSE)
  }))
  result_metrics <- do.call(rbind, metric_rows)
  result_metrics$warmup_acceptance <- mean(diagnostics$warmup_acceptance)
  result_metrics$sample_acceptance <- mean(diagnostics$sample_acceptance)
  result_metrics$adaptation_opportunities <- min(diagnostics$adaptation_opportunities)
  for (ch in seq_len(n_chains))
    result_metrics[[paste0("mode_occupancy_chain_", ch)]] <-
      if (family == "two_mode") mode_occupancy[ch] else NA_real_
  list(metrics = result_metrics, mode_occupancy = mode_occupancy,
       diagnostics = diagnostics, elapsed = elapsed, draws = transformed)
}

rows <- list()
for (family in families) for (seed in seeds) for (method in methods) {
  key <- sprintf("%s-seed%d-%s", family, seed, method)
  cache <- file.path(outdir, paste0(key, ".rds"))
  result <- if (file.exists(cache)) readRDS(cache) else NULL
  if (is.null(result) || !identical(result$code_fingerprint, fingerprint)) {
    result <- testthat::with_mocked_bindings(
      run_case(family, seed, method), calc_ll_manager = target_ll,
      .package = "EMC2")
    result$code_fingerprint <- fingerprint
    saveRDS(result, cache)
  }
  rows[[length(rows) + 1L]] <- result$metrics
  all_columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned_rows <- lapply(rows, function(x) {
    for (name in setdiff(all_columns, names(x))) x[[name]] <- NA_real_
    x[all_columns]
  })
  write.csv(do.call(rbind, aligned_rows), file.path(outdir, "metrics.csv"),
            row.names = FALSE)
  cat(key, "worst Rhat", round(max(result$metrics$rhat), 3),
      "minimum bulk ESS", round(min(result$metrics$bulk_ess)),
      "variance range",
      paste(round(range(result$metrics$variance_standardized), 2), collapse = "–"),
      "seconds", round(result$elapsed, 1),
      "mode occupancy",
      if (family == "two_mode") paste(round(result$mode_occupancy, 2), collapse = ",") else "n/a",
      "\n")
  flush.console()
}
writeLines(c(capture.output(sessionInfo()),
             capture.output(system("git rev-parse HEAD", intern = TRUE))),
           file.path(outdir, "session.txt"))

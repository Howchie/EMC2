# Compare a factor-four rank-one scale on the allocator's slow LBA shear
# direction with a no-collapsed control. Charge the exact pilot to x4's time.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
n_draws <- as.integer(Sys.getenv("COLL_DIRECTION_DRAWS", "200"))
n_pilots <- as.integer(Sys.getenv("COLL_DIRECTION_PILOTS", "64"))
factor <- 4
state <- readRDS("/tmp/lba_sampling_audit/shared_allocation_current.rds")
base <- EMC2:::restore_custom_kernel_pointers(
  structure(state$chains, class = "emc"), quiet = TRUE)
options(emc2.subject_hmc = FALSE, emc2.collapsed_move = TRUE)
for (j in seq_along(base)) {
  base[[j]]$production_options$emc2.subject_hmc <- FALSE
  base[[j]]$production_options$emc2.collapsed_move <- TRUE
  base[[j]]$sample_proposals_next <- Inf
  base[[j]]$collapsed$plan$rounds <- 0
}

directional <- unserialize(serialize(base, NULL))
chain <- directional[[2L]]
plan <- chain$collapsed$plan
ctx <- EMC2:::.emc_collapsed_context(chain)
p <- chain$n_pars
k <- 1L
D <- chain$shear$plan$dirs
d <- D[, k]
d[startsWith(chain$par_names, "t0")] <- 0
stopifnot(sum(d^2) > 0)
d <- d / sqrt(sum(d^2))
current_C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
  plan$t0_ratio, EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, ctx$chart),
  plan$t0_hyper_ratio,
  EMC2:::.emc_collapsed_t0_mean_chart_idx(chain$par_names),
  plan$mean_idx, plan$mean_ratio)
T <- diag(ctx$chart$dim)
T[seq_len(p), seq_len(p)] <- diag(p) + (factor - 1) * tcrossprod(d)
Sigma <- T %*% current_C %*% t(current_C) %*% t(T)
C_directional <- t(chol((Sigma + t(Sigma)) / 2))
idx <- chain$samples$idx
phi <- EMC2:::.emc_collapsed_to_phi(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$chart)
plan$C <- C_directional
plan$scale <- 1
plan$t0_ratio <- 1
plan$t0_mean_ratio <- 1
plan$t0_hyper_ratio <- 1
plan$mean_idx <- integer()
plan$mean_ratio <- 1
plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi, C_directional,
  plan$eps, ctx, plan$length)

# Measure the one-time exact-pilot cost on private proposals. Restore the
# caller's RNG and leave the retained chain streams untouched.
had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
on.exit_rng <- function() {
  if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
    suppressWarnings(rm(".Random.seed", envir = globalenv()))
}
set.seed(20261010L + factor)
pilot_seeds <- sample.int(.Machine$integer.max, n_pilots)
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
pilot_base <- list(cur = cur,
  alpha = matrix(chain$samples$alpha[, , idx], p),
  ll = chain$samples$subj_ll[, idx], rounds = 0L, wpool = NULL,
  wpool_part = NULL)
pilot_acceptance <- numeric(n_pilots)
pilot_target_esjd <- numeric(n_pilots)
Dinv <- solve(D)
target_var <- stats::var(drop(Dinv[k, , drop = FALSE] %*%
  chain$samples$theta_mu[, EMC2:::.emc_alloc_draws(chain, 200L), drop = FALSE]))
pilot_started <- proc.time()[["elapsed"]]
for (i in seq_len(n_pilots)) {
  set.seed(pilot_seeds[i])
  x <- pilot_base
  x$moved <- FALSE
  x <- EMC2:::.emc_collapsed_round(x, plan, chain, ctx)
  pilot_acceptance[i] <- x$acc
  jump <- sum(Dinv[k, ] * x$jump[seq_len(p)])
  pilot_target_esjd[i] <- x$acc * jump^2
}
pilot_seconds <- proc.time()[["elapsed"]] - pilot_started
on.exit_rng()

plan$rounds <- 4
plan$calibrate <- FALSE
plan$calibration_screen <- FALSE
plan$last <- EMC2:::.emc_collapsed_record(NULL, p)
plan$pool <- EMC2:::.emc_collapsed_record(NULL, p)
plan$block <- EMC2:::.emc_collapsed_record(NULL, p)
directional[[2L]]$collapsed$plan <- plan

summarise_fit <- function(fit, arm, sample_seconds) {
  nch <- length(fit)
  p <- fit[[1L]]$n_pars
  nsub <- fit[[1L]]$n_subjects
  par_names <- fit[[1L]]$par_names
  target_Dinv <- solve(fit[[2L]]$shear$plan$dirs)
  group <- array(NA_real_, c(n_draws, nch, p),
    dimnames = list(NULL, NULL, par_names))
  target <- array(NA_real_, c(n_draws, nch, 1L),
    dimnames = list(NULL, NULL, "slow_shear_1"))
  logsd <- array(NA_real_, c(n_draws, nch, p),
    dimnames = list(NULL, NULL, paste0("logsd:", par_names)))
  subject <- array(NA_real_, c(n_draws, nch, p * nsub))
  ids <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  subject_id <- match("52091", ids)
  subject_target <- array(NA_real_, c(n_draws, nch, p),
    dimnames = list(NULL, NULL, par_names))
  for (j in seq_len(nch)) {
    ch <- fit[[j]]
    ii <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"),
               n_draws)
    stopifnot(length(ii) == n_draws)
    gm <- t(ch$samples$theta_mu[, ii, drop = FALSE])
    group[, j, ] <- gm
    target[, j, 1L] <- drop(gm %*% t(target_Dinv[1L, , drop = FALSE]))
    logsd[, j, ] <- t(vapply(ii, function(t)
      EMC2:::.emc_alloc_lsd(ch$samples$theta_var[, , t]), numeric(p)))
    alpha <- array(ch$samples$alpha[, , ii], c(p, nsub, n_draws))
    for (t in seq_len(n_draws)) subject[t, j, ] <- as.vector(alpha[, , t])
    if (!is.na(subject_id)) subject_target[, j, ] <- t(matrix(
      ch$samples$alpha[, subject_id, ii], nrow = p, ncol = n_draws))
  }
  get_summaries <- function(a) {
    z <- posterior::summarise_draws(posterior::as_draws_array(a),
      posterior::ess_bulk, posterior::ess_tail, posterior::rhat)
    ess_col <- grep("ess_bulk", names(z), value = TRUE)[1L]
    rhat_col <- grep("rhat", names(z), value = TRUE)[1L]
    list(ess = z[[ess_col]], rhat = z[[rhat_col]], variable = z$variable)
  }
  families <- list(group_mu = get_summaries(group),
    slow_shear_1 = get_summaries(target), group_logsd = get_summaries(logsd),
    subject_52091 = get_summaries(subject_target), all_subject = get_summaries(subject))
  total <- sample_seconds + if (arm == "direction_x4") pilot_seconds else 0
  block <- lapply(fit, function(ch) ch$collapsed$plan$block)
  family_summary <- do.call(rbind, lapply(names(families), function(nm) {
    z <- families[[nm]]
    data.frame(arm = arm, family = nm,
      min_bulk_ess = min(z$ess, na.rm = TRUE),
      median_bulk_ess = median(z$ess, na.rm = TRUE),
      mean_bulk_ess = mean(z$ess, na.rm = TRUE),
      min_ess_per_second = min(z$ess, na.rm = TRUE) / total,
      median_ess_per_second = median(z$ess, na.rm = TRUE) / total,
      max_rhat = max(z$rhat, na.rm = TRUE),
      frac_rhat_over_1_05 = mean(z$rhat > 1.05, na.rm = TRUE))
  }))
  list(summary = family_summary,
    per_variable = lapply(families, function(z)
      data.frame(variable = z$variable, ess_bulk = z$ess, rhat = z$rhat)),
    timing = data.frame(arm = arm, sample_seconds = sample_seconds,
      pilot_seconds = if (arm == "direction_x4") pilot_seconds else 0,
      total_seconds_charged = total,
      collapsed_rounds = paste(vapply(fit, function(ch)
        ch$collapsed$plan$rounds, numeric(1)), collapse = ";"),
      collapsed_acceptance = if (sum(vapply(block, `[[`, 0L, "n")) > 0)
        sum(vapply(block, function(b) b$acc, 0)) /
          sum(vapply(block, `[[`, 0L, "n")) else NA_real_))
}

arms <- list(direction_x4 = directional, no_collapsed = base)
results <- list()
for (arm in names(arms)) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(arms[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = n_draws,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  sample_seconds <- proc.time()[["elapsed"]] - started
  results[[arm]] <- summarise_fit(fit, arm, sample_seconds)
  cat(arm, "sample seconds", round(sample_seconds, 1),
    "pilot seconds", round(results[[arm]]$timing$pilot_seconds, 1),
    "rounds", results[[arm]]$timing$collapsed_rounds,
    "acceptance", round(results[[arm]]$timing$collapsed_acceptance, 3), "\n")
  rm(fit)
  gc()
}
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
timing <- do.call(rbind, lapply(results, `[[`, "timing"))
per_variable <- lapply(results, `[[`, "per_variable")
write.csv(summary, file.path(artifact_dir,
  "allocator_direction_x4_summary.csv"), row.names = FALSE)
write.csv(timing, file.path(artifact_dir,
  "allocator_direction_x4_timing.csv"), row.names = FALSE)
saveRDS(list(factor = factor, n_draws = n_draws, n_pilots = n_pilots,
  pilot_seconds = pilot_seconds, pilot_acceptance = mean(pilot_acceptance),
  pilot_target_esjd_per_variance = mean(pilot_target_esjd) / target_var,
  summary = summary, timing = timing, per_variable = per_variable),
  file.path(artifact_dir, "allocator_direction_x4.rds"))
print(summary, row.names = FALSE)
print(timing, row.names = FALSE)

pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
state <- readRDS("/tmp/lba_sampling_audit/shared_allocation_current.rds")
base <- EMC2:::restore_custom_kernel_pointers(
  structure(state$chains, class = "emc"), quiet = TRUE)
n_draws <- 200L
pilot_seconds <- 22.881 # 64 exact proposals for the x2 direction pilot, one core
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
d <- d / sqrt(sum(d^2))
current_C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
  plan$t0_ratio, EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, ctx$chart),
  plan$t0_hyper_ratio,
  EMC2:::.emc_collapsed_t0_mean_chart_idx(chain$par_names),
  plan$mean_idx, plan$mean_ratio)
T <- diag(ctx$chart$dim)
T[seq_len(p), seq_len(p)] <- diag(p) + tcrossprod(d)
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
plan$rounds <- 4
plan$calibrate <- FALSE
plan$calibration_screen <- FALSE
plan$last <- EMC2:::.emc_collapsed_record(NULL, p)
plan$pool <- EMC2:::.emc_collapsed_record(NULL, p)
plan$block <- EMC2:::.emc_collapsed_record(NULL, p)
chain$collapsed$plan <- plan
directional[[2L]] <- chain

summarise_fit <- function(fit, arm, elapsed) {
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
  }
  ess <- function(a) posterior::summarise_draws(
    posterior::as_draws_array(a), posterior::ess_bulk, posterior::rhat)
  ge <- ess(group)
  te <- ess(target)
  le <- ess(logsd)
  se <- ess(subject)
  ec <- grep("ess_bulk", names(ge), value = TRUE)[1L]
  tc <- grep("ess_bulk", names(te), value = TRUE)[1L]
  lc <- grep("ess_bulk", names(le), value = TRUE)[1L]
  sc <- grep("ess_bulk", names(se), value = TRUE)[1L]
  rc <- grep("rhat", names(ge), value = TRUE)[1L]
  max_rhat <- if (length(rc) && !is.na(rc))
    max(unlist(lapply(list(ge, te, le, se), `[[`, rc), use.names = FALSE),
        na.rm = TRUE) else NA_real_
  total <- elapsed + if (arm == "direction_x2") pilot_seconds else 0
  block <- lapply(fit, function(ch) ch$collapsed$plan$block)
  data.frame(arm = arm, sample_seconds = elapsed,
    directional_pilot_seconds = if (arm == "direction_x2") pilot_seconds else 0,
    total_seconds_charged = total,
    collapsed_rounds = paste(vapply(fit, function(ch)
      ch$collapsed$plan$rounds, numeric(1)), collapse = ";"),
    collapsed_acceptance = if (sum(vapply(block, `[[`, 0L, "n")) > 0)
      sum(vapply(block, function(b) b$acc, 0)) /
        sum(vapply(block, `[[`, 0L, "n")) else NA_real_,
    target_ess = te[[tc]][1L], target_ess_per_second = te[[tc]][1L] / total,
    median_group_mu_ess = median(ge[[ec]]),
    median_group_mu_ess_per_second = median(ge[[ec]]) / total,
    median_group_logsd_ess = median(le[[lc]]),
    median_group_logsd_ess_per_second = median(le[[lc]]) / total,
    median_subject_ess = median(se[[sc]]),
    median_subject_ess_per_second = median(se[[sc]]) / total,
    max_rhat = max_rhat)
}

arms <- list(direction_x2 = directional, no_collapsed = base)
results <- list()
for (arm in c("direction_x2", "no_collapsed")) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(arms[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  results[[arm]] <- list(summary = summarise_fit(fit, arm, elapsed))
  saveRDS(lapply(results, `[[`, "summary"), file.path(artifact_dir,
    "allocator_direction_ess_partial.rds"))
  cat(arm, "seconds", round(elapsed, 1), "rounds",
    results[[arm]]$summary$collapsed_rounds,
    "acceptance", round(results[[arm]]$summary$collapsed_acceptance, 3),
    "target ESS/sec", round(results[[arm]]$summary$target_ess_per_second, 3),
    "group mean ESS/sec", round(results[[arm]]$summary$median_group_mu_ess_per_second, 3),
    "subject ESS/sec", round(results[[arm]]$summary$median_subject_ess_per_second, 3), "\n")
  rm(fit)
  gc()
}
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
write.csv(summary, file.path(artifact_dir,
  "allocator_direction_ess.csv"), row.names = FALSE)
saveRDS(summary, file.path(artifact_dir, "allocator_direction_ess.rds"))
print(summary, row.names = FALSE)

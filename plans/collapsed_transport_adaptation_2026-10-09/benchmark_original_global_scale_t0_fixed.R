pkgload::load_all(".", quiet = TRUE, compile = FALSE)
options(emc2.failure_policy = "strict")

fit_dir <- "/tmp/collapsed_transport_followup"
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
chain <- readRDS(file.path(fit_dir, "safe_fit_ess_start.rds"))$start[[2L]]
selected <- readRDS(file.path(fit_dir, "mean_scaled_calibrated.rds"))

ctx <- EMC2:::.emc_collapsed_context(chain)
idx <- chain$samples$idx
p <- chain$n_pars
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
stopifnot(!is.null(cur), is.finite(cur$lp))
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
base_plan <- chain$collapsed$plan
base_C <- base_plan$C
base <- list(cur = cur, alpha = matrix(chain$samples$alpha[, , idx], p),
  ll = chain$samples$subj_ll[, idx], rounds = 0L, wpool = NULL,
  wpool_part = NULL)

selected_plan <- selected$collapsed$plan
t0_idx <- base_plan$t0_idx
t0_mean_idx <- base_plan$t0_mean_idx
t0_hyper_idx <- setdiff(t0_idx, t0_mean_idx)
target_idx <- match("v.g", chain$par_names)
stopifnot(!is.na(target_idx), target_idx <= chain$n_pars,
  !target_idx %in% t0_mean_idx)

# Preserve the independently selected effective t0 scales in every arm.
# The current frozen plan has scale=.044194, t0 mean=.044194, and t0
# hierarchy=.011049. The unit-scale arm therefore changes only non-t0 rows.
t0_mean_scale <- selected_plan$scale * selected_plan$t0_ratio
t0_hyper_scale <- selected_plan$scale * selected_plan$t0_hyper_ratio
target_scale <- selected_plan$scale * selected_plan$mean_ratio
arms <- list(
  unit_global_t0_fixed = list(scale = 1, mean_scale = 1, mean_idx = integer()),
  old_global_t0_fixed = list(scale = base_plan$scale,
    mean_scale = base_plan$scale, mean_idx = integer()),
  selected_directional = list(scale = selected_plan$scale,
    mean_scale = target_scale, mean_idx = target_idx)
)

n_pilot <- 64L
seed0 <- 20261011L
set.seed(seed0)
pilot_seed_start <- get(".Random.seed", envir = globalenv())
posterior_draws <- EMC2:::.emc_alloc_draws(chain, 200L)
target_var <- stats::var(chain$samples$theta_mu[target_idx, posterior_draws])

result <- vector("list", length(arms))
for (i in seq_along(arms)) {
  arm_name <- names(arms)[i]
  arm <- arms[[i]]
  plan <- base_plan
  plan$scale <- arm$scale
  plan$t0_ratio <- t0_mean_scale / arm$scale
  plan$t0_mean_ratio <- plan$t0_ratio
  plan$t0_hyper_ratio <- t0_hyper_scale / arm$scale
  plan$mean_idx <- arm$mean_idx
  plan$mean_ratio <- if (length(arm$mean_idx)) arm$mean_scale / arm$scale else 1
  plan$length <- base_plan$length
  C_candidate <- EMC2:::.emc_collapsed_preconditioner(base_C, plan$scale,
    plan$t0_ratio, t0_idx, plan$t0_hyper_ratio, t0_mean_idx,
    plan$mean_idx, plan$mean_ratio)

  assign(".Random.seed", pilot_seed_start, envir = globalenv())
  plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C_candidate,
    base_plan$eps, ctx, plan$length)
  pilot_seed <- get(".Random.seed", envir = globalenv())
  acceptance <- total_esjd <- target_esjd <- target_norm_esjd <- numeric(n_pilot)
  start <- proc.time()[["elapsed"]]
  for (k in seq_len(n_pilot)) {
    assign(".Random.seed", pilot_seed, envir = globalenv())
    x <- base
    x$moved <- FALSE
    x <- EMC2:::.emc_collapsed_round(x, plan, chain, ctx)
    pilot_seed <- get(".Random.seed", envir = globalenv())
    acceptance[k] <- x$acc
    total_esjd[k] <- x$acc * sum(forwardsolve(base_C, x$phi_jump)^2)
    target_esjd[k] <- x$acc * x$jump[target_idx]^2
    target_norm_esjd[k] <- target_esjd[k] / target_var
  }
  elapsed <- proc.time()[["elapsed"]] - start
  result[[i]] <- data.frame(
    arm = arm_name, global_scale = plan$scale,
    t0_mean_scale = plan$scale * plan$t0_ratio,
    t0_hyper_scale = plan$scale * plan$t0_hyper_ratio,
    target_mean_scale = if (length(plan$mean_idx))
      plan$scale * plan$mean_ratio else NA_real_,
    length = plan$length, eps = plan$eps,
    expected_exact_acceptance = mean(acceptance),
    total_C_metric_esjd = mean(total_esjd),
    v_g_esjd = mean(target_esjd),
    v_g_esjd_per_posterior_variance = mean(target_norm_esjd),
    seconds = elapsed, proposals = n_pilot,
    proposals_per_second = n_pilot / elapsed
  )
  cat(arm_name, "scale", plan$scale, "t0 scales", t0_mean_scale,
      t0_hyper_scale, "eps", plan$eps, "acceptance", mean(acceptance),
      "total ESJD", mean(total_esjd), "v.g ESJD/var",
      mean(target_norm_esjd), "seconds", elapsed, "\n")
  flush.console()
}

result <- do.call(rbind, result)
write.csv(result,
  file.path(artifact_dir, "original_global_scale_t0_fixed.csv"),
  row.names = FALSE)
saveRDS(result,
  file.path(artifact_dir, "original_global_scale_t0_fixed.rds"))
print(result, row.names = FALSE)

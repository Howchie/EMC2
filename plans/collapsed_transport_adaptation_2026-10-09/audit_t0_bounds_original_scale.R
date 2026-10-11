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
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
base_plan <- chain$collapsed$plan
base_C <- base_plan$C
base <- list(cur = cur, alpha = matrix(chain$samples$alpha[, , idx], p),
  ll = chain$samples$subj_ll[, idx], rounds = 0L, wpool = NULL,
  wpool_part = NULL)
selected_plan <- selected$collapsed$plan
t0_mean_scale <- selected_plan$scale * selected_plan$t0_ratio
t0_hyper_scale <- selected_plan$scale * selected_plan$t0_hyper_ratio
target_idx <- match("v.g", chain$par_names)
arms <- list(
  unit_global_t0_rows_fixed = list(scale = 1, mean_idx = integer()),
  unit_global_t0_group_frozen = list(scale = 1, mean_idx = integer(),
    freeze_t0 = TRUE),
  old_global_t0_rows_fixed = list(scale = base_plan$scale, mean_idx = integer()),
  selected_directional = list(scale = selected_plan$scale, mean_idx = target_idx)
)

# Count proposed, rather than accepted, subject-level t0 support crossings.
# The wrapper delegates to the package likelihood evaluator unchanged.
bound_rows <- list()
capture_bounds <- function(rows, sampler) {
  lower <- upper <- logical(length(rows))
  for (s in seq_along(rows)) {
    pars <- EMC2:::get_pars_oo(rows[[s]], sampler$data[[s]], sampler$model)
    t0 <- pars[, "t0"]
    rt <- sampler$data[[s]]$rt
    stopifnot(length(t0) == length(rt))
    lower[s] <- any(is.finite(t0) & t0 < 0.05)
    upper[s] <- any(is.finite(t0) & is.finite(rt) & t0 >= rt)
  }
  bound_rows[[length(bound_rows) + 1L]] <<- c(lower = sum(lower), upper = sum(upper))
  invisible(NULL)
}
ns <- asNamespace("EMC2")
binding <- ".emc_wpool_subject_ll"
old_ll <- get(binding, envir = ns, inherits = FALSE)
unlockBinding(binding, ns)
assign(binding, function(pool, rows, ctx, part = NULL, source) {
  if (identical(source, "collapsed")) capture_bounds(rows, ctx)
  old_ll(pool, rows, ctx, part, source)
}, envir = ns)
lockBinding(binding, ns)

n_pilot <- 64L
set.seed(20261011L)
seed0 <- get(".Random.seed", envir = globalenv())
result <- list()
for (arm_name in names(arms)) {
  arm <- arms[[arm_name]]
  plan <- base_plan
  plan$scale <- arm$scale
  plan$t0_ratio <- if (isTRUE(arm$freeze_t0)) 1 else t0_mean_scale / arm$scale
  plan$t0_mean_ratio <- plan$t0_ratio
  plan$t0_hyper_ratio <- if (isTRUE(arm$freeze_t0)) 1 else
    t0_hyper_scale / arm$scale
  if (isTRUE(arm$freeze_t0)) plan$C[base_plan$t0_idx, ] <- 0
  plan$mean_idx <- arm$mean_idx
  plan$mean_ratio <- if (length(arm$mean_idx))
    (selected_plan$scale * selected_plan$mean_ratio) / arm$scale else 1
  plan$length <- base_plan$length
  C_candidate <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
    plan$t0_ratio, base_plan$t0_idx, plan$t0_hyper_ratio,
    base_plan$t0_mean_idx, plan$mean_idx, plan$mean_ratio)
  assign(".Random.seed", seed0, envir = globalenv())
  plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C_candidate,
    base_plan$eps, ctx, plan$length)
  seed_pilots <- get(".Random.seed", envir = globalenv())
  bound_rows <- list()
  acc <- numeric(n_pilot)
  for (k in seq_len(n_pilot)) {
    assign(".Random.seed", seed_pilots, envir = globalenv())
    x <- base
    x <- EMC2:::.emc_collapsed_round(x, plan, chain, ctx)
    seed_pilots <- get(".Random.seed", envir = globalenv())
    acc[k] <- x$acc
  }
  b <- do.call(rbind, bound_rows)
  result[[arm_name]] <- data.frame(
    arm = arm_name, proposals = n_pilot,
    proposals_any_lower_violation = sum(b[, "lower"] > 0),
    proposals_any_rt_violation = sum(b[, "upper"] > 0),
    mean_subjects_below_005 = mean(b[, "lower"]),
    mean_subjects_t0_ge_rt = mean(b[, "upper"]),
    mean_acceptance_on_violation = mean(acc[(b[, "lower"] > 0) |
      (b[, "upper"] > 0)]),
    mean_acceptance_without_violation = mean(acc[(b[, "lower"] == 0) &
      (b[, "upper"] == 0)]),
    mean_acceptance = mean(acc)
  )
  cat(arm_name, "lower-cross proposals",
      result[[arm_name]]$proposals_any_lower_violation,
      "RT-cross proposals", result[[arm_name]]$proposals_any_rt_violation,
      "mean acceptance", mean(acc), "\n")
}
result <- do.call(rbind, result)
write.csv(result, file.path(artifact_dir, "t0_bounds_original_scale_audit.csv"),
  row.names = FALSE)
print(result, row.names = FALSE)

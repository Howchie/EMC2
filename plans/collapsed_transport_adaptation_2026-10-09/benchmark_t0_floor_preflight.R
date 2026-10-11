pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
target <- Sys.getenv("EMC2_PREFLIGHT_TARGET", "LBA")
if (identical(target, "LBA")) {
  prepared <- readRDS("/tmp/collapsed_transport_followup/safe_fit_ess_start.rds")
  chain <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(prepared$start), quiet = TRUE)[[2L]]
} else if (identical(target, "T30")) {
  archive <- Sys.getenv("EMC2_FULLC_RDS",
    ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final/T30_fullc_1_emc.rds")
  chain <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(readRDS(archive)), quiet = TRUE)[[1L]]
} else stop("EMC2_PREFLIGHT_TARGET must be LBA or T30")
ctx <- EMC2:::.emc_collapsed_context(chain)
idx <- chain$samples$idx
p <- chain$n_pars
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
x0 <- list(cur = cur, alpha = matrix(chain$samples$alpha[, , idx], p),
  ll = chain$samples$subj_ll[, idx], rounds = 0L, wpool = NULL,
  wpool_part = NULL)

plan <- chain$collapsed$plan
plan$scale <- 1
plan$t0_ratio <- plan$t0_mean_ratio <- plan$t0_hyper_ratio <- 1
plan$mean_idx <- integer()
plan$mean_ratio <- 1
plan$length <- if (isTRUE(is.finite(plan$length)) && plan$length > 0)
  plan$length else 3
plan$t0_idx <- EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, ctx$chart)
plan$t0_mean_idx <- EMC2:::.emc_collapsed_t0_mean_chart_idx(chain$par_names)
C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
  plan$t0_ratio, plan$t0_idx, plan$t0_hyper_ratio, plan$t0_mean_idx,
  plan$mean_idx, plan$mean_ratio)
plan$eps <- EMC2:::.emc_with_preserved_rng(
  EMC2:::.emc_collapsed_tune_eps(phi0, C, plan$eps, ctx, plan$length))

ns <- asNamespace("EMC2")
rounds_per_arm <- as.integer(Sys.getenv("EMC2_PREFLIGHT_REPS", "12"))
stopifnot(is.finite(rounds_per_arm), rounds_per_arm > 0L)
preflight_binding <- ".emc_collapsed_t0_below_bound"
ll_binding <- ".emc_wpool_subject_ll"
old_preflight <- get(preflight_binding, envir = ns, inherits = FALSE)
old_ll <- get(ll_binding, envir = ns, inherits = FALSE)
enabled <- TRUE
preflight_count <- 0L
skip_count <- 0L
exact_rows <- 0L
exact_elapsed <- exact_cpu <- 0
unlockBinding(preflight_binding, ns)
unlockBinding(ll_binding, ns)
assign(preflight_binding, function(alpha, sampler, support = NULL) {
  preflight_count <<- preflight_count + 1L
  rejected <- if (enabled) old_preflight(alpha, sampler, support) else FALSE
  skip_count <<- skip_count + as.integer(rejected)
  rejected
}, envir = ns)
assign(ll_binding, function(pool, rows, ctx, part = NULL, source) {
  started <- proc.time()
  result <- old_ll(pool, rows, ctx, part, source)
  elapsed <- proc.time() - started
  if (identical(source, "collapsed")) {
    exact_rows <<- exact_rows + length(rows)
    exact_elapsed <<- exact_elapsed + unname(elapsed[["elapsed"]])
    exact_cpu <<- exact_cpu + unname(elapsed[["user.self"]] + elapsed[["sys.self"]])
  }
  result
}, envir = ns)
lockBinding(preflight_binding, ns)
lockBinding(ll_binding, ns)

run_arm <- function(label, use_preflight, seed = 20261010L) {
  enabled <<- use_preflight
  arm_plan <- plan
  arm_plan$t0_floor_preflight <- use_preflight
  preflight_count <<- skip_count <<- exact_rows <<- 0L
  exact_elapsed <<- exact_cpu <<- 0
  set.seed(seed)
  outputs <- vector("list", rounds_per_arm)
  total_cpu <- 0
  started <- proc.time()[["elapsed"]]
  for (k in seq_along(outputs)) {
    x <- x0
    x$moved <- FALSE
    cpu0 <- proc.time()
    out <- EMC2:::.emc_collapsed_round(x, arm_plan, chain, ctx)
    cpu <- proc.time() - cpu0
    total_cpu <- total_cpu + unname(cpu[["user.self"]] + cpu[["sys.self"]])
    outputs[[k]] <- list(acc = out$acc, phi_jump = out$phi_jump,
      jump = out$jump, rng = get(".Random.seed", envir = globalenv()))
  }
  list(result = data.frame(arm = label, rounds = length(outputs),
    subjects = chain$n_subjects, preflight_checks = preflight_count,
    preflight_skips = skip_count, exact_subject_rows = exact_rows,
    exact_eval_wall_seconds = exact_elapsed, exact_eval_cpu_seconds = exact_cpu,
    full_round_wall_seconds = proc.time()[["elapsed"]] - started,
    full_round_cpu_seconds = total_cpu,
    mean_exact_acceptance = mean(vapply(outputs, `[[`, 0.0, "acc"))),
    outputs = outputs)
}

invisible(run_arm("warm_exact", FALSE))
invisible(run_arm("warm_preflight", TRUE))
reference <- run_arm("exact_all", FALSE)
screened <- run_arm("floor_preflight", TRUE)
unlockBinding(preflight_binding, ns)
unlockBinding(ll_binding, ns)
assign(preflight_binding, old_preflight, envir = ns)
assign(ll_binding, old_ll, envir = ns)
lockBinding(preflight_binding, ns)
lockBinding(ll_binding, ns)
same <- all(vapply(seq_along(reference$outputs), function(k)
  identical(reference$outputs[[k]]$acc, screened$outputs[[k]]$acc) &&
    identical(reference$outputs[[k]]$phi_jump, screened$outputs[[k]]$phi_jump) &&
    identical(reference$outputs[[k]]$jump, screened$outputs[[k]]$jump) &&
    identical(reference$outputs[[k]]$rng, screened$outputs[[k]]$rng), logical(1)))
if (!same) stop("preflight changed proposal, acceptance, or RNG results")
summary <- rbind(reference$result, screened$result)
summary$same_rng_proposals_and_decisions <- same
write.csv(summary, file.path(artifact_dir,
  paste0("t0_floor_preflight_cost_", tolower(target), ".csv")),
  row.names = FALSE)
print(summary, row.names = FALSE)

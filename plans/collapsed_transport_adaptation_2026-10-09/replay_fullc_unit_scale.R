pkgload::load_all(".", quiet = TRUE, compile = FALSE)

archive <- Sys.getenv("EMC2_FULLC_RDS",
  "/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final/T30_fullc_1_emc.rds")
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
n_proposals <- as.integer(Sys.getenv("EMC2_FULLC_PROPOSALS", "32"))
seed <- as.integer(Sys.getenv("EMC2_FULLC_SEED", "20261012"))
chain_id <- as.integer(Sys.getenv("EMC2_FULLC_CHAIN", "1"))

emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(archive)), quiet = TRUE)
chain <- emc[[chain_id]]
iteration <- chain$samples$idx
ctx <- EMC2:::.emc_collapsed_context(chain)
current <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, iteration],
  as.matrix(chain$samples$theta_var[, , iteration]), ctx$fits,
  chain$prior, ctx$spec)
stopifnot(!is.null(current), is.finite(current$lp))
x0 <- list(cur = current,
  alpha = matrix(chain$samples$alpha[, , iteration], chain$n_pars),
  ll = chain$samples$subj_ll[, iteration], rounds = 0L,
  wpool = NULL, wpool_part = NULL)

# The archive predates scalar calibration. Make its unshrunk metric explicit
# under the current plan schema, while retaining the archived eps and C.
plan <- chain$collapsed$plan
plan$scale <- 1
plan$t0_ratio <- 1
plan$t0_mean_ratio <- 1
plan$t0_hyper_ratio <- 1

violations <- list()
ns <- asNamespace("EMC2")
binding <- ".emc_wpool_subject_ll"
old_ll <- get(binding, envir = ns, inherits = FALSE)
unlockBinding(binding, ns)
assign(binding, function(pool, rows, ctx, part = NULL, source) {
  if (identical(source, "collapsed")) {
    lower <- upper <- logical(length(rows))
    for (s in seq_along(rows)) {
      pars <- EMC2:::get_pars_oo(rows[[s]], ctx$data[[s]], ctx$model)
      t0 <- pars[, "t0"]
      rt <- ctx$data[[s]]$rt
      stopifnot(length(t0) == length(rt))
      lower[s] <- any(is.finite(t0) & t0 < 0.05)
      upper[s] <- any(is.finite(t0) & is.finite(rt) & t0 >= rt)
    }
    violations[[length(violations) + 1L]] <<-
      c(lower = sum(lower), upper = sum(upper))
  }
  old_ll(pool, rows, ctx, part, source)
}, envir = ns)
lockBinding(binding, ns)

set.seed(seed)
acceptance <- numeric(n_proposals)
started <- proc.time()[["elapsed"]]
for (k in seq_len(n_proposals)) {
  x <- EMC2:::.emc_collapsed_round(x0, plan, chain, ctx)
  acceptance[k] <- x$acc
}
elapsed <- proc.time()[["elapsed"]] - started
unlockBinding(binding, ns)
assign(binding, old_ll, envir = ns)
lockBinding(binding, ns)
crossings <- do.call(rbind, violations)

result <- data.frame(case = "T30", seed = 1L, chain = chain_id,
  subjects = chain$n_subjects, group_parameters = chain$n_pars,
  chart_dim = plan$dim, scale = 1, eps = plan$eps,
  proposals = n_proposals, mean_expected_exact_acceptance = mean(acceptance),
  proposals_any_t0_below_005 = sum(crossings[, "lower"] > 0),
  proposals_any_t0_ge_observed_rt = sum(crossings[, "upper"] > 0),
  elapsed_seconds = elapsed)
write.csv(result,
  file.path(artifact_dir, "fullc_current_source_replay.csv"), row.names = FALSE)
print(result, row.names = FALSE)

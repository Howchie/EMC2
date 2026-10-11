pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
prepared <- readRDS("/tmp/collapsed_transport_followup/safe_fit_ess_start.rds")
make_start <- function() EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(prepared$start), quiet = TRUE)[[2L]]
ns <- asNamespace("EMC2")
preflight_binding <- ".emc_collapsed_t0_below_bound"
ll_binding <- ".emc_wpool_subject_ll"
old_preflight <- get(preflight_binding, envir = ns, inherits = FALSE)
old_ll <- get(ll_binding, envir = ns, inherits = FALSE)
enabled <- TRUE
preflight_checks <- preflight_skips <- exact_rows <- 0L
exact_cpu <- 0
unlockBinding(preflight_binding, ns)
unlockBinding(ll_binding, ns)
assign(preflight_binding, function(alpha, sampler, support = NULL) {
  preflight_checks <<- preflight_checks + 1L
  rejected <- if (enabled) old_preflight(alpha, sampler, support) else FALSE
  preflight_skips <<- preflight_skips + as.integer(rejected)
  rejected
}, envir = ns)
assign(ll_binding, function(pool, rows, ctx, part = NULL, source) {
  started <- proc.time()
  result <- old_ll(pool, rows, ctx, part, source)
  elapsed <- proc.time() - started
  if (identical(source, "collapsed")) {
    exact_rows <<- exact_rows + length(rows)
    exact_cpu <<- exact_cpu + unname(elapsed[["user.self"]] + elapsed[["sys.self"]])
  }
  result
}, envir = ns)
lockBinding(preflight_binding, ns)
lockBinding(ll_binding, ns)

run_arm <- function(label, use_preflight) {
  enabled <<- use_preflight
  preflight_checks <<- preflight_skips <<- exact_rows <<- 0L
  exact_cpu <<- 0
  chain <- make_start()
  chain$collapsed$plan$calibrate <- TRUE
  chain$collapsed$plan$calibration_screen <- FALSE
  if (is.null(chain$rng)) chain$rng <- list()
  set.seed(20261042L)
  chain$rng$gibbs <- get(".Random.seed", envir = globalenv())
  started <- proc.time()
  chain <- EMC2:::.emc_collapsed_calibrate(chain, n_cores = 1L,
    validation_batches = 16L)
  elapsed <- proc.time() - started
  list(summary = data.frame(arm = label, pilot_rounds = preflight_checks,
    preflight_skips = preflight_skips, exact_subject_rows = exact_rows,
    exact_eval_cpu_seconds = exact_cpu,
    calibration_wall_seconds = unname(elapsed[["elapsed"]]),
    calibration_cpu_seconds = unname(elapsed[["user.self"]] + elapsed[["sys.self"]]),
    scale = chain$collapsed$plan$scale,
    t0_mean_ratio = chain$collapsed$plan$t0_mean_ratio,
    t0_hyper_ratio = chain$collapsed$plan$t0_hyper_ratio,
    exact_acceptance = chain$collapsed$plan$calibration$acceptance),
    plan = chain$collapsed$plan, rng = chain$rng$gibbs)
}

if (identical(Sys.getenv("EMC2_PREFLIGHT_ORDER", "exact_first"), "screen_first")) {
  screened <- run_arm("floor_preflight", TRUE)
  reference <- run_arm("exact_all", FALSE)
} else {
  reference <- run_arm("exact_all", FALSE)
  screened <- run_arm("floor_preflight", TRUE)
}
unlockBinding(preflight_binding, ns)
unlockBinding(ll_binding, ns)
assign(preflight_binding, old_preflight, envir = ns)
assign(ll_binding, old_ll, envir = ns)
lockBinding(preflight_binding, ns)
lockBinding(ll_binding, ns)

same <- identical(reference$plan$scale, screened$plan$scale) &&
  identical(reference$plan$t0_mean_ratio, screened$plan$t0_mean_ratio) &&
  identical(reference$plan$t0_hyper_ratio, screened$plan$t0_hyper_ratio) &&
  identical(reference$plan$calibration$acceptance,
    screened$plan$calibration$acceptance) &&
  identical(reference$rng, screened$rng)
if (!same) stop("preflight changed calibration results or RNG state")
summary <- rbind(reference$summary, screened$summary)
summary$same_calibration_and_rng <- same
output <- if (identical(Sys.getenv("EMC2_PREFLIGHT_ORDER", "exact_first"),
    "screen_first")) "t0_floor_preflight_calibration_repeat2.csv" else
      "t0_floor_preflight_calibration.csv"
write.csv(summary, file.path(artifact_dir, output), row.names = FALSE)
print(summary, row.names = FALSE)

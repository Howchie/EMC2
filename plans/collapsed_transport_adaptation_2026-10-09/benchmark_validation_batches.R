pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
start_file <- Sys.getenv("VALIDATION_START",
  "/tmp/collapsed_transport_followup/length_adaptive_allocation_start.rds")
raw <- readRDS(start_file)
base <- if (is.list(raw$start)) raw$start else raw
base <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(base), quiet = TRUE)
options(emc2.subject_hmc = FALSE, emc2.failure_policy = "strict")
if (identical(Sys.getenv("VALIDATION_INITIAL_REBUILD", "0"), "1")) {
  for (j in seq_along(base))
    base[[j]]$sample_proposals_next <- base[[j]]$sample_transitions
  base <- EMC2:::add_proposals(base, stage = "sample", n_cores = 12L,
    n_blocks = 1L, block_iter = 200L)
}
production_rng <- lapply(base, `[[`, "rng")
production_rng <- lapply(production_rng, `[[`, "gibbs")

run_batch_count <- function(n_batches) {
  emc <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(unserialize(serialize(base, NULL))), quiet = TRUE)
  elapsed <- numeric(length(emc))
  reserved <- numeric(length(emc))
  for (j in seq_along(emc)) {
    old <- emc[[j]]$collapsed$plan
    emc[[j]] <- EMC2:::.emc_collapsed_set_plan(emc[[j]], old$C, "sample")
    stopifnot(isTRUE(emc[[j]]$collapsed$plan$calibration_screen))
    started <- proc.time()[["elapsed"]]
    emc[[j]] <- EMC2:::.emc_collapsed_calibrate(emc[[j]], n_cores = 4L,
      validation_batches = n_batches)
    elapsed[j] <- proc.time()[["elapsed"]] - started
    reserved[j] <- EMC2:::.emc_reserved_cost(elapsed[j], 4L)
  }

  Vbar <- lapply(emc, function(ch) {
    k <- EMC2:::.emc_alloc_draws(ch,
      if (is.null(ch$alloc$window)) 200L else max(2L, ch$alloc$window$n))
    apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean)
  })
  B <- 200L
  rows <- vector("list", length(emc))
  for (j in seq_along(emc)) {
    ms <- EMC2:::.emc_alloc_measure(emc, j, NULL, Vbar, "sample")
    ms$c_calibration <- EMC2:::.emc_alloc_amortized_cost(reserved[j], B)
    u <- EMC2:::.emc_alloc_units(ms, NULL, ms$shear, optimistic = FALSE)
    sol <- EMC2:::.emc_alloc_solve_units(u)
    p <- emc[[j]]$collapsed$plan
    emc[[j]]$collapsed$plan$rounds <- sol$joint
    rows[[j]] <- data.frame(
      batches = n_batches, chain = j, seconds = elapsed[j],
      proposals = p$calibration$pilots,
      scale = p$scale,
      t0_mean_scale = p$scale * p$t0_ratio,
      t0_hyper_scale = p$scale * p$t0_hyper_ratio,
      mean_scale = if (length(p$mean_idx)) p$scale * p$mean_ratio else NA_real_,
      path_length = p$length, eps = p$eps,
      validation_acceptance = p$calibration$acceptance,
      validation_score = p$calibration$log_score,
      marginal_cost = p$calibration$cost,
      setup_cost = p$calibration$fixed_cost,
      mean_collapsed_factor = mean(ms$c_fac),
      mean_saturation_factor = mean(ms$c_sat),
      rounds = sol$joint,
      objective = sol$objective,
      objective_without_collapsed = sol$objective_without_joint)
  }
  save_dir <- Sys.getenv("VALIDATION_SAVE_DIR", "")
  if (nzchar(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    for (j in seq_along(emc)) {
      emc[[j]]$rng$gibbs <- production_rng[[j]]
      emc[[j]]$sample_proposals_next <- Inf
    }
    saveRDS(emc, file.path(save_dir, paste0("screen", n_batches, ".rds")))
  }
  do.call(rbind, rows)
}

result <- rbind(run_batch_count(8L), run_batch_count(16L))
label <- if (grepl("T30", start_file, ignore.case = TRUE)) "t30_" else "lba_"
write.csv(result, file.path(artifact_dir,
  paste0(label, "validation_batch_comparison.csv")),
          row.names = FALSE)
print(result, row.names = FALSE)

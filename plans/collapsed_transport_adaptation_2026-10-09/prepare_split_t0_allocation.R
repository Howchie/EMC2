pkgload::load_all(".", quiet = TRUE, compile = FALSE)
fit_dir <- "/tmp/collapsed_transport_followup"
x <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
options(emc2.subject_hmc = FALSE)
for (j in seq_along(x)) {
  x[[j]]$production_options$emc2.subject_hmc <- FALSE
  x[[j]]$sample_proposals_next <- x[[j]]$sample_transitions
}
.captured <- list()
trace(".emc_alloc_solve_units", where = asNamespace("EMC2"), print = FALSE,
  exit = quote({
    .GlobalEnv$.captured[[length(.GlobalEnv$.captured) + 1L]] <-
      list(units = u, solution = returnValue())
  }))
started <- proc.time()[["elapsed"]]
x <- EMC2:::add_proposals(x, "sample", n_cores = 12L, n_blocks = 1L,
                          block_iter = 200L)
rebuild_seconds <- proc.time()[["elapsed"]] - started
untrace(".emc_alloc_solve_units", where = asNamespace("EMC2"))
stopifnot(length(.captured) == length(x))
for (j in seq_along(x)) {
  p <- x[[j]]$collapsed$plan
  a <- x[[j]]$alloc$last
  cat("chain", j, "scale", p$scale, "t0 mean ratio", p$t0_ratio,
      "t0 hyper ratio", p$t0_hyper_ratio, "length", p$length,
      "eps", p$eps, "pilot seconds", p$calibration$elapsed,
      "acceptance", p$calibration$acceptance,
      "collapsed rate", p$rounds, "objective", a$objective,
      "without collapsed", a$objective_without_collapsed, "\n")
}
cat("rebuild seconds", rebuild_seconds, "\n")
saveRDS(list(start = x, captured = .captured,
             rebuild_seconds = rebuild_seconds),
        file.path(fit_dir, "split_t0_allocation_start.rds"))

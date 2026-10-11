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
summary <- do.call(rbind, lapply(seq_along(x), function(j) {
  p <- x[[j]]$collapsed$plan
  a <- x[[j]]$alloc$last
  data.frame(chain = j, scale = p$scale,
    scale_step_capped = isTRUE(p$scale_step_capped),
    t0_mean_scale = p$scale * p$t0_ratio,
    t0_hyper_scale = p$scale * p$t0_hyper_ratio,
    slow_mean = if (length(p$mean_idx)) x[[j]]$par_names[p$mean_idx] else NA_character_,
    slow_mean_scale = p$scale * p$mean_ratio,
    length = p$length, eps = p$eps,
    pilot_seconds = p$calibration$elapsed,
    validation_acceptance = p$calibration$acceptance,
    proposed_rounds = p$rounds, objective = a$objective,
    objective_without_collapsed = a$objective_without_collapsed,
    objective0 = a$objective0)
}))
print(summary, row.names = FALSE)
cat("rebuild_seconds", rebuild_seconds, "\n")
saveRDS(list(start = x, captured = .captured,
             rebuild_seconds = rebuild_seconds, summary = summary),
        file.path(fit_dir, "capped_allocator_probe.rds"))

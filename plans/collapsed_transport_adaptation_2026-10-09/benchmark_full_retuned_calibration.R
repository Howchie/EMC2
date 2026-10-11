pkgload::load_all(".", quiet = TRUE, compile = FALSE)
prepared <- readRDS("/tmp/collapsed_transport_followup/t0_adaptive_allocation_start.rds")
emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(prepared$start), quiet = TRUE)
ch <- emc[[1L]]
options(emc2.failure_policy = "strict")
old <- ch$collapsed$plan
stopifnot(is.list(old$calibration), identical(old$stage, "sample"))
# Mark the previous plan as having run so this benchmark exercises the full
# grid; the retained state and starting metric remain those in the checkpoint.
old$rounds <- 1L
ch$collapsed$plan <- old
ch <- EMC2:::.emc_collapsed_set_plan(ch, old$C, "sample")
stopifnot(!isTRUE(ch$collapsed$plan$calibration_screen))
started <- proc.time()[["elapsed"]]
ch <- EMC2:::.emc_collapsed_calibrate(ch, n_cores = 4L)
elapsed <- proc.time()[["elapsed"]] - started
plan <- ch$collapsed$plan
out <- list(scale = plan$scale, t0_ratio = plan$t0_ratio, eps = plan$eps,
            acceptance = plan$calibration$acceptance, elapsed = elapsed,
            global_search = plan$calibration$search,
            t0_search = plan$calibration$t0_search)
saveRDS(out, "plans/collapsed_transport_adaptation_2026-10-09/full_retuned_calibration.rds")
write.csv(out$global_search,
          "plans/collapsed_transport_adaptation_2026-10-09/full_retuned_global_search.csv",
          row.names = FALSE)
write.csv(out$t0_search,
          "plans/collapsed_transport_adaptation_2026-10-09/full_retuned_t0_search.csv",
          row.names = FALSE)
cat("mode full scale", out$scale, "t0 ratio", out$t0_ratio,
    "eps", out$eps, "acceptance", out$acceptance,
    "elapsed", out$elapsed, "seconds\n")
print(out$global_search, row.names = FALSE)
print(out$t0_search, row.names = FALSE)

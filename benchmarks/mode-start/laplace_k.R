# Pareto k of the Laplace importance weights from the real burn mode start,
# for every fixture with a saved preburn state.
#   Rscript benchmarks/mode-start/laplace_k.R
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- list()
for (case in c("corr", "freeA", "wide", "lba71", "joint")) {
  f <- file.path(out, paste0("preburn_", case, ".rds"))
  emc <- EMC2:::restore_duplicates(readRDS(f)); set.seed(1)
  t <- system.time(emc <- suppressMessages(EMC2:::.emc_burn_mode_start(emc, n_cores = 12)))[["elapsed"]]
  st <- lapply(emc, `[[`, "burn_mode_stats")
  k <- unlist(lapply(st, `[[`, "pareto_k")); o <- unlist(lapply(st, `[[`, "outside"))
  out[[case]] <- data.frame(case = case, p = emc[[1]]$n_pars, jobs = length(k), hessians = sum(is.finite(k)),
    k_med = round(median(k, na.rm = TRUE), 2), k_gt_0.7 = sum(k > 0.7, na.rm = TRUE), k_gt_1 = sum(k > 1, na.rm = TRUE),
    outside_med = round(median(o, na.rm = TRUE), 2), outside_gt_0.5 = sum(o > 0.5, na.rm = TRUE), secs = round(t))
}
print(do.call(rbind, out), row.names = FALSE)

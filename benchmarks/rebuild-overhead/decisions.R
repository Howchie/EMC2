# Strict identity on a shared input, including immutable model environments.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2"); candidate <- new.env(parent = ns)
for (file in c("R/fitting.R", "R/objects.R", "R/rebuild_queue.R", "R/statistics.R"))
  sys.source(file, candidate)
keys <- c("split_along_dim", "get_summary_stat", "check_gd", "check_progress")
original <- setNames(lapply(keys, get, envir = ns), keys)
rows <- list()
for (path in commandArgs(TRUE)) {
  e <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(path)), quiet = TRUE)
  stage <- ns$get_last_stage(e)
  for (rep in 1:2) {
    values <- list()
    for (build in c("base", "new")) {
      for (key in keys) assignInNamespace(key,
        if (build == "base") original[[key]] else get(key, candidate), ns)
      set.seed(713); started <- proc.time()
      ans <- ns$check_progress(e, stage, iter = NULL,
        stop_criteria = list(max_gd = 1.1, min_es = 1, selection = c("mu", "sigma2", "alpha")),
        max_tries = 20, step_size = 100, n_cores = 9, verbose = FALSE, n_blocks = 1)
      values[[build]] <- list(ans, seed = .Random.seed)
      rows[[length(rows) + 1L]] <- data.frame(fixture = basename(path), rep, build,
        wall = (proc.time() - started)[["elapsed"]])
    }
    for (key in keys) assignInNamespace(key, original[[key]], ns)
    stopifnot(identical(values$base, values$new))
    cat(basename(path), "complete check result and RNG identical", rep, "\n")
  }
}
saveRDS(do.call(rbind, rows), Sys.getenv("DECISIONS_OUT"))

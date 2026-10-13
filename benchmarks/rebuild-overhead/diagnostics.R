# Paired diagnostic timings with only the array-splitting implementation changed.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2"); candidate <- new.env(parent = ns)
sys.source("R/objects.R", candidate)
original <- ns$split_along_dim
rows <- list()
for (path in commandArgs(TRUE)) {
  e <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(path)), quiet = TRUE)
  stage <- ns$get_last_stage(e)
  selection <- c("mu", "sigma2", "alpha")
  for (rep in 1:2) {
    values <- list()
    for (build in c("base", "new")) {
      assignInNamespace("split_along_dim", if (build == "base") original else candidate$split_along_dim, ns)
      set.seed(713)
      t <- proc.time()
      result <- ns$check_progress(e, stage, iter = NULL,
        stop_criteria = list(max_gd = 1.1, min_es = 1, selection = selection),
        max_tries = 20, step_size = 100, n_cores = 9, verbose = FALSE, n_blocks = 1)
      spent <- proc.time() - t
      values[[build]] <- list(result, seed = .Random.seed)
      rows[[length(rows) + 1L]] <- data.frame(fixture = basename(path), rep, build,
        wall = spent[["elapsed"]], master_cpu = sum(spent[c("user.self", "sys.self")]))
    }
    assignInNamespace("split_along_dim", original, ns)
    stopifnot(identical(values$base, values$new))
    cat(basename(path), "diagnostics identical", rep, "\n")
  }
}
rows <- do.call(rbind, rows)
print(rows, row.names = FALSE)
saveRDS(rows, Sys.getenv("DIAGNOSTICS_OUT"))

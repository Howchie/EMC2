# Compare full stopping decisions in separate installed baseline/candidate runs.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
cores <- as.integer(Sys.getenv("CORES", "9"))
repeats <- as.integer(Sys.getenv("REPEATS", "2"))
keep <- Sys.getenv("KEEP_RESULTS", "1") == "1"
collect_memory <- Sys.getenv("MEMORY", "0") == "1"
rows <- list(); results <- list(); dispatch <- list(); worker_cpu <- 0
rss <- function() {
  lines <- readLines("/proc/self/status", warn = FALSE)
  line <- grep("^VmRSS:", lines, value = TRUE)
  as.numeric(sub("^VmRSS:[[:space:]]*([0-9]+).*$", "\\1", line))
}
original <- ns$auto_mclapply
assignInNamespace("auto_mclapply", function(X, FUN, mc.cores, ...,
    mc.preschedule = TRUE, mc.set.seed = TRUE) {
  if (mc.cores <= 1L) return(original(X, FUN, mc.cores, ...,
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed))
  vals <- original(X, function(...) {
    started <- proc.time()
    value <- FUN(...)
    list(value = value, pid = Sys.getpid(), rss = rss(),
      cpu = sum((proc.time() - started)[c("user.self", "sys.self")]))
  }, mc.cores, ..., mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed)
  worker_cpu <<- worker_cpu + sum(vapply(vals, `[[`, numeric(1), "cpu"))
  dispatch[[length(dispatch) + 1L]] <<- data.frame(fixture, rep,
    pid = vapply(vals, `[[`, integer(1), "pid"),
    rss_kb = vapply(vals, `[[`, numeric(1), "rss"))
  lapply(vals, `[[`, "value")
}, ns)
for (path in commandArgs(TRUE)) {
  fixture <- basename(path)
  e <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(path)), quiet = TRUE)
  stage <- ns$get_last_stage(e)
  expected <- NULL
  for (rep in seq_len(repeats)) {
    if (collect_memory) invisible(gc())
    set.seed(713); started <- proc.time(); cpu0 <- worker_cpu
    ans <- ns$check_progress(e, stage, iter = NULL,
      stop_criteria = list(max_gd = 1.1, min_es = 1,
                           selection = c("mu", "sigma2", "alpha")),
      max_tries = 20, step_size = 100, n_cores = cores, verbose = FALSE, n_blocks = 1)
    elapsed <- proc.time() - started
    saved <- list(ans, seed = .Random.seed)
    if (is.null(expected)) expected <- saved else stopifnot(identical(saved, expected))
    rows[[length(rows) + 1L]] <- data.frame(fixture, rep,
      wall = elapsed[["elapsed"]], master_cpu = sum(elapsed[c("user.self", "sys.self")]),
      worker_cpu = worker_cpu - cpu0, master_rss_kb = rss(),
      fds = length(list.files("/proc/self/fd")))
    if (keep) results[[paste(fixture, rep)]] <- saved
  }
}
timings <- do.call(rbind, rows)
timings$idle_core_seconds <- with(timings, pmax(0, cores * wall - master_cpu - worker_cpu))
saveRDS(list(timings = timings, results = results, dispatch = do.call(rbind, dispatch)),
        Sys.getenv("DIAG_OUT"))
print(timings, row.names = FALSE)

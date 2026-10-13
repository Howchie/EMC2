# Paired component timings and identity checks with both implementations in
# one process. Use nine cores and pass the baseline private library.
lib <- Sys.getenv("EMC_LIB")
suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
rss <- function(pid = Sys.getpid()) {
  lines <- tryCatch(readLines(file.path("/proc", pid, "status"), warn = FALSE),
    error = function(e) character())
  line <- grep("^VmRSS:", lines, value = TRUE)
  if (length(line)) as.numeric(sub("^VmRSS:[[:space:]]*([0-9]+).*$", "\\1", line[1L])) * 1024 else NA_real_
}
candidate <- new.env(parent = ns)
for (file in c("R/fitting.R", "R/mode_start.R", "R/rebuild_queue.R", "R/shear_move.R", "R/collapsed_move.R", "R/objects.R"))
  sys.source(file, envir = candidate)
worker_total <- 0; dispatch <- list()
original_apply <- ns$auto_mclapply
assignInNamespace("auto_mclapply", function(X, FUN, mc.cores, ...,
    mc.preschedule = TRUE, mc.set.seed = TRUE) {
  if (mc.cores <= 1L) return(original_apply(X, FUN, mc.cores, ...,
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed))
  res <- original_apply(X, function(...) {
    t <- proc.time(); nested <- worker_total
    value <- FUN(...)
    list(value = value, cpu = sum((proc.time() - t)[c("user.self", "sys.self")]) +
      worker_total - nested, pid = Sys.getpid(), rss = rss())
  }, mc.cores, ..., mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed)
  worker_total <<- worker_total + sum(vapply(res, `[[`, numeric(1), "cpu"))
  dispatch[[length(dispatch) + 1L]] <<- data.frame(fixture, component, stage, rep, build,
    workers = length(unique(vapply(res, `[[`, integer(1), "pid"))),
    rss_max = max(vapply(res, `[[`, numeric(1), "rss"), na.rm = TRUE))
  lapply(res, `[[`, "value")
}, ns)
# The source overlay defines its own auto_mclapply in fitting.R. Route it
# through the same telemetry wrapper as the installed baseline namespace.
candidate$auto_mclapply <- ns$auto_mclapply
paths <- commandArgs(TRUE)
cores <- as.integer(Sys.getenv("CORES", "9"))
rows <- list(); skipped <- character()
for (path in paths) {
  fixture <- basename(path)
  e <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(path)), quiet = TRUE)
  e <- ns$.emc_prepare_chain_rng(e)
  # Raw data fixtures need a fit before they have a proposal window.
  if (!isTRUE(e[[1L]]$init) || is.null(e[[1L]]$samples$idx)) {
    cat("SKIP", fixture, "has no initialized sampler history\n"); skipped <- c(skipped, fixture)
    next
  }
  for (component in c(".emc_laplace_refresh", ".emc_check_likelihood_fits", ".emc_shear_atlas_add", "create_chain_proposals", ".emc_shear_orbit_steps", ".emc_collapsed_fits", "split_along_dim")) {
    original <- get(component, ns); replacement <- get(component, candidate)
    if (component %in% c(".emc_shear_orbit_steps", ".emc_shear_atlas_add") && !any(vapply(e, function(ch)
        is.list(ch$shear$plan), logical(1)))) next
    state <- e
    window <- ns$.emc_chain_proposal_window(e, NULL, "burn")
    if (component == ".emc_shear_atlas_add") {
      for (j in seq_along(state)) state[[j]] <- ns$.emc_shear_set_plan(state[[j]], window, "burn")
    }
    for (stage in if (component == "create_chain_proposals") c("burn", "adapt", "sample") else "refresh") {
      for (rep in 1:2) {
        answers <- list()
        for (build in c("base", "new")) {
          set.seed(713)
          fun <- if (build == "base") original else replacement
          before <- proc.time(); cpu0 <- worker_total
          ans <- if (component %in% c(".emc_laplace_refresh", ".emc_check_likelihood_fits")) fun(state, cores)
            else if (component == ".emc_shear_atlas_add") fun(state, window, cores, "burn") else if (component == ".emc_shear_orbit_steps" || component == ".emc_collapsed_fits")
            if (component == ".emc_collapsed_fits" && build == "base")
              fun(state, window) else
              fun(state, window, cores) else if (component == "split_along_dim") fun(e[[1L]]$samples$alpha, 2L)
            else fun(e, stage = stage)
          elapsed <- proc.time() - before
          answers[[build]] <- list(state = ans, seed = .Random.seed)
          rows[[length(rows) + 1L]] <- data.frame(fixture, component, stage, rep, build,
            wall = elapsed[["elapsed"]], master_cpu = sum(elapsed[c("user.self", "sys.self")]),
            worker_cpu = worker_total - cpu0)
        }
        same <- identical(answers$base, answers$new)
        cat(fixture, component, stage, rep, "identical:", same, "\n")
        if (!same) stop("Component state or RNG changed")
      }
    }
  }
}
out <- do.call(rbind, rows)
out$idle_core_seconds <- with(out, pmax(0, cores * wall - master_cpu - worker_cpu))
print(out, row.names = FALSE)
saveRDS(list(timings = out, dispatch = do.call(rbind, dispatch), skipped = skipped,
  cores = cores, telemetry_complete = TRUE), Sys.getenv("SCREEN_OUT"))

# Profile saved-state rebuilds or a short continuation in an isolated library.
lib <- Sys.getenv("EMC_LIB")
stopifnot(nzchar(lib))
suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
rss <- function(pid = Sys.getpid()) {
  lines <- tryCatch(readLines(file.path("/proc", pid, "status"), warn = FALSE),
    error = function(e) character())
  line <- grep("^VmRSS:", lines, value = TRUE)
  if (length(line)) as.numeric(sub("^VmRSS:[[:space:]]*([0-9]+).*$", "\\1", line[1L])) * 1024 else NA_real_
}
options(emc2.sampler_profile = TRUE)
rows <- list(); stack <- list(); jobs <- list(); worker_total <- 0
clock <- function() {
  t <- proc.time()
  c(wall = unname(t["elapsed"]), master = sum(t[c("user.self", "sys.self")]),
    worker = worker_total)
}
components <- c("run_emc", "reset_pm_settings", "fit_remove_samples", ".emc_burn_mode_start",
  ".emc_laplace_refresh", "add_proposals", "create_chain_proposals",
  ".emc_check_likelihood_fits", ".emc_shear_set_plans", ".emc_shear_atlas_add",
  ".emc_shear_build_curves", ".emc_shear_curve_validate", ".emc_shear_orbit_steps",
  ".emc_collapsed_set_plans", "sub_blocking", "create_eff_proposals", ".emc_alloc_rebuild",
  "subset.emc", "concat_emc", "check_progress", "strip_duplicates", "restore_duplicates")
for (name in components) local({
  original <- get(name, ns); label <- name
  wrapper <- function(...) {
    started <- clock(); level <- length(stack) + 1L
    stack[[level]] <<- list(label = label, child = c(wall = 0, master = 0, worker = 0))
    on.exit({
      inclusive <- clock() - started
      exclusive <- inclusive - stack[[level]]$child
      stack[[level]] <<- NULL
      if (level > 1L) stack[[level - 1L]]$child <<- stack[[level - 1L]]$child + inclusive
      rows[[length(rows) + 1L]] <<- data.frame(component = label, depth = level,
        wall = exclusive[["wall"]], master_cpu = exclusive[["master"]],
        worker_cpu = exclusive[["worker"]], inclusive_wall = inclusive[["wall"]])
    })
    original(...)
  }
  assignInNamespace(label, wrapper, ns)
})
# proc.time() child CPU is charged when a fork is reaped, which can happen
# after mclapply has returned. Return each job's own CPU with its result so
# costs cannot drift into the following component.
original_apply <- ns$auto_mclapply
assignInNamespace("auto_mclapply", function(X, FUN, mc.cores, ...,
    mc.preschedule = TRUE, mc.set.seed = TRUE) {
  if (mc.cores <= 1L) return(original_apply(X, FUN, mc.cores, ...,
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed))
  label <- if (length(stack)) stack[[length(stack)]]$label else "unattributed"
  sampling <- identical(label, "run_emc")
  started <- clock()
  if (sampling) {
    level <- length(stack) + 1L
    stack[[level]] <<- list(label = "sampling", child = c(wall = 0, master = 0, worker = 0))
  }
  res <- original_apply(X, function(...) {
    cpu0 <- sum(proc.time()[c("user.self", "sys.self")]); child0 <- worker_total
    value <- FUN(...)
    list(value = value, pid = Sys.getpid(),
      cpu = sum(proc.time()[c("user.self", "sys.self")]) - cpu0 + worker_total - child0,
      rss = rss())
  }, mc.cores, ..., mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed)
  cpu <- sum(vapply(res, `[[`, numeric(1), "cpu"))
  worker_total <<- worker_total + cpu
  jobs[[length(jobs) + 1L]] <<- data.frame(component = if (sampling) "sampling" else label,
    requested_cores = mc.cores, workers_used = length(unique(vapply(res, `[[`, integer(1), "pid"))),
    wall = clock()[["wall"]] - started[["wall"]], worker_cpu = cpu,
    rss_max = max(vapply(res, `[[`, numeric(1), "rss"), na.rm = TRUE))
  if (sampling) {
    inclusive <- clock() - started
    stack[[level]] <<- NULL
    stack[[level - 1L]]$child <<- stack[[level - 1L]]$child + inclusive
    rows[[length(rows) + 1L]] <<- data.frame(component = "sampling", depth = level,
      wall = inclusive[["wall"]], master_cpu = inclusive[["master"]],
      worker_cpu = inclusive[["worker"]], inclusive_wall = inclusive[["wall"]])
  }
  lapply(res, `[[`, "value")
}, ns)
memory <- list()
original_record <- ns$.emc_block_record
assignInNamespace(".emc_block_record", function(...) {
  args <- list(...)
  pids <- ns$.emc_profile_tree_pids()
  memory[[length(memory) + 1L]] <<- data.frame(block = args$block, pid = pids,
    rss = vapply(pids, rss, numeric(1)))
  original_record(...)
}, ns)
input <- Sys.getenv("EMC_IN"); stage <- Sys.getenv("STAGE", "burn")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(input)), quiet = TRUE)
if (!isTRUE(emc[[1L]]$init) || is.null(emc[[1L]]$samples$idx))
  stop("Input has no initialized sampler history; run preburn first")
rows <- list(); jobs <- list(); worker_total <- 0
seed <- as.integer(Sys.getenv("SEED", "1")); cores <- as.integer(Sys.getenv("CORES", "9"))
set.seed(seed)
mode <- Sys.getenv("MODE", "rebuild")
# Force a saved-state sample rebuild for profiling without collecting new draws.
# "first" uses the stage-entry history; "next" uses stored sample draws.
sample_rebuild <- Sys.getenv("SAMPLE_REBUILD", "")
if (mode == "rebuild" && stage == "sample" && nzchar(sample_rebuild)) {
  stopifnot(sample_rebuild %in% c("first", "next"))
  for (j in seq_along(emc)) {
    if (sample_rebuild == "first") {
      emc[[j]]$sample_proposals_next <- NULL
    } else {
      at <- emc[[j]]$sample_proposals_next
      if (is.null(at)) at <- 100L
      emc[[j]]$sample_transitions <- at
      emc[[j]]$sample_proposals_next <- at
    }
  }
}
if (mode == "run") {
  n <- as.integer(Sys.getenv("ITER", "20")); step <- as.integer(Sys.getenv("STEP", "2"))
  result <- run_emc(emc, stage, list(iter = n, selection = c("mu", "sigma2", "alpha")),
    step_size = step, max_tries = 20, cores_for_chains = min(3L, length(emc)),
    cores_per_chain = max(1L, cores %/% min(3L, length(emc))), trim = FALSE, verbose = TRUE)
} else {
  result <- ns$.emc_laplace_refresh(emc, cores)
  result <- ns$add_proposals(result, stage, cores, NULL, block_iter = 100L)
}
profile <- do.call(rbind, rows)
profile$busy_cores <- with(profile, ifelse(wall > 0, (master_cpu + worker_cpu) / wall, NA_real_))
profile$idle_core_seconds <- with(profile, pmax(0, cores * wall - master_cpu - worker_cpu))
print(profile, row.names = FALSE)
out <- Sys.getenv("PROF_OUT")
if (nzchar(out)) saveRDS(list(profile = profile, blocks = ns$.emc_block_profile(),
  state = result, seed = .Random.seed, library = find.package("EMC2"), input = input,
  stage = stage, cores = cores, jobs = do.call(rbind, jobs), memory = do.call(rbind, memory)), out)
